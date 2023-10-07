// Support fragment additive model.

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;

static int verbose = 0;

static int molecules_read = 0;

/*
  Some queries can be applied after implicit Hydrogens are made explicit. We
  look for the appropriate keyword in the control file
*/

static IWString hydrogen_keyword;

static IWString no_explicit_hydrogen_keyword;

static int queries_to_be_run_after_hydrogen_addition = 0;

static int make_implicit_hydrogens_explicit = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fp_tag;
static int bit_replicates = 5;

static int function_as_tdt_filter = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static IWString descriptor_name_prefix;

class Substituent_Model_Query : public Substructure_Hit_Statistics {
 private:
  typedef float contribution_t;

  Set_or_Unset<contribution_t>* _contribution;

  int _is_ring;

  int _apply_after_hydrogen_addition;

  //  We can speed things up by invalidating matched atoms

  int _invalidate_matched_atoms;

  int _must_match_unmatched_atoms;

  int _isotope;

  //  We want to be able temporarily turn on and turn off queries depending
  //  on the after hydrogen business

  int _active;

  //  There are multiple places from which our "name" can come, so we get it once

  IWString _comment;

  //  In the scenario where we add explicit hydrogens, we can specify whether or
  //  not the atoms matched by a query get explicit hydrogens or not

  int _attach_explicit_hydrogens;

  //  One of the components in the _comment field may be designated as the official name

  IWString _name;

  // private functions

  int _get_comment();
  void _default_values();

 public:
  Substituent_Model_Query();
  Substituent_Model_Query(const const_IWSubstring&);  // just needed for process_queries
  ~Substituent_Model_Query();

  int debug_print(IWString_and_File_Descriptor&) const;

  int initialise_query_specific_things();

  const IWString&
  name() const {
    return _name;
  }

  void
  set_name(const IWString& s) {
    _name = s;
  }

  int
  apply_after_hydrogen_addition() const {
    return _apply_after_hydrogen_addition;
  }

  int
  attach_explicit_hydrogens() const {
    return _attach_explicit_hydrogens;
  }

  int
  must_match_unmatched_atoms() const {
    return _must_match_unmatched_atoms;
  }

  int
  active() const {
    return _active;
  }

  void
  set_active(int s) {
    _active = s;
  }

  int initialise_this_many_models(int);

  int get_contribution_for_model(const IWString& model_prefix, int ndx);
  int get_text_token(const IWString& model_prefix, IWString& s);

  int
  isotope() const {
    return _isotope;
  }

  int
  is_ring() const {
    return _is_ring;
  }

  int
  invalidate_atoms_matched() const {
    return _invalidate_matched_atoms;
  }

  void
  set_invalidate_atoms_matched(int i) {
    _invalidate_matched_atoms = i;
  }

  int increment_results(int, float*) const;

  int
  has_value_for_model(int i) const {
    return _contribution[i].is_set();
  }

  float value_for_model(int) const;
};

void
Substituent_Model_Query::_default_values() {
  _contribution = nullptr;

  _is_ring = 0;

  _invalidate_matched_atoms = 0;

  _must_match_unmatched_atoms = 1;

  _isotope = 0;

  _apply_after_hydrogen_addition = 0;

  _active = 1;

  _attach_explicit_hydrogens = 1;

  return;
}

Substituent_Model_Query::Substituent_Model_Query() {
  _default_values();

  return;
}

Substituent_Model_Query::Substituent_Model_Query(const const_IWSubstring& notused) {
  _default_values();

  return;
}

Substituent_Model_Query::~Substituent_Model_Query() {
  if (NULL != _contribution) {
    delete[] _contribution;
  }

  return;
}

int
Substituent_Model_Query::_get_comment() {
  _comment = Substructure_Hit_Statistics::comment();

  if (0 == _comment.length()) {
    _comment = _things[0]->comment();
    if (0 == _comment.length()) {
      cerr << "Substituent_Model_Query::_get_comment:no name!\n";
      return 0;
    }
  }

  return 1;
}

int
Substituent_Model_Query::initialise_this_many_models(int n) {
  assert(n > 0);
  assert(NULL == _contribution);

  if (0 == _comment.length()) {
    if (!_get_comment()) {
      return 0;
    }
  }

  _contribution = new Set_or_Unset<contribution_t>[n];

  if (NULL == _contribution) {
    cerr << "Substituent_Model_Query::initialise_this_many_models:bad news, cannot "
            "allocate "
         << n << " contributions\n";
    return 0;
  }

  return 1;
}

int
Substituent_Model_Query::initialise_query_specific_things() {
  if (0 == _comment.length()) {
    if (!_get_comment()) {
      return 0;
    }
  }

  int i = 0;
  const_IWSubstring token;
  while (_comment.nextword(token, i)) {
    if ("INVALIDATE" == token) {
      _invalidate_matched_atoms = 1;
    } else if ("IS_RING" == token) {
      _is_ring = -523332;  // some unlikely number, but the same for all objects, and
                           // unlikely to conflict with any IS_RING=nnn directive
    } else if (token.starts_with("IS_RING=")) {
      token.remove_leading_chars(8);
      if (!token.numeric_value(_is_ring)) {
        cerr << "Substituent_Model_Query::initialise_query_specific_things:invalid ring "
                "specification '"
             << _comment << "'\n";
        return 0;
      }
    } else if (hydrogen_keyword == token) {
      _apply_after_hydrogen_addition = 1;
    } else if (no_explicit_hydrogen_keyword == token) {
      _attach_explicit_hydrogens = 0;
    } else if (token.starts_with("NAME=")) {
      token.remove_leading_chars(5);
      _name = token;
    } else if (token.starts_with("ISO=")) {
      token.remove_leading_chars(4);
      if (!token.numeric_value(_isotope) || _isotope <= 0) {
        cerr << "Substituent_Model_Query::initialise_query_specific_things:invalid "
                "isotope '"
             << _comment << "'\n";
        return 0;
      }
    } else if ("OKOVERLAP" == token) {
      _must_match_unmatched_atoms = 0;
    }
  }

  return 1;
}

int
Substituent_Model_Query::get_contribution_for_model(const IWString& model_prefix,
                                                    int ndx) {
  if (0 == _comment.length()) {
    if (!_get_comment()) {
      return 0;
    }
  }

  int i = 0;
  const_IWSubstring token;

  while (_comment.nextword(token, i)) {
    if (!token.starts_with(model_prefix)) {
      continue;
    }

    token.remove_up_to_first('=');

    contribution_t d;
    if (!token.numeric_value(d)) {
      cerr << "Substituent_Model_Query::get_contribution_for_model:invalid numeric for "
              "model '"
           << model_prefix << "' '" << token << "'\n";
      return 0;
    }

    _contribution[ndx].set(d);

    return 1;
  }

  cerr << "Substituent_Model_Query::get_contribution_for_model:no contribution to model '"
       << model_prefix << "'\n";
  cerr << _comment << '\n';

  return 0;
}

int
Substituent_Model_Query::get_text_token(const IWString& model_prefix, IWString& s) {
  if (0 == _comment.length()) {
    if (!_get_comment()) {
      return 0;
    }
  }

  int i = 0;
  const_IWSubstring token;

  while (_comment.nextword(token, i)) {
    if (!token.starts_with(model_prefix)) {
      continue;
    }

    token.remove_up_to_first('=');

    s = token;

    return 1;
  }

  cerr << "Substituent_Model_Query::get_text_token:nothing matching '" << model_prefix
       << "' in '" << _comment << "'\n";
  return 0;
}

int
Substituent_Model_Query::increment_results(int number_models, float* zresult) const {
  contribution_t f;  // scope here for efficiency

  for (int i = 0; i < number_models; i++) {
    if (_contribution[i].value(f)) {
      zresult[i] += static_cast<float>(f);
    }
  }

  return 1;
}

float
Substituent_Model_Query::value_for_model(int ndx) const {
  contribution_t f;
  if (_contribution[ndx].value(f)) {
    return static_cast<float>(f);
  }

  return static_cast<float>(0.0);
}

static resizable_array_p<Substituent_Model_Query> queries;

static int nq = 0;  // same as queries.number_elements();

/*
  For each model, we have text in the query name with the contribution for each model
*/

static int number_models = 0;

static IWString* model_name = nullptr;

static int isotopes_present = 0;

static Accumulator<float>* acc = nullptr;

/*
  We may also insist on different queries hitting completely different atoms

  Jan 2003. Ran into the problem of benzene as a query. What happens when
  we get to Napthalene. If we insist on no previously matched atoms
  being hit, then the first benzene match will prevent the second one
  from being perceived.

  Therefore if 1 == queries_must_not_hit_previously_matched_atoms then
  we make sure that each new embedding does not overlap with any
  atoms matched by previous queries

  If 2 == queries_must_not_hit_previously_matched_atoms then we check
  for any previously matched atom - including atoms matched by the
  current query
*/

static int queries_must_not_hit_previously_matched_atoms = 0;

/*
  It will be convenient to label the matched atoms
*/

static Molecule_Output_Object stream_for_labelled_atoms;

static IWString_and_File_Descriptor stream_for_atomic_values;

/*
  Normal descriptors are the number of hits to each query
*/

static IWString_and_File_Descriptor stream_for_descriptors;

/*
  Chavali wanted the descriptor to be the contribution made by each query
*/

static IWString_and_File_Descriptor stream_for_per_query_values;

static int* hits_to_query = nullptr;

static int write_to_stdout = 1;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off

  cerr << "Implements one or more substituent models\n";
  cerr << "  -M <descr>     produce output for model <descr> - must be defined in query file\n";
  cerr << "  -q <query>     specify one or more queries that define the substituents\n";
  cerr << "  -q S:<file>    file of smarts\n";
  cerr << "  -q F:<file>    file of queries\n";
  cerr << "  -I             invalidate matched atoms - speeds things up\n";
  cerr << "  -u             each query must hit atoms not matched by any previous query\n";
  cerr << "                 two -u's mean embeddings of the current query cannot overlap\n";
  cerr << "  -L <fname>     write labelled atoms to <fname>\n";
  cerr << "  -W <fname>     dump atomic contributions to <fname>\n";
  cerr << "  -D <fname>     create descriptor file - descriptors are counts of times hit\n";
  cerr << "  -B <fname>     create descriptor file - descriptors are per query contributions\n";
  cerr << "  -h             make implicit Hydrogens explicit before any processing\n";
  cerr << "  -H <string>    queries with <string> in them will only be run after explicit Hydrogens are added\n";
  cerr << "  -H noxh=<string> atoms matched by queries with <string> in them will not get explicit Hydrogens\n";
  cerr << "  -J <tag>       convert results to fingerprint\n";
  cerr << "  -f             function as a tdt filter\n";
  cerr << "  -P <prefix>    prefix for descriptor names\n";
  cerr << "  -l             reduce to largest fragment\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -E <sym>       create element with symbol <sym>\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

template <typename T>
int
do_write_multi_component_values(const IWString& mname, int n, float zresult, const T* f,
                                IWString_and_File_Descriptor& output) {
  append_first_token_of_name(mname, output);

  output << ' ' << zresult;

  for (int i = 0; i < nq; i++) {
    output << ' ' << f[i];
  }

  output << '\n';

  return output.good();
}

static int
do_write_per_query_contribution(const IWString& mname, float zresult,
                                const float* per_query_contribution,
                                IWString_and_File_Descriptor& output) {
  return do_write_multi_component_values(mname, nq, zresult, per_query_contribution,
                                         output);
}

static int
do_write_descriptors(const IWString& mname, float zresult, const int* hits_to_query,
                     IWString_and_File_Descriptor& output) {
  return do_write_multi_component_values(mname, nq, zresult, hits_to_query, output);
}

static int
write_atomic_contributions(Molecule& m, int* already_hit, const float* zresult,
                           IWString_and_File_Descriptor& os) {
  const IWString& mname = m.name();

  int matoms = m.natoms();

  os << mname << ' ' << matoms << " atoms\n";

  for (int i = 0; i < matoms; i++) {
    os << i << ' ';

    const Atom* a = m.atomi(i);

    os << a->atomic_symbol() << ' ' << a->ncon() << " connections ";

    int j = already_hit[i];

    if (0 == j) {
      os << ". . .";
    } else {
      const Substituent_Model_Query* q = queries[j - 1];
      os << j << ' ' << zresult[i] << ' ' << q->name();
    }

    os << '\n';
  }

  os << "|\n";

  os.write_if_buffer_holds_more_than(32768);

  return os.good();
}

static int
write_labelled_atoms(Molecule& m, const int* already_hit, const float* zresult,
                     Molecule_Output_Object& output) {
  if (isotopes_present) {
    int matoms = m.natoms();
    for (int i = 0; i < matoms; i++) {
      if (0 == already_hit[i]) {
        continue;
      }

      int query_number = already_hit[i] - 1;
      m.set_isotope(i, queries[query_number]->isotope());
    }
  }

  IWString initial_name = m.name();  // we need to temporarily switch the name

  IWString buffer(initial_name);

  for (int i = 0; i < number_models; i++) {
    buffer << ' ' << model_name[i] << ' ' << zresult[i];
  }

  m.set_name(buffer);

  int rc = output.write(m);

  m.set_name(initial_name);

  return rc;
}

static int
match_rejected_for_some_reason_is_ring(const Set_of_Atoms& embedding,
                                       const int* already_hit, int f, int is_ring) {
  if (1 == queries_must_not_hit_previously_matched_atoms) {
    int unmatched_atoms_hit =
        0;  // we must hit at least some atoms not previously matched

    int n = embedding.number_elements();
    for (int i = 0; i < n; i++) {
      atom_number_t j = embedding[i];

      if (0 == already_hit[j]) {
        unmatched_atoms_hit++;
        continue;
      }

      if (already_hit[j] ==
          f) {  // hit another instance of the current query, ignore these
        continue;
      }

      int k = already_hit[j];

      const Substituent_Model_Query* q = queries[k - 1];

      if (q->is_ring() == is_ring) {
        continue;
      }

      return 1;  // hit another query that was not a ring of the same number
    }

    return 0 == unmatched_atoms_hit;
  }

  assert(2 == queries_must_not_hit_previously_matched_atoms);  // we must not hit even our
                                                               // own query

  // Do something a little tricky here. We temporarily switch the value of
  // queries_must_not_hit_previously_matched_atoms and pass a value of F that will not be
  // present in the ALREADY_HIT array

  queries_must_not_hit_previously_matched_atoms = 1;

  int rc = match_rejected_for_some_reason_is_ring(embedding, already_hit, f + 1, is_ring);

  queries_must_not_hit_previously_matched_atoms = 2;

  return rc;
}

static int
match_rejected_for_some_reason(const Set_of_Atoms& embedding, const int* already_hit,
                               int f, int isring) {
  if (0 == queries_must_not_hit_previously_matched_atoms) {  // no checking
    return 0;
  }

  if (isring) {
    return match_rejected_for_some_reason_is_ring(embedding, already_hit, f, isring);
  }

  if (1 == queries_must_not_hit_previously_matched_atoms) {
    int n = embedding.number_elements();
    for (int i = 0; i < n; i++) {
      atom_number_t j = embedding[i];

      if (already_hit[j] > 0 &&
          already_hit[j] <
              f) {  // hit by a previous query - no check on hit by the current query
        return 1;
      }
    }

    return 0;
  }

  assert(2 == queries_must_not_hit_previously_matched_atoms);

  if (embedding.any_members_set_in_array(already_hit)) {
    return 1;  // rejected
  } else {
    return 0;  // not rejected
  }
}

// template class resizable_array_p<Substituent_Model_Query>;
// template class resizable_array_base<Substituent_Model_Query *>;
// template int process_queries<Substituent_Model_Query>(Command_Line &,
// resizable_array_p<Substituent_Model_Query> &, int, char); template int
// process_cmdline_token<Substituent_Model_Query>(const_IWSubstring &,
// resizable_array_p<Substituent_Model_Query> &, int); template int
// queries_from_file<Substituent_Model_Query>(const_IWSubstring const &,
// resizable_array_p<Substituent_Model_Query> &, int, int);

static int
run_queries(Molecule_to_Match& target, float* zresult, int* already_hit, float* each_atom,
            float* per_query_contribution) {
  for (int i = 0; i < nq; i++) {
    if (!queries[i]->active()) {
      continue;
    }

    Substructure_Results sresults;
    sresults.set_remove_invalid_atom_numbers_from_new_embeddings(1);

    int nhits = queries[i]->substructure_search(target, sresults);
    if (0 == nhits) {
      continue;
    }

    if (verbose > 2) {
      cerr << nhits << " hits to query " << i << " '" << queries[i]->comment() << "'\n";
    }

    int hits_used = 0;

    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* ej = sresults.embedding(j);

      if (!queries[i]->must_match_unmatched_atoms()) {
        ;
      } else if (match_rejected_for_some_reason(*ej, already_hit, i + 1,
                                                queries[i]->is_ring())) {
        continue;
      }

      hits_used++;

#ifdef DEBUG_RUN_QUERIES
      cerr << "Before updating already_hit: ring? " << queries[i]->is_ring() << '\n';
      for (int k = 0; k < target.natoms(); k++) {
        cerr << " atom " << k << " already_hit " << already_hit[k] << '\n';
      }
#endif

      ej->set_vector(already_hit, i + 1);

      queries[i]->increment_results(number_models, zresult);

      if (NULL != each_atom) {
        ej->set_vector(each_atom, queries[i]->value_for_model(0));
      }

      //    cerr << "Query " << i << " value " << queries[i]->value_for_model(0) << '\n';

      if (NULL != per_query_contribution) {
        per_query_contribution[i] += queries[i]->value_for_model(0);
      }

      if (queries[i]->invalidate_atoms_matched()) {
        target.invalidate(*ej);
      }
    }

    hits_to_query[i] = hits_used;

    if (verbose > 2) {
      cerr << hits_used << " of " << nhits << " embeddings used\n";
    }
  }

  return 1;
}

static int
add_explicit_hydrogens(Molecule& m, const int* already_hit) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    int j = already_hit[i];

    if (0 == j) {  // atom not matched by any query, will add hydrogens
      ;
    } else  // was hit. Does the matching query allow addition of an explicit Hydrogen
    {
      const Substituent_Model_Query* q = queries[j - 1];
      if (!q->attach_explicit_hydrogens()) {
        continue;
      }
    }

    m.make_implicit_hydrogens_explicit(i);
    rc++;
  }

  return rc;
}

static void
set_active(int after_h, int before_h) {
  for (int i = 0; i < nq; i++) {
    Substituent_Model_Query* q = queries[i];

    if (q->apply_after_hydrogen_addition()) {
      q->set_active(after_h);
    } else {
      q->set_active(before_h);
    }
  }

  return;
}

static int
substituent_model(Molecule& m, float* zresult, int* already_hit, float* each_atom,
                  float* per_query_contribution) {
  set_vector(hits_to_query, nq, 0);

  Molecule_to_Match target(&m);

  if (hydrogen_keyword.length() > 0) {
    set_active(0, 1);
  }

  run_queries(target, zresult, already_hit, each_atom, per_query_contribution);

  if (0 == hydrogen_keyword.length()) {
    ;
  } else if (0 == queries_to_be_run_after_hydrogen_addition) {
    ;
  } else {
    set_active(1, 0);

    add_explicit_hydrogens(m, already_hit);

    Molecule_to_Match t2(&m);

    run_queries(t2, zresult, already_hit, each_atom, per_query_contribution);
  }

  if (stream_for_descriptors.is_open()) {
    do_write_descriptors(m.name(), zresult[0], hits_to_query, stream_for_descriptors);
  }

  if (stream_for_labelled_atoms.active()) {
    write_labelled_atoms(m, already_hit, zresult, stream_for_labelled_atoms);
  }

  if (NULL != each_atom) {
    write_atomic_contributions(m, already_hit, each_atom, stream_for_atomic_values);
  }

  if (NULL != per_query_contribution) {
    do_write_per_query_contribution(m.name(), zresult[0], per_query_contribution,
                                    stream_for_per_query_values);
  }

  return 1;
}

/*
  This is hard coded to do the sa and sb models.
  sa appears to run between 0 and 5
  sb appears to run between -1 and 9
*/

static int
do_fingerprint_output(Molecule& m, int number_models, const float* zresult,
                      IWString_and_File_Descriptor& output) {
  if (!function_as_tdt_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < number_models; i++) {
    int v;
    if (0 == i) {
      v = static_cast<int>(zresult[0] * 2.0 + 0.4999) + 1;
    } else {
      v = static_cast<int>(zresult[0] + 1.4999) + 1;
    }

    if (v < 1) {
      v = 1;
    }

    for (int j = 0; j < bit_replicates; j++) {
      sfc.hit_bit(j * number_models + i, v);
    }
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(fp_tag, tmp);

  output << tmp << '\n';

  if (!function_as_tdt_filter) {
    output << "|\n";
  }

  return 1;
}

static int
do_output(Molecule& m, int number_models, const float* zresult,
          IWString_and_File_Descriptor& output) {
  if (fp_tag.length()) {
    return do_fingerprint_output(m, number_models, zresult, output);
  }

  IWString buffer;
  if (write_to_stdout) {
    if (0 == m.name().length()) {
      buffer << "M" << molecules_read;
    } else {
      append_first_token_of_name(m.name(), buffer);
    }
  }

  for (int i = 0; i < number_models; i++) {
    if (buffer.length()) {
      buffer += ' ';
    }

    buffer.append_number(zresult[i], 4);

    acc[i].extra(zresult[i]);
  }

  if (write_to_stdout) {
    buffer += '\n';
    output << buffer;

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

static void
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  return;
}

static int
substituent_model(Molecule& m, IWString_and_File_Descriptor& output) {
  molecules_read++;

  preprocess(m);

  int matoms = m.natoms();

  if (hydrogen_keyword.length() > 0) {
    matoms += m.implicit_hydrogens();
  }

  int* tmp = new_int(matoms, 0);
  std::unique_ptr<int[]> free_tmp(tmp);

  float* zresult = new_float(number_models);
  std::unique_ptr<float[]> free_zresult(zresult);

  float* each_atom;

  if (stream_for_atomic_values.is_open()) {
    each_atom = new_float(matoms);
  } else {
    each_atom = nullptr;
  }

  float* per_query_contribution;

  if (stream_for_per_query_values.is_open()) {
    per_query_contribution = new_float(nq);
  } else {
    per_query_contribution = nullptr;
  }

  int rc = substituent_model(m, zresult, tmp, each_atom, per_query_contribution);

  if (rc) {
    do_output(m, number_models, zresult, output);
  }

  if (NULL != each_atom) {
    delete[] each_atom;
  }

  if (NULL != per_query_contribution) {
    delete[] per_query_contribution;
  }

  return rc;
}

static int
substituent_model(data_source_and_type<Molecule>& input,
                  IWString_and_File_Descriptor& output) {
  Molecule* m;

  while (NULL != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (0 == m->natoms()) {
      cerr << "Skipping molecule with no atoms\n";
      continue;
    }

    if (!substituent_model(*m, output)) {
      return 0;
    }
  }

  return output.good();
}

static int
substituent_model_record(const_IWSubstring buffer,  // local copy
                         IWString_and_File_Descriptor& output) {
  buffer.remove_leading_chars(smiles_tag.length());
  assert(buffer.ends_with('>'));
  buffer.chop();

  Molecule m;

  if (!m.build_from_smiles(buffer)) {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return substituent_model(m, output);
}

static int
substituent_model(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!buffer.starts_with(smiles_tag)) {
      output << buffer << '\n';
    } else if (!substituent_model_record(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
substituent_model(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return substituent_model(input, output);
}

static int
substituent_model(const char* fname, FileType input_type,
                  IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  return substituent_model(input, output);
}

/*
  Examine the query names and determine the contribution to each model
*/

static int
determine_contributions(const IWString* model_prefix, int number_models,
                        const resizable_array_p<Substituent_Model_Query>& queries) {
  assert(nq = queries.number_elements());

  int rc = 1;

  for (int i = 0; i < number_models; i++) {
    const IWString& m = model_prefix[i];

    for (int j = 0; j < nq; j++) {
      if (!queries[j]->get_contribution_for_model(m, i)) {
        cerr << "Query " << j << " '" << queries[j]->comment()
             << "' cannot initialise model '" << m << "'\n";
        rc = 0;
      }
    }
  }

  return rc;
}

static int
write_header(const IWString* h, int number_models, IWString_and_File_Descriptor& output) {
  IWString header;
  for (int i = 0; i < number_models; i++) {
    const_IWSubstring hi = h[i];
    assert(hi.ends_with('='));

    hi.chop();

    if (descriptor_name_prefix.length()) {
      header << ' ' << descriptor_name_prefix << hi;
    } else {
      header.append_with_spacer(hi);
    }
  }

  if (write_to_stdout) {
    output << "Name " << header << '\n';
  }

  return output.good();  // even if not used
}

static int
determine_contributions(Command_Line& cl, char flag, int number_models,
                        const resizable_array_p<Substituent_Model_Query>& queries,
                        IWString_and_File_Descriptor& output) {
  assert(number_models > 0);

  model_name = new IWString[number_models];

  for (int i = 0; i < number_models; i++) {
    model_name[i] = cl.string_value(flag, i);
    if (!model_name[i].ends_with('=')) {
      model_name[i] += '=';
    }
  }

  if (!determine_contributions(model_name, number_models, queries)) {
    return 0;
  }

  if (0 == fp_tag.length()) {
    write_header(model_name, number_models, output);
  }

  return 1;
}

static int
write_descriptor_header(const const_IWSubstring& model_name, int n,
                        IWString_and_File_Descriptor& output) {
  output << "Name ";

  const_IWSubstring tmp(model_name);
  if (tmp.ends_with('=')) {
    tmp.chop();
  }

  output << tmp;

  for (int i = 0; i < n; i++) {
    const Substituent_Model_Query* q = queries[i];

    output << ' ' << q->name();
  }

  output << '\n';

  return output.good();
}

static int
open_descriptor_file(IWString& fname, const const_IWSubstring& model_name, int n,
                     IWString_and_File_Descriptor& output) {
  if (number_models > 1) {
    cerr << "Sorry, writing a descriptor file only works with one model present, see "
            "Ian\n";
    return 4;
  }

  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open descriptor file '" << fname << "'\n";
    return 0;
  }

  return write_descriptor_header(model_name, n, output);
}

static int
substituent_model(int argc, char** argv) {
  Command_Line cl(argc, argv, "vq:A:E:o:i:M:L:nIuD:W:B:H:hJ:fg:lP:p:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  (void)process_standard_aromaticity_options(cl, verbose > 1);

  (void)process_elements(cl, 'E', verbose > 1);

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will strip to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      usage(6);
    }
  }

  if (cl.option_present('P')) {
    cl.value('P', descriptor_name_prefix);

    if (verbose) {
      cerr << "Descriptors produced with prefix '" << descriptor_name_prefix << "'\n";
    }
  }

  if (cl.option_present('H')) {
    int i = 0;
    const_IWSubstring h;
    while (cl.value('H', h, i++)) {
      if (h.starts_with("noxh=")) {
        h.remove_leading_chars(5);
        no_explicit_hydrogen_keyword = h;

        if (verbose) {
          cerr << "Queries marked with the '" << no_explicit_hydrogen_keyword
               << "' keyword will not get explicit hydrogens\n";
        }
      } else {
        hydrogen_keyword = h;
        if (verbose) {
          cerr << "Queries containing the keyword '" << hydrogen_keyword
               << "' will be run after explicit hydrogen addition\n";
        }
      }
    }
  }

  if (cl.option_present('h')) {
    if (cl.option_present('H')) {
      cerr << "The -h and -H options are mutually inconsistent\n";
      usage(5);
    }

    make_implicit_hydrogens_explicit = 1;

    if (verbose) {
      cerr << "Implicit hydrogens will be made explicit\n";
    }
  }

  if (!cl.option_present('q')) {
    cerr << "Must specify one or more substitutent queries via the -q option\n";
    usage(2);
  }

  if (!process_queries(cl, queries, verbose, 'q')) {
    cerr << "Cannot read queries(-q option)\n";
    return 8;
  }

  nq = queries.number_elements();

  if (verbose) {
    cerr << "Read " << nq << " queries\n";
  }

  queries_to_be_run_after_hydrogen_addition = 0;

  for (int i = 0; i < nq; i++) {
    Substituent_Model_Query* q = queries[i];

    const IWString& c = q->comment();

    if (0 == c.length()) {
      cerr << "Zero length query comment, impossible\n";
      return i + 1;
    }

    q->set_find_unique_embeddings_only(1);

    q->initialise_query_specific_things();

    if (q->apply_after_hydrogen_addition()) {
      queries_to_be_run_after_hydrogen_addition++;
    }

    if (0 == q->name().length()) {
      IWString tmp(c);
      tmp.truncate_at_first(' ');
      q->set_name(tmp);
    }
  }

  if (verbose) {
    cerr << queries_to_be_run_after_hydrogen_addition
         << " queries to be run after Hydrogen addition\n";
  }

  number_models = cl.option_count('M');

  if (0 == number_models) {
    cerr << "Must specify one or more models via the -M option\n";
    usage(16);
  }

  for (int i = 0; i < nq; i++) {
    queries[i]->initialise_this_many_models(number_models);
  }

  if (cl.option_present('J')) {
    cl.value('J', fp_tag);

    if (!fp_tag.ends_with('<')) {
      fp_tag << '<';
    }

    if (verbose) {
      cerr << "Results written as fingerprint with tag '" << fp_tag << "'\n";
    }

    if (cl.option_present('f')) {
      function_as_tdt_filter = 1;

      if (verbose) {
        cerr << "Will function as a TDT filter\n";
      }
    }

    if (cl.option_present('p')) {
      if (!cl.value('p', bit_replicates) || bit_replicates < 1) {
        cerr << "The number of bit replicates (-p) option must be a whole +ve number\n";
        usage(2);
      }

      if (verbose) {
        cerr << "Will produce " << bit_replicates << " bit replicates\n";
      }
    }
  }

  IWString_and_File_Descriptor output(1);

  acc = new Accumulator<float>[number_models];

  if (!determine_contributions(cl, 'M', number_models, queries, output)) {
    cerr << "Cannot determine contribution of each query to each model\n";
    return 4;
  }

  for (int i = 0; i < number_models; i++) {
    int number_non_zero = 0;

    for (int j = 0; j < nq; j++) {
      if (queries[j]->has_value_for_model(i)) {
        number_non_zero++;
      }
    }

    if (0 == number_non_zero) {
      cerr << "Yipes, no queries have a value for model '" << model_name[i] << "'\n";
      return 3;
    }
  }

  if (cl.option_present('I')) {
    for (int i = 0; i < nq; i++) {
      queries[i]->set_invalidate_atoms_matched(1);
    }

    if (verbose) {
      cerr << "Atoms matched by a query will be prevented from matched other queries\n";
    }
  }

  if (cl.option_present('u')) {
    queries_must_not_hit_previously_matched_atoms = cl.option_count('u');

    if (0 == verbose) {
      ;
    } else if (1 == queries_must_not_hit_previously_matched_atoms) {
      cerr << "Each query must hit atoms not previously matched by any other query\n";
    } else if (2 == queries_must_not_hit_previously_matched_atoms) {
      cerr << "Each query must hit atoms not previously matched by any embedding\n";
    } else {
      cerr << "Can have at most two -u options\n";
      usage(4);
    }
  }

  if (cl.option_present('n')) {
    write_to_stdout = 0;

    if (verbose) {
      cerr << "Normal output suppressed\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (function_as_tdt_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 5;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(7);
  }

  if (cl.option_present('L')) {
    if (!cl.option_present('o')) {
      stream_for_labelled_atoms.add_output_type(FILE_TYPE_SMI);
    } else if (!stream_for_labelled_atoms.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s) for -L file\n";
      usage(6);
    }

    const IWString& fname = cl.string_value('L');

    if ('-' == fname) {
      ;
    } else if (stream_for_labelled_atoms.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) '" << fname << "'\n";
      return 6;
    }

    if (!stream_for_labelled_atoms.new_stem(fname)) {
      cerr << "Cannot open stream for isotopically labelled atoms '" << fname << "'\n";
      return 8;
    }

    if (verbose) {
      cerr << "Labelled molecules written to '" << fname << "'\n";
    }

    isotopes_present = 0;

    for (int i = 0; i < nq; i++) {
      if (queries[i]->isotope() > 0) {
        isotopes_present++;
      }
    }

    if (0 == isotopes_present) {
      cerr << "Sorry, none of the queries have isotopes specified\n";
      return 3;
    }

    if (verbose) {
      cerr << isotopes_present << " of " << nq
           << " queries have isotopic label specifications\n";
    }
  }

  if (cl.option_present('W')) {
    if (number_models > 1) {
      cerr << "Sorry, dumping atomic contributions only works with one model present. "
              "See Ian\n";
      return 6;
    }

    IWString d = cl.string_value('W');

    if (!stream_for_atomic_values.open(d.null_terminated_chars())) {
      cerr << "Bad news, cannot open stream for atomic contributions '" << d << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Will write atomic contributions to '" << d << "'\n";
    }
  }

  if (cl.option_present('D')) {
    IWString d = cl.string_value('D');

    if (!open_descriptor_file(d, model_name[0], nq, stream_for_descriptors)) {
      return 5;
    }

    if (verbose) {
      cerr << "Descriptors as counts written to '" << d << "'\n";
    }
  }

  if (cl.option_present('B')) {
    IWString b = cl.string_value('B');

    if (!open_descriptor_file(b, model_name[0], nq, stream_for_per_query_values)) {
      return 6;
    }

    if (verbose) {
      cerr << "Descriptors as contributions written to '" << b << "'\n";
    }
  }

  hits_to_query = new int[nq];
  std::unique_ptr<int[]> free_hits_to_query(hits_to_query);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (function_as_tdt_filter) {
      if (!substituent_model(cl[i], output)) {
        rc = i + 1;
        break;
      }
    } else if (!substituent_model(cl[i], input_type, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < number_models; i++) {
      const Accumulator<float>& acci = acc[i];

      cerr << "For model " << i << " values between " << acci.minval() << " and "
           << acci.maxval() << '\n';
    }

    if (verbose > 1) {
      for (int i = 0; i < nq; i++) {
        queries[i]->report(cerr, verbose > 1);
      }
    }
  }

  delete[] acc;

  return rc;
}

int
main(int argc, char** argv) {
  int rc = substituent_model(argc, argv);

  return rc;
}
