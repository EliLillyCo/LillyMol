/*
  Compute degree of overlap in a substructure match
*/

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static Molecule_to_Query_Specifications mqs;

static int allow_self_matches = 1;

extending_resizable_array<int> queries_hitting;

static float upper_fraction_matched = static_cast<float>(1.0);

static float lower_fraction_matched = static_cast<float>(0.0);

static int matches_below_lower_fraction_matched = 0;

static int matches_above_upper_fraction_matched = 0;

static Accumulator<float> highest_fraction;

static int write_molecules_with_no_neighbours = 0;

static int output_distances = 1;

static int compute_fraction_on_first_embedding_only = 1;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -q <fname>    query file specification\n";
  cerr << "  -f            only find one embedding of the query (default is to find all)\n";
  cerr << "  -h            avoid self matches\n";
  cerr << "  -T <fraction> upper distance cutoff\n";
  cerr << "  -t <fraction> lower distance cutoff\n";
  cerr << "  -z            write needle molecules that have no matches\n";
  cerr << "  -r r          ring atoms preserve their ring membership in matches\n";
  cerr << "  -r c          chain atoms will only match chain atoms in matches\n";
  cerr << "  -a            use all embeddings to determine query coverage\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

class Smiles_ID_and_Fraction
{
 private:
  const IWString _smiles;
  const IWString _id;
  const float _fraction;

 public:
  Smiles_ID_and_Fraction(const IWString& s, const IWString& i, float f)
      : _smiles(s), _id(i), _fraction(f){};

  const IWString&
  smiles() const
  {
    return _smiles;
  }

  const IWString&
  id() const
  {
    return _id;
  }

  float
  fraction_atoms_matched() const
  {
    return _fraction;
  }
};

class Needle : public Substructure_Hit_Statistics
{
 private:
  IWString _smiles;

  resizable_array_p<Smiles_ID_and_Fraction> _neighbours;

 public:
  void
  set_smiles(const IWString& s)
  {
    _smiles = s;
  }

  int perform_search(Molecule& m, Molecule_to_Match& target);

  int do_write(IWString_and_File_Descriptor&);
};

static float
fraction_for_union_of_all_matches(const Molecule& m, const Substructure_Results& sresults)
{
  int matoms = m.natoms();

  int* atom_hit = new_int(matoms);
  std::unique_ptr<int[]> free_atom_hit(atom_hit);

  int nhits = sresults.number_embeddings();

  for (int i = 0; i < nhits; i++) {
    const Set_of_Atoms* e = sresults.embedding(i);

    e->set_vector(atom_hit, 1);
  }

  int n = count_non_zero_occurrences_in_array(atom_hit, matoms);

  return static_cast<float>(n) / static_cast<float>(matoms);
}

int
Needle::perform_search(Molecule& m, Molecule_to_Match& target)
{
  Substructure_Results sresults;

  int nhits = this->Substructure_Hit_Statistics::substructure_search(target, sresults);

  if (0 == nhits) {
    return 0;
  }

  // cerr << "Got " << nhits << " matches, " << sresults.number_embeddings() << "
  // embeddings\n";

  float fraction;

  if (compute_fraction_on_first_embedding_only || 1 == nhits) {
    const Set_of_Atoms* e = sresults.embedding(0);
    fraction = static_cast<float>(e->number_elements()) / static_cast<float>(m.natoms());
  } else {
    fraction = fraction_for_union_of_all_matches(m, sresults);
  }

  if (fraction < lower_fraction_matched) {
    matches_below_lower_fraction_matched++;
    return 1;  // hard to know what the return code should be
  } else if (fraction > upper_fraction_matched) {
    matches_above_upper_fraction_matched++;
    return 1;
  }

  Smiles_ID_and_Fraction* t = new Smiles_ID_and_Fraction(m.smiles(), m.name(), fraction);

  _neighbours.add(t);

  return 1;
}

class Fraction_Matched_Comparator
{
 private:
 public:
  int operator()(const Smiles_ID_and_Fraction*, const Smiles_ID_and_Fraction*) const;
};

int
Fraction_Matched_Comparator::operator()(const Smiles_ID_and_Fraction* q1,
                                        const Smiles_ID_and_Fraction* q2) const
{
  if (q1->fraction_atoms_matched() < q2->fraction_atoms_matched()) {
    return 1;
  }

  if (q1->fraction_atoms_matched() > q2->fraction_atoms_matched()) {
    return -1;
  }

  return 0;
}

static Fraction_Matched_Comparator fmc;

int
Needle::do_write(IWString_and_File_Descriptor& output)
{
  int n = _neighbours.number_elements();

  if (n > 0) {
    if (verbose) {
      highest_fraction.extra(1.0F - _neighbours[0]->fraction_atoms_matched());
    }
  } else if (!write_molecules_with_no_neighbours) {
    return 1;
  }

  // cerr << "Writing " << n << " nbrs\n";

  output << smiles_tag << _smiles << ">\n";
  output << identifier_tag << comment() << ">\n";

  if (n > 1) {
    _neighbours.iwqsort(fmc);
  }

  for (int i = 0; i < n; i++) {
    const Smiles_ID_and_Fraction* s = _neighbours[i];

    output << smiles_tag << s->smiles() << ">\n";
    output << identifier_tag << s->id() << ">\n";
    output << distance_tag;
    if (output_distances) {
      output << (1.0f - s->fraction_atoms_matched()) << ">\n";
    } else {
      output << s->fraction_atoms_matched() << ">\n";
    }
  }

  output << "|\n";

  return 1;
}

static int
substructure_match_fraction(Molecule& m, const resizable_array_p<Needle>& queries)
{
  Molecule_to_Match target(&m);

  int nq = queries.number_elements();

  int queries_hitting_this_molecule = 0;

  for (int i = 0; i < nq; i++) {
    //  cerr << "Compare " << m.name() << "' and '" << queries[i]->comment() << "'\n";
    if (allow_self_matches) {
      ;
    } else if (m.name() == queries[i]->comment()) {
      continue;
    }

    if (queries[i]->perform_search(m, target)) {
      queries_hitting_this_molecule++;
    }
  }

  queries_hitting[queries_hitting_this_molecule]++;

  return 1;
}

static int
substructure_match_fraction(data_source_and_type<Molecule>& input,
                            resizable_array_p<Needle>& queries)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!substructure_match_fraction(*m, queries)) {
      return 0;
    }
  }

  return 1;
}

static int
read_query_molecule(Molecule& m, resizable_array_p<Needle>& queries)
{
  preprocess(m);

  Needle* q = new Needle;

  if (!q->create_from_molecule(m, mqs)) {
    cerr << "Cannot create query from '" << m.name() << "'\n";
    return 0;
  }

  q->set_smiles(m.smiles());

  queries.add(q);

  return 1;
}

static int
read_query_molecules(data_source_and_type<Molecule>& input,
                     resizable_array_p<Needle>& queries)
{
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (!read_query_molecule(*m, queries)) {
      cerr << "Fatal error processing query molecule '" << m->name() << "'\n";
      return 0;
    }
  }

  return queries.number_elements();
}

static int
read_query_molecules(const const_IWSubstring& fname, resizable_array_p<Needle>& queries)
{
  FileType input_type = discern_file_type_from_name(fname);

  if (input_type == FILE_TYPE_INVALID) {
    cerr << "Cannot discern input type '" << fname << "'\n";
    return 0;
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  return read_query_molecules(input, queries);
}

static int
substructure_match_fraction(const char* fname, FileType input_type,
                            resizable_array_p<Needle>& queries)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return substructure_match_fraction(input, queries);
}

static int
substructure_match_fraction(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lq:fht:T:zr:a");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('r'))  // BEOFRE the -q option
  {
    int i = 0;
    const_IWSubstring r;
    while (cl.value('r', r, i++)) {
      if ('r' == r) {
        mqs.set_atoms_conserve_ring_membership(1);

        if (verbose) {
          cerr << "Query atoms will preserve their ring membership in matches\n";
        }
      } else if ('c' == r) {
        mqs.set_non_ring_atoms_become_nrings_0(1);

        if (verbose) {
          cerr << "Chain atoms will only match chain atoms in matches\n";
        }
      } else {
        cerr << "Unrecognised -r qualifier '" << r << "'\n";
        usage(22);
      }
    }
  }

  if (cl.option_present('a')) {
    compute_fraction_on_first_embedding_only = 0;

    if (verbose) {
      cerr << "All query matches will be used to determine atom coverage\n";
    }
  }

  resizable_array_p<Needle> queries;

  if (cl.option_present('q')) {
    int i = 0;
    const_IWSubstring q;

    while (cl.value('q', q, i++)) {
      if (q.starts_with("M:")) {
        q.remove_leading_chars(2);
      }

      if (!read_query_molecules(q, queries)) {
        cerr << "Cannot read query molecules from '" << q << "'\n";
        return 3;
      }
    }
  }

  int nq = queries.number_elements();

  if (0 == nq) {
    cerr << "Use the -q option to specify queries\n";
    usage(3);
  }

  if (verbose) {
    cerr << "Read " << nq << " query molecules\n";
  }

  for (int i = 0; i < nq; i++) {
    queries[i]->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('f')) {
    for (int i = 0; i < nq; i++) {
      queries[i]->set_max_matches_to_find(1);
    }

    if (verbose) {
      cerr << "Will find a maximum of one query match\n";
    }
  }

  if (cl.option_present('h')) {
    allow_self_matches = 0;

    if (verbose) {
      cerr << "Will avoid self matches\n";
    }
  }

  // Slight complication here because internally we operate on fraction atoms matched, but
  // to the external world we present distances

  if (cl.option_present('t')) {
    float t;
    if (!cl.value('t', t) || t < 0.0f || t > 1.0f) {
      cerr << "The lower distance cutoff option (-t) must be a valid distance\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will discard matches where the distance is below " << t << endl;
    }

    upper_fraction_matched = 1.0f - t;  // convert to fraction matched
  }

  if (cl.option_present('T')) {
    float t;
    if (!cl.value('T', t) || t < 0.0f || t > 1.0f) {
      cerr << "The upper distance cutoff option (-T) must be a valid distance\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will discard matches where the distance is above " << t << endl;
    }

    lower_fraction_matched = 1.0f - t;  // convert to fraction matched
  }

  if (lower_fraction_matched > upper_fraction_matched) {
    cerr << "Inconsistent lower (" << lower_fraction_matched << ") and upper ("
         << upper_fraction_matched << ") distance cutoffs\n";
    usage(4);
  }

  if (cl.option_present('z')) {
    write_molecules_with_no_neighbours = 1;

    if (verbose) {
      cerr << "Will write molecules with no neighbours\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!substructure_match_fraction(cl[i], input_type, queries)) {
      rc = i + 1;
      break;
    }
  }

  if (0 != rc) {
    return rc;
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < nq; i++) {
    queries[i]->do_write(output);

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < queries_hitting.number_elements(); i++) {
      if (queries_hitting[i]) {
        cerr << queries_hitting[i] << " molecules hit " << i << " queries\n";
      }
    }

    if (cl.option_present('T')) {
      float t;
      cl.value('T', t);
      cerr << matches_below_lower_fraction_matched
           << " molecules produced distances above " << t << endl;
    }

    if (cl.option_present('t')) {
      float t;
      cl.value('t', t);
      cerr << matches_above_upper_fraction_matched
           << " molecules produced distances below " << t << endl;
    }

    if (highest_fraction.n() > 1) {
      cerr << "Shortest distances between " << highest_fraction.minval() << " and "
           << highest_fraction.maxval() << " ave "
           << static_cast<float>(highest_fraction.average()) << endl;
    }
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = substructure_match_fraction(argc, argv);

  return rc;
}
