/*
  Do a series of transformations on molecules
*/

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

static int verbose = 0;

static uint64_t molecules_read = 0;

static uint64_t molecules_written = 0;

static uint64_t molecules_changed = 0;

static int enumerate_scaffold_hits_individually = 0;

// static int mmx = -1;

static int ignore_scaffold_multiple_substucture_matches = 0;

static int discard_molecules_not_reacting = 0;

static int ignore_molecules_not_reacting = 0;

static int ignore_sidechains_causing_increased_core_hit_count = 0;

static int write_molecules_not_reacting = 0;

static IWString append_to_changed_molecules;

static int append_reaction_names = 0;
static int append_reagent_names = 0;

static int write_original_molecule_to_output_stream = 0;

static int write_parent_of_changed_molecules = 0;

static int break_after_first_match = 0;

static int make_implicit_hydrogens_explicit = 0;

static int strip_products_to_largest_fragment = 0;

static int strip_input_molecules_to_largest_fragment = 0;

/*
  We can choose to do a specific embedding number - why?
*/

static int mdo = -1;

static int or_depth = 0;

static int or_matches_can_overlap = 0;

static Element_Transformations element_transformations;

/*
  We can also be part of a pipeline
*/

static IWString smiles_in, smiles_out;

static IWString identifier_tag("PCN<");

static Molecule_Output_Object output_stream;

static Molecule_Output_Object stream_for_molecules_not_reacting;

static Elements_to_Remove elements_to_remove;

static Chemical_Standardisation chemical_standardisation;

static int convert_isotopes = 0;

static int one_embedding_per_start_atom = 0;

/*
  When suppressing duplicates, we adopt a two-phase approach. First
  compare the smiles, then the unique smiles
*/

static int suppress_duplicate_molecules = 0;

static IW_STL_Hash_Set smiles_generated_current_molecule;
static IW_STL_Hash_Set unique_smiles_generated_current_molecule;

static int preserve_initial_aromaticity = 0;

static uint64_t molecules_with_recovered_aromaticity = 0;
static uint64_t molecules_rejected_for_aromaticity_loss = 0;

static IWString name_separator = " + ";
static IWString reaction_reagent_separator;

static int multiple_reaction_schemes = 0;

static int one_reaction_scheme_per_input_molecule = 0;

static int ok_atom_removal_reaction = 0;

static std::random_device rd;

static std::mt19937_64 rng(rd());

static double probability_sample_reactions = 1.0;

static uint64_t molecules_to_produce = 1;

static int unique_products_only = 0;

static IW_STL_Hash_Set seen;

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
  cerr << "  -R <rxn>       specify reaction file(s)\n";
  cerr << "  -T <fname>     read reaction and sidechain information from <fname>\n";
  cerr << "  -z i           ignore molecules not reacting, passed to subsequent reactions\n";
  cerr << "  -z d           discard molecules not reacting\n";
  cerr << "  -z w           write molecules not reacting\n";
  //cerr << "  -m <number>    the maximum number of scaffold reaction sites to process\n";
  cerr << "  -m do=number   process site number <number> in the scaffold\n";
  cerr << "  -m each        enumerate each scaffold hit separately together with the combined set\n";
  cerr << "  -m all         synomym for -m each\n";
  //cerr << "  -m each1       enumerate each scaffold hit separately\n";
  cerr << "  -m RMX         ignore any scaffold with multiple substructure matches\n";
  cerr << "  -r             when in -m each mode, ignore any further expansion when a sidechain addition adds a core matc (infinite loop) \n";

  cerr << "  -Z             ignore sidechains not reacting\n";
  cerr << "  -M each        generate all regio-isomers from multiple sidechain matches\n";
  cerr << "  -M all         synomym for -M each\n";
  cerr << "  -M do=number   process site number <number> in the sidechains\n";
//  cerr << "  -M mskip=text  append <text> to names where just one possible sidechain attachment chosen\n";
  cerr << "  -M write=file  write non-reacting sidechains to <file>\n";
  cerr << "  -M RMX         ignore any sidechains with multiple substructure matches\n";

  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
  cerr << "  -I             change isotopes to natural form in product molecules\n";
  cerr << "  -C <string>    append <string> to the name of all changed molecules\n";
  cerr << "  -C RXN         append the reaction name(s) of reactions that change the molecule\n";
  cerr << "  -j             write original molecule to output stream\n";
  cerr << "  -k             write original molecule to output stream only if changed\n";
  cerr << "  -F IN=<tag>    work as a filter, read smiles in <tag>\n";
  cerr << "  -F OUT=<tag>   work as a filter, write changed molecule smiles in <tag>\n";
  cerr << "  -b             break after first reaction matches\n";
  cerr << "  -h             make implicit Hydrogens explicit\n";
  cerr << "  -p <number>    perform transformations in parallel. Do all <number> combinations\n";
  cerr << "  -p ovok        embeddings can overlap - watch out, dangerous\n";
//cerr << "  -P <prob>      sample active reactions with probability <prob>\n";
  cerr << "  -d             suppress duplicate molecules - only checks current molecule\n";
  cerr << "  -u             one embedding per start atom\n";
  cerr << "  -L             strip products to largest fragment\n";
  cerr << "  -O ...         miscellaneous options, enter '-O help' for info\n";
  cerr << "  -W <string>    token put between names of products (default \" + \")\n";
//cerr << "  -s <col>       single reagent smiles is column <col> of compound name \n";
  cerr << "  -U <fname>     stream for molecules not reacting\n";
  cerr << "  -e             unique products only\n";
  cerr << "  -a rej         try to preserve initial aromaticity, reject molecules losing aromaticity\n";
  cerr << "  -a write       try to preserve initial aromaticity, write  molecules losing aromaticity\n";
  cerr << "  -t <ele1=ele2> standard element transformation options\n";
  cerr << "  -i <type>      specify input type\n";
  cerr << "  -o <type>      specify output type(s)\n";
  cerr << "  -S <filename>  output filename stem (should not include an extention - see \"-o\")\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  We need a small class to hold any AND or OR relationships within the set of reactions.
  They will be stored in the user specified void ptr of the reaction objects
*/

#define NO_GROUPING 0
#define AND_GROUPING 1
#define OR_GROUPING 2

class AndORProb {
 private:
  const int _group;
  const int _andor;
  const int _regioIsomer;
  double _probability;

  int _active;

 public:
  AndORProb(const int g, const int ao, const double p = 0.0, const int ri = 0)
      : _group(g),
        _andor(ao),
        _regioIsomer(ri),
        _probability(p),
        _active(1){
            //     _active = 1;
        };

  int
  andor() const {
    return _andor;
  }

  int
  group() const {
    return _group;
  }

  int
  regioIsomer() const {
    return _regioIsomer;
  }

  int
  active() const {
    return _active;
  }

  void
  set_active(const int s) {
    _active = s;
  }

  void
  set_probability(const double p) {
    _probability = p;
  }

  double
  probability() const {
    return _probability;
  }

  int operator==(const AndORProb& rhs) const;
};

int
AndORProb::operator==(const AndORProb& rhs) const {
  if (_group != rhs._group) {
    return 0;
  }

  return _andor == rhs._andor;
}

#define NO_OP_RXN "NOOP"

class Set_of_Reactions {
 private:
  resizable_array_p<IWReaction> _reaction;
  int _number_reactions;  // same as _reaction.number_elements()

  //  For sampling we can turn of various reactions

  int* _active;

  int* _is_noop;  // to avoid checking names

  double* _probability;
  std::uniform_real_distribution<double> _u;
  std::mt19937_64 _rng;

  Molecule _start_molecule;

  iwstring_data_source _input;

  int* _hits_to_query;

  //  private functions

  int _maybe_prepend_directory_name(IWString& file) const;
  int _reactions_have_starting_molecules();
  void _free_reaction_array();
  int _read_file_of_smirks(const const_IWSubstring& fname);
  int _read_file_of_smirks(iwstring_data_source& input);
  int _read_start_molecule(const IWString&);
  int _read_reaction(const IWString& buffer, Sidechain_Match_Conditions& smc, int& group,
                     const int andor, const double current_group_probability,
                     int& regioIsomer);
  int _read_file_of_reactions(const const_IWSubstring& fname,
                              Sidechain_Match_Conditions& smc);
  int _read_file_of_reactions(iwstring_data_source& input, const IWString& dirname,
                              Sidechain_Match_Conditions& smc);
  int _add_reaction(const_IWSubstring& fname, Sidechain_Match_Conditions& smc);
  int _add_reaction_from_smirks(const_IWSubstring& smirks,
                                Sidechain_Match_Conditions& smc);
  int _add_reaction_from_proto_file(const_IWSubstring& fname,
                                    Sidechain_Match_Conditions& smc);
  int _add_reaction_from_file_of_proto_reactions(iwstring_data_source& input,
                                                 const IWString& dirname,
                                                 Sidechain_Match_Conditions& smc);
  int _add_reaction_from_file_of_proto_reactions(const_IWSubstring& fname,
                                                 Sidechain_Match_Conditions& smc);
  int _build(const resizable_array_p<IWString>& zdata, Sidechain_Match_Conditions& smc);
  int prepOneSideChain(IWString& rxn_string,
                       std::vector<std::vector<IWString>>& smilesLists,
                       int sideChainIndex,
                       std::vector<Molecule_and_Embedding*>& currentMolAndEmbeddings,
                       Sidechain_Match_Conditions& smc, int& thisGroup,
                       const int thisAndor, const double current_group_probability,
                       int& regioIsomer);

 public:
  Set_of_Reactions();
  ~Set_of_Reactions();

  int debug_print(std::ostream& output) const;

  int open(const char* fname);
  int build(Sidechain_Match_Conditions& smc, int& fatal);
  int next(Sidechain_Match_Conditions& smc, int& fatal);
  int build(Command_Line& cl, const char opt, Sidechain_Match_Conditions& smc,
            const int verbose);

  int set_probability(const double);
  int set_probability_group(const int reaction_number);
  int reaction_is_active(const int rno);

  int is_start_of_and_group(const int reaction_number) const;
  int members_in_group(const int reaction_number) const;
  int lastRegioIsomerRxnId(const int reaction_number) const;

  int
  number_reactions() const {
    return _number_reactions;
  }

  int reactions_in_group(resizable_array<int>& rpg) const;

  int multiple_reaction_schemes_in_input();

  int identify_random_reactions(resizable_array<int>& do_rxn) const;

  void
  another_match_to_query(const int s) {
    _hits_to_query[s]++;
  }

  int report_matches(std::ostream& output) const;

  int
  has_starting_molecule() const {
    return _start_molecule.natoms();
  }

  int check_for_atom_removals(const int ok_atom_removal_reaction) const;

  void set_find_unique_embeddings_only(const int s);
  void set_max_matches_to_find(const int s);
  void set_find_one_embedding_per_start_atom(const int s);
  void set_append_names(const int s);

  int reactions_have_starting_molecules();

  int append_reagent_names(Molecule& m, const int* changed_by) const;
  int grouping_is_exhausted(const int reaction_number,
                            const AndORProb* current_andor) const;
  int lastOfRegioIsomers(const int reaction_number, const AndORProb* current_andor) const;
  void set_active_status(int reaction_number, const int a);

  Molecule&
  start_molecule() {
    return _start_molecule;
  }  // probably should be const

  IWReaction&
  reaction(const int i) {
    return *(_reaction[i]);
  }

  const IWReaction&
  reaction(const int i) const {
    return *_reaction[i];
  }

  resizable_array_p<IWReaction>&
  reactions() {
    return _reaction;
  }

  int shuffle_or_group(const int reaction_number);
};

Set_of_Reactions::Set_of_Reactions() : _u(0.0, 1.0), _rng(rd()) {
  _number_reactions = 0;

  _active = nullptr;

  _probability = nullptr;

  _is_noop = nullptr;

  _hits_to_query = nullptr;

  return;
}

Set_of_Reactions::~Set_of_Reactions() {
  _free_reaction_array();

  return;
}

void
Set_of_Reactions::_free_reaction_array() {
  for (int i = 0; i < _reaction.number_elements(); ++i) {
    void* v = _reaction[i]->user_specified_void_ptr();

    if (nullptr == v) {
      continue;
    }

    AndORProb* a = reinterpret_cast<AndORProb*>(v);

    delete a;
  }

  // cerr << "_reaction has " << _reaction.number_elements() << " items\n";
  _reaction.resize_keep_storage(0);
  _number_reactions = 0;

  if (nullptr != _hits_to_query) {
    delete[] _hits_to_query;
    _hits_to_query = nullptr;
  }

  if (nullptr != _probability) {
    delete[] _probability;
    _probability = nullptr;
  }

  if (nullptr != _active) {
    delete[] _active;
    _active = nullptr;
  }

  if (nullptr != _is_noop) {
    delete[] _is_noop;
    _is_noop = nullptr;
  }

  return;
}

int
Set_of_Reactions::debug_print(std::ostream& output) const {
  output << "Set_of_Reactions:debug_print: " << _number_reactions << " reactions\n";
  for (int i = 0; i < _number_reactions; ++i) {
    output << ' ' << i << ' ' << _reaction[i]->comment();

    const auto ao =
        reinterpret_cast<const AndORProb*>(_reaction[i]->user_specified_void_ptr());
    if (nullptr != ao) {
      output << " group " << ao->group() << " regioIsomer " << ao->regioIsomer()
             << " and/or " << ao->andor();
    }

    output << '\n';
  }

  return 1;
}

int
Set_of_Reactions::report_matches(std::ostream& output) const {
  if (nullptr == _hits_to_query) {
    return 0;
  }

  for (int i = 0; i < _number_reactions; i++) {
    output << _hits_to_query[i] << " hits to query " << i;
    const IWString& c = _reaction[i]->comment();

    if (c.length()) {
      output << " '" << c << "'";
    }
    output << '\n';
  }

  return 1;
}

int
Set_of_Reactions::multiple_reaction_schemes_in_input() {
  return _input.count_records_starting_with("$$$$");
}

int
Set_of_Reactions::reactions_have_starting_molecules() {
  const auto o = _input.tellg();

  const int rc = _reactions_have_starting_molecules();

  if (!_input.seekg(o)) {
    cerr << "Set_of_Reactions::reactions_have_starting_molecules:cannot seek back to "
         << o << '\n';
    return 0;
  }

  return rc;
}

int
Set_of_Reactions::_reactions_have_starting_molecules() {
  const_IWSubstring buffer;

  int found_start_molecule = 0;
  int foundNewRecord = 0;

  while (_input.next_record(buffer)) {
    if (buffer.starts_with("START:")) {
      found_start_molecule = 1;
    } else if ("$$$$" == buffer) {
      if (!found_start_molecule) {
        return 0;
      }

      found_start_molecule = 0;
      foundNewRecord = 0;
    } else {
      buffer.strip_leading_blanks();
      buffer.strip_trailing_blanks();
      foundNewRecord = 1;
    }
  }

  if (foundNewRecord && !found_start_molecule) {
    return 0;
  }

  return 1;
}

int
Set_of_Reactions::_maybe_prepend_directory_name(IWString& file) const {
  if (dash_s(file.null_terminated_chars())) {
    return 1;
  }

  IWString tmp(_input.fname());

  // cerr << tmp << '\n';

  int i = tmp.rindex('/');

  if (i < 0) {  // no directory component
    return 0;
  }

  tmp.iwtruncate(i);

  tmp << '/' << file;

  if (dash_s(tmp.null_terminated_chars())) {
    file = tmp;
    return 1;
  }

  return 0;
}

static int
random_number_between(const int zmin, const int zmax) {
  std::uniform_int_distribution<int> u(zmin, zmax);

  return u(rng);
}

int
Set_of_Reactions::reactions_in_group(resizable_array<int>& rpg) const {
  rpg.resize_keep_storage(0);

  for (int i = 0; i < _number_reactions; ++i) {
    const int g = members_in_group(i);
    rpg.add(g);
    i += g - 1;
  }

  return rpg.number_elements();
}

int
Set_of_Reactions::identify_random_reactions(resizable_array<int>& do_rxn) const {
  do_rxn.resize_keep_storage(0);

  for (int i = 0; i < _number_reactions; ++i) {
    const int g = members_in_group(i);
    if (1 == g) {
      do_rxn.add(i);
      continue;
    }

    const auto j = random_number_between(i, i + g - 1);
    do_rxn.add(j);
    i += g - 1;
  }

  return do_rxn.number_elements();
}

int
Set_of_Reactions::is_start_of_and_group(const int reaction_number) const {
  const AndORProb* curr_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number]->user_specified_void_ptr());  // may be null

  if (AND_GROUPING != curr_andor->andor()) {
    return 0;
  }

  if (0 == reaction_number) {  // must be the start
    return 1;
  }

  const AndORProb* prev_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number - 1]->user_specified_void_ptr());  // may be null

  if (curr_andor->group() == prev_andor->group()) {  // continuing already started group
    return 0;
  }

  return 1;
}

int
Set_of_Reactions::lastRegioIsomerRxnId(const int reaction_number) const {
  const AndORProb* curr_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number]->user_specified_void_ptr());  // may be null
  if (nullptr == curr_andor) {                                 // should not happen
    return reaction_number;
  }

  if (0 == curr_andor->andor() || 0 == curr_andor->regioIsomer()) {
    return reaction_number;
  }

  const int currentRegioIsomer = curr_andor->regioIsomer();

  for (int ndx = reaction_number; ndx < _number_reactions; ++ndx) {
    if (ndx + 1 >= _number_reactions) {
      return ndx;
    }

    const AndORProb* next_andor = reinterpret_cast<const AndORProb*>(
        _reaction[ndx + 1]->user_specified_void_ptr());  // may be null
    if (nullptr == next_andor ||
        next_andor->regioIsomer() != currentRegioIsomer) {  // should not happen
      return ndx;
    }
  }

  return _number_reactions;  // should never happen
}

int
Set_of_Reactions::members_in_group(const int reaction_number) const {
  const AndORProb* curr_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number]->user_specified_void_ptr());  // may be null
  if (nullptr == curr_andor) {                                 // should not happen
    return 1;
  }

  if (0 == curr_andor->andor()) {
    return 1;
  }

  const int current_group = curr_andor->group();

#ifdef DEBUG_MEMBERS_IN_GROUP
  cerr << "begin and group " << current_group << " with reaction " << reaction_number
       << '\n';
#endif

  int ndx = reaction_number +
            1;  // the index of the last member of the group, decrement at end of loop
  for (; ndx < _number_reactions; ++ndx) {
    const AndORProb* next_andor = reinterpret_cast<const AndORProb*>(
        _reaction[ndx]->user_specified_void_ptr());  // may be null
    if (nullptr == next_andor) {                     // should not happen
      break;
    }

#ifdef DEBUG_MEMBERS_IN_GROUP
    cerr << " reaction " << ndx << " group " << next_andor->group() << '\n';
#endif

    if (next_andor->group() != current_group) {
      break;
    }
  }

  return ndx - reaction_number;
}

/*
  We must return the number of items in the group

  Note that we do not check that we are in fact part of an AND group
*/

// #define DEBUG_SET_PROBABILITY_GROUP

int
Set_of_Reactions::set_probability_group(const int reaction_number) {
  const int gsize = members_in_group(reaction_number);
  if (1 == gsize) {
    return 1;
  }

#ifdef DEBUG_SET_PROBABILITY_GROUP
  cerr << " from " << reaction_number << " gsize " << gsize << '\n';
#endif

  int still_active = 0;

  const int gend = reaction_number + gsize - 1;

  for (int j = reaction_number; j < gend; ++j) {
    if (_is_noop[j]) {
      _active[j] = 1;
      continue;
    }

    const AndORProb* a =
        reinterpret_cast<const AndORProb*>(_reaction[j]->user_specified_void_ptr());

    if (1.0 == a->probability()) {
      _active[j] = 1;
      still_active++;
    } else if (_u(_rng) > _probability[j]) {
      _active[j] = 0;
    } else {
      _active[j] = 1;
      still_active++;
    }
  }

#ifdef DEBUG_SET_PROBABILITY_GROUP
  cerr << still_active << " still_active\n";
#endif

  if (0 == still_active) {
    std::uniform_int_distribution<int> mygroup(reaction_number, gend + 1);
    const int j = mygroup(_rng);
    _active[j] = 1;
  }

  return gsize;
}

int
Set_of_Reactions::set_probability(const double p) {
  if (nullptr == _probability) {
    _probability = new double[_number_reactions];
  }

  if (nullptr == _active) {
    _active = new int[_number_reactions];
  }

  if (nullptr == _is_noop) {
    _is_noop = new int[_number_reactions];
  }

  std::fill_n(_probability, _number_reactions, p);
  std::fill_n(_active, _number_reactions, 1);

  for (int i = 0; i < _number_reactions; ++i) {
    if (NO_OP_RXN == _reaction[i]->comment()) {
      _is_noop[i] = 1;
      _active[i] = 1;
    } else {
      _is_noop[i] = 0;
    }
  }

  // need to make sure we have at least one item from each AND group

  for (int i = 0; i < _number_reactions; ++i) {
    const AndORProb* curr_andor = reinterpret_cast<const AndORProb*>(
        _reaction[i]->user_specified_void_ptr());  // may be null
    if (nullptr == curr_andor) {                   // should not happen
      continue;
    }

    if (NO_GROUPING == curr_andor->andor()) {
      continue;
    }

    if (AND_GROUPING != curr_andor->andor()) {  // not sure what to do about OR groups
      continue;
    }

    const int gsize = set_probability_group(i);
    i += gsize - 1;  // don't forget the loop increment above
  }

  for (int i = 0; i < _number_reactions; ++i) {
    const AndORProb* ao =
        reinterpret_cast<const AndORProb*>(_reaction[i]->user_specified_void_ptr());
    if (nullptr == ao) {
      continue;
    }

    if (OR_GROUPING != ao->andor()) {
      continue;
    }

    const auto items_in_group = shuffle_or_group(i);
    i += items_in_group - 1;
  }

#ifdef DEBUG_SET_PROBABILITY
  cerr << "Set_of_Reactions::set_probability\n";

  for (int i = 0; i < _number_reactions; ++i) {
    cerr << ' ' << i << " active " << _active[i] << ' ' << _reaction[i]->comment()
         << '\n';
  }
#endif

  return 1;
}

int
Set_of_Reactions::reaction_is_active(const int rno) {
  if (nullptr == _active) {  // no probability in effect, will be done
    return 1;
  }

  return _active[rno];
}

static void
reset_duplicate_hash_sets() {
  smiles_generated_current_molecule.clear();
  unique_smiles_generated_current_molecule.clear();

  return;
}

static int
molecule_is_duplicate(Molecule& m) {
  if (unique_products_only) {  // compare across all molecules...
    if (seen.contains(m.unique_smiles())) {
      return 1;
    }

    seen.insert(m.unique_smiles());

    return 0;
  }

  if (smiles_generated_current_molecule.contains(m.smiles())) {
    return 1;
  }

  smiles_generated_current_molecule.insert(m.smiles());

  if (unique_smiles_generated_current_molecule.contains(m.unique_smiles())) {
    return 1;
  }

  unique_smiles_generated_current_molecule.insert(m.unique_smiles());

  return 0;
}

#ifdef NOT_BEING_USED
static int
do_append(IWString& new_name, const IWString& name_separator, const IWString& zextra) {
  if (new_name.length() > 0) {
    new_name << name_separator;
  }

  new_name << zextra;

  return zextra.length();
}
#endif

static int
do_append(Molecule& m, const IWString& name_separator, const IWString& zextra) {
  IWString tmp(m.name());
  tmp.append_with_spacer(zextra, name_separator);
  m.set_name(tmp);

  return 1;
}

/*
  There is ambiguity between this functionality and the property of the reactions
  themselves. The reactions will append their comment attribute, but we append
  the names of the substituents
*/

#ifdef NO_LONG_ER_NEEDEDDDD
int
Set_of_Reactions::append_reagent_names(Molecule& m, const int* changed_by) const {
  IWString new_name(m.name());
  cerr << "Appending reagent names " << _number_reactions << '\n';

  for (int i = 0; i < _number_reactions; ++i) {
    //  cerr << changed_by[i] << " changed_by " << i << '\n';

    if (!changed_by[i]) {
      continue;
    }

    if (NO_OP_RXN == _reaction[i]->comment()) {
      continue;
    }

    if (0 == _reaction[i]->number_sidechains()) {
      new_name.append_with_spacer(_reaction[i]->comment, name_separator);
      continue;
    }

    for (int j = 0; j < _reaction[i]->number_sidechains(); ++j) {
      const auto& s = _reaction[i]->sidechain(j);

      if (1 == s->number_reagents()) {
        new_name.append_with_spacer(s->reagent(0)->name(), name_separator);
      } else {
        new_name.append_with_spacer(_reaction[i]->comment(), name_separator);
      }
    }
  }

  m.set_name(new_name);

  return 1;
}
#endif

static void
do_append_to_name(Molecule& m) {
  IWString tmp = m.name();

  tmp.append_with_spacer(append_to_changed_molecules, name_separator);

  m.set_name(tmp);

  return;
}

void
Set_of_Reactions::set_find_unique_embeddings_only(const int s) {
  if (0 == _number_reactions) {
    return;
  }

  for (int i = 0; i < _number_reactions; ++i) {
    if (NO_OP_RXN != _reaction[i]->comment()) {
      _reaction[i]->set_find_unique_embeddings_only(s);
    }
  }

  return;
}

void
Set_of_Reactions::set_max_matches_to_find(const int s) {
  if (0 == _number_reactions) {
    return;
  }

  for (int i = 0; i < _number_reactions; ++i) {
    if (NO_OP_RXN != _reaction[i]->comment()) {
      printf("Setting max per reaction to %d\n", s);
      _reaction[i]->set_max_matches_to_find(s);
    }
  }

  return;
}

void
Set_of_Reactions::set_append_names(const int s) {
  for (int i = 0; i < _number_reactions; ++i) {
    _reaction[i]->set_append_names(s);
  }

  return;
}

void
Set_of_Reactions::set_find_one_embedding_per_start_atom(const int s) {
  // cerr << "Set_of_Reactions::set_find_one_embedding_per_start_atom:setting " << s << "
  // for " << _number_reactions << " reactions\n";
  for (int i = 0; i < _number_reactions; i++) {
    if (NO_OP_RXN == _reaction[i]->comment()) {
      continue;
    }

    _reaction[i]->set_one_embedding_per_start_atom(1);
    _reaction[i]->set_find_one_embedding_per_atom(1);

    //  cerr << "Reaction " << i << " has " << _reaction[i]->number_sidechains() << "
    //  sidechains\n";
    for (int j = 0; j < _reaction[i]->number_sidechains(); ++j) {
      _reaction[i]->sidechain(j)->set_find_one_embedding_per_atom(1);
    }
  }

  return;
}

int
Set_of_Reactions::check_for_atom_removals(const int ok_atom_removal_reaction) const {
  for (int i = 0; i < _number_reactions; i++) {
    if (!_reaction[i]->will_remove_atoms()) {
      continue;
    }

    cerr << "//t_of_Reactions::check_for_atom_removals:molecular_transformations may not "
            "be able to handle reactions that remove atoms, especially with multiple "
            "matches\n";

    if (ok_atom_removal_reaction) {
      continue;
    }

    cerr << "If you think this will be OK, use the undocumented -H option to give it a "
            "try\n";
    cerr << "Failures are likely to be catstrophic in the case of multiple substructure "
            "matches\n";
    return 0;
  }

  return 1;
}

static int
write_the_molecule(Molecule& m, std::ostream& output) {
  if (strip_products_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (suppress_duplicate_molecules && molecule_is_duplicate(m)) {
    return 1;
  }

  molecules_written++;

  if (output_stream.active()) {
    return output_stream.write(m);
  }

  if (smiles_out.length()) {
    output << smiles_out << m.smiles() << ">\n";

    if (0 == smiles_in.length()) {
      output << "|\n";
    }
  }

  return output.good();
}

static int
final_processing(Molecule& m, int changes_this_molecule, const Set_of_Reactions& sor,
                 const int* changed_by) {
  if (elements_to_remove.active()) {
    elements_to_remove.process(m);
  }

  if (convert_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  // cerr << "final_processing  changes_this_molecule " << changes_this_molecule << '\n';
  if (changes_this_molecule) {
    //  sor.append_reagent_names(m, changed_by);

    if (append_to_changed_molecules.length()) {
      do_append_to_name(m);
    }

    return write_the_molecule(m, std::cout);
  }

  if (write_molecules_not_reacting) {
    return write_the_molecule(m, std::cout);
  } else if (stream_for_molecules_not_reacting.is_open()) {
    stream_for_molecules_not_reacting.write(m);
  }

  return 1;
}

static int
do_write_original_molecule_to_output_stream(Molecule& m, std::ostream& output) {
  if (output_stream.active()) {
    return output_stream.write(m);
  }

  if (0 == smiles_out.length()) {  // hard to imagine this happening
    return 1;
  }

  if (smiles_in.length()) {
    output << smiles_in << m.smiles() << ">\n";
  } else {
    output << "$SMI<" << m.smiles() << ">\n";
  }

  output << identifier_tag << m.name() << ">\n";

  return output.good();
}

static int
do_name_appends(Molecule& m, const IWReaction& rxn) {
  if (append_reaction_names && append_reagent_names) {  // handled below
    ;
  } else if (append_reaction_names) {
    return do_append(m, name_separator, rxn.comment());
  } else if (append_reagent_names) {  // handled below
    ;
  } else {  // nothing to do
    return 1;
  }

  if (0 == rxn.number_sidechains()) {  // no reagents to worry about
    if (append_reaction_names) {
      do_append(m, name_separator, rxn.comment());
    }

    return 1;
  }

  for (int i = 0; i < rxn.number_sidechains(); ++i) {
    const auto& s = rxn.sidechain(i);
    if (1 == s->number_reagents()) {
      do_append(m, reaction_reagent_separator, s->reagent(0)->name());
    } else {
      //    do_append(m, name_separator, rxn.comment());    not sure what to do
    }
  }

  if (append_reaction_names) {
    do_append(m, name_separator, rxn.comment());
  }

  return 1;
}

static int
in_place_transformation_maybe_append_to_name(IWReaction& rxn, Molecule& m,
                                             const Set_of_Atoms* embedding) {
  if (!rxn.in_place_transformation(m, embedding)) {  // no reaction, nothing to do
    return 0;
  }

  return do_name_appends(m, rxn);
}

#ifdef NOT_BEING_USED
static int
in_place_transformations_maybe_append_to_name(IWReaction& rxn, Molecule& m,
                                              Substructure_Results& sresults) {
  if (!rxn.in_place_transformations(m, sresults)) {
    return 0;
  }

  return do_name_appends(m, rxn);
}
#endif

int
Set_of_Reactions::grouping_is_exhausted(const int reaction_number,
                                        const AndORProb* current_andor) const {
  if (reaction_number == (_number_reactions - 1)) {
    return 1;  // no more reactions to check.
  }

  const AndORProb* next_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number + 1]->user_specified_void_ptr());  // may be null

  if (nullptr == next_andor) {  // next reaction is not part of group, exhausted
    return 1;
  }

  if (next_andor->group() !=
      current_andor->group()) {  // different groups, current group is exhausted
    return 1;
  }

  // if they are not equal, then grouping is exhausted
  return current_andor->andor() != next_andor->andor();
}

int
Set_of_Reactions::lastOfRegioIsomers(const int reaction_number,
                                     const AndORProb* current_andor) const {
  if (reaction_number == (_number_reactions - 1)) {
    return 1;  // no more reactions to check.
  }

  const AndORProb* next_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number + 1]->user_specified_void_ptr());  // may be null

  if (nullptr == next_andor) {  // next reaction is not part of group, exhausted
    return 1;
  }

  if (next_andor->group() !=
      current_andor->group()) {  // different groups, current regioIsomer is exhausted
    return 1;
  }

  if (next_andor->regioIsomer() !=
      current_andor->regioIsomer()) {  // different regioIsomers, current one is exhausted
    return 1;
  }

  // if they are not equal, then grouping is exhausted
  return current_andor->andor() != next_andor->andor();
}

/*
  It is important that we return the number of items in the group
*/

int
Set_of_Reactions::shuffle_or_group(const int reaction_number) {
  const AndORProb* curr_andor = reinterpret_cast<const AndORProb*>(
      _reaction[reaction_number]->user_specified_void_ptr());  // may be null
  if (nullptr == curr_andor) {
    return 0;
  }

  if (OR_GROUPING != curr_andor->andor()) {
    return 0;
  }

  int ndx = reaction_number;

  for (; ndx < _number_reactions; ++ndx) {
    const AndORProb* next_andor = reinterpret_cast<const AndORProb*>(
        _reaction[ndx + 1]->user_specified_void_ptr());  // may be null
    if (nullptr == next_andor) {
      break;
    }

    if (next_andor->group() != curr_andor->group()) {
      break;
    }

    if (next_andor->andor() !=
        curr_andor->andor()) {  // not really necesssary as long as groups match
      break;
    }
  }

  if (ndx == reaction_number) {  // an OR group with just one member???
    return 1;
  }

  auto r = _reaction.rawdata();

  std::shuffle(r + reaction_number, r + ndx, _rng);

  return ndx - reaction_number + 1;
}

void
Set_of_Reactions::set_active_status(int reaction_number, const int a) {
  const AndORProb* ao =
      reinterpret_cast<AndORProb*>(_reaction[reaction_number]->user_specified_void_ptr());

  const int group = ao->group();

  for (; reaction_number < _number_reactions; ++reaction_number) {
    AndORProb* ao = reinterpret_cast<AndORProb*>(
        _reaction[reaction_number]->user_specified_void_ptr());

    if (ao->group() != group) {
      return;
    }

    ao->set_active(a);
  }

  return;
}

static uint64_t
compute_hash(const resizable_array<int>& do_rxn, const int number_reactions) {
  uint64_t rc = do_rxn[0];

  const int n = do_rxn.number_elements();

  for (int i = 1; i < n; ++i) {
    rc = rc * number_reactions + do_rxn[i];
  }

  return rc;
}

/*
  Made this a class just to keep the number of arguments being passed around controllable.
*/

class Molecular_Transformations_Random_Set_Generator {
 private:
  Set_of_Reactions& _sor;

  const int _number_reactions;

  //  We detect dups based on a hash of reaction numbers

  std::unordered_set<uint64_t> _seen;

  //  How many times do we randomly generate the same random reactions

  int _duplicate_items_skipped_top_level;

  //  How many times do we only use a subset of reactions and end up with a dup

  int _duplicate_items_skipped_final;

  //  A global array of how many times a reaction changes something

  int* _changed_by;

  //  private functions

  int _already_seen(const resizable_array<int>& do_rxn);

  template <typename O>
  int _final_processing(Molecule& m, const int changes_this_molecule, O& output);
  template <typename O>
  int _process(Molecule& m, int& changes_this_molecule, const int ndx,
               const resizable_array<int>& do_rxn, resizable_array<int>& actually_made,
               O& output);

 public:
  Molecular_Transformations_Random_Set_Generator(Set_of_Reactions& s, int* c);

  template <typename T>
  int process(Molecule& m, int& changes_this_molecule, T& output);
};

Molecular_Transformations_Random_Set_Generator::
    Molecular_Transformations_Random_Set_Generator(Set_of_Reactions& s, int* c)
    : _sor(s), _number_reactions(s.number_reactions()), _changed_by(c) {
  _duplicate_items_skipped_top_level = 0;
  _duplicate_items_skipped_final = 0;

  return;
}

int
Molecular_Transformations_Random_Set_Generator::_already_seen(
    const resizable_array<int>& do_rxn) {
  const auto h = compute_hash(do_rxn, _number_reactions);

  const auto f = _seen.find(h);

#ifdef DEBUG_ALREADY_SEEN
  cerr << "Q";
  for (int i = 0; i < do_rxn.number_elements(); ++i) {
    cerr << '_' << do_rxn[i];
  }
  cerr << " N " << do_rxn.size() << " hash value " << h << " seen " << (f != seen.end())
       << '\n';
#endif

  if (f != _seen.end()) {  // seen before, is dup
    _duplicate_items_skipped_final++;
    return 1;
  }

  _seen.insert(h);

  return 0;  // not seen before
}

template <typename O>
int
Molecular_Transformations_Random_Set_Generator::process(Molecule& m,
                                                        int& changes_this_molecule,
                                                        O& output) {
  changes_this_molecule = 0;

  resizable_array<int> do_rxn, actually_made;

  uint64_t i;
  for (i = 0; i < molecules_to_produce * 10; ++i) {
    _sor.identify_random_reactions(do_rxn);

    if (_already_seen(do_rxn)) {
      _duplicate_items_skipped_top_level++;
      //    do_write(do_rxn, cerr);
      continue;
    }

    actually_made.resize_keep_storage(0);

    if (!_process(m, changes_this_molecule, 0, do_rxn, actually_made, output)) {
      continue;
    }

    if (molecules_written == molecules_to_produce) {
      break;
    }
  }

  if (verbose) {
    cerr << "Ran " << i << " iterations to produce " << molecules_written
         << " random molecules. Skipped " << _duplicate_items_skipped_top_level
         << " duplicates at top level, " << _duplicate_items_skipped_final
         << " subsequent dups\n";
  }

  return 1;
}

/*
  If we did not use all of the selected reactions, we need to re-check for uniqueness at
  the end. 10 20 30 40 and 10 21 30 40

  are duplicates if the second reaction did not happen
*/

template <typename O>
int
Molecular_Transformations_Random_Set_Generator::_process(
    Molecule& m, int& changes_this_molecule, const int ndx,
    const resizable_array<int>& do_rxn, resizable_array<int>& actually_made, O& output) {
  if (ndx >= do_rxn.number_elements()) {
// #define SHOW_FINAL_MOLECULE
#ifdef SHOW_FINAL_MOLECULE
    cerr << "final_processing " << m.name() << ' ' << actually_made.number_elements()
         << " reactions";
    for (int i = 0; i < actually_made.number_elements(); ++i) {
      cerr << ' ' << actually_made[i];
    }
    cerr << " hash " << compute_hash(actually_made, sor.number_reactions()) << '\n';
    cerr << "actually_made cotains " << actually_made.size() << " items\n";
#endif

    if (actually_made.number_elements() < do_rxn.number_elements() &&
        _already_seen(actually_made)) {  // check partial set
#ifdef SHOW_FINAL_MOLECULE
      cerr << "DISC DUP";
      for (auto x : actually_made) {
        cerr << ' ' << x;
      }
      cerr << '\n';
#endif
      return 1;
    }

    return _final_processing(m, changes_this_molecule, output);
  }

  const int reaction_number = do_rxn[ndx];

  IWReaction& ri = _sor.reaction(reaction_number);

  if (NO_OP_RXN == ri.comment()) {
    _changed_by[reaction_number] += 1;
    _sor.another_match_to_query(reaction_number);
    actually_made.add(reaction_number);
    return _process(m, changes_this_molecule, ndx + 1, do_rxn, actually_made, output);
  }

  Substructure_Results sresults;

  const int nhits = ri.determine_matched_atoms(m, sresults);

  if (0 == nhits) {
    return _process(m, changes_this_molecule, ndx + 1, do_rxn, actually_made, output);
  }

  actually_made.add(reaction_number);
  if (verbose > 2) {
    cerr << nhits << " hits for reaction " << reaction_number << " level " << ndx << '\n';
  }

  _sor.another_match_to_query(reaction_number);
  changes_this_molecule++;
  _changed_by[reaction_number] = 1;

  if (1 == changes_this_molecule && write_parent_of_changed_molecules) {
    do_write_original_molecule_to_output_stream(m, output);
  }

  int hdo = 0;
  if (nhits > 1) {
    hdo = random_number_between(0, nhits - 1);
  }

  Molecule mcopy(m);

  in_place_transformation_maybe_append_to_name(ri, mcopy, sresults.embedding(hdo));
  return _process(mcopy, changes_this_molecule, ndx + 1, do_rxn, actually_made, output);
}

template <typename O>
int
Molecular_Transformations_Random_Set_Generator::_final_processing(
    Molecule& m, int changes_this_molecule, O& output) {
  if (elements_to_remove.active()) {
    elements_to_remove.process(m);
  }

  if (convert_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  // cerr << "final_processing  changes_this_molecule " << changes_this_molecule << '\n';
  if (changes_this_molecule) {
    //  sor.append_reagent_names(m, changed_by);

    if (append_to_changed_molecules.length()) {
      do_append_to_name(m);
    }

    return write_the_molecule(m, output);
  }

  if (write_molecules_not_reacting) {
    return write_the_molecule(m, output);
  } else if (stream_for_molecules_not_reacting.is_open()) {
    stream_for_molecules_not_reacting.write(m);
  }

  return 1;
}

static int
molecular_transformations_random_set(Molecule& m, int& changes_this_molecule,
                                     int* changed_by, Set_of_Reactions& sor) {
  if (verbose) {
    resizable_array<int> rpg;
    sor.reactions_in_group(rpg);

    for (int i = 0; i < rpg.number_elements(); ++i) {
      cerr << rpg[i] << " reactions in group " << i << '\n';
    }
  }

  Molecular_Transformations_Random_Set_Generator mtrsg(sor, changed_by);
  if (!mtrsg.process(m, changes_this_molecule, std::cout)) {
    return 0;
  }

  return 1;
}

static int molecular_transformations(Molecule& m, int& changes_this_molecule,
                                     int* changed_by, Set_of_Reactions& sor,
                                     int reaction_number);

static int
_doAllRegioIsomersAndSitesRecurs(Molecule& m, int& changes_this_molecule, int* changed_by,
                                 Set_of_Reactions& sor, int firstReactionNumber,
                                 int lastReactionNumber, Substructure_Results& sresults,
                                 int firstHit, int lastHit,
                                 std::set<IWString>* smilesAlreadyDone, int level) {
  // loop over the the hits to process

  int hitsToDo = lastHit - firstHit + 1;
  // cerr << "In _doAllRegioIsomersAndSitesRecurs do hits btw " << firstHit << " and " <<
  // lastHit << " hitsToDo " << hitsToDo << '\n';

  for (int thisHit = firstHit; thisHit <= lastHit; ++thisHit) {
    // loop over the  sidechain regioIsomers to be done (each one is a reaction)

    for (int reactionNumber = firstReactionNumber; reactionNumber <= lastReactionNumber;
         ++reactionNumber) {
      if (verbose > 1) {
        cerr << "level: " << level << '\n';
        cerr << " thisHit=" << thisHit << " of " << firstHit << " to " << lastHit << '\n';
        cerr << "reactionNumber=" << reactionNumber << " of " << firstReactionNumber
             << " to " << lastReactionNumber << '\n';
      }

      IWReaction& ri = sor.reaction(reactionNumber);
      // const AndORProb * ao = reinterpret_cast<const AndORProb
      // *>(ri.user_specified_void_ptr());    // may be null

      // do the reaction

      Molecule mcopy(m);

      in_place_transformation_maybe_append_to_name(ri, mcopy,
                                                   sresults.embedding(thisHit));
      int hitsLeftToDo = hitsToDo - 1;

      Substructure_Results newSresults;
      if (hitsLeftToDo) {  // we think we have more to do
        // we must re-do the search because the transformation above will have
        // (potentiallly) changed the atom mappings, and could have removed or added
        // matches to the query.  If it added hits, we bail out - this could mean that the
        // sidechain had a match to the core part of the reaction, and this would/could be
        // an infinate loop if we tried to react the new matches

        const int nhits = ri.determine_matched_atoms(mcopy, newSresults);

        // nHits should be one (or more) less than the hitsToDo passed in

        if (nhits > hitsLeftToDo) {
          cerr << "molecular_transformations::_doAllRegioIsomersAndSites, number of hits "
                  "has INCREASED.\n";
          cerr << "The side chain might have added a core match. " << reactionNumber
               << " '" << ri.comment() << '\n'
               << "  molecule was: " << mcopy.smiles() << '\n'
               << "  previous Molecule was: " << m.smiles() << '\n'
               << " expected " << hitsLeftToDo << " left to do, but found " << nhits
               << '\n';

          if (ignore_sidechains_causing_increased_core_hit_count) {
            return 1;
          }

          return 0;
        }
        hitsLeftToDo = nhits;  // should be at least one less hit
      }

      // See if we have already seen this one - if so do not repeat it

      Molecule umcopy(
          mcopy);  // so the uniqing does not alter the one we are really working with
      const IWString& uniqSmiles = umcopy.unique_smiles();
      if (smilesAlreadyDone->find(uniqSmiles) != smilesAlreadyDone->end()) {
        if (verbose > 1) {
          cerr << "Skipping a duplicate: " << uniqSmiles << '\n';
        }
      } else {
        smilesAlreadyDone->insert(uniqSmiles);

        if (enumerate_scaffold_hits_individually || hitsLeftToDo == 0) {
          // if doing all intermediate conversions , or this is the final level for this
          // hit/regio set, proceed with this one

          //        smilesAlreadyDone->insert(uniqSmiles);  IAW done above
          sor.set_active_status(reactionNumber, 0);
          if (!molecular_transformations(mcopy, changes_this_molecule, changed_by, sor,
                                         reactionNumber + 1)) {
            return 0;
          }
          sor.set_active_status(reactionNumber, 1);
        }

        // recursively call this routine to do the next mapping on the core

        if (hitsLeftToDo) {
          if (!_doAllRegioIsomersAndSitesRecurs(
                  mcopy, changes_this_molecule, changed_by, sor, firstReactionNumber,
                  lastReactionNumber, newSresults, 0, hitsLeftToDo - 1, smilesAlreadyDone,
                  level + 1)) {
            return 0;
          }
        }
      }
    }
  }

  return 1;
}

static int
_doAllRegioIsomersAndSites(Molecule& m, int& changes_this_molecule, int* changed_by,
                           Set_of_Reactions& sor, int firstReactionNumber,
                           int lastReactionNumber, Substructure_Results& sresults,
                           int firstHit, int lastHit) {
  std::set<IWString>
      smilesAlreadyDone;  // this is local so it will go away when this routine
                          // terminates.  It MAY not be used, if one was passed in

  return _doAllRegioIsomersAndSitesRecurs(
      m, changes_this_molecule, changed_by, sor, firstReactionNumber, lastReactionNumber,
      sresults, firstHit, lastHit, &smilesAlreadyDone, 1);
}

static int
molecular_transformations(Molecule& m, int& changes_this_molecule, int* changed_by,
                          Set_of_Reactions& sor, int reaction_number) {
// #define DEBUG_MOLECULAR_TRANDFORMATIONS
#ifdef DEBUG_MOLECULAR_TRANDFORMATIONS
  cerr << "molecular_transformations processing reaction number " << reaction_number
       << " of " << sor.number_reactions() << " start " << m.smiles() << ' ' << m.name()
       << '\n';
#endif

  if (reaction_number >= sor.number_reactions()) {
    return final_processing(m, changes_this_molecule, sor, changed_by);
  }

  if (break_after_first_match && changes_this_molecule) {
    return final_processing(m, 1, sor, changed_by);
  }

  if (sor.is_start_of_and_group(reaction_number) && probability_sample_reactions < 1.0) {
    sor.set_probability_group(reaction_number);
  }

  if (!sor.reaction_is_active(reaction_number)) {
    return molecular_transformations(m, changes_this_molecule, changed_by, sor,
                                     reaction_number + 1);
  }

  IWReaction& ri = sor.reaction(reaction_number);

  const AndORProb* ao =
      reinterpret_cast<const AndORProb*>(ri.user_specified_void_ptr());  // may be null

#ifdef DEBUG_MOLECULAR_TRANDFORMATIONS
  cerr << "Begin reaction " << reaction_number << " ao active " << ao->active() << '\n';
#endif

  if (!ao->active()) {
    return molecular_transformations(m, changes_this_molecule, changed_by, sor,
                                     reaction_number + 1);
  }

  const int andor = ao->andor();

  if (NO_OP_RXN == ri.comment()) {
    changed_by[reaction_number] += 1;
    sor.another_match_to_query(reaction_number);
    return molecular_transformations(m, changes_this_molecule, changed_by, sor,
                                     reaction_number + 1);
  }

  Substructure_Results sresults;

  const int nhits = ri.determine_matched_atoms(m, sresults);

#ifdef DEBUG_MOLECULAR_TRANDFORMATIONS
  cerr << nhits << " searching reaction " << reaction_number << " in " << m.smiles()
       << ' ' << m.name() << '\n';
#endif

  if (0 == nhits) {
    //  cerr << "zero hits, andor " << andor << " rxn " << reaction_number << " in " <<
    //  m.smiles() << ' ' << m.name() << '\n';
    if (OR_GROUPING == andor) {
      if (!sor.grouping_is_exhausted(reaction_number, ao)) {
        return molecular_transformations(
            m, changes_this_molecule, changed_by, sor,
            reaction_number + 1);  // try next member of OR group
      }
    }

    if (discard_molecules_not_reacting) {
      return 1;
    }

    if (!ignore_molecules_not_reacting) {
      cerr << "molecular_transformations::yipes, zero hits for reaction "
           << reaction_number << " '" << ri.comment() << "', molecule " << m.smiles()
           << ' ' << m.name() << "'\n";
      return 0;
    }

    if (verbose > 2) {
      cerr << m.smiles() << ' ' << m.name() << " no hits to reaction " << reaction_number
           << '\n';
    }

    return molecular_transformations(m, changes_this_molecule, changed_by, sor,
                                     reaction_number + 1);
  }

  if (verbose > 2) {
    cerr << nhits << " hits for reaction " << reaction_number << '\n';
  }

  if (ignore_scaffold_multiple_substucture_matches && nhits > 1) {
    if (0 == verbose) {
      cerr << "'" << m.name() << "' "
           << " " << nhits << " hits in scaffold, ignored" << '\n';
    }
    return 1;
  }

  sor.another_match_to_query(reaction_number);
  changes_this_molecule++;
  changed_by[reaction_number] = 1;

  if (1 == changes_this_molecule && write_parent_of_changed_molecules) {
    do_write_original_molecule_to_output_stream(m, std::cout);
  }

  // Can we process a specific hit number?
  int hit_number = -1;

  if (1 == nhits) {
    hit_number = 0;
  } else if (mdo >= 0) {
    if (static_cast<uint32_t>(mdo) >= sresults.number_embeddings()) {
      cerr << "Reaction '" << ri.comment() << "' only " << nhits << " hits, so embedding "
           << mdo << " invalid\n";
      return 0;
    }

    hit_number = mdo;
  }

  // see how many sidechain regioisomres there are come from the same sidechain
  // There are actually reactions in with the same regioIsomer id

  int last_reaction_number = sor.lastRegioIsomerRxnId(reaction_number);

  // now the recursive calls to generate all reactions at the sites specified (one or all)
  // and the regioIsomers (could be just one!)

  if (hit_number >= 0) {
    if (!_doAllRegioIsomersAndSites(m, changes_this_molecule, changed_by, sor,
                                    reaction_number, last_reaction_number, sresults,
                                    hit_number, hit_number)) {
      return 0;
    }
  } else {
    if (!_doAllRegioIsomersAndSites(m, changes_this_molecule, changed_by, sor,
                                    reaction_number, last_reaction_number, sresults, 0,
                                    nhits - 1)) {
      return 0;
    }
  }

  if (AND_GROUPING == andor) {
    if (sor.grouping_is_exhausted(last_reaction_number, ao)) {
      return 1;
    }

    return molecular_transformations(m, changes_this_molecule, changed_by, sor,
                                     last_reaction_number + 1);
  } else if (OR_GROUPING == andor) {
    return 1;
  } else {
    return 1;
  }
}

/*
  We need to pass a lot of info from one level to the next. Easier to put it in a class
*/

class MT_OR_Data {
 private:
  int _number_reactions;
  int _matoms;
  Substructure_Results* _sresults;
  int* _changed_by;
  int* _atoms_hit;
  int* _initial_aromaticity;
  int _depth;
  int _istart;
  int _jstart;

  //  We need the ability to push our state

  resizable_array<int> _depth_stack, _istart_stack, _jstart_stack;

 public:
  MT_OR_Data(int nr, int na);
  ~MT_OR_Data();

  int identify_matched_atoms(Set_of_Reactions&, Molecule& m);

  int discern_initial_aromaticity(Molecule& m);

  Substructure_Results*
  sresults() {
    return _sresults;
  }

  int*
  changed_by() {
    return _changed_by;
  }

  int
  depth() const {
    return _depth;
  }

  int
  istart() const {
    return _istart;
  }

  int
  jstart() const {
    return _jstart;
  }

  void
  set_depth(int s) {
    _depth = s;
  }

  void
  set_istart(int s) {
    _istart = s;
  }

  void
  set_jstart(int s) {
    _jstart = s;
  }

  void
  change_depth(int s) {
    _depth += s;
  }

  void set_atoms_hit(const int s);
  void set_atoms_hit(const Set_of_Atoms& e, const int s);

  int overlapping_embeddings(const Set_of_Atoms& e) const;

  void push_state();
  int pop_state();

  int next_state();

  int aromaticity_preserved(Molecule& m);
};

MT_OR_Data::MT_OR_Data(int nr, int na) {
  _number_reactions = nr;
  _matoms = na;

  _sresults = new Substructure_Results[nr];
  _changed_by = new_int(nr);

  if (or_matches_can_overlap) {
    _atoms_hit = new int[na];
  } else {
    _atoms_hit = nullptr;
  }

  if (preserve_initial_aromaticity) {
    _initial_aromaticity = new int[na];
  } else {
    _initial_aromaticity = nullptr;
  }

  _depth = 1;
  _istart = 0;
  _jstart = 0;

  return;
}

MT_OR_Data::~MT_OR_Data() {
  if (nullptr != _atoms_hit) {
    delete[] _atoms_hit;
  }
  if (nullptr != _initial_aromaticity) {
    delete[] _initial_aromaticity;
  }

  return;
}

void
MT_OR_Data::set_atoms_hit(const int s) {
  if (nullptr == _atoms_hit) {
    return;
  }

  set_vector(_atoms_hit, _matoms, s);

  return;
}

void
MT_OR_Data::set_atoms_hit(const Set_of_Atoms& e, const int s) {
  if (nullptr == _atoms_hit) {
    return;
  }

  e.set_vector(_atoms_hit, s);

  return;
}

int
MT_OR_Data::identify_matched_atoms(Set_of_Reactions& sor, Molecule& m) {
  resizable_array_p<IWReaction>& reaction = sor.reactions();

  int changes_this_molecule = 0;

  for (int i = 0; i < _number_reactions; i++) {
    const auto nhits = reaction[i]->determine_matched_atoms(m, _sresults[i]);

    if (0 == nhits) {
      if (discard_molecules_not_reacting) {
        continue;
      }

      if (!ignore_molecules_not_reacting) {
        cerr << "Yipes, zero hits for reaction '" << reaction[i]->comment()
             << "', molecule " << m.name() << "'\n";
        return 0;
      }

      continue;
    }

    if (verbose > 2) {
      cerr << nhits << " hits for reaction " << i << " '" << reaction[i]->comment()
           << "'\n";
    }

    _changed_by[i] = 1;

    sor.another_match_to_query(i);

    changes_this_molecule++;
  }

  return changes_this_molecule;
}

int
MT_OR_Data::discern_initial_aromaticity(Molecule& m) {
  assert(nullptr != _initial_aromaticity);

  return m.aromaticity(_initial_aromaticity);
}

int
MT_OR_Data::overlapping_embeddings(const Set_of_Atoms& e) const {
  if (nullptr == _atoms_hit) {
    return 0;
  }

  return e.any_members_set_in_array(_atoms_hit);
}

void
MT_OR_Data::push_state() {
  _depth_stack.add(_depth);
  _istart_stack.add(_istart);
  _jstart_stack.add(_jstart);

  return;
}

int
MT_OR_Data::pop_state() {
  if (0 == _depth_stack.size()) {
    cerr << "MT_OR_Data::pop_state:no state\n";
    return 0;
  }

  _depth = _depth_stack.pop();
  _istart = _istart_stack.pop();
  _jstart = _jstart_stack.pop();

  return 1;
}

/*
  Can we move to the next embedding within the current set of matches, or the next
  reaction
*/

int
MT_OR_Data::next_state() {
  const int nhits = _sresults[_istart].number_embeddings();

  if (_jstart < nhits - 1) {  // just move to the next embedding for the current reaction
    push_state();
    _jstart++;
    _depth++;
    return 1;
  }

  // Need to move to the next reaction

  if (_istart == _number_reactions - 1) {  // these are done
    return 0;
  }

  push_state();
  _istart++;
  _depth++;
  _jstart = 0;

  return 1;
}

int
MT_OR_Data::aromaticity_preserved(Molecule& m) {
  if (0 == preserve_initial_aromaticity) {  // we don't care
    return 1;
  }

  m.compute_aromaticity_if_needed();

  int aromaticity_mismatches = 0;

// #define DEBUG_AROMATICITY_PRESEVRATION
#ifdef DEBUG_AROMATICITY_PRESEVRATION
  cerr << "New molecule " << m.smiles() << '\n';
#endif

  for (int i = 0; i < _matoms; ++i) {
#ifdef DEBUG_AROMATICITY_PRESEVRATION
    if (_initial_aromaticity[i]) {
      cerr << "Initially atom " << i << " aromatic " << m.smarts_equivalent_for_atom(i)
           << " now " << m.is_aromatic(i) << '\n';
    }
#endif

    if (_initial_aromaticity[i] && m.is_aromatic(i)) {
      continue;
    }

    aromaticity_mismatches++;
  }

  if (0 == aromaticity_mismatches) {  // great
    return 1;
  }

  int* a = new int[_matoms];
  std::unique_ptr<int[]> free_a(a);

  std::copy_n(_initial_aromaticity, _matoms, a);

  const int nr = m.nrings();

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (ri->is_aromatic()) {  // ring is aromatic after the change, great
      continue;
    }

    if (!ri->all_members_set_in_array(_initial_aromaticity,
                                      1)) {  // ring was not aromatic to begin with
      continue;
    }

    ri->set_vector(a, 2);  // magic number

    if (ri->is_fused()) {
      for (int j = i + 1; j < nr; ++j) {
        const Ring* rj = m.ringi(j);
        if (ri->fused_system_identifier() != rj->fused_system_identifier()) {
          continue;
        }

        if (!rj->all_members_set_in_array(_initial_aromaticity, 1)) {
          continue;
        }

        rj->set_vector(a, 2);  // magic number
      }
    }
  }

  for (int i = 0; i < _matoms; ++i) {
    if (2 == a[i]) {  // magic number
      ;
    } else {
      a[i] = 0;
    }
  }

  for (int i = 0; i < _matoms; ++i) {
    if (a[i] && !m.is_ring_atom(i)) {
      cerr << "HUH, non ring atom " << i << " in " << m.smiles() << ' ' << m.name() << ' '
           << m.smarts_equivalent_for_atom(i) << '\n';
    }

    if (a[i] && m.ncon(i) > 3) {
      cerr << "HUH, 4 connected aromatic atom " << i << " in " << m.smiles() << ' '
           << m.name() << ' ' << m.smarts_equivalent_for_atom(i) << '\n';
    }
  }

  if (m.find_kekule_form(a)) {
    molecules_with_recovered_aromaticity++;
    return 1;
  }

  if (1 == preserve_initial_aromaticity) {
    return 1;
  }

  molecules_rejected_for_aromaticity_loss++;

  return 0;
}

static int
molecular_transformations_or_operation(Molecule& m, MT_OR_Data& mtord,
                                       Set_of_Reactions& sor) {
  if (1 == mtord.depth()) {
    mtord.set_atoms_hit(0);
  }

  const int number_reactions = sor.number_reactions();

  resizable_array_p<IWReaction>& reaction = sor.reactions();

  // cerr << "Starting recursion depth " << depth << " istart " << istart << " jstart " <<
  // jstart << '\n';

  Substructure_Results* sresults = mtord.sresults();
  int* changed_by = mtord.changed_by();

  for (int i = mtord.istart(); i < number_reactions; i++) {
    int nhits = sresults[i].number_embeddings();

    if (0 == nhits) {
      continue;
    }

    changed_by[i] = 1;

    for (int j = mtord.jstart(); j < nhits; j++) {
      const Set_of_Atoms* e = sresults[i].embedding(j);

      //    cerr << "Hit " << j << '\n';

      if (mtord.overlapping_embeddings(*e)) {
        continue;
      }

      Molecule mcopy(m);

      reaction[i]->in_place_transformation(mcopy, e);

      if (!mtord.aromaticity_preserved(mcopy)) {
        continue;
      }

      final_processing(mcopy, 1, sor, mtord.changed_by());

      //    cerr << " i = " << i << " j = " << j << " depth " << depth << '\n';

      if (or_depth == mtord.depth()) {
        continue;
      }

      if (!mtord.next_state()) {
        continue;
      }

      mtord.set_atoms_hit(*e, 1);

      molecular_transformations_or_operation(mcopy, mtord, sor);

      mtord.pop_state();
    }

    changed_by[i] = 0;
  }

  return 1;
}

static int
molecular_transformations_or_operation(Molecule& m, Set_of_Reactions& sor) {
  MT_OR_Data mtord(sor.number_reactions(), m.natoms());

  const int changes_this_molecule = mtord.identify_matched_atoms(sor, m);

  if (0 == changes_this_molecule) {
    if (write_molecules_not_reacting) {
      write_the_molecule(m, std::cout);
    }

    return 1;
  }

  if (preserve_initial_aromaticity) {
    mtord.discern_initial_aromaticity(m);
  }

  molecules_changed++;
  return molecular_transformations_or_operation(m, mtord, sor);
}

#ifdef MAYBE_NOT_SUCH_A_GOOD_IDEA
static int
molecular_transformations_or_operation_SS_each_stem(Molecule& m, int d) {
  Substructure_Results sresults;

  const auto nhits = reaction[i]->determine_matched_atoms(m, sresults);

  if (0 == nhits) {
    if (d < or_depth) {
      return molecular_transformations_or_operation_SS_each_stem(m, d + 1);
    } else {
      return 0;
    }
  }

  for (int i = 0; i < nhits; ++i) {
    const Set_of_Atoms* e = sresults.embedding(i);

    Molecule mcopy(m);
  }
}
#endif

static int
molecular_transformations(Molecule& m, Set_of_Reactions& sor) {
  int changes_this_molecule = 0;

  int* changed_by = new_int(sor.number_reactions());
  std::unique_ptr<int[]> free_changed_by(changed_by);

  int rc = 0;
  if (1 == molecules_to_produce) {
    rc = molecular_transformations(m, changes_this_molecule, changed_by, sor, 0);
  } else {
    rc = molecular_transformations_random_set(m, changes_this_molecule, changed_by, sor);
  }

  if (changes_this_molecule) {
    molecules_changed++;
  }

  return rc;
}

static int
molecular_transformations(Set_of_Reactions& sor, Sidechain_Match_Conditions& smc) {
  int reaction_schemes_read = 0;

  int fatal = 0;
  while (sor.next(smc, fatal)) {
    if (verbose > 1) {
      sor.debug_print(cerr);
    }

    Molecule& s = sor.start_molecule();
    if (0 == s.natoms()) {
      cerr << "molecular_transformations:no start molecule\n";
      return 0;
    }

    reaction_schemes_read++;

    //  cerr << "after NEXT: start molecule " << s.smiles() << ' ' << s.name() << '\n';
    if (!molecular_transformations(s, sor)) {
      return 0;
    }
  }

  if (0 == reaction_schemes_read) {
    cerr << "No reaction schemes fetched from input\n";
  }

  return reaction_schemes_read;
}

static void
preprocess(Molecule& m) {
  if (strip_input_molecules_to_largest_fragment) {
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
molecular_transformations(data_source_and_type<Molecule>& input, Set_of_Reactions& sor,
                          Sidechain_Match_Conditions& smc) {
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1) {
      cerr << molecules_read << " processing '" << m->name() << "'\n";
    }

    if (suppress_duplicate_molecules) {
      reset_duplicate_hash_sets();
    }

    preprocess(*m);

    if (one_reaction_scheme_per_input_molecule || 0 == sor.number_reactions()) {
      int fatal;
      if (!sor.next(smc, fatal)) {
        cerr << "molecular_transformations::cannot read next set of reactions\n";
        return 0;
      }

      sor.set_find_unique_embeddings_only(
          1);  // not sure if this is the right idea or not
      if (one_embedding_per_start_atom) {
        sor.set_find_one_embedding_per_start_atom(1);
      }
      if (append_reaction_names) {
        sor.set_append_names(1);
      }
    }

    if (write_original_molecule_to_output_stream) {
      do_write_original_molecule_to_output_stream(*m, std::cout);
    }

    int rc;

    if (or_depth) {
      rc = molecular_transformations_or_operation(*m, sor);
    } else {
      rc = molecular_transformations(*m, sor);
    }

    if (0 == rc) {
      return 0;
    }
  }

  return 1;
}

static int
molecular_transformations_tdt(const const_IWSubstring& buffer, Set_of_Reactions& sor) {
  assert(buffer.ends_with('>'));

  const_IWSubstring smi(buffer);
  smi.chop();
  smi.remove_up_to_first('<');

  Molecule m;
  if (!m.build_from_smiles(smi)) {
    cerr << "Cannot create molecule from smiles '" << smi << "'\n";
    return 0;
  }

  preprocess(m);

  return molecular_transformations(m, sor);
}

static int
molecular_transformations_tdt(iwstring_data_source& input, Set_of_Reactions& sor) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!buffer.starts_with(smiles_in)) {
      std::cout << buffer << '\n';
      continue;
    }

    molecules_read++;

    if (write_original_molecule_to_output_stream) {
      std::cout << buffer << '\n';
    }

    if (!molecular_transformations_tdt(buffer, sor)) {
      return 0;
    }
  }

  return std::cout.good();
}

static int
molecular_transformations(const char* fname, FileType input_type, Set_of_Reactions& sor,
                          Sidechain_Match_Conditions& smc) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  return molecular_transformations(input, sor, smc);
}

static int
report_results(const Set_of_Reactions& sor, std::ostream& os) {
  os << "Read " << molecules_read << " molecules, changed " << molecules_changed << '\n';

  sor.report_matches(os);

  cerr << "Wrote " << molecules_written << " molecules\n";

  if (preserve_initial_aromaticity) {
    cerr << molecules_with_recovered_aromaticity
         << " molecules that were able to recover otherwise lost aromaticity\n";
    cerr << molecules_rejected_for_aromaticity_loss
         << " molecules discarded for loss of aromaticity\n";
  }

  return os.good();
}

int
Set_of_Reactions::_read_file_of_smirks(iwstring_data_source& input) {
  const_IWSubstring buffer;

  int rc = 0;

  while (input.next_record(buffer)) {
    if (0 == buffer.length()) {
      continue;
    }

    if (buffer.starts_with('#')) {
      continue;
    }

    IWReaction* r = new IWReaction();

    if (!r->construct_from_smirks(buffer)) {
      cerr << "Invalid smirks '" << buffer << "'\n";
      delete r;
      return 0;
    }

    _reaction.add(r);
    rc++;
  }

  return rc;
}

int
Set_of_Reactions::_read_file_of_smirks(const const_IWSubstring& fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return _read_file_of_smirks(input);
}

int
Set_of_Reactions::_read_file_of_reactions(iwstring_data_source& input,
                                          const IWString& dirname,
                                          Sidechain_Match_Conditions& smc) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    buffer.strip_leading_blanks();
    buffer.strip_trailing_blanks();

    if (0 == buffer.length()) {
      continue;
    }

    if (buffer.starts_with('#')) {
      continue;
    }

    int rc;

    IWReaction* r = new IWReaction();
    set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(
        0);  // do not show the message

    if (buffer.starts_with('/')) {
      rc = r->do_read(buffer, smc);
    } else {
      IWString fname;
      fname << dirname << '/' << buffer;

      rc = r->do_read(fname.null_terminated_chars(), smc);
    }

    if (0 == rc) {
      cerr << "Cannot read '" << buffer << "'\n";
      delete r;
      return 0;
    }

    _reaction.add(r);
  }

  return 1;
}

// Return the directory name of `fname`.
template <typename T>
IWString
GetDirName(const T& fname) {
  IWString dirname(fname);

  int i = dirname.rindex(std::filesystem::path::preferred_separator);

  if (i < 0) {
    dirname = ".";
  } else {
    dirname.iwtruncate(i);
  }

  return dirname;
}

int
Set_of_Reactions::_read_file_of_reactions(const const_IWSubstring& fname,
                                          Sidechain_Match_Conditions& smc) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Set_of_Reactions::_read_file_of_reactions:cannot open '" << fname << "'\n";
    return 0;
  }

  IWString dirname(fname);

  int i = dirname.rindex('/');

  if (i < 0) {
    dirname = ".";
  } else {
    dirname.iwtruncate(i);
  }

  // cerr << "Reading reactions, dirname '" << dirname << "'\n";

  return _read_file_of_reactions(input, dirname, smc);
}

int
Set_of_Reactions::prepOneSideChain(
    IWString& rxn_string, std::vector<std::vector<IWString>>& smilesLists,
    int sideChainIndex, std::vector<Molecule_and_Embedding*>& currentMolAndEmbeddings,
    Sidechain_Match_Conditions& smc, int& thisGroup, const int thisAndor,
    const double current_group_probability, int& regioIsomer) {
  if (sideChainIndex >= static_cast<int>(smilesLists.size())) {
    // done with side chains - so add the new reaction with the chosen mappings

    // get a local copy of Sidechain_Match_Conditions - We manually handle the settings
    // for multiple embeddings, so we always send a particular embedding into each
    // reaction

    Sidechain_Match_Conditions smcLocal(smc);
    smcLocal.set_process_hit_number(0);  // only one embedding per query

    // make sure we have a new rxn

    IWReaction* r = new IWReaction();
    if (!r->do_read(rxn_string, smcLocal)) {
      cerr << "Set_of_Reactions::prepOneSideChain:cannot read '" << rxn_string << "'\n";
      delete r;
      return 0;
    }

    set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(
        0);  // do not show the message

    if (verbose > 1) {
      cerr << "adding new reaction - thisGroup: " << thisGroup
           << "  thisAndor: " << thisAndor << "   regioIsomer: " << regioIsomer << '\n';
    }
    r->set_user_specified_void_ptr(
        new AndORProb(thisGroup, thisAndor, current_group_probability, regioIsomer));
    _reaction.add(r);

    if (static_cast<int>(currentMolAndEmbeddings.size()) != r->number_sidechains()) {
      cerr << "Set_of_Reactions::prepOneSideChain:number of sidechains in rxn does not "
              "match directive from config file "
           << '\n';
      return 0;
    }

    r->setup_to_skip_multiple_embeddings_involving_non_changing_atoms();
    for (int sideIndex = 0; sideIndex < r->number_sidechains(); ++sideIndex) {
      if (r->sidechain(sideIndex)->single_reagent_only()) {
        continue;
      }

      // make a deep copy of the Molecule and embedding - the reaction thinks it owns this
      // and will try to delete it, so there must be a copy for each one

      // make a copy of the mol and embedding  - these are handed to the reaction

      Molecule* mCopy(currentMolAndEmbeddings[sideIndex]);
      Set_of_Atoms* eCopy =
          new Set_of_Atoms(*(currentMolAndEmbeddings[sideIndex]->embedding()));

      Molecule_and_Embedding* qCopy = new Molecule_and_Embedding();
      qCopy->add_molecule(mCopy);
      qCopy->set_name(mCopy->name());
      qCopy->set_embedding(*eCopy);

      if (!r->sidechain(sideIndex)->add_reagent(qCopy, smc)) {
        return 0;
      }
    }

    return 1;
  }

  std::vector<IWString>& smilesList = smilesLists[sideChainIndex];

  // still more sidechains to prep  if the next side chain does not take smiles (it has
  // its own smiles defined), just go to the next one

  if (smilesList.size() == 0) {
    currentMolAndEmbeddings[sideChainIndex] = nullptr;

    // process the next side chain recursively

    if (!prepOneSideChain(rxn_string, smilesLists, sideChainIndex + 1,
                          currentMolAndEmbeddings, smc, thisGroup, thisAndor,
                          current_group_probability, regioIsomer)) {
      return 0;
    }
    return 1;
  }

  for (std::vector<IWString>::const_iterator smiIter = smilesList.begin();
       smiIter != smilesList.end(); ++smiIter) {
    regioIsomer++;

    IWString smiles;
    IWString id;
    (*smiIter).split(smiles, ' ', id);

    if (verbose > 1) {
      cerr << "smiles: " << smiles << "    id: " << id << '\n';
    }

    Molecule m;

    if (!m.build_from_smiles(smiles)) {
      cerr << "Set_of_Reactions::prepOneSideChain:invalid sidechain smiles " << smiles
           << "\n";
      return 0;
    }
    m.set_name(id);

    // we need a  reaction to do the search

    IWReaction rxn;
    set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(
        0);  // do not show the message
    if (!rxn.do_read(rxn_string, smc)) {
      cerr << "Set_of_Reactions::prepOneSideChain:cannot read '" << rxn_string << "'\n";
      return 0;
    }
    // IAW: note this is a memory leak...
    rxn.set_user_specified_void_ptr(
        new AndORProb(thisGroup, thisAndor, current_group_probability));

    // see how many mappings for this reaction  there are

    Substructure_Results sresults;
    if (verbose > 1) {
      cerr << "smc.process_hit_number(): " << smc.process_hit_number() << '\n';
    }
    const int nhits = rxn.sidechain(sideChainIndex)->substructure_search(m, sresults);

    if (verbose > 1) {
      cerr << "nhits: " << nhits << '\n';
    }

    if (0 == nhits) {
      cerr << "no hits to query 0 " << rxn.comment() << ' ' << m.smiles() << '\n';
      if (smc.ignore_not_reacting()) {
        return 1;
      }
      return 0;
    }

    int thisHit = smc.process_hit_number();

    if (nhits > 1 &&
        !(smc.make_new_reagent_for_each_hit() ||
          thisHit >= 0)) {  // multiple hits, but not doing all hits and do= not specified
      cerr << " smc.suppress_if_more_than_this_many_substructure_search_hits()= "
           << smc.suppress_if_more_than_this_many_substructure_search_hits() << '\n';
      if (smc.suppress_if_more_than_this_many_substructure_search_hits() < nhits) {
        cerr << "multiple hits to query 0 " << ' ' << rxn.comment() << ' ' << m.smiles()
             << ".  Ignored." << '\n';
        return 1;
      }
      cerr << "Error: multiple hits to query 0 " << ' ' << rxn.comment() << ' '
           << m.smiles() << '\n';
      return 0;
    }

    if (nhits == 1) {
      thisHit = 0;
    }

    if (thisHit >= nhits) {
      cerr << "Error: hit requested (" << thisHit << ") is out of range.  #hits=" << nhits
           << " for " << rxn.comment() << ' ' << m.smiles() << '\n';
      return 0;
    }
    for (int i = 0; i < nhits; ++i) {
      if (thisHit >= 0 && i != thisHit) {
        continue;
      }

      Molecule_and_Embedding q;
      q.add_molecule(&m);
      q.set_name(m.name());
      q.set_embedding(*sresults.embedding(i));
      currentMolAndEmbeddings[sideChainIndex] = &q;

      // process the next side chain recursively

      if (!prepOneSideChain(rxn_string, smilesLists, sideChainIndex + 1,
                            currentMolAndEmbeddings, smc, thisGroup, thisAndor,
                            current_group_probability, regioIsomer)) {
        return 0;
      }
    }
  }

  return 1;
}

/*
  Buffer will look like
  file.rxn|directives|sidechain id|sidechain id
*/

int
Set_of_Reactions::_read_reaction(const IWString& buffer, Sidechain_Match_Conditions& smc,
                                 int& group, const int andor,
                                 const double current_group_probability,
                                 int& regioIsomer) {
  static const char sep = '|';

  if (!buffer.contains(sep)) {
    if (verbose > 1) {
      cerr << "Found reaction without sidechain information: " << buffer << '\n';
    }

    IWReaction* r = new IWReaction();
    set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(
        0);  // do not show the message

    if (!r->do_read(buffer, smc)) {
      cerr << "Set_of_Reactions::_read_reaction:cannot read '" << buffer << "'\n";
      delete r;
      return 0;
    }

    if (NO_GROUPING != andor) {
      r->set_user_specified_void_ptr(new AndORProb(group, andor));
    } else {
      r->set_user_specified_void_ptr(new AndORProb(++group, NO_GROUPING));
    }

    _reaction.add(r);

    return 1;
  }

  int i = 0;
  IWString rxn_string;

  buffer.nextword_single_delimiter(rxn_string, i, sep);

  // Bug #19912: quick fix to trim trailing spaces
  while (rxn_string.back() == ' ') {
    rxn_string.chop();
  }
  // Bug #19912: end

  if (verbose > 1) {
    cerr << "Reaction: " << rxn_string << '\n';
  }

  _maybe_prepend_directory_name(rxn_string);

  IWString directives;

  buffer.nextword_single_delimiter(directives, i, sep);

  std::vector<std::vector<IWString>>
      smilesLists;  // one vector of smiles for each sidechain

  // get a reaction for checking things

  IWReaction rxn;
  set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(
      0);  // do not show the message
  if (!rxn.do_read(rxn_string, smc)) {
    cerr << "Set_of_Reactions::prepOneSideChain:cannot read '" << rxn_string << "'\n";
    return 0;
  }
  rxn.set_user_specified_void_ptr(new AndORProb(group, andor, current_group_probability));

  for (int ndx = 0;; ++ndx) {
    const_IWSubstring token;
    if (!buffer.nextword_single_delimiter(token, i, sep)) {
      break;
    }

    // make sure the number of sidechains in the reactin is at least as big as the ndx
    // index we are processing now

    if (rxn.number_sidechains() <= ndx) {
      cerr << "Set_of_Reactions::_read_reaction:too many |'s (tokens) in " << buffer
           << ".  Rxn only has " << rxn.number_sidechains() << "sidechains\n";
      return 0;
    }

    smilesLists.push_back(std::vector<IWString>());
    std::vector<IWString>& thisSmilesList = smilesLists.back();

    token.strip_leading_blanks();

    Sidechain_Reaction_Site* s = rxn.sidechain(ndx);

    // see if the reaction needs smiles from the config file.  If not, the token should be
    // blank

    if (s->single_reagent_only()) {
      if (token.length() > 0) {
        if (verbose > 1) {
          cerr << "Set_of_Reactions::_read_reaction:in '" << buffer
               << "', the smiles specifiecation for sidechain " << ndx + 1
               << " should be blank - the rxn does not require a regaent\n";
        }
        return 0;
      }
      continue;  // go the next sidechain
    }

    if (0 == token.length()) {
      cerr << "Error:Set_of_Reactions::_read_reaction:warning: zero length sidechain "
              "specification '"
           << buffer << "'\n";
      return 0;
    }

    if (token.starts_with("file=")) {
      if (verbose > 1) {
        cerr << "found a file specification for the sidechain structures:" << token
             << '\n';
      }

      IWString thisFile = token;
      thisFile.remove_leading_chars(5);
      _maybe_prepend_directory_name(thisFile);

      iwstring_data_source input(thisFile.null_terminated_chars());

      if (!input.good()) {
        cerr << "Set_of_Reactions::_read_reaction:cannot open '" << thisFile << "'\n";
        return 0;
      }

      const_IWSubstring thisSmilesAndName;
      while (input.next_record(thisSmilesAndName)) {
        if (verbose > 1) {
          cerr << "Read smiles: " << thisSmilesAndName << '\n';
        }
        thisSmilesList.push_back(thisSmilesAndName);
      }

      input.do_close();
    } else {
      if (verbose > 1) {
        cerr << "found an explicit smiles definition: " << token << '\n';
      }

      thisSmilesList.push_back(token);
    }
  }

  // collect all the smiles to process - it might be only one from a direct smiles spec,
  // or several from a file
  //  In any case each smiles might have multiple mappings
  //  we split each mapping or each smiles into separate reactions
  //
  //  Yes, if there are 300 molecules in the input file, we make 300 (or more) copies of
  //  the reaction

  // now get ready for the andor grouping.  If there is already one active, we use it.
  // If it is an AND grouping, each smiles and each mapping of the smiles becomes one more
  // reaction in the AND group, which means each one will generate one output record. if
  // the current grouping is an OR grouping, each mapping of each smiles is part of that
  // OR, and if any one actually matches, the rest of the OR will not be processed
  //
  // If there is no current grouping, we make An AND grouping.  This soes not hurt even if
  // there is only one mapping of one smiles because an AND with only one thing in it is
  // the same as no mapping group.

  int thisAndor;
  int thisGroup;

  if (NO_GROUPING != andor) {
    thisGroup = group;
    thisAndor = andor;
  } else {
    thisGroup = ++group;
    ++group;  // after this file of reactants, the group must start anew
    thisAndor = AND_GROUPING;
  }

  if (verbose > 1) {
    cerr << "thisGroup: " << thisGroup << "  thisAndor: " << thisAndor << '\n';
  }

  // now process each smiles for each sidechain - recursively

  std::vector<Molecule_and_Embedding*> currentMolAndEmbeddings(smilesLists.size());
  if (!prepOneSideChain(rxn_string, smilesLists,
                        0 /* first sidechain - rest are done recursively*/
                        ,
                        currentMolAndEmbeddings, smc, thisGroup, thisAndor,
                        current_group_probability, regioIsomer)) {
    return 0;
  }

  return 1;
}

int
Set_of_Reactions::_read_start_molecule(const IWString& s) {
  char sep = ' ';
  if (s.contains(',')) {
    sep = ',';
  }

  // cerr << "Set_of_Reactions::_read_start_molecule:start " << s << '\n';

  int i = 0;
  const_IWSubstring token;

  if (!s.nextword(token, i, sep)) {
    return 0;
  }

  if (!_start_molecule.build_from_smiles(token)) {
    cerr << "Invalid smiles " << token << '\n';
    return 0;
  }

  if (!s.nextword(token, i, sep)) {
    return 1;
  }

  _start_molecule.set_name(token);

  return 1;
}

int
Set_of_Reactions::_build(const resizable_array_p<IWString>& zdata,
                         Sidechain_Match_Conditions& smc) {
  if (_number_reactions > 0) {
    _free_reaction_array();
  }

  const int n = zdata.number_elements();

  int andor = NO_GROUPING;

  int group = 0;
  int regioIsomer = 0;

  double current_group_probability = 0.0;

  for (int i = 0; i < n; ++i) {
    IWString& s = *zdata[i];  // non const is OK

    if (s.starts_with("|")) {
      andor = 0;
    } else if (s.starts_with("OR:")) {
      andor = OR_GROUPING;
      current_group_probability = 0.0;
      group++;
    } else if (s.starts_with("AND")) {
      andor = AND_GROUPING;
      current_group_probability = 0.0;
      group++;
    } else if (s.starts_with("PROBABILITY:")) {
      s.remove_leading_chars(12);
      s.strip_leading_blanks();
      if (!s.numeric_value(current_group_probability) ||
          current_group_probability <= 0.0 || current_group_probability > 1.0) {
        cerr << "Set_of_Reactions::_build:invalid probability '" << s << "'\n";
        return 0;
      }
    } else if (s.starts_with("K:")) {
      s.remove_leading_chars(2);
      IWReaction* r = new IWReaction();
      // do not show the message
      set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(0);

      if (!r->construct_from_smirks(s)) {
        cerr << "Invalid smirks '" << s << "'\n";
        delete r;
        return 0;
      }

      if (NO_GROUPING != andor) {
        r->set_user_specified_void_ptr(
            new AndORProb(group, andor, current_group_probability));  // note never freed
      } else {
        r->set_user_specified_void_ptr(new AndORProb(
            ++group, NO_GROUPING, current_group_probability));  // group number is unique
      }

      _reaction.add(r);
    } else if (s.starts_with("START:")) {
      s.remove_leading_chars(6);
      if (!_read_start_molecule(s)) {
        cerr << "Set_of_Reactions::_build:cannot read starting molecule '" << s << "'\n";
        return 0;
      }
    } else if (NO_OP_RXN == s) {
      IWReaction* r = new IWReaction();
      // do not show the message
      set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(0);

      r->set_comment(NO_OP_RXN);

      if (NO_GROUPING != andor) {
        r->set_user_specified_void_ptr(
            new AndORProb(group, andor, current_group_probability));
      } else {
        r->set_user_specified_void_ptr(
            new AndORProb(++group, NO_GROUPING, current_group_probability));
      }
      _reaction.add(r);
    } else if (!_read_reaction(s, smc, group, andor, current_group_probability,
                               regioIsomer)) {
      cerr << "Set_of_Reactions::_build:cannot process control file record '" << s
           << "'\n";
      return 0;
    }
  }

  _number_reactions = _reaction.number_elements();

  if (verbose > 1) {
    cerr << "Read " << _number_reactions << " reactions\n";
  }

  for (int i = 0; i < _number_reactions; ++i) {
    if (nullptr != _reaction[i]->user_specified_void_ptr()) {
      continue;
    }

    cerr << "Reaction " << i << " no andor\n";
  }

  if (one_embedding_per_start_atom) {
    set_find_one_embedding_per_start_atom(1);
  }

  _hits_to_query = new_int(_number_reactions);

  return _number_reactions;
}

int
Set_of_Reactions::build(Sidechain_Match_Conditions& smc, int& fatal) {
  resizable_array_p<IWString> zdata;
  const_IWSubstring buffer;

  fatal = 1;

  while (_input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if ("$$$$" == buffer) {
      break;
    }

    if (0 == buffer.length()) {
      continue;
    }

    zdata.add(new IWString(buffer));
  }

  if (verbose) {
    cerr << "Read " << zdata.number_elements() << " records from control file\n";
  }

  if (zdata.empty()) {
    cerr << "Set_of_Reactions::build:no data\n";
    return 0;
  }

  return _build(zdata, smc);
}

int
Set_of_Reactions::next(Sidechain_Match_Conditions& smc, int& fatal) {
  fatal = 0;

  if (!_input.is_open()) {
    return 0;
  }

  if (_input.eof()) {
    return 0;
  }

  return build(smc, fatal);
}

int
Set_of_Reactions::open(const char* fname) {
  if (!_input.open(fname)) {
    cerr << "Set_of_Reactions::open:cannot open control file '" << fname << "'\n";
    return 0;
  }

  _input.set_skip_blank_lines(1);
  _input.set_strip_trailing_blanks(1);

  return 1;
}

int
Set_of_Reactions::_add_reaction_from_proto_file(const_IWSubstring& fname,
                                                Sidechain_Match_Conditions& smc) {
  IWString tmp(fname);  // Call below needs an IWString.
  std::optional<ReactionProto::Reaction> maybe_rxn =
      iwmisc::ReadTextProto<ReactionProto::Reaction>(tmp);
  if (!maybe_rxn) {
    cerr << "Set_of_Reactions::_add_reaction_from_proto_file:cannot read '" << fname
         << "'\n";
    return 0;
  }
  std::unique_ptr<IWReaction> rxn = std::make_unique<IWReaction>();
  if (!rxn->ConstructFromProto(*maybe_rxn, fname)) {
    cerr << "Set_of_Reactions::_add_reaction_from_proto_file:cannot parse proto\n";
    return 0;
  }
  _reaction.add(rxn.release());

  return 1;
}

int
Set_of_Reactions::_add_reaction_from_file_of_proto_reactions(
    const_IWSubstring& fname, Sidechain_Match_Conditions& smc) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Set_of_Reactions::_add_reaction_from_file_of_proto_reactions:cannot open '"
         << fname << "'\n";
    return 0;
  }
  const IWString dirname = GetDirName(fname);
  return _add_reaction_from_file_of_proto_reactions(input, dirname, smc);
}

int
Set_of_Reactions::_add_reaction_from_file_of_proto_reactions(
    iwstring_data_source& input, const IWString& dirname,
    Sidechain_Match_Conditions& smc) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with("#")) {
      continue;
    }
    IWString fname;
    fname << dirname << std::filesystem::path::preferred_separator << buffer;
    const_IWSubstring tmp(fname);
    if (!_add_reaction_from_proto_file(tmp, smc)) {
      cerr << "Set_of_Reactions::_add_reaction_from_file_of_proto_reactions:cannot "
              "process '"
           << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Set_of_Reactions::_add_reaction_from_smirks(const_IWSubstring& smirks,
                                            Sidechain_Match_Conditions& smc) {
  IWReaction* r = new IWReaction();

  if (!r->construct_from_smirks(smirks)) {
    cerr << "Set_of_Reactions::_add_reaction_from_smirks;invalid smirks " << smirks
         << '\n';
    delete r;
    return 0;
  }

  _reaction.add(r);

  return 1;
}

int
Set_of_Reactions::_add_reaction(const_IWSubstring& fname,
                                Sidechain_Match_Conditions& smc) {
  IWReaction* r = new IWReaction();

  if (!r->do_read(fname, smc)) {
    cerr << "Set_of_Reactions::_add_reaction:cannot read reaction " << fname << '\n';
    delete r;
    return 0;
  }

  _reaction.add(r);

  return 1;
}

int
Set_of_Reactions::build(Command_Line& cl, const char opt, Sidechain_Match_Conditions& smc,
                        const int verbose) {
  if (_number_reactions > 0) {
    _free_reaction_array();
  }

  _number_reactions = count_number_of_reactions(cl, opt);

  if (0 == _number_reactions) {
    cerr << "Set_of_Reactions::build:no reactions on command line -" << opt << '\n';
    return 0;
  }

  int i = 0;
  const_IWSubstring r;

  while (cl.value(opt, r, i++)) {
    if (r.starts_with("F:")) {
      r.remove_leading_chars(2);

      if (!_read_file_of_reactions(r, smc)) {
        cerr << "Invalid reaction specification in '" << r << "'\n";
        return 0;
      }
    } else if (r.starts_with("SMIRKS:")) {
      r.remove_leading_chars(7);
      if (!_add_reaction_from_smirks(r, smc)) {
        cerr << "Invalid smirks '" << r << "'\n";
        return 0;
      }
    } else if (r.starts_with("K:")) {
      r.remove_leading_chars(2);
      if (!_read_file_of_smirks(r)) {
        cerr << "Cannot read file of smirks '" << r << "'\n";
        return 0;
      }
    } else if (r.starts_with("proto:")) {
      r.remove_leading_chars(6);
      if (!_add_reaction_from_proto_file(r, smc)) {
        cerr << "Invalid proto file '" << r << "'\n";
        return 0;
      }
    } else if (r.starts_with("PROTO:")) {
      r.remove_leading_chars(6);
      if (!_add_reaction_from_file_of_proto_reactions(r, smc)) {
        cerr << "Invalid proto file '" << r << "'\n";
        return 0;
      }
    } else if (!_add_reaction(r, smc)) {
      cerr << "Cannot create reaction from '" << r << "'\n";
      return 0;
    }
  }

  _number_reactions = _reaction.number_elements();

  // cerr << "LINE " << __LINE__ << " j " << j << " _number_reactions " <<
  // _number_reactions << '\n';

  // Set a default AndORProb if needed.

  const int grp = 1;
  for (int i = 0; i < _number_reactions; ++i) {
    if (nullptr != _reaction[i]->user_specified_void_ptr()) {  //
      continue;
    }

    AndORProb* a = new AndORProb(grp, 0);
    _reaction[i]->set_user_specified_void_ptr(a);
  }

  _hits_to_query = new_int(_number_reactions);

  return _number_reactions;
}

static int
display_miscellaneous_options(std::ostream& output, const char flag) {
  output << "  -" << flag << " append      append sidechain names if present\n";
  output << "  -" << flag << " apprxn      append reaction names if present\n";
  output << "  -" << flag << " apprgnt     append reagent  names if present\n";
  output << "  -" << flag
         << " aeuao       Atom Environments match Unmatched Atoms Only\n";
  output << "  -" << flag << " orspm       one reaction scheme per input molecule\n";

  exit(1);
}

static void
display_dash_W_options(std::ostream& output) {
  output << " -W rgsep=<...>     separator between reaction name and reagent \n";
  output << " -W rxsep=<...>     separator between reaction\n";

  exit(0);
}

static int
molecular_transformations(int argc, char** argv) {
  Command_Line cl(argc, argv,
                  "vA:E:R:t:i:o:S:m:z:M:ZX:C:jkF:g:bhp:J:Y:uO:lLdU:a:HT:W:P:N:es:r");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot parse element specifications (-E)\n";
      return 5;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialis chemical standardisation\n";
      return 4;
    }
  }

  if (cl.option_present('t')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 't')) {
      usage(8);
    }
  }

  if (cl.option_present('I')) {
    convert_isotopes = 1;
    if (verbose) {
      cerr << "Isotopes in product molecules will be converted\n";
    }
  }

  if (cl.option_present('C')) {
    int i = 0;
    const_IWSubstring c;

    while (cl.value('C', c, i++)) {
      if ("RXN" == c) {
        append_reaction_names = 1;

        if (verbose) {
          cerr << "Will append the name of the reactions changing the molecule\n";
        }
      } else {
        append_to_changed_molecules = cl.string_value('C');

        if (verbose) {
          cerr << "Will append '" << append_to_changed_molecules
               << "' and reaction number to changed molecules\n";
        }
      }
    }
  }

  if (cl.option_present('j')) {
    write_original_molecule_to_output_stream = 1;

    if (verbose) {
      cerr << "Will write initial molecules to output stream as well\n";
    }
  }

  if (cl.option_present('k')) {
    write_parent_of_changed_molecules = 1;

    if (verbose) {
      cerr << "Will write initial molecules to output stream if it is changed\n";
    }
  }

  if (cl.option_present('l')) {
    strip_input_molecules_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Input molecules stripped to largest fragment\n";
    }
  }

  if (cl.option_present('L')) {
    strip_products_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Product molecules stripped to largest fragment\n";
    }
  }

  if (cl.option_present('W')) {
    IWString w;
    for (int i = 0; cl.value('W', w, i); ++i) {
      if (w.starts_with("rgsep=")) {
        reaction_reagent_separator = w;
        reaction_reagent_separator.remove_leading_chars(6);
        append_reagent_names = 1;
      } else if (w.starts_with("rxsep=")) {
        name_separator = w;
        name_separator.remove_leading_chars(6);
        name_separator.unhtml();
        append_reaction_names = 1;
      } else if ("help" == w) {
        display_dash_W_options(cerr);
      } else {
        name_separator = w;
        name_separator.unhtml();

        if (verbose) {
          cerr << "Name separator set to '" << name_separator << "'\n";
        }
        set_component_separator(name_separator);
      }
    }
  }

  if (cl.option_present('O')) {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('O', m, i++)) {
      if ("append" == m) {
        append_reagent_names = 1;

        if (verbose) {
          cerr << "Will concatenate reagent names\n";
        }

      } else if ("apprxn" == m) {
        append_reaction_names = 1;
        if (verbose) {
          cerr << "Will append reaction names to changed molecules\n";
        }
      } else if ("apprgnt" == m) {
        append_reagent_names = 1;
        if (verbose) {
          cerr << "Will append reagent names to changed molecules\n";
        }
      } else if ("aeuao" == m) {
        set_atom_environment_only_matches_unmatched_atoms(1);
      } else if ("orspm" == m) {
        one_reaction_scheme_per_input_molecule = 1;
        if (verbose) {
          cerr << "One reaction scheme per input molecule\n";
        }
      } else if ("help" == m) {
        display_miscellaneous_options(cerr, 'O');
      } else {
        cerr << "Unrecognised -O qualifier '" << m << "'\n";
        display_miscellaneous_options(cerr, 'O');
      }
    }
  }

  if (cl.option_present('a')) {
    const_IWSubstring a = cl.string_value('a');
    {
      if ("rej" == a) {
        preserve_initial_aromaticity = 2;
        if (verbose) {
          cerr << "Will reject any molecules losing aromaticity\n";
        }
      } else if ("write" == a) {
        preserve_initial_aromaticity = 1;
        if (verbose) {
          cerr << "Molecules losing aromaticity will be written\n";
        }
      } else {
        cerr << "Unrecognised -a qualifier '" << a << "'\n";
        usage(1);
      }
    }
  }

  if (cl.option_present('F') && cl.option_present('i')) {
    cerr << "Pipe functionality specified (-F), input option (-i) ignored\n";
  }

  if (cl.option_present('F')) {
    int i = 0;
    const_IWSubstring f;
    while (cl.value('F', f, i++)) {
      if (f.starts_with("IN=")) {
        f.remove_leading_chars(3);
        smiles_in = f;
        if (!smiles_in.ends_with('<')) {
          smiles_in += '<';
        }

        if (verbose) {
          cerr << "Will work as a pipe, smiles input as '" << smiles_in << "'\n";
        }
      } else if (f.starts_with("OUT=")) {
        f.remove_leading_chars(4);
        smiles_out = f;
        if (!smiles_out.ends_with('<')) {
          smiles_out += '<';
        }

        if (verbose) {
          cerr << "Will work as a pipe, smiles output as '" << smiles_out << "'\n";
        }
      } else {
        cerr << "Unrecognised -F qualifier '" << f << "'\n";
        usage(4);
      }

      if (smiles_in.length() && 0 == smiles_out.length()) {
        smiles_out = smiles_in;
      }
    }
  }

  if (cl.option_present('b')) {
    break_after_first_match = 1;

    if (verbose) {
      cerr << "Will stop processing after the first reaction matches\n";
    }
  }

  if (cl.option_present('h')) {
    make_implicit_hydrogens_explicit = 1;

    if (verbose) {
      cerr << "Implicit Hydrogens will be made explicit\n";
    }
  }

  if (cl.option_present('p')) {
    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++)) {
      if ("ovok" == p) {
        or_matches_can_overlap = 1;

        if (verbose) {
          cerr << "Parallel embeddings can overlap DANGEOUS if atoms removed\n";
        }
      } else if (!p.numeric_value(or_depth) || or_depth < 1) {
        cerr << "Invalid OR depth '" << p << "'\n";
        usage(8);
      } else if (verbose) {
        cerr << "Will perform all " << or_depth << " combinations of changes\n";
      }
    }
  }

  if (cl.option_present('X')) {
    if (!elements_to_remove.construct_from_command_line(cl, verbose, 'X')) {
      cerr << "Cannot discern elements to remove from -X switch\n";
      usage(18);
    }
  }

  if (cl.option_present('d')) {
    suppress_duplicate_molecules = 1;

    if (verbose) {
      cerr << "Will suppress duplicate molecules - within current generation\n";
    }
  }

  if (cl.option_present('s')) {
    unsigned int seed;
    if (!cl.value('s', seed)) {
      cerr << "Invaid seed (-s)\n";
      return 1;
    }
    rng.seed(seed);
  }

  if (cl.option_present('U')) {
    const_IWSubstring u = cl.string_value('U');

    stream_for_molecules_not_reacting.add_output_type(FILE_TYPE_SMI);

    if (stream_for_molecules_not_reacting.would_overwrite_input_files(cl, u)) {
      cerr << "Cannot overwrite input file(s) '" << u << "'\n";
      return 2;
    }

    if (!stream_for_molecules_not_reacting.new_stem(u)) {
      cerr << "Cannot open stream for non reacting molecules '" << u << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Non reacting molecules written to '" << u << "'\n";
    }
  }

  set_copy_name_in_molecule_copy_constructor(1);

  FileType input_type = FILE_TYPE_INVALID;

  if (smiles_in.length()) {
    ;
  } else if (!cl.option_present('i')) {
    if (all_files_recognised_by_suffix(cl)) {
      ;
    } else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
      input_type = FILE_TYPE_SMI;
    } else {
      cerr << "Cannot discern file types from names\n";
      return 4;
    }
  } else {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (cl.option_present('R') && cl.option_present('T')) {
    cerr << "Sorry, the -R and -T options cannot be used together, see Ian\n";
    return 1;
  }

  Set_of_Reactions sor;

  if (append_reaction_names) {
    sor.set_append_names(1);
  }

  Sidechain_Match_Conditions sidechain_match_conditions;
  sidechain_match_conditions.set_verbose(verbose);

  if (cl.option_present('M')) {
    const_IWSubstring m;

    cl.value('M', m);

    if ("all" == m || "each" == m) {
      sidechain_match_conditions.set_make_new_reagent_for_each_hit(1);
      if (verbose) {
        cerr << "All regio-isomers will be generated\n";
      }
    } else if (m.starts_with("do=")) {
      m.remove_leading_chars(3);

      int mm;

      if (!m.numeric_value(mm) || mm < 0) {
        cerr << "The '-M do=' option must be followed by a whole non-negative number\n";
        usage(99);
      }

      sidechain_match_conditions.set_process_hit_number(mm);
    }
    //    else if (m.starts_with("mskip="))
    //    {
    //      const_IWSubstring ms = m;
    //      ms.remove_leading_chars(6);
    //      sidechain_match_conditions.set_multiple_match_string(ms);
    //      if (verbose)
    //        cerr << "Where just one of several possible sidechain matches used, will
    //        append '" << ms << "'\n";
    //    }
    else if (m.starts_with("write=")) {
      m.remove_leading_chars(6);
      if (!set_stream_for_sidechains_not_matching_query(m)) {
        cerr << "Cannot initialise stream for non reacting sidechains '" << m << "'\n";
        return 0;
      }

      if (verbose) {
        cerr << "Sidechains not reacting written to '" << m << ".smi\n";
      }

      sidechain_match_conditions.set_ignore_not_reacting(1);
    } else if ("RMX" == m) {
      sidechain_match_conditions.set_ignore_multiple_substucture_matches(1);
    } else {
      cerr << "Unrecognised -M qualifier '" << m << "'\n";
      usage(23);
    }
  }

  if (cl.option_present('Z')) {
    sidechain_match_conditions.set_ignore_not_reacting(1);

    if (verbose) {
      cerr << "Reagents not matching the query will be ignored\n";
    }
  }

  if (cl.option_present('R')) {
    if (!sor.build(cl, 'R', sidechain_match_conditions, verbose)) {
      cerr << "Cannot read reactions (-R)\n";
      return 1;
    }

    if (verbose) {
      cerr << "Defined " << sor.number_reactions() << " reactions\n";
    }
  } else if (cl.option_present('T')) {
    const char* fname = cl.option_value('T');

    if (!sor.open(fname)) {
      cerr << "Cannot open set of reactions specification '" << fname << "'\n";
      return 1;
    }

    multiple_reaction_schemes = sor.multiple_reaction_schemes_in_input();
    if (verbose && multiple_reaction_schemes) {
      cerr << "Found " << multiple_reaction_schemes
           << " reaction schemes in control file\n";
    }
  } else {
    cerr << "Must specify reactions via the -R option, or a control file via the -T "
            "option\n";
    usage(1);
  }

  if (cl.option_present('H')) {
    ok_atom_removal_reaction = 1;

    if (verbose) {
      cerr << "With OR reactions, it will be OK to have atom removals. Beware!\n";
    }
  }

  if (or_depth) {  // we cannot have atom removals.
    if (!sor.check_for_atom_removals(ok_atom_removal_reaction)) {
      return 1;
    }
  }

  sor.set_find_unique_embeddings_only(1);  // not sure if this is the right idea or not

  // cerr << cl.option_present('u') << " setting " << number_reactions << " reactions\n";

  if (cl.option_present('u')) {
    one_embedding_per_start_atom = 1;
    sor.set_find_one_embedding_per_start_atom(1);
    sidechain_match_conditions.set_one_embedding_per_start_atom(1);
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    int i = 0;
    while (cl.value('z', z, i++)) {
      if ('w' == z) {
        write_molecules_not_reacting = 1;
        ignore_molecules_not_reacting = 1;
        if (verbose) {
          cerr << "Molecules not reacting will be written to the output stream\n";
        }
      } else if ('d' == z) {
        discard_molecules_not_reacting = 1;
        if (verbose) {
          cerr << "Molecules not reacting will be discarded\n";
        }
      } else if ('i' == z) {
        ignore_molecules_not_reacting = 1;
        if (verbose) {
          cerr << "Molecules not reacting will be ignored\n";
        }
      } else {
        cerr << "Unrecognised -z directive '" << z << "'\n";
        usage(19);
      }
    }
  }

  if (cl.option_present('r')) {
    ignore_sidechains_causing_increased_core_hit_count = 1;

    if (verbose) {
      cerr << "Reagents not matching the query will be ignored\n";
    }
  }

  // Apply -m after reading from file so we override any value in the file
  // Note that the variable mdo is global to this file

  if (cl.option_present('m')) {
    int i = 0;
    const_IWSubstring m;
    // int mmx = -1;
    mdo = -1;  // global variable
    while (cl.value('m', m, i++)) {
      if (m.starts_with("do=")) {
        m += 3;

        if (!m.numeric_value(mdo) || mdo < 0) {
          cerr << "The '-m do=' option must be followed by a whole non-negative number\n";
          usage(29);
        }
      } else if ("each" == m || "all" == m) {
        enumerate_scaffold_hits_individually = 1;
        // enumerate_scaffold_hits_individually = 2;
        if (verbose) {
          cerr << "Will enumerate multiple scaffold hits individually and collectively\n";
        }
      } else if ("RMX" == m) {
        ignore_scaffold_multiple_substucture_matches = 1;
        if (verbose) {
          cerr << "Will discard scaffolds that show multiple substructure search hits\n";
        }
      }
      //      else if ("each1" == m)
      //      {
      //        enumerate_scaffold_hits_individually = 1;
      //        if (verbose)
      //          cerr << "Will enumerate multiple scaffold hits individually\n";
      //      }
      //      else if (m.numeric_value(mmx))
      //      {
      //        if (! m.numeric_value(mmx) || mmx < 1)
      //        {
      //          cerr << "The '-m <number>' option must be followed by a whole positive
      //          number\n"; usage(91);
      //        }
      //        enumerate_scaffold_hits_individually = 2;
      //      }
      else {
        cerr << "Unrecognised -m qualifier '" << m << "'\n";
        usage(37);
      }
    }

    //  If multiple forms are specified, they must be compatible

    if (enumerate_scaffold_hits_individually && mdo >= 0) {
      cerr << "Enumerating multiple scaffold hits (-m each) is incompatible with '-m "
              "do=nn'\n";
      usage(51);
    }

    //    if (mmx > 0 && mdo >= 0)
    //    {
    //      if (mdo >= mmx)
    //      {
    //        cerr << "You asked to process scaffold hit number " << mdo << '\n';
    //        cerr << "You asked to process a maximum of " << mmx << " scaffold hits\n";
    //        cerr << "Impossible\n";
    //        usage(28);
    //      }
    //    }
    //    if (mmx > 0)
    //    {
    //      //sor.set_max_matches_to_find(mmx);   //sor is not yet loaded
    //      if (verbose)
    //        cerr << "Max matches to find is " << m << '\n';
    //    }

    if (mdo >= 0) {
      if (verbose) {
        cerr << "Will process hit number " << mdo << " of scaffold hits\n";
      }
    }
  }

  if (cl.option_present('e')) {
    unique_products_only = 1;

    if (verbose) {
      cerr << "Will suppress duplicate products\n";
    }
  }

  if (cl.option_present('P')) {
    double p;
    if (!cl.value('P', p) || p <= 0.0 || p > 1.0) {
      cerr << "Invalid reaction probability (-P)\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will sample reactions with probability " << p << '\n';
    }

    probability_sample_reactions = p;
  }

  if (cl.option_present('N')) {
    if (!cl.value('N', molecules_to_produce) || molecules_to_produce < 1) {
      cerr << "The number of molecules to produce (-N) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will produce " << molecules_to_produce << " random selections\n";
    }
  }

  if (smiles_in.length()) {  // special processing for TDT type input
    int fatal = 0;
    if (!sor.next(sidechain_match_conditions, fatal)) {
      cerr << "Cannot read reactions\n";
      return 1;
    }

    iwstring_data_source input("-");
    if (!input.good()) {
      cerr << "Yipes, cannot open stdin\n";
      return 4;
    }

    if (!molecular_transformations_tdt(input, sor)) {
      return 4;
    }

    if (verbose) {
      report_results(sor, cerr);
    }

    return 0;
  }

  if (smiles_out.length()) {  // output written in TDT form
    ;
  } else if (!cl.option_present('o')) {
    output_stream.add_output_type(FILE_TYPE_SMI);
  } else if (!output_stream.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    return 9;
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    // if the string contains an extention, remove it - it should be  a stem

    for (int i = s.length() - 1; i >= 0; --i) {
      const_IWSubstring thisChar = s.substr(i, 1);
      if (thisChar == ".") {
        s = s.substr(0, i);
        break;
      }

      if (thisChar == "/") {
        break;
      }
    }

    if (output_stream.would_overwrite_input_files(cl, s)) {
      cerr << "Cannot overwrite input file(s)\n";
      return 4;
    }

    if (!output_stream.new_stem(s, 1)) {
      cerr << "Cannot initialise output stem '" << s << "'\n";
      return 8;
    }
  } else if (smiles_out.length()) {
    ;
  } else {
    if (!output_stream.new_stem("-")) {
      cerr << "Cannot initialise stdout\n";
      return 3;
    }
  }

  int rc = 0;
  if (cl.empty()) {
    if (!sor.reactions_have_starting_molecules()) {
      cerr << "Insufficient arguments\n";

      usage(1);
    }

    if (!molecular_transformations(sor, sidechain_match_conditions)) {
      rc = 1;
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!molecular_transformations(cl[i], input_type, sor,
                                     sidechain_match_conditions)) {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose) {
    report_results(sor, cerr);
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = molecular_transformations(argc, argv);

  return rc;
}
