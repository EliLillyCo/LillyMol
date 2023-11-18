/*
  Performs retrosynthesis based on templates
*/

#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <atomic>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/rxn_file.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/retrosynthesis.pb.h"

typedef unsigned int atype_t;

using std::cerr;
using std::endl;

using iw_tf_data_record::TFDataWriter;

struct Error : std::exception {
  char text[1000];

  Error(char const* thisMessage)
  {
    sprintf(text, "%s", thisMessage);
  }

  char const*
  what() const throw()
  {
    return text;
  }
};

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_deconstructed = 0;

static int molecules_produced = 0;

static Chemical_Standardisation input_molecule_chemical_standardisation;
static Chemical_Standardisation product_molecule_chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int strip_product_to_largest_fragment = 0;

static int lower_atom_count_products = 0;

static int add_output_record_with_small_fragments_removed = 0;  // Jibo

static Atom_Typing_Specification atom_typing_specification;

static int* molecules_deconstructed_at_radius = nullptr;

static int* deconstructions_at_radius = nullptr;

// static resizable_array_p<IWReaction> * rxn;

// static int * reactions_performed

static int take_first_of_multiple_reagents = 1;

static int skip_reactions_with_multiple_ragents = 1;

static int echo_initial_molecule_before_each_product = 1;

static int discard_existing_atom_map = 0;

static resizable_array<int> changing_atoms_should_be;

static int ignore_molecules_not_reacting = 1;

static int molecules_not_reacting = 0;

static int remove_widows = 0;

static IWString_and_File_Descriptor stream_for_non_reacting_molecules;
static IWString_and_File_Descriptor reactantsStream;
static IWString_and_File_Descriptor reactant1Stream;
static IWString_and_File_Descriptor reactant2Stream;
static IWString_and_File_Descriptor reactantOnly1Stream;
static IWString outputFilename = "";
static IWString reactantFilenameBase = "";

// static int fragments_produced = 0;

// static int found_in_db = 0;

static int split_products_into_separate_molecules = 1;

static int break_after_first_deconstruction = 0;

static Report_Progress report_progress;

static int multi_threaded = 0;

static int reaction_name_is_file_name = 0;
static int reaction_name_is_path_name = 0;

static extending_resizable_array<int> acc_nhits;

static int regressionFlag = 0;

/*
  Various things for Molecule_to_Query_Specifications
*/

// static int fix_connectivity = 1;
static int aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic = 0;
static int atoms_conserve_ring_membership = 0;
static int preserve_saturation = 0;
static int preserve_ring_size = 0;

std::ofstream* queryEchoStreamPtr = nullptr;
// std::ofstream *msiEchoStreamPtr = nullptr;
static IWString rxn_name_header("");

static int skip_molecules_with_multiple_scaffold_embeddings = 0;
static int molecules_skipped_for_multiple_scaffold_embeddings = 0;
static IWString_and_File_Descriptor
    stream_for_molecules_with_multiple_scaffold_embeddings;

static int preserve_kekule_forms = 0;

static int input_is_reaction_smiles = 0;

static int max_matches_each_reaction = std::numeric_limits<int>::max();

static int molecules_failing = 0;
static int keep_going_after_test_failure = 0;
static int ignore_no_matching_reaction_test = 0;
static int molecules_with_no_corresponding_reaction = 0;

static int unfix_implicit_hydrogens_in_products = 0;
static int unfix_implicit_hydrogens_in_reagents = 0;

static IW_STL_Hash_Map<IWString, Molecule*> test_results;

static int remove_isotopes_from_incoming_molecules = 0;

static int ignore_reactions_with_no_changing_atoms = 0;

static int ignore_bad_reactions = 0;

static int write_results_as_proto = 0;

/*
  For each reaction, we need to keep track of how many reactions
  it is able to perform, as a function of radius
  Definitely no thread safety here
*/

class Reaction_With_Stats : public IWReaction
{
 private:
  int _searches_done;
  int _matches_found;

 public:
  Reaction_With_Stats();

  void
  another_search()
  {
    _searches_done++;
  }

  void
  another_match()
  {
    _matches_found++;
  }

  int
  searches_done() const
  {
    return _searches_done;
  }

  int
  matches_found() const
  {
    return _matches_found;
  }
};

Reaction_With_Stats::Reaction_With_Stats()
{
  _searches_done = 0;
  _matches_found = 0;

  return;
}

static void
preprocess(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (input_molecule_chemical_standardisation.active()) {
    input_molecule_chemical_standardisation.process(m);
  }

  if (remove_isotopes_from_incoming_molecules) {
    m.transform_to_non_isotopic_form();
    m.unset_unnecessary_implicit_hydrogens_known_values();
  }

  if (unfix_implicit_hydrogens_in_reagents) {
    m.unset_unnecessary_implicit_hydrogens_known_values();
  }

  return;
}

static int
handle_non_reacting_molecule(Molecule& m)
{
  molecules_not_reacting++;

  if (stream_for_non_reacting_molecules.is_open()) {
    stream_for_non_reacting_molecules << m.smiles() << ' ' << m.name() << '\n';
    stream_for_non_reacting_molecules.write_if_buffer_holds_more_than(4096);
  }

  return ignore_molecules_not_reacting;
}

static int
add_to_test_results(const const_IWSubstring& buffer,
                    IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  IWString smiles, id;

  if (!buffer.split(smiles, ' ', id)) {
    cerr << "add_to_test_results:must have smiles and ID '" << buffer << "' invalid\n";
    return 0;
  }

  std::unique_ptr<Molecule> m = std::make_unique<Molecule>();

  if (!m->build_from_smiles(smiles, 1)) {
    cerr << "add_to_test_results:invalid smiles " << buffer << endl;
    return 0;
  }

  m->set_name(id);

  m->reset_all_atom_map_numbers();
  m->unset_unnecessary_implicit_hydrogens_known_values();
  // cerr << "After unset " << m->smiles() << endl;

  if (test_results.contains(id)) {
    cerr << "add_to_test_results:duplicate results for '" << id << "', ignored\n";
    return 1;
  }

  test_results[id] = m.release();

  return test_results.size();
}

static int
read_test_results(iwstring_data_source& input,
                  IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!add_to_test_results(buffer, test_results)) {
      cerr << "Fatal erorr processing test results " << buffer << "\n";
      return 0;
    }
  }

  return test_results.size();
}

static int
read_test_results(const_IWSubstring& fname,
                  IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open test results '" << fname << "'\n";
    return 0;
  }

  return read_test_results(input, test_results);
}

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
  cerr << "  -R <rad>      maximum radius around changing atoms\n";
  cerr << "  -F <fname>    file containing single set of reaction files\n";
  cerr << "  -G <fname>    file containing names of sets of reaction files\n";
  cerr << "  -T <fname>    individual reactions\n";
  cerr << "  -I <fname>    reaction smiles input\n";
  cerr << "  -s            input reactions are in rxn smiles format\n";
  cerr << "  -o            output filename - if omitted, output is to stderr.\n";
  cerr << "                If specified, the file is only created if hits are found\n";
  cerr << "  -O            filename base for reactants - if specified, all reactants are written\n";
  cerr << "                to <base>_all.smi\n";
  cerr << "                If there are two or more reactants, the first two are\n";
  cerr << "                written to <base>_1.smi and <base.2.smi.\n";
  cerr << "                If there is only 1 reactant, it is writen to <base>_only1.smi\n";
  cerr << "  -P ...        atom typing - determine changing atoms (default=UST:AZUCORS)\n";
  cerr << "  -q p          use the reaction path name as the reaction name\n";
  cerr << "  -q f          use the reaction file name as the reaction name\n";
  cerr << "  -b            only allow each molecule to be deconstructed by one reaction\n";
  cerr << "  -m .          skip molecules with multiple scaffold embeddings\n";
  cerr << "  -m <fname>    write molecules with multiple scaffold embeddings to <fname>\n";
  cerr << "  -Q <fname>    file name into which the queries are saved (optional)\n";
  cerr << "  -n <header>   When queries are written out, this string is pre-pended to the name\n";
  cerr << "  -U <fname>    write non reacting molecules to <fname>\n";
  cerr << "  -a <natoms>   do NOT write products with  fewer than <natoms> atoms\n";
//cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -L            strip products to largest fragment\n";
  cerr << "  -r <number>   report progress every <number> molecules processed\n";
  cerr << "  -u prod       unfix redundant implicit hydrogens in product molecules\n";
  cerr << "  -u reac       remove redundant square brackets on input molecules (for testing)\n";
  cerr << "  -z            ignore reactions with no changing atoms\n";
  cerr << "  -Z            ignore reactions that cannot be constructed\n";
  cerr << "  -S            suppress files names in output (for regression tests)\n";
  cerr << "  -X ...        miscellaneous options, enter '-X help' for info\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        input molecule chemical standardisation options\n";
  cerr << "  -Y ...        output molecule chemical standardisation options\n";
//cerr << "  -E ...        standard element specifications\n";
//cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  The atoms that are really part of the reaction core will have 1 == changed,
  while those that are just Kekule additions will have 2 == changed.
  We want to expand to all atoms that are within RADIUS or the 1 atoms
*/

static int
identify_atoms_within_range(Molecule& m, const int radius, int* changed)
{
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (1 == changed[i]) {  // we only add things relative to the actual reaction core
      changed[i] = matoms;
    }
  }

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (matoms != changed[i]) {
      continue;
    }

    const int i_fragment_membership = m.fragment_membership(i);

    for (int j = 0; j < matoms; ++j) {
      if (changed[j]) {
        continue;
      }

      if (i_fragment_membership != m.fragment_membership(j)) {
        continue;
      }

      const int d = m.bonds_between(i, j);

      if (d <= radius) {
        changed[j] = 1;
        rc++;
      }
    }
  }

  for (int i = 0; i < matoms; ++i) {
    if (matoms == changed[i]) {
      changed[i] = 1;
      rc++;
    }
  }

  // the following code was taken from rxn_signature (which was taken
  // from rxn_fingerprint) to include any atoms that are double or
  // triple bonded to the core+radius atoms.  This was added here to get
  // the reaction seaching in retrosynthesis to agree with the the
  // production of the template clusters which is based on the
  // signatures Tad Hurst - 20180621

  // TODO: ianwatson. I think this could happen by looping over bonds
  // and looking for multiple bonds where both ends change. Much more efficient.

  m.compute_aromaticity_if_needed();

  Set_of_Atoms extras;

  for (int i = 0; i < matoms; ++i) {
    if (0 == changed[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    const int acon = a->ncon();

    if (acon == m.nbonds(i)) {  // fully saturated, no double bonds here
      continue;
    }

    for (int j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      if (b->is_single_bond() || b->is_aromatic()) {
        continue;
      }

      const int k = b->other(i);

      if (changed[k]) {
        continue;
      }

      extras.add(k);
      rc++;
    }
  }

  extras.set_vector(changed, 1);

  return rc;
}

static int
identify_changing_atoms_reagent(RXN_File& rxn, const int reagent, int* changed)
{
  Changing_Atom_Conditions cac;
  cac.set_is_changing_if_different_neighbours(1);
  cac.set_ignore_lost_atom_if_isolated(1);

  cac.set_include_changing_bonds_in_changing_atom_count(1);
  cac.set_consider_aromatic_bonds(0);

  return rxn.identify_atoms_changing_reagent(reagent, atom_typing_specification, changed,
                                             cac);
}

class Set_of_Reactions
{
 private:
  resizable_array_p<Reaction_With_Stats>* _rxn;
  int _max_radius;
  int _number_reactions;
  IWString _name;

  //  Some things to keep track of molecules we process. No thread safety here

  int* _molecules_deconstructed_at_radius;

  //   During testing it can be handy to know where a given reaction might
  //   be.  One per radius

  IW_STL_Hash_Map_int* _id_to_ndx;

  //  private functions

  int _determine_number_reactions();

  void _initialise_radius_dependent_arrays(const resizable_array<int>& radius);

  int _build(const const_IWSubstring& buffer, const resizable_array<int>& radius);

  int _process(Molecule& m, Molecule_to_Match& target, int* already_reacted,
               const int radius, IWString_and_File_Descriptor& output);
  int _process(Molecule& m, Molecule_to_Match& target, const int radius,
               Reaction_With_Stats* rxn, std::vector<IWString>& uniqResultSmi,
               IWString_and_File_Descriptor& output);
  int _process(Molecule& m, Molecule_to_Match& target, const int radius,
               Reaction_With_Stats* rxn, IW_STL_Hash_Set& usmi_here,
               IWString_and_File_Descriptor& output);

  void _reportResultDuplicates(const IWString& nameOfReaction,
                               const IWString& nameOfMolecule,
                               const IW_STL_Hash_Set& newResultSmi,
                               IW_STL_Hash_Set& all_smiles_formed);

 public:
  Set_of_Reactions();
  ~Set_of_Reactions();

  const IWString&
  name() const
  {
    return _name;
  }

  void
  set_name(const const_IWSubstring& s) {
    _name = s;
    cerr << "name now '" << _name << "'\n";
  }

  int
  number_reactions() const {
    return _number_reactions;
  }

  int report() const;
  int reaction_statistics() const;

  int build(const char* fname, const resizable_array<int>& radius);
  int build(iwstring_data_source& input, const resizable_array<int>& radius);
  int build_from_reaction_files(const const_IWSubstring&,
                                const resizable_array<int>& radius);

  int process(Molecule& m, Molecule_to_Match& target,
              IWString_and_File_Descriptor& output);
  int process(Molecule& m, Molecule_to_Match& target, const int radius,
              IWString_and_File_Descriptor& output);
  int process(Molecule& m, const int radius, IWString_and_File_Descriptor& output);

  Reaction_With_Stats* reaction_with_id(const IWString& id, const int radius) const;

  // Write reaction protos for each reaction in _rxn to
  // a TfDataRecord file `fname`.
  int WriteReactionProto(IWString& fname, int radius) const;
  int WriteReactionProto(iw_tf_data_record::TFDataWriter& writer, int radius) const;

  // Write the query implied by each LHS to a TfDataRecord
  // file `fname`.
  int WriteLhsQueryProto(IWString& name_stem, int radius) const;
  int WriteLhsQueryProto(iw_tf_data_record::TFDataWriter& writer, int radius) const;
};

Set_of_Reactions::Set_of_Reactions()
{
  _rxn = nullptr;
  _max_radius = -1;
  _number_reactions = 0;

  _molecules_deconstructed_at_radius = nullptr;

  _id_to_ndx = nullptr;

  return;
}

Set_of_Reactions::~Set_of_Reactions()
{
  if (nullptr != _rxn) {
    delete[] _rxn;
  }

  if (nullptr != _molecules_deconstructed_at_radius) {
    delete[] _molecules_deconstructed_at_radius;
  }

  if (nullptr != _id_to_ndx) {
    delete[] _id_to_ndx;
  }

  return;
}

int
Set_of_Reactions::build(const char* fname, const resizable_array<int>& radius)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Set_of_Reactions::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input, radius);
}

void
Set_of_Reactions::_initialise_radius_dependent_arrays(const resizable_array<int>& radius)
{
  _max_radius = radius.last_item();

  _rxn = new resizable_array_p<Reaction_With_Stats>[_max_radius + 1];

  _id_to_ndx = new IW_STL_Hash_Map_int[_max_radius + 1];

  _molecules_deconstructed_at_radius = new_int(_max_radius + 1);

  return;
}

int
Set_of_Reactions::_determine_number_reactions()
{
  _number_reactions = 0;

  for (int i = 0; i <= _max_radius; ++i) {
    if (_rxn[i].empty()) {
      continue;
    }

    _number_reactions = _rxn[i].number_elements();
    break;
  }

  return _number_reactions;
}

int
Set_of_Reactions::build(iwstring_data_source& input, const resizable_array<int>& radius)
{
  assert(radius.number_elements() > 0);
  assert(nullptr == _rxn);

  _initialise_radius_dependent_arrays(radius);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    // cerr << "Set_of_Reactions::build:processing '" << buffer << "',
    // input_is_reaction_smiles " << input_is_reaction_smiles << endl;
    if (!_build(buffer, radius)) {
      cerr << "Set_of_Reactions:build:cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return _determine_number_reactions();
}

int
Set_of_Reactions::build_from_reaction_files(const const_IWSubstring& buffer,
                                            const resizable_array<int>& radius)
{
  _initialise_radius_dependent_arrays(radius);

  int i = 0;
  const_IWSubstring token;

  while (buffer.nextword(token, i, ',')) {
    //  cerr << "build_from_reaction_files:building '" << token << "'\n";
    if (!_build(token, radius)) {
      cerr << "Set_of_Reactions::build_from_reaction_files:cannot process '" << token
           << "'\n";
      return 0;
    }
  }

  return _determine_number_reactions();
}

#ifdef IF_WE_EVER_SUPPORT_SMIRKS
// Not sure if supporting smirks is a good idea of not.
// Many smirks have a lot of information on the LHS to fully specify
// the reagents, but not much information on the RHS. That will not
// work in a reversed reaction context.
int
Set_of_Reactions::build_from_smirks(IWString& fname, int radius)
{
  iwstring_data_source input(fname.null_terminated_chars());
  if (!input.good()) {
    cerr << "Set_of_Reactions::build:cannot open smirks file '" << fname << "'\n";
    return 0;
  }

  return build_from_smirks(input, radius);
}

// Read smirks from `input`. Note that build_rxn_from_smirks adds
// the resulting reactions to _rxn.
int
Set_of_Reactions::build_from_smirks(iwstring_data_source& input, int radius)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!build_rxn_from_smirks(buffer, radius)) {
      cerr << "Set_of_Reactions::build_from_smirks:invalid smirks '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

// Parse `buffer` as a smirks, and if successful add to _rxn.
int
Set_of_Reactions::build_rxn_from_smirks(const const_IWSubstring& buffer, int radius)
{
  if (buffer.empty()) {
    return 1;
  }

  const_IWSubstring smirks, name;
  if (!buffer.split(smirks, ' ', name) || smirks.empty()) {
    cerr << "Set_of_Reactions::build_rxn_from_smirks:no data\n";
  }

  std::unique_ptr<Reaction_With_Stats> rxn = std::make_unique<Reaction_With_Stats>();
  if (!rxn->construct_from_smirks(smirks)) {
    cerr << "Set_of_Reactions::build_rxn_from_smirks:invalid smirks '" << smirks << "'\n";
    return 0;
  }
  rxn->set_comment(name);
  _rxn[radius] << rxn.release();

  return 1;
}
#endif

#ifdef JUST_FOR_DEBUGGING___
static void
set_all_atom_numbers(const ISIS_RXN_FILE_Molecule& m, int* x, const int flag)
{
  const auto matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const auto j = m.atom_map(i);
    x[j] = flag;
  }

  return;
}
#endif

// Put isotopic labels on the atoms and write the array
// and smiles.
void
for_debugging_write_molecule(const Molecule& m, const int* some_array,
                             const char* message)
{
  const int matoms = m.natoms();
  for (int q = 0; q < matoms; ++q) {
    cerr << ' ' << some_array[q];
  }
  cerr << endl;
  Molecule mcopy(m);
  mcopy.set_isotopes(some_array);
  set_include_atom_map_with_smiles(0);
  cerr << message << mcopy.smiles() << "\n";
  set_include_atom_map_with_smiles(1);
  return;
}

// TODO ianwatson: use c++ filesystem
const_IWSubstring
FileName(const IWString& path_name)
{
  const_IWSubstring result(path_name);
  const int last_dir_separator = result.rindex('/');
  if (last_dir_separator > 0) {
    result += (last_dir_separator + 1);
  }
  return result;
}

int
Set_of_Reactions::_build(const const_IWSubstring& buffer,
                         const resizable_array<int>& radius)
{
  // cerr << "Set_of_Reactions::_build: '" << buffer << endl << "'
  // input_is_reaction_smiles " << input_is_reaction_smiles << endl;
  RXN_File rxn;
  rxn.set_auto_fix_orphans(1);

  if (input_is_reaction_smiles) {
    if (!rxn.build_from_reaction_smiles(buffer, 1)) {
      cerr << "Set_of_Reactions::_build:invalid reaction smiles " << buffer << "\n";
      return 0;
    }
  } else if (!rxn.do_read(buffer)) {
    cerr << "read_reaction:cannot read reaction '" << buffer << "'\n";
    return 0;
  }

  if (rxn_name_header != "") {
    IWString newName = rxn_name_header + rxn.name();
    rxn.set_name(newName);
  }

  rxn.check_for_widows_and_orphans();
  if (remove_widows) {
    rxn.remove_widow_atoms();
  }

  // const int initial_nr = rxn.number_reagents();

#ifdef STUFF_NOW_DONE_BY_RXN2RXNSMILES
  if (rxn.remove_fragments_not_participating() && verbose > 1) {
    cerr << "Removed fragments not reacting from '" << rxn.name() << ", now "
         << rxn.number_reagents() << " reagents\n";
  }

  if (rxn.eliminate_reagents_not_participating() && verbose > 1) {
    cerr << "Removed reagents not reacting from '" << rxn.name() << "', now "
         << rxn.number_reagents() << " reagents\n";
  }

  if (rxn.remove_non_participating_fragments() && verbose > 1) {
    cerr << "Removed one or more non participating fragments from '" << rxn.name()
         << "\n";
  }
#endif

  // probably need to turn this on as well
  // if (remove_duplicate_reagents_atom_maps_scrambled &&
  // rxn.remove_duplicate_reagents_atom_maps_scrambled())
  /// cerr << "Removed duplicate reagents, but with different atom maps from '" <<
  /// rxn.name() << "' now " << rxn.number_reagents() << " reagents\n";

  // rxn.do_write("FOO.rxn");

  if (1 == rxn.number_reagents()) {
    ;
  } else if (take_first_of_multiple_reagents) {
    ;
  } else if (rxn.all_reagents_the_same()) {
    cerr << "Reaction " << rxn.name() << " all reagents the same\n";
  } else if (skip_reactions_with_multiple_ragents) {
    if (verbose) {
      cerr << "Set_of_Reactions::_build:skipping multi reagent reaction '" << rxn.name()
           << "' in " << _name << endl;
    }
    return 1;
  } else {
    cerr << "multiple reagents in '" << rxn.name() << "', cannot process\n";
    return 0;
  }

  if (preserve_kekule_forms) {
    rxn.set_preserve_kekule_forms(1);
  } else {
    rxn.set_aromatic_bonds_lose_kekule_identity(1);
  }

  rxn.set_do_automatic_atom_mapping(0);
  // rxn.set_auto_fix_orphans(1);

  if (input_is_reaction_smiles) {
    ;
  } else if (reaction_name_is_path_name) {
    rxn.set_name(rxn_name_header + buffer);
  } else if (reaction_name_is_file_name) {
    const const_IWSubstring fname = FileName(buffer);
    rxn.set_name(rxn_name_header + fname);
  } else if (rxn.name().nwords() > 1) {
    IWString tmp(rxn.name());
    tmp.truncate_at_first(' ');
    rxn.set_name(tmp);
  }

  if (discard_existing_atom_map) {
    rxn.discard_atom_map();
  }

  const int highest_atom_map_number = rxn.highest_atom_map_number();

// #define DEBUG_RETROSYNTHETIC_Q
#ifdef DEBUG_RETROSYNTHETIC_Q
  cerr << __LINE__ << " line\n";
  rxn.debug_print();
#endif

  if (!rxn.prepare_for_reaction_construction()) {
    cerr << "prepare_for_reaction_construction failed '" << rxn.name() << "'\n";
    return 0;
  }

  // set up the controls for how to create the rxn query

  Molecule_to_Query_Specifications mqs;
  mqs.set_set_element_hits_needed_during_molecule_to_query(0);
  mqs.set_aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic(
      aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic);
  mqs.set_atoms_conserve_ring_membership(atoms_conserve_ring_membership);
  mqs.set_bonds_preserve_ring_membership(1);

  if (atoms_conserve_ring_membership) {
    mqs.set_non_ring_atoms_become_nrings_0(1);
  }
  mqs.set_preserve_saturation(preserve_saturation);
  if (preserve_ring_size) {
    mqs.set_preserve_smallest_ring_size(1);
  }

  int* include_atom_map = new int[highest_atom_map_number + 1];
  std::unique_ptr<int[]> free_include_atom_map(include_atom_map);

  for (const auto r : radius) {
    std::fill_n(include_atom_map, highest_atom_map_number + 1, 0);

    for (int reagentIndex = 0; reagentIndex < rxn.number_reagents(); ++reagentIndex) {
      auto& thisReagent = rxn.reagent(reagentIndex);

      // preprocess(thisReagent);    no, we preprocess the molecules we are reacting,
      // doing it here would be dangerous...

      const int matoms = thisReagent.natoms();

      int* changed = new_int(matoms + matoms);
      std::unique_ptr<int[]> free_changed(changed);

      const int c = identify_changing_atoms_reagent(rxn, reagentIndex, changed);

      if (verbose > 1) {
        cerr << "Reaction '" << rxn.name() << "' has " << c << " changing atoms\n";
      }

// #define ECHO_CHANGED 1
#ifdef ECHO_CHANGED
      for_debugging_write_molecule(thisReagent, changed, "Initial Changed: ");
#endif

      if (0 == c &&
          reagentIndex ==
              0)  // after the first one, some reagents may not have any changing atoms
      {
        cerr << "Reaction '" << thisReagent.name()
             << "' has no changing atoms, cannot process\n";
        return ignore_reactions_with_no_changing_atoms;
      }

      for (int i = 0; i < matoms; ++i) {
        auto* a = thisReagent.mdl_atom_data(i);
        a->set_exact_change(1);
      }

      identify_atoms_within_range(thisReagent, r, changed);  // updates the changed array

#ifdef ECHO_CHANGED
      for_debugging_write_molecule(thisReagent, changed, "Changed after radius applied:");
#endif

      for (int i = 0; i < matoms; ++i) {
        //    cerr << " atom " << i << " amap " << thisReagent.atom_map(i) << " " <<
        //    thisReagent.smarts_equivalent_for_atom(i) << " changed " << changed[i] <<
        //    endl;

        if (0 == changed[i]) {
          continue;
        }

        const int amap = thisReagent.atom_map(i);
        // if (amap <= 0)
        if (amap < 0) {
          cerr << "missing amap for atom " << i << endl;
        }

        thisReagent.set_substitution(i, thisReagent.ncon(i));
        //      if (1 == changed[i])
        //        thisReagent.set_substitution(i, thisReagent.ncon(i));
        //      else
        //      {
        //        thisReagent.set_substitution(i, -3);   // IAW extention to substitution,
        //        that means ignore thisReagent.set_ring_bond(i, -3);     // IAW extention
        //        to ring bond count, ignore
        //      }

        include_atom_map[amap] = 1;
      }

      //   for (int i = 0; i < matoms; ++i)
      //   {
      //     if (changed[i])
      //       thisReagent.set_substitution(i, 0);
      //   }
      //

#ifdef ECHO_CHANGED
      cerr << "atoms maps in frag:";
      for (int q = 1; q < highest_atom_map_number + 1; ++q) {
        cerr << " " << include_atom_map[q];
      }
      cerr << endl;
#endif
    }

    // if we have added orphans, they will be additional reagents, make sure to include
    // them

    //  for (int i = 1; i < rxn.number_reagents(); ++i)
    //  {
    //    const auto & y = rxn.reagent(i);
    //    set_all_atom_numbers(y, include_atom_map, 1);
    //  }

#ifdef ECHO_CHANGED
    iw_write_array(include_atom_map, highest_atom_map_number,
                   "atoms maps in frag (final):");
#endif

    for (int reagentIndex = 0; reagentIndex < rxn.number_reagents(); ++reagentIndex) {
      ISIS_RXN_FILE_Molecule& thisMol = rxn.reagent(reagentIndex);
      atype_t* atype2 = new atype_t[thisMol.natoms()];
      std::unique_ptr<atype_t[]> free_atype2(atype2);

      atom_typing_specification.assign_atom_types(thisMol, atype2);

      for (int atomIndex = 0; atomIndex != thisMol.natoms(); ++atomIndex) {
        thisMol.set_userAtomType(atomIndex, atype2[atomIndex]);
      }
      // cerr << "thisMol (reactant): " << thisMol << endl;
    }

    for (int productIndex = 0; productIndex < rxn.number_products(); ++productIndex) {
      ISIS_RXN_FILE_Molecule& thisMol = rxn.product(productIndex);
      atype_t* atype2 = new atype_t[thisMol.natoms()];
      std::unique_ptr<atype_t[]> free_atype2(atype2);

      atom_typing_specification.assign_atom_types(thisMol, atype2);

      for (int atomIndex = 0; atomIndex != thisMol.natoms(); ++atomIndex) {
        thisMol.set_userAtomType(atomIndex, atype2[atomIndex]);
      }
      // cerr << "thisMol (Product): " << *thisMol << endl;
    }

    RXN_File_Create_Reaction_Options rxnfcro;
    rxnfcro.set_only_create_query_from_first_reagent(1);

    rxn.setQueryOutStream(queryEchoStreamPtr);

    rxn.set_auto_fix_orphans(1);

    std::unique_ptr<Reaction_With_Stats> t = std::make_unique<Reaction_With_Stats>();

    if (!rxn.create_reaction(*t, rxnfcro, mqs, include_atom_map)) {
      cerr << "Reaction " << rxn.name() << " cannot create reaction object\n";
      return ignore_bad_reactions;
    }

    t->set_do_not_perceive_symmetry_equivalent_matches(1);

    _rxn[r].add(t.release());

    _id_to_ndx[r][rxn.name()] = _rxn[r].number_elements() - 1;
  }

  return 1;
}

int
Set_of_Reactions::process(Molecule& m, Molecule_to_Match& target,
                          IWString_and_File_Descriptor& output)
{
  int* already_reacted = new_int(_number_reactions);
  std::unique_ptr<int[]> free_already_reacted(already_reacted);

  for (int r = _max_radius; r >= 0; --r) {
    if (_rxn[r].empty()) {
      continue;
    }

    if (_process(m, target, already_reacted, r, output)) {
      _molecules_deconstructed_at_radius[r]++;
      return 1;
    }
  }

  return 0;
}

int
Set_of_Reactions::process(Molecule& m, Molecule_to_Match& target, const int radius,
                          IWString_and_File_Descriptor& output)
{
  const int nr = _rxn[radius].number_elements();

  int rc = 0;

  IW_STL_Hash_Set all_smiles_formed;

  for (int i = 0; i < nr; ++i) {
    IW_STL_Hash_Set usmi_here;

    if (_process(m, target, radius, _rxn[radius][i], usmi_here, output)) {
      _reportResultDuplicates(_rxn[radius][i]->name(), m.name(), usmi_here,
                              all_smiles_formed);

      _molecules_deconstructed_at_radius[radius]++;
      rc++;
      if (break_after_first_deconstruction) {
        break;
      }
    }
  }

  return rc;
}

void
Set_of_Reactions::_reportResultDuplicates(const IWString& nameOfReaction,
                                          const IWString& molName,
                                          const IW_STL_Hash_Set& newResultSmi,
                                          IW_STL_Hash_Set& all_smiles_formed)
{
  for (const auto& smiles : newResultSmi) {
    if (all_smiles_formed.contains(smiles)) {
      cerr << "Reaction " << nameOfReaction << " returned a duplicate result for "
           << molName << endl;
    } else {
      all_smiles_formed.insert(smiles);
    }
  }

  return;
}

int
Set_of_Reactions::_process(Molecule& m, Molecule_to_Match& target, int* already_reacted,
                           const int radius, IWString_and_File_Descriptor& output)
{
  const int number_reactions = _rxn[radius].number_elements();

  int rc = 0;
  IW_STL_Hash_Set all_smiles_formed;

  for (int i = 0; i < number_reactions; ++i) {
    if (already_reacted[i]) {
      continue;
    }

    if (_rxn[radius][i]->matches_found() >= max_matches_each_reaction) {  // testing
      continue;
    }

    IW_STL_Hash_Set newUniqResultSmi;

    if (_process(m, target, radius, _rxn[radius][i], newUniqResultSmi, output)) {
      already_reacted[i] = 1;

      _reportResultDuplicates(_rxn[radius][i]->name(), m.name(), newUniqResultSmi,
                              all_smiles_formed);

      rc++;
      if (break_after_first_deconstruction) {
        return 1;
      }
    }
  }

  return rc;
}

// Set smiles and name in `proto` based on `parent`.
void
SetParentAttributes(Molecule& parent, retros::Output& proto)
{
  const IWString& usmi = parent.unique_smiles();
  proto.set_smiles(usmi.data(), usmi.length());
  const IWString& name = parent.name();
  proto.set_name(name.data(), name.length());
}

void
SetSmiles(const IWString& usmi, retros::Result& result)
{
  result.add_smiles(usmi.data(), usmi.length());
}

void
SetRxnName(const IWString& rxn_name, retros::Result& result)
{
  if (rxn_name.empty()) {
    return;
  }
  result.set_reaction(rxn_name.data(), rxn_name.length());
}

void
SetReactionFamily(const IWString& rxn_family, retros::Result& result)
{
  if (rxn_family.empty()) {
    return;
  }

  result.set_reaction_family(rxn_family.data(), rxn_family.length());
}

void
PlaceTrimmedForm(const IWString& usmi, retros::Output& result)
{
  const auto n = result.result_size();
  if (n == 0) {
    cerr << "PlaceTrimmedForm:no results in proto\n";
    return;
  }

  retros::Result* res = result.mutable_result(n - 1);
  res->set_small_fragments_removed(usmi.data(), usmi.length());
}

// Append `smi` to the latest result in `proto`.
void
AppendProductSmiles(const IWString& smi, retros::Output& proto)
{
  const uint32_t nresult = proto.result_size();
  if (nresult == 0) {
    return;
  }

  retros::Result* res = proto.mutable_result(nresult - 1);
  res->add_smiles(smi.data(), smi.length());
}

// Write `smiles` followed by a text_format representation
// of `proto` to `output`.
int
WriteProto(const IWString& smiles, const retros::Output& proto,
           IWString_and_File_Descriptor& output)
{
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (!printer.PrintToString(proto, &buffer)) {
    cerr << "WriteProto:cannot print " << proto.ShortDebugString() << '\n';
    return 0;
  }

  // Don't write the smiles, too confusing.
  // output << smiles << ' ' << buffer << '\n';
  output << buffer << '\n';
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
do_output(Molecule& parent, Molecule& result, const IWString& reaction_family,
          Reaction_With_Stats& rxn, const int radius,
          IWString_and_File_Descriptor& output)
{
  IWString desc;

  retros::Output proto;

  desc << parent.name() << " via " << rxn.comment();
  if (reaction_family.length() && !regressionFlag) {
    desc << ' ' << reaction_family;
  }
  desc << " R " << radius;

  if (reactantFilenameBase != "") {
    if (!reactantsStream.is_open()) {
      if (!reactantsStream.open(
              (reactantFilenameBase + "_all.smi").null_terminated_chars())) {
        cerr << "Cannot open reactantsStream '" << reactantFilenameBase << "_all.smi"
             << "'\n";
        throw Error("Cannot open reactantsStream");
      }
    }
  }

  if (write_results_as_proto) {
    SetParentAttributes(parent, proto);
    retros::Result* res = proto.add_result();
    SetSmiles(result.unique_smiles(), *res);
    res->set_radius(radius);
    SetRxnName(rxn.comment(), *res);
    SetReactionFamily(reaction_family, *res);
  } else {
    output << result.unique_smiles() << ' ' << desc << " ALL\n";
  }

  const int nf = result.number_fragments();

  IWString thisUniqSmiles =
      "";  // TODO ianwatson:initialisation not needed, here and below?
  if (!add_output_record_with_small_fragments_removed) {
    ;
  } else if (nf > 1) {
    const int matoms = result.natoms();

    int* remove_atom = new_int(matoms);
    std::unique_ptr<int[]> free_remove_atom(remove_atom);

    int nremove = 0;
    for (int i = 0; i < matoms; ++i) {
      const int f = result.fragment_membership(i);

      if (result.atoms_in_fragment(f) >= lower_atom_count_products) {
        continue;
      }

      remove_atom[i] = 1;
      nremove++;
    }

    if (nremove) {
      Molecule mcopy(result);
      mcopy.remove_atoms(remove_atom);

      thisUniqSmiles = mcopy.unique_smiles();
      if (write_results_as_proto) {
        PlaceTrimmedForm(thisUniqSmiles, proto);
      } else {
        output << thisUniqSmiles << ' ' << desc << " SPFRM.1\n";
      }
    } else {
      thisUniqSmiles = result.unique_smiles();
      if (!write_results_as_proto) {
        output << thisUniqSmiles << ' ' << desc << " SPFRM.0\n";
      }
    }
  } else {
    thisUniqSmiles = result.unique_smiles();
    if (!write_results_as_proto) {
      output << thisUniqSmiles << ' ' << desc << " SPFRM.0\n";
    }
  }

  if (nf > 1 && split_products_into_separate_molecules) {
    resizable_array_p<Molecule> f;
    result.create_components(f);

    IW_STL_Hash_Set seen;
    for (int i = 0; i < nf; ++i) {
      if (f[i]->natoms() < lower_atom_count_products) {
        continue;
      }

      IWString partSmiles = f[i]->unique_smiles();
      if (seen.contains(partSmiles)) {
        continue;
      }
      // TODO ianwatson: update seen?

      if (write_results_as_proto) {
        AppendProductSmiles(partSmiles, proto);
      } else {
        output << partSmiles << ' ' << desc << '\n';
      }
    }
  } else if (!write_results_as_proto) {
    output << thisUniqSmiles << ' ' << desc << '\n';
  }

  if (write_results_as_proto) {
    WriteProto(result.unique_smiles(), proto, output);
  }

  IWString reactant1 = "";
  IWString reactant2 = "";

  if (reactantFilenameBase != "") {
    IW_STL_Hash_Set seen;

    // split the parts and record the reactants

    iwaray<const_IWSubstring> reactantSmilesParts;
    int partCount = thisUniqSmiles.split(reactantSmilesParts, '.');

    for (int i = 0; i < partCount; ++i) {
      const_IWSubstring thisPart = reactantSmilesParts[i];

      if (seen.contains(thisPart)) {
        continue;
      }
      reactantsStream << thisPart << '\n';
      if (reactant1 == "") {
        reactant1 = thisPart;
      } else if (reactant2 == "") {
        reactant2 = thisPart;
      }

      seen.insert(thisPart);
    }

    if (reactant2 != "") {
      if (!reactant1Stream.is_open()) {
        if (!reactant1Stream.open(
                (reactantFilenameBase + "_R1.smi").null_terminated_chars())) {
          cerr << "Cannot open reactant1Stream '" << reactantFilenameBase << "_R1.smi"
               << "'\n";
          throw Error("Cannot open reactant1Stream");
        }
        if (!reactant2Stream.open(
                (reactantFilenameBase + "_R2.smi").null_terminated_chars())) {
          cerr << "Cannot open reactant2Stream '" << reactantFilenameBase + "_R2.smi"
               << "'\n";
          throw Error("Cannot open reactant2Stream");
        }
      }

      reactant1Stream << reactant1 << "\n";
      reactant2Stream << reactant2 << "\n";
    } else if (reactant1 != "") {
      if (!reactantOnly1Stream.is_open()) {
        if (!reactantOnly1Stream.open(
                (reactantFilenameBase + "_Only1.smi").null_terminated_chars())) {
          cerr << "Cannot open reactantOnly1Stream '"
               << reactantFilenameBase + "_Only1.smi"
               << "'\n";
          throw Error("Cannot open reactantOnly1Stream");
        }
      }
      reactantOnly1Stream << reactant1 << "\n";
    }
  }

  return 1;
}

static void
do_echo_initial_molecule_before_each_product(Molecule& m,
                                             IWString_and_File_Descriptor& output)
{
  if (!write_results_as_proto) {
    output << m.smiles() << ' ' << m.name() << " PARENT\n";
  }

  return;
}

/*
  Molecule was found to have multiple scaffold embeddings
*/

static int
handle_molecules_skipped_for_multiple_scaffold_embeddings(Molecule& m)
{
  molecules_skipped_for_multiple_scaffold_embeddings++;
  if (stream_for_molecules_with_multiple_scaffold_embeddings.is_open()) {
    stream_for_molecules_with_multiple_scaffold_embeddings << m.smiles() << ' '
                                                           << m.name() << '\n';
    stream_for_molecules_with_multiple_scaffold_embeddings
        .write_if_buffer_holds_more_than(4096);
  }

  return 1;  // always
}

int
Set_of_Reactions::_process(Molecule& m, Molecule_to_Match& target, const int radius,
                           Reaction_With_Stats* rxn, IW_STL_Hash_Set& usmi_here,
                           IWString_and_File_Descriptor& output)
{
  Substructure_Results sresults;

  // cerr << "Reaction has " << rxn->number_sidechains() << " sidechains\n";
  int nhits = rxn->substructure_search(target, sresults);

  acc_nhits[nhits]++;

  rxn->another_search();
  if (verbose > 1) {
    cerr << "at Radius " << radius << " got " << nhits << " hits with " << rxn->comment()
         << " in " << m.name() << ", matched " << rxn->max_query_atoms_matched_in_search()
         << " atoms\n";
  }

  if (0 == nhits) {
    return 0;
  }

  if (verbose > 1 && !multi_threaded) {
    cerr << rxn->comment() << " has " << nhits << " hits in " << m.name() << " "
         << m.smiles() << " radius " << radius << ' ' << _name << endl;
  }

  if (nhits != sresults.number_embeddings()) {
    cerr << "HUH, nhits " << nhits << " number_embeddings "
         << sresults.number_embeddings() << " mismatch, " << rxn->comment()
         << " searching " << m.name() << endl;
    nhits = sresults.number_embeddings();
    //  abort();
  }

  if (nhits > 1 && skip_molecules_with_multiple_scaffold_embeddings) {
    return handle_molecules_skipped_for_multiple_scaffold_embeddings(m);
  }

  rxn->another_match();

  _molecules_deconstructed_at_radius[radius]++;

  for (int i = 0; i < nhits; ++i) {
    Molecule result;
    //  cerr << " processing embedding " << *sresults.embedding(i) << endl;
    if (!rxn->perform_reaction(&m, sresults.embedding(i), result)) {
      cerr << "Reacting " << m.name() << " via " << rxn->comment() << " failed\n";
      return 0;
    }

    const IWString& thisSmi = result.unique_smiles();
    if (usmi_here.contains(thisSmi)) {
      continue;
    }

    usmi_here.insert(thisSmi);

    if (echo_initial_molecule_before_each_product) {
      do_echo_initial_molecule_before_each_product(m, output);
    }

    if (strip_product_to_largest_fragment) {
      result.reduce_to_largest_fragment_carefully();
    }

    if (unfix_implicit_hydrogens_in_products) {
      result.unset_unnecessary_implicit_hydrogens_known_values();
    }

    do_output(m, result, _name, *rxn, radius, output);
  }

  return nhits;
}

Reaction_With_Stats*
Set_of_Reactions::reaction_with_id(const IWString& id, const int radius) const
{
  const auto f = _id_to_ndx[radius].find(id);

  if (f == _id_to_ndx[radius].end()) {
    return nullptr;
  }

  return _rxn[radius][f->second];
}

int
Set_of_Reactions::report() const
{
  if (regressionFlag) {
    cerr << "Set_of_Reactions: "
         << " with " << _number_reactions << " reactions\n";
  } else {
    cerr << "Set_of_Reactions: " << _name << " with " << _number_reactions
         << " reactions\n";
  }

  int molecules_deconstructed = 0;

  for (int i = 0; i <= _max_radius; ++i) {
    if (_rxn[i].empty()) {
      continue;
    }

    cerr << _molecules_deconstructed_at_radius[i] << " molecules deconstructed at radius "
         << i << '\n';
    molecules_deconstructed += _molecules_deconstructed_at_radius[i];
  }

  cerr << molecules_deconstructed << " molecules deconstructed\n";

  return 1;
}

int
Set_of_Reactions::reaction_statistics() const
{
  if (regressionFlag) {
    cerr << "Set_of_Reactions: "
         << " with " << _number_reactions << " reactions\n";
  } else {
    cerr << "Set_of_Reactions: " << _name << " with " << _number_reactions
         << " reactions\n";
  }

  int molecules_deconstructed = 0;

  for (int i = 0; i <= _max_radius; ++i) {
    const int nri = _rxn[i].number_elements();

    if (0 == nri) {
      continue;
    }

    for (int j = 0; j < nri; ++j) {
      const Reaction_With_Stats* r = _rxn[i][j];

      cerr << ' ' << i << ' ' << r->comment() << ' ' << r->searches_done()
           << " searches, " << r->matches_found() << " matches found\n";
    }

    cerr << _molecules_deconstructed_at_radius[i] << " molecules deconstructed at radius "
         << i << '\n';
    molecules_deconstructed += _molecules_deconstructed_at_radius[i];
  }

  cerr << molecules_deconstructed << " molecules deconstructed\n";

  return 1;
}

int
Set_of_Reactions::WriteReactionProto(IWString& fname, int radius) const
{
  iw_tf_data_record::TFDataWriter writer;
  if (!writer.Open(fname)) {
    cerr << "Set_of_Reactions::WriteReactionProto:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteReactionProto(writer, radius);
}

int
Set_of_Reactions::WriteReactionProto(iw_tf_data_record::TFDataWriter& writer,
                                     int radius) const
{
  for (const Reaction_With_Stats* rxn : _rxn[radius]) {
    cerr << "Set_of_Reactions::WriteReactionProto:write reaction not implemented "
         << rxn->comment() << '\n';
  }

  return 1;
}

int
Set_of_Reactions::WriteLhsQueryProto(IWString& name_stem, int radius) const
{
  // as part of debugging also generate a .txtproto output.
  // turn this off once we are convinced the binary protos work.
  IWString txtproto_fname;
  txtproto_fname << name_stem << '_' << _name << '.' << radius << ".txtproto";
  IWString_and_File_Descriptor txtproto_output;
  if (!txtproto_output.open(txtproto_fname.null_terminated_chars())) {
    cerr << "Set_of_Reactions::WriteLhsQueryProto:cannot open '" << txtproto_fname
         << "'\n";
    return 0;
  }

  for (auto* rxn : _rxn[radius]) {
    Substructure_Query& q = *rxn;
    SubstructureSearch::SubstructureQuery proto = q.BuildProto();
    std::string buffer;
    google::protobuf::TextFormat::PrintToString(proto, &buffer);
    txtproto_output << buffer;
  }

  IWString fname;
  fname << name_stem << '_' << _name << '.' << radius << ".dat";
  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Set_of_Reactions::WriteLhsQueryProto:cannot open '" << fname << "'\n";
    return 0;
  }
  iw_tf_data_record::TFDataWriter writer;
  if (!writer.Open(fname)) {
    cerr << "Set_of_Reactions::WriteLhsQueryProto:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteLhsQueryProto(writer, radius);
}

int
Set_of_Reactions::WriteLhsQueryProto(iw_tf_data_record::TFDataWriter& writer,
                                     int radius) const
{
  cerr << "Set_of_Reactions::WriteLhsQueryProto:radius " << radius << " has "
       << _rxn[radius].size() << " reactions\n";
  for (const Reaction_With_Stats* rxn : _rxn[radius]) {
    const Substructure_Query& qry = *rxn;
    SubstructureSearch::SubstructureQuery proto = qry.BuildProto();
    if (!writer.WriteSerializedProto<SubstructureSearch::SubstructureQuery>(proto)) {
      cerr << "Set_of_Reactions::WriteLhsQueryProto:cannot write\n";
      return 0;
    }
  }

  return 1;
}

/*
  Partof retrosynthetic pipeline.
*/

static int
retrosynthesis(Molecule& m, const int* atom_map,
               const resizable_array_p<Set_of_Reactions>& rxn, const int max_radius,
               IWString_and_File_Descriptor& output)
{
  atype_t* atype = new atype_t[m.natoms()];
  std::unique_ptr<atype_t[]> free_atype(atype);

  atom_typing_specification.assign_atom_types(m, atype);

  for (int atomIndex = 0; atomIndex != m.natoms(); ++atomIndex) {
    m.set_userAtomType(atomIndex, atype[atomIndex]);
  }

  Molecule_to_Match target(&m);

  const int nr = rxn.number_elements();

  int* deconstucted_by = new_int(nr);
  std::unique_ptr<int[]> free_deconstucted_by(deconstucted_by);

  int this_molecule_deconstructed = 0;
  for (int r = max_radius; r >= 0; --r) {
    int deconstructed_at_r = 0;

    for (int j = 0; j < nr; ++j) {
      if (deconstucted_by[j]) {
        continue;
      }

      if (!rxn[j]->process(m, target, r, output)) {
        continue;
      }

      deconstucted_by[j] = 1;

      this_molecule_deconstructed++;

      deconstructed_at_r++;

      if (break_after_first_deconstruction) {
        break;
      }
    }

    if (deconstructed_at_r) {
      molecules_deconstructed_at_radius[r]++;
      deconstructions_at_radius[r] += deconstructed_at_r;
    }
  }

  if (this_molecule_deconstructed) {
    molecules_deconstructed++;
    return 1;
  }

  return handle_non_reacting_molecule(m);
}

static int
contains_atom_map_values(const Molecule& m)
{
  int istop = m.natoms();
  if (istop > 5) {  // arbitrary choice
    istop = 5;
  }

  for (int i = 0; i < istop; ++i) {
    if (m.atomi(i)->atom_map()) {
      return 1;
    }
  }

  return 0;
}

static int
retrosynthesis(Molecule& m, const resizable_array_p<Set_of_Reactions>& rxn,
               const int max_radius, IWString_and_File_Descriptor& output)
{
  if (!contains_atom_map_values(m)) {
    return retrosynthesis(m, nullptr, rxn, max_radius, output);
  }

  const int matoms = m.natoms();

  int* atom_map = new int[matoms];
  std::unique_ptr<int[]> free_atom_map(atom_map);

  for (int i = 0; i < matoms; ++i) {
    atom_map[i] = m.atomi(i)->atom_map();
    m.set_atom_map_number(i, 0);
  }

  m.unset_unnecessary_implicit_hydrogens_known_values();

  return retrosynthesis(m, atom_map, rxn, max_radius, output);
}

static int
retrosynthesis(data_source_and_type<Molecule>& input,
               const resizable_array_p<Set_of_Reactions>& rxn, const int max_radius,
               IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!retrosynthesis(*m, rxn, max_radius, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);

    if (queryEchoStreamPtr != nullptr && queryEchoStreamPtr->is_open()) {
      queryEchoStreamPtr->close();
      delete queryEchoStreamPtr;
      queryEchoStreamPtr = nullptr;
    }

    if (report_progress()) {
      cerr << "Read " << molecules_read << " molecules, " << molecules_deconstructed
           << " deconstructed, " << molecules_produced << " molecules produced\n";
    }
  }

  return 1;
}

static int
read_reactions_buffer(const const_IWSubstring& s, const resizable_array<int>& radius,
                      resizable_array_p<Set_of_Reactions>& rxn)
{
  IWString fname, name;

  if (s.contains(' ')) {
    int i = 0;
    s.nextword(fname, i);
    s.nextword(name, i);

  } else {
    fname = s;
  }

  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "read_reactions:cannot open '" << fname << "'\n";
    return 0;
  }

  std::unique_ptr<Set_of_Reactions> r = std::make_unique<Set_of_Reactions>();

  r->set_name(rxn_name_header + name);

  if (!r->build(input, radius)) {
    cerr << "Cannot read reactions '" << s << "'\n";
    return 0;
  }

  rxn << r.release();

  return rxn.number_elements();
}

static int
read_reactions(iwstring_data_source& input, const resizable_array<int>& radius,
               resizable_array_p<Set_of_Reactions>& rxn)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!read_reactions_buffer(buffer, radius, rxn)) {
      cerr << "Fatal error processing reaction specification '" << buffer << "'\n";
      return 0;
    }
  }

  return rxn.number_elements();
}

static int
read_reactions(const const_IWSubstring& s, const resizable_array<int>& radius,
               resizable_array_p<Set_of_Reactions>& rxn)
{
  iwstring_data_source input(s);

  if (!input.good()) {
    cerr << "Cannot open '" << s << "'\n";
    return 0;
  }

  return read_reactions(input, radius, rxn);
}

static int
molecule_found(Molecule& needle, const resizable_array_p<Molecule>& frags)
{
  const IWString& s = needle.unique_smiles();

  const int nfrags = frags.number_elements();

  for (int i = 0; i < nfrags; ++i) {
    if (s == frags[i]->unique_smiles()) {
      return 1;
    }
  }

  return 0;
}

static int
run_self_react_test(Molecule& m, Reaction_With_Stats& r, const int max_radius,
                    int& nfailures,
                    const IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  Substructure_Results sresults;

  const int nhits = r.substructure_search(m, sresults);
  acc_nhits[nhits]++;
  if (nhits != sresults.number_embeddings()) {
    cerr << "HUH " << nhits << " vs " << sresults.number_embeddings() << endl;
  }

  if (0 == nhits) {
    cerr << "Molecule " << m.name() << " no substructure match to '" << r.name()
         << "', only matched " << sresults.max_query_atoms_matched_in_search()
         << " query atoms\n";
    nfailures++;
    return 0;
  }

  // cerr << nhits << " hits with " << m.smiles() << ' ' << m.name() << endl;

  if (0 == test_results.size()) {
    return 1;
  }

  const auto f = test_results.find(m.name());

  if (f == test_results.end()) {
    cerr << "run_self_react_test:no result for " << m.name() << ", ignored\n";
    return 1;
  }

  Molecule* correct = f->second;
  if (product_molecule_chemical_standardisation.active()) {
    product_molecule_chemical_standardisation.process(*correct);
  }

  // cerr << nhits << " hits to " << m.name() << " for reaction " << r.name() << endl;

  resizable_array_p<Molecule> results;

  for (int i = 0; i < nhits; ++i) {
    Molecule* result = new Molecule;

    results.add(result);

    if (!r.perform_reaction(&m, sresults.embedding(i), *result)) {
      cerr << "Reacting " << m.name() << " with " << r.comment() << " failed\n";
      return 0;
    }
    //  cerr << "Result is " << result->smiles() << endl;

    result->reset_all_atom_map_numbers();
    result->unset_unnecessary_implicit_hydrogens_known_values();

    product_molecule_chemical_standardisation.process(*result);

    result->transform_to_non_isotopic_form();

    if (1 == result->number_fragments() &&
        1 == correct->number_fragments())  // the easy case
    {
      if (result->unique_smiles() == correct->unique_smiles()) {
        return 1;
      }

      continue;
    }

    // Now the more difficult case of one or more fragments either side

    if (1 == result->number_fragments()) {
      resizable_array_p<Molecule> frags;
      correct->create_components(frags);

      if (molecule_found(*result, frags)) {
        return 1;
      }

      continue;
    }

    if (1 == correct->number_fragments()) {
      resizable_array_p<Molecule> frags;
      result->create_components(frags);

      if (molecule_found(*correct, frags)) {
        return 1;
      }

      continue;
    }

    // both are multi fragment

    resizable_array_p<Molecule> rfrags, cfrags;
    result->create_components(rfrags);
    correct->create_components(cfrags);

    for (int j = 0; j < rfrags.number_elements(); ++j) {
      //    cerr << "Checking " << rfrags[j]->unique_smiles() << endl;
      if (molecule_found(*rfrags[j], cfrags)) {
        return 1;
      }
    }
  }

  cerr << "Result mismatch, " << m.name() << ' ' << " expected " << correct->smiles()
       << " but produced " << nhits << "\n";
  for (int i = 0; i < results.number_elements(); ++i) {
    cerr << "  " << i << ' ' << results[i]->unique_smiles() << endl;
  }

  nfailures++;
  return 0;
}

static int
run_self_react_test(Molecule& m, const resizable_array_p<Set_of_Reactions>& rxn,
                    const int max_radius, int& nfailures,
                    const IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  const IWString& id = m.name();

  for (int i = 0; i < rxn.number_elements(); ++i) {
    Reaction_With_Stats* r = rxn[i]->reaction_with_id(id, max_radius);

    if (nullptr == r) {
      continue;
    }

    return run_self_react_test(m, *r, max_radius, nfailures, test_results);
  }

  cerr << "run_self_react_test:no reactions for '" << id << "'";
  molecules_with_no_corresponding_reaction++;

  if (ignore_no_matching_reaction_test) {
    cerr << ", ignored\n";
    return 1;
  }

  cerr << endl;

  return 0;
}

static int
run_self_react_test(data_source_and_type<Molecule>& input,
                    const resizable_array_p<Set_of_Reactions>& rxn, const int max_radius,
                    int& nfailures,
                    const IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  int rc = 1;

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (run_self_react_test(*m, rxn, max_radius, nfailures, test_results)) {
      continue;
    }

    rc = 0;

    molecules_failing++;

    if (!keep_going_after_test_failure) {
      return 0;
    }
  }

  return rc;
}

static int
run_self_react_test(const char* fname, FileType input_type,
                    const resizable_array_p<Set_of_Reactions>& rxn, const int max_radius,
                    int& nfailures,
                    const IW_STL_Hash_Map<IWString, Molecule*>& test_results)
{
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return run_self_react_test(input, rxn, max_radius, nfailures, test_results);
}

static int
retrosynthesis(const char* fname, FileType input_type,
               const resizable_array_p<Set_of_Reactions>& rxn, const int max_radius,
               IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return retrosynthesis(input, rxn, max_radius, output);
}

static int
get_single_value(const const_IWSubstring& buffer, resizable_array<int>& radius)
{
  int r;
  if (!buffer.numeric_value(r) || r < 0) {
    return 0;
  }

  radius.add_if_not_already_present(r);

  return 1;
}

static int
get_range_of_values(const const_IWSubstring& buffer, resizable_array<int>& radius)
{
  const_IWSubstring s1, s2;
  if (!buffer.split(s1, '-', s2) || 0 == s1.length() || 0 == s2.length()) {
    return 0;
  }

  int r1, r2;
  if (!s1.numeric_value(r1) || r1 < 0 || !s2.numeric_value(r2) || r2 < r1) {
    return 0;
  }

  for (int r = r1; r <= r2; ++r) {
    radius.add_if_not_already_present(r);
  }

  return 1;
}

static int
gather_radii(const const_IWSubstring& buffer, resizable_array<int>& radius)
{
  int i = 0;
  const_IWSubstring s;
  while (buffer.nextword(s, i, ',')) {
    int rc;

    if (s.contains('-')) {
      rc = get_range_of_values(s, radius);
    } else {
      rc = get_single_value(s, radius);
    }
    if (!rc) {
      cerr << "Invalid radius specification '" << s << "'\n";
      return 0;
    }
  }

  return radius.number_elements();
}

static void
display_misc_dash_X_options()
{
  cerr << " -X ersfrm       add an extra output with small fragments (-a) removed\n";
  cerr << " -X nsmfp        do NOT split multi fragment products into individual "
          "components\n";
  cerr << " -X kekule       preserve Kekule forms\n";
  cerr << " -X test         test self reaction - must provide reaction as both Reaction "
          "and input\n";
  cerr << " -X kg           keep going after a test failure\n";
  cerr << " -X oknr         just skip cases where we cannot find a reaction to match the "
          "name of the molecule\n";
  cerr << " -X maxr=<n>     once a reaction has reacted <n> times, do not try any more. "
          "Just for run times during testing\n";
  cerr << " -X fastexit     immediate exit - do not de-allocate data structures\n";
  cerr << " -X rmiso        remove isotopes from incoming molecules\n";
  cerr << " -X proto        write results in retros::Output proto form\n";
  cerr << " -X rmwidow      remove widow atoms\n";

  exit(1);
}

static void
report_acc_nhits(const extending_resizable_array<int>& acc_nhits)
{
  for (int i = 0; i < acc_nhits.number_elements(); ++i) {
    if (acc_nhits[i]) {
      cerr << acc_nhits[i] << " reactions had " << i << " hits\n";
    }
  }

  return;
}

static int
WriteReactionQueries(const Set_of_Reactions& rxn, int radius, IWString& name_stem)
{
  return rxn.WriteLhsQueryProto(name_stem, radius);
}

void
SetName(const IWString& fname, Set_of_Reactions& rxn)
{
  if (reaction_name_is_path_name) {
    rxn.set_name(fname);
  } else if (reaction_name_is_file_name) {
    const IWString filename = FileName(fname);
    rxn.set_name(filename);
  } else {
    cerr << "Not sure what the name should be\n";
    rxn.set_name(fname);
  }

  cerr << "SetName:name is " << rxn.name() << "'\n";
}

static int
retrosynthesis(int argc, char** argv)
{
  Command_Line cl(argc, argv,
                  "vA:E:i:g:lR:F:P:xc:D:fbG:r:q:Q:U:TSm:a:X:I:Lu:Y:zZSO:o:sn:");

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
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  } else {
    set_auto_create_new_elements(1);
  }

  if (cl.option_present('g')) {
    if (!input_molecule_chemical_standardisation.construct_from_command_line(
            cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('Y')) {
    if (!product_molecule_chemical_standardisation.construct_from_command_line(
            cl, verbose > 1, 'Y')) {
      cerr << "Cannot process chemical standardisation options (-Y)\n";
      usage(32);
    }
  }

  if (cl.option_present('z')) {
    ignore_reactions_with_no_changing_atoms = 1;
    if (verbose) {
      cerr << "Will ignore reactions with no changing atoms\n";
    }
  }

  if (cl.option_present('Z')) {
    ignore_bad_reactions = 1;
    if (verbose) {
      cerr << "Will ignore reactions that cannot be constructed\n";
    }
  }

  if (cl.option_present('S')) {
    regressionFlag = 1;
    if (verbose) {
      cerr << "Will supress filenames in output (for regression testing)\n";
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('L')) {
    strip_product_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce product molecules to largest fragment\n";
    }
  }

  if (cl.option_present('u')) {
    const_IWSubstring u;
    for (int i = 0; cl.value('u', u, i); ++i) {
      if (u.starts_with("prod")) {
        unfix_implicit_hydrogens_in_products = 1;

        if (verbose) {
          cerr << "Will unfix implicit hydrogens in product molecules\n";
        }
      } else if (u.starts_with("reac")) {
        unfix_implicit_hydrogens_in_reagents = 1;

        if (verbose) {
          cerr << "Will unfix implicit hydrogens in input molecules\n";
        }
      } else {
        cerr << "Unrecognised -u qualifier '" << u << "'\n";
        return 1;
      }
    }
  }

  if (cl.option_present('a')) {
    if (!cl.value('a', lower_atom_count_products) || lower_atom_count_products < 0) {
      cerr << "The atom count cutoff for product molecules (-a) must be a whole +ve "
              "number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will not write product fragments with fewer than "
           << lower_atom_count_products << " atoms\n";
    }
  }

  int test_mode = 0;
  int immediate_exit = 0;
  // if we are writing protos of either the reactions or the
  // queries, those names will be in these variables.
  // Note that reactions differ by data, so these file names
  // are stems, and will get the radius and a '.dat' suffix.
  IWString reaction_fname, reaction_queries_fname;

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if ("ersfrm" == x) {
        add_output_record_with_small_fragments_removed = 1;
      } else if ("nsmfp" == x) {
        split_products_into_separate_molecules = 0;
        if (verbose) {
          cerr << "Will NOT split products into separate molecules\n";
        }
      } else if ("kekule" == x || "Kekule" == x) {
        preserve_kekule_forms = 1;
        if (verbose) {
          cerr << "Will preserve Kekule forms\n";
        }
      } else if (x.starts_with("maxr=")) {
        x.remove_leading_chars(5);
        if (!x.numeric_value(max_matches_each_reaction) ||
            max_matches_each_reaction < 1) {
          cerr << "The max number of reactions to perform (marxr=) just be a whole +ve "
                  "number\n";
          return 1;
        }
        if (verbose) {
          cerr << "Once a reaction has performed " << max_matches_each_reaction
               << " matches at each radius, it will no longer try\n";
        }
      } else if ("test" == x) {
        test_mode = 1;
        if (verbose) {
          cerr << "WIll test reactions for ability to react themselves\n";
        }
      } else if ("kg" == x) {
        keep_going_after_test_failure = 1;
        if (verbose) {
          cerr << "Will keep going after an individual test failure\n";
        }
      } else if ("oknr" == x) {
        ignore_no_matching_reaction_test = 1;
        if (verbose) {
          cerr << "Will skip cases where no matching reaction can be found\n";
        }
      } else if (x.starts_with("testp=")) {
        x.remove_leading_chars(6);
        if (!read_test_results(x, test_results)) {
          cerr << "Cannot read test results from '" << x << "'\n";
          return 1;
        }

        if (verbose) {
          cerr << "Read " << test_results.size() << " test results from '" << x << "'\n";
        }
      } else if ("fastexit" == x) {
        immediate_exit = 1;
      } else if ("rmiso" == x) {
        remove_isotopes_from_incoming_molecules = 1;
        if (verbose) {
          cerr << "Will remove isotopes from incoming molecules\n";
        }
      } else if (x == "proto") {
        write_results_as_proto = 1;
        if (verbose) {
          cerr << "Results will be written in proto form\n";
        }
      } else if (x.starts_with("rproto=")) {
        x.remove_leading_chars(7);
        reaction_fname = x;
        if (verbose) {
          cerr << "Reaction protos written to file name stem '" << reaction_fname
               << "'\n";
        }
      } else if (x.starts_with("qproto=")) {
        x.remove_leading_chars(7);
        reaction_queries_fname = x;
        if (verbose) {
          cerr << "Reaction query protos written to '" << reaction_queries_fname << "'\n";
        }
      } else if (x == "rmwidow") {
        remove_widows = 1;
        if (verbose) {
          cerr << "Widow type atoms removed from LHS\n";
        }
      } else if ("help" == x) {
        display_misc_dash_X_options();
      } else {
        cerr << "Unrecognised -X directive '" << x << "'\n";
        display_misc_dash_X_options();
      }
    }
  }

  // msiEchoStreamPtr = new std::ofstream("MsiRxn.smi", std::ofstream::out);

  if (cl.option_present('Q')) {
    IWString queryFilename("");
    cl.value('Q', queryFilename);
    queryEchoStreamPtr = new std::ofstream(queryFilename, std::ofstream::out);

    if (verbose) {
      cerr << "Queries echo'd to the files '" << queryFilename << "'\n";
    }
  }

  if (cl.option_present('n')) {
    cl.value('n', rxn_name_header);

    if (verbose) {
      cerr << "Reaction queries will be saved with names that start with: "
           << rxn_name_header << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialis progress reporting (-r)\n";
      return 1;
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "INvalid atom typing specification '" << p << "'\n";
      return 1;
    }
  } else  // set the default atom typing
  {
    if (!atom_typing_specification.build("UST:AZUCORS")) {
      cerr << "Error setting the defalt UST:AZUCORS atom typing specification\n";
      return 1;
    }
  }

  if (cl.option_present('c')) {
    int c;
    for (int i = 0; cl.value('c', c, i); ++i) {
      changing_atoms_should_be.add(c);
      if (verbose) {
        cerr << "If the number of changing atoms is not " << c
             << " will recompute the atom map\n";
      }
    }
  }

  if (cl.option_present('b')) {
    break_after_first_deconstruction = 1;

    if (verbose) {
      cerr << "Will only do one reaction for any input molecule\n";
    }
  }

  if (cl.option_present('q')) {
    const_IWSubstring q = cl.string_value('q');

    if ('f' == q) {
      reaction_name_is_file_name = 1;
      if (verbose) {
        cerr << "Will use the reaction file name as the reaction name\n";
      }
    } else if ('p' == q) {
      reaction_name_is_path_name = 1;
      if (verbose) {
        cerr << "Will use the reaction path name as the reaction name\n";
      }
    } else {
      cerr << "Unrecognised reaction name modifier '" << q << "'\n";
      usage(1);
    }
  }

  resizable_array<int> radius;

  if (!cl.option_present('R')) {
    cerr << "Must specify one of more radii via the -R option\n";
    usage(1);
  }

  if (cl.option_present('R')) {
    const_IWSubstring r;
    for (int i = 0; cl.value('R', r, i); ++i) {
      if (!gather_radii(r, radius)) {
        cerr << "Invalid radius specification '" << r << "'\n";
        return 1;
      }
    }

    radius.iwqsort_lambda([](int r1, int r2) {
      if (r1 < r2) {
        return -1;
      }
      if (r1 > r2) {
        return 1;
      }
      return 0;
    });
  }

  int max_radius = radius.last_item();

  molecules_deconstructed_at_radius = new_int(max_radius + 1);
  std::unique_ptr<int[]> free_molecules_deconstructed_at_radius(
      molecules_deconstructed_at_radius);
  deconstructions_at_radius = new_int(max_radius + 1);
  std::unique_ptr<int[]> free_deconstructions_at_radius(deconstructions_at_radius);

  resizable_array_p<Set_of_Reactions> rxn;

  // resizable_array_p<Reaction_With_Stats> * rxn = new
  // resizable_array_p<Reaction_With_Stats>[max_radius + 1];
  // std::unique_ptr<resizable_array_p<Reaction_With_Stats>[]> free_rxn(rxn);

  if (!cl.option_present('F') && !cl.option_present('G') && !cl.option_present('T') &&
      !cl.option_present('I')) {
    cerr << "Must specify file of reactions via the -F, -G, -T or -I options\n";
    usage(1);
  }

  const auto t0 = time(nullptr);

  set_iwreaction_display_no_atoms_in_query_message(0);

  if (cl.option_present('F')) {
    const_IWSubstring f;
    for (int i = 0; cl.value('F', f, i); ++i) {
      IWString fname, name;

      if (fname.contains(',')) {
        f.split(fname, ',', name);
      } else {
        fname = f;
      }

      std::unique_ptr<Set_of_Reactions> r = std::make_unique<Set_of_Reactions>();

      if (!r->build(fname, radius)) {
        cerr << "Cannot read reactions from '" << f << "'\n";
        return 1;
      }

      if (verbose) {
        cerr << "read " << r->number_reactions() << " from " << f << '\n';
      }
      r->set_name(rxn_name_header + name);
      rxn << r.release();
    }
  }

  if (cl.option_present('G')) {
    const_IWSubstring g;
    for (int i = 0; cl.value('G', g, i); ++i) {
      if (!read_reactions(g, radius, rxn)) {
        cerr << "Cannot read reactions from '" << g << "'\n";
        return 1;
      }
    }
  }

  if (cl.option_present('T')) {
    const_IWSubstring t;
    for (int i = 0; cl.value('T', t, i); ++i) {
      std::unique_ptr<Set_of_Reactions> r = std::make_unique<Set_of_Reactions>();

      if (!r->build_from_reaction_files(t, radius)) {
        cerr << "Cannot read reactions in '" << t << "'\n";
        return 1;
      }
      SetName(t, *r);
      rxn << r.release();
    }
  }

#ifdef IF_WE_EVER_SUPPORT_SMIRKS
  if (cl.option_present('K')) {
    IWString k;
    for (int i = 0; cl.value('K', k, i); ++i) {
      std::unique_ptr<Set_of_Reactions> r = std::make_unique<Set_of_Reactions>();
      if (!r->build_from_smirks(k, radius)) {
        cerr << "Cannot read smirks in '" << k << "'\n";
        return 1;
      }
      r->set_name(k);
      rxn << r.release();
    }
  }
#endif

  if (cl.option_present('s')) {
    input_is_reaction_smiles = 1;
  }

  if (cl.option_present('I')) {
    input_is_reaction_smiles = 1;

    IWString s;
    for (int i = 0; cl.value('I', s, i); ++i) {
      std::unique_ptr<Set_of_Reactions> r = std::make_unique<Set_of_Reactions>();

      if (!r->build(s.null_terminated_chars(), radius)) {
        cerr << "Cannot read file of reaction smiles '" << s << "'\n";
        return 1;
      }

      SetName(s, *r);
      rxn << r.release();
    }
  }

  if (verbose) {
    const auto t1 = time(NULL);
    cerr << "Reading reactions took " << (t1 - t0) << " seconds\n";
  }

  if (!reaction_queries_fname.empty()) {
    cerr << "have " << rxn.size() << " sets of reactions\n";
    for (const Set_of_Reactions* sor : rxn) {
      if (!WriteReactionQueries(*sor, max_radius, reaction_queries_fname)) {
        cerr << "Cannot write reaction queries\n";
        return 1;
      }
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

  std::unique_ptr<IWString_and_File_Descriptor> output;
  if (cl.option_present('o')) {
    IWString outputFilename;
    cl.value('o', outputFilename);
    output = std::make_unique<IWString_and_File_Descriptor>();
    if (!output->open(outputFilename.null_terminated_chars())) {
      cerr << "Cannot open output file '" << outputFilename << "'\n";
    }

    if (verbose) {
      cerr << "Output written to " << outputFilename << endl;
    }
  } else {
    output = std::make_unique<IWString_and_File_Descriptor>(1);
  }

  if (cl.option_present('O')) {
    cl.value('O', reactantFilenameBase);
  }

  if (cl.option_present('x')) {
    discard_existing_atom_map = 1;

    if (verbose) {
      cerr << "Will discard existing atom map data\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (test_mode) {
    int rc = 0;
    int nfailures = 0;
    for (int i = 0; i < cl.number_elements(); ++i) {
      if (!run_self_react_test(cl[i], input_type, rxn, max_radius, nfailures,
                               test_results)) {
        rc = 1;
      }
    }

    cerr << "encountered " << nfailures << " failures in " << molecules_read
         << " molecules ";
    if (0 == molecules_read) {
      cerr << " 0\n";
    } else {
      cerr << static_cast<float>(nfailures) / static_cast<float>(molecules_read) << endl;
    }

    cerr << molecules_with_no_corresponding_reaction
         << " input molecules had no corresponding reaction (name matching failure)\n";

    if (0 == rc) {
      cerr << "All tests successful\n";
    }

    report_acc_nhits(acc_nhits);

    return rc;
  }

  if (cl.option_present('U')) {
    IWString u = cl.string_value('U');
    if (!u.ends_with(".smi")) {
      u << ".smi";
    }

    if (!stream_for_non_reacting_molecules.open(u.null_terminated_chars())) {
      cerr << "Cannot open stream for non reacting molecules '" << u << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Non reacting molecules written to '" << u << "'\n";
    }
  }

  if (cl.option_present('m')) {
    IWString m = cl.string_value('m');
    if ('.' == m) {
      ;
    } else {
      if (!stream_for_molecules_with_multiple_scaffold_embeddings.open(
              m.null_terminated_chars())) {
        cerr << "Cannot open stream for molecules with multiple scaffold embeddings '"
             << m << "'\n";
        return 1;
      }

      if (verbose) {
        cerr << "MOlecules with multiple scaffold embeddings written to "
             << stream_for_molecules_with_multiple_scaffold_embeddings << endl;
      }
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!retrosynthesis(cl[i], input_type, rxn, max_radius, *output.get())) {
      cerr << "Failure processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  if (queryEchoStreamPtr != nullptr && queryEchoStreamPtr->is_open()) {
    queryEchoStreamPtr->flush();
  }
  if (reactantsStream.is_open()) {
    reactantsStream.flush();
  }
  if (reactant1Stream.is_open()) {
    reactant1Stream.flush();
  }
  if (reactant2Stream.is_open()) {
    reactant2Stream.flush();
  }
  if (reactantOnly1Stream.is_open()) {
    reactantOnly1Stream.flush();
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules, " << molecules_deconstructed
         << " deconstructed\n";
    if (skip_molecules_with_multiple_scaffold_embeddings) {
      cerr << molecules_skipped_for_multiple_scaffold_embeddings
           << " molecules skipped for multiple scaffold substructure matches\n";
    }
    for (int i = 0; i <= max_radius; ++i) {
      if (!radius.contains(i)) {
        continue;
      }

      cerr << molecules_deconstructed_at_radius[i]
           << " molecules deconstructed at radius " << i << endl;
    }

    for (int i = 0; i <= max_radius; ++i) {
      cerr << deconstructions_at_radius[i] << " deconstructions done at radius " << i
           << endl;
    }

    const int nr = rxn.number_elements();

    for (int i = 0; i < nr; ++i) {
      rxn[i]->report();
    }
  }

  if (verbose && molecules_read > 0) {
    for (int i = 0; i < rxn.number_elements(); ++i) {
      rxn[i]->reaction_statistics();
    }

    report_acc_nhits(acc_nhits);
  }

  if (immediate_exit) {
    if (stream_for_molecules_with_multiple_scaffold_embeddings.is_open()) {
      stream_for_molecules_with_multiple_scaffold_embeddings.flush();
    }
    if (stream_for_non_reacting_molecules.is_open()) {
      stream_for_non_reacting_molecules.flush();
    }
    if (queryEchoStreamPtr != nullptr && queryEchoStreamPtr->is_open()) {
      queryEchoStreamPtr->flush();
    }

    if (reactantsStream.is_open()) {
      reactantsStream.flush();
      reactant1Stream.flush();
      reactant2Stream.flush();
      reactantOnly1Stream.flush();
    }

    cerr.flush();

    _exit(0);
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = retrosynthesis(argc, argv);

  return rc;
}
