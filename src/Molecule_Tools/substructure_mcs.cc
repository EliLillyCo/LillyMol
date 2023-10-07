// Look for one set of molecules in another as substructures, by applying
// various molecular transformations.

#include <iostream>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/msi_object.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Lib/reaction.pb.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace substructure_mcs {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on

  cerr << "Does substructure searches between sets of molecules\n";
  cerr << "Takes two arguments, haystack and needles. Task is to find needles in haystack molecules\n";
  cerr << "If matches are not found, transformations are applied until a match is found\n";
  cerr << " -R <fname>   file containing proto form transformation reactions\n";
  cerr << " -r <fname>   file containing msi   form transformation reactions\n";
  cerr << " -a           perform reactions on needles - as well as the haystack\n";
  cerr << " -T e1=d2     element transformations\n";
  cerr << " -G ...       graph reduction specifications, suggest '-G def'\n";
  cerr << " -c           remove all chirality on input\n";
  cerr << " -l           strip to largest fragment\n";
  cerr << " -U <fname>   write unmatched molecules to <fname>\n";
  cerr << " -g ...       chemical standardisation\n";
  cerr << " -s           perform scaffold matches\n";
  cerr << " -y           smiles output (rather than the default tabular)\n";
  cerr << " -i ...       specify input type\n";
  cerr << " -v           verbose output\n";
  ::exit(rc);
}

// A reaction class that remembers its results.
// The tracking of success is not robust.
class ReactionAndStats {
  private:
    IWReaction _rxn;
    int _successful_substructure_searches = 0;
    int _successful_lookups = 0;
  public:
    ReactionAndStats();

    int construct_from_msi_object(msi_object& msi, const Sidechain_Match_Conditions& smc) {
      return _rxn.construct_from_msi_object(msi, smc);
    }

    const IWString& name() const {
      return _rxn.name();
    }

    int ConstructFromProto(ReactionProto::Reaction& proto, IWString& fname);
    int substructure_search(Molecule_to_Match& target, Substructure_Results& sresults);
    int substructure_search(Molecule& m, Substructure_Results& sresults);

    int in_place_transformation(Molecule& m, const Set_of_Atoms* embedding) {
      return _rxn.in_place_transformation(m, embedding);
    }
    int in_place_transformations(Molecule& m, Substructure_Results& sresults)  {
      return _rxn.in_place_transformations(m, sresults);
    }

    void GotSuccessfulLookup() {
      ++_successful_lookups;
    }

    int Report(std::ostream& output) const;
};

ReactionAndStats::ReactionAndStats() {
}

int
ReactionAndStats::ConstructFromProto(ReactionProto::Reaction& proto,
                        IWString& fname) {
  return _rxn.ConstructFromProto(proto, fname);
}

int
ReactionAndStats::substructure_search(Molecule_to_Match& target,
                                      Substructure_Results& sresults) {
  const int nhits = _rxn.substructure_search(target, sresults);
  if (nhits == 0) {
    return 0;
  }

  ++_successful_substructure_searches;
  return nhits;
}

int
ReactionAndStats::substructure_search(Molecule& m,
                                      Substructure_Results& sresults) {
  const int nhits = _rxn.substructure_search(m, sresults);
  if (nhits == 0) {
    return 0;
  }

  ++_successful_substructure_searches;
  return nhits;
}

int
ReactionAndStats::Report(std::ostream& output) const {
  output << _rxn.name() << ' ' << _successful_substructure_searches <<
          " successful substructure searches and " <<
          _successful_lookups << " successful lookups";
  if (_successful_substructure_searches) {
    output << " %";
  }
  if ( _successful_lookups) {
    output << "%";
  }
  output << '\n';
  return 1;
}

enum class RxnFileType {
  kMsi,
  kProto,
};

class Options {
  private:
    int _verbose = 0;

    FileType _input_type = FILE_TYPE_INVALID;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    Element_Transformations _element_transformations;

    Mol2Graph _mol2graph;

    Molecule_to_Query_Specifications _mqs;

    resizable_array_p<ReactionAndStats> _reaction;

    int _check_scaffolds = 0;

    int _molecules_read = 0;

    int _perform_reactions_on_needles = 0;

    IWString_and_File_Descriptor stream_for_unmatched_molecules;

    Sidechain_Match_Conditions _smc;

    int _tabular_output = 1;

    char _output_separator = ' ';

  // private functions

    int ReadReactions(RxnFileType rxn_file_type, IWString& fname);
    int ReadReactions(RxnFileType rxn_file_type, IWString& fname, iwstring_data_source& input);
    int ReadReaction(RxnFileType rxn_file_type, const IWString& fname, const IWString& line);
    int ReadReactionProto(IWString& fname);
    int ReadReactionMsi(IWString& fname);

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int MaybeDiscernInputType(const char * fname);

    FileType input_type() const {
      return _input_type;
    }

    void ReadAnotherMolecule() {
      ++_molecules_read;
    }

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }

    char output_separator() const {
      return _output_separator;
    }

    int perform_reactions_on_needles() const {
      return _perform_reactions_on_needles;
    }

    int check_scaffolds() const {
      return _check_scaffolds;
    }

    int tabular_output() const {
      return _tabular_output;
    }

    Mol2Graph& mol2graph() {
      return _mol2graph;
    }

    Molecule_to_Query_Specifications& mqs() {
      return _mqs;
    }

    resizable_array_p<ReactionAndStats>& reactions() {
      return _reaction;
    }

    int MaybeWriteNonMatch(Molecule& m);
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _perform_reactions_on_needles = 0;
  _molecules_read = 0;
  _tabular_output = 1;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  set_warn_duplicate_atom_deletion(0);
  set_warn_non_bonded_breaks(0);
  set_display_strange_chemistry_messages(0);

  _mqs.set_make_embedding(1);

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }
  _chemical_standardisation.set_verbose(0);

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('G')) {
    if (! _mol2graph.construct(cl, 'G', _verbose)) {
      cerr << "Options::Initialise:cannot initialise mol2graph (-G) specifications\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('a')) {
    _perform_reactions_on_needles = 1;
    if (_verbose) {
      cerr << "Will perform reactions on needles\n";
    }
  }

  if (cl.option_present('s')) {
    _check_scaffolds = 1;
    if (_verbose) {
      cerr << "Will try to match scaffolds\n";
    }
  }

  if (cl.option_present('y')) {
    _tabular_output = 0;
    if (_verbose) {
      cerr << "Will output a smiles file\n";
    }
  }

  if (cl.option_present('R')) {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      if (! ReadReactions(RxnFileType::kProto, fname)) {
        cerr << "Cannot read file of reactions " << fname << "\n";
        return 0;
      }
    }
  }

  if (cl.option_present('r')) {
    IWString fname;
    for (int i = 0; cl.value('r', fname, i); ++i) {
      if (! ReadReactions(RxnFileType::kMsi, fname)) {
        cerr << "Cannot read file of reactions " << fname << "\n";
        return 0;
      }
    }
  }

  if (_verbose) {
    cerr << "Read " << _reaction.size() << " reactions\n";
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 0;
  }

  if (cl.option_present('U')) {
    IWString fname = cl.option_value('U');
    if (! fname.ends_with(".smi")) {
      fname << ".smi";
    }

    if (! stream_for_unmatched_molecules.open(fname.null_terminated_chars())) {
      cerr << "Options::initialise:cannot open stream for unmatched molecules '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Unmatched molecules written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::ReadReactions(RxnFileType rxn_file_type, IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadReactions:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReactions(rxn_file_type, fname, input);
}

int
Options::ReadReactions(RxnFileType rxn_file_type, IWString& fname, iwstring_data_source& input) {
  input.set_skip_blank_lines(1);
  input.set_ignore_pattern("^#");

  IWString line;
  while (input.next_record(line)) {
    if (! ReadReaction(rxn_file_type, fname, line)) {
      cerr << "Options::ReadReaction:cannot read reaction '" << line << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadReaction(RxnFileType rxn_file_type, 
                      const IWString& fname,
                      const IWString& line) {
  std::optional<IWString> path_name = iwmisc::FileOrPath(fname, line);
  if (! path_name) {
    cerr << "Options::ReadReaction:cannot discern '" << line << "'\n";
    return 0;
  }
  switch (rxn_file_type) {
    case RxnFileType::kProto:
      return ReadReactionProto(*path_name);
    case RxnFileType::kMsi:
      return ReadReactionMsi(*path_name);
    default:
      cerr << "What kind of file is this?\n";
      return 0;
  }
}

int
Options::ReadReactionProto(IWString& fname) {
  std::optional<ReactionProto::Reaction> proto = 
    iwmisc::ReadTextProtoCommentsOK<ReactionProto::Reaction>(fname);
  if (! proto) {
    cerr << "Options::ReadReaction:cannot read proto file " << fname << '\n';
    return 0;
  }

  std::unique_ptr<ReactionAndStats> rxn = std::make_unique<ReactionAndStats>();
  if (! rxn->ConstructFromProto(*proto, fname)) {
    cerr << "Options::ReadReaction:cannot parse " << proto->ShortDebugString() << '\n';
    return 0;
  }

  _reaction << rxn.release();

  return 1;
}

int
Options::ReadReactionMsi(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr <<"Options::ReadReactionMsi:cannot open '" << fname << "'\n";
    return 0;
  }
  msi_object msi;
  if (! msi.read(input)) {
    cerr << "Options::ReadReactionMsi:cannot read msi data " << fname << "\n";
    return 0;
  }

  std::unique_ptr<ReactionAndStats> rxn = std::make_unique<ReactionAndStats>();
  if (! rxn->construct_from_msi_object(msi, _smc)) {
    cerr << "Options::ReadReactionMsi:cannot parse " << msi << '\n';
    return 0;
  }

  _reaction << rxn.release();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  for (int i = 0; i < _reaction.number_elements(); ++i) {
    output << " rxn " << i << ' ';
    _reaction[i]->Report(output);
  }
  return 1;
}

int
Options::MaybeDiscernInputType(const char * fname) {
  if (_input_type == FILE_TYPE_INVALID) {
    _input_type = discern_file_type_from_name(fname);
  }
  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    if (_element_transformations.process(m)) {
      IWString mname(m.name());
      mname << ".T";
      m.set_name(mname);
    }
  }

  return 1;
}

int
Options::MaybeWriteNonMatch(Molecule& m) {
  if (! stream_for_unmatched_molecules.is_open()) {
    return 1;
  }

  stream_for_unmatched_molecules << m.smiles() << _output_separator
      << m.name() << '\n';
  stream_for_unmatched_molecules.write_if_buffer_holds_more_than(4096);

  return 1;
}

class HaystackMolecule : public Molecule {
};

// Unique smiles stored in hashes come fro haystack molecules via various
// transformations. These are the transformations.
enum class VariantType {
  kUnchanged,
  kSubstructure,
  kNochiral,
  kGraph,
  kSkeleton,
  kReacted,
  kScaffold,
};

// We make various changes to the molecules in the haystack and store
// those variants in a hash - keyed by unique smiles.
// For the value we store a number for the haystack member that generated
// this unique smiles, and the kind of transformation.
struct HaystackMember {
  // Which molecule generated this entry.
  int molecule_number;
  // the kind of transformation that lead to this entry.
  VariantType type;
};

class Haystack {
  private:
    resizable_array_p<HaystackMolecule> _haystack;

    // All variants get loaded into this map.
    IW_STL_Hash_Map<IWString, HaystackMember> _usmi;

    // For each molecule, a substructure search target.
    std::unique_ptr<Molecule_to_Match[]> _target;

    // For each molecule, a Substructure_Query
    resizable_array_p<Substructure_Query> _query;

    int _molecules_examined = 0;
    int _found_as_exact_match = 0;
    int _found_as_non_chiral = 0;
    int _found_as_substructure = 0;
    int _found_as_skeleton = 0;
    int _found_as_graph = 0;
    int _found_as_reaction_variant = 0;
    int _found_as_scaffold = 0;
    int _molecules_with_no_matches = 0;

  // private functions

    int BuildInternalDataStructure(Options& options);
    int EnumerateReactions(Options& options);
    int EnumerateReactions(Options& options, HaystackMolecule& m);
    int EnumerateReactions(Options& options, int molecule_number);
    int EnumerateHits(Options& options,
                        ReactionAndStats& rxn,
                        int molecule_number,
                        int recursion_depth,
                        int max_recursion_depth,
                        const Set_of_Atoms* embedding);
    int StoreReactedForm(Options& options,
                       ReactionAndStats& rxn,
                       int molecule_number,
                       const Set_of_Atoms* embedding);
    int StoreScaffolds(Options& options);
    int StoreScaffold(Options& options, int molecule_number,
                        int * storage);

    int DoLookup(Options& options,
                   Molecule& initial_molecule,
                   Molecule& m,
                   IWString_and_File_Descriptor& output);
    int DoLookupInner(Options& options,
                   Molecule& initial_molecule,
                   Molecule& m,
                   IWString_and_File_Descriptor& output);

    int MatchAsNeedleVariant(Options& options,
                             Molecule& m,
                             IWString_and_File_Descriptor& output);
    int MatchAsNeedleVariant(Options& options,
                        Molecule& m,
                        ReactionAndStats& rxn,
                        const Set_of_Atoms* embedding,
                        IWString_and_File_Descriptor& output);
    int MatchAsSubstructure(Options& options,
                               Molecule& m,
                               IWString_and_File_Descriptor& output);
    int TryMatchSimultaneousChanges(Options& options,
                        Molecule& m,
                        ReactionAndStats& rxn,
                        Substructure_Results& sresults,
                        IWString_and_File_Descriptor& output);

    int GotMatch(Options& options,
                   Molecule& initial_molecule,
                   Molecule& m,
                   const HaystackMember& match,
                   IWString_and_File_Descriptor& output);

  public:
    Haystack();
    ~Haystack() {
    }

    int ReadHaystack(Options& options, const char * fname);
    int ReadHaystack(Options& options, 
                        data_source_and_type<HaystackMolecule>& input);

    int Process(Options& options,
                  Molecule& m,
                  IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Haystack::Haystack() {
  _molecules_examined = 0;
  _found_as_exact_match = 0;
  _found_as_non_chiral = 0;
  _found_as_substructure = 0;
  _found_as_skeleton = 0;
  _found_as_graph = 0;
  _found_as_reaction_variant = 0;
  _found_as_scaffold = 0;
  _molecules_with_no_matches = 0;
}

int
Haystack::ReadHaystack(Options& options,
                       const char * fname) {
  FileType input_type = discern_file_type_from_name(fname);
  if (input_type == FILE_TYPE_INVALID) {
    cerr << "Haystack::ReadHaystack:cannot discern type '" << fname << "'\n";
    return 0;
  }
  data_source_and_type<HaystackMolecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Haystack::ReadHaystack:cannot open " << fname << '\n';
    return 0;
  }

  return ReadHaystack(options, input);
}

int
Haystack::ReadHaystack(Options& options, 
                        data_source_and_type<HaystackMolecule>& input) {
  HaystackMolecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    if (! options.Preprocess(*m)) {
      delete m;
      return 0;
    }
    _haystack << m;
  }

  if (options.verbose()) {
    cerr << "Haystack::ReadHaystack:read " << _haystack.size() << " haystack molecules\n";
  }

  return BuildInternalDataStructure(options);
}

Molecule
TransformToSkeleton(const Molecule& m) {
  Molecule result(m);
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    result.set_atomic_number(i, 6);
  }
  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = result.bondi(i);
    result.set_bond_type_between_atoms(b->a1(), b->a2(), SINGLE_BOND);
  }
  result.remove_all_chiral_centres();
  return result;
}

Molecule
AsGraph(Mol2Graph& mol2graph, const Molecule& m) {
  Molecule result(m);
  result.remove_all_chiral_centres();
  result.change_to_graph_form(mol2graph);
  return result;
}

std::unique_ptr<Substructure_Query>
AsQuery(Molecule& m,
        Molecule_to_Query_Specifications& mqs) {
  std::unique_ptr<Substructure_Query> result = std::make_unique<Substructure_Query>();

  auto chiral_centres = m.ReleaseChiralCentres();
  int ok = result->create_from_molecule(m, mqs);
  if (! ok) {
    cerr << "AsQuery:warning, result may be invalid\n";
  }
  m.SetChiralCentres(std::move(chiral_centres));

  return result;
}

int
Haystack::BuildInternalDataStructure(Options& options) {
  const int haystack_size = _haystack.number_elements();
  if (haystack_size == 0) {
    return 0;
  }

  _target = std::make_unique<Molecule_to_Match[]>(haystack_size);

  for (int i = 0; i < haystack_size; ++i) {
    HaystackMolecule* m = _haystack[i];
    _target[i].initialise_molecule(m);

    HaystackMember value{i, VariantType::kUnchanged};
    _usmi.emplace(m->unique_smiles(), value);

    if (m->chiral_centres()) {
      Molecule mcopy(*m);
      mcopy.remove_all_chiral_centres();
      HaystackMember value{i, VariantType::kNochiral};
      _usmi.try_emplace(mcopy.unique_smiles(), value);
    }
  }

  _query.resize(haystack_size);
  for (HaystackMolecule * m : _haystack) {
    std::unique_ptr<Substructure_Query> qry = AsQuery(*m, options.mqs());
    _query << qry.release();
  }

  if (options.verbose()) {
    cerr << "Haystack::BuildInternalDataStructure:found " <<
            _usmi.size() << " unique smiles\n";
  }

  if (options.verbose() > 1) {
    for (const auto& [key, value] : _usmi) {
      cerr << key << '\n';
    }
  }

  for (int i = 0; i < haystack_size; ++i) {
    const HaystackMolecule * m = _haystack[i];
    Molecule skeleton = TransformToSkeleton(*m);
    HaystackMember value{i, VariantType::kSkeleton};
    _usmi.try_emplace(skeleton.unique_smiles(), value);
  }

  if (options.mol2graph().active()) {
    Mol2Graph& mol2graph = options.mol2graph();
    for (int i = 0; i < haystack_size; ++i) {
      HaystackMolecule * m = _haystack[i];
      Molecule as_graph = AsGraph(mol2graph, *m);
      HaystackMember value{i, VariantType::kGraph};
      _usmi.try_emplace(as_graph.unique_smiles(), value);
    }
  }

  if (options.check_scaffolds()) {
    StoreScaffolds(options);
  }

  if (options.reactions().size() > 0) {
    EnumerateReactions(options);
  }

  return haystack_size;
}

int
Haystack::Report(std::ostream& output) const {
  output << _molecules_examined << " molecules examined\n";
  output << _found_as_exact_match << " molecules found as unique smiles matches\n";
  output << _found_as_non_chiral << " molecules found as non chiral unique smiles matches\n";
  output << _found_as_substructure << " molecules found as substructures\n";
  output << _found_as_reaction_variant << " molecules found as reaction variants\n";
  output << _found_as_skeleton << " molecules found as skeleton matches\n";
  output << _found_as_graph << " molecules found as graph unique smiles\n";
  output << _found_as_scaffold << " molecules found as scaffold matches\n";
  output << _molecules_with_no_matches << " molecules did not match anything\n";
  return 1;
}

int
Haystack::StoreScaffolds(Options& options) {
  const int haystack_size = _haystack.number_elements();

  int max_natoms = 0;
  for (const HaystackMolecule* m : _haystack) {
    if (m->natoms() > max_natoms) {
      max_natoms = m->natoms();
    }
  }

  std::unique_ptr<int[]> storage = std::make_unique<int[]>(max_natoms);

  for (int i = 0; i < haystack_size; ++i) {
    StoreScaffold(options, i, storage.get());
  }

  return 1;
}

Molecule
ToScaffold(const Molecule& m,
           int * storage) {
  Molecule result;
  result.identify_spinach(storage);
  result.remove_atoms(storage);
  return result;
}

Molecule
ToScaffold(const Molecule& m) {
  std::unique_ptr<int[]> storage = std::make_unique<int[]>(m.natoms());
  return ToScaffold(m, storage.get());
}

int
Haystack::StoreScaffold(Options& options, int molecule_number,
                        int * storage) {
  Molecule * m = _haystack[molecule_number];
  Molecule scaffold = ToScaffold(*m, storage);
  HaystackMember value{molecule_number, VariantType::kScaffold};
  _usmi.try_emplace(scaffold.unique_smiles(), value);

  return 1;
}

int
Haystack::EnumerateReactions(Options& options) {
  const int haystack_size = _haystack.number_elements();

  for (int i = 0; i < haystack_size; ++i) {
    EnumerateReactions(options, i);
  }

  if (options.verbose()) {
    cerr << "Haystack::EnumerateReactions:after generationg reaction variants " << _usmi.size() << " haystack members\n";
    if (options.verbose() > 1) {
      for (const auto& [key, value] : _usmi) {
        cerr << key << '\n';
      }
    }
  }

  return 1;
}

int
Haystack::EnumerateReactions(Options& options, int molecule_number) {
  for (ReactionAndStats * rxn : options.reactions()) {
    Substructure_Results sresults;
    const int nhits = rxn->substructure_search(_target[molecule_number], sresults);
    if (nhits == 0) {
      continue;
    }

    for (const Set_of_Atoms * embedding : sresults.embeddings()) {
      StoreReactedForm(options, *rxn, molecule_number, embedding);
    }

    if (nhits > 1) {
      Molecule mcopy(*_haystack[molecule_number]);
      rxn->in_place_transformations(mcopy, sresults);
      HaystackMember value{molecule_number, VariantType::kReacted};
      _usmi.try_emplace(mcopy.unique_smiles(), value);
    }
  }

  return 1;
}

int
Haystack::StoreReactedForm(Options& options,
                       ReactionAndStats& rxn,
                       int molecule_number,
                       const Set_of_Atoms* embedding) {
  Molecule mcopy(*_haystack[molecule_number]);
  if (! rxn.in_place_transformation(mcopy, embedding)) {
    return 0;
  }

  HaystackMember value{molecule_number, VariantType::kReacted};
  auto [_, inserted] = _usmi.try_emplace(mcopy.unique_smiles(), value);

  return inserted;
}

const char *
MatchType(VariantType type) {
  switch (type) {
    case VariantType::kUnchanged:
      return "parent";
    case  VariantType::kNochiral:
      return "nochiral";
    case VariantType::kGraph:
      return "graph";
    case VariantType::kSkeleton:
      return "skeleton";
    case VariantType::kReacted:
      return "rxn";
    case VariantType::kSubstructure:
        return "substructure";
    default:
      cerr << "What kind of transformation\n";
      return "??";
  }
}

int
Haystack::GotMatch(Options& options,
                   Molecule& initial_molecule,
                   Molecule& m,
                   const HaystackMember& match,
                   IWString_and_File_Descriptor& output) {
  const char sep = options.output_separator();

  int molecule_number = match.molecule_number;

  if (options.tabular_output()) {
  output << initial_molecule.smiles() << sep << 
            initial_molecule.name() << sep <<
            MatchType(match.type) << sep <<
            m.smiles() << sep <<
            _haystack[molecule_number]->smiles() << sep <<
            _haystack[molecule_number]->name() << '\n';
  } else {
    output << initial_molecule.smiles() << sep <<
              initial_molecule.name() << sep << MatchType(match.type) << '\n';
    output << m.smiles() << sep << '\n';
    output << _haystack[molecule_number]->smiles() << sep <<
              _haystack[molecule_number]->name() << '\n';
  }

  output.write_if_buffer_holds_more_than(4096);

  switch (match.type) {
    case VariantType::kUnchanged:
      ++_found_as_exact_match;
      break;
    case VariantType::kNochiral:
      ++_found_as_non_chiral;
      break;
    case VariantType::kGraph:
      ++_found_as_graph;
      break;
    case VariantType::kSkeleton: 
      ++_found_as_skeleton;
      break;
    case VariantType::kReacted: 
      ++_found_as_reaction_variant;
      break;
    case VariantType::kSubstructure: 
      ++_found_as_substructure;
      break;
    case VariantType::kScaffold: 
      ++_found_as_scaffold;
      break;
    default:
      cerr << "What kind of transformation\n";
      return 0;
  }
  return 1;
}

// If `m` is found in _usmi, return true.
// Then, if `m` has chirality, remove it and try that.
int
Haystack::DoLookup(Options& options,
                   Molecule& initial_molecule,
                   Molecule& m,
                   IWString_and_File_Descriptor& output) {
  if (DoLookupInner(options, initial_molecule, m, output)) {
    return 1;
  }

  if (m.chiral_centres() == 0) {
    return 0;
  }

  resizable_array_p<Chiral_Centre> chiral_centres = m.ReleaseChiralCentres();
  int rc = DoLookupInner(options, initial_molecule, m, output);
  m.SetChiralCentres(std::move(chiral_centres));

  return rc;
}

int
Haystack::DoLookupInner(Options& options,
                   Molecule& initial_molecule,
                   Molecule& m,
                   IWString_and_File_Descriptor& output) {
  auto iter = _usmi.find(m.unique_smiles());

  if (iter == _usmi.cend()) {
    return 0;
  }

  return GotMatch(options, initial_molecule, m, iter->second, output);
}

int
Haystack::Process(Options& options,
                  Molecule& m,
                  IWString_and_File_Descriptor& output) {
  ++_molecules_examined;

  if (DoLookup(options, m, m, output)) {
    return 1;
  }

  if (m.chiral_centres()) {
    Molecule mcopy(m);
    mcopy.remove_all_chiral_centres();
    if (DoLookup(options, m, mcopy, output)) {
      return 1;
    }
  }

  if (options.perform_reactions_on_needles()) {
    if (MatchAsNeedleVariant(options, m, output)) {
      return 1;
    }
  }

  if (MatchAsSubstructure(options, m, output)) {
    return 1;
  }

  if (options.mol2graph().active()) {
    Molecule as_graph = AsGraph(options.mol2graph(), m);
    if (DoLookup(options, m, as_graph, output)) {
      return 1;
    }
  }

  Molecule skeleton = TransformToSkeleton(m);
  if (DoLookup(options, m, skeleton, output)) {
    return 1;
  }

  if (options.check_scaffolds()) {
    Molecule scaffold = ToScaffold(m);
    if (DoLookup(options, m, scaffold, output)) {
      return 1;
    }
  }

  ++_molecules_with_no_matches;

  options.MaybeWriteNonMatch(m);

  return 1;
}

int
Haystack::MatchAsNeedleVariant(Options& options,
                               Molecule& m,
                               IWString_and_File_Descriptor& output) {
  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  for (ReactionAndStats * rxn : options.reactions()) {
    const int nhits = rxn->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      if (MatchAsNeedleVariant(options, m, *rxn, embedding, output)) {
        return 1;
      }
    }

    if (nhits > 1) {
      if (TryMatchSimultaneousChanges(options, m, *rxn, sresults, output)) {
        return 1;
      }
    }
  }

  return 0;
}

int
Haystack::TryMatchSimultaneousChanges(Options& options,
                        Molecule& m,
                        ReactionAndStats& rxn,
                        Substructure_Results& sresults,
                        IWString_and_File_Descriptor& output) {
  Molecule mcopy(m);
  if (! rxn.in_place_transformations(mcopy, sresults)) {
    return 0;
  }

  if (DoLookup(options, m, mcopy, output)) {
    rxn.GotSuccessfulLookup();
    return 1;
  }

  if (mcopy.chiral_centres()) {
    mcopy.remove_all_chiral_centres();
    if (DoLookup(options, m, mcopy, output)) {
      rxn.GotSuccessfulLookup();
      return 1;
    }
  }

  return 0;
}

int
Haystack::MatchAsNeedleVariant(Options& options,
                        Molecule& m,
                        ReactionAndStats& rxn,
                        const Set_of_Atoms* embedding,
                        IWString_and_File_Descriptor& output) {
  Molecule mcopy(m);
  if (! rxn.in_place_transformation(mcopy, embedding)) {
    return 0;
  }

  if (DoLookup(options, m, mcopy, output)) {
    rxn.GotSuccessfulLookup();
    return 1;
  }

  mcopy.remove_all_chiral_centres();
  if (mcopy.chiral_centres()) {
    if (DoLookup(options, m, mcopy, output)) {
      rxn.GotSuccessfulLookup();
      return 1;
    }
  }

  // See if any of the haystack are substructures.
  Molecule_to_Match target(&mcopy);
  const int haystack_size = _haystack.number_elements();
  Substructure_Results sresults;
  for (int i = 0; i < haystack_size; ++i) {
    int nhits = _query[i]->substructure_search(target, sresults);
    if (nhits) {
      HaystackMember mtype{i, VariantType::kSubstructure};
      return GotMatch(options, m, mcopy, mtype, output);
    }
  }

  std::unique_ptr<Substructure_Query> qry = AsQuery(mcopy, options.mqs());
  for (int i = 0; i < haystack_size; ++i) {
    if (qry->substructure_search(_target[i], sresults)) {
      HaystackMember mtype{i, VariantType::kSubstructure};
      return GotMatch(options, m, mcopy, mtype, output);
    }
  }

  return 0;
}

int
Haystack::MatchAsSubstructure(Options& options,
                              Molecule& m,
                              IWString_and_File_Descriptor& output) {
  std::unique_ptr<Substructure_Query> qry = AsQuery(m, options.mqs());

  Substructure_Results sresults;
  const int haystack_size = _haystack.number_elements();
  for (int i = 0; i < haystack_size; ++i) {
    const int nhits = qry->substructure_search(_target[i], sresults);
    if (nhits == 0) {
      continue;
    }
    HaystackMember mtype{i, VariantType::kSubstructure};
    return GotMatch(options, m, m, mtype, output);
  }

  // Now try the other way.
  Molecule_to_Match target(&m);
  for (int i = 0; i < haystack_size; ++i) {
    const int nhits = _query[i]->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    HaystackMember mtype{i, VariantType::kSubstructure};
    return GotMatch(options, m, m, mtype, output);
  }

  return 0;
}

int
SubstructureMcs(Options& options,
                Haystack& haystack,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return haystack.Process(options, m, output);
}

int
SubstructureMcs(Options& options,
                Haystack& haystack,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    options.ReadAnotherMolecule();

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! SubstructureMcs(options, haystack, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
SubstructureMcs(Options& options,
                Haystack& haystack,
                const char * fname,
                IWString_and_File_Descriptor& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.good()) {
    cerr << "SubstructureMcs:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return SubstructureMcs(options, haystack, input, output);
}

int
SubstructureMcs(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:R:r:G:aU:sy");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.number_elements() == 1) {
    cerr << "SubstructureMcs:requires at least two files, the first is the haystack\n";
    Usage(1);
  }

  Haystack haystack;
  if (! haystack.ReadHaystack(options, cl[0])) {
    cerr << "Cannot read haystack '" << cl[0] << "'\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 1; i < cl.number_elements(); ++i) {
    if (! SubstructureMcs(options, haystack, cl[i], output)) {
      cerr << "SubstructureMcs::fatal error processing '" << cl[i] << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
    haystack.Report(cerr);
  }

  return 0;
}

}  // namespace substructure_mcs

int
main(int argc, char ** argv) {

  int rc = substructure_mcs::SubstructureMcs(argc, argv);

  return rc;
}
