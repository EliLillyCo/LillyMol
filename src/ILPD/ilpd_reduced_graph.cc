// Convert full connection tables of iomisable lipid molecules to
// reduced graph forms.

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "ILPD/ilpd_reduced_graph.pb.h"

namespace ilpd_reduced_graph {

using std::cerr;

// By convention the Usage function tells how to use the tool.
void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
cerr << R"(Transforms Ionisable Lipid molecules into reduced grah forms
-C <fname>              read paramaters from <fname> - always needed
-S <fname>              write descriptors to <fname>
-F <fname>              write fragmented lipids to <fname>
-v                      verbose output\n";
)";
// clang-format on

  ::exit(rc);
}

// when scanning a chain, we can keep track of various things found in the region.
struct ChainData {
  int unsaturation = 0;
  int ring_bond_count = 0;
};

// Values assigned to the _classification array in PerMoleculeData objects.
static constexpr int kUnclassified = 0;
static constexpr int kHead = 1;
static constexpr int kDivider = 2;
static constexpr int kBranch = 3;
static constexpr int kCarbonJoin = 4;
static constexpr int kBetweenLinker = 5;
static constexpr int kCarbonChain = 6;
static constexpr int kLinker = 7;
static constexpr int kOtherHead = 8;
// chain regions get numbers [5...]

static constexpr int kNatoms = 0;
static constexpr int kMaxDist = 1;
static constexpr int kAveDist = 2;
static constexpr int kSymmetricAtoms = 3;
static constexpr int kHeadCarbon = 4;
static constexpr int kAtomsInHead = 5;
static constexpr int kAtomsInTrimmedHead = 6;

static constexpr int kCarbonChainCount = 7;
static constexpr int kCarbonChainBranchPoints = 8;
static constexpr int kNumberDividers = 9;
static constexpr int kNumberBranches = 10;

static constexpr int kCarbonChains = 11;
static constexpr int kAtomsInCarbonChains = 12;

static constexpr int kMaxSymmSize = 13;
static constexpr int kAveSymmSize = 14;
static constexpr int kMaxSymmDist = 15;
static constexpr int kSymmClass = 16;

static constexpr int kNchains = 17;
static constexpr int kMaxChainLength = 18;
static constexpr int kInstancesMaxChainLength = 19;
static constexpr int kMinChainLength = 20;
static constexpr int kAveChainLength = 21;

static constexpr int kNumberFragments = 22;
static constexpr int kNumberUniqueFragments = 23;

// This must always be the last descriptor.
static constexpr int kNextToAssign = 24;

// A container for holding descriptors.
class IlpdDescriptors {
  private:
    int _n;
    float* _d;

  public:
    IlpdDescriptors();
    ~IlpdDescriptors();

    void set(int ndx, float v) {
      _d[ndx] = v;
    }
    void set(int ndx, int v) {
      _d[ndx] = static_cast<float>(v);
    }
    void set(int ndx, uint32_t v) {
      _d[ndx] = static_cast<float>(v);
    }
    void set(int ndx, uint64_t v) {
      _d[ndx] = static_cast<float>(v);
    }
    void set(int ndx, double v) {
      _d[ndx] = static_cast<float>(v);
    }

    const float* values() const {
      return _d;
    }

    int WriteDescriptorHeader(IWString_and_File_Descriptor& output, const IWString& sep) const;
};

IlpdDescriptors::IlpdDescriptors() {
  _n = kNextToAssign;
  _d = new float[kNextToAssign];
  std::fill_n(_d, kNextToAssign, 0);
}

IlpdDescriptors::~IlpdDescriptors() {
  _n = 0;
  delete [] _d;
}

int
IlpdDescriptors::WriteDescriptorHeader(IWString_and_File_Descriptor& output, const IWString& sep) const {
  IWString dname[kNextToAssign];
  for (int i = 0; i < kNextToAssign; ++i) {
    dname[i] = ".";
  }

  output << "Name";
  dname[kNatoms] = "natoms";
  dname[kMaxDist] = "maxdist";
  dname[kAveDist] = "ave_dist";
  dname[kSymmetricAtoms] = "symmetric_atoms";
  dname[kHeadCarbon] = "head_carbon";
  dname[kAtomsInHead] = "atoms_in_head";
  dname[kAtomsInTrimmedHead] = "atoms_in_trimmed_head";
  dname[kCarbonChainCount] = "carbon_chain_count";
  dname[kCarbonChainBranchPoints] = "carbon_chain_branch_points";
  dname[kNumberDividers] = "number_dividers";
  dname[kNumberBranches] = "number_branches";
  dname[kCarbonChains] = "carbon_chains";
  dname[kAtomsInCarbonChains] = "atoms_in_carbon_chains";
  dname[kMaxSymmSize] = "atoms_in_largest_symmetry_group";
  dname[kMaxSymmDist]  = "max_distance_btw_symmetric_atoms";
  dname[kSymmClass]  = "number_symmetry_classes";
  dname[kNchains]  = "number_of_chains";
  dname[kMaxChainLength]  = "length_of_longest_chain";
  dname[kInstancesMaxChainLength]  = "instances_length_of_longest_chain";
  dname[kMinChainLength]  = "shortest_chain_length";
  dname[kAveChainLength]  = "ave_chain_length";
  dname[kNumberFragments]  = "number_fragments";
  dname[kNumberUniqueFragments]  = "number_unique_fragments";

  for (int i = 0; i < kNextToAssign; ++i) {
    output << sep << dname[i];
  }

  output << '\n';

  return 1;
}

class PerMoleculeData {
  private:
    Molecule& _m;

    int* _attached_heteroatom_count;


    extending_resizable_array<int> _carbon_chain_length;

    int _carbon_chain_branch_points;
    int _carbon_chain_count;

    // Used as a temporary during region finding.
    int* _visited;

    // Will be set for atoms that are the join points between chains.
    int* _chain_divider;

    int* _classification;

    IlpdDescriptors _descriptors;

  public:
    PerMoleculeData(Molecule& m);
    ~PerMoleculeData();

#ifdef ADASDASD
    int* chain_dividers() {
      return _chain_divider;
    }
    int chain_divider(int i) const {
      return _chain_divider[i];
    }
    int& chain_divider(int i) {
      return _chain_divider[i];
    }
#endif

    template <typename T>
    int WriteIsotopicallyLabelled(bool add_newline, T& output) const;

    int* visited() {
      return _visited;
    }
    int visited(atom_number_t a) const {
      return _visited[a];
    }
    int& visited(atom_number_t a) {
      return _visited[a];
    }

    void ResetVisited() {
      std::fill_n(_visited, _m.natoms(), 0);
    }

    int* classification() {
      return _classification;
    }
    int& classification(atom_number_t a) {
      return _classification[a];
    }
    int classification(atom_number_t a) const {
      return _classification[a];
    }

    int attached_heteroatom_count(atom_number_t a) const {
      return _attached_heteroatom_count[a];
    }

    void increment_carbon_chain_length(int len) {
      ++_carbon_chain_length[len];
    }
    const extending_resizable_array<int>& carbon_chain_length() const {
      return _carbon_chain_length;
    }

    void increment_carbon_chain_branch_points() {
      ++_carbon_chain_branch_points;
    }
    int carbon_chain_branch_points() const {
      return _carbon_chain_branch_points;
    }

    void increment_carbon_chain_count() {
      ++_carbon_chain_count;
    }
    int carbon_chain_count() const {
      return _carbon_chain_count;
    }

    const IlpdDescriptors& descriptors() const {
      return _descriptors;
    }
    IlpdDescriptors& descriptors() {
      return _descriptors;
    }
};

PerMoleculeData::PerMoleculeData(Molecule& m) : _m(m) {
  const int matoms = _m.natoms();

  _chain_divider = new_int(matoms);

  _visited = new int[matoms];

  _classification = new_int(matoms);

  _attached_heteroatom_count = new_int(matoms);

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (m.atomic_number(a1) != 6) {
      ++_attached_heteroatom_count[a2];
    }
    if (m.atomic_number(a2) != 6) {
      ++_attached_heteroatom_count[a1];
    }
  }

  _carbon_chain_branch_points = 0;

  _carbon_chain_count = 0;
}

PerMoleculeData::~PerMoleculeData() {
  delete [] _chain_divider;
  delete [] _visited;
  delete [] _classification;
  delete [] _attached_heteroatom_count;
}

template <typename T>
int
PerMoleculeData::WriteIsotopicallyLabelled(bool add_newline, T& output) const {
  Molecule mcopy(_m);
  mcopy.set_isotopes(_classification);
  output << mcopy.smiles() << ' ' << mcopy.name();
  if (add_newline) {
    output << '\n';
  }
  return 1;
}

// A chain region is characterised by the atoms that comprise it.
// Some
class Chain : public Set_of_Atoms {
  private:
    // A unique id for this chain
    int _uid;

    uint32_t _unsaturation;
    uint32_t _ring_atoms;

    std::vector<int> _join_point_indices;

  public:
    Chain();

    int Build(Molecule& m, atom_number_t zatom, PerMoleculeData& pmd);
};

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    Charge_Assigner _charge_assigner;

    // This should be an option, but for now we just ignore them.
    int _ignore_no_charged_atoms = 1;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    // There are molecules that are not hit by the charge assigner. We can set up
    // queries, that are applied to molecules without formal charges.
    resizable_array_p<Substructure_Query> _if_not_charged;

    // Substructure queries that define groups known to be the functional
    // groups where sidechains are attached.
    resizable_array_p<Substructure_Query> _chain_dividers;

    // Queries that define short sections between chain_dividers queries.
    // These are usually CC or CCC groupings. The query should include the two
    // chain divider atoms, O=C-O-C-C-O-C=O for example.
    resizable_array_p<Substructure_Query> _between_linker_query;

    // Queries for identifying a carbon that at the end of long chains.
    resizable_array_p<Substructure_Query> _carbon_chain_join;

    // Keep track of how many times molecules have 0,1,2,... hits to 
    // chain divider queries.
    extending_resizable_array<int> _number_divider_query_hits;

    // Keep track of how many times each query hits a molecule.
    int* _divider_query_hits;

    resizable_array_p<Substructure_Query> _branch_query;

    // See proto;
    resizable_array_p<Substructure_Query> _break_tail_at_attachment;

    // See proto
    resizable_array_p<Substructure_Query> _break_head_at;

    // See proto
    int _label_chains_from_head;

    // See proto.
    int _head_linker_tail;


    int _min_chain_length;

    int _molecules_read = 0;
    int _uncharged_molecules = 0;

    IWString_and_File_Descriptor _stream_for_descriptors;
    IWString_and_File_Descriptor _stream_for_fragments;

    std::unique_ptr<ipld::Options> _proto;

    IWString _output_separator;

  // Private functions
    int Process(Molecule& m,
                 PerMoleculeData& pmd,
                 IWString_and_File_Descriptor& output);
    int ToReducedGraph(Molecule& m, PerMoleculeData& pmd);
    int ProcessChargeState(Molecule& m);
    int SatisfiesBranchSizeConstraint(Molecule& m, atom_number_t zatom, PerMoleculeData& pmd);
    int IdentifyCarbonChainBranchPoints(Molecule& m, Molecule_to_Match& target, PerMoleculeData& pmd);
    int IdentifyAllCarbonChainRegion(Molecule& m, atom_number_t zatom,
                        PerMoleculeData& pmd, Set_of_Atoms& three_connected);
    int Identify2Chain(Molecule& m, const Set_of_Atoms& embedding, PerMoleculeData& pmd);
    int MaybeSplitTails(resizable_array_p<Molecule>& fragments);
    int MaybeWriteFragments(Molecule& parent, PerMoleculeData& pmd);
    int WriteDescriptors(Molecule& m, const PerMoleculeData& pmd);
    int MaybeTruncateHead(Molecule& m, PerMoleculeData& pmd);
    int TruncateHead(Molecule& m, atom_number_t zatom, int* dtb,
             PerMoleculeData& pmd);
    int HeadLinkerTail(Molecule& m, PerMoleculeData& pmd,
                        IWString_and_File_Descriptor& output);

  public:
    Options();
    ~Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _uncharged_molecules = 0;
  _divider_query_hits = nullptr;
  _min_chain_length = 0;
  _head_linker_tail = 0;
  _label_chains_from_head = 0;
  _output_separator = ' ';
}

Options::~Options() {
  delete [] _divider_query_hits;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
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

  if (cl.option_present('N')) {
    if (! _charge_assigner.construct_from_command_line(cl, _verbose, 'N')) {
      cerr << "Cannot initialise charge assigner (-N)\n";
      return 0;
    }
  }

  if (! cl.option_present('C')) {
    cerr << "Must specify ipld::Options textproto file via the -C option\n";
    Usage(1);
  }

  if (cl.option_present('C')) {
    IWString fname = cl.option_value('C');
    _proto = std::move(iwmisc::ReadTextProtoPtr<ipld::Options>(fname));
    if (! _proto) {
      cerr << "Options::Initialise:cannot read textproto -C '" << fname << "'\n";
      return 0;
    }
  }

  //cerr << "Proto is\n";
  //cerr << _proto->ShortDebugString() << '\n';
  for (const std::string& s : _proto->applied_to_neutral()) {
    IWString tmp(s);
    if (! process_cmdline_token('*', tmp, _if_not_charged, _verbose)) {
      cerr << "Invalid if not charged query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _if_not_charged.size() << 
              " charge queries for uncharged molecules\n";
    }
  }

  cerr << "Have " << _proto->chain_dividers_size() << " chain divider queries\n";
  for (const std::string& s : _proto->chain_dividers()) {
    IWString tmp(s);
    tmp.ExpandEnvironmentVariablesInPlace();
    if (! process_cmdline_token('*', tmp, _chain_dividers, _verbose)) {
      cerr << "Invalid chain dividers query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _chain_dividers.size() << " chain divider queries\n";
    }

    _divider_query_hits = new_int(_chain_dividers.size());
  }

  for (const std::string& s : _proto->branch_query()) {
    IWString tmp(s);
    tmp.ExpandEnvironmentVariablesInPlace();
    if (! process_cmdline_token('*', tmp, _branch_query, _verbose)) {
      cerr << "Invalid chain branch query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _branch_query.size() << " chain branch queries\n";
    }
  }

  if (_chain_dividers.empty()) {
    cerr << "Options::Initialise:no chain divider queries, cannot continue\n";
    return 0;
  }

  for (const std::string& s : _proto->between_linker_query()) {
    IWString tmp(s);
    tmp.ExpandEnvironmentVariablesInPlace();
    if (! process_cmdline_token('*', tmp, _between_linker_query, _verbose)) {
      cerr << "Invalid between linker query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _between_linker_query.size() << " between linker queries\n";
    }
  }

  for (const std::string& s : _proto->carbon_chain_join()) {
    IWString tmp(s);
    tmp.ExpandEnvironmentVariablesInPlace();
    if (! process_cmdline_token('*', tmp, _carbon_chain_join, _verbose)) {
      cerr << "Invalid carbon chain join query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _carbon_chain_join.size() << " carbon chain join queries\n";
    }

    for (Substructure_Query* q : _carbon_chain_join) {
      q->set_find_one_embedding_per_atom(1);
    }
  }

  for (const std::string& s : _proto->break_tail_at_attachment()) {
    IWString tmp(s);
    tmp.ExpandEnvironmentVariablesInPlace();
    if (! process_cmdline_token('*', tmp, _break_tail_at_attachment, _verbose)) {
      cerr << "Invalid tail break at head query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _break_tail_at_attachment.size() << " break carbon at head join queries\n";
    }

    for (Substructure_Query* q : _break_tail_at_attachment) {
      q->set_find_one_embedding_per_atom(1);
    }
  }

  for (const std::string& s : _proto->break_head_at()) {
    IWString tmp(s);
    tmp.ExpandEnvironmentVariablesInPlace();
    if (! process_cmdline_token('*', tmp, _break_head_at, _verbose)) {
      cerr << "Invalid break head at query specification '" << s << "'\n";
      cerr << _proto->ShortDebugString() << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _break_head_at.size() << " break head at queries\n";
    }

    for (Substructure_Query* q : _break_head_at) {
      q->set_find_one_embedding_per_atom(1);
    }
  }

  if (_proto->has_min_chain_length()) {
    _min_chain_length = _proto->min_chain_length();
  }

  if (_proto->has_head_linker_tail()) {
    _head_linker_tail = _proto->head_linker_tail();
  }
  if (_proto->has_label_chains_from_head()) {
    _label_chains_from_head = _proto->label_chains_from_head();
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify name of descriptor file output via the -S option\n";
    return 0;
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! _stream_for_descriptors.open(fname)) {
      cerr << "Cannot open stream for molecular descriptors '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Descriptors written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    fname.EnsureEndsWith(".smi");
    if (! _stream_for_fragments.open(fname)) {
      cerr << "Cannot open stream for fragments '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Fragments written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  // Other information about what has happened.

  for (uint32_t i = 0; i < _chain_dividers.size(); ++i) {
    output << "chain dividier " << i << " query hit " << 
              _divider_query_hits[i] << " molecules\n";
  }

  for (int i = 0; i < _number_divider_query_hits.number_elements(); ++i) {
    if (_number_divider_query_hits[i]) {
      output << _number_divider_query_hits[i] << " molecules had " << i <<
                " substructure matches to chain divider queries\n";
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

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
    _element_transformations.process(m);
  }

  return 1;
}

uint32_t
IdentifyChain(Molecule& m, atom_number_t zatom, ChainData& chain_data, PerMoleculeData& pmd) {

  pmd.visited(zatom) = 1;
  if (m.ring_bond_count(zatom)) {
    ++chain_data.ring_bond_count;
  }

  uint32_t rc = 1;
  for (const Bond* b : m[zatom]) {
    const atom_number_t o = b->other(zatom);
    if (pmd.visited(o)) {
      continue;
    }

    const Atom& ao = m[o];
    if (ao.atomic_number() != 6) {
      continue;
    }

    if (! b->is_single_bond()) {
      ++chain_data.unsaturation;
    }

    if (pmd.classification(o) == kDivider) {
      continue;
    }

    rc += IdentifyChain(m, o, chain_data, pmd);
  }

  return rc;
}

// This version stops at a 3 connected atom. But we accommodate CH3 substituents on 
// the ring.
uint32_t
Identify2ChainInner(Molecule& m, atom_number_t zatom, ChainData& chain_data, PerMoleculeData& pmd) {

  pmd.visited(zatom) = 1;
  if (m.ring_bond_count(zatom)) {
    ++chain_data.ring_bond_count;
  }

  uint32_t rc = 1;
  atom_number_t next_atom = kInvalidAtomNumber;
  for (const Bond* b : m[zatom]) {
    const atom_number_t o = b->other(zatom);
    if (pmd.visited(o)) {
      continue;
    }

    const Atom& ao = m[o];
    if (ao.atomic_number() != 6) {
      continue;
    }

    if (b->is_single_bond() && ao.ncon() == 1) {
      pmd.visited(o) = 1;
      ++rc;
      continue;
    }

    if (! b->is_single_bond()) {
      ++chain_data.unsaturation;
    }

    if (pmd.classification(o) == kDivider) {
      continue;
    }

    if (next_atom == kInvalidAtomNumber) {
      next_atom = o;
    } else {
      return 0;
    }
  }

  if (next_atom != kInvalidAtomNumber) {
    rc += Identify2ChainInner(m, next_atom, chain_data, pmd);
  }

  return rc;
}

uint32_t
IdentifyHead(Molecule& m, atom_number_t zatom, PerMoleculeData& pmd) {

  pmd.classification(zatom) = kHead;
  // cerr << "Head continues to " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << '\n';

  uint32_t rc = 1;
  for (const Bond* b : m[zatom]) {
    const atom_number_t o = b->other(zatom);
    if (pmd.classification(o)) {
      continue;
    }

    rc += IdentifyHead(m, o, pmd);
  }

  return rc;
}

int
Options::ToReducedGraph(Molecule& m, PerMoleculeData& pmd) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> visited = std::make_unique<int[]>(matoms);
  // As we identify regions of the molecule we assign a unique id.
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.ncon() != 1 || a.atomic_number() != 6) {
      continue;
    }

    ChainData chain_data;
    std::fill_n(pmd.visited(), matoms, 0);

    uint32_t chain_length = IdentifyChain(m, i, chain_data, pmd);
    if (chain_length < _proto->min_chain_length()) {
      continue;
    }
  }

  return 1;
}

int
Options::ProcessChargeState(Molecule& m) {

  if (_charge_assigner.active()) {
    _charge_assigner.process(m);
  }

  // This is not entirely robust, hopefully nothing in the non-head region gets charged.
  if (m.has_formal_charges()) {
    return 1;
  }

  if (_if_not_charged.empty()) {
    cerr << "No formal charges in " << m.name() << '\n';
    return 0;
  }

  Molecule_to_Match target(&m);

  for (Substructure_Query* q : _if_not_charged) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }
    m.set_formal_charge(sresults.embedding(0)->item(0), +1);
    return 1;
  }

  cerr << "None of " << _if_not_charged.size() << " charge queries matched " <<
          m.name() << '\n';
  cerr << m.smiles() << '\n';
  return 0;
}

int
DiscernSymmetryRelated(Molecule& m, PerMoleculeData& pmd) {
  int nsymm = m.number_symmetry_classes();

  pmd.descriptors().set(kSymmClass, nsymm);

  const int matoms = m.natoms();
  const int* symmetry_class = m.symmetry_classes();

  m.recompute_distance_matrix();

  std::unique_ptr<int[]> sdone = std::make_unique<int[]>(nsymm + 1);
  std::fill_n(sdone.get(), nsymm + 1, 0);

  int atoms_in_largest_symmetry_class = 0;
  Accumulator_Int<uint32_t> acc_size;
  Accumulator_Int<uint32_t> acc_longest_distance;
  Accumulator_Int<uint32_t> acc_multiplicity;  // not really implemented yet, need to think more about this.
  int symmetric_atoms = 0;

  for (int i = 0; i < matoms; ++i) {
    const int s = symmetry_class[i];
    if (sdone[s]) {
      continue;
    }

    sdone[s] = 1;

    int atoms_in_class = 0;
    int max_distance = 0;
    int multiplicity = 0;

    for (int j = i + 1; j < matoms; ++j) {
      if (s != symmetry_class[j]) {
        continue;
      }
      ++atoms_in_class;
      const int d = m.bonds_between(i, j);
      if (d > max_distance) {
        max_distance = d;
        multiplicity = 1;
      } else if (d == max_distance) {
        ++multiplicity;
      }
    }

    if (atoms_in_class > atoms_in_largest_symmetry_class) {
      atoms_in_largest_symmetry_class = atoms_in_class;
    }
    acc_size.extra(atoms_in_class);
    acc_longest_distance.extra(max_distance);
    acc_multiplicity.extra(multiplicity);
    if (atoms_in_class > 1) {
      ++symmetric_atoms;
    }
  }

  pmd.descriptors().set(kSymmetricAtoms, symmetric_atoms);
  pmd.descriptors().set(kMaxSymmSize, acc_size.maxval());
  pmd.descriptors().set(kAveSymmSize, acc_size.average());
  pmd.descriptors().set(kMaxSymmDist, acc_longest_distance.maxval());
  pmd.descriptors().set(kMaxDist, m.longest_path());

  Accumulator_Int<uint32_t> acc_all_distances;

  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      acc_all_distances.extra(m.bonds_between(i, j));
    }
  }

  pmd.descriptors().set(kAveDist, acc_all_distances.average());

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  if (! ProcessChargeState(m)) {
    ++_uncharged_molecules;
    return _ignore_no_charged_atoms;
  }

  PerMoleculeData pmd(m);
  pmd.descriptors().set(kNatoms, m.natoms());

  static bool first_call = true;

  if (first_call && _stream_for_descriptors.is_open()) {
    pmd.descriptors().WriteDescriptorHeader(_stream_for_descriptors, _output_separator);
    first_call = false;
  }

  DiscernSymmetryRelated(m, pmd);

  Molecule_to_Match target(&m);

  cerr << "Have " << _chain_dividers.size() << " chain divider queries\n";
  if (_chain_dividers.size() > 0) {
    int nhits = 0;
    int queries_matching = 0;
    for (uint32_t i = 0; i < _chain_dividers.size(); ++i) {
      Substructure_Results sresults;
      if (_chain_dividers[i]->substructure_search(target, sresults) == 0) {
        continue;
      }

      ++_divider_query_hits[i];

      sresults.each_embedding_set_vector(pmd.classification(), kDivider);
      nhits += sresults.number_embeddings();
      ++queries_matching;
    }

    if (nhits == 0) {
      cerr << "Options::Process:no hits to " << _chain_dividers.size() <<
              " chain divider queries\n";
    }

    ++_number_divider_query_hits[nhits];
    if (_verbose) {
      cerr << m.name() << " had " << nhits << " to chain divider queries\n";
    }

    pmd.descriptors().set(kNumberDividers, nhits);

    pmd.WriteIsotopicallyLabelled(false, cerr);
    cerr << " %chain_divider " << nhits << '\n';
  }

  if (_branch_query.size() > 0) {
    int nhits = 0;
    int queries_matching = 0;
    for (uint32_t i = 0; i < _branch_query.size(); ++i) {
      Substructure_Results sresults;
      if (_branch_query[i]->substructure_search(target, sresults) == 0) {
        continue;
      }
      cerr << "Got branch query match\n";
      ++queries_matching;

      for (const Set_of_Atoms* e : sresults.embeddings()) {
        for (atom_number_t a : *e) {
          if (pmd.classification(a) != 0) {
            continue;
          }
          if (! SatisfiesBranchSizeConstraint(m, a, pmd)) {
            continue;
          }

          pmd.classification(a) = kBranch;
          ++nhits;
        }
      }
    }

    pmd.descriptors().set(kNumberBranches, nhits);

    pmd.WriteIsotopicallyLabelled(false, cerr);
    cerr << " %chain_splitter " << nhits << '\n';
  }

  if (_between_linker_query.size() > 0) {
    int nhits = 0;
    for (uint32_t i = 0; i < _between_linker_query.size(); ++i) {
      Substructure_Results sresults;
      if (_between_linker_query[i]->substructure_search(target, sresults) == 0) {
        continue;
      }
      for (const Set_of_Atoms* e : sresults.embeddings()) {
        for (atom_number_t a : *e) {
          if (pmd.classification(a)) {
            continue;
          }
          pmd.classification(a) = kBetweenLinker;
          ++nhits;
        }
      }
    }
    pmd.WriteIsotopicallyLabelled(false, cerr);
    cerr << " %between_linker " << nhits << '\n';
  }

  // Atoms where Carbon chains join are too hard to identify via
  // queries.
  IdentifyCarbonChainBranchPoints(m, target, pmd);

  pmd.descriptors().set(kCarbonChainCount, pmd.carbon_chain_count());
  pmd.descriptors().set(kCarbonChainBranchPoints, pmd.carbon_chain_branch_points());
  Accumulator_Int<uint32_t> acc_chain_length;
  const extending_resizable_array<int>& carbon_chain_length = pmd.carbon_chain_length();
  int multiplicity_longest_chain = 1;
  for (int i = 0; i < carbon_chain_length.number_elements(); ++i) {
    if (carbon_chain_length[i] == 0) {
      continue;
    }
    acc_chain_length.extra(i, carbon_chain_length[i]);
    multiplicity_longest_chain = carbon_chain_length[i];
  }
  pmd.descriptors().set(kNchains, acc_chain_length.n());
  pmd.descriptors().set(kMaxChainLength, acc_chain_length.maxval());
  pmd.descriptors().set(kInstancesMaxChainLength, multiplicity_longest_chain);
  pmd.descriptors().set(kMinChainLength, acc_chain_length.minval());
  pmd.descriptors().set(kAveChainLength, acc_chain_length.average());

  return Process(m, pmd, output);
}

std::optional<int>
ChainSize(const Molecule& m, atom_number_t zatom, PerMoleculeData& pmd, int flag) {
  pmd.visited(zatom) = flag;

  int rc = 1;
  for (const Bond * b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (pmd.visited(o) || pmd.classification(o)) {
      continue;
    }
    if (m.atomic_number(o) != 6) {
      return std::nullopt;
    }

    std::optional<int> q = ChainSize(m, o, pmd, flag);
    if (! q) {
      return std::nullopt;
    }
    rc += *q;
  }

  return rc;
}

// Atom `zatom` matches the query for a branch in a chain. Is this
// consistent with requirements for min_chain_length
int
Options::SatisfiesBranchSizeConstraint(Molecule& m, atom_number_t zatom,
                PerMoleculeData& pmd) {
  if (_min_chain_length == 0) {
    return 1;
  }

  pmd.ResetVisited();
  int flag = 1;
  pmd.visited(zatom) = flag;

  int number_ok_chains = 0;
  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      continue;
    }
    const atom_number_t o = b->other(zatom);
    if (pmd.classification(o) != 0) {
      continue;
    }
    if (m.atomic_number(o) != 6 || m.ncon(o) == 1) {
      continue;
    }

    std::optional<int> csize = ChainSize(m, o, pmd, flag);
    if (! csize) {
      continue;
    }
    if (*csize >= _min_chain_length) {
      ++number_ok_chains;
      ++flag;
    }
  }
 
  if (number_ok_chains < 2) {
    return 0;
  }

  pmd.classification(zatom) = kBranch;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
  }

  return 1;
}

// By convention matched atom 0 must be a 3 connected carbon, and matched
// atoms 1 and 2 must be [CH2], each going down a different chain.
int
Options::Identify2Chain(Molecule& m, const Set_of_Atoms& embedding,
               PerMoleculeData& pmd) {
  cerr << "Identify2Chain atoms " << embedding << '\n';
  atom_number_t c0 = embedding[0];
  if (m.atomic_number(c0) != 6 || m.ncon(c0) != 3) {
    cerr << "Identify2Chain:Not [CD3] " << m.smarts_equivalent_for_atom(c0) << '\n';
    return 0;
  }
  const atom_number_t c1 = embedding[1];
  if (m.atomic_number(c1) != 6 || m.ncon(c1) != 2) {
    cerr << "Identify2Chain:Not [CD2] " << m.smarts_equivalent_for_atom(c1) << '\n';
    return 0;
  }
  const atom_number_t c2 = embedding[2];
  if (m.atomic_number(c2) != 6 || m.ncon(c2) != 2) {
    cerr << "Identify2Chain:Not [CD2] " << m.smarts_equivalent_for_atom(c2) << '\n';
    return 0;
  }

  const int matoms = m.natoms();

  pmd.ResetVisited();
  pmd.visited(c0) = 1;
  ChainData chain_data1;
  int chain1 = Identify2ChainInner(m, c1, chain_data1, pmd);
  if (chain1 < _min_chain_length) {
    return 0;
  }

  std::unique_ptr<int[]> copy_visited = std::make_unique<int[]>(matoms);
  pmd.visited(c0) = 0;
  std::copy_n(pmd.visited(), matoms, copy_visited.get());

  pmd.ResetVisited();
  pmd.visited(c0) = 1;
  ChainData chain_data2;
  int chain2 = Identify2ChainInner(m, c2, chain_data2, pmd);
  if (chain2 < _min_chain_length) {
    return 0;
  }
  pmd.visited(c0) = 0;

  for (int i = 0; i < matoms; ++i) {
    if (pmd.visited(i) || copy_visited[i]) {
      pmd.classification(i) = kCarbonChain;
    }
  }

  pmd.classification(c0) = kCarbonJoin;
  cerr << "3 chain " << c0 << " atom " << c1 << ' ' << chain1 << " atom " << c2 << " " << chain2 << '\n';

//pmd.increment_carbon_chain_length(chain1);
//pmd.increment_carbon_chain_length(chain2);
  pmd.increment_carbon_chain_branch_points();
  pmd.increment_carbon_chain_count();

  return 1;
}

int
Options::IdentifyCarbonChainBranchPoints(Molecule& m, 
                Molecule_to_Match& target, PerMoleculeData& pmd) {
  cerr << " Testing " << _carbon_chain_join.size() << " _carbon_chain_join queries\n";
  for (Substructure_Query* q : _carbon_chain_join) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      cerr << "no match to _carbon_chain_join query\n";
      continue;
    }
    cerr << " got " << sresults.embeddings().size() << " matches\n";
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      Identify2Chain(m, *e, pmd);
    }
  }

  int carbon_chains = 0;
  int atoms_in_carbon_chains = 0;
  int carbon_atom_join_points = 0;

  pmd.ResetVisited();

  const int matoms = m.natoms();

  // When we get here some atoms have been assigned classification 6,
  // but they might be better assigned here. So save all the assignments,
  // reset any that are 6, then when we are done, restore a 6 to any
  // that do not get classified here.
  // No this is not what we want long term....
  std::unique_ptr<int[]> save6 = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    save6[i] = pmd.classification(i);
    if (pmd.classification(i) == 6) {
      pmd.classification(i) = 0;
    }
  }

  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << " b4 chain\n";
  for (int i = 0; i < matoms; ++i) {
    // if (pmd.classification(i) || pmd.visited(i)) {  not sure why this was this way
    if (pmd.classification(i)) {
      // cerr << "Skipping atom " << i << " because classified " << m.smarts_equivalent_for_atom(i) << ' ' << pmd.classification(i) << '\n';
      continue;
    }

    const Atom& a = m[i];
    if (a.atomic_number() != 6 || a.ncon() != 1) {
      continue;
    }

    Set_of_Atoms three_connected;
    pmd.ResetVisited();
    uint32_t ncarbon = IdentifyAllCarbonChainRegion(m, i, pmd, three_connected);
    cerr << m.name() << " from atom " << i << " ncarbon " << ncarbon << " 3con " << three_connected <<  ' ' << m.smarts_equivalent_for_atom(i) << '\n';
    if (ncarbon < _proto->min_atoms_in_all_carbon_chain()) {
      continue;
    }
    ++carbon_chains;
    atoms_in_carbon_chains += ncarbon;
    pmd.increment_carbon_chain_length(ncarbon);
    carbon_atom_join_points += three_connected.number_elements();
    for (atom_number_t a : three_connected) {
      pmd.classification(a) = kCarbonJoin;
    }
    pmd.increment_carbon_chain_count();

    for (int i = 0; i < matoms; ++i) {
      if (! pmd.visited(i)) {
        continue;
      }
      if (pmd.classification(i)) {
        continue;
      }
      if (m.saturated(i)) {
      } else if (pmd.attached_heteroatom_count(i) == 0) {
      } else {
        continue;
      }
      pmd.classification(i) = kCarbonChain;
    }
  }

  pmd.descriptors().set(kCarbonChains, carbon_chains);
  pmd.descriptors().set(kAtomsInCarbonChains, atoms_in_carbon_chains);
  pmd.descriptors().set(kNumberBranches, carbon_atom_join_points);

  // Restore the classification of 6 types if they were not assigned here.
  for (int i = 0; i < matoms; ++i) {
    if (save6[i] == 6 && pmd.classification(i) == 0) {
      pmd.classification(i) = 6;
    }
  }

  pmd.WriteIsotopicallyLabelled(false, cerr);
  cerr << " %carbon_chain\n";
  return 1;
}

// Return true if `zatom` is bonded to a singly connected carbon atom.
// Return false if any bond is non a single bond.
int
BondedToSinglyConnectedCarbon(const Molecule& m, atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      return 0;
    }
    atom_number_t o = b->other(zatom);
    const Atom& ao = m[o];
    if (ao.atomic_number() == 6 && ao.ncon() == 1) {
      return 1;
    }
  }

  return 0;
}

// Return true if `zatom` is a candidate for a 3 connected, fully
// saturated atom.
// Return false if it is bonded to any singly connected atom.
// Return false if any bond is non a single bond.
int
IsThreeConnectedCarbonChainJoin(Molecule& m, atom_number_t zatom) {
  if (m.ring_bond_count(zatom) > 0) {
    return 0;
  }

  const Atom& a = m[zatom];
  if (a.ncon() < 3) {
    return 0;
  }

  // cerr << "IsThreeConnectedCarbonChainJoin checking " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << '\n';
  // cerr << "Checking bonds and connections\n";
  for (const Bond* b : a) {
    if (! b->is_single_bond()) {
      return 0;
    }
    const atom_number_t o = b->other(zatom);
    const Atom& ao = m[o];
    if (ao.atomic_number() != 6) {
      return 0;
    }
    if (ao.ncon() == 1) {
      return 0;
    }
  }

  return 1;
}

int
Options::IdentifyAllCarbonChainRegion(Molecule& m, atom_number_t zatom,
                        PerMoleculeData& pmd, Set_of_Atoms& three_connected) {
  pmd.visited(zatom) = 1;
  // cerr << "IdentifyAllCarbonChainRegion visits atom " << zatom << '\n';

  if (IsThreeConnectedCarbonChainJoin(m, zatom)) {
    three_connected << zatom;
  }

  int rc = 1;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (pmd.visited(o) || pmd.classification(o)) {
      continue;
    }

    if (pmd.attached_heteroatom_count(o) == 0) {
    } else if (m.ncon(o) == 2) {  // end of chain O-CCCCCC
      pmd.visited(o) = 1;
      ++rc;
      continue;
    } else {
      continue;
    }

    if (m.atomic_number(o) == 6) {
    } else {
      continue;
    }

    rc += IdentifyAllCarbonChainRegion(m, o, pmd, three_connected);
  }

  return rc;
}

int
Options::Process(Molecule& m,
                 PerMoleculeData& pmd,
                 IWString_and_File_Descriptor& output) {
  Set_of_Atoms positive_charge;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.formal_charge(i) > 0) {
      positive_charge << i;
    }
  }

  if (positive_charge.empty()) {
    cerr << "No positively charged atoms " << m.name() << "\n";
    ++_uncharged_molecules;
    return 0;
  }

  if (positive_charge.size() != 1) {
    cerr << "Warning " << m.name() << " has multiple positive charges " <<
            positive_charge.size() << '\n';
  }

  for (atom_number_t a : positive_charge) {
    uint32_t hsize = IdentifyHead(m, a, pmd);
    cerr << "Head size " << hsize << '\n';
    pmd.WriteIsotopicallyLabelled(false, cerr);
    cerr << " %head " << hsize << "\n";
    pmd.descriptors().set(kAtomsInHead, hsize);
  }

  if (_break_head_at.size()) {
    MaybeTruncateHead(m, pmd);
  }

  if (_head_linker_tail) {
    return HeadLinkerTail(m, pmd, output);
  }

  MaybeWriteFragments(m, pmd);

  if (_stream_for_descriptors.is_open()) {
    WriteDescriptors(m, pmd);
  }

  return 1;
}

// The head has been identified, the tail has been identified, everything else
// is linker.
int
Options::HeadLinkerTail(Molecule& m, PerMoleculeData& pmd,
                        IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (pmd.classification(i) == kHead) {
      continue;
    }
    if (pmd.classification(i) == kCarbonChain)  {
      continue;
    }
    if (pmd.classification(i) == kCarbonJoin)  {
      continue;
    }
    if (pmd.classification(i) == kDivider) {
      continue;
    }
    if (pmd.classification(i) == kOtherHead) {
      continue;
    }

    pmd.classification(i) = kLinker;
  }

  MaybeWriteFragments(m, pmd);

  pmd.WriteIsotopicallyLabelled(false, cerr);
  cerr << " %headlinkertail\n";

  WriteDescriptors(m, pmd);

  return 1;
}

int
MarkConnectedAtoms(const Molecule& m, atom_number_t zatom, atom_number_t prev, 
                  int flag, int* visited) {
  int rc = 1;
  visited[zatom] = flag;
  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (visited[o] == flag) {
      continue;
    }
    if (o == prev) {
      continue;
    }
    const Atom& ao = m[o];
    if (ao.atomic_number() != 6) {
      continue;
    }
    if (!ao.fully_saturated()) {
      continue;
    }
    rc += MarkConnectedAtoms(m, o, zatom, flag, visited);
  }

  return rc;
}

int
Options::TruncateHead(Molecule& m, atom_number_t zatom, int* dtb,
             PerMoleculeData& pmd) {
  Set_of_Atoms conn;
  m.connections(zatom, conn);
  if (conn.size() != 3) {
    return 0;
  }

  if (m.ring_bond_count(zatom)) {
    return 0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<int[]> myvisited = std::make_unique<int[]>(matoms);

  // We rely on there being 2 connections with a lot of atoms and 1 with a small number.
  int lowest_dtb = matoms;
  int connection_with_lowest_dtb = -1;
  for (int i = 0; i < 3; ++i) {
    std::optional<int> d = m.DownTheBond(zatom, conn[i], dtb, i + 1);
    if (! d) {
      return 0;
    }
    if (*d < lowest_dtb) {
      lowest_dtb = *d;
      connection_with_lowest_dtb = i;
    }
  }

  int connected_carbon = 0;

  for (int i = 0; i < 3; ++i) {
    if (i == connection_with_lowest_dtb) {
      continue;
    }
    std::fill_n(myvisited.get(), matoms, 0);
    myvisited[zatom] = 1;
    if (_label_chains_from_head && m.atomic_number(conn[i]) == 6) {
      connected_carbon += MarkConnectedAtoms(m, conn[i], zatom, kOtherHead, myvisited.get());
      for (int q = 0; q < matoms; ++q) {
        if (myvisited[q] == kOtherHead) {
          pmd.classification(q) = kOtherHead;
        }
      }
    }

    for (int j = 0; j < matoms; ++j) {
      if (dtb[j] != i + 1) {
        continue;
      }

      if (pmd.classification(j) == kCarbonJoin ||
          pmd.classification(j) == kCarbonChain ||
          pmd.classification(j) == kDivider ||
          pmd.classification(j) == kOtherHead ||
          pmd.classification(j) == kBranch) {
      } else {
        pmd.classification(j) = kLinker;
      }
    }
  }

  pmd.descriptors().set(kHeadCarbon, connected_carbon);
  int atoms_in_head = 0;
  for (int i = 0; i < matoms; ++i) {
    if (pmd.classification(i) == kHead) {
      ++atoms_in_head;
    }
  }

  pmd.descriptors().set(kAtomsInTrimmedHead, atoms_in_head);

  return 1;
}

// We have assigned the head region by default rules, but if there is a query
// that truncates the head at a given atom (probably a 3 connected Nitrogen)
// remove those atoms from the head region.
int
Options::MaybeTruncateHead(Molecule& m, PerMoleculeData& pmd) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> dtb = std::make_unique<int[]>(matoms);

  for (Substructure_Query* q : _break_head_at) {
    Substructure_Results sresults;
    if (! q->substructure_search(&m, sresults)) {
      continue;
    }
    atom_number_t n = sresults.embedding(0)->front();
    if (TruncateHead(m, n, dtb.get(), pmd)) {
      return 1;
    }
  }

  return 0;
}

int
Options::MaybeWriteFragments(Molecule& parent, PerMoleculeData& pmd) {
  if (!_stream_for_fragments.is_open()) {
    return 0;
  }

  for (const Bond* b : parent.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (pmd.classification(a1) == kHead && pmd.classification(a2) != kHead) {
    } else if (pmd.classification(a1) != kHead && pmd.classification(a2) == kHead) {
    } else {
      continue;
    }

    parent.set_isotope(a1, 1);
    parent.set_isotope(a2, 1);
  }

  Molecule head;
  parent.create_subset(head, pmd.classification(), kHead);

  Molecule mcopy(parent);
  for (int i = mcopy.natoms() - 1; i >= 0; --i) {
    if (pmd.classification(i) == kHead) {
      mcopy.remove_atom(i);
    }
  }

  resizable_array_p<Molecule> fragments;
  if (! mcopy.create_components(fragments)) {
    fragments << new Molecule(mcopy);
  }

  MaybeSplitTails(fragments);

  IW_STL_Hash_Map<IWString, int> usmi_count;
  for (Molecule* frag : fragments) {
    auto iter = usmi_count.find(frag->unique_smiles());
    if (iter == usmi_count.end()) {
      usmi_count[frag->unique_smiles()] = 1;
    } else {
      ++iter->second;
    }
  }

  pmd.descriptors().set(kNumberFragments, fragments.size());
  pmd.descriptors().set(kNumberUniqueFragments, usmi_count.size());

  static constexpr char kSep = ' ';
  _stream_for_fragments << parent.smiles() << kSep << parent.name() << kSep <<
                           "parent" << kSep << fragments.size() << kSep << usmi_count.size() << '\n';
  pmd.WriteIsotopicallyLabelled(false, _stream_for_fragments);
  _stream_for_fragments << kSep << parent.name() << " parent_iso\n";

  _stream_for_fragments << head.smiles() << kSep << parent.name() << " head\n";

  if (fragments.empty()) {
    _stream_for_fragments << mcopy.smiles() << kSep << parent.name() << kSep << "frag" << kSep << 1 << '\n';
    return 1;
  }

  for (const auto& [usmi, count] : usmi_count) {
    _stream_for_fragments << usmi << kSep << parent.name() << kSep << "frag" << kSep << count << '\n';
  }

  return 1;
}

// if `zatom` is a saturated, 2 connected atom, create fragments by
// breaking each of the bonds attached to `zatom` and adding those
// fragments to `destination`. Note that `zatom` will be included with
// both frgaments added to `destination`.
int
MakeBothSides(Molecule& m,
              atom_number_t zatom,
              resizable_array_p<Molecule>& destination) {
  cerr << "MakeBothSides " << m.smiles() << '\n';
  if (m.ncon(zatom) != 2) {
    return 0;
  }

  Set_of_Atoms conn;
  conn.reserve(2);

  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      return 0;
    }
    conn << b->other(zatom);
  }

  if (conn.size() != 2) {
    return 0;
  }

  for (atom_number_t c : conn) {
    std::unique_ptr<Molecule> lhs = std::make_unique<Molecule>(m);
    lhs->remove_bond_between_atoms(zatom, c);
    lhs->remove_fragment_containing_atom(c);
    destination << lhs.release();
  }

  return 1;
}

int
Options::MaybeSplitTails(resizable_array_p<Molecule>& fragments) {
  cerr << " have " << _break_tail_at_attachment.size() << " _break_tail_at_attachment queries\n";
  if (_break_tail_at_attachment.empty()) {
    return 0;
  }

  int rc = 0;
  for (int i = fragments.number_elements() - 1; i >= 0; --i) {
    for (Substructure_Query* q : _break_tail_at_attachment) {
      Substructure_Results sresults;
      if (! q->substructure_search(fragments[i], sresults)) {
        continue;
      }
      // Don't check for number_embeddings, look for one that works.
      atom_number_t break_point = sresults.embedding(0)->front();

      if (MakeBothSides(*fragments[i], break_point, fragments))  {
        fragments.remove_item(i);
        ++rc;
      }
    }
  }

  return rc;
}

int
Options::WriteDescriptors(Molecule& m, const PerMoleculeData& pmd) {
  _stream_for_descriptors << m.name();

  const float* values = pmd.descriptors().values();

  for (int i = 0; i < kNextToAssign; ++i) {
    _stream_for_descriptors << _output_separator << values[i];
  }
  _stream_for_descriptors << '\n';

  _stream_for_descriptors.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
IlpdReducedGraph(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
IlpdReducedGraph(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! IlpdReducedGraph(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
IlpdReducedGraph(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "IlpdReducedGraph:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return IlpdReducedGraph(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:H:N:T:A:lcg:i:C:S:F:z:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  int ignore_errors = 0;
  int errors_ignored = 0;
  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == "i") {
        ignore_errors = 0;
        if (verbose) {
          cerr << "Will ignore otherwise fatal errors\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 1;
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (IlpdReducedGraph(options, fname, input_type, output)) {
    } else if (ignore_errors) {
      ++errors_ignored;
    } else {
      cerr << "IlpdReducedGraph::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    if (ignore_errors) {
      cerr << errors_ignored << " otherwise fatal errors ignored\n";
    }
    options.Report(cerr);
  }

  return 0;
}

}  // namespace ilpd_reduced_graph

int
main(int argc, char ** argv) {

  int rc = ilpd_reduced_graph::Main(argc, argv);

  return rc;
}
