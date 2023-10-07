// Make minor changes to incoming molecules.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>

#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/minor_changes.h"

#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools/minor_changes.pb.h"

namespace minor_changes {

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
  // This is redundant with the function in minor_changes_main.cc, consolidate...
  cerr << "Makes a potentially large number of small changes to a molecule.\n";
  cerr << " -C <fname>    textproto file with MinorChangesData options\n";
  cerr << " -F <fname>    one or more DicerFragment textproto files (from get_substituents)\n";
  cerr << " -P <atype>    atom typing specification -needed if atom types in the get_substituents data\n";
  cerr << " -u <support>  support level for inclusion from fragment libraries\n";
  cerr << " -M <max>      maximum number of products per starting molecule\n";
  cerr << " -v            verbose output\n";
// clang-format on

  ::exit(rc);
}

Fragment::Fragment() {
  _attach = INVALID_ATOM_NUMBER;
  _number_exemplars = 0;

  _atype = 0;
}

// get_substituents writes protos with the key as well as
// the proto, so building must remove the first token from
// `buffer`.
int
Fragment::Build(const_IWSubstring buffer, uint32_t min_support) {
  if (buffer.empty()) {
    return 0;
  }

  buffer.remove_leading_words(1);

  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());

  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&input, &proto)) {
    cerr << "Fragment::Build:invalid input '" << buffer << "'\n";
    return 0;
  }

  if (min_support > 0 && proto.n() < min_support) {
    return 0;
  }

  return BuildFromSmiles(proto);
}

// Build the molecule and attachment points from `proto`.
// We should be looking at the `iso` field in `proto`.
// Note too that Sulphur containing fragments can be problematic.
// As a fragment, the S may appear to have no hydrogens, but it
// would be quite happy being attached, thereby adopting an alternate
// valence. S(=O)(=NC(=O)OCC)(N)C1=CC=CC=C1 CHEMBL1556501
int
Fragment::BuildFromSmiles(const dicer_data::DicerFragment& proto) {
  if (! _m.build_from_smiles(proto.smi())) {
    cerr << "Fragment::BuildFromSmiles:cannot parse smiles '" << proto.smi() << "'\n";
    return 0;
  }

  _m.set_name(proto.par());
  _number_exemplars = proto.n();

  // If there is an isotope, use that as the atom type.

  for (int i = 0; i < _m.natoms(); ++i) {
    isotope_t iso = _m.isotope(i);
    if (iso == 0) {
      continue;
    }
    // Must have an implicit Hydrogen in order to join.
    if (_m.hcount(i) == 0) {
      continue;
    }

    if (_attach == INVALID_ATOM_NUMBER) {
      _attach = i;
      _atype = iso;
    } else {
      cerr << "Fragment::BuildFromSmiles:more than 1 isotope '" << _m.smiles() << "'\n";
      return 0;
    }
  }

  if (_attach != INVALID_ATOM_NUMBER) {
    return 1;
  }

  // No isotopes, the first atom with an implicit Hydrogen.

  for (int i = 0; i < _m.natoms(); ++i) {
    if (_m.hcount(i)) {
      _attach = i;
      return 1;
    }
  }

  cerr << "Fragment::BuildFromSmiles:no attachment points in '" << proto.smi() << "'\n";
  return 0;
}

int
Fragment::OkForJoining() const {
  return 1;
}

BivalentFragment::BivalentFragment() {
  _attach1 = INVALID_ATOM_NUMBER;
  _attach2 = INVALID_ATOM_NUMBER;

  _bonds_between = -1;

  _atype1 = 0;
  _atype2 = 0;
  _number_exemplars = 0;
}

int
BivalentFragment::Build(const_IWSubstring buffer, uint32_t min_support) {
  if (buffer.empty()) {
    return 0;
  }

  buffer.remove_leading_words(1);

  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());

  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&input, &proto)) {
    cerr << "Fragment::Build:invalid input '" << buffer << "'\n";
    return 0;
  }

  if (min_support > 0 && proto.n() < min_support) {
    return 0;
  }

  return BuildFromSmiles(proto);
}

// Build the molecule and attachment points from `proto`.
// We should be looking at the `iso` field in `proto`.
int
BivalentFragment::BuildFromSmiles(const dicer_data::DicerFragment& proto) {
  if (! _m.build_from_smiles(proto.smi())) {
    cerr << "Fragment::BuildFromSmiles:cannot parse smiles '" << proto.smi() << "'\n";
    return 0;
  }

  _m.set_name(proto.par());
  _number_exemplars = proto.n();

  const int matoms = _m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso = _m.isotope(i);
    if (iso == 0) {
      continue;
    }

    // Must have an implicit Hydrogen in order to join.
    if (_m.hcount(i) == 0) {
      continue;
    }

    if (_attach1 == INVALID_ATOM_NUMBER) {
      _attach1 = i;
      _atype1 = iso;
    } else if (_attach2 == INVALID_ATOM_NUMBER) {
      _attach2 = i;
      _atype2 = iso;
    } else {
      cerr << "BivalentFragment::BuildFromSmiles:too many isotopes " << proto.smi() << "'\n";
      return 0;
    }
  }

  if (_attach1 == INVALID_ATOM_NUMBER) {
    cerr << "BivalentFragment::BuildFromSmiles:no isotopes '" << proto.smi() << "'\n";
    return 0;
  }

  // If only 1 isotope, this will be a simple insertion, with maybe
  // a sidechain. But now this join point will need to accept two
  // bonds, so there must be at least two implicit hydrogens.
  if (_attach2 == INVALID_ATOM_NUMBER) {
    if (_m.hcount(_attach1) < 2) {
      return 0;
    }

    _attach2 = _attach1;
    _atype2 = _atype1;
    _bonds_between = 0;
  } else {
    _bonds_between = _m.bonds_between(_attach1, _attach2);
  }

  return 1;
}

using fixed_bit_vector::FixedBitVector;


// A class that holds information frequently used during computations.
class MoleculeData {
  private:
    // the number of atoms in the molecule.
    int _matoms;

    // The atomic numbers.
    atomic_number_t* _z;

    // Is the atom unsaturated (nbonds - ncon);
    int* _unsaturated;
    // Is the atom bonded to a heteroatom.
    int* _attached_heteroatom_count;
    // Is the atom a 5 valent nitrogen.
    int* _five_valent_nitrogen;
    // Save time by storing the hydrogen count.
    int* _hcount;
    // Save time by storing the fused system size.
    int* _fused_system_size;

    // For each atom, the atomic numbers of atoms that are alpha
    // and beta.
    FixedBitVector* _alpha;
    FixedBitVector* _beta;

    // If the proto config specifies sizes of fragments that can
    // be lost, those will be transferred to this array.
    int* _remove_fragment;

    // Maybe there are atom types.
    uint32_t* _atype;

    // If queries are present for specifying atoms that can change.
    int* _process_these_atoms;

  public:
    MoleculeData(Molecule& m);
    ~MoleculeData();

    int StoreAtomTypes(Molecule& m, Atom_Typing_Specification& atype_spec);

    // If `proto` contains remove_fragment entries, copy those to _remove_fragment.
    int StoreLoseFragment(const minor_changes_data::MinorChangesData& proto);

    uint32_t atype(atom_number_t zatom) const {
      return _atype[zatom];
    }

    // return true if atom `zatom` is alpha to atomic number `atn`.
    int alpha(atom_number_t zatom, atomic_number_t atn) const {
      return _alpha[zatom].is_set(atn);
    }
    int beta(atom_number_t zatom, atomic_number_t atn) const {
      return _beta[zatom].is_set(atn);
    }

    int unsaturated(atom_number_t zatom) const {
      return _unsaturated[zatom];
    }

    int attached_heteroatom_count(atom_number_t zatom) const {
      return _attached_heteroatom_count[zatom];
    }

    int five_valent_nitrogen(atom_number_t zatom) const {
      return _five_valent_nitrogen[zatom];
    }

    int hcount(atom_number_t zatom) const {
      return _hcount[zatom];
    }

    int fused_system_size(atom_number_t zatom) const {
      return _fused_system_size[zatom];
    }

    int remove_fragment(int frag_size) {
      return _remove_fragment[frag_size];
    }

    int* process_these_atoms() {
      return _process_these_atoms;
    }

    int process_atom(atom_number_t a) {
      return _process_these_atoms[a];
    }

    bool CanChange(atom_number_t a1) const {
      return _process_these_atoms[a1];
    }

    // Returns true if _process_these_atoms is set for both `a1` and `a2`.
    bool CanChange(atom_number_t a1, atom_number_t a2) const {
      return _process_these_atoms[a1] && _process_these_atoms[a2];
    }

    bool AllAtomsCanChange(const Set_of_Atoms& s) const {
      return s.all_members_set_in_array(_process_these_atoms, 1);
    }
    
    // True if every item set in `values` is also set in `_process_these_atoms`.
    bool AllAtomsCanChange(const int* values) const;
};

// Set bits for the atomic numbers alpha and beta to `zatom`.
void
FillAlphaBeta(const Molecule& m,
         atom_number_t zatom,
         FixedBitVector* alpha,
         FixedBitVector* beta,
         int* attached_heteroatom_count)  {
  for (const Bond* b1 : m.atom(zatom)) {
    atom_number_t j = b1->other(zatom);
    alpha[zatom].set_bit(m.atomic_number(j));
    if (m.atomic_number(j) != 6) {
      ++attached_heteroatom_count[zatom];
    }
    for (const Bond* b2 : m.atom(j)) {
      atom_number_t k = b2->other(j);
      if (k == zatom) {
        continue;
      }
      beta[zatom].set_bit(m.atomic_number(k));
    }
  }
}

MoleculeData::MoleculeData(Molecule& m) {
  _matoms = m.natoms();

  _z = new atomic_number_t[_matoms];
  _unsaturated = new int[_matoms];
  _attached_heteroatom_count = new_int(_matoms);
  _five_valent_nitrogen = new_int(_matoms);
  _hcount = new_int(_matoms);
  _fused_system_size = new_int(_matoms);
  _alpha = new FixedBitVector[_matoms];
  _beta = new FixedBitVector[_matoms];
  _process_these_atoms = new_int(_matoms, 1);  // initialised so all processed.

  // 128 leaves room for all known atomic numbers.
  for (int i = 0; i < _matoms; ++i) {
    _alpha[i].resize(128);
    _beta[i].resize(128);
  }

  for (int i = 0; i < _matoms; ++i) {
    const Atom& a = m.atom(i);
    _z[i] = a.atomic_number();
    _unsaturated[i] = m.nbonds(i) - a.ncon();
    _hcount[i] = m.hcount(i);
    _fused_system_size[i] = m.fused_system_size(i);
    if (_z[i] == 7 && a.formal_charge() == 0 && m.nbonds(i) > 4) {
      _five_valent_nitrogen[i] = 1;
    }
    FillAlphaBeta(m, i, _alpha, _beta, _attached_heteroatom_count);
  }

  _remove_fragment = nullptr;
  _atype = nullptr;
}

MoleculeData::~MoleculeData() {
  delete [] _z;
  delete [] _unsaturated;
  delete [] _attached_heteroatom_count;
  delete [] _five_valent_nitrogen;
  delete [] _hcount;
  delete [] _alpha;
  delete [] _beta;
  delete [] _process_these_atoms;
  delete [] _fused_system_size;

  if (_atype != nullptr) {
    delete [] _atype;
  }
  if (_remove_fragment != nullptr) {
    delete [] _remove_fragment;
  }
}

int
MoleculeData::StoreLoseFragment(const minor_changes_data::MinorChangesData& proto) {
  if (proto.remove_fragment_size() == 0) {
    return 0;
  }

  _remove_fragment = new_int(_matoms);
  for (int f : proto.remove_fragment()) {
    if (f < _matoms) {
      _remove_fragment[f] = 1;
    }
  }

  return 1;
}

int
MoleculeData::StoreAtomTypes(Molecule& m,
        Atom_Typing_Specification& atype_spec) {
  _atype = new uint32_t[m.natoms()];

  return atype_spec.assign_atom_types(m, _atype);
}

bool
MoleculeData::AllAtomsCanChange(const int* values) const {
  for (int i = 0; i < _matoms; ++i) {
    if (values[i] && _process_these_atoms[i] == 0) {
      return false;
    }
  }

  return true;
}

Options::Options() {
  _verbose = 0;
  _first_call = true;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _molecules_generated = 0;
  _invalid_valence = 0;
  _remove_chirality = 0;
  _acc = new Accumulator_Int<uint32_t>[TransformationType::kHighest];
}

Options::~Options() {
  delete [] _acc;
}

// Sort `lib` by number of exemplars, and if larger than `max_size`
// truncate to that length.
template <typename F>
void
SortAndTrim(resizable_array_p<F>& lib, uint32_t max_size) {
  lib.iwqsort_lambda([](const F* f1, const F* f2) {
    if (f1->number_exemplars() < f2->number_exemplars()) {
      return 1;
    }
    if (f1->number_exemplars() > f2->number_exemplars()) {
      return -1;
    }

    return 0;
  });

  if (lib.size() <= max_size) {
    return;
  }

  lib.resize(max_size);

  return;
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
  } else {
    _reduce_to_largest_fragment = 1;
    cerr << "Default is to reduce to largest fragment\n";
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  // Parse the proto first, values will then be overwritten by command
  // line values.
  if (cl.option_present('C')) {
    const_IWSubstring c = cl.string_value('C');
    if (! ReadOptions(c)) {
      cerr << "Cannot read MinorChangesData proto from '" << c << "'\n";
      return 0;
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! _atom_typing.build(p)) {
      cerr << "Invalid atom typing directive '" << p << "'\n";
      return 0;
    }
  }

  // Parse support level before reading fragments.
  if (cl.option_present('u')) {
    uint32_t support;
    if (! cl.value('u', support)) {
      cerr << "Invalid support specification (-u)\n";
      return 0;
    }

    _config.set_fragment_support(support);
    if (_verbose) {
      cerr << "Will only include fragments with at least " << support << " instances\n";
    }
  }

  if (cl.option_present('M')) {
    uint32_t max_variants;
    if (! cl.value('M', max_variants)) {
      cerr << "Invalid max products (-M)\n";
      return 0;
    }
    _config.set_max_variants(max_variants);
  }

  if (cl.option_present('F')) {
    const_IWSubstring f;
    for (int i = 0; cl.value('F', f, i); ++i) {
      if (! ReadFragments(f)) {
        cerr << "Cannot read fragments from '" << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Read " << _fragment.size() << " fragments\n";
    }

    if (_remove_chirality) {
      for (Fragment* f : _fragment) {
        f->mol().remove_all_chiral_centres();
      }
    }

    if (_fragment.empty()) {
      cerr << "Options::Initialise:no fragments read\n";
      return 0;
    }

    if (_config.max_fragment_lib_size() > 0) {
      SortAndTrim(_fragment, _config.max_fragment_lib_size());

      if (_verbose) {
        cerr << "After SortAndTrim have " << _fragment.size() << " fragments\n";
      }
    }
  }

  if (cl.option_present('B')) {
    const_IWSubstring b;
    for (int i = 0; cl.value('B', b, i); ++i) {
      if (! ReadBivalentFragments(b)) {
        cerr << "Cannot read bivalent fragments from '" << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Read " << _bivalent_fragment.size() << " bivalent fragments\n";
    }

    if (_bivalent_fragment.empty()) {
      cerr << "Options::Initialise:no bivalent fragments read\n";
      return 0;
    }

    if (_remove_chirality) {
      for (BivalentFragment* f : _bivalent_fragment) {
        f->mol().remove_all_chiral_centres();
      }
    }

    if (_config.max_bivalent_fragment_lib_size() > 0) {
      SortAndTrim(_bivalent_fragment, _config.max_bivalent_fragment_lib_size());

      if (_verbose) {
        cerr << "After SortAndTrim have " << _bivalent_fragment.size() << " bivalent fragments\n";
      }
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (!q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }

      _only_process_query << q.release();
    }
  }

  if (cl.option_present('q')) {
    static constexpr int kFlag = ' ';  // not used.
    if (!process_queries(cl, _only_process_query, _verbose, kFlag)) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  return 1;
}

// Return 1 if any of the optional features have been enabled.
int
Options::AnythingSpecified() const {
  if (_config.has_add_fragments()) { return 1;}  // 2
  if (_config.has_replace_terminal_fragments()) { return 1;}  // 3
  if (_config.has_single_to_double_bond()) { return 1;}  // 4
  if (_config.has_double_to_single_bond()) { return 1;}  // 5
  if (_config.has_unspiro()) { return 1;}  // 6
  if (_config.has_make_three_membered_rings()) { return 1;}  // 7
  if (_config.has_change_carbon_to_nitrogen()) { return 1;}  // 8
  if (_config.has_change_carbon_to_oxygen()) { return 1;}  // 9
  if (_config.has_change_nitrogen_to_carbon()) { return 1;}  // 10
  if (_config.has_insert_ch2()) { return 1;}  // 11
  if (_config.has_remove_ch2()) { return 1;}  // 12
  if (_config.has_destroy_aromatic_rings()) { return 1;}  // 13
  if (_config.has_destroy_aromatic_ring_systems()) { return 1;}  // 14
  if (_config.has_swap_adjacent_atoms()) { return 1;}  // 15
  if (_config.remove_fragment().size()) { return 1;}   // 16
  if (_config.has_insert_fragments()) { return 1;}  // 23
  if (_config.has_replace_inner_fragments()) { return 1;}  // 24

  return 0;
}

// Either _config has been reset, or this is the first run with
// default settings. Make sure settings are OK.
int
Options::CheckConditions() {
  if (_config.max_variants() == 0) {
    _config.set_max_variants(std::numeric_limits<uint32_t>::max());
  }

  if ((_config.fragment_size() && _fragment.empty()) ||
       (_config.bivalent_fragment_size() && _bivalent_fragment.empty())) {
    if (! ReadFragmentsFromConfig()) {
      return 0;
    }
  }

  if (_config.only_process_query_size()) {
    if (! SetupOnlyProcessQueries()) {
      cerr << "Options::CheckConditions:process only queries bad\n";
      cerr << _config.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (_config.has_atype()) {
    const IWString tmp(_config.atype());
    if (! _atom_typing.build(tmp)) {
      cerr << "Options::CheckConditions:cannot initialise atom typing '" << _config.atype() << "'\n";
      return 0;
    }
  }

  if (_config.fragment_support() > 0 &&
      (_config.max_fragment_lib_size() > 0 || _config.max_bivalent_fragment_lib_size())) {
    cerr << "Options::CheckConditions:cannot specify support level and fragment library size\n";
    return 0;
  }

  if (AnythingSpecified()) {
    return 1;
  }

  // If nothing has been specified, turn on everything.

  // only if we have fragments.
  _config.set_add_fragments(true);  // 2

  _config.set_replace_terminal_fragments(true);  // 3
  _config.set_replace_inner_fragments(true);  // 24
  // Set a reasonable default for max atoms lost.
  if (! _config.has_max_atoms_lost()) {
    _config.set_max_atoms_lost(3);  // 19
  }

  _config.set_single_to_double_bond(true);  // 4
  _config.set_double_to_single_bond(true);  // 5
  _config.set_unspiro(true);  // 6
  _config.set_make_three_membered_rings(true);  // 7
  _config.set_change_carbon_to_nitrogen(true);  // 8
  _config.set_change_carbon_to_oxygen(true);  // 9
  _config.set_change_nitrogen_to_carbon(true); // 10
  _config.set_insert_ch2(true); // 11
  _config.set_remove_ch2(true); // 12
  _config.set_destroy_aromatic_rings(true); // 13
  _config.set_destroy_aromatic_ring_systems(true); // 14
  _config.set_swap_adjacent_atoms(true); // 15
  _config.set_insert_fragments(true); // 23

  return 1;
}

// Return true if results.size > max_variants
int
Options::GeneratedEnough(const resizable_array_p<Molecule>& results) const {
  return results.size() >= _config.max_variants();
}

// There are fragment smiles in _config.fragment().
// Transfer to the _fragment array.
int
Options::ReadFragmentsFromConfig() {
  for (const std::string& smi : _config.fragment()) {
    dicer_data::DicerFragment proto;
    proto.set_smi(smi);
    if (_atom_typing.active()) {
      proto.set_iso(dicer_data::ATYPE);
    } else {
      proto.set_iso(dicer_data::ATT);
    }

    std::unique_ptr<Fragment> frag = std::make_unique<Fragment>();
    if (! frag->BuildFromSmiles(proto)) {
      cerr << "Options::ReadFragmentsFromConfig:cannot parse '" << smi << "'\n";
      return 0;
    }
    _fragment << frag.release();
  }

  for (const std::string& smi : _config.bivalent_fragment()) {
    dicer_data::DicerFragment proto;
    proto.set_smi(smi);
    if (_atom_typing.active()) {
      proto.set_iso(dicer_data::ATYPE);
    } else {
      proto.set_iso(dicer_data::ATT);
    }

    std::unique_ptr<BivalentFragment> frag = std::make_unique<BivalentFragment>();
    if (! frag->BuildFromSmiles(proto)) {
      cerr << "Options::ReadFragmentsFromConfig:cannot parse bivalent '" << smi << "'\n";
      return 0;
    }
    _bivalent_fragment << frag.release();
  }

  if (_config.max_fragment_lib_size() > 0) {
    SortAndTrim(_fragment, _config.max_fragment_lib_size());
    if (_verbose) {
      cerr << "After SortAndTrim have " << _fragment.size() << " fragments\n";
    }
  }
  if (_config.max_bivalent_fragment_lib_size() > 0) {
    SortAndTrim(_bivalent_fragment, _config.max_bivalent_fragment_lib_size());
    if (_verbose) {
      cerr << "After SortAndTrim have " << _bivalent_fragment.size() << " bivalent fragments\n";
    }
  }

  return 1;
}

int
Options::SetupOnlyProcessQueries() {
  for (const std::string& q : _config.only_process_query()) {
    if (! SetupOnlyProcessQuery(q)) {
      cerr << "Options::SetupOnlyProcessQueries:invalid query specification '" << q << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::SetupOnlyProcessQuery(const std::string& qry) {
  static constexpr char kOption = ' ';   // not used.
  static constexpr int kVerbose = 0;     // not used.
  const_IWSubstring tmp(qry);
  return process_cmdline_token(kOption, tmp, _only_process_query, kVerbose);
}

int
Options::SetConfig(const minor_changes_data::MinorChangesData& proto) {
  _config = proto;

  _first_call = false;

  return CheckConditions();
}

int
Options::ReadFragments(const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFragments:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragments(input);
}

int
Options::ReadFragments(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.empty() || buffer.starts_with('#')) {
      continue;
    }

    if (! ReadFragment(buffer)) {
      cerr << "Options::ReadFragments:cannot process '" << buffer << "'\n";
      cerr << "Ignored\n";
    }
  }

  return _fragment.size();
}

int
Options::ReadBivalentFragments(const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadBivalentFragments:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadBivalentFragments(input);
}

int
Options::ReadBivalentFragments(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.empty() || buffer.starts_with('#')) {
      continue;
    }

    if (! ReadBivalentFragment(buffer)) {
      cerr << "Options::ReadBivalentFragments:cannot process '" << buffer << "'\n";
      cerr << "ignored\n";  // might just be a support failure.
      continue;
    }
  }

  return _bivalent_fragment.size();
}

int
Options::ReadFragment(const_IWSubstring& buffer) {
  std::unique_ptr<Fragment> frag = std::make_unique<Fragment>();
  if (! frag->Build(buffer, _config.fragment_support())) {
    return 0;
  }

  _fragment << frag.release();

  return 1;
}

int
Options::ReadBivalentFragment(const_IWSubstring& buffer) {
  std::unique_ptr<BivalentFragment> frag = std::make_unique<BivalentFragment>();
  if (! frag->Build(buffer, _config.fragment_support())) {
    return 0;
  }

  _bivalent_fragment << frag.release();

  return 1;
}

#ifdef NOT_USED_WANT_COMMENTS
int
Options::ReadOptions(const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadOptions(input);
}

int
Options::ReadOptions(iwstring_data_source& input) {

  // Copied from proto_support.h

  using google::protobuf::io::ZeroCopyInputStream;
  using google::protobuf::io::FileInputStream;
  std::unique_ptr<FileInputStream> zero_copy_input(new FileInputStream(input.fd()));

  if (! google::protobuf::TextFormat::Parse(zero_copy_input.get(), &_config)) {
    cerr << "Options::ReadOptions:cannot read '" << input.fname() << "'\n";
    return 0;
  }

  return 1;
}
#endif

int
Options::ReadOptions(const_IWSubstring& fname) {
  IWString tmp(fname);
  std::optional<minor_changes_data::MinorChangesData> proto =
      iwmisc::ReadTextProtoCommentsOK<minor_changes_data::MinorChangesData>(tmp);
  if (! proto) {
    cerr << "optional::ReadOptions:cannot read '" << fname << "'\n";
    return 0;
  }

  _config = std::move(*proto);

  return 1;
}

// Return true if we have already encountered `m`.
int
Options::Seen(Molecule& m) {
  if (_seen.contains(m.unique_smiles())) {
    return 1;
  }

  _seen.insert(m.unique_smiles());

  return 0;
}

// If this is an ok new molecule, add it to `results`, and release
// the pointer.
// Returns 1 if `m` was added to `results`.
int
Options::AddToResultsIfNew(std::unique_ptr<Molecule>& m,
                           resizable_array_p<Molecule>& results) {
  ++_molecules_generated;

  if (! m->valence_ok()) {
    ++_invalid_valence;
    if (_config.echo_bad_valence()) {
      cerr << "Bad valence " << m->smiles() << "\n";
    }
    return 0;
  }

  if (Seen(*m)) {
    return 0;
  }

  results << m.release();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules, generated " << _molecules_generated << '\n';
  output << _invalid_valence << " molecules with invalid valence discarded\n";
  // Other information about what has happened.

  Accumulator_Int<uint32_t> acc;
  for (int i = 0; i < _variants_generated.number_elements(); ++i) {
    if (_variants_generated[i]) {
      output << _variants_generated[i] << " molecules generated " << i << " variants\n";
      acc.extra(i, _variants_generated[i]);
    }
  }

  output << "Btw " << acc.minval() << " and " << acc.maxval() << " mean " << acc.average() << '\n';

  std::unordered_map<TransformationType, std::string> name {
    { TransformationType::kAddFragments, "Add Fragments"},
    { TransformationType::kReplaceTerminalFragments, "Replace Terminal Fragments"},
    { TransformationType::kSingleToDoubleBond, "Single to Double"},
    { TransformationType::kDoubleToSingleBond, "Double to Single"},
    { TransformationType::kUnspiro, "Unspiro"},
    { TransformationType::kMakeThreeMemberedRings, "Make 3 ring"},
    { TransformationType::kChangeCarbonToNitrogen, "C->N"},
    { TransformationType::kChangeCarbonToOxygen, "C->O"},
    { TransformationType::kChangeNitrogenToCarbon, "N->C"},
    { TransformationType::kInsertCH2, "Insert CH2"},
    { TransformationType::kRemoveCH2, "Remove CH2"},
    { TransformationType::kDestroyAromaticRings, "Destroy arom ring"},
    { TransformationType::kDestroyAromaticRingSystems, "Destroy arom ring system"},
    { TransformationType::kSwapAdjacentAtoms, "Swap adjacent"},
    { TransformationType::kRemoveFragment, "Remove Fragment"},
    { TransformationType::kInsertFragments, "Insert Fragments"},
    { TransformationType::kReplaceInnerFragments, "Replace Inner Fragments"},
    { TransformationType::kFuseBiphenyls, "Fuse biphenyls"}
  };

#ifdef NOT_USED_ASDASDASD
  resizable_array<TransformationType> types {
    TransformationType::kAddFragments,
    TransformationType::kReplaceTerminalFragments,
    TransformationType::kSingleToDoubleBond,
    TransformationType::kDoubleToSingleBond,
    TransformationType::kUnspiro,
    TransformationType::kMakeThreeMemberedRings,
    TransformationType::kChangeCarbonToNitrogen,
    TransformationType::kChangeCarbonToOxygen,
    TransformationType::kChangeNitrogenToCarbon,
    TransformationType::kInsertCH2,
    TransformationType::kRemoveCH2,
    TransformationType::kDestroyAromaticRings,
    TransformationType::kDestroyAromaticRingSystems,
    TransformationType::kSwapAdjacentAtoms,
    TransformationType::kRemoveFragment,
    TransformationType::kInsertFragments,
    TransformationType::kReplaceInnerFragments
  };
#endif

  for (const auto& iter : name) {
    const Accumulator_Int<uint32_t>& acc = _acc[iter.first];
    output << iter.second << ' ' << acc.n() << " btw " << acc.minval()
           << ' ' << acc.maxval() << " tot " << acc.sum();
    if (acc.n() > 0) {
      output << " ave " << ' ' << static_cast<float>(acc.average());
    }

    output << '\n';
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

int
Options::Process(Molecule& m,
                 resizable_array_p<Molecule>& results) {
  if (_first_call) {
    if (! CheckConditions()) {
      return 0;
    }

    _first_call = false;
  }

  results.reserve(1000);  // arbitrary guess.

  ++_molecules_read;

  MoleculeData molecule_data(m);
  molecule_data.StoreLoseFragment(_config);

  if (_atom_typing.active()) {
    molecule_data.StoreAtomTypes(m, _atom_typing);
  }

  if (_only_process_query.size()) {
    if (! DetermineChangingAtoms(m, molecule_data)) {
      return 0;
    }
  }

  // Add the parent to the hash of what has been generated.
  Seen(m);

  const int rc = Process(m, molecule_data, results);

  ++_variants_generated[results.size()];

  return rc;
}

// Run the _only_process_query queries on `m` and store
// the results in `molecule_data`
int
Options::DetermineChangingAtoms(Molecule& m,
                                MoleculeData& molecule_data) {
  // Turn off all atoms. To be selectively enabled below.
  std::fill_n(molecule_data.process_these_atoms(), m.natoms(), 0);

  Molecule_to_Match target(&m);

  int rc = 0;
  for (Substructure_Query* q : _only_process_query) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    sresults.each_embedding_set_vector(molecule_data.process_these_atoms(), 1);
    ++rc;
  }

  return rc;
}

int
Options::Process(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {

  uint32_t rc = 0;

  rc += SingleToDoubleBond(m, molecule_data, results);  // 4
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "SingleToDoubleBond " << rc << '\n';
  rc += DoubleToSingleBond(m, molecule_data, results);  // 5
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "DoubleToSingleBond " << rc << '\n';
  rc += Unspiro(m, molecule_data, results);  // 6
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "Unspiro " << rc << '\n';
  rc += MakeThreeRing(m, molecule_data, results);  // 7
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "MakeThreeRing " << rc << '\n';
  rc += ChangeCarbonToNitrogen(m , molecule_data, results);  // 8
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "ChangeCarbonToNitrogen " << rc << '\n';
  rc += ChangeCarbonToOxygen(m , molecule_data, results);  // 9
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "ChangeCarbonToOxygen " << rc << '\n';
  rc += ChangeNitrogenToCarbon(m , molecule_data, results);  // 10
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "ChangeNitrogenToCarbon " << rc << '\n';
  rc += InsertCh2(m, molecule_data, results);  // 11
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "InsertCh2 " << rc << '\n';
  rc += RemoveCh2(m, molecule_data, results);  // 12
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "RemoveCh2 " << rc << '\n';
  rc += DestroyAromaticRings(m, molecule_data, results);  // 13
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "DestroyAromaticRings " << rc << '\n';
  rc += DestroyAromaticRingSystems(m, molecule_data, results);  // 14
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "DestroyAromaticRingSystems " << rc << '\n';
  rc += SwapAdjacentAtoms(m, molecule_data, results);  // 15

  rc += ReplaceTerminalFragments(m, molecule_data, results);  // 3
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "ReplaceTerminalFragments " << rc << '\n';

  rc += RemoveFragments(m, molecule_data, results);  // 3
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "RemoveFragment " << rc << '\n';

  rc += InsertBivalentFragments(m, molecule_data, results);  // 23
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "InsertBivalentFragments " << rc << '\n';

  rc += ReplaceInnerFragments(m, molecule_data, results);  // 24
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "ReplaceInnerFragments " << rc << '\n';

  rc += AddFragments(m, molecule_data, results);  // 2
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "AddFragments " << rc << '\n';

  rc += FuseBiphenyls(m, molecule_data, results);  // 2
  if (rc > _config.max_variants()) {
    return rc;
  }
  // q cerr << "FuseBiphenyls " << rc << '\n';

  // q cerr << "Options::Process:returning " << rc << '\n';

  return rc;
}

int
Options::SingleToDoubleBond(Molecule& m,
                            const MoleculeData& molecule_data,
                            resizable_array_p<Molecule>& results) {
  if (! _config.single_to_double_bond()) {
    return 0;
  }

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (const Bond* b : m.bond_list()) {
    if (! b->is_single_bond() || b->is_aromatic()) {
      continue;
    }
    atom_number_t a1 = b->a1();
    if (m.hcount(a1) == 0) {
      continue;
    }
    atom_number_t a2 = b->a2();
    if (m.hcount(a2) == 0) {
      continue;
    }

    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    if (molecule_data.unsaturated(a1) || molecule_data.unsaturated(a2)) {
      continue;
    }

    // cerr << a1 << ' ' << m.smarts_equivalent_for_atom(a1) << " to " << a2 << ' ' << m.smarts_equivalent_for_atom(a2) << '\n';

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    mcopy->set_bond_type_between_atoms(a1, a2, DOUBLE_BOND);
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kSingleToDoubleBond].extra(rc);

  return rc;
}

// Return true if the a1=a2 bond is part of an amide.
int
IsAmide(const Molecule& m,
        const MoleculeData& molecule_data,
        atom_number_t a1,
        atom_number_t a2) {
  const Atom* atom1 = m.atomi(a1);
  const Atom* atom2 = m.atomi(a2);

  // Make sure atom1 is the Carbon or Sulphur.
  if (atom1->atomic_number() == 8 &&
     (atom2->atomic_number() == 6 || atom2->atomic_number() == 16)) {
    std::swap(atom1, atom2);
    std::swap(a1, a2);
  } else if (atom2->atomic_number() == 8 &&
    (atom1->atomic_number() == 6 || atom1->atomic_number() == 16)) {
  } else {
    return 0;
  }

  return molecule_data.alpha(a1, 7);

#ifdef SLOW
  // atom1 must be singly bonded to a Nitrogen.
  for (const Bond* b : *atom1) {
    if (b->is_double_bond()) {
      continue;
    }
    if (m.atomic_number(b->other(a1)) == 7) {
      return 1;
    }
  }

  return 0;
#endif
}

// There is a double or triple bond between `a1` and `a2` in `m`. Is it
// OK to lower that bond order?
int
Options::OkLowerBondOrder(Molecule& m,
                  const MoleculeData& molecule_data,
                  const Bond* b) const {
  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();

  if (! _config.destroy_amide() && IsAmide(m, molecule_data, a1, a2)) {
    return 0;
  }

  // We want to make sure we do not lower a bond order in a guanidine or amidine
  // or acid or nitro.
  const Atom* atom1 = m.atomi(a1);
  const Atom* atom2 = m.atomi(a2);

  // Both carbon, always lower.
  if (atom1->atomic_number() == 6 && atom2->atomic_number() == 6) {
    return 1;
  }

  // Both heteroatoms, very likely NO2 or S=O type. Do not destroy.
  if (atom1->atomic_number() != 6 && atom2->atomic_number() != 6) {
    return 0;
  }

  // Make sure `a1` is the carbon.
  if (atom1->atomic_number() == 6) {
  } else {
    std::swap(a1, a2);
    std::swap(atom1, atom2);
  }

  if (molecule_data.attached_heteroatom_count(a1) > 1) {
    return 0;
  }

  return 1;
}

int
Options::DoubleToSingleBond(Molecule& m,
                            const MoleculeData& molecule_data,
                            resizable_array_p<Molecule>& results) {
  if (! _config.double_to_single_bond()) {
    return 0;
  }

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (const Bond* b : m.bond_list()) {
    // We allow triple bonds to be lowered.
    if (b->is_single_bond() || b->is_aromatic()) {
      continue;
    }

    // always OK to lower a triple bond to a double bond.
    if (b->is_triple_bond()) {
    } else if (! OkLowerBondOrder(m, molecule_data, b)) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    if (molecule_data.five_valent_nitrogen(a1) ||
        molecule_data.five_valent_nitrogen(a2)) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    if (b->is_triple_bond()) {
      mcopy->set_bond_type_between_atoms(a1, a2, DOUBLE_BOND);
    } else {
      mcopy->set_bond_type_between_atoms(a1, a2, SINGLE_BOND);
    }
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kDoubleToSingleBond].extra(rc);

  return rc;
}

#ifdef NOT_USED_ASDASD
// Return true if atom `zatom` in `m` is bonded
// to a heteroatom.
bool
BondedToHeteroatom(const Molecule& m,
                   atom_number_t zatom) {
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    if (m.atomic_number(j) != 6) {
      return 1;
    }
  }

  return 0;
}
#endif

int
Options::ChangeCarbonToNitrogen(Molecule& m,
                                const MoleculeData& molecule_data,
                                resizable_array_p<Molecule>& results) {
  if (! _config.change_carbon_to_nitrogen()) {
    return 0;
  }

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    if (a.atomic_number() != 6) {
      continue;
    }

    if (! molecule_data.CanChange(i)) {
      continue;
    }

    if (a.ncon() == 4 || a.nbonds() == 4) {
      continue;
    }

    if (molecule_data.attached_heteroatom_count(i)) {
      continue;
    }

    if (molecule_data.beta(i, 7) || molecule_data.beta(i, 8)) {
      continue;
    }

    if (m.in_ring_of_given_size(i, 3)) {
      return 0;
    }

    if (m.is_aromatic(i) && a.ncon() != 2) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    mcopy->set_atomic_number(i, 7);
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kChangeCarbonToNitrogen].extra(rc);

  return rc;
}

int
Options::ChangeCarbonToOxygen(Molecule& m,
                              MoleculeData& molecule_data,
                              resizable_array_p<Molecule>& results) {
  if (! _config.change_carbon_to_oxygen()) {
    return 0;
  }

  const int matoms = m.natoms();
  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    if (a.atomic_number() != 6) {
      continue;
    }

    if (! molecule_data.CanChange(i)) {
      continue;
    }

    if (molecule_data.unsaturated(i)) {
      continue;
    }

    if (a.ncon() > 2) {
      continue;
    }

    if (m.is_aromatic(i)) {
      continue;
    }

    if (molecule_data.attached_heteroatom_count(i)) {
      continue;
    }

    if (molecule_data.beta(i, 7) || molecule_data.beta(i, 8)) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    mcopy->set_atomic_number(i, 8);
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kChangeCarbonToOxygen].extra(rc);

  return rc;
}

// Return true if `zatom` is singly connected and adjacent
// to an aromatic ring
int
ExocyclicToAromatic(Molecule& m,
                    atom_number_t zatom) {
  const Atom& a = m.atom(zatom);
  if (a.ncon() != 1) {
    return 0;
  }

  const atom_number_t j = a.other(zatom, 0);

  return m.is_aromatic(j);
}

int
Options::ChangeNitrogenToCarbon(Molecule& m,
                              MoleculeData& molecule_data,
                              resizable_array_p<Molecule>& results) {
  if (! _config.change_nitrogen_to_carbon()) {
    return 0;
  }

  const int matoms = m.natoms();
  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    if (a.atomic_number() != 7) {
      continue;
    }

    if (! molecule_data.CanChange(i)) {
      continue;
    }

    // We could possibly check for amides...

    // Avoid aromatics.
    if (m.is_aromatic(i) && a.ncon() == 3) {
      continue;
    }
    if (m.is_aromatic(i) && a.ncon() == 2 && m.hcount(i)) {
      continue;
    }
    
    // Do not destroy aromatic rings.
    if (ExocyclicToAromatic(m, i)) {
      continue;
    }

    if (molecule_data.five_valent_nitrogen(i)) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    mcopy->set_atomic_number(i, 6);
    if (a.formal_charge()) {
      mcopy->set_formal_charge(i, 0);
    }
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kChangeNitrogenToCarbon].extra(rc);

  return rc;
}

int
Options::RemoveCh2(Molecule& m,
                   const MoleculeData& molecule_data,
                   resizable_array_p<Molecule>& results) {
  if (! _config.remove_ch2()) {
    return 0;
  }

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    if (6 != a.atomic_number()) {
      continue;
    }
    if (a.ncon() != 2) {
      continue;
    }
    if (! molecule_data.CanChange(i)) {
      continue;
    }
    if (molecule_data.unsaturated(i)) {
      continue;
    }
    if (molecule_data.attached_heteroatom_count(i) == 2) {
      continue;
    }
    if (m.is_ring_atom(i)) {
      continue;
    }

    Set_of_Atoms connected;
    a.connections(i, connected);
    assert(connected.size() == 2);

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    mcopy->remove_atom(i);
    if (connected[0] > i) {
      --connected[0];
    }
    if (connected[1] > i) {
      --connected[1];
    }
    mcopy->add_bond(connected[0], connected[1], SINGLE_BOND);
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kRemoveCH2].extra(rc);

  return rc;
}

int
Options::InsertCh2(Molecule& m,
                   MoleculeData& molecule_data,
                   resizable_array_p<Molecule>& results) {
  if (! _config.insert_ch2()) {
    return 0;
  }

  int rc = 0;

  static const Element* carbon = get_element_from_atomic_number(6);

  for (const Bond* b : m.bond_list()) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    if (b->nrings() && m.is_aromatic(a1) && m.is_aromatic(a2)) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    const int initial_matoms = mcopy->natoms();
    mcopy->remove_bond_between_atoms(a1, a2);
    mcopy->add(carbon);
    mcopy->add_bond(a1, initial_matoms, SINGLE_BOND);
    mcopy->add_bond(a2, initial_matoms, SINGLE_BOND);
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kInsertCH2].extra(rc);

  return rc;
}

// For each bond where both atoms are `in_system`, set that
// bond to single.
int
AllInSystemBondsSingle(Molecule& m,
                       const MoleculeData& molecule_data,
                       const int * in_system) {
  // cerr << "Checking " << m.smiles() << '\n';
  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    if (! b->is_double_bond()) {
      continue;
    }

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    if (! in_system[a1] || ! in_system[a2]) {
      continue;
    }

    if (molecule_data.five_valent_nitrogen(a1) ||
        molecule_data.five_valent_nitrogen(a2)) {
      return 0;
    }

    // Try to deal with fused rings - this is not correct, but
    // close enough. Will fail if rings fused at each end of the bond.
    if (m.ring_bond_count(a1) == 3 && m.ring_bond_count(a2) == 3) {
      continue;
    }

    m.set_bond_type_between_atoms(b->a1(), b->a2(), SINGLE_BOND);
    ++rc;
  }

  // cerr << "With single bonds " << m.smiles() << '\n';

  // Neutralize formal charges.
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (! in_system[i]) {
      continue;
    }

    const Atom& atom = m.atom(i);
    if (atom.formal_charge() == 0) {
      continue;
    }

    // cerr << "Atom " << i << " has formal " << atom.atomic_number() << " ncon " << atom.ncon() << " nbonds " << m.nbonds(i) << '\n';
    // C1=C[N+]2=C(N=C1)N(CCCC)[C@](O)(C1=CC=CC=C1)C2 CHEMBL1851918
    if (atom.atomic_number() == 7 && atom.ncon() == 3 && m.nbonds(i) == 4) {
      continue;
    }

    if (atom.atomic_number() == 7 && atom.ncon() == 4) {
      continue;
    }

    m.set_formal_charge(i, 0);
  }

  // Fail if any n=O groups found.

  for (int i = 0; i < matoms; ++i) {
    if (! in_system[i]) {
      continue;
    }

    const Atom& atom = m.atom(i);
    // cerr << " z " << atom.atomic_number() << " ncon " << atom.ncon() << '\n';
    if (atom.atomic_number() != 7) {
      continue;
    }
    if (atom.ncon() != 3) {
      continue;
    }
    // cerr << "Atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " maybe\n";
    for (const Bond* b : atom) {
      if (! b->is_double_bond()) {
        continue;
      }

      atom_number_t j = b->other(i);
      if (m.atomic_number(j) == 8) {
        return 0;
      }
    }
  }

  return rc;
}


int
Options::DestroyAromaticRings(Molecule& m,
                              MoleculeData& molecule_data,
                              resizable_array_p<Molecule>& results) {
  if (! _config.destroy_aromatic_rings()) {
    return 0;
  }

  int rc = 0;

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  std::unique_ptr<int[]> in_system = std::make_unique<int[]>(matoms);

  for (const Ring* r : m.sssr_rings()) {
    // cerr << "CHecking ring " << *r << '\n';
    if (! r->is_aromatic()) {
      continue;
    }

    if (! molecule_data.AllAtomsCanChange(*r)) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    std::fill_n(in_system.get(), matoms, 0);
    r->set_vector(in_system.get(), 1);
    if (! AllInSystemBondsSingle(*mcopy, molecule_data, in_system.get())) {
      continue;
    }
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kDestroyAromaticRings].extra(rc);

  return rc;
}

// In each ring system, destroy all aromatic rings.
int
Options::DestroyAromaticRingSystems(Molecule& m,
                              MoleculeData& molecule_data,
                              resizable_array_p<Molecule>& results) {
  if (! _config.destroy_aromatic_ring_systems()) {
    return 0;
  }

  const int nr = m.nrings();
  if (nr == 0) {
    return 0;
  }

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  std::unique_ptr<int[]> ring_already_done(new_int(nr));
  std::unique_ptr<int[]> in_system = std::make_unique<int[]>(matoms);

  int rc = 0;

  for (int i = 0; i < nr; ++i) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);
    if (! ri->is_aromatic()) {
      continue;
    }

    ring_already_done[i] = 1;

    std::fill_n(in_system.get(), matoms, 0);
    ri->set_vector(in_system.get(), 1);

    for (int j = 0; j < nr; ++j) {
      if (ring_already_done[j]) {
        continue;
      }

      const Ring* rj = m.ringi(j);
      if (ri->fused_system_identifier() != rj->fused_system_identifier()) {
        continue;
      }
      ring_already_done[j] = 1;
      rj->set_vector(in_system.get(), 1);
    }

    if (! molecule_data.AllAtomsCanChange(in_system.get())) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    if (! AllInSystemBondsSingle(*mcopy, molecule_data, in_system.get())) {
      continue;
    }
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kDestroyAromaticRingSystems].extra(rc);

  return rc;
}

// If there is an atom that is bonded to both `a1` and `a2`
// return it. Otherwise INVALID_ATOM_NUMBER.
atom_number_t
AtomBetween(const Molecule& m,
        atom_number_t a1,
        atom_number_t a2) {
  const Atom& atom2 = m.atom(a2);

  atom_number_t result = INVALID_ATOM_NUMBER;

  for (const Bond* b : m.atom(a1)) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t j = b->other(a1);
    // Fail if the atoms are bonded.
    if (j == a2) {
      return INVALID_ATOM_NUMBER;
    }

    if (atom2.is_bonded_to(j)) {
      result = j;
    }
  }

  return result;
}

// Make three membered rings.
int
Options::MakeThreeRing(Molecule& m,
                        MoleculeData& molecule_data,
                        resizable_array_p<Molecule>& results) {
  if (! _config.make_three_membered_rings()) {
    return 0;
  }

  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (m.atomic_number(i) != 6 || molecule_data.unsaturated(i)) {
      continue;
    }
    if (m.hcount(i) == 0) {
      continue;
    }
    if (! molecule_data.CanChange(i)) {
      continue;
    }

    for (int j = i + 1; j < matoms; ++j) {
      if (m.atomic_number(j) != 6 || molecule_data.unsaturated(j)) {
        continue;
      }
      if (! molecule_data.CanChange(j)) {
        continue;
      }

      if (m.hcount(j) == 0) {
        continue;
      }
      const atom_number_t k = AtomBetween(m, i, j);
      // cerr << "Atoms " << i << " and " << j << " btw " << k << '\n';
      if (k == INVALID_ATOM_NUMBER || m.atomic_number(k) != 6) {
        continue;
      }
      if (! molecule_data.CanChange(k)) {
        continue;
      }

      if (m.are_bonded(i, j)) {
        cerr << "atoms " << i << " " << m.smarts_equivalent_for_atom(i) << " and " << j << ' ' << m.smarts_equivalent_for_atom(j) << " already bonded\n";
        cerr << " btw atom " << k << '\n';
        write_isotopically_labelled_smiles(m, false, cerr);
        cerr << '\n';
        continue;
      }
      std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
      mcopy->add_bond(i, j, SINGLE_BOND);
      rc += AddToResultsIfNew(mcopy, results);
    }
  }

  _acc[TransformationType::kMakeThreeMemberedRings].extra(rc);

  return rc;
}

// Make a new molecule by breaking the bond btw `a1` and `a2`
// and discarding the fragment containing `a2`.
int
Options::RemoveFragment(Molecule& m,
               atom_number_t a1,
               atom_number_t a2,
               resizable_array_p<Molecule>& results) {
  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
  mcopy->remove_bond_between_atoms(a1, a2);
  mcopy->remove_fragment_containing_atom(a2);
  return AddToResultsIfNew(mcopy, results);
}

int
Options::RemoveFragments(Molecule& m,
                         MoleculeData& molecule_data,
                         resizable_array_p<Molecule>& results) {
  if (_config.remove_fragment().empty()) {
    return 0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<int[]> frag = std::make_unique<int[]>(matoms);

  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    // Note that we do not check the bond type. Maybe we should: sulfonamide?
    if (b->nrings()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    // If either of these is terminal and we are removing terminal
    // groups, do it.
    if (molecule_data.remove_fragment(1)) {
      if (m.ncon(a1) == 1) {
        rc += RemoveFragment(m, a2, a1, results);
      } else if (m.ncon(a2) == 1) {
        rc += RemoveFragment(m, a1, a2, results);
      }
      continue;
    }

    std::fill_n(frag.get(), matoms, 0);
    m.identify_side_of_bond(frag.get(), a1, 1, a2);
    int frag_size = std::count(frag.get(), frag.get() + matoms, 1);
    if (frag_size > matoms / 2) {
      std::swap(a1, a2);
      frag_size = matoms - frag_size;
    }
    if (! molecule_data.remove_fragment(frag_size)) {
      continue;
    }

    if (! molecule_data.AllAtomsCanChange(frag.get())) {
      continue;
    }

    std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
    mcopy->remove_bond_between_atoms(a1, a2);
    mcopy->remove_fragment_containing_atom(a2);
    rc += AddToResultsIfNew(mcopy, results);
  }

  _acc[TransformationType::kRemoveFragment].extra(rc);

  return rc;
}

// copied from random_molecular_permuations.
// Common code for swapping adjacent atoms. It can be called by the do and undo transformation
//     A1-J-A2-A3
//  becomes
//     A1-A2-J-A3

void
DoSwapAdjacentAtoms(Molecule & m,
                    const atom_number_t a1,
                    const atom_number_t j,
                    const atom_number_t a2,
                    const atom_number_t a3) {
  if (a1 != INVALID_ATOM_NUMBER) {
    m.remove_bond_between_atoms(a1, j);
  }

  m.remove_bond_between_atoms(j, a2);

  if (a3 != INVALID_ATOM_NUMBER) {
    m.remove_bond_between_atoms(a2, a3);
  }

  if (a1 != INVALID_ATOM_NUMBER) {
    m.add_bond(a1, a2, SINGLE_BOND);
  }

  m.add_bond(a2, j, SINGLE_BOND);

  if (a3 != INVALID_ATOM_NUMBER) {
    m.add_bond(j, a3, SINGLE_BOND);
  }

  return;
}

// We are contemplating joining atoms A1 and A2

int
AnyUndesirableAdjacencies(Molecule & m,
                          atom_number_t a1,
                          atom_number_t a2) {
  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  atomic_number_t z1 = aa1->atomic_number();
  atomic_number_t z2 = aa2->atomic_number();

  if (6 != z1 && 6 != z2)
    return 1;    // never join heteroatoms

// Avoid forming alkyl halides

  if (aa1->is_halogen()) {
    if (! m.is_aromatic(a2)) {
      return 1;
    }
  }
  else if (aa2->is_halogen()) {
    if (! m.is_halogen(a1)) {
      return 1;
    }
  }

// Never make adjacent carbonyls or such...

  if (m.multiple_bond_to_heteroatom(a1) && m.multiple_bond_to_heteroatom(a2))
    return 1;

  return 0;
}


// We are just about to swap atoms A1-J-A2-A3 to A1-A2-J-A3
// Return true if that is problematic.

int
AnyUndesirableAdjacencies(Molecule & m,
                          const atom_number_t a1,
                          const atom_number_t j,
                          const atom_number_t a2,
                          const atom_number_t a3)
{
  if (a1 == INVALID_ATOM_NUMBER) {
  } else if (AnyUndesirableAdjacencies(m, a1, a2)) {
    return 1;
  }

  if (a3 == INVALID_ATOM_NUMBER) {
  } else if (AnyUndesirableAdjacencies(m, j, a3)) {
    return 1;
  }

  return 0;
}

// Within `m` the atoms look like a0-a1-a2-a3 with any kind of bond.
int
Options::SwapAdjacentAtoms(Molecule& m,
                           MoleculeData& molecule_data,
                           atom_number_t a0,
                           atom_number_t a1,
                           atom_number_t a2,
                           atom_number_t a3,
                           resizable_array_p<Molecule>& results) {
  // Avoid three membered rings.
  if (a0 == INVALID_ATOM_NUMBER) {
  } else if (m.are_bonded(a0, a2)) {
    return 0;
  }

  if (a3 == INVALID_ATOM_NUMBER) {
  } else if (m.are_bonded(a1, a3)) {
    return 0;
  }

  if (AnyUndesirableAdjacencies(m, a0, a1, a2, a3)) {
    return 0;
  }

  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);

  DoSwapAdjacentAtoms(*mcopy, a0, a1, a2, a3);

  return AddToResultsIfNew(mcopy, results);
}

// Get an atom that is connected to `zatom` that is not `exclude`.
atom_number_t
GetConnectedAtom(const Atom& atom,
                 atom_number_t zatom,
                 atom_number_t exclude) {
  for (const Bond* b : atom) {
    atom_number_t j = b->other(zatom);
    if (j == exclude) {
      continue;
    }

    return j;
  }

  return INVALID_ATOM_NUMBER;
}

// Adapted from random_molecular_permutations.
int
Options::SwapAdjacentAtoms(Molecule& m,
                           MoleculeData& molecule_data,
                           resizable_array_p<Molecule>& results) {
  if (! _config.swap_adjacent_atoms()) {
    return 0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<int[]> aromatic;
  if (_config.swap_adjacent_aromatic_atoms()) {
    aromatic.reset(new_int(matoms));
    m.aromaticity(aromatic.get());
  }

  int rc = 0;

  for (const Bond* b : m.bond_list()) {
    if (b->is_aromatic()) {
      if (_config.swap_adjacent_aromatic_atoms()) {
      } else {
        continue;
      }
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    if (_config.swap_adjacent_aromatic_atoms()) {
    } else if (molecule_data.unsaturated(a1) || molecule_data.unsaturated(a2)) {
      continue;
    }

    const Atom& at1 = m.atom(a1);
    const Atom& at2 = m.atom(a2);

    // Skip if nothing would change.
    if ((at1.atomic_number() == at2.atomic_number()) &&
         (at1.ncon() == at2.ncon()) &&
         (molecule_data.unsaturated(a1) == molecule_data.unsaturated(a2))) {
      continue;
    }
      
    // Get the connected atoms outside the bond. We only
    // process two connected atoms for now.

    const atom_number_t a0 = GetConnectedAtom(at1, a1, a2);
    const atom_number_t a3 = GetConnectedAtom(at2, a2, a1);
    // cerr << "Atoms are " << a0 << ' ' << a1 << ' ' << a2 << ' ' << a3 << '\n';

    rc += SwapAdjacentAtoms(m, molecule_data, a0, a1, a2, a3, results);
  }

  _acc[TransformationType::kSwapAdjacentAtoms].extra(rc);

  return rc;
}

// `atom` which is atom `zatom` is a member of a ring. Ring members
// are set in `in_ring`.
// Identify the two atoms attached to `zatom` which are also in the ring.
int
GetRingMembers(const Atom& atom, atom_number_t zatom, const int* in_ring,
               atom_number_t& r1a, atom_number_t& r1b) {
  assert(in_ring[zatom]);

  r1a = INVALID_ATOM_NUMBER;
  r1b = INVALID_ATOM_NUMBER;

  for (const Bond* b : atom) {
    atom_number_t j = b->other(zatom);
    if (! in_ring[j]) {
      continue;
    }
    if (r1a == INVALID_ATOM_NUMBER) {
      r1a = j;
    } else if (r1b == INVALID_ATOM_NUMBER) {
      r1b = j;
    } else {  // not sure this can happen.
      return 0;
    }
  }

  return r1b != INVALID_ATOM_NUMBER;
}

// If there is just a single atom in `ring` that is marked in `in_ring`,
// return that atom. If not, return INVALID_ATOM_NUMBER.
atom_number_t 
SingleAtomShared(const Ring& ring, const int* in_ring) {
  atom_number_t result = INVALID_ATOM_NUMBER;

  const int ring_size = ring.number_elements();
  for (int i = 0; i < ring_size; ++i) {
    atom_number_t a = ring[i];

    if (in_ring[a] == 0) {
      continue;
    }

    // If we have already found a shared atom, the rings are fused. Fail.
    if (result != INVALID_ATOM_NUMBER) {
      return INVALID_ATOM_NUMBER;
    }

    result = a;
  }

  return result;
}

// Change a spiro ring fusion to 
int
Options::Unspiro(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {
  const int nrings = m.nrings();
  if (nrings < 2) {
    return 0;
  }

  if (! _config.unspiro()) {
    return 0;
  }

  static const Element* carbon = get_element_from_atomic_number(6);

  const int matoms = m.natoms();

  std::unique_ptr<int[]> in_ring = std::make_unique<int[]>(matoms);

  int rc = 0;

  for (int i = 0; i < nrings; ++i) {
    const Ring* ri = m.ringi(i);
    std::fill_n(in_ring.get(), matoms, 0);
    ri->set_vector(in_ring.get(), 1);

    for (int j = i + 1; j < nrings; ++j) {
      const Ring* rj = m.ringi(j);
      atom_number_t spiro = SingleAtomShared(*rj, in_ring.get());
      if (spiro == INVALID_ATOM_NUMBER) {
        continue;
      }

      const Atom* a = m.atomi(spiro);
      if (a->ncon() != 4) {  // Not sure this can happen.
        continue;
      }

      if (! molecule_data.CanChange(spiro)) {
        continue;
      }

      atom_number_t r1a, r1b;
      if (! GetRingMembers(*a, spiro, in_ring.get(), r1a, r1b)) {
        continue;
      }

      // cerr << "Spiro " << spiro <<  " and " << r1a << ',' << r1b << '\n';

      std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
      mcopy->remove_bond_between_atoms(r1a, spiro);
      mcopy->remove_bond_between_atoms(r1b, spiro);
      mcopy->add(carbon);
      mcopy->add_bond(matoms, r1a, SINGLE_BOND);
      mcopy->add_bond(matoms, r1b, SINGLE_BOND);
      mcopy->add_bond(matoms, spiro, SINGLE_BOND);

      rc += AddToResultsIfNew(mcopy, results);
    }
  }

  _acc[TransformationType::kUnspiro].extra(rc);

  return rc;
}
  
int
AtomTypeMatch(Molecule& m,
              const MoleculeData& molecule_data,
              atom_number_t zatom,
              const Fragment& frag) {
  const uint32_t iso = frag.mol().isotope(frag.attachment_point());
  // cerr << "AtomTypeMatch iso " << iso << " cmp " << molecule_data.atype(zatom) << '\n';
  return iso == molecule_data.atype(zatom);
}

int
Options::AddFragments(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {
  if (! _config.add_fragments()) {
    return 0;
  }

  if (_fragment.empty()) {
    return 0;
  }

  int rc = 0;

  const int matoms = m.natoms();

  for (const Fragment* frag : _fragment) {
    for (int i = 0; i < matoms; ++i) {
      if (molecule_data.hcount(i) == 0) {
        continue;
      }
      if (! molecule_data.CanChange(i)) {
        continue;
      }

      if (! _atom_typing.active()) {
      } else if (! AtomTypeMatch(m, molecule_data, i, *frag)) {
        continue;
      }
      if (m.atomic_number(i) != 6 &&
          frag->mol().atomic_number(frag->attachment_point()) != 6) {
        continue;
      }

      std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
      mcopy->add_molecule(&frag->mol());
      mcopy->add_bond(i, matoms + frag->attachment_point(), SINGLE_BOND);
      mcopy->unset_all_implicit_hydrogen_information(matoms + frag->attachment_point());
      rc += AddToResultsIfNew(mcopy, results);
    }
  }

  _acc[TransformationType::kAddFragments].extra(rc);

  return rc;
}

// We are considering joining `frag` to atoms `a1` and `a2` in `m`.
// Does this result in any adjacent heteroatoms?
int
OkAdjacencies(Molecule& m,
                const MoleculeData& molecule_data,
                atom_number_t a1,
                atom_number_t a2,
                const BivalentFragment& frag) {
  const atomic_number_t z1 = m.atomic_number(a1);
  const atomic_number_t z2 = m.atomic_number(a2);

  const atomic_number_t f1 = frag.mol().atomic_number(frag.attach1());
  const atomic_number_t f2 = frag.mol().atomic_number(frag.attach2());

  // cerr << "OkAdjacencies " << z1 << '-' << f1 << ' ' << z2 << '-' << f2 << '\n';

  // If everything is carbon we are good.
  if (z1 == 6 && f1 == 6 && z2 == 6 && f2 == 6) {
    // cerr << "All carbon, adjacencies ok\n";
    return 1;
  }

  // Do not join heteroatoms.
  if (z1 != 6 || f1 != 6) {
    return 0;
  }
  if (z2 != 6 || f2 != 6) {
    return 0;
  }

  if (f1 != 6 && molecule_data.attached_heteroatom_count(a1)) {
    return 0;
  }
  if (f2 != 6 && molecule_data.attached_heteroatom_count(a2)) {
    return 0;
  }

  // cerr << "Adjacencies OK\n";
  return 1;
}

int
Options::InsertBivalentFragments(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {
  if (! _config.insert_fragments()) {
    return 0;
  }

  if (_bivalent_fragment.empty()) {
    return 0;
  }

  int rc = 0;

  // cerr << "InsertBivalentFragments:have " << _bivalent_fragment.size() << " fragments\n";
  // Not const because the hcount() method is called.
  for (BivalentFragment* frag : _bivalent_fragment) {
    for (const Bond* b : m.bond_list()) {
      atom_number_t a1 = b->a1();
      atom_number_t a2 = b->a2();
      if (! molecule_data.CanChange(a1, a2)) {
        continue;
      }

      if (OkAdjacencies(m, molecule_data, a1, a2, *frag)) {
        rc += InsertBivalentFragment(m, molecule_data, a1, a2, *frag, results);
      }

      if (OkAdjacencies(m, molecule_data, a2, a1, *frag)) {
        rc += InsertBivalentFragment(m, molecule_data, a2, a1, *frag, results);
      }
    }
  }

  _acc[TransformationType::kInsertFragments].extra(rc);

  return rc;
}

// insert `frag` between atoms `a1` and `a2` in `m`.
int
Options::InsertBivalentFragment(const Molecule& m,
                 const MoleculeData& molecule_data,
                 atom_number_t a1,
                 atom_number_t a2,
                 BivalentFragment& frag,
                 resizable_array_p<Molecule>& results) {
  const int matoms = m.natoms();

  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
  // cerr << "InsertBivalentFragment " << frag.mol().smiles() << " into " << mcopy->smiles() << '\n';

  mcopy->add_molecule(&frag.mol());
  mcopy->remove_bond_between_atoms(a1, a2);

  const atom_number_t f1 = matoms + frag.attach1();
  const atom_number_t f2 = matoms + frag.attach2();

  // If atom typing active, make sure things match.
  if (_atom_typing.active()) {
    uint32_t frag_atype1 = frag.atype1();
    uint32_t frag_atype2 = frag.atype2();
    // cerr << "Fragment atom types " << frag_atype1 << "," << frag_atype2 << " scaffold " << molecule_data.atype(a1) << " and " << molecule_data.atype(a2) << '\n';
    if (frag_atype1 != molecule_data.atype(a1)) {
      return 0;
    }
    if (frag_atype2 != molecule_data.atype(a2)) {
      return 0;
    }
  }

  if (f1 == f2 && frag.mol().hcount(frag.attach1()) < 2) {
    // cerr << "Insufficient hcount " << frag.mol().hcount(frag.attach1()) << '\n';
    return 0;
  }

  mcopy->add_bond(a1, f1, SINGLE_BOND);
  mcopy->add_bond(a2, f2, SINGLE_BOND);
  mcopy->unset_all_implicit_hydrogen_information(f1);
  mcopy->unset_all_implicit_hydrogen_information(f2);

  return AddToResultsIfNew(mcopy, results);
}

// We are considering inserting `frag` between atoms `a1` and `a2` in `m`.
// Return OK if the atom types match.
int
Options::OkAtomTypes(Molecule& m,
            const MoleculeData& molecule_data,
            atom_number_t a1,
            atom_number_t a2,
            const BivalentFragment& frag) const {
  if (! _atom_typing.active()) {
    return 1;
  }

  uint32_t atype1 = molecule_data.atype(a1);
  uint32_t atype2 = molecule_data.atype(a2);

  if (atype1 != frag.mol().isotope(frag.attach1())) {
    // cerr << "Fail1\n";
    return 0;
  }
  if (atype2 != frag.mol().isotope(frag.attach2())) {
    // cerr << "Fail2\n";
    return 0;
  }

  // cerr << "OkAtomTypes alright\n";

  return 1;
}

// Identify a path from `zatom` to `destination` in `m`,
// marking the atoms visited in `destination`.
int
IdentifyPath(Molecule& m,
             atom_number_t zatom, 
             atom_number_t destination, 
             int * shortest_path) {
  const int current_distance = m.bonds_between(zatom, destination);

  const Atom& a = m.atom(zatom);
  for (const Bond*b : a) {
    atom_number_t j = b->other(zatom);
    if (j == destination) {
      return 1;
    }

    if (m.bonds_between(j, destination) >= current_distance) {
      continue;
    }
    shortest_path[j] = 1;
    IdentifyPath(m, j, destination, shortest_path);
  }

  return 1;
}

// When we find two atoms that are separated by the same bond separation
// as a fragment, excise the intervening atoms and inser the fragment.
int
Options::ReplaceInnerFragments(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {
  if (! _config.replace_inner_fragments()) {
    return 0;
  }

  if (_bivalent_fragment.empty()) {
    return 0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<int[]> shortest_path = std::make_unique<int[]>(matoms);

  int rc = 0;

  for (const BivalentFragment* frag : _bivalent_fragment) {
    for (int i = 0; i < matoms; ++i) {
      if (! molecule_data.CanChange(i)) {
        continue;
      }
      for (int j = i + 1; j < matoms; ++j) {
        if (! molecule_data.CanChange(j)) {
          continue;
        }

        if (m.bonds_between(i, j) != frag->bonds_between() + 2) {
          continue;
        }

        std::fill_n(shortest_path.get(), matoms, 0);
        if (! IdentifyPath(m, i, j, shortest_path.get())) {
          continue;
        }
        if (OkAtomTypes(m, molecule_data, i, j, *frag)) {
          rc += ReplaceInnerFragment(m, i, j, shortest_path.get(), *frag, results);
        }
        if (OkAtomTypes(m, molecule_data, j, i, *frag)) {
          rc += ReplaceInnerFragment(m, j, i, shortest_path.get(), *frag, results);
        }
      }
    }

    if (GeneratedEnough(results)) {
      break;
    }
  }

  _acc[TransformationType::kReplaceInnerFragments].extra(rc);

  return rc;
}

// Remove any of the first `matoms` in `m` that are set in `remove`.
int
RemoveAtoms(Molecule& m,
            int matoms,
            const int* remove) {
  int rc = 0;
  for (int i = matoms - 1; i >= 0; --i) {
    if (remove[i]) {
      m.remove_atom(i);
      ++rc;
    }
  }

  return rc;
}

int
Options::ReplaceInnerFragment(const Molecule& m,
                              atom_number_t a1,
                              atom_number_t a2,
                              int* shortest_path,
                              const BivalentFragment& frag,
                              resizable_array_p<Molecule>& results) {
  const int matoms = m.natoms();

  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);

  mcopy->add_molecule(&frag.mol());
  // cerr << "After adding frag " << mcopy->smiles() << '\n';

  atom_number_t f1 = matoms + frag.attach1();
  atom_number_t f2 = matoms + frag.attach2();

  mcopy->add_bond(a1, f1, SINGLE_BOND);
  mcopy->add_bond(a2, f2, SINGLE_BOND);
  mcopy->unset_all_implicit_hydrogen_information(f1);
  mcopy->unset_all_implicit_hydrogen_information(f2);
  RemoveAtoms(*mcopy, matoms, shortest_path);

  return AddToResultsIfNew(mcopy, results);
}

int
Count(const int* haystack,
      int n,
      int needle) {
  return std::count(haystack, haystack + n, needle);
}

// Invert the `n` values in `values`.
void
Invert(int* values,
       uint32_t n) {
  std::transform(values, values + n, values, [](int i) {
    return !i;
  });
}

// Return true if the atomic number of `zatom` in `m` is the same
// as the attachment atom in `frag`.
bool
AtomicNumbersMatch(const Molecule& m,
                   atom_number_t zatom,
                   const Fragment& frag) {
  return m.atomic_number(zatom) == 
         frag.mol().atomic_number(frag.attachment_point());
}

int
Options::ReplaceTerminalFragments(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {
  if (!_config.replace_terminal_fragments()) {
    return 0;
  }

  if (_fragment.empty()) {
    return 0;
  }

  const uint32_t matoms = m.natoms();

  std::unique_ptr<int[]> either_side = std::make_unique<int[]>(matoms);

  int rc = 0;

  // Force ring membership.
  (void) m.sssr_rings();

  for (const Bond* b : m.bond_list()) {
    if (! b->is_single_bond() || b->nrings()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    std::fill_n(either_side.get(), matoms, 0);

#ifdef DEBUG_REPLACE_TERMINAL_GROUP
    cerr << " Atoms " << a1 << ' ' << m.smarts_equivalent_for_atom(a1) << " and " << a2 << ' ' << m.smarts_equivalent_for_atom(a2) << '\n';
    write_isotopically_labelled_smiles(m, false, cerr);
    cerr << '\n';
#endif
    m.identify_side_of_bond(either_side.get(), a1, 1, a2);
#ifdef DEBUG_REPLACE_TERMINAL_GROUP
    for (int i = 0; i < matoms; ++i) {
      cerr << " ES " << i << ' ' << m.smarts_equivalent_for_atom(i) << " side " << either_side[i] << '\n';
    }
#endif

    uint32_t frag_size = Count(either_side.get(), matoms, 1);
    if (frag_size > matoms / 2) {
      frag_size = matoms - frag_size;
      Invert(either_side.get(), matoms);
      std::swap(a1, a2);
    }
#ifdef DEBUG_REPLACE_TERMINAL_GROUP
    '' cerr << " a1 " << a2 << " a2 " << a2 << " frag_size " << frag_size << '\n';
    for (int i = 0; i < matoms; ++i) {
      cerr << i << " side " << either_side[i] << '\n';
    }
#endif

    if (_config.max_atoms_lost() > 0 && frag_size > _config.max_atoms_lost()) {
      continue;
    }

    for (Fragment* frag : _fragment) {
      // Only do the replacement if the atomic numbers match.
      if (! AtomicNumbersMatch(m, a2, *frag)) {
        continue;
      }

      std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
      mcopy->add_molecule(&frag->mol());
      mcopy->remove_bond_between_atoms(a1, a2);
      mcopy->add_bond(a2, matoms + frag->attachment_point(), SINGLE_BOND);
      mcopy->remove_fragment_containing_atom(a1);
      mcopy->unset_all_implicit_hydrogen_information(matoms + frag->attachment_point());
      // cerr << "Added " << frag->mol().smiles() << " to form "  << mcopy->smiles() << '\n';
      rc += AddToResultsIfNew(mcopy, results);
    }

    if (GeneratedEnough(results)) {
      break;
    }
  }

  _acc[TransformationType::kReplaceTerminalFragments].extra(rc);

  return rc;
}

int
Options::FuseBiphenyls(Molecule& m,
                 MoleculeData& molecule_data,
                 resizable_array_p<Molecule>& results) {
  // not implemented yet, TODO:ianwatson
  return 0;
  if (! _config.fuse_biphenyls()) {
    return 0;
  }

  m.compute_aromaticity_if_needed();

  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    if (b->is_aromatic()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (! molecule_data.CanChange(a1, a2)) {
      continue;
    }

    // Do not form large fused ring systems.
    if (molecule_data.fused_system_size(a1) != 1 ||
        molecule_data.fused_system_size(a2) != 1) {
      continue;
    }

    // Make things more likely to succeed. Also the most common case.
    if (m.atomic_number(a1) != 6 || m.atomic_number(a2) != 6) {
      continue;
    }

    if (! m.is_aromatic(a1) || ! m.is_aromatic(a2)) {
      continue;
    }

    rc += FuseBiphenyls(m, molecule_data, a1, a2, results);
  }

  _acc[TransformationType::kFuseBiphenyls].extra(rc);

  return rc;
}

// Identify any aromatic neighbours of `zatom` that have a hydrogen atom.
Set_of_Atoms
GetAromaticNbrs(Molecule& m,
                const MoleculeData& molecule_data,
                atom_number_t zatom) {
  Set_of_Atoms result;
  for (const Bond* b : m.atom(zatom)) {
    if (! b->is_aromatic()) {
      continue;
    }

    atom_number_t j = b->other(zatom);
    if (molecule_data.hcount(j) == 0) {
      continue;
    }

    result << j;
  }

  return result;
}

int
Options::FuseBiphenyls(Molecule& m,
                 MoleculeData& molecule_data,
                 atom_number_t a1,
                 atom_number_t a2,
                 resizable_array_p<Molecule>& results) {
  Set_of_Atoms nbrs1 = GetAromaticNbrs(m, molecule_data, a1);
  if (nbrs1.empty()) {
    return 0;
  }
  Set_of_Atoms nbrs2 = GetAromaticNbrs(m, molecule_data, a2);
  if (nbrs2.empty()) {
    return 0;
  }

  int rc = 0;

  for (atom_number_t n1 : nbrs1) {
    for (atom_number_t n2 : nbrs2) {
      rc += FuseBiphenyl(m, molecule_data, a1, n1, a2, n2, results);
      rc += FuseBiphenyl(m, molecule_data, a2, n2, a1, n1, results);
    }
  }

  return rc;
}

// Return a neighbour of 3 connected atom `atom` that is not
// either of the excluded nbrs.
// That nbr must have 2 connections
atom_number_t
GetRemainingNbr(const Molecule& m,
                atom_number_t zatom,
                atom_number_t exclude1,
                atom_number_t exclude2) {
  assert(m.ncon(zatom) == 3);

  for (const Bond*b : m.atom(zatom)) {
    atom_number_t j = b->other(zatom);
    if (j == exclude1 || j == exclude2) {
      continue;
    }

    return j;
  }

  return INVALID_ATOM_NUMBER;
}

/*
   Looks like:

    - n1         n2 -
        \       /
         a1 - a2
        /       \
    - x1         x2 -

   Fuse the rings by forming


*/
int
Options::FuseBiphenyl(Molecule& m,
                 MoleculeData& molecule_data,
                 atom_number_t a1,
                 atom_number_t n1,
                 atom_number_t a2,
                 atom_number_t n2,
                 resizable_array_p<Molecule>& results) {
  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
  const atomic_number_t z1 = m.atomic_number(a1);
  const atomic_number_t z2 = m.atomic_number(a2);

  // n1 and n2 must merge. If one is a hetoeratom, use that. n2 will
  // be removed;
  if (z1 == z2) {
    // same, does not matter.
  } else if (z1 == 6) {
    std::swap(a1, a2);
    std::swap(n1, n2);
  } else {
    // n1 is a heteroatom, n2 will be lost.
  }

  atom_number_t x1 = GetRemainingNbr(m, a1, a2, n1);
  if (x1 == INVALID_ATOM_NUMBER) {
    return 0;
  }
  atom_number_t x2 = GetRemainingNbr(m, a2, a1, n2);
  if (x2 == INVALID_ATOM_NUMBER) {
    return 0;
  }

  // TODO:ianwatson finish this...
  // get the bond types, and merge across a common bond type.
  // mcopy->add_bond(

  return 1;
}

}  // namespace minor_changes
