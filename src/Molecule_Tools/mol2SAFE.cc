// Convert molecules to SAFE form.
// Gotta be SAFE: A New Framework for Molecular Design
// Emanuel Noutahi, Christian Gabellini, Michael Craig
// Jonathan Lim, Prudencio Tosou. Valence Labs.
// https://arxiv.org/pdf/2310.10773.pdf

#include <algorithm>
#include <iostream>
#include <memory>

#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"

#include "fragment_molecule.h"
#include "highest_ring_number.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace mol2safe {

using lillymol::RingNumberControl;
using std::cerr;

void
Usage(int rc) {
  cerr << R"(Transform smiles to SAFE representations.
Input must be a smiles file.
 -p             write the parent molecule as well
 -c             remove chirality
 -l             strip to largest fragment
 -y             re-use ring numbers in generated ring numbers
 -I <iso>       place a fixed isotope on all attachment points
 -P <atype>     atom typing specification. If specified the isotope will be the atom type.
 -S <fname>     write fragment statistics in dicer_data::DicerFragment textproto form
 -M ...         constraints on max fragment size
 -M <n>         fragments can contain no more than <n> atoms
 -M maxnr=<n>   the maximum number of non-ring atoms that can be in a fragment
 -g ...         standardisation options
 -z             ignore connection table errors on input
 -T ...         Element transformations - suggest '-T Br=Cl -T I=Cl'
 -t             test the unique smiles of the parent against that formed from the SAFE representation.
                This only works if chirality is removed and isotopes NOT added.
 -v             verbose output
  )";

  ::exit(rc);
}

// In order to avoid passing around many arguments.
struct PerMoleculeArrays {
 public:
  int atoms_in_parent;

  // For each atom across a broken bond, the aton(s) on the other side.
  Set_of_Atoms* nbrs;

  // Each component needs an array of strings for the atomic smarts.
  IWString* atom_smarts;

  // an natoms*natoms array that holds what is going on between two
  // atom numbers (indexed into the starting molecule). Initialised
  // to zero, which means nothing set. If the two atoms end up in
  // different fragments, this entry will be set.
  int* ring_number_status;

  // when forming the per fragment smiles  we need an include_atom array.
  // If we size this to the size of the starting molecule, it will work for
  // every fragment
  int* include_atom;

  uint32_t* atype = nullptr;

  IWString parent_smiles;
  IWString parent_name;

 public:
  PerMoleculeArrays(Molecule& m);
  ~PerMoleculeArrays();

  int
  RingStatusIndex(atom_number_t a1, atom_number_t a2) const {
    return a1 * atoms_in_parent + a2;
  }

  void
  SetRingNumber(atom_number_t a1, atom_number_t a2, int value) {
    int ndx = RingStatusIndex(a1, a2);
    ring_number_status[ndx] = value;
    ndx = RingStatusIndex(a2, a1);
    ring_number_status[ndx] = value;
  }

  int StoreAtomTypes(Molecule& m, Atom_Typing_Specification& atom_typing);

  uint32_t
  AtomType(atom_number_t zatom) const {
    return atype[zatom];
  }
};

class Options {
 private:
  int _verbose = 0;
  int _molecules_read = 0;
  int _bad_smiles = 0;
  int _too_many_ring_numbers = 0;

  // Normally an invalid smiles being read stops processing.
  int _ignore_invalid_input_smiles = 0;
  int _invalid_input_smiles = 0;

  // The fragmentation rules are incorporated from here.
  fragment_molecule::MoleculeFragmenter _fragmenter;

  // We can put a fixed isotope at each attachment point.
  // If atom typing is specified, the atom type will be the isotope
  isotope_t _isotope = 0;

  // For diagnostic work it can be helpful to see the parent.
  int _write_parent = 0;

  // Chirality us often not helpful. Much of it may be destroyed on fragmentation.
  int _remove_chirality = 0;

  int _reduce_to_largest_fragment = 0;

  Element_Transformations _etrans;

  Chemical_Standardisation _chemical_standardisation;

  // If set, the atom type number will be the isotope.
  Atom_Typing_Specification _atom_typing;

  // By default, ring numbers are not re-used. This is faster, but results in more
  // ring numbers being used, which may be undesirable.
  int _minimise_ring_numbers_used = 0;

  int _test_unique_smiles = 0;
  int _invalid_safe_smiles = 0;
  int _successful_test = 0;
  int _test_failure = 0;

  extending_resizable_array<int> _breakable_bond_count;

  IWString_and_File_Descriptor _stream_for_fragment_stats;

  // When accumulating fragments, we can impose a maximum atom count limit.
  int _max_atoms = 0;
  // when accumulating fragments, we can impose a limit on the number of
  // non ring atoms.
  int _max_non_ring_atoms = 0;

  absl::flat_hash_map<IWString, dicer_data::DicerFragment> _fragment;

  // Private functions
  int Preprocess(Molecule& m);
  int PerformTest(Molecule& m, PerMoleculeArrays& data, const IWString& SAFE_smiles);
  void AccumulateFragmentStatistics(const resizable_array_p<Molecule>& components,
                                    const IWString& parent_name);
  void AccumulateFragmentStatistics(Molecule& m, const IWString& mname);

 public:
  int Initialise(Command_Line& cl);

  int Process(Molecule& m, int hring, IWString_and_File_Descriptor& output);

  int WriteFragmentStatistics();

  int Report(std::ostream& output) const;
};

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1)) {
      cerr << "Cannot parse -g option\n";
      return 0;
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (!_atom_typing.build(p)) {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 1;
    }
  }

  if (cl.option_present('I')) {
    if (!cl.value('I', _isotope)) {
      cerr << "Options::Initialise:invalid isotope specification (-I)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will label break points with isotope " << _isotope << '\n';
    }
  }

  if (cl.option_present('p')) {
    _write_parent = 1;
    if (_verbose) {
      cerr << "Will write parent molecule\n";
    }
  }

  if (cl.option_present('y')) {
    _minimise_ring_numbers_used = 1;
    if (_verbose) {
      cerr << "Will re-use ring numbers\n";
    }
  }

  if (cl.option_present('z')) {
    _ignore_invalid_input_smiles = 1;
    if (_verbose) {
      cerr << "Will ignore connection table errors on read\n";
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to the largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('T')) {
    if (! _etrans.construct_from_command_line(cl, _verbose, 'T')) {
      cerr << "Cannot initialise element transformations (-T)\n";
      return 0;
    }
  }

  if (cl.option_present('t')) {
    _test_unique_smiles = 1;
    if (_verbose) {
      cerr << "Will check for unique smiles match\n";
    }
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (!_stream_for_fragment_stats.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for fragment statustics '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Fragment statistics written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('M')) {
    const_IWSubstring m;
    for (int i = 0; cl.value('M', m, i); ++i) {
      if (m.numeric_value(_max_atoms) && _max_atoms > 1) {
        if (_verbose) {
          cerr << "Will not accumulate fragments with more than " << _max_atoms
               << " atoms\n";
        }
      } else if (m.starts_with("maxnr=")) {
        m.remove_leading_chars(6);
        if (!m.numeric_value(_max_non_ring_atoms) || _max_non_ring_atoms < 1) {
          cerr << "Invalid maxnr= specification '" << m << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will not accumulate fragments with more than " << _max_non_ring_atoms
               << " non ring atoms\n";
        }
      } else {
        cerr << "Unrecognised -M specification '" << m << "'\n";
        return 0;
      }
    }
  }

  return 1;
}

class NatomsComparitor {
 private:
 public:
  int
  operator()(const Molecule* m1, const Molecule* m2) {
    if (m1->natoms() < m2->natoms()) {
      return -1;
    }
    if (m1->natoms() > m2->natoms()) {
      return 1;
    }
    return 0;
  }
};

PerMoleculeArrays::PerMoleculeArrays(Molecule& m) {
  atoms_in_parent = m.natoms();
  nbrs = new Set_of_Atoms[atoms_in_parent];
  atom_smarts = new IWString[atoms_in_parent];
  ring_number_status = new_int(atoms_in_parent * atoms_in_parent);
  include_atom = new_int(atoms_in_parent, 1);
}

PerMoleculeArrays::~PerMoleculeArrays() {
  delete[] nbrs;
  delete[] atom_smarts;
  delete[] ring_number_status;
  delete[] include_atom;
  if (atype != nullptr) {
    delete[] atype;
  }
}

int
PerMoleculeArrays::StoreAtomTypes(Molecule& m, Atom_Typing_Specification& atom_typing) {
  atype = new uint32_t[m.natoms()];
  return atom_typing.assign_atom_types(m, atype);
}

// We need to keep track of the initial atom numbers
struct Initial {
  atom_number_t initial_atom_number;
};

void
PlaceSmilesSymbol(Molecule& m, atom_number_t zatom, IWString& smi) {
  static constexpr char kOpenSquareBracket = '[';
  static constexpr char kCloseSquareBracket = ']';

  const Atom& a = m[zatom];

  const Element* e = a.element();

  bool square_bracket = false;
  if (e->needs_square_brackets() || a.isotope() || a.formal_charge()) {
    square_bracket = true;
    smi << kOpenSquareBracket;
  }

//  e->append_smiles_symbol(smi, m.is_aromatic(zatom), m.isotope(zatom));
  e->append_smiles_symbol(smi, NOT_AROMATIC, m.isotope(zatom));
  if (square_bracket) {
    if (a.formal_charge() == 0) {
    } else if (a.formal_charge() > 0) {
      smi << '+';
    } else {
      smi << '-';
    }

    smi << kCloseSquareBracket;
  }
}

atom_number_t
InitialAtomNumber(const Molecule& m, atom_number_t zatom) {
  const Initial* ini =
      reinterpret_cast<const Initial*>(m.user_specified_atom_void_ptr(zatom));
  if (ini == nullptr) {
    cerr << "InitialAtomNumber:no void prt " << m.smarts_equivalent_for_atom(zatom)
         << '\n';
    abort();
    return INVALID_ATOM_NUMBER;
  }
  return ini->initial_atom_number;
}

int
AppendSmiles(Molecule& m, PerMoleculeArrays& data, RingNumberControl& rnc,
             IWString& smiles) {
  if (!smiles.empty()) {
    smiles << '.';
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    data.atom_smarts[i].resize_keep_storage(0);
  }

  for (int i = 0; i < matoms; ++i) {
    PlaceSmilesSymbol(m, i, data.atom_smarts[i]);
    const atom_number_t initial_atom_number = InitialAtomNumber(m, i);
    const Set_of_Atoms& nbrsi = data.nbrs[initial_atom_number];
    if (nbrsi.empty()) {
      continue;
    }
    // Join these neighbours to this atom.
    for (atom_number_t nbr : nbrsi) {
      const int ndx = data.RingStatusIndex(nbr, initial_atom_number);

      int ring_number = data.ring_number_status[ndx];
      if (ring_number == 0) {
        ring_number = rnc.GetRing();
        data.SetRingNumber(initial_atom_number, nbr, ring_number);
      } else {
        rnc.OkToReuse(ring_number);
      }

      // cerr << "Between " << initial_atom_number << " and " << nbr << " ring number " <<
      // data.ring_number_status[ndx] << " ndx " << ndx << '\n';
      if (ring_number > 9) {
        data.atom_smarts[i] << '%';
      }
      data.atom_smarts[i] << ring_number;
    }
  }

  Smiles_Information smiles_information(m.natoms());
  smiles_information.allocate_user_specified_atomic_smarts();

  for (int i = 0; i < matoms; ++i) {
    smiles_information.set_user_specified_atomic_smarts(i, data.atom_smarts[i]);
  }

  m.invalidate_smiles();
  const IWString& smt = m.smarts(smiles_information, data.include_atom);

  smiles << smt;

  return 1;
}

int
AppendSmiles(Molecule& m, PerMoleculeArrays& data, int& ring_number, IWString& smiles) {
  if (!smiles.empty()) {
    smiles << '.';
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    data.atom_smarts[i].resize_keep_storage(0);
  }

  for (int i = 0; i < matoms; ++i) {
    PlaceSmilesSymbol(m, i, data.atom_smarts[i]);
    const atom_number_t initial_atom_number = InitialAtomNumber(m, i);
    const Set_of_Atoms& nbrsi = data.nbrs[initial_atom_number];
    if (nbrsi.empty()) {
      continue;
    }
    // Join these neighbours to this atom.
    for (atom_number_t nbr : nbrsi) {
      const int ndx = data.RingStatusIndex(nbr, initial_atom_number);

      if (data.ring_number_status[ndx] == 0) {
        data.SetRingNumber(initial_atom_number, nbr, ring_number);
        ++ring_number;
      }

      // cerr << "Between " << initial_atom_number << " and " << nbr << " ring number " <<
      // data.ring_number_status[ndx] << " ndx " << ndx << '\n';
      if (data.ring_number_status[ndx] > 9) {
        data.atom_smarts[i] << '%';
      }
      data.atom_smarts[i] << data.ring_number_status[ndx];
    }
  }

  Smiles_Information smiles_information(m.natoms());
  smiles_information.allocate_user_specified_atomic_smarts();

  for (int i = 0; i < matoms; ++i) {
    smiles_information.set_user_specified_atomic_smarts(i, data.atom_smarts[i]);
  }

  m.invalidate_smiles();
  const IWString& smt = m.smarts(smiles_information, data.include_atom);

  smiles << smt;

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
  if (_etrans.active()) {
    _etrans.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m, int hring, IWString_and_File_Descriptor& output) {
  Preprocess(m);

  // Added rings are always >= 10
  if (hring < 10) {
    hring = 10;
  }

  resizable_array<int> breakable_bonds;
  int nbreakable = _fragmenter.IdentifyBreakableBonds(m, breakable_bonds);

  ++_breakable_bond_count[nbreakable];

  if (nbreakable == 0) {
    return 1;
  }

  const int matoms = m.natoms();

  PerMoleculeArrays data(m);
  data.parent_name = m.name();
  if (_write_parent) {
    if (_test_unique_smiles) {
      data.parent_smiles = m.unique_smiles();
    } else {
      data.parent_smiles = m.smiles();
    }
  }

  if (_atom_typing.active()) {
    data.StoreAtomTypes(m, _atom_typing);
  }

  std::unique_ptr<Initial[]> initial = std::make_unique<Initial[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    initial[i].initial_atom_number = i;
    m.set_user_specified_atom_void_ptr(i, initial.get() + i);
  }

  for (int bnum : breakable_bonds) {
    // cerr << "Examining bond " << bnum << '\n';
    const Bond* b = m.bond_list()[bnum];
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    data.nbrs[a1] << a2;
    data.nbrs[a2] << a1;
    m.remove_bond_between_atoms(a1, a2);
    if (_isotope) {
      m.set_isotope(a1, _isotope);
      m.set_isotope(a2, _isotope);
    } else if (_atom_typing.active()) {
      // Note we deliberately store the atom type of the connected atom, not the atom
      // itself.
      m.set_isotope(a1, data.AtomType(a2));
      m.set_isotope(a1, data.AtomType(a1));
    }
  }

  set_copy_user_specified_atom_void_ptrs_during_create_subset(1);

  resizable_array_p<Molecule> components;
  m.create_components(components);
  // cerr << "Molecule generates " << components.size() << " components\n";
  if (_stream_for_fragment_stats.is_open()) {
    AccumulateFragmentStatistics(components, data.parent_name);
  }

  static NatomsComparitor nac;
  components.iwqsort(nac);

  // Now create a cross reference from initial atom number to new atoms number
  std::unique_ptr<int[]> xref = std::make_unique<int[]>(matoms);
  int ndx = 0;
  for (const Molecule* c : components) {
    const int catoms = c->natoms();
    for (int i = 0; i < catoms; ++i, ++ndx) {
      const Initial* ini =
          reinterpret_cast<const Initial*>(c->user_specified_atom_void_ptr(i));
      xref[ini->initial_atom_number] = ndx;
    }
  }

#ifdef FOOABA
  for (Molecule * f : components) {
    const int matoms = f->natoms();
    for (int i = 0; i < matoms; ++i) {
      f->set_implicit_hydrogens_known(i, 0);
      f->recompute_implicit_hydrogens(i);
    }
  }
#endif

  // When the smiles is being formed, and a fragment sees that it was bonded
  // to another atom, it does not know if that atom has already been processed,
  // possibly creating a ring number.

  IWString smiles;
  int ring_number = hring;

  if (_minimise_ring_numbers_used) {
    RingNumberControl rnc(ring_number, components.size());
    for (Molecule* c : components) {
      AppendSmiles(*c, data, rnc, smiles);
    }
  } else {
    for (Molecule* c : components) {
      AppendSmiles(*c, data, ring_number, smiles);
    }
  }

  if (_write_parent) {
    output << data.parent_smiles << ' ' << data.parent_name << '\n';
  }

  output << smiles;
  if (!_write_parent) {
    output << ' ' << data.parent_name;
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  if (_test_unique_smiles) {
    return PerformTest(m, data, smiles);
  }

  return 1;
}

int
Options::PerformTest(Molecule& m, PerMoleculeArrays& data, const IWString& SAFE_smiles) {
  Molecule m2;
  if (!m2.build_from_smiles(SAFE_smiles)) {
    ++_invalid_safe_smiles;
    cerr << "Options::PerformTest:invalid SAFE smiles " << SAFE_smiles << '\n';
    ++_invalid_input_smiles;
    return _ignore_invalid_input_smiles;
  }

  if (m2.unique_smiles() == data.parent_smiles) {
    ++_successful_test;
    return 1;
  }

  cerr << "SAFE smiles mismatch, parent " << data.parent_smiles << " SAFE "
       << m2.unique_smiles() << '\n';
  ++_test_failure;

  return 1;
}

void
Options::AccumulateFragmentStatistics(const resizable_array_p<Molecule>& components,
                                      const IWString& parent_name) {
  for (Molecule* m : components) {
    AccumulateFragmentStatistics(*m, parent_name);
  }
}

// Return true if `m` has more than `max_non_ring_atoms`.
int
TooManyNonRingAtoms(Molecule& m, int max_non_ring_atoms) {
  int count = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i)) {
      continue;
    }

    ++count;
    if (count > max_non_ring_atoms) {
      return 1;
    }
  }

  return 0;
}

void
Options::AccumulateFragmentStatistics(Molecule& m, const IWString& mname) {
  if (_max_atoms > 0 && m.natoms() > _max_atoms) {
    return;
  }
  if (_max_non_ring_atoms > 0 && TooManyNonRingAtoms(m, _max_non_ring_atoms)) {
    return;
  }

  // Save a copy.
  IWString smi = m.smiles();

  const IWString& usmi = m.unique_smiles();
  auto iter = _fragment.find(usmi);
  if (iter == _fragment.end()) {
    dicer_data::DicerFragment proto;
    proto.set_n(1);
    proto.set_smi(smi.AsString());
    proto.set_par(mname.AsString());
    if (_isotope) {
      proto.set_iso(dicer_data::ATT);
    }
    _fragment.emplace(usmi, std::move(proto));
  } else {
    auto n = iter->second.n();
    iter->second.set_n(n + 1);
  }
}

int
Options::Report(std::ostream& output) const {
  output << "Processed " << _molecules_read << " molecules\n";
  for (int i = 0; i < _breakable_bond_count.number_elements(); ++i) {
    if (_breakable_bond_count[i]) {
      output << _breakable_bond_count[i] << " molecules had " << i
             << " breakable bonds\n";
    }
  }

  if (_test_unique_smiles) {
    output << _invalid_safe_smiles << " invalid SAFE smiles generated\n";
    output << _successful_test << " successful unique smiles matches\n";
    output << _test_failure << " smiles mismatches\n";
  }
  return 1;
}

int
Options::WriteFragmentStatistics() {
  if (!_stream_for_fragment_stats.is_open()) {
    cerr << "Options::WriteFragmentStatistics:stream not open\n";
    return 0;
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  for (const auto& [usmi, proto] : _fragment) {
    if (!printer.PrintToString(proto, &buffer)) {
      cerr << "Options::WriteFragmentStatistics write '" << proto.ShortDebugString()
           << "'\n";
      return 0;
    }

    _stream_for_fragment_stats << buffer;
    _stream_for_fragment_stats << '\n';
    _stream_for_fragment_stats.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Mol2SAFE(Options& options, Molecule& m, const int hring,
         IWString_and_File_Descriptor& output) {
  return options.Process(m, hring, output);
}

int
Mol2SAFEInner(Options& options, const IWString& buffer,
              IWString_and_File_Descriptor& output) {
  Molecule m;
  if (!m.build_from_smiles(buffer)) {
    cerr << "Mol2SAFEInner:invalid smiles\n";
    return 0;
  }

  std::optional<int> hring = lillymol::HighestRingNumber(buffer);
  if (!hring) {
    return 1;
  }

  return Mol2SAFE(options, m, *hring, output);
}

int
Mol2SAFE(Options& options, iwstring_data_source& input,
         IWString_and_File_Descriptor& output) {
  IWString buffer;
  while (input.next_record(buffer)) {
    if (!Mol2SAFEInner(options, buffer, output)) {
      cerr << "Mol2SAFEInner:error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Mol2SAFE(Options& options, const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Mol2SAFE:cannot open '" << fname << "'\n";
    return 0;
  }

  return Mol2SAFE(options, input, output);
}

int
Mol2SAFE(int argc, char** argv) {
  Command_Line cl(argc, argv, "vI:pcltg:S:M:P:yzT:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (!Mol2SAFE(options, fname, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (cl.option_present('S')) {
    options.WriteFragmentStatistics();
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace mol2safe

int
main(int argc, char** argv) {
  int rc = mol2safe::Mol2SAFE(argc, argv);

  return rc;
}
