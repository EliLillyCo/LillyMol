// Systematic enumeration around a set of molecles.
// Takes a set of molecules and a fragment library,
// and adds all members of the fragment library to
// available sites in the molecules.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <random>

#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools/enumeration.pb.h"

namespace substituent_enumeration {

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
  cerr << R"(
Performs systematic enumeration of all available sites in a molecule.
Uses a library of fragments to add. These must be either smiles or
textproto represetation of DicerFragment protos.
Fragments added must be marked with an isotope. Only the first isotopic atom
in each fragment is used.
 -C <fname>             substituent_enumeration.EnumerationConfig textproto file with options.
 -F <fname>             Read file of substituents as smiles.
 -F PROTO:<fname>       Read file of substituents as DicerFragment textproto.
 -F                     If any -F directive ends with ',arom' those fragments are only added to aroamtic sites.
 -m <natoms>            Only use substituents with at least <matoms>
 -M <natoms>            Only use substituents with at most <matoms>
 -y <query>             Specify the atom(s) in the starting molecule to which fragments can    be added.
 -n <query>             Specify the atom(s) in the starting molecule to which fragments cannot be added.
                        only one of either -y -n can be specified.
 -p                     write the parent molecule before any variants
 -Y ...                 other options, enter '-Y help' for info

 -v                     verbose output
)";
// clang-format on

  ::exit(rc);
}

void
DisplayDashYOptions(std::ostream& output) {
  // clang-format off
  output << " -Y xiso           remove isotopic labels from product molecules\n";
  output << " -Y dmv=<n>        do not process any molecule that has more than <n> reactive sites\n";
  output << " -Y dns=<n>        if a molecule has more than <n> sites, randomly downsample to limit to approximately <n>\n";
  output << " -Y valence        check for OK valences in product molecules, and discard bad\n";
  output << " -Y stop=<n>       stop processing a molecule once it has generated <n> products\n";
  output << "                   note that this will be a function of the size of the fragment library as well as sites in the molecule\n";
  // clang-format on

  ::exit(0);
}

// Keep track of the sidechains being added.
// A molecule and the attachment point, and
// possibly with the number of times this is found in
// the originating collection.
class Fragment {
  private:
     Molecule* _frag = nullptr;

  // The atom by which this fragment is attached.
     atom_number_t _zatom = kInvalidAtomNumber;

  // What kind of atom is it, and is it a halogen
     atomic_number_t _atomic_number = kInvalidAtomicNumber;
     bool _is_halogen = false;

   // Some atomic properties useful for avoiding forming
   // undesirable adjacencies.
     int _aromatic = 0;
     int _ring_bond_count = 0;
     int _single_bond_count = 0;
     int _double_bond_count = 0;
     int _triple_bond_count = 0;
     int _aryl = 0;
     int _vinyl = 0;
     int _attached_carbon_count = 0;
     int _attached_heteroatom_count = 0;
     int _attached_halogen_count = 0;
     int _attached_nitrogen_count = 0;
     int _attached_oxygen_count = 0;
     int _saturated = 1;
     int _saturated_oxygen_neighbour = 0;
     int _saturated_nitrogen_neighbour = 0;

     // If this is the carbon atom of an amide, both
     // of these will be set.
     int _carbon_of_amide = 0;
     int _nitrogen_of_amide = 0;

     // Will be true if the labelled atom is the O of an acid.
     int _is_acid = 0;

     // These are whole molecule properties.
     int _positive_charge = 0;
     int _negative_charge = 0;
     int _halogen_count = 0;

     int _rotbonds = 0;

     int _nitro_count = 0;

  // If this comes from a DicerFragment proto, the number
  // of instances.
     uint32_t _count = 0;

     // If set, this fragment will only attach to aromatic atoms in the
     // starting molecule.
     int _aromatic_only = 0;

  // private functions
    void EstablishAtomicProperties(Molecule& m, atom_number_t zatom);
    int IsAcid(Molecule& m, atom_number_t zatom);

  public:
    Fragment();
    ~Fragment();

    // Any caller should not change the molecule since that could
    // invalidate stored values. But methods like nrings() are non const.
    Molecule* frag() const {
      return _frag;
    }

    int SetFragment(Molecule* f);

    int DetermineFormalCharges(Charge_Assigner& charge_assigner);

    int is_halogen() const {
      return _is_halogen;
    }

    atom_number_t zatom() const {
      return _zatom;
    }

    atomic_number_t atomic_number() const {
      return _atomic_number;
    }

    uint32_t count() const {
      return _count;
    }
    void set_count(uint32_t s) {
      _count = s;
    }

    int rotbonds() const {
      return _rotbonds;
    }
    void set_rotbonds(int s) {
      _rotbonds = s;
    }

    int is_acid() const {
      return _is_acid;
    }

    int aromatic_only() const {
      return _aromatic_only;
    }
    void set_aromatic_only(int s) {
      _aromatic_only = s;
    }

    int saturated() const {
      return _saturated;
    }

    int aromatic() const {
      return _aromatic;
    }
    int aryl() const {
      return _aryl;
    }
    int vinyl() const {
      return _vinyl;
    }

    int saturated_heteroatom_neighbour() const {
      if (_saturated_oxygen_neighbour) {
        return 1;
      }
      if (_saturated_nitrogen_neighbour) {
        return 1;
      }

      return 0;
    }

    int attached_carbon_count() const {
      return _attached_carbon_count;
    }
    int attached_nitrogen_count() const {
      return _attached_nitrogen_count;
    }
    int attached_oxygen_count() const {
      return _attached_oxygen_count;
    }

    int positive_charge() const {
      return _positive_charge;
    }
    int negative_charge() const {
      return _negative_charge;
    }
    int halogens() const {
      return _halogen_count;
    }
    int nitro_count() const {
      return _nitro_count;
    }
    int carbon_of_amide() const {
      return _carbon_of_amide;
    }
    int nitrogen_of_amide() const {
      return _nitrogen_of_amide;
    }

    // Return true if the atom is a charged aromatic nitrogen.
    int IsQuaternaryAryl();
};

Fragment::Fragment() {
}

Fragment::~Fragment() {
  if (_frag != nullptr) {
    delete _frag;
  }
}

bool
IsHalogem(atomic_number_t z) {
  if (z == 6) {
    return 0;
  }

  if (z == 7) {
    return 0;
  }

  if (z == 8) {
    return 0;
  }

  if (z == 9) {
    return 1;
  }
  if (z == 17) {
    return 1;
  }
  if (z == 35) {
    return 1;
  }
  if (z == 53) {
    return 1;
  }

  return 0;
}


// If we can identify an attachment point in `f` take ownership
// of `f`, otherwise return 0 and it will be up to the caller to
// dispose of `f`.
int
Fragment::SetFragment(Molecule* f) {
  const int matoms = f->natoms();

  for (int i = 0; i < matoms; ++i) {
    if (f->isotope(i) == 0) {
      continue;
    }

    _zatom = i;
    _atomic_number = f->atomic_number(i);
    _is_halogen = IsHalogem(f->atomic_number(i));
    EstablishAtomicProperties(*f, i);
    _frag = f;
    return 1;
  }

  return 0;
}

int
CountHalogens(const Molecule& m) {
  const int matoms = m.natoms();
  int halogens = 0;
  for (int i = 0; i < matoms; ++i) {
    atomic_number_t z = m.atomic_number(i);
    if (z < 9) {
      continue;
    }
    if (z == 9 || z == 17 || z == 35 || z == 53) {
      ++halogens;
    }
  }

  return halogens;
}

int
CountNitro(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& n = m[i];
    if (n.atomic_number() != 7) {
      continue;
    }
    if (n.ncon() != 3) {
      continue;
    }
    if (m.ring_bond_count(i)) {
      continue;
    }

    int doubly_bonded_oxygen_count = 0;
    for (const Bond* b : n) {
      if (! b->is_double_bond()) {
        continue;
      }
      const atom_number_t o = b->other(i);
      if (m.atomic_number(o) == 8) {
        ++doubly_bonded_oxygen_count;
      }
    }
    if (doubly_bonded_oxygen_count == 2) {
      return 1;
    }
  }

  return 0;
}

int
Fragment::DetermineFormalCharges(Charge_Assigner& charge_assigner) {
  const int matoms = _frag->natoms();

  std::unique_ptr<formal_charge_t[]> q = std::make_unique<formal_charge_t[]>(matoms);
  std::fill_n(q.get(), matoms, 0);

  if (charge_assigner.process(*_frag, q.get()) == 0) {
    return 1;
  }

  for (int i = 0; i < matoms; ++i) {
    if (q[i] == 0) {
    } else if (q[i] > 0) {
      ++_positive_charge;
    } else {
      ++_negative_charge;
    }
  }

  _halogen_count = CountHalogens(*_frag);

  return 1;
}

// Return true if `zatom` is the oxygen of an acid.
int
Fragment::IsAcid(Molecule& m, atom_number_t zatom) {
  if (m.atomic_number(zatom) != 8) {
    return 0;
  }

  if (m.ncon(zatom) != 1) {
    return 0;
  }
  if (0 == _vinyl) {
    return 0;
  }

  const atom_number_t carbon = m.other(zatom, 0);
  if (m.atomic_number(carbon) == 6) {
  } else if (m.atomic_number(carbon) == 16) {
  } else {
    return 0;
  }

  for (const Bond* b : m[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(carbon);
    if (m.atomic_number(o) == 8) {
      return 1;
    }
    if (m.atomic_number(o) == 16) {
      return 1;
    }
  }

  return 0;
}

// Examine the atoms attached to `zatom` and if any of them are saturated
// nitrogen or oxygen atoms, update the _saturated_*_neighbour class variable.
void
Fragment::EstablishAtomicProperties(Molecule& m, atom_number_t zatom) {
  // Only compute this for some atom types.
  m.compute_aromaticity_if_needed();

  _aromatic = m.is_aromatic(zatom);

  // these help with amide perception.
  int doubly_bonded_oxygen = kInvalidAtomNumber;
  int singly_bonded_nitrogen = kInvalidAtomNumber;

  for (const Bond* b : m[zatom]) {
    if (b->nrings()) {
      ++_ring_bond_count;
    }

    if (b->is_aromatic()) {
      continue;
    }

    const atom_number_t o = b->other(zatom);
    const atomic_number_t z = m.atomic_number(o);

    if (b->is_single_bond()) {
      ++_single_bond_count;
      if (z == 7) {
        singly_bonded_nitrogen = o;
      }
    } else if (b->is_double_bond()) {
      _saturated = 0;
      ++_double_bond_count;
      if (z == 8) {
        doubly_bonded_oxygen = o;
      }
    } else {
      _saturated = 0;
      ++_triple_bond_count;
    }

    const int saturated = m.saturated(o);

    if (m.is_aromatic(o)) {
      ++_aryl;
    } else if (! saturated) {
      ++_vinyl;
    }

    if (z == 6) {
      ++_attached_carbon_count;
      continue;
    } 
    
    ++_attached_heteroatom_count;

    if (z == 7) {
      ++_attached_nitrogen_count;
      if (saturated) {
        ++_saturated_nitrogen_neighbour;
      }
    } else if (z == 8) {
      ++_attached_oxygen_count;
      if (saturated) {
        ++_saturated_oxygen_neighbour;
      }
    } else if (IsHalogem(z)) {
      ++_attached_halogen_count;
    }
  }

  if (singly_bonded_nitrogen != kInvalidAtomNumber &&
      doubly_bonded_oxygen != kInvalidAtomNumber) {
    _carbon_of_amide = 1;
  }

  _nitro_count = CountNitro(m);

  _is_acid = IsAcid(m, zatom);
}

int
Fragment::IsQuaternaryAryl() {
  if (_frag->atomic_number(_zatom) != 7) {
    return 0;
  }
  if (_frag->formal_charge(_zatom) != 1) {
    return 0;
  }
  if (! _frag->is_aromatic(_zatom)) {
    return 0;
  }

  return 1;
}

// Data we harvest from a molecule.
struct MoleculeData {
  int positive = 0;
  int negative = 0;
  int aromatic_ring_count = 0;
  int nrings = 0;

  // includes F.
  int halogens = 0;
  // only set if a limit on rotatable bonds is specified.
  int rotbonds = 0;

  int nitro_count = 0;

  int Initialise(Molecule& m, Charge_Assigner& charge_assigner);
};

int
MoleculeData::Initialise(Molecule& m, Charge_Assigner& charge_assigner) {
  m.compute_aromaticity_if_needed();

  nrings = m.nrings();
  for (const Ring* r : m.sssr_rings()) {
    if (r->is_aromatic()) {
      ++aromatic_ring_count;
    }
  }

  if (charge_assigner.active()) {
    const int matoms = m.natoms();

    std::unique_ptr<formal_charge_t[]> q = std::make_unique<formal_charge_t[]>(matoms);
    charge_assigner.process(m, q.get());
    for (int i = 0; i < matoms; ++i) {
      if (q[i] == 0) {
      } else if (q[i] > 0) {
        ++positive;
      } else {
        ++negative;
      }
    }
  }

  halogens = CountHalogens(m);

  nitro_count = CountNitro(m);

  return 1;
}


int
AttachedTo(const Molecule& m,
           atom_number_t zatom,
           atomic_number_t z) {
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (m.atomic_number(o) == z) {
      return 1;
    }
  }

  return 0;
}

int
AttachedToSaturatedON(const Molecule& m,
                atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);

    const atomic_number_t zo = m.atomic_number(o);
    if ((zo == 7 || zo == 8) && m.saturated(o)) {
      return 1;
    }
  }

  return 0;
}

// Return true if joining `fragment` to `m` via `zatom` would likely
// result in formation of an NCN motif. Note that this is stricter than
// an aminal, where the Nitrogen atoms must be saturated.
// We need to check two cases.
// molecule NC - N fragment
// molecule N  - CN fragment
// Also look for O-C-N
//
int
IsNCN(Molecule& m, atom_number_t zatom, const Fragment& fragment) {
  if (m.atomic_number(zatom) == 6 &&
      (fragment.atomic_number() == 7 || fragment.atomic_number() == 8)) {
    if (! m.saturated(zatom)) {
      return 0;
    }
    return AttachedToSaturatedON(m, zatom);
  } else if ((m.atomic_number(zatom) == 7 || m.atomic_number(zatom) == 8) &&
             fragment.atomic_number() == 6) {
    if (! fragment.saturated()) {
      return 0;
    }
    if (fragment.attached_nitrogen_count()) {
      return 1;
    }
    if (fragment.attached_oxygen_count()) {
      return 1;
    }
  }

  return 0;
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    // A list of fragments that are to be added.
    resizable_array_p<Fragment> _fragment;

    EnumerationConfig _config;

    // We can specify atoms where substitution is OK or
    // atoms where substitution is not allowed.
    resizable_array_p<Substructure_Query> _ok_attach;
    resizable_array_p<Substructure_Query> _no_attach;

    quick_rotbond::QuickRotatableBonds _rotbond;

    // Once the new molecule is formed, we can do substructure
    //searches to discard known problems.
    resizable_array_p<Substructure_Query> _discard_products_containing;
    int _product_discarded_by_query = 0;

    // By default, we do not check the valence of the products.
    int _products_with_bad_valence = 0;

    // Triggered by _config.stop_generating_if_more_than()
    int _stopped_for_too_many_products = 0;

    // Do not generate duplicates.
    IW_STL_Hash_Set _seen;

    // Triggered by _config.discard_if_too_many_sites()
    int _discarded_for_too_many_sites = 0;

    // triggered by _config.max_atoms_in_product();
    int _product_too_many_atoms = 0;

    // Triggered by _config.max_rings_in_product()
    int _product_too_many_rings = 0;

    // Triggered by _config.max_halogens_in_product()
    int _product_too_many_halogens = 0;

    // Triggered by _config.max_rotbonds_in_product()
    int _product_too_many_rotatable_bonds = 0;

    std::default_random_engine _generator;

    int _molecules_read = 0;

    Charge_Assigner _charge_assigner;

    // The number of sites in each starting molecule.
    extending_resizable_array<int> _nsites;
    // The number of molecules generated per starting molecule.
    // Generally this will be roughly nsites * number of fragments,
    // although that is an upper value.
    extending_resizable_array<int> _variants_generated;

    // Private functions

    int InitialiseFromProto();
    int BuildQueries(resizable_array_p<Substructure_Query>& queries,
                      const google::protobuf::RepeatedPtrField<std::string>& tokens);

    int IdentifyAttachmentPoints(Molecule& m, int * attachment_point);
    int Seen(Molecule& m);
    int GenerateVariants(Molecule& m, const MoleculeData& mdata, atom_number_t zatom,
                          IWString_and_File_Descriptor& output);
    int OkProperties(Molecule& m);

    int ReadFragmentFromTextProto(const dicer_data::DicerFragment& proto, int aromatic_only);
    int ReadFragmentsFromTextProto(IWString& fname, int aromatic_only);
    int ReadFragmentsFromTextProto(iwstring_data_source& input, int aromatic_only);
    int ReadFragmentFromTextProto(const const_IWSubstring& buffer, int aromatic_only);

    int ReadFragmentsAsSmiles(IWString& fname, int aromatic_only);
    int ReadFragmentsAsSmiles(data_source_and_type<Molecule>& input, int aromatic_only);

    int OkFragment(Molecule& frag) const;
    int ContainsRejectedSubstructure(Molecule& m);

    int OkToBeCombined(Molecule& m,
                const MoleculeData& mdata,
                atom_number_t zatom,
                Fragment& frag);

  public:
    Options();

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
  _molecules_read = 0;
  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    std::optional<EnumerationConfig> proto =
        iwmisc::ReadTextProtoCommentsOK<EnumerationConfig>(fname);
    if (! proto) {
      cerr << "Options::Initialise:cannot read config proto '" << fname << "'\n";
      return 0;
    }

    _config = std::move(*proto);

    if (! InitialiseFromProto()) {
      cerr << "Options::Initialise:cannot initialise state from proto\n";
      cerr << proto->ShortDebugString() << '\n';
      return 0;
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    _config.set_reduce_to_largest_fragment(true);
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _config.set_remove_chirality(true);
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('m')) {
    int m;
    if (! cl.value('m', m)) {
      cerr << "Invalid min atoms in fragment (-c)\n";
      return 0;
    }

    _config.set_min_atoms_in_fragment(m);

    if (_verbose) {
      cerr << "Will only use fragments with at least " << m << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    int m;
    if (! cl.value('M', m)) {
      cerr << "Invalid max atoms in fragment (-c)\n";
      return 0;
    }

    _config.set_max_atoms_in_fragment(m);
    if (_verbose) {
      cerr << "Will only use fragments with at most " << m << " atoms\n";
    }
  }

  if (cl.option_present('F')) {
    IWString fname;
    for (int i = 0; cl.value('F', fname, i); ++i) {
      int only_aromatic = 0;
      if (fname.ends_with(",arom")) {
        fname.chop(5);
        only_aromatic = 1;
      }

      if (fname.starts_with("PROTO:")) {
        fname.remove_leading_chars(6);
        if (! ReadFragmentsFromTextProto(fname, only_aromatic)) {
          cerr << "Options::Initialise:cannot read testproto file '" << fname << "'\n";
          return 0;
        }
      } else {
        if (! ReadFragmentsAsSmiles(fname, only_aromatic)) {
          cerr << "Options::Initialise:cannot read smiles file '" << fname << "'\n";
          return 0;
        }
      }
    }
  }

  if (_fragment.empty()) {
    cerr << "Options::Initialise:no fragments - use the -F option to specify\n";
    Usage(0);
  }

  if (cl.option_present('y')) {
    if (! process_queries(cl, _ok_attach, _verbose, 'y')) {
      cerr << "Options::Initialise:cannot assemble ok attachment point queries (-y)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _ok_attach.size() << " ok attachment point queries\n";
    }
  }

  if (cl.option_present('n')) {
    if (! process_queries(cl, _no_attach, _verbose, 'n')) {
      cerr << "Options::Initialise:cannot assemble not ok attachment point queries (-n)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _no_attach.size() << " queries for excluding attachment points\n";
    }
  }

  if (cl.option_present('p')) {
    _config.set_write_parent(true);
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('Y')) {
    IWString y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "xiso") {
        _config.set_remove_isotopes(true);
        if (_verbose) {
          cerr << "Will remove isotopic labels from products\n";
        }
      } else if (y.starts_with("dmv=")) {
        y.remove_leading_chars(4);
        int d;
        if (! y.numeric_value(d)) {
          cerr << "Invalid discard if too many sites directive 'dmv=" << y << "'\n";
          return 0;
        }
        _config.set_discard_if_too_many_sites(d);
      } else if (y.starts_with("dns=")) {
        y.remove_leading_chars(4);
        uint32_t d;
        if (! y.numeric_value(d)) {
          cerr << "Invalid downsample threshold 'dns=" << y << "'\n";
          return 0;
        }
        _config.set_downsample_threshold(d);
        if (_verbose) {
          cerr << "Will downsample molecules with more than " << d << " sites\n";
        }
      } else if (y == "valence") {
        _config.set_check_valences(true);
        if (_verbose) {
          cerr << "Will discard products with bad valences\n";
        }
      } else if (y.starts_with("stop=")) {
        y.remove_leading_chars(5);
        uint32_t s;
        if (! y.numeric_value(s)) {
          cerr << "Invalid stop= directive '" << y << "'\n";
          return 0;
        }
        _config.set_stop_generating_if_more_than(s);
        if (_verbose) {
          cerr << "Will abandon an input molecule if it generates more than " << s << " products\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  return 1;
}

// The proto contains various things that we need to convert into our own
// data structures, charge assigner, queries...
int
Options::InitialiseFromProto() {
  if (_config.ok_attach().size()) {
    if (! BuildQueries(_ok_attach, _config.ok_attach())) {
      cerr << "Options::InitialiseFromProto:cannot initialise ok_attach queries\n";
      return 0;
    }
  }

  if (_config.no_attach().size()) {
    if (! BuildQueries(_no_attach, _config.no_attach())) {
      cerr << "Options::InitialiseFromProto:cannot initialise no_attach queries\n";
      return 0;
    }
  }

  if (_ok_attach.size() && _no_attach.size()) {
    cerr << "Options::InitialiseFromProto:cannot have both ok attach and no attach queries\n";
    return 0;
  }

  if (_config.has_charge_assigner()) {
    IWString tmp(_config.charge_assigner());
    if (! _charge_assigner.build(tmp)) {
      cerr << "Options::InitialiseFromProto:cannot initialise charge assigner '" << tmp << "'\n";
      return 0;
    }
    _charge_assigner.set_apply_charges_to_molecule(0);
  }

  if (_config.discard_products_containing_size()) {
    if (! BuildQueries(_discard_products_containing, _config.discard_products_containing())) {
      cerr << "Options::Initialise:cannot read queries for discards\n";
      return 0;
    }
  }

  return 1;
}

int
Options::BuildQueries(resizable_array_p<Substructure_Query>& queries,
                      const google::protobuf::RepeatedPtrField<std::string>& tokens) {
  for (const std::string& token : tokens) {
    IWString tmp(token);
    if (! process_cmdline_token('*', tmp, queries, 0)) {
      cerr << "Options::BuildQueries:cannot process '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_discarded_for_too_many_sites) {
    output << _discarded_for_too_many_sites << " input molecules discarded for more than " << _config.discard_if_too_many_sites() << " sites\n";
  }
  if (_config.check_valences()) {
    output << "Discarded " << _products_with_bad_valence << " products with bad valence\n";
  }
  if (_config.stop_generating_if_more_than() > 0) {
    output << _stopped_for_too_many_products << " enumerations stopped for generating " <<  _config.stop_generating_if_more_than() << " products\n";
  }
  if (_discard_products_containing.size()) {
    output << _product_discarded_by_query << " products discarded for hitting rejection query\n";
  }
  if (_config.has_max_atoms_in_product()) {
    output << _product_too_many_atoms << " products with more than " << _config.max_atoms_in_product() << " atoms\n";
  }
  if (_config.has_max_rings_in_product()) {
    output << _product_too_many_rings << " products with more than " << _config.max_rings_in_product() << " rings\n";
  }
  if (_config.has_max_halogens_in_product()) {
    output << _product_too_many_halogens << " products with more than " << _config.max_halogens_in_product() << " halogen atoms\n";
  }
  if (_config.has_max_rotbonds_in_product()) {
    output << _product_too_many_rotatable_bonds << " products with more than " << _config.max_rotbonds_in_product() << " rotatable bonds\n";
  }

  Accumulator_Int<int> acc;
  for (int i = 0; i < _nsites.number_elements(); ++i) {
    if (_nsites[i]) {
      output << _nsites[i] << " molecules had " << i << " sites for substitution\n";
      acc.extra(i, _nsites[i]);
    }
  }
  output << " ave " << acc.average() << " sites per starting molecule\n";
  for (int i = 0; i < _variants_generated.number_elements(); ++i) {
    if (_variants_generated[i]) {
      output << _variants_generated[i] << " molecules generated " << i << " variants\n";
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_config.reduce_to_largest_fragment()) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_config.remove_chirality()) {
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

// We don't put extra bonds on Carbon atoms that already have 3 connection,
// even if they do have an available H.
int
Is3ConnectedCarbon(Molecule& m, atom_number_t zatom) {
  if (m.is_aromatic(zatom)) {
    return 0;
  }

  const Atom& a = m[zatom];

  if (a.atomic_number() != 6) {
    return 0;
  }

  return a.ncon() == 3;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  const int matoms = m.natoms();

  if (matoms == 0) {
    return 1;
  }

  if (_config.write_parent()) {
    static constexpr char kSep = ' ';

    output << m.smiles() << kSep << m.name() << '\n';
  }

  MoleculeData mdata;
  mdata.Initialise(m, _charge_assigner);
  if (_config.has_max_rotbonds_in_product()) {
    mdata.rotbonds = _rotbond.Process(m);
  }

  std::unique_ptr<int[]> attachment_point = std::make_unique<int[]>(matoms);
  if (_ok_attach.size() || _no_attach.size()) {
    if (! IdentifyAttachmentPoints(m, attachment_point.get())) {
      cerr << "Options::Process:cannot identify attachment point atoms\n";
      return 0;
    }
  } else {
    std::fill_n(attachment_point.get(), matoms, 1);
  }

  for (int i = 0; i < matoms; ++i) {
    if (m.formal_charge(i)) {
      attachment_point[i] = 0;
    } else if (m.hcount(i) == 0) {
      attachment_point[i] = 0;
    } else if (Is3ConnectedCarbon(m, i)) {
      attachment_point[i] = 0;
    }
  }

  const uint32_t nsites = std::count(attachment_point.get(), attachment_point.get() + matoms, 1);

  ++_nsites[nsites];

  if (nsites == 0) {
    return 1;
  }

  if (_config.discard_if_too_many_sites() > 0) {
    if (nsites > _config.discard_if_too_many_sites()) {
      ++_discarded_for_too_many_sites;
      return 0;
    }
  }
  // cerr << "Processing '" << m.name() << "' nsites " << nsites << "\n";

  std::unique_ptr<std::bernoulli_distribution> rng;
  if (nsites > _config.downsample_threshold()) {
    float fraction = iwmisc::Fraction<double>(_config.downsample_threshold(), nsites);
    rng = std::make_unique<std::bernoulli_distribution>(fraction);
  }

  uint32_t rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (! attachment_point[i]) {
      continue;
    }

    if (rng && ! (*rng)(_generator)) {
      continue;
    }

    rc += GenerateVariants(m, mdata, i, output);

    if(_config.has_stop_generating_if_more_than() && rc > _config.stop_generating_if_more_than()) {
      ++_stopped_for_too_many_products;
      break;
    }
  }

  output.write_if_buffer_holds_more_than(4092);

  ++_variants_generated[rc];

  return rc;
}

int
IdentifyAtoms(Molecule_to_Match& target,
              resizable_array_p<Substructure_Query>& queries,
              int* attachment_point,
              int flag) {
  int rc = 0;
  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults)) {
      sresults.each_embedding_set_vector(attachment_point, flag);
      ++rc;
    }
  }

  return rc;
}

int
Options::IdentifyAttachmentPoints(Molecule& m, int * attachment_point) {
  const int matoms = m.natoms();

  Molecule_to_Match target(&m);

  if (_ok_attach.size()) {
    std::fill_n(attachment_point, matoms, 0);
    return IdentifyAtoms(target, _ok_attach, attachment_point, 1);
  } else {
    std::fill_n(attachment_point, matoms, 1);
    return IdentifyAtoms(target, _no_attach, attachment_point, 0);
  }
}

// Return true if it is OK to join `zatom` in `m` to `frag`.
int
OkBondFormation(Molecule& m,
                atom_number_t zatom,
                const Fragment& frag) {
  if (frag.is_halogen() && ! m.is_aromatic(zatom)) {
    return 0;
  }

  if (frag.atomic_number() == 6) {
    return 1;
  }

  const atomic_number_t z = m.atomic_number(zatom);
  if (z == 6) {
    return 1;
  }

  // Must be two heteratoms, which we do not allow.
  return 0;
}

int
Options::GenerateVariants(Molecule& m, 
                          const MoleculeData& mdata,
                          atom_number_t zatom,
                          IWString_and_File_Descriptor& output) {
  const int initial_matoms = m.natoms();

  static constexpr char kSep = ' ';

  static IWString percents(" %% ");

  // cerr << "For atom " << zatom << " in " << m.name() << " will generate " << _fragment.size() << " variants\n";
  int rc = 0;
  // `frag` is actually const, but some underlying methods may call some non-cost
  // methods on frag.frag().
  for (Fragment* frag : _fragment) {
    if (! OkBondFormation(m, zatom, *frag)) {
      continue;
    }

    if (frag->aromatic_only() && ! m.is_aromatic(zatom)) {
      continue;
    }

    if (!OkToBeCombined(m, mdata, zatom, *frag)) {
      continue;
    }

    Molecule mcopy(m);
    mcopy += *frag->frag();
    mcopy.add_bond(zatom, initial_matoms + frag->zatom(), SINGLE_BOND);
    mcopy.unset_all_implicit_hydrogen_information(zatom);
    mcopy.unset_all_implicit_hydrogen_information(initial_matoms + frag->zatom());
    if (_config.has_isotope_at_join()) {
      mcopy.set_isotope(zatom, _config.isotope_at_join());
    }

    if (_config.remove_isotopes()) {
      mcopy.unset_isotopes();
    }

    if (Seen(mcopy)) {
      continue;
    }

    if (! OkProperties(mcopy)) {
      continue;
    }

    mcopy.invalidate_smiles();

    if (_config.check_valences() && ! mcopy.valence_ok()) {
      ++_products_with_bad_valence;
      cerr << "Bad valence " << mcopy.smiles() << ' ' << m.name() << '\n';
      continue;
    }

    if (! _discard_products_containing.empty() &&
        ContainsRejectedSubstructure(mcopy)) {
      continue;
    }

    output << mcopy.smiles() << kSep << m.name() << percents << frag->frag()->name();
    if (frag->count() > 0) {
      output << kSep << frag->count();
    }
    output << '\n';
    ++rc;
  }

  output.write_if_buffer_holds_more_than(8192);

  return rc;
}

// If the fragment is an acid and the molecule is aromatic.
int
WouldFormPhenolicEster(Molecule& m, 
                const MoleculeData& mdata,
                atom_number_t zatom,
                Fragment& frag) {
  if (! m.is_aromatic(zatom)) {
    return 0;
  }

  return frag.is_acid();
}

int
WouldFormAdjacentRings(Molecule& m, atom_number_t zatom,
                  Fragment& frag) {
  if (m.ring_bond_count(zatom) == 0) {
    return 0;
  }

  return frag.frag()->ring_bond_count(frag.zatom());
}

int
WouldFormBiphenyl(Molecule& m, atom_number_t zatom,
                  Fragment& frag) {
  if (! m.is_aromatic(zatom)) {
    return 0;
  }

  return frag.aromatic();
}

int
WouldContainTooManyNitros(Molecule& m, 
                          const MoleculeData& mdata,
                          Fragment& frag) {
  int nnitro = mdata.nitro_count + frag.nitro_count();
  return nnitro > 1;
}

// Looking for a Nitrogen atom being added to an aromatic
// where the N atom has no nearby unsaturation.
// can be either molecule -N  - a- fragment
// can be either molecule -a  - N- fragment
int
WouldFormAnilineNh(Molecule& m,
                 const MoleculeData& mdata,
                 atom_number_t zatom,
                 Fragment& frag) {
  const Molecule* f = frag.frag();

  // THe case of molecule-a -  N-fragment
  if (m.is_aromatic(zatom) &&
      frag.atomic_number() == 7 && f->ncon(frag.zatom()) == 1 &&
      frag.saturated() && 
      frag.aryl() == 0 && frag.vinyl() == 0) {
    return 1;
  }

  // The case of molecule-N a- fragment
  if (! frag.aromatic()) {
    return 0;
  }

  if (m.atomic_number(zatom) != 7) {
    return 0;
  }
  if (m.ncon(zatom) != 1) {
    return 0;
  }

  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      return 0;
    }

    atom_number_t o = b->other(zatom);
    if (m.atomic_number(o) != 6) {
      return 0;
    }
    if (! m.saturated(o)) {
      return 0;
    }
  }

  return 1;
}


// Return true if `zatom` is doubly bonded to an aliphatic carbon atom.
int
IsEnamineCarbon(Molecule& m, atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(zatom);
    if (m.is_aromatic(o)) {
      continue;
    }
    if (m.atomic_number(o) == 6) {
      return 1;
    }
  }

  return 0;
}

int
IsCarbonOfAmide(Molecule& m, atom_number_t zatom) {
  const Atom& a = m[zatom];
  if (a.atomic_number() != 6) {
    return 0;
  }

  atom_number_t singly_bonded_nitrogen = kInvalidAtomNumber;
  atom_number_t doubly_bonded_oxygen = kInvalidAtomNumber;
  for (const Bond* b : a) {
    atom_number_t o = b->other(zatom);
    if (b->is_single_bond()) {
      if (m.atomic_number(o) == 7) {
        singly_bonded_nitrogen = o;
      }
    } else if (b->is_double_bond()) {
      if (m.atomic_number(o) == 8) {
        doubly_bonded_oxygen = o;
      }
    }
  }

  if (singly_bonded_nitrogen == kInvalidAtomNumber) {
    return 0;
  }
  if (doubly_bonded_oxygen == kInvalidAtomNumber) {
    return 0;
  }

  return 1;
}

// There are several possibilities.
// Molecule is the carbon of an amide, fragment is a nitrogen in an amide.
// Molecule is the nitrogen in an amide, and fragment is the carbon in an amide.
int
AdjacentAmides(Molecule& m,
                 const MoleculeData& mdata,
                 atom_number_t zatom,
                 Fragment& frag) {
  if (frag.carbon_of_amide() && m.atomic_number(zatom) == 7) {
    return 1;
  } else if (frag.nitrogen_of_amide() && IsCarbonOfAmide(m, zatom)) {
    return 1;
  }

  return 0;
}

// If zatom is doubly bonded to another carbon atom and frag
// is a nitrogen.
// Or the other way around
int
WouldFormEnamine(Molecule& m,
                 const MoleculeData& mdata,
                 atom_number_t zatom,
                 Fragment& frag) {
  if (m.atomic_number(zatom) == 6) {
    if (m.ring_bond_count(zatom) == 0 &&
        frag.atomic_number() == 7 &&
        IsEnamineCarbon(m, zatom)) {
      return 1;
    }
    return 0;
  }

  if (m.atomic_number(zatom) == 7 && frag.atomic_number() == 6) {
    if (frag.saturated() || frag.aromatic()) {
      return 0;
    }
    return IsEnamineCarbon(*frag.frag(), frag.zatom());
  }

  return 0;
}

// Return true if there is a positive charge in both `mdata` and `frag`.
int
MultiplePositiveCharges(const Molecule& m,
                        const MoleculeData& mdata,
                        const Fragment& frag) {
  if (mdata.positive == 0) {
    return 0;
  }
  if (frag.positive_charge() == 0) {
    return 0;
  }

  return 1;
}

int
Options::OkToBeCombined(Molecule& m,
                const MoleculeData& mdata,
                atom_number_t zatom,
                Fragment& frag) {
  if (_config.has_max_atoms_in_product()) {
    uint32_t atoms_in_product = m.natoms() + frag.frag()->natoms();
    if (atoms_in_product > _config.max_atoms_in_product()) {
      ++_product_too_many_atoms;
      return 0;
    }
  }
  if (_config.has_max_rings_in_product()) {
    uint32_t rings_in_product = m.nrings() + frag.frag()->nrings();
    if (rings_in_product > _config.max_rings_in_product()) {
      ++_product_too_many_rings;
      return 0;
    }
  }
  if (_config.has_max_halogens_in_product()) {
    uint32_t halogens_in_product = mdata.halogens + frag.halogens();
    if (halogens_in_product > _config.max_halogens_in_product()) {
      ++_product_too_many_halogens;
      return 0;
    }
  }
  if (_config.has_max_rotbonds_in_product()) {
    uint32_t rotbonds_in_product = mdata.rotbonds + frag.rotbonds();
    if (rotbonds_in_product > _config.max_rotbonds_in_product()) {
      ++_product_too_many_rotatable_bonds;
      return 0;
    }
  }

  // Rotatable bonds not handled yet.

  if (WouldFormBiphenyl(m, zatom, frag)) {
    return 0;
  }

  if (WouldContainTooManyNitros(m, mdata, frag)) {
    return 0;
  }

  // a more stringent test than biphenyl - if we are
  // using this, we could omit the biphenyl test.
  if (WouldFormAdjacentRings(m, zatom, frag)) {
    return 0;
  }

  if (IsNCN(m, zatom, frag)) {
    return 0;
  }

  if (WouldFormEnamine(m, mdata, zatom, frag)) {
    return 0;
  }

  if (AdjacentAmides(m, mdata, zatom, frag)) {
    return 0;
  }

  if (WouldFormAnilineNh(m, mdata, zatom, frag)) {
    return 0;
  }

  if (MultiplePositiveCharges(m, mdata, frag)) {
    return 0;
  }

  if (WouldFormPhenolicEster(m, mdata, zatom, frag)) {
    return 0;
  }

  return 1;
}

int
Options::ContainsRejectedSubstructure(Molecule& m) {
  Molecule_to_Match target(&m);

  for (Substructure_Query* q : _discard_products_containing) {
    ++_product_discarded_by_query;
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

int
Options::Seen(Molecule& m) {
  if (_seen.contains(m.unique_smiles())) {
    return 1;
  }

  _seen.insert(m.unique_smiles());

  return 0;
}

int
Options::ReadFragmentsFromTextProto(IWString& fname, int only_aromatic) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadFragmentsFromTextproto:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragmentsFromTextProto(input, only_aromatic);
}

int
Options::ReadFragmentsFromTextProto(iwstring_data_source& input, int only_aromatic) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadFragmentFromTextProto(buffer, only_aromatic)) {
      cerr << "Options::ReadFragmentsFromTextProto:cannot process " << buffer << '\n';
      return 0;
    }
  }

  return _fragment.size();
}

int
Options::ReadFragmentFromTextProto(const const_IWSubstring& buffer, int only_aromatic) {
  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());

  dicer_data::DicerFragment proto;
  if (! google::protobuf::TextFormat::Parse(&input, &proto)) {
    cerr << "DicerFragmentLookupImpl::Lookup:invalid db contents '";
    cerr << buffer << '\n';
    return 0;
  }

  return ReadFragmentFromTextProto(proto, only_aromatic);
}

int
Options::ReadFragmentFromTextProto(const dicer_data::DicerFragment& proto, int only_aromatic) {
  std::unique_ptr<Fragment> fragment = std::make_unique<Fragment>();

  Molecule* f = new Molecule();
  if (! f->build_from_smiles(proto.smi())) {
    cerr << "Options::ReadFragmentFromTextProto:invalid smiles\n";
    cerr << proto.ShortDebugString() << '\n';
    delete f;
    return 0;
  }

  f->set_name(proto.par());

  if (! fragment->SetFragment(f)) {
    delete f;
    return 0;
  }

  if (! OkFragment(*f)) {
    return 0;
  }

  if (_charge_assigner.active()) {
    fragment->DetermineFormalCharges(_charge_assigner);
  }

  fragment->set_count(proto.n());
  fragment->set_aromatic_only(only_aromatic);

  if (fragment->IsQuaternaryAryl()) {
    // Non fatal error.
    return 1;
  }

  _fragment << fragment.release();
  
  return 1;
}

int
Options::OkFragment(Molecule& frag) const {
  const uint32_t matoms = frag.natoms();
  if (matoms < _config.min_atoms_in_fragment()) {
    return 0;
  }

  if (! _config.has_max_atoms_in_fragment()) {
  } else if (matoms > _config.max_atoms_in_fragment()) {
    return 0;
  }

  return 1;
}

int
Options::ReadFragmentsAsSmiles(IWString& fname, int only_aromatic) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadFragmentsAsSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragmentsAsSmiles(input, only_aromatic);
}

int
Options::ReadFragmentsAsSmiles(data_source_and_type<Molecule>& input, int only_aromatic) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Fragment> fragment = std::make_unique<Fragment>();
    if (! fragment->SetFragment(m)) {
      cerr << "Options::ReadFragmentsAsSmiles:cannot identify attachment in " <<
              m->smiles() << ' ' << m->name() << '\n';
      return 0;
    }

    if (! OkFragment(*m)) {
      return 0;
    }

    fragment->set_aromatic_only(only_aromatic);

    _fragment << fragment.release();
  }

  return _fragment.size();
}

// This could be expanded to include molecule_filter when/if that
// gets turned into an API.
int
Options::OkProperties(Molecule& m) {
  return 1;
}

int
SubstituentEnumeration(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
SubstituentEnumeration(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    SubstituentEnumeration(options, *m, output);
  }

  return 1;
}

int
SubstituentEnumeration(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "SubstituentEnumeration:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return SubstituentEnumeration(options, input, output);
}

int
SubstituentEnumeration(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:F:m:M:y:n:pY:C:");

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
    if (! SubstituentEnumeration(options, fname, input_type, output)) {
      cerr << "SubstituentEnumeration::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace substituent_enumeration

int
main(int argc, char ** argv) {

  int rc = substituent_enumeration::SubstituentEnumeration(argc, argv);

  return rc;
}
