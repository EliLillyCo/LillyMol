// Generate a set of ReplacementRing protos from a set of molecules.
// Protos are written as text_format since there are never that many of them.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/replacement_ring.pb.h"

namespace ring_extraction {

using std::cerr;

constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

// Used in file name construction for aromatic and aliphatic rings.
// Note that some file systems may have trouble with the presence of
// files 6a.smi and 6A.smi.
IWString arom_suffix = 'a';
IWString aliph_suffix = 'A';

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
  cerr << "Extracts rings and ring systems creating ReplacementRing protos that can be used by ring_replacement\n";
  cerr << " -S <stem>         create ring data files with <stem>\n";
  cerr << " -R <rsize>        max ring size to process (def 7)\n";
  cerr << " -Z <size>         max ring system size to process (def 3)\n";
  cerr << " -k                also generate smarts with connectivity not specified\n";
  cerr << " -a                transform within ring aliphatic double bonds to type any\n";
  cerr << " -P <atype>        label atoms by atom type of exocyclic attached atom\n";
  cerr << " -X ...            more options\n";
  cerr << " -c                remove chirality\n";
  cerr << " -g ...            chemical standardisation - enter '-g help' for info\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";
  // clang-format on

  ::exit(rc);
}

// An exocyclic doubly bonded atom is characterised by its
//  the atom number in the ring
//  the atom number of the exocyclic atom
//  the atomic number of the exocyclic atom
//  the bond type - which will always be double.
struct Exocyclic {
  atom_number_t inside_ring;
  atom_number_t outside_ring;
  atomic_number_t atomic_number;
  bond_type_t btype;
};

enum class ExocyclicStatus {
  kNotInvolved = 0,
  kIsBase = 1,
  kIsExocyclic = 2
};

struct PerMoleculeArrays {
  public:
    int* ring_sys;
    int* xref;
    int* include_atom;
    uint32_t* atype;
    ExocyclicStatus* exocyclic_status;

    // True if the atom is aromatic in the parent
    int* aromatic;

    // True it an atom in the child is aromatic
    int* aromatic_in_child;

  public:
    PerMoleculeArrays(int matoms);
    ~PerMoleculeArrays();

    uint32_t* AssignAtomTypes(int matoms) {
      atype = new uint32_t[matoms];
      std::fill_n(atype, matoms, 0);
      return atype;
    }
};

PerMoleculeArrays::PerMoleculeArrays(int matoms) {
  ring_sys = new int[matoms];
  xref = new int[matoms];
  include_atom = new_int(matoms, 1);
  atype = nullptr;
  exocyclic_status = new ExocyclicStatus[matoms];
  std::fill_n(exocyclic_status, matoms, ExocyclicStatus::kNotInvolved);
  aromatic = new int[matoms];
  aromatic_in_child = new int[matoms];
}

PerMoleculeArrays::~PerMoleculeArrays() {
  delete [] ring_sys;
  delete [] xref;
  delete [] include_atom;
  delete [] atype;
  delete [] exocyclic_status;
  delete [] aromatic;
  delete [] aromatic_in_child;
}

class ExtractRings {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    // We can optionally mark the attachment point.
    isotope_t _isotope;

    // We can also have the isotope mean the atom type of what used to be
    // attached. Note that is creates ambiguity if there were multiple atom
    // types attached to the ring at a given atom. For now we ignore those rings.
    // TODO:ianwatson revisit this maybe.
    Atom_Typing_Specification _atype;

    // We can ignore rings that are too large.
    uint32_t _max_ring_size;

    // We can ignore ring systems containing too many rings.
    uint32_t _max_ring_system_size;

    int _transform_ring_double_bonds;

    // For every ring/system type we have a mapping from unique
    // smiles to protos.
    IW_STL_Hash_Map<IWString, IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>> _ring;

    // Outputs will be written to a set of files with prefix `_stem`.
    IWString _stem;

    // As a check, once the smarts for a ring is generated, we can do a search
    // in the starting molecule.
    int _substructure_search_starting_molecule;
    int _substructure_search_failures;

    // Optionally we can generate raw rings, with no substituent information.
    int _generate_substitution_not_specified;

    Chemical_Standardisation _chemical_standardisation;

    Report_Progress _report_progress;

  // Private functions.
    int GenerateRing(Molecule& parent, Molecule& m, const IWString& usmi,
                     const IWString& exocyclic_smiles,
                     const IWString& label, PerMoleculeArrays& data, int include_d);
    int GenerateSmarts(Molecule& parent, Molecule& m, PerMoleculeArrays& data,
                       int include_d, IWString& result) const;

    int ExtendToDoublyBonded(Molecule& m, int * ring_sys,
                resizable_array_p<Exocyclic>& exocyclic);
    int ExtendToDoublyBonded(Molecule& m,
                int* ring_sys,
                ExocyclicStatus* status);

    int MakeUniqueSmiles(Molecule& m,
                 int sys,
                 PerMoleculeArrays& data,
                 IWString& usmi);

    int LabelAttachmentPoints(Molecule& parent,
                              Molecule& r,
                              int sys_num,
                              const PerMoleculeArrays& data) const;

    int LabelAttachmentPoints(Molecule& parent,
                              Molecule& r,
                              const int* ring_sys,
                              int sys_num,
                              const int* xref,
                              const uint32_t* atypes) const;
    isotope_t IsotopeOfExocyclicAtom(Molecule& m,
                atom_number_t zatom,
                const int* ring_sys,
                int sys_num,
                const uint32_t* atypes) const;

    std::optional<IWString> CanonicalRingName(Molecule& m,
                        const int * ring_sys,
                        int sys) const;

    void ChangeRingDoubleBonds(Molecule& m,
                      IWString& smt) const;
    int MaybeCheckSubstructureMatch(Molecule& m, const IWString& smt);

    int WriteRings(IWString& fname,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const;
    int WriteRings(IWString_and_File_Descriptor& output,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const;
    isotope_t IsotopeForAtom(Molecule& m, atom_number_t zatom,
                             int sys_num,
                             const PerMoleculeArrays& data) const;


  public:
    ExtractRings();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    // Extract the rings from `m` and accumulate in our internal data structures.
    int Process(Molecule& m);

    // At the completion of processing, report summary information.
    int Report(std::ostream& output) const;

    // Once data is assembled, write protos.
    int WriteRings() const;
};

ExtractRings::ExtractRings() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _isotope = 1;
  _max_ring_size = 7;
  _max_ring_system_size = 3;
  _transform_ring_double_bonds = 0;
  _substructure_search_starting_molecule = 0;
  _substructure_search_failures = 0;
  _generate_substitution_not_specified = 0;
}

void
DisplayDashXOption(std::ostream& output) {
  output << " -X sss     substructure search the starting molecule as a check\n";
  output << " -X a:A     letters assigned to aromatic and aliphatic rings\n";

  ::exit(0);
}

int
ExtractRings::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    return 0;
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Output files created with name step '" << _stem << "'\n";
    }
  }

  if (cl.option_present('R')) {
    if (! cl.value('R', _max_ring_size) || _max_ring_size < 3) {
      cerr << "The max ring size (-R) option must be a valid ring size\n";
      return 0;
    }
    if (_verbose) {
      cerr << "WIll skip rings with more than " << _max_ring_size << " atoms\n";
    }
  }

  if (cl.option_present('Z')) {
    if (! cl.value('Z', _max_ring_system_size) || _max_ring_system_size < 1) {
      cerr << "The max ring system size (-Z) option must be a valid ring system size\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only extract ring systems containing as many as "
           << _max_ring_system_size << " rings\n";
    }
  }

  if (cl.option_present('a')) {
    _transform_ring_double_bonds = 1;
    if (_verbose) {
      cerr << "Will transform within ring aliphatic double bonds to type any\n";
    }
  }

  if (cl.option_present('k')) {
    _generate_substitution_not_specified = 1;
    if (_verbose) {
      cerr << "Will also generate a smarts variant without connectivity\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! _atype.build(p)) {
      cerr << "Invalid atom type '" << p << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Atom typing initialised '" << p << "'\n";
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "sss") {
        _substructure_search_starting_molecule = 1;
        if (_verbose) {
          cerr << "Will perform a substructure search vs the starting molecule\n";
        }
      } else if (x.contains(':')) {
        if (! x.split(arom_suffix, ':', aliph_suffix)
            || arom_suffix.empty() || aliph_suffix.empty()) {
          cerr << "Invalid aromatic:aliphatic file name stub '" << x << "' must be 'a:A' form\n";
          return 1;
        }
        if (_verbose) {
          cerr << "aromatic suffix '" << arom_suffix << "' aliphatic suffix '" << aliph_suffix << "'\n";
        }
      } else if (x == "help") {
        DisplayDashXOption(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOption(cerr);
      }
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  return 1;
}

int
ExtractRings::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.empty()) {
    return 0;
  }

  // Do not process phosphorus.
  if (m.natoms(15) > 0) {
    return 0;
  }

  return 1;
}

int
ExtractRings::Report(std::ostream& output) const {
  output << "ExtractRings:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  if (_substructure_search_starting_molecule) {
    output << _substructure_search_failures << " substructure search failures\n";
  }

  return 1;
}

// For each atom, unset implicit hydrogen information.
void
UnsetImplicitHydrogenInformation(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
  }
}


// within `m`, `ring_sys` designates a set of atoms comprising a
// ring system. Extend that to any doubly bonded extensions to the
// ring system.
// Update `status` for exocyclic motifs.
int
ExtractRings::ExtendToDoublyBonded(Molecule& m,
                int* ring_sys,
                ExocyclicStatus* status) {
  // Force sssr if needed.
  m.ring_membership();

  // The atoms that get added here - the =O atoms
  Set_of_Atoms added_here;
  // For each added atom, the ring systen identifier associated with it.
  resizable_array<int> sys;

  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_sys[i] == 0) {
      continue;
    }
    if (m.ncon(i) < 3) {
      continue;
    }

    const Atom& a = m.atom(i);
    if (a.nbonds() == a.ncon()) {
      continue;
    }
    for (const Bond* b : a) {
      if (b->is_single_bond()) {
        continue;
      }
      if (b->nrings()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (ring_sys[o] == ring_sys[i]) {
        continue;
      }

      if (ring_sys[o] > 0) {
        continue;
      }
      added_here << o;
      sys << ring_sys[i];

      status[i] = ExocyclicStatus::kIsBase;
      status[o] = ExocyclicStatus::kIsExocyclic;

      ++rc;
    }
  }

  for (int i = 0; i < added_here.number_elements(); ++i) {
    atom_number_t atom = added_here[i];
    ring_sys[atom] = sys[i];
    //cerr << "Atom " << atom << " added as exocyclic, sys " << sys[i] << '\n';
  }

  return rc;
}

// For each member of `ring_sys` that is equal to `sys`, and if that atom
// is an exocyclic atom, set ring_sys to 0, and an isotopic atom to `exocyclic_smiles`, and
// place the same isotope on the ring atom to which it is attached.
void
TurnOffExocyclic(Molecule& m, int* ring_sys, int sys,
                 ExocyclicStatus* exocyclic_status,
                 IWString& exocyclic_smiles,
                 isotope_t isostart) {
  const int matoms = m.natoms();
#ifdef DEBUG_TURN_OFF_EXOCYCLIC
  cerr << "Turning off exoclic from " << m.smiles() << '\n';

  for (int i = 0; i < matoms; ++i) {
    cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " to ";
    for (const Bond* b : m[i]) {
      cerr << ' ' << b->other(i);
    }
    cerr << '\n';
  }
#endif

  for (int i = 0; i < matoms; ++i) {
    // cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " status " << '\n';
    // cerr << " status " << (exocyclic_status[i] != ExocyclicStatus::kIsExocyclic) << '\n';
    if (exocyclic_status[i] != ExocyclicStatus::kIsExocyclic) {
      continue;
    }
    // cerr << "Exocyclic atom " << i << " sys " << ring_sys[i] << " sys " << sys << '\n';
    if (ring_sys[i] != sys) {
      continue;
    }

    ring_sys[i] = 0;
    if (exocyclic_smiles.size() > 0) {
      exocyclic_smiles << '.';
    }
    exocyclic_smiles << '[' << isostart << m.atomic_symbol(i) << ']';
    // cerr << "Update exocyclic smiles to " << exocyclic_smiles << '\n';
    // Place corresponding isotope on attached atom in system as an atom map number.
    bool placed = false;
    for (const Bond* b : m[i]) {
      if (! b->is_double_bond()) {
        continue;
      }
      atom_number_t j = b->other(i);
      // This probably never happens.
      // cerr << " atom " << i << " bonded to " << j << " ring_sys " << ring_sys[j] << " cmp " << sys << '\n';
      if (ring_sys[j] != sys) {
        continue;
      }
      m.set_atom_map_number(j, isostart);
      ++isostart;
      placed = true;
      break;
    }
    if (! placed) {
      cerr << "TurnOffExocyclic:No ring sys atom for atom " << i << ' '
           << m.smarts_equivalent_for_atom(i) << ' ' << m.name() << '\n';
      return;
    }
  }
}

// Create a subset of `m` consisting of those atoms in `data.ring_sys` equal
// to `sys`. Place the unique smiles of the subset in `usmi`.
int
ExtractRings::MakeUniqueSmiles(Molecule& m,
                 int sys,
                 PerMoleculeArrays& data,
                 IWString& usmi) {
  Molecule subset;
#ifdef DEBUG_MAKE_UNIQUE_SMILES
  cerr << "system " << sys << '\n';
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << data.ring_sys[i] << '\n';
  }
#endif

  if (! m.create_subset(subset, data.ring_sys, sys, data.xref)) {
    return 0;
  }

  if (! LabelAttachmentPoints(m, subset, sys, data)) {
    return 0;
  }

  usmi = subset.unique_smiles();
  return 1;
}

// Return true if any atom that is in `sys` is joined to the
// rest of the molecule via a double bond.
int
JoinIsViaDoubleBond(Molecule& m,
                    int sys,
                    PerMoleculeArrays& data) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (data.ring_sys[i] != sys) {
      continue;
    }

    const Atom& a = m[i];
    if (a.ncon() == 2) {
      continue;
    }

    for (const Bond* b : a) {
      if (! b->is_double_bond()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (data.ring_sys[o] == sys) {
        continue;
      }

      if (m.ncon(o) > 1) {
        return 1;
      }
    }
  }

  return 0;
}

int
ExtractRings::Process(Molecule& m) {
  ++_molecules_read;

  // cerr << "Begin processing " << m.name() << " nrings " << m.nrings() << "\n";

  if (m.nrings() == 0) {
    return 1;
  }
  static std::time_t tzero = std::time(nullptr);
  if (_report_progress()) {
    std::time_t tnow = std::time(nullptr);
    cerr << "ExtractRings::Processed " << _molecules_read <<
            " molecules " << (tnow - tzero) << " seconds\n";
  }

  const int matoms = m.natoms();
  PerMoleculeArrays data(matoms);

  for (int i = 0; i < matoms; ++i) {
    if (m.is_aromatic(i)) {
      data.aromatic[i] = 1;
    } else {
      data.aromatic[i] = 1;
    }
  }

  if (_atype.active()) {
    _atype.assign_atom_types(m, data.AssignAtomTypes(matoms));
  }

  m.label_atoms_by_ring_system_including_spiro_fused(data.ring_sys);

  int exocyclic_present = ExtendToDoublyBonded(m, data.ring_sys, data.exocyclic_status);

#ifdef DEBUG_PROCESS
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << '\n';
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " i " << i << " ring_sys " << data.ring_sys[i] << '\n';
  }
#endif

  // Handling exocyclic atoms is very problematic. We form an initial unique
  // smiles of all the atoms in the subset, but then, if there are exocyclic
  // atoms, form a different subset for the smarts.

  IWString exocyclic_smiles;
  static constexpr isotope_t kIsotope = 70;

  // For each ring system identified....
  for (int sys = 1; ; ++sys) {
    // cerr << "Processing ring sys " << sys << '\n';
    if (JoinIsViaDoubleBond(m, sys, data)) {
      continue;
    }

    // The subset used for uniqueness.
    IWString usmi;
    if (!MakeUniqueSmiles(m, sys, data, usmi)) {
      break;
    }

    // We want exocyclic atoms not present for the query, but present for the smiles.
    // We remove the atoms from ring_sys, place an isotope on the attachment points,
    // and get a list of the exocyclic atoms temporarily removed.
    // exocyclic_smiles gets a list of the atoms to be attached.
    exocyclic_smiles.resize_keep_storage(0);
    if (exocyclic_present) {
      TurnOffExocyclic(m, data.ring_sys, sys, data.exocyclic_status,
                       exocyclic_smiles, kIsotope);
      --exocyclic_present;
    }

    int* xref = data.xref;
    Molecule r;
    if (! m.create_subset(r, data.ring_sys, sys, xref)) {
      break;
    }

    for (int i = 0; i < matoms; ++i) {
      if (xref[i] >= 0){
        data.aromatic_in_child[xref[i]] = data.aromatic[i];
      }
    }

    UnsetImplicitHydrogenInformation(r);
    // cerr << "Subset created " << m.smiles() << ' ' << r.smiles() << '\n';

    if (! LabelAttachmentPoints(m, r, sys, data)) {
      continue;
    }

    std::optional<IWString> label = CanonicalRingName(m, data.ring_sys, sys);
    if (! label) {
      continue;
    }
    // cerr << "Label '" << label << "' " << r.unique_smiles() << '\n';

    r.set_name(m.name());
    r << '.' << *label;

    GenerateRing(m, r, usmi, exocyclic_smiles, *label, data, 1);
    //Maybe also compute the variant with no substition information.
    // This does not work if we have added exocyclic isotopes. TODO:ianwatson fix.
    if (_generate_substitution_not_specified && exocyclic_smiles.empty()) {
      r.transform_to_non_isotopic_form();
      GenerateRing(m, r, usmi, exocyclic_smiles, *label, data, 0);
    }
  }

  return 1;
}

// Atom `zatom` is part of a ring system where ring_sys[zatom] == sys_num.
// If it is bonded to an exocyclic atom, return the `atypes` value for
// that attached atom. If there are no attached atoms outside the ring
// system, return 0;
isotope_t
ExtractRings::IsotopeOfExocyclicAtom(Molecule& m,
                atom_number_t zatom,
                const int* ring_sys,
                int sys_num,
                const uint32_t* atypes) const {
  const Atom& atom = m.atom(zatom);
  // If 2 connected, no exocyclic bonds. Check is redundant, has been
  // checked before here.
  if (atom.ncon() == 2) {
    return 0;
  }

  for (const Bond* b : atom) {
    atom_number_t j = b->other(zatom);
    if (ring_sys[j] == sys_num) {
      continue;
    }

    return atypes[j];
  }

  // All attached atoms are part of the ring system.
  return 0;
}

// We are applying some kind of isotope to an atom. If we have atom types
// use that, otherwise _isotope.
isotope_t
ExtractRings::IsotopeForAtom(Molecule& m, atom_number_t zatom, 
                             int sys_num,
                             const PerMoleculeArrays& data) const {
//cerr << "IsotopeForAtom " << atom_number << " atypes " << atypes.get() << " value " << atypes[atom_number] << '\n';

  if (data.atype) {
    return IsotopeOfExocyclicAtom(m, zatom, data.ring_sys, sys_num, data.atype);
  }
  
  return _isotope;
}

// Return 1 if there any atoms where ring_sys[i] == sys_num that correspond
// to more than chain bonds attached to the ring atom.
int
AnyMultiplyConnectedAtoms(Molecule& m,
                const int* ring_sys,
                int sys_num) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_sys[i] != sys_num) {
      continue;
    }

    if (m.is_aromatic(i)) {
      continue;
    }

    const Atom& atom = m.atom(i);
    if (atom.ncon() <= 3) {
      continue;
    }

    // There can be just 1 connection that is not in a ring.
    if (m.ring_bond_count(i) + 1 != atom.ncon()) {
      return 1;
    }
  }

  return 0;  // none detected.
}

int
ExtractRings::LabelAttachmentPoints(Molecule& parent,
                              Molecule& r,
                              int sys_num,
                              const PerMoleculeArrays& data) const {
  const int matoms = parent.natoms();

  // If we are applying isotopes for the attachment point, make sure we do not have
  // any ring atoms attached to more than 1 exocyclic chain atom.
  if (data.atype) {
    if (AnyMultiplyConnectedAtoms(parent, data.ring_sys, sys_num)) {
      return 0;
    }
  }

  // Return the number of isotopes applied.
  int rc = 0;
#ifdef DEBUG_LABEL_ATTACHMENT_POINTS
  cerr << "LabelAttachmentPoints::processing ring " << sys_num << '\n';
#endif
  for (int i = 0; i < matoms; ++i) {
    // cerr << "parent atom " << i << " ring_sys " << data.ring_sys[i] << '\n';
    if (data.ring_sys[i] != sys_num) {
      continue;
    }

    const Atom& a = parent.atom(i);
    // cerr << "  acon " << a.ncon() << '\n';
    if (a.ncon() < 3) {
      continue;
    }

#ifdef DEBUG_LABEL_ATTACHMENT_POINTS
    cerr << "exocyclic_status base " << (data.exocyclic_status[i] == ExocyclicStatus::kIsBase) << '\n';
    cerr << "exocyclic_status exo " << (data.exocyclic_status[i] == ExocyclicStatus::kIsExocyclic) << '\n';
    cerr << "exocyclic_status not " << (data.exocyclic_status[i] == ExocyclicStatus::kNotInvolved) << '\n';
#endif
    // No isotopes at exocyclic ring atoms.
    if (data.exocyclic_status[i] == ExocyclicStatus::kIsBase) {
      continue;
    }

    for (const Bond* b : a) {
      atom_number_t o = b->other(i);

      // cerr << "   attached to " << o << " ring_sys " << data.ring_sys[o] << '\n';
      if (data.ring_sys[o] == data.ring_sys[i]) {
        continue;
      }

      r.set_isotope(data.xref[i], IsotopeForAtom(parent, i, sys_num, data));
      ++rc;
      //cerr << "ring set isotope on atom " << xref[i] << " value " << _isotope << '\n';
    }
  }

  return rc;
}

// When building a canonical label, information needed about the
// component rings in a ring system.
struct RType {
  int rsize;
  int aromatic;

  RType() {
    rsize = 0;
    aromatic = 0;
  }
  RType(int s, int a) {
    rsize = s;
    aromatic = a;
  }
};

IWString
operator<< (IWString& buffer, const RType& rtype) {
  buffer << rtype.rsize;
  if (rtype.aromatic) {
    buffer << arom_suffix;
  } else {
    buffer << aliph_suffix;
  }

  return buffer;
}

// For sorting ring types.
//First on ring size then aromaticity
struct
CompareRType {
  int operator()(const RType& rt1, const RType& rt2) const {
    if (rt1.rsize < rt2.rsize) {
      return -1;
    }
    if (rt1.rsize > rt2.rsize) {
      return 1;
    }
    if (rt1.aromatic > rt2.aromatic) {
      return -1;
    }
    if (rt1.aromatic < rt2.aromatic) {
      return 1;
    }
    return 0;
  }
};

// Maybe return a canonical name for the ring system defined
// by the atoms `ring_sys[i] == sys`.
std::optional<IWString>
ExtractRings::CanonicalRingName(Molecule& m,
                        const int * ring_sys,
                        int sys) const {
  m.compute_aromaticity_if_needed();

  // Gather the rings in the system.
  std::vector<RType> rtype;
  for (const Ring* r : m.sssr_rings()) {
    if (r->count_members_set_in_array(ring_sys, sys) == 0) {
      continue;
    }
    if (r->size() > _max_ring_size) {
      return std::nullopt;
    }
    rtype.emplace_back(RType(r->number_elements(), r->is_aromatic()));
  }

  if (rtype.size() > 1) {
    if (rtype.size() > _max_ring_system_size) {
      return std::nullopt;
    }

    static CompareRType cmp;
    iwqsort(rtype.data(), rtype.size(), cmp);
  }

  IWString result;
  for (const RType& r: rtype) {
    result << r;
  }

  return result;
}

int
GetPosition(const resizable_array<atom_number_t>& order, int n) {
  return order.index(n);
}

// A smarts `smt` has been formed, if we are doing substructure checks, return
// true if `smt` matches `m`.
int
ExtractRings::MaybeCheckSubstructureMatch(Molecule& m, const IWString& smt) {
  if (! _substructure_search_starting_molecule) {
    return 1;
  }

  Substructure_Query query;
  if (! query.create_from_smarts(smt)) {
    cerr << "ExtractRings::MaybeCheckSubstructureMatch:invalid smarts '" << smt << "'\n";
    return 0;
  }

  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  if (query.substructure_search(target, sresults)) {
    return 1;
  }

  cerr << "No substructure match " << m.smiles() << " smt '" << smt << "' matched " <<
      sresults.max_query_atoms_matched_in_search() << " query atoms\n";
  cerr << "Contains " << m.aromatic_atom_count() << " aromatic atoms\n";
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << '\n';
  ++_substructure_search_failures;

  return 0;
}

// `smt` is a unique smiles for `m`. If there are ring bonds in `smt`
// that are type double, change to type any.
void
ExtractRings::ChangeRingDoubleBonds(Molecule& m,
                      IWString& smt) const {
  if (! _transform_ring_double_bonds) {
    return;
  }

  if (! smt.contains('=')) {
    return;
  }

  const resizable_array<atom_number_t> & atom_order_in_smiles = m.atom_order_in_smiles();

  IWString new_smt;
  new_smt.reserve(smt.size());

  int atom_number = -1;
  for (int i = 0; i < smt.number_elements(); ++i) {
    const char c = smt[i];
    if (c == kOpenSquareBracket) {
      ++atom_number;
      new_smt << c;
      continue;
    }
    if (smt[i] != '=') {
      new_smt << c;
      continue;
    }
    int previous_atom = atom_order_in_smiles.index(atom_number);
    if (m.ncon(previous_atom) == 1) {
      new_smt << '=';
      continue;
    }
    int next_atom = atom_order_in_smiles.index(atom_number + 1);
    // Trailing ring closures will not have a next atom
    if (next_atom < 0) {
      new_smt << '=';
      continue;
    }

    if (m.is_aromatic(previous_atom) && m.is_aromatic(next_atom)) {
      new_smt << ':';
      continue;
    }
    if (m.ncon(next_atom) == 1) {
      new_smt << '=';
      continue;
    }
    if (m.in_same_ring(previous_atom, next_atom)) {
      new_smt << '~';
    }
  }
  // cerr << "Convert " << smt << " to " << new_smt << '\n';
  smt = new_smt;
}

// Given a subset of atoms in `m` indicated by `include_atom`,
// generate a smarts.
// If `include_d` is set, each atomic smarts will include the D directive,
// either as a fixed value or as a D> directive.
// This fails for molecules like
//  O1CC[N+](=C2C=C(C)OC(=C2)C)CC1 74332-96-0
// where the exocyclic heteroatom is needed in order to make the ring
// aromatic, but that information is lost here.
int
ExtractRings::GenerateSmarts(Molecule& parent,
               Molecule& m,
               PerMoleculeArrays& data,
               int include_d,
               IWString& result) const {

  m.compute_aromaticity_if_needed();
  // cerr << "Smiles after arom " << m.aromatic_smiles() << '\n';

  const int matoms = m.natoms();

  Smiles_Information smiles_information(matoms);
  smiles_information.allocate_user_specified_atomic_smarts();
  smiles_information.set_smiles_is_smarts(1);

  int aromatic_atoms_encountered = 0;

  for (int i = 0; i < matoms; ++i) {
    // cerr << "Atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " included , arom " << m.is_aromatic(i) << '\n';

    IWString smt;
    smt << kOpenSquareBracket;
    // this is not quite correct. If the atom is exocyclic and aliphatic here
    // it prevents matching an aromatic later. Ignore for now.
    if (data.aromatic_in_child[i]) {
      smt << 'a';
      ++aromatic_atoms_encountered;
    } else {
      smt << 'A';
    }
    if (m.ring_bond_count(i)) {
      smt << 'x' << m.ring_bond_count(i);
    }
    // This is inefficient, we could precomute the string ring membership(s) for each atom.
    for (const Ring* r : m.sssr_rings()) {
      if (r->contains(i)) {
        smt << 'r' << r->size();
      }
    }

    if (! include_d) {
    } else if (_atype.active() && m.isotope(i) > 0) {
      smt << "D>" << m.ncon(i);
    }else if (_isotope == m.isotope(i)) {
      if (m.is_aromatic(i)) {
        smt << "D3";
      } else {
        smt << "D>" << m.ncon(i);
      }
    } else if (m.ncon(i) > 1) {
      smt << "D" << m.ncon(i);
    }
    smt << kCloseSquareBracket;
    smiles_information.set_user_specified_atomic_smarts(i, smt);
  }

  m.ring_membership();

  if (aromatic_atoms_encountered) {
    for (const Bond* b : m.bond_list()) {
      if (b->nrings() == 0) {
        continue;
      }

      const atom_number_t a1 = b->a1();
      const atom_number_t a2 = b->a2();
      if (data.aromatic_in_child[a1] && data.aromatic_in_child[a2]) {
        m.set_bond_permanent_aromatic(b->a1(), b->a2());
      }
    }
  }

  set_write_smiles_aromatic_bonds_as_colons(1);
  m.smiles(smiles_information);
  set_write_smiles_aromatic_bonds_as_colons(0);

  result = smiles_information.smiles();

  ChangeRingDoubleBonds(m, result);

  return 1;
}

// `m` is a subset of `parent`, governed by `include_atom`.
// `usmi` is the unique smiles of `m`, although possibly different in
// the presence of exocyclic double bonds.
// Create smiles and smarts from `m` and add to `_hash`.
// Currently exocyclic double bonds are not handled properly. Ignoring for now.
// O1CC[N+](=C2C=C(C)OC(=C2)C)CC1 74332-96-0 does not work properly...
int
ExtractRings::GenerateRing(Molecule& parent,
                           Molecule& m,
                           const IWString& usmi,
                           const IWString& exocyclic_smiles,
                           const IWString& label,
                           PerMoleculeArrays& data,
                           int include_d) {

  const IWString smi = m.smiles();

  IWString smt;
  GenerateSmarts(parent, m, data, include_d, smt);

  MaybeCheckSubstructureMatch(parent, smt);

  auto iter_label = _ring.find(label);
  if (iter_label == _ring.end()) {
    IW_STL_Hash_Map<IWString, RplRing::ReplacementRing> r;
    auto [iter, _] = _ring.emplace(std::make_pair(label, std::move(r)));
    iter_label = iter;
  }

  auto iter_usmi = iter_label->second.find(usmi);
  if (iter_usmi == iter_label->second.end()) {
    RplRing::ReplacementRing r;
    r.set_smi(smi.data(), smi.length());
    r.set_smt(smt.data(), smt.length());
    r.set_id(m.name().AsString());
    r.set_usmi(usmi.data(), usmi.length());
    r.set_conn(include_d);
    if (! exocyclic_smiles.empty()) {
      r.set_exo(exocyclic_smiles.data(), exocyclic_smiles.size());
    }
    r.set_n(1);

    iter_label->second.emplace(std::make_pair(usmi, std::move(r)));
  } else {
    const auto n = iter_usmi->second.n();
    iter_usmi->second.set_n(n + 1);
  }
  
  return 1;
}

// For each stored ring write the protos to a file
// with prefix `_stem`.
int
ExtractRings::WriteRings() const {
  for (const auto& [label, rings] : _ring) {
    IWString fname;
    fname << _stem << '_' << label << ".smi";
    if (! WriteRings(fname, rings)) {
      cerr << "ExtractRings::WriteRings:cannot write '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

// Write a set of rings `rings` to `fname` as text_proto form.
int
ExtractRings::WriteRings(IWString& fname,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "ExtractRings::WriteRings:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteRings(output, rings);
}

// Write a set of rings `rings` to `output` as text_proto form.
// Optionally sort them.
int
ExtractRings::WriteRings(IWString_and_File_Descriptor& output,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  resizable_array<const RplRing::ReplacementRing*> for_sorting;
  for_sorting.resize(rings.size());
  for (const auto& [_, r] : rings) {
    for_sorting << &r;
  }

  // Sort by number of exemplars, and secondarily by the unique smiles
  // in order to generate a canonical order, which helps unit tests.
  for_sorting.iwqsort_lambda([](const RplRing::ReplacementRing* r1,
                                const RplRing::ReplacementRing* r2) {
      if (r1->n() > r2->n()) {
        return -1;
      } else if (r1->n() < r2->n()) {
        return 1;
      }
      const std::string& usmi1 = r1->usmi();
      const std::string& usmi2 = r2->usmi();
      if (usmi1.size() < usmi2.size()) {
        return -1;
      } else if (usmi1.size() > usmi2.size()) {
        return 1;
      }

      return std::strncmp(usmi1.data(), usmi2.data(), usmi1.size());
    }
  );

  std::string buffer;

  for (const RplRing::ReplacementRing* r : for_sorting) {
    printer.PrintToString(*r, &buffer);
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }
  return 1;

  for (const auto& [usmi, r] : rings) {
    printer.PrintToString(r, &buffer);
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
ReplaceRings(ExtractRings& extract_rings,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return extract_rings.Process(m);
}

int
ReplaceRings(ExtractRings& extract_rings,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! extract_rings.Preprocess(*m)) {
      continue;
    }

    if (! ReplaceRings(extract_rings, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceRings(ExtractRings& extract_rings,
            const char * fname,
            FileType input_type,
            IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceRings(extract_rings, input, output);
}


int
ReplaceRings(ExtractRings& extract_rings,
            iwstring_data_source& input,
            FileType input_type,
            IWString_and_File_Descriptor& output) {
  IWString fname;
  while (input.next_record(fname)) {
    if (! ReplaceRings(extract_rings, fname.null_terminated_chars(), input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcS:aR:P:kr:X:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  ExtractRings extract_rings;
  if (! extract_rings.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      return 1;
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    IWString tmp(fname);
    if (tmp.starts_with("F:")) {
      tmp.remove_leading_chars(2);
      iwstring_data_source input;
      if (! input.open(tmp.null_terminated_chars())) {
        cerr << "ReplaceRings:cannot open '" << tmp << "'\n";
        return 0;
      }
      if (! ReplaceRings(extract_rings, input, input_type, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }

      continue;
    }

    if (! ReplaceRings(extract_rings, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  extract_rings.WriteRings();

  if (verbose) {
    extract_rings.Report(cerr);
  }

  return 0;
}

}  // namespace ring_extraction

int
main(int argc, char** argv) {
  int rc = ring_extraction::Main(argc, argv);

  return rc;
}
