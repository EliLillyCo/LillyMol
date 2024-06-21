#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_set>

#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/path.h"

#include "Molecule_Tools/scaffolds.h"

#include "Molecule_Tools/scaffolds.pb.h"

namespace scaffolds {

using std::cerr;

// Used during region assignments.
constexpr int kLinker = -1;
constexpr int kSpinach = -2;

// When two rings are adjacent, like biphenyl, there is
// an empty linker between them
constexpr int kEmptyLinker = -1;

// `in_sys` identifies atons in the different ring systems in `m`.
// Identify =* groups and add them to `in_sys`.
int
ExtendToSinglyAttachedDoublyBonded(Molecule& m,
                                   int* in_sys) {
  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    // cerr << " atom " << i << " " << m.smarts_equivalent_for_atom(i) << " in_sys " << in_sys[i] << '\n';
    if (in_sys[i]) {
      continue;;
    }

    const Atom& a = m.atom(i);
    if (a.ncon() != 1) {
      continue;
    }

    const Bond* b = a[0];
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(i);
    // cerr << "  atom " << j << " at other end of double bond\n";
    if (in_sys[j] == 0) {
      continue;
    }

    in_sys[i] = in_sys[j];
    ++rc;
  }

  return rc;
}

int
MarkAtoms(Molecule& m,
        atom_number_t zatom,
        int* visited,
        int* in_sys) {
  visited[zatom] = 1;
  in_sys[zatom] = 1;
  int rc = 1;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (visited[o] || in_sys[o]) {
      continue;
    }

    rc += MarkAtoms(m, o, visited, in_sys);
  }

  return rc;
}

// `in_set` has all the scaffold atoms marked.
// For any atom that joins the scaffold at a ring atom, add
// that atom to `in_set`.
// Start by looping over all atoms in the scaffold and starting
// a search at all >2 connected ring atoms.
int
AddBackFromRing(Molecule& m, int* in_set) {
  int rc = 0;
  const int matoms = m.natoms();
  std::unique_ptr<int[]> visited(new_int(matoms));
  // Loop over atoms in the scaffold.
  for (int i = 0; i < matoms; ++i) {
    if (in_set[i] == 0) {
      continue;
    }

    if (m.ring_bond_count(i) == 0) {
      continue;
    }
    // We have a ring atom that is in the scaffold.
    const Atom& a = m[i];
    if (a.ncon() == 2) {
      continue;
    }
    visited[i] = 1;
    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      // Ignore other members of the scaffold.
      if (in_set[o] == 1) {
        continue;
      }
      rc += MarkAtoms(m, o, visited.get(), in_set);
    }
  }

  return rc;
}

// Exactly the same as UnsetFromRing except for the check on ring
// membership. If this ever gets changed, it would make sense to combine them
// with a functor that did the comparison with m.ring_bond_count(i).
// For now, we avoid that complexity with code duplication.

int
AddBackFromLinker(Molecule& m, int* in_set) {
  int rc = 0;
  const int matoms = m.natoms();
  std::unique_ptr<int[]> visited(new_int(matoms));
  // Loop over atoms in the scaffold.
  for (int i = 0; i < matoms; ++i) {
    if (in_set[i] == 0) {
      continue;
    }

    // Only difference from previous function.
    if (m.ring_bond_count(i)) {
      continue;
    }
    // We have a ring atom that is in a ring and in the scaffold.
    const Atom& a = m[i];
    if (a.ncon() == 2) {
      continue;
    }
    visited[i] = 1;
    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      // Ignore other members of the scaffold.
      if (in_set[o] == 1) {
        continue;
      }
      rc += MarkAtoms(m, o, visited.get(), in_set);
    }
  }

  return rc;
}

// Update `spinach` to include any singly connected, doubly bonded atoms.
int
AddDoublyBonded(const Molecule& m,
                int* spinach) {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (spinach[i]) {
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
      if (m.ncon(o) == 1) {
        spinach[o] = 0;
        ++rc;
      }
    }
  }

  return rc;
}

PerMoleculeData::PerMoleculeData(Molecule& m, int classify_spinach) {
  _matoms = m.natoms();
  _ring_sys = new int[_matoms];
  _in_system = new int[_matoms];
  _region = new int[_matoms];

  _nsys = m.label_atoms_by_ring_system_including_spiro_fused(_ring_sys);

  // This array gets indexed by the region number, which starts at _nsys + 1
  // 3*_nsys is more than what is needed...
  _atoms_in_region = new_int(_nsys + _nsys + _nsys);

  ExtendToSinglyAttachedDoublyBonded(m, _ring_sys);

  _atoms_in_subset = new int[_matoms];

  _tmp = new int[_matoms];

  _between = new_int((_nsys + 1) * (_nsys + 1));

  m.identify_spinach(_tmp);

  AddDoublyBonded(m, _tmp);

  ClassifyRegions(m, _tmp, classify_spinach);
}

PerMoleculeData::~PerMoleculeData() {
  delete [] _ring_sys;
  delete [] _in_system;
  delete [] _tmp;
  delete [] _region;
  delete [] _atoms_in_subset;
  delete [] _between;
  delete [] _atoms_in_region;
}

int
PerMoleculeData::Seen(Molecule& m) {
  const auto iter = _seen.find(m.unique_smiles());
  // Has been seen before
  if (iter != _seen.end()) {
    return 1;
  }

  // Never seen before.
  _seen.emplace(m.unique_smiles());

  // Do not output unique smiles.
  m.invalidate_smiles();

  return 0;
}

void
PerMoleculeData::SetBetween(int r1, int r2, int region_id) {
  // cerr << "between " << r1 << " and " << r2 << " region " << region_id << '\n';
  _between[r1 * _nsys + r2] = _between[r2 * _nsys + r1] = region_id;
}

// The purpose is to establish the _region array.
// for atoms in a ring system their _region value is their _ring_sys value.
// If `classify_spinach` is set, then we need to classify the spinach atoms.
// By default, only the rings and linkers are retained, but if parts of
// the spinach are being retained, then we need to classify those regions.
// Classify spinach atoms as either attached to a ring or attached
// to a linker.
// Note that `spinach` is actually `_tmp` so nothing here touches _tmp.
int
PerMoleculeData::ClassifyRegions(Molecule& m, const int* spinach, int classify_spinach) {
  for (int i = 0; i < _matoms; ++i) {
    if (_ring_sys[i]) {
      _region[i] = _ring_sys[i];
    } else if (spinach[i]) {
      _region[i] = kSpinach;
    } else {
      _region[i] = kLinker;
    }
  }

  // First identify the linker regions.
  int next_region = _nsys + 1;
  for (int i = 0; i < _matoms; ++i) {
    if (_ring_sys[i] == 0) {
      continue;
    }
    const Atom& a = m[i];
    if (a.ncon() == 2) {
      continue;
    }
    for (const Bond* b : a) {
      const atom_number_t o = b->other(i);
      if (_ring_sys[o] == _ring_sys[i]) {
        continue;
      }
      if (_region[o] == kSpinach) {
        continue;
      }
      if (_region[o] > 0) {
        continue;
      }
      if (_ring_sys[o] > 0 && _ring_sys[i] != _ring_sys[o]) {  // biphenyl
        // cerr << "biphenyl linkage " << i << '-' << o << '\n';
        SetBetween(_ring_sys[i], _ring_sys[o], kEmptyLinker);
        continue;
      }
      int other_end = AssignLinkerRegion(m, i, o, spinach, next_region);
      // cerr << "From " << i << " got to " << other_end << " linker\n";
      SetBetween(_ring_sys[i], other_end, next_region);
      _atoms_in_region[next_region] = std::count(_region, _region + _matoms, next_region);
      ++next_region;
    }
  }

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "After classifying linkers\n";
  for (int i = 0; i < _matoms; ++i) {
    cerr << i << " region " << _region[i] << '\n';
  }
#endif

  if (! classify_spinach) {
    return 1;
  }

  // Now that the linker regions are identified, look at sidechains - spinach.

  // Loop over the atoms in the scaffold, ring atoms and linkers
  for (int i = 0; i < _matoms; ++i) {
    if (spinach[i]) {
      continue;
    }

    const Atom& a = m[i];
    if (a.ncon() == 2) {
      continue;
    }

    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      // Skip if part of the scaffold.
      if (spinach[o] == 0) {
        continue;
      }
      // We now have an atom outside the scaffold, mark the attached atoms
      // according to whether we are in a ring or a linker between rings.
      // Things attached to a ring become part of that ring system.
      if (m.ring_bond_count(i)) {
        ExtendRingSystem(m, o, _ring_sys[i]);
      } else {
        PropagateRegion(m, o, _region[i]);
      }
    }
  }

  return 1;
}

// Mark all the atoms that comprise the linker atoms between two rings.
// `zatom` is a chain atom in a linker region.
// Update the _region array. Follow attached atoms that
// are NOT set in `spinach` till we get an atom for which _ring_sys is set.
// Return the value of _ring_sys.
int 
PerMoleculeData::AssignLinkerRegion(const Molecule& m, atom_number_t previous,
                 atom_number_t zatom,
                 const int *spinach, int flag) {
  // Allow for biphenyls on the first call.
  if (_ring_sys[zatom]) {
    return _ring_sys[zatom];
  }

  _region[zatom] = flag;

  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (o == previous) {
      continue;
    }
    if (b->is_double_bond() && m.ncon(o) == 1) {
      _region[o] = flag;
      continue;
    }
    if (spinach[o]) {
      continue;
    }
    return AssignLinkerRegion(m, zatom, o, spinach, flag);
  }

  return kInvalidAtomNumber;
}

int
PerMoleculeData::PropagateRegion(Molecule& m,
                                 atom_number_t zatom,
                                 int flag) {
  _region[zatom] = flag;

  int rc = 1;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (_ring_sys[o] || _region[o] == flag) {
      continue;
    }
    rc += PropagateRegion(m, o, flag);
  }

  return rc;
}

// We are adding a set of atoms which are to become part of a ring system.
// `flag` is the ring_system number to which we are adding these atoms.
int
PerMoleculeData::ExtendRingSystem(const Molecule& m,
                atom_number_t zatom, int flag) {
  _ring_sys[zatom] = flag;

  int rc = 1;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (_ring_sys[o]) {
      continue;
    }
    rc += ExtendRingSystem(m, o, flag);
  }

  return rc;
}

ScaffoldFinder::ScaffoldFinder() {
  _max_systems_lost = 0;
  _min_systems_in_subset = 0;

  return;
}

// Use Euler's formula to figure out if the combination of ring systems
// implied by `state` is connected or not.
int
PerMoleculeData::StateIsDisconnected(const std::vector<int>& state) const {
  assert(state.size() == static_cast<int>(_nsys));

  const uint32_t state_size = state.size();
  int nodes = 0;
  int nedges = 0;
  for (uint32_t i = 0; i < state_size; ++i) {
    if (state[i] == 0) {
      continue;
    }

    ++nodes;

    for (uint32_t j = i + 1; j < state_size; ++j) {
      if (state[j] == 0) {
        continue;
      }
      // cerr << " items " << i << " and " << j << " in stage, btw " << _between[(i + 1) * _nsys + (j + 1)] << '\n';

      if (_between[(i + 1) * _nsys + (j + 1)] != 0) {
        ++nedges;
      }
    }
  }

  int ncircuits = nedges - nodes + 1;
  // cerr << "State has " << nodes << " nodes and " << nedges << " edges, ncircuits " << ncircuits << '\n';
  if (ncircuits > 0) {
    return 1;
  }

  return 0;
}

// Given that ring systems `r1` and `r2` are in a subset, add any atoms
// that are in a linker between them.
int
PerMoleculeData::AddInterRingAtoms(int r1, int r2, int* atoms_in_subset, int flag) const {
  int region = _between[r1 * _nsys + r2];
  // cerr << "Betw " << r1 << " and " << r2 << " find " << region << '\n';
  if (region == 0) {
    return 0;
  }

  int rc = 0;
  for (int i = 0; i < _matoms; ++i) {
    if (_region[i] == region) {
      atoms_in_subset[i] = flag;
      ++rc;
    }
  }

  return rc;
}

int
PerMoleculeData::ApplyIsotopicLabels(Molecule& m, isotope_t linker,
                                     isotope_t substituent) const {
  int rc = 0;
  if (substituent) {
    rc = ApplySubstituentIsotope(m, substituent);
  }
  if (linker) {
    rc += ApplyLinkerIsotope(m, linker);
  }

  return rc;
}

int
PerMoleculeData::ApplySubstituentIsotope(Molecule& m, isotope_t substituent) const {
  m.ring_membership();

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    const Atom& a = m[i];
    if (a.ncon() < 2) {
      continue;
    }

    for (const Bond* b : a) {
      if (b->nrings()) {
        continue;
      }
      // If we have a chain bond to something marked as same ring system, that
      // is a substituent.
      atom_number_t o = b->other(i);
      if (_ring_sys[i] != _ring_sys[o]) {
        continue;
      }
      if (b->is_double_bond() && m.ncon(o) == 1) {
        continue;
      }

      m.set_isotope(i, substituent);
      ++rc;
      break;
    }
  }

  return rc;
}

int
PerMoleculeData::ApplyLinkerIsotope(Molecule& m, isotope_t linker) const {
  m.ring_membership();

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    for (const Bond* b : m[i]) {
      if (b->nrings()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (_ring_sys[i] != _ring_sys[o]) {
        m.set_isotope(i, linker);
        ++rc;
        break;
      }
    }
  }

  return rc;
}

int
ScaffoldFinder::Initialise(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return 1;
  }

  IWString fname = cl.string_value(flag);

  std::optional<scaffolds::ScaffoldsOptions> maybe_options =
       iwmisc::ReadTextProto<scaffolds::ScaffoldsOptions>(fname);
   if (! maybe_options) {
     cerr << "ScaffoldFinder::Initialise:cannot read proto '" << fname << "'\n";
     return 0;
   }

  return Initialise(*maybe_options);
}

int 
ScaffoldFinder::Initialise(const scaffolds::ScaffoldsOptions& proto) {
  _config = proto;
  // cerr << _config.ShortDebugString() << '\n';

  return 1;
}

// Return true if a molecule that starts with `in_molecule` ring systems,
// and a subset that contains `in_subset` ring systems is OK wrt any
// constraints on ring systems lost or a minimum count.
int
ScaffoldFinder::OkSubsetSize(uint32_t in_molecule, uint32_t in_subset) const {
  if (_config.max_systems_lost() > 0 &&
      (in_molecule - in_subset) > _config.max_systems_lost()) {
    return 0;
  }

  if (_config.min_systems_in_subset() > 0 && in_subset < _config.min_systems_in_subset()) {
    return 0;
  }

  if (_config.max_systems_in_subset() > 0 && in_subset > _config.max_systems_in_subset()) {
    return 0;
  }

  return 1;
}

int
IsCyclopropyl(Molecule& m, const Ring& ring, Set_of_Atoms& to_be_removed) {
  assert(ring.size() == 3);

  atom_number_t two_connected1 = kInvalidAtomNumber;
  atom_number_t two_connected2 = kInvalidAtomNumber;
  atom_number_t three_connections = kInvalidAtomNumber;
  for (const atom_number_t a : ring) {
    const Atom& atom = m[a];
    if (atom.atomic_number() != kCarbon) {
      continue;
    }

    if (atom.ncon() == 2) {
      if (two_connected1 == kInvalidAtomNumber) {
        two_connected1 = a;
      } else if (two_connected2 == kInvalidAtomNumber) {
        two_connected2 = a;
      }
    } else if (atom.ncon() == 3) {
      if (three_connections == kInvalidAtomNumber) {
        three_connections = a;
      } else {
        return false;
      }
    }
  }

  if (three_connections == kInvalidAtomNumber) {
    return 0;
  }
  if (two_connected1 == kInvalidAtomNumber || two_connected2 == kInvalidAtomNumber) {
    return 0;
  }

  to_be_removed << two_connected1 << two_connected2;

  return 1;
}

int
DiscardCycloPropylRings(Molecule& m) {
  // The atoms at the ends of the bonds that will be removed.
  Set_of_Atoms to_be_removed;

  int rc = 0;

  for (const Ring* r : m.sssr_rings()) {
    if (r->size() != 3) {
      continue;
    }
    
    rc += IsCyclopropyl(m, *r, to_be_removed);
  }
 
  if (rc == 0) {
    return 0;
  }

  int n = to_be_removed.number_elements();
  for (int i = 0; i < n; i += 2) {
    m.remove_bond_between_atoms(to_be_removed[i], to_be_removed[i + 1]);
  }

  return rc;
}

//#define DEBUG_SCAFFOLD_FINDER

int
ScaffoldFinder::MakeScaffolds(Molecule& m,
                              scaffolds::ScaffoldData& result) {
  const int nr = m.nrings();
  if (nr == 0) {
    return 1;
  }

  // By default, all non scaffold atoms are removed.
  bool remove_ring_based_non_scaffold_atoms = true;
  bool remove_linker_based_non_scaffold_atoms = true;

  if (_config.has_remove_ring_based_non_scaffold_atoms()) {
    remove_ring_based_non_scaffold_atoms = _config.remove_ring_based_non_scaffold_atoms();
  }
  if (_config.has_remove_linker_based_non_scaffold_atoms()) {
    remove_linker_based_non_scaffold_atoms = _config.remove_linker_based_non_scaffold_atoms();
  }

  if (_config.discard_cyclopropyl_ring()) {
    DiscardCycloPropylRings(m);
  }

  std::unique_ptr<int[]> in_sys = std::make_unique<int[]>(m.natoms());
  m.identify_spinach(in_sys.get());

  MaybeApplyIsotopicLabels(m, in_sys.get());

  // Invert to get the scaffold and count
  const int matoms = m.natoms();
  int atoms_in_subset = 0;
  for (int i = 0; i < matoms; ++i) {
    if (in_sys[i]) {
      in_sys[i] = 0;
    } else {
      in_sys[i] = 1;
      ++atoms_in_subset;
    }
  }
#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "MakeScaffolds begin\n";
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << "Atom " << i << " scafold " << in_sys[i] << '\n';
  }
  cerr << "remove_ring_based_non_scaffold_atoms " << remove_ring_based_non_scaffold_atoms << '\n';
  cerr << "remove_linker_based_non_scaffold_atoms " << remove_linker_based_non_scaffold_atoms << '\n';

#endif

  if (remove_ring_based_non_scaffold_atoms &&
      remove_linker_based_non_scaffold_atoms) {
    // Nothing to add back in, in_sys is already scaffold form.
  } else {
    if (! remove_ring_based_non_scaffold_atoms) {
      atoms_in_subset += AddBackFromRing(m, in_sys.get());
    } 
    if (! remove_linker_based_non_scaffold_atoms) {
      atoms_in_subset += AddBackFromLinker(m, in_sys.get());
    }
  }

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "After adding back sidechains\n";
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << " scafold " << in_sys[i] << '\n';
  }
  cerr << atoms_in_subset << " atoms_in_subset , cmp " << m.natoms() << '\n';
#endif

  if (atoms_in_subset == matoms) {
    return Process(m, result);
  }

  // Hmmm, these shortcuts do not update _generated...
  if (atoms_in_subset + ExtendToSinglyAttachedDoublyBonded(m, in_sys.get()) ==
      matoms) {
    return Process(m, result);
  }

  Molecule subset;
  m.create_subset(subset, in_sys.get(), 1);

  subset.set_name(m.name());

  int rc = Process(subset, result);

  ++_generated[result.subset_size()];

  return rc;
}

// `m` is the scaffold of an incoming molecule.
int
ScaffoldFinder::Process(Molecule& m,
                        scaffolds::ScaffoldData& result) {
  // If everything is being removed, no need to classify the scaffold atoms
  int classify_spinach = 0;
  if (_config.remove_ring_based_non_scaffold_atoms() &&
      _config.remove_linker_based_non_scaffold_atoms()) {
    classify_spinach = 0;
  } else {
    classify_spinach = 1;
  }

  PerMoleculeData pmd(m, classify_spinach);

  uint32_t nsys = pmd.nsys();
#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Scaffold molecule " << m.smiles() << " has " << m.nrings() << " rings and " << nsys << " ring systems\n";
#endif

  result.set_ring_sys(nsys);

  ++_nsys[nsys];

  // If just 1 ring system in the molecule, process it directly.
  if (nsys == 1) {
    if (! OkSubsetSize(1, 1)) {
      return 1;
    }
    scaffolds::ScaffoldSubset* s = result.add_subset();
    s->set_smi(m.smiles().AsString());
    s->set_ring_sys(1);
    return 1;
  }

  // Multiple ring systems present.

#ifdef DEBUG_SCAFFOLD_FINDER
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " sys " << pmd.ring_sys(i) << '\n';
  }
#endif

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Region assignments\n";
  for (int q = 0; q < m.natoms(); ++q) {
    cerr << " atom " << q << ' ' << m.smarts_equivalent_for_atom(q) << " region " << pmd.region()[q] << '\n';
  }
#endif

  // Set up arrays needed for combinations.
  std::vector<int> count(nsys);
  std::fill(count.begin(), count.end(), 2);
  std::vector<int> state(nsys);

  combinations::Combinations comb(count);

  while (comb.Next(state)) {
    const uint32_t nset = std::count(state.begin(), state.end(), 1);
#ifdef DEBUG_SCAFFOLD_FINDER
    for (int a: state) {
      cerr << ' ' << a;
    }
    cerr << " nset " << nset << '\n';
#endif
    if (nset == 0) {
      continue;
    }

    if (! OkSubsetSize(pmd.nsys(), nset)) {
      // cerr << "Fails subset size\n";
      continue;
    }

    if (nset < _config.min_systems_in_subset()) {
      continue;
    }
    // cerr << "Will be processed\n";

    Process(m, state, pmd, result);
  }

  return 1;
}


// Process a subset of the ring systems, as in `state`.
// `region` holds, for each non ring atom, a label that
// indicates a particular set of inter-ring atoms.
// `systems_reachable` holds, for each region, a list of
// the ring systems that can be reached from atoms in that region.
int
ScaffoldFinder::Process(const Molecule& m,
                 const std::vector<int>& state,
                 PerMoleculeData& pmd,
                 scaffolds::ScaffoldData& result) {
  assert(state.size() > 1);

  if (pmd.StateIsDisconnected(state)) {
    // cerr << "Disconnected state\n";
    return 0;
  }

  const int matoms = m.natoms();

  const int* ring_sys = pmd.ring_sys();

  int* atoms_in_subset = pmd.atoms_in_subset();
  std::fill_n(atoms_in_subset, matoms, 0);

  // Set this value in `atoms_in_subset` for the retained atoms.
  static constexpr int kFlag = 1;

  // First add the atoms in the ring systems present.
  // `state` is numbered 0,1,2, but ring_sys starts at 1 for atoms that are in a ring.
#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Just before enumeration\n";
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " ring_sys " << ring_sys[i] << " region " << pmd.region()[i] << '\n';
  }
#endif

  // First add the ring atoms.
  const int n = state.size();
  int nsys = 0;
  for (int i = 0; i < n; ++i) {
    if (state[i] == 0) {
      continue;
    }

    ++nsys;
    for (int j = 0; j < matoms; ++j) {
      if (ring_sys[j] == i + 1) {
        atoms_in_subset[j] = kFlag;
      }
    }
  }

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "At end of ring system assignments\n";
  for (int q = 0; q < m.natoms(); ++q) {
    cerr << q << " atoms_in_subset " << atoms_in_subset[q] << '\n';
  }
#endif

  // Now each region that is between any pair of ring systems included.
  for (int i = 0; i < n; ++i) {
    if (state[i] == 0) {
      continue;
    }
    for (int j = i + 1; j < n; ++j) {
      if (state[j] == 0) {
        continue;
      }

      // cerr << " state " << i  << ',' << j << " look between\n";
      int atoms_added = pmd.AddInterRingAtoms(i + 1, j + 1, atoms_in_subset, kFlag);
      // cerr << "Added " << atoms_added << " atoms as part of linker\n";
      if (! OkLinkerSize(atoms_added)) {
        return 1;
      }
    }
  }

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Subset atoms ";
  for (int i = 0; i < matoms; ++i) {
    if (atoms_in_subset[i]) {
      cerr << ' ' << i;
    }
  }
  cerr << '\n';
#endif

  Molecule subset;
  m.create_subset(subset, atoms_in_subset, kFlag);
  Molecule mcopy(m);
//  ExtendToSinglyAttachedDoublyBonded(mcopy, atoms_in_subset);
  // cerr << "Subset " << subset->smiles() << '\n';
  if (subset.number_fragments() > 1) {
    return 0;
  }

  return AddToResultsIfUnique(subset, pmd, nsys, result);
}

// Return true if `atoms_in_linker` does not violate atoms in
// linker constraints.
int
ScaffoldFinder::OkLinkerSize(uint32_t atoms_in_linker) const {
  if (_config.min_length_linker() == 0 &&
      _config.max_length_linker() == 0) {
    return 1;
  }

  if (atoms_in_linker < _config.min_length_linker()) {
    return 0;
  }

  if (_config.max_length_linker() > 0 &&
      atoms_in_linker > _config.max_length_linker()) {
    return 0;
  }

  return 1;
}

int
ScaffoldFinder::MaybeApplyIsotopicLabels(Molecule& m, const int* spinach) {
  isotope_t linker = _config.linker_isotope();
  isotope_t substituent = _config.substituent_isotope();
  if (linker == 0 && substituent == 0) {
    return 0;
  }

  if (_config.has_remove_linker_based_non_scaffold_atoms() && !
      _config.remove_linker_based_non_scaffold_atoms()) {
    linker = 0;
  }

  if (_config.has_remove_ring_based_non_scaffold_atoms() &&
      ! _config.remove_ring_based_non_scaffold_atoms()) {
    substituent = 0;
  }

  return ApplyIsotopicLabels(m, linker, substituent, spinach);
}

int
ScaffoldFinder::ApplyIsotopicLabels(Molecule& m, isotope_t linker,
                isotope_t substituent, const int* spinach) const {
  int rc = 0;
  if (linker) {
    rc += ApplyLinkerIsotope(m, linker, spinach);
  }
  if (substituent) {
    rc += ApplySubstituentIsotope(m, substituent, spinach);
  }

  return rc;
}

int
ScaffoldFinder::ApplyLinkerIsotope(Molecule& m, isotope_t linker, const int* spinach) const {
  m.ring_membership();

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    const Atom& a = m[i];
    if (a.ncon() <= 2) {
      continue;
    }

    for (const Bond* b : a) {
      if (b->nrings()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (spinach[o]) {
        continue;
      }

      m.set_isotope(i, linker);
      ++rc;
      break;
    }
  }

  return rc;
}

int
ScaffoldFinder::ApplySubstituentIsotope(Molecule& m, isotope_t substituent, const int* spinach) const {
  m.ring_membership();

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    const Atom& a = m[i];
    if (a.ncon() <= 2) {
      continue;
    }

    for (const Bond* b : m[i]) {
      if (b->nrings()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (spinach[o] == 0) {
        continue;
      }
      if (b->is_double_bond() && m.ncon(o) == 1) {
        continue;
      }

      m.set_isotope(i, substituent);
      ++rc;
      break;
    }
  }

  return rc;
}

int
ScaffoldFinder::AddToResultsIfUnique(Molecule& candidate,
                PerMoleculeData& pmd,
                int nsys,
                scaffolds::ScaffoldData& result) {
  if (pmd.Seen(candidate)) {
    return 0;
  }

  scaffolds::ScaffoldSubset* subset = result.add_subset();
  subset->set_smi(candidate.unique_smiles().AsString());
  subset->set_ring_sys(nsys);

  // cerr << "Added unique " << candidate.smiles() << " result now contains " << result.subset_size() << " items\n";

  return 1;
}

int
ScaffoldFinder::Report(std::ostream& output) const {
  for (int i = 0; i < _nsys.number_elements(); ++i) {
    if (_nsys[i]) {
      output << _nsys[i] << " molecules had " << i << " ring systems\n";
    }
  }

  for (int i = 0; i < _generated.number_elements(); ++i) {
    if (_generated[i]) {
      output << _generated[i] << " molecules generated " << i << " scaffold subsets\n";
    }
  }

  return 1;
}

}  // namespace scaffolds
