#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>

#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwmisc/misc.h"

#include "misc2.h"
#include "path.h"
#include "substructure.h"
#include "target.h"

using std::cerr;
using std::endl;

namespace ss_ring_sys {
void
FillMatchedAtomsArray(const int* atoms_in_system, int value, const int matoms,
                      std::unique_ptr<int[]>& matched_by_global_specs) {
  if (!matched_by_global_specs) {
    matched_by_global_specs.reset(new_int(matoms));
  }
  for (int i = 0; i < matoms; ++i) {
    if (atoms_in_system[i] > 0) {
      matched_by_global_specs[i] = value;
    }
  }
}

void
SetTargetGlobalIds(const int* atoms_in_system, int value, Molecule_to_Match& target) {
  const int matoms = target.molecule()->natoms();
  for (int i = 0; i < matoms; ++i) {
    if (atoms_in_system[i]) {
      target[i].set_global_id(value);
      // cerr << "Atom " << i << " set gid " << value << '\n';
    }
  }
}

}  // namespace ss_ring_sys

Substructure_Ring_System_Specification::Substructure_Ring_System_Specification() {
  _need_per_atom_array = -1;

  _ring_systems_extend_across_spiro_fusions = false;

  _every_group_matches_atoms_in_spinach_group = false;
  _every_group_matches_length_of_spinach_group = false;

  return;
}

Substructure_Ring_System_Specification::~Substructure_Ring_System_Specification() {
  if (-5 == _match_as_match_or_rejection) {
    cerr << "Deleting already deleted Substructure_Ring_System_Specification\n";
  }

  _match_as_match_or_rejection = -5;

  return;
}

int
Substructure_Ring_System_Specification::ok() const {
  if (_match_as_match_or_rejection < 0) {
    return 0;
  }

  return 1;
}

int
Substructure_Ring_System_Specification::debug_print(std::ostream& os,
                                                    const IWString& indentation) const {
  IWString ind = indentation;
  ind += "  ";

  os << ind << "Details on Substructure_Ring_System_Specification";

  if (0 == _match_as_match_or_rejection) {
    os << ", REJECTION";
  }

  os << endl;

  if (_comment.length()) {
    os << ind << "  " << _comment << endl;
  }

  os << ind << "  Hits needed " << _hits_needed << endl;
  os << ind << "  System size " << _rings_in_system << endl;
  os << ind << "  Ncon        " << _ncon << endl;
  os << ind << "  Fused       " << _degree_of_fusion << endl;
  os << ind << "  Heteroatoms " << _heteroatom_count << endl;
  os << ind << "  Aromatic    " << _aromatic_ring_count << endl;
  os << ind << "  Non arom    " << _non_aromatic_ring_count << endl;
  os << ind << "  Max Fused N " << _largest_number_of_bonds_shared_with_another_ring
     << endl;
  os << ind << "  Strongly Fs " << _strongly_fused_ring_neighbours << endl;
  os << ind << "  Nbr Spch Gp " << _number_spinach_groups << endl;
  os << ind << "  Nbr NSPCH Gp " << _number_non_spinach_groups << endl;
  os << ind << "  Atm in Spch " << _atoms_in_spinach_group << endl;
  os << ind << "  Len of Spch " << _length_of_spinach_group << endl;
  os << ind << "  Dist to Rng " << _distance_to_another_ring << endl;
  os << ind << "  Set global id " << _set_global_id << endl;

  return os.good();
}

int
Substructure_Ring_System_Specification::terse_details(std::ostream& os,
                                                      const IWString& indentation) const {
  IWString ind = indentation;
  ind += "  ";

  os << ind << "Details on Substructure_Ring_System_Specification";

  if (0 == _match_as_match_or_rejection) {
    os << ", REJECTION";
  }

  os << endl;

  if (_hits_needed.is_set()) {
    os << ind << "  Hits needed " << _hits_needed << endl;
  }
  if (_rings_in_system.is_set()) {
    os << ind << "  Ring size   " << _rings_in_system << endl;
  }
  if (_ncon.is_set()) {
    os << ind << "  Ncon        " << _ncon << endl;
  }
  if (_degree_of_fusion.is_set()) {
    os << ind << "  Fused       " << _degree_of_fusion << endl;
  }
  if (_heteroatom_count.is_set()) {
    os << ind << "  Heteroatoms " << _heteroatom_count << endl;
  }
  if (_aromatic_ring_count.is_set()) {
    os << ind << "  Aromatic    " << _aromatic_ring_count << endl;
  }
  if (_non_aromatic_ring_count.is_set()) {
    os << ind << "  Non arom    " << _non_aromatic_ring_count << endl;
  }
  if (_number_spinach_groups.is_set()) {
    os << ind << "  number spinach groups " << _number_spinach_groups << endl;
  }
  if (_number_non_spinach_groups.is_set()) {
    os << ind << "  number non-spinach groups " << _number_non_spinach_groups << endl;
  }
  if (_atoms_in_spinach_group.is_set()) {
    os << ind << "  atoms in spinach " << _atoms_in_spinach_group << endl;
  }
  if (_length_of_spinach_group.is_set()) {
    os << ind << "  length of spinach " << _length_of_spinach_group << endl;
  }
  if (_distance_to_another_ring.is_set()) {
    os << ind << "  dist to other ring " << _distance_to_another_ring << endl;
  }
  if (_set_global_id >= 0) {
    os << ind << "  set global id " << _set_global_id << endl;
  }

  return os.good();
}

/*
  ATOMS_IN_SYSTEM is an array over target.natoms().
  It will be 1 for each atom in the ring system
*/

int
Substructure_Ring_System_Specification::_check_ncon(const int* atoms_in_system,
                                                    Molecule_to_Match& target) const {
  int atoms_in_molecule = target.natoms();

  int ncon = 0;
  int ahc = 0;  // attached heteroatoms(outside the ring system)

  for (int i = 0; i < atoms_in_molecule; i++) {
    if (0 == atoms_in_system[i]) {
      continue;
    }

    const Atom* a = target[i].atom();
    int acon = a->ncon();
    if (2 == acon) {  // can only be in 1 ring
      continue;
    }

    for (int j = 0; j < acon; j++) {
      atom_number_t k = a->other(i, j);
      if (atoms_in_system[k]) {  // atom K is in the ring system
        continue;
      }

      ncon++;
      if (!_attached_heteroatom_count.is_set()) {  // no need to determine heteroatoms
        continue;
      }

      if (6 != target[k].atomic_number()) {
        ahc++;
      }
    }
  }

  if (_ncon.is_set() && !_ncon.matches(ncon)) {
    return 0;
  }

  if (_attached_heteroatom_count.is_set() && !_attached_heteroatom_count.matches(ahc)) {
    return 0;
  }

  return 1;
}

#ifdef NOT_BEING_USED
static int
compute_atoms_with_pi_electrons(const atom_number_t* in_ring, Molecule_to_Match& target) {
  assert(nullptr != in_ring);

  int rc = 0;

  int matoms = target.natoms();
  for (int i = 0; i < matoms; i++) {
    if (0 == in_ring[i]) {
      continue;
    }

    Target_Atom& a = target[i];

    int acon = a.ncon();

    if (a.nbonds() > acon) {  // unsaturation so we assume pi electrons
      rc++;
    } else {
      Atom* a1 = const_cast<Atom*>(a.atom());

      int pi;
      if (a1->pi_electrons(pi) && pi > 0) {
        rc++;
      }
    }
  }

  // cerr << "ring system has " << rc << " atoms with pi electrons\n";

  return rc;
}
#endif

// Take a simpler view of saturation and unsaturation.
// If there is any multiple bond, an atom is considered to have pi
// electrons.
static std::tuple<int, int>
compute_pi_electrons(const int* in_ring, Molecule_to_Match& target) {
  assert(nullptr != in_ring);

  int awpe = 0;  // atoms with pi electrons
  int fsat = 0;  // fully saturated atoms

  int matoms = target.natoms();
  for (int i = 0; i < matoms; i++) {
    if (0 == in_ring[i]) {
      continue;
    }

    Target_Atom& a = target[i];

    if (a.nbonds() > a.ncon()) {  // unsaturation so we assume pi electrons
      ++awpe;
    } else {
      ++fsat;
    }
  }

  // cerr << "ring system has " << rc << " atoms with pi electrons\n";

  return std::tuple<int, int>(awpe, fsat);
}

int
Substructure_Ring_System_Specification::_check_heteroatoms(
    const int* atoms_in_system, Molecule_to_Match& target) const {
  int atoms_in_molecule = target.natoms();

  int hac = 0;
  for (int i = 0; i < atoms_in_molecule; i++) {
    if (atoms_in_system[i] && 6 != target[i].atomic_number()) {
      hac++;
    }
  }

  if (!_heteroatom_count.matches(hac)) {
    return 0;
  }

  return 1;
}

// Ring number `ring_number` is being added to a fused system.
// Update the various accumulators...
void
AddToSystem(const Ring* r, int ring_number, int* ring_already_done,
            resizable_array<int>& rings_in_system, int* atoms_in_system) {
  ring_already_done[ring_number] = 1;
  rings_in_system << ring_number;
  if (atoms_in_system != nullptr) {
    r->set_vector(atoms_in_system, 1);
  }
}

// #define DEBUG_IDENTIFY_RINGS_IN_SYSTEM

int
Substructure_Ring_System_Specification::IdentifyRingsInSystem(
    Molecule_to_Match& target, int* ring_already_done,
    resizable_array<int>& rings_in_system, int* atoms_in_system) const {
  rings_in_system.resize_keep_storage(0);
  if (atoms_in_system != nullptr) {
    std::fill_n(atoms_in_system, target.natoms(), 0);
  }

#ifdef DEBUG_IDENTIFY_RINGS_IN_SYSTEM
  cerr << "Substructure_Ring_System_Specification::IdentifyRingsInSystem:looking for "
          "systems among "
       << target.nrings() << " rings\n";
  for (int i = 0; i < target.nrings(); ++i) {
    cerr << "  ring " << i << " done " << ring_already_done[i] << '\n';
  }
#endif

  const int nrings = target.nrings();
  for (int i = 0; i < nrings; ++i) {
#ifdef DEBUG_IDENTIFY_RINGS_IN_SYSTEM
    cerr << "Check ring " << i << " (" << target.ringi(i)->size()
         << " atoms) ring_already_done " << ring_already_done[i] << '\n';
#endif
    if (ring_already_done[i]) {
      continue;
    }
    const Ring* ri = target.ringi(i);
    AddToSystem(ri, i, ring_already_done, rings_in_system, atoms_in_system);
    if (!ri->is_fused()) {
      break;
    }
    for (int j = i + 1; j < nrings; ++j) {
      if (ring_already_done[j]) {
        continue;
      }
      const Ring* rj = target.ringi(j);
      if (rj->fused_system_identifier() != ri->fused_system_identifier()) {
        continue;
      }
#ifdef DEBUG_IDENTIFY_RINGS_IN_SYSTEM
      cerr << "Ring " << j << " (" << rj->size() << "atoms) has same fsid\n";
#endif
      AddToSystem(rj, j, ring_already_done, rings_in_system, atoms_in_system);
    }
    break;
  }

#ifdef DEBUG_IDENTIFY_RINGS_IN_SYSTEM
  cerr << " Now have " << rings_in_system.size() << " rings in this system\n";
  for (int i = 0; i < target.nrings(); ++i) {
    cerr << "Ring " << i << " size " << target.ringi(i)->size() << " ring_already_done "
         << ring_already_done[i] << " fsid " << target.ringi(i)->fused_system_identifier()
         << '\n';
  }
  for (int r : rings_in_system) {
    cerr << " system contains ring number " << r << " size " << target.ringi(r)->size()
         << '\n';
  }
#endif

  if (rings_in_system.size() == 0) {
    return 0;
  }

  if (!_ring_systems_extend_across_spiro_fusions) {
    return rings_in_system.number_elements();
  }

  // Scan all not-yet done rings looking for spiro fusions.

  const int nr = target.nrings();

  assert(atoms_in_system != nullptr);

  while (true) {
    int found_ring = 0;
#ifdef DEBUG_IDENTIFY_RINGS_IN_SYSTEM
    cerr << "Checking " << nr << " rings for spiro fused, currently "
         << rings_in_system.size() << " rings in system\n";
#endif
    for (int i = 0; i < nr; ++i) {
      if (ring_already_done[i]) {
        continue;
      }
      const Ring* ri = target.ringi(i);
      // cerr << " in sys " << ri->any_members_set_in_array(atoms_in_system, 1) << '\n';
      if (!ri->any_members_set_in_array(atoms_in_system, 1)) {
        continue;
      }
      AddToSystem(ri, i, ring_already_done, rings_in_system, atoms_in_system);
      found_ring = 1;
      if (!ri->is_fused()) {
        continue;
      }
      // Capture all other rings with same fsid.
      for (int j = 0; j < nr; ++j) {
        // cerr << " inner loop j = " << j << " ring_already_done " <<
        // ring_already_done[j] << '\n';
        if (ring_already_done[j]) {
          continue;
        }
        const Ring* rj = target.ringi(j);
        if (!rj->is_fused()) {
          continue;
        }
        if (rj->fused_system_identifier() == ri->fused_system_identifier()) {
          AddToSystem(rj, j, ring_already_done, rings_in_system, atoms_in_system);
        }
      }
    }
    if (!found_ring) {
      break;
    }
  }

  return 1;
}

void
UpdateAromaticity(const Ring* r, int& non_arom, int& arom) {
  if (r->is_aromatic()) {
    ++arom;
  } else {
    ++non_arom;
  }
}

// #define DEBUG_RING_SYS_MATCHES

int
Substructure_Ring_System_Specification::_matches(
    Molecule_to_Match& target, atom_number_t* atoms_in_system,
    std::unique_ptr<int[]>& matched_by_global_specs) {
#ifdef DEBUG_RING_SYS_MATCHES
  cerr << "Substructure_Ring_System_Specification::_matches: checking " << target.nrings()
       << " rings in target\n";
#endif

  int* ring_already_done = new_int(target.nrings());
  std::unique_ptr<int[]> free_rad(ring_already_done);

  int nhits = 0;

  int largest_ring_size = target.molecule()->LargestRingSize();
  std::unique_ptr<int[]> ring_sizes_encountered =
      std::make_unique<int[]>(largest_ring_size + 1);

  // When _all_hits_in_same_fragment is set, we need to keep track of the number
  // of hits in each fragment

  extending_resizable_array<int> hits_in_fragment;

  // If we have Substituents we need an array of natoms, and an array to keep track
  // of which substituents have matched.
  std::unique_ptr<int[]> storage;
  std::unique_ptr<int[]> substituent_matched;

  if (_substituent.size() > 0) {
    storage.reset(new int[target.natoms()]);
    substituent_matched.reset(new_int(_substituent.size()));
  }

  resizable_array<int> rings_in_system;
  while (IdentifyRingsInSystem(target, ring_already_done, rings_in_system,
                               atoms_in_system)) {
#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Found ring sys with " << rings_in_system.size() << " rings";
    for (int x : rings_in_system) {
      cerr << ' ' << x;
    }
    cerr << '\n';
    for (int k = 0; k < target.natoms(); ++k) {
      cerr << " atom " << k << " atoms_in_system " << atoms_in_system[k] << '\n';
    }
#endif

    // Info about this group of rings.

    int system_rejected = 0;
    int ring_sizes_matched = 0;
    int rings_with_strongly_fused_neighbours = 0;
    iwmax<int> max_fused(0);
    iwmax<int> largest_number_of_bonds_shared_with_another_ring(0);
    iwmax<int> strongly_fused_ring_neighbours(0);
    std::fill_n(ring_sizes_encountered.get(), largest_ring_size + 1, 0);
    int arom = 0;
    int non_arom = 0;

    for (int ring_number : rings_in_system) {
      const Ring* r = target.ringi(ring_number);
#ifdef DEBUG_RING_SYS_MATCHES
      cerr << "Processing ring " << ring_number << " size " << r->size()
           << " fused = " << r->is_fused() << '\n';
#endif

      const int rsize = r->number_elements();

      ring_sizes_encountered[rsize]++;
      if (_ring_sizes.matches(rsize)) {
        ++ring_sizes_matched;
      } else {
        system_rejected = 1;
      }

      max_fused.extra(r->fused_ring_neighbours());
      largest_number_of_bonds_shared_with_another_ring.extra(
          r->largest_number_of_bonds_shared_with_another_ring());
      strongly_fused_ring_neighbours.extra(r->strongly_fused_ring_neighbours());
      rings_with_strongly_fused_neighbours += r->strongly_fused_ring_neighbours() > 0;
      UpdateAromaticity(r, non_arom, arom);
    }

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Found system with " << rings_in_system.size()
         << " rings, rejected = " << system_rejected << endl;
    cerr << "_rings_that_must_match_ring_sizes? "
         << _rings_that_must_match_ring_sizes.is_set() << " matched "
         << ring_sizes_matched << " rings\n";
    for (int i = 3; i <= largest_ring_size; ++i) {
      if (ring_sizes_encountered[i]) {
        cerr << ring_sizes_encountered[i] << " rings of size " << i << '\n';
      }
    }
#endif

    if (!_rings_that_must_match_ring_sizes.is_set()) {  // every ring must match
      if (system_rejected) {
        continue;
      }
    } else if (!_rings_that_must_match_ring_sizes.matches(ring_sizes_matched)) {
      continue;
    }

    if (_ring_size_requirement.number_elements() > 0 &&
        !_ring_size_requirements_matched(ring_sizes_encountered.get(),
                                         largest_ring_size)) {
      system_rejected = 1;
      continue;
    }

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Now checking _rings_in_system "
         << _rings_in_system.matches(rings_in_system.number_elements()) << '\n';
#endif

    if (!_rings_in_system.matches(rings_in_system.number_elements())) {
      continue;
    }

    if (!_degree_of_fusion.matches(max_fused.maxval())) {
      continue;
    }

    if (!_aromatic_ring_count.matches(arom)) {
      continue;
    }

    if (!_non_aromatic_ring_count.matches(non_arom)) {
      continue;
    }

    if (!_largest_number_of_bonds_shared_with_another_ring.matches(
            largest_number_of_bonds_shared_with_another_ring.maxval())) {
      continue;
    }

    if (!_strongly_fused_ring_neighbours.matches(
            strongly_fused_ring_neighbours.maxval())) {
      continue;
    }

    if (!_strongly_fused_ring_count.matches(rings_with_strongly_fused_neighbours)) {
      continue;
    }

    //  Looking good, count heteroatoms in this system if needed. These
    //  cases must be handled specially because they require identification
    //  of every atom in the system

    if (_heteroatom_count.is_set()) {
      if (!_check_heteroatoms(atoms_in_system, target)) {
        continue;
      }
    }

    if (_ncon.is_set() || _attached_heteroatom_count.is_set()) {
      if (!_check_ncon(atoms_in_system, target)) {
        continue;
      }
    }

    if (_atoms_in_system.is_set()) {
      int ais = count_occurrences_of_item_in_array(1, target.natoms(), atoms_in_system);
      if (!_atoms_in_system.matches(ais)) {
        continue;
      }
    }

    if (_atoms_with_pi_electrons.is_set() || _fully_saturated_atoms.is_set()) {
      const auto [awpe, fsat] = compute_pi_electrons(atoms_in_system, target);
      if (_atoms_with_pi_electrons.is_set() && !_atoms_with_pi_electrons.matches(awpe)) {
        continue;
      }
      if (_fully_saturated_atoms.is_set() && !_fully_saturated_atoms.matches(fsat)) {
        continue;
      }
    }

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Checking environment " << _environment_atom.number_elements() << '\n';
    if (atoms_in_system != nullptr) {
      for (int i = 0; i < target.natoms(); ++i) {
        cerr << i << " atoms_in_system " << atoms_in_system[i] << '\n';
      }
    }
#endif
    if (_environment_atom.number_elements()) {
      if (!_environment_matches(target, atoms_in_system)) {
        continue;
      }
    }

    if (_number_spinach_groups.is_set() || _atoms_in_spinach_group.is_set() ||
        _length_of_spinach_group.is_set() || _number_non_spinach_groups.is_set()) {
      if (!_spinach_matches(target, atoms_in_system)) {
        continue;
      }
    }

    if (_distance_to_another_ring.is_set()) {
      if (!_match_distance_to_another_ring(target, atoms_in_system)) {
        continue;
      }
    }

    if (!_substituent.empty()) {
      bool got_match = false;
      for (int i = 0; i < _substituent.number_elements(); ++i) {
        if (!_substituent[i]->Matches(target, atoms_in_system, storage.get(),
                                      matched_by_global_specs)) {
          continue;
        }
        ++substituent_matched[i];
        got_match = true;
        break;
      }
      if (!got_match) {
        continue;
      }
    }

    //  We have a ring system which matches!

    nhits++;
    if (_all_hits_in_same_fragment) {
      const Ring* r = target.ringi(rings_in_system[0]);
      atom_number_t j = r->item(0);
      Molecule* m = target.molecule();
      hits_in_fragment[m->fragment_membership(j)]++;
    }
    if (_set_global_id >= 0) {
      if (_ring_extends_to_carbonyl) {
        ExtendToCarbonyl(*target.molecule(), atoms_in_system);
      }
      ss_ring_sys::FillMatchedAtomsArray(atoms_in_system, _set_global_id, target.natoms(),
                                         matched_by_global_specs);
      ss_ring_sys::SetTargetGlobalIds(atoms_in_system, _set_global_id, target);
    }
  }

#ifdef DEBUG_RING_SYS_MATCHES
  cerr << "Found " << nhits
       << " matching fused systems, matches nhits = " << _hits_needed.matches(nhits)
       << '\n';
#endif

  if (0 == nhits) {
    return !_match_as_match_or_rejection;
  }

  if (substituent_matched) {
    for (int i = 0; i < _substituent.number_elements(); ++i) {
      if (substituent_matched[i] == 0) {
        return 0;
      }
    }
  }

  if (_all_hits_in_same_fragment) {
    for (int i = 0; i < hits_in_fragment.number_elements(); i++) {
      if (_hits_needed.matches(hits_in_fragment[i])) {
        return _match_as_match_or_rejection;
      }
    }

    return !_match_as_match_or_rejection;
  }

  if (_hits_needed.matches(nhits)) {
    return _match_as_match_or_rejection;
  }

  return !_match_as_match_or_rejection;
}

// #define DEBUG_RING_SIZE_REQUIREMENTS_MATCHED

int
Substructure_Ring_System_Specification::_ring_size_requirements_matched(
    const int* ring_sizes_encountered, int largest_ring_size) {
#ifdef DEBUG_RING_SIZE_REQUIREMENTS_MATCHED
  cerr << "Checking " << _ring_size_requirement.size() << " ring size requirements\n";
#endif

  for (const RingSizeRequirement* rsc : _ring_size_requirement) {
    if (!rsc->Matches(ring_sizes_encountered, largest_ring_size)) {
      return 0;
    }
  }

  return 1;
}

// Work out the total number of rings in range of _ring_size and see if that
// is compatible with _count.
int
RingSizeRequirement::Matches(const int* ring_sizes_encountered,
                             int largest_ring_size) const {
  int count = 0;
  for (int rsize = 3; rsize <= largest_ring_size; ++rsize) {
    if (_ring_size.matches(rsize)) {
      count += ring_sizes_encountered[rsize];
    }
  }
#ifdef DEBUG_RING_SIZE_REQUIREMENT_MATCHES
  cerr << "Count " << count << " ring sizes in range, matches " << _count.matches(count)
       << '\n';
#endif

  return _count.matches(count);
}

int
Substructure_Ring_System_Specification::matches(
    Molecule_to_Match& target, std::unique_ptr<int[]>& matched_by_global_specs) {
  const int nr = target.nrings();

  if (0 == nr) {
    return !_match_as_match_or_rejection;
  }

  if (_aromatic_ring_count.is_set() || _non_aromatic_ring_count.is_set()) {
    target.molecule()->compute_aromaticity_if_needed();
  }

  // If we will be examining the individual atoms in the system, we need an array for them

  if (_need_per_atom_array < 0) {
    if (_heteroatom_count.is_set() || _ncon.is_set() || _atoms_in_system.is_set() ||
        _environment_atom.number_elements() || _number_spinach_groups.is_set() ||
        _atoms_in_spinach_group.is_set() || _length_of_spinach_group.is_set() ||
        _distance_to_another_ring.is_set() || _number_non_spinach_groups.is_set() ||
        _ring_systems_extend_across_spiro_fusions || !_substituent.empty() ||
        _set_global_id >= 0) {
      _need_per_atom_array = 1;
    } else {
      _need_per_atom_array = 0;
    }
  }

  std::unique_ptr<int[]> atmp;
  if (_need_per_atom_array) {
    atmp.reset(new int[target.natoms()]);
  }

  return _matches(target, atmp.get(), matched_by_global_specs);
}

int
Substructure_Ring_System_Specification::_spinach_matches(Molecule_to_Match& target,
                                                         const int* in_system) const {
  assert(nullptr != in_system);

  const int natoms = target.natoms();

  const Molecule* m = target.molecule();

  int number_spinach_groups = 0;

  int between_ring_connections = 0;

  int got_match_to_atoms_in_group = 0;
  int got_match_to_length_of_spinach = 0;

  for (int i = 0; i < natoms; i++) {
    if (!in_system[i]) {
      continue;
    }

    //  cerr << "System contains atom " << i << '\n';

    const Atom* a = m->atomi(i);

    int acon = a->ncon();

    if (2 == acon) {  // 2 connections, definitely no spinach
      continue;
    }

    for (int j = 0; j < acon; j++) {
      const Bond* b = a->item(j);

      atom_number_t k = b->other(i);

      if (in_system[k]) {
        continue;
      }

      const int spinach = target.is_spinach(k);
      // cerr << " ATom " << k << " spinach " << spinach << '\n';

      if (spinach <= 0) {
        between_ring_connections++;
        continue;
      }

      if (b->is_double_bond() && 1 == target[k].ncon()) {  // considered part of the ring
        continue;
      }

      number_spinach_groups++;

      if (! _atoms_in_spinach_group.is_set()) {
      } else if (_atoms_in_spinach_group.matches(spinach)) {
        got_match_to_atoms_in_group++;
      } else if (_every_group_matches_atoms_in_spinach_group) {
        return 0;
      }

      // 
      if (! _length_of_spinach_group.is_set()) {
        continue;
      }

      // Evaluating the atoms in the spinach group is expensive. Only look once
      // if we can.
      if (! _every_group_matches_length_of_spinach_group && got_match_to_length_of_spinach) {
        continue;
      }
      
      if (_check_length_of_spinach(*m, in_system, i, k)) {
        ++got_match_to_length_of_spinach;
      } else if (_every_group_matches_length_of_spinach_group) {
        return 0;
      }
    }
  }

  // cerr << "Checking spinach constraints, spinach " << number_spinach_groups << "
  // between ring " << between_ring_connections << '\n';
  if (_number_spinach_groups.is_set() &&
      !_number_spinach_groups.matches(number_spinach_groups)) {
    return 0;
  }

  if (_number_non_spinach_groups.is_set() &&
      !_number_non_spinach_groups.matches(between_ring_connections)) {
    return 0;
  }

  // Maybe we have the case where there

  if (!_atoms_in_spinach_group.is_set()) {  // nothing to check
    ;
  } else if (got_match_to_atoms_in_group) {
    ;
  } else if (number_spinach_groups > 0) {
    return 0;
  } else if (_atoms_in_spinach_group.number_elements() > 0) {
    return 0;
  } else  // if the only specification is a max, and we never got anything, that's OK
  {
    int notused;
    if (_atoms_in_spinach_group.min(
            notused)) {  // if min was set, and we didn't test anything, that's a fail
      return 0;
    }
  }

  if (_length_of_spinach_group.is_set() && !got_match_to_length_of_spinach) {
    return 0;
  }

  // cerr << "Spinach specification returning match\n";

  return 1;
}

static int
max_length_of_spinach(const Molecule& m, int* already_done, atom_number_t zatom) {
  already_done[zatom] = 1;

  const Atom* a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(zatom, i);

    if (already_done[j]) {
      continue;
    }

    int tmp = max_length_of_spinach(m, already_done, j);

    if (tmp > rc) {
      rc = tmp;
    }
  }

  return rc + 1;
}

int
Substructure_Ring_System_Specification::_check_length_of_spinach(
    const Molecule& m, const int* in_system, atom_number_t atom_in_ring,
    atom_number_t first_spinach_atom) const {
  int* already_done = new_int(m.natoms());

  std::unique_ptr<int[]> free_already_done(already_done);

  already_done[atom_in_ring] = 1;

  int maxdist = max_length_of_spinach(m, already_done, first_spinach_atom);

  return _length_of_spinach_group.matches(maxdist);
}

// #define DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING

static int
shortest_distance_to_another_ring(Molecule_to_Match& target, int* already_done,
                                  atom_number_t zatom) {
#ifdef DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING
  cerr << "shortest_distance_to_another_ring atom " << zatom << '\n';
#endif

  already_done[zatom] = 1;

  const Molecule* m = target.molecule();

  int rc = 0;

  const Atom* a = m->atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    const atom_number_t j = a->other(zatom, i);

    if (already_done[j]) {
      continue;
    }

    if (target[j].is_ring_atom()) {
      rc = 1;
      continue;
    }

    if (1 == target[j].ncon()) {
      continue;
    }

    if (target.is_spinach(j)) {
      continue;
    }

    int tmp = shortest_distance_to_another_ring(target, already_done, j);

    if (0 == rc) {
      rc = tmp;
    } else if (tmp < rc) {
      rc = tmp;
    }
  }

  return 1 + rc;
}

// #define DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING

int
Substructure_Ring_System_Specification::_match_distance_to_another_ring(
    Molecule_to_Match& target, const int* in_ring_system) const {
  const int matoms = target.natoms();

  int* already_done = new int[matoms];
  std::unique_ptr<int[]> free_already_done(already_done);

  const Molecule* m = target.molecule();

  int rc = 0;

#ifdef DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING
  cerr << "_match_distance_to_another_ring:checking " << matoms << " atoms\n";
#endif

  for (int i = 0; i < matoms; i++) {
    if (!in_ring_system[i]) {
      continue;
    }

    const Atom* a = m->atomi(i);

    int acon = a->ncon();

    if (2 == acon) {  // no branches outside the ring
      continue;
    }

    for (int j = 0; j < acon; j++) {
      const Bond* b = a->item(j);

      const atom_number_t k = b->other(i);

      if (in_ring_system[k]) {
        continue;
      }

      int d;
      if (0 == b->nrings() && target[k].nrings()) {
        d = 1;
      } else {
        if (!target.is_between_rings(k)) {
          continue;
        }

        set_vector(already_done, matoms, 0);

        already_done[i] = 1;

        d = shortest_distance_to_another_ring(target, already_done, k);
      }

#ifdef DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING
      cerr << "From atom " << i << " shortest dist to another ring " << d << '\n';
#endif

      const auto mm = _distance_to_another_ring.matches(d);

      if (_distance_to_another_ring
              .empty())  // we have a min and/or max specification only
      {
        if (!mm) {  // did not match, violated constraint, we are done
          return 0;
        }

        rc = 1;         // we found a ring that did NOT violate the constraint
      } else if (mm) {  // found a specific distance that matches
        return 1;
      }

      if (3 == acon) {  // just one bond outside the ring, hopefully the most common case
        break;
      }
    }
  }

  return rc;
}
