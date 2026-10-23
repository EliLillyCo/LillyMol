// IWDescr owns descriptor selection, descriptor storage, assigners, and
// calculation-only parameters. iwdescr_main owns tests, filters, descriptor
// range scaling/output policy, output modes, preprocessing, and all I/O.

#include "Molecule_Tools/iwdescr_lib.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/path_around_ring.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/partial_symmetry.h"
#include "Molecule_Tools/xlogp.h"

#include "Molecule_Tools/iwdescr_internal.h"

namespace {

using std::cerr;

// We store atom status flags for distance-matrix descriptors.
static constexpr uint32_t kDistanceMatrixIsRing = 1;
static constexpr uint32_t kDistanceMatrixIsAromatic = 2;
static constexpr uint32_t kDistanceMatrixIsHeteroatom = 4;
static constexpr uint32_t kDistanceMatrixIsONS = 8;


std::unique_ptr<uint32_t[]>
AssignDistanceMatrixAtomTypes(Molecule& m, const atomic_number_t* z) {
  const int matoms = m.natoms();

  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);

  for (int i = 0; i < matoms; ++i) {
    atype[i] = 0;
    if (m.ring_bond_count(i) > 0) {
      atype[i] |= kDistanceMatrixIsRing;
      if (m.is_aromatic(i)) {
        atype[i] |= kDistanceMatrixIsAromatic;
      }
    }
    if (z[i] == 6) {
      continue;
    }
    atype[i] |= kDistanceMatrixIsHeteroatom;
    if (z[i] == 7 || z[i] == 8 || z[i] == 16) {
      atype[i] |= kDistanceMatrixIsONS;
    }
  }

  return atype;
}

int
LongestContiguousRotatableBondGroup(Molecule& m, const int* bond_rotatable) {
  const int nedges = m.nedges();
  if (nedges == 0) {
    return 0;
  }

  std::vector<int> visited(nedges, 0);
  std::vector<int> to_process;
  to_process.reserve(nedges);

  int longest = 0;
  for (int bond_number = 0; bond_number < nedges; ++bond_number) {
    if (! bond_rotatable[bond_number] || visited[bond_number]) {
      continue;
    }

    int bonds_in_component = 0;
    to_process.resize(0);
    to_process.push_back(bond_number);
    visited[bond_number] = 1;

    while (! to_process.empty()) {
      const int current_bond_number = to_process.back();
      to_process.pop_back();
      ++bonds_in_component;

      const Bond* current_bond = m.bondi(current_bond_number);
      const atom_number_t atoms[2] = {current_bond->a1(), current_bond->a2()};
      for (atom_number_t atom : atoms) {
        const Atom& a = m.atom(atom);
        for (const Bond* b : a) {
          const int other_bond_number = b->bond_number();
          if (other_bond_number == current_bond_number ||
              ! bond_rotatable[other_bond_number] || visited[other_bond_number]) {
            continue;
          }
          visited[other_bond_number] = 1;
          to_process.push_back(other_bond_number);
        }
      }
    }

    if (bonds_in_component > longest) {
      longest = bonds_in_component;
    }
  }

  return longest;
}

int
AllFlexibleBonds(const Molecule& m, atom_number_t zatom, atom_number_t destination,
                 int dist, const int* dm) {
  const int matoms = m.natoms();

  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      return 0;
    }
    if (b->nrings()) {
      continue;
    }

    atom_number_t o = b->other(zatom);
    if (o == destination) {
      return 1;
    }

    if (dm[zatom * matoms + o] != (dist - 1)) {
      continue;
    }

    if (AllFlexibleBonds(m, o, destination, dist - 1, dm)) {
      return 1;
    }
  }

  return 0;
}

// Andrews-Craik-Martin atom contribution helpers. These are file-local because
// they do not depend on IWDescrImpl state; the member function owns descriptor
// storage and calls them in the legacy order.
float
IdentifyAcmN(const Molecule& m, int matoms, const atomic_number_t* z,
             const Atom* const* atom, int* already_done) {
  (void)m;

  int nplus = 0;
  int n = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] != 7) {
      continue;
    }

    if (atom[i]->formal_charge() == 1) {
      ++nplus;
    } else {
      ++n;
    }

    already_done[i] = 1;
  }

  return nplus * 11.5f + n * 1.2f;
}

float
IdentifyAcmC(const Molecule& m, int matoms, const atomic_number_t* z,
             const Atom* const* atom, int* already_done) {
  (void)m;

  int csp2 = 0;
  int csp3 = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] != 6) {
      continue;
    }

    const Atom* ac = atom[i];

    if (ac->ncon() < ac->nbonds()) {
      ++csp2;
    } else {
      ++csp3;
    }

    already_done[i] = 1;
  }

  return csp3 * 0.8f + csp2 * 0.7f;
}

float
IdentifyAcmHalogen(const Molecule& m, int matoms, const atomic_number_t* z,
                  const Atom* const* atom, int* already_done) {
  (void)m;
  (void)atom;

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] != 9 && z[i] != 17 && z[i] != 35 && z[i] != 53) {
      continue;
    }

    already_done[i] = 1;
    ++rc;
  }

  return 1.3f * rc;
}

float
IdentifyAcmOS(const Molecule& m, int matoms, const atomic_number_t* z,
              const Atom* const* atom, int* already_done) {
  (void)m;
  (void)atom;

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] == 8 || z[i] == 16) {
      already_done[i] = 1;
      ++rc;
    }
  }

  return 10.0f * rc;
}

float
IdentifyAcmCO(const Molecule& m, int matoms, const atomic_number_t* z,
              const Atom* const* atom, int* already_done) {
  (void)m;

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] != 8) {
      continue;
    }

    if (atom[i]->ncon() != 1) {
      continue;
    }

    const Bond* b = atom[i]->item(0);

    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t c = b->other(i);

    if (already_done[c]) {
      continue;
    }

    already_done[c] = 1;
    already_done[i] = 1;
    ++rc;
  }

  return rc * 3.5f;
}

float
IdentifyAcmOH(const Molecule& m, int matoms, const atomic_number_t* z,
              const Atom* const* atom, int* already_done) {
  (void)m;

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] != 8) {
      continue;
    }

    if (atom[i]->ncon() != 1) {
      continue;
    }

    if (const_cast<Atom*>(atom[i])->implicit_hydrogens() == 1) {
      ++rc;
    }
  }

  return rc * 2.5f;
}

float
IdentifyAcmAcids(const Molecule& m, int matoms, const atomic_number_t* z,
                 const Atom* const* atom, int* already_done) {
  (void)m;

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    if (z[i] != 8) {
      continue;
    }

    const Atom* ai = atom[i];

    if (ai->ncon() != 1) {
      continue;
    }

    const Bond* b = ai->item(0);

    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t c = b->other(i);

    if (already_done[c]) {
      continue;
    }

    if (z[c] != 6) {
      continue;
    }

    const Atom* ac = atom[c];

    if (ac->ncon() != 3) {
      continue;
    }

    atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 3; ++j) {
      const Bond* b = ac->item(j);

      if (b->is_double_bond()) {
        continue;
      }

      atom_number_t o = b->other(c);

      if (z[o] != 8) {
        continue;
      }

      if (already_done[o]) {
        continue;
      }

      if (atom[o]->ncon() != 1) {
        continue;
      }

      singly_bonded_oxygen = o;
      break;
    }

    if (singly_bonded_oxygen == INVALID_ATOM_NUMBER) {
      continue;
    }

    already_done[i] = 1;
    already_done[c] = 1;
    already_done[singly_bonded_oxygen] = 1;
    ++rc;
  }

  return rc * 8.2f;
}

float
IdentifyAcmPhosphoricAcids(const Molecule& m, int matoms, const atomic_number_t* z,
                           const Atom* const* atom, int* already_done) {
  (void)m;

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 15) {
      continue;
    }

    const Atom* ai = atom[i];

    if (ai->ncon() != 4) {
      continue;
    }

    atom_number_t o1 = INVALID_ATOM_NUMBER;
    atom_number_t o2 = INVALID_ATOM_NUMBER;
    atom_number_t o3 = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 4; ++j) {
      atom_number_t k = ai->other(i, j);

      if (z[k] != 8) {
        continue;
      }

      if (atom[k]->ncon() != 1) {
        continue;
      }

      if (o1 == INVALID_ATOM_NUMBER) {
        o1 = k;
      } else if (o2 == INVALID_ATOM_NUMBER) {
        o2 = k;
      } else if (o3 == INVALID_ATOM_NUMBER) {
        o3 = k;
      }
    }

    if (o3 == INVALID_ATOM_NUMBER) {
      continue;
    }

    already_done[i] = 1;
    already_done[o1] = 1;
    already_done[o2] = 1;
    already_done[o3] = 1;

    ++rc;
  }

  return rc * 10.0f;
}


int
LookForInterRingAtomsForSpinach(const Molecule& m, atom_number_t avoid,
                                atom_number_t zatom, int* inter_ring,
                                const int* ring_membership,
                                const Atom* const* atom) {
  const Atom* a = atom[zatom];

  if (a->ncon() == 1) {
    return 0;
  }

  int rc = 0;

  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);

    if (j == avoid) {
      continue;
    }

    if (ring_membership[j]) {
      inter_ring[zatom] = 1;
      rc = 1;
    } else if (LookForInterRingAtomsForSpinach(m, zatom, j, inter_ring,
                                               ring_membership, atom)) {
      inter_ring[zatom] = 1;
      rc = 1;
    }
  }

  return rc;
}

void
IdentifyInterRingAtomsForSpinach(const Molecule& m, int* inter_ring,
                                 const int* ring_membership,
                                 const Atom* const* atom) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (ring_membership[i] == 0) {
      continue;
    }

    const Atom* a = atom[i];
    if (a->ncon() == 2) {
      continue;
    }

    for (const Bond* b : *a) {
      if (b->nrings()) {
        continue;
      }

      atom_number_t k = b->other(i);

      if (ring_membership[k]) {
        continue;
      }

      if (inter_ring[k]) {
        continue;
      }

      LookForInterRingAtomsForSpinach(m, i, k, inter_ring, ring_membership, atom);
    }
  }
}

int
IdentifySpinachAtoms(const Molecule& m, int* spinach, const int* ring_membership,
                     const Atom* const* atom) {
  IdentifyInterRingAtomsForSpinach(m, spinach, ring_membership, atom);

  const int matoms = m.natoms();

  // Add doubly bonded, singly connected atoms to the ring/scaffold side before
  // inverting the selection.
  for (int i = 0; i < matoms; ++i) {
    if (! ring_membership[i] && ! spinach[i]) {
      continue;
    }

    const Atom* a = atom[i];
    const int acon = a->ncon();
    if (acon <= 2 || acon >= a->nbonds()) {
      continue;
    }

    for (const Bond* b : *a) {
      if (! b->is_double_bond()) {
        continue;
      }

      atom_number_t k = b->other(i);

      if (atom[k]->ncon() != 1) {
        continue;
      }

      spinach[k] = 1;
    }
  }

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (ring_membership[i] || spinach[i]) {
      spinach[i] = 0;
    } else {
      spinach[i] = 1;
      ++rc;
    }
  }

  return rc;
}

int
WithinScaffoldBranchedForSpinach(const Molecule& m, atom_number_t zatom) {
  const Atom& a = m[zatom];
  if (a.ncon() <= 2) {
    return 0;
  }

  int branches = 0;
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);
    if (b->is_double_bond() && m.ncon(j) == 1) {
      continue;
    }
    ++branches;
  }

  assert(branches >= 2);

  return branches - 2;
}

void
CountSpinachConnections(const Molecule& m, const Ring& r, const int* spinach,
                        int& spinach_connections,
                        int& non_spinach_connections) {
  spinach_connections = 0;
  non_spinach_connections = 0;

  const int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; ++i) {
    atom_number_t j = r[i];

    const Atom* a = m.atomi(j);

    if (a->ncon() == 2) {
      continue;
    }

    for (const Bond* b : *a) {
      if (b->nrings()) {
        continue;
      }

      atom_number_t k = b->other(j);

      if (spinach[k]) {
        ++spinach_connections;
      } else {
        ++non_spinach_connections;
      }
    }
  }
}

}  // namespace


class IWDescr::IWDescrImpl {
 public:
  IWDescrImpl();
  ~IWDescrImpl();

  IWDescrImpl(const IWDescrImpl&) = delete;
  IWDescrImpl& operator=(const IWDescrImpl&) = delete;

  int Initialise(Command_Line& cl);
  int InitialiseAll();
  int Process(Molecule& m, float* results);

  int number_descriptors_value() const {
    return number_descriptors;
  }

  const IWString& descriptor_name(int i) const;
  const Descriptor& descriptor_ref(int i) const;
  const Descriptor* descriptor_data() const { return descriptor; }
  Descriptor* GetDescriptor(const std::string& name) const;
  std::span<Descriptor> Descriptors() const {
    return std::span<Descriptor>(descriptor, number_descriptors);
  }

  int valid_descriptor_index(int i) const {
    return descriptor != nullptr && i >= 0 && i < number_descriptors;
  }

  DescriptorsToCompute& mutable_descriptors_to_compute() {
    return descriptors_to_compute;
  }
  const DescriptorsToCompute& descriptors_to_compute_options() const {
    return descriptors_to_compute;
  }

  int ReportDescriptorStatistics(std::ostream& output) const;
 private:

  // Per-molecule scratch arrays used by the legacy descriptor driver. Keeping
  // this as a private IWDescrImpl detail makes the Process(Molecule&, float*)
  // boundary clean while preserving the old calculation ordering.
  struct PerMoleculeData {
    explicit PerMoleculeData(Molecule& m);

    atomic_number_t* atomic_numbers() const { return z.get(); }
    int* connections() const { return ncon.get(); }
    const int* ring_membership_data() const { return ring_membership.get(); }
    atom_number_t* connection_buffer() const { return conn.get(); }
    const Atom** atoms() const { return atom.get(); }
    const int* aromaticity() const { return is_aromatic.get(); }
    int* already_done_data() const { return already_done.get(); }
    int* donor_acceptor_results_data() const { return donor_acceptor_results.get(); }

    // Some legacy calculations, notably charge assignment, can alter the
    // molecule's connection state. Refresh this array before later descriptor
    // families that depend on ncon.
    void RefreshConnections(Molecule& m);

    // Scratch visited array used by several migrated descriptor families. The
    // old implementation repeatedly allocated/reset local arrays. Keeping this
    // with the other per-molecule arrays makes the Process() data flow clearer
    // and avoids adding another loose local to ComputeTopologicalDescriptors().
    void ResetAlreadyDone();

    // Allocate/reset donor/acceptor scratch only when donor_acceptor descriptors
    // are enabled. Keeping this here avoids local allocations in
    // ComputeDonorAcceptorDescriptors() and mirrors the other per-molecule arrays.
    int EnsureDonorAcceptorResults();

    int matoms = 0;
    int maxcon = 0;

    std::unique_ptr<atomic_number_t[]> z;
    std::unique_ptr<int[]> ncon;
    std::unique_ptr<int[]> ring_membership;
    std::unique_ptr<atom_number_t[]> conn;
    std::unique_ptr<const Atom*[]> atom;
    std::unique_ptr<int[]> is_aromatic;
    std::unique_ptr<int[]> already_done;
    std::unique_ptr<int[]> donor_acceptor_results;
  };

  // Accumulators for the large first loop in legacy
  // compute_topological_descriptors().  These are deliberately grouped so the
  // loop can be migrated without adding dozens of temporary locals to
  // ComputeTopologicalDescriptors().  This is per-call scratch state, not
  // IWDescr persistent state.
  struct TopologicalDescriptorCounts {
    int heavy_atom_count = 0;
    int two_connected_chain_atom = 0;
    int rotatable_bonds = 0;
    int ring_atom_count = 0;
    int heteroatom_count = 0;
    int ring_heteroatom_count = 0;
    int fused_ring_atom_count = 0;
    int singly_connected_oxygen_or_sulphur_count = 0;
    int carboxylic_acid_count = 0;
    int aromatic_atom_count = 0;
    int aromatic_heteroatom_count = 0;
    int chain_multiple_bonds = 0;
    int ring_multiple_bonds = 0;
    int carbon_hydrogen_count = 0;
    int ch2 = 0;
    int d2sp3 = 0;
    int unsaturation = 0;
    int atoms_with_pi_electrons = 0;
    int electron_rich_sections = 0;
    int atoms_in_electron_rich_sections = 0;
    int largest_electron_rich_section = 0;
    int halogen_count = 0;
    int bigatom_count = 0;
    int aliphatic_carbon_count = 0;
    int non_ring_non_halogen_heteroatoms = 0;
    int csp3 = 0;
    int csp3_chain = 0;
    int hcount = 0;

    int connected[6] = {0, 0, 0, 0, 0, 0};

    int total_connectivity = 0;
    int total_aliphatic_connectivity = 0;
    int total_chain_connectivity = 0;
  };

 private:
  // Split Initialise into small computation-owned pieces. These are private
  // implementation hooks; they make it easier to migrate option blocks from
  // the old iwdescr(argc, argv) without moving caller-owned output logic.
  int InitialiseDescriptorSelection(Command_Line& cl);
  int InitialiseAssigners(Command_Line& cl);
  int InitialiseCalculationOptions(Command_Line& cl);
  int InitialiseDescriptorStorage(Command_Line& cl);
  int AllocateDescriptors();
  void SetDescriptorName(int descriptor_number, const char* name);
  void MarkBestFingerprint(int descriptor_number);
  int VerifyDescriptorNames() const;
  void InitialiseDescriptorDefaults();
  void ResetDescriptors();
  int ComputeDescriptorValues(Molecule& m);
  int ComputeDescriptorValues(Molecule& m, PerMoleculeData& data);

  // Core of the legacy compute_topological_descriptors() function. This is the
  // large atom/bond counting pass before the optional descriptor-family tail.
  int ComputeCoreTopologicalDescriptors(Molecule& m, PerMoleculeData& data, int* already_done);
  int ComputeOptionalTopologicalDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalChainDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalLongCarbonChainDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalSaturatedChainDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalSpecificAndHbondDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalSpecificGroupDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalSimpleHbondDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalLogPAndConjugationDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalLogPOnlyDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalConjugationDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalRingAndDistanceDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalCoreRingDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeOptionalRingSubstitutionDescriptors(Molecule& m, PerMoleculeData& data);
  int PopulateTopologicalCounts(Molecule& m, PerMoleculeData& data, int* already_done,
                                TopologicalDescriptorCounts& counts);
  int ComputeDistancesFromLongestPath(Molecule& m, const int* dm,
                                const int longest_path, int* in_path);
  int ComputeDistancesFromLongestPathInstance(Molecule& m, const int* dm,
                                                   atom_number_t a1, atom_number_t a2,
                                                   int* in_path, const int longest_path);
  int ComputeAtomDistributionAlongLongestPath(Molecule& m, const int* dm,
                                             atom_number_t a1, atom_number_t a2);
  int StoreTopologicalCounts(Molecule& m, PerMoleculeData& data,
                             const TopologicalDescriptorCounts& counts);
  int StoreCoreTopologicalCountDescriptors(Molecule& m, PerMoleculeData& data,
                                           const TopologicalDescriptorCounts& counts);
  int StoreAtomCompositionDescriptors(Molecule& m, PerMoleculeData& data,
                                      const TopologicalDescriptorCounts& counts);
  int StoreBasicAtomDescriptors(Molecule& m, PerMoleculeData& data,
                                const TopologicalDescriptorCounts& counts);
  int StoreCarbonAndHeteroatomDescriptors(Molecule& m, PerMoleculeData& data,
                                          const TopologicalDescriptorCounts& counts);
  int StoreCarbonHybridisationDescriptors(Molecule& m, PerMoleculeData& data,
                                          const TopologicalDescriptorCounts& counts);
  int StoreHeteroatomDescriptors(Molecule& m, PerMoleculeData& data,
                                 const TopologicalDescriptorCounts& counts);
  int StoreBondAndRingDescriptors(Molecule& m, PerMoleculeData& data,
                                  const TopologicalDescriptorCounts& counts);
  int StoreConnectionCountDescriptors(Molecule& m, PerMoleculeData& data,
                                      const TopologicalDescriptorCounts& counts);
  int StoreMultipleBondDescriptors(Molecule& m, PerMoleculeData& data,
                                   const TopologicalDescriptorCounts& counts);
  int StoreRotatableBondDescriptors(Molecule& m, PerMoleculeData& data,
                                    const TopologicalDescriptorCounts& counts);
  int StoreRingAndAromaticDescriptors(Molecule& m, PerMoleculeData& data,
                                      const TopologicalDescriptorCounts& counts);
  int StoreRingDescriptors(Molecule& m, PerMoleculeData& data,
                           const TopologicalDescriptorCounts& counts);
  int StoreAromaticDescriptors(Molecule& m, PerMoleculeData& data,
                               const TopologicalDescriptorCounts& counts);
  int StorePiElectronDescriptors(Molecule& m, PerMoleculeData& data,
                                 const TopologicalDescriptorCounts& counts);
  int StoreConnectivityDescriptors(Molecule& m, PerMoleculeData& data,
                                   const TopologicalDescriptorCounts& counts);
  int ComputeHalogenAttachmentDescriptors(Molecule& m, PerMoleculeData& data,
                                          int halogen_count);
  int ComputeRadhaEntropyDescriptors(Molecule& m, PerMoleculeData& data);

  // Legacy helper logic used by PopulateTopologicalCounts(). These are private
  // implementation details because they only affect descriptor values. The
  // first pass keeps conservative stubs for the helpers whose bodies will be
  // migrated from the legacy file next.
  int IsCarboxylicAcid(const Atom** atom, atom_number_t oxygen,
                       const atomic_number_t* z, const int* ncon) const;
  int ComputeElectronRichSection(Molecule& m, atom_number_t start,
                                 const Atom** atom, int* already_done) const;
  int TripleBondAtEitherEnd(Molecule& m, const Bond* b) const;
  int PartOfOtherwiseNonRotatableEntity(Molecule& m, atom_number_t a1,
                                        atom_number_t a2) const;
  int Guanidine(Molecule& m, atom_number_t doubly_bonded_nitrogen,
                int* nitrogens, const Atom** atom, const atomic_number_t* z,
                const int* ncon, const int* ring_membership) const;
  int IsSpiroFused(Molecule& m, atom_number_t a) const;
  int AtomsNotAlreadyMarked(const Ring& r, int* atom_already_done) const;
  int ComputeExocyclicBonds(Molecule& m, const Ring& r, const atomic_number_t* z,
                            const int* ncon, int& double_bond_attachments,
                            int& singly_connected_attachments,
                            int& singly_connected_heteroatoms,
                            int& singly_connected_donors) const;
  int JustTerminalGroupsOutsideRing(const Molecule& m, const Ring& r, const int* ncon,
                                    atom_number_t a) const;
  float ComputeRingIsolation(const Molecule& m, const Ring& r, const int* ncon) const;
  int HeteroatomsInRing(const Ring* r, const atomic_number_t* z) const;

  int StronglyFused(const Ring& r1, const Ring& r2, int matoms, int* tmp) const;
  int GrowFusedSystem(Molecule& m, int ring_index, int* ring_already_done, int flag,
                      int* atmp, int growing_strongly_fused_system) const;
  int JoinedRingsTooDifferentInSize(Molecule& m, const Ring& r, int* tmp) const;
  int ComputeNonPlanarFusedRings(Molecule& m, const atomic_number_t* z,
                                 int* ring_already_done, int* atmp);
  int ComputePlanarFusedRings(Molecule& m, const atomic_number_t* z,
                              int* ring_already_done, int* atmp);

  // Top-level descriptor-family hooks, listed in the same order as the legacy
  // iwdescriptors(Molecule&, ...) driver. These stubs make the migration path
  // explicit: each old file-scope helper that touches descriptor or calculation
  // state should become one of these private IWDescrImpl methods, or a helper
  // called by one of them.
  int ComputeChiralityDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeTopologicalDescriptors(Molecule& m, PerMoleculeData& data);

  // Sub-families historically called from compute_topological_descriptors().
  // Keeping these as separate hooks makes the next migration steps smaller and
  // avoids putting more local scratch arrays into PerMoleculeData prematurely.
  int ComputeLongCarbonChainDescriptors(Molecule& m, PerMoleculeData& data, int* already_done);
  int ComputeSaturatedChainDescriptors(Molecule& m, PerMoleculeData& data, int* already_done);
  int ComputeSpecificGroupDescriptors(Molecule& m, PerMoleculeData& data, int* already_done);
  int ComputeSimpleHbondDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeRingDescriptors(Molecule& m, PerMoleculeData& data);
  int compute_ring_fusion_descriptors(Molecule& m, const int* ncon,
                                const int* ring_membership);
  int ComputeBetweenRingDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeSpiroFusionDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeRameyDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeLogPDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeRingSubstitutionDescriptors(Molecule& m, PerMoleculeData& data);
  int compute_4_connected_carbon_stuff(const Molecule& m, const atomic_number_t* z,
                                 const int* ncon, const int* ring_membership);
  int ComputeExtendedConjugationDescriptors(Molecule& m, PerMoleculeData& data);

  int ComputeCrowdingDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeSpinachDescriptors(Molecule& m, PerMoleculeData& data);
//int ComputeRingChainDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeComplexityDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeNovartisPsaDescriptor(Molecule& m, PerMoleculeData& data);
  int ComputeMolarRefractivityDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeRuleOfFiveDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeDistanceMatrixDescriptors(Molecule& m, PerMoleculeData& data);
  int do_compute_mean_shell_occupancies(const Molecule& m, const int* dm);
  int compute_double_bond_substitution(Molecule& m, const atomic_number_t* z, const int* ncon);
  int do_compute_descriptors_related_to_centroid_atom(Molecule& m, const int* dm,
                                                const int longest_path,
                                                const atomic_number_t* z, const int* ncon,
                                                int* in_shell);
  int compute_polar_bond_descriptors(Molecule& m, const atomic_number_t* z, const int* ncon);
  int ComputeDistanceMatrixDescriptors(Molecule& m, PerMoleculeData& data,
                const int* dm, int* eccentricity);
  int ComputeSymmetryDescriptors(Molecule& m, PerMoleculeData& data);
  int NearSymmetricDescriptors(Molecule& m);
  int ComputeChargeDescriptors(Molecule& m, PerMoleculeData& data);
  int StoreChargeAssignerResults(Molecule& m, PerMoleculeData& data);
  int ComputeAndrewsCraikMartinDescriptors(Molecule& m, PerMoleculeData& data);
  int ComputeDonorAcceptorDescriptors(Molecule& m, PerMoleculeData& data);
  int StoreDonorAcceptorResults(Molecule& m, PerMoleculeData& data);
  int ComputeHbondSeparationDescriptors(Molecule& m, PerMoleculeData& data);

  void ZeroAllRingRelatedDescriptors();
  void ZeroRunFusionDescriptors();

  int CopyDescriptorsToResults(float* results) const;

  // IWDescr-owned computation state. Names intentionally preserve legacy
  // file-scope names during the mechanical conversion.
  Charge_Assigner charge_assigner;
  Donor_Acceptor_Assigner donor_acceptor_assigner;
  alogp::ALogP alogp_engine;
  xlogp::XLogPCalc xlogp_engine;
  nvrtspsa::NovartisPolarSurfaceArea tpsa_calculator;

  DescriptorsToCompute descriptors_to_compute;

  // Central descriptor result/metadata storage. Raw pointer is retained for
  // the first compile-oriented conversion; convert to unique_ptr once tests pass.
  Descriptor* descriptor = nullptr;

  // Keep the legacy member name during migration, but expose it through
  // number_descriptors_value() inside IWDescrImpl to avoid colliding with the
  // public IWDescr::number_descriptors() API.
  int number_descriptors = 0;

  int min_hbond_feature_separation = 0;
  int max_hbond_feature_separation = std::numeric_limits<int>::max();

  int compute_spiro_fusions = 1;
  int molecules_with_no_rings = 0;
  int fsdrng_descriptors_consider_just_two_ring_systems = 1;
  uint32_t max_difference_in_ring_size_for_strongly_fused = 0;

  Molecular_Weight_Control mwc;

  int ignore_molecules_with_no_atoms = 0;

  int shortest_internal_hydrogen_bond_separation = 3;
  int longest_internal_hydrogen_bond_separation = 6;

  // Descriptor calculation diagnostics/statistics that are coupled to the
  // descriptor engine rather than output formatting. These can later be exposed
  // through a Report(std::ostream&) method if needed.
  uint64_t molecules_processed = 0;
};

IWDescr::IWDescrImpl::IWDescrImpl() {
  // Even if alogp is not being computed, set some useful default values.
  alogp_engine.set_use_alcohol_for_acid(1);
  alogp_engine.set_rdkit_phoshoric_acid_hydrogen(1);
}

IWDescr::IWDescrImpl::~IWDescrImpl() {
  delete[] descriptor;
}

IWDescr::IWDescr() : _impl(std::make_unique<IWDescrImpl>()) {}
IWDescr::~IWDescr() = default;
IWDescr::IWDescr(IWDescr&&) noexcept = default;
IWDescr& IWDescr::operator=(IWDescr&&) noexcept = default;

int
IWDescr::Initialise(Command_Line& cl) {
  return _impl->Initialise(cl);
}

int
IWDescr::InitialiseAll() {
  return _impl->InitialiseAll();
}

int
IWDescr::Process(Molecule& m, float* results) {
  return _impl->Process(m, results);
}

int
IWDescr::number_descriptors() const {
  return _impl->number_descriptors_value();
}

DescriptorsToCompute&
IWDescr::mutable_descriptors_to_compute() {
  return _impl->mutable_descriptors_to_compute();
}

const DescriptorsToCompute&
IWDescr::descriptors_to_compute() const {
  return _impl->descriptors_to_compute_options();
}

const IWString&
IWDescr::descriptor_name(int i) const {
  return _impl->descriptor_name(i);
}

int
IWDescr::descriptor_active(int ndx) const {
  return _impl->descriptor_ref(ndx).active();
}

const Descriptor&
IWDescr::descriptor(int i) const {
  return _impl->descriptor_ref(i);
}

const Descriptor*
IWDescr::descriptor_data() const {
  return _impl->descriptor_data();
}

std::span<Descriptor>
IWDescr::Descriptors() const {
  return _impl->Descriptors();
}

int
IWDescr::valid_descriptor_index(int i) const {
  return _impl->valid_descriptor_index(i);
}

Descriptor*
IWDescr::GetDescriptor(const std::string& name) const {
  return _impl->GetDescriptor(name);
}

int
IWDescr::ReportDescriptorStatistics(std::ostream& output) const {
  return _impl->ReportDescriptorStatistics(output);
}

IWDescr::IWDescrImpl::PerMoleculeData::PerMoleculeData(Molecule& m) {
  matoms = m.natoms();

  z = std::make_unique<atomic_number_t[]>(matoms);
  m.atomic_numbers(z.get());

  ncon = std::make_unique<int[]>(matoms);
  maxcon = m.ncon(ncon.get());

  ring_membership = std::make_unique<int[]>(matoms);
  m.ring_membership_including_non_sssr_rings(ring_membership.get());

  conn = std::make_unique<atom_number_t[]>(maxcon);

  atom = std::make_unique<const Atom*[]>(matoms);
  m.atoms(atom.get());

  is_aromatic = std::make_unique<int[]>(matoms);
  m.compute_aromaticity_if_needed();
  for (int i = 0; i < matoms; ++i) {
    is_aromatic[i] = m.is_aromatic(i);
  }

  already_done = std::make_unique<int[]>(matoms);
  ResetAlreadyDone();
}

void
IWDescr::IWDescrImpl::PerMoleculeData::RefreshConnections(Molecule& m) {
  if (matoms == 0 || ncon == nullptr) {
    return;
  }

  maxcon = m.ncon(ncon.get());
}

void
IWDescr::IWDescrImpl::PerMoleculeData::ResetAlreadyDone() {
  if (matoms == 0 || already_done == nullptr) {
    return;
  }

  std::fill_n(already_done.get(), matoms, 0);
}

int
IWDescr::IWDescrImpl::PerMoleculeData::EnsureDonorAcceptorResults() {
  if (matoms == 0) {
    return 1;
  }

  if (donor_acceptor_results == nullptr) {
    donor_acceptor_results = std::make_unique<int[]>(matoms);
  }

  std::fill_n(donor_acceptor_results.get(), matoms, 0);
  return 1;
}

const IWString&
IWDescr::IWDescrImpl::descriptor_name(int i) const {
  static IWString empty;
  if (descriptor == nullptr || i < 0 || i >= number_descriptors) {
    return empty;
  }

  return descriptor[i].name();
}

const Descriptor&
IWDescr::IWDescrImpl::descriptor_ref(int i) const {
  // Invalid descriptor access should be rare. The public valid_descriptor_index()
  // method lets callers check before accessing. For now return a process-lifetime
  // empty Descriptor for invalid access so legacy output paths can be migrated
  // without introducing undefined behaviour.
  static Descriptor empty_descriptor;
  if (descriptor == nullptr || i < 0 || i >= number_descriptors) {
    return empty_descriptor;
  }

  return descriptor[i];
}

Descriptor*
IWDescr::IWDescrImpl::GetDescriptor(const std::string& name) const {
  for (int i = 0; i < number_descriptors; ++i) {
    if (descriptor[i].name() == name) {
      return descriptor + i;
    }
  }

  return nullptr;
}

int
IWDescr::IWDescrImpl::Initialise(Command_Line& cl) {
  const int verbose = cl.option_present('v');

  if (! InitialiseDescriptorSelection(cl)) {
    return 0;
  }

  if (! InitialiseAssigners(cl)) {
    return 0;
  }

  if (! InitialiseCalculationOptions(cl)) {
    return 0;
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', min_hbond_feature_separation) ||
        min_hbond_feature_separation < 1) {
      cerr << "The minimum Hydrogen bond feature separation (-b) option must be followed "
              "by a whole number\n";
      return 0;
    }

    if (cl.option_present('v')) {
      cerr << "Will perceive H-Bond features separated by more than "
           << min_hbond_feature_separation << " bonds or more\n";
    }

    descriptors_to_compute.hbond_descriptors = 1;
  }

  if (cl.option_present('S')) {
    tpsa_calculator.SetRDKitCompatibility();
    if (cl.option_present('v')) {
      cerr << "TPSA computation done with maximum RDKit compatability\n";
    }
  }

  // Note that these options may be set based on verbose, but may get 
  // reversed by the -B option.
  if (verbose) {
    tpsa_calculator.set_display_psa_unclassified_atom_messages(1);
  } else {
    tpsa_calculator.set_display_psa_unclassified_atom_messages(0);
  }

  if (verbose == 0) {
    alogp_engine.set_display_error_messages(0);
  }

  if (cl.option_present('B')) {
    const_IWSubstring b;
    for (int i = 0; cl.value('B', b, i); ++i) {
      if (b == "quiet") {
        tpsa_calculator.set_display_psa_unclassified_atom_messages(0);
        alogp_engine.set_display_error_messages(0);
        xlogp_engine.SetIssueUnclassifiedAtomMessages(false);
        if (verbose) {
          cerr << "Will not report unclassified atoms\n";
        }
      } else if (b.starts_with("mxsfsdf=")) {
        // This option is not documented, and should be deprecated.
        b.remove_leading_chars(8);
        if (!b.numeric_value(max_difference_in_ring_size_for_strongly_fused) ||
            max_difference_in_ring_size_for_strongly_fused < 1) {
          cerr << "The max difference in ring size in a strongly fused system must be "
                  "positive\n";
          return 0;
        }

        if (verbose) {
          cerr
              << "Strongly fused ring systems will only grow if the ring sizes differ by "
              << max_difference_in_ring_size_for_strongly_fused << " atoms or less\n";
        }
      }
    }
  }

  if (! InitialiseDescriptorStorage(cl)) {
    return 0;
  }

  return 1;
}

int
IWDescr::IWDescrImpl::InitialiseAll() {
  if (getenv("LILLYMOL_HOME") == nullptr) {
    cerr << "IWDescr::InitialiseAll:LILLYMOL_HOME is not defined\n";
    return 0;
  }

  descriptors_to_compute.SetAll(1);
  min_hbond_feature_separation = 5;

  set_aromatic_bonds_lose_kekule_identity(0);

  charge_assigner.set_verbose(0);
  if (! charge_assigner.BuildFromDefaultEnvs()) {
    cerr << "IWDescr::InitialiseAll:cannot initialise charge assigner\n";
    return 0;
  }

  IWString empty;
  if (! donor_acceptor_assigner.BuildFromDir(empty, 0)) {
    cerr << "IWDescr::InitialiseAll:cannot initialise donor/acceptor assigner\n";
    return 0;
  }

  tpsa_calculator.set_display_psa_unclassified_atom_messages(0);
  alogp_engine.set_display_error_messages(0);
  xlogp_engine.SetIssueUnclassifiedAtomMessages(false);

  if (! AllocateDescriptors()) {
    return 0;
  }
  InitialiseDescriptorDefaults();

  return VerifyDescriptorNames();
}

int
IWDescr::IWDescrImpl::InitialiseDescriptorSelection(Command_Line& cl) {
  return descriptors_to_compute.Initialise(cl);
}

int
IWDescr::IWDescrImpl::InitialiseAssigners(Command_Line& cl) {
  const int verbose = cl.option_present('v');

  // Not sure why this is important.
  set_aromatic_bonds_lose_kekule_identity(0);

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, (verbose > 1), 'N')) {
      return 0;
    }
  }

  if (cl.option_present('H')) {
    if (!donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose)) {
      cerr << "Cannot initialise donor/acceptor assignment object\n";
      return 0;
    }
  }

  return 1;
}

int
IWDescr::IWDescrImpl::InitialiseCalculationOptions(Command_Line& cl) {
  // Move computation-only option blocks from the old driver here, including:
  //   * alogp_engine and xlogp setup
  //   * min/max hbond separation
  //   * ring fusion/spiro settings
  //   * TPSA calculator compatibility setup
  //   * molecular weight control
  //   * ignore_molecules_with_no_atoms, if it controls calculation handling
  // Mixed options should be split: calculation pieces here, output/filter/test
  // pieces in iwdescr_main.cc.
  return 1;
}

int
IWDescr::IWDescrImpl::InitialiseDescriptorStorage(Command_Line& cl) {
  if (! AllocateDescriptors()) {
    return 0;
  }

  InitialiseDescriptorDefaults();

  return 1;
}

int
IWDescr::IWDescrImpl::AllocateDescriptors() {
  // Next compile-oriented step: move the existing allocate_descriptors() body
  // here. The descriptor index enum and kNumberDescriptors now live in
  // iwdescr_internal.h, so this method can allocate the IWDescr-owned storage
  // without relying on a file-scope descriptor pointer.
  if (descriptor != nullptr) {
    delete[] descriptor;
    descriptor = nullptr;
    number_descriptors = 0;
  }

  descriptor = new Descriptor[kNumberDescriptors];
  if (descriptor == nullptr) {
    return 0;
  }

  number_descriptors = kNumberDescriptors;

  SetDescriptorName(iwdescr_natoms, "natoms");
  SetDescriptorName(iwdescr_nrings, "nrings");
  SetDescriptorName(iwdescr_nelem, "nelem");
  SetDescriptorName(iwdescr_amw, "amw");
  if (descriptors_to_compute.ncon_descriptors) {
    SetDescriptorName(iwdescr_ncon1, "ncon1");
    SetDescriptorName(iwdescr_fncon1, "fncon1");
    SetDescriptorName(iwdescr_ncon2, "ncon2");
    SetDescriptorName(iwdescr_fncon2, "fncon2");
    SetDescriptorName(iwdescr_ncon3, "ncon3");
    SetDescriptorName(iwdescr_fncon3, "fncon3");
    SetDescriptorName(iwdescr_ncon4, "ncon4");
    SetDescriptorName(iwdescr_fncon4, "fncon4");
  }
  SetDescriptorName(iwdescr_frhc, "frhc");
  SetDescriptorName(iwdescr_mltbd, "mltbd");
  SetDescriptorName(iwdescr_fmltbd, "fmltbd");
  SetDescriptorName(iwdescr_chmltbd, "chmltbd");
  SetDescriptorName(iwdescr_fchmltbd, "fchmltbd");
  SetDescriptorName(iwdescr_rgmltbd, "rgmltbd");
  SetDescriptorName(iwdescr_frgmltbd, "frgmltbd");
  SetDescriptorName(iwdescr_dcca, "dcca");
  SetDescriptorName(iwdescr_fdcca, "fdcca");
  if (descriptors_to_compute.distance_matrix_descriptors) {
    SetDescriptorName(iwdescr_mxdst, "mxdst");
    SetDescriptorName(iwdescr_fmxdst, "fmxdst");
    SetDescriptorName(iwdescr_mxsdlp, "mxsdlp");
    SetDescriptorName(iwdescr_avsdlp, "avsdlp");
    SetDescriptorName(iwdescr_mxsdlprl, "mxsdlprl");
    SetDescriptorName(iwdescr_mdallp, "mdallp");
    SetDescriptorName(iwdescr_fmdallp, "fmdallp");
    SetDescriptorName(iwdescr_fdiffallp, "fdiffallp");
    SetDescriptorName(iwdescr_harary, "harary");
  }
  SetDescriptorName(iwdescr_rotbond, "rotbond");
  SetDescriptorName(iwdescr_frotbond, "frotbond");
  SetDescriptorName(iwdescr_ringatom, "ringatom");
  SetDescriptorName(iwdescr_rhacnt, "rhacnt");
  SetDescriptorName(iwdescr_rhaf, "rhaf");
  SetDescriptorName(iwdescr_frafus, "frafus");
  SetDescriptorName(iwdescr_rngatmf, "rngatmf");
  SetDescriptorName(iwdescr_aroma, "aroma");
  SetDescriptorName(iwdescr_aromha, "aromha");
  SetDescriptorName(iwdescr_fraromha, "fraromha");
  SetDescriptorName(iwdescr_aromdens, "aromdens");
  SetDescriptorName(iwdescr_ch2, "ch2");
  SetDescriptorName(iwdescr_d2sp3, "d2sp3");
  SetDescriptorName(iwdescr_ch, "ch");
  SetDescriptorName(iwdescr_htroatom, "htroatom");
  SetDescriptorName(iwdescr_htroaf, "htroaf");
  SetDescriptorName(iwdescr_nrgnhlht, "nrgnhlht");
  SetDescriptorName(iwdescr_ohsh, "ohsh");
  SetDescriptorName(iwdescr_co2h, "co2h");

  if (descriptors_to_compute.specific_groups) {
    SetDescriptorName(iwdescr_amine, "amine");
    SetDescriptorName(iwdescr_pyridine, "pyridine");
    SetDescriptorName(iwdescr_pyrrole, "pyrrole");
  }

  if (descriptors_to_compute.simple_hbond_descriptors) {
    SetDescriptorName(iwdescr_hacts, "hacts");
    SetDescriptorName(iwdescr_hdons, "hdons");
    SetDescriptorName(iwdescr_hduals, "hduals");
  }
  SetDescriptorName(iwdescr_mhr, "mhr");
  SetDescriptorName(iwdescr_mxhrf, "mxhrf");
  SetDescriptorName(iwdescr_mnhrf, "mnhrf");
  SetDescriptorName(iwdescr_lrsysz, "lrsysz");
  SetDescriptorName(iwdescr_srsz, "srsz");
  SetDescriptorName(iwdescr_lrsz, "lrsz");
  SetDescriptorName(iwdescr_rng7atoms, "rng7atoms");
  SetDescriptorName(iwdescr_nrsyscmr, "nrsyscmr");
  SetDescriptorName(iwdescr_mars, "mars");
  if (descriptors_to_compute.spinach_descriptors) {
    SetDescriptorName(iwdescr_frspch, "frspch");
    SetDescriptorName(iwdescr_spchtro, "spchtro");
    SetDescriptorName(iwdescr_rbfrspch, "rbfrspch");
    SetDescriptorName(iwdescr_scafrotb, "scafrotb");
    SetDescriptorName(iwdescr_spchrotb, "spchrotb");
    SetDescriptorName(iwdescr_maxrotbgrp, "maxrotbgrp");
    SetDescriptorName(iwdescr_satspcha, "satspcha");
    SetDescriptorName(iwdescr_unsatspcha, "unsatspcha");
    SetDescriptorName(iwdescr_fsatspcha, "fsatspcha");
    SetDescriptorName(iwdescr_scaffoldbranches, "scaffoldbranches");
    SetDescriptorName(iwdescr_nrnspch, "nrnspch");
    SetDescriptorName(iwdescr_fnrnspc, "fnrnspc");

    SetDescriptorName(iwdescr_trmnlrng, "trmnlrng");
    SetDescriptorName(iwdescr_intrnlrng, "intrnlrng");
    SetDescriptorName(iwdescr_rng2spch, "rng2spch");
    SetDescriptorName(iwdescr_rng2bridge, "rng2bridge");
  }

  if (descriptors_to_compute.ring_chain_descriptors) {
    SetDescriptorName(iwdescr_rcj, "rcj");
    SetDescriptorName(iwdescr_rchj, "rchj");
    SetDescriptorName(iwdescr_amrcj, "amrcj");
    SetDescriptorName(iwdescr_alrcj, "alrcj");
    SetDescriptorName(iwdescr_biphenyl, "biphenyl");
  }

  if (descriptors_to_compute.polar_bond_descriptors) {
    SetDescriptorName(iwdescr_pbcount, "pbcount");
    SetDescriptorName(iwdescr_frpbond, "frpbond");
    SetDescriptorName(iwdescr_nonpbond, "nonpbond");
    SetDescriptorName(iwdescr_pbarom, "pbarom");
    SetDescriptorName(iwdescr_npbarom, "npbarom");
    SetDescriptorName(iwdescr_pbunset, "pbunset");
    SetDescriptorName(iwdescr_dvinylb, "dvinylb");
  }
  SetDescriptorName(iwdescr_ringsys, "ringsys");
  SetDescriptorName(iwdescr_arring, "arring");
  SetDescriptorName(iwdescr_alring, "alring");
  SetDescriptorName(iwdescr_excybond, "excybond");
  SetDescriptorName(iwdescr_excydbond, "excydbond");
  SetDescriptorName(iwdescr_excydscon, "excydscon");
  SetDescriptorName(iwdescr_excydsconh, "excydsconh");
  SetDescriptorName(iwdescr_excydscondon, "excydscondon");
  // SetDescriptorName(iwdescr_scra, "scra");
  // SetDescriptorName(iwdescr_scrha, "scrha");
  // SetDescriptorName(iwdescr_scrd, "scrd");
  SetDescriptorName(iwdescr_atmpiele, "atmpiele");
  SetDescriptorName(iwdescr_fratmpie, "fratmpie");
  SetDescriptorName(iwdescr_unsatura, "unsatura");
  SetDescriptorName(iwdescr_funsatura, "funsatura");
  SetDescriptorName(iwdescr_ringisol, "ringisol");
  SetDescriptorName(iwdescr_isolrc, "isolrc");
  SetDescriptorName(iwdescr_isolhtrc, "isolhtrc");
  SetDescriptorName(iwdescr_erichsct, "erichsct");
  SetDescriptorName(iwdescr_aiercsct, "aiercsct");
  SetDescriptorName(iwdescr_lercsct, "lercsct");
  SetDescriptorName(iwdescr_faiercst, "faiercst");
  SetDescriptorName(iwdescr_avcon, "avcon");
  SetDescriptorName(iwdescr_avchcon, "avchcon");
  SetDescriptorName(iwdescr_avalcon, "avalcon");
  SetDescriptorName(iwdescr_platt, "platt");
  SetDescriptorName(iwdescr_weiner, "weiner");
  SetDescriptorName(iwdescr_internalhbd, "internalhbd");
  if (descriptors_to_compute.crowding_descriptors) {
    SetDescriptorName(iwdescr_crowding, "crowding");
    SetDescriptorName(iwdescr_fcrowdng, "fcrowdng");
  }
  SetDescriptorName(iwdescr_halogen, "halogen");
  SetDescriptorName(iwdescr_halogena, "halogena");
  SetDescriptorName(iwdescr_bigatom, "bigatom");
  SetDescriptorName(iwdescr_fbigatom, "fbigatom");
  SetDescriptorName(iwdescr_csp3, "csp3");
  SetDescriptorName(iwdescr_fcsp3, "fcsp3");
  SetDescriptorName(iwdescr_fccsp3, "fccsp3");
  SetDescriptorName(iwdescr_csp3_chain, "csp3_chain");
  SetDescriptorName(iwdescr_aromc, "aromc");
  SetDescriptorName(iwdescr_aliphc, "aliphc");
  SetDescriptorName(iwdescr_numcdb, "numcdb");
  SetDescriptorName(iwdescr_totdbsub, "totdbsub");
  SetDescriptorName(iwdescr_avcdbsub, "avcdbsub");
  SetDescriptorName(iwdescr_nflxchn, "nflxchn");
  SetDescriptorName(iwdescr_atflxchn, "atflxchn");
  SetDescriptorName(iwdescr_faflxchn, "faflxchn");
  SetDescriptorName(iwdescr_fnflxchn, "fnflxchn");
  SetDescriptorName(iwdescr_lflxchn, "lflxchn");
  SetDescriptorName(iwdescr_avflxchn, "avflxchn");
  SetDescriptorName(iwdescr_rkentrpy, "rkentrpy");
  SetDescriptorName(iwdescr_nconjgsc, "nconjgsc");
  SetDescriptorName(iwdescr_atincnjs, "atincnjs");
  SetDescriptorName(iwdescr_mxcnjscz, "mxcnjscz");
  SetDescriptorName(iwdescr_cinconjs, "cinconjs");
  if (descriptors_to_compute.charge_descriptors) {
    SetDescriptorName(iwdescr_brunsneg, "brunsneg");
    SetDescriptorName(iwdescr_brunspos, "brunspos");
    SetDescriptorName(iwdescr_formal_charge, "formal_charge");
  }
  if (descriptors_to_compute.donor_acceptor) {
    SetDescriptorName(iwdescr_brunsacc, "brunsacc");
    SetDescriptorName(iwdescr_brnsdual, "brnsdual");
    SetDescriptorName(iwdescr_brunsdon, "brunsdon");
    SetDescriptorName(iwdescr_brunshbdsum, "brunshbdsum");
    SetDescriptorName(iwdescr_nplus, "nplus");
    SetDescriptorName(iwdescr_nminus, "nminus");
  }

  if (descriptors_to_compute.distance_matrix_descriptors) {
    SetDescriptorName(iwdescr_muldiam, "muldiam");
    SetDescriptorName(iwdescr_rad, "rad");
    SetDescriptorName(iwdescr_mulrad, "mulrad");
    SetDescriptorName(iwdescr_tm, "tm");
    SetDescriptorName(iwdescr_tg3, "tg3");
    SetDescriptorName(iwdescr_ishape, "ishape");
    SetDescriptorName(iwdescr_maxdrng, "maxdrng");
    SetDescriptorName(iwdescr_maxdarom, "maxdarom");
    SetDescriptorName(iwdescr_maxdhtro, "maxdhtro");
    SetDescriptorName(iwdescr_maxdons, "maxdons");
    SetDescriptorName(iwdescr_avebbtwn, "avebbtwn");
    SetDescriptorName(iwdescr_normbbtwn, "normbbtwn");
    SetDescriptorName(iwdescr_compact, "compact");
    SetDescriptorName(iwdescr_nolp, "nolp");
    SetDescriptorName(iwdescr_avdcentre, "avdcentre");
    SetDescriptorName(iwdescr_stddcentre, "stddcentre");
    SetDescriptorName(iwdescr_centre3, "centre3");
    SetDescriptorName(iwdescr_centre3h, "centre3h");
    SetDescriptorName(iwdescr_mh3b, "mh3b");
    SetDescriptorName(iwdescr_cntrdgncy, "cntrdgncy");
    SetDescriptorName(iwdescr_cntrdshell1, "cntrdshell1");
    SetDescriptorName(iwdescr_cntrdshell2, "cntrdshell2");
    SetDescriptorName(iwdescr_cntrdshell3, "cntrdshell3");

    SetDescriptorName(iwdescr_aveshell1, "aveshell1");
    SetDescriptorName(iwdescr_aveshell2, "aveshell2");
    SetDescriptorName(iwdescr_aveshell3, "aveshell3");
    SetDescriptorName(iwdescr_maxshell3, "maxshell3");
  }

  SetDescriptorName(iwdescr_nnsssrng, "nnsssrng");
  SetDescriptorName(iwdescr_nrings3, "nrings3");
  SetDescriptorName(iwdescr_nrings4, "nrings4");
  SetDescriptorName(iwdescr_nrings5, "nrings5");
  SetDescriptorName(iwdescr_nrings6, "nrings6");
  SetDescriptorName(iwdescr_nrings7, "nrings7");
  SetDescriptorName(iwdescr_nrings8, "nrings8");

  if (descriptors_to_compute.ring_substitution_descriptors) {
    SetDescriptorName(iwdescr_rsarom1, "rsarom1");
    SetDescriptorName(iwdescr_rsarom2, "rsarom2");
    SetDescriptorName(iwdescr_rsarom3, "rsarom3");

    SetDescriptorName(iwdescr_rsaliph1, "rsaliph1");
    SetDescriptorName(iwdescr_rsaliph2, "rsaliph2");
    SetDescriptorName(iwdescr_rsaliph3, "rsaliph3");
    SetDescriptorName(iwdescr_rsaliph4, "rsaliph4");

    SetDescriptorName(iwdescr_rssys1, "rssys1");
    SetDescriptorName(iwdescr_rssys2, "rssys2");
    SetDescriptorName(iwdescr_rssys3, "rssys3");
    SetDescriptorName(iwdescr_rssys4, "rssys4");
    SetDescriptorName(iwdescr_rssys5, "rssys5");
    SetDescriptorName(iwdescr_rssys6, "rssys6");
    SetDescriptorName(iwdescr_rssys7, "rssys7");
    SetDescriptorName(iwdescr_rssys8, "rssys8");
    SetDescriptorName(iwdescr_rssys9, "rssys9");
  }

  if (descriptors_to_compute.ring_fusion_descriptors) {
    SetDescriptorName(iwdescr_ar5, "ar5");
    SetDescriptorName(iwdescr_ar6, "ar6");
    SetDescriptorName(iwdescr_al5, "al5");
    SetDescriptorName(iwdescr_al6, "al6");

    SetDescriptorName(iwdescr_fsdrng5l5l, "fsdrng5l5l");
    SetDescriptorName(iwdescr_fsdrng5l5r, "fsdrng5l5r");
    SetDescriptorName(iwdescr_fsdrng5r5r, "fsdrng5r5r");
    SetDescriptorName(iwdescr_fsdrng5l6l, "fsdrng5l6l");
    SetDescriptorName(iwdescr_fsdrng5l6r, "fsdrng5l6r");
    SetDescriptorName(iwdescr_fsdrng5r6l, "fsdrng5r6l");
    SetDescriptorName(iwdescr_fsdrng5r6r, "fsdrng5r6r");
    SetDescriptorName(iwdescr_fsdrng6r6r, "fsdrng6r6r");
    SetDescriptorName(iwdescr_fsdrng6l6r, "fsdrng6l6r");
    SetDescriptorName(iwdescr_fsdrng6l6l, "fsdrng6l6l");

    SetDescriptorName(iwdescr_fsdrngarar, "fsdrngarar");
    SetDescriptorName(iwdescr_fsdrngalar, "fsdrngalar");
    SetDescriptorName(iwdescr_fsdrngalal, "fsdrngalal");
  }

  SetDescriptorName(iwdescr_nchiral, "nchiral");
  if (descriptors_to_compute.psa) {
    SetDescriptorName(iwdescr_nvrtspsa, "nvrtspsa");
  }
  if (descriptors_to_compute.long_carbon_chains) {
    SetDescriptorName(iwdescr_mxlencchain2, "mxlencchain2");
    SetDescriptorName(iwdescr_mxlencchain3, "mxlencchain3");
  }
  if (descriptors_to_compute.charge_descriptors) {
    SetDescriptorName(iwdescr_acmbe, "acmbe");
  }
  SetDescriptorName(iwdescr_cmr, "cmr");
  SetDescriptorName(iwdescr_cd4ring, "cd4ring");
  SetDescriptorName(iwdescr_cd4chain, "cd4chain");
  if (descriptors_to_compute.ring_substitution_ratio_descriptors) {
    SetDescriptorName(iwdescr_frsub, "frsub");
    SetDescriptorName(iwdescr_frssub, "frssub");
    SetDescriptorName(iwdescr_arorthoring, "arorthoring");
    SetDescriptorName(iwdescr_alorthoring, "alorthoring");
  }
  if (descriptors_to_compute.bonds_between_rings) {
    SetDescriptorName(iwdescr_bbr1, "bbr1");
    SetDescriptorName(iwdescr_bbr2, "bbr2");
    SetDescriptorName(iwdescr_bbr3, "bbr3");
    SetDescriptorName(iwdescr_bbr4, "bbr4");
    SetDescriptorName(iwdescr_bbr5, "bbr5");
    SetDescriptorName(iwdescr_bbr6, "bbr6");
  }
  if (descriptors_to_compute.adjacent_ring_fusion_descriptors) {
    SetDescriptorName(iwdescr_sboradjf, "sboradjf");
    SetDescriptorName(iwdescr_dboradjf, "dboradjf");
  }
  SetDescriptorName(iwdescr_hcount, "hcount");
  SetDescriptorName(iwdescr_hperatom, "hperatom");
  SetDescriptorName(iwdescr_ro5_ohnh, "ro5_ohnh");
  SetDescriptorName(iwdescr_ro5_on, 
      "ro5_on");  // the last one which will always be computed

  if (descriptors_to_compute.donor_acceptor && min_hbond_feature_separation > 0) {
    SetDescriptorName(iwdescr_aamind, "aamind");
    SetDescriptorName(iwdescr_aa2mind, "aa2mind");
    SetDescriptorName(iwdescr_aaave, "aaave");
    SetDescriptorName(iwdescr_admind, "admind");
    SetDescriptorName(iwdescr_ad2mind, "ad2mind");
    SetDescriptorName(iwdescr_adave, "adave");
    SetDescriptorName(iwdescr_ddmind, "ddmind");
    SetDescriptorName(iwdescr_dd2mind, "dd2mind");
    SetDescriptorName(iwdescr_ddave, "ddave");
  }

  if (descriptors_to_compute.complexity_descriptors) {
    SetDescriptorName(iwdescr_nspiro, "nspiro");
    SetDescriptorName(iwdescr_nsfsdsys, "nsfsdsys");
    SetDescriptorName(iwdescr_rnginsfs, "rnginsfs");
    SetDescriptorName(iwdescr_lgstrfsy, "lgstrfsy");
    SetDescriptorName(iwdescr_htrcsfsy, "htrcsfsy");
    SetDescriptorName(iwdescr_mxhtsfsy, "mxhtsfsy");
    SetDescriptorName(iwdescr_npfsdsys, "npfsdsys");
    SetDescriptorName(iwdescr_rnginpfs, "rnginpfs");
    SetDescriptorName(iwdescr_lgplnfsy, "lgplnfsy");
    SetDescriptorName(iwdescr_htrcpfsy, "htrcpfsy");
    SetDescriptorName(iwdescr_mxhtpfsy, "mxhtpfsy");
  }

  if (descriptors_to_compute.symmetry_descriptors) {
    SetDescriptorName(iwdescr_symmatom, "symmatom");
    SetDescriptorName(iwdescr_fsymmatom, "fsymmatom");
    SetDescriptorName(iwdescr_lsepsymatom, "lsepsymatom");
    SetDescriptorName(iwdescr_flsepsymatom, "flsepsymatom");
    SetDescriptorName(iwdescr_maxsymmclass, "maxsymmclass");
  }

  if (descriptors_to_compute.partial_symmetry_descriptors) {
    SetDescriptorName(iwdescr_maxpsymd, "maxpsymd");
    SetDescriptorName(iwdescr_fmaxpsymd, "fmaxpsymd");
    SetDescriptorName(iwdescr_maxpsymdmean, "maxpsymdmean");
    SetDescriptorName(iwdescr_psymdnumzero, "psymdnzero");
  }

  if (descriptors_to_compute.ramey_descriptors) {
    SetDescriptorName(iwdescr_obalance, "obalance");
    SetDescriptorName(iwdescr_rmync, "rmync");
    SetDescriptorName(iwdescr_rmynn, "rmynn");
    SetDescriptorName(iwdescr_rmyno, "rmyno");
    SetDescriptorName(iwdescr_rmynf, "rmynf");
    SetDescriptorName(iwdescr_rmyns, "rmyns");
    SetDescriptorName(iwdescr_rmyncl, "rmyncl");
    SetDescriptorName(iwdescr_rmynbr, "rmynbr");
    SetDescriptorName(iwdescr_rmyni, "rmyni");
    SetDescriptorName(iwdescr_rmy_heavy_halogen, "heavy_halogen");
  }
  if (descriptors_to_compute.compute_alogp) {
    SetDescriptorName(iwdescr_alogp, "alogp");
  }
  if (descriptors_to_compute.compute_xlogp) {
    SetDescriptorName(iwdescr_xlogp, "xlogp");
  }

  if (descriptors_to_compute.distance_matrix_descriptors) {
    MarkBestFingerprint(iwdescr_maxdarom);
  }
  if (descriptors_to_compute.psa) {
    MarkBestFingerprint(iwdescr_nvrtspsa);
  }
  if (descriptors_to_compute.long_carbon_chains) {
    MarkBestFingerprint(iwdescr_mxlencchain2);
    MarkBestFingerprint(iwdescr_mxlencchain3);
  }
  if (descriptors_to_compute.saturated_chains) {
    SetDescriptorName(iwdescr_nsatchain, "nsatchain");
    SetDescriptorName(iwdescr_mxsatchain, "mxsatchain");
    SetDescriptorName(iwdescr_fsatchain, "fsatchain");
  }
  MarkBestFingerprint(iwdescr_natoms);
  MarkBestFingerprint(iwdescr_frafus);
  if (descriptors_to_compute.distance_matrix_descriptors) {
    MarkBestFingerprint(iwdescr_maxdrng);
  }
  MarkBestFingerprint(iwdescr_aromc);
  MarkBestFingerprint(iwdescr_ro5_ohnh);
  MarkBestFingerprint(iwdescr_rmync);
  MarkBestFingerprint(iwdescr_fraromha);
  MarkBestFingerprint(iwdescr_ringatom);
  MarkBestFingerprint(iwdescr_rhaf);
  MarkBestFingerprint(iwdescr_rgmltbd);
  MarkBestFingerprint(iwdescr_ncon3);
  MarkBestFingerprint(iwdescr_brunsdon);
  MarkBestFingerprint(iwdescr_brunshbdsum);
  MarkBestFingerprint(iwdescr_avebbtwn);
  MarkBestFingerprint(iwdescr_amine);
  MarkBestFingerprint(iwdescr_npbarom);
  MarkBestFingerprint(iwdescr_hacts);
  MarkBestFingerprint(iwdescr_fchmltbd);
  MarkBestFingerprint(iwdescr_fnrnspc);
  MarkBestFingerprint(iwdescr_atmpiele);
  MarkBestFingerprint(iwdescr_ncon2);
  MarkBestFingerprint(iwdescr_funsatura);
  MarkBestFingerprint(iwdescr_csp3);
  MarkBestFingerprint(iwdescr_ch);
  MarkBestFingerprint(iwdescr_avsdlp);


#ifdef TEST_FOR_NAMES
  return VerifyDescriptorNames();
#else
  return 1;
#endif
}

void
IWDescr::IWDescrImpl::SetDescriptorName(int descriptor_number, const char* name) {
  if (descriptor == nullptr || descriptor_number < 0 ||
      descriptor_number >= number_descriptors) {
    return;
  }

  descriptor[descriptor_number].set_name(name);
}

void
IWDescr::IWDescrImpl::MarkBestFingerprint(int descriptor_number) {
  if (descriptor == nullptr || descriptor_number < 0 ||
      descriptor_number >= number_descriptors) {
    return;
  }

  descriptor[descriptor_number].set_best_fingerprint(1);
}

int
IWDescr::IWDescrImpl::VerifyDescriptorNames() const {
  if (descriptor == nullptr) {
    return 0;
  }

  std::unordered_set<std::string> names;
  names.reserve(number_descriptors);
  for (int i = 0; i < number_descriptors; ++i) {
    if (descriptor[i].descriptor_name().empty()) {
      cerr << "IWDescr::VerifyDescriptorNames:descriptor " << i << " has no name\n";
      return 0;
    }
    const std::string name = descriptor[i].descriptor_name().AsString();
    if (! names.insert(name).second) {
      cerr << "IWDescr::VerifyDescriptorNames:duplicate descriptor name '" << name << "'\n";
      return 0;
    }
  }

  return 1;
}

// Currently turned off. Descriptors have initialised with Nan values
// elsewhere.
void
IWDescr::IWDescrImpl::InitialiseDescriptorDefaults() {
  if (descriptor == nullptr) {
    return;
  }

//static constexpr float kZero = 0.0f;
//for (int i = 0; i < number_descriptors; ++i) {
//  descriptor[i].set_default_value(kZero);
//}
}

int
IWDescr::IWDescrImpl::Process(Molecule& m, float* results) {
  if (results == nullptr) {
    return 0;
  }

  ResetDescriptors();

  // All descriptor computation happens here. This is deliberately separated
  // from CopyDescriptorsToResults so the legacy descriptor functions can be
  // migrated one family at a time into IWDescrImpl methods.
  if (! ComputeDescriptorValues(m)) {
    return 0;
  }

  ++molecules_processed;
  CopyDescriptorsToResults(results);

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeDescriptorValues(Molecule& m) {
  // This mirrors the old iwdescriptors(Molecule&, output_separator, output)
  // wrapper, except it builds calculation scratch data and then delegates to a
  // computation-only method. There is no output object and no output separator
  // in the IWDescr path.
  PerMoleculeData data(m);
  return ComputeDescriptorValues(m, data);
}

int
IWDescr::IWDescrImpl::ComputeDescriptorValues(Molecule& m, PerMoleculeData& data) {
  // This is the migrated shape of the legacy
  // iwdescriptors(Molecule&, output, z, ncon, ring_membership, conn, atom,
  // is_aromatic, output_separator) driver. It deliberately preserves the old
  // descriptor-family ordering, but removes all output concerns from the
  // computation path.

  ComputeChiralityDescriptors(m, data);

  // Historical behaviour: once chirality descriptors are computed, chiral
  // centres are removed. This means Process() is not currently non-mutating.
  m.remove_all_chiral_centres();

  // cerr << "Going to ComputeTopologicalDescriptors\n";

  ComputeTopologicalDescriptors(m, data);

  if (descriptors_to_compute.crowding_descriptors) {
      ComputeCrowdingDescriptors(m, data);
  }

  if (descriptors_to_compute.spinach_descriptors) {
      ComputeSpinachDescriptors(m, data);
  }

//if (descriptors_to_compute.ring_chain_descriptors) {
//    ComputeRingChainDescriptors(m, data);
//}

  if (descriptors_to_compute.complexity_descriptors) {
      ComputeComplexityDescriptors(m, data);
  }

  // The Novartis PSA descriptor historically precedes charge assignment.
  ComputeNovartisPsaDescriptor(m, data);

  ComputeMolarRefractivityDescriptors(m, data);

  ComputeRuleOfFiveDescriptors(m, data);

  if (descriptors_to_compute.distance_matrix_descriptors) {
      ComputeDistanceMatrixDescriptors(m, data);
  }

  if (descriptors_to_compute.symmetry_descriptors) {
      ComputeSymmetryDescriptors(m, data);
  }

  if (descriptors_to_compute.partial_symmetry_descriptors) {
    NearSymmetricDescriptors(m);
  }

  if (descriptors_to_compute.charge_descriptors && charge_assigner.active()) {
    ComputeChargeDescriptors(m, data);

    // Preserve the legacy refresh after charge assignment.
    data.RefreshConnections(m);
  }

  if (descriptors_to_compute.donor_acceptor && donor_acceptor_assigner.active()) {
      ComputeDonorAcceptorDescriptors(m, data);
  }

  // Partial-charge descriptor families remain intentionally omitted, matching
  // the legacy `if (0 && m.has_charges())` block. If that is resurrected later,
  // it should become a private IWDescrImpl computation hook, not main-side I/O.
  return 1;
}

/*
  Counts the number of 4 connected carbons that are in chains and rings
  asymc is the number of 4 connected atoms with 4 different atom types connected.
*/

int
IWDescr::IWDescrImpl::compute_4_connected_carbon_stuff(const Molecule& m, const atomic_number_t* z,
                                 const int* ncon, const int* ring_membership) {
  int cd4_ring = 0;
  int cd4_chain = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    if (4 != ncon[i]) {
      continue;
    }

    if (6 != z[i]) {
      continue;
    }

    if (ring_membership[i]) {
      cd4_ring++;
    } else {
      cd4_chain++;
    }
  }

  descriptor[iwdescr_cd4ring].set(static_cast<float>(cd4_ring));
  descriptor[iwdescr_cd4chain].set(static_cast<float>(cd4_chain));

  return 1;
}


bool
AnyAtomsMoreThanTwoRings(Molecule& m, const Ring& r) {
  for (atom_number_t a : r) {
    // cerr << " atom " << a << " in " << m.nrings(a) << '\n';
    if (m.nrings(a) > 2) {
      return true;
    }
  }

  return false;
}

// Turn off any atom that is in a strongly fused ring.
// This is unstable wrt SSSR ring determination.
// Perhaps not as precise as I would like, but it is not dependent on 
// SSSR ring perception.
void
TurnOffStronglyFused(Molecule& m, const atomic_number_t* z,
                        const int* ncon, int* ok_to_check) {
  if (m.nrings() < 2) {
    return;
  }

  m.ring_membership();

  // First turn off all ring atoms.
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) > 1) {
      ok_to_check[i] = 0;
    }
  }

  // Identify ring atoms at which we can do expensive chirality computations.
  auto possible_chiral = [&](atom_number_t zatom)->bool {
    if (z[zatom] != 6) {
      return false;
    }
    if (ncon[zatom] != 3 && ncon[zatom] != 4) {
      return false;
    }
    if (m.is_aromatic(zatom)) {
      return false;
    }

    const int rbc = m.ring_bond_count(zatom);

    if (rbc == 0) {
      return false;
    }

    if (rbc > 2) {
      return false;
    }

    return true;
  };

  // All ring atoms have been turned off. Turn back on those that
  // are eligible - two [D2] connections.
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    if (! possible_chiral(i)) {
      continue;
    }

    int adjacent_D2 = 0;  // within the ring
    for (const Bond* b : m[i]) {
      if (b->nrings() == 0) {
        continue;
      }
      if (! b->is_single_bond()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (ncon[o] == 2) {
        ++adjacent_D2;
      }
    }

    if (adjacent_D2 == 2) {
      ok_to_check[i] = 1;
    }
  }

#ifdef UNSTABLE_WRT_SSR
  for (int i = 0; i < nrings; ++i) {
    const Ring* ri = m.ringi(i);
    //  cerr << "What abt " << *ri << '\n';
    //  cerr << "Bonds shared " << ri->largest_number_of_bonds_shared_with_another_ring()
    //  << '\n';
    if (!ri->is_fused()) {
      continue;
    }

    if (!AnyAtomsMoreThanTwoRings(m, *ri)) {
      continue;
    }

    ri->set_vector(ok_to_check, 0);

    const int rsize = ri->number_elements();
    atom_number_t previous_atom = ri->back();

    for (int i = 0; i < rsize; ++i) {
      atom_number_t a = (*ri)[i];
      if (ncon[a] == 2 || m.nrings_including_non_sssr_rings(a) > 1) {
        previous_atom = a;
        continue;
      }

      // Allow chiral centres on "uncluttered" atoms.
      int j = i;
      atom_number_t next_atom = ri->next_after_wrap(j, 1);
      if (ncon[previous_atom] == 2 && ncon[next_atom] == 2) {
        ok_to_check[a] = 1;
      }

      previous_atom = a;
    }
  }
#endif
}


int
IWDescr::IWDescrImpl::ComputeChiralityDescriptors(Molecule& m, PerMoleculeData& data) {
  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();

  std::unique_ptr<int[]> ok_to_check;
  if (descriptors_to_compute.perform_expensive_chirality_perception) {
    ok_to_check.reset(new_int(matoms, 1));
    TurnOffStronglyFused(m, z, ncon, ok_to_check.get());
  }

  int chiral_centres = 0;
#ifdef DEBUG_CHIRALITY
  Set_of_Atoms cc;
#endif

  for (int i = 0; i < matoms; ++i) {
    if (ncon[i] < 3) {
      continue;
    }

    if (6 == z[i]) {
      if (4 == ncon[i]) {
      } else if (3 == ncon[i] && 1 == m.implicit_hydrogens(i)) {
      } else {
        continue;
      }
    } else if (16 == z[i]) {
      if (3 == ncon[i] && 4 == m.nbonds(i)) {
      } else {
        continue;
      }
    } else {
      continue;
    }

    if (nullptr != m.chiral_centre_at_atom(i)) {
      ++chiral_centres;
      continue;
    }

    if (!descriptors_to_compute.perform_expensive_chirality_perception) {
      continue;
    }

    if (!ok_to_check[i]) {
      continue;
    }

    // is_actually_chiral() is expensive, so try to avoid if possible.
    // Careful what we check: C1C2C3C1N1C2C31 has 5 chiral centers
    // otherwise. By using non-SSSR rings, all atoms in cubane will not be
    // checked.
    if (is_actually_chiral(m, i)) {
      ++chiral_centres;
#ifdef DEBUG_CHIRALITY
      cc << i;
#endif
    }
  }

#ifdef DEBUG_CHIRALITY
  if (cc.size() > 0) {
    cerr << "Unlabelled chiral centres " << cc << '\n';
    write_isotopically_labelled_smiles(m, true, cerr);
  }
#endif

  descriptor[iwdescr_nchiral].set(static_cast<float>(chiral_centres));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeCoreTopologicalDescriptors(Molecule& m, PerMoleculeData& data,
                                                       int* already_done) {
  TopologicalDescriptorCounts counts;

  // cerr << "Goingto PopulateTopologicalCounts\n";
  if (! PopulateTopologicalCounts(m, data, already_done, counts)) {
    return 0;
  }

  //  cerr << "Going to StoreTopologicalCounts\n";
  return StoreTopologicalCounts(m, data, counts);
}

int
IWDescr::IWDescrImpl::PopulateTopologicalCounts(Molecule& m, PerMoleculeData& data,
                                                int* already_done,
                                                TopologicalDescriptorCounts& counts) {
  // Migrated outline of the large atom/bond loop from legacy
  // compute_topological_descriptors().  This method owns accumulation only;
  // StoreTopologicalCounts() remains responsible for writing Descriptor values.
  if (already_done == nullptr || data.matoms == 0) {
    return 0;
  }

  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const Atom** atom = data.atoms();
  const int* is_aromatic = data.aromaticity();

  // Preserve the legacy force of SSSR perception before the atom loop.
  m.ring_membership();

  for (int i = 0; i < matoms; ++i) {
    counts.total_connectivity += ncon[i];

    if (ncon[i] <= 4) {
      ++counts.connected[ncon[i]];
    } else {
      ++counts.connected[5];
    }

    if (ncon[i] == 2 && m.saturated(i)) {
      ++counts.d2sp3;
    }

    const atomic_number_t zi = z[i];
    if (zi == 1) {
      ++counts.hcount;
    } else {
      ++counts.heavy_atom_count;
      counts.hcount += m.hcount(i);
    }

    if (zi != 6) {
      ++counts.heteroatom_count;
    }

    if (zi > 10) {
      ++counts.bigatom_count;
    }

    if ((zi == 8 || zi == 16) && ncon[i] == 1) {
      ++counts.singly_connected_oxygen_or_sulphur_count;
    }

    // Loss of const OK.
    Atom* ai = const_cast<Atom*>(atom[i]);

    if (zi == 8 && ncon[i] == 1 && ai->nbonds() == 2) {
      counts.carboxylic_acid_count += IsCarboxylicAcid(atom, i, z, ncon);
    }

    if (zi == 6 && ncon[i] == ai->nbonds()) {
      ++counts.csp3;
      if (ring_membership[i] == 0) {
        ++counts.csp3_chain;
      }
    }

    const int is_halogen = (zi == 9 || zi == 17 || zi == 35 || zi == 53);
    if (is_halogen) {
      ++counts.halogen_count;
    }

    if (zi != 6 && ring_membership[i] == 0 && ! is_halogen) {
      ++counts.non_ring_non_halogen_heteroatoms;
    }

    bool aromatic = false;
    if (ring_membership[i]) {
      ++counts.ring_atom_count;

      if (is_aromatic[i]) {
        ++counts.aromatic_atom_count;
        aromatic = true;
      }

      if (zi != 6) {
        ++counts.ring_heteroatom_count;
        if (aromatic) {
          ++counts.aromatic_heteroatom_count;
        }
      }

      if (ring_membership[i] > 1) {
        ++counts.fused_ring_atom_count;
      }
    } else {
      counts.total_chain_connectivity += ncon[i];

      if (ncon[i] == 2) {
        ++counts.two_connected_chain_atom;
      }
    }

    if (zi == 6 && ! aromatic) {
      ++counts.aliphatic_carbon_count;
    }

    if (zi == 6 && ai->implicit_hydrogens()) {
      ++counts.carbon_hydrogen_count;
      if (ai->implicit_hydrogens() == 2) {
        ++counts.ch2;
      }
    }

    int lone_pairs = 0;
    bool electron_rich = false;

    if (aromatic) {
      ++counts.atoms_with_pi_electrons;
      electron_rich = true;
    } else {
      counts.total_aliphatic_connectivity += ncon[i];

      if (ncon[i] < ai->nbonds()) {
        ++counts.atoms_with_pi_electrons;
        ++counts.unsaturation;
        electron_rich = true;
      } else if (ai->lone_pair_count(lone_pairs) && lone_pairs) {
        ++counts.atoms_with_pi_electrons;
        electron_rich = true;
      }
    }

    if (electron_rich && already_done[i] == 0) {
      ++counts.electron_rich_sections;
      const int section_size = ComputeElectronRichSection(m, i, atom, already_done);
      if (section_size > counts.largest_electron_rich_section) {
        counts.largest_electron_rich_section = section_size;
      }
      counts.atoms_in_electron_rich_sections += section_size;
    }

    for (int j = 0; j < ncon[i]; ++j) {
      const Bond* b = ai->item(j);
      const atom_number_t k = b->other(i);

      if (k < i) {
        continue;
      }

      if (b->nrings()) {
        if (! b->is_aromatic() && ! b->is_single_bond()) {
          ++counts.ring_multiple_bonds;
        }
        continue;
      }

      if (b->is_single_bond() && ncon[i] > 1 && ncon[k] > 1 &&
          ! TripleBondAtEitherEnd(m, b) &&
          ! part_of_otherwise_non_rotabable_entity(m, i, k)) {
        ++counts.rotatable_bonds;
      } else if (! b->is_single_bond()) {
        ++counts.chain_multiple_bonds;
      }
    }
  }

  return 1;
}

int
IWDescr::IWDescrImpl::IsCarboxylicAcid(const Atom** atom, atom_number_t oxygen,
                                        const atomic_number_t* z, const int* ncon) const {
  // Migrated from legacy carboxylic_acid(atom, oxygen, z, ncon). The caller
  // has already identified `oxygen` as a singly-connected doubly bonded oxygen.
  const atom_number_t carbon = atom[oxygen]->other(oxygen, 0);

  if (ncon[carbon] != 3 || z[carbon] != 6) {
    return 0;
  }

  const Atom* ac = atom[carbon];

  int found_singly_bonded_oxygen = 0;
  for (int i = 0; i < 3; ++i) {
    const Bond* b = ac->item(i);
    if (! b->is_single_bond()) {
      continue;
    }

    const atom_number_t j = b->other(carbon);
    if (z[j] == 8 && ncon[j] == 1) {
      ++found_singly_bonded_oxygen;
    }
  }

  return found_singly_bonded_oxygen == 1;
}

int
IWDescr::IWDescrImpl::ComputeElectronRichSection(Molecule& m, atom_number_t start,
                                                  const Atom** atom,
                                                  int* already_done) const {
  // Migrated from legacy compute_electron_rich(m, start, atom, already_done).
  // `m` is not used by the legacy traversal, but remains in the signature to
  // keep this helper aligned with the other descriptor-family hooks.
  (void)m;

  already_done[start] = 1;

  int rc = 1;
  const Atom* ai = atom[start];

  for (const Bond* b : *ai) {
    const atom_number_t j = b->other(start);
    if (already_done[j]) {
      continue;
    }

    // lone_pair_count is non const.
    Atom* aj = const_cast<Atom*>(atom[j]);

    bool electron_rich = false;
    int lone_pairs = 0;

    if (aj->ncon() < aj->nbonds()) {
      electron_rich = true;
    } else if (aj->lone_pair_count(lone_pairs) && lone_pairs) {
      electron_rich = true;
    }

    if (electron_rich) {
      rc += ComputeElectronRichSection(m, j, atom, already_done);
    }
  }

  return rc;
}

int
IWDescr::IWDescrImpl::TripleBondAtEitherEnd(Molecule& m, const Bond* b) const {
  // A single bond adjacent to a triple bond is not counted as freely rotatable
  // by the legacy rotatable-bond heuristic. Keep this as an implementation
  // helper because it only affects topological descriptor values.
  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  for (const Bond* b1 : m[a1]) {
    if (b1 == b) {
      continue;
    }
    if (b1->is_triple_bond()) {
      return 1;
    }
  }

  for (const Bond* b2 : m[a2]) {
    if (b2 == b) {
      continue;
    }
    if (b2->is_triple_bond()) {
      return 1;
    }
  }

  return 0;
}

int
IWDescr::IWDescrImpl::PartOfOtherwiseNonRotatableEntity(Molecule& m, atom_number_t a1,
                                                        atom_number_t a2) const {
  // Conservative implementation of the legacy non-rotatable single-bond
  // exclusions. The important case is an amide-like single bond where one
  // endpoint is carbon and the other is N/O/S, and the carbon has a multiple
  // bond to a hetero atom.
  auto is_amide_like_pair = [&m](atom_number_t carbon, atom_number_t hetero) -> int {
    if (m.atomic_number(carbon) != 6) {
      return 0;
    }

    const atomic_number_t zhetero = m.atomic_number(hetero);
    if (zhetero != 7 && zhetero != 8 && zhetero != 16) {
      return 0;
    }

    const Atom& c = m[carbon];
    if (c.ncon() == c.nbonds()) {
      return 0;
    }

    for (const Bond* b : c) {
      const atom_number_t j = b->other(carbon);
      if (j == hetero) {
        continue;
      }
      if (! b->is_single_bond() && m.atomic_number(j) != 6) {
        return 1;
      }
    }

    return 0;
  };

  return is_amide_like_pair(a1, a2) || is_amide_like_pair(a2, a1);
}

int
IWDescr::IWDescrImpl::StoreCoreTopologicalCountDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // This method is now primarily a dispatcher for topological descriptor
  // storage. Keep new descriptor-store blocks in small private Store* helpers.

  // cerr << "StoreCoreTopologicalCountDescriptors\n";
  if (descriptor == nullptr || data.matoms == 0) {
    return 0;
  }

  StoreAtomCompositionDescriptors(m, data, counts);

  StoreBondAndRingDescriptors(m, data, counts);

  StorePiElectronDescriptors(m, data, counts);

  StoreConnectivityDescriptors(m, data, counts);

  ComputeHalogenAttachmentDescriptors(m, data, counts.halogen_count);

  ComputeRadhaEntropyDescriptors(m, data);


  return 1;
}

int
IWDescr::IWDescrImpl::StoreAtomCompositionDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Split from StoreTopologicalCounts(). Keep the top-level atom-composition
  // grouping, but delegate the basic atom/molecular-weight descriptors and the
  // carbon/heteroatom descriptors to smaller private helpers.
  // cerr << "Going to StoreBasicAtomDescriptors\n";
  StoreBasicAtomDescriptors(m, data, counts);

  // cerr << "Goint to StoreCarbonAndHeteroatomDescriptors\n";
  StoreCarbonAndHeteroatomDescriptors(m, data, counts);

  return 1;
}

int
IWDescr::IWDescrImpl::StoreBasicAtomDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Atom count, molecular weight, elemental diversity, hydrogens, and coarse
  // non-organic atom counts.
  const int matoms = data.matoms;

  descriptor[iwdescr_natoms].set(static_cast<float>(matoms));

  // An isotopic label must not change the weight, and must not be read as a
  // statement that the atom has no Hydrogens. Agrees with molecule_filter,
  // fileconv and the python bindings.
  descriptor[iwdescr_amw].set(static_cast<float>(lillymol::MolecularWeightIsotopesAsLabels(m)));

  descriptor[iwdescr_bigatom].set(static_cast<float>(counts.bigatom_count));
  descriptor[iwdescr_fbigatom].set(static_cast<float>(counts.bigatom_count) /
                                   static_cast<float>(matoms));
  descriptor[iwdescr_halogen].set(static_cast<float>(counts.halogen_count));
  descriptor[iwdescr_nelem].set(static_cast<float>(m.number_different_elements()));

  descriptor[iwdescr_hcount].set(static_cast<float>(counts.hcount));
  descriptor[iwdescr_hperatom].set(static_cast<float>(counts.hcount) /
                                   static_cast<float>(matoms));

  return 1;
}

int
IWDescr::IWDescrImpl::StoreCarbonAndHeteroatomDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  StoreCarbonHybridisationDescriptors(m, data, counts);

  StoreHeteroatomDescriptors(m, data, counts);

  return 1;
}

int
IWDescr::IWDescrImpl::StoreCarbonHybridisationDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  (void)m;
  const int matoms = data.matoms;

  descriptor[iwdescr_csp3].set(static_cast<float>(counts.csp3));
  descriptor[iwdescr_fcsp3].set(static_cast<float>(counts.csp3) /
                                static_cast<float>(matoms));

  if (counts.heteroatom_count == matoms) {
    // Preserve the legacy convention for all-heteroatom molecules.
    descriptor[iwdescr_fccsp3].set(0.0f);
  } else {
    descriptor[iwdescr_fccsp3].set(static_cast<float>(counts.csp3) /
                                   static_cast<float>(matoms - counts.heteroatom_count));
  }

  descriptor[iwdescr_csp3_chain].set(static_cast<float>(counts.csp3_chain));

  descriptor[iwdescr_ch2].set(static_cast<float>(counts.ch2));
  descriptor[iwdescr_d2sp3].set(static_cast<float>(counts.d2sp3));
  descriptor[iwdescr_ch].set(static_cast<float>(counts.carbon_hydrogen_count));

  return 1;
}

int
IWDescr::IWDescrImpl::StoreHeteroatomDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  (void)m;
  (void)data;

  descriptor[iwdescr_htroatom].set(static_cast<float>(counts.heteroatom_count));
  if (counts.heavy_atom_count > 0) {
    descriptor[iwdescr_htroaf].set(static_cast<float>(counts.heteroatom_count) /
                                   static_cast<float>(counts.heavy_atom_count));
  } else {
    descriptor[iwdescr_htroaf].set(0.0f);
  }

  descriptor[iwdescr_nrgnhlht].set(
      static_cast<float>(counts.non_ring_non_halogen_heteroatoms));
  descriptor[iwdescr_ohsh].set(
      static_cast<float>(counts.singly_connected_oxygen_or_sulphur_count));
  descriptor[iwdescr_co2h].set(static_cast<float>(counts.carboxylic_acid_count));

  return 1;
}

int
IWDescr::IWDescrImpl::StoreBondAndRingDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Split from StoreTopologicalCounts(). These assignments cover connection
  // counts, multiple bonds, rotatable bonds, and ring/aromatic atom counts.
  // cerr << "Inside StoreBondAndRingDescriptors\n";
  StoreConnectionCountDescriptors(m, data, counts);

  StoreMultipleBondDescriptors(m, data, counts);

  StoreRotatableBondDescriptors(m, data, counts);

  StoreRingAndAromaticDescriptors(m, data, counts);

  return 1;
}

int
IWDescr::IWDescrImpl::StoreConnectionCountDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Connection-count descriptors and highly connected atom fraction.
  const int matoms = data.matoms;

  if (descriptors_to_compute.ncon_descriptors) {
    const int offset = iwdescr_ncon1;
    for (int i = 1; i < 5; ++i) {
      descriptor[offset + 2 * (i - 1)].set(static_cast<float>(counts.connected[i]));
      descriptor[offset + 2 * (i - 1) + 1].set(
          static_cast<float>(counts.connected[i]) / static_cast<float>(matoms));
    }
  }

  int highly_connected = 0;
  for (int i = 2; i < 5; ++i) {
    highly_connected += counts.connected[i];
  }

  descriptor[iwdescr_frhc].set(static_cast<float>(highly_connected) /
                               static_cast<float>(matoms));

  return 1;
}

int
IWDescr::IWDescrImpl::StoreMultipleBondDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Multiple-bond counts and fractions.
  int mbonds = m.nedges();
  if (mbonds == 0) [[unlikely]] {
    // Preserve the legacy convention for single atom molecules.
    mbonds = 1;
  }

  const int multiple_bonds = counts.chain_multiple_bonds + counts.ring_multiple_bonds;
  descriptor[iwdescr_mltbd].set(static_cast<float>(multiple_bonds));
  descriptor[iwdescr_fmltbd].set(static_cast<float>(multiple_bonds) /
                                 static_cast<float>(mbonds));
  descriptor[iwdescr_chmltbd].set(static_cast<float>(counts.chain_multiple_bonds));
  descriptor[iwdescr_fchmltbd].set(static_cast<float>(counts.chain_multiple_bonds) /
                                   static_cast<float>(mbonds));
  descriptor[iwdescr_rgmltbd].set(static_cast<float>(counts.ring_multiple_bonds));
  descriptor[iwdescr_frgmltbd].set(static_cast<float>(counts.ring_multiple_bonds) /
                                   static_cast<float>(mbonds));

  return 1;
}

int
IWDescr::IWDescrImpl::StoreRotatableBondDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Chain two-connected atoms and rotatable-bond counts/fractions.
  const int matoms = data.matoms;

  descriptor[iwdescr_dcca].set(static_cast<float>(counts.two_connected_chain_atom));
  descriptor[iwdescr_fdcca].set(static_cast<float>(counts.two_connected_chain_atom) /
                                static_cast<float>(matoms));

  descriptor[iwdescr_rotbond].set(static_cast<float>(counts.rotatable_bonds));
  if (m.nedges() == 0) {
    descriptor[iwdescr_frotbond].set(0.0f);
  } else {
    descriptor[iwdescr_frotbond].set(
        iwmisc::Fraction<float>(counts.rotatable_bonds, m.nedges()));
  }

  return 1;
}

int
IWDescr::IWDescrImpl::StoreRingAndAromaticDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Ring/fused-ring and aromatic descriptors are separate legacy store blocks
  // but are kept behind this top-level grouping to preserve the current
  // StoreBondAndRingDescriptors() call structure.
  StoreRingDescriptors(m, data, counts);

  StoreAromaticDescriptors(m, data, counts);

  return 1;
}

int
IWDescr::IWDescrImpl::StoreRingDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Ring and fused-ring atom descriptors.
  const int matoms = data.matoms;

  descriptor[iwdescr_ringatom].set(static_cast<float>(counts.ring_atom_count));
  if (counts.ring_atom_count) {
    descriptor[iwdescr_nrings].set(static_cast<float>(m.nedges() - matoms + 1));
    descriptor[iwdescr_nnsssrng].set(static_cast<float>(m.non_sssr_rings()));
    descriptor[iwdescr_rhacnt].set(static_cast<float>(counts.ring_heteroatom_count));
    descriptor[iwdescr_rhaf].set(static_cast<float>(counts.ring_heteroatom_count) /
                                 static_cast<float>(counts.ring_atom_count));
    descriptor[iwdescr_frafus].set(static_cast<float>(counts.fused_ring_atom_count) /
                                   static_cast<float>(counts.ring_atom_count));
  } else {
    descriptor[iwdescr_nrings].set(0.0f);
    descriptor[iwdescr_nnsssrng].set(0.0f);
    descriptor[iwdescr_rhacnt].set(0.0f);
  }

  if (counts.heavy_atom_count > 0) {
    descriptor[iwdescr_rngatmf].set(static_cast<float>(counts.ring_atom_count) /
                                    static_cast<float>(counts.heavy_atom_count));
  } else {
    descriptor[iwdescr_rngatmf].set(0.0f);
  }

  return 1;
}

int
IWDescr::IWDescrImpl::StoreAromaticDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Aromatic atom descriptors.
  (void)m;
  const int matoms = data.matoms;

  descriptor[iwdescr_aroma].set(static_cast<float>(counts.aromatic_atom_count));
  descriptor[iwdescr_aromha].set(static_cast<float>(counts.aromatic_heteroatom_count));
  descriptor[iwdescr_aromc].set(
      static_cast<float>(counts.aromatic_atom_count - counts.aromatic_heteroatom_count));
  descriptor[iwdescr_aliphc].set(static_cast<float>(counts.aliphatic_carbon_count));

  if (counts.ring_atom_count) {
    descriptor[iwdescr_fraromha].set(static_cast<float>(counts.aromatic_heteroatom_count) /
                                     static_cast<float>(counts.ring_atom_count));
    descriptor[iwdescr_aromdens].set(static_cast<float>(counts.aromatic_atom_count) /
                                     static_cast<float>(matoms));
  } else {
    descriptor[iwdescr_aromdens].set(0.0f);
  }

  return 1;
}


int
IWDescr::IWDescrImpl::StorePiElectronDescriptors(Molecule& m, PerMoleculeData& data,
                                                  const TopologicalDescriptorCounts& counts) {
  // Pulled out of StoreTopologicalCounts so the pi-electron/electron-rich
  // descriptor block is independently testable and easier to compare against
  // the legacy compute_topological_descriptors() implementation.
  (void)m;

  const int matoms = data.matoms;
  if (matoms == 0) {
    return 0;
  }

  descriptor[iwdescr_atmpiele].set(static_cast<float>(counts.atoms_with_pi_electrons));
  descriptor[iwdescr_fratmpie].set(static_cast<float>(counts.atoms_with_pi_electrons) /
                                   static_cast<float>(matoms));
  descriptor[iwdescr_unsatura].set(static_cast<float>(counts.unsaturation));
  descriptor[iwdescr_funsatura].set(static_cast<float>(counts.unsaturation) /
                                    static_cast<float>(matoms));
  descriptor[iwdescr_erichsct].set(static_cast<float>(counts.electron_rich_sections));
  descriptor[iwdescr_aiercsct].set(
      static_cast<float>(counts.atoms_in_electron_rich_sections));
  descriptor[iwdescr_lercsct].set(
      static_cast<float>(counts.largest_electron_rich_section));
  descriptor[iwdescr_faiercst].set(
      static_cast<float>(counts.atoms_in_electron_rich_sections) / static_cast<float>(matoms));

  return 1;
}

int
IWDescr::IWDescrImpl::StoreConnectivityDescriptors(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {
  // Split from StoreTopologicalCounts(). These are pure connectivity-derived
  // descriptor assignments and do not depend on output policy.
  const int matoms = data.matoms;

  descriptor[iwdescr_platt].set(static_cast<float>(counts.total_connectivity));
  descriptor[iwdescr_avcon].set(static_cast<float>(counts.total_connectivity) /
                                static_cast<float>(matoms));
  if (matoms - counts.aromatic_atom_count > 0) {
    descriptor[iwdescr_avalcon].set(static_cast<float>(counts.total_aliphatic_connectivity) /
                                    static_cast<float>(matoms - counts.aromatic_atom_count));
  }

  if (matoms - counts.ring_atom_count > 0) {
    descriptor[iwdescr_avchcon].set(static_cast<float>(counts.total_chain_connectivity) /
                                    static_cast<float>(matoms - counts.ring_atom_count));
  }

  return 1;
}

int
IWDescr::IWDescrImpl::StoreTopologicalCounts(
    Molecule& m, PerMoleculeData& data, const TopologicalDescriptorCounts& counts) {

  // cerr << "StoreTopologicalCounts calling StoreBasicAtomDescriptors\n";
  StoreBasicAtomDescriptors(m, data, counts);

  StoreCoreTopologicalCountDescriptors(m, data, counts);

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeHalogenAttachmentDescriptors(Molecule& m,
                                                          PerMoleculeData& data,
                                                          int halogen_count) {
  if (halogen_count <= 1) {
    descriptor[iwdescr_halogena].set(static_cast<float>(halogen_count));
    return 1;
  }

  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const Atom* const* atom = data.atoms();

  int attached_to_halogen = 0;
  for (int i = 0; i < matoms; ++i) {
    if (z[i] == 9 || z[i] == 17 || z[i] == 35 || z[i] == 53) {
      continue;
    }

    if (ncon[i] == 1) {
      continue;
    }

    const Atom* ai = atom[i];
    for (int j = 0; j < ncon[i]; ++j) {
      const atom_number_t k = ai->other(i, j);
      if (z[k] == 9 || z[k] == 17 || z[k] == 35 || z[k] == 53) {
        ++attached_to_halogen;
        break;
      }
    }
  }

  descriptor[iwdescr_halogena].set(static_cast<float>(attached_to_halogen));

  return 1;
}

namespace {

int
RadhaEntropyChain(Molecule& m, atom_number_t zatom, const int* ncon,
                  const int* ring_membership, const Atom* const* atom,
                  int* already_done) {
  const Atom* ai = atom[zatom];

  already_done[zatom] = 1;

  int rc = 0;
  for (int i = 0; i < ncon[zatom]; ++i) {
    const Bond* b = ai->item(i);
    const atom_number_t j = b->other(zatom);

    if (already_done[j]) {
      continue;
    }

    if (ring_membership[j]) {
      continue;
    }

    if (ncon[j] > 2) {
      ++rc;
    } else if (ncon[j] == 1) {
      continue;
    } else {
      ++rc;
      rc += RadhaEntropyChain(m, j, ncon, ring_membership, atom, already_done);
    }
  }

  return rc;
}

}  // namespace

int
IWDescr::IWDescrImpl::ComputeRadhaEntropyDescriptors(Molecule& m,
                                                     PerMoleculeData& data) {
  const int matoms = data.matoms;
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const Atom* const* atom = data.atoms();

  data.ResetAlreadyDone();
  int* already_done = data.already_done_data();

  int longest_chain = 0;
  int total_chain_atoms = 0;
  int nchains = 0;
  int non_ring_atoms = 0;

  double radha_entropy = 1.0;

  for (int i = 0; i < matoms; ++i) {
    if (ring_membership[i]) {
      continue;
    }

    ++non_ring_atoms;

    if (already_done[i] || ncon[i] != 2) {
      continue;
    }

    const int chain_length = RadhaEntropyChain(m, i, ncon, ring_membership, atom,
                                               already_done);

    if (chain_length == 0) {
      continue;
    }

    if (chain_length == 1) {
      radha_entropy += 3.0;
    } else if (chain_length == 2) {
      radha_entropy += 9.0;
    } else {
      radha_entropy += std::pow(3.0, static_cast<double>(chain_length));
    }

    ++nchains;
    total_chain_atoms += chain_length;
    if (chain_length > longest_chain) {
      longest_chain = chain_length;
    }
  }

  descriptor[iwdescr_nflxchn].set(static_cast<float>(nchains));
  descriptor[iwdescr_lflxchn].set(static_cast<float>(longest_chain));
  descriptor[iwdescr_atflxchn].set(static_cast<float>(total_chain_atoms));

  if (nchains) {
    descriptor[iwdescr_avflxchn].set(static_cast<float>(total_chain_atoms) /
                                     static_cast<float>(nchains));
  } else {
    descriptor[iwdescr_avflxchn].set(0.0f);
  }

  descriptor[iwdescr_faflxchn].set(static_cast<float>(total_chain_atoms) /
                                   static_cast<float>(matoms));

  if (non_ring_atoms) {
    descriptor[iwdescr_fnflxchn].set(static_cast<float>(total_chain_atoms) /
                                     static_cast<float>(non_ring_atoms));
  } else {
    descriptor[iwdescr_fnflxchn].set(0.0f);
  }

  descriptor[iwdescr_rkentrpy].set(static_cast<float>(std::log(radha_entropy)));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeTopologicalDescriptors(Molecule& m, PerMoleculeData& data) {
  // Core topological counts are computed first. Optional topological descriptor
  // families then run with a clean already_done array, preserving the legacy
  // sequencing from compute_topological_descriptors().
  data.ResetAlreadyDone();

  // cerr << "Going to ComputeCoreTopologicalDescriptors\n";
  ComputeCoreTopologicalDescriptors(m, data, data.already_done_data());

  ComputeOptionalTopologicalDescriptors(m, data);

  compute_polar_bond_descriptors(m, data.atomic_numbers(), data.connections());

  if (descriptors_to_compute.complexity_descriptors) {
    ComputeSpiroFusionDescriptors(m, data);
  }

  if (descriptors_to_compute.ramey_descriptors) {
    ComputeRameyDescriptors(m, data);
  }

  return 1;
}

int
IWDescr::IWDescrImpl::compute_polar_bond_descriptors(Molecule& m, const atomic_number_t* z, const int* ncon) {
  int nb = m.nedges();

  int polar_bond_count = 0;
  int aromatic_polar_bond_count = 0;
  int aromatic_non_polar_bond_count = 0;
  int polar_multiple_bonds = 0;
  int di_vinyl_bonds = 0;  // matches the single bond in  *=*-*=*, not in an aromatic ring

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (z[a1] != z[a2]) {
      polar_bond_count++;
    }

    int isar = m.in_same_aromatic_ring(a1, a2);

    if (isar) {  // in same aromatic ring(s)
      if (z[a1] == z[a2]) {
        aromatic_non_polar_bond_count++;
      } else {
        aromatic_polar_bond_count++;
      }
    } else if (z[a1] == z[a2]) {  // a non-polar bond
      ;
    } else if (!b->is_single_bond()) {
      polar_multiple_bonds++;
    }

    if (!isar && b->is_single_bond() && (ncon[a1] < m.nbonds(a1)) &&
        (ncon[a2] < m.nbonds(a2))) {
      di_vinyl_bonds++;
    }
  }

  descriptor[iwdescr_pbcount].set(static_cast<float>(polar_bond_count));
  if (nb > 0) {
    descriptor[iwdescr_frpbond].set(static_cast<float>(polar_bond_count) /
                                    static_cast<float>(nb));
  }
  descriptor[iwdescr_nonpbond].set(static_cast<float>(nb - polar_bond_count));
  descriptor[iwdescr_pbarom].set(static_cast<float>(aromatic_polar_bond_count));
  descriptor[iwdescr_npbarom].set(static_cast<float>(aromatic_non_polar_bond_count));
  descriptor[iwdescr_pbunset].set(static_cast<float>(polar_multiple_bonds));
  descriptor[iwdescr_dvinylb].set(static_cast<float>(di_vinyl_bonds));

  return 1;
}

int
IWDescr::IWDescrImpl::compute_double_bond_substitution(Molecule& m, const atomic_number_t* z, const int* ncon) {
  int nb = m.nedges();

  int number_doubly_bonded_carbons = 0;

  int total_double_bond_substituents = 0;

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    if (!b->is_double_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (m.in_same_aromatic_ring(a1, a2)) {
      continue;
    }

    if (6 == z[a1]) {
      number_doubly_bonded_carbons++;
      total_double_bond_substituents += ncon[a1] - 1;
    }

    if (6 == z[a2]) {
      number_doubly_bonded_carbons++;
      total_double_bond_substituents += ncon[a2] - 1;
    }
  }

  descriptor[iwdescr_numcdb].set(static_cast<float>(number_doubly_bonded_carbons));
  descriptor[iwdescr_totdbsub].set(static_cast<float>(total_double_bond_substituents));

  if (number_doubly_bonded_carbons) {
    descriptor[iwdescr_avcdbsub].set(static_cast<float>(total_double_bond_substituents) /
                                     static_cast<float>(number_doubly_bonded_carbons));
  } else {
    descriptor[iwdescr_avcdbsub].set(0.0f);
  }

  return 1;
}


int
IWDescr::IWDescrImpl::ComputeOptionalTopologicalDescriptors(Molecule& m,
                                                           PerMoleculeData& data) {
  // Preserve the legacy convention that the optional tail families get a clean
  // already_done array after the electron-rich-section traversal used by the
  // core topological pass.
  // cerr << "ComputeOptionalTopologicalDescriptors\n";
  data.ResetAlreadyDone();

  // Keep these calls in the legacy order unless there is a verified reason to change
  // them. Some optional families use and reset `already_done`.
  ComputeOptionalChainDescriptors(m, data);

  ComputeOptionalSpecificAndHbondDescriptors(m, data);

  if (m.nrings() == 0) {
    ZeroAllRingRelatedDescriptors();
  } else {
    ComputeOptionalRingAndDistanceDescriptors(m, data);
  }

  ComputeOptionalLogPAndConjugationDescriptors(m, data);

  return 1;
}

void
IWDescr::IWDescrImpl::ZeroAllRingRelatedDescriptors() {
  static constexpr float kZero = 0.0f;

  descriptor[iwdescr_ar5].set(kZero);
  descriptor[iwdescr_al5].set(kZero);
  descriptor[iwdescr_ar6].set(kZero);
  descriptor[iwdescr_al6].set(kZero);
  descriptor[iwdescr_npfsdsys].set(kZero);
  descriptor[iwdescr_rnginpfs].set(kZero);
  descriptor[iwdescr_lgplnfsy].set(kZero);
  descriptor[iwdescr_htrcpfsy].set(kZero);
  descriptor[iwdescr_mxhtpfsy].set(kZero);

  return;
}

void
IWDescr::IWDescrImpl::ZeroRunFusionDescriptors() {
  static constexpr float kZero = 0.0f;

  descriptor[iwdescr_npfsdsys].set(kZero);
  descriptor[iwdescr_rnginpfs].set(kZero);
  descriptor[iwdescr_lgplnfsy].set(kZero);
  descriptor[iwdescr_htrcpfsy].set(kZero);
  descriptor[iwdescr_mxhtpfsy].set(kZero);
}

int
IWDescr::IWDescrImpl::ComputeOptionalChainDescriptors(Molecule& m,
                                                      PerMoleculeData& data) {
  // cerr << "in ComputeOptionalChainDescriptors \n";
  ComputeOptionalLongCarbonChainDescriptors(m, data);

  // cerr << "Going to ComputeOptionalSaturatedChainDescriptors\n";
  return ComputeOptionalSaturatedChainDescriptors(m, data);
}

int
IWDescr::IWDescrImpl::ComputeOptionalLongCarbonChainDescriptors(Molecule& m,
                                                               PerMoleculeData& data) {
  if (descriptors_to_compute.long_carbon_chains) {
    ComputeLongCarbonChainDescriptors(m, data, data.already_done_data());
    data.ResetAlreadyDone();
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeOptionalSaturatedChainDescriptors(Molecule& m,
                                                              PerMoleculeData& data) {
  if (descriptors_to_compute.saturated_chains) {
    ComputeSaturatedChainDescriptors(m, data, data.already_done_data());
    data.ResetAlreadyDone();
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeOptionalSpecificAndHbondDescriptors(
    Molecule& m, PerMoleculeData& data) {
  ComputeOptionalSpecificGroupDescriptors(m, data);

  return ComputeOptionalSimpleHbondDescriptors(m, data);
}

int
IWDescr::IWDescrImpl::ComputeOptionalSpecificGroupDescriptors(
    Molecule& m, PerMoleculeData& data) {

  if (descriptors_to_compute.specific_groups &&
      ! ComputeSpecificGroupDescriptors(m, data, data.already_done_data())) {
    return 0;
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeOptionalSimpleHbondDescriptors(
    Molecule& m, PerMoleculeData& data) {
  if (descriptors_to_compute.simple_hbond_descriptors &&
      ! ComputeSimpleHbondDescriptors(m, data)) {
    return 0;
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeOptionalRingAndDistanceDescriptors(Molecule& m,
                                                                PerMoleculeData& data) {
  // Ring substitution descriptors are kept after the core ring descriptors to
  // preserve the legacy optional-tail ordering.
  ComputeOptionalCoreRingDescriptors(m, data);

  return ComputeOptionalRingSubstitutionDescriptors(m, data);
}

int
IWDescr::IWDescrImpl::ComputeOptionalCoreRingDescriptors(Molecule& m,
                                                         PerMoleculeData& data) {
  ComputeRingDescriptors(m, data);

  compute_4_connected_carbon_stuff(m, data.atomic_numbers(), data.connections(), data.ring_membership_data());

  if (descriptors_to_compute.bonds_between_rings) {
    ComputeBetweenRingDescriptors(m, data);
  }

  if (descriptors_to_compute.ring_fusion_descriptors && m.nrings() > 1) {
    compute_ring_fusion_descriptors(m, data.connections(), data.ring_membership_data());
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeOptionalRingSubstitutionDescriptors(
    Molecule& m, PerMoleculeData& data) {
  return ComputeRingSubstitutionDescriptors(m, data);
}

#define MAX_RING_SIZE 9

/*
  We count the number of types of ring fusions.
  These are laid out in the descriptor array as:

  5l5l
  5l5r
  5r5r
  5l6l
  5l6r
  5r6l
  5r6r
  6l6l
  6l6r
  6r6r

  Where l means aLiphatic and r means aRomatic
*/

int
IWDescr::IWDescrImpl::compute_ring_fusion_descriptors(Molecule& m, const int* ncon,
                                const int* ring_membership) {
  int ar5 = 0;
  int ar6 = 0;
  int al5 = 0;
  int al6 = 0;

  int fsdrng5l5l = 0;
  int fsdrng5l5r = 0;
  int fsdrng5r5r = 0;
  int fsdrng5l6l = 0;
  int fsdrng5l6r = 0;
  int fsdrng5r6l = 0;
  int fsdrng5r6r = 0;
  int fsdrng6l6l = 0;
  int fsdrng6l6r = 0;
  int fsdrng6r6r = 0;

  m.compute_aromaticity_if_needed();

  int nr = m.nrings();

  if (nr < 2) {
    return 1;
  }

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    int rsi = ri->number_elements();

    if (rsi < 5 || rsi > 6) {
      continue;
    }

    int n = ri->fused_ring_neighbours();

    if (0 == n) {  // isolated ring
      if (ri->is_aromatic()) {
        if (6 == rsi) {
          ar6++;
        } else if (5 == rsi) {
          ar5++;
        }
      } else {
        if (6 == rsi) {
          al6++;
        } else if (5 == rsi) {
          al5++;
        }
      }
      continue;
    }

    if (fsdrng_descriptors_consider_just_two_ring_systems && n > 1) {
      continue;
    }

    for (int j = 0; j < n; j++) {
      const Ring* rj = ri->fused_neighbour(j);

      if (rj->ring_number() < ri->ring_number()) {
        continue;
      }

      int rsj = rj->number_elements();

      if (rsj < 5 || rsj > 6) {
        continue;
      }

      //    cerr << "lbswar " << ri->largest_number_of_bonds_shared_with_another_ring() <<
      //    '\n'; cerr << "Shared " << ri->compute_bonds_shared_with(*rj) << '\n';

      if (ri->largest_number_of_bonds_shared_with_another_ring() > 1 &&  // strongly fused
          ri->compute_bonds_shared_with(*rj) > 1) {
        continue;
      }

      if (fsdrng_descriptors_consider_just_two_ring_systems &&
          rj->fused_ring_neighbours() > 1) {
        continue;
      }

      if (ri->is_aromatic() && rj->is_aromatic()) {
        if (5 == rsi && 5 == rsj) {
          fsdrng5r5r++;
        } else if (6 == rsi && 6 == rsj) {
          fsdrng6r6r++;
        } else {
          fsdrng5r6r++;
        }
      } else if (!ri->is_aromatic() && !rj->is_aromatic()) {  // both aliphatic
        if (5 == rsi && 5 == rsj) {
          fsdrng5l5l++;
        } else if (6 == rsi && 6 == rsj) {
          fsdrng6l6l++;
        } else {
          fsdrng5l6l++;
        }
      } else  // rings have different aromaticity
      {
        if (5 == rsi && 5 == rsj) {  // both of size 5
          fsdrng5l5r++;
        } else if (6 == rsi && 6 == rsj) {  // both of size 6
          fsdrng6l6r++;
        } else if (ri->is_aromatic()) {  // therefore rj is aliphatic
          if (5 == rsi) {
            fsdrng5r6l++;
          } else {
            fsdrng5l6r++;
          }
        } else  // rj is aromatic and ri is aliphatic
        {
          if (5 == rsi) {
            fsdrng5l6r++;
          } else {
            fsdrng5r6l++;
          }
        }
      }
    }
  }

  descriptor[iwdescr_ar5].set(ar5);
  descriptor[iwdescr_ar6].set(ar6);
  descriptor[iwdescr_al5].set(al5);
  descriptor[iwdescr_al6].set(al6);

  descriptor[iwdescr_fsdrng5l5l].set(fsdrng5l5l);
  descriptor[iwdescr_fsdrng5l5r].set(fsdrng5l5r);
  descriptor[iwdescr_fsdrng5r5r].set(fsdrng5r5r);
  descriptor[iwdescr_fsdrng5l6l].set(fsdrng5l6l);
  descriptor[iwdescr_fsdrng5l6r].set(fsdrng5l6r);
  descriptor[iwdescr_fsdrng5r6l].set(fsdrng5r6l);
  descriptor[iwdescr_fsdrng5r6r].set(fsdrng5r6r);
  descriptor[iwdescr_fsdrng6l6l].set(fsdrng6l6l);
  descriptor[iwdescr_fsdrng6l6r].set(fsdrng6l6r);
  descriptor[iwdescr_fsdrng6r6r].set(fsdrng6r6r);

  descriptor[iwdescr_fsdrngarar].set(fsdrng5r5r + fsdrng5r6r + fsdrng6r6r);
  descriptor[iwdescr_fsdrngalar].set(fsdrng5l5r + fsdrng5l6r + fsdrng6l6r);
  descriptor[iwdescr_fsdrngalal].set(fsdrng5l5l + fsdrng5l6l + fsdrng6l6l);

  return 1;
}


static int
bonds_to_nearest_ring(const Molecule& m, atom_number_t previous_atom,
                      atom_number_t current_atom, const Atom* const* atom,
                      const int* fsid) {
  const Atom* a = atom[current_atom];

  int acon = a->ncon();

  int matoms = m.natoms();

  int rc = matoms;  // we look for the shortest path

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(current_atom, i);

    if (previous_atom == j) {
      continue;
    }

    if (fsid[j] > 0) {
      return 1;
    }

    int tmp = bonds_to_nearest_ring(m, current_atom, j, atom, fsid);

    if (tmp <= 0) {
      continue;
    }

    tmp++;

    if (tmp < rc) {
      rc = tmp;
    }
  }

  if (matoms == rc) {  // didn't find anything
    return -1;
  }

  return rc;
}


int
IWDescr::IWDescrImpl::ComputeBetweenRingDescriptors(Molecule& m, 
                                                    PerMoleculeData& data) {
  const int* ncon = data.connections();
  const Atom** atom = data.atoms();

  int nr = m.nrings();


  if (nr < 2) {  // we deal with between ring information
    descriptor[iwdescr_bbr1].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr2].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr3].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr4].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr5].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr6].set(static_cast<float>(0.0));

    return 1;
  }

  int matoms = m.natoms();

  int* fsid = new int[matoms];
  std::unique_ptr<int[]> free_fsid(fsid);

  if (1 == m.label_atoms_by_ring_system(fsid)) {  // more than 1 ring, but 1 ring system
    descriptor[iwdescr_bbr1].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr2].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr3].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr4].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr5].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr6].set(static_cast<float>(0.0));

    return 1;
  }

  extending_resizable_array<int> bonds_between_rings;

  for (int i = 0; i < matoms; i++) {
    if (0 == fsid[i]) {  // not in a ring
      continue;
    }

    if (2 == ncon[i]) {
      continue;
    }

    const Atom* ai = atom[i];

    for (int j = 0; j < ncon[i]; j++) {
      const Bond* b = ai->item(j);

      if (b->nrings()) {
        continue;
      }

      atom_number_t k = b->other(i);

      if (fsid[k]) {
        bonds_between_rings[1]++;
      } else {
        int tmp = bonds_to_nearest_ring(m, i, k, atom, fsid);
        if (tmp > 0) {
          bonds_between_rings[tmp + 1]++;
        }
      }
    }
  }

  descriptor[iwdescr_bbr1].set(static_cast<float>(bonds_between_rings[1] / 2));
  descriptor[iwdescr_bbr2].set(static_cast<float>(bonds_between_rings[2] / 2));
  descriptor[iwdescr_bbr3].set(static_cast<float>(bonds_between_rings[3] / 2));
  descriptor[iwdescr_bbr4].set(static_cast<float>(bonds_between_rings[4] / 2));
  descriptor[iwdescr_bbr5].set(static_cast<float>(bonds_between_rings[5] / 2));
  descriptor[iwdescr_bbr6].set(static_cast<float>(bonds_between_rings[6] / 2));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeOptionalLogPAndConjugationDescriptors(
    Molecule& m, PerMoleculeData& data) {
  ComputeOptionalLogPOnlyDescriptors(m, data);

  return ComputeOptionalConjugationDescriptors(m, data);
}

int
IWDescr::IWDescrImpl::ComputeOptionalLogPOnlyDescriptors(Molecule& m,
                                                         PerMoleculeData& data) {
  return ComputeLogPDescriptors(m, data);
}

int
IWDescr::IWDescrImpl::ComputeOptionalConjugationDescriptors(Molecule& m,
                                                            PerMoleculeData& data) {
  return ComputeExtendedConjugationDescriptors(m, data);
}

int
IWDescr::IWDescrImpl::ComputeLongCarbonChainDescriptors(Molecule& m, PerMoleculeData& data,
                                                        int* already_done) {
  // Migrated from legacy ComputeLongCarbonChains(). This descriptor family is
  // pure computation: it consumes PerMoleculeData and writes IWDescr-owned
  // Descriptor values.
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const int matoms = data.matoms;

  auto traverse = [&](auto&& self, atom_number_t zatom, int max_ncon) -> int {
    already_done[zatom] = 1;
    int rc = 1;

    for (const Bond* b : m[zatom]) {
      if (! b->is_single_bond()) {
        return rc;
      }

      const atom_number_t o = b->other(zatom);
      if (already_done[o]) {
        continue;
      }
      if (z[o] != 6) {
        continue;
      }
      if (ring_membership[o]) {
        continue;
      }
      if (ncon[o] > max_ncon) {
        continue;
      }
      if (! m.saturated(o)) {
        continue;
      }

      if (ncon[o] == 1) {
        already_done[o] = 1;
        ++rc;
      } else {
        rc += self(self, o, max_ncon);
      }
    }

    return rc;
  };

  Accumulator_Int<uint32_t> acc2;

  int max_ncon = 2;
  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }
    if (z[i] != 6 || ring_membership[i]) {
      continue;
    }
    if (! m.saturated(i)) {
      continue;
    }
    if (ncon[i] > max_ncon) {
      continue;
    }

    acc2.extra(traverse(traverse, i, max_ncon));
  }

  if (acc2.n()) {
    descriptor[iwdescr_mxlencchain2].set(static_cast<float>(acc2.maxval()));
  } else {
    descriptor[iwdescr_mxlencchain2].set(0.0f);
  }

  std::fill_n(already_done, matoms, 0);

  Accumulator_Int<uint32_t> acc3;
  max_ncon = 3;
  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }
    if (z[i] != 6 || ring_membership[i]) {
      continue;
    }
    if (! m.saturated(i)) {
      continue;
    }
    if (ncon[i] > max_ncon) {
      continue;
    }

    acc3.extra(traverse(traverse, i, max_ncon));
  }

  if (acc3.n()) {
    descriptor[iwdescr_mxlencchain3].set(static_cast<float>(acc3.maxval()));
  } else {
    descriptor[iwdescr_mxlencchain3].set(0.0f);
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeSaturatedChainDescriptors(Molecule& m, PerMoleculeData& data,
                                                       int* already_done) {
  // Migrated from legacy ComputeSaturatedChains(). The arbitrary legacy cutoff
  // of three atoms is preserved.
  // cerr << " in ComputeSaturatedChainDescriptors \n";
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const int matoms = data.matoms;

  auto traverse = [&](auto&& self, atom_number_t zatom) -> int {
    int rc = 1;
    already_done[zatom] = 1;

    for (const Bond* b : m[zatom]) {
      const atom_number_t o = b->other(zatom);
      if (ring_membership[o]) {
        continue;
      }
      if (already_done[o]) {
        continue;
      }
      if (ncon[o] != 2) {
        continue;
      }
      if (! m.saturated(o)) {
        continue;
      }

      rc += self(self, o);
    }

    return rc;
  };

  int number_regions = 0;
  int atoms_in_regions = 0;
  int largest_region = 0;

  for (int i = 0; i < matoms; ++i) {
    if (ring_membership[i]) {
      continue;
    }
    if (ncon[i] != 2) {
      continue;
    }
    if (already_done[i]) {
      continue;
    }
    if (! m.saturated(i)) {
      continue;
    }

    const int tmp = traverse(traverse, i);
    if (tmp < 3) {
      continue;
    }

    ++number_regions;
    atoms_in_regions += tmp;
    if (tmp > largest_region) {
      largest_region = tmp;
    }
  }

  if (number_regions == 0) {
//  descriptor[iwdescr_nsatchain].set(0.0f);
//  descriptor[iwdescr_mxsatchain].set(0.0f);
//  descriptor[iwdescr_fsatchain].set(0.0f);
//  return 1;
  }

  descriptor[iwdescr_nsatchain].set(static_cast<float>(number_regions));
  descriptor[iwdescr_mxsatchain].set(static_cast<float>(largest_region));
  descriptor[iwdescr_fsatchain].set(static_cast<float>(atoms_in_regions) /
                                    static_cast<float>(matoms));

  return 1;
}


int
IWDescr::IWDescrImpl::ComputeExtendedConjugationDescriptors(Molecule& m,
                                                           PerMoleculeData& data) {
  const int matoms = data.matoms;
  const Atom** atom = data.atoms();

  std::vector<int> already_done(matoms, 0);

  auto back_to_zero = [&](int from, int to) {
    for (int i = 0; i < matoms; ++i) {
      if (already_done[i] == from) {
        already_done[i] = to;
      }
    }
  };

  auto discern_conjugated_section = [&](auto&& self, atom_number_t zatom,
                                        int flag) -> int {
    already_done[zatom] = flag;

    int rc = 1;
    const Atom* a = atom[zatom];

    for (const Bond* b : *a) {
      atom_number_t j = b->other(zatom);

      if (already_done[j]) {
        continue;
      }

      int pe = 0;
      if (b->is_double_bond() || b->is_triple_bond()) {
        ;  // definitely extend section across these.
      } else if (const_cast<Atom*>(atom[j])->pi_electrons(pe) && pe) {
        ;
      } else {
        continue;
      }

      rc += self(self, j, flag);
    }

    return rc;
  };

  int number_conjugated_sections = 0;
  int max_conjugated_section_size = 0;
  int atoms_in_conjugated_sections = 0;

  // Preserve the legacy call. Molecule computes ring/aromaticity lazily, but
  // this makes the intent explicit and keeps behaviour aligned with the old code.
  (void) m.ring_membership();

  for (int i = 0; i < matoms; ++i) {
    if (already_done[i]) {
      continue;
    }

    int pe = 0;
    if (! const_cast<Atom*>(atom[i])->pi_electrons(pe) || pe == 0) {
      continue;
    }

    const int flag = number_conjugated_sections + 1;
    const int conjugated_section_size = discern_conjugated_section(
        discern_conjugated_section, i, flag);

    // Just an isolated double bond, or a section that did not extend.
    if (conjugated_section_size <= 2) {
      back_to_zero(flag, 0);
      continue;
    }

    atoms_in_conjugated_sections += conjugated_section_size;
    ++number_conjugated_sections;

    if (conjugated_section_size > max_conjugated_section_size) {
      max_conjugated_section_size = conjugated_section_size;
    }
  }

  int carbons_in_conjugated_sections = 0;
  if (number_conjugated_sections) {
    for (int i = 0; i < matoms; ++i) {
      if (already_done[i] > 0 && atom[i]->atomic_number() == 6) {
        ++carbons_in_conjugated_sections;
      }
    }
  }

  descriptor[iwdescr_nconjgsc].set(static_cast<float>(number_conjugated_sections));
  descriptor[iwdescr_atincnjs].set(static_cast<float>(atoms_in_conjugated_sections));
  descriptor[iwdescr_mxcnjscz].set(static_cast<float>(max_conjugated_section_size));
  descriptor[iwdescr_cinconjs].set(static_cast<float>(carbons_in_conjugated_sections));

  return 1;
}

int
IWDescr::IWDescrImpl::Guanidine(Molecule& m, atom_number_t doubly_bonded_nitrogen,
                                int* nitrogens, const Atom** atom,
                                const atomic_number_t* z, const int* ncon,
                                const int* ring_membership) const {
  atom_number_t c = atom[doubly_bonded_nitrogen]->other(doubly_bonded_nitrogen, 0);

  if (3 != ncon[c] || ring_membership[c]) {
    return 0;
  }

  const Atom* carbon = atom[c];

  if (4 != carbon->nbonds()) {
    return 0;
  }

  atom_number_t terminal_atom = INVALID_ATOM_NUMBER;
  atom_number_t linker_atom = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 3; ++i) {
    const Bond* b = carbon->item(i);

    if (b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(c);

    if (7 != z[j] || ring_membership[j]) {
      return 0;
    }

    if (1 == ncon[j] && INVALID_ATOM_NUMBER == terminal_atom) {
      terminal_atom = j;
    } else if (2 == ncon[j] && INVALID_ATOM_NUMBER == linker_atom) {
      linker_atom = j;
    } else {
      return 0;
    }
  }

  if (INVALID_ATOM_NUMBER == terminal_atom || INVALID_ATOM_NUMBER == linker_atom) {
    return 0;
  }

  nitrogens[doubly_bonded_nitrogen] = 1;
  nitrogens[terminal_atom] = 1;
  nitrogens[linker_atom] = 1;

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeSpecificGroupDescriptors(Molecule& m, PerMoleculeData& data,
                                                      int* already_done) {
  // Migrated from legacy compute_amine_count() and
  // compute_pyridine_pyrrole(). These descriptors are calculation-only and
  // use PerMoleculeData plus IWDescr-owned descriptor storage.
  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const Atom** atom = data.atoms();

  int amine_count = 0;

  // First perceive guanidines: CNC(=N)N. Mark all guanidine nitrogens in
  // already_done so they are not also counted as ordinary amines.
  for (int i = 0; i < matoms; ++i) {
    if (7 == z[i] && 1 == ncon[i] && 2 == atom[i]->nbonds()) {
      amine_count += Guanidine(m, i, already_done, atom, z, ncon, ring_membership);
    }
  }

  for (int i = 0; i < matoms; ++i) {
    if (7 == z[i] && 0 == already_done[i] && m.hcount(i)) {
      ++amine_count;
    }
  }

  descriptor[iwdescr_amine].set(static_cast<float>(amine_count));

  int pyridine = 0;
  int pyrrole = 0;

  for (int i = 0; i < matoms; ++i) {
    if (7 != z[i]) {
      continue;
    }

    if (0 == ring_membership[i]) {
      continue;
    }

    if (! m.is_aromatic(i)) {
      continue;
    }

    if (3 == ncon[i]) {
      continue;
    }

    if (0 == m.implicit_hydrogens(i)) {
      ++pyridine;
    } else {
      ++pyrrole;
    }
  }

  descriptor[iwdescr_pyridine].set(static_cast<float>(pyridine));
  descriptor[iwdescr_pyrrole].set(static_cast<float>(pyrrole));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeSimpleHbondDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_hbond_descriptors(m, z, ncon, atom).
  // This is pure descriptor computation and writes only IWDescr-owned descriptor
  // values, so it belongs behind the IWDescrImpl::Process() boundary.
  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const Atom** atom = data.atoms();

  double acceptor_score = 0.0;
  double donor_score = 0.0;
  double dual_score = 0.0;

  for (int i = 0; i < matoms; ++i) {
    Atom* ai = const_cast<Atom*>(atom[i]);

    const int nbi = ai->nbonds();
    const int hci = ai->implicit_hydrogens();

    if (z[i] == 8 && ncon[i] == 1 && nbi == 2) {  // carbonyl
      acceptor_score += 1.0;
    } else if (z[i] == 8 && ncon[i] == 1 && nbi == 1) {  // hydroxy
      dual_score += 1.0;
    } else if (z[i] == 8 && ncon[i] == 2) {  // ether
      acceptor_score += 0.5;
    } else if (z[i] == 7 && hci == 2) {  // primary amine
      donor_score += 1.0;
    } else if (z[i] == 7 && hci == 1) {  // secondary amine
      donor_score += 1.0;
    } else if (z[i] == 7 && ncon[i] == 2 && nbi == 3) {  // =N-
      acceptor_score += 1.0;
    } else if (z[i] == 7 && hci == 0) {  // tertiary amine
      acceptor_score += 0.2;
    } else if (z[i] == 16 && hci == 1) {
      donor_score += 1.0;
    } else if (z[i] == 16 && ncon[i] == 2) {
      acceptor_score += 0.2;
    }
  }

  descriptor[iwdescr_hacts].set(static_cast<float>(acceptor_score));
  descriptor[iwdescr_hdons].set(static_cast<float>(donor_score));
  descriptor[iwdescr_hduals].set(static_cast<float>(dual_score));

  return 1;
}

int
IWDescr::IWDescrImpl::AtomsNotAlreadyMarked(const Ring& r, int* atom_already_done) const {
  int rc = 0;

  const int ring_size = r.number_elements();
  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t j = r[i];

    if (atom_already_done[j]) {
      continue;
    }

    atom_already_done[j] = 1;
    ++rc;
  }

  return rc;
}

#ifdef NO_LONGER_USED
unstable WRT non-sssr rings
int
IWDescr::IWDescrImpl::ComputeExocyclicBonds(Molecule& m, const Ring& r,
                                            const atomic_number_t* z,
                                            const int* ncon,
                                            int& double_bond_attachments,
                                            int& singly_connected_attachments,
                                            int& singly_connected_heteroatoms,
                                            int& singly_connected_donors) const {
  int rc = 0;

  const int ring_size = r.number_elements();
  for (atom_number_t j : r) {
    if (ncon[j] == 2) {
      continue;
    }

    for (const Bond* b : m[j]) {
      if (b->nrings()) {
        continue;
      }

      ++rc;

      const atom_number_t o = b->other(j);
      if (ncon[o] != 1) {
        continue;
      }

      if (b->is_double_bond()) {
        ++double_bond_attachments;
      }

      ++singly_connected_attachments;
      if (z[o] == 6) {
        continue;
      }

      ++singly_connected_heteroatoms;
      if (m.hcount(o)) {
        ++singly_connected_donors;
      }
      // Cannot break here because there might be two attachments in an aliphatic ring
    }
  }

  return rc;
}
#endif

int
IWDescr::IWDescrImpl::JustTerminalGroupsOutsideRing(const Molecule& m, const Ring& r,
                                                    const int* ncon,
                                                    atom_number_t a) const {
  const Atom* ai = m.atomi(a);

  for (int i = 0; i < ncon[a]; ++i) {
    const atom_number_t j = ai->other(a, i);

    if (ncon[j] == 1) {
      continue;
    }

    if (r.contains(j)) {
      continue;
    }

    if (ncon[j] > 1) {
      return 0;
    }
  }

  return 1;
}

float
IWDescr::IWDescrImpl::ComputeRingIsolation(const Molecule& m, const Ring& r,
                                           const int* ncon) const {
  const int ring_size = r.number_elements();

  int branches = 0;
  int terminal_groups_found_here = 0;

  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t j = r[i];
    if (ncon[j] == 2) {
      continue;
    }

    ++branches;

    if (JustTerminalGroupsOutsideRing(m, r, ncon, j)) {
      ++terminal_groups_found_here;
    }
  }

  branches -= terminal_groups_found_here;

  if (branches == 0) {
    return 2.0f;
  }

  return 1.0f / static_cast<float>(branches);
}

int
IWDescr::IWDescrImpl::HeteroatomsInRing(const Ring* r, const atomic_number_t* z) const {
  int rc = 0;

  const int ring_size = r->number_elements();
  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t j = r->item(i);
    if (z[j] != 6) {
      ++rc;
    }
  }

  return rc;
}

int
IWDescr::IWDescrImpl::ComputeRingDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_ring_descriptors(m, z, ncon). This method now
  // owns the main ring-size, ring-system, ring-isolation, heterocycle, and
  // exocyclic-bond descriptor calculations. More specialised ring families
  // remain split into their own IWDescrImpl methods.
  constexpr int kMaxRingSize = 9;  // Must remain aligned with iwdescr_nrings3..8.

  const int nr = m.nrings();
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const int* is_aromatic_atom = data.aromaticity();

  descriptor[iwdescr_trmnlrng].set(0.0f);
  descriptor[iwdescr_intrnlrng].set(0.0f);
  descriptor[iwdescr_rng2spch].set(0.0f);
  descriptor[iwdescr_rng2bridge].set(0.0f);

  if (nr < 2) {
    descriptor[iwdescr_fsdrng5l5l].set(0.0f);
    descriptor[iwdescr_fsdrng5l5r].set(0.0f);
    descriptor[iwdescr_fsdrng5r5r].set(0.0f);
    descriptor[iwdescr_fsdrng5l6l].set(0.0f);
    descriptor[iwdescr_fsdrng5l6r].set(0.0f);
    descriptor[iwdescr_fsdrng5r6l].set(0.0f);
    descriptor[iwdescr_fsdrng5r6r].set(0.0f);
    descriptor[iwdescr_fsdrng6r6r].set(0.0f);
    descriptor[iwdescr_fsdrng6l6r].set(0.0f);
    descriptor[iwdescr_fsdrng6l6l].set(0.0f);

    descriptor[iwdescr_fsdrngarar].set(0.0f);
    descriptor[iwdescr_fsdrngalar].set(0.0f);
    descriptor[iwdescr_fsdrngalal].set(0.0f);

    descriptor[iwdescr_nsfsdsys].set(0.0f);
    descriptor[iwdescr_rnginsfs].set(0.0f);
    descriptor[iwdescr_lgstrfsy].set(0.0f);
    descriptor[iwdescr_htrcsfsy].set(0.0f);
    descriptor[iwdescr_mxhtsfsy].set(0.0f);
  }

  if (nr == 0) {
    descriptor[iwdescr_mhr].set(0.0f);
    descriptor[iwdescr_mxhrf].set(0.0f);
    descriptor[iwdescr_mnhrf].set(0.0f);
    descriptor[iwdescr_lrsysz].set(0.0f);
    descriptor[iwdescr_lrsz].set(0.0f);
    descriptor[iwdescr_rng7atoms].set(0.0f);
    descriptor[iwdescr_mars].set(0.0f);
    descriptor[iwdescr_ringsys].set(0.0f);
    descriptor[iwdescr_ringisol].set(0.0f);
    descriptor[iwdescr_isolrc].set(0.0f);
    descriptor[iwdescr_isolhtrc].set(0.0f);
    descriptor[iwdescr_arring].set(0.0f);
    descriptor[iwdescr_alring].set(0.0f);
    descriptor[iwdescr_excybond].set(0.0f);
    descriptor[iwdescr_excydbond].set(0.0f);
    descriptor[iwdescr_excydscon].set(0.0f);
    descriptor[iwdescr_excydsconh].set(0.0f);
    descriptor[iwdescr_excydscondon].set(0.0f);

    descriptor[iwdescr_nrings3].set(0.0f);
    descriptor[iwdescr_nrings4].set(0.0f);
    descriptor[iwdescr_nrings5].set(0.0f);
    descriptor[iwdescr_nrings6].set(0.0f);
    descriptor[iwdescr_nrings7].set(0.0f);
    descriptor[iwdescr_nrings8].set(0.0f);

    descriptor[iwdescr_rsarom1].set(0.0f);
    descriptor[iwdescr_rsarom2].set(0.0f);
    descriptor[iwdescr_rsarom3].set(0.0f);

    descriptor[iwdescr_rsaliph1].set(0.0f);
    descriptor[iwdescr_rsaliph2].set(0.0f);
    descriptor[iwdescr_rsaliph3].set(0.0f);
    descriptor[iwdescr_rsaliph4].set(0.0f);

    descriptor[iwdescr_rssys1].set(0.0f);
    descriptor[iwdescr_rssys2].set(0.0f);
    descriptor[iwdescr_rssys3].set(0.0f);
    descriptor[iwdescr_rssys4].set(0.0f);
    descriptor[iwdescr_rssys5].set(0.0f);
    descriptor[iwdescr_rssys6].set(0.0f);
    descriptor[iwdescr_rssys7].set(0.0f);
    descriptor[iwdescr_rssys8].set(0.0f);
    descriptor[iwdescr_rssys9].set(0.0f);

    descriptor[iwdescr_ar5].set(0.0f);
    descriptor[iwdescr_ar6].set(0.0f);
    descriptor[iwdescr_al5].set(0.0f);
    descriptor[iwdescr_al6].set(0.0f);

    descriptor[iwdescr_rhaf].set(0.0f);
    descriptor[iwdescr_frafus].set(0.0f);
    descriptor[iwdescr_fraromha].set(0.0f);
    descriptor[iwdescr_srsz].set(0.0f);
    descriptor[iwdescr_nrsyscmr].set(0.0f);

    descriptor[iwdescr_spchtro].set(0.0f);
    descriptor[iwdescr_rbfrspch].set(0.0f);
    descriptor[iwdescr_satspcha].set(0.0f);
    descriptor[iwdescr_unsatspcha].set(0.0f);
    descriptor[iwdescr_fsatspcha].set(0.0f);
    descriptor[iwdescr_scaffoldbranches].set(0.0f);
    descriptor[iwdescr_nrnspch].set(0.0f);
    descriptor[iwdescr_fnrnspc].set(0.0f);

    descriptor[iwdescr_avcdbsub].set(0.0f);

    descriptor[iwdescr_npfsdsys].set(0.0f);
    descriptor[iwdescr_rnginpfs].set(0.0f);
    descriptor[iwdescr_lgplnfsy].set(0.0f);
    descriptor[iwdescr_htrcpfsy].set(0.0f);
    descriptor[iwdescr_mxhtpfsy].set(0.0f);

    return 1;
  }

  int nrings[kMaxRingSize];
  std::fill_n(nrings, kMaxRingSize, 0);

  std::unique_ptr<int[]> ring_already_done(new_int(nr));
  std::unique_ptr<int[]> atom_already_done(new_int(data.matoms));

  int rings_in_largest_system = 1;
  int atoms_in_largest_system = 0;
  int ring_systems_containing_multiple_rings = 0;

  int isolated_rings = 0;
  int isolated_heterocycles = 0;

  int max_heteroatoms_in_ring = 0;
  float max_ring_heteroatom_fraction = 0.0f;
  float min_ring_heteroatom_fraction = 1.0f;

  float ring_isolation_score = 0.0f;

  int aromatic_rings = 0;
  int aliphatic_rings = 0;

  int smallest_ring_size = 0;
  int largest_ring_size = 0;

  int exocyclic_bonds = 0;
  int double_bond_attachments = 0;
  int singly_connected_attachments = 0;
  int singly_connected_heteroatoms = 0;
  int singly_connected_donors = 0;

  int biphenyl = 0;  // arom_single_arom
  int arom_single_aliph = 0;
  int aliph_single_aliph = 0;

  int ring_chain_heteroatom_join_count = 0;
  int aromatic_ring_chain_join_count = 0;
  int aliphatic_ring_chain_join_count = 0;


  // Originally I had a function that traversed each ring looking for singly
  // connected attachments, but that is unstable when different SSSR ring
  // determinations are done, so just look for singly connected atoms that are
  // attached to a ring.
  // THis is probably more efficient than looking at all bonds attached to
  // all ring atoms.
  for (const Bond* b : m.bond_list()) {
    if (b->nrings()) {
      continue;
    }

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    if (ring_membership[a1] == 0 && ring_membership[a2] == 0) {
      continue;
    }

    // non ring bond bonded to at least one ring.
    ++exocyclic_bonds;

    if (b->is_double_bond()) {
      ++double_bond_attachments;
      continue;
    }

    if (is_aromatic_atom[a1] && is_aromatic_atom[a2]) {
      ++biphenyl;
      continue;
    } 
    
    if (is_aromatic_atom[a1] || is_aromatic_atom[a2]) {
      ++arom_single_aliph;
    } else {
      ++aliph_single_aliph;
    }

    if (ncon[a1] == 1 || ncon[a2] == 1) {
      ++singly_connected_attachments;
      if (z[a1] != 6 || z[a2] != 6) {  // ignore heteroatom-heteratom case.
        ++singly_connected_heteroatoms;
        if ((z[a1] != 6 && m.hcount(a1)) || 
            (z[a2] != 6 && m.hcount(a2))) {
          ++singly_connected_donors;
        }
      }
    } else if (z[a1] != 6 || z[a2] != 6) {
      ++ring_chain_heteroatom_join_count;
    }

    if (is_aromatic_atom[a1] || is_aromatic_atom[a2]) {
      ++ aromatic_ring_chain_join_count;
    } else if (! is_aromatic_atom[a1] && ! is_aromatic_atom[a2]) {
      ++aliphatic_ring_chain_join_count;
    }
  }

  const int ring_chain_join_count = arom_single_aliph + aliph_single_aliph;

  int more_than_7_atoms = 0;

  m.compute_aromaticity_if_needed();

  int ring_systems = nr;

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (ri->is_aromatic()) {
      ++aromatic_rings;
    } else {
      ++aliphatic_rings;
    }

    if (!ri->is_fused()) {
      ring_isolation_score += ComputeRingIsolation(m, *ri, ncon);
    }

    const int hac = HeteroatomsInRing(ri, z);
    if (hac > max_heteroatoms_in_ring) {
      max_heteroatoms_in_ring = hac;
    }

//  Unstable when non-sssr rings present - computed above.
//  exocyclic_bonds += ComputeExocyclicBonds(m, *ri, z, ncon, double_bond_attachments,
//                                           singly_connected_attachments,
//                                           singly_connected_heteroatoms,
//                                           singly_connected_donors);

    if (!ri->is_fused()) {
      ++isolated_rings;
      if (hac > 0) {
        ++isolated_heterocycles;
      }
    }

    const int rs = ri->number_elements();

    if (smallest_ring_size == 0) {
      smallest_ring_size = rs;
    }
    if (rs > largest_ring_size) {
      largest_ring_size = rs;
    }
    if (rs > 7) {
      ++more_than_7_atoms;
    }

    const float tmp = static_cast<float>(hac) / static_cast<float>(rs);
    if (tmp > max_ring_heteroatom_fraction) {
      max_ring_heteroatom_fraction = tmp;
    }
    if (tmp < min_ring_heteroatom_fraction) {
      min_ring_heteroatom_fraction = tmp;
    }

    if (rs < kMaxRingSize) {
      ++nrings[rs];
    }

    if (ring_already_done[i]) {
      continue;
    }

    int system_size = 1;
    int atoms_in_system = ri->number_elements();

    if (ri->is_fused()) {
      std::fill_n(atom_already_done.get(), data.matoms, 0);
      ri->set_vector(atom_already_done.get(), 1);

      for (int j = i + 1; j < nr; ++j) {
        const Ring* rj = m.ringi(j);
        if (rj->fused_system_identifier() == ri->fused_system_identifier()) {
          ++system_size;
          atoms_in_system += AtomsNotAlreadyMarked(*rj, atom_already_done.get());
          ring_already_done[j] = 1;
          --ring_systems;
        }
      }
    }

    if (system_size > rings_in_largest_system) {
      rings_in_largest_system = system_size;
    }
    if (system_size > 1) {
      ++ring_systems_containing_multiple_rings;
    }
    if (atoms_in_system > atoms_in_largest_system) {
      atoms_in_largest_system = atoms_in_system;
    }
  }

  for (int i = 0; i < m.non_sssr_rings(); ++i) {
    const Ring* ri = m.non_sssr_ring(i);

    const int hac = HeteroatomsInRing(ri, z);
    if (hac > max_heteroatoms_in_ring) {
      max_heteroatoms_in_ring = hac;
    }

    const int rs = ri->number_elements();
    const float tmp = static_cast<float>(hac) / static_cast<float>(rs);
    if (tmp > max_ring_heteroatom_fraction) {
      max_ring_heteroatom_fraction = tmp;
    }
  }

  descriptor[iwdescr_mhr].set(static_cast<float>(max_heteroatoms_in_ring));
  descriptor[iwdescr_mxhrf].set(max_ring_heteroatom_fraction);
  descriptor[iwdescr_mnhrf].set(min_ring_heteroatom_fraction);
  descriptor[iwdescr_lrsysz].set(static_cast<float>(rings_in_largest_system));
  descriptor[iwdescr_srsz].set(static_cast<float>(smallest_ring_size));
  descriptor[iwdescr_lrsz].set(static_cast<float>(largest_ring_size));
  descriptor[iwdescr_rng7atoms].set(static_cast<float>(more_than_7_atoms));

  descriptor[iwdescr_nrsyscmr].set(
      static_cast<float>(ring_systems_containing_multiple_rings));
  descriptor[iwdescr_mars].set(static_cast<float>(atoms_in_largest_system));
  descriptor[iwdescr_ringsys].set(static_cast<float>(ring_systems));
  descriptor[iwdescr_ringisol].set(ring_isolation_score);
  descriptor[iwdescr_isolrc].set(static_cast<float>(isolated_rings));
  descriptor[iwdescr_isolhtrc].set(static_cast<float>(isolated_heterocycles));

  descriptor[iwdescr_arring].set(static_cast<float>(aromatic_rings));
  descriptor[iwdescr_alring].set(static_cast<float>(aliphatic_rings));

  descriptor[iwdescr_excybond].set(static_cast<float>(exocyclic_bonds));
  descriptor[iwdescr_excydbond].set(static_cast<float>(double_bond_attachments));
  descriptor[iwdescr_excydscon].set(static_cast<float>(singly_connected_attachments));
  descriptor[iwdescr_excydsconh].set(static_cast<float>(singly_connected_heteroatoms));
  descriptor[iwdescr_excydscondon].set(static_cast<float>(singly_connected_donors));

  int offset = iwdescr_nrings3;
  for (int i = 3; i < kMaxRingSize; ++i, ++offset) {
    descriptor[offset].set(static_cast<float>(nrings[i]));
  }

  descriptor[iwdescr_rcj].set(static_cast<float>(ring_chain_join_count));
  descriptor[iwdescr_rchj].set(static_cast<float>(ring_chain_heteroatom_join_count));
  descriptor[iwdescr_amrcj].set(static_cast<float>(aromatic_ring_chain_join_count));
  descriptor[iwdescr_alrcj].set(static_cast<float>(aliphatic_ring_chain_join_count));
  return 1;
}

int
IWDescr::IWDescrImpl::IsSpiroFused(Molecule& m, atom_number_t a) const {
  // Migrated from legacy is_spiro_fused().
  // We look for atoms that are in two rings. If either ring is non-fused, or
  // the two fused rings are in different fused systems, this atom represents
  // a spiro fusion. The historical code deliberately ignored spiro-like atoms
  // wholly within one fused ring system.
  const int nr = m.nrings();

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (! ri->contains(a)) {
      continue;
    }

    if (! ri->is_fused()) {
      return 1;
    }

    for (int j = i + 1; j < nr; ++j) {
      const Ring* rj = m.ringi(j);

      if (! rj->contains(a)) {
        continue;
      }

      if (! rj->is_fused()) {
        return 1;
      }

      if (ri->fused_system_identifier() != rj->fused_system_identifier()) {
        return 1;
      }

      return 0;
    }
  }

  return 0;
}

int
IWDescr::IWDescrImpl::ComputeSpiroFusionDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy do_compute_spiro_fusions(m, ncon, ring_membership).
  // This is descriptor computation state because it writes iwdescr_nspiro and
  // depends on the compute_spiro_fusions option owned by IWDescrImpl.
  const int nr = m.nrings();

  if (! compute_spiro_fusions || nr < 2) {
    descriptor[iwdescr_nspiro].set(0.0f);
    return 1;
  }

  const int matoms = data.matoms;
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();

  int spiro_fusions = 0;

  for (int i = 0; i < matoms; ++i) {
    if (ncon[i] != 4) {
      continue;
    }

    if (ring_membership[i] < 2) {
      continue;
    }

    if (IsSpiroFused(m, i)) {
      ++spiro_fusions;
    }
  }

  descriptor[iwdescr_nspiro].set(static_cast<float>(spiro_fusions));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeRameyDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_ramey_descriptors(m, z).
  if (m.contains_non_periodic_table_elements()) {
    descriptor[iwdescr_obalance].set(0.0f);
    return 1;
  }

  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();

  int nc = 0;
  int nn = 0;
  int no = 0;
  int nf = 0;
  int ns = 0;
  int ncl = 0;
  int nbr = 0;
  int ni = 0;

  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t zi = z[i];

    if (zi == 6) {
      ++nc;
    } else if (zi == 7) {
      ++nn;
    } else if (zi == 8) {
      ++no;
    } else if (zi == 9) {
      ++nf;
    } else if (zi == 16) {
      ++ns;
    } else if (zi == 17) {
      ++ncl;
    } else if (zi == 35) {
      ++nbr;
    } else if (zi == 53) {
      ++ni;
    }
  }

  descriptor[iwdescr_rmync].set(static_cast<float>(nc));
  descriptor[iwdescr_rmynn].set(static_cast<float>(nn));
  descriptor[iwdescr_rmyno].set(static_cast<float>(no));
  descriptor[iwdescr_rmynf].set(static_cast<float>(nf));
  descriptor[iwdescr_rmyns].set(static_cast<float>(ns));
  descriptor[iwdescr_rmyncl].set(static_cast<float>(ncl));
  descriptor[iwdescr_rmynbr].set(static_cast<float>(nbr));
  descriptor[iwdescr_rmyni].set(static_cast<float>(ni));
  descriptor[iwdescr_rmy_heavy_halogen].set(static_cast<float>(ncl + nbr + ni));

  const int nh = m.implicit_hydrogens();

  const double oxygen_balance =
      -1600.0 * (2 * nc + static_cast<double>(nh) / 2.0 - no) /
      m.molecular_weight_ignore_isotopes();

  descriptor[iwdescr_obalance].set(static_cast<float>(oxygen_balance));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeLogPDescriptors(Molecule& m, PerMoleculeData& data) {
  (void)data;

  if (descriptors_to_compute.compute_xlogp) {
    std::optional<double> x = xlogp_engine.LogP(m);
    if (x) {
      descriptor[iwdescr_xlogp].set(*x);
    }
  }

  if (descriptors_to_compute.compute_alogp) {
    std::optional<double> x = alogp_engine.LogP(m);
    if (x) {
      descriptor[iwdescr_alogp].set(*x);
    }
  }

  return 1;
}



namespace {

// Ring-substitution helper functions migrated from the legacy file. They remain
// file-local because they do not own IWDescr state; the member function below
// owns descriptor storage and the computation-option gates.
void
AccumulateRingSubstitutionDistance(extending_resizable_array<int>& diffs, int n,
                                   int n2, int d) {
  assert(d > 0);

  if (d > n2) {
    d = n - d;
  }

  diffs[d]++;
}

int
AccumulateRingSubstitutionDistances(Molecule& m, const int* ncon,
                                    const Set_of_Atoms& par,
                                    extending_resizable_array<int>& diffs) {
  const int n = par.number_elements();

  resizable_array<int> offset_of_branch;

  for (int i = 0; i < n; ++i) {
    const atom_number_t j = par[i];

    if (ncon[j] == 2) {
      continue;
    }

    if (m.nrings(j) > 1) {
      continue;
    }

    offset_of_branch.add(i);
  }

  const int noffset = offset_of_branch.number_elements();
  if (noffset < 2) {
    return 1;
  }

  const int n2 = n / 2;

  if (noffset == 2) {
    AccumulateRingSubstitutionDistance(diffs, n, n2,
                                       offset_of_branch[1] - offset_of_branch[0]);
    return 1;
  }

  for (int i = 0; i < noffset; ++i) {
    const int di = offset_of_branch[i];

    for (int j = i + 1; j < noffset; ++j) {
      const int dj = offset_of_branch[j];
      AccumulateRingSubstitutionDistance(diffs, n, n2, dj - di);
    }
  }

  return 1;
}

int
BondedToAnotherRingForSubstitution(Molecule& m, const int* ring_membership,
                                   atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    if (b->nrings()) {
      continue;
    }

    const atom_number_t o = b->other(zatom);
    if (ring_membership[o]) {
      return 1;
    }
  }

  return 0;
}

int
OrthoRingCountForSubstitution(Molecule& m, const int* ncon,
                              const int* ring_membership, const Ring& ring) {
  int rc = 0;
  const int ring_size = ring.number_elements();

  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t a1 = ring[i];
    if (ncon[a1] == 2) {
      continue;
    }

    const atom_number_t a2 = ring[(i + 1) % ring_size];
    if (ncon[a2] == 2) {
      continue;
    }

    if (! BondedToAnotherRingForSubstitution(m, ring_membership, a1)) {
      continue;
    }

    if (BondedToAnotherRingForSubstitution(m, ring_membership, a2)) {
      ++rc;
    }
  }

  return rc;
}

int
CountRingBondsForAdjacentFusion(const Atom& a) {
  int rc = 0;

  for (const Bond* b : a) {
    if (b->nrings() > 0) {
      ++rc;
    }
  }

  return rc;
}

}  // namespace

int
IWDescr::IWDescrImpl::ComputeRingSubstitutionDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from the legacy ring-substitution, adjacent-fusion, and
  // ring-substitution-ratio helpers. The low-level accumulation helpers remain
  // file-local and stateless; this method owns the descriptor writes and option
  // gates.
  const int nr = m.nrings();
  const int matoms = data.matoms;
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const Atom** atom = data.atoms();

  if (descriptors_to_compute.ring_substitution_descriptors) {
    if (nr > 0) {
      extending_resizable_array<int> bonds_between_ring_substitutions_aromatic_ring;
      extending_resizable_array<int> bonds_between_ring_substitutions_aliphatic_ring;

      int fused_rings_present = 0;

      for (int i = 0; i < nr; ++i) {
        const Ring* ri = m.ringi(i);

        if (ri->is_fused()) {
          ++fused_rings_present;
          continue;
        }

        if (ri->is_aromatic()) {
          AccumulateRingSubstitutionDistances(
              m, ncon, *ri, bonds_between_ring_substitutions_aromatic_ring);
        } else {
          AccumulateRingSubstitutionDistances(
              m, ncon, *ri, bonds_between_ring_substitutions_aliphatic_ring);
        }
      }

      extending_resizable_array<int> bonds_between_ring_substitutions_ring_system;

      if (fused_rings_present) {
        std::unique_ptr<int[]> fsid(new_int(matoms, -1));
        std::unique_ptr<int[]> ring_already_done(new_int(nr));

        for (int i = 0; i < nr; ++i) {
          if (ring_already_done[i]) {
            continue;
          }

          const Ring* ri = m.ringi(i);
          const int f = ri->fused_system_identifier();

          if (f < 0) {
            continue;
          }

          int strongly_fused_ring_found = 0;
          if (ri->largest_number_of_bonds_shared_with_another_ring() > 1) {
            strongly_fused_ring_found = 1;
          }

          ri->set_vector(fsid.get(), f + 1);

          for (int j = i + 1; j < nr; ++j) {
            if (ring_already_done[j]) {
              continue;
            }

            const Ring* rj = m.ringi(j);
            if (f != rj->fused_system_identifier()) {
              continue;
            }

            ring_already_done[j] = 1;

            if (strongly_fused_ring_found) {
              continue;
            }

            if (rj->largest_number_of_bonds_shared_with_another_ring() > 1) {
              strongly_fused_ring_found = 1;
            } else {
              rj->set_vector(fsid.get(), f + 1);
            }
          }

          if (strongly_fused_ring_found) {
            continue;
          }

          Set_of_Atoms s;
          if (! path_around_edge_of_ring_system(m, fsid.get(), f + 1, s)) {
            continue;
          }

          AccumulateRingSubstitutionDistances(
              m, ncon, s, bonds_between_ring_substitutions_ring_system);
        }
      }

      descriptor[iwdescr_rsarom1].set(bonds_between_ring_substitutions_aromatic_ring[1]);
      descriptor[iwdescr_rsarom2].set(bonds_between_ring_substitutions_aromatic_ring[2]);
      descriptor[iwdescr_rsarom3].set(bonds_between_ring_substitutions_aromatic_ring[3]);

      descriptor[iwdescr_rsaliph1].set(bonds_between_ring_substitutions_aliphatic_ring[1]);
      descriptor[iwdescr_rsaliph2].set(bonds_between_ring_substitutions_aliphatic_ring[2]);
      descriptor[iwdescr_rsaliph3].set(bonds_between_ring_substitutions_aliphatic_ring[3]);

      if (bonds_between_ring_substitutions_aliphatic_ring.number_elements() < 4) {
        descriptor[iwdescr_rsaliph4].set(0);
      } else {
        int n4 = 0;
        for (int i = 4;
             i < bonds_between_ring_substitutions_aliphatic_ring.number_elements(); ++i) {
          n4 += bonds_between_ring_substitutions_aliphatic_ring[i];
        }
        descriptor[iwdescr_rsaliph4].set(n4);
      }

      descriptor[iwdescr_rssys1].set(bonds_between_ring_substitutions_ring_system[1]);
      descriptor[iwdescr_rssys2].set(bonds_between_ring_substitutions_ring_system[2]);
      descriptor[iwdescr_rssys3].set(bonds_between_ring_substitutions_ring_system[3]);
      descriptor[iwdescr_rssys4].set(bonds_between_ring_substitutions_ring_system[4]);
      descriptor[iwdescr_rssys5].set(bonds_between_ring_substitutions_ring_system[5]);
      descriptor[iwdescr_rssys6].set(bonds_between_ring_substitutions_ring_system[6]);
      descriptor[iwdescr_rssys7].set(bonds_between_ring_substitutions_ring_system[7]);
      descriptor[iwdescr_rssys8].set(bonds_between_ring_substitutions_ring_system[8]);
      descriptor[iwdescr_rssys9].set(bonds_between_ring_substitutions_ring_system[9]);
    }
  }

  if (descriptors_to_compute.adjacent_ring_fusion_descriptors) {
    int single_bonds_found = 0;
    int double_bonds_found = 0;

    (void)m.ring_membership();

    for (int i = 0; i < matoms; ++i) {
      if (ring_membership[i] < 2) {
        continue;
      }

      if (ncon[i] != 3) {
        continue;
      }

      const Atom* ai = m.atomi(i);

      // Very important, otherwise C1CC2C1C(C2CCCNC)C1=CC=CC=C1 fails.
      if (CountRingBondsForAdjacentFusion(*ai) != 3) {
        continue;
      }

      for (const Bond* b : *ai) {
        const atom_number_t k = b->other(i);
        const Atom* ak = m.atomi(k);

        if (CountRingBondsForAdjacentFusion(*ak) != 2) {
          continue;
        }

        // No exocyclic bonds here.
        if (ak->ncon() == 2) {
          continue;
        }

        for (const Bond* b2 : *ak) {
          if (b2->nrings()) {
            continue;
          }

          if (b2->is_single_bond()) {
            ++single_bonds_found;
          } else {
            ++double_bonds_found;
          }
        }
      }
    }

    descriptor[iwdescr_sboradjf].set(static_cast<float>(single_bonds_found));
    descriptor[iwdescr_dboradjf].set(static_cast<float>(double_bonds_found));
  }

  if (descriptors_to_compute.ring_substitution_ratio_descriptors) {
    int ring_atoms = 0;
    int substituted_outside_ring = 0;
    int simple_substituted_outside_ring = 0;

    if (nr == 0) {
      descriptor[iwdescr_frsub].set(0.0f);
      descriptor[iwdescr_frssub].set(0.0f);
      descriptor[iwdescr_arorthoring].set(0.0f);
      descriptor[iwdescr_alorthoring].set(0.0f);
      return 1;
    }

    for (int i = 0; i < matoms; ++i) {
      if (ring_membership[i] == 0) {
        continue;
      }

      ++ring_atoms;

      if (ncon[i] == 2) {
        continue;
      }

      const Atom* a = atom[i];

      for (int j = 0; j < ncon[i]; ++j) {
        const atom_number_t k = a->other(i, j);

        if (ring_membership[k]) {
          continue;
        }

        ++substituted_outside_ring;
        if (ncon[k] == 1) {
          ++simple_substituted_outside_ring;
        }
      }
    }

    descriptor[iwdescr_frsub].set(static_cast<float>(substituted_outside_ring) /
                                  static_cast<float>(ring_atoms));
    descriptor[iwdescr_frssub].set(static_cast<float>(simple_substituted_outside_ring) /
                                   static_cast<float>(ring_atoms));

    m.compute_aromaticity_if_needed();

    int aromatic_ortho_ring_count = 0;
    int aliphatic_ortho_ring_count = 0;

    for (const Ring* r : m.sssr_rings()) {
      if (r->is_aromatic()) {
        aromatic_ortho_ring_count +=
            OrthoRingCountForSubstitution(m, ncon, ring_membership, *r);
      } else {
        aliphatic_ortho_ring_count +=
            OrthoRingCountForSubstitution(m, ncon, ring_membership, *r);
      }
    }

    descriptor[iwdescr_arorthoring].set(aromatic_ortho_ring_count);
    descriptor[iwdescr_alorthoring].set(aliphatic_ortho_ring_count);
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeCrowdingDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_crowding_descriptors(m, ncon, atom).
  // Crowding is pure descriptor computation and consumes only per-molecule
  // connection/atom arrays plus IWDescr-owned descriptor storage.
  const int matoms = data.matoms;
  const int* ncon = data.connections();
  const Atom** atom = data.atoms();

  float rc = 0.0f;

  for (int i = 0; i < matoms; ++i) {
    if (ncon[i] < 3) {
      continue;
    }

    const Atom* ai = atom[i];
    for (int j = 0; j < ncon[i]; ++j) {
      const atom_number_t k = ai->other(i, j);

      if (ncon[k] > 2) {
        rc += 1.0f;
      } else if (ncon[k] == 2) {
        const Atom* ak = atom[k];

        for (int l = 0; l < 2; ++l) {
          const atom_number_t n = ak->other(k, l);

          if (n == i) {
            continue;
          }

          if (ncon[n] > 2) {
            rc += 0.5f;
          }

          break;
        }
      }
    }
  }

  rc = rc * 0.5f;  // The loop above double counts.

  descriptor[iwdescr_crowding].set(rc);
  descriptor[iwdescr_fcrowdng].set(rc / static_cast<float>(matoms));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeSpinachDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_spinach_descriptors(). Spinach is a
  // descriptor-family calculation, while the atom classification helpers are
  // file-local and stateless.
  const int matoms = data.matoms;
  if (matoms == 0) {
    return 1;
  }

  std::vector<int> spinach(matoms, 0);

  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const int* ring_membership = data.ring_membership_data();
  const Atom** atom = data.atoms();

  const int nedges = m.nedges();
  std::vector<int> bond_rotatable(nedges, 0);
  m.assign_bond_numbers_to_bonds_if_needed();

  int total_rotatable_bonds = 0;
  for (int bond_number = 0; bond_number < nedges; ++bond_number) {
    const Bond* b = m.bondi(bond_number);
    if (b->nrings() || ! b->is_single_bond()) {
      continue;
    }

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    if (ncon[a1] <= 1 || ncon[a2] <= 1 || TripleBondAtEitherEnd(m, b) ||
        part_of_otherwise_non_rotabable_entity(m, a1, a2)) {
      continue;
    }

    bond_rotatable[bond_number] = 1;
    ++total_rotatable_bonds;
  }

  descriptor[iwdescr_maxrotbgrp].set(
      static_cast<float>(LongestContiguousRotatableBondGroup(m, bond_rotatable.data())));

  if (m.nrings() == 0) {
    ++molecules_with_no_rings;
    descriptor[iwdescr_frspch].set(1.0f);
    descriptor[iwdescr_scafrotb].set(0.0f);
    descriptor[iwdescr_spchrotb].set(static_cast<float>(total_rotatable_bonds));
    return 1;
  }

  const int spinach_atoms = IdentifySpinachAtoms(m, spinach.data(), ring_membership, atom);

  descriptor[iwdescr_frspch].set(static_cast<float>(spinach_atoms) /
                                 static_cast<float>(matoms));

  int non_ring_non_spinach_atoms = 0;
  int heteroatoms_in_spinach = 0;
  int bonds_in_spinach = 0;
  int saturated_spinach_atoms = 0;
  int unsaturated_spinach_atoms = 0;

  for (int i = 0; i < matoms; ++i) {
    if (spinach[i] == 0) {
      if (ring_membership[i] == 0) {
        ++non_ring_non_spinach_atoms;
      }
    } else if (z[i] != 6) {
      ++heteroatoms_in_spinach;
    }

    const Atom* ai = atom[i];
    const int icon = ai->ncon();

    if (spinach[i]) {
      if (icon == ai->nbonds()) {
        ++saturated_spinach_atoms;
      } else {
        ++unsaturated_spinach_atoms;
      }

      for (int j = 0; j < icon; ++j) {
        const Bond* b = ai->item(j);
        atom_number_t k = b->other(i);

        if (k < i || spinach[k] == 0) {
          continue;
        }

        ++bonds_in_spinach;
      }
    }
  }

  int rotatable_bonds_in_scaffold = 0;
  for (int bond_number = 0; bond_number < nedges; ++bond_number) {
    if (! bond_rotatable[bond_number]) {
      continue;
    }
    const Bond* b = m.bondi(bond_number);
    if (spinach[b->a1()] == 0 && spinach[b->a2()] == 0) {
      ++rotatable_bonds_in_scaffold;
    }
  }
  const int rotatable_bonds_in_spinach =
      total_rotatable_bonds - rotatable_bonds_in_scaffold;

  descriptor[iwdescr_scafrotb].set(static_cast<float>(rotatable_bonds_in_scaffold));
  descriptor[iwdescr_spchrotb].set(static_cast<float>(rotatable_bonds_in_spinach));

  if (spinach_atoms > 0) {
    descriptor[iwdescr_spchtro].set(
        iwmisc::Fraction<float>(heteroatoms_in_spinach, spinach_atoms));
  } else {
    descriptor[iwdescr_spchtro].set(0.0f);
  }

  descriptor[iwdescr_nrnspch].set(static_cast<float>(non_ring_non_spinach_atoms));
  if (matoms == spinach_atoms) {
    descriptor[iwdescr_fnrnspc].set(1.0f);
  } else {
    descriptor[iwdescr_fnrnspc].set(
        static_cast<float>(non_ring_non_spinach_atoms) /
        static_cast<float>(matoms - spinach_atoms));
  }

  if (bonds_in_spinach) {
    descriptor[iwdescr_rbfrspch].set(
        static_cast<float>(rotatable_bonds_in_spinach) /
        static_cast<float>(bonds_in_spinach));
  } else {
    descriptor[iwdescr_rbfrspch].set(0.0f);
  }

  descriptor[iwdescr_satspcha].set(static_cast<float>(saturated_spinach_atoms));
  descriptor[iwdescr_unsatspcha].set(static_cast<float>(unsaturated_spinach_atoms));
  if (spinach_atoms > 0) {
    descriptor[iwdescr_fsatspcha].set(
        iwmisc::Fraction<float>(saturated_spinach_atoms, spinach_atoms));
  } else {
    descriptor[iwdescr_fsatspcha].set(0.0f);
  }

  int branches_in_scaffold = 0;
  for (int i = 0; i < matoms; ++i) {
    if (spinach[i] || ring_membership[i]) {
      continue;
    }

    branches_in_scaffold += WithinScaffoldBranchedForSpinach(m, i);
  }

  descriptor[iwdescr_scaffoldbranches].set(static_cast<float>(branches_in_scaffold));

  int terminal_rings = 0;
  int internal_rings = 0;
  int spinach_connections = 0;
  int non_spinach_connections = 0;

  const int nr = m.nrings();

  resizable_array<int> strongly_fused;
  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);
    if (ri->strongly_fused_ring_neighbours()) {
      strongly_fused.add_if_not_already_present(ri->fused_system_identifier());
    }
  }

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (strongly_fused.contains(ri->fused_system_identifier())) {
      continue;
    }

    int s = 0;
    int ns = 0;
    CountSpinachConnections(m, *ri, spinach.data(), s, ns);

    if (ns == 1) {
      ++terminal_rings;
    } else if (ns > 1) {
      ++internal_rings;
    }

    spinach_connections += s;
    non_spinach_connections += ns;
  }

  descriptor[iwdescr_trmnlrng].set(static_cast<float>(terminal_rings));
  descriptor[iwdescr_intrnlrng].set(static_cast<float>(internal_rings));
  descriptor[iwdescr_rng2spch].set(static_cast<float>(spinach_connections));
  descriptor[iwdescr_rng2bridge].set(static_cast<float>(non_spinach_connections));

  return 1;
}

#ifdef NOW_HANDLED_ELSEWHERE
int
IWDescr::IWDescrImpl::ComputeRingChainDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_ring_chain_descriptors(m, z, ncon, atom,
  // ring_membership).  These descriptors count ring atoms with chain
  // substituents, split by heteroatom substituents and aromatic/aliphatic
  // ring atoms. This is pure descriptor computation and belongs in IWDescrImpl.
  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();
  const Atom** atom = data.atoms();
  const int* ring_membership = data.ring_membership_data();
  const int* is_aromatic_atom = data.aromaticity();

  int ring_chain_join_count = 0;
  int ring_chain_heteroatom_join_count = 0;
  int aromatic_ring_chain_join_count = 0;
  int aliphatic_ring_chain_join_count = 0;

  for (int i = 0; i < matoms; ++i) {
    if (ring_membership[i] == 0) {
      continue;
    }

    // Two connected ring atoms cannot have an external chain attached.
    if (ncon[i] == 2) {
      continue;
    }

    const Atom* ai = atom[i];

    for (int j = 0; j < ncon[i]; ++j) {
      const atom_number_t k = ai->other(i, j);

      if (ring_membership[k]) {
        continue;
      }

      ++ring_chain_join_count;

      if (z[k] != 6) {
        ++ring_chain_heteroatom_join_count;
      }
      
      if (is_aromatic_atom[i]) {
        ++aromatic_ring_chain_join_count;
      } else {
        ++aliphatic_ring_chain_join_count;
      }
    }
  }

  descriptor[iwdescr_rcj].set(static_cast<float>(ring_chain_join_count));
  descriptor[iwdescr_rchj].set(static_cast<float>(ring_chain_heteroatom_join_count));
  descriptor[iwdescr_amrcj].set(static_cast<float>(aromatic_ring_chain_join_count));
  descriptor[iwdescr_alrcj].set(static_cast<float>(aliphatic_ring_chain_join_count));

  return 1;
}
#endif


// `r1` and `r2` are fused. Determine whether or not they are strongly fused.
int
IWDescr::IWDescrImpl::StronglyFused(const Ring& r1, const Ring& r2, int matoms,
                                    int* tmp) const {
  std::fill_n(tmp, matoms, 0);
  r1.set_vector(tmp, 1);
  r2.increment_vector(tmp);

  // If there are only two atoms in common, then these cannot be strongly fused.
  if (count_occurrences_of_item_in_array(2, matoms, tmp) == 2) {
    return 0;
  }

  // There could be a case where two rings had 3 atoms in common but not two
  // bonds. Imagine two rings fused somewhere and joined somewhere else in a
  // spiro fusion. Preserve the historical assumption that this is too rare to
  // worry about.

  if (max_difference_in_ring_size_for_strongly_fused == 0) {
    return 1;
  }

  // This works because number_elements() is a signed quantity.
  if (std::abs(r1.number_elements() - r2.number_elements()) >
      max_difference_in_ring_size_for_strongly_fused) {
    return 0;
  }

  return 1;
}

int
IWDescr::IWDescrImpl::GrowFusedSystem(Molecule& m, int ring_index,
                                      int* ring_already_done, int flag, int* atmp,
                                      int growing_strongly_fused_system) const {
  const Ring* ri = m.ringi(ring_index);

  const int max_bonds_shared = ri->largest_number_of_bonds_shared_with_another_ring();

  ring_already_done[ring_index] = flag;

  int rc = 1;

  const int matoms = m.natoms();
  const int number_neighbours = ri->fused_ring_neighbours();

  for (int j = 0; j < number_neighbours; ++j) {
    const Ring* rj = ri->fused_neighbour(j);

    const int k = rj->ring_number();

    // Decide whether we are growing a strong or weakly fused system.
    if (max_bonds_shared == 1) {  // growing a weakly fused system
      if (rj->largest_number_of_bonds_shared_with_another_ring() > 1) {
        continue;
      }
    } else if (rj->largest_number_of_bonds_shared_with_another_ring() < 2) {
      continue;  // pair of rings not strongly fused
    }

    if (ring_already_done[k]) {
      continue;
    }

    const int sf = StronglyFused(*ri, *rj, matoms, atmp);

    if (growing_strongly_fused_system && sf) {
    } else if (! growing_strongly_fused_system && ! sf) {
    } else {
      continue;
    }

    rc += GrowFusedSystem(m, k, ring_already_done, flag, atmp,
                          growing_strongly_fused_system);
  }

  return rc;
}

int
IWDescr::IWDescrImpl::JoinedRingsTooDifferentInSize(Molecule& m, const Ring& r,
                                                    int* tmp) const {
  const int matoms = m.natoms();

  for (int i = 0; i < r.fused_ring_neighbours(); ++i) {
    const Ring* rj = r.fused_neighbour(i);

    if (StronglyFused(r, *rj, matoms, tmp)) {
      return 0;
    }
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ComputePlanarFusedRings(Molecule& m, const atomic_number_t* z,
                                              int* ring_already_done, int* atmp) {
  const int matoms = m.natoms();
  const int nr = m.nrings();

  int number_planar_fused_systems_found = 0;
  int rings_in_planar_fused_systems = 0;
  int largest_planar_fused_system_size = 0;
  int max_heteroatoms_in_system = 0;
  int heterocycles_in_planar_fused_systems = 0;

  for (int i = 0; i < nr; ++i) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    if (! ri->is_fused()) {
      continue;
    }

    if (ri->largest_number_of_bonds_shared_with_another_ring() > 1) {
      continue;  // in a strongly fused system
    }

    const int flag =
        matoms + 2 * i;  // unique number that won't collide with a system id.

    const int planar_system_size = GrowFusedSystem(m, i, ring_already_done, flag,
                                                   atmp, 0);

    std::fill_n(atmp, matoms, 0);
    ri->set_vector(atmp, 1);

    ++number_planar_fused_systems_found;
    rings_in_planar_fused_systems += planar_system_size;

    if (planar_system_size > largest_planar_fused_system_size) {
      largest_planar_fused_system_size = planar_system_size;
    }

    if (HeteroatomsInRing(ri, z)) {
      ++heterocycles_in_planar_fused_systems;
    }

    for (int j = i + 1; j < nr; ++j) {
      if (flag != ring_already_done[j]) {
        continue;
      }

      const Ring* rj = m.ringi(j);
      rj->set_vector(atmp, 1);

      if (HeteroatomsInRing(rj, z)) {
        ++heterocycles_in_planar_fused_systems;
      }
    }

    int heteroatoms_in_system = 0;
    for (int j = 0; j < matoms; ++j) {
      if (atmp[j] == 0) {
        continue;
      }

      if (z[j] != 6) {
        ++heteroatoms_in_system;
      }
    }

    if (heteroatoms_in_system > max_heteroatoms_in_system) {
      max_heteroatoms_in_system = heteroatoms_in_system;
    }
  }

  descriptor[iwdescr_npfsdsys].set(static_cast<float>(number_planar_fused_systems_found));
  descriptor[iwdescr_rnginpfs].set(static_cast<float>(rings_in_planar_fused_systems));
  descriptor[iwdescr_lgplnfsy].set(static_cast<float>(largest_planar_fused_system_size));
  descriptor[iwdescr_htrcpfsy].set(static_cast<float>(heterocycles_in_planar_fused_systems));
  descriptor[iwdescr_mxhtpfsy].set(static_cast<float>(max_heteroatoms_in_system));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeNonPlanarFusedRings(Molecule& m,
                                                 const atomic_number_t* z,
                                                 int* ring_already_done,
                                                 int* atmp) {
  const int matoms = m.natoms();
  const int nr = m.nrings();

  int number_strongly_fused_systems_found = 0;
  int rings_in_strongly_fused_systems = 0;
  int largest_strongly_fused_system_size = 0;  // number of rings in it.
  int max_heteroatoms_in_system = 0;
  int heterocycles_in_strongly_fused_systems = 0;

  for (int i = 0; i < nr; ++i) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    if (! ri->is_fused()) {
      continue;
    }

    if (ri->largest_number_of_bonds_shared_with_another_ring() < 2) {
      continue;
    }

    if (max_difference_in_ring_size_for_strongly_fused == 0) {
    } else if (JoinedRingsTooDifferentInSize(m, *ri, atmp)) {
      continue;
    }

    const int flag = ri->fused_system_identifier() + 1;  // unique identifier.

    const int strongly_fused_size = GrowFusedSystem(m, i, ring_already_done,
                                                    flag, atmp, 1);

    ++number_strongly_fused_systems_found;
    rings_in_strongly_fused_systems += strongly_fused_size;

    if (strongly_fused_size > largest_strongly_fused_system_size) {
      largest_strongly_fused_system_size = strongly_fused_size;
    }

    std::fill_n(atmp, matoms, 0);

    for (int j = 0; j < nr; ++j) {
      if (flag != ring_already_done[j]) {
        continue;
      }

      const Ring* rj = m.ringi(j);
      rj->set_vector(atmp, 1);

      if (HeteroatomsInRing(rj, z)) {
        ++heterocycles_in_strongly_fused_systems;
      }
    }

    int heteroatoms_in_system = 0;
    for (int j = 0; j < matoms; ++j) {
      if (atmp[j] && z[j] != 6) {
        ++heteroatoms_in_system;
      }
    }

    if (heteroatoms_in_system > max_heteroatoms_in_system) {
      max_heteroatoms_in_system = heteroatoms_in_system;
    }
  }

  descriptor[iwdescr_nsfsdsys].set(static_cast<float>(number_strongly_fused_systems_found));
  descriptor[iwdescr_rnginsfs].set(static_cast<float>(rings_in_strongly_fused_systems));
  descriptor[iwdescr_lgstrfsy].set(static_cast<float>(largest_strongly_fused_system_size));
  descriptor[iwdescr_htrcsfsy].set(static_cast<float>(heterocycles_in_strongly_fused_systems));
  descriptor[iwdescr_mxhtsfsy].set(static_cast<float>(max_heteroatoms_in_system));

  return 1;
}


int
IWDescr::IWDescrImpl::ComputeComplexityDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_fused_rings(m, z). The fused-ring helpers
  // write Descriptor values and depend on max_difference_in_ring_size_for_strongly_fused,
  // so they belong inside IWDescrImpl rather than remaining file-scope helpers.
  if (! descriptors_to_compute.complexity_descriptors) {
    return 1;
  }

  const int nr = m.nrings();
  if (nr < 2) {
    ZeroRunFusionDescriptors();
    return 1;
  }

  std::unique_ptr<int[]> ring_already_done = std::make_unique<int[]>(nr);
  std::unique_ptr<int[]> atmp = std::make_unique<int[]>(data.matoms);

  ComputeNonPlanarFusedRings(m, data.atomic_numbers(), ring_already_done.get(),
                                   atmp.get());

  std::fill_n(ring_already_done.get(), nr, 0);

  return ComputePlanarFusedRings(m, data.atomic_numbers(), ring_already_done.get(),
                                 atmp.get());
}

int
IWDescr::IWDescrImpl::ComputeNovartisPsaDescriptor(Molecule& m, PerMoleculeData& data) {
  // Migrated from the legacy Novartis polar surface area assignment. This is
  // computation-only and intentionally runs before charge assignment, matching
  // the historical iwdescriptors() ordering.
  if (! descriptors_to_compute.psa) {
    return 1;
  }

  std::optional<double> psa = tpsa_calculator.PolarSurfaceArea(
      m, data.atomic_numbers(), data.atoms(), data.aromaticity());
  descriptor[iwdescr_nvrtspsa].set(static_cast<float>(psa.value_or(0.0)));

  return 1;
}


namespace {

// Molar refractivity atom contributions migrated from the legacy
// compute_molar_refractivity() helper. These remain file-local helpers because
// they do not own IWDescr state; the member function below owns descriptor
// storage and calls them in the historical order.
double
IdentifyMrHalogen(const Molecule& m, const atomic_number_t* z, int* already_done) {
  const int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; ++i) {
    if (z[i] == 9) {
      rc += 1.0632;
      already_done[i] = 15;
    } else if (z[i] == 17) {
      rc += 5.6105;
      already_done[i] = 16;
    } else if (z[i] == 35) {
      rc += 8.6782;
      already_done[i] = 17;
    } else if (z[i] == 53) {
      rc += 13.8741;
      already_done[i] = 18;
    }
  }

  return rc;
}

double
IdentifyMrOxygen(const Molecule& m, const atomic_number_t* z, const Atom* const* atom,
                 int* already_done) {
  const int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 8 || already_done[i]) {
      continue;
    }

    const Atom* ai = atom[i];

    if (ai->ncon() == 2) {
      rc += 1.6351;
      already_done[i] = 8;
      continue;
    }

    if (ai->ncon() == 0) {
      continue;
    }

    const Bond* b = ai->item(0);

    if (b->is_single_bond()) {
      rc += 1.6351;
      already_done[i] = 7;
      continue;
    }

    const atom_number_t n = b->other(i);

    if (z[n] == 7) {
      rc += 2.1407;
      already_done[i] = 9;
    } else {
      rc += 1.7956;
      already_done[i] = 8;
    }
  }

  return rc;
}

double
IdentifyMrSulphur(const Molecule& m, const atomic_number_t* z, const Atom* const* atom,
                  int* already_done) {
  double rc = 0.0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 16 || already_done[i]) {
      continue;
    }

    const Atom* ai = atom[i];
    const int icon = ai->ncon();

    if (icon == 1) {
      if (ai->nbonds() == 1) {
        rc += 7.3190;
        already_done[i] = 19;
      } else {
        rc += 9.1680;
        already_done[i] = 20;
      }
      continue;
    }

    if (icon == ai->nbonds()) {
      rc += 7.3190;
      already_done[i] = 19;
      continue;
    }

    atom_number_t o1 = INVALID_ATOM_NUMBER;
    atom_number_t o2 = INVALID_ATOM_NUMBER;

    for (int j = 0; j < icon; ++j) {
      const Bond* b = ai->item(j);
      if (! b->is_double_bond()) {
        continue;
      }

      const atom_number_t k = b->other(i);
      if (z[k] != 8 || already_done[k]) {
        continue;
      }

      if (o1 == INVALID_ATOM_NUMBER) {
        o1 = k;
      } else {
        o2 = k;
      }
    }

    if (o2 != INVALID_ATOM_NUMBER) {
      rc += 5.3321;
      already_done[i] = 22;
    } else if (o1 != INVALID_ATOM_NUMBER) {
      rc += 6.0762;
      already_done[i] = 21;
    } else {
      rc += 9.1680;
      already_done[i] = 20;
    }
  }

  return rc;
}

double
IdentifyMrPhosphorus(const Molecule& m, const atomic_number_t* z,
                     const Atom* const* atom, int* already_done) {
  (void) atom;
  (void) already_done;

  double rc = 0.0;
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (z[i] == 15) {
      ++rc;
    }
  }

  return 5.3 * rc;
}

double
IdentifyMrNitro(const Molecule& m, const atomic_number_t* z, const Atom* const* atom,
                int* already_done) {
  const int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 8 || already_done[i]) {
      continue;
    }

    const Atom* ai = atom[i];
    if (ai->ncon() != 1) {
      continue;
    }

    const Bond* b = ai->item(0);
    if (! b->is_double_bond()) {
      continue;
    }

    const atom_number_t n = b->other(i);
    if (z[n] != 7 || already_done[n]) {
      continue;
    }

    const Atom* an = atom[n];
    if (an->ncon() != 3 || an->nbonds() != 5) {
      continue;
    }

    atom_number_t other_oxygen = INVALID_ATOM_NUMBER;
    for (int j = 0; j < 3; ++j) {
      const Bond* bn = an->item(j);
      if (! bn->is_double_bond()) {
        continue;
      }

      const atom_number_t k = bn->other(n);
      if (z[k] == 8 && ! already_done[k]) {
        other_oxygen = k;
        break;
      }
    }

    if (other_oxygen == INVALID_ATOM_NUMBER) {
      continue;
    }

    already_done[n] = 13;
    rc += 3.5054;
  }

  return rc;
}

double
IdentifyMrNitrogen(Molecule& m, const atomic_number_t* z, const Atom* const* atom,
                   int* already_done) {
  const int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 7 || already_done[i]) {
      continue;
    }

    if (m.is_aromatic(i)) {
      rc += 2.7662;
      already_done[i] = 12;
      continue;
    }

    const Atom* n = atom[i];
    const int ncon = n->ncon();

    if (n->nbonds() == ncon) {
      rc += 3.0100;
      already_done[i] = 10;
      continue;
    }

    if (ncon == 1) {
      rc += 3.2009;
      already_done[i] = 11;
      continue;
    }

    atom_number_t aromatic_neighbour = INVALID_ATOM_NUMBER;
    atom_number_t doubly_bonded_heteroatom = INVALID_ATOM_NUMBER;
    atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;
    atom_number_t singly_bonded_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon; ++j) {
      const Bond* b = n->item(j);
      const atom_number_t k = b->other(i);

      if (b->is_single_bond()) {
        if (m.is_aromatic(k)) {
          aromatic_neighbour = k;
        } else if (z[k] == 8) {
          singly_bonded_oxygen = k;
        } else if (z[k] == 7) {
          singly_bonded_nitrogen = k;
        }
        continue;
      }

      if (b->is_double_bond() && z[k] != 6) {
        doubly_bonded_heteroatom = k;
      }
    }

    if (doubly_bonded_heteroatom != INVALID_ATOM_NUMBER &&
        aromatic_neighbour != INVALID_ATOM_NUMBER) {
      already_done[i] = 14;
      rc += 3.8095;
    } else if (singly_bonded_oxygen != INVALID_ATOM_NUMBER) {
      already_done[i] = 23;
      rc += 4.2109;
    } else if (singly_bonded_nitrogen != INVALID_ATOM_NUMBER) {
      already_done[i] = 23;
      rc += 3.9809;
    } else {
      already_done[i] = 11;
      rc += 3.2009;
    }
  }

  return rc;
}

double
IdentifyMrCarbon(Molecule& m, const atomic_number_t* z, const Atom* const* atom,
                 int* already_done) {
  const int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 6 || already_done[i]) {
      continue;
    }

    if (m.is_aromatic(i)) {
      already_done[i] = 4;
      rc += 3.5090;
      continue;
    }

    const Atom* c = atom[i];
    const int ncon = c->ncon();

    if (ncon == c->nbonds()) {
      already_done[i] = 1;
      rc += 2.8158;
      continue;
    }

    atom_number_t doubly_bonded_heteroatom = INVALID_ATOM_NUMBER;
    int found_triple_bond = 0;

    for (int j = 0; j < ncon; ++j) {
      const Bond* b = c->item(j);
      if (b->is_single_bond()) {
        continue;
      }

      if (b->is_triple_bond()) {
        found_triple_bond = 1;
        break;
      }

      const atom_number_t k = b->other(i);
      if (z[k] != 6) {
        doubly_bonded_heteroatom = k;
      }
    }

    if (doubly_bonded_heteroatom != INVALID_ATOM_NUMBER) {
      rc += 3.0887;
      already_done[i] = 5;
    } else if (found_triple_bond) {
      already_done[i] = 3;
      rc += 3.8974;
    } else {
      rc += 3.8278;
      already_done[i] = 2;
    }
  }

  return rc;
}

}  // namespace

int
IWDescr::IWDescrImpl::ComputeMolarRefractivityDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_molar_refractivity(). The atom contribution
  // helpers remain file-local because they do not own descriptor state; this
  // method owns the descriptor write and uses PerMoleculeData for atom arrays.
  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const Atom* const* atom = data.atoms();

  std::unique_ptr<int[]> already_done(new_int(matoms));

  double rc = 0.0;
  rc += IdentifyMrHalogen(m, z, already_done.get());
  rc += IdentifyMrSulphur(m, z, atom, already_done.get());
  rc += IdentifyMrPhosphorus(m, z, atom, already_done.get());
  rc += IdentifyMrNitro(m, z, atom, already_done.get());
  rc += IdentifyMrOxygen(m, z, atom, already_done.get());
  rc += IdentifyMrNitrogen(m, z, atom, already_done.get());
  rc += IdentifyMrCarbon(m, z, atom, already_done.get());

  for (int i = 0; i < matoms; ++i) {
    rc += m.implicit_hydrogens(i) * 0.9155;
  }

  descriptor[iwdescr_cmr].set(static_cast<float>(rc / 10.0));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeRuleOfFiveDescriptors(Molecule& m, PerMoleculeData& data) {
  // The counts themselves live in Molecule_Lib. This used to be a local loop,
  // and it disagreed with the copies in molecule_filter and the python
  // bindings over whether a donor was a hydrogen or a heteroatom carrying one.
  // ro5_ohnh is "the sum of OHs and NHs", so it counts hydrogens.
  int on = 0;
  int ohnh = 0;
  m.LipinskiHbaHbd(on, ohnh);

  descriptor[iwdescr_ro5_ohnh].set(static_cast<float>(ohnh));
  descriptor[iwdescr_ro5_on].set(static_cast<float>(on));

  return 1;
}

void
AppendShell(const Molecule& m, atom_number_t zatom, int radius, const int* dm,
            Set_of_Atoms& shell) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == zatom) {
      continue;
    }
    if (dm[zatom * matoms + i] == radius) {
      shell << i;
    }
  }
}

int
IWDescr::IWDescrImpl::do_compute_mean_shell_occupancies(const Molecule& m, const int* dm) {
  const int matoms = m.natoms();
  const int max_radius = 3;

  std::unique_ptr<Accumulator_Int<int>[]> acc(new Accumulator_Int<int>[max_radius + 1]);

  Set_of_Atoms shell;  // Scope here for efficiency.
  for (int i = 0; i < matoms; ++i) {
    shell.resize_keep_storage(0);
    for (int r = 1; r <= max_radius; ++r) {
      AppendShell(m, i, r, dm, shell);
      acc[r].extra(shell.number_elements());
    }
  }

  descriptor[iwdescr_aveshell1].set(acc[1].average());
  descriptor[iwdescr_aveshell2].set(acc[2].average());
  descriptor[iwdescr_aveshell3].set(acc[3].average());
  descriptor[iwdescr_maxshell3].set(acc[3].maxval());

  return 1;
}

#define CURRENT_SHELL 30
#define NEXT_SHELL 99
#define SHELL_COMPLETED -5

/*
  We determine the atoms that have the smallest total distance to
  other atoms in the molecule.
  We count them.

  Then we do a couple of shell expansions to count the number of atoms
  encountered.
*/

int
IWDescr::IWDescrImpl::do_compute_descriptors_related_to_centroid_atom(Molecule& m, const int* dm,
                                                const int longest_path,
                                                const atomic_number_t* z, const int* ncon,
                                                int* in_shell) {
  const auto matoms = m.natoms();

  Set_of_Atoms atoms_at_min;
  int mind = longest_path * longest_path;

  for (auto i = 0; i < matoms; ++i) {
    int tot_dist = std::accumulate(dm + i * matoms, dm + i * matoms + matoms, 0);

    if (tot_dist < mind) {
      atoms_at_min.resize_keep_storage(0);
      atoms_at_min.add(i);
      mind = tot_dist;
    } else if (tot_dist == mind) {
      atoms_at_min.add(i);
    }
  }

  std::fill_n(in_shell, matoms, 0);

  int atoms_in_next_shell = 0;

  for (auto i = 0; i < atoms_at_min.number_elements(); ++i) {
    const auto j = atoms_at_min[i];
    // m.set_isotope(j, 1);
    in_shell[j] = CURRENT_SHELL;

    const auto a = m.atomi(j);
    for (auto k = 0; k < ncon[j]; ++k) {
      const atom_number_t s = a->other(j, k);
      if (atoms_at_min.contains(s)) {
        continue;
      }

      atoms_in_next_shell++;
      in_shell[s] = NEXT_SHELL;
    }
  }

  descriptor[iwdescr_cntrdgncy].set(atoms_at_min.number_elements());
  descriptor[iwdescr_cntrdshell1].set(atoms_in_next_shell);

  const int max_centroid_bond_radius = 2;  // change if ever needed

  for (auto r = 1; r <= max_centroid_bond_radius; ++r) {
    for (auto i = 0; i < matoms; ++i) {
      if (0 == in_shell[i]) {
        continue;
      }

      if (CURRENT_SHELL == in_shell[i]) {
        in_shell[i] = SHELL_COMPLETED;
      } else if (NEXT_SHELL == in_shell[i]) {
        in_shell[i] = CURRENT_SHELL;
      }
    }

    atoms_in_next_shell = 0;

    for (auto i = 0; i < matoms; ++i) {
      if (CURRENT_SHELL != in_shell[i]) {
        continue;
      }

      const auto a = m.atomi(i);

      // cerr << "Processing atom " << i << " ncon " << a->ncon() << " array " << ncon[i]
      // << '\n';

      for (auto j = 0; j < ncon[i]; ++j) {
        //      cerr << "from atom " << i << " ncon " << ncon[i] << " connection " << j <<
        //      " type " << m.smarts_equivalent_for_atom(i) << '\n';
        atom_number_t k = a->other(i, j);

        if (0 != in_shell[k]) {
          continue;
        }

        atoms_in_next_shell++;
        in_shell[k] = NEXT_SHELL;
      }
    }
    descriptor[iwdescr_cntrdshell1 + r].set(atoms_in_next_shell);
  }

  // cerr << m.smiles() << ' ' << m.name() << '\n';

  return 1;
}


static void
compute_shortest_distance_from_longest_path(Molecule& m, const int* dm,
                                            const int* in_path,
                                            int* shortest_distance_from_longest_path) {
  int matoms = m.natoms();
#ifdef DEBUG_DISTANCE_MATRIX
  for (int i = 0; i < matoms; ++i) {
    cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " in_path "
         << in_path[i] << '\n';
  }
#endif

  for (int i = 0; i < matoms; i++) {  // atom I is outside the longest path
    if (in_path[i]) {
      continue;
    }

    for (int j = 0; j < matoms; j++) {  // atom J is in the longest path
      if (!in_path[j]) {
        continue;
      }

      int b = dm[i * matoms + j];

      if (b < shortest_distance_from_longest_path[i]) {
        shortest_distance_from_longest_path[i] = b;
      }
    }
  }

  return;
}

// There are molecules for which enumerating all paths between two most separated
// atoms blows up. Imagine a long series of para substituted rings. At each ring
// there is a choice of which way to go, so it explodes as 2^N. If we detect that
// problem, we abort the calculation.
static constexpr uint32_t kMaxVisited = 16384;  // an arbitrary choice, 2^14

static void
do_compute_distances_from_longest_path_descriptors(
    Molecule& m, const int* dm, atom_number_t astart, atom_number_t astop, int* in_path,
    int* shortest_distance_from_longest_path,
    uint32_t& visited) {
  if (visited > kMaxVisited) {
    return;
  }

  in_path[astart] = 1;

  const int matoms = m.natoms();

  const Atom* a = m.atomi(astart);

  // const int current_distance = m.bonds_between(astart, astop);
  const int current_distance = dm[astart * matoms + astop];

  for (const Bond* b : *a) {
    const atom_number_t j = b->other(astart);

    if (j == astop) {
      compute_shortest_distance_from_longest_path(m, dm, in_path,
                                                  shortest_distance_from_longest_path);
      ++visited;
      if (visited > kMaxVisited) {
        return;
      }
    } else if (in_path[j]) {  // no doubling back
      ;
    } else if (current_distance - 1 == dm[j * matoms + astop]) {
      do_compute_distances_from_longest_path_descriptors(
          m, dm, j, astop, in_path, shortest_distance_from_longest_path, visited);
    }
  }

  in_path[astart] = 0;

  return;
}

// Given separated atoms `a1` and `a2` compute the mean distance of
// all other atoms to these.
int
IWDescr::IWDescrImpl::ComputeAtomDistributionAlongLongestPath(Molecule& m, const int* dm,
                                             atom_number_t a1, atom_number_t a2) {
  Accumulator_Int<int> to_a1, to_a2;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == a1 || i == a2) {
      continue;
    }
    int d = dm[a1 * matoms + i];
    to_a1.extra(d);
    d = dm[a2 * matoms + i];
    to_a2.extra(d);
  }

  // cerr << "Between atoms " << a1 << ' ' << to_a1 << '\n';
  // cerr << "Between atoms " << a2 << ' ' << to_a2 << '\n';

  if (to_a1.n() == 0) {
    descriptor[iwdescr_mdallp] = 0.0f;
    descriptor[iwdescr_fmdallp] = 0.0f;
    descriptor[iwdescr_fdiffallp] = 0.0f;
    return 1;
  }

  const float mean1 = to_a1.average();
  const float mean2 = to_a2.average();
  float minval = std::min(mean1, mean2);
  descriptor[iwdescr_mdallp] = minval;
  descriptor[iwdescr_fmdallp] = minval / static_cast<float>(dm[a1 * matoms + a2]);
  const float diff = abs(mean1 - mean2);
  descriptor[iwdescr_fdiffallp] = diff / static_cast<float>(dm[a1 * matoms + a2]);

  return 1;
}

// `longest_path_atom` will be set for all atoms on any of the longest
// paths. For each atom that is NOT on longest path atom, compute
// the shortest distance to an atom that is on a longest path.
// `shortest_distance_from_longest_path[i]` comes in initialised
// to `matoms`.
void
MaxDistLongestPaths(const Molecule& m,
                    const int* dm,
                    const int* longest_path_atom,
                    int* shortest_distance_from_longest_path) {
  const int matoms = m.natoms();
  // Loop over all atoms not on a longest path.
  for (int i = 0; i < matoms; ++i) {
    if (longest_path_atom[i]) {
      continue;
    }

    // Examine distances from `i` to all atoms on the longest path.
    for (int j = 0; j < matoms; ++j) {
      if (! longest_path_atom[j]) {
        continue;
      }

      const int d = dm[i * matoms + j];
      if (d < shortest_distance_from_longest_path[i]) {
        shortest_distance_from_longest_path[i] = d;
      }
    }
  }
}

int
IWDescr::IWDescrImpl::ComputeDistancesFromLongestPathInstance(Molecule& m, const int* dm,
                                                   atom_number_t a1, atom_number_t a2,
                                                   int* in_path, const int longest_path) {
  ComputeAtomDistributionAlongLongestPath(m, dm, a1, a2);

  int matoms = m.natoms();

  std::fill_n(in_path, matoms, 0);
  in_path[a1] = 1;
  in_path[a2] = 1;

  int* shortest_distance_from_longest_path = new_int(matoms, matoms);
  std::unique_ptr<int[]> free_sdlp(shortest_distance_from_longest_path);

  uint32_t visited = 0;
  do_compute_distances_from_longest_path_descriptors(m, dm, a1, a2, in_path,
                shortest_distance_from_longest_path, visited);

  // If we have detected a combinatorial blow-up switch to a different computation.
  // Note that it does not compute the same values, but a reasonable approximation.
  if (visited >= kMaxVisited) {
    Set_of_Atoms longest_path_atoms;
    m.AllAtomsBetween(a1, a2, longest_path_atoms);
    longest_path_atoms.set_vector(in_path, 1);
    MaxDistLongestPaths(m, dm, in_path, shortest_distance_from_longest_path);
  }

  Accumulator_Int<int> acc;
  for (int i = 0; i < matoms; i++) {
    if (matoms != shortest_distance_from_longest_path[i]) {
      acc.extra(shortest_distance_from_longest_path[i]);
    }
  }

  if (0 == acc.n()) {  // linear chain type molecule
    descriptor[iwdescr_mxsdlp].set(0.0);
    descriptor[iwdescr_avsdlp].set(0.0);
    descriptor[iwdescr_mxsdlprl].set(0.0);
  } else {
    descriptor[iwdescr_mxsdlp].set(acc.maxval());
    descriptor[iwdescr_avsdlp].set(acc.average());
    descriptor[iwdescr_mxsdlprl].set(static_cast<float>(acc.maxval()) /
                                     static_cast<float>(longest_path));
  }

  return 1;
}

/*
  This is complicated by the fact that there can be many paths of
  maximum length. To break ties, we use the molecular canonicalisation
  function
*/

int
IWDescr::IWDescrImpl::ComputeDistancesFromLongestPath(Molecule& m, const int* dm,
                                                   const int longest_path, int* in_path) {
  assert(nullptr != dm);

  const int matoms = m.natoms();

  Set_of_Atoms a1, a2;

  for (int i = 0; i < matoms; i++) {
    for (int j = i + 1; j < matoms; j++) {
      if (longest_path != dm[i * matoms + j]) {
        continue;
      }

      a1.add(i);
      a2.add(j);
    }
  }

  // cerr << "longest_path " <<longest_path << " atoms " << a1 << " and " << a2 << '\n';

  // Single atom molecule, or all atoms disconnected.
  if (a1.empty()) {
    return 1;
  }

  int n = a1.number_elements();

  if (1 == n) {
    return ComputeDistancesFromLongestPathInstance(m, dm, a1[0], a2[0],
                                                   in_path, longest_path);
  }

  int highest_score = -1;
  int pair_with_highest_score = -1;

  for (int i = 0; i < n; i++) {
    int c1 = m.canonical_rank(a1[i]);
    int c2 = m.canonical_rank(a2[i]);

    int canonical_sum;

    if (c1 > c2) {
      canonical_sum = c1 * (matoms + matoms) + c2;
    } else {
      canonical_sum = c2 * (matoms + matoms) + c1;
    }

    if (canonical_sum > highest_score) {
      highest_score = canonical_sum;
      pair_with_highest_score = i;
    }
  }

  return ComputeDistancesFromLongestPathInstance(
      m, dm, a1[pair_with_highest_score], a2[pair_with_highest_score], in_path,
      longest_path);
}

int
IWDescr::IWDescrImpl::ComputeDistanceMatrixDescriptors(Molecule& m, PerMoleculeData& data) {
  const int matoms = data.matoms;
  if (matoms <= 1) {
    return 1;
  }

  std::unique_ptr<int[]> working_storage = std::make_unique<int[]>(matoms + matoms);

  const int* dm = m.distance_matrix_warning_may_change();

  const int longest_path = ComputeDistanceMatrixDescriptors(m, data, dm, 
                working_storage.get());

  ComputeDistancesFromLongestPath(m, dm, longest_path, working_storage.get());

  do_compute_descriptors_related_to_centroid_atom(
      m, dm, longest_path, data.atomic_numbers(), data.connections(),
      working_storage.get());

  do_compute_mean_shell_occupancies(m, dm);

  compute_double_bond_substitution(m, data.atomic_numbers(), data.connections());

  return longest_path;
}

int
IWDescr::IWDescrImpl::ComputeDistanceMatrixDescriptors(Molecule& m, PerMoleculeData& data,
                const int* dm,
                int* eccentricity) {
  const int matoms = m.natoms();

  const atomic_number_t* z = data.atomic_numbers();
  const int* ncon = data.connections();

  int* totd = eccentricity + matoms;

  int terminal_methyl_groups = 0;
  int terminal_groups_separated_by_three_bonds = 0;
  int max_distance_between_ring_atoms = 0;
  int max_distance_between_aromatic_atoms = 0;
  int max_distance_between_heteratoms = 0;
  int max_distance_between_ons = 0;
  double harary = 0.0;

  Accumulator_Int<int> bonds_between_stats;

  std::unique_ptr<uint32_t[]> atype = AssignDistanceMatrixAtomTypes(m, z);

  for (int i = 0; i < matoms; ++i) {
    const int terminal_group = (ncon[i] == 1);

    const int i_is_ring = atype[i] & kDistanceMatrixIsRing;
    const int i_is_aromatic = atype[i] & kDistanceMatrixIsAromatic;
    const int i_is_heteroatom = atype[i] & kDistanceMatrixIsHeteroatom;
    const int i_is_ons = atype[i] & kDistanceMatrixIsONS;

    if (terminal_group && z[i] == 6 && m.nbonds(i) == 1) {
      ++terminal_methyl_groups;
    }

    for (int j = i + 1; j < matoms; ++j) {
      const int d = dm[i * matoms + j];

      bonds_between_stats.extra(d);
      harary += 1.0 / static_cast<float>(d);

      totd[i] += d;
      totd[j] += d;

      if (terminal_group && ncon[j] == 1 && d == 3) {
        ++terminal_groups_separated_by_three_bonds;
      }

      if (d > eccentricity[i]) {
        eccentricity[i] = d;
      }
      if (d > eccentricity[j]) {
        eccentricity[j] = d;
      }

      if (i_is_ring && d > max_distance_between_ring_atoms &&
          (atype[j] & kDistanceMatrixIsRing)) {
        max_distance_between_ring_atoms = d;
      }
      if (i_is_aromatic && d > max_distance_between_aromatic_atoms &&
          (atype[j] & kDistanceMatrixIsAromatic)) {
        max_distance_between_aromatic_atoms = d;
      }
      if (i_is_heteroatom && d > max_distance_between_heteratoms &&
          (atype[j] & kDistanceMatrixIsHeteroatom)) {
        max_distance_between_heteratoms = d;
      }
      if (i_is_ons && d > max_distance_between_ons &&
          (atype[j] & kDistanceMatrixIsONS)) {
        max_distance_between_ons = d;
      }
    }
  }

  descriptor[iwdescr_harary].set(harary);
  descriptor[iwdescr_weiner].set(static_cast<float>(bonds_between_stats.sum()));

  int max_eccentricity = 0;  // same as longest path.
  int muldiam = 0;
  int min_eccentricity = matoms;
  int mulrad = 0;

  int min_tot_d = totd[0];
  Set_of_Atoms atoms_at_min_tot_d;

  int max_tot_d = totd[0];

  for (int i = 0; i < matoms; ++i) {
    const int ecc = eccentricity[i];

    if (ecc > max_eccentricity) {
      max_eccentricity = ecc;
      muldiam = 1;
    } else if (ecc == max_eccentricity) {
      ++muldiam;
    }

    if (ecc < min_eccentricity) {
      min_eccentricity = ecc;
      mulrad = 1;
    } else if (ecc == min_eccentricity) {
      ++mulrad;
    }

    if (totd[i] > max_tot_d) {
      max_tot_d = totd[i];
    } else if (totd[i] < min_tot_d) {
      min_tot_d = totd[i];
      atoms_at_min_tot_d.resize_keep_storage(0);
      atoms_at_min_tot_d.add(i);
    } else if (totd[i] == min_tot_d) {
      atoms_at_min_tot_d.add(i);
    }
  }

  descriptor[iwdescr_mxdst].set(static_cast<float>(max_eccentricity));
  descriptor[iwdescr_muldiam].set(static_cast<float>(muldiam));
  descriptor[iwdescr_rad].set(static_cast<float>(min_eccentricity));
  descriptor[iwdescr_mulrad].set(static_cast<float>(mulrad));

  descriptor[iwdescr_tm].set(static_cast<float>(terminal_methyl_groups));
  descriptor[iwdescr_tg3].set(static_cast<float>(terminal_groups_separated_by_three_bonds));

  descriptor[iwdescr_mxdst].set(static_cast<float>(max_eccentricity));
  descriptor[iwdescr_fmxdst].set(iwmisc::Fraction<float>(max_eccentricity, matoms));

  if (max_eccentricity > 0) {
    descriptor[iwdescr_ishape].set(static_cast<float>(max_eccentricity - min_eccentricity) /
                                   static_cast<float>(max_eccentricity));
  }

  descriptor[iwdescr_maxdrng].set(max_distance_between_ring_atoms);
  descriptor[iwdescr_maxdarom].set(max_distance_between_aromatic_atoms);
  descriptor[iwdescr_maxdhtro].set(max_distance_between_heteratoms);
  descriptor[iwdescr_maxdons].set(max_distance_between_ons);
  descriptor[iwdescr_avebbtwn].set(
      bonds_between_stats.average_if_available_minval_if_not());
  descriptor[iwdescr_normbbtwn].set(
      bonds_between_stats.average_if_available_minval_if_not() /
      static_cast<float>(matoms));

  descriptor[iwdescr_compact].set(1.0f - static_cast<float>(max_eccentricity) /
                                             static_cast<float>(matoms));
  descriptor[iwdescr_nolp].set(matoms - max_eccentricity - 1);

  descriptor[iwdescr_avdcentre].set(static_cast<float>(min_tot_d) /
                                    static_cast<float>(matoms - 1));

  Accumulator_Int<int> from_middle;

  int atoms_within_3_bonds = 0;
  int heteroatoms_within_3_bonds = 0;

  for (int i = 0; i < atoms_at_min_tot_d.number_elements(); ++i) {
    const atom_number_t j = atoms_at_min_tot_d[i];
    for (int k = 0; k < matoms; ++k) {
      if (k == j) {
        continue;
      }

      const int d = dm[matoms * j + k];

      from_middle.extra(d);

      if (d > 3) {
        continue;
      }

      ++atoms_within_3_bonds;

      if (z[k] != 6) {
        ++heteroatoms_within_3_bonds;
      }
    }
  }

  descriptor[iwdescr_stddcentre].set(sqrt(from_middle.variance()));
  descriptor[iwdescr_centre3].set(static_cast<float>(atoms_within_3_bonds));
  descriptor[iwdescr_centre3h].set(static_cast<float>(heteroatoms_within_3_bonds));

  // Historical behaviour: despite the variable name, this counts heteroatoms
  // at distances >= 3. The old comment says this is wrong but useful; preserve it.
  int largest_heteroatoms_within_three_bonds = 0;

  for (int i = 0; i < matoms; ++i) {
    int heteroatoms = (z[i] != 6);

    for (int j = 0; j < matoms; ++j) {
      if (j == i) {
        continue;
      }

      const int d = dm[i * matoms + j];

      if (d < 3) {
        continue;
      }

      if (z[j] != 6) {
        ++heteroatoms;
      }
    }

    if (heteroatoms > largest_heteroatoms_within_three_bonds) {
      largest_heteroatoms_within_three_bonds = heteroatoms;
    }
  }

  descriptor[iwdescr_mh3b].set(
      static_cast<float>(largest_heteroatoms_within_three_bonds));

  m.ring_membership();

  int internal_hydrogen_bond_possibilities = 0;
  for (int i = 0; i < matoms; ++i) {
    if (z[i] == 6) {
      continue;
    } else if (z[i] == 7) {
    } else if (z[i] == 8 || z[i] == 16) {
    } else {
      continue;
    }

    const int hi = m.hcount(i);

    for (int j = i + 1; j < matoms; ++j) {
      const int d = dm[i * matoms + j];
      if (d > longest_internal_hydrogen_bond_separation) {
        continue;
      }
      if (d < shortest_internal_hydrogen_bond_separation) {
        continue;
      }

      if (z[j] == 6) {
        continue;
      } else if (z[j] == 7 || z[j] == 8 || z[j] == 16) {
      } else {
        continue;
      }

      const int hj = m.hcount(j);
      if (hi == 0 && hj == 0) {
        continue;
      }
      if (hi && hj) {  // omit consideration of dual donor/acceptor types
        continue;
      }

      if (AllFlexibleBonds(m, i, j, d, dm)) {
        ++internal_hydrogen_bond_possibilities;
      }
    }
  }

  descriptor[iwdescr_internalhbd].set(
      static_cast<float>(internal_hydrogen_bond_possibilities));

  return max_eccentricity;
}

// Partial symmetry descriptors.
int
IWDescr::IWDescrImpl::NearSymmetricDescriptors(Molecule& m) {
  partial_symmetry::PartialSymmetry psim(m);
  const int* symmetry = psim.SymmetricAtRadius();
#ifdef DEBUG_PARTIAL_SYMMETRY
  Molecule tmp(m);
  write_isotopically_labelled_smiles(tmp, false, cerr);
  cerr << '\n';
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << " " << tmp.smarts_equivalent_for_atom(i) << " value "
         << symmetry[i] << '\n';
  }
#endif

  Accumulator_Int<int> acc;
  const int matoms = m.natoms();
  acc.extra(symmetry, matoms);
  descriptor[iwdescr_maxpsymd] = acc.maxval();
  descriptor[iwdescr_fmaxpsymd] = iwmisc::Fraction<float>(acc.maxval(), matoms);
  descriptor[iwdescr_maxpsymdmean] = acc.average();

  const int nzero = std::count(symmetry, symmetry + matoms, 0);
  descriptor[iwdescr_psymdnumzero] = nzero;

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeSymmetryDescriptors(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy compute_symmetry_related_descriptors(m).
  // Molecule lazily computes symmetry classes and the distance matrix as needed;
  // this descriptor family does not impose ordering dependencies on other
  // descriptor families.
  const int matoms = data.matoms;
  if (matoms == 0) {
    return 1;
  }

  const int* dm = m.distance_matrix_warning_may_change();
  const int* symmetry_class = m.symmetry_classes();

  int asymmetric_atoms = 0;
  int furthest_separated_symmetry_related_atoms = 0;
  int max_equivalent_atoms = 0;
  int longest_path = 0;

  extending_resizable_array<int> class_done;

  for (int i = 0; i < matoms; ++i) {
    const int si = symmetry_class[i];
    if (class_done[si] > 0) {
      continue;
    }

    class_done[si] = 1;

    int atoms_in_class = 1;
    int max_separation_in_class = 0;
    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }

      const int d = dm[i * matoms + j];
      if (d > longest_path) {
        longest_path = d;
      }

      if (si != symmetry_class[j]) {
        continue;
      }

      ++atoms_in_class;
      if (d > max_separation_in_class) {
        max_separation_in_class = d;
      }
    }

    if (atoms_in_class == 1) {
      ++asymmetric_atoms;
      continue;
    }

    if (atoms_in_class > max_equivalent_atoms) {
      max_equivalent_atoms = atoms_in_class;
    }

    if (max_separation_in_class > furthest_separated_symmetry_related_atoms) {
      furthest_separated_symmetry_related_atoms = max_separation_in_class;
    }
  }

  descriptor[iwdescr_symmatom].set(static_cast<float>(matoms - asymmetric_atoms));
  descriptor[iwdescr_fsymmatom].set(
      iwmisc::Fraction<float>(matoms - asymmetric_atoms, matoms));
  descriptor[iwdescr_lsepsymatom].set(
      static_cast<float>(furthest_separated_symmetry_related_atoms));
  if (longest_path > 0) {
    descriptor[iwdescr_flsepsymatom].set(
        iwmisc::Fraction<float>(furthest_separated_symmetry_related_atoms,
                                longest_path));
  } else {
    descriptor[iwdescr_flsepsymatom].set(0.0f);
  }
  descriptor[iwdescr_maxsymmclass].set(static_cast<float>(max_equivalent_atoms));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeChargeDescriptors(Molecule& m, PerMoleculeData& data) {
  if (! descriptors_to_compute.charge_descriptors || ! charge_assigner.active()) {
    return 1;
  }

  // Legacy order from iwdescriptors(): assign charges, store charge-derived
  // descriptors, compute Andrews-Craik-Martin descriptors, then refresh ncon
  // because charge assignment can alter molecule connection state.
  (void) charge_assigner.process(m);

  StoreChargeAssignerResults(m, data);

  ComputeAndrewsCraikMartinDescriptors(m, data);
  data.RefreshConnections(m);
  return 1;
}

int
IWDescr::IWDescrImpl::StoreChargeAssignerResults(Molecule& m, PerMoleculeData& data) {
  // Migrated from legacy store_charge_assigner_results(m, z, atom).
  // This remains IWDescr-owned because it writes Descriptor values.
  (void)m;

  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const Atom** atom = data.atoms();

  int npos = 0;
  int nneg = 0;
  int positive_nitrogen = 0;
  int negative_nitrogen = 0;

  for (int i = 0; i < matoms; ++i) {
    const formal_charge_t fc = atom[i]->formal_charge();
    if (fc == 0) {
      continue;
    }

    if (fc < 0) {
      ++nneg;
      if (z[i] == 7) {
        ++negative_nitrogen;
      }
    } else {
      ++npos;
      if (z[i] == 7) {
        ++positive_nitrogen;
      }
    }
  }

  descriptor[iwdescr_brunsneg].set(static_cast<float>(nneg));
  descriptor[iwdescr_brunspos].set(static_cast<float>(npos));
  descriptor[iwdescr_formal_charge].set(static_cast<int>(nneg + npos));
  descriptor[iwdescr_nplus].set(static_cast<float>(positive_nitrogen));
  descriptor[iwdescr_nminus].set(static_cast<float>(negative_nitrogen));

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeAndrewsCraikMartinDescriptors(Molecule& m,
                                                           PerMoleculeData& data) {
  // Migrated from legacy andrews_craik_martin(m, z, atom). The ACM helper
  // functions are stateless and remain file-local; this member owns the
  // descriptor storage and preserves the historical helper ordering.
  (void)m;

  const int matoms = data.matoms;
  const atomic_number_t* z = data.atomic_numbers();
  const Atom** atom = data.atoms();

  std::vector<int> already_done(matoms, 0);

  float rc = 0.0f;

  rc += IdentifyAcmPhosphoricAcids(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmAcids(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmOH(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmCO(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmOS(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmHalogen(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmC(m, matoms, z, atom, already_done.data());
  rc += IdentifyAcmN(m, matoms, z, atom, already_done.data());

  descriptor[iwdescr_acmbe].set(rc);

  return 1;
}

int
IWDescr::IWDescrImpl::ComputeDonorAcceptorDescriptors(Molecule& m, PerMoleculeData& data) {
  if (! descriptors_to_compute.donor_acceptor || ! donor_acceptor_assigner.active()) {
    return 1;
  }

  if (! data.EnsureDonorAcceptorResults()) {
    return 0;
  }

  donor_acceptor_assigner.process(m, data.donor_acceptor_results_data());

  StoreDonorAcceptorResults(m, data);

  if (min_hbond_feature_separation > 0) {
    ComputeHbondSeparationDescriptors(m, data);
  }

  return 1;
}

int
IWDescr::IWDescrImpl::StoreDonorAcceptorResults(Molecule& m, PerMoleculeData& data) {
  const int matoms = data.matoms;
  const int* da_results = data.donor_acceptor_results_data();

  int acc = 0;
  int dual = 0;
  int don = 0;

  for (int i = 0; i < matoms; ++i) {
    if (da_results[i] == 0) {
      continue;
    } else if (da_results[i] == 1) {
      ++acc;
    } else if (da_results[i] == 2) {
      ++dual;
    } else if (da_results[i] == 3) {
      ++don;
    }
  }

  acc += dual;
  don += dual;

  descriptor[iwdescr_brunsacc].set(static_cast<float>(acc));
  descriptor[iwdescr_brnsdual].set(static_cast<float>(dual));
  descriptor[iwdescr_brunsdon].set(static_cast<float>(don));
  descriptor[iwdescr_brunshbdsum].set(static_cast<float>(acc + don - dual));

  return 1;
}

namespace {

void
ComputeAverage(Accumulator_Int<int>& acc, double& result) {
  if (acc.n() == 0) {
    result = 0.0;
  } else if (acc.n() == 1) {
    result = acc.sum();
  } else {
    result = acc.average();
  }
}

// Track the shortest and next-to-shortest distances in the same way as the
// legacy compute_bond_separation_parameters() helper.
void
CheckAgainstTwo(int tocheck, int& shortest, int& nshortest) {
  if (tocheck < shortest) {
    nshortest = shortest;
    shortest = tocheck;
  } else if (tocheck < nshortest) {
    nshortest = tocheck;
  }
}

void
CheckForNoValues(int& shortest, int& nshortest, int matoms) {
  if (shortest == matoms) {
    shortest = 0;
    nshortest = 0;
  } else if (nshortest == matoms) {
    nshortest = 0;
  }
}

}  // namespace

int
IWDescr::IWDescrImpl::ComputeHbondSeparationDescriptors(Molecule& m,
                                                        PerMoleculeData& data) {
  const int matoms = data.matoms;
  const int* da_results = data.donor_acceptor_results_data();
  if (da_results == nullptr) {
    return 1;
  }

  int shortest_aa = matoms;
  int shortest_ad = matoms;
  int shortest_dd = matoms;

  int nshortest_aa = matoms;
  int nshortest_ad = matoms;
  int nshortest_dd = matoms;

  Accumulator_Int<int> aa;
  Accumulator_Int<int> ad;
  Accumulator_Int<int> dd;

  for (int i = 0; i < matoms; ++i) {
    if (da_results[i] == 0) {
      continue;
    }

    const int i_is_acceptor = da_results[i] >= 2;

    for (int j = i + 1; j < matoms; ++j) {
      if (da_results[j] == 0) {
        continue;
      }

      const int d = m.bonds_between(i, j);
      if (d < min_hbond_feature_separation || d > max_hbond_feature_separation) {
        continue;
      }

      const int j_is_acceptor = da_results[j] >= 2;

      if (i_is_acceptor && j_is_acceptor) {
        aa.extra(d);
        CheckAgainstTwo(d, shortest_aa, nshortest_aa);
      }
      if (! i_is_acceptor && ! j_is_acceptor) {
        dd.extra(d);
        CheckAgainstTwo(d, shortest_dd, nshortest_dd);
      } else {
        ad.extra(d);
        CheckAgainstTwo(d, shortest_ad, nshortest_ad);
      }
    }
  }

  CheckForNoValues(shortest_aa, nshortest_aa, matoms);
  CheckForNoValues(shortest_ad, nshortest_ad, matoms);
  CheckForNoValues(shortest_dd, nshortest_dd, matoms);

  double ave_aa = 0.0;
  double ave_ad = 0.0;
  double ave_dd = 0.0;
  ComputeAverage(aa, ave_aa);
  ComputeAverage(ad, ave_ad);
  ComputeAverage(dd, ave_dd);

  descriptor[iwdescr_aamind].set(static_cast<float>(shortest_aa));
  descriptor[iwdescr_aa2mind].set(static_cast<float>(nshortest_aa));
  descriptor[iwdescr_aaave].set(static_cast<float>(ave_aa));

  descriptor[iwdescr_admind].set(static_cast<float>(shortest_ad));
  descriptor[iwdescr_ad2mind].set(static_cast<float>(nshortest_ad));
  descriptor[iwdescr_adave].set(static_cast<float>(ave_ad));

  descriptor[iwdescr_ddmind].set(static_cast<float>(shortest_dd));
  descriptor[iwdescr_dd2mind].set(static_cast<float>(nshortest_dd));
  descriptor[iwdescr_ddave].set(static_cast<float>(ave_dd));

  return 1;
}

void
IWDescr::IWDescrImpl::ResetDescriptors() {
  if (descriptor == nullptr) {
    return;
  }

  for (int i = 0; i < number_descriptors; ++i) {
    descriptor[i].reset();
  }
}

int
IWDescr::IWDescrImpl::CopyDescriptorsToResults(float* results) const {
  if (descriptor == nullptr || results == nullptr) {
    return 0;
  }

  for (int i = 0; i < number_descriptors; ++i) {
    float v;
    if (descriptor[i].value(v)) {
      results[i] = v;
    } else {
      results[i] = std::numeric_limits<float>::quiet_NaN();
    }
  }

  return 1;
}

int
IWDescr::IWDescrImpl::ReportDescriptorStatistics(std::ostream& output) const {
  for (int i = 0; i < number_descriptors; ++i) {
    if (! descriptor[i].active()) {
      continue;
    }

    descriptor[i].report_statistics(output);
  }

  return output.good();
}
