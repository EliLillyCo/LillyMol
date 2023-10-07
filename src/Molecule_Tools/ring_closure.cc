#include <algorithm>
#include <limits>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/ring_closure.h"

namespace ring_closure {

using std::cerr;

// Count the number of spiro fusions in `m`.
// `ring_atoms` is a temporary array used here.
int
SpiroFusionCount(Molecule& m,
        int* ring_atoms) {
  const int nr = m.nrings();
  if (nr < 2) {
    return 0;
  }

  const int matoms = m.natoms();

  static constexpr int kOne = 1;

  int rc = 0;

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);
    std::fill_n(ring_atoms, matoms, 0);
    ri->set_vector(ring_atoms, kOne);

    for (int j = i + 1; j < nr; ++j) {
      if (m.ringi(j)->count_members_set_in_array(ring_atoms, kOne) == 1) {
        ++rc;
      }
    }
  }

  return rc;
}

MoleculeData::MoleculeData() {
  initial_aromatic_ring_count = 0;
  initially_aromatic = nullptr;
  spiro_count_tmp_array = nullptr;
}

MoleculeData::~MoleculeData() {
  if (initially_aromatic != nullptr) {
    delete [] initially_aromatic;
  }
  if (visited != nullptr) {
    delete [] visited;
  }

  if (spiro_count_tmp_array != nullptr) {
    delete [] spiro_count_tmp_array;
  }
}

int
MoleculeData::Initialise(Molecule& m) {
  initial_aromatic_ring_count = m.aromatic_ring_count();
  const int matoms = m.natoms();
  initially_aromatic = new int[matoms];

  if (initial_aromatic_ring_count == 0) {
    std::fill_n(initially_aromatic, matoms, 0);
  } else {
    for (int i = 0; i < matoms; ++i) {
      initially_aromatic[i] = m.is_aromatic(i);
    }
  }

  spiro_count_tmp_array = new int[matoms];
  initial_spiro_ring_count = SpiroFusionCount(m, spiro_count_tmp_array);

  visited = new int[matoms];

  static constexpr int kHydrogen = 1;

  if (m.natoms(kHydrogen) > 0) {
    temp_detach_atoms.detach_atoms(m, kHydrogen);
  }

  return 1;
}

RingClosureOptions::RingClosureOptions() {
  _min_bonds_between = 2;
  // Let's have a sensible default.
  // _max_bonds_between = std::numeric_limits<int>::max();
  _max_bonds_between = 6;

  _max_variants_formed = std::numeric_limits<uint32_t>::max();

  _molecules_processed = 0;

  _min_distance = 0.0f;
  _max_distance = 0.0f;
    
  _min_bond_angle = 0.0f;
  _max_bond_angle = 0.0f;

  _pairs_failing_distance_constraint = 0;

  _isotope = 0;
  _allow_caged_aromatic_rings = false;
  _allow_aromatic_3_membered_ring_fusions = false;
  _allow_cyclopropene = false;
  _allow_spiro_fusions = true;
  _attempt_aromatization = false;
  _aromatized_forms_found = 0;

  _nd3x2_to_c = 0;
  _od2x2_to_n = 0;

  _spiro_ring_fusions_created = 0;
}

int
TokensToQuery(const google::protobuf::RepeatedPtrField<std::basic_string<char> >& tokens,
              resizable_array_p<Substructure_Query>& queries,
              const int verbose) {
  for (const std::string& token : tokens) {
    if (! process_cmdline_token('C', token, queries, verbose)) {
      cerr << "TokensToQuery:cannot process '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
RingClosureOptions::Initialise(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return 0;
  }

  const int verbose = cl.option_present('v');

  IWString fname = cl.string_value(flag);

  std::optional<RingClosure::ring_closure> proto =
     iwmisc::ReadTextProtoCommentsOK<RingClosure::ring_closure>(fname);
  if (! proto) {
    cerr << "RingClosureOptions::Initialise:cannot open/read '" << fname << "'\n";
    return 0;
  }

  if (proto->has_min_bonds_between()) {
    _min_bonds_between = proto->min_bonds_between();
  }

  if (proto->has_max_bonds_between()) {
    _max_bonds_between = proto->max_bonds_between();
  }

  if (proto->has_max_variants_formed()) {
    _max_variants_formed = proto->max_variants_formed();
  }

  if (! proto->has_allow_spiro_fusions()) {
  } else {
    _allow_spiro_fusions = proto->allow_spiro_fusions();
  }

  if (proto->has_attempt_aromatization()) {
    _attempt_aromatization = proto->attempt_aromatization();
  }

  if (proto->btype_size() > 0) {
    for (uint32_t btype : proto->btype()) {
      switch (btype) {
        case RingClosure::SS_DOUBLE_BOND: {
          _btype << static_cast<bond_type_t>(DOUBLE_BOND);
          break;
        }
        case RingClosure::SS_SINGLE_BOND: {
          _btype << static_cast<bond_type_t>(SINGLE_BOND);
          break;
        }
        default: {
        }
      }
    }
  } else {
    _btype << static_cast<bond_type_t>(SINGLE_BOND);
  }

  _allow_caged_aromatic_rings = proto->allow_caged_aromatic_rings();

  if (! TokensToQuery(proto->connection(), _connection, verbose)) {
    cerr << "RingClosureOptions::Initialise:cannot parse connection specifications\n";
    return 0;
  }

  if (! TokensToQuery(proto->product_must_have(), _product_must_have, verbose)) {
    cerr << "RingClosureOptions::Initialise:cannot parse product_must_have specifications\n";
    return 0;
  }

  if (! TokensToQuery(proto->product_must_not_have(), _product_must_not_have, verbose)) {
    cerr << "RingClosureOptions::Initialise:cannot parse product_must_not_have specifications\n";
    return 0;
  }

  // Should handle the cases where either min or max is unspecified.
  if (proto->has_distance_constraint()) {
    _min_distance = proto->distance_constraint().min();
    _max_distance = proto->distance_constraint().max();

    if (_min_distance > _max_distance) {
      cerr << "RingClosureOptions::Initialise:inconsistent min/max distance\n";
      return 0;
    }
  }

  if (_max_distance > 0.0f) {
  }

  if (proto->has_bond_angle_constraint()) {
    _min_bond_angle = proto->bond_angle_constraint().min() * DEG2RAD;
    _max_bond_angle = proto->bond_angle_constraint().max() * DEG2RAD;

    if (_min_bond_angle > _max_bond_angle) {
      cerr << "RingClosureOptions::Initialise:inconsistent min/max bond angle\n";
      return 0;
    }
  }

  if (proto->has_isotope()) {
    _isotope = proto->isotope();
  }

  if (! proto->has_nd3x2_to_c()) {
  } else {
    _nd3x2_to_c = proto->nd3x2_to_c();
  }

  if (! proto->has_od2x2_to_n()) {
  } else {
    _od2x2_to_n = proto->od2x2_to_n();
  }

  return 1;
}

int
RingClosureOptions::IdentifyRingEnds(Molecule& m,
                        int* ok_atom) {
  Molecule_to_Match target(&m);
  int queries_matching = 0;
  for (Substructure_Query* q : _connection) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    sresults.each_embedding_set_vector(ok_atom, 1);
    ++queries_matching;
  }

  if (queries_matching) {
    return queries_matching;
  }

  return 0;
}

// Return true if `zatom` can get a bond of type `btype`
int
RingClosureOptions::OkImplicitHydrogens(Molecule& m, atom_number_t zatom,
                        bond_type_t btype) const {
  if (btype == SINGLE_BOND) {
    return m.implicit_hydrogens(zatom) > 0;
  }
  if (btype == DOUBLE_BOND) {
    return m.implicit_hydrogens(zatom) > 1;
  }

  return 0;
}

// Return true of a bond separation of `dist` is consistent with
// both _min_bonds_between and _max_bonds_between.
int
RingClosureOptions::OkBondSeparation(const int dist) const {
  if (dist < _min_bonds_between) {
    return 0;
  }

  if (dist > _max_bonds_between) {
    return 0;
  }

  return 1;
}

// Return true if distance `d` is consistent with both
// _min_distance and _max_distance.
// the only reason this is not const is that we update
// _pairs_failing_distance_constraint if we fail.
int
RingClosureOptions::OkDistance(float d) {
  if (_max_distance == 0.0f) {
    return 1;
  }

  if (d < _min_distance) {
    ++_pairs_failing_distance_constraint;
    return 0;
  }

  if (d > _max_distance) {
    ++_pairs_failing_distance_constraint;
    return 0;
  }

  return 1;
}

// Return true if all the bonds defined by *-a1-a2-* are
// consistent with _min_bond_angle and _max_bond_angle.
int
RingClosureOptions::OkBondAngles(const Molecule& m,
                                 atom_number_t a1,
                                 atom_number_t a2) const {
#ifdef DEBUG_OK_BOND_ANGLES
  cerr << "OkBondAngles: btw " << _min_bond_angle << " and " << _max_bond_angle << " atoms " << a1 << " and " << a2 << '\n';
#endif

  if (_max_bond_angle == 0.0f) {
    return 1;
  }

  if (! OkBondAnglesDirectional(m, a1, a2)) {
    return 0;
  }

  return OkBondAnglesDirectional(m, a2, a1);
}

// Return true of all the bond angles defined *-a1-a2
// are consistent with _min_bond_angle and _max_bond_angle.
int
RingClosureOptions::OkBondAnglesDirectional(const Molecule& m,
                                 atom_number_t a1,
                                 atom_number_t a2) const {
#ifdef DEBUG_OK_BOND_ANGLES
  cerr << "Checking bonds defined by " << a1 << " to " << a2 << '\n';
#endif
  for (const Bond* b : m.atom(a1)) {
    atom_number_t j = b->other(a1);
    if (j == a2) {
      continue;
    }

    float angle = m.bond_angle(j, a1, a2);
    if (angle < 0.0f) {  // not sure this is necessary.
      angle = - angle;
    }
#ifdef DEBUG_OK_BOND_ANGLES
    cerr << "Angle is " << angle << " in degrees " << (angle * RAD2DEG) << " atoms " << j << ' ' << a1 << " " << a2 << '\n';
#endif

    if (angle < _min_bond_angle || angle > _max_bond_angle) {
      //cerr << "should reject\n";
      return 0;
    }
  }

  return 1;
}

int
RingClosureOptions::OkQueryRequirements(Molecule& m) {
  Molecule_to_Match target(&m);

  if (! ContainsRequiredMotif(target)) {
    return 0;
  }

  if (_product_must_not_have.empty()) {
    return 1;
  }

  for (Substructure_Query* q : _product_must_not_have) {
    if (q->substructure_search(target)) {
      return 0;
    }
  }

  return 1;
}

int
RingClosureOptions::ContainsRequiredMotif(Molecule_to_Match& target) {
  if (_product_must_have.empty()) {
    return 1;
  }

  for (Substructure_Query* q : _product_must_have) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

resizable_array_p<Molecule>
RingClosureOptions::Process(Molecule& m) {
  resizable_array_p<Molecule> result;

  const int matoms = m.natoms();

  if (matoms < 3) {
    return result;
  }

  // Unless we do this, we run into complications in temp_detach_atoms
  // and replacing implicit hydrogens.
  m.move_hydrogens_to_end_of_connection_table(1);

  // write_isotopically_labelled_smiles(m, false, std::cerr);

  ++_molecules_processed;

  MoleculeData molecule_data;
  molecule_data.Initialise(m);

  std::unique_ptr<int[]> ok_atom = std::make_unique<int[]>(matoms);
  if (_connection.empty()) {
    std::fill_n(ok_atom.get(), matoms, 1);
  } else if (IdentifyRingEnds(m, ok_atom.get())) {
  } else {
    return result;
  }

#ifdef DEBUG_QUERY_MATCHES
  for (int i = 0; i < matoms; ++i) {
    cerr << "Atom " << i << " ok_atom " << ok_atom[i] << ' ' << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

  for (int i = 0; i < matoms; ++i) {
    if (m.implicit_hydrogens(i) == 0) {
      ok_atom[i] = 0;
    }
  }

  // If we are just placing double bonds, we could further constrain ok_atom...

  std::optional<atom_number_t> only_process;
  for (bond_type_t btype : _btype) {
    FormRings(m, ok_atom.get(), only_process, btype, molecule_data, result);
    if (result.size() >= _max_variants_formed) {
      break;
    }
  }

  if (_nd3x2_to_c) {
    for (int i = 0; i < matoms; ++i) {
      if (m.atomic_number(i) == 7 && m.ring_bond_count(i) == 2) {
        Molecule mcopy(m);
        mcopy.set_atomic_number(i, 6);
        mcopy.set_isotope(i, 2);
        only_process = i;
        const auto oksave = ok_atom[i];
        ok_atom[i] = 1;
        FormRings(mcopy, ok_atom.get(), only_process, SINGLE_BOND, molecule_data, result);
        ok_atom[i] = oksave;
        if (result.size() >= _max_variants_formed) {
          break;
        }
      }
    }
  }

  if (_od2x2_to_n) {
    for (int i = 0; i < matoms; ++i) {
      if (m.atomic_number(i) == 8 && m.ring_bond_count(i) == 2) {
        Molecule mcopy(m);
        mcopy.set_atomic_number(i, 7);
        mcopy.set_isotope(i, 2);
        only_process = i;
        const auto oksave = ok_atom[i];
        ok_atom[i] = 1;
        FormRings(mcopy, ok_atom.get(), only_process, SINGLE_BOND, molecule_data, result);
        ok_atom[i] = oksave;
        if (result.size() >= _max_variants_formed) {
          break;
        }
      }
    }
  }

  ++_variants_formed[result.size()];

  return result;
}

// Form as many rings as we can with bond type `btype`.
// `ok_atom` marks the atoms that are OK to participate.
// If `only_process` is set, then we only process that atom in the
// outer loop. That gets used when a given atom gets transformed
// and we want to explore new possibilities.
int
RingClosureOptions::FormRings(Molecule& m,
                        const int* ok_atom,
                        std::optional<atom_number_t> only_process,
                        bond_type_t btype,
                        MoleculeData& molecule_data,
                        resizable_array_p<Molecule>& result) {

  const int matoms = m.natoms();
#ifdef DEBUG_FORM_RINGS
  cerr << "Begin processing " << m.name() << " with " << matoms << " atoms\n";
#endif

  for (int i = 0; i < matoms && result.size() < _max_variants_formed; ++i) {
    if (only_process && i != *only_process) {
      continue;
    }
#ifdef DEBUG_FORM_RINGS
    cerr << "Processing atom " << i << " ok_atom " << ok_atom[i] << " ih " << OkImplicitHydrogens(m, i, btype) << '\n';
#endif
    if (! ok_atom[i] || ! OkImplicitHydrogens(m, i, btype)) {
      continue;
    }
    for (int j = i + 1; j < matoms; ++j) {
      if (! ok_atom[j] || ! OkImplicitHydrogens(m, j, btype)) {
        continue;
      }

      const int bonds_between = m.bonds_between(i, j);

      if (! OkBondSeparation(bonds_between)) {
        continue;
      }

      if (! OkDistance(m.distance_between_atoms(i, j))) {
        continue;
      }

      if (! OkChemistry(m, i, j, btype)) {
        continue;
      }

#ifdef DEBUG_FORM_RINGS
      cerr << "Processing atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " and " << 
              j << ' ' << m.smarts_equivalent_for_atom(j) << " dist " << m.bonds_between(i, j) << '\n';
#endif

      // It is an open question how this should be treated.
      // For now, disallow it, but perhaps revisit. But some
      // of the other filtering functions will need to be changed.
      if (m.in_same_ring(i, j)) {
        continue;
      }
      if (m.in_same_ring_system(i, j)) {
        continue;
      }

      if (InSameAromaticSystem(m, i, j)) {
        continue;
      }

      if (WouldFormDisallowedFused3MemberedRing(m, i, j, bonds_between)) {
        continue;
      }

      if (WouldFormBridgedAromatic(m, i, j)) {
        continue;
      }

      if (WouldFormDisallowedCyloPropene(m, i, j, bonds_between, btype)) {
        continue;
      }

      if (WouldFormDisallowedSpiroFusion(m, i, j, bonds_between, btype)) {
        continue;
      }

      std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);
      mcopy->add_bond(i, j, btype);

      if (! OkQueryRequirements(*mcopy)) {
        continue;
      }

      if (molecule_data.initial_aromatic_ring_count > 0 &&
          mcopy->aromatic_ring_count() < molecule_data.initial_aromatic_ring_count) {
        continue;
      }
#ifdef DEBUG_FORM_RINGS
      cerr << "Aromatic ring count ok\n";
#endif

      if (! _allow_caged_aromatic_rings && ContainsCagedAromaticRing(*mcopy)) {
        continue;
      }

      if (StronglyFusedRingIncludesAromatic(*mcopy, i, j, molecule_data)) {
        continue;
      }

      if (TripleBondInRing(*mcopy)) {
        continue;
      }

      if (! OkBondAngles(*mcopy, i, j)) {
        continue;
      }

      UpdateRingsFormed(bonds_between + 1, btype);

      MaybePlaceIsotopes(*mcopy, i, j);

      MaybeAttemptAromatization(*mcopy, i, j, molecule_data, result);

      MaybeReattachHydrogens(*mcopy, i, j, molecule_data);
      if (AlreadySeen(*mcopy, molecule_data.seen)) {
        continue;
      }

      TrackSpiroFusions(*mcopy, molecule_data);

      result << mcopy.release();
    }
  }

  return 1;
}

void
RingClosureOptions::TrackSpiroFusions(Molecule& m,
                        MoleculeData& molecule_data) {
  if (! _allow_spiro_fusions) {
    return;
  }

  if (m.nrings() < 2) {
    return;
  }

  if (SpiroFusionCount(m, molecule_data.spiro_count_tmp_array) >
      molecule_data.initial_spiro_ring_count) {
    ++_spiro_ring_fusions_created;
  }
}

// A new ring has been formed, update counters.
void
RingClosureOptions::UpdateRingsFormed(int nbonds,
                        bond_type_t btype) {
  ++_ring_size_formed[nbonds];

  switch (btype) {
    case SINGLE_BOND:
      ++_formed_with_bond_type[1];
      break;
    case DOUBLE_BOND:
      ++_formed_with_bond_type[2];
      break;
  }
}

// Return true if `a1` and `a2` are both aromatic and in the
// same ring system.
int
RingClosureOptions::InSameAromaticSystem(Molecule& m,
                        atom_number_t a1,
                        atom_number_t a2) const {
  if (! m.is_aromatic(a1) || ! m.is_aromatic(a2)) {
    return 0;
  }

  // Both are aromatic, are they in the same ring or ring system.

  const int fsid1 = m.fused_system_identifier(a1);
  const int fsid2 = m.fused_system_identifier(a2);

  if (fsid1 == fsid2) {
    return 1;
  }

  // Not in same aromatic system.
  return 0;
}

// Return true if joining `a1` and `a2` would create a bridge bond across
// an aromatic ring or ring system.
// For example, a1 might be in a phenyl ring, and a2 might be a phenol
// atom nearby. They cannot be joined.
// Note that we assume the atoms are not in the same ring.
int
RingClosureOptions::WouldFormBridgedAromatic(Molecule& m,
                        atom_number_t a1,
                        atom_number_t a2) const {
  // If neither is aromatic, nothing to worry about.
  if (! m.is_aromatic(a1) && ! m.is_aromatic(a2)) {
    return 0;
  }

  // If both are aromatic, this does not apply.
  if (m.is_aromatic(a1) && m.is_aromatic(a2)) {
    return 0;
  }

  // If they are both in the same ring system, this test does
  // not apply.
  if (m.fused_system_identifier(a1) == m.fused_system_identifier(a2)) {
    return 0;
  }

  // Make sure a1 is the aromatic atom.
  if (m.is_aromatic(a2)) {
    std::swap(a1, a2);
  }

  // Is a2 bonded to an atom in the same ring as a1.
  for (const Bond* b : m.atom(a2)) {
    const atom_number_t j = b->other(a2);
    if (m.in_same_aromatic_ring(a1, j)) {
      return 1;
    }
  }

  // No problems detected.
  return 0;
}

// If _isotope is set, place that isotope at `a1` and `a2`.
int
RingClosureOptions::MaybePlaceIsotopes(Molecule& m,
                                       atom_number_t a1,
                                       atom_number_t a2) const {
  if (_isotope == 0) {
    return 0;
  }

  m.set_isotope(a1, _isotope);
  m.set_isotope(a2, _isotope);

  return 1;
}

int
RingClosureOptions::OkChemistry(Molecule& m, atom_number_t atom1, atom_number_t atom2, 
                        bond_type_t btype) const {
  const atomic_number_t z1 = m.atomic_number(atom1);
  const atomic_number_t z2 = m.atomic_number(atom2);

  // No peroxides.
  if (z1 == 8 && z2 == 8) {
    return 0;
  }

  if (z1 == 16 && z2 == 16) {
    return 0;
  }

  if (z1 == 15 && z2 == 16) {
    return 0;
  }

  // No N-O bonds
  if (z1 == 7 && z2 == 8) {
    return 0;
  }

  if (z1 == 8 && z2 == 7) {
    return 0;
  }

  if (z1 == 7 && z2 == 7) {
    return 0;
  }

  if (btype == SINGLE_BOND) {
    return 1;
  }

  // If both fully saturated, a double bond is OK.
  if (m.saturated(atom1) && m.saturated(atom2)) {
    return 1;
  }

  // One of the atoms is unsaturated, no double bonds.

  return 0;
}

// Return true if there is a triple bond in a ring.
int
RingClosureOptions::TripleBondInRing(Molecule& m) const {
  for (const Bond* b : m.bond_list()) {
    if (! b->is_triple_bond()) {
      continue;
    }
    if (m.ring_bond_count(b->a1()) > 0 ||
        m.ring_bond_count(b->a2()) > 0) {
      return 1;
    }
  }

  return 0;
}

int
RingClosureOptions::AlreadySeen(Molecule& m,
                IW_STL_Hash_Set& seen) const {
  if (seen.contains(m.unique_smiles())) {
    return 1;
  }

  seen.insert(m.unique_smiles());
  return 0;
}

// Return true if `m` contains any strongly fused aromatic.
int
RingClosureOptions::ContainsCagedAromaticRing(Molecule& m) {
  m.compute_aromaticity_if_needed();

  for (const Ring* r : m.sssr_rings()) {
    if (! r->is_aromatic()) {
      continue;
    }
#ifdef DEBUG_CONTAINS_CATED_AROMATIC_RING
    cerr << " arom ring " << r->largest_number_of_bonds_shared_with_another_ring() << '\n';
#endif
    if (r->largest_number_of_bonds_shared_with_another_ring() > 1) {
#ifdef DEBUG_CONTAINS_CATED_AROMATIC_RING
      cerr << "Caged aromatic found\n";
#endif
      return 1;
    }
  }
#ifdef DEBUG_CONTAINS_CATED_AROMATIC_RING
  cerr << "NO caged aromatic rings\n";
#endif

  return 0;
}

// For all unset members of `visited` that have the same fused
// system id as `zatom`, set those array members to `flag`.
// Return the number of atoms set.
int
FillSameFsid(Molecule& m, atom_number_t zatom, int* visited, int flag) {
  const int fsid = m.fused_system_identifier(zatom);

  int rc = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (visited[i] > 0) {
      continue;
    }

    if (m.fused_system_identifier(i) != fsid) {
      continue;
    }

    visited[i] = flag;
    ++rc;
  }

  return rc;
}

// #define DEBUG_ENCOUNTERS_TARGET
// Returns true if starting with `zatom` we can find a path
// to an atom where `visited` == target.
// `target` should not be 1 since we use that here.

int
EncountersTarget(const Molecule& m,
                atom_number_t zatom,
                int* visited,
                int target) {
#ifdef DEBUG_ENCOUNTERS_TARGET
  cerr << "EncountersTarget from atom " << zatom << " value " << visited[zatom] << " target " << target << '\n';
#endif
  for (const Bond* b : m.atom(zatom)) {
    atom_number_t j = b->other(zatom);
#ifdef DEBUG_ENCOUNTERS_TARGET
    cerr << " from atom " << zatom << " to " << j << " visited " << visited[j] << '\n';
#endif

    if (visited[j] == target) {
      return 1;
    }

    if (visited[j]) {
      continue;
    }

    visited[j] = 1;
    if (EncountersTarget(m, j, visited, target)) {
      return 1;
    }
    visited[j] = 0;
  }

  return 0;
}

// We have placed a bond btw `a1` and `a2`. If either one of those
// atoms is aromatic, have we created a strongly fused system involving
// an aromatic ring.
int
RingClosureOptions::StronglyFusedRingIncludesAromatic(Molecule& m,
                atom_number_t a1,
                atom_number_t a2,
                MoleculeData& molecule_data) {
  // If neither is aromatic, no problem.
  if (! m.is_aromatic(a1) && ! m.is_aromatic(a2)) {
    return 0;
  }

  // At least one is aromatic, look for a strongly fused ring involving a1 and a2.
  for (const Ring* r : m.sssr_rings()) {
    // Ignore isolated or singly fused.
    if (r->largest_number_of_bonds_shared_with_another_ring() < 2) {
      continue;
    }

    if (r->contains_bond(a1, a2)) {
      return 1;
    }
  }

  // No problems found.
  return 0;
}

// Return true if a1 and a2 are two bonds apart, one of them is in an aromatic
// ring, and joining them would form a 3 membered, fused ring.
int
RingClosureOptions::WouldFormDisallowedFused3MemberedRing(Molecule& m,
                atom_number_t a1, atom_number_t a2, int bonds_between) const {
  if (_allow_aromatic_3_membered_ring_fusions) {
    return 0;
  }

  if (bonds_between != 2) {
    return 0;
  }

  const int arom1 = m.is_aromatic(a1);
  const int arom2 = m.is_aromatic(a2);

  if (arom1 && arom2) {  // not sure if this could happen.
    return 0;
  }

  if (! arom1 && ! arom2) {
    return 0;
  }

  // Make sure a1 is the aromatic atom.
  if (arom2) {
    std::swap(a1, a2);
  }

  // Identify the atom that is in between a1 and a2.
  // If it is aromatic, we reject this bonding.
  for (const Bond* b : m.atom(a1)) {
    atom_number_t j = b->other(a1);
    if (! m.are_bonded(j, a2)) {
      continue;
    }
    if (m.is_aromatic(j)) {
      return 1;
    }
  }

  return 0;
}

int
RingClosureOptions::WouldFormDisallowedCyloPropene(Molecule& m,
                atom_number_t a1, atom_number_t a2, int bonds_between,
                bond_type_t btype) const {
  if (_allow_cyclopropene) {
    return 0;
  }

  if (bonds_between != 2) {
    return 0;;
  }

  if (btype == DOUBLE_BOND) {
    return 1;
  }

  // We are placing a single bond between a1 and a2. Both must be
  // fully saturated.
  if (! m.saturated(a1) || ! m.saturated(a2)) {
    return 1;
  }

  return 0;
}

// Return 1 if forming a bond between `a1` and `a2` would result in a
// spiro fusion.
int
WouldBeSpiroFusion(Molecule& m,
                   atom_number_t a1,
                   atom_number_t a2) {
  int rbc1 = m.ring_bond_count(a1);
  int rbc2 = m.ring_bond_count(a2);

  // If neither is in a ring, no spiro.
  if (rbc1 == 0 && rbc2 == 0) {
    return 0;
  }

  // If either one already in more than 1 ring, no spiro.
  if (rbc1 > 2 || rbc2 > 2) {
    return 0;
  }

  if (m.in_same_ring(a1, a2)) {
    return 0;
  }

  // One must be 3 connected. Ensure a1 is 3 connected. A2 might also be...
  if (m.ncon(a1) == 3 && rbc1 == 2) {
  } else if (m.ncon(a2) == 3 && rbc2 == 2) {
    std::swap(a1, a2);
    std::swap(rbc1, rbc2);
  } else {
    return 0;
  }

  return 1;
}

// Return 1 if there would be a disallowed spiro form generated
int
RingClosureOptions::WouldFormDisallowedSpiroFusion(Molecule& m,
                atom_number_t a1, atom_number_t a2, int bonds_between,
                bond_type_t btype) const {
  // If we are not forming spiro fusions, notify it this looks like a spiro form.
  if (! _allow_spiro_fusions) {
    return WouldBeSpiroFusion(m, a1, a2);
  }

  // No restrictions active.
  return 0;
}

int
RingClosureOptions::MaybeAttemptAromatization(Molecule& m, atom_number_t a1, atom_number_t a2,
                MoleculeData& molecule_data,
                resizable_array_p<Molecule>& result) {
  if (! _attempt_aromatization) {
    return 0;
  }

  if (m.is_aromatic(a1) || m.is_aromatic(a2)) {
    return 0;
  }

  // We only process isolated rings.
  if (m.fused_system_size(a1) != 1) {
    return 0;
  }

  if (m.ring_bond_count(a1) != 2 ||
      m.ring_bond_count(a2) != 2) {
    return 0;
  }

  const Ring* r = m.ring_containing_atom(a1);
  assert(r->contains(a2));
  if (r->size() == 3) {
    return 0;
  }

  for (int i = 0; i < r->number_elements(); ++i) {
    if (m.ncon(r->item(i)) == 4) {
      return 0;
    }
  }
#ifdef FIX_RING_ITERATOR_PROBLEM
  for (atom_number_t i : *r) {
    if (m.ncon(i) == 4) {
      return 0;
    }
  }
#endif

  std::fill_n(molecule_data.visited, m.natoms(), 0);
  r->set_vector(molecule_data.visited, 1);

  std::unique_ptr<Molecule> mcopy = std::make_unique<Molecule>(m);

  mcopy->set_implicit_hydrogens_known(a1, 0);
  mcopy->set_implicit_hydrogens_known(a2, 0);

  if (! mcopy->find_kekule_form(molecule_data.visited)) {
    return 0;
  }

  MaybeReattachHydrogens(*mcopy, a1, a2, molecule_data);
  if (AlreadySeen(*mcopy, molecule_data.seen)) {
    return 0;
  }

  ++_aromatized_forms_found;

  result << mcopy.release();

  return 1;
}

// Note that this only works because we have moved Hydrogens to the
// end of the connection table. If we don't do that, the atom numbers
// may change during ReattachAtoms, since it removes Hydrogen atoms no
// longer needed.
int
RingClosureOptions::MaybeReattachHydrogens(Molecule& m, 
                atom_number_t a1, 
                atom_number_t a2,
                MoleculeData& molecule_data) {
  // cerr << "temp_detach_atoms active " << molecule_data.temp_detach_atoms.active() << '\n';
  if (! molecule_data.temp_detach_atoms.active()) {
    return 0;
  }

  // First reattach, but do not attach anything to a1 or a2.
  std::fill_n(molecule_data.visited, m.natoms(), 1);
  molecule_data.visited[a1] = 0;
  molecule_data.visited[a2] = 0;

  molecule_data.temp_detach_atoms.ReattachAtoms(m, molecule_data.visited);
  // cerr << "After reattach " << m.smiles() << '\n';

  Make_Implicit_Hydrogens_Explicit mihe;
  if (m.highest_coordinate_dimensionality() == 3) {
    mihe.set_dimensionality(3);
  }

  mihe.set_atom(a1);
  m.make_implicit_hydrogens_explicit(mihe);
  mihe.set_atom(a2);
  m.make_implicit_hydrogens_explicit(mihe);

  // Now reattach, excluding a1 and a2.

  return 1;
}

int
RingClosureOptions::Report(std::ostream& output) const {
  output << "RingClosureOptions::Report:processed " << _molecules_processed << " molecules\n";
  int molecules_formed = 0;
  for (int i = 0; i < _variants_formed.number_elements(); ++i) {
    if (_variants_formed[i] > 0) {
      output << _variants_formed[i] << " molecules formed " << i << " variants\n";
      molecules_formed += (_variants_formed[i] * i);
    }
  }

  output << "Generated " << molecules_formed << " molecules\n";
  for (int i = 2; i < _ring_size_formed.number_elements(); ++i) {
    if (_ring_size_formed[i]) {
      output << _ring_size_formed[i] << " rings of size " << i << " formed\n";
    }
  }

  for (int i = 1; i < _formed_with_bond_type.number_elements(); ++i) {
    if (_formed_with_bond_type[i]) {
      output << _formed_with_bond_type[i] << " rings formed with bond type " << i << '\n';
    }
  }

  if (_attempt_aromatization) {
    output << _aromatized_forms_found << " aromatic forms found\n";
  }

  if (_max_distance > 0.0f) {
    output << _pairs_failing_distance_constraint << " pairs of atoms failed 3D constraints\n";
  }

  output << _spiro_ring_fusions_created << " spiro ring fusions formed\n";

  return output.good();
}

}  // namespace ring_closure
