#include "ring_substitution.h"

#include <algorithm>
#include <vector>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/path.h"
#include "Molecule_Lib/path_around_ring.h"

namespace ring_substitution {

RingSubstitutionGenerator::RingSubstitutionGenerator() {
  _positional_information_only = true;
  _simple_atom_types = false;
  _full_atom_types = false;
  _max_path_length = 24;
  _place_single_feature_bits = false;
}

/*
  Make sure the ring atom types have gaps in them because we
  incrememt by 1 if the ring atom is aromatic
*/

namespace {

static constexpr int kRsTypeDoubleBondOutsideRing = 1;
static constexpr int kRsTypeSpiro = 2;
static constexpr int kRsTypeRingJoin = 3;
static constexpr int kRsTypeTwoSubstituents = 4;
static constexpr int kRsTypeSubstituted = 5;
static constexpr int kRsTypeAnotherRing = 7;

static constexpr int kRsTypeSaturatedCarbon = 9;
static constexpr int kRsTypeUnsaturatedCarbon = 11;
static constexpr int kRsTypeSaturatedHeteroatom = 13;
static constexpr int kRsTypeUnsaturatedHeteroatom = 15;

static constexpr int kRsTypeMethyl = 17;
static constexpr int kRsTypeTerminalNitrogen = 19;
static constexpr int kRsTypeNitro = 21;
static constexpr int kRsTypeSaturatedNitrogen = 23;
static constexpr int kRsTypeSp2Nitrogen = 25;
static constexpr int kRsTypeHydroxy = 27;
static constexpr int kRsTypeEther = 29;
static constexpr int kRsTypeSulphur = 31;
static constexpr int kRsTypeFluorine = 33;
static constexpr int kRsTypeHalogen = 35;

}  // namespace

bool
RingSubstitutionGenerator::Generate(const IWString& mname,
                                    const resizable_array<int>& abstract_path, int* tmp,
                                    Sparse_Fingerprint_Creator& sfpc) const {
  int n = abstract_path.number_elements();

// #define DEBUG_RING_SUBSTITUTION
#ifdef DEBUG_RING_SUBSTITUTION
  cerr << mname << '\n';
  for (int i = 0; i < n; i++) {
    cerr << " i = " << i << " abstract_path " << abstract_path[i] << '\n';
  }
#endif

  assert(n > 2);

  copy_vector(tmp, abstract_path.rawdata(), n);
  copy_vector(tmp + n, abstract_path.rawdata(), n);

  int jstop = n;
  if (jstop > _max_path_length) {
    jstop = _max_path_length;
  }

  for (int i = 0; i < n; i++) {
    int ta = tmp[i];

    if (0 == ta) {
      continue;
    }

    for (int j = 1; j < jstop; j++) {
      int k = j + i;  // index of next atom

      int tb = tmp[k];
      if (0 == tb) {
        continue;
      }

      unsigned int b = 40 + 30 + j * 30 * 30 * _max_path_length +
                       (ta * tb);  // 30 is larger than all atom types assigned

#ifdef DEBUG_RING_SUBSTITUTION
      cerr << " i = " << i << " ar " << ta << " j = " << j << " k = " << k << " ar " << tb
           << " b = " << b << '\n';
#endif
      sfpc.hit_bit(b);
    }
  }

  return true;
}

bool
RingSubstitutionGenerator::Generate(Molecule& m, const int* atype,
                                    const Set_of_Atoms& par, int* tmp,
                                    Sparse_Fingerprint_Creator& sfpc) const {
  resizable_array<int> abstract_path;

  const int n = par.number_elements();

  int number_features = 0;
  int the_feature = -1;

  for (int i = 0; i < n; i++) {
    atom_number_t j = par[i];

    abstract_path.add(atype[j]);

    if (0 == atype[j]) {
      continue;
    }

    if (0 == number_features) {
      the_feature = atype[j];
    }

    number_features++;

    if (_place_single_feature_bits) {
      sfpc.hit_bit(atype[j]);
    }
  }

  // cerr << "Abstract ring contains " << number_features << " features\n";

  if (1 == number_features) {
    sfpc.hit_bit(1383 * n + the_feature);
    return true;
  }

  assert(abstract_path.number_elements() == n);

  return Generate(m.name(), abstract_path, tmp, sfpc);
}

bool
RingSubstitutionGenerator::Generate(Molecule& m, const int* atype,
                                    Sparse_Fingerprint_Creator& sfpc) const {
  const int nr = m.nrings();

  const int matoms = m.natoms();

  std::vector<int> ring_already_done(nr, 0);
  std::vector<int> tmp(matoms + matoms, 0);

  m.compute_aromaticity_if_needed();

  // Do all the non-fused rings first. Includes spiro fused

  int rings_processed = 0;

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    //  cerr << "Ring " << i << " has " <<
    //  ri->largest_number_of_bonds_shared_with_another_ring() << " shared bonds\n";

    if (ri->largest_number_of_bonds_shared_with_another_ring() > 0) {
      continue;
    }

    Generate(m, atype, *ri, tmp.data(), sfpc);

    ring_already_done[i] = 1;
    rings_processed++;
  }

  if (nr == rings_processed) {
    return true;
  }

  // Now any fused rings

  std::vector<int> in_ring_system(matoms, 0);

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i]) {
      continue;
    }

    std::fill(in_ring_system.begin(), in_ring_system.end(), 0);

    const Ring* ri = m.ringi(i);

    ri->set_vector(in_ring_system.data(), 1);

    rings_processed++;

    for (int j = i + 1; j < nr; j++) {
      if (ring_already_done[j]) {
        continue;
      }

      const Ring* rj = m.ringi(j);

      if (rj->fused_system_identifier() != ri->fused_system_identifier()) {
        continue;
      }

      rj->set_vector(in_ring_system.data(), 1);
      ring_already_done[j] = 1;

      rings_processed++;
    }

    Set_of_Atoms s;
    if (!path_around_edge_of_ring_system(m, in_ring_system.data(), 1,
                                         s))  // strongly fused
    {
      sfpc.hit_bit(count_occurrences_of_item_in_array(
          1, matoms,
          in_ring_system.data()));  // unprocessed ring, hit bit according to size
      continue;
    }

    Generate(m, atype, s, tmp.data(), sfpc);

    if (nr == rings_processed) {
      break;
    }
  }

  return true;
}

int
RingSubstitutionGenerator::determine_substitution_type(Molecule& m, atom_number_t zatom,
                                                       const Atom& a) const {
  int acon = a.ncon();
  // cerr << "determine_substitution_type for atom type " <<
  // m.smarts_equivalent_for_atom(zatom) << '\n';

  if (2 == acon) {  // unsubstituted
    return 0;
  }

  if (m.nrings(zatom) > 1) {
    return kRsTypeRingJoin;
  }

  if (4 == acon) {
    return kRsTypeTwoSubstituents;
  }

  for (int i = 0; i < acon; i++) {
    const Bond* b = a[i];

    if (b->nrings()) {
      continue;
    }

    if (b->is_double_bond()) {
      return kRsTypeDoubleBondOutsideRing;
    }

    if (_positional_information_only) {
      return kRsTypeSubstituted;
    }

    atom_number_t o = b->other(zatom);

    if (m.is_ring_atom(o)) {
      return kRsTypeAnotherRing;
    }

    atomic_number_t zo = m.atomic_number(o);

    int ocon = m.ncon(o);

    if (_simple_atom_types) {
      int unsaturation = m.nbonds(o) - ocon;
      if (6 == zo) {
        if (0 == unsaturation) {
          return kRsTypeSaturatedCarbon;
        } else {
          return kRsTypeUnsaturatedCarbon;
        }
      } else if (1 == ocon && (9 == zo || 17 == zo || 35 == zo || 53 == zo)) {
        return kRsTypeHalogen;
      } else if (0 == unsaturation) {
        return kRsTypeSaturatedHeteroatom;
      } else {
        return kRsTypeUnsaturatedHeteroatom;
      }

      continue;
    }

    //  Full atom typing

    if (6 == zo) {
      if (1 == ocon) {
        return kRsTypeMethyl;
      }

      int nbonds = m.nbonds(o);
      if (ocon == nbonds) {
        return kRsTypeSaturatedCarbon;
      } else {
        return kRsTypeUnsaturatedCarbon;
      }
    } else if (7 == zo) {
      if (1 == ocon) {
        return kRsTypeTerminalNitrogen;
      }

      int nbonds = m.nbonds(o);

      if (3 == ocon && 5 == nbonds) {
        return kRsTypeNitro;
      } else if (ocon == nbonds) {
        return kRsTypeSaturatedNitrogen;
      } else {
        return kRsTypeSp2Nitrogen;
      }
    } else if (8 == zo) {
      if (1 == ocon) {
        return kRsTypeHydroxy;
      } else {
        return kRsTypeEther;
      }
    } else if (16 == zo) {
      if (1 == ocon) {
        return kRsTypeHydroxy;
      } else if (2 == ocon) {
        return kRsTypeEther;
      } else {
        return kRsTypeSulphur;  // some other kind of state, who knows
      }
    } else if (9 == zo) {
      return kRsTypeFluorine;
    } else if (17 == zo || 35 == zo || 53 == zo) {
      return kRsTypeHalogen;
    }

    return kRsTypeSubstituted;  // of not classified above
  }

  assert(nullptr == "Could not classify non-ring bond???");

  return 0;
}

namespace {

bool
IsSpiroFused(Molecule& m, atom_number_t zatom, const Atom& a) {
  int acon = a.ncon();

  if (4 != acon) {
    return false;
  }

  if (2 != m.nrings(zatom)) {  // too hard and rare otherwise
    return false;
  }

  for (int i = 0; i < acon; i++) {
    const Bond* b = a[i];

    if (0 == b->nrings()) {
      return false;
    }
  }

  return true;
}

bool
HasDoubleBondOutsideRing(Molecule& m, atom_number_t zatom, const Atom& a) {
  int acon = a.ncon();

  if (3 != acon) {
    return false;
  }

  for (int i = 0; i < acon; i++) {
    const Bond* b = a[i];
    //  cerr << " from atom " << zatom << " have bond to " << b->other(zatom) << " ring "
    //  << b->nrings() << " single? " << b->is_single_bond() << '\n';

    if (b->nrings()) {
      continue;
    }

    if (b->is_double_bond()) {
      return true;
    }
  }

  return false;
}

}  // namespace

bool
RingSubstitutionGenerator::assign_atom_types(Molecule& m, int* atype) const {
  m.ring_membership();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int nr = m.nrings(i);

    if (0 == nr) {
      continue;
    }

    const Atom* a = m.atomi(i);

    if (HasDoubleBondOutsideRing(m, i, *a)) {
      atype[i] = kRsTypeDoubleBondOutsideRing;
      if (m.is_aromatic(i)) {
        atype[i]++;
      }
      continue;
    }

    if (IsSpiroFused(m, i, *a)) {
      atype[i] = kRsTypeSpiro;
      continue;
    }

    atype[i] = determine_substitution_type(m, i, *a);

    if (0 == atype[i])
      ;
    else if (m.is_aromatic(i)) {
      atype[i]++;
    }
  }

  return true;
}

bool
RingSubstitutionGenerator::Generate(Molecule& m, Sparse_Fingerprint_Creator& sfc) const {
  std::vector<int> atype(m.natoms(), 0);

  assign_atom_types(m, atype.data());

  return Generate(m, atype.data(), sfc);
}

}  // namespace ring_substitution
