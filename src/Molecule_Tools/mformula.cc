#include <algorithm>

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"

#include "mformula.h"

namespace mformula {

constexpr int kMFCarbon = 0;
constexpr int kMFArCarbon = 1;
constexpr int kMFNitrogen = 2;
constexpr int kMFArNitrogen = 3;
constexpr int kMFOxygen = 4;
constexpr int kMFArOxygen = 5;
constexpr int kMFFluorine = 6;
constexpr int kMFPhosphorus = 7;
constexpr int kMFSulphur = 8;
constexpr int kMFArSulphur = 9;
constexpr int kMFChlorine = 10;
constexpr int kMFBromine = 11;
constexpr int kMFIodine = 12;
constexpr int kHydrogenOnHeteroatom = 13;
constexpr int kHydrogenAromatic = 14;
constexpr int kHydrogen = 15;
constexpr int kRingAtom = 16;
// constexpr int kMFOther = 17;

void
MFormula::ZeroCountArray() {
  std::fill_n(_count, kMFOther + 1, 0);
}

MFormula::MFormula() {
  ZeroCountArray();

  _initialised = 0;
}

int
MFormula::Build(Molecule& m) {
  m.compute_aromaticity_if_needed();

  ZeroCountArray();

  for (int i = 0; i < m.natoms(); ++i) {
    Build(m, i);
  }

  _initialised = 1;

  return 1;
}

// Build a particular atom
int
MFormula::Build(Molecule& m, atom_number_t i) {
  atomic_number_t z = m.atomic_number(i);
  if (z == 6) {
    if (m.is_aromatic(i)) {
      ++_count[kMFArCarbon];
    } else {
      ++_count[kMFCarbon];
    }
  } else if (z == 7) {
    if (m.is_aromatic(i)) {
      ++_count[kMFArNitrogen];
    } else {
      ++_count[kMFNitrogen];
    }
  } else if (z == 8) {
    if (m.is_aromatic(i)) {
      ++_count[kMFArOxygen];
    } else {
      ++_count[kMFOxygen];
    }
  } else if (z == 9) {
    ++_count[kMFFluorine];
  } else if (z == 15) {
    ++_count[kMFPhosphorus];
  } else if (z == 16) {
    if (m.is_aromatic(i)) {
      ++_count[kMFArSulphur];
    } else {
      ++_count[kMFSulphur];
    }
  } else if (z == 17) {
    ++_count[kMFChlorine];
  } else if (z == 35) {
    ++_count[kMFBromine];
  } else if (z == 53) {
    ++_count[kMFIodine];
  } else {
    ++_count[kMFOther];
  }

  if (int hcount = m.hcount(i); hcount > 0) {
    if (z != 6) {
      ++_count[kHydrogenOnHeteroatom];
    } else if (m.is_aromatic(i)) {
      ++_count[kHydrogenAromatic];
    } else {
      ++_count[kHydrogen];
    }
  }

  if (int rbc = m.ring_bond_count(i); rbc > 0) {
    ++_count[kRingAtom];
  }

  return 1;
}

int
MFormula::Build(Molecule& m, const Set_of_Atoms& embedding) {
  ZeroCountArray();
  for (atom_number_t a : embedding) {
    Build(m, a);
  }

  _initialised = 1;

  return 1;
}

uint32_t
MFormula::Diff(const MFormula& rhs) const {
  assert(_initialised);
  assert(rhs._initialised);

  uint32_t rc = 0;
  for (int i = 0; i <= kMFOther; ++i) {
    if (_count[i] == rhs._count[i]) {
    } else if (_count[i] < rhs._count[i]) {
      rc += rhs._count[i] - _count[i];
    } else {
      rc += _count[i] - rhs._count[i];
    }
  }

  return rc;
}

int
MFormula::ToSparseFingerprint(IWString& destination) const {
  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i <= kMFOther; ++i) {
    if (_count[i] > 0) {
      sfc.hit_bit(i, _count[i]);
    }
  }

  return sfc.daylight_ascii_form_with_counts_encoded(destination);
}

}  // namespace mformula
