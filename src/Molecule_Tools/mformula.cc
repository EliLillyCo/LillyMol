#include <algorithm>
#include <cmath>
#include <limits>

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
// constexpr int kMFOther = 17;  in header.

// in some cases, we want to just check element counts.
// for (int i = 0; i < kLastElement; ++i)
constexpr int kLastElement = 13;

void
MFormula::ZeroCountArray() {
  std::fill_n(_count, kMFOther + 1, 0);

  _natoms = 0;
}

MFormula::MFormula() {
  ZeroCountArray();

  _initialised = 0;

  _consider_aromatic = 1;
}

int
MFormula::Build(Molecule& m) {
  if (_consider_aromatic) {
    m.compute_aromaticity_if_needed();
  }

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
  if (z == 6){
    if (! _consider_aromatic || ! m.is_aromatic(i)) {
      ++_count[kMFCarbon];
    } else {
      ++_count[kMFArCarbon];
    }
  } else if (z == 7) {
    if (! _consider_aromatic || ! m.is_aromatic(i)) {
      ++_count[kMFNitrogen];
    } else {
      ++_count[kMFArNitrogen];
    }
  } else if (z == 8) {
    if (! _consider_aromatic || ! m.is_aromatic(i)) {
      ++_count[kMFOxygen];
    } else {
      ++_count[kMFArOxygen];
    }
  } else if (z == 9) {
    ++_count[kMFFluorine];
  } else if (z == 15) {
    ++_count[kMFPhosphorus];
  } else if (z == 16) {
    if (! _consider_aromatic || ! m.is_aromatic(i)) {
      ++_count[kMFSulphur];
     } else {
      ++_count[kMFArSulphur];
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
    } else if (_consider_aromatic && m.is_aromatic(i)) {
      ++_count[kHydrogenAromatic];
    } else {
      ++_count[kHydrogen];
    }
  }

  if (! _consider_aromatic) {
    return 1;
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

void
MFormula::ApplyLogCountScaling(Sparse_Fingerprint_Creator& sfc) const {
  const float denom = std::log(_log_d);

  for (int i = 0; i < kMFOther; ++i) {
    if (_count[i] == 0) {
      continue;
    }
    int c = static_cast<int>(_log_s * std::log(_count[i] + 1) / denom);
    // We don't worry about overflow in `c`, `sfc` takes care of that.
    sfc.hit_bit(i, c);
  }
}

int
MFormula::ToSparseFingerprint(IWString& destination) const {
  Sparse_Fingerprint_Creator sfc;
  if (_log_s > 0.0f) {
    ApplyLogCountScaling(sfc);
  } else {
    for (int i = 0; i <= kMFOther; ++i) {
      if (_count[i] > 0) {
        sfc.hit_bit(i, _count[i]);
      }
    }
  }

  return sfc.daylight_ascii_form_with_counts_encoded(destination);
}

int
MFormula::ToFixedCountedFingerprint(int* destination, int n) const {
  if (destination == nullptr || n < kMFOther + 1) {
    return 0;
  }

  std::fill_n(destination, n, 0);

  if (_log_s > 0.0f) {
    const float denom = std::log(_log_d);
    for (int i = 0; i <= kMFOther; ++i) {
      if (_count[i] == 0) {
        continue;
      }
      int c = static_cast<int>(_log_s * std::log(_count[i] + 1) / denom);
      if (c > std::numeric_limits<uint8_t>::max()) {
        c = std::numeric_limits<uint8_t>::max();
      }
      destination[i] = c;
    }
  } else {
    for (int i = 0; i <= kMFOther; ++i) {
      destination[i] = std::min(_count[i], static_cast<int>(std::numeric_limits<uint8_t>::max()));
    }
  }

  return 1;
}

constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

constexpr char kOpenParen = '(';
constexpr char kCloseParen = ')';

static bool skip_char[256];
static bool need_to_initialise_skip_char = true;

void
InitialiseSkipChar() {
  std::fill_n(skip_char, 256, false);

  skip_char['('] = true;
  skip_char[')'] = true;
  skip_char['+'] = true;
  skip_char['-'] = true;
  skip_char['='] = true;
  skip_char['#'] = true;
  skip_char['@'] = true;
  // Explicit H or H inside square brackets always skipped.
  skip_char['H'] = true;

  for (int i = 0; i < 9; ++i) {
    skip_char['0' + i] = true;
  }

  need_to_initialise_skip_char = false;
}

int
MFormula::Build(const IWString& smiles) {
  if (need_to_initialise_skip_char) {
    InitialiseSkipChar();
  }

  // We do not know anything about aromaticity or ring membership.
  _consider_aromatic = 0;

  bool in_square_bracket = 0;
  bool got_element_in_square_bracket = false;
  const int nchars = smiles.length();

  for (int i = 0; i < nchars; ++i) {
    char c = smiles[i];
    if (skip_char[static_cast<int>(c)]) {
      continue;
    }

    if (kOpenSquareBracket == c) {
      in_square_bracket = true;
      got_element_in_square_bracket = false;
    } else if (kCloseSquareBracket == c) {
      in_square_bracket = false;
    } else if (in_square_bracket && got_element_in_square_bracket) {
    } else {
      if (c == 'C') {
        if (i != nchars - 1 && smiles[i + 1] == 'l') {
          ++_count[kMFChlorine];
          ++i;
        } else {
          ++_count[kMFCarbon];
        }
      } else if (c == 'c') {
        ++_count[kMFCarbon];
      } else if (c == 'N' || c == 'n') {
        ++_count[kMFNitrogen];
      } else if (c == 'O' || c == 'o') {
        ++_count[kMFOxygen];
      } else if (c == 'F') {
        ++_count[kMFFluorine];
      } else if (c == 'P') {
        ++_count[kMFPhosphorus];
      } else if (c == 'S' || c == 's') {
        ++_count[kMFSulphur];
      } else if (c == 'B' && i != nchars - 1 && smiles[i + 1] == 'r') {
        ++_count[kMFBromine];
      } else if (c == 'I') {
        ++_count[kMFIodine];
      } else {
        continue;
      }
      if (in_square_bracket && ! got_element_in_square_bracket) {
        got_element_in_square_bracket = true;
      }
    }
  }

  _natoms = 0;
  for (int i = 0; i < kMFOther; ++i) {
    _natoms += _count[i];
  }

  if (_natoms) {
    _initialised = 1;
  }

  return _natoms;
}

bool
MFormula::IsSubset(const MFormula& rhs) const {
  for (int i = 0; i < kMFOther; ++i) {
    if (_count[i] > rhs._count[i]) {
      return false;
    }
  }

  return true;
}
bool
MFormula::IsElementCountSubset(const MFormula& rhs) const {
  for (int i = 0; i <kLastElement; ++i) {
    if (_count[i] > rhs._count[i]) {
      return false;
    }
  }

  return true;
}

int
MFormula::Carbon() const {
  return _count[kMFCarbon];
}
int
MFormula::Nitrogen() const {
  return _count[kMFNitrogen];
}
int
MFormula::Oxygen() const {
  return _count[kMFOxygen];
}
int
MFormula::Fluorine() const {
  return _count[kMFFluorine];
}
int
MFormula::Phosphorus() const {
  return _count[kMFPhosphorus];
}
int
MFormula::Sulphur() const {
  return _count[kMFSulphur];
}
int
MFormula::Chlorine() const {
  return _count[kMFChlorine];
}
int
MFormula::Bromine() const {
  return _count[kMFBromine];
}
int
MFormula::Iodine() const {
  return _count[kMFIodine];
}

}  // namespace mformula
