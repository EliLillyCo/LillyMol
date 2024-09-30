#include <iostream> 
#include <random> 

#include "Molecule_Lib/molecule.h"

#include "safe_generate_lib.h"

namespace safe_generate {

using std::cerr;

constexpr char kPercent = '%';
constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';


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
// constexpr int kMFOther = 13;

void
MFormula::ZeroCountArray() {
  std::fill_n(_count, kMFOther, 0);
}

MFormula::MFormula() {
  ZeroCountArray();
}

int
MFormula::Build(Molecule& m) {
  m.compute_aromaticity_if_needed();

  ZeroCountArray();

  for (int i = 0; i < m.natoms(); ++i) {
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
  }

  return 1;
}

int
MFormula::Diff(const MFormula& rhs) const {
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
IsSulfonamide(const Molecule& m, atom_number_t s, atom_number_t n) {
  assert(m.atomic_number(s) == kSulphur);

  for (const Bond* b : m[s]) {
    atom_number_t o = b->other(s);
    if (o == n) {
      continue;
    }

    if (b->is_double_bond() && m.atomic_number(o) == 8) {
      return 1;
    }
  }

  return 0;
}

int
BondsOk(Molecule& m) {
  m.compute_aromaticity_if_needed();

  for (const Bond* b : m.bond_list()) {
    // We don't form new rings, so ignore all ring bonds.
    if (b->nrings()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    // Bonds are only formed with isotopes
    if (m.isotope(a1) == 0 || m.isotope(a2) == 0) {
      continue;
    }

    atomic_number_t z1 = m.atomic_number(a1);
    atomic_number_t z2 = m.atomic_number(a2);
    // Carbon-Carbon always OK.
    if (z1 == 6 && z2 == 6) {
      continue;
    }

    // All NN or OO bonds fail
    if (z1 == 7 && z2 == 7) {
      return 0;
    }
    if (z1 == 8 && z2 == 8) {
      return 0;
    }

    // One atom must be a carbon.
    if (z1 > z2) {
      std::swap(a1, a2);
      std::swap(z1, z2);
    }

    // All #6-F bonds are OK.
    if (z2 == 9) {
      if (z1 == 6) {
        continue;
      } else {
        return 0;
      }
    }

    // #6-#7 and #6-#8 are OK.
    if (z2 == 7 || z2 == 8) {
      continue;
    }

    // reject alkyl halides.
    if (z2 == 17 || z2 == 35 || z2 == 53) {
      if (m.is_aromatic(a1)) {  // aromatic attachments OK.
        continue;
      } else {
        return 0;
      }
    }

    if (z2 == kSulphur) {
      if (z1 == 7) {
        if (IsSulfonamide(m, a2, a1)) {
          continue;
        } else {
          return 0;
        }
      }

      if (z1 == 6 || z1 == 8) {
        continue;
      } else {
        return 0;
      }
    }

    // P-O bonds are OK.
    if (z2 == 15) {
      if (z1 == 8 || z1 == 6) {
        continue;
      } else {
        return 0;
      }
    }

    return 0;
  }

  return 1;
}

SafeFragment::SafeFragment() {
  _natoms = 0;
  _nrings = 0;
  _ncon = 0;
  _distance = 0;
  _ok_to_select = 1;
}

int
SafeFragment::DebugPrint(std::ostream& output) const {
  output << "SafeFragment smiles " << _smiles << " natoms " << _natoms << " ncon " << _ncon << '\n';
  output << "first";;
  for (int f : _first_digit) {
    output << ' ' << f << " '" << _smiles[f] << "'";
  }
  output << '\n';

  return output.good();
}

int
SafeFragment::Build(const const_IWSubstring& smi) {
  const int nchars = smi.length();
  // cerr << "SafeFragment::Build from '" << smi << "'\n";
  // Must have at least an atomic symbol followed by % and a two digit ring number.
  if (nchars < 4) {
    cerr << "SafeFragment::Build:too short '" << smi << "'\n";
    return 0;
  }

  // Construct a smiles without the extra %nn ring openings.
  IWString tmp;
  tmp.reserve(smi.length());

  for (int i = 0; i < nchars; ++i) {
    const char c = smi[i];
    if (c == kPercent) {
      int n = (smi[i + 1] - '0') * 10 + smi[i + 2] - '0';
      _ring.push_back(n);
      _first_digit << (i + 1);
      i += 2;
    } else if (c == ' ') {
      break;
    } else {
      tmp << c;
    }
  }

  if (! _m.build_from_smiles(tmp)) {
    cerr << "SafeFragment::Build:invalid smiles '" << tmp << "' from '" << smi << "'\n";
    return 0;
  }

  _smiles = smi;
  _smiles.truncate_at_first(' ');
  _natoms = _m.natoms();
  _nrings = _m.nrings();
  _ncon = _first_digit.number_elements();

  const int niso = _m.number_isotopic_atoms();
  if (niso == 0) {
    cerr << "SafeFragment::Build:no isotopes '" << _smiles << "' from '" << smi << "'\n";
    return 0;
  }

  if (niso > 1) {
    _distance = _natoms + 1;
    for (int i = 0; i < _natoms; ++i) {
      if (_m.isotope(i) == 0) {
        continue;
      }
      for (int j = i + 1; j < _natoms; ++j) {
        if (_m.isotope(j) == 0) {
          continue;
        }
        const int d = _m.bonds_between(i, j);
        if (d < _distance) {
          _distance = d;
        }
      }
    }
  }

  _mformula.Build(_m);

  return 1;
}

int
SafeFragment::ProcessSquareBracket(const IWString& smi, int &i, int& next_ring,
                     IWString& destination) {
  assert(smi[i] == kOpenSquareBracket);
  destination << kOpenSquareBracket;
  ++i;

  // Should not happen.
  if (i == smi.length()) {
    return 0;
  }

  int got_isotope = 0;
  if (std::isdigit(smi[i])) {
    got_isotope = 1;
  }

  for ( ;i < smi.length(); ++i) {
    destination << smi[i];
    if (smi[i] == kCloseSquareBracket) {
      break;
    }
  }

  if (got_isotope) {
    // Skip over any ring closures already present.
    for ( ; (i + 1) < smi.length() && std::isdigit(smi[i + 1]); ++i) {
      destination << smi[i + 1];
    }
    destination << '%';
    _first_digit << destination.length();
    destination << next_ring;
    ++next_ring;
  }

  return 1;
}

int
SafeFragment::Build(const dicer_data::DicerFragment& proto) {
  IWString smi = proto.smi();
  if (smi.empty()) {
    cerr << "SafeFragment::Build:empty smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (! _m.build_from_smiles(smi)) {
    cerr << "SafeFragment::Build:invalid smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _m.set_name(proto.par());

  _mformula.Build(_m);

  _natoms = _m.natoms();
  // cerr << "SafeFragment::build:smiles " << smi << " natoms " << _natoms << '\n';

  _ncon = _m.number_isotopic_atoms();
  if (_ncon == 0) {
    cerr << "SafeFragment::Build:no isotopes '" << smi << "' from '" <<
            proto.ShortDebugString() << "'\n";
    return 0;
  }

  if (_ncon > 1) {
    _distance = _natoms + 1;
    for (int i = 0; i < _natoms; ++i) {
      if (_m.isotope(i) == 0) {
        continue;
      }
      for (int j = i + 1; j < _natoms; ++j) {
        if (_m.isotope(j) == 0) {
          continue;
        }
        const int d = _m.bonds_between(i, j);
        if (d < _distance) {
          _distance = d;
        }
      }
    }
  }

  // The smiles that is read does not have the %nn ring opening/closing motifs.
  // We insert those after the isotopic labels.

  IWString tmp;
  tmp.reserve(smi.length() + 9);
  int next_ring = 11;
  for (int i = 0; i < smi.length(); ++i) {
    const char c = smi[i];
    if (c == kOpenSquareBracket) {
      ProcessSquareBracket(smi, i, next_ring, tmp);
    } else {
      tmp << c;
    }
  }

  // cerr << "from " << smi << " transform to '" << tmp << " atoms " << _natoms << '\n';
  _smiles = tmp;

  return 1;
}

int
SafeFragment::SameNumbers(const SafeFragment& f2, IWString& new_smiles) const {
  new_smiles = f2._smiles;

  // Now copy our ring numbers to `new_smiles`.
  const int nf = _first_digit.number_elements();
  for (int i = 0; i < nf; ++i) {
    new_smiles[f2._first_digit[i]] = _smiles[_first_digit[i]];
    new_smiles[f2._first_digit[i] + 1] = _smiles[_first_digit[i] + 1];
  }

  return 1;
}

}  // namespace safe_generate
