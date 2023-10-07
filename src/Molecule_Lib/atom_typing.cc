#include <stdlib.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <optional>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;

static int assign_arbitrary_values_to_unclassified_atoms = 0;

void
set_assign_arbitrary_values_to_unclassified_atoms(int s) {
  assign_arbitrary_values_to_unclassified_atoms = s;
}

static int use_version_2_augmented_atom_algorithm = 0;

void
set_use_version_2_augmented_atom_algorithm(int s) {
  use_version_2_augmented_atom_algorithm = s;
}

Atom_Typing_Specification::Atom_Typing_Specification() {
  _type = 0;
  _user_specified_type = 0;
  _differentiate_rings = 0;
  _perform_shell_iteration = 0;

  _assign_other_type = 0;
  _combine_donor_and_acceptor = 0;

  return;
}

int
Atom_Typing_Specification::active() const {
  if (0 != _type) {
    return 1;
  }

  if (0 != _user_specified_type) {
    return 1;
  }

  if (_query_and_value.size() > 0) {
    return 1;
  }

  return 0;
}

void
Atom_Typing_Specification::set_user_specified_type(unsigned int s) {
  _type = 0;

  _user_specified_type = s;

  return;
}

/*int
Atom_Typing_Specification::determine_atom_type_and_set_tag (const const_IWSubstring & s,
                                IWString & tag,
                                int verbose)
{
}*/

int
determine_atom_type_and_set_tag(const const_IWSubstring& s, IWString& tag, int verbose) {
  if (0 == tag.length()) {
    cerr << "determine_atom_type_and_set_tag:tag must be initialised\n";
    return 0;
  }

  if (!tag.starts_with("NC")) {
    cerr << "determine_atom_type_and_set_tag:warning, tag not NC type '" << tag << "'\n";
  }

  const_IWSubstring mys(s);

  int rc = 0;

  int differentiate_rings = 0;
  int perform_shell_expansion = 0;

  IWString p(s);

  if (p.starts_with('R'))  // note that 'R2' works but '2R' does not...
  {
    differentiate_rings = 1;
    p.remove_leading_chars(1);
  }

  if (isdigit(p.last_item())) {
    perform_shell_expansion = (p.last_item() - '0');
    p.pop();
  }

  IWString psave(p);  // might contain file names which must be preserved

  p.to_lowercase();
  if ('c' == p || "complex" == p) {
    rc = IWATTYPE_COMPLEX;
    tag << "C<";
    if (verbose) {
      cerr << "Will produce complex atom type, tag '" << tag << "'\n";
    }
  } else if ('z' == p) {
    rc = IWATTYPE_Z;
    tag << "Z<";
    if (verbose) {
      cerr << "Will produce atomic number atom type, tag '" << tag << "'\n";
    }
  } else if ("zp" == p) {
    rc = IWATTYPE_ZP;
    tag << "Z<";
    if (verbose) {
      cerr << "Will produce atomic number (prime) atom type, tag '" << tag << "'\n";
    }
  } else if ("tt" == p) {
    rc = IWATTYPE_TT;
    tag << "TT<";
    if (verbose) {
      cerr << "Will produce topotorsion atom type, tag '" << tag << "'\n";
    }
  } else if ("bs" == p) {
    rc = IWATTYPE_BASIC;
    tag << "BS<";
    if (verbose) {
      cerr << "Will produce basic atom type, tag '" << tag << "'\n";
    }
  } else if ("sb" == p) {
    rc = IWATTYPE_SYBYL;
    tag << "SB<";
    if (verbose) {
      cerr << "Will produce sibyl atom type, tag '" << tag << "'\n";
    }
  } else if (p.starts_with("pp")) {
    rc = IWATTYPE_PPHORE;
    tag << "PP<";
    if (verbose) {
      cerr << "Will produce pharmacaphore atom type, tag '" << tag << "'\n";
    }
  } else if ("expt" == p) {
    rc = IWATTYPE_EXPT;
    tag << "EXPT<";
    if (verbose) {
      cerr << "Will produce experimental atom type, tag '" << tag << "'\n";
    }
  } else if ("HB" == p) {
    rc = IWATTYPE_HB;
    tag << "HB<";
    if (verbose) {
      cerr << "Will produce HB atom type, tag '" << tag << "'\n";
    }
  } else if ("sf" == p) {
    rc = IWATTYPE_SF;
    differentiate_rings = 0;  // this typing includes ring membership
  } else if ("sfx" == p) {
    rc = IWATTYPE_SFX;
    differentiate_rings = 0;  // this typing includes ring membership
  } else if ("ch" == p) {
    rc = IWATTYPE_CH;
    tag << "CH<";
    differentiate_rings = 0;
  } else if (p.starts_with("expt=")) {
    p.remove_leading_chars(5);
  } else if ("none" == p) {
    rc = IWATTYPE_NONE;
    tag << "NON<";
    if (verbose) {
      cerr << " Will produce NONE atom type, tag '" << tag << "'\n";
    }
  } else if ("cc" == p) {
    rc = IWATTYPE_CC;
    tag << "CC";
    if (verbose) {
      cerr << "Will produce CC atom type, tag '" << tag << "'\n";
    }
  } else if ("za" == p) {
    rc = IWATTYPE_ZA;
    tag << "ZA";
    if (verbose) {
      cerr << "Will produce ZA atom type, tag '" << tag << "'\n";
    }
  } else if ("nox" == p) {
    rc = IWATTYPE_NOX;
    tag << "NOX";
    if (verbose) {
      cerr << "Will produce NOX (Nitrogen Oxygen and other) type, tag '" << tag << "'\n";
    }
  } else {
    cerr << "Unrecognised atom typing '" << p << "'\n";
    return 0;
  }

  if (differentiate_rings) {
    rc |= DIFFERENTIATE_RINGS;
  }

  if (perform_shell_expansion) {
    rc |= PERFORM_SHELL_ITERATION;
  }

  return rc;
}

int
Atom_Typing_Specification::_parse_user_specified_type(const const_IWSubstring& p) {
  IWString myp(p);
  myp.to_lowercase();

  _user_specified_type = 0;

  for (int i = 0; i < p.length(); i++) {
    char c = tolower(myp[i]);

    if ('z' == c) {
      _user_specified_type |= IWATTYPE_USP_Z;
    } else if ('y' == c) {
      _user_specified_type |= IWATTYPE_USP_Y;
    } else if ('h' == c) {
      _user_specified_type |= IWATTYPE_USP_H;
    } else if ('o' == c) {
      _user_specified_type |= IWATTYPE_USP_O;
    } else if ('p' == c) {
      _user_specified_type |= IWATTYPE_USP_P;
    } else if ('q' == c) {
      _user_specified_type |= IWATTYPE_USP_Q;
    } else if ('a' == c) {
      _user_specified_type |= IWATTYPE_USP_A;
    } else if ('c' == c) {
      _user_specified_type |= IWATTYPE_USP_C;
    } else if ('r' == c) {
      _user_specified_type |= IWATTYPE_USP_R;
    } else if ('e' == c) {
      _user_specified_type |= IWATTYPE_USP_E;
    } else if ('n' == c) {
      _user_specified_type |= IWATTYPE_USP_N;
    } else if ('u' == c) {
      _user_specified_type |= IWATTYPE_USP_U;
    } else if ('b' == c) {
      _user_specified_type |= IWATTYPE_USP_B;
    } else if ('i' == c) {
      _user_specified_type |= IWATTYPE_USP_I;
    } else if ('f' == c) {
      _user_specified_type |= IWATTYPE_USP_F;
    } else if ('m' == c) {
      _user_specified_type |= IWATTYPE_USP_M;
    } else if ('t' == c) {
      _user_specified_type |= IWATTYPE_USP_T;
    } else if ('x' == c) {
      _user_specified_type |= IWATTYPE_USP_X;
    } else if ('s' == c) {
      _user_specified_type |= IWATTYPE_USP_S;
    } else if ('l' == c) {
      _user_specified_type |= IWATTYPE_USP_L;
    } else if ('k' == c) {
      _user_specified_type |= IWATTYPE_USP_K;
    } else if ('g' == c) {
      _user_specified_type |= IWATTYPE_USP_G;
    } else {
      cerr << "Atom_Typing_Specification::_parse_user_specified_type:unrecognised '" << c
           << "'\n";
      return 0;
    }
  }

  return 1;
}

int
iwattype_convert_to_string_form(int atype, IWString& s) {
  if (IWATTYPE_Z == atype) {
    s = 'z';
  } else if (IWATTYPE_COMPLEX == atype) {
    s = 'c';
  } else if (IWATTYPE_TT == atype) {
    s = "tt";
  } else if (IWATTYPE_SYBYL == atype) {
    s = "sb";
  } else if (IWATTYPE_PPHORE == atype) {
    s = "pp";
  } else if (IWATTYPE_BASIC == atype) {
    s = "bs";
  } else if (IWATTYPE_SF == atype) {
    s = "sf";
  } else if (IWATTYPE_SFX == atype) {
    s = "sfx";
  } else if (IWATTYPE_CH == atype) {
    s = "ch";
  } else if (IWATTYPE_CC == atype) {
    s = "cc";
  } else if (IWATTYPE_NONE == atype) {
    s = "none";
  } else if (IWATTYPE_ZA == atype) {
    s = "za";
  } else if (IWATTYPE_ZP == atype) {
    s = "zp";
  } else if (IWATTYPE_NOX == atype) {
    s = "nox";
  } else {
    cerr << "iwattpe_convert_to_string_form:unrecognised form " << atype << '\n';
    return 0;
  }

  return 1;
}

int
determine_atom_type(const IWString& s) {
  IWString mys(s);
  mys.to_lowercase();

  // cerr << "Checking atom type '" << mys << "'\n";

  if ('z' == mys) {
    return IWATTYPE_Z;
  } else if ('y' == mys) {
    return IWATTYPE_USP_Y;
  } else if ('c' == mys || "complex" == s) {
    return IWATTYPE_COMPLEX;
  } else if ("tt" == mys) {
    return IWATTYPE_TT;
  } else if ("sb" == mys) {
    return IWATTYPE_SYBYL;
  } else if ("pp" == mys) {
    return IWATTYPE_PPHORE;
  } else if ("bs" == mys) {
    return IWATTYPE_BASIC;
  } else if ("sf" == mys) {
    return IWATTYPE_SF;
  } else if ("sfx" == mys) {
    return IWATTYPE_SFX;
  } else if ("ch" == mys) {
    return IWATTYPE_CH;
  } else if ("cc" == mys) {
    return IWATTYPE_CC;
  } else if ("none" == mys) {
    return IWATTYPE_NONE;
  } else if ("za" == mys) {
    return IWATTYPE_ZA;
  } else if ("zp" == mys) {
    return IWATTYPE_ZP;
  } else if ("nox" == mys) {
    return IWATTYPE_NOX;
  } else {
    cerr << "determine_atom_type:unrecognised type '" << s << "'\n";
    return 0;
  }
}

template <typename T>
int
assign_atom_types_za(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    invariant[i] = m.atomic_number(i);

    if (m.is_aromatic(i)) {
      invariant[i] = 2 * invariant[i] + 1;
    } else {
      invariant[i] = 2 * invariant[i];
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_none(Molecule& m, T* invariant) {
  std::fill_n(invariant, m.natoms(), 1);

  return 1;
}

template <typename T>
int
assign_atom_types_ch(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (6 == m.atomic_number(i)) {
      invariant[i] = 1;
    } else {
      invariant[i] = 11;
    }
  }

  return 1;
}

/*
  Jibo's initial implementation of atomic number atom
  typing
*/

template <typename T>
int
assign_atom_types_z_prime_numbers(Molecule& m, T* atom_constant) {
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    switch (m.atomic_number(i)) {
      case 6:
        atom_constant[i] = 13;
        break;
      case 7:
        atom_constant[i] = 17;
        break;
      case 8:
        atom_constant[i] = 19;
        break;
      case 9:
        atom_constant[i] = 23;
        break;
      case 15:
        atom_constant[i] = 29;
        break;
      case 16:
        atom_constant[i] = 31;
        break;
      case 17:
        atom_constant[i] = 37;
        break;
      case 35:
        atom_constant[i] = 41;
        break;
      case 53:
        atom_constant[i] = 43;
        break;
      default:
        atom_constant[i] = 83 + m.atomic_number(i);
        break;
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_z(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    invariant[i] = m.atomic_number(i);
  }

  return 1;
}

template <typename T>
int
assign_atom_types_z_ring(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    invariant[i] = m.atomic_number(i);

    if (m.is_ring_atom(i)) {
      invariant[i] += HIGHEST_ATOMIC_NUMBER;
    }
  }

  return 1;
}

/*
  Two variants of the same function, one where the ncon array is
  available, and one where it is not. Just for efficiency
*/

template <typename T>
int
assign_atom_types_complex(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int hcount = m.hcount(i);

    int pi;
    (void)m.pi_electrons(i, pi);

    if (pi > 4) {
      pi = 4;
    }

    int arom;
    // Cannot have a breaking change, but I like this better - needs testing.
#ifdef NEW_INCOMPATIBLE_VERSION
    if (m.is_aromatic(i)) {
      arom = 8;
    } else {
      arom = 0;
    }
#else
    arom = m.is_aromatic(i);
#endif

    invariant[i] = 100 * hcount + 10 * pi + 4 * arom + m.ncon(i);
  }

  return 1;
}

template <typename T>
int
assign_atom_types_complex(Molecule& m, const int* ncon, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const int hcount = m.hcount(i);

    int pi;
    (void)m.pi_electrons(i, pi);

    if (pi > 4) {
      pi = 4;
    }

    const int arom = m.is_aromatic(i);

    invariant[i] = 100 * hcount + 10 * pi + 4 * arom + ncon[i];
  }

  return 1;
}

template <typename T>
int
assign_atom_types_complex_ring(Molecule& m, T* invariant) {
  assign_atom_types_complex(m, invariant);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_ring_atom(i)) {
      invariant[i] = invariant[i] + 400;
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_complex_ring(Molecule& m, const int* ncon, T* invariant) {
  assign_atom_types_complex(m, ncon, invariant);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_ring_atom(i)) {
      invariant[i] = invariant[i] + 400;
    }
  }

  return 1;
}

/*
  Capture all the important elements, the rest are all equivalent

  Note that I have an off by 1 error here, and because of that
  all the halogens get the same atom type. Seems that is beneficial
*/

// clang-format off
static int z2inv[HIGHEST_ATOMIC_NUMBER + 1] = {
//0, //    this line initially missing!! Still missing! Never touch this again
  0,
  0,
  0,
  0,
  0,
  1,    // C
  2,    // N
  3,    // O
  4,    // F
  0,
  0,
  0,
  0,
  0,
  5,    // P
  6,    // S
  7,    // Cl
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  8,     // Br
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  9,      // I
  0,
  0,
  0,
  0,
  0,
  0,
  0,   
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0};

// clang-format on

template <typename T>
int
assign_atom_types_tt(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* ai = m.atomi(i);

    int pi = 0;
    (void)m.pi_electrons(i, pi);

    if (pi > 4) {
      pi = 4;
    }

    //  cerr << "Atom " << i << " type " << m.atomic_number(i) << " invariant " <<
    //  z2inv[ai->atomic_number()] << '\n';

    invariant[i] = 50 * (ai->ncon() - 1) + 10 * pi + z2inv[ai->atomic_number()];
  }

  return 1;
}

template <typename T>
int
assign_atom_types_tt(Molecule& m, const int* ncon, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int pi = 0;
    (void)m.pi_electrons(i, pi);

    if (pi > 4) {
      pi = 4;
    }

    const Atom* ai = m.atomi(i);

    invariant[i] = 50 * (ncon[i] - 1) + 10 * pi + z2inv[ai->atomic_number()];
  }

  return 1;
}

template <typename T>
int
assign_atom_types_tt_ring(Molecule& m, T* invariant) {
  assign_atom_types_tt(m, invariant);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_ring_atom(i)) {
      invariant[i] = invariant[i] + 400;
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_tt_ring(Molecule& m, const int* ncon, T* invariant) {
  assign_atom_types_tt(m, ncon, invariant);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_ring_atom(i)) {
      invariant[i] = invariant[i] + 400;
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_sybyl(Molecule& m, T* invariant) {
  int rc = m.find_simplified_sybyl_atom_type(reinterpret_cast<atom_type_t*>(invariant));

  if (rc) {
    return rc;
  }

  if (assign_arbitrary_values_to_unclassified_atoms) {
    return 1;
  }

  return 0;
}

template <typename T>
int
assign_atom_types_expt(Molecule& m, T* invariant) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
//  invariant[i] = m.ncon(i);
//  invariant[i] = 1;
#ifdef AT_HETEROATOM_ONLY
    if (6 == m.atomic_number(i)) {
      invariant[i] = 1;
    } else {
      invariant[i] = 2;
    }
#endif

#ifdef AT_PI_ELECTRONS
    int pi;
    m.pi_electrons(i, pi);

    if (pi > 0) {
      pi = 1;
    }
    invariant[i] = pi + 1;
#endif

#ifdef AT_HCOUNT
    if (m.hcount(i)) {
      invariant[i] = 1;
    } else {
      invariant[i] = 2;
    }
#endif

    if (m.is_aromatic(i)) {
      invariant[i] = 1;
    } else {
      invariant[i] = 2;
    }

    //  invariant[i] = 100 * m.atomic_number(i) + 20 * m.ncon(i);
    //  invariant[i] = 100 * m.atomic_number(i) + 20 * m.is_aromatic(i);

    atomic_number_t z = m.atomic_number(i);

    if (z <= 17) {
      invariant[i] = 1000 * z;
    } else {
      invariant[i] = 20000;
    }

    int pi;
    m.pi_electrons(i, pi);
    if (pi > 4) {
      pi = 4;
    }

    invariant[i] += pi;
  }

  return 1;
}

template <typename T>
int
assign_atom_types_basic(Molecule& m, T* atype) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      atype[i] = 1;
    } else if (6 == m.atomic_number(i)) {
      if (m.ncon(i) == m.nbonds(i)) {
        atype[i] = 2;
      } else {
        atype[i] = 3;
      }
    } else if (m.ncon(i) == m.nbonds(i)) {
      if (m.hcount(i)) {
        atype[i] = 4;
      } else {
        atype[i] = 5;
      }
    } else if (m.hcount(i)) {
      atype[i] = 6;
    } else {
      atype[i] = 7;
    }
  }

  return 1;
}

#ifdef OLDER_VERSION_SFX

static void
count_ring_type(Molecule& m, atom_number_t zatom, int& aromatic_rings,
                int& aliphatic_rings) {
  aromatic_rings = 0;
  aliphatic_rings = 0;

  m.compute_aromaticity_if_needed();

  int nr = m.nrings();

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    if (!ri->contains(zatom)) {
      continue;
    }

    if (ri->is_aromatic()) {
      aromatic_rings++;
    } else {
      aliphatic_rings++;
    }
  }

  nr = m.non_sssr_rings();

  if (0 == nr) {
    return;
  }

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.non_sssr_ring(i);

    if (!ri->contains(zatom)) {
      continue;
    }

    aliphatic_rings++;  // there are no non-sssr aromatic rings
  }

  return;
}

/*
  For ring types we have 11 cases
    0 rings
    1 rings
       Ar       100
       Al       200
    2 rings
       Ar Ar      300
       Al Al      400
       Ar Al      500

    3 rings
       Ar Ar Ar   600
       Ar Ar Al   700
       Ar Al Al   800
       Al Al Al   900

    more than 3 rings   1000
*/

template <typename T>
int
assign_atom_types_synthetic_feasibility_expt(Molecule& m, T* atype) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    atomic_number_t z = a->atomic_number();

    int t = 200000 * z + 40000 * a->ncon();

    t += 8000 * m.hcount(i);

    int nr = m.nrings(i);
    if (0 == nr) {
      ;
    } else if (1 == nr) {
      if (m.is_aromatic(i)) {
        t += 100;
      } else {
        t += 200;
      }
    } else if (nr > 3) {
      t += 1000;
    } else {
      int aromatic_rings, aliphatic_rings;
      count_ring_type(m, i, aromatic_rings, aliphatic_rings);
      if (2 == aromatic_rings && 0 == aliphatic_rings) {
        t += 300;
      } else if (0 == aromatic_rings && 2 == aliphatic_rings) {
        t += 400;
      } else if (1 == aromatic_rings && 2 == aliphatic_rings) {
        t += 500;
      } else if (3 == aromatic_rings) {
        t += 600;
      } else if (2 == aromatic_rings && 1 == aliphatic_rings) {
        t += 700;
      } else if (1 == aromatic_rings && 2 == aliphatic_rings) {
        t += 800;
      } else if (0 == aromatic_rings && 3 == aliphatic_rings) {
        t += 900;
      }
    }

    t += static_cast<int>(a->formal_charge());

    atype[i] = t;
  }

  return 1;
}

/*
  We want to make a compressed form of atom types.
*/

static int atomic_number_to_ndx[HIGHEST_ATOMIC_NUMBER + 1] = {
    0,
    1,   // H  1
    2,   // He 2
    2,   // Li 3
    2,   // Be 4
    5,   // B  5
    6,   // C  6
    7,   // N  7
    8,   // O  8
    9,   // F  9
    2,   // Ne 10
    2,   // Na 11
    2,   // Mg 12
    2,   // Al 13
    2,   // Si 14
    15,  // P  15
    16,  // S  16
    17,  // Cl 17
    2,   // Ar 18
    2,   // K  19
    2,   // Ca 20
    2,   // Sc 21
    2,   // Ti 22
    2,   // V  23
    2,   // Cr 24
    2,   // Mn 25
    2,   // Fe 26
    2,   // Co 27
    2,   // Ni 28
    2,   // Cu 29
    2,   // Zn 30
    2,   // Ga 31
    2,   // Ge 32
    2,   // As 33
    2,   // Se 34
    35,  // Br 35
    2,   // Kr 36
    2,   // Rb 37
    2,   // Sr 38
    2,   // Y  39
    2,   // Zr 40
    2,   // Nb 41
    2,   // Mo 42
    2,   // Tc 43
    2,   // Ru 44
    2,   // Rh 45
    2,   // Pd 46
    2,   // Ag 47
    2,   // Cd 48
    2,   // In 49
    2,   // Sn 50
    2,   // Sb 51
    2,   // Te 52
    53,  // I  53
    2,   // Xe 54
    2,   // Cs 55
    2,   // Ba 56
    2,   // La 57
    2,   // Ce 58
    2,   // Pr 59
    2,   // Nd 60
    2,   // Kr 61
    2,   // Kr 62
    2,   // Kr 63
    2,   // Kr 64
    2,   // Kr 65
    2,   // Kr 66
    2,   // Kr 67
    2,   // Kr 68
    2,   // Kr 69
    2,   // Kr 70
    2,   // Kr 71
    2,   // Kr 72
    2,   // Kr 73
    2,   // Kr 74
    2,   // Kr 75
    2,   // Kr 76
    2,   // Kr 77
    2,   // Kr 78
    2,   // Kr 79
    2,   // Kr 80
    2,   // Kr 81
    2,   // Kr 82
    2,   // Kr 83
    2,   // Kr 84
    2,   // Kr 85
    2,   // Kr 86
    2,   // Kr 87
    2,   // Kr 88
    2,   // Kr 89
    2,   // Kr 90
    2,   // Kr 91
    2,   // Kr 92
    2,   // Kr 93
    2,   // Kr 94
    2,   // Kr 95
    2,   // Kr 96
    2,   // Kr 97
    2,   // Kr 98
    2,   // Kr 99
    2,   // Kr 100
    2,   // Kr 101
    2,   // Kr 102
    2,   // Kr 103
    2,   // Kr 104
    2,   // Kr 105
    2,   // Kr 106
    2,   // Kr 107
    2,   // Kr 108
    2,   // Kr 109
    2,   // Kr 110
    2,   // Kr 111
    2,   // Kr 112
    2,   // Kr 113
    2,   // Kr 114
    2,   // Kr 115
    2    // Kr 116
};
#endif

/*
  For each atom set the corresponding item in RING_STATUS
    0 not a ring
    1 in an isolated aliphatic ring
    2 in an isolated aromatic ring
    3 in an aromatic ring
    4 in an aliphatic ring
    5 in both an aliphatic and an aromatic ring
*/

static void
determine_ring_status(Molecule& m, int* ring_status) {
  const int nr = m.nrings();

  if (0 == nr) {
    return;
  }

  m.compute_aromaticity_if_needed();

  if (1 == nr) {
    const Ring* r = m.ringi(0);
    if (r->is_aromatic()) {
      r->set_vector(ring_status, 2);
    } else {
      r->set_vector(ring_status, 1);
    }

    return;
  }

  // Now the case of more than one ring

  int* ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (ri->is_fused()) {
      continue;
    }

    if (ri->is_aromatic()) {
      ri->set_vector(ring_status, 2);
    } else {
      ri->set_vector(ring_status, 1);
    }

    ring_already_done[i] = 1;
  }

  // Set all atoms in fused aromatic rings to 3

  for (int i = 0; i < nr; ++i) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    if (!ri->is_aromatic()) {
      continue;
    }

    ri->set_vector(ring_status, 3);
    ring_already_done[i] = 1;
  }

  // Now look at the remaining aliphatic rings

  for (int i = 0; i < nr; ++i) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    if (ri->is_aromatic()) {  // should not happen
      continue;
    }

    const int ring_size = ri->number_elements();

    for (int j = 0; j < ring_size; ++j) {
      const auto k = ri->item(j);

      if (3 == ring_status[k]) {  // has been identified as aromatic, so it is in both
                                  // aliphatic and aromatic rings
        ring_status[k] = 5;
      } else {
        ring_status[k] = 4;
      }
    }
  }

  return;
}

/*
  Here we fix the problem with the halogen indexing above
*/

// clang-format off
static int z2inv_sfx[HIGHEST_ATOMIC_NUMBER + 1] = {
  0,
  0,
  0,
  0,
  0,
 10,    // B
  1,    // C
  2,    // N
  3,    // O
  4,    // F
  0,
  0,
  0,
  0,
 11,    // Si
  5,    // P
  6,    // S
  7,    // Cl
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  8,     // Br
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  9,      // I
  0,
  0,
  0,
  0,
  0,
  0,
  0,   
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0};

// clang-format on

template <typename T>
int
assign_atom_types_synthetic_feasibility_expt(Molecule& m, T* atype) {
  const int matoms = m.natoms();

  int* ring_status = new_int(matoms);
  std::unique_ptr<int[]> free_ring_status(ring_status);

  determine_ring_status(m, ring_status);

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    atomic_number_t z = a->atomic_number();

    int fc = a->formal_charge() + 1;

    atype[i] = (10 * 5 * 4 * 3) * z2inv_sfx[z] + (5 * 4 * 3) * ring_status[i] +
               (4 * 3) * a->ncon() + (3) * m.hcount(i) + fc;
    //  cerr << "Atom " << i << " type " << m.smarts_equivalent_for_atom(i) << " type " <<
    //  atype[i] << " z2inv " << z2inv_sfx[z] << '\n';
  }

  return 1;
}

template <typename T>
int
assign_atom_types_synthetic_feasibility(Molecule& m, T* atype) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    atomic_number_t z = a->atomic_number();

    if (m.is_aromatic(i)) {
      atype[i] = 5 * z;
      continue;
    }

    int r = m.is_ring_atom(i);

    if (a->nbonds() > a->ncon()) {
      atype[i] = 5 * z + 1 + r;
    } else {
      atype[i] = 5 * z + 3 + r;
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_basic(Molecule& m, const int* ncon, T* atype) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      atype[i] = 1;
    } else if (6 == m.atomic_number(i)) {
      if (ncon[i] == m.nbonds(i)) {
        atype[i] = 2;
      } else {
        atype[i] = 3;
      }
    } else if (ncon[i] == m.nbonds(i)) {
      if (m.hcount(i)) {
        atype[i] = 4;
      } else {
        atype[i] = 5;
      }
    } else if (m.hcount(i)) {
      atype[i] = 6;
    } else {
      atype[i] = 7;
    }
  }

  return 1;
}

template <typename T>
int
assign_atom_types_cc(Molecule& m, T* atype) {
  const int matoms = m.natoms();

  if (!assign_atom_types_complex(m, atype)) {
    return 0;
  }

  if (!m.has_formal_charges()) {
    return 1;
  }

  for (int i = 0; i < matoms; i++) {
    formal_charge_t fc = m.formal_charge(i);

    if (0 == fc) {
      ;
    } else if (fc > 0) {
      atype[i] = 703;
    } else {
      atype[i] = 809;
    }
  }

  return 2;
}

/*
  Differentiate Oxygen and Nitrogen from all others
*/

template <typename T>
int
assign_atom_types_nox(const Molecule& m, T* atype) {
  const int matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    const atomic_number_t z = m.atomic_number(i);

    if (6 == z) {
      atype[i] = 1;
    } else if (7 == z || 8 == z) {
      atype[i] = 7;
    } else {
      atype[i] = 1;
    }
  }

  return 1;
}

/*
  Need a number that does not collide with any of the donor/acceptor isotopes
*/

#define IWATTYPE_HYDROPHOBE_ISOTOPE 5

#define IWATTYPE_PPORE_POS 1
#define IWATTYPE_PPORE_NEG 2
#define IWATTYPE_PPORE_DON 4
#define IWATTYPE_PPORE_ACC 8
#define IWATTYPE_PPORE_HYP 16
#define IWATTYPE_PPORE_OTH 32

template <typename T>
int
Atom_Typing_Specification::_assign_atom_types_pharmacaphore(Molecule& m, T* atype) {
  const auto matoms = m.natoms();

  // cerr << "Atom_Typing_Specification::_assign_atom_types_pharmacaphore, ust " <<
  // _user_specified_type << '\n';

  if (0 != _user_specified_type) {
    (void)_assign_user_specified_type(m, atype, nullptr);
  }

  // cerr << "Assigned " << std::count_if(atype, atype + matoms, [] (int t) { return 0 !=
  // t;}) << '\n';

  (void)_charge_assigner.process(m);

  int* da = new_int(matoms);
  std::unique_ptr<int[]> free_da(da);

  (void)_donor_acceptor_assigner.process(m, da);

  // cerr << "After donor acceptor " << m.smiles() << '\n';

  Molecule_to_Match target(&m);

  for (auto i = 0; i < _hydrophobe.number_elements(); ++i) {
    Substructure_Results sresults;

    auto nhits = _hydrophobe[i]->substructure_search(target, sresults);

    for (auto j = 0; j < nhits; ++j) {
      const Set_of_Atoms* e = sresults.embedding(j);

      atom_number_t k = e->first();

      if (0 != da[k]) {  // do not allow multi-valued for hydrophobes
        continue;
      }

      da[k] = IWATTYPE_HYDROPHOBE_ISOTOPE;
    }
  }

  for (auto i = 0; i < matoms; ++i) {
    const Atom* a = m.atomi(i);

    const auto dai = da[i];

    const auto fc = a->formal_charge();

    if (0 == dai && 0 == fc) {
      if (_assign_other_type) {
        atype[i] = IWATTYPE_PPORE_OTH;
      }
      continue;
    }

    atype[i] = 0;

    if (0 == fc) {
      ;
    } else if (fc > 0) {
      atype[i] |= IWATTYPE_PPORE_POS;
    } else {
      atype[i] |= IWATTYPE_PPORE_NEG;
    }

    if (0 == dai) {
      ;
    } else if (_combine_donor_and_acceptor && dai <= 3) {
      atype[i] |= IWATTYPE_PPORE_ACC;  // just chosen at random
    } else if (1 == dai) {
      atype[i] |= IWATTYPE_PPORE_ACC;
    } else if (3 == dai) {
      atype[i] |= IWATTYPE_PPORE_DON;
    } else if (2 == dai) {
      atype[i] |= (IWATTYPE_PPORE_DON | IWATTYPE_PPORE_ACC);
    } else if (IWATTYPE_HYDROPHOBE_ISOTOPE == dai) {
      atype[i] |= IWATTYPE_PPORE_HYP;
    }
  }

  if (_write_isotopically_labelled.active()) {
    Molecule mcopy(m);
    for (auto i = 0; i < matoms; ++i) {
      mcopy.set_isotope(i, atype[i]);
    }
    _write_isotopically_labelled.write(mcopy);
  }

  return 1;
}

template <typename T>
int
_assign_atom_types(Molecule& m, int typing_to_use, T* atype, const int* ncon) {
  int differentiate_rings = 0;
  if (DIFFERENTIATE_RINGS & typing_to_use) {
    differentiate_rings = 1;
    typing_to_use ^= DIFFERENTIATE_RINGS;
  }

  // cerr << "Assigning, DIFFERENTIATE_RINGS " << differentiate_rings << " type " <<
  // typing_to_use << '\n';

  const int matoms = m.natoms();

  if (0 == matoms) {
    return 0;
  }

  if (IWATTYPE_Z == typing_to_use) {
    if (differentiate_rings) {
      return assign_atom_types_z_ring(m, atype);
    } else {
      return assign_atom_types_z(m, atype);
    }
  }

  if (IWATTYPE_COMPLEX == typing_to_use) {
    if (nullptr == ncon) {
      if (differentiate_rings) {
        return assign_atom_types_complex_ring(m, atype);
      } else {
        return assign_atom_types_complex(m, atype);
      }
    } else {
      if (differentiate_rings) {
        return assign_atom_types_complex_ring(m, ncon, atype);
      } else {
        return assign_atom_types_complex(m, ncon, atype);
      }
    }
  }

  if (IWATTYPE_ZP == typing_to_use) {
    return assign_atom_types_z_prime_numbers(m, atype);
  }

  if (IWATTYPE_TT == typing_to_use) {
    if (nullptr == ncon) {
      if (differentiate_rings) {
        return assign_atom_types_tt_ring(m, atype);
      } else {
        return assign_atom_types_tt(m, atype);
      }
    } else {
      if (differentiate_rings) {
        return assign_atom_types_tt_ring(m, ncon, atype);
      } else {
        return assign_atom_types_tt(m, ncon, atype);
      }
    }
  }

#ifdef IWHB
  if (IWATTYPE_HB == typing_to_use) {
    if (nullptr == ncon) {
      return assign_atom_types_hb(m, atype);
    } else {
      return assign_atom_types_hb(m, ncon, atype);
    }
  }
#endif

  if (IWATTYPE_SYBYL == typing_to_use) {
    return assign_atom_types_sybyl(m, atype);
  }

  if (IWATTYPE_BASIC == typing_to_use) {
    assert(!differentiate_rings);  // not implemented

    if (nullptr == ncon) {
      return assign_atom_types_basic(m, atype);
    } else {
      return assign_atom_types_basic(m, ncon, atype);
    }
  }

  if (IWATTYPE_SF == typing_to_use) {
    return assign_atom_types_synthetic_feasibility(m, atype);
  }

  if (IWATTYPE_SFX == typing_to_use) {
    return assign_atom_types_synthetic_feasibility_expt(m, atype);
  }

  if (IWATTYPE_EXPT == typing_to_use) {
    return assign_atom_types_expt(m, atype);
  }

  if (IWATTYPE_NONE == typing_to_use) {
    return assign_atom_types_none(m, atype);
  }

  if (IWATTYPE_ZA == typing_to_use) {
    return assign_atom_types_za(m, atype);
  }

  if (IWATTYPE_CH == typing_to_use) {
    return assign_atom_types_ch(m, atype);
  }

  if (IWATTYPE_CC == typing_to_use) {
    return assign_atom_types_cc(m, atype);
  }

  if (IWATTYPE_NOX == typing_to_use) {
    return assign_atom_types_nox(m, atype);
  }

  cerr << "What to do for atom type " << typing_to_use << '\n';
  return 0;
}

static int
bond_constant(const Bond* b) {
  if (b->is_aromatic()) {
    return 1;
  }

  if (b->is_single_bond()) {
    return 2;
  }

  if (b->is_double_bond()) {
    return 3;
  }

  if (b->is_triple_bond()) {
    return 4;
  }

  return 5;  // perhaps a coordination bond
}

template <typename T>
int
assign_atom_types(Molecule& m, int typing_to_use, T* atype, const int* ncon) {
  if (!_assign_atom_types(m, typing_to_use, atype, ncon)) {
    return 0;
  }

  if (0 == (PERFORM_SHELL_ITERATION & typing_to_use)) {
    return 1;
  }

  typing_to_use = (PERFORM_SHELL_ITERATION ^ typing_to_use);

  const int matoms = m.natoms();

  if (matoms < 2) {  // no expansion possible
    return 1;
  }

  m.compute_aromaticity_if_needed();

  int* updated_value = new_int(matoms);
  std::unique_ptr<int[]> free_updated_value(updated_value);

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    int acon = a->ncon();

    updated_value[i] = 792 * atype[i];

    for (int j = 0; j < acon; j++) {
      const Bond* b = a->item(j);

      int bc = bond_constant(b);

      atom_number_t k = b->other(i);

      updated_value[i] += bc * atype[k];
    }
  }

  copy_vector(atype, updated_value, matoms);

  return 1;
}

static void
display_atom_type_specifications(std::ostream& os) {
  os << " -P z              atomic number atom type\n";
  os << " -P c              \"c\" atom type\n";
  os << " -P tt             topological torsion atom type\n";
  os << " -P syb            Sibyl atom type\n";
  os << " -P pp             Pharmacaphore type\n";
  os << " -P UST:achprtyzne User Specified Type\n";
  os << "        A          aromatic or not\n";
  os << "        C          number of connections\n";
  os << "        E          carbon atoms get one type, all heteroatoms a different "
        "type\n";
  os << "        F          ring fusion\n";
  os << "        G          atomic symbol hash\n";
  os << "        H          implicit hydrogen count\n";
  os << "        I          isotope\n";
  os << "        K          atomic numbers of connected atoms\n";
  os << "        L          largest ring\n";
  os << "        M          all aromatic atoms the same\n";
  os << "        N          all atoms get the same type\n";
  os << "        O          formal charge\n";
  os << "        P          pi electron count\n";
  os << "        Q          presence or absence of a pi electron\n";
  os << "        R          ring bond count\n";
  os << "        S          smallest ring\n";
  os << "        T          atomic number, but possibly tautomeric nitrogens grouped\n";
  os << "        U          unsaturation\n";
  os << "        X          centrality of atom (expensive, uses distance matrix)\n";
  os << "        Y          atomic number, heavy halogens compressed\n";
  os << "        Z          atomic numbers\n";
  os << "              end specification with number for shell expansion\n";

  exit(1);
}

/*
  Syntax is
  type fractal

  C0         Means C atom type, zero fractal iterations
  Z1         Means atomic number type, one fractal iteration
*/

int
Atom_Typing_Specification::build(const const_IWSubstring& s) {
  _built_from = s;

  const_IWSubstring mys(s);

  if (mys.starts_with("R:")) {
    _differentiate_rings = 1;

    mys += 2;
  }

  if (0 == mys.length()) {
    cerr << "Atom_Typing_Specification::build:empty specification\n";
    return 0;
  }

  if (isdigit(mys[0])) {
    cerr << "Atom_Typing_Specification::build:cannot start with digit '" << s << "'\n";
    return 0;
  }

  // If it ends with a digit, extract the fractal expansion value. Note
  // we don't allow for multi-digit fractal atom types

  if (isdigit(mys.last_item())) {
    _perform_shell_iteration = mys.last_item() - '0';
    mys.chop();
    //  cerr << "Detected shell iteration " << _perform_shell_iteration << '\n';
  }

  if (mys.starts_with("EXT:")) {
    mys.remove_leading_chars(4);
    IWString tmp(mys);
    if (! ParseExternalQuery(tmp)) {
      return 0;
    }

    return 1;
  }

  if (mys.starts_with("UST:") || mys.starts_with("ust:")) {
    mys.remove_leading_chars(4);
    if (!_parse_user_specified_type(mys)) {
      return 0;
    }

    return 1;
  }

  if ("help" == mys) {
    display_atom_type_specifications(cerr);
  }

  // Some types might have extra information. For example the PP type might have
  // a qualifier indicating the file to read.

  const_IWSubstring stype, qualifier;

  if (mys.split(stype, '=', qualifier) && stype.length() > 0 && qualifier.length() > 0) {
    _type = determine_atom_type(stype);
  } else {
    _type = determine_atom_type(mys);
  }

  if (0 == _type) {
    return 0;
  }

  if (IWATTYPE_PPHORE == _type) {
    return _build_pharmacaphore_specification(mys);
  }

  // cerr << "TYpe " << _type << '\n';

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::assign_atom_types(Molecule& m, T* atype, const int* ncon) {
  std::fill_n(atype, m.natoms(), 0);

  if (IWATTYPE_PPHORE ==
      _type)  // because this owns the charge assigner etc queries. Must be done first
  {
    if (!_assign_atom_types_pharmacaphore(m, atype)) {  // we allow
      return 0;
    }
  } else if (0 != _user_specified_type) {
    if (!_assign_user_specified_type(m, atype, ncon)) {
      return 0;
    }
  } else if (_query_and_value.size() > 0) {
    if (! AssignExternalAtomTypes(m, atype)) {
      return 0;
    }
  } else if (!_assign_atom_types(m, _type, atype, ncon)) {
    return 0;
  }

  if (0 == _perform_shell_iteration) {
    return 1;
  }

  if (use_version_2_augmented_atom_algorithm) {
    _perform_shell_expansion_v2(m, atype);
  } else {
    _perform_shell_expansion_v1(m, atype);
  }

  return 1;
}

/*
  This really does not work properly for anything but radius 1
*/

template <typename T>
void
Atom_Typing_Specification::_perform_shell_expansion_v1(Molecule& m, T* atype) const {
  int need_to_compute_bond_constants = 1;

  if (IWATTYPE_NONE == _type) {
    need_to_compute_bond_constants = 0;
  } else if (IWATTYPE_CH == _type) {
    need_to_compute_bond_constants = 0;
  }

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  T* tmp = new T[matoms];
  std::unique_ptr<T[]> free_tmp(tmp);

  for (int i = 0; i < _perform_shell_iteration; i++) {
    for (int j = 0; j < matoms; j++) {
      tmp[j] = 792 * atype[j];

      const Atom* a = m.atomi(j);

      int acon = a->ncon();

      for (int k = 0; k < acon; k++) {
        const Bond* b = a->item(k);

        int bc;
        if (need_to_compute_bond_constants) {
          bc = bond_constant(b);
        } else {
          bc = 1;
        }

        atom_number_t l = b->other(j);

        tmp[j] += bc * atype[l];
      }
    }

    copy_vector(atype, tmp, matoms);
  }

  return;
}

template <typename T>
void
Atom_Typing_Specification::_perform_shell_expansion_v2(Molecule& m, T* atype) const {
  int need_to_compute_bond_constants = 1;

  if (IWATTYPE_NONE == _type) {
    need_to_compute_bond_constants = 0;
  } else if (IWATTYPE_CH == _type) {
    need_to_compute_bond_constants = 0;
  }

  m.compute_aromaticity_if_needed();
  const auto matoms = m.natoms();

  assert(sizeof(T) == sizeof(int));

  int* complete = new int[matoms + matoms];
  std::unique_ptr<int[]> free_complete(
      complete);  // hopefully one larger allocation faster than two smaller
  T* newtype = reinterpret_cast<T*>(complete + matoms);

  for (auto i = 0; i < _perform_shell_iteration; ++i) {
    for (auto j = 0; j < matoms; j++) {
      std::fill_n(complete, matoms, 0);

      newtype[j] = _perform_shell_iteration_v2(m, j, atype, complete,
                                               need_to_compute_bond_constants,
                                               _perform_shell_iteration);
    }
  }

  // std::copy_n(newtype, matoms, atype);
  copy_vector(atype, newtype, matoms);

  return;
}

/*
  Aug 2013.
  Radius 2 shell expansion never worked anywhere. Perhaps it was
  because the typing just got too messed up above, due to the poor
  algorithm used. This version gets the augmented atom type more
  carefully - atoms are not double counted
*/

template <typename T>
int
Atom_Typing_Specification::_perform_shell_iteration_v2(Molecule& m,
                                                       const atom_number_t zatom,
                                                       const T* atype, int* complete,
                                                       int need_to_compute_bond_constants,
                                                       int radius) const {
  complete[zatom] = 1;

  int rc = 792 * atype[zatom] + radius -
           _perform_shell_iteration;  // we want the first iteration of expansion to
                                      // reflect the proper atom type

  const Atom* a = m.atomi(zatom);

  auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i) {
    const Bond* b = a->item(i);
    auto j = b->other(zatom);

    if (complete[j]) {
      continue;
    }

    int bc;
    if (need_to_compute_bond_constants) {
      bc = bond_constant(b);
    } else {
      bc = 1;
    }

    rc += bc * atype[j];

    if (radius > 0) {
      rc = rc + _perform_shell_iteration_v2(m, j, atype, complete,
                                            need_to_compute_bond_constants, radius - 1);
    }
  }

  return rc;
}

static constexpr int kCarbon = 3001;
static constexpr int kNitrogen = 6007;
static constexpr int kOxygen = 9001;

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_z_prime_numbers(
    const Molecule& m, T* atype, int compress_halogens) const {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    switch (m.atomic_number(i)) {
      case 6:
        atype[i] = kCarbon;
        break;
      case 7:
        atype[i] = kNitrogen;
        break;
      case 8:
        atype[i] = kOxygen;
        break;
      case 9:
        atype[i] = 12007;
        break;
      case 14:
        atype[i] = 15013;
        break;
      case 15:
        atype[i] = 18013;
        break;
      case 16:
        atype[i] = 21001;
        break;
      case 17:
        atype[i] = 24001;
        break;
      case 35:
        if (compress_halogens) {
          atype[i] = 24001;
        } else {
          atype[i] = 27011;
        }
        break;
      case 53:
        if (compress_halogens) {
          atype[i] = 24001;
        } else {
          atype[i] = 30011;
        }
        break;
      default:
        atype[i] = 33013 + m.atomic_number(i);  // Maybe should be optional
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_hcount(Molecule& m, T* atype) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    const int h = m.hcount(i);

    if (0 == h) {
      ;
    } else if (1 == h) {
      atype[i] += 743;
    } else if (2 == h) {
      atype[i] += 1511;
    } else if (3 == h) {
      atype[i] += 1801;
    } else {
      atype[i] += 2111;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_pi(Molecule& m, T* atype) const {
  int pi;
  for (int i = m.natoms() - 1; i >= 0; i--) {
    m.pi_electrons(i, pi);

    if (0 == pi) {
      ;
    } else if (1 == pi) {
      atype[i] += 167;
    } else if (2 == pi) {
      atype[i] += 317;
    } else if (3 == pi) {
      atype[i] += 641;
    } else {
      atype[i] += 691;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_pi_boolean(Molecule& m,
                                                             T* atype) const {
  int pi;
  for (int i = m.natoms() - 1; i >= 0; i--) {
    m.pi_electrons(i, pi);

    if (pi > 0) {
      atype[i] += 167;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_isotope(const Molecule& m,
                                                          T* atype) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    atype[i] += m.isotope(i);
  }

  return 1;
}

/*
  Note that we truncate ring system size considerations at 1. Hard to know what will be
  most useful...
*/

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_ring_fusion(Molecule& m,
                                                              T* atype) const {
  if (0 == m.nrings()) {
    return 1;
  }

  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (!m.is_ring_atom(i)) {
      continue;
    }

    const int f = m.fused_system_size(i);

    if (0 == f) {  // not sure this can happen
      ;
    } else if (1 == f) {
      atype[i] = atype[i] * 3 + 311;
    } else if (f > 1) {
      atype[i] = atype[i] * 4 + 17;
    }
  }

  return 1;
}

// just some random number

#define TYPE_M_COMMON_AROMATIC_TYPE 12737

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_aromatic_all_the_same(Molecule& m,
                                                                        T* atype) const {
  // cerr << "_ust_assign_atom_types_aromatic_all_the_same\n";
  _ust_assign_atom_types_z_prime_numbers(m, atype, 1);  // regular Y type

  m.compute_aromaticity_if_needed();

  if (0 == m.nrings()) {
    return 1;
  }

  if (0 == m.aromatic_ring_count()) {
    return 1;
  }

  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    if (m.is_aromatic(i)) {
      atype[i] = TYPE_M_COMMON_AROMATIC_TYPE;
    }
  }

  return 1;
}

/*
  We increment the existing atom type
  Don't like this because it will distinguish quite similar
  atoms when they should not be differentiated.
*/

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_x(Molecule& m, T* atype) const {
  m.recompute_distance_matrix();

  const int matoms = m.natoms();

  if (1 == matoms) {
    atype[0] += 5;
    return 1;
  }

  for (int i = 0; i < matoms; ++i) {
    Accumulator_Int<int> acc;

    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }

      acc.extra(m.bonds_between(i, j));
    }

    atype[i] += 10 * static_cast<int>(sqrt(static_cast<float>(acc.maxval()))) +
                static_cast<int>(acc.average());
  }

  return 1;
}

template <typename T>
int
ust_assign_atom_types_t_ring(Molecule& m, const Ring& r, T* atype) {
  const int ring_size = r.number_elements();

  Set_of_Atoms two_connected_nitrogen_count, exocyclic_oxygen;

  for (auto i = 0; i < ring_size; ++i) {
    const auto j = r[i];

    const auto a = m.atomi(j);

    const auto z = a->atomic_number();

    if (7 == z && 2 == a->ncon()) {
      two_connected_nitrogen_count.add(j);
      continue;
    }

    if (6 != z) {
      continue;
    }

    if (3 != a->ncon()) {
      continue;
    }

    for (auto k = 0; k < 3; ++k) {
      const Bond* b = a->item(k);

      if (b->nrings()) {
        continue;
      }

      const auto o = b->other(j);

      if (1 == m.ncon(i) && 8 == m.atomic_number(o)) {  // =O or -O is fine
        exocyclic_oxygen.add(o);
      }

      break;
    }
  }

  if (two_connected_nitrogen_count.empty()) {
    return 0;
  }

  if (two_connected_nitrogen_count.number_elements() > 1) {
    ;
  } else if (exocyclic_oxygen.number_elements() > 0) {
    ;
  } else {
    return 0;
  }

  two_connected_nitrogen_count.set_vector(atype, static_cast<T>(37013));
  if (exocyclic_oxygen.number_elements()) {
    exocyclic_oxygen.set_vector(atype, static_cast<T>(1282));
  }

  return 1;
}

// Return true if `zatom` is the sulphur of a sulfonamide.
static int
IsSulfonamide(const Molecule& m,
              atom_number_t zatom) {
  const Atom& sulphur = m.atom(zatom);
  if (sulphur.ncon() != 4) {
    return 0;
  }

  int doubly_bonded_oxygen = 0;
  atom_number_t nitrogen = INVALID_ATOM_NUMBER;
  for (const Bond* b : sulphur) {
    const atom_number_t nbr = b->other(zatom);

    if (b->is_double_bond() && m.atomic_number(nbr) == 8) {
      ++doubly_bonded_oxygen;
    } else if (b->is_single_bond() && m.atomic_number(nbr) == 7) {
      nitrogen = nbr;
    }
  }

  if (doubly_bonded_oxygen < 2 ||
      nitrogen == INVALID_ATOM_NUMBER) {
    return 0;
  }

  return 1;
}

// Nove 2022. Modified to make low connected sulphur equivalent to oxygen.
// Do the same for sulfonamides, but that is only slighly effective because
// of the different connectivity of a sulfonamide (4) vs an amide (3).
template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_t(Molecule & m,
                                                    T * atype) const
{
  _ust_assign_atom_types_z_prime_numbers(m, atype, 1);    // last arg means compress halogens

  // Make 1 and 2 connected Sulphur same as Oxygen, and
  // Sulfonamides same as amides.
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.atomic_number(i) != 16) {
      continue;
    }
    const Atom& s = m.atom(i);
    if (s.ncon() <= 2) {
      atype[i] = kOxygen;
      continue;
    }
    if (IsSulfonamide(m, i)) {
      atype[i] = kCarbon;
    }
  }

  const int nr = m.nrings();

  if (0 == nr)
    return 1;

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    ust_assign_atom_types_t_ring(m, *ri, atype);
  }

  return 1;
}


#ifdef OLD_VERSION_PRIOR_NOV_22
template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_t(Molecule& m, T* atype) const {
  _ust_assign_atom_types_z_prime_numbers(m, atype,
                                         1);  // last arg means compress halogens

  const int nr = m.nrings();

  if (0 == nr) {
    return 1;
  }

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < nr; ++i) {
    const Ring* ri = m.ringi(i);

    if (!ri->is_aromatic()) {
      continue;
    }

    ust_assign_atom_types_t_ring(m, *ri, atype);
  }

  return 1;
}
#endif

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_arom(Molecule& m, T* atype) const {
  if (0 == m.nrings()) {
    return 1;
  }

  for (int i = m.natoms() - 1; i >= 0; i--) {
    if (m.is_aromatic(i)) {
      atype[i] += 37;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_ncon(const Molecule& m,
                                                       T* atype) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    switch (m.ncon(i)) {
      case 0:
        break;
      case 1:
        atype[i] += 7;
        break;
      case 2:
        atype[i] += 11;
        break;
      case 3:
        atype[i] += 13;
        break;
      case 4:
        atype[i] += 17;
        break;
    }

    atype[i] += m.ncon(
        i);  // probably just a bug this still being here. It is not hurting anything...
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_ncon(const Molecule& m, T* atype,
                                                       const int* ncon) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    switch (ncon[i]) {
      case 0:
        break;
      case 1:
        atype[i] += 7;
        break;
      case 2:
        atype[i] += 11;
        break;
      case 3:
        atype[i] += 13;
        break;
      case 4:
        atype[i] += 17;
        break;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_heteroatom(const Molecule& m,
                                                             T* atype) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    if (6 == m.atomic_number(i)) {
      atype[i] += 761;
    } else {
      atype[i] += 2477;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_unsaturated(const Molecule& m,
                                                              T* atype) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    const Atom* a = m.atomi(i);
    const int acon = a->ncon();
    const int nbonds = a->nbonds();

    atype[i] += 1025 + (nbonds - acon);
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_unsaturated_x_aromatic(Molecule & m,
                                                T * atype) const
{
  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    if (m.is_aromatic(i)) {
      continue;
    }

    const Atom * a = m.atomi(i);
    const int acon = a->ncon();
    const int nbonds = a->nbonds();

    atype[i] += 7878 + (9 * nbonds - acon);
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_none(const Molecule& m,
                                                       T* atype) const {
  for (int i = m.natoms() - 1; i >= 0; i--) {
    atype[i] = 763;
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_nrings(Molecule& m, T* atype) const {
  if (0 == m.nrings()) {
    return 1;
  }

  for (int i = m.natoms() - 1; i >= 0; i--) {
    atype[i] += 232 * m.ring_bond_count(i);
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_smallest_ring(Molecule& m,
                                                                T* atype) const {
  const int nr = m.nrings();

  if (0 == nr) {
    return 1;
  }

  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (!m.is_ring_atom(i)) {
      continue;
    }

    for (int j = 0; j < nr; ++j) {
      const Ring* rj = m.ringi(j);
      if (rj->contains(i)) {
        atype[i] += (73 * rj->number_elements());
        break;
      }
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_largest_ring(Molecule& m,
                                                               T* atype) const {
  const int nr = m.nrings();

  if (0 == nr) {
    return 1;
  }

  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (!m.is_ring_atom(i)) {
      continue;
    }

    for (int j = nr - 1; j >= 0; --j)  // scan from largest ring size downwards
    {
      const Ring* rj = m.ringi(j);
      if (rj->contains(i)) {
        atype[i] += (19 * rj->number_elements());
        break;
      }
    }
  }

  return 1;
}

static int
compute_connected_atoms_hash(resizable_array<int>& z) {
  const int n = z.number_elements();

  if (0 == n) {
    return 0;
  }

  if (1 == n) {
    return z[0];
  }

  if (2 == n) {
    if (z[0] > z[1]) {
      z.swap_elements(0, 1);
    }
  } else {
    Int_Comparator_Larger icl;
    z.iwqsort(icl);
  }

  int rc = z[0];

  for (int i = 1; i < n; ++i) {
    rc = 100 * rc + z[i];
  }

  return rc;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_connected_atoms_hash(const Molecule& m,
                                                                       T* atype) const {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const Atom* a = m.atomi(i);

    const int acon = a->ncon();

    resizable_array<int> z;

    for (int j = 0; j < acon; ++j) {
      const atom_number_t k = a->other(i, j);
      z.add(m.atomic_number(k));
    }

    atype[i] += compute_connected_atoms_hash(z);
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_atomic_symbol_hash_value(
    const Molecule& m, T* atype) const {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const Element* e = m.elementi(i);
    atype[i] = e->atomic_symbol_hash_value();
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_ust_assign_atom_types_formal_charge(const Molecule& m,
                                                                T* atype) const {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const formal_charge_t q = m.formal_charge(i);

    if (0 == q) {
      continue;
    } else if (q > 0) {
      atype[i] += q * 113;
    } else {
      atype[i] += q * 204;
    }
  }

  return 1;
}

template <typename T>
int
Atom_Typing_Specification::_assign_user_specified_type(Molecule& m, T* atype,
                                                       const int* ncon) const {
  const int matoms = m.natoms();

  std::fill_n(atype, matoms, 1);  // don't start with zero, hmmm, why not?

  unsigned int t = _user_specified_type;
  while (0 != t) {
    if (t & IWATTYPE_USP_Z)  // the two atomic number types set atype
    {
      _ust_assign_atom_types_z_prime_numbers(m, atype, 0);
      t = t ^ IWATTYPE_USP_Z;
    } else if (t & IWATTYPE_USP_Y)  // the two atomic number types set atype
    {
      _ust_assign_atom_types_z_prime_numbers(m, atype, 1);
      t = t ^ IWATTYPE_USP_Y;
    } else if (t & IWATTYPE_USP_G) {
      _ust_assign_atom_types_atomic_symbol_hash_value(m, atype);
      t = t ^ IWATTYPE_USP_G;
    } else if (t & IWATTYPE_USP_N) {
      _ust_assign_atom_types_none(m, atype);
      t = t ^ IWATTYPE_USP_N;
    } else if (t & IWATTYPE_USP_T) {
      _ust_assign_atom_types_t(m, atype);
      t = t ^ IWATTYPE_USP_T;
    } else if (t & IWATTYPE_USP_M) {
      _ust_assign_atom_types_aromatic_all_the_same(m, atype);
      t = t ^ IWATTYPE_USP_M;
    } else if (t & IWATTYPE_USP_X) {
      _ust_assign_atom_types_x(m, atype);
      t = t ^ IWATTYPE_USP_X;
    } else if (t & IWATTYPE_USP_H) {
      _ust_assign_atom_types_hcount(m, atype);
      t = t ^ IWATTYPE_USP_H;
    } else if (t & IWATTYPE_USP_O) {
      _ust_assign_atom_types_formal_charge(m, atype);
      t = t ^ IWATTYPE_USP_O;
    } else if (t & IWATTYPE_USP_P) {
      _ust_assign_atom_types_pi(m, atype);
      t = t ^ IWATTYPE_USP_P;
    } else if (t & IWATTYPE_USP_Q) {
      _ust_assign_atom_types_pi_boolean(m, atype);
      t = t ^ IWATTYPE_USP_Q;
    } else if (t & IWATTYPE_USP_I) {
      _ust_assign_atom_types_isotope(m, atype);
      t = t ^ IWATTYPE_USP_I;
    } else if (t & IWATTYPE_USP_F) {
      _ust_assign_atom_types_ring_fusion(m, atype);
      t = t ^ IWATTYPE_USP_F;
    } else if (t & IWATTYPE_USP_A) {
      _ust_assign_atom_types_arom(m, atype);
      t = t ^ IWATTYPE_USP_A;
    } else if (t & IWATTYPE_USP_C) {
      if (ncon) {
        _ust_assign_atom_types_ncon(m, atype, ncon);
      } else {
        _ust_assign_atom_types_ncon(m, atype);
      }
      t = t ^ IWATTYPE_USP_C;
    } else if (t & IWATTYPE_USP_R) {
      _ust_assign_atom_types_nrings(m, atype);
      t = t ^ IWATTYPE_USP_R;
    } else if (t & IWATTYPE_USP_E) {
      _ust_assign_atom_types_heteroatom(m, atype);
      t = t ^ IWATTYPE_USP_E;
    } else if (t & IWATTYPE_USP_U) {
      _ust_assign_atom_types_unsaturated(m, atype);
      t = t ^ IWATTYPE_USP_U;
    } else if (t & IWATTYPE_USP_B) {
      _ust_assign_atom_types_unsaturated_x_aromatic(m, atype);
      t = t ^ IWATTYPE_USP_B;
    } else if (t & IWATTYPE_USP_S) {
      _ust_assign_atom_types_smallest_ring(m, atype);
      t = t ^ IWATTYPE_USP_S;
    } else if (t & IWATTYPE_USP_L) {
      _ust_assign_atom_types_largest_ring(m, atype);
      t = t ^ IWATTYPE_USP_L;
    } else if (t & IWATTYPE_USP_K) {
      _ust_assign_atom_types_connected_atoms_hash(m, atype);
      t = t ^ IWATTYPE_USP_K;
    }

#ifdef DEBUG_UST_ATOM_TYPE
    cerr << "t = " << std::hex << t << std::dec;
    for (int i = 0; i < matoms; ++i) {
      cerr << " " << m.smarts_equivalent_for_atom(i) << ' ' << atype[i];
    }
    cerr << '\n';
#endif
  }

  if (_perform_shell_iteration) {
  }

  return 1;
}

int
Atom_Typing_Specification::_append_tag_for_user_specified_type(IWString& tag) const {
  if (IWATTYPE_USP_Z & _user_specified_type) {
    tag << 'Z';
  }
  if (IWATTYPE_USP_Y & _user_specified_type) {
    tag << 'Y';
  }
  if (IWATTYPE_USP_H & _user_specified_type) {
    tag << 'H';
  }
  if (IWATTYPE_USP_O & _user_specified_type) {
    tag << 'O';
  }
  if (IWATTYPE_USP_P & _user_specified_type) {
    tag << 'P';
  }
  if (IWATTYPE_USP_Q & _user_specified_type) {
    tag << 'Q';
  }
  if (IWATTYPE_USP_A & _user_specified_type) {
    tag << 'A';
  }
  if (IWATTYPE_USP_C & _user_specified_type) {
    tag << 'C';
  }
  if (IWATTYPE_USP_R & _user_specified_type) {
    tag << 'R';
  }
  if (IWATTYPE_USP_U & _user_specified_type) {
    tag << 'U';
  }
  if (IWATTYPE_USP_B & _user_specified_type) {
    tag << 'B';
  }
  if (IWATTYPE_USP_I & _user_specified_type) {
    tag << 'I';
  }
  if (IWATTYPE_USP_F & _user_specified_type) {
    tag << 'F';
  }
  if (IWATTYPE_USP_M & _user_specified_type) {
    tag << 'M';
  }
  if (IWATTYPE_USP_T & _user_specified_type) {
    tag << 'T';
  }
  if (IWATTYPE_USP_X & _user_specified_type) {
    tag << 'X';
  }
  if (IWATTYPE_USP_S & _user_specified_type) {
    tag << 'S';
  }
  if (IWATTYPE_USP_L & _user_specified_type) {
    tag << 'L';
  }
  if (IWATTYPE_USP_K & _user_specified_type) {
    tag << 'K';
  }
  if (IWATTYPE_USP_G & _user_specified_type) {
    tag << 'G';
  }

  return 1;
}

int
Atom_Typing_Specification::append_to_tag(IWString& tag) const {
  if (0 != _user_specified_type) {
    return _append_tag_for_user_specified_type(tag);
  }

  if (0 == _type) {
    cerr << "Atom_Typing_Specification::append_to_tag:no type specified\n";
    return 0;
  }

  IWString tmp;

  if (!iwattype_convert_to_string_form(_type, tmp)) {
    cerr << "Atom_Typing_Specification::append_to_tag:unrecognised type\n";
    return 0;
  }

  tag << tmp;

  tag.to_uppercase();

  tag << _perform_shell_iteration;

  return 1;
}

int
Atom_Typing_Specification::string_representation(IWString& s) const {
  if (0 == _type) {
    cerr << "Atom_Typing_Specification::string_representation:no type specified\n";
    return 0;
  }

  return iwattype_convert_to_string_form(_type, s);
}

/*
  Pharmacaphore type specified. We need to go and open the file
*/

int
Atom_Typing_Specification::_build_pharmacaphore_specification(
    const const_IWSubstring& s) {
  IWString myfname(s);

  if ("PP" == myfname || 0 == myfname.length()) {
    const char* from_env = getenv("IW_PHARMACAPHORE");
    if (nullptr != from_env) {
      myfname = from_env;
    } else {
      // The default private home directory is removed.
      // A new directory for the queries may be provided in the future
      cerr << "Could not find required queries \n";
      return 0;
    }
  } else if (myfname.starts_with("PP=") || myfname.starts_with("PP:")) {
    myfname.remove_leading_chars(3);
  }

  iwstring_data_source input(myfname.null_terminated_chars());

  if (!input.good()) {
    cerr << "Atom_Typing_Specification::_build_pharmacaphore_specification:cannot open '"
         << myfname << "'\n";
    return 0;
  }

  return _build_pharmacaphore_specification(myfname, input);
}

/*
  In order to achieve installation independence, we must do string interpolation
  on any of the records read.
  We look for any of the tokens that can are shell variables
*/

static int
string_interpolation(const const_IWSubstring& starting_string, IWString& expanded) {
  auto dollar = starting_string.index("${");
  if (dollar < 0) {
    expanded = starting_string;
    return 1;
  }

  const auto n = starting_string.length();

  int closing_brace = -1;

  for (auto i = dollar + 2; i < n; ++i) {
    if ('}' == starting_string[i]) {
      closing_brace = i;
      break;
    }
  }

  if (closing_brace < 0) {  // should we issue a warning if not found
    return 0;
  }

#ifdef DEBUG_STRING_INTERPOLATION
  cerr << "Dollar at " << dollar << " closing_brace " << closing_brace << '\n';
#endif
  const_IWSubstring before_dollar, after_dollar;

  IWString interior;

  starting_string.from_to(0, dollar - 1, before_dollar);
  starting_string.from_to(dollar + 2, closing_brace - 1, interior);
  starting_string.from_to(closing_brace + 1, n - 1, after_dollar);

#ifdef DEBUG_STRING_INTERPOLATION
  cerr << "From '" << starting_string << "'\nget:  '" << before_dollar << "' '"
       << interior << "' and '" << after_dollar << "'\n";
#endif

  // At this point, it might be a single token ${FOO} or maybe a set of possibilities

  const char* e = getenv(interior.null_terminated_chars());

  if (nullptr != e) {
    expanded << before_dollar << e << after_dollar;
    IWString tmp;
    string_interpolation(expanded, tmp);
    expanded = tmp;
    return 1;
  }

  if (!interior.contains(
          '|')) {  // not a single shell variable, no alternatives avaialble.
    return 0;
  }

  IWString token, last_token;
  int i = 0;

  while (interior.nextword(token, i, '|')) {
    e = getenv(token.null_terminated_chars());
#ifdef DEBUG_STRING_INTERPOLATION
    cerr << "Checking shell variable '" << token << "' ? " << (nullptr != e) << '\n';
#endif
    if (nullptr != e) {
      expanded << before_dollar << e << after_dollar;
      return 1;
    }
    last_token = token;
  }

  // We just assume that the last alternative is a default

  expanded << before_dollar << last_token << after_dollar;

  IWString tmp;
  string_interpolation(expanded, tmp);
  expanded = tmp;

  return 1;
}

int
Atom_Typing_Specification::_build_pharmacaphore_specification(
    const IWString& fname, iwstring_data_source& input) {
  input.set_translate_tabs(1);
  input.set_strip_leading_blanks(1);
  input.set_strip_trailing_blanks(1);
  input.set_skip_blank_lines(1);

  auto echo_inputs = false;

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (0 == buffer.length()) {
      continue;
    }

    //  cerr << "Building pp specification '" << buffer << "'\n";
    IWString tmp;

    string_interpolation(buffer, tmp);

    buffer = tmp;

    if (echo_inputs) {
      cerr << "Atom_Typing_Specification::_build_pharmacaphore_specification:reading '"
           << buffer << "'\n";
    }

    //  cerr << " expanded               '" << buffer << "'\n";

    if (buffer == "echo") {
      echo_inputs = true;
    } else if (buffer.starts_with("charge_assigner")) {
      buffer.remove_leading_words(1);

      if (!_charge_assigner.build(buffer)) {
        cerr << "Atom_Typing_Specification::_build_pharmacaphore_specification, cannot "
                "build charge assigner '"
             << buffer << "'\n";
        return 0;
      }
    } else if (buffer.starts_with("donor_acceptor")) {
      buffer.remove_leading_words(1);

      if (!_donor_acceptor_assigner.build(buffer)) {
        cerr << "Atom_Typing_Specification::_build_pharmacaphore_specification, cannot "
                "build donor/acceptor assigner '"
             << buffer << "'\n";
        return 0;
      }
      _donor_acceptor_assigner.set_apply_isotopic_labels(
          0);  // we allocate our own array, do not want to perturb the incoming molecule
    } else if (buffer.starts_with("hydrophobe")) {
      buffer.remove_leading_words(1);
      if (!process_cmdline_token('*', buffer, _hydrophobe, 0)) {
        cerr << "Atom_Typing_Specification::_build_pharmacaphore_specification:cannot "
                "process pharmacaphore '"
             << buffer << "'\n";
        return 0;
      }
    } else if (buffer.starts_with("write")) {
      buffer.remove_leading_words(1);

      _write_isotopically_labelled.add_output_type(FILE_TYPE_SMI);
      if (!_write_isotopically_labelled.new_stem(buffer)) {
        cerr << "Atom_Typing_Specification::_build_pharmacaphore_specification:cannot "
                "initialise write= stream '"
             << buffer << "'\n";
        return 0;
      }
    } else if ("other" ==
               buffer) {  // no directive for what to assign, will get a default type
      _assign_other_type = 1;
    } else if (buffer.starts_with("other")) {
      buffer.remove_leading_words(1);

      buffer.remove_leading_chars(4);

      if (!_parse_user_specified_type(buffer)) {
        return 0;
      }
    } else if ("combine" == buffer) {
      _combine_donor_and_acceptor = 1;
    } else {
      cerr << "Atom_Typing_Specification::_build_hydrophobe_specification:unrecognised "
              "directive '"
           << buffer << "'\n";
      return 0;
    }
  }

  set_global_aromaticity_type(Daylight);

  return 1;
}

QueryAndValue::QueryAndValue() {
  _value = 0;
}

int
QueryAndValue::Build(const atom_typing_spec::QueryAndValue& proto) {
  if (proto.query().empty()) {
    cerr << "QueryAndValue::Build:no query " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! proto.has_value()) {
    cerr << "QueryAndValue::Build:no value " << proto.ShortDebugString() << '\n';
    return 0;
  }

  static constexpr int kVerbose = 0;
  // There is no command line flag here, use arbitrary value.
  static constexpr char kFlag = 'x';

  for (const std::string& query_string : proto.query()) {
    const_IWSubstring token(query_string);
    if (! process_cmdline_token(kFlag, token, _query, kVerbose)) {
      cerr << "QueryAndValue::Build:cannot process '" << token << "'\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  _value = proto.value();

  return _query.size();
}

// Perform all substructure searches against `target`
// and for those atoms that match, which have not been
// previously assigned a value, set the corresponding
// value in `atype`.
template <typename T>
int
QueryAndValue::AssignAtomTypes(Molecule_to_Match& target,
                        T* atype) {
  int rc = 0;
  // cerr << "QueryAndValue has " << _query.size() << " queries\n";
  for (Substructure_Query* q : _query) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      for (atom_number_t a : *e) {
        if (a < 0) {  // Atoms excluded from the embedding.
          continue;
        }
        if (atype[a] != 0) {  // already set.
          continue;
        }

        atype[a] = _value;
        ++rc;
      }
    }
  }

  return 1;
}

int
Atom_Typing_Specification::ParseExternalQuery(IWString& proto_fname) {
  std::optional<atom_typing_spec::External> proto =
        iwmisc::ReadTextProto<atom_typing_spec::External>(proto_fname);

  if (! proto) {
    cerr << "Atom_Typing_Specification::ParseExternalQuery:cannot read '" << proto_fname << '\n';
    return 0;
  }

  return BuildExternal(*proto);
}

int
Atom_Typing_Specification::BuildExternal(const atom_typing_spec::External& proto) {
  for (const atom_typing_spec::QueryAndValue& qv : proto.query_and_value()) {
    std::unique_ptr<QueryAndValue> item = std::make_unique<QueryAndValue>();
    if (! item->Build(qv)) {
      cerr << "Atom_Typing_Specification::BuildExternal:cannot parse ";
      cerr << qv.ShortDebugString() << '\n';
      return 0;
    }

    _query_and_value << item.release();
  }

  return _query_and_value.size();
}

// Run each of the _query_and_value queries against `m` and for
// each of the matched atoms, set the corresponding numeric value
// in `atype`.
// Note that once a value in `atype` is set, it will not be
// overwritten by a subsequent query, so order matters.
template <typename T>
int
Atom_Typing_Specification::AssignExternalAtomTypes(Molecule& m,
        T* atype) {
  std::fill_n(atype, m.natoms(), 0);

  Molecule_to_Match target(&m);
  // cerr << "AssignExternalAtomTypes:have " << _query_and_value.size() << " queries\n";
  for (QueryAndValue* qv : _query_and_value) {
    qv->AssignAtomTypes(target, atype);
  }

  return 1;
}
                        

template int Atom_Typing_Specification::assign_atom_types<int>(Molecule&, int*,
                                                               int const*);
template int Atom_Typing_Specification::assign_atom_types<unsigned int>(Molecule&,
                                                                        unsigned int*,
                                                                        int const*);
template int Atom_Typing_Specification::assign_atom_types<uint64_t>(Molecule&, uint64_t*,
                                                                    int const*);

template int assign_atom_types(Molecule&, int, int*, int const*);
