#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#define COMPILING_IS_ACTUALLY_CHIRAL_CC

#include "chiral_centre.h"
#include "is_actually_chiral.h"
#include "molecule.h"
#include "path_scoring.h"
#include "smiles.h"

using std::cerr;

static int max_iterations = std::numeric_limits<int>::max();

void
set_max_iterations(int m) {
  assert(m > 0);

  max_iterations = m;
}

static int allow_unsaturated_atoms_to_be_chiral = 0;

void
set_allow_unsaturated_atoms_to_be_chiral(int s) {
  allow_unsaturated_atoms_to_be_chiral = s;
}

/*
  To determine if an atom is chiral or not, we need to perform path tracing
  from that atom.
*/

static int
is_actually_chiral(Molecule& m, atom_number_t zatom, resizable_array_p<Path_Scoring>& ps,
                   int* claimed, Atom* const* atom) {
  const Atom* a = atom[zatom];

  int acon = a->ncon();

  if (ps.number_elements()) {
    ps.resize_keep_storage(0);
  }

  ps.resize(acon);

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    atom_number_t j = b->other(zatom);

    Path_Scoring* p = new Path_Scoring;

    const Atom* aj = atom[j];

    p->initialise(j, aj);

    if (b->is_single_bond()) {
      p->set_first_bond(1);
    } else if (!allow_unsaturated_atoms_to_be_chiral || a->atomic_number() < 14) {
      delete p;
      return 0;
    } else if (b->is_double_bond()) {  // tetrahedral Sulphur types
      p->set_first_bond(2);
    }

    claimed[j] = 1;

    ps.add(p);
  }

  int stopped;
  if (resolved(ps, stopped)) {
    return 1;
  }

  // cerr << "Begin path expansions around atom " << zatom << ' ' <<
  // m.smarts_equivalent_for_atom(zatom) << '\n';

  for (int iterations = 0; iterations < max_iterations; ++iterations) {
    for (int i = 0; i < acon; i++) {
      if (ps[i]->active()) {
        ps[i]->advance(atom, claimed);
      }
    }

    int stopped;
    if (resolved(ps, stopped)) {
      return 1;
    }

    if (stopped) {  // not resolved, but cannot go any further
      return 0;
    }

    int number_active = 0;
    for (int i = 0; i < acon; i++) {
      if (!ps[i]->active()) {
        continue;
      }

      ps[i]->update_claimed(claimed);
      number_active++;
    }

    if (number_active < 2) {
      return 0;
    }
  }

  return 1;
}

/*
  The query for an asymmetric carbon atom will hit things like the
  carbon in t-butyl. We need to examine the neighbours to make sure
  that this atom actually is an asymmetric centre

  Thought about putting in a more aggressive check on the number of
  connections, but too dangerous. Even 2 == ncon is problematic because
  you could have an atom with a lone-pair and an implicit Hydrogen
*/

int
is_actually_chiral(Molecule& m, atom_number_t zatom) {
  resizable_array_p<Path_Scoring> ps;

  return is_actually_chiral(m, zatom, ps);
}

int
is_actually_chiral(Molecule& m, atom_number_t zatom,
                   resizable_array_p<Path_Scoring>& ps) {
  const Atom* a = m.atomi(zatom);

  const int acon = a->ncon();

  if (acon < 2 || acon > 4) {
    return 0;
  }

  const int hcount = m.hcount(zatom);

  if (hcount > 1) {  // what if isotopic Hydrogen???
    return 0;
  }

  int lp;
  if (m.lone_pair_count(zatom, lp) && lp > 1) {
    return 0;
  }

  if (1 == hcount && 1 == lp && 7 == a->atomic_number()) {  // never
    return 0;
  }

  if (acon < a->nbonds() && !allow_unsaturated_atoms_to_be_chiral) {
    return 0;
  }

  m.compute_aromaticity_if_needed();  // so bonds get aromatic character

  if (m.is_aromatic(zatom)) {
    return 0;
  }

  const int matoms = m.natoms();

  int* claimed = new_int(matoms);
  std::unique_ptr<int[]> free_claimed(claimed);

  claimed[zatom] = 1;

  Atom* const* atoms = new Atom*[matoms];
  std::unique_ptr<Atom* const[]> free_atoms(atoms);

  m.atoms((const Atom**)atoms);

  // cerr << "Detailed calculation on " << m.smarts_equivalent_for_atom(zatom) << '\n';

  return is_actually_chiral(m, zatom, ps, claimed, atoms);
}

namespace lillymol {

int
do_remove_invalid_chiral_centres(Molecule& m) {
  int nc = m.chiral_centres();
  if (0 == nc) {
    return 0;
  }

  // Removing a chiral centre while we are scanning the set would mess things up,
  // so we make a list of the atoms with invalid chiral centres and remove them later

  Set_of_Atoms centres_to_be_removed;

  for (int i = 0; i < nc; i++) {
    Chiral_Centre* c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    atom_number_t a = c->a();

    // cerr << "Atom chiral? " << is_actually_chiral(m, a) << ' ' <<
    // m.smarts_equivalent_for_atom(a) << '\n';

    if (!is_actually_chiral(m, a)) {
      centres_to_be_removed.add(a);
    }
  }

  if (centres_to_be_removed.number_elements()) {
    for (int i = 0; i < centres_to_be_removed.number_elements(); i++) {
      m.remove_chiral_centre_at_atom(centres_to_be_removed[i]);
    }
  }

  return centres_to_be_removed.number_elements();
}

} // namespace lillymol

std::ostream&
operator<<(std::ostream& os, const CahnIngoldPrelog& cip) {
  switch (cip) {
    case CahnIngoldPrelog::kNeither:
      os << "Neither";
      return os;
    case CahnIngoldPrelog::R:
      os << 'R';
      return os;
    case CahnIngoldPrelog::S:
      os << 'S';
      return os;
    case CahnIngoldPrelog::kUnspecified:
      os << "Unspecified";
      return os;
    default:
      os << '?';
      return os;
  }
}

std::optional<CahnIngoldPrelog>
Molecule::CahnIngoldPrelogValue(atom_number_t zatom) {
  const Chiral_Centre* c = chiral_centre_at_atom(zatom);
  if (c == nullptr) {
    return std::nullopt;
  }

  return CahnIngoldPrelogValue(c);
}

// Convert zatom to an atomic number equivalent for CIP perception.
// `zatom` must be either
//   a non-atom part of a chiral centre
//   a valid atom number.
// In the case of an implicit hydrogen, return 1.
// In the case of a lone pair, return 0;
std::optional<uint32_t>
Molecule::ChiralCentreMemberToCipInt(int zatom) const {
  if (zatom >= 0 && zatom < _number_elements) {
    return _things[zatom]->atomic_number();
  }

  if (zatom == kChiralConnectionIsImplicitHydrogen) {
    return 1;
  }

  if (zatom == kChiralConnectionIsLonePair) {
    return 0;
  }

  return std::nullopt;
}

CahnIngoldPrelog
Molecule::CahnIngoldPrelogValue(const Chiral_Centre* c) {
  CahnIngoldPrelog rc = CahnIngoldPrelog::kNeither;

  auto top_front = ChiralCentreMemberToCipInt(c->top_front());
  if (!top_front) {
    return rc;
  }

  auto top_back = ChiralCentreMemberToCipInt(c->top_back());
  if (!top_back) {
    return rc;
  }

  auto left_down = ChiralCentreMemberToCipInt(c->left_down());
  if (!left_down) {
    return rc;
  }

  auto right_down = ChiralCentreMemberToCipInt(c->right_down());
  if (!right_down) {
    return rc;
  }

  return CahnIngoldPrelogValue(c, *top_front, *top_back, *left_down, *right_down);
}

// clang-format off
enum Position { kTopFront = 0,
                kTopBack = 1,
                kLeftDown = 2,
                kRightDown = 3
};
// clang-format on

// An arbitrary value indicating that during an expansion, an atom is the centre
// of the Chiral_Centre being resolved.
static constexpr int kCentre = -11;
// an arbitrary negative value indicating that an atom has already been visited
// as part of a shell expansion.
static constexpr int kVisited = -12;
static constexpr int kThisExpansion = -13;

// `top_front` etc are all atomic number equivalents.
CahnIngoldPrelog
Molecule::CahnIngoldPrelogValue(const Chiral_Centre* c, int top_front, int top_back,
                                int left_down, int right_down) {
#ifdef DEBUG_CAHN_INGOLD_PRELOG
  cerr << "top_front " << top_front << " top_back " << top_back << " left_down "
       << left_down << " right_down " << right_down << '\n';
#endif

  // See if resolved by the atoms directly attached.
  if (top_front == top_back || top_front == left_down || top_front == right_down ||
      top_back == left_down || top_back == right_down || left_down == right_down) {
    // a least two atoms the same, resolve by expansion below.
  } else if (top_front < top_back && top_front < left_down && top_front < right_down) {
    return CahnIngoldPrelogValue(top_back, left_down, right_down);
  } else if (top_back < top_front && top_back < left_down && top_back < right_down) {
    return CahnIngoldPrelogValue(top_front, right_down, left_down);
  } else if (right_down < top_front && right_down < top_back && right_down < left_down) {
    return CahnIngoldPrelogValue(top_front, left_down, top_back);
  } else if (left_down < top_front && left_down < top_back && left_down < right_down) {
    return CahnIngoldPrelogValue(top_front, top_back, right_down);
  } else {
    cerr << "Molecule::CahnIngoldPrelogValue:internal error\n";
    return CahnIngoldPrelog::kNeither;
  }

  cerr << "Molecule::CahnIngoldPrelogValue:shell expansion not implemented\n";
  CahnIngoldPrelog rc = CahnIngoldPrelog::kNeither;
  return rc;

  // Atom needs to be resolved by expansion.
  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(4 * _number_elements);
  std::fill_n(tmp.get(), 4 * _number_elements, 0);
  // For each expansion, mark the centre atom.
  for (int i = 0; i < 4; ++i) {
    tmp.get()[i * _number_elements + c->a()] = kCentre;
  }

  // TODO:ianwatson implement shell expansion
}

// The low priority connection has been identified, and we can set up
// an ordering of the remaining connections.
//       north                  *
//         |                    *
//         *                    *
//        / \                   *
//      /    \                  *
//    sw      se                *
// Examine the relative values and return R or S.
CahnIngoldPrelog
Molecule::CahnIngoldPrelogValue(int north, int se, int sw) const {
  if (north > se && se > sw) {
    return CahnIngoldPrelog::R;
  } else if (se > sw && sw > north) {
    return CahnIngoldPrelog::R;
  } else if (sw > north && north > se) {
    return CahnIngoldPrelog::R;
  } else if (north > sw && sw > se) {
    return CahnIngoldPrelog::S;
  } else if (sw > se && se > north) {
    return CahnIngoldPrelog::S;
  } else if (se > north && north > sw) {
    return CahnIngoldPrelog::S;
  } else {
    cerr << "CahnIngoldPrelogValue:internal error, north " << north << " se " << se
         << " sw " << sw << '\n';
    return CahnIngoldPrelog::kUnspecified;
  }
}

namespace lillymol {

// return true if any two atoms attached to `zatom` are symmetric.
int
AnySymmetry(Molecule& m,
            atom_number_t zatom,
            const int* symm) {
  resizable_array<int> symm_values;

  // cerr << "AnySymmetry from atom " << zatom << "\n";
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    // cerr << "From atom " << zatom << " sym " << symm[zatom]  << " to " << o << " symm " << symm[o] << ' ' << m.smarts_equivalent_for_atom(o) << '\n';
    // If we cannot add this symmetry flag to the array of existing values,
    // this connection must be the same as a previous connection, so there
    // is symmetry at `zatom`.
    // cerr << " to atom " << o << " symm " << symm[o] << '\n';
    if (!symm_values.add_if_not_already_present(symm[o])) {
      return 1;
    }
  }

  return 0;
}

// The default symmetry calculation takes into account chirality, so
// we need to temporarily turn off consideration of chirality in
// forming the symmetry classes. This begs the question what is actually
// correct. A lot of this is motivated by molecules like
// N1C[C@@H](O)[C@@H](O)C(O)C1 CHEMBL313657
// Is the chiral centre opposite the Nitrogen atom valid or not?
// Since there is an unmarked chiral centre, we really do not know.
// But does this decision compromise cases where the chirality is known?
int
RemoveInvalidChiralCentresUsingSymmetry(Molecule & m) {
  const int nc = m.chiral_centres();
  if (nc == 0) {
    return 0;
  }

  const int matoms = m.natoms();
  set_include_chiral_info_in_smiles(0);
  m.invalidate_smiles();
  std::unique_ptr<int[]> symm = std::make_unique<int[]>(matoms);
  std::copy_n(m.symmetry_classes(), matoms, symm.get());
  set_include_chiral_info_in_smiles(1);
  m.invalidate_smiles();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Chiral_Centre* c = m.chiral_centre_at_atom(i);
    if (c == nullptr) {
      continue;
    }

    if (AnySymmetry(m, i, symm.get())) {
      m.remove_chiral_centre_at_atom(i);
      ++rc;
    }
  }

  return rc;
}

int
IsActuallyChiralBySymmetry(Molecule& m, atom_number_t zatom) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> symm = std::make_unique<int[]>(matoms);
  std::copy_n(m.symmetry_classes(), matoms, symm.get());

  if (AnySymmetry(m, zatom, symm.get())) {
    return 0;
  }

  return 1;
}

std::unique_ptr<ChiralStatus[]>
ChiralityStatus(Molecule& m) {
  const int matoms = m.natoms();
  std::unique_ptr<ChiralStatus[]> result = std::make_unique<ChiralStatus[]>(matoms);

  set_include_chiral_info_in_smiles(0);
  m.invalidate_smiles();
  std::unique_ptr<int[]> symm = std::make_unique<int[]>(matoms);
  std::copy_n(m.symmetry_classes(), matoms, symm.get());
  set_include_chiral_info_in_smiles(1);
  m.invalidate_smiles();

  for (int i = 0; i < matoms; ++i) {
    result[i] = ChiralStatus::kNotChiral;

    const Atom& a = m[i];
    if (a.ncon() < 3) {
      continue;
    }

    if (a.ncon() == 4) {
      if (allow_unsaturated_atoms_to_be_chiral) {
      } else if (! m.saturated(i)) {
        continue;
      }
    } else if (a.ncon() == 3 && m.hcount(i) == 1) {
    } else if (allow_unsaturated_atoms_to_be_chiral && ! m.saturated(i)) {
    } else {
      continue;
    }

    const Chiral_Centre* c = m.chiral_centre_at_atom(i);
    if (AnySymmetry(m, i, symm.get())) {
      if (c) {
        result[i] = ChiralStatus::kInvalid;
      }
    } else {
      if (c) {
        result[i] = ChiralStatus::kChiral;
      } else {
        result[i] = ChiralStatus::kUnmarked;
      }
    }
  }

  return result;
}

}  // namespace lillymol
