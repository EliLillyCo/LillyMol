#include <algorithm>

// This is not working properly. This molecule generates different
// results depending on the smiles ordering. Might be something
// conceptually wrong about the algorithm.
// TODO:ianwatson resolve this sometime.
// O=C(NC)[C@H]1NC(=O)CC1 CHEMBL1892080
// N([1CH3])[2C](=[3O])[4CH]1[5CH2][6CH2][7C](=[8O])[9NH]1 CHEMBL1892080

#include "Foundational/iwmisc/misc.h"

#include "partial_symmetry.h"

namespace partial_symmetry {

constexpr int kRingFormed = 46;

void
AssignIfGreater(int new_value, int& existing_value) {
  if (new_value > existing_value) {
    existing_value = new_value;
  }
}

// Return a score associated with atom `zatom` in `m` being found
// in the next layer via bond `bond`.
// No allowance made for non periodic table elements.
int
ComputeScore(const Molecule& m,
             atom_number_t zatom,
             const Bond& bond) {
  if (bond.is_aromatic()) {
    return m.atomic_number(zatom);
  }
  if (bond.is_single_bond()) {
    return 120 + m.atomic_number(zatom);
  }
  if (bond.is_double_bond()) {
    return 240 + m.atomic_number(zatom);
  }
  if (bond.is_triple_bond()) {
    return 360 + m.atomic_number(zatom);
  }

  return 1;
}

PartialSymmetry::PartialSymmetry(Molecule& m) : _m(m) {
  _matoms = m.natoms();
  if (_matoms == 0) {
    _score = nullptr;
    _symmetric_at_radius = nullptr;
    _computation_done = true;
    return;
  }

  _score = new InLayer[_matoms];
  _symmetric_at_radius = new_int(_matoms);

  _current_shell.resize(_matoms);
  _next_shell.resize(_matoms);
    
  _computation_done = false;
}

// Deleting a nullptr is harmless.
PartialSymmetry::~PartialSymmetry() {
  delete [] _score;
  delete [] _symmetric_at_radius;
}

void
PartialSymmetry::AllScoresZero() {
  for (int i = 0; i < _matoms; ++i) {
    _score[i].score = 0;
    _score[i].atom = i;
  }
}

// #define DEBUG_PARTIAL_SYMMETRY

// Examine the _score array and for every pair that consists of a
// bonded pair of atoms, adjust their score.
void
PartialSymmetry::AdjustForRings(const int * dm) {

  for (int i = 0; i < _matoms; ++i) {
    if (_score[i].score == 0) {
      continue;
    }
    for (int j = i + 1; j < _matoms; ++j) {
      if (_score[j].score == 0) {
        continue;
      }

      if (dm[i * _matoms + j] == 1) {
        _score[i].score += kRingFormed + 1;
        _score[j].score += kRingFormed + 1;
      }
    }
  }
}

// Expand the shell, that started at `starting_atom` into atoms at
// `radius` bonds from `starting_atom`.
int
PartialSymmetry::Expand(atom_number_t starting_atom, int radius) {
  if (_current_shell.empty()) {
    return 1;
  }

  _next_shell.resize_keep_storage(0);
#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "PartialSymmetry::Expand:";
  for (int i : _current_shell) {
    std::cerr << ' ' << i;
  }
  std::cerr << '\n';
#endif
  AllScoresZero();

  int atoms_in_next_shell = 0;
  const int * dm = _m.distance_matrix_warning_may_change();
  for (const atom_number_t i : _current_shell) {
#ifdef DEBUG_PARTIAL_SYMMETRY
    std::cerr << " atom " << i << " is in the current shell\n";
#endif
    const Atom& atom = _m.atom(i);
    for (const Bond * b : atom) {
      const atom_number_t j = b->other(i);
#ifdef DEBUG_PARTIAL_SYMMETRY
      std::cerr << " is atom " << j << " in radius " << radius << " " << dm[starting_atom * _matoms + j] << '\n';
#endif
      if (dm[starting_atom * _matoms + j] != radius) {
        continue;
      }
      _score[j].score += ComputeScore(_m, j, *b);
      atoms_in_next_shell++;
    }
  }

#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "atoms_in_next_shell " << atoms_in_next_shell << '\n';
#endif

  // Finding 0 or 1 atoms in the next shell, means we are done.
  if (atoms_in_next_shell < 2) {
    return 1;
  }

  AdjustForRings(dm);

  // Also fills _next_shell.
  SortAndAssign(radius);

#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "Next shell contains " << _next_shell.size() << " items\n";
#endif
  if (_next_shell.empty()) {
    return 1;
  }

  _current_shell.swap_contents(_next_shell);

  return Expand(starting_atom, radius + 1);
}

// Expand from a bond defined by atoms `starting_atom1` 
// and `starting_atom2`.
int
PartialSymmetry::Expand(atom_number_t starting_atom1, atom_number_t starting_atom2, int radius) {
  if (_current_shell.empty()) {
    return 1;
  }

  _next_shell.resize_keep_storage(0);
#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "PartialSymmetry::Expand2:";
  for (int i : _current_shell) {
    std::cerr << ' ' << i;
  }
  std::cerr << '\n';
#endif
  AllScoresZero();

  int atoms_in_next_shell = 0;
  const int * dm = _m.distance_matrix_warning_may_change();
  for (const atom_number_t i : _current_shell) {
#ifdef DEBUG_PARTIAL_SYMMETRY
    std::cerr << " atom " << i << " is in the current shell\n";
#endif
    const Atom& atom = _m.atom(i);
    for (const Bond * b : atom) {
      const atom_number_t j = b->other(i);
      if (radius == 0 && (j == starting_atom1 || j == starting_atom2)) {
        continue;
      }
#ifdef DEBUG_PARTIAL_SYMMETRY
      std::cerr << " is atom " << j << " in radius " << radius << " " << dm[starting_atom1 * _matoms + j] << " and " << dm[starting_atom2 * _matoms + j] << '\n';
#endif
      if (dm[starting_atom1 * _matoms + j] != radius &&
          dm[starting_atom2 * _matoms + j] != radius) {
        continue;
      }
      _score[j].score += ComputeScore(_m, j, *b);
      atoms_in_next_shell++;
    }
  }

#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "atoms_in_next_shell " << atoms_in_next_shell << '\n';
#endif

  // Finding 0 or 1 atom in the next shell, means we are done.
  if (atoms_in_next_shell < 2) {
    return 1;
  }

  // Also fills _next_shell.
  SortAndAssign(radius);

#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "Next shell contains " << _next_shell.size() << " items\n";
#endif
  if (_next_shell.empty()) {
    return 1;
  }

  _current_shell.swap_contents(_next_shell);

  return Expand(starting_atom1, starting_atom2, radius + 1);
}


const int *
PartialSymmetry::SymmetricAtRadius() {
  if (_computation_done) {
    return _symmetric_at_radius;
  }

  if (_matoms == 0) {
    return nullptr;
  }

  _m.compute_aromaticity_if_needed();

  for (int i = 0; i < _matoms; ++i) {
    _current_shell.resize_keep_storage(0);
    _current_shell << i;
#ifdef DEBUG_PARTIAL_SYMMETRY
    std::cerr << "Begin expansion from atom " << i << '\n';
#endif
    Expand(i, 1);
  }

  for (const Bond * b : _m.bond_list()) {
    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    if (_m.atomic_number(a1) != _m.atomic_number(a2)) {
      continue;
    }
    _current_shell.resize_keep_storage(0);
    _current_shell << a1 << a2;
    Expand(a1, a2, 0);
  }

  _computation_done = true;

  return _symmetric_at_radius;
}

// The _score array is filled with atoms at radius `radius`.
// Compress out the values not in use, sort the _score array
// and if there are duplicate values present, update _symmetric_at_radius.
int
PartialSymmetry::SortAndAssign(int radius) {
  // Remove atoms that are not part of this shell.
  int ndx = 0;
  for (int i = 0; i < _matoms; ++i) {
    if (_score[i].score == 0) {
      continue;
    }
    if (i > ndx) {
      _score[ndx] = _score[i];
    }
    ndx++;
  }

#ifdef DEBUG_PARTIAL_SYMMETRY
  std::cerr << "Actual count next shell " << ndx << '\n';
#endif
  assert(ndx > 0);

  // Just one atom, cannot be any duplicate values, nothing to do.
  if (ndx == 1) {
    return 1;
  }

  // Handle a common case.
  if (ndx == 2) {
    if (_score[0].score != _score[1].score) {
        return 1;
    }
    const atom_number_t a0 = _score[0].atom;
    const atom_number_t a1 = _score[1].atom;
    AssignIfGreater(radius, _symmetric_at_radius[a0]);
    AssignIfGreater(radius, _symmetric_at_radius[a1]);
    _next_shell << a0 << a1;
    return 1;
  }

  std::sort(_score, _score + ndx, [](const InLayer& s1, const InLayer& s2) {
    return s1.score < s2.score;
  });

  // All items that have the same score are symmetric at this radius.

  if (_score[0].score == _score[1].score) {
    AssignIfGreater(radius, _symmetric_at_radius[_score[0].atom]);
    _next_shell << _score[0].atom;
  }
  for (int i = 1; i < ndx; ++i) {
    if (_score[i].score == _score[i-1].score) {
      AssignIfGreater(radius, _symmetric_at_radius[_score[i].atom]);
      _next_shell << _score[i].atom;
    }
  }

  return 1;
}

}  // namespace partial_symmetry
