#include <iostream>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/contiguous_atoms.h"

namespace contiguous_atoms {
int
LargestContiguousCarbonGroup(Molecule& m) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms);

  return LargestContiguousCarbonGroup(m, tmp.get());
}

int
SizeOfCarbonGroup(const Molecule& m, atom_number_t zatom,
                int* visited) {
  visited[zatom] = 1;
  int rc = 1;

  const Atom& a = m[zatom];
  for (const Bond* b : a) {
    atom_number_t o = b->other(zatom);
    if (visited[o]) {
      continue;
    }

    if (m.atomic_number(o) != 6) {
      continue;
    }

    rc += SizeOfCarbonGroup(m, o, visited);
  }

  return rc;
}

int
HeteroatomCount(const Molecule& m, const Ring& r) {
  int rc = 0;
  for (atom_number_t a : r) {
    if (m.atomic_number(a) != 6) {
      ++rc;
    }
  }

  return rc;
}

int
LargestContiguousCarbonGroup(Molecule& m, int* visited) {
  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  // Mark heteroatoms and O=C carbons as visited.
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() == 6) {
      continue;
    }

    visited[i] = 1;

    if (a.ncon() == 1 &&
        (a.atomic_number() == 7 || a.atomic_number() == 8)) {
      const Bond* b = a[0];
      if (! b->is_single_bond()) {
        visited[b->other(i)] = 1;
      }
    }
  }

  for (const Ring* r : m.sssr_rings()) {
    if (! r->is_aromatic()) {
      continue;
    }

    if (HeteroatomCount(m, *r) > 0) {
      r->set_vector(visited, 1);
    }
  }

  // anything attached to a heteroatom is excluded.
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() == 6) {
      continue;
    }

    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      if (m.is_halogen(o)) {
        continue;
      }
      visited[b->other(i)] = 1;
    }
  }

  int largest_size = 0;

  for (int i = 0; i < matoms; ++i) {
    if (visited[i]) {
      continue;
    }

    if (m.atomic_number(i) != 6) {
      continue;
    }

    if (int s = SizeOfCarbonGroup(m, i, visited); s > largest_size) {
      largest_size = s;
    }
  }

  return largest_size;
}


}  // namespace contiguous_atoms
