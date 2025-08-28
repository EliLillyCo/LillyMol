#include <iostream>
#include <memory>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/set_of_atoms.h"

#include "subset_of_interest.h"

namespace subset_of_interest {

using std::cerr;

SubsetOfInterest::SubsetOfInterest() {
  _isotope_at_attachment_points = 0;
  _break_bonds_to_aromatic_rings = 0;
}

int
SubsetOfInterest::ReduceToFragmentOfInterest(Molecule& m, atom_number_t zatom){
  // The atoms at the ends of the bonds to be removed.
  Set_of_Atoms lhs, rhs;
  Set_of_Atoms aromatic_atom;

  m.ring_membership();

#ifdef DEBUG_REDUCE_TO_FRAGMENT_OF_INTEREST
  cerr << "_break_bonds_to_aromatic_rings " << _break_bonds_to_aromatic_rings <<  '\n';
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << '\n';
  cerr << "Must keep atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << '\n';
#endif


  for (const Bond* b : m.bond_list()) {
    if (! b->is_single_bond()) {
      continue;
    }
    if (b->nrings()) {
      continue;
    }
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (a1 == zatom|| a2 == zatom) [[unlikely]] {
      continue;
    }

    if (m.atomic_number(a1) != 6 || m.atomic_number(a2) != 6) {
      continue;
    }

    if (m.saturated(a1) && m.saturated(a2)) {
      // bond will be removed.
    } else if (_break_bonds_to_aromatic_rings == 0) {
      continue;
    } else if (m.is_aromatic(a1)) {
      aromatic_atom << a2;
    } else if (m.is_aromatic(a2)) {
      aromatic_atom << a1;
    } else {
      continue;
    }

    lhs << a1;
    rhs << a2;
  }

  if (lhs.empty()) {
    return 0;
  }

  //cerr << lhs << " lhs being removed\n";
  for (int i = lhs.number_elements() - 1; i >= 0; --i) {
    m.remove_bond_between_atoms(lhs[i], rhs[i]);
    if (_isotope_at_attachment_points) {
      m.set_isotope(lhs[i], _isotope_at_attachment_points);
      m.set_isotope(rhs[i], _isotope_at_attachment_points);
    }
  }

  //cerr << "ARomatic " << aromatic_atom << '\n';
  if (aromatic_atom.size() > 0) {
    for (atom_number_t a : aromatic_atom) {
      const isotope_t iso = m.isotope(a);
      m.set_isotope(a, iso + _break_bonds_to_aromatic_rings);
    }
  }

  int keep_frag = m.fragment_membership(zatom);

  const int matoms = m.natoms();

  std::unique_ptr<int[]> to_remove = std::make_unique<int[]>(matoms);
  //cerr << m.smiles() << " before removing atoms\n";

  for (int i = 0; i < matoms; ++i) {
    if (m.fragment_membership(i) == keep_frag) {
      to_remove[i] = 0;
    } else {
      to_remove[i] = 1;
    }
    //cerr << i << " to_remove " << to_remove[i] << '\n';
  }

  m.remove_atoms(to_remove.get());

  return 1;
}

};  // namespace subset_of_interest

