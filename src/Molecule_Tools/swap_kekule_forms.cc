#include <stdlib.h>

#include "molecule.h"
#include "path.h"
#include "toggle_kekule_form.h"

int
find_adjacent_atoms_in_common_between_two_rings (const Ring & ri, 
                      const Ring & rj,
                      atom_number_t & a1,
                      atom_number_t & a2)
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  for (Ring_Bond_Iterator i(ri); i != ri.end(); ++i)
  {
    if (! rj.contains_bond(i.a1(), i.a2()))
      continue;

    a1 = i.a1();
    a2 = i.a2();

    return 1;
  }

  return 0;   // hard to imagine
}

int
arrange_kekule_forms_in_fused_rings (Molecule & m)
{
  int nr = m.nrings();

  if (nr < 2)
    return 0;

//cerr << "arrange_kekule_forms_in_fused_rings: checking " << nr << " rings\n";

  m.compute_aromaticity_if_needed();

  Set_of_Atoms f1, f2;    // we will store the atom pairs at the ends of fused aromatic rings

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    int nf = ri->fused_ring_neighbours();

    if (0 == nf)
      continue;

    for (int j = 0; j < nf; j++)
    {
      const Ring * rj = ri->fused_neighbour(j);

      atom_number_t a1, a2;

      if (! find_adjacent_atoms_in_common_between_two_rings(*ri, *rj, a1, a2))
        continue;

      const Bond * b = m.bond_between_atoms(a1, a2);

      if (b->is_double_bond())
        continue;

      if (7 == m.atomic_number(a1) || 7 == m.atomic_number(a2))  // these cannot be toggled
        continue;

      f1.add(a1);
      f2.add(a2);
    }
  }

// Now we have the bonds identified, start switching things

  int n = f1.number_elements();
  assert (n == f2.number_elements());

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t a1 = f1[i];
    atom_number_t a2 = f2[i];

#ifdef DEBUG_ARRANGE_KEKULE_FORMS_IN_FUSED_RINGS
    cerr << "Atoms in common " << a1 << " and " << a2 <<endl;
#endif

    const Bond * b = m.bond_between_atoms(a1, a2);

    if (b->is_double_bond())   // great, already set
      continue;

    Toggle_Kekule_Form tkf;

    int changed = 0;

#ifdef DEBUG_ARRANGE_KEKULE_FORMS_IN_FUSED_RINGS
      IWString ismi = m.smiles();
#endif

    if (2 != m.nrings(a1) || 2 != m.nrings(a2))
      continue;

    if (! tkf.process(m, a1, a2, DOUBLE_BOND, changed))
      continue;

    if (! m.valence_ok())
    {
      cerr << "Fatal error, bad valence in Kekule switch " << m.smiles() << ' ' << m.name() << endl;
    }

#ifdef DEBUG_ARRANGE_KEKULE_FORMS_IN_FUSED_RINGS
    cerr << "Switched " << ismi << " to " << m.smiles() << endl;
#endif

    rc++;
  }

  if (rc)
    m.compute_aromaticity_if_needed();

  return rc;
}

