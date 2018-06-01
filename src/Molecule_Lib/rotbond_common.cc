/*
  Rotatable bonds are computed in iwdescr and in the rotatable_bonds programme.
  In order to keep these in sync as much as possible, we include as many common
  functions here as possible
*/

#include <stdlib.h>

#include "misc.h"

#include "molecule.h"
#include "rotbond_common.h"

/*
  A bond might be considered rotatable. Is there a triple bond at either end?
*/

int
triple_bond_at_either_end(const Molecule & m,
                          const Bond * b)
{
  Atom * a1 = const_cast<Atom *>(m.atomi(b->a1()));     // nbonds() is non-const

  if (2 == a1->ncon() && 4 == a1->nbonds())
    return 1;

  Atom * a2 = const_cast<Atom *>(m.atomi(b->a2()));     // nbonds() is non-const

  if (2 == a2->ncon() && 4 == a2->nbonds())
    return 1;

  return 0;        // nope, no triple bonds here
}

int
is_non_rotatable_amide(Molecule & m,
                       atom_number_t n,
                       atom_number_t c)
{
  const Atom * nitrogen = m.atomi(n);
  const Atom * carbon = m.atomi(c);

// Swap things around to get the pointers pointing to the right atoms

  if (7 == nitrogen->atomic_number() && 6 == carbon->atomic_number())    // great, got it first time
    ;
  else if (6 == nitrogen->atomic_number() && 7 == carbon->atomic_number())
  {
    std::swap(n, c);
    std::swap(nitrogen, carbon);
  }
  else
    return 0;

  if (2 != nitrogen->ncon())    // only O=C-[NH]- is non-rotatable
    return 0;

  if (1 != m.hcount(n))    // probably not necessary to test
    return 0;

  int ccon = carbon->ncon();
  int cbonds = carbon->nbonds();

  if (3 == ccon && 4 == cbonds)      // [NH]-C(=O)-*
    ;
  else if (2 == ccon && 3 == cbonds)    // [NH]-[CH]=O
    ;
  else
    return 0;

  int doubly_bonded_oxygens = 0;

  for (int i = 0; i < ccon; i++)
  {
    const Bond * b = carbon->item(i);
    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(c);

    atomic_number_t z = m.atomic_number(o);
    if (8 == z)
      ;
//  else if (16 == z)    // allow sulphur
//    ;
    else
      return 0;

    doubly_bonded_oxygens++;
  }

  return 1 == doubly_bonded_oxygens;
}

static int
is_cf3_or_t_butyl(Molecule & m,
                  const atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int n = 0;    // the number of F or C's attached

  atomic_number_t z = INVALID_ATOMIC_NUMBER;

  for (int i = 0; i < 4; i++)
  {
    atom_number_t j = a->other(zatom, i);

    const Atom * aj = m.atomi(j);

    if (1 != aj->ncon())
      continue;

    if (INVALID_ATOMIC_NUMBER == z)
    {
      z = aj->atomic_number();
      n = 1;
    }
    else if (z == aj->atomic_number())
      n++;
  }

  return 3 == n;    // we need 3 singly connected C or F neighbours
}

static int
is_cf3_or_t_butyl (Molecule & m,
                   atom_number_t zatom1,
                   atom_number_t zatom2)
{
  const Atom * a1 = m.atomi(zatom1);

  if (4 == a1->ncon() && 4 == a1->nbonds())
  {
    if (is_cf3_or_t_butyl(m, zatom1))
      return 1;
  }

  const Atom * a2 = m.atomi(zatom2);

  if (4 == a2->ncon() && 4 == a2->nbonds())
  {
    if (is_cf3_or_t_butyl(m, zatom2))
      return 1;
  }


  return 0;
}

int
is_non_rotatable_sulphonamide(Molecule & m,
                              atom_number_t zatom1,
                              atom_number_t zatom2)
{
  const Atom * a1 = m.atomi(zatom1);
  const Atom * a2 = m.atomi(zatom2);

  if (16 == a1->atomic_number() && 7 == a2->atomic_number())
    ;
  else if (7 == a1->atomic_number() && 16 == a2->atomic_number())
  {
    std::swap(zatom1, zatom2);
    std::swap(a1, a2);
  }
  else
    return 0;

#ifdef NOT_SURE_IF_THESE_ARE_NEEDED_OR_NOT
  if (2 != a2->ncon())    // only O=C-[NH]- is non-rotatable
    return 0;

  if (1 != m.hcount(zatom2))    // probably not necessary to test
    return 0;
#endif

  if (4 != a1->ncon())
    return 0;

  if (6 != a1->nbonds())
    return 0;

  int doubly_bonded_oxygen = 0;

  for (int i = 0; i < 4; i++)
  {
    const Bond * b = a1->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(zatom1);

    if (1 != m.ncon(o))
      continue;

    if (8 != m.atomic_number(o))
      continue;

    doubly_bonded_oxygen++;
  }

//cerr << "doubly_bonded_oxygen " << doubly_bonded_oxygen  << endl;
  return 2 == doubly_bonded_oxygen;
}

/*
  The bond between A1 and A2 might be rotatable. Is it part of a CF3, t-Butyl or O=C-[NH]-
*/

int
part_of_otherwise_non_rotabable_entity(Molecule & m,
                                       atom_number_t a1,
                                       atom_number_t a2)
{
//cerr << " atom " << m.smarts_equivalent_for_atom(a1) << ' ' << m.smarts_equivalent_for_atom(a2) << " is cf3 " << is_cf3_or_t_butyl(m, a1, a2) << endl;

  if (is_cf3_or_t_butyl(m, a1, a2))
    return 1;

  if (is_non_rotatable_amide(m, a1, a2))    // O=C-[NH]- is non-rotatable
    return 1;

  if (is_non_rotatable_sulphonamide(m, a1, a2))    // O=S(=O)-[NH]- is non-rotatable
    return 1;

  return 0;
}
