#include "Molecule_Tools/dicer_lib.h"

namespace dicer_lib {

USPVPTR::USPVPTR() {
  _atom_number_in_parent = -1;
  _atom_type_in_parent = 0;
}

std::optional<atom_number_t>
InitialAtomNumber(const Atom* a) {
  const void * v = a->user_specified_void_ptr();
  if (v == nullptr) {
    return std::nullopt;
  }

  const USPVPTR* uspvptr = reinterpret_cast<const USPVPTR*>(v);
  return uspvptr->_atom_number_in_parent;
}


int
initialise_atom_pointers(Molecule & m,
                         const uint32_t * atype,
                         USPVPTR * uspvptr)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    uspvptr[i].set_atom_number_in_parent(i);
    m.set_user_specified_atom_void_ptr(i, (void*)(uspvptr + i));
  }

  if (nullptr != atype) {
    for (int i = 0; i < matoms; ++i) {
      uspvptr[i].set_atom_type_in_parent(atype[i]);
    }
  }

  return 1;
}

/*
   Detect CX3 or t-butyl
*/

static int
is_cf3_like (const Molecule & m, atom_number_t zatom)
{
  const Atom * c = m.atomi(zatom);

  if (4 != c->ncon())
    return 0;

  int halogens_attached = 0;
  int ch3_attached = 0;

  for (int i = 0; i < 4; i++)
  {
    atom_number_t j = c->other(zatom, i);

    if (m.is_halogen(j))
      halogens_attached++;
    else if (1 == m.ncon(j) && 6 == m.atomic_number(j))
      ch3_attached++;
  }

  if (3 == halogens_attached)
    return 1;

  return 3 == ch3_attached;
}


/*
   Identify bonds to a CF3 or t-Butyl group
   We must be the bond leading to the CF3, not within the CF3
*/

int
is_cf3_like (const Molecule & m, const Bond & b) {
  atom_number_t a1 = b.a1();
  atom_number_t a2 = b.a2();

  if (1 == m.ncon(a1) || 1 == m.ncon(a2))
    return 0;

  if (4 == m.ncon(a1) && is_cf3_like(m, a1))
    return 1;

  if (4 == m.ncon(a2) && is_cf3_like(m, a2))
    return 1;

  return 0;
}

int
is_amide(const Molecule & m,
         atom_number_t a1,
         atom_number_t a2) {
  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  int z1 = aa1->atomic_number();
  int z2 = aa2->atomic_number();

  // At least one must be a heteroatom

  if (6 != z1)
    ;
  else if (6 != z2)
    ;
  else
    return 0;

  // And at least one must be unsaturated

  if (aa1->nbonds() > aa1->ncon())
    ;
  else if (aa2->nbonds() > aa2->ncon())
    ;
  else
    return 0;

  return 1;
}

}  // namespace dicer_lib
