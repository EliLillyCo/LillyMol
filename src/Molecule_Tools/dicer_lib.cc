#include "Molecule_Tools/dicer_lib.h"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "Molecule_Lib/chiral_centre.h"

namespace dicer_lib {

namespace {
bool
HasChiralDummy(const Molecule& m, atom_number_t centre) {
  for (const Bond* bond : m[centre]) {
    const atom_number_t neighbour = bond->other(centre);
    if (m.atomic_number(neighbour) == 0 && !InitialAtomNumber(m.atomi(neighbour))) {
      return true;
    }
  }
  return false;
}

atom_number_t
AddChiralDummy(Molecule& m) {
  const atom_number_t result = m.add(new Atom(0));
  m.set_implicit_hydrogens(result, 0, 1);
  return result;
}
}  // namespace

bool
BreakBondWithChiralDummies(Molecule& m, atom_number_t a1, atom_number_t a2) {
  const Bond* bond = m.bond_between_atoms(a1, a2);
  if (bond == nullptr) {
    return false;
  }
  Chiral_Centre* c1 = m.chiral_centre_at_atom(a1);
  Chiral_Centre* c2 = m.chiral_centre_at_atom(a2);
  if (!c1 && !c2) {
    m.remove_bond_between_atoms(a1, a2);
    return true;
  }
  // Preflight both ends before adding atoms or changing connections.
  if (!bond->is_single_bond() || (c1 && HasChiralDummy(m, a1)) ||
      (c2 && HasChiralDummy(m, a2))) {
    return false;
  }
  if (!c1) {
    std::swap(a1, a2);
    std::swap(c1, c2);
  }
  const atom_number_t dummy1 = AddChiralDummy(m);
  if (!m.stereo_preserving_substitute(a1, a2, dummy1)) {
    return false;
  }
  if (c2) {
    // The substitution has already detached the shared bond at its other end.
    // Replace that end's saved neighbour slot directly, then add its own cap.
    const atom_number_t dummy2 = AddChiralDummy(m);
    // add_bond would otherwise replace an implicit-H stereo slot. Detach the
    // centre temporarily, then restore the explicitly substituted ordering.
    std::unique_ptr<Chiral_Centre> second(m.remove_no_delete_chiral_centre_at_atom(a2));
    second->change_atom_number(a1, dummy2);
    m.add_bond(a2, dummy2, SINGLE_BOND);
    m.add_chiral_centre(second.release());
  } else {
    m.unset_all_implicit_hydrogen_information(a2);
  }
  return true;
}

bool
CreateChiralDummySubsets(Molecule& parent, const int* subset, Molecule& fragment,
                         Molecule& complement, int* xref) {
  Molecule capped(parent);
  // Only capped is modified; the parent's bond list remains safe to traverse.
  for (const Bond* bond : parent.bond_list()) {
    if (subset[bond->a1()] == subset[bond->a2()]) {
      continue;
    }
    if (!BreakBondWithChiralDummies(capped, bond->a1(), bond->a2())) {
      return false;
    }
  }
  std::vector<int> membership(capped.natoms());
  std::copy_n(subset, parent.natoms(), membership.begin());
  for (int i = parent.natoms(); i < capped.natoms(); ++i) {
    membership[i] = membership[capped.other(i, 0)];
  }
  std::vector<int> full_xref(capped.natoms());
  capped.create_subsets(membership.data(), complement, fragment, full_xref.data());
  std::copy_n(full_xref.begin(), parent.natoms(), xref);
  return true;
}

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
