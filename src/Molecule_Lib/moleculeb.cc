/*
  Molecule functions for when there is a MOLECULE_BOND_LIST
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"

#define COMPILING_CTB
#define COMPILING_MOLECULE_MAIN

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "chiral_centre.h"
#include "mdl.h"
#include "misc2.h"
#include "molecule.h"
#include "path.h"

using std::cerr;
using std::endl;

static int copy_atom_based_user_specified_void_pointers_during_add_molecule = 0;

void
set_copy_atom_based_user_specified_void_pointers_during_add_molecule(int s)
{
  copy_atom_based_user_specified_void_pointers_during_add_molecule = s;
}

/*int
Molecule::nbonds() const
{
  assert (ok());

  return _bond_list.number_elements();
}*/

int
Molecule::nedges() const
{
  assert(ok());

  return _bond_list.number_elements();
}

const Bond *
Molecule::bondi(int i) const
{
  assert(ok());
  assert(i >= 0);

  return _bond_list[i];
}

/*
  Return the J'th bond to atom I
*/

const Bond *
Molecule::bondi(atom_number_t i, int j) const
{
  assert(ok_atom_number(i));

  const Atom * a = _things[i];
  assert(a->ok_index(j));

  return a->item(j);
}

int
Molecule::saturated(atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->fully_saturated();
}

int
Molecule::unsaturation(atom_number_t zatom) const {
  assert(ok_atom_number(zatom));

  return _things[zatom]->nbonds() - _things[zatom]->ncon();
}

/*
  Return the atom number of the J'th connection to atom I in
  molecule THIS
*/

atom_number_t
Molecule::other(atom_number_t i, int j) const
{
  assert(ok_atom_number(i));

  return _things[i]->other(i, j);
}

/*
  Returns the bond type of the J'th bond to atom I in THIS molecule
*/

bond_type_t
Molecule::btype_to_connection(atom_number_t i, int j) const
{
  assert(ok_atom_number(i));
  assert(j >= 0 && j < _number_elements);

  const Atom * a1 = _things[i];
  return a1->btype_to_connection(j);
}

bond_type_t
Molecule::btype_between_atoms(atom_number_t myAtomId, atom_number_t otherAtomId) const
{
  assert(ok_2_atoms(myAtomId, otherAtomId));

  return _things[myAtomId]->btype_to_atom(myAtomId, otherAtomId);
}

int
Molecule::other_and_type(atom_number_t i, int j, atom_number_t & a, bond_type_t & bt) const
{
  assert(ok_atom_number(i));

  return _things[i]->other_and_type(i, j, a, bt);
}

/*
  Return all connections to atom A
*/

int
Molecule::connections(atom_number_t a, atom_number_t * others, bond_type_t * bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections(a, others, bt);
}

int
Molecule::connections(atom_number_t a, Set_of_Atoms & others) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections(a, others);
}

Set_of_Atoms
Molecule::connections(atom_number_t a) const {
  assert(ok_atom_number(a));
  return _things[a]->connections(a);
}

int
Molecule::connections_and_types(atom_number_t a, Set_of_Atoms & others,
                                resizable_array<bond_type_t> & bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections_and_types(a, others, bt);
}

int
Molecule::connections_and_types(atom_number_t a, atom_number_t * others, bond_type_t * bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections_and_types(a, others, bt);
}

int
Molecule::bond_types(atom_number_t a, resizable_array<bond_type_t> & bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->bond_types(bt);
}

int
Molecule::bond_types(atom_number_t a, bond_type_t * bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->bond_types(bt);
}

const Bond *
Molecule::bond_between_atoms(atom_number_t a1, atom_number_t a2) const
{

  assert(ok_2_atoms(a1, a2));

  const Bond * rc = _things[a1]->bond_to_atom(a2);

  if (nullptr != rc)
    return rc;

  cerr << "Molecule::bond_between_atoms: atoms " << a1 << " and " << a2 << " not bonded\n";
  assert(0);
  return nullptr;
}

/*
  A version that is designed to also function as a query as to whether or not
  the two atoms are bonded
*/

const Bond *
Molecule::bond_between_atoms_if_present(const atom_number_t a1, const atom_number_t a2) const
{
  if (a1 == a2)
    return nullptr;

  assert(ok_2_atoms(a1, a2));

  return _things[a1]->bond_to_atom(a2);
}

bool
Molecule::bond_endpoints(int e, atom_number_t & a1, atom_number_t & a2) const
{
  if (e > nedges())
    return false;

  const Bond * b = _bond_list[e];

  a1 = b->a1();
  a2 = b->a2();

  assert(ok_2_atoms(a1, a2));

  return true;
}

/*
  This function returns true if atoms a1 and a2 (in molecule m)
  are bonded to each other.
*/

int
Molecule::are_bonded(atom_number_t a1, atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  return _things[a1]->is_bonded_to(a2);
}

int
Molecule::are_bonded(atom_number_t a1, atom_number_t a2, bond_type_t & bt) const
{
  assert(ok_2_atoms(a1, a2));

  return _things[a1]->is_bonded_to(a2, bt);
}

int
Molecule::are_adjacent(atom_number_t a1, atom_number_t a2) const {
  return _things[a1]->connections(a1).contains(a2);
}

/*
  Return the bond number of the bond between atoms A1 and A2
*/

int
Molecule::which_bond(atom_number_t a1, atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  return _bond_list.which_bond(a1, a2);
}

//#define DEBUG_ADD_MOLECULE

Molecule&
Molecule::operator+=(const Molecule& rhs) {
  this->add_molecule(&rhs);
  return *this;
}

Molecule
Molecule::operator+(const Molecule& rhs) const {
  Molecule result(*this);
  result.add_molecule(&rhs);
  return result;
}

/*
  Function to add all the atoms of molecule m2 to molecule m1
  There is no attempt to actually join the fragments, we
  just copy over connectivity, geometry and perhaps charges.
*/

int
Molecule::add_molecule(const Molecule * m2)
{
  assert(ok());

  const int m2atoms = m2->natoms();    // will call OK_MOLECULE (m2)

#ifdef DEBUG_ADD_MOLECULE
  cerr << "add_molecule : add " << m2atoms << " atoms to " << _number_elements << endl;
  cerr << _molecule_name << endl << m2->molecule_name() << endl;
#endif

  if (m2atoms <= 0)    // wierd
    return 1;

  int m1atoms = _number_elements;

  if (_elements_allocated < m1atoms + m2atoms)
    resize(m1atoms + m2atoms);

  int m2_has_charges = m2->has_charges();
  if (m2_has_charges)
    (void)allocate_charges();

  for (int i = 0; i < m2atoms; i++)
  {
    const Atom * a2 = m2->_things[i];

   // Use copy constructor to preserve as many properties as possible, including coordinates
    Atom * tmp = new Atom(a2);

    add(tmp, 1);    // extra arg means partial molecule

    if (m2_has_charges)
      set_charge(_number_elements - 1, m2->charge_on_atom(i));

    //  cerr << "built atom " << i << ' ';
    //  tmp->debug_print(cerr);
    //  a2->debug_print(cerr);

    //  if (copy_atom_based_user_specified_void_pointers_during_add_molecule)   // redundant, the Atom copy constructor does this always
    //    tmp->set_user_specified_void_ptr(const_cast<void *>(a2->user_specified_void_ptr()));
  }

  const Bond_list & m2blist = m2->_bond_list;
  int m2bonds = m2blist.nbonds();
  _bond_list.resize(_bond_list.number_elements() + m2bonds);

  for (int i = 0; i < m2bonds; i++)
  {
    const Bond * oldb = m2blist[i];
    atom_number_t a1 = m1atoms + oldb->a1();
    atom_number_t a2 = m1atoms + oldb->a2();
    bond_type_t btype = oldb->btype();

    add_bond(a1, a2, btype, 1);    // last argument means partial molecule, avoid expensive checking

    //  Bond * b = const_cast<Bond *> (bond_between_atoms (a1, a2));   // loss of const OK
    Bond * b = _bond_list.last_item();

    b->copy_directionality_specifications(oldb);

    if (oldb->is_permanent_aromatic())
      b->set_permanent_aromatic(1);
  }

  // Copy any chiral centres

  int nc = m2->_chiral_centres.number_elements();
  for (int i = 0; i < nc; i++)
  {
    const Chiral_Centre * c2 = m2->_chiral_centres[i];
    assert(c2->ok());

    atom_number_t a2 = c2->a();

    Chiral_Centre * c1 = new Chiral_Centre(m1atoms + a2);
    c1->make_copy(c2);

    _chiral_centres.add(c1);
  }

  _set_modified();

  if (! _check_chiral_centres())
  {
    cerr << "add_molecule: invalid chiral centres\n";
    return 0;
  }

  if (! ok())
  {
    cerr << "Molecule::add_molecule failed\n";
    debug_print(cerr);
    return 0;
  }

  for (int i = 0; i < m2->_text_info.number_elements(); i++)
  {
    const IWString * inf = m2->_text_info[i];

    IWString * n = new IWString(*inf);
    _text_info.add(n);
  }

  return 1;
}

/*
  During growth, we often need to add one molecule to another
  all except for a given atom - typically an atom being eliminated
  during a condensation
  except is an atom number within molecule m2.
*/

int
Molecule::add_molecule_without(const Molecule * m2, atom_number_t except)
{
  assert(ok());
  assert(OK_ATOM_NUMBER(m2, except));

  int m1atoms = _number_elements;
  int m2atoms = m2->_number_elements;
  resize(m1atoms + m2atoms - 1);

  // First add atoms from m2 to THIS

  for (int i = 0; i < m2atoms; i++)
    if (i != except)
    {
      Atom * a = m2->_things[i];

      Atom * tmp = new Atom(a);    // zero connected

      formal_charge_t fc2 = a->formal_charge();
      if (fc2)
        tmp->set_formal_charge(fc2);

      Molecule::add(tmp);
      atom_number_t j = ADJUSTED_ATOM_NUMBER_1(i, except);
      set_charge(j, m2->charge_on_atom(i));
    }

  // Now the bonds. Be careful about atom numbering, as EXCEPT is gone

  int m2bonds = m2->nedges();
  for (int i = 0; i < m2bonds; i++)
  {
    const Bond * m2bond = m2->bondi(i);
    atom_number_t a1 = m2bond->a1();
    atom_number_t a2 = m2bond->a2();

    if ((except != a1) && (except != a2))
      add_bond(m1atoms + ADJUSTED_ATOM_NUMBER_1(a1, except),
               m1atoms + ADJUSTED_ATOM_NUMBER_1(a2, except), m2bond->btype());
  }

  _set_modified();
  check_bonding();

  return 1;
}

int
Molecule::add_molecule_without(const Molecule * m2, atom_number_t except1, atom_number_t except2)
{
  assert(OK_2_ATOMS(m2, except1, except2));

  int m1atoms = natoms();
  int m2atoms = m2->natoms();
  resize(m1atoms + m2atoms - 2);    // 2 atoms lost this time

  // First add atoms from m2 to THIS

  for (int i = 0; i < m2atoms; i++)
    if (i != except1 && i != except2)
    {
      Atom * a = m2->_things[i];
      Atom * tmp = new Atom(a);    // zero connected

      formal_charge_t fc2 = a->formal_charge();
      if (fc2)
        tmp->set_formal_charge(fc2);

      Molecule::add(tmp);
      atom_number_t j = ADJUSTED_ATOM_NUMBER_2(i, except1, except2);
      set_charge(j, m2->charge_on_atom(i));
    }

  // Now the bonds. Be careful about atom numbering, as EXCEPT is gone

  int m2bonds = m2->nedges();
  for (int i = 0; i < m2bonds; i++)
  {
    const Bond * m2bond = m2->bondi(i);
    atom_number_t a1 = m2bond->a1();
    atom_number_t a2 = m2bond->a2();

    if ((except1 != a1) && (except1 != a2) && (except2 != a1) && (except2 != a2))
      add_bond(m1atoms + ADJUSTED_ATOM_NUMBER_2(a1, except1, except2),
               m1atoms + ADJUSTED_ATOM_NUMBER_2(a2, except1, except2), m2bond->btype());
  }

  _set_modified();
  check_bonding();

  return 1;
}

/*
  When add_bond is called with the flag to indicate a partially built molecule,
  the _ncon field for each atom is not updated (this is because bonds to atoms
  not yet present in the molecule may be added).
*/

int
Molecule::finished_bond_addition()
{
  assert(ok());

  if (! _partially_built)
  {
    cerr << "finished_bond_addition: molecule already complete\n";
    return 1;
  }

  _partially_built = 0;

  if (0 == _number_elements)
    return 1;

  return 1;
}

int
Molecule::remove_bond(int which_bond)
{
  const Bond * b = _bond_list[which_bond];

  return remove_bond_between_atoms(b->a1(), b->a2());

  return 1;
}

int
Molecule::remove_bond_between_atoms(atom_number_t a1, atom_number_t a2)
{
  assert(ok_2_atoms(a1, a2));

  _things[a1]->remove_bonds_to_atom(a2);
  _things[a2]->remove_bonds_to_atom(a1);

  _bond_list.remove_bond_between_atoms(a1, a2);

  _set_modified();

  // If either of these two atoms were involved in specifying a chiral atom,
  // we are in trouble.

  int nc = _chiral_centres.number_elements();

  if (0 == nc)    // great, none to worry about
    return 1;

  _check_chirality_after_loss_of_bond(a1, a2);

  return 1;
}

/*
  We have removed a bond from the molecule.
  Do any of the chiral centres need to be updated
*/

int
Molecule::_check_chirality_after_loss_of_bond(const atom_number_t a1, const atom_number_t a2)
{
  for (int i = _chiral_centres.number_elements() - 1; i >= 0; i--)
  {
    Chiral_Centre * c = _chiral_centres[i];

    atom_number_t ca = c->a();

    if (ca == a1 && c->involves(a2))
      ;
    else if (ca == a2 && c->involves(a1))
      ;
    else
      continue;

    if (_things[ca]->ncon() < 2)    // can happen with chiral C[P@](=O)(-O)-O...
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    //  cerr << "Chiral centre at " << c->a() << " has " << implicit_hydrogens(c->a()) << " implicit hydrogens\n";
    //  c->debug_print(cerr);

    atomic_number_t z = _things[ca]->atomic_number();

    //  N1C(=S)N(C)C2=C(C1=O)C=CC=C2 PBCHM915229

    // if ((7 == z) && 2 == _things[ca]->ncon()) change jul 2023
    if ((z == 7 || z == 6) && 2 == _things[ca]->ncon())
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    // after breaking the bond, centre now has an implicit hydrogen
    if (1 == implicit_hydrogens(c->a())) {
      if (ca == a1)
        c->atom_is_now_implicit_hydrogen(a2);
      else if (ca == a2)
        c->atom_is_now_implicit_hydrogen(a1);
      else
      {
        cerr << "Molecule::remove_bond_between_atoms: huh IH " << a1 << " and " << a2 << endl;
        cerr << "Chiral Centre ";
        c->debug_print(cerr);
      }

      continue;
    }

    //  Ran into problems with things like this. Just skip over unsaturated N and P

    //  N1=C(N(C)C=C1)N(=O)=O 14835 PBCHM15477

    // Also problematic

    // C1=CC=CC=C1[S@@](=O)(=N)C PBCHM25036286

    if (((7 == z) || (15 == z)) && _things[ca]->ncon() < _things[ca]->nbonds())
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    //  cerr << "Atomic number " << z << " ncon " << _things[ca]->ncon() << endl;

    if (16 == z && _things[ca]->ncon() < 4)
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    int lp;

    // after breaking the bond, the centre atom has a lone pair
    if (0 == implicit_hydrogens(ca) && lone_pair_count(ca, lp) && 1 == lp)
    {
      if (ca == a1)
        c->atom_is_now_lone_pair(a2);
      else if (ca == a2)
        c->atom_is_now_lone_pair(a1);
      else
      {
        cerr << "Molecule::remove_bond_between_atoms: huh LP " << a1 << " and " << a2 << endl;
        cerr << "Chiral Centre ";
        c->debug_print(cerr);
      }

      if (! c->ok())
      {
        cerr << "Messed up chiral centre during removal\n";
        c->debug_print(cerr);
      }
      assert(c->ok());
    }
    else
      _chiral_centres.remove_item(i);
  }

  // Should do something about cis-trans bonds too

  return 1;
}

//#define DEBUG_REMOVE_BONDS_TO_ATOM

/*
  Function used by remove_atom to get rid of all the bonds going
  to the atom being removed
*/

#ifdef DEBUG_REMOVE_BONDS_TO_ATOM
extern int
BAD_BONDS(const Molecule &);
#endif

int
Molecule::_remove_bonds_to_atom(atom_number_t zatom, int adjust_atom_numbers)
{
  assert(ok_atom_number(zatom));

#ifdef DEBUG_REMOVE_BONDS_TO_ATOM
  cerr << "Removing bonds to atom " << zatom << " adjust = " << adjust_atom_numbers << endl;
#endif

  Atom * a = _things[zatom];

  const int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    atom_number_t j = b->other(zatom);

    _things[j]->remove_bonds_to_atom(zatom);
  }

  for (int i = _bond_list.number_elements() - 1; i >= 0; i--)
  {
    Bond * b = _bond_list[i];
    if (b->involves(zatom))
      _bond_list.remove_item(i);
    else if (adjust_atom_numbers)
      b->adjust_for_loss_of_atom(zatom);
  }

#ifdef DEBUG_PLACE_WEDGE_BOND
  cerr << "Molecule::_remove_bonds_to_atom: atom " << zatom << endl;
  BAD_BONDS(*this);
#endif

  a->resize(0);    // remove all the bonds

  a->set_modified();

  _set_modified();

  return 1;
}

int
Molecule::remove_bonds_to_atom(atom_number_t zatom, int preserve_chirality)
{
  _atom_being_unbonded_check_directional_bonds(zatom, preserve_chirality);

  if (! _remove_bonds_to_atom(zatom))
    return 0;

  for (int i = _chiral_centres.number_elements() - 1; i >= 0; i--)
  {
    Chiral_Centre * c = _chiral_centres[i];

    if (! c->involves(zatom))
      continue;

    if (! preserve_chirality)
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    if (c->implicit_hydrogen_count())    // already got an implicit Hydrogen, cannot have another one
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    c->atom_is_now_implicit_hydrogen(zatom);
  }

  return 1;
}

int
Molecule::remove_bonds_involving_these_atoms(const int * rm, int check_chirality)
{
  for (int i = 0; i < _number_elements; ++i)
  {
    if (0 == rm[i])
      continue;

    Atom * a = _things[i];    // all bonds to atom I will be removed

    for (int j = a->ncon() - 1; j >= 0; --j)
    {
      const atom_number_t k = a->other(i, j);

      if (! rm[k])    // all bonds to atom K is not being removed, so we need to carefully remove the bond to atom I
        _things[k]->remove_bonds_to_atom(i, 0);    // last arg means do NOT adjust atom numbers

      if (check_chirality)
        _check_chirality_after_loss_of_bond(i, k);
    }

    a->resize(0);
    a->set_modified();
  }

  _bond_list.remove_bonds_involving_these_atoms(rm);

  _set_modified();    // always

  return 1;
}

int
Molecule::set_bond_type_between_atoms(atom_number_t a1, atom_number_t a2, bond_type_t bt)
{
  assert(ok_2_atoms(a1, a2));
  assert(OK_BOND_TYPE(bt));

  Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));    // loss of const OK
  if (nullptr == b)
  {
    cerr << "Molecule::set_bond_type_between_atoms: atoms " << a1 << " and " << a2
         << ", no bond found\n";
    abort();
  }

  // was directional, now no longer a single bond
  if (0 == (bt & SINGLE_BOND) && b->is_directional())
    b->set_not_directional();

  // We don't check whether any cis-trans bond groupings should be broken up...

  // was part of cis/trans, now cannot be
  if (0 == (bt & DOUBLE_BOND) && b->part_of_cis_trans_grouping())
    b->set_part_of_cis_trans_grouping(0);

  b->set_bond_type(bt);

  // Perhaps we don't need to discard everything. Invalidate smiles only? Aromaticity?
  _set_modified();

  _things[a1]->set_modified();
  _things[a2]->set_modified();

  // Jan 99. If we have just placed a double bond, we need to remove any
  // chiral centres which may have been at either end.

  if (! b->is_single_bond())
  {
    if (nullptr != chiral_centre_at_atom(a1))
      remove_chiral_centre_at_atom(a1);

    if (nullptr != chiral_centre_at_atom(a2))
      remove_chiral_centre_at_atom(a2);
  }

  // Check to see if any cis-trans bonds need to be adjusted

  // double bond in middle of cis trans group was made into a single bond
  if (! b->is_double_bond() && b->part_of_cis_trans_grouping())
  {
    b->set_part_of_cis_trans_grouping(0);
  }

  return 1;
}

int
Molecule::set_all_bonds_to_type(bond_type_t bt)
{
  int rc = _bond_list.set_all_bond_types(bt);

  if (rc)
    invalidate_smiles();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);

  for (int i = 0; i < _sssr_rings.number_elements(); i++)
  {
    _sssr_rings[i]->invalidate_aromaticity();
  }

  return rc;
}

int
Molecule::set_wedge_bond_between_atoms(atom_number_t a1, atom_number_t a2, int dir)
{
  assert(ok_2_atoms(a1, a2));

  Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));    // loss of const OK
  if (nullptr == b)
  {
    cerr << "Molecule::set_bond_type_between_atoms: atoms " << a1 << " and " << a2
         << ", no bond found\n";
    abort();
  }

  if (! b->is_single_bond())
  {
    cerr << "Molecule::set_wedge_bond_between_atoms: only single bonds can be directional\n";
    cerr << "Bond between atoms " << a1 << " and " << a2 << " is not a single bond\n";
    return 0;
  }

#ifdef DEBUG_PLACE_WEDGE_BOND
  cerr << "SWBBA: Processing bond " << b->a1() << " to " << b->a2() << endl;
#endif

  if (0 == dir)
  {
    b->set_not_directional();

    return 1;
  }

  // Gets a little more complicated here, because we don't know whether a1 corresponds to a1 in the bond

#ifdef DEBUG_PLACE_WEDGE_BOND
  if (a1 != b->a1())    // reverse the direction
    cerr << "Direction reversed, now " << dir << endl;
#endif

  if (a1 != b->a1())    // reverse the direction
    dir = -dir;

  if (dir > 0)
    b->set_wedge_up();
  else
    b->set_wedge_down();

  // Don't bother with set modified, this doesn't change anything

  return 1;
}

int
Molecule::all_atoms_connected(atom_number_t zatom, Set_of_Atoms & connected_atoms) const
{
  assert(ok_atom_number(zatom));
  connected_atoms.add(zatom);

  const Atom * a = _things[zatom];

  int rc = 0;
  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);
    if (connected_atoms.contains(j))
      continue;

    rc++;
    rc += all_atoms_connected(j, connected_atoms);
  }

  return rc;
}

int
Molecule::chop(int items_removed)
{
  assert(ok());

  if (items_removed <= 0) {
    return 0;
  }

  if (items_removed > _number_elements) {
    cerr << "Molecule::chop:natoms " << _number_elements << " cannot remove " << items_removed << '\n';
    return 0;
  }

  return resize(_number_elements - items_removed);
}

int
Molecule::number_up_or_down_wedge_bonds() const
{
  int rc = 0;

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    if (_bond_list[i]->is_wedge_definitive())
      rc++;
  }

  return rc;
}

int
Molecule::assign_bond_numbers_to_bonds()
{
  return _bond_list.assign_bond_numbers(0);
}

int
Molecule::assign_bond_numbers_to_bonds_if_needed()
{
  return _bond_list.assign_bond_numbers_to_bonds_if_needed();
}

int
Bond_list::assign_bond_numbers(int istart)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_bond_number(i + istart);
  }

  return 1;
}

int
Bond_list::assign_bond_numbers_to_bonds_if_needed()
{
  if (0 == _number_elements)
    return 0;

  // If first and last items are assigned, assume OK

  if (_things[0]->bond_number_assigned() && _things[_number_elements - 1]->bond_number_assigned())
    return 0;

  return assign_bond_numbers(0);
}

int
Molecule::unset_all_permanent_aromatic_bonds()
{
  return _bond_list.unset_all_permanent_aromatic_bonds();
}

int
Bond_list::unset_all_permanent_aromatic_bonds()
{
  assert(ok());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_permanent_aromatic(0);
  }

  return 1;
}

atom_number_t
Molecule::atom_with_isotope(const isotope_t iso) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (iso == _things[i]->isotope())
      return i;
  }

  return INVALID_ATOM_NUMBER;
}

int
Molecule::remove_all_bonds()
{
  _set_modified();

  _bond_list.resize(0);

  _chiral_centres.resize(0);

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->resize(0);
    _things[i]->set_modified();
  }

  return 1;
}

int
Molecule::remove_all_bonds_keep_storage()
{
  _set_modified();

  _bond_list.resize_keep_storage(0);

  _chiral_centres.resize_keep_storage(0);

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->resize_keep_storage(0);
    _things[i]->set_modified();
  }

  return 1;
}

int
Molecule::get_bond_types(bond_type_t * bt) const
{
  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    bt[i] = BOND_TYPE_ONLY_MASK & _bond_list[i]->btype();
  }

  return 1;
}

int
Molecule::set_bond_types_no_set_modified(const bond_type_t * bt)
{
  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    _bond_list[i]->set_bond_type(bt[i]);
  }

  return 1;
}

static uint64_t
do_or(uint64_t x, int bits_left, const uint64_t maxval, const uint64_t rc) {
  if (x > maxval) {
    x = maxval;
  }

  return rc | (x << bits_left);
}

uint64_t
Molecule::quick_bond_hash() {
  if (_bond_list.empty()) {
    return 0;
  }

  compute_aromaticity_if_needed();

  uint64_t c_single_c = 0;
  uint64_t c_double_c = 0;
  uint64_t c_triple_c = 0;    // 2 bits
  uint64_t c_arom_c = 0;

  uint64_t c_single_n = 0;
  uint64_t c_double_n = 0;
  uint64_t c_triple_n = 0;    // 3 bits
  uint64_t c_arom_n = 0;

  uint64_t c_single_o = 0;
  uint64_t c_double_o = 0;
  uint64_t c_arom_o = 0;    // 3 bits

  uint64_t c_single_s = 0;
  uint64_t c_double_s = 0;    // 2 bits
  uint64_t c_arom_s = 0;      // 3 bits

  uint64_t n_arom_n = 0;      // 2 bits
  uint64_t n_single_n = 0;    // bug, we are overlooking this
  uint64_t n_double_n = 0;    // bug, we are overlooking this
  uint64_t n_triple_n = 0;    // bug, we are overlooking this

  uint64_t n_single_o = 0;    // bug, overlooked
  uint64_t n_double_o = 0;
  uint64_t n_arom_o = 0;    // bug, overlooked

  uint64_t o_single_s = 0;
  uint64_t o_double_s = 0;

  uint64_t aliphatic_ring_bond_count = 0;

  const int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = _bond_list[i];

    if (b->nrings() && ! b->is_aromatic())
      aliphatic_ring_bond_count++;

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    atomic_number_t z1 = _things[a1]->atomic_number();

    atomic_number_t z2 = _things[a2]->atomic_number();

    if (z1 > z2)
      std::swap(z1, z2);

    if (6 == z1)
    {
      if (6 == z2)
      {
        if (b->is_aromatic())
          c_arom_c++;
        else if (b->is_single_bond())
          c_single_c++;
        else if (b->is_double_bond())
          c_double_c++;
        else if (b->is_triple_bond())
          c_triple_c++;
      }
      else if (7 == z2)
      {
        if (b->is_aromatic())
          c_arom_n++;
        else if (b->is_single_bond())
          c_single_n++;
        else if (b->is_double_bond())
          c_double_n++;
        else if (b->is_triple_bond())
          c_triple_n++;
      }
      else if (8 == z2)
      {
        if (b->is_aromatic())
          c_arom_o++;
        else if (b->is_single_bond())
          c_single_o++;
        else if (b->is_double_bond())
          c_double_o++;
      }
      else if (16 == z2)
      {
        if (b->is_aromatic())
          c_arom_s++;
        else if (b->is_single_bond())
          c_single_s++;
        else if (b->is_double_bond())
          c_double_s++;
      }
    }
    else if (7 == z1)
    {
      if (7 == z2)
      {
        if (b->is_aromatic())
          n_arom_n++;
        else if (b->is_single_bond())
          n_single_n++;
        else if (b->is_double_bond())
          n_double_n++;
        else if (b->is_triple_bond())
          n_triple_n++;
      }
      else if (8 == z2)
      {
        if (b->is_aromatic())
          n_arom_o++;
        else if (b->is_single_bond())
          n_single_o++;
        else if (b->is_double_bond())
          n_double_o++;
      }
    }
    else if (8 == z1)
    {
      if (16 == z2)
      {
        if (b->is_single_bond())
          o_single_s++;
        else if (b->is_double_bond())
          o_double_s++;
      }
    }
  }

  uint64_t rc = 0;

  // clang-format off
  int bshift = 64-4;

  rc = do_or(c_single_c, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_double_c, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_triple_c, bshift, 3, rc); bshift -= 2;
  rc = do_or(c_arom_c/4, bshift, 15, rc); bshift -= 4;

  rc = do_or(c_single_n, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_double_n, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_triple_n, bshift, 7, rc); bshift -= 3;
  rc = do_or(c_arom_n,   bshift, 15, rc); bshift -= 4;

  rc = do_or(c_single_o, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_double_o, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_arom_o,   bshift, 1, rc); bshift -= 1;

  rc = do_or(c_single_s, bshift, 15, rc); bshift -= 4;
  rc = do_or(c_double_s, bshift, 3, rc); bshift -= 2;
  rc = do_or(c_arom_s,   bshift, 7, rc); bshift -= 3;

//rc = do_or(n_double_o, bshift, 15, rc); bshift -= 4

  rc = do_or(o_single_s, bshift, 15, rc); bshift -= 4;
  rc = do_or(o_double_s, bshift, 15, rc); bshift -= 4;
  // clang-format on

  rc = do_or(aliphatic_ring_bond_count / 2, bshift, 15, rc);
  bshift -= 4;

  assert(bshift >= 0);
  //if (bshift < 0)
  //  cerr << "negative bshift\n";

  return rc;
}

// Shift `destination` left by nbits and OR `value`.
// If `value` is too large for that many bits, truncate and increment `overflows`.
static void
ShiftOr(uint64_t value, int nbits, uint64_t &overflows, uint64_t& destination) {

  static uint64_t maxval[] = {0, 1, 7, 16, 31, 63, 127, 255, 511, 1023};

  if (value > maxval[nbits]) {
    value = maxval[nbits];
    ++overflows;
  }

  destination = (destination << nbits) | value;
}

void
Molecule::QuickBondHash(uint64_t hash[2]) {
  hash[0] = 0;
  hash[1] = 0;

  if (_bond_list.empty()) {
    return;
  }

  compute_aromaticity_if_needed();

  uint64_t c_single_c = 0;
  uint64_t c_single_ring_c = 0;
  uint64_t c_double_c = 0;
  uint64_t c_triple_c = 0;
  uint64_t c_arom_c = 0;
  uint64_t carom_single_c = 0;
  uint64_t carom_single_n = 0;
  uint64_t carom_single_o = 0;
  uint64_t biphenyl = 0;
  uint64_t carom_double_n = 0;
  uint64_t carom_double_o = 0;
  uint64_t carom_ring_single_n = 0;
  uint64_t carom_ring_single_o = 0;

  uint64_t c_single_n = 0;
  uint64_t c_double_n = 0;
  uint64_t c_triple_n = 0;
  uint64_t c_arom_n = 0;

  uint64_t c_single_o = 0;
  uint64_t c_double_o = 0;
  uint64_t c_arom_o = 0;

  uint64_t c_single_s = 0;
  uint64_t c_double_s = 0;
  uint64_t c_arom_s = 0;

  uint64_t carom_f = 0;
  uint64_t c_single_f = 0;
  uint64_t c_single_halogen = 0;
  uint64_t c_other = 0;

  uint64_t n_arom_n = 0;
  uint64_t n_single_n = 0;
  uint64_t n_double_n = 0;
  uint64_t n_triple_n = 0;

  uint64_t n_single_o = 0;
  uint64_t n_double_o = 0;
  uint64_t n_arom_o = 0;
  uint64_t n_arom_s = 0;
  uint64_t narom_single_c = 0;

  uint64_t o_single_o = 0;
  uint64_t o_single_s = 0;
  uint64_t o_double_s = 0;

  uint64_t other_bond = 0;

  uint64_t aliphatic_ring_bond_count = 0;

  // The number of times something overflows the maximum
  uint64_t overflow_count = 0;

  for (const Bond* b : _bond_list) {
    if (b->nrings() && ! b->is_aromatic())
      aliphatic_ring_bond_count++;

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    atomic_number_t z1 = _things[a1]->atomic_number();

    atomic_number_t z2 = _things[a2]->atomic_number();

    int arom1 = is_aromatic(a1);
    int arom2 = is_aromatic(a2);

    if (z1 > z2) {
      std::swap(z1, z2);
      std::swap(arom1, arom2);
    }

    if (6 == z1) {
      if (6 == z2)
      {
        if (b->is_aromatic()) {
          c_arom_c++;
        } else if (b->is_single_bond()) {
          const uint64_t arom = arom1 + arom2;
          if (arom == 0) {
            if (b->nrings()) {
              c_single_ring_c++;
            } else {
              c_single_c++;
            }
          } else if (arom == 1) {
            ++carom_single_c;
          } else {
            ++biphenyl;
          }
        }
        else if (b->is_double_bond())
          c_double_c++;
        else if (b->is_triple_bond())
          c_triple_c++;
      }
      else if (7 == z2)
      {
        if (b->is_aromatic()) {
          c_arom_n++;
        } else if (b->is_single_bond()) {
          if (arom1) {
            if (b->nrings()) {
              ++carom_ring_single_n;
            } else {
              ++carom_single_n;
            }
          } else if (arom2) {
            ++narom_single_c;
          } else {
            ++c_single_n;
          }
        }
        else if (b->is_double_bond()) {
          if (arom1 == 0) {
            c_double_n++;
          } else {
            ++carom_double_n;
          }
        }
        else if (b->is_triple_bond())
          c_triple_n++;
      }
      else if (8 == z2)
      {
        if (b->is_aromatic())
          c_arom_o++;
        else if (b->is_single_bond()) {
          if (arom1) {
            if (b->nrings()) {
              ++carom_ring_single_o;
            } else {
              ++carom_single_o;
            }
          } else {
            c_single_o++;
          }
        }
        else if (b->is_double_bond()) {
          if (arom1) {
            ++carom_double_o;
          } else {
            c_double_o++;
          }
        }
      }
      else if (16 == z2)
      {
        if (b->is_aromatic())
          c_arom_s++;
        else if (b->is_single_bond())
          c_single_s++;
        else if (b->is_double_bond())
          c_double_s++;
      }
      else if (z2 == 9) {
        if (arom1) {
          ++carom_f;
        } else {
          ++c_single_f;
        }
      } else if (z2 == 17 || z2 == 35 || z2 == 53) {
        ++c_single_halogen;
      } else {
        ++c_other;
      }
    }

    else if (7 == z1)
    {
      if (7 == z2)
      {
        if (b->is_aromatic())
          n_arom_n++;
        else if (b->is_single_bond())
          n_single_n++;
        else if (b->is_double_bond())
          n_double_n++;
        else if (b->is_triple_bond())
          n_triple_n++;
      }
      else if (8 == z2)
      {
        if (b->is_aromatic())
          n_arom_o++;
        else if (b->is_single_bond())
          n_single_o++;
        else if (b->is_double_bond())
          n_double_o++;
      }
      else if (z2 == 16) {
        if (b->is_double_bond()) {
          ++n_arom_s;
        }
      }
    }
    else if (8 == z1)
    {
      if (z2 == 8) {
        ++o_single_o;
      } else if (16 == z2) {
        if (b->is_single_bond())
          o_single_s++;
        else if (b->is_double_bond())
          o_double_s++;
      }
    } else {
      ++other_bond;
    }
  }

  // clang-format off

  // Handy table of number of bits vs max value accommodated.
  // 9 1 2 3  4  5
  // 0 1 3 7 15 31

  ShiftOr(c_single_c, 6, overflow_count, hash[0]);  // total = 6
  ShiftOr(c_double_c, 3, overflow_count, hash[0]);  // total = 9
  ShiftOr(c_triple_c, 2, overflow_count, hash[0]);  // total = 11
  ShiftOr(c_arom_c, 5, overflow_count, hash[0]); // total = 16

  ShiftOr(c_single_n, 4, overflow_count, hash[0]); // total = 20
  ShiftOr(c_double_n, 3, overflow_count, hash[0]); // total = 23
  ShiftOr(c_triple_n, 2, overflow_count, hash[0]); // total = 25
  ShiftOr(c_arom_n, 4, overflow_count, hash[0]); // total = 31

  ShiftOr(c_single_o, 4, overflow_count, hash[0]); // total = 35
  ShiftOr(c_double_o, 3, overflow_count, hash[0]); // total = 38
  ShiftOr(c_arom_o, 3, overflow_count, hash[0]); // total = 41

  ShiftOr(c_single_s, 3, overflow_count, hash[0]); // total = 44
  ShiftOr(c_double_s, 2, overflow_count, hash[0]); // total = 46
  ShiftOr(c_arom_s, 3, overflow_count, hash[0]); // total = 49

  ShiftOr(n_arom_n, 3, overflow_count, hash[0]); // total = 52
  ShiftOr(n_single_n, 2, overflow_count, hash[0]); // total = 54
  ShiftOr(n_double_n, 1, overflow_count, hash[0]); // total = 55
  ShiftOr(n_triple_n, 1, overflow_count, hash[0]); // total = 56
  ShiftOr(n_arom_o, 2, overflow_count, hash[0]); // total = 58
  ShiftOr(n_arom_s, 1, overflow_count, hash[0]); // total = 59
  ShiftOr(n_single_o, 2, overflow_count, hash[0]); // total = 61
  ShiftOr(n_double_o, 3, overflow_count, hash[0]); // total = 64

  ShiftOr(o_single_o, 1, overflow_count, hash[1]); // total = 1
  ShiftOr(o_single_s, 3, overflow_count, hash[1]); // total = 4
  ShiftOr(o_double_s, 3, overflow_count, hash[1]); // total = 7

  ShiftOr(aliphatic_ring_bond_count, 5, overflow_count, hash[1]);  // total = 14

  ShiftOr(carom_single_c, 4, overflow_count, hash[1]); // total = 18
  ShiftOr(carom_single_n, 3, overflow_count, hash[1]); // total = 21
  ShiftOr(carom_single_o, 3, overflow_count, hash[1]); // total = 24
  ShiftOr(biphenyl, 2, overflow_count, hash[1]); // total = 26
  ShiftOr(carom_double_n, 2, overflow_count, hash[1]); // total = 28
  ShiftOr(carom_double_o, 2, overflow_count, hash[1]); // total = 30

  ShiftOr(carom_ring_single_o, 2, overflow_count, hash[1]); // total = 32
  ShiftOr(carom_ring_single_n, 2, overflow_count, hash[1]); // total = 34
  ShiftOr(narom_single_c, 2, overflow_count, hash[1]); // total = 36
  ShiftOr(c_single_ring_c, 4, overflow_count, hash[1]); // total = 40

  // Must be last.
  ShiftOr(c_other, 1, overflow_count, hash[1]); // total = 41
  ShiftOr(other_bond, 2, overflow_count, hash[1]); // total = 43

  // clang-format on
}

const_BondIterator::const_BondIterator(const Atom * a, const atom_number_t zatom) : _atnum(zatom)
{
  _b = a->cbegin();

  return;
}
int
Molecule::CanonicaliseBondList() {
  for (Bond * b : _bond_list) {
    if (b->a1() > b->a2()) {
      b->set_a1a2(b->a2(), b->a1());
    }
  }
  _bond_list.iwqsort_lambda([](const Bond * b1, const Bond * b2) {
    if (b1->a1() < b2->a1()) {
      return -1;
    }
    if (b1->a1() > b2->a1()) {
      return 1;
    }
    if (b1->a2() < b2->a2()) {
      return -1;
    }
    return 1;
  });

  for (int i = 0; i < _number_elements; ++i) {
    _things[i]->CanonicaliseBonds();
  }

  _set_modified();

  return 1;
}
