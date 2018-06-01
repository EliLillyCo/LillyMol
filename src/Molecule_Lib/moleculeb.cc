/*
  Molecule functions for when there is a MOLECULE_BOND_LIST
*/

#include <stdlib.h>
#include <memory>

#include "cmdline.h"
#include "misc.h"

#define COMPILING_CTB
#define COMPILING_MOLECULE_MAIN

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "molecule.h"
#include "mdl.h"
#include "chiral_centre.h"
#include "path.h"
#include "misc2.h"

static int copy_atom_based_user_specified_void_pointers_during_add_molecle = 0;

void
set_copy_atom_based_user_specified_void_pointers_during_add_molecle (int s)
{
  copy_atom_based_user_specified_void_pointers_during_add_molecle = s;
}

static int exclude_triple_bonds_from_graph_reduction = 0;

void
set_exclude_triple_bonds_from_graph_reduction(int s)
{
  exclude_triple_bonds_from_graph_reduction = s;
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
  assert (ok());

  return _bond_list.number_elements();
}

const Bond *
Molecule::bondi (int i) const
{
  assert(ok());
  assert(i >= 0);

  return _bond_list[i];
}

/*
  Return the J'th bond to atom I
*/

const Bond *
Molecule::bondi (atom_number_t i, int j) const
{
  assert(ok_atom_number(i));

  const Atom * a = _things[i];
  assert (a->ok_index(j));

  return a->item(j);
}

int
Molecule::saturated(atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->fully_saturated();
}

/*
  Return the atom number of the J'th connection to atom I in
  molecule THIS
*/

atom_number_t
Molecule::other (atom_number_t i, int j) const
{
  assert(ok_atom_number(i));

  return _things[i]->other(i, j);
}

/*
  Returns the bond type of the J'th bond to atom I in THIS molecule
*/

bond_type_t
Molecule::btype_to_connection (atom_number_t i, int j) const
{
  assert(ok_atom_number(i));
  assert(j >= 0 && j < _number_elements);

  const Atom * a1 = _things[i];
  return a1->btype_to_connection(j);
}

bond_type_t
Molecule::btype_between_atoms (atom_number_t a1, atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  return _things[a1]->btype_to_atom(a2);
}

int
Molecule::other_and_type (atom_number_t i, int j, atom_number_t & a, bond_type_t & bt) const
{
  assert(ok_atom_number(i));

  return _things[i]->other_and_type(i, j, a, bt);
}

/*
  Return all connections to atom A
*/

int
Molecule::connections (atom_number_t a, atom_number_t * others,
                       bond_type_t * bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections(a, others, bt);
}

int
Molecule::connections (atom_number_t a, Set_of_Atoms & others) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections(a, others);
}

int
Molecule::connections_and_types (atom_number_t a,
                       Set_of_Atoms & others,
                       resizable_array<bond_type_t> & bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections_and_types(a, others, bt);
}

int
Molecule::connections_and_types (atom_number_t a,
                       atom_number_t * others,
                       bond_type_t * bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->connections_and_types(a, others, bt);
}

int
Molecule::bond_types (atom_number_t a,
                      resizable_array<bond_type_t> & bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->bond_types(bt);
}

int
Molecule::bond_types (atom_number_t a,
                      bond_type_t * bt) const
{
  assert(ok_atom_number(a));

  const Atom * aa = _things[a];
  return aa->bond_types(bt);
}

const Bond *
Molecule::bond_between_atoms (atom_number_t a1, atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  const Bond * rc = _things[a1]->bond_to_atom(a2);

  if (NULL != rc)
    return rc;

  cerr << "Molecule::bond_between_atoms: atoms " << a1 << " and " << a2 << " not bonded\n";
  assert(0);
  return NULL;
}

/*
  A version that is designed to also function as a query as to whether or not
  the two atoms are bonded
*/

const Bond * 
Molecule::bond_between_atoms_if_present(const atom_number_t a1, const atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  return _things[a1]->bond_to_atom(a2);
}

bool
Molecule::bond_endpoints (int e, atom_number_t & a1, atom_number_t & a2 ) const
{
  if( e > nedges() ) 
    return false;

  const Bond * b = _bond_list[e];

  a1 = b->a1();
  a2 = b->a2();

  assert(ok_2_atoms( a1, a2 ) );

  return true;
}
  
/*
  This function returns true if atoms a1 and a2 (in molecule m)
  are bonded to each other.
*/

int
Molecule::are_bonded (atom_number_t a1, atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  return _things[a1]->is_bonded_to(a2);
}

int
Molecule::are_bonded (atom_number_t a1, atom_number_t a2,
                      bond_type_t & bt) const
{
  assert(ok_2_atoms(a1, a2));

  return _things[a1]->is_bonded_to(a2, bt);
}

/*
  Return the bond number of the bond between atoms A1 and A2
*/

int
Molecule::which_bond (atom_number_t a1, atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  return _bond_list.which_bond(a1, a2);
}

//#define DEBUG_ADD_MOLECULE

/*
  Function to add all the atoms of molecule m2 to molecule m1
  There is no attempt to actually join the fragments, we
  just copy over connectivity, geometry and perhaps charges.
*/

int
Molecule::add_molecule(const Molecule *m2)
{
  assert(ok());

  const int m2atoms = m2->natoms();   // will call OK_MOLECULE (m2)

#ifdef DEBUG_ADD_MOLECULE
  cerr << "add_molecule : add " << m2atoms << " atoms to " << _number_elements << endl;
  cerr << _molecule_name << endl << m2->molecule_name() << endl;
#endif

  if (m2atoms <= 0)       // wierd
    return 1;

  int m1atoms = _number_elements;

  if (_elements_allocated < m1atoms + m2atoms)
    resize(m1atoms + m2atoms);

  int m2_has_charges = m2->has_charges();
  if (m2_has_charges)
    (void) allocate_charges();

  for (int i = 0; i < m2atoms; i++)
  {
    const Atom * a2 = m2->_things[i];

    Atom * tmp = new Atom(a2);     // using copy constructor preserves as many properties as possible, including coordinates

    add(tmp, 1);    // extra arg means partial molecule

    if (m2_has_charges)
      set_charge(_number_elements - 1, m2->charge_on_atom(i));

//  cerr << "built atom " << i << ' ';
//  tmp->debug_print(cerr);
//  a2->debug_print(cerr);

//  if (copy_atom_based_user_specified_void_pointers_during_add_molecle)   // redundant, the Atom copy constructor does this always
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
    bond_type_t   btype = oldb->btype();

    add_bond(a1, a2, btype, 1);   // last argument means partial molecule, avoid expensive checking

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
Molecule::add_molecule_without (const Molecule *m2, atom_number_t except)
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
      Atom *a = m2->_things[i];

      Atom *tmp = new Atom(a);     // zero connected

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
    const Bond *m2bond = m2->bondi(i);
    atom_number_t a1 = m2bond->a1();
    atom_number_t a2 = m2bond->a2();

    if ((except != a1) && (except != a2) )
      add_bond(m1atoms + ADJUSTED_ATOM_NUMBER_1(a1, except),
                m1atoms + ADJUSTED_ATOM_NUMBER_1(a2, except), m2bond->btype());
  }

  _set_modified();
  check_bonding();

  return 1;
}

int
Molecule::add_molecule_without (const Molecule *m2, atom_number_t except1,
                                                    atom_number_t except2)
{
  assert(OK_2_ATOMS(m2, except1, except2));

  int m1atoms = natoms();
  int m2atoms = m2->natoms();
  resize(m1atoms + m2atoms - 2);    // 2 atoms lost this time

// First add atoms from m2 to THIS

  for (int i = 0; i < m2atoms; i++)
    if (i != except1 && i != except2)
    {
      Atom *a = m2->_things[i];
      Atom *tmp = new Atom(a);     // zero connected

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
    const Bond *m2bond = m2->bondi(i);
    atom_number_t a1 = m2bond->a1();
    atom_number_t a2 = m2bond->a2();

    if ( (except1 != a1) && (except1 != a2) &&
         (except2 != a1) && (except2 != a2) )
      add_bond (m1atoms + ADJUSTED_ATOM_NUMBER_2(a1, except1, except2),
                m1atoms + ADJUSTED_ATOM_NUMBER_2(a2, except1, except2),
                m2bond->btype());
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
Molecule::remove_bond (int which_bond)
{
  const Bond * b = _bond_list[which_bond];

  return remove_bond_between_atoms(b->a1(), b->a2());

  return 1;
}

int
Molecule::remove_bond_between_atoms (atom_number_t a1, atom_number_t a2)
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
Molecule::_check_chirality_after_loss_of_bond (const atom_number_t a1,
                                               const atom_number_t a2)
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

    if (_things[ca]->ncon() < 2)   // can happen with chiral C[P@](=O)(-O)-O...
    {
      _chiral_centres.remove_item(i);
      continue;
    }

//  cerr << "Chiral centre at " << c->a() << " has " << implicit_hydrogens(c->a()) << " implicit hydrogens\n";
//  c->debug_print(cerr);

    atomic_number_t z = _things[ca]->atomic_number();

//  N1C(=S)N(C)C2=C(C1=O)C=CC=C2 PBCHM915229

    if ((7 == z) && 2 == _things[ca]->ncon())
    {
      _chiral_centres.remove_item(i);
      continue;
    }

    if (1 == implicit_hydrogens(c->a()))   // after breaking the bond, centre now has an implicit hydrogen
    {
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

    if (0 == implicit_hydrogens(ca) && lone_pair_count(ca, lp) && 1 == lp)    // after breaking the bond, the centre atom has a lone pair
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
extern int BAD_BONDS(const Molecule &);
#endif

int
Molecule::_remove_bonds_to_atom (atom_number_t zatom, int adjust_atom_numbers)
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

  a->resize(0);     // remove all the bonds

  a->set_modified();

  _set_modified();

  return 1;
}

int
Molecule::remove_bonds_to_atom (atom_number_t zatom,
                                int preserve_chirality)
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
Molecule::remove_bonds_involving_these_atoms (const int * rm, int check_chirality)
{
  for (int i = 0; i < _number_elements; ++i)
  {
    if (0 == rm[i])
      continue;

    Atom * a = _things[i];    // all bonds to atom I will be removed

    for (int j = a->ncon() - 1; j >= 0; --j)
    {
      const atom_number_t k = a->other(i, j);

      if (! rm[k])     // all bonds to atom K is not being removed, so we need to carefully remove the bond to atom I
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
Molecule::set_bond_type_between_atoms (atom_number_t a1, atom_number_t a2,
                                       bond_type_t bt)
{
  assert(ok_2_atoms(a1, a2));
  assert(OK_BOND_TYPE(bt));

  Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));   // loss of const OK
  if (NULL == b)
  {
    cerr << "Molecule::set_bond_type_between_atoms: atoms " << a1 << " and " << a2 << ", no bond found\n";
    abort();
  }

  if (0 == (bt | SINGLE_BOND) && b->is_directional())   // was directional, now no longer a single bond
    b->set_not_directional();

// We don't check whether any cis-trans bond groupings should be broken up...

  if (0 == (bt | DOUBLE_BOND) && b->part_of_cis_trans_grouping())  // was part of cis/trans, now cannot be
    b->set_part_of_cis_trans_grouping(0);

  b->set_bond_type(bt);

  _set_modified();

  _things[a1]->set_modified();
  _things[a2]->set_modified();

// Jan 99. If we have just placed a double bond, we need to remove any
// chiral centres which may have been at either end.

  if (! b->is_single_bond())
  {
    if (NULL != chiral_centre_at_atom(a1))
      remove_chiral_centre_at_atom(a1);

    if (NULL != chiral_centre_at_atom(a2))
      remove_chiral_centre_at_atom(a2);
  }

// Check to see if any cis-trans bonds need to be adjusted

  if (! b->is_double_bond() && b->part_of_cis_trans_grouping())   // double bond in middle of cis trans group was made into a single bond
  {
    b->set_part_of_cis_trans_grouping(0);
  }

  return 1;
}

int
Molecule::set_all_bonds_to_type (bond_type_t bt)
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
Molecule::set_wedge_bond_between_atoms (atom_number_t a1, atom_number_t a2,
                                        int dir)
{
  assert (ok_2_atoms(a1, a2));

  Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));   // loss of const OK
  if (NULL == b)
  {
    cerr << "Molecule::set_bond_type_between_atoms: atoms " << a1 << " and " << a2 << ", no bond found\n";
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
  if (a1 != b->a1())     // reverse the direction
    cerr << "Direction reversed, now " << dir << endl;
#endif

  if (a1 != b->a1())     // reverse the direction
    dir = - dir;

  if (dir > 0)
    b->set_wedge_up();
  else
    b->set_wedge_down();

// Don't bother with set modified, this doesn't change anything

  return 1;
}

int
Molecule::all_atoms_connected (atom_number_t zatom,
                               Set_of_Atoms & connected_atoms) const 
{
  assert (ok_atom_number(zatom));
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
    rc+= all_atoms_connected(j, connected_atoms);
  }

  return rc;
}

int
Molecule::chop (int new_size)
{
  assert(ok());

  assert(new_size > 0);

  if (new_size == _number_elements)
    return 0;

  assert (new_size < _number_elements);

  return resize(new_size);
}

int
Molecule::change_to_graph_form()
{
  int changes = 0;   // changes which require recomputation of nbonds & implicit hydrogens

// We remove all explicit hydrogens

  Set_of_Atoms atoms_to_be_removed;

  int nb = _bond_list.number_elements();
  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (b->is_single_bond())
      continue;

    if (b->is_triple_bond() && exclude_triple_bonds_from_graph_reduction)
      continue;

    b->set_bond_type(SINGLE_BOND);
    b->set_permanent_aromatic(0);
    changes++;
  }

  if (_chiral_centres.number_elements())
  {
    _chiral_centres.resize(0);
    changes++;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (1 == a->atomic_number())
    {
      atoms_to_be_removed.add(i);
      changes++;
      continue;
    }

//  Remove charges unless it would create valence problems

    if (a->formal_charge())
    {
      if (7 == a->atomic_number() && 4 == a->ncon())
        ;
      else
      {
        a->set_formal_charge(0);
        changes++;
      }
    }

    if (a->implicit_hydrogens_known())
    {
      a->set_implicit_hydrogens_known(0);
      changes++;
    }

    if (a->is_isotope())
    {
      a->set_isotope(0);
    }
  }

  if (atoms_to_be_removed.number_elements())
    remove_atoms(atoms_to_be_removed);

  if (revert_all_directional_bonds_to_non_directional())
    changes++;

  if (changes)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->recompute_nbonds();
      int unused;
      _things[i]->recompute_implicit_hydrogens(unused);
    }
  }

  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);

  return 1;
}

/*
  Are all the connections from atom ZATOM fully saturated
*/

int
Molecule::_all_connections_saturated (const atom_number_t zatom,
                                      const atom_number_t ignore) const
{
  const Atom * a = _things[zatom];

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);

    if (j == ignore)
      continue;

    const Atom * aj = _things[j];

    if (aj->ncon() < aj->nbonds())    // unsaturated
      return 0;
  }

  return 1;
}

int
Molecule::_double_bond_needs_changing_for_graph_form (const Bond & b,
                                               const Mol2Graph & mol2graph) const
{
  assert (b.is_double_bond());

  if (! mol2graph.some_kind_of_double_bond_preservation_active())
    return 1;

// Now we have non-aromatic double bonds

  const Atom * a1 = _things[b.a1()];
  if (6 != a1->atomic_number())    // must be changed
    return 1;

  const Atom * a2 = _things[b.a2()];
  if (6 != a2->atomic_number())    // must be changed
    return 1;

  if (mol2graph.preserve_cc_double_bonds_saturated())
  {
    if (! _all_connections_saturated(b.a1(), b.a2()))    // must be changed
      return 1;

    if (! _all_connections_saturated(b.a2(), b.a1()))    // must be changed
      return 1;

    return 0;     // all neighbours fully saturated, does not need changing
  }

  if (mol2graph.preserve_cc_double_bonds_no_heteroatoms())
  {
    if (! _all_connections_saturated(b.a1(), b.a2()))    // must be changed
      return 1;

    if (! _all_connections_saturated(b.a2(), b.a1()))    // must be changed
      return 1;

    if (attached_heteroatom_count(b.a1() > 0) || attached_heteroatom_count(b.a2() > 0))
      return 1;

    return 0;
  }

  return 1;
}

int
Molecule::change_to_graph_form (const Mol2Graph & mol2graph)
{
  int changes = 0;   // changes which require recomputation of nbonds & implicit hydrogens

  const int nb = _bond_list.number_elements();

  int * change_to_single = new_int(nb, 1); std::unique_ptr<int[]> free_change_to_single(change_to_single);

  if (mol2graph.some_kind_of_double_bond_preservation_active())
    compute_aromaticity_if_needed();

  for (int i = 0; i < nb; ++i)
  {
    Bond * b = _bond_list[i];

    if (b->is_single_bond())    // no need to change
      change_to_single[i] = 0;
    else if (b->is_aromatic())    // non single, aromatic bond, definitely change
      b->set_permanent_aromatic(0);
    else if (b->is_triple_bond())
      change_to_single[i] = ! mol2graph.exclude_triple_bonds_from_graph_reduction();
    else if (b->is_double_bond() && ! mol2graph.some_kind_of_double_bond_preservation_active())
      ;
    else if (! _double_bond_needs_changing_for_graph_form(*b, mol2graph))
      change_to_single[i] = 0;
  }

  for (int i = 0; i < nb; i++)
  {
    if (! change_to_single[i])
      continue;

    Bond * b = _bond_list[i];

    b->set_bond_type(SINGLE_BOND);
    b->set_permanent_aromatic(0);
    changes++;
  }

  if (mol2graph.remove_chiral_centres() && _chiral_centres.number_elements())
  {
    _chiral_centres.resize(0);
    changes++;
  }

  Set_of_Atoms atoms_to_be_removed;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (1 == a->atomic_number())
    {
      atoms_to_be_removed.add(i);
      changes++;
      continue;
    }

//  Remove charges unless it would create valence problems

    if (a->formal_charge())
    {
      if (7 == a->atomic_number() && 4 == a->ncon())
        ;
      else
      {
        a->set_formal_charge(0);
        changes++;
      }
    }

    if (a->implicit_hydrogens_known())
    {
      a->set_implicit_hydrogens_known(0);
      changes++;
    }

    if (a->is_isotope())
      a->set_isotope(0);
  }

  if (atoms_to_be_removed.number_elements())
    remove_atoms(atoms_to_be_removed);

  if (mol2graph.revert_all_directional_bonds_to_non_directional())
  {
    if (revert_all_directional_bonds_to_non_directional())
      changes++;
  }

  if (changes)
  {
    int unused;
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->recompute_nbonds();
      _things[i]->recompute_implicit_hydrogens(unused);
    }
  }

  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);

  return 1;
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
Molecule::write_set_of_bonds_as_mdl_v30_collection (const resizable_array<int> & b,
                                                    const const_IWSubstring & zname,
                                                    const const_IWSubstring & subname,
                                                    std::ostream & output) const
{
  output << "M  V30 BEGIN COLLECTION\n";
  output << "M  V30 " << zname;
  if (subname.length() > 0)
    output << '/' << subname;
  output << endl;

  IWString output_buffer;
  output_buffer.resize(200);

  int n = b.number_elements();

  assert(n <= _bond_list.number_elements());

  output_buffer << "M  V30 BONDS=" << n;

  for (int i = 0; i < n; i++)
  {
    output_buffer << ' ' << (b[i] + 1);
  }

  write_v30_record(output_buffer, output);

  output << "M  V30 END COLLECTION\n";

  return output.good();
}

int
Molecule::write_set_of_bonds_as_mdl_v30_collection (const int * b,
                                                    const const_IWSubstring & zname,
                                                    const const_IWSubstring & subname,
                                                    std::ostream & output) const
{
  output << "M  V30 BEGIN COLLECTION\n";
  output << "M  V30 " << zname;
  if (subname.length() > 0)
    output << '/' << subname;
  output << endl;

  IWString output_buffer;
  output_buffer.resize(200);

  int nb = _bond_list.number_elements();

  int n = count_non_zero_occurrences_in_array(b, nb);

  output_buffer << "M  V30 BONDS=" << n;

  for (int i = 0; i < nb; i++)
  {
    if (0 != b[i])
      output_buffer << ' ' << (i + 1);
  }

  write_v30_record(output_buffer, output);

  output << "M  V30 END COLLECTION\n";

  return output.good();
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
Bond_list::assign_bond_numbers (int istart)
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

  if (_things[0]->bond_number_assigned() &&
      _things[_number_elements - 1]->bond_number_assigned())
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
  assert (ok());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_permanent_aromatic(0);
  }

  return 1;
}

atom_number_t
Molecule::atom_with_isotope(int iso) const
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
Molecule::get_bond_types (bond_type_t * bt) const
{
  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    bt[i] = BOND_TYPE_ONLY_MASK & _bond_list[i]->btype();
  }
  
  return 1;
}

int
Molecule::set_bond_types_no_set_modified (const bond_type_t * bt)
{
  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    _bond_list[i]->set_bond_type(bt[i]);
  }

  return 1;
}

Mol2Graph::Mol2Graph ()
{
  _exclude_triple_bonds_from_graph_reduction = 0;
  _revert_all_directional_bonds_to_non_directional = 1;
  _preserve_cc_double_bonds_saturated = 0;
  _preserve_cc_double_bonds_no_heteroatoms = 0;
    
  _remove_chiral_centres = 1;

  _append_molecular_formula = 0;

  _aromatic_distinguishing_formula = 0;

  return;
}

int
Mol2Graph::debug_print (std::ostream & output) const
{
  output << "Mol2Graph:debug_print:\n";
  output << " _remove chiral centres " << _remove_chiral_centres << '\n';
  output << " _exclude_triple_bonds_from_graph_reduction " << _exclude_triple_bonds_from_graph_reduction << '\n';
  output << " _revert_all_directional_bonds_to_non_directional " << _revert_all_directional_bonds_to_non_directional << '\n';
  output << " _preserve_cc_double_bonds_saturated " << _preserve_cc_double_bonds_saturated << '\n';
  output << " _preserve_cc_double_bonds_no_heteroatoms " << _preserve_cc_double_bonds_no_heteroatoms << '\n';
  output << " _append_molecular_formula " << _append_molecular_formula << '\n';
  output << " _aromatic_distinguishing_formula " << _aromatic_distinguishing_formula << '\n';

  return output.good();
}

int
Mol2Graph::construct (Command_Line & cl, 
                      const char flag,
                      const int verbose)
{
  const_IWSubstring h;
  for (int i = 0; cl.value(flag, h, i); ++i)
  {
    if ('s' == h)
    {
      _preserve_cc_double_bonds_saturated = 1;
      if (verbose)
        cerr << "During graph reduction, will preserve double bonds adjacent to fully saturated atoms\n";
    }
    else if ('c' == h)
    {
      _preserve_cc_double_bonds_no_heteroatoms = 1;
      if (verbose)
        cerr << "During graph reduction, will preserve double bonds adjacent to only carbon atoms\n";
    }
    else if ("chiral" == h)
    {
      _remove_chiral_centres = 0;
      if (verbose)
        cerr << "Will NOT remove chiral centres during graph reduction\n";
    }
    else if ("rmchiral" == h)
    {
      _remove_chiral_centres = 1;
      if (verbose)
        cerr << "Will remove chiral centres during graph reduction\n";
    }
    else if ("help" == h)
    {
      cerr << "Mol2Graph::construct:the following options are recognised\n";
      cerr << " -" << flag << " s          preserve C=C bonds that are adjacent to fully saturated atoms\n";
      cerr << " -" << flag << " c          preserve C=C bonds that are adjacent to only carbon atoms\n";
      cerr << " -" << flag << " chiral     preserve chiral centres during graph reduction\n";
      cerr << " -" << flag << " rmchiral   remove chiral centres during graph reduction\n";
      exit(0);
    }
    else
    {
      cerr << "Mol2Graph::construct:unrecognised -" << flag << " qualifier '" << h << "'\n";
      return 0;
    }
  }

  return 1;
}

static uint64_t
do_or(uint64_t x, int bits_left, 
      const uint64_t maxval,
      const uint64_t rc)
{
  if (x > maxval)
    x = maxval;

  return rc | (x << bits_left);
}

uint64_t
Molecule::quick_bond_hash()
{
  if (0 == _number_elements)
    return 0;

  compute_aromaticity_if_needed();

  uint64_t c_single_c = 0;
  uint64_t c_double_c = 0;
  uint64_t c_triple_c = 0;   // 2 bits
  uint64_t c_arom_c = 0;

  uint64_t c_single_n = 0;
  uint64_t c_double_n = 0;
  uint64_t c_triple_n = 0;    // 3 bits
  uint64_t c_arom_n = 0;

  uint64_t c_single_o = 0;
  uint64_t c_double_o = 0;
  uint64_t c_arom_o   = 0;       // 3 bits

  uint64_t c_single_s = 0;
  uint64_t c_double_s = 0;     // 2 bits
  uint64_t c_arom_s = 0;       // 3 bits

  uint64_t n_arom_n = 0;         // 2 bits
  uint64_t n_single_n = 0;         // bug, we are overlooking this
  uint64_t n_double_n = 0;         // bug, we are overlooking this
  uint64_t n_triple_n = 0;         // bug, we are overlooking this

  uint64_t n_single_o = 0;    // bug, overlooked
  uint64_t n_double_o = 0;
  uint64_t n_arom_o = 0;      // bug, overlooked

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

  int bshift = 64-4;

  rc = do_or(c_single_c, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_double_c, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_triple_c, bshift, 2, rc); bshift -= 2;
  rc = do_or(c_arom_c/4, bshift, 8, rc); bshift -= 4;

  rc = do_or(c_single_n, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_double_n, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_triple_n, bshift, 4, rc); bshift -= 3;
  rc = do_or(c_arom_n,   bshift, 8, rc); bshift -= 4;

  rc = do_or(c_single_o, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_double_o, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_arom_o,   bshift, 1, rc); bshift -= 1;

  rc = do_or(c_single_s, bshift, 8, rc); bshift -= 4;
  rc = do_or(c_double_s, bshift, 2, rc); bshift -= 2;
  rc = do_or(c_arom_s,   bshift, 4, rc); bshift -= 3;

//rc = do_or(n_double_o, bshift, 8, rc); bshift -= 4

  rc = do_or(o_single_s, bshift, 8, rc); bshift -= 4;
  rc = do_or(o_double_s, bshift, 8, rc); bshift -= 4;

  rc = do_or(aliphatic_ring_bond_count/2, bshift, 8, rc); bshift -= 4;

  assert (bshift >= 0);
//if (bshift < 0)
//  cerr << "negative bshift\n";

  return rc;
}

const_BondIterator::const_BondIterator(const Atom * a, const atom_number_t zatom) : _atom(&a), _atnum(zatom)
{
  _b = a->cbegin();

  return;
}
