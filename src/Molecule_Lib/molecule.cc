#include <iostream>
#include <iomanip>
#include <utility>
#include <memory>
using std::cerr;
using std::endl;

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#define COMPILING_CTB
#define COMPILING_MOLECULE_MAIN

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

/*
  A molecule.
  The object model is complex, because the decision was to base
  most manipulations on atom numbers. Atom objects do not know
  very much at all, the molecule must be consulted for almost
  everything.
*/

/*
  The charge array is initialised to NULL. As a molecule
  is built, if we encounter a non-zero charge, we allocate
  an array for it, and fill all members with 0.0 (except
  for the non zero charge just encountered).
*/

#include "iwminmax.h"
#include "misc.h"

#include "molecule.h"
#include "misc2.h"
#include "chiral_centre.h"
#include "path.h"
#include "pearlman.h"
#include "smiles.h"

static int display_already_bonded_error_message = 1;

void 
set_display_already_bonded_error_message(int s)
{
  display_already_bonded_error_message = s;
}

static char useThisCharAsDotInSmiles = '.';

void 
setUseThisCharAsDotInSmiles(char thisChar)
{
  useThisCharAsDotInSmiles = thisChar;
}

char
getUseThisCharAsDotInSmiles(void)
{
  return useThisCharAsDotInSmiles ;
}

static charge_t _max_reasonble_atomic_partial_charge_value = static_cast<charge_t>(3.0);
static charge_t _min_reasonble_atomic_partial_charge_value = static_cast<charge_t>(-3.0);

int
set_max_reasonble_atomic_partial_charge_value(charge_t s)
{
  if (s > _min_reasonble_atomic_partial_charge_value)
  {
    _max_reasonble_atomic_partial_charge_value = s;
    return 1;
  }

  cerr << "Invalid max reasonable charge " << s << " min is " << _min_reasonble_atomic_partial_charge_value << endl;
  return 0;
}

int
set_min_reasonble_atomic_partial_charge_value (charge_t s)
{
  if (s < _max_reasonble_atomic_partial_charge_value)
  {
    _min_reasonble_atomic_partial_charge_value = s;
    return 1;
  }

  cerr << "Invalid min reasonable charge " << s << " max is " << _max_reasonble_atomic_partial_charge_value << endl;
  return 0;
}

int
set_reasonable_atomic_partial_charge_range(charge_t mn, charge_t mx)
{
  if (mn < mx)
  {
    _min_reasonble_atomic_partial_charge_value = mn;
    _max_reasonble_atomic_partial_charge_value = mx;

    return 1;
  }

  cerr << "INvalid range of valid partial atomic charges " << mn << " and " << mx << "\n";
  return 0;
}

int
reasonable_atomic_partial_charge_value (charge_t q)
{
  if (q < _min_reasonble_atomic_partial_charge_value)
    return 0;

  if (q > _max_reasonble_atomic_partial_charge_value)
    return 0;

  return 1;
}

static int invalidate_bond_list_ring_info_during_invalidate_ring_info = 1;

void
set_invalidate_bond_list_ring_info_during_invalidate_ring_info (int s)
{
  invalidate_bond_list_ring_info_during_invalidate_ring_info = s;
}

void
Molecule::_default_values (int atoms_in_new_molecule)
{
  assert(atoms_in_new_molecule >= 0);

  _magic   = MOLECULE_MAGIC_NUMBER;

  _charges = NULL;
  _atom_type = NULL;

  _distance_matrix = NULL;

  _partially_built = 0;

  _nrings = IW_NRINGS_NOT_COMPUTED;
  _number_sssr_rings = IW_NRINGS_NOT_COMPUTED;
  _ring_membership = NULL;

  _aromaticity = NULL;

  _fragment_information.invalidate();

  if (atoms_in_new_molecule > 0)
    resize(atoms_in_new_molecule);

  _user_specified_void_ptr = NULL;
  
  doNotComputeAromaticity = false; 

  return;
}

static int copy_name_in_molecule_copy_constructor = 0;

void
set_copy_name_in_molecule_copy_constructor (int s)
{
  copy_name_in_molecule_copy_constructor = s;
}

Molecule::Molecule (int atoms_in_new_molecule)
{
  assert(atoms_in_new_molecule >= 0);

  _default_values(atoms_in_new_molecule);

  return;
}

Molecule::Molecule (const Molecule & rhs)
{
  _default_values(0);

  add_molecule(&rhs);

  if (copy_name_in_molecule_copy_constructor)
    _molecule_name = rhs._molecule_name;

  return;
}

/*
  Both the destructor and delete_all_atoms_and_bonds need to free all dynamically allocated arrays.
*/

int
Molecule::_free_all_dynamically_allocated_things()
{
  DELETE_IF_NOT_NULL(_charges);

  DELETE_IF_NOT_NULL(_atom_type);

  if (NULL != _distance_matrix)
  {
    delete [] _distance_matrix;
    _distance_matrix = NULL;
  }

  if (NULL != _aromaticity)
  {
    delete [] _aromaticity;
    _aromaticity = NULL;
  }

  if (NULL != _ring_membership)
  {
    delete [] _ring_membership;
    _ring_membership = NULL;
  }

  return 1;
}

Molecule::~Molecule()
{
  assert(ok());

#ifndef NDEBUG
  if (-7373 == _magic)
    cerr << "Deleting an already deleted molecule '" << _molecule_name << "'\n";
#endif

  _free_all_dynamically_allocated_things();

  _smiles_information.invalidate();

#ifndef NDEBUG
  _magic = -7373;
#endif

  return;
}

int
Molecule::delete_all_atoms_and_bonds()
{
  assert(ok());

  invalidate_smiles();

  _invalidate_ring_info();

  _set_modified();

  _chiral_centres.resize(0);

  resize(0);

  _bond_list.resize(0);

  _symmetry_class_and_canonical_rank.invalidate();

  _free_all_dynamically_allocated_things();

  return 1;
}

Molecule &
Molecule::operator= (const Molecule & rhs)
{
  delete_all_atoms_and_bonds();

  _molecule_name = rhs._molecule_name;

  _text_info.resize(0);

  add_molecule(&rhs);

  return *this;
}


Molecule &
Molecule::operator= (Molecule && rhs)
{
//cerr << "Molecule::operator move\n";

  delete_all_atoms_and_bonds();

  rhs.invalidate_smiles();

  _molecule_name = std::move(rhs._molecule_name);

  resizable_array_p<Atom> * myatoms = this;
  resizable_array_p<Atom> * rhsatoms = &rhs;
  *myatoms = std::move(*rhsatoms);

  _bond_list = std::move(rhs._bond_list);

  _chiral_centres = std::move(rhs._chiral_centres);

//cerr << "We have " << _number_elements << " atoms, RHS has " << rhs._number_elements << " atoms\n";

//rhs.delete_all_atoms_and_bonds();

  return *this;
}

bool
Molecule::operator== ( Molecule &rhs )
{
  if (_number_elements != rhs._number_elements)
    return false;

  if (_bond_list.number_elements() != rhs._bond_list.number_elements())
    return false;

  if (number_fragments() != rhs.number_fragments())
    return false;

  if (nrings() != rhs.nrings())
    return false;

  if (_chiral_centres.number_elements() != rhs._chiral_centres.number_elements())
    return false;

  if( unique_smiles() == rhs.unique_smiles() )
    return true;
  else
    return false;
}

int
Molecule::debug_print (std::ostream & os) const
{
  assert(os.good());
  os << "Molecule::debug_print " << this << ", information, " << _number_elements << " atoms";
  os << ' ' << _bond_list.number_elements() << " bonds";
  os << endl;

  if (! ok())
    os << "Warning, OK failed\n";
  if (! _ok_ring_info())
    os << "Warning, OK RING INFO failed\n";

  if (_molecule_name != "")
    os << "Molecule name '" << _molecule_name << "'\n";
  else
    os << "No name\n";

  _smiles_information.debug_print(os);

  if (NULL != _symmetry_class_and_canonical_rank.symmetry_class())
    os << "Symmetry class array allocated\n";

  if (NULL != _symmetry_class_and_canonical_rank.canonical_rank())
    os << "Canonical order function allocated\n";

  if (NULL != _aromaticity)
    os << "Aromaticity data is available\n";

  if (_charges)
    os << _charges->number_elements() << " charges, type '" << _charges->ztype() << "'\n";

  if (_atom_type)
    os << _atom_type->number_elements() << " atom types, type '" << _atom_type->ztype() << "'\n";

  if (_fragment_information.contains_valid_data())
    _fragment_information.debug_print(os);

  int hcd = highest_coordinate_dimensionality();

  os << "Highest coordinate dimensionality " << hcd << endl;

  charge_t net_charge = static_cast<charge_t>(0.0);

  const int * canonical_rank = _symmetry_class_and_canonical_rank.canonical_rank();
  const int * symmetry_class = _symmetry_class_and_canonical_rank.symmetry_class();

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];    // not const because a->implicit_hydrogens is non const

    os << std::setw(4) << i << " " << std::setw(2) << a->atomic_symbol();
    if (a->atom_map())
      os << ' ' << a->atom_map();

    os << ' ';

    if (! a->implicit_hydrogens_computed())
      os << "?  implicit H";
    else
    {
      os << a->implicit_hydrogens();
      if (a->implicit_hydrogens_known())
        os << '*';
      else
        os << ' ';
      os << " implicit H";
    }

    if (a->number_elements())
      os << " (ncon " << a->number_elements() << ')';

    if (NULL != canonical_rank)
      os << " canon = " << canonical_rank[i];

    if (NULL != symmetry_class && IW_SYMMETRY_CLASS_UNDEFINED != symmetry_class[i])
      os << " sym = " << symmetry_class[i];

    if (_fragment_information.contains_valid_data())
      os << " (frag " << _fragment_information.fragment_membership()[i] << ")";

    if (a->formal_charge())
      os << " (fc " << a->formal_charge() << ')';
    if (_charges)
    {
      os << " (q " << std::setw(7) << _charges->item(i) << ')';
      net_charge += _charges->item(i);
    }

    if (hcd > 1)
      os << " (" << a-> x() << "," << a->y() << "," << a->z() << ") ";

    int icon = ncon(i);
    if (icon)
    {
      os << " bonded to";
      for (int j = 0; j < ncon(i); j++)
        os << ' ' << other(i, j);
    }
    else
      os << " unconnected";

    if (NULL == _aromaticity)
      ;
    else if (IS_AROMATIC_ATOM(_aromaticity[i]))
      os << " aromatic";
    else
      os << " aliph";

    if (a->isotope())
      os << " ISO " << a->isotope();

    os << endl;


    if (! a->audit())
      os << "Warning, audit function fails for this atom\n";
  }

  if (_charges)
    os << "Total net charge " << net_charge << endl;

  int nca = _chiral_centres.number_elements();
  for (int i = 0; i < nca; i++)
  {
    const Chiral_Centre * c = _chiral_centres[i];
    os << "Chiral Center " << i << ' ';
    print_chiral_centre_details(c, os);
  }

  _bond_list.debug_print(os);

  print_ring_info(os);

  return 1;
}

int
Molecule::ok_atom_number(atom_number_t a) const
{
  if (! ok())
    return 0;

  if (a < 0 || a >= _number_elements)
    return 0;

  return 1;
}

int
Molecule::ok_2_atoms (atom_number_t a1, atom_number_t a2) const
{
  if (! ok())
    return 0;

  if (a1 < 0 || a1 >= _number_elements)
    return 0;
  if (a2 < 0 || a2 >= _number_elements)
    return 0;
  if (a1 == a2)
    return 0;

  return 1;
}

int
Molecule::ok_3_atoms (atom_number_t a1, atom_number_t a2, atom_number_t a3) const
{
  if (! ok_2_atoms(a1, a2))
    return 0;

  if (a3 < 0 || a3 >= _number_elements)
    return 0;

  if (a3 == a1 || a3 == a2)
    return 0;

  return 1;
}

int
Molecule::ok_4_atoms (atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4) const
{
  if (! ok_3_atoms(a1, a2, a3))
    return 0;

  if (a4 < 0 || a4 >= _number_elements)
    return 0;

  if (a1 == a4 || a2 == a4 || a3 == a4)
    return 0;

  return 1;
}

int
Molecule::check_bonding() const
{
  assert(ok());

  if (0 == _number_elements)
  {
    cerr << "check_bonding: warning, empty molecule encountered\n";
    assert(0 == _bond_list.number_elements());
  }

//  Check to make sure that all bonded atoms are within range,
//  and also that bonding info is symmetric

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];
    if (! a->audit())
    {
      cerr << "check bonding: bad atom found, at address " << &a << "\n";
      a->debug_print(cerr);
      iwabort();
    }
    
    const int icon = a->ncon();

    Set_of_Atoms connected(icon);

    for (int j = 0; j < icon; j++)
    {
      const Bond * b = a->item(j);
      atom_number_t k = b->other(i);

      if (INVALID_ATOM_NUMBER == k || (! ok_index(k)))
      {
        cerr << "check bonding: bad connection " << i << ' ' << j << ' ' << k << "\n";
	return 0;
      }

      if (! _things[k]->is_bonded_to(i))    // check symmetry of bonding
      {
        cerr << "check bonding: asymetric bond, atoms " << i << " and " << k << "\n";
        return 0;
      }

      if (i == b->a1())
      {
        if (0 == connected.number_elements())
          connected.add(b->a2());
        else if (0 == connected.add_if_not_already_present(b->a2()))
        {
          cerr << "Molecule::check_bonding:atom " << i << " has multiple bonds to " << b->a2() << endl;
          return 0;
        }
      }
      else if (i == b->a2())
      {
        if (0 == connected.number_elements())
          connected.add(b->a1());
        else if (0 == connected.add_if_not_already_present(b->a1()))
        {
          cerr << "Molecule::check_bonding:atom " << i << " has multiple bonds to " << b->a1() << endl;
          return 0;
        }
      }
      else
      {
        cerr << "Molecule::check_bonding:bond not involving owner, atom " << i << ' ' << *b << endl;
        return 0;
      }
    }
  }

  if (! check_ring_info())
    return 0;

  return 1;
}

/*
  Audits a molecule for chemical reasonableness.
  Pretty crude now...
  We return the number of definite problems
*/

int
Molecule::check_chemistry() const
{
  assert(ok());

  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
//  atomic_number_t z = atomic_number(i);
//  int icon = ncon(i);
//  int ibonds = nbonds(i);
//  formal_charge_t fcharge = formal_charge(i);
  }

  return rc;
}

/*
  There is a special function for adding a new atom, because the charge
  array (if present) must be kept in sync with the number of atoms.
*/

int
Molecule::add (Atom * a, int partial_molecule)
{
  assert(ok());
  assert(OK_ATOM(a));

  resizable_array_p<Atom>::add(a);

  if (_charges)
    _charges->add(static_cast<charge_t>(0.0));

  if (_atom_type)
    _atom_type->add(static_cast<atom_type_t>(0));

  if (! partial_molecule)
    _set_modified();

  return 1;
}

int
Molecule::add (const Element * e)
{
  Atom * a = new Atom(e);

  return add(a);
}

int
Molecule::attached_heteroatom_count (atom_number_t zatom) const
{
//assert(ok_atom_number(zatom));   // used within pearlman.cc, no check..

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int rc = 0;
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);
    atomic_number_t zj = _things[j]->atomic_number();
    if (1 != zj && 6 != zj)
      rc++;
  }

  return rc;
}

/*
  In Nov 96 change this from counting the number of multiple bonds
  to heteroatoms, to just reporting presence or absence.
  If EXCLUDE is specified, the bond from A to EXCLUDE will be ignored
*/

int
Molecule::multiple_bond_to_heteroatom (atom_number_t zatom,
                                       atom_number_t exclude) const
{
  assert(INVALID_ATOM_NUMBER == exclude ? (OK_ATOM_NUMBER(this, zatom)) :
                          (OK_2_ATOMS(this, zatom, exclude)));

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (exclude == j)
      continue;

    atomic_number_t zj = _things[j]->atomic_number();
    if (6 == zj || 1 == zj)
      continue;

    if (b->is_aromatic())
      continue;

    if (! b->is_single_bond())
      return 1;
  }

  return 0;
}

/*
  Similar to the above, but now there is a list of excluded atoms
  in an array
*/

int
Molecule::multiple_bond_to_heteroatom (atom_number_t zatom,
                                       const int * exclude) const
{
  assert(ok_atom_number(zatom));

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int rc = 0;
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (exclude[j])
      continue;

    const atomic_number_t zj = _things[j]->atomic_number();
    if (6 == zj || 1 == zj)
      continue;

    if (b->is_aromatic())
      continue;

    if (! b->is_single_bond())
      rc++;
  }

  return rc;
}

int
Molecule::doubly_bonded_oxygen_count (atom_number_t zatom) const
{
  Atom * a = _things[zatom];    // not const because it might compute _nbonds

  const int acon = a->ncon();

  const int max_available = a->nbonds() - acon;

  if (0 == max_available)
    return 0;

  int rc = 0;

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    const atom_number_t o = b->other(zatom);

    if (8 == _things[o]->atomic_number())
    {
      rc++;
      if (rc == max_available)
        return rc;
    }
  }

  return rc;
}

void
Molecule::allocate_charges()
{
  assert(ok());
  assert(NULL == _charges);

  _charges = new Set_of_Charges;

  assert(NULL != _charges);

  if (_number_elements > 0)
    _charges->extend(_number_elements, static_cast<charge_t>(0.0));

  return;
}

int
Molecule::set_charge (atom_number_t i, charge_t qq)
{
  assert(ok_atom_number(i));

  if (! reasonable_atomic_partial_charge_value(qq))
  {
    cerr << "Molecule::set_charge: invalid charge " << qq << " for atom " << i << endl;
    return 0;
  }

  if (NULL == _charges)
    allocate_charges();

  _charges->seti(i, qq);

  return 1;
}

void
Molecule::set_charges (const charge_t q [], const const_IWSubstring & qt)
{
  assert(ok());

  if (NULL == _charges)
    allocate_charges();

  for (int i = 0; i < _number_elements; i++)
  {
    _charges->seti(i, q[i]);
  }

  _charges->set_type(qt);

  return;
}

void
Molecule::_set_partial_charge_type (const const_IWSubstring & qtype)
{
  if (NULL == _charges)
    allocate_charges();

  _charges->set_type(qtype);

  return;
}

static IWString empty_string;

const IWString &
Molecule::partial_charge_type() const 
{
  if (NULL == _charges)
    return empty_string;

  return _charges->ztype();
}

/*
  Special care is needed with the sticky bit for hcounts. We overwrite it.
*/

void
Molecule::set_formal_charge (atom_number_t zatom, formal_charge_t qq)
{
  assert(ok_atom_number(zatom));
  assert(reasonable_formal_charge_value(qq));

  _things[zatom]->set_formal_charge(qq);

  _set_modified(zatom);

// If we are removing a positive charge from a nitrogen, remove any chiral center on that atom.
// Should this be an optional behaviour??

  if (0 != qq)   // no, we set some kind of formal charge, we are not resetting a Nitrogen
    return;

  const Atom * a = _things[zatom];

  if (7 != a->atomic_number())    // not on a nitrogen
    return;

  const int acon = a->ncon();

  if (acon < 3)   // cannot be a chiral centre there
    return;

// We may have an N atom with 3 connections and an implicit Hydrogen, or 4 connections including an explicit H
// But we really do not know what is going on, so let's just remove any chiral centre that might be there..

  const int nchiral = _chiral_centres.number_elements();

  for (int i = 0; i < nchiral; ++i)
  {
    if (zatom == _chiral_centres[i]->a())
    {
      _chiral_centres.remove_item(i);
      return;
    }
  }

  return;
}

int
Molecule::set_formal_charge_if_different (atom_number_t zatom, formal_charge_t qq)
{
  assert(ok_atom_number(zatom));
  assert(reasonable_formal_charge_value(qq));

  if (qq == _things[zatom]->formal_charge())
    return 0;

  _things[zatom]->set_formal_charge(qq);

  _set_modified(zatom);

// should do the check for chirality on nitrogen atoms...

  return 1;
}

/*
*/

int
Molecule::resize (int new_size)
{
  assert(ok());
  assert(new_size >= 0);

  if (new_size > _elements_allocated)
  {
    resizable_array_p<Atom>::resize(new_size);

    if (NULL != _charges)
      _charges->resize(new_size);

    if (NULL != _atom_type)
      _atom_type->resize(new_size);

    return 1;
  }

  _set_modified();

  resizable_array_p<Atom>::resize(new_size);

  if (NULL != _charges)
    _charges->resize(new_size);

  if (NULL != _atom_type)
    _atom_type->resize(new_size);

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->molecule_being_resized(new_size);
  }

  int nb = _bond_list.number_elements();
  for (int i = nb - 1; i >= 0; i--)
  {
    const Bond * b = _bond_list[i];
    if (b->a1() >= new_size || b->a2() >= new_size)
      _bond_list.remove_item(i);
  }

  int ncc = _chiral_centres.number_elements();
  for (int i = ncc - 1; i >= 0; i--)
  {
    const Chiral_Centre * c = _chiral_centres[i];
    if (c->left_down() >= new_size || c->right_down() >= new_size ||
        c->top_front() >= new_size || c->top_back()   >= new_size ||
        c->a() >= new_size)
      _chiral_centres.remove_item(i);
  }

  _remove_directionality_from_bonds_not_actually_directional();

  return new_size;
}

/*
  The private _resize just does a dumb resize! Use with care!
*/

/*void
Molecule::_resize (int new_size)
{
  assert(new_size >= 0);

  resizable_array_p<Atom>::resize(new_size);

  if (NULL != _charges)
    _charges->resize(new_size);

  if (NULL != _atom_type)
    _atom_type->resize(new_size);

  if (0 == new_size)
    _bond_list.resize(0);

  if (new_size < _number_elements)
    _symmetry_class_and_canonical_rank.invalidate();

  return;
}*/

//#define SHOW_OK 1

int
Molecule::ok() const
{
#ifdef SHOW_OK
  if (NULL == this)
    cerr << "Null molecule\n";
  cerr << "Checking Molecule " <<   " array " << resizable_array_p<Atom>::ok() << 
          " magic " << (MOLECULE_MAGIC_NUMBER == _magic) << endl;
#endif

  if (0 == resizable_array_p<Atom>::ok())
  {
#ifdef SHOW_OK
    cerr << "resizable_array_p<Atom>::ok() failed\n";
#endif
    return 0;
  }

  if (MOLECULE_MAGIC_NUMBER != _magic)
  {
#ifdef SHOW_OK
    cerr << "MAGIC NUMBER IS WRONG\n";
#endif
    return 0;
  }
  if (NULL == _charges)
    ;
  else if (_charges->number_elements() > 0 &&
           _number_elements != _charges->number_elements())
  {
#ifdef SHOW_OK
    cerr << "Charge mis-match, " << _number_elements << " atoms, and " <<
            _charges->number_elements() << " charges\n";
#endif

    return 0;
  }

  if (NULL == _atom_type)
    ;
  else if (_atom_type->number_elements() > 0 &&
           _number_elements != _atom_type->number_elements())
  {
#ifdef SHOW_OK
    cerr << "Atom type mis-match, " << _number_elements << " atoms, and " <<
            _atom_type->number_elements() << " atom types\n";
#endif

    return 0;
  }

#ifdef SHOW_OK
  cerr << "Checking ring info " << _nrings << " sssr = " << _sssr_rings.number_elements() << endl;
#endif

  if (IW_NRINGS_NOT_COMPUTED == _nrings)
    ;
  else if (_nrings < 0)
    return 0;
  else if (_nrings > _bond_list.number_elements())
    return 0;

#ifdef SHOW_OK
  cerr << "Checking fragment info " << _fragment_information.number_fragments() << endl;
#endif

//if (_fragment_information. < 0)    // not computed
//  ;
//else if (_bonds_in_fragment.number_elements() != _number_fragments)
//  return 0;
//else if (_atoms_in_fragment.number_elements() != _number_fragments)
//  return 0;

  if (! _ok_ring_info())
    return 0;

#ifdef SHOW_OK
  if (! _bond_list.ok())
    cerr << "Bond list is bad\n";
#endif

  if (! _bond_list.ok())
    return 0;

#ifdef SHOW_OK
  cerr << "Checking chiral centres\n";
#endif

//if (! _check_chiral_centres())
//  return 0;

#ifdef SHOW_OK
  cerr << "Molecule is OK\n";
#endif

  return 1;
}

int
Molecule::has_charges() const
{
  assert(ok());

  return (NULL != _charges);
}


void
Molecule::allocate_atom_types()
{
  assert(ok());
  assert(NULL == _atom_type);

  _atom_type = new Atom_Types;
  _atom_type->extend(_number_elements, static_cast<atom_type_t>(INVALID_ATOM_TYPE));

  return;
}

int
Molecule::has_atom_types() const
{
  assert(ok());

  return (NULL != _atom_type);
}

int
Molecule::copy_atom_types (const Molecule & m2)
{
  assert(ok());
  assert(m2.ok());

  if (_number_elements != m2._number_elements)
  {
    cerr << "molecule::copy_atom_types: atom count mismatch " << _number_elements << " vs " <<
            m2._number_elements << "\n";
    return 0;
  }

  if (! m2.has_atom_types())
  {
    if (NULL != _atom_type)
      invalidate_atom_types();

    return 1;
  }

  if (! has_atom_types())
    allocate_atom_types();

  *_atom_type = *(m2._atom_type);

  return 1;
}

void
Molecule::invalidate_atom_types()
{
  if (NULL != _atom_type)
  {
    delete _atom_type;
    _atom_type = NULL;
  }

  return;
}

atom_type_t
Molecule::atom_type (atom_number_t i) const
{
  assert(ok_atom_number(i));

  if (! has_atom_types())
    return static_cast<atom_type_t>(INVALID_ATOM_TYPE);

  return _atom_type->item(i);
}

Atom_Types &
Molecule::atom_types()
{
  assert(ok());

  if (NULL == _atom_type)
    allocate_atom_types();

  return * _atom_type;
}

void
Molecule::set_atom_type (atom_number_t a, atom_type_t t)
{
  assert(ok_atom_number(a));

  if (! has_atom_types())
    allocate_atom_types();

  _atom_type->seti(a, t);

  return;
}


int
Molecule::has_formal_charges() const
{
  assert(ok());
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 != _things[i]->formal_charge())
      return 1;
  }

  return 0;
}

int
Molecule::has_no_formal_charges() const
{
  return ! has_formal_charges();
}

int
Molecule::number_formally_charged_atoms() const
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->formal_charge())
      rc++;
  }

  return rc;
}

int
Molecule::net_formal_charge() const
{
  int rc = 0;
  for (int i = 0; i < _number_elements; ++i)
  {
    rc += _things[i]->formal_charge();
  }

  return rc;
}

/*
  many times we need the address of the i'th atom in a molecule
*/

const Atom *
Molecule::atomi (int i) const
{
  assert(ok_atom_number(i));

  return _things[i];
}

int
Molecule::atoms (const Atom ** a) const
{
  for (int i = 0; i < _number_elements; i++)
    a[i] = _things[i];

  return _number_elements;
}

const Atom &
Molecule::atom (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return *(_things[i]);
}

const Element *
Molecule::elementi (int i) const
{
  assert(ok_atom_number(i));

  return _things[i]->element();
}

int
Molecule::set_element (atom_number_t a, const Element * e)
{
  assert(ok_atom_number(a));
  assert(OK_ELEMENT(e));

  Atom * aa = _things[a];
  aa->set_element(e);

  _set_modified(a);

  return 1;
}

int
Molecule::set_atomic_number (atom_number_t a, atomic_number_t z)
{
  assert(ok_atom_number(a));

  const Element * e = get_element_from_atomic_number(z);
  assert(NULL != e);

  _things[a]->set_element(e);

  _set_modified(a);

  return 1;
}

const Element &
Molecule::element (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->elementq();
}

atomic_number_t
Molecule::atomic_number (int i) const
{
  assert(ok_atom_number(i));

  return _things[i]->atomic_number();
}

void
Molecule::atomic_numbers (atomic_number_t * z) const
{
  assert(ok());
  assert(NULL != z);

  for (int i = 0; i < _number_elements; i++)
  {
    z[i] = _things[i]->atomic_number();
  }

  return;
}

/*
  Returns the number of connections for the I'th atom in a molecule
  Note this public interface does an OK_ATOM_NUMBER check
*/

int
Molecule::ncon (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->ncon();
}

/*
  Fill an array with the ncon values for each atom. Often used for
  efficiency...
  Return the highest connectivity found.
*/

int
Molecule::ncon (int * con) const
{
  assert(ok());

  iwmax<int> maxcon(0);
  for (int i = 0; i < _number_elements; i++)
  {
    con[i] = _things[i]->ncon();
//  cerr << " atom " << i << " value " << con[i] << " type " << _things[i]->element()->symbol() << endl;
    maxcon.extra(con[i]);
  }

  return maxcon.maxval();
}

int
Molecule::ncon (resizable_array<int> & con) const
{
  assert(ok());
  con.extend(_number_elements);

  iwmax<int> maxcon(0);
  for (int i = 0; i < _number_elements; i++)
  {
    int ic = _things[i]->ncon();
    con[i] = ic;
    maxcon.extra(ic);
  }

  return maxcon.maxval();
}

int
Molecule::maximum_connectivity() const
{
  iwmax<int> rc(0);
  for (int i = 0; i < _number_elements; i++)
  {
    rc.extra(_things[i]->number_elements());
  }

  return rc.maxval();
}

int
Molecule::nbonds (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->nbonds();
}

int
Molecule::nbonds (int * bonds) const
{
  assert(NULL != bonds);

  for (int i = 0; i < _number_elements; i++)
    bonds[i] = _things[i]->nbonds();

  return _number_elements;
}

/*
  Something about atom A has been modified. Invalidate all computed
  properties which are tied to atom A
*/

int
Molecule::_set_modified (atom_number_t a)
{
  invalidate_smiles();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);

  _symmetry_class_and_canonical_rank.invalidate();

  _things[a]->set_modified();

//invalidate_fragment_membership();   not needed, only some property of the atom has changed

// We must notify all rings that aromaticity is now unknown

  return 1;
}

/*
  The molecule has been modified. We need to invalidate many whole
  molecule computed properties.
  Note that this does NOT call set_modified for each atom!

  Jun 97: At one stage I had a return if _ring_membership was NULL,
  but then iwfrag died. Even if _ring_membership is NULL, we still
  have to check the _aromaticity array, because the _aromaticity
  array can exist even if the molecule has no rings
*/

int
Molecule::_set_modified()
{
  assert(ok());

  return _set_modified_no_ok();
}

/*
  Due to unknown reasons, the internal state may be inconsistent. Recover
  if possible
*/

int
Molecule::invalidate_from_possibly_invalid_state()
{
  return _set_modified_no_ok();
}

class Bond_Invalidator
{
  private:
  public:
    int operator() (Bond * b) const;
};

/*
  During set_modified, we need to do a bunch of things to the bond list.
*/

int
Bond_Invalidator::operator() (Bond * b) const
{
  b->set_non_aromatic();
  b->invalidate_bond_number();
  b->invalidate_nrings();

  return 1;
}

int
Molecule::_set_modified_no_ok()
{
  invalidate_smiles();

  DELETE_IF_NOT_NULL_ARRAY(_aromaticity);    // must be after invalidate_smiles()

  _symmetry_class_and_canonical_rank.invalidate();

  if (NULL != _ring_membership)
    _invalidate_ring_info();

  DELETE_IF_NOT_NULL_ARRAY(_distance_matrix);

  _bond_list.invalidate_bond_numbers();

  if (invalidate_bond_list_ring_info_during_invalidate_ring_info)
    _bond_list.invalidate_ring_info();

  _nrings = IW_NRINGS_NOT_COMPUTED;
  _number_sssr_rings = IW_NRINGS_NOT_COMPUTED;

  invalidate_fragment_membership();

  return 1;
}

int
Molecule::invalidate_smiles()
{
// No call to ok() on purpose

  _smiles_information.invalidate();

// If any bonds had been assigned aromaticity, we must reset them.
// Sept 98, not sure why this is here - it creates problems with testing iwfp.
// For now, I'll leave it...

/* Dec 2009. the non aromatic setting is now done in _invalidate_ring_info
  if (NULL != _aromaticity && locate_item_in_array (AROMATIC, _number_elements, _aromaticity) >= 0)
  {
    int nb = _bond_list.number_elements();
    for (int i = 0; i < nb; i++)
    {
      _bond_list[i]->set_non_aromatic();
    }
  }*/

  if (NULL != _ring_membership || IW_NRINGS_NOT_COMPUTED != _nrings || _number_sssr_rings > 0)
    _invalidate_ring_info();

  _number_sssr_rings = IW_NRINGS_NOT_COMPUTED;

  return 1;
}

int
Molecule::_invalidate_for_changed_isotope()
{
  _smiles_information.invalidate();

  if (! include_isotopic_information_in_unique_smiles())
    return 1;

  _symmetry_class_and_canonical_rank.invalidate();

  return 1;
}

int
Molecule::invalidate_canonical_ordering_information()
{
  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  return 1;
}

int
Molecule::_ok_ring_info() const
{
  if (IW_NRINGS_NOT_COMPUTED == _nrings)
  {
    if (0 == _sssr_rings.number_elements())
      return 1;
    else
      return 0;
  }

  if (_sssr_rings.number_elements() == _nrings)     // which may include the case of _nrings = 0
    return 1;

// If no ring determinations have yet been made, that's OK.

  if (0 == _sssr_rings.number_elements())
    return 1;

  if (perceive_sssr_rings() && _sssr_rings.number_elements() > _nrings)
  {
    cerr << "Molecule::_ok_ring_info:too many SSSR rings " << _sssr_rings.number_elements() << " expect " << _nrings << endl;
    for (int i = 0; i < _sssr_rings.number_elements(); i++)
    {
      cerr << *(_sssr_rings[i]) << endl;
    }
    return 0;
  }

// Not sure what to do with the esssr case. _nrings is based on the sssr formula

  if (! perceive_sssr_rings() && _sssr_rings.number_elements() > _nrings)
    return 1;

  return 1;    // for now, fix later
  
//if (_experimental_sssr_rings.number_elements() + _raw_rings.number_elements() != _nrings)
//  return 0;

  return 0;
}

int
Molecule::_invalidate_ring_info()
{
  _nrings = IW_NRINGS_NOT_COMPUTED;
  _number_sssr_rings = IW_NRINGS_NOT_COMPUTED;
  if (NULL != _ring_membership)
  {
    delete [] _ring_membership;
    _ring_membership = NULL;
  }

  _sssr_rings.resize(0);
  _raw_rings.resize(0);
  _non_sssr_rings.resize(0);

  _experimental_raw_rings.resize(0);
  _experimental_sssr_rings.resize(0);

  if (invalidate_bond_list_ring_info_during_invalidate_ring_info)
    _bond_list.invalidate_ring_info();

  if (NULL != _aromaticity)
  {
    delete [] _aromaticity;
    _aromaticity = NULL;
  }

  return 1;
}

/*int
Molecule::_invalidate_ring_aromaticity_info()
{
  for (int i = 0; i < _sssr_rings.number_elements(); i++)
  {
    Ring * r = _sssr_rings[i];
    r->set_aromaticity_to_not_determined();
  }

  return 1;
}*/

/*
  Function to return a const pointer to a molecule's name
*/

const char *
Molecule::molecule_name() const
{
  assert(ok());

  return ((IWString &)(_molecule_name)).chars();
//return _molecule_name.chars();
}

const IWString &
Molecule::name() const
{
  assert(ok());

  return _molecule_name;
}

int
Molecule::natoms() const
{
  assert(ok());

  return _number_elements;
}

/*
  Note that we do not include implicit hydrogens (even if 1 == z)
*/

int
Molecule::natoms (atomic_number_t z) const
{
  assert(ok());
  assert(z >= 0);

  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
    if (z == _things[i]->atomic_number())
      rc++;

  return rc;
}

int
Molecule::natoms (const Element *e) const
{
  assert(ok());
  assert(OK_ELEMENT(e));

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (e == _things[i]->element())
      rc++;
  }

  return rc;
}

int
Molecule::natoms (const char *s) const
{
  assert(ok());
  assert(NULL != s);

  const Element *e;

  if (NULL == (e = get_element_from_symbol_no_case_conversion(s)))
  {
    cerr << "molecule::natoms: unrecognised element '" << s << "\n";
    return -1;
  }

  return natoms(e);
}

/*
  Often we need charges copied from one molecule to another
*/

int
Molecule::copy_charges (const Molecule & m2)
{
  assert(ok());
  assert(m2.ok());

  if (_number_elements != m2._number_elements)
  {
    cerr << "molecule::copy_charges: atom count mismatch " << _number_elements << " vs " <<
            m2._number_elements << "\n";
    return 0;
  }

  if (! m2.has_charges())
  {
    if (NULL != _charges)
      invalidate_charges();

    return 1;
  }

  if (! has_charges())
    allocate_charges();

  *_charges = *(m2._charges);

  return 1;
}

void
Molecule::invalidate_charges()
{
  if (NULL != _charges)
  {
    delete _charges;
    _charges = NULL;
  }

  return;
}

/*
  produce the vector which goes from atom N1 to atom N2
*/

int
Molecule::vector_between_atoms (atom_number_t n1,
                                atom_number_t n2,
                                Coordinates & v) const
{
  assert(ok_2_atoms(n1, n2));

  const Atom * a1 = _things[n1];
  const Atom * a2 = _things[n2];

  v = *a2;
  v -= *a1;

  return 1;
}

void
Molecule::_standardise_name()
{
  if (0 == _molecule_name.length())
    return;

  _molecule_name.strip_leading_blanks();
  _molecule_name.strip_trailing_blanks();
  
  return;
}

void
Molecule::set_name (const char *new_name)
{
  assert(ok());

  if (NULL == new_name)
  {
    _molecule_name = "";
    return;
  }

  _molecule_name = new_name;

  _standardise_name();

  return;
}

void
Molecule::set_name (const char * new_name, int lens)
{
  _molecule_name.set(new_name, lens);

  _standardise_name();

  return;
}

void
Molecule::set_name (const IWString & new_name)
{
  _molecule_name = new_name;

  _standardise_name();

  return;
}

void
Molecule::append_to_name (const IWString & zextra)
{
  _molecule_name += zextra;

  _standardise_name();

  return;
}

static void
append_formula_symbol (IWString & formula,
                       const const_IWSubstring & symbol,
                       int count)
{
  formula += symbol;
  if (count > 1)
    formula << count;

  return;
}

IWString
Molecule::molecular_formula() 
{
  IWString f;

  (void) molecular_formula(f);

  return f;
}

/*
  For things line molecular weight determinations, formula and exact mass, 
  having a count of each element type is good.
*/

void
Molecule::_compute_element_count (int * element_count,
                                  int & highest_atomic_number,
                                  int & isotopes_present,
                                  int & non_periodic_table_atoms_present) const
{
  set_vector(element_count, HIGHEST_ATOMIC_NUMBER + 1, 0);

  non_periodic_table_atoms_present = 0;
  isotopes_present = 0;

  highest_atomic_number = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];    // not const

    atomic_number_t z = a->atomic_number();

    if (z < 0)
    {
      non_periodic_table_atoms_present++;
      continue;
    }

    if (a->isotope())
      isotopes_present++;

    element_count[z]++;

    if (z > highest_atomic_number)
      highest_atomic_number = z;

    element_count[1] += a->implicit_hydrogens();
  }

  return;
}

void
Molecule::_compute_element_count (int * element_count,
                                  const int * include_atom,
                                  int & highest_atomic_number,
                                  int & isotopes_present,
                                  int & non_periodic_table_atoms_present) const
{
  set_vector(element_count, HIGHEST_ATOMIC_NUMBER + 1, 0);

  non_periodic_table_atoms_present = 0;
  isotopes_present = 0;

  highest_atomic_number = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (! include_atom[i])
      continue;

    Atom * a = _things[i];

    atomic_number_t z = a->atomic_number();

    if (z < 0)
    {
      non_periodic_table_atoms_present++;
      continue;
    }

    if (a->isotope())
      isotopes_present++;

    element_count[z]++;

    if (z > highest_atomic_number)
      highest_atomic_number = z;

    element_count[1] += a->implicit_hydrogens();
  }

  return;
}

/*
  Only include those atoms for which FLAG == ATOM_FLAG[i]
*/

void
Molecule::_compute_element_count (int * element_count,
                                  const int * atom_flag,
                                  int flag,
                                  int & highest_atomic_number,
                                  int & isotopes_present,
                                  int & non_periodic_table_atoms_present) const
{
  set_vector(element_count, HIGHEST_ATOMIC_NUMBER + 1, 0);

  non_periodic_table_atoms_present = 0;
  isotopes_present = 0;

  highest_atomic_number = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (flag != atom_flag[i])
      continue;

    Atom * a = _things[i];

    atomic_number_t z = a->atomic_number();

    if (z < 0)
    {
      non_periodic_table_atoms_present++;
      continue;
    }

    if (a->isotope())
      isotopes_present++;

    element_count[z]++;

    if (z > highest_atomic_number)
      highest_atomic_number = z;

    element_count[1] += a->implicit_hydrogens();
  }

  return;
}

/*
  Several methods need an array for holding element counts.
  No thread safety here!
*/

static int element_count[HIGHEST_ATOMIC_NUMBER + 1];

int
Molecule::molecular_formula (IWString & f) const
{
  f = "";

  if (0 == _number_elements)
    return 1;

  if (f.nchars() < 100)     // should be wide enough!
    f.resize(100);

  int highest_atomic_number, non_periodic_table_atoms_present, isotopes_present;
  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  int atoms_counted = 0;

  if (element_count[6])
  {
    append_formula_symbol(f, "C", element_count[6]);
    atoms_counted += element_count[6];
  }

  if (element_count[7])
  {
    append_formula_symbol(f, "N", element_count[7]);
    atoms_counted += element_count[7];
  }

  if (element_count[8])
  {
    append_formula_symbol(f, "O", element_count[8]);
    atoms_counted += element_count[8];
  }

  if (element_count[15])
  {
    append_formula_symbol(f, "P", element_count[15]);
    atoms_counted += element_count[15];
  }

  if (element_count[16])
  {
    append_formula_symbol(f, "S", element_count[16]);
    atoms_counted += element_count[16];
  }

  if (element_count[9])
  {
    append_formula_symbol(f, "F", element_count[9]);
    atoms_counted += element_count[9];
  }

  if (element_count[17])
  {
    append_formula_symbol(f, "Cl", element_count[17]);
    atoms_counted += element_count[17];
  }

  if (element_count[35])
  {
    append_formula_symbol(f, "Br", element_count[35]);
    atoms_counted += element_count[35];
  }

  if (element_count[53])
  {
    append_formula_symbol(f, "I", element_count[53]);
    atoms_counted += element_count[53];
  }

  if (element_count[1])
  {
    append_formula_symbol(f, "H", element_count[1]);
    atoms_counted += element_count[1];
  }

  if (atoms_counted == _number_elements)
    return 1;

// Now we have to do all the other periodic table elements - ignore the others

  for (int i = 0; i <= highest_atomic_number; i++)
  {
    if (1 == i || 6 == i || 7 == i || 8 == i || 9 == i || 15 == i || 16 == i || 17 == i || 35 == i || 53 == i)
      continue;

    int j = element_count[i];
    if (j)
    {
      const Element * e = get_element_from_atomic_number(i);
      assert(NULL != e);

      append_formula_symbol(f, e->symbol(), j);

      atoms_counted += j;

      if (atoms_counted == _number_elements)
        return 1;
    }
  }

  return 1;
}

/*
  The elements ordered in alphabetic order by symbol
*/

static int alphabetic_element_symbol_order [] = {
                    0,     /*  *   0  */
                   89,     /*  Ac  1  */
                   47,     /*  Ag  2  */
                   13,     /*  Al  3  */
                   95,     /*  Am  4  */
                   18,     /*  Ar  5  */
                   33,     /*  As  6  */
                   85,     /*  At  7  */
                   79,     /*  Au  8  */
                    5,     /*   B  9  */
                   56,     /*  Ba  10 */
                    4,     /*  Be  11 */
                  107,     /*  Bh  12 */
                   83,     /*  Bi  13 */
                   97,     /*  Bk  14 */
                   35,     /*  Br  15 */
                    6,     /*   C  16 */
                   20,     /*  Ca  17 */
                   48,     /*  Cd  18 */
                   58,     /*  Ce  19 */
                   98,     /*  Cf  20 */
                   17,     /*  Cl  21 */
                   96,     /*  Cm  22 */
                  112,     /*  Cn  23 */
                   27,     /*  Co  24 */
                   24,     /*  Cr  25 */
                   55,     /*  Cs  26 */
                   29,     /*  Cu  27 */
                  105,     /*  Db  28 */
                  110,     /*  Ds  29 */
                   66,     /*  Dy  30 */
                   68,     /*  Er  31 */
                   99,     /*  Es  32 */
                   63,     /*  Eu  33 */
                    9,     /*   F  34 */
                  114,     /*  Fl  35 */
                  100,     /*  Fm  36 */
                   87,     /*  Fr  37 */
                   31,     /*  Ga  38 */
                   64,     /*  Gd  39 */
                   32,     /*  Ge  40 */
                    1,     /*   H  41 */
                    2,     /*  He  42 */
                   72,     /*  Hf  43 */
                   80,     /*  Hg  44 */
                   67,     /*  Ho  45 */
                  108,     /*  Hs  46 */
                   53,     /*   I  47 */
                   49,     /*  In  48 */
                   77,     /*  Ir  49 */
                   26,     /*  Fe  50 */
                   19,     /*   K  51 */
                   36,     /*  Kr  52 */
                   57,     /*  La  53 */
                    3,     /*  Li  54 */
                  103,     /*  Lr  55 */
                   71,     /*  Lu  56 */
                  116,     /*  Lv  57 */
                  115,     /*  Mc  58 */
                  101,     /*  Md  59 */
                   12,     /*  Mg  60 */
                   25,     /*  Mn  61 */
                   42,     /*  Mo  62 */
                  109,     /*  Mt  63 */
                    7,     /*   N  64 */
                   11,     /*  Na  65 */
                   41,     /*  Nb  66 */
                   60,     /*  Nd  67 */
                   10,     /*  Ne  68 */
                  113,     /*  Nh  69 */
                   28,     /*  Ni  70 */
                   93,     /*  Np  71 */
                  102,     /*  No  72 */
                    8,     /*   O  73 */
                  118,     /*  Og  74 */
                   76,     /*  Os  75 */
                   15,     /*   P  76 */
                   91,     /*  Pa  77 */
                   82,     /*  Pb  78 */
                   46,     /*  Pd  79 */
                   84,     /*  Po  80 */
                   61,     /*  Pm  81 */
                   59,     /*  Pr  82 */
                   78,     /*  Pt  83 */
                   94,     /*  Pu  84 */
                   88,     /*  Ra  85 */
                   37,     /*  Rb  86 */
                   75,     /*  Re  87 */
                  104,     /*  Rf  88 */
                  111,     /*  Rg  89 */
                   45,     /*  Rh  90 */
                   86,     /*  Rn  91 */
                   44,     /*  Ru  92 */
                   16,     /*   S  93 */
                   51,     /*  Sb  94 */
                   21,     /*  Sc  95 */
                   34,     /*  Se  96 */
                  106,     /*  Sg  97 */
                   14,     /*  Si  98 */
                   62,     /*  Sm  99 */
                   50,     /*  Sn  100 */
                   38,     /*  Sr  101 */
                   73,     /*  Ta  102 */
                   65,     /*  Tb  103 */
                   43,     /*  Tc  104 */
                   52,     /*  Te  105 */
                   90,     /*  Th  106 */
                   22,     /*  Ti  107 */
                   81,     /*  Tl  108 */
                   69,     /*  Tm  109 */
                  117,     /*  Ts  110 */
                   92,     /*   U  111 */
                   23,     /*   V  112*/
                   74,     /*   W  113 */
                   54,     /*  Xe  114 */
                   39,     /*   Y  115 */
                   70,     /*  Yb  116 */
                   30,     /*  Zn  117 */
                   40};    /*  Zr  118 */

int
Molecule::isis_like_molecular_formula_dot_between_fragments (IWString & f)
{
  f = "";

  if (0 == _number_elements)
    return 1;

  int nf = number_fragments();

  if (1 == nf)
    return isis_like_molecular_formula(f);

  f.make_room_for_extra_items(24 * nf);

  resizable_array_p<Molecule> fragments;
  create_components(fragments);

  for (int i = 0; i < nf; i++)
  {
    if (i > 0)
      f += '.';

    IWString tmp;
    fragments[i]->isis_like_molecular_formula(tmp);

    f += tmp;
  }

  return 1;
}

int
Molecule::isis_like_molecular_formula (IWString & f)
{
  f = "";
  if (0 == _number_elements)
    return 1;

  f.make_room_for_extra_items(32);

  int highest_atomic_number, isotopes_present, non_periodic_table_atoms_present;     // not used here

  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  if (element_count[6])
    append_formula_symbol(f, "C", element_count[6]);

  if (element_count[1])
    append_formula_symbol(f, "H", element_count[1]);

  int completed = element_count[6];    // the number of atoms completed. Note that explicit Hydrogens are not counted.

// Loop through all the other elements in correct order

  for (int i = 0; i <= HIGHEST_ATOMIC_NUMBER && completed < _number_elements; i++)
  {
    int j = alphabetic_element_symbol_order[i];
    if (j < 0)
      continue;

    if (6 == j || 1 == j)     // did those above
      continue;

    if (element_count[j])
    {
      const Element * e = get_element_from_atomic_number(j);
      assert(NULL != e);

      append_formula_symbol(f, e->symbol(), element_count[j]);

      completed += element_count[j];
    }
  }

// Don't forget any non-periodic table elements. We don't handle multiple instances gracefully
  
  for (int i = 0; i < _number_elements && completed < _number_elements; i++)
  {
    const Element * e = _things[i]->element();

    if (e->is_in_periodic_table())
      continue;

    append_formula_symbol(f, e->symbol(), 1);
    completed++;
  }

  return 1;
}

static int
all_carbon_atoms (const Molecule & m,
                  const Set_of_Atoms & s)
{
  int n = s.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (6 != m.atomic_number(s[i]))
      return 0;
  }

  return 1;
}

static void
append_atomic_symbol (IWString & f,
                      const IWString & s,
                      int count)
{
  if (0 == count)
    return;

  f << s;

  if (count > 1)
    f << count;

  return;
}

/*
  We identify hydrogen atoms associated with aromatic carbon atoms and those
  associated with other carbon atoms.

  All hydrogens associated with heteroatoms just get counted with the Molecule's
  overall hcount.

  The only hydrogens that are completely safe to associate with an aromatic ring
  are those on aromatic carbon rings only.
*/

int
Molecule::formula_distinguishing_aromatic (IWString & f)
{
  f.resize_keep_storage(0);

  if (0 == _number_elements)
    return 1;

  compute_aromaticity_if_needed();

  int highest_atomic_number, non_periodic_table_atoms_present, isotopes_present;
  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  int * aromatic_carbon = new_int(_number_elements); std::unique_ptr<int[]> free_aromatic_carbon(aromatic_carbon);

  IWString arom_string;     // concatenation of aromatic ring sizes

  int nr = nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    if (! ri->is_aromatic())
      continue;

    arom_string << ri->number_elements();

    if (! all_carbon_atoms(*this, *ri))
      continue;

    ri->set_vector(aromatic_carbon, 1);
  }

  int molecular_hcount = 0;

  for (int i = 0; i <= highest_atomic_number; i++)
  {
    if (0 == element_count[i])
      continue;

    if (1 == i)
      continue;

    int arom_count = 0;
    int arom_hcount = 0;
    int aliph_count = 0;
    int aliph_hcount = 0;

    if (6 == i)    // handle carbon separately
    {
      int aromatic_hydrogen_count = 0;

      for (int j = 0; j < _number_elements; j++)
      {
        if (6 != _things[j]->atomic_number())
          continue;

        if (aromatic_carbon[j])
        {
          aromatic_hydrogen_count += hcount(j);
          arom_count++;
        }
        else if (IS_AROMATIC_ATOM(_aromaticity[j]))
        {
          arom_count++;
          molecular_hcount += hcount(j);
        }
        else
        {
          aliph_count++;
          molecular_hcount += hcount(j);
        }
      }

      if (arom_count)
      {
        f << 'c';
        if (arom_count > 1)
          f << arom_count;
        if (aromatic_hydrogen_count)
          append_atomic_symbol(f, 'H', aromatic_hydrogen_count);
      }
      if (aliph_count)
      {
        f << 'C';
        if (aliph_count > 1)
          f << aliph_count;
        if (aliph_hcount)
          append_atomic_symbol(f, 'H', aliph_hcount);
      }

      continue;
    }

//  cerr << "Before atomic number " << i << " molecular_hcount " << molecular_hcount << endl;

//  elements other than carbon

    for (int j = 0; j < _number_elements; j++)
    {
      if (i != _things[j]->atomic_number())
        continue;

      if (IS_AROMATIC_ATOM(_aromaticity[j]))
      {
        arom_count++;
        arom_hcount += implicit_hydrogens(j);
      }
      else
      {
        aliph_count++;
        aliph_hcount += implicit_hydrogens(j);
      }
    }

    const Element * e = get_element_from_atomic_number(i);

    if (arom_count)
      append_atomic_symbol(f, e->aromatic_symbol(), arom_count);
    if (aliph_count)
      append_atomic_symbol(f, e->symbol(), aliph_count);

    molecular_hcount += (arom_count + aliph_count);
//  cerr << "After atomic number " << i << " molecular_hcount " << molecular_hcount << endl;
  }

  if (molecular_hcount)
    append_atomic_symbol(f, 'H', molecular_hcount);

  if (arom_string.length())
    f << 'a' << arom_string;

  return 1;
}

/*
  The number of hydrogen atoms in a molecule
*/

int
Molecule::number_hydrogens() const
{
  assert(ok());

  return natoms("H");
}

int
Molecule::transform_atoms (const Element *efrom, const Element *eto)
{
  assert(ok());
  assert(OK_ELEMENT(efrom) && OK_ELEMENT(eto));

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
    if (efrom == _things[i]->element())
    {
      _things[i]->set_element(eto);
      rc++;
    }

  _set_modified();

  return rc;
}

charge_t
Molecule::charge_on_atom (atom_number_t i) const
{
  assert(ok_atom_number(i));

  if (! has_charges())
    return static_cast<charge_t>(0.0);

  return _charges->item(i);
}
formal_charge_t
Molecule::formal_charge (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->formal_charge();
}

int
Molecule::formal_charge() const
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i]->formal_charge();
  }

  return rc;
}

/*
  Normally adding a second bond between any two atoms is fatal.
  Sphinx had an application where they wanted to check molecules
  for syntax, so we need to be able to make this a non-fatal
  error.
*/

static int add_same_bond_twice_fatal = 1;

void
set_add_same_bond_twice_fatal (int f)
{
  add_same_bond_twice_fatal = f;
}

/*
  We need to put a bond between a pair of atoms

  The atoms must already be present in the molecule, and they must not
  already be bonded.

  if PARTIAL_MOLECULE, then this is being called from one of the input
  functions. This info is used to decide which Atom::add() function to
  call. The overloaded add() function will reset _implicit_hcount,
  which we don't want when we are reading molecules with a known
  _implicit_hcount
*/

int
Molecule::add_bond (atom_number_t a1, atom_number_t a2, bond_type_t bt, 
                    int partial_molecule)
{
  assert(ok_2_atoms(a1, a2));

#ifdef DEBUG_ADD_BOND
  cerr << "Adding bond between atoms " << a1 << " and " << a2 << endl;
#endif

  if (_things[a1]->is_bonded_to(a2))
  {
    if (display_already_bonded_error_message)
      cerr << "Molecule::add_bond: atoms " << a1 << " and " << a2 << " are already bonded\n";
    if (! add_same_bond_twice_fatal)
      return 0;

    debug_print(cerr);
    // IL, no need to abort, just skip this mol 
    // iwabort();
    return 0;
  }

  assert(OK_BOND_TYPE(bt));

  if (_bond_list.elements_allocated() < 30)
    _bond_list.resize(30);

  Bond * b = new Bond(a1, a2, bt);

//cerr << "Adding bond of type " << bt << endl;

  _bond_list.add(b);

  if (partial_molecule)
  {
    ((resizable_array<Bond *> *) _things[a1])->add(b);  // avoid overloaded function
    ((resizable_array<Bond *> *) _things[a2])->add(b);  // avoid overloaded function
  }
  else     // use the overloaded version to allow recomputation of necessary
  {
    _things[a1]->add(b);
    _things[a2]->add(b);

    _set_modified();    // not for partial molecule
  }

  if (partial_molecule)
    return 1;

// Dec 97, when doing reactions I ran into a problem with making a bond to an
// atom which has a chiral centre.

  int nc = _chiral_centres.number_elements();
  if (0 == nc)
    return 1;

  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * c = _chiral_centres[i];

    atom_number_t zatom;
    atom_number_t zother;

    if (a1 == c->a())
    {
      zatom = a1;
      zother = a2;
    }
    else if (a2 == c->a())
    {
      zatom = a2;
      zother = a1;
    }
    else         // C does not involve either A1 or A2
      continue;

#ifdef DEBUG_ADD_BOND
    cerr << "Must deal with chiral centre at atom " << zatom << endl;
#endif

    _things[zatom]->set_implicit_hydrogens_known(0);
    _things[zatom]->set_modified();

    if (_things[zatom]->ncon() > 4 || (! b->is_single_bond()))    // ZATOM cannot be a chiral centre any more
    {
      _chiral_centres.remove_item(i);
      i--;
      nc--;
      continue;
    }

    if (1 == c->implicit_hydrogen_count())
    {
      c->implicit_hydrogen_is_now_atom_number(zother);
      continue;
    }

    if (1 == c->lone_pair_count())
    {
      c->lone_pair_is_now_atom_number(zother);
      continue;
    }

    if (0 == c->number_connections_specified())   // Oct 2007. Reading a Kekule sdf with explicit Hydrogens
      continue;

    cerr << "Molecule::add_bond: very strange, atom " << zatom << " type " << _things[zatom]->atomic_symbol() << " ncon " << _things[zatom]->ncon() << endl;
    cerr << "Adding bond between " << a1 << " and " << a2 << endl;
    c->debug_print(cerr);
  }

  return 1;
}

static int
int_comparitor_shorter (const int * p1, const int * p2)
{
  assert(NULL != p1);
  assert(NULL != p2);

  if (*p1 < *p2)
    return 1;
  else if (*p1 == *p2)
    return 0;
  else
    return -1;
}

void
Molecule::remove_atom_from_charge_arrays (const atom_number_t atom_to_remove)
{
  if (NULL != _charges)
    _charges->remove_item(atom_to_remove);

  if (_atom_type)
    _atom_type->remove_item(atom_to_remove);

  return;
}

/*
  This core functionality is used by remove_atom and by remove_atoms.
  Note that it does no checking, and does not call set_modified.

  August 2015.
  Put in extra processing for dealing with Hydrogen atoms being removed,
  and the implicit Hydrogen properties of any attached atoms

  Jan 2016. If we are removing an explicit H on a non-organic,
  we need to save that info and put it back when done. Because
  of how the Atom methods behave...
*/

int
Molecule::_remove_atom (const atom_number_t atom_to_remove)
{
  const atomic_number_t z = _things[atom_to_remove]->atomic_number();

  atom_number_t explicit_hydrogen_attached_to_non_organic = INVALID_ATOM_NUMBER;
  
  const int acon = _things[atom_to_remove]->ncon();

  if (1 == z && 1 == acon)    // singly connected explicit H
  {
    const atom_number_t x = _things[atom_to_remove]->other(atom_to_remove, 0);
    if (! _things[x]->element()->organic())     // likely will not know anything about implicit hydrogens. Bug: what if the unknown value is non-zero? Kind of hard to know what is best here
      explicit_hydrogen_attached_to_non_organic = x;
    else if (_things[x]->implicit_hydrogens_known())    // organic elements will be able to recompute in hopefully all cases. Easy to imagine cases where this would be wrong. Hopefully unusual...
      _things[x]->set_implicit_hydrogens_known(0);
  }

#ifdef DEBUG_REMOVE_ATOM
  cerr << "Molecule:: removing atom " << atom_to_remove << " acon " << acon << endl;
  if (INVALID_ATOM_NUMBER != explicit_hydrogen_attached_to_non_organic)
  {
    cerr << "Our atom has " << implicit_hydrogens(explicit_hydrogen_attached_to_non_organic) << " IH\n";
    debug_print(cerr);
  }
#endif

  if (acon > 0)
  {
    _atom_being_unbonded_check_directional_bonds(atom_to_remove);

    int initial_implicit_hydrogen_count = 0;
    if (INVALID_ATOM_NUMBER != explicit_hydrogen_attached_to_non_organic)
      initial_implicit_hydrogen_count = _things[explicit_hydrogen_attached_to_non_organic]->implicit_hydrogens();

    (void) _remove_bonds_to_atom(atom_to_remove, 1);

    if (INVALID_ATOM_NUMBER != explicit_hydrogen_attached_to_non_organic)
    {
      _things[explicit_hydrogen_attached_to_non_organic]->set_implicit_hydrogens(initial_implicit_hydrogen_count+1, 1);
      _things[explicit_hydrogen_attached_to_non_organic]->set_implicit_hydrogens_known(1);
    }
  }
  else
    _bond_list.adjust_atom_numbers_for_loss_of_atom(atom_to_remove);

  remove_atom_from_charge_arrays(atom_to_remove);

// We must tell _adjust_chiral_centres.. whether or not this 
// was a hydrogen, as it handles hydrogens specially.

  _adjust_chiral_centres_for_loss_of_atom(atom_to_remove, (1 == z));

  (void) remove_item(atom_to_remove);

  return 1;
}

int
Molecule::remove_atom (atom_number_t atom_to_remove)
{
  assert(ok_atom_number(atom_to_remove));

  int rc = _remove_atom(atom_to_remove);

  _set_modified();

  return rc;
}

//#define DEBUG_REMOVE_ATOMS

int
Molecule::remove_atoms (Set_of_Atoms & atoms_to_remove)
{
  assert(ok());

  int nr = atoms_to_remove.number_elements();
  if (0 == nr)
    return 0;

  if (nr > 1)
    atoms_to_remove.sort(int_comparitor_shorter);

#ifdef DEBUG_REMOVE_ATOMS
  cerr << "Molecule::remove_atoms: will remove atoms " << atoms_to_remove << endl;
#endif

  if (! ok_index(atoms_to_remove[0]) || ! ok_index(atoms_to_remove.last_item()))
  {
    cerr << "One or more invalid atom numbers encountered. Molecule has " << _number_elements << " atoms\n";
    cerr << atoms_to_remove << endl;

    return 0;
  }

  int rc = 0;

// Note that we ordered the atoms above, so that as we remove atoms, the
// atom numbers in ATOMS_TO_REMOVE do not change.

  for (int i = 0; i < nr; i++)
  {
    rc += _remove_atom(atoms_to_remove[i]);
  }

  _set_modified();

  return rc;
}

int
Molecule::remove_atoms (const int * to_remove)
{
  assert(ok());

#ifdef DEBUG_REMOVE_ATOMS
  for (auto i = 0; i < _number_elements; ++i)
  {
    if (to_remove[i])
      cerr << "Molecule will remove atom " << i << endl;
  }
#endif

  int rc = 0;
  for (int i = _number_elements - 1; i >= 0; i--)
  {
    if (to_remove[i])
    {
      _remove_atom(i);
      rc++;
    }
  }

  if (rc)
    _set_modified();

  return rc;
}

int
Molecule::remove_many_atoms (const int * to_remove)
{
  if (_chiral_centres.number_elements())
    return remove_atoms(to_remove);

  for (auto i = 0; i < _number_elements; ++i)
  {
    _things[i]->remove_connections_to_any_of_these_atoms(to_remove);
  }

  _bond_list.remove_bonds_involving_these_atoms(to_remove, 1);

  if (NULL != _atom_type)
    _atom_type->remove_items(to_remove);

  if (NULL != _charges)
    _charges->remove_items(to_remove);

  resizable_array_p<Atom>::remove_items(to_remove);

  _set_modified();

  return 1;
}

int
Molecule::remove_fragment_containing_atom (atom_number_t zremove)
{
  int f = fragment_membership(zremove);

  if (1 == number_fragments())
  {
    cerr << "Molecule::remove_fragment_containing_atom: molecule contains only one fragment\n";
    return 0;
  }

  return delete_fragment(f);
}

int
Molecule::number_isotopic_atoms() const
{
  assert(ok());

  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom *a = _things[i];
    if (a->is_isotope())
      rc++;
  }

  return rc;
}

int
Molecule::number_isotopic_atoms (int iso) const
{
  assert(ok());

  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom *a = _things[i];
    if (iso == a->isotope())
      rc++;
  }

  return rc;
}

/*
  Special function for removing an isotopic specification
*/

int
Molecule::_set_isotope_zero(atom_number_t zatom)
{
  Atom * a = _things[zatom];

  if (0 == a->isotope())
    return 1;

  a->set_isotope(0);
  if (a->implicit_hydrogens_known())
  {
    int ih;
    a->compute_implicit_hydrogens(ih);
    if (ih == a->implicit_hydrogens())
      a->set_implicit_hydrogens_known(0);
  }

  return 1;
}

int
Molecule::transform_to_non_isotopic_form()
{
  assert(ok());

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];
    if (a->isotope())
    {
      _set_isotope_zero(i);

      rc++;
    }
  }

  if (rc)
    _invalidate_for_changed_isotope();

  return rc;
}

template <typename T> 
int
Molecule::set_isotopes(const T * iso)
{
  assert(NULL != iso);

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (iso[i] < 0)    // never needed if T is unsigned
      continue;

    if (iso[i] > 0)
      _things[i]->set_isotope(iso[i]);
    else
      _set_isotope_zero(i);

    rc++;
  }

   if (rc)
    _invalidate_for_changed_isotope();

  return rc;
}

template int Molecule::set_isotopes(const int *);
template int Molecule::set_isotopes(const unsigned int *);

int
Molecule::unset_isotopes (const int * iso)
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (iso[i] > 0)
    {
      _set_isotope_zero(i);
      rc++;
    }
  }

  if (rc)
    _invalidate_for_changed_isotope();

  return rc;
}

int
Molecule::set_isotope (const Set_of_Atoms & s,
                       int iso)
{
  int n = s.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = s[i];

    if (iso > 0)
      _things[j]->set_isotope(iso);
    else if (0 == iso)
      _set_isotope_zero(j);
  }

  _invalidate_for_changed_isotope();

  return 1;
}

void
Molecule::get_isotopes (int * iso) const
{
  assert(NULL != iso);

  for (int i = 0; i < _number_elements; i++)
  {
    iso[i] = _things[i]->isotope();
  }

  return;
}

int
Molecule::set_isotope (atom_number_t a, int iso)
{
  assert(ok_atom_number(a));

  if (iso > 0)
    _things[a]->set_isotope(iso);
  else if (0 == iso)
    _set_isotope_zero(a);

  _invalidate_for_changed_isotope();

  return 1;
}

int
Molecule::set_userAtomType (atom_number_t a, int atomType)
{
  assert(ok_atom_number(a));

  _things[a]->set_userAtomType(atomType);
  
  return 1;
}


void
Molecule::set_isotope_to_atom_number_no_perturb_canonical_ordering()
{
  for (int i = 0; i < _number_elements; ++i)
  {
    _things[i]->set_isotope(i);
  }

  _smiles_information.invalidate();

  return;
}

int
Molecule::set_isotope_no_perturb_canonical_ordering (atom_number_t a,
                                int iso)
{
  assert(ok_atom_number(a));

  _things[a]->set_isotope(iso);

  _smiles_information.invalidate();

  return 1;
}

int
Molecule::isotope (atom_number_t a) const
{
  assert(ok_atom_number(a));

  return _things[a]->isotope();
}

int
Molecule::userAtomType (atom_number_t a) const
{
  assert(ok_atom_number(a));

  return _things[a]->userAtomType();
}

int
Molecule::maximum_isotope( ) const
{
  int maxi = _things[0]->isotope();

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->isotope() > maxi)
      maxi = _things[i]->isotope();
  }

  return maxi;
}

int
Molecule::increment_isotope (atom_number_t zatom,
                             int incr) 
{
  assert(ok_atom_number(zatom));

//int current_isotope = _things[zatom]->isotope();

  int new_isotope = _things[zatom]->isotope() + incr;

  if (new_isotope < 0)
  {
    cerr << "Molecule::increment_isotope:out of range, from " << _things[zatom]->isotope() << " increment " << incr << endl;
    return 0;
  }

  _things[zatom]->set_isotope(new_isotope);

  _invalidate_for_changed_isotope();

  return 1;
}

static int
issue_non_periodic_table_molecular_weight_warning = 1;

void 
set_issue_non_periodic_table_molecular_weight_warning (int s)
{
  issue_non_periodic_table_molecular_weight_warning = s;
}

molecular_weight_t
Molecule::molecular_weight() const
{
  assert(ok());

  int highest_atomic_number = 0;
  int non_periodic_table_atoms_present;
  int isotopes_present;
  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  if (non_periodic_table_atoms_present)
  {
    if (issue_non_periodic_table_molecular_weight_warning)
      cerr << "Molecule::molecular_weight: " << non_periodic_table_atoms_present << " non periodic table elements present\n";
    return static_cast<molecular_weight_t>(0.0);
  }

  if (isotopes_present)
  {
    cerr << "Molecule::molecular_weight: " << isotopes_present << " isotopic atoms present\n";
    return static_cast<molecular_weight_t>(0.0);
  }

  molecular_weight_t rc = static_cast<molecular_weight_t>(0.0);

  int atoms_encountered = 0;

  for (int i = 0; i <= highest_atomic_number; i++)
  {
    if (0 == element_count[i])
      continue;

    const Element * e = get_element_from_atomic_number(i);

    atoms_encountered += element_count[i];

    rc += static_cast<molecular_weight_t>((element_count[i]) * e->atomic_mass());
  }

  return rc;
}

molecular_weight_t
Molecule::molecular_weight_ignore_isotopes() const
{
  assert(ok());

  int highest_atomic_number = 0;
  int non_periodic_table_atoms_present;
  int isotopes_present;
  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  if (non_periodic_table_atoms_present)
  {
    if (issue_non_periodic_table_molecular_weight_warning)
      cerr << "Molecule::molecular_weight: " << non_periodic_table_atoms_present << " non periodic table elements present\n";
    return static_cast<molecular_weight_t>(0.0);
  }

  molecular_weight_t rc = static_cast<molecular_weight_t>(0.0);

  int atoms_encountered = 0;

  for (int i = 0; i <= highest_atomic_number; i++)
  {
    if (0 == element_count[i])
      continue;

    const Element * e = get_element_from_atomic_number(i);

    atoms_encountered += element_count[i];

    rc += static_cast<molecular_weight_t>((element_count[i]) * e->atomic_mass());
//  cerr << "got " << element_count[i] << " of atomic number " << i << " sum now " << rc << endl;
  }

  return rc;
}

static const Element * hydrogen = NULL;

molecular_weight_t
Molecule::molecular_weight_count_isotopes() const
{
  assert(ok());

  if (NULL == hydrogen)
    hydrogen = get_element_from_atomic_number(1);

  molecular_weight_t rc = static_cast<molecular_weight_t>(0.0);

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = const_cast<Atom *>(_things[i]);    // loss of const for the implicit_hydrogen computation

    const Element * e = a->element();

    if (! e->is_in_periodic_table())
    {
      cerr << "Molecule::molecular_weight_count_isotopes:non periodic table elements present\n";
      return static_cast<molecular_weight_t>(0.0);
    }

    if (a->isotope())
      rc += a->isotope();
    else
      rc += e->atomic_mass() + a->implicit_hydrogens() * hydrogen->atomic_mass();
  }

  return rc;
}

int
Molecule::exact_mass (exact_mass_t & zresult) const
{
  int highest_atomic_number = 0;
  int non_periodic_table_atoms_present = 0;
  int isotopes_present = 0;

  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  return _exact_mass(element_count, highest_atomic_number, non_periodic_table_atoms_present, zresult);
}

exact_mass_t
Molecule::exact_mass() const
{
  exact_mass_t rc;
  if (! exact_mass(rc))
    return static_cast<exact_mass_t>(0.0);

  return rc;
}

int
Molecule::exact_mass (const int * include_atom, exact_mass_t & zresult) const
{
  int highest_atomic_number = 0;
  int non_periodic_table_atoms_present;
  int isotopes_present;

  _compute_element_count(element_count, include_atom, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  return _exact_mass(element_count, highest_atomic_number, non_periodic_table_atoms_present, zresult);
}

int
Molecule::exact_mass (const int * atom_flag, int flag, exact_mass_t & zresult) const
{
  int highest_atomic_number = 0;
  int non_periodic_table_atoms_present;
  int isotopes_present;

  _compute_element_count(element_count, atom_flag, flag, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  return _exact_mass(element_count, highest_atomic_number, non_periodic_table_atoms_present, zresult);
}

int
Molecule::_exact_mass (const int * element_count,
                       int highest_atomic_number,
                       int non_periodic_table_atoms_present,
                       exact_mass_t & zresult) const
{
  zresult = 0.0;

  int rc = 1;      // let's hope the result is OK

  double tmp = 0.0;     // maximum accuracy

  for (int i = 0; i <= highest_atomic_number; i++)
  {
    if (0 == element_count[i])    // none of this type present
      continue;

    const Element * e = get_element_from_atomic_number(i);

    exact_mass_t x = e->exact_mass();
    if (x > 0.0)
      tmp += element_count[i] * x;
    else
      rc = 0;
  }

  zresult = static_cast <exact_mass_t>(tmp);

  if (non_periodic_table_atoms_present)
    cerr << "Molecule::_exact_mass: " << non_periodic_table_atoms_present << " non periodic table atoms present\n";

  return rc;
}

int
Molecule::molecular_weight (const Molecular_Weight_Control & mwc,
                            Molecular_Weight_Calculation_Result & mwcr) const
{
  mwcr.reset();

  double amw = 0.0;

  int ih = 0;     // implicit Hydrogens

  int nc = 0;    // carbon
  int nn = 0;    // nitrogen
  int no = 0;    // oxygen

  for (int i = 0; i < _number_elements; ++i)
  {
    const int iso = _things[i]->isotope();

    if (iso > 0)
    {
      mwcr._isotopes_found++;
      if (mwc._ignore_isotopes)
        ;
      else
        amw += static_cast<double>(iso);

      continue;
    }

    const Element * e = _things[i]->element();

    if (e->is_in_periodic_table())
    {
      const atomic_number_t z = e->atomic_number();
      if (6 == z)
        nc++;
      else if (7 == z)
        nn++;
      else if (8 == z)
        no++;
      else
        amw += e->atomic_mass();

      ih += _things[i]->implicit_hydrogens();
    }
    else
    {
      mwcr._non_periodic_table_elements_found++;
      if (! mwc._ignore_non_periodic_table_elements)
        return 0;
    }
  }

  if (ih && ! mwc._ignore_hydrogens)
  {
    const Element * h = get_element_from_atomic_number(1);
    amw += ih * h->atomic_mass();
  }

  if (nc)
  {
    const double amw_carbon = get_element_from_atomic_number(6)->atomic_mass();
    amw += nc * amw_carbon;
  }
  if (nn)
  {
    const double amw_nitrogen = get_element_from_atomic_number(7)->atomic_mass();
    amw += nn * amw_nitrogen;
  }
  if (no)
  {
    const double amw_oxygen = get_element_from_atomic_number(8)->atomic_mass();
    amw += no * amw_oxygen;
  }

  mwcr._amw = amw;

  return 1;
}

Molecular_Weight_Control::Molecular_Weight_Control()
{
  _ignore_isotopes = false;
  _ignore_non_periodic_table_elements = false;
  _ignore_hydrogens = false;

  return;
}

void
Molecular_Weight_Calculation_Result::_default_values()
{
  _isotopes_found = 0;
  _non_periodic_table_elements_found = 0;
  _amw = 0.0;

  return;
}

Molecular_Weight_Calculation_Result::Molecular_Weight_Calculation_Result()
{
  _default_values();
}

void
Molecular_Weight_Calculation_Result::reset ()
{
  _default_values();

  return;
}

int
Molecule::number_different_elements() const
{
  int highest_atomic_number = 0;
  int non_periodic_table_atoms_present;
  int isotopes_present;

  _compute_element_count(element_count, highest_atomic_number, isotopes_present, non_periodic_table_atoms_present);

  int rc = 0;
  for (int i = 0; i <= highest_atomic_number; i++)
  {
    if (element_count[i])
      rc++;
  }

  return rc;
}

atomic_mass_t
Molecule::atomic_mass (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->element()->atomic_mass();
}

void
Molecule::translate_atoms (coord_t x, coord_t y, coord_t z)
{
  assert(ok());

  for (int i = 0; i < _number_elements; i++)
  {
    Coordinates * t = _things[i];
    t->add(x, y, z);
  }

  return;
}

void
Molecule::translate_atoms (coord_t x, coord_t y, coord_t z,
                           const Set_of_Atoms & atoms_to_move)
{
  assert(ok());

  int moving_atoms = atoms_to_move.number_elements();
  for (int i = 0; i < moving_atoms; i++)
  {
    atom_number_t atom_to_move = atoms_to_move[i];
    assert(ok_index(atom_to_move));

    Coordinates * t = _things[atom_to_move];

    t->add(x, y, z);
  }

  return;
}

void
Molecule::translate_atoms (const Coordinates & whereto,
                           const Set_of_Atoms & atoms_to_move)
{
  assert(ok());

  int moving_atoms = atoms_to_move.number_elements();

  for (int i = 0; i < moving_atoms; i++)
  {
    atom_number_t atom_to_move = atoms_to_move[i];
    assert(ok_index(atom_to_move));

    Coordinates * t = _things[atom_to_move];
    *t += whereto;
  }

  return;
}

void
Molecule::translate_atoms (const Coordinates & whereto)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->translate(whereto);
  }

  return;
}

void
Molecule::translate_atoms (const Coordinates & whereto,
                           const int * to_move,
                           int flag)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (flag == to_move[i])
      _things[i]->translate(whereto);
  }

  return;
}

int
Molecule::rotate_atoms (const Coordinates & axis, angle_t theta)
{
  assert(ok());

  if (static_cast<angle_t>(0.0) == theta)
    return 1;

  if (0 == _number_elements)
    return 0;

// The direction cosines for the vector
 
  const coord_t dc1 = axis.x();
  const coord_t dc2 = axis.y();
  const coord_t dc3 = axis.z();
//cerr << "Dc's are " << dc1 << "," << dc2 << "," << dc3 << " sum = " <<
//     dc1 * dc1 + dc2 * dc2 + dc3 * dc3 << ", angle " << theta << endl;
  
// Rather than deal properly with a matrix, just use individual variables

  coord_t rotmat11 = static_cast<coord_t>(cos(theta) + dc1 * dc1 * (1.0 - cos(theta)) );
  coord_t rotmat12 = static_cast<coord_t>(dc1 * dc2 * (1.0 - cos(theta)) - dc3 * sin(theta) );
  coord_t rotmat13 = static_cast<coord_t>(dc1 * dc3 * (1.0 - cos(theta)) + dc2 * sin(theta) );
  coord_t rotmat21 = static_cast<coord_t>(dc1 * dc2 * (1.0 - cos(theta)) + dc3 * sin(theta) );
  coord_t rotmat22 = static_cast<coord_t>(cos(theta) + dc2 * dc2 * (1.0 - cos(theta)) );
  coord_t rotmat23 = static_cast<coord_t>(dc2 * dc3 * (1.0 - cos(theta)) - dc1 * sin(theta) );
  coord_t rotmat31 = static_cast<coord_t>(dc3 * dc1 * (1.0 - cos(theta)) - dc2 * sin(theta) );
  coord_t rotmat32 = static_cast<coord_t>(dc3 * dc2 * (1.0 - cos(theta)) + dc1 * sin(theta) );
  coord_t rotmat33 = static_cast<coord_t>(cos(theta) + dc3 * dc3 * (1.0 - cos(theta)) );

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];
  
    coord_t x0 = a->x();
    coord_t y0 = a->y();
    coord_t z0 = a->z();

    coord_t xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    coord_t yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    coord_t zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    a->setxyz(xx, yy, zz);
  }

  return 1;
}

//#define DEBUG_ROTATE_ATOMS

template <typename T>
int
Molecule::rotate_atoms (const Space_Vector<T> & axis, T theta,
                        const Set_of_Atoms & atoms_to_move)
{
  assert(ok());

  if (static_cast<T>(0.0) == theta)
    return 1;

  int moving_atoms = atoms_to_move.number_elements();
  if (0 == moving_atoms)
    return 0;

// The direction cosines for the vector
 
  const T dc1 = axis.x();
  const T dc2 = axis.y();
  const T dc3 = axis.z();
  
#ifdef DEBUG_ROTATE_ATOMS
  cerr << "Dc's are " << dc1 << "," << dc2 << "," << dc3 << " sum = " <<
       dc1 * dc1 + dc2 * dc2 + dc3 * dc3 << "\n";
#endif

// Rather than deal properly with a matrix, just use individual variables

  T rotmat11 = static_cast<T>(cos(theta) + dc1 * dc1 * (1.0 - cos(theta)) );
  T rotmat12 = static_cast<T>(dc1 * dc2 * (1.0 - cos(theta)) - dc3 * sin(theta) );
  T rotmat13 = static_cast<T>(dc1 * dc3 * (1.0 - cos(theta)) + dc2 * sin(theta) );
  T rotmat21 = static_cast<T>(dc1 * dc2 * (1.0 - cos(theta)) + dc3 * sin(theta) );
  T rotmat22 = static_cast<T>(cos(theta) + dc2 * dc2 * (1.0 - cos(theta)) );
  T rotmat23 = static_cast<T>(dc2 * dc3 * (1.0 - cos(theta)) - dc1 * sin(theta) );
  T rotmat31 = static_cast<T>(dc3 * dc1 * (1.0 - cos(theta)) - dc2 * sin(theta) );
  T rotmat32 = static_cast<T>(dc3 * dc2 * (1.0 - cos(theta)) + dc1 * sin(theta) );
  T rotmat33 = static_cast<T>(cos(theta) + dc3 * dc3 * (1.0 - cos(theta)) );

#ifdef DEBUG_ROTATE_ATOMS
  cerr << "Rotation matrix is\n" << rotmat11 << " " << rotmat12 << " " << rotmat13 << "\n";
  cerr << rotmat21 << " " << rotmat22 << " " << rotmat23 << "\n";
  cerr << rotmat31 << " " << rotmat32 << " " << rotmat33 << "\n";
#endif

  for (int i = 0; i < moving_atoms; i++)
  {
    atom_number_t j = atoms_to_move[i];

    assert(j >= 0 && j < _number_elements);

    Atom *a = _things[j];

//  cerr << "Initial coordinates for atom " << j << " " << *a << endl;
  
    T x0 = a->x();
    T y0 = a->y();
    T z0 = a->z();

    T xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    T yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    T zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    a->setxyz(static_cast<coord_t>(xx), static_cast<coord_t>(yy), static_cast<coord_t>(zz) );
//  cerr << "New coordinates for atom " << j << " " << *a << endl;
  }

  return 1;
}

template int Molecule::rotate_atoms(const Space_Vector<double> &, double, const Set_of_Atoms &);
template int Molecule::rotate_atoms(const Space_Vector<coord_t> &, angle_t, const Set_of_Atoms &);

void
Molecule::rotate_to_longest_distance_along_x (atom_number_t & left, atom_number_t & right)
{
  left  = INVALID_ATOM_NUMBER;
  right = INVALID_ATOM_NUMBER;

  coord_t maxd = static_cast<coord_t>(0.0);

  if (2 == _number_elements)
  {
    left = 0;
    right = 1;

    maxd = _things[0]->distance(*(_things[1]));
  }
  else
  {
    for (auto i = 0; i < _number_elements; i++)
    {
      const Atom * ai = _things[i];

      for (int j = i + 1; j < _number_elements; j++)
      {
        const Atom * aj = _things[j];

        coord_t d = ai->distance(*aj);

        if (d > maxd)
        {
          maxd = d;
          left = i;
          right = j;
        }
      }
    }
  }

// In order to enforce more uniform behaviour, canonicalise left and right. 

  if (_things[left]->x() > _things[right]->x())
    std::swap(left, right);

  assert(INVALID_ATOM_NUMBER != left);
  assert(INVALID_ATOM_NUMBER != right);

// Translate the molecule so that atom LEFT is at the origin

  translate_atoms(- *(_things[left]));

  assert(static_cast<coord_t>(0.0) == _things[left]->x() && static_cast<coord_t>(0.0) == _things[left]->y() && static_cast<coord_t>(0.0) == _things[left]->z());

// We now want to rotate the molecule so that RIGHT is along the X axis

  Coordinates x(static_cast<coord_t>(1.0), static_cast<coord_t>(0.0), static_cast<coord_t>(0.0));
  Coordinates r(_things[right]->x(), _things[right]->y(), _things[right]->z());
  r.normalise();
  angle_t theta = x.angle_between_unit_vectors(r);

  r.cross_product(x);
  r.normalise();

//#define DEBUG_ROTATE_TO_LOGNEST_DISTANCE
#ifdef DEBUG_ROTATE_TO_LOGNEST_DISTANCE
  cerr << "Right starts " << atom[right]->x() << ' ' << atom[right]->y() << ' ' << atom[right]->z() << " angle " << (theta * RAD2DEG) << endl;
#endif

  rotate_atoms(r, theta);

#ifdef DEBUG_ROTATE_TO_LOGNEST_DISTANCE
  cerr << "Right now    " << atom[right]->x() << ' ' << atom[right]->y() << ' ' << atom[right]->z() << " angle " << (theta * RAD2DEG) << endl;
#endif

// At this stage, the molecule should have LEFT at the origin and RIGHT somwhere along the X axis

  assert(static_cast<coord_t>(0.0) == _things[left]->x() && static_cast<coord_t>(0.0) == _things[left]->y() && static_cast<coord_t>(0.0) == _things[left]->z());
  assert(fabs(_things[right]->y()) < static_cast<coord_t>(0.02));
  assert(fabs(_things[right]->z()) < static_cast<coord_t>(0.02));

#ifdef DEBUG_ROTATE_TO_LOGNEST_DISTANCE
  report_average_position(m, "At end of rotate_2", cerr);
  cerr << "left is " << left << " (" << _things[left]->x() << ',' << _things[left]->y() << "), right " << right << " (" << _things[right]->x() << ',' << _things[right]->y() << ")\n";
#endif

  return;
}

Coordinates 
Molecule::get_coords (atom_number_t i) const
{
  assert(ok_atom_number(i));

  const Atom *a = _things[i];

  return Coordinates(*a);
}

int
Molecule::get_coords (atom_number_t i, Coordinates & v) const
{
  assert(ok_atom_number(i));

  const Atom * a = _things[i];

  v.setxyz(a->x(), a->y(), a->z());

  return 1;
}

int
Molecule::get_coords (Coordinates * c) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];

    c[i] = *a;
  }

  return _number_elements;
}

coord_t
Molecule::x (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->x();
}
coord_t
Molecule::y (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->y();
}
coord_t
Molecule::z (atom_number_t i) const
{
  assert(ok_atom_number(i));

  return _things[i]->z();
}

void
Molecule::setx (atom_number_t a, coord_t newx)
{
  assert(ok_atom_number(a));

  _things[a]->x() = newx;

  return;
}

void
Molecule::sety (atom_number_t a, coord_t newy)
{
  assert(ok_atom_number(a));

  _things[a]->y() =  newy;

  return;
}

void
Molecule::setz (atom_number_t a, coord_t newz)
{
  assert(ok_atom_number(a));

  _things[a]->z() = newz;

  return;
}

void
Molecule::setxyz (atom_number_t a, coord_t newx, coord_t newy, coord_t newz)
{
  assert(ok_atom_number(a));

  _things[a]->setxyz(newx, newy, newz);

  return;
}

void
Molecule::setxyz (const Coordinates * ca)
{
  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    const Coordinates & c = ca[i];
    a->setxyz(c.x(), c.y(), c.z());
  }

  return;
}

distance_t
Molecule::bond_length (atom_number_t a1, atom_number_t a2, int not_bonded_ok) const
{
  assert(ok_2_atoms(a1, a2));
  if (! not_bonded_ok)
  {
    assert(are_bonded(a1, a2));
  }

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];

  return aa1->distance(*aa2);
}

int
Molecule::set_bond_length (atom_number_t a1, atom_number_t a2,
                           distance_t d,
                           atom_number_t atom_to_move)
{
  if (! are_bonded(a1, a2))
  {
    cerr << "Molecule::set_bond_length: atoms " << a1 << " and " << a2 << " are not bonded\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER == atom_to_move)
    atom_to_move = a2;
  else if (atom_to_move == a1)    // swap them, we move the atoms attached to A2
  {
    a1 = a2;
    a2 = atom_to_move;
  }
  else if (atom_to_move == a2)
    ;
  else
  {
    cerr << "Molecule::set_bond_length: setting bond between " << a1 << " and " << a2 << " move " << atom_to_move << endl;
    return 0;
  }

  int * moving_atoms = new_int(_number_elements); std::unique_ptr<int[]> free_moving_atoms(moving_atoms);

  return _set_bond_length(a1, a2, d, moving_atoms);
}

//#define DEBUG_SET_BOND_LENGTH

int
Molecule::_set_bond_length (atom_number_t a1, atom_number_t a2,
                            distance_t d,
                            int * moving_atoms)
{
  moving_atoms[a1] = 2;    // special flag - if this value is encountered, in _determine_moving_atoms, we abort
  moving_atoms[a2] = 1;

  const Atom * aa2 = _things[a2];

  int acon = aa2->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = aa2->other(a2, i);
    if (j == a1)
      continue;

    if (! _determine_moving_atoms(j, moving_atoms))
    {
      cerr << "Molecule::set_bond_length:cannot identify atoms to move, atoms " << a1 << " and " << a2 << endl;
      return 0;
    }
  }

  moving_atoms[a1] = 0;    // atom a1 does not move

  Coordinates c12 = *(_things[a2]) - *(_things[a1]);

  distance_t current_distance = c12.length();

#ifdef DEBUG_SET_BOND_LENGTH
  cerr << "Setting bond between atoms " << a1 << " '" << smarts_equivalent_for_atom(a1) << "' and " << a2 << " '" << smarts_equivalent_for_atom(a2) << "'\n";
  cerr << "Current distance " << current_distance << " vector " << c12 << endl;
#endif

  if (fabs(current_distance - d) < 0.00001)    // highly unlikely
    return 1;

// Need to handle the case where the atoms are on top of each other already

  if (c12.norm() < 1.0e-03)
    c12.setxyz(1.0, 0.0, 0.0);    // random - could theoretically do better by going and looking at bonded atoms
  else
    c12.normalise();

  c12 *= (d - current_distance);

#ifdef DEBUG_SET_BOND_LENGTH
  int atoms_moving = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (moving_atoms[i])
    {
      atoms_moving++;
      cerr << "Moving atom " << i << " '" << smarts_equivalent_for_atom(i) << endl;
    }
  }

  cerr << "Moving " << atoms_moving << " atoms along " << c12 << endl;
  cerr << *(_things[a1]) << " and " << *(_things[a2]) << endl;
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == moving_atoms[i])
      continue;

    Atom * a = _things[i];

    a->setxyz(a->x() + c12.x(), a->y() + c12.y(), a->z() + c12.z());
  }

#ifdef DEBUG_SET_BOND_LENGTH
  cerr << "After moving, distance " << distance_between_atoms(a1, a2) << endl;
  cerr << *(_things[a1]) << " and " << *(_things[a2]) << endl;
#endif

  return 1;
}

angle_t
Molecule::bond_angle (atom_number_t a1, atom_number_t a2, atom_number_t a3, int not_bonded_ok) const
{
  assert(ok_3_atoms(a1, a2, a3));
  if (! not_bonded_ok)
  {
    assert(are_bonded(a1, a2));
    assert(are_bonded(a2, a3));
  }

  const Atom * aa1 = _things[a1];
  const Atom * aa3 = _things[a3];

  return _things[a2]->angle_between(*aa1, *aa3);
}

int
Molecule::remove_all(atomic_number_t to_remove)
{
  assert(ok());

  const Element * e = get_element_from_atomic_number(to_remove);

  if (NULL == e)
  {
    cerr << "Molecule::remove_all: what element is this " << to_remove << endl;
    abort();
    return 0;
  }


  return remove_all(e);
}

int
Molecule::remove_all_atoms_with_isotope (int iso)
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (iso != _things[i]->isotope())
      continue;

    remove_atom(i);
    i--;
    rc++;
  }

  return rc;
}

int
Molecule::remove_all(const Element * to_remove)
{
  assert(ok());
  assert(OK_ELEMENT(to_remove));
   
  int rc = 0;

  if (1 == to_remove->atomic_number())
    return remove_explicit_hydrogens();

  if (_number_elements < 100 || _chiral_centres.number_elements())
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (_things[i]->element() != to_remove)
        continue;

      rc++;
      (void) remove_atom(i);
      i--;
    }
  }
  else
  {
    int * to_be_removed = new_int(_number_elements); std::unique_ptr<int[]> free_to_be_removed(to_be_removed);

    for (auto i = 0; i < _number_elements; ++i)
    {
      if (_things[i]->element() == to_remove)
        to_be_removed[i] = 1;
    }

    return remove_atoms(to_be_removed);
  }

  return rc;
}

/*
  The reason this is special is if there are non-organic atoms attached to multiple explicit hydrogens.
  We need to get the H count of the remaining atom correct.

  S1C(=CC(=C1)B)C=O PBCHM3255522

  particular case in point
*/

int
Molecule::remove_explicit_hydrogens ()
{
  if (0 == _number_elements)
    return 0;

  if (1 == _number_elements && 1 == _things[0]->atomic_number())    // do not disappear H
    return 0;

  int * hcount = new_int(_number_elements + _number_elements + _number_elements);std::unique_ptr<int[]> free_hcount(hcount);
  int * is_hydrogen = hcount + _number_elements;
  int * xref = hcount + _number_elements + _number_elements;
  std::fill_n(xref, _number_elements, -1);

//#define DEBUG_REMOVE_EXPLICIT_HYDROGENS
#ifdef DEBUG_REMOVE_EXPLICIT_HYDROGENS
  cerr << "Molecule::remove_explicit_hydrogens molecule has " << _number_elements << " atoms\n";
  debug_print(cerr);
#endif

  int ndx = 0;
  for (int i = 0; i < _number_elements; ++i)
  {
    xref[i] = ndx;
    ndx++;

    const Atom * a = _things[i];

    if (1 != a->atomic_number())
      continue;

    if (0 != a->isotope() || 0 != a->formal_charge())
      continue;

    const int acon = a->ncon();
    if (acon > 1)
      continue;

    is_hydrogen[i] = 1;
    xref[i] = -1;     // in the new molecule, what is the new atom number for atom I
    ndx--;

    if (0 == acon)
      continue;

    const atom_number_t j = a->other(i, 0);

    hcount[j]++;

    _things[j]->remove_bonds_to_atom(i); 
    _things[i]->remove_bonds_to_atom(j);

#ifdef DEBUG_REMOVE_EXPLICIT_HYDROGENS
    cerr << "Atoms " << i << " and " << j << " are across a Hydrogen bond\n";
#endif

    Chiral_Centre * c = chiral_centre_at_atom(j);

    if (NULL == c)
      continue;

    c->atom_is_now_implicit_hydrogen(i);

    _atom_being_unbonded_check_directional_bonds(i);
  }

#ifdef DEBUG_REMOVE_EXPLICIT_HYDROGENS
  for (int i = 0; i < _number_elements; ++i)
  {
    cerr << " i = " << i << ' ' << _things[i]->atomic_symbol() << " xref " << xref[i] << endl;
  }

  cerr << "NDX " << ndx << " cmp " << _number_elements << endl;
  debug_print(cerr);
#endif

  if (_number_elements == ndx)
    return 0;

  if (2 == _number_elements && 2 == ndx)    // hydrogen molecule??
    return 0;

#ifdef DEBUG_REMOVE_EXPLICIT_HYDROGENS
  cerr << "Molecule::remove_explicit_hydrogens:will remove explicit Hydrogen atoms\n";
#endif

// First deal with atoms to which our Hydrogens were attached

  for (int i = 0; i < _number_elements; ++i)
  {
    if (0 == hcount[i])   // atom had no explicit hyddrogens
      continue;

    if (! _things[i]->element()->organic())    // non organic, must set the known flag
    {
      _things[i]->set_implicit_hydrogens(hcount[i], 1);
      _things[i]->set_implicit_hydrogens_known(1);
    }
    else if (_things[i]->implicit_hydrogens_known())   // wow, organic, implicit Hydrogens known, but had an explicit H. 
    {
      _things[i]->set_implicit_hydrogens(hcount[i], 1);
      _things[i]->set_implicit_hydrogens_known(1);
    }
    else                               // let it go free
      _things[i]->set_implicit_hydrogens_known(0);
  }

  _bond_list.remove_bonds_involving_these_atoms(is_hydrogen);

  _bond_list.new_atom_numbers(xref);

  const int nc = _chiral_centres.number_elements();

  for (int i = 0; i < nc; ++i)
  {
    _chiral_centres[i]->new_atom_numbers(xref);
  }

  for (int i = _number_elements - 1; i >= 0; --i)
  {
    if (! is_hydrogen[i])
      continue;

    remove_item(i);
    remove_atom_from_charge_arrays(i);
  }

  _set_modified();

   return ndx;
}

int
Molecule::remove_all_non_natural_elements()
{
  assert(ok());

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (NOT_AN_ELEMENT == _things[i]->atomic_number())
    {
      rc++;
      (void) remove_atom(i);
      i--;
    }
  }

  return rc;
}

const IWString &
Molecule::atomic_symbol (atom_number_t a) const
{
  assert(ok_atom_number(a));

  const Element * e = elementi(a);

  return e->symbol();
}

int
Molecule::organic_only() const
{
  assert(ok());

  for (int i = 0; i < _number_elements; i++)
  {
    const Element * e = _things[i]->element();
    if (! e->organic())
      return 0;
  }

  return 1;
}

int
Molecule::swap_atoms (int i1, int i2,
                      int call_set_modified) 
{
  assert(ok_2_atoms(i1, i2));

  Atom * a1 = _things[i1];
  Atom * a2 = _things[i2];

  _things[i1] = a2;
  _things[i2] = a1;

  if (NULL != _charges)
    _charges->swap_elements(i1, i2);

  if (NULL != _atom_type)
    _atom_type->swap_elements(i1, i2);

  for (int i = 0; i < _chiral_centres.number_elements(); i++)
  {
    Chiral_Centre * c = _chiral_centres[i];

    if (c->involves(i1) || c->involves(i2))
      c->atom_numbers_are_swapped(i1, i2);
  }

// Should do cis-trans bonds too

  for (auto i = 0; i < _things[i1]->number_elements(); ++i)
  {
    _things[i1]->item(i)->swap_atoms(i1, i2);
  }

  for (auto i = 0; i < _things[i2]->number_elements(); ++i)
  {
    _things[i2]->item(i)->swap_atoms(i1, i2);
  }

  if (call_set_modified)
    _set_modified();

  return 1;
}

int
Molecule::move_atom_to_end_of_atom_list (atom_number_t zatom)
{
  assert(ok_atom_number(zatom));

  if (_number_elements - 1 == zatom)    // already the last atom
    return 0;

  Atom * a = _things[zatom];

  for (int i = zatom; i < _number_elements - 1; i++)
  {
    _things[i] = _things[i + 1];
  }

  _things[_number_elements - 1] = a;   

  for (int i = 0; i < _chiral_centres.number_elements(); i++)
  {
    Chiral_Centre * c = _chiral_centres[i];

    c->move_atom_to_end_of_atom_list(zatom, _number_elements);
  }

// Should do cis-trans bonds too

  _set_modified();

  return _bond_list.move_atom_to_end_of_atom_list(zatom, _number_elements);
}

int
Molecule::is_halogen (atom_number_t a) const
{
  assert(ok_atom_number(a));

  const Element * e = _things[a]->element();

  return e->is_halogen();
}

int
Molecule::delete_fragment (int frag)
{
  assert(ok());

//int nf = number_fragments();

  assert(frag >= 0 && frag < number_fragments());

  Set_of_Atoms atoms_to_be_removed;

  _fragment_information.atoms_in_fragment(_number_elements, frag, atoms_to_be_removed);

  assert(atoms_to_be_removed.number_elements() > 0);

  return remove_atoms(atoms_to_be_removed);
}

int
Molecule::delete_fragments (const resizable_array<int> & to_be_deleted)
{
  (void) number_fragments();

  const int * fragment_membership = _fragment_information.fragment_membership();

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    if (to_be_deleted.contains(fragment_membership[i]))
      atoms_to_be_removed.add(i);
  }

  return remove_atoms(atoms_to_be_removed);
}

int
Molecule::delete_fragments (const int * fragments_to_be_deleted)
{
  (void) number_fragments();

  int * atoms_to_be_deleted = new_int(_number_elements); std::unique_ptr<int[]> free_atoms_to_be_deleted(atoms_to_be_deleted);

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    int f = fragment_membership[i];

    if (fragments_to_be_deleted[f])
      atoms_to_be_deleted[i] = 1;
  }

  return remove_atoms(atoms_to_be_deleted);
}

int
Molecule::delete_all_fragments_except (int frag)
{
  assert(ok());

  assert(frag >= 0 && frag < number_fragments());

  const int * fragment_membership = _fragment_information.fragment_membership();

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    if (frag != fragment_membership[i])
      atoms_to_be_removed.add(i);
  }

  assert(atoms_to_be_removed.number_elements());

  return remove_atoms(atoms_to_be_removed);
}

distance_t
Molecule::distance_between_atoms (atom_number_t a1,
                                  atom_number_t a2) const
{
  assert(ok_2_atoms(a1, a2));

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];

  return aa1->distance(*aa2);
}

distance_t
Molecule::longest_intra_molecular_distance() const
{
  distance_t rc = 0.0;

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * ai = _things[i];

    for (int j = 0; j < _number_elements; j++)
    {
      const Atom * aj = _things[j];

      distance_t d = ai->distance(*aj);

      if (d > rc)
        rc = d;
    }
  }

  return rc;
}

void
Molecule::compute_centre (const Set_of_Atoms * s, Coordinates & result) const
{
  coord_t x = 0.0;
  coord_t y = 0.0;
  coord_t z = 0.0;

  int ns = s->number_elements();
  for (int i = 0; i < ns; i++)
  {
    const Atom * a = _things[s->item(i)];
    x += a->x();
    y += a->y();
    z += a->z();
  }

  result.setxyz(x / float(ns), y / float(ns), z / float(ns));

  return;
}

int
Molecule::add_extra_text_info (IWString * extra)
{
  return _text_info.add(extra);
}

int
Molecule::add_extra_text_info (const IWString & extra)
{
  IWString * tmp = new IWString(extra);
  return _text_info.add(tmp);
}

int
Molecule::add_extra_text_info (const char * extra)
{
  IWString * tmp = new IWString(extra);

  return _text_info.add(tmp);
}

int
Molecule::copy_extra_text_info_to (Molecule & rhs) const
{
  int ninfo = _text_info.number_elements();
  for (int i = 0; i < ninfo; i++)
  {
    const IWString & infi = *(_text_info[i]);

    IWString * tmp = new IWString(infi);

    rhs._text_info.add(tmp);
  }

  return ninfo;
}

void
Molecule::discard_extra_text_info()
{
  _text_info.resize(0);

  return;
}

int
Molecule::centroid (Coordinates & result, int frag)
{
  result.setxyz(0.0, 0.0, 0.0);
  
  assert(frag >= 0 && frag < number_fragments());

  const int * fragment_membership = _fragment_information.fragment_membership();

  int atoms_included = 0;

  for (int i = 0; i < _number_elements; i++)   // i is atom number
  {
    if (frag == fragment_membership[i])
    {
      result += *(_things[i]);
      atoms_included++;
    }
  }

  result /= static_cast<coord_t>(atoms_included);

  return 1;
}

/*
  compute centroid for each fragment.
*/

int
Molecule::centroids (resizable_array_p<Coordinates> & result)
{
  assert(0 == result.number_elements());

  int nf = number_fragments();
  result.resize(nf);

// If only one fragment call the method which does not compute fragment membership

  if (1 == nf)
  {
    Coordinates * c = new Coordinates;
    centroid(*c);

    result.add(c);

    return 1;
  }

  for (int i = 0; i < nf; i++)
  {
    Coordinates * c = new Coordinates;

    centroid(*c, i);

    result.add(c);
  }

  return 1;
}

//#define DEBUG_STEREO_PRESERVING_SUBSTITUTE

/*
  We are substituting atom A2 for atom A1. Atom A1 is bonded to atom C.
  We care about preserving any stereochemistry associated with atom C
*/

int
Molecule::stereo_preserving_substitute (atom_number_t c,
                                        const atom_number_t a1,
                                        const atom_number_t a2)
{
#ifdef DEBUG_STEREO_PRESERVING_SUBSTITUTE
  cerr << "Molecule::stereo_preserving_substitute: c = " << c << " a1 = " << a1 << " a2 = " << a2 << endl;
#endif

  assert(ok_3_atoms(c, a1, a2));

  Atom * ac = _things[c];

  Bond * bca1 = NULL;

  for (auto b : *ac)
  {
    const auto x = b->other(c);

    if (a1 == x)
      bca1 = const_cast<Bond *>(b);
    else if (a2 == x)     // almost certainly too hard
    {
      cerr << "Molecule::stereo_preserving_substitute:atom " << a2 << " alread bonded to " << c << " replace " << a1 << endl;
      return 0;
    }
  }

  if (NULL == bca1)
  {
    cerr << "Molecule::stereo_preserving_substitute:atoms " << c << " and " << a1 << " not bonded\n";
    return 0;
  }

// Change the atoms in the bond

  if (a1 == bca1->a1())    
    bca1->set_a1(a2);
  else
    bca1->set_a2(a2);

#ifdef DEBUG_STEREO_PRESERVING_SUBSTITUTE
  cerr << "Final bond " << bca1->a1() << " to " << bca1->a2() << endl;
#endif

// Tell A1 that he is no longer bonded to C

  _things[a1]->remove_bonds_to_atom(c);

// Now tell A2 that he is now bonded to C

  _things[a2]->add(bca1);

#ifdef DEBUG_STEREO_PRESERVING_SUBSTITUTE
  assert(are_bonded(c, a2));
  assert(! are_bonded(c, a1));
#endif

  Chiral_Centre * cc = chiral_centre_at_atom(c);

  if (NULL != cc)
  {
    if (! cc->change_atom_number(a1, a2))
    {
      cerr << "Molecule::stereo_preserving_substitute: cannot change atom numbers for chiral center on atom " << c << endl;
      cerr << "new atoms " << a1 << " and " << a2 << endl;

      return 0;
    }
  }

  _set_modified();

  return 1;
}

/*
  Replace A1 with A2
*/

int
Molecule::stereo_preserving_substitute (atom_number_t a1,
                                        atom_number_t a2)
{
  assert(ok_2_atoms(a1, a2));

#ifdef DEBUG_STEREO_PRESERVING_SUBSTITUTE
  cerr << "Molecule::stereo_preserving_substitute: begin a1 = " << a1 << " a2 = " << a2 << " natoms = " << _number_elements << " nb = " << _bond_list.number_elements() << endl;
  cerr << "Molecule::stereo_preserving_substitute:initial nrings " << nrings() << endl;
  debug_print(cerr);
#endif

  Atom * at1 = _things[a1];
  Atom * at2 = _things[a2];

  if (at1->is_bonded_to(a2))
  {
    cerr << "Molecule::stereo_preserving_substitute:atoms " << a1 << " and " << a2 << " bonded, cannot process\n";
    return 0;
  }

  for (int i = _chiral_centres.number_elements() - 1; i >= 0; i--)
  {
    Chiral_Centre * c = _chiral_centres[i];

    int inva1 = c->involves(a1);
    int inva2 = c->involves(a2);

    if (! inva1 && ! inva2)
      continue;

    if (inva1 && inva2)    // too wierd, how could this happen, maybe some kind of chiral-rearrangent...
      continue;

    int is_central_atom1 = (c->a() == a1);   // special case if A1 is the centre of a chiral centre
    int is_central_atom2 = (c->a() == a2);   // special case if A2 is the centre of a chiral centre

    int rc;
    if (inva1)
      rc = c->change_atom_number(a1, a2);
    else
      rc = c->change_atom_number(a2, a2);

    if (0 == rc)
    {
      cerr << "Molecule::stereo_preserving_substitute: cannot change atom numbers for chiral center on atom " << c->a() << endl;
      cerr << "new atoms " << a1 << " and " << a2 << endl;
      debug_print(cerr);

      return 0;
    }

    if (! is_central_atom1 && ! is_central_atom2)
      continue;

//  Chiral centre has A1 or A2 as its centre. Unless there is an empty slot on the chiral centre object,
//  we must delete it

    if (0 == c->implicit_hydrogen_count() && 0 == c->lone_pair_count())
    {
      _chiral_centres.remove_item(i);
      continue;
    }

//  Unless the replacement atom has just one connection, we must delete the chiral centre

    atom_number_t o;
    if (is_central_atom1 && 1 == at2->ncon())
      o = at2->other(a2, 0);
    else if (is_central_atom2 && 1 == at1->ncon())
      o = at1->other(a1, 0);
    else
    {
      _chiral_centres.remove_item(i);
      continue;
    }

//  We replace the lone pair or Hydrogen with the one atom bonded to the replacement atom

    if (c->implicit_hydrogen_count())
      c->implicit_hydrogen_is_now_atom_number(o);
    else
      c->lone_pair_is_now_atom_number(o);
  }

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    atom_number_t o;    // the other atom involved in the bond - not A1

    if (a1 == b->a1())
    {
      b->set_a1(a2);
      o = b->a2();
    }
    else if (a1 == b->a2())
    {
      b->set_a2(a2);
      o = b->a1();
    }
    else
      continue;

    if (! at2->is_bonded_to(o))
      at2->add(b);
  }

  at2->set_modified();

  at1->resize(0);   // get rid of all the bonds
  

  at1->set_implicit_hydrogens_known(0);

  at1->set_modified();

  _set_modified();

#ifdef DEBUG_STEREO_PRESERVING_SUBSTITUTE
  cerr << "After stereo_preserving_substitute between " << a1 << " and " << a2 << " molecule is\n";
  debug_print(cerr);
#endif

  return 1;
}

int
Molecule::highest_coordinate_dimensionality() const
{
  int rc = 1;    

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * ai = _things[i];

    if (ai->z() != static_cast<coord_t>(0.0))
      return 3;

    if (ai->x() != static_cast<coord_t>(0.0) || ai->y() != static_cast<coord_t>(0.0))
      rc = 2;
  }

  return rc;
}

int
Molecule::convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s,
                                                int * barray)
{
  assign_bond_numbers_to_bonds_if_needed();

  return _convert_set_of_atoms_to_bond_numbers(s, barray);
}

int
Molecule::convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s,
                                                int * barray) const
{
  if (0 == _bond_list.number_elements())
    return 0;

  assert(_bond_list[0]->bond_number_assigned());

  return _convert_set_of_atoms_to_bond_numbers(s, barray);
}

int
Molecule::_convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s,
                                                 int * barray) const
{
  int rc = 0;

  int n = s.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = s[i];

    const Atom * aj = _things[j];

    for (int k = i + 1; k < n; k++)
    {
      atom_number_t j2 = s[k];

      const Bond * b = aj->bond_to_atom(j,j2);

      if (NULL == b)
        continue;

      int bn = b->bond_number();

      barray[bn] = 1;

      rc++;
    }
  }

  return rc;
}

int
Molecule::contains_non_periodic_table_elements() const
{
  assert(ok());

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->atomic_number() <= 0)
      return 1;
  }

  return 0;
}

void *
Molecule::user_specified_atom_void_ptr (atom_number_t zatom) const
{
  assert(ok_atom_number(zatom));

  return _things[zatom]->user_specified_void_ptr();
}

void
Molecule::set_user_specified_atom_void_ptr (atom_number_t zatom, void * v)
{
  assert(ok_atom_number(zatom));

  _things[zatom]->set_user_specified_void_ptr(v);

  return;
}

void
Molecule::clear_all_user_specified_atom_pointers()
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_user_specified_void_ptr(NULL);
  }

  return;
}

const Atom *
Molecule::atom_with_user_specified_void_ptr (const void * v) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (v == _things[i]->user_specified_void_ptr())
      return _things[i];
  }

  return NULL;
}

void
Molecule::spatial_extremeties (coord_t & xmin, coord_t & xmax, coord_t & ymin, coord_t & ymax) const
{
  if (0 == _number_elements)
    return;

  const Atom * a = _things[0];

  xmin = a->x();
  xmax = a->x();
  ymin = a->y();
  ymax = a->y();

  coord_t c;

  for (int i = 1; i < _number_elements; i++)
  {
    a = _things[i];

    c = a->x();

    if (c < xmin)
      xmin = c;
    else if (c > xmax)
      xmax = c;

    c = a->y();
    if (c < ymin)
      ymin = c;
    else if (c > ymax)
      ymax = c;
  }

  return;
}
void
Molecule::spatial_extremeties (coord_t & xmin, coord_t & xmax,
                               coord_t & ymin, coord_t & ymax,
                               coord_t & zmin, coord_t & zmax) const
{
  if (0 == _number_elements)
    return;

  const Atom * a = _things[0];

  xmin = a->x();
  xmax = a->x();
  ymin = a->y();
  ymax = a->y();
  zmin = a->z();
  zmax = a->z();

  coord_t c;

  for (int i = 1; i < _number_elements; i++)
  {
    a = _things[i];

    c = a->x();

    if (c < xmin)
      xmin = c;
    else if (c > xmax)
      xmax = c;

    c = a->y();
    if (c < ymin)
      ymin = c;
    else if (c > ymax)
      ymax = c;

    c = a->z();
    if (c < zmin)
      zmin = c;
    else if (c > zmax)
      zmax = c;
  }

  return;
}

void
Molecule::spatial_extremeties_x (coord_t & xmin, coord_t & xmax) const
{
  if (0 == _number_elements)
    return;

  xmin = _things[0]->x();
  xmax = xmin;

  for (auto i = 0; i < _number_elements; ++i)
  {
    const auto x = _things[i]->x();

    if (x < xmin)
      xmin = x;
    else if (x > xmax)
      xmax = x;
  }

  return;
}

void
Molecule::spatial_extremeties_x (atom_number_t & left,
                                 atom_number_t & right) const
{
  if (0 == _number_elements)
    return;

  coord_t xmin = _things[0]->x();
  coord_t xmax = xmin;
  left = 0;
  right = 0;

  for (auto i = 0; i < _number_elements; ++i)
  {
    const auto x = _things[i]->x();

    if (x < xmin)
    {
      xmin = x;
      left = i;
    }
    else if (x > xmax)
    {
      xmax = x;
      right = i;
    }
  }

  return;
}

int
Molecule::count_heteroatoms (const Set_of_Atoms & r) const
{
  int n = r.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    if (6 != _things[j]->atomic_number())
      rc++;
  }

  return rc;
}

template <typename T>
int
write_isotopically_labelled_smiles (Molecule & m, const bool uniq, T & output)
{
  const int matoms = m.natoms();
  int * isosave = new int[matoms]; std::unique_ptr<int[]> free_isosave(isosave);

  m.get_isotopes(isosave);

  for (int i = 0; i < matoms; ++i)
  {
    m.set_isotope(i, i);
  }

  if (uniq)
    output << m.unique_smiles();
  else
    output << m.smiles();

  output << ' ' << m.name();

  m.set_isotopes(isosave);

  return 1;
}

template int write_isotopically_labelled_smiles(Molecule & m, const bool, std::ostream &);
template int write_isotopically_labelled_smiles(Molecule & m, const bool, IWString_and_File_Descriptor &);

void
Molecule::reset_all_atom_map_numbers()
{
  int changes = 0;
  for (int i = 0; i < _number_elements; ++i)
  {
    if (0 == _things[i]->atom_map())
      continue;

    _things[i]->set_atom_map(0);
    changes = 1;
  }

  if (changes)
    unset_unnecessary_implicit_hydrogens_known_values();

  return;
}

void
Molecule::set_atom_map_number(const atom_number_t zatom, const int s)
{
  assert(ok_atom_number(zatom));

  _things[zatom]->set_atom_map(s);

  invalidate_smiles();

  return;
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

/*
  Pack some info into a 64 bit word

  We keep track of various things that can be discerned from a single pass through the atom array
  Mostly we allocate 4 bits to each.
  Note that we are not too careful about this. We do not consider charge for example

  Bits Max
  1    1
  2    2
  3    4
  4    8
  5   16
  6   32
  7   64
*/

uint64_t
Molecule::quick_atom_hash() const
{
  if (0 == _number_elements)
    return 0;

  uint64_t cd1 = 0;   // 4 bits
  uint64_t cd2 = 0;   // 6 bits
  uint64_t cd3 = 0;   // 5 bits
  uint64_t cd4 = 0;   // 4 bits

  uint64_t nd1 = 0;   // 4 bits
  uint64_t nd2 = 0;   // 4 bits
  uint64_t nd3 = 0;   // 4 bits
  uint64_t nd4 = 0;   // 1 bits

  uint64_t od1 = 0;   // 4 bits
  uint64_t od2 = 0;   // 4 bits

  uint64_t f = 0;     // 3 bits
  uint64_t sd1 = 0;   // 2 bits
  uint64_t sd2 = 0;   // 3 bits
  uint64_t sdx = 0;   // 3 bits
  uint64_t cl = 0;    // 3 bits
  uint64_t x = 0;    // 3 bits   

  uint64_t h = 0;   // 6 bits

  uint64_t fc = 0;  // 1 bit


  for (int i = 0; i < _number_elements; ++i)
  {
    Atom * a = const_cast<Atom *>(_things[i]);        // beware, loss of const not really a good thing to do...
    const atomic_number_t z = a->atomic_number();
    const int acon = a->ncon();

    if (6 == z)
    {
      if (1 == acon)
      {
        cd1++;
        h += 3;
      }
      else if (2 == acon)
      {
        cd2++;
        h += 2;
      }
      else if (3 == acon)
      {
        cd3++;
        h += 1;
      }
      else
        cd4++;
    }
    else if (7 == z)
    {
      if (1 == acon)
      {
        nd1++;
        h += 2;
      }
      else if (2 == acon)
      {
        nd2++;
        h++;
      }
      else if (3 == acon)
        nd3++;
      else   
      {
        nd4++;
        fc++;
      }

    }
    else if (8 == z)
    {
      if (1 == acon)
      {
        od1++;
        h++;
      }
      else
        od2++;
    }
    else if (9 == z)
      f++;
    else if (16 == z)
    {
      if (1 == acon)
      {
        sd1++;
        h++;
      }
      else if (2 == acon)
        sd2++;
      else
        sdx++;
    }
    else if (17 == z)
      cl++;
    else if (1 == z)
      h++;
    else
    {
      x++;
      h += a->implicit_hydrogens();
      if (a->formal_charge())
        fc++;
    }
  }

// Now put all that back into the final result

  uint64_t rc = 0;

  int bshift = 64;

  bshift -= 4; rc = do_or(cd1, bshift, 8, rc);
  bshift -= 6; rc = do_or(cd2, bshift, 32, rc);
  bshift -= 5; rc = do_or(cd3, bshift, 16, rc);
  bshift -= 4; rc = do_or(cd4, bshift, 8, rc);

  bshift -= 4; rc = do_or(nd1, bshift, 8, rc);
  bshift -= 4; rc = do_or(nd2, bshift, 8, rc);
  bshift -= 4; rc = do_or(nd3, bshift, 8, rc);
  bshift -= 1; rc = do_or(nd4, bshift, 2, rc);

  bshift -= 4; rc = do_or(od1, bshift, 8, rc);
  bshift -= 4; rc = do_or(od2, bshift, 8, rc);

  bshift -= 3; rc = do_or(f, bshift, 4, rc);

  bshift -= 2; rc = do_or(sd1, bshift, 2, rc);
  bshift -= 3; rc = do_or(sd2, bshift, 4, rc);
  bshift -= 3; rc = do_or(sdx, bshift, 4, rc);

  bshift -= 3; rc = do_or(cl, bshift, 4, rc);
  bshift -= 3; rc = do_or(x, bshift, 4, rc);
  bshift -= 6; rc = do_or(h, bshift, 32, rc);

  bshift -= 1; rc = do_or(fc, bshift, 1, rc);

//cerr << " from " << cd1 << ' ' << cd2 << ' ' << cd3 << ' ' << cd4 << ' ' << nd1 << ' ' << nd2 << ' ' << nd3 << ' ' << nd4 << ' ' << od1 << ' ' << od2 << ' ' << f << ' ' << sd1 << ' ' << sd2 << ' ' << sdx << ' ' << cl << ' ' << x << ' ' << fc << " get " << rc << endl;

  return rc;
}

atom_number_t
Molecule::atom_with_atom_map_number(const int n) const
{
  for (int i = 0; i < _number_elements; ++i)
  {
    if (n == _things[i]->atom_map())
      return i;
  }

  return INVALID_ATOM_NUMBER;
}
