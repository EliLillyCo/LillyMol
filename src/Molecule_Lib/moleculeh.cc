#include <iostream>
#include <iomanip>
using std::cerr;
using std::endl;

#define COMPILING_MOLECULEH_H

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "molecule.h"

#include "misc.h"
#include "misc2.h"

//#define DEBUG_COMPUTE_IMPLICIT_HYDROGENS

static int _display_messages_about_unable_to_compute_implicit_hydgogens = 1;

void 
set_display_messages_about_unable_to_compute_implicit_hydgogens (int s)
{
  _display_messages_about_unable_to_compute_implicit_hydgogens = s;

  return;
}

int
Molecule::_compute_and_store_implicit_hydrogens (atom_number_t zatom)
{
  Atom * a = _things[zatom];

  if (a->implicit_hydrogens_known())
    return 1;

  int ih;
  if (! a->compute_implicit_hydrogens(ih))
  {
    if (_display_messages_about_unable_to_compute_implicit_hydgogens)
    {
      cerr << "Molecule::_compute_and_store_implicit_hydrogens: cannot determine atom " << zatom << 
              " (" << a->atomic_symbol() << ")\n";
      cerr << _molecule_name << endl;
    }
    a->set_implicit_hydrogens(0);
    return 0;
  }

//cerr << "Atom " << zatom << " '" << smarts_equivalent_for_atom (zatom) << "' has " << ih << " implicit hydrogens\n";

  a->set_implicit_hydrogens(ih);
  return 1;
}

int
Molecule::_compute_implicit_hydrogens (atom_number_t i, int & result) const
{
  return _things[i]->compute_implicit_hydrogens(result);
}

int
Molecule::_compute_implicit_hydrogens (atom_number_t a)
{
  int ih;
  if (! _compute_implicit_hydrogens(a, ih))
    return 0;

  return ih;
}

/*
  Looking at this, I'm not respecting the setting of the sticky bit.
  Just why did I have a sticky bit in the first place?
*/

/*int
Molecule::compute_implicit_hydrogens()
{
  int rc = 1;
  for (int i = 0; i < _number_elements; i++)
    (void) _things[i]->implicit_hydrogens();

#ifdef DEBUG_COMPUTE_IMPLICIT_HYDROGENS
  cerr << "After computing implicit hydrogens for " << _number_elements << " atoms\n";
  for (int i = 0; i < _number_elements; i++)
    cerr << "Atom " << i << " has " << _things[i]->implicit_hydrogens() << " implicit hydrogens\n";
#endif

  return rc;
}*/

int
Molecule::recompute_implicit_hydrogens (atom_number_t a)
{
  return _compute_and_store_implicit_hydrogens(a);
}

int
Molecule::recompute_implicit_hydrogens()
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (! _compute_and_store_implicit_hydrogens(i))
      rc = 0;
  }

  return rc;
}

int
Molecule::implicit_hydrogens()
{
  int rc = 0;
  
  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i]->implicit_hydrogens();
  }

  return rc;
}

/*
  The implicit hydrogen count is somewhat complex.
  Sometimes a known value is available, typically from the input
  file. Most other times, this value is derived.

  Derived values are stored in the array as regular numbers.
  Values which are considered immutable should be entered with
  OVERRIDE = 1. These will be stored with the hcount_sticky_bit
  set can only be changed with another call with OVERRIDE = 1
*/

int
Molecule::set_implicit_hydrogens (atom_number_t i, int ih,
                                  int override)
{
  assert(ih < 6);     // hard to imagine > 5 implicit H's

  int previous_value = implicit_hydrogens(i);   // calls OK_ATOM_NUMBER

  if (previous_value == ih)    // unchanged
    return 1;

  _smiles_information.invalidate();

  _symmetry_class_and_canonical_rank.invalidate();

  return _things[i]->set_implicit_hydrogens(ih, override);
}

int
Molecule::set_implicit_hydrogens_known (atom_number_t a, int known)
{
  assert(ok_atom_number(a));

  _things[a]->set_implicit_hydrogens_known(known);

  _smiles_information.invalidate();

  return 1;
}

//#define DEBUG_IMPLICIT_HYDROGENS

/*
  Returning the number of implicit hydrogens is complicated by
  the need to remove the implicit_hcount_sticky_bit if present.
*/

int
Molecule::implicit_hydrogens (atom_number_t zatom)
{
  return _things[zatom]->implicit_hydrogens();
}

int
Molecule::implicit_hydrogens (atom_number_t a, int & sticky_bit_set)
{
  assert(ok_atom_number(a));

  int rc = implicit_hydrogens(a);

  if (_things[a]->implicit_hydrogens_known())
    sticky_bit_set = 1;
  else
    sticky_bit_set = 0;

  return rc;
}

int
Molecule::unset_all_implicit_hydrogen_information (atom_number_t zatom)
{
  _things[zatom]->unset_all_implicit_hydrogen_information();

  return 1;
}

/*
  this function returns the number of hydrogen atom connections to a given
  atom number.
*/

int
Molecule::explicit_hydrogens (atom_number_t i) const
{
  assert(ok_atom_number(i));

  int hcount = 0;

  const Atom * a = _things[i];
  int acon = a->number_elements();
  for (int j = 0; j < acon; j++)
  {
    const Bond * b = a->item(j);

    atom_number_t k = b->other(i);

    if (1 == _things[k]->atomic_number())
      hcount++;
  }
  
  return hcount;
}

#include "cmdline.h"

/*
  Lots of choices for the coordinates to be assigned to newly created H
  atoms.
*/

#define NEW_H_GET_ZERO_COORDINATES 1
#define NEW_H_GET_ZOFFSET_COORDINATES 2
#define NEW_H_GET_CAREFUL_COORDINATES 3

static int new_h_coordinates = NEW_H_GET_CAREFUL_COORDINATES;

int
process_explicit_h_options (Command_Line & cl, const int verbose, const char cflag)
{
  int i = 0;
  IWString tmp;
  while (cl.value(cflag, tmp, i++))
  {
    if ('0' == tmp || "0.0" == tmp)
    {
      new_h_coordinates = NEW_H_GET_ZERO_COORDINATES;
      if (verbose)
        cerr << "Newly created Hydrogens get zero coordinates\n";
    }
    else if ('Z' == tmp)
    {
      new_h_coordinates = NEW_H_GET_ZOFFSET_COORDINATES;
      if (verbose)
        cerr << "Newly created Hydrogens get Z axis offset coordinates only\n";
    }
    else if ("careful" == tmp)
    {
      new_h_coordinates = NEW_H_GET_CAREFUL_COORDINATES;
      if (verbose)
        cerr << "Newly created Hydrogens get carefully generated coordinates only\n";
    }
    else
    {
      cerr << "Process_explicit_h_options: unrecognised '-" << cflag << "' option '" << tmp << "'\n";
      return 0;
    }
  }

  return 1;
}

static coord_t default_h_bond_length = static_cast<coord_t> (1.08);

void
set_default_h_bond_length (coord_t d)
{
  assert(d > 0.0);

  default_h_bond_length = d;

  return;
}

#include "chiral_centre.h"

/*
  Placing an extra H on a chiral centre is difficult.
  H is the newly created Hydrogen atom, and AH is its atom number
*/

/*int
Molecule::_place_chiral_h_atom (Chiral_Centre * c,
                                Atom * h,
                                atom_number_t ah)
{
  assert(c->involves (ah));

  atom_number_t anchor = c->a();

  c->make_top_front(ah);
  assert(ah == c->top_front());

  Coordinates va;
  get_coords(anchor, va);

  Coordinates v_top_back, v_left_down, v_right_down;
  get_coords(c->top_back(), v_top_back);
  get_coords(c->left_down(), v_left_down);
  get_coords(c->right_down(), v_right_down);

  v_top_back -= va;
  v_left_down -= va;
  v_right_down -= va;

  v_top_back.normalise();
  v_left_down.normalise();
  v_right_down.normalise();

#ifdef DEBUG_PLACE_CHIRAL_H_ATOM
  cerr << "Top back vector is " << v_top_back << endl;
  cerr << "Left down vector is " << v_left_down << endl;
  cerr << "right down vector is " << v_right_down << endl;
#endif

  Coordinates bl, br;
  bl = v_top_back;
  br = v_top_back;

  bl.cross_product(v_left_down);
  br.cross_product(v_right_down);

#ifdef DEBUG_PLACE_CHIRAL_H_ATOM
  cerr << "bl is " << bl << endl;
  cerr << "br is " << br << endl;
#endif

// This is easy if we have the straightfoward case of
//
//         |
//        / \

// In which case the Z component of one cross product vector is
// negative and the other is positive. What happens when both
// signs are the same? Basically we do something random, it is
// too hard to get right.

  coord_t z;
  if (bl.z() < 0.0 && br.z() > 0.0)
    z = default_h_bond_length;
  else if (bl.z() > 0.0 && br.z() < 0.0)
    z = - default_h_bond_length;
  else if (bl.z() > 0.0 && br.z() > 0.0)
  {
    z = default_h_bond_length;
  }
  else    // both are negative
  {
    z = -default_h_bond_length;
  }

  const Atom * a = _things[anchor];

  h->setxyz (a->x(), a->y(), z);

  return 1;
}*/

/*
  We need to set the coordinates for newly created Hydrogen atom H
  It is bonded to atom ANCHOR
*/

int
Molecule::_place_1_hydrogen (const Make_Implicit_Hydrogens_Explicit & mihe)
{
  atom_number_t anchor = mihe.a();

  const Atom * a = _things[anchor];

  int acon = a->ncon();

// Well, I suppose some bozo could have a chiral centre where 1 Hydrogen was
// implicit and one was explicit. Should be shot.

  Chiral_Centre * c = chiral_centre_at_atom(anchor);
  if (NULL != c && 1 != c->implicit_hydrogen_count())
  {
    cerr << "Molecule::_place_1_hydrogen: Hmmm, adding one H to chiral centre at atom " << anchor << endl;
    cerr << "Chiral centre is ";
    c->debug_print(cerr);

    return 0;
  }

  Atom * h = mihe.new_atom();

  add(h);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  if (NULL != c && c->implicit_hydrogen_count())
    c->implicit_hydrogen_is_now_atom_number(_number_elements - 1);

  int dimensionality = mihe.dimensionality();

  if (dimensionality < 2)    // no need to worry about coordinates
    return 1;

  if (2 == dimensionality || 0 == acon || acon > 3)
  {
    coord_t z;
    if ((anchor & 1))       // slightly random
      z = - default_h_bond_length;
    else
      z = default_h_bond_length;

    h->setxyz(a->x(), a->y(), z);

    return 1;
  }

// Needs to be done carefully. Sum all the bond vectors coming from the attachment point
// and then place the extra atom on the other side

  Coordinates va = (*a);

// avoid straight bonds

//cerr << "Just about to put 1 Hydrogen, acon " << acon << " nbonds " << a->nbonds() << endl;

  coord_t blen = default_h_bond_length;
  if (6 != a->atomic_number())
    blen = static_cast<coord_t>(blen*0.9);

  if (1 == acon && 4 != a->nbonds())      // straight bond OK for acetylene
  {
    atom_number_t j = a->other(anchor, 0);
    Coordinates v = *(_things[j]);

    v -= va;
    v.normalise();

    Coordinates vperp( static_cast<coord_t>(a->x() + 1.0), static_cast<coord_t>(a->y() + 1.0), static_cast<coord_t>(a->z() + 1.0) );
    vperp -= va;
    vperp.cross_product(v);
    vperp.normalise();

	vperp = vperp * static_cast<coord_t>(blen/sqrt(5.0) );
	v = v * static_cast<coord_t>( blen * 2 / sqrt(5.0) );

//  cerr << "v is " << v << ", vperp " << vperp << endl;

    h->setxyz(a->x() - v.x() + vperp.x(),
               a->y() - v.y() + vperp.y(), 
               a->z() - v.z() + vperp.z());
//  cerr << "Distance is " << va.distance(*h) << endl;
  }
  else
  {
    Coordinates vsum;
    for (int i = 0; i < acon; i++)
    {
      atom_number_t j = a->other(anchor, i);
      Coordinates v = *(_things[j]);
  
      v -= va;
      v.normalise();
      vsum += v;
    }

    vsum.normalise();
    vsum *= blen;

    h->setxyz(a->x() - vsum.x(), a->y() - vsum.y(), a->z() - vsum.z());
  }

  return 1;
}

int
Molecule::set_coordinates_of_singly_connected_atom (atom_number_t zatom,
                                                    coord_t default_bond_length)
{
  assert(1 == _things[zatom]->ncon());

  atom_number_t anchor = _things[zatom]->other(zatom, 0);
  const Atom * a = _things[anchor];

  int acon = a->ncon();

  if (1 == acon)    // anchor was previously disconnected
  {
    _things[zatom]->setxyz(a->x() + default_bond_length, a->y(), a->z());   // randomly choose the X axis
    return 1;
  }

  if (acon > 4)
  {
    cerr << "Molecule::set_coordinates_of_singly_connected_atom:do not know how to handle " << acon << " connected atom\n";
    return 0;
  }

// Needs to be done carefully. Sum all the bond vectors coming from the attachment point
// and then place the extra atom on the other side

  Coordinates va = (*a);

  Coordinates vsum;
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(anchor, i);
    if (j == zatom)
      continue;

    Coordinates v = *(_things[j]);
  
    v -= va;
    v.normalise();
    vsum += v;
  }

  vsum.normalise();
  vsum *= default_bond_length;

  _things[zatom]->setxyz(a->x() - vsum.x(), a->y() - vsum.y(), a->z() - vsum.z());

// If 2 == acon, we want to avoid straight bonds unless there is a triple bond present

  if (2 == acon && 4 != a->nbonds())
  {
    atom_number_t j = a->other(anchor, 0);
    if (j == zatom)
      j = a->other(anchor, 1);

    set_bond_angle(j, anchor, zatom, static_cast<angle_t>(120.0 * DEG2RAD) );
  }

  return 1;
}

/*
  For now, two extra hydrogens is always one above and one below
*/

int
Molecule::_place_2_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe)
{
  atom_number_t anchor = mihe.a();
  int dimensionality = mihe.dimensionality();

  const Atom * a = _things[anchor];
  int acon = a->ncon();

  Atom * h1 = mihe.new_atom();
  Atom * h2 = mihe.new_atom();

  add(h1);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  add(h2);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);

// take care of water before checking the dimensionality

  if (0 == acon)    // adding hydrogens to a water
  {
    h1->setxyz(a->x() + default_h_bond_length, a->y(), a->z());
    coord_t dx = 0.95 * default_h_bond_length;
    coord_t dy = 0.31 * default_h_bond_length;
    h2->setxyz(a->x() - dx, a->y() + dy, a->z());

    return 1;
  }

  if (dimensionality < 2)
    return 1;

// The 2D case will never come up so just do something silly

  if (2 == dimensionality)
  {
    h1->setxyz(a->x(), a->y(), default_h_bond_length);
    h2->setxyz(a->x(), a->y(), - default_h_bond_length);

    return 1;
  }

  if (acon > 2 || 0 == acon)
  {
    cerr << "Molecule::_place_2_hydrogens: sorry, don't know how to place 2 Hydrogens on something with " << acon << " existing connections\n";
    return 0;
  }

// We need to set two vectors. vsum we just add all the vectors departing from anchor and take the opposite direction
// vperp must be a vector normal to vsum

  Coordinates va = *(a);

  Coordinates v0;
  get_coords(other(anchor, 0), v0);
  v0 -= va;
  v0.normalise();

  Coordinates vsum;
  Coordinates vperp;

  if (2 == acon)
  {
    Coordinates v1;
    get_coords(other(anchor, 1), v1);
    v1 -= va;
    v1.normalise();

    vsum = v0 + v1;

    vsum.normalise();

    vperp = v0;
    vperp *= v1;
  }
  else      // adding 2 hydrogens to -N or =C for example. The orientation is random
  {
    vsum = v0;

    vperp.setxyz( static_cast<coord_t>(va.x() + 1.0), static_cast<coord_t>(va.y() + 1.0), static_cast<coord_t>(va.z() + 1.0) );
    vperp.cross_product(v0);
  }

  vperp.normalise();

  coord_t blen = default_h_bond_length;
  if ( 6 != a->atomic_number() ) 
  {
	  blen *= static_cast<coord_t>(0.9);
  }

  float factor = static_cast<float>( blen / sqrt(2.0) );
  vsum *= factor;
  vperp *= factor;

//cerr << "Angle between vectors " << vsum.angle_between(vperp) << " acon = " << acon << endl;

  h1->setxyz(a->x() - vsum.x() + vperp.x(),
              a->y() - vsum.y() + vperp.y(),
              a->z() - vsum.z() + vperp.z());
  h2->setxyz(a->x() - vsum.x() - vperp.x(),
              a->y() - vsum.y() - vperp.y(),
              a->z() - vsum.z() - vperp.z());

  return 1;
}

int
Molecule::set_coordinates_of_singly_connected_atoms (atom_number_t a1,
                                                     atom_number_t a2,
                                                     coord_t default_h_bond_length)
{
  assert(1 == _things[a1]->ncon());
  assert(1 == _things[a2]->ncon());

  atom_number_t anchor = _things[a1]->other(a1, 0);

  assert(anchor == _things[a2]->other(a2, 0));   // both a1 and a2 must be bonded to the same atom

  const Atom * a = _things[anchor];

  int acon = a->ncon();

  if (acon > 4)
  {
    cerr << "Molecule::set_coordinates_of_singly_connected_atoms: sorry, don't know how to place 2 atoms o n something with " << acon << " connections\n";
    return 0;
  }

// There are really two cases here. If (4 == acon), we create a 90 degree cross. If (3 == ncon)
// then we create a Y shape

  const Coordinates va = *(a);

  if (4 == acon)   // cross shape
  {
    atom_number_t j = a->other(anchor, 0);
    Coordinates vj = *(_things[j]);
    vj -= va;
    vj.normalise();
    vj *= default_h_bond_length;
    _things[a1]->setxyz(a->x() + vj.x(), a->y() + vj.y(), a->z() + vj.z());

    j = a->other(anchor, 1);
    vj = *(_things[j]);
    vj -= va;
    vj.normalise();
    vj *= default_h_bond_length;
    _things[a1]->setxyz(a->x() + vj.x(), a->y() + vj.y(), a->z() + vj.z());

    return 1;
  }

  if (3 == acon)    // Y shape
  {
    atom_number_t j = a->other(anchor, 0);
    if (j == a1 || j == a2)
      j = a->other(anchor, 1);
    if (j == a1 || j == a1)
      j = a->other(anchor, 2);

    Coordinates vj = *(_things[j]);
    vj -= va;
    vj.normalise();
    vj *= default_h_bond_length;

    _things[a1]->setxyz(a->x() + vj.x(), a->y() + vj.y(), a->z() + vj.z());
    _things[a2]->setxyz(a->x() + vj.x(), a->y() + vj.y(), a->z() + vj.z());

    set_bond_angle( j, anchor, a1, static_cast<angle_t>(120 * DEG2RAD) );

    vj.normalise();

    Coordinates vh = *(_things[a1]);
    vh -= va;
    vh.normalise();

    vj *= -1.0;
    vh *= -1.0;

    vj += vh;
    vj.normalise();

    vj *= default_h_bond_length;

    _things[a2]->setxyz(a->x() + vj.x(), a->y() + vj.y(), a->z() + vj.z());

    return 1;
  }

  cerr << "Molecule::set_coordinates_of_singly_connected_atoms:cannot process ncon " << acon << endl;
  return 0;
}


int
Molecule::_place_3_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe)
{
  atom_number_t anchor = mihe.a();
  int dimensionality = mihe.dimensionality();

  const Atom * a = _things[anchor];
  int acon = a->ncon();

  Atom * atoms[3];
  for (int i = 0; i < 3; i++)
  {
    atoms[i] = mihe.new_atom();
    add(atoms[i]);
    add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  }

  if (dimensionality < 2)
    return 1;

  if (2 == dimensionality)
  {
    atoms[0]->setxyz(a->x(), a->y(), - default_h_bond_length);
    atoms[1]->setxyz(static_cast<coord_t>( a->x() - 0.25 * default_h_bond_length ), a->y(), default_h_bond_length);
    atoms[2]->setxyz(static_cast<coord_t>( a->x() + 0.25 * default_h_bond_length ), a->y(), default_h_bond_length);

    return 1;
  }

  coord_t blen = default_h_bond_length;
  if (6 != a->atomic_number()) 
    blen *= static_cast<coord_t>(0.9);

  if (1 == acon)
  {
    Coordinates va = *a;
//  get_coords(anchor, va);

    Coordinates v0 = *(_things[a->other(anchor, 0)]);
    get_coords(other(anchor, 0), v0);

    va -= v0;
    va.normalise();

//  We need a vector normal to the plane

    Coordinates perp(0.0, 0.0, 1.0);
    perp.cross_product(va);

    atoms[0]->setxyz(static_cast<coord_t>( a->x() + va.x() * blen * 0.25 ),
                      static_cast<coord_t>( a->y() + va.y() * blen * 0.25 ),
                      static_cast<coord_t>( a->z() - blen * 0.5) );

    atoms[1]->setxyz(static_cast<coord_t>( a->x() + 0.5 * va.x() * blen - perp.x() * 0.4 ),
                      static_cast<coord_t>( a->y() + 0.5 * va.y() * blen - perp.y() * 0.4 ),
                      static_cast<coord_t>( a->z() + blen * 0.5) );

    atoms[2]->setxyz(static_cast<coord_t>( a->x() + 0.5 * va.x() * blen + perp.x() * 0.4 ),
                      static_cast<coord_t>( a->y() + 0.5 * va.y() * blen + perp.y() * 0.4 ),
                      static_cast<coord_t>( a->z() +  blen * 0.5) );

    return 1;
  }

  return 1;
}

int
Molecule::_place_4_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe)
{
  atom_number_t anchor = mihe.a();
  int dimensionality = mihe.dimensionality();

  const Atom * a = _things[anchor];

  Atom * atoms[4];
  for (int i = 0; i < 4; i++)
  {
    atoms[i] = mihe.new_atom();
    add (atoms[i]);
    add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  }

  if (dimensionality < 2)
    return 1;

  if (2 == dimensionality)
  {
    atoms[0]->setxyz(static_cast<coord_t>( a->x() - 0.25 * default_h_bond_length ),
                      a->y(),
                      static_cast<coord_t>( 0.60 * default_h_bond_length) );

    atoms[1]->setxyz(static_cast<coord_t>( a->x() - 0.25 * default_h_bond_length ),
                      a->y(),
                      static_cast<coord_t>( -0.60 * default_h_bond_length) );

    atoms[2]->setxyz(static_cast<coord_t>( a->x() + 0.25 * default_h_bond_length ),
                      a->y(),
                      static_cast<coord_t>( 0.60 * default_h_bond_length) );

    atoms[3]->setxyz(static_cast<coord_t>( a->x() + 0.25 * default_h_bond_length ),
                      a->y(),
                      static_cast<coord_t>( -0.60 * default_h_bond_length) );

    return 1;
  }

  cerr << "Molecule::_place_4_hydrogens: don't know how to place 4 Hydrogens\n";

  return 1;
}

int
Molecule::_place_lots_of_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe,
                                    int ih)
{
  atom_number_t anchor = mihe.a();
  int dimensionality = mihe.dimensionality();

  resizable_array<Atom *> atoms;
  atoms.resize(ih);
  for (int i = 0; i < ih; i++)
  {
    Atom * h = mihe.new_atom();
    add(h);
    add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  }

  if (dimensionality < 2)
    return 1;

  cerr << "Molecule::_place_lots_of_hydrogens: don't know how to place " << ih << " hydrogens\n";

  return 1;
}

int
Molecule::make_implicit_hydrogens_explicit()
{
  Make_Implicit_Hydrogens_Explicit mihe;

  return make_implicit_hydrogens_explicit(mihe);
}

int 
Molecule::make_implicit_hydrogens_explicit (atom_number_t a)
{
  Make_Implicit_Hydrogens_Explicit mihe;
  mihe.set_atom(a);

  return make_implicit_hydrogens_explicit(mihe);
}

int
Molecule::valence_ok (atom_number_t a)
{
  assert(ok_atom_number(a));

  return _things[a]->valence_ok();
}

int
Molecule::valence_ok()
{
  for (int i = 0; i < _number_elements; i++)
  {
//  cerr << "Checking atom " << i << " type " << smarts_equivalent_for_atom(i) << endl;

    if (! _things[i]->valence_ok())
      return 0;

// Aug 2005. The valence check does not catch errors
// For now, just forget it...

  }

  return 1;
}

Make_Implicit_Hydrogens_Explicit::Make_Implicit_Hydrogens_Explicit()
{
  _a = INVALID_ATOM_NUMBER;
  _isotope = 0;
  _dimensionality = -1;

  return;
}

Atom *
Make_Implicit_Hydrogens_Explicit::new_atom() const
{
  Atom * rc = new Atom(1);

  if (_isotope)
    rc->set_isotope(_isotope);

  return rc;
}

int
Molecule::make_implicit_hydrogens_explicit (Make_Implicit_Hydrogens_Explicit & mihe)
{
  int na = _number_elements;

  if (mihe.dimensionality() < 0)
    mihe.set_dimensionality(highest_coordinate_dimensionality());

  atom_number_t a = mihe.a();

  if (INVALID_ATOM_NUMBER != a)
  {
    if (a < 0 || a >= _number_elements)
    {
      cerr << "Molecule::make_implicit_hydrogens_explicit: atom number " << a << " is invalid\n";
      return 0;
    }

    int ih = implicit_hydrogens(a);
    if (0 == ih)
      return 1;

    int rc;

    if (1 == ih)
      rc = _place_1_hydrogen(mihe);
    else if (2 == ih)
      rc = _place_2_hydrogens(mihe);
    else if (3 == ih)
      rc = _place_3_hydrogens(mihe);
    else if (4 == ih)
      rc = _place_4_hydrogens(mihe);
    else
      rc = _place_lots_of_hydrogens(mihe, ih);

    _things[a]->set_implicit_hydrogens(0, 1);    // 1 means override known value

    return rc;
  }

  int rc = 1;

  for (int i = 0; i < na; i++)
  {
    mihe.set_atom(i);

    if (! make_implicit_hydrogens_explicit(mihe))    // recursive call
      rc = 0;
  }

  for (int i = na; i < _number_elements; i++)
  {
    const Atom * ai = _things[i];
    if (1 != ai->atomic_number())
      continue;

    atom_number_t j = ai->other(i, 0);
    Atom * aj = _things[j];

    if (ai->distance(*aj) < default_h_bond_length * 1.5)
      continue;

    cerr << "Hydrogen " << i << " bonded to atom " << j << " distance " << ai->distance(*aj) << ". Coordinates changed\n";

    aj->setxyz(static_cast<coord_t>( ai->x() + 0.5 ),
                static_cast<coord_t>( ai->y() + 0.5 ),
                static_cast<coord_t>( ai->z() + 0.5) );
  }

  return rc;
}

// We try to preserve the ordering as much as possible. There are probably a bunch of
// Hydrogens already at the end of the list of atoms

int
Molecule::move_hydrogens_to_end_of_connection_table (atomic_number_t z)
{
  int rc = 0;

  int last_non_hydrogen = _number_elements - 1;

  for (int i = _number_elements - 1; i >= 0; i--)
  {
    if (z != _things[i]->atomic_number())
      break;

    last_non_hydrogen--;
  }

  for (int i = 0; i < last_non_hydrogen; i++)
  {
    const Atom * ai = _things[i];

    if (z != ai->atomic_number())
      continue;

    for (int j = i; j < last_non_hydrogen; j++)
    {
      swap_atoms(j, j + 1, 0);
    }
    i--;
    last_non_hydrogen--;
    rc++;
  }

  if (rc)
    _set_modified();

  return rc;
}

/*
  Someone has a smiles like CC[C]CC where the only thing wrong
  is the implicit hydrogens known attribute
*/

int
Molecule::remove_hydrogens_known_flag_to_fix_valence_errors()
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (! a->implicit_hydrogens_known())
      continue;

    if (a->valence_ok())
      continue;

    int current_value = a->implicit_hydrogens();

    a->unset_all_implicit_hydrogen_information();

    if (a->valence_ok())
      rc++;
    else
    {
      a->set_implicit_hydrogens(current_value);
      a->set_implicit_hydrogens_known(1);
    }
  }

  return rc;
}

int
Molecule::unset_unnecessary_implicit_hydrogens_known_values()
{
  int rc = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    Atom * a = _things[i];

    if (! a->implicit_hydrogens_known())
      continue;

    if (a->isotope() || a->formal_charge())
      continue;

    if (a->element()->organic())
      ;
    else if (5 == a->element()->atomic_number())    // Boron is OK
      ;
    else
      continue;

    int ih;
    if (! a->compute_implicit_hydrogens(ih))
      continue;

    if (ih != a->implicit_hydrogens())
      continue;

    if (! valence_ok(i))    // being very careful
      continue;

    a->set_implicit_hydrogens_known(0);
    rc++;
  }

  return rc;
}
int
Molecule::implicit_hydrogens_known(const atom_number_t zatom) const
{
  return _things[zatom]->implicit_hydrogens_known();
}
