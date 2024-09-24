#include <iostream>
#include <iomanip>
using std::cerr;
using std::endl;

#define COMPILING_MOLECULEH_H

#include "Foundational/cmdline/cmdline.h"

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "molecule.h"

#include "misc2.h"

//#define DEBUG_COMPUTE_IMPLICIT_HYDROGENS

static int _display_messages_about_unable_to_compute_implicit_hydgogens = 1;

namespace moleculeh {

// Disallow things like =[SH]- where the implicit hydrogen is used to
// satisfy an alternate valence.
int implicit_hydrogens_cannot_satisfy_alternate_valence = 1;

void
set_implicit_hydrogens_cannot_satisfy_alternate_valence(int s) {
  implicit_hydrogens_cannot_satisfy_alternate_valence = s;
}

// Return true if the molecule contains things like =S-
// where the implicit Hydrogen on the S makes up the valence.
// This is more likely a structure error.
int
ImplicitHydrogenSatisfiesAlternateValence(Molecule& m)
{
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    int ih = m.implicit_hydrogens(i);
    if (ih == 0) {
      continue;
    }
    const Atom& a = m.atom(i);
    const Element* e = a.element();
    const resizable_array<int>& alternate_valences = e->alternate_valences();
    if (alternate_valences.empty()) {
      continue;
    }

    int existing_connections = a.nbonds() + ih - a.formal_charge();
    if (existing_connections == e->normal_valence()) {
      continue;
    }

    // Special case for Phorphorus.  P(=O)(OC)OC CHEMBL3183964
    if (e->atomic_number() == 15 && ih == 1) {
      continue;
    }
    // CC1=CC=C(C=C1)N1N=C(C=C1NC(=O)C1=C2N=CC=CN2N=C1)C1=CC=C(NS(=O)=O)C=C1 CHEMBL3917931
    if (e->atomic_number() == 16 && a.ncon() == 3 && existing_connections == 6 && ih == 1) {
      continue;
    }

    return 1;
  }

  return 0;  // No instances found
}

}  // namespace moleculeh

void 
set_display_messages_about_unable_to_compute_implicit_hydgogens(int s)
{
  _display_messages_about_unable_to_compute_implicit_hydgogens = s;

  return;
}

int
Molecule::_compute_and_store_implicit_hydrogens(atom_number_t zatom)
{
  Atom * a = _things[zatom];

#ifdef DEBUG_COMPUTE_AND_STORE_IMPLICIT_HYDROGENS
  cerr << "Molecule::_compute_and_store_implicit_hydrogens:implicit_hydrogens_known " << a->implicit_hydrogens_known() << ' ' << smarts_equivalent_for_atom(zatom) << '\n';
#endif
  if (a->implicit_hydrogens_known()) {
    return 1;
  }

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

#ifdef DEBUG_COMPUTE_AND_STORE_IMPLICIT_HYDROGENS
  cerr << "Atom " << zatom << " '" << smarts_equivalent_for_atom (zatom) << "' has " << ih << " implicit hydrogens\n";
#endif

  a->set_implicit_hydrogens(ih);
  return 1;
}

int
Molecule::_compute_implicit_hydrogens(atom_number_t i, int & result) const
{
  return _things[i]->compute_implicit_hydrogens(result);
}

int
Molecule::_compute_implicit_hydrogens(atom_number_t a)
{
  int ih;
  if (! _compute_implicit_hydrogens(a, ih)) {
    return 0;
  }

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
Molecule::recompute_implicit_hydrogens(atom_number_t a)
{
  int rc = _compute_and_store_implicit_hydrogens(a);
  invalidate_smiles();

  return rc;
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

  invalidate_smiles();

  return rc;
}

int
Molecule::implicit_hydrogens()
{
  int rc = 0;
  
  for (int i = 0; i < _number_elements; i++) {
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
Molecule::set_implicit_hydrogens(atom_number_t i, int ih,
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

int
Molecule::GeometryIsSp2(atom_number_t zatom) const {
  const int acon = _things[zatom]->ncon();

  if (acon < 2) {
    return -1;
  }
  if (acon >= 4) {
    return 0;
  }

  // Make the decision to NOT check for unsaturation. This is
  // a geometry based method.
#ifdef CONSIDER_UNSATURATION
  if (_things[zatom]->nbonds() > acon) {
    return 0;
  }
#endif

  Set_of_Atoms conn;
  _things[zatom]->connections(zatom, conn);

  // If any bond angle is less than this, we assume `zatom` is tetrahedral
  static constexpr float kTetrahedral = 112.0f;

  // The angles must all be larger than this
  static constexpr float kFlat = 118.0f;

  const float a01 = bond_angle(conn[0], zatom, conn[1]) * RAD2DEG;
  if (a01 < kTetrahedral) {
    return 0;
  }
  if (acon == 2) {
    return a01 > kFlat;
  }

  const float a02 = bond_angle(conn[0], zatom, conn[2]) * RAD2DEG;
  if (a02 < kTetrahedral) {
    return 0;
  }

  const float a12 = bond_angle(conn[1], zatom, conn[2]) * RAD2DEG;
  if (a12 < kTetrahedral) {
    return 0;
  }

  if (a01 < kFlat || a02 < kFlat || a12 < kFlat) {
    return 0;
  }

  return 1;
}

#ifdef NOT_NEEDED_HANDLED_IN_EXISTING_FN
// Adding an H to something like *=C-*
int
Molecule::PlaceOneHydrogenSp2(const Make_Implicit_Hydrogens_Explicit& mihe) {
  atom_number_t anchor = mihe.a();

  const Atom * a = _things[anchor];

  atom_number_t o = a->other(anchor, 0);
  Coordinates v1(*a);
  v1 -= *_things[o];

  o = a->other(anchor, 1);
  Coordinates v2(*a);
  v2 =- *_things[o];

  v1.normalise();
  v2.normalise();

  v1 += v2;
  v1.normalise();
  v1 *= 1.09;
  v1 += *a;

  Atom * h = mihe.new_atom();
  h->set_x(v1.x());
  h->set_y(v1.y());
  h->set_z(v1.z());

  add(h);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);

#ifdef DEBUG_PLACEONEHYDROGENSP2
  cerr << "Anchor " << a->x() << ',' << a->y() << ',' << a->z() << '\n';
  const atom_number_t o0 = a->other(anchor, 0);
  const Atom* q = _things[o0];
  cerr << "0      " << q->x() << ',' << q->y() << ',' << q->z() << '\n';
  const atom_number_t o1 = a->other(anchor, 1);
  q = _things[o1];
  cerr << "1      " << q->x() << ',' << q->y() << ',' << q->z() << '\n';
  q = _things[_number_elements - 1];
  cerr << "H      " << q->x() << ',' << q->y() << ',' << q->z() << '\n';
  cerr << "dists " << distance_between_atoms(_number_elements - 1, anchor) << ' '
       << distance_between_atoms(_number_elements - 1, o0) << ' '
       << distance_between_atoms(_number_elements - 1, o1) << ' '
       << distance_between_atoms(anchor, o0) << ' '
       << distance_between_atoms(anchor, o1) << '\n';
#endif

  return 1;
}
#endif

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
  if (nullptr != c && 1 != c->implicit_hydrogen_count()) {
    cerr << "Molecule::_place_1_hydrogen: Hmmm, adding one H to chiral centre at atom " << anchor << endl;
    cerr << "Chiral centre is ";
    c->debug_print(cerr);

    return 0;
  }

  Atom * h = mihe.new_atom();

  add(h);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  if (nullptr != c && c->implicit_hydrogen_count()) {
    c->implicit_hydrogen_is_now_atom_number(_number_elements - 1);
  }

  const int dimensionality = mihe.dimensionality();
  // cerr << "_place_1_hydrogen dimensionality " << dimensionality << '\n';

  if (dimensionality < 2) {    // no need to worry about coordinates
    return 1;
  }

  if (2 == dimensionality || 0 == acon || acon > 3) {
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

  const Coordinates va = (*a);

// avoid straight bonds

//cerr << "Just about to put 1 Hydrogen, acon " << acon << " nbonds " << a->nbonds() << endl;

  coord_t blen = default_h_bond_length;
  if (6 != a->atomic_number()) {
    blen = static_cast<coord_t>(blen*0.9);
  }

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
Molecule::set_coordinates_of_singly_connected_atom(atom_number_t zatom,
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
    if (j == zatom) {
      continue;
    }

    Coordinates v = *(_things[j]);
  
    v -= va;
    v.normalise();
    vsum += v;
  }

  vsum.normalise();
  vsum *= default_bond_length;

  _things[zatom]->setxyz(a->x() - vsum.x(), a->y() - vsum.y(), a->z() - vsum.z());

// If 2 == acon, we want to avoid straight bonds unless there is a
// triple bond present

  if (2 == acon && 4 != a->nbonds()) { 
    atom_number_t j = a->other(anchor, 0);
    if (j == zatom) {
      j = a->other(anchor, 1);
    }

    set_bond_angle(j, anchor, zatom, static_cast<angle_t>(120.0 * DEG2RAD) );
  }

  return 1;
}

/*
  For now, two extra hydrogens is always one above and one below
*/

int
Molecule::_place_2_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe) {

  atom_number_t anchor = mihe.a();

  const int dimensionality = mihe.dimensionality();

  const Atom * a = _things[anchor]; 
  const int acon = a->ncon();

  Atom * h1 = mihe.new_atom();
  Atom * h2 = mihe.new_atom();

  add(h1);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  add(h2);
  add_bond(anchor, _number_elements - 1, SINGLE_BOND);

// take care of water before checking the dimensionality

  // Adding Hydrogens to water.
  if (0 == acon) {
    h1->setxyz(a->x() + default_h_bond_length, a->y(), a->z());
    coord_t dx = 0.95 * default_h_bond_length;
    coord_t dy = 0.31 * default_h_bond_length;
    h2->setxyz(a->x() - dx, a->y() + dy, a->z());

    return 1;
  }

  if (dimensionality < 2) {
    return 1;
  }

// The 2D case will never come up so just do something silly

  if (2 == dimensionality) {
    h1->setxyz(a->x(), a->y(), default_h_bond_length);
    h2->setxyz(a->x(), a->y(), - default_h_bond_length);

    return 1;
  }

  if (acon > 2 || 0 == acon) {
    cerr << "Molecule::_place_2_hydrogens: sorry, don't know how to place 2 Hydrogens on something with " <<
    acon << " existing connections\n"; return 0;
 }

// We need to set two vectors.  vsum we just add all the vectors
// departing from anchor and take the opposite direction vperp must be
// a vector normal to vsum

  Coordinates va = *(a);

  Coordinates v0;
  get_coords(other(anchor, 0), v0);
  v0 -= va;
  v0.normalise();

  Coordinates vsum;
  Coordinates vperp;

  if (2 == acon) {
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

    vperp.setxyz(static_cast<coord_t>(va.x() + 1.0),
                 static_cast<coord_t>(va.y() + 1.0),
                 static_cast<coord_t>(va.z() + 1.0));
    vperp.normalise();  // is this needed?
    vperp.cross_product(v0);
  }

  vperp.normalise();

  coord_t blen = default_h_bond_length;
  if (6 != a->atomic_number()) {
    blen *= static_cast<coord_t>(0.9);
  }

  float factor = static_cast<float>(blen / sqrt(2.0));
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
Molecule::set_coordinates_of_singly_connected_atoms(atom_number_t a1,
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
Molecule::_place_3_hydrogens(const Make_Implicit_Hydrogens_Explicit & mihe)
{
  atom_number_t anchor = mihe.a();
  const int dimensionality = mihe.dimensionality();

  const Atom * a = _things[anchor];
  const int acon = a->ncon();

  int first_h = _number_elements;

  Atom * atoms[3];
  for (int i = 0; i < 3; i++)
  {
    atoms[i] = mihe.new_atom();
    add(atoms[i]);
    add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  }

  if (dimensionality < 2) {
    return 1;
  }

  if (2 == dimensionality) {
    atoms[0]->setxyz(a->x(), a->y(), - default_h_bond_length);
    atoms[1]->setxyz(static_cast<coord_t>( a->x() - 0.25 * default_h_bond_length ), a->y(), default_h_bond_length);
    atoms[2]->setxyz(static_cast<coord_t>( a->x() + 0.25 * default_h_bond_length ), a->y(), default_h_bond_length);

    return 1;
  }

  coord_t blen = default_h_bond_length;
  if (6 != a->atomic_number())  {
    blen *= static_cast<coord_t>(0.9);
  }

  if (acon != 1) {
    cerr << "Molecule::_place_3_hydrogens:cannot place onto atom with " << acon <<
            " connections " << smarts_equivalent_for_atom(anchor) << '\n';
    // Should randomly place coords near the atom.
    return 0;
  }

  Set_of_Atoms h{first_h, first_h + 1, first_h + 2};

  const atom_number_t previous_atom = a->other(anchor, 0);

  // Place the added Hydrogens in the xy plane, centered on `a`.
  _things[first_h]->setxyz(0.9, 0.5196, 0.0);
  _things[first_h + 1]->setxyz(-0.9, 0.5196, 0.0);
  _things[first_h + 2]->setxyz(0, -1.0392, 0.0);
#ifdef DEBUG_PLACE_3_HYDROGENS
  cerr << "Btw hydrogens " << distance_between_atoms(first_h, first_h + 1) << '\n';
  cerr << "Btw hydrogens " << distance_between_atoms(first_h, first_h + 2) << '\n';
  cerr << "Btw hydrogens " << distance_between_atoms(first_h + 1, first_h + 2) << '\n';
#endif

  // Get a unit vector along the bond.

  Coordinates va = *a;
  Coordinates v0 = *(_things[previous_atom]);
  va -= v0;
  va.normalise();

  // Find a rotation axis that is perpendicular to the Z axis and the bond, and rotate around that.
  Coordinates vperp(0.0, 0.0, 1.0);
  const angle_t angle = va.angle_between_unit_vectors(vperp);
  vperp.cross_product(va);
  vperp.normalise();
  rotate_atoms(vperp, angle, h);
#ifdef DEBUG_PLACE_3_HYDROGENS
  cerr << "After rotation\n";
  cerr << "btw " << first_h << " and " << anchor << " " << distance_between_atoms(first_h, anchor) << '\n';
  cerr << "btw " << (first_h+1) << " and " << anchor << " "  << distance_between_atoms(first_h+1, anchor) << '\n';
  cerr << "btw " << (first_h+2) << " and " << anchor << " "  << distance_between_atoms(first_h+2, anchor) << '\n';
  cerr << vperp << '\n';
#endif

  Coordinates& h0 = *_things[first_h];
  h0 += *a;
  Coordinates& h1 = *_things[first_h + 1];
  h1 += *a;
  Coordinates& h2 = *_things[first_h + 2];
  h2 += *a;
  va *= 0.40;  // Empirically derived.
  h0 += va;
  h1 += va;
  h2 += va;
#ifdef DEBUG_PLACE_3_HYDROGENS
  cerr << "After tralsnation\n";
  cerr << "btw " << first_h << " and " << anchor << " " << distance_between_atoms(first_h, anchor) << '\n';
  cerr << "btw " << (first_h+1) << " and " << anchor << " "  << distance_between_atoms(first_h+1, anchor) << '\n';
  cerr << "btw " << (first_h+2) << " and " << anchor << " "  << distance_between_atoms(first_h+2, anchor) << '\n';
#endif
  return 1;
}

int
Molecule::_place_4_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe)
{
  const atom_number_t anchor = mihe.a();
  const int dimensionality = mihe.dimensionality();

  const int initial_matoms = _number_elements;

  const Atom * a = _things[anchor];

  Atom * atoms[4];
  for (int i = 0; i < 4; i++) {
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

  // Wikipedia
  _things[initial_matoms + 0]->setxyz(1.0, 0.0, -1.0 / sqrt(2.0));
  _things[initial_matoms + 1]->setxyz(-1.0, 0.0, -1.0 / sqrt(2.0));
  _things[initial_matoms + 2]->setxyz(0.0, 1.0, 1.0 / sqrt(2.0));
  _things[initial_matoms + 3]->setxyz(0.0, -1.0, 1.0 / sqrt(2.0));

  return 1;
}

int
Molecule::_place_lots_of_hydrogens (const Make_Implicit_Hydrogens_Explicit & mihe,
                                    int ih)
{
  atom_number_t anchor = mihe.a();
  const int dimensionality = mihe.dimensionality();

  resizable_array<Atom *> atoms;
  atoms.resize(ih);
  for (int i = 0; i < ih; i++)
  {
    Atom * h = mihe.new_atom();
    add(h);
    add_bond(anchor, _number_elements - 1, SINGLE_BOND);
  }

  if (dimensionality < 2) {
    return 1;
  }

  cerr << "Molecule::_place_lots_of_hydrogens:don't know how to place " << ih << " hydrogens\n";
  return 1;
}

int
Molecule::make_implicit_hydrogens_explicit()
{
  Make_Implicit_Hydrogens_Explicit mihe;

  return make_implicit_hydrogens_explicit(mihe);
}

int 
Molecule::make_implicit_hydrogens_explicit(atom_number_t a)
{
  Make_Implicit_Hydrogens_Explicit mihe;
  mihe.set_atom(a);

  return make_implicit_hydrogens_explicit(mihe);
}

int
Molecule::valence_ok(atom_number_t a)
{
  assert(ok_atom_number(a));

  return _things[a]->valence_ok();
}

int
Molecule::valence_ok()
{
  bool need_to_check_alternative_valences = false;
  for (int i = 0; i < _number_elements; i++)
  {
    if (! _things[i]->valence_ok()) {
      return 0;
    }

    if (_things[i]->element()->number_alternate_valences()) {
      need_to_check_alternative_valences = true;
    }
  }

  if (! need_to_check_alternative_valences) {
    return 1;
  }

  if (moleculeh::implicit_hydrogens_cannot_satisfy_alternate_valence &&
      moleculeh::ImplicitHydrogenSatisfiesAlternateValence(*this)) {
      return 0;
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
Molecule::make_implicit_hydrogens_explicit(Make_Implicit_Hydrogens_Explicit & mihe) {
  int initial_atoms = _number_elements;

  // Assume that the first molecule encountered indicates the dimensionality of
  // everything that follows... Almost always correct.
  if (mihe.dimensionality() < 0) {
    mihe.set_dimensionality(highest_coordinate_dimensionality());
  }

  if (const atom_number_t a = mihe.a(); a != INVALID_ATOM_NUMBER) {
    if (a < 0 || a >= _number_elements) {
      cerr << "Molecule::make_implicit_hydrogens_explicit: atom number " << a << " is invalid\n";
      return 0;
    }

    const int ih = implicit_hydrogens(a);
    // cerr << "ATom " << a << " has " << ih << " implicit hydrogens " << aromatic_smiles() << '\n';
    if (0 == ih) {
      return 1;
    }

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

  for (int i = 0; i < initial_atoms; ++i) {
    const int ih = implicit_hydrogens(i);
    if (ih == 0) {
      continue;
    }

    mihe.set_atom(i);

    int rc;

    if (1 == ih) {
      rc = _place_1_hydrogen(mihe);
    } else if (2 == ih) {
      rc = _place_2_hydrogens(mihe);
    } else if (3 == ih) {
      rc = _place_3_hydrogens(mihe);
    } else if (4 == ih) {
      rc = _place_4_hydrogens(mihe);
    } else {
      rc = _place_lots_of_hydrogens(mihe, ih);
    }

    _things[i]->set_implicit_hydrogens(0, 1);    // 1 means override known value
    if (rc == 0) {
      return 0;
    }
  }

  if (mihe.dimensionality() < 3) {
    return _number_elements - initial_atoms;
  }

  // A funbling effort to clear any major geometry problems. This is kind of
  // futile. Delete or fix...

  for (int i = initial_atoms; i < _number_elements; i++) {
    const Atom * ai = _things[i];
    // This test should never happen, they are all H atoms out here.
    if (1 != ai->atomic_number()) {
      continue;
    }

    atom_number_t j = ai->other(i, 0);
    Atom * aj = _things[j];

    if (ai->distance(*aj) < default_h_bond_length * 1.5) {
      continue;
    }

//  cerr << "Hydrogen " << i << " bonded to atom " << j << " distance " << ai->distance(*aj) << ". Coordinates changed\n";

    aj->setxyz(static_cast<coord_t>( ai->x() + 0.5 ),
                static_cast<coord_t>( ai->y() + 0.5 ),
                static_cast<coord_t>( ai->z() + 0.5) );
  }

  return _number_elements - initial_atoms;
}

#ifdef OLD_VERSION_MIHE
int
Molecule::make_implicit_hydrogens_explicit(Make_Implicit_Hydrogens_Explicit & mihe)
{
  // The initial number of atoms.
  int na = _number_elements;

  if (mihe.dimensionality() < 0) {
    mihe.set_dimensionality(highest_coordinate_dimensionality());
  }

  atom_number_t a = mihe.a();

  if (INVALID_ATOM_NUMBER != a) {
    if (a < 0 || a >= _number_elements) {
      cerr << "Molecule::make_implicit_hydrogens_explicit: atom number " << a << " is invalid\n";
      return 0;
    }

    int ih = implicit_hydrogens(a);
    // cerr << "ATom " << a << " " << smarts_equivalent_for_atom(a) << " has " << ih << " imp H\n";
    if (0 == ih) {
      return 1;
    }

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

  for (int i = 0; i < na ; i++) {
    mihe.set_atom(i);

    if (! make_implicit_hydrogens_explicit(mihe)) {    // recursive call
      rc = 0;
    }
  }

  for (int i = na; i < _number_elements; i++) {
    const Atom * ai = _things[i];
    if (1 != ai->atomic_number()) {
      continue;
    }

    atom_number_t j = ai->other(i, 0);
    Atom * aj = _things[j];

    if (ai->distance(*aj) < default_h_bond_length * 1.5) {
      continue;
    }

//  cerr << "Hydrogen " << i << " bonded to atom " << j << " distance " << ai->distance(*aj) << ". Coordinates changed\n";

    aj->setxyz(static_cast<coord_t>( ai->x() + 0.5 ),
                static_cast<coord_t>( ai->y() + 0.5 ),
                static_cast<coord_t>( ai->z() + 0.5) );
  }

  return rc;
}
#endif

// We try to preserve the ordering as much as possible.
// There are may already be Hydrogens at the end of the list of atoms.

int
Molecule::move_hydrogens_to_end_of_connection_table(atomic_number_t z)
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

int
Molecule::MoveToEndOfConnectionTable(atomic_number_t z) {
  int first_H = -1;
  int last_H = -1;
  int number_hydrogens = 0;
  for (int i = 0; i < _number_elements; ++i) {
    if (_things[i]->atomic_number() != z) {
      continue;
    }
    ++number_hydrogens;
    if (first_H < 0) {
      first_H = i;
    }
    last_H = i;
  }
  // If none found we are done.
  if (number_hydrogens == 0) {
    return 0;
  }
  // If all the H atoms are already at the end, we are done.
  if (_number_elements - first_H == number_hydrogens) {
    return 0;
  }
  // If just one found, easy...
  if (first_H == last_H) {
    move_atom_to_end_of_atom_list(first_H);
    return 1;
  }
  // Multiple Hydrogens are found, preserve their order.
  int * xref = new int[_number_elements]; std::unique_ptr<int[]> free_xref(xref);
  Atom ** original_order = new Atom *[_number_elements]; std::unique_ptr<Atom *[]> free_original_order(original_order);
  int heavy_atom_index = 0;
  int hydrogen_index = _number_elements - number_hydrogens;
  for (int i = 0; i < _number_elements; ++i) {
    original_order[i] = _things[i];
    if (_things[i]->atomic_number() == z) {
      xref[i] = hydrogen_index;
      ++hydrogen_index;
    } else {
      xref[i] = heavy_atom_index;
      ++heavy_atom_index;
    }
  }
  for (int i = 0; i < _number_elements; ++i) {
    _things[xref[i]] = original_order[i];
//  _things[i] = original_order[xref[i]];
  }
  for (Bond * b : _bond_list) {
    b->new_atom_numbers(xref);
  }
  for (int i = 0; i < _chiral_centres.number_elements(); i++) {
    _chiral_centres[i]->adjust_atom_numbers(xref);
  }
  // cis-trans bonds?
  _set_modified();

  return 1;
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

    if (! a->implicit_hydrogens_known()) {
      continue;
    }

    if (a->valence_ok()) {
      continue;
    }

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

    if (! a->implicit_hydrogens_known()) {
      continue;
    }

    if (a->isotope() || a->formal_charge()) {
      continue;
    }

    if (a->element()->organic())
      ;
    else if (5 == a->element()->atomic_number())    // Boron is OK
      ;
    else
      continue;

    int ih;
    if (! a->compute_implicit_hydrogens(ih)) {
      continue;
    }

    if (ih != a->implicit_hydrogens()) {
      continue;
    }

    if (! valence_ok(i)) {   // being very careful
      continue;
    }

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
