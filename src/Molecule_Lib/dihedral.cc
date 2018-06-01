/*
  Fucntions for the Dihedral_Atoms class.
*/

#include <iostream>
#include <iostream>
#include <assert.h>
#include <memory>
using std::cerr;
using std::endl;

#include "misc.h"

#include "molecule.h"

angle_t
Molecule::dihedral_angle (atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4,
                          int not_bonded_ok) const
{
  assert (ok_4_atoms(a1, a2, a3, a4));

  if (! not_bonded_ok)
  {
    assert (are_bonded(a1, a2));
    assert (are_bonded(a2, a3));
    assert (are_bonded(a3, a4));
  }

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];
  const Atom * aa3 = _things[a3];
  const Atom * aa4 = _things[a4];

  return angle_between_atoms(*aa1, *aa2, *aa3, *aa4);
}

//#define DEBUG_SET_DIHEDRAL

int
Molecule::set_dihedral(atom_number_t a1, 
                       atom_number_t a2,
                       atom_number_t a3,
                       atom_number_t a4,
                       angle_t theta)
{
  assert (ok());
  assert (natoms() >= 4);

  assert (ok_4_atoms(a1, a2, a3, a4));

  angle_t current_angle = dihedral_angle(a1, a2, a3, a4, 1);

  if (current_angle == theta)     // nothing to do
    return 1;

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "Molecule::set_dihedral:atoms " << a1 << ',' << a2 << ',' << a3 << ',' << a4 << " angle " << theta << " (" << (theta * RAD2DEG) << " deg)\n";
  cerr << "current value " << current_angle << " (" << (current_angle * RAD2DEG) << " deg)\n";
#endif

  int * atoms_to_move = new_int(_number_elements); std::unique_ptr<int[]> free_atoms_to_move(atoms_to_move);

  atoms_to_move[a2] = 2;
  atoms_to_move[a3] = 1;

  const Atom * aa3 = _things[a3];

  const int acon = aa3->ncon();

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "A3 has " << acon << " connections\n";
#endif

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = aa3->other(a3, i);

#ifdef DEBUG_SET_DIHEDRAL
    cerr << " bonded to " << j << ' ' << smarts_equivalent_for_atom(j) << endl;
#endif

    if (j == a2)
      continue;

    if (0 == _determine_moving_atoms(j, atoms_to_move))
    {
      cerr << "Molecule::set_dihedral:possible ring structure, a1 " << a2 << " a2 " << a2 << " a3 " << a3 << " a4 " << a4 << endl;
      return 0;
    }
  }

  atoms_to_move[a3] = 0;

  double rot = current_angle - theta;

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "Rotating " << rot << " (" << (rot * RAD2DEG) << " deg)\n";
#endif

  double r = _things[a2]->distance(*(_things[a3]));
  if (0.0 == r)
  {
#ifdef DEBUG_SET_DIHEDRAL
    cerr << "Molecule::set_dihedral:zero distance between atoms " << a2 << " and " << a3 << endl;
    _things[a2]->debug_print(cerr);
    _things[a3]->debug_print(cerr);
    write_molecule_mdl(cerr, "");
#endif
    return 0;
  }

//  The direction cosines of the a2-a3 bond

  double dc1 = (_things[a3]->x() - _things[a2]->x()) / r;
  double dc2 = (_things[a3]->y() - _things[a2]->y()) / r;
  double dc3 = (_things[a3]->z() - _things[a2]->z()) / r;

  double rotmat11 = cos(rot) + dc1 * dc1 * (1.0 - cos(rot));
  double rotmat12 = dc1 * dc2 * (1.0 - cos(rot)) - dc3 * sin(rot);
  double rotmat13 = dc1 * dc3 * (1.0 - cos(rot)) + dc2 * sin(rot);
  double rotmat21 = dc1 * dc2 * (1.0 - cos(rot)) + dc3 * sin(rot);
  double rotmat22 = cos(rot) + dc2 * dc2 * (1.0 - cos(rot));
  double rotmat23 = dc2 * dc3 * (1.0 - cos(rot)) - dc1 * sin(rot);
  double rotmat31 = dc3 * dc1 * (1.0 - cos(rot)) - dc2 * sin(rot);
  double rotmat32 = dc3 * dc2 * (1.0 - cos(rot)) + dc1 * sin(rot);
  double rotmat33 = cos(rot) + dc3 * dc3 * (1.0 - cos(rot));

  double x0 = _things[a3]->x();
  double y0 = _things[a3]->y();
  double z0 = _things[a3]->z();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == atoms_to_move[i])
      continue;

#ifdef DEBUG_SET_DIHEDRAL
    cerr << "Molecule::set_dihedral::moving atom " << i << endl;
#endif

    Atom *a = _things[i];

    double xx = a->x() - x0;
    double yy = a->y() - y0;
    double zz = a->z() - z0;
    double xnew = rotmat11 * xx + rotmat12 * yy + rotmat13 * zz + x0;
    double ynew = rotmat21 * xx + rotmat22 * yy + rotmat23 * zz + y0;
    double znew = rotmat31 * xx + rotmat32 * yy + rotmat33 * zz + z0;

    a->setxyz( static_cast<coord_t>(xnew), static_cast<coord_t>(ynew), static_cast<coord_t>(znew) );

//  cerr << i << " from (" << xold << ',' << yold << ',' << zold << ") to (" << a->x() << ',' << a->y() << ',' << a->z() << ")\n";
  }

#ifdef DEBUG_SET_DIHEDRAL
  angle_t tmp = dihedral_angle(a1, a2, a3, a4);
  cerr << "Molecule::set_dihedral:after moving " << tmp << " (" << (tmp * RAD2DEG) << " deg)\n";
#endif
  
  return 1;
}

int
Molecule::bump_check (atom_number_t a1,
                      atom_number_t a2,
                      atom_number_t a3,
                      atom_number_t a4,
                      distance_t d) const
{
  int * either_side = new int[_number_elements]; std::unique_ptr<int[]> free_either_side(either_side);

  return _bump_check(a1, a2, a3, a4, d, either_side);
}

int
Molecule::_bump_check (atom_number_t a1,
                       atom_number_t a2,
                       atom_number_t a3,
                       atom_number_t a4,
                       distance_t too_close, int * either_side) const
{
  if (! _determine_either_side_of_bond(a2, a3, either_side))
    return 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == either_side[i])     // one of the atoms forming the bond
      continue;

    const Atom * ai = _things[i];

    for (int j = i + 1; j < _number_elements; j++)
    {
      if (either_side[i] == either_side[j])    // on same side of the bond
        continue;

      if (0 == either_side[j])     // one of the atoms forming the bond
        continue;

      distance_t dij = ai->distance(*(_things[j]));

      if (dij <= too_close)
      {
        cerr << "Molecule::_bump_check: atoms " << i << " and " << j << " too close " << dij << " min is " << too_close << endl;
        return 0;
      }
    }
  }

  return 1;
}

int
Molecule::_determine_either_side_of_bond (atom_number_t a1,
                                          atom_number_t a2,
                                          int * either_side) const
{
#ifdef DEBUG_DETERMINE_EITHER_SIDE_OF_BOND
  cerr << "_determine_either_side_of_bond with atoms " << a1 << " and " << a2 << endl;
#endif

  set_vector(either_side, _number_elements, -9);

  identify_side_of_bond(either_side, a1, -1, a2);

#ifdef DEBUG_DETERMINE_EITHER_SIDE_OF_BOND
  cerr << "Ring? " << (-1 == either_side[a2]) << endl;
#endif

  if (-1 == either_side[a2])      // bond must be in a ring
    return 0;

  identify_side_of_bond(either_side, a2, 1, a1);    // will definitely work

  either_side[a1] = 0;
  either_side[a2] = 0;

  return 1;
}

/*
  There are two identify_side_of_bond functions. One is called at the
  top level, and here is it OK if it encounters the atom at the end of
  the initial bond. In the second instance, it is a fatal error to 
  encounter that atom.
*/

int
Molecule::identify_side_of_bond (int * either_side,
                                 atom_number_t astart,
                                 int flag,
                                 atom_number_t avoid) const
{
  int rc = 0;

  const Atom * a = _things[astart];

  either_side[astart] = flag;

  const int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(astart, i);

    if (j == avoid)     // just looking at the other end of our starting bond
      continue;

#ifdef DEBUG_IDENTIFY_SIDE_OF_BOND
    cerr << "Top level, astart " << astart << " to " << j << " must avoid " << avoid << ", flag " << flag << endl;
#endif

    int tmp = __identify_side_of_bond(either_side, j, flag, avoid);
    if (0 == tmp)
      return 0;

    rc += tmp;
  }

#ifdef DEBUG_IDENTIFY_SIDE_OF_BOND
  cerr << "From atom " << astart << " returning " << rc << endl;
#endif

  return rc;
}

int
Molecule::__identify_side_of_bond (int * either_side,
                                   atom_number_t astart,
                                   int flag,
                                   atom_number_t avoid) const
{
#ifdef DEBUG_IDENTIFY_SIDE_OF_BOND
  cerr << "Molecule::__identify_side_of_bond: at atom " << astart << " type [" << _things[astart]->atomic_symbol() << 'D' << _things[astart]->ncon() << "], flag " << flag << "\n";
#endif

  either_side[astart] = flag;

  int rc = 1;

  const Atom * a = _things[astart];

  for (int i = 0; i < a->ncon(); i++)
  {
    atom_number_t j = a->other(astart, i);

//  cerr << "Molecule::__identify_side_of_bond:how about atom " << j << endl;

    if (flag == either_side[j])    // already done this atom - a ring somewhere
      continue;

    if (avoid == j)      // we have doubled back to the other end of the initial bond - must have been in a ring
      return 0;

    int tmp = __identify_side_of_bond(either_side, j, flag, avoid);
    if (0 == tmp)
      return 0;

    rc += tmp;
  }

  return rc;
}

/*
  This determines whether or not it is possible to perform a 
  conformational search around a given bond. If either of the
  two atoms which define the twisting bond have only one connection,
  then a conformational search will not be possible.
*/

/*int
Molecule::conf_search_possible (const Dihedral_Atoms & dihedral) const
{
  assert (ok());
  if (! valid_dihedral(dihedral))
    return 0;

  atom_number_t a2 = dihedral.a2();
  atom_number_t a3 = dihedral.a3();

  const Atom *aa2 = atomi(a2);
  const Atom *aa3 = atomi(a3);

  if (1 == aa2->ncon())
    return 0;
  if (1 == aa3->ncon())
    return 0;

  return 1;
}*/

int
Molecule::bump_check(const int * tmp,
                     distance_t tolerance) const
{
  int rc = 1;

  for (int i = 0; i < _number_elements; i++)
  {
    for (int j = i + i; j < _number_elements; j++)
    {
      if (tmp[i] == tmp[j])    // same grouping, don't compare them
        continue;

      if (_things[i]->distance(*(_things[j])) < tolerance)
        rc = 0;
    }
  }

  return rc;
}

//#define DEBUG_SET_BOND_ANGLE

int
Molecule::_determine_moving_atoms(atom_number_t zatom,
                                  int * moving_atom) const
{
  moving_atom[zatom] = 1;

  int rc = 1;

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

//  if (2 == moving_atom[j])    // gack, came back to the non-moving atom, a1
//    cerr << "Gack, from " << zatom << " got to " << j << endl;
    if (2 == moving_atom[j])    // gack, came back to the non-moving atom, a1
      return 0;

    if (moving_atom[j])
      continue;

    int tmp = _determine_moving_atoms(j, moving_atom);
    if (0 == tmp)
      return 0;

    rc += tmp;
  }

  return rc;
}

/*
  New implementation, oct 2005.
  Need to handle several cases:
    a1 outside a ring, a2 and a3 in a ring
    multiple bonds from a2
*/

int
Molecule::set_bond_angle (atom_number_t a1,
                          atom_number_t a2,
                          atom_number_t a3,
                          angle_t theta)
{
  assert (ok_3_atoms(a1, a2, a3));

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];
  const Atom * aa3 = _things[a3];

  if (! aa1->is_bonded_to(a2))
  {
    cerr << "Molecule::set_bond_angle:atom " << a1 << " not bonded to " << a2 << endl;
    return 0;
  }

  if (! aa2->is_bonded_to(a3))
  {
    cerr << "Molecule::set_bond_angle:atom " << a2 << " not bonded to " << a3 << endl;
    return 0;
  }

  if (aa1->is_bonded_to(a3))    // /yipes a 3 membered ring
  {
    cerr << "Molecule::set_bond_angle:atom " << a2 << " bonded to " << a3 << " 3 membered ring!\n";
    return 0;
  }

// If ring info is available use it. Don't force a ring determination.

  if (_sssr_rings.number_elements() && in_same_ring(a1, a3))
  {
    cerr << "Molecule::set_bond_angle:cannot change ring bond, atoms " << a2 << " and " << a3 << endl;
    return 0;
  }

  int * moving_atoms = new_int(_number_elements); std::unique_ptr<int[]> free_moving_atoms(moving_atoms);

  moving_atoms[a1] = 2;    // special flag - if this value is encountered, in _determine_moving_atoms, we abort
  moving_atoms[a2] = 1;

  const int acon = aa2->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = aa2->other(a2, i);

    if (a1 == j)
      continue;

    if (! _determine_moving_atoms(j, moving_atoms))
    {
      cerr << "Molecule::set_bond_angle:cannot identify atoms to move, atoms " << a1 << ',' << a2 << ',' << a3 << "\n";
      return 0;
    }
  }

  moving_atoms[a1] = 0;    // atom a1 does not move
  moving_atoms[a2] = 0;    // atom a2 does not move

#ifdef DEBUG_SET_BOND_ANGLE
  for (int i = 0; i < _number_elements; i++)
  {
    if (moving_atoms[i])
      cerr << "Atom " << i << " moving_atoms " << moving_atoms[i] << ' ' << smarts_equivalent_for_atom(i) << endl;
  }
#endif

  Space_Vector<double> v12(aa1->x() - aa2->x(), aa1->y() - aa2->y(), aa1->z() - aa2->z()); //  = *aa1 - *aa2;
  v12.normalise();

  Space_Vector<double> v32(aa3->x() - aa2->x(), aa3->y() - aa2->y(), aa3->z() - aa2->z()); // *aa3 - *aa2;
  v32.normalise();

  double current_angle = v12.angle_between_unit_vectors(v32);

#ifdef DEBUG_SET_BOND_ANGLE
  cerr << "Angle currently " << current_angle << " (" << (current_angle * RAD2DEG) << " deg), desired angle " << theta << " (" << (theta * RAD2DEG) << " deg)\n";
  write_molecule_mdl(cerr, "Within dihedral routine");
#endif

// We'll use v32 as the axis

  if (fabs(current_angle - M_PI) < 0.01)    // yipes almost a straight line, V12 can be anything
    v12.setxyz(v32.y(), v32.z(), v32.x());

// The cross product should be normal to the plane.

  v32.cross_product(v12);

// But if there was a 90 degree angle between the two vectors, the
// cross product will be zero

  if (v32.norm() < 1.0e-03)
  {
    v32.setxyz(aa3->x() - aa2->x(), aa3->y() - aa2->y(), aa3->z() - aa2->z());
    v32.normalise();

    v32 += v12;
    v32.normalise();
    v32.cross_product(v12);
  }

  v32.normalise();

  double rot = current_angle - theta;

#ifdef DEBUG_SET_BOND_ANGLE
  cerr << "Rotating " << rot << " (" << (rot * RAD2DEG) << " degrees) around " << v32 << endl;
#endif

  Set_of_Atoms atoms_to_move;
  atoms_to_move.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    if (1 == moving_atoms[i])
      atoms_to_move.add(i);
  }

  if (0 == atoms_to_move.size())
  {
    cerr << "Molecule::set_bond_angle:no atoms being moved!\n";
    return 0;
  }

#ifdef DEBUG_SET_BOND_ANGLE
  cerr << "Moving " << atoms_to_move.size() << " atoms around " << v32 << endl;
  cerr << atoms_to_move << endl;
#endif

// Need to shift so that atom A2 is the origin

  coord_t a2x = _things[a2]->x();
  coord_t a2y = _things[a2]->y();
  coord_t a2z = _things[a2]->z();

  translate_atoms(-a2x, -a2y, -a2z, atoms_to_move);

  int rc = rotate_atoms(v32, rot, atoms_to_move);

  translate_atoms(a2x, a2y, a2z, atoms_to_move);

#ifdef DEBUG_SET_BOND_ANGLE
  angle_t tmp = bond_angle(a1, a2, a3);
  cerr << "After moving angle is " << tmp << " or " << (RAD2DEG * tmp) << " degrees\n";
#endif

  return rc;
}
