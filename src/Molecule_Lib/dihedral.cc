/*
  Fucntions for the Dihedral_Atoms class.
*/

#include <assert.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "molecule.h"

using std::cerr;

angle_t
Molecule::dihedral_angle(atom_number_t a1, atom_number_t a2, atom_number_t a3, atom_number_t a4,
                         BondedStatus bonded_status) const
{
  assert(ok_4_atoms(a1, a2, a3, a4));

  if (bonded_status == BondedStatus::kMustBeBonded) {
    if (! are_bonded(a1, a2) ||
        ! are_bonded(a2, a3) ||
        ! are_bonded(a3, a4)) {
      cerr << "Molecule::dihedral_angle:atoms not bonded\n";
      cerr << "a1 " << a1 << " a2 " << a2 << " bonded " << are_bonded(a1, a2) << '\n';
      cerr << "a2 " << a2 << " a3 " << a3 << " bonded " << are_bonded(a2, a3) << '\n';
      cerr << "a3 " << a3 << " a4 " << a4 << " bonded " << are_bonded(a3, a4) << '\n';
      return 0.0;
    }
  }

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];
  const Atom * aa3 = _things[a3];
  const Atom * aa4 = _things[a4];

  return angle_between_atoms(*aa1, *aa2, *aa3, *aa4);
}

angle_t
Molecule::signed_dihedral_angle(atom_number_t a1, atom_number_t a2,
                                atom_number_t a3, atom_number_t a4) const {
  assert(ok_4_atoms(a1, a2, a3, a4));
  const Atom* aa1 = _things[a1];
  const Atom* aa2 = _things[a2];
  const Atom* aa3 = _things[a3];
  const Atom* aa4 = _things[a4];

  Space_Vector<double> v21(aa1->x() - aa2->x(), aa1->y() - aa2->y(), aa1->z() - aa2->z());
  Space_Vector<double> v32(aa2->x() - aa3->x(), aa2->y() - aa3->y(), aa2->z() - aa3->z());
  Space_Vector<double> v43(aa3->x() - aa4->x(), aa3->y() - aa4->y(), aa3->z() - aa4->z());

  v21.normalise();
  v32.normalise();
  v43.normalise();

  v21.cross_product(v32);
  v43.cross_product(v32);
  v43.negate();

  // The unsigned result.
  double result = v21.angle_between(v43);
#ifdef DEBUG_SIGNED_DIHEDRAL
  cerr << "Unsigned result " << result << '\n';
#endif

  // Use notation from https://mathworld.wolfram.com/Plane.html
  // v21 is normal to the plane defined by a1, a2, a3.
  // The coefficients for the equation for the a1,a2,a3 plane.

  const double a = v21.x();
  const double b = v21.y();
  const double c = v21.z();

  // s1 always evaluates to near zero. By construction, it is
  // in the plane, so only consider s2.
  // double s1 = a * aa1->x() + b * aa1->y() + c * aa1->z();

  double s2 = a * aa4->x() + b * aa4->y() + c * aa4->z();

  // When we have things that are close to planar, s2 may assume
  // very small values, and become unstable. Introduce a convention that if
  // the value is small, we return a positive result.
#ifdef DEBUG_SIGNED_DIHEDRAL
  cerr << "Wrt plane " << s2 << '\n';
#endif
  if (abs(s2) < std::numeric_limits<float>::epsilon()) {
    return result;
  }

  if (s2 < 0.0) {
    return result;
  } else {
    return -result;
  }
}

//#define DEBUG_SET_DIHEDRAL

int
Molecule::set_dihedral(atom_number_t a1, 
                       atom_number_t a2,
                       atom_number_t a3,
                       atom_number_t a4,
                       angle_t theta)
{
  assert(ok());
  assert(natoms() >= 4);

  assert(ok_4_atoms(a1, a2, a3, a4));

  angle_t current_angle = signed_dihedral_angle(a1, a2, a3, a4);

  if (abs(current_angle - theta) < 0.005) {     // nothing to do 
    return 1;
  }

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "Molecule::set_dihedral:atoms " << a1 << ',' << a2 << ',' << a3 << ',' << a4 << " angle " << theta << " (" << (theta * RAD2DEG) << " deg)\n";
  cerr << "current value " << current_angle << " (" << (current_angle * RAD2DEG) << " deg)\n";
#endif

  int * atoms_to_move = new_int(_number_elements); std::unique_ptr<int[]> free_atoms_to_move(atoms_to_move);

  atoms_to_move[a2] = 2;
  atoms_to_move[a3] = 1;

  const Atom * aa3 = _things[a3];

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "A3 has " << aa3->ncon() << " connections " << smarts_equivalent_for_atom(a3) << '\n';
  cerr << smiles() << '\n';
#endif

  for (const Bond* b : *aa3) {
    const atom_number_t j = b->other(a3);

    if (j == a2)
      continue;

    if (0 == _determine_moving_atoms(j, atoms_to_move))
    {
      cerr << "Molecule::set_dihedral:possible ring structure, a1 " << a2 << " a2 " << a2 << " a3 " << a3 << " a4 " << a4 << '\n';
      return 0;
    }
  }

  atoms_to_move[a3] = 0;

  double rot = current_angle - theta;

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "Rotating " << rot << " (" << (rot * RAD2DEG) << " deg)\n";
  write_isotopically_labelled_smiles(*this, false, cerr);
  cerr << '\n';
#endif

  const double r = _things[a2]->distance(*(_things[a3]));
  if (0.0 == r)
  {
#ifdef DEBUG_SET_DIHEDRAL
    cerr << "Molecule::set_dihedral:zero distance between atoms " << a2 << " and " << a3 << '\n';
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

#ifdef DEBUG_SET_DIHEDRAL
  cerr << rotmat11 << ' ' << rotmat12 << ' ' << rotmat13 << '\n';
  cerr << rotmat21 << ' ' << rotmat22 << ' ' << rotmat23 << '\n';
  cerr << rotmat31 << ' ' << rotmat32 << ' ' << rotmat33 << '\n';
#endif

  double x0 = _things[a3]->x();
  double y0 = _things[a3]->y();
  double z0 = _things[a3]->z();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == atoms_to_move[i])
      continue;

#ifdef DEBUG_SET_DIHEDRAL
    cerr << "Molecule::set_dihedral::moving atom " << i << '\n';
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
  angle_t tmp = signed_dihedral_angle(a1, a2, a3, a4);
  cerr << "Molecule::set_dihedral:after moving " << tmp << " (" << (tmp * RAD2DEG) << " deg)\n";
#endif
  
  return 1;
}

// Return true if there are no pairs of atoms with different atoms_to_move values
// that are within `bump_check`.
// Atoms `a1` and `a2` are the central bond, so ignore them in any distance calculation
static int
OkBumpCheck(const Molecule& m, const int* atoms_to_move,
            atom_number_t a1, atom_number_t a2, float bump_check) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == a1 || i == a2) {
      continue;
    }
    for (int j = i + 1; j < matoms; ++j) {
      if (j == a1 || j == a2) {
        continue;
      }
      if (atoms_to_move[i] == atoms_to_move[j]) {
        continue;
      }
      const float d = m.distance_between_atoms(i, j);
      if (d > bump_check) {
        continue;
      }
      if (m.are_bonded(i, j)) {
        continue;
      }
      return 0;
    }
  }

  // No violations found.
  return 1;
}

std::vector<std::unique_ptr<float[]>>
Molecule::DihedralScan(atom_number_t a2,
                atom_number_t a3, angle_t delta, float bump_check) {
  std::unique_ptr<float[]> starting_coords = GetCoords();
  
  delta = std::abs(delta);
  if (delta <= 1.0f || delta >= 180.0f) {
    cerr << "Molecule::DihedralScan:invalid delta " << delta << "\n";
    return std::vector<std::unique_ptr<float[]>>();
  }

  // Remember, we do not return the current coordinates, so decrement by 1.
  const int nconf = static_cast<int>(360.0f / delta) - 1;
  assert(nconf >= 0);

  // We could guard against a delta of 44.5 degrees, which see the last
  // conformation very similar to the starting values. Assume well behaved.

  std::vector<std::unique_ptr<float[]>> result;
  result.reserve(nconf);

  if (nconf <= 1) {  // cannot happen, paranoid...
    cerr << "Molecule::DihedralScan:no conformers for " << delta << " degrees\n";
    return result;
  }

  std::unique_ptr<int[]> atoms_to_move(new_int(_number_elements));

  // Thought about making a function to consolidate with set_dihedral, but
  // I think the extra complexity is not worth it.

  atoms_to_move[a2] = 2;
  atoms_to_move[a3] = 1;

  const Atom * aa3 = _things[a3];

#ifdef DEBUG_SET_DIHEDRAL
  cerr << "A3 has " << aa3->ncon() << " connections " << smarts_equivalent_for_atom(a3) << '\n';
  cerr << smiles() << '\n';
#endif

  for (const Bond* b : *aa3) {
    const atom_number_t j = b->other(a3);

    if (j == a2)
      continue;

    if (0 == _determine_moving_atoms(j, atoms_to_move.get())) {
      cerr << "Molecule::DihedralScan:possible ring structure, a1 " << a2 << " a2 " << a2 << " a3 " << a3 << '\n';
      result.resize(0);
      return result;
    }
  }

  atoms_to_move[a3] = 0;

  const double r = _things[a2]->distance(*(_things[a3]));

//  The direction cosines of the a2-a3 bond

  const double dc1 = (_things[a3]->x() - _things[a2]->x()) / r;
  const double dc2 = (_things[a3]->y() - _things[a2]->y()) / r;
  const double dc3 = (_things[a3]->z() - _things[a2]->z()) / r;

  const double rot = delta * DEG2RAD;

  const double rotmat11 = cos(rot) + dc1 * dc1 * (1.0 - cos(rot));
  const double rotmat12 = dc1 * dc2 * (1.0 - cos(rot)) - dc3 * sin(rot);
  const double rotmat13 = dc1 * dc3 * (1.0 - cos(rot)) + dc2 * sin(rot);
  const double rotmat21 = dc1 * dc2 * (1.0 - cos(rot)) + dc3 * sin(rot);
  const double rotmat22 = cos(rot) + dc2 * dc2 * (1.0 - cos(rot));
  const double rotmat23 = dc2 * dc3 * (1.0 - cos(rot)) - dc1 * sin(rot);
  const double rotmat31 = dc3 * dc1 * (1.0 - cos(rot)) - dc2 * sin(rot);
  const double rotmat32 = dc3 * dc2 * (1.0 - cos(rot)) + dc1 * sin(rot);
  const double rotmat33 = cos(rot) + dc3 * dc3 * (1.0 - cos(rot));

#ifdef DEBUG_SET_DIHEDRAL
  cerr << rotmat11 << ' ' << rotmat12 << ' ' << rotmat13 << '\n';
  cerr << rotmat21 << ' ' << rotmat22 << ' ' << rotmat23 << '\n';
  cerr << rotmat31 << ' ' << rotmat32 << ' ' << rotmat33 << '\n';
#endif

  const double x0 = _things[a3]->x();
  const double y0 = _things[a3]->y();
  const double z0 = _things[a3]->z();

  for (int i = 0; i < nconf; ++i) {
    for (int j = 0; j < _number_elements; ++j) {
      if (atoms_to_move[j] == 0) {
        continue;
      }

      Atom *a = _things[j];

      double xx = a->x() - x0;
      double yy = a->y() - y0;
      double zz = a->z() - z0;
      double xnew = rotmat11 * xx + rotmat12 * yy + rotmat13 * zz + x0;
      double ynew = rotmat21 * xx + rotmat22 * yy + rotmat23 * zz + y0;
      double znew = rotmat31 * xx + rotmat32 * yy + rotmat33 * zz + z0;

      a->setxyz(static_cast<coord_t>(xnew), static_cast<coord_t>(ynew), static_cast<coord_t>(znew));
    }

    if (bump_check > 0.0f && ! OkBumpCheck(*this, atoms_to_move.get(), a2, a3, bump_check)) {
      continue;
    }

    result.emplace_back(GetCoords());
#ifdef DEBUG_DIHEDRAL_SCAN
    cerr << "In c++ these are the coodinates\n";
    for (int q = 0; q < _number_elements; ++q) {
      cerr << ' ' << result.back()[q * 3 + 0];
      cerr << ' ' << result.back()[q * 3 + 1];
      cerr << ' ' << result.back()[q * 3 + 2];
    }
    cerr << '\n';
#endif
  }

  SetXyz(starting_coords.get());

  return result;
}

int
Molecule::IncrementDihedral(atom_number_t a2, atom_number_t a3, angle_t angle) {
  assert(ok_2_atoms(a2, a3));

  const atom_number_t a1 = _things[a2]->other(a2, 0);
  const atom_number_t a4 = _things[a3]->other(a3, 0);

  const angle_t current_angle = signed_dihedral_angle(a1, a2, a3, a4);

  return set_dihedral(a1, a2, a3, a4, current_angle + angle);
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
        cerr << "Molecule::_bump_check: atoms " << i << " and " << j << " too close " << dij << " min is " << too_close << '\n';
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
  cerr << "_determine_either_side_of_bond with atoms " << a1 << " and " << a2 << '\n';
#endif

  set_vector(either_side, _number_elements, -9);

  identify_side_of_bond(either_side, a1, -1, a2);

#ifdef DEBUG_DETERMINE_EITHER_SIDE_OF_BOND
  cerr << "Ring? " << (-1 == either_side[a2]) << '\n';
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
    cerr << "Top level, astart " << astart << " to " << j << " must avoid " << avoid << ", flag " << flag << '\n';
#endif

    int tmp = __identify_side_of_bond(either_side, j, flag, avoid);
    if (0 == tmp)
      return 0;

    rc += tmp;
  }

#ifdef DEBUG_IDENTIFY_SIDE_OF_BOND
  cerr << "From atom " << astart << " returning " << rc << '\n';
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

//  cerr << "Molecule::__identify_side_of_bond:how about atom " << j << '\n';

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

int
Molecule::bump_check(distance_t dist) const {
  int rc = 0;
  for (int i = 0; i < _number_elements; ++i) {
    for (int j = i + 1; j < _number_elements; ++j) {
      if (are_bonded(i, j)) {
        continue;
      }
      if (_things[i]->distance(*_things[j]) < dist) {
        ++rc;
      }
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
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);

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
Molecule::set_bond_angle(atom_number_t a1,
                         atom_number_t a2,
                         atom_number_t a3,
                         angle_t theta)
{
  assert (ok_3_atoms(a1, a2, a3));

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];
  const Atom * aa3 = _things[a3];

  if (! aa1->is_bonded_to(a2)) {
    cerr << "Molecule::set_bond_angle:atom " << a1 << " not bonded to " << a2 << '\n';
    return 0;
  }

  if (! aa2->is_bonded_to(a3)) {
    cerr << "Molecule::set_bond_angle:atom " << a2 << " not bonded to " << a3 << '\n';
    return 0;
  }

  if (aa1->is_bonded_to(a3))    // /yipes a 3 membered ring
  {
    cerr << "Molecule::set_bond_angle:atom " << a2 << " bonded to " << a3 << " 3 membered ring!\n";
    return 0;
  }

// If ring info is available use it. Don't force a ring determination.

  if (_sssr_rings.number_elements() && in_same_ring(a1, a3)) {
    cerr << "Molecule::set_bond_angle:cannot change ring bond, atoms " << a2 << " and " << a3 << '\n';
    return 0;
  }

#ifdef DEBUG_SET_BOND_ANGLE
  cerr << "set_bond_angle atoms " << a1 << ' ' << a2 << ' ' << a2 << ' ' << a3 << " angle " << theta << '\n';
#endif

  std::unique_ptr<int[]> moving_atoms = std::make_unique<int[]>(_number_elements);
  std::fill_n(moving_atoms.get(), _number_elements, 0);

  // special flag - if this value is encountered, in _determine_moving_atoms fail.
  moving_atoms[a1] = 2;
  moving_atoms[a2] = 1;

  if (! _determine_moving_atoms(a3, moving_atoms.get())) {
    cerr << "Molecule::set_bond_angle:cannot identify moving atoms " << a1 << ' ' << a2 << ' ' << a3 << '\n';
    return 0;
  }

#ifdef DEBUG_SET_BOND_ANGLE
  for (int i = 0; i < _number_elements; i++)
  {
    if (moving_atoms[i]) {
      cerr << "Atom " << i << " moving_atoms " << moving_atoms[i] << ' ' << smarts_equivalent_for_atom(i) << '\n';
    }
  }
#endif

  //  Generate *aa1 - *aa2;
  Space_Vector<double> v12(aa1->x() - aa2->x(), aa1->y() - aa2->y(), aa1->z() - aa2->z());
  v12.normalise();

  // Generate *aa3 - *aa2;
  Space_Vector<double> v32(aa3->x() - aa2->x(), aa3->y() - aa2->y(), aa3->z() - aa2->z());
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

  if (v32.norm() < 1.0e-03) {
    v32.setxyz(aa3->x() - aa2->x(), aa3->y() - aa2->y(), aa3->z() - aa2->z());
    v32.normalise();

    v32 += v12;
    v32.normalise();
    v32.cross_product(v12);
  }

  v32.normalise();

  double rot = current_angle - theta;

#ifdef DEBUG_SET_BOND_ANGLE
  cerr << "Rotating " << rot << " (" << (rot * RAD2DEG) << " degrees) around " << v32 << '\n';
#endif

  Set_of_Atoms atoms_to_move;
  atoms_to_move.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    if (1 == moving_atoms[i])
      atoms_to_move.add(i);
  }

  if (atoms_to_move.empty()) {
    cerr << "Molecule::set_bond_angle:no atoms being moved!\n";
    return 0;
  }

#ifdef DEBUG_SET_BOND_ANGLE
  cerr << "Moving " << atoms_to_move.size() << " atoms around " << v32 << '\n';
  cerr << atoms_to_move << '\n';
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

angle_t
DihedralAngle(const Atom* a1, const Atom* a2, const Atom* a3, const Atom* a4) {
  return angle_between_atoms(a1, a2, a3, a4);
}

int
Molecule::LocationOfSubstituent(atom_number_t zatom, distance_t dist,
                      Coordinates& result) const {
  if (_number_elements == 0) {
    return 0;
  }

  if (_number_elements == 1) {
    result = *_things[0];
    return 0;
  }

  const int acon = _things[zatom]->ncon();

  if (acon == 1) {
    return LocationOfSubstituent1(zatom, dist, result);
  } 
  if (acon == 2) {
    return LocationOfSubstituent2(zatom, dist, result);
  } 
  if (acon == 3) {
    return LocationOfSubstituent3(zatom, dist, result);
  } 

  return 0;
}

int
Molecule::LocationOfSubstituent1(atom_number_t zatom, distance_t dist,
                      Coordinates& result) const {
  const Atom* a1 = _things[zatom];
  assert(a1->ncon() == 1);

  const atom_number_t o = a1->other(zatom, 0);

  const Atom* a2 = _things[o];

  result = *a1;
  result -= *a2;
  result.normalise();
  result *= dist;
  result += *a1;

  return 1;
}

int
Molecule::LocationOfSubstituent2(atom_number_t zatom, distance_t dist,
                      Coordinates& result) const {
  const Atom* a1 = _things[zatom];
  assert(a1->ncon() == 2);

  const atom_number_t o1 = a1->other(zatom, 0);
  const atom_number_t o2 = a1->other(zatom, 1);

  const Atom* a2 = _things[o1];
  const Atom* a3 = _things[o2];

  result = *a1;
  result -= *a2;
  result.normalise();
  Coordinates v3(*a1);
  v3 -= *a3;
  v3.normalise();

  result += v3;
  result.normalise();
  result *= dist;
  result += *a1;

  return 1;
}

int
Molecule::LocationOfSubstituent3(atom_number_t zatom, distance_t dist,
                      Coordinates& result) const {
  const Atom* a1 = _things[zatom];
  assert(a1->ncon() == 3);

  const atom_number_t o1 = a1->other(zatom, 0);
  const atom_number_t o2 = a1->other(zatom, 1);
  const atom_number_t o3 = a1->other(zatom, 2);

  const Atom* a2 = _things[o1];
  const Atom* a3 = _things[o2];
  const Atom* a4 = _things[o3];

  result = *a1;
  result -= *a2;
  result.normalise();
  Coordinates v3(*a1);
  v3 -= *a3;
  v3.normalise();
  Coordinates v4(*a1);
  v4 -= *a4;
  v4.normalise();

  result += v3;
  result += v4;
  result.normalise();
  result *= dist;
  result += *a1;

  return 1;
}

int
Position3D(const Molecule& m1, atom_number_t atom1,
                 float distance,
                 Molecule& m2, atom_number_t atom2) {
  Coordinates ext1, ext2;
  if (! m1.LocationOfSubstituent(atom1, distance, ext1) ||
      ! m2.LocationOfSubstituent(atom2, distance, ext2)) {
    cerr << "Position3D:cannot get location for substituent\n";
    return 0;
  }
  const Atom& a2 = m2[atom2];

  // We must orient `m2` so it is pointing towards `m1`.
  // First translate it so `a2` is at `ext1`.
  m2.translate_atoms(ext1 - a2);

  // Redetermine ext2.
  m2.LocationOfSubstituent(atom2, distance, ext2);
  ext2.cross_product(ext1);
  angle_t angle = ext1.angle_between(ext2);
  m2.rotate_atoms(ext1, -angle);

  return 1;
}


namespace lillymol {

//#define DEBUG_3D_POSITIONING

// try to position the fragment defined by `atom2` so it can join with
// `atom1`. `hydrogen1` and `hydrogen2` are coordinats of what would be
// explicit Hydrogens attached to `atom1` and `atom2` and they define
//  the vectors where the two atoms can grow.
float
TryPlacing(Molecule& m,
           const int initial_natoms,
           atom_number_t atom1,
           const Coordinates& hydrogen1,
           atom_number_t atom2,
           const Coordinates hydrogen2,
           float distance,
           const int* fragment_membership) {
  const Atom& a1 = m[atom1];
  const Atom& a2 = m[atom2];

#ifdef DEBUG_3D_POSITIONING
  cerr << " h1dist " << m[atom1].distance(hydrogen1) << " h2dist " << m[atom2].distance(hydrogen2) << '\n';
#endif

  int moving_fragment = fragment_membership[atom2];

  Coordinates v1 = hydrogen1 - a1;
  Coordinates v2 = hydrogen2 - a2;

  // translate the moving fragment to the origin.
  m.translate_atoms(-a2, fragment_membership, moving_fragment);

  const angle_t angle = v1.angle_between(v2);
#ifdef DEBUG_3D_POSITIONING
  cerr << "Initial angle " << (angle * RAD2DEG) << '\n';
#endif

  v1.normalise();
  v2.normalise();
  v2.cross_product(v1);

  m.rotate_atoms(v2, -static_cast<coord_t>(M_PI - angle),
                 fragment_membership, moving_fragment);

  // Now translate the moving fragment back to where atom1 wants it to be
  // what about distance...
  m.translate_atoms(hydrogen1, fragment_membership, moving_fragment);

#ifdef DEBUG_3D_POSITIONING
  cerr << "atom 2 now at " << m.x(atom2) << ',' << m.y(atom2) << ',' << m.z(atom2) << '\n';
  cerr << "Distance between " << m.distance_between_atoms(atom1, atom2) << '\n';
  v1 = hydrogen1 - a1;
  v2 = hydrogen2 - a2;
  cerr << "Angle now " << (v1.angle_between(v2) * RAD2DEG) << '\n';
#endif

  // Now set the bond length.
  float current_distance = m.distance_between_atoms(atom1, atom2);
  // The vector along which we will move.
  v1 = a2 - a1;
  v1.normalise();
  v1 *= (distance - current_distance);
  m.translate_atoms(v1, fragment_membership, moving_fragment);

  // Return the shortest distance between fragments.
  // Note that we do not specifically check that the atoms are
  // in fragment_membership[atom1] and fragment_membership[atom2]
  // since it would slow things down for likely no benefit.
  float rc = std::numeric_limits<float>::max();
  for (int i = 0; i < initial_natoms; ++i) {
    for (int j = i + 1; j < initial_natoms; ++j) {
      if (fragment_membership[i] == fragment_membership[j]) {
        continue;
      }
      float d = m.distance_between_atoms(i, j);
      if (d < rc) {
        rc = d;
      }
    }
  }

  return rc;
}

#ifdef NO_LONGER_NEEDED_BUT_GOOD_WAY_OF_DOING_THIS
// Remove any explicit Hydrogens attached to `zatom`.
// If any atom number below `other` is removed, decrement `other`.
int
RemoveExplicitHydrogens(Molecule& m,
                        atom_number_t zatom,
                        atom_number_t& other) {
  int rc = 0;
  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (m.atomic_number(i) != 1) {
      continue;
    }
    if (! m.are_bonded(zatom, i)) {
      continue;
    }

    m.remove_atom(i);
    if (i < other) {
      --other;
    }
    ++rc;
  }

  return rc;
}
#endif

struct Position3DAtom {
  atom_number_t atom;
  int aromatic;
  int sp2;
};

int
AddExplicitHydrogens(Molecule& m,
                     Position3DAtom& atom1,
                     Position3DAtom& atom2,
                     Set_of_Atoms& hydrogens) {
  atom_number_t zatom = atom1.atom;

  const int initial_matoms = m.natoms();

  if (! m.make_implicit_hydrogens_explicit(zatom)) {
    cerr << "AddExplicitHydrogens:make_implicit_hydrogens_explicit failed\n";
    return 0;
  }

  // cerr << "Molecule starts with " << initial_matoms << " now has " << m.natoms() << '\n';
  for (int i = initial_matoms; i < m.natoms(); ++i) {
    m.set_atomic_number(i, 9);
    hydrogens << i;
  }
  // m.write_molecule_mdl<std::ostream>(cerr, "");

  return hydrogens.size();
}

// Return an array of the coordindates of the atoms in `atoms`.
static std::unique_ptr<Coordinates[]>
AtomsToCoordinates(const Molecule& m, const Set_of_Atoms& atoms) {
  const int n = atoms.size();

  std::unique_ptr<Coordinates[]> result = std::make_unique<Coordinates[]>(n);
  for (uint32_t i = 0; i < atoms.size(); ++i) {
    result[i] = m[atoms[i]];
  }

  return result;
}

int
Position3DInner(Molecule& m,
                atom_number_t atom1, float distance, atom_number_t atom2) {
  const int initial_natoms = m.natoms();

  Position3DAtom pa1, pa2;
  pa1.atom = atom1;
  pa1.aromatic = m.is_aromatic(atom1);
  pa1.sp2 = m.GeometryIsSp2(atom1);

  pa2.atom = atom2;
  pa2.aromatic = m.is_aromatic(atom2);
  pa2.sp2 = m.GeometryIsSp2(atom2);

  Set_of_Atoms hydrogen1, hydrogen2;

  if (! AddExplicitHydrogens(m, pa1, pa2, hydrogen1)) {
    cerr << "Position3D:cannot place explicit Hydrogen atom " <<
            m.smarts_equivalent_for_atom(atom1) << '\n';
    return 0;
  }
  if (! AddExplicitHydrogens(m, pa2, pa1, hydrogen2)) {
    cerr << "Position3D:cannot place explicit Hydrogen atom " << 
            m.smarts_equivalent_for_atom(atom2) << '\n';
    return 0;
  }

  std::unique_ptr<Coordinates[]> hcoords1 = AtomsToCoordinates(m, hydrogen1);
  std::unique_ptr<Coordinates[]> hcoords2 = AtomsToCoordinates(m, hydrogen2);

  // Done with the explicit Hydrogens added.
  m.resize(initial_natoms);

#ifdef DEBUG_POSITION3D
  for (atom_number_t h : hydrogen1) {
    cerr << " h1 " << m.distance_between_atoms(atom1, h) << '\n';
    for (const Bond* b : m[atom1]) {
      atom_number_t o = b->other(atom1);
      cerr << h << ' ' << atom1 << ' ' << o << " angle " << (m.bond_angle(h, atom1, o) * RAD2DEG) << '\n';
    }
  }

  for (atom_number_t h : hydrogen2) {
    cerr << " h2 " << m.distance_between_atoms(atom2, h) << '\n';
    for (const Bond* b : m[atom1]) {
      atom_number_t o = b->other(atom1);
      cerr << h << ' ' << atom1 << ' ' << o << " angle " << (m.bond_angle(h, atom2, o) * RAD2DEG) << '\n';
    }
  }
#endif

  std::unique_ptr<int[]> fragment_membership = m.fragment_membership();

  if (fragment_membership[atom1] == fragment_membership[atom2]) {
    cerr << "Position3D:atoms in same fragment, cannot process " << atom1 << ' ' << atom2 << '\n';
    return 0;
  }

  std::unique_ptr<float[]> starting_coords = m.GetCoordinates();

  // On each of `atom1` and `atom2` we have implicit Hydrogen atom positions.
  // We need to decide which combination to use.
  // TryPlacing will return the closest approach between atoms in the two
  // fragments, and we want this to be as large as possible.
  float furthest_separation = 0.0f;
  std::unique_ptr<float[]> best_coords;

  const int nh1 = hydrogen1.size();
  const int nh2 = hydrogen2.size();

  for (int i = 0; i < nh1; ++i) {
    for (int j = 0; j < nh2; ++j) {
      float score = TryPlacing(m, initial_natoms, atom1, hcoords1[i], atom2, hcoords2[j],
                               distance, fragment_membership.get());
      if (score > furthest_separation) {
        furthest_separation = score;
        best_coords = m.GetCoordinates();
      }

      m.SetCoordinates(starting_coords.get());
    }
  }

  m.SetCoordinates(best_coords.get());

  return 1;
}

int
Position3D(Molecule& m, atom_number_t atom1, float distance, atom_number_t atom2) {
  if (! m.ok_2_atoms(atom1, atom2)) {
    cerr << "Position3D:invalid atom numbers " << atom1 << ' ' << atom2 << '\n';
    return 0;
  }

  if (m.hcount(atom1) == 0 || m.hcount(atom2) == 0) {
    cerr << "Position3D:No Hydrogen atoms on " << atom1 << ' ' << 
            m.smarts_equivalent_for_atom(atom1) << " or " << atom2 << 
            m.smarts_equivalent_for_atom(atom2) << '\n';
    return 0;
  }

#ifdef DEBUG_POSITION3D
  cerr << "in c++ " << m.smiles() << " atoms " << atom1 << ' ' << m.smarts_equivalent_for_atom(atom1) << " and " << m.smarts_equivalent_for_atom(atom2) << '\n';
#endif

  return Position3DInner(m, atom1, distance, atom2);
}

}  // namespace lillymol
