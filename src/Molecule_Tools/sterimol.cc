/*
  Compute sterimol parameters

  The strategy is to align the longest axis of the molecule along the X
  axis. Then compute the distance of each atom from the X axis. Then, those
  atoms with a chance of being the most distant are in turn, placed in the
  X/Y plane and the atom with the largest spread in the X/Y plane identified.
  Once that if done, we go back to that atom to do a computation which
  includes averages.
*/

#include "sterimol.h"

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/iw_vdw.h"
#include "Molecule_Lib/molecule.h"

using std::cerr;

static angle_t rotation_angle = 5.0 * M_PI / 180.0;

angle_t
default_sterimol_rotation_angle() {
  return rotation_angle;
}

void
set_sterimol_rotation_angle(angle_t a) {
  rotation_angle = a;

  return;
}

static int default_vdw_type = IW_VDW_SAVOL;

/*
  Once we have checked all atoms, we move the molecule so that the longest
  extremity is in the X/Y plane, we can compute the average displacements
  and determine the angle extrema
*/

static void
compute_average_offsets(const Molecule& m, const int* process_atom, const Atom** atom,
                        const vdw_radius_t* vdw, Sterimol& b) {
  const coord_t zero = static_cast<coord_t>(0.0);  // try to minimise casts

  Accumulator<coord_t> yave, zave;

  coord_t xmin = zero;
  atom_number_t atom_with_shortest_x_coord = INVALID_ATOM_NUMBER;  // not on the X axis

  coord_t highest_yx_ratio = zero;
  coord_t highest_zx_ratio = zero;

  coord_t ymax = zero;
  coord_t x_at_ymax = zero;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (!process_atom[i]) {
      continue;
    }

    const Atom* a = atom[i];

    assert(a->x() >= zero);

    coord_t y = a->y();
    if (y < zero) {
      y = -y;
    }

    coord_t z = a->z();
    if (z < zero) {
      z = -z;
    }

    yave.extra(y + vdw[i]);
    zave.extra(z + vdw[i]);

    if (zero == y || zero == z) {  // doesn't subtend an angle to the X axis
      continue;
    }

    if (y + vdw[i] > ymax) {
      ymax = y + vdw[i];
      x_at_ymax = a->x();
    }

    coord_t r = (y + vdw[i]) / a->x();
    if (r > highest_yx_ratio) {
      highest_yx_ratio = r;
    }

    r = (z + vdw[i]) / a->x();
    if (r > highest_zx_ratio) {
      highest_zx_ratio = r;
    }

    if (zero == a->x()) {  // in the Y/Z plane
      ;
    } else if (INVALID_ATOM_NUMBER == atom_with_shortest_x_coord)  // first atom checked
    {
      atom_with_shortest_x_coord = i;
      xmin = a->x();
    } else if (a->x() < xmin) {
      atom_with_shortest_x_coord = i;
      xmin = a->x();
    }
  }

  b[STERIMOL_AVE_Y] = yave.average();
  b[STERIMOL_AVE_Z] = zave.average();

  angle_t theta = atan(highest_yx_ratio);

  b[STERIMOL_Y_ANGLE] = theta * RAD2DEG;

  theta = atan(highest_zx_ratio);

  b[STERIMOL_Z_ANGLE] = theta * RAD2DEG;

  b[STERIMOL_ANGLE_RATIO] = b[STERIMOL_Y_ANGLE] / b[STERIMOL_Z_ANGLE];

  if (zero != x_at_ymax) {  // could it ever be 0.0?
    b[STERIMOL_OVERALL_CONE] = RAD2DEG * atan(ymax / x_at_ymax);
  }

  if (INVALID_ATOM_NUMBER == atom_with_shortest_x_coord) {
    return;
  }

  const Atom* a = atom[atom_with_shortest_x_coord];

  double ryx = fabs(a->y() / a->x());
  double rzx = fabs(a->z() / a->x());

  if (ryx > rzx) {
    b[STERIMOL_LARGE_CONE] = atan(ryx) * RAD2DEG;
  } else {
    b[STERIMOL_LARGE_CONE] = atan(rzx) * RAD2DEG;
  }

  return;
}

/*
  The molecule has been rotated to a given alignment. Look for the
  longest Z coordinate.
*/

static void
sterimol_config(Molecule& m, const Atom** atom, const vdw_radius_t* vdw, Sterimol& b) {
  // We can initialise these to 0.0 because we know the molecule is along the X axis

  coord_t ymax = 0.0;
  coord_t ymin = 0.0;
  coord_t zmax = 0.0;
  coord_t zmin = 0.0;


  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    const Atom* a = atom[i];

    if (a->y() + vdw[i] > ymax) {
      ymax = a->y() + vdw[i];
    }

    if (a->y() - vdw[i] < ymin) {
      ymin = a->y() - vdw[i];
    }

    if (a->z() + vdw[i] > zmax) {
      zmax = a->z() + vdw[i];
    }

    if (a->z() - vdw[i] < zmin) {
      zmin = a->z() - vdw[i];
    }
  }

  coord_t ycross = (ymax - ymin);
  coord_t zcross = (zmax - zmin);

  if (fabs(ymax) > fabs(ymin)) {
    b[STERIMOL_LONGEST_PERP_DIST] = fabs(ymax);
  } else {
    b[STERIMOL_LONGEST_PERP_DIST] = fabs(ymin);
  }

  b[STERIMOL_YCROSS] = ycross;
  b[STERIMOL_ZCROSS] = zcross;

  b[STERIMOL_YZPROD] = ycross * zcross;

  return;
}

static Coordinates x_axis(1.0, 0.0, 0.0);
static Coordinates y_axis(0.0, 1.0, 0.0);

/*
  Rotate the molecule so that atom A is in the X/Y plane (z = 0)
*/

static void
rotate_atom_to_z_zero(Molecule& m, const Atom* a) {
  Coordinates tmp(0.0, a->y(), a->z());  // project into Y/Z plane
  tmp.normalise();

  angle_t theta = tmp.angle_between_unit_vectors(y_axis);

  // cerr << "X = " << a->x () << " Y = " << a->y() << " Z = " << a->z() << " angle is "
  // << theta << '\n';

  if (a->z() < 0.0) {
    m.rotate_atoms(x_axis, theta);
  } else if (a->z() > 0.0) {
    m.rotate_atoms(x_axis, -theta);
  }

  if (fabs(a->z()) >= 0.005) {
    cerr << "Yipes, Z coordinate bad " << a->x() << ' ' << a->y() << ' ' << a->z()
         << '\n';
  }
  assert(fabs(a->z()) < 0.005);

  return;
}

// #define WRITE_INTERMEDIATE_FILES

#ifdef WRITE_INTERMEDIATE_FILES
static int files_written = 0;
#endif

static int
sterimol(Molecule& m, const int* process_atom, const Atom** atom, vdw_radius_t* vdw,
         coord_t* d, Sterimol& b) {
  if (!assign_vdw_radii(m, default_vdw_type, vdw)) {
    return 0;
  }

  atom_number_t left = INVALID_ATOM_NUMBER;
  atom_number_t right = INVALID_ATOM_NUMBER;

  coord_t maxd = 0.0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (!process_atom[i]) {
      continue;
    }

    const Atom* ai = atom[i];

    for (int j = i + 1; j < matoms; j++) {
      if (!process_atom[j]) {
        continue;
      }

      //    don't exclude this in case just a single atom is the substituent
      //    if (ai->is_bonded_to (j))
      //      continue;

      const Atom* aj = atom[j];

      coord_t d = ai->distance(*aj);

      if (d > maxd) {
        maxd = d;
        left = i;
        right = j;
      }
    }
  }

  assert(INVALID_ATOM_NUMBER != left);

  b[STERIMOL_LONGEST_DISTANCE] = maxd;

  // Translate the molecule so that atom LEFT is at the origin

  m.translate_atoms(-*(atom[left]));

  assert(0.0 == atom[left]->x() && 0.0 == atom[left]->y() && 0.0 == atom[left]->z());

  // We now want to rotate the molecule so that RIGHT is along the X axis

  Coordinates x(1.0, 0.0, 0.0);
  Coordinates r(atom[right]->x(), atom[right]->y(), atom[right]->z());
  r.normalise();
  angle_t theta = x.angle_between_unit_vectors(r);

  r.cross_product(x);
  r.normalise();
  m.rotate_atoms(r, theta);

#ifdef WRITE_INTERMEDIATE_FILES
  // cerr << "Left = " << left << " right = " << right << '\n';
  IWString fname = "foo";
  fname += files_written;
  fname += ".sdf";
  IWString comment;
  comment << "Left = " << left << " right = " << right;
  m.write_molecule_mdl(fname, comment.null_terminated_chars());
#endif

  // At this stage, the molecule should have LEFT at the origin and RIGHT somwhere along
  // the X axis

  assert(fabs(atom[right]->y()) < 0.002);
  assert(fabs(atom[right]->z()) < 0.002);

  // Determine the furthest atom from the X axis

  coord_t dmax = static_cast<coord_t>(0.0);

  for (int i = 0; i < matoms; i++) {
    if (!process_atom[i]) {
      continue;
    }

    if (i == left) {
      continue;
    }

    const Atom* a = atom[i];

    d[i] = sqrt((a->y() * a->y()) + (a->z() * a->z())) + vdw[i];
    if (d[i] > dmax) {
      dmax = d[i];
    }
  }

  // Only those atoms which extend >= 0.5 dmax need be considered

  coord_t dmax2 = dmax * 0.5;

  // For each atom being processed, rotate the molecule so that atom is along the Y axis

  atom_number_t extremity_atom = left;

  for (int i = 0; i < matoms; i++) {
    if (!process_atom[i]) {
      continue;
    }

    if (i == left) {
      continue;
    }

    if (d[i] <= dmax2) {
      continue;
    }

    const Atom* a = atom[i];

    rotate_atom_to_z_zero(m, a);  // put atom A in X/Y plane

    //  cerr << "Atom " << i << " moved to " << a->x () << ',' << a->y () << ',' << a->z() << '\n';

#ifdef WRITE_INTERMEDIATE_FILES
    IWString fname = "foo";
    fname += files_written;
    fname += ".sdf";
    IWString comment;
    comment << "Rotated to atom " << i;
    m.write_molecule_mdl(fname, comment.null_terminated_chars());
#endif

    assert(0.0 == atom[left]->x() && 0.0 == atom[left]->y() &&
           0.0 == atom[left]->z());  // left must be at the origin

    assert(0.0 != atom[right]->x() && fabs(atom[right]->y()) < 0.002 &&
           fabs(atom[right]->z()) < 0.002);

    Sterimol bb;

    sterimol_config(m, atom, vdw, bb);

    if (b.copy_if_longer(bb))  // returns 1 if updated successfully
    {
      extremity_atom = i;
      b[STERIMOL_LONGEST_DISTANCE] =
          maxd;  // doesn't change as the molecule is rotated and isn't set above
    }
  }

  assert(INVALID_ATOM_NUMBER != extremity_atom);

  rotate_atom_to_z_zero(m, atom[extremity_atom]);

  compute_average_offsets(m, process_atom, atom, vdw, b);

  return 1;
}

int
sterimol(Molecule& m, Sterimol& b) {
  b.reset();

  int matoms = m.natoms();

  if (matoms < 3) {
    cerr << "Sterimol parameters not defined for " << matoms << " atoms\n";
    return 0;
  }

  const Atom** atoms = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atoms(atoms);

  m.atoms(atoms);

  vdw_radius_t* vdw = new vdw_radius_t[matoms];
  std::unique_ptr<vdw_radius_t[]> free_vdw(vdw);

  coord_t* d = new coord_t[matoms];
  std::unique_ptr<coord_t[]> free_d(d);  // distance from X axis

  int* tmp = new_int(m.natoms(), 1);
  std::unique_ptr<int[]> free_tmp(tmp);

  return sterimol(m, tmp, atoms, vdw, d, b);
}

Sterimol::Sterimol() {
  _do_partial_charge_descriptors = 0;

  reset();

  return;
}

void
Sterimol::reset() {
  set_vector(_b, NSTERIMOL, static_cast<sterimol_t>(0.0));

  set_vector(_partial_charge_descriptor, NSTERIMOL_PARTIAL_CHARGE,
             static_cast<charge_t>(0.0));
  _from_base_atom.reset();

  return;
}

sterimol_t
Sterimol::operator[](int i) const {
  assert(i >= 0 && i < NSTERIMOL);

  return _b[i];
}

sterimol_t&
Sterimol::operator[](int i) {
  assert(i >= 0 && i < NSTERIMOL);

  return _b[i];
}

/*
  A frequent operation is to compare two Sterimol configurations and save
  the largest one. We base our comparison on the longest perpendicular distance
*/

int
Sterimol::copy_if_longer(const Sterimol& rhs) {
  if (rhs._b[STERIMOL_LONGEST_PERP_DIST] <= _b[STERIMOL_LONGEST_PERP_DIST]) {
    return 0;  // no change, rhs is shorter/smaller
  }

  copy_vector(_b, rhs._b, NSTERIMOL);

  // Note that we don't copy the _from_base_atom accumulator

  return 1;
}

int
Sterimol::do_computation(Molecule& m, const int* don_acc, const int* process,
                         atom_number_t a1) {
  reset();

  int matoms = m.natoms();

  if (matoms < 3) {
    cerr << "Sterimol parameters not defined for " << matoms << " atoms\n";
    return 0;
  }

  _atoms_in_fragment = count_non_zero_occurrences_in_array(process, matoms);

  if (_atoms_in_fragment < 2) {
    cerr << "Sterimol::do_computation: only " << _atoms_in_fragment
         << " atoms flagged as being in fragment, impossible\n";
    return 0;
  }

  const Atom** atoms = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atoms(atoms);

  m.atoms(atoms);

  vdw_radius_t* vdw = new vdw_radius_t[matoms];
  std::unique_ptr<vdw_radius_t[]> free_vdw(vdw);

  coord_t* d = new coord_t[matoms];
  std::unique_ptr<coord_t[]> free_d(d);  // distance from X axis

  if (INVALID_ATOM_NUMBER != a1) {
    _compute_distances_from_a1(matoms, process, a1, atoms);
  }

  return sterimol(m, process, atoms, vdw, d, *this);
}

/*
  We need to provide two arrays.
  An array of donor acceptor values, and an array of which atoms to process
*/

int
Sterimol::do_computation(Molecule& m) {
  const int matoms = m.natoms();

  int* tmp = new int[matoms + matoms];
  std::unique_ptr<int[]> free_tmp(tmp);

  std::fill_n(tmp, matoms, 0);           // donor acceptor
  std::fill_n(tmp + matoms, matoms, 1);  // atoms to proces

  return do_computation(m, tmp, tmp + matoms);
}

int
Sterimol::write_descriptors(std::ostream& os) const {
  for (int i = 0; i < NSTERIMOL; i++) {
    os << ' ' << _b[i];
  }

  if (_from_base_atom.n()) {
    os << ' ' << _from_base_atom.maxval();

    if (_from_base_atom.n() > 1) {
      os << ' ' << _from_base_atom.average();
    } else {
      os << ' ' << _from_base_atom.minval();
    }
  }

  if (_do_partial_charge_descriptors) {
    for (int i = 0; i < NSTERIMOL_PARTIAL_CHARGE; i++) {
      os << ' ' << _partial_charge_descriptor[i];
    }
  }

  return os.good();
}

int
Sterimol::_compute_distances_from_a1(int matoms, const int* process, atom_number_t a1,
                                     const Atom** atom) {
  assert(0 == _from_base_atom.n());

  const Atom* aa1 = atom[a1];

  for (int i = 0; i < matoms; i++) {
    if (!process[i]) {
      continue;
    }

    if (a1 == i) {
      continue;
    }

    if (2 == _atoms_in_fragment) {
      ;
    } else if (aa1->is_bonded_to(i)) {
      continue;
    }

    coord_t d = aa1->distance(*(atom[i]));

    _from_base_atom.extra(d);
  }

  return 1;
}

int
write_sterimol_descriptor_headers(std::ostream& os, const IWString& prefix,
                                  int include_distance_from_base_atom) {
  os << ' ' << prefix << "longd";
  os << ' ' << prefix << "perpd";
  os << ' ' << prefix << "ycross";
  os << ' ' << prefix << "zcross";
  os << ' ' << prefix << "yave";
  os << ' ' << prefix << "zave";
  os << ' ' << prefix << "yzprod";

  os << ' ' << prefix << "yangle";
  os << ' ' << prefix << "zangle";
  os << ' ' << prefix << "yzrtio";
  os << ' ' << prefix << "lgcone";
  os << ' ' << prefix << "ovcone";

  if (include_distance_from_base_atom) {
    os << ' ' << prefix << "bsmaxd" << ' ' << prefix << "bsaved";
  }

  return os.good();
}
