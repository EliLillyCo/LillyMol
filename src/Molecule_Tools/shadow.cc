/*
  Determine the projection of a molecule onto the Y/Z plane

  We do this by allocating a grid and projecting each atom onto that grid.
*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>

#include "misc.h"

#include "shadow.h"

static double default_grid_resolution = 0.1;

Shadow_Area::Shadow_Area ()
{
  _resolution = default_grid_resolution;

  _ny = _nz = 0;

  _grid = NULL;

  return;
}

Shadow_Area::~Shadow_Area ()
{
  if (NULL != _grid)
    delete [] _grid;

  _grid = NULL;

  return;
}

int
Shadow_Area::ok () const
{
  if (_resolution <= 0.0)
    return 0;

  if (NULL == _grid)
    return 1;

  if (_ny <= 0 || _nz <= 0)
    return 0;

  return 1;
}

/*
  Note, this isn't useful - see the constructor
*/

int
Shadow_Area::active () const
{
  return _resolution > 0.0;
}

int
Shadow_Area::_nearest_y_grid_point (coord_t y) const
{
  double tmp = (y - _ymin) / _resolution;

  if (tmp < 0.0)
    return 0;

  tmp = rint (tmp);

  int rc = static_cast<int> (rint (tmp) + 0.001);

  if (rc >= _ny)
    rc = _ny - 1;

  return rc;
}

area_t
Shadow_Area::shadow_area (const Molecule & m,
                          const vdw_radius_t * vdw)
{
  assert (ok ());

  if (NULL != _grid)
  {
    delete [] _grid;
    _grid = NULL;
  }

  int matoms = m.natoms ();

  if (0 == matoms)
    return 0.0;

  const Atom ** atom = new const Atom * [matoms];

  m.atoms (atom);

  area_t rc = _shadow_area (matoms, atom, vdw);

  delete [] atom;

  return rc;
}

//#define DEBUG_SHADOW_AREA

#ifdef DEBUG_SHADOW_AREA
#include "accumulator.h"
#endif

area_t 
Shadow_Area::_shadow_area (int matoms,
                           const Atom * const * atom,
                           const vdw_radius_t vdw [])
{
  assert (matoms > 0);

  const Atom * a = atom[0];

// Determine the extremeties.

#ifdef DEBUG_SHADOW_AREA
  Accumulator<float> X, Y, Z;
  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = atom[i];
    X.extra (a->x ());
    Y.extra (a->y ());
    Z.extra (a->z ());
  }

  cerr << "X between " << X.minval () << " and " << X.maxval () << endl;
  cerr << "Y between " << Y.minval () << " and " << Y.maxval () << endl;
  cerr << "Z between " << Z.minval () << " and " << Z.maxval () << endl;
#endif

  _ymin = a->y () - vdw[0];
  _zmin = a->z () - vdw[0];

  double ymax = a->y () + vdw[0];
  double zmax = a->z () + vdw[0];

  for (int i = 1; i < matoms; i++)
  {
    const Atom * a = atom[i];

    if (a->y () + vdw[i] > ymax)
      ymax = a->y () + vdw[i];
    if (a->y () - vdw[i] < _ymin)
      _ymin = a->y ()- vdw[i];

    if (a->z () + vdw[i] > zmax)
      zmax = a->z () + vdw[i];
    if (a->z () - vdw[i] < _zmin)
      _zmin = a->z ()- vdw[i];
  }

  _ny = static_cast<int> ((ymax - _ymin) / _resolution) + 2;
  _nz = static_cast<int> ((zmax - _zmin) / _resolution) + 2;

  _grid = new_int (_ny * _nz);

#ifdef DEBUG_SHADOW_AREA
  cerr << matoms << " atoms between y = (" << _ymin << ',' << ymax << ") and z = (" << _zmin << ',' << zmax << ")\n";
  cerr << "_ny = " << _ny << " _nz = " << _nz << " product " << (_ny * _nz) << endl;
#endif

  _project_atoms_to_grid (matoms, atom, vdw);

  area_t rc = _determine_area ();

  return rc;
}

area_t
Shadow_Area::_determine_area () const
{
  int n = _ny * _nz;

  int rc = 0;      // the number of grid points set

  for (int i = 0; i < n; i++)
  {
    assert (_grid[i] >= 0);

    if (_grid[i])
      rc++;
  }

//return static_cast<area_t> (rc);

  return _resolution * _resolution * static_cast<area_t> (rc);
}

void
Shadow_Area::_project_atoms_to_grid (int matoms,
                                     const Atom * const * atom,
                                     const vdw_radius_t * vdw)
{
  for (int i = 0; i < matoms; i++)
  {
    _project_atom (atom[i], vdw[i]);
  }

  return;
}

/*
  At first I tried to write a unified project atom, which would work
  regardless of the orientation of the grid. Turned out to be too hard
  so I gave up and we now have a different function for each.
*/

void
Shadow_Area::_project_atom (const Atom * a,
                            vdw_radius_t vdw)
{
//cerr << "Projecting atom " << a->x() << "," << a->y() << "," << a->z() << ", ny " << _ny << " nz " << _nz << endl;

  if (_ny > _nz)
    _project_atom_y (a, vdw);
  else
    _project_atom_z (a, vdw);

  return;
}

void
Shadow_Area::_project_atom_y (const Atom * a,
                              vdw_radius_t vdw)
{
  double y0 = a->y ();
  double z0 = a->z ();

  double miny = y0 - vdw;
  double maxy = y0 + vdw;

  int ystart = static_cast<int> ((miny - _ymin) / _resolution);
  assert (ystart >= 0);

  int ystop  = static_cast<int> ((maxy - _ymin) / _resolution) + 1;
  if (ystop >= _ny)
    ystop = _ny - 1;

//cerr << "Mol " << _ymin << " Z " << _zmin << " resolution " << _resolution << endl;
//cerr << "Box " << ystart << ',' << ystop << " Z " << zstart << ',' << zstop << endl;

#ifdef DEBUG_SHADOW_AREA
  cerr << "Atom at " << y0 << ',' << z0 << " (vdw " << vdw << ") between Y=" << ystart << " and " << ystop << endl;
#endif

  double vdw2 = vdw * vdw;

  for (int i = ystart; i <= ystop; i++)
  {
    double yg = _ymin + static_cast<coord_t> (i) * _resolution;

    double tmp = vdw2 - (yg - y0) * (yg - y0);
    if (tmp < 0.0)     // doesn't actually intersect
      continue;

    tmp = sqrt (tmp);

    int zstart = static_cast<int> ((z0 - tmp - _zmin) / _resolution);
    if (zstart < 0)
      zstart = 0;

    int zstop  = static_cast<int> ((z0 + tmp - _zmin) / _resolution) + 1;
    if (zstop >= _nz)
      zstop = _nz - 1;

#ifdef DEBUG_SHADOW_AREA
    cerr << "I = " << i << " yg = " << yg << " Z between " << zstart << " and " << zstop << endl;
#endif

    for (int j = zstart; j <= zstop; j++)
    {
//    cerr << i << ',' << j << " at " << zg << " hit " << _grid[i * _ny + j] << endl;

      _grid[i * _nz + j] = 1;
    }
  }

  return;
}

void
Shadow_Area::_project_atom_z (const Atom * a,
                              vdw_radius_t vdw)
{
  double y0 = a->y ();
  double z0 = a->z ();

  double minz = z0 - vdw;
  double maxz = z0 + vdw;

  int zstart = static_cast<int> ((minz - _zmin) / _resolution);
  assert (zstart >= 0);

  int zstop  = static_cast<int> ((maxz - _zmin) / _resolution) + 1;
  if (zstop >= _nz)
    zstop = _nz - 1;

//cerr << "Mol " << _ymin << " Z " << _zmin << " resolution " << _resolution << endl;
//cerr << "Box " << zstart << ',' << zstop << endl;

#ifdef DEBUG_SHADOW_AREA
  cerr << "Atom at " << a->y () << ',' << z0 << " (vdw " << vdw << ") between Z = " << zstart << " and " << zstop << endl;
#endif

  double vdw2 = vdw * vdw;

  for (int i = zstart; i <= zstop; i++)
  {
    double zg = _zmin + static_cast<coord_t> (i) * _resolution;

    double tmp = vdw2 - (zg - z0) * (zg - z0);
    if (tmp < static_cast<coord_t> (0.0))     // doesn't actually intersect
      continue;

    tmp = sqrt (tmp);

    int ystart = static_cast<int> ((y0 - tmp - _ymin) / _resolution);
    if (ystart < 0)
      ystart = 0;

    int ystop  = static_cast<int> ((y0 + tmp - _ymin) / _resolution) + 1;
    if (ystop >= _ny)
      ystop = _ny - 1;

#ifdef DEBUG_SHADOW_AREA
    cerr << "I = " << i << " zg = " << zg << " Y between " << ystart << " and " << ystop << endl;
#endif

    for (int j = ystart; j <= ystop; j++)
    {
//    cerr << i << ',' << j << " at " << zg << " hit " << _grid[i * _ny + j] << endl;

      _grid[i * _ny + j] = 1;
    }
  }

  return;
}

int
Shadow_Area::print_grid (std::ostream & os,
                         char on,
                         char off) const
{
  assert (ok ());

  assert (NULL != _grid);

  os << "Grid " << _ny << " x " << _nz << " total " << (_ny * _nz) << endl;

  if (_ny > _nz)
    return _print_grid_y (os, on, off);
  else
    return _print_grid_z (os, on, off);
}

int
Shadow_Area::_print_grid_y (std::ostream & os,
                            char on,
                            char off) const
{
  int grid_points_hit = 0;

  IWString buffer;
  buffer.resize (_nz);

  for (int i = 0; i < _ny; i++)
  {
    for (int j = 0; j < _nz; j++)
    {
      int k = _grid[i * _nz + j];
      if (k)
      {
        grid_points_hit++;
        buffer += on;
      }
      else
        buffer += off;
    }

    os << buffer << endl;
    buffer.resize_keep_storage (0);
  }

  os << grid_points_hit << " grid points hit\n";

  return os.good ();
}

int
Shadow_Area::_print_grid_z (std::ostream & os,
                            char on,
                            char off) const
{
  int grid_points_hit = 0;

  IWString buffer;
  buffer.resize (_ny);

#ifdef PRINT_LANDSCAPE
  os << "Landscape\n";
  for (int i = 0; i < _ny; i++)
  {
    for (int j = 0; j < _nz; j++)
    {
      int k = _grid[i + j * _ny];
      if (k)
        buffer += on;
      else
        buffer += off;
    }

    os << buffer << endl;
    buffer.resize_keep_storage (0);
  }

  os << "Portrait\n";
#endif

  for (int i = 0; i < _nz; i++)
  {
    for (int j = 0; j < _ny; j++)
    {
      int k = _grid[i * _ny + j];
      if (k)
      {
        grid_points_hit++;
        buffer += on;
      }
      else
        buffer += off;
    }

    os << buffer << endl;

    buffer.resize_keep_storage (0);
  }

  os << grid_points_hit << " grid points hit\n";

  return os.good ();
}
