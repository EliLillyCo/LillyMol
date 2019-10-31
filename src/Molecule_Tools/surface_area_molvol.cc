/*
  interface to molvol
*/

#include <stdlib.h>

#include "iwrandom.h"
#include "molecule.h"
#include "iw_vdw.h"

#include "surface_area_molvol.h"
#include "jwrandom_preset_array.h"

extern "C" void volume_ (int &);

// n_Max must be the same as in volume.f

#define n_Max 301

struct RawData
{
  double x_Chain[n_Max];
  double y_Chain[n_Max];
  double z_Chain[n_Max];
  double Radii[n_Max];
  double RefRadius;
  int    n_Bonds;
};

extern struct RawData rawdata_;


/*
  We get our results from this common block
*/

struct MyResults
{
  double atomic_volume[n_Max];
  double atomic_area[n_Max];
};

extern struct MyResults volumeoutput_;

Surface_Area_Molvol::Surface_Area_Molvol ()
{
  _vdw_radius_type = IW_VDW_MOLVOL;

  _random_number_counter = 0;

  _probe_radius = 1.40;

  _molecules_processed = 0;
  _molecules_wobbled = 0;

  _max_wobble = 20;

  _wobble_amplitude = 0.05;

  _wobble_first = 1;

  _wobbled.resize (_max_wobble + 2);

  _save_coordinates = 1;

  return;
};

int
Surface_Area_Molvol::report (std::ostream & os) const
{
  os << "Surface_Area_Molvol::report: " <<  _molecules_processed << " molecules processed, " << _molecules_wobbled << " wobbled\n";

  if (_wobble_first)
    os << "All molecules wobbled before computation\n";

  for (int i = 0; i < _wobbled.number_elements (); i++)
  {
    int w = _wobbled.item (i);      // operator [] is non-const

    if (w)
        os << w << " molecules wobbled " << i << " times\n";
  }

  return os.good ();
}

/*
  We want to ensure that when we perturb the molecule, we actually move it.
  This function generates an offset that is somewhere between 0.1 and 1.0 of
  amplitude, or the same as a negative
*/

static coord_t
delta (coord_t amplitude, int counter)
{
  //  random_number r = iwrandom ();

  random_number_t r = random_number_preset_array [counter];
  //  cout<<r<<" ";

  if (r < 0.0)
  {
    return -( amplitude * 0.2 - r * 0.8 * amplitude );
  }
  else
  {
    return  ( amplitude * 0.2 + r * 0.8 * amplitude );
  }
}

void
Surface_Area_Molvol::_wobble_structure (Molecule & m)
{
//cerr << "Wobbling structure, amplitude " << wobble_amplitude << endl;
  int matoms = m.natoms ();
  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi (i);

    coord_t x = a->x () + delta (_wobble_amplitude, _random_number_counter*3);
    coord_t y = a->y () + delta (_wobble_amplitude, _random_number_counter*3 +1);
    coord_t z = a->z () + delta (_wobble_amplitude, _random_number_counter*3 +2);

    _random_number_counter++;
    if (_random_number_counter>=MAX_RANDOM_NUMBER_TRIPLE)
      _random_number_counter = 0;

    m.setxyz (i, x, y, z);
  }

  return;
}

extern int verbose;

static void
copy_coordinates_to_common_blocks (const Molecule & m, double RMax_Radius)
{
  int matoms = m.natoms ();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi (i);

    rawdata_.x_Chain[i] = a->x () * RMax_Radius;
    rawdata_.y_Chain[i] = a->y () * RMax_Radius;
    rawdata_.z_Chain[i] = a->z () * RMax_Radius;
  }

  return;
}

int
Surface_Area_Molvol::_molvol (Molecule & m,
                              area_t * area,
                              area_t & total_area,
                              volume_t & volume)
{
  int matoms = m.natoms ();

// molvol uses n_bonds where we would think of atoms. Not sure why...
// also, the n_Spheres variable in MOLVOL seems to be unused...

  rawdata_.RefRadius = 0.0;
  rawdata_.n_Bonds = matoms;
  //  rawdata_.n_Bondsm1 = matoms - 1;
  //  rawdata_.n_Spheres = 0;

  double MaxRadius = 0.0;
  for (int i = 0; i < matoms; i++)
  {
    rawdata_.Radii[i] += _probe_radius;

    if (rawdata_.Radii[i] > MaxRadius)
      MaxRadius = rawdata_.Radii[i];
  }

// Rescale radii with the maximum radius

  double RMax_Radius = 1.0 / MaxRadius;
  for (int i = 0; i < matoms; i++)
  {
    rawdata_.Radii[i]   *= RMax_Radius;
  }

  rawdata_.RefRadius = MaxRadius;

  copy_coordinates_to_common_blocks (m, RMax_Radius);

  int wobbles_this_molecule = 0;
  while (wobbles_this_molecule <= _max_wobble)
  {
    int error_encountered=0;
    volume_ (error_encountered);

    if (! error_encountered)     // great, we are done
      break;

    wobbles_this_molecule++;

    if (verbose)
      cerr << "Molecule needs to be wobbled " << wobbles_this_molecule << endl;

    _wobble_structure (m);
    copy_coordinates_to_common_blocks(m, RMax_Radius);
  }

  if (wobbles_this_molecule > _max_wobble)
    return 0;

  _wobbled[wobbles_this_molecule]++;

  volume = 0.0;
  total_area = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    volume += volumeoutput_.atomic_volume[i];
    area[i] = volumeoutput_.atomic_area[i];
    total_area   += area[i];
  }

  return 1;
}

int
Surface_Area_Molvol::surface_area (Molecule & m,
                                   area_t * area,
                                   area_t & total_area,
                                   volume_t & volume)
{
  int matoms = m.natoms ();

  if (matoms > n_Max)
  {
    cerr << "Surface_Area_Molvol::compute_surface_area: too many atoms " << matoms << " in '" << m.name () << "'\n";
    return 0;
  }

  if (! assign_vdw_radii (m, _vdw_radius_type, rawdata_.Radii))
  {
    cerr << "Cannot assign Van der Walls radii, type " << _vdw_radius_type << endl;
    return 0;
  }

  _molecules_processed++;

  Coordinates * coords;

  if (_save_coordinates)
  {
    coords = new Coordinates[matoms];
    m.get_coords (coords);
  }
  else
    coords = NULL;

  if (_wobble_first)
    _wobble_structure (m);    // found things work better if we do a first wobble

  int rc = _molvol (m, area, total_area, volume);

  if (NULL != coords)
  {
    m.setxyz (coords);
    delete [] coords;
  }

  return rc;
}
