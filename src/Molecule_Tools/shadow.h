#ifndef IWSHADOW_H
#define IWSHADOW_H

/*
  Several places in the code I need to compute shadow areas
*/

#include "molecule.h"
#include "iw_vdw.h"

class Shadow_Area 
{
  private:
    double _resolution;

    int * _grid;

//  We need to keep track of how many divisions along each axis

    int _ny, _nz;

//  The extremeties of the molecule on the Y/Z plane

    double _ymin, _zmin;

//  private functions

    area_t _shadow_area (int matoms, const Atom * const * atom, const vdw_radius_t vdw []);
    void _project_atoms_to_grid (int matoms,
                                 const Atom * const * atom,
                                 const vdw_radius_t * vdw);

    void _project_atom (const Atom * a, vdw_radius_t vdw);
    void _project_atom_y (const Atom * a, vdw_radius_t vdw);
    void _project_atom_z (const Atom * a, vdw_radius_t vdw);

    area_t _determine_area () const;

    int _print_grid_y (std::ostream &, char, char) const;
    int _print_grid_z (std::ostream &, char, char) const;

    int _nearest_y_grid_point (coord_t) const;
    int _nearest_z_grid_point (coord_t) const;

  public:
    Shadow_Area ();
    ~Shadow_Area ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int active () const;

    void set_resolution (float r) { _resolution = r;}

    area_t shadow_area (const Molecule &, const vdw_radius_t *);

    int print_grid (std::ostream &, char, char) const;
};

#endif
