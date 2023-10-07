#ifndef SAF_MOLVOL_H
#define SAF_MOLVOL_H

#include "Foundational/iwaray/iwaray.h"

#include "Molecule_Lib/iw_vdw.h"
#include "Molecule_Lib/iwmtypes.h"

class Surface_Area_Molvol
{
  private:
    vdw_radius_t _probe_radius;
    int _random_number_counter;

//  What type of Van der Waals radii to use

    int _vdw_radius_type;

//  We can impose a maximum number of times a structure can be wobbled

    int _max_wobble;

    double _wobble_amplitude;

//  I've found that things work better if you wobble first
    
    int _wobble_first;

//  For our report, we keep track of how much wobbling we do

    int _molecules_processed;

    int _molecules_wobbled;

    extending_resizable_array<int> _wobbled;

//  It can be desirable to return the molecule unperturbed

    int _save_coordinates;

//  private functions

    void _wobble_structure (Molecule & m);

    int _molvol (Molecule & m,
                 area_t * area,
                 area_t & total_area,
                 volume_t & volume);

  public:
    Surface_Area_Molvol ();

    void set_probe_radius (vdw_radius_t r) { _probe_radius = r;}

    void set_vdw_type (int v) { _vdw_radius_type = v;}

    void set_max_wobble (int m) { _max_wobble = m;}

    void set_wobble_first (int w) { _wobble_first = w;}

    int report (std::ostream &) const;

    int surface_area (Molecule &, area_t *, area_t &, volume_t &);

    //    int reset_random_number_counter () {_random_number_counter = 0;}
};

#endif
