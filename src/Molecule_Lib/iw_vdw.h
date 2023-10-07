#ifndef MOLECULE_LIB_IWVDW_H_
#define MOLECULE_LIB_IWVDW_H_

#include <optional>

#define IW_VDW_SHRAKE_AND_RUPLEY 1
#define IW_VDW_SAVOL 2
#define IW_VDW_MOLVOL 3
#define IW_VDW_SYBYL63 4
#define IW_VDW_WIKI 5

#include "Molecule_Lib/molecule.h"

typedef double vdw_radius_t;


extern int assign_vdw_radii(Molecule &, int, vdw_radius_t *);

class Command_Line;

extern int set_default_van_der_waals_radius_type(Command_Line &, char, int &, int = 0);
extern int display_standard_vdw_radius_types(std::ostream &, char, int = 0);

namespace vdw {
  // Return the radius for atom `zatom` given `vdw_type`.
  std::optional<vdw_radius_t> vdw_radius(Molecule& m, atom_number_t zatom, vdw::VdwType vdw_type);
  // Assign all vdw radii in
  int AssignVdwRadii(Molecule& m, vdw::VdwType vdw_type, vdw_radius_t* vdw);
}  // namespace vdw


#endif  // MOLECULE_LIB_IWVDW_H_
