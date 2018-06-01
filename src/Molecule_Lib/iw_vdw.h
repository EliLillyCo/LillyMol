#ifndef IWVDW_H
#define IWVDW_H

#define IW_VDW_SHRAKE_AND_RUPLEY 1
#define IW_VDW_SAVOL 2
#define IW_VDW_MOLVOL 3
#define IW_VDW_SYBYL63 4
#define IW_VDW_WIKI 5

typedef double vdw_radius_t;

class Molecule;

extern int assign_vdw_radii (Molecule &, int, vdw_radius_t *);

class Command_Line;

extern int set_default_van_der_waals_radius_type (Command_Line &, char, int &, int = 0);
extern int display_standard_vdw_radius_types (std::ostream &, char, int = 0);

#endif
