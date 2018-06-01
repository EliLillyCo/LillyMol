#ifndef IS_ACTUALLY_CHIRAL_H
#define IS_ACTUALLY_CHIRAL_H

#include "iwmtypes.h"
#include "path_scoring.h"

class Molecule;

extern int is_actually_chiral (Molecule & m, atom_number_t);
extern int is_actually_chiral (Molecule & m, atom_number_t, resizable_array_p<Path_Scoring> &);

// Returns the number of invalid chiral centres removed

extern int do_remove_invalid_chiral_centres (Molecule & m);

extern void set_allow_unsaturated_atoms_to_be_chiral (int s);

extern void set_max_iterations (int m);

#endif
