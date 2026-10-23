#ifndef E_STATE_COMPUTATION_H
#define E_STATE_COMPUTATION_H

#include "Molecule_Lib/molecule.h"

namespace lillymol_estate {
int determine_atom_e_state_index (Molecule &m, double [], const Atom * const *, const atomic_number_t *, const int * hcount);
int determine_atom_e_state_index (Molecule &, double []);
int determine_hydrogen_e_state_index (Molecule &, double []);
int determine_atom_e_state_index (Molecule &m, double e_state_index [], const double * i_state, const atomic_number_t * z);
int number_of_bond_in_conjugated_system (Molecule &, const int);
int is_triple_bonded (Molecule &, const int);
double value_of_Kier_and_Hall_atom_intrinsic_state (Molecule &, const int, const Atom * const *, const atomic_number_t *, const int *);
}  // namespace lillymol_estate

#endif
