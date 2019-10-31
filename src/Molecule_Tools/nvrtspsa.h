#ifndef NVRTSPSA_H
#define NVRTSPSA_H

extern double novartis_polar_surface_area (Molecule & m);

extern double novartis_polar_surface_area (Molecule & m, 
                             const atomic_number_t * z,
                             const Atom ** atom,
                             const int * is_aromatic);
extern void   set_display_psa_unclassified_atom_mesages(int s);
extern void   set_return_zero_for_unclassified_atoms (int s);
extern void   set_non_zero_constribution_for_SD2(int s);
#endif
