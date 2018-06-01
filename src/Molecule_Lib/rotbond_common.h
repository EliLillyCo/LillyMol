#ifndef ROTBOND_COMMON_H
#define ROTBOND_COMMON_H

extern int part_of_otherwise_non_rotabable_entity (Molecule & m,
                                        atom_number_t a1,
                                        atom_number_t a2);

extern int triple_bond_at_either_end (const Molecule & m,
                           const Bond * b);

extern int is_non_rotatable_amide (Molecule & m,
                        atom_number_t a1,
                        atom_number_t a2);

extern int is_non_rotatable_sulphonamide (Molecule & m,
                               atom_number_t zatom1,
                               atom_number_t zatom2);
#endif

/* arch-tag: 28fd8c6d-c5cb-450f-9aed-6a78cdaa3400

*/
