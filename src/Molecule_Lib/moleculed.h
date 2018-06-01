#ifndef COMPILING_MOLECULED
  THIS FILE SHOULD ONLY BE INCLUDED IN MOLECULED.CC
#else

    void _compute_distance_matrix ();
    int _initialise_distance_matrix ();
    int _bonds_between (atom_number_t, atom_number_t);
    int _recompute_distance_matrix (int (Molecule::*identify_first_atom) (const int *, atom_number_t &),
               int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &));
    void _compute_row_of_distance_matrix (int * row_of_distance_matrix,
                                 atom_number_t current_atom,
                                 int distance);
    void _compute_row_of_distance_matrix (CRDM_args & crdm,
               int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &));
    void _compute_row_of_distance_matrix (int * row_of_distance_matrix,
                                 int & distance,
                                 int * atom_stack,
                                 int stack_ptr,
                                 int * ring_atom);

    int _atoms_between (atom_number_t a1,
                          atom_number_t a2,
                          int d,
                          Set_of_Atoms & s);

#endif
