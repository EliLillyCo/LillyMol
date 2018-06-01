    int _read_molecule_mol2_ds (iwstring_data_source & input);
    int _read_molecule_mol2_ds (iwstring_data_source & input, int na, int nb, int * aromatic_atom, int * aromatic_bonds, Tripos_Residue_Information * tri);

    int _parse_tripos_atom_record (const const_IWSubstring & buffer, atom_number_t, int &, Tripos_Residue_Information * tri);
    int _parse_tripos_bond_record (const const_IWSubstring & buffer, int * aromatic_atoms, int, int * aromatic_bond);
    int _mol2_assign_default_formal_charges ();
    int _doubly_bonded_to_oxygen (atom_number_t zatom) const;
    int _tripos_atom_type_from_string (atom_number_t, const const_IWSubstring &);
    int _place_formal_charges_on_quat_n_from_mol2 ();
