    int _write_molecule_mrv (std::ostream & os) const;
    int _write_atoms_mrv (std::ostream & os) const;
    int _write_bonds_mrv (std::ostream & os) const;

    int read_molecule_mrv_mchemical (XMLNode & cml);
    int _read_atom_array_mrv (XMLNode & cml);
    int _read_bond_array_mrv (XMLNode & cml, int *);
    int _read_atom_array_mrv_individual_attributes (const XMLNode & xml);
