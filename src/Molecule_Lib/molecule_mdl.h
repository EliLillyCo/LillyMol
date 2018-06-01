#ifndef COMPILING_MDL_CC
  THIS FILE SHOULD ONLY BE INCLUDED WHEN COMPILING MDL.CC
#endif

  private:

    template <typename T> int  _read_molecule_mdl_ds (T &, int, MDL_File_Supporting_Material &);

    template <typename T> int  _read_mdl_atom_connection_table (T &, int &, int &, MDL_File_Supporting_Material & mdlfos);
    template <typename T> int  _read_mdl_bond_list (T &, int, int *, int *, MDL_File_Supporting_Material & mdlfos);
    template <typename T> int  _read_molecule_mdl_trailing_records (T &, int, MDL_File_Supporting_Material & mdlfos);

    int _write_M_RGP_records (const MDL_File_Supporting_Material &, std::ostream &) const;
  protected:
    template <typename T> int  _write_m_chg_records (T & os, int) const;
    template <typename T> int  _write_m_iso_records (T & os, int) const;
    int _common_parse_M_record (const const_IWSubstring & buffer,
                                int & fatal);

    int _mdl_atom_is_chiral_centre (atom_number_t zatom, int cfg, MDL_File_Supporting_Material &);

  private:
    int _contains_isotope_above (int) const;


    int _has_delocalised_neighbours (atom_number_t zatom, const int * aromatic_atoms, const int * aromatic_bonds,
                                     Set_of_Atoms & s, const int * ab) const;
    int _unset_aromatic_bond_entry(atom_number_t a1, atom_number_t a2, int * aromatic_bonds) const;
    int _unset_aromatic_bond_entry(atom_number_t a1, atom_number_t a2, int * aromatic_bonds, const int * ab) const;

    int _identify_atoms_at_ends_of_aromatic_bonds(const int *, int *, int *) const;

    int _read_molecule_rdf_ds (iwstring_data_source & input, IWString & possible_name, MDL_File_Supporting_Material & mdlfos);

//  May 98, stuff for Version 3 sd files

    int _read_mdl_V3 (iwstring_data_source &);
    int _write_mdl_V3 (iwstring_data_source &);
    template <typename T> int _read_mdl_atom_connection_table_v30 (T & input, int & nb, MDL_File_Supporting_Material & mdlfos);
    template <typename T> int _process_v30_composite_records(T & input, const MDL_File_Supporting_Material & mdlfsm);
    int _convert_sgroups_to_elements(const resizable_array_p<IWString> & sgroup);
    int _convert_sgroup_to_elements(const IWString & s);
  protected:
    int _parse_v30_atom_record (const IWString & buffer, int, MDL_File_Supporting_Material & mdlfos);
  private:
    template <typename T> int _read_v30_bond_list (T & input, int nb, int *, int *);

    template <typename T> int _write_molecule_atom_list_v30 (T & os, MDL_File_Supporting_Material & mdlfos) const;
    template <typename T> int _write_molecule_bond_list_v30 (T & os) const;

//  Sept 2000, stuff for discerning chirality from wedge bonds

    int _discern_chirality_from_wedge_bond (atom_number_t a1, atom_number_t a2,
                          int direction);
    int _discern_chirality_from_wedge_bond_4 (atom_number_t zatom, atom_number_t a2,
                          int direction);
    int _create_chiral_centre (atom_number_t zatom,
                                 atom_number_t a1,
                                 atom_number_t a2,
                                 atom_number_t a3,
                                 atom_number_t a4,
                                 int direction);
    int _create_unspecified_chirality_object (atom_number_t zatom);

    int _mdl_atom_stereo_value (atom_number_t a, const MDL_File_Supporting_Material & mdlfos) const;

    int _multiple_wedge_bonds_to_atom (atom_number_t a) const;
    int _discern_chirality_from_multiple_wedge_bonds (int bstart);

    int _parse_M_RGP_record (const const_IWSubstring & buffer);

    template <typename T> int _read_mdl_data_following_tag (T & input, const MDL_File_Supporting_Material & mdlfos);

    int _set_elements_based_on_atom_aliases (const resizable_array_p<Atom_Alias> & a);


    int _count_elements (const atomic_number_t z, const Set_of_Atoms & e) const;
  public:
    int  mdl_add_m_formal_charge  (int, const Aprop * atom_properties);
    int  mdl_add_m_radical (int, const Aprop * atom_properties);
    int  mdl_add_m_isotope (int, const Aprop * atom_properties);

//  Shared with the function that writes ISIS reaction files

    int _write_mdl_atom_record_element_charge_and_chirality (atom_number_t, IWString & output_buffer, MDL_File_Supporting_Material & mdlfos) const;

    template <typename T> int _mdl_write_atoms_and_bonds_record (T &, int nfc, int iat, MDL_File_Supporting_Material & mdlfos) const;

    int _process_mdl_g_record (const IWString &, const const_IWSubstring & buffer);
    int _assign_strange_atomic_symbol_to_atom (atom_number_t zatom, const_IWSubstring s);

  private:
