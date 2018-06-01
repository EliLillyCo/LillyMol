#ifndef MOLECULE_MAIN_H
#define MOLECULE_MAIN_H

  private:
    int _convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s, int * barray) const;

    int  _ok_ring_info () const;
    int  _invalidate_ring_info ();
    int  _invalidate_ring_aromaticity_info ();

    void _compute_element_count (int * element_count, int & highest_atomic_number, int & isotopes_present, int & non_periodic_table_elements_present) const;
    void _compute_element_count (int * element_count, const int * include_atom, int & highest_atomic_number, int & isotopes_present, int & non_periodic_table_elements_present) const;
    void _compute_element_count (int * element_count, const int * atom_flag, int flag, int & highest_atomic_number, int & isotopes_present, int & non_periodic_table_elements_present) const;

    int  _remove_atom (atom_number_t);

    void remove_atom_from_charge_arrays (const atom_number_t atom_to_remove);

    int _invalidate_for_changed_isotope ();
    int _exact_mass (const int * element_count, int highest_atomic_number,
                     int non_periodic_table_elements_present,
                     exact_mass_t & result) const;

    int _set_bond_length (atom_number_t a1, atom_number_t a2,
                            distance_t d, int * either_side);

    int _set_isotope_zero(atom_number_t zatom);

    int _all_connections_saturated (const atom_number_t zatom, const atom_number_t ignore) const;
    int _double_bond_needs_changing_for_graph_form (const Bond & b, const Mol2Graph & mol2graph) const;

#endif
