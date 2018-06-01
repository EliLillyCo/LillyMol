//  Temporary stuff for cis-trans bonds within the molecule class

  private:
//  Cis_Trans_Bond * _create_coupled_cis_trans_bond (atom_number_t previous_atom, int);
//  Cis_Trans_Bond * _create_new_cis_trans_bond (atom_number_t previous_atom, int);

//  int _handle_right_side_of_existing_cis_trans_bond (const char *smiles,
//                                                     Cis_Trans_Bond *ctb);
//  int _handle_left_side_of_existing_cis_trans_bond (const char *smiles,
//                                                     Cis_Trans_Bond *ctb);
//  int _handle_part_of_existing_cis_trans_bond (atom_number_t previous_atom,
//                                               const char *smiles,
//                                               Cis_Trans_Bond *ctb);

//  int _check_compatibility_with_existing_directional_bonds (Cis_Trans_Bond * ctb);

//  int _cis_trans_bonds_consistent (const Cis_Trans_Bond * ctbi,
 //                                  const Cis_Trans_Bond *ctbj) const;
    int _check_all_cis_trans_bonds () const;

//  int _assign_directional_bonds (const Cis_Trans_Bond * ctb);
//  int _assign_directional_bonds ();

    int _process_directional_bond_for_smiles (IWString & smiles,
                                              const Bond *,
                                              atom_number_t anchor);

    int _adjacent_directional_bonds_ok (const Bond &) const;
    int _adjacent_directional_bonds_ok (atom_number_t) const;

    int _adjacent_directional_bonds_mutually_consistent (atom_number_t zatom);
    int _can_be_end_of_directional_bond (atom_number_t zatom) const;
    int _mark_adjacent_double_bond_with_directional_atoms_at_other_end (atom_number_t a3, int *);

    int _unset_directional_bonds_not_adjacent_to_a_double_bond();;

//  int _got_other_end_of_directional_bond (atom_number_t, int);
    int _check_for_cis_trans_bond (atom_number_t previous_atom);
//  Cis_Trans_Bond * _cis_trans_bond_anchored_at (atom_number_t a) const;

//  int _next_available_cis_trans_coupled_identifier () const;

    int _smiles_write_directional_bond (atom_number_t a,
                     atom_number_t anchor, IWString & smiles);
//  int _discern_cis_trans_couplings (Cis_Trans_Bond * ctb);
//  int _set_directional_bond (Cis_Trans_Bond *, Bond *);

//  int _adjust_cis_trans_for_loss_of_atom (atom_number_t);

//  int _transfer_cis_trans_bonds (const Molecule & m2);

    int _discern_cis_trans_bond_from_depiction (Bond * b);
    int _extend_cis_trans_system (atom_number_t zatom);

    int _discern_cis_trans_bond_from_depiction (atom_number_t zatom);
    int _invalidate_directional_bonds_involving (atom_number_t zatom);
    int _atom_being_unbonded_check_directional_bonds (atom_number_t, int preserve_chirality = 0);
    int _finished_reading_smiles_assign_and_check_directional_bonds ();
    int __finished_reading_smiles_assign_and_check_directional_bonds ();
    int _invalidate_directional_double_bond (Bond &);
    int _invalidate_directional_bonds_at_end_of_double_bond (atom_number_t zatom);

//                                                resizable_array<const Bond *> & already_processed);
//  int _discern_cis_trans_bond_couplings_from_directional_bonds ();
//  int _discern_cis_trans_bond_couplings_from_directional_bonds (atom_number_t astart,
//                                              int * at_end_of_directional_bond);
    int _find_atoms_attached_to_root (atom_number_t root,
                                        atom_number_t & down_atom,
                                        atom_number_t & up_atom,
                                        atom_number_t & doubly_bonded) const;
    int _append_bad_cis_trans_input_text_to_name ();

    int _cis_trans_bond_has_been_invalidated (atom_number_t zatom);

    int _adjust_cis_trans_bonds_to_canonical_form(const int *);
    const Bond * _identify_double_bond (atom_number_t zatom) const;
    int _identify_linked_cis_trans_bonds(resizable_array<const Bond *> & bonds_to_be_flipped,
                                           const Bond * current_bond,
                                           int * bond_already_done) const;
    int _identify_linked_cis_trans_bonds(resizable_array<const Bond *> & bonds_to_be_flipped,
                                         atom_number_t previous_atom,
                                         atom_number_t zatom,
                                         int * bond_already_done) const;
    int _canonicalise_linked_group_of_cis_trans_bonds(const resizable_array<const Bond *> &bonds_in_grouping,
                                        const int * canonical_rank);

    int _bond_is_no_longer_directional (const Bond * b);
    int _remove_directional_bonding_associated_with_bond (const Bond * b);
    int _process_directional_system (atom_number_t lhs1,
                                       atom_number_t db1,
                                       int & valid_cis_trans_form_found);
    int _identify_directional_bonds_across_double_bonds (const Bond * b,
                                resizable_array<const Bond *> & coupled) const;
    int _identify_directional_bonds_across_double_bonds (atom_number_t zatom,
                                                         resizable_array<const Bond *> & coupled) const;

    int _set_any_unambiguous_unset_directional_bonds(resizable_array<const Bond *> & directional_bonds);
    int _get_single_bonds (atom_number_t zatom,
                             resizable_array<const Bond *> & sb) const;
    int _fill_in_missing_directional_bond_specification (atom_number_t zatom,
                                        resizable_array<const Bond *> & directional_bonds);
  public:
//  const Cis_Trans_Bond * cis_trans_bond_anchored_at (atom_number_t a) const;
//  const Cis_Trans_Bond * part_of_cis_trans_bond (atom_number_t, atom_number_t) const;
//  int extra_cis_trans_bond (Cis_Trans_Bond *, int = 1);

