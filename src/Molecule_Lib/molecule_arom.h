//  Include file for all the private molecule functions in aromatic.cc

#ifndef COMPILING_AROMATIC_CC
    YIPES, THIS SHOULD ONLY BE INCLUDED FROM AROMATIC.CC
#endif

//  Stuff added jan 2004 for dealing with purely aromatic structures

    int _a_atoms_found_this_set_of_rings (const Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings);
    int _kekule_all_atoms_and_bonds_in_system_aromatic (Kekule_Temporary_Arrays & kta,
                                const resizable_array<const Ring *> & rings);

    int _allocate_aromaticity ();

    int _count_pi_electrons_in_any_fused_but_kekule_form_rings (Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings);
    int _sort_possible_continuations_by_ring_size (atom_number_t, Set_of_Atoms & p);
    int _bond_only_in_small_carbocycle (atom_number_t a0, atom_number_t a1);

    int __compute_aromaticity_for_ring (const Ring &, aromaticity_type_t &,
                                        int &, int *, int &, int);
    int _compute_aromaticity_for_ring (int, aromaticity_type_t &,
                                       int &, int *, int &, int);

    int _compute_aromaticity_with_fused_rings (int *, const int *, int *, int);

    int _compute_aromaticity (int *, int *, int *, int *, int);

//  Update aromaticity values for a ring number, or for a group of atoms (most probably
//  the group is part of a ring.

    int _update_aromaticity (int, aromaticity_type_t);
    int _update_aromaticity (const Ring &, aromaticity_type_t);

//  These two functions do practically the same thing

    int _make_bonds_aromatic (const Ring *);
    int _make_all_bonds_aromatic (const Ring & r);
    int _make_all_unshared_bonds_non_aromatic (const Ring & r);

    int _count_aromatic_bonds_in_just_one_ring (const Ring & r);

    int _determine_aromaticity_by_single_fusions_to_ring (int * pi_electron_count,
                                                    int aromaticity_rules,
                                                    int * rings_in_fused_system,
                                                    int * ok_to_include,
                                                    const int * impossible_aromatic,
                                                    const Ring * r);
    int _determine_aromaticity_by_single_fusions (int * ring_already_done, 
                                                    int * pi_electron_count,
                                                    const int * impossible_aromatic,
                                                    int aromaticity_rules);

    int _assemble_super_ring (Set_of_Atoms & super_ring,
                                int * to_be_processed,
                                resizable_array<int> & rings_used);
    int _whole_fused_system_is_aromatic (int * already_done,
                                         int * pi_electron_count,
                                         int aromaticity_rules,
                                         int * rings_in_fused_system);
    int _whole_fused_system_is_aromatic (int * already_done,
                                         int fused_system_id,
                                         int * pi_electron_count,
                                         const int * impossible_aromatic,
                                         int aromaticity_rules);

    int _unshared_pi_electrons (const Ring * r, const int * pi_electrons);

    int _combine_non_arom_ring (int zring,
                                  int system_identifier,
                                  int * already_done,
                                  int * pi_electron_count,
                                  const int * impossible_aromaticic,
                                  const int aromaticity_rules,
                                  int * rings_in_ring_system);

    int _determine_aromaticity_of_fused_systems (int * already_done,
                                                 int * pi_electron_count,
                                                 const int * impossible_aromatic,
                                                 const int aromaticity_rules,
                                                 int & fused_system_identifier);
    int _determine_aromaticity_of_fused_systems (int * already_done,
                                                 int * pi_electron_count,
                                                 const int * impossible_aromatic,
                                                 const int aromaticity_rules,
                                                 int *);
    int _determine_aromaticity (const Set_of_Atoms &,
                                aromaticity_type_t &,
                                int &,
                                int *, 
                                int);
    int _determine_aromaticity (const int * process_these,
                                aromaticity_type_t & result,
                                int & impossible_aromatic,
                                int * pi_electron_count,
                                int);

//  All the private functions needed for kekule determinations

    int _count_pi_electrons (atom_number_t zatom,
                              const int * process_these,
                              int & result);

    int _do_obvious_hcount_adjustments (Kekule_Temporary_Arrays &);
    int _do_obvious_bond_order_settings (Kekule_Temporary_Arrays &);

    int _kekule_find_ring_system (Kekule_Temporary_Arrays &,
                                  resizable_array<const Ring *> & rings);

    int __kekule_arom_test_rings (const int * process_these,
                                   const resizable_array<const Ring *> & rings,
                                   aromaticity_type_t * ring_aromaticity,
                                   aromaticity_type_t * atom_aromaticity,
                                   int * pi_electrons,
                                   int * unshared_pi_electrons);
    int _kekule_arom_test_rings (Kekule_Temporary_Arrays & kta,
                                 const resizable_array<const Ring *> & rings);

    int _kekule_identify_next_atom (const int * process_these,
                                    atom_number_t & a) const;
    int _doubly_bonded_to_something_outside_ring(atom_number_t zatom) const;
    int _bonded_to_heteroatom_outside_system (atom_number_t a, const int * aromatic_system);

    int __find_kekule_form_current_config (const resizable_array<const Ring *> & rings,
                                Kekule_Temporary_Arrays & kta,
                               atom_number_t a,
                               int pi_electron_count);
    int _find_kekule_form_current_config (const resizable_array<const Ring *> & rings,
                                Kekule_Temporary_Arrays & kta,
                               atom_number_t a,
                               int pi_electron_count);
    int _find_kekule_form (const resizable_array<const Ring *> & rings,
                           Kekule_Temporary_Arrays &,
                           atom_number_t a,
                           int pi_electron_count);
    int _find_kekule_form_ring_system (Kekule_Temporary_Arrays &,
                           resizable_array<const Ring *> & rings);

    int _find_kekule_form (Kekule_Temporary_Arrays &);


    int _kekule_suppress_non_aromatic_rings (Kekule_Temporary_Arrays &);
    int _kekule_check_rings_containing_aliphatic_bonds (Kekule_Temporary_Arrays & kta);

    int _kekule_identify_non_arom_rings (Kekule_Temporary_Arrays &);

    int _kekule_cannot_be_aromatic (int *, int &, int &);

    int _convert_chain_aromatic_bonds (int * aromatic_atoms, const int * rm);
    int _convert_any_chain_aromatic_atoms_to_permanent_aromatic (const int * process_these_atoms);
    int _maybe_all_atoms_are_permanent_aromatic (Kekule_Temporary_Arrays &, const resizable_array<const Ring *> & rings);
    int _convert_atom_to_permanent_aromatic (atom_number_t zatom);

    int _explicit_hydrogen_bonded_to_aromatic_atom (const int * aromatic_atoms) const;

    int _find_kekule_form_positive_nitrogen (Kekule_Temporary_Arrays &);
    int _find_kekule_form_positive_nitrogen (Set_of_Atoms & possible_nplus, Kekule_Temporary_Arrays &);

//  int _find_kekule_form_positive_nitrogen_single_ring (atom_number_t n,
//                           const Ring * rn,
//                           Kekule_Temporary_Arrays &);
    int _find_kekule_form_positive_nitrogen (int * possible_nplus,
                             atom_number_t n,
                             Kekule_Temporary_Arrays &);

    void _grow_fused_system (const Ring * rn,
                             Kekule_Temporary_Arrays &,
                              const int * possible_nplus);

    const Ring * _single_active_ring_containing_atom (atom_number_t n, const int * process_these_rings);

    int _has_aromatic_bonds (const Atom * a) const;

    int _find_kekule_form_positive_nitrogen_ring_system (Kekule_Temporary_Arrays &,
                                        const resizable_array<const Ring *> &,
                                        const Set_of_Atoms &);
    int _find_kekule_form_positive_nitrogen_ring_system (Kekule_Temporary_Arrays & kta,
                                                        const resizable_array<const Ring *> & rings,
                                                        const Set_of_Atoms & possible_nplus_in_system,
                                                        Set_of_Atoms & s,
                                                        int maxdepth);
    int _find_kekule_form_positive_nitrogen_ring_system_place_charges (Kekule_Temporary_Arrays & kta,
                                                        const resizable_array<const Ring *> & rings,
                                                        const Set_of_Atoms & possible_nplus_in_system);

    int _all_bonds_aromatic (const Ring & r) const;

    int _only_process_rings_with_all_aromatic_atoms (Kekule_Temporary_Arrays & kta);

    int _determine_hydrogen_counts (Kekule_Temporary_Arrays & kta);
    int _bonds_to_atom_cannot_vary (Kekule_Temporary_Arrays & kta, atom_number_t a);

    int _group_rings_into_groups_to_be_processed_together (Kekule_Temporary_Arrays & kta);
    int _grow_fused_system (const Ring & r, int flag, int * in_system);

    int _find_next_ring_or_ring_system (Kekule_Temporary_Arrays & kta,
                                         int flag,
                                         resizable_array<const Ring *> & rings);

    int _identify_possible_nplus_in_system (Set_of_Atoms & possible_nplus,
                                            Set_of_Atoms & possible_nplus_in_system,
                                            const int * process_these_atoms) const;
    int _identify_kekule_search_starting_atom (Kekule_Temporary_Arrays & kta,
                                                 atom_number_t & astart) const;
    void _set_all_bonds_in_system_to_single (const int * process_these_atoms);
    int _restore_existing_bond_types (const bond_type_t * bsave);

    int _atom_is_adjacent_to_aromatic_carbonyl (Kekule_Temporary_Arrays & kta, atom_number_t zatom);
    int _is_nitrogen_double_bond_to_something_outside_ring (atom_number_t zatom);
    int _is_aromatic_carbonyl (atom_number_t zatom,
                                 atom_number_t & a1,
                                 atom_number_t & a2) const;
    int _is_furan_or_thiophene (atom_number_t zatom,
                                  atom_number_t & a1,
                                  atom_number_t & a2) const;
    int _is_pyrrole_type_nitrogen (atom_number_t zatom,
                                     atom_number_t & a1,
                                     atom_number_t & a2,
                                     atom_number_t & a3,
                                     const int * process_these_atoms) const;
    int _must_be_single_bond_between (Kekule_Temporary_Arrays & kta,
                                      atom_number_t stop_atom,
                                      atom_number_t zatom,
                                      atom_number_t aprev);
    int _identify_continuation_atom (atom_number_t stop_atom, atom_number_t zatom,
                                       atom_number_t aprev,
                                       const int * aromatic_atom) const;

    int _process_without_kekule_perception (const int * aromatic_atoms,
                                            const int * aromatic_bonds);

    int _bond_in_ring_with_all_members_permanent_aromatic (atom_number_t a1,  
                                                             atom_number_t a2,
                                                             const int * is_permanent_aromatic) const;
    int _both_double_bonds_in_set(const Set_of_Atoms & p, const atom_number_t zatom) const;

    int _has_kekule_forms(const Set_of_Atoms & r) const;

/* arch-tag: 37e305b3-dfa5-44ea-8f3c-c036cab21fa5 */
