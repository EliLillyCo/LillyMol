#ifndef COMPILING_MOLECULER_CC
   YIPES, THIS FILE SHOULD ONLY BE INCLUDED WHEN COMPILING MOLECULER_CC
#endif

    int _make_ring (const Beep * bp, Ring * r, int * atom_in_ring);
    int _make_rings (const resizable_array_p<Beep> & beeps,
                       resizable_array<Ring *> & rings,
                       int * tmp);
    int _make_rings (const resizable_array_p<Beep> & beeps,
                       resizable_array_p<Ring> & rings,
                       int * tmp);
    int _make_rings (Rings_Found & ring_info,
                       resizable_array_p<Ring> & sssr_rings,
                       resizable_array_p<Ring> & non_sssr_rings,
                       int fid,
                       int * tmp);
    int _unused_fused_system_identifier () const;

    int _make_rings (Rings_Found & rings_found,
                       resizable_array_p<Ring> & sssr_rings,
                       resizable_array_p<Ring> & non_sssr_rings);

    int _convert_fused_raw_rings_to_sssr_form (int fused_sys_id, int * tmp);

    int _sssr_for_all_raw_rings (int * tmp);
    int _pearlman_sssr (const int * process_these, int id, Tnode ** tnodes, int);
    int _pearlman_sssr (const int * process_these, int id, Tnode ** tnodes, const int * pi, const int * sac);
    int _pearlman_sssr (const int * process_these, int id, Tnode ** tnodes);
    int _pearlman_sssr (const int *, int);

    int _initialise_tnode (atom_number_t zatom, const int * process_these, int id, Tnode ** tnodes, const int * pi, const int * sac);

    int _transfer_from_non_sssr_to_sssr_ring_set (int nssr_ndx, int sssr_ndx);

    int _identify_spinach (int * spinach, Set_of_Atoms & stack) const;

    void _sort_by_ring_size ();

  public:
    int experimental_sssr ();
    int experimental_nrings (atom_number_t);
    const int * experimental_ring_membership ();
    int experimental_print_ring_info (std::ostream &) const;

    int experimental_bonds_between (atom_number_t, atom_number_t);

  private:

    int _bonds_to_ring (const int *, Ring *);
    int _set_bonds_in_ring (IW_Bits_Base *, const Ring * r);
    int _set_bonds_in_ring (int * zbonds, const Ring * r);
    int _convert_to_sssr_form (const resizable_array_p<int> & rings_as_bonds, int * tmp, int ** matrix_form);
    int _convert_to_sssr_form (const resizable_array_p<int> & rings_as_bonds, int * tmp);
    int _convert_to_sssr_form (resizable_array_p<Ring> & raw_rings);

    int _merge_fused_system_identifiers (resizable_array<Ring *> & rings,
                                         int rstart,
                                         int fused_system_identifier,
                                         resizable_array<int> & fused_sys_ids_to_be_changed);


    int _unused_fused_ring_system_identifier ();
    int _find_sssr_for_these_fused_raw_rings (int fused_sys_id, int * tmp);
    int _convert_all_raw_rings_to_sssr_form ();
    int _force_complete_sssr_determination ();
    int _determine_sssr_ring_membership (atom_number_t a);
    int _determine_sssr_for_fragment (int f);

//  All the stuff from iwring

    int _determine_unique_rings (int fused_sys_id,
                                 int max_rings_to_find_this_iteration,
                                 IW_Bits_Base ** matrix_form,
                                 resizable_array_p<Ring> & raw_rings_found_this_iteration);
    int _determine_unique_rings (int fused_sys_id,
                                 int max_rings_to_find_this_iteration,
                                 int ** matrix_form,
                                 resizable_array_p<Ring> & raw_rings_found_this_iteration);
    int _go_from_to (atom_number_t afrom, atom_number_t ato,
                     atom_number_t aexclude,
                     int id, const int * process_these,
                     Ring * r);
    int _construct_new_raw_ring (int id, const int * process_these,
                            atom_number_t i1,
                            atom_number_t i2,
                            atom_number_t i3,
                            Ring * r);
    int _found_ring (int id, const int * process_these,
                      atom_number_t i1,
                      atom_number_t i2,
                      atom_number_t i3,
                      resizable_array_p<Ring> & raw_rings);
    int _advance_sssr_determination (int id, const int * process_these,
                                     const int * con,
                                     int distance,
                                     resizable_array_p<Ring> & rings_found_this_iteration);
    int _form_sssr (int id, const int * process_these,
                    const int * con);
    int _form_sssr (int id, const int * process_these);
    int _next_available_fused_system_identifier () const;

    void _assign_ring_numbers (int);

    int _handle_spiro_between_isolated_and_fused (const resizable_array<Ring *> & rings, int, int *);

    void _add_ring_to_sssr (Ring * ri);

    int _easy_case_two_rings (int fused_sys_id, int * tmp);

    int _compute_number_sssr_rings_by_eulers_formula();

  private:

    int _just_one_unclassified_spinach_connection (atom_number_t zatom, const int * spinach) const;
    int _pi_electrons_in_ring(const atom_number_t zatom, int & result) const;
