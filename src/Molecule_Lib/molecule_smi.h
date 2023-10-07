#ifndef COMPILING_SMILES_CC
  YIPES! THIS FILE IS ONLY SUPPOSED TO BE INCLUDED IN SMILES.CC
#else
  private:

    int _smiles_add_bond (atom_number_t previous_atom,
                            atom_number_t current_atom,
                            bond_type_t bt);

    int _smiles_choose_next_atom (const int * zorder, 
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom);
    int _smiles_choose_nice_next_atom (const int * zorder, 
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom);
    int _smiles_choose_random_next_atom (const int * zorder, 
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom);
    int _smiles_choose_unique_next_atom (const int * zorder, 
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom);

    int _smiles_choose_first_atom (const int * zorder, atom_number_t & first_atom, const int * include_atom);
    int _smiles_choose_nice_first_atom (const int * zorder, atom_number_t & first_atom, const int * include_atom);
    int _smiles_choose_unique_first_atom (const int * zorder, atom_number_t & first_atom);
    int _smiles_choose_unique_first_atom (const int * zorder, atom_number_t & first_atom, const int * include_atom);
    int _smiles_choose_random_first_atom (const int * zorder, atom_number_t & first_atom, const int * include_atom);

    template <typename F> int _smiles_choose_first_atom_fctr (F chooser, const int * zorder, atom_number_t & first_atom, const int * include_atom);

//  Various functions for building smiles

    const Bond * _identify_first_smiles_bond (atom_number_t zatom,
                                       const int * zorder);

    int _construct_smiles_for_fragment (Smiles_Formation_Info &, Smiles_Information & smi_info);
    int _construct_smiles (const int *, int *, Smiles_Information &, const int *);
    int _construct_smiles (const Fragment_Information &, Smiles_Information &, const int * include_atom);

    template <typename N> int _build_smiles_ordering_fctr (N next_atom_selector,
                                       const int * include_atom,
                                       Smiles_Information & smi_info);
    template <typename N> int _build_smiles_ordering_fctr (N next_atom_selector,
                                  atom_number_t previous_atom,
                                  atom_number_t zatom,
                                  int & icounter,
                                  const int * include_atom,
                                  Smiles_Information & smi_info);

    int _build_smiles_ordering (Smiles_Information & smi_info, const int * include_atom);


    int _append_smarts_equivalent (Smiles_Formation_Info & smi_info, IWString & s);

    int   _process_atom_for_smiles (Smiles_Formation_Info &, IWString & smiles);
    int   _process_atom_for_smiles (Smiles_Formation_Info &,
                  const int * zorder,
                  const resizable_array<const Bond *> &,
                  const resizable_array<atom_number_t> &,
                  const Chiral_Centre * c,
                  IWString &);
    int _update_ring_membership (const Ring * r);

    const IWString & _unique_smiles (const Fragment_Information & frag_info, Smiles_Information & smi_info, Symmetry_Class_and_Canonical_Rank & sccr, const int * include_atom);


    int _smiles_choose_first_atom (const int * zorder,
                                   Smiles_First_Atom & smfa,
                                   atom_number_t & first_atom,
                                   const int *);

//  int _all_atoms_are_chain_atoms (const int * process_these_atoms);

    int _determine_ring_closure_bonds (const int * zorder,
                                        const int * include_atom);

    atom_number_t _choose_highest_canonical_order_in_fragment (int f,
                                                       const int * zorder,
                                                       const int * include_atom) const;
    int _determine_ring_closure_bonds (atom_number_t aprev,
                                         atom_number_t zatom,
                                         const int * zorder,
                                         const int * include_atom,
                                         int * already_done);

//  in frag.cc

    int _recursive_fragment_membership(Fragment_Information& fragment_information);
    int _recursive_fragment_membership(Fragment_Information& fragment_information,
                atom_number_t zatom,
                int fragment_number,
                int & bonds_in_fragment);
    int _compute_fragment_information_subset(Fragment_Information & fragment_information,
                                               const int * include_atom) const;
    int _create_bond_subset_starting_with (Molecule & subset,
                                             atom_number_t zatom,
                                             int * bond_lookup_table,
                                             int * xref) const;

    int _compute_fragment_information (Fragment_Information &, int = 1);

//  Functions for transferring directional bonds in subsets

    int _transfer_directional_bond_info (Molecule & subset, const int * xref) const;
    int _directional_bond_at_other_end_in_subset(atom_number_t zatom, const int * xref) const;
    void _set_directional_bonds (Molecule & subset,
                                  atom_number_t zatom,
                                  const int * xref) const;
    int _count_directional_attachments (atom_number_t zatom) const;

    int _do_unfix_implicit_hydrogens_on_aromatic_isotopic_atoms(const int * aromatic_atoms);

    void _invalidate_attatched_directional_bonds (Molecule & subset,
                                                   const int * xref,
                                                   const Bond * b) const;
    int _transfer_wedge_bond_info (Molecule & subset, const int * xref) const;

    int _assign_fsid_values_to_isolated_rings ();

    int _common_fragment_extraction(int * tmp,
                                      const int initial_natoms,
                                      const int atoms_in_residual,
                                      const int atoms_in_fragment,
                                      Molecule & f);

    int SmilesSetName(const char * s, int nchars, int processing_quoted_smiles);
    int MaybeParseAsChemaxonExtension(const_IWSubstring& name, int processing_quoted_smiles);
    int ParseChemaxonExtension(const const_IWSubstring& chemaxon);
    int ParseCoords(const const_IWSubstring& chemaxon, int * claimed);
    int ParseSpecialAtoms(const const_IWSubstring& chemaxon, int * claimed);
#endif
