#ifndef MDL_MOLECULE_H
#define MDL_MOLECULE_H

#include "iwstring_data_source.h"

/*
  An MDL_Molecule carries with it all the query information read from an MDL file
*/

#include "substructure.h"
#include "mdl_file_data.h"
#include "mdl.h"
#include "ematch.h"
#include "atom_alias.h"

class MDL_Molecule : public Molecule, public MDL_File_Data
{
  private:

//  If we have something like explicit R atoms, or a query that
//  defines the atoms at which substitutions can be, or someone is
//  doing just substitute at isotopic atoms, then we may know the
//  id's of the atoms at which substitution can happen

    Set_of_Atoms _substitution_points;

//  private functions

    int _parse_link_record (const IWString & buffer, ::resizable_array_p<ISIS_Link_Atom> & ltmp);

    int _parse_M_record (iwstring_data_source & input,
                         const const_IWSubstring & buffer,
                         ::resizable_array_p<ISIS_Link_Atom> & link_atom,
                         int & fatal);
    int _common_set_from_aprop (const Aprop * atom_properties,
                                                int ntokens,
                                                int * dest);
    int _parse_atom_list (const IWString & buffer);
    int _parse_atom_alias(iwstring_data_source & input, const const_IWSubstring & buffer);
    int _set_unsaturation_specifications (const Aprop * atom_properties, int tokens);
    int _set_substitution_specifications (const Aprop * atom_properties, int tokens);
    int _set_ring_bond_specifications (const Aprop * atom_properties, int tokens);
    int _remove_explicit_link_atoms ();

    int _read_v3000 (iwstring_data_source & input);
    int _parse_v3000_atom_record (const const_IWSubstring &);
    int _read_v3000 (iwstring_data_source & input, int na, int nb);
    int _convert_symbol_to_element (int ndx, const IWString & s) ;
    int _look_for_atom_query_directives (int ndx, const IWString & buffer);
    int _look_for_bond_query_directives (int ndx, const IWString & buffer);

  public:
    MDL_Molecule();
    MDL_Molecule(const Molecule &);
    MDL_Molecule(const MDL_Molecule &);
    ~MDL_Molecule();

    const ISIS_Atom_List * atom_list_for_atom (atom_number_t a) const;

    int read_molecule_ds (iwstring_data_source &, int notused);
    int read_molecule_mdl_ds (iwstring_data_source &, int return_on_m_end = 0);

    int initialise_mqs (Molecule_to_Query_Specifications & mqs) const;

//  In reality, I should overload every function that can change the atom ordering, but
//  that would be too much work. Just beware and implement what you need... None of
//  this would be necessary if the info in an MDL_File_Data object were atom properties...

    int remove_atom (atom_number_t zatom);
    int remove_atoms (const int *);

    int remove_explicit_hydrogens(atomic_number_t = 1);
    
    int change_R_groups_to_substitutions (Element_Matcher & rgroup, int enable_hydrogen_substituent);
    int change_R_groups_to_match_any_atom (Element_Matcher & rgroup, int only_substituents_at_matched_atoms);

    int only_allow_substitutions_at_isotopic_atoms(const Molecule_to_Query_Specifications & mqs);
    int only_allow_substitutions_at_non_isotopic_atoms();

    void set_substitution(const atom_number_t a, const int s);
    void set_ring_bond(const atom_number_t a, const int s);

    int determine_attachment_points_by_query(Molecule_to_Query_Specifications & mqs);
    
    int not_atom_lists_present() const;

    int swap_atoms (atom_number_t, atom_number_t);

    const MDL_Bond_Data * mdl_bond_between_atoms(atom_number_t, atom_number_t) const;

    int add (const Element *);

    int add_atom_based_on_symbol (const IWString & s);

    int add_bond (atom_number_t a1, atom_number_t a2, bond_type_t bond_for_molecule,
                        bond_type_t query_bond, int bond_topology);

    int remove_bond_between_atoms (atom_number_t, atom_number_t);
    
    int set_mdl_atom_data (int, const MDL_Atom_Data *);

    int compute_aromaticity_handle_atom_lists();   // handles atom lists in rings
};

/*
  An optional behaviour is to automatically convert A and Q atoms to atom lists with just
  the organic subset. This is for Valhalla
*/

extern void set_convert_a_and_q_atoms_to_atom_lists(int);
extern void set_convert_not_atom_lists_to_organic_lists(int);
void reset_mdl_molecule_file_scope_variables();
extern void set_mdl_molecule_discard_chirality(const int s);

#endif
