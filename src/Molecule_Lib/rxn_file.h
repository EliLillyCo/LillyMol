#ifndef IWRXN_FILE_H
#define IWRXN_FILE_H

#include "sparse_fp_creator.h"
#include "iw_stl_hash_map.h"

#include "iwreaction.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include "mdl_molecule.h"

class Atom_Typing_Specification;
class ISIS_RXN_FILE_Molecule;
class Product_Atom_Types;

/*
  When creating a reaction from a subset of an input file, there
  is quite a lot of book-keeping to keeping track of which atom
  goes where.
*/

class Reaction_Subset
{
  private:
    const int _nr;

    int * _xref;
    int * _reagent_offset;

//  Turns out that sometimes we need a cross reference array, and sometimes we need an include/exclude
//  flag, so we have two variables. 

    int * _include_atom;

    int _total_atoms;

  public:
    Reaction_Subset(const int nr);
    ~Reaction_Subset();

    int debug_print(std::ostream &) const;

    const int * include_atom (const int i) const { if (NULL == _include_atom)
                                                     return NULL;
                                                   return _include_atom + _reagent_offset[i];}

    void convert_cross_reference_to_boolean_include();

    int build(const ISIS_RXN_FILE_Molecule * reagent, const int * reagent_locator,
                       const int * include_these_atom_map_numbers, const int highest_atom_map);

    atom_number_t atom_number_in_subset (const atom_number_t a, const int r) const;
};

struct Aprop;
class Atom_in_Fragment
{
  private:
    atom_number_t _a;
    int           _fragment;
  public:
    Atom_in_Fragment (atom_number_t zatom, int f) : _a (zatom), _fragment (f) {};

    atom_number_t atom () const { return _a;}
    int fragment () const { return _fragment;}

    void set_atom (atom_number_t s) { _a = s;}
    void set_fragment (int s) { _fragment = s;}
};

class RXN_File_Create_Reaction_Options
{
  public:
    int _only_create_query_from_first_reagent;
    int _write_agent;

  public:
    RXN_File_Create_Reaction_Options();
    ~RXN_File_Create_Reaction_Options(){};

    void set_only_create_query_from_first_reagent(const int s) { _only_create_query_from_first_reagent = s;}

    int only_create_query_from_first_reagent() const { return _only_create_query_from_first_reagent;}
   
    
};

class Reaction_Site;

class ISIS_RXN_FILE_Molecule : public MDL_Molecule
{
  private:
    int * _atom_map;

//  Because of the various bonding possibilities (like double or aromatic), it can
//  sometimes be challenging to compute the number of attached hydrogens. We do it once
//  and store the result

    int * _explicit_hydrogen_count;
    int * _computed_implicit_hydrogen_count;

//  We remove any explicit Hydrogens found in the input file. We need to
//  know how many explicit Hydrogens have been removed from each atom so we
//  can make inferences about hcount in the query

    int * _explicit_hydrogen_atom_removed;

//  With link atoms, we need to keep track of the atoms that have lost a connection

    int * _connections_lost;

//  Some of the bonds will be valid bonds, but some will be OR bonds. Store all
//  bonding info

    bond_type_t * _btype;

//  NOT atom lists need special care so check whether any are present or not

    int _not_atom_lists_present;

//  If we have OR aromatic bonds present, some of the rings drawn may be aromatic,
//  but since we've put single bonds into our Molecule, we have lost the aromaticity.
//  The funciton _look_for_rings_that_are_supposed_to_be_aromatic discerns hidden aromaticity

    int * _aromatic;

//  People expect to enter A atoms in aromatic rings, but since A is a non-periodic table
//  element, such a ring cannot be aromatic. We temporarily convert to Carbon to get aromaticity

    int _convert_A_to_C_for_aromaticity;

//  When deciding whether we need to toggle Kekule forms, we need to know whether
//  ring membership changes or not

    resizable_array_p<Bond> _toggle_kekule_form;

//  For some reactions, it is better to forget Kekule forms.

    int _aromatic_bonds_lose_kekule_identity;

//  Aug 2016. Sometimes we really do want to keep the existing Kekule form

    int _preserve_kekule_forms;

//  Since I strip off explicit Hydrogens, it becomes much easier to compute the number of
//  bonds before any of that happens

    int * _nbonds;

//  When doing automatic atom mapping, we need to know symmetry. If an atom
//  is set in this array, it is either unique (value 1) or a representative
//  of a symmetry class, with the value being the number of symmetrically equivalent atoms
    
    int * _use_as_representative;

//  We can speed things up by putting a heteroatom first

    int _swap_atoms_to_put_rare_atoms_first;

//  Sometimes we will have a reaction where the thing drawn is the (only) reagent,
//  rather than a query

    int _molecule_is_sole_reagent;

    int _remove_explicit_hydrogens;

//  If we have identical molecules as products, then _use_as_representative needs to span
//  molecules. Fragments that have the same structure will have the same non-zero structure_value

    int _structure_value;

    int _is_duplicate_fragment;

    IWString _unique_smiles;

    int _interpret_atom_alias_as_smarts;   // sept 2014

//  It can be convenient to do atom typing across a reaction - does not come from a file

    int * _atype;

//  private functions

    int _compute_implicit_hydrogens ();

    int _aromatic_bonds_attached (atom_number_t zatom) const;

    int _number_atoms_with_substitution_directives() const;
    int _number_atoms_with_unsaturation_directives() const;
    int _number_atoms_with_ring_bond_directives() const;

    int _adjust_computed_implicit_hydrogens_for_unstautration (atom_number_t zatom);

    int _set_unsaturation_specifications (const Aprop * atom_properties, int ntokens);
    int _set_substitution_specifications (const Aprop * atom_properties, int ntokens);
    int _set_ring_bond_specifications (const Aprop * atom_properties, int ntokens);

    int _grow_symmetry_class (atom_number_t zatom, int * symmetry_equivalent_atoms);

    int _common_set_from_aprop (const Aprop * atom_properties, int ntokens, int * dest);

    int _parse_atom_list (const IWString &);
    int _parse_link_record (const IWString &);

    int _check_non_periodic_table_elements_and_atom_lists ();
    int _create_query (Reaction_Site & r,
                           Molecule_to_Query_Specifications & mqs,
                           int ndx);
                           
                         
    int _create_query (Reaction_Site & r,
                           Molecule_to_Query_Specifications & mqs, 
                           const int * include_these_atoms);

//  int _do_remove_explicit_hydrogens (Set_of_Atoms & contains_explicit_hydrogen);
    int _identify_explicit_hydrogens_to_be_removed (Set_of_Atoms & to_be_removed, int * xref);
    int _identify_link_atoms_to_be_removed (Set_of_Atoms & to_be_removed, int * xref);

    int _do_atom_removals (Set_of_Atoms & to_be_removed, const int * xref);

    int _look_for_rings_that_are_supposed_to_be_aromatic ();
    int __look_for_rings_that_are_supposed_to_be_aromatic (int *);

    int _back_to_single_bonds (bond_type_t);

    int _change_explicit_kekule_forms_to_aromatic (const Ring & r);
    int _change_explicit_kekule_forms_to_aromatic ();

    int _do_aromatic_bonds_lose_kekule_identity ();

    int _do_swap_atoms_to_put_rare_atoms_first ();

    int _write_m_sub_records (int n, std::ostream & output) const;
    int _write_m_uns_records (int n, std::ostream & output) const;
    int _write_m_rbc_records (int n, std::ostream & output) const;
		
  public:
    ISIS_RXN_FILE_Molecule ();
    ~ISIS_RXN_FILE_Molecule ();

    ISIS_RXN_FILE_Molecule & operator=(ISIS_RXN_FILE_Molecule&& rhs);

    int do_read (iwstring_data_source &);

    int allocate_arrays (int na, int nb);

    int highest_atom_map_number () const;

    int remove_atoms (const int * to_remove);
    int remove_atom (const int * to_remove);

    int identify_which_atoms_we_have (int * locator_array, int mark, int * atom_map_to_atom_number) const;

    void recompute_unique_smiles ();

    int add (const Element *);

    int identify_symmetry_classes ();

    int use_atom_as_symmetry_representative (atom_number_t a) const { return _use_as_representative[a];}
    void set_use_as_symmetry_representative (atom_number_t a, int s) { _use_as_representative[a] = s;}

    int fill_bonds_between_mapped_atoms_array (atom_type_t * b, const int);

    int break_symmetry (atom_number_t);

    void discard_atom_map ();

//  Sometimes callers will need to know whether a particular atom is a list

    int is_an_atom_list (atom_number_t) const;

    int structure_value () const { return _structure_value;}
    void set_structure_value (int s) {  _structure_value = s;}

    int is_duplicate_fragment () const { return _is_duplicate_fragment;}
    void set_is_duplicate_fragment (int s) { _is_duplicate_fragment = s;}

    void set_interpret_atom_alias_as_smarts(int s) { _interpret_atom_alias_as_smarts = s;}

    const IWString & unique_smiles () const { return _unique_smiles;}

    bond_type_t mapped_atoms_are_bonded (int, int) const;

    const int * atom_map () const { return _atom_map;}
    int atom_map (atom_number_t a) const { return _atom_map[a];}

    int number_mapped_atoms () const;

    void set_atom_map (atom_number_t a, int m) {assert (ok_atom_number (a)); _atom_map[a] = m;}

    int transfer_atom_map_data_to_global_array();

    int allocate_mdl_atom_array();

    void set_aromatic_bonds_lose_kekule_identity (int s) { _aromatic_bonds_lose_kekule_identity = s;}
    void set_preserve_kekule_forms (int s) { _preserve_kekule_forms= s;}
    void set_swap_atoms_to_put_rare_atoms_first (int s) { _swap_atoms_to_put_rare_atoms_first = s;}
    void set_molecule_is_sole_reagent (int s) { _molecule_is_sole_reagent = s;}
    void set_remove_explicit_hydrogens (int s) { _remove_explicit_hydrogens = s;}
    void set_convert_A_to_C_for_aromaticity (int s) { _convert_A_to_C_for_aromaticity = s;}

    int single_reagent_only () const { return _molecule_is_sole_reagent;}

    int which_is_mapped_atom (int) const;

    int identify_unmapped_neighbour  (int m, atom_number_t & n, int &);
    int identify_unmapped_neighbours (int m, Set_of_Atoms & neighbour);

    int identify_mapped_neighbours (atom_number_t zatom, resizable_array<int> & neighbour) const;

    int identify_singly_bonded_mapped_neighbours (atom_number_t zatom, resizable_array<int> & neighbour) const;

    int identify_connected_mapped_atoms (atom_number_t zatom, resizable_array<int> & connected_to) const;

//  int transfer_isis_info_to_mqs (Molecule_to_Query_Specifications &) const;

    int swap_atoms_to_put_rare_atoms_first ();

    int create_query (Reaction_Site & r, 
    									const int * include_these_atoms, 
    									std::ofstream *queryOutStream);
    int create_query (Reaction_Site & r, 
    									const int * include_these_atoms,
    									Molecule_to_Query_Specifications & mqs, 
    									std::ofstream *queryOutStream);

    int add_chiral_centres_to_be_inverted (Reaction_Site &) const;

    int nbonds_for_mapped_atom (int);

    int discern_initial_conditions (int m, int final_nbonds);

    int do_write (std::ostream &) const;

    int add_toggle_kekule_form (atom_number_t, atom_number_t, bond_type_t);

    int add_toggle_kekule_forms (Reaction_Site & rxn, const Reaction_Subset & subset, const int ndx,
    												const RXN_File_Create_Reaction_Options & rxnfcro) const;

    int has_inter_fragment_changes () const;

    int transfer_atom_alias_to_isotope ();

    int assign_atom_types (Atom_Typing_Specification &);

    const int * atom_type () const { return _atype;}

    int gather_mapped_neighbours (const atom_number_t zatom, resizable_array<int> & nbrs) const;
    
};

/*
  When working out the number of changing atoms in a reaction, we may want to impose
  some conditions on what gets counted
*/

class Changing_Atom_Conditions
{
  private:

// feb 2016. An atom is changing if the mapped atoms attached to it change during the reaction

    int _is_changing_if_different_neighbours;

// may 2016. Sometimes we want to be able to ignore small fragments

    int _ignore_lost_atom_if_isolated;

// jul 2016. We can have atoms in small fragments that really do not matter

    int _only_consider_largest_reagent_fragment;

//  Sep 2016. Sometimes it is useful to include changing bonds in determining changing atoms,
//  but frequently not. Especially in the case where the only difference might be a Kekule 
//  representation

    int _include_changing_bonds_in_changing_atom_count;

//  Oct 2016. If we have previously added orphan atoms to the reagents, we may want to
//  exclude those atoms when computing the changing atom count

    int _discern_changing_atoms_only_in_first_fragment;

//  Nov 2017. Sometimes we only consider bond aromaticity

    int _consider_aromatic_bonds;

  public:
    Changing_Atom_Conditions();

    void set_is_changing_if_different_neighbours (int s) { _is_changing_if_different_neighbours = s;}
    void set_ignore_lost_atom_if_isolated (const int s) { _ignore_lost_atom_if_isolated = s;}
    void set_only_consider_largest_reagent_fragment (const int s) { _only_consider_largest_reagent_fragment = s;}
    void set_include_changing_bonds_in_changing_atom_count(const int s) { _include_changing_bonds_in_changing_atom_count = s;}
    void set_discern_changing_atoms_only_in_first_fragment(const int s) { _discern_changing_atoms_only_in_first_fragment = s;}
    void set_consider_aromatic_bonds(const int s) { _consider_aromatic_bonds = s;}

    int is_changing_if_different_neighbours () const { return _is_changing_if_different_neighbours;}
    int ignore_lost_atom_if_isolated () const { return _ignore_lost_atom_if_isolated;}
    int only_consider_largest_reagent_fragment () const { return _only_consider_largest_reagent_fragment;}
    int include_changing_bonds_in_changing_atom_count() const { return _include_changing_bonds_in_changing_atom_count;}
    int discern_changing_atoms_only_in_first_fragment() const { return _discern_changing_atoms_only_in_first_fragment;}
    int consider_aromatic_bonds() const { return _consider_aromatic_bonds;}
};

/*
  When looking across from reagents to products, we need a data structure to make those
  cross references faster. Otherwise we have to use the atom map to find the product locator,
  then ask that product to fetch the atom number with that atom map.

  For each atom map encountered in reagents, we store a 32 bit number. The first 16 bits are
  the product number, and the last 16 bits are the atom number
*/

class Atom_Locations
{
  private:
    std::unordered_map<int, unsigned int> _reagent_to_product;

  public:
    Atom_Locations();
    ~Atom_Locations();

    int initialise(ISIS_RXN_FILE_Molecule * reagents, const int nr, ISIS_RXN_FILE_Molecule * products, const int np, const int * product_locator);

    int get_product_and_atom_number(const int amap, int & product, atom_number_t & zatom) const;
};

class Reaction_Smiles_Options
{
  public:
    int _reagent_product_plus_rather_than_dot;
    //int _orphan_plus_rather_than_dot;
    int _write_reaction_name;
    int _write_agent;
    //IWString _output_separator;

  public:
    Reaction_Smiles_Options();

    void set_reagent_product_plus_rather_than_dot(const bool s) { _reagent_product_plus_rather_than_dot = s;}
    //void set_orphan_plus_rather_than_dot(const int s) { _orphan_plus_rather_than_dot = s;}
    void set_write_reaction_name(const int s) { _write_reaction_name = s;}
    //void set_output_separator(const const_IWSubstring & s) { _output_separator = s;}
    void set_write_agent (const int s) { _write_agent = s;}

    bool reagent_product_plus_rather_than_dot() const { return _reagent_product_plus_rather_than_dot;}
    //int orphan_plus_rather_than_dot() const { return _orphan_plus_rather_than_dot;}
    int write_reaction_name() const { return _write_reaction_name;}
    int write_agent() const { return _write_agent;}
    //const IWString & output_separator() const { return _output_separator;}
};

class RXN_File
{
  private:
    IWString _fname;

    IWString _comment;

    int _nr;
    ISIS_RXN_FILE_Molecule * _reagent;

    int _np;
    ISIS_RXN_FILE_Molecule * _product;

    int _na;
    ISIS_RXN_FILE_Molecule * _agent;

//  If we have a reaction that produces atoms that are not present in the LHS
//  we need somewhere to store those atoms. I do not want to store those
//  atoms with any reagent since they would become part of the query

    ISIS_RXN_FILE_Molecule _orphan_atoms;

//  It is convenient to know the reagent for each mapped atom. 

    int * _reagent_locator;

    int * _product_locator;

//  for a given atom map number, what is the corresponding atom number. Note that
//  the atom number will be in whatever reagent/product has that atom map

    int * _atom_map_to_reagent_atom_number;
    int * _atom_map_to_product_atom_number;

//  Saves time to know the bonding status both before and after. But because we may extend the atom
//  map, we need to keep the dimensionality of these arrays

    int _btype_dim;
    bond_type_t * _initial_bond_type;
    bond_type_t * _final_bond_type;

//  We can optionally remove all product fragments except the first

    int _remove_product_fragments;

//  These variables are only used during initial construction. We put them in the
//  object just to avoid passing them around as arguments

    int _unmapped_atoms_in_reagents;
    int _unmapped_atoms_in_products;

//  What if there is an unmapped atom of a given kind on the LHS but absent on the RHS

    int _remove_unmapped_atoms_that_disappear;

    int _aromatic_bonds_lose_kekule_identity;
    int _preserve_kekule_forms;

    int _swap_atoms_to_put_rare_atoms_first;

//  Since we get told which reagent is a single reagent before we have read the .rxn file,
//  we don't know the max number of reagents in advance. Therefore an extending array...

    extending_resizable_array<int> _molecule_is_sole_reagent;

//  If we need to echo the reaction after atom mapping

    IWString _fname_for_echo;

    int _do_automatic_atom_mapping;

//  make removal of explicit Hydrogens optional

    int _remove_explicit_hydrogens;

//  Need to identify cases where a rearrangement across two bonds takes place

    int * _involved_in_square_bond;

//  We can optionally unconnect unmapped atoms that, which if left in place,
//  would violate the valence of the product.
  
    int _unconnect_unmapped_atoms_that_exceed_product_valence;

    int _convert_A_to_C_for_aromaticity;

    int _convert_atom_aliases_to_isotopes;

//  Sept 2014. When processing files from Marvin, we need to NOT interpret
//  atom aliases as smarts

    int _interpret_atom_alias_as_smarts;

//  Sept 2014. What do we do when we see a reaction with an atom list on the LHS
//  and a carbon on the RHS. Did they mean to transform whatever matched in the
//  reagent into a carbon atom, or was the carbon atom just a placeholder.

    int _interpret_carbon_on_rhs_of_list_as_no_change;

    int _auto_fix_orphans;     // jan 2016

// Mar 2016. When creating subsets, we need to pass info to the molecule_to_query object that
// controls generation of the substructure query. But the Molecule_to_Query_Specifications object is 
// not visible outside. Rather than changing a whole bunch of signatures to allow it, we have
// a bit of a kludge

    int _mol2qry_isotope_special_meaning;
    
// Mar 2018.  retrosynthetic_quick had a core dump when atoms that are parts of a bond that only changed
// its kekule form were NOT considered changed atoms. However, this code is used in many programs, and could
// cause changes in expected hehavior.  A complete test of regression tests did NOT reveal any issues, so the default
// has been set to make these atoms as changed.
//  Tad Hurst

    int _mark_atoms_changed_when_kekule_form_of_bond_changes;

// 
    std::ofstream *_queryOutStream;
    	
//  private functions

    int _identify_square_bonding_changes (int highest_atom_map);
    int _identify_square (int m1, int m2, int m3) const;

    int _create_reaction (IWReaction & rxn,
                          const RXN_File_Create_Reaction_Options & rxnfcro,
                          Molecule_to_Query_Specifications & mqs,
                          int highest_atom_map, const int * include_these_atoms);

    int _looks_like_atom_replacement (int m,
                                       resizable_array<int> & initial_connections,
                                       resizable_array<int> & final_connections) const;

    int _look_for_unmapped_atoms_that_exceed_product_valence (int highest_atom_map, IWReaction & rxn) const;

    int _identify_bonds_to_be_broken (int highest_atom_map,
                                      resizable_array_p<Bond> & bonds_to_be_broken,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions,
                                      const int * include_these_atom_map_numbers) const;
    int _identify_bonds_to_be_made (int highest_atom_map,
                                      resizable_array_p<Bond> & bonds_to_be_made,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions,
                                      const int * include_these_atom_map_numbers);
    int _identify_isotopes_to_be_placed (int highest_atom_map, resizable_array_p<Reaction_Place_Isotope> & isotopes) const;

    int _identify_atom_substitutions (int highest_atom_map,
                                         resizable_array<Replace_Atom *> & atom_substitutions,
                                         resizable_array_p<Bond> & connections_involved_in_substitutions) const;
    void _specify_stereo_centre_component (Stereo_Centre_Component & s,
                                            int centre_atom_reagent,
                                            int our_reagent,
                                            int mapped) const;
		int _remove_unmapped_components (ISIS_RXN_FILE_Molecule * component,  int & n);

#ifdef COMPILING_RXN_FILE
    int  _number_reagent_atoms () const;
    int  _number_product_atoms () const;

    int _parse_ChemAxon_Extensions(const_IWSubstring buffer, IWString & smiles);

    int _establish_atom_mapping (int & highest_atom_map_number);

    int _establish_atom_mapping (const IW_Hash_Map<int, int> & unmapped_occurrences_in_reagents,
                                 const IW_Hash_Map<int, int> & unmapped_occurrences_in_products,
                                 int & highest_atom_map_number);
    int _mapped_atoms_bonded_in_products (int m1, int m2) const;

    int _reestablish_reagent_locator_array();
    int _reestablish_locator_array(ISIS_RXN_FILE_Molecule * x, const int n, int * locator);

    int _all_atoms_still_aromatic(ISIS_RXN_FILE_Molecule & m, const Ring & r);

    int _fix_kekule_differences_reagent(const int rgnt);
    template <typename CMP> int _identify_atoms_with_different_bonding (const ISIS_RXN_FILE_Molecule & m,
                                                  const Ring & r,
                                                  atom_number_t & a1,
                                                  atom_number_t & a2,
                                                  CMP compare,
                                                  const Set_of_Atoms & already_found) const;

//  int _map_corresponding_atom_list (int highest_atom_map_number, const ISIS_Atom_List & als);
    int _identify_corresponding_list_in_products (const ISIS_Atom_List & als, int &, atom_number_t &) const;
    int _map_any_atom_lists_on_either_side (int & highest_atom_map_number);

    int _map_symmetry_equivalent_atoms (int & highest_atom_map_number,
                                          int r,
                                          atom_number_t ar,
                                          int p,
                                          atom_number_t ap);
    int _set_atom_map (int & highest_atom_map_number,
                                          int r,
                                          atom_number_t ar,
                                          int p,
                                          atom_number_t ap);

    int _identify_kekule_forms_to_be_toggled (const int * include_these_atom_numbers);

    int _changes_in_ring_membership (ISIS_RXN_FILE_Molecule & mfrom, const Ring & r) const;

    void _fill_reagent_and_product_locator_arrays ();
    int _setup_reagent_product_locator_arrays();

    int _look_for_unmapped_atoms_that_disappear (int & highest_atom_map,  
    																			IWReaction & rxn, 
    																			const Reaction_Subset & subset,
    																			const RXN_File_Create_Reaction_Options & rxnfcro);
    																			
    int _identify_small_fragments_showing_up_in_products(const int rgnt, int * not_in_largest_fragment) const;
    int _all_atoms_in_fragment_in_products(ISIS_RXN_FILE_Molecule & r, const int f, int * not_in_largest_fragment) const;

    int _add_atom_removals (IWReaction & rxn, int h) const;
    int _add_atom_removal (Reaction_Site & r, int h) const;

    int _different_mapped_atoms_attached(const ISIS_RXN_FILE_Molecule & r, const atom_number_t zatom, const int zmap) const;

    int _map_unmapped_surrounded_by_mapped (int & highest_atom_map_number);
    int _identify_product_with_same_mapped_neighbours (int h,
                                                       const resizable_array<int> & mapped_neighbours,
                                                       int & pm,
                                                       atom_number_t & pa) const;

    int _do_read_v3000(iwstring_data_source & input);
    int _fix_orphan_condition (const Set_of_Atoms & orphans, const int p);
    int _identify_bonding_changes_involving_matched_atom (const int mstart,
                                      ISIS_RXN_FILE_Molecule & reagent,
                                      ISIS_RXN_FILE_Molecule & product,
                                      const int highest_atom_map,
                                      resizable_array_p<Bond> & bonds_to_be_made,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions,
                                      const int * include_these_atom_map_numbers) const;
    int _identify_bonds_to_be_made_involving_orphan (const resizable_array<int> & orphan,
                                      resizable_array_p<Bond> & bonds_to_be_made,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions);
    int _identify_bonds_to_be_made_involving_orphan(resizable_array_p<Bond> & bonds_to_be_made,
                                                      const resizable_array_p<Bond> & connections_involved_in_substitutions) const;
    int _identify_bonds_to_be_made_involving_orphan_atom(const int oi,     // atom map number of the orphan atom
                                                      resizable_array_p<Bond> & bonds_to_be_made,
                                                      const resizable_array_p<Bond> & connections_involved_in_substitutions) const;
    int _all_atom_maps_only_here (const ISIS_RXN_FILE_Molecule & m, const int * locator) const;

    int _remove_common_fragments(ISIS_RXN_FILE_Molecule & r, ISIS_RXN_FILE_Molecule & p);

    int _same_atoms_and_bonds(ISIS_RXN_FILE_Molecule & r, const int rf, ISIS_RXN_FILE_Molecule & p, const int pf) const;
    int _atoms_the_same(ISIS_RXN_FILE_Molecule & m1, const atom_number_t a1, ISIS_RXN_FILE_Molecule & m2, const atom_number_t a2) const;
    int _remove_non_participating_fragments (ISIS_RXN_FILE_Molecule & r);

    int _identify_atoms_changing_by_bond_reagent(ISIS_RXN_FILE_Molecule & r, int * changed, const Changing_Atom_Conditions & cac) const;;
    int _identify_just_changed_kekule_form(ISIS_RXN_FILE_Molecule & r, const atom_number_t zatom) const;
    int _identify_just_changed_kekule_form(ISIS_RXN_FILE_Molecule & r, const int * changed, int * just_changed_kekule_form) const;

    int _identify_atoms_in_rings_separated_from_changing_atoms(ISIS_RXN_FILE_Molecule & r, const int * changed, int * just_changed_kekule_form) const;
    int _identify_atoms_in_rings_separated_from_changing_atoms(ISIS_RXN_FILE_Molecule & reagent, const Ring & r, const int * changed, int * just_changed_kekule_form) const;

    int _involves_exocyclic_bonding_changes(const int rgnt, const Ring & r) const;
    int _bond_the_same_in_product(const ISIS_RXN_FILE_Molecule & r, const Bond & b) const;
    int _same_connectivity_and_bonds_in_product(const ISIS_RXN_FILE_Molecule & r, const atom_number_t x1) const;
    int _loss_or_gain_of_singly_bonded_atom (const int rgnt, const int am) const;
    template <typename T>
    int _reaction_fingerprint_rev(const int p,       // processing product molecule P
                                    const T * reagent_atom_type,
                                    const Product_Atom_Types & product_atom_type,
                                    const Set_of_Atoms & changing_in_reagent,
                                    const Set_of_Atoms & changing_in_product,
                                    Sparse_Fingerprint_Creator & sfc);
    template <typename T> int _form_possible_lost_or_changed_bond_fingerprint(const int m1, const int m2,
                                      const Bond * b1, const T * reagent_atom_type,      // in the reagent
                                      const Bond * b2, const Product_Atom_Types & pat,  // in the product
                                      Sparse_Fingerprint_Creator & sfc) const;
/*  template <typename T>
    int _form_possible_lost_or_changed_bond_fingerprint(const int m1, const int m2,
                                      const Bond * b1, const T * reagent_atom_type, // in the reagent
                                      const Bond * b2, const Product_Atom_Types & pat,  // in the product
                                      Sparse_Fingerprint_Creator & sfc) const;*/
    template <typename T>
    int _reaction_fingerprint_bonds(const T * reagent_atom_type,
                                      const Product_Atom_Types & product_atom_type,
                                      Sparse_Fingerprint_Creator & sfc) const;
#endif

    int _create_query (Reaction_Site & r, 
    										ISIS_RXN_FILE_Molecule & m, 
    										Molecule_to_Query_Specifications & mqs, 
    										const int * include_these_atoms);
                         
    int _look_for_stereo_centres_made (IWReaction &);

    int _map_unmapped_atoms (int & highest_atom_map_number);
    int __map_unmapped_atoms (int & highest_atom_map_number);

    int _aggressively_map_unmapped_atoms (int highest_atom_map_number);

    int _highest_atom_map_number_or_atoms_in_reagents () const;
    int _max_rings_in_any_reagent() const;
    int _reagent_ring_contains_adjacent_mapped_atoms (const ISIS_RXN_FILE_Molecule & m,
                                const Ring & r, atom_number_t & a1, atom_number_t & a2,
                                const int * include_these_atom_map_numbers) const;
    int _identify_changing_aromatic_systems(const int reagent_number, int * changed, const Changing_Atom_Conditions & cac) const;
    int _all_bonds_unchanged(const ISIS_RXN_FILE_Molecule & m, const Ring & r) const;
    void _update_aromatic_ring_system(ISIS_RXN_FILE_Molecule & m, const Ring & r, int * changed, const int flag, int * ring_already_done) const;

  public:
    RXN_File ();
    ~RXN_File ();

    const IWString & name() const { return _comment;}
    void set_name (const IWString & s) { _comment = s;}

    int debug_print (std::ostream &) const;
    int write_mapped_reaction (std::ostream & os) const;
    int print_atom_map_into (std::ostream & output) const;

    const IWString & fname () const { return _fname;}
    void set_fname (const char * s) { _fname = s;}

    int number_reagents() const { return _nr;}
    int number_products() const { return _np;}
    int number_agents() const { return _na;}

		void setQueryOutStream(std::ofstream *thisStream){_queryOutStream = thisStream;}
    //std::ofstream *queryOutStream() { return _queryOutStream;}
    	
    void set_remove_product_fragments (int s) { _remove_product_fragments = s;}
    void set_remove_unmapped_atoms_that_disappear (int s) { _remove_unmapped_atoms_that_disappear = s;}
    void set_aromatic_bonds_lose_kekule_identity (int s);
    void set_preserve_kekule_forms (int s);
    void set_swap_atoms_to_put_rare_atoms_first (int s) { _swap_atoms_to_put_rare_atoms_first = s;}
    void set_molecule_is_sole_reagent (int, int);
    void set_do_automatic_atom_mapping (int s) { _do_automatic_atom_mapping = s;}
    void set_remove_explicit_hydrogens (int s) { _remove_explicit_hydrogens = s;}
    void set_unconnect_unmapped_atoms_that_exceed_product_valence (int s) { _unconnect_unmapped_atoms_that_exceed_product_valence = s;}
    void set_convert_A_to_C_for_aromaticity (int s) { _convert_A_to_C_for_aromaticity = s;}
    void set_convert_atom_aliases_to_isotopes (int s) { _convert_atom_aliases_to_isotopes = s;}
    void set_interpret_atom_alias_as_smarts (int s) { _interpret_atom_alias_as_smarts = s;}
    void set_interpret_carbon_on_rhs_of_list_as_no_change (int s) { _interpret_carbon_on_rhs_of_list_as_no_change = s;}
    void set_auto_fix_orphans (int s) { _auto_fix_orphans = s;}
//  void set_is_changing_if_different_neighbours (int s) { _is_changing_if_different_neighbours = s;}
    void set_mol2qry_isotope_special_meaning (int s) { _mol2qry_isotope_special_meaning = s;}
    void set_mark_atoms_changed_when_kekule_form_of_bond_changes (int s) { _mark_atoms_changed_when_kekule_form_of_bond_changes = s;}

    int contains_orphan_atoms() const { return _orphan_atoms.natoms();}
    int check_for_widows_and_orphans ();

    int all_reagents_the_same();   // identify duplicated reagents (drawing errors we have encountered). Will reset _nr to 1 if everything the same

    int remove_duplicate_reagents_ignore_atom_map();
    int remove_duplicate_products_ignore_atom_map();
    int remove_duplicate_agents();
    int remove_all_agents();
    int reduce_to_largest_product();
    int reduce_to_largest_reactant();
    int reduce_to_largest_component(ISIS_RXN_FILE_Molecule *component, int &n);
    
    int eliminate_reagents_not_participating();
    int remove_fragments_not_participating ();
    int remove_non_participating_fragments ();
    int remove_unchanging_fragments();
    int remove_unchanging_components();

    
    
    int move_small_counterions_to_orphan_status();

    int remove_cis_trans_bonding();

    int contains_duplicate_atom_map_numbers() const;
    int unmap_duplicate_atom_map_numbers();
    void assign_unmapped_atoms();

    int max_atom_in_any_reagent() const;
    int max_atom_in_any_product() const;

    int fix_kekule_differences();

    int remove_duplicate_reagents_atom_maps_scrambled ();

    void discard_atom_map ();

    void set_fname_for_echo (const const_IWSubstring & s) { _fname_for_echo = s;}

    int highest_atom_map_number () const;

    int do_read (iwstring_data_source &);
    int do_read (const const_IWSubstring &);

    int build_from_reaction_smiles (const const_IWSubstring & buffer);

    int do_write (std::ostream &) const;
    int do_write (const char * fname) const;

    template <typename T> int write_rxn_smiles(const Reaction_Smiles_Options &, T & output);

    int prepare_for_reaction_construction ();   // allocates arrays, filly reagent/product locator arrays, does atom mapping

    int create_reaction (IWReaction &, const RXN_File_Create_Reaction_Options & rxnfcro, const int * include_these_atoms = NULL);
    int create_reaction (IWReaction &, const RXN_File_Create_Reaction_Options & rxnfcro, Molecule_to_Query_Specifications & mqs, const int * include_these_atoms = NULL);

    int transfer_atom_aliases_to_isotopes ();

    int convert_unchanging_fragments_to_agents ();

    ISIS_RXN_FILE_Molecule & reagent (const int i) { return _reagent[i];}
    ISIS_RXN_FILE_Molecule & product (const int i) { return _product[i];}

    int identify_atoms_changing_reagent (const int r, int *, const Changing_Atom_Conditions & cac);
    int identify_atoms_changing_product (const int p, int *, const Changing_Atom_Conditions & cac);
    int identify_atoms_changing_reagent (const int r, Atom_Typing_Specification & ats, int *, const Changing_Atom_Conditions & cac);
    int identify_atoms_changing_product (const int p, Atom_Typing_Specification & ats, int *, const Changing_Atom_Conditions & cac);

    int contains_isotopic_reagent_atoms() const;
    int at_least_some_mapped_atoms_common_btw_reagents_and_products();

    int largest_fragment_is_changed() const;

    template <typename T> int reaction_fingerprint(const int * changing_atoms,
                               const T * reagent_atom_type,
                               const Product_Atom_Types & pat,
                               const int expand,
                               Sparse_Fingerprint_Creator & sfc);
                               
    int remove_unmapped_components ();

};

/*
  Whereas the reagents will be a single molecule, the products may contain multiple different molecules.
  We need a means of keeping track of their atom types
*/

class Product_Atom_Types
{
  private:
    int * _pstart;
    int * _atom_type;
  public:
    Product_Atom_Types();
    ~Product_Atom_Types();

    int initialise(RXN_File &, Atom_Typing_Specification & ats);

    int atom_type(const int f, const atom_number_t a) const;
};

extern void set_warn_no_mapped_atoms(const int s);

extern int
parse_isis_rxn_file_options (Command_Line & cl,
                             char flag,
                             RXN_File & ISIS_rxn);

extern int
write_isis_reaction_file_header (int nr,
                                 int np,
                                 std::ostream &);

#endif
