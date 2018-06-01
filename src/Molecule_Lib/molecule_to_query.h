#ifndef MOL2QRY_SPEC_H
#define MOL2QRY_SPEC_H

#include "msi_object.h"

#include "mdl_file_data.h"
#include "substructure.h"


/*
  To support the ability to specify substructures at atoms, we need a
  relationship between an element and a smarts
*/

class Element_to_Smarts
{
  private:
    const Element * _e;
    IWString _smarts;

  public:
    Element_to_Smarts ();

    int build (const const_IWSubstring &);

    const Element * element () const { return _e;}
    const IWString & smarts () const { return _smarts;}
};

/*
  We have the ability to create various substructure things from a molecule.
  There are lots of options controlling that process, so rather than long
  argument lists, we have an object that holds all the options

  The atom map stuff is for reading ISIS reaction files

  _highest_atom_number is needed because a temporary array of Substructure_Atom objects
  is used. We need to know the largest number in that array

  Jan 2005. Note there are serious problems with use of this class. It was built
  initially for trxn, but since then, I've used it in a variety of places. The 
  problem is that it while it works well for processing a single molecule, it
  is not really designed for being re-used for multiple molecules - there is no
  dynamic resizing of arrays and such. Just beware how you use this...
  I'll need to fix this!
*/

#define MQS_ISOTOPE_MEANS_NCON 1
#define MQS_ISOTOPE_MEANS_RING_BOND_COUNT 2
#define MQS_ISOTOPE_MEANS_AROMATIC 4

class Molecule_to_Query_Specifications : public MDL_File_Data
{
  private:
    int _make_embedding;
    int _all_ring_bonds_become_undefined;
    int _non_ring_atoms_become_nrings_0;
    int _atoms_conserve_ring_membership;

//  May 2015. Preserve ring membership for ring atoms, but leave others undefined

    int _ring_atoms_conserve_ring_membership;

    int _copy_bond_attributes;

//  Mar 2005. By default, when creating the query, I only insist on aromatic if
//  the atoms are aromatic in the molecule. Chain atoms are matched by atomic number only

    int _only_aromatic_atoms_match_aromatic_atoms;

// Apr 2005.  We need different behaviour when we come from an ISIS
// reaction file than when we come from a molecule.

    int _built_from_isis_reaction_file;

    int _atoms_in_molecule;

// But what if we want to specify that the unsaturation is actually 0. 
// Normally, 0 means not specified.  Use a kludge

#define SPECIAL_MEANING_FULLY_SATURATED -9321

// The _ncon attribute of a single_substructure_query object - number
// of connections to matched atoms

    int _ncon;
    int _max_ncon;
    int _min_ncon;

// Aug 2003.  We want the ability to specify just those atoms that
// allow substitutions, and things like that.  Will be derived from a
// string following the file name. 

    Substructure_Query _substitutions_only_at;

//  Jan 2005. I want to be able to specify an environment around the _substitutions_only_at.

    msi_object _environment_near_substitutions_only_at;
    msi_object _environment_no_match_near_substitutions_only_at;

// The atom numbers where the _substitutions_only_at query matches. 
// note that the values in here will be a function of which molecule
// is being examined

    Set_of_Atoms _substitution_points;

//  When converting R atoms to substitution points on the adjacent atom, I
//  need a means of telling it to use those atoms as attachment points

    Set_of_Atoms _externally_specified_substitution_points;
    
//  We can change all of a certain element to an atomic smarts

    resizable_array_p<Element_to_Smarts> _element_to_smarts;

// Mar 2005.  Want to be able to remove explicit hydrogens and have
// their presence be an attribute of the anchor atoms

    int _condense_explicit_hydrogens_to_anchor_atoms;

// Apr 2005.  Jeff Sutherland wants to be able to match structures
// from crystal.  There, all we really know is atomic numbers and
// whether connected or not

    int _just_atomic_number_and_connectivity;

//  May 2007. VDOM project. I want to be able to NOT set hcount

    int _discern_hcount;

//  Dec 2007. For VDOM. ISIS uses the ring_bond_count designator for
//  most ring stuff. We can change the behaviour of query construction
//  so that the nrings attribute is filled with the ring bond count

    int _nrings_is_ring_bond_count;

//  Mar 2009. Want to be able to specify a minimum number of extra
//  atoms in the target

    int _min_extra_atoms_in_target;
    int _max_extra_atoms_in_target;

    float _min_fraction_atoms_matched;
    float _max_fraction_atoms_matched;

//  Mar 2009. I think I can force matches to one particular atom in
//  symmetry derived queries by using preference values

    int _use_preference_values_to_distinguish_symmetry;

//  Mar 2009. Explicit hydrogens

    int _convert_explicit_hydrogens_to_match_any_atom;

//  Sept 2010. Sometimes we might want to convert all aromatic atoms to just "aromatic"
//  rather than remembering their element type. Optionally, that can be restricted to
//  having the same number of heteroatoms as the starting molecule.

    int _convert_all_aromatic_atoms_to_generic_aromatic;

//  Aug 2013. If we are to create a query that will match a molecule
//  with either explicit or implicit hydrogens, we need to be careful
//  not sure if I ever implemented this...

    int _query_must_match_both_explicit_and_implicit_hydrogens;

//  We might want to enforce only matches to similar saturation

    int _preserve_saturation;

//  Nov 2013. Ignore all Hydrogen information

    int _ignore_molecular_hydrogen_information;

//  Sept 2014. Make interpretation of atom aliases as smarts optional

    int _interpret_atom_alias_as_smarts;

//  Apr 2015 for Jibo

    int _convert_explicit_hydrogens_to_match_any_atom_including_hydrogen;

//  May 2015. Must start getting some of the global variables into this object

    int _substituents_only_at_isotopic_atoms;

//  Mar 2016. It can be convenient to have some per-atom meaning passed in via
//  the isotopic label. Kind of a kludge, but flexible.

    int _isotopic_label_means;

//  Jun 2016. Make setting the element hits needed attribute optional - wastes time in many cases

    int _set_element_hits_needed_during_molecule_to_query;
    int _aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic;

//  Jul 2016. Preserve smallest ring size

    int _preserve_smallest_ring_size;

//  aug 2016. With the existing attributes we cannot accomplish this very simple request

    int _bonds_preserve_ring_membership;

//  private functions

//  template <typename T>
//  int _common_array_copy (const T * a, int n, T * & v);

    int _parse_directive (const_IWSubstring directives);

    int _parse_onlysub_directive (const const_IWSubstring smarts);

  public:
    Molecule_to_Query_Specifications ();
    ~Molecule_to_Query_Specifications ();

    int built_from_isis_reaction_file () const { return _built_from_isis_reaction_file;}
    void set_built_from_isis_reaction_file (int s) { _built_from_isis_reaction_file = s;}

    void set_ncon (int s) { _ncon = s;}
    void set_min_ncon (int s) { _min_ncon = s;}
    void set_max_ncon (int s) { _max_ncon = s;}

    int  ncon () const { return _ncon;}
    int  min_ncon () const { return _min_ncon;}
    int  max_ncon () const { return _max_ncon;}

    int initialise_sub_array_from_only_sub_query (MDL_Molecule &);

    int  just_atomic_number_and_connectivity () const { return _just_atomic_number_and_connectivity;}
    void set_just_atomic_number_and_connectivity (int s) { _just_atomic_number_and_connectivity = s;}

    int  discern_hcount() const { return _discern_hcount;}
    void set_discern_hcount(int s) { _discern_hcount = s;}

    int nrings_is_ring_bond_count () const { return _nrings_is_ring_bond_count;}
    void set_nrings_is_ring_bond_count (int s) { _nrings_is_ring_bond_count = s;}

    int make_embedding () const { return _make_embedding;}
    int all_ring_bonds_become_undefined () const { return _all_ring_bonds_become_undefined;}
    int non_ring_atoms_become_nrings_0 () const { return _non_ring_atoms_become_nrings_0;}
    int atoms_conserve_ring_membership() const { return _atoms_conserve_ring_membership;}
    int ring_atoms_conserve_ring_membership() const { return _ring_atoms_conserve_ring_membership;}
    int copy_bond_attributes () const { return _copy_bond_attributes;}
    int min_extra_atoms_in_target () const { return _min_extra_atoms_in_target;}
    int max_extra_atoms_in_target () const { return _max_extra_atoms_in_target;}
    int use_preference_values_to_distinguish_symmetry () const { return _use_preference_values_to_distinguish_symmetry;}
    int convert_explicit_hydrogens_to_match_any_atom () const { return _convert_explicit_hydrogens_to_match_any_atom;}
    int convert_explicit_hydrogens_to_match_any_atom_including_hydrogen () const { return _convert_explicit_hydrogens_to_match_any_atom_including_hydrogen;}
    float min_fraction_atoms_matched () const { return _min_fraction_atoms_matched;}
    float max_fraction_atoms_matched () const { return _max_fraction_atoms_matched;}
    int convert_all_aromatic_atoms_to_generic_aromatic () const { return _convert_all_aromatic_atoms_to_generic_aromatic;}
    int query_must_match_both_explicit_and_implicit_hydrogens () const { return _query_must_match_both_explicit_and_implicit_hydrogens;}
    int preserve_saturation () const { return _preserve_saturation;}
    int ignore_molecular_hydrogen_information () const { return _ignore_molecular_hydrogen_information;}
    int interpret_atom_alias_as_smarts () const { return _interpret_atom_alias_as_smarts;}
    int substituents_only_at_isotopic_atoms () const { return _substituents_only_at_isotopic_atoms;}
    int isotopic_label_means () const { return _isotopic_label_means;}
    int set_element_hits_needed_during_molecule_to_query () const { return _set_element_hits_needed_during_molecule_to_query;}
    int aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic () const { return _aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic;}
    int preserve_smallest_ring_size() const { return _preserve_smallest_ring_size;}
    int bonds_preserve_ring_membership () const { return _bonds_preserve_ring_membership;}

    void set_make_embedding (int s) { _make_embedding = s;}
    void set_all_ring_bonds_become_undefined (int s) { _all_ring_bonds_become_undefined = s;}
    void set_non_ring_atoms_become_nrings_0 (int s) { _non_ring_atoms_become_nrings_0 = s;}
    void set_atoms_conserve_ring_membership(int s) { _atoms_conserve_ring_membership = s;}
    void set_ring_atoms_conserve_ring_membership(int s) { _ring_atoms_conserve_ring_membership = s;}
    void set_min_extra_atoms_in_target (int s) { _min_extra_atoms_in_target = s;}
    void set_max_extra_atoms_in_target (int s) { _max_extra_atoms_in_target = s;}
    void set_use_preference_values_to_distinguish_symmetry(int s) { _use_preference_values_to_distinguish_symmetry = s;}
    void set_convert_explicit_hydrogens_to_match_any_atom (int s) { _convert_explicit_hydrogens_to_match_any_atom = s;}
    void set_convert_explicit_hydrogens_to_match_any_atom_including_hydrogen (int s) { _convert_explicit_hydrogens_to_match_any_atom_including_hydrogen = s;}
    void set_min_fraction_atoms_matched(float s) { _min_fraction_atoms_matched = s;}
    void set_max_fraction_atoms_matched(float s) { _max_fraction_atoms_matched = s;}
    void set_convert_all_aromatic_atoms_to_generic_aromatic (int s) { _convert_all_aromatic_atoms_to_generic_aromatic = s;}
    void set_query_must_match_both_explicit_and_implicit_hydrogens (int s) { _query_must_match_both_explicit_and_implicit_hydrogens = s;}
    void set_preserve_saturation (int s) { _preserve_saturation = s;}
    void set_ignore_molecular_hydrogen_information (int s) { _ignore_molecular_hydrogen_information = s;}
    void set_interpret_atom_alias_as_smarts (int s) { _interpret_atom_alias_as_smarts = s;}
    void set_substituents_only_at_isotopic_atoms (int s) { _substituents_only_at_isotopic_atoms = s;}
    void set_isotopic_label_means (int s) { _isotopic_label_means = s;}
    void set_set_element_hits_needed_during_molecule_to_query (int s) { _set_element_hits_needed_during_molecule_to_query = s;}
    void set_aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic (int s) { _aromatic_only_matches_aromatic_aliphatic_only_matches_aliphatic = s;}
    void set_preserve_smallest_ring_size(int s) { _preserve_smallest_ring_size = s;}
    void set_bonds_preserve_ring_membership(int s) { _bonds_preserve_ring_membership = s;}

    int parse_directives (const const_IWSubstring &);

    Substructure_Query & substitutions_only_at () { return _substitutions_only_at;}
    const Substructure_Query & substitutions_only_at () const { return _substitutions_only_at;}

    int read_environment_specification (iwstring_data_source & input) { return _environment_near_substitutions_only_at.read (input);}
    int read_environment_no_match_specification (iwstring_data_source & input) { return _environment_no_match_near_substitutions_only_at.read (input);}

    int environment_near_substitution_points_specified () const { return _environment_near_substitutions_only_at.active ();}
    int environment_no_match_near_substitution_points_specified () const { return _environment_no_match_near_substitutions_only_at.active ();}

    const msi_object & environment_near_substitution_points () const { return _environment_near_substitutions_only_at;}
    const msi_object & environment_no_match_near_substitution_points () const { return _environment_no_match_near_substitutions_only_at;}

    int only_aromatic_atoms_match_aromatic_atoms () const { return _only_aromatic_atoms_match_aromatic_atoms;}
    void set_only_aromatic_atoms_match_aromatic_atoms (int s) { _only_aromatic_atoms_match_aromatic_atoms = s;}

    const Set_of_Atoms & substitution_points () const { return _substitution_points;}

    Set_of_Atoms & externally_specified_substitution_points () { return _externally_specified_substitution_points;}

    int set_smarts_for_atom (const const_IWSubstring &);

    int smarts_for_element (const Element *, IWString &) const;

    int condense_explicit_hydrogens_to_anchor_atoms () const { return _condense_explicit_hydrogens_to_anchor_atoms;}
    void set_condense_explicit_hydrogens_to_anchor_atoms (int s) { _condense_explicit_hydrogens_to_anchor_atoms = s;}
};

/*
  When creating query objects from Molecules, sometimes isotopes have
  special meaning
*/

extern void set_substituents_only_at_isotopic_atoms (int);
extern void set_must_have_substituent_at_every_isotopic_atom (int s);
extern void set_isotope_count_means_extra_connections (int s);
extern void set_substitutions_only_at_non_isotopic_atoms (int s);
extern void set_only_include_isotopically_labeled_atoms (int s);
extern void set_only_aromatic_atoms_match_aromatic_atoms (int s);

extern int substituents_only_at_isotopic_atoms();
extern int must_have_substituent_at_every_isotopic_atom();
extern int isotope_count_means_extra_connections();
extern int substitutions_only_at_non_isotopic_atoms();
extern int only_include_isotopically_labeled_atoms();

extern void set_respect_ring_membership(int);

extern void set_molecule_to_query_always_condense_explicit_hydrogens_to_anchor_atoms (int);

#endif
