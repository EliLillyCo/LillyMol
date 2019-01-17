#ifndef IW_MOLECULE_H
#define IW_MOLECULE_H 1

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

class XMLNode;

/*
  Header file for Molecule objects
*/

#include "iwstring.h"
#include "iwmtypes.h"
class iwstring_data_source;

// forward declaration

class Molecule;
class Ring;
class Beep;
class Path_Scoring;
class Smiles_First_Atom;
class Smiles_Formation_Info;
class MDL_File_Supporting_Material;

class CMarkup;    // xml reader used by Marvin

#ifdef COMPILING_MOLECULER_CC
class Rings_Found;
#endif

#ifdef COMPILING_CAREFUL_FRAG
class Fragment_Data;
#endif

#ifdef COMPILING_AROMATIC_CC
class Kekule_Temporary_Arrays;
#endif

class MDL_File_Artifacts_Last_Molecule_Read;

/*
  We need a "smiles" for the empty molecule
*/

#define EMPTY_MOLECULE_SMILES "."

#include "atom.h"

#include "collection_template.h"

class Set_of_Charges : public Collection_Template <charge_t>
{
};

typedef unsigned int atom_type_t;

class Atom_Types : public Collection_Template <atom_type_t>
{
};

/*
  There are lots of choices when we produce a molecular graph
*/

class Mol2Graph
{
  private:
    int _exclude_triple_bonds_from_graph_reduction;
    int _revert_all_directional_bonds_to_non_directional;
    int _preserve_cc_double_bonds_saturated;
    int _preserve_cc_double_bonds_no_heteroatoms;
    int _remove_chiral_centres;
    
    int _append_molecular_formula;

    int _aromatic_distinguishing_formula;

  public:
    Mol2Graph();

    int construct (Command_Line & cl, const char flag, const int verbose);

    int debug_print (std::ostream &) const;

    int exclude_triple_bonds_from_graph_reduction () const { return _exclude_triple_bonds_from_graph_reduction;}
    int revert_all_directional_bonds_to_non_directional() const { return _revert_all_directional_bonds_to_non_directional;}
    int preserve_cc_double_bonds_saturated () const { return _preserve_cc_double_bonds_saturated;}
    int preserve_cc_double_bonds_no_heteroatoms () const { return _preserve_cc_double_bonds_no_heteroatoms;}
    int append_molecular_formula () const { return _append_molecular_formula;}
    int aromatic_distinguishing_formula () const {return _aromatic_distinguishing_formula;}
    int remove_chiral_centres () const { return _remove_chiral_centres;}

    int some_kind_of_double_bond_preservation_active () const { return _preserve_cc_double_bonds_saturated + _preserve_cc_double_bonds_no_heteroatoms;}

    void set_exclude_triple_bonds_from_graph_reduction(int s) { _exclude_triple_bonds_from_graph_reduction = s;}
    void set_revert_all_directional_bonds_to_non_directional(int s) { _revert_all_directional_bonds_to_non_directional = s;}
    void set_preserve_cc_double_bonds_no_heteroatoms (int s) { _preserve_cc_double_bonds_no_heteroatoms = s;}
    void set_preserve_cc_double_bonds_saturated (int s) { _preserve_cc_double_bonds_saturated = s;}
    void set_append_molecular_formula (int s) { _append_molecular_formula = s;}
    void set_aromatic_distinguishing_formula(int s) { _aromatic_distinguishing_formula = s;}
    void set_remove_chiral_centres (int s) { _remove_chiral_centres = s;}
};

/*
  We sometimes want to be able to preserve various atom attributes
  from a pdb file
*/

class PDB_Stored_Atom_Information
{
  private:
    int _atom_number;
    IWString _residue_name;
    IWString _atom_name;
    IWString _occupancy_factor;
    IWString _temperature_factor;

//  private functions

    void _do_copy (const PDB_Stored_Atom_Information & rhs);

  public:
    PDB_Stored_Atom_Information();
    PDB_Stored_Atom_Information(const PDB_Stored_Atom_Information &);

    PDB_Stored_Atom_Information & operator = (const PDB_Stored_Atom_Information & rhs);

    void set_atom_number (int s) { _atom_number = s;}
    void set_residue_name (const IWString & s) { _residue_name = s;}
    void set_residue_name (const char * s) { _residue_name = s;}
    void set_atom_name (const IWString & s) { _atom_name = s;}
    void set_occupancy (const IWString & s) { _occupancy_factor = s;}
    void set_temperature (const IWString & s) { _temperature_factor = s;}

    int atom_number () const { return _atom_number;}
    const IWString & atom_name () const { return _atom_name;}
    const IWString & residue_name () const { return _residue_name;}
    const IWString & occupancy_factor () const { return _occupancy_factor;}
    const IWString & temperature_factor () const { return _temperature_factor;}
};


#include "bond_list.h"

/*
  The charge array is initialised to NULL. As a molecule
  is built, if we encounter a non-zero charge, we allocate
  an array for it, and fill all members with 0.0 (except
  for the non zero charge just encountered).
*/

/*
  The flag _partially_built is set whenever add_bond is called, indicating
  a partially built molecule. Typically, this means that bonds have been
  added for which the atoms may not yet be in the molecule.
*/

#include "coordinates.h"

#define MOLECULE_MAGIC_NUMBER -7215237

#define IW_NRINGS_NOT_COMPUTED -87654
#define IW_RING_MEMBERSHIP_NOT_COMPUTED -41871
#define IW_RING_MEMBERSHIP_IS_A_RING_ATOM   -76

#define FRAGMENT_MEMBERSHIP_NOT_SET -1

/*
  The implicit hydrogen count is somewhat complex, because most times, it is
  a derived quantity. Sometimes (as with reading an aromatic smiles for example)
  there is a known, immutable value. In those cases, we set the HCOUNT_STICKY_BIT
  to indicate that the value should not be changed blithely
*/

#define IMPLICIT_HCOUNT_STICKY_BIT 0x10000000

/*
  Make sure that this value does not collide with the implicit_hydrogen_sticky_bit,
*/

#define IMPLICIT_HYDROGENS_NOT_COMPUTED 0x01000000

// Values for _smiles_order_type

#define INVALID_SMILES_ORDER_TYPE 0
#define UNIQUE_SMILES_ORDER_TYPE 2
#define RANDOM_SMILES_ORDER_TYPE 3
#define DEFAULT_SMILES_ORDER_TYPE 4
#define SUBSET_SMILES_ORDER_TYPE 5
#define USER_SPECIFIED_SMILES_ORDER 6


#define REASONABLE_RING_SIZE(r) ((r) > 2)

class IW_Bits_Base;
class Path;
class List_of_Ring_Sizes : public resizable_array<int>
{
};

class Fragment_Information
{
  protected:
    int _number_fragments;

    int * _fragment_membership;

    resizable_array<int> _atoms_in_fragment;

    resizable_array<int> _bonds_in_fragment;

  public:
    Fragment_Information ();
    ~Fragment_Information ();

    int debug_print (std::ostream &) const;

    int  contains_valid_data () const;

    int initialise (int);

    void invalidate ();

    int number_fragments () const {return _number_fragments;}

    int set_number_fragments (int);

    int * fragment_membership () { return _fragment_membership;}
    const int * fragment_membership () const { return _fragment_membership;}

    int fragment_membership (int a) const { return _fragment_membership[a];}

    resizable_array<int> & atoms_in_fragment () {return _atoms_in_fragment;}
    resizable_array<int> & bonds_in_fragment () {return _bonds_in_fragment;}

    int atoms_in_fragment (int f) const { return _atoms_in_fragment[f];}
    int bonds_in_fragment (int f) const { return _bonds_in_fragment[f];}

    const resizable_array<int> & atoms_in_fragment () const {return _atoms_in_fragment;}
    const resizable_array<int> & bonds_in_fragment () const {return _bonds_in_fragment;}

    int atoms_in_fragment (int matoms, int f, Set_of_Atoms &) const;

    int rings_in_fragment (int f) const { return _bonds_in_fragment[f] - _atoms_in_fragment[f] + 1;}

    int all_atoms_in_one_fragment (int natoms, int nbonds);   // initialise the structure
};

#include "iwrcb.h"

/*
  Once fragments, rings and aromaticity are known, we can generate a smiles ordering

  Apr 2004. Need to make sure that

  m.smarts ()
  m.smiles ()

  forces recomputation of the smiles
*/

class Smiles_Information
{
  protected:
    int _natoms;

    int * _smiles_order;

    int _smiles_order_type;

//  Save the atom numbers where the smiles for each fragment starts.
//  This is for efficiency.

    Set_of_Atoms _smiles_start_atom;

//  During smiles determination, we keep track of any ring closing bonds found

    Ring_Closure_Bonds _ring_closure_bonds;

// Whereas _smiles_order is a guide for making decisions at each atom
// about how to assemble the smiles, sometimes we want to keep track
// of the order that each atom shows up in the smiles.

    resizable_array<atom_number_t> _atom_order_in_smiles;

    IWString _smiles;          

    int _smiles_is_smarts;

//  If we are creating a smarts, we can specify, on a per atom basis, whether each
//  atom is represented as it appears in the molecule, or as an embedding.

    int * _create_smarts_embedding;

//  We can also specify any arbitrary string to be the smarts for any atom

    IWString * _user_specified_atomic_smarts;

//  private functions

    void _default_values ();

  public:
    Smiles_Information ();
    Smiles_Information (int);
    ~Smiles_Information ();

    int debug_print (std::ostream &) const;

    void invalidate ();

    int prepare_to_build_ordering (int matoms);

    int prepare_to_build_smiles (int matoms);

    const IWString & smiles () const { return _smiles;}
    IWString & smiles () { return _smiles;}

    int allocate_user_specified_atomic_smarts();
    IWString & user_specified_atomic_smarts(atom_number_t);
    const IWString & user_specified_atomic_smarts(atom_number_t) const;
    void set_user_specified_atomic_smarts(atom_number_t, const IWString &);
    const IWString * user_specified_atomic_smarts() const { return _user_specified_atomic_smarts;}

    int contains_smiles () const { return _smiles.length ();}
    int contains_valid_ordering () const { return INVALID_SMILES_ORDER_TYPE != _smiles_order_type;}

    int smiles_is_smarts () const { return _smiles_is_smarts;}
    void set_smiles_is_smarts (int s) { _smiles_is_smarts = s;}

    void make_empty ();

    const IWString & set_error () { _smiles = "ERROR"; return _smiles;}

    void add_atom (atom_number_t a) { _atom_order_in_smiles.add (a);}

    int * smiles_order () {return _smiles_order;}
    const int * smiles_order () const {return _smiles_order;}

    int smiles_order_type () const { return _smiles_order_type;}
    void set_smiles_order_type (int s) { _smiles_order_type = s;}

    int contains_ring_closure_bond (atom_number_t a1, atom_number_t a2) const { return _ring_closure_bonds.contains (a1, a2);}
    int add_ring_closure_bond (atom_number_t a1, atom_number_t a2) { return _ring_closure_bonds.add (a1, a2);}

    int add_start_atom (atom_number_t a) { return _smiles_start_atom.add (a);}
    const Set_of_Atoms & smiles_start_atom () const { return _smiles_start_atom;}

    const resizable_array<atom_number_t> & atom_order_in_smiles () const { return _atom_order_in_smiles;}

    const int * create_smarts_embedding () const { return _create_smarts_embedding;}
    int create_smarts_embedding (atom_number_t) const;

//  We have the ability to set embedding mode on or off for the whole molecule,
//  or for individual atoms

    int set_create_smarts_embedding (int);
    int set_create_smarts_embedding (atom_number_t, int);
};

/*
  We break these out as a separate class because of subsets.
*/

#define IW_SYMMETRY_CLASS_UNDEFINED -1

class Symmetry_Class_and_Canonical_Rank
{
  private:
    int * _symmetry_class;

    int * _canonical_rank;

  public:
    Symmetry_Class_and_Canonical_Rank ();
    ~Symmetry_Class_and_Canonical_Rank ();

    int debug_print (std::ostream &) const;

    int invalidate ();

    int allocate_arrays (int);

    int arrays_allocated () const { return NULL != _symmetry_class;}

    const int * symmetry_class () const { return _symmetry_class;}
    const int * canonical_rank () const { return _canonical_rank;}

    int * symmetry_class () { return _symmetry_class;}
    int * canonical_rank () { return _canonical_rank;}

    int symmetry_class (int a) { return _symmetry_class[a];}
    int canonical_rank (int a) { return _canonical_rank[a];}

    int store_values_from (const Symmetry_Class_and_Canonical_Rank &, int);   // like operator=, but it needs to know the number of atoms
};

class Smiles_Ring_Status;
class Tnode;
class Ring_Number_Manager;

class Chiral_Centre;

#include "set_of_atoms.h"

/*
  We want various degrees of control over how implicit Hydrogens are added
*/

class Make_Implicit_Hydrogens_Explicit
{
  private:
    atom_number_t _a;

    int _isotope;
    int _dimensionality;

  public:
    Make_Implicit_Hydrogens_Explicit ();

    void reset () { _a = INVALID_ATOM_NUMBER;}

    int isotope () const { return _isotope;}
    int dimensionality () const { return _dimensionality;}

    void set_isotope (int s) { _isotope = s;}
    void set_dimensionality (int s) { _dimensionality = s;}

    atom_number_t a () const { return _a;}
    void set_atom (atom_number_t s) { _a = s;}

    Atom * new_atom () const;
};

/*
  When computing molecular weights in the thread safe version, we use this class
  to control the calculation
*/

class Molecular_Weight_Control
{
  public:
    bool _ignore_isotopes;
    bool _ignore_non_periodic_table_elements;
    bool _ignore_hydrogens;

  public:
    Molecular_Weight_Control();

    void set_ignore_isotopes(int s) { _ignore_isotopes = s;}
};

class Molecular_Weight_Calculation_Result
{
  public:
    int _isotopes_found;
    int _non_periodic_table_elements_found;
    double _amw;

  private:
    void _default_values();

  public:
    Molecular_Weight_Calculation_Result();

    void reset();

    double amw() const { return _amw;}
};

class const_BondIterator
{
  private:
    const Atom * const * _atom;
    const atom_number_t _atnum;
    const Bond * const * _b;

  public:
    const_BondIterator(const Atom *, const atom_number_t);

    int operator!=(const Bond * const * b) const { return b != _b;}

    void operator++() { _b++;}

    atom_number_t operator*() const { return (*_b)->other(_atnum);}
};

class Molecule : private resizable_array_p<Atom>
{
  private:

    Bond_list _bond_list;

    resizable_array_p<Chiral_Centre> _chiral_centres;

//  We may also have information about the fragment characteristics of the molecule.
//  Only computed if needed

    Fragment_Information _fragment_information;

    int       _nrings;
    int       _number_sssr_rings;
    int *     _ring_membership;       // within the SSSR set
    resizable_array_p<Ring> _raw_rings;
    resizable_array_p<Ring> _sssr_rings;

//  During uniqueness determination, we need to know true ring membership,
//  rather than the SSSR ring membership. For example, in cubane, we want
//  to know that each atom is in three four membered rings, whereas the
//  SSSR ring set has four atoms in three rings, and four in just two rings.

//  We keep any rings which were rejected by the SSSR process.

    resizable_array_p<Ring> _non_sssr_rings;

    resizable_array_p<Ring> _experimental_raw_rings;
    resizable_array_p<Ring> _experimental_sssr_rings;

//  We store aromaticity for the atoms. Aromaticity of rings
//  is stored with the Ring objects

    aromaticity_type_t * _aromaticity;
    
    bool doNotComputeAromaticity;   // dangerous to use - this prevent the aromatcity routines from functioning - needed for usmiles of reaction cores.

    Smiles_Information _smiles_information;

    Symmetry_Class_and_Canonical_Rank _symmetry_class_and_canonical_rank;

//  Distance_Matrix * _distmat;
    int * _distance_matrix;

//  Fractional atomic charges

    Set_of_Charges * _charges;    // not yet implemented, using old way

//  The set of atom types can be used to hold any information

    Atom_Types * _atom_type;

//  When reading an MDL MOLFILE or a TDT file we can optionally bring
//  along the text info from the file

    resizable_array_p<IWString> _text_info;

    void * _user_specified_void_ptr;

    IWString _molecule_name;
    int _partially_built;
    magic_number_t _magic;

//  Private functions

    void _resize (int);

    int  _set_modified (atom_number_t);
    int  _set_modified ();
    int _set_modified_no_ok ();
    int  _remove_bonds_to_atom (atom_number_t, int = 0);    // optional arg is whether or not to renumber things for loss of the atom
    int  _free_all_dynamically_allocated_things ();

    void _add_bond (Bond *, int);

    int _build_smiles ();

    void _default_values (int);

//  Lots of functions for dealing with MDL files

#ifdef COMPILING_MDL_CC
#include "molecule_mdl.h"
#endif

#ifdef COMPILING_CAREFUL_FRAG
#include "careful_frag.h"
#endif

#ifdef COMPILING_MOLECULE_SMARTS
#include "molecule_smarts.h"        
#endif

#ifdef COMPILING_MOLECULE_CIF
#include "molecule_cif.h"        
#endif

//  Stuff needed for finding kekule forms:

    int  _find_kekule_form (resizable_array<Bond *> &);
    int  _find_kekule_form (resizable_array<Bond *> &, int);
    int  _identify_possible_aromatic_rings (resizable_array<Bond *> & , int *);

    int _do_unconnect_covalently_bonded_non_organics ();

//  private because it could be used to build partially correct cis-trans groupings

    int _set_bond_directionality (atom_number_t, atom_number_t, int);

#ifdef COMPILING_MOLECULER_CC
#include "tmpsssr.h"
#endif

#ifdef COMPILING_MOLECULE_MAIN
#include "molecule_main.h"
#endif

#ifdef COMPILING_MOLECULE_B
#include "moleculeb.h"
#endif

#ifdef COMPILING_SMILES_CC
#include "molecule_smi.h"
#endif

//  New distance matrix stuff

#ifdef COMPILING_MOLECULED
#include "moleculed.h"
#endif

#ifdef COMPILING_TRIPOS_CC
#include "molecule_tripos.h"
#endif

#ifdef COMPILING_MARVIN_CC
#include "molecule_marvin.h"
#endif

    int  _build_from_smiles     (const char *, int, Smiles_Ring_Status &, int *);
    int  _build_from_smiles     (const char *, int, int *);
    int  _build_from_smiles     (const char *, int);

    int _read_mrk_atom_record (const const_IWSubstring & buffer);
    int _read_mrk_bond_record (const const_IWSubstring & buffer, int na);

    int _final_processing_of_aromatic_mdl_input(int * aromatic_atoms, int * aromatic_bonds);

//  By convention, molecule names are standardised, no leading spaces, etc..

    void _standardise_name ();

//  If we get something like C[CH2]C we don't need the implicit hydrogens
//  known attribute on the central carbon atom.

    int _unset_implicit_hydrogens_known_if_computed_matches ();

    int _unset_all_implicit_hydrogens_known_attributes ();

    int _compute_implicit_hydrogens (atom_number_t, int &) const;
    int _compute_and_store_implicit_hydrogens (atom_number_t);
    int _compute_implicit_hydrogens (atom_number_t);
    int _compute_implicit_hydrogens ();

//  Things associated with rings.

    void _initialise_ring_membership ();
    void _determine_ring_or_non_ring (atom_number_t a);

    int  _compute_aromaticity ();

#ifdef COMPILING_AROMATIC_CC
#include "molecule_arom.h"
#endif

//  Stuff related to symmetry determinations

    void _new_symmetry_class (atom_number_t);
    void _add_symmetry_class (const resizable_array<atom_number_t> &);

    void _ensure_symmetry_perceived ();

    int _find_raw_rings (const atom_number_t previous_atom,
                         const atom_number_t current_atom,
                         resizable_array<Ring *> & rings,
                         resizable_array<atom_number_t> & active_rings,
                         int * already_done);
    int _find_raw_rings_for_fragment (int id, int * already_done);
    int _find_raw_rings (int * already_done);
    int _find_raw_rings ();

//  Functions dealing with chiral centres

    int  _complete_chiral_centre_from_mdl_files (Chiral_Centre *, const MDL_File_Supporting_Material &);
    int  _write_mdl_atom_stereo_info (std::ostream & os, atom_number_t a) const;

    int _remove_directionality_from_bonds_not_actually_directional ();
    int _remove_directionality_from_bonds_not_actually_directional (atom_number_t);

    int _discern_chirality_from_3d_structure(atom_number_t zatom);

    void _print_atom_and_type (std::ostream & os, const char * s, atom_number_t centre, atom_number_t a) const;   // only used in chiral_centre code

    int _stereo_centre_hydrogens_become_implicit (Chiral_Centre * c);

    int _add_chiral_centre_checking_for_duplicate (Chiral_Centre * c);

    int _check_chirality_after_loss_of_bond (const atom_number_t a1, const atom_number_t a2);

    int  _write_molecule_tdt_pcn (std::ostream & os, const IWString &) const;

    int _write_molecule_mol2 (std::ostream & os, const int * atype);

    int  _adjust_chiral_centres_for_loss_of_atom (atom_number_t a, int = 0);

    int  _smi_process_new_chiral_centre (Chiral_Centre * c, int hcount) const;
    int  _smi_last_atom_is_part_of_chiral_centre (atom_number_t previous_atom,
                                                 int previous_atom_chiral_count);
    int _smi_atom_bonded_to_chiral_centre (atom_number_t previous_atom,
                                           int previous_atom_chiral_count,
                                           atom_number_t atom_bonded_to_chiral_centre);
    int  _check_for_incomplete_chiral_specifications (Chiral_Centre * c);
    int  _check_for_incomplete_chiral_specifications ();
    int  _check_chiral_centres () const;

    int _score_ez_bond (int * already_done, Path_Scoring & ps1, Path_Scoring & ps2) const;

    void _set_partial_charge_type (const const_IWSubstring &);

    int  _renumber_atoms (const int *, int *);

    int _bump_check (atom_number_t, atom_number_t, atom_number_t, atom_number_t, distance_t too_close, int * either_side) const;

    int _determine_moving_atoms(atom_number_t zatom, int * moving_atoms) const;
    int _determine_either_side_of_bond (atom_number_t a1, atom_number_t a2, int * either_side) const;
    int __identify_side_of_bond (int * either_side, atom_number_t astart, int flag, atom_number_t avoid) const;

    int _do_put_formal_charges_on_neutral_ND3v4 ();

  protected:

    int  _mdl_set_bond_directionality (atom_number_t, atom_number_t, int);
    int  _complete_chiral_centres_from_mdl_files (const MDL_File_Supporting_Material &);
    int _parse_v30_bond_record (const const_IWSubstring & buffer, int * aromatic_atom, int & aromatic_bond, int = 0);
    int _fill_empty_molecule_with_null_atoms (int na);

#include "moleculej.h"

  public:
    Molecule (int = 0);
    Molecule (const Molecule &);
    ~Molecule ();

    Molecule & operator = (const Molecule &);
    Molecule & operator = (Molecule &&);
    bool       operator == ( Molecule & );

    int  add (Atom *, int = 0);   // optional arg means partial molecule, no invalidation
    int  add (const Element *);
    int  resize (int);

    int  ok () const;          // quick audit function
    int  debug_print (std::ostream &) const;

    int  check_bonding   () const;     // detailed audit function
    int  check_ring_info () const;     // checks rings for compatibility w/ this
    int  check_chemistry () const;     // detailed audit function

    int  print_ring_info (std::ostream &) const;

    int ok_atom_number (atom_number_t) const;
    int ok_2_atoms   (atom_number_t, atom_number_t) const;
    int ok_3_atoms   (atom_number_t, atom_number_t, atom_number_t) const;
    int ok_4_atoms   (atom_number_t, atom_number_t, atom_number_t, atom_number_t) const;

    int invalidate_from_possibly_invalid_state();

    coord_t   x (atom_number_t) const;        // returns the x coord
    coord_t   y (atom_number_t) const;        // returns the y coord
    coord_t   z (atom_number_t) const;        // returns the z coord

    int delete_all_atoms_and_bonds ();  // remove all contents from a molecule.

    int   add_bond (atom_number_t, atom_number_t, bond_type_t, int = 0);
    int   are_bonded (atom_number_t, atom_number_t) const;
    int   are_bonded (atom_number_t, atom_number_t, bond_type_t &) const;
    int   finished_bond_addition ();    // call after add_bond with partial mol
    int   assign_bond_numbers_to_bonds ();
    int   assign_bond_numbers_to_bonds_if_needed ();
    int   which_bond (atom_number_t, atom_number_t) const;

//  Which bonds are set by the atoms in a Set_of_Atoms

    int   convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s, int * barray);
    int   convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s, int * barray) const;

//  In reactions, we want a substitution to be able to preserve stereochemistry

    int stereo_preserving_substitute (atom_number_t, atom_number_t, atom_number_t);
    int stereo_preserving_substitute (atom_number_t, atom_number_t);

    int has_formal_charges () const;
    int has_no_formal_charges () const;     // could be ! has_formal_charges ()
    int number_formally_charged_atoms () const;
    int net_formal_charge() const;

    int natoms () const;
    int natoms (atomic_number_t) const;     // number atoms with this atomic_number
    int natoms (const Element *) const;     // number atoms of this element type
    int natoms (const char *) const;        // number atoms of this atomic symbol

    int number_different_elements () const;

    const Atom * atomi  (atom_number_t) const;                 // the i'th atom
    const Atom & atom   (atom_number_t) const;
    atom_number_t which_atom (const Atom * a) const { return index((Atom *) a);}
    int    atoms  (const Atom **) const;                 // dangerous, use sparingly
    int    ncon   (atom_number_t) const;                 // number connections to i'th atom
    int    nbonds (atom_number_t) const;                 // number of bonds to i'th atom
    atom_number_t other (atom_number_t, int) const;       // atom # of J'th connection to atom I

    bond_type_t   btype_to_connection (atom_number_t, int) const;       // btype of J'th connection to atom I
    bond_type_t   btype_between_atoms (atom_number_t, atom_number_t) const;  // btype between two atom numbers

    int   other_and_type (atom_number_t, int, atom_number_t &, bond_type_t &) const;

    int   connections (atom_number_t, atom_number_t *, bond_type_t * = NULL) const;
    int   connections (atom_number_t, Set_of_Atoms &) const;

    int   connections_and_types (atom_number_t,
                Set_of_Atoms &,
                resizable_array<bond_type_t> &) const;
    int   connections_and_types (atom_number_t, atom_number_t *, bond_type_t *) const;

    int   bond_types (atom_number_t, resizable_array<bond_type_t> &) const;
    int   bond_types (atom_number_t, bond_type_t *) const;

    int all_atoms_connected (atom_number_t, Set_of_Atoms &) const;

    template <typename F> void each_atom(F &) const;
    template <typename F> void each_atom_lambda(F) const;
    template <typename F> void each_bond(F &) const;
    template <typename F> void each_ring(F &) const;

    int maximum_connectivity () const;

    int   set_bond_type_between_atoms (atom_number_t a1, atom_number_t a2,
                                       bond_type_t bt);
    int set_bond_type_between_atoms_no_check_btype (atom_number_t a1, atom_number_t a2,
                                       bond_type_t bt);

    int   set_wedge_bond_between_atoms (atom_number_t a1, atom_number_t a2, int);
    int   number_up_or_down_wedge_bonds () const;

    const Element * elementi (atom_number_t) const;     // the element of the I'th atom
    const Element & element  (atom_number_t) const;     // the element of the I'th atom

    int   set_element (atom_number_t a, const Element * e);
    int   set_atomic_number (atom_number_t a, atomic_number_t);

    int   number_isotopic_atoms () const;
    int   number_isotopic_atoms (int) const;
    int   transform_to_non_isotopic_form ();
    int   set_isotope  (atom_number_t, int);
    int   set_isotope_no_perturb_canonical_ordering (atom_number_t a, int iso);
    int   set_userAtomType(atom_number_t a, int atomType);
    template <typename T> int   set_isotopes (const T *);     // changed 2016. The set the isotope to the value in the array
    int   unset_isotopes (const int *);   // anything > 0 will be set to 0
    int   unset_isotopes () { return transform_to_non_isotopic_form ();}
    int   set_isotope (const Set_of_Atoms &, int);
    void  get_isotopes (int *) const;
    int   isotope (atom_number_t) const;
    int   userAtomType (atom_number_t) const;
    
    int   increment_isotope (atom_number_t, int);
    int   maximum_isotope() const;
    atom_number_t atom_with_isotope(int) const;
    void  set_isotope_to_atom_number_no_perturb_canonical_ordering();

    const IWString & atomic_symbol (atom_number_t) const;    // atomic symbol of I'th atom
    atomic_number_t atomic_number (atom_number_t) const; // atomic number of I'th atom
    void atomic_numbers (atomic_number_t *) const;     // all atomic numbers
    int  ncon (int *) const;                           // all ncon values
    int  ncon (resizable_array<int> &) const;
    int  nbonds (int *) const;                         // all nbonds values

//  Once we introduce ESSSR rings, the nrings() method is ambiguous. If the
//  ESSSR rings have been determined, it will report the number of ESSSR rings,
//  but if ring determination hasn't been done, it will report the number
//  of SSSR rings computed by Euler's formula.
//  
    int nrings ();                  // number of rings in molecule
    int nrings_no_compute() const;    // number of rings in molecule without performing ring perception (assumes has already been performed)
    int number_sssr_rings ();

//  If the 2nd argument is specified, we report only rings of that size.
//  If the first argument is specified, we report data for that atom.

    int nrings (atom_number_t);     // number of rings for atom I
    int nrings (atom_number_t, int);  // rings for atom I of a given size
    int nrings_size (int);          // number of rings of size N

    int ring_bond_count (atom_number_t);

//  Given a set of atoms with uncertain bonding, can a valid Kekule representation be found
//  First argument is set for each aromatic atom. The optional 2nd arg is for aromatic bonds

    int find_kekule_form (int *, const int * = NULL);

    int generate_switched_kekule_forms (resizable_array_p<Molecule> & variant);   // not exhaustive

//  Called by mdl and tripos reading functions. Sometimes carboxyllic acids, nitros, sulf* and phosp* acids
//  come in with aromatic bonds

    int process_delocalised_carbonyl_bonds (int * aromatic_atoms, int * aromatic_bonds = NULL);

    int is_part_of_fused_ring_system (atom_number_t a);
    int fused_system_size (atom_number_t a);
    int fused_system_size_no_compute (atom_number_t a) const;
    int rings_with_fused_system_identifier (int);
    int fused_system_identifier (atom_number_t);

//  in_same_ring returns 1 if the two atoms are in the same ring.
//  in_same_rings returns a count of the number of rings which contain the two atoms

    int in_same_ring_no_compute (const atom_number_t &, const atom_number_t & ) const;
    int in_same_ring          (atom_number_t, atom_number_t);
    int in_same_rings         (atom_number_t, atom_number_t);
    int in_same_aromatic_ring (atom_number_t, atom_number_t);

    int in_same_ring_system   (atom_number_t, atom_number_t);

    const int * ring_membership ();
    int ring_membership (int *);

    int is_ring_atom (atom_number_t);
    int is_non_ring_atom (atom_number_t);

    int ring_or_non_ring (int *);

//  The ring sizes of the rings for a given atom;

    int ring_sizes_for_atom (atom_number_t, List_of_Ring_Sizes &);
    int ring_sizes_for_all_atoms (resizable_array_p<List_of_Ring_Sizes> &);
    int in_ring_of_given_size (atom_number_t a, int ring_size);

    const Ring * ringi_no_compute (int i) const;    // the I'th ring  (assuming rings already computed)
    const Ring * ringi (int);        // the I'th ring
    const Ring * ringi (int, int);   // the I'th ring of size J
    const Ring * ring_containing_atom (atom_number_t);  // first ring containing atom I

    int count_heteroatoms (const Set_of_Atoms & r) const;

    int ok_ring (const Ring *) const;

    int get_fused_system (int, Set_of_Atoms &);

    int label_atoms_by_ring_system_no_compute (int *) const;
    int label_atoms_by_ring_system (int *);
    int label_atoms_by_ring_system_including_spiro_fused(int * r);

//  Functions which deal with the non sssr rings

    int nrings_including_non_sssr_rings (atom_number_t);
    int ring_membership_including_non_sssr_rings (int *);
    int ring_sizes_for_non_sssr_rings (atom_number_t, List_of_Ring_Sizes &, int = 0);
    int non_sssr_rings ();
    int non_sssr_rings_no_compute() const;
    const Ring * non_sssr_ring (int);
    const Ring * non_sssr_ring_no_compute (int) const;

    int rings_with_strongly_fused_ring_neighbours();

    int is_spiro_fused (atom_number_t);

    int is_halogen (atom_number_t) const;

    int nedges () const;

    const IWString & smiles ();

//  If you want a subset, you need to pass a Smiles_Information object to hold the subset into. It
//  is never stored in the molecule itself

    const IWString & smiles (Smiles_Information &, const int *);
    const IWString & unique_smiles ();
    const IWString & unique_smiles (Smiles_Information &, const int *);
    const IWString & non_aromatic_unique_smiles ();   // use only in special circumstances
    const IWString & random_smiles ();
    const IWString & smiles_starting_with_atom (atom_number_t);
    const IWString & smiles_starting_with_atom (atom_number_t, Smiles_Information &, const int *);
    const IWString & smiles_using_order (const int * zorder);

    const Smiles_Information & smiles_information () const { return _smiles_information;}

    IWString         isotopically_labelled_smiles();

    int  change_to_graph_form ();
    int  change_to_graph_form (const Mol2Graph &);
    int  set_all_bonds_to_type (bond_type_t);

    int  invalidate_smiles ();
    int  invalidate_fragment_membership ();
    int  invalidate_canonical_ordering_information ();

    int  smiles_atom_order (int *);

    IWString molecular_formula ();
    int  molecular_formula (IWString &) const;
    int  isis_like_molecular_formula_dot_between_fragments (IWString &);
    int  isis_like_molecular_formula (IWString &);    // does all components present
    int  formula_distinguishing_aromatic (IWString &);   // C8H2c6H3O2H
    void set_name (const char *);
    void set_name (const char *, int);
    void set_name (const IWString &);

    void append_to_name (const IWString &);

    const char * molecule_name () const;
    const IWString & name () const;

    int write_molecule       (std::ostream &, int, const IWString & = "");

    int write_molecule_smi   (std::ostream &, const IWString &);

    int write_molecule_tdt   (std::ostream &, const IWString &);
    int write_molecule_tdt_unique (std::ostream &, const IWString &);
    int write_molecule_tdt_nausmi (std::ostream & os, const IWString & comment);

    int write_molecule_usmi    (std::ostream &, const IWString &);
    int write_molecule_nausmi  (std::ostream &, const IWString &);
    int write_molecule_rsmi    (std::ostream &, const IWString &);

//  template <typename T> int write_molecule_mdl     (const char *, const char *) const;
    int write_molecule_mdl (const char *, const char *) const;
    template <typename T> int write_molecule_mdl     (T &, const IWString &) const;
    template <typename T> int write_molecule_mdl_v30 (T &, const IWString &, int) const;

    int write_molecule_pdb   (const char *, const IWString &);
    int write_molecule_pdb   (std::ostream &,    const IWString &);

    int write_molecule_mmod  (std::ostream &) const;

    int write_molecule_msi   (std::ostream &, const IWString & = "") const;

    int write_molecule_bfile (std::ostream &);

    int write_molecule_mol2  (std::ostream &);

    int write_molecule_crd   (std::ostream &);
    int write_molecule_psf   (std::ostream &);

    int write_molecule_mrk   (std::ostream &);
    int write_molecule_wchm  (std::ostream &);

    int write_molecule_cif (std::ostream &);

    int write_molecule_smarts (std::ostream &);

    int write_molecule_mrv (std::ostream &);
    int write_molecule_inchi (std::ostream &);

    int read_molecule_mrv_molecule (XMLNode & cml);

    template <typename T> int write_connection_table_mdl (T &) const;
    int write_connection_table_pdb (std::ostream &);

    int write_set_of_bonds_as_mdl_v30_collection (const resizable_array<int> & b,
                                                    const const_IWSubstring & zname,
                                                    const const_IWSubstring & subname,
                                                    std::ostream & output) const;
    int write_set_of_bonds_as_mdl_v30_collection (const int * b,
                                                    const const_IWSubstring & zname,
                                                    const const_IWSubstring & subname,
                                                    std::ostream & output) const;

    int read_molecule_ds       (iwstring_data_source &, int);
    int read_molecule_pdb_ds   (iwstring_data_source &);
    template <typename T> int read_molecule_mdl_ds   (T &, int = 0);
    int read_molecule_rdf_ds   (iwstring_data_source &);
    int read_molecule_smi_ds   (iwstring_data_source &);
    int read_molecule_tdt_ds   (iwstring_data_source &);
    int read_molecule_mmod_ds  (iwstring_data_source &);
    int read_molecule_msi_ds   (iwstring_data_source &);
    int read_molecule_mol2_ds  (iwstring_data_source &);
    int read_molecule_mrk_ds   (iwstring_data_source &);
    int read_molecule_mrv_ds   (iwstring_data_source &);
    int read_molecule_inchi_ds (iwstring_data_source &);
    int read_molecule_cif_ds   (iwstring_data_source &);

    int build_from_smiles     (const char *);
    int build_from_smiles     (const char *, int);
    int build_from_smiles     (const IWString &);
    int build_from_smiles     (const const_IWSubstring &);

    int build_from_inchi      (const const_IWSubstring &);

    int InChI                 (IWString &);

    const Bond * bondi (int) const;                // pointer to I'th bond in molecule
    const Bond * bondi (atom_number_t, int) const;   // pointer to J'th bond to atom I
    const Bond * bond_between_atoms (atom_number_t, atom_number_t) const;    // pointer to bond between two atoms
    const Bond * bond_between_atoms_if_present(const atom_number_t a1, const atom_number_t a2) const;   // will not complain if they are not bonded
    bool  bond_endpoints (int ndx, atom_number_t & a1, atom_number_t & a2) const;   // bond endpoint indices

    const Bond_list & bond_list() const {return _bond_list;}   // dangerous, but we return it const

    int compute_canonical_ranking (Symmetry_Class_and_Canonical_Rank &, const int *);
    int compute_canonical_ranking ();
    int canonical_rank (atom_number_t);
    int canonical_ranks (int *);
    const int * canonical_ranks ();

    const resizable_array<atom_number_t> & atom_order_in_smiles () const { return _smiles_information.atom_order_in_smiles ();}

    int symmetry_class (atom_number_t);
    int number_symmetry_classes ();
    int symmetry_equivalents (atom_number_t, Set_of_Atoms &);

    const int * symmetry_classes ();

    int bond_symmetry_class_large_memory (int * s);
    int bond_symmetry_class_small_memory (int * s);

    int number_hydrogens () const;

    int attached_heteroatom_count   (atom_number_t) const;
    int multiple_bond_to_heteroatom (atom_number_t, atom_number_t = INVALID_ATOM_NUMBER) const;
    int multiple_bond_to_heteroatom (atom_number_t, const int *) const;
    int doubly_bonded_oxygen_count  (atom_number_t) const;

    int identify_side_of_bond  (int * either_side, atom_number_t astart, int flag, atom_number_t avoid) const;

//  All these functions work the same. If the flag (last argument) is set, then there
//  will be no requirement for the atoms to be bonded.

    distance_t    bond_length    (atom_number_t, atom_number_t, int = 0) const;
    angle_t       bond_angle     (atom_number_t, atom_number_t, atom_number_t, int = 0) const;
    angle_t       dihedral_angle (atom_number_t, atom_number_t, atom_number_t, atom_number_t, int = 0) const;

//  After changing a dihedral, we may want to do a bump check

    int bump_check   (atom_number_t, atom_number_t, atom_number_t, atom_number_t, distance_t) const;

//  if different parts of the molecule have been shifted, we can do a bump check between all
//  the '0' atoms and all the '1' atoms

    int bump_check (const int *, distance_t) const;

    int    set_dihedral (atom_number_t, atom_number_t, atom_number_t, atom_number_t, angle_t);

//  A1 and A2 stay fixed. A3 is moved

    int     set_bond_angle (atom_number_t a1, atom_number_t a2, atom_number_t a3, angle_t theta);

//  When setting a bond length the option parameter is the identity of the atom to move

    int     set_bond_length (atom_number_t, atom_number_t, distance_t, atom_number_t = INVALID_ATOM_NUMBER);

//  These functions add one molecule to another, perhaps excluding certain atoms

    int add_molecule         (const Molecule *);        
    int add_molecule_without (const Molecule *, atom_number_t);
    int add_molecule_without (const Molecule *, atom_number_t, atom_number_t);

//  Partial charges and atom types behave in similar ways
                  
    int     has_charges () const;
    void    allocate_charges ();
    void    invalidate_charges ();

//  Set_of_Charges * partial_charges () const;

    const IWString & partial_charge_type () const;

    charge_t charge_on_atom (atom_number_t) const;
    int    copy_charges     (const Molecule &);      // copy charges from another MOL
    int    set_charge       (atom_number_t, charge_t);
    void   set_charges      (const charge_t [], const const_IWSubstring &);

    int    has_atom_types () const;
    void   allocate_atom_types ();
    void   invalidate_atom_types ();

    Atom_Types & atom_types ();

    atom_type_t atom_type  (atom_number_t) const;
    int    copy_atom_types (const Molecule &);      // copy atom types from another molecule
    void   set_atom_type   (atom_number_t, atom_type_t);

    int    compute_Abraham_partial_charges ();
    int    compute_Gasteiger_partial_charges ();
    int    compute_Huckel_partial_charges ();
    int    compute_Gasteiger_Huckel_partial_charges ();
    int    compute_Del_Re_partial_charges ();
    int    compute_Pullman_partial_charges ();

    formal_charge_t formal_charge (atom_number_t) const;
//  formal_charge_t formal_charge_on_atom (atom_number_t) const;
    void   set_formal_charge (atom_number_t, formal_charge_t);

//  We don't want to bother first asking if the charge is already the value we want

    int    set_formal_charge_if_different (atom_number_t, formal_charge_t);

    int    formal_charge () const;    // whole molecule charge

    int transform_atoms (const Element *, const Element *);    // all H's become C

    int    remove_atom (atom_number_t);
    int    remove_atoms (Set_of_Atoms &);
    int    remove_atoms (const int *);
    int    remove_many_atoms (const int *);      // more efficient version - with molecules having lots of atoms
    int    remove_all_atoms_with_isotope (int);
    int    delete_fragment (int);
    int    delete_fragments (const resizable_array<int> &);
    int    delete_fragments (const int *);
    int    remove_fragment_containing_atom (atom_number_t);
    int    delete_all_fragments_except (int);
    int    remove_all  (atomic_number_t);
    int    remove_all  (const Element *);
    int    remove_all_non_natural_elements ();
    int    remove_explicit_hydrogens ();    // need to be treated specially because of H property of adjacent atoms

    int    chop (int);

    int    remove_bonds_to_atom (atom_number_t, int = 0);
    int    remove_bond (int);
    int    remove_bond_between_atoms (atom_number_t, atom_number_t);
    int    remove_all_bonds();
    int    remove_bonds_involving_these_atoms (const int *, int check_chirality = 1);

    molecular_weight_t molecular_weight () const;
    molecular_weight_t molecular_weight_count_isotopes () const;
    molecular_weight_t molecular_weight_ignore_isotopes () const;

    int molecular_weight(const Molecular_Weight_Control &, Molecular_Weight_Calculation_Result &) const;

    Coordinates get_coords (atom_number_t) const;
    int get_coords (atom_number_t, Coordinates &) const;
    int get_coords (Coordinates *) const;
    int vector_between_atoms (atom_number_t, atom_number_t, Coordinates &) const;

    void setx (atom_number_t, coord_t);
    void sety (atom_number_t, coord_t);
    void setz (atom_number_t, coord_t);
    void setxyz (atom_number_t, coord_t, coord_t, coord_t);
    void setxyz (atom_number_t, const coord_t *);
    void setxyz (const Coordinates *);

    void spatial_extremeties (coord_t & xmin, coord_t & xmax, coord_t & ymin, coord_t & ymax) const;
    void spatial_extremeties (coord_t & xmin, coord_t & xmax, coord_t & ymin, coord_t & ymax, coord_t & zmin, coord_t & zmax) const;
    void spatial_extremeties_x (coord_t & xmin, coord_t & xmax) const;
    void spatial_extremeties_x (atom_number_t & left, atom_number_t & right) const;

//  Sometimes we want to know if the molecule has 3D, 2D or "1D" coordinates
//  Function will return 3, 2 or 1 by examining the coordinates. Note that it will
//  be fooled by a molecule with all Z coordinates set to 1.0 - that is 3D

    int  highest_coordinate_dimensionality () const;

    atomic_mass_t atomic_mass (atom_number_t) const;

    exact_mass_t exact_mass () const;
    int exact_mass (exact_mass_t &) const;
    int exact_mass (const int *, exact_mass_t &) const;
    int exact_mass (const int *, int, exact_mass_t &) const;

    void translate_atoms (coord_t, coord_t, coord_t);
    void translate_atoms (const Coordinates &);
    void translate_atoms (coord_t, coord_t, coord_t, const Set_of_Atoms &);
    void translate_atoms (const Coordinates &, const Set_of_Atoms &);
    void translate_atoms (const Coordinates & whereto, const int * to_move, int flag);

    int rotate_atoms (const Coordinates &, angle_t);

    void rotate_to_longest_distance_along_x (atom_number_t & left, atom_number_t & right);

//  Found that I lost accurace with the coord_t version, so make available 
//  versions of rotate_atoms with either Space_Vector<coord_t> or Space_Vector<double>

    template <typename T>
    int rotate_atoms (const Space_Vector<T> &, T, const Set_of_Atoms &);

    int number_fragments ();
    int fragment_membership (atom_number_t);
    int fragment_membership (int *);
    int atoms_in_fragment (int);
    int atoms_in_largest_fragment ();
    int largest_fragment();     // fragment number

//  Fragment data for a subset

    int compute_fragment_information (Fragment_Information &, const int *) const;

    int identify_spinach (int *);      // molecules outside and not between rings
    int identify_spinach_preset (int * spinach) const;      // same, but anything already set in SPINACH will also be included

    int atoms_in_fragment (Set_of_Atoms &, int);               // Set_of_Atoms must start empty
    int add_atoms_in_fragment (Set_of_Atoms &, int);           // Set_of_Atoms may initially contain atoms
    atom_number_t first_atom_in_fragment (int);
    int rings_in_fragment (int);

    template <typename T> int create_components (resizable_array_p<T> &);
    int create_components (int bond_number, Molecule & m1, Molecule & m2);
    int create_components (const int * frag, resizable_array_p<Molecule> &) const;
    int create_components_across_bonds (const int * bonds_to_remove, resizable_array_p<Molecule> &);

//  Retain the fragment containing ZATOM and put all others into FRAGS

    int split_off_fragments (const atom_number_t zatom, Molecule & frags);

//  kind of the opposite

    int excise_fragment (const atom_number_t zatom, Molecule & frag);

//  Different calls to create subset depending on whether the caller
//  provides an array of length natoms () or not

    int create_subset (Molecule &, const int *, int = 1) const;
    int create_subset (Molecule &, const int *, int, int *) const;
    int create_subset_by_bond (Molecule & subset,
                                 const int * these_bonds_only,
                                 int flag) const;

    int reduce_to_largest_fragment ();
    int reduce_to_largest_organic_fragment ();
    int reduce_to_largest_fragment_carefully ();

    int identify_largest_organic_fragment (Set_of_Atoms & atoms_to_be_removed,
                                           int & fragments_same_size_as_largest_organic);

    int organic_only () const;

    int contains_non_periodic_table_elements() const;

    int swap_atoms (atom_number_t, atom_number_t, int call_set_modified = 1);
    int move_atom_to_end_of_atom_list (atom_number_t a);

//  Things using the distance matrix

    int longest_path ();
    int bonds_between (atom_number_t, atom_number_t);
    int atoms_between (atom_number_t, atom_number_t, Set_of_Atoms &);

//  Asking for the implicit hydrogens attached to an atom may force
//  computation of all implicit hydrogens, so no const here.

    int implicit_hydrogens (atom_number_t);

//  This variant also passes back the value of the implicit hcount sticky bit

    int implicit_hydrogens (atom_number_t, int &);

    int explicit_hydrogens (atom_number_t) const;
    int set_implicit_hydrogens (atom_number_t, int, int = 0);
    int set_implicit_hydrogens_known (atom_number_t, int);
    int implicit_hydrogens_known(const atom_number_t) const;
    int unset_all_implicit_hydrogen_information (atom_number_t);
    int hcount (atom_number_t a) { return implicit_hydrogens(a) + explicit_hydrogens(a);}
    int compute_implicit_hydrogens ();
    int implicit_hydrogens ();
    int recompute_implicit_hydrogens (atom_number_t);
    int recompute_implicit_hydrogens ();
    int make_implicit_hydrogens_explicit ();
    int make_implicit_hydrogens_explicit (atom_number_t);
    int make_implicit_hydrogens_explicit (Make_Implicit_Hydrogens_Explicit &);

//  Works with any atomic number

    int move_hydrogens_to_end_of_connection_table (atomic_number_t z = 1);

//  It is often handy to be able to place atoms around another atom. Both A1 and A2 must
//  be singly connected atoms that are bonded to a common anchor. This is used by
//  the make_implicit_hydrogens_explicit () functions

    int set_coordinates_of_singly_connected_atom (atom_number_t zatom,
                                                  coord_t default_bond_length);
    int set_coordinates_of_singly_connected_atoms (atom_number_t a1,
                                                   atom_number_t a2,
                                                   coord_t default_h_bond_length);

//  We may need to know whether or not atomic valences are valid.

    int valence_ok ();
    int valence_ok (atom_number_t);

//  Things for aromaticity

    int compute_aromaticity ();
    int compute_aromaticity_if_needed ();
    int aromaticity (atom_number_t, aromaticity_type_t &);
    int is_aromatic (atom_number_t);
    int is_permanent_aromatic (const atom_number_t &) const;
    int is_aromatic_no_computation (atom_number_t) const;
    int aromaticity (int *);
    int pi_electrons (atom_number_t, int &);
    int lone_pair_count (atom_number_t, int &);
    int contains_aromatic_atoms ();
    int aromatic_atom_count ();
    int aromatic_ring_count ();
    int aromaticity_computed () const;
    int all_rings_containing_atom_are_kekule(const atom_number_t zatom);
    void setDoNotComputeAromaticity(bool valueToSet)
    {
        doNotComputeAromaticity = valueToSet;
    }

//  Optional parameter means don't invalidate the smiles

    int set_aromaticity (atom_number_t, aromaticity_type_t, int = 1);

//  Set the aromaticity of two atoms and the bond between them.

    int set_aromaticity_two_atoms (atom_number_t, atom_number_t, const aromaticity_type_t, int = 1);

//  Somewhat dangerous function. set_modified is not called. The optional
//  argument is passed to the corresponding Atom method

    int set_permanent_aromatic (atom_number_t, int = 1);

    int unset_all_permanent_aromatic_bonds ();
    int unset_all_permanent_aromatic_atoms ();

//  optional arg is for calling set modified. By default it is called

    int change_double_bonds_between_permanent_aromatic_to_single (int call_set_modified = 1);
    int change_ring_bonds_between_permanent_aromatic_to_aromatic (int call_set_modified = 1);

//  Various things for distances

    distance_t distance_between_atoms (atom_number_t, atom_number_t) const;
    distance_t longest_intra_molecular_distance () const;

    int recompute_distance_matrix ();

    const int * distance_matrix_warning_may_change()
    { if (NULL == _distance_matrix)
        recompute_distance_matrix();
      return _distance_matrix;
    }

//  By default, bond types are not updated to aromatic types. If this function
//  is called, that happens.

//  int add_aromaticity_to_bonds ();

    int saturated (atom_number_t) const;

//  Construct a molecule based on a given subset of a molecule

//  int construct_subset (const Path *, Molecule &) const;
//  int construct_subset (const Path *, Molecule &, int *) const;

//  Functions for stereo centres

    int chiral_centres () const { return _chiral_centres.number_elements ();}

    int stereo_centre_hydrogens_become_implicit ();

//  Fetch the Chiral_Centre associated with atom number A

    Chiral_Centre * chiral_centre_at_atom (atom_number_t a) const;
    int valid_chiral_centre (const Chiral_Centre * c) const;

//  Fetch the I'th chiral centre in the molecule

    Chiral_Centre * chiral_centre_in_molecule_not_indexed_by_atom_number (int i) const;

    Chiral_Centre * create_chiral_centre (atom_number_t a, int = 0);
    int print_chiral_centre_details (const Chiral_Centre * c, std::ostream & os) const;

//  The optional parameter governs whether or not we check for an existing chiral centre
//  on that atom. If present, it is replaced with the new chiral centre

    int  add_chiral_centre (Chiral_Centre *, int = 0);
    int  remove_chiral_centre_at_atom (atom_number_t);

    Chiral_Centre * remove_no_delete_chiral_centre_at_atom (atom_number_t a);    // does not delete the Chiral_Centre object - caller must free it

    int  remove_all_chiral_centres ();

//  Set all atoms involved in a chiral centre

    int at_centre_of_chiral_centre (int *, int) const;
    int involved_in_chiral_centre (int *, int) const;

    int revert_all_directional_bonds_to_non_directional ();

    int remove_invalid_directional_bonds ();

    int  invert_chirality_on_atom (atom_number_t);

    int  discern_cis_trans_bonds_from_depiction ();

    int discern_chirality_from_3d_structure();

    int cis_trans_bonds_present() const;

//  during unique smiles determination, we may need to flip the directionality
//  on a group of related cis-trans bonds

    int flip_cis_trans_bond_and_all_related_directional_bonds(const Bond * b);

    int discern_chirality_from_wedge_bonds ();

#ifdef COMPILING_CTB
#include "molecule_ctb.h"
#endif

    int number_records_text_info () const { return _text_info.number_elements ();}
    const IWString & text_info(int i) const { return *(_text_info.item(i));}
    IWString & text_info(int i) { return *(_text_info.item(i));}
    template <typename T> int write_extra_text_info (T &) const;
    int write_extra_text_info (IWString &) const;

    int add_extra_text_info (IWString * extra);
    int copy_extra_text_info_to (Molecule &) const;
    void discard_extra_text_info ();
    int add_extra_text_info (const IWString &);
    int add_extra_text_info (const char *);

    void compute_centre (const Set_of_Atoms *, Coordinates &) const;

    int  centroid  (Coordinates &, int = -1);
    int  centroids (resizable_array_p<Coordinates> &);

    int identify_ez_atoms (atom_number_t a3, atom_number_t a4, atom_number_t & alhs, atom_number_t & arhs) const;
    int ez_by_geometry (atom_number_t a1,
                        atom_number_t a2,
                        atom_number_t a3,
                        atom_number_t a4,
                        angle_t theta = -1.0) const;

#ifdef COMPILING_MOLECULEH_H
#include "moleculeh.h"
#endif

//  This should be a template function, but the implementation goes through
//  qsort, and that seems to make it impossible to use templates. Therefore
//  we have multiple signatures

    int sort (const int *, int = 1);
//  int sort (const float *, int = 1);
//  int sort (const double *, int = 1);

    int renumber_atoms (const int *);

//  const and non-const versions

    IWString smarts_equivalent_for_atom (atom_number_t zatom);
    IWString smarts_equivalent_for_atom (atom_number_t zatom) const;

    const IWString & smarts ();
    const IWString & smarts (Smiles_Information &, const int *);

    const IWString & smarts_starting_with_atom (atom_number_t, Smiles_Information & smi_info, const int *);

//  When we want to make sure that we don't force anything unwanted

    IWString const_smarts_equivalent_for_atom (atom_number_t zatom) const;

    uint64_t quick_atom_hash() const;
    uint64_t quick_bond_hash();    // not const because aromaticity is perceived

//  We are making a smarts of a subset and we want all the D and v directives to specify just
//  the atoms in the subset or a minimum requirement

    int append_smarts_equivalent_for_atom (atom_number_t zatom, IWString & s, const int * include_atom);
    int append_smarts_equivalent_for_atom (atom_number_t zatom, IWString & s);
    int append_smarts_equivalent_for_atom (atom_number_t zatom, IWString & s) const;
    int append_smarts_equivalent_for_atom (atom_number_t zatom, IWString &, const Smiles_Formation_Info &);

    int get_bond_types (bond_type_t *) const;
    int set_bond_types_no_set_modified (const bond_type_t *);
    int set_bond_types_for_isis_aromaticity_matching ();

    int remove_hydrogens_known_flag_to_fix_valence_errors();

    void set_user_specified_void_ptr (void * p) { _user_specified_void_ptr = p;}
    void * user_specified_void_ptr () const { return _user_specified_void_ptr;}

    template <typename F> int while_true_each_user_specified_atom_void_ptr (F &);

    void set_user_specified_atom_void_ptr (atom_number_t, void *);
    void * user_specified_atom_void_ptr (atom_number_t) const;
    void clear_all_user_specified_atom_pointers();
    const Atom * atom_with_user_specified_void_ptr (const void *) const;

    void reset_all_atom_map_numbers();
    void set_atom_map_number(const atom_number_t zatom, const int s);
    int  atom_map_number(const atom_number_t zatom) const { return _things[zatom]->atom_map();}
    atom_number_t atom_with_atom_map_number(const int n) const;

    int unset_unnecessary_implicit_hydrogens_known_values();

    const Atom * const * cbegin() const { return resizable_array_p<Atom>::cbegin();}
    const Atom * const * cend()   const { return resizable_array_p<Atom>::cend();}
    const Atom * const * begin() const { return resizable_array_p<Atom>::cbegin();}
    const Atom * const * end()   const { return resizable_array_p<Atom>::cend();}
    const Bond * const * cbeginBond()  const { return _bond_list.cbegin();}
    const Bond * const * cendBond()    const { return _bond_list.cend();}
    const Ring * const * cbeginRing();
    const Ring * const * cendRing();
    const Chiral_Centre * const * cbeginChiral() const { return _chiral_centres.cbegin();}
    const Chiral_Centre * const * cendChiral() const { return _chiral_centres.cend();}
};

/*
  When writing a Molecule to an std::ostream, we can have it write different types
*/

extern void set_type_to_write_with_operator(int s);

extern std::ostream & operator << (std::ostream & os, Molecule & m);

extern int string_to_file_type (const const_IWSubstring &);
extern int discern_file_type_from_name (const IWString &);
extern const char * suffix_for_file_type (int);
extern int valid_file_type (int);
extern int create_file_with_appropriate_name (const const_IWSubstring &, IWString &, int, int = 0);
extern int append_appropriate_suffix (IWString &, int);
extern int all_files_recognised_by_suffix (const Command_Line &);

/*
  In those cases when two atoms are removed from a molecule, we need a means
  of determining the new atom numbers for those which remain.
*/

#define ADJUSTED_ATOM_NUMBER_1(i, except) \
         ((i) < (except) ? i : i - 1)

#define ADJUSTED_ATOM_NUMBER_2(i, except1, except2) \
         ((i) < (except1) ? ADJUSTED_ATOM_NUMBER_1((i), (except2)) : \
                            ADJUSTED_ATOM_NUMBER_1((i), (except1)) - 1)

/*
  Molecule objects should be non-NULL and also pass their OK function
*/

#define OK_MOLECULE(m) ( (NULL != (m)) && ((Molecule *) (m))->ok ())

/*
  Frequently we need to check whether or not a molecule and a given
  atom within it are OK.
*/

#define OK_ATOM_NUMBER(m, a) ( (m)->ok_atom_number ((a)))

/*
  Some macro's for checking more than one atom within the same
  molecule. Note that we impose the restriction of the atoms being
  distinct.
*/

#define OK_2_ATOMS(m, a1, a2)         (NULL != (m) && (m)->ok_2_atoms ((a1), (a2)))
#define OK_3_ATOMS(m, a1, a2, a3)     (NULL != (m) && (m)->ok_3_atoms ((a1), (a2), (a3)))
#define OK_4_ATOMS(m, a1, a2, a3, a4) (NULL != (m) && (m)->ok_4_atoms ((a1), (a2), (a3), (a4)))

extern void set_mdl_ignore_self_bonds       (int);
extern void set_write_mdl_charges_as_m_chg  (int);
extern void set_mdl_write_aromatic_bonds    (int);
extern void set_mdl_write_aromatic_atoms    (int);
extern void set_mdl_allow_deuterium         (int);
extern void set_mdl_allow_tritium           (int);
extern void set_write_mdl_chiral_flags      (int);
extern void set_write_v30_mdl_files         (int);
extern void set_write_mdl_dollars           (int);
extern void set_store_pdb_atom_information  (int);
extern void set_write_mdl_m_end_record      (int);

extern const resizable_array_p<PDB_Stored_Atom_Information> & stored_pdb_atom_information_last_molecule_read();
// pass any int to distinguish from the const version
extern       resizable_array_p<PDB_Stored_Atom_Information> & stored_pdb_atom_information_last_molecule_read(int);

//extern void set_mdl_write_m_end             (int);
extern void set_include_chiral_info_in_mdl_outputs (int);
//extern void set_include_chirality_in_mdl_outputs   (int);

extern void set_mdl_read_isotopes_as_numbers_rather_than_differences_from_normal (int s);
extern void set_mdl_read_M_isotopes_as_numbers_rather_than_differences_from_normal (int s);
extern void set_mdl_write_isotopes_as_numbers_rather_than_differences_from_normal (int s);
extern void set_mdl_write_M_isotopes_as_numbers_rather_than_differences_from_normal (int s);
extern void set_mdl_insert_between_sdf_name_tokens (const const_IWSubstring &);
extern void set_write_isis_standard    (int);
extern void set_prepend_sdf_identifier (int);
extern int  add_rdfile_identifier (const IWString &);
extern int  set_rdfile_start_of_record (const const_IWSubstring & s);
extern void set_ignore_unrecognised_mdl_m_records (int);
extern void set_die_on_erroneous_m_input (int);
extern void set_mdl_report_unrecognised_records   (int);

extern void set_write_pdb_files_in_fragment_order(int);
extern void set_pdb_number_by_element_count(int);
extern void set_pdb_number_within_sequence(int);
extern void set_use_stored_atom_information_when_writing_pdb_files (int s);
extern void set_replace_first_sdf_tag (const IWString & s);

/*
  When [ND3v4+0] is read, we can optionally automatically convert it to N+
*/

extern void set_put_formal_charges_on_neutral_ND3v4 (int);
extern int  put_formal_charges_on_neutral_ND3v4 ();

extern void set_mdl_accumulate_mdl_chirality_features (int);
extern void set_mdl_truncate_long_elements (int s);
extern void set_mdl_change_long_symbols_to (const const_IWSubstring & s);
extern void set_mdl_discern_chirality_from_wedge_bonds (int);
extern void set_mdl_read_h_correct_chiral_centres (int);
extern void set_mdl_write_h_correct_chiral_centres (int);
extern int  set_mdl_input_bond_type_translation (int zfrom, int zto);
extern int  mdl_discern_chirality_from_wedge_bonds ();

extern int  ignore_self_bonds ();

template <typename T> int read_next_v30_record (T & input, IWString & buffer);

extern off_t seek_to_from_command_line ();
extern void  set_seek_to (off_t);

extern off_t max_offset_from_command_line ();
extern void  set_max_offset_from_command_line (off_t);

extern void set_mol2_assign_default_formal_charges (int);
extern void set_mol2_write_assigned_atom_types (int s);
extern void set_place_mol2_residue_information_in_user_specified_void_ptr (int s);
extern void set_mol2_write_formal_charge_as_partial_charge (int s);
extern void set_mol2_read_charge_column_contains_formal_charges(int s);

extern int  ignore_all_chiral_information_on_input ();
extern void set_ignore_all_chiral_information_on_input (int);

extern void set_ignore_tdts_with_no_smiles (int);
extern void set_smiles_tag (const const_IWSubstring &);

extern int  flush_files_after_writing_each_molecule ();
extern void set_flush_files_after_writing_each_molecule (int);

extern int  ignore_incorrect_chiral_input ();
extern void set_ignore_incorrect_chiral_input (int);
extern void set_automatically_add_implicit_hydrogen_to_incomplete_chiral_centre (int s);

extern int  set_sdf_identifier (const const_IWSubstring &);
extern void set_extract_isis_extregno (int e);
extern void set_mdl_name_in_m_tag (const const_IWSubstring &);
extern void set_fetch_all_sdf_identifiers (int);
extern void set_take_first_tag_as_name (int s);
//extern void set_isis_output_name_tag (const const_IWSubstring &);
extern void set_discard_sdf_molecule_name (int);
extern void set_display_non_organic_chirality_messages (int);
extern void set_mdl_display_invalid_chiral_connectivity (int);

extern void set_tdt_identifier_dataitem (const const_IWSubstring &);
extern void set_tdt_append_dataitem (const const_IWSubstring &);

extern void set_discern_cis_trans_bonds (int);
extern int  discern_cis_trans_bonds ();

extern void set_discern_chirality_from_3d_coordinates(int);
extern int  discern_chirality_from_3d_coordinates();

extern void set_ignore_bad_cis_trans_input (int s);
extern int ignore_bad_cis_trans_input ();
extern void set_discard_directional_bonds_on_input (int);

#ifdef BONDS_KNOW_RING_MEMBERSHIP

// It is optional whether or not bonds know about their ring membership.

//extern int set_bonds_get_ring_info (int);

#endif

extern void set_read_extra_text_info (int);
extern int  read_extra_text_info ();
extern void set_write_extra_text_info (int);
extern int  write_extra_text_info ();
extern char input_file_delimiter ();
extern int  input_is_dos_mode ();

extern void set_write_DOS_records (int s);
extern int  write_DOS_records ();
extern const IWString & newline_string();

extern void set_skip_first_molecules (int);
extern int  skip_first_molecules     ();
extern void set_do_only_n_molecules  (int);
extern int  do_only_n_molecules      ();

//extern Molecule * next_molecule (iwstring_data_source &, int);

//extern int next_molecule (Molecule &, iwstring_data_source &, int);

//extern int construct_ring (resizable_array<const Bond *> & bonds, Ring * r);

class Command_Line;

extern int process_file_types  (const Command_Line &, int &, int &);
extern int process_input_type  (const Command_Line &, int &);
extern int process_output_type (const Command_Line &, int &);
//extern int determine_output_types (const Command_Line &, resizable_array<int> &);

//extern int process_file_types  (const char *, int &, const char *, int &);

//extern int identify_fused_rings (Molecule * m, resizable_array<int> & rings, int count);

extern int construct_all_paths_between (Molecule & m, atom_number_t a1, atom_number_t a2,
                             resizable_array_p<Path> & paths);
extern int construct_all_paths_between (Molecule & m, atom_number_t a1, atom_number_t a2,
                             resizable_array_p<Set_of_Atoms> & paths);

// When adding Hydrogens, we can set the bond length used

extern void set_default_h_bond_length (coord_t d);

// For Sphinx, we make it optional whether or not adding the same
// bond twice to a molecule is fatal

extern void set_add_same_bond_twice_fatal (int);

/*
  Do we issue a printed warning when we are asked for the molecular
  weight of a molecule with non-periodic table elements in it
*/

extern void set_issue_non_periodic_table_molecular_weight_warning (int s);

/*
  Do we issue warnings when we encounter a likely invalid valence? This
  is particularly important with trxn, where a great many invalid
  intermediates get formed.
*/

extern void set_display_messages_about_unable_to_compute_implicit_hydgogens (int s);

extern void set_display_already_bonded_error_message(int s);

/*
  We can optinally compute only the upper half of the distance matrix - by default full
*/

extern void set_full_distance_matrix(const int s);


/* 
Setting the char to be used for the "dot" in 	SMILES - this is used by LIlly reaction stuff where + is a separator
*/

extern void setUseThisCharAsDotInSmiles(char thisChar);
extern char getUseThisCharAsDotInSmiles(void);

/*
  We can temporarily disable such messages with an object. Messages get
  turned back on when this goes out of scope
*/

class Temporarily_Disable_Messages_About_Unable_to_Compute_Implicit_Hydrogens
{
  private:
  public:
    Temporarily_Disable_Messages_About_Unable_to_Compute_Implicit_Hydrogens ()
    {
      set_display_abnormal_valence_messages(0);
      set_display_messages_about_unable_to_compute_implicit_hydgogens(0);
    }
    ~Temporarily_Disable_Messages_About_Unable_to_Compute_Implicit_Hydrogens () 
    {
      set_display_abnormal_valence_messages(1);
      set_display_messages_about_unable_to_compute_implicit_hydgogens(1);
    }
};

int is_actually_chiral (Molecule & m, atom_number_t zatom);

extern int unconnect_covalently_bonded_non_organics_on_read ();
extern void set_unconnect_covalently_bonded_non_organics_on_read (int);

/*
  When asking for smarts from a molecule, we need to decide what kind of
  smarts do we want. Do we want a smarts that as completely describies the
  atom(s) in the current molecule, or are we looking for something that
  can be used for substructure searches?
*/

extern void set_make_smarts_embedding (int);
extern int  make_smarts_embedding ();

/*
  When making a smarts embedding, do we include the unsaturation of the atom.
  This comes up when we, for example, pass through a carbonyl group, without
  including the Oxygen atom
*/

extern void set_include_unsaturation_in_smarts_embeddings (int);
extern int  include_unsaturation_in_smarts_embeddings ();

/*
  When making a smarts out of a subset, do we make the ring membership
  that of the whole parent molecule, or just the subset that is
  part of the subset
*/

extern void set_include_whole_molecule_ring_membership_in_subset_smarts (int);
extern int  include_whole_molecule_ring_membership_in_subset_smarts ();

extern void set_accumulate_non_sssr_rings (int);
extern int  accumulate_non_sssr_rings ();

extern void set_include_D_in_smarts(int);
extern int  include_D_in_smarts();

/*
  When doing substructure searches with queries that depend on having implicit rather than
  explicit hydrogens, we can temporarily break the bonds between Hydrogens and heavy atoms.

  First call unattach_hydrogens (), then call reattach_atoms ()
*/

/*extern int unattach_atoms (Molecule & m,
                    int * connection,
                    const atomic_number_t z = 1);*/

/*extern int reattach_atoms (Molecule & m,
                const int * connection,
                const bond_type_t bt = SINGLE_BOND);*/

/*
  Mar 2004. We are getting files with arbitrary partial charge values.
  Convert the REASONBLE_CHARGE_VALUE macro to a function so we can
  alter the threshold
*/

extern int reasonable_atomic_partial_charge_value (charge_t);

extern int set_max_reasonble_atomic_partial_charge_value (charge_t);
extern int set_min_reasonble_atomic_partial_charge_value (charge_t);
extern int set_reasonable_atomic_partial_charge_range (charge_t, charge_t);

extern int number_connection_table_errors_to_skip ();
extern void set_number_connection_table_errors_to_skip (int);

/*
  Determine the number of atoms in a smiles just by examining text
*/

extern int count_atoms_in_smiles(const const_IWSubstring & smiles);

extern void set_copy_name_in_molecule_copy_constructor(int);

extern int is_fixed_kekule_form (const Molecule & m, const Ring & r);

extern int inchi_to_inchi_key(const char * inchi,IWString & key);

template <typename F>
void
Molecule::each_atom (F & f) const
{
  assert (ok());

  resizable_array_p<Atom>::each(f);
}

template <typename F>
void
Molecule::each_atom_lambda (F f) const
{
  assert (ok());

  resizable_array_p<Atom>::each_lambda(f);
}

template <typename F>
void
Molecule::each_bond (F & f) const
{
  assert (ok());

  _bond_list.each(f);
}

template <typename M>
int
Molecule::while_true_each_user_specified_atom_void_ptr (M & f)
{
  Atom * const * a = _things;

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == f(a[i]->user_specified_void_ptr()))
      return i;
  }

  return _number_elements;
}

extern void reset_all_file_scope_variables();
extern void reset_mdl_file_scope_variables();
extern void reset_unique_file_scope_variables();
extern void reset_chiral_centre_file_scope_variables();

extern void set_copy_user_specified_atom_void_ptrs_during_create_subset (int);
extern void set_copy_atom_based_user_specified_void_pointers_during_add_molecle (int);

extern void set_exclude_triple_bonds_from_graph_reduction (int);

/*
  Forces a double bond to be placed at all ring fusion points - helps when pulling
  rings apart
*/

extern int arrange_kekule_forms_in_fused_rings (Molecule & m);

extern void set_invalidate_bond_list_ring_info_during_invalidate_ring_info (int s);

template <typename T>
int write_isotopically_labelled_smiles (Molecule &, const bool uniq, T &);


#endif

/* arch-tag: 2137d4c2-7ac1-4028-8347-e12e3515b13f

*/
