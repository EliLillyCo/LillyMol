#ifndef MDL_FILE_DATA_H
#define MDL_FILE_DATA_H

#include "substructure.h"
#include "mdl_atom_record.h"

class ISIS_Atom_List
{
  private:
    int _normal_list;     // the T or F value

    resizable_array<const Element *> _element;

  public:
    ISIS_Atom_List ();

    int debug_print (std::ostream &) const;

    int active() const { return _element.number_elements();}
    int number_elements() const { return _element.number_elements();}

    int write_M_ALS (atom_number_t, std::ostream &) const;

    int normal_list () const { return _normal_list;}
    void set_normal_list (int s) { _normal_list = s;}

    int create_from_ALS_record (const IWString &);
    int initialise_from_mdl_A_symbol ();
    int initialise_from_mdl_AH_symbol ();
    int initialise_from_mdl_Q_symbol ();
    int initialise_from_mdl_QH_symbol ();
    int initialise_atom_list_from_symbol (const const_IWSubstring &);

    int invert_to_normal_list_based_on_organics ();

    bool operator == (const ISIS_Atom_List &) const;
    bool operator != (const ISIS_Atom_List & rhs) const { return ! operator== (rhs);}

    const Element * elementi(int i) const { return _element[i];}

    void swap_atoms (atom_number_t, atom_number_t) { return;} // nothing to do

    int convert_not_atom_lists_to_organic_lists ();
};

/*
  There is a bunch of stuff that can be read from MDL files.
  This later gets converted to query and reaction objects.

  These objects perform two functions:
    First, read and store the data from an MDL file, that's easy.

  Then as the associated molecule is examined, various things
  can be discerned about how the atoms need to be translated
  to substructure query specifications, so there are also 
  variables that describe this. 
*/

class MDL_Atom_Data
{
  private:

    atom_number_t _atom_number;

//  These pieces of information come from the atom records

    int _hcount;
    int _h0designator;
    int _valence;

//  These come from M records

    int _unsaturated;
    int _substitution;     // note that I extend this with the value -3, which means ignore. I should hvae done that with 0, but too late...
    int _ring_bond;       // a value of -3 means ignore

    IWString _alias;

    ISIS_Atom_List _atom_list;

//  Some things specific to reactions - read from the atom records

    int _atom_map;
    int _inversion;
    int _exact_change;

//  Now derived information that can be imposed by any programme on
//  the way to forming a substructure query

//  Sometimes we have just a minimum number of substitutions specified

    int _min_ncon;
    int _max_ncon;

// When parsing an ISIS reaction file we may be able to discern that
// various atoms need a minimum number of Hydrogen atoms

    int _min_hcount;

//  Usually when a molecule is read, we will remove any explicit hydrogen
//  atoms. We need to keep track of that

    int _explicit_hydrogen_atoms_removed;

//  When we remove link atoms, the atoms to which they were connected
//  need to know that they lost a connection.

    int _connections_lost;

// An ISIS query may have OR aromatic bonds.  The molecule will have
// single bonds, but aromaticity may in fact be known

    int _aromatic;

//  Atoms get read in as things like A, Q and L. We need to keep track of
//  that information

    IWString _initial_atomic_symbol;

//  One of the hardest things is to keep the association between
//  the ordering of the atoms in a Molecule and these items. So,
//  keep a pointer

    const Atom * _atom;

//  Jul 2007. Not read from an MDL file, but sometimes you know something
//  about the number of bonds to an atom. Initial use in remove_and_label.cc

    int _nbonds;

//  don't forget the copy constructor if you add things here

//  private functions

    void _default_values();
    void _do_copy (const MDL_Atom_Data & rhs);

  public:
    MDL_Atom_Data ();
    MDL_Atom_Data (const MDL_Atom_Data &);

    MDL_Atom_Data & operator= (const MDL_Atom_Data &);

    int write_M_ALS (atom_number_t, std::ostream &) const;

    void set_atom_number (atom_number_t a) { _atom_number = a;}
    atom_number_t atom_number () const { return _atom_number;}

    int extract_info_from_mdl_file_record (const MDL_Atom_Record & mar);

    const ISIS_Atom_List & atom_list () const { return _atom_list;}

    int initialise_atom_list_from_symbol (const const_IWSubstring &);

    void set_alias (const const_IWSubstring & s) { _alias = s;}
    const IWString & alias () const { return _alias;}

    void set_atom_map (const int s) { _atom_map = s;}    // dangerous

    const IWString & initial_atomic_symbol () const { return _initial_atomic_symbol;}

    int build_atom_list (const const_IWSubstring &);
    int convert_a_or_q_atoms_to_atom_list (const IWString &);

    void set_unsaturation (int s) { _unsaturated = s;}
    void set_substitution (int s) { _substitution = s;}
    void set_ring_bond (int s) { _ring_bond = s;}
    void set_exact_change (int s) { _exact_change = s;}

    void increment_connections_lost () { _connections_lost++;}
    void increment_min_hcount () { _min_hcount++;}
    void set_min_hcount (int s) { _min_hcount = s;}
    void increment_explicit_hydrogen_atoms_removed () { _explicit_hydrogen_atoms_removed++;}
    void decrement_substitution () { _substitution--;}

    void set_hcount (int s) { _hcount = s;}

    int min_hcount () const { return _min_hcount;}

    int explicit_hydrogen_atoms_removed () const { return _explicit_hydrogen_atoms_removed;}

    int connected_atom_is_being_removed(atomic_number_t);

    int aromatic () const {return _aromatic;}
    void set_aromatic (int s) { _aromatic = s;}

    void set_min_ncon (int s) { _min_ncon = s;}
    void set_max_ncon (int s) { _max_ncon = s;}   // no check for consistency
    int  min_ncon () const { return _min_ncon;}
    int  max_ncon () const { return _max_ncon;}

    void set_nbonds (int s) { _nbonds = s;}
    int  nbonds() const { return _nbonds;}

    int atom_map () const { return _atom_map;}

    int hcount () const { return _hcount;}
    int inversion () const { return _inversion;}
    int substitution () const { return _substitution;}
    int unsaturated () const { return _unsaturated;}
    int ring_bond () const { return _ring_bond;}
    int h0designator () const { return _h0designator;}
    int valence () const { return _valence;}

    int exact_change () const { return _exact_change;}

    int swap_atoms (atom_number_t, atom_number_t);

    int convert_not_atom_lists_to_organic_lists ();
};

class MDL_Bond_Data
{
  private:
    int _bond_type_read_in;

//  the _bond_type_read_in may contain OR'd combinations of regular bond types

    bond_type_t _btype;

    int _bond_topology;

    int _reacting_centre_status;   // not used anywhere

    int _cfg;

  public:
    MDL_Bond_Data();
    MDL_Bond_Data (const MDL_Bond_Data &);

    int extract_info_from_mdl_file_record (const MDL_Bond_Record & mbr);

    bond_type_t btype () const {return _btype;}
    void set_btype (bond_type_t s) { _btype = s;}

    bond_type_t bond_type_read_in () const { return _bond_type_read_in;}

    int is_or_aromatic () const;

    int bond_topology () const { return _bond_topology;}
    void set_bond_topology (int s) { _bond_topology = s;}

    int reacting_centre_status () const { return _reacting_centre_status;}

    void swap_atoms (atom_number_t, atom_number_t) { return;}   // don't think there is anything to do
};

class MDL_File_Data
{
  protected:
    IWString _third_line_of_input_sdf_file;

    resizable_array_p<MDL_Atom_Data> _mdl_atom;
    resizable_array_p<MDL_Bond_Data> _mdl_bond;

    resizable_array_p<Link_Atom> _link_atom;

//  private functions

    int _do_copy (const MDL_File_Data &);

  public:
    MDL_File_Data();
    MDL_File_Data(const MDL_File_Data &);
    ~MDL_File_Data();

    MDL_File_Data & operator= (MDL_File_Data && rhs);

//  We make things much easier if we allocate all the arrays the moment na and nb are known

    int allocate_arrays(int na, int nb);

//  Sometimes we really don't have any info, but we need to have an object
//  that matches a given molecule

    int build (Molecule &);

    const IWString & third_line_of_input_sdf_file() const { return _third_line_of_input_sdf_file;}

    int arrays_allocated() const { return _mdl_atom.number_elements();}
    int active() const { return _mdl_atom.number_elements();}

    void discard_atom_map ();

    int number_atoms() const { return _mdl_atom.number_elements();}
    int number_bonds() const { return _mdl_bond.number_elements();}

    const MDL_Atom_Data * mdl_atom_data(int i) const { return _mdl_atom[i];}
    MDL_Atom_Data * mdl_atom_data(int i) { return _mdl_atom[i];}
    const MDL_Atom_Data * mdl_atom(int i) const { return _mdl_atom[i];}
    MDL_Atom_Data * mdl_atom(int i) { return _mdl_atom[i];}

    const MDL_Bond_Data * mdl_bond_data(int i) const;
    MDL_Bond_Data * mdl_bond_data(int i);
    const MDL_Bond_Data * mdl_bond(int i) const;
    MDL_Bond_Data * mdl_bond(int i);

//  Do any of the bonds involve an aromatic bond in an OR operation

    int or_aromatic_bonds_present() const;

    int number_link_atoms() const { return _link_atom.number_elements();}

    const Link_Atom * link_atom(int ndx) const { return _link_atom[ndx];}

    void remove_link_atoms() { _link_atom.resize(0);}

    const resizable_array_p<Link_Atom> & link_atoms() const { return _link_atom;}
    int transfer_link_atoms (resizable_array_p<Link_Atom> & s) { return s.transfer_in(_link_atom);}

    int swap_atoms (atom_number_t, atom_number_t);

    int add (const Element * e);
};

/*
  When we read in bond types from query files, they can be things like "single or double".
  From this we need two things:
    a bond that actually encodes all that info
    a bond that can be used for building a molecule - defaults to single
*/

extern int convert_mdl_bond_type_read_in_to_query (int, bond_type_t & for_query, bond_type_t & for_building_a_molecule);

#endif
