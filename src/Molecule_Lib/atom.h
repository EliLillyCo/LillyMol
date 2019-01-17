/*
  Atom class. 
  Note that connection information can be held by atoms.
*/

#ifndef IW_ATOM_H
#define IW_ATOM_H 1

#include <iostream>

#include "iwmtypes.h"
#include "element.h"
#include "coordinates.h"

#define ATOM_PROPERTY_UNKNOWN -476

#include "bond.h"
#include "set_of_atoms.h"

#include "iwaray.h"

class Atom : public resizable_array <Bond *>, public Coordinates
{
  private:
    const Element *_element;

//  Normally the implicit hydrogen count is a computed value. Sometimes
//  however a known hcount value is read from a file.

    short _implicit_hydrogens;
    short _implicit_hydrogens_known;    // known value from file

    short _formal_charge;
    short _radical;
    short _nbonds;

//  Some atoms may be classified as permanently aromatic

    short _permanent_aromatic;

    int _isotope;
    int _userAtomType;

    void * _user_specified_void_ptr;

    int _atom_map;

//  private functions

    void _default_values (const Element *);
    void _constructor_copy_atom_attributes (const Atom & other_atom);
    int _compute_implicit_hydrogens (int & result);

    template <typename T> int _common_saturation (const T & comparitor) const;

  public:
    Atom (const char *);
    Atom (const Element *);
    Atom (atomic_number_t);
    Atom (const Atom *);                // make a copy of an atom - connections are not copied
    Atom (const Atom &);
    ~Atom ();

    int  debug_print (std::ostream &) const;

    int ok () const { return NULL != _element;}

    int audit () const;

    int  isotope () const { return _isotope;}
    void set_isotope (int i) { _isotope = i;}
    int  is_isotope () const;
    
    int  userAtomType () const { return _userAtomType;}
    void set_userAtomType (int i) { _userAtomType = i;}
 
    int add (Bond *);

    int ncon () const { return _number_elements;}    // how many atoms connected

    int ncon (atom_number_t, const int * include_atom) const;   // how many of the atoms set in the array are connected

    int nbonds ();
    int nbonds () const;
    int recompute_nbonds ();

    int permanent_aromatic () const { return _permanent_aromatic;}
    void set_permanent_aromatic (int s) { _permanent_aromatic = s;}

    int molecule_being_resized (int);

    const Bond * bond_to_atom (atom_number_t myAtomId,atom_number_t otherAtomId) const;

//  int is_hydrogen () const;
//  int is_carbon () const;
//  int is_asterisk () const;
//  int is_hydrogen_or_asterisk () const;
    int is_halogen () const;

    const IWString & atomic_symbol () const;
    const Element * element () const { return _element;}
    const Element & elementq () const { return *_element;}
    atomic_number_t atomic_number () const;

    int atomic_symbol_hash_value () const { return _element->atomic_symbol_hash_value ();}

    void set_element (const Element *);    // dangerous to have public

    int implicit_hydrogens_computed () const { return _implicit_hydrogens >= 0;}
    int implicit_hydrogens ();       // not const as it stores the result
    int recompute_implicit_hydrogens (int &);
    int compute_implicit_hydrogens (int &);   // note this does not store the result, it is non const because it may compute nbonds
    int set_implicit_hydrogens (int, int = 0);

    void unset_all_implicit_hydrogen_information();

//  during Kekule determinations I may change the bonds to an atom and need a quick means
//  of having that atom reset anything that may be dependent on its bonding

    void invalidate_computed_values_after_bond_change ();

    int  implicit_hydrogens_known () const { return _implicit_hydrogens_known;}
    void set_implicit_hydrogens_known (int i);

    formal_charge_t formal_charge () const { return formal_charge_t (_formal_charge);}
    int set_formal_charge (formal_charge_t);

//  int nrings_computed () const { return ATOM_PROPERTY_UNKNOWN != _nrings;}
//  int nrings () const;
//  int set_nrings (int r);
//  int set_is_ring_atom ();
//  int set_is_non_ring_atom ();
//  void invalidate_nrings () { _nrings = ATOM_PROPERTY_UNKNOWN;}
//  int in_another_ring ();

    int  is_radical() const { return _radical;}
    void set_radical(const int s ) { _radical = s;}

    void set_modified ();

//  void set_cb (int, atom_number_t, bond_type_t);
//  int  add_con (atom_number_t, bond_type_t);
    int  remove_bonds_to_atom (atom_number_t, int = 0);

    atom_number_t other (atom_number_t, int) const;    // atom number of i'th connection
    bond_type_t   btype_to_connection (int) const;    // bond order  to i'th connection
    bond_type_t   btype_to_atom (atom_number_t myAtomId, atom_number_t otherAtomId) const;    // bond order of bond to otherAtomId

    int other_and_type (atom_number_t, int, atom_number_t &, bond_type_t &) const;

//  void set_bond_type_to_connection (int, bond_type_t);    // set bond order of i'th connection
    int  set_bond_type_to_atom (atom_number_t, bond_type_t);

//  void set_atom_number (int, int);        // set atom number of i'th connection
//  void set_directional_bond (int, int);    // set the I'th bond to be directional

    int is_centre_of_cis_trans_bond () const;

    int number_directional_bonds_attached () const;

    int connections (atom_number_t, atom_number_t *, bond_type_t * = NULL) const;
    int connections (atom_number_t, Set_of_Atoms &) const;

    int connections_and_types (atom_number_t, Set_of_Atoms &,
                               resizable_array<bond_type_t> &) const;
    int connections_and_types (atom_number_t, atom_number_t *, bond_type_t *) const;

    int bond_types (resizable_array<bond_type_t> &) const;
    int bond_types (bond_type_t *) const;

    int is_bonded_to (atom_number_t) const;       // are we bonded to atom number J
    int is_bonded_to (atom_number_t, bond_type_t &) const;

    atomic_mass_t atomic_weight () const;

    int exact_mass (exact_mass_t &) const;

    int lone_pair_count (int &);
    int pi_electrons (int &);

    int valence_ok ();

//  When making a smiles, it is often convenient to have multiple bonds first

    void multiple_bonds_first ();

//  various things for smiles

    int next_atom_for_smiles (atom_number_t my_atom_number,
                            int * already_done,
                            int * canonical_order,
                            atom_number_t & next_atom) const;
    int next_atom_for_unique_smiles (atom_number_t my_atom_number,
                            int * already_done,
                            int * canonical_order,
                            atom_number_t & next_atom) const;
    int next_atom_for_random_smiles (atom_number_t my_atom_number,
                            int * already_done,
                            int * canonical_order,
                            atom_number_t & next_atom) const;

//  Several molecule formats need coordinates written in a common format

    int write_coordinates (std::ostream &, int = 0) const;

    int fully_saturated () const;     // nbonds () == ncon ()
    int unsaturated () const;         // nbonds () < ncon ()

    const void * user_specified_void_ptr () const { return _user_specified_void_ptr;}
    void * user_specified_void_ptr () { return _user_specified_void_ptr;}
    void set_user_specified_void_ptr (void * v) {_user_specified_void_ptr = v;}

    int remove_connections_to_any_of_these_atoms (const int *);

    int atom_map () const { return _atom_map;}
    void set_atom_map (const int s) { _atom_map = s;}
};
extern void reset_atom_file_scope_variables();
extern angle_t angle_between_atoms (const Atom &, const Atom &, const Atom &, const Atom &);

#define OK_ATOM(a) (NULL != (a) && (a)->ok () )

extern int how_many_atoms ();

extern Coordinates & form_unit_vector (const Atom &, const Atom &);

extern void set_display_abnormal_valence_messages (int);
extern int  display_abnormal_valence_messages ();

extern void set_copy_implicit_hydrogen_count_in_atom_copy_constructor (int s);

extern int set_reasonable_formal_charge_range(formal_charge_t, formal_charge_t);
extern int reasonable_formal_charge_value(formal_charge_t c);
extern void set_reset_implicit_hydrogens_known_on_bond_removal(const int s);
extern void set_four_connected_neutral_nitrogen_has_h(const int s);

extern void set_alternate_valences_give_hcount(const int s);
#endif
