/*
  Elements are mostly defined by their atomic numbers.
*/

#ifndef ELEMENT_H
#define ELEMENT_H

// Make sure you adjust the alphabetic_element_symbol_order array in molecule.cc when adding new elements

#define HIGHEST_ATOMIC_NUMBER 118

#define REASONABLE_ATOMIC_NUMBER(z) ((z) >= 0 && (z) <= HIGHEST_ATOMIC_NUMBER)

class IWString;
class const_IWSubstring;

#include <iostream>

#include "iwstring.h"

#include "iwmtypes.h"

#define LEN_AT_SYMBOL 4

#define ELEMENT_MAGIC -71773

#define NOT_AN_ELEMENT -1

#define VALENCE_NOT_DEFINED -5

class Element {
  private:
    atomic_number_t _atomic_number;
    IWString _symbol;
    IWString _aromatic_symbol;        // save doing run-time conversions
    int _normal_isotope;              // count of protons + neutrons for most abundant
    atomic_mass_t _atomic_mass;       // closely related to _normal_isotope
    int _organic;
    int _metal;
    int _normal_valence;
    int _outer_shell_electrons;                // only when in normal valence state

//  Alternate valences should be stored in increasing order. Currently there is no
//  check for this.

    resizable_array<int> _alternate_valence;

//  Rick wanted exact masses for Mass Spec

    exact_mass_t  _exact_mass;

    int _needs_square_brackets;

//  Oct 2002. Store the hash value - based on the symbol

    int _atomic_symbol_hash_value;

//  Jan 2004. Some elements can be designated as always aromatic

    int _permanent_aromatic;

//  June 2016. Need some quick means of uniquely identifying elements.

    int _unique_id;

// October 2016. When we are dealing with atoms as peptides, handy to know
// if a peptide is a natural one or not.

    int _natural_peptide;

//  private functions

    void _non_periodic_table_element_constructor (const char * s, int nchars);

    void _set_symbol (const char *);
    void _set_symbol (const char *, int);
    void _set_symbol (const const_IWSubstring &);
    void _default_values (atomic_number_t);

  public:
    Element (atomic_number_t);
    Element (const char *);
    Element (const char *, int nchars);
    Element (const const_IWSubstring &);
    ~Element ();

    int ok () const;      // audit function
    int debug_print (std::ostream &) const;

    atomic_number_t atomic_number () const;

//  When we create permanent aromatic elements, like 'c' and 'n', they are just non-periodic
//  table entities. We need a means ot letting them know that they are proper elements

    void copy_element_data (const Element *);

    const IWString & symbol () const { return _symbol;}
    const IWString & aromatic_symbol () const { return _aromatic_symbol;}
    void set_aromatic_symbol (const char * s) { _aromatic_symbol = s;}

    int atomic_symbol_hash_value () const { return _atomic_symbol_hash_value;}

    int unique_id() const { return _unique_id;}
    void set_unique_id(int s) { _unique_id = s;}   // should never be set by the user

    int read_ptable_record (const const_IWSubstring &);

    int append_smiles_symbol (IWString & smiles, aromaticity_type_t arom, int isotope) const;

    int normal_isotope () const { return _normal_isotope;}

    int organic () const { return _organic;}
    void set_organic (int o) { _organic = o;}
    int is_in_periodic_table () const;
    int is_metal () const { return _metal;}

    int normal_valence () const { return _normal_valence;}

//  May 2002. Competing requirements for pathd and sample ID. SampleID wants things like Ca
//  to have a valence, pathd doesn't.

    void set_normal_valence (int v) { _normal_valence = v;}

    int number_alternate_valences () const { return _alternate_valence.number_elements ();}
    int alternate_valence (int) const;

//  Is a valence either the normal or one of the alternate ones?

    int is_valid_valence (int) const;

    int add_alternate_valence (int av) { return _alternate_valence.add (av);}  // should check for legal values of av

    int needs_square_brackets () const { return _needs_square_brackets;}
    void set_needs_square_brackets (int n) { _needs_square_brackets = n;}

    atomic_mass_t atomic_mass () const;
    void set_atomic_mass (atomic_mass_t a) { _atomic_mass = a;}

    exact_mass_t exact_mass () const { return _exact_mass;}

    int is_halogen () const;

    int pi_electrons (int, formal_charge_t, int &) const;
    int lone_pairs   (int, formal_charge_t, int &) const;

    int outer_shell_electrons (int &) const;

    int permanent_aromatic () const { return _permanent_aromatic;}
    void set_permanent_aromatic (int s) { _permanent_aromatic = s;}

    int natural_peptide() const { return _natural_peptide;}
    void set_natural_peptide(const int s) { _natural_peptide = s;}
};

extern void  debug_print_all_elements (std::ostream &);

/*
  All the get_element_from_symbol variants convert the first char to uppercase
*/

extern const Element * get_element_from_symbol (const char *, int, int &);
extern const Element * get_element_from_symbol (const char *, int &);
extern const Element * get_element_from_symbol (const const_IWSubstring &, int &);
extern const Element * get_element_from_symbol (const IWString &, int &);
extern const Element * get_element_from_symbol (char);

extern const Element * get_element_from_symbol_no_case_conversion (const char * s, int nchars);
extern const Element * get_element_from_symbol_no_case_conversion (const char * s);
extern const Element * get_element_from_symbol_no_case_conversion (const const_IWSubstring & s);
extern const Element * get_element_from_symbol_no_case_conversion (const IWString & s);

extern const Element * get_element_from_atomic_number (atomic_number_t);

extern const Element * create_element_with_symbol (const char *, int);
extern const Element * create_element_with_symbol (const char *);
extern const Element * create_element_with_symbol (const IWString &);
extern const Element * create_element_with_symbol (const const_IWSubstring &);

extern int element_from_smiles_string (const char * smiles, int nchars, const Element * & result);
extern int element_from_smarts_string (const char * smiles, int nchars, const Element * & result);

extern int element_from_long_smiles_string (const char * smiles, int nchars, const Element * & result);

#define OK_ELEMENT(e) ( (NULL != (e)) && (e)->ok () )

extern int set_auto_create_new_elements (int);
extern int auto_create_new_elements ();

extern int symbol_for_atomic_symbol_hash_value (int, IWString &);

class Command_Line;

extern int process_elements    (const Command_Line &, int = 0, char = 'E');

extern void set_include_isotopes_in_smiles (int);

extern int display_standard_element_options (std::ostream &);

/*
  Daylight insists that explicit Hydrogens have square brackets. We
  can make that optional
*/

extern void set_explicit_hydrogens_need_square_brackets_in_smiles (int);

extern int print_element_hash_table (std::ostream & os);

extern void set_atomic_symbols_can_have_arbitrary_length (int s);
extern int  atomic_symbols_can_have_arbitrary_length ();

extern void set_display_strange_chemistry_messages(int);

extern void de_allocate_periodic_table();

/*
  Mar 2010, extend the meaning of D and T to input formats other than mdl
*/

extern int interpret_d_as_deuterium ();
extern int interpret_t_as_tritium ();

extern void reset_element_file_scope_variables ();

#endif
