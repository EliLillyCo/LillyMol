#ifndef IW_SMILES_H
#define IW_SMILES_H 1

#include "molecule.h"
#include "iwrandom.h"

/*
  This contains functions used internally by the smiles routines.
*/


extern boolean include_atom_in_smiles (const Molecule *, atom_number_t);

/*
  Usually, aromaticity is never written. 
  Aromaticity can be written if these functions are used.
*/

extern void set_include_aromaticity_in_smiles (int);    // smiles only
extern int  get_include_aromaticity_in_smiles ();       // smiles only

extern void set_include_cis_trans_in_smiles (int);
extern int  include_cis_trans_in_smiles ();

extern void set_include_chiral_info_in_smiles (int);
extern int  include_chiral_info_in_smiles ();

/*
  Sets include_chiral_info_in_smiles to the constructor's value and
  restores initial value on destruction
*/

class Temporarily_Set_Include_Chiral_Info_in_Smiles
{
  private:
    int _initial_value;

  public:
    Temporarily_Set_Include_Chiral_Info_in_Smiles(int);
    ~Temporarily_Set_Include_Chiral_Info_in_Smiles();
};

extern void set_smiles_random_number_seed (random_number_seed_t);
extern random_number_seed_t set_smiles_random_number_seed_random ();

extern int smiles_process_atom (Molecule *, IWString &, atom_number_t, bond_type_t,
                                 atom_number_t, chiral_type_t = NON_CHIRAL);

class Command_Line;

extern int display_standard_smiles_options (std::ostream &);

extern int process_standard_smiles_options (Command_Line &, int = 0, const char = 'K');

extern void set_smiles_reuse_ring_closure_numbers (int);
extern int  smiles_reuse_ring_closure_numbers ();

extern void set_append_coordinates_after_each_atom (int);

extern void set_smiles_native_ordering (int);

extern int set_datatype_name_for_structure_in_tdt_files (const char *);

extern void set_tdt_append_dataitem_content (int s);

extern int smiles_error_message (const char *, int, int, const char *);

extern int set_default_unique_smiles_aromaticity (int a);

extern void set_unset_implicit_hydrogens_known_if_possible (int s);
extern void set_unset_all_implicit_hydrogens_known_attributes (int s);

extern void set_include_isotopic_information_in_unique_smiles (int s);
extern int  include_isotopic_information_in_unique_smiles ();

extern void set_include_directional_bonding_information_in_unique_smiles (int s);

extern void set_unique_determination_version(const int s);

extern void set_include_implicit_hydrogens_on_aromatic_n_and_p (int);

extern void set_smiles_ring_number_offset (int s);

extern void set_display_smiles_interpretation_error_messages(int s);

extern int display_smiles_interpretation_error_messages();

#define CANON_IMPH_CONSIDER_ALL 86
#define CANON_IMPH_CONSIDER_JUST_HETEROATOMS 73
#define CANON_IMPH_CONSIDER_NONE 61

extern void set_consider_implicit_hydrogens_in_unique_smiles(int s);
extern void set_consider_implicit_hydrogens_known_in_unique_smiles(int s);

// Parse_smiles_token is external because it is also used by the smarts routines

extern int
parse_smiles_token (const char * smiles,
                    int characters_to_process,
                    const Element * &    e,
                    aromaticity_type_t & aromatic,
                    formal_charge_t &    fc,
                    int             &    hcount,
                    int             &    chiral_count,
                    int             &    atomic_mass);
int
parse_smiles_token (const char * smiles,
                    int characters_to_process,
                    const Element * &    e,
                    int & aromatic);

/*
  To avoid passing around a lot of arguments when generating a smiles, we bundle
  all the information into an object
*/

#include "iwrnm.h"

class Smiles_Formation_Info
{
  private:
    int _natoms;

    int * _already_done;

    Ring_Number_Manager _rnm;

    atom_number_t _previous_atom;

    atom_number_t _zatom;

    const int * _include_atom;

    int _write_smiles;    // are we writing smiles or smarts

//  The per-atom create embedding information comes from
//  the Smiles_Information object

    const int * _make_smarts_embedding;

//  Similarly, any user specified per-atom smarts information comes from
//  the Smiles_Information object

    const IWString * _user_specified_atomic_smarts;

  public:
    Smiles_Formation_Info (int na, int nr);
    ~Smiles_Formation_Info ();

    void set_already_done (int * s) { _already_done = s;}
    int * already_done () { return _already_done;}

    Ring_Number_Manager & rnm () { return _rnm;}

    int ok () const;

    void set_zatom (atom_number_t s) { _previous_atom = _zatom; _zatom = s;}
    void set_previous_atom (atom_number_t s) { _previous_atom = s;}
    void set_prev_and_zatom(const atom_number_t p, const atom_number_t z) {_previous_atom = p; _zatom = z;}

    atom_number_t zatom () { return _zatom;}
    atom_number_t previous_atom () { return _previous_atom;}

    const int * include_atom () const { return _include_atom;}

    void set_include_atom (const int * s) { _include_atom = s;}

    void set_write_smiles (int s) { _write_smiles = s;}
    int  write_smiles () const { return _write_smiles;}

    void set_make_smarts_embedding (const int * s) { _make_smarts_embedding = s;}
    int make_smarts_embedding (atom_number_t) const;

    void set_user_specified_atomic_smarts(const IWString * const s) { _user_specified_atomic_smarts = s;}
    const IWString * user_specified_atomic_smarts() const { return _user_specified_atomic_smarts;}
};

extern void set_write_smiles_with_smarts_atoms (int);
extern int  write_smiles_with_smarts_atoms ();

extern void set_resolve_unique_smiles_ties_by_geometry(int);

extern void set_consider_isotopes_as_zero_and_non_zero (int s);
extern int  consider_isotopes_as_zero_and_non_zero (int s);

extern int  add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens();
extern void set_add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens(int);

extern void set_write_single_bonds_in_smiles(int);
extern int  write_single_bonds_in_smiles ();

extern void set_include_hcount_in_smiles(int);

extern void set_write_smiles_aromatic_bonds_as_colons(int);
extern int  write_smiles_aromatic_bonds_as_colons ();

extern  void set_display_unusual_hcount_warning_messages(int);
extern  int  display_unusual_hcount_warning_messages ();

extern  void set_include_directionality_in_ring_closure_bonds (int);

extern void reset_smi_file_scope_variables ();
extern void reset_smiles_support_file_scope_variables ();

extern void set_include_atom_map_with_smiles(const int s);
extern int  include_atom_map_with_smiles();

/*
  Used for parsing leading numeric qualifiers
*/

extern int smarts_fetch_numeric (const char * string, int & value, int & qualifier);

/*
  DO NOT CHANGE THE ORDERING - will kill the ok_adjacency matrix
*/

enum SMILESSMARTS_COMPONENT
{
  SSC_ATOM,
  SSC_BOND,
  SSC_RING,
  SSC_OPAREN,
  SSC_CPAREN,
  SSC_DOT
};

enum SS_SPECIFIC
{
  SSC_SINGLE_BOND,
  SSC_DOUBLE_BOND,
  SSC_TRIPLE_BOND,
  SSC_AROMATIC_BOND,
  SSC_DIRECTIONAL_UP,
  SSC_DIRECTIONAL_DOWN,
  SSC_LEADING_NUMERIC,
  SSC_SMARTS_BOND
};

// found in smiles and smarts

#ifdef MAYBE_SOME_DAY_
enum class SMILESSMARTS_COMPONENT {
  ATOM,
  OPAREN,
  CPAREN,
  BOND,
  DOT,
  RING,

// found in smarts

  DOTS,
  SS_SINGLE_BOND,
  SS_DOUBLE_BOND,
  SS_TRIPLE_BOND,
  SS_AROMATIC_BOND,
  SS_DIRECTIONAL_UP_BOND,
  DIRECTIONAL_DOWN_BOND,
  LEADING_NUMERIC,
};
#endif


/*
  the _ISA variable will be set to one of the #defines above.
*/

class SmilesSmarts_Component : public const_IWSubstring
{
  private:
  private:
    int _isa;

//  within the types above, there can be many specific values (atomic number, bond type, ring number...)

    int _specific;

//  atoms get an atom number

    int _atom_number;

  public:
    SmilesSmarts_Component();

    int isa() const { return _isa;}

    int specific() const { return _specific;}

    template <typename T> int debug_print(T & output) const;

    int set(const char * s, const int nchars, const int t);    // last arg will set _isa
    int set_element(const char * s, const int nchars, const int z);
    int set_bond(const char * s, const int nchars, const int btype);
    int set_ring(const char * s, const int nchars, const int r);
    int set_leading_numeric(const char * s, const int nchars);

    void set_atom_number(const int s) { _atom_number = s;}
    int  atom_number() const { return _atom_number;}

    int is_atom() const { return SSC_ATOM == _isa;}
    int is_ring() const { return SSC_RING == _isa;}
    int is_bond() const { return SSC_BOND == _isa;}
    int is_dot() const { return SSC_DOT == _isa && 1 == this->length();}
    int is_dots() const { return SSC_DOT == _isa && this->length() > 1;}
    int is_leading_numeric() const { return SSC_LEADING_NUMERIC == _isa;}
    int is_open_paren() const { return SSC_OPAREN == _isa;}
    int is_close_paren() const { return SSC_CPAREN == _isa;}
};

class Smiles_Text
{
  private:
    IWString _smiles;
    int _ntokens;
    SmilesSmarts_Component * _token;

//  private functions

    void _free_array();
    int _index_of_closing(const int x, const char copen, const char cclose) const;
    int _error(const int pos, const char * msg);
    int _consume_digits(int & cstart);
    int _consume_digit(int & cstart);
    int _consume_leading_numeric(int & i);
    int _handle_dot(int &);
    int _consume_bond(int & i, const int is_smiles);
    int _consume_element(int & i);
    int _consume_smarts_bond(int & i);

    int _ok_adjacency_smarts() const;

  public:
    Smiles_Text();
    ~Smiles_Text();

    template <typename T> int debug_print(T & output) const;

    int build(const char * s, const int nchars, const int is_smiles = 1);

    int ntokens() const { return _ntokens;}

    int establish_atom_numbers();

    const SmilesSmarts_Component & operator[](int s) { return _token[s];}

    int ok_adjacency(const int is_smiles) const;
};

#endif
