#ifndef IW_ATOM_TYPING_H
#define IW_ATOM_TYPING_H

#include "iwbits.h"

#include "charge_assigner.h"
#include "donor_acceptor.h"

#define IWATTYPE_Z 1
#define IWATTYPE_COMPLEX 2
#define IWATTYPE_TT 4
#define IWATTYPE_SYBYL 8
#define IWATTYPE_BASIC 16
#define IWATTYPE_EXPT 32
#define IWATTYPE_HB 64
#define IWATTYPE_SF 128
#define IWATTYPE_NONE 256
#define IWATTYPE_SFX 512
#define IWATTYPE_CH 1024
#define IWATTYPE_CC 2048
#define IWATTYPE_PP 4096
#define IWATTYPE_ZA 8192
#define IWATTYPE_ZP 16384

// Realised that all the choices above are mutually exclusive, so no need to keep
// them as separate bits

#define IWATTYPE_NOX 73

#define DIFFERENTIATE_RINGS 32768
#define PERFORM_SHELL_ITERATION 65536

/*
  Oct 2009. We introduce an entirely flexible atom typing that can come
  from the command line. 
  There will be two things,
    an atom type that indicates this is the type to use
    and an int that holds the components
*/

#define IWATTYPE_USP 131072

// atomic numbers
#define IWATTYPE_USP_Z one_bit_32[0]
// atomic numbers halogens equivalent
#define IWATTYPE_USP_Y one_bit_32[1]
// Hydrogen count
#define IWATTYPE_USP_H one_bit_32[2]
// pi electrons
#define IWATTYPE_USP_P one_bit_32[3]
// aromatic
#define IWATTYPE_USP_A one_bit_32[4]
// ncon
#define IWATTYPE_USP_C one_bit_32[5]
// ring bond count
#define IWATTYPE_USP_R one_bit_32[6]
// no type
#define IWATTYPE_USP_N one_bit_32[7]
// carbon or heteroatom
#define IWATTYPE_USP_E one_bit_32[8]
// unsaturated - for now, includes arom
#define IWATTYPE_USP_U one_bit_32[9] 
// include isotope information in the atom type
#define IWATTYPE_USP_I one_bit_32[10] 
// isolated or fused ring
#define IWATTYPE_USP_F one_bit_32[11]

// Y atom typing, but all aromatic atoms the same
#define IWATTYPE_USP_M one_bit_32[12]

// Same as Y, but all possibly tautomeric aromatic nitrogens get same type

#define IWATTYPE_USP_T one_bit_32[13]

// Concept of centrality of an atom
#define IWATTYPE_USP_X one_bit_32[14]

// Smallest ring size containing the atom
#define IWATTYPE_USP_S one_bit_32[15]
// Largest ring size containing the atom
#define IWATTYPE_USP_L one_bit_32[16]

// the atomic numbers of the connected atoms
#define IWATTYPE_USP_K one_bit_32[17]

// Boolean pi electron
#define IWATTYPE_USP_Q one_bit_32[18]

// formal charge

#define IWATTYPE_USP_O one_bit_32[19]

#define IWATTYPE_PPHORE 262144

// The atomic symbol hash value - useful for non periodic table elements

#define IWATTYPE_USP_G one_bit_32[17]

class Molecule;

class Atom_Typing_Specification
{
  private:
    int _type;
    unsigned int _user_specified_type;
    int _differentiate_rings;
    int _perform_shell_iteration;

    IWString _built_from;

//  For Pharmacaphore types

    Charge_Assigner _charge_assigner;
    Donor_Acceptor_Assigner _donor_acceptor_assigner;
    resizable_array_p<Substructure_Query> _hydrophobe;
    int _assign_other_type;
    int _combine_donor_and_acceptor;

    Molecule_Output_Object _write_isotopically_labelled;

//  private functions

    int _parse_user_specified_type(const const_IWSubstring & p);
    int _append_tag_for_user_specified_type (IWString & tag) const;
    template <typename T> int _assign_user_specified_type    (Molecule & m, T * atype, const int * ncon) const;
    template <typename T> int _ust_assign_atom_types_z_prime_numbers (const Molecule & m, T * atype, int compress_halogens) const;
    template <typename T> int _ust_assign_atom_types_x       (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_t       (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_hcount  (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_pi      (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_arom    (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_ncon    (const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_ncon    (const Molecule & m, T * atype, const int * ncon) const;
    template <typename T> int _ust_assign_atom_types_nrings  (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_none    (const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_isotope (const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_heteroatom (const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_unsaturated (const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_pi_boolean (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_ring_fusion (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_aromatic_all_the_same (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_smallest_ring (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_largest_ring (Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_connected_atoms_hash(const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_atomic_symbol_hash_value (const Molecule & m, T * atype) const;
    template <typename T> int _ust_assign_atom_types_formal_charge(const Molecule & m, T * atype) const;

    template <typename T> void _perform_shell_expansion_v1 (Molecule & m, T *) const;
    template <typename T> void _perform_shell_expansion_v2 (Molecule & m, T * atype) const;
    template <typename T> int _perform_shell_iteration_v2 (Molecule & m, const atom_number_t zatom,
                                const T * atype,
                                int * complete,
                                int need_to_compute_bond_constants,
                                int radius) const;

    int _build_pharmacaphore_specification (const const_IWSubstring &);
    int _build_pharmacaphore_specification (const IWString &, iwstring_data_source &);
    template <typename T> int _assign_atom_types_pharmacaphore (Molecule & m, T *atype);

  public:
    Atom_Typing_Specification();

    int active () const;

    void set_atom_type (int s) { _type = s;}   // should probably check validity
    int  atom_type () const { return _type;}

    void set_user_specified_type (unsigned int);

    int append_to_tag (IWString &) const;
    int string_representation (IWString &) const;

    int build (const const_IWSubstring & s);

    template <typename T> int assign_atom_types (Molecule & m, T * atype, const int * ncon = NULL);   // non const because the pharmacaphore type contains queries

//  ec fingerprints started with Jibo's prime number atomic number types,
//  but other programmes do not use that convention.

    void swap_atomic_number_atom_type_to_atomic_number_prime() {
      if (_type == IWATTYPE_Z)
        _type = IWATTYPE_ZP;
    }
};

/*
  Most programmes want to examine a string and determine the atom typing - this is the
  return code.
  It also appends a standard suffix to an existing tag.
  for example, map would pass NCMAP in TAG and it would come out as something like NCMAPC<
*/

extern int determine_atom_type_and_set_tag(const const_IWSubstring & p, IWString & tag, int verbose = 0);

extern int iwattype_convert_to_string_form (int atype, IWString & s);
extern int determine_atom_type (const IWString & s);

/*
  Some programmes have an array of ncon values available. Use it if possible
*/

template <typename T> int assign_atom_types(Molecule & m, int typing_to_use, T * atype, const int * ncon = NULL);

extern void set_assign_arbitrary_values_to_unclassified_atoms(int s);

extern void set_use_version_2_augmented_atom_algorithm (int s);

#endif

