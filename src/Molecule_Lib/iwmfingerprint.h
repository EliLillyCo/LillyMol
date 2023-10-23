#ifndef MOLECULE_LIB_IWMFINGERPRINT_H_
#define MOLECULE_LIB_IWMFINGERPRINT_H_

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwstring/iwstring.h"

#include "atom_typing.h"

class Molecule;
class Possible_Fingerprint_Addition;
class MFingerprint;

//#define COUNT_TIMES_ATOM_IN_PATH

class IWMFingerprint : public IW_Bits_Base
{
  friend
    int operator == (const IWMFingerprint &, const IWMFingerprint &);

  private:

//  during generation of the fingerprint we keep the bits in various vectors.

    int * _bvector;

    int * _auxiliary_bvector;

    int _min_heteroatoms_at_path_ends;

    Atom_Typing_Specification _atom_typing_specification;

#ifdef COUNT_TIMES_ATOM_IN_PATH
    int _cycles_adjust_times_used;
#endif

//  private functions

    void _default_values ();
    int _construct_fingerprint (Molecule &, Atom_Typing_Specification & ats, const int * = nullptr, const int * = nullptr);

#ifdef COUNT_TIMES_ATOM_IN_PATH
    int _extra_bits_to_equalise_times_in_path_3(Molecule & m, MFingerprint & mfp, const int * include_these_atoms, const int niter);
    int _extra_bits_to_equalise_times_in_path_4(Molecule & m, MFingerprint & mfp, const int * include_these_atoms, const int niter);
    int _all_3_paths_with_centre(const Molecule & m, const atom_number_t zatom, const int * include_these_atoms, const int * times_in_path,
                                         resizable_array_p<Possible_Fingerprint_Addition> & possible_additions) const;
    int _all_4_paths_from_bond(const Molecule & m, const Bond * b, const int * include_these_atoms, const int * times_in_path,
                                       resizable_array_p<Possible_Fingerprint_Addition> & possible_additions) const;
#endif

  public:
    IWMFingerprint ();
    IWMFingerprint (int);
    ~IWMFingerprint ();

    int nset () const;

    int initialise_atom_typing_specification (const const_IWSubstring &);

    int construct_fingerprint (Molecule &);
    int construct_fingerprint (Molecule &, Atom_Typing_Specification & ats, const int * inc);
    int construct_fingerprint (Molecule &, const int * atype, const int * inc);
    int construct_fingerprint (Molecule &, const int * atype);

#ifdef COUNT_TIMES_ATOM_IN_PATH
    int construct_fingerprint_label_by_times_in_path (Molecule & m);
    void set_cycles_adjust_times_used(const int s) { _cycles_adjust_times_used = s;}
#endif

    int delete_vector_representation ();     // to save space

    const int * vector () const { return _bvector;}

    int level2_fingerprint (int *) const;

    int write_as_zero_and_one (std::ostream &) const;
    int write_as_zero_and_one (IWString &) const;

    int write_count (std::ostream &) const;
    int write_count (IWString &) const;

    void set_min_heteroatoms_at_path_ends (int);

    int  truncate_to_max_hits (int);
};

extern void set_iwmfingerprint_nbits (int);
extern int  iwmfingerprint_nbits ();

extern int  set_min_path_length (int);
extern int  set_max_path_length (int);
extern void set_omit_ch2 (int);
extern void set_include_formal_charge_bits (int);
extern void set_do_atom_pair_bits (int);
extern void set_isotopic_paths (int);
extern void set_unsaturated_includes_aromatic (int);
extern void set_maximal_daylight_compatibility(int);
extern void set_formally_charged_atoms_lose_identity (int);
extern void set_topotorsion_atom_types (int);
extern void set_bits_for_hydrogen_attachments (int);
extern void set_include_bits_for_rings (int);
extern void set_include_unsaturation_in_atom_type (int);

/*
  The value we set for max_path_length_isotopic_bits is the max path
  length for which we want isotopic bits set
*/

extern void set_max_path_length_isotopic_bits (int);

extern int operator == (const IWMFingerprint &, const IWMFingerprint &);
extern int parse_path_length_options (Command_Line & cl, char cflag, int verbose);

extern int parse_misc_fingerprint_options (Command_Line &, char flag, int verbose);
extern int display_misc_fingerprint_options (std::ostream &, char);

/*
  When we create "non-colliding" fingerprints, we just use a very wide 
  creation vector. 
*/

extern void set_form_bit_vector_during_fingerprint_creation(int);

extern void set_only_bits_preserving_substructure_perception(int s);

#endif  // MOLECULE_LIB_IWMFINGERPRINT_H_
