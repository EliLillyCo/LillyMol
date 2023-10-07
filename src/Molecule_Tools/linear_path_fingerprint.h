#ifndef LINEAR_FP_H
#define LINEAR_FP_H

//#define COMPILE_WITH_WATCH_BITS

#include <unordered_map>

#ifdef COMPILE_WITH_WATCH_BITS
#include <unordered_set>
#endif

class Molecule;
class Command_Line;
class IWDigits;
#include "Molecule_Tools/ct_cache.h"

namespace LFP 
{

class Linear_Fingerprint_Defaults
{
  public:
    int _min_path_length;
    int _max_path_length;
    int _min_heteroatoms_at_path_ends;
    int _only_bits_preserving_substructure_perception;
    int _formally_charged_atoms_lose_identity;
    int _omit_ch2;
    int _set_formal_charge_bits;
    int _isotopic_paths;
    int _max_path_length_isotopic_bits;
    int _bits_for_hydrogen_attachments;
    int _place_bits_for_rings;
    int _bits_for_atoms_in_multiple_rings;
    int _terminal_path_bits;
    int _maximal_daylight_compatibility;
    int _differentiate_ring_and_chain_bonds;

  public:
    Linear_Fingerprint_Defaults();

    int build (Command_Line & cl, char flag, const int verbose);

    int min_path_length () const { return _min_path_length;}
    int max_path_length () const { return _max_path_length;}

    int bits_for_hydrogen_attachments () const { return _bits_for_hydrogen_attachments;}
    int max_path_length_isotopic_bits () const { return _max_path_length_isotopic_bits;}

    int bits_for_atoms_in_multiple_rings () const { return _bits_for_atoms_in_multiple_rings;}
    int only_bits_preserving_substructure_perception () const { return _only_bits_preserving_substructure_perception;}
    int set_formal_charge_bits () const { return _set_formal_charge_bits;}
    int maximal_daylight_compatibility () const { return _maximal_daylight_compatibility;}
    int isotopic_paths () const { return _isotopic_paths;}
    int omit_ch2 () const { return _omit_ch2;}
    int place_bits_for_rings () const { return _place_bits_for_rings;}
    int formally_charged_atoms_lose_identity () const { return _formally_charged_atoms_lose_identity;}
    int terminal_path_bits () const { return _terminal_path_bits;}
    int min_heteroatoms_at_path_ends () const { return _min_heteroatoms_at_path_ends;}
    int differentiate_ring_and_chain_bonds() const { return _differentiate_ring_and_chain_bonds;}

    void set_min_path_length(int s) { _min_path_length = s;}
    void set_max_path_length(int s) { _max_path_length = s;}
    void set_differentiate_ring_and_chain_bonds(const int s) { _differentiate_ring_and_chain_bonds = s;}

    void set_max_path_length_isotopic_bits(int s) { _max_path_length_isotopic_bits = s;}

    void set_place_bits_for_rings(int s) { _place_bits_for_rings = s;}

    int numeric_bond_code(const Bond * b) const;
};

class Linear_Fingerprint_Creator
{
  private:
    const Linear_Fingerprint_Defaults * _lfpd;   // must be kept in scope somewhere else

  public:
    Linear_Fingerprint_Creator();
    Linear_Fingerprint_Creator(const Linear_Fingerprint_Defaults *);

    void set_defaults (const Linear_Fingerprint_Defaults * s) { _lfpd = s;}

    template <typename T, typename B> int process (Molecule & m, const T * atype, const int * include_these_atoms, B & bits);
    template <typename T, typename B> int process (Molecule & m, const T * atype, B & bits);

    template <typename C, typename T, typename B> int process (const CT_Cache<C> &, const T * atype, const int * include_these_atoms, B & bits);
};

typedef unsigned int fwf_bits_t;

class Fixed_Width_Fingerprint_Bits
{
  private:
    const int _nbits;
    fwf_bits_t * _bits;

  public:
    Fixed_Width_Fingerprint_Bits();
    Fixed_Width_Fingerprint_Bits(int);
    ~Fixed_Width_Fingerprint_Bits();

    void clear();

    int truncate_to_max_hits(const int);
    int identify_differences (const Fixed_Width_Fingerprint_Bits & rhs, std::ostream &) const;
    int append_fingerprint (const IWString & tag, IWString &) const;

    fwf_bits_t & operator[](const unsigned int);
    fwf_bits_t & operator[](const unsigned int) const;

    int nbits() const;
    int nset() const;

    unsigned int bit_number (unsigned int b) { return b % _nbits;}

    const fwf_bits_t * cbegin() const { return _bits;}
    fwf_bits_t * begin() const { return _bits;}
    const fwf_bits_t * cend() const { return _bits + _nbits;}
    fwf_bits_t * end() const { return _bits + _nbits;}

    template <typename T> int write_as_zero_and_one (T &) const;
    template <typename T> int write_count (const IWDigits & d, T & output) const;
    template <typename T> int write_as_feature_count(const IWString & sep, T & output) const;
    int write_as_md5_sum(IWString & output) const;
};

class Sparse_Fingerprint_Bits : public std::unordered_map<unsigned int, int>
{
  typedef unsigned int sfp_bits_type;

  private:
  public:
//  void clear();

    int truncate_to_max_hits(const int);
    int identify_differences (const Sparse_Fingerprint_Bits & rhs, std::ostream &) const;

    int nbits() const { return this->size();}
//  int nset() const { return this->size();}
    int nset() const { return std::unordered_map<sfp_bits_type, int>::size();}

    sfp_bits_type bit_number (unsigned int b) { return b;}

    int append_fingerprint (const IWString & tag, IWString &) const;
    template <typename T> int write_as_zero_and_one (T &) const;
    template <typename T> int write_count (const IWDigits & d, T & output) const;
    template <typename T> int write_as_feature_count(const IWString & sep, T & output) const;
    int write_as_md5_sum (IWString & output) const;
};

extern void set_default_linear_path_nbits (const int s);

extern void set_default_min_path_length (const int s);
extern void set_default_max_path_length (const int s);
extern void set_only_bits_preserving_substructure_perception(const int s);
extern void set_default_set_formal_charge_bits(const int s);
extern void set_default_omit_ch2(const int s);
extern void set_default_min_heteroatoms_at_path_ends (const int s);
extern void set_default_bits_for_hydrogen_attachments(int s);
extern void set_default_differentiate_ring_and_chain_bonds(const int s);

extern int parse_misc_fingerprint_options (Command_Line & cl, char flag, int verbose);

#if defined(COMPILE_WITH_WATCH_BITS)
extern void add_bit_to_watch(const unsigned int);
extern int read_bits_to_watch (const IWString & fname);
#endif

}    // end namespace LFP

#endif
