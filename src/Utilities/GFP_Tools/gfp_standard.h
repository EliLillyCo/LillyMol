#ifndef GFP_STANDARD_H
#define GFP_STANDARD_H

#include <optional>

/*
  We so often deal with MPR MK MK2 IW we want to make a very efficient
  implementation
*/

#include "gfp.h"

inline constexpr int kMkFp = 1;
inline constexpr int kMk2Fp = 2;
inline constexpr int kIwFp = 0;

class GFP_Standard
{
  private:
    int _molecular_properties[8];
    unsigned char _iw[256];
    unsigned char _mk[32];
    unsigned char _mk2[32];

    int _nset_mk;
    int _nset_mk2;
    int _nset_iw;

  // private functions.

    void CheckOkMkNBits(uint32_t nbits);

  public:
    void build_molecular_properties (const Molecular_Properties_Integer &);
    void build_iw  (IWDYFP &);
    void build_mk  (IWDYFP &);   // only non const because it calls nset()
    void build_mk2 (IWDYFP &);

    int DebugPrint(std::ostream& output) const;

    int * molecular_properties () { return _molecular_properties;}

    void build_iwfp(const unsigned char *, int nset);             // note last param is bits set
    template <typename T> void build_iwfp(const T *, int nset);   // note last param is bits set

    // Externally specified maccs keys.
    void build_mk(const int *, uint32_t nbits);
    void build_mk2(const int *, uint32_t nbits);

    int natoms () const { return _molecular_properties[0];}
    int nrings () const { return _molecular_properties[1];}
    int aromatic_atoms () const { return _molecular_properties[4];}

    float tanimoto (const GFP_Standard &) const;
    float tanimoto_distance (const GFP_Standard & rhs) const { return 1.0f - tanimoto(rhs);}
    void tanimoto_distance_2 (const GFP_Standard & rhs, const GFP_Standard & rhs2, float *) const;

    // Only return a result if the value will be <= `must_be_closer_than`.
    std::optional<float> tanimoto_distance_if_less(const GFP_Standard& rhs, float must_be_closer_than) const;
    // Only return a result if the value will be >= `must_be_further_than`.
    std::optional<float> tanimoto_distance_if_greater(const GFP_Standard& rhs, float must_be_further_than) const;
};

extern int can_be_compared (const GFP_Standard & fp1, const GFP_Standard & fp2);
extern int standard_fingerprints_present ();
extern void set_bits_in_mk (int s);

#endif
