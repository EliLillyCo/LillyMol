#ifndef GFP_STANDARD_H
#define GFP_STANDARD_H

/*
  We so often deal with MPR MK MK2 IW we want to make a very efficient
  implementation
*/

#include "gfp.h"

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

  public:
    void build_molecular_properties (const Molecular_Properties_Integer &);
    void build_iw  (IWDYFP &);
    void build_mk  (IWDYFP &);   // only non const because it calls nset()
    void build_mk2 (IWDYFP &);

    int * molecular_properties () { return _molecular_properties;}

    void build_iwfp(const unsigned char *, int nset);             // note last param is bits set
    template <typename T> void build_iwfp(const T *, int nset);   // note last param is bits set

    void build_mk(const int *);
    void build_mk2(const int *);

    int natoms () const { return _molecular_properties[0];}
    int nrings () const { return _molecular_properties[1];}
    int aromatic_atoms () const { return _molecular_properties[4];}

    float tanimoto (const GFP_Standard &) const;
    float tanimoto_distance (const GFP_Standard & rhs) const { return 1.0f - tanimoto(rhs);}
    void tanimoto_distance_2 (const GFP_Standard & rhs, const GFP_Standard & rhs2, float *) const;
};

extern int can_be_compared (const GFP_Standard & fp1, const GFP_Standard & fp2);
extern int standard_fingerprints_present ();
extern void set_bits_in_mk (int s);

#endif
