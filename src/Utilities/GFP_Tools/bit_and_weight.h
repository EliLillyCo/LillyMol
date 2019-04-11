#ifndef BIT_AND_WEIGHT_H
#define BIT_AND_WEIGHT_H

#include <math.h>

#include "iw_stl_hash_map.h"

#include "gfp.h"

/*
  Horrendous HACK here.

  It is so much easier to assume a fixed size for these arrays. Otherwise
  we need to construct these in the presence of the first fingerprint.
*/

#define NBC_FIXED_ARRAY_SIZE 10

class Bit_and_Weight
{
  protected:
    typedef IW_Hash_Map<unsigned int, double> BWHT;

    BWHT _weight;

  public:

    void set_weight(unsigned int b, double w) { _weight[b] = log10(w);}

    void set_count(unsigned int b, double c) { _weight[b] = c;}

    void increment(unsigned int b, double e) { _weight[b] += e;}

    unsigned int size() const { return _weight.size();}

    void increment_weight(unsigned int b, double & rc) const;

    int build_scoring_model (Bit_and_Weight & bw, int population) const;

    int write_scoring_model (const IWString & prefix, IWString_and_File_Descriptor & output) const;

};

class Fixed_and_Sparse_Bit_and_Weight
{
  public:    // public!!!
    int _nfp;
    Bit_and_Weight * _fp;
    int _nsfp;
    Bit_and_Weight * _sparsefp;

  public:
    Fixed_and_Sparse_Bit_and_Weight();
    ~Fixed_and_Sparse_Bit_and_Weight();

    int status(std::ostream &) const;

    int arrays_allocated() const { return _nfp > 0;}

    int initialise(const IW_General_Fingerprint & fp);

    int build_scoring_model (Fixed_and_Sparse_Bit_and_Weight & rhs,
                                        int population) const;

    int write_scoring_model (IWString_and_File_Descriptor & output) const;

    void set_fixed_count(int f, int b, double c) { _fp[f].set_count(b, c);}
    void set_sparse_count(int f, int b, double c) { _sparsefp[f].set_count(b, c);}

    void increment_fixed(int f, int b, double e) { _fp[f].increment(b, e);}
    void increment_sparse(int f, unsigned int b, double e) { _sparsefp[f].increment(b, e);}
};

#endif
