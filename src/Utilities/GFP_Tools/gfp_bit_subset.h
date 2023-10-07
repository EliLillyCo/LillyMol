#ifndef UTILITIES_GFP_TOOLS_GFP_BIT_SUBSET_H
#define UTILITIES_GFP_TOOLS_GFP_BIT_SUBSET_H
#include <iostream>

/*
  These symbols are shared between gfp_to_svm_lite and gfp_bit_subset
*/

#define HEADER_RECORD "# written by gfp_to_svm_lite"
#define COUNT_FIXED "count_fixed"
#define COUNT_SPARSE "count_sparse"
#define PROPERTIES_TAG "PROPERTIES"
#define FIXED_TAG "FIXED"
#define SPARSE_TAG "SPARSE"

#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwbits/iwbits.h"
#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/sparsefp.h"

/*
  To subset a fixed width fingerprint, we just mask out the bits not used
  and whether or not each bit is included or not
*/

class Fixed_Width_Fingerprint_Subset
{
  private:
    int _active;   // maybe all bits are set, don't need to do anything

    IW_Bits_Base  _mask;

  public:
    Fixed_Width_Fingerprint_Subset();

    int do_read (int nb, iwstring_data_source & input);

    int active() const { return _active;}

    int create_subset (IWDYFP & fp) const;

    int nbits() const { return _mask.nbits();}

    int bits_used_in_fixed_width_fingerprint() const;
};

class GFP_Bit_Subset
{
  private:
    int _active;

//  for properties, a flag as to whether each property is included

    int _number_properties;
    int * _properties;

//  for fixed size fingerprints, each gets an array of which bits are in

    int _number_fixed;

    Fixed_Width_Fingerprint_Subset * _fixed;

//  Each sparse fingerprint gets its own hash set

    int _number_sparse;

    IW_Hash_Set<unsigned int> * _sparse;

//  private functions

    int _read_property_mapping (const const_IWSubstring & buffer,
                                        iwstring_data_source & input);
    int _read_fixed_fingerprint (const const_IWSubstring & buffer,
                                         Fixed_Width_Fingerprint_Subset & fp,
                                         iwstring_data_source & input);

    int _do_property_subset (Molecular_Properties_Integer & p) const;
    int _create_sparse_subset (const IW_Hash_Set<unsigned int> & h,
                                       Sparse_Fingerprint & fp) const;

  public:
    GFP_Bit_Subset();
    ~GFP_Bit_Subset();

    int debug_print (std::ostream &) const;

    int active() const { return _active;}

//  Reads files created by gfp_to_svm_lite

    int do_read (const char *);
    int do_read (iwstring_data_source &);

    int number_fixed_width_fingerprints() const { return _number_fixed;}
    int number_sparse_fingerprints() const { return _number_sparse;}

    int number_properties() const { return _number_properties;}
    int bits_in_fixed_width_fingerprint(int i) const { return _fixed[i].nbits();}
    int bits_used_in_fixed_width_fingerprint(int i) const { return _fixed[i].bits_used_in_fixed_width_fingerprint();}
    int bits_in_sparse_fingerprint (int i) const { return _sparse[i].size();}

    int reduce_to_subset (IW_General_Fingerprint &) const;
};

#endif  // UTILITIES_GFP_TOOLS_GFP_BIT_SUBSET_H
