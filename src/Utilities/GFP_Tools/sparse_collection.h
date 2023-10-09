#ifndef SPARSE_COLLECTION_H
#define SPARSE_COLLECTION_H

/*
  When we are dealing with a collection of sparse fingerprints,
  we can convert them to fixed width fingerprints if we have a
  translation table from the arbitrary fingerprint numbers to
  some fixed number
*/

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "gfp.h"

/*
  We need a profile for each sparse fingerprint found in the gfp object
*/

class Sparse_Fingerprint_Collection_Profile
{
  private:

    typedef IW_STL_Hash_Map<unsigned int, unsigned int> xref_t;

    xref_t _xref; 

//  We also keep track of the number of items that have a given bit set.
//  We won't bother converting singletons into the xref array

    extending_resizable_array<int> _count;

//  It can be informative to give into on the average number of bits set in the fingerprints
//  we profile

    Accumulator_Int<int> _nset;

  public:
    Sparse_Fingerprint_Collection_Profile ();

    int build_profile (const Sparse_Fingerprint &);

    int report (std::ostream &) const;

    int nbits () const { return _xref.size ();}

    int remove_singletons (int = 1, int = 0);

    int convert_to_fixed_width (const Sparse_Fingerprint & fpfrom,
                                IWDYFP & fpto, int & nextra) const;

    int convert_to_fixed_width (const Sparse_Fingerprint & fpfrom,
                                Fixed_Size_Counted_Fingerprint_uchar & fpto,
                                int & extra_bits,
                                int & extra_count) const;
    int convert_to_fixed_width (const Sparse_Fingerprint & fpfrom,
                                Fixed_Size_Counted_Fingerprint_uint & fpto,
                                int & extra_bits,
                                int & extra_count) const;
                                        
};

/*
  An IW_General_Fingerprint object can have any number of sparse fingerprints
*/

class IW_General_Fingerprint;

class Set_of_Sparse_Fingerprint_Collection_Profile
{
  private:
    int _number_sparse_fingerprints;

    Sparse_Fingerprint_Collection_Profile * _sfcp;
  public:
    Set_of_Sparse_Fingerprint_Collection_Profile ();
    ~Set_of_Sparse_Fingerprint_Collection_Profile ();

    int report (std::ostream &) const;

    int resize (int);

    int active () const { return _number_sparse_fingerprints;}

    int build_profile (const IW_General_Fingerprint &);

    int finished_profiling (int verbose);

    int convert_to_fixed_width (int which_collection, const Sparse_Fingerprint & fpfrom, IWDYFP & fpto, int & nextra) const {
      return _sfcp[which_collection].convert_to_fixed_width(fpfrom, fpto, nextra);
    }
    int convert_to_fixed_width (int which_collection, const Sparse_Fingerprint & fpfrom, Fixed_Size_Counted_Fingerprint_uchar & fpto, int & extra_bits, int & extra_count) const {
      return _sfcp[which_collection].convert_to_fixed_width(fpfrom, fpto, extra_bits, extra_count);
    }
    int convert_to_fixed_width (int which_collection, const Sparse_Fingerprint & fpfrom, Fixed_Size_Counted_Fingerprint_uint & fpto, int & extra_bits, int & extra_count) const {
      return _sfcp[which_collection].convert_to_fixed_width(fpfrom, fpto, extra_bits, extra_count);
    }
};

class Command_Line;

extern int parse_sparse_to_dense_fingerprint_specifications (Command_Line &, char, int verbose);

#endif
