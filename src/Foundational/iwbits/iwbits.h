#ifndef IW_BITS_BASE_H
#define IW_BITS_BASE_H

#include <iostream>

using std::cerr;
using std::endl;

#define IW_BITS_PER_WORD 32
#define IW_BITS_PER_BYTE 8
#define IW_BYTES_PER_WORD (sizeof (unsigned int) / sizeof (unsigned char))

class IWString;
class const_IWSubstring;
/*
  Byte values in which different bits are turned on
*/

static const unsigned char one_bit_8[8] = {128, 64, 32, 16, 8, 4, 2, 1};

#if defined(__i386__) || defined(__x86_64)
static const unsigned int one_bit_32[] = {
  0x00000080,
  0x00000040,
  0x00000020,
  0x00000010,
  0x00000008,
  0x00000004,
  0x00000002,
  0x00000001,
  0x00008000,
  0x00004000,
  0x00002000,
  0x00001000,
  0x00000800,
  0x00000400,
  0x00000200,
  0x00000100,
  0x00800000,
  0x00400000,
  0x00200000,
  0x00100000,
  0x00080000,
  0x00040000,
  0x00020000,
  0x00010000,
  0x80000000,
  0x40000000,
  0x20000000,
  0x10000000,
  0x08000000,
  0x04000000,
  0x02000000,
  0x01000000
};
#else
static const unsigned int one_bit_32[] = {
  0x80000000,
  0x40000000,
  0x20000000,
  0x10000000,
  0x08000000,
  0x04000000,
  0x02000000,
  0x01000000,
  0x00800000,
  0x00400000,
  0x00200000,
  0x00100000,
  0x00080000,
  0x00040000,
  0x00020000,
  0x00010000,
  0x00008000,
  0x00004000,
  0x00002000,
  0x00001000,
  0x00000800,
  0x00000400,
  0x00000200,
  0x00000100,
  0x00000080,
  0x00000040,
  0x00000020,
  0x00000010,
  0x00000008,
  0x00000004,
  0x00000002,
  0x00000001 };
#endif

static const char hex_char [] = {
  '0',
  '1',
  '2',
  '3',
  '4',
  '5',
  '6',
  '7',
  '8',
  '9',
  'a',
  'b',
  'c',
  'd',
  'e',
  'f'
};

class IW_Bits_Base
{
//friend
//  int operator == (const IW_Bits_Base &, const IW_Bits_Base &);
  friend
    int operator != (const IW_Bits_Base &, const IW_Bits_Base &);
  friend
    IW_Bits_Base
      operator & (const IW_Bits_Base &, const IW_Bits_Base &);
  friend
    IW_Bits_Base &
      operator | (const IW_Bits_Base &, const IW_Bits_Base &);
  friend
    IW_Bits_Base &
      operator ^ (const IW_Bits_Base &, const IW_Bits_Base &);

  protected:
    int _nbits;
    unsigned char * _bits;

//  For bit vectors which are not a multiple of 8 bits,
//  we keep track of any odd bits off the end.
//  _nbits == _whole_bytes * bits_per_byte + _extra_bits;

    int _whole_bytes;
    int _extra_bits;

//  private functions

    void _default_values ();

    int _increase_size_for_bits (int);
    int _allocate_space_for_bits (int);

    int _left_shift (int);
    int _right_shift (int);

    void _copy_bits (const IW_Bits_Base & rhs);

    int _construct_from_tdt_record (const char * ascii,
                                         int nchars,
                                         int check_nset);

  public:
    IW_Bits_Base ();
    IW_Bits_Base (int);
    IW_Bits_Base (const IW_Bits_Base &);
    ~IW_Bits_Base ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int allocate_space_for_bits (int);

    template <typename T> int construct_from_array_of_ints (const T *, int);
    int construct_from_array_of_bits (const unsigned char *, const int);
    int construct_from_array_of_chars (const unsigned char *, int, char = '1', char = '0');
    int construct_from_daylight_ascii_bit_rep (const char *, const int);
    int construct_from_daylight_ascii_representation (const const_IWSubstring &);
    int construct_from_sparse_representation (const const_IWSubstring & bit_representation);
    int construct_from_ascii_01_representation (const char *, int);
    int construct_from_ascii_01_representation_with_spaces (const char * s, int n);
    int construct_from_tdt_record    (const IWString &, int = 0);
    int construct_from_tdt_record    (const const_IWSubstring &, int = 0);
    int construct_from_tdt_record_nset (const const_IWSubstring & tdt_record, int & nset);
    int construct_from_hex (const const_IWSubstring &);

    int hex_form (IWString &) const;

    IW_Bits_Base & operator = (const IW_Bits_Base &);

    int nbits  () const { return _nbits;}
    const unsigned char * bits () const { return _bits;}

    int is_empty () const;    // nothing allocated

    void set (int, int = 1);
    int  is_set (int) const;

    void clear ();     // set all bits to 0
    void set_all ();   // set all bits to 1

    void set_first (int, int);   // set the first N bits to a value
    void set_all_bits (int, int, int);  // set all bits between N1 and N2 to a value

    void flip ();     // invert all bits

    void swap_byte_order ();     // little endian/big endian conversion

    int  nset () const;

    bool any_bits_set () const;   // doesn't count, just checks whole words for non zero values - fast

    int  printon (std::ostream &, const char = '1', const char = '0', int = 0) const;

    int  printon_fast (std::ostream & os, int include_space) const;    // uses '1' and '0' only

    int  append_string_form (IWString &, const char = '1', const char = '0', int = 0) const;
    int  append_string_form_fast (IWString & buffer, int include_space) const;
    int  append_ascii_01_representation (IWString & s) const;

    char *  daylight_ascii_representation (int & nchars) const;
    int daylight_ascii_representation (IWString &) const;
    int daylight_ascii_representation_including_nset_info (IWString &) const;

    int write_daylight_ascii_representation (std::ostream & output, const IWString & tag) const;
    int write_daylight_ascii_representation (std::ostream & output) const;

    int  set_vector (int *) const;
    int  set_vector (int *, int) const;
    int  increment_vector (int *) const;
    int  increment_vector (int *, int) const;

    int  rearrange (const int * xref, IW_Bits_Base & b2) const;

    void iwand (const IW_Bits_Base &);
    void iwand (const IW_Bits_Base &, int &);

    void iwor  (const IW_Bits_Base &);
    void iwor  (const IW_Bits_Base &, int &);

    void iwxor (const IW_Bits_Base &);
    void iwxor (const IW_Bits_Base &, int &);

    int  bits_in_common (const IW_Bits_Base &) const;

    void unset_bits_in_rhs(const IW_Bits_Base & f2);

    int first_bit () const;

//  This next function used by druglike: produce sum of weights for each bit set.

    float compute_weight (const float *) const;
    int   compute_weight (const int   *) const;

//  When dealing with clusters and such, we can compute a weight which also
//  includes a negative component for missing values

    float compute_weight_inc_missing (const float *, int) const;
    int   compute_weight_inc_missing (const int *, int) const;

    void operator += (const IW_Bits_Base &);

    int operator == ( const IW_Bits_Base &) const;

    int  equal (const IW_Bits_Base *) const;

    int fold (int);

//  With substructure screening we need to know whether or not one fingerprint
//  is a subset of another

    int is_subset (const IW_Bits_Base &) const;

    int shift (int);

    void * copy_to_contiguous_storage (void *) const;
    void * copy_to_contiguous_storage_gpu (void *) const;
    const void * build_from_contiguous_storage (const void *, int);

//  Sequentially fetch all the 1 bits. Pass 0 at first, Returns -1 when done

    int next_on_bit (int &) const;
};

typedef float similarity_type_t;

//extern similarity_type_t tanimoto         (IW_Bits_Base *, IW_Bits_Base *);
extern similarity_type_t fraction_matched (IW_Bits_Base *, IW_Bits_Base *);

extern int print_bits (std::ostream & os, const void * bits, int nbits,
                       const char t = '1', const char f = '0', int = 0);

#ifdef _WIN32
extern int bits_in_common (const unsigned int *, const unsigned int *, int);
#else
extern "C" int bits_in_common (const unsigned int *, const unsigned int *, int);
#endif
extern "C" int count_bits_set (const unsigned char *, int);

extern int write_fixed_size_counted_fingerprint (const int * c,
                                      int n,
                                      int nset,
                                      int bits_in_original_fingerprint,
                                      std::ostream & output);
extern int append_fixed_size_counted_fingerprint (const int * c,
                                      int n,
                                      int nset,
                                      int bits_in_original_fingerprint,
                                      IWString & output_buffer);

#endif
