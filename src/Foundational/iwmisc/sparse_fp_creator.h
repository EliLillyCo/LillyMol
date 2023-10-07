#ifndef SPARSE_FP_CREATOR_H
#define SPARSE_FP_CREATOR_H

#include <algorithm>
#include <ostream>

#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "timsort.hpp"

class IW_Unsigned_Int_Hash_Fn
{
  private:
  public:
    size_t operator () (const unsigned int & u) const { return u;}
};

class Sparse_Fingerprint_Creator
{
  public:

    using FPHash = IW_Hash_Map<unsigned int, int>;

  private:
    FPHash _fp;

//  private functions

    int _write_constant_width_fingerprint (unsigned int, int *, const IWString &, std::ostream &) const;
    inline void _convert_to_unsigned_char (unsigned int b, unsigned char & count) const;

    int _daylight_ascii_form_with_counts_encoded (const unsigned int * s,
                                  unsigned int * tmp,
                                  IWString & dyascii) const;

  public:
    Sparse_Fingerprint_Creator ();

    int operator == (const Sparse_Fingerprint_Creator & rhs) const;
    int operator != (const Sparse_Fingerprint_Creator & rhs) const { return ! operator == (rhs);}

    void hit_bit (unsigned int b);
    void hit_bit (unsigned int b, int c);

    unsigned int nbits () const { return static_cast<unsigned int>(_fp.size ());}
    unsigned int nset  () const { return static_cast<unsigned int>(_fp.size ());}

    int debug_print (std::ostream &) const;
    template <typename O> int to_svml(O &) const;

    void clear() { _fp.clear();}

    void copy_bits_to_unsigned_int_array (unsigned int * b, int & ndx) const;
    int  fill_count_array (const unsigned int * b, int * c, int n) const;

    int write_fingerprint (const IWString &, std::ostream &) const;
    int write_fingerprint (const IWString &, IWString_and_File_Descriptor &) const;
    int write_as_descriptors (int, std::ostream &) const;
    int write_constant_width_fingerprint (unsigned int, const IWString &, std::ostream &) const;
    int daylight_ascii_form_with_counts_encoded (IWString & dyascii) const;
    int daylight_ascii_form_with_counts_encoded (const const_IWSubstring & tag, IWString & dyascii) const;    // first argument is the TDT tag
    int append_daylight_ascii_form_with_counts_encoded (const const_IWSubstring & tag, IWString & dyascii) const;    // first argument is the TDT tag
    int create_from_array_of_ints (const int *, int);

    // Return a Daylight encoded string of a dense vector of the bits hashed to `nbits`.
    IWString FixedWidthFingerprint(int nbits) const;

    // A Daylight encoded string of the sorted bit numbers.
    IWString BitsWithoutCounts() const;

    template <typename O> int write_as_feature_count(const char sep, O &) const;
    template <typename O> int write_as_md5_sum(O & output) const;

    int write_in_svml_form (IWString &) const;

    int increment_vector (int *) const;

    const FPHash & bits_found () const { return _fp;}

    int flatten_to_01();
};

/*
  Many programmes (maccskeys, tnass) generate arrays of hit counts. We can
  write them as non-colliding counted fingerprints too
*/

extern int non_colliding_counted_fingerprint_daylight_representation (int nb, const int * k, IWString & dyascii);

extern int words_needed_for_counted_form (int);

extern int form_sparse_fingerprint (int ndx, unsigned int * tmp, IWString & dyascii);

template <typename C, typename O> int unordered_map_to_md5(const std::unordered_map<unsigned int, C> &, O & output);


#ifdef IW_IMPLEMENTATIONS_EXPOSED

#include <memory>

#define SPARSE_FP_MD5
#ifdef SPARSE_FP_MD5
#include "md5.h"
#else
#include "sha2.h"
#endif

template <typename C, typename O> 
int
unordered_map_to_md5(const std::unordered_map<unsigned int, C> & h, O & output)
{
  const auto sz = h.size();

  if (0 == sz)    // not sure what to do 
    return 0;

  const int bytes_needed = sz * sizeof(unsigned int) + sz * sizeof(C);

  unsigned char * x = new unsigned char[bytes_needed]; std::unique_ptr<unsigned char[]> free_x(x);

  unsigned int * u = reinterpret_cast<unsigned int *>(x);

  unsigned int ndx = 0;
  for (auto f : h)
  {
    u[ndx] = f.first;
    ndx++;
  }

  if (ndx != sz)
    std::cerr << "Yipes, ndx " << ndx << " vs size " << sz << '\n';

  assert (sz == static_cast<unsigned int>(ndx));

//std::sort(u, u + ndx);
  gfx::timsort(u, u + ndx);

  C * y = reinterpret_cast<C*>((x + ndx * sizeof(unsigned int)));    // position where we write the counts

  for (unsigned int i = 0; i < ndx; i++)
  {
    const unsigned int b = u[i];

    const auto f = h.find(b);

    *y = static_cast<C>(f->second);
    y++;
  }

  int nbytes;

#ifdef SPARSE_FP_MD5
  ::MD5_CTX ctx;
  ::MD5Init(&ctx);
  ::MD5Update(&ctx, x, bytes_needed);

// Valgrind will report unititialised values within digest, but I cannot see how that would happen

//unsigned char digest[MD5_DIGEST_LENGTH];   // if using system version
  unsigned char out[16];

  ::MD5Final(out, &ctx);

  nbytes = sizeof(out);

#else
  Sha256 sha2;

  sha2.update(x, bytes_needed);

  std::vector<Byte> out;

  sha2.terminate(out);

  nbytes = out.size();
#endif

  static const char * letters = "1234567890abcdef";

  for (int q = 0; q < nbytes; ++q)
  {
    const auto c = out[q];

    output << letters[(c >> 4) & 0x0f] << letters[c & 0x0f];
  }

  return 1;
}

#endif

#endif
