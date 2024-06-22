#ifndef FOUNDATIONAL_IWMISC_SPARSE_FP_CREATOR_H_
#define FOUNDATIONAL_IWMISC_SPARSE_FP_CREATOR_H_

#include <algorithm>
#include <memory>
#include <ostream>

#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwbits/iwbits.h"

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

    template <typename O>
    int _write_constant_width_fingerprint(unsigned int, int *, const IWString &, O& output) const;

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
    int daylight_ascii_form_with_counts_encoded (IWString & dyascii) const;
    int daylight_ascii_form_with_counts_encoded (const const_IWSubstring & tag, IWString & dyascii) const;    // first argument is the TDT tag
    int append_daylight_ascii_form_with_counts_encoded (const const_IWSubstring & tag, IWString & dyascii) const;    // first argument is the TDT tag
    int create_from_array_of_ints (const int *, int);

    // Return a Daylight encoded string of a dense vector of the bits hashed to `nbits`.
    IWString FixedWidthFingerprint(int nbits) const;

    // Write a fixed width binary fingerprint with `nbits` bits to `output`.
    template <typename O>
    int write_constant_width_fingerprint(uint32_t nbits, const IWString & tag, O&  output) const;

    // A Daylight encoded string of the sorted bit numbers.
    IWString BitsWithoutCounts() const;

    // Write as `ncols` space separated values to `output`.
    template <typename O>
    int WriteAsDescriptors(int ncols, O & output) const;

    // Write as `ncols` space separated values to `output`.
    // The method needs an array to work, and it can be provided.
    // Note that it is NOT initialised, so if `count` already contains
    // that will be incremented and written.
    template <typename O>
    int WriteAsDescriptors(int ncols, int* count, O& output) const;

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


#if defined(SPARSE_FP_CREATOR_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

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

template <typename O>
int
Sparse_Fingerprint_Creator::WriteAsDescriptors(int ncols, O & output) const {
  std::unique_ptr<int[]> count = std::make_unique<int[]>(ncols);
  std::fill_n(count.get(), ncols, 0);
  return WriteAsDescriptors(ncols, count.get(), output);
}

template <typename O>
int
Sparse_Fingerprint_Creator::WriteAsDescriptors(int ncols, int* count, O & output) const {
  for (auto& [b, c] : _fp) {
    int col = b % ncols;
    count[col] += c;
  }

  for (int i = 0; i < ncols; ++i) {
    output << ' ' << count[i];
  }

  return 1;
}

template <typename O>
int
Sparse_Fingerprint_Creator::write_constant_width_fingerprint(unsigned int nb,
                                 const IWString & tag,
                                 O& output) const
{
  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(nb);
  std::fill_n(tmp.get(), nb, 0);

  return _write_constant_width_fingerprint(nb, tmp.get(), tag, output);
}

template <typename O>
int
Sparse_Fingerprint_Creator::_write_constant_width_fingerprint(unsigned int nb,
                                 int * tmp,
                                 const IWString & tag,
                                 O& output) const
{
  for (FPHash::const_iterator i = _fp.begin(); i != _fp.end(); i++) {
    unsigned int b = (*i).first;

    b = b % nb;

    tmp[b]++;
  }

  IW_Bits_Base dyfp;

  (void) dyfp.construct_from_array_of_ints(reinterpret_cast<const int *>(tmp), nb);

  IWString dy_ascii;
  dyfp.daylight_ascii_representation_including_nset_info(dy_ascii);

  output << tag << dy_ascii << ">\n";

  return output.good();
}

#endif  // SPARSE_FP_CREATOR_IMPLEMENTATION

#endif  // FOUNDATIONAL_IWMISC_SPARSE_FP_CREATOR_H_
