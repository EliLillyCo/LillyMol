#ifndef SPARSEFP_H
#define SPARSEFP_H

#include <cstdint>

#include "dyfp.h"
#include "tversky.h"

class Sparse_Fingerprint_Creator;

class Sparse_Fingerprint
{
  protected:
    int _nbits;
    int _nset;
    unsigned int * _bit;
    int * _count;
    uint32_t _sum_squared;
    double _norm;

//  Private functions

    void _recompute_sums();

    similarity_type_t _tanimoto_with_counts(const Sparse_Fingerprint & rhs) const;

    int _counted_form_construct_from_daylight_ascii_representation(const const_IWSubstring & daylight);
    int _counted_form_construct_from_array_of_bits(const void * b);

    int _build_bit(int ndx, const const_IWSubstring & s);

    int _index_for_bit(unsigned int bit) const;

  public:
    Sparse_Fingerprint();
    Sparse_Fingerprint(const Sparse_Fingerprint &);
    ~Sparse_Fingerprint();

    Sparse_Fingerprint & operator = (const Sparse_Fingerprint &);

    int nbits() const { return _nbits;}
    int nset () const { return _nset;}

    int resize(int);    // not sure if this being public is a good idea or not

    // Aborts if the _bit array is out of order.
    int check_sorted() const;

    int debug_print(std::ostream &) const;

    // Like svm_lite: bit:count bit:count ...
    int construct_from_sparse_ascii_representation(const const_IWSubstring &);

    int construct_from_daylight_ascii_representation(const const_IWSubstring &);
    int construct_from_daylight_ascii_representation_uncounted(const const_IWSubstring& daylight);
    int construct_from_tdt_record(const const_IWSubstring &);

    int append_daylight_ascii_form_with_counts_encoded(IWString & s) const;

    similarity_type_t distance(const Sparse_Fingerprint &) const;
    similarity_type_t tanimoto(const Sparse_Fingerprint &) const;
    similarity_type_t tanimoto_with_unlimited_counts(const Sparse_Fingerprint & rhs) const;   // counts not limited to 255
    similarity_type_t continuous_tanimoto(const Sparse_Fingerprint &) const;
    similarity_type_t tversky(const Sparse_Fingerprint &, const Tversky &) const;
    similarity_type_t tversky_distance(const Sparse_Fingerprint &, const Tversky &) const;
    similarity_type_t optimistic_distance(const Sparse_Fingerprint &, const Tversky &) const;
    similarity_type_t tversky_distance01(const Sparse_Fingerprint &, const Tversky &) const;
    similarity_type_t manhattan_distance(const Sparse_Fingerprint & rhs) const;
    similarity_type_t soergel_similarity(const Sparse_Fingerprint & rhs) const;
    similarity_type_t soergel_variant_similarity(const Sparse_Fingerprint & rhs) const;
    similarity_type_t fvb_modified_tanimoto(const Sparse_Fingerprint & rhs) const;
    similarity_type_t tanimoto_binary(const Sparse_Fingerprint &) const;

    // Currently these are both the same.
    similarity_type_t cosine_measure(const Sparse_Fingerprint &) const;
    double cosine_coefficient(const Sparse_Fingerprint & rhs) const;

    uint32_t dot_product(const Sparse_Fingerprint & rhs) const;

    int is_set(unsigned int) const;

    int next_bit_set(int &, unsigned int &, int &) const;     // istart, ibit, icount

    int bits_in_common(const Sparse_Fingerprint &) const;

    // Given any T that supports find(bitnum), reduce `this` to a subset
    // containing only the bits in `hash`.
    // Returns the number of bits removed.
    template <typename T> int ReduceToSubset(const T& hash);

    int truncate_counts_at(int);

    template <typename T> void each_count(T & t);

    void vector_difference(const Sparse_Fingerprint &, const Sparse_Fingerprint &);

    int  count_for_bit(unsigned int b) const;   // will return 0 if not present

    int build_from_sparse_fingerprint_creator(Sparse_Fingerprint_Creator & sfc);

    int remove_bit(unsigned int b);
    int set_count(unsigned int b, int c);

    void * copy_to_contiguous_storage(void *) const;
    void * copy_to_contiguous_storage_gpu(void *) const;
    const void * build_from_contiguous_storage(const void *, int);
};

extern std::ostream & operator << (std::ostream &, const Sparse_Fingerprint &);

extern void set_sparsefp_warn_empty_data(int);
extern void set_continuous_tanimoto_exponent(double s);
// The separator used in sparse representations, bit:count
extern void set_sparse_ascii_separator(char s);

template <typename T>
void
Sparse_Fingerprint::each_count(T & t)
{
  for (int i = 0; i < _nbits; i++)
  {
    t(_count[i]);
  }

  return;
}

// Deliberate decision to not resize the arrays. Probably
// would not make any difference to anything.
template <typename T>
int
Sparse_Fingerprint::ReduceToSubset(const T& hash) {
  int ndx = 0;
  _nset = 0;
  for (int i = 0; i < _nbits; ++i) {
    const auto iter = hash.find(_bit[i]);

    if (iter == hash.cend())
      continue;

    _bit[ndx] = _bit[i];
    _count[ndx] = _count[i];
    _nset += _count[ndx];
    ndx++;
  }

  const int rc = _nbits - ndx;  // Bits removed.

  _nbits = ndx;

  return rc;
}

#endif
