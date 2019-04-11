#include <limits>
#include <stdlib.h>
#include <stdint.h>
#ifdef _WIN32_
#include <stdafx.h>
#endif

#include <nmmintrin.h>

#include "misc.h"

#include "gfp_standard.h"

extern float * precomputed_ratio;

void
GFP_Standard::build_molecular_properties (const Molecular_Properties_Integer & mpr)
{
  copy_vector(_molecular_properties, mpr.rawdata(), 8);
}

// 2017: it is safer to set to max unsigned int, than to -1, also is
// more logical to compare with, rather than < 0 - it is an unsigned variable
// after all.  The previous way failed the gc3tk tests when build with C++-11
static const unsigned int max_uint = std::numeric_limits<unsigned int>::max();
static unsigned int bits_in_mk = max_uint;

void
set_bits_in_mk (int s)
{
  bits_in_mk = s;

  return;
}

void
GFP_Standard::build_mk (IWDYFP & fp)
{
  if (bits_in_mk == max_uint)
  {
    bits_in_mk = fp.nbits();
    if (bits_in_mk > (sizeof(_mk) * IW_BITS_PER_BYTE) || 0 != bits_in_mk % 32)
    {
      cerr << "GFP_Standard::build_mk:maccskeys must be multiple of 32 bits, " << bits_in_mk << " invalid\n";
      abort();
    }
  }
  else if (fp.nbits() != static_cast<int>(bits_in_mk))
  {
    cerr << "GFP_Standard::build_mk:expecting " << bits_in_mk << " bits, got " << fp.nbits() << " in input\n";
    abort();
  }

  _nset_mk = fp.nset();
  memcpy(_mk, fp.bits(), bits_in_mk / 8);

  for (int i = bits_in_mk / 8; i < 32; i++)
  {
    _mk[i] = static_cast<unsigned char>(0);
  }

  return;
}

void
GFP_Standard::build_mk2 (IWDYFP & fp)
{
  if (fp.nbits() != static_cast<decltype(fp.nbits())>(bits_in_mk))
  {
    cerr << "GFP_Standard::build_mk2:expecting " << bits_in_mk << " bits, got " << fp.nbits() << " in input\n";
    abort();
  }

  _nset_mk2 = fp.nset();
  memcpy(_mk2, fp.bits(), bits_in_mk / 8);

  for (int i = bits_in_mk / 8; i < 32; i++)
  {
    _mk2[i] = static_cast<unsigned char>(0);
  }

  return;
}

/*
  Make sure the array is zero'd before calling this
*/

static void
bits_from_array (const int * b,
                 int n,
                 unsigned char * destination,
                 int & nset)
{
  nset = 0;

  for (auto i = 0; i < n/IW_BITS_PER_BYTE; ++i)
  {
    for (auto j = 0; j < IW_BITS_PER_BYTE; ++j)
    {
      if (*b)
      {
        destination[i] |= one_bit_8[j];
        nset++;
      }

      b++;
    }
  }

  return;
}

void
GFP_Standard::build_iwfp(const unsigned char * b,
                         int nset)    // not last param is bits set
{
  _nset_iw = nset;
  memcpy(_iw, b, sizeof(_iw));

  return;
}

template <typename T>
void
GFP_Standard::build_iwfp(const T * b, const int nset)
{
  _nset_iw = nset;

  std::fill_n(_iw, sizeof(_iw), 0);

  for (unsigned int i = 0; i < sizeof(_iw); ++i)
  {
    unsigned char u = _iw[i];

    for (int j = 0; j < IW_BITS_PER_BYTE; ++j)
    {
      if (*b)
        u |= one_bit_8[j];
      b++;
    }
    _iw[i] = u;
  }

  return;
}

template void GFP_Standard::build_iwfp(const int *, const int);
template void GFP_Standard::build_iwfp(const unsigned int *, const int);

void 
GFP_Standard::build_mk (const int * b)
{
  assert (bits_in_mk > 0);

  _nset_mk2 = 0;

  for (size_t i = 0; i < sizeof(_mk); ++i)
  {
    _mk[i] = static_cast<unsigned char>(0);
  }

  bits_from_array(b, bits_in_mk, _mk, _nset_mk);

  return;
}

void 
GFP_Standard::build_mk2 (const int * b)
{
  assert (bits_in_mk > 0);

  _nset_mk2 = 0;

  for (size_t i = 0; i < sizeof(_mk2); ++i)
  {
    _mk2[i] = static_cast<unsigned char>(0);
  }

  bits_from_array(b, bits_in_mk, _mk2, _nset_mk2);

  return;
}

void
GFP_Standard::build_iw  (IWDYFP & fp)
{
  _nset_iw = fp.nset();
  memcpy(_iw, fp.bits(), 256);

  return;
}

/*
  Code from 
 Written by Imran S. Haque (ihaque@cs.stanford.edu)
*/

static inline int
popcount_2fp (const unsigned* bufA,const unsigned* bufB,const int nwords)
{
    int count = 0;
    assert(nwords % 8 == 0);
     
#if defined(__x86_64__)
    int nquads = nwords/2;
    uint64_t* a64 = (uint64_t*)bufA;
    uint64_t* b64 = (uint64_t*)bufB;
    for (int i = 0; i < nquads; i += 4) {
        count +=  _mm_popcnt_u64(a64[i]&b64[i])     + _mm_popcnt_u64(a64[i+1]&b64[i+1])
                + _mm_popcnt_u64(a64[i+2]&b64[i+2]) + _mm_popcnt_u64(a64[i+3]&b64[i+3]);
    }
#else
    const uint32_t * a32 = (const uint32_t *)bufA;
    const uint32_t * b32 = (const uint32_t *)bufB;
    for (int i = 0; i < nwords; i += 4)
    {
      count +=  _mm_popcnt_u32(a32[i]&b32[i])     + _mm_popcnt_u32(a32[i+1]&b32[i+1]) +
                _mm_popcnt_u32(a32[i+2]&b32[i+2]) + _mm_popcnt_u32(a32[i+3]&b32[i+3]);
    }
#endif
    return count;
}

static inline int
popcount (const unsigned char * b, const int nwords)
{
  const int nquads = nwords / 2;
  int count = 0;

  const uint64_t * b64 = reinterpret_cast<const uint64_t *>(b);

  for (int i = 0; i < nquads; i += 4)
  {
    count += _mm_popcnt_u64(b64[i]) + _mm_popcnt_u64(b64[i+1]) + _mm_popcnt_u64(b64[i+2]) + _mm_popcnt_u64(b64[i+3]);
  }

  return count;
}

float
GFP_Standard::tanimoto (const GFP_Standard & rhs) const
{
  float rc = static_cast<float>(0.0);
  
  for (int i = 0; i < 8; ++i)
  {
    const int j = _molecular_properties[i] * 256 + rhs._molecular_properties[i];
    rc += precomputed_ratio[j];
  }

//cerr << "Properties " << rc << " _nset_iw " << _nset_iw << ' ' << rhs._nset_iw << " mk " << _nset_mk << ' ' << rhs._nset_mk << " mk2 " << _nset_mk2 << ' ' << rhs._nset_mk2 << endl;

  if (_nset_iw || rhs._nset_iw)
  {
    int bic = popcount_2fp((const unsigned *) _iw, (const unsigned *)rhs._iw, 64);
    rc += (static_cast<float>(bic) / static_cast<float>(_nset_iw + rhs._nset_iw - bic));
    assert (rc <= 2.0f);
  }

//cerr << "AFter iw " << rc << endl;

  if (_nset_mk || rhs._nset_mk)
  {
    int bic = popcount_2fp((const unsigned *) _mk, (const unsigned *) rhs._mk, 8);

    rc += (static_cast<float>(bic) / static_cast<float>(_nset_mk + rhs._nset_mk - bic));
    assert (rc <= 3.0f);
  }

  if (_nset_mk2 || rhs._nset_mk2)
  {
    int bic = popcount_2fp((const unsigned *) _mk2, (const unsigned *) rhs._mk2, 8);

    rc += (static_cast<float>(bic) / static_cast<float>(_nset_mk2 + rhs._nset_mk2 - bic));
    assert (rc <= 4.0f);
  }

  assert (rc <= 4.0f);

  return rc * 0.25f;
}

void
GFP_Standard::tanimoto_distance_2(GFP_Standard const & fp1,
                                  GFP_Standard const & fp2,
                                  float* d) const
{
  abort();
  float rc1 = static_cast<float>(0.0);
  float rc2 = static_cast<float>(0.0);
  
  int j;
  for (int i = 0; i < 8; i++)
  {
    j = _molecular_properties[i] * 256 + fp1._molecular_properties[i];
    rc1 += precomputed_ratio[j];

    j = _molecular_properties[i] * 256 + fp2._molecular_properties[i];
    rc2 += precomputed_ratio[j];
  }

  if (_nset_iw || fp1._nset_iw)
  {
    int bic = popcount_2fp((const unsigned *) _iw, (const unsigned *)fp1._iw, 64);
    rc1 += (static_cast<float>(bic) / static_cast<float>(_nset_iw + fp1._nset_iw - bic));
  }

  if (_nset_iw || fp2._nset_iw)
  {
    int bic = popcount_2fp((const unsigned *) _iw, (const unsigned *)fp2._iw, 64);
    rc2 += (static_cast<float>(bic) / static_cast<float>(_nset_iw + fp2._nset_iw - bic));
  }

  if (_nset_mk || fp1._nset_mk)
  {
    int bic = popcount_2fp((const unsigned *) _mk, (const unsigned *) fp1._mk, 8);
    rc1 += (static_cast<float>(bic) / static_cast<float>(_nset_mk + fp1._nset_mk - bic));
  }

  if (_nset_mk || fp2._nset_mk)
  {
    int bic = popcount_2fp((const unsigned *) _mk, (const unsigned *) fp2._mk, 8);
    rc2 += (static_cast<float>(bic) / static_cast<float>(_nset_mk + fp2._nset_mk - bic));
  }

  if (_nset_mk2 || fp1._nset_mk2)
  {
    int bic = popcount_2fp((const unsigned *) _mk2, (const unsigned *) fp1._mk2, 8);
    rc1 += (static_cast<float>(bic) / static_cast<float>(_nset_mk2 + fp1._nset_mk2 - bic));
  }

  if (_nset_mk2 || fp2._nset_mk2)
  {
    int bic = popcount_2fp((const unsigned *) _mk2, (const unsigned *) fp2._mk2, 8);
    rc2 += (static_cast<float>(bic) / static_cast<float>(_nset_mk2 + fp2._nset_mk2 - bic));
  }

  d[0] = rc1 * 0.25f;
  d[1] = rc2 * 0.25f;

  return;
}

int
standard_fingerprints_present ()
{
  if (3 != number_fingerprints())
  {
    cerr << "standard_fingerprints_present::must be 3 fingerprints in the file, -IW -MK -MK2 and -MPR, found " << number_fingerprints() << " fingerprints\n";
    return 0;
  }

  int seen_fpiw = 0;
  int seen_fpmk = 0;
  int seen_fpmk2 = 0;
  for (auto i = 0; i < number_fingerprints(); ++i)
  {
    const auto & tag = fixed_fingerprint_tag(i);

    if (tag.starts_with("FPIW"))
      seen_fpiw++;
    else if (tag.starts_with("FPMK2"))
      seen_fpmk2++;
    else if (tag.starts_with("FPMK"))
      seen_fpmk++;
    else
    {
      cerr << "standard_fingerprints_present::unrecognised fingerprint '" << tag << "'\n";
      return 0;
    }
  }

  if (1 != seen_fpiw || 1 != seen_fpmk || 1 != seen_fpmk2)
  {
    cerr << "standard_fingerprints_present::not all fingerprints present, only processes -MPR -IW -MK -MK2\n";
    return 0;
  }

  const IWString & mpr = property_tag();

  if (! mpr.starts_with("MPR"))
  {
    cerr << "standard_fingerprints_present::property tag invalid '" << mpr << "'\n";
    return 0;
  }

  return 1;
}
