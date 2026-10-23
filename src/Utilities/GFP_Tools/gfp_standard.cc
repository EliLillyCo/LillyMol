#include <stdint.h>
#include <stdlib.h>
#ifdef _WIN32_
#include <stdafx.h>
#endif

#include "Foundational/iwbits/popcount.h"
#include <algorithm>

#include <iostream>
#include <limits>

#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"
#include "gfp_standard.h"
#include "molecular_property_ratio.h"

using std::cerr;

int
GFP_Standard::DebugPrint(std::ostream& output) const {
  output << "GFP_Standard::DebugPrint\n";
  const uint64_t * x = _iw;
  output << "IW\n";
  for (int i = 0; i < 32; ++i) {
    output << i << ' ' << x[i] << '\n';
  }
  x = _mk;
  output << "MK\n";
  for (int i = 0; i < 4; ++i) {
    output << i << ' ' << x[i] << '\n';
  }
  x = _mk2;
  output << "MK2\n";
  for (int i = 0; i < 4; ++i) {
    output << i << ' ' << x[i] << '\n';
  }

  return output.good();
}

void
GFP_Standard::build_molecular_properties(const Molecular_Properties_Integer &mpr) {
  build_molecular_properties(mpr.rawdata(), mpr.nproperties());
}

void
GFP_Standard::build_molecular_properties(const uint8_t* properties, int nproperties) {
  if (nproperties != 8) {
    cerr << "GFP_Standard::build_molecular_properties:expected 8 properties, got "
         << nproperties << '\n';
    abort();
  }

  memcpy(_molecular_properties, properties, 8 * sizeof(uint8_t));
}

void
GFP_Standard::build_molecular_properties(const int* properties, int nproperties) {
  if (nproperties != 8) {
    cerr << "GFP_Standard::build_molecular_properties:expected 8 properties, got "
         << nproperties << '\n';
    abort();
  }

  for (int i = 0; i < 8; ++i) {
    if (properties[i] < 0 || properties[i] > 255) {
      cerr << "GFP_Standard::build_molecular_properties:property " << i
           << " out of range " << properties[i] << '\n';
      abort();
    }
    _molecular_properties[i] = static_cast<uint8_t>(properties[i]);
  }
}

static const unsigned int max_uint = std::numeric_limits<unsigned int>::max();
static uint32_t bits_in_mk = std::numeric_limits<uint32_t>::max();

// A new maccs keys fingerprint is being constructed. Make sure the number
// of bits in that fingerprint is consistent with what is in `bits_in_mk`.
// If it is not initialised, set to `nbits`, otherwise check to see if
// `nbits` is the same as `bits_in_mk`.
// Abort on failure.
void
GFP_Standard::CheckOkMkNBits(uint32_t nbits) {
  // If not initialised, must be OK size.
  if (bits_in_mk == std::numeric_limits<uint32_t>::max()) {
    if (nbits > (sizeof(_mk) * IW_BITS_PER_BYTE) || 0 != nbits % 32) {
      cerr << "GFP_Standard::set_bits_in_mk:maccskeys must be multiple of 32 bits, "
           << nbits << " invalid\n";
      abort();
    }
    bits_in_mk = nbits;
    return;
  } 

  // Already initialised, must be the same.
  if (nbits != bits_in_mk) {
    cerr << "GFP_Standard::build_mk2:expecting " << bits_in_mk << " bits, got "
         << nbits << " in input\n";
    abort();
  }
}

void
GFP_Standard::build_mk(IWDYFP &fp) {
  CheckOkMkNBits(fp.nbits());

  _nset_mk = fp.nset();
  memcpy(_mk, fp.bits(), bits_in_mk / 8);

  // Byte access is valid for the object representation of word storage.
  auto* bytes = reinterpret_cast<unsigned char*>(_mk);
  std::fill(bytes + bits_in_mk / 8, bytes + sizeof(_mk), 0);

  return;
}

void
GFP_Standard::build_mk2(IWDYFP &fp) {
  CheckOkMkNBits(fp.nbits());

  _nset_mk2 = fp.nset();
  memcpy(_mk2, fp.bits(), bits_in_mk / 8);

  // Byte access is valid for the object representation of word storage.
  auto* bytes = reinterpret_cast<unsigned char*>(_mk2);
  std::fill(bytes + bits_in_mk / 8, bytes + sizeof(_mk2), 0);

  return;
}

/*
  Make sure the array is zero'd before calling this
*/

static void
bits_from_array(const int *b, int n, unsigned char *destination, int &nset) {
  nset = 0;

  for (auto i = 0; i < n / IW_BITS_PER_BYTE; ++i) {
    for (auto j = 0; j < IW_BITS_PER_BYTE; ++j) {
      if (*b) {
        destination[i] |= one_bit_8[j];
        nset++;
      }

      b++;
    }
  }

  return;
}

void
GFP_Standard::build_iwfp(const unsigned char *b,
                         int nset)  // not last param is bits set
{
  _nset_iw = nset;
  memcpy(_iw, b, sizeof(_iw));

  return;
}

template <typename T>
void
GFP_Standard::build_iwfp(const T *b, const int nset) {
  _nset_iw = nset;

  std::fill_n(_iw, 32, uint64_t{0});
  auto* bytes = reinterpret_cast<unsigned char*>(_iw);

  for (unsigned int i = 0; i < sizeof(_iw); ++i) {
    unsigned char u = bytes[i];

    for (int j = 0; j < IW_BITS_PER_BYTE; ++j) {
      if (*b) {
        u |= one_bit_8[j];
      }
      b++;
    }
    bytes[i] = u;
  }

  return;
}

template void GFP_Standard::build_iwfp(const int *, const int);
template void GFP_Standard::build_iwfp(const unsigned int *, const int);

void
GFP_Standard::build_mk(const int *b, uint32_t nbits) {
  CheckOkMkNBits(nbits);

  _nset_mk = 0;

  std::fill_n(_mk, 4, uint64_t{0});

  bits_from_array(b, nbits, reinterpret_cast<unsigned char*>(_mk), _nset_mk);

  return;
}

void
GFP_Standard::build_mk2(const int *b, uint32_t nbits) {
  CheckOkMkNBits(nbits);

  _nset_mk2 = 0;

  std::fill_n(_mk2, 4, uint64_t{0});

  bits_from_array(b, nbits, reinterpret_cast<unsigned char*>(_mk2), _nset_mk2);

  return;
}

void
GFP_Standard::build_iw(IWDYFP &fp) {
  _nset_iw = fp.nset();
  memcpy(_iw, fp.bits(), 256);

  return;
}

/*
  Code from
 Written by Imran S. Haque (ihaque@cs.stanford.edu)
*/

static inline int
popcount_word(uint64_t word) {
#ifdef __i386__
  return POPCOUNT(static_cast<uint32_t>(word)) +
         POPCOUNT(static_cast<uint32_t>(word >> 32));
#else
  return POPCOUNT(word);
#endif
}

static inline int
popcount_2fp(const uint64_t* lhs, const uint64_t* rhs, int nwords) {
  int count = 0;
  // Callers specify the historical number of 32-bit words.
  assert(nwords % 8 == 0);
  const int nquads = nwords / 2;
  for (int i = 0; i < nquads; i += 4) {
    count += popcount_word(lhs[i] & rhs[i]) +
             popcount_word(lhs[i + 1] & rhs[i + 1]) +
             popcount_word(lhs[i + 2] & rhs[i + 2]) +
             popcount_word(lhs[i + 3] & rhs[i + 3]);
  }
  return count;
}

//#define DEBUG_GFP_STANDARD_TANIMOTO

float
GFP_Standard::tanimoto(const GFP_Standard &rhs) const {
  float rc = static_cast<float>(0.0);

  int ndiv = 1;
  rc += gfp_internal::molecular_property_ratio_table().similarity(_molecular_properties,
                                                         rhs._molecular_properties, 8);

#ifdef DEBUG_GFP_STANDARD_TANIMOTO
  cerr << "Properties " << rc << " _nset_iw " << _nset_iw << ' ' << rhs._nset_iw <<
          " mk " << _nset_mk << ' ' << rhs._nset_mk << " mk2 " << _nset_mk2 << ' ' << rhs._nset_mk2 << '\n';
#endif

  // Using AVX does not seem to make any difference.
  if (_nset_iw || rhs._nset_iw) {
    int bic = popcount_2fp(_iw, rhs._iw, 64);
//  int bic = iwbits::BitsInCommonAvx((const uint64_t *)_iw, (const uint64_t *)rhs._iw, 32);
    rc += (static_cast<float>(bic) / static_cast<float>(_nset_iw + rhs._nset_iw - bic));
    ++ndiv;
    assert(rc <= 2.0f);
  }

#ifdef DEBUG_GFP_STANDARD_TANIMOTO
  cerr << "AFter iw " << rc << '\n';
#endif

  if (_nset_mk || rhs._nset_mk) {
    int bic = popcount_2fp(_mk, rhs._mk, 8);

    rc += (static_cast<float>(bic) / static_cast<float>(_nset_mk + rhs._nset_mk - bic));
    ++ndiv;
    assert(rc <= 3.0f);
  }
#ifdef DEBUG_GFP_STANDARD_TANIMOTO
  cerr << "After MK " << rc << '\n';
#endif

  if (_nset_mk2 || rhs._nset_mk2) {
    int bic = popcount_2fp(_mk2, rhs._mk2, 8);

    rc += (static_cast<float>(bic) / static_cast<float>(_nset_mk2 + rhs._nset_mk2 - bic));
    assert(rc <= 4.0f);
    ++ndiv;
  }

  assert(rc <= ndiv);

#ifdef DEBUG_GFP_STANDARD_TANIMOTO
  cerr << "tanimoto " << rc << " ndiv " << ndiv << " returning " << static_cast<float>(rc) / static_cast<float>(ndiv) << '\n';
#endif

  return static_cast<float>(rc) / static_cast<float>(ndiv);
}

void
GFP_Standard::tanimoto_distance_2(GFP_Standard const &fp1, GFP_Standard const &fp2,
                                  float *d) const {
  abort();
  float rc1 = static_cast<float>(0.0);
  float rc2 = static_cast<float>(0.0);

  rc1 += gfp_internal::molecular_property_ratio_table().similarity(_molecular_properties,
                                                           fp1._molecular_properties, 8);
  rc2 += gfp_internal::molecular_property_ratio_table().similarity(_molecular_properties,
                                                           fp2._molecular_properties, 8);

  int ndiv1 = 1;
  int ndiv2 = 1;

  if (_nset_iw || fp1._nset_iw) {
    int bic = popcount_2fp(_iw, fp1._iw, 64);
    rc1 += (static_cast<float>(bic) / static_cast<float>(_nset_iw + fp1._nset_iw - bic));
    ++ndiv1;
  }

  if (_nset_iw || fp2._nset_iw) {
    int bic = popcount_2fp(_iw, fp2._iw, 64);
    rc2 += (static_cast<float>(bic) / static_cast<float>(_nset_iw + fp2._nset_iw - bic));
    ++ndiv2;
  }

  if (_nset_mk || fp1._nset_mk) {
    int bic = popcount_2fp(_mk, fp1._mk, 8);
    rc1 += (static_cast<float>(bic) / static_cast<float>(_nset_mk + fp1._nset_mk - bic));
    ++ndiv1;
  }

  if (_nset_mk || fp2._nset_mk) {
    int bic = popcount_2fp(_mk, fp2._mk, 8);
    rc2 += (static_cast<float>(bic) / static_cast<float>(_nset_mk + fp2._nset_mk - bic));
    ++ndiv2;
  }

  if (_nset_mk2 || fp1._nset_mk2) {
    int bic = popcount_2fp(_mk2, fp1._mk2, 8);
    rc1 +=
        (static_cast<float>(bic) / static_cast<float>(_nset_mk2 + fp1._nset_mk2 - bic));
    ++ndiv1;
  }

  if (_nset_mk2 || fp2._nset_mk2) {
    int bic = popcount_2fp(_mk2, fp2._mk2, 8);
    rc2 +=
        (static_cast<float>(bic) / static_cast<float>(_nset_mk2 + fp2._nset_mk2 - bic));
    ++ndiv2;
  }

  d[0] = rc1 / static_cast<float>(ndiv1);
  d[1] = rc2 / static_cast<float>(ndiv2);

  return;
}

std::optional<float>
GFP_Standard::tanimoto_distance_if_less(const GFP_Standard &rhs,
                                        float must_be_closer_than) const {
  // The calculations are done in similarity space, so convert one time.
  const float similarity_needed = 1.0 - must_be_closer_than;

  float rc = static_cast<float>(0.0);

  rc += gfp_internal::molecular_property_ratio_table().similarity(_molecular_properties,
                                                         rhs._molecular_properties, 8);

  if ((rc + 3.0f) / 4.0f < similarity_needed) {
    return std::nullopt;
  }

  if (_nset_mk2 || rhs._nset_mk2) {
    int bic = popcount_2fp(_mk2, rhs._mk2, 8);
    rc += iwmisc::Fraction<float>(bic, _nset_mk2 + rhs._nset_mk2 - bic);
    if ((rc + 2.0) / 4.0 < similarity_needed) {
      return std::nullopt;
    }
  }

  if (_nset_mk || rhs._nset_mk) {
    int bic = popcount_2fp(_mk, rhs._mk, 8);
    rc += iwmisc::Fraction<float>(bic, _nset_mk + rhs._nset_mk - bic);
    if ((rc + 1.0) / 3.0 < similarity_needed) {
      return std::nullopt;
    }
  }

  int bic = popcount_2fp(_iw, rhs._iw, 64);
  rc += iwmisc::Fraction<float>(bic, _nset_iw + rhs._nset_iw - bic);

  rc = rc * 0.25;

  if (rc < similarity_needed) {
    return std::nullopt;
  }

  return 1.0f - rc;
}

int
standard_fingerprints_present() {
  if (3 != number_fingerprints()) {
    cerr << "standard_fingerprints_present::must be 3 fingerprints in the file, -IW -MK "
            "-MK2 and -MPR, found "
         << number_fingerprints() << " fingerprints\n";
    return 0;
  }

  int seen_fpiw = 0;
  int seen_fpmk = 0;
  int seen_fpmk2 = 0;
  for (auto i = 0; i < number_fingerprints(); ++i) {
    const auto &tag = fixed_fingerprint_tag(i);

    if (tag.starts_with("FPIW")) {
      seen_fpiw++;
    } else if (tag.starts_with("FPMK2")) {
      seen_fpmk2++;
    } else if (tag.starts_with("FPMK")) {
      seen_fpmk++;
    } else {
      cerr << "standard_fingerprints_present::unrecognised fingerprint '" << tag << "'\n";
      return 0;
    }
  }

  if (1 != seen_fpiw || 1 != seen_fpmk || 1 != seen_fpmk2) {
    cerr << "standard_fingerprints_present::not all fingerprints present, only processes "
            "-MPR -IW -MK -MK2\n";
    return 0;
  }

  const IWString &mpr = property_tag();

  if (!mpr.starts_with("MPR")) {
    cerr << "standard_fingerprints_present::property tag invalid '" << mpr << "'\n";
    return 0;
  }

  return 1;
}
