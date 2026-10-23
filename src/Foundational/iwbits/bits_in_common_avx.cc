#include <cstdint>
#include <cstddef>

// immintrin.h is x86 only - on other architectures it is a hard #error, not an
// empty header. The AVX512 functions below are already guarded, so this file
// contributes nothing off x86 and the include must be skipped too. Without this
// the whole library fails to compile on arm64, Apple silicon included.
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#include <immintrin.h>
#endif

#if __has_include(<bit>)
  #include <bit>
#endif

namespace iwbits {

// ------------------------------------------------------------
// Thanks ChatGPT
// Unfortunately it looks like this does not make any difference between
// the 64 bit version.

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)

uint64_t
intersection_popcount_avx512(const uint64_t* a,
                             const uint64_t* b,
                             std::size_t n_words) {
    std::size_t i = 0;
    __m512i vtotal = _mm512_setzero_si512();

    for (; i + 8 <= n_words; i += 8) {
        __m512i va = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(a + i));
        __m512i vb = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(b + i));
        __m512i vand = _mm512_and_si512(va, vb);
        vtotal = _mm512_add_epi64(vtotal, _mm512_popcnt_epi64(vand));
    }
    uint64_t total = _mm512_reduce_add_epi64(vtotal);

    // Tail
    for (; i < n_words; ++i) {
        total += _mm_popcnt_u64(a[i] & b[i]);
    }

    return total;
}

uint64_t
IntersectionPopcount(const uint64_t* a, const uint64_t* b, std::size_t nwords) {
  uint64_t rc = 0;

  for (uint64_t i = 0; i < nwords; ++i) {
    rc += _mm_popcnt_u64(a[i] & b[i]);
  }

  return rc;
}

uint64_t
BitsInCommonAvx(const uint64_t* a, const uint64_t* b, std::size_t nwords) {
  uint64_t bic;
  if (nwords < 8) {
    bic = IntersectionPopcount(a, b, nwords);
  } else {
     bic = intersection_popcount_avx512(a, b, nwords);
  }

  return bic;
}

#else
#endif

}  // namespace iwbits
