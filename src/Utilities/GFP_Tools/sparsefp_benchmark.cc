// Benchmarking to examine performance characteristics of things supporting Sparse_Fingerprints.

#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Utilities/GFP_Tools/sparsefp.h"

namespace {

void
BM_BinarySearch(benchmark::State& state) {
  std::default_random_engine rng;
  int nbits = state.range(0);
  std::uniform_int_distribution<int> count_distribution(1, 255);

  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < nbits; ++i) {
    uint32_t b = rng();
    int count = count_distribution(rng);
    sfc.hit_bit(b, count);
  }
  Sparse_Fingerprint sfp;
  sfp.build_from_sparse_fingerprint_creator(sfc);
  int missing = 0;
  for (auto _ : state) {
    uint32_t b = rng();
    int c = sfp.is_set(b);
    if (c == 0) {
      missing++;
    }
  }
}
BENCHMARK(BM_BinarySearch)
->DenseRange(20, 180, 10);

Sparse_Fingerprint
MakeFingerprint(int nbits) {
  std::default_random_engine rng;
  std::uniform_int_distribution<uint32_t> bits(0, std::numeric_limits<uint32_t>::max());

  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < nbits; ++i) {
    sfc.hit_bit(bits(rng), 1);
  }
  Sparse_Fingerprint result;
  result.build_from_sparse_fingerprint_creator(sfc);
  return result;
}

void
BM_Tanimoto(benchmark::State& state) {
  std::default_random_engine rng;
  const int nbits = state.range(0);
  std::uniform_int_distribution<int> count_distribution(40, nbits + 10);
  const int nbits1 = count_distribution(rng);
  const int nbits2 = count_distribution(rng);

  Sparse_Fingerprint sfp1 = MakeFingerprint(nbits1);
  Sparse_Fingerprint sfp2 = MakeFingerprint(nbits2);
  float maxsim = 0.0f;  // Useless operation to guard against compiler optimization.
  for (auto _ :state) {
    const float t = sfp1.tanimoto_binary(sfp2);
    if (t > maxsim) {
      maxsim = t;
    }
  }
}
BENCHMARK(BM_Tanimoto)->DenseRange(40, 140, 10);

}  // namespace

BENCHMARK_MAIN();
