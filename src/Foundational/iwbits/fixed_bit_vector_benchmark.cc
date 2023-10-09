// Benchmarking FixedBitVector BitsInCommon function.

#include <random>

#include <benchmark/benchmark.h>

#include "fixed_bit_vector.h"

namespace {
void
BM_BitsInCommon(benchmark::State& state) {
  const int nbits = state.range(0);

  fixed_bit_vector::FixedBitVector bitvector(nbits);

  std::default_random_engine rng;
  std::uniform_int_distribution<int> nset_distribution(1, nbits);

  const int on_bits = nset_distribution(rng);
  // No attempt to ensure different bits selected.
  for (int i = 0; i < on_bits; ++i) {
    int b = nset_distribution(rng);
    bitvector.set_bit(b);
  }

  int total = 0;
  for (auto _ : state) {
    int bits_in_common = bitvector.BitsInCommon(bitvector);
    total += bits_in_common;
  }

  // Useless thing to try and ensure the compiler does not
  // optimize things away.
  if (total == -3) {
    exit(1);
  }
}
BENCHMARK(BM_BitsInCommon)->RangeMultiplier(2)->Range(64, 2048);

}  // namespace

BENCHMARK_MAIN();
