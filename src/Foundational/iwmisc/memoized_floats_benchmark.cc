
#include <random>
#include <benchmark/benchmark.h>
#include "Foundational/iwmisc/memoized_floats.h"

namespace memoized_floats {
namespace {

static void BM_WithoutMemoization(benchmark::State& state) {
  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(state.range(0), state.range(1));

  for (auto _ : state) {
    IWString buffer;
    float value = distribution(generator);
    buffer << value;
  }
}
BENCHMARK(BM_WithoutMemoization)
->Args({0,1})
->Args({-10, 10});

static void BM_WithMemoization(benchmark::State& state) {
  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(state.range(0), state.range(1));
  MemoizedFloats mf;
  mf.Build(10, 4);

  for (auto _ : state) {
    IWString buffer;
    float value = distribution(generator);
    mf.Representation(value);
  }
}
BENCHMARK(BM_WithMemoization)
->Args({0,1})
->Args({-10, 10});


}  // namespace
}  // namespace memoized_floats
