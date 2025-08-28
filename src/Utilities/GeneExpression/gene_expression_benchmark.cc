// Benchmarking for the GeneProfile class

#include <algorithm>
#include <numeric>
#include <random>
#include <string>

#include "benchmark/benchmark.h"

#include "gene_expression.h"

namespace gene_expression {

namespace {
// FIll the `ngenes` rank values in `proto` with random values
void
FillRank(int ngenes, int* tmp, std::mt19937_64& rng, Profile& proto) {
  std::iota(tmp, tmp + ngenes, 1);
  std::shuffle(tmp, tmp + ngenes, rng);

  std::bernoulli_distribution bernouilli(0.50);

  for (int i = 0; i < ngenes; ++i) {
    if (bernouilli(rng)) {
      tmp[i] = - tmp[i];
    }
  }

  proto.mutable_rank()->Reserve(ngenes);
  for (int i = 0; i < ngenes; ++i) {
    proto.add_rank(tmp[i]);
  }
}

void
BM_Association(benchmark::State& state) {
  const int ngenes = state.range(0);

  GeneProfile p1;
  GeneProfile p2;

  Profile proto1;
  Profile proto2;

  proto1.set_name("id1");
  proto2.set_name("id2");

  std::random_device rd;
  std::mt19937_64 rng(rd());

  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(ngenes);

  FillRank(ngenes, tmp.get(), rng, proto1);
  FillRank(ngenes, tmp.get(), rng, proto2);

  p1.Build(proto1);
  p2.Build(proto2);

  double total = 0.0;
  for (auto _ : state) {
    total += p1.Association(p2);
  }
}
BENCHMARK(BM_Association)->RangeMultiplier(2)->Range(64, 512);

}  // namespace

}  // namespace gene_expression

BENCHMARK_MAIN();
