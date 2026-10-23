// Benchmark molecular-property similarity implementations.
//
// Real input:
//   bazel-bin/Molecule_Tools/temperature -J MPR file.smi > /tmp/file.mpr.gfp
//   MPR_BENCHMARK_INPUT=/tmp/file.mpr.gfp bazel-bin/Utilities/GFP_Tools/mpr_similarity_benchmark
//
// Without MPR_BENCHMARK_INPUT, synthetic data is generated. Optional controls:
//   MPR_BENCHMARK_MOLECULES, MPR_BENCHMARK_PROPERTIES, MPR_BENCHMARK_PAIRS.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#include "benchmark/benchmark.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Utilities/GFP_Tools/dyfp.h"

namespace {

struct BenchmarkData {
  int nproperties = 8;
  std::vector<int> properties;
  std::vector<uint8_t> packed_properties;
  std::vector<std::pair<uint32_t, uint32_t>> pairs;
  std::vector<float> ratio_table;

  uint32_t nmolecules() const {
    if (nproperties == 0) {
      return 0;
    }
    return properties.size() / nproperties;
  }

  const int* row(uint32_t ndx) const {
    return properties.data() + static_cast<size_t>(ndx) * nproperties;
  }

  const uint8_t* packed_row(uint32_t ndx) const {
    return packed_properties.data() + static_cast<size_t>(ndx) * nproperties;
  }
};

BenchmarkData g_data;

float
Ratio(int lhs, int rhs) {
  if (lhs == rhs) {
    return 1.0f;
  }

  if (lhs == 0 || rhs == 0) {
    return 0.5f;
  }

  if (lhs < rhs) {
    return static_cast<float>(lhs) / static_cast<float>(rhs);
  }

  return static_cast<float>(rhs) / static_cast<float>(lhs);
}

void
InitialiseRatioTable(BenchmarkData& data) {
  data.ratio_table.resize(256 * 256);
  const float scale = 1.0f / static_cast<float>(data.nproperties);

  for (int i = 0; i < 256; ++i) {
    for (int j = 0; j < 256; ++j) {
      data.ratio_table[i * 256 + j] = Ratio(i, j) * scale;
    }
  }
}

float
TableSimilarity(const int* lhs, const int* rhs, int nproperties,
                const std::vector<float>& table) {
  float result = 0.0f;
  for (int i = 0; i < nproperties; ++i) {
    result += table[lhs[i] * 256 + rhs[i]];
  }

  return result;
}

float
PackedTableSimilarity(const uint8_t* lhs, const uint8_t* rhs, int nproperties,
                      const std::vector<float>& table) {
  float result = 0.0f;
  for (int i = 0; i < nproperties; ++i) {
    result += table[static_cast<int>(lhs[i]) * 256 + static_cast<int>(rhs[i])];
  }

  return result;
}

float
ScalarSimilarity(const int* lhs, const int* rhs, int nproperties) {
  float result = 0.0f;
  for (int i = 0; i < nproperties; ++i) {
    result += Ratio(lhs[i], rhs[i]);
  }

  return result / static_cast<float>(nproperties);
}

#if defined(__AVX2__)
float
Avx2Similarity(const int* lhs, const int* rhs, int nproperties) {
  const __m256i zero_i = _mm256_setzero_si256();
  const __m256 one = _mm256_set1_ps(1.0f);
  const __m256 half = _mm256_set1_ps(0.5f);
  __m256 sum = _mm256_setzero_ps();

  for (int i = 0; i < nproperties; i += 8) {
    const __m256i l = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(lhs + i));
    const __m256i r = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(rhs + i));

    const __m256i eq = _mm256_cmpeq_epi32(l, r);
    const __m256i l_zero = _mm256_cmpeq_epi32(l, zero_i);
    const __m256i r_zero = _mm256_cmpeq_epi32(r, zero_i);
    const __m256i either_zero = _mm256_or_si256(l_zero, r_zero);

    const __m256i mn = _mm256_min_epi32(l, r);
    const __m256i mx = _mm256_max_epi32(l, r);
    const __m256 ratio = _mm256_div_ps(_mm256_cvtepi32_ps(mn), _mm256_cvtepi32_ps(mx));

    __m256 contribution = _mm256_blendv_ps(ratio, half, _mm256_castsi256_ps(either_zero));
    contribution = _mm256_blendv_ps(contribution, one, _mm256_castsi256_ps(eq));
    sum = _mm256_add_ps(sum, contribution);
  }

  alignas(32) float tmp[8];
  _mm256_store_ps(tmp, sum);
  float result = tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];

  return result / static_cast<float>(nproperties);
}
#endif

bool
ReadMprRecord(const const_IWSubstring& buffer, BenchmarkData& data) {
  IWDYFP fp;
  if (!fp.construct_from_tdt_record(buffer)) {
    std::cerr << "Cannot parse MPR record '" << buffer << "'\n";
    return false;
  }

  const int nproperties = fp.nbits() / 8;
  if (nproperties == 0 || fp.nbits() % 8 != 0) {
    std::cerr << "Invalid MPR bit count " << fp.nbits() << '\n';
    return false;
  }

  if (data.properties.empty()) {
    data.nproperties = nproperties;
    if (data.nproperties % 8 != 0) {
      std::cerr << "Number of properties must be a multiple of 8, got "
                << data.nproperties << '\n';
      return false;
    }
  } else if (nproperties != data.nproperties) {
    std::cerr << "Inconsistent property count, got " << nproperties
              << " expected " << data.nproperties << '\n';
    return false;
  }

  const unsigned char* bits = reinterpret_cast<const unsigned char*>(fp.bits());
  for (int i = 0; i < data.nproperties; ++i) {
    data.properties.push_back(static_cast<int>(bits[i]));
    data.packed_properties.push_back(bits[i]);
  }

  return true;
}

bool
ReadMprFile(const char* fname, BenchmarkData& data) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    std::cerr << "Cannot open '" << fname << "'\n";
    return false;
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!buffer.starts_with("MPR<")) {
      continue;
    }

    if (!ReadMprRecord(buffer, data)) {
      return false;
    }
  }

  if (data.properties.empty()) {
    std::cerr << "No MPR records in '" << fname << "'\n";
    return false;
  }

  return true;
}

int
EnvironmentValue(const char* env, int default_value) {
  const char* value = std::getenv(env);
  if (value == nullptr || *value == '\0') {
    return default_value;
  }

  return std::max(1, std::atoi(value));
}

void
GenerateSyntheticData(BenchmarkData& data) {
  data.nproperties = EnvironmentValue("MPR_BENCHMARK_PROPERTIES", 8);
  if (data.nproperties % 8 != 0) {
    std::cerr << "MPR_BENCHMARK_PROPERTIES must be a multiple of 8\n";
    std::exit(1);
  }

  const int nmolecules = EnvironmentValue("MPR_BENCHMARK_MOLECULES", 10000);
  data.properties.reserve(static_cast<size_t>(nmolecules) * data.nproperties);
  data.packed_properties.reserve(static_cast<size_t>(nmolecules) * data.nproperties);

  std::mt19937 rng(8675309);
  std::poisson_distribution<int> small_count(3);
  std::uniform_int_distribution<int> occasional_large(0, 40);
  for (int i = 0; i < nmolecules; ++i) {
    for (int j = 0; j < data.nproperties; ++j) {
      int value = small_count(rng);
      if ((rng() & 0x1f) == 0) {
        value += occasional_large(rng);
      }
      value = std::min(value, 255);
      data.properties.push_back(value);
      data.packed_properties.push_back(static_cast<uint8_t>(value));
    }
  }
}

void
GeneratePairs(BenchmarkData& data) {
  const uint32_t nmolecules = data.nmolecules();
  const int npairs = EnvironmentValue("MPR_BENCHMARK_PAIRS", 1000000);
  data.pairs.reserve(npairs);

  std::mt19937 rng(424242);
  std::uniform_int_distribution<uint32_t> molecule(0, nmolecules - 1);
  for (int i = 0; i < npairs; ++i) {
    data.pairs.emplace_back(molecule(rng), molecule(rng));
  }
}


bool
ApproximatelyEqual(float lhs, float rhs) {
  const float delta = lhs > rhs ? lhs - rhs : rhs - lhs;
  return delta < 1.0e-6f;
}

void
CheckImplementationsAgree() {
  const size_t ncheck = std::min<size_t>(g_data.pairs.size(), 10000);
  for (size_t i = 0; i < ncheck; ++i) {
    const auto& [lhs_ndx, rhs_ndx] = g_data.pairs[i];
    const int* lhs = g_data.row(lhs_ndx);
    const int* rhs = g_data.row(rhs_ndx);
    const float table = TableSimilarity(lhs, rhs, g_data.nproperties, g_data.ratio_table);
    const float packed = PackedTableSimilarity(g_data.packed_row(lhs_ndx), g_data.packed_row(rhs_ndx),
                                               g_data.nproperties, g_data.ratio_table);
    if (!ApproximatelyEqual(table, packed)) {
      std::cerr << "Table/packed mismatch " << table << " vs " << packed << '\n';
      std::exit(1);
    }
    const float scalar = ScalarSimilarity(lhs, rhs, g_data.nproperties);
    if (!ApproximatelyEqual(table, scalar)) {
      std::cerr << "Table/scalar mismatch " << table << " vs " << scalar << '\n';
      std::exit(1);
    }
#if defined(__AVX2__)
    const float avx2 = Avx2Similarity(lhs, rhs, g_data.nproperties);
    if (!ApproximatelyEqual(table, avx2)) {
      std::cerr << "Table/AVX2 mismatch " << table << " vs " << avx2 << '\n';
      std::exit(1);
    }
#endif
  }
}

void
InitialiseData() {
  const char* fname = std::getenv("MPR_BENCHMARK_INPUT");
  if (fname != nullptr && *fname != '\0') {
    if (!ReadMprFile(fname, g_data)) {
      std::exit(1);
    }
  } else {
    GenerateSyntheticData(g_data);
  }

  InitialiseRatioTable(g_data);
  GeneratePairs(g_data);
  CheckImplementationsAgree();

  std::cerr << "MPR benchmark loaded " << g_data.nmolecules() << " molecules, "
            << g_data.nproperties << " properties, " << g_data.pairs.size()
            << " comparison pairs\n";
}

void
BM_TableSimilarity(benchmark::State& state) {
  float total = 0.0f;
  for (auto _ : state) {
    for (const auto& [lhs, rhs] : g_data.pairs) {
      total += TableSimilarity(g_data.row(lhs), g_data.row(rhs), g_data.nproperties,
                               g_data.ratio_table);
    }
  }
  benchmark::DoNotOptimize(total);
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(g_data.pairs.size()));
}
BENCHMARK(BM_TableSimilarity);


void
BM_PackedTableSimilarity(benchmark::State& state) {
  float total = 0.0f;
  for (auto _ : state) {
    for (const auto& [lhs, rhs] : g_data.pairs) {
      total += PackedTableSimilarity(g_data.packed_row(lhs), g_data.packed_row(rhs),
                                     g_data.nproperties, g_data.ratio_table);
    }
  }
  benchmark::DoNotOptimize(total);
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(g_data.pairs.size()));
}
BENCHMARK(BM_PackedTableSimilarity);

void
BM_ScalarSimilarity(benchmark::State& state) {
  float total = 0.0f;
  for (auto _ : state) {
    for (const auto& [lhs, rhs] : g_data.pairs) {
      total += ScalarSimilarity(g_data.row(lhs), g_data.row(rhs), g_data.nproperties);
    }
  }
  benchmark::DoNotOptimize(total);
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(g_data.pairs.size()));
}
BENCHMARK(BM_ScalarSimilarity);

#if defined(__AVX2__)
void
BM_Avx2Similarity(benchmark::State& state) {
  float total = 0.0f;
  for (auto _ : state) {
    for (const auto& [lhs, rhs] : g_data.pairs) {
      total += Avx2Similarity(g_data.row(lhs), g_data.row(rhs), g_data.nproperties);
    }
  }
  benchmark::DoNotOptimize(total);
  state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(g_data.pairs.size()));
}
BENCHMARK(BM_Avx2Similarity);
#endif

}  // namespace

int
main(int argc, char** argv) {
  benchmark::Initialize(&argc, argv);
  InitialiseData();
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
