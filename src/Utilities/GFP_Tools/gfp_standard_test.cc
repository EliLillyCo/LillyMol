#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "gfp_standard.h"

namespace {
// Match the additional state in spread items: array stride must preserve the
// base fingerprint's alignment, rather than just aligning the first element.
struct SpreadItem : GFP_Standard {
  int state[4];
};
static_assert(alignof(GFP_Standard) >= alignof(uint64_t));
static_assert(alignof(SpreadItem) >= alignof(uint64_t));

enum class Builder { kIntegers, kBytes, kFingerprints };

void Build(GFP_Standard& fp, bool second, Builder builder) {
  const int properties[8] = {10, 2, 3, 4, 5, 6, 7, 8};
  fp.build_molecular_properties(properties, 8);
  std::array<int, 2048> iw{};
  std::array<int, 192> mk{}, mk2{};
  iw[63] = mk[63] = 1;
  iw[second ? 2047 : 0] = 1;
  mk[second ? 191 : 0] = 1;
  mk2[second ? 191 : 0] = 1;
  if (builder == Builder::kFingerprints) {
    IWDYFP fiw, fmk, fmk2;
    ASSERT_TRUE(fiw.construct_from_array_of_ints(iw.data(), iw.size()));
    ASSERT_TRUE(fmk.construct_from_array_of_ints(mk.data(), mk.size()));
    ASSERT_TRUE(fmk2.construct_from_array_of_ints(mk2.data(), mk2.size()));
    fp.build_iw(fiw);
    fp.build_mk(fmk);
    fp.build_mk2(fmk2);
  } else {
    if (builder == Builder::kBytes) {
      std::array<unsigned char, 256> bytes{};
      for (int i = 0; i < 2048; ++i) {
        if (iw[i]) {
          bytes[i / 8] |= one_bit_8[i % 8];
        }
      }
      fp.build_iwfp(bytes.data(), 2);
    } else {
      fp.build_iwfp(iw.data(), 2);
    }
    fp.build_mk(mk.data(), mk.size());
    fp.build_mk2(mk2.data(), mk2.size());
  }
}

void CheckPair(const GFP_Standard& lhs, const GFP_Standard& rhs) {
  // Properties contribute 1, IW and MK each 1/3, MK2 contributes zero.
  constexpr float similarity = 5.0f / 12.0f;
  EXPECT_NEAR(lhs.tanimoto(rhs), similarity, 1.0e-7f);
  EXPECT_FLOAT_EQ(lhs.tanimoto(rhs), rhs.tanimoto(lhs));
  EXPECT_FLOAT_EQ(lhs.tanimoto(lhs), 1.0f);
  EXPECT_NEAR(lhs.tanimoto_distance(rhs), 1.0f - similarity, 1.0e-7f);
  const auto distance = lhs.tanimoto_distance_if_less(rhs, 0.7f);
  ASSERT_TRUE(distance.has_value());
  EXPECT_NEAR(*distance, 1.0f - similarity, 1.0e-7f);
  EXPECT_FALSE(lhs.tanimoto_distance_if_less(rhs, 0.5f).has_value());

}

TEST(GfpStandard, StackArray) {
  std::array<GFP_Standard, 4> pool;
  for (int i = 0; i < 4; ++i) {
    Build(pool[i], i % 2, Builder::kIntegers);
  }
  CheckPair(pool[0], pool[1]);
  CheckPair(pool[2], pool[3]);
}

TEST(GfpStandard, HeapArray) {
  auto pool = std::make_unique<GFP_Standard[]>(4);
  for (int i = 0; i < 4; ++i) {
    Build(pool[i], i % 2, Builder::kFingerprints);
  }
  CheckPair(pool[0], pool[1]);
  CheckPair(pool[2], pool[3]);
}

TEST(GfpStandard, DerivedVector) {
  std::vector<SpreadItem> pool(4);
  for (int i = 0; i < 4; ++i) {
    Build(pool[i], i % 2, Builder::kBytes);
  }
  CheckPair(pool[0], pool[1]);
  CheckPair(pool[2], pool[3]);
}

TEST(GfpStandard, BuildersAndCopiesAgree) {
  GFP_Standard integers, bytes, fingerprints;
  Build(integers, false, Builder::kIntegers);
  Build(bytes, false, Builder::kBytes);
  Build(fingerprints, false, Builder::kFingerprints);
  EXPECT_FLOAT_EQ(integers.tanimoto(bytes), 1.0f);
  EXPECT_FLOAT_EQ(integers.tanimoto(fingerprints), 1.0f);
  GFP_Standard copy = integers;
  EXPECT_FLOAT_EQ(copy.tanimoto(integers), 1.0f);
  Build(copy, true, Builder::kFingerprints);
  CheckPair(copy, integers);
  copy = integers;
  EXPECT_FLOAT_EQ(copy.tanimoto(integers), 1.0f);
}
}  // namespace
