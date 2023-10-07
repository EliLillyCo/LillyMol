// Tester for various Sparse_Fingerprint things.

#include <array>
#include <cstdint>
#include <limits>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Utilities/GFP_Tools/sparsefp.h"

namespace {

class TestSparseFingerprint : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    Sparse_Fingerprint_Creator _sfc;
};

Sparse_Fingerprint
ASparseFingerprint(const int nbits,
                   uint32_t bit,
                   const uint32_t bit_delta,
                   int count,
                   const int count_delta) {
  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < nbits; ++i) {
    sfc.hit_bit(bit, count);
    bit += bit_delta;
    count += count_delta;
    if (count > std::numeric_limits<uint8_t>::max()) {
      count = 10;   // An arbitrary choice.
    }
  }

  Sparse_Fingerprint result;
  result.build_from_sparse_fingerprint_creator(sfc);
  return result;
}

Sparse_Fingerprint
ASparseFingerprint(const std::vector<uint32_t>& bits,
                   const std::vector<uint32_t>& count) {
  Sparse_Fingerprint_Creator sfc;
  for (uint32_t i = 0; i < bits.size(); ++i) {
    sfc.hit_bit(bits[i], count[i]);
  }

  Sparse_Fingerprint result;
  result.build_from_sparse_fingerprint_creator(sfc);
  return result;
}

TEST(TestSparseFingerprint, TestEmpty) {
  Sparse_Fingerprint_Creator sfc;
  Sparse_Fingerprint sfp1;
  ASSERT_EQ(sfp1.build_from_sparse_fingerprint_creator(sfc), 1);
  EXPECT_EQ(sfp1.nbits(), 0);
  EXPECT_EQ(sfp1.nset(), 0);

  Sparse_Fingerprint sfp2;
  ASSERT_EQ(sfp2.build_from_sparse_fingerprint_creator(sfc), 1);
  EXPECT_EQ(sfp2.nbits(), 0);
  EXPECT_EQ(sfp2.nset(), 0);

  EXPECT_EQ(sfp1.bits_in_common(sfp2), 0);
}

TEST(TestSparseFingerprint, TestBuild) {
  Sparse_Fingerprint_Creator sfc;
  constexpr int nbits = 1081;
  std::array<uint32_t, nbits> bits;
  std::array<uint32_t, nbits> counts;
  uint32_t b = 2;
  int count = 2;
  uint32_t nset = 0;
  int trunc50 = 0;
  for (int i = 0; i < nbits; ++i) {
    b += 133;
    bits[i] = b;
    count += 5;
    if (count > std::numeric_limits<uint8_t>::max()) {
      count = 4;
    }
    sfc.hit_bit(b, count);
    counts[i] = count;
    nset += count;
    if (count > 50) {
      trunc50 += 50;
    } else {
      trunc50 += count;
    }
  }

  Sparse_Fingerprint sfp;
  ASSERT_TRUE(sfp.build_from_sparse_fingerprint_creator(sfc));

  for (int i = 0; i < nbits; ++i) {
    EXPECT_TRUE(sfp.is_set(bits[i]));
    EXPECT_EQ(sfp.count_for_bit(bits[i]), counts[i]);
  }

  EXPECT_EQ(sfp.nbits(), nbits);
  EXPECT_EQ(sfp.nset(), nset);

  sfp.truncate_counts_at(50);
  EXPECT_EQ(sfp.nset(), trunc50);
}

TEST(TestSparseFingerprint, TestOperatorEquals) {
  Sparse_Fingerprint_Creator sfc;
  constexpr int nbits = 10;
  std::array<uint32_t, nbits> bits;
  std::array<uint32_t, nbits> counts;
  int b = 13;
  int count = 2;
  int nset = 0;
  for (int i = 0; i < nbits; ++i) {
    bits[i] = b;
    counts[i] = count;
    sfc.hit_bit(b, count);
    nset += count;

    b += 157;
    count += 17;
    if (count > std::numeric_limits<uint8_t>::max()) {
      count = 3;
    }
  }

  Sparse_Fingerprint sfp1;
  ASSERT_TRUE(sfp1.build_from_sparse_fingerprint_creator(sfc));

  const Sparse_Fingerprint sfp2 = sfp1;

  EXPECT_EQ(sfp1.nbits(), sfp2.nbits());
  EXPECT_EQ(sfp1.nset(), sfp2.nset());
  EXPECT_EQ(sfp1.bits_in_common(sfp2), sfp2.bits_in_common(sfp1));
  EXPECT_EQ(sfp1.bits_in_common(sfp2), nset);
}

TEST(TestSparseFingerprint, TestConstructor) {
  constexpr int nbits = 78;
  const Sparse_Fingerprint sfp1 = ASparseFingerprint(nbits, 200000, 110000, 6, 7);
  const Sparse_Fingerprint sfp2(sfp1);
  EXPECT_EQ(sfp1.nbits(), nbits);
  EXPECT_EQ(sfp2.nbits(), nbits);
  EXPECT_EQ(sfp1.nset(), sfp2.nset());

  int i = 0;
  unsigned int bit;
  int count;
  while (sfp1.next_bit_set(i, bit, count)) {
    EXPECT_EQ(sfp2.count_for_bit(bit), count);
  }
}

TEST(TestSparseFingerprint, TestRemoveBit) {
  Sparse_Fingerprint_Creator sfc;
  constexpr int nbits = 10;
  std::array<uint32_t, nbits> bits;
  std::array<uint32_t, nbits> counts;
  int b = 13;
  int count = 2;
  int nset = 0;
  for (int i = 0; i < nbits; ++i) {
    bits[i] = b;
    counts[i] = count;
    sfc.hit_bit(b, count);
    nset += count;

    b += 157;
    count += 17;
    if (count > std::numeric_limits<uint8_t>::max()) {
      count = 3;
    }
  }

  Sparse_Fingerprint sfp;
  ASSERT_TRUE(sfp.build_from_sparse_fingerprint_creator(sfc));

  ASSERT_EQ(sfp.nbits(), nbits);
  ASSERT_EQ(sfp.nset(), nset);
  for (int i = 0; i < nbits; ++i) {
    EXPECT_EQ(sfp.count_for_bit(bits[i]), counts[i]);
  }

  EXPECT_EQ(sfp.remove_bit(bits[3]), 1);

  EXPECT_EQ(sfp.nbits(), nbits - 1);
  EXPECT_EQ(sfp.nset(), nset - counts[3]);

  // Remove the last bit.
  EXPECT_EQ(sfp.remove_bit(bits[nbits - 1]), 1);
  EXPECT_EQ(sfp.nbits(), nbits - 2);
  EXPECT_EQ(sfp.nset(), nset - counts[3] - counts[nbits - 1]);
}

TEST(TestSparseFingerprint, TestSubset) {
  Sparse_Fingerprint_Creator sfc;
  constexpr int nbits = 7;
  std::array<uint32_t, nbits> bits {32768, 0, 19282, 11, 1, 
                             std::numeric_limits<uint32_t>::max(),
                             104};
  std::array<uint32_t, nbits> count {1, 5, 3, 1, 4, 255, 1};

  for (int i = 0; i < nbits; ++i) {
    sfc.hit_bit(bits[i], count[i]);
  }
  
  Sparse_Fingerprint sfp;
  ASSERT_TRUE(sfp.build_from_sparse_fingerprint_creator(sfc));
  EXPECT_EQ(sfp.nbits(), nbits);

  for (int i = 0; i < nbits; ++i) {
    int c = sfp.count_for_bit(bits[i]);
    EXPECT_EQ(c, count[i]);
  }

  std::unordered_set<uint32_t> to_keep {bits[2], bits[4], bits[5], bits[6]};
  EXPECT_EQ(sfp.ReduceToSubset(to_keep), 3);
  EXPECT_EQ(sfp.nbits(), 4);
  EXPECT_EQ(sfp.nset(), count[2] + count[4] + count[5] + count[6]);
  EXPECT_EQ(sfp.count_for_bit(bits[2]), count[2]);
  EXPECT_EQ(sfp.count_for_bit(bits[4]), count[4]);
  EXPECT_EQ(sfp.count_for_bit(bits[5]), count[5]);
  EXPECT_EQ(sfp.count_for_bit(bits[6]), count[6]);

  EXPECT_EQ(sfp.count_for_bit(bits[0]), 0);
  EXPECT_EQ(sfp.count_for_bit(bits[1]), 0);
  EXPECT_EQ(sfp.count_for_bit(bits[3]), 0);
}

TEST(TestSparseFingerprint, TestIsSet) {
  Sparse_Fingerprint_Creator sfc;
  constexpr int nbits = 97;
  std::array<uint32_t, nbits> bits;
  std::array<uint32_t, nbits> counts;
  uint32_t b = 8;
  int count = 23;
  for (int i = 0; i < nbits; ++i) {
    bits[i] = b;
    counts[i] = count;
    sfc.hit_bit(b, count);
    b += 201;
    count += 1;

  }

  Sparse_Fingerprint sfp;
  ASSERT_TRUE(sfp.build_from_sparse_fingerprint_creator(sfc));

  for (int i = 0; i < nbits; ++i) {
    EXPECT_TRUE(sfp.is_set(bits[i]));
    EXPECT_FALSE(sfp.is_set(bits[i] - 1));
    EXPECT_FALSE(sfp.is_set(bits[i] + 1));
  }
}

TEST(TestSparseFingerprint, TestTanimoto1) {
  const std::vector<uint32_t> bits1  {0, 1, 2, 3};
  const std::vector<uint32_t> count1 {1, 1, 1, 1};
  const std::vector<uint32_t> bits2        {2, 3, 4, 5};
  const std::vector<uint32_t> count2       {1, 1, 4, 5};

  const Sparse_Fingerprint sfp1 = ASparseFingerprint(bits1, count1);
  const Sparse_Fingerprint sfp2 = ASparseFingerprint(bits2, count2);
  EXPECT_EQ(sfp1.bits_in_common(sfp2), 2);
  EXPECT_EQ(sfp2.bits_in_common(sfp1), 2);
  // bits_in_common 2
  // Just1 = 2
  // Just2 = 9
  float expected = 2.0 / static_cast<float>(2 + 9 + 2);
  EXPECT_FLOAT_EQ(sfp1.tanimoto(sfp2), expected);
  EXPECT_FLOAT_EQ(sfp2.tanimoto(sfp1), expected);
  EXPECT_FLOAT_EQ(sfp1.distance(sfp2), 1.0f - expected);

  EXPECT_FLOAT_EQ(sfp1.dot_product(sfp2), 2);
  EXPECT_FLOAT_EQ(sfp2.dot_product(sfp1), 2);
}

TEST(TestSparseFingerprint, TestTanimoto2) {
  const std::vector<uint32_t> bits1  {0, 2, 3,      10, 13, 15, 17};
  const std::vector<uint32_t> count1 {1, 1, 1,       1,  1,  2,  1};
  const std::vector<uint32_t> bits2  {   2,   4, 6, 10,     15};
  const std::vector<uint32_t> count2 {   1,   1, 4,  2,      2};

  // bits_in_common 4
  // just1 = 4
  // just2 = 6
  const Sparse_Fingerprint sfp1 = ASparseFingerprint(bits1, count1);
  const Sparse_Fingerprint sfp2 = ASparseFingerprint(bits2, count2);
  EXPECT_EQ(sfp1.bits_in_common(sfp2), 4);
  EXPECT_EQ(sfp2.bits_in_common(sfp1), 4);
  float expected = 4.0 / static_cast<float>(4 + 6 + 4);
  EXPECT_FLOAT_EQ(sfp1.tanimoto(sfp2), expected);
  EXPECT_FLOAT_EQ(sfp2.tanimoto(sfp1), expected);
  EXPECT_FLOAT_EQ(sfp1.distance(sfp2), 1.0f - expected);

  EXPECT_FLOAT_EQ(sfp1.dot_product(sfp2), 7);
  EXPECT_FLOAT_EQ(sfp2.dot_product(sfp1), 7);
}

// The cosine coefficient depends on _sum_squared being set
// properly, so test the operations that update it.
TEST(TestSparseFingerprint, TestSumSquared) {
  const std::vector<uint32_t> bits1  {0, 1, 2, 3};
  const std::vector<uint32_t> count1 {1, 1, 1, 1};
  const std::vector<uint32_t> bits2        {2, 3, 4, 5};
  const std::vector<uint32_t> count2       {1, 1, 4, 5};

  Sparse_Fingerprint sfp1 = ASparseFingerprint(bits1, count1);
  Sparse_Fingerprint sfp2 = ASparseFingerprint(bits2, count2);
  EXPECT_FLOAT_EQ(sfp1.cosine_coefficient(sfp1), 1.0);

  sfp1.set_count(bits1[1], count1[1] + 1);
  EXPECT_FLOAT_EQ(sfp1.cosine_coefficient(sfp1), 1.0);

  EXPECT_FLOAT_EQ(sfp2.cosine_coefficient(sfp2), 1.0);
  sfp2.remove_bit(bits2[3]);
  EXPECT_FLOAT_EQ(sfp2.cosine_coefficient(sfp2), 1.0);
  sfp2.truncate_counts_at(2);
  EXPECT_FLOAT_EQ(sfp2.cosine_coefficient(sfp2), 1.0);
}

TEST(TestSparseFingerprint, TestCosineCoefficient) {
  const std::vector<uint32_t> bits1  {0, 1, 2, 3};
  const std::vector<uint32_t> count1 {1, 1, 1, 2};
  const std::vector<uint32_t> bits2        {2, 3, 4, 5};
  const std::vector<uint32_t> count2       {1, 1, 4, 5};

  Sparse_Fingerprint sfp1 = ASparseFingerprint(bits1, count1);
  Sparse_Fingerprint sfp2 = ASparseFingerprint(bits2, count2);
  EXPECT_FLOAT_EQ(sfp1.cosine_coefficient(sfp1), 1.0);
  EXPECT_FLOAT_EQ(sfp2.cosine_coefficient(sfp2), 1.0);
  EXPECT_FLOAT_EQ(sfp1.cosine_coefficient(sfp2),
                  sfp2.cosine_coefficient(sfp1));
  // Sum of products 3
  // Sum squares 7 and 43
  const float expected_result = 3.0 / sqrt(7 * 43);
  EXPECT_FLOAT_EQ(sfp1.cosine_coefficient(sfp2), expected_result);
  EXPECT_FLOAT_EQ(sfp2.cosine_coefficient(sfp1), expected_result);
}

TEST(TestSparseFingerprint, TestTanimotoBinary) {
  const std::vector<uint32_t> bits1  {0, 1, 2, 3, 4};
  const std::vector<uint32_t> count1 {1, 1, 1, 2, 3};
  const std::vector<uint32_t> bits2        {2, 3, 4, 5, 6};
  const std::vector<uint32_t> count2       {1, 1, 4, 5, 1};

  Sparse_Fingerprint sfp1 = ASparseFingerprint(bits1, count1);
  Sparse_Fingerprint sfp2 = ASparseFingerprint(bits2, count2);
  EXPECT_EQ(sfp1.bits_in_common(sfp1), sfp1.nset());
  EXPECT_EQ(sfp2.bits_in_common(sfp2), sfp2.nset());
  EXPECT_EQ(sfp1.bits_in_common(sfp2), 5);
  EXPECT_FLOAT_EQ(sfp1.tanimoto(sfp1), 1.0f);
  EXPECT_FLOAT_EQ(sfp2.tanimoto(sfp2), 1.0f);
  EXPECT_FLOAT_EQ(sfp2.tanimoto(sfp2), sfp1.tanimoto(sfp1));
  // bits_in_common 5
  // just1 3 just2 7
  float expected_result = 5.0 / (3 + 7 + 5);
  EXPECT_FLOAT_EQ(sfp1.tanimoto(sfp2), expected_result);
  EXPECT_FLOAT_EQ(sfp2.tanimoto(sfp1), expected_result);
  expected_result = 3.0 / (2 + 2 + 3);
  EXPECT_FLOAT_EQ(sfp1.tanimoto_binary(sfp1), 1.0f);
  EXPECT_FLOAT_EQ(sfp2.tanimoto_binary(sfp2), 1.0f);
  EXPECT_FLOAT_EQ(sfp1.tanimoto_binary(sfp2), sfp2.tanimoto_binary(sfp1));
  EXPECT_FLOAT_EQ(sfp1.tanimoto_binary(sfp2), expected_result);
}

TEST(TestSparseFingerprint, TestManhattan) {
  const std::vector<uint32_t> bits1  {0, 1, 2, 3, 4};
  const std::vector<uint32_t> count1 {1, 1, 1, 2, 3};
  const std::vector<uint32_t> bits2        {2, 3, 4, 5, 6};
  const std::vector<uint32_t> count2       {1, 1, 4, 5, 1};

  Sparse_Fingerprint sfp1 = ASparseFingerprint(bits1, count1);
  Sparse_Fingerprint sfp2 = ASparseFingerprint(bits2, count2);
  EXPECT_EQ(sfp1.bits_in_common(sfp1), sfp1.nset());
  EXPECT_EQ(sfp2.bits_in_common(sfp2), sfp2.nset());
  EXPECT_EQ(sfp1.bits_in_common(sfp2), 5);
  EXPECT_EQ(sfp1.manhattan_distance(sfp1), 1.0f);
  EXPECT_EQ(sfp2.manhattan_distance(sfp2), 1.0f);
  EXPECT_EQ(sfp1.manhattan_distance(sfp2), 1.0f / (1 + 10));
  EXPECT_EQ(sfp2.manhattan_distance(sfp1), 1.0f / (1 + 10));
}

TEST(TestSparseFingerprint, DaylightRoundTrip) {
  const std::vector<uint32_t> bits  {1049, 801, 10000000, 4, 32769, 7777777, 5011};
  const std::vector<uint32_t> count {1,    5,   255,      1, 251,   1,       189};
  Sparse_Fingerprint sfp1 = ASparseFingerprint(bits, count);
  IWString ascii;
  sfp1.append_daylight_ascii_form_with_counts_encoded(ascii);

  Sparse_Fingerprint sfp2;
  ASSERT_TRUE(sfp2.construct_from_daylight_ascii_representation(ascii));
  EXPECT_FLOAT_EQ(sfp1.tanimoto(sfp2), 1.0f);
}

TEST(TestSparseFingerprint, TdtRecordRoundTrip) {
  const std::vector<uint32_t> bits  {1049, 801, 10000000, 4, 32769, 7777777, 5011};
  const std::vector<uint32_t> count {1,    5,   255,      1, 251,   1,       189};
  Sparse_Fingerprint sfp1 = ASparseFingerprint(bits, count);
  IWString ascii;
  ascii << "NCFOOBAR<";
  sfp1.append_daylight_ascii_form_with_counts_encoded(ascii);
  ascii << ">";

  Sparse_Fingerprint sfp2;
  ASSERT_TRUE(sfp2.construct_from_tdt_record(ascii));
  EXPECT_FLOAT_EQ(sfp1.tanimoto(sfp2), 1.0f);
}

TEST(TestSparseFingerprint, FromSvmLite) {
  set_sparse_ascii_separator(':');
  IWString ascii("33:1 901:504 14002:1000000 32769:2000000 45000:1");

  Sparse_Fingerprint sfp;
  ASSERT_TRUE(sfp.construct_from_sparse_ascii_representation(ascii));
  EXPECT_EQ(sfp.count_for_bit(33), 1);
  EXPECT_EQ(sfp.count_for_bit(901), 504);
  EXPECT_EQ(sfp.count_for_bit(14002), 1000000);
  EXPECT_EQ(sfp.count_for_bit(32769), 2000000);
  EXPECT_EQ(sfp.count_for_bit(45000), 1);
}

TEST(TestSparseFingerprint, TestUncounted) {
  const std::vector<uint32_t> bits{98000, 12, 13, 65000, 4, 88};
  Sparse_Fingerprint_Creator sfc;
  for (auto b : bits) {
    sfc.hit_bit(b);
  }

  const IWString daylight = sfc.BitsWithoutCounts();

  Sparse_Fingerprint sfp;

  ASSERT_TRUE(sfp.construct_from_daylight_ascii_representation_uncounted(daylight));

  EXPECT_EQ(sfp.nbits(), bits.size());
  for (auto b : bits) {
    EXPECT_TRUE(sfp.is_set(b));
  }

  EXPECT_FLOAT_EQ(sfp.tanimoto_binary(sfp), 1.0f);
}

} // namespace
