#include "Utilities/GFP_Tools/fixed_size_counted_fingerprint.h"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

Fixed_Size_Counted_Fingerprint_uchar
Fingerprint(std::initializer_list<unsigned char> values) {
  Fixed_Size_Counted_Fingerprint_uchar result;
  EXPECT_TRUE(result.construct_from_array(values.begin(), values.size()));
  return result;
}

TEST(TestFixedSizeCountedFingerprint, CopyConstructor) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar copy(fp);

  EXPECT_EQ(copy.size(), fp.size());
  EXPECT_EQ(copy.nset(), fp.nset());
  EXPECT_FLOAT_EQ(copy.tanimoto(fp), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, CopyAssignment) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar copy = Fingerprint({1, 0});

  copy = fp;

  EXPECT_EQ(copy.size(), fp.size());
  EXPECT_EQ(copy.nset(), fp.nset());
  EXPECT_FLOAT_EQ(copy.tanimoto(fp), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, MoveConstructor) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar moved(std::move(fp));

  Fixed_Size_Counted_Fingerprint_uchar expected = Fingerprint({1, 1, 0, 1});
  EXPECT_EQ(moved.size(), expected.size());
  EXPECT_EQ(moved.nset(), expected.nset());
  EXPECT_FLOAT_EQ(moved.tanimoto(expected), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, MoveAssignment) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar moved = Fingerprint({1, 0});

  moved = std::move(fp);

  Fixed_Size_Counted_Fingerprint_uchar expected = Fingerprint({1, 1, 0, 1});
  EXPECT_EQ(moved.size(), expected.size());
  EXPECT_EQ(moved.nset(), expected.nset());
  EXPECT_FLOAT_EQ(moved.tanimoto(expected), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, VectorReallocation) {
  std::vector<Fixed_Size_Counted_Fingerprint_uchar> fingerprints;
  fingerprints.reserve(1);

  fingerprints.push_back(Fingerprint({1, 0, 1, 0}));
  fingerprints.push_back(Fingerprint({0, 1, 0, 1}));
  fingerprints.push_back(Fingerprint({1, 0, 0, 1}));

  Fixed_Size_Counted_Fingerprint_uchar expected0 = Fingerprint({1, 0, 1, 0});
  Fixed_Size_Counted_Fingerprint_uchar expected1 = Fingerprint({0, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar expected2 = Fingerprint({1, 0, 0, 1});

  EXPECT_FLOAT_EQ(fingerprints[0].tanimoto(expected0), 1.0f);
  EXPECT_FLOAT_EQ(fingerprints[1].tanimoto(expected1), 1.0f);
  EXPECT_FLOAT_EQ(fingerprints[2].tanimoto(expected2), 1.0f);
}

}  // namespace
