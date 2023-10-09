// Tester for memoized floats.

#include <random>

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"
//#include <benchmark/benchmark.h>


#include "Foundational/iwmisc/memoized_floats.h"

namespace memoized_floats {
namespace {

TEST(MemoizedFloats, TestPositiveWholeNumbers) {
  MemoizedFloats mf;
  constexpr int kMax = 10;
  ASSERT_TRUE(mf.Build(kMax, 2));
  for (int i = 0; i < kMax; ++i) {
    const IWString s = mf.Representation(static_cast<float>(i));
    IWString expected;
    expected << i << '.' << "0";
    EXPECT_EQ(s, expected);
  }
}
TEST(MemoizedFloats, TestNegativeWholeNumbers) {
  MemoizedFloats mf;
  constexpr int kMax = 10;
  ASSERT_TRUE(mf.Build(kMax, 2));
  for (int i = 1; i < kMax; ++i) {
    const IWString s = mf.Representation(static_cast<float>(-i));
    IWString expected;
    expected << '-' << i << '.' << "0";
    EXPECT_EQ(s, expected);
  }
}

TEST(MemoizedFloats, TestPositiveFractions1) {
  MemoizedFloats mf;
  ASSERT_TRUE(mf.Build(5, 1));

  for (int i = 0; i < 10; ++i) {
    const IWString s = mf.Representation(static_cast<float>(i) / 10.0f);
    IWString expected;
    expected << "0." << i;
    EXPECT_EQ(s, expected);
  }
}
TEST(MemoizedFloats, TestNegativeFractions1) {
  MemoizedFloats mf;
  ASSERT_TRUE(mf.Build(5, 1));

  for (int i = 1; i < 10; ++i) {
    const IWString s = mf.Representation(static_cast<float>(-i) / 10.0f);
    IWString expected;
    expected << "-0." << i;
    EXPECT_EQ(s, expected);
  }
}

TEST(MemoizedFloats, TestpositiveFractionRoundDown) {
  MemoizedFloats mf;
  ASSERT_TRUE(mf.Build(5, 1));

  float value = 0.01f;
  for (int i = 0; i < 10; ++i, value += 0.1f) {
    const IWString s = mf.Representation(value);
    IWString expected;
    expected << "0." << i;
    EXPECT_EQ(s, expected);
  }
}

TEST(MemoizedFloats, TestpositiveFractionRoundUp) {
  MemoizedFloats mf;
  ASSERT_TRUE(mf.Build(5, 1));

  float value = 0.099f;
  for (int i = 0; i < 10; ++i, value += 0.1f) {
    const IWString s = mf.Representation(value);
    IWString expected;
    expected << "0." << (i + 1);
    EXPECT_EQ(s, expected);
  }
}

}  // namespace
}  // namespace memoized_floats
