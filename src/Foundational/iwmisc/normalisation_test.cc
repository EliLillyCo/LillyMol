// Tester for Normalisation

// Acknowledge, never a good idea to have random numbers in tests.
#include <random>
#include <vector>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "Foundational/accumulator/accumulator.h"

#include "normalisation.h"

#ifdef BUILD_BAZEL
#include "Foundational/iwmisc/normalisation.pb.h"
#else
#include "normalisation.pb.h"
#endif

namespace {

using normalisation::NormalisationData;
using normalisation::Normalisation;

// Return a vector of `n` random numbers between `minval` and `maxval`.
std::vector<double>
VectorOfRandomNumbers(int n, double minval, double maxval) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dist(minval, maxval);

  std::vector<double>result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = dist(gen);
  }

  return result;
}

// Return an Accumulator built on `values`.
Accumulator<double>
GetAccumulator(const std::vector<double>& values) {
  Accumulator<double> result;
  for (double v : values) {
    result.extra(v);
  }

  return result;
}


TEST(Test01, TestInRange01) {
  const int npoints = 100;

  std::vector<double> values = VectorOfRandomNumbers(npoints, 50.0, 100.0);

  Accumulator<double> acc = GetAccumulator(values);

  Normalisation mynorm;
  mynorm.SetScalingType(normalisation::NormalisationType::k01);
  ASSERT_TRUE(mynorm.Build(acc));

  for (double v : values) {
    double scaled = v;
    EXPECT_TRUE(mynorm.Scale(scaled));
    double unscaled = scaled;
    EXPECT_TRUE(mynorm.Unscale(unscaled));
    EXPECT_NEAR(v, unscaled, 0.001) << "original " << v << " scaled " << scaled << " unscaled " << unscaled;
  }
}

TEST(Test01, TestInRange11) {
  const int npoints = 100;

  std::vector<double> values = VectorOfRandomNumbers(npoints, -100.0, -50.0);

  Accumulator<double> acc = GetAccumulator(values);

  Normalisation mynorm;
  mynorm.SetScalingType(normalisation::NormalisationType::k11);
  ASSERT_TRUE(mynorm.Build(acc));

  for (double v : values) {
    double scaled = v;
    EXPECT_TRUE(mynorm.Scale(scaled));
    double unscaled = scaled;
    EXPECT_TRUE(mynorm.Unscale(unscaled));
    EXPECT_NEAR(v, unscaled, 0.001) << "original " << v << " scaled " << scaled << " unscaled " << unscaled;
  }
}


TEST(Test01, TestInRange11NoTruncation) {
  const int npoints = 100;

  std::vector<double> values = VectorOfRandomNumbers(npoints, -100.0, -50.0);

  Accumulator<double> acc = GetAccumulator(values);

  Normalisation mynorm;
  mynorm.SetScalingType(normalisation::NormalisationType::k11);
  ASSERT_TRUE(mynorm.Build(acc));
  mynorm.set_truncate_out_of_range_scale_requests(0);
  mynorm.set_truncate_out_of_range_unscale_requests(0);

  double x = -101.0;
  double scaled = x;
  EXPECT_TRUE(mynorm.Scale(scaled));
  double unscaled = scaled;
  EXPECT_TRUE(mynorm.Unscale(unscaled));
  EXPECT_NEAR(x, unscaled, 0.001) << "begin " << x << " scaled " << scaled << " unscaled " << unscaled;

  x = -45.0;
  scaled = x;
  EXPECT_TRUE(mynorm.Scale(scaled));
  unscaled = scaled;
  EXPECT_TRUE(mynorm.Unscale(unscaled));
  EXPECT_NEAR(x, unscaled, 0.001) << "begin " << x << " scaled " << scaled << " unscaled " << unscaled;
}

TEST(Test01, TestInRange11WithTruncation) {
  const int npoints = 100;

  std::vector<double> values = VectorOfRandomNumbers(npoints, -100.0, 50.0);

  Accumulator<double> acc = GetAccumulator(values);

  Normalisation mynorm;
  mynorm.SetScalingType(normalisation::NormalisationType::k11);
  ASSERT_TRUE(mynorm.Build(acc));

  double x = -101.0;
  double scaled = x;
  EXPECT_TRUE(mynorm.Scale(scaled));
  double unscaled = scaled;
  EXPECT_TRUE(mynorm.Unscale(unscaled));
  EXPECT_NEAR(acc.minval(), unscaled, 0.001) << "begin " << x << " scaled " << scaled << " unscaled " << unscaled;

  x = 55.0;
  scaled = x;
  EXPECT_TRUE(mynorm.Scale(scaled));
  unscaled = scaled;
  EXPECT_TRUE(mynorm.Unscale(unscaled));
  EXPECT_NEAR(acc.maxval(), unscaled, 0.001) << "begin " << x << " scaled " << scaled << " unscaled " << unscaled;
}

TEST(Test01, TestMeanZero) {
  const int npoints = 100;

  std::vector<double> values = VectorOfRandomNumbers(npoints, -100.0, 50.0);

  Accumulator<double> acc = GetAccumulator(values);

  Normalisation mynorm;
  mynorm.SetScalingType(normalisation::NormalisationType::kZeroMean);
  ASSERT_TRUE(mynorm.Build(acc));

  Accumulator<double> acc_rescaled;

  for (double v : values) {
    double scaled = v;
    mynorm.Scale(scaled);
    acc_rescaled.extra(scaled);
  }

  EXPECT_NEAR(acc_rescaled.average(), 0.0, 0.0001);
}

}  // namespace
