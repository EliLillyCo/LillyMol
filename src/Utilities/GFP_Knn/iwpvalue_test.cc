#include "Utilities/GFP_Knn/iwpvalue.h"

#include "gtest/gtest.h"

namespace {

TEST(IwPValue, ZeroTValueHasUnitProbability) {
  EXPECT_DOUBLE_EQ(iwpvalue(10, 0.0), 1.0);
}

TEST(IwPValue, SymmetricInTValue) {
  EXPECT_DOUBLE_EQ(iwpvalue(7, 2.0), iwpvalue(7, -2.0));
}

TEST(IwPValue, KnownTwoSidedCriticalValues) {
  EXPECT_NEAR(iwpvalue(1, 12.7062047361747), 0.05, 1.0e-10);
  EXPECT_NEAR(iwpvalue(10, 2.22813885196494), 0.05, 1.0e-10);
  EXPECT_NEAR(iwpvalue(30, 2.04227245630124), 0.05, 1.0e-10);
}

TEST(IwPValue, CauchySpecialCase) {
  // Student t with one degree of freedom is the Cauchy distribution.
  EXPECT_NEAR(iwpvalue(1, 1.0), 0.5, 1.0e-12);
}

}  // namespace
