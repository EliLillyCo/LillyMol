// Tests for ladderane and related

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "ladderane.h"

namespace {

struct Data {
  IWString smiles;
  int result;
};

class LadderaneTest : public testing::TestWithParam<Data> {
  protected:
    Molecule _m;
};

TEST_P(LadderaneTest, LadderaneTest) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  EXPECT_EQ(ladderane::CountLadderane(_m), params.result) << params.smiles << ' ' << _m.name()
                << " expect " << params.result;
}
INSTANTIATE_TEST_SUITE_P(LadderaneTest, LadderaneTest, testing::Values(
  Data{ "C1CCC1 t1", 0},
  Data{ "C1CC2CCC12 t2", 1},
  Data{ "C1CCCC2CCC12 t3", 0},
  Data{ "C1CCCC2CCC12.C1CCC1 t4", 0},
  Data{ "C1CCCC2C3CCC3C12 t5", 1},
  Data{ "C1CCCC2C3C4CCC4C3C21 t6", 2},
  Data{ "C1CCCC2C3C4CC5CCC45C3C21 t7", 3},
  Data{ "C12C3C4C1C5C2C3C45 cubane t8", 0},
  Data{ "C1CC2CCC12.C1CC2C3CCC3C12 t9", 2},
  Data{ "C1CCC1.C1CCC1 t10", 0},
  Data{ "CC12CC(C1)(C2)C", 0}
));

class PolySpiroCycloPropaneTest : public testing::TestWithParam<Data> {
  protected:
    Molecule _m;
};

TEST_P(PolySpiroCycloPropaneTest, PolySpiroCycloPropaneTest) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  EXPECT_EQ(ladderane::CountPolySpiroCycloPropane(_m), params.result) << params.smiles << ' ' << _m.name()
                << " expect " << params.result;
}
INSTANTIATE_TEST_SUITE_P(PolySpiroCycloPropaneTest, PolySpiroCycloPropaneTest, testing::Values(
  Data{ "C1CC1 t1", 0},
  Data{ "C1CC12CC2 t2", 0},
  Data{ "C1CC12CC23CC3 t3", 3},
  Data{ "C1CC12CC23CC3.C1CC12CC2 t3", 3}
));

}  // namespace

