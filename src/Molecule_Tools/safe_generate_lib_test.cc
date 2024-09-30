// Tests for safe_generate_library

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/safe_generate_lib.h"

namespace {

using safe_generate::BondsOk;

TEST(TestBondsOk, Test1) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("N#[1C][1CH]1[1CH]([1NH][1CH3])[1CH]([1S][1CH2][1CH]2NC3C(N2)CCCC3)[1CH]([1CH]2CCCCC2)C=C1"));
  EXPECT_TRUE(BondsOk(m));
}

struct SmilesResult {
  IWString smiles;
  int result;
};

class TestBondsOk: public testing::TestWithParam<SmilesResult> {
  protected:
    Molecule _m;
};
TEST_P(TestBondsOk, Tests) {
  auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(BondsOk(_m), params.result) << params.smiles;
}
INSTANTIATE_TEST_SUITE_P(TestBondsOk, TestBondsOk, testing::Values(
  SmilesResult{"[1Cl]-[1CH3]", 0},
  SmilesResult{"[1Br]-[1CH3]", 0},
  SmilesResult{"[1I]-[1CH3]", 0},
  SmilesResult{"[1F]-[1CH3]", 1},
  SmilesResult{"[1F]-[1c]1ccccc1", 1},
  SmilesResult{"[1Cl]-[1c]1ccccc1", 1},
  SmilesResult{"[1Br]-[1c]1ccccc1", 1},
  SmilesResult{"[1I]-[1c]1ccccc1", 1},
  SmilesResult{"C[1NH][1F]", 0},
  SmilesResult{"C[1NH][1Cl]", 0},
  SmilesResult{"C[1NH][1Br]", 0},
  SmilesResult{"C[1NH][1I]", 0},
  SmilesResult{"C[1O][1F]", 0},
  SmilesResult{"C[1O][1Cl]", 0},
  SmilesResult{"C[1O][1Br]", 0},
  SmilesResult{"C[1O][1I]", 0},
  SmilesResult{"C[1O]-[1O]C", 0},
  SmilesResult{"C[1NH]-[1NH]C", 0},
  SmilesResult{"C[1S]-[1S]C", 0},
  SmilesResult{"C[1S]-[1CH2]C", 1},
  SmilesResult{"C[1S](=O)-[1NH]C", 1},
  SmilesResult{"C[1S][1NH]C", 0},
  SmilesResult{"C[1S]-[1P](=O)OC", 0},
  SmilesResult{"C[1O][1O]C", 0},
  // In a ring is not checked
  SmilesResult{"C1[1O][1O]C1", 1},
  SmilesResult{"C[1NH][1NH]C", 0},
  SmilesResult{"C1[1NH][1NH]C1", 1},
  SmilesResult{"N#[1C][1P]1(=O)O[1P](=O)([1OH])OC[1CH]([1OH])[1CH]([1CH2][1OH])O1", 1},
  SmilesResult{"[1Cl][1CH]1[1N]([1C]2=NO[1C]([1CH]3NCCCC3)([1CH3])O2)CCC1", 0}
));



TEST(TestBondsOk, Test2) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("N#CC1C(NC)C(SCC2NC3C(N2)CCCC3)C(C2CCCCC2)C=C1"));
  EXPECT_TRUE(BondsOk(m));
}



}  // namespace
