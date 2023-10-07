// Tester for PartialSymmetry.

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "partial_symmetry.h"

namespace {

using partial_symmetry::PartialSymmetry;

class TestPartialSymmetry : public testing::Test {
  protected:
    Molecule _m;
};

struct SmilesAndExpected {
  IWString smiles;
  resizable_array<int> expected;
};

class TestPartialSymmetryP : public testing::TestWithParam<SmilesAndExpected> {
  protected:
    Molecule _m;
};

TEST_P(TestPartialSymmetryP, TestMolecule) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  PartialSymmetry p(_m);
  const int * symmetric = p.SymmetricAtRadius();
  EXPECT_TRUE(symmetric);
//cerr << "Processing " << params.smiles << '\n';
  for (int i = 0; i < params.expected.number_elements(); ++i) {
//  cerr << "Checking atom " << i << " expected " << params.expected[i] << " got " << symmetric[i] << '\n';
    EXPECT_EQ(params.expected[i], symmetric[i]);
  }
}
INSTANTIATE_TEST_SUITE_P(TestPartialSymmetryP, TestPartialSymmetryP, testing::Values(
   SmilesAndExpected{"C", {0}},
   SmilesAndExpected{"CN", {0, 0}},
   SmilesAndExpected{"CC", {0, 0}},
   SmilesAndExpected{"CCC", {1, 0, 1}},
   SmilesAndExpected{"OCN", {0, 0, 0}},
   SmilesAndExpected{"C1CC1", {1, 1, 1}},
   SmilesAndExpected{"CCCC", {1, 1, 1, 1}},
   SmilesAndExpected{"C1CCC1", {1, 1, 1, 1}},
   SmilesAndExpected{"CC(C)(C)C", {1, 0, 1, 1, 1}},
   SmilesAndExpected{"CCCCC", {2, 1, 1, 1, 2}},
   SmilesAndExpected{"C1CCCC1", {2, 2, 2, 2, 2}},
   SmilesAndExpected{"CCCCCCC", {3, 2, 2, 1, 2, 2, 3}},
   SmilesAndExpected{"CCCCCCCCC", {4, 3, 3, 2, 2, 2, 3, 3, 4}},
   SmilesAndExpected{"NCCCCCO", {0, 2, 1, 1, 1, 2, 0}},
   SmilesAndExpected{"NCCC(CCC)CCO", {0, 3, 2, 1, 2, 3, 2, 2, 3, 0}},
   SmilesAndExpected{"NCCC(CCC)(CCO)CCF", {0, 3, 2, 1, 2, 3, 2, 2, 3, 0, 2, 3, 0}},
   SmilesAndExpected{"CCC1CCCC1", {1, 1, 2, 2, 3, 3, 2}}
));

TEST_F(TestPartialSymmetry, TestEmptyMolecule) {
  PartialSymmetry p(_m);
  EXPECT_FALSE(p.SymmetricAtRadius());
}

}  // namespace
