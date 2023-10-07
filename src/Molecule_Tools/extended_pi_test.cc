// Tester for extended_pi.

#include <algorithm>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "extended_pi.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

namespace {

struct ExtendedPiArgs {
  IWString smiles;
  int center_atom;
  int radius;
  int expected_count;
};

class TestExtendedPiTest : public testing::TestWithParam<ExtendedPiArgs> {
};

TEST_P(TestExtendedPiTest, LinearMolecules) {
  const ExtendedPiArgs test_params = GetParam();
  radius_pi_extension::RadiusPiExtensionParams params;
  params.radius = test_params.radius;
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(test_params.smiles));
  int * in_system = new_int(m.natoms());
  std::unique_ptr<int[]> free_in_system(in_system);
  in_system[test_params.center_atom] = 1;
  EXPECT_EQ(radius_pi_extension::RadiusPiExtension(m, in_system, params), test_params.expected_count);
}

// Could do something clever with multiple iterators through values...
INSTANTIATE_TEST_SUITE_P(TestLinearMolecules, TestExtendedPiTest, testing::Values(
  ExtendedPiArgs{"CCC", 1, 0, 1},   // 0
  ExtendedPiArgs{"CCC", 1, 1, 3},   // 1
  ExtendedPiArgs{"CCC", 0, 1, 2},   // 2
  ExtendedPiArgs{"CCC", 1, 2, 3},   // 3
  ExtendedPiArgs{"CCCCC", 2, 1, 3},   // 4
  ExtendedPiArgs{"CCCCC", 2, 2, 5},   // 5
  ExtendedPiArgs{"CCCCC", 2, 3, 5},   // 6
  ExtendedPiArgs{"CCCCC", 1, 3, 5},   // 7
  ExtendedPiArgs{"CCCCC", 0, 3, 4},   // 8
  ExtendedPiArgs{"CCCCCCC", 0, 1, 2},   // 9
  ExtendedPiArgs{"CCCCCCC", 0, 2, 3},   // 10
  ExtendedPiArgs{"CCCCCCC", 0, 3, 4},   // 11
  ExtendedPiArgs{"CCCCCCC", 0, 4, 5},   // 12
  ExtendedPiArgs{"CCCCCCC", 0, 5, 6},   // 13
  ExtendedPiArgs{"CCCCCCC", 0, 6, 7},   // 14
  ExtendedPiArgs{"CCCCCCC", 3, 1, 3},   // 15
  ExtendedPiArgs{"CCCCCCC", 3, 2, 5}   // 16
));

TEST(test_extended_pi, TestNothingToDo) {
  radius_pi_extension::RadiusPiExtensionParams params;
  params.radius = 0;
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC"));
  int * in_system = new_int(m.natoms()); std::unique_ptr<int[]> free_in_system(in_system);
  in_system[1] = 1;
  EXPECT_EQ(radius_pi_extension::RadiusPiExtension(m, in_system, params), 1);
}

INSTANTIATE_TEST_SUITE_P(TestBranched, TestExtendedPiTest, testing::Values(
  ExtendedPiArgs{"CC(C)(C)C", 0, 1, 2},
  ExtendedPiArgs{"CC(C)(C)C", 1, 1, 5},
  ExtendedPiArgs{"CC(C)(C)CC", 1, 1, 5},
  ExtendedPiArgs{"CC(C)(C)CC", 0, 3, 6}
));

// Could make this all one struct, but seems desirable to separate tests.
struct ExtendedPiArgsArom {
  IWString smiles;
  int center_atom;
  int radius;
  int extend_aromatic;
  int expected_count;
};

class TestExtendedPiArom : public testing::TestWithParam<ExtendedPiArgsArom> {
};

TEST_P(TestExtendedPiArom, TestArom) {
  const ExtendedPiArgsArom test_params = GetParam();
  radius_pi_extension::RadiusPiExtensionParams params;
  params.radius = test_params.radius;
  params.aromatic_bond_extend = test_params.extend_aromatic;
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(test_params.smiles));
  int * in_system = new_int(m.natoms());
  std::unique_ptr<int[]> free_in_system(in_system);
  in_system[test_params.center_atom] = 1;
  EXPECT_EQ(radius_pi_extension::RadiusPiExtension(m, in_system, params), test_params.expected_count);
}

INSTANTIATE_TEST_SUITE_P(TestARom, TestExtendedPiArom, testing::Values(
  ExtendedPiArgsArom{"CC(C)(C)C", 0, 1, 0, 2},
  ExtendedPiArgsArom{"Cc1ccccc1", 0, 1, 0, 2},
  ExtendedPiArgsArom{"Cc1ccccc1", 0, 1, 1, 4},
  ExtendedPiArgsArom{"Cc1ccccc1", 0, 1, 2, 6},
  ExtendedPiArgsArom{"Cc1ccccc1", 0, 1, 3, 7}
));

// Could make this all one struct, but seems desirable to separate tests.
struct ExtendedPiArgsConjugated {
  IWString smiles;
  int center_atom;
  int radius;
  int extend_conjugated;
  int expected_count;
};

class TestExtendedPiConjugated : public testing::TestWithParam<ExtendedPiArgsConjugated> {
};

TEST_P(TestExtendedPiConjugated, TestConjugated) {
  const ExtendedPiArgsConjugated test_params = GetParam();
  radius_pi_extension::RadiusPiExtensionParams params;
  params.radius = test_params.radius;
  params.conjugated_bond_extend = test_params.extend_conjugated;
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(test_params.smiles));
  int * in_system = new_int(m.natoms());
  std::unique_ptr<int[]> free_in_system(in_system);
  in_system[test_params.center_atom] = 1;
  EXPECT_EQ(radius_pi_extension::RadiusPiExtension(m, in_system, params), test_params.expected_count);
}

INSTANTIATE_TEST_SUITE_P(TestConjugated, TestExtendedPiConjugated, testing::Values(
  ExtendedPiArgsConjugated{"CC(C)(C)C", 0, 1, 0, 2},
  ExtendedPiArgsConjugated{"Cc1ccccc1", 0, 1, 0, 2},
  ExtendedPiArgsConjugated{"Cc1ccccc1", 0, 1, 3, 7},
  ExtendedPiArgsConjugated{"CC(=O)OC", 0, 1, 1, 4},
  ExtendedPiArgsConjugated{"CC(=O)OC", 0, 1, 2, 4},
  ExtendedPiArgsConjugated{"CC(=O)Oc1ccccc1", 0, 1, 2, 7},
  ExtendedPiArgsConjugated{"CC#N", 0, 1, 1, 3}
));

TEST(AromHidesConjugated, SimpleCase) {
  radius_pi_extension::RadiusPiExtensionParams params;
  params.radius = 1;
  params.conjugated_bond_extend = 5;
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccc(C(=O)NC)cc1"));
  int * in_system = new_int(m.natoms());
  std::unique_ptr<int[]> free_in_system(in_system);
  in_system[0] = 1;
  EXPECT_EQ(radius_pi_extension::RadiusPiExtension(m, in_system, params), 10);

  std::fill_n(in_system, m.natoms(), 0);
  in_system[0] = 1;
  params.aromatic_bond_extend =  3;
  EXPECT_EQ(radius_pi_extension::RadiusPiExtension(m, in_system, params), 7);
}

}   //  namespace
