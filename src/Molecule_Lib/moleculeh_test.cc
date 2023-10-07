// Tests for functions in moleculeh.cc

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"

namespace {

class TestMoveHydrogensLast : public testing::Test {
  protected:
    IWString _smiles;
    Molecule _m;
    IWString _initial_smiles, _final_smiles;
};

TEST_F(TestMoveHydrogensLast, TestNoHydrogens) {
  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  _initial_smiles = _m.unique_smiles();
  EXPECT_EQ(_m.move_hydrogens_to_end_of_connection_table(), 0);
  EXPECT_EQ(_m.smiles(), _initial_smiles);
}

TEST_F(TestMoveHydrogensLast, AlreadyAtEnd) {
  _smiles = "CCCCO[H]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  const int natoms = _m.natoms();
  ASSERT_EQ(_m.atomic_number(natoms - 1), 1);
  EXPECT_EQ(_m.move_hydrogens_to_end_of_connection_table(), 0);
  ASSERT_EQ(_m.atomic_number(natoms - 1), 1);
}


TEST_F(TestMoveHydrogensLast, AlreadyLast1) {
  _smiles = "CCCCO[H]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  const int matoms = _m.natoms();
  ASSERT_EQ(_m.atomic_number(matoms - 1), 1);
  EXPECT_EQ(_m.move_hydrogens_to_end_of_connection_table(), 0);
  EXPECT_EQ(_m.atomic_number(matoms - 1), 1);
}

TEST_F(TestMoveHydrogensLast, AlreadyLast2) {
  _smiles = "CCCCO[1H][2H]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  const int matoms = _m.natoms();
  ASSERT_EQ(_m.atomic_number(matoms - 1), 1);
  ASSERT_EQ(_m.atomic_number(matoms - 2), 1);
  EXPECT_EQ(_m.move_hydrogens_to_end_of_connection_table(), 0);
  EXPECT_EQ(_m.atomic_number(matoms - 1), 1);
  EXPECT_EQ(_m.atomic_number(matoms - 2), 1);
  EXPECT_EQ(_m.isotope(matoms - 2), 1);
  EXPECT_EQ(_m.isotope(matoms - 1), 2);
}

class TestMoveHydrogensLastP : public testing::TestWithParam<IWString> {
  protected:
    Molecule _m;
};

TEST_P(TestMoveHydrogensLastP, TestHMoved) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params));
  const int matoms = _m.natoms();
  EXPECT_EQ(_m.move_hydrogens_to_end_of_connection_table(), 1);
  EXPECT_EQ(_m.natoms(), matoms);
  EXPECT_EQ(_m.atomic_number(matoms - 1), 1);
}
INSTANTIATE_TEST_SUITE_P(TestMoveHydrogensLastP, TestMoveHydrogensLastP, testing::Values(
  IWString{"[H]CCCC"},
  IWString{"C[H]CCC"},
  IWString{"CC[H]CC"},
  IWString{"CCC[H]C"}
));

class TestMoveHydrogensLastOrderPreserved : public testing::TestWithParam<IWString> {
  protected:
    Molecule _m;
};

TEST_P(TestMoveHydrogensLastOrderPreserved, TestOrderPreserved) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params));
  int hydrogens = 0;
  const int matoms = _m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (_m.atomic_number(i) != 1) {
      continue;
    }
    ++hydrogens;
    _m.set_isotope(i, hydrogens);
  }
  _m.MoveToEndOfConnectionTable(1);
  int next_isotope_expected = 1;
  for (int i = 0; i < matoms; ++i) {
    if (_m.atomic_number(i) != 1) {
      continue;
    }
    EXPECT_EQ(_m.isotope(i), next_isotope_expected);
    ++next_isotope_expected;
  }
  for (int i = matoms - hydrogens; i < matoms; ++i) {
    EXPECT_EQ(_m.atomic_number(i), 1);
  }
}
INSTANTIATE_TEST_SUITE_P(TestMoveHydrogensLastOrderPreserved, TestMoveHydrogensLastOrderPreserved, testing::Values(
  IWString{"[H][H]CC"},
  IWString{"[H]C[H]CC"},
  IWString{"[H]CC[H]C"},
  IWString{"[H]CCC[H]"},
  IWString{"C[H]CC[H]"},
  IWString{"CC[H]C[H]"},
  IWString{"CC[H][H]C"}
));

class TestHydrogensLastChirality : public testing::TestWithParam<IWString> {
  protected:
    IWString _smiles;
    Molecule _m;
};

// This tests that the order of the Hydrogens is indeed
// preserved.
TEST_P(TestHydrogensLastChirality, WithChirality) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params));
  const IWString initial_usmi = _m.unique_smiles();
  _m.MoveToEndOfConnectionTable(1);
  EXPECT_EQ(_m.unique_smiles(), initial_usmi);
  const int matoms = _m.natoms();
  bool first_hydrogen_encountered = false;
  int previous_isotope = -1;
  for (int i = 0; i < matoms; ++i) {
    if (_m.atomic_number(i) != 1) {
      EXPECT_FALSE(first_hydrogen_encountered) << "Heavy atom after H";
      continue;
    }
    first_hydrogen_encountered = true;
    EXPECT_GT(_m.isotope(i), previous_isotope);
    previous_isotope = _m.isotope(i);
  }

}
INSTANTIATE_TEST_SUITE_P(TestHydrogensLastChirality, TestHydrogensLastChirality, testing::Values(
  IWString{"[1H][C@](C)(N)O"},
  IWString{"[1H][C@](C[2H])(N[3H])O[4H]"},
  IWString{"[1H][C@@](C[2H])(N[3H])O[4H]"},
  IWString{"[1H][2H][3H][4H]"},
  IWString{"[1H][2H][3H][4H]C"},
  IWString{"C[1H][2H][3H][4H]"},
  IWString{"C[1H]C[2H]C[3H]C[4H]C"},
  IWString{"CCC[1H]CCC[2H]CCC[3H]CCC[4H]CCC"}
));

TEST(TestImplicitHSatisfiesAltValence, Test1) {
  const IWString smiles("C=[SH]-C");
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));
  EXPECT_FALSE(m.ok_valence());
}
}  // namespace
