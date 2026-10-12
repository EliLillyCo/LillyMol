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
  EXPECT_EQ(_m.isotope(matoms - 2), static_cast<isotope_t>(1));
  EXPECT_EQ(_m.isotope(matoms - 1), static_cast<isotope_t>(2));
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
  isotope_t next_isotope_expected = 1;
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
  isotope_t previous_isotope = 0;
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
  EXPECT_FALSE(m.valence_ok());
}

struct SmilesHbaHbd {
  IWString smiles;
  int hba;
  int hbd;
};

class TestHbaHbd : public testing::TestWithParam<SmilesHbaHbd> {
  protected:
    Molecule _m;
};

TEST_P(TestHbaHbd, Default) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  int hba = _m.LipinskiNumHAcceptors();
  int hbd = _m.LipinskiNumHDonors();
  EXPECT_EQ(params.hba, hba);
  EXPECT_EQ(params.hbd, hbd);
}
// Donors are the sum of OHs and NHs - a count of hydrogens, not of the
// heteroatoms carrying them. Acceptors are every N and O, with no exclusions.
INSTANTIATE_TEST_SUITE_P(TestHbaHbd, TestHbaHbd, testing::Values(
  SmilesHbaHbd{"C", 0, 0},
  SmilesHbaHbd{"CN", 1, 2},
  SmilesHbaHbd{"COC", 1, 0},
  SmilesHbaHbd{"CO", 1, 1},
  SmilesHbaHbd{"O", 1, 2},
  SmilesHbaHbd{"N", 1, 3},
  SmilesHbaHbd{"NCCN", 2, 4},
  SmilesHbaHbd{"NCCO", 2, 3},
  // Explicit hydrogens must count the same as implicit ones.
  SmilesHbaHbd{"[NH2]C", 1, 2},
  // No acceptor exclusions - the amide N and both nitro oxygens all count.
  SmilesHbaHbd{"CC(N)=O", 2, 2},
  SmilesHbaHbd{"c1ccccc1[N+](=O)[O-]", 3, 0},
  // Pyrrole NH is one donor, pyridine N is an acceptor with no donor.
  SmilesHbaHbd{"c1cc[nH]c1", 1, 1},
  SmilesHbaHbd{"c1ccncc1", 1, 0},
  // Sulfur is neither, no matter how it is bonded.
  SmilesHbaHbd{"CS", 0, 0},
  SmilesHbaHbd{"CS(=O)(=O)N", 3, 2}
));

// The RDKit compatible counts. Values here are what RDKit 2026.03.3
// rdMolDescriptors.CalcNumHBA and CalcNumHBD return. These are deliberately
// NOT the Lipinski numbers above - note that sulfur participates, that amide
// nitrogen is not an acceptor, and that a primary amine is one donor rather
// than two.
class TestRDKitHbaHbd : public testing::TestWithParam<SmilesHbaHbd> {
  protected:
    Molecule _m;
};

TEST_P(TestRDKitHbaHbd, Default) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  EXPECT_EQ(params.hba, _m.RDKitNumHAcceptors()) << params.smiles;
  EXPECT_EQ(params.hbd, _m.RDKitNumHDonors()) << params.smiles;
}
INSTANTIATE_TEST_SUITE_P(TestRDKitHbaHbd, TestRDKitHbaHbd, testing::Values(
  SmilesHbaHbd{"C", 0, 0},
  // Donors count atoms here, not hydrogens - contrast the Lipinski cases.
  SmilesHbaHbd{"CN", 1, 1},
  SmilesHbaHbd{"NCCN", 2, 2},
  SmilesHbaHbd{"CO", 1, 1},
  SmilesHbaHbd{"COC", 1, 0},
  // Amide and sulfonamide nitrogen are excluded as acceptors.
  SmilesHbaHbd{"CC(N)=O", 1, 1},
  SmilesHbaHbd{"CNC(C)=O", 1, 1},
  SmilesHbaHbd{"CS(=O)(=O)N", 2, 1},
  // The hydroxyl of a carboxylic acid is not an acceptor, the carbonyl is.
  SmilesHbaHbd{"OC(=O)c1ccccc1", 1, 1},
  // A plain alcohol is.
  SmilesHbaHbd{"OCc1ccccc1", 1, 1},
  // Aromatic nitrogen must have no hydrogen and exactly two connections.
  SmilesHbaHbd{"c1ccncc1", 1, 0},
  SmilesHbaHbd{"c1cc[nH]c1", 0, 1},
  SmilesHbaHbd{"Cn1cccc1", 0, 0},
  // Sulfur participates, unlike the Lipinski counts.
  SmilesHbaHbd{"CSC", 1, 0},
  SmilesHbaHbd{"CS", 1, 1},
  SmilesHbaHbd{"c1ccsc1", 1, 0},
  // Nitro nitrogen is tetravalent, so not an acceptor. Both oxygens count.
  SmilesHbaHbd{"c1ccccc1[N+](=O)[O-]", 2, 0},
  // Explicit hydrogens must give the same answer as implicit ones.
  SmilesHbaHbd{"OCc1ccccc1", 1, 1},
  SmilesHbaHbd{"[OH]Cc1ccccc1", 1, 1},
  SmilesHbaHbd{"O([H])Cc1ccccc1", 1, 1},
  SmilesHbaHbd{"OC(=O)C", 1, 1},
  SmilesHbaHbd{"O([H])C(=O)C", 1, 1}
));
}  // namespace
