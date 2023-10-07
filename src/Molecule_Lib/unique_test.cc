// Tester for unique smiles generation.

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "aromatic.h"
#include "molecule.h"
#include "smiles.h"

namespace {

// A smiles and an expected unique smiles.
struct SmilesAndExpected {
  IWString smiles;
  IWString expected;
};

class TestUniqueSmiles: public testing::TestWithParam<SmilesAndExpected> {
  protected:
    void SetUp();
};

void
TestUniqueSmiles::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestUniqueSmiles, Lots) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
//std::cerr << "usmi " << m.unique_smiles() << std::endl;
  EXPECT_EQ(m.unique_smiles(), params.expected);
}

INSTANTIATE_TEST_SUITE_P(RandomlySelectedMolecules, TestUniqueSmiles, testing::Values(
   SmilesAndExpected{"C", "C"},
   SmilesAndExpected{"CC", "CC"},
   SmilesAndExpected{"OC", "OC"},
   SmilesAndExpected{"CO", "OC"},
   SmilesAndExpected{"O=C(O)C12CC3(CC(C1)CC(C2)C3)", "O=C(O)C12CC3CC(CC(C3)C1)C2"}
));

class TestUniqueKekuleSmiles: public testing::TestWithParam<SmilesAndExpected> {
  protected:
    Molecule _m;

    void SetUp();
};

void
TestUniqueKekuleSmiles::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestUniqueKekuleSmiles, Kekule) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.UniqueKekuleSmiles(), params.expected);
}

INSTANTIATE_TEST_SUITE_P(AromaticMolecules, TestUniqueKekuleSmiles, testing::Values(
  SmilesAndExpected{"C", "C"},
  SmilesAndExpected{"CC", "CC"},
  SmilesAndExpected{"C1CC1", "C1CC1"},
  SmilesAndExpected{"c1ccccc1", "c1=cc=cc=c1"},
  SmilesAndExpected{"N1=CC=CC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=NC=CC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CN=CC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CC=NC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CC=CN=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CC=CC=N1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=C2C=CC=CC2=NC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C12=CC=CN=C2C=CC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C1=C2C=CC=CC2=CN=C1", "[n]1=cc2=c(c=c1)c=cc=c2"},
  SmilesAndExpected{"C1=CC=CC2=CC=NC=C12", "[n]1=cc2=c(c=c1)c=cc=c2"}
));

class TestUniqueKekuleSmilesRandom: public testing::TestWithParam<SmilesAndExpected> {
  protected:
    Molecule _m;

  void SetUp();
};

void
TestUniqueKekuleSmilesRandom::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestUniqueKekuleSmilesRandom, Kekule) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.UniqueKekuleSmiles(), params.expected);
  constexpr int kNtest = 10;
  for (int i = 0; i < kNtest; ++i) {
    const IWString& smi = _m.random_smiles();
    Molecule mcopy;
    ASSERT_TRUE(mcopy.build_from_smiles(smi));
    ASSERT_EQ(mcopy.unique_smiles(), _m.unique_smiles());
    EXPECT_EQ(mcopy.UniqueKekuleSmiles(), params.expected);
  }
}

INSTANTIATE_TEST_SUITE_P(KekuleForms, TestUniqueKekuleSmilesRandom, testing::Values(
  SmilesAndExpected{"C1=C2C=CC=CC2=NC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C12=CC=CN=C2C=CC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C1=C2C=CC=CC2=CN=C1", "[n]1=cc2=c(c=c1)c=cc=c2"},
  SmilesAndExpected{"C1=CC=CC2=CC=NC=C12", "[n]1=cc2=c(c=c1)c=cc=c2"},
  SmilesAndExpected{"N1C=CC2=C1C=CC=C2", "[nH]1c=cc2=cc=cc=c12"},  // indole
  SmilesAndExpected{"N1=CN=CC=C1", "c1=c[n]=c[n]=c1"},   // pyrimidine
  SmilesAndExpected{"N1=CNC2=C1C=NC=N2", "[n]1=cc2=c([n]=c1)[nH]c=[n]2"},
  SmilesAndExpected{"C1=CC=CC2=C1CC3=C(O2)C=CC=C3", "O1c2=c(Cc3=cc=cc=c13)c=cc=c2"},  // xanthene
  SmilesAndExpected{"C1=CC=CC2=C1NC3=C2C=CC=C3", "[nH]1c2=c(c=cc=c2)c2=c1c=cc=c2"}  // carbazole
));

}  // namespace
