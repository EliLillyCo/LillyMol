// Tests for Substructure_Bond

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "aromatic.h"
#include "substructure.h"

namespace {

class TestSSBond : public testing::Test {
  protected:
    Substructure_Bond bond;
    const char* smarts;
    int chars_to_process;
    int characters_processed;
    IWString returned_smarts;
};

TEST_F(TestSSBond, FailIfInputEmpty) {
  smarts = "";
  chars_to_process = 0;
  EXPECT_FALSE(bond.construct_from_smarts(smarts, 0, characters_processed));
}

TEST_F(TestSSBond, SingleBond) {
  smarts = "-";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, 1);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}
TEST_F(TestSSBond, DoubleBond) {
  smarts = "=";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, 1);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, TripleBond) {
  smarts = "#";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, 1);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, AromaticBond) {
  smarts = ":";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, 1);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleOrDouble) {
  smarts = "-,=";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleOrArom) {
  smarts = "-,:";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleOrDoubleOrTriple) {
  smarts = "-,=,#";
  chars_to_process = strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, Ring) {
  smarts = "@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, "@");
}

TEST_F(TestSSBond, NotRing) {
  smarts = "!@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleAndRingImplicitAnd) {
  smarts = "-@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, "-&@");
}

TEST_F(TestSSBond, SingleAndRingExplicitAnd) {
  smarts = "-&@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleAndRingLowPriorityAnd) {
  smarts = "-;@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleNotRing) {
  smarts = "-!@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, "-&!@");
}

TEST_F(TestSSBond, DoubleTripleOrRing) {
  smarts = "=,#,@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, DoubleTripleHpAndRing) {
  smarts = "=,#&@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, DoubleTripleLpAndRing) {
  smarts = "=,#;@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts);
}

TEST_F(TestSSBond, SingleOrNotDoubleOrNotTriple) {
  smarts = "-,!=,!#";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts) << " mismatch '" << returned_smarts << "'";
}

TEST_F(TestSSBond, SingleLpAndNotRing) {
  smarts = "-;!@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts) << " mismatch '" << returned_smarts << "'";
}

TEST_F(TestSSBond, SingleHpAndNotRing) {
  smarts = "-&!@";
  chars_to_process = ::strlen(smarts);
  EXPECT_TRUE(bond.construct_from_smarts(smarts, chars_to_process, characters_processed));
  EXPECT_EQ(characters_processed, chars_to_process);
  returned_smarts = bond.Smarts();
  EXPECT_EQ(returned_smarts, smarts) << " mismatch '" << returned_smarts << "'";
}


}  // namespace
