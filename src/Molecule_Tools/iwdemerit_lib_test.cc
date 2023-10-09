// Tests for iwdemerit_lib

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "iwdemerit_lib.h"

namespace {

using lilly_medchem_rules::MCDemerit;

class MCDemeritTest : public testing::Test {
  protected:
    Molecule _m;

    MCDemerit _mcdemerit;
};

TEST_F(MCDemeritTest, TestHardLowerAtomCountCutoffLower) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=4", "-c", "smax=8", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_TRUE(demerit.rejected());
  EXPECT_EQ(demerit.demerits().size(), 1);
  EXPECT_EQ(demerit.demerits()[0]->reason(), "too_few_atoms");
}

TEST_F(MCDemeritTest, TestHardLowerAtomCountCutoffEquals) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=4", "-c", "smax=8", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_TRUE(demerit.rejected());
}
TEST_F(MCDemeritTest, TestHardLowerAtomCountCutoffGreaterZero) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=4", "-c", "smax=8", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_FALSE(demerit.rejected());
  EXPECT_FLOAT_EQ(demerit.score(), 0.0);
}
TEST_F(MCDemeritTest, TestHardLowerAtomCountCutoffGreaterDemerited) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=8", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_FALSE(demerit.rejected());
  EXPECT_EQ(demerit.score(), 50);
}
TEST_F(MCDemeritTest, TestSoftLowerAtomCountCutoffGreaterNoDemerit) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=8", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_FALSE(demerit.rejected());
  EXPECT_EQ(demerit.score(), 0);
}
TEST_F(MCDemeritTest, TestMidRangeNoDemerit) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=8", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_FALSE(demerit.rejected());
  EXPECT_EQ(demerit.score(), 0);
}
TEST_F(MCDemeritTest, TestSoftUpperNoDemerit) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=7", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCCCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_FALSE(demerit.rejected());
  EXPECT_EQ(demerit.score(), 0);
}
TEST_F(MCDemeritTest, TestDemeritedUpper) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=7", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCCCCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_FALSE(demerit.rejected());
  EXPECT_EQ(demerit.score(), 50);
}
TEST_F(MCDemeritTest, TestRejectedAtUpper) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=7", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCCCCCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_TRUE(demerit.rejected());
  EXPECT_EQ(demerit.demerits().size(), 1);
  EXPECT_EQ(demerit.demerits()[0]->reason(), "too_many_atoms");
}
TEST_F(MCDemeritTest, TestAboveUpperThreshold) {
  const char* opts[]{"notused", "-c", "hmin=3", "-c", "smin=5", "-c", "smax=7", "-c", "hmax=9", "-v"};
  Command_Line cl(10, const_cast<char**>(opts), "vc:");
  ASSERT_TRUE(_mcdemerit.Build(cl));
  ASSERT_TRUE(_m.build_from_smiles("CCCCCCCCCCCCCC"));
  Demerit demerit = _mcdemerit.Process(_m, 0);
  EXPECT_TRUE(demerit.rejected());
  EXPECT_EQ(demerit.demerits().size(), 1);
  EXPECT_EQ(demerit.score(), 130);
  EXPECT_EQ(demerit.demerits()[0]->reason(), "too_many_atoms");
}
}  // namespace
