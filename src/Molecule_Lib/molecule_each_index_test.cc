// Tester for the each_index member function

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"

namespace {

TEST(TestEachIndex, test_ncon) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CC"));
  Molecule::const_member_fn<int> fn = &Molecule::ncon;
  EXPECT_EQ(m.each_index<int>(fn), 2);
}

TEST(TestEachIndex, test_lambda_ncon) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CC"));
  int sum = 0;
  m.each_index_lambda([&m, &sum](int i) {sum += m.ncon(i);});
  EXPECT_EQ(sum, 2);
}

TEST(TestEachIndex, test_lambda_charge) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("NN"));
  m.each_index_lambda([&m](int i) {m.set_formal_charge(i, 1);});
  EXPECT_EQ(m.unique_smiles(), "[NH3+][NH3+]");
}


}  // namespace
