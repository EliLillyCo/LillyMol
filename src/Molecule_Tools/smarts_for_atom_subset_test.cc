#include <cstdint>
#include <vector>

#include "gtest/gtest.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/smiles.h"
#include "smarts_for_atom_subset.h"

namespace {

struct TestCase {
  const char* smiles;
  uint32_t atom_type;
  std::vector<int> include_atom;
  const char* expected;
};

class SmartsForAtomSubsetTest : public testing::TestWithParam<TestCase> {
 protected:
  void SetUp() override {
    set_global_aromaticity_type(Daylight);
    set_include_bond_aromaticity_in_smiles(0);
  }
};

TEST_P(SmartsForAtomSubsetTest, GeneratesExpectedSmarts) {
  const TestCase& test = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(test.smiles));
  ASSERT_EQ(test.include_atom.size(), static_cast<size_t>(m.natoms()));

  const IWString smarts = lillymol::SmartsForAtomSubset(
      m, test.atom_type, test.include_atom.data());

  EXPECT_EQ(smarts.AsString(), test.expected);
  EXPECT_EQ(include_bond_aromaticity_in_smiles(), 0);
}

const uint32_t kZac =
    IWATTYPE_USP_Z | IWATTYPE_USP_A | IWATTYPE_USP_C;

INSTANTIATE_TEST_SUITE_P(
    Subsets, SmartsForAtomSubsetTest,
    testing::Values(
        TestCase{"CCO", kZac, {1, 1, 1}, "[#6;D1][#6;D2][#8;D1]"},
        TestCase{"CC(C)O", kZac, {1, 1, 1, 1},
                 "[#6;D1][#6;D3]([#6;D1])[#8;D1]"},
        TestCase{"CC(C)O", kZac, {1, 0, 1, 1},
                 "[#6;D1].[#6;D1].[#8;D1]"},
        TestCase{"C=C", kZac, {1, 1}, "[#6;D1]=[#6;D1]"},
        TestCase{"c1ccccc1", kZac, {1, 1, 1, 1, 1, 1},
                 "[c;D2]1[c;D2][c;D2][c;D2][c;D2][c;D2]1"},
        TestCase{"c1ccccc1", kZac, {1, 1, 1, 0, 0, 0},
                 "[c;D2][c;D2][c;D2]"},
        TestCase{"C1CC1", kZac, {1, 1, 1},
                 "[#6;D2]1[#6;D2][#6;D2]1"}));

TEST(SmartsForAtomSubsetEdgeCases, EmptyAndNullInput) {
  Molecule empty;
  EXPECT_TRUE(
      lillymol::SmartsForAtomSubset(empty, IWATTYPE_USP_Z, nullptr).empty());

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC"));
  const int include_none[] = {0, 0};

  EXPECT_TRUE(
      lillymol::SmartsForAtomSubset(m, IWATTYPE_USP_Z, nullptr).empty());
  EXPECT_TRUE(lillymol::SmartsForAtomSubset(
                  m, IWATTYPE_USP_Z, include_none)
                  .empty());
}

TEST(SmartsForAtomSubsetEdgeCases, ResetsAromaticBondOutputToDefault) {
  set_global_aromaticity_type(Daylight);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));
  const std::vector<int> include_atom(m.natoms(), 1);

  set_include_bond_aromaticity_in_smiles(1);
  EXPECT_FALSE(
      lillymol::SmartsForAtomSubset(m, kZac, include_atom.data()).empty());
  EXPECT_EQ(include_bond_aromaticity_in_smiles(), 0);
}

}  // namespace
