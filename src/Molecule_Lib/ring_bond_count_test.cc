// Tester for ring bond count

#include <unordered_map>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"

struct MoleculeExpected
{
  IWString smiles;
  int ring_bond_count;
};

namespace {
class TestRingBondCountAllAtomsSame : public testing::TestWithParam<MoleculeExpected> {
};

TEST_P(TestRingBondCountAllAtomsSame, TestAllSame) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    EXPECT_EQ(m.ring_bond_count(i), params.ring_bond_count);
  }
}
INSTANTIATE_TEST_SUITE_P(TestAllSame, TestRingBondCountAllAtomsSame, testing::Values(
  MoleculeExpected{"C", 0},
  MoleculeExpected{"C.C", 0},
  MoleculeExpected{"C.C.C.C.C.C.C.C.C.C.C.C.CCCC", 0},
  MoleculeExpected{"CC", 0},
  MoleculeExpected{"CCC", 0},
  MoleculeExpected{"CC(C)C", 0},
  MoleculeExpected{"CC(C)(C)C", 0},
  MoleculeExpected{"CC(C)(C)CCCCCCC(C)(C)CC", 0},

  MoleculeExpected{"C1CC1", 2},
  MoleculeExpected{"C1CCC1", 2},
  MoleculeExpected{"C1CCCC1", 2},
  MoleculeExpected{"C1CCCCC1", 2},
  MoleculeExpected{"C1CCCCCC1", 2},
  MoleculeExpected{"C1CC1C2CC2", 2},
  MoleculeExpected{"C12C3C4C1C5C2C3C45", 3},
  MoleculeExpected{"C1CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC1", 2},
  MoleculeExpected{"C1CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC1C1CC1", 2},

  MoleculeExpected{"C1CC1.C1CC1.C1CC1.C1CC1.C1CC1.C1CC1.C1CC1.C1CC1.C1CC1.C1CC1.C1CC1", 2}
));

struct MoleculeExpectedSparse
{
  IWString smiles;
  std::unordered_map<int, int> ring_bond_count;
};

class TestRingBondCount : public testing::TestWithParam<MoleculeExpectedSparse> {
};

TEST_P(TestRingBondCount, TestUnchangedMolecule) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  for (const auto& [a, b] : params.ring_bond_count) {
    EXPECT_EQ(m.ring_bond_count(a), b);
  }
}

INSTANTIATE_TEST_SUITE_P(TestSparse, TestRingBondCount, testing::Values(
  MoleculeExpectedSparse{"C", {{0, 0}}},
  MoleculeExpectedSparse{"CC1CC1", {{0, 0}, {1, 2}, {2, 2}, {3,2}}},
  MoleculeExpectedSparse{"CC1CC1CC1CC1", {{0, 0}, {1, 2}, {2, 2}, {3,2}, {4, 0}, {5, 2}, {6, 2}, {7, 2}}},
  MoleculeExpectedSparse{"CC1CC1CCC1CC1", {{0, 0}, {1, 2}, {2, 2}, {3,2}, {4, 0}, {5, 0}, {6, 2}, {7, 2}, {8, 2}}},
  MoleculeExpectedSparse{"C1CC12CC2", {{0, 2}, {1, 2}, {2, 4}, {3,2}, {4, 2}}},
  MoleculeExpectedSparse{"C1CC12CC2C", {{0, 2}, {1, 2}, {2, 4}, {3,2}, {4, 2}, {5, 0}}},
  MoleculeExpectedSparse{"C1C2CC12", {{0, 2}, {1, 3}, {2, 2}, {3, 3}}},
  MoleculeExpectedSparse{"C1C2CCC12", {{0, 2}, {1, 3}, {2, 2}, {3, 2}, {4, 3}}},
  MoleculeExpectedSparse{"C12C3C4C1C243", {{0, 3}, {1, 3}, {2, 3}, {3, 3}, {4, 4}}},
  MoleculeExpectedSparse{"C1C2C(C2)C2C(C2)C1", {{0, 2}, {1, 3}, {2, 3}, {3, 2}, {4, 3}, {5, 3}, {6, 2}, {7, 2}}},
  MoleculeExpectedSparse{"C1C2CC3CC1CC(C2)C3", {{0, 2}, {1, 3}, {2, 2}, {3, 3}, {4, 2}, {5, 3}, {6, 2}, {7, 3}, {8, 2}}},
  MoleculeExpectedSparse{"C1C2C(C2)C2C(C2)C1", {{0, 2}, {1, 3}, {2, 3}, {3, 2}, {4, 3}, {5, 3}, {6, 2}, {7, 2}}},
  MoleculeExpectedSparse{"C12C3C4C5C6C7C8C9C%10C1[U]23456789%10", {{0, 3}, {1, 3}, {2, 3}, {9, 3}, {10, 10}}},
  MoleculeExpectedSparse{"C1C23CC1(C2)C3", {{0, 2}, {1, 4}, {2, 2}, {3, 4}, {4, 2}, {5, 2}}},
  MoleculeExpectedSparse{"C1C23CC1(C2)C3.C", {{0, 2}, {1, 4}, {2, 2}, {3, 4}, {4, 2}, {5, 2}, {6, 0}}}
));

struct AfterAtomRemoval {
  IWString smiles;
  std::unordered_map<int, int> before;
  atom_number_t atom_to_remove;
  std::unordered_map<int, int> after;
};

class TestRemoveAtom : public testing::TestWithParam<AfterAtomRemoval> {
};

TEST_P(TestRemoveAtom, TestAfterAtomRemoval) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  for (const auto& [a, b] : params.before) {
    EXPECT_EQ(m.ring_bond_count(a), b);
  }
  m.remove_atom(params.atom_to_remove);
  for (const auto& [a, b] : params.after) {
    EXPECT_EQ(m.ring_bond_count(a), b);
  }
}

INSTANTIATE_TEST_SUITE_P(TestRemove, TestRemoveAtom, testing::Values(
  AfterAtomRemoval{"CC", {{0, 0}, {1, 0}}, 1, {{0, 0}}},
  AfterAtomRemoval{"C1CC1", {{0, 2}, {1, 2}, {2,2}}, 1, {{0, 0}, {1, 0}}},
  AfterAtomRemoval{"C1CC1C", {{0, 2}, {1, 2}, {2,2}}, 3, {{0, 2}, {1, 2}, {2, 2}}}
));


}  // namespace
