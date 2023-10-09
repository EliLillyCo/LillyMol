// Tests for symmetry

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"

namespace {

using testing::UnorderedElementsAreArray;

class TestSymmetry : public testing::Test {
  protected:
    Molecule _m;
};

TEST_F(TestSymmetry, TestEmptyMolecule) {
  const int * my_null = nullptr;
  EXPECT_EQ(_m.symmetry_classes(), my_null);
}

TEST_F(TestSymmetry, TestOneAtom) {
  ASSERT_TRUE(_m.build_from_smiles("C"));
  EXPECT_EQ(_m.symmetry_class(0), 1);
  EXPECT_EQ(_m.number_symmetry_classes(), 1);
}

TEST_F(TestSymmetry, TestTwoAtoms) {
  ASSERT_TRUE(_m.build_from_smiles("CC"));
  EXPECT_EQ(_m.symmetry_class(0), 1);
  EXPECT_EQ(_m.symmetry_class(1), 1);
  EXPECT_EQ(_m.number_symmetry_classes(), 1);
}

TEST_F(TestSymmetry, TestBenzene) {
  ASSERT_TRUE(_m.build_from_smiles("C1=CC=CC=C1"));
  EXPECT_EQ(_m.number_symmetry_classes(), 1);
  for (int i = 0; i < _m.natoms(); ++i) {
    EXPECT_EQ(_m.symmetry_class(i), 1);
  }
}

TEST_F(TestSymmetry, TestTButyl) {
  ASSERT_TRUE(_m.build_from_smiles("C(C)(C)(C)C"));
  EXPECT_EQ(_m.number_symmetry_classes(), 2);
  EXPECT_EQ(_m.symmetry_class(0), 1);
  for (int i = 1; i < _m.natoms(); ++i) {
    EXPECT_EQ(_m.symmetry_class(i), 2);
  }
}

/*
//       0   6
//      / \ / \
//     1   5   7
//     |   |   |
//     2   4   8
//       \ /\ /
//        3  9
*/

TEST_F(TestSymmetry, TestNaphthalene) {
  ASSERT_TRUE(_m.build_from_smiles("c1cccc2c1cccc2"));
  EXPECT_EQ(_m.number_symmetry_classes(), 3);
  std::vector<int> expected_equilvalent {3, 3, 3, 3, 1, 1, 3, 3, 3, 3};
  std::vector<std::vector<int>> equivalent_atoms {
    {3, 6, 9},
    {2, 7, 8},
    {1, 7, 8},
    {0, 6, 9},
    {5},
    {4},
    {0, 3, 9},
    {1, 2, 8},
    {1, 2, 7},
    {0, 3, 6}};

  for (int i = 0; i < _m.natoms(); ++i) {
    Set_of_Atoms equivalent;
    _m.symmetry_equivalents(i, equivalent);
    EXPECT_EQ(equivalent.number_elements(), expected_equilvalent[i]);
    EXPECT_THAT(equivalent, UnorderedElementsAreArray(equivalent_atoms[i]));
  }
}

TEST_F(TestSymmetry, TestSymmetryClasses) {
  ASSERT_TRUE(_m.build_from_smiles("Fc1c(F)cccc1"));
  const int * symm = _m.symmetry_classes();
  EXPECT_EQ(symm[0], symm[3]);
  EXPECT_EQ(symm[1], symm[2]);
}

}  // namespace
