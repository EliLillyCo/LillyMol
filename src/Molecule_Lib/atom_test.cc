// Tester for the atom class

//#include "googlemock/include/gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"

namespace {

// Used with EXPECT_NEAR comparisons.
constexpr float kAbsTol = 0.0001f;

TEST(TestAtom, TestDihedralZero) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, 1.0, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, 1.0, 0.0);
  EXPECT_FLOAT_EQ(m.dihedral_angle(0, 1, 2, 3), 0.0f);
}

TEST(TestAtom, TestDihedralSmallPositive) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, 1.0, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, 1.0, 0.1);
  EXPECT_NEAR(m.dihedral_angle(0, 1, 2, 3), 0.09966857f, kAbsTol);
}

TEST(TestAtom, TestDihedralSmallNegative) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, 1.0, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, 1.0, -0.1);
  EXPECT_NEAR(m.dihedral_angle(0, 1, 2, 3), 0.09966857f, kAbsTol);
}

TEST(TestAtom, TestDihedral90) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, 1.0, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, 0.0, 1.0);
  EXPECT_NEAR(m.dihedral_angle(0, 1, 2, 3), 0.5 * M_PI, kAbsTol);
}

TEST(TestAtom, TestDihedralneg90) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, 1.0, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, 0.0, -1.0);
  EXPECT_NEAR(m.dihedral_angle(0, 1, 2, 3), 0.5 * M_PI, kAbsTol);
}

TEST(TestAtom, TestDihedralneg180) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  m.setxyz(0, -0.5, 1.0, 0.0);
  m.setxyz(1, 0.0, 0.0, 0.0);
  m.setxyz(2, 1.0, 0.0, 0.0);
  m.setxyz(3, 1.5, -1.0, 0.0);
  EXPECT_NEAR(m.dihedral_angle(0, 1, 2, 3), M_PI, kAbsTol);
}

TEST(TestAtom, TestBondsAndAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C(F)(N)(O)C"));
  const Atom& a = m.atom(0);
  int atoms_found = 0;
  for (const auto [bond, atom] : a.BondsAndConnections(0)) {
    EXPECT_NE(atom, 0);
    ++atoms_found;
  }
  EXPECT_EQ(atoms_found, m.ncon(0));
}

}  //  namespace
