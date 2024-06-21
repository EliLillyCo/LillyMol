// Tests for fragment_molecule

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/fragment_molecule.h"

namespace {


TEST(TestFragmentMolecule, TestNoBreakCC) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CC"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 0);
  EXPECT_EQ(bonds.size(), 0);
}

TEST(TestFragmentMolecule, TestBreakCC) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CCC(C)C"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 0);
  fragmenter.set_break_carbon_carbon_bonds(1);
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 4);
  EXPECT_EQ(bonds.size(), 4);
}

TEST(TestFragmentMolecule, TestCAmide) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CNC(=O)C"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 2);
  EXPECT_EQ(bonds.size(), 2);
}

TEST(TestFragmentMolecule, TestSAmide1) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CNS(=O)C"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 2);
  EXPECT_EQ(bonds.size(), 2);
}

}  // namespace
