// Tests for is_actually_chiral

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/is_actually_chiral.h"
namespace {

// A problematic molecule
// This is not completely resolved. is_actually_chiral is causing significant
// failures in iwdescr and needs to be resolved...
// TODO:ianwatson fix bug in is_actually_chiral.
TEST(TestIsActuallyChiral, DISABLED_Test1) {
  Molecule m;
  // This works
  ASSERT_TRUE(m.build_from_smiles("ClC1=CC(=[2C]([1CH](OCCN2CC(C(=O)O)CCC2)[3C]2=C(C)C=C(Cl)C=C2)C=C1)C CHEMBL553153"));
  // This fails
  ASSERT_TRUE(m.build_from_smiles("ClC1=CC=[2C](C(=C1)C)[1CH](OCCN1CCCC(C(=O)O)C1)[3C]1=C(C)C=C(Cl)C=C1 CHEMBL553153"));
  m.compute_aromaticity();
  atom_number_t centre = m.atom_with_isotope(1);
  atom_number_t arom1 = m.atom_with_isotope(2);
  atom_number_t arom2 = m.atom_with_isotope(3);
  ASSERT_NE(centre, kInvalidAtomNumber);
  ASSERT_NE(arom1, kInvalidAtomNumber);
  ASSERT_NE(arom2, kInvalidAtomNumber);

  m.transform_to_non_isotopic_form();

  const int * symm = m.symmetry_classes();
  EXPECT_EQ(symm[arom1], symm[arom2]);

  // No other atom should have the same symmetry class as centre
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == centre) {
      continue;
    }
    EXPECT_NE(symm[centre], symm[i]);
  }

  EXPECT_FALSE(is_actually_chiral(m, centre)) << m.smarts_equivalent_for_atom(centre);
  std::cerr << "atom is " << m.smarts_equivalent_for_atom(centre) << '\n';
}

}  // namespace
