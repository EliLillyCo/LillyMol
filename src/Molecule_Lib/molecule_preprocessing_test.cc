#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"

namespace {

using molecule_processing::MoleculePreprocessing;

Molecule
MolFromSmiles(const char* smiles) {
  Molecule result;
  EXPECT_TRUE(result.build_from_smiles(smiles)) << smiles;
  return result;
}

TEST(MoleculePreprocessing, InitiallyInactive) {
  MoleculePreprocessing preprocessing;

  EXPECT_FALSE(preprocessing.active());
}

TEST(MoleculePreprocessing, ReduceToLargestFragment) {
  MoleculePreprocessing preprocessing;
  preprocessing.set_reduce_to_largest_fragment(true);

  EXPECT_TRUE(preprocessing.active());

  Molecule m = MolFromSmiles("C.CC");
  ASSERT_EQ(m.number_fragments(), 2);

  EXPECT_GT(preprocessing.Process(m), 0);

  EXPECT_EQ(m.number_fragments(), 1);
  EXPECT_EQ(m.natoms(), 2);
}

TEST(MoleculePreprocessing, RemoveChirality) {
  MoleculePreprocessing preprocessing;
  preprocessing.set_remove_chirality(true);

  EXPECT_TRUE(preprocessing.active());

  Molecule m = MolFromSmiles("C[C@H](N)F");
  ASSERT_GT(m.chiral_centres(), 0);

  EXPECT_GT(preprocessing.Process(m), 0);

  EXPECT_EQ(m.chiral_centres(), 0);
}

TEST(MoleculePreprocessing, RemoveCisTransBonds) {
  MoleculePreprocessing preprocessing;
  preprocessing.set_remove_cis_trans_bonds(true);

  EXPECT_TRUE(preprocessing.active());

  Molecule m = MolFromSmiles("C\\C=C\\C");

  const IWString& before = m.smiles();
  EXPECT_TRUE(before.contains('\\') || before.contains('/')) << before;

  EXPECT_GT(preprocessing.Process(m), 0);

  const IWString& after = m.smiles();
  EXPECT_FALSE(after.contains('\\')) << after;
  EXPECT_FALSE(after.contains('/')) << after;
}

TEST(MoleculePreprocessing, RemoveIsotopes) {
  MoleculePreprocessing preprocessing;
  preprocessing.set_remove_isotopes(true);

  EXPECT_TRUE(preprocessing.active());

  Molecule m = MolFromSmiles("[1CH4]");
  ASSERT_EQ(m.isotope(0), static_cast<isotope_t>(1));

  EXPECT_GT(preprocessing.Process(m), 0);

  EXPECT_EQ(m.isotope(0), static_cast<isotope_t>(0));
}

TEST(MoleculePreprocessing, MultipleOperations) {
  MoleculePreprocessing preprocessing;
  preprocessing.set_reduce_to_largest_fragment(true);
  preprocessing.set_remove_chirality(true);
  preprocessing.set_remove_isotopes(true);

  EXPECT_TRUE(preprocessing.active());

  Molecule m = MolFromSmiles("[1CH3].[2CH]([C@H](N)F)C");
  ASSERT_EQ(m.number_fragments(), 2);

  EXPECT_GT(preprocessing.Process(m), 0);

  EXPECT_EQ(m.number_fragments(), 1);
  EXPECT_EQ(m.chiral_centres(), 0);

  for (int i = 0; i < m.natoms(); ++i) {
    EXPECT_EQ(m.isotope(i), static_cast<isotope_t>(0)) << "atom " << i;
  }
}

}  // namespace
