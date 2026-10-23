// Tests for some of the functions in rxnfile.cc

#include "rxn_file.h"
#include "element.h"

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

class TestRxnFile : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    RXN_File _rxn;
};

void
TestRxnFile::SetUp()
{
  _rxn.set_do_automatic_atom_mapping(0);

  return;
}

TEST_F(TestRxnFile, TestremoveNonParticipatingFragmentsMapped)
{
  const_IWSubstring buffer = "[C:1][C:2](=[O:3])[O:4]+[N:5][C:6]>>[C:1][C:2](=[O:3])[N:5][C:6].[C:8][C:9]";

  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 1));
  EXPECT_EQ(_rxn.number_reagents(), 2);
  EXPECT_EQ(_rxn.number_products(), 1);
  const int removed = _rxn.remove_non_participating_fragments();
  EXPECT_EQ(removed, 1);
  EXPECT_EQ(_rxn.product(0).natoms(), 5);
}

TEST_F(TestRxnFile, TestremoveNonParticipatingFragmentsUnMapped)
{
  const_IWSubstring buffer = "[C:1][C:2](=[O:3])[O:4]+[N:5][C:6]>>CC.[C:1][C:2](=[O:3])[N:5][C:6]";

  EXPECT_TRUE(_rxn.build_from_reaction_smiles(buffer, 1));
  EXPECT_EQ(_rxn.number_reagents(), 2);
  EXPECT_EQ(_rxn.number_products(), 1);
  const int removed = _rxn.remove_non_participating_fragments();
  EXPECT_EQ(removed, 1);
  EXPECT_EQ(_rxn.product(0).natoms(), 5);
}

// MDL/reaction add wrappers retain success/failure returns, whereas Molecule::add
// returns a zero-based index. The first atom must not short-circuit metadata setup.
TEST(MdlMoleculeAdd, FirstAtomReportsSuccess) {
  MDL_Molecule molecule;
  ASSERT_TRUE(molecule.add(get_element_from_atomic_number(6)));
  EXPECT_EQ(molecule.natoms(), 1);
  EXPECT_EQ(molecule.atomic_number(0), 6);
  ASSERT_TRUE(molecule.add(get_element_from_atomic_number(8)));
  EXPECT_EQ(molecule.natoms(), 2);
}

TEST(MdlMoleculeAdd, ExistingMetadataGrowsWithAtoms) {
  Molecule source;
  ASSERT_TRUE(source.build_from_smiles("C"));
  MDL_Molecule molecule(source);
  ASSERT_EQ(molecule.number_atoms(), 1);
  ASSERT_TRUE(molecule.add(get_element_from_atomic_number(8)));
  ASSERT_EQ(molecule.number_atoms(), molecule.natoms());
  ASSERT_NE(molecule.mdl_atom(1), nullptr);
}

TEST(ReactionMoleculeAdd, FirstAtomAllocatesMapsAndLaterAtomsPreserveThem) {
  ISIS_RXN_FILE_Molecule molecule;
  ASSERT_TRUE(molecule.add(get_element_from_atomic_number(12)));
  ASSERT_EQ(molecule.natoms(), 1);
  ASSERT_NE(molecule.atom_map(), nullptr);
  EXPECT_EQ(molecule.atom_map(0), -1);
  molecule.set_atom_map(0, 26);
  ASSERT_TRUE(molecule.add(get_element_from_atomic_number(35)));
  ASSERT_EQ(molecule.natoms(), 2);
  EXPECT_EQ(molecule.atom_map(0), 26);
  EXPECT_EQ(molecule.atom_map(1), -1);
  molecule.set_atom_map(1, 27);
  ASSERT_TRUE(molecule.add_bond(0, 1, SINGLE_BOND, SINGLE_BOND, 0));
  EXPECT_EQ(molecule.nedges(), 1);
  EXPECT_EQ(molecule.atom_map(1), 27);
}

TEST_F(TestRxnFile, RepairsMappedProductWithBondedOrphanAtoms) {
  const_IWSubstring reaction = "[CH4:1]>>[CH3:1][Mg:2][Br:3]";
  ASSERT_TRUE(_rxn.build_from_reaction_smiles(reaction, 1));
  _rxn.set_auto_fix_orphans(1);
  ASSERT_TRUE(_rxn.check_for_widows_and_orphans());
  EXPECT_EQ(_rxn.contains_orphan_atoms(), 2);
}

} // namespace
