// Tests for some of the functions in rxnfile.cc

#include "rxn_file.h"

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

} // namespace
