// Tester for rxn_to_openrxn

#include "Molecule_Tools/rxn_to_openrxn.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace rxn_to_openrxn {
namespace {
TEST(TestRxnToOpenRxn, Works1) {
  RXN_File rxn;
  const const_IWSubstring buffer = "[Br:1][CH2:2][CH2:3][OH:4].[CH2:5]([S:7](Cl)(=[O:9])=[O:8])[CH3:6].CCOCC>C(N(CC)CC)C>[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]     US03930839";
  ASSERT_TRUE(rxn.build_from_reaction_smiles(buffer, false));

  JobOptions job_options;

  const ord::Reaction ord_proto = BuildOrdReaction(job_options, buffer, rxn);
//EXPECT_THAT(ord_proto, testing::EqualsProto(R"(...)"));
}

}  //  namespace
}  //  namespace rxn_to_openrxn
