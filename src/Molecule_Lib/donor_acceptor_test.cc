// Tester for Donor_Acceptor

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/qry_wstats.h"

namespace {
struct ProtoMolResult {
  // BrunsDonorAcceptor proto
  // Note that shell variable expansion is performed.
  std::string proto;
  // Smiles of the input molecule.
  IWString smiles;
  // Numeric return code from process()
  int numeric_result;
  // The smiles of the resulting molecule.
  IWString result;
};

class TestHbonds: public testing::TestWithParam<ProtoMolResult> {
  protected:
    BrunsDonorAcceptor::BrunsDonorAcceptor _proto;
    Molecule _mol;
    Donor_Acceptor_Assigner _bruns;
};

const std::string default_proto = R"pb(
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/carbonyl.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/cyano.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/imine.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/aminunch.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/hydroxyl.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/ether.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/hydroxam.textproto"
  acceptor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/qminus.textproto"

  donor: "${C3TK_DATA_PERSISTENT}/queries/hbonds/donor.textproto"
)pb";

TEST_P(TestHbonds, TestBuilding) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(_bruns.BuildFromProto(_proto, ""));
  EXPECT_EQ(_bruns.process(_mol), params.numeric_result);
  EXPECT_EQ(_mol.smiles(), params.result);
}
INSTANTIATE_TEST_SUITE_P(TestHbonds, TestHbonds, testing::Values(
  ProtoMolResult{default_proto, "CC", 0, "CC"}
));

} // namespace
