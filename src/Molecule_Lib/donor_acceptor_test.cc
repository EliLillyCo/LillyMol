// Tester for Donor_Acceptor

#include <filesystem>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/qry_wstats.h"

namespace {

#ifdef USES_TEXTPROTO_NOT_IMPLEMENTED_YET
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

class TestHbondsProto: public testing::TestWithParam<ProtoMolResult> {
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

TEST_P(TestHbondsProto, TestBuilding) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(_bruns.BuildFromProto(_proto, ""));
  EXPECT_EQ(_bruns.process(_mol), params.numeric_result);
  EXPECT_EQ(_mol.smiles(), params.result);
}
INSTANTIATE_TEST_SUITE_P(TestHbonds, TestHbondsProto, testing::Values(
  ProtoMolResult{default_proto, "CC", 0, "CC"}
));

#endif

// Traditional style file specifications.
struct MolResult {
  IWString smiles;
  int numeric_result;
  IWString result;
};

class TestHbonds: public testing::TestWithParam<MolResult> {
  protected:
    Molecule _mol;
    Donor_Acceptor_Assigner _bruns;

  // protected functions.
    void SetUp();
};

void
TestHbonds::SetUp() {
  const char* test_srcdir = getenv("TEST_SRCDIR");
  IWString queries_dir(test_srcdir);
  queries_dir << "/../donor_acceptor_test.runfiles/donor_acceptor/";

#ifdef LIST_DIRECTORY
  std::string qq(test_srcdir);
  qq += "/../donor_acceptor_test.runfiles/donor_acceptor/";
  for (auto const& dir_entry : std::filesystem::directory_iterator{qq}) {
    std::cerr << dir_entry << '\n';
  }
#endif

  IWString cmd;
  cmd << "a=F:" << queries_dir << "acceptor d=" << queries_dir << "donor.qry";
  if (_bruns.build(cmd)) {
    // std::cerr << "DOnor acceptor initialised " << queries_dir << '\n';
  } else {
    std::cerr << "Could not build donor acceptor from '" << queries_dir << "'\n";
  }

  _bruns.set_apply_isotopic_labels(1);
}

TEST_P(TestHbonds, TestBuilding) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  EXPECT_EQ(_bruns.process(_mol), params.numeric_result);
  EXPECT_EQ(_mol.aromatic_smiles(), params.result) << _mol.aromatic_smiles() <<
        " expected " << params.result;
}
INSTANTIATE_TEST_SUITE_P(TestHbonds, TestHbonds, testing::Values(
  MolResult{"CC", 0, "CC"},
  MolResult{"CCN", 2, "CC[2NH2]"},
  MolResult{"O=CC1CCC1 CHEMBL18475", 1, "[1O]=CC1CCC1"},
  MolResult{"CCC(O)CC CHEMBL47100", 2, "CCC([2OH])CC"},
  MolResult{"N1C(=O)CNC1 CHEMBL30446 CHEMBL30446", 4, "[3NH]1C(=[1O])C[2NH]C1"},
  MolResult{"C(=[1S])=NCC=C CHEMBL233248", 1, "C(=[1S])=NCC=C"}
));

} // namespace
