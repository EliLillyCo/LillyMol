// Tests for matched atoms match.

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "aromatic.h"
#include "substructure.h"

namespace {

struct ProtoSmilesResult {
  std::string proto_string;
  IWString smiles;
  int expected;
};

class TestMam : public testing::TestWithParam<ProtoSmilesResult> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _m;
};

std::string nitrogen = R"pb(
  query {
    smarts: "N"
    matched_atom_must_be {
      atom: 0
      smarts: "N"
    }
  }
)pb";

std::string nD2 = R"pb(
  query {
    smarts: "N"
    matched_atom_must_be {
      atom: 0
      smarts: "[ND2]"
    }
  }
)pb";

std::string nD23 = R"pb(
  query {
    smarts: "N"
    matched_atom_must_be {
      atom: 0
      smarts: "[ND2]"
      smarts: "[ND3]"
    }
  }
)pb";

std::string not_amide = R"pb(
  query {
    smarts: "N"
    matched_atom_must_be {
      atom: 0
      smarts: "!NC=O"
    }
  }
)pb";

std::string is_amide = R"pb(
  query {
    smarts: "C-N"
    matched_atom_must_be {
      atom: 1
      smarts: "NC=O"
    }
  }
)pb";

std::string two_ethers = R"pb(
  query {
    smarts: "O-c:c-O"
    matched_atom_must_be {
      atom: [0, 3]
      smarts: "[OD2](-c)-C"
    }
  }
)pb";

std::string primary_amine = R"pb(
  query {
    smarts: "N"
    matched_atom_must_be {
      atom: 0
      smarts: "![ND0]"
      smarts: "![NG>0]"
      smarts: "!N-C=O"
      smarts: "!N-S=O"
      smarts: "!N-C=S"
      smarts: "!N-C=S"
      smarts: "![ND2]"
      smarts: "!N-a"
      smarts: "!N-[!#6]"
      smarts: "!N...{<3}N"
    }
  }
)pb";


std::string multiple_conditions = R"pb(
  query {
    smarts: "[!#6]c1ccc([!#6])cc1"
    matched_atom_must_be {
      atom: 0
      smarts: "[F,Cl,Br,I]"
      smarts: "!S"
    }
    matched_atom_must_be {
      atom: 5
      smarts: "[O,N]"
      smarts: "!S"
    }
  }
)pb";

TEST_P(TestMam, Test1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto_string, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  //cerr << "TestingH '" << params.smiles << "' smarts '" << params.proto_string << " xpt " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(&_m), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestMam, TestMam, testing::Values(
  ProtoSmilesResult{nitrogen, "C", 0},
  ProtoSmilesResult{nitrogen, "N", 1},
  ProtoSmilesResult{nD2, "N", 0},
  ProtoSmilesResult{nD2, "NC", 0},
  ProtoSmilesResult{nD2, "CNC", 1},
  ProtoSmilesResult{nD2, "CN(C)C", 0},
  ProtoSmilesResult{nD23, "N", 0},
  ProtoSmilesResult{nD23, "CNC", 1},
  ProtoSmilesResult{nD23, "CN(C)C", 1},
  ProtoSmilesResult{nD23, "C[N+](C)(C)C", 0},
  ProtoSmilesResult{not_amide, "N", 1},
  ProtoSmilesResult{not_amide, "NC", 1},
  ProtoSmilesResult{not_amide, "NC", 1},
  ProtoSmilesResult{not_amide, "NC", 1},
  ProtoSmilesResult{not_amide, "CN(C)C", 1},
  ProtoSmilesResult{not_amide, "CNC=O", 0},
  ProtoSmilesResult{not_amide, "CNC(=O)C", 0},
  ProtoSmilesResult{is_amide, "CNC(=O)C", 2},
  ProtoSmilesResult{is_amide, "CNC(-O)C", 0},

  ProtoSmilesResult{two_ethers, "CNC(-O)C", 0},
  ProtoSmilesResult{two_ethers, "c1ccccc1", 0},
  ProtoSmilesResult{two_ethers, "c1(O)ccccc1",  0},
  ProtoSmilesResult{two_ethers, "c1(O)cc(O)ccc1", 0},
  ProtoSmilesResult{two_ethers, "c1(O)c(O)cccc1", 0},
  ProtoSmilesResult{two_ethers, "c1(OC)c(OC)cccc1", 2},
  ProtoSmilesResult{two_ethers, "c1(O)c(OC)cccc1", 0},

  ProtoSmilesResult{primary_amine, "C", 0},
  ProtoSmilesResult{primary_amine, "N", 0},
  ProtoSmilesResult{primary_amine, "NC", 1},
  ProtoSmilesResult{primary_amine, "CNC", 0},
  ProtoSmilesResult{primary_amine, "NC(O)C", 1},
  ProtoSmilesResult{primary_amine, "NC(=O)C", 0},
  ProtoSmilesResult{primary_amine, "Nc1ccccc1", 0},
  ProtoSmilesResult{primary_amine, "NC(=S)C", 0},
  ProtoSmilesResult{primary_amine, "NS(=O)C", 0},
  ProtoSmilesResult{primary_amine, "NCN", 0},
  ProtoSmilesResult{primary_amine, "NCCN", 0},
  ProtoSmilesResult{primary_amine, "NCCCN", 2},

  ProtoSmilesResult{multiple_conditions, "c1ccccc1", 0},
  ProtoSmilesResult{multiple_conditions, "Fc1ccccc1", 0},
  ProtoSmilesResult{multiple_conditions, "Fc1ccc(F)cc1", 0},
  ProtoSmilesResult{multiple_conditions, "Nc1ccc(N)cc1", 0},
  ProtoSmilesResult{multiple_conditions, "Nc1ccc(F)cc1", 2}
));

}  // namespace
