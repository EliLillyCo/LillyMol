// tests for DownTheBond functionality

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "substructure.h"

namespace {

struct SmilesSmartsNhits {
  IWString smiles;
  IWString smarts;
  int nhits;
};

class TestDTB : public testing::TestWithParam<SmilesSmartsNhits> {
  protected:
    Substructure_Query _query;
    Molecule _m;
};

TEST_P(TestDTB, Test1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  //std::cerr << "TestingH '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestDTB, TestDTB, testing::Values(
  SmilesSmartsNhits{"OC", "OC", 1},
  SmilesSmartsNhits{"OC", "O-{a0}C", 0},
  SmilesSmartsNhits{"OCC", "O-{a1}C", 0},
  SmilesSmartsNhits{"OCC", "O-{a2}CC", 1},
  SmilesSmartsNhits{"OCC", "O-{a3}CC", 0},
  SmilesSmartsNhits{"OCC", "O-{a>1}CC", 1},
  SmilesSmartsNhits{"OCC", "O-{a<3}CC", 1},
  SmilesSmartsNhits{"OCC(C)C", "O-{a{1-3}}C", 0},
  SmilesSmartsNhits{"OCC(C)C", "O-{a4}C", 1},
  SmilesSmartsNhits{"OCC(C)CC", "O-{a{2-7}}C", 1},
  SmilesSmartsNhits{"O1CCCC1", "O-{a{1-3}}CC", 0}  // contains ring
));

const std::string zero_atoms = R"pb(
query {
  smarts: "O-{a0}C"
  down_the_bond {
    a1: 0
    a2: 1
    natoms: 2
  }
}
)pb";

const std::string one_atom = R"pb(
query {
  smarts: "O-{a1}C"
  down_the_bond {
    a1: 0
    a2: 1
    natoms: 1
  }
}
)pb";

const std::string two_atoms = R"pb(
query {
  smarts: "O-{a2}C"
  down_the_bond {
    a1: 0
    a2: 1
    natoms: 2
  }
}
)pb";

struct SmilesProtoNhits {
  IWString smiles;
  std::string proto;
  int nhits;
};

class TestDTBProto : public testing::TestWithParam<SmilesProtoNhits> {
  protected:
    Substructure_Query _query;
    Molecule _m;
};

TEST_P(TestDTBProto, Test1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  SubstructureSearch::SubstructureQuery proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &proto));
  ASSERT_TRUE(_query.ConstructFromProto(proto));

  // std::cerr << "TestingH '" << params.smiles << "' smarts '" << params.proto << " xpt " << params.nhits << '\n';
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestDTBProto, TestDTBProto, testing::Values(
  SmilesProtoNhits{"OC", zero_atoms, 0},
  SmilesProtoNhits{"OC", one_atom, 1},
  SmilesProtoNhits{"OCC", one_atom, 0},
  SmilesProtoNhits{"OC", two_atoms, 0},
  SmilesProtoNhits{"OCC", two_atoms, 1},
  SmilesProtoNhits{"OCCC", two_atoms, 0}
));

}  // namespace
