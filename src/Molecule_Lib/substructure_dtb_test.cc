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
  SmilesSmartsNhits{"O1CCCC1", "O-{a{1-3}}CC", 0},  // contains ring
  SmilesSmartsNhits{"S1N=C(CC(N)C(=O)O)C(=O)N1 CHEMBL1094324",
                "[OHD1]-[CD3R0](=O)-[CD3R0](-{[CD2]1;[$(O=c1nsnc1)]1;d3;m5;r5;u1;a7}*)-[ND1H2]", 1},
  SmilesSmartsNhits{"c1ccccc1O", "c1ccccc1-{a1}O", 2},
  SmilesSmartsNhits{"CN", "C-{h1}*", 1},
  SmilesSmartsNhits{"CN", "C-{a1}*", 1},
  SmilesSmartsNhits{"CN", "C-{a1;h1;r0;m0;u0;d0}*", 1}
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

std::ostream&
operator<<(std::ostream& output, const SmilesProtoNhits& proto)  {
  output << proto.smiles << ' ' << proto.proto << ' ' << proto.nhits;

  return output;
}

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
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits) << "failed " << params;
}
INSTANTIATE_TEST_SUITE_P(TestDTBProto, TestDTBProto, testing::Values(
  SmilesProtoNhits{"OC", zero_atoms, 0},
  SmilesProtoNhits{"OC", one_atom, 1},
  SmilesProtoNhits{"OCC", one_atom, 0},
  SmilesProtoNhits{"OC", two_atoms, 0},
  SmilesProtoNhits{"OCC", two_atoms, 1},
  SmilesProtoNhits{"OCCC", two_atoms, 0},
  SmilesProtoNhits{"OCCCC",
    R"pb(
query {
  smarts: "OC-{a3;r0;h0;m0}C"
}
)pb", 1},
  SmilesProtoNhits{"OCCCC",
    R"pb(
query {
  smarts: "OC-{m>0}C"
}
)pb", 0},
  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{m5}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{m{5}}C"
},
)pb", 1},

  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{m{5-}}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{m>4}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{u>0}C"
}
)pb", 0},

  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{r5}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCCn1cccc1",
    R"pb(
query {
  smarts: "OC-{r5;m>4;a7;u0}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCC(=O)N",
    R"pb(
query {
  smarts: "OC-{u1}C"
}
)pb", 0},

  SmilesProtoNhits{"OCCC(=O)N",
    R"pb(
query {
  smarts: "OC-{u2}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCC(=O)N",
    R"pb(
query {
  smarts: "OC-{a4;h2;u2;r0;m0}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCC#C",
    R"pb(
query {
  smarts: "OC-{u2}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCC#C",
    R"pb(
query {
  smarts: "OC-{u2;h0;a3;r0;m0}C"
}
)pb", 1},

  SmilesProtoNhits{"OCCC(=O)N(CCC)CO",
    R"pb(
query {
  smarts: "[CD3T2](=O)N"
  down_the_bond {
    a1: 0
    a2: 2
    heteroatom_count: 1
  }
}
)pb", 0},

  SmilesProtoNhits{"OCCC(=O)N(CCC)CO",
    R"pb(
query {
  smarts: "[CD3T2](=O)N"
  down_the_bond {
    a1: 0
    a2: 2
    heteroatom_count: 1
    match_individual_substituent: true
  }
}
)pb", 1},

  SmilesProtoNhits{"OCCC(=O)N(CCC)CO",
    R"pb(
query {
  smarts: "[CD3T2](=O)N"
  down_the_bond {
    a1: 0
    a2: 2
    heteroatom_count: 1
    match_individual_substituent: true
    no_other_substituents_allowed: true
  }
}
)pb", 0},

  SmilesProtoNhits{"P(=O)(O)(O)CCCC(N)C(O)=O CHEMBL28862",
    R"pb(
query {
  smarts: "[ND1H2]-[CD3x0](-[R0])-[CD3](=O)-[OD1]",
  down_the_bond {
    a1: 1
    a2: 2
    query_matches {
      smarts: "[$(P(=O)-([OH])-[OH])]"
      hits_needed: 0
    }
  }
}
)pb", 0},

  SmilesProtoNhits{"C(N)(C(=O)O)C(O)C CHEMBL30037",
    R"pb(
query {
  smarts: "[ND1H2]-[CD3x0](-[R0])-[CD3](=O)-[OD1]",
  down_the_bond {
    a1: 1
    a2: 2
    query_matches {
      smarts: "[$([OD1]-[CD3])]"
      hits_needed: 1
    }
  }
}
)pb", 1},

  SmilesProtoNhits{"C(=O)(O)[C@@H](N)[C@H](O)[C@H](O)C(=O)O CHEMBL28259",
    R"pb(
query {
  smarts: "[ND1H2]-[CD3x0](-[R0])-[CD3](=O)-[OD1]",
  down_the_bond {
    a1: 1
    a2: 2
    query_matches {
      smarts: "[$([OD1]-[CD3G0])]"
      min_hits_needed: 2
    }
  }
},
)pb", 1},

  SmilesProtoNhits{"C(=O)(O)C(N)CCCCN CHEMBL28328",
    R"pb(
query {
  smarts: "[ND1H2]-[CD3x0](-[R0])-[CD3](=O)-[OD1]",
  down_the_bond {
    a1: 1
    a2: 2
    query_matches {
      smarts: "[CD2]",
      hits_needed: 4
    }
  }
}
)pb", 1},

  SmilesProtoNhits{"C(=O)(O)C(N)CCCCN CHEMBL28328",
    R"pb(
query {
  smarts: "[ND1H2]-[CD3x0](-[R0])-[CD3](=O)-[OD1]",
  down_the_bond {
    a1: 1
    a2: 2
    query_matches {
      smarts: "[CD2]",
      max_hits_needed: 3
    }
  }
}
)pb", 0},

  SmilesProtoNhits{"c1ccccc1CC(C)(C)C",
    R"pb(
query {
  smarts: "c-[CH2][CD4]"
  down_the_bond {
    a1: 1
    a2: 2
    natoms: 4
  }
}
)pb", 1},

  SmilesProtoNhits{"c1ccccc1CC(C)(C)C",
    R"pb(
query {
  smarts: "c-[CH2][CD4]"
  down_the_bond {
    a1: 1
    a2: 2
    match_individual_substituent: true
    natoms: 2
    heteroatom_count: 0
    ring_atom_count: 0
    unsaturation_count: 0
    aromatic_count: 0
  }
}
)pb", 1},

  SmilesProtoNhits{"FC(C)N",
    R"pb(
query {
  smarts: "F[CD3]",
  down_the_bond {
    match_individual_substituent: true
    a1: 0
    a2: 1
    natoms: 2
    heteroatom_count: 1
  }
}
)pb", 1},

  SmilesProtoNhits{"NCCC(F)(F)F",
    R"pb(
query {
  smarts: "[ND1H2]-[CD2]",
  down_the_bond {
    a1: 0
    a2: 1
    query_matches {
      smarts: "[$([CD4](F)(F)F)]",
      hits_needed: 1
    }
  }
}
)pb", 1}

));

}  // namespace
