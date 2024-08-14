// Tester for the separated_atoms query

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "google/protobuf/text_format.h"

#include "aromatic.h"
#include "substructure.h"

namespace {

class TestSeparatedAtoms : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    Substructure_Query _query;

    Molecule _m;

    SubstructureSearch::SubstructureQuery _proto;
};

void
TestSeparatedAtoms::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestSeparatedAtoms, invalid_atoms) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 0
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_FALSE(_query.ConstructFromProto(_proto));
}

TEST_F(TestSeparatedAtoms, invalid_min_max) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 0
        min_bonds_between: 2
        max_bonds_between: 1
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_FALSE(_query.ConstructFromProto(_proto));
}

TEST_F(TestSeparatedAtoms, TooFarApart) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        bonds_between: 1
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("COC"));
  EXPECT_FALSE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, Tooclose) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        min_bonds_between: 3
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("COC"));
  EXPECT_FALSE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, TwoBondsOK) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        bonds_between: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("COC"));
  EXPECT_TRUE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, TwoBondsMinOk) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        min_bonds_between: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("COC"));
  EXPECT_TRUE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, TwoBondsMaxOk) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        max_bonds_between: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("COC"));
  EXPECT_TRUE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, TestRing) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 7
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        max_bonds_between: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("COON"));
  EXPECT_FALSE(_query.substructure_search(&_m));
  ASSERT_TRUE(_m.build_from_smiles("C1OON1"));
  EXPECT_TRUE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, TestMultiple) {
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 16
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 7
        }
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number : 8
        }
      }
      separated_atoms {
        a1: 0
        a2: 1
        max_bonds_between: 2
      }
      separated_atoms {
        a1: 0
        a2: 2
        max_bonds_between: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("S(CO)CN"));
  EXPECT_TRUE(_query.substructure_search(&_m));
  ASSERT_TRUE(_m.build_from_smiles("S(CC)CN"));
  EXPECT_FALSE(_query.substructure_search(&_m));
}

TEST_F(TestSeparatedAtoms, TestRespectInitialNumbering) {
  _string_proto = R"(query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 4
        atom_properties {
          atomic_number : 16
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 7
        }
      }
      query_atom {
        id: 8
        atom_properties {
          atomic_number : 8
        }
      }
      separated_atoms {
        a1: 4
        a2: 1
        max_bonds_between: 2
      }
      separated_atoms {
        a1: 4
        a2: 8
        max_bonds_between: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_query.ConstructFromProto(_proto));

  ASSERT_TRUE(_m.build_from_smiles("S(CO)CN"));
  EXPECT_TRUE(_query.substructure_search(&_m));
  ASSERT_TRUE(_m.build_from_smiles("S(CC)CN"));
  EXPECT_FALSE(_query.substructure_search(&_m));
}

struct ProtoSmilesResult {
  std::string proto;
  IWString smiles;
  int expected;
};

std::ostream&
operator<< (std::ostream& output, const ProtoSmilesResult& psr) {
  output << psr.proto << ' ' << psr.smiles << ' ' << psr.expected;
  return output;
}

class TestSepma: public testing::TestWithParam<ProtoSmilesResult> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Substructure_Results _sresults;
    Molecule _mol;
};
TEST_P(TestSepma, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  // cerr << "Testing " << params.smiles << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected) << params;
}
INSTANTIATE_TEST_SUITE_P(TestSepma, TestSepma, testing::Values(
  ProtoSmilesResult{
    R"pb(
query {
  smarts: "O.O"
  unique_embeddings_only: true
  separated_atoms {
    a1: 0
    a2: 1
    rotbond: 3
  }
}
)pb", "OCCCCO", 1},

  ProtoSmilesResult{
    R"pb(
query {
  smarts: "O.O"
  separated_atoms {
    a1: 0
    a2: 1
    rotbond: 3
  }
}
)pb", "OCCCCO", 2},

  ProtoSmilesResult{
    R"pb(
query {
  smarts: "O.O"
  separated_atoms {
    a1: 0
    a2: 1
    bonds_between: 5
  }
}
)pb", "OC1CCC(O)CC1", 2},

  ProtoSmilesResult{
    R"pb(
query {
  smarts: "O.O"
  separated_atoms {
    a1: 0
    a2: 1
    bonds_between: 5
    rotbond: 0
  }
}
)pb", "OC1CCC(O)CC1", 2},

  ProtoSmilesResult{
    R"pb(
query {
  smarts: "O.O"
  separated_atoms {
    a1: 0
    a2: 1
    min_rotbond: 2
    max_rotbond: 2
  }
}
)pb", "OCC(=O)NCO", 2}
));

}  // namespace
