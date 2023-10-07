// Tester for Element_Hits_Needed
#include <iostream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/substructure.h"

namespace {

using std::cerr;
using std::endl;

const std::string q1 = R"pb(
query{
  smarts: "FC(F)(Cl)CCO"
  element_hits_needed {
    atomic_number: 6
    hits_needed: 0
  }
}
)pb";

struct Data {
  std::string proto;
  IWString smiles;
  int nhits;
};

class TestRingSys: public testing::TestWithParam<Data> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Molecule _mol;
    Substructure_Query _query;
};

TEST_P(TestRingSys, TestOperators) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  //cerr << "Testing " << params.smiles << ' ' << params.proto << "\n expecting " << params.nhits << '\n';
  EXPECT_EQ(_query.substructure_search(&_mol), params.nhits) << params.smiles << " query " << params.proto << " expect " << params.nhits;
}
INSTANTIATE_TEST_SUITE_P(TestRingSys, TestRingSys, testing::Values(
  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: 6
    hits_needed: 0
  }
}
)pb", "FC(F)(Cl)CCO", 0},

  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: 6
    hits_needed: 1
  }
}
)pb", "FC(F)(Cl)CCO", 0},

  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: 6
    hits_needed: 2
  }
}
)pb", "FC(F)(Cl)CCO", 1},

  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: 6
    hits_needed: 3
  }
}
)pb", "FC(F)(Cl)CCO", 0},

  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: 6
    min_hits_needed: 2
  }
}
)pb", "FC(F)(Cl)CCO", 1},

  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: 6
    max_hits_needed: 2
  }
}
)pb", "FC(F)(Cl)CCO", 1},

  Data{
    R"pb(
query{
  smarts: "[F,Cl]C([F,Cl])([F,Cl])C"
  element_hits_needed {
    atomic_number: [9, 17]
    max_hits_needed: 2
    multiple_values_operator: SS_OR
  }
}
)pb", "FC(F)(Cl)CCO", 0},

  Data{
    R"pb(
query{
  smarts: "[F,Cl]C([F,Cl])([F,Cl])C"
  element_hits_needed {
    atomic_number: [9, 17]
    multiple_values_operator: SS_OR
    hits_needed: 3
  }
}
)pb", "FC(F)(Cl)CCO", 6},

  Data{
    R"pb(
query{
  smarts: "[F,Cl]C([F,Cl])([F,Cl])C"
  element_hits_needed {
    atomic_number: [6, 9]
    hits_needed: 4
  }
}
)pb", "FC(F)(Cl)CCO", 6},

  Data{
    R"pb(
query{
  smarts: "[F,Cl]C([F,Cl])([F,Cl])C"
  element_hits_needed {
    atomic_number: [6, 9]
    multiple_values_operator: SS_XOR
  }
}
)pb", "FC(F)(Cl)CCO", 0},

  Data{
    R"pb(
query{
  smarts: "CCO"
  element_hits_needed {
    atomic_number: [6, 7]
    multiple_values_operator: SS_XOR
  }
}
)pb", "CCO", 1},

// also test global condtion elements_needed
  Data{
    R"pb(
query{
  smarts: "C"
  required_molecular_properties {
    elements_needed {
      atomic_number: [6, 9]
      min_hits_needed: 4
    }
  }
}
)pb", "FC(F)(Cl)CCO", 3},

  Data{
    R"pb(
query{
  smarts: "C"
  required_molecular_properties {
    elements_needed {
      atomic_number: [6, 9]
      multiple_values_operator: SS_OR
    }
  }
}
)pb", "CN", 1}

));

}  // namespace
