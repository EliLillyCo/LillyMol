// Tests for substructure ring specifiers

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/substructure.h"

namespace {

struct Data {
  IWString smiles;
  std::string proto;
  int nhits;
};

class TestRingSys: public testing::TestWithParam<Data> {
  protected:
    IWString _smiles;
    Molecule _mol;
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
};

// Repeatedly used as parameter inputs.
const std::string f_env_query = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 1
  ring_specifier {
    base {
      environment: "2a-[$(N(=O)=O),$(C#N)]"
      set_global_id: 8
    }
  }
}
)pb";

const std::string f_no_operator_env_query = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 1
  ring_specifier {
    base {
      environment: "a-N(=O)=O"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specifier {
    base {
      environment: "a-C#N"
      set_global_id: 8
    }
    aromatic: 1
  }
}
)pb";

const std::string f_or_operator_env_query = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 1
  ring_specifier {
    base {
      environment: "a-N(=O)=O"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specifier {
    base {
      environment: "a-C#N"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specification_logexp: SS_OR
}
)pb";

const std::string f_xor_operator_env_query = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 1
  ring_specifier {
    base {
      environment: "a-N(=O)=O"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specifier {
    base {
      environment: "a-C#N"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specification_logexp: SS_XOR
}
)pb";

const std::string f_lpand_operator_env_query = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 1
  ring_specifier {
    base {
      environment: "a-N(=O)=O"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specifier {
    base {
      environment: "a-C#N"
      set_global_id: 8
    }
    aromatic: 1
  }
  ring_specification_logexp: SS_LP_AND
}
)pb";

const std::string f_complex_environment = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 1
  ring_specifier {
    base {
      environment: "1a-[O,N]&&2a-C"
      set_global_id: 8
    }
    aromatic: 1
  }
}
)pb";

const std::string f_arom = R"pb(
query {
  smarts: "F"
  ring_specifier {
    aromatic: 1
  }
}
)pb";

const std::string ring_size = R"pb(
query {
  ring_specifier {
    ring_size: 3
    ring_size: 5
  }
}
)pb";

const std::string fused_0 = R"pb(
query {
  ring_specifier {
    fused: 0
  }
}
)pb";

const std::string min_fused = R"pb(
query {
  ring_specifier {
    min_fused: 1
  }
}
)pb";

const std::string f_ring_size = R"pb(
query {
  smarts: "F-[/IWgid8c]"
  hits_needed: 2
  ring_specifier {
    base {
      set_global_id: 8
    }
    aromatic: 1
    ring_size: 6
  }
}
)pb";

const std::string heteroatom_count = R"pb(
query {
  ring_specifier {
    base {
      heteroatom_count: 1
    }
  }
}
)pb";

const std::string attached_heteroatom_count = R"pb(
query {
  ring_specifier {
    base {
      attached_heteroatom_count: 2
    }
  }
}
)pb";

const std::string ncon = R"pb(
query {
  ring_specifier {
    base {
      ncon: 2
    }
  }
}
)pb";

const std::string all_hits_in_same_fragment = R"pb(
query {
  ring_specifier {
    base {
      hits_needed: 2
      all_hits_in_same_fragment: true
    }
    ring_size: 3
  }
}
)pb";

const std::string atoms_with_pi_electrons = R"pb(
query {
  ring_specifier {
    base {
      atoms_with_pi_electrons: 1
    }
    ring_size: 3
  }
}
)pb";

const std::string fully_saturated_atoms = R"pb(
query {
  ring_specifier {
    base {
      fully_saturated_atoms: 1
    }
  }
}
)pb";

const std::string fully_saturated_atoms_8 = R"pb(
query {
  ring_specifier {
    base {
      fully_saturated_atoms: 8
    }
  }
}
)pb";

const std::string two_ring_sys = R"pb(
query {
  name: "Cl/F -> 2 EWD"
  smarts: "[Cl,F]-[/IWgid8c].[/IWgid8a]!@[/IWgid9a]"
  unique_embeddings_only: true

  ring_system_specifier {
    base {
      set_global_id: 8
      environment: "1a-[F,Cl]"
    }
  }

  ring_system_specifier {
    base {
      set_global_id: 9
      environment: "a-!@a&&0[/IWgid8a]"
    }
    min_aromatic_ring_count: 1
  }
}
)pb";

TEST_P(TestRingSys, TestOperators) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  //std::cerr << "Testing " << params.smiles << ' ' << params.proto << "\n expecting " << params.nhits << '\n';
  EXPECT_EQ(_query.substructure_search(&_mol), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRingSys, TestRingSys, testing::Values(
  Data{"Fc1ccccc1", f_env_query, 0},
  Data{"Fc1c(N(=O)=O)cccc1", f_env_query, 0},
  Data{"Fc1cc(N(=O)=O)ccc1", f_env_query, 0},
  Data{"Fc1ccc(N(=O)=O)cc1", f_env_query, 0},
  Data{"Fc1c(C#N)cccc1", f_env_query, 0},
  Data{"Fc1cc(C#N)ccc1", f_env_query, 0},
  Data{"Fc1ccc(C#N)cc1", f_env_query, 0},
  Data{"Fc1cc(O)c(C#N)cc1", f_env_query, 0},
  Data{"Fc1cc(C#N)c(C#N)cc1", f_env_query, 1},
  Data{"Fc1cc(N(=O)=O)c(C#N)cc1", f_env_query, 1},
  Data{"Fc1cc(N(=O)=O)c(N(=O)=O)cc1", f_env_query, 1},
  Data{"Fc1c(F)c(N(=O)=O)c(N(=O)=O)cc1", f_env_query, 0},
  Data{"Fc1c(C#N)c(N(=O)=O)c(N(=O)=O)cc1", f_env_query, 0},

  Data{"Fc1ccccc1", f_no_operator_env_query, 0},
  Data{"Fc1c(N(=O)=O)cccc1", f_no_operator_env_query, 0},
  Data{"Fc1c(C#N)cccc1", f_no_operator_env_query, 0},
  Data{"Fc1c(C#N)c(C#N)ccc1", f_no_operator_env_query, 0},
  Data{"Fc1c(C#N)c(N(=O)=O)ccc1", f_no_operator_env_query, 1},
  Data{"Fc1c(C#N)c(N(=O)=O)cc(C#N)c1", f_no_operator_env_query, 1},

  Data{"Fc1ccccc1", f_or_operator_env_query, 0},
  Data{"Fc1ccc(C#N)cc1", f_or_operator_env_query, 1},
  Data{"Fc1ccc(N(=O)=O)cc1", f_or_operator_env_query, 1},
  Data{"Fc1c(N(=O)=O)cc(N(=O)=O)cc1", f_or_operator_env_query, 1},
  Data{"Fc1c(N(=O)=O)cc(C#N)cc1", f_or_operator_env_query, 1},

  Data{"Fc1ccccc1", f_xor_operator_env_query, 0},
  Data{"Fc1ccc(C#N)cc1", f_xor_operator_env_query, 1},
  Data{"Fc1ccc(N(=O)=O)cc1", f_xor_operator_env_query, 1},
  Data{"Fc1c(N(=O)=O)cc(N(=O)=O)cc1", f_xor_operator_env_query, 1},
  Data{"Fc1c(C#N)cc(C#N)cc1", f_xor_operator_env_query, 1},
  Data{"Fc1c(N(=O)=O)cc(C#N)cc1", f_xor_operator_env_query, 0},

  Data{"Fc1ccccc1", f_lpand_operator_env_query, 0},
  Data{"Fc1ccc(C#N)cc1", f_lpand_operator_env_query, 0},
  Data{"Fc1ccc(N(=O)=O)cc1", f_lpand_operator_env_query, 0},
  Data{"Fc1c(N(=O)=O)cc(N(=O)=O)cc1", f_lpand_operator_env_query, 0},
  Data{"Fc1c(C#N)cc(N(=O)=O)cc1", f_lpand_operator_env_query, 1},
  Data{"Fc1c(C#N)c(C#N)c(N(=O)=O)cc1", f_lpand_operator_env_query, 1},

  Data{"Fc1ccccc1", f_complex_environment, 0},
  Data{"Fc1c(O)cccc1", f_complex_environment, 0},
  Data{"Fc1c(O)cc(N)cc1", f_complex_environment, 0},
  Data{"Fc1c(O)c(N)c(F)cc1", f_complex_environment, 0},
  Data{"Fc1c(O)cc(C)cc1", f_complex_environment, 0},
  Data{"Fc1c(O)c(C)c(C)cc1", f_complex_environment, 1},
  Data{"Fc1c(N)c(C)c(C)cc1", f_complex_environment, 1},
  Data{"Fc1c(C)c(C)c(C)cc1", f_complex_environment, 0},

  Data{"C1CC1", ring_size, 1},
  Data{"C1CCC1", ring_size, 0},
  Data{"C1CCCC1", ring_size, 1},
  Data{"C1CCCCC1", ring_size, 0},
  Data{"C1CCCCCC1", ring_size, 0},

  Data{"C1CC1", fused_0, 1},
  Data{"C12CC1C2", fused_0, 0},

  Data{"C1CC1", heteroatom_count, 0},
  Data{"C1NC1", heteroatom_count, 1},
  Data{"C1NO1", heteroatom_count, 0},

  Data{"C1CC1", attached_heteroatom_count, 0},
  Data{"C1CC1F", attached_heteroatom_count, 0},
  Data{"FC1CC1F", attached_heteroatom_count, 1},
  Data{"FC1C(F)C1F", attached_heteroatom_count, 0},

  Data{"C1CC1", ncon, 0},
  Data{"C1CC1F", ncon, 0},
  Data{"FC1CC1F", ncon, 1},
  Data{"FC1C(F)C1F", ncon, 0},

  Data{"C1CC1", atoms_with_pi_electrons, 0},
  Data{"C1CC1=O", atoms_with_pi_electrons, 1},
  Data{"C1CCC1=O", atoms_with_pi_electrons, 0},
  Data{"O=C1CC1=O", atoms_with_pi_electrons, 0},

  Data{"O=C1CC1=O", fully_saturated_atoms, 1},
  Data{"O=C1CCC1=O", fully_saturated_atoms, 0},
  Data{"C1=CCCCCCCCC1", fully_saturated_atoms_8, 1},

  Data{"C1CC1", all_hits_in_same_fragment, 0},
  Data{"C1CC1CC1CC1", all_hits_in_same_fragment, 1},
  Data{"C12CC1C2", all_hits_in_same_fragment, 1},
  Data{"C1CC1.C1CC1", all_hits_in_same_fragment, 0},

  Data{"C1CC1", min_fused, 0},
  Data{"C12CC1C2", min_fused, 1},
  // Note that global conditions do NOT report numeric matches. Change sometime?
  Data{"C123CC1C2CC3", min_fused, 1},

  Data{"Fc1ccccc1", f_ring_size, 0},
  Data{"Fc1cc(F)ccc1", f_ring_size, 2},
  Data{"FC1CC(F)CCC1", f_ring_size, 0},
  Data{"FC1CC(F)C1", f_ring_size, 0},

  Data{"OC(=O)c1ccc(cc1F)c2ccccc2[N+](=O)[O-]", two_ring_sys, 1}
));

const std::string ring_includes_carbonyl = R"pb(
query {
  smarts: "[/IWgid1]"
  compress_embeddings: true
  ring_system_specifier {
    base {
      environment: "C"
      set_global_id: 1
      ring_extends_to_carbonyl: true
    }
  }
}
)pb";

const std::string ring_systems_span_spiro2 = R"pb(
query {
  smarts: "[/IWgid1]"
  compress_embeddings: true
  ring_system_specifier {
    rings_in_system: 2
    ring_systems_extend_across_spiro_fusions: true
    base {
      environment: "C"
      set_global_id: 1
    }
  }
}
)pb";

const std::string ring_systems_span_spiro3 = R"pb(
query {
  smarts: "[/IWgid1]"
  compress_embeddings: true
  ring_system_specifier {
    rings_in_system: 3
    ring_systems_extend_across_spiro_fusions: true
    base {
      environment: "C"
      set_global_id: 1
    }
  }
}
)pb";

const std::string ring_systems_substituent_len_2 = R"pb(
query {
  smarts: "[/IWgid1]"
  ring_system_specifier {
    rings_in_system: 2
    base {
      substituent {
        length: 2
        set_global_id: 1
      }
    }
  }
}
)pb";

const std::string ring_spiro_count0 = R"pb(
query {
  ring_specifier {
    spiro_fusion_count: 0
    base {
      set_global_id: 1
    }
  }
  smarts: "[/IWgid1]"
}
)pb";

const std::string ring_spiro_count1 = R"pb(
query {
  ring_specifier {
    spiro_fusion_count: 1
    base {
      set_global_id: 1
    }
  }
  smarts: "[/IWgid1]"
}
)pb";


class TestRingSysCarbonyl: public testing::TestWithParam<Data> {
  protected:
    IWString _smiles;
    Molecule _mol;
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Substructure_Results _sresults;
};

TEST_P(TestRingSysCarbonyl, TestMatches) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  // std::cerr << "Testing " << params.smiles << ' ' << params.proto << "\n expecting " << params.nhits << '\n';
  ASSERT_GT(_query.substructure_search(&_mol, _sresults), 0);
  EXPECT_EQ(_sresults.embedding(0)->size(), params.nhits);
}
// Note that only tests that generate matches can be tested here.
INSTANTIATE_TEST_SUITE_P(TestRingSysCarbonyl, TestRingSysCarbonyl, testing::Values(
  Data{"C1CC1", ring_includes_carbonyl, 3},
  Data{"O=C1CC1", ring_includes_carbonyl, 4},
  Data{"O=C1CC1=N", ring_includes_carbonyl, 5},
  Data{"O=C1C(C)C1=N", ring_includes_carbonyl, 5},
  Data{"C1C2CC12", ring_includes_carbonyl, 4},
  Data{"O=C1C2CC12", ring_includes_carbonyl, 5},
  Data{"O=C1C2C(=O)C12", ring_includes_carbonyl, 6},
  Data{"S=C1C2C(=O)C12", ring_includes_carbonyl, 6},

  //Data{"C1CC1C1CC1", ring_systems_span_spiro2, 0}
  Data{"CC1C2CC12", ring_systems_span_spiro2, 4},  // not spiro fused
  Data{"c1ccccc1C1CC12CC2", ring_systems_span_spiro2, 5},
  Data{"c1ccccc1C1CC12CC23CC3", ring_systems_span_spiro3, 7},
  Data{"c1ccccc1C1CC12C3CCC23", ring_systems_span_spiro3, 7},

  Data{"C1CC1C1CC1", ring_spiro_count0, 1},
  Data{"C12(CC1)CC2", ring_spiro_count1, 1}
));


class TestRingSysOkZero: public testing::TestWithParam<Data> {
  protected:
    IWString _smiles;
    Molecule _mol;
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Substructure_Results _sresults;
};

TEST_P(TestRingSysOkZero, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  //Substructure_Results qq;
  //std::cerr << "Testing " << params.smiles << ' ' << params.proto << "\n expecting " << params.nhits << '\n';
  //std::cerr << "Query result " << _query.substructure_search(&_mol, qq) << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.nhits);
}
// Generally tests that produce no matches should be tested here.
INSTANTIATE_TEST_SUITE_P(TestRingSysOkZero, TestRingSysOkZero, testing::Values(
  Data{"C1CC1CC", ring_systems_substituent_len_2, 0},
  Data{"C1CC1CCC", ring_systems_substituent_len_2, 0},
  Data{"C1CC1CCC", ring_systems_substituent_len_2, 0},
  Data{"C12CC1C2C", ring_systems_substituent_len_2, 0},
  Data{"C12CC1C2CC", ring_systems_substituent_len_2, 2},
  Data{"C12CC1C2CCC", ring_systems_substituent_len_2, 0},
  Data{"C123CC1C2C3CC", ring_systems_substituent_len_2, 0},

  Data{"C1CC1C1CC1", ring_spiro_count1, 0},
  Data{"C12(CC1)CC2", ring_spiro_count0, 0}
));

}  // namespace

