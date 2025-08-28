// Tests for donor acceptor.
#include <filesystem>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/qry_wstats.h"

namespace {

struct MoleculeResults {
  IWString smiles;
  std::vector<ChargeAndQuery> results;
};

std::ostream&
operator <<(std::ostream& output, const MoleculeResults& mr) {
  output << "MoleculeResults:smiles " << mr.smiles << " with " << mr.results.size() << " results\n";
  for (const ChargeAndQuery& afq : mr.results) {
    output << afq << '\n';
  }

  return output;
}

class TestChargeAssigner: public testing::TestWithParam<MoleculeResults> {
  protected:
    Molecule _mol;
    Charge_Assigner _charge_assigner;
    std::vector<ChargeAndQuery> _results;

  // protected functions.
    void SetUp();
};

void
TestChargeAssigner::SetUp() {
  const char* test_srcdir = getenv("TEST_SRCDIR");
  if (test_srcdir == NULL) {
    std::cerr << "TEST_SRCDIR not set, catastrophe!\n";
    return;
  }

  IWString queries_dir(test_srcdir);
  queries_dir << "/data+/queries/charges/";

// Helps to figure out the directory structure. Should not be this hard
// but the recipes I have found do not work.
//#define LIST_DIRECTORY
#ifdef LIST_DIRECTORY
  std::string qq(test_srcdir);
  std::cerr << "Copied to '" << qq << "'\n";
  qq += "/data+/queries/charges";
  for (auto const& dir_entry : std::filesystem::directory_iterator{qq}) {
    std::cerr << dir_entry << '\n';
  }
#endif

  IWString cmd;
  cmd << "F:" << queries_dir << '/' << "queries";
  if (_charge_assigner.build(cmd)) {
    // std::cerr << "Charge Assigner initialised " << queries_dir << '\n';
  } else {
    std::cerr << "Could not build charge assigner from '" << queries_dir << "'\n";
  }
}
TEST_P(TestChargeAssigner, TestBuilding) {
  ASSERT_TRUE(_charge_assigner.active());

  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));

  _charge_assigner.Process(_mol, _results);

  EXPECT_EQ(_results.size(), params.results.size()) << "failed " << params;
  EXPECT_THAT(_results, testing::ContainerEq(params.results)) << "failed " << params;
}
INSTANTIATE_TEST_SUITE_P(TestChargeAssigner, TestChargeAssigner, testing::Values(
  MoleculeResults {
    "C",
    {}
  },
  MoleculeResults {
    "NCCCC",
    {
      ChargeAndQuery(0, 1, 12)
    }
  },
  MoleculeResults {
    "C(=O)(O)C=O CHEMBL3986754 (1 matches to 'carboxylic_acid')",
    {
      ChargeAndQuery(2, -1, 0)
    }
  },
  MoleculeResults {
    "CS(O)(=O)=O CHEMBL3949263 (1 matches to 'sulfonic_acid')",
    {
      ChargeAndQuery(2, -1, 1)
    }
  },
  MoleculeResults {
    "NC(C)P(O)(O)=O CHEMBL296494 (1 matches to 'phosphonate')",
    {
      // Kind of arbitrary as to whether atom 4 or atom 5 is hit.
      ChargeAndQuery(5, -1, 2)
    }
  },
  MoleculeResults {
    "C1NCC1OC1=NN=NN1 CHEMBL1165160 (1 matches to 'tetrazole')",
    {
      ChargeAndQuery(1, 1, 12),
      ChargeAndQuery(9, -1, 3)
    }
  },
  MoleculeResults {
    "C1=CC=C2C(=O)NS(=O)(=O)NC2=C1 CHEMBL503420 (2 matches to 'N_acylsulfonamide')",
    {
      ChargeAndQuery(6, -1, 4),
    }
  },
  MoleculeResults {
    "S1C(=CC2=CC=CO2)C(=O)NC1=S CHEMBL2007469 (1 matches to 'thiazolidinedione')",
    {
      ChargeAndQuery(10, -1, 5),
    }
  },
  MoleculeResults {
    "O=C1N(NC(=O)C1(CC)CC)C(=O)C1=CC=CC=C1 CHEMBL1966387 (1 matches to 'N_acylpyrazolidinone')",
    {
      ChargeAndQuery(3, -1, 6),
    }
  },
  MoleculeResults {
    "C1=C(Cl)C(=CC(=C1Cl)Cl)O CHEMBL109095 (1 matches to 'phenolate')",
    {
      ChargeAndQuery(9, -1, 7),
    }
  },
  MoleculeResults {
    "C(N)(=N)N1CCCCC1 CHEMBL103102 (1 matches to 'guanidine')",
    {
      ChargeAndQuery(2, 1, 8),
    }
  },
  MoleculeResults {
    "C1(=NCCC1C)N CHEMBL359703 (1 matches to 'amidine')",
    {
      ChargeAndQuery(1, 1, 9),
    }
  },
  MoleculeResults {
    "C1(=C(N)C=CN=C1)C1N(C)CCC1 CHEMBL193763 (1 matches to '4_amino_pyridine')",
    {
      ChargeAndQuery(5, 1, 10),
      ChargeAndQuery(8, 1, 12),
    }
  },
  MoleculeResults {
    "NCC(=O)N1CCCC1 CHEMBL97025 (1 matches to 'amino_terminal_restricted')",
    {
      ChargeAndQuery(0, 1, 11),
    }
  },
  MoleculeResults {
    "C(N)C(N)CS CHEMBL3247435 (2 matches to 'aliphatic_amine_restricted')",
    {
      ChargeAndQuery(1, 1, 12),
    }
  },
  MoleculeResults {
    "N1C(=CN=C1N)CCCC(N)C(O)=O CHEMBL1099170 (1 matches to 'imidazole_basic')",
    {
      ChargeAndQuery(3, 1, 13),
      ChargeAndQuery(10, 1, 12),
      ChargeAndQuery(12, -1, 0),
    }
  },
  // Other interesting molecules
  MoleculeResults {
    "N1(CCN(CC1)CCOC(C1=CC=CC=C1)C1=CC=CC=C1)CC=CC1=CC=CC(=C1)O CHEMBL422499",
    {
      ChargeAndQuery(0, 1, 12),
    }
  }
));

// This class tests the API signature with an array. In this case the query_number
// attribute will not be specified.
class TestChargeAssignerArray: public testing::TestWithParam<MoleculeResults> {
  protected:
    Molecule _mol;
    Charge_Assigner _charge_assigner;
    std::unique_ptr<formal_charge_t[]> _charges;

  // protected functions.
    void SetUp();
};
void
TestChargeAssignerArray::SetUp() {
  const char* test_srcdir = getenv("TEST_SRCDIR");
  if (test_srcdir == NULL) {
    std::cerr << "TEST_SRCDIR not set, catastrophe!\n";
    return;
  }

  IWString queries_dir(test_srcdir);
  queries_dir << "/data+/queries/charges/";

// Helps to figure out the directory structure. Should not be this hard
// but the recipes I have found do not work.
//#define LIST_DIRECTORY
#ifdef LIST_DIRECTORY
  std::string qq(test_srcdir);
  std::cerr << "Copied to '" << qq << "'\n";
  qq += "/data+/queries/charges";
  for (auto const& dir_entry : std::filesystem::directory_iterator{qq}) {
    std::cerr << dir_entry << '\n';
  }
#endif

  IWString cmd;
  cmd << "F:" << queries_dir << '/' << "queries";
  if (_charge_assigner.build(cmd)) {
    // std::cerr << "Charge Assigner initialised " << queries_dir << '\n';
  } else {
    std::cerr << "Could not build charge assigner from '" << queries_dir << "'\n";
  }
}
TEST_P(TestChargeAssignerArray, TestBuilding) {
  ASSERT_TRUE(_charge_assigner.active());

  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));

  const int matoms = _mol.natoms();
  _charges = std::make_unique<formal_charge_t[]>(matoms);

  _charge_assigner.process(_mol, _charges.get());

  for (const ChargeAndQuery& mr : params.results) {
    atom_number_t a = mr.atom;
    EXPECT_EQ(_charges[a], mr.formal_charge);
  }

  for (int i = 0; i < matoms; ++i) {
    if (_charges[i] == 0) {
      continue;
    }

    // Look for this atom in `params`.
    bool found_match = false;
    for (const ChargeAndQuery& mr : params.results) {
      if (i == mr.atom) {
        found_match = true;
        break;
      }
    }
    EXPECT_TRUE(found_match) << _mol.smiles() << " atom " << i << " charged in array";
  }
}
INSTANTIATE_TEST_SUITE_P(TestChargeAssignerArray, TestChargeAssignerArray, testing::Values(
  MoleculeResults {
    "C",
    {}
  },
  MoleculeResults {
    "N1(CCN(CC1)CCOC(C1=CC=CC=C1)C1=CC=CC=C1)CC=CC1=CC=CC(=C1)O CHEMBL422499",
    {
      ChargeAndQuery(0, 1, 12),
    }
  }
));

}  // namespace
