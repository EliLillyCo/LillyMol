// Unit tests for toggle_kekule_form

#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwstring/iwstring.h"
#include "molecule.h"
#include "set_of_atoms.h"
#include "toggle_kekule_form.h"
#include "Molecule_Lib/toggle_kekule_form.pb.h"

namespace {
struct ToggleTest {
  std::string string_proto;
  IWString smiles;
  Set_of_Atoms embedding;
  int change_expected;
  // Some tests start with invalid molecules.
  int ok_check_usmi;
};

class TestToggleKekuleFormP: public testing::TestWithParam<ToggleTest> {
  protected:
    ToggleKekuleForm::ToggleKekuleForm _proto;

    Toggle_Kekule_Form _toggle_kekule_form;

    Molecule _mol;
};

static constexpr int kOkCheckUsmi = 1;
static constexpr int kNoCheckUsmi = 1;

static constexpr int kUnChanged = 0;
static constexpr int kChanged = 1;

// After toggling has been done, `found` is the bond found in the
// molecule.
// `directed` is what was requested.
// Return true if they match.
bool
BondIsCorrect(const Bond& found, const ToggleKekuleForm::KekuleBond& directed) {
  if (found.is_single_bond() && directed.btype() == SubstructureSearch::BondType::SS_SINGLE_BOND) {
    return true;
  }
  if (found.is_double_bond() && directed.btype() == SubstructureSearch::BondType::SS_DOUBLE_BOND) {
    return true;
  }

  return false;
}

TEST_P(TestToggleKekuleFormP, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.string_proto, &_proto));
  ASSERT_TRUE(_toggle_kekule_form.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));

  Molecule m(_mol);
  int changed = 0;
  _toggle_kekule_form.process(m, params.embedding, changed);
  EXPECT_EQ(changed, params.change_expected) << "change mismatch " << changed << ' '
       << _mol.smiles() << ' ' << m.smiles();

  // Nothing happened and nothing was expected. Good.
  if (! changed && ! params.change_expected) {
    return;
  }

  if (params.ok_check_usmi) {
    EXPECT_EQ(_mol.unique_smiles(), m.unique_smiles()) << "unique smiles mismatch " <<
         _mol.unique_smiles() << " changed " << m.unique_smiles();
    EXPECT_EQ(_mol.aromatic_atom_count(), m.aromatic_atom_count());
  }

  for (const auto& bond : _proto.bond()) {
    int i1 = bond.a1();
    int i2 = bond.a2();
    ASSERT_TRUE(params.embedding.ok_index(i1));
    ASSERT_TRUE(params.embedding.ok_index(i2));
    atom_number_t a1 = params.embedding[i1];
    atom_number_t a2 = params.embedding[i2];
    const Bond* b = m.bond_between_atoms(a1, a2);
    ASSERT_NE(b, nullptr);
    EXPECT_TRUE(BondIsCorrect(*b, bond)) << "bad bond " << *b << " exected " <<
                bond.ShortDebugString();
  }
}
INSTANTIATE_TEST_SUITE_P(TestToggleKekuleFormP, TestToggleKekuleFormP, testing::Values(
  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
)pb",
    "C1=CC=CC=C1",
    {0, 1},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
)pb",
    "C1C=CC=CC=1",
    {0, 1},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
)pb",
    "C12=CC=NC=C2CN1 CHEMBL286232",
    {0, 1},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
)pb",
    "C12=CC=NC=C2CN1 CHEMBL286232",
    {0, 1},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
)pb",
    "C12=CN=NC=C1C=CC=C2 CHEMBL1650268",
    {0, 1},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
)pb",
    "C12=CN=NC=C1C=CC=C2 CHEMBL1650268",
    {0, 1},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
)pb",
    "N=C1N=CN=C2N(C)NC=C12 CHEMBL4572446",
    {1, 2},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
    display_error_messages: false
)pb",
    "N=C1N=CN=C2N(C)NC=C12 CHEMBL4572446",
    {1, 2},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
    display_error_messages: false
)pb",
    "CC1=CC2=CN=CC=C2NC1=O CHEMBL342320",
    {1, 2},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
    display_error_messages: false
)pb",
    "CC1=CC2=CN=CC=C2NC1=O CHEMBL342320",
    {3, 4},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
    display_error_messages: false
)pb",
    "CC1=CC2=CN=CC=C2NC1=O CHEMBL342320",
    {3, 4},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
    display_error_messages: false
)pb",
    "NC1=C2C=NNC2=NC=N1",
    {1, 2},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
    display_error_messages: false
)pb",
    "NC1=C2C=NNC2=NC=N1",
    {1, 2},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
    bond {
      a1: 2
      a2: 3
      btype: SS_DOUBLE_BOND
    }
    display_error_messages: false
)pb",
    "NC1=NC=NC2C3=CC=CN=C3SC=21 CHEMBL457313",
    {1, 2, 6, 7},
    kUnChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
    bond {
      a1: 2
      a2: 3
      btype: SS_DOUBLE_BOND
    }
    display_error_messages: false
)pb",
    "NC1=NC=NC2C3=CC=CN=C3SC=21 CHEMBL457313",
    {1, 2, 6, 7},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
    bond {
      a1: 2
      a2: 3
      btype: SS_SINGLE_BOND
    }
    display_error_messages: false
)pb",
    "NC1=NC=NC2C3=CC=CN=C3SC=21 CHEMBL457313",
    {1, 2, 6, 7},
    kChanged, kOkCheckUsmi
  },

  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_DOUBLE_BOND
    }
    bond {
      a1: 2
      a2: 3
      btype: SS_SINGLE_BOND
    }
    bond {
      a1: 4
      a2: 5
      btype: SS_SINGLE_BOND
    }
    display_error_messages: false
)pb",
    "OC1=C2CCCC2=NC3=CC(C)=NN13 CHEMBL1331686",
    {1, 2, 2, 6, 8, 13},
    kUnChanged, kOkCheckUsmi
  },

  // This one cannot be changed because of the 3 connected nitrogen.
  ToggleTest{ R"pb(
    bond {
      a1: 0
      a2: 1
      btype: SS_SINGLE_BOND
    }
    bond {
      a1: 2
      a2: 3
      btype: SS_DOUBLE_BOND
    }
    bond {
      a1: 4
      a2: 5
      btype: SS_DOUBLE_BOND
    }
    bond {
      a1: 6
      a2: 7
      btype: SS_SINGLE_BOND
    }
    display_error_messages: false
    check_all_bonds_aromatic: false
    unset_unnecessary_implicit_hydrogens_known_values: true
)pb",
    "O=C1=C2CCCC2=NC3=CC(C)=NN13 CHEMBL1331686",
    {1, 2, 2, 6, 8, 13, 6, 7},
    kUnChanged, kNoCheckUsmi
  }
));



}  // namespace
