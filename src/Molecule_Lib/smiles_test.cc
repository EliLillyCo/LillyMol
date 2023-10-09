// Tests for smiles reading and generation

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "coordinate_box.h"
#include "smiles.h"

namespace {

using testing::ElementsAre;
using testing::FloatNear;

// For tests that need an input smiles and an expected output.
struct SmilesAndExpected {
  IWString input_smiles;
  IWString expected_smiles;
};

// For tests that need an input, a parameter to set, and an expected output.
struct SmilesAndParam {
  IWString input_smiles;
  int int_param;
  IWString expected_smiles;
};

class TestSmilesCharges : public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestSmilesCharges, MultiplePlusCharges) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_write_formal_charge_as_consecutive_signs(params.int_param);
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestMultipleCharges, TestSmilesCharges, testing::Values(
   SmilesAndParam{"C", 0, "C"},
   SmilesAndParam{"C", 1, "C"},
   SmilesAndParam{"[Fe+++]", 1, "[Fe+++]"},
   SmilesAndParam{"[Fe+++]", 0, "[Fe+3]"},
   SmilesAndParam{"[P---]", 1, "[P---]"},
   SmilesAndParam{"[P---]", 0, "[P-3]"},
   SmilesAndParam{"[Zn--]", 1, "[Zn--]"},
   SmilesAndParam{"[Zn--]", 0, "[Zn-2]"}
));

class TestAromaticSmiles : public testing::TestWithParam<SmilesAndExpected> {
  public:
    void SetUp() {
      set_include_aromaticity_in_smiles(1);
    }
};

TEST_P(TestAromaticSmiles, AromaticSmiles) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestAromaticSmiles, TestAromaticSmiles, testing::Values(
  SmilesAndExpected{"C", "C"},
  SmilesAndExpected{"C1=CC=CC=C1", "c1ccccc1"}
));

class TestHOnAromaticNP: public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestHOnAromaticNP, TestHOnAromaticNP) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_include_implicit_hydrogens_on_aromatic_n_and_p(params.int_param);
  EXPECT_EQ(m.unique_smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestHOnAromaticNP, TestHOnAromaticNP, testing::Values(
  SmilesAndParam{"C", 1, "C"},
  SmilesAndParam{"C1=CC=CN1", 1, "[nH]1cccc1"},
  SmilesAndParam{"C1=CC=CN1", 0, "n1cccc1"}
));

class TestAddHToIsotopes : public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestAddHToIsotopes, TestAddHToIsotopes) {
  const auto params = GetParam();
  set_unset_implicit_hydrogens_known_if_possible(params.int_param);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  EXPECT_EQ(m.unique_smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestAddHToIsotopes, TestAddHToIsotopes, testing::Values(
  SmilesAndParam{"C", 1, "C"},
  SmilesAndParam{"[2C]", 1, "[2C]"},
  SmilesAndParam{"[2C]C", 1, "[2C]C"},
  SmilesAndParam{"C[2C]C", 1, "C[2C]C"}
//SmilesAndParam{"[9C]1=CC=CC=C1", 1, "[9cH]1ccccc1"}  not sure why this does not work, investigate.
));


class TestBoxedCoordinates : public testing::Test {
    void SetUp() {
      set_append_coordinate_box_after_each_atom(1);
    }
    void TearDown() {
      set_append_coordinate_box_after_each_atom(0);
    }
};

TEST_F(TestBoxedCoordinates, TestBoxedCoordinates0) {
  const_IWSubstring smiles("C{{B0:0}}");
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));
  EXPECT_EQ(m.x(0), 0.0f);
  EXPECT_EQ(m.y(0), 0.0f);
  EXPECT_EQ(m.z(0), 0.0f);
}

TEST_F(TestBoxedCoordinates, TestBoxedCoordinates1) {
  Coordinates coords(1.0f, 2.0f, 0.5f);
  coordinate_box::ConcentricBox box;
  const coordinate_box::LayerPosition layer_position = box.Position(coords);

  IWString smiles("CC");
  smiles << "{{B" << layer_position << "}}";

  constexpr float kAbsDiff = 0.001;

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));
  EXPECT_THAT(m.atom(1).ToResizableArray(), 
              ElementsAre(FloatNear(coords.x(), kAbsDiff),
              FloatNear(coords.y(), kAbsDiff),
              FloatNear(coords.z(), kAbsDiff)));
}

TEST_F(TestBoxedCoordinates, TestBoxedCoordinates2) {
  constexpr int matoms = 10;
  std::vector<Coordinates> coords(matoms);
  for (int i = 0; i < matoms; ++i) {
    coords[i].setxyz(i, 2 * i, i + 8);
  }

  coordinate_box::ConcentricBox box;

  IWString smiles;
  for (int i = 0; i < matoms; ++i) {
    smiles << 'C';
    const coordinate_box::LayerPosition layer_position = box.Position(coords[i]);
    smiles << "{{B" << layer_position << "}}";
  }

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));

  constexpr float kAbsDiff = 0.001;

  ASSERT_EQ(m.natoms(), matoms);

  for (int i = 0; i < matoms; ++i) {
    EXPECT_THAT(m.atom(i).ToResizableArray(),
                ElementsAre(FloatNear(coords[i].x(), kAbsDiff),
                            FloatNear(coords[i].y(), kAbsDiff),
                            FloatNear(coords[i].z(), kAbsDiff)));
  }
}

struct RemoveHTest {
  IWString starting_smiles;
  IWString expected_smiles;
};

// Tests to see if we can convert an explicit Hydrogen to an
// implicit Hydrogen when a chiral centre is involved.
class TestRemoveH : public testing::TestWithParam<RemoveHTest> {
  protected:
    Molecule _mol;
};

TEST_P(TestRemoveH, Test) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.starting_smiles));
  _mol.remove_all(1);
  EXPECT_EQ(_mol.unique_smiles(), params.expected_smiles) << "Mismatch " << _mol.unique_smiles() << " vs " << params.expected_smiles << '\n';
}

INSTANTIATE_TEST_SUITE_P(TestRemoveH, TestRemoveH, testing::Values(
  RemoveHTest{"C", "C"},
  RemoveHTest{"[H][C@](F)(N)C", "F[C@@H](N)C"},
  RemoveHTest{"F[C@]([H])(C)N", "F[C@@H](N)C"},
  RemoveHTest{"N[C@](C)([H])F", "F[C@@H](N)C"},
  RemoveHTest{"C[C@](N)(F)[H]", "F[C@@H](N)C"},

  RemoveHTest{"[H][C@@](F)(N)C", "F[C@H](N)C"},
  RemoveHTest{"F[C@@]([H])(C)N", "F[C@H](N)C"},
  RemoveHTest{"N[C@@](C)([H])F", "F[C@H](N)C"},
  RemoveHTest{"C[C@@](N)(F)[H]", "F[C@H](N)C"}

));

}  // namespace
