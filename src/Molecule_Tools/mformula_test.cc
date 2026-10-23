// Tests for MFormula

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwstring/iwstring.h"

#include "mformula.h"

namespace {

struct SmilesFormula {
  IWString smiles;
  int carbon;
  int nitrogen;
  int oxygen;
  int fluorine;
  int phosphorus;
  int sulphur;
  int chlorine;
  int bromine;
  int iodine;

  int natoms;

//SmilesFormula();
};

#ifdef HAS_CONSTRUCTOR
SmilesFormula::SmilesFormula() {
  carbon = 0;
  nitrogen = 0;
  oxygen = 0;
  fluorine = 0;
  phosphorus = 0;
  sulphur = 0;
  chlorine = 0;
  bromine = 0;
  iodine = 0;
  natoms = 0;
}
#endif

class TestMFormula: public testing::TestWithParam<SmilesFormula> {
  protected:
    mformula::MFormula _formula;
};

TEST_P(TestMFormula, Tests) {
  const auto params = GetParam();
  EXPECT_EQ(_formula.Build(params.smiles), params.natoms);
  EXPECT_EQ(_formula.Carbon(), params.carbon);
  EXPECT_EQ(_formula.Nitrogen(), params.nitrogen) << params.smiles;
  EXPECT_EQ(_formula.Oxygen(), params.oxygen);
  EXPECT_EQ(_formula.Fluorine(), params.fluorine);
  EXPECT_EQ(_formula.Phosphorus(), params.phosphorus);
  EXPECT_EQ(_formula.Sulphur(), params.sulphur) << params.smiles;
  EXPECT_EQ(_formula.Chlorine(), params.chlorine) << params.smiles;
  EXPECT_EQ(_formula.Bromine(), params.bromine);
  EXPECT_EQ(_formula.Iodine(), params.iodine);

  EXPECT_TRUE(_formula.initialised());
}
INSTANTIATE_TEST_SUITE_P(TestMFormula, TestMFormula, testing::Values(
  SmilesFormula{"C", 1, 0, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"c", 1, 0, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"N", 0, 1, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"n", 0, 1, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"O", 0, 0, 1, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"o", 0, 0, 1, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"F", 0, 0, 0, 1, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"P", 0, 0, 0, 0, 1, 0, 0, 0, 0,   1},
  SmilesFormula{"S", 0, 0, 0, 0, 0, 1, 0, 0, 0,   1},
  SmilesFormula{"s", 0, 0, 0, 0, 0, 1, 0, 0, 0,   1},
  SmilesFormula{"Cl", 0, 0, 0, 0, 0, 0, 1, 0, 0,   1},
  SmilesFormula{"Br", 0, 0, 0, 0, 0, 0, 0, 1, 0,   1},
  SmilesFormula{"I", 0, 0, 0, 0, 0, 0, 0, 0, 1,   1},

  SmilesFormula{"CC", 2, 0, 0, 0, 0, 0, 0, 0, 0,   2},
  SmilesFormula{"CCC", 3, 0, 0, 0, 0, 0, 0, 0, 0,   3},
  SmilesFormula{"CN", 1, 1, 0, 0, 0, 0, 0, 0, 0,   2},
  SmilesFormula{"CCl", 1, 0, 0, 0, 0, 0, 1, 0, 0,   2},
  SmilesFormula{"ClC", 1, 0, 0, 0, 0, 0, 1, 0, 0,   2},
  SmilesFormula{"[Cl][C]", 1, 0, 0, 0, 0, 0, 1, 0, 0,   2},
  SmilesFormula{"[4NH2+][12C]", 1, 1, 0, 0, 0, 0, 0, 0, 0,   2}
));

}  // namespace
