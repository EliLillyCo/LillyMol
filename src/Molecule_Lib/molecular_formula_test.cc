// Tests for molecular formula

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecular_formula.h"

namespace {

using molecular_formula::MolecularFormula;

struct SmilesMf {
  IWString smiles;
  int expected;
  int mf[molecular_formula::kNTypes];
};

class TestMfFail: public testing::TestWithParam<SmilesMf> {
  protected:
    MolecularFormula<uint32_t> _mf;
};

#define ALLZERO {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

TEST_P(TestMfFail, TestFailures) {
  const auto params = GetParam();
  EXPECT_EQ(_mf.Build(params.smiles), params.expected) << params.smiles;
}
INSTANTIATE_TEST_SUITE_P(TestFailures, TestMfFail, testing::Values(
  SmilesMf{"", 0, ALLZERO},
  SmilesMf{"[", 0, ALLZERO},
  SmilesMf{"C[", 0, ALLZERO},
  SmilesMf{"C[ ", 0, ALLZERO},
  SmilesMf{"C[1]", 0, ALLZERO},
  SmilesMf{"C[12]", 0, ALLZERO},
  SmilesMf{"Q", 0, ALLZERO},
  SmilesMf{"Qq", 0, ALLZERO}
));

class TestMf: public testing::TestWithParam<SmilesMf> {
  protected:
    MolecularFormula<uint32_t> _mf;
};

TEST_P(TestMf, TestWorks) {
  const auto params = GetParam();
  // std::cerr << "Testing '" << params.smiles << "'\n";
  ASSERT_EQ(_mf.Build(params.smiles), params.expected) << params.smiles;

  const uint32_t* count = _mf.mf();

  for (int i = 0; i < 10; ++i) {
    EXPECT_EQ(count[i], params.mf[i]) << "MF diff " << i << " calc " << count[i] <<
        " expect " << params.mf[i] << ' ' << params.smiles;
  }
}

INSTANTIATE_TEST_SUITE_P(TestWorks, TestMf, testing::Values(
  SmilesMf{"C", 1, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"c", 1, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"N", 1, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"n", 1, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"O", 1, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"o", 1, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"F", 1, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"P", 1, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"p", 1, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}},
  SmilesMf{"S", 1, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}},
  SmilesMf{"s", 1, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}},
  SmilesMf{"Cl", 1, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}},
  SmilesMf{"Br", 1, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}},
  SmilesMf{"I", 1, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}},
  SmilesMf{"B", 1, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}},
  SmilesMf{"[Q]", 1, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
  SmilesMf{"C[X]", 2, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
  SmilesMf{"N[13X]", 2, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
  // Note that this works - even though it is invalid.
  // Would slow things down to intercept it, and if the
  // input is valid smiles, this will never happen.
  SmilesMf{"Cr", 1, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},  // invalid smiles
  SmilesMf{"[Cr]", 1, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
  SmilesMf{"[By]", 1, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}
));

}  // namespace
