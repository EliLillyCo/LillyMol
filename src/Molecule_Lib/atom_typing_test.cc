#include <cstdint>

#include "gtest/gtest.h"

#include "aromatic.h"
#include "atom_typing.h"

namespace {

struct SmartsForAtomTypeTestCase {
  const char* smiles;
  atom_number_t atom;
  uint32_t atom_type;
  const char* expected;
};

class SmartsForAtomTypeTest
    : public testing::TestWithParam<SmartsForAtomTypeTestCase> {
 protected:
  void SetUp() override {
    set_global_aromaticity_type(Daylight);
  }
};

TEST_P(SmartsForAtomTypeTest, GeneratesExpectedSmarts) {
  const SmartsForAtomTypeTestCase& test = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(test.smiles));
  ASSERT_GE(test.atom, 0);
  ASSERT_LT(test.atom, m.natoms());

  const IWString smarts =
      lillymol::SmartsForAtomType(m, test.atom, test.atom_type);

  EXPECT_EQ(smarts.AsString(), test.expected);
}

INSTANTIATE_TEST_SUITE_P(
    AtomProperties, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_Z, "[#6]"},
        SmartsForAtomTypeTestCase{"N", 0, IWATTYPE_USP_Z, "[#7]"},
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_A, "[A]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0, IWATTYPE_USP_A, "[a]"},
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_Z | IWATTYPE_USP_A,
                                  "[#6]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A, "[c]"},
        SmartsForAtomTypeTestCase{"Br", 0, IWATTYPE_USP_Y, "[#17]"},
        SmartsForAtomTypeTestCase{"I", 0, IWATTYPE_USP_Y, "[#17]"},
        SmartsForAtomTypeTestCase{"F", 0, IWATTYPE_USP_Y, "[#9]"},
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_N, "[*]"}));

INSTANTIATE_TEST_SUITE_P(
    Connections, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"CC(C)(C)C", 1, IWATTYPE_USP_C, "[D4]"},
        SmartsForAtomTypeTestCase{"CC(C)(C)C", 1,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A |
                                      IWATTYPE_USP_C,
                                  "[#6;D4]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A |
                                      IWATTYPE_USP_C,
                                  "[c;D2]"},
        SmartsForAtomTypeTestCase{"Cc1ccccc1", 1,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A |
                                      IWATTYPE_USP_C,
                                  "[c;D3]"}));

INSTANTIATE_TEST_SUITE_P(
    RingsAndUnsaturation, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"C1CC1", 0, IWATTYPE_USP_R, "[x2]"},
        SmartsForAtomTypeTestCase{"C1CC1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_R, "[#6;x2]"},
        SmartsForAtomTypeTestCase{"C1CC1", 0, IWATTYPE_USP_S, "[r3]"},
        SmartsForAtomTypeTestCase{"C1CC1", 0, IWATTYPE_USP_L, "[r3]"},
        SmartsForAtomTypeTestCase{"C1CC2(C1)CCCC2", 2, IWATTYPE_USP_S, "[r4]"},
        SmartsForAtomTypeTestCase{"C1CC2(C1)CCCC2", 2, IWATTYPE_USP_L, "[r5]"},
        SmartsForAtomTypeTestCase{"CC", 0, IWATTYPE_USP_U, "[G0]"},
        SmartsForAtomTypeTestCase{"C=C", 0, IWATTYPE_USP_U, "[G1]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0, IWATTYPE_USP_U, "[G1]"},
        SmartsForAtomTypeTestCase{"C=C", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_B, "[#6;G1]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_B, "[#6]"}));

INSTANTIATE_TEST_SUITE_P(
    ChargeAndCompositeTypes, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"[NH4+]", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_O, "[#7;+]"},
        SmartsForAtomTypeTestCase{"[O-]C=O", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_O, "[#8;-]"},
        SmartsForAtomTypeTestCase{"O", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_O, "[#8]"},
        SmartsForAtomTypeTestCase{
            "C1CC1", 0,
            IWATTYPE_USP_Z | IWATTYPE_USP_A | IWATTYPE_USP_C |
                IWATTYPE_USP_R | IWATTYPE_USP_U | IWATTYPE_USP_S |
                IWATTYPE_USP_L | IWATTYPE_USP_O,
            "[#6;D2;x2;G0;r3;r3]"}));

}  // namespace
