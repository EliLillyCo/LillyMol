// Tests for functions in smi.cc

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"
#include "smiles.h"

namespace {

struct NameAndResult {
  int discern_chemaxon_smiles_extensions;
  int process_quoted_smiles;
  // the full smiles record, including name
  IWString smiles;

  // True if build_from_smiles is expected to pass.
  bool expected_tf;

  IWString expected_name;
};

class TestSetName : public testing::TestWithParam<NameAndResult> {
  protected:
    Molecule _m;
};

TEST_P(TestSetName, TestSetName) {
  const auto params = GetParam();

  smiles::SetDiscernChemaxonSmilesExtensions(params.discern_chemaxon_smiles_extensions);
  smiles::SetProcessQuotedSmiles(params.process_quoted_smiles);

  bool rc = _m.build_from_smiles(params.smiles) > 0;
  EXPECT_EQ(rc, params.expected_tf);
  if (rc) {
    EXPECT_EQ(_m.name(), params.expected_name);
  }
}
INSTANTIATE_TEST_SUITE_P(TestSetName, TestSetName, testing::Values(
  NameAndResult{0, 0, "C", true, ""},
  NameAndResult{1, 0, "C", true, ""},
  NameAndResult{0, 0, "C mm93", true, "mm93"},
  NameAndResult{1, 0, "C mm93", true, "mm93"},
  NameAndResult{0, 0, "C  mm93 ", true, "mm93"},
  NameAndResult{1, 0, "C  mm93 ", true, "mm93"},
  NameAndResult{0, 0, "C |mm93 ", true, "|mm93"},
  NameAndResult{1, 0, "C |mm93 ", true, "|mm93"},
  NameAndResult{0, 0, "C |mm93| ", true, "|mm93|"},
  NameAndResult{1, 0, "C |...| ", true, ""},
  NameAndResult{1, 0, "C |...| mm93 ", true, "mm93"},
  NameAndResult{1, 1, "C |...| mm93 ", true, "mm93"},
  NameAndResult{1, 0, "\"C |...|\" mm93 ", false, "mm93"},
  NameAndResult{1, 1, "\"C |...|\" mm93 ", true, "mm93"},
  NameAndResult{1, 9, "COc1cccnc1 |(-2.5355,0.254433,-0.166002;-1.59542,-0.742704,0.102521;-0.227576,-0.517763,0.163097;0.334801,0.720379,-0.0376926;1.70145,0.934598,0.0249043;2.50314,-0.155939,0.301217;1.9648,-1.37757,0.500599;0.621235,-1.56955,0.43539),atomProp:0.atom_charge.0.37360690016521403:0.atom_idx.0:1.atom_charge.-0.73275654250267053:1.atom_idx.1:2.atom_charge.0.45541822770253759:2.atom_idx.2:3.atom_charge.0.042214187674204595:3.atom_idx.7:4.atom_charge.-0.05385077472243438:4.atom_idx.6:5.atom_charge.0.31280460106793484:5.atom_idx.5:6.atom_charge.-0.81047997202506394:6.atom_idx.4:7.atom_charge.0.29428079025576892:7.atom_idx.3| vr46", true, "vr46"}
));

TEST(TestChemaxonCoords, Test1) {
  smiles::SetDiscernChemaxonSmilesExtensions(1);
  Molecule m;
  const IWString smiles = "COc1cccnc1 |(-2.5355,0.254433,-0.166002;-1.59542,-0.742704,0.102521;-0.227576,-0.517763,0.163097;0.334801,0.720379,-0.0376926;1.70145,0.934598,0.0249043;2.50314,-0.155939,0.301217;1.9648,-1.37757,0.500599;0.621235,-1.56955,0.43539),atomProp:0.atom_charge.0.37360690016521403:0.atom_idx.0:1.atom_charge.-0.73275654250267053:1.atom_idx.1:2.atom_charge.0.45541822770253759:2.atom_idx.2:3.atom_charge.0.042214187674204595:3.atom_idx.7:4.atom_charge.-0.05385077472243438:4.atom_idx.6:5.atom_charge.0.31280460106793484:5.atom_idx.5:6.atom_charge.-0.81047997202506394:6.atom_idx.4:7.atom_charge.0.29428079025576892:7.atom_idx.3|";
  EXPECT_TRUE(m.build_from_smiles(smiles));
  EXPECT_FLOAT_EQ(m.x(0), -2.5355);
  EXPECT_FLOAT_EQ(m.y(0), 0.254433);
  EXPECT_FLOAT_EQ(m.z(0), -0.166002);

  EXPECT_FLOAT_EQ(m.x(1), -1.59542);
  EXPECT_FLOAT_EQ(m.y(1), -0.742704);
  EXPECT_FLOAT_EQ(m.z(1), 0.102521);

  EXPECT_FLOAT_EQ(m.x(2), -0.227576);
  EXPECT_FLOAT_EQ(m.y(2), -0.517763);
  EXPECT_FLOAT_EQ(m.z(2), 0.163097);

  EXPECT_FLOAT_EQ(m.x(3), 0.334801);
  EXPECT_FLOAT_EQ(m.y(3), 0.720379);
  EXPECT_FLOAT_EQ(m.z(3), -0.0376926);

  EXPECT_FLOAT_EQ(m.x(4), 1.70145);
  EXPECT_FLOAT_EQ(m.y(4), 0.934598);
  EXPECT_FLOAT_EQ(m.z(4), 0.0249043);

  EXPECT_FLOAT_EQ(m.x(5), 2.50314);
  EXPECT_FLOAT_EQ(m.y(5), -0.155939);
  EXPECT_FLOAT_EQ(m.z(5), 0.301217);

  EXPECT_FLOAT_EQ(m.x(6), 1.9648);
  EXPECT_FLOAT_EQ(m.y(6), -1.37757);
  EXPECT_FLOAT_EQ(m.z(6), 0.500599);

  EXPECT_FLOAT_EQ(m.x(7), 0.621235);
  EXPECT_FLOAT_EQ(m.y(7), -1.56955);
  EXPECT_FLOAT_EQ(m.z(7), 0.43539);
}


struct SmilesInOut {
  // The starting smiles
  IWString input;
  // The expected smiles
  IWString output;
};

class TestChemaxonSpecialAtoms : public testing::TestWithParam<SmilesInOut> {
  protected:
    Molecule m;
};

TEST_P(TestChemaxonSpecialAtoms, Test1) {
  smiles::SetDiscernChemaxonSmilesExtensions(1);
  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);

  const auto params = GetParam();
  EXPECT_TRUE(m.build_from_smiles(params.input));
  EXPECT_EQ(m.smiles(), params.output);
}
INSTANTIATE_TEST_SUITE_P(TestChemaxonSpecialAtoms, TestChemaxonSpecialAtoms, testing::Values(
  SmilesInOut{"*C(*)CC(*)CC(*)* |$;;Pol_p;;;Q_e;;;star_e;M_p$|", "[A]C([Pol])CC([Q])CC(*)[M]"}
));

}  // namespace
