// Tester for tp_first_pass_lib

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "tp_first_pass_lib.h"

namespace {

using lilly_medchem_rules::MCFirstPass;
using lilly_medchem_rules::MCFirstPassCounter;
using lilly_medchem_rules::InterestingAtoms;

TEST(TestInterestingAtoms, NotInterestingC) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCC"));
  EXPECT_EQ(InterestingAtoms(m), 0);
}
TEST(TestInterestingAtoms, NotInterestingN) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("NN"));
  EXPECT_EQ(InterestingAtoms(m), 0);
}
TEST(TestInterestingAtoms, NotInterestingO) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("O"));
  EXPECT_EQ(InterestingAtoms(m), 0);
}
TEST(TestInterestingAtoms, NotInterestingX) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("FBS"));
  EXPECT_EQ(InterestingAtoms(m), 0);
}
TEST(TestInterestingAtoms, InterestingO) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));
  EXPECT_GT(InterestingAtoms(m), 0);
}
TEST(TestInterestingAtoms, InterestingN) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCN"));
  EXPECT_GT(InterestingAtoms(m), 0);
}
TEST(TestInterestingAtoms, InterestingBoth) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("OCCN"));
  EXPECT_GT(InterestingAtoms(m), 0);
}

class TestTpFirstPass : public testing::Test {
  protected:
    Molecule _m;
    MCFirstPass _mc_first_pass;
    MCFirstPassCounter _counter;
    IWString _reason;

  void SetUp() {
    set_auto_create_new_elements(1);
  }
};

TEST_F(TestTpFirstPass, TestEmptyMolecule) {
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.empty_molecule, 1);
  EXPECT_EQ(_reason, "Empty molecule");
}

TEST_F(TestTpFirstPass, TestTooFewAtomsFails) {
  ASSERT_TRUE(_m.build_from_smiles("CC"));
  _mc_first_pass.set_lower_atom_count_cutoff(3);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_reason, "too_few_atoms");
}

TEST_F(TestTpFirstPass, TestTooFewAtomsOK) {
  ASSERT_TRUE(_m.build_from_smiles("CC"));
  _mc_first_pass.set_lower_atom_count_cutoff(2);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestTooManyAtomsFails) {
  ASSERT_TRUE(_m.build_from_smiles("CCCC"));
  _mc_first_pass.set_upper_atom_count_cutoff(3);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_reason, "too_many_atoms");
}

TEST_F(TestTpFirstPass, TestTooManyAtomsOK) {
  ASSERT_TRUE(_m.build_from_smiles("CCC"));
  _mc_first_pass.set_upper_atom_count_cutoff(3);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestNonPeriodicTableRejected) {
  ASSERT_TRUE(_m.build_from_smiles("[Fx].CCC"));
  _mc_first_pass.set_exclude_molecules_containing_non_periodic_table_elements(1);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_containing_non_periodic_table_elements, 1);
  EXPECT_EQ(_reason, "non_periodic_table_atom");
}

TEST_F(TestTpFirstPass, TestNonPeriodicTableOK) {
  ASSERT_TRUE(_m.build_from_smiles("[Fx].CCC"));
  _mc_first_pass.set_exclude_molecules_containing_non_periodic_table_elements(0);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestNonPeriodicTableDisconnectedOK) {
  ASSERT_TRUE(_m.build_from_smiles("[Fx].CCC"));
  _mc_first_pass.set_exclude_molecules_containing_non_periodic_table_elements(1);
  _mc_first_pass.set_allow_non_periodic_table_elements_if_not_connected(1);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestNonPeriodicTableDisconnectedNotOk) {
  ASSERT_TRUE(_m.build_from_smiles("[Fx].CCC"));
  _mc_first_pass.set_exclude_molecules_containing_non_periodic_table_elements(1);
  _mc_first_pass.set_allow_non_periodic_table_elements_if_not_connected(0);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_containing_non_periodic_table_elements, 1);
  EXPECT_EQ(_reason, "non_periodic_table_atom");
}

TEST_F(TestTpFirstPass, TestOKElementsDefault) {
  ASSERT_TRUE(_m.build_from_smiles("[Fe].OC"));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_containing_non_allowed_atom_types, 1);
  EXPECT_EQ(_reason, "non_allowed_atom");
}

TEST_F(TestTpFirstPass, TestOKElementsIronOk) {
  ASSERT_TRUE(_m.build_from_smiles("[Fe].OC"));
  _mc_first_pass.set_ok_element(get_element_from_symbol_no_case_conversion("Fe"));
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestMinFractionInterestingAtomsZero) {
  ASSERT_TRUE(_m.build_from_smiles("SC"));
  _mc_first_pass.set_min_fraction_interesting_atoms(0.0);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestMinFractionInterestingAtomsBelow5) {
  ASSERT_TRUE(_m.build_from_smiles("NC"));
  _mc_first_pass.set_min_fraction_interesting_atoms(0.49);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestMinFractionInterestingAtomsSNotInteresting) {
  ASSERT_TRUE(_m.build_from_smiles("SC"));
  _mc_first_pass.set_min_fraction_interesting_atoms(0.1);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_too_few_interesting_atoms, 1);
  EXPECT_EQ(_reason, "too few interesting atoms");
}

TEST_F(TestTpFirstPass, TestMinFractionInterestingAtomsOK) {
  ASSERT_TRUE(_m.build_from_smiles("OCN"));
  _mc_first_pass.set_min_fraction_interesting_atoms(0.6);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestMaxFragSizeDoesNotApplyToSingleFragmentMolecules) {
  ASSERT_TRUE(_m.build_from_smiles("OCN"));
  _mc_first_pass.set_reject_if_any_fragment_larger_than(1);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestMaxFragSizeOkSize) {
  ASSERT_TRUE(_m.build_from_smiles("OCN.C"));
  _mc_first_pass.set_reject_if_any_fragment_larger_than(2);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_fragment_too_large, 1);
  EXPECT_THAT(_reason.AsString(), testing::HasSubstr("fragment too large"));
}

TEST_F(TestTpFirstPass, TestIsMixture) {
  ASSERT_TRUE(_m.build_from_smiles("OCN.CCC"));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.mixtures_rejected, 1);
  EXPECT_EQ(_reason, "mixture");
}

TEST_F(TestTpFirstPass, TestLowerRingCountCutoffBelow) {
  ASSERT_TRUE(_m.build_from_smiles("OCN"));
  _mc_first_pass.set_lower_ring_count_cutoff(1);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_too_few_rings, 1);
  EXPECT_EQ(_reason, "not_enough_rings");
}

TEST_F(TestTpFirstPass, TestLowerRingCountCutoffAbove) {
  ASSERT_TRUE(_m.build_from_smiles("O1CN1"));
  _mc_first_pass.set_lower_ring_count_cutoff(1);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestUpperRingCountCutoffAbove) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC1CC1CC1"));
  _mc_first_pass.set_upper_ring_count_cutoff(1);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_too_many_rings, 1);
  EXPECT_EQ(_reason, "too_many_rings");
}

TEST_F(TestTpFirstPass, TestUpperRingCountCutoffBelow) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC1CC1CC1"));
  _mc_first_pass.set_upper_ring_count_cutoff(2);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestIsotopeReject) {
  ASSERT_TRUE(_m.build_from_smiles("[1CH4]"));
  _mc_first_pass.set_exclude_isotopes(1);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_containing_isotopes, 1);
  EXPECT_EQ(_reason, "isotopes");
}

TEST_F(TestTpFirstPass, TestIsotopeAllowed) {
  ASSERT_TRUE(_m.build_from_smiles("[1CH4]"));
  _mc_first_pass.set_exclude_isotopes(0);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_m.unique_smiles(), "[1CH4]");
}

TEST_F(TestTpFirstPass, TestConvertIsotopes) {
  ASSERT_TRUE(_m.build_from_smiles("[1CH4]"));
  _mc_first_pass.set_convert_isotopes(1);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_m.unique_smiles(), "C");
}

TEST_F(TestTpFirstPass, TestInvalidValence) {
  ASSERT_TRUE(_m.build_from_smiles("CC(C)(C)(C)C"));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_abnormal_valences, 1);
  EXPECT_EQ(_reason, "abnormal_valence");
}

TEST_F(TestTpFirstPass, TestMaxRingBondRatioFails) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC1CC"));
  _mc_first_pass.set_max_ring_bond_ratio(0.5);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_bad_max_ring_bond_ratios, 1);
  EXPECT_EQ(_reason, "bad_ring_bond_ratio");
}

TEST_F(TestTpFirstPass, TestMaxRingBondRatioPasses) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC1CCCC"));
  _mc_first_pass.set_max_ring_bond_ratio(0.5);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestRingSizeBad) {
  ASSERT_TRUE(_m.build_from_smiles("C1CCC1"));
  _mc_first_pass.set_upper_ring_size_cutoff(3);
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_counter.molecules_with_ring_sizes_out_of_range, 1);
  EXPECT_EQ(_reason, "ring size");
}

TEST_F(TestTpFirstPass, TestRingSizeOk) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC1CCCC"));
  _mc_first_pass.set_upper_ring_size_cutoff(4);
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestElementCountNotRejected) {
  ASSERT_TRUE(_m.build_from_smiles("CCCN"));
  const char argc = 3;
  const char * argv[] = {"notused", "-F", "7:1"};
  Command_Line cl(argc, const_cast<char**>(argv), "F:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestElementCountRejected) {
  ASSERT_TRUE(_m.build_from_smiles("CCCNN"));
  const char argc = 3;
  const char * argv[] = {"notused", "-F", "7:1"};
  Command_Line cl(argc, const_cast<char **>(argv), "F:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_reason, "element_fraction");
}

TEST_F(TestTpFirstPass, TestElementFractionNotRejected) {
  ASSERT_TRUE(_m.build_from_smiles("CCCN"));
  const char argc = 3;
  const char * argv[] = {"notused", "-F", "7:0.26"};
  Command_Line cl(argc, const_cast<char**>(argv), "F:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestElementFractionRejected) {
  ASSERT_TRUE(_m.build_from_smiles("CCCNN"));
  const char argc = 3;
  const char * argv[] = {"notused", "-F", "7:0.26"};
  Command_Line cl(argc, const_cast<char**>(argv), "F:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_reason, "element_fraction");
}

TEST_F(TestTpFirstPass, TestRingSystemTooLargePasses) {
  ASSERT_TRUE(_m.build_from_smiles("C1C23CC3CC12"));
  const char argc = 3;
  const char * argv[] = {"notused", "-H", "3"};
  Command_Line cl(argc, const_cast<char**>(argv), "H:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_FALSE(_mc_first_pass.Rejected(_m, _counter, _reason));
}

TEST_F(TestTpFirstPass, TestRingSystemTooLargeFails) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC2C3C4CCC4C3C21"));
  const char argc = 3;
  const char * argv[] = {"notused", "-H", "3"};
  Command_Line cl(argc, const_cast<char**>(argv), "H:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_reason, "ring_system_size");
}

TEST_F(TestTpFirstPass, TestRejectPhosphorus) {
  ASSERT_TRUE(_m.build_from_smiles("c1nccccc1P(=O)O"));
  const char argc = 3;
  const char * argv[] = {"notused", "-n", "P"};
  Command_Line cl(argc, const_cast<char**>(argv), "n:");
  ASSERT_TRUE(_mc_first_pass.Build(cl));
  EXPECT_TRUE(_mc_first_pass.Rejected(_m, _counter, _reason));
  EXPECT_EQ(_reason, "covalent_non-organic") << _reason;
}


}  // namespace
