#include "Molecule_Tools/molecule_filter_lib.h"

#include "gtest/gtest.h"

namespace {

using molecule_filter_lib::Feature;
using molecule_filter_lib::FeatureCalculators;
using molecule_filter_lib::FeatureFromName;
using molecule_filter_lib::FeatureValues;
using molecule_filter_lib::FeatureName;
using molecule_filter_lib::MoleculeFilter;
using molecule_filter_lib::RejectionReason;
using molecule_filter_lib::Utility;



TEST(FeatureTest, RecognisesCanonicalFeatureNames) {
  EXPECT_EQ(FeatureFromName("natoms"), Feature::kNatoms);
  EXPECT_EQ(FeatureFromName("nrings"), Feature::kNrings);
  EXPECT_EQ(FeatureFromName("heteroatom_count"), Feature::kHeteroatomCount);
  EXPECT_EQ(FeatureFromName("heteroatom_fraction"), Feature::kHeteroatomFraction);
  EXPECT_EQ(FeatureFromName("aromatic_ring_count"), Feature::kAromaticRingCount);
  EXPECT_EQ(FeatureFromName("aliphatic_ring_count"), Feature::kAliphaticRingCount);
  EXPECT_EQ(FeatureFromName("rotatable_bonds"), Feature::kRotatableBonds);
  EXPECT_EQ(FeatureFromName("max_ring_system_size"), Feature::kMaxRingSystemSize);
  EXPECT_EQ(FeatureFromName("aromatic_rings_in_system"), Feature::kAromaticRingsInSystem);
  EXPECT_EQ(FeatureFromName("tpsa"), Feature::kTpsa);
  EXPECT_EQ(FeatureFromName("alogp"), Feature::kAlogp);
  EXPECT_EQ(FeatureFromName("xlogp"), Feature::kXlogp);
  EXPECT_EQ(FeatureFromName("hba"), Feature::kHba);
  EXPECT_EQ(FeatureFromName("hbd"), Feature::kHbd);
  EXPECT_EQ(FeatureFromName("largest_ring_size"), Feature::kLargestRingSize);
  EXPECT_EQ(FeatureFromName("halogen_count"), Feature::kHalogenCount);
  EXPECT_EQ(FeatureFromName("max_distance"), Feature::kMaxDistance);
  EXPECT_EQ(FeatureFromName("sp3_carbon"), Feature::kSp3Carbon);
  EXPECT_EQ(FeatureFromName("sp3_carbon_fraction"), Feature::kSp3CarbonFraction);
  EXPECT_EQ(FeatureFromName("aromatic_density"), Feature::kAromaticDensity);
  EXPECT_EQ(FeatureFromName("chiral"), Feature::kChiral);
  EXPECT_EQ(FeatureFromName("number_fragments"), Feature::kNumberFragments);
  EXPECT_EQ(FeatureFromName("qed"), Feature::kQed);
}

TEST(FeatureTest, RecognisesLimitedAliases) {
  EXPECT_EQ(FeatureFromName("heteroatoms"), Feature::kHeteroatomCount);
  EXPECT_EQ(FeatureFromName("aromatic_rings"), Feature::kAromaticRingCount);
  EXPECT_EQ(FeatureFromName("aliphatic_rings"), Feature::kAliphaticRingCount);
  EXPECT_EQ(FeatureFromName("rotbond"), Feature::kRotatableBonds);
  EXPECT_EQ(FeatureFromName("ring_system_size"), Feature::kMaxRingSystemSize);
  EXPECT_EQ(FeatureFromName("halogens"), Feature::kHalogenCount);
  EXPECT_EQ(FeatureFromName("longest_path"), Feature::kMaxDistance);
  EXPECT_EQ(FeatureFromName("fraction_csp3"), Feature::kSp3CarbonFraction);
  EXPECT_EQ(FeatureFromName("fcsp3"), Feature::kSp3CarbonFraction);
  EXPECT_EQ(FeatureFromName("nfrag"), Feature::kNumberFragments);
  EXPECT_EQ(FeatureFromName("fragments"), Feature::kNumberFragments);
}

TEST(FeatureTest, RejectsUnknownFeatureName) {
  EXPECT_EQ(FeatureFromName("too_many_atoms"), std::nullopt);
  EXPECT_EQ(FeatureFromName("planar"), std::nullopt);
}

TEST(FeatureTest, ReturnsCanonicalFeatureName) {
  EXPECT_EQ(FeatureName(Feature::kRotatableBonds), "rotatable_bonds");
  EXPECT_EQ(FeatureName(Feature::kSp3CarbonFraction), "sp3_carbon_fraction");
  EXPECT_EQ(FeatureName(Feature::kNumberFragments), "number_fragments");
  EXPECT_EQ(FeatureName(Feature::kQed), "qed");
}

molecule_filter_data::Utility* AddUtility(molecule_filter_data::Requirements& proto,
                                          const char* name) {
  molecule_filter_data::Utility* utility = proto.add_utility();
  utility->set_name(name);
  return utility;
}

molecule_filter_data::Point* AddPoint(molecule_filter_data::Utility& proto,
                                      double x, double y) {
  molecule_filter_data::Point* point = proto.add_point();
  point->set_x(x);
  point->set_y(y);
  return point;
}

void AddLinearUtility(molecule_filter_data::Requirements& proto, const char* name,
                      double x1, double y1, double x2, double y2,
                      float weight = 1.0f) {
  molecule_filter_data::Utility* utility = AddUtility(proto, name);
  if (weight != 1.0f) {
    utility->set_weight(weight);
  }
  AddPoint(*utility, x1, y1);
  AddPoint(*utility, x2, y2);
}


TEST(UtilityTest, SortsPointsAndInterpolates) {
  molecule_filter_data::Utility proto;
  proto.set_name("natoms");
  proto.set_weight(2.5f);
  AddPoint(proto, 10.0, 1.0);
  AddPoint(proto, 0.0, 0.0);
  AddPoint(proto, 5.0, 0.5);

  Utility utility;
  ASSERT_TRUE(utility.BuildFromProto(proto));

  EXPECT_EQ(utility.name(), "natoms");
  EXPECT_EQ(utility.feature(), Feature::kNatoms);
  EXPECT_FLOAT_EQ(utility.weight(), 2.5f);
  EXPECT_EQ(utility.npoints(), 3);
  EXPECT_DOUBLE_EQ(utility.Value(0.0), 0.0);
  EXPECT_DOUBLE_EQ(utility.Value(2.5), 0.25);
  EXPECT_DOUBLE_EQ(utility.Value(5.0), 0.5);
  EXPECT_DOUBLE_EQ(utility.Value(7.5), 0.75);
  EXPECT_DOUBLE_EQ(utility.Value(10.0), 1.0);
}

TEST(UtilityTest, ClampsOutsideRange) {
  molecule_filter_data::Utility proto;
  proto.set_name("alogp");
  AddPoint(proto, -1.0, 0.2);
  AddPoint(proto, 3.0, 0.8);

  Utility utility;
  ASSERT_TRUE(utility.BuildFromProto(proto));

  EXPECT_FLOAT_EQ(utility.weight(), 1.0f);
  EXPECT_NEAR(utility.Value(-10.0), 0.2, 1.0e-6);
  EXPECT_NEAR(utility.Value(10.0), 0.8, 1.0e-6);
}

TEST(UtilityTest, AllowsDecreasingUtility) {
  molecule_filter_data::Utility proto;
  proto.set_name("natoms");
  AddPoint(proto, 10.0, 1.0);
  AddPoint(proto, 20.0, 0.0);

  Utility utility;
  ASSERT_TRUE(utility.BuildFromProto(proto));

  EXPECT_DOUBLE_EQ(utility.Value(15.0), 0.5);
}


TEST(UtilityTest, RejectsUnknownFeatureName) {
  molecule_filter_data::Utility proto;
  proto.set_name("too_many_atoms");
  AddPoint(proto, 0.0, 0.0);
  AddPoint(proto, 1.0, 1.0);

  Utility utility;
  EXPECT_FALSE(utility.BuildFromProto(proto));
}

TEST(UtilityTest, RejectsMissingName) {
  molecule_filter_data::Utility proto;
  AddPoint(proto, 0.0, 0.0);
  AddPoint(proto, 1.0, 1.0);

  Utility utility;
  EXPECT_FALSE(utility.BuildFromProto(proto));
}

TEST(UtilityTest, RejectsTooFewPoints) {
  molecule_filter_data::Utility proto;
  proto.set_name("natoms");
  AddPoint(proto, 0.0, 0.0);

  Utility utility;
  EXPECT_FALSE(utility.BuildFromProto(proto));
}

TEST(UtilityTest, RejectsDuplicateXValues) {
  molecule_filter_data::Utility proto;
  proto.set_name("natoms");
  AddPoint(proto, 0.0, 0.0);
  AddPoint(proto, 0.0, 1.0);

  Utility utility;
  EXPECT_FALSE(utility.BuildFromProto(proto));
}

TEST(UtilityTest, RejectsYOutsideUnitInterval) {
  molecule_filter_data::Utility proto;
  proto.set_name("natoms");
  AddPoint(proto, 0.0, 0.0);
  AddPoint(proto, 1.0, 1.1);

  Utility utility;
  EXPECT_FALSE(utility.BuildFromProto(proto));
}

TEST(UtilityTest, RejectsNonPositiveWeight) {
  molecule_filter_data::Utility proto;
  proto.set_name("natoms");
  proto.set_weight(0.0f);
  AddPoint(proto, 0.0, 0.0);
  AddPoint(proto, 1.0, 1.0);

  Utility utility;
  EXPECT_FALSE(utility.BuildFromProto(proto));
}



TEST(FeatureValuesTest, ComputesSimpleMolecularFeatures) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  xlogp::XLogPCalc xlogp;
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues values(m, m.natoms(), m.nrings(), calculators);

  ASSERT_TRUE(values.Value(Feature::kNatoms));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kNatoms), 3.0);
  ASSERT_TRUE(values.Value(Feature::kNrings));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kNrings), 0.0);
  ASSERT_TRUE(values.Value(Feature::kHeteroatomCount));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHeteroatomCount), 1.0);
  ASSERT_TRUE(values.Value(Feature::kHeteroatomFraction));
  EXPECT_NEAR(*values.Value(Feature::kHeteroatomFraction), 1.0 / 3.0, 1.0e-6);
  ASSERT_TRUE(values.Value(Feature::kLargestRingSize));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kLargestRingSize), 0.0);
  ASSERT_TRUE(values.Value(Feature::kHalogenCount));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHalogenCount), 0.0);
  ASSERT_TRUE(values.Value(Feature::kSp3Carbon));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kSp3Carbon), 2.0);
  ASSERT_TRUE(values.Value(Feature::kSp3CarbonFraction));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kSp3CarbonFraction), 1.0);
  ASSERT_TRUE(values.Value(Feature::kNumberFragments));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kNumberFragments), 1.0);
}

TEST(FeatureValuesTest, Sp3CarbonFractionUsesCarbonAtomDenominator) {
  Molecule ethene;
  ASSERT_TRUE(ethene.build_from_smiles("C=C"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  xlogp::XLogPCalc xlogp;
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues ethene_values(ethene, ethene.natoms(), ethene.nrings(), calculators);
  ASSERT_TRUE(ethene_values.Value(Feature::kSp3CarbonFraction));
  EXPECT_DOUBLE_EQ(*ethene_values.Value(Feature::kSp3CarbonFraction), 0.0);

  Molecule nitrogen;
  ASSERT_TRUE(nitrogen.build_from_smiles("N"));
  FeatureValues nitrogen_values(nitrogen, nitrogen.natoms(), nitrogen.nrings(), calculators);
  ASSERT_TRUE(nitrogen_values.Value(Feature::kSp3CarbonFraction));
  EXPECT_DOUBLE_EQ(*nitrogen_values.Value(Feature::kSp3CarbonFraction), 0.0);
}

TEST(FeatureValuesTest, ComputesRingFeatures) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccc2ccccc2c1"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  xlogp::XLogPCalc xlogp;
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues values(m, m.natoms(), m.nrings(), calculators);

  ASSERT_TRUE(values.Value(Feature::kNrings));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kNrings), 2.0);
  ASSERT_TRUE(values.Value(Feature::kAromaticRingCount));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kAromaticRingCount), 2.0);
  ASSERT_TRUE(values.Value(Feature::kAliphaticRingCount));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kAliphaticRingCount), 0.0);
  ASSERT_TRUE(values.Value(Feature::kLargestRingSize));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kLargestRingSize), 6.0);
  ASSERT_TRUE(values.Value(Feature::kAromaticDensity));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kAromaticDensity), 1.0);
  ASSERT_TRUE(values.Value(Feature::kMaxRingSystemSize));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kMaxRingSystemSize), 2.0);
  ASSERT_TRUE(values.Value(Feature::kAromaticRingsInSystem));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kAromaticRingsInSystem), 2.0);
}

TEST(FeatureValuesTest, ComputesGroupedRuleOfFiveFeatures) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("NCCO"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  xlogp::XLogPCalc xlogp;
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues values(m, m.natoms(), m.nrings(), calculators);

  ASSERT_TRUE(values.Value(Feature::kHba));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHba), 2.0);
  // NCCO - the NH2 contributes two donors and the OH one. Donors are a count
  // of hydrogens, not of the heteroatoms carrying them.
  ASSERT_TRUE(values.Value(Feature::kHbd));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHbd), 3.0);
}

// The RDKit compatible features are a different definition and must not be
// confused with the Lipinski ones above. Thioanisole exercises both
// differences at once - sulfur is an acceptor for RDKit but not for Lipinski,
// and the amide nitrogen of the acetamide case is an acceptor for Lipinski but
// not for RDKit. The counts themselves are tested in moleculeh_test, this is
// checking that the feature names reach the right calculation.
TEST(FeatureValuesTest, ComputesRdkitHbaHbdFeatures) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CSC"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  xlogp::XLogPCalc xlogp;
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues values(m, m.natoms(), m.nrings(), calculators);

  // Lipinski counts nitrogen and oxygen only, so a thioether is nothing.
  ASSERT_TRUE(values.Value(Feature::kHba));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHba), 0.0);
  ASSERT_TRUE(values.Value(Feature::kHbaRdkit));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHbaRdkit), 1.0);
  ASSERT_TRUE(values.Value(Feature::kHbdRdkit));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHbdRdkit), 0.0);
}

TEST(FeatureValuesTest, RdkitHbaExcludesAmideNitrogen) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(N)=O"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  xlogp::XLogPCalc xlogp;
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues values(m, m.natoms(), m.nrings(), calculators);

  // Lipinski counts both the nitrogen and the oxygen, RDKit only the oxygen.
  ASSERT_TRUE(values.Value(Feature::kHba));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHba), 2.0);
  ASSERT_TRUE(values.Value(Feature::kHbaRdkit));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHbaRdkit), 1.0);
  // Lipinski counts the two hydrogens, RDKit counts the one nitrogen.
  ASSERT_TRUE(values.Value(Feature::kHbd));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHbd), 2.0);
  ASSERT_TRUE(values.Value(Feature::kHbdRdkit));
  EXPECT_DOUBLE_EQ(*values.Value(Feature::kHbdRdkit), 1.0);
}

// min_ and max_ on a float valued feature are inclusive. Whether a fraction
// that is conceptually equal to the threshold lands just inside or just
// outside must not depend on binary representation - the threshold arrives as
// a 32 bit float while the value is computed as a double, and the direction of
// that error varies with the number.
struct SmilesThresholdOutcome {
  const char* smiles;
  float threshold;
  bool passes_min;
  bool passes_max;
};

class TestFractionThreshold : public testing::TestWithParam<SmilesThresholdOutcome> {
};

TEST_P(TestFractionThreshold, InclusiveAtTheBoundary) {
  const auto params = GetParam();

  molecule_filter_data::Requirements req_min;
  req_min.set_min_sp3_carbon_fraction(params.threshold);
  MoleculeFilter filter_min;
  ASSERT_TRUE(filter_min.Build(req_min));

  molecule_filter_data::Requirements req_max;
  req_max.set_max_sp3_carbon_fraction(params.threshold);
  MoleculeFilter filter_max;
  ASSERT_TRUE(filter_max.Build(req_max));

  Molecule m1, m2;
  ASSERT_TRUE(m1.build_from_smiles(params.smiles));
  ASSERT_TRUE(m2.build_from_smiles(params.smiles));

  EXPECT_EQ(static_cast<bool>(filter_min.Ok(m1)), params.passes_min)
      << params.smiles << " min " << params.threshold;
  EXPECT_EQ(static_cast<bool>(filter_max.Ok(m2)), params.passes_max)
      << params.smiles << " max " << params.threshold;
}
INSTANTIATE_TEST_SUITE_P(TestFractionThreshold, TestFractionThreshold, testing::Values(
  // CCCCc1ccccc1 - 4 sp3 carbons of 10, so exactly 0.4. Both bounds keep it.
  SmilesThresholdOutcome{"CCCCc1ccccc1", 0.4f, true, true},
  SmilesThresholdOutcome{"CCCCc1ccccc1", 0.35f, true, false},
  SmilesThresholdOutcome{"CCCCc1ccccc1", 0.45f, false, true},
  // CCc1ccccc1 - 2 of 8, exactly 0.25. As a float this rounds up, as a double
  // it is exact, so a bare comparison would reject it from min_.
  SmilesThresholdOutcome{"CCc1ccccc1", 0.25f, true, true},
  // Cc1ccccc1 - 1 of 7, a non terminating fraction either way.
  SmilesThresholdOutcome{"Cc1ccccc1", 0.1f, true, false},
  SmilesThresholdOutcome{"Cc1ccccc1", 0.2f, false, true}
));

// Integer requirements are exact - no tolerance is applied.
TEST(MoleculeFilterTest, IntegerThresholdsAreExact) {
  molecule_filter_data::Requirements req;
  req.set_max_natoms(3);
  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(req));

  Molecule m3, m4;
  ASSERT_TRUE(m3.build_from_smiles("CCO"));
  ASSERT_TRUE(m4.build_from_smiles("CCCO"));
  EXPECT_TRUE(filter.Ok(m3));
  EXPECT_FALSE(filter.Ok(m4));
}

TEST(FeatureFromNameTest, RdkitHbondNames) {
  EXPECT_EQ(FeatureFromName("hba_rdkit"), Feature::kHbaRdkit);
  EXPECT_EQ(FeatureFromName("hbd_rdkit"), Feature::kHbdRdkit);
  EXPECT_EQ(FeatureName(Feature::kHbaRdkit), "hba_rdkit");
  EXPECT_EQ(FeatureName(Feature::kHbdRdkit), "hbd_rdkit");
}

TEST(FeatureValuesTest, ComputesOptionalContinuousFeatures) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  quick_rotbond::QuickRotatableBonds rotbond;
  rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  alogp::ALogP alogp;
  alogp.set_use_alcohol_for_acid(1);
  alogp.set_apply_zwitterion_correction(1);
  xlogp::XLogPCalc xlogp;
  xlogp.SetIssueUnclassifiedAtomMessages(false);
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  qed::Qed qed;
  FeatureCalculators calculators{rotbond, alogp, xlogp, tpsa, qed};

  FeatureValues values(m, m.natoms(), m.nrings(), calculators);

  EXPECT_TRUE(values.Value(Feature::kRotatableBonds));
  EXPECT_TRUE(values.Value(Feature::kTpsa));
  EXPECT_TRUE(values.Value(Feature::kAlogp));
  EXPECT_TRUE(values.Value(Feature::kXlogp));
  EXPECT_TRUE(values.Value(Feature::kMaxDistance));
}

TEST(MoleculeFilterTest, BuildsUtilitiesFromProto) {
  molecule_filter_data::Requirements proto;
  molecule_filter_data::Utility* natoms = AddUtility(proto, "natoms");
  AddPoint(*natoms, 0.0, 0.0);
  AddPoint(*natoms, 10.0, 1.0);

  molecule_filter_data::Utility* alogp = AddUtility(proto, "alogp");
  alogp->set_weight(3.0f);
  AddPoint(*alogp, -1.0, 0.2);
  AddPoint(*alogp, 3.0, 0.8);

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));
  EXPECT_TRUE(filter.active());
  ASSERT_EQ(filter.number_utilities(), 2);
  EXPECT_EQ(filter.utility(0).name(), "natoms");
  EXPECT_EQ(filter.utility(0).feature(), Feature::kNatoms);
  EXPECT_DOUBLE_EQ(filter.utility(0).Value(5.0), 0.5);
  EXPECT_EQ(filter.utility(1).name(), "alogp");
  EXPECT_EQ(filter.utility(1).feature(), Feature::kAlogp);
  EXPECT_FLOAT_EQ(filter.utility(1).weight(), 3.0f);
  EXPECT_NEAR(filter.utility(1).Value(1.0), 0.5, 1.0e-6);
}

TEST(MoleculeFilterTest, RebuildClearsUtilities) {
  molecule_filter_data::Requirements proto;
  molecule_filter_data::Utility* natoms = AddUtility(proto, "natoms");
  AddPoint(*natoms, 0.0, 0.0);
  AddPoint(*natoms, 10.0, 1.0);

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));
  ASSERT_EQ(filter.number_utilities(), 1);

  molecule_filter_data::Requirements empty_proto;
  ASSERT_TRUE(filter.Build(empty_proto));
  EXPECT_TRUE(filter.active());
  EXPECT_EQ(filter.number_utilities(), 0);
}

TEST(MoleculeFilterTest, InvalidUtilityFailsBuildAndClearsUtilities) {
  molecule_filter_data::Requirements proto;
  molecule_filter_data::Utility* natoms = AddUtility(proto, "natoms");
  AddPoint(*natoms, 0.0, 0.0);
  AddPoint(*natoms, 10.0, 1.0);

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));
  ASSERT_EQ(filter.number_utilities(), 1);

  molecule_filter_data::Requirements invalid_proto;
  molecule_filter_data::Utility* invalid = AddUtility(invalid_proto, "bad");
  AddPoint(*invalid, 0.0, 0.0);

  EXPECT_FALSE(filter.Build(invalid_proto));
  EXPECT_FALSE(filter.active());
  EXPECT_EQ(filter.number_utilities(), 0);
}

TEST(MoleculeFilterTest, AppliesNewMinMaxComplementFields) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  molecule_filter_data::Requirements max_heteroatoms;
  max_heteroatoms.set_max_heteroatom_count(0);
  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(max_heteroatoms));
  RejectionReason rejection_reason;
  EXPECT_FALSE(filter.Ok(m, rejection_reason));
  EXPECT_EQ(rejection_reason, RejectionReason::kTooManyHeteroatoms);

  molecule_filter_data::Requirements min_halogens;
  min_halogens.set_min_halogen_count(1);
  ASSERT_TRUE(filter.Build(min_halogens));
  EXPECT_FALSE(filter.Ok(m, rejection_reason));
  EXPECT_EQ(rejection_reason, RejectionReason::kTooFewHalogen);

  molecule_filter_data::Requirements min_largest_ring;
  min_largest_ring.set_min_largest_ring_size(3);
  ASSERT_TRUE(filter.Build(min_largest_ring));
  EXPECT_FALSE(filter.Ok(m, rejection_reason));
  EXPECT_EQ(rejection_reason, RejectionReason::kRingTooSmall);

  molecule_filter_data::Requirements min_sp3_carbon_fraction;
  min_sp3_carbon_fraction.set_min_sp3_carbon_fraction(1.1f);
  ASSERT_TRUE(filter.Build(min_sp3_carbon_fraction));
  EXPECT_FALSE(filter.Ok(m, rejection_reason));
  EXPECT_EQ(rejection_reason, RejectionReason::kSp3CarbonFractionTooLow);

  molecule_filter_data::Requirements min_aromatic_density;
  min_aromatic_density.set_min_aromatic_density(0.1f);
  ASSERT_TRUE(filter.Build(min_aromatic_density));
  EXPECT_FALSE(filter.Ok(m, rejection_reason));
  EXPECT_EQ(rejection_reason, RejectionReason::kAromaticDensityTooLow);
}


TEST(MoleculeFilterTest, EvaluateUtilitiesWeightedAverageByDefault) {
  molecule_filter_data::Requirements proto;
  AddLinearUtility(proto, "natoms", 0.0, 0.0, 10.0, 1.0, 2.0f);
  AddLinearUtility(proto, "heteroatom_count", 0.0, 0.0, 2.0, 1.0);

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  std::vector<double> values;
  double overall;
  ASSERT_TRUE(filter.EvaluateUtilities(m, m.natoms(), m.nrings(), values, overall));

  ASSERT_EQ(values.size(), 2u);
  EXPECT_DOUBLE_EQ(values[0], 0.3);
  EXPECT_DOUBLE_EQ(values[1], 0.5);
  EXPECT_NEAR(overall, (2.0 * 0.3 + 0.5) / 3.0, 1.0e-6);
}

TEST(MoleculeFilterTest, EvaluateUtilitiesWeightedSum) {
  molecule_filter_data::Requirements proto;
  proto.set_utility_combination(molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_SUM);
  AddLinearUtility(proto, "natoms", 0.0, 0.0, 10.0, 1.0, 2.0f);
  AddLinearUtility(proto, "heteroatom_count", 0.0, 0.0, 2.0, 1.0);

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  std::vector<double> values;
  double overall;
  ASSERT_TRUE(filter.EvaluateUtilities(m, m.natoms(), m.nrings(), values, overall));
  EXPECT_NEAR(overall, 2.0 * 0.3 + 0.5, 1.0e-6);
}

TEST(MoleculeFilterTest, EvaluateUtilitiesProduct) {
  molecule_filter_data::Requirements proto;
  proto.set_utility_combination(molecule_filter_data::UTILITY_COMBINATION_PRODUCT);
  AddLinearUtility(proto, "natoms", 0.0, 0.0, 10.0, 1.0);
  AddLinearUtility(proto, "heteroatom_count", 0.0, 0.0, 2.0, 1.0);

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  std::vector<double> values;
  double overall;
  ASSERT_TRUE(filter.EvaluateUtilities(m, m.natoms(), m.nrings(), values, overall));
  EXPECT_NEAR(overall, 0.3 * 0.5, 1.0e-6);
}

TEST(MoleculeFilterTest, EvaluateUtilitiesMinMax) {
  molecule_filter_data::Requirements proto_min;
  proto_min.set_utility_combination(molecule_filter_data::UTILITY_COMBINATION_MIN);
  AddLinearUtility(proto_min, "natoms", 0.0, 0.0, 10.0, 1.0);
  AddLinearUtility(proto_min, "heteroatom_count", 0.0, 0.0, 2.0, 1.0);

  MoleculeFilter filter_min;
  ASSERT_TRUE(filter_min.Build(proto_min));

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  std::vector<double> values;
  double overall;
  ASSERT_TRUE(filter_min.EvaluateUtilities(m, m.natoms(), m.nrings(), values, overall));
  EXPECT_DOUBLE_EQ(overall, 0.3);

  molecule_filter_data::Requirements proto_max;
  proto_max.set_utility_combination(molecule_filter_data::UTILITY_COMBINATION_MAX);
  AddLinearUtility(proto_max, "natoms", 0.0, 0.0, 10.0, 1.0);
  AddLinearUtility(proto_max, "heteroatom_count", 0.0, 0.0, 2.0, 1.0);

  MoleculeFilter filter_max;
  ASSERT_TRUE(filter_max.Build(proto_max));
  ASSERT_TRUE(filter_max.EvaluateUtilities(m, m.natoms(), m.nrings(), values, overall));
  EXPECT_DOUBLE_EQ(overall, 0.5);
}

TEST(MoleculeFilterTest, EvaluateUtilitiesNoUtilitiesIsNoOp) {
  molecule_filter_data::Requirements proto;

  MoleculeFilter filter;
  ASSERT_TRUE(filter.Build(proto));
  EXPECT_FALSE(filter.has_utilities());

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  std::vector<double> values = {1.0};
  double overall = 99.0;
  ASSERT_TRUE(filter.EvaluateUtilities(m, m.natoms(), m.nrings(), values, overall));
  EXPECT_TRUE(values.empty());
  EXPECT_DOUBLE_EQ(overall, 0.0);
}

}  // namespace
