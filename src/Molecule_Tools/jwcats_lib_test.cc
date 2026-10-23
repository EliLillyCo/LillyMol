#include "Molecule_Tools/jwcats_lib.h"

#include "gtest/gtest.h"

namespace {

TEST(JWCatsLib, PairNumbers) {
  EXPECT_EQ(jwcats::PairNumber(0, 0), 0);
  EXPECT_EQ(jwcats::PairNumber(0, 4), 4);
  EXPECT_EQ(jwcats::PairNumber(4, 0), 4);
  EXPECT_EQ(jwcats::PairNumber(4, 4), 14);
}

TEST(JWCatsLib, FeatureNames) {
  jwcats::JWCats cats;
  ASSERT_TRUE(cats.Initialise());
  std::vector<std::string> names = cats.FeatureNames();
  ASSERT_EQ(names.size(), 150u);
  EXPECT_EQ(names[0], "jwc_B1PAA");
  EXPECT_EQ(names[14], "jwc_B1PHH");
  EXPECT_EQ(names.back(), "jwc_B10PHH");

  cats.SetIncludeHydrophobicPairs(0);
  ASSERT_TRUE(cats.Initialise());
  names = cats.FeatureNames();
  ASSERT_EQ(names.size(), 140u);
  EXPECT_EQ(names[0], "jwc_B1PAA");
  EXPECT_EQ(names[13], "jwc_B1PNH");
  EXPECT_EQ(names[14], "jwc_B2PAA");
}

TEST(JWCatsLib, EthaneHasOneHydrophobeHydrophobePair) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC"));

  jwcats::JWCats cats;
  ASSERT_TRUE(cats.Initialise());
  jwcats::Result result;
  EXPECT_EQ(cats.Compute(m, result), jwcats::ComputeStatus::kOk);

  ASSERT_EQ(result.scaled_counts.size(), 150u);
  EXPECT_EQ(result.number_heavy_atoms, 2);
  EXPECT_EQ(result.property_count[4], 2);
  EXPECT_FLOAT_EQ(result.scaled_counts[14], 0.5f);
}

TEST(JWCatsLib, RequiresInitialise) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC"));

  jwcats::JWCats cats;
  jwcats::Result result;
  EXPECT_EQ(cats.Compute(m, result), jwcats::ComputeStatus::kNotInitialised);
  EXPECT_TRUE(cats.FeatureNames().empty());

  ASSERT_TRUE(cats.Initialise());
  EXPECT_TRUE(cats.initialised());
  EXPECT_EQ(cats.Compute(m, result), jwcats::ComputeStatus::kOk);

  cats.SetMaxBondSeparation(5);
  EXPECT_FALSE(cats.initialised());
  EXPECT_EQ(cats.Compute(m, result), jwcats::ComputeStatus::kNotInitialised);
}

}  // namespace
