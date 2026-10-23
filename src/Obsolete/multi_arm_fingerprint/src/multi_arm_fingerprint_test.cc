#include "multi_arm_fingerprint.h"

#include <vector>

#include "Molecule_Lib/aromatic.h"

#include "gtest/gtest.h"

namespace multi_arm_fingerprint {
namespace {

class MultiArmTest : public testing::Test {
 protected:
  Generator _generator;
  Options _options;
  Molecule _m;
  Sparse_Fingerprint_Creator _fp;
  Statistics _stats;

  void
  SetUp() override {
    set_global_aromaticity_type(Daylight);
  }

  void
  Generate(const char* smiles) {
    ASSERT_TRUE(_m.build_from_smiles(smiles));
    ASSERT_TRUE(_generator.SetOptions(_options));
    ASSERT_TRUE(_generator.Fingerprint(_m, _fp, nullptr, &_stats));
  }
};

TEST_F(MultiArmTest, LinearMoleculeHasNoFeatures) {
  Generate("CCCC");
  EXPECT_EQ(_stats._centres, 0u);
  EXPECT_EQ(_fp.nbits(), 0u);
}

TEST_F(MultiArmTest, ThreeArmsAreMandatory) {
  _options._max_total_length = 3;
  Generate("CC(C)C");
  EXPECT_EQ(_stats._centres, 1u);
  EXPECT_EQ(_stats._arm_paths, 3u);
  EXPECT_EQ(_stats._combinations, 1u);
  EXPECT_EQ(_fp.nbits(), 1u);
}

TEST_F(MultiArmTest, FourConnectedCentreHasOneFourArmFeature) {
  _options._max_total_length = 4;
  Generate("CC(C)(C)C");
  EXPECT_EQ(_stats._centres, 1u);
  EXPECT_EQ(_stats._combinations, 1u);
  _options._max_total_length = 3;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp));
  EXPECT_EQ(_fp.nbits(), 0u);
}

TEST_F(MultiArmTest, LengthTuplesAndMultiplicity) {
  _options._max_total_length = 5;
  Generate("CCC(CC)CC");  // Three arms, each of length one or two.
  EXPECT_EQ(_stats._arm_paths, 6u);
  EXPECT_EQ(_stats._combinations, 7u);  // (1,1,1), three permutations each of 112/122.
  int total = 0;
  for (const auto& [bit, count] : _fp.bits_found()) {
    total += count;
  }
  EXPECT_EQ(total, 7);
  EXPECT_EQ(_fp.nbits(), 3u);
  _options._presence_only = true;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp));
  for (const auto& [bit, count] : _fp.bits_found()) {
    EXPECT_EQ(count, 1);
  }
  _options._presence_only = false;
  _options._min_total_length = 5;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, nullptr, &_stats));
  EXPECT_EQ(_stats._combinations, 3u);
}

TEST_F(MultiArmTest, IndividualArmLimit) {
  _options._max_arm_length = 1;
  Generate("CCC(CC)CC");
  EXPECT_EQ(_stats._combinations, 1u);
}

TEST_F(MultiArmTest, RingClosureAndIndependentOverlappingArms) {
  _options._max_total_length = 7;
  Generate("CC1CC1");
  EXPECT_EQ(_stats._arm_paths, 7u);  // Methyl: one. Each ring direction: 1,2,closure.
  EXPECT_EQ(_stats._combinations, 9u);
  _options._ring_closures = false;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, nullptr, &_stats));
  EXPECT_EQ(_stats._arm_paths, 5u);
  EXPECT_EQ(_stats._combinations, 4u);
}

TEST_F(MultiArmTest, RootsRestrictCentresNotTraversalAndReplaceOutput) {
  ASSERT_TRUE(_m.build_from_smiles("CC(C)C(C)C"));
  std::vector<int> roots(_m.natoms(), 0);
  roots[1] = 1;
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, roots.data(), &_stats));
  EXPECT_EQ(_stats._centres, 1u);
  EXPECT_GT(_stats._combinations, 1u);
  roots[1] = 0;
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, roots.data(), &_stats));
  EXPECT_EQ(_fp.nbits(), 0u);
  EXPECT_EQ(_stats._combinations, 0u);
}

TEST_F(MultiArmTest, DegreeSelection) {
  _options._centre_degree = 3;
  Generate("CC(C)(C)C");
  EXPECT_EQ(_fp.nbits(), 0u);
}

TEST_F(MultiArmTest, InvalidConfigurationIsRejected) {
  _options._max_total_length = 2;
  EXPECT_FALSE(_generator.SetOptions(_options));
  _options = Options{};
  _options._centre_degree = 4;
  _options._max_total_length = 3;
  EXPECT_FALSE(_generator.SetOptions(_options));
  _options = Options{};
  _options._min_total_length = 5;
  _options._max_arm_length = 1;
  EXPECT_FALSE(_generator.SetOptions(_options));
}

TEST_F(MultiArmTest, SuppliedAtomTypesAndEmptyMolecule) {
  ASSERT_TRUE(_m.build_from_smiles("CC(C)C"));
  std::vector<uint32_t> types(_m.natoms(), 42);
  ASSERT_TRUE(_generator.Fingerprint(_m, types.data(), _fp));
  EXPECT_EQ(_fp.nbits(), 1u);
  EXPECT_FALSE(_generator.Fingerprint(_m, static_cast<const uint32_t*>(nullptr), _fp));
  EXPECT_EQ(_fp.nbits(), 0u);
  Molecule empty;
  EXPECT_TRUE(_generator.Fingerprint(empty, _fp));
  EXPECT_EQ(_fp.nbits(), 0u);
}

TEST_F(MultiArmTest, RandomSmilesInvariance) {
  for (const char* smiles :
       {"CC(C)(N)C", "CC1CC1", "Cc1ccccc1", "CC1CCC2(C)CCCCC2C1", "CC(C)C(C)C"}) {
    Molecule original;
    ASSERT_TRUE(original.build_from_smiles(smiles));
    Sparse_Fingerprint_Creator expected;
    ASSERT_TRUE(_generator.Fingerprint(original, expected));
    for (int i = 0; i < 30; ++i) {
      Molecule reordered;
      ASSERT_TRUE(reordered.build_from_smiles(original.random_smiles()));
      ASSERT_TRUE(_generator.Fingerprint(reordered, _fp));
      EXPECT_TRUE(_fp == expected) << smiles;
    }
  }
}

TEST_F(MultiArmTest, AromaticRepresentationsAgree) {
  Generate("Cc1ccccc1");
  Sparse_Fingerprint_Creator expected = _fp;
  Molecule kekule;
  ASSERT_TRUE(kekule.build_from_smiles("CC1=CC=CC=C1"));
  ASSERT_TRUE(_generator.Fingerprint(kekule, _fp));
  EXPECT_TRUE(_fp == expected);
}

TEST_F(MultiArmTest, FixedHashesWithExplicitAtomTypes) {
  for (const auto& item : {std::pair<const char*, uint32_t>{"CC(C)C", 3526788589u},
                           {"CC(C)(C)C", 829943902u}}) {
    ASSERT_TRUE(_m.build_from_smiles(item.first));
    std::vector<uint32_t> types(_m.natoms(), 42);
    ASSERT_TRUE(_generator.Fingerprint(_m, types.data(), _fp));
    ASSERT_EQ(_fp.nbits(), 1u);
    EXPECT_EQ(_fp.bits_found().at(item.second), 1);
  }
}

TEST_F(MultiArmTest, ClosureAwayFromRootUsesCycleSizeNotArmLength) {
  ASSERT_TRUE(_m.build_from_smiles("CC(C)CC1CC1"));
  _options._min_total_length = 7;
  _options._max_total_length = 7;
  ASSERT_TRUE(_generator.SetOptions(_options));
  std::vector<uint32_t> types(_m.natoms(), 42);
  std::vector<int> roots(_m.natoms(), 0);
  roots[1] = 1;
  ASSERT_TRUE(_generator.Fingerprint(_m, types.data(), _fp, roots.data(), &_stats));
  EXPECT_EQ(_stats._centres, 1u);
  EXPECT_EQ(_stats._arm_paths, 10u);
  EXPECT_EQ(_stats._combinations, 2u);
  ASSERT_EQ(_fp.nbits(), 1u);
  // Two directions close a three-membered ring at arm step five.
  EXPECT_EQ(_fp.bits_found().at(2196532180u), 2);
}

TEST_F(MultiArmTest, ExplicitHydrogenIsRejectedByBothEntryPoints) {
  ASSERT_TRUE(_m.build_from_smiles("[H]C(C)(C)C"));
  ASSERT_EQ(_m.natoms(1), 1);
  std::vector<uint32_t> types(_m.natoms(), 42);
  // Rejection clears a previous fingerprint and statistics as well.
  for (bool supplied_types : {false, true}) {
    _fp.hit_bit(123);
    _stats._centres = _stats._arm_paths = _stats._combinations = _stats._distinct_bits =
        99;
    if (supplied_types) {
      EXPECT_FALSE(_generator.Fingerprint(_m, types.data(), _fp, nullptr, &_stats));
    } else {
      EXPECT_FALSE(_generator.Fingerprint(_m, _fp, nullptr, &_stats));
    }
    EXPECT_EQ(_fp.nbits(), 0u);
    EXPECT_EQ(_stats._centres, 0u);
    EXPECT_EQ(_stats._arm_paths, 0u);
    EXPECT_EQ(_stats._combinations, 0u);
    EXPECT_EQ(_stats._distinct_bits, 0u);
  }
}

}  // namespace

TEST_F(MultiArmTest, SynchronisedBalancedArms) {
  _options._synchronised = true;
  Generate("CCC(CC)CC");
  EXPECT_EQ(_stats._combinations, 2u);  // Radii one and two, not eight length tuples.
  EXPECT_EQ(_fp.nbits(), 2u);
  _options._max_total_length = 5;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, nullptr, &_stats));
  EXPECT_EQ(_stats._combinations, 1u);  // Radius two needs six bonds.
}

TEST_F(MultiArmTest, SynchronisedTerminalArmsAreNotRepeated) {
  _options._synchronised = true;
  _options._max_total_length = 12;
  Generate("CC(C)C");
  EXPECT_EQ(_stats._combinations, 1u);
  Generate("CC(C)(C)C");
  EXPECT_EQ(_stats._combinations, 1u);
}

TEST_F(MultiArmTest, SynchronisedRetainsEachTerminalVariant) {
  ASSERT_TRUE(_m.build_from_smiles("CC(CC)C(C)CC"));
  std::vector<int> roots(_m.natoms(), 0);
  roots[1] = 1;
  _options._synchronised = true;
  _options._max_total_length = 12;
  _options._max_arm_length = 3;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, roots.data(), &_stats));
  EXPECT_EQ(_stats._centres, 1u);
  EXPECT_EQ(_stats._combinations, 4u);  // One at r1, two at r2, one new at r3.
  _options._max_arm_length = 1;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, roots.data(), &_stats));
  EXPECT_EQ(_stats._combinations, 1u);
}

TEST_F(MultiArmTest, SynchronisedRingClosuresAndOmittedClosures) {
  _options._synchronised = true;
  _options._max_total_length = 12;
  Generate("CC1CC1");
  EXPECT_EQ(_stats._combinations, 3u);  // Ring arms at lengths 1, 2, closure 3.
  _options._ring_closures = false;
  ASSERT_TRUE(_generator.SetOptions(_options));
  ASSERT_TRUE(_generator.Fingerprint(_m, _fp, nullptr, &_stats));
  EXPECT_EQ(_stats._combinations, 2u);
}

TEST_F(MultiArmTest, SynchronisedRandomSmilesInvariance) {
  _options._synchronised = true;
  _options._max_total_length = 12;
  _options._max_arm_length = 3;
  ASSERT_TRUE(_generator.SetOptions(_options));
  for (const char* smiles :
       {"CC(CC)C(C)CC", "CC1CC1", "Cc1ccccc1", "CC1CCC2(C)CCCCC2C1", "CC(C)(N)C"}) {
    Molecule original;
    ASSERT_TRUE(original.build_from_smiles(smiles));
    Sparse_Fingerprint_Creator expected;
    ASSERT_TRUE(_generator.Fingerprint(original, expected));
    for (int i = 0; i < 30; ++i) {
      Molecule reordered;
      ASSERT_TRUE(reordered.build_from_smiles(original.random_smiles()));
      ASSERT_TRUE(_generator.Fingerprint(reordered, _fp));
      EXPECT_TRUE(_fp == expected) << smiles;
    }
  }
}

}  // namespace multi_arm_fingerprint
