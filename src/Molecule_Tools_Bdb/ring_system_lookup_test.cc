#include "ring_system_lookup.h"

#include <unistd.h>

#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "google/protobuf/text_format.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/smiles.h"

#include "db_cxx.h"

#include "Molecule_Tools_Bdb/smi2rings.pb.h"
#include "gtest/gtest.h"

namespace smi2rings_bdb {
namespace {
class RingSystemLookupTest : public testing::Test {
 protected:
  std::string _file =
      testing::TempDir() + "ring_system_lookup." + std::to_string(getpid());
  std::vector<RingSystemMatch> _results;

  void
  SetUp() override {
    set_global_aromaticity_type(Daylight);
  }

  void
  TearDown() override {
    std::remove(_file.c_str());
  }

  IWString
  Key(const char* smiles) {
    Molecule m;
    EXPECT_TRUE(m.build_from_smiles(smiles));
    return m.unique_smiles();
  }

  std::string
  Record(uint32_t count) {
    Smi2Rings::Ring proto;
    proto.set_n(count);
    proto.set_ex("exemplar");
    std::string result;
    EXPECT_TRUE(google::protobuf::TextFormat::PrintToString(proto, &result));
    return result;
  }

  void
  Write(const std::vector<std::pair<IWString, std::string>>& records) {
    Db db(nullptr, DB_CXX_NO_EXCEPTIONS);
    ASSERT_EQ(
        db.open(nullptr, _file.c_str(), nullptr, DB_BTREE, DB_CREATE | DB_TRUNCATE, 0600),
        0);
    for (const auto& [key, value] : records) {
      Dbt k(const_cast<char*>(key.rawchars()), key.length());
      Dbt v(const_cast<char*>(value.data()), value.size());
      ASSERT_EQ(db.put(nullptr, &k, &v, 0), 0);
    }
    ASSERT_EQ(db.close(0), 0);
  }

  void
  BareOnly(RingSystemLookup& lookup) {
    lookup.extractor().options().set_label_attachment_points(0);
    lookup.extractor().options().set_always_process_unsubstituted_ring(1);
  }
};

TEST_F(RingSystemLookupTest, SpecificMatchThenFallbackAndInputPreserved) {
  Write({{Key("[1CH2]1CCCCC1"), Record(1)}, {Key("C1CCCCC1"), Record(9)}});
  RingSystemLookup lookup;
  lookup.extractor().options().set_always_process_unsubstituted_ring(1);
  ASSERT_TRUE(lookup.Open({_file}));
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC1CCCCC1 query"));
  const IWString before = m.unique_smiles();
  ASSERT_TRUE(lookup.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_TRUE(_results[0]._found);
  EXPECT_EQ(_results[0]._database._count, 1u);
  EXPECT_EQ(_results[0]._representation, RingKeyType::kAttachment);
  EXPECT_EQ(m.unique_smiles(), before);
  EXPECT_EQ(m.name(), "query");
  ASSERT_TRUE(m.build_from_smiles("CC1CCC(C)CC1 fallback"));
  ASSERT_TRUE(lookup.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_EQ(_results[0]._database._count, 9u);
  EXPECT_EQ(_results[0]._representation, RingKeyType::kBare);
}

TEST_F(RingSystemLookupTest, MissingExcludedAndAcyclicMolecules) {
  Write({{Key("C1CCCCC1"), Record(9)}});
  RingSystemLookup lookup;
  BareOnly(lookup);
  ASSERT_TRUE(lookup.Open({_file}));
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CCOC1 missing"));
  ASSERT_TRUE(lookup.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_FALSE(_results[0]._found);
  EXPECT_FALSE(_results[0]._representation);
  ASSERT_TRUE(m.build_from_smiles("CC acyclic"));
  ASSERT_TRUE(lookup.Lookup(m, _results));
  EXPECT_TRUE(_results.empty());
  ASSERT_TRUE(m.build_from_smiles("C1CCCCC1 excluded"));
  lookup.extractor().options().set_max_ring_size_to_consider(5);
  ASSERT_TRUE(lookup.Lookup(m, _results));
  EXPECT_TRUE(_results.empty());
}

TEST_F(RingSystemLookupTest, RemoveInputIsotopesIsExplicitAndPreserving) {
  Write({{Key("C1CCCCC1"), Record(9)}});
  RingSystemLookup lookup;
  BareOnly(lookup);
  ASSERT_TRUE(lookup.Open({_file}));
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("[13CH2]1CCCCC1 labelled"));
  const IWString before = m.unique_smiles();
  ASSERT_TRUE(lookup.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_FALSE(_results[0]._found);
  lookup.preprocessor().options().set_remove_isotopes(true);
  ASSERT_TRUE(lookup.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_TRUE(_results[0]._found);
  EXPECT_EQ(m.unique_smiles(), before);
  ASSERT_TRUE(lookup.LookupInPlace(m, _results));
  EXPECT_EQ(m.number_isotopic_atoms(), 0);
}

TEST_F(RingSystemLookupTest, LargestFragmentIsOptional) {
  Write({{Key("C1CCCCC1"), Record(9)}});
  RingSystemLookup lookup;
  BareOnly(lookup);
  ASSERT_TRUE(lookup.Open({_file}));
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CCCCC1.C1CCOC1 mixture"));
  ASSERT_TRUE(lookup.Lookup(m, _results));
  EXPECT_EQ(_results.size(), 2u);
  lookup.preprocessor().options().set_largest_fragment(true);
  ASSERT_TRUE(lookup.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_TRUE(_results[0]._found);
  EXPECT_EQ(m.number_fragments(), 2);
}

TEST_F(RingSystemLookupTest, InstancesHaveIndependentExtractionAndPreprocessing) {
  Write({{Key("[1CH2]1CCCCC1"), Record(1)}, {Key("C1CCCCC1"), Record(9)}});
  RingSystemLookup specific, bare;
  BareOnly(bare);
  ASSERT_TRUE(specific.Open({_file}));
  ASSERT_TRUE(bare.Open({_file}));
  ASSERT_TRUE(specific.preprocessor().SetStandardisation(IWString("all")));
  EXPECT_TRUE(specific.preprocessor().standardisation().active());
  EXPECT_FALSE(bare.preprocessor().standardisation().active());
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC1CCCCC1"));
  ASSERT_TRUE(specific.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_EQ(_results[0]._database._count, 1u);
  ASSERT_TRUE(bare.Lookup(m, _results));
  ASSERT_EQ(_results.size(), 1u);
  EXPECT_EQ(_results[0]._database._count, 9u);
}

TEST_F(RingSystemLookupTest, InvalidSetupAndErrorsClearResults) {
  RingSystemLookup unopened;
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC"));
  _results.emplace_back();
  EXPECT_FALSE(unopened.Lookup(m, _results));
  EXPECT_TRUE(_results.empty());
  Write({{Key("C1CCCCC1"), Record(9)}, {Key("C1CCOC1"), "invalid ring proto"}});
  RingSystemLookup lookup;
  BareOnly(lookup);
  ASSERT_TRUE(lookup.Open({_file}));
  ASSERT_TRUE(m.build_from_smiles("C1CCCCC1.C1CCOC1"));
  EXPECT_FALSE(lookup.Lookup(m, _results));
  EXPECT_TRUE(_results.empty());
  lookup.extractor().options().set_always_process_unsubstituted_ring(0);
  EXPECT_FALSE(lookup.Lookup(m, _results));
  EXPECT_TRUE(_results.empty());
}

TEST_F(RingSystemLookupTest, KeyGenerationRestoresCoordinateSmilesSetting) {
  Write({{Key("C1CCCCC1"), Record(9)}});
  RingSystemLookup lookup;
  BareOnly(lookup);
  ASSERT_TRUE(lookup.Open({_file}));
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CCCCC1"));
  const int saved = lillymol::include_coordinates_with_smiles();
  lillymol::set_include_coordinates_with_smiles(1);
  EXPECT_TRUE(lookup.Lookup(m, _results));
  EXPECT_EQ(lillymol::include_coordinates_with_smiles(), 1);
  lillymol::set_include_coordinates_with_smiles(saved);
}
}  // namespace
}  // namespace smi2rings_bdb
