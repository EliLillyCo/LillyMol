#include "ring_database_lookup.h"

#include <unistd.h>

#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "google/protobuf/text_format.h"

#include "db_cxx.h"

#include "Molecule_Tools_Bdb/smi2rings.pb.h"
#include "gtest/gtest.h"

namespace smi2rings_bdb {
namespace {
class RingDatabaseLookupTest : public testing::Test {
 protected:
  std::string _one = testing::TempDir() + "ring_lookup_one." + std::to_string(getpid());
  std::string _two = testing::TempDir() + "ring_lookup_two." + std::to_string(getpid());

  void
  TearDown() override {
    std::remove(_one.c_str());
    std::remove(_two.c_str());
  }

  void
  Write(const std::string& filename,
        const std::vector<std::pair<std::string, std::string>>& entries) {
    Db database(nullptr, DB_CXX_NO_EXCEPTIONS);
    ASSERT_EQ(database.open(nullptr, filename.c_str(), nullptr, DB_BTREE,
                            DB_CREATE | DB_TRUNCATE, 0600),
              0);
    for (const auto& [key, value] : entries) {
      Dbt dkey(const_cast<char*>(key.data()), key.size());
      Dbt data(const_cast<char*>(value.data()), value.size());
      ASSERT_EQ(database.put(nullptr, &dkey, &data, 0), 0);
    }
    ASSERT_EQ(database.close(0), 0);
  }

  std::string
  Record(uint32_t count, const std::string& exemplar) {
    Smi2Rings::Ring proto;
    proto.set_n(count);
    proto.set_ex(exemplar);
    std::string result;
    EXPECT_TRUE(google::protobuf::TextFormat::PrintToString(proto, &result));
    return result;
  }
};

TEST_F(RingDatabaseLookupTest, SingleDatabaseAndMissClearResult) {
  Write(_one, {{"ring", Record(7, "example")}});
  RingDatabaseReader reader;
  ASSERT_TRUE(reader.Open({_one}));
  RingLookupResult result;
  EXPECT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._count, 7u);
  EXPECT_EQ(result._exemplar, "example");
  EXPECT_EQ(result._databases_matched, 1u);
  EXPECT_EQ(reader.Lookup(IWString("missing"), result), LookupStatus::kNotFound);
  EXPECT_EQ(result._count, 0u);
  EXPECT_TRUE(result._exemplar.empty());
  EXPECT_EQ(result._databases_matched, 0u);
}

TEST_F(RingDatabaseLookupTest, MaximumAcrossDatabasesNotSum) {
  Write(_one, {{"ring", Record(2, "first")}, {"only_one", Record(1, "one")}});
  Write(_two, {{"ring", Record(5, "second")}});
  RingDatabaseReader reader;
  ASSERT_TRUE(reader.Open({_one, _two}));
  RingLookupResult result;
  EXPECT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._count, 5u);
  EXPECT_EQ(result._exemplar, "second");
  EXPECT_EQ(result._databases_matched, 2u);
  EXPECT_EQ(reader.Lookup(IWString("only_one"), result), LookupStatus::kFound);
  EXPECT_EQ(result._databases_matched, 1u);
}

TEST_F(RingDatabaseLookupTest, EqualCountsUseLastExemplar) {
  Write(_one, {{"ring", Record(3, "first")}});
  Write(_two, {{"ring", Record(3, "second")}});
  RingDatabaseReader reader;
  ASSERT_TRUE(reader.Open({_one, _two}));
  RingLookupResult result;
  ASSERT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._exemplar, "second");
  ASSERT_TRUE(reader.Open({_two, _one}));
  ASSERT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._exemplar, "first");
}

TEST_F(RingDatabaseLookupTest, IndependentReadersAndFailedReopen) {
  Write(_one, {{"ring", Record(1, "first")}});
  Write(_two, {{"ring", Record(2, "second")}});
  RingDatabaseReader first, second;
  ASSERT_TRUE(first.Open({_one}));
  ASSERT_TRUE(second.Open({_two}));
  EXPECT_FALSE(first.Open({_two, _one + ".missing"}));
  RingLookupResult result;
  ASSERT_EQ(first.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._count, 1u);
  ASSERT_EQ(second.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._count, 2u);
  EXPECT_FALSE(first.Open({}));
  ASSERT_EQ(first.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._count, 1u);
}

TEST_F(RingDatabaseLookupTest, ErrorIsNotAMissOrPartialSuccess) {
  Write(_one, {{"ring", Record(1, "valid")}});
  Write(_two, {{"ring", "not a ring proto"}});
  RingDatabaseReader reader;
  ASSERT_TRUE(reader.Open({_one, _two}));
  RingLookupResult result;
  EXPECT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kError);
  EXPECT_EQ(result._count, 0u);
  EXPECT_EQ(result._databases_matched, 0u);
  EXPECT_TRUE(result._exemplar.empty());
}

TEST_F(RingDatabaseLookupTest, ZeroCountRecordIsStillFound) {
  Write(_one, {{"ring", Record(0, "zero")}});
  RingDatabaseReader reader;
  ASSERT_TRUE(reader.Open({_one}));
  RingLookupResult result;
  EXPECT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kFound);
  EXPECT_EQ(result._count, 0u);
  EXPECT_EQ(result._databases_matched, 1u);
}

TEST_F(RingDatabaseLookupTest, UnopenedAndMissingDatabasesAreErrors) {
  RingDatabaseReader reader;
  RingLookupResult result;
  EXPECT_EQ(reader.Lookup(IWString("ring"), result), LookupStatus::kError);
  const std::string missing = _one + ".missing";
  EXPECT_FALSE(reader.Open({missing}));
  EXPECT_NE(access(missing.c_str(), F_OK), 0);
}
}  // namespace
}  // namespace smi2rings_bdb
