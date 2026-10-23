#include "ring_database_builder.h"

#include <unistd.h>

#include <cstdio>
#include <string>

#include "db_cxx.h"
#include "ring_database_lookup.h"

#include "gtest/gtest.h"

namespace smi2rings_bdb {
namespace {
class RingDatabaseBuilderTest : public testing::Test {
 protected:
  std::string _file = testing::TempDir() + "ring_builder." + std::to_string(getpid());

  void
  TearDown() override {
    std::remove(_file.c_str());
  }

  RingLookupResult
  Lookup(const char* key) {
    RingDatabaseReader reader;
    EXPECT_TRUE(reader.Open({_file}));
    RingLookupResult result;
    EXPECT_EQ(reader.Lookup(IWString(key), result), LookupStatus::kFound);
    return result;
  }
};

TEST_F(RingDatabaseBuilderTest, DirectAndCachedBuildsPreserveCountsAndFirstExemplar) {
  for (bool cached : {false, true}) {
    std::remove(_file.c_str());
    {
      RingDatabaseBuildOptions options;
      options.set_cache(cached);
      RingDatabaseBuilder builder;
      ASSERT_TRUE(builder.Open(_file.c_str(), options));
      ASSERT_TRUE(builder.Store(IWString("ring"), IWString("first")));
      ASSERT_TRUE(builder.Store(IWString("ring"), IWString("second")));
      ASSERT_TRUE(builder.Store(IWString("other"), IWString("third")));
      EXPECT_EQ(builder.new_ring_systems(), 2u);
      EXPECT_EQ(builder.pending_ring_systems(), cached ? 2u : 0u);
      ASSERT_TRUE(builder.Finish());
      EXPECT_EQ(builder.pending_ring_systems(), 0u);
      EXPECT_TRUE(builder.Finish());
      EXPECT_FALSE(builder.Store(IWString("ring"), IWString("too late")));
    }
    const auto result = Lookup("ring");
    EXPECT_EQ(result._count, 2u);
    EXPECT_EQ(result._exemplar, "first");
    EXPECT_EQ(Lookup("other")._count, 1u);
  }
}

TEST_F(RingDatabaseBuilderTest, DirectAppendAndLegacyCachedReplacement) {
  {
    RingDatabaseBuilder builder;
    ASSERT_TRUE(builder.Open(_file.c_str()));
    ASSERT_TRUE(builder.Store(IWString("ring"), IWString("original")));
    ASSERT_TRUE(builder.Finish());
  }
  {
    RingDatabaseBuilder builder;
    ASSERT_TRUE(builder.Open(_file.c_str()));
    ASSERT_TRUE(builder.Store(IWString("ring"), IWString("later")));
    EXPECT_EQ(builder.new_ring_systems(), 0u);
    ASSERT_TRUE(builder.Finish());
  }
  EXPECT_EQ(Lookup("ring")._count, 2u);
  EXPECT_EQ(Lookup("ring")._exemplar, "original");
  {
    RingDatabaseBuildOptions options;
    options.set_cache(true);
    RingDatabaseBuilder builder;
    ASSERT_TRUE(builder.Open(_file.c_str(), options));
    ASSERT_TRUE(builder.Store(IWString("ring"), IWString("cached")));
    ASSERT_TRUE(builder.Finish());
  }
  EXPECT_EQ(Lookup("ring")._count, 1u);
  EXPECT_EQ(Lookup("ring")._exemplar, "cached");
}

TEST_F(RingDatabaseBuilderTest, BorrowedHandleRemainsOpen) {
  Db database(nullptr, DB_CXX_NO_EXCEPTIONS);
  ASSERT_EQ(database.open(nullptr, _file.c_str(), nullptr, DB_BTREE, DB_CREATE, 0600), 0);
  {
    RingDatabaseBuilder builder(database, {});
    ASSERT_TRUE(builder.Store(IWString("ring"), IWString("example")));
    ASSERT_TRUE(builder.Finish());
  }
  IWString key("ring");
  Dbt k(const_cast<char*>(key.rawchars()), key.length());
  Dbt value;
  EXPECT_EQ(database.get(nullptr, &k, &value, 0), 0);
  EXPECT_EQ(database.close(0), 0);
}

TEST_F(RingDatabaseBuilderTest, InvalidSetupAndMalformedStoredValueFail) {
  RingDatabaseBuilder unopened;
  EXPECT_FALSE(unopened.Store(IWString("ring"), IWString("example")));
  EXPECT_FALSE(unopened.Finish());
  EXPECT_FALSE(unopened.Open((_file + ".missing/file").c_str()));
  {
    Db database(nullptr, DB_CXX_NO_EXCEPTIONS);
    ASSERT_EQ(database.open(nullptr, _file.c_str(), nullptr, DB_BTREE, DB_CREATE, 0600),
              0);
    std::string key = "ring", bad = "invalid proto";
    Dbt k(key.data(), key.size()), value(bad.data(), bad.size());
    ASSERT_EQ(database.put(nullptr, &k, &value, 0), 0);
    ASSERT_EQ(database.close(0), 0);
  }
  RingDatabaseBuilder builder;
  ASSERT_TRUE(builder.Open(_file.c_str()));
  EXPECT_FALSE(builder.Open(_file.c_str()));
  EXPECT_FALSE(builder.Store(IWString("ring"), IWString("example")));
  EXPECT_EQ(builder.new_ring_systems(), 0u);
}
}  // namespace
}  // namespace smi2rings_bdb
