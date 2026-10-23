#include "Utilities/General/proto_collator.h"

#include <sstream>
#include <string>

#include "gtest/gtest.h"

namespace {

struct TestRecord {
  std::string key;
  int count;
};

struct TestTraits {
  using Proto = TestRecord;
  using Value = TestRecord;

  bool IsValid(const Proto& proto, std::ostream& err) const {
    if (proto.key.empty()) {
      err << "missing key";
      return false;
    }
    return true;
  }

  const std::string& Key(const Proto& proto) const {
    return proto.key;
  }

  Value MakeValue(Proto&& proto, const std::string&) const {
    return std::move(proto);
  }

  void Merge(Value& destination, Proto&& incoming) const {
    destination.count += incoming.count;
  }

  int Write(const std::string&, const Value&, IWString_and_File_Descriptor&) const {
    return 1;
  }
};

TEST(ProtoCollator, AccumulatesDuplicateKeys) {
  proto_collator::ProtoCollator<TestTraits> collator;
  std::ostringstream err;

  EXPECT_TRUE(collator.Add(TestRecord{"a", 1}, err));
  EXPECT_TRUE(collator.Add(TestRecord{"b", 2}, err));
  EXPECT_TRUE(collator.Add(TestRecord{"a", 3}, err));

  EXPECT_EQ(collator.items_read(), 3u);
  EXPECT_EQ(collator.duplicates(), 1u);
  EXPECT_EQ(collator.size(), 2u);
  EXPECT_EQ(collator.items().at("a").count, 4);
  EXPECT_EQ(collator.items().at("b").count, 2);
}

TEST(ProtoCollator, RejectsInvalidRecords) {
  proto_collator::ProtoCollator<TestTraits> collator;
  std::ostringstream err;

  EXPECT_FALSE(collator.Add(TestRecord{"", 1}, err));
  EXPECT_EQ(collator.errors(), 1u);
  EXPECT_TRUE(collator.empty());
}

}  // namespace
