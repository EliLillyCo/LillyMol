// Test for RecordIO
#include <unistd.h>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string_view>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "google/protobuf/text_format.h"

#include "Foundational/data_source/iwrecordio.h"
#include "Foundational/data_source/proto_for_testing.pb.h"

#include "Foundational/data_source/proto_for_testing.pb.h"

namespace {

using just_for_testing::TestMessage;

class TestReadWrite : public ::testing::Test {
  protected:
    IWString _fname;

  void SetUp();
  void TearDown();
};

void
TestReadWrite::SetUp() {
  auto pid = getpid();
  _fname << "/tmp/iwrecordio_test" << pid << ".dat";
}

void
TestReadWrite::TearDown() {
  std::string_view sv(_fname.data(), _fname.length());
  std::filesystem::path mypath(sv);
  std::error_code ec;
  if (! std::filesystem::remove(mypath, ec)) {
    std::cerr << "TestReadWrite::TearDown:cannot remove '" << _fname << "' ec " << ec << '\n';
  }
}

TEST_F(TestReadWrite, TestOne) {
  iwrecordio::IWRecordIoWriter writer;
  // std::cerr << "Writing to '" <<  fname << "'\n";
  ASSERT_TRUE(writer.Open(_fname));

  TestMessage proto;
  proto.set_d(1.0);
  proto.set_f(2.0f);
  proto.set_i(3);
  proto.set_s("four");
  proto.set_u(5);

  ASSERT_TRUE(writer.WriteSerializedProto<TestMessage>(proto));

  writer.Close();

  iwrecordio::IWRecordIoReader reader;
  ASSERT_TRUE(reader.Open(_fname));

  std::optional<TestMessage> maybe_proto = reader.ReadProto<TestMessage>();
  EXPECT_THAT(maybe_proto, testing::Ne(std::nullopt));
  // EqualsProto seems not to exist outside Google, the CookBook mentions it, not sure....
  //EXPECT_THAT(proto, testing::EqualsProto(*maybe_proto));
  EXPECT_DOUBLE_EQ(proto.d(), maybe_proto->d());
  EXPECT_FLOAT_EQ(proto.f(), maybe_proto->f());
  EXPECT_EQ(proto.i(), maybe_proto->i());
  EXPECT_EQ(proto.s(), maybe_proto->s());
  EXPECT_EQ(proto.u(), maybe_proto->u());
}

TEST_F(TestReadWrite, TestMany) {
  iwrecordio::IWRecordIoWriter writer;
  // std::cerr << "Writing to '" <<  _fname << "'\n";
  ASSERT_TRUE(writer.Open(_fname));

  constexpr int kNiter = 1000;

  TestMessage proto;
  for (int i = 0; i < kNiter; ++i) {
    proto.set_d(10.0 * i);
    proto.set_f(10.0f * i + 1);
    proto.set_i(10 * i + 2);
    proto.set_s("hello world");
    proto.set_u(10 * i + 5);

    ASSERT_TRUE(writer.WriteSerializedProto<TestMessage>(proto));
  }

  writer.Close();

  iwrecordio::IWRecordIoReader reader;
  ASSERT_TRUE(reader.Open(_fname));

  for (int i = 0; i < kNiter; ++i) {
    std::optional<TestMessage> maybe_proto = reader.ReadProto<TestMessage>();
    EXPECT_THAT(maybe_proto, testing::Ne(std::nullopt));
    EXPECT_DOUBLE_EQ(maybe_proto->d(), 10.0 * i);
    EXPECT_FLOAT_EQ(maybe_proto->f(), 10.0f * i + 1);
    EXPECT_EQ(maybe_proto->i(), 10 * i + 2);
    EXPECT_EQ(maybe_proto->s(), "hello world");
    EXPECT_EQ(maybe_proto->u(), 10 * i + 5);
  }

  EXPECT_TRUE(reader.seek_zero());
  int found = 0;
  while (true) {
    std::optional<const_IWSubstring> maybe_string = reader.Next();
    if (! maybe_string) {
      break;
    }
    ++found;
  }
  EXPECT_EQ(found, kNiter);
}

}  // namespace
