// Unit tests

#include <optional>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwmisc/toml_support.h"

#include "Foundational/iwmisc/for_testing.pb.h"

namespace {

using testing::ElementsAre;
using testing::Eq;

using for_testing::TestMessage;

TEST(Test, TestStrings) {
  const std::string toml = R"(
str1 = "hello"
str2 = "world"
)";
  std::optional<TestMessage> msg = iwmisc::ParseFromToml<TestMessage>(toml);
  EXPECT_NE(msg, std::nullopt);

  EXPECT_EQ(msg->str1(), "hello");
  EXPECT_EQ(msg->str2(), "world");
}


TEST(Test, TestInt) {
  const std::string toml = R"(
i1 = -12
)";
  std::optional<TestMessage> msg = iwmisc::ParseFromToml<TestMessage>(toml);
  EXPECT_NE(msg, std::nullopt);

  EXPECT_EQ(msg->i1(), -12);
}

TEST(Test, TestUInt) {
  const std::string toml = R"(
ui1 = 12
)";
  std::optional<TestMessage> msg = iwmisc::ParseFromToml<TestMessage>(toml);
  EXPECT_NE(msg, std::nullopt);

  EXPECT_EQ(msg->ui1(), 12);
}

TEST(Test, TestFloat) {
  const std::string toml = R"(
x = 12.0
)";
  std::optional<TestMessage> msg = iwmisc::ParseFromToml<TestMessage>(toml);
  EXPECT_NE(msg, std::nullopt);

  EXPECT_FLOAT_EQ(msg->x(), static_cast<float>(12.0));
}

TEST(Test, TestIntArray) {
  const std::string toml = R"(
int_array = [1,2,3,4,5]
)";
  std::optional<TestMessage> msg = iwmisc::ParseFromToml<TestMessage>(toml);
  EXPECT_NE(msg, std::nullopt);

  EXPECT_THAT(msg->int_array(), ElementsAre(1,2,3,4,5));
}

TEST(Test, TestSubMessage) {
  const std::string toml = R"(
repeated_string = ["hello", "world"]
[sub_message]
i1 = 3
str1 = "vr46"
)";
  std::optional<TestMessage> msg = iwmisc::ParseFromToml<TestMessage>(toml);
  EXPECT_NE(msg, std::nullopt);

  EXPECT_THAT(msg->repeated_string(), ElementsAre("hello", "world"));
  EXPECT_THAT(msg->sub_message().i1(), Eq(3));
  EXPECT_THAT(msg->sub_message().str1(), Eq("vr46"));
}

}  // namespace
