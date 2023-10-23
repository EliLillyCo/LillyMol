// Tester for the string classes
#include <stdlib.h>

#include <optional>
#include <unordered_set>
#include <unordered_set>
#include <string>

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "iwstring.h"

namespace {
TEST(TestIWString, TestAsString) {
  IWString s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}
TEST(TestConstIWSubstring, TestAsString) {
  const_IWSubstring s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}
TEST(TestIWString, TestExpandEnvironmentVariablesNothing) {
  const IWString hello("hello world");
  EXPECT_EQ(hello.ExpandEnvironmentVariables(), hello);
}
TEST(TestIWString, TestExpandEnvironmentVariablesTooShort) {
  const IWString hello("${}");
  EXPECT_EQ(hello.ExpandEnvironmentVariables(), hello);
}
TEST(TestIWString, TestExpandEnvironmentVariablesNotSet) {
  const IWString hello("${NOTSET}");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}
TEST(TestIWString, TestExpandEnvironmentVariablesEmpty) {
  const IWString hello("xx${}foobar");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}
TEST(TestIWString, TestExpandEnvironmentVariablesNotClosed) {
  const IWString hello("xx${xyfoobar");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}

struct EnvData {
  // Pairs of SHELL_VAR=>value
  std::unordered_map<std::string, std::string> to_set;
  // something like "hello ${world}"
  IWString input;
  // something like "hello world"
  IWString expected;
};

class TestExpandEnv: public testing::TestWithParam<EnvData> {
  protected:
    // Set of shell variables to be unset
    std::unordered_set<std::string> _to_clear;

    void TearDown();
};

void
TestExpandEnv::TearDown() {
  for (const auto& vname : _to_clear) {
    ::unsetenv(vname.c_str());
  }
}

TEST_P(TestExpandEnv, TestExpandEnv) {
  const auto params = GetParam();
  for (const auto& [vname, value] : params.to_set) {
    _to_clear.emplace(vname);
    ::setenv(vname.c_str(), value.c_str(), 1);
  }

  std::optional<IWString> expanded = params.input.ExpandEnvironmentVariables();
  EXPECT_EQ(expanded, params.expected) << " expected " << params.expected <<
        " got '" << expanded->AsString();
}
INSTANTIATE_TEST_SUITE_P(TestExpandEnv, TestExpandEnv, testing::Values(
  EnvData{{}, "hello", "hello"},
  EnvData{{{"world", "world"}}, "hello ${world}", "hello world"},
  EnvData{{{"hello", "hello"}}, "${hello} world", "hello world"},
  EnvData{{{"hello", "welcome"}}, "${hello} world", "welcome world"},
  EnvData{{{"a", "abcdefghi"}}, "${a} world", "abcdefghi world"},
  EnvData{{{"a", "abcdefghi"}}, "xxx ${a} world", "xxx abcdefghi world"},
  EnvData{{{"abcdefghi", "a"}}, "${abcdefghi} world", "a world"},
  EnvData{{{"abcdefghi", "a"}}, "xxx ${abcdefghi} y", "xxx a y"},
  EnvData{{{"mm93", "marc marquez"}, 
           {"vr46", "valentino rossi"}},
           "motogp ${mm93} and ${vr46} greats",
           "motogp marc marquez and valentino rossi greats"},
  EnvData{{{"mm93", "marc marquez"}}, "hello $mm93 motogp", "hello $mm93 motogp"},
  EnvData{{{"mm93", "marquez marquez"}}, "hello $mm93}", "hello $mm93}"}
));

}  // namespace
