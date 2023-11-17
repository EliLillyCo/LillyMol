// Tester for the string classes
#include <stdlib.h>

#include <optional>
#include <string>
#include <type_traits>
#include <unordered_set>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

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

TEST(TestIWString, TestEnsureEndsWithEmpty) {
  IWString s;
  EXPECT_TRUE(s.EnsureEndsWith('a'));
  EXPECT_EQ(s, 'a');
}

template <typename T>
class EnsureEndsWithTest : public testing::Test {
 protected:
  void TestEmpty();
  void TestAlreadyEndsWith();
  void TestMustBeAdded();
};

template <typename T>
void
EnsureEndsWithTest<T>::TestEmpty() {
  IWString s;
  T extra("a");
  EXPECT_EQ(s.EnsureEndsWith(extra), 1);
  EXPECT_EQ(s, "a");
}

template <typename T>
void
EnsureEndsWithTest<T>::TestAlreadyEndsWith() {
  IWString s('a');
  T extra("a");
  EXPECT_EQ(s.EnsureEndsWith(extra), 0);
  EXPECT_EQ(s, 'a');
}

template <typename T>
void
EnsureEndsWithTest<T>::TestMustBeAdded() {
  IWString s('a');
  T extra("b");
  EXPECT_EQ(s.EnsureEndsWith(extra), 1);
  EXPECT_EQ(s, "ab");
}

using MyTypes = ::testing::Types<const char*, const IWString&, const const_IWSubstring&>;
TYPED_TEST_SUITE_P(EnsureEndsWithTest);

TYPED_TEST_P(EnsureEndsWithTest, StartsEmpty) {
  // Inside a test, refer to TypeParam to get the type parameter.
  // TypeParam n = 0;
  // Maybe something could be done to combine the const char* type into the template?
  // std::is_pointer...
  // if (std::is_integral<TypeParam>::value) {
  // }

  // You will need to use `this` explicitly to refer to fixture members.
  this->TestEmpty();
}

TYPED_TEST_P(EnsureEndsWithTest, AlreadyEndsWith) {
  // TypeParam n = 0;

  this->TestAlreadyEndsWith();
}

TYPED_TEST_P(EnsureEndsWithTest, MustBeAdded) { 
  // TypeParam n = 0;

  this->TestMustBeAdded();
}

REGISTER_TYPED_TEST_SUITE_P(EnsureEndsWithTest,
                            StartsEmpty, AlreadyEndsWith, MustBeAdded);

INSTANTIATE_TYPED_TEST_SUITE_P(My, EnsureEndsWithTest, MyTypes);



}  // namespace
