#include "Foundational/iwmisc/iwre2.h"

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
TEST(TEST_FullMatch, test_IWstring) {
  std::unique_ptr<RE2> rx(new RE2("^hello world$"));
  const IWString s("hello world");
  ASSERT_TRUE(iwre2::RE2FullMatch(s, *rx));
}
TEST(TEST_FullMatch, test_const_IWSubstring) {
  std::unique_ptr<RE2> rx(new RE2("^hello.worl[d-g]$"));
  const const_IWSubstring s("hello world");
  ASSERT_TRUE(iwre2::RE2FullMatch(s, *rx));
}

TEST(TestReset, Test_IWString) {
  std::unique_ptr<RE2> rx;
  const IWString s("^hello");
  ASSERT_TRUE(iwre2::RE2Reset(rx, s));
  const IWString hello_world("hello world");
  ASSERT_TRUE(iwre2::RE2PartialMatch(hello_world, *rx));
}
TEST(TestReset, Test_const_IWSubstring) {
  std::unique_ptr<RE2> rx;
  const const_IWSubstring s("^hello.");
  ASSERT_TRUE(iwre2::RE2Reset(rx, s));
  const IWString hello_world("hello world");
  ASSERT_TRUE(iwre2::RE2PartialMatch(hello_world, *rx));
}
}  // namespace
