// Test code for the Matcher object.

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "matcher.h"

namespace {

TEST(TestMatcher, TestMatchAny) {
  iwmatcher::Matcher<int> matcher;
  ASSERT_TRUE(matcher.ok());
  EXPECT_TRUE(matcher.matches(0));
  EXPECT_TRUE(matcher.matches(1));
}

TEST(TestMatcher, TestMatchSingleValue) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(5);
  ASSERT_TRUE(matcher.ok());
  EXPECT_FALSE(matcher.matches(0));
  EXPECT_TRUE(matcher.matches(5));
  EXPECT_FALSE(matcher.matches(1));
  EXPECT_TRUE(matcher.matches(5));
}

TEST(TestMatcher, TestMatchMin) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(2);
  ASSERT_TRUE(matcher.ok());
  EXPECT_FALSE(matcher.matches(0));
  EXPECT_FALSE(matcher.matches(1));
  EXPECT_TRUE(matcher.matches(2));
  EXPECT_TRUE(matcher.matches(3));
}

TEST(TestMatcher, TestMatchMax) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(2);
  ASSERT_TRUE(matcher.ok());
  EXPECT_TRUE(matcher.matches(0));
  EXPECT_TRUE(matcher.matches(1));
  EXPECT_TRUE(matcher.matches(2));
  EXPECT_FALSE(matcher.matches(3));
}

TEST(TestMatcher, TestMatchBadBuildMinMax) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(1);
  matcher.set_min(2);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild1MaxViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(1);
  matcher.add(2);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild2MaxViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(1);
  matcher.add(0);
  EXPECT_TRUE(matcher.ok());
  matcher.add(2);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild1MinViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(2);
  matcher.add(1);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild2MinViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(1);
  matcher.add(2);
  EXPECT_TRUE(matcher.ok());
  matcher.add(0);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestAddIfNotAlreadPresent) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  EXPECT_EQ(matcher.add_if_not_already_present(3), 0);
  matcher.add(4);
  EXPECT_EQ(matcher.add_if_not_already_present(4), 0);
  EXPECT_TRUE(matcher.ok());
}

TEST(TestMatcher, TestFromSeveralValues) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.add(4);
  matcher.add(5);

  EXPECT_FALSE(matcher.matches(1));
  EXPECT_FALSE(matcher.matches(2));
  EXPECT_TRUE(matcher.matches(3));
  EXPECT_TRUE(matcher.matches(4));
  EXPECT_TRUE(matcher.matches(5));
  EXPECT_FALSE(matcher.matches(6));
  EXPECT_FALSE(matcher.matches(7));
}

TEST(TestMatcher, TestMatchesMinOnly) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.set_min(1);
  ASSERT_TRUE(matcher.ok());
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatcher, TestMatchesMaxOnly) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.set_max(5);
  ASSERT_TRUE(matcher.ok());
  EXPECT_FALSE(matcher.matches(6));
  EXPECT_TRUE(matcher.matches(3));
}

TEST(TestMatchesMaxOnly, TestAdjustToAccommodateIndividualValue) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  EXPECT_FALSE(matcher.matches(2));
  EXPECT_TRUE(matcher.matches(3));
  matcher.adjust_to_accommodate(2);
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatchesMaxOnly, TestAdjustToAccommodateMin) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(3);
  EXPECT_FALSE(matcher.matches(2));
  matcher.adjust_to_accommodate(2);
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatchesMaxOnly, TestAdjustToAccommodateMax) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(3);
  EXPECT_FALSE(matcher.matches(4));
  matcher.adjust_to_accommodate(2);
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatches, TestAnyValueFunctor) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.add(4);
  EXPECT_FALSE(matcher.AnyValue([] (const int v) {
    return v == 2;
  }));
  EXPECT_TRUE(matcher.AnyValue([] (const int v) {
    return v == 3;
  }));
}

TEST(TestMatcher, TestUPdateValues) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.add(4);
  EXPECT_FALSE(matcher.matches(5));
  matcher.UpdateValues([](int v) {
    return v + 1;
  });
  EXPECT_FALSE(matcher.matches(3));
  EXPECT_TRUE(matcher.matches(4));
  EXPECT_TRUE(matcher.matches(5));
  EXPECT_FALSE(matcher.matches(6));
}

}  // namespace
