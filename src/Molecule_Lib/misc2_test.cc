// Tests for misc2

#include "misc2.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

class TestSplitOnPlusses : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    resizable_array_p<const_IWSubstring> _tokens;
};

void
TestSplitOnPlusses::SetUp() {
}

TEST_F(TestSplitOnPlusses, Empty)
{
  EXPECT_EQ(SplitOnPlusses("", _tokens, '+'), 0);
}

TEST_F(TestSplitOnPlusses, One)
{
  EXPECT_EQ(SplitOnPlusses("One", _tokens, '+'), 1);
  EXPECT_EQ(_tokens.size(), 1);
  EXPECT_EQ(*_tokens[0], "One");
}

TEST_F(TestSplitOnPlusses, EmptyFirst)
{
  EXPECT_EQ(SplitOnPlusses("+One", _tokens, '+'), 2);
  EXPECT_EQ(_tokens.size(), 2);
  EXPECT_TRUE(_tokens[0]->empty());
  EXPECT_EQ(*_tokens[1], "One");
}

TEST_F(TestSplitOnPlusses, EmptyLast)
{
  EXPECT_EQ(SplitOnPlusses("One+Two+", _tokens, '+'), 3);
  EXPECT_EQ(_tokens.size(), 3);
  EXPECT_EQ(*_tokens[0], "One");
  EXPECT_EQ(*_tokens[1], "Two");
  EXPECT_TRUE(_tokens[2]->empty());
}

TEST_F(TestSplitOnPlusses, EmptyMiddle)
{
  EXPECT_EQ(SplitOnPlusses("One++Two", _tokens, '+'), 3);
  EXPECT_EQ(_tokens.size(), 3);
  EXPECT_EQ(*_tokens[0], "One");
  EXPECT_TRUE(_tokens[1]->empty());
  EXPECT_EQ(*_tokens[2], "Two");
}

TEST_F(TestSplitOnPlusses, VeryLong)
{
  EXPECT_EQ(SplitOnPlusses("One+++++Two++Three+Four++", _tokens, '+'), 11);
  EXPECT_EQ(*_tokens[8], "Four");
}

}  // namespace
