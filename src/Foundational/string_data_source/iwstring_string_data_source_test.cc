// Tests for iwstring_data_source

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "iwstring_string_data_source.h"

namespace {
TEST(TestStringDataSource, EmptyString) {
  String_Data_Source ds("");
  EXPECT_GT(ds.good(), 0);
  EXPECT_GT(ds.ok(), 0);
  EXPECT_EQ(ds.records_remaining(), 0);
  EXPECT_TRUE(ds.eof());
  EXPECT_EQ(ds.tellg(), 0);
}

TEST(TestStringDataSource, MultipleZeroLength) {
  const char * s = "\n\n\n\n\n";
  String_Data_Source ds(s);
  const int n = 5;
  EXPECT_EQ(ds.records_remaining(), n);
  const_IWSubstring buffer;
  for (int i = 0; i < 5; ++i) {
    ASSERT_TRUE(ds.next_record(buffer));
    EXPECT_EQ(buffer.length(), 0);
    EXPECT_EQ(ds.tellg(), i + 1);
    EXPECT_EQ(ds.records_remaining(), n - i - 1);
  }
  EXPECT_EQ(ds.lines_read(), n);
  EXPECT_TRUE(ds.eof());
}

TEST(TestStringDataSource, TestRecordsRemainingAllEmpty) {
  const char * s = "\n\n\n\n\n";
  String_Data_Source ds(s);
  EXPECT_EQ(ds.records_remaining(), 5);
}

TEST(TestStringDataSource, StripTrailingBlanks) {
  const char * s = "a\na \na   \n a \n     a      \n";
  String_Data_Source ds(s);
  ds.set_strip_trailing_blanks(1);
  const int n = 5;
  EXPECT_EQ(ds.records_remaining(), n);
  const_IWSubstring buffer;
  for (int i = 0; i < 5; ++i) {
    ASSERT_TRUE(ds.next_record(buffer));
    EXPECT_TRUE(buffer.ends_with('a'));
    const_IWSubstring stripped_here(buffer);
    stripped_here.strip_trailing_blanks();
    EXPECT_EQ(stripped_here, buffer);
  }
  EXPECT_EQ(ds.lines_read(), n);
  EXPECT_TRUE(ds.eof());
}

TEST(TestStringDataSource, TestSkipNTooFewRecords) {
  const char * s = "1\n2\n3\n4\n5\n6\n7\n8\n9\n";
  String_Data_Source ds(s);
  EXPECT_EQ(ds.records_remaining(), 9);
  EXPECT_FALSE(ds.skip_records(10));
  EXPECT_TRUE(ds.eof());
  EXPECT_EQ(ds.lines_read(), 9);
}

TEST(TestStringDataSource, TestSkipN) {
  const char * s = "abc\nabcdef\n0\n1\n2\n";
  String_Data_Source ds(s);
  ASSERT_TRUE(ds.skip_records(2));
  const_IWSubstring buffer;
  for (int i = 0; i < 3; ++i) {
    EXPECT_TRUE(ds.next_record(buffer));
    IWString as_string;
    as_string << i;
    EXPECT_EQ(buffer, as_string);
  }
  EXPECT_EQ(ds.lines_read(), 5);
  EXPECT_TRUE(ds.eof());
}

TEST(TestStringDataSource, TestSkipRegexTooFew) {
  const char * s = "vr46\nvr46\nmm93\n";
  String_Data_Source ds(s);
  re2::RE2 rx("vr46");
  ASSERT_FALSE(ds.skip_records(rx, 3));
  EXPECT_TRUE(ds.eof());
  EXPECT_EQ(ds.lines_read(), 3);
}

TEST(TestStringDataSource, TestSkipRegex) {
  const char * s = "vr46\nvr46\nmm93\nvr46\nvr46\nmm93\n";
  String_Data_Source ds(s);
  re2::RE2 rx("vr46");
  ASSERT_TRUE(ds.skip_records(rx, 3));
  const_IWSubstring buffer;
  EXPECT_TRUE(ds.next_record(buffer));
  EXPECT_EQ(buffer, "vr46");
  EXPECT_TRUE(ds.next_record(buffer));
  EXPECT_EQ(buffer, "mm93");
  EXPECT_TRUE(ds.eof());
}

TEST(TestStringDataSource, TestMostRecentRecord) {
  const char * s = "mm93\nvr46\nfoo\n\nbar  \n";
  String_Data_Source ds(s);
  ds.set_strip_trailing_blanks(1);
  const_IWSubstring buffer1;
  const_IWSubstring buffer2;
  EXPECT_TRUE(ds.next_record(buffer1));
  EXPECT_EQ(buffer1, "mm93");
  EXPECT_TRUE(ds.most_recent_record(buffer2));
  EXPECT_EQ(buffer2, "mm93");

  EXPECT_TRUE(ds.next_record(buffer1));
  EXPECT_EQ(buffer1, "vr46");
  EXPECT_TRUE(ds.most_recent_record(buffer2));
  EXPECT_EQ(buffer2, "vr46");

  EXPECT_TRUE(ds.next_record(buffer1));
  EXPECT_EQ(buffer1, "foo");

  ds.push_record();
  EXPECT_TRUE(ds.next_record(buffer2));
  EXPECT_EQ(buffer2, "foo");

}

}  // namespace
