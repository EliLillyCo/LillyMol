// Tester for iwstring_data_source

#include <filesystem>
#include <fstream>

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "iwstring_data_source.h"

namespace {

using std::endl;
using std::cerr;

int CreateTempFile(const IWString& contents, IWString& fname)
{
  const auto f = std::filesystem::temp_directory_path();

  const std::string f_string(f);
  IWString path(f_string);
  path << '/' << fname;
  fname = path;

  const std::string string_fname(fname.data(), fname.length());

  std::ofstream output;
  output.open(string_fname);
  if (! output.good()) {
    cerr << "Cannot open " << fname << "\n";
    return 0;
  }
  
  output.write(contents.data(), contents.length());
  if (! output.good()) {
    cerr << "CreateTempFile:output failed\n";
    return 0;
  }

  return 1;
}

TEST(TestDataSource, Read1)
{
  const IWString contents("hello world\ngoodbye world");
  IWString fname("foo");
  CreateTempFile(contents, fname);
  cerr << "Contents written to " << fname << endl;

  iwstring_data_source input(fname);
  ASSERT_TRUE(input.ok());

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_TRUE("hello world" == buffer);
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_TRUE("goodbye world" == buffer);
}

TEST(TestDataSource, Filter)
{
  const IWString contents("hello world\ngoodbye world\nhello again\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);
  cerr << "Contents written to " << fname << endl;

  iwstring_data_source input(fname);
  ASSERT_TRUE(input.set_filter_pattern("hello"));
  ASSERT_TRUE(input.ok());

  int records_found = 0;
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    records_found++;
  }
  EXPECT_EQ(records_found, 2);
}

TEST(TestDataSource, Exclude)
{
  const IWString contents("hello world\ngoodbye world\nhello again\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);
  cerr << "Contents written to " << fname << endl;

  iwstring_data_source input(fname);
  ASSERT_TRUE(input.set_ignore_pattern("hello"));
  ASSERT_TRUE(input.ok());

  int records_found = 0;
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    records_found++;
  }
  EXPECT_EQ(records_found, 1);
}

TEST(TestDataSource, Grep)
{
  const IWString contents("hello world\ngoodbye world\nhello again\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);
  cerr << "Contents written to " << fname << endl;

  iwstring_data_source input(fname);
  EXPECT_EQ(input.grep("hello"), 2);

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(input.grep("hello"), 1);
}

TEST(TestDataSource, LeadingBlanks)
{
  const IWString contents("  hello world\n    goodbye world\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  input.set_strip_leading_blanks(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    EXPECT_FALSE(buffer.starts_with(' '));
  }
}

TEST(TestDataSource, TrailingBlanks)
{
  const IWString contents("hello world \ngoodbye world    \n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  input.set_strip_trailing_blanks(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    EXPECT_FALSE(buffer.ends_with(' '));
  }
}

TEST(TestDataSource, RecordsRemaining) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);

  const_IWSubstring buffer;
  for (int i = 0; i <= 4; ++i) {
    ASSERT_TRUE(input.next_record(buffer));
    EXPECT_EQ(buffer.length(), 1);
    EXPECT_EQ(buffer[0] - '0', i);
    EXPECT_EQ(input.records_remaining(), 4 - i);
  }
}

TEST(TestDataSource, SkipRecordsNumeric) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  RE2 rx("^[0-9]");
  ASSERT_TRUE(input.skip_records(1));

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "1");

  ASSERT_TRUE(input.skip_records(2));
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "4");
}

TEST(TestDataSource, SkipRecordsRx) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  RE2 rx("^[0-9]");
  ASSERT_TRUE(input.skip_records(rx, 1));

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "1");

  ASSERT_TRUE(input.skip_records(rx, 2));
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "4");
}

TEST(TestDataSource, CountRecordsStartingWith) 
{
  const IWString contents("0\n1\n0\n0\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);

  EXPECT_EQ(input.count_records_starting_with("0"), 3);
  EXPECT_EQ(input.count_records_starting_with("1"), 1);
  EXPECT_EQ(input.count_records_starting_with("4"), 1);

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));

  EXPECT_EQ(input.count_records_starting_with("0"), 2);
  EXPECT_EQ(input.count_records_starting_with("1"), 1);

  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(input.count_records_starting_with("0"), 2);
  EXPECT_EQ(input.count_records_starting_with("4"), 1);
}

TEST(TestDataSource, SkipTo) 
{
  const IWString contents("0\n1\n0\n0\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  ASSERT_TRUE(input.skip_to("^4"));

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "4");
}

TEST(TestDataSource, SkipPast) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  ASSERT_TRUE(input.skip_past("^1"));

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "2");

  ASSERT_TRUE(input.skip_past("^4"));
  EXPECT_FALSE(input.next_record(buffer));
}

TEST(TestDataSource, NextRecordMatches) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);
  EXPECT_TRUE(input.next_record_matches("^0"));

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "0");
  EXPECT_TRUE(input.next_record_matches("^1$"));
}

TEST(TestDataSource, PushRecord) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);

  const_IWSubstring buffer;

  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "0");
  input.push_record();
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "0");

  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "1");
  input.push_record();
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "1");
}

TEST(TestDataSource, SeekG) 
{
  const IWString contents("0\n1\n2\n3\n4\n");
  IWString fname("foo");
  CreateTempFile(contents, fname);

  iwstring_data_source input(fname);

  input.seekg(2);

  const_IWSubstring buffer;
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "1");

  input.push_record();
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "1");

  input.seekg(0);
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "0");

  input.seekg(8);
  ASSERT_TRUE(input.next_record(buffer));
  EXPECT_EQ(buffer, "4");
}


}  // namespace
