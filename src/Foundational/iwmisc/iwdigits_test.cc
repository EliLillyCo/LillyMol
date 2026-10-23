// Tester for the IWDigits class

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "iwdigits.h"

namespace {
TEST(TestIWDigits, TestTrailingSpace) {
  IWDigits iwdigits;
  iwdigits.initialise(10);
  IWString buffer;

  for (int i = 0; i < 3; ++i) {
    iwdigits.append_number(buffer, i);
  }
  EXPECT_EQ(buffer, "012");
  iwdigits.append_number(buffer, 93);
  EXPECT_EQ(buffer, "01293");

  iwdigits.append_to_each_stored_string(".");
  buffer.resize_keep_storage(0);
  for (int i = 0; i < 3; ++i) {
    iwdigits.append_number(buffer, i);
  }
  EXPECT_EQ(buffer, "0.1.2.");
  iwdigits.append_number(buffer, 2033);
  EXPECT_EQ(buffer, "0.1.2.2033.");
}


TEST(TestFractionAsString, FixedPrecisionTrimmedCachedFractions) {
  Fraction_as_String fraction;
  ASSERT_TRUE(fraction.initialise(0.0f, 1.0f, 4));

  EXPECT_EQ(fraction.string_for_fraction(0.0f), "0");
  EXPECT_EQ(fraction.string_for_fraction(0.1f), "0.1");
  EXPECT_EQ(fraction.string_for_fraction(0.2f), "0.2");
  EXPECT_EQ(fraction.string_for_fraction(0.0001f), "0.0001");
  EXPECT_EQ(fraction.string_for_fraction(1.0f), "1");
}

TEST(TestFractionAsString, FixedPrecisionTrimmedFallback) {
  Fraction_as_String fraction;
  ASSERT_TRUE(fraction.initialise(0.0f, 1.0f, 4));

  IWString buffer;
  fraction.append_number(buffer, 30.069f);
  EXPECT_EQ(buffer, "30.069");

  buffer.resize_keep_storage(0);
  fraction.append_number(buffer, -0.00001f);
  EXPECT_EQ(buffer, "0");
}

TEST(TestFractionAsString, LeadingStringPreserved) {
  Fraction_as_String fraction;
  ASSERT_TRUE(fraction.initialise(0.0f, 1.0f, 4));
  fraction.set_leading_string(' ');

  IWString buffer;
  fraction.append_number(buffer, 0.2f);
  EXPECT_EQ(buffer, " 0.2");
}


}  // namespace
