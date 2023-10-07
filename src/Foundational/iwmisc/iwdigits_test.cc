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

}  // namespace
