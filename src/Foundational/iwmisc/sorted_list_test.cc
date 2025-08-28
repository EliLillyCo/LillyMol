// Tests for sorted list class

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "sorted_list.h"

namespace {

struct Data {
  resizable_array<uint64_t> data;
  uint64_t extra;
  std::vector<uint64_t> expected;
};

class TestSortedList: public testing::TestWithParam<Data>  {
};

TEST_P(TestSortedList, Tests) {
  const auto params = GetParam();

  resizable_array<uint64_t> mydata(params.data);
  sorted_list::InsertIntoList(mydata, params.extra);
  EXPECT_THAT(params.expected, testing::ElementsAreArray(mydata.rawdata(), mydata.size()));
}
INSTANTIATE_TEST_SUITE_P(TestSortedList, TestSortedList, testing::Values(
  Data{ {1, 3}, 2, {1, 2, 3}},
  Data{ {1, 3, 4}, 2, {1, 2, 3, 4}},
  Data{ {1, 2, 4}, 3, {1, 2, 3, 4}},
  Data{ {1, 3, 4}, 2, {1, 2, 3, 4}},
  Data{ {1, 3, 4, 5}, 2, {1, 2, 3, 4, 5}},
  Data{ {1, 2, 4, 5}, 3, {1, 2, 3, 4, 5}},
  Data{ {1, 2, 3, 5}, 4, {1, 2, 3, 4, 5}},
  Data{ {1, 3, 4, 5, 6}, 2, {1, 2, 3, 4, 5, 6}},
  Data{ {1, 2, 4, 5, 6}, 3, {1, 2, 3, 4, 5, 6}},
  Data{ {1, 2, 3, 5, 6}, 4, {1, 2, 3, 4, 5, 6}},
  Data{ {1, 2, 3, 4, 6}, 5, {1, 2, 3, 4, 5, 6}},

  // Uses binary search.
  Data{ {1, 3, 4, 5, 6, 7}, 2, {1, 2, 3, 4, 5, 6, 7}},
  Data{ {1, 2, 4, 5, 6, 7}, 3, {1, 2, 3, 4, 5, 6, 7}},
  Data{ {1, 2, 3, 5, 6, 7}, 4, {1, 2, 3, 4, 5, 6, 7}},
  Data{ {1, 2, 3, 4, 6, 7}, 5, {1, 2, 3, 4, 5, 6, 7}},
  Data{ {1, 2, 3, 4, 5, 7}, 6, {1, 2, 3, 4, 5, 6, 7}},

  Data{ {1, 3, 4, 5, 6, 7, 8}, 2, {1, 2, 3, 4, 5, 6, 7, 8}},
  Data{ {1, 2, 4, 5, 6, 7, 8}, 3, {1, 2, 3, 4, 5, 6, 7, 8}},
  Data{ {1, 2, 3, 5, 6, 7, 8}, 4, {1, 2, 3, 4, 5, 6, 7, 8}},
  Data{ {1, 2, 3, 4, 6, 7, 8}, 5, {1, 2, 3, 4, 5, 6, 7, 8}},
  Data{ {1, 2, 3, 4, 5, 7, 8}, 6, {1, 2, 3, 4, 5, 6, 7, 8}},
  Data{ {1, 2, 3, 4, 5, 6, 8}, 7, {1, 2, 3, 4, 5, 6, 7, 8}}
));

}  // namespace
