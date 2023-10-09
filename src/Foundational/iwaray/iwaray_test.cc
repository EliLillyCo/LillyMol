// Tests for iwaray.

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "iwaray.h"

namespace {

using testing::ElementsAreArray;

TEST(iwaray, test_reserve) {
  constexpr int n = 10;
  resizable_array<float> v(n);
  EXPECT_EQ(v.size(), 0);
  for (int i = 0; i < n; ++i) {
    v.add(static_cast<float>(i));
  }
  EXPECT_EQ(v.size(), n);
  v.reserve(n);
  EXPECT_EQ(v.size(), n);
  EXPECT_EQ(v.capacity(), n);
  v.reserve(n / 2);
  EXPECT_EQ(v.size(), n);
  EXPECT_EQ(v.capacity(), n);
}

TEST(iwaray, test_make_room) {
  constexpr int n = 10;
  resizable_array<int> v(n);
  EXPECT_EQ(v.size(), 0);
  for (int i = 0; i < n; ++i) {
    v.add(i);
  }
  EXPECT_EQ(v.size(), n);
  v.make_room_for_extra_items(1);
  EXPECT_GE(v.capacity(), n + 1);
  v.add(0);
  EXPECT_GE(v.capacity(), n + 1);
}

class ForTesting {
  private:
    int _value;
  public:
    ForTesting() {
      _value = {};
    }
    int get_value() const { return _value;}
    void set_value(int s) { _value = s;}
};

template <typename T>
resizable_array<T> ArrayOfScalars(const int n) {
  resizable_array<T> to_be_returned(n);
  for (int i = 0; i < n; ++i) {
    to_be_returned.add(static_cast<T>(i));
  }

  return to_be_returned;
}

template <typename T>
resizable_array_p<T> ArrayOfObjects(const int n) {
  resizable_array_p<T> to_be_returned(n);
  for (int i = 0; i < n; ++i) {
    to_be_returned.add(new T());
  }

  return to_be_returned;
}

TEST(iwaray, test_swap) {
  constexpr int n = 17;
  resizable_array<int> v = ArrayOfScalars<int>(n);
  EXPECT_EQ(v[2], 2);
  EXPECT_EQ(v[3], 3);
  v.swap_elements(2, 3);
  EXPECT_EQ(v[2], 3);
  EXPECT_EQ(v[3], 2);
}

TEST(iwaray, test_add_non_duplicated_elements) {
  constexpr int n = 11;
  resizable_array<int> v1 = ArrayOfScalars<int>(n);
  resizable_array<int> v2 = v1;
  EXPECT_EQ(v1.add_non_duplicated_elements(v2), 0);
  v2.resize(n / 2);
  EXPECT_EQ(v1.add_non_duplicated_elements(v2), 0);
  for (int& i : v2) {
    i = n + i;
  }
  EXPECT_EQ(v1.add_non_duplicated_elements(v2), v2.size());
  v1.resize(n);
  EXPECT_EQ(v1.size(), n);

  for (int& i : v2) {
    i = n;
  }
  EXPECT_EQ(v1.add_non_duplicated_elements(v2), 1);
}

TEST(iwaray, test_add_if_not_already_present) {
  constexpr int n = 8;
  resizable_array<int> v = ArrayOfScalars<int>(n);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v.add_if_not_already_present(i), 0);
  }
  EXPECT_GT(v.add_if_not_already_present(n), 0);
}

TEST(iwaray, Pointers) {
  constexpr int n = 50;
  resizable_array_p<ForTesting> v;

  for (int i = 0;i < n; ++i) {
    ForTesting* f = new ForTesting();
    f->set_value(i);
    v.add(f);
  }
  EXPECT_EQ(v.number_elements(), n);
  EXPECT_GE(v.capacity(), n);

  v.resize(n / 2);
  EXPECT_EQ(v.number_elements(), n / 2);
}

TEST(iwaray, test_remove_no_delete) {
  constexpr int n = 10;
  resizable_array_p<ForTesting> v = ArrayOfObjects<ForTesting>(n);
  EXPECT_EQ(v.size(), n);

  ForTesting* removed = v.remove_no_delete(0);
  delete removed;
}

TEST(resizable_array, operator_ltlt) {
  resizable_array<int> x;
  for (int i = 0; i < 5; ++i) {
    x << i;
    EXPECT_EQ(x.back(), i);
    EXPECT_EQ(x.size(), i + 1);
  }

  x << 99 << 100;
  EXPECT_EQ(x.size(), 5 + 2);
  EXPECT_EQ(x.back(), 100);
}

TEST(TestResizableArray, TestInitializerList) {
  resizable_array<int> foo {3, 2, 1, 5};
  EXPECT_THAT(foo, ElementsAreArray({3, 2, 1, 5}));
}

TEST(TestResizableArray, TestIndexRelativeToPositive) {
  resizable_array<int> foo {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  const int n = foo.number_elements();
  int first_ndx = 3;
  int expected = first_ndx;
  for (int delta = 0; delta < 100; ++delta) {
    int pos = foo.index_relative_to(first_ndx, delta);
    EXPECT_EQ(pos, expected);
    expected++;
    if (expected == n) {
      expected = 0;
    }
  }
}

TEST(TestResizableArray, TestIndexRelativeToNegative) {
  resizable_array<int> foo {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  const int n = foo.number_elements();
  int first_ndx = 3;
  int expected = first_ndx;
  for (int delta = 0; delta < 100; ++delta) {
    int pos = foo.index_relative_to(first_ndx, -delta);
    EXPECT_EQ(pos, expected);
    expected--;
    if (expected < 0) {
      expected = n - 1;
    }
  }
}

TEST(TestResizableArray, TestRemoveFromToEmpty) {
  resizable_array<int> foo;
  EXPECT_LT(foo.remove_from_to(0, 1), 0);
}

TEST(TestResizableArray, TestRemoveOnlyElement) {
  resizable_array<int> foo {1};
  EXPECT_EQ(foo.remove_from_to(0, 1), 1);
  EXPECT_TRUE(foo.empty());
}

TEST(TestResizableArray, TestRemoveFirstElement) {
  resizable_array<int> foo {1, 2, 3, 4};
  EXPECT_EQ(foo.remove_from_to(0, 1), 1);
  EXPECT_THAT(foo, ElementsAreArray({2, 3, 4}));
}

TEST(TestResizableArray, TestRemoveFirstTwoElements) {
  resizable_array<int> foo {1, 2, 3, 4};
  EXPECT_EQ(foo.remove_from_to(0, 2), 2);
  EXPECT_THAT(foo, ElementsAreArray({3, 4}));
}

TEST(TestResizableArray, TestRemoveLastElement) {
  resizable_array<int> foo {1, 2, 3, 4};
  EXPECT_EQ(foo.remove_from_to(3, 4), 1);
  EXPECT_THAT(foo, ElementsAreArray({1, 2, 3}));
}

TEST(TestResizableArray, TestRemoveTwoLastElements) {
  resizable_array<int> foo {1, 2, 3, 4};
  EXPECT_EQ(foo.remove_from_to(2, 4), 2);
  EXPECT_THAT(foo, ElementsAreArray({1, 2}));
}

TEST(TestResizableArray, TestRemoveAllElements) {
  resizable_array<int> foo {1, 2, 3, 4};
  EXPECT_EQ(foo.remove_from_to(0, 4), 4);
  EXPECT_TRUE(foo.empty());
}

TEST(TestResizableArray, TestRemoveOneMiddle) {
  resizable_array<int> foo {1, 2, 3, 4};
  EXPECT_EQ(foo.remove_from_to(2, 3), 1);
  EXPECT_THAT(foo, ElementsAreArray({1, 2, 4}));
}

TEST(TestResizableArray, TestRemoveTwoMiddle) {
  resizable_array<int> foo {1, 2, 3, 4, 5};
  EXPECT_EQ(foo.remove_from_to(2, 4), 2);
  EXPECT_THAT(foo, ElementsAreArray({1, 2, 5}));
}

}  // namespace
