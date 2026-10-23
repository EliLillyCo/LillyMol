#include <cstdint>
#include <limits>
#include <string>
#include <tuple>

#include "gtest/gtest.h"

#include "string_with_commas.h"

class WithCommasTest
    : public ::testing::TestWithParam<std::tuple<long long, std::string>> {};

TEST_P(WithCommasTest, FormatsIntegerWithCommas) {
    const auto& [input, expected] = GetParam();
    EXPECT_EQ(iwstring::with_commas(input), expected);
}

INSTANTIATE_TEST_SUITE_P(
    IntegerFormatting,
    WithCommasTest,
    ::testing::Values(
        std::make_tuple(-1LL, "-1"),
        std::make_tuple(0LL, "0"),
        std::make_tuple(1LL, "1"),

        std::make_tuple(9LL, "9"),
        std::make_tuple(10LL, "10"),
        std::make_tuple(99LL, "99"),
        std::make_tuple(100LL, "100"),

        std::make_tuple(999LL, "999"),
        std::make_tuple(1000LL, "1,000"),
        std::make_tuple(1001LL, "1,001"),
        std::make_tuple(9999LL, "9,999"),
        std::make_tuple(10000LL, "10,000"),
        std::make_tuple(99999LL, "99,999"),
        std::make_tuple(100000LL, "100,000"),
        std::make_tuple(999999LL, "999,999"),
        std::make_tuple(1000000LL, "1,000,000"),

        std::make_tuple(-999LL, "-999"),
        std::make_tuple(-1000LL, "-1,000"),
        std::make_tuple(-1001LL, "-1,001"),
        std::make_tuple(-999999LL, "-999,999"),
        std::make_tuple(-1000000LL, "-1,000,000"),

        std::make_tuple(1234567890LL, "1,234,567,890"),
        std::make_tuple(-1234567890LL, "-1,234,567,890"),

        std::make_tuple(
            std::numeric_limits<long long>::max(),
            "9,223,372,036,854,775,807"
        )

        // Do not test std::numeric_limits<long long>::min()
        // with the original implementation if it does `value < 0`
        // and then relies on std::to_string(value), unless the code
        // avoids negating the value. The version I gave earlier is okay
        // because it does not negate, but keep this note if you later
        // change the implementation.
    )
);
