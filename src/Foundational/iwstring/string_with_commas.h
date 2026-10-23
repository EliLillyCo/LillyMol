#ifndef FOUNDATIONAL_IWSTRING_STRING_WITH_COMMAS_H_
#define FOUNDATIONAL_IWSTRING_STRING_WITH_COMMAS_H_

#include <string>
#include <type_traits>

namespace iwstring {

// Returns the decimal representation of `value` with commas separating
// groups of three digits.
template <typename Int>
std::string
with_commas(Int value) {
  static_assert(std::is_integral_v<Int>, "with_commas requires an integer type");

  std::string result = std::to_string(value);
  const std::size_t first_digit =
      !result.empty() && result[0] == '-' ? 1 : 0;

  for (std::size_t position = result.size();
       position > first_digit + 3; position -= 3) {
    result.insert(position - 3, 1, ',');
  }

  return result;
}

}  // namespace iwstring

#endif  // FOUNDATIONAL_IWSTRING_STRING_WITH_COMMAS_H_
