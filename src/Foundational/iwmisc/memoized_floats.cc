#include "Foundational/iwmisc/memoized_floats.h"

namespace memoized_floats {
MemoizedFloats::MemoizedFloats() {
}

int
MemoizedFloats::Build(int max_int, int  float_digits) {
  _positive_ints.initialise(max_int);
  _positive_ints.append_to_each_stored_string(".");

  _negative_ints.initialise(max_int);
  _negative_ints.set_leading_string("-");
  _negative_ints.append_to_each_stored_string(".");

  _fraction_max_value = 1;  // Assume at least one digit of precision.
  for (int i = 0; i < float_digits; ++i) {
    _fraction_max_value *= 10;
  }

  _fractions.initialise(_fraction_max_value);

  return 1;
}

IWString
MemoizedFloats::Representation(float value) const {
  IWString to_be_returned;

  float fraction;

  if (value < 0.0f) {
    int int_part = static_cast<int>(-value);
    _negative_ints.append_number(to_be_returned, int_part);
    fraction = -value - int_part;
  } else {
    int int_part = static_cast<int>(value);
    _positive_ints.append_number(to_be_returned, int_part);
    fraction = value - int_part;
  }

  // fraction now holds a value in [0,1). Convert to an int.

  int ndx = static_cast<int>(fraction * _fraction_max_value + 0.49999f);
  _fractions.append_number(to_be_returned, ndx);

  return to_be_returned;
}

}
