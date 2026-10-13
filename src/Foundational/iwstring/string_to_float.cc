// Those IWString member functions associated with string to float conversions

// #define USE_DRAGONBOX
#ifdef USE_DRAGONBOX
#include "dragonbox/dragonbox_to_chars.h"
#endif

#include "iwstring.h"
#include <vector>
#include <limits>

namespace {
// Common formats stay on the stack; allocate only when the output exceeds it.
void AppendFormatted(IWString& output, double value, const char* format, int precision = -1) {
  char buffer[64];
  int length = precision < 0 ? snprintf(buffer, sizeof(buffer), format, value)
                            : snprintf(buffer, sizeof(buffer), format, precision, value);
  if (length < 0) {
    return;
  }
  if (length < static_cast<int>(sizeof(buffer))) {
    output.strncat(buffer, length);
    return;
  }
  std::vector<char> expanded(static_cast<size_t>(length) + 1);
  if (precision < 0) {
    snprintf(expanded.data(), expanded.size(), format, value);
  } else {
    snprintf(expanded.data(), expanded.size(), format, precision, value);
  }
  output.strncat(expanded.data(), length);
}
}  // namespace

static int float_precision = 7;

void
set_default_iwstring_float_concatenation_precision (int s)
{
  float_precision = s;
}

static int double_precision = 10;

void
set_default_iwstring_double_concatenation_precision (int s)
{
  double_precision = s;
}

void
IWString::append_number(float f, int fprecision)
{
  // Preserve gcvt output on the usual path, with a safe significant-digit bound.
  if (fprecision > 0 && fprecision <= std::numeric_limits<double>::max_digits10) {
    char buffer[32];
    const char* formatted = gcvt(static_cast<double>(f), fprecision, buffer);
    if (formatted == nullptr) return;
    strncat(buffer, static_cast<int>(::strlen(buffer)));
  } else {
    AppendFormatted(*this, f, "%.*g", std::max(1, fprecision));
  }

}

void
IWString::append_number(float f) {
  append_number(f, float_precision);
}

void
IWString::append_number(double d)
{
  append_number(d, double_precision);
}

void
IWString::append_number(double d, int dprecision)
{
  // Preserve gcvt output on the usual path, with a safe significant-digit bound.
  if (dprecision > 0 && dprecision <= std::numeric_limits<double>::max_digits10) {
    char buffer[32];
    const char* formatted = gcvt(static_cast<double>(d), dprecision, buffer);
    if (formatted == nullptr) return;
    strncat(buffer, static_cast<int>(::strlen(buffer)));
  } else {
    AppendFormatted(*this, d, "%.*g", std::max(1, dprecision));
  }

}

void
IWString::append_number(float f, const char * fformat)
{
  AppendFormatted(*this, f, fformat);

}

#ifdef USE_DRAGONBOX
void
IWString::append_number_dragonbox(float f) {
  if (f == 0.0f) {
    resizable_array<char>::add('0');
    return;
  }

  constexpr int kBufferLength = 1 + // for '\0'
    jkj::dragonbox::max_output_string_length<jkj::dragonbox::ieee754_binary64>;

  this->make_room_for_extra_items(kBufferLength);

  int initial_number_elements = _number_elements;
  char* end_ptr = jkj::dragonbox::to_chars_n(f, _things + _number_elements);

  _number_elements = end_ptr - _things;

  for (int i = _number_elements - 2; i > initial_number_elements; --i) {
    if (_things[i] == 'E') {
      _things[i] = 'e';
      break;
    }
  }

}
#endif
