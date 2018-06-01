#include <stdlib.h>

#include "iwstring.h"

int
IWString::numeric_value (long long & result) const
{
  return string_class_is_int(_things, _number_elements, result);
}

int
const_IWSubstring::numeric_value (long long & result) const
{
  return string_class_is_int(_data, _nchars, result);
}
