#include <stdlib.h>

#include "iwstring.h"

static int
common_remove_suffix (const char * s,
                      int & nchars,
                      char dot,
                      char directory_separator)
{
  for (int i = nchars - 1; i >= 0; i--)
  {
    if (s[i] == directory_separator)    // found directory separator before dot
      return 0;

    if (dot != s[i])
      continue;

//  Found a dot character

    nchars = i;
    return 1;
  }

  return 0;
}

int
const_IWSubstring::remove_suffix (char dot, char directory_separator)
{
  return common_remove_suffix (_data, _nchars, dot, directory_separator);
}

int
IWString::remove_suffix (char dot, char directory_separator)
{
  return common_remove_suffix (_things, _number_elements, dot, directory_separator);
}

