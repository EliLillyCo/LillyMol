#include <stdlib.h>

#ifdef TEST_RXPOSIX
#include "rxposix.h"

#include "iwstring.h"

int
IWString::matches(regex_t & rx) const
{
  return 0 == regnexec(&rx, _things, _number_elements, 0, NULL, 0);
}

int
const_IWSubstring::matches(regex_t & rx) const
{
  return 0 == regnexec(&rx, _data, _nchars, 0, NULL, 0);
}

#endif
