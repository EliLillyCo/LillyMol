#include <stdlib.h>

#define IWCREX_IMPLEMENTATION
#include "iwcrex.h"
#include "iwlregex.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class IW_Regular_Expression_Template<IW_lregex>;

#else

static void
unused ()
{
  IW_Regular_Expression_Template<IW_lregex> foo1;

  const_IWSubstring s2;
  IW_Regular_Expression_Template<IW_lregex> foo2 (s2);

  IW_Regular_Expression_Template<IW_lregex> foo3 ("");

  IWString iwstring;
  foo1.set_pattern (iwstring);
  (void) foo2.matches (iwstring);
  (void) foo2.matches ("");
  (void) foo2.matches ("", 0);

  IW_Regular_Expression_Template<IW_lregex> foo4 (iwstring);
}

#endif
