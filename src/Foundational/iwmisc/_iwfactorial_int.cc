#include <stdlib.h>

#define IWFACTORIAL_IMPLEMENTATION
#include "iwfactorial.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class IW_Factorial<int>;

#else

static void
unused ()
{
  IW_Factorial<int> foo;
}

#endif
