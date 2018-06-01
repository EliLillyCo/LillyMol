#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmin<int>;

#else

static void
unused ()
{
  iwmin<int> foo (9);
}

#endif
