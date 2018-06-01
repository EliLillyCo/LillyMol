#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwminmax<int>;

#else

static void
unused ()
{
  iwminmax<int> foo (0, 0);
}

#endif
