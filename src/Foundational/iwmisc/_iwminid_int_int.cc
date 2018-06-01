#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwminid<int, int>;

#else

static void
unused ()
{
  iwminid<int, int> foo (1, 0);
}

#endif
