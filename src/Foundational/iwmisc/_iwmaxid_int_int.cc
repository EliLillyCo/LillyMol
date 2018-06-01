#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmaxid<int, int>;

#else

static void
unused ()
{
  iwmaxid<int, int> foo (0, 0);
}

#endif
