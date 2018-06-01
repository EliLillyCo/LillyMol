#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwminmax<double>;

#else

static void
unused ()
{
  iwminmax<double> foo (0.0, 8.0);
}

#endif
