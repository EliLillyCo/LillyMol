#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmax<double>;

#else

static void
unused ()
{
  iwmax<double> f(0.0);
}

#endif

