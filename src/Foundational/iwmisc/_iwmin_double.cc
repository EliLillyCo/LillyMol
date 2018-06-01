#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmin<double>;

#else

static void
unused ()
{
  iwmin<double> foo (9.0);
}

#endif
