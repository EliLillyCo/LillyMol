#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmin<float>;

#else

static void
unused ()
{
  iwmin<float> foo (9.0);
}

#endif
