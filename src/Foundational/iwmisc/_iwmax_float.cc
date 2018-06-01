#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmax<float>;

#else

static void
unused ()
{
  iwmax<float> f(0.0);
}

#endif
