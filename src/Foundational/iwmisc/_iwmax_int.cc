#define IWMINMAX_IMPLEMENTATION
#include "iwminmax.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwmax<int>;

#else

static void
unused ()
{
  iwmax<int> f(0);
}
#endif
