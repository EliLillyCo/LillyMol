#include "misc.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template float sum_vector (const float *, int);

#else

static void
unused ()
{
  double x[2];
  (void) sum_vector (x, 0);
}

#endif
