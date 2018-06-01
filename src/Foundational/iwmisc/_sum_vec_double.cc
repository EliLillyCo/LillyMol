#include "misc.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template double sum_vector (const double *, int);

#else

static void
unused ()
{
  double x[2];
  (void) sum_vector (x, 0);
}

#endif
