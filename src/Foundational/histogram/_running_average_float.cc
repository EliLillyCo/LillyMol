#include <stdlib.h>

#define RUNNING_AVERAGE_IMPLEMENTATION
#include "running_average.h"

#if defined (__GNUG__) || defined (__SUNPRO_CC)

template class Running_Average<float>;

#else

static void
unused ()
{
  Running_Average<float> foo;
}

#endif
