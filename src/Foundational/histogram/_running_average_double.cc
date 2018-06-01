#include <stdlib.h>

#define RUNNING_AVERAGE_IMPLEMENTATION
#include "running_average.h"

#ifdef __GNUG__

template class Running_Average<double>;

#else

static void
unused ()
{
  Running_Average<double> foo;
}

#endif
