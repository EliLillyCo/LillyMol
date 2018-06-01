#include <stdlib.h>

#include "new_array_.h"

double *
new_double (int size, double initial_value)
{
  assert (size > 0);

  return iw_new_array (size, initial_value);
}

template double * iw_new_array (int, double);
