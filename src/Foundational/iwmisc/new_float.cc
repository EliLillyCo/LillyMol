#include <stdlib.h>

#include "new_array_.h"

float *
new_float (int size, float initial_value)
{
  assert (size > 0);

  return iw_new_array (size, initial_value);
}

template float * iw_new_array (int, float);
