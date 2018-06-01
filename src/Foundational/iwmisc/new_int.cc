#include <stdlib.h>

#include "new_array_.h"

int *
new_int (int size, int initial_value)
{
  assert (size > 0);

  return iw_new_array (size, initial_value);
}

template int * iw_new_array (int, int);
