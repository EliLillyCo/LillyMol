#include <stdlib.h>

#include "misc.h"

#include "new_array_.h"

unsigned int *
new_unsigned_int(int size, unsigned int initial_value)
{
  assert (size > 0);

  return iw_new_array(size, initial_value);
}

template unsigned int * iw_new_array(int, unsigned int);
