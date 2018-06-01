#include <stdlib.h>

#include "misc.h"

/*
  No, this is dangerous. unsigned char will overflow
*/

template unsigned int sum_vector (const unsigned char *, int);
