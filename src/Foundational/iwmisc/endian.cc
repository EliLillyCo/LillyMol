#include <stdlib.h>

static int little_endian = -1;

int
iw_little_endian ()
{
  if (little_endian >= 0)    // already been determined
    return little_endian;

  int i = 1;
  unsigned char * p = (unsigned char *) &i;

  if (0 == p[0])
    little_endian = 0;
  else
    little_endian = 1;

  return little_endian;
}
