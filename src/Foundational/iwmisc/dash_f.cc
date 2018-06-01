#include <stdlib.h>

#include "misc.h"

int
dash_f (const char * fname)
{
  struct stat stbuf;

  if (0 == stat (fname, &stbuf))
    return 1;

  return 0;
}
