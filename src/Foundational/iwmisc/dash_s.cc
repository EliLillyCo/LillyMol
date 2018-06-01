#include <stdlib.h>

#include "misc.h"

off_t
dash_s (const char * fname)
{
  struct stat stbuf;

  int rc = stat (fname, &stbuf);

  if (0 != rc)    // file not present, definitely not executable
    return 0;

  return stbuf.st_size;
}
