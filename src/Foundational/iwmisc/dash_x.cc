#include <stdlib.h>

#include "misc.h"

int
dash_x (const char * fname)
{
  struct stat stbuf;

  int rc = stat (fname, &stbuf);

  if (0 != rc)    // file not present, definitely not executable
    return 0;

  if (0 == stbuf.st_size)    // an executable file must have finite size
    return 0;

  if (0 != ((stbuf.st_mode) & S_IEXEC))
    return 1;

  return 0;
}
