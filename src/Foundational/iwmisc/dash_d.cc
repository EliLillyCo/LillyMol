#include <stdlib.h>

#include "misc.h"

int
dash_d (const char * fname)
{
  struct stat stbuf;

  int rc = stat (fname, &stbuf);

  if (0 != rc)    // file not present, definitely not executable
    return 0;

#ifdef _WIN32
  return ( _S_IFDIR & stbuf.st_mode );
#else
  return S_ISDIR (stbuf.st_mode);
#endif
}
