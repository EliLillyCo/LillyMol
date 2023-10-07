#include <stdlib.h>
#ifdef _WIN32
#else
#include <unistd.h>
#endif
#include <sys/types.h>
#include <iostream>

#include "iwstring.h"

static int
common_write(int fd, const char * s, int nchars)
{
  assert (fd >= 0);

  if (0 == nchars)
    return 1;

  int rc = IW_FD_WRITE(fd, s, nchars);

  if (rc == nchars)
    return 1;

  std::cerr << "iwstring::common_write: cannot write " << nchars << " bytes to fd " << fd << '\n';

  return 0;
}

int
const_IWSubstring::write(int fd) const
{
  return common_write(fd, _data, _nchars);
}

int
IWString::write(int fd) const
{
  return common_write(fd, _things, _number_elements);
}
