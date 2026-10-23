#include <stdlib.h>
#include <cerrno>
#include <sys/types.h>
#ifdef _WIN32
#else
#include <unistd.h>
#endif

#include <iostream>

using std::cerr;
using std::endl;

#include "Foundational/iwmisc/iwconfig.h"
#include "iwstring.h"

#define IWWRITE_BLKSIZE 4096

int
IWString::write_whole_blocks_shift_unwritten (int fd)
{
  if (_number_elements < IWWRITE_BLKSIZE)
    return 1;


  int blocks_to_write = _number_elements / IWWRITE_BLKSIZE;

  const int requested = blocks_to_write * IWWRITE_BLKSIZE;
  int chars_written = 0;
  while (chars_written < requested) {
    const ssize_t rc = IW_FD_WRITE(fd, _things + chars_written, requested - chars_written);
    if (rc < 0 && errno == EINTR) {
      continue;
    }
    if (rc <= 0) {
      if (chars_written > 0) {
        ::memmove(_things, _things + chars_written, _number_elements - chars_written);
        _number_elements -= chars_written;
      }
      return 0;
    }
    chars_written += rc;
  }
  ::memmove(_things, _things + chars_written, _number_elements - chars_written);
  _number_elements -= chars_written;

  return chars_written;
}
