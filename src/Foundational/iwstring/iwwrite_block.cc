#include <stdlib.h>
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

  int chars_written = IW_FD_WRITE (fd, _things, blocks_to_write * IWWRITE_BLKSIZE);

  if (chars_written != blocks_to_write * IWWRITE_BLKSIZE)
  {
    cerr << "IWString::write_whole_blocks_shift_unwritten:cannot write " << blocks_to_write << " blocks to " << fd << endl;
    return 0;
  }

  ::memcpy (_things, _things + chars_written, _number_elements - chars_written);

  _number_elements -= chars_written;

  return chars_written;
}
