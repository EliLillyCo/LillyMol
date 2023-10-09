#include <stdlib.h>
#ifdef _WIN32
#else
#include <unistd.h>
#endif
#include <sys/types.h>
#include <iostream>

#include "Foundational/iwmisc/iwconfig.h"
#include "iwstring.h"

//#define DEBUG_GETLINE

int
IWString::getline(std::istream & is, char record_delimiter)
{
  if (! is.good())
  {
    if (is.eof())
      return 0;
    std::cerr << "IWString::getline: input stream not good\n";
    return 0;
  }

  if (is.eof())
    return 0;

  if (0 == _elements_allocated)
    resize(256);

  _number_elements = 0;

  while (1)
  {
    int bytes_available = _elements_allocated - _number_elements;
    is.get (&(_things[_number_elements]), bytes_available, record_delimiter);
    int nchars = is.gcount();
    _number_elements += nchars;

#ifdef DEBUG_GETLINE
    std::cerr << "Just read " << nchars << " characters from stream\n";
    std::cerr << "n = " << _number_elements << " alloc = " << _elements_allocated << '\n';
#endif

    if (is.eof())
    {
      if (_number_elements)
        std::cerr << "Incomplete last record, have " << _number_elements << " chars\n";

      return _number_elements;
    }

//  Important change gcc-3.0.1 
//  When a blank line is encountered, is.good() is set to bad....

    if (0 == nchars && ! is.good())
    {
      is.clear();
      if (record_delimiter == is.get())
        return 1;
      
      return 0;
    }

    if (_number_elements < _elements_allocated - 1)   // whole record fit into buffer
    {
      (void) is.get();     // fetch the next character just to check for EOF
      if (is.eof())
      {
        std::cerr << "IWString::getline: record not terminated\n";
        return _number_elements;
      }

      return _number_elements;
    }

    resize(_elements_allocated * 2);
//  std::cerr << "Buffer expanded to " << _elements_allocated << '\n';
  }
}

int
operator >> (std::istream & input, IWString & rhs)
{
  return rhs.getline(input);
}

//#ifndef linux

static int
_bytes_to_char(int fd, char zchar)
{
  char rdbuf[4096];

  ssize_t chars_read;

  int rc = 0;

  while ( (chars_read = IW_FD_READ(fd, rdbuf, sizeof(rdbuf))) > 0)
  {
    char * c = static_cast<char *>(memchr(rdbuf, zchar, chars_read));

    if (nullptr != c)
      return static_cast<int>(rc + (c - rdbuf + 1));   // we include ZCHAR itself

     rc += chars_read;
  }

  return rc;
}

static int
bytes_to_char(int fd, char zchar)
{
//off_t current_position = tell(fd);
  off_t current_offset = IW_FD_LSEEK(fd, 0, SEEK_CUR);

  int rc = _bytes_to_char(fd, zchar);

  IW_FD_LSEEK(fd, current_offset, 0);

//std::cerr << "At " << current_offset << " need " << rc << " bytes to newline\n";

  return rc;
}

int
IWString::getline(int fd, char newline)
{
  resize_keep_storage(0);

  int bytes_to_newline = bytes_to_char(fd, newline);

  if (0 == bytes_to_newline)
    return 0;

  if (_elements_allocated < bytes_to_newline)
  {
    if (! resize(bytes_to_newline))
    {
      std::cerr << "IWString::getline: cannot allocate " << bytes_to_newline << " bytes\n";
      return 0;
    }
  }

  ssize_t chars_read = IW_FD_READ(fd, _things, bytes_to_newline);

  if (chars_read <= 0)
  {
    std::cerr << "IWString::getline: cannot read " << bytes_to_newline << " bytes from fd\n";
    return 0;
  }

  if (newline != _things[bytes_to_newline - 1])
  {
    std::cerr << "IWString::getline_fd: record not terminated\n";
    _number_elements = bytes_to_newline;
  }
  else
    _number_elements = bytes_to_newline - 1;    // we don't want the newline separator in our buffer

  return 1;
}

//#endif
