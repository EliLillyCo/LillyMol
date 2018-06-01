#include <stdlib.h>
#include "iwconfig.h"
#ifdef _WIN32
#else
#include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "zlib.h"

#include "iwstring.h"

IWString_and_File_Descriptor::IWString_and_File_Descriptor()
{
  _fd = -9;

  _gzfile = NULL;

  return;
}

IWString_and_File_Descriptor::IWString_and_File_Descriptor(int f) : _fd (f)
{
  _gzfile = NULL;

  return;
}

IWString_and_File_Descriptor::~IWString_and_File_Descriptor()
{
  if (_gzfile)
  {
    _compress_and_write();
    gzclose(_gzfile);
    _gzfile = NULL;
  }
  else if (_fd < 0)
    ;
  else if (_fd ==1)
  {
	IWString::write(_fd);
    return;
  }
  else if (0 == IWString::length())
  {
    IW_FD_CLOSE(_fd);
    _fd = -9;
  }
  else
  {
    IWString::write(_fd);
    IW_FD_CLOSE(_fd);
    _fd = -9;
  }

  return;
}

int
IWString_and_File_Descriptor::open(const char * fname)
{
  if (_fd > 0)
  {
    cerr << "IWString_and_File_Descriptor::open:already open\n";
    return 0;
  }

// I haven't allowed for '>>foo.gz'. Tough

  int strlen_fname = strlen(fname);

  if (1 == strlen_fname && '-' == fname[0])
  {
    _fd = 1;
    return 1;
  }

  if (strlen_fname <= 3)   // cannot be of the form 'xxx.gz'
    ;
  else if (0 == ::strcmp(fname + strlen_fname - 3, ".gz"))
    return _open_gzipd_stream(fname);

  int mode;

  if (strlen_fname > 2 && 0 == ::strcmp(fname, ">>"))
  {
    fname += 2;
    mode = O_WRONLY | O_APPEND | O_CREAT;
  }
  else
    mode = O_WRONLY | O_TRUNC | O_CREAT;

#ifdef _WIN32
  int flags = 0;
#else
  int flags = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
#endif

  _fd = IW_FD_OPEN(fname, mode, flags);

  if (_fd < 0)
  {
    cerr << "IWString_and_File_Descriptor::open:cannot open '" << fname << "'\n";
    return 0;
  }

  return _fd;
}

int
IWString_and_File_Descriptor::_open_gzipd_stream(const char * fname)
{
#if (__GNUC__ >= 6)
  _gzfile = gzopen64(fname, "wb");
#else
  _gzfile = gzopen(fname, "wb");
#endif

  if (NULL != _gzfile)
    return 1;

  cerr << "IWString_and_File_Descriptor::_open_gzipd_stream:cannot open '" << fname << "'\n";
  return 0;
}

int
IWString_and_File_Descriptor::_do_resize(int keep_storage)
{
  if (IWSFD_KEEP_STORAGE_ON_FLUSH == keep_storage)
    IWString::resize_keep_storage(0);
  else if (0 == keep_storage)
    IWString::resize(0);
  else
  {
    cerr << "IWString_and_File_Descriptor::_do_resize:unrecognised option " << keep_storage << '\n';
    return 0;
  }

  return 1;
}

int
IWString_and_File_Descriptor::flush(int keep_storage)
{
  if (0 == IWString::length())
    return 1;

  if (_gzfile)
    _compress_and_write();
  else if (_fd < 0)
  {
    cerr << "IWString_and_File_Descriptor::flush:no file descriptor\n";
    return 0;
  }
  else
    IWString::write(_fd);

  return _do_resize(IWSFD_KEEP_STORAGE_ON_FLUSH);
}

int
IWString_and_File_Descriptor::write_whole_blocks_shift_unwritten()
{
  if (0 == IWString::length())
    return 0;

  if (_gzfile)
  {
    if (! _compress_and_write())
      return 0;

    return _do_resize(IWSFD_KEEP_STORAGE_ON_FLUSH);
  }

  if (_fd < 0)
  {
    cerr << "IWString_and_File_Descriptor::write_whole_blocks_shift_unwritten:invalid file descriptor\n";
    return 0;
  }

  return IWString::write_whole_blocks_shift_unwritten(_fd);
}

int
IWString_and_File_Descriptor::close()
{
  if (_gzfile)
  {
    _compress_and_write();
    gzclose (_gzfile);
    _gzfile = NULL;
    _do_resize(0);
    return 1;
  }

  if (_fd < 0)
  {
    cerr << "IWString_and_File_Descriptor::close:not open\n";
    return 0;
  }

  if (_number_elements)
    flush();

  IW_FD_CLOSE(_fd);

  _fd = -1;

  return 1;
}

int
IWString_and_File_Descriptor::_compress_and_write()
{
  if (0 == _number_elements)
    return 1;

  int bytes_written = gzwrite(_gzfile, _things, static_cast<unsigned>(_number_elements));

  if (bytes_written > 0)
    return bytes_written;

  cerr << "IWString_and_File_Descriptor::_compress_and_write:cannot write\n";
  return 0;
}
ssize_t
IWString_and_File_Descriptor::write (const char * s, size_t nchars)
{
  if (nchars + _number_elements < 32768)
    return IWString::strncat(s, nchars);

  if (! is_open())
  {
    cerr << "IWString_and_File_Descriptor::write:not open\n";
    return 0;
  }

  if (! flush(1))    // 1 means keep storage
    return 0;

  ssize_t rc = 0;

  while (nchars)
  {
    if (nchars < 32768)
      return IWString::strncat(s, nchars);

    int chars_to_write;
    if (nchars >= 32768)
      chars_to_write = 32768;
    else
      chars_to_write = nchars;

    ssize_t bytes_written = IW_FD_WRITE(_fd, s, chars_to_write);
    if (bytes_written != static_cast<ssize_t>(chars_to_write))
    {
      cerr << "IWString_and_File_Descriptor::write:cannot write " << chars_to_write << " to " << _fd << '\n';
      return 0;
    }

    rc += bytes_written;

    s += chars_to_write;
    nchars -= chars_to_write;
  }

  return rc;
}
