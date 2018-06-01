#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

#include "Large_File.h"

void
IW_Large_File::_default_values ()
{
  _fd = -1;

  _pstart = 0;

  _pend = 0;

  _records_read = 0;

  _ok = 1;

  return;
}

IW_Large_File::IW_Large_File ()
{
  _default_values ();

  return;
}

IW_Large_File::IW_Large_File (const char * fname)
{
  _default_values ();

  (void) open (fname);

  return;
}

IW_Large_File::~IW_Large_File ()
{
  if (_fd >= 0)
    ::close (_fd);

  return;
}

int 
IW_Large_File::open (const char * fname)
{
  if (is_open ())
  {
    cerr << "IW_Large_File::open: already open, fd = " << _fd << endl;
    return 0;
  }

  _fd = ::open64 (fname, O_RDONLY);

  if (_fd >= 0)
  {
    _ok = 1;
    return 1;
  }

  cerr << "IW_Large_File::open: cannot open '" << fname << "'\n";

  _ok = 0;

  return 0;
}

int
IW_Large_File::_index_of_newline (char newline) const
{
  for (int i = _pstart; i < _pend; i++)
  {
    if (newline == _pending[i])
      return i;
  }

  return -1;
}


int
IW_Large_File::next_record (IWString & buffer, char newline)
{
  if (! is_open ())
  {
    cerr << "IW_Large_File::next_record: file not open\n";
    return 0;
  }

  buffer.resize_keep_storage (0);

  while (1)
  {
    cerr << "start = " << _pstart << " end is " << _pend << endl;
     
    if (_pstart < _pend)    // there are some characters to process
    {
      int p = _index_of_newline (newline);
      cerr << "newline found at " << p << endl;

      if (p >= 0)     // great, found a newline in _pending
      {
        buffer.strncpy (_pending + _pstart, p + 1);

        _pstart = p + 1;

        _records_read++;

        return 1;
      }

//    no newline in _pending. Copy everything left to buffer

      buffer.strncat (_pending + _pstart, _pend - _pstart);
    }

    _pstart = 0;
    _pend = read (_fd, _pending, sizeof (_pending));
    cerr << "Read " << _pend << " bytes\n";

    if (_pend < 0)
    {
      cerr << "IW_Large_File::next_record: error reading file\n";
      return 0;
    }

    if (0 == _pend)
    {
      cerr << "IW_Large_File::next_record: missing EOF\n";
      return 0;
    }
  }

  return 1;
}
