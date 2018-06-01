#include <stdlib.h>
#include <assert.h>
#include <iostream>
//#include <cstring>

#include "iwzlib.h"

#ifdef NO_ZLIB
#else

#include "iwstring.h"

#define IWZLIB_BUF_SIZE 1024

IW_ZLib_Wrapper::IW_ZLib_Wrapper ()
{
  _gzfile = NULL;

 _record_delimiter = '\n';

 _buffer = NULL;

 _chars_in_buffer = 0;

 _start_of_next_record = 0;

  return;
}

IW_ZLib_Wrapper::~IW_ZLib_Wrapper ()
{
  if (NULL != _gzfile)
    gzclose (_gzfile);

  if (NULL != _buffer)
    delete _buffer;
  
  return;
}

int
IW_ZLib_Wrapper::open_file (const char * fname)
{
  if (NULL != _gzfile)
    close_file ();

#if (__GNUC__ >= 6)
  _gzfile = gzopen64 (fname, "r");
#else
  _gzfile = gzopen (fname, "r");
#endif

  if (NULL == _gzfile)
  {
    int errnum;
    cerr << "IW_ZLib_Wrapper::open_file:cannot open '" << fname << "' " << gzerror (_gzfile, &errnum) << endl;
    return 0;
  }

  _start_of_next_record = 1;

  _chars_in_buffer = 0;

  _buffer = new char[IWZLIB_BUF_SIZE];

  if (NULL == _buffer)
  {
    cerr << "IW_ZLib_Wrapper::open_file:cannot allocate buffer " << IWZLIB_BUF_SIZE << " bytes\n";
    return 0;
  }

  return 1;
}

int
IW_ZLib_Wrapper::close_file ()
{
  if (NULL == _gzfile)
  {
    cerr << "IW_ZLib_Wrapper::close_file:no file open\n";
    return 0;
  }

  gzclose (_gzfile);
  _gzfile = NULL;

  return 1;
}

//#define DEBUG_NEXT_RECORD

size_t
IW_ZLib_Wrapper::next_record (IWString & destination)
{
  assert (NULL != _gzfile);

  destination.resize_keep_storage (0);

#ifdef DEBUG_NEXT_RECORD
  cerr << "Fetching record\n";

  cerr << "_chars_in_buffer " << _chars_in_buffer << " EOF? " << gzeof(_gzfile) << endl;
#endif

  if (0 == _chars_in_buffer && gzeof (_gzfile))
    return 0;

  while (1)
  {
#ifdef DEBUG_NEXT_RECORD
#if (__GNUC__ >= 6)
    cerr << "At start of loop _chars_in_buffer = " << _chars_in_buffer << " start " << _start_of_next_record << " gztell64 " << gztell64(_gzfile) << endl;
#else
    cerr << "At start of loop _chars_in_buffer = " << _chars_in_buffer << " start " << _start_of_next_record << " gztell " << gztell(_gzfile) << endl;
#endif
#endif

    if (_start_of_next_record > _chars_in_buffer)   // read some more data
    {
      _chars_in_buffer = gzread (_gzfile, _buffer, IWZLIB_BUF_SIZE);

#ifdef DEBUG_NEXT_RECORD
      cerr << "Read " << _chars_in_buffer << " bytes\n";
#endif

      if (_chars_in_buffer < 0)
      {
        int errnum;
        cerr << "IW_ZLib_Wrapper::next_record:fatal error '" << gzerror (_gzfile, &errnum) << "'\n";
        close_file ();
        return 0;
      }

      if (0 == _chars_in_buffer)
      {
        if (destination.length () > 0)
        {
          cerr << "IW_ZLib_Wrapper::next_record:no record terminator\n";
          return 1;
        }

        return 0;    // normal eof
      }

      _start_of_next_record = 0;
    }

//  Now that our buffer has some data, copy it to DESTINATION

    const void * c = ::memchr (static_cast<const void *> (_buffer + _start_of_next_record), _record_delimiter, _chars_in_buffer - _start_of_next_record);

    if (NULL == c)    // no record delimiter found, fetch another record
    {
      destination.strncat (_buffer + _start_of_next_record, _chars_in_buffer - _start_of_next_record);
      _start_of_next_record = _chars_in_buffer + 1;    // force reading more data
      continue;
    }

    int chars_to_copy = static_cast<const char *> (c) - (_buffer + _start_of_next_record);

    if (0 == chars_to_copy)     // blank line
    {
      _start_of_next_record++;
      return 1;
    }

#ifdef DEBUG_NEXT_RECORD
    cerr << "Found " << chars_to_copy << " characters to copy\n";
#endif

    destination.strncat (_buffer + _start_of_next_record, chars_to_copy);
    _start_of_next_record += chars_to_copy + 1;

    return 1;
  }
}

int
IW_ZLib_Wrapper::seekg (z_off_t o)
{
  assert (NULL != _gzfile);

//cerr << "Seeking to " << o << ", eof? " << gzeof(_gzfile) << endl;

  if (gzeof(_gzfile))
    gzclearerr(_gzfile);

#if (__GNUC__ >= 6)
  if (gzseek64(_gzfile, o, 0) >= 0)
#else
  if (gzseek(_gzfile, o, 0) >= 0)
#endif
  {
    _chars_in_buffer = 0;
    return 1;
  }

  int errnum;

  cerr << "IW_ZLib_Wrapper::seekg:cannot seek '" << gzerror (_gzfile, &errnum) << "'\n";
  return 0;
}

z_off_t
IW_ZLib_Wrapper::tellg () const
{
  assert (NULL != _gzfile);

#if (__GNUC__ >= 6)
  z_off_t o = gztell64(_gzfile);
#else
  z_off_t o = gztell(_gzfile);
#endif

//cerr << "gztell yields " << o << ", buffer contains " << _chars_in_buffer << " start " << _start_of_next_record << endl;
//cerr << "gztell yields " << static_cast<unsigned int>(o) << ", buffer contains " << _chars_in_buffer << " start " << _start_of_next_record << endl;

// Dec 2006. I saw very wierd things. gztell was yielding huge numbers,
// so, I just cast it back to something reasonable. This is obviously
// quite wrong, but it seems to work. But seeking in a gzip'd file isn't
// a good idea anyway...

  o = static_cast<unsigned int>(o);

  o -= _chars_in_buffer;

  o += _start_of_next_record;

  return o;
}

int
IW_ZLib_Wrapper::eof () const
{
  assert (NULL != _gzfile);

  return gzeof (_gzfile);
}

int
IW_ZLib_Wrapper::read_bytes (void * destination, int nbytes)
{
  assert (NULL != _gzfile);

  if (gzeof (_gzfile))
    return 0;

  int rc = gzread (_gzfile, destination, nbytes);

  if (rc >= 0)
    return rc;

  int errnum;
  cerr << "IW_ZLib_Wrapper::read_bytes:fatal error '" << gzerror (_gzfile, &errnum) << "'\n";

  return 0;
}

#endif
