#include "iwconfig.h"

#include <stdlib.h>
#include <assert.h>
#include <limits>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef _WIN32
#else
#include <unistd.h>
#endif

#include <iostream>

#include "iwstring_data_source.h"

static int
iw_open_file (const char * fname)
{
#ifdef IW_SUN
        int rc = open64 (fname, O_RDONLY);
#else
        int rc = IW_FD_OPEN (fname, O_RDONLY);
#endif

        return rc;
}

void
iwstring_data_source::_default_values (int lrecl)
{
        _isstringbuffer = false;

        _fd = -1;
        _record_buffered = 0;
        _longest_record = 0;
        _record_delimiter = '\n';
        _lines_read = 0;
        _lines_which_are_returned = 0;
        _strip_leading_blanks = 0;
        _strip_trailing_blanks = 0;
        _compress_spaces  = 0;
        _skip_blank_lines = 0;
        _convert_to_lowercase = 0;
        _convert_to_uppercase = 0;
        _open = 0;
        _good = 1;
        _eof = 0;
        _dos = 0;
        _translate_tabs = 0;
        _echo_returned_records_stream = NULL;

        _buffer.resize(lrecl + 1);

        //cerr << "Allocate _read_buffer with " << lrecl << " bytes\n";
        _read_buffer = new char[lrecl];

        if (NULL == _read_buffer)
        {
                cerr << "iwstring_data_source::_default_values:cannot allocate " << lrecl << " read buffer\n";
                _good = 0;
        }

        _chars_in_read_buffer = 0;
        _next_char_in_read_buffer_to_transfer_to_buffer = 0;

        _lrecl = lrecl;

        assert(ok());

        return;
}

void
iwstring_data_source::_setup_stream (const char * fname)
{
        _open = 0;
        _good = 0;

        const int strlen_fname = ::strlen(fname);

        if (1 == strlen_fname && '-' == fname[0])
        {
                _fd = 0;
                _open = 1;
                _good = 1;
        }
        else if (0 == strncmp(fname + strlen_fname - 3, ".gz", 3))
        {
                if (_gzfile.open_file(fname))
                {
                        _open = 1;
                        _good = 1;
                        _gzfile.set_record_delimiter(_record_delimiter);
                }
                else
                {
                        cerr << "iwstring_data_source::_setup_stream:cannot open '" << fname << "'\n";
                        _good = 0;
                }
        }
        else
        {
                _fd = iw_open_file(fname);

                //  March 2006. In some circumstances Windows may use file descriptor 0

#ifdef IWCYGWIN
                if (_fd >= 0)
#else
                if (_fd > 0)
#endif
                {
                        _open = 1;
                        _good = 1;
                }
                else
                        cerr << "iwstring_data_source:cannot open '" << fname << "'\n";
        }

        if (_good)
          _fname.set(fname, strlen_fname);

        return;
}

/*
  A data source is created with a suggested starting size for the buffer
*/

iwstring_data_source::iwstring_data_source (const char * fname, int lrecl) 
{
        assert (lrecl > 1);
        assert (NULL != fname);

        _default_values(lrecl);

        _setup_stream(fname);

        return;

}

iwstring_data_source::iwstring_data_source (const IWString & fname, int lrecl)
{
        assert (lrecl > 1);
        assert (fname.length());

        _default_values(lrecl);

        IWString tmp = fname;    // fname is const, so we cannot use it's .chars() function
        _setup_stream(tmp.null_terminated_chars());

        return;
}

iwstring_data_source::iwstring_data_source()
{
        int lrecl = 132;     // probably OK for many situations

        _default_values(lrecl);

        return;
}

iwstring_data_source::iwstring_data_source (int f)
{
  _default_values(STRING_DEFAULT_BUF_SIZE);

  _fd = f;
  _open = 1;
  _good = 1;

  return;
}

iwstring_data_source::iwstring_data_source (bool isstringbuffer, const char *stringbuffer, int stringbuffer_size)
{
  _default_values(stringbuffer_size);
  _isstringbuffer = isstringbuffer;
  _stringbuffer = stringbuffer;
  _stringbuffer_size = stringbuffer_size;
  _isstringbuffer_loaded = false;

  return;
}

iwstring_data_source::~iwstring_data_source()
{
        // TODO: for stringbuffer testing.              
  if(_isstringbuffer)
  {
  }
  else
  {
    if (_gzfile.active())    // destructor takes care of things
    ;
    else if (_open)
    {
    assert (_fd >= 0);
    if (0 != _fd)    // closing stdin messes up all subsequent file openings!!
      IW_FD_CLOSE(_fd);
    }

    _open = 0;

    if (NULL != _read_buffer)
      delete [] _read_buffer;

    return;
  }
}

int
iwstring_data_source::good() const
{
  if (_gzfile.active())
    return 1;

  return _good;
}

void
iwstring_data_source::set_record_delimiter (char s)
{
  _record_delimiter = s;

  if (_gzfile.active())
    _gzfile.set_record_delimiter(s);

  return;
}

//#define DEBUG_OK

int
iwstring_data_source::ok() const
{
  if (_gzfile.active())
    return 1;

  if (! _open)
    return 1;

  return _good;
}

int
iwstring_data_source::debug_print (std::ostream & os) const
{
  assert (os.good());

  if (_open)
    os << "Info on open string data source\n";
  else
    os << "Info on string data source not open\n";

  if (_gzfile.active())
    os << "gzip'd input\n";
  else if (_fd > 0)
    os << "File descriptor " << _fd << endl;
        //else if (_fd > 0)
        //  cerr << "File descriptor " << _fd << " currently at " << tellg() << endl;

  os << _lines_read << " lines have been read\n";
  os << _lines_which_are_returned << " lines passed the filters\n";

  if (_strip_leading_blanks)
    os << "Leading blanks will be stripped\n";
  if (_strip_trailing_blanks)
    os << "trailing blanks will be stripped\n";
  if (_convert_to_lowercase)
    os << "Records converted to lowercase\n";
  if (_convert_to_uppercase)
    os << "Records converted to uppercase\n";

  if (_filter_pattern.active())
    os << "Will filter lines '" << _filter_pattern.source() << "'\n";
  else
    os << "No filter pattern\n";

  if (_ignore_pattern.active())
    os << "Will ignore lines '" << _ignore_pattern.source() << "'\n";
  else
    os << "No ignore pattern\n";

  if (_record_buffered)
    os << "There is a record buffered\n" << _buffer << endl;
  else
    os << "No record buffered\n";

  os << _chars_in_read_buffer << " characters in read buffer\n";

  return 1;
}

int
iwstring_data_source::open (const char * fname)
{
  assert (ok());

  if (_open)
  {
    cerr << "iwstring_data_source::open: attempt to open already open file\n";
    return 0;
  }

  _setup_stream(fname);

  return _good;
}

int
iwstring_data_source::do_close()
{
  assert (_open);

  if (_gzfile.active())
    _gzfile.close_file();
  else
  {
    assert (_fd >= 0);
    if (_fd > 0)
      IW_FD_CLOSE(_fd);
  }

  _open = 0;

  return 1;
}

int
iwstring_data_source::push_record()
{
  assert (ok());

  _record_buffered = 1;

  return 1;    // why does this function return anything ?
}

/*
  In the interests of clarity (rather than efficiency), we always
  delete any existing pattern - even if the new one could fit in
  that space.
*/

int
iwstring_data_source::set_ignore_pattern (const const_IWSubstring & pattern)
{
  return _ignore_pattern.set_pattern(pattern);
}

int
iwstring_data_source::set_filter_pattern (const const_IWSubstring & pattern)
{
  return _filter_pattern.set_pattern(pattern);
}

/*
  This function applies all filters which are set to _buffer
  If it passes all of them, we return 1
  The order in which these filters are applied can be significant.
  Change this order and break existing code.
*/

#include <ctype.h>

int
iwstring_data_source::_apply_all_filters()
{
  //cerr << "_apply_all_filters:dos " << _dos << endl;
  if (_dos && _buffer.length())
  {
    //  if (static_cast<char> (13) == _buffer.last_item())   // could not put in ^M character because it would not compile under Linux
    //    cerr << "Doing chop\n";
    if (static_cast<char>(13) == _buffer.last_item())   // could not put in ^M character because it would not compile under Linux
      _buffer.chop();
  }
  //else
  //  cerr << "Not chopped\n";

  if (_strip_trailing_blanks)
    _buffer.strip_trailing_blanks();

  if (_strip_leading_blanks && _buffer.length() > 0 && isspace(_buffer[0]))
    _buffer.strip_leading_blanks();

  if (_skip_blank_lines && 0 == _buffer.length())
    return 0;

  if (_compress_spaces && _buffer.length())
    _buffer.compress_blanks();

  // Do case conversion before we apply string based filters.

  if (_convert_to_lowercase)
    _buffer.to_lowercase();

  if (_convert_to_uppercase)
    _buffer.to_uppercase();

  // If this record matches the ignore pattern, reject this line

  if (_ignore_pattern.active() && _ignore_pattern.matches(_buffer))
    return 0;

  // If we don't match the filter, reject this record.

  if (_filter_pattern.active() && ! _filter_pattern.matches(_buffer))
    return 0;

  _lines_which_are_returned++;
  return 1;
}

/*
  Frequently a programme will need to read down until it gets a record
  which matches a pattern. This function does that, and leaves the record
  buffered.
*/

int
iwstring_data_source::next_record_matches (const char * pattern)
{
  if (! ok())
    return 0;

  if (0 == strlen(pattern))    // what about EOF?
    return 0;

  if (! _record_buffered)
  {
    _fetch_record();

    if (_eof)
      return 0;
  }

  push_record();

  IW_Regular_Expression regexp(pattern);

  return regexp.matches(_buffer);
}

//#define DEBUG_COPY_NEXT_RECORD_FROM_READ_BUFFER_TO_BUFFER

/*
  We have some data in _read_buffer. Copy that data to _buffer.
  If we find a record delimiter, we return 1. Otherwise we just
  append the data to _buffer and return 0;
*/

int
iwstring_data_source::_copy_next_record_from_read_buffer_to_buffer()
{
  assert (_chars_in_read_buffer > 0);

  const char * s = _read_buffer + _next_char_in_read_buffer_to_transfer_to_buffer;   // where the new data starts

  const char * c = reinterpret_cast<const char *>(memchr(s, _record_delimiter, _chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer));

  //cerr << "Did we find a delimiter '" << (NULL != c) << "'\n";

  if (NULL != c)   // there is a record delimiter character in _read_buffer
  {
    int chars_to_copy = c - s;

    _buffer.strncat(s, chars_to_copy);

    _next_char_in_read_buffer_to_transfer_to_buffer += chars_to_copy + 1;

#ifdef DEBUG_COPY_NEXT_RECORD_FROM_READ_BUFFER_TO_BUFFER
    cerr << "copied " << chars_to_copy << " to buffer making complete record, buffer contains " << _buffer.length() << " chars\n";
#endif

    return 1;
  }

  // No record delimiter found, copy everything

  if (_buffer.number_elements() + _lrecl > _buffer.elements_allocated())
    _buffer.resize(_buffer.elements_allocated() + _buffer.elements_allocated());

  _buffer.strncat(s, _chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer);

  _chars_in_read_buffer = 0;
  _next_char_in_read_buffer_to_transfer_to_buffer = 0;

#ifdef DEBUG_COPY_NEXT_RECORD_FROM_READ_BUFFER_TO_BUFFER
  cerr << "Not a complete record yet, buffer " << _buffer.size() << ", allocated " << _buffer.elements_allocated() << endl;
#endif

  return 0;    // we did not find a whole record
}

//#define DEBUG_READ_MORE_DATA_INTO_READ_BUFFER

int
iwstring_data_source::_read_more_data_into_read_buffer()
{
  // for stringbuffer.
  if(_isstringbuffer)
  {
    if(!_isstringbuffer_loaded)
    {
      // copy string from _stringbuffer into _read_buffer, and set the _isstringbuffer_loaded flag.
      strncpy(_read_buffer, _stringbuffer,_stringbuffer_size);//_stringbuffer is NOT null terminated! use strncpy, strcpy will cause serious issue.
      _chars_in_read_buffer = _stringbuffer_size;
      _next_char_in_read_buffer_to_transfer_to_buffer = 0;
      _isstringbuffer_loaded = true;
      return 1;
    }
    else
    {
      _eof = 1;
      return 0;
    }
  }
  else
  {
    _chars_in_read_buffer = IW_FD_READ(_fd, _read_buffer, _lrecl);
    _next_char_in_read_buffer_to_transfer_to_buffer = 0;

#ifdef DEBUG_READ_MORE_DATA_INTO_READ_BUFFER
    cerr << "_read_more_data_into_read_buffer:read " << _chars_in_read_buffer << " bytes\n";
#endif

    if (_chars_in_read_buffer > 0)
      return 1;

    if (0 == _chars_in_read_buffer)
    {
      _eof = 1;
      if (_buffer.length())
      {
        cerr << "iwstring_data_source::_read_more_data_into_read_buffer:unterminated record\n";
        return 1;
      }

      return 0;
    }
    else
    {
      cerr << "iwstring_data_source::_read_more_data_into_read_buffer:fatal error\n";
      _good = 0;
    }

    return 0;
  }
}

int
iwstring_data_source::_fetch_record_into_buffer()
{
  _buffer.resize_keep_storage(0);

  // If we have some data buffered, start copying it to _buffer.

//#define DEBUG_FETCH_RECORD_INTO_BUFFER
#ifdef DEBUG_FETCH_RECORD_INTO_BUFFER
  cerr << "On entry to _fetch_record_into_buffer, _chars_in_read_buffer " << _chars_in_read_buffer << endl;
#endif

  if (_chars_in_read_buffer > 0 && _next_char_in_read_buffer_to_transfer_to_buffer < _chars_in_read_buffer)
  {
    if (_copy_next_record_from_read_buffer_to_buffer())
      return 1;
  }

  // _read_buffer is exhausted

  while (1)
  {
    if (! _read_more_data_into_read_buffer())
    {
      if (! _good)
        return 0;
      if (_eof)
        return 0;
      if (_buffer.length())
      {
        cerr << "iwstring_data_source::_fetch_record_into_buffer:unterminated record\n";
        _eof = 1;
        return 1;
      }
      return 0;
    }

    if (_eof)    // unterminated record at end of file
      return 1;

    if (_copy_next_record_from_read_buffer_to_buffer())
      return 1;
  }
}

/*
  this function fills _buffer with the next record
  Returns 1 if successful
*/

int
iwstring_data_source::_fetch_record()
{
  if(_isstringbuffer)
  {
    int nchars;
    if (! _fetch_record_into_buffer())
      return 0;

    nchars = _buffer.length();

    _lines_read++;

    if (nchars > _longest_record)
      _longest_record = nchars;

    if (_translate_tabs)
      _buffer.gsub('\t', ' ');

    return 1;
  }
  else
  {
    if (_eof)
     return 0;

    if (! _good)
      return 0;

    int nchars;

    if (_gzfile.active())
    {
      nchars = _gzfile.next_record(_buffer);
      if (nchars < 0)
      {
        cerr << "iwstring_data_source:_fetch_record:error reading gzip'd file\n";
        _good = 0;
        return 0;
      }
      else if (0 == nchars)
      {
        if (_gzfile.eof())
        {
          _eof = 1;
          return 0;
        }

        return 1;   // read zero length record, that's just fine
      }
    }
    else
    {
      if (! _fetch_record_into_buffer())
        return 0;

      nchars = _buffer.length();
    }

    _lines_read++;

    if (nchars > _longest_record)
      _longest_record = nchars;

    if (_translate_tabs)
      _buffer.gsub('\t', ' ');

    return 1;
  }
}

int
iwstring_data_source::most_recent_record (IWString & buffer)
{
  assert (_lines_read > 0);

  buffer = _buffer;

  return 1;
}

template <typename T>
int
iwstring_data_source::next_record (T & buffer)
{
        // TODO: for stringbuffer testing.
  if(_isstringbuffer)
  {
    if (_record_buffered)
    {
      _record_buffered = 0;
      if (_apply_all_filters())
      {
        buffer = _buffer;
        return 1;
      }
    }

    //#define DEBUG_NEXT_RECORD2

    // This infinite loop will loop until either we get a record that
    // matches the filters, or we hit EOF

    _record_buffered = 0;
    _buffer.resize_keep_storage(0);

    while (1)
    {
#ifdef DEBUG_NEXT_RECORD2
      cerr << "In infinite loop fetching records " << _good << endl;
#endif
      if (! _fetch_record())
        return 0;

#ifdef DEBUG_NEXT_RECORD2
      cerr << "Fetched a record, good " << _good << ", " << _buffer.length() << " chars\n";
#endif

      if (_apply_all_filters())
      {
        buffer = _buffer;
#ifdef DEBUG_NEXT_RECORD2
        cerr << "Got a record that passed all filters " << _good << endl;
#endif
        return 1;
      }
    }

    assert (NULL == "should never come here");
  }
  else
  {

    if (! _open)
      return 0;

    if (_eof)
      return 0;

    if (! _good)
    {
      cerr << "iwstring_data_source::next_record: is not good\n";
      return 0;
    }

    if (_record_buffered)
    {
      _record_buffered = 0;
      if (_apply_all_filters())
      {
        buffer = _buffer;
        return 1;
      }
    }

    //#define DEBUG_NEXT_RECORD2

    // This infinite loop will loop until either we get a record that
    // matches the filters, or we hit EOF

    _record_buffered = 0;
    _buffer.resize_keep_storage(0);

    while (1)
    {
#ifdef DEBUG_NEXT_RECORD2
      cerr << "In infinite loop fetching records " << _good << endl;
#endif
      if (! _fetch_record())
        return 0;

#ifdef DEBUG_NEXT_RECORD2
      cerr << "Fetched a record, good " << _good << ", " << _buffer.length() << " chars\n";
#endif

      if (_apply_all_filters())
      {
        buffer = _buffer;
#ifdef DEBUG_NEXT_RECORD2
        cerr << "Got a record that passed all filters " << _good << endl;
#endif
        return 1;
      }
    }

    assert (NULL == "should never come here");
  }
}

template int iwstring_data_source::next_record(IWString &);
template int iwstring_data_source::next_record(const_IWSubstring &);

/*
  Reads records until we get to one containing pattern.
  If pattern starts with '^', it must match at the beginning
  of the record.
  We return 1 if successful, 0 otherwise
  The next record read after this call, will be the record after
  the one which contains the pattern
*/

int
iwstring_data_source::skip_past (const char * pattern)
{
        assert (ok());

        if (0 == strlen(pattern))
                return 0;

        IW_Regular_Expression regexp(pattern);

        if (_record_buffered && regexp.matches(_buffer))
        {
          _record_buffered = 0;
          return 1;
        }

        while (_fetch_record())
        {
                if (regexp.matches(_buffer))
                        return 1;
        }

        return 0;
}

/*
Very similar to skip_past, but we push the record when we are done
so that the next record read will be one with the pattern
*/

int
iwstring_data_source::skip_to (const char * pattern)
{
        assert (ok());

        if (skip_past(pattern))
        {
                push_record();
                return 1;
        }
        else
                return 0;
}

/*
tellg is complicated because of the buffering we do.

When someone requests the current offset, they expect it to
be at the end of _buffer.
*/

off_t
iwstring_data_source::tellg() const
{
  assert (ok());

  if (_gzfile.active())
    return _gzfile.tellg();

  if(_isstringbuffer)
  {
    off_t current_offset = 1;
    off_t rc = current_offset - static_cast<off_t>(_chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer);
    return rc;
  }
  else
  {
    assert (_fd >= 0);

    off_t current_offset = IW_FD_LSEEK(_fd, 0, SEEK_CUR);

#ifdef DEBUG_TELLG
    cerr << "iwstring_data_source::tellg:current_offset " << current_offset << " _chars_in_read_buffer " << _chars_in_read_buffer << " _next_char_in_read_buffer_to_transfer_to_buffer " << _next_char_in_read_buffer_to_transfer_to_buffer << endl;
#endif

    if (-1 == current_offset && 0 == _fd)
    {
      cerr << "iwstring_data_source::tellg:not supported on pipe read\n";
      return static_cast<off_t>(0);
    }

    if (current_offset < static_cast<off_t>(0))   // cannot happen, fix sometime
    {
      cerr << "iwstring_data_source::tellg:negative current offset " << current_offset << ", file " << _fd << endl;
      assert (NULL == "This is very bad");
      return current_offset;
    }

    assert (current_offset >= static_cast<off_t>(_chars_in_read_buffer));

    off_t rc = current_offset - static_cast<off_t>(_chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer);

    if (rc < 0)          // silly test, cannot happen. Fix sometime
    {
      cerr << "iwstring_data_source::tellg:negative offset " << rc << endl;
      cerr << _chars_in_read_buffer << " chars in read buffer\n";
      cerr << _next_char_in_read_buffer_to_transfer_to_buffer << " next to transfer\n";
      cerr << current_offset << " current offset\n";
    }

    assert (rc >= static_cast<off_t>(0));
    return rc;
  }

}

int
iwstring_data_source::seekg (off_t zoffset, int whence)
{
//cerr << "iwstring_data_source::seekg:seeking " << zoffset << " by " << whence << endl;

  assert (ok());

  _record_buffered = 0;

  if (_gzfile.active())
  {
    _eof = 0;
    return _gzfile.seekg(zoffset);
  }

  if (! _good)
  {
    cerr << "iwstring_data_source::seek:file descriptor is bad\n";
    return 0;
  }

  if (_fd <= 0 || ! _open)
  {
    cerr << "iwstring_data_source::seekg: NULL file pointer\n";
    return 0;
  }

  if (SEEK_CUR == whence)
  {
    off_t current_offset = tellg();
    if (current_offset < 0)    // cannot be negative, fix sometime
    {
      cerr << "iwstring_data_source::seekg:cannot determine current offset\n";
      return 0;
    }

    zoffset = current_offset + zoffset;
  }

  //if (_is->eof())
  //  cerr << "File is at EOF\n";

  if (_eof)    // the data source is still valid, although at EOF
    _eof = 0;

  off_t rc = IW_FD_LSEEK(_fd, zoffset, SEEK_SET);

  if (rc < 0)
  {
    cerr << "iwstring_data_source::seekg:cannot seek to " << zoffset << endl;
    _good = 0;
    return 0;
  }

//off_t current_offset = IW_FD_LSEEK(_fd, 0, SEEK_CUR);

  //cerr << "After seek " << whence << " offset " << current_offset << " rc " << rc << endl;

  _chars_in_read_buffer = 0;
  _record_buffered = 0;
  _next_char_in_read_buffer_to_transfer_to_buffer = 0;

  //cerr << "After seek to " << zoffset << " good is " << _is->good() << endl;

  return 1;
}

class IWSDS_State
{
private:
        off_t _offset;
        int _save_lines_read;
        int _save_longest_record;
        int _save_record_buffered;

        IWString _save_buffer;

public:
        IWSDS_State();
        ~IWSDS_State();

        void set_offset (off_t p) { _offset = p;};
        off_t offset() const { return _offset;}

        void set_lines_read (int lr) { _save_lines_read = lr;}
        int  lines_read() const { return _save_lines_read;}

        void set_longest_record (int lr) { _save_longest_record = lr;}
        int  longest_record() const { return _save_longest_record;}

        void set_record_buffered (int lr) { _save_record_buffered = lr;}
        int  record_buffered() const { return _save_record_buffered;}

        void set_buffer (const IWString & b) { _save_buffer = b;}
        IWString & buffer() { return _save_buffer;}
};

IWSDS_State::IWSDS_State()
{
        _save_lines_read = -5;

        return;
}

IWSDS_State::~IWSDS_State()
{
}

off_t
iwstring_data_source::file_size()
{
        if (_gzfile.active())
        {
                cerr << "iwstring_data_source::file_size:unsupported on gzip'd input\n";
                return 0;
        }

        if (! _good || ! _open)
                return 0;

        if (1 == _fd)
        {
                cerr << "iwstring_data_source::file_size:undefined for stdin\n";
                return 0;
        }

        const off_t current_position = IW_FD_LSEEK(_fd, 0, SEEK_CUR);

        const off_t rc = IW_FD_LSEEK(_fd, 0, SEEK_END);

        const off_t back_again = IW_FD_LSEEK(_fd, current_position, SEEK_SET);

        assert (back_again == current_position);

        return rc;
}

/*
The optional parameter, STOP_COUNTING_WHEN, is just used by
at_least_X_records_remaining.
*/

int
iwstring_data_source::records_remaining (int stop_counting_when)
{
        if(_isstringbuffer)
        {
                // TODO: 
                return _stringbuffer_size; 
        }
        else
        {
                if (_gzfile.active())
                {
                        cerr << "iwstring_data_source::records_remaining:unsupported on gzip'd input\n";
                        return 0;
                }

                if (! _open || ! _good)
                {
                        cerr << "iwstring_data_source::records_remaining:file not open, or bad\n";
                        return 0;
                }

                if (_eof)
                        return 0;

                if (0 == _fd)
                {
                        cerr << "iwstring_data_source::records_remaining:undefined for stdin\n";
                        return 0;
                }

                assert (_fd >= 0);

                if (_eof)
                        return 0;

                // Get the offset before we save the state because that destroys _read_more_data_into_read_buffer

                off_t current_offset = tellg();

                IWSDS_State zstate;

                _save_state(zstate);

                if (0 == stop_counting_when)
                        stop_counting_when = std::numeric_limits<int>::max();

                int rc = 0;

                //#define DEBUG_RECORDS_REMAINING
        #ifdef DEBUG_RECORDS_REMAINING
                cerr << "iwstring_data_source::records_remaining:seek to " << current_offset << endl;
        #endif

                if (IW_FD_LSEEK(_fd, current_offset, SEEK_SET) < 0)
                {
                        cerr << "iwstring_data_source::records_remaining:cannot seek to " << current_offset << endl;
                        return 0;
                }

        #define RRBUFSIZE 4096

                char tmpbuffer[RRBUFSIZE];

                int last_char_was_record_terminator = 0;

                // If we are at the end of the file, we will read zero characters
                // and last_char_was_record_terminator will never get set

                int read_some_data = 0;

                while (1)
                {
                        int chars_read = IW_FD_READ(_fd, tmpbuffer, RRBUFSIZE);

        #ifdef DEBUG_RECORDS_REMAINING
                        cerr << "Read " << chars_read << " characters, rc = " << rc << endl;
        #endif

                        if (0 == chars_read)
                                break;

                        read_some_data = 1;

                        last_char_was_record_terminator = 0;

                        for (int i = 0; i < chars_read; i++)
                        {
                                if (_record_delimiter != tmpbuffer[i])
                                        continue;

                                rc++;
                                if (i == chars_read - 1)
                                        last_char_was_record_terminator = 1;

                                if (rc >= stop_counting_when)
                                        break;
                        }

                        if (rc >= stop_counting_when)
                                break;
                }

                // Check for an unterminated final record

                if (rc >= stop_counting_when)
                        ;
                else if (! read_some_data)   // there cannot have been an unterminated record
                        ;
                else if (! last_char_was_record_terminator)
                {
                        cerr << "iwstring_data_source::records_remaining:unterminated last record\n";
                        rc++;
                }

                _open = 0;     // look at _restore_state to see why this is important

                _restore_state(zstate);

                _eof = 0;

        #ifdef DEBUG_RECORDS_REMAINING
                cerr << "After records remaining, offset " << tellg() << endl;
        #endif

                return rc;
        }
}

int
iwstring_data_source::at_least_X_records_remaining (int n)
{
        assert (n > 0);

        return records_remaining(n);
}

/*
Note the many possible problems with grep.
For example, if they want strip_leading_blanks(), then should
grep() do that before attempting a match? Well, it doesn't
*/

int
iwstring_data_source::grep (IW_Regular_Expression & rx)
{
        assert (ok());
        assert (rx.ok());

        IWSDS_State zstate;

        if (! _save_state(zstate))
                return 0;

        int rc = 0;
        while (_fetch_record())
        {
                if (rx.matches(_buffer))
                        rc++;
        }

        _restore_state(zstate);

        return rc;
}

/*
  Note that we do NOT initialise the count array
*/

int
iwstring_data_source::grep (int n,
                            IW_Regular_Expression * rx,
                            int * count)
{
        assert (ok());

        IWSDS_State zstate;

        if (! _save_state(zstate))
                return 0;

        int rc = 0;

        while (_fetch_record())
        {
                for (int i = 0; i < n; i++)
                {
                        if (rx[i].matches(_buffer))
                        {
                                count[i]++;
                                rc++;
                        }
                }
        }

        return rc;
}

int
iwstring_data_source::grep (const const_IWSubstring & c)
{
        IW_Regular_Expression rx(c);

        if (! rx.ok())
        {
                cerr << "iwstring_data_source::grep: cannot parse regexp '" << c << "'\n";
                return 0;
        }

        return grep(rx);
}

int
iwstring_data_source::count_records_starting_with(const const_IWSubstring & s)
{
  assert (ok());

  IWSDS_State zstate;

  if (! _save_state(zstate))
    return 0;

  int rc = 0;

  while (_fetch_record())
  {
    if (_buffer.starts_with(s))
      rc++;
  }

  _restore_state(zstate);

  return rc;
}

/*
  When we save our state, we discard _read_buffer
*/

int 
iwstring_data_source::_save_state (IWSDS_State & zstate)
{
  assert (ok());

  if (! _open)
          return 0;

  if (_gzfile.active())
  {
          cerr << "iwstring_data_source::_save_state:save state not available with gzip'd data\n";
          return 0;
  }

// because we have discarded the read buffer, we must position the file descriptor to the apparent position

  const off_t o = tellg();

  zstate.set_offset(o);

  off_t s = IW_FD_LSEEK(_fd, o, SEEK_SET);

  if (s < 0)
  {
    cerr << "iwstring_data_source::_save_state:cannot seek to " << o << endl;
    return 0;
  }

  _chars_in_read_buffer = 0;
  _next_char_in_read_buffer_to_transfer_to_buffer = 0;

  zstate.set_lines_read(_lines_read);
  zstate.set_longest_record(_longest_record);
  zstate.set_record_buffered(_record_buffered);

  zstate.set_buffer(_buffer);

  return 1;
}

/*
  Might not be entirely robust
*/

int
iwstring_data_source::_restore_state (IWSDS_State & zstate)
{
  if (_fd <= 0)
  {
          cerr << "iwstring_data_source::_restore_state: no file pointer\n";
          return 0;
  }

  if (_eof)
    _eof = 0;

  if (! _open)
    _open = 1;

  assert (ok());

#ifdef DEBUG_RESTORE_STATE
  cerr << "Seeking back to " << zstate.offset() << endl;
#endif

  off_t s = IW_FD_LSEEK(_fd, zstate.offset(), SEEK_SET);

  if (s < 0)
  {
    cerr << "iwstring_data_source::_restore_state:cannot seek back to " << zstate.offset() << endl;
    _good = 0;
    return 0;
  }

  _lines_read = zstate.lines_read();
  _longest_record = zstate.longest_record();
  _record_buffered = zstate.record_buffered();

  _buffer = zstate.buffer();

  _chars_in_read_buffer = 0;
  _next_char_in_read_buffer_to_transfer_to_buffer = 0;

  return 1;
}

int
iwstring_data_source::_write_read_buffer(std::ostream & output,
                                         size_t & nbytes)
{
  int chars_to_write;
  if (static_cast<size_t>(_chars_in_read_buffer) >= nbytes)
    chars_to_write = nbytes;
  else
    chars_to_write = _chars_in_read_buffer;

  output.write(_read_buffer, chars_to_write);

  _chars_in_read_buffer -= chars_to_write;
  _next_char_in_read_buffer_to_transfer_to_buffer -= chars_to_write;

  nbytes -= chars_to_write;

  return output.good();
}

int
iwstring_data_source::_write_read_buffer(IWString_and_File_Descriptor & output,
                                         size_t & nbytes)
{
  int chars_to_write;
  if (static_cast<size_t>(_chars_in_read_buffer) >= nbytes)
    chars_to_write = nbytes;
  else
    chars_to_write = _chars_in_read_buffer;

  output.write(_read_buffer, chars_to_write);

  _chars_in_read_buffer -= chars_to_write;
  _next_char_in_read_buffer_to_transfer_to_buffer -= chars_to_write;

  nbytes -= chars_to_write;

  return output.good();
}

/*
  Echo a certain number of bytes to an ostream.

  Note that we don't keep track of records read or anything like that,
  this is usually done as a post-processing step: read the file,
  determine what gets read, and then echo certain parts of the input.
*/

int
iwstring_data_source::echo(std::ostream & output, size_t nbytes)
{
  assert (ok());

  if (! _open)
  {
    cerr << "iwstring_data_source::echo: not open\n";
    return 0;
  }

  _buffer.resize_keep_storage(0);
  _record_buffered = 0;

  if (_chars_in_read_buffer)
  {
    if (! _write_read_buffer(output, nbytes))
      return 0;

    if (0 == nbytes)
      return 1;
  }

  while (1)
  {
    if (! _read_more_data_into_read_buffer())
    {
      if (! _good)
        return 0;
    }

    if (! _write_read_buffer(output, nbytes))
      return 0;

    if (0 == nbytes)
      return output.good();
  }

  return output.good();
}

int
iwstring_data_source::echo(IWString_and_File_Descriptor & output, size_t nbytes)
{
  assert (ok());

  if (! _open)
  {
    cerr << "iwstring_data_source::echo: not open\n";
    return 0;
  }

  _buffer.resize_keep_storage(0);
  _record_buffered = 0;

  if (_chars_in_read_buffer)
  {
    if (! _write_read_buffer(output, nbytes))
      return 0;

    if (0 == nbytes)
      return 1;
  }

  while (1)
  {
    if (! _read_more_data_into_read_buffer())
    {
      if (! _good)
        return 0;
    }

    if (! _write_read_buffer(output, nbytes))
      return 0;

    if (0 == nbytes)
      return output.good();
  }

  return output.good();
}

int
iwstring_data_source::skip_records (int nskip)
{
        assert (ok());
        assert (nskip >= 0);

        const_IWSubstring buffer;

        for (int i = 0; i < nskip; i++)
        {
                if (! next_record(buffer))
                {
                        cerr << "iwstring_data_source::skip_records: EOF at record " << i << " of " << nskip << endl;
                        return 0;
                }
        }


        return 1;
}

int
iwstring_data_source::skip_records (IW_Regular_Expression & rx, int nskip)
{
        assert (ok());
        assert (nskip >= 0);

        const_IWSubstring buffer;

        int nfound = 0;

        while (nfound < nskip)
        {
                if (! next_record(buffer))
                {
                        cerr << "iwstring_data_source::skip_records: EOF at record " << _lines_read << " of " << nskip << endl;
                        return 0;
                }

                if (rx.matches(buffer))
                {
                        nfound++;
                        if (nfound >= nskip)
                                return 1;
                }
        }

        return 1;
}

int
iwstring_data_source::echo_records (std::ostream & os, int necho)
{
        assert (ok());
        assert (necho >= 0);

        const_IWSubstring buffer;

        for (int i = 0; i < necho && os.good(); i++)
        {
                if (! next_record(buffer))
                {
                        cerr << "iwstring_data_source::echo_records: EOF at record " << i << " of " << necho << endl;
                        return 0;
                }

                os << buffer << '\n';
        }

        return os.good();
}

int
iwstring_data_source::echo_records (IWString_and_File_Descriptor & os, int necho)
{
        assert (ok());
        assert (necho >= 0);

        const_IWSubstring buffer;

        for (int i = 0; i < necho && os.good(); i++)
        {
                if (! next_record(buffer))
                {
                        cerr << "iwstring_data_source::echo_records: EOF at record " << i << " of " << necho << endl;
                        return 0;
                }

                os << buffer << '\n';

                if (os.size() > 32768)
                        os.write_whole_blocks_shift_unwritten();
        }

        return os.good();
}

/*
  NBYTES will be set to the number of bytes copied
*/

int
iwstring_data_source::_copy_read_buffer_to_destination(void * destination,
                                                       int & nbytes)
{
  assert (_chars_in_read_buffer > 0);

  int bytes_to_copy;

  int bytes_available = _chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer;

  if (bytes_available >= nbytes)
    bytes_to_copy = nbytes;
  else
    bytes_to_copy = bytes_available;

  ::memcpy(destination, _read_buffer + _next_char_in_read_buffer_to_transfer_to_buffer, bytes_to_copy);

  _chars_in_read_buffer -= bytes_to_copy;
  _next_char_in_read_buffer_to_transfer_to_buffer -= bytes_to_copy;

  nbytes = bytes_to_copy;

  return 1;
}

//#define DEBUG_READ_BYTES

int
iwstring_data_source::read_bytes (void * destination, 
                                  size_t bytes_requested)
{
#ifdef DEBUG_READ_BYTES
  cerr << "read_bytes:looking for " << bytes_requested << " bytes\n";
#endif

  assert (ok());
  assert (bytes_requested > 0);

  _buffer.resize_keep_storage(0);

  if (_gzfile.active())
    return _gzfile.read_bytes(destination, bytes_requested);

  if (! _open || ! _good)
  {
    cerr << "iwstring_data_source::read_bytes: not open or bad\n";
    return 0;
  }

  size_t bytes_written = 0;

  if (_chars_in_read_buffer)
  {
#ifdef DEBUG_READ_BYTES
    cerr << "There are " << _chars_in_read_buffer << " chars in read buffer\n";
#endif

    int nbytes = bytes_requested;
    _copy_read_buffer_to_destination(destination, nbytes);

#ifdef DEBUG_READ_BYTES
    cerr << nbytes << " bytes copied from read buffer\n";
#endif

    bytes_requested -= nbytes;
    bytes_written = nbytes;

    if (0 == bytes_requested)   // satisfied the entire request from the buffer
      return bytes_written;
  }

  while (bytes_requested > 0)
  {
#ifdef DEBUG_READ_BYTES
    cerr << bytes_written << " written so far, requested " << bytes_requested << endl;
#endif

    if (! _read_more_data_into_read_buffer())
    {
      //    cerr << "Could not read more data, good " << _good << endl;
      if (! _good)
        return 0;
      if (_buffer.length())
        cerr << "iwstring_data_source::_fetch_record_into_buffer:unterminated record\n";
      return 0;
    }

    int nbytes = bytes_requested;     // nbytes will be set to how many bytes get copied

    _copy_read_buffer_to_destination(reinterpret_cast<char *>(destination) + bytes_written, nbytes);

#ifdef DEBUG_READ_BYTES
    cerr << "Copied " << nbytes << " bytes from read buffer\n";
#endif

    if (0 == nbytes)
      return bytes_written;

    bytes_written += nbytes;
    bytes_requested -= nbytes;

#ifdef DEBUG_READ_BYTES
    cerr << "At end of loop, written " << bytes_written << " still requrested " << bytes_requested << endl;
#endif
  }

  return bytes_written;
}

#define IWSTRDS_BUF_SIZE 4096

size_t
iwstring_data_source::copy_raw_bytes (void * destination, const size_t bytes_to_copy)
{
  IWSDS_State zstate;

  if (! _save_state(zstate))
    return 0;

  const size_t rc = _copy_raw_bytes(destination, bytes_to_copy);

  _restore_state(zstate);

  return rc;
}

size_t 
iwstring_data_source::_copy_raw_bytes(void * destination, const size_t bytes_to_copy)
{
  unsigned char buf[IWSTRDS_BUF_SIZE];

  size_t ncopy = bytes_to_copy;

//cerr << "iwstring_data_source::_copy_raw_bytes:copying " << bytes_to_copy << " bytes, read buffer " << _chars_in_read_buffer << ", _next_char_in_read_buffer_to_transfer_to_buffer " << _next_char_in_read_buffer_to_transfer_to_buffer << endl;

  if (_chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer > 0)
  {
    size_t copy_from_read_buffer = _chars_in_read_buffer - _next_char_in_read_buffer_to_transfer_to_buffer;

    if (copy_from_read_buffer > bytes_to_copy)
      copy_from_read_buffer = bytes_to_copy;

    ::memcpy(destination, _read_buffer + _next_char_in_read_buffer_to_transfer_to_buffer, copy_from_read_buffer);

    if (copy_from_read_buffer == bytes_to_copy)   // we are done
      return bytes_to_copy;

    ncopy -= copy_from_read_buffer;
    destination += copy_from_read_buffer;
  }

// At this stage the read buffer is empty, so it is just about transferring data

  while (ncopy > 0)
  {
    int bytes_read = IW_FD_READ(_fd, buf, IWSTRDS_BUF_SIZE);
//  cerr << "Read " << bytes_read << endl;

    if (0 == bytes_read)
      break;

    ::memcpy(destination, buf, bytes_read);
    destination += bytes_read;
    ncopy -= bytes_read;

    if (bytes_read != IWSTRDS_BUF_SIZE)
      break;
  }

  return bytes_to_copy;
}

/*
  Might be safe enough to just check for 0 == _fd
*/

int
iwstring_data_source::is_pipe() const
{
  if (! is_open())
    return 0;

  struct stat s;

  const auto rc = fstat(_fd, &s);

  return S_ISFIFO(s.st_mode);
}
