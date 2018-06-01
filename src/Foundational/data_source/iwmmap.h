#ifndef IWSTRING_DS_MMAP_H
#define IWSTRING_DS_MMAP_H

#include<algorithm>

#include "iwcrex.h"

class IW_MMapd_File
{
  protected:
    const unsigned char * _p;
    size_t _len;

    int _good;

//  protected functions

    int  _close_file();

  private:
    int _fd;

    IWString _fname;

//  private functions

    int  _open_file(const char * fname);

  public:
    IW_MMapd_File();
    IW_MMapd_File(const char * fname);
    ~IW_MMapd_File();

    const IWString & fname () const { return _fname;}

    int open (const char * fname);
    int do_close ();

    int is_pipe () const { return 0 == _fd;}

    size_t file_size() const { return _len;}

    int is_open () const { return nullptr != _p;}

    int do_madvise (const int);    // passed directly to madvise
};

class IW_File_Contents
{
  protected:
    const unsigned char * _p;
    size_t _len;

    int _good;

    IWString _fname;
    resizable_array<unsigned char> _acc;

//  private functions

    int _read_the_data(const int fd);

  public:
    IW_File_Contents();
    IW_File_Contents(const char *);
    ~IW_File_Contents();

    int open (const char * fname);
    int do_close ();

    const IWString & fname () const { return _fname;}

    size_t file_size() const { return _len;}

    int is_open () const { return nullptr != _p;}
};

/*
  Conceptually, we have an interface that knows how to read records, seek and so forth.
  It has a pointer to some data.
  this can come from either a memory mapped file, or from someone who might have read
  the contents of stdin into memory
*/

template <typename T>
class IW_Storage_Reader : public T
{ 
  protected:
    const unsigned char * _myptr;

    int _lines_read;
    int _dos;
    int _record_buffered;

    int _strip_leading_blanks;
    int _strip_trailing_blanks;
    int _skip_blank_lines;    // not implemented yet

    IW_Regular_Expression _ignore_pattern;
    IW_Regular_Expression _filter_pattern;

//  private functions

    void _default_values();
    int  _find_next_record (const_IWSubstring & buffer);
    int  _fetch_previous_record (const_IWSubstring & buffer);
    int _next_record (const_IWSubstring & buffer);
    int _next_record_inner (const_IWSubstring & buffer);

  public:
    IW_Storage_Reader();
    IW_Storage_Reader(const char * fname);
    ~IW_Storage_Reader();

    int eof () const;

    int open (const char * fname);

    int  records_remaining (const int = 0) const;

    void push_record () { _record_buffered = 1;}

    int next_record (const_IWSubstring &);
    int next_record (IWString &);

    int lines_read() const { return _lines_read;}
    void set_dos (int d) { _dos = d;}
    void set_strip_leading_blanks (int s=1) { _strip_leading_blanks = s; }
    void set_strip_trailing_blanks (int s=1) { _strip_trailing_blanks = s; }
    void set_record_delimiter (char d) { return;}   // noop for now
    void set_skip_blank_lines (int s=1) { _skip_blank_lines = s; }
//  int  set_filter_pattern (const const_IWSubstring &);
    int  set_ignore_pattern (const const_IWSubstring &);

    int skip_records (const int nskip);
    int skip_records (IW_Regular_Expression & rx, int nskip);

    int skip_to   (const char *);
    int skip_past (const char *);

    size_t tellg() const;
    int    seekg (size_t, const int = 0);

    int  echo (std::ostream &, size_t);    // echo's bytes
    int  echo (IWString_and_File_Descriptor &, size_t);    // echo's bytes

    int echo_records (std::ostream & os, int necho);    // echo's records
    int echo_records (IWString_and_File_Descriptor & os, int necho);    // echo's records

    int grep (IW_Regular_Expression & rx);

    int most_recent_record (IWString &);

    int  count_records_starting_with(const const_IWSubstring &);   // the most common thing we do with regexps

    size_t copy_raw_bytes (void *, const size_t) const;    // does a save and restore of the state, so it will not advance the file pointer - just got too complicated...
};

class IWString_Data_Source_MMAP : public IW_Storage_Reader<IW_MMapd_File>
{
  private:
  public:
    IWString_Data_Source_MMAP();
    IWString_Data_Source_MMAP(const char * fname);

    int good () const { return _good;}
    int ok () const { return _good;}
};

class IWString_Data_Source_RAM : public IW_Storage_Reader<IW_File_Contents>
{
  private:
  public:
    IWString_Data_Source_RAM();
    IWString_Data_Source_RAM(const char * fname);

    int good () const { return _good;}
    int ok () const { return _good;}
};

#if defined(STORAGE_READER_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

template <typename T>
void
IW_Storage_Reader<T>::_default_values()
{
  _myptr = nullptr;
  _dos = 1;
  _record_buffered = 0;
  _lines_read = 0;
  _strip_leading_blanks = 0;
  _strip_trailing_blanks = 0;
  _skip_blank_lines = 0;

  return;
}

template <typename T>
IW_Storage_Reader<T>::IW_Storage_Reader()
{
  _default_values();
}

template <typename T>
IW_Storage_Reader<T>::IW_Storage_Reader (const char * fname) : T(fname)
{
  _default_values();

  _myptr = T::_p;

  return;
}

template <typename T>
IW_Storage_Reader<T>::~IW_Storage_Reader()
{
  _myptr = nullptr;
  return;
}

template <typename T>
int
IW_Storage_Reader<T>::_find_next_record (const_IWSubstring & buffer)
{
  const unsigned char * p = T::_p;

  for (size_t i = _myptr - p; i < T::_len; ++i)
  {
    if ('\n' != p[i])
      continue;

    int nchars = i - (_myptr - p);

//  cerr << "Found newline at i = " << i << ", nchars " << nchars << endl;

    if (_dos && nchars > 0 && static_cast<unsigned char>(13) == p[i-1])
      nchars--;

    buffer.set(reinterpret_cast<const char *>(_myptr), nchars);

    _myptr = p + i + 1;

    return 1;
  }
  
  return 0;
}

template <typename T>
int
IW_Storage_Reader<T>::open(const char * fname)
{
  if (! T::open(fname))
    return 0;

  _myptr = T::_p;

  return 1;
}

template <typename T>
int
IW_Storage_Reader<T>::eof () const
{
  if (! T::_good || nullptr == T::_p)
  {
    cerr << "IW_Storage_Reader:eof:closed or invalid state\n";
    return -1;
  }

  assert(_myptr >= T::_p);

  return static_cast<size_t>(_myptr - T::_p) == T::_len;
}

template <typename T>
int
IW_Storage_Reader<T>::next_record (const_IWSubstring & buffer)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  if (! _next_record(buffer))
    return 0;

  if (_strip_leading_blanks)
    buffer.strip_trailing_blanks();

  if (_strip_trailing_blanks)
    buffer.strip_trailing_blanks();

  return 1;
}

template <typename T>
int
IW_Storage_Reader<T>::_next_record_inner (const_IWSubstring & buffer)
{

//#define DEBUG_NEXT_RECORD
#ifdef DEBUG_NEXT_RECORD
  cerr << _lines_read << " offset " << (_myptr - _p) << " len " << _len << endl;
#endif

  assert (_myptr >= T::_p);

  if (static_cast<size_t>(_myptr - T::_p) == T::_len)    // eof
    return 0;

#ifdef DEBUG_NEXT_RECORD
  cerr << "Not EOF\n";
#endif

  if (! _find_next_record(buffer))
    return 0;

  _lines_read++;

  return 1;
}

template <typename T>
int
IW_Storage_Reader<T>::set_ignore_pattern (const const_IWSubstring & s)
{
  return _ignore_pattern.set_pattern(s);
}

template <typename T>
int
IW_Storage_Reader<T>::_next_record (const_IWSubstring & buffer)
{
  if (_record_buffered)
  {
    const unsigned char * msave = _myptr;

    const int rc = _fetch_previous_record(buffer);

    _myptr = msave;

    _record_buffered = 0;

    return rc;
  }

  while (_next_record_inner(buffer))
  {
    if (_strip_trailing_blanks)
      buffer.strip_trailing_blanks();
    if (_strip_leading_blanks)
      buffer.strip_leading_blanks();
    if (_skip_blank_lines && 0 == buffer.length())
      continue;
    if (_ignore_pattern.active() && _ignore_pattern.matches(buffer))
      continue;
    return 1;
  }

  buffer.make_empty();

  return 0;
}

template <typename T>
int
IW_Storage_Reader<T>::next_record (IWString & buffer)
{
  const_IWSubstring tmp;
  if (! next_record(tmp))
    return 0;

  buffer = tmp;

  return 1;
}

/*
   Look at the string below, 'abc' with a trailing dos carriage return q.
   _myptr will be 8, 7 will be the previous record's newline
   6 might be a dos carriage return

   2345678
   .abcq.
   ^    ^
*/

template <typename T>
int
IW_Storage_Reader<T>::_fetch_previous_record (const_IWSubstring & buffer)
{
  const unsigned char * p = T::_p;

  if (p == _myptr)    // there is no previous record
    return 0;

  assert (_myptr > p);

//#define DEBUG_FETCH_PREVIOUS_RECORD
#ifdef DEBUG_FETCH_PREVIOUS_RECORD
  cerr << "_fetch_previous_record:offset " << (_myptr - p) << endl;
#endif

  for (int i = _myptr - 2 - p; i >= 0; --i)     // skip back past previous record newline
  {
    if ('\n' != p[i])
      continue;

    int nchars = (_myptr - p) - 2 - i;
    assert (nchars >= 0);

#ifdef DEBUG_FETCH_PREVIOUS_RECORD
    cerr << "previous record contains " << nchars << " characters possible dos " << (_myptr - _p - 2) << endl;
#endif

    if (_dos && nchars > 0 && static_cast<unsigned char>(13) == p[_myptr - p - 2])
      nchars--;

    buffer.set(reinterpret_cast<const char *>(p + i + 1), nchars);

    return 1;
  }

// If we come to here, the previous record was the first one, myptr = 5
//     abcq.
//         ^

  int nchars = _myptr - p - 1;
  if (_dos && nchars > 0 && static_cast<unsigned char>(13) == p[_myptr - p - 2])
    nchars--;

  buffer.set(reinterpret_cast<const char *>(p), nchars);

  return 1;
}

template <typename T>
int
IW_Storage_Reader<T>::records_remaining (const int stop_counting_when) const
{
  const unsigned char * p = T::_p;

  if (! T::_good || nullptr == p)
    return 0;

  int rc = 0;

  for (size_t i = _myptr - p; i < T::_len; ++i)
  {
    if ('\n' != p[i])
      continue;

    rc++;
    if (stop_counting_when > 0 && rc >= stop_counting_when)
      return stop_counting_when;
  }

  return rc;
}

template <typename T>
int
IW_Storage_Reader<T>::seekg (size_t o,
                                  const int flag)
{
  const unsigned char * p = T::_p;

  if (! T::_good || nullptr == p)
    return 0;

  if (o > T::_len)
  {
    cerr << "IWString_Data_Source_MMAP::seekg:offset " << o << " out of range, file size " << T::_len << endl;
    return 0;
  }

  if (SEEK_CUR == flag && (_myptr-p) + o > T::_len)
  {
    cerr << "IWString_Data_Source_MMAP::seekg:beyond end of file " << o << " file size " << T::_len << endl;
    return 0;
  }

  if (SEEK_CUR == flag)
    _myptr = _myptr + o;
  else
    _myptr = p + o;

  return 1;
}

template <typename T>
size_t
IW_Storage_Reader<T>::tellg() const
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  return _myptr - T::_p;
}

template <typename T>
int
IW_Storage_Reader<T>::echo_records (std::ostream & output, int necho)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const_IWSubstring buffer;
  for (int i = 0; i < necho; ++i)
  {
    if (! _next_record(buffer))
      return 0;

    output << buffer << '\n';
  }

  return necho;
}

template <typename T>
int
IW_Storage_Reader<T>::echo_records (IWString_and_File_Descriptor & output, int necho)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const_IWSubstring buffer;
  for (int i = 0; i < necho; ++i)
  {
    if (! _next_record(buffer))
      return 0;

    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return necho;
}

template <typename T>
int
IW_Storage_Reader<T>::echo (IWString_and_File_Descriptor & output, size_t necho)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const int chunk_size = 8192;

  const auto o = tellg();

  const unsigned char * e = T::_p + T::_len;

  if (_myptr + necho > e)
  {
    cerr << "IW_Storage_Reader::echo:cannot echo " << necho << " bytes\n";
    return 0;
  }

  while (necho > 0)
  {
    int t = chunk_size;
    if (_myptr + t > e)
      t = e - _myptr;

    output.strncat(reinterpret_cast<const char *>(_myptr), t);
    output.write_if_buffer_holds_more_than(chunk_size);
    _myptr += t;
  }

  return 1;
}

template <typename T>
int
IW_Storage_Reader<T>::echo (std::ostream & output, size_t nbytes)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const auto o = tellg();

  const unsigned char * e = T::_p + T::_len;

  if (_myptr + nbytes > e)
  {
    cerr << "IW_Storage_Reader::echo:cannot echo " << nbytes << " bytes\n";
    return 0;
  }

  output.write(reinterpret_cast<const char *>(_myptr), static_cast<std::streamsize>(nbytes));

  _myptr += nbytes;

  return 1;

}

template <typename T>
int
IW_Storage_Reader<T>::count_records_starting_with(const const_IWSubstring & s)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const unsigned char * msave = _myptr;
  const int lsave = _lines_read;

  const_IWSubstring buffer;
  int rc = 0;
  while (_next_record(buffer))
  {
    if (buffer.starts_with(s))
      rc++;
  }

  _myptr = msave;
  _lines_read = lsave;

  return rc;
}

template <typename T>
int
IW_Storage_Reader<T>::grep (IW_Regular_Expression & rx)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const unsigned char * msave = _myptr;
  const int lsave = _lines_read;

  const_IWSubstring buffer;
  int rc = 0;
  while (_next_record(buffer))
  {
    if (rx.matches(buffer))
      rc++;
  }

  _myptr = msave;
  _lines_read = lsave;

  return rc;
}

template <typename T>
int
IW_Storage_Reader<T>::skip_records (const int nskip)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const_IWSubstring buffer;

  for (int i = 0; i < nskip; ++i)
  {
    if (! _next_record(buffer))
      return 0;
  }

  return 1;
}

template <typename T>
int
IW_Storage_Reader<T>::skip_records (IW_Regular_Expression & rx, const int nskip)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const_IWSubstring buffer;

  int rc = 0;

  while (_next_record(buffer))
  {
    if (! rx.matches(buffer))
      continue;

    rc++;
    if (rc == nskip)
      return nskip;
  }

  return 0;
}

template <typename T>
int
IW_Storage_Reader<T>::skip_past (const char * pattern)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  // IL added T::
  // assert (T::ok());

  if (0 == strlen(pattern))
    return 0;

  IW_Regular_Expression regexp(pattern);

  if (_record_buffered)
  {
    const_IWSubstring buffer;
    if (! _fetch_previous_record(buffer))
      return 0;
    _record_buffered = 0;

    if ( regexp.matches(buffer))
      return 1;
  }

  const_IWSubstring buffer;

  int rc = 0;

  while (_next_record(buffer))
  {
    rc++;
    if (regexp.matches(buffer))
      return rc;
  }

  return 0;
}

template <typename T>
int
IW_Storage_Reader<T>::skip_to (const char * pattern)
{
  const int rc = skip_past(pattern);

  if (0 == rc)
    return 0;

  _record_buffered = 1;

  return rc;
}

template <typename T>
int
IW_Storage_Reader<T>::most_recent_record (IWString & s)
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const_IWSubstring buffer;

  if (! _fetch_previous_record(buffer))
    return 0;

  s = buffer;

  return 1;
}

template <typename T>
size_t
IW_Storage_Reader<T>::copy_raw_bytes (void * dest, const size_t nbytes) const
{
  if (! T::_good || nullptr == T::_p)
    return 0;

  const unsigned char * e = T::_p + T::_len;

  if (_myptr + nbytes > e)
  {
    cerr << "IW_Storage_Reader::copy_raw_bytes:cannot copy " << nbytes << "\n";
    return 0;
  }

  std::copy_n(_myptr, nbytes, reinterpret_cast<unsigned char *>(dest));

  return nbytes;
}

#endif

#endif
