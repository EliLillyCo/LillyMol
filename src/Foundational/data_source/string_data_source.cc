#include <stdlib.h>
#include <iostream>

#include "string_data_source.h"

void
String_Data_Source::_default_values()
{
  _lines_read = 0;

  _iptr = 0;

  _strip_trailing_blanks = 0;
  _dos = 0;

  return;
}

String_Data_Source::String_Data_Source(const char * s)
{
  assert (nullptr != s);

  _default_values();

  _src = s;

  return;
}

String_Data_Source::String_Data_Source(const char * s,
                                       int notused)
{
  assert (nullptr != s);

  _default_values();

  _src = s;

  return;
}

/*
  Really not sure what to do about good()
*/

int
String_Data_Source::good() const
{
  return 1;
}

template <typename T>
int
String_Data_Source::next_record(T & buffer)
{
  if ('\0' == _src[_iptr])
    return 0;

  int newline = _iptr;

  while ('\n' != _src[newline])
  {
    newline++;
  }

  if (_strip_trailing_blanks)
  {
    int zend = newline - 1;

    while (isspace(_src[zend]) && zend >= _iptr)
      zend--;

    buffer.set(_src + _iptr, zend - _iptr);
  }
  else
    buffer.set(_src + _iptr, newline - _iptr);

  _iptr = newline + 1;

  _lines_read++;

  return 1;
}

template int String_Data_Source::next_record(const_IWSubstring &);
template int String_Data_Source::next_record(IWString &);

int
String_Data_Source::push_record()
{
  if (0 == _iptr)
  {
    std::cerr << "String_Data_Source::push_record:at beginning\n";
    return 0;
  }

  assert ('\n' == _src[_iptr] - 1);

  _iptr = _iptr - 2;    // skip back past newline

  while (1)
  {
    if (0 == _iptr)
      return 1;

    if ('\n' == _src[_iptr])
    {
      _iptr++;
      return 1;
    }

    _iptr--;
  }
}

int
String_Data_Source::seekg(off_t o)
{
  const size_t s = strlen(_src);

  if (o > static_cast<off_t>(s))
  {
    std::cerr << "String_Data_Source::seekg:cannot seek to " << o << ", len " << strlen(_src) << '\n';
    return 0;
  }

  _iptr = o;

  return 1;
}

int
String_Data_Source::most_recent_record(IWString & buffer) const
{
  if ('\0' == _src[_iptr])
    return 0;

  int newline = _iptr;

  while ('\n' != _src[newline])
  {
    newline++;
  }

  if (_strip_trailing_blanks)
  {
    int zend = newline - 1;

    while (isspace(_src[zend]) && zend >= _iptr)
      zend--;

    buffer.set(_src + _iptr, zend - _iptr);
  }
  else
    buffer.set(_src + _iptr, newline - _iptr);

  return 1;
}
