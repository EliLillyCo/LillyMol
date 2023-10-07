#include "Foundational/iwmisc/iwre2.h"

#include "iwstring_string_data_source.h"

using std::cerr;

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

    while (zend >= _iptr && isspace(_src[zend])) {
      --zend;
    }

    buffer.set(_src + _iptr, zend - _iptr + 1);
  }
  else
    buffer.set(_src + _iptr, newline - _iptr);

  _iptr = newline + 1;

  _lines_read++;

  _most_recent_record = buffer;

  return 1;
}

template int String_Data_Source::next_record(const_IWSubstring &);
template int String_Data_Source::next_record(IWString &);

int
String_Data_Source::push_record()
{
  if (0 == _iptr)
  {
    cerr << "String_Data_Source::push_record:at beginning\n";
    return 0;
  }

  assert ('\n' == _src[_iptr - 1]);

  _iptr = _iptr - 2;    // skip back past newline

  for ( ; _iptr > 0; --_iptr) {
    if ('\n' == _src[_iptr])
    {
      _iptr++;
      return 1;
    }
  }

  return 1;
}

int
String_Data_Source::seekg(off_t o)
{
  const size_t s = strlen(_src);

  if (o > static_cast<off_t>(s))
  {
    cerr << "String_Data_Source::seekg:cannot seek to " << o << ", len " << strlen(_src) << '\n';
    return 0;
  }

  _iptr = o;

  return 1;
}

int
String_Data_Source::good() const {
  return 1;
}

int
String_Data_Source::records_remaining() const {
  int result = 0;
  for (int myptr = _iptr; _src[myptr] != '\0'; ++myptr) {
    if (_src[myptr] == '\n') {
      ++result;
    }
  }

  return result;
}

#ifdef TEMPLATES_DID_NOT_WORK_WITH_CNAME
template <typename T>
int
String_Data_Source::most_recent_record(T& buffer) const {
  buffer = _most_recent_record;
  return 1;
}

template int String_Data_Source::most_recent_record<IWString>(IWString&) const;
template int String_Data_Source::most_recent_record<const_IWSubstring>(const_IWSubstring&) const;
#endif

int
String_Data_Source::most_recent_record(IWString& buffer) const {
  buffer = _most_recent_record;
  return 1;
}

int
String_Data_Source::most_recent_record(const_IWSubstring& buffer) const {
  buffer = _most_recent_record;
  return 1;
}

// Note that _most_recent_record is not saved, so after skipping
// records, most_recent_record() will return the last record
// skipped.
int
String_Data_Source::skip_records(int nskip) {
  const_IWSubstring buffer;
  for (int i = 0; i < nskip; ++i) {
    if (! next_record(buffer)) {
      return 0;
    }
  }

  return 1;
}

int
String_Data_Source::skip_records(re2::RE2& rx, int nskip) {
  int found = 0;
  const_IWSubstring buffer;
  while(true) {
    if (! next_record(buffer)) {
      return 0;
    }
    if (iwre2::RE2PartialMatch(buffer, rx)) {
      ++found;
      if (found == nskip) {
        return 1;
      }
    }
  }

  return 0;
}
