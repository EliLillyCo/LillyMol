#include "Foundational/iwmisc/iwconfig.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>

#include <charconv>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

#ifdef USE_FAST_FLOAT
// Somewhat surprisingly fast_float is slower than
// what we have. That does not make sense, debug sometime.
#include "fast_float.h"
#endif

#include "iwstring.h"

using std::cerr;
using std::endl;

/*
  The philosophy of this class is that everything is inherited from
  resizable_array<char>.
*/

void
IWString::_default_values()
{
}

IWString::IWString()
{
  _default_values();

  return;
}

IWString::IWString(const char * s)
{
  _default_values();

  if (nullptr == s)
    return;

  IWString::strncpy(s, static_cast<int>(::strlen(s)));

  return;
}

IWString::IWString(const IWString & s)
{
  _default_values();

  IWString::strncpy(s._things, s._number_elements);

  return;
}

IWString::IWString(const const_IWSubstring & s)
{
  _default_values();

  IWString::strncpy(s._data, s._nchars);

  return;
}

IWString::IWString(char c)
{
  _default_values();

  resizable_array<char>::add(c);

  return;
}

IWString::IWString(const char * s, int n)
{
  _default_values();

  IWString::strncpy(s, n);

  return;
}

IWString::IWString(IWString&& rhs) {
  _default_values();
  _number_elements = rhs._number_elements;
  _elements_allocated = rhs._elements_allocated;
  _things = rhs._things;

  rhs._number_elements = 0;
  rhs._elements_allocated = 0;
  rhs._things = nullptr;

  return;
}

int
IWString::ok() const
{
  if (! resizable_array<char>::ok())
    return 0;

  return 1;
}

/*
  Jul 99, changed the definition so that IWString and const_IWSubstring have
  the same signature for this function
*/

int 
IWString::set(const char * s, int n)
{
  return IWString::strncpy(s, n);
}

int 
IWString::set(const char * s, size_t n)
{
  return IWString::strncpy(s, static_cast<int>(n));
}

int 
IWString::set(const char * s, u_int32_t n)
{
  return IWString::strncpy(s, static_cast<int>(n));
}

int
IWString::set_and_assume_ownership(char * s, int n)
{
  if (_things)
    delete [] _things;

  _things = s;
  _elements_allocated = n;
  _number_elements = n;

  assert(ok());

  return 1;
}

void
IWString::strip_leading_blanks()
{
  int i;
  for (i = 0; i < _number_elements; i++)
  {
    if (! isspace(_things[i]))
      break;
  }

  if (0 == i)     // no leading spaces
    return;

  _number_elements -= i;

  for (int j = 0; j < _number_elements; j++)
  {
    _things[j] = _things[i + j];
  }

  _things[_number_elements] = '\0';

  return;
}

int
const_IWSubstring::debug_print(std::ostream & os) const
{
  os << "const_IWSubstring::debug_print: " << _nchars << " at " << std::hex << &(_data[0]) << " " << _data << endl;

  return os.good();
}

void
const_IWSubstring::strip_leading_blanks()
{
  while (_nchars && isspace(_data[0]))
  {
    _data++;
    _nchars--;
  }

  if (0 == _nchars)
    _data = nullptr;

  return;
}

void
IWString::strip_trailing_blanks()
{
  if (0 == _number_elements)
    return;

  for (int i = _number_elements - 1; i >= 0; i--)
  {
    if (! isspace(_things[i]))
    {
      _number_elements = i + 1;
      if (_number_elements < _elements_allocated)
        _things[_number_elements] = '\0';

      return;
    }
  }

  _number_elements = 0;
  _things[0] = '\0';

  return;
}

void
const_IWSubstring::strip_trailing_blanks()
{
  while (_nchars && isspace(_data[_nchars - 1]))
  {
    _nchars--;
  }

  if (0 == _nchars)
    _data = nullptr;

  return;
}

/*
  Note that we don't increase the number of elements. This is really
  just so that callers can safely use rawchars().
*/

int
IWString::null_terminate()
{
  if (_number_elements == _elements_allocated)
    resize(_elements_allocated + 1);

//if ('\0' == _things[_number_elements])     drives valgrind NUTS
//  return 0;

  _things[_number_elements] = '\0';

  return 1;
}

/*
  This is pretty awful...
*/

void
IWString::_const_appearing_null_terminate() const
{
  ((IWString *)this)->null_terminate();

  return;
}

/*
  Some operation has occurred which shortened the string representation.
  Find the '\0' and adjust _number_elements;
*/

/*void
IWString::_recompute_length()
{
  for (int i = 0; i < _number_elements; i++)
  {
    if ('\0' == _things[i])
    {
      _number_elements = i;
      return;
    }
  }

// If we pass to here, we did not find a null terminator. For now, let's assume that's OK

  return;
}*/

const char *
IWString::chars()
{
  assert(ok());

  null_terminate();

  return _things;
}

const char *
IWString::null_terminated_chars()
{
  assert(ok());

  null_terminate();

  return _things;
}
void
IWString::chop(int nchars)
{
  assert (nchars >= 0 && nchars <= _number_elements);
  if (0 == nchars)
    return;

  int i = _number_elements - nchars;
  _number_elements = i;
  _things[i] = '\0';    

  return;
}


void
IWString::iwtruncate(int nchars)
{
  assert (nchars >= 0 && nchars <= _number_elements);
  if (nchars == _number_elements)    // already the appropriate length
    return;

  _number_elements = nchars;
  if (_number_elements < _elements_allocated)
    _things[_number_elements] = '\0';

  return;
}

static int
common_truncate_at_first(const char * s, int & nchars, char c)
{
  for (int i = 0; i < nchars; i++)
  {
    if (s[i] == c)
    {
      nchars = i;
      return 1;
    }
  }

  return 0;
}

int
IWString::truncate_at_first(char c)
{
  return common_truncate_at_first(_things, _number_elements, c);
}

int
const_IWSubstring::truncate_at_first(char c)
{
  return common_truncate_at_first(_data, _nchars, c);
}

static int 
common_truncate_at_last(const char * s, int & nchars, char c)
{
  for (int i = nchars - 1; i >= 0; i--)
  {
    if (c == s[i])
    {
      nchars = i;
      return 1;
    }
  }

  return 0;      // did not find the character
}

int
IWString::truncate_at_last(char c)
{
  return common_truncate_at_last(_things, _number_elements, c);
}

int
const_IWSubstring::truncate_at_last(char c)
{        
  return common_truncate_at_last(_data, _nchars, c);
}

void
IWString::remove_leading_chars(int nremove)
{
  assert (nremove >= 0 && nremove <= _number_elements);

  for (int i = nremove; i < _number_elements; i++)
  {
    _things[i - nremove] = _things[i];
  }

  _number_elements = _number_elements - nremove;

  return;
}

/*
  this variant keeps the string the same length and fills the empty spaces with 'PAD'
*/

void
IWString::remove_leading_chars(int nremove, char pad)
{
  assert (nremove >= 0 && nremove <= _number_elements);

  for (int i = nremove; i < _number_elements; i++)
  {
    _things[i - nremove] = _things[i];
  }

  for (int i = _number_elements - nremove; i < _number_elements; i++)
  {
    _things[i] = pad;
  }

  return;
}


static int
common_count_leading_chars(const char * s, int nchars, char toremove)
{
  int rc = 0;

  for (int i = 0; i < nchars; i++)
  {
    if (toremove != s[i])
      break;

    rc++;
  }

  return rc;
}

int
IWString::remove_leading_chars(char toremove)
{
  assert (ok());

  int chars_removed = common_count_leading_chars(_things, _number_elements, toremove);

  if (0 == chars_removed)
    return 0;

  _number_elements = _number_elements - chars_removed;

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = _things[i + chars_removed];
  }

  return chars_removed;
}

int
const_IWSubstring::remove_leading_chars(char toremove)
{
  int chars_removed = common_count_leading_chars(_data, _nchars, toremove);

  if (0 == chars_removed)
    return 0;

  _data += chars_removed;

  _nchars = _nchars - chars_removed;

  return chars_removed;
}

void
const_IWSubstring::remove_line_terminators()
{
  if (0 >= _nchars)
    return;


  // make sure the buffer does not end with a \r or \n
  int i = _nchars - 1;
  while(_data[i] == '\r' || _data[i] == '\n')
  {
     _nchars = i;
    --i;
    if (i <= 0)
      break;
  }

  return;
}


int
IWString::shift(int nshift, char pad)
{
  if (nshift < 0)
  {
    remove_leading_chars(-nshift, pad);
    return 1;
  }

  if (0 == nshift)
    return 0;

  if (_number_elements + nshift > _elements_allocated)
    resize(_number_elements + nshift + 1);      // +1 just in case a null termination follows

  for (int i = _number_elements + nshift - 1; i >= nshift; i--)
  {
    _things[i] = _things[i - nshift];
  }

  for (int i = 0; i < nshift; i++)
  {
    _things[i] = pad;
  }

  _number_elements += nshift;

  return 1;
}

int
IWString::remove_all_these(const char * remove_these)
{
  int rc = 0;
  int j = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    char c = _things[i];
    if (NULL == ::strchr(remove_these, c))
    {
      _things[j] = c;
      j++;
    }
    else
      rc++;
  }

  _number_elements = j;

  if (rc)
    _things[j] = '\0';

  return rc;
}

int
IWString::compress_blanks()
{
  int istop = _number_elements;

  _number_elements = 0;
  for (int i = 0; i < istop; i++)
  {
    _things[_number_elements] = _things[i];
    _number_elements++;

    if (! isspace(_things[i]))
      continue;

//  We found a space. Skip over any other spaces which follow

    do
    {
      i++;

      if (i == istop)
        return 1;
    }
    while (isspace(_things[i]));

    _things[_number_elements] = _things[i];
    _number_elements++;
  }

  return 1;
}

IWString &
IWString::operator=(const char * rhs)
{
  int nchars = static_cast<int>(::strlen(rhs));

  if (0 == nchars)
  {
    resize(0);
    return *this;
  }

  resize(nchars + 1);
  assert(_things);

  IW_STRCPY(_things, rhs);
  _number_elements = nchars;

  assert(ok());

  return *this;
}

const_IWSubstring &
const_IWSubstring::operator=(const char & rhs)
{
  _data = & rhs;

  _nchars = 1;

  return *this;
}

const_IWSubstring &
const_IWSubstring::operator=(const char * rhs)
{
  _data = rhs;

  if (nullptr == rhs)
    _nchars = 0;
  else
    _nchars = static_cast<int>(::strlen(rhs));

  return *this;
}

const_IWSubstring &
const_IWSubstring::operator=(const IWString & rhs)
{
  _nchars = rhs.length();

  _data = rhs.rawchars();

  return *this;
}

void
IWString::operator +=(const char * rhs)
{
  int nchars = static_cast<int>(::strlen(rhs));
  if (0 == nchars)
    return;

  (void) add(rhs, nchars);

  return;
}

void
IWString::operator +=(const IWString & rhs)
{
  if (0 == rhs._number_elements)
    return;

  add(rhs._things, rhs._number_elements);
  
  return;
}

void
IWString::operator +=(const const_IWSubstring & rhs)
{
  if (0 == rhs.nchars())
    return;

  (void) add(rhs.rawchars(), rhs.nchars());

  return;
}

void
IWString::to_lowercase()
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = tolower(_things[i]);
  }

  return;
}

void
IWString::to_uppercase()
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = toupper(_things[i]);
  }

  return;
}

/*
  We need to turn strcmp into a strncmp because neither of the strings may
  be null terminated.
*/

static int
_strcmp (const char * s1, int l1, const char * s2, int l2)
{
  uint32_t ncomp = l1;
  if (l2 < l1) {
    ncomp = l2;
  }

  return ::strncmp(s1, s2, ncomp);
}

static int
_strcasecmp(const char * s1,
             int l1,
             const char * s2,
             int l2)
{
  if (l1 < l2)    // s1 is shorter
  {
    int rc = iwstrncasecmp(s1, s2, l1);
    if (0 != rc)
      return rc;

    return -1;     // s1 is less than s2
  }
  else if (l1 == l2)
    return iwstrncasecmp(s1, s2, l1);
  else    // s2 is shorter
  {
    int rc = iwstrncasecmp(s1, s2, l2);
    if (0 != rc)
      return rc;

    return 1;    // s1 is greater than s2
  }
}

int
IWString::strcmp(const IWString & rhs) const
{
  return _strcmp(_things, _number_elements, rhs._things, rhs._number_elements);
}

int
IWString::strcmp(const char * rhs) const
{
  return _strcmp(_things, _number_elements, rhs, static_cast<int>(::strlen(rhs)));
}

int
const_IWSubstring::strcmp(const const_IWSubstring & rhs) const
{
  return _strcmp(_data, _nchars, rhs._data, rhs._nchars);
}       

int
IWString::strcasecmp(const IWString & rhs) const
{
  return _strcasecmp(_things, _number_elements, rhs._things, rhs._number_elements);
}

int
IWString::strcasecmp(const char * rhs) const
{
  return _strcasecmp(_things, _number_elements, rhs, static_cast<int>(strlen(rhs)));
}

int
const_IWSubstring::strcasecmp(const const_IWSubstring & rhs) const
{
  return _strcasecmp(_data, _nchars, rhs._data, rhs._nchars);
}       

/*
  Because neither of the string may be null terminated, we need to be
  careful in doing a strncmp
*/

static int
_strncmp (const char * s1, int len1, const char * s2, int len2,
          int chars_to_compare)
{
#ifdef DEBUG_STRNCMP_Q
  cerr << "_strncmp '";
  cerr.write(s1, len1);
  cerr << "' " << len1 << " with '";
  cerr.write(s2, len2);
  cerr << "' " << len2 << endl;
#endif

  if (len1 < len2)
    return -1;

  if (len1 > len2)
    return 1;

// lengths are equal

  if (chars_to_compare > len1)    // 
  {
    cerr << "_strncmp:comparison length " << chars_to_compare << " longer than string length " << len1 << " impossible\n";
    return -1;    // an arbitrary choice
  }

  return strncmp(s1, s2, chars_to_compare);
}

int
IWString::strncmp(const IWString & rhs, int chars_to_compare) const
{
  return _strncmp(_things, _number_elements, rhs._things, rhs._number_elements, chars_to_compare);
}

int
IWString::strncmp(const char * rhs, int chars_to_compare) const
{
  return _strncmp(_things, _number_elements, rhs, chars_to_compare, chars_to_compare);
}

int
const_IWSubstring::strncmp(const char * rhs, int chars_to_compare) const
{
  return _strncmp(_data, _nchars, rhs, chars_to_compare, chars_to_compare);
}

static int
_common_matches_ignore_case (const char * needle,
                             size_t len_needle,
                             const char * haystack,
                             size_t len_haystack,
                             int offset)
{
  if (offset + len_needle > len_haystack)
    return 0;

  return 0 == iwstrncasecmp(haystack + offset, needle, static_cast<int>(len_needle));
}

int
IWString::matches_ignore_case(const char * s,
                               size_t lens,
                               int offset) const
{
  return _common_matches_ignore_case(s, lens, _things, _number_elements, offset);
}

int
IWString::matches_ignore_case(const IWString & s,
                               int offset) const
{
  return _common_matches_ignore_case(s._things, s._number_elements, _things, _number_elements, offset);
}

int
const_IWSubstring::matches_ignore_case(const char * s,
                       size_t lens,
                       int offset) const
{
  return _common_matches_ignore_case(s, lens, _data, _nchars, offset);
}

static int
common_equals_ignore_case (const char * s1,
                           size_t lens1,
                           const char * s2,
                           size_t lens2)
{
  if (lens1 != lens2)
    return 0;

  return 0 == iwstrncasecmp(s1, s2, static_cast<int>(lens1));
}

int
const_IWSubstring::equals_ignore_case(const char * s, size_t lens) const
{
  return common_equals_ignore_case(_data, _nchars, s, lens);
}

int
IWString::equals_ignore_case(const char * s, size_t lens) const
{
  return common_equals_ignore_case(_things, _number_elements, s, lens);
}

int
const_IWSubstring::equals_ignore_case(const IWString & s) const
{
  return common_equals_ignore_case(s.rawchars(), s.length(), _data, _nchars);
}

static int
do_find (const char * haystack, int size_of_haystack, 
         const char * needle, int size_of_needle)
{
  assert (needle);

  if (0 == size_of_haystack)
    return -1;

  if (size_of_needle > size_of_haystack)
    return -1;

  int last_i = size_of_haystack - size_of_needle;

  for (int i = 0; i <= last_i; i++)
  {
    if (haystack[i] != needle[0])
      continue;

    if (0 == ::strncmp(haystack + i, needle, size_of_needle))
      return i;
/*  int found_match = 1;
    for (int j = 1; j < size_of_needle; j++)
    {
      if (haystack[i + j] != needle[j])
      {
        found_match = 0;
        break;
      }
    }

    if (found_match)
      return i;*/
  }

  return -1;
}

int
IWString::find(const char * needle) const
{
  return do_find(_things, _number_elements, needle, static_cast<int>(::strlen(needle)));
}

int
const_IWSubstring::find(const char * needle) const
{
  return do_find(_data, _nchars, needle, static_cast<int>(::strlen(needle)));
}

int
IWString::find(const IWString & needle) const
{
  return do_find(_things, _number_elements, needle._things, needle._number_elements);
}

int
const_IWSubstring::find(const IWString & needle) const
{
  return do_find(_data, _nchars, needle.rawchars(), needle.nchars());
}

int
IWString::find(const const_IWSubstring & needle) const
{
  return do_find(_things, _number_elements, needle._data, needle._nchars);
}

int
const_IWSubstring::find(const const_IWSubstring & needle) const
{
  return do_find(_data, _nchars, needle._data, needle._nchars);
}

/*
  Unsigned int conversion
  This is not robust
*/

template <typename T>
int
string_class_is_unsigned_int_4 (const char * s, int nchars, T & result)
{
  if (nchars < 10)
    ;
  else if (nchars > 10)
    return 0;
  else if (s[0] > '4')    // very rough check, fix this mess sometime
    return 0;

  if (nchars < 10)
  {
    result = 0;
    for (int i = 0; i < nchars; i++)
    {
      int j = s[i] - '0';
      if (j < 0)
        return 0;
      if (j > 9)
        return 0;
  
      result = 10 * result + j;
    }

    return 1;
  }

  T tmp;
  if (! string_class_is_unsigned_int_4(s, nchars - 1, tmp))
    return 0;

  result = tmp * 10 + s[9] - '0';

  if (result > tmp)
    return 1;

  return 0;
}

/*
  max unsigned long is 18446744073709551615
*/

template <typename T>
int
string_class_is_unsigned_int_8 (const char * s, int nchars, T & result)
{
  if (nchars > 20)
    return 0;

  if (nchars < 20)
    ;
  else if (s[0] > '1')    // very rough check, fix this mess sometime
    return 0;

  result = 0;
  if (nchars < 20)
  {
    for (int i = 0; i < nchars; i++)
    {
      int j = s[i] - '0';
      if (j < 0)
        return 0;
      if (j > 9)
        return 0;

      result = 10 * result + j;
    }

    return 1;
  }

  T tmp;
  if (! string_class_is_unsigned_int_8(s, nchars - 1, tmp))
    return 0;

  result = tmp * 10 + s[19] - '0';

  if (result > tmp)
    return 1;

//cerr << "Overflowed int 8\n";

  return 0;
}

//#ifdef __GNUG__
template int string_class_is_unsigned_int_4 (const char *, int, unsigned int &);
template int string_class_is_unsigned_int_8 (const char *, int, unsigned long int &);
template int string_class_is_unsigned_int_8 (const char *, int, unsigned long long int &);
template int string_class_is_unsigned_int_4<unsigned long>(char const*, int, unsigned long&);
//#endif

/*
  Pretty dumb int conversion. Note that this does not deal adequately with
  large negative values and such.
*/

template <typename T>
int
string_class_is_int_4 (const char * s, int nchars, T & result)
{
  if (0 == nchars)
    return 0;

  int sign = 0;
  if ('-' == *s)
  {
    sign = 1;
    s++;
    nchars--;
  }
  else if ('+' == *s)
  {
    s++;
    nchars--;
  }

  if (0 == nchars)
    return 0;

  if (nchars > 10) // must be erroneous
    return 0;

  if (nchars < 10)    // the most common case
  {
    result = 0;
    for (int i = 0; i < nchars; i++)
    {
      if (s[i] >= '0' && s[i] <= '9')
        result = 10 * result + (s[i] - '0');
      else
        return 0;
    }

    if (sign)
      result = - result;

    return 1;
  }

// Now the case of 10 digits. Let the system function do this

  if (*s > '3')     // int max is 2147483647, so any 10 digit number starting with 3 is invalid
    return 0;

//cerr << "Converting long string\n";

  unsigned int tmp;
  if (! string_class_is_unsigned_int_4(s, nchars, tmp))
    return 0;

  if (0 == sign)
  {
    if (tmp > static_cast<unsigned int>(std::numeric_limits<T>::max()))
      return 0;
    result = static_cast<T>(tmp);
  }
  else if (tmp > 2147483648u)
    return 0;
  else
    result = - tmp;

  return 1;

#ifdef OLD_STUFFQ
// Strings must be null terminated for the system functions

  char buffer[13];
  IW_STRNCPY (buffer, s, nchars);
  buffer[nchars] = '\0';

  char *p;
  unsigned long tmp = strtoul (buffer, &p, 10);

  if ('\0' != *p)
    return 0;

  if (tmp > static_cast<unsigned long>(std::numeric_limits<int>::max()))
    return 0;

  result = T (tmp);
  if (sign)
    result = -result;

  return 1;
#endif
}

/*
  This is not robust, should be fixed sometime.
  Highest 8 byte signed int is 9223372036854775807
*/

template <typename T>
int
string_class_is_int_8 (const char * s, int nchars, T & result)
{
  if (0 == nchars)
    return 0;

  int sign = 0;
  if ('-' == *s)
  {
    sign = 1;
    s++;
    nchars--;
  }
  else if ('+' == *s)
  {
    s++;
    nchars--;
  }

  if (0 == nchars)
    return 0;

  if (nchars > 19) // must be erroneous
    return 0;

  if (nchars < 19)    // the most common case
  {
    result = 0;
    for (int i = 0; i < nchars; i++)
    {
      if (s[i] >= '0' && s[i] <= '9')
        result = 10 * result + (s[i] - '0');
      else
        return 0;
    }

    if (sign)
      result = - result;

    return 1;
  }

// Now the case of 19 digits. just need to be careful. Use an unsigned value to get the number
// then check to see if it is OK
  
  unsigned long int tmp;
  if (! string_class_is_unsigned_int_8(s, nchars, tmp))
    return 0;

  if (0 == sign)
  {
    if (tmp > static_cast<unsigned long int>(std::numeric_limits<long int>::max()))
      return 0;
  }
  else if (tmp > 9223372036854775808U)
    return 0;
  
  result = static_cast<T>(tmp);
  if (sign)
    result = - result;

  return 1;
}

//#if defined (__GNUG__) || defined (__SUNPRO_CC)
template int string_class_is_int_4 (const char *, int, int &);
template int string_class_is_int_8 (const char *, int, long &);
template int string_class_is_int_8 (const char *, int, long long &);
template int string_class_is_int_4<long>(char const*, int, long&);
template int string_class_is_unsigned_int_8<unsigned int>(char const*, int, unsigned int&);
template int string_class_is_int_8<int>(char const*, int, int&);
//#endif

#define MAX_FRACTIONAL_DIGITS 15

static double reciprocal_power_of_10[MAX_FRACTIONAL_DIGITS] = {1.0,
                               0.1,
                               0.01,
                               0.001,
                               0.0001,
                               0.00001,
                               0.000001,
                               0.0000001,
                               0.00000001,
                               0.000000001,
                               0.0000000001,
                               0.00000000001,
                               0.000000000001,
                               0.0000000000001,
                               0.00000000000001};

#define MAX_LEADING_DIGITS 15

static double power_of_10[MAX_LEADING_DIGITS] = {1.0,
                               10.0,
                               100.0,
                               1000.0,
                               10000.0,
                               100000.0,
                               1000000.0,
                               10000000.0,
                               100000000.0,
                               1000000000.0,
                               10000000000.0,
                               100000000000.0,
                               1000000000000.0,
                               10000000000000.0,
                               100000000000000.0};

static double double_digit [] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

static int
iw_parse_numeric_pos (const char * s, int nchars, double & result)
{
  if (1 == nchars)
  {
    int i = *s - '0';
    if (i >= 0 && i <= 9)
    {
      result = double_digit[i];
      return 1;
    }

    return 0;
  }

  while (nchars && isdigit (*s))
  {
    result = 10.0 * result + static_cast<double>(*s - '0');
    s++;
    nchars--;
  }

  if (0 == nchars)
    return 1;

// In parsing the fraction, we treat leading 0's specially. That way we can
// parse '0.0000000001234567891234'

  if ('.' == *s)
  {
    s++;
    nchars--;

    int leading_zero = 0;
    while (nchars && '0' == *s)
    {
      leading_zero++;
      s++;
      nchars--;
    }

    if (0 == nchars)
      return 1;

    double fraction = 0.0;

    for (int i = 1; nchars && i < MAX_FRACTIONAL_DIGITS && isdigit(*s); i++, nchars--, s++)
    {
      fraction = fraction + reciprocal_power_of_10[i] * static_cast<double>(*s - '0');
    }

    if (0 == leading_zero)
      ;
    else if (leading_zero <= 10)
      fraction = fraction * reciprocal_power_of_10[leading_zero];
    else
      fraction = fraction / pow(10.0, static_cast<double>(leading_zero));

    result = result + fraction;

//  ignore any extra digits

    while (nchars && isdigit(*s))
    {
      s++;
      nchars--;
    }

    if (0 == nchars)
      return 1;
  }

  if ('e' == *s)
    ;
  else if ('E' == *s)
    ;
  else 
    return 0;

  s++;
  nchars--;

  if (0 == nchars)     // cannot end with 'e'
    return 0;

  if (nchars > 4)      // max 4 digits in the exponent
    return 0;

  int exponent_sign;
  if ('-' == *s)
  {
    exponent_sign = -1;
    s++;
    nchars--;
    if (0 == nchars)
      return 0;
  }
  else if ('+' == *s)
  {
    exponent_sign = 1;
    s++;
    nchars--;
    if (0 == nchars)
      return 0;
  }
  else
    exponent_sign = 1;

  int exponent = 0;

  while (nchars && isdigit(*s))
  {
    exponent = 10 * exponent + (*s - '0');
    s++;
    nchars--;
  }

  if (nchars)     // some strange character out there
    return 0;

  if (exponent_sign > 0 && exponent < 10)
  {
    result = result * power_of_10[exponent];
    return 1;
  }
  if (exponent_sign < 0 && exponent < 10)
  {
    result = result * reciprocal_power_of_10[exponent];
    return 1;
  }

  if (exponent_sign > 0)
    result = result * pow(10.0, static_cast<double>(exponent));
  else
    result = result * pow(10.0, static_cast<double>(-exponent));

  return 1;
}

static int
iw_parse_numeric(const char * s, int nchars, double & result)
{
  if (0 == nchars)
    return 0;

  int zsign;

  if ('-' == *s)
  {
    zsign = -1;
    s++;
    nchars--;

    if (0 == nchars)
      return 0;
  }
  else
    zsign = 1;

  result = 0.0;

  int rc = iw_parse_numeric_pos(s, nchars, result);

  if (0 == rc)
    return 0;

  if (zsign < 0)
    result = -result;

  return rc;
}

int
IWString::is_int(int & result) const
{
  return string_class_is_int_4(_things, _number_elements, result);
}

int
IWString::numeric_value(unsigned char & result) const
{
  int tmp;
  if (! string_class_is_int_4(_things, _number_elements, tmp))
    return 0;

  if (tmp < 0 || tmp > std::numeric_limits<unsigned char>::max())
    return 0;

  result = static_cast<unsigned char>(tmp);

  return 1;
}

int
IWString::numeric_value(unsigned short & result) const
{
  int tmp;
  if (! string_class_is_int_4(_things, _number_elements, tmp))
    return 0;

  if (tmp < 0 || tmp > std::numeric_limits<unsigned short>::max())
    return 0;

  result = static_cast<unsigned short>(tmp);

  return 1;
}

int
IWString::numeric_value(int & result) const
{
  if (4 == sizeof(result))
    return string_class_is_int_4(_things, _number_elements, result);
  else
    return string_class_is_int_8(_things, _number_elements, result);
}

int
IWString::numeric_value(long & result) const
{
  if (4 == sizeof(result))
    return string_class_is_int_4(_things, _number_elements, result);
  else
    return string_class_is_int_8(_things, _number_elements, result);
}

int
IWString::numeric_value(long long & result) const
{
  return string_class_is_int_8(_things, _number_elements, result);
}

int
IWString::numeric_value(unsigned int & result) const
{
  if (4 == sizeof(result))
    return string_class_is_unsigned_int_4(_things, _number_elements, result);
  else
    return string_class_is_unsigned_int_8(_things, _number_elements, result);
}

int
IWString::numeric_value(unsigned long & result) const
{
  if (4 == sizeof(result))
    return string_class_is_unsigned_int_4(_things, _number_elements, result);
  else
    return string_class_is_unsigned_int_8(_things, _number_elements, result);
}

/*
  Somewhat of a kludge for the hex value
*/

int
IWString::is_hex(unsigned int & result) const
{
  if (_number_elements < 3)
    return 0;

  _const_appearing_null_terminate();


  char *c;
  long tmp = ::strtol(_things, &c, 16);

  if (c == _things)
    return 0;

  if ('\0' == *c)
    ;     // good
  else if (! isspace(*c))
    return 0;

  result = tmp;

  return 1;
}

int
IWString::is_double(double & result) const
{
  _const_appearing_null_terminate();

  return ::is_double(_things, &result);
}

int
IWString::numeric_value(float & result) const
{
#ifdef USE_FAST_FLOAT
  auto conv = fast_float::from_chars(_things, _things + _number_elements, result);
  if (conv.ec != std::errc()) {
    return 0;
  }
  return 1;
#endif
  double tmp;

  if (! numeric_value(tmp))
    return 0;

  // Check for overflow.
  static const double maxfloat = static_cast<double>(std::numeric_limits<float>::max());

  if (fabs(tmp) > maxfloat)
    return 0;

  result = tmp;

  return 1;
}

int
IWString::numeric_value(double & result) const
{
#ifdef USE_FAST_FLOAT
  auto conv = fast_float::from_chars(_things, _things + _number_elements, result);
  if (conv.ec != std::errc()) {
    return 0;
  }
  return 1;
#endif
  return iw_parse_numeric(_things, _number_elements, result);
}

static int
common_nwords (const char * s,
               int nchars)
{
  int rc = 0;
  int in_word = 0;
  for (int i = 0; i < nchars; i++)
  {
    if (' ' == s[i])      // changed 6 Sept from isspace (s[i])
    {
      if (in_word)
        in_word = 0;
    }
    else
    {
      if (! in_word)
      {
        rc++;
        in_word = 1;
      }
    }
  }

  return rc;
}

int
IWString::nwords() const
{
  return common_nwords (_things, _number_elements);
}

int
const_IWSubstring::nwords() const
{
  return common_nwords (_data, _nchars);
}

#if defined (__SUNPRO_CC)
#else
static
#endif
int
common_nwords_with_separator (const char * s,
               int nchars,
               char word_separator)
{
  int rc = 0;
  int in_word = 0;
  for (int i = 0; i < nchars; i++)
  {
    if (s[i] == word_separator)
    {
      if (in_word)
        in_word = 0;
    }
    else
    {
      if (! in_word)
      {
        rc++;
        in_word = 1;
      }
    }
  }

  return rc;
}

int
IWString::nwords(char word_separator) const
{
  return common_nwords_with_separator (_things, _number_elements, word_separator);
}

int
const_IWSubstring::nwords(char word_separator) const
{
  return common_nwords_with_separator (_data, _nchars, word_separator);
}

/*
  Isn't there some way we can share the code between the word locating functions - maybe
  some clever templating?
*/

static int
locate_whitespace_word_boundaries (const char * s, int nchars,
                               int which_word,
                               int & word_start, int & word_stop)
{
  int in_word = 0;
  int word_count = 0;
  for (int i = 0; i < nchars; i++)
  {
    if (isspace (s[i]))
    {
      if (in_word)
        in_word = 0;
    }
    else
    {
      if (! in_word)
      {
        if (word_count == which_word)
        {
          word_start = i;
          for (int j = i + 1; j < nchars; j++)
          {
            if (isspace (s[j]))
            {
              word_stop = j - 1;
              return 1;
            }
          }

//        If we get to here, we ran into the end of the string

          word_stop = nchars - 1;

          return 1;
        }
        word_count++;
      }
      in_word = 1;
    }
  }

  return 0;
}

/*
  Someone has entered something like word -1
*/

static int
common_locate_word_boundaries_rev (const char * s, int nchars,
                                   int which_word,
                                   char word_separator,
                                   int & word_start, int & word_stop)
{
  assert (which_word < 0);

  which_word = - which_word;   // convert to positive
  which_word--;

  int in_word = 0;
  int word_count = 0;

  for (int i = nchars - 1; i >= 0; i--)
  {
    if (s[i] == word_separator)
    {
      if (in_word)
        in_word = 0;
    }
    else if (! in_word)   // encountered start of new word
    {
      if (word_count == which_word)
      {
        word_stop = i;
        for (int j = i - 1; j >= 0; j--)
        {
          if (s[j] == word_separator)
          {
            word_start = j + 1;
            return 1;
          }
        }

//      If we get to here, we ran into the beginning of the string

        word_start = 0;

        return 1;
      }
      word_count++;
      in_word = 1;
    }
  }

  return 0;
}

static int
common_locate_word_boundaries(const char * s, int nchars,
                              int which_word,
                              char word_separator,
                              int & word_start, int & word_stop)
{
  if (which_word < 0)
    return common_locate_word_boundaries_rev(s, nchars, which_word, word_separator, word_start, word_stop);

  int in_word = 0;
  int word_count = 0;
  for (int i = 0; i < nchars; i++)
  {
    if (s[i] == word_separator)
    {
      if (in_word)
        in_word = 0;
    }
    else
    {
      if (! in_word)
      {
        if (word_count == which_word)
        {
          word_start = i;
          for (int j = i + 1; j < nchars; j++)
          {
            if (s[j] == word_separator)
            {
              word_stop = j - 1;
              return 1;
            }
          }

//        If we get to here, we ran into the end of the string

          word_stop = nchars - 1;

          return 1;
        }
        word_count++;
      }
      in_word = 1;
    }
  }

  return 0;
}

int
IWString::_locate_word_boundaries(int which_word, char word_separator,
                                   int & word_start, int & word_stop) const
{
  int in_word = 0;
  int word_count = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i] == word_separator)
    {
      if (in_word)
        in_word = 0;
    }
    else
    {
      if (! in_word)
      {
        word_count++;
        if (word_count == which_word)
        {
          word_start = i;
          for (int j = i + 1; j < _number_elements; j++)
          {
            if (_things[j] == word_separator)
            {
              word_stop = j - 1;
              return 1;
            }
          }

//        If we get to here, we ran into the end of the string

          word_stop = _number_elements - 1;

          return 1;
        }
      }
      in_word = 1;
    }
  }

  return 0;
}

const_IWSubstring
IWString::word(int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
    return const_IWSubstring ("", 0);

  return const_IWSubstring (&(_things[word_start]), word_stop - word_start + 1);
}

int
IWString::word(int which_word, IWString & result, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
  {
    result = "";
    return 0;
  }

  result = from_to (word_start, word_stop);

  return 1;
}

int
IWString::word(int which_word, const_IWSubstring & result, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
  {
    result = "";
    return 0;
  }

  from_to (word_start, word_stop, result);

  return 1;
}

int
IWString::wordindex(int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
    return -1;

  return word_start;
}

const_IWSubstring
const_IWSubstring::word(int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_data, _nchars,
                 which_word, word_separator, word_start, word_stop))
    return const_IWSubstring ("", 0);

  return const_IWSubstring (&(_data[word_start]), word_stop - word_start + 1);
}

int
const_IWSubstring::word(int which_word, IWString & result, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_data, _nchars,
                 which_word, word_separator, word_start, word_stop))
  {
    result = "";
    return 0;
  }

//from_to(word_start, word_stop, result);
  result = from_to (word_start, word_stop);

  return 1;
}

int
const_IWSubstring::word(int which_word, const_IWSubstring & result, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_data, _nchars,
                 which_word, word_separator, word_start, word_stop))
  {
    result = "";
    return 0;
  }

  from_to (word_start, word_stop, result);

  return 1;
}

int
const_IWSubstring::whitespace_delimited_word(int which_word, IWString & result) const
{
  int word_start, word_stop;

  if (! locate_whitespace_word_boundaries (_data, _nchars,
                 which_word, word_start, word_stop))
  {
    result = "";
    return 0;
  }

  from_to (word_start, word_stop, result);

  return 1;
}

int
const_IWSubstring::wordindex(int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_data, _nchars,
                 which_word, word_separator, word_start, word_stop))
    return -1;

  return word_start;
}

static int
_locate_word_beginnings (const char * s,
                         int nchars,
                         resizable_array<int> & word_beginnings,
                         char word_separator)
{
  assert (0 == word_beginnings.number_elements ());

  int in_word = 0;

  for (int i = 0; i < nchars; i++)
  {
    if (s[i] == word_separator)
    {
      if (in_word)
        in_word = 0;
    }
    else
    {
      if (! in_word)
        word_beginnings.add (i);

      in_word = 1;
    }
  }

  return word_beginnings.number_elements ();
}

int
const_IWSubstring::locate_word_beginnings(resizable_array<int> & word_beginnings,
                                           char separator) const
{
  return _locate_word_beginnings(_data, _nchars, word_beginnings, separator);
}

int
IWString::locate_word_beginnings(resizable_array<int> & word_beginnings,
                                           char separator) const
{
  return _locate_word_beginnings(_things, _number_elements, word_beginnings, separator);
}

static int
common_locate_word_beginnings_single_delimiter(const char * s,
                                               int nchars,
                                               resizable_array<int> & word_beginnings,
                                               char separator)
{
  word_beginnings.resize_keep_storage(0);

  if (0 == nchars)
    return 0;

  word_beginnings.add(0);   // always

  for (int i = 0; i < nchars; i++)
  {
    if (separator != s[i])
      continue;

//  if (nchars - 1 != i)       removed April 2010, what follows is an empty token
    word_beginnings.add(i + 1);
  }

  return word_beginnings.number_elements();
}

int
const_IWSubstring::locate_word_beginnings_single_delimiter(resizable_array<int> & word_beginnings,
                        char separator) const
{
  return common_locate_word_beginnings_single_delimiter(_data, _nchars, word_beginnings, separator);
}

int
IWString::locate_word_beginnings_single_delimiter(resizable_array<int> & word_beginnings,
                        char separator) const
{
  return common_locate_word_beginnings_single_delimiter(_things, _number_elements, word_beginnings, separator);
}

/*
  The code for nextword is common between IWString and constIWSubstring.
*/

static int
internal_nextword(const char * s, int nchars, int & i,
                  char separator,
                  int & istart, int & istop)
{
  if (i >= nchars)
    return 0;

// strip over any whitespace we may be in

  while (separator == s[i])
  {
    i++;
    if (i >= nchars)
      return 0;
  }

// We are in a 'word'

  istart = i;

  i++;

  for ( ; i < nchars; ++i)
  {
    if (separator == s[i])
    {
      istop = i-1;
      return 1;
    }
  }

  istop = nchars - 1;
  return 1;
}
  
static int
internal_nextword(IWString & result, int & i,
                  char separator,
                  const char * s, int nchars)
{
  result.resize_keep_storage(0);

  int istart, istop;

  if (! internal_nextword(s, nchars, i, separator, istart, istop))
    return 0;

//cerr << "After STRING determination, istart = " << istart << " istop = " << istop << endl;

  result.resize(istop - istart + 1);

  result.set(s + istart, istop - istart + 1);

  return 1;
}

#if defined (__SUNPRO_CC)
#else
static
#endif
int
internal_nextword(const_IWSubstring & result, int & i,
                  char separator,
                  const char * s, int nchars)
{
  result.set(nullptr, 0);

  int istart, istop;

  if (! internal_nextword(s, nchars, i, separator, istart, istop))
    return 0;

//cerr << "After determination, istart = " << istart << " istop = " << istop << endl;

  result.set(&s[istart], istop - istart + 1);

  return 1;
}

/*
  Used when scanning for words backwards in a string.
  Note that we must scan back to the start of the word, and then
  add all the letters up to where we started.
*/

static int
internal_prevword(IWString & result, int & i,
                  char separator,
                  const char * s)
{
  if (i <= 0)
    return 0;

  result.resize_keep_storage (0);

// strip over any whitespace we may be in

  while (separator == s[i])
  {
    i--;
    if (i <= 0)
      return 0;
  }

  int end_of_word = i;
  int start_of_word = 0;

  i--;

  for ( ; i >= 0; --i)
  {
    if (separator == s[i])
    {
      start_of_word = i + 1;
      break;
    }
  }


// Now build the word from i to end_of_word

  for (int j = start_of_word; j <= end_of_word; j++)
  {
    result += s[j];
  }

  return 1;
}

int
IWString::nextword(IWString & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword(result, i, separator, _things, _number_elements);
}

int
const_IWSubstring::nextword(IWString & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword(result, i, separator, _data, _nchars);
}

int
IWString::nextword(const_IWSubstring & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword(result, i, separator, _things, _number_elements);
}

int
const_IWSubstring::nextword(const_IWSubstring & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword(result, i, separator, _data, _nchars);
}

int
IWString::prevword(IWString & result, int & i, char separator) const
{
  if (i <= 0)
    return 0;

  assert (i < _number_elements);

  return internal_prevword(result, i, separator, _things);
}

int
const_IWSubstring::prevword(IWString & result, int & i, char separator) const
{
  if (i <= 0)
    return 0;

  assert (i < _nchars);

  return internal_prevword(result, i, separator, _data);
}

int 
IWString::nwords_single_delimiter(char separator) const
{
  if (0 == _number_elements)
    return 0;

  return ccount(separator) + 1;
}

int 
const_IWSubstring::nwords_single_delimiter(char separator) const
{
  if (0 == _nchars)
    return 0;

  return ccount(separator) + 1;
}

const_IWSubstring
IWString::substr(int istart, int nchars) const
{
  assert (istart >= 0 && istart < _number_elements);
  if (nchars < 0)
    nchars = _number_elements - istart;

  assert (nchars >= 0 && istart + nchars <= _number_elements);

  return const_IWSubstring(&(_things[istart]), nchars);
}

const_IWSubstring
IWString::from_to(int istart, int istop) const
{
  assert (istart >= 0);
  assert (istop >= istart && istop < _number_elements);

  return const_IWSubstring(&(_things[istart]), istop - istart + 1);
}

const_IWSubstring
const_IWSubstring::from_to(int istart, int istop) const
{
  assert (istart >= 0);
  assert (istop >= istart && istop < _nchars);

  return const_IWSubstring(&(_data[istart]), istop - istart + 1);
}

/*
  THIS DOES NOW WORK, from_to(int,int) GETS CALLED INSTEAD.
  return all the characters from a given position up to a given
  character.
  right now, behaviour when the char is not found is undefined.
  Should we return the entire rest of the string?
  Should we return an empty string?
*/

const_IWSubstring
IWString::from_to(int istart, const char * s) const
{
  assert (istart >= 0 && istart < _number_elements);

  if (0 == _number_elements)
    return const_IWSubstring(nullptr, 0);

  int lens = static_cast<int>(strlen(s));

  if (0 == lens)
    return const_IWSubstring(nullptr, 0);

  int i = do_find(_things, _number_elements, s, lens);

  if (i < 0)
    return const_IWSubstring(nullptr, 0);

  return const_IWSubstring(&_things[istart], i - istart + 1);
}

void
IWString::from_to(int istart, int istop, char * destination) const
{
  assert (istart >= 0);
  assert (istop >= istart && istart < _number_elements);

  int nchars = istop - istart + 1;
  memcpy(destination, &(_things[istart]), static_cast<size_t>(nchars));
  destination[nchars] = '\0';

  return;
}

void
IWString::from_to(int istart, int istop, const_IWSubstring & destination) const
{
  assert (ok_index(istart));
  assert (istop >= istart && istop < _number_elements);

  destination._data = &_things[istart];
  destination._nchars = istop - istart + 1;

  return;
}

void
IWString::from_to(int istart, int istop, IWString & destination) const
{
  assert (ok_index(istart));
  assert (istop >= istart && istop < _number_elements);

  destination.strncpy(&_things[istart], istop - istart + 1);

  return;
}

void
const_IWSubstring::from_to(int istart, int istop, const_IWSubstring & destination) const
{
  assert (istart >= 0 && istart < _nchars);
  assert (istop >= istart && istop < _nchars);

  destination._data = &_data[istart];
  destination._nchars = istop - istart + 1;

  return;
}

void
const_IWSubstring::from_to(int istart, int istop, IWString & destination) const
{
  assert (istart >= 0 && istart < _nchars);
  assert (istop >= istart && istop < _nchars);

  destination.strncpy(&_data[istart], istop - istart + 1);

  return;
}


void
const_IWSubstring::from_to(int istart, IWString & zresult) const
{
  assert (istart >= 0 && istart < _nchars);

  zresult.strncpy(_data + istart, _nchars - istart);

  return;
}

IWString &
IWString::operator=(const const_IWSubstring & rhs)
{
  if (_elements_allocated < rhs._nchars)
    resize(rhs._nchars + 1);

  IW_STRNCPY(_things, rhs._data, rhs._nchars);
  _number_elements = rhs._nchars;

  return *this;
}

IWString &
IWString::operator=(char c)
{
  if (_elements_allocated < 1)
    resize(2);

  _things[0] = c;
  if (_elements_allocated > 1)
    _things[1] = '\0';

  _number_elements = 1;

  return *this;
}

const_IWSubstring::const_IWSubstring()
{
  _data = nullptr;
  _nchars = 0;

  return;
}

const_IWSubstring::const_IWSubstring(const char * data, int nchars) :
          _data(data), _nchars(nchars)
{
  return;
}

const_IWSubstring::const_IWSubstring(const char * data) :
          _data(data)
{
  if (nullptr != data)
    _nchars = static_cast<int>(::strlen(data));
  else
    _nchars = 0;

  return;
}


const_IWSubstring::const_IWSubstring(const IWString & s) :
          _data(s._things), _nchars(s._number_elements)
{
  return;
}

const_IWSubstring::~const_IWSubstring()
{
  _nchars = 0;
}

const_IWSubstring &
const_IWSubstring::operator =(const const_IWSubstring & other)
{
  _data = other._data;
  _nchars = other._nchars;

  return *this;
}

IWString &
IWString::operator =(const IWString & other)
{
  if (_elements_allocated < other._number_elements)
    resize(other._number_elements);

#ifdef _WIN32
  IW_STRNCPY(_things, other._things, other._number_elements);
#else
  ::memcpy(_things, other._things, other._number_elements);
#endif

//for (int i = 0; i < other._number_elements; i++)
//  _things[i] = other._things[i];

  _number_elements = other._number_elements;

  return *this;
}

/*
  Maybe there is some way to make this more elegant...
*/

IWString &
IWString::operator =(IWString && other)
{
  resizable_array<char> & t = *this;
  resizable_array<char> & o = other;
  t = std::move(o);

//cerr << " After move '" << other << "'\n";

  return *this;
}

int
const_IWSubstring::copy_to_char_array(char * s) const
{
  IW_STRNCPY(s, _data, _nchars);
  s[_nchars] = '\0';

  return _nchars;
}

int
IWString::copy_to_char_array(char * s) const
{
  if (_number_elements)
    IW_STRNCPY(s, _things, _number_elements);

  s[_number_elements] = '\0';

  return _number_elements;
}

const_IWSubstring
const_IWSubstring::substr(int cstart, int nchars) const
{
  if (nchars < 0)
    nchars = _nchars - cstart;

  assert (nchars >= 0 && nchars < _nchars);

  return const_IWSubstring(&(_data[cstart]), nchars);
}

int
const_IWSubstring::starts_with(char s) const
{
  if (0 == _nchars)
    return 0;

  return s == _data[0];
}

int
const_IWSubstring::starts_with(const char * s) const
{
  int lens = static_cast<int>(strlen(s));
  if (_nchars < lens)
    return 0;

  return 0 == ::strncmp(_data, s, lens);
}

int
const_IWSubstring::starts_with(const const_IWSubstring & s) const
{
  if (s._nchars > _nchars)
    return 0;

  return 0 == ::strncmp(_data, s._data, s._nchars);
}

int
const_IWSubstring::starts_with(const IWString & s) const
{
  if (s._number_elements > _nchars)
    return 0;

  return 0 == ::strncmp(_data, s._things, s._number_elements);
}

int
const_IWSubstring::ends_with(char s) const
{
  if (0 == _nchars)
    return 0;

  return s == _data[_nchars - 1];
}

int
const_IWSubstring::ends_with(const char * s) const
{
  int lens = static_cast<int>(strlen(s));
  if (_nchars < lens)
    return 0;

  return 0 == ::strncmp(_data + _nchars - lens, s, lens);
}

int
const_IWSubstring::ends_with(const const_IWSubstring & s) const
{
  if (s.length() > _nchars)
    return 0;

  return 0 == ::strncmp(_data + _nchars - s._nchars, s._data, s._nchars);
}

int
const_IWSubstring::ends_with(const IWString & s) const
{
  if (s.length() > _nchars)
    return 0;

  return 0 == ::strncmp(_data + _nchars - s.length(), s.rawdata(), s.length());
}

int
const_IWSubstring::remove_leading_chars(int nremove)
{
  assert (nremove > 0 && nremove <= _nchars);

  _data += nremove;
  _nchars -= nremove;

  return _nchars;
}

/*
  We return the number of characters to be chopped
  note that when asked to chop the first word, we
  also get rid of separators after the first word
*/

static int
common_remove_leading_words (const char * s, int nchars,
                             int nremove, char word_separator)
{
  int word_start, word_stop;

// Most of the time, we will strip off everything before the first
// word to be retained

  if (common_locate_word_boundaries(s, nchars,
                 nremove, word_separator, word_start, word_stop))

    return word_start;

// Well, that didn't work. Why? Well, if the string contains just
// one word and we were asked to remove the first word, we should
// be able to do that

  if (nremove > 1)
    return 0;

  return nchars;
}

int
const_IWSubstring::remove_leading_words(int nremove, char separator)
{
  int rc = common_remove_leading_words(_data, _nchars, nremove, separator);

  if (rc <= 0)
    return rc;

  remove_leading_chars(rc);

  return rc;
}

int
IWString::remove_leading_words(int nremove, char separator)
{
  int rc = common_remove_leading_words(_things, _number_elements, nremove, separator);

  if (rc <= 0)
    return rc;

  remove_leading_chars(rc);

  return rc;
}

int
IWString::remove_word(int which_word, char word_separator)
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries(_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
  {
    return 0;
  }

// Special case for spaces - remove trailing spaces from words other than
// the last, leading spaces from the last word

  if (' ' == word_separator)
  {
    if (word_stop == _number_elements - 1 && word_start > 0)
    {
      while (word_start > 0 && ' ' == _things[word_start - 1])
      {
        word_start--;
      }
    }
    else
    {
      while (word_stop < _number_elements - 1 && ' ' == _things[word_stop + 1])
      {
        word_stop++;
      }
    }
  }

  erase(word_start, word_stop);

  return 1;
}

int
const_IWSubstring::chop(int nchop)
{
  assert (nchop > 0 && nchop <= _nchars);

  return _nchars -= nchop;
}

int
const_IWSubstring::iwtruncate(int n)
{
  assert (n >= 0 && n <= _nchars);

  return _nchars = n;
}

/*
  Cannot think of an easy way to do this with a lib routine.
  We do it the dumb (and broken for MAX_INT) way
*/

int
const_IWSubstring::is_int(int & result) const
{
  return string_class_is_int_4(_data, _nchars, result);
}

int
const_IWSubstring::numeric_value(float & result) const
{
  double tmp;

  if (! numeric_value(tmp))
    return 0;

  static const double maxfloat = static_cast<double>(std::numeric_limits<float>::max());

  if (fabs(tmp) > maxfloat)
    return 0;

  result = static_cast<float>(tmp);

  return 1;
}

int
const_IWSubstring::numeric_value(double & result) const
{
  if (0 == _nchars)
    return 0;

  if (_nchars > 25)
    return 0;

  return iw_parse_numeric(_data, _nchars, result);

//char buffer[30];

//strncpy (buffer, _data, _nchars);
//buffer[_nchars] = '\0';

//char * endptr;

//result = strtod (buffer, &endptr);
//if ('\0' == *endptr)
//  return 1;

//return 0;
}

int
const_IWSubstring::numeric_value(int & result) const
{
  return string_class_is_int_4(_data, _nchars, result);
}

int
const_IWSubstring::numeric_value(unsigned char & result) const
{
  int tmp;
  if (! string_class_is_int_4(_data, _nchars, tmp))
    return 0;

  if (tmp < 0 || tmp > std::numeric_limits<unsigned char>::max())
    return 0;

  result = static_cast<unsigned char>(tmp);

  return 1;
}

int
const_IWSubstring::numeric_value(unsigned short & result) const
{
  int tmp;
  if (! string_class_is_int_4(_data, _nchars, tmp))
    return 0;

  if (tmp < 0 || tmp > std::numeric_limits<unsigned short>::max())
    return 0;

  result = static_cast<unsigned char>(tmp);

  return 1;
}

int
const_IWSubstring::numeric_value(unsigned int & result) const
{
  return string_class_is_unsigned_int_4(_data, _nchars, result);
}

int
const_IWSubstring::numeric_value(long & result) const
{
  if (4 == sizeof(result))
    return string_class_is_int_4(_data, _nchars, result);
  else
    return string_class_is_int_8(_data, _nchars, result);
}

int
const_IWSubstring::numeric_value(long long & result) const
{
  return string_class_is_int_8(_data, _nchars, result);
}

int
const_IWSubstring::numeric_value(unsigned long & result) const
{
  return string_class_is_unsigned_int_8(_data, _nchars, result);
}

int
const_IWSubstring::numeric_value(unsigned long long & result) const
{
  return string_class_is_unsigned_int_8(_data, _nchars, result);
}

int
const_IWSubstring::index(const char c) const
{
  for (int i = 0; i < _nchars; i++)
  {
    if (c == _data[i])
      return i;
  }

  return -1;
}

int
const_IWSubstring::rindex(const char c) const
{
  for (int i = _nchars - 1; i >= 0; i--)
  {
    if (c == _data[i])
      return i;
  }

  return -1;
}

static int
do_next (const char * s, int nchars, const char c, int & istart)
{
  while (istart < nchars)
  {
    if (c == s[istart])
      return 1;

    istart++;
  }

  return 0;
}

int
const_IWSubstring::next(const char c, int & istart) const
{
  assert (istart >= 0 && istart < _nchars);

  return do_next(_data, _nchars, c, istart);
}

/*ostream &
operator << (ostream & os, const const_IWSubstring & ss)
{
  assert (os.good ());

  os.write (ss.rawchars (), ss.nchars ());

  return os;
}*/

/*ostream &
operator << (ostream & os, const const_IWSubstring * ss)
{
  assert (os.good ());

  os.write (ss->rawchars(), ss->nchars());

  return os;
}*/

void
append_int (IWString & s, int n)
{
  char buffer[24];

  (void) IW_SPRINTF(buffer, "%d", n);

  s += buffer;

  return;
}

static void
_cat (const char * lhs, int l1, const char * rhs, int l2,
     IWString & result)
{
  result.resize(l1 + l2 + 1);   // leave room for newline

  result.add(lhs, l1);
  result.add(rhs, l2);

  return;
}

IWString
operator + (const char * lhs, const IWString & rhs)
{
  int l1 = static_cast<int>(strlen(lhs));

  IWString rc;

  _cat(lhs, l1, rhs._things, rhs._number_elements, rc);

  return rc;
}

IWString
operator + (const IWString & lhs, const char * rhs)
{
  IWString rc;

  int l2 = static_cast<int>(strlen(rhs));

  _cat(lhs._things, lhs._number_elements, rhs, l2, rc);

  return rc;
}

IWString
operator + (const IWString & lhs, const IWString & rhs)
{
  IWString rc;

  _cat(lhs._things, lhs._number_elements, rhs._things, rhs._number_elements, rc);

  return rc;
}

void
append_digit (IWString & s, int digit)
{
  char buffer[24];

  IW_SPRINTF(buffer, "%d", digit);

  s += buffer;

  return;
}

void
IWString::append_number(int znumber)
{
  if (znumber < 0)
  {
    resizable_array<char>::add('-');
    znumber = - znumber;    // WRONG, won't work for largest negative number
  }

  _append_int_form(znumber);

  return;
}

void
IWString::append_number(unsigned int znumber)
{
  _append_int_form(znumber);

  return;
}

void
IWString::append_number(long znumber)
{
  if (znumber < 0)
  {
    resizable_array<char>::add('-');
    znumber = -znumber;        // breaks for largest negative number
  }

  _append_int_form(znumber);

  return;
}

void
IWString::append_number(long long znumber)
{
  if (znumber < 0)
  {
    resizable_array<char>::add('-');
    znumber = -znumber;                 // breaks for largest negative number
  }

  _append_int_form(znumber);

  return;
}
void
IWString::append_number(unsigned long long znumber)
{
  _append_int_form(znumber);

  return;
}

void
IWString::append_number(unsigned long znumber)
{
  _append_int_form(znumber);
  return;
}

static int float_precision = 7;

void
set_default_iwstring_float_concatenation_precision (int s)
{
  float_precision = s;
}

void
IWString::append_number(float f, int fprecision)
{
  char buffer[32];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  gcvt(static_cast<double>(f), fprecision, buffer);
#pragma GCC diagnostic pop
  resizable_array<char>::add(buffer, static_cast<int>(::strlen(buffer)));
  return;
}

void
IWString::append_number(float f) {
  append_number(f, float_precision);
}

static int double_precision = 10;

void
set_default_iwstring_double_concatenation_precision (int s)
{
  double_precision = s;
}

void
IWString::append_number(double d)
{
  append_number(d, double_precision);
}

void
IWString::append_number(double d, int dprecision)
{
  char buffer[32];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  gcvt(d, dprecision, buffer);
#pragma GCC diagnostic pop
  resizable_array<char>::add(buffer, static_cast<int>(::strlen(buffer)));
  return;
}

void
IWString::append_number(float f, const char * fformat)
{
  char buffer[32];

  int nchars = IW_SPRINTF(buffer, fformat, f);

  resizable_array<char>::add(buffer, nchars);

  return;
}

const_IWSubstring
substr (const IWString &s, int cstart, int nchars)
{ 
  if (nchars < 0)
    nchars = s.number_elements() - cstart;

  return s.substr(cstart, nchars);
}

const_IWSubstring
substr (const const_IWSubstring & s, int cstart, int nchars)
{
  if (nchars < 0)
    nchars = s._nchars - cstart;

  return s.substr(cstart, nchars);
}

const_IWSubstring
from_to (const IWString & s, int cstart, int cstop)
{
  return s.from_to(cstart, cstop);
}

/*
  Note that this is a non-robust version. Things larger than MAX_INT will
  just die.
*/

static int
_is_int (const char * s, int nchars, int & result)
{
// skip over leading blanks

  int istart = 0;
  while (isspace(s[istart]) && istart < nchars)
  {
    istart++;
  }

  if (nchars == istart)    // blank token
    return 0;

  int tmp = 0;
  int zsign = 1;
  if ('-' == s[istart])
  {
    zsign = -1;
    istart = 1;
  }
  else if ('+' == s[istart])
    istart = 1;

  for (int i = istart; i < nchars; i++)
  {
    int j = s[i] - '0';
    if (j < 0 || j > 9)
      return 0;

    tmp = 10 * tmp + j;
    if (tmp < 0)    // feeble attempt at detecting overflow
      return 0;
  }

  result = zsign * tmp;
  return 1;
}

int
is_int (const const_IWSubstring & s, int & result)
{
  return _is_int(s.rawchars(), s.nchars(), result);
}

int
change_suffix (IWString & fname, const IWString & new_suffix)
{
  int len_fname = fname.length();

  if (0 == len_fname)
    return 0;

  int last_period = -1;
  for (int i = 0; i < len_fname; i++)
  {
    if ('.' == fname[i])
      last_period = i;
  }

  if (last_period < 0)
    fname += '.';
  else
    fname.shorten(last_period + 1);

  fname += new_suffix;

  return 1;
}

/*
  Shorten is a convenient way to shorten a string. If you are
  interested in storage, use resize ();
*/

int
IWString::shorten(int new_size)
{
  assert (new_size >= 0 && new_size <= _number_elements);

  _number_elements = new_size;

  return 1;
}

int
IWString::starts_with(char s) const
{
  if (0 == _number_elements)
    return 0;

  return (s == _things[0]);
}

int
IWString::starts_with(const char * s, int lens) const
{
  if (0 == _number_elements)
    return 0;

  // This seems very ugly and fragile, TODO:ianwatson get rid of this?
  if (lens < 0) {
    lens = static_cast<int>(::strlen(s));
  }

  if (_number_elements < lens) {
    return 0;
  }

  return (0 == ::strncmp(_things, s, lens));
}

int
IWString::starts_with(const IWString & s) const
{
  int lens = s._number_elements;

  if (lens > _number_elements)
    return 0;

  return (0 == ::strncmp(_things, s._things, lens));
}

int
IWString::starts_with(const const_IWSubstring & s) const
{
  int lens = s._nchars;

  if (lens > _number_elements)
    return 0;

  return (0 == ::strncmp(_things, s._data, lens));
}

int
IWString::translate(char cfrom, char cto)
{
  assert (ok());

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (cfrom == _things[i])
    {
      _things[i] = cto;
      rc++;
    }
  }

  return rc;
}

int
IWString::ends_with(char s) const
{
  if (0 == _number_elements)
    return 0;

  return (s == _things[_number_elements - 1]);
}

int
IWString::ends_with(const char * s, int lens) const
{
  assert(s);

  if (lens < 0)
    lens = static_cast<int>(::strlen(s));

  if (lens > _number_elements)
    return 0;

  if (0 == _number_elements)     // catch the case of both blank
    return 0;

  int i = _number_elements - lens;

#ifdef MANUALLY_CHECK_ENDS_WITH
  for (int q = 0; q < lens; ++q)
  {
    if (s[q] == _things[i+q])
      cerr << "Character " << q << " '" << _things[i+q] << "' the same\n";
    else
      cerr << "Character " << q << " '" << _things[i+q] << "' different\n";
  }
#endif

  return (0 == ::strncmp(s, _things + i, static_cast<uint32_t>(lens)));
}

int
IWString::ends_with(const IWString & s) const
{
  if (s._number_elements > _number_elements)
    return 0;

  if (0 == _number_elements)     // catch the case of both being blank
    return 0;

  int i = _number_elements - s._number_elements;
  return (0 == ::strncmp(s._things, &(_things[i]), s._number_elements));
}

int
IWString::ends_with(const const_IWSubstring & s) const
{
  return ends_with(s.rawchars(), s.length());
}

int
IWString::looks_like(const char * s, int min_chars_needed_for_match) const
{
  assert (ok());
  assert (s);
  assert (min_chars_needed_for_match);

  if (_number_elements < min_chars_needed_for_match)
    return 0;

  int lens = static_cast<int>(::strlen(s));

  if (_number_elements > lens) 
    return 0;

  return 0 == ::strncmp(_things, s, _number_elements);
}

/*
  For now, we limit ourselves to the $VARIABLE at the start of OLD_NAME
*/

int
expand_environment_variables (const char * old_name, IWString & expanded_name)
{
  if ('$' != *old_name)
    return 0;

  char buffer[256];
  old_name++;
  int iptr = 0;
  while (isalnum(*old_name) || '_'==*old_name)
  {
    buffer[iptr] = *old_name;
    old_name++;
    iptr++;
  }
  buffer[iptr] = '\0';

  const char * env = getenv(buffer);

  if (nullptr == env)
  {
    cerr << "expand_environment_variables: no env value for '" << buffer << "'\n";
    return 0;
  }

//cerr << "Value for environment variable '" << buffer << "' is '" << env << "'\n";

  expanded_name = env;
  expanded_name += old_name;

//cerr << "After expansion, name is '" << expanded_name << "'\n";
  return 1;
}

int
IWString::insert(char c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  if (_number_elements == _elements_allocated)
    resize(_number_elements + 1);

  for (int i = _number_elements; i > pos; i--)
  {
    _things[i] = _things[i - 1];
  }

  _things[pos] = c;

  _number_elements++;

  return 1;
}

int
IWString::insert(const char * c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  int lens = static_cast<int>(::strlen(c));

  return _insert(pos, c, lens);
}

int
IWString::insert(const IWString & c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  return _insert(pos, c._things, c._number_elements);
}

int
IWString::insert(const const_IWSubstring & c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  return _insert(pos, c._data, c._nchars);
}

int
IWString::_insert(int pos, const char * s, int lens)
{
  if (_number_elements + lens > _elements_allocated)
    resize(_number_elements + lens);

  for (int i = _number_elements + lens - 1; i >= pos + lens; i--)
  {
    _things[i] = _things[i - lens];
  }

  for (int i = 0; i < lens; i++)
  {
    _things[pos + i] = s[i];
  }

  _number_elements += lens;

  return 1;
}


int
IWString::overwrite(char c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  _things[pos] = c;

  return 1;
}

int
IWString::overwrite(const char * s, int pos)
{
  int lens = static_cast<int>(::strlen(s));

  return _overwrite(s, lens, pos);
}

int
IWString::overwrite(const const_IWSubstring & s, int pos)
{
  return _overwrite(s.rawchars(), s.length(), pos);
}

int
IWString::overwrite(const IWString & s, int pos)
{
  return _overwrite(s.rawchars(), s.length(), pos);
}

int
IWString::_overwrite(const char * s, int nchars, int pos)
{
  assert (pos >= 0 && pos < _number_elements);
  assert (pos + nchars - 1 < _number_elements);

  for (int i = 0; i < nchars; i++)
  {
    _things[pos + i] = s[i];
  }

  return 1;
}

/*
  Locate the next occurrence of C. Increment istart to allow repeated calls.
*/

int
IWString::next(const char c, int & istart) const
{
  assert (ok_index(istart));

  return do_next(_things, _number_elements, c, istart);
}

/*
  This looks to be exactly the same as translate
*/

int
IWString::gsub(char cfrom, char cto, int how_many)
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (cfrom == _things[i])
    {
      _things[i] = cto;
      rc++;
      if (how_many > 0 && rc >= how_many)
        return rc;
    }
  }

  return rc;
}

static int
_internal_ccount (const char * haystack,
                  int nchars,
                  char needle)
{
  int rc = 0;
  for (int i = 0; i < nchars; i++)
  {
    if (haystack[i] == needle)
      rc++;
  }

  return rc;
}

int
const_IWSubstring::ccount(char c) const
{
  return _internal_ccount(_data, _nchars, c);
}

int
IWString::ccount(char c) const
{
  return _internal_ccount(_things, _number_elements, c);
}

int
IWString::strncpy(const char * s, int nchars)
{
  if (_elements_allocated < nchars + 1)
    resize(nchars + 1);

#ifdef _WIN32
  IW_STRNCPY(_things, s, nchars);
#else
  ::memcpy(_things, s, static_cast<uint32_t>(nchars));
#endif

  _number_elements = nchars;

  return 1;
}

int
IWString::strncat(const char * s, int nchars)
{
  make_room_for_extra_items(nchars);

  ::memcpy(_things + _number_elements, s, nchars);

  for (int i = 0; i < nchars; i++)
  {
    _things[_number_elements + i] = s[i];
  }

  _number_elements += nchars;

  return _number_elements;
}

int
IWString::strncat(const IWString & s, int nchars)
{
  assert (nchars <= s.length());

  return strncat(s.rawchars(), nchars);
}

int
IWString::strncat(const const_IWSubstring & s, int nchars)
{
  assert (nchars <= s.length());

  return strncat(s.rawchars(), nchars);
}

int
IWString::append_with_spacer(const const_IWSubstring zextra, char spacer)
{
  if (0 == _number_elements)
  {
    operator= (zextra);
    return 1;
  }

  if (_elements_allocated < _number_elements + 1 + zextra.nchars())
    resize(_number_elements + 1 + zextra.nchars());

  *this << spacer << zextra;

  return _number_elements;
}

int
IWString::append_with_spacer(const const_IWSubstring zextra, const IWString & spacer)
{
  if (0 == _number_elements)
  {
    operator= (zextra);
    return 1;
  }

  make_room_for_extra_items(zextra.length() + spacer.length());

  *this << spacer << zextra;

  return _number_elements;
}

/*
  A post-fix increment operator
*/

char
const_IWSubstring::operator++(int)
{
//cerr << "Substring operator ++ postfix, " << _nchars << " characters\n";
  if (0 == _nchars)
    return '\0';

//cerr << "Currently pointing at '" << *_data << "', will return '" << (*(_data + 1)) << "'\n";
  _data++;
  _nchars--;

  if (0 == _nchars)
    return '\0';

  return *_data;
}

/*
  A pre-fix increment operator
*/

char
const_IWSubstring::operator++()
{
  if (0 == _nchars)
    return '\0';

  char rc = *_data;
  _data++;
  _nchars--;
  return rc;
}

void
const_IWSubstring::operator+=(int howfar)
{
  assert (howfar >= 0 && howfar <= _nchars);

  _data += howfar;
  _nchars -= howfar;

  return;
}

int
IWString::strspn(const char * s2)
{
  null_terminate ();

  return static_cast<int>(::strspn(_things, s2));
}

/*
  We frequently want the part of a string before a character as well
  as the part after.

  If the character is not found, we return 0
  If the character is found at the beginning or end of the string,
  we return 1 and fill in either before or after
  If it is somewhere in the middle, we return 2
*/

template <typename S>
int
common_split (const char * s, int nchars,
              S & before,
              char separator,
              S & after)
{
  before = "";
  after = "";

  int sindex = -1;
  for (int i = 0; i < nchars; i++)
  {
    if (separator == s[i])
    {
      sindex = i;
      break;
    }
  }

  if (sindex < 0)
    return 0;

  if (0 == sindex)       // before is empty
  {
    after.set(s + 1, nchars - 1);
    return 1;
  }

  if (sindex == nchars - 1)    // after is empty
  {
    before.set(s, nchars - 1);
    return 1;
  }

  before.set(s, sindex);
  after.set(s + sindex + 1, nchars - sindex - 1);

  return 2;
}

//#if defined (__GNUG__) || defined (__SUNPRO_CC)
template int common_split(const char *, int, const_IWSubstring &, char, const_IWSubstring &);
template int common_split(const char *, int, IWString &, char, IWString &);
//#endif

int
IWString::split(const_IWSubstring & before,
                 char separator,
                 const_IWSubstring & after) const
{
  return common_split(_things, _number_elements, before, separator, after);
}

int
IWString::split(IWString & before,
                 char separator,
                 IWString & after) const
{
  return common_split(_things, _number_elements, before, separator, after);
}

int
const_IWSubstring::split(const_IWSubstring & before,
                          char separator,
                          const_IWSubstring & after) const
{
  return common_split(_data, _nchars, before, separator, after);
}

int
const_IWSubstring::split(IWString & before,
                          char separator,
                          IWString & after) const
{
  return common_split(_data, _nchars, before, separator, after);
}

IWString &
operator << (IWString & s1, const IWString & s2)
{
  s1 += s2;

  return s1;
}

IWString &
operator << (IWString & s1, const const_IWSubstring & s2)
{
  s1 += s2;

  return s1;
}

IWString &
IWString::operator <<(const const_IWSubstring & s2)
{
  operator += (s2);

  return *this;
}

IWString & 
IWString::operator <<(const IWString & s2)
{
  operator += (s2);

  return *this;
}

IWString &
operator << (IWString & s1, const char * s2)
{
  s1 += s2;

  return s1;
}

IWString &
operator << (IWString & s1, char s2)
{
//s1.add(s2);
  s1 += s2;

  return s1;
}

IWString &
IWString::operator <<(char s2)
{
  operator += (s2);

  return *this;
}

IWString &
IWString::operator <<(const char * s2)
{
  add (s2, static_cast<int>(::strlen(s2)));

  return *this;
}

IWString &
operator << (IWString & s1, int s2)
{
  s1.append_number(s2);

  return s1;
}

IWString &
IWString::operator <<(int i)
{
  append_number(i);

  return *this;
}

IWString &
IWString::operator <<(long i)
{
  append_number(i);

  return *this;
}

IWString &
IWString::operator <<(long long i)
{
  append_number(i);

  return *this;
}

IWString &
IWString::operator <<(unsigned long long i)
{
  append_number(i);

  return *this;
}

IWString &
IWString::operator <<(unsigned long i)
{
  append_number(i);

  return *this;
}

IWString &
IWString::operator <<(unsigned int i)
{
  append_number(i);

  return *this;
}

IWString &
operator << (IWString & s1, float s2)
{
  s1.append_number(s2);

  return s1;
}

IWString &
IWString::operator <<(float s2)
{
  append_number(s2);

  return *this;
}

IWString &
IWString::operator <<(double s2)
{
  append_number(s2);

  return *this;
}

static int
common_split (resizable_array_p<const_IWSubstring> & tokens,
              char separator,
              const char * zdata,
              int nchars)
{
  int rc = 0;

  int i = 0;
  const_IWSubstring token;
  while (internal_nextword(token, i, separator, zdata, nchars))
  {
    const_IWSubstring * tmp = new const_IWSubstring(zdata + i - token.nchars(), token.nchars());
    tokens.add(tmp);

    rc++;
  }

  return rc;
}

static int
common_split (resizable_array_p<IWString> & tokens,
              char separator,
              const char * zdata,
              int nchars)
{
  int rc = 0;

  int i = 0;
  const_IWSubstring token;
  while (internal_nextword(token, i, separator, zdata, nchars))
  {
    IWString * tmp = new IWString(zdata + i - token.nchars(), token.nchars());
    tokens.add(tmp);

    rc++;
  }

  return rc;
}

int
IWString::split(resizable_array_p<const_IWSubstring> & tokens,
                 char separator) const
{
  tokens.resize_keep_storage(0);
  int nw = nwords(separator);
  if (tokens.elements_allocated() < nw)
    tokens.resize(nw);

  return common_split(tokens, separator, _things, _number_elements);
}

int
IWString::split(resizable_array_p<IWString> & tokens,
                 char separator) const
{
  tokens.resize_keep_storage(0);

  int nw = nwords(separator);
  if (tokens.elements_allocated() < nw)
    tokens.resize(nw);

  return common_split(tokens, separator, _things, _number_elements);
}

int
const_IWSubstring::split(resizable_array_p<const_IWSubstring> & tokens,
                          char separator) const
{
  tokens.resize_keep_storage(0);
  int nw = nwords(separator);
  if (tokens.elements_allocated() < nw)
    tokens.resize(nw);

  return common_split(tokens, separator, _data, _nchars);
}

template <typename T>
int
common_split (const char * s, int nchars,
              iwaray<T> & tokens,
              char separator)
{
  int nw = common_nwords_with_separator(s, nchars, separator);

  if (! tokens.resize(nw))
  {
    cerr << "Memory failure in split\n";
    return 0;
  }

  int rc = 0;

  int i = 0;
  const_IWSubstring token;
  while (internal_nextword(token, i, separator, s, nchars))
  {
    tokens[rc] = token;

    rc++;
  }

  assert (rc == nw);

  return rc;
}

//#if defined (__GNUG__) || defined (__SUNPRO_CC)
template int common_split(const char *, int, iwaray<IWString> &, char);
template int common_split(const char *, int, iwaray<const_IWSubstring> &, char);
//#endif

int
IWString::split(iwaray<const_IWSubstring> & tokens, char separator) const
{
  return common_split(_things, _number_elements, tokens, separator);
}

int
const_IWSubstring::split(iwaray<const_IWSubstring> & tokens, char separator) const
{
  return common_split(_data, _nchars, tokens, separator);
}

int
IWString::split(iwaray<IWString> & tokens, char separator) const
{
  return common_split(_things, _number_elements, tokens, separator);
}

static int
common_split (const char * s,
              int nchars,
              IWString * tokens,
              char separator)
{
  assert (nullptr != tokens);

  int i = 0;
  const_IWSubstring token;
  int rc = 0;

  while (internal_nextword(token, i, separator, s, nchars))
  {
    tokens[rc] = token;

    rc++;
  }

  return rc;
}

int
const_IWSubstring::split(IWString * tokens, char separator) const
{
  return common_split(_data, _nchars, tokens, separator);
}


int
IWString::split(IWString * tokens, char separator) const
{
  return common_split(_things, _number_elements, tokens, separator);
}

int
const_IWSubstring::split (std::vector<std::string>& tokens, std::string separators) const
{
  // I have given up trying to do this in IW world - std::string is easier to understand
  // create tokens for both strings and separators, i.e. break down
  // 'reactant1.qry+reactant2.qry>>product.qry' to
  // 'reactant1.qry', '+', 'reactant2.qry', '>>', 'product.qry'
  // assuming + and > are separators
  // returns the number of tokens
  std::string s(_data, _nchars);  // essential to use the length here
  tokens.clear();
  std::string::size_type last_pos = s.find_first_not_of(separators, 0);
  int cnt = 0;
  if (last_pos > 0)
  {
    tokens.push_back(s.substr(0, last_pos));
    cnt++;
  }
  std::string::size_type pos = s.find_first_of(separators, last_pos);

  while (std::string::npos != pos || std::string::npos != last_pos)
  {
    tokens.push_back(s.substr(last_pos, pos - last_pos));
    last_pos = s.find_first_not_of(separators, pos); // string token
    cnt++;
    if (std::string::npos != pos) {
      if (std::string::npos == last_pos) {  // ends with a separator
        tokens.push_back(s.substr(pos, s.size() - pos));
        cnt++;
      }
      else
      {
        tokens.push_back(s.substr(pos, last_pos - pos));  // separator token
        cnt++;
      }
    }
    pos = s.find_first_of(separators, last_pos);
  }
  return cnt;
}

int
IWString::remove_up_to_first(char target)
{
  int index_of_target = -1;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i] == target)
    {
      index_of_target = i;
      break;
    }
  }

  if (index_of_target < 0)
    return 0;

  if (index_of_target == _number_elements - 1)
  {
    resize(0);
    return index_of_target + 1;
  }

  int jptr = 0;
  for (int i = index_of_target + 1; i < _number_elements; i++)
  {
    _things[jptr] = _things[i];
    jptr++;
  }

  _number_elements = jptr;

  return index_of_target + 1;    // the number of characters removed
}

int
const_IWSubstring::remove_up_to_first(char target)
{
  int index_of_target = -1;
  for (int i = 0; i < _nchars; i++)
  {
    if (target == _data[i])
    {
      index_of_target = i;
      break;
    }
  }

  if (index_of_target < 0)
    return 0;

  _data += index_of_target + 1;

  _nchars -= (index_of_target + 1);

  return (index_of_target + 1);
}

/*
  All these operators with scalars need template specialisation.
  Once we move to 2.8.xxx change this
*/

int 
const_IWSubstring::operator <(int rhs) const
{
  int intme;
  if (! numeric_value(intme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return intme < rhs;
}

int
const_IWSubstring::operator ==(int rhs) const
{
  int intme;
  if (! numeric_value(intme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return intme == rhs;
}

int
const_IWSubstring::operator <=(int rhs) const
{
  int intme;
  if (! numeric_value(intme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return intme <= rhs;
}

int
const_IWSubstring::operator >=(int rhs) const
{
  int intme;
  if (! numeric_value(intme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return intme >= rhs;
}

int
const_IWSubstring::operator !=(int rhs) const
{
  int intme;
  if (! numeric_value(intme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return intme != rhs;
}

int 
const_IWSubstring::operator <(float rhs) const
{
  float floatme;
  if (! numeric_value(floatme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return floatme < rhs;
}

int
const_IWSubstring::operator ==(float rhs) const
{
  float floatme;
  if (! numeric_value(floatme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return floatme == rhs;
}

int
const_IWSubstring::operator <=(float rhs) const
{
  float floatme;
  if (! numeric_value(floatme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return floatme <= rhs;
}

int
const_IWSubstring::operator >=(float rhs) const
{
  float floatme;
  if (! numeric_value(floatme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return floatme >= rhs;
}

int
const_IWSubstring::operator !=(float rhs) const
{
  float floatme;
  if (! numeric_value(floatme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return floatme != rhs;
}

int 
const_IWSubstring::operator <(double rhs) const
{
  double doubleme;
  if (! numeric_value(doubleme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return doubleme < rhs;
}

int
const_IWSubstring::operator ==(double rhs) const
{
  double doubleme;
  if (! numeric_value(doubleme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return doubleme == rhs;
}

int
const_IWSubstring::operator <=(double rhs) const
{
  double doubleme;
  if (! numeric_value(doubleme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return doubleme <= rhs;
}

int
const_IWSubstring::operator >=(double rhs) const
{
  double doubleme;
  if (! numeric_value(doubleme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return doubleme >= rhs;
}

int
const_IWSubstring::operator !=(double rhs) const
{
  double doubleme;
  if (! numeric_value(doubleme))
  {
    cerr << "Non numeric string value '";
    cerr.write(_data, _nchars) << "'\n";
    abort();
  }

  return doubleme != rhs;
}

static int
common_balance (char open, char close, const char * s, int nchars)
{
  int rc = 0;
  for (int i = 0; i < nchars; i++)
  {
    if (open == s[i])
      rc++;
    else if (close == s[i])
      rc--;
    else if ('\\' == s[i])
      i++;
  }

  return rc;
}

int
const_IWSubstring::balance(char open, char close) const
{
  return common_balance(open, close, _data, _nchars);
}

int
IWString::balance(char open, char close) const
{
  return common_balance(open, close, _things, _number_elements);
}

/*
  Append N copies of S to the end
*/

void
IWString::append(int n, const char * s)
{
  _append(n, s, static_cast<int>(::strlen(s)));

  return;
}

void
IWString::append(int n, const IWString & s)
{
  _append(n, s.rawchars(), s.length());

  return;
}

void
IWString::_append(int n, const char * s, int nchars)
{
  int storage_needed = _number_elements + n * nchars;

  if (_elements_allocated < storage_needed)
    resize(storage_needed);

  for (int i = 0; i < n; i++)
  {
    IW_STRNCPY(_things + _number_elements, s, nchars);
    _number_elements += nchars;
  }

  return;
}

const_IWSubstring
const_IWSubstring::before(char c) const
{
  int i = index(c);
  if (i <= 0)
    return "";

  return const_IWSubstring(_data, i);
}

const_IWSubstring
const_IWSubstring::after(char c) const
{
  int i = index(c);
  if (i <= 0)
    return "";

  return const_IWSubstring(_data + i + 1, _nchars - i - 1);
}

const_IWSubstring
IWString::before(char c) const
{
  int i = index(c);
  if (i <= 0)
    return "";

  return const_IWSubstring(_things, i);
}

const_IWSubstring
IWString::after(char c) const
{
  int i = index(c);
  if (i <= 0)
    return "";

  return const_IWSubstring(_things + i + 1, _number_elements - i - 1);
}

template <typename T>
int
common_basename (const char * s, int nchars,
                 T & zresult,
                 char separator)

{
  if (0 == nchars)
  {
    zresult = "";
    return 1;
  }

  while (separator == s[nchars - 1])     // strip trailing 'separator' characters
  {
    nchars--;
    if (0 == nchars)
    {
      zresult = "";
      return 1;
    }
  }

  int last_separator = nchars + 1;

  for (int i = nchars - 1; i >= 0; i--)
  {
    if (separator == s[i])
    {
      last_separator = i;
      break;
    }
  }

  if (last_separator > nchars)    // no separator found
    zresult.set(s, nchars);
  else
    zresult.set(s + last_separator + 1, nchars - last_separator - 1);

  return 1;
}

//#if defined (__GNUG__) || defined(__SUNPRO_CC)
template int common_basename(const char *, int, const_IWSubstring &, char);
template int common_basename(const char *, int, IWString &, char);
//#endif

void
const_IWSubstring::iwbasename(const_IWSubstring & zresult, char separator) const
{
  common_basename(_data, _nchars, zresult, separator);

  return;
}

void
IWString::iwbasename(const_IWSubstring & zresult, char separator) const
{
  common_basename(_things, _number_elements, zresult, separator);

  return;
}

void
const_IWSubstring::iwbasename(IWString & zresult, char separator) const
{
  common_basename(_data, _nchars, zresult, separator);

  return;
}

void
IWString::iwbasename(IWString & zresult, char separator) const
{
  common_basename(_things, _number_elements, zresult, separator);

  return;
}

IWString &
IWString::append_to_space_separated_list(const IWString & rhs)
{
  if (_number_elements)
  {
    resizable_array<char>::make_room_for_extra_items(1 + rhs._number_elements);
    resizable_array<char>::add(' ');
  }
  else
    resizable_array<char>::make_room_for_extra_items(rhs._number_elements);

  resizable_array<char>::add(rhs._things, rhs._number_elements);

  return *this;
}

int
IWString::remove_chars(int istart, int nchars)
{
  assert (istart >= 0 && istart < _number_elements);
  assert (nchars >= 0);
  assert (istart + nchars <= _number_elements);

  int istop = _number_elements - nchars;
  assert (istop >= istart);

  for (int i = istart; i < istop; i++)
  {
    _things[i] = _things[i + nchars];
  }

  _number_elements -= nchars;

  return nchars;
}

int
IWString::remove_from_to(int zfrom, int zto)
{
  return remove_chars(zfrom, zto - zfrom + 1);
}

int 
IWString::_common_gsub(const char * zfrom,
                        int len_from,
                        const char * zto,
                        int len_to)
{
  if (0 == _number_elements)
    return 0;

  if (len_from > _number_elements)
    return 0;

  if (len_from > len_to)
    return _common_gsub_getting_smaller(zfrom, len_from, zto, len_to);
  else if (len_from < len_to)
    return _common_gsub_getting_larger(zfrom, len_from, zto, len_to);
  else
    return _common_gsub_same_size(zfrom, len_from, zto);
}

/*
  Expanding the string is hard. We need some heuristics to make sure we aren't too
  inefficient, without being wasteful
*/

int
IWString::_common_gsub_getting_larger(const char * zfrom,
                        int len_from,
                        const char * zto,
                        int len_to)
{
  int max_possible_matches = _number_elements / len_from;

  int delta = len_to - len_from;

  assert (delta > 0);

  int extra_allocation = 0;

  if (0 == max_possible_matches)   // only one resize possible
    extra_allocation = 0;
  else if (max_possible_matches < 10)
    extra_allocation = 5 * delta;
  else
    extra_allocation = 15 * delta;

  int rc = 0;

  int iptr = 0;

/*int increment_iptr;
  if (len_from > len_to)
    increment_iptr = len_to;
  else 
    increment_iptr = len_from;*/

  while (iptr < _number_elements - len_from + 1)
  {
    if (_things[iptr] != zfrom[0] || 0 != ::strncmp(_things + iptr, zfrom, len_from))
    {
      iptr++;
      continue;
    }

    if (_number_elements + delta > _elements_allocated)
    {
      if (rc > 1 && 0 == rc % 10)           // looks like lots of matches
        extra_allocation = extra_allocation * 2;
      
      make_room_for_extra_items(extra_allocation);
    }

    ::memmove(_things + iptr + len_to, _things + iptr + len_from, _number_elements - iptr - len_from);
    _number_elements += delta;

    ::memcpy(_things + iptr, zto, len_to);

    iptr = iptr + len_to;

    rc++;
  }

  return rc;
}

int
IWString::_common_gsub_same_size(const char * zfrom,
                                  int len_from,
                                  const char * zto)
{
  int rc = 0;

  int iptr = 0;

  while (iptr < _number_elements - len_from + 1)
  {
    if (_things[iptr] != zfrom[0] || 0 != ::strncmp(_things + iptr, zfrom, len_from))
    {
      iptr++;
      continue;
    }

    ::memcpy(_things + iptr, zto, len_from);    // copy replacement string in place

    iptr += len_from;

    rc++;
  }

  return rc;
}

int
IWString::_common_gsub_getting_smaller(const char * zfrom,
                                        int len_from,
                                        const char * zto,
                                        int len_to)
{
  assert (len_from > len_to);

  int rc = 0;

  int readfrom = 0;
  int putback = 0;

  while (readfrom < _number_elements - len_from + 1)
  {
    if (_things[readfrom] != zfrom[0] || 0 != ::strncmp(_things + readfrom, zfrom, len_from))
    {
      if (readfrom != putback)
        _things[putback] = _things[readfrom];

      readfrom++;
      putback++;

      continue;
    }

    ::memcpy(_things + putback, zto, len_to);    // copy replacement string in place

    readfrom += len_from;
    putback += len_to;

    rc++;
  }

// Because we don't run the loop above all the way out to _number_elements, we will have
// items still to be copied

  if (rc > 0)  
  {
    if (readfrom < _number_elements)
      ::memmove(_things + putback, _things + readfrom, _number_elements - readfrom);

    _number_elements = _number_elements - rc * (len_from - len_to);
  }

  return rc;
}

int
IWString::gsub(const char * zfrom, const char * zto)
{
  return _common_gsub(zfrom, static_cast<int>(::strlen(zfrom)), zto, static_cast<int>(::strlen(zto)));
}

int
IWString::gsub(const const_IWSubstring & zfrom, const const_IWSubstring & zto)
{
  return _common_gsub(zfrom.rawchars(), zfrom.length(), zto.rawchars(), zto.length());
}

int
IWString::gsub(char zfrom, const IWString & zto)
{
  return _common_gsub(&zfrom, 1, zto.rawchars(), zto.length());
}

int
IWString::gsub(char zfrom, const char * zto)
{
  return _common_gsub(&zfrom, 1, zto, static_cast<int>(::strlen(zto)));
}

static int
common_matches_at_position (const char * s1,
                            int lens1,
                            int offset,
                            const char * s2,
                            int lens2)
{
  if (offset + lens2 > lens1)
    return 0;

  return 0 == ::strncmp(s1 + offset, s2, lens2);
}

int
const_IWSubstring::matches_at_position(int o,
                                        const char * s,
                                        int lens) const
{
  return common_matches_at_position(_data, _nchars, o, s, static_cast<int>(lens));
}

int
const_IWSubstring::matches_at_position(int o,
                                        const IWString & s) const
{
  return common_matches_at_position(_data, _nchars, o, s.rawchars(), s.length());
}

int
IWString::matches_at_position(int o,
                               const char * s,
                               int lens) const
{
  return common_matches_at_position(_things, _number_elements, o, s, static_cast<int>(lens));
}

template <typename T>
int
common_nextword_single_delimiter(const char * s,
                                 int nchars,
                                 int & ndx,
                                 char separator,
                                 T & zresult)
{
  if (ndx >= nchars)
    return 0;

  if (separator == s[ndx])   // consecutive separators, possibly a common scenario
  {
    zresult.make_empty();
    ndx++;
    return 1;
  }

  const int initial_ndx = ndx;

  ndx++;

  while (ndx < nchars)
  {
    if (separator == s[ndx])
    {
      zresult.set(s + initial_ndx, ndx - initial_ndx);
      ndx++;
      return 1;
    }
    else
      ndx++;
  }

// If we get to here, we came off the end

  zresult.set(s + initial_ndx, ndx - initial_ndx);

  return 1;
}


template <typename T>
int
const_IWSubstring::nextword_single_delimiter(T & zresult, int & i,  char separator) const
{
  return common_nextword_single_delimiter(_data, _nchars, i, separator, zresult);
}

template <typename T>
int
IWString::nextword_single_delimiter(T & zresult, int & i,  char separator) const
{
  return common_nextword_single_delimiter(_things, _number_elements, i, separator, zresult);
}


template int IWString::nextword_single_delimiter<const_IWSubstring>(const_IWSubstring&, int&, char) const;
template int IWString::nextword_single_delimiter<IWString>(IWString&, int&, char) const;
template int const_IWSubstring::nextword_single_delimiter<const_IWSubstring>(const_IWSubstring&, int&, char) const;
template int const_IWSubstring::nextword_single_delimiter<IWString>(IWString&, int&, char) const;

template int common_nextword_single_delimiter<const_IWSubstring>(char const*, int, int&, char, const_IWSubstring&);
template int common_nextword_single_delimiter<IWString>(char const*, int, int&, char, IWString&);

#if defined (IW_STD_STRING_DEFINED)
const_IWSubstring::const_IWSubstring(const std::string & rhs)
{
  set(rhs.data(), rhs.length());

  return;
}

const_IWSubstring &
const_IWSubstring::operator =(const std::string & rhs)
{
  set(rhs.data(), rhs.length());

  return *this;
}

IWString::IWString(const std::string & rhs)
{
  _default_values();

  IWString::strncpy(rhs.data(), static_cast<int>(rhs.length()));
}

IWString::IWString(int s)
{
  _default_values();

  resize(s);

  return;
}

IWString::IWString(unsigned int s)
{
  _default_values();

  resize(s);

  return;
}

IWString &
IWString::operator =(const std::string & rhs)
{
  IWString::strncpy(rhs.data(), static_cast<int>(rhs.length()));

  return *this;
}

void
IWString::operator +=(const std::string & rhs)
{
  const_IWSubstring tmp(rhs);

  operator += (tmp);

  return;
}

#endif

int
IWString::EnsureEndsWith(const char* s) {
  const_IWSubstring mycopy(s, ::strlen(s));

  return this->EnsureEndsWith<const_IWSubstring&>(mycopy);
}

template <typename T>
int
IWString::EnsureEndsWith(const T& c) {
  if (_number_elements == 0) {
    operator<<(c);
    return 1;
  }

  if (ends_with(c)) {
    return 0;
  }

  operator<<(c);

  return 1;
}

template int IWString::EnsureEndsWith(const char&);
template int IWString::EnsureEndsWith(const IWString&);
template int IWString::EnsureEndsWith(const const_IWSubstring&);

/*
  Some names taken from 

  https://blog.codinghorror.com/ascii-pronunciation-rules-for-programmers/
*/

void char_name_to_char_usage(const IWString & optionTag)
{
  cerr << "    " << optionTag << " tab" << endl;

  cerr << "    " << optionTag << " comma -  for a comma (\",\") separator" << endl;
  cerr << "    " << optionTag << " space -  for a space (\" \") separator" << endl;
  cerr << "    " << optionTag << " vbar or bar -  for a vertical bar (\"|\") separator" << endl;
  cerr << "    " << optionTag << " squote -  for a single quote  (\"'\") separator" << endl;
  cerr << "    " << optionTag << " dquote -  for a comma ('\"') separator" << endl;
  cerr << "    " << optionTag << " semic -  for a semicolon (\";\") separator" << endl;
  cerr << "    " << optionTag << " colon -  for a comma (\",\") separator" << endl;
  cerr << "    " << optionTag << " and or amp -  for a apersand (\"&\") separator" << endl;
  cerr << "    " << optionTag << " dollar -  for a dollar sign (\"$\") separator" << endl;
  cerr << "    " << optionTag << " hash -  for a hash tag or number sign (\"#\") separator" << endl;
  cerr << "    " << optionTag << " at -  for an at sign (\"@\") separator" << endl;
  cerr << "    " << optionTag << " excl or bang -  for a exclmation point (\"!\") separator" << endl;
  cerr << "    " << optionTag << " squiggle or tilde -  for a tilda (\"~\") separator" << endl;
  cerr << "    " << optionTag << " star or asterisk -  for an asterisk (\"*\") separator" << endl;
  cerr << "    " << optionTag << " hat or caret -  for a caret (\"^\") separator" << endl;
  cerr << "    " << optionTag << " oparen or lparen -  for an open parenthesis (\"(\") separator" << endl;
  cerr << "    " << optionTag << " cparen or rparen -  for a close parenthesis (\")\") separator" << endl;
  cerr << "    " << optionTag << " uscore -  for an underscore (\"_\") separator" << endl;
  cerr << "    " << optionTag << " plus -  for a plus sign (\"+\") separator" << endl;
  cerr << "    " << optionTag << " minus or dash -  for a dash (\"-\") separator" << endl;
  cerr << "    " << optionTag << " ocbrace or lcbrace -  for an open curly brace (\"{\") separator" << endl;
  cerr << "    " << optionTag << " ccbrace or rcbrace -  for a close curly brace (\"}\") separator" << endl;
  cerr << "    " << optionTag << " osqb or lsqb -  for an opeen square bracket (\"[\") separator" << endl;
  cerr << "    " << optionTag << " csqb or rsqb -  for a close square bracket (\"]\") separator" << endl;
  cerr << "    " << optionTag << " bslash -  for a backslash (\"\\\") separator" << endl;
  cerr << "    " << optionTag << " fslash or slash -  for a slash (\"/\") separator" << endl;
  cerr << "    " << optionTag << " dot or period -  for a period (\".\") separator" << endl;
  cerr << "    " << optionTag << " oangle or less -  for a less than sign (\"<\") separator" << endl;
  cerr << "    " << optionTag << " cangle or greater -  for a greater than sign  (\">\") separator" << endl;
  cerr << "    " << optionTag << " qmark -  for a question mark (\"?\") separator" << endl;
  cerr << "    " << optionTag << " pct -  for a percent sign (\"%\") separator" << endl;
  cerr << "    " << optionTag << " grave or backquote -  for a back quote sign (\"`\") separator" << endl;
}


int
char_name_to_char(IWString & s, int message_if_unrecognised) {
  if (1 == s.length())
    return 1;

  s.to_lowercase();
  if ("tab" == s)
  {
    s = '\t';
    return 1;
  }

  if ("comma" == s)
  {
    s = ',';
    return 1;
  }

  if ("space" == s)
  {
    s = ' ';
    return 1;
  }

  if ("vbar" == s || "bar" == s)
  {
    s = '|';
    return 1;
  }

  if ("squote" == s)
  {
    s = '\'';
    return 1;
  }

  if ("dquote" == s)
  {
    s = '"';
    return 1;
  }

  if ("semic" == s)
  {
    s = ';';
    return 1;
  }

  if ("colon" == s)
  {
    s = ':';
    return 1;
  }

  if ("and" == s || "amp" == s)
  {
    s = '&';
    return 1;
  }

  if ("dollar" == s)
  {
    s = '$';
    return 1;
  }

  if ("hash" == s)
  {
    s = '#';
    return 1;
  }

  if ("at" == s)
  {
    s = '@';
    return 1;
  }

  if ("excl" == s || "bang" == s)
  {
    s = '!';
    return 1;
  }

  if ("squiggle" == s || "tilde" == s)
  {
    s = '~';
    return 1;
  }

  if ("star" == s || "asterisk" == s)
  {
    s = '*';
    return 1;
  }

  if ("hat" == s || "caret" == s)
  {
    s = '^';
    return 1;
  }

  if ("oparen" == s || "lparen" == s)
  {
    s = '(';
    return 1;
  }

  if ("cparen" == s || "rparen" == s)
  {
    s = ')';
    return 1;
  }

  if ("uscore" == s)
  {
    s = '_';
    return 1;
  }

  if ("plus" == s)
  {
    s = '+';
    return 1;
  }

  if ("minus" == s || "dash" == s)
  {
    s = '-';
    return 1;
  }

  if ("ocbrace" == s || "lcbrace" == s)
  {
    s = '{';
    return 1;
  }

  if ("ccbrace" == s || "rcbrace" == s)
  {
    s = '}';
    return 1;
  }

  if ("osqb" == s || "lsqb" == s)
  {
    s = '[';
    return 1;
  }

  if ("csqb" == s || "rsqb" == s)
  {
    s = ']';
    return 1;
  }

  if ("bslash" == s)
  {
    s = '\\';
    return 1;
  }

  if ("fslash" == s || "slash" == s)
  {
    s = '/';
    return 1;
  }

  if ("dot" == s || "period" == s)
  {
    s = '.';
    return 1;
  }

  if ("oangle" == s || "less" == s)
  {
    s = '<';
    return 1;
  }

  if ("cangle" == s || "greater" == s)
  {
    s = '>';
    return 1;
  }

  if ("qmark" == s)
  {
    s = '?';
    return 1;
  }

  if ("pct" == s)
  {
    s = '%';
    return 1;
  }

  if ("grave" == s || "backquote" == s)
  {
    s = '`';
    return 1;
  }

  else if (s == "empty") {
    s.resize(0);
    return 1;
  }
  
  if ("help" == s)
  {
    char_name_to_char_usage( "");
    return (0);
  }

  if (message_if_unrecognised) {
    cerr << "Unrecognised character name '" << s << "'\n";
  }
  return 0;
}


int
const_IWSubstring::expand_environment_variables(IWString & destination) const
{
  int needs_close_brace = 0;

  int rc = 0;

  for (int i = 0; i < _nchars; ++i)
  {
    if ('$' != _data[i])
    {
      destination += _data[i];
      continue;
    }

    if (i == _nchars - 1)    // ends with $. Should this be a failure instead?
    {
      destination += '$';
      return rc;
    }

    ++i;
    if ('{' == _data[i])
    {
      ++i;
      needs_close_brace = 1;
    }
    else
      needs_close_brace = 0;

    IWString v;
    for ( ; i < _nchars; ++i)
    {
      if (isalnum(_data[i]) || '_' == _data[i])
        v += _data[i];
      else
        break;
    }

    if (needs_close_brace)
    {
      if (_data[i] != '}')
      {
        cerr << "const_IWSubstring::expand_environment_variables:unclosed open brace, got '" << v << "'\n";
        return 0;
      }
    }

    const char * env = getenv(v.null_terminated_chars());
    if (nullptr == env)
    {
      cerr << "const_IWSubstring::expand_environment_variables:no value for '" << env << "'\n";
      return 0;
    }
//  cerr << "Expanded '" << v << "' to '" << env << "'\n";

    destination += env;
    rc++;
  }

  return rc;
}

static constexpr char kOpenBrace = '{';
static constexpr char kCloseBrace = '}';

std::optional<IWString>
IWString::ExpandEnvironmentVariables() const {
  IWString result;
  // Must contain at least ${x}.
  if (_number_elements < 4) {
    result = *this;
    return result;
  }

  result.resize(_number_elements + 50);  // just a guess, might even shrink...

  for (int i = 0; i < _number_elements; ++i) {
    if ('$' != _things[i]) {
      result += _things[i];
      continue;
    }

    if (i == _number_elements - 1) {  // ends with $.
      result += '$';
      return result;
    }

    if (_things[i + 1] != kOpenBrace) {  // $ without following braces, not expanded.
      result += '$';
      continue;
    }

    i += 2;  // skip over $ and opening brace

    int start_of_variable = i;
    int closing_brace = -1;
    for ( ; i < _number_elements; ++i) {
      if (_things[i] == kCloseBrace) {
        closing_brace = i;
        break;
      }
    }
    if (closing_brace < 0) {
      cerr << "IWString::ExpandEnvironmentVariables:no closing brace '" << *this << "'\n";
      return std::nullopt;
    }

    if (start_of_variable == closing_brace) {
      cerr << "IWString::ExpandEnvironmentVariables:empty variable spec '" << *this << "'\n";
      return std::nullopt;
    }

    // cerr << "varname from " << start_of_variable << " to " << closing_brace << '\n';
    IWString vname(_things + start_of_variable, closing_brace - start_of_variable);
    // cerr << "Variable is '" << vname << "'\n";

    const char * env = getenv(vname.null_terminated_chars());
    if (nullptr == env)
    {
      cerr << "IWString::ExpandEnvironmentVariables:no value for '" << vname << "'\n";
      return std::nullopt;
    }
    // cerr << "Expanded '" << vname << "' to '" << env << "'\n";

    result += env;
    i = closing_brace;
  }

  return result;
}

std::string
IWString::AsString() const {
  if (_number_elements == 0) {
    return std::string();
  }

  return std::string(_things, _number_elements);
}

std::string
const_IWSubstring::AsString() const {
  if (_nchars == 0) {
    return std::string();
  }

  return std::string(_data, _nchars);
}
