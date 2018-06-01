#include "stdlib.h"
#include "strings.h"
#include "string.h"
#include "ctype.h"
#include "math.h"

#include "iostream.h"
#include "strstream.h"
#include "fstream.h"
#include "iwstring.h"

/*
  The philosophy of this class is that everything is inherited from
  resizable_array<char>.


*/

void
IWString::_default_values ()
{
}

IWString::IWString ()
{
  _default_values ();

  return;
}

IWString::IWString (const char * s)
{
  _default_values ();

  *this = s;
  return;
}

IWString::IWString (const IWString & s)
{
  _default_values ();

  resize (s._elements_allocated);
  memcpy (_things, s._things, s._number_elements);

//for (int i = 0; i < s._number_elements; i++)
//  _things[i] = s._things[i];

  _number_elements = s._number_elements;

//*this = s;

  return;
}

IWString::IWString (const const_IWSubstring & s)
{
  _default_values ();

  *this = s;

  return;
}

IWString::IWString (char c)
{
  _default_values ();

  *this = c;

  return;
}

int
IWString::ok () const
{
  if (! resizable_array<char>::ok ())
    return 0;

  return 1;
}

void
IWString::strip_leading_blanks ()
{
  int i;
  for (i = 0; i < _number_elements; i++)
    if (! isspace (_things[i]))
      break;

  if (0 == i)     // no leading spaces
    return;

  _number_elements -= i;

  int j;     // BARF
  for (j = 0; j < _number_elements; j++)
    _things[j] = _things[i + j];

  _things[_number_elements] = '\0';

  return;
}

void
const_IWSubstring::strip_leading_blanks ()
{
  while (_nchars && isspace (_data[0]))
  {
    _data++;
    _nchars--;
  }

  if (0 == _nchars)
    _data = NULL;

  return;
}

void
IWString::strip_trailing_blanks ()
{
  int i;     // BARF
  for (i = _number_elements - 1; i >= 0; i--)
  {
    if (! isspace (_things[i]))
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
const_IWSubstring::strip_trailing_blanks ()
{
  while (_nchars && isspace (_data[_nchars - 1]))
  {
    _nchars--;
  }

  if (0 == _nchars)
    _data = NULL;

  return;
}

/*
  Note that we don't increase the number of elements. This is really
  just so that callers can safely use rawchars ().
*/

int
IWString::null_terminate ()
{
  if (_number_elements == _elements_allocated)
    resize (_elements_allocated + 1);

  if ('\0' == _things[_number_elements])
    return 0;

  _things[_number_elements] = '\0';

  return 1;
}

/*
  This is pretty awful, but seemed like the easiest way of getting ostream << to work.
*/

void
IWString::_const_appearing_null_terminate () const
{
  ((IWString *)this)->null_terminate ();

  return;
}

/*
  Some operation has occurred which shortened the string representation.
  Find the '\0' and adjust _number_elements;
*/

void
IWString::_recompute_length ()
{
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
  {
    if ('\0' == _things[i])
    {
      _number_elements = i;
      return;
    }
  }

// If we pass to here, we did not find a null terminator. For now, let's assume that's OK

  return;
}

const char *
IWString::chars ()
{
  assert (ok ());

  null_terminate ();

  return _things;
}

const char *
IWString::null_terminated_chars ()
{
  assert (ok ());

  null_terminate ();

  return _things;
}
void
IWString::chop (int nchars)
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
IWString::truncate (int nchars)
{
  assert (nchars >= 0 && nchars <= _number_elements);
  if (nchars == _number_elements)    // already the appropriate length
    return;

  _number_elements = nchars;
  if (_number_elements < _elements_allocated)
    _things[_number_elements] = '\0';

  return;
}

int
IWString::truncate_at_first (char c)
{
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
  {
    if (c == _things[i])
    {
      _number_elements = i;
      return 1;
    }
  }

  return 0;
}
void
IWString::remove_leading_chars (int nremove)
{
  assert (nremove >= 0 && nremove <= _number_elements);

  int i;     // BARF
  for (i = nremove; i < _number_elements; i++)
    _things[i - nremove] = _things[i];

  _number_elements = _number_elements - nremove;

  return;
}

int
IWString::remove_from_to (int cfrom, int cto)
{
  assert (cfrom >= 0 && cfrom <= cto && cto < _number_elements);

  int j = cfrom;
  int i;     // BARF
  for (i = cto + 1; i < _number_elements; i++)
    _things[j++] = _things[i];

  _number_elements = j;

  return 1;
}

int
IWString::remove_all_these (const char * remove_these)
{
  int rc = 0;
  int j = 0;
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
  {
    char c = _things[i];
    if (NULL == ::strchr (remove_these, c))
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
IWString::compress_blanks ()
{
  null_terminate ();

  int rc = ::compress_blanks (_things);

  _recompute_length ();

  return rc;
}

IWString &
IWString::operator= (const char * rhs)
{
  int nchars = strlen (rhs);

  if (0 == nchars)
  {
    resize (0);
    return *this;
  }

  resize (nchars + 1);
  assert (_things);

  strcpy (_things, rhs);
  _number_elements = nchars;

  assert (ok ());

  return *this;
}

const_IWSubstring &
const_IWSubstring::operator= (const char * rhs)
{
  _data = rhs;

  if (NULL == rhs)
    _nchars = 0;
  else
    _nchars = strlen (rhs);

  return *this;
}

const_IWSubstring &
const_IWSubstring::operator= (const IWString & rhs)
{
  _nchars = rhs.length ();

  _data = rhs.rawchars ();

  return *this;
}

void
IWString::operator += (const char * rhs)
{
  int nchars = strlen (rhs);
  if (0 == nchars)
    return;

  (void) add (rhs, nchars);

  return;
}

void
IWString::operator += (const IWString & rhs)
{
  if (0 == rhs._number_elements)
    return;

  add (rhs._things, rhs._number_elements);
  
  return;
}

void
IWString::operator += (const const_IWSubstring & rhs)
{
  if (0 == rhs._nchars)
    return;

  (void) add (rhs._data, rhs._nchars);

  return;
}

void
IWString::to_lowercase ()
{
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
    _things[i] = tolower (_things[i]);

  return;
}

void
IWString::to_uppercase ()
{
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
    _things[i] = toupper (_things[i]);

  return;
}

/*
  We need to turn strcmp into a strncmp because neither of the strings may
  be null terminated.
*/

static int
_strcmp (const char * s1, int l1, const char * s2, int l2)
{
  int ncomp = l1;
  if (l2 < l1)
    ncomp = l2;

  return strncmp (s1, s2, ncomp);
}

int
IWString::strcmp (const IWString & rhs) const
{
  int l2 = rhs._number_elements;

  return _strcmp (_things, _number_elements, rhs._things, l2);
}

int
IWString::strcmp (const char * rhs) const
{
  int l2 = strlen (rhs);

  return _strcmp (_things, _number_elements, rhs, l2);
}

/*
  Because neither of the string may be null terminated, we need to be
  careful in doing a strncmp
*/

static int
_strncmp (const char * s1, int l1, const char * s2, int l2,
          int chars_to_compare)
{
  int ncomp = chars_to_compare;
  if (l1 < ncomp)
    ncomp = l1;
  if (l2 < ncomp)
    ncomp = l2;

  return strncmp (s1, s2, ncomp);
}

int
IWString::strncmp (const IWString & rhs, int chars_to_compare) const
{
  return _strncmp (_things, _number_elements, rhs._things, rhs._number_elements, chars_to_compare);
}

int
IWString::strncmp (const char * rhs, int chars_to_compare) const
{
  return _strncmp (_things, _number_elements, rhs, strlen (rhs), chars_to_compare);
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

  int i;     // BARF
  for (i = 0; i <= last_i; i++)
  {
    if (haystack[i] != needle[0])
      continue;

    int found_match = 1;
    for (int j = 1; j < size_of_needle; j++)
    {
      if (haystack[i + j] != needle[j])
      {
        found_match = 0;
        break;
      }
    }

    if (found_match)
      return i;
  }

  return -1;
}

int
IWString::find (const char * needle) const
{
  return do_find (_things, _number_elements, needle, strlen (needle));
}

int
const_IWSubstring::find (const char * needle) const
{
  return do_find (_data, _nchars, needle, strlen (needle));
}

int
IWString::find (const IWString & needle) const
{
  return do_find (_things, _number_elements, needle._things, needle._number_elements);
}

int
const_IWSubstring::find (const IWString & needle) const
{
  return do_find (_data, _nchars, needle.rawchars (), needle.nchars ());
}

int
IWString::find (const const_IWSubstring & needle) const
{
  return do_find (_things, _number_elements, needle._data, needle._nchars);
}

int
const_IWSubstring::find (const const_IWSubstring & needle) const
{
  return do_find (_data, _nchars, needle._data, needle._nchars);
}

/*
  Pretty dumb int conversion. Note that this does not deal adequately with
  large negative values and such.
*/

static int
string_class_is_int (const char * s, int nchars, int & result)
{
  if (nchars > 9)
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

  result = 0;
  int i;     // BARF
  for (i = 0; i < nchars; i++)
  {
    if (s[i] >= '0' && s[i] <= '9')
    {
      result = 10 * result + (s[i] - '0');
    }
    else
      return 0;
  }

  if (sign)
    result = - result;

  return 1;
}

int
IWString::is_int (int & result) const
{
  if (0 == _number_elements)
    return 0;

  return string_class_is_int (_things, _number_elements, result);
}

int
IWString::numeric_value (int & result) const
{
  if (0 == _number_elements)
    return 0;

  return string_class_is_int (_things, _number_elements, result);
}

/*
  Somewhat of a kludge for the hex value
*/

int
IWString::is_hex (unsigned int & result) const
{
  if (_number_elements < 3)
    return 0;

  _const_appearing_null_terminate ();


  char *c;
  long tmp = strtol (_things, &c, 16);

  if (c == _things)
    return 0;

  if ('\0' == *c)
    ;     // good
  else if (! isspace (*c))
    return 0;

  result = tmp;

  return 1;
}

int
IWString::is_double (double & result) const
{
  _const_appearing_null_terminate ();

  return ::is_double (_things, &result);
}

int
IWString::numeric_value (double & result) const
{
  _const_appearing_null_terminate ();

  return ::is_double (_things, &result);
}

static int
common_nwords (const char * s,
               int nchars)
{
  int rc = 0;
  int in_word = 0;
  int i;     // BARF
  for (i = 0; i < nchars; i++)
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
IWString::nwords () const
{
  return common_nwords (_things, _number_elements);
}

int
const_IWSubstring::nwords () const
{
  return common_nwords (_data, _nchars);
}

static int
common_nwords_with_separator (const char * s,
               int nchars,
               char word_separator)
{
  int rc = 0;
  int in_word = 0;
  int i;     // BARF
  for (i = 0; i < nchars; i++)
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
IWString::nwords (char word_separator) const
{
  return common_nwords_with_separator (_things, _number_elements, word_separator);
}

int
const_IWSubstring::nwords (char word_separator) const
{
  return common_nwords_with_separator (_data, _nchars, word_separator);
}

static int
common_locate_word_boundaries (const char * s, int nchars,
                               int which_word,
                               char word_separator,
                               int & word_start, int & word_stop)
{
  int in_word = 0;
  int word_count = 0;
  int i;     // BARF
  for (i = 0; i < nchars; i++)
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
IWString::_locate_word_boundaries (int which_word, char word_separator,
                                   int & word_start, int & word_stop) const
{
  int in_word = 0;
  int word_count = 0;
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
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
IWString::word (int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
    return const_IWSubstring ("", 0);

  return const_IWSubstring (&(_things[word_start]), word_stop - word_start + 1);
}

int
IWString::word (int which_word, IWString & result, char word_separator) const
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
IWString::word (int which_word, const_IWSubstring & result, char word_separator) const
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
IWString::wordindex (int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_things, _number_elements,
                 which_word, word_separator, word_start, word_stop))
    return -1;

  return word_start;
}

const_IWSubstring
const_IWSubstring::word (int which_word, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_data, _nchars,
                 which_word, word_separator, word_start, word_stop))
    return const_IWSubstring ("", 0);

  return const_IWSubstring (&(_data[word_start]), word_stop - word_start + 1);
}

int
const_IWSubstring::word (int which_word, IWString & result, char word_separator) const
{
  int word_start, word_stop;

  if (! common_locate_word_boundaries (_data, _nchars,
                 which_word, word_separator, word_start, word_stop))
  {
    result = "";
    return 0;
  }

  result = from_to (word_start, word_stop);

  return 1;
}

int
const_IWSubstring::word (int which_word, const_IWSubstring & result, char word_separator) const
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
const_IWSubstring::wordindex (int which_word, char word_separator) const
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

  int i;     // BARF
  for (i = 0; i < nchars; i++)
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
const_IWSubstring::locate_word_beginnings (resizable_array<int> & word_beginnings,
                                           char separator) const
{
  return _locate_word_beginnings (_data, _nchars, word_beginnings, separator);
}

int
IWString::locate_word_beginnings (resizable_array<int> & word_beginnings,
                                           char separator) const
{
  return _locate_word_beginnings (_things, _number_elements, word_beginnings, separator);
}

/*
  The code for nextword is common between IWString and constIWSubstring.
*/

static int
internal_nextword (const char * s, int nchars, int & i,
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

  istart = i;
  istop = i;

  while (i < nchars)
  {
    i++;
    if (separator == s[i])
      return 1;

    istop++;
  }

// If we ran off the end of the string, decrement istop

  istop--;

  return 1;
}
  
static int
internal_nextword (IWString & result, int & i,
                   char separator,
                   const char * s, int nchars)
{
  result.resize_keep_storage (0);

  int istart, istop;

  if (! internal_nextword (s, nchars, i, separator, istart, istop))
    return 0;

//cerr << "After STRING determination, istart = " << istart << " istop = " << istop << endl;

  result.resize (istop - istart + 1);

  int j;     // BARF
  for (j = istart; j <= istop; j++)
    result += s[j];

  return 1;
}

static int
internal_nextword (const_IWSubstring & result, int & i,
                   char separator,
                   const char * s, int nchars)
{
  result.const_IWSubstring (NULL, 0);

  int istart, istop;

  if (! internal_nextword (s, nchars, i, separator, istart, istop))
    return 0;

//cerr << "After determination, istart = " << istart << " istop = " << istop << endl;

  result.const_IWSubstring (&s[istart], istop - istart + 1);

  return 1;
}

/*
  Used when scanning for words backwards in a string.
  Note that we must scan back to the start of the word, and then
  add all the letters up to where we started.
*/

static int
internal_prevword (IWString & result, int & i,
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

  while (i >= 0)
  {
    i--;
    if (separator == s[i])
    {
      start_of_word = i + 1;
      break;
    }
  }

// Now build the word from i to end_of_word

  int j;     // BARF
  for (j = start_of_word; j <= end_of_word; j++)
    result += s[j];

  return 1;
}

int
IWString::nextword (IWString & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword (result, i, separator, _things, _number_elements);
}

int
const_IWSubstring::nextword (IWString & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword (result, i, separator, _data, _nchars);
}

int
IWString::nextword (const_IWSubstring & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword (result, i, separator, _things, _number_elements);
}

int
const_IWSubstring::nextword (const_IWSubstring & result, int & i, char separator) const
{
  assert (i >= 0);

  return internal_nextword (result, i, separator, _data, _nchars);
}

int
IWString::prevword (IWString & result, int & i, char separator) const
{
  if (i <= 0)
    return 0;

  assert (i < _number_elements);

  return internal_prevword (result, i, separator, _things);
}

int
const_IWSubstring::prevword (IWString & result, int & i, char separator) const
{
  if (i <= 0)
    return 0;

  assert (i < _nchars);

  return internal_prevword (result, i, separator, _data);
}

const_IWSubstring
IWString::substr (int istart, int nchars) const
{
  assert (istart >= 0 && istart < _number_elements);
  if (nchars < 0)
    nchars = _number_elements - istart;

  assert (nchars >= 0 && istart + nchars <= _number_elements);

  return const_IWSubstring (&(_things[istart]), nchars);
}

const_IWSubstring
IWString::from_to (int istart, int istop) const
{
  assert (istart >= 0);
  assert (istop >= istart && istop < _number_elements);

  return const_IWSubstring (&(_things[istart]), istop - istart + 1);
}

const_IWSubstring
const_IWSubstring::from_to (int istart, int istop) const
{
  assert (istart >= 0);
  assert (istop >= istart && istop < _nchars);

  return const_IWSubstring (&(_data[istart]), istop - istart + 1);
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
IWString::from_to (int istart, const char * s) const
{
  assert (istart >= 0 && istart < _number_elements);

  if (0 == _number_elements)
    return const_IWSubstring (NULL, 0);

  int lens = strlen (s);

  if (0 == lens)
    return const_IWSubstring (NULL, 0);

  int i = do_find (_things, _number_elements, s, lens);

  if (i < 0)
    return const_IWSubstring (NULL, 0);

  return const_IWSubstring (&_things[istart], i - istart + 1);
}

void
IWString::from_to (int istart, int istop, char * destination) const
{
  assert (istart >= 0);
  assert (istop >= istart && istart < _number_elements);

  int nchars = istop - istart + 1;
  ::strncpy (destination, &(_things[istart]), nchars);
  destination[nchars] = '\0';

  return;
}

void
IWString::from_to (int istart, int istop, const_IWSubstring & destination) const
{
  assert (ok_index (istart));
  assert (istop >= istart && istop < _number_elements);

  destination._data = &_things[istart];
  destination._nchars = istop - istart + 1;

  return;
}

void
const_IWSubstring::from_to (int istart, int istop, const_IWSubstring & destination) const
{
  assert (istart >= 0 && istart < _nchars);
  assert (istop >= istart && istop < _nchars);

  destination._data = &_data[istart];
  destination._nchars = istop - istart + 1;

  return;
}

IWString &
IWString::operator= (const const_IWSubstring & rhs)
{
  if (_elements_allocated < rhs._nchars)
    resize (rhs._nchars + 1);

  ::strncpy (_things, rhs._data, rhs._nchars);
  _number_elements = rhs._nchars;

  return *this;
}

IWString &
IWString::operator= (char c)
{
  if (_elements_allocated < 1)
    resize (2);

  _things[0] = c;
  if (_elements_allocated > 1)
    _things[1] = '\0';

  _number_elements = 1;

  return *this;
}

int
IWString::getline (istream & is, char record_delimiter)
{
  if (! is.good ())
  {
    if (is.eof ())
      return 0;
    cerr << "IWString::getline: input stream not good\n";
    return 0;
  }

  if (is.eof ())
    return 0;

  if (0 == _elements_allocated)
    resize (256);

  _number_elements = 0;

  while (1)
  {
    int bytes_available = _elements_allocated - _number_elements;
    is.get (&(_things[_number_elements]), bytes_available, record_delimiter);
    int nchars = is.gcount ();
//  cerr << "Just read " << nchars << " characters from stream\n";
    _number_elements += nchars;
//  cerr << "n = " << _number_elements << " alloc = " << _elements_allocated << endl;

    if (is.eof ())
    {
      if (_number_elements)
        cerr << "Incomplete last record\n";

      return _number_elements;
    }

    if (_number_elements < _elements_allocated - 1)   // whole record fit into buffer
    {
      int nextchar = is.get ();
      if (is.eof ())
      {
        cerr << "IWString::getline: record not terminated\n";
        return _number_elements;
      }

      return _number_elements;
    }

    resize (_elements_allocated * 2);
//  cerr << "Buffer expanded to " << _elements_allocated << endl;
  }
}

const_IWSubstring::const_IWSubstring ()
{
  _data = NULL;
  _nchars = 0;

  return;
}

const_IWSubstring::const_IWSubstring (const char * data, int nchars) :
          _data (data), _nchars (nchars)
{
  return;
}

const_IWSubstring::const_IWSubstring (const char * data) :
          _data (data)
{
  _nchars = ::strlen (data);

  return;
}


const_IWSubstring::const_IWSubstring (const IWString & s) :
          _data (s._things), _nchars (s._number_elements)
{
  return;
}

const_IWSubstring::~const_IWSubstring ()
{
  _nchars = 0;
}

const_IWSubstring &
const_IWSubstring::operator = (const const_IWSubstring & other)
{
  _data = other._data;
  _nchars = other._nchars;

  return *this;
}

IWString &
IWString::operator = (const IWString & other)
{
  if (_elements_allocated < other._number_elements)
    resize (other._number_elements);

  bcopy (other._things, _things, other._number_elements);

//for (int i = 0; i < other._number_elements; i++)
//  _things[i] = other._things[i];

  _number_elements = other._number_elements;

  return *this;
}

IWString &
IWString::operator = (ostrstream & zstream)
{
  IWString::strncpy (zstream.str (), zstream.pcount ());

  zstream.freeze (0);   // allow zstream to delete the buffer, we have a copy

  return * this;
}

int
const_IWSubstring::copy_to_char_array (char * s) const
{
  ::strncpy (s, _data, _nchars);
  s[_nchars] = '\0';

  return _nchars;
}

int
IWString::copy_to_char_array (char * s) const
{
  if (_number_elements)
    ::strncpy (s, _things, _number_elements);

  s[_number_elements] = '\0';

  return _number_elements;
}

const_IWSubstring
const_IWSubstring::substr (int cstart, int nchars) const
{
  if (nchars < 0)
    nchars = _nchars - cstart;

  assert (nchars >= 0 && nchars < _nchars);

  return const_IWSubstring (&(_data[cstart]), nchars);
}

int
const_IWSubstring::starts_with (char s) const
{
  if (0 == _nchars)
    return 0;

  return s == _data[0];
}

int
const_IWSubstring::starts_with (const char * s) const
{
  int lens = strlen (s);
  if (_nchars < lens)
    return 0;

  return 0 == ::strncmp (_data, s, lens);
}

int
const_IWSubstring::ends_with (char s) const
{
  if (0 == _nchars)
    return 0;

  return s == _data[_nchars - 1];
}

int
const_IWSubstring::ends_with (const char * s) const
{
  int lens = strlen (s);
  if (_nchars < lens)
    return 0;

  return 0 == ::strncmp (_data + _nchars - lens, s, lens);
}

int
const_IWSubstring::remove_leading_chars (int nremove)
{
  assert (nremove > 0 && nremove <= _nchars);

  _data += nremove;
  _nchars -= nremove;

  return _nchars;
}

int
const_IWSubstring::chop (int nchop)
{
  assert (nchop > 0 && nchop <= _nchars);

  return _nchars -= nchop;
}

int
const_IWSubstring::truncate (int n)
{
  assert (n >= 0 && n <= _nchars);

  return _nchars = n;
}

int
const_IWSubstring::truncate_at_first (char c)
{
  int i;     // BARF
  for (i = 0; i < _nchars; i++)
  {
    if (c == _data[i])
    {
      _nchars = i;
      return 1;
    }
  }

  return 0;
}

/*
  Cannot think of an easy way to do this with a lib routine.
  We do it the dumb (and broken for MAX_INT) way
*/

int
const_IWSubstring::is_int (int & result) const
{
  if (0 == _nchars)
    return 0;

  return string_class_is_int (_data, _nchars, result);
}

int
const_IWSubstring::numeric_value (float & result) const
{
  double tmp;

  if (! numeric_value (tmp))
    return 0;

// Check for overflow!!!

  result = tmp;

  return 1;
}

int
const_IWSubstring::numeric_value (double & result) const
{
  if (0 == _nchars)
    return 0;

  if (_nchars > 25)
    return 0;

  char buffer[30];

  strncpy (buffer, _data, _nchars);
  buffer[_nchars] = '\0';

  char * endptr;

  result = strtod (buffer, &endptr);
  if ('\0' == *endptr)
    return 1;

  return 0;
}

int
const_IWSubstring::numeric_value (int & result) const
{
  return string_class_is_int (_data, _nchars, result);
}

int
const_IWSubstring::index (const char c) const
{
  int i;     // BARF
  for (i = 0; i < _nchars; i++)
    if (c == _data[i])
      return i;

  return -1;
}

int
const_IWSubstring::rindex (const char c) const
{
  int i;     // BARF
  for (i = _nchars - 1; i > 0; i--)
    if (c == _data[i])
      return i;

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
const_IWSubstring::next (const char c, int & istart) const
{
  assert (istart >= 0 && istart < _nchars);

  return do_next (_data, _nchars, c, istart);
}

ostream &
operator << (ostream & os, const const_IWSubstring & ss)
{
  assert (os.good ());

  os.write (ss._data, ss._nchars);

  return os;
}

#include "stdio.h"

void
append_int (IWString & s, int n)
{
  char buffer[32];

  (void) sprintf (buffer, "%d", n);

  s += buffer;

  return;
}

static void
_cat (const char * lhs, int l1, const char * rhs, int l2,
     IWString & result)
{
  result.resize (l1 + l2 + 1);   // leave room for newline

  result.add (lhs, l1);
  result.add (rhs, l2);

  return;
}

IWString
operator + (const char * lhs, const IWString & rhs)
{
  int l1 = strlen (lhs);

  IWString rc;

  _cat (lhs, l1, rhs._things, rhs._number_elements, rc);

  return rc;
}

IWString
operator + (const IWString & lhs, const char * rhs)
{
  IWString rc;

  int l2 = strlen (rhs);

  _cat (lhs._things, lhs._number_elements, rhs, l2, rc);

  return rc;
}

IWString
operator + (const IWString & lhs, const IWString & rhs)
{
  IWString rc;

  _cat (lhs._things, lhs._number_elements, rhs._things, rhs._number_elements, rc);

  return rc;
}

void
append_digit (IWString & s, int digit)
{
  char buffer[32];

  sprintf (buffer, "%d", digit);

  s += buffer;

  return;
}

void
IWString::append_number (int znumber)
{
  int sign = 1;
  if (znumber < 0)
  {
    sign = -1;
    znumber = - znumber;
  }

  if (0 == znumber)
  {
    if (_number_elements == _elements_allocated)
      resize (_number_elements + 1);

    _things[_number_elements] = '0';
    _number_elements++;
    return;
  }

  char buffer[32];
  int ndigits = 0;
  while (znumber)
  {
    int j = znumber % 10;
    buffer[ndigits++] = '0' + j;
    znumber /= 10;
  }

  int chars_needed = ndigits + (sign < 0);

  if (_number_elements + chars_needed > _elements_allocated)
    resize (_number_elements + chars_needed);

  if (sign < 0)
    _things[_number_elements++] = '-';

  int i;     // BARF
  for (i = 0; i < ndigits; i++)
  {
    _things[_number_elements++] = buffer[ndigits - i - 1];
  }

  return;
}

void
IWString::append_number (float f)
{
  char buffer[20];

  sprintf (buffer, "%f", f);

  operator += (buffer);

  return;
}

void
IWString::append_number (float f, const char * fformat)
{
  char buffer[20];

  sprintf (buffer, fformat, f);

  operator += (buffer);

  return;
}

const_IWSubstring
substr (const IWString &s, int cstart, int nchars)
{ 
  if (nchars < 0)
    nchars = s.number_elements () - cstart;

  return s.substr (cstart, nchars);
}

const_IWSubstring
substr (const const_IWSubstring & s, int cstart, int nchars)
{
  if (nchars < 0)
    nchars = s._nchars - cstart;

  return s.substr (cstart, nchars);
}

const_IWSubstring
from_to (const IWString & s, int cstart, int cstop)
{
  return s.from_to (cstart, cstop);
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
  while (isspace (s[istart]) && istart < nchars)
    istart++;

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

  int i;     // BARF
  for (i = istart; i < nchars; i++)
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
  return _is_int (s._data, s._nchars, result);
}

/*int
is_int (const IWString & s, int & result)
{
  return _is_int (s._things, s._number_elements, result);
}*/

int
change_suffix (IWString & fname, const IWString & new_suffix)
{
  int len_fname = fname.length ();

  if (0 == len_fname)
    return 0;

  int last_period = -1;
  int i;     // BARF
  for (i = 0; i < len_fname; i++)
  {
    if ('.' == fname[i])
      last_period = i;
  }

  if (last_period < 0)
    fname += '.';
  else
    fname.shorten (last_period + 1);

  fname += new_suffix;

  return 1;
}

/*
  Shorten is a convenient way to shorten a string. If you are
  interested in storage, use resize ();
*/

int
IWString::shorten (int new_size)
{
  assert (new_size >= 0 && new_size <= _number_elements);

  _number_elements = new_size;

  return 1;
}

int
IWString::starts_with (char s) const
{
  if (0 == _number_elements)
    return 0;

  return (s == _things[0]);
}

int
IWString::starts_with (const char * s, int lens) const
{
  if (0 == _number_elements)
    return 0;

  if (lens < 0)
    lens = strlen (s);

  return (0 == ::strncmp (_things, s, lens));
}

int
IWString::starts_with (const IWString & s) const
{
  int lens = s._number_elements;

  return (0 == ::strncmp (_things, s._things, lens));
}

int
IWString::starts_with (const const_IWSubstring & s) const
{
  int lens = s._nchars;

  return (0 == ::strncmp (_things, s._data, lens));
}

int
IWString::translate (char cfrom, char cto)
{
  assert (ok ());

  int rc = 0;
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
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
IWString::ends_with (char s) const
{
  if (0 == _number_elements)
    return 0;

  return (s == _things[_number_elements - 1]);
}

int
IWString::ends_with (const char * s, int lens) const
{
  assert (s);

  if (lens < 0)
    lens = strlen (s);

  if (lens > _number_elements)
    return 0;

  if (0 == _number_elements)     // catch the case of both blank
    return 0;

  int i = _number_elements - lens;
  return (0 == ::strncmp (s, &(_things[i]), lens));
}

int
IWString::ends_with (const IWString & s) const
{
  if (s._number_elements > _number_elements)
    return 0;

  if (0 == _number_elements)     // catch the case of both being blank
    return 0;

  int i = _number_elements - s._number_elements;
  return (0 == ::strncmp (s._things, &(_things[i]), s._number_elements));
}

int
IWString::looks_like (const char * s, int min_chars_needed_for_match) const
{
  assert (ok ());
  assert (s);
  assert (min_chars_needed_for_match);

  if (_number_elements < min_chars_needed_for_match)
    return 0;

  int lens = strlen (s);

  if (_number_elements > lens) 
    return 0;

  return 0 == ::strncmp (_things, s, _number_elements);
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
  while (isalpha (*old_name))
  {
    buffer[iptr] = *old_name;
    old_name++;
    iptr++;
  }
  buffer[iptr] = '\0';

  const char * env = getenv (buffer);

  if (NULL == env)
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
IWString::insert (char c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  if (_number_elements == _elements_allocated)
    resize (_number_elements + 1);

  int i;     // BARF
  for (i = _number_elements; i >= pos; i--)
    _things[i] = _things[i - 1];

  _things[pos] = c;

  _number_elements++;

  return 1;
}

int
IWString::insert (const char * c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  int lens = strlen (c);

  return _insert (pos, c, lens);
}

int
IWString::insert (const IWString & c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  return _insert (pos, c._things, c._number_elements);
}

int
IWString::insert (const const_IWSubstring & c, int pos)
{
  assert (pos >= 0 && pos < _number_elements);

  return _insert (pos, c._data, c._nchars);
}

int
IWString::_insert (int pos, const char * s, int lens)
{
  if (_number_elements + lens > _elements_allocated)
    resize (_number_elements + lens);

  int i;     // BARF
  for (i = _number_elements + lens - 1; i >= pos + lens; i--)
    _things[i] = _things[i - lens];

  for (i = 0; i < lens; i++)
    _things[pos + i] = s[i];

  _number_elements += lens;

  return 1;
}

/*
  Locate the next occurrence of C. Increment istart to allow repeated calls.
*/

int
IWString::next (const char c, int & istart) const
{
  assert (ok_index (istart));

  return do_next (_things, _number_elements, c, istart);
}

/*
  This looks to be exactly the same as translate
*/

int
IWString::gsub (char cfrom, char cto, int how_many)
{
  int rc = 0;
  int i;     // BARF
  for (i = 0; i < _number_elements; i++)
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
  int i;     // BARF
  for (i = 0; i < nchars; i++)
    if (haystack[i] == needle)
      rc++;

  return rc;
}

int
const_IWSubstring::ccount (char c) const
{
  return _internal_ccount (_data, _nchars, c);
}

int
IWString::ccount (char c) const
{
  return _internal_ccount (_things, _number_elements, c);
}

int
IWString::strncpy (const char * s, int nchars)
{
  if (_elements_allocated < nchars + 1)
    resize (nchars);

  ::strncpy (_things, s, nchars);

  _number_elements = nchars;

  return 1;
}

int
IWString::strncat (const char * s, int nchars)
{
  if (_number_elements + nchars >= _elements_allocated)
    resize (_number_elements + nchars + 1);

  int i;     // BARF
  for (i = 0; i < nchars; i++)
    _things[_number_elements + i] = s[i];

  _number_elements += nchars;

  return _number_elements;
}

/*
  A post-fix increment operator
*/

char
const_IWSubstring::operator++ (int)
{
//cerr << "Substring operator ++ postfix, " << _nchars << " characters\n";
  if (0 == _nchars)
    return '\0';

//cerr << "Currently pointing at '" << *_data << "', will return '" << (*(_data + 1)) << "'\n";
  _data++;
  _nchars--;
  return *_data;
}

/*
  A pre-fix increment operator
*/

char
const_IWSubstring::operator++ ()
{
  if (0 == _nchars)
    return '\0';

  char rc = *_data;
  _data++;
  _nchars--;
  return rc;
}

void
const_IWSubstring::operator+= (int howfar)
{
  assert (howfar >= 0 && howfar <= _nchars);

  _data += howfar;
  _nchars -= howfar;

  return;
}

int
IWString::strspn (const char * s2)
{
  null_terminate ();

  return ::strspn (_things, s2);
}

/*
  Append the text in ZEXTRA
*/

IWString &
operator << (IWString & s, ostrstream & zextra)
{
  assert (s.ok ());

  s.strncat (zextra.str (), zextra.pcount ());

  zextra.freeze (0);

  return s;
}

/*
  We frequently want the part of a string before a character as well
  as the part after.

  If the character is not found, we return 0
  If the character is found at the beginning or end of the string,
  we return 1 and fill in either before or after
  If it is somewhere in the middle, we return 2
*/

int
IWString::split (const_IWSubstring & before,
                 char s,
                 const_IWSubstring & after) const
{
  before = "";
  after = "";

  int i = index (s);
  if (i < 0)
    return 0;

  if (0 == i)
  {
    from_to (1, _number_elements - 1, after);
    return 1;
  }

  if (i == _number_elements - 1)
  {
    from_to (0, _number_elements - 2, before);
    return 1;
  }

  from_to (0, i - 1, before);
  from_to (i + 1, _number_elements - 1, after);

  return 2;
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
operator << (IWString & s1, const char * s2)
{
  s1 += s2;

  return s1;
}
