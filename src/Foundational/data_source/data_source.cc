#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
//#include <String.h>

#include "data_source.h"
#include "iwstring.h"

void
data_source::_default_values ()
{
  _fp = NULL;
  _buffer = NULL;
  _buf_size = 0;
  _record_buffered = 0;
  _longest_record = 0;
  _lines_read = 0;
  _at_eof = 0;
  _strip_newlines = 0;
  _strip_leading_blanks = 0;
  _strip_trailing_blanks = 0;
  _skip_blank_lines = 0;
  _ignore_pattern = NULL;
  _strlen_ignore_pattern = 0;
  _filter_pattern = NULL;
  _strlen_filter_pattern = 0;
  _match_filter_at_start = 0;
  _filter_match_position = NULL;
  _created_from_fp = 0;
  _convert_to_lowercase = 0;
  _convert_to_uppercase = 0;
}

data_source *
new_data_source (FILE *fp, int lrecl)
{
  assert (NULL != fp);
  assert (lrecl > 1);

  data_source *tmp = new data_source (lrecl);
  tmp->_fp = fp;
  tmp->_created_from_fp = 1;

  return tmp;
}

data_source *
new_data_source (const char *fname, int lrecl)
{
  assert (NULL != fname);
  assert (lrecl > 1);

  FILE *fp = fopen (fname, "r");
  if (NULL == fp)
  {
    (void) fprintf (stderr, "new data source: cannot open '%s'\n", fname);
    perror ("");
    return NULL;
  }

  return new_data_source (fp, lrecl);
}

/*
  A data source is created with a suggested starting size for the buffer
*/

data_source::data_source (int lrecl)
{
  assert (lrecl > 1);

  _default_values ();

  _buffer = new char[lrecl];
  _buf_size = lrecl;
  _buffer[0] = '\0';

  return;
}

#define DELETE_IF_NOT_NULL(p) { if (NULL != (p)) delete (p); }

data_source::~data_source ()
{
  DELETE_IF_NOT_NULL(_buffer);

  if (! _created_from_fp)
    fclose (_fp);

  DELETE_IF_NOT_NULL(_ignore_pattern);
  DELETE_IF_NOT_NULL(_filter_pattern);
}

int
data_source::_line_is_blank () const
{
  assert (NULL != _buffer);

  char *c = _buffer;
  while (*c)
  {
    if (! isspace (*c))
      return 0;
    c++;
  }

  return 1;
}

int
data_source::push_record ()
{
  assert (NULL != _buffer);

  _record_buffered = 1;

  return 1;    // why does this function return anything ?
}

/*
  In the interests of clarity (rather than efficiency), we always
  delete any existing pattern - even if the new one could fit in
  that space.
*/

void
data_source::set_ignore_pattern (const char *pattern)
{
  DELETE_IF_NOT_NULL (_ignore_pattern);

  _ignore_pattern = NULL;
  _strlen_ignore_pattern = 0;

  if (NULL == pattern)
    return;

  int i = strlen (pattern);
  if (0 == i)
    return;

  _ignore_pattern = new char[i + 1];
  _strlen_ignore_pattern = i;
  strcpy (_ignore_pattern, pattern);
}

void
data_source::set_filter_pattern (const char *pattern)
{
  if (NULL != _filter_pattern)
    delete _filter_pattern;

  _filter_pattern = NULL;
  _strlen_filter_pattern = 0;
  _match_filter_at_start = 0;
  _filter_match_position = NULL;

  if (NULL == pattern)
    return;

  int i = strlen (pattern);

  if (0 == i)
    return;

  if ('^' == *pattern)    // match pattern at start of line only
  {
    assert (i > 1);      // "^" by itself is not accepted
    _match_filter_at_start = 1;
    _filter_pattern = new char [i];
    _strlen_filter_pattern = i - 1;
    strcpy (_filter_pattern, pattern + 1);
  }
  else
  {
    _filter_pattern = new char [i + 1];
    _strlen_filter_pattern = i;
    strcpy (_filter_pattern, pattern);
  }
}

/*
  This function applies all filters which are set to _buffer
  If it passes all of them, we return 1
  The order in which these filters are applied can be significant.
  Change this order and break existing code.
*/

int
data_source::_apply_all_filters ()
{
  assert (NULL != _buffer);

  if (_strip_trailing_blanks)
    strip_trailing_blanks (_buffer);

// Here is where we can get an order dependency. This is quite
// ambiguous. For now, the code defines the behaviour.

  if (_strip_leading_blanks)
    strip_leading_blanks (_buffer);

  if (_skip_blank_lines && _line_is_blank ())
    return 0;

// Do case conversion before we apply string based filters.

  if (_convert_to_lowercase)
    to_lowercase (_buffer);

  if (_convert_to_uppercase)
    to_uppercase (_buffer);

  if (_strlen_ignore_pattern && 
      0 == strncmp (_buffer, _ignore_pattern, _strlen_ignore_pattern))
    return 0;

  if (0 == _strlen_filter_pattern)     // no filter pattern
    ;
  else if (_match_filter_at_start &&
           0 != strncmp (_buffer, _filter_pattern, _strlen_filter_pattern))
    return 0;
  else if (NULL == (_filter_match_position = strstr (_buffer, _filter_pattern)))
    return 0;

  return 1;
}

/*
  Frequently a programme will need to read down until it gets a record
  which matches a pattern. This function does that, and leaves the record
  buffered.
*/

int
data_source::next_record_matches (const char *pattern)
{
  assert (NULL != pattern);

  if (! _record_buffered)
  {
    const char *c = next_record ();
    if (NULL == c)     // reached EOF
      return 0;
  }

  _record_buffered = 1;

  int strlen_pattern = strlen (pattern);
  return 0 == strncmp (_buffer, pattern, strlen_pattern);
}

/*
  Whenever a complete record is obtained, there are some tasks to be
  performed.
  c is a pointer to the newline character in _buffer
*/

int
data_source::_got_complete_record (char *c)
{
  assert ('\n' == *c);

  _lines_read++;
  if (c - _buffer > _longest_record)
    _longest_record = c - _buffer;

  if (_strip_newlines)
    *c = '\0';

  return 1;
}

/*
  Whenever an EOF is encountered, this function is called to clean up.
  We return 1 if there is any data in the buffer
*/

int
data_source::_handle_eof_processing ()
{
  _at_eof = 1;

  if ('\0' == _buffer[0])
    return 0;

  (void) fprintf (stderr, "data source: unterminated last line\n");
  int strlen_buffer = strlen (_buffer);
  if (strlen_buffer > _longest_record)
    _longest_record = strlen_buffer;

  return 1;
}

/*
  this function fills _buffer with the next record
  Returns 1 if successful
*/

int
data_source::_fetch_record ()
{
  assert (NULL != _buffer);

  _buffer[0] = '\0';
  if (NULL == fgets (_buffer, _buf_size - 1, _fp))
    return _handle_eof_processing ();

  char *c = strchr (_buffer, '\n');
  if (NULL != c)
    return _got_complete_record (c);

/*
  We did not get a complete record at our first attempt, we must
  grow the buffer until we can fit this record.
*/

  while (1)
  {
    char * nb = new char[2 * _buf_size];
    strcpy (nb, _buffer);

    char *rc = fgets (_buffer, _buf_size - 1, _fp);
    strcat (nb, _buffer);

    delete _buffer;
    _buffer = nb;
    _buf_size = 2 * _buf_size;

    if (NULL == rc)
      return _handle_eof_processing ();

    if (NULL != (c = strchr (_buffer, '\n')))
      return _got_complete_record (c);
  }

  abort ();     // should not come to here
  return 0;
}

const char *
data_source::most_recent_record ()
{
  assert (_lines_read > 0);

  return _buffer;
}

const char *
data_source::next_record ()
{
  assert (NULL != _buffer);

  if (_record_buffered)
  {
    _record_buffered = 0;
    if (_apply_all_filters ())
      return _buffer;
  }

/*
  This infinite loop will loop until either we get a record that
  matches the filters, or we hit EOF
*/

  _record_buffered = 0;
  while (1)
  {
    if (! _fetch_record ())
      return NULL;

    if (_apply_all_filters ())
      return _buffer;
  }

  abort ();     // should never come here
  return NULL;
}

/*
  Reads records until we get to one containing pattern.
  If pattern starts with '^', it must match at the beginning
  of the record.
  We return 1 if successful, 0 otherwise
  The next record read after this call, will be the record after
  the one which contains the pattern
*/

int
data_source::skip_past (const char *pattern)
{
  assert (NULL != pattern);

  if (_at_eof)
    return 0;

// local copy of pattern.

  const char *lpattern;

// do we match at the beginning or anywhere

  int match_anywhere;

  if ('^' == *pattern)
  {
    lpattern = pattern + 1;
    match_anywhere = 0;
  }
  else
  {
    match_anywhere = 1;
    lpattern = pattern;
  }

  assert (0 != strlen (lpattern));

  const char *buffer;

  while (NULL != (buffer = next_record ()))
  {
    const char *c = strstr (buffer, lpattern);
    if (NULL != c)
    {
      if (match_anywhere)
        return 1;

      if (c == buffer)     // match at beginnign only
        return 1;

    }
  }

  return 0;
}

/*
  Very similar to skip_past, but we push the record when we are done
  so that the next record read will be one with the pattern
*/

int
data_source::skip_to (const char *pattern)
{
  if (1 == skip_past (pattern))
  {
    push_record ();
    return 1;
  }
  else
    return 0;
}
