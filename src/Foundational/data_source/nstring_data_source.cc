#include <stdlib.h>
#include <libgen.h>
#include <assert.h>

#include "String_Data_Source.h"
#include "iwstring.h"

void
String_Data_Source::_default_values (int lrecl)
{
  _record_buffered = 0;
  _longest_record = 0;
  _record_delimiter = '\n';
  _lines_read = 0;
  _lines_which_are_returned = 0;
  _strip_leading_blanks = 0;
  _strip_trailing_blanks = 0;
  _compress_spaces  = 0;
  _skip_blank_lines = 0;
  _ignore_pattern_string = "";
  _ignore_pattern_regexp = NULL;
  _filter_pattern_string = "";
  _filter_pattern_regexp = NULL;
  _filter_match_position = -1;
  _convert_to_lowercase = 0;
  _convert_to_uppercase = 0;

  _buffer.alloc (lrecl);
  _buf_size = lrecl;
  _input_buffer = new char[_buf_size];

  return;
}

String_Data_Source::String_Data_Source ()
{
  assert (NULL != fname);

  _default_values ();

  return;
}

String_Data_Source::String_Data_Source (const char * fname, int lrecl)
{
  _default_values (lrecl);

  ifstream::open (fname, ios::in);

  return;
}

#define DELETE_IF_NOT_NULL(p) { if (NULL != (p)) delete (p); }

String_Data_Source::~String_Data_Source ()
{
  DELETE_IF_NOT_NULL (_ignore_pattern_regexp);
  DELETE_IF_NOT_NULL (_filter_pattern_regexp);
  DELETE_IF_NOT_NULL (_input_buffer);
}

int
String_Data_Source::ok () const
{
  if (NULL == _input_buffer)
    return 0;

  if (eof ())    // the data source is still valid, although at EOF
    return 1;

  if (! good ())
    return 0;

  return 1;
}

int
String_Data_Source::debug_print (ostream & os) const
{
  assert (os.good ());

  if (good ())
    os << "OSTREAM is good, currently at " << tellg () << endl;

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

  if (NULL != _filter_pattern_regexp)
    os << "Will filter lines '" << _filter_pattern_string << "\n";
  else
    os << "No filter pattern\n";

  if (NULL != _ignore_pattern_regexp)
    os << "Will ignore lines '" << _ignore_pattern_string << "\n";
  else
    os << "No ignore pattern\n";

  if (_record_buffered)
    os << "There is a record buffered\n" << _buffer << endl;
  else
    os << "No record buffered\n";

  return 1;
}

int
String_Data_Source::push_record ()
{
  assert (ok ());

  if (_record_buffered)
    return 0;    // unlikely anyone will ever check this

  _record_buffered = 1;

  return 1;
}

static char *
compile_regular_expression (const char * pattern)
{
  char * rc = regcmp (pattern, (char *) 0);
  if (rc)
    return rc;

  cerr << "compile_regular_expression:: regcmp failed '" << pattern << "'\n";
  assert (NULL == "Bad luck");
  return NULL;
}

/*
  In the interests of clarity (rather than efficiency), we always
  delete any existing pattern - even if the new one could fit in
  that space.
*/

void
String_Data_Source::set_ignore_pattern (const char * pattern)
{
  DELETE_IF_NOT_NULL (_ignore_pattern_regexp);

  _ignore_pattern_string = pattern;
  _ignore_pattern_regexp = compile_regular_expression (pattern);
}

void
String_Data_Source::set_filter_pattern (const char * pattern)
{
  DELETE_IF_NOT_NULL (_filter_pattern_regexp);

  _filter_pattern_string = pattern;
  _filter_pattern_regexp = compile_regular_expression (pattern);
}

/*
  This function applies all filters which are set to _buffer
  If it passes all of them, we return 1
  The order in which these filters are applied can be significant.
  Change this order and break existing code.
*/

#include <ctype.h>

int
String_Data_Source::_apply_all_filters ()
{
  if (_strip_trailing_blanks)
    strip_trailing_blanks (_buffer);

  if (_strip_leading_blanks && _buffer.length () && isspace (_buffer[0]))
    strip_leading_blanks (_buffer);

  if (_skip_blank_lines && 0 == _buffer.length ())
    return 0;

  if (_compress_spaces && _buffer.length ())
    compress_blanks (_buffer);

// Do case conversion before we apply string based filters.

  if (_convert_to_lowercase)
    to_lowercase (_buffer);

  if (_convert_to_uppercase)
    to_uppercase (_buffer);

// If this record matches the ignore pattern, reject this line

  if (_ignore_pattern_regexp)
  {
    if (NULL != regex (_ignore_pattern_regexp, _buffer.data ()))
      return 0;
  }

// If we don't match the filter, reject this record.

  if (NULL != _filter_pattern_regexp)
  {
    if (0 == _buffer.length ())
      return 0;

    if (NULL == regex (_filter_pattern_regexp, _buffer.data ()))
      return 0;
  }

  _lines_which_are_returned++;
  return 1;
}

/*
  Frequently a programme will need to read down until it gets a record
  which matches a pattern. This function does that, and leaves the record
  buffered.
*/

int
String_Data_Source::next_record_matches (const char * pattern)
{
  if (! ok () || ! good ())
    return 0;

  if (0 == strlen (pattern))   // what about EOF?
    return 0;

  if (! _record_buffered)
  {
    _fetch_record ();
    if (bad ())
      return 0;
  }

  push_record ();

  char * compiled_regular_expression = compile_regular_expression (pattern);
  assert (compiled_regular_expression);

  int rc = 0;
  if (NULL == regex (compiled_regular_expression, _buffer.data ()))
    rc = 0;
  else
    rc = 1;

  delete compiled_regular_expression;

  return rc;
}

/*
  Return non zero if we successfully read a record.
*/

int
String_Data_Source::_fetch_complete_record (int first)
{
  if (! good ())
    return 0;

  get (_input_buffer, _buf_size, _record_delimiter);
  int nchars = gcount ();
//cerr << "Just read " << nchars << " characters\n";

  if (nchars)
  {
    if (first)
      _buffer = _input_buffer;
    else
      _buffer += _input_buffer;
  }
  else if (first)
    _buffer = "";

  if (eof ())
  {
    if (nchars > 0)
      cerr << "Incomplete last record\n";
    return nchars;
  }

  if (_record_delimiter == peek ())
  {
    seekg (tellg () + 1);
    return 1;       // we have a complete record.
  }

  return nchars + _fetch_complete_record (0);
}

/*
  this function fills _buffer with the next record
  Returns 1 if successful
*/

int
String_Data_Source::_fetch_record ()
{
  int nchars = _fetch_complete_record (1);

  if (0 == nchars)
  {
    _open = 0;
    _buffer = "";
    return 0;
  }

  _lines_read++;

  if (nchars > _longest_record)
    _longest_record = nchars;

  return 1;
}

int
String_Data_Source::most_recent_record (string & buffer)
{
  assert (_lines_read > 0);

  buffer = _buffer;

  return 1;
}

int
String_Data_Source::next_record (string & buffer)
{
  if (! good ())
  {
    cerr << "next_record: istream is not good\n";
    return 0;
  }

  if (eof ())
    return 0;

  if (_record_buffered)
  {
    _record_buffered = 0;
    if (_apply_all_filters ())
    {
      buffer = _buffer;
      return 1;
    }
  }

// This infinite loop will loop until either we get a record that
// matches the filters, or we hit EOF

  _record_buffered = 0;
  while (1)
  {
    if (! _fetch_record ())
      return 0;

    if (_apply_all_filters ())
    {
      buffer = _buffer;
      return 1;
    }
  }

  assert (NULL == "should never come here");
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
string_data_source::skip_past (const char * pattern)
{
  assert (ok ());

  if (0 == strlen (pattern))
    return 0;

  char * compiled_regular_expression = compile_regular_expression (pattern);

  int rc = _skip_past (compiled_regular_expression);

  delete compiled_regular_expression;

  return rc;
}

int
string_data_source::_skip_past (const char * compiled_regular_expression)
{
  while (_fetch_record ())
    if (regex (compiled_regular_expression, _buffer.data ()))
      return 1;

  return 0;
}

/*
  Very similar to skip_past, but we push the record when we are done
  so that the next record read will be one with the pattern
*/

int
String_Data_Source::skip_to (const char * pattern)
{
  assert (ok ());

  if (skip_past (pattern))
  {
    push_record ();
    return 1;
  }
  else
    return 0;
}
