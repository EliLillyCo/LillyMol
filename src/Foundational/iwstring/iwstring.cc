/*
  Strips leading spaces
*/

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring.h"

void
strip_leading_blanks (char *cc)
{
  char *c = cc;
  while (isspace (*c) && '\0' != *c)
    c++;

  if (c == cc)    /* no whitespace found */
    return;
  while (*cc)
  {
    *cc++ = *c++;
  }
  *cc = '\0';
}

void
strip_trailing_blanks (char *cc)
{
  size_t i = strlen (cc);
  char *c = cc + i - 1;

  if (0 == i)
    return;

  while (isspace (*c))
  {
    if (c == cc)
    {
      *c = '\0';
      return;
    }
    c--;
  }
  
  c++;
  *c = '\0';
  return;
}

void
no_newline (char *cc)
{
  char *c = strchr (cc, '\n');
  if (nullptr != c)
    *c = '\0';
  return;
}

#ifdef CONFLICTS_WITH_STDLIB
char *
basename(char *path_name)
{
  char *c, *file_name;

  c = file_name = path_name;
  while (*c)
  {
    if (*c == '/')
      file_name = c + 1;
    c++;
  }

  return (file_name);
}
#endif

int
remove_blanks (char *cc)
{
  char *c1, *c2;
  int nblanks = 0;

  c1 = c2 = cc;

  while (*c1)
  {
    if (' ' == *c1)
    {
      nblanks++;
      c1++;
    }
    else
    {
      *c2 = *c1;
      c1++;
      c2++;
    }
  }
  *c2 = '\0';
  return nblanks;
}

void
to_lowercase (char *cc, int characters_to_convert)
{
  int characters_processed = 0;
  while (*cc)
  {
    *cc = tolower (*cc);
    cc++;
    characters_processed++;
    if (characters_to_convert && characters_processed >= characters_to_convert)
      return;
  }

  return;
}

/*void
to_lowercase (resizable_array_p<char> *buffers)
{
  for (int i = 0; i < buffers->number_elements (); i++)
    to_lowercase (buffers->item (i));
}*/

void
to_uppercase (char *cc, int characters_to_convert)
{
  int characters_processed = 0;
  while (*cc)
  {
    *cc = toupper (*cc);
    cc++;
    characters_processed++;
    if (characters_to_convert > 0 && characters_processed >= characters_to_convert)
      return;
  }

  return;
}

/*void
to_uppercase (resizable_array_p<char> *buffers)
{
  for (int i = 0; i < buffers->number_elements (); i++)
    to_uppercase (buffers->item (i));
}*/

int
compress_blanks (char * zstring)
{
  int from = 0;
  int to = 0;

  int rc = 0;
  int previous_was_space = 0;
  while (zstring[from])
  {
    int copy_char = 1;

    if (isspace (zstring[from]))
    {
      if (previous_was_space)
      {
        copy_char = 0;
        rc++;
      }
      previous_was_space = 1;
    }
    else
      previous_was_space = 0;

    if (copy_char)
    {
      zstring[to] = zstring[from];
      to++;
    }

    from++;
  }

  zstring[to] = '\0';

  return rc;
}

/*
  We return the length of the string
*/

int
compress_blanks (char * zstring, int nchars)
{
  int from = 0;
  int to = 0;

  int previous_was_space = 0;
  for (int i = 0; i < nchars; i++)
  {
    int copy_char = 1;

    if (isspace (zstring[from]))
    {
      if (previous_was_space)
        copy_char = 0;
      previous_was_space = 1;
    }
    else
      previous_was_space = 0;

    if (copy_char)
    {
      zstring[to] = zstring[from];
      to++;
    }

    from++;
  }

  return to;
}

int
is_int (const char *buffer, int *i)
{
  char *c;
  int tmp = strtol (buffer, &c, 10);

  if (c == buffer)
    return 0;

  if ('\0' == *c)
    ;     // good
  else if (! isspace (*c))
    return 0;

  *i = tmp;
  return 1;
}

int
is_double (const char *buffer, double *x)
{
  assert (nullptr != buffer);
  assert (nullptr != x);
  char *c;
  double tmp = strtod (buffer, &c);

//cerr << "Parsing '" << buffer << "' as double\n";

  if (c == buffer)
    return 0;

// c must be pointing to a string terminator, or whitespace

  if ('\0' == *c)
    ;     // good
  else if (! isspace (*c))
    return 0;

  *x = tmp;
  return 1;
}

/*
  Counts the number of occurrences of NEEDLE in HAYSTACK
*/

int
ccount (const char *haystack, char needle)
{
  assert (nullptr != haystack);

  int count = 0;

  while (*haystack) 
  {
    if (needle == *haystack)
      count++; 
    haystack++;
  }
    
  return count;
} 

char *
make_copy (const char *c)
{
  assert (nullptr != c);
  int lenc = static_cast<int>(strlen (c));

  char * copyc = new char[lenc + 1];
  memcpy (copyc, c, lenc);

  return copyc;
}

/*#ifdef __GNUG__

template class resizable_array_p<char>;   // instantiate the template
template class resizable_array_base<char *>;   // instantiate the template

#endif*/

/*
  This function returns a resizable array of char *, consisting of
  all the tokens in buffer. strtok is used with PATTERN to divide
  the line.
*/

/*resizable_array_p<char> *
tokenise_char (const char *buffer, const char *pattern)
{
  assert (nullptr != buffer);
  assert (nullptr != pattern);

  int lenbuf = strlen (buffer);
  if (0 == lenbuf)
    return nullptr;

  char * copy_buf = new char[lenbuf + 1];
  strcpy (copy_buf, buffer);

  char *c = strtok (copy_buf, pattern);

  if (nullptr == c)
    return nullptr;

  resizable_array_p<char> *a = new resizable_array_p<char> (make_copy (c));

  while (nullptr != (c = strtok (nullptr, pattern)))
  {
    a->add (make_copy (c));
  }

  delete copy_buf;

  return a;
}*/

/*
  Tokenise as double and tokenise as int are identical except
  that one calls is_int and the other calls is_double.
  This template can instantiate either one.
*/

/*template <typename T>
resizable_array_p<T> *
tokenise_as (const char *buffer, const char *pattern,
             int (* is_t) (const char *, T *))
{
  assert (nullptr != buffer);
  assert (nullptr != pattern);

  resizable_array_p<char> * tk = tokenise (buffer, pattern);

  if (nullptr == tk)
    return nullptr;

  if (0 == tk->number_elements ())
  {
    delete tk;
    return nullptr;
  }

  resizable_array_p<T> * r = new resizable_array_p<T>;

  for (int i = 0; i < tk->number_elements (); i++)
  {
    T j;
    if (is_t (tk->item (i), &j))
    {
      T *jp = new int;
      *jp = j;
      r->add (jp);
    }
    else
    {
      r->add (nullptr);
    }
  }

  return r;
}*/

/*
  This doesn't work!!!
  Leave commented out until it does.

resizable_array_p<int> *
tokenise_as_int (const char *buffer, const char *pattern)
{
  int i;
  return tokenise_as (buffer, pattern, i, is_int);
}

resizable_array_p<double> *
tokenise_as_int (const char *buffer, const char *pattern)
{
  return tokenise_as (buffer, pattern, is_double);
}
*/

/*int
find_string (const resizable_array_p<char> * strings, const char *pattern)
{
  for (int i = 0; i < strings->number_elements (); i++)
    if (0 == strcmp (pattern, strings->item (i)))
      return i;

  return -1;
}*/

int
words (const char *string, const char separator)
{

  int in_word;

  if ('\0' == *string)
    return 0;        // null string has no words

  if (*string == separator)
    in_word = 0;
  else
    in_word = 1;

  int words = 0;

  while (1)
  {
    string++;

    if ('\0' == *string)
    {
      if (in_word)
        words++;
      return words;
    }

    if (separator == *string)
    {
      if (in_word)
      {
        words++;
        in_word = 0;
      }
    }
    else
      in_word = 1;
  }
}

const char *
iwbasename(const char * fname)
{
  const char * rc = fname;
  int previous_char_was_slash = 0;

  while (*fname)
  {
    if ('/' == *fname)
      previous_char_was_slash = 1;
    else
    {
      if (previous_char_was_slash)
        rc = fname;
      previous_char_was_slash = 0;
    }
    fname++;
  }

  return rc;
}
