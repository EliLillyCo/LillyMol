// Tester for iwstring library functions.
// This should be run with a debugging malloc

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring.h"
//#include "iwmalloc.h"
using std::cerr;
using std::cout;
using std::endl;

static int verbose = 0;
static int abort_on_error = 0;

void
handle_failure ()
{
  if (abort_on_error)
  {
    abort ();
  }
  else
    exit (1);
}

void
should_match (const char *result, const char *expected_result,
              const char *invoker)
{
  if (0 == strcmp (result, expected_result))
    return;

  cerr << invoker <<  " failed, actual '" << result << "', expected '" << expected_result << "'\n";
  handle_failure ();
}

void
numeric_match (int result, int expected_result, const char *buffer,
              const char *invoker)
{
  if (result == expected_result)
    return;

  cerr << "Failure on numeric match, result = " << result << " expected = "
       << expected_result << " buffer = '" << buffer << "', function '" <<
       invoker << "'\n";

  handle_failure ();
}

static int
test_is_double ()
{
  const char * buffer = "9";
  int i;
  double x;

  i = is_double (buffer, &x);
  numeric_match (i, 1, buffer, "is_double");
  
  buffer = "    .1";
  i = is_double (buffer, &x);
  numeric_match (i, 1, buffer, "is_double");

  buffer = "  3.9    ";
  i = is_double (buffer, &x);
  numeric_match (i, 1, buffer, "is_double");

  buffer = "  3.9000u    ";
  i = is_double (buffer, &x);
  numeric_match (i, 0, buffer, "is_double");

  return 0;
}

static int
test_ccount ()
{

  const char * buffer = "      a b c de f";
  char needle = 'g';
  int i = ccount (buffer, needle);
  numeric_match (i, 0, buffer, "ccount");

  needle = 'A';
  buffer = "FLDJFAAA:LLAALKLKA";
  i = ccount (buffer, needle);
  numeric_match (i, 6, buffer, "ccount");

  return 0;
}

int
test_remove_blanks ()
{
  char buffer[80];
  strcpy (buffer, "abcdef");
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  strcpy (buffer, "      abcdef");
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  strcpy (buffer, "   abcdef  ");
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  strcpy (buffer, "a b c d      ef ");
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  strcpy (buffer, "  a b cd      ef ");
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  return 0;
}

void
do_words_test (int test_number, const char *buffer, 
               const char separator, int expected_result)
{
  int result = words (buffer, separator);

  if (result == expected_result)
    return;

  cerr << "Words test " << test_number << " failed\n";

  (void) fprintf (stderr, "Got %d, should have been %d\n", result, expected_result);

  (void) fprintf (stderr, "String is '%s', separator is '%c'\n", buffer, separator);

  return;
}

void
test_words ()
{
  char buffer[512];
  buffer[0] = '\0';

  do_words_test (1, "", ' ', 0);
  do_words_test (2, " ", ' ', 0);

  do_words_test (3, " ", '1', 1);
  do_words_test (4, "$", '$', 0);

  do_words_test (5, "    ", '@', 1);

  do_words_test (6, "  $  ", '$', 2);

  do_words_test (7, "!   a  !", '!', 1);

  return;
}

/*
  Driver for string testing routines.
  -v option prints verbose output
  -d option puts iwmalloc into debug mode
*/
int
main (int argc, char **argv)
{

  int o;
  while ((o = getopt(argc, argv, "vd")) != EOF)
    switch (o) {
      case 'v':
        verbose = 1;
        break;

//    case 'd':
//      iwmalloc_set_debug (1);
//      break;

      case 'a':
        abort_on_error = 1;
        break;

      default:
        cerr << "Unrecognised option '" << o << "'\n";
        return 1;
    }

  if (verbose)
    cout << "String tester begins\n";

  char buffer[512];
  strcpy (buffer, "      abcdef");

  if (verbose)
    cout << "Testing strip leading blanks on '" << buffer << "'\n";

  strip_leading_blanks (buffer);
  should_match (buffer, "abcdef", "strip_leading_blanks");

  strcpy (buffer, "      ");
  if (verbose)
    cerr << "Testing strip_leading_blanks on '" << buffer << "'\n";
  strip_leading_blanks (buffer);
  should_match (buffer, "", "strip_leading_blanks");

  strcpy (buffer, "     1 23   4   ");
  if (verbose)
    cout << "Testing strip leading blanks on '" << buffer << "'\n";
  strip_leading_blanks (buffer);
  should_match (buffer, "1 23   4   ", "strip_leading_blanks");

  strcpy (buffer, "      ");
  if (verbose)
    cerr << "Testing strip_trailing_blanks on '" << buffer << "'\n";
  strip_trailing_blanks (buffer);
  should_match (buffer, "", "strip_trailing_blanks");

  strcpy (buffer, "  abcdef      ");
  if (verbose)
    cerr << "Test strip_trailing_blanks on '" << buffer << "'\n";
  strip_trailing_blanks (buffer);
  should_match (buffer, "  abcdef", "strip_trailing_blanks");

  strcpy (buffer, "a  bdef");
  if (verbose)
    cerr << "Test strip_trailing_blanks on '" << buffer << "\n";
  strip_trailing_blanks (buffer);
  should_match (buffer, "a  bdef", "strip_trailing_blanks");

  strcpy (buffer, "foo");
  if (verbose)
    cout << "Testing basename on '" << buffer << "'\n";
  const char *c = iwbasename (buffer);
  should_match (c, "foo", "basename");

  strcpy (buffer, "foo/bar");
  if (verbose)
    cout << "Testing basename on '" << buffer << "'\n";
  c = iwbasename (buffer);
  should_match (c, "bar", "basename");

  strcpy (buffer, "/usr/local/");
  if (verbose)
    cout << "testing basename on '" << buffer << "'\n";
  c = iwbasename (buffer);
  should_match (c, "local/", "basename");

  strcpy (buffer, " ABCdEF123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?");
  if (verbose)
    cout << "Testing lowercase on '" << buffer << "'\n";
  to_lowercase (buffer);
  should_match (buffer, " abcdef123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?", "to_lowercase");

  strcpy (buffer, " abcDef123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?");
  if (verbose)
    cout << "Testing uppercase on '" << buffer << "'\n";
  to_uppercase (buffer);
  should_match (buffer, " ABCDEF123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?", "to_uppercase");

  test_remove_blanks ();

  if (verbose)
    cout << "Testing ccount\n";
  test_ccount ();

  if (verbose)
    cout << "Testing is double\n";
  test_is_double ();

  if (verbose)
    cout << "Testing words\n";
  test_words ();

  cerr << "All tests successful\n";

  return 0;
}
