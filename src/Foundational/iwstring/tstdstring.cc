// Tester for iwstring library functions.
// This should be run with a debugging malloc

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ostream>

#include <string>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring.h"
#include "iwmalloc.h"

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
should_match (const string & result, const char *expected_result,
              const char *invoker)
{
  if (result == expected_result)
    return;

  cerr << invoker << " failed, actual '" << result << "', expected '" <<
          expected_result << "'\n";

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
test_expand_env ()
{
  if (0 != putenv ("FOO=barf"))
  {
    cerr << "putenv failed\n";
    handle_failure ();
  }

  string expanded;

  if (! expand_environment_variables ("$FOO/bar", expanded))
  {
    cerr << "expand_environment_variables failed\n";
    handle_failure ();
  }

  should_match (expanded, "barf/bar", "test_expand_env");

  return 0;
}

static int
test_tokenise ()
{
  const char *buffer = "";
  resizable_array_p<string> *tk = tokenise (buffer, " ");
  if (NULL != tk)
  {
    cerr << "tokenise failed on empty string\n";
    handle_failure ();
  }

  buffer = "                    ";
  tk = tokenise (buffer, " ");
  if (NULL != tk)
  {
    cerr << "tokenise failed on blank string\n";
    handle_failure ();
  }

  buffer = "asdfasdf";
  if (verbose)
    cout << "Testing tokenise on '" << buffer << "'\n";
  tk = tokenise (buffer, " ");
  if (NULL == tk)
  {
    cerr << "tokenise failed on '" << buffer << "'\n";
    handle_failure ();
  }
  if (1 != tk->number_elements ())
  {
    cerr << "Bad tokenise token count " << tk->number_elements () << ", '" <<
             buffer << "'\n";
    handle_failure ();
  }

  delete tk;

  buffer = "abc def,ghi ";
  if (verbose)
    cout << "Testing tokenise on '" << buffer << "'\n";

  tk = tokenise (buffer, " ,");
  if (3 != tk->number_elements ())
  {
    cerr << "Bad tokenise count " << tk->number_elements () <<
            " string was '" << buffer << "'\n";
    handle_failure ();
  }

  if (0 == tk->item (0)->compare ("abc") &&
      0 == tk->item (1)->compare ("def") &&
      0 == tk->item (2)->compare ("ghi") )
    ;
  else
  {
    cerr << "Bad tokenisation of '" << buffer << "'" << *(tk->item (0)) <<
            "' '" << *(tk->item (1)) << "' '" << *(tk->item (2)) << "'\n";
    handle_failure ();
  }
  delete tk;

  return 0;
}

int
test_wordindex ()
{
  string buffer ("abc def");
  int i = wordindex (buffer, 1);
  numeric_match (i, 0, buffer.data (), "wordindex");

  i = wordindex (buffer, 2);
  numeric_match (i, 4, buffer.data (), "wordindex");

  i = wordindex (buffer, 8);
  numeric_match (i, -1, buffer.data (), "wordindex");

  buffer = " ";
  i = wordindex (buffer, 1);
  numeric_match (i, -1, buffer.data (), "wordindex");

  buffer = "  abc      def g";
  i = wordindex (buffer, 1);
  numeric_match (i, buffer.find ('a'), buffer.data (), "wordindex");

  i = wordindex (buffer, 2);
  numeric_match (i, buffer.find ('d'), buffer.data (), "wordindex");

  i = wordindex (buffer, 3);
  numeric_match (i, buffer.find ('g'), buffer.data (), "wordindex");

  return 1;
}

int
test_get_word ()
{
  string buffer ("abc");

  string word;

  int i = get_word (buffer, word, 1);
  should_match (word, "abc", "test_get_word");

  i = get_word (buffer, word, 2);
  if (i)
    cerr << "get word failure\n";

  buffer = "abc def";
  i = get_word (buffer, word, 1);
  should_match (word, "abc", "test_get_word");

  i = get_word (buffer, word, 2);
  should_match (word, "def", "test_get_word");

  buffer = " abc def ";
  i = get_word (buffer, word, 1);
  should_match (word, "abc", "test_get_word");

  i = get_word (buffer, word, 2);
  should_match (word, "def", "test_get_word");

  return 1;
}


int
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

int
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
  char * buffer = "abcdef";
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  buffer = "      abcdef";
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  buffer = "   abcdef  ";
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  buffer = "a b c d      ef ";
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  buffer = "  a b cd      ef ";
  remove_blanks (buffer);
  should_match (buffer, "abcdef", "remove_blanks");

  return 0;
}

void
test_change_suffix ()
{
  string foo = "bar";
  change_suffix (foo, "bar");

  should_match (foo, "bar.bar", "change suffix");

  change_suffix (foo, "foo");
  should_match (foo, "bar.foo", "change suffix");

  return;
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

      case 'd':
        iwmalloc_set_debug (1);
        break;

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
  string zstring = buffer;

  if (verbose)
    cout << "Testing strip leading blanks on '" << buffer << "'\n";

  strip_leading_blanks (buffer);
  should_match (buffer, "abcdef", "strip_leading_blanks");

  strip_leading_blanks (zstring);
  should_match (zstring, "abcdef", "strip_leading_blanks_string");

  strcpy (buffer, "      ");
  zstring = buffer;
  if (verbose)
    cerr << "Testing strip_leading_blanks on '" << buffer << "'\n";
  strip_leading_blanks (buffer);
  should_match (buffer, "", "strip_leading_blanks");

  strip_leading_blanks (zstring);
  should_match (zstring, "", "strip_leading_blanks_string");
  strcpy (buffer, "1 23   4   ");
  zstring = buffer;

  if (verbose)
    cout << "Testing strip leading blanks on '" << buffer << "'\n";
  strip_leading_blanks (buffer);
  should_match (buffer, "1 23   4   ", "strip_leading_blanks");

  strip_leading_blanks (zstring);
  should_match (zstring, "1 23   4   ", "strip_leading_blanks");

  strcpy (buffer, "      ");
  zstring = buffer;
  if (verbose)
    cerr << "Testing strip_trailing_blanks on '" << buffer << "'\n";
  strip_trailing_blanks (buffer);
  should_match (buffer, "", "strip_trailing_blanks");

  strip_trailing_blanks (zstring);
  should_match (zstring, "", "strip_trailing_blanks_string");
  strcpy (buffer, "  abcdef      ");
  zstring = buffer;
  if (verbose)
    cerr << "Test strip_trailing_blanks on '" << buffer << "'\n";
  strip_trailing_blanks (buffer);
  should_match (buffer, "  abcdef", "strip_trailing_blanks");

  strip_trailing_blanks (zstring);
  should_match (zstring, "  abcdef", "strip_trailing_blanks_string");

  strcpy (buffer, "a  bdef");
  zstring = buffer;
  if (verbose)
    cerr << "Test strip_trailing_blanks on '" << buffer << "\n";
  strip_trailing_blanks (buffer);
  should_match (buffer, "a  bdef", "strip_trailing_blanks");

  strip_leading_blanks (zstring);
  should_match (zstring, "a  bdef", "strip_trailing_blanks_string");

  strcpy (buffer, "foo");
  if (verbose)
    cout << "Testing basename on '" << buffer << "'\n";
  char *c = basename (buffer);
  should_match (c, "foo", "basename");

  strcpy (buffer, "foo/bar");
  if (verbose)
    cout << "Testing basename on '" << buffer << "'\n";
  c = basename (buffer);
  should_match (c, "bar", "basename");

  strcpy (buffer, "/usr/local/");
  if (verbose)
    cout << "testing basename on '" << buffer << "'\n";
  c = basename (buffer);
  should_match (c, "", "basename");

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

  string sbuffer = " ABCdEF123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?";
  if (verbose)
    cout << "Testing lowercase on '" << sbuffer << "'\n";
  to_lowercase (sbuffer);
  should_match (sbuffer, " abcdef123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?", "to_lowercase");

  sbuffer = " abcDef123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?";
  if (verbose)
    cout << "Testing uppercase on '" << sbuffer << "'\n";
  to_uppercase (sbuffer);
  should_match (sbuffer, " ABCDEF123~!@#$%^&*()_-+={}[]\\|'\";:,./<>?", "to_uppercase");

  sbuffer = "abc";
  if (verbose)
    cerr << "testing uppercase on '" << sbuffer << "'\n";
  to_uppercase (sbuffer);
  should_match (sbuffer, "ABC", "to_uppercase");
  if (verbose)
    cout << "Testing remove blanks\n";
  test_remove_blanks ();

  sbuffer = "abc";
  if (verbose)
    cerr << "testing compress_multiple blanks on '" << sbuffer << "'\n";
  int tmp = compress_blanks (sbuffer);
  numeric_match (tmp, 0, sbuffer.data (), "compress_blanks");
  should_match (sbuffer, "abc", "compress_blanks");

  sbuffer = " abc ";
  if (verbose)
    cerr << "testing compress_multiple blanks on '" << sbuffer << "'\n";
  tmp = compress_blanks (sbuffer);
  numeric_match (tmp, 0, sbuffer.data (), "compress_blanks");
  should_match (sbuffer, " abc ", "compress_blanks");

  sbuffer = " a b d e ";
  if (verbose)
    cerr << "testing compress_multiple blanks on '" << sbuffer << "'\n";
  tmp = compress_blanks (sbuffer);
  numeric_match (tmp, 0, sbuffer.data (), "compress_blanks");
  should_match (sbuffer, " a b d e ", "compress_blanks");

  sbuffer = "  abd    ee   ";
  if (verbose)
    cerr << "testing compress_multiple blanks on '" << sbuffer << "'\n";
  tmp = compress_blanks (sbuffer);
  numeric_match (tmp, 6, sbuffer.data (), "compress_blanks");
  should_match (sbuffer, " abd ee ", "compress_blanks");

  sbuffer = "    ";
  if (verbose)
    cerr << "testing compress_multiple blanks on '" << sbuffer << "'\n";
  tmp = compress_blanks (sbuffer);
  numeric_match (tmp, 3, sbuffer.data (), "compress_blanks");
  should_match (sbuffer, " ", "compress_blanks");

  if (verbose)
    cout << "Testing ccount\n";
  test_ccount ();

  if (verbose)
    cout << "Testing is double\n";
  test_is_double ();

  if (verbose)
    cout << "Testing tokenise\n";
  test_tokenise ();

  if (verbose)
    cout << "Testing words\n";
  test_words ();

  if (verbose)
    cout << "Testing new suffix\n";
  test_change_suffix ();

  if (verbose)
    cout << "Testing wordindex\n";
  test_wordindex ();

  if (verbose)
    cout << "Testing get_word\n";
  test_get_word ();

  if (verbose)
    cout << "Testing expand_environment\n";
  test_expand_env ();

  cerr << "All tests successful\n";

  if (verbose)
    malloc_status (stdout);
  else
    check_all_malloced (stdout);

  return 0;
}
