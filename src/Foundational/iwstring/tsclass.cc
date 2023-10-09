// Tester for the string class

// Set if you want to test std::string related things

#define DO_TEST_STD_STRING 0

#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <limits>
#include <math.h>
#include <time.h>
#include <random>

#include <string>

#include <fstream>
#include <sstream>

#include "Foundational/cmdline/cmdline.h"

#include "iwstring.h"

using std::cerr;
using std::cout;
using std::endl;

const char * prog_name = nullptr;

int verbose = 0;     // note, is externally visible

#include "should_match.h"

void
numeric_value(int int_result, int expected_int_result,
              const IWString & result, const char *expected_result,
              const char *invoker)
{
  int die = 0;
  if (int_result != expected_int_result)
  {
    cerr << "Numeric mis-match, got " << int_result << " expected " << 
            expected_int_result << endl;
     die = 1;
  }

  if (result != expected_result)
  {
    cerr << "Result mismatch, got '" << result << "', expected '" <<
            expected_result << "'\n";
    die = 1;
  }

  if (!  die)
    return;

  cerr << "Fatal error for " << invoker << endl;
  exit(2);
}

static void
numeric_value (int int_result, int expected_int_result,
               const IWString & string_input,
               const char * invoker)
{
  if (verbose > 1)
    cerr << "Doing numeric test, got " << int_result << " expected " <<
             expected_int_result << " derived by " << invoker <<
             " from '" << string_input << "'\n";

  if (int_result == expected_int_result)
    return;

  cerr << "Numeric mis-match, got " << int_result << " expected " << 
          expected_int_result << endl;

  cerr << "Fatal error for " << invoker << " on input '" << string_input << "'\n";
  exit(2);
}

static IWString
produce_a_string (const char * s1,
                  const char * s2)
{
  IWString rc(s1);
  rc << s2;

  return rc;
}

static int
test_unhtml()
{
  IWString foo = "abcdef";
  int t = foo.unhtml();
  numeric_value(t, 0, foo, "unhtml");

  foo = "&abcdef";
  t = foo.unhtml();
  numeric_value(t, 0, foo, "unhtml");

  foo = "  &foo;bar";
  t = foo.unhtml();
  numeric_value(t, 0, foo, "unhtml");

  foo = "&nbsp;hello";
  t = foo.unhtml();
  should_match(foo, " hello", "unhtml");

  foo = "hello&nbsp;world";
  t = foo.unhtml();
  should_match(foo, "hello world", "unhtml");

  foo = "hello world&nbsp;";
  foo.unhtml();
  should_match(foo, "hello world ", "unhtml");

  foo = "&nbsp;&lt;&gt;&quot;&amp;&apos;";
  foo.unhtml();
  should_match(foo, " <>\"&'", "unhtml");

  return 1;
}

static int
test_append_hex()
{
  unsigned char y[256];
  for (int i = 0; i < 256; ++i)
  {
    y[i] = i;
  }

  IWString s;
  s.append_hex(y, 256);

  cout << "test_append_hex, all numbers to 255\n";
  cout << s << endl;

  return 1;
}

static int
test_returned_string()
{
  IWString foo = produce_a_string("abc", "def");

  should_match(foo, "abcdef", "returned string");

  foo << "1";

  should_match(foo, "abcdef1", "returned appended");

  return 1;
}

static int
test_operators_with_scalars()
{
  return 1;
}

static int
test_strncasecmp()
{
  IWString foo("hello");
  IWString bar("hello");
  
  if (0 != foo.strcmp(bar))
  {
    cerr << "foo.strcmp failed, compare '" << foo << "' with '" << bar << "'\n";
    exit(1);
  }

  if (0 != bar.strncmp(foo, 4))
  {
    cerr << "foo.strncmp failed, compare '" << foo << "' with '" << bar << "', n = 4\n";
    exit(1);
  }

  bar = "HELLO";

  if (0 != foo.strcasecmp(bar))
  {
    cerr << "foo.strcasecmp failed, compare '" << foo << "' with '" << bar << "'\n";
    cerr << "Value " << foo.strcasecmp(bar) << endl;
    exit(1);
  }

  if (0 != bar.strcasecmp(foo))
  {
    cerr << "foo.strcasecmp failed, compare '" << bar << "' with '" << foo << "'\n";
    exit(1);
  }

  return 1;
}

static int
test_index()
{
  IWString foo("hello");

  auto i = foo.rindex(' ');

  numeric_value(i, -1, foo, "rindex");

  foo = " hello";

  i = foo.rindex(' ');

  numeric_value(i, 0, foo, "rindex");

  const_IWSubstring xx(foo);
  i = xx.rindex(' ');

  numeric_value(i, 0, xx, "rindex");

  foo = "hello";
  i = foo.rindex('l');

  numeric_value(i, 3, foo, "rindex");

  i = foo.index('l');

  numeric_value(i, 2, foo, "index");

  i = foo.index('h');

  numeric_value(i, 0, foo, "index");

  return 1;
}

static int
test_strncat()
{
  IWString foo("hello");

  foo.strncat(" world", 6);

  should_match(foo, "hello world", "strncat");

  return 1;
}

static int
test_remove_suffix()
{
  IWString foo("hello world");

  if (foo.remove_suffix())
  {
    cerr << "Should not be able to remove suffix from 'hello world', result '" << foo << "'\n";
    return 0;
  }

  foo = "hello.world";

  if (! foo.remove_suffix())
  {
    cerr << "Could not remove suffix from '" << foo << "'\n";
    return 0;
  }

  should_match(foo, "hello", "remove_suffix");

  foo = ".smi";

  if (! foo.remove_suffix())
  {
    cerr << "Could not remove suffix from '" << foo << "'\n";
    return 0;
  }

  should_match(foo, "", "remove suffix .smi");

  return 1;
}

static int
test_split_into_directive_and_value()
{
  IWString foo("foo=3");

  const_IWSubstring directive;
  int v;

  if (! foo.split_into_directive_and_value(directive, '=', v) || 3 != v)
  {
    cerr << "split_into_directive_and_value failed '" << foo << "'\n";
    return 0;
  }

  return 1;
}

template <typename T, typename V>
int
test_string_relationals_template (T & s1, V & s2)
{
  s1 = "abc";

  int c = s1 < s1;
  numeric_value(0, c, s1, "less than itself");

  c = s1 > s1;
  numeric_value(0, c, s1, "greater than itself");

  s2 = "abcd";

  c = s1 < s2;
  numeric_value(1, c, s1, "less than longer string");

  c = s2 < s1;
  numeric_value(0, c, s1, "less than shorter string");

  s2 = "b";

  c = s1 < s2;
  numeric_value(1, c, s1, "less than 'b' string");

  s2 = "bbbbb";
  c = s1 < s2;
  numeric_value(1, c, s1, "less than 'bbbbb' string");

//cerr << "Finished string relationals, sizes " << sizeof(T) << " and " << sizeof(V) << endl;

  return 1;
}

//#ifdef __GNUG__
template int test_string_relationals_template(IWString &, IWString &);
template int test_string_relationals_template(const_IWSubstring &, IWString &);
template int test_string_relationals_template(IWString &, const_IWSubstring &);
template int test_string_relationals_template(const_IWSubstring &, const_IWSubstring &);
//#endif

static int
test_compare_without_case()
{
  IWString foo("ABCDEF");

  if (! foo.matches_ignore_case('A'))
  {
    cerr << "matches_ignore_case:cannot exact match single char\n";
    exit(3);
  }

  if (! foo.matches_ignore_case('a'))
  {
    cerr << "matches_ignore_case:cannot match single char\n";
    exit(1);
  }

  if (! foo.matches_ignore_case("ABCDEF"))
  {
    cerr << "matches_ignore_case:cannot exact match string\n";
    exit(1);
  }

  if (! foo.matches_ignore_case("Abcdef"))
  {
    cerr << "matches_ignore_case:cannot match string\n";
    exit(1);
  }

  return 1;
}

static int
test_string_relationals()
{
  IWString s1, s2;
  test_string_relationals_template(s1, s2);

  const_IWSubstring b1, b2;
  test_string_relationals_template(b1, b2);

  test_string_relationals_template(s1, b1);

  test_string_relationals_template(b1, s1);

  return 1;
}

static int
test_change()
{
  IWString foo("hello world");

  foo.change(0, 4, "hello");

  should_match(foo, "hello world", "change 0-4");

  foo.change(6, 10, "world");

  should_match(foo, "hello world", "change 6-10");

  foo.change(3, 3, "q");

  should_match(foo, "helqo world", "change 3");

  foo.change(2, 3, "ll");

  should_match(foo, "hello world", "change 6-10");

  foo.change(6, 10, "potato");

  should_match(foo, "hello potato", "change 6-10");

  foo.change(6, 11, "");

  should_match(foo, "hello ", "change 6-11 to nothing");

  foo.change(5, 5, " world");

  should_match(foo, "hello world", "change 5 back to world");

  foo.change(9, 10, "ld torana");

  should_match(foo, "hello world torana", "change at end - grow");

  foo = "abcdefghijklmnopqrstuvwxyz0123456789";

  foo.change(5, 10, "x");

  should_match(foo, "abcdexlmnopqrstuvwxyz0123456789", "change");

  return 1;
}

static int
test_read_lines (const char * fname)
{
  std::ifstream input(fname, std::ios::in);

  if (! input.good())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 3;
  }

  IWString buffer;

  while (input >> buffer)
  {
    cout << "Read '" << buffer << "'\n";
  }

  return 1;
}

static int
test_basename()
{
  IWString foo("foo");

  const_IWSubstring b;

  foo.iwbasename(b);

  should_match(b, "foo", "basename - nothing");

  foo = "foo/";

  foo.iwbasename(b);

  should_match(b, "foo", "basename - strip 1");

  foo = "glerf//////";

  foo.iwbasename(b);

  should_match(b, "glerf", "basename -strip lots");

  foo = "   ";

  foo.iwbasename(b);

  should_match(b, "   ", "basename - nothing");

  foo = "abc/def";

  foo.iwbasename(b);

  should_match(b, "def", "basename 1");

  const_IWSubstring bar("a/b/c/d/e/f");

  bar.iwbasename(b);

  should_match(b, "f", "basename (substring)");

  return 1;
}

#if (DO_TEST_STD_STRING)

static int
test_std_string()
{
  std::string foo = "hello world";

  IWString s1(foo);

  should_match(s1, "hello world", "IWString constructor from std::string");

  const_IWSubstring s2(foo);

  should_match(s2, "hello world", "const_IWSubstring constructor from std::string");

  foo = " hunga dunga";

  s1 = foo;

  should_match(s1, " hunga dunga", "IWString operator = const std::string");

  s2 = foo;

  should_match(s2, " hunga dunga", "const_IWSubstring operator = const std::string");

  foo = " extra";

  s1 += foo;

  should_match(s1, " hunga dunga extra", "IWString operator += std::string");

  s1 = "hello world";

  foo = s1;

  if (foo != "hello world")
  {
    cerr << "Assignment to std::string from IWString failed '" << foo << "' should be 'hello world'\n";
    return 0;
  }

  s2 = " &#^@ _";

/*foo = s2;

  if (foo != " &#^@ _")
  {
    cerr << "Assignment to std::string from const_IWSubstring failed '" << foo << "' should be ' &#^@ _'\n";
    return 0;
  }*/

  return 1;
}

#endif

static int
test_remove_to_first()
{
  IWString foo("hello world");

  int chars_removed = foo.remove_up_to_first('e');
  numeric_value(chars_removed, 2, foo, "IWString.remove_up_to_first");

  should_match(foo, "llo world", "IWString.remove_up_to_first");

  const_IWSubstring bar(foo);

  chars_removed = bar.remove_up_to_first('z');
  numeric_value(chars_removed, 0, bar, "const_IWSubstring.remove_up_to_first not found");

  chars_removed = bar.remove_up_to_first(' ');
  numeric_value(chars_removed, 4, bar, "const_IWSubstring.remove_up_to_first");
  should_match(bar, "world", "const_IWSubstring.remove_up_to_first");

  return 1;
}


static int
test_perl_split()
{
  IWString foo = "abc d e fghijkl m n OPQRST UVW XYZ";
  resizable_array_p<const_IWSubstring> tokens;

  int ntokens = foo.split(tokens);

  assert (ntokens == tokens.number_elements());

  should_match(*(tokens[0]), "abc", "perl split, first token");
  should_match(*(tokens[1]), "d", "perl split, 2nd token");
  should_match(*(tokens[ntokens - 1]), "XYZ", "perl split, last token");

  foo = "abc/def/ ";
  ntokens = foo.split(tokens, '/');

  assert (ntokens == tokens.number_elements());
  numeric_value(ntokens, 3, foo, "split on /");

  should_match(*(tokens[0]), "abc", "perl split, first token slash");
  should_match(*(tokens[1]), "def", "perl split, 2nd token slash");
  should_match(*(tokens[2]), " ", "perl split, last token slash");

// Since the iwaray<> template cannot be resized, we do those tests within braces

  if (1)
  {
    iwaray<const_IWSubstring> t2;

    foo = "hello world   testing";
  
    ntokens = foo.split(t2, ' ');

    numeric_value(ntokens, 3, foo, "split on ' '");

    should_match(t2[0], "hello", "iwaray split first token");
    should_match(t2[1], "world", "iwaray split first token");
    should_match(t2[2], "testing", "iwaray split first token");
  }

  if (1)
  {
    iwaray<const_IWSubstring> t2;

    foo = "w0;1;hello world;3;444;glerf;last one";

    ntokens = foo.split(t2, ';');

    numeric_value(ntokens, 7, foo, "iwaray split on ';'");

    should_match(t2[0], "w0", "iwaray split first token ;");
    should_match(t2[1], "1", "iwaray split 2nd token ;");
    should_match(t2[2], "hello world", "iwaray split 3rd token ;");
    should_match(t2[3], "3", "iwaray split 4th token ;");
    should_match(t2[4], "444", "iwaray split 5th token ;");
    should_match(t2[5], "glerf", "iwaray split 6th token ;");
    should_match(t2[6], "last one", "iwaray split 7th token ;");
  }

  if (1)
  {
    foo = " hello     world     ";
    iwaray<IWString> tokens;

    int ntokens = foo.split(tokens, ' ');

    numeric_value(ntokens, 2, foo, "iwaray split on ' '");
    should_match(tokens[0], "hello", "iwaray first token");
    should_match(tokens[1], "world", "iwaray second token");
  }

  return 1;
}

static int
test_operators()
{
  IWString foo;

  foo << "Hello " << 3 << '5';

  should_match(foo, "Hello 35", "operator << IWString");

  foo = "Hello ";
  foo << 8 << " world";

  should_match(foo, "Hello 8 world", "operator << IWSring");

  IWString bar = foo;

  if (foo == bar)
    ;
  else
  {
    cerr << "Operator == failed on two IWString objects\n";
  }

  if (foo != bar)
  {
    cerr << "Operator != failed on two identifical IWString objects\n";
  }

  const_IWSubstring g = "hello world";
  const_IWSubstring h = g;

  if (g == h)
    ;
  else
  {
    cerr << "const_IWSubstring::operator== (const const_IWSubstring &) failed\n";
    return 0;
  }

  if (g != h)
  {
    cerr << "const_IWSubstring::operator != (const const_IWSubstring &) failed\n";
    return 0;
  }

#ifdef DOES_NOT_WORK_MEMORY_PROBLEM
  g = 'q';

  if (g == 'q')
    ;
  else
  {
    cerr << "const_IWSubstring::operator != (char) failed, line " << __LINE__ << "\n";
    return 0;
  }
#endif

  const char * rhs = nullptr;
  IWString nn = rhs;

  if (0 != nn.length())
  {
    cerr << "Null assignment failed\n";
    return 0;
  }

  return 1;
}

static int
test_split()
{
  IWString foo = "hello/world";
  const_IWSubstring before, after;

  int rc = foo.split(before, ' ', after);
  numeric_value(rc, 0, foo, "split, separator not present");

  rc = foo.split(before, '/', after);
  numeric_value(rc, 2, foo, "split");

  should_match(before, "hello", "split before");
  should_match(after, "world", "split after");

  foo = "/blerfl";

  rc = foo.split(before, '/', after);
  numeric_value(rc, 1, foo, "split, separator at beginning");
  should_match(before, "", "split, nothing before separator");
  should_match(after, "blerfl", "split, after separator at beginning");

  foo = "gbwerdy?";
  rc = foo.split(before, '?', after);
  numeric_value(rc, 1, foo, "split, separator at end");
  should_match(before, "gbwerdy", "split, before separator at end");
  should_match(after, "", "split, after token at end");

  return 1;
}

/*
  Too hard to make this work, just skip
*/

#ifdef I_CARE_ABOUT_STRSTREAM

#if (__GNUC_MINOR__ == 95)
static int
test_strstream_stuff()
{

  IWString foo = "hello ";

  ostrstream ss;
  ss << 35;

  foo << ss;

  should_match(foo, "hello 35", "operator << strstream");


  ostrstream ss2;
  ss2 << foo;

  foo = ss2;

  should_match(foo, "hello 35", "operator << from assignment");

  return 1;
}

#else

static int
test_strstream_stuff()
{
  IWString foo = "hello ";

  ostringstream ss;
  ss << 35;

  foo << ss;

  should_match(foo, "hello 35", "operator << strstream");


  ostringstream ss2;
  ss2 << foo;

  foo = ss2;

  should_match(foo, "hello 35", "operator << from assignment");

  return 1;
}

#endif

#endif

static int
test_operator_plus()
{
  IWString foo("abc");
  foo += "def";

  should_match(foo, "abcdef", "operator +=");

  IWString bar = foo + "1";

  should_match(bar, "abcdef1", "operator +");

  foo = "abc";
  bar = "111";

  IWString s = foo + "3" + bar;

  should_match(s, "abc3111", "operator + multiple");

  foo='3';
  s = "12" + foo;

  should_match(s, "123", "operator + string");

  return 1;
}

static int
test_append_number_width()
{
  IWString foo;
  append_number(foo, 1, 1);

  should_match(foo, "1", "append_number 11");

  append_number(foo, 1, 2);

  should_match(foo, "1 1", "append_number 1 2");

  append_number(foo, 876, 1);

  should_match(foo, "1 1876", "append_number 1 *");

  append_number(foo, -312, 1);

  should_match(foo, "1 1876-312", "append_number 1 -*");

  foo.resize(0);

  append_number(foo, 1, 4);

  should_match(foo, "   1", "append_number 4");

  append_number(foo, -14, 4);

  should_match(foo, "   1 -14", "append_number -4");
  return 1;
}

static int
test_append_number()
{
  IWString foo;
  foo.append_number(0);

  should_match(foo, "0", "append number 0");

  foo.append_number(-8);

  should_match(foo, "0-8", "append number -8");

  foo.append_number(123456789);

  should_match(foo, "0-8123456789", "append number 123456789");

  foo = "";

  foo.append_number(-9);

  should_match(foo, "-9", "append number -9");

  foo.append_number(1000);

  should_match(foo, "-91000", "append number 1000");

  return 1;
}

static int
test_operator_angle_angle()
{
  IWString foo ("hello");
  IWString bar ("world");

  foo << " " << bar;

  should_match (foo, "hello world", "operator <<");

  const_IWSubstring orl = bar.from_to (1, 3);

  foo << " " << orl;

  should_match (foo, "hello world orl", " operator << substring");

  return 1;
}

static int
test_operator_plus_plus()
{
  const IWString foo = "hello world";

  const_IWSubstring bar;

  foo.from_to (0, foo.length() - 1, bar);

  should_match (bar, "hello world", "substring all");

  IWString result;
  char c;
  while ((c = ++bar))
    result += c;

  should_match (result, "hello world", "substring build");

  bar = foo.from_to (0, foo.length() - 1);

  should_match (bar, "hello world", "substring assignment from from_to");

  result.iwtruncate (0);
  while ((c = bar++))
    result += c;

  should_match (result, "ello world", "substring build via postfix ++");
  
  return 1;
}

static int
test_numeric_value()
{
  const_IWSubstring foo = "hello world";
  int i;
  float x;
  double xx;
  if (foo.numeric_value (i) || foo.numeric_value (x) || foo.numeric_value (xx))
  {
    cerr << "Should not make numeric match '" << foo << "'\n";
    exit (2);
  }

  foo = "1234";
  if (! foo.numeric_value (i) || i != 1234)
  {
    cerr << "Not getting correct numeric value '" << foo << "'\n";
    exit (2);
  }

  if (! foo.numeric_value (x) || fabs (x - 1234) > 0.1)
  {
    cerr << "Incorrect float translation '" << foo << "' result " << x << endl;
    exit (3);
  }

  if (! foo.numeric_value (xx) || fabs (xx - 1234) > 0.1)
  {
    cerr << "Incorrect float translation '" << foo << "' result " << xx << endl;
    exit (4);
  }

  foo = "-1.234";
  double expected_value = -1.234;

  if (foo.numeric_value (i))
  {
    cerr << "Incorrect match as int '" << foo << "'\n";
    exit (5);
  }

  if (! foo.numeric_value (x) || fabs (x - expected_value) > 0.1)
  {
    cerr << "Incorrect float value '" << foo << "' result " << x << endl;
    exit (6);
  }

  if (! foo.numeric_value (xx) || fabs (xx - expected_value) > 0.1)
  {
    cerr << "Incorrect double value '" << foo << "' result " << xx << endl;
    exit (7);
  }

  foo = "1234567890";
  if (! foo.numeric_value (i))
  {
    cerr << "Valid 10 digit number not recognised '" << foo << "'\n";
    exit (8);
  }

  numeric_value (i, 1234567890, foo, "numeric value - 10 digit valid");

  foo = "2147483647";
  if (! foo.numeric_value (i))
  {
    cerr << "Int max not recognised as int '" << foo << "'\n";
    exit (9);
  }

  numeric_value (i, std::numeric_limits<int>::max(), foo, "numeric value, INT_MAX");

  foo = "2147483648";

  if (foo.numeric_value (i))
  {
    cerr << "Number larger than INT_MAX recognised '" << foo << "'\n";
    exit (10);
  }

  foo = "3147483647";

  if (foo.numeric_value (i))
  {
    cerr << "3 billion number recognised as int '" << foo << "'\n";
    exit (11);
  }

  foo = "9223372036854775807";
  long int i8;
  if (! foo.numeric_value(i8) || std::numeric_limits<long int>::max() != i8)
  {
    cerr << "Largest int '" << foo << "' not recognised as " << std::numeric_limits<long int>::max() << endl;
    exit(12);
  }

  foo = "-9223372036854775808";
  if (! foo.numeric_value(i8) || std::numeric_limits<long int>::min() != i8)
  {
    cerr << "Smallest int '" << foo << "' not recognised as " << std::numeric_limits<long int>::min() << endl;
    exit(12);
  }

  foo = "9223372036858775807";     // larger than max()
  if (foo.numeric_value(i8))
  {
    cerr << "Long int larger than max incorrectly recognised '" << foo << "'\n";
    exit(12);
  }

  foo = "-9223372036854775809";     // smaller than min()
  if (foo.numeric_value(i8))
  {
    cerr << "Long int smaller than min incorrectly recognised '" << foo << "'\n";
    exit(12);
  }


  foo = "18446744073709551615";
  unsigned long int ul8;
  if (! foo.numeric_value(ul8) || std::numeric_limits<unsigned long int>::max() != ul8)
  {
    cerr << "Largest int '" << foo << "' not recognised as " << std::numeric_limits<unsigned long int>::max() << endl;
    exit(12);
  }

  foo = "18446744073709551616";    // larger than max
  if (foo.numeric_value(ul8))
  {
    cerr << "String '" << foo << "' larger than " << std::numeric_limits<unsigned long int>::max() << " recognised\n";
    exit(12);
  }

  foo = "2147483647";    // int max
  if (! foo.numeric_value(i) || i != std::numeric_limits<int>::max())
  {
    cerr << "Largest int '" << foo << "' not recognised as " << std::numeric_limits<int>::max() << endl;
    exit(12);
  }

  foo = "2147483648";    // larger than int max
  if (foo.numeric_value(i))
  {
    cerr << "Larger than largest int '" << foo << "' recognised\n";
    exit(12);
  }

  foo = "-2147483648";    // min int
  if (! foo.numeric_value(i) || i != std::numeric_limits<int>::min())
  {
    cerr << "Smallest int '" << foo << "' not recognised as " << std::numeric_limits<int>::min() << endl;
    exit(12);
  }

  foo = "-2147483658";    // smaller than min int
  if (foo.numeric_value(i))
  {
    cerr << "Value smaller than smallest int '" << foo << "' recognised, value " << i << endl;
    exit(12);
  }

  foo = "4020100000";
  unsigned long long ull;

  if (! foo.numeric_value(ull) || 4020100000 != ull)
  {
    cerr << "unsigned long long " << foo << " not recognised\n";
    exit(12);
  }


  foo = "666666.666666667";

  xx = 0.0;

  if (! foo.numeric_value (xx) || fabs (xx - 666666.666666667) > 0.1)
  {
    cerr << "Could not convert '" << foo << "' to double, result " << xx << endl;
    exit (14);
  }

//cerr.precision (10);
//cerr << "Foo '" << foo << "' result " << xx << endl;

  foo = "4245319541";
  unsigned int ui;
  if (! foo.numeric_value (ui) || ui != 4245319541U)
  {
    cerr << "Cound not convert large unsigned int '" << foo << "'\n";
    exit (8);
  }

  std::default_random_engine eng(std::random_device{}());
  std::uniform_int_distribution<> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
  for (auto i = 0; i < 10000; ++i)
  {
    int j = dist(eng);

    IWString s;
    s << j;

    int k;
    if (! s.numeric_value(k) || k != j)
    {
      cerr << "Integer conversion mismatch initial " << j << " string '" << s << "' from string " << k << endl;
      exit(10);
    }
  }

  srand (time(NULL));

  int numeric_value_differences = 0;    // very bad
  int numerically_equivalent = 0;       // close enough

  std::uniform_real_distribution<double> float_dist(std::numeric_limits<double>::min(), std::numeric_limits<double>::max());

  for (int i = 0; i < 100; i++)
  {
    double t1 = float_dist(eng);

    char buffer[32];

    sprintf (buffer, "%.14g", t1);

//  cerr << "Number is '" << buffer << "'\n";

    IWString s = buffer;

    double t2;
    if (! s.numeric_value (t2))
    {
      cerr << "Yipes, cannot read back numeric value '" << s << "'\n";
      exit (4);
    }

    if (t1 == t2)
      continue;

    double largest;
    if (0.0 == t1)
      largest = fabs (t2);
    else if (0.0 == t2)
      largest = fabs (t1);
    else if (fabs (t1) > fabs (t2))
      largest = fabs (t2);
    else
      largest = fabs (t1);

    if (fabs (t1 - t2) / largest < 1.0e-13)
    {
      numerically_equivalent++;
      continue;
    }

    cerr << "Numeric value difference " << t1 << " vs " << t2 << " diff " << (t1 - t2) << endl;
    numeric_value_differences++;
  }

  if (1)
  {
    std::vector<double> v;
    v.push_back(std::numeric_limits<double>::min());
    v.push_back(std::numeric_limits<double>::max());
    v.push_back(std::numeric_limits<double>::epsilon());
    v.push_back(std::numeric_limits<double>::denorm_min());

    char buffer[32];
    for (unsigned int i = 0; i < v.size(); ++i)
    {
      sprintf(buffer, "%.14g", v[i]);

      IWString s = buffer;

      double t2;

      if (! s.numeric_value(t2))
      {
        cerr << "Yipes, cannot reinterpret " << v[i] << " as double, string rep '" << s << "'\n";
        exit(1);
      }

      if (t2 == v[i])
        continue;

      if (fabs(t2 - v[i]) / std::max(fabs(t2), fabs(v[i])) < 1.0e-11)
        continue;

      cerr << "Numeric difference " << v[i] << " string rep '" << s << "' reinterpreted as " << t2 << " diff " << fabs(t2 - v[i]) << endl;
    }
  }

  if (numeric_value_differences)
  {
    cerr << numeric_value_differences << " values were numerically different\n";
    exit (5);
  }

  if (numerically_equivalent)
    cerr << numerically_equivalent << " tests were numerically equivalent\n";

  foo = "1.0e+03";

  double y;
  if (! foo.numeric_value (y))
  {
    cerr << "Cannot determine numeric from '" << foo << "'\n";
    exit (3);
  }

  if (fabs (y - 1.0e+03) > 0.0003)
  {
    cerr << "Incorrect numeric translation from '" << foo << "' got " << y << endl;
  }

  foo = "-2.0E+03";
  y = 0.0;
  if (! foo.numeric_value (y) || fabs (y + 2.0e+03) > 0.003)
  {
    cerr << "INcorrect numeric translation from '" << foo << "', bot " << y << endl;
    exit (6);
  }

  foo = ".5e-3";
  y = 0.0;
  if (! foo.numeric_value (y) || fabs (y - 0.5e-03) > 0.0003)
  {
    cerr << "INcorrect numeric translation from '" << foo << "', bot " << y << endl;
    exit (6);
  }

  foo = "100E-02";
  if (! foo.numeric_value (y) || fabs (y - 1.0) > 0.0003)
  {
    cerr << "INcorrect numeric translation from '" << foo << "', bot " << y << endl;
    exit (6);
  }

  foo = "00003.4";
  y = 0.0;
  if (! foo.numeric_value(y) || fabs(y - 3.4) > 0.00001)

  {
    cerr << "Incorrect numeric translation - leading 0's '" << foo << "' got " << y << endl;
    exit (3);
  }

  foo = "-00000.4";
  y = 0.0;
  if (! foo.numeric_value(y) || fabs(y + 0.4) > 0.00001)

  {
    cerr << "Incorrect numeric translation - leading 0's '" << foo << "' got " << y << endl;
    exit (3);
  }

  return 1;
}

static int
test_find()
{
  IWString foo = "hello world";
  int i = foo.find ("hello");
  numeric_value (i, 0, foo, "find hello");

  i = foo.find ("world");
  numeric_value (i, 6, foo, "find world");

  IWString bar = "hello";
  i = foo.find (bar);
  numeric_value (i, 0, foo, "find IWString hello");

  bar = "world";
  i = foo.find (bar);
  numeric_value (i, 6, foo, "find IWString world");

  const_IWSubstring s = foo.from_to (0, 4);
  i = foo.find (s);
  numeric_value (i, 0, foo, "find const_IWSubstring hello");

  s = foo.from_to (6, 10);
  i = foo.find (s);
  numeric_value (i, 6, foo, "find const_IWSubstring world");


  return 1;
}

static int
test_strncpy()
{
  IWString foo;

  foo.strncpy ("hello world", 11);

  should_match (foo, "hello world", "strncpy 11");

  foo.strncpy ("hello world", 5);

  should_match (foo, "hello", "strncpy 5");

  const_IWSubstring bar;
  
  bar = foo.from_to (0, 4);
  should_match (bar, "hello", "from_to 4");

  char buffer[9];

  bar.copy_to_char_array (buffer);

#ifdef NONONO
// this is invalid! bar points to a subset of foo, so we cannot assign it this way!!
   unfortunately, there is no autmatic guard against this
  foo = bar;

  should_match (foo, "hello", "from const_IWSubstring::operator=");
#endif

  return 1;
}

static int
test_prevword()
{
  IWString foo = "hello world from foo";
  int i = foo.length() - 1;
  IWString bar;
  IWString result;
  while (foo.prevword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  should_match (result, "foo from world hello", "prevword");

  foo = "  hello     world";
  i = foo.length() - 1;
  result = "";
  while (foo.prevword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  should_match (result, "world hello", "prevword leading space");

  foo = "  hello world   ";
  i = foo.length() - 1;
  result = "";
  while (foo.prevword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  should_match (result, "world hello", "prevword trailing space");

  return 1;
}

static int
test_nextword_single_delimiter()
{
  IWString foo = "hello world from foo";
  int i = 0;
  const_IWSubstring token;
  IWString tmp;

  while (foo.nextword_single_delimiter(token, i, ' '))
  {
    tmp.append_with_spacer(token);
  }

  should_match(tmp, foo, "nextword_single_delimiter single space");

  foo = "hello  world";

  i = 0;
  foo.nextword_single_delimiter(token, i, ' ');

  should_match(token, "hello", "nextword_single_delimiter first word");

  foo.nextword_single_delimiter(token, i, ' ');

  should_match(token, "", "nextword_single_delimiter 2 spaces");

  foo.nextword_single_delimiter(token, i, ' ');

  should_match(token, "world", "nextword_single_delimiter last token");

  const_IWSubstring glerf ("  abc   def  ");
  i = 0;

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "", "nextword_single_delimiter leading blank 0");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "", "nextword_single_delimiter leading blank 1");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "abc", "nextword_single_delimiter leading blank 1");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "", "nextword_single_delimiter middle blank 1");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "", "nextword_single_delimiter middle blank 2");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "def", "nextword_single_delimiter 2nd word");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "", "nextword_single_delimiter 1st trailing blank");

  glerf.nextword_single_delimiter(token, i, ' ');
  should_match(token, "", "nextword_single_delimiter 2nd trailing blank");

  glerf = "a,b,c,d,,,e,,g";

  i = 0;
  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "a", "nextword_single_delimiter comma1");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "b", "nextword_single_delimiter comma2");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "c", "nextword_single_delimiter comma3");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "d", "nextword_single_delimiter comma4");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter first empty");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter second empty");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "e", "nextword_single_delimiter after empties");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter next empty");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "g", "nextword_single_delimiter last token");

  glerf = "a,b,,,,";
  i = 0;

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "a", "nextword_single_delimiter first");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "b", "nextword_single_delimiter second");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter first empty");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter 2nd empty");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter 3nd empty");

  glerf.nextword_single_delimiter(token, i, ',');
  should_match(token, "", "nextword_single_delimiter 4th empty");

  return 1;
}

static int
test_nextword()
{
  IWString foo = "hello world from foo";
  int i = 0;
  IWString bar;
  IWString result;
  while (foo.nextword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  should_match (foo, result.null_terminated_chars(), "nextword");

  const_IWSubstring sbar;
  result.iwtruncate (0);

  i = 0;
  while (foo.nextword (sbar, i))
  {
    if (result.length())
      result += ' ';
    result += sbar;
  }

  should_match (result, foo.null_terminated_chars(), "nextword const_IWSubstring");

  foo = "  hello     world";
  i = 0;
  result = "";
  while (foo.nextword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  should_match (result, "hello world", "nextword leading space");

  foo = "  hello world   ";
  i = 0;
  result = "";
  while (foo.nextword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  should_match (result, "hello world", "nextword trailing space");

  foo = "d1";

  IWString d1;
  i = 0;
  foo.nextword (d1, i);

  should_match (d1, "d1", "nextword - one token");

  sbar = "d1";
  i = 0;
  sbar.nextword (d1, i);

  should_match (d1, "d1", "substring nextword - one token");

  foo = "hello world ";
  i = 0;
  foo.nextword (d1, i);
  foo.nextword (d1, i);

  should_match (d1, "world", "substring nextword trailing space");

  return 1;
}

static int
test_gsub()
{
  IWString foo = "hello world";
  int i = foo.gsub (' ', '_');
  numeric_value (i, 1, foo, "hello_world", "gsub space to underscore");

  i = foo.gsub ('o', ' ', 1);
  numeric_value (1, i, foo, "hell _world", "gsub 'o' to space (1)");

  i = foo.gsub ('o', ' ');
  numeric_value (1, i, foo, "hell _w rld", "gsub 'o' to space (2)");

  i = foo.gsub (' ', 'o');
  numeric_value (2, i, foo, "hello_world", "gsub ' ' to 'o'");
  
  return 1;
}

static int
test_gsub_string()
{
  IWString foo ("hello world");

  int i = foo.gsub (" ", "_");
  numeric_value (i, 1, foo, "hello_world", "gsub space string to underscore");

  foo.gsub ('_', ' ');

  i = foo.gsub ("h", "H");
  numeric_value (i, 1, foo, "Hello world", "gsub char string at start");

  i = foo.gsub ("d", "D");
  numeric_value (i, 1, foo, "Hello worlD", "gsub char string at end");

  i = foo.gsub (" ", "");
  numeric_value (i, 1, foo, "HelloworlD", "gsub char string to nothing");

  foo = "hello world";

  i = foo.gsub ("o w", "O W");
  numeric_value (i, 1, foo, "hellO World", "gsub 3 chars");

  foo = "hello world";

  i = foo.gsub ("hello", "HELLO");
  numeric_value (i, 1, foo, "HELLO world", "gsub first 5 chars");

  i = foo.gsub ("world", "WORLD");
  numeric_value (i, 1, foo, "HELLO WORLD", "gsub last 5 chars");

  foo = "hello world";

  i = foo.gsub ("hello", "gday");
  numeric_value (i, 1, foo, "gday world", "gsub 'hello' to 'gday'");

  i = foo.gsub ("world", "mate");
  numeric_value (i, 1, foo, "gday mate", "gsub 'world' to 'mate'");

  i = foo.gsub ("day mat", "e");
  numeric_value (i, 1, foo, "gee", "gsub 'day mat' to 'e'");

  foo = "hello world";
  i = foo.gsub ("hello world", "gday mate");
  numeric_value (i, 1, foo, "gday mate", "gsub 'hello world' to 'gday mate'");

  foo = "gday mate";

  i = foo.gsub ("gday", "hello");
  numeric_value (i, 1, foo, "hello mate", "gsub 'gday' to 'hello'");

  i = foo.gsub ("mate", "world");
  numeric_value (i, 1, foo, "hello world", "gsub 'mate' to 'world'");

  foo = "hello hello";

  i = foo.gsub ("hello", "gdbye");
  numeric_value (i, 2, foo, "gdbye gdbye", "gsub 'hello' to 'gdbye'");

  foo = "hello hello hello hello";

  i = foo.gsub ("hello", "xx");
  numeric_value (i, 4, foo, "xx xx xx xx", "gsub 'hello' to 'xx'");

  i = foo.gsub ("xx", "yyy");
  numeric_value (i, 4, foo, "yyy yyy yyy yyy", "gsub 'hello' to 'yyy'");

  foo = "";

  int n = 100;

  for (int i = 0; i < n; i++)
  {
    foo << " x";
  }

  i = foo.gsub ("x", "yz");
  if (3 * n != foo.length())
  {
    cerr << "gsub'd string wrong length " << foo.length() << " expected " << (3 * n) << endl;
  }

  for (int i = 0; i < 3 * n; i++)
  {
    int j = i % 3;
    if (0 == j && ' ' == foo[i])
      ;
    else if (1 == j && 'y' == foo[i])
      ;
    else if (2 == j && 'z' == foo[i])
      ;
    else
    {
      cerr << "Improper character at i = " << i << " got '" << foo[i] << "'\n";
    }
  }

  return 1;
}

static int
test_next()
{
  IWString foo = "01234012340123401234";

  int nzero = 0;

  int j = 0;
  while (foo.next ('0', j))
  {
    assert ('0' == foo[j]);
    j++;
    nzero++;
  }

  assert (4 == nzero);

  return 1;
}

static int
test_remove_chars()
{
  IWString foo = "hello world";

  foo.remove_chars (0, 1);
  should_match (foo, "ello world", "remove_chars");

  foo.remove_chars (1, 2);
  should_match (foo, "eo world", "remove_chars");

  foo.remove_chars (foo.length() - 1, 1);
  should_match (foo, "eo worl", "remove_chars");

  foo = "hello world";
  foo.remove_chars (0, foo.length());
  should_match (foo, "", "remove_chars all");

  foo = "hello world";

  foo.remove_from_to (0, 10);
  should_match (foo, "", "remove_from_to all");

  foo = "hello world";
  foo.remove_from_to (1, 1);
  should_match (foo, "hllo world", "remove_from_to 1 1");

  foo.remove_from_to (5, 9);
  should_match (foo, "hllo ", "remove_from_to 5 9");

  return 1;
}

static int
test_erase()
{
  IWString foo = "hello world";

  foo.erase (0, 1);
  should_match (foo, "llo world", "erase (0 1)");

  foo.erase (0, 0);
  should_match (foo, "lo world", "remove from to");

  foo.erase (1, 1);
  should_match (foo, "l world", "remove from to");

  foo.erase (1, 4);
  should_match (foo, "lld", "remove from to");

  foo.erase (2, 2);
  should_match (foo, "ll", "remove from to");

  foo = "abcdef";

  foo.erase (0);
  should_match (foo, "bcdef", "erase 0");

  foo.erase (1);
  should_match (foo, "bdef", "erase 1");

  return 1;
}

static int
test_truncate()
{
  IWString foo = "Hello world";

  foo.iwtruncate (foo.length());
  should_match (foo, "Hello world", "truncate no action");

  foo.iwtruncate (5);
  should_match (foo, "Hello", "truncate");

  foo << " world   ";

  foo.truncate_at_first (' ');
  should_match (foo, "Hello", "truncate at first");

  foo = "abc.def.ghi";

  foo.truncate_at_last ('.');
  should_match (foo, "abc.def", "truncate at last");

  return 1;
}

static int
test_insert()
{
  IWString foo ("helloworld");

  foo.insert (' ', 5);
  should_match (foo, "hello world", "insert char");

  foo.insert (" everyone", 5);
  should_match (foo, "hello everyone world", "insert string");

  foo = "world";
  foo.insert ("ello ", 0);
  should_match (foo, "ello world", "insert string at beginning");

  foo.insert ('h', 0);
  should_match (foo, "hello world", "insert char at beginning");

  IWString bar = "hello";
  foo = " world";

  foo.insert (bar, 0);
  should_match (foo, "hello world", "insert iwstring at beginning");

  bar = " the";
  foo.insert (bar, 5);
  should_match (foo, "hello the world", "insert iwstring somewhere");

  foo.insert (substr (bar, 1, 1), 3);
  should_match (foo, "heltlo the world", "insert substring");

  return 1;
}

static int
test_matches_at_position()
{
  IWString foo ("hello world");

  int tmp = foo.matches_at_position (0, "hello");
  if (0 == tmp)
  {
    cerr << "matches_at_position failed '" << foo << "'\n";
    exit (4);
  }

  tmp = foo.matches_at_position (1, "e");
  if (0 == tmp)
  {
    cerr << "matches_at_position single character failed '" << foo << "'\n";
    exit (4);
  }

  tmp = foo.matches_at_position (6, "world");
  if (0 == tmp)
  {
    cerr << "matches_at_position failed, end of variable '" << foo << "'\n";
    exit (6);
  }

  return 1;
}

static int
test_strncmp()
{
  IWString foo ("hello world");

  int tmp = foo.strncmp ("hello world", ::strlen ("hello world"));

  numeric_value (tmp, 0, foo, "hello world", "strncmp");

  const_IWSubstring bar = foo;

  tmp = bar.strncmp ("hello world", ::strlen ("hello world"));

  numeric_value (tmp, 0, bar, "hello world", "strncmp");

  return 1;
}

static int
test_expand_env_mf()
{
  if (0 != putenv (const_cast<char *>("FOO=barf")))
  {
    cerr << "putenv failed\n";
    exit (8);
  }

  if (0 != putenv (const_cast<char *>("FOO2=jarjar")))
  {
    cerr << "putenv failed\n";
    exit (8);
  }

  IWString expanded;

  const_IWSubstring s("hello ${FOO} and $FOO2");

  if (2 != s.expand_environment_variables(expanded))
  {
    cerr << "expand_environment_variables MF failed\n";
    return 0;
  }

  should_match(expanded, "hello barf and jarjar", "test_expand_env_mf");

  return 1;
}

static int
test_expand_env()
{
  if (0 != putenv (const_cast<char *>("FOO=barf")))
  {
    cerr << "putenv failed\n";
    exit (8);
  }

  IWString expanded;

  if (! expand_environment_variables ("$FOO/bar", expanded))
  {
    cerr << "expand_environment_variables failed\n";
    exit (9);
  }

  should_match (expanded, "barf/bar", "test_expand_env");

  return 1;
}

static int
test_looks_like()
{
  IWString foo ("he");

  if (! foo.looks_like ("hello world", 2))
  {
    cerr << "Looks like 2 failed\n";
    exit (2);
  }

  if (foo.looks_like ("h", 1))
  {
    cerr << "Looks like failed on too short string\n";
    exit (4);
  }

  if (! foo.looks_like ("he", 2))
  {
    cerr << "Looks like failed on strings of same length\n";
    exit (9);
  }

  return 1;
}

static int
test_starts_with()
{
  IWString foo;
  foo = "abcdef";

  if (! foo.starts_with ("abc"))
  {
    cerr << "Starts with (string) failed '" << foo << "'\n";
    exit (3);
  }

  IWString bar ("abcd");
  if (! foo.starts_with (bar))
  {
    cerr << "Starts with (IWString) failed '" << foo << "'\n";
    exit (3);
  }

  const_IWSubstring bb = substr (bar, 0, 2);
  if (! foo.starts_with (bb))
  {
    cerr << "Starts with (const_IWSubtring) failed '" << foo << "'\n";
    exit (3);
  }

  return 1;
}

static int
test_ends_with()
{
  IWString foo ("barbleqwerf");

  if (! foo.ends_with ('f'))
  {
    cerr << "Ends with char failed\n";
    exit (5);
  }

  if (! foo.ends_with ("f"))
  {
    cerr << "Ends with char * failed\n";
    exit (6);
  }

  if (! foo.ends_with ("qwerf"))
  {
    cerr << "Ends with char * (n) failed\n";
    exit (7);
  }

  if (! foo.ends_with ("barbleqwerf"))
  {
    cerr << "Ends with char * all failed\n";
    exit (8);
  }

  if (! foo.ends_with (foo))
  {
    cerr << "Ends with IWString all failed\n";
    exit (9);
  }

  IWString bar ("werf");
  if (! foo.ends_with (bar))
  {
    cerr << "Ends with IWString failed\n";
    exit (10);
  }

  return 1;
}

static int
test_is_int()
{
  IWString foo ("err");
  
  int tmp;
  numeric_value (foo.is_int (tmp), 0, foo, "err", "is_int");

  foo = "";
  numeric_value (foo.is_int (tmp), 0, foo, "is_int");

  foo = "      ";
  numeric_value (foo.is_int (tmp), 0, foo, "is_int");

  foo = "+1";
  numeric_value (foo.is_int (tmp), 1, foo, "is_int");
  numeric_value (tmp, 1, foo, "is_int (result)");

  foo = "234";
  numeric_value (foo.is_int (tmp), 1, foo, "is_int");
  numeric_value (tmp, 234, foo, "is_int (result)");

  foo = "-009";
  numeric_value (foo.is_int (tmp), 1, foo, "is_int");
  numeric_value (tmp, -9, foo, "is_int (result)");

  foo = "a1998";
  numeric_value (is_int (substr (foo, 1), tmp), 1, foo, "foo (substr)");
  numeric_value (tmp, 1998, foo, "foo (substr - result)");

  numeric_value (foo.substr (1).is_int (tmp), 1, foo, "substr.is_int");
  numeric_value (tmp, 1998, foo, "substr.is_int");

  foo = "00987";
  numeric_value (foo.is_int (tmp), 1, foo, "foo.is_int");
  numeric_value (tmp, 987, foo, "foo.is_int (value)");

  foo = "0x1";
  unsigned int utmp;
  numeric_value (foo.is_hex (utmp), 1, foo, "is_hex()");
  numeric_value (int (utmp), 1, foo, "is_hex() result");

  foo = "0x1234123";
  int expected_result = 0x1234123;
  numeric_value (foo.is_hex (utmp), 1, foo, "is_hex()");
  numeric_value (int (utmp), expected_result, foo, "is_hex()");

  return 1;
}

static int
test_from_to()
{
  IWString foo = "0123456789";

  should_match (from_to (foo, 0, 0), "0", "from_to");

  should_match (foo.from_to (0, 0), "0", "from_to");

  should_match (from_to (foo, 0, 1), "01", "from_to");

  should_match (foo.from_to (0, 1), "01", "from_to");

  should_match (from_to (foo, 2, 8), "2345678", "from_to");

  should_match (foo.from_to (2, 8), "2345678", "from_to");

  should_match (foo.from_to (2, "6"), "23456", "from_to");

  should_match (foo.before ('0'), "", "IWString::before");
  should_match (foo.before ('1'), "0", "IWString::before");
  should_match (foo.after  ('9'), "", "IWString::after");
  should_match (foo.after  ('3'), "456789", "IWString::after");

  const_IWSubstring bar;

  foo.from_to (1, 4, bar);
  should_match (bar, "1234", "from_to const_IWSubstring &");

  bar = "hello world";

  should_match (bar.before (' '), "hello", "const_IWSubstring::before");
  should_match (bar.after (' '),  "world", "const_IWSubstring::after");

  return 1;
}

static int
test_read_from_cin()
{
  IWString foo;

  cout << "Enter string for interactive test : ";
  int i = foo.getline (std::cin);

  cerr << "Read " << i << " characters, '" << foo << "'\n";

  return 1;
}

static int
test_remove_leading_chars()
{
  IWString foo ("abcdef");

  foo.remove_leading_chars (1);

  should_match (foo, "bcdef", "remove_leading_chars");

  foo.remove_leading_chars (2);

  should_match (foo, "def", "remove leading chars");

  foo = "hello     world   to you and you";

  const_IWSubstring bar = foo;

  bar.remove_leading_words (1);

  should_match (bar, "world   to you and you", "remove leading words");

  bar.remove_leading_words (2);

  should_match (bar, "you and you", "remove_leading_words");

  foo = "abcdef";

  foo.remove_leading_words (1);

  should_match (foo, "", "remove_leading_words all");

  foo = "abcdef";

  foo.shift (3, ' ');

  should_match (foo, "   abcdef", "shift (3)");

  foo.remove_leading_chars (4, 'q');

  should_match (foo, "bcdefqqqq", "remove_leading_chars with pad");

  foo = "a b c d e";

  foo.remove_word (0);

  should_match (foo, "b c d e", "remove word (first)");

  foo.remove_word (3);

  should_match (foo, "b c d", "removeword (last)");

  foo.remove_word (1);

  should_match (foo, "b d", "removeword (middle)");

  return 1;
}

static int
test_substring()
{
  IWString foo = "abcdef";
  if (foo.substr (1, 3) != "bcd")
  {
    cerr << "Substr failed\n";
    exit (9);
  }

  should_match (substr (foo, 0, 1), "a", "substr");

  should_match (foo.substr (0, 1), "a", "substr member");

  should_match (substr (foo, 0, 2), "ab", "substr");

  should_match (foo.substr (0, 2), "ab", "substr member");

  should_match (substr (foo, 2, 2), "cd", "substr");

  should_match (foo.substr (2, 2), "cd", "substr member");

  should_match (foo.substr (1), "bcdef", "substr member rest");

  if (foo.from_to (2, 3) != "cd")
  {
    cerr << "from_to failed\n";
    cerr << "Got '" << foo.from_to (2, 3) << "', expected 'cd'\n";
    exit (10);
  }

  should_match (substr (foo, 3), "def", "substr");

  should_match (substr (substr (foo, 2), 2), "ef", "substr");

  return 1;
}

static int
test_word()
{
  IWString foo = "abc def";
  should_match (foo.word (0), "abc", "word");

  should_match (foo.word (1), "def", "word");

  should_match (foo.word (2), "", "word");

  IWString bar;

  foo = "     abc     def    ";

  foo.word (0, bar);
  should_match (bar, "abc", "word 1");

  foo.word (1, bar);
  should_match (bar, "def", "word 2");

  foo.word (2, bar);
  should_match (bar, "", "word 3");

  foo.word(-1, bar);
  should_match(bar, "def", "word -1");

  foo.word(-2, bar);
  should_match(bar, "abc", "word -2");

  bar = foo.word(-1);
  should_match(bar, "def", "word -1");
  bar = foo.word(-2);
  should_match(bar, "abc", "word -2");

  bar = "hello ";

  bar.word (0, foo);
  should_match (foo, "hello", "word - just 1, trailing blank");

  const_IWSubstring glerf = "	yabba	dabba do";
  IWString w;
  glerf.whitespace_delimited_word (0, w);

  should_match (w, "yabba", "whitespace word 1");

  glerf.whitespace_delimited_word (1, w);
  should_match (w, "dabba", "whitespace word 2");

  glerf.whitespace_delimited_word (2, w);
  should_match (w, "do", "whitespace_delimited_word 3");

  return 1;
}

static int
test_nwords()
{
  IWString foo = "";
  int i = foo.nwords();
  numeric_value (i, 0, foo, "", "nwords");

  foo = " ";
  i = foo.nwords();
  numeric_value (i, 0, foo, " ", "nwords");

  i = foo.nwords ('"');
  numeric_value (i, 1, foo, " ", "nwords");

  foo = "ab";
  i = foo.nwords();
  numeric_value (i, 1, foo, "ab", "nwords");

  foo = " ab";
  i = foo.nwords();
  numeric_value (i, 1, foo, " ab", "nwords");

  foo = "ab ";
  i = foo.nwords();
  numeric_value (i, 1, foo, "ab ", "nwords");

  foo = " a  b  c ";
  i = foo.nwords();
  numeric_value (i, 3, foo, " a  b  c ", "nwords");
  return 1;
}

static int
test_compress_blanks()
{
  IWString foo ("abcdef");
  foo.compress_blanks();
  should_match (foo, "abcdef", "compress_blanks (none)");

  foo = "abc d e ";
  foo.compress_blanks();
  should_match (foo, "abc d e ", "compress_blanks (none w/spaces)");

  foo = " a b c";
  foo.compress_blanks();
  should_match (foo, " a b c", "compress blanks (none, w/spaces 2");

  foo = "      ";
  foo.compress_blanks();
  should_match (foo, " ", "compress_blanks, long blank");

  foo = " a   b";
  foo.compress_blanks();
  should_match (foo, " a b", "compress_blanks");

  foo = "       a            b     zyz      ";
  foo.compress_blanks();
  should_match (foo, " a b zyz ", "compress_blanks");

  return 1;
}

static int
test_getline (std::ifstream & is)
{
  IWString foo;

  while (is.good() && ! is.eof())
  {
    int nchars = foo.getline (is);
    cerr << "Read " << nchars << " nchars '" << foo << "'\n";
  }

  return 1;
}

static int
test_getline (const char * fname)
{
  std::ifstream is;
  is.open (fname, std::ios::in);    // for some reason, cannot get open translated properly
  if (! is.good())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  return test_getline (is);
}

#ifndef linux

static int
test_getline_fd (int fd)
{
  IWString buffer;
  
  while (buffer.getline (fd))
  {
    cout << buffer << endl;
  }

  return cout.good();
}

static int
test_getline_fd (const char * fname)
{
  int fd = IW_FD_OPEN(fname, O_RDONLY);    // cannot get compiler to do this for me

  if (fd < 0)
  {
    cerr << "Cannot open '" << fname << "' for getline test ";
    perror ("open");
    return 0;
  }

  return test_getline_fd (fd);
}
#endif

static int
test_equality()
{
  IWString foo;

  foo = "f";

  if (foo == 'f')
    ;
  else
  {
    cerr << "Equality test for rhs single character failed\n";
    exit (6);
  }

  if ('f' == foo)
    ;
  else
  {
    cerr << "Equality test for lhs single character failed\n";
    exit (7);
  }

  if (foo != 'f')
  {
    cerr << "Inequality test for rhs single character failed\n";
    exit (8);
  }

  if ('f' != foo)
  {
    cerr << "Inequality test for lhs single character failed\n";
    exit (9);
  }

  return 1;
}

static int
test_constructors()
{
  IWString foo;
  
  should_match (foo, "", "constructor");

  if (foo == "")
    ;
  else
  {
    cerr << "Operator == failed\n";
    return 1;
  }

  if (foo != "")
  {
    cerr << "Operator != failed on '" << foo << "'\n";
    return 2;
  }

  if (foo == " ")
  {
    cerr << "Operator == failed for blank string\n";
    return 3;
  }

  IWString bar ("hello world");
  should_match (bar, "hello world", "constructor");

  IWString bg = 'h';
  should_match (bg, "h", "single char assignment");

  const_IWSubstring s1;
  s1 = "hello world";
  should_match (s1, "hello world", "substring char assignment");

  s1 = bar;
  should_match (s1, "hello world", "substring from string");

  const char * hello = "hello";
  IWString h = "hello";

  if (hello == h)
    ;
  else
  {
    cerr << "Operator == (const char *, const IWString *) failed '" << hello << "' vs '" << h << "'\n";
    return 88;
  }

  const_IWSubstring hs = h;
  if (hello == hs)
    ;
  else
  {
    cerr << "Operator == (const char *, const const_IWSubstring *) failed '" << hello << "' vs '" << hs << "'\n";
    return 88;
  }

  foo = "abc def ghi";

  if (foo == "abc def ghi")
    ;
  else
  {
    cerr << "const_IWSubstring::operator== (const char *) failed\n";
    return 0;
  }

  if ("abc def ghi" == foo)
    ;
  else
  {
    cerr << "const_IWSubstring::operator == (const char *) failed\n";
  }

  if ("def" != foo.word (1))
  {
    should_match (foo.word (1), "def", "word and operator 1=");
  }

  if ("ghi" == foo.word (2))
  {
  }
  else
  {
    cerr << "Yipes, operator == seems to have failed\n";
  }

  hs = "blargldy glerf";

  if ("blargldy glerf" == hs)
    ;
  else
  {
    cerr << "operator == (const char *, const const_IWSubstring &) failed\n";
    return 0;
  }

  return 1;
}

static int
test_append()
{
  IWString foo;

  foo += "hello";
  should_match (foo, "hello", "+=");

  foo += ' ';

  should_match (foo, "hello ", "+= char");

  foo += "world";
  should_match (foo, "hello world", "+= char *");

  IWString bar = "hello " + foo;
  should_match (bar, "hello hello world", "+ constructur");

  foo = " hello";
  const_IWSubstring world (" world");

  foo += world;

  should_match (foo, " hello world", "+= string + substring");

  foo = "";
  foo.append_with_spacer ("hello");

  should_match (foo, "hello", "append_with_spacer, empty");

  foo.append_with_spacer ("world");

  should_match (foo, "hello world", "append_with_spacer started");

  return 1;
}

static int
test_chop()
{
  IWString foo ("a;sldkfjas;dlfkj");

  foo.chop();
  should_match (foo, "a;sldkfjas;dlfk", "chop");

  foo.chop (2);
  should_match (foo, "a;sldkfjas;dl", "chop");

  return 1;
}

static int
test_strips()
{
  IWString foo ("hello world");

  foo.strip_leading_blanks();

  should_match (foo, "hello world", "strip_leading_blanks");

  foo.strip_trailing_blanks();

  should_match (foo, "hello world", "strip_trailing_blanks");

  foo = " ";
  foo.strip_leading_blanks();
  should_match (foo, "", "strip leading blanks");

  foo = " ";
  foo.strip_trailing_blanks();
  should_match (foo, "", "strip trailing blanks");

  foo = "          hello world  ";
  foo.strip_leading_blanks();
  should_match (foo, "hello world  ", "strip_leading blanks");

  foo = "   hello   world        ";
  foo.strip_trailing_blanks();
  should_match (foo, "   hello   world", "strip trailing blanks");

  const_IWSubstring bar = "hello world";
  bar.strip_leading_blanks();
  should_match (bar, "hello world", "strip leading blanks (Substr)");

  bar.strip_trailing_blanks();
  should_match (bar, "hello world", "strip trailing blanks (Substr)");

  bar = "   ";
  bar.strip_leading_blanks();
  should_match (bar, "", "strip leading blanks (empty) (substr)");

  bar = "       ";
  bar.strip_trailing_blanks();
  should_match (bar, "", "strip leading blanks (empty) (substr)");

  bar = "    hello world   ";
  bar.strip_leading_blanks();
  should_match (bar, "hello world   ", "strip leading blanks (substr)");

  bar.strip_trailing_blanks();
  should_match (bar, "hello world", "strip leading blanks (substr)");

  return 1;
}

static int
tsclass (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "ivR:");

  verbose = cl.option_count ('v');

  test_constructors();

  test_equality();

  test_strips();

  test_chop();

  test_append();

  test_compress_blanks();

  test_nwords();

  test_word();

  test_substring();

  test_remove_leading_chars();

  test_from_to();

  test_is_int();

  test_starts_with();

  test_ends_with();

  test_looks_like();

  test_expand_env();

  test_expand_env_mf();

  test_insert();

  test_truncate();

  test_erase();

  test_remove_chars();

  test_next();

  test_gsub();

  test_gsub_string();

  test_nextword();

  test_nextword_single_delimiter();

  test_prevword();

  test_strncmp();
  test_strncpy();

  test_index();

  test_matches_at_position();

  test_find();

  test_operator_plus_plus();

  test_operator_angle_angle();

  test_append_number();


  test_append_number_width();

  test_numeric_value();

  test_operator_plus();

#ifdef I_CARE_ABOUT_STRSTREAM
  test_strstream_stuff();
#endif

  test_split();

  test_perl_split();

  test_operators();

  test_remove_to_first();

  test_basename();

#if (DO_TEST_STD_STRING)
  test_std_string();
#endif

  test_string_relationals();

  test_change();

  test_split_into_directive_and_value();

  test_remove_suffix();

  test_compare_without_case();

  test_strncat();

  test_strncasecmp();

  test_returned_string();

  test_unhtml();

  test_append_hex();

  if (cl.option_present ('i'))
    test_read_from_cin();

  if (0 == cl.number_elements())
    cerr << "No arguments specified, getline test not done\n";
  else
  {
    test_getline (cl[0]);
#ifndef linux
    test_getline_fd (cl[0]);
#endif
    test_read_lines (cl[0]);
  }

  cerr << "All tests successful\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

#ifdef USE_IWMALLOC
  iwmalloc_set_scramble_freed_memory (1);
  iwmalloc_initialise_memory_tracking (128);
#endif

  int rc = tsclass (argc, argv);

#ifdef USE_IWMALLOC
  iwmalloc_malloc_status (stderr);
#endif

  return rc;
}
