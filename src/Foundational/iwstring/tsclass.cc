// Tester for the string class

// Set if you want to test std::string related things

#define DO_TEST_STD_STRING 1

#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <limits>
#include <math.h>
#include <time.h>
#include <random>
#include <vector>

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
test_read_from_cin()
{
  IWString foo;

  cout << "Enter string for interactive test : ";
  int i = foo.getline (std::cin);

  cerr << "Read " << i << " characters, '" << foo << "'\n";

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
tsclass (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "ivR:");

  verbose = cl.option_count ('v');

  test_expand_env();

  test_expand_env_mf();

  test_numeric_value();

#ifdef I_CARE_ABOUT_STRSTREAM
  test_strstream_stuff();
#endif

  test_string_relationals();

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

  int rc = tsclass (argc, argv);

  return rc;
}
