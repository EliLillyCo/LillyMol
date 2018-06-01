/*
  tester for iwdigits class
*/

#include <iostream>
#include <random>

#include "cmdline_v2.h"
#include "iwdigits.h"

const char * prog_name = NULL;

static int verbose = 0;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "tester for iwdigits class\n";
  cerr << " -N <n>         number of outer loop iterations to perform (default 1000)\n";
  cerr << " -n <n>         number of inner loop iterations to perform (default 1000)\n";
  cerr << " -s <s>         maximum value to test with iwdigits\n";
  cerr << " -S <fname>     use Fraction_as_String to write a random file\n";
  cerr << " -r <nrows>     number of rows in random file\n";
  cerr << " -c <ncols>     number of cols in random file\n";
  cerr << " -min <v>       min value in random file\n";
  cerr << " -max <v>       max value in random file\n";
  cerr << " -d <d>         number of digits of accuracy in output\\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static std::random_device rd;

class Random_File_Parameters
{
  public:
    int _nrows;
    int _ncols;
    double _min;
    double _max;
    int _ndigits;

  public:
    Random_File_Parameters();

    int debug_print(std::ostream &) const;

    void set_nrows(const int s) { _nrows = s;}
    void set_ncols(const int s) { _ncols = s;}
    void set_min(const double s) { _min = s;}
    void set_max(const double s) { _max = s;}
    void set_ndigits(const int s) { _ndigits = s;}

    int nrows() const { return _nrows;}
    int ncols() const { return _ncols;}
    double min() const { return _min;}
    double max() const { return _max;}
    int ndigits() const { return _ndigits;}
};

Random_File_Parameters::Random_File_Parameters()
{
  _nrows = 100;
  _ncols = 100;
  _min = -1.0;
  _max = 1.0;
  _ndigits = 3;

  return;
}

int
Random_File_Parameters::debug_print(std::ostream & output) const
{
  output << "Random_File_Parameters: " << _nrows << " rows, " << _ncols << " cols , btw " << _min << " and " << _max << " ndigits " << _ndigits << '\n';

  return 1;
}

static int
test_iwdigits_fraction_as_string_speed_test(const Random_File_Parameters & parms,
                         IWString_and_File_Descriptor & output)
{
  std::mt19937_64 rng(rd());

  Fraction_as_String f;
  f.set_leading_string(",");

  if (! f.initialise(parms.min(), parms.max(), parms.ndigits()))
  {
    cerr << "Cannot initialise fraction as string with parms min " << parms.min() << " max " << parms.max() << " digits " <<  parms.ndigits() << endl;
    return 0;
  }

  if (verbose)
    cerr << "Initialised with " << f.nbuckets() << " buckets\n";

  std::uniform_real_distribution<double> u (parms.min(), parms.max());

  const int nrows = parms.nrows();
  const int ncols = parms.ncols();

  for (int r = 0; r < nrows; ++r)
  {
    output << u(rng);
    for (int c = 1; c < ncols; ++c)
    {
      f.append_number(output, u(rng));

      output.write_if_buffer_holds_more_than(16384);
    }
    output << '\n';
  }

  output.flush();

  return 1;
}

static int
test_iwdigits_fraction_as_string_speed_test(const Random_File_Parameters & parms,
                         const char * fname)
{
  IWString_and_File_Descriptor output;

  if (! output.open(fname))
  {
    cerr << "test_iwdigits_fraction_as_string_speed_test:cannot open '" << fname << "'\n";
    return 0;
  }

  return test_iwdigits_fraction_as_string_speed_test(parms, output);
}

static int
test_iwdigits(int n, int s,
               const IWString & sbefore,
               const IWString & safter)
{
  if (verbose > 1)
    cerr << "Testing " << n << " tests across range " << s << endl;

  IWDigits iwdigits1;
  if (sbefore.length() > 0)
    iwdigits1.set_leading_string(sbefore);             // 1

  iwdigits1.initialise(s);                             // 2

  if (safter.length() > 0)
    iwdigits1.append_to_each_stored_string(safter);    // 3

  IWDigits iwdigits2;

  if (sbefore.length() > 0)
    iwdigits2.set_leading_string(sbefore);             // 1

  if (safter.length() > 0)
    iwdigits2.append_to_each_stored_string(safter);    // 3

  iwdigits2.initialise(s);                             // 2

  IWDigits iwdigits3;

  if (safter.length() > 0)
    iwdigits3.append_to_each_stored_string(safter);    // 3

  if (sbefore.length() > 0)
    iwdigits3.set_leading_string(sbefore);             // 1

  iwdigits3.initialise(s);                             // 2

  IWDigits iwdigits4;

  if (safter.length() > 0)
    iwdigits4.append_to_each_stored_string(safter);    // 3

  iwdigits4.initialise(s);                             // 2

  if (sbefore.length() > 0)
    iwdigits4.set_leading_string(sbefore);             // 1

  IWDigits iwdigits5;

  iwdigits5.initialise(s);                             // 2

  if (sbefore.length() > 0)
    iwdigits5.set_leading_string(sbefore);             // 1

  if (safter.length() > 0)
    iwdigits5.append_to_each_stored_string(safter);    // 3

  IWDigits iwdigits6;

  iwdigits6.initialise(s);                             // 2

  if (safter.length() > 0)
    iwdigits6.append_to_each_stored_string(safter);    // 3

  if (sbefore.length() > 0)
    iwdigits6.set_leading_string(sbefore);             // 1

  IWDigits * iwd[6];
  iwd[0] = &iwdigits1;
  iwd[1] = &iwdigits2;
  iwd[2] = &iwdigits3;
  iwd[3] = &iwdigits4;
  iwd[4] = &iwdigits5;
  iwd[5] = &iwdigits6;

  const auto niwd = 6;

  int rc = 1;

  std::mt19937_64 rng (rd());

  for (auto i = 0; i < n; ++i)
  {
    std::uniform_int_distribution<int> u (-s, s);

    int d = u(rng);

    IWString expected;
    if (sbefore.length())
      expected << sbefore;
    expected << d;
    if (safter.length())
      expected << safter;

    for (auto j = 0; j < niwd; ++j)
    {
      IWString got;
      iwd[j]->append_number(got, d);

      if (expected == got)
        continue;

      cerr << "Mismatch, d = " << d << " type " << j << " got '" << got << "' expected '" << expected << "', s = " << s << endl;
      rc = 0;
    }
  }

  return rc;
}


static int
test_iwdigits (int argc, char ** argv)
{
  Command_Line_v2 cl(argc, argv, "-v-N=ipos-n-ipos-s-ipos-S=s-r=ipos-c=ipos-min=float-max=float-d=ipos");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('S') && ! cl.option_present('N'))
  {
    cerr << "Must specify either iwdigits test (-N) or fraction as string speed test (-S)\n";
    usage(1);
  }


  if (cl.option_present('S'))
  {
    const char * fname = cl.option_value('S');

    Random_File_Parameters parms;

    if (cl.option_present('r'))
    {
      int r;
      cl.value('r', r);
      parms.set_nrows(r);
    }

    if (cl.option_present('c'))
    {
      int c;
      cl.value('c', c);
      parms.set_ncols(c);
    }

    if (cl.option_present("min"))
    {
      double v;
      cl.value("min", v);
      parms.set_min(v);
    }

    if (cl.option_present("max"))
    {
      double v;
      cl.value("max", v);
      parms.set_max(v);
    }

    if (cl.option_present("d"))
    {
      int c;
      cl.value("d", c);
      parms.set_ndigits(c);
      set_default_iwstring_double_concatenation_precision(c);
    }

    if (verbose)
      parms.debug_print(cerr);

    if (! test_iwdigits_fraction_as_string_speed_test(parms, fname))
    {
      cerr << "Speed test to " << fname << " failed\n";
      return 1;
    }

    return 0;
  }

  int outer_loop_tests = 1000;

  if (cl.option_present('N'))
  {
    if (! cl.value('N', outer_loop_tests) || outer_loop_tests < 1)
    {
      cerr << "Invalid outer loop ntest value (-N)\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will perform " << outer_loop_tests << " outer loop tests\n";
  }

  int inner_loop_tests = 1000;

  if (cl.option_present('n'))
  {
    if (! cl.value('n', inner_loop_tests) || inner_loop_tests < 1)
    {
      cerr << "Invalid inner loop ntest value (-n)\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will perform " << inner_loop_tests << " inner loop tests\n";
  }

  int max_size_to_test = 1000;
  if (cl.option_present('s'))
  {
    if (! cl.value('s', max_size_to_test) || max_size_to_test < 1)
    {
      cerr << "The max size to test must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will test array sizes as large as " << max_size_to_test << endl;
  }

  std::mt19937_64 rng(rd());

  std::uniform_int_distribution<int> u(2, max_size_to_test);

  int failures = 0;

  for (auto i = 0; i < outer_loop_tests; ++i)
  {
    int s = u(rng);

    if (s > -2 && s < 2)
      continue;

    IWString sbefore, safter;

    if (! test_iwdigits(inner_loop_tests, s, sbefore, safter))
      failures++;

    sbefore = 'q';
    if (! test_iwdigits(inner_loop_tests, s, sbefore, safter))
      failures++;

    sbefore = "";
    safter = '_';
    if (! test_iwdigits(inner_loop_tests, s, sbefore, safter))
      failures++;

    sbefore = '%';
    safter = '_';
    if (! test_iwdigits(inner_loop_tests, s, sbefore, safter))
      failures++;
  }

  if (verbose)
  {
    if (failures)
      cerr << "encountered " << failures << " failures\n";
    return 1;
  }

  return 0;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_iwdigits(argc, argv);

  return rc;
}
