/*
  Need to check my T test stuff
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "iwpvalue.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::numeric_limits;

const char * prog_name = nullptr;

static int verbose = 0;

static int header_records_to_skip = 0;

static int write_h_value = 0;

static int records_to_read = numeric_limits<int>::max();

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes paired T tests on series of numbers\n";
  cerr << "Reads either one or two files\n";
  cerr << " -c <column>    column for values - same number as number of files on cmd line\n";
  cerr << " -s <nskip>     header records to skip\n";
  cerr << " -h             also write h value (1.0 - p/2.0)\n";
  cerr << " -r <nrec>      number records to read\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

#ifdef NOT_BEING_USED
static int
test_t_test (iwstring_data_source & input,
             ostream & output)
{
  Accumulator<double> acc1, acc2;
  Accumulator<double> diff;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (buffer.nwords() < 2)
    {
      cerr << "Must have at least two words on each line '" << buffer << "'\n";
      break;
    }

    const_IWSubstring s1, s2;

    (void) buffer.split (s1, ' ', s2);

    double d1, d2;

    if (! s1.numeric_value (d1) || ! s2.numeric_value (d2))
    {
      cerr << "Invalid numeric values\n";
      return 0;
    }

    acc1.extra (d1);
    acc2.extra (d2);

    diff.extra (d1 - d2);
  }

  int n = diff.n();
  if (n < 2)
    return 0;
  return output.good();
}
#endif


static int
do_the_computation (const Accumulator<double> & acc1,
                    const IWString & colname1,
                    const Accumulator<double> & acc2,
                    const IWString & colname2,
                    const Accumulator<double> & diff,
                    ostream & output)
{
  int n = diff.n();

  if (n <= 2)
  {
    cerr << "No data\n";
    return 0;
  }

  if (verbose)
  {
    cerr << colname1 << ' ' << acc1.minval() << " " << acc1.maxval() << " ave " << acc1.average() << endl;
    cerr << colname2 << ' ' << acc2.minval() << " " << acc2.maxval() << " ave " << acc2.average() << endl;
  }

  double v = sqrt (diff.variance());

  double t = diff.average() / (v / sqrt (static_cast<double> (n)));

  if (verbose)
    cerr << "Ave " << diff.average() << " variance " << v << " n = " << n << endl;

  double p = iwpvalue (n, fabs(t));

  if (verbose)
    output << "P value for " << n << " values t = " << fabs(t) << " " << p << endl;
  else
    output << "P " << p << endl;

  if (write_h_value)
  {
    double h = 1.0 - p * 0.5;

    output << "h = " << h << endl;
  }

  return output.good();
}

static int
get_numeric_value (const const_IWSubstring & buffer,
                   const int col,
                   double & v)
{
  int c = 0;
  int i = 0;
  const_IWSubstring token;

  while (buffer.nextword(token, i))
  {
    if (col == c)
    {
      if (! token.numeric_value(v))
      {
        cerr << "Invalid numeric '" << token << "'\n";
        return 0;
      }

      return 1;
    }

    c++;
  }

  cerr << "No column " << (col + 1) << " in '" << buffer << "'\n";
  return 0;
}

static int
test_t_test (iwstring_data_source & input1,
             int col1,
             iwstring_data_source & input2, 
             int col2,
             ostream & output)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input1.next_record(buffer))
    {
      cerr << "premature EOF on file1\n";
      return 0;
    }
    if (! input2.next_record(buffer))
    {
      cerr << "premature EOF on file2\n";
      return 0;
    }
  }

  Accumulator<double> acc1, acc2;
  Accumulator<double> diff;

  while (input1.next_record(buffer))
  {
    double v1;
    if (! get_numeric_value(buffer, col1, v1))
    {
      cerr << "Invalid input (file 1) '" << buffer << "' in column " << (col1 + 1) << endl;
      return 0;
    }

    if (! input1.next_record(buffer))
    {
      cerr << "Premature EOF on file2\n";
      return 0;
    }

    double v2;
    if (! get_numeric_value (buffer, col2, v2))
    {
      cerr << "Invalid input (file 2) '" << buffer << "' in column " << (col2 + 1) << endl;
      return 0;
    }

    acc1.extra(v1);
    acc2.extra(v2);
    diff.extra (v1 - v2);
  }

  IWString colname1, colname2;

  colname1 << "Col " << col1;
  colname2 << "Col " << col2;

  return do_the_computation (acc1, colname1, acc2, colname2, diff, output);
}

static int
test_t_test (iwstring_data_source & input,
             int col1,
             int col2, 
             ostream & output)
{
  const_IWSubstring buffer;

  IWString d1, d2;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot read " << header_records_to_skip << " header records\n";
      return 0;
    }

    if (0 == i)
    {
      buffer.word(col1, d1);
      buffer.word(col2, d2);
    }
//  cerr << "'" << d1 << " in column " << col1 << " and '" << d2 << "' in column " << col2 << endl;
  }

  Accumulator<double> acc1, acc2;
  Accumulator<double> diff;

  while (input.next_record (buffer))
  {
    double v1, v2;

    if (! get_numeric_value (buffer, col1, v1))
    {
      cerr << "Invalid input '" << buffer << "' in column " << (col1 + 1) << endl;
      cerr << __LINE__ << endl;
      return 0;
    }

    if (! get_numeric_value (buffer, col2, v2))
    {
      cerr << "Invalid input '" << buffer << "' in column " << (col2 + 1) << endl;
      cerr << __LINE__ << endl;
      return 0;
    }

    acc1.extra(v1);
    acc2.extra(v2);
    diff.extra (v1 - v2);

    if (diff.n() >= static_cast<uint32_t>(records_to_read))
      break;
  }

  if (d1.length () && d2.length())
    output << d1 << ' ' << static_cast<float>(acc1.average()) << ' ' << d2 << ' ' << static_cast<float>(acc2.average()) << ' ';

  return do_the_computation (acc1, d1, acc2, d2, diff, output);
}


static int
test_t_test (const char * fname1,
             int col1,
             int col2, 
             ostream & output)
{
  iwstring_data_source input1(fname1);

  if (! input1.good())
  {
    cerr << "Cannot open '" << fname1 << "'\n";
    return 0;
  }

  return test_t_test (input1, col1, col2, output);
}

static int
test_t_test (const char * fname1,
             int col1,
             const char * fname2,
             int col2, 
             ostream & output)
{
  iwstring_data_source input1(fname1);

  if (! input1.good())
  {
    cerr << "Cannot open '" << fname1 << "'\n";
    return 0;
  }

  iwstring_data_source input2(fname2);

  if (! input2.good())
  {
    cerr << "Cannot open '" << fname2 << "'\n";
    return 0;
  }

  return test_t_test (input1, col1, input2, col2, output);
}

#ifdef NOT_BEING_USED
static int
test_t_test (const char * fname,
             ostream & output)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_t_test (input, output);
}
#endif

static int
test_t_test (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:t:c:hs:r:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present('h'))
  {
    write_h_value = 1;

    if (verbose)
      cerr << "Will also write H values\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', header_records_to_skip) || header_records_to_skip < 1)
    {
      cerr << "The number of header records to skip (-s) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will skip the first " << header_records_to_skip << " records of each file\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', records_to_read) || records_to_read < 2)
    {
      cerr << "The number of records to read (-r) option must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will read " << records_to_read << " records\n";
  }

  if (cl.option_present ('n'))
  {
    int j = 0;
    int n;
    while (cl.value ('n', n, j++))
    {
      for (int i = 0; i < 200; i++)
      {
        double t = i / 40.0;
        double p = iwpvalue (n, t);

        double h = 1.0 - p * 0.5;

        cerr << " t = " << t << " N = " << n << " P = " << p << " h = " << h << endl;
      }
    }
  }

  if (cl.option_present ('t'))
  {
    int j = 0;
    float t;
    while (cl.value ('t', t, j++))
    {
      for (int n = 1; n < 30; n++)
      {
        double p = iwpvalue (n, t);

        double h = 1.0 - p * 0.5;

        cerr << " t = " << t << " N = " << n << " P = " << p << " h = " << h << endl;
      }
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;

  int column1 = 0;
  int column2 = 1;

  if (1 == cl.number_elements() && cl.option_count('c') > 1)  // multiple columns in same file
  {
    resizable_array<int> cols;
    int i = 0;
    int c;
    while (cl.value('c', c, i++))
    {
      c--;
//    cerr << "Got column " << c << endl;
      if (c < 1 || cols.contains(c))
      {
        cerr << "Invalid column specification '" << (c+1) << "'\n";
        if (cols.contains(c))
          cerr << "Already specified column " << (c+1) << endl;
        return 3;
      }

      cols.add(c);
    }

    int nc = cols.number_elements();
    for (int i = 1; i < nc && 0 == rc; i++)
    {
      int ci = cols[i];
      if (! test_t_test (cl[0], cols[0], ci, cout))
      {
        rc = 3;
        break;
      }
    }
  }
  else if (1 == cl.number_elements())
  {
    if (0 == cl.option_count('c'))
      ;
    else if (1 == cl.option_count('c'))
    {
      if (! cl.value('c', column1) || column1 < 1)
      {
        cerr << "Invalid column (-c)\n";
        usage(3);
      }

      if (verbose)
        cerr << "Data in column " << column1 << " and  column " << (column1 + 1) << endl;

      column1--;
      column2 = column1 + 1;
    }
    else
    {
      cerr << "Must have same number of -c options as files on command line\n";
      usage(4);
    }
    rc = test_t_test (cl[0], column1, column2, cout);
  }
  else if (2 == cl.number_elements())
  {
    if (0 == cl.option_count('c'))
    {
      column1 = 0;
      column2 = 1;
    }
    else if (1 == cl.option_count('c'))
    {
      if (! cl.value('c', column1) || column1 < 1)
      {
        cerr << "The column number option (-c) must be a whole +ve number\n";
        usage(3);
      }

      if (verbose)
        cerr << "Data in column " << column1 << endl;

      column1 --;
      column2 = column1;
    }
    else if (2 == cl.option_count('c'))
    {
      if (! cl.value('c', column1) || column1 < 1)
      {
        cerr << "Invalid column 1 (-c)\n";
        usage(3);
      }

      if (! cl.value('c', column2) || column2 < 1)
      {
        cerr << "Invalid column 2 (-c)\n";
        usage(3);
      }

      if (verbose)
        cerr << "Data from columns " << column1 << " and " << column2 << endl;

      column1--;
      column2--;

      rc = test_t_test (cl[0], column1, cl[1], column2, cout);
    }
    else
    {
      cerr << "Must specify as many -c options as files on the command line\n";
      usage(3);
    }
  }
  else
  {
    cerr << "Sorry, don't know how to handle " << cl.number_elements() << " files\n";
    return 3;
  }

  if (! rc)
    return 3;

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_t_test (argc, argv);

  return rc;
}
