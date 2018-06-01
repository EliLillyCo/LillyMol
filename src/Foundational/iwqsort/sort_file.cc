/*
  Reads records from a file and sorts them - tester for iwqsort
*/

#include <stdlib.h>

#include "iwstring_data_source.h"
#include "cmdline.h"
#define IWQSORT_IMPLEMENTATION
#include "iwqsort.h"
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << " -v             verbose output\n";

  exit (rc);
}

template <typename T>
class SValue
{
  private:
    T _v;

  public:
    SValue ();

    void set_value (T v) { _v = v;}
    T    zvalue () const { return _v;}
    int iwqsortcompare (SValue &);
};

template <typename T>
SValue<T>::SValue ()
{
  _v = static_cast<T> (0);

  return;
}

template <typename T>
int
SValue<T>::iwqsortcompare (SValue & rhs)
{
  if (_v < rhs._v)
    return -1;

  if (_v > rhs._v)
    return 1;

  return 0;
}

template <typename T>
int
test_iwqsort (SValue<T> * sv,
              int n,
              ostream & output)
{
  iwqsort (sv, n);

  for (int i = 1; i < n; i++)
  {
    int c = sv[i].iwqsortcompare (sv[i - 1]);
    if (c >= 0)
      continue;

    cerr << "Results out of order, i = " << i << " value " << sv[i].zvalue () << " vs " << sv[i - 1].zvalue () << endl;
    for (int j = 0; j < n; j++)
    {
      cerr << " i = " << j << " value " << sv[j].zvalue () << endl;
    }

    return 0;
  }

  return 1;
}

template <typename T>
int
test_iwqsort (const const_IWSubstring & buffer,
              SValue<T> * sv,
              int n,
              ostream & output)
{
  int i = 0;
  const_IWSubstring token;
  int col = 0;
  while (buffer.nextword (token, i))
  {
    T v;
    if (! token.numeric_value (v))
    {
      cerr << "Invalid numeric '" << token << "'\n";
      return 0;
    }

    sv[col].set_value (v);
    col++;
  }

  assert (col == n);

  return test_iwqsort (sv, n, output);
}

static int
test_iwqsort (const const_IWSubstring & buffer,
              ostream & output)
{
  int n = buffer.nwords ();
  if (n < 2)
  {
    cerr << "Need at least two tokens to sort\n";
    return 0;
  }

  SValue<int> * sv = new SValue<int>[n];

  int rc = test_iwqsort (buffer, sv, n, output);

  delete [] sv;

  return rc;
}

static int
test_iwqsort (iwstring_data_source & input,
              ostream & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! test_iwqsort (buffer, output))
    {
      cerr << "Fatal error on line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return output.good ();
}

static int
test_iwqsort (const char * fname, ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_iwqsort (input, output);
}

static int
test_iwqsort (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! test_iwqsort (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_iwqsort (argc, argv);

  return rc;
}
