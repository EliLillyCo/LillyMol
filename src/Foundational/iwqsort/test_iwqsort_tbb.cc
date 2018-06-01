#include <stdlib.h>

#include "tbb/task_scheduler_init.h"

#include "cmdline.h"
#include "iwrandom.h"
#include "iwqsort_tbb.h"

#include "foo_.h"

static int verbose = 0;

static int integer_range = 100;

static tbb::task_scheduler_init init;

static void
usage (int rc)
{
  cerr << " -n <niter>       number of iterations to perform\n";
  cerr << " -s <size>        size of arrays to sort\n";

  exit (rc);
}

template <typename T, typename Compare>
int
check_sorted (const T * f,
              int s,
              const Compare & c)
{
  int rc = 1;

  for (int i = 1; i < s; i++)
  {
    if (c(f[i], f[i - 1]))
    {
      cerr << "Sorted items out of order " << f[i - 1] << " followed by " << f[i] << endl;
      rc = 0;
    }
    if (verbose > 1)
      cerr << " i = " << i << " value " << f[i] << endl;
  }

  return rc;
}

template <typename T, typename Compare, typename InitialValueAssigner>
int
test_tbb_parallel_sort (T * f, 
                        int s,
                        int iterations,
                        Compare & c,
                        InitialValueAssigner & iva)
{
  for (int i = 0; i < iterations; i++)
  {
    for (int j = 0; j < s; j++)
    {
      iva(f[j]);
    }

    parallel_sort(f, f + s, c);
    check_sorted(f, s, c);
  }

  return 1;
}

class Random_Numbers_Int
{
  private:
    const int _maxval;

  public:
    Random_Numbers_Int (int maxval_) : _maxval(maxval_) {}

    int operator() (int & s) const
    {
      s = intbtwij(0, _maxval);
      return 1;
    };
};

class Random_Numbers_Float
{
  private:
    const float _range;
  public:
    Random_Numbers_Float (float range_) : _range(range_) {}

    int operator () (float & s) const
    {
      s = iwrandom() * _range;
      return 1;
    }
};

class Int_Comparator
{
  private:
  public:
    bool operator () (int i1, int i2) const
    {
      return i1 < i2;
    }
};

class Float_Comparator
{
  private:
  public:
    bool operator () (float i1, float i2) const
    {
      return i1 < i2;
    }
};

static int
test_tbb_parallel_sort_int (int s, int iterations)
{
  int * f = new int[s];

  Random_Numbers_Int rni(integer_range);
  Int_Comparator intc;

  int rc = test_tbb_parallel_sort (f, s, iterations, intc, rni);

  delete [] f;

  return rc;
}

static int
test_tbb_parallel_sort_float (int s, int iterations)
{
  float * f = new float[s];

  Random_Numbers_Float rnf(100.0);
  Float_Comparator floatc;

  int rc = test_tbb_parallel_sort (f, s, iterations, floatc, rnf);

  delete [] f;

  return rc;
}

static int
test_tbb_parallel_sort (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:i:fr:h");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present('h'))
  {
    usage(1);
  }

  if (cl.option_present ('r'))
  {
    int seed;
    if (! cl.value ('r', seed))
    {
      cerr << "Invalid value for random number seed\n";
      usage (5);
    }

    if (verbose)
      cerr << "Using seed " << seed << endl;

    iw_set_rnum_seed (seed);
  }
  else
    iw_random_seed ();

  int s;
  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', s) || s < 2)
    {
      cerr << "Must sort at least 2 items\n";
      usage (3);
    }
  }
  else
  {
    s = intbtwij (2, 1000);
  }

  if (verbose)
    cerr << "Will test vectors of length " << s << " items\n";

  int iterations = 1;
  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', iterations) || iterations < 1)
    {
      cerr << "The number of iterations must be a whole number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will perform " << iterations << " iterations\n";
  }

  if (cl.option_present('i'))
  {
    if (! cl.value ('i', integer_range) || integer_range < 1)
    {
      cerr << "Invalid maximum integer range (-i)\n";
      usage (5);
    }

    if (verbose)
      cerr << "Testing int's between 0 and " << integer_range << endl;
  }

  test_tbb_parallel_sort_int (s, iterations);
  test_tbb_parallel_sort_float (s, iterations);

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_tbb_parallel_sort (argc, argv);

  return rc;
}
