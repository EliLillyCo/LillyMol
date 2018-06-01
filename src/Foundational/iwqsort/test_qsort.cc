#include <stdlib.h>

/*
  Test the system qsort
*/

#include "cmdline.h"
#include "iwrandom.h"
#define IWQSORT_IMPLEMENTATION
#include "iwqsort.h"

#include "foo_.h"

static int verbose = 0;

static int integer_range = 100;

static void
usage (int rc)
{
  exit (rc);
}

template <typename F>
int
check_sorted (const F * f,
              int s)
{
  int rc = 0;

  for (int i = 1; i < s; i++)
  {
//  cerr << "i = " << i << " value " << f[i].zvalue () << endl;
    int c = f[i].iwqsortcompare (f[i - 1]);
    if (c < 0)
    {
      cerr << "Sorted items out of order " << f[i - 1].zvalue () << " followed by " << f[i].zvalue () << endl;
      rc = 0;
    }
  }

  return rc;
}

static int
comparitor_function_int (const void * pf1, const void * pf2)
{
  const Foo<int> * f1 = (const Foo<int> *) (pf1);
  const Foo<int> * f2 = (const Foo<int> *) (pf2);

  return f1->iwqsortcompare (*f2);
}

static int
comparitor_function_float (const void * pf1, const void * pf2)
{
  const Foo<float> * f1 = (const Foo<float> *) (pf1);
  const Foo<float> * f2 = (const Foo<float> *) (pf2);

  return f1->iwqsortcompare (*f2);
}

template <typename F, typename T>
int
test_iwqsort (F * f, 
              int s,
              int iterations,
              T * initial_values,
              void (&new_values) (T *, int),
              int (* comparitor_function) (const void *, const void *))
{
  for (int i = 0; i < s; i++)
  {
    f[i].set_value (static_cast<T> (8));
  }

  if (verbose)
    cerr << "Constant values\n";

  qsort (f, s, sizeof (F), comparitor_function);

  check_sorted (f, s);

  for (int i = 0; i < s; i++)
  {
    f[i].set_value (static_cast<T> (i));
  }

  if (verbose)
    cerr << "Increasing values\n";

  qsort (f, s, sizeof (F), comparitor_function);

  check_sorted (f, s);

  for (int i = 0; i < s; i++)
  {
    f[i].set_value (static_cast<T> (s - i));
  }

  if (verbose)
    cerr << "Decreasing values\n";

  qsort (f, s, sizeof (F), comparitor_function);

  check_sorted (f, s);

  for (int i = 0; i < iterations; i++)
  {
    new_values (initial_values, s);
    for (int j = 0; j < s; j++)
    {
      f[j].set_value (initial_values[j]);
    }

    qsort (f, s, sizeof (F), comparitor_function);

    check_sorted (f, s);
  }

  return 1;
}

#ifdef __GNUG__

template int test_iwqsort (Foo<int> *,   int, int, int *,   void (&) (int *,   int), int (*) (const void *, const void *));
template int test_iwqsort (Foo<float> *, int, int, float *, void (&) (float *, int), int (*) (const void *, const void *));

template int check_sorted (const Foo<int> *, int);
template int check_sorted (const Foo<float> *, int);

template void iwqsort (Foo<int> *, int);
template void iwqsort (Foo<float> *, int);

#endif

static void
random_numbers_int (int * f, int s)
{
  for (int i = 0; i < s; i++)
  {
    f[i] = intbtwij (0, integer_range);
  }

  return;
}

static int
test_iwqsort_int (int s, int iterations)
{
  Foo<int> * f = new Foo<int>[s];

  int * initial_values = new int[s];

  int rc = test_iwqsort (f, s, iterations, initial_values, random_numbers_int, comparitor_function_int);

  delete initial_values;
  
  delete [] f;

  return rc;
}

static void
random_numbers_float (float * f, int s)
{
  for (int i = 0; i < s; i++)
  {
    f[i] = static_cast<float> (iwrandom ());
  }

  return;
}

static int
test_iwqsort_float (int s, int iterations)
{
  Foo<float> * f = new Foo<float>[s];

  float * initial_values = new float[s];

  int rc = test_iwqsort (f, s, iterations, initial_values, random_numbers_float, comparitor_function_float);
  
  delete initial_values;
  delete [] f;

  return rc;
}

static int
test_iwqsort (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:i:fr:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

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

  if (cl.option_present ('f') && cl.option_present ('i'))
  {
    cerr << "Can do either floating point or integer tests\n";
    usage (6);
  }

  if (cl.option_present ('i'))
  {
    if (! cl.value ('i', integer_range) || integer_range < 1)
    {
      cerr << "Invalid maximum integer range (-i)\n";
      usage (5);
    }

    if (verbose)
      cerr << "Testing int's between 0 and " << integer_range << endl;

    test_iwqsort_int (s, iterations);
  }
  else
  {
    test_iwqsort_float (s, iterations);
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_iwqsort (argc, argv);

  return rc;
}
