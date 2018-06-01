#include <stdlib.h>

#include "cmdline.h"
#include "iwrandom.h"
#include "mfn.h"

#include "foo_.h"

static int verbose = 0;

static void
usage (int rc)
{
  exit (rc);
}

template <typename F>
int
check_sorted (F * f,
              int s,
              int (*mfn) (F &, F &))
{
  int rc = 0;

  for (int i = 1; i < s; i++)
  {
//  cerr << "i = " << i << " value " << f[i].zvalue () << endl;
    int c = mfn (f[i], f[i - 1]);

    if (c > 1)
    {
      cerr << "Sorted items out of order " << f[i - 1].zvalue () << " followed by " << f[i].zvalue () << endl;
      rc = 0;
    }
  }

  return rc;
}

template <typename F, typename T>
int
test_iwqsort (F * f, 
              int s,
              int iterations,
              T * initial_values,
              int (*mfn) (F &, F &),
              void (&new_values) (T *, int))
{
  for (int i = 0; i < s; i++)
  {
    f[i].set_value (static_cast<T> (8));
  }

  if (verbose)
    cerr << "Constant values\n";

  iwqsort_mfn (f, s, mfn);

  check_sorted (f, s, mfn);

  for (int i = 0; i < s; i++)
  {
    f[i].set_value (static_cast<T> (i));
  }

  if (verbose)
    cerr << "Increasing values\n";

  iwqsort_mfn (f, s, mfn);

  check_sorted (f, s, mfn);

  for (int i = 0; i < s; i++)
  {
    f[i].set_value (static_cast<T> (s - i));
  }

  if (verbose)
    cerr << "Decreasing values\n";

  iwqsort_mfn (f, s, mfn);

  check_sorted (f, s, mfn);

  for (int i = 0; i < iterations; i++)
  {
    new_values (initial_values, s);
    for (int j = 0; j < s; j++)
    {
      f[j].set_value (initial_values[j]);
    }

    iwqsort_mfn (f, s, mfn);

    check_sorted (f, s, mfn);
  }

  return 1;
}

#ifdef __GNUG__q

template int test_iwqsort (Foo<int> *,   int, int, int *,     int (Foo<int>::*) (Foo<int> &, Foo<int> &), void (&) (int *,   int));
template int test_iwqsort (Foo<float> *,   int, int, int *,   int (Foo<float>::*) (Foo<float> &, Foo<float> &), void (&) (int *,   int));

template int check_sorted (const Foo<int> *, int, int (*) (Foo<int> &, Foo<int> &));
template int check_sorted (const Foo<float> *, int, int (*) (Foo<float> &, Foo<float> &));

template void iwqsort_mfn (Foo<int> *, int, int (Foo<int>::*) (const Foo<int> &, const Foo<int> &));
template void iwqsort_mfn (Foo<float> *, int, int (Foo<float>::*) (const Foo<float> &, const Foo<float> &));

#endif

static void
random_numbers_int (int * f, int s)
{
  for (int i = 0; i < s; i++)
  {
    f[i] = intbtwij (-1000, 1000);
  }

  return;
}

static int
test_iwqsort_int (int s, int iterations)
{
  Foo<int> * f = new Foo<int>[s];

  int * initial_values = new int[s];

  int rc = test_iwqsort (f, s, iterations, initial_values, &(Foo<int>::iwqsort_mfn), random_numbers_int);
//int rc = test_iwqsort (f, s, iterations, initial_values, &foo_comparitor_int, random_numbers_int);

  delete [] f;
  delete initial_values;

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

//int rc = test_iwqsort (f, s, iterations, initial_values, &foo_comparitor_float, random_numbers_float);
  int rc = test_iwqsort (f, s, iterations, initial_values, &(Foo<float>::iwqsort_mfn), random_numbers_float);
  
  delete [] f;
  delete initial_values;

  return rc;
}

static int
test_iwqsort (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:if");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

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

  if (cl.option_present ('f'))
  {
    test_iwqsort_float (s, iterations);
  }
  else
  {
    test_iwqsort_int (s, iterations);
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_iwqsort (argc, argv);

  return rc;
}
