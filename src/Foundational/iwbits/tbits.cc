#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwbits.h"
#include "iwrandom.h"
#include "iwstring.h"

const char * prog_name = nullptr;

static int verbose = 0;

/*
  We do ITERATIONS iterations of various tests. To exercise things more, we
  can change the bit sizes used every CHANGE_BIT_SIZE iterations
*/

static int change_bit_size = 1000;

/*
  When doing tests, we use vectors which are of a length
  between lower_vector_size and upper_vector_size.
*/

static int lower_vector_size = 1;
static int upper_vector_size = 1000;

static int iterations = 100;

/*
  Some tests take relatively longer
*/

static int fewer_iterations = 10;
static int fewer_iterations_split = 5;

/*
  There are some tests for which we divide the iterations between an inner and
  an outer loop
*/

static int iterations_split = 10;

static int failures = 0;

static void
usage (int rc)
{
  cerr << "Usage " << prog_name << " <options>\n";
  cerr << " -s <seed>              specify random number seed (* for random seed)\n";
  cerr << " -i <iterations>        for tests which iterate, specify iteration count\n";
  cerr << " -b <nbits>             smallest size for test vectors\n";
  cerr << " -B <nbits>             largest size for test vectors\n";
  cerr << " -p                     print the bit masks\n";
  cerr << " -v                     verbose output\n";

  exit (rc);
}

static int
problem_apr03 ()
{
  IW_Bits_Base foo = IW_Bits_Base (6);
  foo.set (3);
  int previous = 4;
  int current = foo.next_on_bit (previous);

  if (current >= 0)
  {
    cerr << "Next on bit failed ";
    foo.printon (cerr);
    cerr << " got " << current << " expected negative\n";
    failures++;
  }

  IW_Bits_Base foo2 (15);
  foo2.set (12);
  previous = 13;
  current = foo2.next_on_bit (previous);
  if (current >= 0)
  {
    cerr << "Next on bit failed ";
    foo.printon (cerr);
    cerr << " got " << current << " expected negative\n";
    failures++;
  }

  foo2.set (12, 0);
  foo2.set (5);
  previous = 0;
  current = foo2.next_on_bit (previous);
  if (5 != current)
  {
    cerr << "Next on bit failed ";
    foo.printon (cerr);
    cerr << " got " << current << " expected 5\n";
    failures++;
  }

  current = foo2.next_on_bit (previous);
  if (current >= 0)
  {
    cerr << "Next on bit failed ";
    foo.printon (cerr);
    cerr << " got " << current << " expected negative\n";
    failures++;
  }

  return 1;
}

static int
test_reading_daylight_ascii_representation (const const_IWSubstring & buffer)
{
  if (! buffer.starts_with ("FP"))
    return 1;

  if (! buffer.ends_with ('>'))
    return 1;

  IWString mybuffer (buffer);

  mybuffer.remove_up_to_first ('<');
  mybuffer.chop ();

  IWString zrep, zextra;

  if (! mybuffer.split (zrep, ';', zextra))
  {
    cerr << "Cannot split into representation and data\n";
    cerr << mybuffer << endl;
    failures++;
    return 0;
  }

  IW_Bits_Base b;

  if (! b.construct_from_daylight_ascii_bit_rep (zrep.rawchars (), zrep.length ()))
  {
    cerr << "Cannot read Daylight bit representation\n";
    cerr << mybuffer << endl;
    failures++;
    return 0;
  }

  int initial_nbits, initial_nset;
  int final_nbits, final_nset;
  if (4 != sscanf (zextra.null_terminated_chars (), "%d;%d;%d;%d", &initial_nbits, & initial_nset, &final_nbits, &final_nset))
  {
    cerr << "Yipes, cannot extract nset and such\n";
    cerr << mybuffer << endl;
    failures++;
    return 0;
  }

  if (final_nbits != b.nbits ())
  {
    cerr << "test_reading_daylight_ascii_representation: expected " << final_nbits << " got " << b.nbits () << " bits\n";
    cerr << buffer << endl;
    failures++;
    return 0;
  }

  if (final_nset != b.nset ())
  {
    cerr << "test_reading_daylight_ascii_representation: expected " << final_nbits << " got " << b.nbits () << " bits set\n";
    cerr << buffer << endl;
    failures++;
    return 0;
  }

  return 1;
}

static int
test_reading_daylight_ascii_representation (ifstream & input)
{
  IWString buffer;
  while (buffer.getline (input))
  {
    if (! test_reading_daylight_ascii_representation (buffer))
      return 0;
  }

  return 1;
}

static int
test_reading_daylight_ascii_representation_file (const char * fname)
{
  ifstream input;
  input.open (fname, ios::in);

  if (! input.good ())
  {
    cerr << "test_reading_daylight_ascii_representation: cannot open '" << fname << "'\n";
    return 0;
  }

  return test_reading_daylight_ascii_representation (input);
}

static void
set_random_bits (IW_Bits_Base & b)
{
  int nb = b.nbits ();

  for (int i = 0; i < nb; i++)
  {
    if (intbtwij (0, 1))
      b.set (i, 1);
  }

  return;
}

static int
test_iterator ()
{
  for (int i = 0; i < iterations; i++)
  {
    IW_Bits_Base b;

    int nb = intbtwij (lower_vector_size, upper_vector_size);

    set_random_bits (b);

    int nset = b.nset ();
    if (0 == nset)
      continue;

    int * tmp = new int[nb];

    b.set_vector (tmp);

    int j = 0;

    int number_bits_set = 0;

    int bit_set;
    while ((bit_set = b.next_on_bit (j)) >= 0)
    {
      if (0 == tmp[bit_set])
      {
        cerr << "next_set, bit " << bit_set << " not set\n";
        failures++;
      }

      number_bits_set++;
    }

    delete tmp;

    if (number_bits_set != b.nset ())
    {
      cerr << "next_set:bit count mismatch " << number_bits_set << " vs " << b.nset () << endl;
    }
  }
  
  return 1;
}

static int
test_hex_stuff ()
{
  IW_Bits_Base b1;

  int nb = intbtwij (4, 2048);
  if (0 != nb % 8)
    nb = nb / 8 * 8;

  b1.allocate_space_for_bits (nb);

  for (int i = 0; i < iterations; i++)
  {
    set_random_bits (b1);

    IWString zhex;

    b1.hex_form (zhex);

    if (zhex.length () != b1.nbits () / 4)
    {
      cerr << "Hmm, nbits = " << b1.nbits () << " but hex form " << zhex.length () << " chars\n";
      cerr << zhex << endl;
      failures++;
      return 0;
    }

    IW_Bits_Base b2;

    if (! b2.construct_from_hex (zhex))
    {
      cerr << "Cannot parse hex form '" << zhex << "'\n";
      return 0;
    }

    if (b1 != b2)
    {
      cerr << "Hex echo failed '" << zhex << "'\n";
      cerr << "nbits = " << b1.nbits () << " b2 bits = " << b2.nbits () << endl;

      for (int i = 0; i < b1.nbits (); i++)
      {
        int s1 = b1.is_set (i);
        int s2 = b2.is_set (i);
        if (s1 != s2)
          cerr << "Bit " << i << " 1 = " << s1 << " 2 = " << s2 << endl;
      }

      failures++;

      return 0;
    }
  }

  return 1;
}

/*
  We want to check that all the bits in B correspond to a 1 in IS_SET
*/

static int
check_bits_set (const IW_Bits_Base & b,
                int nb, 
                int * is_set)
{
  assert (nb = b.nbits ());

  b.increment_vector (is_set);

  int rc = 1;
  for (int i = 0; i < nb; i++)
  {
    if (0 == is_set[i])
      continue;

    if (2 == is_set[i])
      continue;

    int s = b.is_set (i);
    cerr << "Error on bit " << i << " in object " << s << " in array " << is_set[i] << endl;
    rc = 0;
  }

  return rc;
}

static void
build_string_rep (IWString & string_rep,
                  int nb,
                  int * is_set)
{
  for (int i = 0; i < nb; i++)
  {
    is_set[i] = 0;
  }

  int n1 = intbtwij (0, nb / 2);    // start somewhere low

  while (1)
  {
    is_set[n1] = 1;

    if (string_rep.length ())
      string_rep += ',';

    string_rep << n1;

    int range = intbtwij (0, 1);

    if (1 == range)     // do a range
    {
      int n2 = intbtwij (n1, nb - 1);

      string_rep << '-' << n2;

      for (int i = n1; i <= n2; i++)
      {
        is_set[i] = 1;
      }

      n1 = intbtwij (n2 + 1, nb - 1);
    }
    else
    {
      n1 = intbtwij (n1 + 1, nb - 1);
    }

    if (n1 >= nb - 1)
      break;
  }

  string_rep << ';' << nb;

  return;
}

static int
test_sparse_things (int nb, int * tmp)
{
  for (int i = 0; i < fewer_iterations_split; i++)
  {
    IWString string_rep;

    build_string_rep (string_rep, nb, tmp);

    IW_Bits_Base b;

    b.allocate_space_for_bits (nb);

    if (! b.construct_from_sparse_representation (string_rep))
    {
      cerr << "Cannot parse string representation '" << string_rep << "'\n";
      return 0;
    }

    if (! check_bits_set (b, nb, tmp))
    {
      cerr << "Sparse bits representation '" << string_rep << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
test_sparse_things ()
{
  int * tmp = nullptr;

  int bsize = 0;
  for (int i = 0; i < fewer_iterations_split; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (1, 2048);

      if (verbose)
        cerr << "test_sparse_things bsize " << bsize << endl;

      if (nullptr != tmp)
        delete tmp;

      tmp = new int[bsize];
    }

    (void) test_sparse_things (bsize, tmp);
  }

  if (nullptr != tmp)
    delete tmp;

  return 1;
}

static int
test_ascii_01_things (int nb)
{
  IW_Bits_Base b (nb);

  for (int i = 0; i < fewer_iterations_split; i++)
  {
    set_random_bits (b);

    IWString ascii;

    b.append_ascii_01_representation (ascii);

    if (ascii.length () != b.nbits ())
    {
      cerr << "Yipes, bit vector has " << b.nbits () << " bits, but string has " << ascii.nchars () << " characters\n";
      failures++;
      return 0;
    }

    IW_Bits_Base b2;

    if (! b2.construct_from_ascii_01_representation (ascii.rawchars (), ascii.length ()))
    {
      cerr << "Cannot parse ASCII '" << ascii << "'\n";
      failures++;
      return 0;
    }

    if (b2 != b)
    {
      cerr << "Yipes, mismatch via ASCII 01 conversion\n";
      cerr << ascii << endl;
      b.debug_print (cerr);
      b2.debug_print (cerr);
      failures++;
      return 0;
    }

    b.clear ();
  }

  return 1;
}

static int
test_ascii_01_things ()
{
  int bsize = 0;
  for (int i = 0; i < fewer_iterations_split; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (1, 2048);

      if (verbose)
        cerr << "test_ascii_01_things bsize " << bsize << endl;
    }

    test_ascii_01_things (bsize);
  }

  return 1;
}

static int
test_shift ()
{
  int bsize = 512;

  for (int i = 1; i <= 11; i++)
  {
    IW_Bits_Base b;
    b.allocate_space_for_bits (bsize);

    for (int j = 0; j < bsize; j++)
    {
      if (0 == j % i)
        b.set (j, 1);
    }

    b.shift (i + 1);

    assert (bsize + i + 1 == b.nbits ());

    for (int j = 0; j < i + 1; j++)
    {
      if (b.is_set (j))
      {
        cerr << "Leading whitespace bit " << j << " set\n";
        failures++;
        return 0;
      }
    }

    for (int j = i + 1; j < b.nbits (); j++)
    {
      if (1 == j % i)
      {
        if (! b.is_set (j))
        {
          cerr << "Yipes, bit " << j << " not set in shift " << i << " test\n";
          failures++;
        }
      }
      else if (i > 1)
      {
        if (b.is_set (j))
        {
          cerr << "Yipes, bit " << j << " set in shift " << i << " test\n";
          failures++;
        }
      }
    }
  }

  return 1;
}

static int
test_set_bit_vector1 (const IW_Bits_Base & b,
                      int istart, int istop)
{
  for (int i = 0; i < istart; i++)
  {
    if (b.is_set (i))
    {
      cerr << "After setting range " << istart << " to " << istop << " bit " << i << " is set\n";
      failures++;
    }
  }

  for (int i = istart; i <= istop; i++)
  {
    if (! b.is_set (i))
    {
      cerr << "After setting range " << istart << " to " << istop << " bit " << i << " is not set\n";
      failures++;
    }
  }

  for (int i = istop + 1; i < b.nbits (); i++)
  {
    if (b.is_set (i))
    {
      cerr << "After setting range " << istart << " to " << istop << " bit " << i << " is set\n";
      failures++;
    }
  }

  return 1;
}

static int
test_set_bit_vector0 (const IW_Bits_Base & b,
                      int istart, int istop)
{
  for (int i = 0; i < istart; i++)
  {
    if (! b.is_set (i))
    {
      cerr << "After zeroing range " << istart << " to " << istop << " bit " << i << " is not set\n";
      failures++;
    }
  }

  for (int i = istart; i <= istop; i++)
  {
    if (b.is_set (i))
    {
      cerr << "After zeroing range " << istart << " to " << istop << " bit " << i << " is set\n";
      failures++;
    }
  }

  for (int i = istop + 1; i < b.nbits (); i++)
  {
    if (! b.is_set (i))
    {
      cerr << "After zeroing range " << istart << " to " << istop << " bit " << i << " is not set\n";
      failures++;
    }
  }

  return 1;
}

static int
test_printon ()
{
  for (int i = 0; i < iterations; i++)
  {
    int nb = intbtwij (24, 80);

    IW_Bits_Base b (nb);

    set_random_bits (b);

    IWString s1, s2;
  
    b.append_string_form (s1, '1', '0', 1);
    b.append_string_form_fast (s2, 1);
  
    if (s1 != s2)
    {
      cerr << "Yipes, append string form failure, nb = " << nb << endl;
      cerr << "slow " << s1 << endl;
      cerr << "fast " << s2 << endl;

      return 0;
    }
  }

  return 1;
}

static int
test_set_bit_vector (const IW_Bits_Base & b,
                     int istart, int istop,
                     int value_in_range)
{
  if (value_in_range)
    return test_set_bit_vector1 (b, istart, istop);
  else
    return test_set_bit_vector0 (b, istart, istop);
}

static int
test_construct_from_array_of_ints (int bsize, int * ii)
{
  IW_Bits_Base foo;
  for (int i = 0; i < bsize; i++)
  {
    ii[i] = intbtwij (0, 1);
  }

  foo.construct_from_array_of_ints (ii, bsize);

  if (foo.nbits () != bsize)
  {
    cerr << "Bits constructed from vector has wrong count " << foo.nbits () << endl;
    failures++;
  }

  for (int i = 0; i < bsize; i++)
  {
    if (ii[i] && foo.is_set (i))
      ;
    else if (0 == ii[i] && ! foo.is_set (i))
      ;
    else
    {
      cerr << "Construct from array of bits failed, i = " << i << " array is " <<
              ii[i] << " set is " << foo.is_set (i) << endl;
      failures++;
    }
  }

  return 0;
}

static int
test_set_all (int bsize)
{
  IW_Bits_Base foo (bsize);

  foo.set_all ();

  for (int i = 0; i < bsize; i++)
  {
    if (! foo.is_set (i))
    {
      cerr << "bsize " << bsize << " just set all bits but bit " << i << " not set\n";
      failures++;
    }
  }

  return 1;
}

int 
test_set_all ()
{
  int bsize = 0;

  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "test_set_all using vectors of length " << bsize << endl;
    }

    test_set_all (bsize);
  }

  return 1;
}

static int
test_first_bit ()
{
  for (int i = 0; i < 132; i++)
  {
    for (int j = 0; j < i; j++)
    {
      IW_Bits_Base foo (i);

      foo.set (j);

      if (j != foo.first_bit ())
      {
        cerr << "First bit failure, nbits " << i << " bit " << j << " first found to be " << foo.first_bit() << endl;
        failures++;
      }
    }
  }

  for (int i = 0; i < iterations; i++)
  {
  }

  return 1;
}

static int
test_setting_ranges (int bsize)
{
  for (int i = 0; i < fewer_iterations_split; i++)
  {
    IW_Bits_Base foo (bsize);

    int bstart = intbtwij (0, bsize - 1);
    int bstop  = intbtwij (bstart, bsize - 1);

    if (0 == i % 2)    // we will be filling 0's, so fill bit vector with 1's
      foo.set_all ();

    foo.set_all_bits (bstart, bstop, i % 2);

    test_set_bit_vector (foo, bstart, bstop, i % 2);
  }

// There are some specially difficult cases. Make sure we hit them

  IW_Bits_Base foo (bsize + 32);
  foo.set_all_bits (0, 31, 1);

  test_set_bit_vector (foo, 0, 31, 1);

  foo.set_all ();
  foo.set_all_bits (0, 31, 0);
  test_set_bit_vector (foo, 0, 31, 0);

  foo.clear ();
  foo.set_all_bits (0, 32, 1);
  test_set_bit_vector (foo, 0, 32, 1);

// Make sure we test operations in a single word

  for (int i = 0; i < iterations_split; i++)
  {
    if (i % 2)
      foo.clear ();
    else
      foo.set_all ();

    int i1 = intbtwij (0, foo.nbits () - 2);

    int i2 = i1 + intbtwij (1, 31);
    if (i2 > foo.nbits () - 1)
      i2 = foo.nbits () - 1;

    foo.set_all_bits (i1, i2, i % 2);
    test_set_bit_vector (foo, i1, i2, i % 2);
  }

  IW_Bits_Base bar (bsize + 64);

  bar.set_all_bits (32, 63, 1);

  test_set_bit_vector (bar, 32, 63, 1);

  return 1;
}

static int
test_setting_ranges ()
{
  int bsize = 0;

  for (int i = 0; i < fewer_iterations_split; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "test_setting_ranges using vectors of length " << bsize << endl;
    }

    test_setting_ranges (bsize);
  }

  return 1;
}

static int
test_construct_from_array_of_ints ()
{
  int * ii = nullptr;

  int bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "test_construct_from_array_of_ints using vectors of length " << bsize << endl;

      if (nullptr != ii)
        delete ii;

      ii = new int[bsize];
    }

    test_construct_from_array_of_ints (bsize, ii);
  }

  if (nullptr != ii)
    delete ii;

  return 0;
}

static int
test_compute_weight (const int bsize, float * x)
{
  IW_Bits_Base foo (bsize);

  float expected_result = 0.0;
  for (int i = 0; i < bsize; i++)
  {
    int new_value = intbtwij (0, 1);
    x[i] = iwrandom () * 10.0;    // try to keep reasonable magnitude numbers
    if (new_value)
    {
      expected_result += x[i];
      foo.set (i, 1);
    }
  }

  float result = foo.compute_weight (x);
  if (result == expected_result)
    ;
  else if (fabs (result - expected_result) > 0.1)
  {
    cerr << "Error on compute weight, expected " << expected_result << " got " << result << " diff " << fabs (expected_result - result) << endl;
    cerr << "nbits = " << foo.nbits () << endl;
    failures++;
  }

  return 0;
}
static int
test_compute_weight ()
{  
  float * x = nullptr;

  int bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "test_compute_weight using vectors of length " << bsize << endl;

      if (nullptr != x)
        delete x;

      x = new float[bsize];
    }

    test_compute_weight (bsize, x);
  }

  if (nullptr != x)
    delete x;

  return 0;
}

static int
test_increment_vector (int bsize,
                       int * i1,
                       int * i2)
{
  IW_Bits_Base foo (bsize);

  for (int i = 0; i < bsize; i++)
  {
    int new_value = intbtwij (0, iterations);    // iterations just used as a random number
    i1[i] = i2[i] = new_value;
    if (new_value)
      foo.set (i, 1);
  }

  foo.increment_vector (i1);

  for (int i = 0; i < bsize; i++)
  {
    if (foo.is_set (i) && i1[i] == i2[i] + 1)
      ;
    else if (! foo.is_set (i) && i1[i] == i2[i])
      ;
    else
    {
      cerr << "Increment vector failed, bit " << i << " set = " << foo.is_set (i) <<
              " i1 = " << i1[i] << " i2 = " << i2[i] << endl;
    }
  }

  return 0;
}

static int
test_increment_vector ()
{
  int * i1 = nullptr;
  int * i2 = nullptr;

  int bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "test_increment_vector using vectors of length " << bsize << endl;

      if (nullptr != i1)
        delete i1;
      if (nullptr != i2)
        delete i2;

      i1 = new int[bsize];
      i2 = new int[bsize];
    }

    test_increment_vector (bsize, i1, i2);
  }

  if (nullptr != i1)
    delete i1;
  if (nullptr != i2)
    delete i2;

  return 0;
}

static int
test_set_vector (const int bsize, int * ii,
                 int * expected_result)
{
  IW_Bits_Base foo (bsize);

  for (int i = 0; i < bsize; i++)
  {
    int new_value = intbtwij (0, 1);
    if (new_value)
      foo.set (i, new_value);
    expected_result[i] = new_value;

    if (new_value && ! foo.is_set (i))
    {
      cerr << "Yipes, just set bit " << i << " but now not set\n";
      failures++;
    }
    else if (0 == new_value && foo.is_set (i))
    {
      cerr << "Yipes, didn't set bit " << i << " but it is set\n";
      failures++;
    }
  }

  foo.set_vector (ii);

  for (int i = 0; i < bsize; i++)
  {
    if (ii[i] != expected_result[i])
    {
      cerr << "Set vector failed at bit " << i << " expected " << expected_result[i] <<
              " got " << ii[i] << endl;
      failures++;
    }
  }

  return 0;
}

static int
test_set_vector ()
{
  int * i1 = nullptr;
  int * i2 = nullptr;

  int bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "Test set vector being done with vectors of size " << bsize << endl;

      if (nullptr != i1)
        delete i1;
      if (nullptr != i2)
        delete i2;

      i1 = new int[bsize];
      i2 = new int[bsize];
    }

    (void) test_set_vector (bsize, i1, i2);
  }

  if (nullptr != i1)
    delete i1;
  if (nullptr != i2)
    delete i2;

  return 0;
}

static int
do_one_set_of_operator_tests (int bsize)
{
  IW_Bits_Base foo (bsize);
  IW_Bits_Base bar (bsize);

  int new_value = intbtwij (0, 1);
  int zbit = intbtwij (0, bsize - 1);

  int changed = 0;
  if (foo.is_set (zbit) && 0 == new_value)
    changed = 1;
  else if (! foo.is_set (zbit) && 1 == new_value)
    changed = 1;

  foo.set (zbit, new_value);
  if (changed && foo == bar)
  {
    cerr << "Equality between unequal vectors failed\n";
    failures++;
  }

  if (0 == changed && foo == bar)
    ;
  else if (changed && foo != bar)
    ;
  else
  {
    cerr << "Inequality between unequal vectors not detected\n";
    failures++;
  }

  bar.set (zbit, new_value);
  if (foo == bar)
    ;
  else
  {
    cerr << "Equality between equal vectors not detected\n";
    failures++;
  }

  if (bar != foo)
  {
    cerr << "Inequality between equal vectors found\n";
    failures++;
  }

  if (bar.nset () != foo.nset ())
  {
    cerr << "Equal bit vectors should have equal nset () values\n";
    failures++;
  }

  IW_Bits_Base tmp = foo;
  tmp.iwand (bar);
  if (tmp != foo)
  {
    cerr << "And between two equal vectors should not change anything\n";
    failures++;
  }

  tmp.iwor (bar);
  if (tmp != foo)
  {
    cerr << "OR between two equal vectors should not change anything\n";
    failures++;
  }

  tmp.iwxor (foo);
  if (0 != tmp.nset ())
  {
    cerr << "XOR between two equal vectors should have no bits\n";
    failures++;
  }

  if ((foo & bar) != foo)
  {
    cerr << "And between two equal vectors should not change anything\n";
    failures++;
  }

/*if ((bar | foo) != foo)
  {
    cerr << "OR between two equal vectors should not change anything\n";
    failures++;
  }

  if (0 != (foo ^ bar).nset ())
  {
    cerr << "XOR between two equal vectors should have no bits\n";
    failures++;
  }*/

  set_random_bits (foo);
  set_random_bits (bar);

  IW_Bits_Base zresult (foo);

  zresult.iwor (bar);

  for (int i = 0; i < bsize; i++)
  {
    int fooi = foo.is_set (i);
    int bari = bar.is_set (i);
    int ri = zresult.is_set (i);

    if (0 == fooi && 0 == bari && 0 == ri)    // both unset
      ;
    else if (fooi && bari && ri)              // both set
      ;
    else if (fooi && ri && 0 == bari)         // foo set, bar unset
      ;
    else if (0 == fooi && ri && bari)         // foo unset, bar set
      ;
    else
    {
      cerr << "Bit " << i << " set or not set " << fooi << " vs " << bari << " in iwor operation\n";
      failures++;
    }
  }

  return 0;
}

static int
do_operator_tests ()
{
  int bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "do_operator_tests using vectors of length " << bsize << endl;
    }

    IW_Bits_Base foo (bsize);
    IW_Bits_Base bar (bsize);

    if (foo == bar)
      ;
    else
    {
      cerr << "Equality between newly created vectors failed\n";
      failures++;
    }

    if (foo != bar || bar != foo)
    {
      cerr << "Inequality between newly created vectors failed\n";
      failures++;
    }
  }

  for (int i = 1; i < 132; i++)
  {
    do_one_set_of_operator_tests (i);
  }

  bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "do_operator_tests using vectors of length " << bsize << endl;
    }

    do_one_set_of_operator_tests (bsize);
  }

// Test short vectors

  bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (8, 31);

      if (verbose)
        cerr << "Short vector tests done on bits of length " << bsize << endl;
    }

    do_one_set_of_operator_tests (bsize);
  }

  bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (1, 7);

      if (verbose)
        cerr << "Very short vector tests done on bits of length " << bsize << endl;
    }

    do_one_set_of_operator_tests (bsize);
  }

  return 0;
}

static int
test_flip (int n,
           int * ivec)
{
  IW_Bits_Base foo (n);

  for (int i = 0; i < n; i++)
  {
    ivec[i] = intbtwij (0, 1);
  }

  foo.construct_from_array_of_ints (ivec, n);

  foo.flip ();

  for (int i = 0; i < n; i++)
  {
    int fi = foo.is_set (i);

    if (fi && 0 == ivec[i])
      ;
    else if (! fi && ivec[i])
      ;
    else
    {
      cerr << "Flip failure, bit " << i << " in vector of size " << n << " bit vector " << fi << " array " << ivec[i] << endl;
      failures++;
    }
  }

  return 1;
}

static int
test_flip (int n)
{
  int * tmp = new int[n];

  int rc = test_flip (n, tmp);

  delete tmp;

  return rc;
}

static int
test_flip ()
{
  IW_Bits_Base foo (1);

  foo.flip ();

  if (! foo.is_set (0))
  {
    cerr << "Flip on 1 bit item failed\n";
    failures++;
  }

  foo.flip ();

  if (foo.is_set (0))
  {
    cerr << "Flip back on 1 bit item failed\n";
    failures++;
  }

  for (int i = 1; i < 129; i++)
  {
    test_flip (i);
  }

  return 1;
}

static int
test_turning_a_bit_on (int bsize)
{
  IW_Bits_Base foo (bsize);

  int zbit = intbtwij (0, foo.nbits () - 1);
  foo.set (zbit, 1);

  if (! foo.is_set (zbit))
  {
    cerr << "Cannot set bit " << zbit << endl;
    failures++;
  }

  IW_Bits_Base bar = foo;
  if (! bar.is_set (zbit))
  {
    cerr << "Bit " << zbit << " in copy found not set\n";
    failures++;
  }

  return 0;
}

static int
test_turning_a_bit_off (int bsize)
{
  IW_Bits_Base foo (bsize);

  int zbit = intbtwij (0, foo.nbits () - 1);
  foo.set (zbit, 0);

  if (foo.is_set (zbit))
  {
    cerr << "Cannot unset bit " << zbit << endl;
    failures++;
  }

  IW_Bits_Base bar = foo;
  if (bar.is_set (zbit))
  {
    cerr << "Bit in copy found set\n";
    failures++;
  }

  return 0;
}

static int
do_single_set_test_test (int bsize)
{
  int on_off = intbtwij (0, 1);
  if (on_off)
    return test_turning_a_bit_on (bsize);
  else
    return test_turning_a_bit_off (bsize);
}

static int
do_multiple_bit_setting (int bsize)
{
  IW_Bits_Base foo (bsize);

  int step = intbtwij (1, bsize - 1);
  for (int i = 0; i < bsize; i += step)
  {
    foo.set (i, 1);
  }

  int expected_nset;
  if (0 == foo.nbits () % step)
    expected_nset = foo.nbits () / step;
  else 
    expected_nset = foo.nbits () / step + 1;

//if (0 == expected_nset)
//  expected_nset = 1;

  if (expected_nset != foo.nset ())
  {
    cerr << "nbits = " << foo.nbits () << " step = " << step << " but nset = " << foo.nset () << endl;
    failures++;
  }

  foo.clear ();
  if (0 != foo.nset ())
  {
    cerr << "After clear, a vector should have no bits set\n";
    failures++;
  }

  return 0;
}

static int
do_set_test_tests ()
{
  int bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "do_set_test_tests using vectors of length " << bsize << endl;
    }

    do_single_set_test_test (bsize);
  }

  bsize = 0;
  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "do_set_test_tests using vectors of length " << bsize << endl;
    }

    do_multiple_bit_setting (bsize);
  }

  return 0;
}

static int
do_constructor_tests (int bsize)
{
  IW_Bits_Base foo (bsize);

  if (bsize != foo.nbits ())
  {
    cerr << "Newly constructed bit vector should have " << bsize << " bits\n";
    failures++;
  }

  IW_Bits_Base bar = foo;
  if (bar != foo)
  {
    cerr << "Inqeuality between two equal vectors found\n";
    failures++;
  }

  if (bar == foo && foo == bar)
    ;
  else
  {
    cerr << "Equality between equal vectors of length " << foo.nbits () << " failed\n";
    failures++;
  }

  return 1;
}

static int
do_constructor_tests ()
{
  int bsize = 0;

  for (int i = 0; i < iterations; i++)
  {
    if (0 == i % change_bit_size)
    {
      bsize = intbtwij (lower_vector_size, upper_vector_size);
      if (verbose)
        cerr << "do_constructor_tests using vectors of length " << bsize << endl;
    }

    do_constructor_tests (bsize);
  }

  return 1;
}

static void
do_print_bit_mask (ostream & os)
{
  for (int i = 0; i < IW_BITS_PER_WORD; i++)
  {
  }
}

static int
check_time (int (& f) (), const char * zname)
{
  time_t t0 = time (NULL);

  int rc = f ();

  if (verbose)
  {
    time_t t1 = time (NULL);
    cerr << zname << " took " << (t1 - t0) << " seconds\n";
  }

  return rc;
}

#include "cmdline.h"

static int
tiwbits (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:rb:B:i:phD:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('h'))
  {
    usage (2);
  }

  if (cl.option_present ('b'))
  {
    if (! cl.value ('b', lower_vector_size) || lower_vector_size <= 0)
    {
      cerr << "The -b option must be followed by a whole positive number\n";
      usage (1);
    }
  }

  if (cl.option_present ('B'))
  {
    if (! cl.value ('B', upper_vector_size) || upper_vector_size < lower_vector_size)
    {
      cerr << "The -b option must be followed by a whole number >= " << lower_vector_size << endl;
      usage (2);
    } 
  }

  if (verbose)
    cerr << "Vectors will be created with between " << lower_vector_size << " and " <<
            upper_vector_size << " bits\n";

  if (cl.option_present ('i'))
  {
    if (! cl.value ('i', iterations) || iterations < 1)
    {
      cerr << "The -i option requires a whole positive number\n";
      usage (3);
    }
    if (verbose)
      cerr << "Will perform " << iterations << " of all iterative tests\n";

    iterations_split = iterations / 10;

    if (0 == iterations_split)
      iterations_split = 1;

    fewer_iterations = iterations / 5;

    if (0 == fewer_iterations)
      fewer_iterations = 1;

    fewer_iterations_split = fewer_iterations / 10;
    if (0 == fewer_iterations_split)
      fewer_iterations_split = 1;
  }

  if (cl.option_present ('s'))
  {
    IWString tstring;
    cl.value ('s', tstring);
    int seed;

    if ('*' == tstring)
      iw_random_seed ();
    else if (! tstring.is_int (seed))
    {
      cerr << "Random number seeds must be whole numbers\n";
      usage (8);
    }
    else
      iw_set_rnum_seed (seed);
  }
  else
  {
    iw_random_seed ();
  }

  if (cl.option_present ('p'))
  {
    do_print_bit_mask (cout);
  }

  time_t t0 = time (NULL);

  check_time (do_constructor_tests, "do_constructor_tests");

  check_time (do_set_test_tests, "do_set_test_tests");

  check_time (do_operator_tests, "do_operator_tests");

  check_time (test_set_vector, "test_set_vector");

  check_time (test_increment_vector, "test_increment_vector");

  check_time (test_compute_weight, "test_compute_weight");

  check_time (test_construct_from_array_of_ints, "test_construct_from_array_of_ints");

  check_time (test_set_all, "test_set_all");

  check_time (test_first_bit, "test_first_bit");

  check_time (test_setting_ranges, "test_setting_ranges");

  check_time (test_shift, "test_shift");

  check_time (test_ascii_01_things, "test_ascii_01_things");

  check_time (test_printon, "test_printon");

  check_time (test_hex_stuff, "test_hex_stuff");

  check_time (test_sparse_things, "test_sparse_things");

  check_time (test_iterator, "test_iterator");

  check_time (test_flip, "test_flip");

  check_time (problem_apr03, "apr03");

  if (cl.option_present ('D'))
  {
    IWString d = cl.string_value ('D');

    test_reading_daylight_ascii_representation_file (d.null_terminated_chars ());
  }

  cerr << "Encountered " << failures << " failures\n";

  if (verbose)
  {
    time_t t1 = time (NULL);
    cerr << "All processing took " << (t1 - t0) << " seconds\n";
  }

  return failures;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tiwbits (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
