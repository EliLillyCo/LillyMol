#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

/*
  The basic problem is to take two arrays of unsigned int's
  and compute the number of bits which are set in both
  arrays.
*/

#define NWORDS 64

/*
  bits_in_common is the critical routine
*/

extern "C" int bits_in_common (const unsigned int *, const unsigned int *, int);

int
main ()
{
  unsigned int v1[NWORDS], v2[NWORDS];     // a couple of vectors

// for the first test, turn off all bits

  for (int i = 0; i < NWORDS; i++)
  {
    v1[i] = v2[i] = 0;
  }

  int bic = bits_in_common (v1, v2, NWORDS);

// should be no bits in common

  assert (0 == bic);

// for the second test, turn on all the bits

  for (int i = 0; i < NWORDS; i++)
  {
    v1[i] = v2[i] = 0xffffffff;
  }

  bic = bits_in_common (v1, v2, NWORDS);

  assert (32 * NWORDS == bic);     // all bits in common

// and an intermediate case - just selected at random

  v1[0] = v2[0] = 0xf0f0f0f0;

  bic = bits_in_common (v1, v2, 1);

  assert (16 == bic);

// For timing, fill the arrays with arbitrary numbers

  unsigned int pid = getpid ();
  for (int i = 0; i < NWORDS; i++)
  {
    v1[i] = 17927 * i * pid;
    v2[i] = 37907 * i * pid;
  }

// Now for timing. Adjust ITERATIONS so that timing can be done

  int iterations = 1000000;

  for (int i = 0; i < iterations; i++)
  {
    (void) bits_in_common (v1, v2, NWORDS);
  }

// Using our existing algorithm, this takes 5 seconds on an Ultra 1 (166 MHz)
// and about 3 seconds on the 336 MHz server.

  return 0;
}
