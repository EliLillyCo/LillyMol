/*
  Mostly a timing tester for bits in common
*/

#include <stdlib.h>

#include <values.h>
#ifdef linux
#else
#include <limits.h>
#endif

#include <iomanip>

#include "cmdline.h"
#include "iwbits.h"
#include "iwrandom.h"

static int verbose = 0;

/*
  To replicate bit densities encountered in real life, we may want to 
  increase or lower the bits set randomly
*/

static float desired_density = 0.0;

static int default_nbits = 512;

static int default_iterations = 2000;

static void
usage (int rc)
{
  cerr << "Tester for bits in common\n";
  cerr << " -n <iterations>  perform the tests <iterations> times (default " << default_iterations << ")\n";
  cerr << " -b <nbits>       size of bit vectors to use (default " << default_nbits << ")\n";
  cerr << " -d <density>     desired bit density\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

/*
  We have adjusted the density of a bit vector to the desired range. Copy
  the bits to a vector.
*/

static void
copy_bits (IW_Bits_Base & b, unsigned char * r)
{
  const unsigned char * x = reinterpret_cast<const unsigned char *> (b.bits ());

  int nb = b.nbits () / IW_BITS_PER_BYTE;

  for (int i = 0; i < nb; i++)
  {
    r[i] = x[i];
  }

  return;
}

static void
adjust_density_upwards (IW_Bits_Base & b, unsigned char * r)
{
  int nset = b.nset ();

  if (nset == b.nbits ())
  {
    cerr << "Cannot adjust density upwards, all bits set\n";
    return;
  }

  while (1)
  {
    int i = intbtwij (0, b.nbits () - 1);
    if (b.is_set (i))
      continue;

    b.set (i, 1);
    nset++;

    float density = static_cast<float> (nset) / static_cast<float> (b.nbits ());

    if (density < desired_density)
      continue;

    assert (nset == b.nset ());

    copy_bits (b, r);

    return;
  }
}

static void
adjust_density_downwards (IW_Bits_Base & b, unsigned char * r)
{
  int nset = b.nset ();

  if (0 == nset)
  {
    cerr << "Cannot adjust density downwards, all bits unset\n";
    return;
  }

  while (1)
  {
    int i = intbtwij (0, b.nbits () - 1);
    if (! b.is_set (i))
      continue;

    b.set (i, 0);
    nset--;

    float density = static_cast<float> (nset) / static_cast<float> (b.nbits ());

    if (density > desired_density)
      continue;

    assert (nset == b.nset ());

    copy_bits (b, r);

    return;
  }
}

static void
adjust_density (IW_Bits_Base & b, unsigned char * r)
{
  if (0.0 == desired_density)
    return;

  float density = static_cast<float> (b.nset ()) / static_cast<float> (b.nbits ());

  if (verbose > 1)
    cerr << "Initial nset " << b.nset () << " density " << density << endl;

  if (density == desired_density)
    return;

  if (density < desired_density)
    adjust_density_upwards (b, r);
  else
    adjust_density_downwards (b, r);

  return;
}

static int
tbic (int iterations, int nbits,
      int s,                   // the number of double's in r1 and r2
      random_number_t * r1,
      random_number_t * r2)
{
  iw_random_seed ();

  IW_Bits_Base b1 (nbits), b2 (nbits);

  for (int i = 0; i < iterations; i++)
  {
    for (int j = 0; j < s; j++)
    {
      r1[j] = DBL_MAX * iwrandom ();
      r2[j] = DBL_MAX * iwrandom ();
      if (iwrandom () < 0.5)
        r2[j] = - r2[j];
//    cerr << "j = " << j << " r1 = " << r1[j] << " r2 = " << r2[j] << endl;
    }

    b1.construct_from_array_of_bits (reinterpret_cast<const unsigned char *> (r1), nbits);
    adjust_density (b1, reinterpret_cast<unsigned char *> (r1));

    b2.construct_from_array_of_bits (reinterpret_cast<const unsigned char *> (r2), nbits);
    adjust_density (b2, reinterpret_cast<unsigned char *> (r2));

    for (int j = 0; j < 100; j++)
    {
      if (b1.bits_in_common (b2) != b2.bits_in_common (b1))
      {
        cerr << "Yipes, asymmetry, " << b1.bits_in_common (b2) << " vs " << b2.bits_in_common (b1) << endl;
        return 0;
      }

      if (b1.bits_in_common (b1) != b1.nset ())
      {
        cerr << "Yipes, nset error, bits in common self " << b1.bits_in_common (b1) << " vs nset " << b1.nset () << ", nbits " << b1.nbits() << endl;
        return 0;
      }

//    Shuffle r2

      random_number_t rsave = r2[0];
      for (int k = 0; k < s - 1; k++)
      {
        r2[k] = r2[k + 1];
      }
      r1[s - 1] = rsave;

      b2.construct_from_array_of_bits (reinterpret_cast<const unsigned char *> (r2), nbits);
    }
  }

  return 1;
}

/*
  Allocate two arrays for holding random numbers
*/

static int
tbic (int iterations, int nbits)
{
  int s = (nbits + IW_BITS_PER_BYTE) / IW_BITS_PER_BYTE / sizeof (random_number_t);
  random_number_t * r1 = new random_number_t[s];
  random_number_t * r2 = new random_number_t[s];

  int rc = tbic (iterations, nbits, s, r1, r2);

  delete r1;
  delete r2;

  return rc;
}

int
tbic (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:b:d:h");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('h'))
  {
    usage (3);
  }

  int iterations = default_iterations;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', iterations) || iterations < 1)
    {
      cerr << "The iteration count (-n) option must be followed by a whole positive number\n";
      usage (3);
    }

    if (verbose)
      cerr << "Will do " << iterations << " iterations\n";
  }
  else if (verbose)
    cerr << "Default iteration count " << iterations << endl;

  int nbits = default_nbits;

  if (cl.option_present ('b'))
  {
    if (! cl.value ('b', nbits) || nbits < 1)
    {
      cerr << "The number of bits (-b) option must be followed by a whole positive number\n";
      usage (5);
    }
    if (verbose)
      cerr << "Will test with bit vectors with " << nbits << " bits\n";
  }
  else if (verbose)
    cerr << "Default nbits " << nbits << endl;

  if (cl.option_present ('d'))
  {
    if (! cl.value ('d', desired_density) || desired_density <= 0.0 || desired_density >= 1.0)
    {
      cerr << "The desired density (-d) option must be followed by a valid density\n";
      usage (7);
    }

    if (verbose)
      cerr << "Desired bit density " << desired_density << endl;
  }

  (void) tbic (iterations, nbits);

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = tbic (argc, argv);

  return rc;
}
