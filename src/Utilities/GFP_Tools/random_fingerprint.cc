/*
  We want to generate random fingerprints - for testing
*/

#include <stdlib.h>
#include <iostream>
#include <random>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/accumulator/accumulator.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString fingerprint_tag ("FPRAND<");

static int work_as_filter = 0;

static int nbits = 256;

static Accumulator_Int<int> acc_nset;

static int required_nset = 0;

static int bits_to_change = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Generates random fingerprints\n";
  cerr << " -F <tag>       tag for fingerprints (default '" << fingerprint_tag << ")\n";
  cerr << " -n <nbits>     number of bits (default " << nbits << ")\n";
  cerr << " -f             work as a filter, reading from stdin\n";
  cerr << " -d <fraction>  desired fraction bits set in each fp\n";
  cerr << " -c <nbits>     number of bits to change on each new fingerprint\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
set_this_many_bits (IW_Bits_Base & fp,
                    int bits_to_set)
{
  int nset = 0;


  std::random_device rd;
  std::default_random_engine rng;
  std::uniform_int_distribution<off_t> u(0, fp.nbits() - 1);
  while (nset < bits_to_set)
  {
    const int i = u(rng);

    if (fp.is_set(i))
      continue;

    fp.set(i);
    nset++;
  }

  return nset;
}

static int
random_fingerprint_density_target (ostream & output)
{
  IW_Bits_Base fp (nbits);

  set_this_many_bits (fp, required_nset);

  return fp.write_daylight_ascii_representation (output, fingerprint_tag);
}

static int
random_fingerprint_small_change (ostream & output)
{
  static IW_Bits_Base fp (nbits);

  if (0 == fp.nset ())
    set_this_many_bits (fp, required_nset);
  else
  {
    std::random_device rd;
    std::default_random_engine rng;
    std::uniform_int_distribution<off_t> u(0, nbits - 1);
    for (int i = 0; i < bits_to_change; i++)
    {
      int j = u(rng);
  
      if (fp.is_set(j))
        fp.set(j, 0);
      else
        fp.set(j, 1);
    }
  }

  fp.write_daylight_ascii_representation (output, fingerprint_tag);

  return output.good ();
}

static int
random_fingerprint(ostream & output)
{
  if (bits_to_change)
    return random_fingerprint_small_change(output);

  if (required_nset)
    return random_fingerprint_density_target(output);

  IW_Bits_Base fp(nbits);

  std::random_device rd;
  std::default_random_engine rng;
  std::bernoulli_distribution d(0.50);
  for (int i = 0; i < nbits; i++)
  {
    if (d(rng))
      fp.set(i);
  }

  acc_nset.extra(fp.nset());

  return fp.write_daylight_ascii_representation(output, fingerprint_tag);
}

static int
random_fingerprint (iwstring_data_source & input,
                    ostream & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (work_as_filter)
    {
      if ("|" == buffer)
        random_fingerprint (output);

      output << buffer << endl;
    }
    else
    {
      const_IWSubstring smiles, id;
      if (! buffer.split (smiles, ' ', id))
      {
        cerr << "Cannot split '" << buffer << "'\n";
        return 0;
      }

      output << smiles_tag << smiles << ">\n";
      output << identifier_tag << id << ">\n";

      if (! random_fingerprint (output))
        return 0;

      output << "|\n";
    }
  }

  return output.good ();
}

static int
random_fingerprint (const const_IWSubstring & fname,
                    ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return random_fingerprint (input, output);
}

static int
random_fingerprint (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:n:fd:c:i:");   // the -i option does nothing - gfp_make insists on passing -i smi

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

  if (cl.option_present ('F'))
  {
    cl.value ('F', fingerprint_tag);

    if (verbose)
      cerr << "Fingerprints generated with tag " << fingerprint_tag << "'\n";

    if (! fingerprint_tag.ends_with ('<'))
      fingerprint_tag.add ('<');
  }

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', nbits) || nbits < 1)
    {
      cerr << "Invalid number of bits\n";
      usage (5);
    }

    if (verbose)
      cerr << "Fingerprints will contain " << nbits << " bits\n";
  }

  if (cl.option_present ('d'))
  {
    float required_density;
    if (! cl.value ('d', required_density) || required_density <= 0.0 || required_density >= 1.0)
    {
      cerr << "Invalid density (-d option)\n";
      usage (5);
    }

    required_nset = static_cast<int> (nbits * required_density);

    if (0 == required_nset)
      required_nset = 1;
    else if (required_nset == nbits)
      required_nset = nbits - 1;

    if (verbose)
      cerr << "Bit density of " << required_density << " requires " << required_nset << " bits\n";
  }

  if (cl.option_present ('c'))
  {
    if (! cl.value ('c', bits_to_change) || bits_to_change < 1)
    {
      cerr << "The number of bits to change option must be a +ve whole number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will change " << bits_to_change << " bits on each iteration\n";

    if (0 == required_nset)
      required_nset = nbits / 4;
  }

  if (cl.option_present ('f'))
  {
    work_as_filter = 1;

    if (cl.number_elements () > 1)
    {
      cerr << "Working as filter, but more than 1 command line argument\n";
      return 5;
    }
  }


  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! random_fingerprint (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Generated " << acc_nset.n () << " fingerprints";
    if (acc_nset.n () > 1)
      cerr << ", nset between " << acc_nset.minval () << " and " << acc_nset.maxval () << " ave " << acc_nset.average_if_available_minval_if_not ();

    cerr << endl;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = random_fingerprint (argc, argv);

  return rc;
}
