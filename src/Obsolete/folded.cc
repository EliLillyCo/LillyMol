/*
  When using hashed fingerprints, we often need to know how badly
  collisions are affecting us
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/misc.h"

#include "Utilities/GFP_Tools/dyfp.h"

using std::cout;
using std::ostream;
using std::cerr;
using std::endl;

static int verbose = 0;

static int fingerprints_read = 0;

/*
  We read fingerprints of a given size and determine how bad the collisions
  are when folded down to FINAL_NBITS
*/

static int final_nbits = 0;

static int * bit_ever_set = nullptr;

static int * collisions_on_bit = nullptr;

/*
  The number of bits in the fingerprints as they are read in
*/

static int initial_nbits = 0;

static IWString identifier_tag ("PCN<");

static IWString fingerprint_tag;

static int fingerprints_with_colliding_bits = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Examines collisions in a fingerprint file\n";
  cerr << " -F <tag>         fingerprint to process\n";
  cerr << " -I <tag>         idenifier tag (default " << identifier_tag << ")\n";
  cerr << " -b <nbits>       final folded size - must be power of 2\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
folded (const IWString & id,
        IWDYFP & fp,
        ostream & output)
{
  assert (initial_nbits == fp.nbits ());

  int multiplier = initial_nbits / final_nbits;

  int final_nset = 0;
  int colliding_bits = 0;

  for (int i = 0; i < final_nbits; i++)
  {
    int bits_folding_to_this_bit = 0;
    for (int j = 0; j < multiplier; j++)
    {
      int b = j * multiplier + i;
      if (fp.is_set (b))
        bits_folding_to_this_bit++;
    }

    if (bits_folding_to_this_bit)
    {
      final_nset++;
      bit_ever_set[i]++;
    }

    if (bits_folding_to_this_bit > 1)
    {
      colliding_bits++;

      collisions_on_bit[i]++;

      if (verbose > 2)
        cerr << " bit " << i << " folded " << bits_folding_to_this_bit << endl;
    }
  }

  if (colliding_bits)
    fingerprints_with_colliding_bits++;

  if (verbose > 1)
    cerr << id << ' ' << colliding_bits << " of " << final_nset << " bits colliding\n";

  return output.good ();
}

/*
  final_nbits must be a multiple of initial_nbits
*/

static int
check_proper_relationship_between_initial_and_final_nbits ()
{
  assert (initial_nbits > 0);
  assert (final_nbits > 0);

  if (initial_nbits == (initial_nbits / final_nbits) * final_nbits)
    return 1;

  cerr << "INcorrect relationship between initial nbits " << initial_nbits << " and final " << final_nbits << endl;

  return 0;
}

static int
folded (iwstring_data_source & input, ostream & output)
{
  const_IWSubstring buffer;
  IWString id;
  IWDYFP fp;
  int fingerprint_found = 0;

  while (input.next_record (buffer))
  {
    if (buffer.starts_with (identifier_tag))
    {
      id = buffer;
      assert (id.ends_with ('>'));
      id.chop ();
      id.remove_leading_chars (identifier_tag.length ());
    }
    else if (buffer.starts_with (fingerprint_tag))
    {
      if (! fp.construct_from_tdt_record (buffer))
      {
        cerr << "Cannot parse fingerprint\n";
        cerr << buffer << endl;
        return 0;
      }
      fingerprint_found = 1;

      if (0 == initial_nbits)
      {
        initial_nbits = fp.nbits ();
        check_proper_relationship_between_initial_and_final_nbits ();

        if (verbose)
          cerr << "Initial fingerprints contain " << initial_nbits << " bits\n";
      }
      else if (fp.nbits () != initial_nbits)
      {
        cerr << "Bit count mismatch, expected " << initial_nbits << " got " << fp.nbits () << endl;
        return 0;
      }
    }
    else if ('|' == buffer)
    {
      if (0 == id.length () || ! fingerprint_found)
      {
        cerr << "Very strange, end of TDT but either no id '" << id << "' or fingerprint\n";
        fingerprint_found = 0;
        id.resize_keep_storage (0);
      }
      else
      {
        fingerprints_read++;

        folded (id, fp, output);
        id.resize_keep_storage (0);
        fingerprint_found = 0;
      }
    }
  }

  return output.good ();
}

static int
folded (const char * fname, ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return folded (input, output);
}

static int
folded (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vI:F:b:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present ('F'))
  {
    cerr << "Must specify fingerprint tag via the -F option\n";
    usage (11);
  }

  if (cl.option_present ('F'))
  {
    fingerprint_tag = cl.string_value ('F');

    if (verbose)
      cerr << "Will examine fingerprint '" << fingerprint_tag << "'\n";

    if (! fingerprint_tag.ends_with ('<'))
      fingerprint_tag += '<';
  }

  if (cl.option_present ('I'))
  {
    identifier_tag = cl.string_value ('I');

    if (verbose)
      cerr << "Identifiers in '" << identifier_tag << "'\n";

    if (! identifier_tag.ends_with ('<'))
      identifier_tag += '<';
  }

  if (! cl.option_present ('b'))
  {
    cerr << " Must specify final number of bits via the -b option\n";
    usage (3);
  }

  if (cl.option_present ('b'))
  {
    if (! cl.value ('b', final_nbits) || final_nbits < 1)
    {
      cerr << "The final number of bits option (-b) must be a whole positive number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will report collisions if fingerprints folded to " << final_nbits << " bits\n";
  }

  collisions_on_bit = new_int (final_nbits);
  bit_ever_set = new_int (final_nbits);

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! folded (cl[0], cout))
    {
      rc = i + 1;
      break;
    }
  }

  cerr << "Read " << fingerprints_read << " fingerprints\n";
  cerr << fingerprints_with_colliding_bits << " fingerprints had one or more bit collisions\n";

  if (verbose)
  {
    int bits_ever_set = 0;
    for (int i = 0; i < final_nbits; i++)
    {
      if (bit_ever_set[i])
        bits_ever_set++;
    }

    cerr << bits_ever_set << " of " << final_nbits << " bits were ever set\n";

    int bits_with_collisions = 0;
    for (int i = 0; i < final_nbits; i++)
    {
      if (collisions_on_bit[i])
        bits_with_collisions++;

      if (verbose > 1 && collisions_on_bit[i])
      {
        cerr << "Bit " << (i + 1) << " " << collisions_on_bit[i] << " collisions\n";
      }
    }

    cerr << bits_with_collisions << " of " << final_nbits << " bits had collisions\n";
  }

  if (NULL != collisions_on_bit)
    delete collisions_on_bit;

  if (NULL != bit_ever_set)
    delete bit_ever_set;

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = folded (argc, argv);

  return rc;
}
