/*
  
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iwbits.h"

using std::cout;
using std::ostream;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString fingerprint_tag ("FP<");

static int fingerprints_read = 0;

static int nbits = 0;

static int nbytes = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Converts fingerprints to binary data\n";
  cerr << " -F <tag>       fingerprint to process, default '" << fingerprint_tag << "'\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
tomdoman_fingerprints_to_bits (const const_IWSubstring & buffer,
                               ostream & output)
{
  IW_Bits_Base fp;

  if (! fp.construct_from_tdt_record (buffer))
    return 0;

  fingerprints_read++;

  if (0 == nbits)
  {
    nbits = fp.nbits ();

    if (nbits > 2048)
    {
      cerr << "Too many bits " << nbits << " to store in a byte\n";
      return 0;
    }

    if (nbits != nbits / IW_BITS_PER_BYTE * IW_BITS_PER_BYTE)
    {
      cerr << "Nbits " << nbits << " not a multiple of " << IW_BITS_PER_BYTE << endl;
      return 0;
    }

    nbytes = nbits / IW_BITS_PER_BYTE;

    unsigned char tmp = static_cast<unsigned char> (nbytes);
    output.write (reinterpret_cast<const char *> (&tmp), 1);

    if (verbose)
      cerr << "Fingerprints contain " << nbits << " bits, " << nbytes << " bytes\n";
  }
  else if (nbits != fp.nbits ())
  {
    cerr << "Bit count mismatch, expected " << nbits << " got " << fp.nbits () << endl;
    return 0;
  }

  output.write (reinterpret_cast<const char *> (fp.bits ()), nbytes);

  return output.good ();
}

static int
tomdoman_fingerprints_to_bits (iwstring_data_source & input,
                               ostream & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! buffer.starts_with (fingerprint_tag))
      continue;

    if (! tomdoman_fingerprints_to_bits (buffer, output))
    {
      cerr << "Fatal error processing line '" << input.lines_read () << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return output.good ();
}

static int
tomdoman_fingerprints_to_bits (const char * fname,
                               ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tomdoman_fingerprints_to_bits (input, output);
}

static int
tomdoman_fingerprints_to_bits (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('F'))
  {
    cl.value ('F', fingerprint_tag);

    if (! fingerprint_tag.ends_with ('<'))
      fingerprint_tag.add ('<');

    if (verbose)
      cerr << "Will process fingerprint '" << fingerprint_tag << "'\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! tomdoman_fingerprints_to_bits (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << fingerprints_read << " fingerprints\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tomdoman_fingerprints_to_bits (argc, argv);

  return rc;
}
