/*
  Test reading Daylight ascii representations

  Just a copy of nnetin

  Creates an input file for the neural set programme
  consisting of
    name  bits
*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <iostream>
#include <iomanip>

#include "iwbits.h"
#include "iwstring_data_source.h"
#include "accumulator.h"

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

using std::cerr;
using std::cout;
using std::ostream;
using std::setw;

const char * prog_name = NULL;

static int verbose = 0;

static int write_header = 0;

static int write_01_representation = 1;

/*
  In addition, we can specify the name stem used in the header
*/

static IWString bit_name ("bit");

static IWString fingerprint_tag ("FP<");

static int is_hex = 0;

static int append_nset = 0;

static int tdts_read = 0;
static int fingerprints_read = 0;

static Accumulator_Int<int> nbits;
static Accumulator_Int<int> nset;

static int fingerprints_with_zero_bits_set = 0;

static int
do_write_header (ostream & os, int nbits)
{
  os << "Name";

  for (int i = 1; i <= nbits; i++)
  {
    os << ' ' << bit_name << i;
  }

  if (append_nset)
    os << ' ' << bit_name << "nset";

  os << endl;

  return os.good ();
}

static void
usage (int rc)
{
  cerr << "Usage: " << prog_name << " <options> <tdtfile1>\n";
  cerr << "  -h                write header record\n";
  cerr << "  -s n              compute population counts for fingerprints <=n long\n";
  cerr << "  -N <name>         specify header record name for bits (" << bit_name << ")\n";
  cerr << "  -x                input is hex encoded\n";
  cerr << "  -F <tag>          fingeprint tag (default '" << fingerprint_tag << "')\n";
  cerr << "  -e                append the number of bits set to the descriptors\n";
  cerr << "  -v                verbose output\n";

  exit (rc);
}

#include "iwstring_data_source.h"

static int population_count_vector_length = 0;
static int * population = NULL;
static int fingerprints_summed = 0;

static int bits_in_fingerprints = 0;

static int
nnetin (iwstring_data_source & input, ostream & output)
{
  IWString buffer;
  IWString name;
  IW_Bits_Base * fp = NULL;

  while (input.next_record (buffer))
  {
    if (buffer.starts_with ("PCN<"))
    {
      name = buffer;
      assert (name.ends_with ('>'));
      name.chop ();
      name.remove_leading_chars (4);
      name.strip_trailing_blanks ();
      name.gsub (' ', '_');
      assert (name.length ());
      continue;
    }

    if (buffer.starts_with (fingerprint_tag))
    {
      fp = new IW_Bits_Base;
      if (is_hex)
      {
        const_IWSubstring h = buffer;
        h.remove_up_to_first ('<');
        h.chop ();

        if (! fp->construct_from_hex (h))
        {
          cerr << "Cannot parse hex encoded fingerprint\n";
          cerr << buffer << endl;
          return 0;
        }
      }
      else if (! fp->construct_from_tdt_record (buffer))
      {
        cerr << "Yipes, cannot contstruct fp from record " << input.lines_read () << endl;
        output << "barf\n";
      }

      if (0 == bits_in_fingerprints)
      {
        bits_in_fingerprints = fp->nbits ();
        if (write_header)
          do_write_header (output, bits_in_fingerprints);
      }
      else if (fp->nbits () != bits_in_fingerprints)
      {
        cerr << "Fingerprint bit count mismatch, previous = " << bits_in_fingerprints <<
                " current = " << fp->nbits () << endl;
        return 0;
      }

      fingerprints_read++;

      continue;
    }

    if (buffer == '|')
    {
      if (name.nchars () && NULL != fp)
      {
        if (write_01_representation)
        {
          IWString buffer (name);
          buffer.resize (buffer.nchars () + fp->nbits () * 2 + 9);   // 9 is arbitrary, just some slop

          fp->append_string_form (buffer, '1', '0', 1);

          if (append_nset)
            buffer += ' ' << fp->nset ();

          buffer += '\n';

          output << buffer;
        }

        if (0 == fp->nset ())
          fingerprints_with_zero_bits_set++;

        if (verbose)
        {
          nbits.extra (fp->nbits ());
          nset.extra (fp->nset ());
        }

        if (0 == population)
          ;
        else if (fp->nbits () >= population_count_vector_length)
          ;
        else
        {
          fp->increment_vector (population);
          fingerprints_summed++;
        }
      }
      else
      {
        if (name.nchars () || NULL != fp)
          cerr << "Name or fp missing\n";
      }
      
      if (NULL != fp)
      {
        delete fp;
        fp = NULL;
      }
      name = "";

      tdts_read++;
    }
  }

  return 1;
}

static int
nnetin (const char * fname, ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 1;
  }

  return nnetin (input, output);
}

#include "cmdline.h"

static int
nnetin (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vN:ns:hxF:e");

  if (cl.unrecognised_options_encountered ())
    usage (1);

  verbose = cl.option_count ('v');

  if (cl.option_present ('h'))
  {
    write_header = 1;
    if (verbose)
      cerr << "A header record will be written\n";
  }

  int tmp;
  if (cl.option_present ('s') && cl.value ('s', tmp) && tmp > 0)
  {
    population_count_vector_length = tmp;
    population = new int[population_count_vector_length];

    for (int i = 0; i < population_count_vector_length; i++)
    {
      population[i] = 0;
    }

    if (verbose)
      cerr << "Will compute bit populations for bits <= " << population_count_vector_length << endl;
  }

  if (cl.option_present ('N'))
  {
    cl.value ('N', bit_name);
    if (verbose)
      cerr << "Column names prefaced by '" << bit_name << "'\n";
  }

  if (cl.option_present ('n'))
  {
    write_01_representation = 0;

    if (verbose)
      cerr << "Normal output suppressed\n";
  }

  if (cl.option_present ('F'))
  {
    fingerprint_tag = cl.string_value ('F');

    if (verbose)
      cerr << "Fingerprint tag '" << fingerprint_tag << "'\n";

    if (! fingerprint_tag.ends_with ('<'))
      fingerprint_tag += '<';
  }

  if (cl.option_present ('x'))
  {
    is_hex = 1;

    if (verbose)
      cerr << "Input will be treated as hex encoded\n";
  }

  if (cl.option_present ('e'))
  {
    append_nset = 1;

    if (verbose)
      cerr << "Will include the number of bits set in the descriptors\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "No input files specified\n";
    usage (8);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! nnetin (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << tdts_read << " tdts with ";
    if (fingerprints_read != tdts_read)
      cerr << fingerprints_read << ' ';
    cerr << "fingerprints\n";
  }

  if (0 == rc && population)
  {
    cerr << "Summed " << fingerprints_summed << " fingerprints\n";
    int zero_count = 0;
    int all_count = 0;
    for (int i = 0; i < population_count_vector_length; i++)
    {
      cerr << "Bit " << i << " population " << population[i];
      float tmp = float (population[i]) / float (fingerprints_summed);
      cerr << " " << setw (3) << tmp;

      if (0 == population[i])
      {
        zero_count++;
        cerr << "    -- 0";
      }
      else if (fingerprints_summed == population[i])
      {
        all_count++;
        cerr << "    -- *";
      }

      cerr << endl;
    }

    if (zero_count)
      cerr << zero_count << " bits had no hits\n";
    if (all_count)
      cerr << all_count << " bits were hit for each molecule\n";

    delete population;
  }

  if (verbose)
  {
    cerr << "Read " << nbits.n () << " fingerprints, had between " << nbits.minval () << " and " << nbits.maxval () << " bits";
    if (nbits.n () > 1 && nbits.minval () != nbits.maxval ())
      cerr << " ave " << nbits.average ();
    cerr << endl;
    if (fingerprints_with_zero_bits_set)
      cerr << fingerprints_with_zero_bits_set << " fingerprints had zero bits set\n";
    cerr << "Fingerprints had between " << nset.minval () << " and " << nset.maxval () << " bits set";
    if (nset.n () > 1)
      cerr << " ave " << nset.average ();
    cerr << endl;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = nnetin (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
