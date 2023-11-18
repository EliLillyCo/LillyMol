/*
  Computes distances to each of a set of fingerprints
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString identifier_tag("PCN<");

static IW_General_Fingerprint * pool = nullptr;
static int pool_size = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes distances from a set of fingerprints and produces descriptors\n";
  cerr << " -p <fname>     name of file against which to compare\n";
  cerr << " -F,-P,...      standard fingerprint options\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
build_pool (iwstring_data_source & input)
{
  assert (pool_size > 0);

  int i = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    if (! pool[i].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    i++;

    if (i >= pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }
  }

  if (verbose)
  {
    cerr << i << " fingerprint objects added to pool\n";
  }

  pool_size = i;

  return 1;
}

static int
build_pool (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)     // must grep the file to find out how many
  {
    assert (NULL == pool);

    IWString tmp = '^';
    tmp << identifier_tag;

    pool_size = input.grep (tmp);
    if (0 == pool_size)
    {
      cerr << "Yipes, cannot find any '" << tmp << "' in the input\n";
      return 0;
    }

    if (verbose)
      cerr << "Input contains " << pool_size << " fingerprints\n";

    pool = new IW_General_Fingerprint[pool_size];

    if (NULL == pool)
    {
      cerr << "Yipes, cannot allocate space for " << pool_size << " fingerprints\n";
      return 0;
    }
  }

  return build_pool (input);
}

static int
gfp_proximity_descriptors (IW_General_Fingerprint & fp,
                           IWString_and_File_Descriptor & output)
{
  append_first_token_of_name(fp.id(), output);

  for (int i = 0; i < pool_size; i++)
  {
    similarity_type_t d = pool[i].distance(fp);

    output << ' ' << d;
  }

  output << '\n';

  return 1;
}

static int
gfp_proximity_descriptors (IW_TDT & tdt,
                           IWString_and_File_Descriptor & output)
{
  int fatal;
  IW_General_Fingerprint fp;

  if (! fp.construct_from_tdt (tdt, fatal))
  {
    if (fatal)
    {
      cerr << "Cannot parse tdt " << tdt;
      return 0;
    }

    return 1;
  }

  return gfp_proximity_descriptors (fp, output);
}

static int
gfp_proximity_descriptors (iwstring_data_source & input,
                           IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    IW_General_Fingerprint fp;

    if (! gfp_proximity_descriptors (tdt, output))
    {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

static int
gfp_proximity_descriptors (const char * fname,
     IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_proximity_descriptors(input, output);
}

static int
write_header (const IW_General_Fingerprint * pool,
              int pool_size,
              IWString_and_File_Descriptor & output)
{
  output << "Name";

  for (int i = 0; i < pool_size; i++)
  {
    output << " D_";
    append_first_token_of_name(pool[i].id(), output);
  }

  output << '\n';

  return 1;
}

static int
gfp_proximity_descriptors (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:P:V:p:s:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage (17);
    }
  }
  else if (! initialise_fingerprints (cl[0], verbose))
  {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }

    pool = new IW_General_Fingerprint[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose)
      cerr << "system sized to " << pool_size << endl;
  }

  if (! cl.option_present('p'))
  {
    cerr << "Must specify pool file via the -p option\n";
    usage(3);
  }

  if (cl.option_present('p'))
  {
    const_IWSubstring p = cl.string_value('p');

    if (! build_pool(p) || 0 == pool_size)
    {
      cerr << "Could not build pool from '" << p << "'\n";
      return 9;
    }

    if (verbose)
      cerr << "Read " << pool_size << " fingerprints\n";
  }

  IWString_and_File_Descriptor output(1);

  write_header (pool, pool_size, output);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! gfp_proximity_descriptors(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_proximity_descriptors(argc, argv);

  return rc;
}
