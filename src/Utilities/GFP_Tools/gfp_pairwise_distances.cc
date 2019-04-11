/*
  Reads a set of fingerprints, and reports any pair-wise
  distances requested
*/

#include <stdlib.h>
#include <limits>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_tdt.h"
#include "iw_stl_hash_map.h"
#include "accumulator.h"
#include "iwhistogram.h"
#include "iwrandom.h"

#include "gfp.h"
#include "tversky.h"

const char * prog_name = NULL;

static int verbose = 0;

static IW_General_Fingerprint * pool = NULL;
static int pool_size = 0;

static Tversky tversky;

static IW_STL_Hash_Map_int id_to_ndx;

static int ignore_duplicate_identifiers_in_pool = 0;

static int duplicate_identifiers_encountered = 0;

static Accumulator<similarity_type_t> distance_acc;

static IWHistogram hist;

static float upper_distance_threshold = std::numeric_limits<float>::max();

static int discarded_by_upper_distance_threshold = 0;

static int write_entire_input_records = 0;

static int equal_weight_tanimoto = 0;

static int strip_leading_zeros = 0;

static int ignore_missing_identifiers = 0;

static int missing_identifiers_ignored = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Reports pair-wise distances. Input file has pairs of identifiers on each line\n";
  cerr << " -p <fname>     file of fingerprints - the haystack\n";
  cerr << " -g             ignore duplicate identifiers in the -p file\n";
  cerr << " -K id1,id2     specify identifiers to be processed on the command line\n";
  cerr << " -s <num>       number of fingerprints in -p file (optional)\n";
  cerr << " -T <dist>      upper distance threshold - don't write\n";
  cerr << " -o <ndigits>   output precision\n";
  cerr << " -V ...         Tversky options\n";
  cerr << " -F,-P, ...     standard gfp fingerprint options\n";
  cerr << " -a             write entire input records (not just first two tokens)\n";
  cerr << " -q             use equal weight tanimoto function on composite fingerprints\n";
  cerr << " -r <n>         report <n> randomly selected pairwise distances\n";
  cerr << " -b             ignore identifiers not present in the haystack file\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
read_pool (iwstring_data_source & input)
{
  int items_read = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;

    if (! pool[items_read].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Fatal error reading TDT, line " << input.lines_read () << endl;
        return 0;
      }

      continue;
    }

    const IWString & id = pool[items_read].id();

//  if (strip_leading_zeros)    implement this sometime
//  {
//  }

    if (id_to_ndx.contains(id))
    {
      cerr << "Duplicate id '" << id << "' cannot process\n";

      duplicate_identifiers_encountered++;

      if (ignore_duplicate_identifiers_in_pool)
        continue;
      else
        return 0;
    }

    id_to_ndx[id] = items_read;

    items_read++;
    if (items_read >= pool_size)
    {
      cerr << "Pool is full " << items_read << endl;
      return 1;
    }
  }

  pool_size = items_read;

  if (verbose)
    cerr << "Read " << pool_size << " items\n";

  if (static_cast<unsigned int>(pool_size) != id_to_ndx.size())
  {
    cerr << "Huh, must be duplicate ID's, " << pool_size << " fingerprints, but " << id_to_ndx.size() << " items in ID hash\n";
    return 1;
  }

  return 1;
}

static int
read_pool (const char * fname)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size)
  {
    int n = input.count_records_starting_with("|");
    if (0 == n)
    {
      cerr << "Yipes, no TDT's in the pool file\n";
      return 0;
    }

    if (n < 2)
    {
      cerr << "Not enough fingerprints in file\n";
      return 0;
    }

    pool = new IW_General_Fingerprint[n];
    if (NULL == pool)
    {
      cerr << "Sorry, cannot allocate " << n << " fingerprints\n";
      return 0;
    }

    pool_size = n;

    if (verbose)
      cerr << "Pool sized to " << n << " molecules\n";
  }

  return read_pool(input);
}

static similarity_type_t
compute_the_distance (int ndx1, int ndx2,
                      Tversky & tversky)
{
  if (tversky.active ())
    return static_cast<similarity_type_t>(1.0) - pool[ndx1].tversky(pool[ndx2], tversky);
  else if (equal_weight_tanimoto)
    return static_cast<similarity_type_t>(1.0) - pool[ndx1].equal_weight_tanimoto(pool[ndx2]);
  else
    return static_cast<similarity_type_t>(1.0) - pool[ndx1].tanimoto(pool[ndx2]);
}

static int
gfp_pairwise_distance(const const_IWSubstring & buffer,
                      IWString_and_File_Descriptor & output)
{
  if (buffer.nwords() < 2)
  {
    cerr << "Input records must contain at least two tokens\n";
    return 0;
  }

  IWString id1;
  int i = 0;

  buffer.nextword(id1, i);

  IW_STL_Hash_Map_int::const_iterator f = id_to_ndx.find(id1);

  if (f == id_to_ndx.end())
  {
    cerr << "No ID '" << id1 << "' in ID hash\n";
    return 0;
  }

  int ndx1 = (*f).second;

// And the second token

  IWString id2;

  buffer.nextword(id2, i);

  f = id_to_ndx.find(id2);

  if (f != id_to_ndx.end())
    ;
  else if (ignore_missing_identifiers)
  {
    missing_identifiers_ignored++;
    if (verbose)
      cerr << "No ID '" << id2 << "' in ID hash\n";

    return 1;
  }
  else
  {
    cerr << "No ID '" << id2 << "' in ID hash\n";
    return 0;
  }

  int ndx2 = (*f).second;

  similarity_type_t d = compute_the_distance (ndx1, ndx2, tversky);

  if (d > upper_distance_threshold)
  {
    discarded_by_upper_distance_threshold++;
    return 1;
  }

  output << id1 << ' ' << id2 << ' ' << d;

  if (write_entire_input_records)
  {
    const_IWSubstring token;

    while (buffer.nextword(token, i))
    {
      output << ' ' << token;
    }
  }

  output << '\n';

  if (verbose)
    distance_acc.extra(d);

  return 1;
}

static int
gfp_pairwise_distances (iwstring_data_source & input,
                        IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! gfp_pairwise_distance(buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int
gfp_pairwise_distances (const char * fname,
                        IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_pairwise_distances(input, output);
}

static int
gfp_pairwise_distances_random (IWString_and_File_Descriptor & output)
{
  int i1 = intbtwij (0, pool_size - 1);
  int i2;
  do 
  {
    i2 = intbtwij (0, pool_size - 1);
  }
  while (i2 == i1);

  similarity_type_t d = compute_the_distance (i1, i2, tversky);

  if (verbose)
    distance_acc.extra(d);

  output << pool[i1].id() << ' ' << pool[i2].id() << ' ' << d << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
gfp_pairwise_distances (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vp:s:F:V:P:gK:T:o:aqr:b");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! initialise_fingerprints (cl, verbose))
  {
    cerr << "Cannot initialise fingerprint specifications\n";
    usage (23);
  }

  if (cl.option_present('V'))
  {
    if (! tversky.parse_command_line(cl, 'V', verbose))
    {
      cerr << "Cannot initialise Tversky specifications (-V)\n";
      usage(5);
    }
  }

  if (! cl.option_present('p'))
  {
    cerr << "Must specify a file of fingerprints via the -p option\n";
    usage(3);
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      cerr << "The pool size (-s) must be a whole number > 1\n";
      usage(4);
    }

    pool = new IW_General_Fingerprint[pool_size];

    if (verbose)
      cerr << "Pool sized for " << pool_size << " fingerprints\n";
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', upper_distance_threshold) || upper_distance_threshold <= 0.0)
    {
      cerr << "Invalid upper distance threshold (-T)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will ignore all distances above " << upper_distance_threshold << endl;
  }

  if (cl.option_present('g'))
  {
    ignore_duplicate_identifiers_in_pool = 1;
    if (verbose)
      cerr << "Will ignore duplicate identifiers in the -p file\n";
  }

  if (cl.option_present('K'))
    ;
  else if (cl.option_present('r'))
    ;
  else if (cl.number_elements() > 0)
    ;
  else
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('p'))
  {
    const char * p = cl.option_value('p');

    if (! read_pool(p))
    {
      cerr << "Cannot read fingerpritns from '" << p << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Read " << pool_size << " fingerprints from '" << p << "'\n";
  }

  if (cl.option_present('o'))
  {
    int o;
    if (! cl.value('o', o) || o < 1)
    {
      cerr << "The output precision value must be a valid precision\n";
      usage(3);
    }

    if (verbose)
      cerr << "Default output precision " << o << endl;

    set_default_iwstring_float_concatenation_precision(o);
  }

  if (cl.option_present('a'))
  {
    write_entire_input_records = 1;

    if (verbose)
      cerr << "Will write all tokens in identifier file\n";
  }

  if (cl.option_present('q'))
  {
    equal_weight_tanimoto = 1;

    if (verbose)
      cerr << "Will use equal weight Tanimoto measure\n";
  }

  if (cl.option_present('b'))
  {
    ignore_missing_identifiers = 1;

    if (verbose)
      cerr << "Will ignore requests for identifiers not in the fingerprint file\n";
  }

  IWString_and_File_Descriptor output(1);
  int rc = 0;

  if (cl.option_present('r'))
  {
    int randomly_choose;
    if (! cl.value('r', randomly_choose) || randomly_choose < 1)
    {
      cerr << "The number of random samples to process (-r) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will choose " << randomly_choose << " pairs at random\n";

    iw_set_rnum_seed (random_seed_from_dev_random());

    for (int i = 0; i < randomly_choose; i++)
    {
      gfp_pairwise_distances_random(output);
    }
  }
  else if (cl.option_present('K'))
  {
    int i = 0;
    IWString k;
    while (cl.value('K', k, i++))
    {
      if (! k.contains(','))
      {
        cerr << "Sorry, command line pairs must be comma separated '" << k << "' is invalid\n";
        rc = i + 1;
        break;
      }

      k.gsub(',', ' ');

      if (! gfp_pairwise_distance(k, output))
      {
        cerr << "Cannot process '" << k << "'\n";
        rc = i + 1;
        break;
      }
    }
  }
  
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! gfp_pairwise_distances(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    if (duplicate_identifiers_encountered)
      cerr << "Skipped " << duplicate_identifiers_encountered << " duplicate identifiers\n";

    if (distance_acc.n() > 1)
      cerr << distance_acc.n() << " distances between " << distance_acc.minval() << " and " << distance_acc.maxval() << " ave " << static_cast<float>(distance_acc.average()) << endl;

    if (cl.option_present('T'))
      cerr << discarded_by_upper_distance_threshold << " comparisons discarded by upper distance threshold " << upper_distance_threshold << endl;

    if (hist.active())
      hist.write_terse(cerr);
  }

  if (missing_identifiers_ignored)
    cerr << "Ignored " << missing_identifiers_ignored << " missing identifiers\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_pairwise_distances(argc, argv);

  return rc;
}
