#include <stdlib.h>

#include <iostream>

/*
  Sometimes we want all those molecules from one file which are within
  a given distance of those in another file
*/

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#define IWMINMAX_IMPLEMENTATION
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp_standard.h"

#ifdef __GNUG__
template iwminid<float, int>::iwminid(float, int);
template int iwminid<float, int>::try_this(float, int);
#endif

using std::cerr;
using std::endl;

/*
  Our pool is an array of FP objects
*/

static GFP_Standard * pool = nullptr;

static int pool_size = 0;

static int verbose = 0;

static const char * prog_name;

static Window_Specification<GFP_Standard> wspec;

static unsigned int tdts_read = 0;

static unsigned int tdts_written = 0;

static Report_Progress report_progress;

static unsigned long long computations_done = 0;

static IWString_and_File_Descriptor stream_for_non_matches;

static int ignore_zero_distances = 0;

/*
  By default, we just search the pool until we find something that is
  within the window. If the user wants the closest distance, we need
  to search the whole pool
*/

static int find_shortest_distance = 0;

static IWString shortest_distance_tag;

static int write_smiles = 0;

static IWString smiles_tag("$SMI<");

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag ("PCN<");

/*
  We can optionally write the closest distance (within the threshold ranges)
  and write that to the output.
*/

static const_IWSubstring distance_tag ("DIST<");

static similarity_type_t upper_distance_threshold = -1.0;

static similarity_type_t lower_distance_threshold = -1.0;

/*
  The lower and upper threshold are different.

  The lower threshold is designed for keeping molecules away from
  the pool.

  A molecule will fail if it violates the lower threshold a specified
  number of times. That is, we are trying to make sure the molecule
  does not get too close to anything in the pool.

  The upper threshold is designed for keeping molecules close to
  the pool. For success, there must be some number of molecules that
  come within the upper threshold.
*/

static int lower_threshold_violation_threshold = 1;

static int upper_threshold_success_requirement = 1;

#ifdef __GNUG__
template int IW_TDT::_add_dataitem<float>(char const*, int, float const&, int);
#endif

static void
add_shortest_distance (IW_TDT & tdt,
                       similarity_type_t mindist)
{
  tdt.add_dataitem (shortest_distance_tag, mindist);

  return;
}

static int
build_pool (iwstring_data_source & input)
{
  int items_in_pool = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    IW_General_Fingerprint gfp;

    int fatal;
    if (! gfp.construct_from_tdt(tdt, fatal))
    {
      cerr << "Cannot read fingerprint\n";
      return 0;
    }

    pool[items_in_pool].build_molecular_properties(gfp.molecular_properties_integer());
    pool[items_in_pool].build_iw(gfp[0]);
    pool[items_in_pool].build_mk(gfp[1]);
    pool[items_in_pool].build_mk2(gfp[2]);

    items_in_pool++;

    if (items_in_pool >= pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      return 1;
    }
  }

  pool_size = items_in_pool;

  return 1;
}

static int
build_pool (const const_IWSubstring & fname)
{
  iwstring_data_source input (fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    IWString tmp;
    tmp << '^' << identifier_tag;

    std::unique_ptr<re2::RE2> pcn = std::make_unique<re2::RE2>(tmp);
    pool_size = input.grep(*pcn);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    pool = new GFP_Standard[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool (input);
}

/*
  Return 1 if this fingerprint passes the filter
*/

static int
distance_filter (GFP_Standard & fp,
                 iwminid<similarity_type_t, int> & mindist)
{
// We keep track of the number of times this molecule violates thresholds

  int lower_threshold_violation_count = 0;
  int upper_threshold_success_count = 0;

  int keep_going = 1;

  for (int i = 0; i < pool_size && keep_going; i++)
  {
    if (! wspec.can_be_compared (fp, pool[i]))
      continue;

    computations_done++;

    auto t = fp.tanimoto_distance (pool[i]);

//#define DEBUG_NN
#ifdef DEBUG_NN
    cerr << "Distance between '" << fp.id() << " and pool " << i << " '" << pool[i].id() << "' is " << t << endl;
#endif

    if (t > static_cast<similarity_type_t> (0.0))
      ;
    else if (ignore_zero_distances)
      continue;

    if (t < lower_distance_threshold)
    {
      lower_threshold_violation_count++;
      if (lower_threshold_violation_count >= lower_threshold_violation_threshold)
        keep_going = 0;
    }

    if (t < upper_distance_threshold)
    {
      upper_threshold_success_count++;
      if (upper_threshold_success_count >= upper_threshold_success_requirement)
      {
        mindist.try_this(t, i);
        keep_going = 0;
      }
    }

    if (find_shortest_distance)
    {
      mindist.try_this(t, i);
      keep_going = 1;    // finding shortest distance, must check everything
    }
  }

  if (lower_threshold_violation_count >= lower_threshold_violation_threshold)
    return 0;

  if (upper_distance_threshold > 0.0 && upper_threshold_success_count < upper_threshold_success_requirement)
    return 0;

  return 1;
}

static int
do_write_smiles (const IW_TDT & tdt,
                 const IW_General_Fingerprint & fp,
                 const similarity_type_t mindist,
                 IWString_and_File_Descriptor & output)
{
  const_IWSubstring smi;
  if (! tdt.dataitem_value(smiles_tag, smi))
  {
    cerr << "Yipes, no smiles for '" << fp.id () << "'\n";
    return 0;
  }

  output << smi << ' ' << fp.id();
  
  if (shortest_distance_tag.length())
    output << ' ' << mindist;

  output << '\n';

  return 1;
}

static int
do_output (IW_TDT & tdt,
           const IW_General_Fingerprint & fp,
           const iwminid<similarity_type_t, int> & mindist,
           IWString_and_File_Descriptor & output)
{
  if (write_smiles)
    do_write_smiles(tdt, fp, mindist.minval(), output);
  else
  {
    if (shortest_distance_tag.length())
      add_shortest_distance (tdt, mindist.minval());

    output << tdt;
  }

  tdts_written++;
  if (verbose > 1)
    cerr << fp.id() << " written, total written " <<  tdts_written << endl;

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
distance_filter (iwstring_data_source & input, 
                 IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;
  while (tdt.next (input) && output.good())
  {
    IW_General_Fingerprint gfp;
    int fatal;

    if (! gfp.construct_from_tdt (tdt, fatal))
    {
      if (fatal)
        return 0;

      continue;
    }
    tdts_read++;

    GFP_Standard stdfp;
    stdfp.build_molecular_properties(gfp.molecular_properties_integer());
    stdfp.build_iw(gfp[0]);
    stdfp.build_mk(gfp[1]);
    stdfp.build_mk2(gfp[2]);

    if (report_progress())
      cerr << tdts_read << " fingerprints processed, " << tdts_written << " written\n";

    iwminid<similarity_type_t, int> mindist(1.01, -1);

    if (distance_filter (stdfp, mindist))
      do_output(tdt, gfp, mindist, output);
    else if (stream_for_non_matches.is_open())
    {
      if (write_smiles)
        do_write_smiles(tdt, gfp, mindist.minval(), stream_for_non_matches);
      else
      {
        if (shortest_distance_tag.length())
          add_shortest_distance(tdt, mindist.minval());

        stream_for_non_matches << tdt;
      }
      stream_for_non_matches.write_if_buffer_holds_more_than(8192);
    }
  }

  return output.good();
}

static int
distance_filter (const char * fname, 
                 IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);
  if (! input.ok())
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 0;
  }

  return distance_filter (input, output);
}


static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Filter molecules according to how close they come to members of a pool\n";
  cerr << prog_name << ": usage <options> <input_file>\n";
  cerr << " -p <file>        specify file against which input is to be compared\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -t <dis>         lower distance threshold\n";
  cerr << " -n <n>           reject molecules that violate lower threshold <n> times or more\n";
  cerr << " -T <dis>         upper distance threshold\n";
  cerr << " -N <n>           reject molecules that violate upper threshold <n> times or more\n";
  cerr << "To pass the thresholds, a distance must be >= lower && <= upper\n";
  cerr << "                  must be less <= pool size\n";
  cerr << " -f               write smiles as output\n";
  cerr << " -U <file>        write molecules that fail the filter to <file>\n";
  cerr << " -D <tag>         add distance to output (forces closest distance determination)\n";
  cerr << " -d <tag>         add distance to output (no shortest distance determination)\n";
  cerr << " -z               ignore zero distances\n";
  cerr << " -h               ignore cases with zero distance and the same ID\n";
  cerr << " -r <n>           report progress every <n> items processed\n";
  cerr << " -F ...           standard fingerprint options, enter '-F help'\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
distance_filter (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:p:t:T:F:P:r:D:d:n:N:U:zfO:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present('t') && cl.option_present('T'))
  {
    cerr << "Sorry, the -t and -T options cannot be used together\n";
    usage(3);
  }

  if (cl.option_present('n') && ! cl.option_present('t'))
  {
    cerr << "The lower threshold violation count option (-n) only makes sense with the -t option\n";
    usage(2);
  }

  if (cl.option_present('N') && ! cl.option_present('T'))
  {
    cerr << "The upper threshold success requirement option (-N) only makes sense with the -T option\n";
    usage(2);
  }

  if (cl.option_present ('t'))
  {
    if (! cl.value ('t', lower_distance_threshold) || lower_distance_threshold < 0.0 || lower_distance_threshold > 1.0)
    {
      cerr << "The -t option must be followed by a valid similarity value\n";
      usage (12);
    }

    if (verbose)
      cerr << "Lower distance threshold set to " << lower_distance_threshold << endl;
  }

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', upper_distance_threshold) || upper_distance_threshold < 0.0 || upper_distance_threshold > 1.0)
    {
      cerr << "The -T option must be followed by a valid similarity value\n";
      usage (12);
    }

    if (verbose)
      cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
  }

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', lower_threshold_violation_threshold) || lower_threshold_violation_threshold < 0)
    {
      cerr << "Invalid value for the lower threshold violation count (-n)\n";
      usage (19);
    }

    if (verbose)
      cerr << "Will fail molecules that violate lower threshold " << lower_threshold_violation_threshold << " times or more\n";
  }

  if (cl.option_present ('N'))
  {
    if (! cl.value ('N', upper_threshold_success_requirement) || upper_threshold_success_requirement < 0)
    {
      cerr << "Invalid value for the upper threshold success requirement (-N)\n";
      usage (19);
    }

    if (verbose)
      cerr << "Will fail molecules unless " << lower_threshold_violation_threshold << " times within upper threshold\n";
  }

  if (cl.option_present('f'))
  {
    write_smiles = cl.option_count('f');
    if (verbose)
      cerr << "Will write smiles\n";
  }

  if (cl.option_present ('z'))
  {
    ignore_zero_distances = 1;
    if (verbose)
      cerr << "Will ignore zero distances\n";
  }

  if (cl.option_present('O'))
  {
    if (! wspec.build(cl, 'O', verbose))
    {
      cerr << "Cannot initialise window comparison specification (-O)\n";
      return 2;
    }
  }

// Kludge for Knime so that '-F help' works in all cases

  if (cl.option_present('F') && "help" == cl.string_value('F'))
  {
    display_standard_gfp_options(cerr);
    return 3;
  }

  if (! cl.option_present ('p'))
  {
    cerr << "Must specify a pool file via the -p option\n";
    usage (5);
  }

  if (1 != cl.option_count ('p'))
  {
    cerr << "Only one pool file may be specified\n";
    usage (6);
  }

// We need to be very careful if we are reading from stdin

  const char * p = cl.option_value('p');

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))
  {
    if (need_to_call_initialise_fingerprints (cl))
    {
      if (! initialise_fingerprints (p, verbose))
      {
        cerr << "Cannot initialise GFP options\n";
        usage (23);
      }
    }
    else if (! initialise_fingerprints(p, verbose))
    {
      cerr << "Cannot initialise fingerprints from '" << p << "'\n";
      return 4;
    }
  }
  else
  {
    if (need_to_call_initialise_fingerprints (cl))
    {
      if (! initialise_fingerprints (cl, verbose))
      {
        cerr << "Cannot initialise GFP options\n";
        usage (23);
      }
    }
    else if (! initialise_fingerprints(cl[0], verbose))
    {
      cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
      return 4;
    }
  }

  if (cl.option_present ('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "The -r option must be followed by a whole positive number\n";
      usage (18);
    }
  }

  if (cl.option_present ('D'))
  {
    find_shortest_distance = 1;

    shortest_distance_tag = cl.string_value ('D');

    if (verbose)
      cerr << "Will write shortest distance as '" << shortest_distance_tag << "' dataitem\n";

    if (! shortest_distance_tag.ends_with('<'))
      shortest_distance_tag << '<';
  }

  if (cl.option_present ('d'))
  {
    shortest_distance_tag = cl.string_value ('d');

    if (verbose)
      cerr << "Will write shortest distance as '" << shortest_distance_tag << "' dataitem\n";

    if (! shortest_distance_tag.ends_with('<'))
      shortest_distance_tag << '<';
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (1);
  }

  if (cl.option_present ('U'))
  {
    IWString n = cl.string_value ('U');

    if (! stream_for_non_matches.open (n.null_terminated_chars()))
    {
      cerr << "Cannot open stream for non-matches '" << n << "'\n";
      return 17;
    }

    if (verbose)
      cerr << "Will write non-matches to '" << n << "'\n";
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }

    pool = new GFP_Standard[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose)
      cerr << "system sized to " << pool_size << endl;

    if (0 == lower_threshold_violation_threshold)
      ;
    else if (lower_threshold_violation_threshold > pool_size)
    {
      cerr << "The lower threshold violation count " << lower_threshold_violation_threshold << " must be <= pool size " << pool_size << endl;
      usage (7);
    }

    if (0 == upper_threshold_success_requirement)
      ;
    else if (upper_threshold_success_requirement > pool_size)
    {
      cerr << "The upper threshold success requirement " << upper_threshold_success_requirement << " must be <= pool size " << pool_size << endl;
      usage (7);
    }
  }

  if (cl.option_present ('p'))
  {
    const_IWSubstring fname;
    cl.value ('p', fname);

    if (! build_pool (fname))
    {
      cerr << "Cannot build pool from '" << fname << "'\n";
      return 76;
    }
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < cl.number_elements(); i++)
  {
    (void) distance_filter (cl[i], output);
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << tdts_read << " wrote " << tdts_written << endl;
    unsigned long long possible_computations = static_cast<unsigned long long> (pool_size) * static_cast<unsigned long long> (tdts_read);
    cerr << "Performed " << computations_done << " of " << possible_computations << " possible distance computations " << static_cast<float>(computations_done) / static_cast<float>(possible_computations) << endl;
  }


  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_filter (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
