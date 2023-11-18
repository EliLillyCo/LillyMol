/*
  Spread implementation
*/

#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

#include <cilk/cilk.h>
#include "tbb/scalable_allocator.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/cmdline/cmdline.h"
#include "iwrandom.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"

#include "gfp_standard.h"

static int verbose = 0;

static int brief_output = 0;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");

static int output_at_end = 0;

static GFP_Standard * pool = nullptr;
static int pool_size = 0;
static float * sdistance = nullptr;
static IWString * smiles = nullptr;
static IWString * pcn = nullptr;
static int * sel = nullptr;
static int * nsn = nullptr;
static int * fll = nullptr;

static int
build_pool (iwstring_data_source & input,
            GFP_Standard * & pool,
            int pool_size)
{
  assert (pool_size > 0);
//cerr << "Pool ptr " << poolptr << ", pool size " << pool_size << endl;

  int ndx = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    IW_General_Fingerprint gfp;

    int fatal;
    if (! gfp.construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }
    }
    else
    {
      tdt.dataitem_value(smiles_tag, smiles[ndx]);
      tdt.dataitem_value(identifier_tag, pcn[ndx]);

      pool[ndx].build_molecular_properties(gfp.molecular_properties_integer());
      pool[ndx].build_iw(gfp[0]);
      pool[ndx].build_mk(gfp[1]);
      pool[ndx].build_mk2(gfp[2]);
      ndx++;
      if (ndx >= pool_size)
        break;
    }
  }

  if (0 == ndx)
  {
    cerr << "No fingerprints read\n";
    return 0;
  }

  return 1;
}

static int
build_pool (const char * fname,
            GFP_Standard * & pool,
            int & pool_size)
{
  iwstring_data_source input (fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (pool_size <= 0)
  {
    std::unique_ptr<re2::RE2> pcn_rx ("^PCN<");

    pool_size = input.grep (pcn_rx);

    if (0 == pool_size)
    {
      cerr << "Zero occurrences of '" << pcn_rx.source() << "' in '" << fname << "'\n";
      return 0;
    }
  }

  pool = new GFP_Standard[pool_size];
  smiles = new IWString[pool_size];
  pcn = new IWString[pool_size];
  fll = new int[pool_size];
  for (auto i = 0; i < pool_size; ++i)
  {
    fll[i] = i;
  }

  return build_pool (input, pool, pool_size);
}


static void
do_output (int number_selected,
           int isel,
           float d,
           IWString_and_File_Descriptor & output)
{
  if (output_at_end)
    sel[number_selected] = isel;
  else if (brief_output)
    output << smiles[isel] << ' ' << pcn[isel] << ' ' << d << "\n";
  else
  {
    output << smiles_tag     << smiles[isel] << ">\n";
    output << identifier_tag << pcn[isel] << ">\n";

    if (nsn[isel] >= 0)
    {
      output << smiles_tag     << smiles[nsn[isel]] << ">\n";
      output << identifier_tag << pcn[nsn[isel]] << ">\n";
      output << distance_tag   << d << ">\n";
    }
    else
    {
      output << smiles_tag << "*>\n";
      output << identifier_tag << "*>\n";
      output << distance_tag << "1>\n";
    }
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(4096);

  return;
}

static int
process_selected_item (int isel)
{
  sdistance[isel] = -1.0f;

  const GFP_Standard & fpsel = pool[isel];

  _Cilk_for (auto i = 0; i < pool_size; i++)
//for (auto i = 0; i < pool_size; i++)
  {
    if (sdistance[i] < 0.0f)
      continue;

    const float d = 1.0f - fpsel.tanimoto(pool[i]);
//  cerr << i << " distance currently " << sdistance[i] << " to new fp " << d << ", furthest is " << _furthest_distance << endl;
    if (d < sdistance[i])
    {
      sdistance[i] = d;
      nsn[i] = isel;
//    cerr << i << " distance to nearest selected updated to " << sdistance[i] << endl;
    }
  }

  return 1;
}

static int
fpobj_spread (GFP_Standard * pool,
              int pool_size,
              int number_to_select,
              IWString_and_File_Descriptor & output)
{
  int first_selected = 0;

  sdistance[first_selected] = -1.0f;

  int ichoose = first_selected;

  int number_selected = 0;

  while (1)
  {
    do_output(number_selected, ichoose, sdistance[ichoose], output);

    if (verbose > 1)
      cerr << "Selected " << number_selected << " '" << smiles[ichoose] << "' distance " << sdistance[ichoose] << " NSN '" << pcn[nsn[ichoose]] << "'\n";

    process_selected_item (ichoose);

    number_selected ++;

    if (number_selected >= number_to_select)
      break;

    ichoose = __sec_reduce_max_ind(sdistance[0:pool_size]);
  }

  return 1;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many items to select\n";
  cerr << " -h <number>      number of threads to use\n";
  cerr << " -b               brief output only, 'smiles id dist', no need for nplotnn\n";
  cerr << " -e               postpone output until job is done\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}


static int
destroy_array_of_semaphores (sem_t * sem,
                             int nsem)
{
  int rc = 1;

  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_destroy(&(sem[i])))
    {
      perror ("sem_destroy error");
      rc = 0;
    }
  }

  return rc;
}

static int
fpobj_spread (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:h:bg:e");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (1);
  }

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise GFP options\n";
      usage (23);
    }
  }
  else if (! initialise_fingerprints (cl[0], verbose))
  {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  if (cl.option_present('b'))
  {
    brief_output = 1;

    if (verbose)
      cerr << "Brief output only, 'smiles id dist'\n";
  }

#ifdef CANNOT_DO_THIS_WITH_CILK
  if (cl.option_present('g'))
  {
    if (! cl.value('g', grainsize) || grainsize < 1)
    {
      cerr << "The grain size must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Grain size set to " << grainsize << endl;
  }
#endif

  if (cl.option_present('e'))
  {
    output_at_end = 1;

    if (verbose)
      cerr << "Output postponed until end of job\n";
  }

#ifdef SQUEEZING_IMPLEMENTED
  if (cl.option_present('q'))
  {
    if (! cl.value('q', squeeze_selected) || squeeze_selected < 1)
    {
      cerr << "The squeeze selected items option (-q) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will squeeze out selected items every " << squeeze_selected << " items selected\n";
  }
#endif

#ifdef CAN_SET_CILK_THREADS
  if (cl.option_present('h'))
  {
    int h;
    if (! cl.value('h', h) || h < 1)
    {
      cerr << "The number of threads must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will use " << h << " threads\n";

    init = new tbb::task_scheduler_init(h);
  }
  else
    init = new tbb::task_scheduler_init();
#endif

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, only processes a single input file\n";
    usage(2);
  }

  if (! build_pool(cl[0], pool, pool_size))
  {
    cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
    return 1;
  }

  if (verbose)
    cerr << "Read " << pool_size << " fingerprints from '" << cl[0] << "'\n";

  sdistance = new_float(pool_size, 2.0f);
  nsn = new_int(pool_size, -1);
  sel = new int[pool_size];

  int number_to_select = pool_size;

  if (cl.option_present ('n'))
  {
    int n;
    if (! cl.value ('n', n) || n < 1)
    {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage (13);
    }
    
    if (n > pool_size)
    {
      cerr << "You asked for " << n << " molecules, but pool only contains " << pool_size << ". Shortened\n";
      n = pool_size;
    }

    number_to_select = n;
    if (verbose)
      cerr << number_to_select << " molecules will be selected\n";
  }

  IWString_and_File_Descriptor output(1);

  (void) fpobj_spread (pool, pool_size, number_to_select, output);

  if (output_at_end)
  {
    output_at_end = 0;
    for (auto i = 0; i < number_to_select; ++i)
    {
      do_output(i, sel[i], sdistance[i], output);
    }
  }

  output.flush();

//delete [] pool;     leave this out for efficiency

  cerr << "Output can be processed with nplotnn\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = fpobj_spread (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
