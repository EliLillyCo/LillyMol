/*
  Spread implementation
*/

#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#include <algorithm>
#include <iostream>
#include <memory>

#include "tbb/scalable_allocator.h"

#include "re2/re2.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"

#include "Utilities/GFP_Tools/gfp_standard.h"

using std::cerr;
using std::endl;
using std::unique_ptr;

static int verbose = 0;

static int brief_output = 0;

static int output_at_end = 0;

/*
  Each worker will have a subset of the global set of fingerprints on which
  to work.

  Never implemented any of this C++ stuff, pthreads and objects don't much go together
*/

#ifdef WITH_OBJECTS

class Spread_Fingerprints
{
  private:
    GFP_Standard * _pool;       // we do not own this, it gets passed to us

    int _istart;
    int _istop;

    float * _distance;
    int * _nsn;

    pthread_t _wait1;
    pthread_t _wait2;

    IWString * _smiles;
    IWString * _pcn;

  public:
    Spread_Fingerprints();
    ~Spread_Fingerprints();

    int build (GFP_Standard * f, int, int);

    void object_has_been_selected (const GFP_Standard &);

    int spread (ostream &);
};

Spread_Fingerprints::Spread_Fingerprints ()
{
  _pool = nullptr;
  _istart = 0;
  _istop  = 0;
  _distance = nullptr;
  _nsn = nullptr;

  _smiles = nullptr;
  _pcn = nullptr;
}

Spread_Fingerprints::~Spread_Fingerprints()
{
  if (NULL != _distance)
    delete [] _distance;

  if (NULL != _nsn)
    delete _nsn;

  if (NULL != _smiles)
    delete _smiles;

  if (NULL != _pcn)
    delete _pcn;

  return;
}

int
Spread_Fingerprints::build (GFP_Standard * f, int a, int o)
{
  _pool = f;
  _istart = a;
  _istop = o;

  int n = o - a;

  _distance = new_float(n, 2.0f);
  _nsn = new_int(n, -1);
  _smiles = new IWString[n];
  _pcn = new IWString[n];

  return 1;
}
#endif

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");

static int nthreads = 0;

static pthread_t * child_thread = nullptr;

static sem_t * time_for_child_to_wake = nullptr;

static sem_t * child_done = nullptr;

static GFP_Standard * pool = nullptr;
static int pool_size = 0;
static IWString * smiles = nullptr;
static IWString * pcn = nullptr;
static int * nsn = nullptr;

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
    std::unique_ptr<re2::RE2> pcn_rx = std::make_unique<re2::RE2>("^PCN<");

    pool_size = input.grep (*pcn_rx);

    if (0 == pool_size)
    {
      cerr << "Zero occurrences of '" << pcn_rx->pattern() << "' in '" << fname << "'\n";
      return 0;
    }
  }

  pool = new GFP_Standard[pool_size];
  smiles = new IWString[pool_size];
  pcn = new IWString[pool_size];

  return build_pool (input, pool, pool_size);
}

static int
post_to_array_of_semaphores (sem_t * sem,
                             int nsem)
{
  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_post(&(sem[i])))
    {
      cerr << "Cannot do semaphore post, semphore " << i << endl;
      perror("sem_post");
      return 0;
    }
  }

  return 1;
}

static void
wait_for_children ()
{
  for (int i = 0; i < nthreads; i++)
  {
    sem_wait (child_done + i);
  }

  return;
}

static int
initialise_array_of_semaphores (sem_t * sem,
                                int nsem)
{
  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_init (&(sem[i]), 0, 0))
    {
      perror ("sem_init error");
      return 0;
    }
  }

  return nsem;
}

/*static void
do_object_has_been_selected (int isel)
{
  GFP_Standard & fpsel = pool[isel];

  for (int i = 0; i < pool_size; i++)
  {
    if (sdistance[i] < 0.0f || i == isel)
      continue;

    float d = 1.0f - fpsel.tanimoto(pool[i]);

    if (d < sdistance[i])
    {
      sdistance[i] = d;
      nsn[i] = isel;
//    cerr << "Distance for " << i << " updated to " << d << endl;
    }
  }

  return;
}*/

/*
  We need a structure via which we communicate between threads
*/

struct Thread_Communication
{
  int thread_number;
  int istart;
  int istop;

// In the first phase, we will fill these two values

  int most_distant;
  similarity_type_t furthest_distance;

// In the second phase, we will receive this value from the controlling
// process

  int sel;
};

/*static void
thread_identify_furthest_away (struct Thread_Communication * tc)
{
  if (tc->thread_number < 0)    // special value for finished
  {
    pthread_exit(NULL);
    return;
  }

  similarity_type_t maxdist = static_cast<similarity_type_t>(-1.0);
  int ichoose = -1;

  for (int i = tc->istart; i < tc->istop; i++)
  {
    if (sdistance[i] < 0.0f)
      continue;

//    cerr << "Distance to " << i << " is " << pool[i]->distance() << endl;
    if (sdistance[i] > maxdist)
    {
      maxdist = sdistance[i];
      ichoose = i;
    }
  }

  tc->most_distant = ichoose;
  tc->furthest_distance = maxdist;

  return;
}*/

//#define DEBUG_SPREAD_THREADS

static void
thread_object_has_been_selected (Thread_Communication * tc,
                                 float * sdistance)
{
#ifdef DEBUG_SPREAD_THREADS
  IWString msg;
  msg << "Thread " << tc->thread_number << " processing, sel " << tc->sel << "\n";
  cerr << msg;
#endif

  if (tc->thread_number < 0)    // special value for finished
  {
    pthread_exit(NULL);
    return;
  }

  int ichoose = -1;
  float maxdist = -1.0f;

  int sel = tc->sel;

  int istart = tc->istart;

  if (sel >= istart && sel < tc->istop)
    sdistance[sel - istart] = -1.0f;

//GFP_Standard & fpsel = pool[tc->sel];    // very important to test, should this be a reference or a copy?
  GFP_Standard & fpsel = pool[sel];   // need to test this...

  int istop = tc->istop - tc->istart;

  for (int i = 0; i < istop; i++)
  {
    float sd = sdistance[i];

    if (sd < 0.0f)
      continue;

    float d = 1.0f - fpsel.tanimoto(pool[istart + i]);
    if (d < sd)
    {
      sdistance[i] = d;
      sd = d;
      if (! brief_output)
        nsn[istart + i] = sel;
    }

    if (sd > maxdist)
    {
      maxdist = sd;
      ichoose = i;
    }
  }

  tc->furthest_distance = maxdist;
  tc->most_distant = istart + ichoose;

#ifdef DEBUG_SPREAD_THREADS
  msg.resize_keep_storage(0);
  msg << "Thread " << tc->thread_number << " done, furthest_distance " << tc->furthest_distance << " for iem " << tc->most_distant << "\n";
  cerr << msg;
#endif
  return;
}

/*
  The threads that look for the furthest away are very simple.
  Just wait until it is our turn, then scan and post
*/

static void *
thread_process (void * v)
{
  struct Thread_Communication * tc = reinterpret_cast<Thread_Communication *>(v);

  int mythread = tc->thread_number;

#ifdef DEBUG_SPREAD_THREADS
  IWString msg;
  msg << "Thread " << mythread << " launched, btw " << tc->istart << " and " << tc->istop << '\n';
  cerr << msg;
#endif

  float * d = new_float(tc->istop - tc->istart, 2.0f); unique_ptr<float[]> free_d(d);

  while (1)
  {
    sem_wait (time_for_child_to_wake + mythread);
    thread_object_has_been_selected(tc, d);
    sem_post(child_done + mythread);
    sched_yield();
  }
}

static struct Thread_Communication ** tc = nullptr;

static int
initialise_threads (int nthreads,
                    int first_selected)
{
  assert (nthreads > 0);
  assert (NULL != tc);

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1)
  {
    cerr << "initialise_threads:too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init (&attr);
  pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++)
  {
    Thread_Communication * tci = tc[i];

    tci->thread_number = i;
    tci->istart = i * items_per_chunk;
    if (i == nthreads - 1)
      tci->istop = pool_size;
    else
      tci->istop = tci->istart + items_per_chunk;

    if (verbose > 1)
      cerr << "Thread processes from " << tci->istart << " to " << tci->istop << endl;

    tci->sel = first_selected;

    if (0 != pthread_create (&(child_thread[i]), &attr, thread_process, tci))
    {
      perror ("cannot create child process");
      return 0;
    }

#ifdef DEBUG_SPREAD_THREADS
    cerr << "Launched thread " << i << endl;
#endif
  }

  return 1;
}

static void
do_output (int isel,
           float d,
           IWString_and_File_Descriptor & output)
{
  if (brief_output)
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
fpobj_spread (GFP_Standard * pool,
              int pool_size,
              int number_to_select,
              IWString_and_File_Descriptor & output)
{
  int first_selected = 0;

  int number_selected = 0;

  if (! initialise_threads(nthreads, first_selected))
  {
    cerr << "Cannot initialise threads\n";
    return 0;
  }

  int * sel = nullptr;
  float * sdistance = nullptr;
  if (output_at_end)    // free these sometime
  {
    sdistance = new float[pool_size];
    sel = new int[pool_size];
  }

  int ichoose = first_selected;
  float d = 2.0f;

  while (number_selected < number_to_select)
  {
    if (! post_to_array_of_semaphores(time_for_child_to_wake, nthreads))
      return 0;

//  Do output while children computing

    if (output_at_end)
    {
      sel[number_selected] = ichoose;
      sdistance[number_selected] = d;
    }
    else
      do_output (ichoose, d, output);

    number_selected++;

//  wait_for_children();

    if (number_to_select == number_selected)
      break;

    d = -1.0f;
    ichoose = -1;

    for (int i = 0; i < nthreads; i++)
    {
#ifdef DEBUG_SPREAD_THREADS
      IWString msg;
      msg << "Thread " << i << " found " << tc[i].furthest_distance << "\n";
      cerr << msg;
#endif
      sem_wait(child_done + i);

      if (tc[i]->furthest_distance > d)
      {
        d = tc[i]->furthest_distance;
        ichoose = tc[i]->most_distant;
      }
    }

#ifdef DEBUG_SPREAD_THREADS
    cerr << "Most distant is " << ichoose << " at distance " << d << endl;
#endif

    assert (ichoose >= 0);

    for (int i = 0; i < nthreads; i++)  // tell each of the threads the identity of the selected fp
    {
      tc[i]->sel = ichoose;
    }

    if (verbose > 1)
      cerr << "Selected " << number_selected << " '" << pcn[ichoose] << "' distance " << d << " NSN '" << pcn[nsn[ichoose]] << "'\n";
  }

  if (output_at_end)
  {
    for (auto i = 0; i < number_to_select; ++i)
    {
      do_output(sel[i], sdistance[i], output);
    }
  }

  output.flush();

// tell all the threads to destroy themselves

  for (int i = 0; i < nthreads; i++)
  {
    tc[i]->thread_number = -1;
  }

  post_to_array_of_semaphores(time_for_child_to_wake, nthreads);

  if (verbose)
    cerr << "Returning with " << number_selected << " items selected\n";

  return number_selected;
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
  Command_Line cl (argc, argv, "vs:n:h:be");

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

  if (! cl.option_present('h'))
  {
    cerr << "Must specify number of threads via the -h option\n";
    usage(3);
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', nthreads) || nthreads < 1)
    {
      cerr << "The number of threads to use (-h) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will use " << nthreads << " threads\n";

    time_for_child_to_wake = new sem_t[nthreads];
    child_done = new sem_t[nthreads];

    if (! initialise_array_of_semaphores(time_for_child_to_wake, nthreads))
      return 3;

    if (! initialise_array_of_semaphores(child_done, nthreads))
      return 3;

    tc = new struct Thread_Communication *[nthreads];
    for (auto i = 0; i < nthreads; ++i)
    {
      tc[i] = new struct Thread_Communication;
    }
  }

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

  nsn = new_int(pool_size, -1);

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

  output.flush();

  if (nthreads)
  {
    destroy_array_of_semaphores(time_for_child_to_wake, nthreads);
    destroy_array_of_semaphores(child_done, nthreads);
  }

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
