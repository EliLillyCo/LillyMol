/*
  Leader implementation with fixed fingerprints
*/

#include <stdlib.h>
#include <semaphore.h>
#include <pthread.h>
#include <iostream>
#include <memory>

#include "tbb/scalable_allocator.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iw_tdt/iw_tdt.h"


#include "Utilities/GFP_Tools/gfp_standard.h"
#include "Utilities/GFP_Tools/smiles_id_dist.h"

using std::cerr;
using std::cout;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static Report_Progress report_progress;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");
static IWString mk_tag ("FPMK<");
static IWString mk2_tag ("FPMK2<");
static IWString iw_tag ("FPIW<");

static float threshold = 0.0;

static int reallocate_unselected_items = 0;

static int
initialise_array_of_semaphores (sem_t * sem,
                                int nsem)
{
  assert(nsem > 0);

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

static int
destroy_array_of_semaphores (sem_t * sem,
                             int nsem)
{
  int rc = 1;

  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_destroy(&sem[i]))
    {
      perror("sem_destroy");
      rc = 0;
    }
  }

  return rc;
}

class Leader_Item : public Smiles_ID_Dist
{
  private:
    const GFP_Standard * _gfp;

  public:
    Leader_Item();

    void set_gfp (const GFP_Standard * s) { _gfp = s;}
};

Leader_Item::Leader_Item()
{
  _gfp = nullptr;

  return;
}

static int pool_size = 0;

static Leader_Item * leader_item = nullptr;
static GFP_Standard * fingerprints = nullptr;

struct Thread_Info
{
  int   thread_number;
  const GFP_Standard * leader;
  const GFP_Standard * fp;
  int   istart;
  int   istop;
  float threshold;
  int   cluster_number;
  int   next_leader;
  int   next_leader_index;
  int   * selected;
  float * distance;
};

static int nthreads = 0;
static Thread_Info * thread_info = nullptr;
static sem_t * leader_identified = nullptr;
static sem_t cluster_formed;
static pthread_t * child_thread = nullptr;

static void *
form_clusters_thread (void * v)
{
  Thread_Info * t = reinterpret_cast<Thread_Info *>(v);

  int mythread = t->thread_number;
  const GFP_Standard * fp = t->fp;
  float threshold =  t->threshold;
  int * selected = t->selected;
  float * distance = t->distance;

//cerr << "Spawn thread from " << t->istart << " to " << t->istop << ", thread " << mythread << endl;

  while (1)
  {
    sem_wait(leader_identified + mythread);

    int cluster_number = t->cluster_number;
    if (cluster_number < 0)
    {
      sem_post(&cluster_formed);
      pthread_exit(NULL);
      return NULL;
    }

    int istart = t->istart;
    int istop = t->istop;
    const GFP_Standard leader = *(t->leader);
    int ndx = 0;
    int next_leader = -1;
    int next_leader_index = -1;

    for (int i = istart; i < istop; i++, ndx++)
    {
//    cerr << " ndx " << ndx << " selected " << selected[ndx] << endl;

      if (selected[ndx])
        continue;

      float d = leader.tanimoto_distance(fp[i]);
      if (d <= threshold)
      {
        selected[ndx] = cluster_number;
        distance[ndx] = d;
      }
      else if (-1 == next_leader)
      {
        next_leader = i;
        next_leader_index = ndx;
      }
    }

    t->next_leader = next_leader;
    t->next_leader_index = next_leader_index;

    sem_post(&cluster_formed);
  }

  return NULL;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Leader implementation requiring MPR IW MK MK2\n";
  cerr << " -t <dist>      distance threshold\n";
  cerr << " -A <fname>     file of previously selected fingerprints\n";
  cerr << " -s <size>      size of fingerprint file\n";
  cerr << " -p <n>         process <n> leaders at once (max 2 for now)\n";
  cerr << " -h <n>         maximum number of OMP threads to use\n";
  cerr << " -r <n>         report progress every <n> clusters formed\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
choose_next_leader (const Thread_Info * thread_info,
                    int nthreads,
                    int cluster_id)
{
//cerr << "nk0 " << thread_info[0].next_leader << " ml1 " << thread_info[1].next_leader << endl;
  for (int i = 0; i < nthreads; i++)
  {
//  cerr << "Thread " << i << " next_leader " << thread_info[i].next_leader << " ndx " << thread_info[i].next_leader_index << endl;
    if (thread_info[i].next_leader >= 0)
    {
      int n = thread_info[i].next_leader;

      int ndx = thread_info[i].next_leader_index;

      thread_info[i].selected[ndx] = cluster_id;
      thread_info[i].distance[ndx] = 0.0f;

      return n;
    }
  }

  return -1;
}

#ifdef NOT_BEING_USED
static int
choose_next_leader (const int * selected,
                    int n)
{
  for (int i = 0; i < n; i++)
  {
//      cerr << " CNL i = " << i << " sel " << selected[i] << endl;
    if (! selected[i])
      return i;
  }

  return -1;
}

static int
choose_next_leaders (const int * selected,
                     int n,
                     int * leaders,
                     int leaders_requested)
{
  int leader_count = 0;

  for (int i = 0; i < n; i++)
  {
    if (selected[i])
      continue;

    leaders[leader_count] = i;
    leader_count++;

    if (leader_count == leaders_requested)
      return leaders_requested;
  }

  return leader_count;
}
#endif

static void
do_reallocate_unselected_items (Thread_Info * thread_info,
                                int nthreads,
                                int * global_selected_array,
                                float * global_distance_array)
{
  int unselected = 0;
  int offset = 0;

  for (int i = 0; i < nthreads; i++)
  {
    int n = thread_info[i].istop - thread_info[i].istart;
    const int * s = thread_info[i].selected;
    unselected = count_occurrences_of_item_in_array(0, n, s);
    copy_vector(global_selected_array + offset, s, n);
    copy_vector(global_distance_array + offset, thread_info[i].distance, n);
    offset += n;
  }

  cerr << "Reallocating " << unselected << " unselected items\n";

  if (unselected < 1000)
    return;

  int items_per_chunk = unselected / nthreads;

  if (items_per_chunk < 10)
    return;

  int ndx = 0;

  unselected = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (global_selected_array[i])
      continue;

    unselected++;
    if (1 == unselected)
    {
      thread_info[ndx].istart = i;
      continue;
    }

    if (unselected >= items_per_chunk)
    {
      thread_info[ndx].istop = i + 1;
      ndx++;
      if (ndx == nthreads)    // processing last thread
      {
        while (i < pool_size)
        {
          if (global_selected_array[i])
            thread_info[ndx].istop = i;
          i++;
        }

        return;
      }

      unselected = 0;
    }
  }

  return;
}


static int
leader (Thread_Info * thread_info,
        int nthreads,
        int * global_selected_array,
        float * global_distance_array)
{
  int cluster_id = 1;

  int icentre;

  int * selected = thread_info[0].selected;
  float * distances = thread_info[0].distance;

  thread_info[0].next_leader = 0;
  thread_info[0].next_leader_index = 0;

  while ((icentre = choose_next_leader (thread_info, nthreads, cluster_id)) >= 0)
  {
    if (report_progress())
      cerr << "Formed " << cluster_id << " clusters, " << count_non_zero_occurrences_in_array(selected, pool_size) << " items in clusters\n";

//#define DEBUG_LEADER1
#ifdef DEBUG_LEADER1
    cerr << "Selected " << icentre << " at centre of cluster " << cluster_id << endl;
#endif

    for (int i = 1; i < nthreads; i++)
    {
      thread_info[i].leader = fingerprints + icentre;
      thread_info[i].threshold = threshold;
      thread_info[i].cluster_number = cluster_id;
    }
    for (int i = 1; i < nthreads; i++)
    {
      sem_post(leader_identified + i);
    }

    GFP_Standard leader = fingerprints[icentre];

    int next_leader = -1;

    int n = thread_info[0].istop;
    for (int i = 0; i < n; i++)
    {
      if (selected[i])
        continue;

      float d = leader.tanimoto_distance(fingerprints[i]);
      if (d <= threshold)
      {
        selected[i] = cluster_id;
        distances[i] = d;
      }
      else if (-1 == next_leader)
        next_leader = i;
    }

    thread_info[0].next_leader = next_leader;
    thread_info[0].next_leader_index = next_leader;

    if (reallocate_unselected_items > 0 && 0 == cluster_id % reallocate_unselected_items)
      do_reallocate_unselected_items(thread_info, nthreads, global_selected_array, global_distance_array);

    cluster_id++;

    for (int i = 1; i < nthreads; i++)
    {
      sem_wait(&cluster_formed);
    }
  }

#ifdef DEBUG_LEADER1
  for (int i = 0; i < pool_size; i++)
  {
    cerr << " i = " << i << " selected " << selected[i] << endl;
  }
#endif

  for (int i = 1; i < nthreads; i++)
  {
    thread_info[i].cluster_number = -1;
    sem_post(leader_identified + i);
  }

  for (int i = 1; i < nthreads; i++)
  {
    sem_wait(&cluster_formed);
  }

  destroy_array_of_semaphores(leader_identified, nthreads);
  sem_destroy(&cluster_formed);

  return cluster_id;
}

static int
do_previously_selected (const IW_General_Fingerprint & gfp,
                        float threshold)
{
#ifdef FIX_THIS
  GFP_Standard sgfp;

  sgfp.build_molecular_properties(gfp.molecular_properties_integer());
  sgfp.build_mk(gfp[1]);
  sgfp.build_mk2(gfp[2]);
  sgfp.build_iw(gfp[0]);

  for (int i = 0; i < pool_size; i++)
  {
    if (selected[i])
      continue;

    float d = 1.0f - sgfp.tanimoto(fingerprints[i]);

    if (d > threshold)
      continue;

    selected[i] = -1;
  }

#endif
  return 1;
}

static int
do_previously_selected (iwstring_data_source & input,
                        float threshold)
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    IW_General_Fingerprint gfp;

    int fatal;
    if (! gfp.construct_from_tdt(tdt, fatal))
    {
      cerr << "Cannot read previously selected fingerprint\n";
      return 0;
    }

    do_previously_selected(gfp, threshold);
  }

  return 1;
}

static int
do_previously_selected (const char * fname,
                        float threshold)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected(input, threshold);
}

static int
read_pool (iwstring_data_source & input)
{
  IW_TDT tdt;

  IWString tmp;

  int ndx = 0;

  for (;tdt.next(input); ndx++)
  {
    leader_item[ndx].set_gfp(fingerprints+ndx);

    tdt.dataitem_value(smiles_tag, tmp);
    leader_item[ndx].set_smiles (tmp);

    tdt.dataitem_value(identifier_tag, tmp);
    leader_item[ndx].set_id (tmp);

    IW_General_Fingerprint gfp;

    int fatal;
    if (! gfp.construct_from_tdt(tdt, fatal))
    {
      cerr << "Cannot read fingerprint\n";
      return 0;
    }

    fingerprints[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    fingerprints[ndx].build_mk(gfp[1]);
    fingerprints[ndx].build_mk2(gfp[2]);
    fingerprints[ndx].build_iw(gfp[0]);
  }

  if (ndx < 2)
  {
    cerr << "Yipes, did not read enough fingerprints\n";
    return 0;
  }

  pool_size = ndx;

  return ndx;
}

static int
allocate_pool (int s)
{
  leader_item = new Leader_Item[s];
  fingerprints = new GFP_Standard[s];

  if (NULL == leader_item || NULL == fingerprints)
  {
    cerr << "Yipes, could not allocate " << pool_size << " fingerprints\n";
    return 0;
  }

  pool_size = s;

  return 1;
}
static int
initialise_threads (int nthreads)
{
  assert (nthreads > 0);

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1)
  {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }
  
  thread_info = new Thread_Info[nthreads];

  child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init (&attr);
  pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

  thread_info[0].istart = 0;
  thread_info[0].istop = items_per_chunk;
  thread_info[0].selected = new_int(items_per_chunk);
  thread_info[0].distance = new_float(items_per_chunk, 1.0f);

  for (int i = 1; i < nthreads; i++)
  {
    Thread_Info & tii = thread_info[i];

    tii.thread_number = i;
    tii.fp = fingerprints;
    tii.istart = i * items_per_chunk;

    if (i == nthreads - 1)
      tii.istop = pool_size;
    else
      tii.istop = tii.istart + items_per_chunk;

    tii.threshold = threshold;
    tii.selected = new_int(tii.istop - tii.istart);
    tii.distance = new_float(tii.istop - tii.istart, 1.0f);

    if (0 != pthread_create (&(child_thread[i]), &attr, form_clusters_thread, (void *)(thread_info + i)))
    {
      perror ("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
gfp_leader_standard (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vt:A:s:p:r:h:q:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('t'))
  {
    cerr << "Must specify cluster threshold via the -t option\n";
    usage(2);
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', threshold) || threshold <= 0.0 || threshold >= 1.0)
    {
      cerr << "The threshold value must be a valid distance\n";
      usage(2);
    }

    if (verbose)
      cerr << "Threshold set to " << threshold << endl;
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      cerr << "The pool size specification (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Problem sized for " << pool_size << " fingerprints\n";
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(2);
    }
  }

  if (cl.option_present('q'))
  {
    if (! cl.value('q', reallocate_unselected_items) || reallocate_unselected_items < 100)
    {
      cerr << "The reallocate unselected items option must be a valid value\n";
      usage(2);
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    usage(2);
  }

  const char * fname = cl[0];

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 2;
  }

  if (0 == pool_size)
  {
    IWString tmp;
    tmp << '^' << identifier_tag;

    std::unique_ptr<re2::RE2> pcn;
    iwre2::RE2Reset(pcn, tmp);
    pool_size = input.grep (*pcn);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    if (! allocate_pool(pool_size))
      return 3;

    cerr << "Job automatically sized to " << pool_size << " fingerprints\n";
  }
  else
  {
    if (! allocate_pool(pool_size))
      return 3;
  }

  if (! read_pool(input))
  {
    cerr << "Cannot read fingerprints\n";
    return 0;
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', nthreads) || nthreads < 0)
    {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve number\n";
      usage(2);
    }
  }
  else
  {
    nthreads = 2;
    cerr << "Using 2 threads by default\n";
  }
   
  if (nthreads)
  {
    leader_identified = new sem_t[nthreads];
    if (! initialise_array_of_semaphores (leader_identified, nthreads))
      return 8;

    if (0 != sem_init(&cluster_formed, 0, 0))
    {
      perror("sem init error");
      return 6;
    }

    if (! initialise_threads(nthreads))
    {
      cerr << "Cannot initialise " << nthreads << " threads\n";
      return 3;
    }
  }

  if (cl.option_present('A'))
  {
    const char * fname = cl.option_value('A');

    if (! do_previously_selected (fname, threshold))
    {
      cerr << "Cannot process previously selected file '" << fname << "'\n";
      return 0;
    }
  }

  int * global_selected_array = new int[pool_size]; std::unique_ptr<int[]> free_global_selected_array(global_selected_array);
  float * global_distance_array = new float [pool_size]; std::unique_ptr<float[]> free_distances(global_distance_array);

  int clusters_formed;

  clusters_formed = leader(thread_info, nthreads, global_selected_array, global_distance_array);

  if (verbose)
    cerr << "Clustered " << pool_size << " fingerprints into " << clusters_formed << " clusters\n";

  int ndx = 0;
  for (int i = 0; i < nthreads; i++)
  {
    const Thread_Info & ti = thread_info[i];
    int n = ti.istop - ti.istart;

    const int * s = ti.selected;
    const float * f = ti.distance;

    for (int j = 0; j < n; j++, ndx++)
    {
      global_distance_array[ndx] = f[j];
      global_selected_array[ndx] = s[j];
    }
  }

  int start = 0;
  for (int c = 1; c <= clusters_formed; c++)
  {
    int items_in_cluster = 0;

    for (int j = start; j < pool_size; j++)
    {
      if (global_selected_array[j] != c)
        continue;

      if (0  == items_in_cluster)
        start = j;

      items_in_cluster++;
    }

    cout << smiles_tag << leader_item[start].smiles() << ">\n";
    cout << identifier_tag << leader_item[start].id() << ">\n";
    cout << "CLUSTER<" << (c-1) << ">\n";
    cout << "CSIZE<" << items_in_cluster << ">\n";

    start++;

    if (start == pool_size)
    {
      cout << "|\n";
      break;
    }

    for (int j = start; j < pool_size; j++)
    {
      if (c != global_selected_array[j])
        continue;

      cout << smiles_tag << leader_item[j].smiles() << ">\n";
      cout << identifier_tag << leader_item[j].id() << ">\n";
      cout << distance_tag << global_selected_array[j] << ">\n";
    }

    cout << "|\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_leader_standard(argc, argv);

  return rc;
}
