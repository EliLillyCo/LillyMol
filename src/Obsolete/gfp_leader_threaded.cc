/*
  An implementation of the leader algorithm for fingerprints
  This variant produces output that can be processed by nplotnn
*/

#include <stdlib.h>
#include <iostream>
#include <limits>
#include <semaphore.h>
#include <pthread.h>

#include "tbb/scalable_allocator.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Foundational/iwqsort/iwqsort.h"

#define LEADER_PARALLEL_IMPLEMENTATION

#include "Utilities/GFP_Tools/tversky.h"
#include "leader_parallel.h"

using std::cerr;
using std::endl;

static Tversky tversky;

/*
  We can have the threshold for each item read from the file
*/

static IWString score_tag;

static int score_column = -1;

/*
  May 99. When looking for the next molecule to select, we can
  add a multiple of the distance to the nearest cluster centre
  to the score. That way, we can include some varying function
  of diversity into the next cluster selection
*/

static similarity_type_t cluster_distance_scale_factor = 0.0;

static int verbose = 0;

static similarity_type_t abandon_distance_cutoff = -1.0;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = std::numeric_limits<int>::max();
static int clusters_found = 0;
static int max_cluster_size = 0;

static similarity_type_t threshold = 0.0;

static int threshold_column = -1;

/*
  when dealing with clusters which are decided by the threshold, we can
  optinally sort the cluster members by their distance from the leader
*/

static int sort_by_distance_from_centre = 0;

static int items_selected = 0;

extending_resizable_array<int> cluster_size;

static Accumulator<float> distance_stats;

static int reallocate_every = 1000;

/*
  The identifier tag used in each TDT
*/

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");

static IWString distance_tag("DIST<");

static IWString threshold_from_file_tag;

static IWString max_cluster_size_tag;

static int max_cluster_size_column = -1;

static int window_specification_made = 0;

GFP_PL::GFP_PL()
{
  _selected = 0;

  _score = 0.0;

  _shortest_distance_to_cluster_centre = static_cast<similarity_type_t>(1.0);

  return;
}

int
GFP_PL::construct_from_tdt (IW_TDT & tdt, int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal))
  {
    return 0;
  }

  if (threshold_from_file_tag.length())
  {
    similarity_type_t tmp;
    if (! tdt.dataitem_value(threshold_from_file_tag, tmp) || tmp < 0.0)
    {
      cerr << "GFP_PL::construct_from_tdt: invalid '" << threshold_from_file_tag << "' in tdt\n";
      return 0;
    }

    _threshold.set(tmp);
  }
  else if (threshold_column >= 0)
  {
    const_IWSubstring t;
    if (! _id.word(threshold_column, t))
    {
      cerr << "GFP_PL::construct_from_tdt: no " << threshold_column << " column in '" << _id << "'";
      if (::threshold > 0.0)
        cerr << " using default threshold\n";
      else
      {
        cerr << endl;
        return 0;
      }
    }
    else
    {
      similarity_type_t d;
      if (! t.numeric_value(d) || d < 0.0 || d > 1.0)
      {
        cerr << "Invalid threshold '" << t << "' in '" << _id << "'\n";
        return 0;
      }

      _threshold.set(d);
    }
  }

  if (score_tag.length())
  {
    if (! tdt.dataitem_value(score_tag, _score))
    {
      cerr << "GFP_PL::construct_from_tdt: cannot extract '" << score_tag << "' from tdt\n";
      return 0;
    }
  }

  if (score_column >= 0)
  {
    const_IWSubstring c;
    if (! _id.word(score_column, c))
    {
      cerr << "Cannot extract column " << score_column << " from '" << _id << "'\n";
      return 0;
    }

    if (! c.numeric_value(_score))
    {
      cerr << "Invalid score, column " << score_column << " in '" << _id << "'\n";
      return 0;
    }

    if (verbose > 2)
      cerr << _id << " set score to " << _score << endl;
  }

  if (max_cluster_size_tag.length() > 0)
  {
    int tmp;
    if (! tdt.dataitem_value(max_cluster_size_tag, tmp) || tmp < 1)
    {
      cerr << "GFP_PL::construct_from_tdt: missing or invalid '" << max_cluster_size_tag << "' in tdt\n";
      return 0;
    }

    _max_cluster_size.set(tmp);
  }

  if (max_cluster_size_column >= 0)
  {
    const_IWSubstring c;
    if (! _id.word(max_cluster_size_column, c))
    {
      cerr << "Cannot extract column " << max_cluster_size_column << " from '" << _id << "'\n";
      return 0;
    }

    int tmp;
    if (! c.numeric_value(tmp) || tmp < 1)
    {
      cerr << "Invalid maximum cluster size, column " << max_cluster_size_column << " in '" << _id << "'\n";
      return 0;
    }

//  cerr << "Max cluster size for '" << _id << "' is " << tmp << endl;
    _max_cluster_size.set(tmp);
  }

  if (! tdt.dataitem_value(smiles_tag, _smiles))
  {
    cerr << "No smiles in tdt '" << _id << "'\n";
    return 0;
  }

  return 1;
}

/*
  Our pool is an array of FP objects
*/

//static resizable_array<GFP_PL *> pool;
static GFP_PL ** pool = nullptr;

static int pool_size = 0;

/*
  A cluster is a set of pointers to such objects
*/

typedef resizable_array<GFP_PL *> Cluster;

static int
write_smiles_and_id (const GFP_PL & f,
                     IWString_and_File_Descriptor & output)
{
  output << smiles_tag << f.smiles() << ">\n";
  output << identifier_tag << f.id() << ">\n";

  return 1;
}

static int
write_cluster_data (const GFP_PL & fp,
                    int clusters_found,
                    int id_within_cluster, 
                    IWString_and_File_Descriptor & output)
{
  write_smiles_and_id(fp, output);

  if (fp.distance() >= 0.0)
    output << distance_tag << fp.distance() << ">\n";

  output.write_if_buffer_holds_more_than(32768);

  if (! output.good())
    throw "Bad output stream";

  return output.good();
}

int
distance_comparitor (GFP_PL * const * ppfp1, GFP_PL * const * ppfp2)
{
  const GFP_PL * pfp1 = *ppfp1;
  const GFP_PL * pfp2 = *ppfp2;

  if (pfp1->distance() < pfp2->distance())
    return -1;
  else if (pfp1->distance() > pfp2->distance())
    return 1;
  else
    return 0;
}

class Distance_Comparitor
{
  private:
  public:
    int operator () (GFP_PL * const, GFP_PL * const);
};

int
Distance_Comparitor::operator () (GFP_PL * const p1, GFP_PL * const p2)
{
  if (p1->distance() < p2->distance())
    return -1;

  if (p1->distance() > p2->distance())
    return 1;

  return 0;
}

template void resizable_array_base<GFP_PL*>::iwqsort<Distance_Comparitor>(Distance_Comparitor&);
template void iwqsort<GFP_PL*, Distance_Comparitor>(GFP_PL**, int, Distance_Comparitor&);
template void iwqsort<GFP_PL*, Distance_Comparitor>(GFP_PL**, int, Distance_Comparitor&, void*);
template void compare_two_items<GFP_PL*, Distance_Comparitor>(GFP_PL**, Distance_Comparitor&, void*);
template void move_in_from_left<GFP_PL*, Distance_Comparitor>(GFP_PL**, int&, int&, int, Distance_Comparitor&, void*);
//template void move_in_from_right<GFP_PL, Distance_Comparitor>(GFP_PL**, int&, int&, Distance_Comparitor&);
template void swap_elements<GFP_PL*>(GFP_PL*&, GFP_PL*&, void*);
template void move_in_from_right<GFP_PL*, Distance_Comparitor>(GFP_PL**, int&, int&, Distance_Comparitor&);

static int
process_cluster (Cluster & cluster,
                 IWString_and_File_Descriptor & output)
{
  if (sort_by_distance_from_centre)
  {
    Distance_Comparitor dc;
    cluster.iwqsort(dc);
  }

  int cs = cluster.number_elements();
  cluster_size[cs]++;      // the leader isn't in the cluster

  if (0 == cs)
  {
    cerr << "Very strange, zero size cluster, ignoring\n";
    return 1;
  }

  GFP_PL * centre = cluster[0];

  if (verbose)
  {
    cerr << "Cluster " << clusters_found << ' ' << cs << " items, centre '" << cluster[0]->id() << "', ";
    if (threshold_from_file_tag.length())
    {
      similarity_type_t threshold = 0.0f;
      (void) centre->threshold(threshold);
      cerr << "threshold " << threshold << ", ";
    }
    cerr << items_selected << " items selected\n";
  }

  write_smiles_and_id(*centre, output);

  output << "CLUSTER<" << clusters_found << ">\n";
  output << "CSIZE<" << cs << ">\n";

  for (int i = 1; i < cs && output.good(); i++)    // start at 1, we've already done centre above
  {
    GFP_PL & fp = *(cluster[i]);

    if (! write_cluster_data(fp, clusters_found, i, output))
      return 0;
  }

  output << "|\n";

  return output.good();
}

static int
choose_next_centre (int & first_unselected, int & icentre)
{
//cerr << "Choosing next cluster centre, from " << first_unselected << " among " << pool_size << endl;
  icentre = -1;

  int next_first_unselected = -1;

  if (0 == score_tag.length() && score_column < 0)    // just grab the first unselected item
  {
    for (int i = first_unselected; i < pool_size; i++)
    {
      if (! pool[i]->selected())
      {
        if (next_first_unselected < 0)
          next_first_unselected = i;

        icentre = i;
        return 1;
      }
    }
  }
  else if (cluster_distance_scale_factor > static_cast<similarity_type_t>(0.0))
  {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = first_unselected; i < pool_size; i++)
    {
      if (pool[i]->selected())
        continue;

      if (next_first_unselected < 0)
        next_first_unselected = i;

      score_t s = pool[i]->score() + cluster_distance_scale_factor * pool[i]->shortest_distance_to_cluster_centre();

      if (icentre < 0 || s > max_score)
      {
        max_score = s;
        icentre = i;
      }
    }
  }
  else    // raw scores
  {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = first_unselected; i < pool_size; i++)
    {
      if (pool[i]->selected())
        continue;

      if (next_first_unselected < 0)
        next_first_unselected = i;

      score_t s = pool[i]->score();

      if (icentre < 0 || s > max_score)
      {
        max_score = s;
        icentre = i;
      }
    }
  }

  first_unselected = next_first_unselected;

  return icentre >= 0;
}

static int
compute_the_distance (IW_General_Fingerprint & fp,
                      IW_General_Fingerprint & p,
                      similarity_type_t & d)
{
  if (! window_specification_made)
    ;
  else if (! can_be_compared(fp, p))
    return 0;
   
  if (tversky.active())
  {
    d = static_cast<similarity_type_t>(1.0) - fp.IW_General_Fingerprint::tversky(p, tversky);
    return 1;
  }

  if (abandon_distance_cutoff > static_cast<similarity_type_t>(0.0))
  {
    if (! fp.IW_General_Fingerprint::tanimoto(p, abandon_distance_cutoff, d))
      return 0;

    d = static_cast<similarity_type_t>(1.0) - d;
    return 1;
  }

  d = fp.IW_General_Fingerprint::distance(p);

  return 1;
}

struct Per_Thread_Data
{
  int thread_number;
  GFP_PL ** pstart;
  GFP_PL ** pstop;
  int leader;
  similarity_type_t threshold;
  Cluster sel;
  Accumulator_Int<int> items_processed;
  int cluster_number;
};

static int nthreads = 0;

/*
  First a series of parallel constructs for dealing with the actual
  leader algorithm.
*/

static struct Per_Thread_Data * ptd = nullptr;

static sem_t * start_forming_next_cluster = nullptr;
static sem_t next_cluster_found;

static pthread_t * child_thread = nullptr;

/*
  Another set of parallel constructs for dealing with the previously
  selected file
*/

struct PS_Thread_Data
{
  int thread_number;
  int istart;
  int istop;
  similarity_type_t threshold;
  IW_General_Fingerprint * previously_selected_fp;
  int discarded_this_thread;
};

static PS_Thread_Data * pstd = nullptr;

static sem_t * start_processing_previously_selected = nullptr;
static sem_t finished_processing_previously_selected;

static int previously_selected_fingerprints_processed = 0;
static int report_previously_selected = 0;

static pthread_t * ps_child_thread = nullptr;

static int
form_next_cluster (Per_Thread_Data * ptd)
{
  if (NULL == ptd->pstart)  // all our items already selected
    return 1;

  GFP_PL ** first_unselected_this_chunk = nullptr;

  GFP_PL & ldr = *(pool[ptd->leader]);
//IW_General_Fingerprint & tmp = *(pool[ptd->leader]);

//IW_General_Fingerprint * leader = new IW_General_Fingerprint(tmp); std::unique_ptr<IW_General_Fingerprint> free_leader(leader);

  ptd->sel.resize_keep_storage(0);

  similarity_type_t d;   // scope here for efficiency

  int tested = 0;

  for (GFP_PL ** pfp = ptd->pstart; pfp != ptd->pstop; pfp++)
  {
    GFP_PL * fp = *pfp;

    if (fp->selected())
      continue;

    if (NULL == first_unselected_this_chunk)
      first_unselected_this_chunk = pfp;

    tested++;

    if (! compute_the_distance(ldr, *(fp), d))
      continue;

    if (d > ptd->threshold)
      continue;

    fp->set_selected(ptd->cluster_number);
    fp->set_distance(d);
    ptd->sel.add(fp);
  }

  ptd->pstart = first_unselected_this_chunk;
  ptd->items_processed.extra(tested);

  return 1;
}

static void * 
thread_process (void * v)
{
  Per_Thread_Data * ptd = reinterpret_cast<Per_Thread_Data *>(v);

  int mythread = ptd->thread_number;

  while (1)
  {
    sem_wait(&(start_forming_next_cluster[mythread]));

    if (ptd->thread_number < 0)   // done
    {
      sem_post(&next_cluster_found);
      pthread_exit(NULL);
      return NULL;
    }

    form_next_cluster(ptd);

    sem_post(&next_cluster_found);
  }
}

static void
ps_thread_process (PS_Thread_Data * pstd)
{
  IW_General_Fingerprint * fp = pstd->previously_selected_fp;

  similarity_type_t threshold = pstd->threshold;

  pstd->discarded_this_thread = 0;

//cerr << "Thread going from " << pstd->istart << " to " << pstd->istop << " pool " << pool_size << endl;

  for (int i = pstd->istart; i < pstd->istop; i++)
  {
    if (pool[i]->selected())
      continue;

    similarity_type_t d;
    if (! compute_the_distance(*fp, *(pool[i]), d))
      continue;

    if (d > threshold)
      continue;

    pool[i]->set_selected(1);
    pool[i]->set_distance(d);

    pstd->discarded_this_thread++;
  }

  return;
}

static void *
ps_thread_process (void * v)
{
  PS_Thread_Data * pstd = reinterpret_cast<PS_Thread_Data *>(v);

  int mythread = pstd->thread_number;

  while (1)
  {
    sem_wait(&(start_processing_previously_selected[mythread]));

    if (pstd->thread_number < 0)  // done
    {
      sem_post(&finished_processing_previously_selected);
      pthread_exit(NULL);
      return NULL;
    }

    ps_thread_process(pstd);

    sem_post(&finished_processing_previously_selected);
  }
}

static int
initialise_threads (int nthreads)
{
  assert (nthreads > 0);

  ptd = new Per_Thread_Data[nthreads];

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1)
  {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++)
  {
    Per_Thread_Data & ptdi = ptd[i];

    ptdi.thread_number = i;
    ptdi.pstart = pool + i * items_per_chunk;

    if (i == nthreads - 1)
    {
      ptdi.pstop = pool + pool_size;
    }
    else
    {
      ptdi.pstop = ptdi.pstart + items_per_chunk;
    }

    ptdi.sel.resize(items_per_chunk);
    ptdi.cluster_number = 1;

    if (0 != pthread_create(&(child_thread[i]), &attr, thread_process, &ptdi))
    {
      perror("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
initialise_ps_child_threads (int nthreads)
{
  assert (nthreads > 0);

  pstd = new PS_Thread_Data[nthreads];

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1)
  {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  ps_child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++)
  {
    PS_Thread_Data & pstdi = pstd[i];

    pstdi.thread_number = i;
    pstdi.istart = i * items_per_chunk;
    if (i == nthreads - 1)
      pstdi.istop = pool_size;
    else
      pstdi.istop = pstdi.istart + items_per_chunk;

    if (0 != pthread_create(&(ps_child_thread[i]), &attr, ps_thread_process, &pstdi))
    {
      perror("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
reallocate_unselected_items ()
{
  if (verbose > 1)
    cerr << "Reallocating from " << pool_size << " items\n";

  int unselected_items = 0;
  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i]->selected())
      ;      // initially I had delete pool[i] but that breaks things
    else
    {
      pool[unselected_items] = pool[i];
      unselected_items++;
    }
  }

  pool_size = unselected_items;

/*for (int i = 0; i < pool_size; i++)
  {
    cerr << " i = " << i <<  " fo " << pool[i]->id() << "'\n";
  }*/

  int items_per_chunk = unselected_items / nthreads;

#ifdef DEBUG_REALLOCATE_UNSELECTED
  cerr << "Means " << items_per_chunk << " items per chunk\n";
#endif

  for (int i = 0; i < nthreads; i++)
  {
    ptd[i].pstart = pool + i * items_per_chunk;

    if (i == nthreads - 1)   // last thread
    {
      ptd[i].pstop = pool + pool_size;
    }
    else
    {
      ptd[i].pstop = ptd[i].pstart + items_per_chunk;
    }
  }

#ifdef ECHO_THREAD_LOCATIONS
  for (int i = 0; i < nthreads; i++)
  {
    cerr << " thread " << i << " starts " << (ptd[i].pstart - pool) << " stop " << (ptd[i].pstop - pool) << endl;
  }
#endif

  return 1;
}

int
leader (IWString_and_File_Descriptor & output)
{
  assert (0 == items_selected);
  assert (0 == clusters_found);

  int first_unselected = 0;

  int icentre;
  if (! choose_next_centre(first_unselected, icentre))
  {
    cerr << "Yipes, cannot find initial leader\n";
    return 0;
  }

  Cluster cluster;
  cluster.resize(pool_size);

  if (! initialise_threads(nthreads))
    return 0;

  int next_reallocate = reallocate_every;

  int initial_pool_size = pool_size;

  while (items_selected < initial_pool_size)
  {
    similarity_type_t my_threshold;
    if (pool[icentre]->threshold(my_threshold))     // has come from the file
      ;
    else
      my_threshold = threshold;

    if (verbose > 1)
      cerr << "Start cluster " << clusters_found << ". ndx " << icentre << ", threshold = " << my_threshold << endl;

    pool[icentre]->set_selected(clusters_found + 1);
    pool[icentre]->set_distance(static_cast<similarity_type_t>(0.0));

    for (int i = 0; i < nthreads; i++)
    {
      ptd[i].leader = icentre;
      ptd[i].threshold = my_threshold;
      ptd[i].cluster_number = clusters_found + 1;
      sem_post(&start_forming_next_cluster[i]);
    }

//  While the threads are forming the next cluster, process the cluster
//  we just built

    if (cluster.size() > 0)
    {
      items_selected += cluster.number_elements();

      (void) process_cluster(cluster, output);

      cluster.resize_keep_storage(0);

      clusters_found++;
      if (clusters_found >= max_clusters_to_find)
        break;
    }

    cluster.add(pool[icentre]);

    for (int i = 0; i < nthreads; i++)
    {
      sem_wait(&next_cluster_found);
    }

    for (int i = 0; i < nthreads; i++)
    {
      cluster += ptd[i].sel;
    }

    if (items_selected > next_reallocate)
    {
      reallocate_unselected_items();
      next_reallocate += reallocate_every;
      first_unselected = 0;
    }

    if (! choose_next_centre(first_unselected, icentre))
      break;
  }

  if (cluster.size() > 0)
  {
    items_selected += cluster.number_elements();

    (void) process_cluster(cluster, output);

    clusters_found++;
  }

  return 1;
}

/*
  If we have a previously selected file, we keep track of the number
  of members of the pool that get selected by the previously selected file
*/

static int molecules_selected_by_previously_selected_file = 0;

static int
do_previously_selected_file (iwstring_data_source & input,
                             similarity_type_t t)
{
  IW_TDT tdt;

  if (! tdt.next(input))
  {
    cerr << "No first TDT in previously selected file\n";
    return 0;
  }

  IW_General_Fingerprint fp1;
  int fatal;
  if (! fp1.construct_from_tdt(tdt, fatal))   // make things easy
    return 0;

  for (int i = 0; i < nthreads; i++)
  {
    pstd[i].threshold = t;
    pstd[i].previously_selected_fp = &fp1;
    sem_post(&(start_processing_previously_selected[i]));
  }

  IW_General_Fingerprint fp2;

  int ndx = 0;

// We use two fingerprints. While one is off being compared, the other
// is being built.

  while (tdt.next(input))
  {
    IW_General_Fingerprint * fp;

    if (0 == ndx)
    {
      fp = &fp2;
      ndx = 1;
    }
    else
    {
      fp = &fp1;
      ndx = 0;
    }

    if (! fp->construct_from_tdt(tdt, fatal))
      return 0;

//  cerr << "Built fp '" << fp->id() << "'\n";

    for (int i = 0; i < nthreads; i++)
    {
      sem_wait(&finished_processing_previously_selected);
    }

    for (int i = 0; i < nthreads; i++)
    {
      molecules_selected_by_previously_selected_file += pstd[i].discarded_this_thread;

      pstd[i].previously_selected_fp = fp;
      sem_post(&(start_processing_previously_selected[i]));
    }

    previously_selected_fingerprints_processed++;

    if (report_previously_selected > 0 && 0 == previously_selected_fingerprints_processed % report_previously_selected)
      cerr << "Processed " << previously_selected_fingerprints_processed << " previously selected items, selected " << molecules_selected_by_previously_selected_file <<endl;
  }


  for (int i = 0; i < nthreads; i++)
  {
    sem_wait(&finished_processing_previously_selected);
  }

  return 1;
}

static int
do_previously_selected_file (const IWString & fname,
                             similarity_type_t t)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected_file(input, t);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Performs leader clustering on a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -C <number>      maximum number of clusters to find\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -t <dis>         specify distance threshold\n";
  cerr << " -t col=<nn>      threshold is column <nn> of the name field\n";
  cerr << " -t tag=<TAG>     threshold for each molecule in dataitem <TAG>\n";
  cerr << " -H <TAG>         threshold for each molecule in dataitem <TAG>\n";
  cerr << " -m <number>      maximum cluster size\n";
  cerr << " -M <tag>         max cluster size for each molecule in <TAG>\n";
  cerr << " -M col=nn        max cluster size is column <nn> of the name field\n";
  cerr << " -S <TAG>         score tag\n";
  cerr << " -S col=nn        score is column <nn> in the name field\n";
  cerr << " -I <TAG>         specify identifier tag\n";
  cerr << " -r               sort clusters by distance from leader\n";
//cerr << " -X <distance>    abandon distance computation if any component > distance\n";
//cerr << " -R <factor>      score = score + <factor> * distance to cluster\n";
  cerr << " -A <file>        file(s) of previously selected molecules - discard all within threshold\n";
  cerr << " -a <dist>        use <dist> as the threshold when comparing against the -A file\n";
  cerr << " -p <nn>          report progress when processing the -A file every <nn> items\n";
//cerr << " -D ...           miscellaneous options, enter '-D help' for info\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -h <number>      number of threads to use\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
initialise_array_of_semaphores (sem_t * sem,
                                int nsem)
{
  assert(nsem > 0);

  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_init(&(sem[i]), 0, 0))
    {
      perror("sem_init error");
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

static int
leader (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:I:t:F:P:W:X:rH:S:C:V:R:m:M:A:a:p:h:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }
  }

  if (cl.option_present('V'))
  {
    if (! tversky.parse_command_line(cl, 'V', verbose))
    {
      cerr << "Cannot parse Tversky specifications\n";
      usage(38);
    }
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', max_clusters_to_find) || max_clusters_to_find < 1)
    {
      cerr << "The -C option (max clusters to find) must be followed by a positive integer\n";
      usage(41);
    }

    if (verbose)
      cerr << "Will find a max of " << max_clusters_to_find << " clusters\n";
  }

  if (cl.option_present('r'))
  {
    sort_by_distance_from_centre = 1;
    if (verbose)
      cerr << "Clusters sorted by distance from centre\n";
  }

// We need to be careful with the -i and -I options. Remember
// that the pool is built first

  if (cl.option_present('I'))
  {
    (void) cl.value('I', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose)
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
  }

  if (cl.option_present('W'))
  {
    window_specification_made = 1;

    if (verbose)
      cerr << "Window specification made, will check for comparability\n";
  }

  if (cl.option_present('H'))
  {
    cl.value('H', threshold_from_file_tag);
    if (verbose)
      cerr << "Each threshold from the '" << threshold_from_file_tag << "' dataitem in the input\n";

    if (! threshold_from_file_tag.ends_with('<'))
      threshold_from_file_tag << '<';
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (s.starts_with("col="))
    {
      s.remove_leading_chars(4);
      if (! s.numeric_value(score_column) || score_column < 1)
      {
        cerr << "Invalid column for score '" << s << "'\n";
        usage(14);
      }

      if (verbose)
        cerr << "Score for each item in column " << score_column << endl;

      score_column--;
    }
    else
    {
      score_tag = s;
      if (verbose)
        cerr << "Score tag is " << score_tag << "'\n";
    }
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', cluster_distance_scale_factor) || cluster_distance_scale_factor <= 0.0)
    {
      cerr << "The cluster distance scale factor option (-R) must be followed by a positive number\n";
      usage(19);
    }

    if (verbose)
      cerr << "Scores adjusted by " << cluster_distance_scale_factor << " times distance to nearest cluster centre\n";
  }

// Because we are reading the file(s) in parallel, we need to always
// initialise the fingerprints here

  if (! initialise_fingerprints(cl, verbose))
  {
    cerr << "Cannot initialise general fingerprint options\n";
    usage(17);
  }

  if (cl.option_present('X'))
  {
    if (! cl.value('X', abandon_distance_cutoff) || abandon_distance_cutoff < 0.0 || abandon_distance_cutoff > 1.0)
    {
      cerr << "The -X option must be followed by a valid distance (0.0, 1.0)\n";
      usage(13);
    }

    if (verbose)
      cerr << "Distance compuations abandoned if any component > " << abandon_distance_cutoff << endl;
  }

  if (! cl.option_present('t') && ! cl.option_present('H') && ! cl.option_present('Y'))
  {
    cerr << "Threshold distance must be specified via -t, -H or -V options\n";
    usage(28);
  }

  if (cl.option_present('t'))
  {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++))
    {
      if (t.starts_with("col="))
      {
        t.remove_leading_chars(4);
        if (! t.numeric_value(threshold_column) || threshold_column < 1)
        {
          cerr << "Invalid column for threshold '" << t << "'\n";
          usage(14);
        }
  
        if (verbose)
          cerr << "Threshold for each item in column " << threshold_column << endl;

        threshold_column--;
      }
      else if (t.starts_with("tag="))
      {
        threshold_from_file_tag = t;
        threshold_from_file_tag.remove_leading_chars(4);
        if (verbose)
          cerr << "Threshold in tag " << threshold_from_file_tag << "'\n";

        if (! threshold_from_file_tag.ends_with('<'))
          threshold_from_file_tag.add('<');
      }
      else if (! t.numeric_value(threshold) || threshold < 0.0 || threshold > 1.0)
      {
        cerr << "The -t option must be followed by a valid distance value\n";
        usage(12);

        if (verbose)
          cerr << "Distance threshold set to " << threshold << endl;
      }
    }
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', max_cluster_size) || max_cluster_size < 2)
    {
      cerr << "The -m (max cluster size) option must be followed by a whole number > 1\n";
      usage(43);
    }

    if (verbose)
      cerr << "Max cluster size " << max_cluster_size << endl;
  }

  if (cl.option_present('M'))
  {
    const_IWSubstring m = cl.string_value('M');

    if (m.starts_with("col="))
    {
      m.remove_leading_chars(4);
      if (! m.numeric_value(max_cluster_size_column) || max_cluster_size_column < 1)
      {
        cerr << "The column for the per molecule maximum cluster size must be a whole positive number\n";
        usage(11);
      }

      if (verbose)
        cerr << "The maximum cluster size per molecule will be in column " << max_cluster_size_column << endl;

      max_cluster_size_column--;
    }
    else
    {
      max_cluster_size_tag = m;

      if (verbose)
        cerr << "Max cluster size in '" << max_cluster_size_tag << "' tag\n";
    }
  }

  if (! cl.option_present('h'))
  {
    cerr << "Must specify number of threads via the -h option\n";
    usage(3);
  }

  Pool_Formation_Info<GFP_PL> pfi(verbose);

  pfi.set_number_single_file_parallel_readers(5);

  if (0 == pool_size)
    ;
  else if (pfi.allocate_pool("cl", pool_size))
    ;
  else
  {
    cerr << "Cannot size fingerprint pool (-s option)\n";
    return 8;
  }

  if (! pfi.build(cl))
  {
    cerr << "Cannot read fingerprints\n";
    return 4;
  }

  pool_size = pfi.pool_size();
  pool = pfi.pool();

  if (verbose)
    cerr << "Read " << pool_size << " fingerprints\n";

  if (pool_size < 1)
  {
    cerr << "No fingerprints\n";
    return 5;
  }

  int initial_fingerprint_count = pool_size;

  if (cl.option_present('h'))
  {
    if (! cl.value('h', nthreads) || nthreads < 1)
    {
      cerr << "The number of threads to use (-h) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will use " << nthreads << " threads\n";

    start_forming_next_cluster = new sem_t[nthreads];
    if (! initialise_array_of_semaphores(start_forming_next_cluster, nthreads))
      return 8;

    if (0 != sem_init(&next_cluster_found, 0, 0))
    {
      perror("sem init error");
      return 6;
    }
  }

  if (cl.option_present('A'))
  {
    if (! cl.option_present('t') && ! cl.option_present('a'))
    {
      cerr << "Must have a threshold available with the -A option (use -t or -a)\n";
      usage(11);
    }

    similarity_type_t t;

    if (cl.option_present('a'))
    {
      if (! cl.value('a', t) || t < 0.0 || t >= 1.0)
      {
        cerr << "Invalid value for previously selected threshold (-a option)\n";
        usage(4);
      }

      if (verbose)
        cerr << "Will use " << t << " as the threshold for the previously selected list\n";
    }
    else
      t = threshold;

    if (cl.option_present('p'))
    {
      if (! cl.value('p', report_previously_selected) || report_previously_selected < 1)
      {
        cerr << "The report scanning of previously selected (-p) option must be a whole +ve number\n";
        usage(4);
      }

      if (verbose)
        cerr << "Will report progress against previously selected (-A) every " << report_previously_selected << " items\n";
    }

    start_processing_previously_selected = new sem_t[nthreads];
    if (! initialise_array_of_semaphores(start_processing_previously_selected, nthreads))
      return 5;

    if (! initialise_ps_child_threads(nthreads))
      return 8;

    IWString fname;
    int i = 0;
    while (cl.value('A', fname, i++))
    {
      if (! do_previously_selected_file(fname, t))
      {
        cerr << "Cannot process previously selected file (-A option)\n";
        return 8;
      }
    }

    if (verbose)
      cerr << "Rejected " << molecules_selected_by_previously_selected_file << " molecules by previously selected file(s)\n";

    if (molecules_selected_by_previously_selected_file == pool_size)
    {
      cerr << "Yipes, the previously selected file knocked out the whole pool\n";
      return 1;
    }

//  Clean up all the threads involved

    for (int i = 0; i < nthreads; i++)
    {
      pstd[i].thread_number = -1;
      pstd[i].previously_selected_fp = nullptr;
      sem_post(&(start_processing_previously_selected[i]));
    }

    for (int i = 0; i < nthreads; i++)
    {
      sem_wait(&finished_processing_previously_selected);
    }

    destroy_array_of_semaphores(start_processing_previously_selected, nthreads);

    delete [] ps_child_thread;
    sem_destroy(&finished_processing_previously_selected);

    if (verbose)
      cerr << "Processed " << previously_selected_fingerprints_processed << " previously selected fingerprints, selected " << molecules_selected_by_previously_selected_file << endl;
  }

  IWString_and_File_Descriptor output(1);

  try {
    leader(output);
  }
  catch (const char * err)
  {
    cerr << "Caught '" << err << "' terminated\n";
    return 81;
  }

  output.flush();

  if (verbose)
  {
    cerr << "Clustered " << initial_fingerprint_count << " fingerprints into " << clusters_found << " clusters\n";
    int isum = 0;
    for (int i = 0; i < cluster_size.number_elements(); i++)
    {
      int j = cluster_size[i];
      if (0 == j)
        continue;

      cerr << j << " clusters were of size " << i << " members\n";

      isum += j * i;
    }

    cerr << "In clusters " << isum << endl;
  }

// Tell all the threads to destroy themselves

  if (nthreads)
  {
    for (int i = 0; i < nthreads; i++)
    {
      ptd[i].thread_number = -1;

      sem_post(&(start_forming_next_cluster[i]));
    }

    for (int i = 0; i < nthreads; i++)
    {
      sem_wait(&next_cluster_found);
    }

    if (0 != sem_destroy(&next_cluster_found))
      cerr << "Could not destroy semaphores\n";

    destroy_array_of_semaphores(start_forming_next_cluster, nthreads);
  }

  delete [] start_forming_next_cluster;

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = leader(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
