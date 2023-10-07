/*
  An implementation of the leader algorithm for fingerprints
  Runs parallel with MPI.
*/

#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include <values.h>

#include <iostream>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Utilities/GFP_Tools/leader.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;

// #undef SEEK_SET
// #undef SEEK_END
// #undef SEEK_CUR
#include "mpi.h"

static int my_mpi_id = 0;
static int my_mpi_size = 0;

typedef struct {
  off_t offset;
  double dis;
} MY_MPI_STRUCT;

static MPI_Datatype my_mpi_type;

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

static resizable_array_p<IWString> dataitems_to_echo;

static int verbose = 0;

static similarity_type_t abandon_distance_cutoff = -1.0;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = std::numeric_limits<int>::max();
static int my_mpi_clusters_found = 0;
static int max_cluster_size = 0;

static similarity_type_t threshold = 0.0;

static int threshold_column = -1;

/*
  when dealing with clusters which are decided by the threshold, we can
  optinally sort the cluster members by their distance from the leader
*/

static int sort_by_distance_from_centre = 0;

static int items_selected = 0;
static int my_mpi_items_selected = 0;

extending_resizable_array<int> cluster_size;

static Accumulator<float> distance_stats;

static int reallocate_every = 500;

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag("PCN<");

static IWString distance_tag("DIST<");

GFP_L::GFP_L()
{
  _selected = 0;

  _score = 0.0;

  _shortest_distance_to_cluster_centre = static_cast<similarity_type_t>(1.0);

  return;
}

static IWString threshold_from_file_tag;

static IWString max_cluster_size_tag;

static int max_cluster_size_column = -1;

int
GFP_L::construct_from_tdt(IW_TDT& tdt, int& fatal)
{
  if (!IW_General_Fingerprint::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  if (threshold_from_file_tag.length()) {
    similarity_type_t tmp;
    if (!tdt.dataitem_value(threshold_from_file_tag, tmp) || tmp < 0.0) {
      cerr << "GFP_L::construct_from_tdt: invalid '" << threshold_from_file_tag
           << "' in tdt\n";
      return 0;
    }

    _threshold.set(tmp);
  } else if (threshold_column >= 0) {
    const_IWSubstring t;
    if (!_id.word(threshold_column, t)) {
      cerr << "GFP_L::construct_from_tdt: no " << threshold_column << " column in '"
           << _id << "'";
      if (::threshold > 0.0) {
        cerr << " using default threshold\n";
      } else {
        cerr << endl;
        return 0;
      }
    } else {
      similarity_type_t d;
      if (!t.numeric_value(d) || d < 0.0 || d > 1.0) {
        cerr << "Invalid threshold '" << t << "' in '" << _id << "'\n";
        return 0;
      }

      _threshold.set(d);
    }
  }

  if (score_tag.length()) {
    if (!tdt.dataitem_value(score_tag, _score)) {
      cerr << "GFP_L::construct_from_tdt: cannot extract '" << score_tag
           << "' from tdt\n";
      return 0;
    }
  }

  if (score_column >= 0) {
    const_IWSubstring c;
    if (!_id.word(score_column, c)) {
      cerr << "Cannot extract column " << score_column << " from '" << _id << "'\n";
      return 0;
    }

    if (!c.numeric_value(_score)) {
      cerr << "Invalid score, column " << score_column << " in '" << _id << "'\n";
      return 0;
    }

    if (verbose > 2) {
      cerr << _id << " set score to " << _score << endl;
    }
  }

  if (max_cluster_size_tag.length() > 0) {
    int tmp;
    if (!tdt.dataitem_value(max_cluster_size_tag, tmp) || tmp < 1) {
      cerr << "GFP_L::construct_from_tdt: missing or invalid '" << max_cluster_size_tag
           << "' in tdt\n";
      return 0;
    }

    _max_cluster_size.set(tmp);
  }

  if (max_cluster_size_column >= 0) {
    const_IWSubstring c;
    if (!_id.word(max_cluster_size_column, c)) {
      cerr << "Cannot extract column " << max_cluster_size_column << " from '" << _id
           << "'\n";
      return 0;
    }

    int tmp;
    if (!c.numeric_value(tmp) || tmp < 1) {
      cerr << "Invalid maximum cluster size, column " << max_cluster_size_column
           << " in '" << _id << "'\n";
      return 0;
    }

    //  cerr << "Max cluster size for '" << _id << "' is " << tmp << endl;
    _max_cluster_size.set(tmp);
  }

  return 1;
}

/*
  Our pool is an array of FP objects
*/

static GFP_L* pool = nullptr;

static int pool_size = 0;

static int my_mpi_pool_size = 0;

static GFP_L my_mpi_centre;

/*
  A cluster is a set of pointers to such objects
*/

typedef resizable_array<GFP_L*> Cluster;

static int
echo_selected_dataitems(IW_TDT& tdt, const resizable_array_p<IWString>& items_to_echo,
                        IWString_and_File_Descriptor& output)
{
  int ne = items_to_echo.number_elements();
  for (int i = 0; i < ne; i++) {
    const IWString& tag = *(items_to_echo[i]);

    if (!tdt.echo_dataitem(tag, 0, output)) {
      cerr << "Cannot echo '" << tag << "' from TDT\n";
      throw "Missing dataitem";
      return 0;
    }
  }

  return output.good();
}

static int
echo_all_dataitems(IW_TDT& tdt, IWString_and_File_Descriptor& output)
{
  return tdt.write_all_except_vbar(output);
}

static int
write_cluster_data(IW_TDT& tdt, int clusters_found, int id_within_cluster,
                   similarity_type_t distance_to_centre,
                   IWString_and_File_Descriptor& output)
{
  if (dataitems_to_echo.number_elements()) {
    echo_selected_dataitems(tdt, dataitems_to_echo, output);
  } else {
    echo_all_dataitems(tdt, output);
  }

  if (distance_to_centre >= 0.0) {
    output << distance_tag << distance_to_centre << ">\n";
  }

  if (!output.good()) {
    throw "Bad output stream";
  }

  return output.good();
}

/*
static int
get_tdt (IW_TDT & tdt, iwstring_data_source & input,
         GFP_L & fp)
{
  off_t offset;
  (void) fp.offset (offset);

  if (! input.seekg (offset))
  {
    cerr << "Cannot seek to offset '" << offset << endl;
    return 0;
  }

  return tdt.next (input);
}*/

static int
get_tdt(IW_TDT& tdt, iwstring_data_source& input, off_t& offset)
{
  if (!input.seekg(offset)) {
    cerr << "Cannot seek to offset '" << offset << endl;
    return 0;
  }

  return tdt.next(input);
}

int
distance_comparitor(GFP_L* const* ppfp1, GFP_L* const* ppfp2)
{
  const GFP_L* pfp1 = *ppfp1;
  const GFP_L* pfp2 = *ppfp2;

  if (pfp1->distance() < pfp2->distance()) {
    return -1;
  } else if (pfp1->distance() > pfp2->distance()) {
    return 1;
  } else {
    return 0;
  }
}

class Distance_Comparitor
{
 private:
 public:
  int
  operator()(const MY_MPI_STRUCT&, const MY_MPI_STRUCT&) const;
};

int
Distance_Comparitor::operator()(const MY_MPI_STRUCT& p1, const MY_MPI_STRUCT& p2) const
{
  if (p1.dis < p2.dis) {
    return -1;
  }

  if (p1.dis > p2.dis) {
    return 1;
  }

  return 0;
}

// template void
// resizable_array_base<GFP_L*>::iwqsort<Distance_Comparitor>(Distance_Comparitor&);
template void
iwqsort<MY_MPI_STRUCT, Distance_Comparitor>(MY_MPI_STRUCT*, int, Distance_Comparitor&);
template void
iwqsort<MY_MPI_STRUCT, Distance_Comparitor>(MY_MPI_STRUCT*, int, Distance_Comparitor&,
                                            void*);
template void
compare_two_items<MY_MPI_STRUCT, Distance_Comparitor>(MY_MPI_STRUCT*,
                                                      Distance_Comparitor&, void*);
template void
move_in_from_left<MY_MPI_STRUCT, Distance_Comparitor>(MY_MPI_STRUCT*, int&, int&, int,
                                                      Distance_Comparitor&, void*);
// template void move_in_from_right<GFP_L, Distance_Comparitor>(GFP_L**, int&, int&,
// Distance_Comparitor&);
template void
swap_elements<MY_MPI_STRUCT>(MY_MPI_STRUCT&, MY_MPI_STRUCT&, void*);
template void
move_in_from_right<MY_MPI_STRUCT, Distance_Comparitor>(MY_MPI_STRUCT*, int&, int&,
                                                       Distance_Comparitor&);

static int
process_cluster(MY_MPI_STRUCT* my_mpi_cluster_all, int cs, iwstring_data_source& input,
                IWString_and_File_Descriptor& output)
{
  if (sort_by_distance_from_centre) {
    Distance_Comparitor dc;
    iwqsort(my_mpi_cluster_all, cs, dc);
    //  cluster.iwqsort (dc);
    //  cluster.sort (&distance_comparitor);
  }

  // int cs = cluster.number_elements ();
  cluster_size[cs]++;  // the leader isn't in the cluster

  // GFP_L * centre = cluster[0];

  /* if (verbose)
   {
     cerr << "Cluster " << my_mpi_clusters_found << ' ' << cs << " items, centre '" <<
   my_mpi_centre_old.id () << "', ";
    // cerr << "Cluster " << my_mpi_clusters_found << ' ' << cs << " items, centre '" <<
   my_mpi_centre_old.id () << "', "; if (threshold_from_file_tag.length ())
     {
       similarity_type_t threshold;
       (void) my_mpi_centre_old.threshold (threshold);
       cerr << "threshold " << threshold << ", ";
     }
     cerr << (my_mpi_items_selected + cs) << " items selected\n";
   }
   */
  IW_TDT tdt;
  if (!get_tdt(tdt, input, my_mpi_cluster_all[0].offset)) {
    return 0;
  }

  if (dataitems_to_echo.number_elements()) {
    echo_selected_dataitems(tdt, dataitems_to_echo, output);
  } else {
    echo_all_dataitems(tdt, output);
  }

  output << "CLUSTER<" << my_mpi_clusters_found << ">\n";
  output << "CSIZE<" << cs << ">\n";

  for (int i = 1; i < cs && output.good(); i++) {
    if (!get_tdt(tdt, input, my_mpi_cluster_all[i].offset)) {
      return 0;
    }
    if (!write_cluster_data(tdt, my_mpi_clusters_found, cs, my_mpi_cluster_all[i].dis,
                            output)) {
      return 0;
    }
  }
  output << "|\n";

  return output.good();
}

static int
choose_next_centre(int& first_unselected, int& icentre)
{
  icentre = -1;

  int next_first_unselected = -1;

  if (0 == score_tag.length() && score_column < 0)  // just grab the first unselected item
  {
    for (int i = first_unselected; i < pool_size; i++) {
      if (!pool[i].selected()) {
        if (next_first_unselected < 0) {
          next_first_unselected = i;
        }

        icentre = i;
        return 1;
      }
    }
  } else if (cluster_distance_scale_factor > static_cast<similarity_type_t>(0.0)) {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = first_unselected; i < pool_size; i++) {
      if (pool[i].selected()) {
        continue;
      }

      if (next_first_unselected < 0) {
        next_first_unselected = i;
      }

      score_t s = pool[i].score() + cluster_distance_scale_factor *
                                        pool[i].shortest_distance_to_cluster_centre();

      if (icentre < 0 || s > max_score) {
        max_score = s;
        icentre = i;
      }
    }
  } else  // raw scores
  {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = first_unselected; i < pool_size; i++) {
      if (pool[i].selected()) {
        continue;
      }

      if (next_first_unselected < 0) {
        next_first_unselected = i;
      }

      score_t s = pool[i].score();

      if (icentre < 0 || s > max_score) {
        max_score = s;
        icentre = i;
      }
    }
  }

  first_unselected = next_first_unselected;

  return icentre >= 0;
}

static int
compute_the_distance(IW_General_Fingerprint& fp, IW_General_Fingerprint& p,
                     similarity_type_t& d)
{
  if (!can_be_compared(fp, p)) {
    return 0;
  }

  if (tversky.active()) {
    d = static_cast<similarity_type_t>(1.0) -
        fp.IW_General_Fingerprint::tversky(p, tversky);
    return 1;
  }

  if (abandon_distance_cutoff > static_cast<similarity_type_t>(0.0)) {
    if (!fp.IW_General_Fingerprint::tanimoto(p, abandon_distance_cutoff, d)) {
      return 0;
    }

    d = static_cast<similarity_type_t>(1.0) - d;
    return 1;
  }

  d = fp.IW_General_Fingerprint::distance(p);

  return 1;
}

struct Per_Thread_Data {
  int thread_number;
  int istart;
  int istop;
  // int leader;
  similarity_type_t threshold;
  Cluster sel;
  Accumulator_Int<int> items_processed;
};

static int nthreads = 0;

/*
  First a series of parallel constructs for dealing with the actual
  leader algorithm.
*/

static struct Per_Thread_Data* ptd = nullptr;

static sem_t* start_forming_next_cluster = nullptr;
static sem_t next_cluster_found;

static pthread_t* child_thread = nullptr;

/*
  Another set of parallel constructs for dealing with the previously
  selected file
*/

struct PS_Thread_Data {
  int thread_number;
  int istart;
  int istop;
  similarity_type_t threshold;
  IW_General_Fingerprint* previously_selected_fp;
  int discarded_this_thread;
};

static PS_Thread_Data* pstd = nullptr;

static sem_t* start_processing_previously_selected = nullptr;
static sem_t finished_processing_previously_selected;

static pthread_t* ps_child_thread = nullptr;

static int
form_next_cluster(Per_Thread_Data* ptd)
{
  if (ptd->istart < 0) {  // all our items already selected
    return 1;
  }

  int first_unselected_this_chunk = -1;

  // GFP_L & ldr = pool[ptd->leader];
  GFP_L& ldr = my_mpi_centre;

  ptd->sel.resize_keep_storage(0);

  similarity_type_t d;  // scope here for efficiency

  int tested = 0;

  for (int i = ptd->istart; i < ptd->istop; i++) {
    if (pool[i].selected()) {
      continue;
    }

    if (first_unselected_this_chunk < 0) {
      first_unselected_this_chunk = i;
    }

    tested++;

    if (!compute_the_distance(ldr, pool[i], d)) {
      continue;
    }

    if (d > ptd->threshold) {
      continue;
    }

    pool[i].set_selected(1);
    pool[i].set_distance(d);
    ptd->sel.add(&(pool[i]));
  }

  ptd->istart = first_unselected_this_chunk;
  ptd->items_processed.extra(tested);

  return 1;
}

static void*
thread_process(void* v)
{
  Per_Thread_Data* ptd = reinterpret_cast<Per_Thread_Data*>(v);

  int mythread = ptd->thread_number;

  while (1) {
    sem_wait(&(start_forming_next_cluster[mythread]));

    if (ptd->thread_number < 0)  // done
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
ps_thread_process(PS_Thread_Data* pstd)
{
  IW_General_Fingerprint* fp = pstd->previously_selected_fp;

  similarity_type_t threshold = pstd->threshold;

  pstd->discarded_this_thread = 0;

  for (int i = pstd->istart; i < pstd->istop; i++) {
    if (pool[i].selected()) {
      continue;
    }

    similarity_type_t d;
    if (!compute_the_distance(*fp, pool[i], d)) {
      continue;
    }

    if (d > threshold) {
      continue;
    }

    pool[i].set_selected(1);
    pool[i].set_distance(d);

    pstd->discarded_this_thread++;
  }

  return;
}

static void*
ps_thread_process(void* v)
{
  PS_Thread_Data* pstd = reinterpret_cast<PS_Thread_Data*>(v);

  int mythread = pstd->thread_number;

  while (1) {
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
initialise_threads(int nthreads)
{
  assert(nthreads > 0);

  ptd = new Per_Thread_Data[nthreads];

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1) {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++) {
    Per_Thread_Data& ptdi = ptd[i];

    ptdi.thread_number = i;
    ptdi.istart = i * items_per_chunk;
    if (i == nthreads - 1) {
      ptdi.istop = pool_size;
    } else {
      ptdi.istop = ptdi.istart + items_per_chunk;
    }

    ptdi.sel.resize(items_per_chunk);

    if (0 != pthread_create(&(child_thread[i]), &attr, thread_process, &ptdi)) {
      perror("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
initialise_ps_child_threads(int nthreads)
{
  assert(nthreads > 0);

  pstd = new PS_Thread_Data[nthreads];

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1) {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  ps_child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++) {
    PS_Thread_Data& pstdi = pstd[i];

    pstdi.thread_number = i;
    pstdi.istart = i * items_per_chunk;
    if (i == nthreads - 1) {
      pstdi.istop = pool_size;
    } else {
      pstdi.istop = pstdi.istart + items_per_chunk;
    }

    if (0 != pthread_create(&(ps_child_thread[i]), &attr, ps_thread_process, &pstdi)) {
      perror("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
reallocate_unselected_items(int items_selected)
{
  int unselected_items = pool_size - items_selected;

  if (unselected_items < 1000) {  // not worth doing
    return 0;
  }

#ifdef DEBUG_REALLOCATE_UNSELECTED
  cerr << "Reallocating, " << items_selected << " items selected, unselected "
       << unselected_items << endl;
#endif

  int items_per_chunk = unselected_items / nthreads;

#ifdef DEBUG_REALLOCATE_UNSELECTED
  cerr << "Means " << items_per_chunk << " items per chunk\n";
#endif

  if (items_per_chunk < 300) {
    return 0;
  }

  int istart = ptd[0].istart;

  for (int i = 0; i < nthreads; i++) {
    ptd[i].istart = -1;
    ptd[i].istop = -1;
  }

  int ndx = 0;
  int items_this_chunk = 0;

  for (int i = istart; i < pool_size; i++) {
    if (pool[i].selected()) {
      continue;
    }

    //  cerr << " i = " << i << " ndx = " << ndx << " istart " << ptd[ndx].istart << ", "
    //  << items_this_chunk << " items in chunk, " <<ptd[0].istart << ' ' << ptd[1].istart
    //  << endl;
    if (ptd[ndx].istart < 0)  // starting new chunk
    {
      //    cerr << ndx << " starts at " << i << ", " << items_this_chunk << " items in
      //    previous chunk\n";
      ptd[ndx].istart = i;
      //    cerr << ptd[0].istart << " and " << ptd[1].istart << endl;
      items_this_chunk = 1;
    } else {
      items_this_chunk++;
      if (items_this_chunk == items_per_chunk) {
        //      cerr << ndx << " ends at " << i << " contains " << items_this_chunk <<
        //      endl;
        ptd[ndx].istop = i + 1;
        ndx++;
      }
    }
  }

  // cerr << "ndx = " << ndx << " last chunk contains " << items_this_chunk << endl;

  ptd[nthreads - 1].istop = pool_size;

  if (verbose > 2) {
    cerr << "Reallocated threads\n";
    for (int i = 0; i < nthreads; i++) {
      cerr << "Thread " << i << " start " << ptd[i].istart << " stop " << ptd[i].istop
           << '\n';
    }
  }

  int number_unselected = 0;
  for (int i = 0; i < pool_size; i++) {
    if (!pool[i].selected()) {
      number_unselected++;
    }
  }

  int method2 = 0;
  for (int i = 0; i < nthreads; i++) {
    for (int j = ptd[i].istart; j < ptd[i].istop; j++) {
      if (!pool[j].selected()) {
        method2++;
      }
    }
  }

  if (number_unselected == method2) {
    return 1;
  }

  cerr << "Mismatch on selected items " << number_unselected << " vs " << method2 << endl;

  for (int i = 0; i < nthreads; i++) {
    cerr << "Thread " << i << " start " << ptd[i].istart << " stop " << ptd[i].istop
         << '\n';
  }

  cerr << items_selected << " items selected, should be " << (pool_size - items_selected)
       << " unselected\n";

  return 1;
}

int
leader(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  cerr << my_mpi_id << "'s Pool_Size=" << pool_size << endl;

  assert(pool_size > 1);

  MPI_Allreduce(&pool_size, &my_mpi_pool_size, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);  // sum
  if (my_mpi_id == 0) {
    cerr << "Pool_Size_All=" << my_mpi_pool_size << endl;
  }

  assert(0 == items_selected);
  assert(0 == my_mpi_clusters_found);

  int first_unselected = 0;

  int icentre;

  Cluster cluster;
  cluster.resize(pool_size);

  if (!initialise_threads(nthreads)) {
    return 0;
  }

  int* cs_sum = NULL;
  int* disp = NULL;
  MY_MPI_STRUCT* my_mpi_cluster = NULL;
  MY_MPI_STRUCT* my_mpi_cluster_all = NULL;
  int my_mpi_struct_size_all = 0;  // total size of all process

  if (0 == my_mpi_id) {
    cs_sum = new int[my_mpi_size];
    disp = new int[my_mpi_size];
  }

#define IMPOSSIBLE_OFFSET 1

  int begin_id = 0;  // find centre from this process id
  do {
    off_t offset = IMPOSSIBLE_OFFSET;
    if (my_mpi_id == begin_id) {
      if (choose_next_centre(first_unselected, icentre)) {
        pool[icentre].offset(offset);
        pool[icentre].set_selected(1);
        items_selected += 1;  // for leader
      }
    }
    MPI_Bcast(&offset, 1, MPI_UNSIGNED_LONG_LONG, begin_id, MPI_COMM_WORLD);  //

    if (IMPOSSIBLE_OFFSET ==
        offset)  // didn't find a centre from begin_id, check next process
    {
      begin_id++;
      if (begin_id < my_mpi_size) {
        continue;
      } else {
        break;
      }
    }

    IW_TDT tdt;
    if (!get_tdt(tdt, input, offset)) {
      break;
    }
    int fatal;
    if (!my_mpi_centre.construct_from_tdt(tdt, fatal)) {
      break;
    }
    // my_mpi_centre.set_offset(offset);
    my_mpi_centre.set_distance(static_cast<similarity_type_t>(0.0));

    similarity_type_t my_threshold;
    if (my_mpi_centre.threshold(my_threshold)) {  // has come from the file
      ;
    } else {
      my_threshold = threshold;
    }

    if (verbose > 1) {
      cerr << "Start cluster " << my_mpi_clusters_found << ". ndx " << offset
           << ", threshold = " << my_threshold << endl;
    }

    //   pool[icentre].set_selected(1);
    //   pool[icentre].set_distance(static_cast<similarity_type_t>(0.0));

    for (int i = 0; i < nthreads; i++) {
      //    ptd[i].leader = icentre;
      ptd[i].threshold = my_threshold;
      sem_post(&start_forming_next_cluster[i]);
    }

    //  While the threads are forming the next cluster, process the cluster
    //  we just built

    if (my_mpi_items_selected > 0) {
      if (my_mpi_id == 0) {
        (void)process_cluster(my_mpi_cluster_all, my_mpi_struct_size_all, input, output);
        delete[] my_mpi_cluster_all;
      }
      my_mpi_clusters_found++;
      if (my_mpi_clusters_found >= max_clusters_to_find) {
        break;
      }
    }

    for (int i = 0; i < nthreads; i++) {
      sem_wait(&next_cluster_found);
    }

    cluster.resize_keep_storage(0);

    //  cluster.add(&(pool[icentre]));
    for (int i = 0; i < nthreads; i++) {
      cluster += ptd[i].sel;
    }

    int cs = cluster.number_elements();
    items_selected += cs;

    my_mpi_cluster = new MY_MPI_STRUCT[cs];
    for (int i = 0; i < cs; i++) {
      GFP_L& fp = *(cluster[i]);
      my_mpi_cluster[i].dis = fp.distance();
      fp.offset(my_mpi_cluster[i].offset);
    }

    MPI_Gather(&cs, 1, MPI_INT, cs_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_mpi_id == 0)  // prepare the buffer to accept the data
    {
      my_mpi_struct_size_all = 1;  // for centre
      int disp_sum = 1;            // for centre;
      for (int i = 0; i < my_mpi_size; i++) {
        my_mpi_struct_size_all += cs_sum[i];
        disp[i] = disp_sum;
        disp_sum += cs_sum[i];
      }

      my_mpi_cluster_all = new MY_MPI_STRUCT[my_mpi_struct_size_all];
      my_mpi_cluster_all[0].offset = offset;
      my_mpi_cluster_all[0].dis = 0;
    }
    MPI_Gatherv(my_mpi_cluster, cs, my_mpi_type, my_mpi_cluster_all, cs_sum, disp,
                my_mpi_type, 0, MPI_COMM_WORLD);

    delete[] my_mpi_cluster;

    MPI_Allreduce(&items_selected, &my_mpi_items_selected, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);  // sum the value
    // cerr<<"5"<<endl;

    if (my_mpi_clusters_found > 0 && 0 == my_mpi_clusters_found % reallocate_every) {
      reallocate_unselected_items(items_selected);
    }

  } while (my_mpi_items_selected < my_mpi_pool_size);

  if (my_mpi_struct_size_all > 0) {
    if (0 == my_mpi_id) {
      (void)process_cluster(my_mpi_cluster_all, my_mpi_struct_size_all, input, output);
      delete[] my_mpi_cluster_all;
    }
    my_mpi_clusters_found++;
  }

  if (0 == my_mpi_id) {
    delete[] cs_sum;
    delete[] disp;
  }

  if (0 == my_mpi_clusters_found) {
    return 0;
  }

  return 1;
}

static int
build_pool(iwstring_data_source& input)
{
  off_t offset = input.tellg();

  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    tdts_read++;
    if ((tdts_read - 1) / pool_size != my_mpi_id)  // distribute to different id
    {
      offset = input.tellg();
      continue;
    }

    int fatal;
    if (!pool[items_in_pool].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      offset = input.tellg();
      continue;
    }

    pool[items_in_pool].set_offset(offset);

    items_in_pool++;

    if (items_in_pool == pool_size) {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }

    offset = input.tellg();
  }

  pool_size = items_in_pool;

  if (verbose) {
    cerr << "Read " << tdts_read << " TDT's, pool contains " << pool_size
         << " fingerprints\n";
  }

  return 1;
}

static int
build_pool(const const_IWSubstring& fname, iwstring_data_source& input)
{
  IWString tmp(fname);

  if (!input.open(tmp))  // method is non-const on its argument!
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size) {
    IWString tmp;
    tmp << '^' << identifier_tag;

    std::unique_ptr<re2::RE2> pcn;
    iwre2::RE2Reset(pcn, tmp);
    pool_size = input.grep(*pcn);

    if (0 == pool_size) {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }
    pool_size = pool_size / my_mpi_size + 1;

    pool = new GFP_L[pool_size];
    if (NULL == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool(input);
}

/*
  If we have a previously selected file, we keep track of the number
  of members of the pool that get selected by the previously selected file
*/

static int molecules_selected_by_previously_selected_file = 0;

static int my_mpi_molecules_selected_by_previously_selected_file = 0;

static int
do_previously_selected_file(iwstring_data_source& input, similarity_type_t t)
{
  IW_TDT tdt;

  if (!tdt.next(input)) {
    cerr << "No first TDT in previously selected file\n";
    return 0;
  }

  IW_General_Fingerprint fp1;
  int fatal;
  if (!fp1.construct_from_tdt(tdt, fatal)) {  // make things easy
    return 0;
  }

  for (int i = 0; i < nthreads; i++) {
    pstd[i].threshold = t;
    pstd[i].previously_selected_fp = &fp1;
    sem_post(&(start_processing_previously_selected[i]));
  }

  IW_General_Fingerprint fp2;

  int ndx = 0;

  // We use two fingerprints. While one is off being compared, the other
  // is being built.

  while (tdt.next(input)) {
    IW_General_Fingerprint* fp;

    if (0 == ndx) {
      fp = &fp2;
      ndx = 1;
    } else {
      fp = &fp1;
      ndx = 0;
    }

    if (!fp->construct_from_tdt(tdt, fatal)) {
      return 0;
    }

    //  cerr << "Built fp '" << fp->id() << "'\n";

    for (int i = 0; i < nthreads; i++) {
      sem_wait(&finished_processing_previously_selected);
    }

    for (int i = 0; i < nthreads; i++) {
      molecules_selected_by_previously_selected_file += pstd[i].discarded_this_thread;

      pstd[i].previously_selected_fp = fp;
      sem_post(&(start_processing_previously_selected[i]));
    }
  }

  for (int i = 0; i < nthreads; i++) {
    sem_wait(&finished_processing_previously_selected);
  }

  MPI_Allreduce(&molecules_selected_by_previously_selected_file,
                &my_mpi_molecules_selected_by_previously_selected_file, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);  // sum the value
  return 1;
}

static int
do_previously_selected_file(const IWString& fname, similarity_type_t t)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected_file(input, t);
}

static int
set_thresholds_via_factor(score_t mvieth_factor)
{
  score_t global_max_score = -FLT_MAX;
  for (int i = 0; i < pool_size; i++) {
    score_t s = pool[i].score();
    if (s > global_max_score) {
      global_max_score = s;
    }
  }

  if (verbose) {
    cerr << "Max score in pool is " << global_max_score << endl;
  }

  for (int i = 0; i < pool_size; i++) {
    GFP_L& fp = pool[i];

    similarity_type_t t = (global_max_score - fp.score()) / mvieth_factor;
    fp.set_threshold(t);

    if (verbose > 1) {
      cerr << "i = " << i << " '" << fp.id() << " score " << fp.score()
           << " threshold set to " << t << endl;
    }
  }

  return 1;
}

static int
process_dash_e_option(Command_Line& cl, char e,
                      resizable_array_p<IWString>& items_to_echo)
{
  if (!cl.option_present(e)) {
    items_to_echo.resize(2);
    IWString* t = new IWString("$SMI<");
    items_to_echo.add(t);
    t = new IWString("PCN<");
    items_to_echo.add(t);

    return 1;
  }

  int all_found = 0;

  const_IWSubstring evalue;
  int i = 0;
  while (cl.value(e, evalue, i++)) {
    if ("ALL" == evalue) {
      all_found = 1;
      if (verbose) {
        cerr << "Will echo entire tdt on output\n";
      }
    } else {
      IWString* t = new IWString(evalue);
      items_to_echo.add(t);
      if (verbose) {
        cerr << "Will echo item '" << evalue << "'\n";
      }
    }
  }

  if (all_found && items_to_echo.number_elements()) {
    cerr << "Using '-" << e << " ALL' and other -" << e
         << " options doesn't make sense\n";
    return 0;
  }

  return 1;
}

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
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
  cerr << " -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -E ALL           echo all dataitems from the pool file\n";
//cerr << " -X <distance>    abandon distance computation if any component > distance\n";
//cerr << " -Y <factor>      thresold = (max_score - score) / factor\n";
//cerr << " -R <factor>      score = score + <factor> * distance to cluster\n";
  cerr << " -A <file>        file(s) of previously selected molecules - discard all within threshold\n";
  cerr << " -a <dist>        use <dist> as the threshold when comparing against the -A file\n";
  cerr << " -D ...           miscellaneous options, enter '-D help' for info\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -h <number>      number of threads to use\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
initialise_array_of_semaphores(sem_t* sem, int nsem)
{
  assert(nsem > 0);

  for (int i = 0; i < nsem; i++) {
    if (0 != sem_init(&(sem[i]), 0, 0)) {
      perror("sem_init error");
      return 0;
    }
  }

  return nsem;
}

static int
destroy_array_of_semaphores(sem_t* sem, int nsem)
{
  int rc = 1;

  for (int i = 0; i < nsem; i++) {
    if (0 != sem_destroy(&sem[i])) {
      perror("sem_destroy");
      rc = 0;
    }
  }

  return rc;
}

static int
leader(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vs:I:E:t:F:P:W:X:rH:S:C:O:Y:V:R:m:M:A:a:h:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool_size = pool_size / my_mpi_size + 1;

    pool = new GFP_L[pool_size];
    if (NULL == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose) {
      cerr << "system sized to " << pool_size << endl;
    }
  }

  if (cl.option_present('V')) {
    if (!tversky.parse_command_line(cl, 'V', verbose)) {
      cerr << "Cannot parse Tversky specifications\n";
      usage(38);
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_clusters_to_find) || max_clusters_to_find < 1) {
      cerr << "The -C option (max clusters to find) must be followed by a positive "
              "integer\n";
      usage(41);
    }

    if (verbose) {
      cerr << "Will find a max of " << max_clusters_to_find << " clusters\n";
    }
  }

  if (cl.option_present('r')) {
    sort_by_distance_from_centre = 1;
    if (verbose) {
      cerr << "Clusters sorted by distance from centre\n";
    }
  }

  // We need to be careful with the -i and -I options. Remember
  // that the pool is built first

  if (cl.option_present('I')) {
    (void)cl.value('I', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose) {
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
    }
  }

  if (cl.option_present('H')) {
    cl.value('H', threshold_from_file_tag);
    if (verbose) {
      cerr << "Each threshold from the '" << threshold_from_file_tag
           << "' dataitem in the input\n";
    }

    if (!threshold_from_file_tag.ends_with('<')) {
      threshold_from_file_tag << '<';
    }
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if (s.starts_with("col=")) {
      s.remove_leading_chars(4);
      if (!s.numeric_value(score_column) || score_column < 1) {
        cerr << "Invalid column for score '" << s << "'\n";
        usage(14);
      }

      if (verbose) {
        cerr << "Score for each item in column " << score_column << endl;
      }

      score_column--;
    } else {
      score_tag = s;
      if (verbose) {
        cerr << "Score tag is " << score_tag << "'\n";
      }
    }
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', cluster_distance_scale_factor) ||
        cluster_distance_scale_factor <= 0.0) {
      cerr << "The cluster distance scale factor option (-R) must be followed by a "
              "positive number\n";
      usage(19);
    }

    if (verbose) {
      cerr << "Scores adjusted by " << cluster_distance_scale_factor
           << " times distance to nearest cluster centre\n";
    }
  }

  if (cl.option_present('F') || cl.option_present('P') || cl.option_present('W')) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise general fingerprint options\n";
      usage(17);
    }
  }

  if (cl.option_present('X')) {
    if (!cl.value('X', abandon_distance_cutoff) || abandon_distance_cutoff < 0.0 ||
        abandon_distance_cutoff > 1.0) {
      cerr << "The -X option must be followed by a valid distance (0.0, 1.0)\n";
      usage(13);
    }

    if (verbose) {
      cerr << "Distance compuations abandoned if any component > "
           << abandon_distance_cutoff << endl;
    }
  }

  if (!process_dash_e_option(cl, 'E', dataitems_to_echo)) {
    cerr << "Cannot process -E option\n";
    usage(15);
  }

  if (!cl.option_present('t') && !cl.option_present('H') && !cl.option_present('Y')) {
    cerr << "Threshold distance must be specified via -t, -H or -V options\n";
    usage(28);
  }

  if (cl.option_present('t')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++)) {
      if (t.starts_with("col=")) {
        t.remove_leading_chars(4);
        if (!t.numeric_value(threshold_column) || threshold_column < 1) {
          cerr << "Invalid column for threshold '" << t << "'\n";
          usage(14);
        }

        if (verbose) {
          cerr << "Threshold for each item in column " << threshold_column << endl;
        }

        threshold_column--;
      } else if (t.starts_with("tag=")) {
        threshold_from_file_tag = t;
        threshold_from_file_tag.remove_leading_chars(4);
        if (verbose) {
          cerr << "Threshold in tag " << threshold_from_file_tag << "'\n";
        }

        if (!threshold_from_file_tag.ends_with('<')) {
          threshold_from_file_tag.add('<');
        }
      } else if (!t.numeric_value(threshold) || threshold < 0.0 || threshold > 1.0) {
        cerr << "The -t option must be followed by a valid distance value\n";
        usage(12);

        if (verbose) {
          cerr << "Distance threshold set to " << threshold << endl;
        }
      }
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_cluster_size) || max_cluster_size < 2) {
      cerr << "The -m (max cluster size) option must be followed by a whole number > 1\n";
      usage(43);
    }

    if (verbose) {
      cerr << "Max cluster size " << max_cluster_size << endl;
    }
  }

  if (cl.option_present('M')) {
    const_IWSubstring m = cl.string_value('M');

    if (m.starts_with("col=")) {
      m.remove_leading_chars(4);
      if (!m.numeric_value(max_cluster_size_column) || max_cluster_size_column < 1) {
        cerr << "The column for the per molecule maximum cluster size must be a whole "
                "positive number\n";
        usage(11);
      }

      if (verbose) {
        cerr << "The maximum cluster size per molecule will be in column "
             << max_cluster_size_column << endl;
      }

      max_cluster_size_column--;
    } else {
      max_cluster_size_tag = m;

      if (verbose) {
        cerr << "Max cluster size in '" << max_cluster_size_tag << "' tag\n";
      }
    }
  }

  if (cl.number_elements() > 1) {
    cerr << "Extra arguments ignored\n";
  }

  iwstring_data_source pool_file;

  if (!build_pool(cl[0], pool_file) || 0 == pool_size) {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 21;
  }

  if (!cl.option_present('h')) {
    cerr << "Must specify number of threads via the -h option\n";
    usage(3);
  }

  if (cl.option_present('h')) {
    if (!cl.value('h', nthreads) || nthreads < 1) {
      cerr << "The number of threads to use (-h) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will use " << nthreads << " threads\n";
    }

    start_forming_next_cluster = new sem_t[nthreads];
    if (!initialise_array_of_semaphores(start_forming_next_cluster, nthreads)) {
      return 8;
    }

    if (0 != sem_init(&next_cluster_found, 0, 0)) {
      perror("sem init error");
      return 6;
    }
  }

  if (cl.option_present('A')) {
    if (!cl.option_present('t') && !cl.option_present('a')) {
      cerr << "Must have a threshold available with the -A option (use -t or -a)\n";
      usage(11);
    }

    similarity_type_t t;

    if (cl.option_present('a')) {
      if (!cl.value('a', t) || t < 0.0 || t >= 1.0) {
        cerr << "Invalid value for previously selected threshold (-a option)\n";
        usage(4);
      }

      if (verbose) {
        cerr << "Will use " << t
             << " as the threshold for the previously selected list\n";
      }
    } else {
      t = threshold;
    }

    start_processing_previously_selected = new sem_t[nthreads];
    if (!initialise_array_of_semaphores(start_processing_previously_selected, nthreads)) {
      return 5;
    }

    if (!initialise_ps_child_threads(nthreads)) {
      return 8;
    }

    IWString fname;
    int i = 0;
    while (cl.value('A', fname, i++)) {
      if (!do_previously_selected_file(fname, t)) {
        cerr << "Cannot process previously selected file (-A option)\n";
        return 8;
      }
    }

    if (verbose) {
      cerr << "Rejected " << my_mpi_molecules_selected_by_previously_selected_file
           << " molecules by previously selected file(s)\n";
    }

    //  Clean up all the threads involved

    for (int i = 0; i < nthreads; i++) {
      pstd[i].thread_number = -1;
      pstd[i].previously_selected_fp = nullptr;
      sem_post(&(start_processing_previously_selected[i]));
    }

    for (int i = 0; i < nthreads; i++) {
      sem_wait(&finished_processing_previously_selected);
    }

    destroy_array_of_semaphores(start_processing_previously_selected, nthreads);

    delete[] ps_child_thread;
    sem_destroy(&finished_processing_previously_selected);

    if (my_mpi_molecules_selected_by_previously_selected_file == my_mpi_pool_size) {
      cerr << "Yipes, the previously selected file knocked out the whole pool\n";
      return 1;
    }
  }

  /*
    Michael Vieth wants to be able to set the threshold as

      threshold = (MAX_SCORE - score) / factor

    where factor is a user settable number
  */

  if (cl.option_present('Y')) {
    if (!cl.option_present('S')) {
      cerr << "Scores must be present (-S) in order to use -V\n";
      usage(42);
    }

    if (cl.option_present('H')) {
      cerr << "The -H (threshold in file) and -Y options are mutually exclusive\n";
      usage(31);
    }

    score_t mvieth_factor;

    if (!cl.value('Y', mvieth_factor) || mvieth_factor <= 0.0) {
      cerr << "The Vieth factor (-Y) must be followed by a positive number\n";
      usage(12);
    }

    if (verbose) {
      cerr << "Vieth factor " << mvieth_factor << endl;
    }

    set_thresholds_via_factor(mvieth_factor);
  }
  IWString_and_File_Descriptor output(1);

  try {
    leader(pool_file, output);
  } catch (const char* err) {
    cerr << "Caught '" << err << "' terminated\n";
    return 81;
  }

  if (verbose && 0 == my_mpi_id) {
    cerr << "Clustered " << my_mpi_pool_size << " fingerprints into "
         << my_mpi_clusters_found << " clusters\n";
    int isum = 0;
    for (int i = 0; i < cluster_size.number_elements(); i++) {
      int j = cluster_size[i];
      if (0 == j) {
        continue;
      }

      cerr << j << " clusters were of size " << i << " members\n";

      isum += j * i;
    }

    cerr << "In clusters " << isum << endl;
  }

  // Tell all the threads to destroy themselves

  if (nthreads) {
    for (int i = 0; i < nthreads; i++) {
      ptd[i].thread_number = -1;

      sem_post(&(start_forming_next_cluster[i]));
    }

    for (int i = 0; i < nthreads; i++) {
      sem_wait(&next_cluster_found);
    }

    if (0 != sem_destroy(&next_cluster_found)) {
      cerr << "Could not destroy semaphores\n";
    }

    destroy_array_of_semaphores(start_forming_next_cluster, nthreads);
  }

  delete[] start_forming_next_cluster;

  return 0;
}

int
main(int argc, char** argv)
{
  MPI_Init(0, 0);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_id);
  MPI_Comm_size(MPI_COMM_WORLD, &my_mpi_size);

  MY_MPI_STRUCT example;

  assert(sizeof(off_t) == sizeof(MPI_UNSIGNED_LONG_LONG));

  MPI_Datatype type[2] = {MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE};
  int block[2] = {1, 1};
  MPI_Aint disp[2];

  MPI_Address(&example, disp);
  MPI_Address(&example.dis, disp + 1);
  disp[1] = disp[1] - disp[0];
  disp[0] = 0;
  MPI_Type_struct(2, block, disp, type, &my_mpi_type);
  MPI_Type_commit(&my_mpi_type);

  int rc = leader(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  MPI_Type_free(&my_mpi_type);
  MPI_Finalize();

  return rc;
}
