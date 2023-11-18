/*
  An implementation of the leader algorithm for fingerprints
  This variant produces output that can be processed by nplotnn
*/

#include <pthread.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <limits>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/scalable_allocator.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwqsort/iwqsort.h"

#define LEADER_PARALLEL_IMPLEMENTATION

#include "leader_parallel.h"

#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;

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
//cerr << " -s <number>      specify max pool size\n";
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
  cerr << " -f               don't bother sorting within clusters at all - faster\n";
//cerr << " -X <distance>    abandon distance computation if any component > distance\n";
//cerr << " -R <factor>      score = score + <factor> * distance to cluster\n";
  cerr << " -A <file>        file(s) of previously selected molecules - discard all within threshold\n";
  cerr << " -a <dist>        use <dist> as the threshold when comparing against the -A file\n";
  cerr << " -L <file>        write fingerprints discarded by -A file(s)\n";
  cerr << " -D ...           miscellaneous options, enter '-D help' for info\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

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

static int process_clusters_on_formation = 1;

static similarity_type_t threshold = 0.0;

static int threshold_column = -1;

/*
  when dealing with clusters which are decided by the threshold, we can
  optinally sort the cluster members by their distance from the leader
*/

static int sort_by_distance_from_centre = 0;

/*
  Shell sort is fast for sorting the selected queue, but it does not
  preserve the input order.
*/

static int shell_sort_cluster_number = 1;

static int squeeze_selected_every = 0;

static int parallel_build_pool = 1;

extending_resizable_array<int> cluster_size;

static Accumulator<float> distance_stats;

static std::ofstream stream_for_discarded_by_previously_selected;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

/*
  Mar 2004. Qi Chen wants to do sub-clustering of the clusters, in order to
  sample an SAR within each cluster. Note that the way this is implemented
  now is somewhat strange because rather than marking the discarded molecules
  in some way, they simply do not appear in the output! Change if this ever
  becomes a problem.

  There are two ways of doing the within-cluster clustering.
  A fixed threshold applied to all clusters
  A constant multiple of the threshold used to form the cluster
*/

static similarity_type_t sub_cluster_threshold = static_cast<similarity_type_t>(0.0);

static float sub_cluster_threshold_ratio = static_cast<float>(0.0);

static int grainsize = 50000;

/*
  A single variable that indicates whether or not any form of cluster clustering
  is active
*/

static int leader_on_clusters = 0;

GFP_PL::GFP_PL()
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
GFP_PL::construct_from_tdt(IW_TDT& tdt, int& fatal)
{
  if (!IW_General_Fingerprint::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  if (threshold_from_file_tag.length()) {
    similarity_type_t tmp;
    if (!tdt.dataitem_value(threshold_from_file_tag, tmp) || tmp < 0.0) {
      cerr << "GFP_PL::construct_from_tdt: invalid '" << threshold_from_file_tag
           << "' in tdt\n";
      return 0;
    }

    _threshold.set(tmp);
  } else if (threshold_column >= 0) {
    const_IWSubstring t;
    if (!_id.word(threshold_column, t)) {
      cerr << "GFP_PL::construct_from_tdt: no " << threshold_column << " column in '"
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
      cerr << "GFP_PL::construct_from_tdt: cannot extract '" << score_tag
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
      cerr << "GFP_PL::construct_from_tdt: missing or invalid '" << max_cluster_size_tag
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

  if (!tdt.dataitem_value(smiles_tag, _smiles)) {
    cerr << "No smiles in tdt '" << _id << "'\n";
    return 0;
  }

  return 1;
}

/*
  Our pool is an array of FP objects
*/

static GFP_PL** pool = nullptr;

static int items_selected = 0;

/*
  Where to start inserting things into the unselected array
*/

static int next_slot_selected_queue = 0;

static GFP_PL** selected_queue = nullptr;

/*
  When we squeeze out selected molecules, pool_size will change
*/

static int initial_pool_size = 0;

static int pool_size = 0;

static int first_unselected = 0;

/*
  A cluster is a set of pointers to such objects
*/

typedef resizable_array<GFP_PL*> Cluster;

static int
squeeze_out_unselected_items()
{
  int ndx = 0;
  for (int i = 0; i < pool_size; i++) {
    GFP_PL* pi = pool[i];

    if (pi->selected()) {
      selected_queue[next_slot_selected_queue] = pi;
      next_slot_selected_queue++;
    } else {
      pool[ndx] = pi;
      ndx++;
    }
  }

  pool_size = ndx;

  first_unselected = 0;

  if (verbose > 1) {
    cerr << "Pool squeezed to " << ndx << " active items\n";
  }

  return ndx;
}

/*
  Each time we form a cluster, we move the selected items out
*/

#ifdef NOT_USED_
static int
squeeze_out_unselected_items(int items_in_cluster)
{
  int ndx = 0;
  for (int i = 0; i < pool_size; i++) {
    GFP_PL* pi = pool[i];
    if (pi->selected()) {
      selected_queue[next_slot_selected_queue] = pi;
      next_slot_selected_queue++;
      items_in_cluster--;
      if (0 == items_in_cluster)  // no need for any more checking of selected attribute
      {
        for (int j = i + 1; j < pool_size; j++) {
          pool[ndx] = pool[j];
          ndx++;
        }

        return 1;
      }
    } else {
      pool[ndx] = pi;
      ndx++;
    }
  }

  pool_size = ndx;

  first_unselected = 0;

  if (verbose) {
    cerr << "Pool squeezed to " << ndx << " active items, time " << time(NULL) << endl;
  }

  return ndx;
}
#endif

static int
determine_cluster_size(const GFP_PL* const* selected_queue, int initial_pool_size,
                       int istart)
{
  const int c = selected_queue[istart]->selected();

  // cerr << "Selected item " << istart << " in cluster " << c << endl;

  for (int i = istart + 1; i <= initial_pool_size; i++)  // note <=
  {
    //  cerr << " selected item " << i << " in cluster " << selected_queue[i]->selected()
    //  << endl;
    if (selected_queue[i]->selected() != c) {
      return i - istart;
    }
  }

  return 0;  // should not come here
}

static int
write_smiles_and_id(const GFP_PL& f, IWString_and_File_Descriptor& output)
{
  output << smiles_tag << f.smiles() << ">\n";
  output << identifier_tag << f.id() << ">\n";

  return 1;
}

static int
do_output(IWString_and_File_Descriptor& output)
{
  time_t tstart = time(NULL);

  if (verbose > 2) {
    cerr << "Begin output at " << tstart << endl;
  }

  int cluster = 0;

  int ndx = 0;

#ifdef ECHO_BEFORE_WRITING
  for (int i = 0; i < initial_pool_size; i++) {
    cerr << " selected item " << i << " in cluster " << selected_queue[i]->selected()
         << endl;
  }
#endif

  while (ndx < initial_pool_size) {
    //  if (verbose > 2)
    //    cerr << "Cluster " << cluster << " starts at " << ndx << endl;

    const GFP_PL* leader = selected_queue[ndx];
    write_smiles_and_id(*leader, output);

    output << "CLUSTER<" << cluster << ">\n";

    const int cs = determine_cluster_size(selected_queue, initial_pool_size, ndx);

    //  if (verbose > 2)
    //    cerr << "Contains " << cs << " items\n";

    output << "CSIZE<" << cs << ">\n";

    for (int j = 0; j < (cs - 1); j++) {
      const GFP_PL* c = selected_queue[ndx + j + 1];

      write_smiles_and_id(*c, output);
      if (c->distance() >= 0.0) {
        output << distance_tag << c->distance() << ">\n";
      }

      output.write_if_buffer_holds_more_than(32768);
    }

    output << "|\n";

    cluster++;
    ndx += cs;
  }

  output.flush();

  if (verbose) {
    time_t tnow = time(NULL);
    cerr << "Output took " << (tnow - tstart) << " seconds\n";
  }

  return 1;
}

static void
do_shell_sort_cluster_number(GFP_PL** const A, int size)
{
  int i, j, increment, temp;
  increment = size / 2;
  GFP_PL* save_i;

  while (increment > 0) {
    for (i = increment; i < size; i++) {
      j = i;
      temp = A[i]->selected();
      save_i = A[i];
      while ((j >= increment) && (A[j - increment]->selected() > temp)) {
        A[j] = A[j - increment];
        j = j - increment;
      }
      A[j] = save_i;
    }

    if (increment == 2) {
      increment = 1;
    } else {
      increment = (int)(increment / 2.2);
    }
  }
}

#ifdef NOT_USED_
static void
shell_sort_ndx(GFP_PL** const A, int size)
{
  int i, j, increment;
  int temp;
  GFP_PL* save_i;
  increment = size / 2;

  while (increment > 0) {
    for (i = increment; i < size; i++) {
      j = i;
      temp = A[i]->ndx();
      save_i = A[i];
      while ((j >= increment) && (A[j - increment]->ndx() > temp)) {
        A[j] = A[j - increment];
        j = j - increment;
      }
      A[j] = save_i;
    }

    if (increment == 2) {
      increment = 1;
    } else {
      increment = (int)(increment / 2.2);
    }
  }
}
#endif

/*static int
do_sort_by_initial_position ()
{
  if (do_quicksort_for_selected_queue)  // already done by quick
    return 0;

  int istart = 0;
  for (int i = 0; i < initial_pool_size; i++)
  {
    const int cluster_number = selected_queue[i]->selected();
    int cluster_size = 1;

    for (int j = i + 1; j <= initial_pool_size; j++)   // note the <=, we allocated an
extra item
    {
      if (cluster_number != selected_queue[j]->selected())
      {
        cluster_size = j - i;
        break;
      }
    }

    if (cluster_size > 1)
      shell_sort_ndx(selected_queue + i, cluster_size);
  }

  return 1;
}*/

/*
  The clusters can be sorted three different ways.
    shell sort on cluster number: fast, may destroy input order
    quicksort  on cluster number: slower, but preserves input order
    quicksort  on cluster number + dist:
*/

class Cluster_Number_Comparator
{
 private:
 public:
  int
  operator()(const GFP_PL*, const GFP_PL*) const;
};

class Cluster_Number_Distance_Comparator
{
 private:
 public:
  int
  operator()(const GFP_PL*, const GFP_PL*) const;
};

int
Cluster_Number_Comparator::operator()(const GFP_PL* ppf1, const GFP_PL* ppf2) const
{
  int c1 = (ppf1)->selected();
  int c2 = (ppf2)->selected();

  if (c1 < c2) {
    return -1;
  }

  if (c1 > c2) {
    return 1;
  }

  return 0;
}

int
Cluster_Number_Distance_Comparator::operator()(const GFP_PL* ppf1,
                                               const GFP_PL* ppf2) const
{
  int c1 = (ppf1)->selected();
  int c2 = (ppf2)->selected();

  if (c1 < c2) {
    return -1;
  }

  if (c1 > c2) {
    return 1;
  }

  similarity_type_t d1 = ppf1->distance();
  similarity_type_t d2 = ppf2->distance();

  if (d1 < d2) {
    return -1;
  }

  if (d1 > d2) {
    return 1;
  }

  return 0;
}

template <typename T, typename O, typename C>
void
shell_sort_template(O** const A, int size, C& fcn)
{
  int i, j, increment;
  T temp;
  O* save_i;

  while (increment > 0) {
    for (i = increment; i < size; i++) {
      j = i;
      temp = fcn(A[i]);
      save_i = A[i];
      while ((j >= increment) && (fcn(A[j - increment]) > temp)) {
        A[j] = A[j - increment];
        j = j - increment;
      }
      A[j] = save_i;
    }

    if (increment == 2) {
      increment = 1;
    } else {
      increment = (int)(increment / 2.2);
    }
  }
}

static void
randomise_order(GFP_PL** q, int n)
{
  int nsteps = n / 50;

  if (0 == nsteps) {
    nsteps = 10;
  }

  std::random_device rd;
  std::default_random_engine rng;
  std::uniform_int_distribution<off_t> u(0, n - 1);
  for (int i = 0; i < nsteps; i++) {
    int i1 = u(rng);
    int i2 = u(rng);
    if (i1 == i2) {
      continue;
    }
    std::swap(q[i1], q[i2]);
  }

  return;
}

static int
sort_selected_queue_by_cluster_number()
{
  time_t tstart = time(NULL);

#ifdef ECHO_BEFORE_SORT
  cerr << "Start sorting by cluster number at " << tstart << endl;

  for (int i = 0; i < next_slot_selected_queue; i++) {
    cerr << "Before sort, i = " << i << " cluster " << selected_queue[i]->selected()
         << endl;
  }
#endif

  // The input is almost sorted, so quicksort is horrible!

  if (shell_sort_cluster_number) {
    do_shell_sort_cluster_number(selected_queue, next_slot_selected_queue);
  } else {
    randomise_order(selected_queue, next_slot_selected_queue);

    if (sort_by_distance_from_centre) {
      Cluster_Number_Distance_Comparator cndc;
      iwqsort(selected_queue, next_slot_selected_queue, cndc);
    } else {
      Cluster_Number_Comparator cnc;
      iwqsort(selected_queue, next_slot_selected_queue, cnc);
    }
  }

  time_t tend = time(NULL);

  if (verbose) {
    cerr << "Sorting took " << (tend - tstart) << " seconds\n";
  }

  return 1;
}

#ifdef NOT_USED_
static int
process_clusters(IWString_and_File_Descriptor& output)
{
  int ndx = 0;

  for (int i = 0; i < pool_size; i++) {
    GFP_PL* pi = pool[i];

    if (pi->selected()) {
    } else {
      pool[ndx] = pi;
      ndx++;
    }
  }

  return 1;
}
#endif

static int
write_cluster_data(const GFP_PL& leader, int clusters_found, int id_within_cluster,
                   similarity_type_t distance_to_centre,
                   IWString_and_File_Descriptor& output)
{
  write_smiles_and_id(leader, output);

  if (distance_to_centre >= 0.0) {
    output << distance_tag << distance_to_centre << ">\n";
  }

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

static int
do_leader_on_clusters(Cluster& cluster, similarity_type_t threshold)
{
  similarity_type_t my_threshold;
  if (sub_cluster_threshold > static_cast<similarity_type_t>(0.0)) {
    my_threshold = sub_cluster_threshold;
  } else if (sub_cluster_threshold_ratio > static_cast<float>(0.0)) {
    my_threshold = threshold * sub_cluster_threshold_ratio;
  } else {
    cerr << "Not sure how to sub cluster\n";
    return 0;
  }

  for (int i = 0; i < cluster.number_elements(); i++) {
    GFP_PL* ldr = cluster[i];

    for (int j = i + 1; j < cluster.number_elements(); j++) {
      GFP_PL* cj = cluster[j];

      similarity_type_t d = ldr->IW_General_Fingerprint::distance(*cj);
      if (d <= my_threshold) {
        cluster.remove_item(j);
        j--;
      }
    }
  }

  return cluster.number_elements();
}

int
distance_comparitor(GFP_PL* const* ppfp1, GFP_PL* const* ppfp2)
{
  const GFP_PL* pfp1 = *ppfp1;
  const GFP_PL* pfp2 = *ppfp2;

  if (pfp1->distance() < pfp2->distance()) {
    return -1;
  } else if (pfp1->distance() > pfp2->distance()) {
    return 1;
  } else {
    return 0;
  }
}

static int
process_cluster(Cluster& cluster, similarity_type_t my_threshold,
                IWString_and_File_Descriptor& output)
{
  if (leader_on_clusters) {
    do_leader_on_clusters(cluster, my_threshold);
  }

  int cs = cluster.number_elements();
  cluster_size[cs]++;  // the leader isn't in the cluster

  GFP_PL* centre = cluster[0];

  if (verbose) {
    cerr << "Cluster " << clusters_found << ' ' << cs << " items, centre '"
         << cluster[0]->id() << "', ";
    if (threshold_from_file_tag.length()) {
      similarity_type_t threshold = 0.0f;
      (void)centre->threshold(threshold);
      cerr << "threshold " << threshold << ", ";
    }
    cerr << (items_selected + cs) << " items selected\n";
  }

  write_smiles_and_id(*centre, output);

  output << "CLUSTER<" << clusters_found << ">\n";
  output << "CSIZE<" << cs << ">\n";

  for (int i = 1; i < cs && output.good();
       i++)  // start at 1, we've already done centre above
  {
    GFP_PL* fp = cluster[i];

    if (!write_cluster_data(*fp, clusters_found, i, fp->distance(), output)) {
      return 0;
    }
  }

  output << "|\n";

  return output.good();
}

static int
choose_next_centre(int& icentre)
{
  icentre = -1;
  if (0 == score_tag.length() && score_column < 0)  // just grab the first unselected item
  {
    for (int i = first_unselected; i < pool_size; i++) {
      if (pool[i]->selected()) {
        continue;
      }

      first_unselected = i;
      icentre = i;
      return 1;
    }
  } else if (cluster_distance_scale_factor > static_cast<similarity_type_t>(0.0)) {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = first_unselected; i < pool_size; i++) {
      if (pool[i]->selected()) {
        continue;
      }

      score_t s = pool[i]->score() + cluster_distance_scale_factor *
                                         pool[i]->shortest_distance_to_cluster_centre();

      if (icentre < 0 || s > max_score) {
        max_score = s;
        icentre = i;
      }
    }
  } else  // raw scores
  {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = first_unselected; i < pool_size; i++) {
      if (pool[i]->selected()) {
        continue;
      }

      if (icentre < 0) {
        first_unselected = i;
      }

      score_t s = pool[i]->score();

      if (icentre < 0 || s > max_score) {
        max_score = s;
        icentre = i;
      }
    }
  }

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

class Form_Cluster_Threshold_TBB
{
 private:
  GFP_PL* _leader;
  const similarity_type_t _my_threshold;
  const int _cluster_number;

 public:
  Cluster cluster;

  Form_Cluster_Threshold_TBB(GFP_PL* l, similarity_type_t t, int c)
      : _leader(l), _my_threshold(t), _cluster_number(c)
  {
  }

  Form_Cluster_Threshold_TBB(const Form_Cluster_Threshold_TBB& rhs, tbb::split)
      : _leader(rhs._leader),
        _my_threshold(rhs._my_threshold),
        _cluster_number(rhs._cluster_number)  //, _myleader(fct_tbb._myleader)
  {
  }

  void
  operator()(const tbb::blocked_range<GFP_PL**>& r)
  {
    cluster.resize(r.end() - r.begin());
    for (GFP_PL** pfp = r.begin(); pfp != r.end(); ++pfp) {
      GFP_PL* fp = *pfp;

      if (fp->selected()) {
        continue;
      }

      similarity_type_t d;
      if (!compute_the_distance(*_leader, *fp, d)) {
        continue;
      }

      if (d <= _my_threshold) {
        cluster.add(fp);
        fp->set_selected(_cluster_number);
        fp->set_distance(d);
      }
      //      else if (d < fp->shortest_distance_to_cluster_centre())
      //               fp->set_shortest_distance_to_cluster_centre (d);
    }
  }

  void
  join(const Form_Cluster_Threshold_TBB& fct_tbb)
  {
    cluster += fct_tbb.cluster;
  }
};

static int
identify_first_and_last(GFP_PL** pool, int pool_size, int& first_unselected,
                        int& last_unselected)
{
  first_unselected = -1;
  last_unselected = -1;

  for (int i = 0; i < pool_size; i++) {
    GFP_PL* fpi = pool[i];

    if (fpi->selected()) {
      continue;
    }

    if (first_unselected < 0) {
      first_unselected = i;
    }

    last_unselected = i;
  }

  return first_unselected >= 0;
}

static int
calculate_grainsize(int zrange)
{
  if (zrange <= 0) {  // should not happen
    return 1;
  }

  if (zrange < 2000) {  // too small to multi-thread
    return zrange;
  }

  if (zrange <= grainsize) {
    return zrange / 5;  // 5 is arbitrary
  }

  return grainsize;
}

static int
form_cluster_threshold(int icentre, Cluster& cluster, int cluster_number,
                       const similarity_type_t my_threshold)
{
  GFP_PL* leader = pool[icentre];

  if (verbose > 2) {
    cerr << "Leader is " << leader->id() << " cluster " << cluster_number << endl;
  }

  int first_unselected, last_unselected;
  if (!identify_first_and_last(pool, pool_size, first_unselected, last_unselected)) {
    return 0;
  }

  if (verbose > 2) {
    cerr << "Start parallel, start " << first_unselected << " to " << last_unselected
         << " range " << (last_unselected - first_unselected) << endl;
  }

  Form_Cluster_Threshold_TBB fct_tbb(leader, my_threshold, cluster_number);
  tbb::parallel_reduce(tbb::blocked_range<GFP_PL**>(
                           pool + first_unselected, pool + last_unselected,
                           calculate_grainsize(last_unselected - first_unselected)),
                       fct_tbb);  //,tbb::auto_partitioner());

  cluster += fct_tbb.cluster;

  if (verbose > 2) {
    cerr << "Cluster " << cluster_number << " contains " << cluster.number_elements()
         << " items\n";
  }

#ifdef CHECK_CLUSTER_NUMBERS
  for (int i = 0; i < cluster.number_elements(); i++) {
    if (cluster[i]->selected() != cluster_number) {
      cerr << "Mismatch on cluster number " << cluster[i]->selected() << " vs "
           << cluster_number << endl;
    }
  }
#endif

  return cluster.number_elements();
}

class Form_Cluster_Max_Cluster_Size_TBB
{
 private:
  GFP_PL* _leader;
  similarity_type_t _my_threshold;

 public:
  Cluster cluster;

  Form_Cluster_Max_Cluster_Size_TBB(GFP_PL* l, similarity_type_t t)
      : _leader(l), _my_threshold(t)
  {
  }

  Form_Cluster_Max_Cluster_Size_TBB(Form_Cluster_Max_Cluster_Size_TBB& fct_tbb,
                                    tbb::split)
      : _leader(fct_tbb._leader),
      _my_threshold(fct_tbb._my_threshold) {
  }

  void
  operator()(const tbb::blocked_range<GFP_PL**>& r)
  {
    cluster.resize(r.end() - r.begin());
    for (GFP_PL** pfp = r.begin(); pfp != r.end(); ++pfp) {
      GFP_PL* fp = *pfp;

      if (fp->selected()) {
        continue;
      }

      similarity_type_t d;
      if (!compute_the_distance(*_leader, *fp, d)) {
        continue;
      }

      if (_my_threshold > static_cast<similarity_type_t>(0.0) && d > _my_threshold) {
        continue;
      }

      cluster.add(fp);
      fp->set_distance(d);
      //    if (d < fp->shortest_distance_to_cluster_centre())
      //      fp->set_shortest_distance_to_cluster_centre (d);
    }
  }

  void
  join(const Form_Cluster_Max_Cluster_Size_TBB& fct_tbb)
  {
    cluster += fct_tbb.cluster;
  }
};

static int
form_cluster_max_cluster_size(int icentre, Cluster& cluster, int cluster_number,
                              similarity_type_t my_threshold,
                              int max_cluster_size_this_molecule)
{
  assert(max_cluster_size_this_molecule > 0);

  cluster.resize(pool_size);

  int first_unselected;
  int last_unselected;
  if (!identify_first_and_last(pool, pool_size, first_unselected, last_unselected)) {
    return 0;
  }

  Form_Cluster_Max_Cluster_Size_TBB fct_tbb(pool[icentre], my_threshold);
  tbb::parallel_reduce(tbb::blocked_range<GFP_PL**>(
                           pool + first_unselected, pool + last_unselected,
                           calculate_grainsize(last_unselected - first_unselected)),
                       fct_tbb);  //,tbb::auto_partitioner());

  cluster += fct_tbb.cluster;

  cluster.sort(&distance_comparitor);

  cluster.resize_keep_storage(max_cluster_size_this_molecule);

  int istop;
  if (cluster.number_elements() < max_cluster_size_this_molecule) {
    istop = cluster.number_elements();
  } else {
    istop = max_cluster_size_this_molecule;
  }

  for (int i = 0; i < istop; i++) {
    GFP_PL* p = cluster[i];

    p->set_selected(cluster_number);
  }

  return 1;
}

/*
static int
form_cluster_max_cluster_size (int icentre, Cluster & cluster,
                               similarity_type_t my_threshold,
                               int max_cluster_size_this_molecule)
{
  assert (max_cluster_size_this_molecule > 0);

  cluster.resize (pool_size);

  GFP_PL & fp = pool[icentre];

  for (int i = 0; i < pool_size; i++)
  {
    GFP_PL & p = pool[i];

    if (p.selected())
      continue;

    if (! can_be_compared (fp, p))
      continue;

    similarity_type_t d;
    if (! compute_the_distance (fp, p, d))
      continue;

    if (my_threshold > static_cast<similarity_type_t> (0.0) && d > my_threshold)
      continue;

    cluster.add (&(pool[i]));
    p.set_distance (d);
    if (d < p.shortest_distance_to_cluster_centre())
      p.set_shortest_distance_to_cluster_centre (d);
  }

  cluster.sort (&distance_comparitor);

  cluster.resize_keep_storage (max_cluster_size_this_molecule);

  int istop;
  if (cluster.number_elements() < max_cluster_size_this_molecule)
    istop = cluster.number_elements();
  else
    istop = max_cluster_size_this_molecule;

  for (int i = 0; i < istop; i++)
  {
    GFP_PL * p = cluster[i];

    p->set_selected (1);
  }

  return 1;
}
*/

/*
  The clustering will be limited either by the maximum number of items which
  can be in a cluster, or a threshold
*/

static int
form_cluster(int icentre, Cluster& cluster, int cluster_number,
             const similarity_type_t my_threshold, int max_cluster_size_this_molecule)
{
  cluster.resize_keep_storage(0);

  cluster.add(pool[icentre]);

  pool[icentre]->set_selected(cluster_number);
  pool[icentre]->set_distance(static_cast<similarity_type_t>(0.0));
  if (max_cluster_size_this_molecule) {
    return form_cluster_max_cluster_size(icentre, cluster, cluster_number, my_threshold,
                                         max_cluster_size_this_molecule);
  } else {
    return form_cluster_threshold(icentre, cluster, cluster_number, my_threshold);
  }
}

int
leader(IWString_and_File_Descriptor& output)
{
  assert(pool_size > 1);

  assert(0 == items_selected);
  assert(0 == clusters_found);

  int icentre;
  if (!choose_next_centre(icentre)) {
    cerr << "Yipes, cannot find initial leader\n";
    return 0;
  }

  Cluster cluster;
  if (!cluster.resize(pool_size)) {
    cerr << "Yipes, cannot allocate " << pool_size << " elements in pool\n";
    return 0;
  }

  int next_time_to_squeeze = std::numeric_limits<int>::max();
  if (squeeze_selected_every > 0) {
    next_time_to_squeeze = squeeze_selected_every;
  }

  while (items_selected < initial_pool_size) {
    GFP_PL* centre = pool[icentre];

    similarity_type_t my_threshold;
    if (centre->threshold(my_threshold)) {  // has come from the file
      ;
    } else {
      my_threshold = threshold;
    }

    int max_cluster_size_this_molecule;

    if (max_cluster_size_tag.length() &&
        pool[icentre]->max_cluster_size(max_cluster_size_this_molecule)) {
      ;
    } else if (max_cluster_size_column >= 0 &&
               pool[icentre]->max_cluster_size(max_cluster_size_this_molecule)) {
      ;
    } else {
      max_cluster_size_this_molecule = max_cluster_size;
    }

    if (verbose > 1) {
      cerr << "Start cluster " << clusters_found << ". ndx " << icentre
           << ", threshold = " << my_threshold;
      if (max_cluster_size_this_molecule > 0) {
        cerr << ", max size " << max_cluster_size_this_molecule;
      }
      cerr << endl;
    }

    (void)form_cluster(icentre, cluster, clusters_found + 1, my_threshold,
                       max_cluster_size_this_molecule);

    cluster_size[cluster.number_elements()]++;

    if (process_clusters_on_formation) {
      process_cluster(cluster, my_threshold, output);
    }

    clusters_found++;
    if (clusters_found >= max_clusters_to_find) {
      break;
    }

    items_selected += cluster.number_elements();

    if (items_selected > next_time_to_squeeze) {
      squeeze_out_unselected_items();
      next_time_to_squeeze = items_selected + squeeze_selected_every;
    }

    if (!choose_next_centre(icentre)) {
      break;
    }
  }

  squeeze_out_unselected_items();  // make sure everything on the selected queue

  sort_selected_queue_by_cluster_number();

#ifdef ECHO_SORTED_CLUSTER_NUMBERS
  for (int i = 0; i < initial_pool_size; i++) {
    cerr << "Item " << i << " in cluster " << selected_queue[i]->selected() << endl;
  }
#endif

  if (process_clusters_on_formation) {
    return 1;
  }

  return do_output(output);
}

/*static int
do_allocate_pool (const char * fname,
                  int s)
{
  assert (0 == pool_size);

  if (s < 2)
  {
    cerr << "Cannot gather offsets from '" << fname << "', maybe empty file?\n";
    return 0;
  }

  pool_size = s;
  initial_pool_size = s;

  pool = new GFP_PL *[pool_size];

// Putting an extra item in the selected queue makes it easier to
// identify clusters later

  selected_queue = new GFP_PL *[pool_size + 1];

  if (NULL == pool || NULL == selected_queue)
  {
    cerr << "Cannot allocate " << pool_size << " items for pool\n";
    return 0;
  }

  if (verbose)
    cerr << "Pool sized for " << pool_size << " items\n";

  selected_queue[pool_size] = new GFP_PL;
  selected_queue[pool_size]->set_selected(-1);

  if (squeeze_selected_every > 0)   // user specified something
    ;
  else if (pool_size < 2000)   // not worth doing
    ;
  else
  {
    squeeze_selected_every = pool_size / 100;

    if (0 == squeeze_selected_every)   // hard to imagine
      squeeze_selected_every = 500;
  }

  return 1;
}*/

/*
  If we have a previously selected file, we keep track of the number
  of members of the pool that get selected by the previously selected file
*/

static int molecules_selected_by_previously_selected_file = 0;

class PS_TBB
{
 private:
  similarity_type_t _threshold;
  IW_General_Fingerprint* _prev_sel;
  int _discarded_this_thread;

 public:
  PS_TBB(PS_TBB& ps_tbb, tbb::split)
      : _threshold(ps_tbb._threshold), _prev_sel(ps_tbb._prev_sel)
  {
    _discarded_this_thread = 0;
  }

  PS_TBB(similarity_type_t t, IW_General_Fingerprint* f) : _threshold(t), _prev_sel(f)
  {
    _discarded_this_thread = 0;
  }

  int
  discarded_this_thread() const
  {
    return _discarded_this_thread;
  }

  void
  operator()(const tbb::blocked_range<GFP_PL**>& r)
  {
    for (GFP_PL** pfp = r.begin(); pfp != r.end(); ++pfp) {
      GFP_PL* fp = *pfp;

      if (fp->selected()) {
        continue;
      }

      similarity_type_t d;
      if (!compute_the_distance(*_prev_sel, *fp, d)) {
        continue;
      }

      if (d > threshold) {
        continue;
      }

      fp->set_selected(1);
      fp->set_distance(d);

      _discarded_this_thread++;
    }
  }

  void
  join(const PS_TBB& ps_tbb)
  {
    _discarded_this_thread += ps_tbb._discarded_this_thread;
  }
};

static int
do_previously_selected_file(IW_General_Fingerprint& fp, similarity_type_t t)
{
  int first_unselected;
  int last_unselected;
  if (!identify_first_and_last(pool, pool_size, first_unselected, last_unselected)) {
    return 0;
  }

  PS_TBB ps_tbb(t, &fp);
  tbb::parallel_reduce(tbb::blocked_range<GFP_PL**>(
                           pool + first_unselected, pool + last_unselected,
                           calculate_grainsize(last_unselected - first_unselected)),
                       ps_tbb);  //,tbb::auto_partitioner());

  molecules_selected_by_previously_selected_file = ps_tbb.discarded_this_thread();
  return 1;
}

/*
static int
do_previously_selected_file (IW_General_Fingerprint & fp,
                             similarity_type_t t)
{
  for (int i = 0; i < pool_size; i++)
  {
    GFP_PL & p = pool[i];

    if (p.selected())     // probably hit during an earlier pass
      continue;

    if (! can_be_compared (fp, p))
      continue;

    similarity_type_t d;
    if (! compute_the_distance (p, fp, d))
      continue;

    if (d > t)
      continue;

    p.set_selected (1);

    molecules_selected_by_previously_selected_file++;

    if (verbose > 1)
      cerr << p.id() << " distance " << d << " from previously selected '" << fp.id() <<
"'\n";

    if (stream_for_discarded_by_previously_selected.rdbuf()->is_open())
    {
      stream_for_discarded_by_previously_selected << identifier_tag << p.id() << ">\n";
      stream_for_discarded_by_previously_selected << identifier_tag << fp.id() << ">\n";
      stream_for_discarded_by_previously_selected << distance_tag << d << ">\n";
      stream_for_discarded_by_previously_selected << "|\n";
    }
  }

  return 1;
}*/

static int
do_previously_selected_file(IW_TDT& tdt, similarity_type_t t)
{
  IW_General_Fingerprint fp;

  int fatal;
  if (!fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot construct fingerprint from TDT\n";
    cerr << tdt;

    if (fatal) {
      return 0;
    }

    return 1;
  }

  return do_previously_selected_file(fp, t);
}

static int
do_previously_selected_file(iwstring_data_source& input, similarity_type_t t)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (!do_previously_selected_file(tdt, t)) {
      return 0;
    }
  }

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
display_dash_D_options(std::ostream& output)
{
  output << " -D sbc=<dist>    cluster each cluster, constant threshold <dist>\n";
  output << " -D sbcr=<ratio>  cluster each cluster. Use ratio of cluster threshold\n";
  output << " -D sqze=<n>      squeeze inactives out of pool every <n> items selected\n";
  output << " -D grsz=<n>      suggested problem subdivision (grainsize dividor)\n";

  exit(0);

  return 0;
}

static int
leader(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vs:I:t:F:P:W:X:rfH:S:C:V:R:m:M:A:L:a:D:eh:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('D')) {
    int i = 0;
    const_IWSubstring d;
    while (cl.value('D', d, i++)) {
      const_IWSubstring directive;
      float dvalue;

      if (!d.split_into_directive_and_value(directive, '=', dvalue)) {
        cerr << "Invalid -D qualifier '" << d << "'\n";
        display_dash_D_options(cerr);
      }

      if ("sbc" == directive) {
        sub_cluster_threshold = dvalue;
        if (verbose) {
          cerr << "Fixed sub-cluster threshold " << sub_cluster_threshold << endl;
        }

        leader_on_clusters = 1;
      } else if ("sbcr" == directive) {
        sub_cluster_threshold_ratio = dvalue;
        if (verbose) {
          cerr << "Variable sub-cluster threshold ratio " << sub_cluster_threshold_ratio
               << endl;
        }

        leader_on_clusters = 1;
      } else if ("sqze" == directive) {
        squeeze_selected_every = static_cast<int>(dvalue + 0.1);

        if (verbose) {
          cerr << "Will squeeze out selected items every " << squeeze_selected_every
               << " items selected\n";
        }
      } else if ("grsz" == directive) {
        grainsize = static_cast<int>(dvalue + 0.1);

        if (verbose) {
          cerr << "Grainsize " << grainsize << " items\n";
        }
      } else if ("dio" == directive) {
        process_clusters_on_formation = 0;
        if (verbose) {
          cerr << "Output delayed until programme end\n";
        }
      } else if ("help" == directive) {
        display_dash_D_options(cerr);
      } else {
        cerr << "Unrecognised -D qualifier '" << d << "'\n";
        display_dash_D_options(cerr);
      }
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('e')) {
    parallel_build_pool = 0;

    if (verbose) {
      cerr << "Serial I/O used for reading input\n";
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new GFP_PL*[pool_size];
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

  if (cl.option_present('r') && cl.option_present('f')) {
    cerr << "The -r and -f options are mutually incompatible\n";
    usage(3);
  }

  if (cl.option_present('r')) {
    sort_by_distance_from_centre = 1;
    if (verbose) {
      cerr << "Clusters sorted by distance from centre\n";
    }
  } else if (cl.option_present('f')) {
    shell_sort_cluster_number = 1;
    if (verbose) {
      cerr << "No within cluster sorting, order may be randomised\n";
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

  // Always initialise fingerprints here, we are reading things in parallel

  if (!initialise_fingerprints(cl, verbose)) {
    cerr << "Cannot initialise general fingerprint options\n";
    usage(17);
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
      } else {
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

  if (!initialise_fingerprints(cl[0], verbose)) {
    cerr << "Cannot initialise fingerprints using '" << cl[0] << "'\n";
    return 32;
  }

  Pool_Formation_Info<GFP_PL> pfi(verbose);

  pfi.set_number_single_file_parallel_readers(5);

  cerr << "What is pool_size " << pool_size << endl;
  if (0 == pool_size) {
    ;
  } else if (pfi.allocate_pool("cl", pool_size)) {
    ;
  } else {
    cerr << "Cannot size fingerprint pool (-s option)\n";
    return 8;
  }

  if (!pfi.build(cl)) {
    cerr << "Cannot read fingerprints\n";
    return 4;
  }

  pool_size = pfi.pool_size();
  pool = pfi.pool();

  if (verbose) {
    cerr << "Read " << pool_size << " fingerprints\n";
  }

  if (pool_size < 1) {
    cerr << "No fingerprints\n";
    return 5;
  }

  initial_pool_size = pool_size;

  if (squeeze_selected_every > 0) {  // user specified something
    ;
  } else if (pool_size < 2000) {  // not worth doing
    ;
  } else {
    squeeze_selected_every = pool_size / 100;

    if (0 == squeeze_selected_every) {  // hard to imagine
      squeeze_selected_every = 500;
    }
  }

  selected_queue = new GFP_PL*[pool_size + 1];

  selected_queue[pool_size] = new GFP_PL;
  selected_queue[pool_size]->set_selected(-1);

  if (cl.option_present('A')) {
    if (!cl.option_present('t') && !cl.option_present('a')) {
      cerr << "Must have a threshold available with the -A option (use -t or -a)\n";
      usage(11);
    }

    if (cl.option_present('L')) {
      IWString fname = cl.string_value('L');
      stream_for_discarded_by_previously_selected.open(fname.null_terminated_chars(),
                                                       std::ios::out);
      if (!stream_for_discarded_by_previously_selected.good()) {
        cerr << "Cannot open stream for molecules discarded by previously selected '"
             << fname << "'\n";
        return 4;
      }

      if (verbose) {
        cerr << "Molecules too close to previously selected file(s) written to '" << fname
             << "'\n";
      }
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

    IWString fname;
    int i = 0;
    while (cl.value('A', fname, i++)) {
      if (!do_previously_selected_file(fname, t)) {
        cerr << "Cannot process previously selected file (-A option)\n";
        return 8;
      }
    }

    if (verbose) {
      cerr << "Rejected " << molecules_selected_by_previously_selected_file
           << " molecules by previously selected file(s)\n";
    }

    if (molecules_selected_by_previously_selected_file == pool_size) {
      cerr << "Yipes, the previously selected file knocked out the whole pool\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (!leader(output)) {  // some kind of error
    rc = 4;
  }

  output.flush();

  if (verbose) {
    cerr << "Clustered " << initial_pool_size << " fingerprints into " << clusters_found
         << " clusters\n";
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

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = leader(argc, argv);

  return rc;
}
