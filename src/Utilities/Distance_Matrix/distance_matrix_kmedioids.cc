/*
  k-Medioids implementation
*/

#include <iostream>
#include <limits>
#include <random>
#include <utility>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"
#define IWDISTANCE_MATRIX_IMPLEMENTATION
#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

using std::pair;

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int clusters_to_form = 0;

static extending_resizable_array<int> cluster_size;

static int assign_clusters_sequentially = 0;
static int assign_clusters_randomly = 1;
static int assign_clusters_density = 0;
static int assign_clusters_spread = 0;

static int nopt_outer = 1;
static int max_inner_opt = 100;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString csize_tag ("CSIZE<");
static IWString distance_tag ("DIST<");

static IWString dummy_smiles;

static int min_cluster_size = -1;

static time_t tstop = static_cast<time_t>(0);

static int inner_loop_break_if_no_change = std::numeric_limits<int>::max();

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "k-medioids clustering\n";
  cerr << " -n <clus>      number of clusters to form\n";
  cerr << " -S <fname>     corresponding smiles file, enter 'NONE' for none\n";
  cerr << " -I ...         how to start initial clustering, enter '-I help' for choices\n";
  cerr << " -O <n>         how many outer optimisation steps to perform\n";
  cerr << " -o <n>         max number inner optimisation steps to perform\n";
  cerr << " -t <dist>      items with no neighbours within <dist> are singletons\n";
  cerr << " -m <number>    minimum cluster size\n";
  cerr << " -T <seconds>   stop clustering after <seconds> have elapsed\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  Since the fundamental need here is to be able to know, for each
  object, is it a medioid, or with which medioid is it associated,
  we adopt the convention that values in the arrays >= 0 means
  the object is not a medioid, and the value will be index of it's
  associated medioid

  A value -f -1 means the point is a medioid
*/

#define IS_MEDIOID -1
#define IS_SINGLETON -5

template <typename T>
class Cluster_Data
{
  private:
    int _n;

    const int _clusters_to_form;

    int * _cluster;

    int * _best_cluster;

    int * _previous_clustering;

    T * _longest_distance_in_cluster;

    T * _per_class_score;

//  private functions

    int  _detect_changes_while_copying_cluster_to_prev () const;

    void _do_assign_clusters_sequentially (const IWDistanceMatrixBase<T> & dm);
    void _do_assign_clusters_randomly (const IWDistanceMatrixBase<T> & dm);
    void _do_assign_clusters_density (const IWDistanceMatrixBase<T> & dm);
    void _do_assign_clusters_spread (const IWDistanceMatrixBase<T> & dm);

    void _assign_to_closest_medioid (const IWDistanceMatrixBase<T> & dm, int ndx);
    void _assign_to_closest_medioid (const IWDistanceMatrixBase<T> & dm);

    void _compute_per_class_score (const IWDistanceMatrixBase<T> & dm);
    T    _sum_per_class_scores () const;

    T   _outer_loop_optimise (const IWDistanceMatrixBase<T> & dm, T);

    int _inner_loop_optimise (const IWDistanceMatrixBase<T> & dm, T & score);

    int _something_clustered_with (const IWDistanceMatrixBase<T> & dm, int) const;

  public:
    Cluster_Data(int, int);
    ~Cluster_Data();

    int  identify_singletons (const IWDistanceMatrixBase<T> & dm, T);

    int do_remove_singletons (IWDistanceMatrixBase<T> & dm,
                                       IW_STL_Hash_Map_String & smiles,
                                       IWString_and_File_Descriptor & output);

    int items_in_cluster (int) const;

    int kmedioids (const IWDistanceMatrixBase<T> & dm);

    int write_clusters (const IW_STL_Hash_Map_String & smiles, const IWDistanceMatrixBase<T> & dm, IWString_and_File_Descriptor & output) const;
};

template <typename T>
Cluster_Data<T>::Cluster_Data(int n, int c) : _n(n), _clusters_to_form(c)
{ 
  _cluster = new int[n];
  _best_cluster = new int[n];
  _previous_clustering = new int[n];
  _longest_distance_in_cluster = new T[n];
  _per_class_score = new T[n];

  return;
}

template <typename T>
Cluster_Data<T>::~Cluster_Data()
{
  delete [] _cluster;
  delete [] _best_cluster;
  delete [] _previous_clustering;
  delete [] _longest_distance_in_cluster;
  delete [] _per_class_score;

  return;
}

template <typename T>
void
Cluster_Data<T>::_do_assign_clusters_randomly (const IWDistanceMatrixBase<T> & dm)
{
  for (int i = 0; i < _n; i++)
  {
    if (IS_SINGLETON != _cluster[i])
      _cluster[i] = 1;
  }

  int medioids_found = 0;

  std::random_device rd;
  std::default_random_engine generator{rd()};
  std::uniform_int_distribution<int> u(0, _n-1);

  while (medioids_found < _clusters_to_form)
  {
    int i = u(generator);

    if (_cluster[i] < 0)
      continue;

    _cluster[i] = IS_MEDIOID;
    medioids_found++;
  }

  _assign_to_closest_medioid(dm);
}

template <typename T>
void
Cluster_Data<T>::_do_assign_clusters_sequentially (const IWDistanceMatrixBase<T> & dm)
{
  for (int i = 0; i < _n; i++)
  {
    if (IS_SINGLETON != _cluster[i])
      _cluster[i] = 1;
  }

  int medioids_found = 0;

  for (int i = 0; i < _n && medioids_found < _clusters_to_form; i++)
  {
    if (_cluster[i] < 0)
      continue;

    _cluster[i] = IS_MEDIOID;
    medioids_found++;
  }

  _assign_to_closest_medioid(dm);
}

template <typename T>
void
Cluster_Data<T>::_do_assign_clusters_spread (const IWDistanceMatrixBase<T> & dm)
{
  int p1 = 0;
  int p2 = 0;
  T longest_distance = dm.zvalue_i_less_than_j(0, 1);

  for (int i = 0; i < _n; i++)
  {
    for (int j = i + 1; j < _n; j++)
    {
      T d = dm.zvalue_i_less_than_j(i, j);

      if (d < longest_distance)
        continue;

      longest_distance = d;
      p1 = i;
      p2 = j;
    }
  }

  _cluster[p1] = IS_MEDIOID;
  _cluster[p2] = IS_MEDIOID;

  T * distance_from_medoid = new T[_n]; std::unique_ptr<T[]> free_d(distance_from_medoid);

  set_vector(distance_from_medoid, _n, std::numeric_limits<T>::max());

  for (int i = 0; i < _n; i++)
  {
    if (i == p1)
    {
      distance_from_medoid[i] = static_cast<T>(0.0);
      continue;
    }

    T d = dm.zvalue(i, p1);

    if (distance_from_medoid[i] < d)
      distance_from_medoid[i] = d;

    if (i == p2)
    {
      distance_from_medoid[i] = static_cast<T>(0.0);
      continue;
    }

    d = dm.zvalue(i, p2);

    if (distance_from_medoid[i] < d)
      distance_from_medoid[i] = d;
  }

  int medioids_found = 2;

  while (medioids_found < _clusters_to_form)
  {
    T longest_distance = static_cast<T>(0.0);
    int ichoose = -1;

    for (int i = 0; i < _n; i++)
    {
      if (IS_MEDIOID == _cluster[i])
        continue;

      if (_cluster[i] < 0)
        continue;

      if (distance_from_medoid[i] > longest_distance)
      {
        longest_distance = distance_from_medoid[i];
        ichoose = i;
      }
    }

    if (ichoose < 0)   // how could this happen?
      return;

    _cluster[ichoose] = IS_MEDIOID;
    medioids_found++;
  }

  _assign_to_closest_medioid(dm);
}


template <typename T>
int
Cluster_Data<T>::identify_singletons (const IWDistanceMatrixBase<T> & dm,
                                      T threshold)
{
  int rc = 0;

  for (int i = 0; i < _n; i++)
  {
    int is_singleton = 1;

    for (int j = 0; j < i; j++)
    {
      if (IS_SINGLETON == _cluster[j])
        continue;

      T d = dm.zvalue_i_less_than_j(j, i);

      if (d > threshold)
        continue;

      is_singleton = 0;
      break;
    }

    if (! is_singleton)
      continue;

    for (int j = i + 1; j < _n; j++)
    {
      T d = dm.zvalue_i_less_than_j(i, j);

      if (d > threshold)
        continue;

      is_singleton = 0;
      break;
    }

    if (is_singleton)
    {
      _cluster[i] = IS_SINGLETON;
      rc++;
    }
  }

  return rc;
}

int
write_smiles_and_id (const IW_STL_Hash_Map_String & smiles,
                     const IWString & id,
                     IWString_and_File_Descriptor & output)
{
  IW_STL_Hash_Map_String::const_iterator f = smiles.find(id);

  if (f != smiles.end())
    output << smiles_tag << (*f).second << ">\n";
  else if (dummy_smiles.length())
    output << dummy_smiles;
  else
  {
    cerr << "No smiles for '" << id << "'\n";
    return 0;
  }

  output << identifier_tag << id << ">\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
write_cluster_begin (const IW_STL_Hash_Map_String & smiles,
                     const IWString & id,
                     int csize,
                     IWString_and_File_Descriptor & output)
{
  if (! write_smiles_and_id(smiles, id, output))
    return 0;

  output << csize_tag << csize << ">\n";

  return 1;
}

static int
write_singleton_cluster (const IW_STL_Hash_Map_String &smiles,
                         const IWString & id,
                         IWString_and_File_Descriptor & output)
{
  write_cluster_begin (smiles, id, 1, output);

  output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

template <typename T>
int
Cluster_Data<T>::do_remove_singletons (IWDistanceMatrixBase<T> & dm,
                                       IW_STL_Hash_Map_String & smiles,
                                       IWString_and_File_Descriptor & output) 
{
  int * to_remove = new int[_n]; std::unique_ptr<int[]> free_to_remove(to_remove);

  int items_removed = 0;
  for (int i = 0; i < _n; i++)
  {
    if (IS_SINGLETON == _cluster[i])
    {
      write_singleton_cluster(smiles, dm.id(i), output);
      to_remove[i] = 1;
      items_removed++;
    }
    else
      to_remove[i] = 0;
  }

  if (0 == items_removed)
  {
    cerr << "Cluster_Data::do_remove_singletons:no singletons\n";
    return 0;
  }

  if (! dm.remove_items(to_remove))
    return 0;

  _n = _n - items_removed;

  set_vector(_cluster, _n, 1);   // not sure what to set it to...

  return items_removed;
}

/*
  We assign our point the negative of the c value for it's closest medioid
*/

template <typename T>
void
Cluster_Data<T>::_assign_to_closest_medioid (const IWDistanceMatrixBase<T> & dm,
                          int ndx)
{
  int closest = -1;
  T shortest_distance = std::numeric_limits<T>::max();

  for (int i = 0; i < ndx; i++)
  {
    if (IS_MEDIOID != _cluster[i])   // only consider medioids
      continue;

    T d = dm.zvalue_i_less_than_j(i, ndx);

    if (d < shortest_distance)
    {
      shortest_distance = d;
      closest = i;
    }
  }

  for (int i = ndx + 1; i < _n; i++)
  {
    if (IS_MEDIOID != _cluster[i])   // only consider medioids
      continue;

    T d = dm.zvalue_i_less_than_j(ndx, i);

    if (d < shortest_distance)
    {
      shortest_distance = d;
      closest = i;
    }
  }

  if (closest < 0)   // not sure how that could happen
    return;

  _cluster[ndx] = closest;

  _per_class_score[closest] += shortest_distance;

  if (shortest_distance > _longest_distance_in_cluster[closest])
    _longest_distance_in_cluster[closest] = shortest_distance;

  return;
}

template <typename T>
void
Cluster_Data<T>::_assign_to_closest_medioid (const IWDistanceMatrixBase<T> & dm)
{
  set_vector(_longest_distance_in_cluster, _n, static_cast<T>(0));

  for (int i = 0; i < _n; i++)
  {
    if (_cluster[i] < 0)   // is a medioid or singleton
      continue;

    _assign_to_closest_medioid(dm, i);
  }

  return;
}

/*
  For each item, we keep track of the number of nbrs within the threshold
*/

typedef pair<int, int> Nbrs_Within_Threshold;

class Nbrs_Comparator
{
  private:
  public:
    int operator () (const Nbrs_Within_Threshold *, const Nbrs_Within_Threshold *) const;
};

int
Nbrs_Comparator::operator() (const Nbrs_Within_Threshold * n1,
                             const Nbrs_Within_Threshold * n2) const
{
  if (n1->second < n2->second)
    return 1;

  if (n1->second > n2->second)
    return -1;

  return 0;
}

/*
  Determine which molecules have the largest number of neighbours within
  some threshold
*/

static float threshold = 0.22;   // arbitrary

template <typename T>
void
Cluster_Data<T>::_do_assign_clusters_density (const IWDistanceMatrixBase<T> & dm)
{
  Nbrs_Within_Threshold ** nbrs = new Nbrs_Within_Threshold *[_n]; std::unique_ptr<Nbrs_Within_Threshold *[]> free_nbrs(nbrs);

  for (int i = 0; i < _n; i++)
  {
    Nbrs_Within_Threshold * t = new Nbrs_Within_Threshold;

    t->first = i;
    t->second = 0;

    nbrs[i] = t;

    if (IS_SINGLETON == _cluster[i])
      continue;

    for (int j = 0; j < i; j++)
    {
      T d = dm.zvalue_i_less_than_j(j, i);
      if (d < threshold)
        t->second++;
    }

    for (int j = i + 1; j < _n; j++)
    {
      T d = dm.zvalue_i_less_than_j(i, j);
      if (d < threshold)
        t->second++;
    }
  }

  Nbrs_Comparator nbc;
  iwqsort(nbrs, _n, nbc);

  for (int i = 0; i < _n; i++)
  {
    const Nbrs_Within_Threshold * t = nbrs[i];

    cerr << "Item " << t->first << " has " << t->second << " items within range\n";
  }

  for (int i = 0; i < _n; i++)
  {
    if (IS_SINGLETON != _cluster[i])
      _cluster[i] = _n + 1;    // an otherwise impossible value
  }

// Now that we have the densest molecules at the front, do a sort of leader
// Assign leaders 2 and cluster members 1

  int jstart = 0;
  int clusters_formed = 0;

  for (int i = 0; i < _clusters_to_form; i++)
  {
    for (int j = jstart; j < _n; j++)
    {
      const Nbrs_Within_Threshold * nj = nbrs[j];

      int jndx = nj->first;

      if (_cluster[jndx] < _n)  // already assigned a cluster
        continue;

      _cluster[jndx] = IS_MEDIOID;

      clusters_formed++;

      for (int k = j + 1; k < _n; k++)
      {
        const Nbrs_Within_Threshold * nk = nbrs[k];

        int kndx = nk->first;

        if (_cluster[kndx] < _n)
          continue;

        T d = dm.zvalue(jndx, kndx);

        if (d > threshold)
          continue;

        _cluster[kndx] = jndx;
      }
    }

    jstart++;
  }

//If we did not get enough, just choose some at random

  std::random_device rd;
  std::default_random_engine generator{rd()};
  std::uniform_int_distribution<int> u(0, _n-1);

  while (clusters_formed < _clusters_to_form)
  {
    int i = u(generator);

    if (_cluster[i] < _n)
      continue;

    _cluster[i] = IS_MEDIOID;

    clusters_formed++;
  }

  _assign_to_closest_medioid (dm);

  return;
}

template <typename T>
int
Cluster_Data<T>::items_in_cluster (int c) const
{
  return count_occurrences_of_item_in_array(c, _n, _cluster);
}

template <typename T>
int
Cluster_Data<T>::_inner_loop_optimise (const IWDistanceMatrixBase<T> & dm,
                                       T & score)
{
#ifdef CHECK_MEDOIDS
  int medoids = 0;
  for (int i = 0; i < _n; i++)
  {
    if (IS_MEDIOID == _cluster[i])
      medoids++;
  }

  cerr << "Initially " << medoids << " medoids\n";
#endif

  std::random_device rd;
  std::default_random_engine generator{rd()};
  std::uniform_int_distribution<int> u(0, _n-1);

  int istart = u(generator);

  int old_medoid = -1;

  for (int i = 0; i < _n; i++, istart++)
  {
    if (istart >= _n)
      istart = 0;

    if (IS_MEDIOID != _cluster[istart])
      ;
    else if (items_in_cluster(istart) < min_cluster_size)
      ;
    else
    {
      old_medoid = istart;
      break;
    }
  }

  if (old_medoid < 0)   // can this happen?
    return 0;

  int new_medoid = _something_clustered_with(dm, old_medoid);

  if (new_medoid < 0)   // medioid is singleton
    return 0;

// Now we need to see if any points in other clusters want to be associated
// with this new medioid, and where do the ones associated with the old
// medioid go...

  T old_medoid_to_new_medoid = dm.zvalue(old_medoid, new_medoid);

  _longest_distance_in_cluster[new_medoid] = static_cast<T>(0);

  _per_class_score[new_medoid] = static_cast<T>(0.0);

  int items_moved = 0;

  for (int i = 0; i < _n; i++)
  {
    if (IS_MEDIOID != _cluster[i])   // process the medioids only
      continue;

    if (i == old_medoid)
      continue;

    T d_old_to_i = dm.zvalue(old_medoid, i);

//  Use triangle-like inequality to see if cluster[i] can possibly be changed

    if (old_medoid_to_new_medoid + _longest_distance_in_cluster[i] < d_old_to_i)
      continue;

//  the items in cluster I needs to be examined

    _longest_distance_in_cluster[i] = static_cast<T>(0.0);
    _longest_distance_in_cluster[new_medoid] = static_cast<T>(0.0);

    for (int j = 0; j < _n; j++)
    {
      if (i != _cluster[j])
        continue;

      T dcurrent = dm.zvalue(j, i);
      T dnew     = dm.zvalue(j, new_medoid);

      if (dnew < dcurrent)    // moving out of cluster I
      {
        _cluster[j] = new_medoid;
        _per_class_score[new_medoid] += dnew;
        _per_class_score[i] -= dcurrent;
        if (dnew > _longest_distance_in_cluster[new_medoid])
          _longest_distance_in_cluster[new_medoid] = dnew;
        items_moved++;
      }
      else if (dcurrent > _longest_distance_in_cluster[i])   // staying in I
        _longest_distance_in_cluster[i] = dcurrent;
    }
  }

  if (0 == items_moved)
    return 0;

// Now all the items moving out of the old cluster

  _cluster[new_medoid] = IS_MEDIOID;

  _cluster[old_medoid] = old_medoid;   // so it gets processed in the next loop below

  for (int i = 0; i < _n; i++)
  {
    if (old_medoid == _cluster[i])
      _assign_to_closest_medioid(dm, i);
  }

  score = _sum_per_class_scores();

  return 1;
}

/*
  We preferentially look for things that are well separated from the existing
  medioid.
*/

template<typename T>
int
Cluster_Data<T>::_something_clustered_with (const IWDistanceMatrixBase<T> & dm,
                                            int medioid) const
{
  int rc = -1;

  std::random_device rd;
  std::default_random_engine generator{rd()};
  std::uniform_int_distribution<int> u(0, _n-1);

  int istart = u(generator);

  for (int i = 0; i < _n; i++)
  {
    if (_cluster[istart] != medioid)
      ;
    else if (dm.zvalue(istart, medioid) < 0.05)   // save for later in case needed
      rc = istart;
    else
      return istart;

    istart++;
    if (istart >= _n)
      istart = 0;
  }

// Maybe a singleton, or maybe all the things we tried were close to
// the medioid

  return rc;
}

template <typename T>
void
Cluster_Data<T>::_compute_per_class_score (const IWDistanceMatrixBase<T> & dm)
{
  set_vector(_per_class_score, _n, static_cast<T>(0));

  for (int i = 0; i < _n; i++)
  {
    if (_cluster[i] < 0)    // medioid or singleton
      continue;

    int m = _cluster[i];   // the medioid for item i
    
    assert (m >= 0 && m < static_cast<int>(dm.size()));

    T d = dm.zvalue(i, m);

    _per_class_score[m] += d;
  }

  return;
}

template <typename T>
T
Cluster_Data<T>::_sum_per_class_scores () const
{
  T rc = static_cast<T>(0);

  for (int i = 0; i < _n; i++)
  {
    if (IS_MEDIOID == _cluster[i])
      rc += _per_class_score[i];
  }

  return rc;
}

template <typename T>
void
report_progress (int iter,
                 const T score,
                 const T best_score,
                 const T global_lowest_score,
                 std::ostream & os)
{
  if (verbose > 2)   // always write
    ;
  else if (score > global_lowest_score)  // not interesting
    return;
  else if (score > best_score)  // not interesting
    return;

  os << " inner loop iteration " << iter << " score " << score;
  if (score < best_score)
      os << " *";
  os << '\n';

  return;
}

template <typename T>
int
Cluster_Data<T>::_detect_changes_while_copying_cluster_to_prev () const
{
  int changes_found = 0;

  for (int i = 0; i < _n; i++)
  {
    if (_cluster[i] != _previous_clustering[i])
    {
      _previous_clustering[i] = _cluster[i];
      changes_found = 1;
    }
  }

  return changes_found;
}

template <typename T>
T
Cluster_Data<T>::_outer_loop_optimise (const IWDistanceMatrixBase<T> & dm,
                                       T global_lowest_score)
{
  set_vector(_cluster, _n, _n + 1);

  if (assign_clusters_randomly)
    _do_assign_clusters_randomly(dm);
  else if (assign_clusters_density)
    _do_assign_clusters_density(dm);
  else if (assign_clusters_spread)
    _do_assign_clusters_spread(dm);
  else if (assign_clusters_sequentially)
    _do_assign_clusters_sequentially(dm);
  else
    return 0;

  _compute_per_class_score(dm);

  T lowest_score = _sum_per_class_scores();

  if (verbose > 2)
    cerr << "initial inner loop score " << lowest_score << endl;

  copy_vector (_previous_clustering, _cluster, _n);

  int last_change = max_inner_opt + 1;
  
  if (verbose > 1)
    cerr << "start " << max_inner_opt << " inner loop optimistations\n";

  T score;
  for (int i = 0; i < max_inner_opt; i++)
  {
    if (! _inner_loop_optimise(dm, score))
      continue;

    if (verbose)
      report_progress (i, score, lowest_score, global_lowest_score, cerr);

    if (score < lowest_score)   // great
      ;
    else if (i - last_change > inner_loop_break_if_no_change)
      break;

    if (! _detect_changes_while_copying_cluster_to_prev())
      break;

    lowest_score = score;
    last_change = i;
  }

  return lowest_score;
}

template <typename T>
int
Cluster_Data<T>::kmedioids (const IWDistanceMatrixBase<T> & dm)
{
  T lowest_score = std::numeric_limits<T>::max();

  for (int i = 0; i < nopt_outer; i++)
  {
    T d = _outer_loop_optimise(dm, lowest_score);

    if (d < lowest_score)
    {
      lowest_score = d;
      copy_vector(_best_cluster, _cluster, _n);
    }

    if (tstop > 0 && time(NULL) > tstop)
    {
      cerr << "Computation terminated by time limit, performed " << (i+1) << " outer optimisations\n";
      if (0 == verbose)
        cerr << "Lowest score " << lowest_score << endl;
      break;
    }
  }

  if (verbose)
    cerr << "Lowest score " << lowest_score << endl;

  return 1;
}

template <typename T>
int
Cluster_Data<T>::write_clusters (const IW_STL_Hash_Map_String & smiles,
                                 const IWDistanceMatrixBase<T> & dm,
                                 IWString_and_File_Descriptor & output) const
{
  for (int i = 0; i < _n; i++)
  {
    if (IS_SINGLETON != _cluster[i])
      continue;

    if (! write_cluster_begin (smiles, dm.id(i), 1, output))
      return 0;

    cluster_size[1]++;

    output << "|\n";
  }

  for (int i = 0; i < _n; i++)
  {
    if (IS_MEDIOID != _cluster[i])
      continue;

    int c = count_occurrences_of_item_in_array(i, _n, _cluster);

    cluster_size[c + 1]++;   // don't forget the mediod itself

    write_cluster_begin(smiles, dm.id(i), c + 1, output);

    for (int j = 0; j < _n; j++)
    {
      if (i != _cluster[j])
        continue;

      write_smiles_and_id(smiles, dm.id(j), output);
      output << distance_tag << dm.zvalue(i, j) << ">\n";
    }

    output << "|\n";
  }

  return 1;
}

static int
read_smiles_record (const const_IWSubstring & buffer,
                    IW_STL_Hash_Map_String & smiles)
{
  int i = 0;
  IWString smi, id;

  if (! buffer.nextword(smi, i) || ! buffer.nextword(id, i))
  {
    cerr << "Not enough tokens on smiles input\n";
    return 0;
  }

  if (smiles.contains(id))
  {
    cerr << "Duplicate smiles identifier '" << id << "'\n";
    return 0;
  }

  smiles[id] = smi;

  return 1;
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! read_smiles_record (buffer, smiles))
    {
      cerr << "Bad smiles '" << buffer << "'\n";
      return 0;
    }
  }

  return smiles.size();
}

static int
read_smiles (const char * fname,
             IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open smiles file '" << fname << "'\n";
    return 0;
  }

  return read_smiles (input, smiles);
}


static void
display_dash_i_options (std::ostream & os)
{
  os << " -I DEN          density influenced leader distrbution\n";
  os << " -I RND          random distrbution\n";
  os << " -I SEQ          sequentially assigned distrbution\n";
  os << " -I SPR          spread-like selection\n";
  
  exit(1);
}

static int
distance_matrix_kmedioids (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:I:O:o:S:m:t:T:b:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('m'))
  {
    if (! cl.value('m', min_cluster_size) || min_cluster_size < 1)
    {
      cerr << "The minimum cluster size option (-m) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will not destroy clusters with less than " << min_cluster_size << " items\n";
  }

  if (cl.option_present('I'))
  {
    const const_IWSubstring i = cl.string_value('I');

    if ("SEQ" == i)
    {
      assign_clusters_sequentially = 1;
      assign_clusters_randomly = 0;
      assign_clusters_density = 0;
      assign_clusters_spread = 0;
      if (verbose)
        cerr << "Initial cluster memberships assigned sequentially\n";
    }
    else if ("RND" == i)
    {
      assign_clusters_sequentially = 0;
      assign_clusters_randomly = 1;
      assign_clusters_density = 0;
      assign_clusters_spread = 0;
      if (verbose)
        cerr << "Initial cluster memberships assigned randomly\n";
    }
    else if ("DEN" == i)
    {
      assign_clusters_sequentially = 0;
      assign_clusters_randomly = 0;
      assign_clusters_density = 1;
      assign_clusters_spread = 0;
      if (verbose)
        cerr << "Initial cluster memberships assigned by density\n";
    }
    else if ("SPR" == i)
    {
      assign_clusters_sequentially = 0;
      assign_clusters_randomly = 0;
      assign_clusters_density = 0;
      assign_clusters_spread = 1;
      if (verbose)
        cerr << "Initial cluster memberships assigned by density\n";
    }
    else if ("help" == i)
    {
      display_dash_i_options(cerr);
    }
    else
    {
      cerr << "Unrecognised -I qualifier '" << i << "'\n";
      display_dash_i_options(cerr);
    }
  }

  if (! cl.option_present('n'))
  {
    cerr << "Must specify number of clusters to for via the -n option\n";
    usage(3);
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', clusters_to_form) || clusters_to_form < 2)
    {
      cerr << "The number of clusters (-n) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will form " << clusters_to_form << " clusters\n";
  }

  if (! cl.option_present('S'))
  {
    cerr << "Must specify smiles file with the -S option\n";
    usage(3);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWDistanceMatrixBase<float> dm;

  if (! dm.do_read(cl[0]))
  {
    cerr << "Cannot read distance matrix from '" << cl[0] << "'\n";
    return 4;
  }

  if (verbose)
    cerr << "Read distance matrix with " << dm.size() << " items from '" << cl[0] << "'\n";

  if (clusters_to_form >= dm.number_molecules())
  {
    cerr << "Distance matrix has " << dm.size() << " items, but request " << clusters_to_form << " clusters, impossible\n";
    return 4;
  }

  if (cl.option_present('m'))
  {
    if (min_cluster_size * clusters_to_form > dm.number_molecules())
    {
      cerr << "Requested " << clusters_to_form << " clusters, each at least " << min_cluster_size << " but only " << dm.size() << " items\n";
      return 4;
    }
  }

  if (cl.option_present('O'))
  {
    if (! cl.value('O', nopt_outer) || nopt_outer < 1)
    {
      cerr << "The number of outer optimisation steps (-O) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will perform " << nopt_outer << " outer optimisations\n";
  }
  else
    nopt_outer = dm.size() / 100;

  if (cl.option_present('o'))
  {
    if (! cl.value('o', max_inner_opt) || max_inner_opt < 1)
    {
      cerr << "The max number of inner optimisation steps (-o) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will perform a max of " << max_inner_opt << " inner optimisations\n";
  }
  else
    max_inner_opt = dm.size() / 100;

  if (cl.option_present('b'))
  {
    if (! cl.value('b', inner_loop_break_if_no_change) || inner_loop_break_if_no_change < 1)
    {
      cerr << "The break if no inner loop change option (-b) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will break inner loop optimisation if no improvement within " << inner_loop_break_if_no_change << " steps\n";
  }

  IW_STL_Hash_Map_String smiles;

  if (cl.option_present('S'))
  {
    const char * fname = cl.option_value('S');

    if (0 == strcmp(fname, "NONE"))
    {
      dummy_smiles << smiles_tag << 'C' << ">\n";
      cerr << "No smiles present\n";
    }
    else if (! read_smiles(fname, smiles))
    {
      cerr << "Cannot read smiles from '" << fname << "'\n";
      return 5;
    }
    else if (verbose)
      cerr << "Read " << smiles.size() << " id->smiles relationships from '" << fname << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  Cluster_Data<float> cd(dm.size(), clusters_to_form);

  if (cl.option_present('t'))
  {
    float t;
    if (! cl.value('t', t) || t <= 0.0 || t >= 1.0)
    {
      cerr << "Invalid singleton threshold (-t), must be valid distance\n";
      usage(3);
    }

    if (verbose)
      cerr << "Items with no nbrs within " << t << " put in singleton cluster\n";

    int singletons = cd.identify_singletons(dm, t);

    if (static_cast<unsigned int>(singletons) >= dm.size())
    {
      cerr << "At distance " << t << " all points are singletons\n";
      return 3;
    }

    if (verbose)
      cerr << "Identified " << singletons << " singletons\n";

    if (singletons)
      cd.do_remove_singletons(dm, smiles, output);

    if (static_cast<int>(dm.size()) <= clusters_to_form)
    {
      cerr << "At distance " << t << " there are " << singletons << " singletons, N = " << dm.size() << ", not enough points to produce " << clusters_to_form << " clusters\n";
      return 3;
    }

    if (min_cluster_size < 0)   // not set
      ;
    else if ((clusters_to_form * min_cluster_size) > static_cast<int>(dm.size()))
    {
      cerr << "Found " << singletons << " singletons, leaving " << dm.size() << " non singletons\n";
      cerr << "But request for " << clusters_to_form << " clusters with min size " << min_cluster_size << ", impossible\n";
      return 3;
    }
  }

  if (cl.option_present('T'))
  {
    int s;

    if (! cl.value('T', s) || s < 1)
    {
      cerr << "The time to cluster option (-T) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will perform clusterings for " << s << " seconds\n";

    tstop = time(NULL) + s;

    nopt_outer = std::numeric_limits<int>::max();
  }

  cd.kmedioids(dm);

  cd.write_clusters (smiles, dm, output);

  output.flush();

  if (verbose)
  {
    for (int i = 0; i < cluster_size.number_elements(); i++)
    {
      if (cluster_size[i])
        cerr << cluster_size[i] << " clusters had " << i << " members\n";
    }
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_kmedioids(argc, argv);

  return rc;
}
