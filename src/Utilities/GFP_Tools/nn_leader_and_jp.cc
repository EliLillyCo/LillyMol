/*
  We have run gfp_nearneighbours and can quickly do several clusterings.
*/

#include <stdlib.h>
#include <memory>
#include <limits>
#include <algorithm>
#include <fstream>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "cmdline.h"
#include "accumulator.h"
#include "set_or_unset.h"
#include "iwstring_data_source.h"
#include "iw_auto_array.h"
#include "iw_tdt.h"
#include "iw_stl_hash_map.h"
#include "timsort.hpp"
#include "misc.h"

using std::cout;
using std::endl;

typedef float similarity_type_t;

/*
  We can have the threshold for each item read from the file
*/

static IWString score_tag;

static int scores_present = 0;

/*
  Or the score may be embedded in the id
*/

static int score_column = -1;

static int verbose = 0;

static float maximum_distance = 1.0f;

static int write_tdt_output = 1;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = std::numeric_limits<int>::max();
static int max_cluster_size = 0;

/*
  We can process any number of thresholds
*/

static int number_thresholds = 0;
static similarity_type_t * threshold = nullptr;

static int select_next_leader_randomly = 0;

/*
  Note order dependency of initialisation...
*/

static std::random_device rd;
static std::mt19937_64 rng(rd());

/*
  The identifier tag used in each TDT
*/

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");
static IWString neighbour_tag("NBR<");

static IWString distance_tag("DIST<");

static int number_threshold_tags = 0;
static IWString * threshold_from_file_tag = nullptr;

static IWString max_cluster_size_tag;
static int max_cluster_size_column = -1;

static IWString cluster_number_tag("CLUSTER<");
static IWString cluster_size_tag("CSIZE<");

typedef float score_t;

static int taylor_butina = 0;

static int suppress_self_neighbours = 0;

class JarvisPatrickParameters
{
  private:
    int _nbr_list_size;
    int _nbrs_in_common;

  public:
    JarvisPatrickParameters();

    int active () const { return _nbrs_in_common > 0;}

    int nbr_list_size () const { return _nbr_list_size;}
    int nbrs_in_common () const { return _nbrs_in_common;}

    int build (const const_IWSubstring &);
};

JarvisPatrickParameters::JarvisPatrickParameters()
{
  _nbr_list_size = 0;
  _nbrs_in_common = 0;

  return;
}

int
JarvisPatrickParameters::build (const const_IWSubstring & buffer)
{
  const_IWSubstring snl, snb;
  if (! buffer.split(snl, ',', snb) || 0 == snl.length() || 0 == snb.length())
  {
    cerr << "JarvisPatrickParameters::build:JP parameters must be of the form 'J,K' where J is the neighbour list size\n";
    cerr << "  and K is the number of neighbours in common\n";
    return 0;
  }

  if ('.' == snl)    // leave as default
    ;
  else if (! snl.numeric_value(_nbr_list_size) || _nbr_list_size < 1)
  {
    cerr << "JarvisPatrickParameters::build:invalid neighbour list size '" << buffer << "'\n";
    return 0;
  }

  if (! snb.numeric_value(_nbrs_in_common) || _nbrs_in_common < 1)
  {
    cerr << "JarvisPatrickParameters::build:invalid neighbours in common size '" << buffer << "'\n";
    return 0;
  }

  if (0 == _nbr_list_size)
    ;
  else if (_nbr_list_size < _nbrs_in_common)
  {
    cerr << "JarvisPatrickParameters::build:invalid combination '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

class Leader_Item
{
  private:
    int _nbrs;

    int * _nbr;
    float * _dist;

//  Each item can limit the number of items in its cluster

    Set_or_Unset<int> _max_cluster_size;

//  Each item may have any number of thresholds associated with it

    similarity_type_t * _threshold;

    score_t _score;

    int _ndx;

//  private functions


  public:
    Leader_Item();
    ~Leader_Item();

    int threshold (similarity_type_t & zresult, int which_one) const;

    int neighbour_count() const { return _nbrs;}

    void set_ndx (int s) {_ndx = s;}

    score_t score () const { return _score;}
    void set_score (const score_t s) { _score = s;}

    void update_max_distance (float & d) const;

    int max_cluster_size (int & m) const { return _max_cluster_size.value (m);}

    int build (iwstring_data_source & input, const IW_STL_Hash_Map_int & id_to_ndx, int n);

    int form_cluster_threshold (int * selected, const float threshold, const int cluster_number, const resizable_array_p<IWString> & smiles, const resizable_array_p<IWString> & id, std::ostream & output);
    int form_cluster_max_cluster_size(int * selected, const float threshold, const int max_cluster_size, const int cluster_number, const resizable_array_p<IWString> & smiles, const resizable_array_p<IWString> & id, std::ostream & output);

    int  mark_neighbours (int * n, const int maxn, const int x) const;
    bool neighbours_in_common (const int * n, const JarvisPatrickParameters &, const int x) const;

    int unselected_neighbours (const int * selected) const;
    int unselected_neighbours (const int * selected, const similarity_type_t within) const;
};

Leader_Item::Leader_Item ()
{
  _nbrs = 0;
  _nbr = nullptr;
  _dist = nullptr;
  _threshold = nullptr;
  _score = static_cast<score_t>(0.0);
  _ndx = -1;

  return;
}

Leader_Item::~Leader_Item()
{
  if (nullptr != _nbr)
    delete [] _nbr;

  if (nullptr != _dist)
    delete [] _dist;

  if (nullptr != _threshold)
    delete [] _threshold;

  return;
}

void
Leader_Item::update_max_distance (float & d) const
{
  if (0 == _nbrs)
    return;

  if (_dist[_nbrs - 1] > d)
    d = _dist[_nbrs - 1];

  return;
}

static int
data_from_tdt (const const_IWSubstring buffer,
               resizable_array_p<IWString> & s)
{
  int openangle = buffer.index('<');
  if (openangle < 0)    // huh?
    return 0;

  IWString * tmp = new IWString;

  buffer.from_to(openangle + 1, buffer.length() - 2, *tmp);

  s.add(tmp);

  return 1;
}

static int
read_smiles_and_ids_and_nbr_count (iwstring_data_source & input,
                                   resizable_array_p<IWString> & smiles,
                                   resizable_array_p<IWString> & id,
                                   resizable_array<int> & ncount)
                
{
  const_IWSubstring buffer;

  int got_smiles_this_molecule = 0;
  int got_id_this_molecule = 0;
  int nbrs_this_molecule = 0;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with(smiles_tag))
    {
      if (0 == got_smiles_this_molecule)
      {
        data_from_tdt(buffer, smiles);
        got_smiles_this_molecule = 1;
      }
    }
    else if (buffer.starts_with(identifier_tag))
    {
      if (0 == got_id_this_molecule)
      {
        data_from_tdt(buffer, id);
        got_id_this_molecule = 1;
      }
      else
        nbrs_this_molecule++;
    }
    else if ('|' == buffer)
      break;
  }

  if (! got_smiles_this_molecule && ! got_id_this_molecule)
    return 0;

  ncount.add(nbrs_this_molecule);

  return 1;
}

static int
read_smiles_and_ids_and_nbr_count (iwstring_data_source & input,
                                   const int max_pool_size,
                                   resizable_array_p<IWString> & smiles,
                                   resizable_array_p<IWString> & id,
                                   resizable_array<int> & ncount)
{
  while (read_smiles_and_ids_and_nbr_count(input, smiles, id, ncount))
  {
    if (smiles.size() != id.size() || id.size() != ncount.size())
    {
      cerr << "Mismatch on smiles " << smiles.size() << " ids " << id.size() << " and/or ncount " << ncount.size() << " at line " << input.lines_read() << endl;
      return 0;
    }

    if (smiles.number_elements() >= max_pool_size)
      break;
  }

  if (0 == smiles.size())
  {
    cerr << "read_smiles_and_ids_and_nbr_count:no data\n";
    return 0;
  }

  return smiles.number_elements();
}

static int
read_smiles_and_ids_and_nbr_count (const char * fname,
                                   const int max_pool_size,
                                   resizable_array_p<IWString> & smiles,
                                   resizable_array_p<IWString> & id,
                                   resizable_array<int> & ncount)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_smiles_and_ids_and_nbr_count:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles_and_ids_and_nbr_count(input, max_pool_size, smiles, id, ncount);
}

static int
get_nbr_id (const const_IWSubstring & buffer,
            const IW_STL_Hash_Map_int & id_to_ndx,
            int & nbr)
{
  int openangle = buffer.index('<');
  if (openangle < 0)
    return 0;

  IWString id;
  buffer.from_to(openangle + 1, buffer.length() - 2, id);

  id.truncate_at_first(' ');

  const auto f = id_to_ndx.find(id);

  if (f == id_to_ndx.end())
  {
    cerr << "get_nbr_id:no number for '" << id << "'\n";
    return 0;
  }

  nbr = f->second;

//cerr << "From '" << buffer << "' id '" << id << "' nbr " << nbr << endl;

  return 1;
}

template <typename T>
int
convert_to_float (const const_IWSubstring & buffer,
                  T & d,
                  const T minv,
                  const T maxv)
{
  int openangle = buffer.index('<');

  const_IWSubstring s;
  buffer.from_to(openangle + 1, buffer.length() - 2, s);

  if (! s.numeric_value(d) || d < minv || d > maxv)
  {
    cerr << "convert_to_float:invalid distance '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
Leader_Item::build (iwstring_data_source & input,
                    const IW_STL_Hash_Map_int & id_to_ndx,
                    int n)
{
  _nbrs = n;
  _nbr = new int[_nbrs];
  _dist = new float[_nbrs];

  IWString id;

  int ndx_nbrs = 0;                 // index into _nbrs
  int ndx_dist = 0;               // index into _dist array

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
//  cerr << "Leader_Item::build:reading '" << buffer << "'\n";

    if ('|' == buffer)
      break;

    if (buffer.starts_with(identifier_tag))
    {
      if (0 == id.length())
        buffer.from_to(identifier_tag.length() + 1, buffer.length() - 2, id);
      else if (! get_nbr_id(buffer, id_to_ndx, _nbr[ndx_nbrs]))
      {
        cerr << "Leader_Item::build:invalid neighbour specification '" << buffer << "'\n";
        return 0;
      }
      else
        ndx_nbrs++;
    }
    else if (buffer.starts_with(distance_tag))
    {
      if (! convert_to_float(buffer, _dist[ndx_dist], 0.0f, maximum_distance))
      {
        cerr << "Leader_Item::build:invalid distance '" << buffer << "'\n";
        return 0;
      }

      ndx_dist++;
    }
    else if (score_tag.length() && buffer.starts_with(score_tag))
    {
      if (! convert_to_float(buffer, _score, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max()))
      {
        cerr << "Leader_Item::build:invalid score '" << buffer << "'\n";
        return 0;
      }
    }
  }

#ifdef DEBUG_BUILD
  cerr << "Item " << _ndx << " has " << _nbrs << " neighbours\n";
  for (int i = 0; i < _nbrs; ++i)
  {
    cerr << ' ' << _nbr[i] << endl;
  }
#endif

  if (ndx_nbrs != ndx_dist)
  {
    cerr << "Leader_Item::build:mismatch btw ids " << ndx_nbrs << " and dists " << ndx_dist << endl;
    return 0;
  }

  if (score_column >= 0)
  {
    IWString s;
    if (! id.word(score_column, s))
    {
      cerr << "Leader_Item::build:cannot extract column " << (score_column+1) << " from '" << id << "'\n";
      return 0;
    }

    if (! s.numeric_value(_score))
    {
      cerr << "Leader_Item::build:invalid score '" << id << "', column " << (score_column+1) << "\n";
      return 0;
    }
  }

  if (max_cluster_size_column >= 0)
  {
    IWString s;
    if (! id.word(max_cluster_size_column, s))
    {
      cerr << "Leader_Item::build:cannot extract column " << (score_column+1) << " from '" << id << "'\n";
      return 0;
    }

    int m;
    if (! s.numeric_value(m) || m < 1)
    {
      cerr << "Leader_Item::build:invalid max cluster size " << id << " column " << (max_cluster_size_column+1) << endl;
      return 0;
    }

    _max_cluster_size.set(m);
  }

  return 1;
}

static int
write_cluster_header (const int ndx,
                      const resizable_array_p<IWString> & smiles,
                      const resizable_array_p<IWString> & id,
                      const int cluster_number,
                      const int csize,
                      std::ostream & output)
{
  if (write_tdt_output)
  {
    output << smiles_tag << *(smiles[ndx]) << ">\n";
    output << identifier_tag << *(id[ndx]) << ">\n";
    output << cluster_number_tag << cluster_number << ">\n";
    output << cluster_size_tag << csize << ">\n";
  }
  else
  {
    output << *(smiles[ndx]) << ' ' << *(id[ndx]) << " CLUSTER " << cluster_number << " (" << csize << " members)\n";
  }


  return 1;
}

static void
write_smiles_id_dist (const IWString * smiles,
                      const IWString * id,
                      const float d,
                      std::ostream & output)
{
  if (write_tdt_output)
  {
    output << smiles_tag << *(smiles) << ">\n";
    output << identifier_tag << *(id) << ">\n";
    if (d >= 0.0f)
      output << distance_tag << d << ">\n";
  }
  else
  {
    output << *smiles << ' ' << *id;
    if (d >= 0.0f)
      output << ' ' << d;
    output << '\n';
  }

  return;
}

static int
finish_cluster (const int csize,
                std::ostream & output)
{
  if (write_tdt_output)
    output << "|\n";

  return csize;
}

int
Leader_Item::form_cluster_threshold (int * selected,
                                     const float threshold,
                                     const int cluster_number,
                                     const resizable_array_p<IWString> & smiles,
                                     const resizable_array_p<IWString> & id,
                                     std::ostream & output)
{
  int csize = 1;

  for (int i = 0; i < _nbrs; ++i)
  {
    if (_dist[i] > threshold)
      break;

    const int j = _nbr[i];

    if (selected[j])
      continue;

    if (suppress_self_neighbours && j == _ndx)
      continue;

    csize++;
  }

  write_cluster_header(_ndx, smiles, id, cluster_number, csize, output);

  if (1 == csize)
  {
    selected[_ndx] = 1;
    return finish_cluster(1, output);
  }

  if (suppress_self_neighbours)
    selected[_ndx] = 1;

  for (int i = 0; i < _nbrs; ++i)
  {
    if (_dist[i] > threshold)
      break;

    const int j = _nbr[i];

    if (selected[j])
      continue;

    write_smiles_id_dist(smiles[j], id[j], _dist[i], output);

    selected[j] = 1;
  }

  selected[_ndx] = 1;          

  return finish_cluster(csize, output);
}

int
Leader_Item::form_cluster_max_cluster_size(int * selected,
                                           const float threshold,
                                           const int max_cluster_size,
                                           const int cluster_number,
                                           const resizable_array_p<IWString> & smiles,
                                           const resizable_array_p<IWString> & id,
                                           std::ostream & output)
{
  int csize = 1;

  for (int i = 0; i < _nbrs; ++i)    // first work out the cluster size
  {
    if (_dist[i] > threshold)
      break;

    const auto j = _nbr[i];

    if (selected[j])
      continue;

    csize++;

    if (csize >= max_cluster_size)
      break;
  }

  selected[_ndx] = 1;

  write_cluster_header(_ndx, smiles, id, cluster_number, csize, output);

  if (1 == csize)
    return finish_cluster(1, output);

  csize = 0;

  for (int i = 0; i < _nbrs; ++i)
  {
    if (_dist[i] > threshold)
      break;

    int j = _nbr[i];

    if (selected[j])
      continue;

    write_smiles_id_dist(smiles[j], id[j], _dist[i], output);

    selected[j] = 1;

    csize++;

    if (csize >= max_cluster_size)
      break;
  }

  return finish_cluster(csize, output);
}

int
Leader_Item::threshold(similarity_type_t & zresult, int which_one) const
{
  assert (which_one >= 0);

  similarity_type_t t = _threshold[which_one];

  if (t >= 0.0)     // valid value, therefore we have a threshold
  {
    zresult = t;
    return 1;
  }

  return 0;
}

int
Leader_Item::mark_neighbours (int * n, const int maxn, const int x) const
{
  int istop = _nbrs;
  if (0 == maxn)    // use all neighbours
    ;
  else if (istop > maxn)
    istop = maxn;

  for (int i = 0; i < istop; ++i)
  {
    n[_nbr[i]] = x;
  }

  return _nbrs;
}

bool
Leader_Item::neighbours_in_common (const int * n,
                                   const JarvisPatrickParameters & jp,
                                   const int x) const
{
  int istop = _nbrs;
  if (0 == jp.nbr_list_size())
    ;
  else if (istop > jp.nbr_list_size())
    istop = jp.nbr_list_size();

  int nic = 0;

  for (int i = 0; i < istop; ++i)
  {
    if (x != n[_nbr[i]])
      continue;

    nic++;

    if (nic >= jp.nbrs_in_common())
      return true;
  }

  return false;
}

int
Leader_Item::unselected_neighbours(const int * sel) const
{
  int rc = 0;

  for (int i = 0; i < _nbrs; ++i)
  {
    if (! sel[_nbr[i]])
      rc++;
  }

  return rc;
}

int
Leader_Item::unselected_neighbours(const int * sel, const similarity_type_t within) const
{
  if (0 == _nbrs)
    return 0;

  const int last_nbr = _nbr[_nbrs-1];

  if (_dist[last_nbr] <= within)    // all nbrs within distance, no need to check individually
    return unselected_neighbours(sel);

  int rc = 0;

  for (int i = 0; i < _nbrs; ++i)
  {
    if (sel[_nbr[i]])
      continue;

    if (_dist[i] > within)
      break;

    rc++;
  }

  return rc;
}


/*
  Items c1 and c2 belong in the same cluster. But they had already been assigned cluster id's.
  Merge those two clusters
*/

static void
merge_clusters(int * cluster,
               const int n,
               const int c1,
               const int c2)
{
  int cfrom = cluster[c2];
  int cto = cluster[c1];

  for (int i = 0; i < n; ++i)
  {
    if (cfrom == cluster[i])
      cluster[i] = cto;
  }

  return;
}

static int
Jarvis_Patrick (Leader_Item * pool,
                const int pool_size,
                const JarvisPatrickParameters & jp,
                int * cluster)
{
  if (verbose)
  {
    cerr << "Doing JP clustering, nbrs ";
    if (0 == jp.nbr_list_size())
      cerr << '.';
    else
      cerr << jp.nbr_list_size();
    cerr << " must have " << jp.nbrs_in_common() << " in common to merge\n";
  }

  int * tmp = new_int(pool_size, -1); std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < pool_size; ++i)
  {
    pool[i].mark_neighbours(tmp, jp.nbr_list_size(), i);

    for (int j = i + 1; j < pool_size; ++j)
    {
      if (pool[j].neighbours_in_common(tmp, jp, i))
      {
//      cerr << " i = " << i << " j = " << j << " common " << pool[j].neighbours_in_common(tmp, jp, i) << " merge \n";
        merge_clusters(cluster, pool_size, i, j);
      }
    }
  }

  return 1;
}

static int
Jarvis_Patrick (Leader_Item * pool,
                const int pool_size,
                const JarvisPatrickParameters & jp,
                const resizable_array_p<IWString> & smiles,
                const resizable_array_p<IWString> & id,
                std::ostream & output)
{
  int * cluster = new int[pool_size]; std::unique_ptr<int[]> free_cluster(cluster);

  for (int i = 0; i < pool_size; ++i)
  {
    cluster[i] = i + 1;
  }

  if (! Jarvis_Patrick(pool, pool_size, jp, cluster))
    return 0;

  extending_resizable_array<int> cluster_size;

  int cluster_number = 0;

  for (int i = 0; i < pool_size; ++i)
  {
    if (0 == cluster[i])
      continue;

    int csize = std::count(cluster + i, cluster + pool_size, cluster[i]);

    write_cluster_header(i, smiles, id, cluster_number, csize, output);
    
    cluster_number++;

    for (int j = i + 1; j < pool_size; ++j)
    {
      if (cluster[i] != cluster[j])
        continue;

      cluster[j] = 0;

      write_smiles_id_dist(smiles[j], id[j], -1.0, output);
    }

    if (verbose)
      cluster_size[csize]++;

    finish_cluster(csize, output);
  }

  if (verbose)
  {
    cerr << "Clustered " << smiles.size() << " items into " << cluster_number << " clusters\n";
    int in_clusters = 0;
    for (int i = 0; i < cluster_size.number_elements(); ++i)
    {
      if (cluster_size[i] > 0)
      {
        cerr << cluster_size[i] << " clusters of size " << i << endl;
        in_clusters += (cluster_size[i]*i);
      }
    }
    cerr << in_clusters << " items in clusters\n";
  }

  return output.good();
}

static int
report_clustering (const int * selected, const resizable_array<int> & cluster_size, int pool_size, int clusters_found)
{
  cerr << "Clustered " << pool_size << " fingerprints into " << clusters_found << " clusters\n";

  int isum = 0;
  for (int i = 0; i < cluster_size.number_elements(); i++)
  {
    int j = cluster_size[i];
    if (0 == j)
      continue;

    cerr << j << " clusters of size " << i << " members\n";

    isum += j * i;
  }

  if (isum != pool_size)
    cerr << "In clusters " << isum << ", started with " << pool_size << " fingerprints\n";

  return cerr.good();
}


/*
  Prepare for Taylor Butina clustering. All we need to do is sort the pool
  by size of neighbour list
  Actually it complicated by the fact that the smiles and id arrays are separate
  We keep a separate structure describing each cluster centre and its number of
  unselected neighbours

  We don't need to sort this list, we only want the single best choice
*/

#ifdef NEED_TO_SORT_EVERYTHING
class Unselected_Item
{
  private:
    int _ndx;
    int _unselected_neighbours;
    score_t _score;

  public:

    int ndx() const { return _ndx;}
    int unselected_neighbours () const { return _unselected_neighbours;}
    score_t score () const { return _score;}

    void set_ndx (const int s) {_ndx = s;}
    void set_unselected_neighbours (const int s) {_unselected_neighbours = s;}
    void set_score (score_t s) { _score = s;}
};

class Unselected_Item_Comparator
{
  private:
  public:
    int operator () (const Unselected_Item & ui1, const Unselected_Item & ui2) const
    {
      if (ui1.unselected_neighbours() > ui2.unselected_neighbours())
        return 1;

      if (ui1.score() > ui2.score())
        return 1;

      return 0;
    }
};

static int
identify_next_tb_leader(const Leader_Item * pool,
                        const int pool_size,
                        const int * selected,
                        int & leader)
{
  Unselected_Item * ui = new Unselected_Item[pool_size]; std::unique_ptr<Unselected_Item[]> free_ui(ui);    // very inefficient to put this here, but cleaner

  int ndx;

  for (int i = 0; i < pool_size; ++i)
  {
    if (selected[i])
      continue;

    ui[ndx].set_ndx(i);
    ui[ndx].set_unselected_neighbours(pool[i].unselected_neighbours(selected));
    ui[ndx].set_score(pool[i].score());
    ndx++;
  }

  if (0 == ndx)
    return 0;

  if (1 == ndx)
    return 1;

  Unselected_Item_Comparator uic;

  gfx::timsort(ui, ui + ndx, uic);

  leader = ui[0].ndx();

  return 1;
}
#endif

static int
do_select_next_leader_randomly(const int * selected,
                               const int pool_size,
                               int & icentre)
{
  std::uniform_int_distribution<int> u(0, pool_size - 1);

  const int j = u(rng);

  for (int i = j; i < pool_size; ++i)
  {
    if (selected[i])
      continue;

    icentre = i;
    return 1;
  }

  for (int i = 0; i < j; ++i)
  {
    if (selected[i])
      continue;

    icentre = i;
    return 1;
  }

  return 0;
}

static int
identify_next_tb_leader(const Leader_Item * pool,
                        const int pool_size,
                        const int * selected,
                        const int threshold_number,
                        const similarity_type_t global_threshold,
                        int & leader)
{
  leader = -1;
  int nmax = -1;
  score_t max_score = static_cast<score_t>(0.0);

  similarity_type_t mythreshold;
  if (threshold_number >= 0)
    mythreshold = threshold[threshold_number];
  else
    mythreshold = global_threshold;

  for (int i = 0; i < pool_size; ++i)
  {
    if (selected[i])
      continue;

    const int n = pool[i].unselected_neighbours(selected, mythreshold);

    if (n < nmax)    // not interested
      ;
    else if (n > nmax)
    {
      nmax = n;
      max_score = pool[i].score();
      leader = i;
    }
    else if (! scores_present)
      ;
    else if (pool[i].score() > max_score)    // same size, break tie via score
    {
      max_score = pool[i].score();
      leader = i;
    }
  }

//#define DEBUG_IDENTIFY_NEXT_TB_LEADER
#ifdef DEBUG_IDENTIFY_NEXT_TB_LEADER
  cerr << unselected << " items avaialbel\n";
  cerr << "Max cluster size " << nmax << " item " << leader;
  if (scores_present)
    cerr << " score " << max_score;
  cerr << endl;
#endif

  return leader >= 0;
}

static int
choose_next_centre(const int * selected,
                   const Leader_Item * pool,
                   const int pool_size,
                   const int threshold_number,
                   const similarity_type_t global_threshold,
                   int & icentre)
{
  if (taylor_butina)
    return identify_next_tb_leader(pool, pool_size, selected, threshold_number, global_threshold, icentre);

  if (select_next_leader_randomly)
    return do_select_next_leader_randomly(selected, pool_size, icentre);

  icentre = -1;
  if (! scores_present)    // just grab the first unselected item
  {
    for (int i = 0; i < pool_size; i++)
    {
      if (! selected[i])
      {
        icentre = i;
        return 1;
      }
    }

    return 0;
  }

// Must look at scores

  score_t max_score = 0.0;
  for (int i = 0; i < pool_size; i++)
  {
    if (selected[i])
      continue;

    score_t s = pool[i].score();

    if (icentre < 0 || s > max_score)
    {
      max_score = s;
      icentre = i;
    }
  }

  return icentre >= 0;
}

int
leader (Leader_Item * pool,
        const int pool_size,
        int * selected,
        const resizable_array_p<IWString> & smiles,
        const resizable_array_p<IWString> & id,
        similarity_type_t global_threshold,
        const int threshold_number,
        std::ostream & output)
{
  if (verbose > 1)
    cerr << "Beginning clustering at " << global_threshold << " threshold " << threshold_number << endl;

  int clusters_found = 0;
  std::fill_n(selected, pool_size, 0);

  extending_resizable_array<int> cluster_size;

  int icentre;
  if (! choose_next_centre(selected, pool, pool_size, threshold_number, global_threshold, icentre))
  {
    cerr << "Yipes, cannot find initial leader\n";
    return 0;
  }

  resizable_array<int> cluster;

  while (1)
  {
    Leader_Item & ldr = pool[icentre];

    similarity_type_t my_threshold;
    if (threshold_number < 0)
      my_threshold = global_threshold;
    else if (ldr.threshold(my_threshold, threshold_number))     // perhaps a per item threshold can be used
      ;
    else
      my_threshold = global_threshold;

    int max_cluster_size_this_molecule;

    if (ldr.max_cluster_size(max_cluster_size_this_molecule))
      ;
    else
      max_cluster_size_this_molecule = max_cluster_size;

    if (verbose > 1)
      cerr << "Start cluster " << clusters_found << ". Threshold = " << my_threshold << ", max size " << max_cluster_size_this_molecule << endl;

    int csize;

    if (max_cluster_size_this_molecule)
      csize = pool[icentre].form_cluster_max_cluster_size(selected, my_threshold, max_cluster_size_this_molecule - 1, clusters_found, smiles, id, output);
    else
      csize = pool[icentre].form_cluster_threshold(selected, my_threshold, clusters_found, smiles, id, output);

    clusters_found++;

    if (verbose)
      cluster_size[csize]++;

    if (clusters_found >= max_clusters_to_find)
      break;

    if (! choose_next_centre(selected, pool, pool_size, threshold_number, global_threshold, icentre))
      break;
  }

  if (verbose)
    report_clustering(selected, cluster_size, pool_size, clusters_found);

  return output.good();
}

static int
leader (Leader_Item * pool,
        const int pool_size,
        int * selected,
        const resizable_array_p<IWString> & smiles,
        const resizable_array_p<IWString> & id,
        similarity_type_t global_threshold,
        int threshold_number,
        const const_IWSubstring & stem,
        int ndx)
{
  IWString fname(stem);
  fname << ndx;
  fname << ".ldr";

  std::ofstream output(fname.null_terminated_chars(), std::ios::out);
  if (! output.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }


  if (verbose)
    cerr << "Begin clustering at " << global_threshold << endl;

  return leader(pool, pool_size, selected, smiles, id, global_threshold, threshold_number, output);
}

static int
build_pool (iwstring_data_source & input,
            Leader_Item * pool,
            const int istop,
            int & ndx,
            const IW_STL_Hash_Map_int & id_to_ndx,
            resizable_array<int> & ncount)
{
  for ( ; ndx < istop && ! input.eof(); ++ndx)
  {
    if (! pool[ndx].build(input, id_to_ndx, ncount[ndx]))
    {
      cerr << "Cannot read pool item " << ndx << ", near line " << input.lines_read() << endl;
      return 0;
    }
  }

  return ndx;
}

static int
build_pool (const char * fname,
            Leader_Item * pool,
            const int istop,
            int & ndx,
            const IW_STL_Hash_Map_int & id_to_ndx,
            resizable_array<int> & ncount)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "build_pool:cannot open '" << fname << "'\n";
    return 0;
  }

  return build_pool(input, pool, istop, ndx, id_to_ndx, ncount);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Performs leader clustering, OR Jarvis Patrick clustering on a near neighbour file (from gfp_nearneighbours_single_file for example)\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " L -t <dis>         specify distance threshold(s)\n";
  cerr << " L -C <number>      maximum number of clusters to find\n";
  cerr << " L -H <TAG>         threshold for each molecule in dataitem <TAG>\n";
  cerr << " L -m <number>      maximum cluster size\n";
  cerr << " L -M <tag>         max cluster size for each molecule in <TAG>\n";
  cerr << " L -M col=<col>     max cluster size for each molecule in column <col> of name\n";
  cerr << " L -S <TAG>         score tag (molecules selected highest score to lowest)\n";
  cerr << " L -S col=<col>     scores are in column <col> of the name\n";
  cerr << " L -e            leader selections are random\n";
  cerr << endl;
  cerr << " -T                 do Taylor Butina clustering\n";
  cerr << endl;
  cerr << " J -J <n,c>         Jarvis Patrick: size of nbr list (. for no limit), neighbours in common\n";
  cerr << endl;
//cerr << " -I <TAG>         identifier tag (def PCN)\n";
  cerr << " -B <stem>        file name stem for multiple outputs\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -a               write smiles as output rather than tdt form\n";
  cerr << " -h               suppress self neighbours during Taylor Butina and Leader\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
do_jarvis_patrick (const JarvisPatrickParameters * jp,
                   const int njp,
                   Leader_Item * pool,
                   const int pool_size,
                   const resizable_array_p<IWString> & smiles,
                   const resizable_array_p<IWString> & id,
                   const IWString & stem)
{
  if (1 == njp)
  {
    if (0 == stem.length())
      return Jarvis_Patrick(pool, pool_size, jp[0], smiles, id, cout);

    IWString fname;
    fname << stem << ".jpc";

    std::ofstream output(fname.null_terminated_chars());

    if (! output.good())
    {
      cerr << "Cannot open output file '" << fname << "'\n";
      return 0;
    }

    return Jarvis_Patrick(pool, pool_size, jp[0], smiles, id, output);
  }

  for (int i = 0; i < njp; ++i)
  {
    IWString fname;

    fname << stem << i << ".jpc";

    std::ofstream output(fname.null_terminated_chars());

    if (! output.good())
    {
      cerr << "Cannot open output file '" << fname << "'\n";
      return 0;
    }

    if (! Jarvis_Patrick(pool, pool_size, jp[i], smiles, id, output))
      return 0;
  }

  return 1;
}

static int
do_leader_clustering (Leader_Item * pool,
                      const int pool_size,
                      const resizable_array_p<IWString> & smiles,
                      const resizable_array_p<IWString> & id,
                      const IWString & stem,
                      const float * threshold,
                      const int number_thresholds)
{
  int * selected = new int[pool_size]; std::unique_ptr<int[]> free_selected(selected);

  int computation_number = 0;

  if (1 == number_thresholds && 0 == number_threshold_tags)
    leader(pool, pool_size, selected, smiles, id, threshold[0], -1, cout);
  else
  {
    for (int i = 0; i < number_thresholds; i++)
    {
      if (! leader(pool, pool_size, selected, smiles, id, threshold[i], -1, stem, computation_number))
      {
        cerr << "Clustering at " << threshold[i] << " failed\n";
        return 0;
      }

      computation_number++;
    }
  }

  if (1 == number_threshold_tags && 0 == number_thresholds)
    leader(pool, pool_size, selected, smiles, id, 0.0, 0, cout);
  else
  {
    for (int i = 0; i < number_threshold_tags; i++)
    {
      if (! leader(pool, pool_size, selected, smiles, id, 0.0, i, stem, computation_number))
      {
        cerr << "Clustering on threshold tag " << i << " failed\n";
        return 0;
      }

      computation_number++;
    }
  }

  return 1;
}

static int
read_scores_from_file (Leader_Item * pool,
                       const int pool_size,
                       const IW_STL_Hash_Map_int & id_to_ndx,
                       iwstring_data_source & input)
{
  const_IWSubstring buffer;

  int scores_processed = 0;

  while (input.next_record(buffer))
  {
    IWString id;
    int i = 0;
    buffer.nextword(id, i);

    const auto f = id_to_ndx.find(id);

//  cerr << "Processing '" << id << "' absent? " << (f == id_to_ndx.end()) << endl;

    if (f == id_to_ndx.end())   // extra data, ignore
      continue;

    const_IWSubstring token;
    if (! buffer.nextword(token, i))
    {
      cerr << "read_scores_from_file:cannot read second column from '" << buffer << "'\n";
      return 0;
    }

    score_t s;
    
    if (token.numeric_value(s))
    {
      pool[f->second].set_score(s);
      scores_processed++;
    }
    else if (0 == scores_processed)
      ;
    else
    {
      cerr << "read_scores_from_file:invalid numeric '" << buffer << "'\n";
      return 0;
    }
  }

  if (pool_size == scores_processed)
    return 1;

  cerr << "read_scores_from_file:pool contains " << pool_size << " fingerprints, but read " << scores_processed << " scores\n";
  return 0;
}

static int
read_scores_from_file (Leader_Item * pool,
                       const int pool_size,
                       const IW_STL_Hash_Map_int & id_to_ndx,
                       IWString & fname)
{
  iwstring_data_source input(fname.null_terminated_chars());

  if (! input.good())
  {
    cerr << "read_scores_from_file:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_scores_from_file(pool, pool_size, id_to_ndx, input);
}

static int
leader (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:I:t:H:S:C:m:M:B:J:aThe");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('a'))
  {
    write_tdt_output = 0;

    if (verbose)
      cerr << "Smiles output\n";
  }

  int pool_size = 0;
  Leader_Item * pool = nullptr;

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    if (verbose)
      cerr << "system sized to " << pool_size << endl;
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

  if (cl.option_present('I'))
  {
    (void) cl.value('I', identifier_tag);

    if (verbose)
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
  }

  if (cl.option_present('H'))
  {
    number_threshold_tags = cl.option_count('H');
    threshold_from_file_tag = new IWString[number_threshold_tags];

    for (int i = 0; i < number_threshold_tags; i++)
    {
      threshold_from_file_tag[i] = cl.string_value('H', i);

      if (verbose)
        cerr << "Threshold in tag '" << threshold_from_file_tag[i] << "' tag\n";
    }
  }

  IWString score_file;

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (s.starts_with("col="))
    {
      s.remove_leading_chars(4);
      if (! s.numeric_value(score_column) || score_column < 1)
      {
        cerr << "INvalid score column specifier 'col=" << s << "'\n";
        return 8;
      }

      if (verbose)
        cerr << "Scores taken from column " << score_column << " of the identifier\n";

      score_column--;
    }
    else if (s.starts_with("FILE="))
    {
      s.remove_leading_chars(5);
      score_file = s;
    }
    else
    {
      score_tag = cl.string_value('S');
      if (verbose)
        cerr << "Score tag is " << score_tag << "'\n";
    }

    scores_present = 1;
  }

  if (cl.option_present('J'))
    ;
  else if (! cl.option_present('t') && ! cl.option_present('H') && ! cl.option_present('Y'))
  {
    cerr << "Threshold distance must be specified via -t, -H or -V options\n";
    usage(28);
  }

  number_thresholds = cl.option_count('t');

  if (1 == number_thresholds)    // maybe it is of the form 0.1,0.2,
  {
    const_IWSubstring t = cl.string_value('t');

    number_thresholds = t.nwords(',');
    threshold = new similarity_type_t[number_thresholds];

    int i = 0;
    const_IWSubstring token;
    for (int ndx = 0; t.nextword(token, i, ','); ndx++)
    {
      if (! token.numeric_value(threshold[ndx]) || threshold[ndx] <= 0.0f || threshold[ndx] >= 1.0f)
      {
        cerr << "INvalud threshold '" << token << "' in '" << t << "'\n";
        return 1;
      }

      if (verbose)
        cerr << "Threshold " << ndx << " set to " << threshold[ndx] << endl;
    }
  }
  else if (number_thresholds > 1)
  {
    threshold = new similarity_type_t[number_thresholds];
    for (int i = 0; i < number_thresholds; i++)
    {
      if (! cl.value('t', threshold[i], i) || threshold[i] < 0.0 || threshold[i] > 1.0)
      {
        cerr << "Invalid threshold\n";
        return 6;
      }

      if (verbose)
        cerr << "Threshold " << i << " set to " << threshold[i] << endl;
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

  if (cl.option_present('e'))
  {
    select_next_leader_randomly = 1;

    if (verbose)
      cerr << "Next leader selected randomly\n";
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

  JarvisPatrickParameters * jp = nullptr;
  int njp = 0;

  if (cl.option_present('J'))
  {
    njp = cl.option_count('J');
    jp = new JarvisPatrickParameters[njp];

    const_IWSubstring j;
    for (int i = 0; cl.value('J', j, i); ++i)
    {
      if (! jp[i].build(j))
      {
        cerr << "Cannot initialise Jarvis Patrick parameters '" << j << "'\n";
        return 2;
      }
    }

    if (verbose)
      cerr << "Recognised " << njp << " Jarvis Patrick clustering specifications\n";
  }

  if (cl.option_present('T'))
  {
    if (nullptr != jp)
    {
      cerr << "Cannot do Jarvis Patrick and Taylor Butina\n";
      return 1;
    }

    taylor_butina = 1;

    if (verbose)
      cerr << "Will do Taylor Butina leader clustering\n";
  }

  if (cl.option_present('h'))
  {
    suppress_self_neighbours = 1;

    if (verbose)
      cerr << "Will suppress self neighbours\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  int max_pool_size;
  if (pool_size > 0)
    max_pool_size = pool_size;
  else
    max_pool_size = std::numeric_limits<int>::max();

  resizable_array_p<IWString> smiles, id;
  resizable_array<int> ncount;
  resizable_array<int> istop;

  for (int i = 0; i < cl.number_elements(); ++i)
  {
    if (! read_smiles_and_ids_and_nbr_count(cl[i], max_pool_size, smiles, id, ncount))
    {
      cerr << "Cannot extract identifiers from " << cl[i] << "'\n";
      return i + 1;
    }

    if (verbose > 1)
      cerr << "After reading " << cl[i] << " have " << smiles.size() << " molecules\n";

    istop.add(smiles.number_elements());
  }

  pool_size = smiles.number_elements();

  IW_STL_Hash_Map_int id_to_ndx;

  for (int i = 0; i < pool_size; ++i)
  {
    IWString s = *(id[i]);

    s.truncate_at_first(' ');

    if (id_to_ndx.contains(s))
    {
      cerr << "Duplicate identifier '" << s << "' index " << i << ", cannot process\n";
      return 1;
    }

    id_to_ndx[s] = i;
  }

  pool = new Leader_Item[pool_size]; std::unique_ptr<Leader_Item[]> free_pool(pool);
  for (int i = 0; i < pool_size; ++i)
  {
    pool[i].set_ndx(i);
  }

  int ndx = 0;

  for (int i = 0; i < cl.number_elements(); ++i)
  {
    if (! build_pool(cl[i], pool, istop[i], ndx, id_to_ndx, ncount))
    {
      cerr << "Cannot read fingerprints from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  if (score_file.length())
  {
    if (! read_scores_from_file(pool, pool_size, id_to_ndx, score_file))
    {
      cerr << "Cannot read scores from '" << score_file << "'\n";
      return 2;
    }
  }

  float max_distance = 0.0f;
  Accumulator_Int<int> nbr_count;

  for (int i = 0; i < pool_size; ++i)
  {
    if (verbose)
    {
      pool[i].update_max_distance(max_distance);
      nbr_count.extra(pool[i].neighbour_count());
    }
  }

  if (verbose)
  {
    cerr << "Read " << pool_size << " NN data sets. Max distance " << max_distance << ". ";
    if (nbr_count.minval() == nbr_count.maxval())
      cerr << nbr_count.minval() << " neighbours";
    else
      cerr << "Neighbours btw " << nbr_count.minval() << " and " << nbr_count.maxval() << " ave " << static_cast<float>(nbr_count.average());
    cerr << endl;
  }

  int number_clusterings = number_thresholds + number_threshold_tags;
  if (number_clusterings > 1 && ! cl.option_present('B'))
  {
    cerr << "If doing multiple thresholds, must specify file name stem via the -B option\n";
    usage(12);
  }

  if (njp > 1 && ! cl.option_present('B'))
  {
    cerr << "If doing multiple JP clusters, must specify file name stem via the -B option\n";
    usage(12);
  }

  IWString stem;
  if (cl.option_present('B'))
  {
    cl.value('B', stem);

    if (verbose)
      cerr << "Output file(s) created with stem '" << stem << "'\n";
  }

  int rc = 0;

  if (njp)
  {
    if (! do_jarvis_patrick(jp, njp, pool, pool_size, smiles, id, stem))
      rc = 1;

    delete [] jp;
  }
  else
  {
    if (! do_leader_clustering(pool, pool_size, smiles, id, stem, threshold, number_thresholds))
      rc = 1;

    if (nullptr != threshold)
      delete [] threshold;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = leader(argc, argv);

  return rc;
}
