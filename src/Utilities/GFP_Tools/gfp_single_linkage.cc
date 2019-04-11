/*
  Single linkage clustering of a gfp file
*/

#include <stdlib.h>
#include <limits>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "cmdline.h"
#include "iwqsort.h"
#include "iwstring_data_source.h"
#include "iw_tdt.h"

#include "gfp.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag("DIST<");
static IWString csize_tag("CSIZE<");
static IWString maxdist_tag("MAXDIST<");
static IWString purity_tag("PURE<");

static similarity_type_t threshold;

static int centroid_first = 0;

static extending_resizable_array<int> cluster_size;
static extending_resizable_array<int> iterations;

static int clusters_formed = 0;
static int items_selected = 0;

static int next_time_to_squeeze = std::numeric_limits<int>::max();

static int squeeze_selected_every = 0;

/*
  We check distances every second expansion. If on iteration I
  a point was at least 2xthreshold away from every point in
  that expansion, it cannot be in the next shell.
  But 2x is quite conservative, since distances seldom go
  beyond 0.6. Leave this user definable.
*/

static float expansion_ratio = static_cast<float>(2.0);

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Single linkage clustering\n";
  cerr << " -t <rad>       threshold for grouping\n";
  cerr << " -c             sort clusters by distance to threshold\n";
//cerr << " -p             identify pure clusters\n";
  cerr << " -s <number>    squeeze out already selected items every <n> items selected\n";
  cerr << " -r <float>     per shell radius distance expansion factor (def 2.0)\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Single_Linkage_Object : public IW_General_Fingerprint
{
  private:
    IWString _smiles;

    int _selected;

//  Each item needs to keep track of the sum of distances to other
//  cluster members, and the longest distance encountered

    similarity_type_t _sumd;
    similarity_type_t _maxd;

    similarity_type_t _dist_to_current_cluster;

  public:
    Single_Linkage_Object ();

    int construct_from_tdt (IW_TDT &);

    void set_selected (int s) { _selected = s;}
    int  selected () const { return _selected;}

    void reset () { _sumd = static_cast<float>(0.0); _maxd = static_cast<float>(0.0);}

    void extra (similarity_type_t s) { _sumd += s; if (s > _maxd) _maxd = s;}

    similarity_type_t dist_to_current_cluster () const { return _dist_to_current_cluster;}
    void set_dist_to_current_cluster(similarity_type_t s) { _dist_to_current_cluster = s;}
    void decrement_dist_to_current_cluster(similarity_type_t s) { _dist_to_current_cluster -= s;}

    void notify_dist_to_current_cluster (similarity_type_t s);

    float sumd () const { return _sumd;}
    void  set_sumd (float f) { _sumd = f;}

    float maxd () const { return _maxd;}

    const IWString & smiles () const { return _smiles;}
};

#ifdef __GNUG__
template class resizable_array<Single_Linkage_Object *>;
#endif

Single_Linkage_Object::Single_Linkage_Object ()
{
  _selected = 0;

  return;
}

void
Single_Linkage_Object::notify_dist_to_current_cluster (similarity_type_t s)
{
  if (s < _dist_to_current_cluster)
    _dist_to_current_cluster = s;

  return;
}

class Single_Linkage_Object_Comparator
{
  private:
  public:
    int operator () (const Single_Linkage_Object *, const Single_Linkage_Object *) const;
};

int
Single_Linkage_Object_Comparator::operator () (const Single_Linkage_Object * s1, const Single_Linkage_Object * s2) const
{
  if (s1->sumd() < s2->sumd())
    return -1;
  else if (s1->sumd() > s2->sumd())
    return 1;
  else
    return 0;
}

int
Single_Linkage_Object::construct_from_tdt (IW_TDT & tdt)
{
  int fatal;

  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal))
    return 0;

  if (! tdt.dataitem_value(smiles_tag, _smiles))
  {
    cerr << "Single_Linkage_Object::construct_from_tdt:no smiles\n";
    return 0;
  }

  return 1;
}

static Single_Linkage_Object ** pool = NULL;
static int pool_size = 0;
static int initial_pool_size = 0;

static int
squeeze_out_selected_items()
{
  int ndx = 0;
  for (int i = 0; i < pool_size; i++)
  {
    Single_Linkage_Object * pi = pool[i];

    if (! pi->selected())
    {
      pool[ndx] = pi;
      ndx++;
    }
  }

  pool_size = ndx;

  if (verbose > 1)
    cerr << "Pool squeezed to " << ndx << " active items\n";

  return ndx;
}

static int
write_smiles_and_id (const Single_Linkage_Object * f,
                     IWString_and_File_Descriptor & output)
{
  output << smiles_tag << f->smiles() << ">\n";
  output << identifier_tag << f->id() << ">\n";

  return 1;
}

static int
write_smiles_id_distance (const Single_Linkage_Object * f,
                          IWString_and_File_Descriptor & output)
{
  write_smiles_and_id (f, output);

  output << distance_tag << f->sumd() << ">\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
move_centroid_to_first (resizable_array<Single_Linkage_Object *> & cluster,
                        IWString_and_File_Descriptor & output)
{
  int n = cluster.number_elements();

  for (int i = 0; i < n; i++)
  {
    cluster[i]->reset();
  }

  write_smiles_and_id (cluster[0], output);
  output << csize_tag << n << ">\n";

  if (n > 2)
  {
    similarity_type_t largest_pairwise_distance = static_cast<float>(0.0);

    for (int i = 0; i < n; i++)
    {
      Single_Linkage_Object * pi = cluster[i];

      for (int j = i + 1; j < n; j++)
      {
        Single_Linkage_Object * pj = cluster[j];

        similarity_type_t d = static_cast<similarity_type_t>(1.0) - pi->tanimoto(*pj);

        pi->extra(d);
        pj->extra(d);

        if (d > largest_pairwise_distance)
          largest_pairwise_distance = d;
      }
    }

    int purity_count = 0;
    for (int i = 0; i < n; i++)
    {
      if (cluster[i]->maxd() <= threshold)
        purity_count++;
    }

    Single_Linkage_Object_Comparator sloc;
    cluster.iwqsort(sloc);

    output << purity_tag << purity_count << ">\n";
    output << maxdist_tag << largest_pairwise_distance << ">\n";
  }

  for (int i = 1; i < n; i++)
  {
    Single_Linkage_Object * ci = cluster[i];

    float f = ci->sumd() / static_cast<float>(n);   // normalise to ave dist

    ci->set_sumd(f);

    write_smiles_id_distance (cluster[i], output);
  }

  output << "|\n";

  return 1;
}

static int
process_cluster (resizable_array<Single_Linkage_Object *> & cluster,
                 IWString_and_File_Descriptor & output)
{
  int n = cluster.number_elements();

  clusters_formed++;
  cluster_size[n]++;
  items_selected += n;

  if (verbose > 1)
    cerr << "Formed cluster " << clusters_formed << ", size " << n << ", nsel " << items_selected << endl;

  if (centroid_first)
    return move_centroid_to_first (cluster, output);

  output << cluster[0]->id();

  for (int i = 1; i < n; i++)
  {
    output << ' ' << cluster[i]->id();
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static void
reset_all_distance_to_current_cluster_values (Single_Linkage_Object ** pool,
                                int pool_size)
{
  for (int j = 0; j < pool_size; j++)
  {
    pool[j]->set_dist_to_current_cluster(static_cast<similarity_type_t>(1.0));
  }

  return;
}

static int
gfp_single_linkage (Single_Linkage_Object * ldr,
                    IWString_and_File_Descriptor & output)
{
  static resizable_array<Single_Linkage_Object *> cluster;
  static resizable_array<Single_Linkage_Object *> next_shell;

  if (0 == cluster.number_elements())
  {
    cluster.resize(pool_size);
    next_shell.resize(pool_size);
  }
  else
  {
    cluster.resize_keep_storage(0);
    next_shell.resize_keep_storage(0);
  }

  cluster.add(ldr);
  next_shell.add(ldr);
  ldr->set_selected(clusters_formed + 1);

#define EXPANSION_RATIO_ACTIVE
#ifdef EXPANSION_RATIO_ACTIVE
  for (int i = 0; i < pool_size; i++)
  {
    pool[i]->set_dist_to_current_cluster(static_cast<similarity_type_t>(0.0));
  }
#endif

// We alternate between having to compute all distances and checking them

  int need_to_compute = 1;

  int iterations_this_cluster = 0;

  while (next_shell.number_elements())
  {
    int n = next_shell.number_elements();

#ifdef EXPANSION_RATIO_ACTIVE
    if (need_to_compute)
      reset_all_distance_to_current_cluster_values(pool, pool_size);
#endif

    for (int i = 0; i < n; i++)
    {
      Single_Linkage_Object * ni = next_shell[i];

      for (int j = 0; j < pool_size; j++)
      {
        Single_Linkage_Object * pj = pool[j];

        if (pj->selected())
          continue;

#ifdef EXPANSION_RATIO_ACTIVE
        if (need_to_compute)
          ;
        else if (pj->dist_to_current_cluster() > expansion_ratio * threshold)
          continue;
#endif

        if (! can_be_compared (*ni, *pj))
          continue;

        similarity_type_t d = static_cast<similarity_type_t>(1.0) - ni->tanimoto(*pj);

        if (d > threshold)
        {
          pj->notify_dist_to_current_cluster(d);
          continue;
        }

        cluster.add(pj);
        next_shell.add(pj);
        pj->set_selected(clusters_formed + 1);
      }
    }

    next_shell.erase (0, n - 1);   // remove the first N to leave only what we just added

    iterations_this_cluster++;

    need_to_compute = ! need_to_compute;
  }

  iterations[iterations_this_cluster]++;

  return process_cluster (cluster, output);
}

static int
gfp_single_linkage (IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < pool_size; i++)
  {
    Single_Linkage_Object * pi = pool[i];

    if(pi->selected())
      continue;

    gfp_single_linkage (pi, output);

    if (items_selected > next_time_to_squeeze)
    {
      squeeze_out_selected_items();
      next_time_to_squeeze = items_selected + squeeze_selected_every;
    }
  }

  return 1;
}


static int
build_pool (iwstring_data_source & input)
{
  int ndx = 0;
  IW_TDT tdt;
  while (tdt.next(input))
  {
    Single_Linkage_Object * s = new Single_Linkage_Object;

    if (s->construct_from_tdt (tdt))
      pool[ndx++] = s;
    else
    {
      delete s;
      return 0;
    }
  }

  pool_size = ndx;
  initial_pool_size = pool_size;

  return pool_size;
}

static int
build_pool (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  assert (NULL == pool);

  IWString tmp = '^';
  tmp << identifier_tag;

  pool_size = input.grep (tmp);
  if (0 == pool_size)
  {
    cerr << "Yipes, cannot find any '" << tmp << "' in the input\n";
    return 0;
  }

  pool = new Single_Linkage_Object *[pool_size];

  if (verbose)
    cerr << "Input contains " << pool_size << " fingerprints\n";

  return build_pool (input);
}

static int
gfp_single_linkage (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:P:G:t:W:cs:r:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('t'))
  {
    cerr << "Must specify clustering threshold via the -t option\n";
    usage(3);
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', threshold) || threshold < 0.0 || threshold >= 1.0)
    {
      cerr << "The threshold (-t) must be a valid distance\n";
      return 4;
    }

    if (verbose)
      cerr << "Threshold set to " << threshold << endl;
  }

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise GFP options\n";
      usage (23);
    }
  }

  if (cl.option_present('c'))
  {
    centroid_first = 1;

    if (verbose)
      cerr << "Will move the centroid to the beginning of each cluster\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', squeeze_selected_every) || squeeze_selected_every < 1)
    {
      cerr << "The squeeze every option (-s) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will squeeze out already selected items every " << squeeze_selected_every << " items selected\n";

    next_time_to_squeeze = squeeze_selected_every;
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', expansion_ratio) || expansion_ratio < 1.0)
    {
      cerr << "The shell radius expansion factor (-r) must be a number >= 1.0\n";
      usage(3);
    }

    if (verbose)
      cerr << "Shell radius expansion factor set to " << expansion_ratio << endl;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Only one command line argument allowed\n";
    return 3;
  }

  if (! build_pool (cl[0]))
  {
    cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
    return 3;
  }

  if (verbose)
    cerr << "Read " << pool_size << " fingerprints from '" <<cl[0] << "'\n";

  IWString_and_File_Descriptor output(1);

  gfp_single_linkage (output);

  output.flush();

  if (verbose)
  {
    cerr << "Grouped " << initial_pool_size << " items into " << clusters_formed << " clusters\n";
    for (int i = 0; i < cluster_size; i++)
    {
      if (cluster_size[i])
        cerr << cluster_size[i] << " clusters of size " << i << endl;
    }

    for (int i = 0; i < iterations.number_elements(); i++)
    {
      if (iterations[i])
        cerr << iterations[i] << " clusters took " << i << " expansion iterations\n";
    }
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_single_linkage(argc, argv);

  return rc;
}
