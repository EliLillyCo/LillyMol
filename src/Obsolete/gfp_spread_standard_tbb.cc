/*
  Spread implementation
*/

#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

#include <iostream>
#include <random>

#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for.h"
#include "tbb/scalable_allocator.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/gfp_standard.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int brief_output = 0;

static Report_Progress report_progress;

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

static int squeeze = 0;
static int next_squeeze = 0;

static int grainsize = 20000;

static int
calculate_grainsize (int zrange)
{
  if (zrange <= 0)   // should not happen
    return 1;

  if (zrange < 2000)  // too small to multi-thread
    return zrange;

  if (zrange <= grainsize)
    return zrange / 5;    // 5 is arbitrary

  return grainsize;
}

/*
  When we squeeze the pool, we must update the sdistance and nsn arrays, but
  we do not want to update the smiles and pcn arrays
  So each item must keep track of their initial index into those arrays
*/

class SSpread_Item : public GFP_Standard
{
  private:
    int _initial_ndx;

  public:
    SSpread_Item ();

    void set_initial_ndx (int s) { _initial_ndx = s;}
    int initial_ndx () const { return _initial_ndx;}

    const IWString & smiles (const IWString * smiles) const { return smiles[_initial_ndx];}
    const IWString & pcn (const IWString * pcn) const { return pcn[_initial_ndx];}
};

SSpread_Item::SSpread_Item()
{
  _initial_ndx = -1;
}

struct Selected_Item
{
  int _sel;
  int _nsn;
  float _dist;
};

static void
set_selected_item (Selected_Item & s,
                   const int isel,
                   const int nsn,
                   const float d)
{
  s._sel = isel;
  s._nsn = nsn;
  s._dist = d;

  return;
}

static SSpread_Item * pool = nullptr;
static int pool_size = 0;
static float * sdistance = nullptr;
static IWString * smiles = nullptr;
static IWString * pcn = nullptr;
static int * nsn = nullptr;
static Selected_Item * selected_item = nullptr;

template <typename F>
int
build_pool (iwstring_data_source & input,
            F * & pool,
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

template <typename F>
int
build_pool (const char * fname,
            F * & pool,
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

    pool_size = input.count_records_starting_with("PCN<");

    if (0 == pool_size)
    {
      cerr << "Zero occurrences of '" << "PCN<" << "' in '" << fname << "'\n";
      return 0;
    }
  }

  pool = new F[pool_size];
  smiles = new IWString[pool_size];
  pcn = new IWString[pool_size];
  selected_item = new Selected_Item[pool_size];

  return build_pool(input, pool, pool_size);
}

template <typename F>
void
do_output (const F * pool,
           const int isel,
           IWString_and_File_Descriptor & output)
{
  const int ini = pool[isel].initial_ndx();

  if (brief_output)
    output << smiles[ini] << ' ' << pcn[ini] << ' ' << sdistance[isel] << "\n";
  else
  {
    output << smiles_tag     << smiles[isel] << ">\n";
    output << identifier_tag << pcn[isel] << ">\n";

    if (nsn[isel] >= 0)
    {
      const int n = pool[nsn[isel]].initial_ndx();

      output << smiles_tag     << smiles[n] << ">\n";
      output << identifier_tag << pcn[n] << ">\n";
      output << distance_tag   << sdistance[isel] << ">\n";
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

template <typename F>
class Propagate_Selected_Item_TBB
{
  private:
    const int _isel;
    const F _fpsel;
    float _furthest_distance;
    int _most_distant;

  public:
    Propagate_Selected_Item_TBB (int i, const F & fp) : _isel (i), _fpsel(fp), _furthest_distance(-1.0f), _most_distant(-1) {}

    float furthest_distance () const { return _furthest_distance;}
    int   most_distant      () const { return _most_distant;}

    Propagate_Selected_Item_TBB (const Propagate_Selected_Item_TBB & rhs, tbb::split)  : _isel(rhs._isel), _fpsel(rhs._fpsel)
    {
      _furthest_distance = rhs._furthest_distance;
      _most_distant = rhs._most_distant;
    }

    void operator()( const tbb::blocked_range<int>& r )
    {
//    cerr << "Scanning from " << r.begin() << " to " << r.end() << " due to " << _isel << endl;
      const auto & rend = r.end();
      for (auto i = r.begin(); i != rend; ++i)
      {
        if (sdistance[i] < 0.0f || i == _isel)
          continue;

        float d = 1.0f - _fpsel.tanimoto(pool[i]);
//      cerr << i << " distance currently " << sdistance[i] << " to new fp " << d << ", furthest is " << _furthest_distance << endl;
        if (d < sdistance[i])
        {
          sdistance[i] = d;
          nsn[i] = _isel;
//        cerr << i << " distance to nearest selected updated to " << sdistance[i] << endl;
        }

        if (sdistance[i] > _furthest_distance)
        {
          _furthest_distance = sdistance[i];
          _most_distant = i;
//        cerr << i << " furthest distance updated to " << _furthest_distance << " because of " << i << endl;
        }
      }

//    cerr << "At end of scan due to " << _isel << " most distant is " << _most_distant << " at dist " << _furthest_distance << endl;
    }

    void join (const Propagate_Selected_Item_TBB & rhs) 
    {
//    cerr << "Join " << _furthest_distance << " with " << rhs._furthest_distance << endl;
      if (rhs._furthest_distance > _furthest_distance)
      {
        _furthest_distance = rhs._furthest_distance;
        _most_distant = rhs._most_distant;
      }
    }
};

template <typename F>
void
do_squeeze (F * pool,
            int & pool_size,
            float * sdistance,
            int * nsn)
{
  if (pool_size < 1000)
    return;

  int ndx = 0;
  for (int i = 0; i < pool_size; ++i)
  {
    if (sdistance[i] < 0.0f)
      continue;

    if (i == ndx)
    {
      ndx++;
      continue;
    }

    pool[ndx] = pool[i];
    sdistance[ndx] = sdistance[i];
    nsn[ndx] = nsn[i];
    ndx++;
  }

  pool_size = ndx + 1;

  return;
}

template <typename F>
int
fpobj_spread (F * pool,
              int pool_size,
              int number_to_select,
              Selected_Item * selected_item)
{
  int first_selected = 0;

  sdistance[first_selected] = -1.0f;

  set_selected_item(selected_item[0], first_selected, -1, 2.0f);

  int number_selected = 1;

  int ichoose = first_selected;

  while (number_selected < number_to_select)
  {
    Propagate_Selected_Item_TBB<F> psitbb(ichoose, pool[ichoose]);
    tbb::parallel_reduce(tbb::blocked_range<int>(0, pool_size, calculate_grainsize(pool_size)), psitbb);//,tbb::auto_partitioner());

    ichoose = psitbb.most_distant();

    assert (ichoose >= 0);

    set_selected_item(selected_item[number_selected], ichoose, nsn[ichoose], sdistance[ichoose]);

    number_selected++;

    if (report_progress())
      cerr << "Selected " << number_selected << " '" << smiles[ichoose] << "' distance " << sdistance[ichoose] << " NSN '" << pcn[nsn[ichoose]] << "'\n";

    sdistance[ichoose] = -1.0f;

    if (number_selected == next_squeeze)
    {
      do_squeeze(pool, pool_size, sdistance, nsn);
      next_squeeze += squeeze;
      if (verbose > 1)
        cerr << "Squeezed, pool now " << pool_size << endl;
    }
  }

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
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
fpobj_spread (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:h:bg:r:q:");

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

  tbb::task_scheduler_init * init = nullptr;

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

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(2);
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

  for (int i = 0; i < pool_size; ++i)
  {
    pool[i].set_initial_ndx(i);
  }

  sdistance = new_float(pool_size, 2.0f);
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

  if (cl.option_present('q'))
  {
    if (! cl.value('q', squeeze) || squeeze < 0)
    {
      cerr << "The squeeze every option (-q) must be a whole non negative number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will squeeze selected items every " << squeeze << " items selected\n";

    next_squeeze = squeeze;
  }
  else
  {
    squeeze = 1000;
    next_squeeze = 1000;
  }

  const int nsel = fpobj_spread (pool, pool_size, number_to_select, selected_item);

  IWString_and_File_Descriptor output(1);
  
  if (cl.option_present('b'))
  {
    for (int i = 0; i < nsel; ++i)
    {
      const Selected_Item & s = selected_item[i];
      output << smiles[s._sel] << ' ' << pcn[s._sel] << ' ' << s._dist << '\n';

      output.write_if_buffer_holds_more_than(4096);
     }
  }
  else
  {
    for (int i = 0; i < nsel; ++i)
    {
      const Selected_Item & s = selected_item[i];

      output << smiles_tag     << smiles[s._sel] << ">\n";
      output << identifier_tag << pcn[s._sel]    << ">\n";
      if (s._nsn >= 0)
      {
        output << smiles_tag     << smiles[s._nsn] << ">\n";
        output << identifier_tag << pcn[s._nsn]    << ">\n";
        output << distance_tag   << s._dist       << ">\n";
      }
      else
      {
        output << smiles_tag << "*>\n";
        output << identifier_tag << "*>\n";
        output << distance_tag << "1>\n";
      }
      output << "|\n";

      output.write_if_buffer_holds_more_than(4096);
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
