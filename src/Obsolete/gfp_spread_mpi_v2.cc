#include <stdlib.h>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <fstream>
#include <ctime>
#include <memory>

# include "mpi.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "gfp_standard.h"

using namespace std;

static int verbose = 0;

#define SPREAD_THREAD_DONE -2776704

static Report_Progress report_progress;

static int squeeze_selected = 0;

void timestamp ( );

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "MPI implementation of spread. Input is a single gfp file.\n";
  cerr << " -s <n>         number of fingerprints in file - must be specified\n";
  cerr << " -n <items>     number of items to select\n";
  cerr << " -L <fname>     use a log file\n";
  cerr << " -r <n>         report progress every <n> items selected\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Spread_Item : public GFP_Standard
{
  private:
    int _nbr;           // previously selected nearest nbr

  public:
    Spread_Item();
    ~Spread_Item();

    int build (IW_TDT & tdt);

    int debug_print (ostream &) const;

    int nbr () const { return _nbr;}

    void set_nbr    (int n) { _nbr = n;}

    int update_distance_if_needed (const Spread_Item & sel, int);
};

Spread_Item::Spread_Item()
{
  _nbr = -1;

  return;
}

Spread_Item::~Spread_Item()
{
  return;
}

int
Spread_Item::debug_print (ostream & output) const
{
  return 1;
}

int
Spread_Item::build (IW_TDT & tdt)
{
  IW_General_Fingerprint gfp;

  int fatal;
   if (! gfp.construct_from_tdt(tdt, fatal))
   {
     cerr << "Cannot read fingerprint\n";
     return 0;
   }

  build_molecular_properties(gfp.molecular_properties_integer());
  build_iw(gfp[0]);
  build_mk(gfp[1]);
  build_mk2(gfp[2]);

  return 1;
}

/*
  The files need to pass back to the controlling thread information about
  what their longest distance is, as well as the item that has that longest
  distance
*/

struct Item_and_Distance
{
  int _item;
  float _distance;
};

/*
  The fingerprint pool retains a vector of the active fingerprints, which means
  more efficient memory usage when querying the pool
*/

class Fingerprint_Pool
{
  private:
    Spread_Item * _pool;
    int _n;

    int _fstart;
    int _fstop;

    float * _distance;    // negative distances will signify already selected

    ofstream _logfile;

//  private functions

    int _read_pool (iwstring_data_source &);
    void _identify_unselected_item_with_longest_distance (struct Item_and_Distance & iad) const;

  public:
    Fingerprint_Pool();
    ~Fingerprint_Pool();

    int read_pool (const char * fname, int s, int a, int o);

    void wait_loop ();

    int initialise_log_file (const IWString &);
};

Fingerprint_Pool::Fingerprint_Pool ()
{
  _pool = nullptr;
  _n = 0;

  _distance = nullptr;

  return;
}

Fingerprint_Pool::~Fingerprint_Pool ()
{
  if (NULL != _pool)
    delete [] _pool;

  if (NULL != _distance)
    delete [] _distance;

  return;
}

int
Fingerprint_Pool::_read_pool (iwstring_data_source & input)
{
  int ndx = 0;

  IW_TDT tdt;

  while (tdt.next(input))
  {
    if (! _pool[ndx].build(tdt))
    {
      cerr << "Cannot build fingerprint\n";
      return 0;
    }

    ndx++;
    if (ndx >= _n)
      break;
  }

  if (ndx != _n)
  {
    cerr << "Huh, expected " << _n << " fingerprints, but read " << ndx << endl;
    return 0;
  }

  return 1;
}

int
Fingerprint_Pool::read_pool (const char * fname,
                             int pool_size,
                             int a, 
                             int o)
{
  assert (NULL == _pool);
  assert (a < o);
  assert (o <= pool_size);

  cerr << "Processor building, range " << a << " to " << o << endl;

  _n = pool_size;

  _pool = new Spread_Item[_n];
  _distance = new_float(_n, 2.0f);

  _fstart = a;
  _fstop = o;

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return _read_pool (input);
}

int
Fingerprint_Pool::initialise_log_file (const IWString & stem)
{
  IWString fname(stem);

  const auto r = MPI::COMM_WORLD.Get_rank ( );

  fname << r;

  _logfile.open(fname.null_terminated_chars());

  if (! _logfile.good())
  {
    cerr << "Fingerprint_Pool::initialise_log_file:cannot open log file '" << fname << "'\n";
    return 0;
  }

  return 1;
}

void
Fingerprint_Pool::_identify_unselected_item_with_longest_distance (Item_and_Distance & iad) const
{
  iad._distance = -1.0f;
  iad._item = -1;

  for (auto i = _fstart; i < _fstop; ++i)
  {
    if (_distance[i] < 0.0f)
      continue;

    if (_distance[i] > iad._distance)
    {
      iad._distance = _distance[i];
      iad._item = i;
    }
  }

  return;
}

//#define DEBUG_CLIENT

void
Fingerprint_Pool::wait_loop ()
{
  MPI::Group world_group_id = MPI::COMM_WORLD.Get_group ( );

  const auto p = MPI::COMM_WORLD.Get_size ( );

  int selected = 0;

  int n = (_fstop - _fstart);

  Item_and_Distance iad;   // scope here for efficiency

  while (1)
  {
    if (selected < n)    // we still have things to select
      _identify_unselected_item_with_longest_distance(iad);
    else
      iad._item = -1;

    MPI::COMM_WORLD.Send(reinterpret_cast<const void *>(&iad), sizeof(Item_and_Distance), MPI::UNSIGNED_CHAR, p-1, 0);   // send message to thread p, tell him we are done
#ifdef DEBUG_CLIENT
    cerr << "Client " << r << " selected " << selected << " posting " << iad._item << " at distance " << iad._distance << endl;
#endif

    int s;
    MPI::COMM_WORLD.Bcast(reinterpret_cast<void *>(&s), 1, MPI::INT, p-1);

    if (numeric_limits<int>::max() == s)
      break;

    if (iad._item < 0)    // no need to update anything more, we are all selected
      continue;

    if (s == iad._item)   // we owned the newly selected fingerprint
    {
      _distance[s] = -1.0f;
      selected++;
    }

    for (auto i = _fstart; i < _fstop; ++i)
    {
      if (_distance[i] < 0.0f)
        continue;

      float d = 1.0f - _pool[i].tanimoto(_pool[s]);
      if (d < _distance[i])
      {
        _distance[i] = d;
//      _pool[i].set_nbr(s);
      }
    }

#ifdef DEBUG_CLIENT
    cerr << "ID of selected fingerprint is " << s << endl;
#endif
  }

  return;
}

/*
  Will be called by controlling process
*/

static int
choose_next_item (Item_and_Distance * spr)
{
  auto p = MPI::COMM_WORLD.Get_size ( );

  auto n = p - 1;    // number of files

  float max_distance = -1.0f;
  int file_with_max_distance = -1;

  for (auto i = 0; i < n; ++i)
  {
    MPI::Status status;
    MPI::COMM_WORLD.Recv(reinterpret_cast<unsigned char *>(spr + i), sizeof(Item_and_Distance), MPI::UNSIGNED_CHAR, i, MPI::ANY_TAG, status);    // wait for the ID of the next selected to be known

//  cerr << "From worker " << i << " got distance " << spr[i]._distance << " item " << spr[i]._item << endl;

    if (spr[i]._item < 0)
      continue;

    if (spr[i]._distance > max_distance)
    {
      max_distance = spr[i]._distance;
      file_with_max_distance = i;
    }
  }

  return file_with_max_distance;
}

static int
gfp_spread_mpi ( int argc, char *argv[] )
{
  int rank = MPI::COMM_WORLD.Get_rank ( );

  Command_Line cl (argc, argv, "vn:L:r:s:q:");

  verbose = cl.option_count('v');

  if (cl.unrecognised_options_encountered())
  {
    if (0 == rank)
      cerr << "Unrecognised options encountered\n";
    usage(2);
  }

  int to_select = numeric_limits<int>::max();

  if (cl.option_present('n'))
  {
    if (! cl.value('n', to_select) || to_select < 1)
    {
      cerr << "The number of items to select (-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose && 0 == rank)
      cerr << "Will select " << to_select << " items\n";
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, only processes one file at a time\n";
    return 2;
  }

  IWString logfile_stem;

  if (cl.option_present('L'))
  {
    cl.value('L', logfile_stem);
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise progress report object\n";
      return 2;
    }
  }

  if (cl.option_present('q'))
  {
    if (! cl.value('q', squeeze_selected) || squeeze_selected < 1)
    {
      cerr << "The squeeze out seleted items (-q) option must be a whole +ve number\n";
      usage(2);
    }

    if (verbose && 0 == rank)
      cerr << "Each process will squeeze out already selected fingerprints every " << squeeze_selected << " items selected\n";
  }

  if (! cl.option_present('s'))
  {
    if (0 == rank)
      cerr << "Must specify fingerprints in file via the -s option\n";

    usage(1);
  }

  int pool_size = 0;

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      cerr << "The number of fingerprints in the input file (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose && 0 == rank)
      cerr << "Problem sized for " << pool_size << " fingerprints\n";
  }

//  Get the number of processes.

  auto p = MPI::COMM_WORLD.Get_size ( );

  if (p < 3)
  {
    if (0 == rank)
      cerr << "Must have at least three processes (two to process the file and one to control)\n";
    usage(2);
  }

  int items_per_process = pool_size / (p-1);

  if (0 == items_per_process)
  {
    if (0 == rank)
      cerr << "Processing " << pool_size << " items across " << p << " processes does not work\n";
    return 1;
  }

  items_per_process++;

  if (verbose && 0 == rank)
    cerr << "Processing " << pool_size << " items with " << p << " processes means each one gets " << items_per_process << " to process\n";

  Fingerprint_Pool fp;

  if (0 != rank)
    set_report_fingerprint_status(0);

  if (rank < p-1)    // one of the threads processing a file
  {
    int a = items_per_process * rank;
    int o = a + items_per_process;
    if (o > pool_size)
      o = pool_size;

    if (! fp.read_pool(cl[0], pool_size, a, o))
    {
      cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
      return 2;
    }

    if (logfile_stem.length() > 0)
      fp.initialise_log_file(logfile_stem);
  }

  if (rank < p-1)                // workers start waiting
  {
    fp.wait_loop();
    return 0;
  }

  assert (rank == p-1);

  int items_selected = 0;

  Item_and_Distance * spr = new Item_and_Distance[p-1];

#ifdef DEBUG_CLIENT
  cerr << "Controlling process starts loop, looking to select " << to_select << " items\n";
#endif

  while (1)
  {
    int j = choose_next_item(spr);

#ifdef DEBUG_CLIENT
    cerr << items_selected << " best selection from rank " << j << endl;
#endif

    if (j < 0)    // all threads fully selected, we are done
      break;

    int s = spr[j]._item;
    assert (s >= 0);

    MPI::COMM_WORLD.Bcast(reinterpret_cast<void *>(&s), 1, MPI::INT, p-1);

    cout << spr[j]._item << ' ' << spr[j]._distance << '\n';

    items_selected++;

    report_progress.report("Selected ", "\n", cerr);
    
    if (items_selected >= to_select)
      break;
  }

  cout.flush();

  if (verbose)
    cerr << "Selected " << items_selected << " items\n";

// We need to tell the workers to exit

  int notused = numeric_limits<int>::max();
  MPI::COMM_WORLD.Bcast(reinterpret_cast<void *>(&notused), 1, MPI::INT, p-1);
#ifdef NO_BROADCAST
  for (auto i = 0; i < p-1; ++i)   // change to Bcast
  {
    MPI::COMM_WORLD.Send(reinterpret_cast<const void *>(&notused), 1, MPI::INT, i, 0);
  }
#endif

  return 0;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
int
main ( int argc, char *argv[] )
{
  MPI::Init ( argc, argv);

  int rc = gfp_spread_mpi (argc, argv);

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize ( );

  return rc;
}
