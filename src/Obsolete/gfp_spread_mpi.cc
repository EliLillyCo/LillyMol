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

static int precompute_distances_while_waiting = 0;

static void
usage (int rc,
       int rank)
{
  if (rank > 0)
    exit(rc);

  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "MPI implementation of spread. Input must be a set of pre-split .gfp files with -s fp's in each\n";
  cerr << " -s <items>     number of fingerprints in EACH input file\n";
  cerr << " -n <items>     number of items to select\n";
  cerr << " -L <fname>     use a log file\n";
  cerr << " -r <n>         report progress every <n> items selected\n";
  cerr << " -p <n>         precompute <n> distances while waiting for communication\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Spread_Item : public GFP_Standard
{
  private:
    int _selected;      // have we been selected or not
    float _distance;    // distance to nearest previously selected item
    int _nbr;           // previously selected nearest nbr

    int _id;                // this will be a globally unique number

  public:
    Spread_Item();
    ~Spread_Item();

    int build (IW_TDT & tdt);

    int debug_print (ostream &) const;

    float distance () const { return _distance;}
    void set_distance(float s) { _distance = s;}
    void set_distance(float s, int n) { _distance = s; _nbr = n;}

    int id () const { return _id;}
    int nbr () const { return _nbr;}

    int selected () const { return _selected;}
    void set_selected (int s) { _selected = s;}

    void set_id     (unsigned int s) { _id = s;}
    void set_nbr    (int n) { _nbr = n;}

    int update_distance_if_needed (const Spread_Item & sel);
};

Spread_Item::Spread_Item()
{
  _selected = 0;
  _distance = 2.0f;
  _nbr = -1;

  _id = -1;

  return;
}

Spread_Item::~Spread_Item()
{
  return;
}

int
Spread_Item::debug_print (ostream & output) const
{
  output << "Spread item " << _id << endl;
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

int
Spread_Item::update_distance_if_needed (const Spread_Item & rhs)
{
  if (_selected)
    return 0;

  float d = 1.0f - tanimoto(rhs);

  if (d < _distance)
  {
    _distance = d;
    _nbr = rhs._id;
  }

  return 1;
}

class Fingerprint_Pool
{
  private:
    Spread_Item * _pool;
    int _n;

    float * _distance;
    float * _new_distance;

    ofstream _logfile;

//  private functions

    int _assign_sequential_identifiers (unsigned int id);
    void _identify_unselected_item_with_longest_distance (float & max_distance, int & item_with_max_distance) const;
    void _precompute_new_distances_while_waiting (int);
    int _update_fingerprints_based_on_newly_selected (const Spread_Item &);

  public:
    Fingerprint_Pool();
    ~Fingerprint_Pool();

    int read_pool (const char * fname, int);
    int read_pool (iwstring_data_source &, int);

    int assign_sequential_identifiers ();
    int assign_sequential_identifiers (unsigned int);

    void wait_loop ();

    int initialise_log_file (const IWString &);
};

Fingerprint_Pool::Fingerprint_Pool ()
{
  _pool = nullptr;
  _n = 0;
  _distance = nullptr;
  _new_distance = nullptr;

  return;
}

Fingerprint_Pool::~Fingerprint_Pool ()
{
  if (NULL != _pool)
    delete [] _pool;

  if (NULL != _distance)
    delete [] _distance;

  if (NULL != _new_distance)
    delete [] _new_distance;

  return;
}

/*
  Reading the pool is a little tricky.
  We need to make sure that the whole file gets read - across
  the different processes.
*/

int
Fingerprint_Pool::read_pool (iwstring_data_source & input,
                             int pool_size)
{
  _n = pool_size;

  _pool = new Spread_Item[_n];
  _distance = new_float(_n, 2.0f);

  if (precompute_distances_while_waiting)
  {
    if (precompute_distances_while_waiting >= _n)
    {
      cerr << "Fingerprint_Pool::read_pool:invalid precompute value " << precompute_distances_while_waiting << " for pool with " << _n << " fingerprints\n";
      return 0;
    }
    _new_distance = new float[_n];
  }

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
                             int pool_size)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_pool (input, pool_size);
}

int
Fingerprint_Pool::_assign_sequential_identifiers (unsigned int id)
{
#ifdef DEBUG_ASSIGN_SEQUENTIAL_IDENTIFIERS
  if (verbose)
    cerr << "Assigning sequential identifiers starting with " << id << endl;
#endif

  for (auto i = 0; i < _n; ++i)
  {
    _pool[i].set_id(id);
    id++;
  }

  auto p = MPI::COMM_WORLD.Get_size ( );
  auto r = MPI::COMM_WORLD.Get_rank ( );

#ifdef DEBUG_ASSIGN_SEQUENTIAL_IDENTIFIERS
  if (r < (p-2))
    cerr << "Rank " << r << " sending message to start numbering at " << id << endl;
#endif

  if (r < (p-2))
    MPI::COMM_WORLD.Send (&id, 1, MPI::UNSIGNED, r+1, 0);   // tell next chunk to where to start

  return 1;
}

int
Fingerprint_Pool::assign_sequential_identifiers (unsigned int id)    // thread 0
{
#ifdef DEBUG_ASSIGN_SEQUENTIAL_IDENTIFIERS
  cerr << "Unique identifiers being assigned starting at " << id << endl;
#endif

  return _assign_sequential_identifiers(id);
}

int
Fingerprint_Pool::assign_sequential_identifiers ()    // all other threads, need to wait for previous thread to know start value for index
{
  auto r = MPI::COMM_WORLD.Get_rank ( );
  assert (r > 0);

  unsigned int id;
  MPI::Status status;
  MPI::COMM_WORLD.Recv(&id, 1, MPI::UNSIGNED, r-1, 0 /*tag*/, status);    // wait for the ID of the next selected to be known

#ifdef DEBUG_ASSIGN_SEQUENTIAL_IDENTIFIERS
  cerr << "Process " << r << " identifiers will be assigned starting at " << id << endl;
#endif

  if (1 != status.Get_count(MPI::UNSIGNED))
  {
    cerr << "Assigning unique id, got " << status.Get_count(MPI::UNSIGNED) << " items, id " << id << endl;
    return 0;
  }

  return _assign_sequential_identifiers(id);
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
Fingerprint_Pool::_identify_unselected_item_with_longest_distance (float & max_distance, 
                                                                   int & item_with_max_distance) const
{
  max_distance = -1.0f;
  item_with_max_distance = -1;

  for (auto i = 0; i < _n; ++i)
  {
    if (_distance[i] > max_distance)
    {
      max_distance = _distance[i];
      item_with_max_distance = i;
    }
  }

  return;
}

void
Fingerprint_Pool::_precompute_new_distances_while_waiting(int s)
{
  int pc = 0;
  for (auto i = 0; pc < precompute_distances_while_waiting; ++i)
  {
    _new_distance[i] = -1.0f;

    if (_pool[i].selected())
      continue;

    if (i == s)
      continue;

    _new_distance[i] = 1.0f - _pool[i].tanimoto(_pool[s]);
    pc++;
  }

  _new_distance[pc] = -2.0f;   // will fail if precompute_distances_while_waiting is same as _n

  return;
}

void
Fingerprint_Pool::wait_loop ()
{
  const auto p = MPI::COMM_WORLD.Get_size ( );
  const auto r = MPI::COMM_WORLD.Get_rank ( );

  int selected = 0;

  if (verbose > 1)
    cerr << "Task " << r << " entering loop between " << _pool[0].id() << " to " << _pool[_n-1].id() << " N = " << _n << endl;

  while (1)
  {
    float max_distance = 0.0f;
    int item_with_max_distance = -1;

    if (selected < _n)
      _identify_unselected_item_with_longest_distance(max_distance, item_with_max_distance);

    if (item_with_max_distance < 0)    // just send the first one, mark it selected, will be ignored by master process
    {
      _pool[0].set_selected(1);
      MPI::COMM_WORLD.Send(reinterpret_cast<const void *>(_pool),                          sizeof(Spread_Item), MPI::UNSIGNED_CHAR, p-1, 0);   // send message to thread p, tell him we are done
    }
    else
    {
      _pool[item_with_max_distance].set_distance(_distance[item_with_max_distance]);
      MPI::COMM_WORLD.Send(reinterpret_cast<const void *>(_pool + item_with_max_distance), sizeof(Spread_Item), MPI::UNSIGNED_CHAR, p-1, 0);   // send message to thread p, tell him we are done
    }

    if (item_with_max_distance >= 0 && NULL != _new_distance)
      _precompute_new_distances_while_waiting(item_with_max_distance);

    Spread_Item s;
    MPI::COMM_WORLD.Bcast(reinterpret_cast<void *>(&s), sizeof(Spread_Item), MPI::UNSIGNED_CHAR, p-1);    // wait for the ID of the next selected to be known

    if (numeric_limits<int>::max() == s.id())   // signal from the caller that we are done
      break;

    if (item_with_max_distance < 0)    // no need to update anything more, we are all selected
      continue;

    selected += _update_fingerprints_based_on_newly_selected (s);
  }

  return;
}

/*
  We return 1 or 0 depending on whether or not the newly selected item is part of our set
*/

int
Fingerprint_Pool::_update_fingerprints_based_on_newly_selected (const Spread_Item & sel)
{
  auto id = sel.id();

  int rc;

  if (id >= _pool[0].id() && id <= _pool[_n-1].id())
  {
    int ndx = id - _pool[0].id();
    _distance[ndx] = -1.0f;
    rc = 1;
  }
  else
    rc = 0;

  int istart = 0;

  if (1 == rc && NULL != _new_distance)    // the new fingerprint was ours, and we have precomputed distances, use them
  {
    int pc = 0;
    for (auto i = 0; pc < precompute_distances_while_waiting; ++i)
    {
      if (_distance[i] < 0.0f)
        continue;

      if (_pool[i].distance() > _new_distance[i])
        _pool[i].set_distance(_new_distance[i], id);
    }
    istart = precompute_distances_while_waiting;
  }

  for (auto i = istart; i < _n; ++i)
  {
    if (_distance[i] < 0.0f)
      continue;

    float d = 1.0f - _pool[i].tanimoto(sel);
    if (d < _distance[i])
    {
      _distance[i] = d;
      _pool[i].set_distance(d, id); 
    }
  }

  return rc;
}

/*
  Will be called by controlling process
*/

static int
choose_next_item (Spread_Item * spr,
                  int tag)
{
  auto p = MPI::COMM_WORLD.Get_size ( );

  auto n = p - 1;    // number of files

  float max_distance = -1.0f;
  int item_with_max_distance = -1;

  for (auto i = 0; i < n; ++i)
  {
    Spread_Item s;
    MPI::Status status;
    MPI::COMM_WORLD.Recv(reinterpret_cast<unsigned char *>(spr + i), sizeof(Spread_Item), MPI::UNSIGNED_CHAR, i, MPI::ANY_TAG, status);    // wait for the ID of the next selected to be known

#ifdef DEBUG_CLIENT
    IWString msg;
    msg << "From task " << i << " selected flag " << spr[i].selected() << "\n";
    cerr << msg;
#endif

    if (spr[i].selected())
      continue;

    if (spr[i].distance() > max_distance)
    {
      max_distance = spr[i].distance();
      item_with_max_distance = i;
    }
  }

  return item_with_max_distance;
}

//****************************************************************************80

static int
gfp_spread_mpi ( int argc, char *argv[] )
{
  const auto p = MPI::COMM_WORLD.Get_size ( );

  const auto rank = MPI::COMM_WORLD.Get_rank ( );

  Command_Line cl (argc, argv, "vn:L:r:s:p:");

  verbose = cl.option_count('v');

  int n = cl.number_elements();

  if (n < 2)
  {
    if (0 == rank)
      cerr << "Must specify > 1 pre-split gfp files on the command line\n";
    usage(2, rank);
  }

  if (! cl.option_present('s'))
  {
    if (0 == rank)
      cerr << "Must specify fingerprints in each input file via the -s option\n";
    usage(2, rank);
  }

  if (cl.unrecognised_options_encountered())
  {
    if (0 == rank)
      cerr << "Unrecognised options encountered\n";
    usage(2, rank);
  }

  int to_select = numeric_limits<int>::max();

  if (cl.option_present('n'))
  {
    if (! cl.value('n', to_select) || to_select < 1)
    {
      if (0 == rank)
        cerr << "The number of items to select (-n) must be a whole +ve number\n";
      usage(2, rank);
    }

    if (verbose && 0 == rank)
      cerr << "Will select " << to_select << " items\n";
  }

  IWString logfile_stem;

  if (cl.option_present('L'))
  {
    cl.value('L', logfile_stem);
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', p-1 == rank ? verbose : 0))
    {
      if (0 == rank)
        cerr << "Cannot initialise progress report object\n";
      return 2;
    }
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', precompute_distances_while_waiting) || precompute_distances_while_waiting < 1)
    {
      if (0 == rank)
        cerr << "The precompute distances while waiting option (-p) must be a whole +ve number\n";
    }

    if (verbose && 0 == rank)
      cerr << "Will precompute " << precompute_distances_while_waiting << " fingerprints while waiting\n";
  }

  int pool_size = 0;
  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      if (0 == rank)
        cerr << "The pool size (-s) must be a whole +ve number\n";
      usage(2, rank);
    }
  }

  if (p != n+1)
  {
    if (0 == rank)
      cerr << "The number of processes must be 1 more than the number of files\n";
    return 1;
  }

// by convention, we have a task for each file, and a task for managing the other processes

  Fingerprint_Pool fp;

  if (0 != rank)
    set_report_fingerprint_status(0);

  if (rank < n)    // one of the threads processing a file
  {
    if (! fp.read_pool(cl[rank], pool_size))
    {
      cerr << "Cannot read fingerprints from '" << cl[rank] << "'\n";
      return 2;
    }

    if (logfile_stem.length() > 0)
      fp.initialise_log_file(logfile_stem);
  }

// Every thread has read it's input file

  if (0 == rank)
    fp.assign_sequential_identifiers(0);
  else if (rank < n)
    fp.assign_sequential_identifiers();     // will wait for earlier ones...

#ifdef DEBUG_ASSIGN_SEQUENTIAL_IDENTIFIERS
  cerr << "Identifiers assigned\n";
#endif

  if (rank < p-1)                // workers start waiting
  {
    fp.wait_loop();
    return 0;
  }

  assert (rank == p-1);

  int items_selected = 0;

  Spread_Item * spr = new Spread_Item[n];

#ifdef DEBUG_CLIENT
  cerr << "Controlling process starts loop, looking to select " << to_select << " items\n";
#endif

  while (1)
  {
    int j = choose_next_item(spr, items_selected+1);

#ifdef DEBUG_CLIENT
    IWString msg;
    msg << items_selected << " best selection from rank " << j << "\n";
    cerr << msg;
#endif

    if (j < 0)    // all threads fully selected, we are done
      break;

    MPI::COMM_WORLD.Bcast(reinterpret_cast<void *>(spr + j), sizeof(Spread_Item), MPI::UNSIGNED_CHAR, p-1);

    const Spread_Item & sel = spr[j];
    cout << sel.id() << ' ' << sel.distance() << ' ' << sel.nbr() << '\n';

    items_selected++;

    report_progress.report("Selected ", "\n", cerr);
    
    if (items_selected >= to_select)
      break;
  }

  cout.flush();

  if (verbose)
    cerr << "Selected " << items_selected << " items\n";

// We need to tell the workers to exit

  Spread_Item notused;
  notused.set_id(numeric_limits<int>::max());
  MPI::COMM_WORLD.Bcast(reinterpret_cast<void *>(&notused), sizeof(Spread_Item), MPI::UNSIGNED_CHAR, p-1);

  return 0;
}
//****************************************************************************80

static void
timestamp ( )

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
