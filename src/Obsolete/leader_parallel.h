#ifndef GFP_LEADER_PARA_H
#define GFP_LEADER_PARA_H

#include <iostream>
#include <memory>
#include <random>

#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/set_or_unset.h"

#include "Utilities/GFP_Tools/gfp.h"

/*
  Variant on the leader class for parallel applications
*/

typedef float score_t;

class GFP_PL : public IW_General_Fingerprint
{
  private:
    similarity_type_t _distance;

    int _selected;

    Set_or_Unset<similarity_type_t> _threshold;

//  Each item can limit the number of items in its cluster

    Set_or_Unset<int> _max_cluster_size;

    score_t _score;

//  shortest distance to a previously selected cluster centre

    similarity_type_t _shortest_distance_to_cluster_centre;

//  When we sort the selected list by cluster number, we need
//  that sort to be stable wrt initial input

    int _ndx;

    IWString _smiles;

  public:
    GFP_PL ();

    similarity_type_t   distance () const { return _distance;}
    void set_distance (similarity_type_t d) { _distance = d;}

    int selected () const { return _selected;}
    int & selected () { return _selected;}
    void set_selected (int s) { _selected = s;}
//  void set_selected (similarity_type_t s) { _selected = 1; _distance = s;}

    score_t score () const { return _score;}

    similarity_type_t shortest_distance_to_cluster_centre () const { return _shortest_distance_to_cluster_centre;}
    void set_shortest_distance_to_cluster_centre (similarity_type_t d) { _shortest_distance_to_cluster_centre = d;}

    int ndx () const { return _ndx;}
    void set_ndx (int s) { _ndx = s;}

    int construct_from_tdt (IW_TDT &, int &);

    int threshold (similarity_type_t & t) const { return _threshold.value (t);}
    int set_threshold (similarity_type_t t) { return _threshold.set (t);}

    int max_cluster_size (int & m) const { return _max_cluster_size.value (m);}

    const IWString & smiles () const { return _smiles;}
};

/*
  Struct needed by pthreads
*/

template <typename T>
struct Info_For_Section_Reader
{
  IWString _fname;
  off_t _start;
  off_t _stop;
  resizable_array<T *> _pool;
};

template <typename T>
class Pool_Formation_Info
{
  private:
    const int _verbose;

    int _number_single_file_parallel_readers;

    int _pool_size;

    T ** _pool;

//  private functions

    int _estimate_bytes_per_fingerprint (const char * fname) const;
    int _bytes_in_next_fingerprint (iwstring_data_source & input, off_t o) const;
    off_t _move_past_next_tdt (iwstring_data_source & input, off_t start) const;
    int _place_offsets_at_start_of_tdts (const char * fname,
                                struct Info_For_Section_Reader<T> * ifsr) const;

    int _build_pool_multiple_files (const Command_Line & cl);
    int _build_pool_single_file (const char * fname);
    int _build_pool_single_file (iwstring_data_source & input);
    int _do_parallel_build_pool_single_file (const char * fname);

  public:
    Pool_Formation_Info (int);

    int verbose () const { return _verbose;}

    void set_number_single_file_parallel_readers(int s) { _number_single_file_parallel_readers = s;}
    int number_single_file_parallel_readers () const { return _number_single_file_parallel_readers;}

    int allocate_pool (const char *, int);

    T ** pool () const { return _pool;}
    int pool_size () const { return _pool_size;}

    int build (const Command_Line & cl);
};

#ifdef LEADER_PARALLEL_IMPLEMENTATION

template <typename T>
Pool_Formation_Info<T>::Pool_Formation_Info(int s) : _verbose(s)
{
  _number_single_file_parallel_readers = 0;
  _pool_size = 0;
  _pool = NULL;

  return;
}

template <typename T>
int
Pool_Formation_Info<T>::allocate_pool (const char * caller,
                                       int s)
{
  assert (s > 0);
  assert (NULL == _pool);

  std::cerr << "Allocating pool of size " << s << '\n';
  _pool = new T * [s + 1];

  if (NULL == _pool)
  {
    std::cerr << "Cannot allocate " << (s + 1) << " pointers for pool\n";
    return 0;
  }

  _pool[s] = NULL;
  _pool_size = s;

  return 1;
}

/*
  A couple of structs needed by the pthreads
*/

template <typename T>
struct Data_For_File_Reader_Thread
{
  IWString _fname;

  resizable_array<T *> _pool;

  int _thread_number;
};

template <typename T>
void *
whole_file_reader_thread (void * v)
{
  Data_For_File_Reader_Thread<T> & dffrt = *(reinterpret_cast<Data_For_File_Reader_Thread<T> *>(v));

  // resizable_array<T *> & fp = dffrt._pool;

  iwstring_data_source input(dffrt._fname.null_terminated_chars());

  if (! input.good())
  {
    std::cerr << "Cannot open '" << dffrt._fname << "'\n";
    pthread_exit(NULL);
  }
  
  IW_TDT tdt;

  while (tdt.next(input))
  {
    T * tmp = new T;

    int fatal;
    if (! tmp->construct_from_tdt(tdt, fatal))
    {
      delete tmp;
      dffrt._pool.resize(0);
      pthread_exit(NULL);
    }

    dffrt._pool.add(tmp);
  }

  pthread_exit(NULL);

  return NULL;   // never comes here
}

template <typename T>
int
Pool_Formation_Info<T>::_build_pool_multiple_files (const Command_Line & cl)
{
  int n = cl.number_elements();
  assert (n > 1);

  pthread_t * child_thread = new pthread_t[n];

  Data_For_File_Reader_Thread<T> * dffrt = new Data_For_File_Reader_Thread<T>[n];

  pthread_attr_t attr;

  pthread_attr_init (&attr);
  pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < n; i++)
  {
    Data_For_File_Reader_Thread<T> & dffrti = dffrt[i];

    dffrti._fname = cl[i];
    dffrti._thread_number = i;

    if (0 != pthread_create (&(child_thread[i]), &attr, whole_file_reader_thread<T>, &dffrti))
    {
      perror ("cannot create child process");
      return 0;
    }
  }

  int s = 0;

  int rc = 1;

  for (int i = 0; i < n; i++)
  {
    pthread_join(child_thread[i], NULL);
    if (0 == dffrt[i]._pool.number_elements())
    {
      std::cerr << "Zero fingerprints read from '" << dffrt[i]._fname << "'\n";
      rc = 0;
    }
    else
      s += dffrt[i]._pool.number_elements();
  }

  if (0 == rc)    // one or more threads returned zero
    return 0;

  if (! allocate_pool("MULTIPLE", s))
    return 0;

  int ndx = 0;
  for (int i = 0; i < n; i++)
  {
    const T * const * f = dffrt[i]._pool.rawdata();
    int m = dffrt[i]._pool.number_elements();

    for (int j = 0; j < m; j++)
    {
      _pool[ndx + j] = const_cast<T *>(f[j]);
      _pool[ndx + j]->set_ndx(ndx + j);
    }

    dffrt[i]._pool.resize(0);

    ndx += m;
  }

  delete [] child_thread;
  delete [] dffrt;     // not clear that this deletes the object in the structs

  return 1;
}

template <typename T>
off_t
Pool_Formation_Info<T>::_move_past_next_tdt (iwstring_data_source & input,
                                             off_t start_pos) const
{
  if (! input.seekg(start_pos))
    return 0;

  IW_TDT tdt;
  if (! tdt.next(input))
    return 0;

  return input.tellg();
}

template <typename T>
int
Pool_Formation_Info<T>::_place_offsets_at_start_of_tdts (const char * fname,
                                struct Info_For_Section_Reader<T> * ifsr) const
{
  iwstring_data_source input(fname);
  if (! input.good())
    return 0;

  size_t file_size = dash_s(fname);

  off_t chunk_size = file_size / _number_single_file_parallel_readers;

  assert(chunk_size > 10);   // an absurdly small number

  for (int i = 0; i < _number_single_file_parallel_readers; i++)
  {
    Info_For_Section_Reader<T> & ifsri = ifsr[i];
    if (0 == i)
      ifsri._start = 0;
    else
      ifsri._start = ifsr[i - 1]._stop;

    if (i == (_number_single_file_parallel_readers - 1))
      ifsri._stop = std::numeric_limits<off_t>::max();
    else
      ifsri._stop = _move_past_next_tdt(input, ifsri._start + chunk_size);
  }

  return 1;
}

template <typename T>
void *
partial_file_reader_thread (void * v)
{
  struct Info_For_Section_Reader<T> * ifsr = reinterpret_cast<struct Info_For_Section_Reader<T> *>(v);

  iwstring_data_source input(ifsr->_fname.null_terminated_chars());

  if (! input.good())
  {
    std::cerr << "Cannot open '" << ifsr->_fname << "'\n";
    pthread_exit(NULL);
  }

  input.seekg(ifsr->_start);

  if (ifsr->_start > 0)
    input.skip_past("^\\|");

  IW_TDT tdt;
  while (tdt.next(input))
  {
    T * tmp = new T;

    int fatal;
    if (! tmp->construct_from_tdt(tdt, fatal))
    {
      delete tmp;
      ifsr->_pool.resize(0);
      pthread_exit(NULL);
    }

    ifsr->_pool.add(tmp);
    if (input.tellg() >= ifsr->_stop)
      break;
  }

  pthread_exit(NULL);

  return NULL;
}

template <typename T>
int
Pool_Formation_Info<T>::_bytes_in_next_fingerprint (iwstring_data_source & input,
                           off_t o) const
{
  if (! input.seekg(o))
    return 0;

  if (! input.skip_past("|"))
    return 0;

  off_t starting_position = input.tellg();

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if ('|' == buffer)
      return o - starting_position;

    o = input.tellg();
  }

  return 0;   // should not come ehre
}

template<typename T>
int
Pool_Formation_Info<T>::_estimate_bytes_per_fingerprint (const char * fname) const
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    std::cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  off_t file_size = input.file_size();

  Accumulator_Int<int> bytes;

  std::random_device rd;
  std::default_random_engine rng;
  std::uniform_int_distribution<off_t> u(0, file_size - 1);
  while (bytes.n() < 100)
  {
    const off_t o = u(rng);
    int b = _bytes_in_next_fingerprint (input, o);
    if (b > 0)
      bytes.extra(b);
  }

  return static_cast<int>(bytes.average());
}

template <typename T>
int
Pool_Formation_Info<T>::_do_parallel_build_pool_single_file (const char * fname)
{
  assert (NULL == _pool);

  size_t file_size = dash_s(fname);

  if (0 == file_size)
  {
    std::cerr << "Empty file '" << fname << "'\n";
    return 0;
  }

  int bytes_per_fingerprint = _estimate_bytes_per_fingerprint(fname);

  if (bytes_per_fingerprint < 50)
  {
    std::cerr << "Impossibly small bytes per fingerprint " << bytes_per_fingerprint << '\n';
    return 0;
  }

  // int estimate_fingerprints_present = file_size / bytes_per_fingerprint + 1000000;  // couple extra MB

  off_t bytes_per_chunk = file_size / _number_single_file_parallel_readers;

  pthread_t * child_thread = new pthread_t[_number_single_file_parallel_readers]; std::unique_ptr<pthread_t[]> free_child_thread(child_thread);

  Info_For_Section_Reader<T> * ifsr = new Info_For_Section_Reader<T>[_number_single_file_parallel_readers]; std::unique_ptr<Info_For_Section_Reader<T>[]> free_ifsr(ifsr);

// Tried to put this block in a separate method, but could not get the
// template stuff to compile

  _place_offsets_at_start_of_tdts(fname, ifsr);

  pthread_attr_t attr;

  pthread_attr_init (&attr);
  pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < _number_single_file_parallel_readers; i++)
  {
    Info_For_Section_Reader<T> & ifsri = ifsr[i];

    ifsri._fname = fname;
    ifsri._start = i * bytes_per_chunk;
    ifsri._stop = ifsri._start + bytes_per_chunk;

//  std::cerr << "Section " << i << " start at " << ifsri._start << " end at " << ifsri._stop << '\n';

    if (0 != pthread_create (&(child_thread[i]), &attr, partial_file_reader_thread<T>, &ifsri))
    {
      perror ("cannot create child process");
      return 0;
    }
  }

  int s = 0;

  int rc = 1;
  for (int i = 0; i < _number_single_file_parallel_readers; i++)
  {
    pthread_join(child_thread[i], NULL);
    if (0 == ifsr[i]._pool.number_elements())
    {
      std::cerr << "No fingerprints read by file chunk reader " << i << '\n';
      rc = 0;
    }
    else 
      s += ifsr[i]._pool.number_elements();

    std::cerr << "Pool reader " << i << " read " << ifsr[i]._pool.number_elements() << " items\n";
  }

  if (0 == rc)
    return 0;

  if (! allocate_pool(fname, s))
    return 0;

  int ndx = 0;
  for (int i = 0; i < _number_single_file_parallel_readers; i++)
  {
    Info_For_Section_Reader<T> & ifsri = ifsr[i];

    int number_fingerprints = ifsri._pool.number_elements();

    const T * const * rawdata = ifsri._pool.rawdata();

    for (int k = 0; k < number_fingerprints; k++)
    {
      _pool[ndx + k] = const_cast<T *>(rawdata[k]);
      _pool[ndx + k]->set_ndx(ndx + k);
    }

    ifsri._pool.resize(0);

    ndx += number_fingerprints;
  }

  return 1;
}

template <typename T>
int
Pool_Formation_Info<T>::_build_pool_single_file (iwstring_data_source & input)
{
  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    tdts_read++;

    T * f = new T;

    int fatal;
    if (! f->construct_from_tdt (tdt, fatal))
    {
      delete f;
      if (fatal)
      {
        std::cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    f->set_ndx(items_in_pool);
    _pool[items_in_pool] = f;

    items_in_pool++;

    if (items_in_pool == _pool_size)
    {
      std::cerr << "Pool is full, max " << _pool_size << '\n';
      break;
    }
  }

  _pool_size = items_in_pool;

  return 1;
}

template <typename T>
int
Pool_Formation_Info<T>::_build_pool_single_file (const char * fname)
{
  size_t file_size = dash_s(fname);

  if (file_size < 1000)   // too small to do parallel
    ;
  else if (_number_single_file_parallel_readers > 0)
    return _do_parallel_build_pool_single_file (fname);

  iwstring_data_source input(fname);
  if (! input.good())
  {
    std::cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == _pool_size)
  {
    IWString tmp;
    tmp << "^\\|";

    std::unique_ptr<re2::RE2> pcn = std::make_unique<re2::RE2>(tmp);
    int s = input.grep (*pcn);

    if (0 == s)
    {
      std::cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    if (! allocate_pool(fname, s))
      return 0;
  }

  return _build_pool_single_file (input);
}

template <typename T>
int
Pool_Formation_Info<T>::build (const Command_Line & cl)
{
  assert (NULL == _pool);

  int rc;
  if (1 == cl.number_elements())
    rc = _build_pool_single_file(cl[0]);
  else
    rc = _build_pool_multiple_files (cl);

  if (0 == rc)
    return 0;

  return rc;
}

#endif
#endif
