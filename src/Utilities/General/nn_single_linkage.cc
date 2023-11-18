/*
  Performs single linkage clustering where input is a list
  of neighbours
*/

#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int items_selected = 0;

static int clusters_formed = 0;

static extending_resizable_array<int> cluster_size;
static extending_resizable_array<int> iterations_per_cluster;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Performs single linkage clustering on a neighbour list\n";
  cerr << " -q             fast exit, do not wait for hash deallocation\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static IW_STL_Hash_Map_int id_to_ndx;
static resizable_array_p<IWString> ndx_to_id;
static int * already_selected = NULL;

class Set_of_Neighbours : public resizable_array<int>
{
  private:
    int _selected;

//  private functions

    void _shell_sort_nbrs();

  public:
    Set_of_Neighbours ();

    int build (const const_IWSubstring &);

    int selected() const { return _selected;}
    void set_selected (int s) { _selected = s;}

    int add_neighbours (resizable_array<int> & n, int * already_selected) const;

    int contains (int needle) const;
};

Set_of_Neighbours::Set_of_Neighbours ()
{
  _selected = 0;

  return;
}

void
Set_of_Neighbours::_shell_sort_nbrs()
{
  int i, j, increment;
  increment = _number_elements / 2;
  int save_i, temp;
 
  while (increment > 0)
  {
    for (i = increment; i < _number_elements; i++) 
    {
      j = i;
      temp = _things[i];
      save_i = _things[i];
      while ((j >= increment) && (_things[j-increment] > temp))
      {
        _things[j] = _things[j - increment];
        j = j - increment;
      }
      _things[j] = save_i;
    }
 
    if (increment == 2)
       increment = 1;
    else 
       increment = (int) (increment / 2.2);
  }
}

int
Set_of_Neighbours::build (const const_IWSubstring & buffer)
{
  int nw = buffer.nwords();
  assert (nw > 0);

  resizable_array<int>::resize(nw);

  int i = 0;
  IWString token;

  while (buffer.nextword(token, i))
  {
    IW_STL_Hash_Map_int::const_iterator f = id_to_ndx.find(token);

    if (f == id_to_ndx.end())
    {
      int s = id_to_ndx.size();
      id_to_ndx[token] = s;
      resizable_array<int>::add(s);
      IWString * tmp = new IWString(token);
      ndx_to_id.add(tmp);
    }
    else
      resizable_array<int>::add((*f).second);
  }

  _shell_sort_nbrs();

  return nw;
}

/*
  Take advantage of the fact that the list is sorted
*/

int
Set_of_Neighbours::contains (int needle) const
{
#ifdef DEBUG_BINARY_SEARCH
  cerr << "Begin binary search over " << _number_elements << " items for " << needle << endl;

  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " i = " << i << " alue " << _things[i] << endl;
  }
#endif

  if (needle < _things[0])
    return 0;

  int right = _number_elements - 1;

  if (needle > _things[right])
    return 0;

  int left = 0;
  int middle, nmid;   // scope here for efficiency

  while (right > left)
  {
    middle = (left + right) / 2;
    if (middle == left)
    {
//    cerr << "Middle == left, match" << (needle == _things[middle]) << endl;
      return (needle == _things[middle]);
    }

    nmid = _things[middle];

//  cerr << "Left = " << left << " value " << _things[left] << " middle " << middle << " value " << _things[middle] << " right " << right << " value " << _things[right] << endl;
    if (needle < nmid)
      right = middle;
    else if (needle > nmid)
      left = middle;
    else
    {
//    cerr << "Found match, i = " << middle << endl;
      return 1;
    }
  }

  return 0;
}

int
Set_of_Neighbours::add_neighbours (resizable_array<int> & n,
                                   int * already_selected) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    int j = _things[i];

    if (already_selected[j])
      continue;

    n.add(j);

    already_selected[j] = 1;
  }

  return 1;
}

static int
nn_single_linkage (resizable_array_p<Set_of_Neighbours> & pool,
                   Set_of_Neighbours * ldr,
                   IWString_and_File_Descriptor & output)
{
  static resizable_array<int> next_shell;

  int pool_size = pool.number_elements();

  if (0 == next_shell.number_elements())
    next_shell.resize(pool_size);
  else
    next_shell.resize_keep_storage(0);

  ldr->add_neighbours(next_shell, already_selected);

  int iterations_this_cluster = 0;
  int items_this_cluster = 0;

  while (next_shell.number_elements())
  {
    int n = next_shell.number_elements();

    for (int i = 0; i < n; i++)
    {
      int ni = next_shell[i];

      output << ' ' << (*(ndx_to_id[ni]));

      for (int j = 0; j < pool_size; j++)
      {
        Set_of_Neighbours * pj = pool[j];

        if (pj->selected())
          continue;

        if (! pj->contains(ni))
          continue;

        pj->add_neighbours(next_shell, already_selected);
        pj->set_selected(1);
      }
    }

    items_this_cluster += n;

    next_shell.erase (0, n - 1);   // remove the first N to leave only what we just added

    iterations_this_cluster++;

    output.write_if_buffer_holds_more_than(32768);
  }

  output << '\n';
  output.write_if_buffer_holds_more_than(32768);

  clusters_formed++;
  items_selected += items_this_cluster;

  if (verbose)
  {
    cluster_size[items_this_cluster]++;
    iterations_per_cluster[iterations_this_cluster]++;

    if (verbose > 1)
      cerr << "cluster " << clusters_formed << " contains " << items_this_cluster << " items, selected " << items_selected << endl;
  }

  return 1;
}

static int
nn_single_linkage (resizable_array_p<Set_of_Neighbours> & pool,
                   IWString_and_File_Descriptor & output)
{
  int n = pool.number_elements();

  for (int i = 0; i < n; i++)
  {
    Set_of_Neighbours * pi = pool[i];

    if (pi->selected())
      continue;

    nn_single_linkage(pool, pi, output);
  }

  return 1;
}

static int
build_pool (iwstring_data_source & input,
            resizable_array_p<Set_of_Neighbours> & pool)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    Set_of_Neighbours * s = new Set_of_Neighbours;

    if (! s->build(buffer))
    {
      delete s;
      cerr << "Fatal error reading line " << input.lines_read() << endl;
      return 0;
    }
    else
      pool.add(s);
  }

  return pool.number_elements();
}

static int
build_pool (const char * fname,
            resizable_array_p<Set_of_Neighbours> & pool)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return build_pool(input, pool);
}


static int
nn_single_linkage (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vq");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  resizable_array_p<Set_of_Neighbours> pool;
  pool.resize(100000);

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! build_pool(cl[i], pool))
    {
      cerr << "Cannot read nn items from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  if (verbose)
    cerr << "Read " << pool.number_elements() << " neighbour relations involving " << id_to_ndx.size() << " items\n";

  already_selected = new_int(id_to_ndx.size());

  IWString_and_File_Descriptor output(1);

  nn_single_linkage(pool, output);

  output.flush();

  if (verbose)
  {
    cerr << "Clustered " << id_to_ndx.size() << " items into " << clusters_formed << " clusters\n";
    for (int i = 0; i < cluster_size; i++)
    {
      if (cluster_size[i])
        cerr << cluster_size[i] << " clusters of size " << i << endl;
    }

    for (int i = 0; i < iterations_per_cluster.number_elements(); i++)
    {
      if (iterations_per_cluster[i])
        cerr << iterations_per_cluster[i] << " clusters took " << i << " expansion iterations\n";
    }
  }

  if (cl.option_present('q'))
  {
    cerr << "Quick exit\n";
    _exit(0);
  }

  if (verbose)
    cerr << "Processing complete, hash deallocation...\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = nn_single_linkage(argc, argv);

  return rc;
}
