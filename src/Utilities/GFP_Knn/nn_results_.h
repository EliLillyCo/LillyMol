#ifndef NN_RESULTS_IMPL
#define  NN_RESULTS_IMPL

#include <stdlib.h>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

#include "nn_results.h"
#include "extract_from_tdt_form.h"

template <typename N>
NN_Results_Base<N>::NN_Results_Base ()
{
  _results = nullptr;
  _number_items = 0;

  return;
}

template <typename N>
NN_Results_Base<N>::~NN_Results_Base ()
{
  if (NULL != _results)
    delete [] _results;

  return;
}

template <typename N>
int 
NN_Results_Base<N>::build (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return build (input);
}

template <typename N>
int
NN_Results_Base<N>::_resize ()
{
  if (NULL != _results)
    delete [] _results;

  assert (_number_items > 0);

  _results = new N[_number_items];

  if (NULL == _results)
  {
    cerr << "NN_Results_Base::_results: memory failure, " << _number_items << " items\n";
    return 0;
  }

  return _number_items;
}


template <typename N>
int
NN_Results_Base<N>::_count_items_in_input_and_resize (iwstring_data_source & input)
{
  std::unique_ptr<re2::RE2> rx = std::make_unique<re2::RE2>("^\\|");

  _number_items = input.grep(*rx);

  if (_number_items <= 0)
  {
    cerr << "NN_Results_Base::build: no tdt's in input file\n";
    return 0;
  }

  cerr << "Found " << _number_items << " items in input pool\n";

  if (! _resize ())
  {
    cerr << "NN_Results_Base::build: cannot resize for " << _number_items << " items\n";
    return 0;
  }

  return 1;
}

template <typename N>
int
NN_Results_Base<N>::build (iwstring_data_source & input)
{
  assert (NULL == _results);

  if (! _count_items_in_input_and_resize (input))
    return 0;

  int fatal = 0;
  for (int i = 0; i < _number_items; i++)
  {
    if (_results[i].build (input))
      continue;

//    cerr << "build function error" << endl;
//    return 0;      
/*GH    if (_results[i].build (input, fatal))
      continue;
*/
    if (! fatal)
      return 1;

    cerr << "NN_Results_Base::build: error building pool item " << i << endl;
    return 0;
  }

  return 1;
}

template <typename N>
int
NN_Item_Base<N>::debug_print (ostream & os, int verbose) const
{
  os << "Item '" << _id << "' has " << _number_elements << " neighbours\n";

  if (verbose)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      const N * sid = _things[i];
      os << "  '" << sid->id () << "' at distance " << sid->distance () << endl;
    }
  }

  return os.good ();
}

template <typename N> template <typename F>
void
NN_Item_Base<N>::change_distances (F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->change_distance (f);
  }

  return;
}

static IWString smiles_tag = "$SMI<";
static IWString identifier_tag = "PCN<";
static IWString neighbour_tag = "NBR<";
static IWString distance_tag = "DIST<";

template <typename N>
int
NN_Item_Base<N>::build (iwstring_data_source & input)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (buffer.starts_with (smiles_tag))
      break;

    if ('|' == buffer)     // no more neighbours
    {
      cerr << "NN_Item_Base::build: no neighbours?\n";
      return 0;
    }
  }

  if (0 == buffer.length ())
  {
    if (! input.eof ())
      cerr << "NN_Item_Base::build: very strange, blank line, but not EOF\n";
      
    return 0;
  }

  if (! extract_from_tdt_form (buffer, smiles_tag, _smiles))
  {
    cerr << "NN_Item_Base::build: invalid smiles form '" << buffer << "', line " << input.lines_read () << endl;
    return 0;
  }

// The fingerprints and other stuff may be nearby, search for the identifier record

  while (input.next_record (buffer))
  {
    if ('|' == buffer)
    {
      if (0 == _id.length ())
      {
        cerr << "NN_Item_Base::build: item with no identifier\n";
        return 0;
      }

      return 1;
    }

    if (! buffer.starts_with (identifier_tag))
      continue;

    if (! extract_from_tdt_form (buffer, identifier_tag, _id))
    {
      cerr << "NN_Item_Base::build: invalid identifier form '" << buffer << "', line " << input.lines_read () << endl;
      return 0;
    }

    break;
  }

// Now that we have the smiles and identifier of our molecule, let's go looking for our neighbours

  while (1)
  {
    N * tmp = new N;

    int fatal = 0;
    //cerr << "build function error" << endl;
    //return 0;
/*GH    if (! tmp->build (input, fatal)) */
    if (! tmp->build (input, fatal))
    {
      delete tmp;

      if (fatal)
      {
        cerr << "NN_Item_Base::build: fatal error on line " << input.lines_read () << endl;
        return 0;
      }

      return 1;
    }

    this->add (tmp);
  }

  return 1;
}

template <typename N>
int
NN_Item_Base<N>::build_from_neighbour_list (const pair<int, similarity_type_t> * p,
                                            int n)
{
  if (! resizable_array_p<N>::resize (n))
  {
    cerr << "NN_Results_Base::build: cannot resize for " << n << " neighbours\n";
    return 0;
  }

  for (int i = 0; i < n; i++)
  {
    const pair<int, similarity_type_t> & pi = p[i];

    N * nbr = new N;

    nbr->set_index_in_pool (pi.first);
    nbr->set_distance (pi.second);

    this->add (nbr);
  }

  return 1;
}

template <typename N>
int
NN_Results_Base<N>::debug_print (ostream & os, int verbose) const
{
  os << "Near neighbour results with " << _number_items << " items\n";

  int neighbours_stored = 0;

  for (int i = 0; i < _number_items; i++)
  {
    os << "ITEM " << i << endl;

    const N & ni = _results[i];
    ni.debug_print (os, verbose);

    neighbours_stored += ni.number_elements ();
  }

  os << neighbours_stored << " total neighbours stored\n";

  return os.good ();
}

template <typename N>
int
NN_Results_Base<N>::build_from_distance_matrix (const char * fname,
                                                const DM_to_NN_Conditions<float> & dmc,
                                                const IWString_STL_Hash_Set & only_use)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "NN_Results_Base::build_from_distance_matrix: cannot open '" << fname << "'\n";
    return 0;
  }

  return build_from_distance_matrix (input, dmc, only_use);
}

template <typename N>
int
NN_Results_Base<N>::build_from_distance_matrix (iwstring_data_source & input,
                                                const DM_to_NN_Conditions<float> & dmc,
                                                const IWString_STL_Hash_Set & only_use)
{
  IWDistanceMatrixFloat dm;

  if (! dm.do_read (input))
  {
    cerr << "NN_Results_Base::build_from_distance_matrix: cannot read distance matrix file\n";
    return 0;
  }

  _number_items = dm.number_molecules ();

  if (0 == _number_items)
  {
    cerr << "NN_Results_Base::build_from_distance_matrix: distance matrix contains no molecules\n";
    return 0;
  }

  if (! _resize ())
    return 0;

// If we are doing a subset, we need to compute the index for the used items in subset

  int * ndx = new int[_number_items]; std::unique_ptr<int[]> free_ndx (ndx);

  if (0 == only_use.size ())
  {
    for (int i = 0; i < _number_items; i++)
    {
      ndx[i] = i;
    }
  }
  else
  {
    set_vector (ndx, _number_items, -1);

    int n = 0;
    for (int i = 0; i < _number_items; i++)
    {
      const IWString & s = dm.id (i);

      if (! only_use.contains (s))
        continue;

      ndx[i] = n;
      n++;
    }

    if (0 == n)
    {
      cerr << "NN_Results_Base::build_from_distance_matrix: no items selected s = " << only_use.size () << endl;
      return 0;
    }
  }

  similarity_type_t * closest_distances;

  if (dmc.max_neighbours ())
  {
    closest_distances = new similarity_type_t[dmc.max_neighbours ()];
    cerr << "For cross validation, must keep at least " << dmc.max_neighbours () << " neighbours\n";
  }
  else
    closest_distances = nullptr;

  pair<int, float> * p = new pair<int, float>[_number_items];

  int rc = _build_from_distance_matrix (dm, dmc, ndx, p, closest_distances);

  delete [] p;

  if (NULL != closest_distances)
    delete closest_distances;

  return rc;
}

/*
  If we need a minimum number of neighbours, we keep an array of the shortest distances encountered so far.
  When a new distance is encountered, we return whether or not it is one of the N shortest so far
*/

static int
is_one_of_shortest_distances_so_far (similarity_type_t * zbest,
                                     similarity_type_t d,
                                     int n)
{
  return 1;
  if (d > zbest[n - 1])    
    return 0;

  for (int i = 0; i < n; i++)   // simple linear search - the number of neighbours is usually small
  {
    if (d > zbest[i])
      continue;

    for (int j = n - 1; j > i; j--)
    {
      zbest[j] = zbest[j - 1];
    }

    zbest[i] = d;

    return 1;
  }

  return 1;
}

static int
pair_comparitor_x (const void * v1, const void * v2)
{
  typedef pair<int, float> pif;

  const pair<int, similarity_type_t>  * p1 = (const pif *) (v1);
  const pair<int, similarity_type_t> * p2 = (const pif *) (v2);

  if (p1->second < p2->second)
    return -1;
  if (p1->second > p2->second)
    return 1;

  return 0;
}

template <typename N>
int
NN_Results_Base<N>::_build_from_distance_matrix (IWDistanceMatrixFloat & dm,
                                                const DM_to_NN_Conditions<similarity_type_t> & dmc,
                                                const int * ndx,
                                                pair<int, similarity_type_t> * p,
                                                similarity_type_t * closest_distances)
{
  for (int i = 0; i < _number_items; i++)
  {
//  const IWString s = dm.id (i);

    if (ndx[i] < 0)
      continue;

    set_vector (closest_distances, dmc.max_neighbours (), static_cast<similarity_type_t> (1.0));
     
    N & ni = _results[i];

    ni.set_id (dm.id (i));

    int nbrs = 0;    // number of neighbours for item
    for (int j = 0; j < _number_items; j++)
    {
      if (ndx[j] < 0)    // not being processed
        continue;

      if (j == i)   // no self neighbours
        continue;

      similarity_type_t d = dm.zvalue (i, j);

//#define DEBUG_BUILD_NEIGHBOUR_LIST
#ifdef DEBUG_BUILD_NEIGHBOUR_LIST
      cerr << "Between " << i << " and " << j << " is " << d << endl;
#endif

      if (dmc.max_neighbours ())
      {
        if (! is_one_of_shortest_distances_so_far (closest_distances, d, dmc.max_neighbours ()))
          continue;
      }
      else if (d > dmc.max_distance ())
        continue;

      pair<int, similarity_type_t> & pi = p[nbrs];
      pi.first = ndx[j];
      pi.second = d;
      nbrs++;
    }

    if (0 == nbrs)
      continue;

    if (nbrs > 1)
      qsort (p, nbrs, sizeof (pair<int, similarity_type_t>), pair_comparitor_x);

//  Get rid of all neighbours that are more than dmc.max_neighbours () and further away than dmc.max_distance ()

    if (1)    // had problems with the following, maybe fix later...
      ;
    else if (0 == dmc.max_neighbours ())    // neighbour list based just on distances, is OK
      ;
    else if (static_cast<similarity_type_t> (0.0) == dmc.max_distance ())   // no maximum distance, just number of neibhbours
      ;
    else    // both max_neighbours and max_distance specified. Maybe do some trimming
    {
      for (int j = dmc.max_neighbours (); j < nbrs; j++)
      {
        const pair<int, similarity_type_t> & pj = p[j];
        if (pj.second > dmc.max_distance ())
        {
          nbrs = j;
          break;
        }
      }
    }

    ni.build_from_neighbour_list (p, nbrs);

#ifdef DEBUG_BUILD_NEIGHBOUR_LIST
    cerr << "Built " << ni.id () << " with " << ni.number_elements () << " neighbours, furthest " << ni.last_item ()->distance () << endl;

    similarity_type_t prevd = ni[0]->distance ();
    for (int j = 1; j < nbrs; j++)
    {
      similarity_type_t d = ni[j]->distance ();
      if (d < prevd)
      {
        cerr << "Neighbour list out of order\n";
        abort ();
      }

      prevd = d;
    }
#endif
  }

  return _number_items;
}

template <typename N>
int
NN_Results_Base<N>::build_from_gfp_nearneighbours_dash_o(const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "NN_Results_Base::build_from_gfp_nearneighbours_dash_o: cannot open '" << fname << "'\n";
    return 0;
  }

  return build_from_gfp_nearneighbours_dash_o(input);
}

template <typename N>
int
NN_Results_Base<N>::build_from_gfp_nearneighbours_dash_o(iwstring_data_source & input)
{
  if (_number_items > 0)
    ;
  else if (! _count_items_in_input_and_resize (input))
    return 0;

  int fatal = 0;
  for (int i = 0; i < _number_items; i++)
  {
/*GH    cerr << "build function error" << endl;
    return 0;      
*/
    if (_results[i].build (input))
      continue;

    cerr << "NN_Results_Base::build: error building pool item " << i << endl;
    return 0;
  }

  return 1;
}

template <typename N> template <typename F>
void
NN_Results_Base<N>::change_distances (F & f)
{
  for (int i = 0; i < _number_items; i++)
  {
    _results[i].change_distances (f);
  }

  return;
}

#endif
