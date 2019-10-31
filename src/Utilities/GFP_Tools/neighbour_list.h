#ifndef NEIGHBOUR_LIST_H
#define NEIGHBOUR_LIST_H

#include "smiles_id_dist.h"
//#include "neighbour_list.h"

class SID_Comparator
{
  private:
  public:
    int operator() (const Smiles_ID_Dist *, const Smiles_ID_Dist *) const;
};

/*
  Parameter D is the type for the distance measurement, probably float.

  The type H is the object type for the haystack - the things that will
  be used to determine the distances

  The type N is not being used right now - I think I had in mind that
  this would be the object type for a neighbour, but right now, that is
  hard coded as a Smiles_ID_Dist object
*/

template <typename D, typename N, typename H>
class Neighbour_List
{
  private:
    resizable_array_p<Smiles_ID_Dist> _neighbours;

    int _neighbours_to_find;

    D _maxd;

    int _verbose;

    int _keep_going_after_fatal_error;

    int _fatal_errors_encountered;

    int _sort_neighbour_list_at_end;

//  private functions

    int  _binary_search (D) const;
    void _extra_no_max_number_neighbours (const IWString & smiles, const IWString & id, D distance);
    void _check_insertion ();

  public:
    Neighbour_List ();
    ~Neighbour_List ();

    void set_keep_going_after_fatal_error (int s) { _keep_going_after_fatal_error = s;}

    void set_verbose (int s) { _verbose = s;}
    void set_sort_neighbour_list_at_end (int s) { _sort_neighbour_list_at_end = s;}

    int set_neighbours_to_find (int);

    int number_neighbours () const;

    void sort_neighbour_list ();

    void extra (const H &, D);
    int  extra (const IWString & smiles, const IWString & id, D d);
    int  write (std::ostream &) const;
    int  write (IWString_and_File_Descriptor &) const;

    void remove_distant_neighbours (int, similarity_type_t);

    similarity_type_t distance_of_closest_neighbour () const { return _neighbours[0]->distance ();}
    similarity_type_t distance_of_furthest_neighbour () const { return _neighbours.last_item()->distance ();}

    void  shrink(int s);
};

extern int initialise_string_distances();

#ifdef NEIGHBOUR_LIST_IMPLEMENTATION

template <typename D, typename N, typename H>
Neighbour_List<D, N, H>::Neighbour_List ()
{
  _neighbours_to_find = 0;

  _maxd = static_cast<D>(2.0);

  _verbose = 0;

  _keep_going_after_fatal_error = 0;

  _fatal_errors_encountered = 0;

  _sort_neighbour_list_at_end = 0;

  return;
}

template <typename D, typename N, typename H>
Neighbour_List<D, N, H>::~Neighbour_List ()
{
  _neighbours_to_find = -2;
}

template <typename D, typename N, typename H>
int
Neighbour_List<D, N, H>::set_neighbours_to_find (int s)
{
  assert (s > 0);

  _neighbours_to_find = s;

  if (_neighbours.number_elements ())    // clean out any existing neighbour data
    _neighbours.resize (0);

  _neighbours.resize (s);

  for (int i = 0; i < _neighbours_to_find; i++)
  {
    Smiles_ID_Dist * s = new Smiles_ID_Dist ();
    _neighbours.add (s);
  }

  return 1;
}

/*
  We need to compute the number of neighbours because we allocate items 
  that may not have been filled with anything
*/

template <typename D, typename N, typename H>
int
Neighbour_List<D, N, H>::number_neighbours () const
{
  for (int i = _neighbours.number_elements () - 1; i >= 0; i--)
  {
    if (_neighbours[i]->id ().length () > 0)
      return i + 1;
  }

  return 0;
}

template <typename D, typename N, typename H>
int
Neighbour_List<D, N, H>::write (std::ostream & os) const
{
  IWString output_buffer;

  for (int i = 0; i < _neighbours.number_elements (); i++)
  {
    const Smiles_ID_Dist & sidi = *(_neighbours[i]);

    if (0 == sidi.id ().length ())     // neighbour not initialised
      break;

    if (smiles_tag.length ())
      output_buffer << smiles_tag << sidi.smiles () << ">\n";

    output_buffer << identifier_tag << sidi.id () << ">\n";

    output_buffer << distance_tag;
//  append_string_distance (sidi.distance (), output_buffer);
    output_buffer << sidi.distance();
    output_buffer << ">\n";
  }

  os << output_buffer;

  return os.good ();
}

template <typename D, typename N, typename H>
int
Neighbour_List<D, N, H>::write (IWString_and_File_Descriptor & output) const
{
  int n = _neighbours.number_elements();

  for (int i = 0; i < n; i++)
  {
    const Smiles_ID_Dist & sidi = *(_neighbours[i]);

    if (0 == sidi.id ().length ())     // neighbour not initialised
      break;

    if (smiles_tag.length ())
      output << smiles_tag << sidi.smiles () << ">\n";

    output << identifier_tag << sidi.id () << ">\n";

    output << distance_tag;
//  append_string_distance (sidi.distance (), output);
    output << sidi.distance();
    output << ">\n";
    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good ();
}

template <typename T, typename N, typename H>
void
Neighbour_List<T, N, H>::shrink (int s)
{
  if (_neighbours.number_elements() > s)
    _neighbours.resize_keep_storage(s);

  return;
}

/*
  If we have a minimum number of neighbours and an upper distance
  threshold, we often get neighbour lists with > the minimum number
  of neighbours, but some distances that are too far - the array was filled
  with distant neighbours before it got to the min size
*/

template <typename D, typename N, typename H>
void
Neighbour_List<D, N, H>::remove_distant_neighbours (int min_nbrs,
                                              similarity_type_t max_dist)
{
  int n = _neighbours.number_elements();

  if (n < min_nbrs)
    return;

  if (_neighbours[n - 1]->distance() <= max_dist)
    return;

  for (int i = min_nbrs; i < n; i++)
  {
    if (_neighbours[i]->distance() > max_dist)
    {
      _neighbours.resize(i);
      return;
    }
  }

  return;
}

template <typename D, typename N, typename H>
void
Neighbour_List<D, N, H>::sort_neighbour_list ()
{
  if (_verbose > 2)
    cerr << "Neighbour_List::sorting " << _neighbours.number_elements() << " nbrs\n";

  SID_Comparator sid_comparator;
  _neighbours.iwqsort(sid_comparator);

  if (_verbose > 2)
   cerr << "Sort complete, resizing\n";

  if (_neighbours_to_find > 0)
    _neighbours.resize(_neighbours_to_find);

  if (_verbose > 2)
    cerr << "Resize complete\n";

  return;
}

template <typename D, typename N, typename H>
int
Neighbour_List<D, N, H>::_binary_search (D distance) const
{
  int left = 0;      // _neighbours[left] will hold the new value

  if (distance > _neighbours[0]->distance ())     // goes after neighbour 0
  {
    int right = _neighbours.number_elements () - 1;

    while (right > left)
    {
      int middle = (left + right) / 2;
      if (middle == left)
      {
        left = right;
        break;
      }

      D dmid = _neighbours[middle]->distance ();
//    cerr << "Left = " << left << " d = " << _neighbours[left]->distance () << " middle " << middle << " d = " << _neighbours[middle].distance () << " right " << right << " d = " << _neighbours[right].distance () << endl;
      if (distance < dmid)
        right = middle;
      else if (distance > dmid)
        left = middle;
      else
      {
        left = middle;
        break;
      }
    }
  }

  return left;
}

template <typename D, typename N, typename H>
void
Neighbour_List<D, N, H>::_check_insertion () 
{
  int failure = 0;

  int nn = _neighbours.number_elements ();
  for (int i = 1; i < nn; i++)
  {
    if (_neighbours[i - 1]->distance () > _neighbours[i]->distance ())
    {
      cerr << "Sort/insertion failed, out of order, i = " << i << endl;
      cerr << _neighbours[i - 1]->distance () << " vs " << _neighbours[i]->distance () << endl;
      if (_keep_going_after_fatal_error)
      {
        _neighbours.remove_item (i);
        nn--;
      }
      failure = 1;
    }
  }

  if (failure && _keep_going_after_fatal_error)
  {
    _fatal_errors_encountered++;
    return;
  }

  if (failure)
  {
    for (int i = 0; i < _neighbours.number_elements (); i++)
    {
      cerr << "i = " << i << " distance " << _neighbours[i]->distance () << endl;
    }

    exit (87);
  }

  return;
}


template <typename D, typename N, typename H>
void
Neighbour_List<D, N, H>::extra (const H & rhs, D distance)
{
  if (distance >= static_cast<D> (0.0) && distance <= static_cast<D> (1.0))
    ;
  else
  {
    cerr << "Neighbour_List::extra: fatal error, distance " << distance << " to '" << rhs.id () << "'\n";
    if (_keep_going_after_fatal_error)
    {
      _fatal_errors_encountered++;
      return;
    }

    abort ();
  }

//cerr << "d = " << distance << " at end? " << _sort_neighbour_list_at_end << endl;
  if (_sort_neighbour_list_at_end)
  {
    Smiles_ID_Dist * sidi = new Smiles_ID_Dist (rhs.smiles (), rhs.id (), distance);
    _neighbours.add(sidi);
    return;
  }

//cerr << "d = " << distance << " to find " << _neighbours_to_find << endl;
  if (0 == _neighbours_to_find)
  {
    _extra_no_max_number_neighbours (rhs.smiles(), rhs.id(), distance);
    return;
  }

  if (distance >= _maxd)
    return;

  Smiles_ID_Dist * x = _neighbours.last_item ();

  int left = _binary_search (distance);

// Shuffle everything one slot to the right

  for (int i = _neighbours.number_elements () - 1; i > left; i--)
  {
    _neighbours[i] = _neighbours[i - 1];
  }

  x->set_distance (distance);
  x->set_id (rhs.id ());
  x->set_smiles (rhs.smiles ());

  _neighbours[left] = x;

  _maxd = _neighbours[_neighbours_to_find - 1]->distance();

//#define CHECK_INSERTION
#ifdef CHECK_INSERTION
  _check_insertion ();
#endif

  return;
}

template <typename D, typename N, typename H>
int
Neighbour_List<D, N, H>::extra (const IWString & smiles, const IWString & id, D distance)
{
#ifdef CHECK_DISTANCES_IN_NBR_LIST
  if (distance >= static_cast<D> (0.0) && distance <= static_cast<D> (1.0))
    ;
  else
  {
    cerr << "Neighbour_List::extra: fatal error, distance " << distance << " to '" << id << "'\n";
    if (_keep_going_after_fatal_error)
    {
      _fatal_errors_encountered++;
      return;
    }

    abort ();
  }
#endif

//cerr << "d = " << distance << " at end? " << _sort_neighbour_list_at_end << endl;
  if (_sort_neighbour_list_at_end)
  {
    Smiles_ID_Dist * sidi = new Smiles_ID_Dist (smiles, id, distance);
    _neighbours.add(sidi);
    return 1;
  }

//cerr << "d = " << distance << " to find " << _neighbours_to_find << endl;
  if (0 == _neighbours_to_find)
  {
    _extra_no_max_number_neighbours (smiles, id, distance);
    return 1;
  }

  if (distance >= _maxd)
    return 0;

  Smiles_ID_Dist * x = _neighbours.last_item ();

  int left = _binary_search (distance);

// Shuffle everything one slot to the right

  for (int i = _neighbours.number_elements () - 1; i > left; i--)
  {
    _neighbours[i] = _neighbours[i - 1];
  }

  x->set_distance (distance);
  x->set_id (id);
  x->set_smiles (smiles);

  _neighbours[left] = x;

  _maxd = _neighbours[_neighbours_to_find - 1]->distance();

//#define CHECK_INSERTION
#ifdef CHECK_INSERTION
  _check_insertion ();
#endif

  return 1;
}


/*
  We are accepting all neighbours
*/

template <typename D, typename N, typename H>
void
Neighbour_List<D, N, H>::_extra_no_max_number_neighbours (const IWString & smiles,
                                                          const IWString & id,
                                                          D distance)
{
  Smiles_ID_Dist * sidi = new Smiles_ID_Dist (smiles, id, distance);

  if (0 == _neighbours.number_elements ())
  {
    _neighbours.add (sidi);

    return;
  }

  if (distance >= _neighbours.last_item ()->distance ())
  {
    _neighbours.add (sidi);

    return;
  }

  int ins = _binary_search (distance);

  _neighbours.insert_before (ins, sidi);

#ifdef CHECK_INSERTION
  _check_insertion ();
#endif

  return;
}

#endif
#endif
