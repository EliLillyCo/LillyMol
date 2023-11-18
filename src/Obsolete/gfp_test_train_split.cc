/*
  Divide a set of molecules into disparate sets
  Reads the output from gfp_nearneighbours
*/

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/histogram/iwhistogram.h"

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString neighbour_tag ("NBR<");
static IWString distance_tag ("DIST<");

using std::cerr;
using std::endl;

typedef float similarity_type_t;

/*
  The cost of being a neighbour to someone is a function of how close we are to them
*/

static similarity_type_t initial_threshold = 0.0;

static int
compute_cost (similarity_type_t d)
{
  return static_cast<int> ((1.0 - d) * 20.0001);
}

class Test_Set_Member;

/*
  As neighbours are read in, they have an index of their position in the set
  That will later be converted to a pointer to the appropriate Test_Set_Member
  item
*/

class Neighbour
{
  private:
    int _index_in_set;

    Test_Set_Member * _parent;

    similarity_type_t _distance;

//  A threshold defines a cost

    int _cost;

  public:
    Neighbour ();

    void initialise (int, similarity_type_t);

    int debug_print (std::ostream &) const;

    Test_Set_Member * parent () { return _parent;}
    const Test_Set_Member * parent () const { return _parent;}

    int index_in_set () const { return _index_in_set;}

    similarity_type_t distance () const { return _distance;}
    int cost () const { return _cost;}

    inline int set_membership () const;
    void set_set_membership (int);

    int establish_parent_pointer (Test_Set_Member * pool, int pool_size);
    int build (const const_IWSubstring &, iwstring_data_source &);
};

Neighbour::Neighbour ()
{
  _index_in_set = -1;
  _parent = nullptr;

  return;
}

void
Neighbour::initialise (int ndx,
                       similarity_type_t d)
{
  _index_in_set = ndx;
  _distance = d;
  _cost = compute_cost (_distance);

  return;
}

static int
fetch_tdt_value (const const_IWSubstring & buffer,
                 IWString & zresult)
{
  assert (buffer.ends_with ('>'));

  int open_angle_bracket = buffer.index ('<');

  zresult.strncpy (buffer.rawchars () + open_angle_bracket + 1, buffer.length () - open_angle_bracket - 2);

  return 1;
}

template <typename T>
int
fetch_numeric_tdt_value (const const_IWSubstring & buffer,
                         T & zvalue)
{
  IWString tmp;

  if (! fetch_tdt_value (buffer, tmp))
  {
    cerr << "Cannot extract TDT value from '" << buffer << "'\n";
    return 0;
  }

  if (! tmp.numeric_value (zvalue))
  {
    cerr << "Invalid numeric '" << tmp << "'\n";
    return 0;
  }

  return 1;
}

/*int
Neighbour::build (const const_IWSubstring & buffer,
                  iwstring_data_source & input)
{
  if (! buffer.starts_with (neighbour_tag))
  {
    cerr << "Neighbour::build: first record must be index '" << buffer << "'\n";
    return 0;
  }

  if (! fetch_numeric_tdt_value (buffer, _index_in_set) || _index_in_set < 0)
  {
    cerr << "Invalid identifier '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring mybuffer;
  if (! input.next_record (mybuffer) || ! mybuffer.starts_with (distance_tag))
  {
    cerr << "Neighbour::build: invalid distance specification '" << mybuffer << "'\n";
    return 0;
  }

  if (! fetch_numeric_tdt_value (mybuffer, _distance) || _distance < 0.0 || _distance > 1.0)
  {
    cerr << "Neighbour::build: invalid distance '" << buffer << "'\n";
    return 0;
  }

  _cost = compute_cost (_distance);

  return 1;
}*/

class Test_Set_Member
{
  private:
    IWString _smiles;
    IWString _id;

    int _singleton;

    int _set_membership;

    int _number_neighbours;
    Neighbour * _neighbour;

//  private functions

    int _build (const IW_TDT & tdt);

  public:
    Test_Set_Member ();
    ~Test_Set_Member ();

    const IWString & id () const { return _id;}

    int number_neighbours () const { return _number_neighbours;}

    int set_all_neighbours_to_set (int s, int maxn);

    int determine_singleton_status (similarity_type_t);

    void cost_of_set (int set_number, similarity_type_t & closest, similarity_type_t & total_cost) const;

    int cost_of_changing (int c, similarity_type_t threshold) const;

    const Neighbour & neighbour (int i) const { return _neighbour[i];}

    similarity_type_t closest_distance () const;

    void set_set_membership (int s) { _set_membership = s;}
    int  set_membership () const { return _set_membership;}

    int build (iwstring_data_source &, int &);

    int establish_parent_pointers (Test_Set_Member *, int);

    int distances_to_class (int, similarity_type_t, similarity_type_t &) const;

    int distances_within (similarity_type_t) const;

    int cost_in_class (int, similarity_type_t) const;

    int neighbours_same_class_as_you (int);

    int closest_distance_to_class (int, similarity_type_t, similarity_type_t &) const;
};

Test_Set_Member::Test_Set_Member ()
{
  _set_membership = -1;
  _singleton = -1;
  _number_neighbours = 0;
  _neighbour = nullptr;

  return;
}

Test_Set_Member::~Test_Set_Member ()
{
  if (-3 == _number_neighbours)
    cerr << "Test_Set_Member::~Test_Set_Member: deleting already deleted Test_Set_Member\n";

  _number_neighbours = -3;

  if (NULL != _neighbour)
    delete [] _neighbour;

  return;
}

similarity_type_t
Test_Set_Member::closest_distance () const
{
  if (0 == _number_neighbours)
    return static_cast<similarity_type_t> (1.0);

  return _neighbour[0].distance ();
}

int
Test_Set_Member::build (iwstring_data_source & input, 
                        int & fatal)
{
  IW_TDT tdt;

  if (! tdt.next (input))
  {
    if (input.eof ())
    {
      fatal = 0;
      return 1;
    }

    fatal = 1;
    cerr << "Test_Set_Member::build: cannot construct tdt, now at line " << input.lines_read () << endl;
    return 0;
  }

  fatal = 1;

  int rc = _build (tdt);

  if (0 == rc)
  {
    cerr << "Test_Set_Member::build: invalid data, now at line " << input.lines_read () << endl;
    cerr << tdt;
    return 0;
  }

  return rc;
}

int
Test_Set_Member::_build (const IW_TDT & tdt)
{

  if (! tdt.dataitem_value (smiles_tag, _smiles))
  {
    cerr << "Test_Set_Member::build: cannot extract '" << smiles_tag << "' from TDT\n";
    return 0;
  }

  if (! tdt.dataitem_value (identifier_tag, _id))
  {
    cerr << "Test_Set_Member:build: cannot extract '" << identifier_tag << "' from TDT\n";
    return 0;
  }

  _number_neighbours = tdt.count_dataitems (neighbour_tag);
  if (0 == _number_neighbours)
  {
    _neighbour = nullptr;
    return 1;
  }

  _neighbour = new Neighbour[_number_neighbours];
  assert (NULL != _neighbour);

  for (int i = 0; i < _number_neighbours; i++)
  {
    int nbr;
    if (! tdt.dataitem_value (neighbour_tag, nbr, i) || nbr < 0)
    {
      cerr << "Test_Set_Member::build: missing or invalid '" << neighbour_tag << "' in TDT\n";
      return 0;
    }

    similarity_type_t d;
    if (! tdt.dataitem_value (distance_tag, d, i) || d < static_cast<similarity_type_t> (0.0) || d > static_cast<similarity_type_t> (1.0))
    {
      cerr << "Test_Set_Member::build: missing or invalid '" << distance_tag << "' in TDT\n";
      return 0;
    }

    _neighbour[i].initialise (nbr, d);
  }

  return 1;
};

int
Test_Set_Member::establish_parent_pointers (Test_Set_Member * pool,
                                            int pool_size)
{
  for (int i = 0; i < _number_neighbours; i++)
  {
    if (! _neighbour[i].establish_parent_pointer (pool, pool_size))
    {
      return 0;
    }
  }

  return 1;
}

int
Test_Set_Member::determine_singleton_status (similarity_type_t threshold)
{
  const Neighbour & n0 = _neighbour[0];

  if (n0.distance () >= threshold)
    _singleton = 1;
  else
    _singleton = 0;

  return _singleton;
}

int 
Neighbour::establish_parent_pointer (Test_Set_Member * pool,
                                     int pool_size)
{
  assert (_index_in_set >= 0);

  if (_index_in_set >= pool_size)
  {
    cerr << "Neighbour::establish_parent_pointer: invalid index " << _index_in_set << " vs " << pool_size << "'\n";
    return 0;
  }

  _parent = & (pool[_index_in_set]);

  return 1;
}

int
Neighbour::set_membership () const
{
  return _parent->set_membership ();
}

void
Neighbour::set_set_membership (int s)
{
  _parent->set_set_membership (s);
}

static Test_Set_Member * pool = nullptr;
static int pool_size = 0;

/*
  Because some molecules may have no neighbours, we can often deal with the effective pool size,
  which is the number of molecules that have neighbours within the threshold
*/

static int effective_pool_size = 0;

static int verbose = 0;

static float fraction_in_set_1 = 0.5;

static int items_needed_in_set_1 = 0;

/*
  Once we know the number of molecules within the threshold, we can then determine
  the minimum number of molecules that can be assigned to set 1. Any difference
  between the mininum and the desired number can be made up by using molecules
  with no neighbours within threshold
*/

static int min_items_in_set1 = 0;
static int max_items_in_set1 = 0;

/*
  When doing iterative refinement, we create a number of copies of our original set memberships,
  and then do a number of steps on those
*/

static int iterative_refinement_replicates = 0;

static int iterative_refinement_steps = 0;

/*
  We may write files of histograms of nearest distances
*/

static IWString histogram_name_stem;

static int
compute_items_needed_in_set1 ()
{
  int rc = int (static_cast<float> (pool_size) * fraction_in_set_1);
  if (0 == rc)
    rc = 1;

  return rc;
}

static int
compute_percent (int numerator, int denominator)
{
  int rc = static_cast<int> (numerator * 100 / denominator);

  return rc;
}

/*
  In order to know how we improve things, it can be useful to generate some random
  splits
*/

static int nrandom = 0;

static int
build_pool (iwstring_data_source & input)
{
  assert (pool_size > 0);

  for (int i = 0; i < pool_size; i++)
  {
    int fatal;
    if (! pool[i].build (input, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot build pool member " << i << "\n";
        return 0;
      }

      pool_size = i;
      break;
    }

  }

  return 1;
}

static int
build_pool (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size)
  {
    assert (NULL == pool);

    IWString tmp ("^\\$SMI<");

    pool_size = input.grep (tmp);
    if (0 == pool_size)
    {
      cerr << "Yipes, cannot find any '" << tmp << "' in the input\n";
      return 0;
    }

    if (verbose)
      cerr << "Input contains " << pool_size << " fingerprints\n";

    pool = new Test_Set_Member[pool_size];

    if (NULL == pool)
    {
      cerr << "Yipes, cannot allocate space for " << pool_size << " fingerprints\n";
      return 0;
    }
  }

  return build_pool (input);
}

/*
  for a molecule that hasn't been placed into a set, we can ask what is the
  distance penalty of doing that
*/

void
Test_Set_Member::cost_of_set (int set_number,
                              similarity_type_t & closest,
                              similarity_type_t & total_cost) const
{
  assert (_set_membership < 0);

  total_cost = static_cast<similarity_type_t> (0.0);
  closest = static_cast<similarity_type_t> (1.0);

  for (int i = 0; i < _number_neighbours; i++)
  {
    const Neighbour & ni = _neighbour[i];

    if (ni.distance () >= initial_threshold)
      break;

    const Test_Set_Member * pi = ni.parent ();

    if (pi->set_membership () != set_number)
      continue;

    total_cost += ni.distance ();
    if (static_cast<similarity_type_t> (1.0) == closest)
      closest = ni.distance ();
  }

  return;
}

/*
  Once we know the number of molecules that have neighbours, we can then work out the
  minimum numbers allowed during optimisation
*/

static int
determine_molecules_with_neighbours_within_threshold (similarity_type_t threshold)
{
  int number_molecules_with_neighbours_within_threshold = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i].closest_distance () <= threshold)
      number_molecules_with_neighbours_within_threshold++;
  }

  if (0 == number_molecules_with_neighbours_within_threshold)
  {
    cerr << "Yipes, none of " << pool_size << " molecules have neighbours closer than " << threshold << endl;
    return 0;
  }

  if (verbose)
    cerr << number_molecules_with_neighbours_within_threshold << " of " << pool_size << " molecules have neighbours within the " << threshold << " radius\n";

  min_items_in_set1 = items_needed_in_set_1 - (pool_size - number_molecules_with_neighbours_within_threshold);
  if (min_items_in_set1 < 0)
    min_items_in_set1 = 0;

  if (verbose)
    cerr << "At least " << min_items_in_set1 << " items in set 1\n";

  max_items_in_set1 = items_needed_in_set_1 + (pool_size - number_molecules_with_neighbours_within_threshold);
  if (max_items_in_set1 >= pool_size)
    max_items_in_set1 = pool_size;

  return 1;
}

static int
establish_parent_pointers ()
{
  int rc = pool_size;

  for (int i = 0; i < pool_size; i++)
  {
    if (! pool[i].establish_parent_pointers (pool, pool_size))
    {
      cerr << "Cannot establish neighbour pointers for '" << pool[i].id () << "'\n";
      rc = 0;
    }
  }

  return rc;
}

static int
Test_Set_Member_comparitor_closest_distance (const void * vt1, const void * vt2)
{
  const Test_Set_Member * t1 = (const Test_Set_Member *) vt1;
  const Test_Set_Member * t2 = (const Test_Set_Member *) vt2;

  similarity_type_t d1 = t1->closest_distance ();
  similarity_type_t d2 = t2->closest_distance ();

  if (d1 < d2)
    return -1;
  if (d1 > d2)
    return 1;

  return 0;
}

static int
Test_Set_Member_comparitor_number_neighbours (const void * vt1, const void * vt2)
{
  const Test_Set_Member * t1 = (const Test_Set_Member *) vt1;
  const Test_Set_Member * t2 = (const Test_Set_Member *) vt2;

  int n1 = t1->number_neighbours ();
  int n2 = t2->number_neighbours ();

  if (n1 < n2)
    return 1;
  if (n1 > n2)
    return -1;

  return 0;
}

int
Test_Set_Member::distances_to_class (int c,
                                     similarity_type_t threshold,
                                     similarity_type_t & closest_distance) const
{
  int rc = 0;    // the number of members of the other class within THRESHOLD

  for (int i = 0; i < _number_neighbours; i++)
  {
    const Neighbour & ni = _neighbour[i];

    if (ni.distance () > threshold)
      break;

    if (ni.set_membership () != c)    // we are interested in distances to members of class C
      continue;

    if (0 == rc)    // first one
      closest_distance = ni.distance ();

    rc++;
  }

  return rc;
}

int
Test_Set_Member::closest_distance_to_class (int c,
                                            similarity_type_t threshold,
                                            similarity_type_t & closest_distance_to_class) const
{
  for (int i = 0; i < _number_neighbours; i++)
  {
    const Neighbour & ni = _neighbour[i];

    if (c != ni.set_membership ())
      return 0;

    closest_distance_to_class = ni.distance ();

    return 1;
  }

  return 0;
}

int
Test_Set_Member::distances_within (similarity_type_t threshold) const
{
  int rc = 0;

  for (int i = 0; i < _number_neighbours; i++)
  {
    if (_neighbour[i].distance () <= threshold)
      rc++;
    else
      break;
  }

  return rc;
}

/*
  When deciding whether an item should switch from 1 class to another, we need
  a cost for being in each class.
*/

int
Test_Set_Member::cost_in_class (int c, similarity_type_t threshold) const
{
  int rc = 0;

  for (int i = 0; i < _number_neighbours; i++)
  {
    const Neighbour & ni = _neighbour[i];

    if (c == ni.set_membership ())    // we want to maximise inter-class distances
      continue;

//  similarity_type_t d = ni.distance ();

//  if (d > threshold)
//    break;

    rc += ni.cost ();
  }

  return rc;
}

int
Test_Set_Member::set_all_neighbours_to_set (int s, int maxn)
{
  int istop = _number_neighbours;
  if (istop > maxn)
    istop = maxn;

  for (int i = 0; i < istop; i++)
  {
    Test_Set_Member * p = _neighbour[i].parent ();

    p->set_set_membership (s);
  }

  return istop;
}

/*
  Return the net change in set 1 membership
*/

int
Test_Set_Member::neighbours_same_class_as_you (int n)
{
  assert (n > 0 && n <= _number_neighbours);

  int rc = 0;
  for (int i = 0; i < n; i++)
  {
    if (_set_membership != _neighbour[i].set_membership ())
    {
      _neighbour[i].set_set_membership (_set_membership);
      rc++;
    }
  }

  if (1 == _set_membership)
    return rc;     // we added this number of molecules to set 1
  else
    return -rc;    // we subtracted this number of molecules from set 1
}

/*
  Someone is thinking of changing the set membership of pool member C
  What would be the cost to us of doing that?
*/

int
Test_Set_Member::cost_of_changing (int c, 
                                   similarity_type_t threshold) const
{
  for (int i = 0; i < _number_neighbours; i++)
  {
    const Neighbour & n = _neighbour[i];

    if (n.distance () > threshold)
      return 0;

    if (c != n.index_in_set ())
      continue;

    return n.cost ();
  }

  return 0;   // wasn't one of my neighbours, therefore no cost
}

class Set_Membership
{
  private:
    int _cost;

    int * _set;

//  private functions

    void _compute_cost (similarity_type_t threshold, int * cost);
    void _determine_cost (similarity_type_t);
    void _assign_set_memberships () const;

    int _iterative_refinement_deterministic (int, similarity_type_t, int *, int *);
    int _iterative_refinement_random (int, similarity_type_t, int *, int *, int *);

    int _cost_of_changing (int c, similarity_type_t threshold) const;

  public:
    Set_Membership ();
    ~Set_Membership ();

    Set_Membership & operator = (const Set_Membership &);

    void copy (const Set_Membership &);

    void initialise_randomly ();

    int  adjust_set_memberships_to_desired_numbers (similarity_type_t);

    int cost () const { return _cost;}

    void determine_cost (similarity_type_t);

    void initial_neighbour_groupings (similarity_type_t);

    int items_in_set1 () const { return count_occurrences_of_item_in_array (1, pool_size, _set);}

    int write (similarity_type_t, std::ostream &, std::ostream &) const;
    int write (similarity_type_t, IWString &, IWString &) const;
    int write_close_distance_histogram (similarity_type_t, IWString &) const;
    int write_close_distance_histogram (similarity_type_t, std::ostream &) const;

    int iterative_refinement_deterministic (int, similarity_type_t);
    int iterative_refinement_random (int, similarity_type_t);
};

Set_Membership::Set_Membership ()
{
  _set = nullptr;
  _cost = -1;

  return;
}

void
Set_Membership::initialise_randomly ()
{
  if (NULL == _set)
    _set = new int[pool_size];

  int n1 = 0;

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < pool_size; i++)
  {
    if (u(rng) < fraction_in_set_1) {
      _set[i] = 1;
      n1++;
    }
    else
      _set[i] = 0;
  }

// Switch the classes to make the number of 1's as close to ITEMS_NEEDED_IN_SET_1

  if (n1 == items_needed_in_set_1)    // wow, got it by chance
    return;

  int d1 = n1 - items_needed_in_set_1;

  if (d1 < 0)
    d1 = - d1;

  int n0 = pool_size - n1;

  int d0 = n0 - items_needed_in_set_1;     // the difference between actual and desired

  if (d0 < 0)
    d0 = - d0;

  if (d1 > d0)
  {
    for (int i = 0; i < pool_size; i++)
    {
      if (0 == _set[i])
        _set[i] = 1;
      else
        _set[i] = 0;
    }
  }

  return;
}

Set_Membership &
Set_Membership::operator = (const Set_Membership & rhs)
{
  Set_Membership::copy (rhs);

  return *this;
}

void
Set_Membership::copy (const Set_Membership & rhs)
{
  if (NULL == _set)
    _set = new int[pool_size];

  copy_vector (_set, rhs._set, pool_size);

  _cost = rhs._cost;

  return;
}

Set_Membership::~Set_Membership ()
{
  if (NULL != _set)
    delete _set;

  return;
}

void
Set_Membership::_assign_set_memberships () const
{
  for (int i = 0; i < pool_size; i++)
  {
    pool[i].set_set_membership (_set[i]);
  }

  return;
}

void
Set_Membership::determine_cost (similarity_type_t threshold)
{
  _assign_set_memberships ();

  _determine_cost (threshold);

  return;
}

/*
  This version assumes that the set memberships are already assigned
*/

void
Set_Membership::_determine_cost (similarity_type_t threshold)
{
  _cost = 0;

  for (int i = 0; i < pool_size; i++)
  {
    _cost += pool[i].cost_in_class (_set[i], threshold);
  }

  return;
}

/*
  Sometimes we need an array of costs
*/

void
Set_Membership::_compute_cost (similarity_type_t threshold,
                               int * cost)
{
  for (int i = 0; i < pool_size; i++)
  {
    cost[i] = pool[i].cost_in_class (_set[i], threshold);
  }

  _cost = sum_vector (cost, pool_size);

  return;
}
void
Set_Membership::initial_neighbour_groupings (similarity_type_t threshold)
{
  _assign_set_memberships ();

  int s1 = items_in_set1 ();

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> u(2, pool_size - 1);

  for (int i = 0; i < 20; i++)
  {
    for (int j = 0; j < pool_size; j++)
    {
      Test_Set_Member & pj = pool[j];

      int n = pj.distances_within (threshold);

      if (n < 6)
        continue;

      if (s1 >= max_items_in_set1 && 1 == pj.set_membership ())   // cannot add more items to set 1
        continue;

      if (s1 <= min_items_in_set1 && 0 == pj.set_membership ())    // cannot remove items from set 1
        continue;

      int k = u(rng) / 2;

      s1 += pj.neighbours_same_class_as_you (k);
    }
  }

  for (int i = 0; i < pool_size; i++)
  {
    _set[i] = pool[i].set_membership ();
  }

  _determine_cost (threshold);

  return;
}

int
Set_Membership::write (similarity_type_t threshold,
                       IWString & test_fname,
                       IWString & train_fname) const
{
  std::ofstream test_output (test_fname.null_terminated_chars (), std::ios::out);
  std::ofstream train_output (train_fname.null_terminated_chars (), std::ios::out);

  if (! test_output.good () || ! train_output.good ())
  {
    cerr << "Set_Membership::write: cannot open '" << test_fname << "' and/or '" << train_fname << "'\n";
    return 0;
  }

  return write (threshold, test_output, train_output);
}

int
Set_Membership::adjust_set_memberships_to_desired_numbers (similarity_type_t threshold)
{
  int s1 = items_in_set1 ();

  if (verbose > 1)
    cerr << "Set initially has " << s1 << " items in set 1\n";
  
  if (s1 < items_needed_in_set_1)
  {
    for (int i = 0; i < pool_size && s1 < items_needed_in_set_1; i++)
    {
      if (1 == _set[i])    // we want to move things into set1
        continue;
  
      Test_Set_Member & pi = pool[i];
  
      if (pi.closest_distance () <= threshold)    // has been optimised
        continue;
  
      _set[i] = 1;
      s1++;
    }
  }
  else if (s1 > items_needed_in_set_1)
  {
    for (int i = 0; i < pool_size && s1 > items_needed_in_set_1; i++)
    {
      if (0 == _set[i])    // we want to move things out of set1
        continue;
  
      Test_Set_Member & pi = pool[i];
  
      if (pi.closest_distance () <= threshold)    // has been optimised
        continue;
  
      _set[i] = 0;
      s1--;
    }
  }

  if (verbose > 1)
    cerr << "Set now has " << s1 << " items in set 1\n";

  return 1;
}

int
Set_Membership::write (similarity_type_t threshold,
                       std::ostream & test_output,
                       std::ostream & train_output) const
{
  _assign_set_memberships ();

  for (int i = 0; i < pool_size; i++)
  {
    const Test_Set_Member & pi = pool[i];

    int cost = pi.cost_in_class (_set[i], threshold);

    if (1 == _set[i])
      train_output << pi.id () << ' ' << cost << endl;
    else
      test_output  << pi.id () << ' ' << cost << endl;
  }

  return train_output.good () && test_output.good ();
}

int
Set_Membership::write_close_distance_histogram (similarity_type_t threshold,
                                                IWString & fname) const
{
  std::ofstream output (fname.null_terminated_chars (), std::ios::out);

  if (! output.good ())
  {
    cerr << "Set_Membership::write_close_distance_histogram: cannot open '" << fname << "'\n";
    return 0;
  }

  return write_close_distance_histogram (threshold, output);
}

int
Set_Membership::write_close_distance_histogram (similarity_type_t threshold,
                                                std::ostream & output) const
{
  _assign_set_memberships ();

  IWHistogram h;
  h.initialise (0.0, 1.0, 0.01);

  for (int i = 0; i < pool_size; i++)
  {
    const Test_Set_Member & pi = pool[i];

    assert (pi.set_membership () == _set[i]);

    int other_class;    
    if (_set[i])
      other_class = 0;
    else
      other_class = 1;

    similarity_type_t d;
    if (! pi.closest_distance_to_class (other_class, threshold, d))
      continue;

    h.extra (d);
  }

  return h.write (output);
}

/*
  We are thinking of changing the set membership of item C
*/

int
Set_Membership::_cost_of_changing (int c,
                                   similarity_type_t threshold) const
{

  const Test_Set_Member & p = pool[c];

  int n = p.number_neighbours ();

// Scan the neighbours of C. Note there is a possible flaw with this. Maybe C is a neighbour
// of something that isn't in C's neighbour list.

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    const Neighbour & ni = p.neighbour (i);

    if (ni.distance () > threshold)
      break;

    int ndx = ni.index_in_set ();

    if (_set[c] == _set[ndx])    // in same class, therefore cost increases if we change the class membership
      rc += ni.cost ();
    else                         // now in different classes, cost comes down
      rc -= ni.cost ();
  }

  return rc;
}

static int 
int_comparitor_larger (const void * v1, const void * v2)
{
  int i1 = * ((const int *) v1);
  int i2 = * ((const int *) v2);

  if (i1 < i2)
    return 1;

  if (i1 > i2)
    return -1;

  return 0;
}

static int
random_selection_from_top (const int * tmp, 
                           int n)
{
  int m = n / 10;
  if (m < 2)
    m = 2;

  static std::random_device rd;
  static std::mt19937 rng(rd());

  std::uniform_int_distribution<int> u(0, m);

  return u(rng);
}

static int
a_high_cost_item (const int * cost,
                  int n,
                  const int * already_done,
                  int * tmp)
{
  int m = 0;
  for (int i = 0; i < n; i++)
  {
    if (already_done[i])
      continue;

    tmp[m] = cost[i];
    m++;
  }

  assert (m > 0);

  qsort (tmp, m, sizeof (int), int_comparitor_larger);

// randomly grab one of the top 10%

  int j = random_selection_from_top (tmp, m);

  int jcost = tmp[j];

// Now find the index of something that has this cost

  std::random_device rd;
  std::mt19937 rng(rd());
  std::bernoulli_distribution d(0.50);

  if (d(rng)) {
    for (int i = 0; i < n; i++)
    {
      if (! already_done[i] &&  cost[i] == jcost)
        return i;
    }
  }
  else
  {
    for (int i = n - 1; i >= 0; i--)
    {
      if (! already_done[i] && cost[i] == jcost)
        return i;
    }
  }

// If we come to here, we weren't able to find something what that cost - impossible

  abort ();

  return -1;
}

static int
next_one_to_try (const int * cost,
                 int n,
                 const int * already_done)
{
  int highest_cost = -1;
  int item_with_highest_cost = -1;

  for (int i = 0; i < n; i++)
  {
    if (already_done[i])
      continue;

    if (cost[i] > highest_cost)
    {
      highest_cost = cost[i];
      item_with_highest_cost = i;
    }
  }

  return item_with_highest_cost;
}

int
Set_Membership::_iterative_refinement_deterministic (int nsteps,
                                       similarity_type_t threshold,
                                       int * cost,
                                       int * already_done)
{
  _assign_set_memberships ();

  _compute_cost (threshold, cost);

  int initial_cost = _cost;

  int items_in_set1 = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (_set[i])
      items_in_set1++;
  }

  for (int i = 0; i < nsteps; i++)
  {
    int j = next_one_to_try (cost, pool_size, already_done);
    if (j < 0)
      break;

    already_done[j] = 1;

    Test_Set_Member & pj = pool[j];

    int cost_in_other_class;
    if (_set[j])
      cost_in_other_class = pj.cost_in_class (0, threshold);
    else
      cost_in_other_class = pj.cost_in_class (1, threshold);

    int cost_of_change = _cost_of_changing (j, threshold);

#ifdef DEBUG__ITERATIVE_REFINEMENT_DETERMINISTIC
    cerr << "Highest cost item " << j << " cost " << cost[j] << ". Cost in other " << cost_in_other_class << " to change neighbours " << cost_of_change << endl;
#endif

    if (cost[j] < cost_in_other_class + cost_of_change)
      continue;

    if (1 == _set[j])
    {
      _set[j] = 0;
      items_in_set1--;
    }
    else
    {
      _set[j] = 1;
      items_in_set1++;
    }

    pj.set_set_membership (_set[j]);

    _compute_cost (threshold, cost);

#ifdef DEBUG__ITERATIVE_REFINEMENT_DETERMINISTIC
    cerr << "Total cost now " << _cost << ", " << items_in_set1 << " items in set 1\n";
#endif

    if (items_in_set1 < min_items_in_set1 || items_in_set1 > max_items_in_set1)
      break;
  }

  if (verbose > 1)
    cerr << "Deterministic refinement reduced const from " << initial_cost << " to " << _cost << ", " << items_in_set1 << " items in set1, " << compute_percent (items_in_set1, pool_size) << endl;

  return 1;
}


int
Set_Membership::_iterative_refinement_random (int nsteps,
                                       similarity_type_t threshold,
                                       int * cost,
                                       int * already_done,
                                       int * tmp)
{
  _assign_set_memberships ();

  _compute_cost (threshold, cost);

  int initial_cost = _cost;

  int items_in_set1 = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (_set[i])
      items_in_set1++;
  }

  for (int i = 0; i < nsteps; i++)
  {
    int j = a_high_cost_item (cost, pool_size, already_done, tmp);
    if (j < 0)
      break;

    already_done[j] = 1;

    Test_Set_Member & pj = pool[j];

    int cost_in_other_class;
    if (_set[j])
      cost_in_other_class = pj.cost_in_class (0, threshold);
    else
      cost_in_other_class = pj.cost_in_class (1, threshold);

    int cost_of_change = _cost_of_changing (j, threshold);

#ifdef DEBUG_ITERATIVE_REFINEMENT_RANDOM
    cerr << "Highest cost item " << j << " cost " << cost[j] << ". Cost in other " << cost_in_other_class << " to change neighbours " << cost_of_change << endl;
#endif

    if (cost[j] < cost_in_other_class + cost_of_change)
      continue;

    if (1 == _set[j])
    {
      _set[j] = 0;
      items_in_set1--;
    }
    else
    {
      _set[j] = 1;
      items_in_set1++;
    }

    pj.set_set_membership (_set[j]);

    _compute_cost (threshold, cost);

#ifdef DEBUG_ITERATIVE_REFINEMENT_RANDOM
    cerr << "Total cost now " << _cost << ", " << items_in_set1 << " items in set 1\n";
#endif

    if (items_in_set1 < min_items_in_set1 || items_in_set1 > max_items_in_set1)
      break;
  }

  if (verbose > 1)
    cerr << "Random refinement reduced cost from " << initial_cost << " to " << _cost << ", " << items_in_set1 << " items in set1, " << compute_percent (items_in_set1, pool_size) << "%\n";

  return 1;
}

int
Set_Membership::iterative_refinement_deterministic (int nsteps, 
                                                    similarity_type_t threshold)
{
  int * cost = new_int (pool_size);

  int * already_done = new_int (pool_size);

  int rc = _iterative_refinement_deterministic (nsteps, threshold, cost, already_done);

  delete[] cost;
  delete[] already_done;

  return rc;
}

int
Set_Membership::iterative_refinement_random (int nsteps, 
                                             similarity_type_t threshold)
{
  int * cost = new_int (pool_size);
  int * already_done = new_int (pool_size);
  int * tmp = new int[pool_size];

  int rc = _iterative_refinement_random (nsteps, threshold, cost, already_done, tmp);

  delete [] cost;
  delete [] already_done;
  delete [] tmp;

  return rc;
}

/*
  We sort on closenest to the desired number of items in set 1
*/

static int
Set_Membership_cost_comparitor (const void * psm1, const void * psm2)
{
  const Set_Membership * s1 = (const Set_Membership *) psm1;
  const Set_Membership * s2 = (const Set_Membership *) psm2;

  int c1 = s1->cost ();
  int c2 = s2->cost ();

  if (c1 < c2)
    return -1;

  if (c2 < c1)
    return 1;

  return 0;
}

static int
write_histogram (similarity_type_t threshold,
                 const Set_Membership & sm, 
                 int ndx, 
                 const IWString & histogram_name_stem)
{
  IWString fname (histogram_name_stem);
  fname << ndx;

  return sm.write_close_distance_histogram (threshold, fname);
}

static int
report_close_distances (similarity_type_t threshold,
                        std::ostream & output)
{
  Accumulator<similarity_type_t> close_distances;

  int number_within_threshold = 0;
  int close_distances_violating_threshold = 0;

  int total_cost = 0;

  for (int i = 0; i < pool_size; i++)
  {
    Test_Set_Member & pi = pool[i];

    int s = pi.set_membership ();

    int other_class;
    if (1 == s)
      other_class = 0;
    else
      other_class = 1;

    similarity_type_t closest_distance;
    int c = pi.distances_to_class (other_class, threshold, closest_distance);

    if (c)    // there were neighbours in the other class
    {
      number_within_threshold += c;
      close_distances.extra (closest_distance);
      if (closest_distance < threshold)
        close_distances_violating_threshold++;

      total_cost += pi.cost_in_class (s, threshold);
    }
  }

  output << number_within_threshold << " inter-class pairs within " << threshold << endl;
  output << "Inter set pair closest distances between " << close_distances.minval () << " and " << close_distances.maxval () << " ave " << close_distances.average () << endl;
  output << "Total cost " << total_cost << endl;

  return output.good ();
}

static int
report_costs (Set_Membership * sm,
              int nrandom,
              std::ostream & output)
{
  Accumulator_Int<int> acc;
  for (int i = 0; i < nrandom; i++)
  {
    Set_Membership & smi = sm[i];

    cerr << " split " << i << " cost " << smi.cost () << ", " << smi.items_in_set1 () << " items in set 1, " << compute_percent (smi.items_in_set1 (), pool_size) << "%\n";

    acc.extra (smi.cost ());
  }

  output << nrandom << " splits, costs between " << acc.minval () << " and " << acc.maxval ();
  if (acc.n () > 0)
    output << ", ave " << acc.average ();
  output << endl;

  return output.good ();
}

/*
  Do iterative refinement on a number of replicates. Return the index of one with a lower
*/

static int
gfp_e_val_divide (similarity_type_t threshold,
                  Set_Membership * rsm,
                  int number_copies)
{
  int initial_cost = rsm[0].cost ();

  for (int i = 0; i < number_copies; i++)
  {
    rsm[i].iterative_refinement_random (iterative_refinement_steps, threshold);
  }

  qsort (rsm, number_copies, sizeof (Set_Membership), Set_Membership_cost_comparitor);

  int number_better = 0;
  for (int i = 0; i < number_copies; i++)
  {
    if (rsm[i].cost () <= initial_cost)
      number_better++;
    else
      break;
  }

  if (0 == number_better)
    return -1;

  static std::random_device rd;
  static std::mt19937 rng(rd());
  std::uniform_int_distribution<int> u(0, number_better / 2);
  return u(rng);
}

/*
  See if we can optimise a given set membership. We make a bunch of copies of our initial set membership,
*/

static int
gfp_e_val_divide (similarity_type_t threshold,
                  Set_Membership & sm)
{
  sm.iterative_refinement_deterministic (iterative_refinement_steps * 20, threshold);

  Set_Membership * rsm = new Set_Membership[iterative_refinement_replicates];

  for (int i = 0; i < iterative_refinement_replicates; i++)
  {
    rsm[i] = sm;
  }

  int better = gfp_e_val_divide (threshold, rsm, iterative_refinement_replicates);

  if (better >= 0)      // great, we got a better splitting
    sm = rsm[better];

  delete [] rsm;

  return 1;
}


static int
gfp_e_val_divide (similarity_type_t threshold,
                  Set_Membership * sm,
                  int nrandom,
                  int nsplits,
                  const IWString & test_stem,
                  const IWString & train_stem)
{
  for (int i = 0; i < nrandom; i++)
  {
    Set_Membership & smi = sm[i];

    if (verbose > 1)
      cerr << "Start random replicate " << i << endl;

    gfp_e_val_divide (threshold, smi);
  }

  qsort (sm, nrandom, sizeof (Set_Membership), Set_Membership_cost_comparitor);

  for (int i = 0; i < nsplits; i++)
  {
    IWString fname (test_stem);
    fname << i;

    std::ofstream test_output (fname.null_terminated_chars (), std::ios::out);

    if (! test_output.good ())
    {
      cerr << "Cannot open output test set file '" << fname << "'\n";
      return 0;
    }

    fname = train_stem;
    fname << i;

    std::ofstream train_output (fname.null_terminated_chars (), std::ios::out);

    if (! train_output.good ())
    {
      cerr << "Cannot open output training set file '" << fname << "'\n";
      return 0;
    }

    sm[i].adjust_set_memberships_to_desired_numbers (threshold);

    if (! sm[i].write (threshold, test_output, train_output))
      return 0;

    if (histogram_name_stem.length ())
      write_histogram (threshold, sm[i], i, histogram_name_stem);
  }

  return 1;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << " -t <dist>     similarity threshold\n";
  cerr << " -n <nsplits>  number of splits to make\n";
  cerr << " -r <number>   number of random splits to start with, must be >= -n option\n";
  cerr << " -I <number>   number of iterative refinement replicates\n";
  cerr << " -i <number>   number of iterative refinement steps\n";
  cerr << " -s <number>   number of items in input file to process (normally not needed)\n";
  cerr << " -R <stem>     file name stem for tRaining set molecules\n";
  cerr << " -E <stem>     file name stem for tEst     set molecules\n";
  cerr << " -f <fraction> fraction of the input that are to be put in set 1\n";
  cerr << " -v            verbose output\n";


  exit (rc);
}

static int
gfp_e_val_divide (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vt:s:n:E:R:I:i:r:f:H:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present ('t'))
  {
    cerr << "Must specify initial threshold via the -t option\n";
    usage (2);
  }

  if (cl.option_present ('t'))
  {
    if (! cl.value ('t', initial_threshold) || initial_threshold <= 0.0 || initial_threshold >= 1.0)
    {
      cerr << "The initial threshold value (-T option) must be a valid distance\n";
      usage (3);
    }

    if (verbose)
      cerr << "Initial threshold set to " << initial_threshold << endl;
  }

  if (cl.option_present ('f'))
  {
    if (! cl.value ('f', fraction_in_set_1) || fraction_in_set_1 <= 0.0 || fraction_in_set_1 >= 1.0)
    {
      cerr << "The fraction in set 1 option (-f) must be followed by a valid fraction\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will include " << fraction_in_set_1 << " of the input in set 1\n";
  }

  if (cl.option_present ('I'))
  {
    if (! cl.value ('I', iterative_refinement_replicates) || iterative_refinement_replicates < 1)
    {
      cerr << "Invalid value for iterative refinement replicates (-I option)\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will make " << iterative_refinement_replicates << " copies of each initial distribution for iterative refinement\n";
  }

  if (cl.option_present ('i'))
  {
    if (! cl.value ('i', iterative_refinement_steps) || iterative_refinement_steps < 1)
    {
      cerr << "Invalid value for iterative refinement nsteps (i option)\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will do " << iterative_refinement_steps << " steps of iterative refinement\n";
  }

// Some defaults so they can specify just the -I or -i option

  if (iterative_refinement_steps && 0 == iterative_refinement_replicates)
    iterative_refinement_replicates = 10;

  if (0 == iterative_refinement_steps && iterative_refinement_replicates)
    iterative_refinement_steps = iterative_refinement_replicates * 10;

  if (cl.option_present ('s'))
  {
    int nfp;
    if (! cl.value ('s', nfp) || nfp < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }

    pool = new Test_Set_Member[nfp];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << nfp << endl;
      return 62;
    }

    pool_size = nfp;

    if (verbose)
      cerr << "System sized to " << nfp << endl;
  }

  int nsplits = 1;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', nsplits) || nsplits < 1)
    {
      cerr << "Invalid number of splits (-n) option\n";
      usage (17);
    }

    if (verbose)
      cerr << "Will produce " << nsplits << " splits of the data\n";
  }

  if (cl.option_present ('r'))
  {
    if (! cl.value ('r', nrandom) || nrandom < 1)
    {
      cerr << "The number of random splits option (-r) must be followed by a whole positive number\n";
      usage (18);
    }

    if (nrandom < nsplits)
    {
      cerr << "You have asked for " << nsplits << " splits of the data, but only " << nrandom << " random starting points\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will do " << nrandom << " random splits of the data\n";
  }
  else
  {
    nrandom = nsplits * 10;
  }

  if (! cl.option_present ('E') || ! cl.option_present ('R'))
  {
    cerr << "Must specify test (-E) and training (-R) file names\n";
    usage (14);
  }

  if (cl.option_present ('H'))
  {
    histogram_name_stem = cl.string_value ('H');

    if (verbose)
      cerr << "Will write inter-class closest distance histogram(s) to '" << histogram_name_stem << "'\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (1 != cl.number_elements ())
  {
    cerr << "Sorry, only processes one file at a time\n";
    return 8;
  }

  if (! build_pool (cl[0]))
  {
    cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
    return 0;
  }

  effective_pool_size = pool_size;

  if (! establish_parent_pointers ())
  {
    cerr << "Cannot establish parent pointers\n";
    return 0;
  }

  items_needed_in_set_1 = compute_items_needed_in_set1 ();

  if (verbose)
    cerr << "Ideal split will have " << items_needed_in_set_1 << " items in set 1\n";

  qsort (pool, pool_size, sizeof (Test_Set_Member), Test_Set_Member_comparitor_number_neighbours);

  effective_pool_size = pool_size;
  for (int i = 0; i < pool_size; i++)
  {
    if (0 == pool[i].number_neighbours ())
    {
      effective_pool_size = i + 1;
      break;
    }
  }

  if (verbose)
    cerr << "only " << effective_pool_size << " of " << pool_size << " molecules have neighbours\n";

  qsort (pool, effective_pool_size, sizeof (Test_Set_Member), Test_Set_Member_comparitor_closest_distance);

  if (verbose > 2)
  {
    for (int i = 0; i < effective_pool_size; i++)
    {
      cerr << i << " '" << pool[i].id () << "' dist " << pool[i].closest_distance () << endl;
    }
  }

  (void) determine_molecules_with_neighbours_within_threshold (initial_threshold);

  int distances_less_than_threshold = 0;

  Accumulator<similarity_type_t> closest_distance;

  for (int i = 0; i < effective_pool_size; i++)
  {
    int dw = pool[i].distances_within (initial_threshold);

    distances_less_than_threshold += dw;

    if (dw)
      closest_distance.extra (pool[i].closest_distance ());
  }

  if (verbose)
  {
    cerr << "Without regard to class membership, " << distances_less_than_threshold << " distances less than " << initial_threshold << endl;
    cerr << "Nearest distances between " << closest_distance.minval () << " and " << closest_distance.maxval ();
    if (closest_distance.n () > 1)
      cerr << " ave " << closest_distance.average ();
    cerr << endl;
  }

  if (0 == distances_less_than_threshold)
  {
    cerr << "Huh, no distances less than " << initial_threshold << endl;
    return 3;
  }

  assert (nrandom > 0);

  Set_Membership * sm = new Set_Membership[nrandom];
  for (int i = 0; i < nrandom; i++)
  {
    sm[i].initialise_randomly ();

    sm[i].determine_cost (initial_threshold);
  }

  if (verbose)
    report_close_distances (initial_threshold, cerr);

  qsort (sm, nrandom, sizeof (Set_Membership), Set_Membership_cost_comparitor);

  if (verbose > 1)
  {
    cerr << "Random splits\n";
    report_costs (sm, nrandom, cerr);
  }

  if (verbose > 1)
    cerr << "Beginning initial neighbour groupings\n";

  for (int i = 0; i < nrandom; i++)
  {
    sm[i].initial_neighbour_groupings (initial_threshold);
  }

  qsort (sm, nrandom, sizeof (Set_Membership), Set_Membership_cost_comparitor);

  if (verbose > 1)
  {
    cerr << "After initial groupings\n";
    report_costs (sm, nrandom, cerr);
  }

  if (verbose)
    report_close_distances (initial_threshold, cerr);

  IWString test_fname  =  cl.string_value ('E');
  IWString train_fname =  cl.string_value ('R');
   
  if (0 == iterative_refinement_replicates)
    sm[0].write (initial_threshold, test_fname, train_fname);
  else
    gfp_e_val_divide (initial_threshold, sm, nrandom, nsplits, test_fname, train_fname);

  delete [] sm;

  if (NULL != pool)
    delete [] pool;

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = gfp_e_val_divide (argc, argv);

  return rc;
}
