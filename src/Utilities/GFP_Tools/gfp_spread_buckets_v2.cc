/*
  Spread algorithm with bucket restrictions.

  We have different types of buckets

  $SMI<C>
  BUCKET1<0>
  BUCKET2 <89>
  BUCKET3<4>
  |
  $SMI<CC>
  BUCKET1<7>
  BUCKET2 <8>
  BUCKET3<48>
  |

  We enforce a uniform selection from each bucket type, and within each bucket
  type, a uniform distribution from bucket numbers
*/

#include <stdlib.h>
#include <limits>

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwstring_data_source.h"
#include "iw_tdt.h"
#include "misc.h"
#include "cmdline.h"
#include "accumulator.h"
#include "iwrandom.h"

#include "spread_v2.h"

static int verbose = 0;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");
static IWString bucket_tag ("BCKT<");

static int strict_rotation = 0;

/*
  It is optional as to whether or not the -A file has bucket information
  in it.
*/

static int previously_selected_file_has_bucket_information = 0;

static int convert_bucket_column_data_from_float = 0;

/*
  Each bucket needs to keep track of the number of items which are
  still in it, and the number which have been selected from it
*/

class Bucket
{
  private:
    int _number_unselected;

    int _number_selected;

    int _can_be_selected_from;

  public:
    Bucket ();

    int   number_unselected () const { return _number_unselected;}
    void  another_unselected_item () { _number_unselected++;}

    int   number_selected ()   const { return _number_selected;}

    int   another_item_selected ();

//  During processing the previously selected file, we need to tell buckets
//  that items have been selected from them

    void  increment_number_selected () { _number_selected++;}

//  During selection, people need to know whether we can choose more items
//  from a given bucket or not

    int can_select_from () const { return _can_be_selected_from;}
    void set_can_select_from (int s) { _can_be_selected_from = s;}
};

Bucket::Bucket ()
{
  _number_selected = 0;
  _number_unselected = 0;

  return;
}

int
Bucket::another_item_selected ()
{
  _number_unselected--;
  _number_selected++;

  assert (_number_unselected >= 0);

  return 1;
}

/*
  We can have any number of bucket types;
*/

int number_bucket_types = 0;

/*
  Each bucket type needs info about its members
*/

class Bucket_Type
{
  private:

//  Two ways of specifying the bucket, a tag or a column

    IWString _tag;

    int _bucket_column_in_tag;

    IWString _name;

//  We don't require every fingerprint to have a value for every bucket

    int _items_with_tag_specified;

    int _lowest_occupied_bucket;

    int _highest_occupied_bucket;

    Bucket * _bucket;

//  private functions

    int _parse_column_in_tag (const IWString & t);

    int _discern_bucket_from_tag (const IW_TDT & tdt, int & b) const;
    int _discern_bucket_from_column_in_name (const IW_TDT & tdt, int & b) const;

  public:
    Bucket_Type ();
    ~Bucket_Type ();

    const IWString & tag () const { return _tag;}
    const IWString & name () const { return _name;}

    int   set_tag (const IWString & t);

    int discern_bucket (const IW_TDT & tdt, int & b) const;

    int status (std::ostream &) const;

    int lowest_occupied_bucket () const { return _lowest_occupied_bucket;}
    int highest_occupied_bucket () const { return _highest_occupied_bucket;}

    void bucket_is_occupied (int);

    int initialise_bucket_arrays ();

    int determine_buckets_from_which_we_can_select ();
    int can_items_be_selected_from_bucket (int b) const;

    void item_in_bucket (int);

    void another_item_selected (int);

    void increment_number_selected (int b) { _bucket[b].increment_number_selected ();}
};

Bucket_Type::Bucket_Type () : _tag (bucket_tag)
{
  _items_with_tag_specified = 0;
  _lowest_occupied_bucket = -1;
  _highest_occupied_bucket = -1;

  _bucket = NULL;

  _bucket_column_in_tag = -1;
}

Bucket_Type::~Bucket_Type ()
{
  if (-85 == _lowest_occupied_bucket)
    cerr << "Deleting already freed Bucket_Type\n";

  _lowest_occupied_bucket = -85;

  if (NULL != _bucket)
    delete [] _bucket;

  return;
}

int 
Bucket_Type::status (std::ostream & os) const
{
  os << "Bucket type '" << _name << "', lowest occupied bucket " << _lowest_occupied_bucket << " highest occupied bucket " << _highest_occupied_bucket << endl;

  if (NULL != _bucket)
  {
    for (int i = _lowest_occupied_bucket; i <= _highest_occupied_bucket; i++)
    {
      os << " bucket " << i << " has " << _bucket[i].number_unselected () << " unselected items";
      if (_bucket[i].number_selected ())
        os << ", " << _bucket[i].number_selected () << " selected";
      os << endl;
    }
  }

  return os.good ();
}

int
Bucket_Type::_parse_column_in_tag (const IWString & t)
{
  int i = t.index ("col=");

  if (i < 0)
    return 0;

  const_IWSubstring myt (t);

  if (i > 0)
  {
    _tag.set (myt.rawchars (), i - 1);

    _name = _tag;

    if (! _tag.ends_with ('<'))
      _tag += '<';

    myt.remove_leading_chars (i);
  }

  assert (myt.starts_with ("col="));

  myt.remove_leading_chars (4);

  if (! myt.numeric_value (_bucket_column_in_tag) || _bucket_column_in_tag < 1)
  {
    cerr << "Bucket_Type::_parse_column_in_tag:invalid column '" << t << "'\n";
    return 0;
  }

  _name << ".COL" << _bucket_column_in_tag << '<';

  _bucket_column_in_tag--;

  return 1;
}

int
Bucket_Type::set_tag (const IWString & t)
{
  if (t.contains ("col="))
  {
    if (! _parse_column_in_tag (t))
      return 0;
  }
  else
  {
    _tag = t;

    if (! _tag.ends_with ('<'))
      _tag += '<';

    _name = t;
  }

  return 1;
}

void
Bucket_Type::bucket_is_occupied (int b)
{
  _items_with_tag_specified++;

  if (-1 == _lowest_occupied_bucket)
  {
    _lowest_occupied_bucket = b;
    _highest_occupied_bucket = b;

    return;
  }

  if (b < _lowest_occupied_bucket)
    _lowest_occupied_bucket = b;
  else if (b > _highest_occupied_bucket)
    _highest_occupied_bucket = b;

  return;
}

int
Bucket_Type::initialise_bucket_arrays ()
{
  if (-1 == _lowest_occupied_bucket)
  {
    cerr << "Bucket_Type::initialise_bucket_arrays: no occupied buckets for '" << _tag << "'\n";
    return 0;
  }

  if (_lowest_occupied_bucket == _highest_occupied_bucket)
  {
    cerr << "Bucket_Type::initialise_bucket_arrays: all buckets the same " << _lowest_occupied_bucket << " for '" << _tag << "'\n";
    return 0;
  }

  _bucket = new Bucket[_highest_occupied_bucket + 1];

  return 1;
}

void
Bucket_Type::item_in_bucket (int b)
{
  assert (b >= _lowest_occupied_bucket && b <= _highest_occupied_bucket);

  _bucket[b].another_unselected_item ();

  return;
}

int
Bucket_Type::determine_buckets_from_which_we_can_select ()
{
  int min_selected = std::numeric_limits<int>::max();

  for (int i = _lowest_occupied_bucket; i <= _highest_occupied_bucket; i++)
  {
    const Bucket & b = _bucket[i];

    if (0 == b.number_unselected ())
      continue;

    if (b.number_selected () < min_selected)
      min_selected = b.number_selected ();
  }

  if (std::numeric_limits<int>::max() == min_selected)    // all our buckets exhausted
    return 0;

  if (verbose > 2)
    cerr << "Min selected from bucket " << min_selected << endl;

// All buckets with the smallest number of items selected can now be selected from

  for (int i = _lowest_occupied_bucket; i <= _highest_occupied_bucket; i++)
  {
    Bucket & b = _bucket[i];

    if (b.number_selected () == min_selected)
      b.set_can_select_from (1);
    else
      b.set_can_select_from (0);
  }

  return 1;
}

int
Bucket_Type::can_items_be_selected_from_bucket (int b) const
{
  assert (b >= _lowest_occupied_bucket && b <= _highest_occupied_bucket);

  return _bucket[b].can_select_from ();
}

void
Bucket_Type::another_item_selected (int b)
{
  assert (b >= _lowest_occupied_bucket && b <= _highest_occupied_bucket);

  _bucket[b].another_item_selected ();

  return;
}

int
Bucket_Type::_discern_bucket_from_column_in_name (const IW_TDT & tdt,
                                int & b) const
{
  const_IWSubstring pcn;
  if (! tdt.dataitem_value(identifier_tag, pcn))
  {
    cerr << "Bucket_Type::_discern_bucket_from_column_in_name:cannot extract '" << identifier_tag << "\n";
    return 0;
  }

  if (pcn.nwords() <= _bucket_column_in_tag)
  {
    cerr << "Bucket_Type::_discern_bucket_from_column_in_name:not enough tokens in name '" << pcn << "'\n";
    return 0;
  }

  const_IWSubstring string_bucket;

  pcn.word(_bucket_column_in_tag, string_bucket);

  if (convert_bucket_column_data_from_float)
  {
    float f;
    if (! string_bucket.numeric_value(f))
    {
      cerr << "Bucket_Type::_discern_bucket_from_column_in_name:invalid float '" << pcn << "'\n";
      return 0;
    }

    if (f < 0.0f)
      b = 0;
    else
      b = static_cast<int>(f + 0.4999f);
  }
  else if (! string_bucket.numeric_value(b) || b < 0)
  {
    cerr << "Bucket_Type::_discern_bucket_from_column_in_name:invalid bucket '" << string_bucket << "'\n";
    return 0;
  }

  return 1;
}

int
Bucket_Type::_discern_bucket_from_tag (const IW_TDT & tdt,
                                       int & b) const
{
  const_IWSubstring string_bucket;
  if (! tdt.dataitem_value(_tag, string_bucket))
  {
    cerr << "Bucket_Type::_discern_bucket_from_tag:no '" << _tag << "' in TDT\n";
    return 0;
  }

  if (! string_bucket.numeric_value(b) || b < 0)
  {
    cerr << "Bucket_Type::_discern_bucket_from_tag:invalid bucket '" << string_bucket << "'\n";
    return 0;
  }

  return 1;
}

int
Bucket_Type::discern_bucket (const IW_TDT & tdt,
                             int & b) const
{
  const_IWSubstring string_bucket;

  if (_bucket_column_in_tag >= 0)
    return _discern_bucket_from_column_in_name(tdt, b);
  else
    return _discern_bucket_from_tag (tdt, b);
}

static Bucket_Type * bucket_type = NULL;

class Bucket_Spread_Object: public Spread_Object
{
  private:
    int * _bucket;

// When selecting from multiple buckets, we could have deadlock, so we
// prioritise the unselected items by the number of times their
// selection is blocked by one of their bucket values

    int _number_times_blocked;

  public:
    Bucket_Spread_Object ();
    ~Bucket_Spread_Object ();

    int construct_from_tdt (IW_TDT &, int &);

    int bucket (int b) const { return _bucket[b];}

    int number_times_blocked () const { return _number_times_blocked;}
    int determine_number_times_blocked ();

    int can_be_selected_in_strict_rotation () const;
};

#define BUCKET_VALUE_NOT_SET -938383

Bucket_Spread_Object::Bucket_Spread_Object ()
{
  assert (number_bucket_types > 0);

  _bucket = new_int (number_bucket_types, BUCKET_VALUE_NOT_SET);

  _number_times_blocked = 0;
}

Bucket_Spread_Object::~Bucket_Spread_Object ()
{
  delete _bucket;
}

int
Bucket_Spread_Object::construct_from_tdt (IW_TDT & tdt, int & fatal)
{
  if (! Spread_Object::construct_from_tdt (tdt, fatal))
    return 0;

  for (int i = 0; i < number_bucket_types; i++)
  {
    if (! bucket_type[i].discern_bucket (tdt, _bucket[i]))
    {
      cerr << "Bucket_Spread_Object::construct_from_tdt:missing or invalid bucket in TDT\n";
      cerr << tdt;
      fatal = 1;
      return 0;
    }
  }

  return 1;
}

int
Bucket_Spread_Object::determine_number_times_blocked ()
{
  _number_times_blocked = 0;

  for (int i = 0; i < number_bucket_types; i++)
  {
    if (! bucket_type[i].can_items_be_selected_from_bucket (_bucket[i]))
      _number_times_blocked++;
  }

  return _number_times_blocked;
}

int
Bucket_Spread_Object::can_be_selected_in_strict_rotation () const
{
  for (int i = 0; i < number_bucket_types; i++)
  {
    if (! bucket_type[i].can_items_be_selected_from_bucket (_bucket[i]))
      return 0;
  }

  return 1;
}

/*
  When doing a run with no pre-selected molecules, start with the
  object which is furthest from the first fingerprint.
*/

static int start_with_object_furthest_from_first = 0;

static int choose_first_item_randomly = 0;

static int already_selected_molecules_present = 0;

/*
  Our pool is an array of FP objects
*/

static Bucket_Spread_Object * pool = NULL;

static int pool_size = 0;

static int poolptr = 0;    // the next item in the pool to be filled

static int number_to_select = 0;

static int report_establish_initial_distances = 0;

/*
  We keep track of the distances of the nearest selected items
*/

static Accumulator<similarity_type_t> nearest_selected_neighbour_distance;

/*
  A two-pass process. First we tell the Bucket_Type objects what buckets
  are occupied, so they can initialise their _lowest_occupied_bucket and
  _highest_occupied_bucket values. We then have them allocate their
  arrays based on those maxima and minima, and the second pass fills them
  in.
*/

static int
initialise_bucket_structures ()
{
  if (verbose)
    cerr << "Initialising spread object, pool size = " << pool_size << endl;

  for (int i = 0; i < pool_size; i++)
  {
    for (int j = 0; j < number_bucket_types; j++)
    {
      int b = pool[i].bucket (j);

      bucket_type[j].bucket_is_occupied (b);
    }
  }

  int rc = number_bucket_types;

  for (int i = 0; i < number_bucket_types; i++)
  {
    Bucket_Type & btype = bucket_type[i];

    if (! btype.initialise_bucket_arrays ())
    {
      cerr << "Cannot allocate bucket structures for bucket type '" << btype.tag () << "'\n";
      rc = 0;
    }
  }

  if (0 == rc)
    return 0;

  for (int i = 0; i < pool_size; i++)
  {
    for (int j = 0; j < number_bucket_types; j++)
    {
      int b = pool[i].bucket (j);

      bucket_type[j].item_in_bucket (b);
    }
  }

  if (verbose)
  {
    for (int i = 0; i < number_bucket_types; i++)
    {
      bucket_type[i].status (cerr);
    }
  }

  return 1;
}

/*
  We are just about to do a selection, ask each bucket type to figure out
  which of its buckets are available for selection
*/

static void
determine_which_buckets_are_avaiable_for_selections ()
{
  for (int i = 0; i < number_bucket_types; i++)
  {
    bucket_type[i].determine_buckets_from_which_we_can_select ();
  }

  return;
}

static int
build_pool (iwstring_data_source & input)
{
  assert (pool_size > 0);
  assert (poolptr >= 0 && poolptr < pool_size);

  off_t offset = input.tellg ();

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    if (! pool[poolptr].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      offset = input.tellg ();
      continue;
    }

    pool[poolptr].set_offset (offset);

    poolptr++;

    if (poolptr >= pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }

    offset = input.tellg ();
  }

  poolptr--;

  if (verbose)
    cerr << "Pool now contains " << (poolptr + 1) << " objects\n";

  return 1;
}

static int
allocate_pool ()
{
  assert (pool_size > 0);
  assert (NULL == pool);

  pool = new Bucket_Spread_Object[pool_size];
  if (NULL == pool)
  {
    cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
    return 62;
  }

  if (verbose)
    cerr << "system sized to " << pool_size << endl;

  return 1;
}

static int
build_pool (const char * fname)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    IW_Regular_Expression pcn_rx ("^PCN<");   // should use identifier-tag

    pool_size = input.grep (pcn_rx);

    if (0 == pool_size)
    {
      cerr << "Zero occurrences of '" << pcn_rx.source () << "' in '" << fname << "'\n";
      return 0;
    }

    if (! allocate_pool ())
      return 0;
  }

  return build_pool (input);
}

template <typename T>
void
compare_against_whole_pool (T & fp,
                            int ntdt)
{
  for (int i = 0; i < pool_size; i++)
  {
    if (can_be_compared (pool[i], fp))
      pool[i].object_has_been_selected (fp);
  }

  if (report_establish_initial_distances && 0 == ntdt % report_establish_initial_distances)
    cerr << "Established initial distances " << ntdt << endl;

  return;
}

template void compare_against_whole_pool (Spread_Object &, int);
template void compare_against_whole_pool (Bucket_Spread_Object &, int);

static int
establish_initial_distances_previous_bucket_information (iwstring_data_source & input)
{
  int ntdt = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    Bucket_Spread_Object fp;
    if (! fp.construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot build fingerpint\n" << tdt << endl;
        return 0;
      }

      continue;
    }

    ntdt++;

    compare_against_whole_pool (fp, ntdt);

    for (int i = 0; i < number_bucket_types; i++)
    {
      int b = fp.bucket (i);

      bucket_type[i].increment_number_selected (b);
    }
  }

  return 1;
}

static int
establish_initial_distances_no_previous_bucket_information (iwstring_data_source & input)
{
  int ntdt = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    int fatal;
    Spread_Object fp;
    if (! fp.construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot build fingerpint\n" << tdt << endl;
        return 0;
      }

      continue;
    }

    ntdt++;

    compare_against_whole_pool (fp, ntdt);
  }

  return 1;
}

static int
establish_initial_distances (iwstring_data_source & input)
{
  int rc;

  if (previously_selected_file_has_bucket_information)
    rc = establish_initial_distances_previous_bucket_information (input);
  else
    rc = establish_initial_distances_no_previous_bucket_information (input);

  if (0 == rc)
    return 0;

  return 1;
}

static int
establish_initial_distances (const const_IWSubstring & fname)
{
  for (int i = 0; i < pool_size; i++)
  {
    pool[i].set_distance (1000.0);
  }

  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open already selected file '" << fname << "'\n";
    return 0;
  }

  return establish_initial_distances (input);
}

static void
do_object_has_been_selected (int isel,
                             int number_selected,
                             IWString_and_File_Descriptor & output)
{
  Bucket_Spread_Object & fpsel = pool[isel];

  const Smiles_ID_Dist & nsn = fpsel.nsn ();

  if (verbose)
  {
    cerr << "Selected " << number_selected << " '" << fpsel.id () << "' distance " << fpsel.distance () << ", bucket ";
    for (int i = 0; i < number_bucket_types; i++)
    {
      if (i > 0)
        cerr << ',';

      cerr << fpsel.bucket (i);
    }

    cerr << " Blocked " << fpsel.number_times_blocked () << endl;

    if (verbose)
      nearest_selected_neighbour_distance.extra (fpsel.distance ());
  }

  output << smiles_tag << fpsel.smiles () << ">\n";
  output << identifier_tag << fpsel.id () << ">\n";
  for (int i = 0; i < number_bucket_types; i++)
  {
    const IWString & btag = bucket_type[i].tag ();

    output << btag << fpsel.bucket (i) << ">\n";
  }

  if (fpsel.has_a_nearest_selected_neighbour ())
  {
    output << smiles_tag << nsn.smiles () << ">\n";
    output << identifier_tag << nsn.id () << ">\n";
    output << distance_tag << fpsel.distance () << ">\n";
  }
  output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  fpsel.set_selected ();

  for (int i = 0; i < number_bucket_types; i++)
  {
    int b = fpsel.bucket (i);

    bucket_type[i].another_item_selected (b);
  }

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i].selected () || i == isel)
      continue;

    if (! can_be_compared (pool[i], fpsel))
      continue;

    pool[i].object_has_been_selected (fpsel);
  }

  return;
}

static int
do_start_with_object_furthest_from_first (int & istart)
{
  resizable_array<int> already_done;
  already_done.resize (start_with_object_furthest_from_first);

  for (int i = 0; i < start_with_object_furthest_from_first; i++)
  {
    already_done.add (istart);

    Bucket_Spread_Object & fp0 = pool[istart];

    int furthest_away = -1;

    similarity_type_t d0 = 0.0;
    for (int j = 1; j < pool_size; j++)
    {
      similarity_type_t d = pool[j].IW_General_Fingerprint::distance (fp0);
      if (d <= d0)
        continue;

      if (j == istart || already_done.contains (j))
        continue;

      d0 = d;
      furthest_away = j;
    }

    assert (furthest_away > 0);

    if (verbose)
      cerr << "Furthest from first fingerprint is " << furthest_away << " '" << pool[istart].id () << "', distance " << d0 << endl;

    istart = furthest_away;
  }

  return 1;
}

/*
  We have a linear penalty against blocked buckets
*/

static similarity_type_t bucket_violation_factor = 0.5;

static similarity_type_t
distance_penalised_for_blockages (similarity_type_t d,
                                  int number_bucket_types,
                                  int number_times_blocked)
{
  similarity_type_t scale = bucket_violation_factor +
                                (1.0 - bucket_violation_factor) / static_cast<similarity_type_t> (number_bucket_types) *
                        static_cast<similarity_type_t> (number_bucket_types - number_times_blocked);

  return scale * d;
}

static int
spread_buckets (IWString_and_File_Descriptor & output)
{
  int first_selected;

  if (choose_first_item_randomly)
    first_selected = intbtwij (0, pool_size - 1);
  else
    first_selected = 0;

  if (start_with_object_furthest_from_first)
    do_start_with_object_furthest_from_first (first_selected);

//#define CHECK_PENALTY_FUNCTION
#ifdef CHECK_PENALTY_FUNCTION
    for (int i = 0; i <= number_bucket_types; i++)
    {
      for (similarity_type_t d = 0.0; d < 1.0; d+= 0.2)
      {
        similarity_type_t w = distance_penalised_for_blockages (d, number_bucket_types, i);
        cerr << "Distance " << d << " blocked " << i << " times: " << w << endl;
      }
    }
#endif

  int number_selected = 1;

  do_object_has_been_selected (first_selected, number_selected, output);

  while (number_selected < number_to_select)
  {
#ifdef USE_IWMALLOC
    iwmalloc_check_all_malloced (stderr);
#endif

    determine_which_buckets_are_avaiable_for_selections ();

    similarity_type_t maxdist = -1.0;
    int ichoose = -1;

    for (int i = 0; i < pool_size; i++)
    {
      Bucket_Spread_Object & p = pool[i];

      if (p.selected ())
        continue;

      if (strict_rotation && ! p.can_be_selected_in_strict_rotation())
        continue;

      similarity_type_t d = p.distance ();

      int number_times_blocked = p.determine_number_times_blocked ();

      if (number_times_blocked > 0)
        d = distance_penalised_for_blockages (d, number_bucket_types, number_times_blocked);

      if (d > maxdist)
      {
        maxdist = d;
        ichoose = i;
      }
    }

    assert (ichoose >= 0);

    number_selected++;

    do_object_has_been_selected (ichoose, number_selected, output);
  }

  return number_selected;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Usage <options> <input_file>\n";
  cerr << " -B <tag>         tag containing bucket information\n";
  cerr << "                  buckets must be integer values\n";
  cerr << " -B col=<col>     bucket value in column <col> of name field\n";
  cerr << " -f               convert column values from float to int to get bucket number\n";
  cerr << " -b <scale>       penalty factor for violating bucket constraints, range (0.0,1.0)\n";
  cerr << "                  high numbers mean a larger penalty for violating bucket constraints\n";
  cerr << " -o               select equally from all buckets before\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many items to select\n";
  cerr << " -A <file>        specify file of already selected items\n";
  cerr << " -a               already selected items have bucket information\n";
  cerr << " -i <dataitem>    specify identifier dataitem in pool\n";
  cerr << " -I <dataitem>    specify identifier dataitem in input file\n";
  cerr << " -p <dataitem>    specify demerit/scale dataitem\n";
  cerr << " -p col=c         demerit/scale factor is column c of name\n";
  cerr << " -r <number>      report progress of initial distance assignments\n";
  cerr << " -C <number>      first selected is <number> times furthest from first fingerprint (DC)\n";
  cerr << " -R               first selected is random\n";
  cerr << " -F -P -Q -W      standard gfp options, enter '-F help' for info\n";
#ifdef USE_IWMALLOC
  (void) display_standard_debug_options (cerr, 'Q');
#endif
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
spread_buckets (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:I:i:A:ar:p:C:RF:P:W:B:b:w:of");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (1);
  }

  if (cl.option_present ('F') || cl.option_present ('P') || cl.option_present ('W'))
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

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }

    if (! allocate_pool ())
      return 83;
  }

  if (cl.option_present('f'))
  {
    convert_bucket_column_data_from_float = 1;

    if (verbose)
      cerr << "Columnar bucket data converted from float to nearest int\n";
  }

  number_bucket_types = cl.option_count ('B');

  if (0 == number_bucket_types)
  {
    cerr << "Must specify bucket tag(s) via the -B option\n";
    usage (7);
  }

  bucket_type = new Bucket_Type[number_bucket_types];

  if (cl.option_present ('B'))
  {
    for (int i = 0; i < number_bucket_types; i++)
    {
      Bucket_Type & btype = bucket_type[i];

      IWString btag = cl.string_value ('B', i);

      btype.set_tag (btag);

      if (verbose)
        cerr << "Bucket tag " << i << " '" << bucket_type[i].tag () << "'\n";
    }
  }

// We need to be careful with the -i and -I options. Remember
// that the pool is built first

  if (cl.option_present ('i'))
  {
    const_IWSubstring id = cl.string_value ('i');

    set_identifier_tag (id);

    if (verbose)
      cerr << "Identifiers in pool tagged as '" << id << "'\n";

    if (! identifier_tag.ends_with ('<'))
      identifier_tag += '<';
  }

  if (cl.option_present ('p'))
  {
    const_IWSubstring p = cl.string_value ('p');
    if (p.starts_with("col="))
    {
      p.remove_leading_chars(4);
      int c;
      if (! p.numeric_value(c) || c < 1)
      {
        cerr << "The scaling factor column must be a valid column number\n";
        usage(3);
      }
      set_scaling_factor_column(c);
      if (verbose)
        cerr << "Scaling factors in column " << c << " of each name\n";
    }
    else
    {
      (void) set_scale_tag (p);

      if (verbose)
        cerr << "The scale factor will be the '" << p << "' dataitem\n";
    }
  }

  if (cl.option_present ('b'))
  {
    if (! cl.value ('b', bucket_violation_factor) || bucket_violation_factor <= 0.0 || bucket_violation_factor >= 1.0)
    {
      cerr << "Invalid bucket violation specifier (-b option)\n";
      usage (6);
    }

    if (verbose)
      cerr << "Bucket violation penalty factor " << bucket_violation_factor << endl;

    bucket_violation_factor = 1.0 - bucket_violation_factor;
  }

// build the pool. Note the problems with multiple files and without the -s option.

  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! build_pool (cl[i]))
    {
      cerr << "Yipes, cannot build pool from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced (stderr);
#endif

  pool_size = poolptr + 1;

  if (! initialise_bucket_structures ())
  {
    cerr << "Yipes, could not initialise buckets\n";
    return 27;
  }

  if (cl.option_present('b') && cl.option_present('o'))
  {
    cerr << "Sorry, the -b and -o options are mutually incompatible\n";
    usage(3);
  }

  if (verbose)
  {
    const Accumulator<float> & sfs = scale_factor_statistics ();
    cerr << sfs.n () << " pool objects had demerit/scale factors\n";
    if (sfs.n ())
    {
      cerr << "Scale factors between " << sfs.minval () << " and " << sfs.maxval ();
      if (sfs.n () > 1)
        cerr << ", ave " << sfs.average ();
      cerr << endl;
    }
  }

  if (cl.option_present('o'))
  {
    strict_rotation = 1;
    if (verbose)
      cerr << "Will select structures in strict rotation\n";
  }

// Now that the pool is built, we can switch identifiers if needed

  if (cl.option_present ('I'))
  {
    const_IWSubstring id;
    cl.value ('I', id);

    set_identifier_tag (id);

    if (verbose)
      cerr << "Identifiers in input tagged as '" << id << "'\n";

    if (! identifier_tag.ends_with ('<'))
      identifier_tag += '<';
  }

  if (cl.option_present ('r'))
  {
    if (! cl.value ('r', report_establish_initial_distances) || report_establish_initial_distances < 1)
    {
      cerr << "The -r option must be followed by a whole positive number\n";
      usage (18);
    }

    if (verbose)
      cerr << "Will report initial neighbour assignments every " << report_establish_initial_distances << " fingerprints\n";
  }

  if (cl.option_present ('A'))
  {
    if (cl.option_present ('a'))
    {
      previously_selected_file_has_bucket_information = 1;
      if (verbose)
        cerr << "Previously selected items have bucket information\n";
    }

    const_IWSubstring fname;
    int i = 0;
    while (cl.value ('A', fname, i++))
    {
      if (! establish_initial_distances (fname))
      {
        cerr << "Cannot establish initial distances from '" << fname << "'\n";
        return 54;
      }
    }

    already_selected_molecules_present = 1;

    if (verbose > 1)
    {
      cerr << "INitial distances established\n";
      for (int i = 0; i < pool_size; i++)
      {
        cerr << "Pool " << i << " is " << pool[i].id () << " dist " << pool[i].distance () << endl;
      }

      for (int i = 0; i < number_bucket_types; i++)
      {
        bucket_type[i].status (cerr);
      }
    }
  }
  else
  {
    already_selected_molecules_present = 0;
  }

  if (cl.option_present ('C'))
  {
    if (already_selected_molecules_present)
    {
      cerr << "The -C option only works in the absence of the -A option\n";
      usage (18);
    }

    if (! cl.value ('C', start_with_object_furthest_from_first) || start_with_object_furthest_from_first < 1)
    {
      cerr << "The -C option must be followed by a positive whole number\n";
      usage (91);
    }

    if (verbose)
      cerr << "Will start with the molecule furthest from the first molecule in the set\n";
  }

  if (cl.option_present ('R'))
  {
    if (already_selected_molecules_present)
    {
      cerr << "The -R option only works in the absence of the -A option\n";
      usage (18);
    }

    choose_first_item_randomly = 1;
    if (verbose)
      cerr << "Will choose first fingerprint at random\n";

    iw_random_seed ();
  }

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
  else
    number_to_select = pool_size;

  IWString_and_File_Descriptor output(1);

  (void) spread_buckets (output);

  output.flush();

  if (verbose)
  {
    cerr << "Nearest previously selected item distances between " << nearest_selected_neighbour_distance.minval () << " and " << nearest_selected_neighbour_distance.maxval ();
    if (nearest_selected_neighbour_distance.n () > 1)
      cerr << " ave " << nearest_selected_neighbour_distance.average ();
    cerr << endl;

    for (int i = 0; i < number_bucket_types; i++)
    {
      bucket_type[i].status (cerr);
    }
  }

  delete [] pool;

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = spread_buckets (argc, argv);

#ifdef USE_IWMALLOC
  iwmalloc_terse_malloc_status (stderr);
#endif

  return rc;
}
