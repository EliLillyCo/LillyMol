/*
  Scans a descriptor file for similarity to a given vector
*/

#include <stdlib.h>
#include <math.h>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWARAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;

#define SET_OF_DESCRIPTORS_IMPLEMENTATION

#include "set_of_descriptors.h"
#include "iwdescriptor.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int keep_going_after_fatal_error = 0;

static int fatal_errors_encountered = 0;

static IWString smiles_tag;
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");
static IWString neighbour_tag ("NBR<");

static IWString dummy_smiles;

static IW_STL_Hash_Map_String smiles_for_haystack;
static IW_STL_Hash_Map_String smiles_for_needles;

//no of ids stored for a single compound
//equivalent of 15 FF's
static int MAX_IDs = 32768;
//static int MAX_IDs = 4096;

/*
  Some distance metrics may generate negative distances
*/

static int negative_distances_allowed = 0;

static int haystack_file_contains_header_record = 1;
static int distance_computation = DISTANCE_CARTESIAN;

static int do_not_compare_molecules_with_themselves = 0;

static int stop_when_all_molecules_have_a_zero_distance_neighbour = 0;

static int molecules_with_zero_distance_neighbours = 0;

static int strip_leading_zeros_from_identifiers = 0;

static int write_neighbours_as_index_numbers = 0;

static int haystack_records_read = 0;

static int report_progress = 0;

static int min_neighbours_to_find = -1;

static int add_number_nbrs = 0;

/*
  If we are writing neighbours as indices rather than ids', we need
  a means of knowing the index of each item
*/

static IW_STL_Hash_Map_int id_to_ndx;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Compares two descriptor files - input file is the \"haystack\"\n";
  cerr << " -p <file>        file for which neighbours are to be found 'needles'\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      number of neighbours to keep\n";
  cerr << " -m <number>      minimum number of neighbours to keep\n";
  cerr << " -T <dis>         discard distances longer than <dis>\n";
  display_distance_metrics (cerr, 'X');
  cerr << " -h               discard nbrs with zero distance and same ID as target\n";
  cerr << " -o               write neighbours as indes numbers rather than id's\n";
  cerr << " -P <file>        smiles file for the \"haystack\" molecules\n";
  cerr << " -P p:<file>      smiles file for -p (\"needle\")descriptors\n";
  cerr << " -D <smiles>      use a dummy smiles wherever needed\n";
  cerr << " -z               strip leading 0's from identifiers in matching up smiles\n";
  cerr << " -Z <number>      stop once every 'needle' has <number> zero distance neighbour(s)\n";
  cerr << " -r <number>      report progress every <number> haystack items processed\n";
  cerr << " -y               haystack file has no header, assume same column order\n";
  cerr << " -a               add the number of neighbours to each needle\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

template <typename D>
class ID_Dist_Template
{
  protected:
    IWString _id;
    D _distance;

  public:
    ID_Dist_Template ();
    ID_Dist_Template (const IWString & s, D d);
    ~ID_Dist_Template ();

    D distance () const { return _distance;}
    const IWString & id () const { return _id;}

    void set_distance (D d) { _distance = d;}
    void set_id (const IWString & s) { _id = s;}
};

template <typename D>
ID_Dist_Template<D>::ID_Dist_Template ()
{
  _distance = static_cast<D> (9.99e+09);

  return;
}

template <typename D>
ID_Dist_Template<D>::ID_Dist_Template (const IWString & s, D d) : _id (s), _distance (d)
{
}

template <typename D>
ID_Dist_Template<D>::~ID_Dist_Template ()
{
}

typedef float distance_t;

static Accumulator<distance_t> stats;

#ifdef __GNUG__
template class ID_Dist_Template<distance_t>;
#endif

//static distance_t max_distance = static_cast<distance_t> (0.0);
static distance_t max_distance = std::numeric_limits<distance_t>::max();

class ID_Distance : public ID_Dist_Template<distance_t>
{
  private:
  public:
    ID_Distance ();
    ID_Distance (const IWString &, distance_t);
};

ID_Distance::ID_Distance ()
{
}

ID_Distance::ID_Distance (const IWString & s, distance_t d) :
      ID_Dist_Template<distance_t> (s, d)
{
}

template <typename D>
class Neighbour_List
{
  private:
    resizable_array_p<ID_Distance> _neighbours;

    int _neighbours_to_find;

    int _number_zero_distance_neighbours;

//  private functions

    int  _binary_search (D) const;
    void _extra_no_max_number_neighbours (const IWString & rhs, D distance);
    void _check_insertion ();

  public:
    Neighbour_List ();
    ~Neighbour_List ();

    int set_neighbours_to_find (int);

    int number_neighbours () const { return _neighbours.number_elements ();}

    void extra (const IWString &, D);
    int  do_write (IWString_and_File_Descriptor &) const;

        //nik
        int scan_for_closer (const IWString &, D);
        //endnik

    D distance_of_closest_neighbour () const;
    D distance_of_furthest_neighbour () const;
};

template class resizable_array_p<ID_Distance>;
template class resizable_array_base<ID_Distance *>;

template <typename D>
Neighbour_List<D>::Neighbour_List ()
{
  _neighbours_to_find = 0;

  _number_zero_distance_neighbours = 0;

  return;
}

template <typename D>
Neighbour_List<D>::~Neighbour_List ()
{
  _neighbours_to_find = -2;
}

template <typename D>
int
Neighbour_List<D>::set_neighbours_to_find (int s)
{
  assert (s > 0);

  _neighbours_to_find = s;

  if (_neighbours.number_elements ())    // clean out any existing neighbour data
    _neighbours.resize (0);

  _neighbours.resize (s);

  for (int i = 0; i < _neighbours_to_find; i++)
  {
    ID_Distance * s = new ID_Distance ();
    _neighbours.add (s);
  }

  return 1;
}

static int
echo_smiles_and_id (const IWString & id,
                    const IW_STL_Hash_Map_String & smiles_for_id,
                    int write_neighbours_as_index_numbers,
                    int nbrs,
                    IWString_and_File_Descriptor & output)
{
  if (write_neighbours_as_index_numbers)    // no smiles
  {
    IW_STL_Hash_Map_int::const_iterator f = id_to_ndx.find (id);
    if (f == id_to_ndx.end ())
    {
      cerr << "Writing neighbours as indices, but no data for '" << id << "'\n";
      return 0;
    }
    output << neighbour_tag << (*f).second << ">\n";

    return output.good ();
  }

  if (smiles_tag.length ())
  {
    output << smiles_tag;
    if (dummy_smiles.length ())
      output << dummy_smiles;
    else
    {
      IW_STL_Hash_Map_String::const_iterator f;

      if (strip_leading_zeros_from_identifiers && id.starts_with ('0'))
      {
        IWString tmp(id);
        tmp.remove_leading_chars ('0');
        f = smiles_for_id.find (tmp);
      }
      else
        f = smiles_for_id.find (id);

      if (f == smiles_for_id.end ())
      {
        if (smiles_for_id.size ())
          cerr << "Warning, no smiles for '" << id << "'\n";

        output << 'C';
      }
      else
        output << (*f).second;
    }

    output << ">\n";
  }

  output << identifier_tag << id;
  if (add_number_nbrs && nbrs >= 0)
    output << ' ' << nbrs;
  output << ">\n";

  return output.good ();
}

template <typename D>
int
Neighbour_List<D>::do_write (IWString_and_File_Descriptor & os) const
{
  for (int i = 0; i < _neighbours.number_elements (); i++)
  {
    const ID_Distance & sidi = *(_neighbours[i]);

    const IWString & id = sidi.id ();

    if (0 == id.length ())     // neighbour not initialised
      break;

    echo_smiles_and_id (id, smiles_for_haystack, write_neighbours_as_index_numbers, -1, os);
    os << distance_tag << sidi.distance () << ">\n";
    os.write_if_buffer_holds_more_than(16384);
  }

  return os.good ();
}

template <typename D>
int
Neighbour_List<D>::_binary_search (D distance) const
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

template <typename D>
void
Neighbour_List<D>::_check_insertion ()
{
  int failure = 0;

  int nn = _neighbours.number_elements ();
  for (int i = 1; i < nn; i++)
  {
    if (_neighbours[i - 1]->distance () > _neighbours[i]->distance ())
    {
      cerr << "Sort/insertion failed, out of order, i = " << i << endl;
      cerr << _neighbours[i - 1]->distance () << " vs " << _neighbours[i]->distance () << endl;
      if (keep_going_after_fatal_error)
      {
        _neighbours.remove_item (i);
        nn--;
      }
      failure = 1;
    }
  }

  if (failure && keep_going_after_fatal_error)
  {
    fatal_errors_encountered++;
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

template <typename D>
void
Neighbour_List<D>::extra (const IWString & rhs, D distance)
{
  if (distance > static_cast<D> (0.0))    // hopefully the most common case
    ;
  else if (static_cast<D> (0.0) == distance)
  {
    _number_zero_distance_neighbours++;

    if (_number_zero_distance_neighbours == stop_when_all_molecules_have_a_zero_distance_neighbour)
      molecules_with_zero_distance_neighbours++;
  }
  else if (negative_distances_allowed)
    ;
  else if (fabs(distance) < 1.0e-07)    // just roundoff
    distance = static_cast<D>(0.0);
  else
  {
    cerr << "Neighbour_List::extra: fatal error, distance " << distance << " to '" << rhs << "'\n";
    if (keep_going_after_fatal_error)
    {
      fatal_errors_encountered++;
      return;
    }

    abort ();
  }

  if (0 == _neighbours_to_find)
  {
    _extra_no_max_number_neighbours (rhs, distance);
    return;
  }

//cerr << "_neighbours array has " << _neighbours.number_elements() << " items, last " << _neighbours.last_item()->distance() << endl;

  if (distance >= _neighbours.last_item ()->distance ())
    return;

  ID_Distance * x = _neighbours.last_item ();

  int left = _binary_search (distance);

  // Shuffle everything one slot to the right
  for (int i = _neighbours.number_elements () - 1; i > left; i--)
  {
    _neighbours[i] = _neighbours[i - 1];
  }

  x->set_distance (distance);
  x->set_id (rhs);

  _neighbours[left] = x;

//#define CHECK_INSERTION
#ifdef CHECK_INSERTION
  _check_insertion ();
#endif

  return;
}

/*
  We are accepting all neighbours
*/

template <typename D>
void
Neighbour_List<D>::_extra_no_max_number_neighbours (const IWString & rhs, D distance)
{
  ID_Distance * sidi = new ID_Distance (rhs, distance);

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

template class Neighbour_List<distance_t>;

class Descriptors_and_Neighbours : public IWDescriptors<float, distance_t>,
                                   public Neighbour_List<distance_t>
{
  private:
  public:
    Descriptors_and_Neighbours() {
    }
    Descriptors_and_Neighbours(const Descriptors_and_Neighbours& rhs);
    Descriptors_and_Neighbours operator=(const Descriptors_and_Neighbours& rhs);
};

Descriptors_and_Neighbours::Descriptors_and_Neighbours(const Descriptors_and_Neighbours& rhs) {
  cerr << "Descriptors_and_Neighbours copy constructor not implemented\n";
  ::abort();
}

Descriptors_and_Neighbours
Descriptors_and_Neighbours::operator=(const Descriptors_and_Neighbours& rhs) {
  cerr << "Descriptors_and_Neighbours assignment operator not implemented\n";
  ::abort();
}

template class iwaray<Descriptors_and_Neighbours>;

template class Set_of_Descriptors<Descriptors_and_Neighbours>;

static distance_t zero_distance = static_cast<distance_t> (0.0);
/*old
static int
descriptor_similarity_erg (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       const IWDescriptors<float, distance_t> & d,
                       int number_descriptors,
                       const int * xref)
{
  int pool_size = pool.number_elements ();

  for (int i = 0; i < pool_size; i++)
  {
    Descriptors_and_Neighbours & pi = pool[i];

    distance_t dist = d.distance (pi, distance_computation, number_descriptors, xref);

    if (zero_distance != max_distance && dist > max_distance)
      continue;

    if (zero_distance == dist && do_not_compare_molecules_with_themselves && pi.id () == d.id ())
      continue;

    pi.extra (d.id (), dist);

    if (verbose)
      stats.extra (dist);
  }

  return 1;
}
endold*/

//nik
template <typename D>
int
Neighbour_List<D>::scan_for_closer (const IWString & rhs, D distance)
{
  int i,j;

  for (i = 0; i < _neighbours.number_elements(); i++)
  {
    if ((_neighbours[i]->id() == rhs) && (distance >= _neighbours[i]->distance()))
      return(0);
    else if ((_neighbours[i]->id() == rhs) && (distance < _neighbours[i]->distance()))
    {
      //since the address of _neighbours[i] will be needed for the last field of
      //_neighbours[j] save this in tmp
      ID_Distance * tmp = _neighbours[i];

      //delete _neighbours[i] entry so that the better one can be included
      //shuffle everything one left
      for (j = i; j < _neighbours.number_elements()-1; j++)
            _neighbours[j] = _neighbours[j+1];

      //currently the last field and last-1 field is pointing in the same
      //memory space! reset this for last into where i is pointing
      _neighbours[j] = tmp;
      //get the values of the last field correct
      (_neighbours[j])->set_distance(_neighbours[j-1]->distance());
      (_neighbours[j])->set_id(_neighbours[j-1]->id());

      return (1);
    }
  }

  //can reach here only if previousID was not stored in _neighbours

  return(1);
}

template <typename D>
D
Neighbour_List<D>::distance_of_closest_neighbour () const
{
  int n = _neighbours.number_elements();

  for (int i = 0; i < n; i++)
  {
    const ID_Distance * sidi = _neighbours[i];

    if (sidi->id().length() > 0)
      return sidi->distance();
  }

//cerr << "Neighbour_List::distance_of_closest_neighbour:no valid nbrs among " << n << " entries\n";

  return 1.0;
}

template <typename D>
D
Neighbour_List<D>::distance_of_furthest_neighbour () const
{
  for (int i = _neighbours.number_elements() - 1; i >= 0; i--)
  {
    const ID_Distance * sidi = _neighbours[i];

    if (sidi->id().length() > 0)
      return sidi->distance();
  }

  return 1.0;
}

static int
descriptor_similarity_erg (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                           resizable_array_p < IWDescriptors<float, distance_t> > & pd,
                          int number_descriptors,
                          const int * xref)
{
  int pool_size = pool.number_elements ();

  int n = pd.number_elements();

  if (DISTANCE_CONTINUOUS_TANIMOTO == distance_computation)
  {
    for (int i = 0; i < n; i++)
    {
      pd[i]->compute_product_of_non_zero_items(number_descriptors);
    }
  }


  const IWString & haystackID = pd[0]->id();

  for (int i = 0; i < pool_size; i++)
  {
    // since we can have multiple pool descriptors from FF and
    // multiple haystack descriptors, first loop over all haystacks
    // and then check on following pool id's and compute distance for
    // each of the pool IDs then use the best on first pool entry
    // (others will not be written!!)

    Descriptors_and_Neighbours & pi = pool[i];

    distance_t best_dist(2);

    for (int pd_i = 0; pd_i < n; pd_i++)
    {
      distance_t dist = (pd[pd_i])->distance (pi, distance_computation, number_descriptors, xref);
      if (dist < best_dist)
        best_dist = dist;
    }

    int j = i+1;
    while ((j < pool_size) && (pool[i].id() == pool[j].id()))
    {
      const Descriptors_and_Neighbours & next_pi = pool[j];
      for (int pd_i = 0; pd_i < n; pd_i++)
      {
        distance_t dist = (pd[pd_i])->distance (next_pi, distance_computation, number_descriptors, xref);
        if (dist < best_dist)
          best_dist = dist;
      }
      j++;
    }
    //caution!! we are changing i within the i for loop!!
    //necessary to jump over the duplicates
    i = j-1;

    if (min_neighbours_to_find > 0 && pi.number_neighbours() < min_neighbours_to_find)
      ;
//  else if (zero_distance != max_distance && best_dist > max_distance)
    else if (best_dist > max_distance)
      continue;

    if (zero_distance == best_dist && do_not_compare_molecules_with_themselves && pi.id () == haystackID)
      continue;

//  cerr << " i = " << i << " id " << pi.id() << " best_dist " << best_dist <<  " '" << haystackID << "'\n";

    pi.extra (haystackID, best_dist);

    if (verbose)
      stats.extra (best_dist);
  }

  return 1;
}

static int
compute_number_neighbours (const Set_of_Descriptors<Descriptors_and_Neighbours> & pool)
{
  int rc = 0;

  int n = pool.number_elements();

  for (int i = 0; i < n; i++)
  {
    rc += pool[i].number_neighbours();
  }

  return rc;
}

static int
descriptor_similarity_erg (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       iwstring_data_source & input,
                       int number_descriptors,
                       const int * xref)
{
  IWString previousID = "";

  int write_error = 1;

  resizable_array_p < IWDescriptors<float, distance_t> > pd;
  pd.resize(MAX_IDs);

  const_IWSubstring buffer;

  while (input.next_record (buffer))   // get first entry
  {
    IWDescriptors<float, distance_t> * d = new IWDescriptors<float, distance_t>;
    int fatal;
    if (! d->build (buffer, number_descriptors, fatal))
    {
      delete d;
      cerr << "Invalid descriptor record at line " << input.lines_read () << endl;
      if (fatal)
      {
        cerr << buffer << endl;
        return 0;
      }
      continue;
    }
    pd.add(d);
    previousID = d->id();
    break;
  }

  while (input.next_record (buffer))
  {
    haystack_records_read++;

    if (report_progress > 0 && 0 == haystack_records_read % report_progress)
    {
      int nbrs = compute_number_neighbours(pool);
      cerr << "Processed " << haystack_records_read << " haystack records, storing " << nbrs << " neighbours\n";
    }

    if (buffer.word(0) == previousID)
    {
     if (pd.number_elements() < MAX_IDs)
     {
        IWDescriptors<float, distance_t> * d = new IWDescriptors<float, distance_t>;
        int fatal;
        if (! d->build (buffer, number_descriptors, fatal))
        {
          delete d;
          cerr << "Invalid descriptor record at line " << input.lines_read () << endl;
          if (fatal)
          {
            cerr << buffer << endl;
            return 0;
          }
          continue;
        }

        pd.add(d);
      }
      else if (write_error)
      {
        cerr << "Compound " << previousID << " maximum number of entries reached, skipping other enumerated forms\n";
        write_error = 0;
      }
    }
    else
    {
      if (0 == pd.number_elements())
        cerr << "Ignoring identifier with no data '" << previousID << "'\n";
      else if (! descriptor_similarity_erg (pool, pd, number_descriptors, xref))
      {
         cerr << "Fatal error processing '" << previousID << "' at line " << input.lines_read () << endl;
         return 0;
      }

      //reset the stuff
      pd.resize_keep_storage(0);

      IWDescriptors<float, distance_t> * d = new IWDescriptors<float, distance_t>;
      int fatal;
      if (! d->build (buffer, number_descriptors, fatal))
      {
        delete d;

        cerr << "Invalid descriptor record at line " << input.lines_read () << endl;
        if (fatal)
        {
          cerr << buffer << endl;
          return 0;
        }
        continue;
      }

      pd.add(d);
      previousID = d->id();
      write_error = 1;
    }
  }

  if (stop_when_all_molecules_have_a_zero_distance_neighbour && pool.number_elements () == molecules_with_zero_distance_neighbours)
  {
    //have to do the last writing!
    if (! descriptor_similarity_erg (pool, pd, number_descriptors, xref))
    {
      cerr << "Fatal error processing '" << previousID << "' at line " << input.lines_read () << endl;
      return 0;
    }

    if (verbose)
      cerr << "All " << pool.number_elements () << " molecules have a zero distance neighbour\n";

    return 1;
  }

  //get the last entry
  if (! descriptor_similarity_erg (pool, pd, number_descriptors, xref))
  {
    cerr << "Fatal error processing '" << previousID << "' at line " << input.lines_read () << endl;
    return 0;
  }

  return 1;
}

/*
  We need to make sure that every descriptor in HEADER is also in MYHEADER.
  We also need to establish the cross reference array which will map descriptors
  in MYHEADER to the corresponding column in HEADER
  pool to reflect the same ordering as the input file
*/

static int
determine_cross_reference (const IWString & header,
                           const IWString & myheader,
                           int * xref)
{
  if (! haystack_file_contains_header_record)
  {
    int nw = header.nwords();

    if (nw != myheader.nwords())
    {
      cerr << "Token count mismatch for no header files, " << nw << " vs " << myheader.nwords() << endl;
      return 0;
    }

    for (int i = 0; i < nw; i++)
    {
      xref[i] = i;
    }

    return 1;
  }

  int i = 0;
  IWString token;
  int col = 0;

  IW_STL_Hash_Map_int descriptor_to_column;

  while (header.nextword (token, i))
  {
    descriptor_to_column[token] = col;
    col++;
  }

  int columns_assigned = 0;

  i = 0;
  col = 0;
  while (myheader.nextword (token, i))
  {
    IW_STL_Hash_Map_int::const_iterator f = descriptor_to_column.find (token);
    if (f != descriptor_to_column.end ())     // descriptor in HEADER
    {
      xref[col] = (*f).second;
      columns_assigned++;
    }

    if (verbose > 1)
      cerr << "Descriptor '" << token << "' in column " << xref[col] << " found in column " << col << endl;

    col++;
  }

  if (columns_assigned == col)
    ;
  else if (0 == columns_assigned)
  {
    cerr << "No column cross reference, cannot continue\n";
    cerr << header << endl;
    cerr << myheader << endl;
    return 0;
  }
  else
  {
    cerr << "Warning, only " << columns_assigned << " of " << col << " columns assigned cross reference values\n";
  }

  return 1;
}

static int
descriptor_similarity_erg (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       iwstring_data_source & input,
                       const IWString & header,
                       const IWString & myheader,
                       int number_descriptors,
                       int * xref)
{
  if (verbose)
    cerr << number_descriptors << " descriptors\n";

  if (! determine_cross_reference (header, myheader, xref))
  {
    cerr << "Cannot reorder the input descriptors to be the same as the pool\n";
    return 0;
  }

  return descriptor_similarity_erg (pool, input, number_descriptors, xref);
}

static int
descriptor_similarity_erg (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       iwstring_data_source & input,
                       const IWString & header)
{
  IWString myheader;
  if (! input.next_record (myheader))
  {
    cerr << "Cannot read header from input\n";
    return 0;
  }

  myheader.remove_leading_words (1);

  int myhn = myheader.nwords ();

  if (myhn < header.nwords ())
  {
    cerr << "Pool file contains " << header.nwords () << " but my file contains only " << myhn << ", impossible\n";
    return 0;
  }

  if (myhn <= 0)
  {
    cerr << "Gack, header is empty!!\n";
    return 0;
  }

  int * xref = new_int (myhn, -1); std::unique_ptr<int[]> free_xref (xref);

  return descriptor_similarity_erg (pool, input, header, myheader, myhn, xref);
}

static int
descriptor_similarity_erg (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       const char * fname,
                       const IWString & header)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return descriptor_similarity_erg (pool, input, header);
}

static int
read_smiles (const const_IWSubstring & buffer,
             IW_STL_Hash_Map_String & id_to_smiles)
{
  const_IWSubstring smiles, id;
  if (! buffer.split (smiles, ' ', id))
  {
    cerr << "Cannot tokenise as smiles and id\n";
    return 0;
  }

  id.truncate_at_first (' ');

  if (strip_leading_zeros_from_identifiers)
    id.remove_leading_chars ('0');

  id_to_smiles[id] = smiles;

  return 1;
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! read_smiles (buffer, smiles))
    {
      cerr << "Fatal error processing smiles file on line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return smiles.size ();
}

static int
read_smiles_from_file (const const_IWSubstring & fname,
                       IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles (input, smiles);
}

static int
fill_id_to_ndx (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                                   IW_STL_Hash_Map_int & id_to_ndx)
{
  int n = pool.number_elements ();

  for (int i = 0; i < n; i++)
  {
    const Descriptors_and_Neighbours & dni = pool[i];

    const IWString & id = dni.id ();

    int j = id_to_ndx.size ();

    id_to_ndx[id] = j;

//  cerr << "Index for " << id << "' at " << id_to_ndx[id] << endl;
  }

  if (0 == id_to_ndx.size ())
  {
    cerr << "Yipes, no index data\n";
    return 0;
  }

  return id_to_ndx.size ();
}

static int
descriptor_similarity_erg (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "p:n:s:vT:m:X:hoEP:D:Z:zr:ya");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', max_distance) || max_distance <= 0.0)
    {
      cerr << "Invalid maximum distance value (-T option)\n";
      usage (3);
    }

    if (verbose)
      cerr << "Will ignore distances > " << max_distance << endl;
  }

  if (cl.option_present ('E'))
  {
    set_descriptors_may_contain_big_E (1);

    if (verbose)
      cerr << "Will translate E to e for numeric input\n";
  }

  if (cl.option_present('y'))
  {
    haystack_file_contains_header_record = 0;
    if (verbose)
      cerr << "Will assume the haystack file does not contain a header record\n";
  }

  if (cl.option_present('a'))
  {
    add_number_nbrs = 1;

    if (verbose)
      cerr << "Will add the number of neighbours to each needle\n";
  }

  Set_of_Descriptors<Descriptors_and_Neighbours> pool;

  if (cl.option_present ('s'))
  {
    int pool_size;
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "Invalid -s option\n";
      usage (4);
    }

    if (! pool.resize (pool_size))
    {
      cerr << "Bad news, cannot size descriptors for " << pool_size << " descriptors\n";
      return 3;
    }

    if (verbose)
      cerr << "Pool sized to " << pool_size << endl;
  }

  if (cl.option_present ('D'))
  {
    if (cl.option_present ('P'))
    {
      cerr << "Cannot specify a dummy smiles (-D) with either -p or -P\n";
      usage (5);
    }

    dummy_smiles = cl.string_value ('D');

    if (verbose)
      cerr << "Will use '" << dummy_smiles << "' as the smiles\n";

    smiles_tag = "$SMI<";
  }

  IWString header;     // the header in the -p file

  if (! cl.option_present ('p'))
  {
    cerr << "Must specify the target descriptor file via the -p option\n";
    usage (3);
  }
  else
  {
    IWString p = cl.string_value ('p');
    if (! pool.build (p, header))
    {
      cerr << "Cannot initialise descriptor pool file '" << p << "'\n";
      return 12;
    }
  }

  int neighbours_to_find = -1;
  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', neighbours_to_find) || neighbours_to_find < 1)
    {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage (13);
    }
  }

  if (neighbours_to_find > 0)
  {
    for (int i = 0; i < pool.number_elements (); i++)
    {
      pool[i].set_neighbours_to_find (neighbours_to_find);
    }

    if (verbose)
      cerr << "A maximum of " << neighbours_to_find << " neighbours of each molecule will be found\n";
  }

  if (cl.option_present ('m'))
  {
    if (! cl.value('m', min_neighbours_to_find) || min_neighbours_to_find < 1)
    {
      cerr << "the -m option must be followed by a whole positive number\n";
      usage (13);
    }

    if (verbose)
      cerr << "Will find at least " << min_neighbours_to_find << " neighbours\n";
  }

  if (cl.option_present ('Z'))
  {
    if (! cl.value ('Z', stop_when_all_molecules_have_a_zero_distance_neighbour) || stop_when_all_molecules_have_a_zero_distance_neighbour < 1)
    {
      cerr << "The stop when zero distance neighbours count option (-Z) must be a whole positive number\n";
      usage (5);
    }

    if (neighbours_to_find > 0 && stop_when_all_molecules_have_a_zero_distance_neighbour > neighbours_to_find)
    {
      cerr << "Only finding " << neighbours_to_find << " neighbours, will never have " << stop_when_all_molecules_have_a_zero_distance_neighbour << " neighbours\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will stop when every 'needle' molecule has " << stop_when_all_molecules_have_a_zero_distance_neighbour << " zero distance neighbours\n";
  }

  if (cl.option_present ('z'))
  {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Leading 0 characters stripped from identifiers\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', report_progress) || report_progress < 1)
    {
      cerr << "The report progress option (-r) must be a whole +ve number\n";
      usage(8);
    }

    if (verbose)
      cerr << "Will report progress every " << report_progress << " haystack items processed\n";
  }

  if (cl.option_present ('X'))
  {
    const_IWSubstring x = cl.string_value ('X');

    distance_computation = determine_distance_type (x, verbose);
    if (0 == distance_computation)
    {
      cerr << "Cannot determine distance type\n";
      usage (4);
    }
  }

  if (cl.option_present ('h'))
  {
    do_not_compare_molecules_with_themselves = 1;

    if (verbose)
      cerr << "Will discard neighbours with zero distance and the same id as the target\n";
  }

  if (cl.option_present ('o'))
  {
    write_neighbours_as_index_numbers = 1;

    if (verbose)
      cerr << "Neighbours written as index numbers\n";

    if (! fill_id_to_ndx (pool, id_to_ndx))
      return 5;
  }

  if (DISTANCE_CONTINUOUS_TANIMOTO == distance_computation)
  {
    int pool_size = pool.number_elements();
    int number_descriptors = header.nwords();

    for (int i = 0; i < pool_size; i++)
    {
      pool[i].compute_product_of_non_zero_items(number_descriptors);
    }
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (0 == cl.option_count ('P'))
    ;
  else if (cl.option_count ('P') > 2)
  {
    cerr << "A maximum of two -P options are possible\n";
    usage (2);
  }
  else
  {
    int gotneedles = 0;
    int gothaystack = 0;
    int i = 0;
    const_IWSubstring p;
    while (cl.value ('P', p, i++))
    {
      if (p.starts_with ("p:"))
      {
        gotneedles++;
        p.remove_leading_chars (2);
        if (! read_smiles_from_file (p, smiles_for_needles))
        {
          cerr << "Cannot read -p smiles from '" << p << "'\n";
          return 3;
        }
      }
      else
      {
        gothaystack++;
        if (! read_smiles_from_file (p, smiles_for_haystack))
        {
          cerr << "Cannot read haystack smiles\n";
          return 4;
        }
      }
    }

    smiles_tag = "$SMI<";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! descriptor_similarity_erg (pool, cl[i], header))
    {
      rc = i + 1;
      break;
    }
  }

  set_default_iwstring_float_concatenation_precision(4);

  IWString_and_File_Descriptor output(1);

//nik
  for (int i = 0; i < pool.number_elements (); i++)
  {
    int j;
    const IWString & id = pool[i].id ();

    echo_smiles_and_id (id, smiles_for_needles, 0, pool[i].number_neighbours(), output);   // write PCN for the target ID

    pool[i].do_write (output);

    for (j = i+1; j < pool.number_elements (); j++)
    {
      if (pool[i].id() != pool[j].id())
        break;
    }

    output << "|\n";
    output.write_if_buffer_holds_more_than(16384);
    // get the i correct (jump over duplicates)
    i = j-1;
  }

  if (verbose)
  {
    cerr << "Distances encountered between " << stats.minval () << " and " << stats.maxval ();
    if (stats.n () > 1)
      cerr << " average " << stats.average ();
    cerr << endl;

    Accumulator<distance_t> mind, maxd;
    for (int i = 0; i < pool.number_elements (); i++)
    {
      const Descriptors_and_Neighbours & dni = pool[i];

      if (i > 0 && dni.id() == pool[i - 1].id())
        continue;

      mind.extra (dni.distance_of_closest_neighbour ());
      maxd.extra (dni.distance_of_furthest_neighbour ());
    }

    cerr << "Of neighbours stored, closest distances between " << mind.minval () << " and " << mind.maxval ();
    if (mind.n () > 1)
      cerr << " ave " << mind.average ();
    cerr << endl;

    cerr << maxd.n () << " furthest neighbours between " << maxd.minval () << " and " << maxd.maxval ();
    if (maxd.n () > 1)
      cerr << " ave " << maxd.average ();
    cerr << endl;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_similarity_erg (argc, argv);

  return rc;
}
