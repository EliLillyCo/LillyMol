/*
  Scans a descriptor file for similarity to a given vector
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWARAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwqsort.h"
#include "iw_stl_hash_map.h"
#include "cmdline.h"
#include "misc.h"
#include "iwstring_data_source.h"
#include "accumulator.h"

#define SET_OF_DESCRIPTORS_IMPLEMENTATION
#include "set_of_descriptors.h"

#include "iwdescriptor.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

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

/*
  Some distance metrics may generate negative distances
*/

static int negative_distances_allowed = 0;

static int distance_computation = DISTANCE_CARTESIAN;

static int do_not_compare_molecules_with_themselves = 0;

static int stop_when_all_molecules_have_a_zero_distance_neighbour = 0;

static int molecules_with_zero_distance_neighbours = 0;

static int strip_leading_zeros_from_identifiers = 0;

static int write_neighbours_as_index_numbers = 0;

static int haystack_records_read = 0;

static int report_progress = 0;

static float subtract_from = static_cast<float>(0.0);

static int min_neighbours_to_find = -1;

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
  cerr << " -p <file>        file of molecules for which neighbours are to be found \"needles\"\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      number of neighbours to keep\n";
  cerr << " -m <number>      minimum number of neighbours to keep\n";
  cerr << " -T <dis>         discard distances longer than <dis>\n";
  display_distance_metrics (cerr, 'X');
  cerr << " -h               discard neighbours with zero distance and the same ID as the target\n";
  cerr << " -o               write neighbours as indes numbers rather than id's\n";
  cerr << " -P <file>        smiles file for the \"haystack\" molecules\n";
  cerr << " -P p:<file>      smiles file for -p (\"needle\")descriptors\n";
  cerr << " -D <smiles>      use a dummy smiles wherever needed\n";
  cerr << " -O <dist>        distance offset, computed distances subtracted from <dist>\n";
  cerr << " -z               strip leading 0's from identifiers in matching up smiles\n";
  cerr << " -W ...           comparison window specification\n";
  cerr << " -Z <number>      stop once every 'needle' has <number> zero distance neighbour(s)\n";
  cerr << " -r <number>      report progress every <number> haystack items processed\n";
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
    ID_Dist_Template();
    ID_Dist_Template (const IWString & s, D d);
    ~ID_Dist_Template();

    D distance() const { return _distance;}
    const IWString & id() const { return _id;}

    void set_distance (D d) { _distance = d;}
    void set_id (const IWString & s) { _id = s;}
};

template <typename D>
ID_Dist_Template<D>::ID_Dist_Template()
{
  _distance = numeric_limits<D>::max();

  return;
}

template <typename D>
ID_Dist_Template<D>::ID_Dist_Template (const IWString & s, D d) : _id (s), _distance (d)
{
}

template <typename D>
ID_Dist_Template<D>::~ID_Dist_Template()
{
}

typedef float distance_t;

static Accumulator<distance_t> stats;

#ifdef __GNUG__
template class ID_Dist_Template<distance_t>;
#endif

static distance_t max_distance = numeric_limits<distance_t>::max();

class ID_Distance : public ID_Dist_Template<distance_t>
{
  private:
  public:
    ID_Distance();
    ID_Distance (const IWString &, distance_t);
};

ID_Distance::ID_Distance()
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
    void _check_insertion();

  public:
    Neighbour_List();
    ~Neighbour_List();

    int set_neighbours_to_find (int);

    int number_neighbours() const { return _neighbours.number_elements();}

    void sort_neighbours ();

    void extra (const IWString &, D);
    int  do_write (IWString_and_File_Descriptor &) const;

    D distance_of_closest_neighbour() const { return _neighbours[0]->distance();}
    D distance_of_furthest_neighbour() const { return _neighbours.last_item()->distance();}
};

template class resizable_array_p<ID_Distance>;
template class resizable_array_base<ID_Distance *>;

template <typename D>
Neighbour_List<D>::Neighbour_List()
{
  _neighbours_to_find = 0;

  _number_zero_distance_neighbours = 0;

  return;
}

template <typename D>
Neighbour_List<D>::~Neighbour_List()
{
  _neighbours_to_find = -2;
}

template <typename D>
int
Neighbour_List<D>::set_neighbours_to_find (int s)
{
  assert (s > 0);

  _neighbours_to_find = s;

  if (_neighbours.number_elements())    // clean out any existing neighbour data
    _neighbours.resize (0);

  _neighbours.resize (s);

  for (int i = 0; i < _neighbours_to_find; i++)
  {
    ID_Distance * s = new ID_Distance();
    _neighbours.add (s);
  }

  return 1;
}

int
echo_smiles_and_id (const IWString & id,
                    const IW_STL_Hash_Map_String & smiles_for_id,
                    int write_neighbours_as_index_numbers,
                    IWString_and_File_Descriptor & output)
{
  if (write_neighbours_as_index_numbers)    // no smiles
  {
    IW_STL_Hash_Map_int::const_iterator f = id_to_ndx.find (id);
    if (f == id_to_ndx.end())
    {
      cerr << "Writing neighbours as indices, but no data for '" << id << "'\n";
      return 0;
    }
    output << neighbour_tag << (*f).second << ">\n";

    return output.good();
  }

  if (smiles_tag.length())
  {
    output << smiles_tag;
    if (dummy_smiles.length())
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

      if (f == smiles_for_id.end())
      {
        if (smiles_for_id.size())
          cerr << "Warning, no smiles for '" << id << "'\n";
    
        output << 'C';
      }
      else
        output << (*f).second;
    }

    output << ">\n";
  }

  output << identifier_tag << id << ">\n";

  output.write_if_buffer_holds_more_than(8192);

  return output.good();
}

template <typename D>
int
Neighbour_List<D>::do_write (IWString_and_File_Descriptor & os) const
{
  for (int i = 0; i < _neighbours.number_elements(); i++)
  {
    const ID_Distance & sidi = *(_neighbours[i]);

    const IWString & id = sidi.id();

    if (0 == id.length())     // neighbour not initialised
      break;

    echo_smiles_and_id (id, smiles_for_haystack, write_neighbours_as_index_numbers, os);
    os << distance_tag << sidi.distance() << ">\n";
  }

  return os.good();
}

template <typename D>
int
Neighbour_List<D>::_binary_search (D distance) const
{
  int left = 0;      // _neighbours[left] will hold the new value

  if (distance > _neighbours[0]->distance())     // goes after neighbour 0
  {
    int right = _neighbours.number_elements() - 1;

    while (right > left)
    {
      int middle = (left + right) / 2;
      if (middle == left)
      {
        left = right;
        break;
      }

      D dmid = _neighbours[middle]->distance();
//    cerr << "Left = " << left << " d = " << _neighbours[left]->distance() << " middle " << middle << " d = " << _neighbours[middle].distance() << " right " << right << " d = " << _neighbours[right].distance() << endl;
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
Neighbour_List<D>::_check_insertion() 
{
  int failure = 0;

  int nn = _neighbours.number_elements();
  for (int i = 1; i < nn; i++)
  {
    if (_neighbours[i - 1]->distance() > _neighbours[i]->distance())
    {
      cerr << "Sort/insertion failed, out of order, i = " << i << endl;
      cerr << _neighbours[i - 1]->distance() << " vs " << _neighbours[i]->distance() << endl;
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
    for (int i = 0; i < _neighbours.number_elements(); i++)
    {
      cerr << "i = " << i << " distance " << _neighbours[i]->distance() << endl;
    }

    exit (87);
  }

  return;
}

template <typename D>
void
Neighbour_List<D>::extra (const IWString & rhs, D distance)
{
  if (distance > static_cast<D>(0.0))    // hopefully the most common case
    ;
  else if (static_cast<D>(0.0) == distance)
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

    abort();
  }

  if (0 == _neighbours_to_find)
  {
    _extra_no_max_number_neighbours (rhs, distance);
    return;
  }

  if (distance >= _neighbours.last_item()->distance())
    return;

  ID_Distance * x = _neighbours.last_item();

  int left = _binary_search (distance);

// Shuffle everything one slot to the right

  for (int i = _neighbours.number_elements() - 1; i > left; i--)
  {
    _neighbours[i] = _neighbours[i - 1];
  }

  x->set_distance (distance);
  x->set_id (rhs);

  _neighbours[left] = x;

//#define CHECK_INSERTION
#ifdef CHECK_INSERTION
  _check_insertion();
#endif

  return;
}

/*
  We are accepting all neighbours. Just add them and sort at the end
*/

template <typename D>
void
Neighbour_List<D>::_extra_no_max_number_neighbours (const IWString & rhs, D distance)
{
  ID_Distance * sidi = new ID_Distance (rhs, distance);

  _neighbours.add(sidi);

  return;
}

/*template <typename D>
void
Neighbour_List<D>::_extra_no_max_number_neighbours (const IWString & rhs, D distance)
{
  ID_Distance * sidi = new ID_Distance (rhs, distance);

  if (0 == _neighbours.number_elements())
  {
    _neighbours.add (sidi);

    return;
  }

  if (distance >= _neighbours.last_item()->distance())
  {
    _neighbours.add (sidi);

    return;
  }

  int ins = _binary_search (distance);

  _neighbours.insert_before (ins, sidi);

#ifdef CHECK_INSERTION
  _check_insertion();
#endif

  return;
}*/

class Neighbour_distance_Comparator
{
  private:
  public:
    int operator() (const ID_Distance *, const ID_Distance *) const;
};

int
Neighbour_distance_Comparator::operator () (const ID_Distance * idd1, 
                                            const ID_Distance * idd2) const
{
  distance_t d1 = idd1->distance();
  distance_t d2 = idd2->distance();

  if (d1 < d2)
    return -1;

  if (d1 > d2)
    return 1;

  return 0;
}

static Neighbour_distance_Comparator ndc;

template <typename D>
void
Neighbour_List<D>::sort_neighbours()
{
  _neighbours.iwqsort(ndc);
}

template class Neighbour_List<distance_t>;

class Descriptors_and_Neighbours : public IWDescriptors<float, distance_t>,
                                   public Neighbour_List<distance_t>
{
  private:
  public:
};

template class iwaray<Descriptors_and_Neighbours>;

template class Set_of_Descriptors<Descriptors_and_Neighbours>;

/*
  We can speed things up a lot by skipping computations where some
  property might differ by a large amount between two items.
  The window consists of a column number and a delta

  Syntax looks like

  col:delta

  or if no col:, then it assumed to be column 1
*/

class Comparison_Window
{
  private:
    int _ndx;
    distance_t _delta;

  public:
    Comparison_Window();

    int build (const const_IWSubstring &);

    int column() const { return _ndx;}

    int too_large_a_difference(const distance_t s) const { return s > _delta;}
};

Comparison_Window::Comparison_Window()
{
  _ndx = -1;
  _delta = static_cast<distance_t>(0.0);

  return;
}

int
Comparison_Window::build(const const_IWSubstring & s)
{
  if (! s.contains(':'))
  {
    if (! s.numeric_value(_delta) || _delta <= 0.0)
    {
      cerr << "Comparison_Window::build:invalid delta '" << s << "'\n";
      return 0;
    }
  
    _ndx = 0;

    return 1;
  }

  const_IWSubstring c, d;

  s.split(c, ':', d);

  if (0 == c.length() || 0 == d.length())
  {
    cerr << "Comparison_Window::build:cannot split directive '" << s << "'\n";
    return 0;
  }

  if (! c.numeric_value(_ndx) || _ndx < 1)
  {
    cerr << "Comparison_Window::build:invalid column '" << s << "'\n";
    return 0;
  }

  if (! d.numeric_value(_delta) || _delta <= 0.0)
  {
      cerr << "Comparison_Window::build:invalid delta '" << s << "'\n";
      return 0;
  }

  _ndx--;

  return 1;
}

static Comparison_Window * window = NULL;
static int nwindow = 0;

static int computations_avoided_because_of_window = 0;

static distance_t zero_distance = static_cast<distance_t> (0.0);

static int
can_be_compared(const IWDescriptors<float, distance_t> & d1,
                const IWDescriptors<float, distance_t> & d2,
                const int * xref)
{
  const float * r1 = d1.rawdata();
  const float * r2 = d2.rawdata();

  for (int i = 0; i < nwindow; i++)
  {
    const Comparison_Window & ci = window[i];

    int col = ci.column();

    if (col[xref] < 0)
      continue;

    float v1 = r1[col];
    float v2 = r2[xref[col]];

    if (ci.too_large_a_difference(fabs(v1 - v2)))
    {
      computations_avoided_because_of_window++;
      return 0;
    }
  }

  return 1;
}

class Request_Average_STD_Computation
{
  private:
    int _number_descriptors;

  public:
    Request_Average_STD_Computation(int n) : _number_descriptors(n) {};

    int operator() (IWDescriptors<float, distance_t> &) const;
};

int
Request_Average_STD_Computation::operator() (IWDescriptors<float, distance_t> & d) const
{
  d.compute_pearson_data(_number_descriptors);

  return 1;
}

static int
descriptor_similarity (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       IWDescriptors<float, distance_t> & d,
                       int number_descriptors,
                       const int * xref)
{
  if (needs_ave_and_variance(distance_computation))
    d.compute_pearson_data(number_descriptors);
  else if (DISTANCE_TANIMOTO == distance_computation)
    d.compute_nset(number_descriptors);
  else if (DISTANCE_CONTINUOUS_TANIMOTO == distance_computation)
    d.compute_product_of_non_zero_items(number_descriptors);

  int pool_size = pool.number_elements();

  for (int i = 0; i < pool_size; i++)
  {
    Descriptors_and_Neighbours & pi = pool[i];

    if (nwindow && ! can_be_compared(d, pi, xref))
      continue;

    distance_t dist = d.distance(pi, distance_computation, number_descriptors, xref);

    if (min_neighbours_to_find > 0 && pi.number_neighbours() < min_neighbours_to_find)
      ;
    else if (dist > max_distance)
      continue;

    if (zero_distance == dist && do_not_compare_molecules_with_themselves && pi.id() == d.id())
      continue;

    if (subtract_from > zero_distance)
    {
      dist = subtract_from - dist;
      if (dist < zero_distance)
        dist = zero_distance;
    }

    pi.extra(d.id(), dist);

    if (verbose)
      stats.extra (dist);
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
descriptor_similarity (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       iwstring_data_source & input,
                       int number_descriptors,
                       const int * xref)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    haystack_records_read++;

    if (report_progress > 0 && 0 == haystack_records_read % report_progress)
    {
      int nbrs = compute_number_neighbours (pool);
      cerr << "Processed " << haystack_records_read << " haystack records, storing " << nbrs << " neighbours\n";
    }

    IWDescriptors<float, distance_t> d;
    int fatal;
    if (! d.build (buffer, number_descriptors, fatal))
    {
      cerr << "Invalid descriptor record at line " << input.lines_read() << endl;
      cerr << buffer << endl;
      if (fatal)
        return 0;
      continue;
    }

    if (! descriptor_similarity(pool, d, number_descriptors, xref))
    {
      cerr << "Fatal error processing '" << d.id() << "' at line " << input.lines_read() << endl;
      return 0;
    }

    if (stop_when_all_molecules_have_a_zero_distance_neighbour && pool.number_elements() == molecules_with_zero_distance_neighbours)
    {
      if (verbose)
        cerr << "All " << pool.number_elements() << " molecules have a zero distance neighbour\n";
      return 1;
    }
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
    if (f != descriptor_to_column.end())     // descriptor in HEADER
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
descriptor_similarity (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
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

  return descriptor_similarity (pool, input, number_descriptors, xref);
}

static int
descriptor_similarity (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
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

  int myhn = myheader.nwords();

  if (myhn < header.nwords())
  {
    cerr << "Pool file contains " << header.nwords() << " but my file contains only " << myhn << ", impossible\n";
    return 0;
  }

  if (myhn <= 0)
  {
    cerr << "Gack, header is empty!!\n";
    return 0;
  }

  int * xref = new_int (myhn, -1); std::unique_ptr<int[]> free_xref (xref);

  return descriptor_similarity (pool, input, header, myheader, myhn, xref);
}

static int
descriptor_similarity (Set_of_Descriptors<Descriptors_and_Neighbours> & pool,
                       const char * fname,
                       const IWString & header)
{
  iwstring_data_source input (fname);
  if (! input.ok())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return descriptor_similarity (pool, input, header);
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
      cerr << "Fatal error processing smiles file on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return smiles.size();
}

static int
read_smiles_from_file (const const_IWSubstring & fname,
                       IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input (fname);
  if (! input.ok())
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
  int n = pool.number_elements();

  for (int i = 0; i < n; i++)
  {
    const Descriptors_and_Neighbours & dni = pool[i];

    const IWString & id = dni.id();

    int j = id_to_ndx.size();

    id_to_ndx[id] = j;

//  cerr << "Index for " << id << "' at " << id_to_ndx[id] << endl;
  }

  if (0 == id_to_ndx.size())
  {
    cerr << "Yipes, no index data\n";
    return 0;
  }

  return id_to_ndx.size();
}

static int
descriptor_similarity (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "p:n:s:vT:n:m:X:hoEP:D:Z:zr:O:W:eq");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  set_default_iwstring_float_concatenation_precision(4);

  if (cl.option_present ('T'))
  {
    if (! cl.value('T', max_distance) || max_distance <= 0.0)
    {
      cerr << "Invalid maximum distance value(-T option)\n";
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

  if (cl.option_present('O'))
  {
    if (! cl.value('O', subtract_from) || subtract_from <= 0.0)
    {
      cerr << "The subtract from option (-O) must be a valid +ve distance\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will subtract all distances from " << subtract_from << '\n';
  }

  if (cl.option_present('W'))
  {
    nwindow = cl.option_count('W');

    window = new Comparison_Window[nwindow];

    for (int i = 0; i < nwindow; i++)
    {
      const_IWSubstring w = cl.string_value('W', i);

      if (! window[i].build(w))
      {
        cerr << "Cannot parse comparison window specification '" << w << "'\n";
        return 5;
      }
    }

    if (verbose)
      cerr << "Defined " << nwindow << " comparison windows\n";
  }

  Set_of_Descriptors<Descriptors_and_Neighbours> pool;

  if (cl.option_present ('s'))
  {
    int pool_size;
    if (! cl.value('s', pool_size) || pool_size < 1)
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

    dummy_smiles = cl.string_value('D');

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
    IWString p = cl.string_value('p');
    if (! pool.build (p, header))
    {
      cerr << "Cannot initialise descriptor pool file '" << p << "'\n";
      return 12;
    }
  }

  int neighbours_to_find = -1;
  if (cl.option_present ('n'))
  {
    if (! cl.value('n', neighbours_to_find) || neighbours_to_find < 1)
    {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage (13);
    }
  }

  if (neighbours_to_find > 0)
  {
    for (int i = 0; i < pool.number_elements(); i++)
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

  if (cl.option_present('e'))
  {
    negative_distances_allowed = 1;
    if (verbose)
      cerr << "Negative distances allowed\n";
  }

  if (cl.option_present ('Z'))
  {
    if (! cl.value('Z', stop_when_all_molecules_have_a_zero_distance_neighbour) || stop_when_all_molecules_have_a_zero_distance_neighbour < 1)
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
    const_IWSubstring x = cl.string_value('X');

    distance_computation = determine_distance_type (x, verbose);
    if (0 == distance_computation)
    {
      cerr << "Cannot determine distance type\n";
      usage (4);
    }
  }

  if (needs_ave_and_variance(distance_computation))
  {
    Request_Average_STD_Computation rastd(header.nwords());
    pool.each(rastd);
  }

  if (DISTANCE_TANIMOTO == distance_computation)
  {
    int pool_size = pool.number_elements();
    int number_descriptors = header.nwords();
    for (int i = 0; i < pool_size; i++)
    {
      pool[i].compute_nset(number_descriptors);
//    cerr << " i = " << i <<  " nset " << pool[i].nset() << endl;
    }
  }
  else if (DISTANCE_CONTINUOUS_TANIMOTO == distance_computation)
  {
    int pool_size = pool.number_elements();
    int number_descriptors = header.nwords();
    for (int i = 0; i < pool_size; i++)
    {
      pool[i].compute_product_of_non_zero_items(number_descriptors);
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

  if (0 == cl.number_elements())
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
    while (cl.value('P', p, i++))
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
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! descriptor_similarity (pool, cl[i], header))
    {
      rc = i + 1;
      break;
    }
  }

  if (neighbours_to_find < 0)
  {
    int pool_size = pool.number_elements();
    for (int i = 0; i < pool_size; i++)
    {
      pool[i].sort_neighbours();
    }
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < pool.number_elements(); i++)
  {
    const IWString & id = pool[i].id();

    echo_smiles_and_id (id, smiles_for_needles, 0, output);   // write PCN for the target ID

    pool[i].do_write (output);
    output << "|\n";

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  if (verbose)
  {
    cerr << "Distances encountered between " << stats.minval() << " and " << stats.maxval();
    if (stats.n() > 1)
      cerr << " average " << stats.average();
    cerr << endl;

    Accumulator<distance_t> mind, maxd;
    for (int i = 0; i < pool.number_elements(); i++)
    {
      const Descriptors_and_Neighbours & dni = pool[i];

      int n = dni.number_neighbours();
      if (0 == n)
        continue;

      mind.extra (dni.distance_of_closest_neighbour());
      if (n > 1)
        maxd.extra (dni.distance_of_furthest_neighbour());
    }

    cerr << "Of neighbours stored, closest distances between " << mind.minval() << " and " << mind.maxval();
    if (mind.n() > 1)
      cerr << " ave " << mind.average();
    cerr << endl;

    cerr << maxd.n() << " furthest neighbours between " << maxd.minval() << " and " << maxd.maxval();
    if (maxd.n() > 1)
      cerr << " ave " << maxd.average();
    cerr << endl;

    if (nwindow)
      cerr << computations_avoided_because_of_window << " computations avoided by presence of " << nwindow << " windows\n";
  }

  if (cl.option_present('q'))
    _exit(0);

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_similarity (argc, argv);

  return rc;
}
