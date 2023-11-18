/*
  Generates difference fingerprints from an activity set
*/

#include <stdlib.h>
#include <iostream>
#include <limits>
#include <math.h>
#include <unordered_map>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/gfp.h"

using std::cerr;
using std::endl;


const char * prog_name = nullptr;

static IW_STL_Hash_Map_float activity;

static int activity_column = 1;

static int verbose = 0;

static int strip_leading_zeros = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");
static IWString activity_difference_tag("ADIFF<");
static IWString tag("NCDIF<");

/*
  We want to keep track of the number of items that get to appear in
  the output
*/

static IW_STL_Hash_Map_int times_in_output;

/*
  We can impose a requirement that at least one data point
  be above some kind of activity threshold
*/

static float activity_threshold = static_cast<float>(0.0);
static int   number_below_threshold = 0;

/*
  We can constrain which differences we examine
  We can constrain based on the observed activity difference, or the size of
  the vector produced
*/

static float min_activity_difference = static_cast<float>(0.0);
static float max_activity_difference = std::numeric_limits<float>::max();

static int max_nset = std::numeric_limits<int>::max();

static unsigned int pairs_skipped_for_too_many_bits_in_differencefp = 0;

static float max_distance_between_fingerprints = static_cast<float>(1.0);

static unsigned int pairs_skipped_for_distance_too_far = 0;

static Accumulator<float> distance_accumulator;
static Accumulator<float> distances_written;

static int fingerprints_written = 0;

static Report_Progress report_progress;

static int include_activity_difference_in_name_field = 1;

static int include_distance_in_output = 1;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Generates difference fingerprints\n";
  cerr << " -A <fname>     activity data\n";
  cerr << " -a <cutoff>    activity threshold - at least one member of a pair must have > cutoff activity\n";
  cerr << " -d <diff>      min activity difference for forming a pair\n";
  cerr << " -D <diff>      max activity difference for forming a pair\n";
  cerr << " -s <nbits>     max allowable bits in difference fingerprint\n";
  cerr << " -T <dist>      only compute difference fingerprints for FP's closer than <dist>\n";
  cerr << " -r <number>    report progress every <number> fingerprints processed\n";
  cerr << " -n <number>    output the <number> shortest distances\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Fingerprint_and_Activity : public IW_General_Fingerprint
{
  private:
    float _activity;
    IWString _smiles;

  public:
    Fingerprint_and_Activity();

    void set_activity (float s) { _activity = s;}
    float activity () const { return _activity;}

    void set_smiles (const IWString & s) { _smiles = s;}
    const IWString & smiles() const { return _smiles;}
};

#ifdef __GNUG__
template class resizable_array_p<Fingerprint_and_Activity>;
template class resizable_array_base<Fingerprint_and_Activity *>;
#endif

Fingerprint_and_Activity::Fingerprint_and_Activity()
{
  _activity = static_cast<float>(0.0);
}

typedef IW_STL_Hash_Map<IWString, Fingerprint_and_Activity *> ID_to_FP;

class Fingerprint_and_Activity_Comparator
{
  private:
  public:
    int operator () (const Fingerprint_and_Activity *, const Fingerprint_and_Activity *) const;
};

int
Fingerprint_and_Activity_Comparator::operator () (const Fingerprint_and_Activity * fpa1,
                                                  const Fingerprint_and_Activity * fpa2) const
{
  float a1 = fpa1->activity();
  float a2 = fpa2->activity();

  if (a1 < a2)
    return -1;

  if (a2 > a2)
    return 1;

  return 0;

}

/*
  Sometimes we might want to just process all the shortest distances in a file
*/

class Distance_Result
{
  private:
    float _d;
    Fingerprint_and_Activity * _fp1;
    Fingerprint_and_Activity * _fp2;

  public:
    Distance_Result();

    Fingerprint_and_Activity * fp1 () const { return _fp1;}
    Fingerprint_and_Activity * fp2 () const { return _fp2;}

    float distance() const { return _d;}

    void set (Fingerprint_and_Activity * f1, Fingerprint_and_Activity * f2, float d) { _fp1 = f1; _fp2 = f2; _d = d;}
};

Distance_Result::Distance_Result()
{
  _d = static_cast<float>(0.0);

  _fp1 = nullptr;
  _fp2 = nullptr;

  return;
}

/*
  In a couple of different contexts we need to keep track of
  the N closest distances
*/

class Set_of_Shortest_Distances
{
  private:
    Distance_Result ** _dr;
    int _max_store;
    int _nstored;

//  private functions

    void _check_shortest_list_ok () const;
    int _binary_search (float d, int nsearch) const;

  public:
    Set_of_Shortest_Distances();
    ~Set_of_Shortest_Distances();

    int initialise (int);

    int extra (Fingerprint_and_Activity * fp1, Fingerprint_and_Activity * fp2,
                                float activity_difference);
    
    int nstored() const { return _nstored;}

    const Distance_Result * const* distance_result() const { return _dr;}
};

Set_of_Shortest_Distances::Set_of_Shortest_Distances()
{
  _dr = nullptr;
  _max_store = 0;
  _nstored = 0;

  return;
}

Set_of_Shortest_Distances::~Set_of_Shortest_Distances()
{
  assert (_nstored >= 0);

  if (NULL != _dr)
  {
    for (int i = 0; 0 < _nstored; i++)
    {
      delete _dr[i];
    }

    delete [] _dr;
  }

  _dr = nullptr;

  _nstored = -3;

  return;
}

int
Set_of_Shortest_Distances::initialise (int s)
{
  assert (0 == _max_store);
  assert (s > 0);

  _dr = new Distance_Result * [s];

  for (int i = 0; i < s; i++)
  {
    _dr[i] = new Distance_Result;
  }

  _max_store = s;

  _nstored = 0;

  return 1;
}

/*
  Identify the insertion point for d.
  We return the index at which it should be inserted
*/

int
Set_of_Shortest_Distances::_binary_search (float d,
                                         int nsearch) const
{
  float dright = _dr[nsearch - 1]->distance();   // check RHS first, expected to be most powerful filter

  if (d >= dright)
    return nsearch;

  float dleft = _dr[0]->distance();
  if (d <= dleft)
    return 0;

  int left = 0;
  int right = nsearch - 1;

  while (right > left)
  {
    int middle = (left + right) / 2;
    if (middle == left)
    {
      left = right;
      break;
    }

    float dmid = _dr[middle]->distance ();
//    cerr << "Left = " << left << " d = " << _neighbours[left]->distance () << " middle " << middle << " d = " << _neighbours[middle].distance () << " right " << right << " d = " << _neighbours[right].distance () << endl;
    if (d < dmid)
      right = middle;
    else if (d > dmid)
      left = middle;
    else
    {
      left = middle;
      break;
    }
  }

  return left;
}

void
Set_of_Shortest_Distances::_check_shortest_list_ok () const
{
  int n = _max_store;

  if (_nstored < _max_store)
    n = _nstored;

  int rc = 1;

  for (int i = 1; i < n; i++)
  {
    if (_dr[i]->distance() < _dr[i - 1]->distance())
    {
      cerr << "Shortest distance array out of order, i = " << i << " prev " << _dr[i-1]->distance() << " compare " << _dr[i]->distance() << ", n = " << n << endl;
      rc = 0;
    }
  }

  if (0 == rc)
    abort();

  return;
}

int
Set_of_Shortest_Distances::extra (Fingerprint_and_Activity * fp1,
                                Fingerprint_and_Activity * fp2,
                                float activity_difference)
{
  float d = fp1->distance(*fp2);

  if (_nstored < 2)
  {
    if (0 == _nstored)
      _dr[0]->set(fp1, fp2, d);
    else
      _dr[1]->set(fp1, fp2, d);
    _nstored++;
    return 1;
  }

  int nsearch;
  if (_nstored < _max_store)
    nsearch = _nstored;
  else
    nsearch = _max_store;

  int ndx = _binary_search(d, nsearch);

#ifdef DEBUG_INSERT_INTO_SHORTEST_LIST
  cerr << "Insert " << d << " between " << _dr[0]->distance() << " and " << _dr[nsearch - 1]->distance() << ", n = " << nsearch << ", insert at " << ndx << endl;
  for (int i = 0; i < _nstored; i++)
  {
    cerr << " i = " << i << " dist " << _dr[i]->distance() << endl;
  }
#endif

  if (ndx == _max_store)   // array is full and we are larger than largest stored
    return 0;

  if (ndx == nsearch)   // array not full, we are larger than largest stored
  {
    assert (nsearch < _max_store);
    _dr[nsearch]->set(fp1, fp2, d);
    _nstored++;
    return 0;
  }

// Gets inserted into the list. This is tricky. We need to identify the
// object that will be lost and put the new contents into it.

  assert (_nstored <= _max_store);

  int istart;
  if (_nstored < _max_store)
    istart = _nstored;
  else
    istart = _nstored - 1;

  Distance_Result * holds_new_value = _dr[istart];

  holds_new_value->set(fp1, fp2, d);

#ifdef DEBUG_INSERT_INTO_SHORTEST_LIST
  cerr << "istart " << istart << endl;
#endif

  for (int i = istart; i > ndx; i--)
  {
    _dr[i] = _dr[i - 1];
  }

  _dr[ndx] = holds_new_value;

  if (_nstored < _max_store)
    _nstored++;

#ifdef DEBUG_INSERT_INTO_SHORTEST_LIST
  cerr << "After insertion\n";
  for (int i = 0; i < _nstored; i++)
  {
    cerr << " i = " << i << " dist " << _dr[i]->distance() << endl;
  }
#endif

  _check_shortest_list_ok();

  return 1;
}

static int
process_fingerprint (IW_TDT & tdt,
                     const IW_STL_Hash_Map_float & activity,
                     resizable_array_p<Fingerprint_and_Activity> & fps)
{
  Fingerprint_and_Activity * fp = new Fingerprint_and_Activity;

  int fatal;
  if (! fp->construct_from_tdt(tdt, fatal))
  {
    delete fp;

    if (fatal)
      return 0;

    return 1;
  }

  if (fp->id().contains(' '))
  {
    IWString id = fp->id();

    id.truncate_at_first(' ');

    fp->set_id(id);
  }

  const IWString & id = fp->id();

  IW_STL_Hash_Map_float::const_iterator f = activity.find(id);

  if (f == activity.end())
  {
    if (id.nwords() > 1)
    {
      IWString tmp(id);
      tmp.truncate_at_first(' ');
      f = activity.find(tmp);
    }

    if (f == activity.end())
    {
      cerr << "No activity data for '" << id << "'\n";
      delete fp;
      return 0;
    }
  }

  fp->set_activity((*f).second);

  IWString smiles;

  if (tdt.dataitem_value(smiles_tag, smiles))
    fp->set_smiles(smiles);

  fps.add(fp);

  return 1;
}

static int
read_fingerprints (iwstring_data_source & input,
                   const IW_STL_Hash_Map_float & activity,
                   resizable_array_p<Fingerprint_and_Activity> & fps)
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    if (! process_fingerprint (tdt, activity, fps))
      return 0;
  }

  return fps.number_elements();
}

static int
read_fingerprints (const char * fname,
                   const IW_STL_Hash_Map_float & activity,
                   resizable_array_p<Fingerprint_and_Activity> & fps)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  return read_fingerprints (input, activity, fps);
}

static int
read_activity_data_record (const const_IWSubstring & buffer,
                           IW_STL_Hash_Map_float & activity)
{
  IWString id, token;
  int i = 0;

  char separator = ' ';

  if (buffer.nextword(id, i) && buffer.nextword(token, i))
    ;
  else if (buffer.nextword(id, i, '\t') && buffer.nextword(token, i, '\t'))
    separator = '\t';
  else if (buffer.nextword(id, i, ',') && buffer.nextword(token, i, ','))
    separator = ',';
  else
  {
    cerr << "Cannot separate into identifier and activity '" << buffer << "'\n";
    return 0;
  }

  if (1 != activity_column)
  {
    if (! buffer.word(activity_column, token, separator))
    {
      cerr << "Cannot extract column '" << (activity_column + 1) << " from record\n";
      return 0;
    }
  }

  float a;
  if (token.numeric_value(a))
    ;
  else if (0 == activity.size())   // header record
    ;
  else
    cerr << "Warning, non numeric activity value '" << token << "', id '" << id << "'\n";
  
  if (strip_leading_zeros)
    id.remove_leading_chars('0');

  activity[id] = a;

//cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token << "'\n";

  return 1;
}

static int
read_activity_data (iwstring_data_source & input,
                    IW_STL_Hash_Map_float & activity)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (! read_activity_data_record(buffer, activity))
    {
      cerr << "Cannot read activity data, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return activity.size();
}

static int
read_activity_data (const const_IWSubstring & fname,
                    IW_STL_Hash_Map_float & activity)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return read_activity_data (input, activity);
}

#ifdef NOT_BEING_USED
static int
do_svml_output (const Fingerprint_and_Activity & fp1,
                const Fingerprint_and_Activity & fp2,
                float d,
                Sparse_Fingerprint_Creator & sfc,
                IWString_and_File_Descriptor & output)
{
  output << d;

  sfc.write_in_svml_form(output);

  output << " # " << fp1.id() << '.' << fp2.id() << '\n';

  return 1;
}
#endif

static int
do_gfp_output (const Fingerprint_and_Activity & fp1,
               const Fingerprint_and_Activity & fp2,
               float activity_difference,
               float distance,
               Sparse_Fingerprint_Creator & sfc,
               IWString_and_File_Descriptor & output)
{
  distances_written.extra(distance);

  if (fp1.smiles().length() && fp2.smiles().length())
    output << smiles_tag << fp1.smiles() << '.' << fp2.smiles() << ">\n";
  
  output << identifier_tag << fp1.id() << " - " << fp2.id();

  if (include_activity_difference_in_name_field)
    output << ' ' << activity_difference;

  output << ">\n";

  if (include_distance_in_output)
    output << distance_tag << distance << ">\n";

  output << activity_difference_tag << activity_difference << ">\n";

  sfc.write_fingerprint(tag, output);

  output << "|\n";

  return 1;
}

static int
do_output (Fingerprint_and_Activity & fp1,
           Fingerprint_and_Activity & fp2,
           float activity_difference,
           float distance,
           Sparse_Fingerprint_Creator & sfc,
           IWString_and_File_Descriptor & output)
{
#ifdef ECHO_FP_BEFORE_WRITING
  const Sparse_Fingerprint_Creator::FPHash & bits_found = sfc.bits_found();

  cerr << "Writing sparse fingerprint with " << bits_found.size() << " bits\n";
  for (Sparse_Fingerprint_Creator::FPHash::const_iterator i = bits_found.begin(); i != bits_found.end(); ++i)
  {
    cerr << "Found " << (*i).second << " occurrences of bit " << (*i).first << endl;
  }
#endif

  int rc = do_gfp_output(fp1, fp2, activity_difference, distance, sfc, output);

  if (0 == rc)
    return 0;

  fingerprints_written++;

  output.write_if_buffer_holds_more_than(32768);

  return rc;
}

/*
  The difference is fp1-fp2
*/

static int
gfp_difference_fingerprint (Fingerprint_and_Activity & fp1,
                            Fingerprint_and_Activity & fp2,
                            float activity_difference,
                            IWString_and_File_Descriptor & output)
{
  const Sparse_Fingerprint & sfp1 = fp1.sparse_fingerprint(0);
  const Sparse_Fingerprint & sfp2 = fp2.sparse_fingerprint(0);

  Sparse_Fingerprint_Creator sfc;

  int i1 = 0;
  unsigned int b1;
  int c1;

  int i2 = 0;
  unsigned int b2;
  int c2;
  (void) sfp1.next_bit_set(i1, b1, c1);
  (void) sfp2.next_bit_set(i2, b2, c2);

  int bic = 0;

  std::unordered_map<unsigned int, int> times_set;

//cerr << "Processing bits with " << sfp1.nbits() << " and " << sfp2.nbits() << " bits\n";
  while (1)
  {
    if (b2 < b1)
    {
      sfc.hit_bit(b2, 128-c2);
      if (times_set[b2])
        cerr << "2 Duplicate set " << b2 << endl;
      times_set[b2]++;
//    cerr << "1 " << b1 << " " << (128-c1) << endl;
      if (! sfp2.next_bit_set(i2, b2, c2))
        break;
    }

    if (b1 < b2)
    {
      sfc.hit_bit(b1, 128+c1);
      if (times_set[b1])
        cerr << "1 Duplicate set " << b1 << endl;
      times_set[b1]++;
      if (! sfp1.next_bit_set(i1, b1, c1))
        break;
    }

    if (b1 == b2)
    {
      if (times_set[b1] > 0)
        cerr << "= Duplicate set " << b1 << endl;
      times_set[b1]++;
      if (c1 > c2)
      {
        sfc.hit_bit(b1, 128 + c1 - c2);
//      cerr << "= " << b1 << ' ' << (128 + c1 - c2) << endl;
        bic += c2;
      }
      else
      {
        sfc.hit_bit(b1, 128 + c1 - c2);
//      cerr << "= " << b1 << ' ' << (128 + c1 - c2) << endl;
        bic += c1;
      }
      if (! sfp1.next_bit_set(i1, b1, c1))
        break;
      if (! sfp2.next_bit_set(i2, b2, c2))
        break;
    }
  }

  while (sfp1.next_bit_set(i1, b1, c1))
  {
    sfc.hit_bit(b1, 128+c1);
//  cerr << "F " << b1 << (128+c1) << endl;
  }

  while (sfp2.next_bit_set(i2, b2, c2))
  {
    sfc.hit_bit(b2, 128-c2);
//  cerr << "G " << b2 << (128+c2) << endl;
  }

  float tanimoto = static_cast<float>(1.0) - static_cast<float>(bic) / static_cast<float>(sfp1.nset() + sfp2.nset() - bic);
//cerr << "bic " << bic << " tanimoto " << tanimoto << endl;
//cerr << "Computed tanimoto " << (1.0 -sfp1.tanimoto(sfp2)) << endl;

  distance_accumulator.extra(tanimoto);

  if (tanimoto > max_distance_between_fingerprints)
  {
    pairs_skipped_for_distance_too_far++;
    return 1;
  }

  int bits_in_difference_fingerprint = 0;

  const Sparse_Fingerprint_Creator::FPHash bits_found = sfc.bits_found();

  for (Sparse_Fingerprint_Creator::FPHash::const_iterator i = bits_found.begin(); i != bits_found.end(); ++i)
  {
    int c = (*i).second;
    if (c > 128)
      bits_in_difference_fingerprint += c - 128;
    else if (c < 128)
      bits_in_difference_fingerprint += 128 - c;
  }

  if (bits_in_difference_fingerprint > max_nset)
  {
    pairs_skipped_for_too_many_bits_in_differencefp++;
    return 1;
  }

  if (verbose)
  {
    times_in_output[fp1.id()]++;
    times_in_output[fp2.id()]++;
  }

  return do_output (fp1, fp2, activity_difference, tanimoto, sfc, output);
}

/*static void
insert_into_shortest_list (Fingerprint_and_Activity * fp1,
                           Fingerprint_and_Activity * fp2,
                           int nsmall,
                           int & nstored,
                           Distance_Result ** dr)
{
  float d = fp1->distance(*fp2);

  if (nstored < 2)
  {
    if (0 == nstored)
      dr[0]->set(fp1, fp2, d);
    else
      dr[1]->set(fp1, fp2, d);
    nstored++;
    return;
  }

  int nsearch;
  if (nstored < nsmall)
    nsearch = nstored;
  else
    nsearch = nsmall;

  int ndx = binary_search(d, nsearch);

#ifdef DEBUG_INSERT_INTO_SHORTEST_LIST
  cerr << "Insert " << d << " between " << dr[0]->distance() << " and " << dr[nsearch - 1]->distance() << ", n = " << nsearch << ", insert at " << ndx << endl;
  for (int i = 0; i < nstored; i++)
  {
    cerr << " i = " << i << " dist " << dr[i]->distance() << endl;
  }
#endif

  if (ndx == nsmall)   // array is full and we are larger than largest stored
    return;

  if (ndx == nsearch)   // array not full, we are larger than largest stored
  {
    assert (nsearch < nsmall);
    dr[nsearch]->set(fp1, fp2, d);
    nstored++;
    return;
  }

// Gets inserted into the list. This is tricky. We need to identify the
// object that will be lost and put the new contents into it.

  assert (nstored <= nsmall);

  Distance_Result * holds_new_value;

  int istart;
  if (nstored < nsmall)
  {
    istart = nstored;
    holds_new_value = dr[nstored];
  }
  else
  {
    istart = nstored - 1;
    holds_new_value = dr[nstored - 1];
  }

  holds_new_value->set(fp1, fp2, d);

#ifdef DEBUG_INSERT_INTO_SHORTEST_LIST
  cerr << "istart " << istart << endl;
#endif

  for (int i = istart; i > ndx; i--)
  {
    dr[i] = dr[i - 1];
  }

  dr[ndx] = holds_new_value;

  if (nstored < nsmall)
    nstored++;

#ifdef DEBUG_INSERT_INTO_SHORTEST_LIST
  cerr << "After insertion\n";
  for (int i = 0; i < nstored; i++)
  {
    cerr << " i = " << i << " dist " << dr[i]->distance() << endl;
  }
#endif

  return;
}*/

static int
gfp_difference_fingerprint (const resizable_array_p<Fingerprint_and_Activity> & fps,
                            int nsmall,
                            IWString_and_File_Descriptor & output)
{
  Set_of_Shortest_Distances sosd;

  sosd.initialise(nsmall);

  int n = fps.number_elements();

  for (int i = 0; i < n; i++)
  {
    Fingerprint_and_Activity * fpi = fps[i];

    float ai = fpi->activity();

    if (report_progress())
      cerr << "Processed " << i << " fingerprints\n";

    int jstart = i + 1;
    if (number_below_threshold > jstart)
      jstart = number_below_threshold;

    for (int j = jstart; j < n; j++)
    {
      Fingerprint_and_Activity * fpj = fps[j];

      float aj = fpj->activity();

      float d = aj - ai;

      if (d < min_activity_difference)
        continue;

      if (d > max_activity_difference)
        break;

      sosd.extra(fpi, fpj, d);
    }
  }

  const Distance_Result * const * dr = sosd.distance_result();

  for (int i = 0; i < sosd.nstored(); i++)
  {
    const Distance_Result * dri = dr[i];

    Fingerprint_and_Activity * fp1 = dri->fp1();
    Fingerprint_and_Activity * fp2 = dri->fp2();

    float a1 = fp1->activity();
    float a2 = fp2->activity();

    gfp_difference_fingerprint(*fp1, *fp2, fabs(a1 - a2), output);
  }

  return 1;
}

static int
gfp_difference_fingerprint (const resizable_array_p<Fingerprint_and_Activity> & fps,
                            IWString_and_File_Descriptor & output)
{
  int nfp = fps.number_elements();

  for (int i = 0; i < nfp; i++)
  {
    Fingerprint_and_Activity * fpi = fps[i];

    float ai = fpi->activity();

    if (report_progress())
      cerr << "Processed " << i << " fingerprints, wrote " << fingerprints_written << " difference fingerprints\n";

    int jstart = i + 1;
    if (number_below_threshold > jstart)
      jstart = number_below_threshold;

    for (int j = jstart; j < nfp; j++)
    {
      Fingerprint_and_Activity * fpj = fps[j];

      float aj = fpj->activity();

      float d = aj - ai;

      if (d < min_activity_difference)
        continue;

      if (d > max_activity_difference)
        break;

      gfp_difference_fingerprint (*fpi, *fpj, d, output);
    }
  }

  return 1;
}

/*
  The list is assumed sorted by activity
*/

static int
compute_number_below_threshold (const resizable_array_p<Fingerprint_and_Activity> & fps,
                                float c)
{
  int n = fps.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (fps[i]->activity() > c)
    {
      return i;
    }
  }

  return n;    // all are below the threshold
}

static int
gfp_difference_fingerprint (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:a:n:d:D:s:T:r:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('a'))
  {
    if (! cl.value('a', activity_threshold))
    {
      cerr << "The activity threshold value (-a) must be a valid floating point number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Activity activity threshold set to " << activity_threshold << endl;
  }

  if (cl.option_present('d'))
  {
    if (! cl.value('d', min_activity_difference) || min_activity_difference < static_cast<float>(0.0))
    {
      cerr << "The minimum activity difference option (-d) must be a valid numeric\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only compare two fingerprints if the activity difference is at least " << min_activity_difference << endl;
  }

  if (cl.option_present('D'))
  {
    if (! cl.value('D', max_activity_difference) || max_activity_difference < static_cast<float>(0.0))
    {
      cerr << "The maximum activity difference option (-D) must be a valid numeric\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only compare two fingerprints if the activity difference is less than " << max_activity_difference << endl;

    if (cl.option_present('d') && max_activity_difference < min_activity_difference)
    {
      cerr << "The max activity difference (-D) must be compatible with the min activity difference (-d)\n";
      return 2;
    }
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', max_nset) || max_nset < 1)
    {
      cerr << "The maximum bits in difference fingerprint (-s) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard difference fingerprints with more than " << max_nset << " bits set\n";
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', max_distance_between_fingerprints) || max_distance_between_fingerprints <= 0.0 || max_distance_between_fingerprints > 1.0)
    {
      cerr << "The max distance between fingerprints option (-T) must be a valid distance\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only produce a difference fingerprint if two fingerprints are within " << max_distance_between_fingerprints << endl;
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "The report progress option (-r) must be a whole +ve number\n";
      usage(3);
    }
  }

  if (cl.option_present('A'))
  {
    const_IWSubstring a = cl.string_value('A');

    if (! read_activity_data(a, activity))
    {
      cerr << "Cannot read activity data from '" << a << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Read " << activity.size() << " activity values from '" << a << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, can only handle one input file\n";
    return 3;
  }

  resizable_array_p<Fingerprint_and_Activity> fps;

  if (! read_fingerprints(cl[0], activity, fps))
  {
    cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
    return 3;
  }

  int nfp = fps.number_elements();

  if (nfp < 2)
  {
    cerr << "Must have more than one fingerprint to process\n";
    return 3;
  }

  if (1 != number_sparse_fingerprints())
  {
    cerr << "Tool is set up to process just one sparse fingerprint\n";
    return 2;
  }

  Fingerprint_and_Activity_Comparator fpac;

  fps.iwqsort (fpac);

  if (verbose)
    cerr << "Read " << nfp << " fingerprints from '" << cl[0] << "'\n";

  if (cl.option_present('a'))
  {
    number_below_threshold = compute_number_below_threshold(fps, activity_threshold);
    if (number_below_threshold <= 0 || number_below_threshold >= nfp)
    {
      cerr << "No fingerprints on one side of threshold " << activity_threshold << ", cannot continue\n";
      return 3;
    }
  }

  if (verbose)
  {
    Accumulator<float> acc;

    for (int i = 0; i < nfp; i++)
    {
      acc.extra(fps[i]->activity());
    }

    cerr << "Activity values between " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;
  }

  IWString_and_File_Descriptor output(1);

  int rc;

  if (cl.option_present('n'))
  {
    int nsmall;
    if (! cl.value('n', nsmall))
    {
      cerr << "The number of shortest distances to process (-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Output will consist of the " << nsmall << " shortest distances\n";

    rc = gfp_difference_fingerprint (fps, nsmall, output);
  }
  else
    rc = gfp_difference_fingerprint (fps, output);

  if (0 == rc)
  {
    cerr << "Fatal error generating difference fingerprints\n";
    return 3;
  }

  output.flush();

  if (verbose)
  {
    cerr << "Wrote " << fingerprints_written << " fingerprints\n";
    if (pairs_skipped_for_distance_too_far)
      cerr << pairs_skipped_for_distance_too_far << " pairs skipped for distance too far\n";
    if (pairs_skipped_for_too_many_bits_in_differencefp > 0)
      cerr << pairs_skipped_for_too_many_bits_in_differencefp << " pairs skipped for too many bits set\n";

    cerr << times_in_output.size() << " of " << fps.number_elements() << " items, " << (static_cast<float>(times_in_output.size()) / static_cast<float>(fps.number_elements())) << " were output as part of a difference fingerprint\n";

    if (distance_accumulator.n() > 1)
      cerr << distance_accumulator.n() << " distances between " << distance_accumulator.minval() << " and " << distance_accumulator.maxval() << " ave " << static_cast<float>(distance_accumulator.average()) << endl;
    if (distances_written.n() > 1)
      cerr << distances_written.n() << " distances written between " << distances_written.minval() << " and " << distances_written.maxval() << " ave " << static_cast<float>(distances_written.average()) << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_difference_fingerprint(argc, argv);

  return rc;
}
