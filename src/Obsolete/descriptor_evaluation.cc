/*
  We have a descriptor file and need to figure out which
  descriptors do the best job of discriminating activity
*/

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <ostream>
#include <unordered_map>
#include <unordered_set>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/sparsefp.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int activity_column = -1;

static IWString fingerprint_tag;

static int show_zero_hits = 1;

static int * process_column = nullptr;

/*
  Fred sometimes doesn't want the data on multiple records.
*/

static int newline_in_output = 1;

/*
  When names have come from a descriptor computation, they may have
  underscores in them
*/

static int gsub_underscores_in_name = 0;

/*
  We can read sparse fingerprints in one of two forms. 
  They can be ascii representation,   1,3,15-20;2048
  Or they can be non-colliding counted form
*/

static int sparse_ascii_representation = 0;

/*
  We may not want to print out all the bits
*/

static int nprint = 0;

/*
  The MACCS keys start at 0, but mk0 is actually maccs key 6
*/

static int print_offset = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;

  cerr << "Ranks descriptors by their ability to differentiate activity\n";
  cerr << " -E col=n         activities are column <N> of the name\n";
  cerr << " -E file=fname    id/activity file\n";
  cerr << " -s <size>        number of molecules to be processed\n";
  cerr << " -p <number>      number of the best descriptors to print\n";
  cerr << " -c <col>         process column(s). Use '1-10' for a range\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

class Descriptor
{
  private:
    typedef float descriptor_t;

    IWString _name;
    descriptor_t * _d;

  public:
    Descriptor ();
    ~Descriptor ();

    const IWString & descriptor_name () const { return _name;}
    void set_name (const IWString & n) { _name = n;}

    int allocate_descriptor_array (int);

    int assign_value (const const_IWSubstring &, int);
};

Descriptor::Descriptor ()
{
  _d = nullptr;

  return;
}

Descriptor::~Descriptor ()
{
  if (NULL != _d)
    delete _d;

  return;
}

int
Descriptor::allocate_descriptor_array (int s)
{
  assert (s > 0);

  if (NULL != _d)
    delete _d;

  _d = new descriptor_t[s];

  return NULL != _d;
}

int
Descriptor::assign_value (const const_IWSubstring & token,
                          int r)
{
  if (! token.numeric_value (_d[r]))
    return 0;

  return 1;
}

class ID_Activity
{
  public:
    ID_Activity ();

    typedef float activity_t;

    void set_id (IWString & s) { _id = s;}
    const IWString & id () const { return _id;}

    activity_t activity () const { return _activity;}

  private:
    IWString _id;
    activity_t _activity;

};

ID_Activity::ID_Activity ()
{
  _activity = static_cast<activity_t> (0.0);

  return;
}

class EResults
{
  private:
    unsigned int _zbit;
    int _hit_1;
    int _hit_2;
    int _nothit_1;
    int _nothit_2;

    float _iab;

    float _explained_by_descriptor_split;

//  private functions

    void _default_values ();

  public:
    EResults ();
    EResults (unsigned int);

//  We rely on the copy operator being a bit-wise copy

    int  set_bit (unsigned int b);

    int  bit_hit_in_class (int c);

    int report (int, std::ostream &) const;

    float iab () const { return _iab;}

    int set_results (int, double, double);
};

EResults::EResults ()
{
  _zbit = 0;
  _default_values ();
}

EResults::EResults (unsigned int b)
{
  _zbit = b;
  _default_values ();
}

void
EResults::_default_values ()
{
  _hit_1 = 0;
  _hit_2 = 0;
  _nothit_1 = 0;
  _nothit_2 = 0;

  _iab = static_cast<float> (0.0);

  return;
}

int
EResults::bit_hit_in_class (int c)
{
  if (c1 == c)
    _hit_1++;
  else if (c2 == c)
    _hit_2++;
  else
  {
    cerr << "EResults::bit_hit_in_class: invalid class " << c << " must be " << c1 << " or " << c2 << endl;
    abort ();
  }

  return 1;
}

int
EResults::set_bit (unsigned int b)
{
  _zbit = b;

  _default_values ();

  return 1;
}


/*
  We have read in all our molecules, and now know how many there are in the whole set
*/

int 
EResults::set_results (int pool_size,
                       double e1,
                       double e2)
{
//cerr << "Bit " << _zbit << ", " << pool_size << " molecules, " << _hit_1 << " hits in class 1, " << _hit_2 << " hits in class 2\n";
  _nothit_1 = n1 - _hit_1;
  assert (_nothit_1 >= 0);

  _nothit_2 = n2 - _hit_2;
  assert (_nothit_2 >= 0);

  double T = static_cast<double> (pool_size);

  double eh = compute_entropy (_hit_1 + _hit_2, T);               // entropy of the hits
  double em = compute_entropy (_nothit_1 + _nothit_2, T);         // entropy of the misses

  double e = compute_entropy (_hit_1, T) + compute_entropy (_hit_2, T) + compute_entropy (_nothit_1, T) + compute_entropy (_nothit_2, T);

// The output variables are all float 

  _iab = (e1 + e2) + (eh + em) - e;

  assert (_iab >= 0.0);

  _explained_by_descriptor_split = _iab / (e1 + e2);

  return 1;
}

int
EResults::report (int pool_size, 
                  std::ostream & os) const
{
  double T = static_cast<double> (pool_size);

  double eh = compute_entropy (_hit_1 + _hit_2, T);               // entropy of the hits
  double em = compute_entropy (_nothit_1 + _nothit_2, T);         // entropy of the misses

  os << "Bit " << _zbit << " hit " << (_hit_1 + _hit_2) << " times. hit1 " << _hit_1 << " hit2 " << _hit_2 << ", nothit1 " << _nothit_1 << ", nothit2 " << _nothit_2 << ". EH " << static_cast<float> (eh) << " EM " << static_cast<float> (em) << " IAB " << _iab;
  if (newline_in_output)
    os << endl;

  float bits_per_compound = _iab / T;
  float descriptor_info_that_predicts_activity = _iab / (eh + em);

  os << " bits/cpd " << bits_per_compound << " explained by split " << _explained_by_descriptor_split << " predicts " << descriptor_info_that_predicts_activity << endl;

  return os.good ();
}

/*static int
eresults_comparitor (const void * p1, const void * p2)
{
  const EResults * e1 = (const EResults *) p1;
  const EResults * e2 = (const EResults *) p2;

  float iab1 = e1->iab ();
  float iab2 = e2->iab ();

  if (iab1 < iab2)
    return 1;
  if (iab1 > iab2)
    return -1;

  return 0;
}*/

static int
increment_class_counters (int c)
{
  if (c == c1)
  {
    n1++;
    return 1;
  }

  if (c == c2)
  {
    n2++;
    return 1;
  }

  if (c1 < 0)
  {
    c1 = c;
    n1 = 1;
    return 1;
  }

  if (c2 < 0)
  {
    c2 = c;
    n2 = 1;
    return 1;
  }

  cerr << "Invalid class " << c << " classes are " << c1 << " and " << c2 << endl;

  return 0;
}

/*
  We need to keep track of all the bits encountered
*/

static resizable_array_p<EResults> bits_found;

#ifdef __GNUG__
template class resizable_array_p<EResults>;
template class resizable_array_base<EResults *>;
#endif

/*
  We need a quick means of getting the results for a particular bit
*/

static std::unordered_map<unsigned int, EResults *> bit_hash;

static int
eresults_comparitor (EResults * const * v1, EResults * const * v2)
{
  const EResults * e1 = *v1;
  const EResults * e2 = *v2;

  float iab1 = e1->iab ();
  float iab2 = e2->iab ();

  if (iab1 < iab2)
    return 1;
  if (iab1 > iab2)
    return -1;

  return 0;
}

static int fingerprints_read = 0;

static int
InformationContent (std::ostream & output)
{
  int nb = bits_found.number_elements ();

  if (verbose)
  {
    cerr << "Input contains " << nb << " different bits in " << fingerprints_read << " fingerprints\n";
    cerr << "N1 = " << n1 << ", n2 = " << n2 << endl;
  }

  double e1 = compute_entropy (n1, static_cast<double> (fingerprints_read));
  double e2 = compute_entropy (n2, static_cast<double> (fingerprints_read));

  if (verbose)
    cerr << "E1 = " << static_cast<float> (e1) << " E2 = " << static_cast<float> (e2) << endl;

  Accumulator<float> acc;
  extending_resizable_array<int> bucket;
  bucket.resize (int (bits_found[0]->iab ()) + 1);

  for (int i = 0; i < nb; i++)
  {
    EResults * eri = bits_found[i];

    eri->set_results (fingerprints_read, e1, e2);
    if (verbose)
      acc.extra (eri->iab ());

    int b = static_cast<int> (eri->iab ());

    bucket[b]++;
  }

  bits_found.sort (eresults_comparitor);

  if (0 == nprint)
    nprint = nb;

  for (int i = 0; i < nprint; i++)
  {
    const EResults * eri = bits_found[i];

    eri->report (fingerprints_read, std::cout);
  }

  if (verbose)
  {
    cerr << "IAB values between " << acc.minval () << " and " << acc.maxval ();
    if (acc.n () > 1)
      cerr << " ave " << acc.average ();
    cerr << endl;

    for (int i = 0; i < bucket.number_elements (); i++)
    {
      if (bucket[i])
        cerr << bucket[i] << " bits with IAB values above " << i << endl;
    }
  }

  return output.good ();
}

template <typename T>
int
InformationContent (T & fp)
{
  fingerprints_read++;

  int c = fp.activity_class ();

  if (! increment_class_counters(c))
    return 0;

  int i = 0;
  unsigned int zbit;
  int hits;
  while (fp.next_bit_set (i, zbit, hits))
  {
    EResults * er;

    std::unordered_map<unsigned int, EResults *>::const_iterator f = bit_hash.find (zbit);
    if (f == bit_hash.end ())
    {
      er = new EResults (zbit);
      bits_found.add (er);
      bit_hash[zbit] = er; 

      f = bit_hash.find(zbit);
    }
    else
      er = (*f).second;

    int c = fp.activity_class();

    er->bit_hit_in_class (c);
  }

  return 1;
}

static int
InformationContent_sparse (iwstring_data_source & input,
                           std::ostream & output)
{
  IW_TDT tdt;
  while (tdt.next (input))
  {
    FBFingerprint_sparse fp;
    int fatal;
    if (! fp.construct_from_tdt (tdt, fatal))
    {
      if (! fatal)
        continue;

      cerr << "Invalid TDT, line " << input.lines_read () << endl;
      cerr << tdt;
      return 0;
    }

    if (! InformationContent (fp))
    {
      cerr << "Fatal error, " << input.lines_read () << " lines read\n";
      return 0;
    }
  }

  return InformationContent (output);
}

static int
InformationContent_dense (iwstring_data_source & input,
                          std::ostream & output)
{
  IW_TDT tdt;
  while (tdt.next (input))
  {
    FBFingerprint_dense fp;
    int fatal;
    if (! fp.construct_from_tdt (tdt, fatal))
    {
      cerr << "INvalid TDT, line " << input.lines_read () << endl;
      cerr << tdt;
      return 0;
    }

    if (! InformationContent(fp))
    {
      cerr << "Fatal error, " << input.lines_read () << " lines read\n";
      return 0;
    }
  }

  return InformationContent (output);
}

static int
InformationContent_sparse (const char * fname,
                          std::ostream & output)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return InformationContent_sparse (input, output);
}

static int
InformationContent_dense (const char * fname,
                          std::ostream & output)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return InformationContent_dense (input, output);
}

static int
get_activities_from_file (const IWString & fname,
                          int number_rows,
                          ID_Activity * idact)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open activity file '" << fname << "'\n";
    return 0;
  }

  return get_activities_from_file (input, number_rows, idact);
}

static int
read_descriptor_file (const const_IWSubstring & buffer,
                      ID_Activity * idact,
                      Descriptor * descriptors,
                      int r)
{
  int i = 0;
  (void) buffer.nextword (idact[r].id(), i);

  const_IWSubstring token;

  int col = 0;
  while (buffer.nextword (token, i))
  {
    if (! descriptors[i].assign_value (token, r))
    {
      cerr << "Invalid value '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_descriptor_file (iwstring_data_source & input,
                      ID_Activity * idact,
                      Descriptor * descriptors)
{
  int r = 0;   // whcih record are we processing

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! read_descriptor_file (buffer, idact, descriptors, r))
    {
      cerr << "Fatal error on line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }

    r++;
  }

  return 1;
}

static int
InformationContent (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:A:o:zT:p:yuc:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  IWString activity_file;

  if (cl.option_present ('E'))
  {
    IWString expt = cl.string_value ('E');

    if (verbose)
      cerr << "Fingerprints in '" << fingerprint_tag << "' dataitem processed\n";
  }

  if (cl.option_present ('A'))
  {
    const_IWSubstring a = cl.string_value ('A');

    if (a.starts_with ("col="))
    {
      a.remove_leading_chars (4);
      if (! a.numeric_value (activity_column) || activity_column < 1)
      {
        cerr << "INvalid activity column '" << a << "'\n";
        usage (3);
      }

      if (verbose)
        cerr << "The activity is the " << activity_column << " token in the name\n";

      activity_column--;
    }
    else if (a.starts_with ("file="))
    {
      a.remove_leading_chars (5);
      activity_file = a;

      if (verbose)
        cerr << "Will read activities from '" << activity_file << "'\n";
    }
    else
    {
      cerr << "Unrecognised -E qualifier '" << a << "'\n";
      usage (6);
    }
  }

  if (cl.option_present ('p'))
  {
    if (! cl.value ('p', nprint) || nprint < 0)
    {
      cerr << "Invalid number of bits to print (-p option)\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will print only the best " << nprint << " bits\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (1);
  }

  iwstring_data_source input (cl[0]);

  if (! input.good ())
  {
    cerr << "Cannot open '" << cl[0] << "'\n";
    return 5;
  }

  const_IWSubstring header;
  if (! input.next_record (header))
  {
    cerr << "Cannot read header record from '" << cl[0] << "'\n";
    return 6;
  }

  number_columns = header.nwords () - 1;
  if (number_columns <= 0)
  {
    cerr << "Bad header record '" << header << "'\n";
    return 6;
  }

  // This is not working, TODO: ianwatson investigate if it ever matters.
  if (cl.option_present ('c'))
  {
    process_column = new_int (number_columns);

    const_IWSubstring c = cl.string_value ('c');
    if (c.contains ('-'))
    {
    }
    else
    {
      int zcol;
      if (! c.numeric_value (zcol) || zcol < 2 || zcol > number_columns)   // check this - upper limit may be wrong
      {
        cerr << "INvalid column specification " << c << "' must be less than " << number_columns << endl;
        return 6;
      }

      if (verbose)
        cerr << "Will only process column " << zcol << endl;
      process_column[zcol - 1] = 1;
    }
  }
  else
  {
    number_descriptors = number_columns;
  }

  Descriptor * descriptors  = new Descriptor[number_descriptors];
  if (NULL == descriptors)
  {
    cerr << "Bad news, cannot allocate " << number_descriptors << " descriptors\n";
    return 7;
  }

  if (verbose)
    cerr << "Input contains " << number_descriptors << " descriptors\n";

  number_rows = input.records_remaining ();
  if (0 == number_rows)
  {
    cerr << "No records in the file\n";
    return 6;
  }

  ID_Activity * idact = new ID_Activity[number_rows];
  if (NULL == idact)
  {
    cerr << "Memory failure, " << number_rows << " rows\n";
    return 6;
  }

  for (int i = 0; i < number_descriptors; i++)
  {
    if (! descriptors[i].allocate_descriptor_array (number_records))
    {
      cerr << "Memory failure\n";
      return 11;
    }
  }

  if (! read_descriptor_file (input, idact, descriptors))
  {
    cerr << "Cannot read descriptors from '" << cl[0] << "'\n";
    return 7;
  }

  int rc;
  if (activity_file.length ())
    rc = get_activities_from_file (activity_file, number_rows, descriptors);
  else
    rc = get_activities_from_name (descriptors, number_rows, descriptors);

  if (0 == rc)
  {
    cerr << "Fatal error establishing activities\n";
    return 5;
  }
  int rc;

  if (! cl.option_present ('T'))
    rc = InformationContent_dense (cl[0], std::cout);
  else
  {
    const_IWSubstring t = cl.string_value ('T');
    if ("NC" == t)
      rc = InformationContent_sparse (cl[0], std::cout);
    else if ("sparse" == t)
    {
      sparse_ascii_representation = 1;
      rc = InformationContent_dense (cl[0], std::cout);
    }
    else
    {
      cerr << "Unrecognised fingerprint type '" << t << "'\n";
      usage (5);
    }
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = InformationContent (argc, argv);

  return rc;
}
