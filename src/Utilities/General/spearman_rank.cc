/*
  Computes Spearman Rank coefficients
*/

#include <iostream>
#include <math.h>
using std::cerr;
using std::endl;


#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

const char * prog_name = NULL;

static int verbose = 0;

static int strip_leading_zeros_from_identifiers = 0;

static int header_records_to_skip = 0;

static int is_descriptor_file = 0;    // not implemented

static IWString * descriptor_name = NULL;

typedef IW_STL_Hash_Map_float ID_Activity_Hash;

static int ignore_missing_identifiers = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes Spearman Rank coefficients\n";
  cerr << " -e <col>       column for experimental (measured) values\n";
  cerr << " -E <file>      activities are in a different file <file>\n";
  cerr << " -z             strip leading 0's from identifiers when using -E\n";
  cerr << " -p <col>       column for predicted values\n";
//cerr << " -s <number>    skip <number> records at the top of each file\n";
  cerr << " -y             if data not present in both sets, ignore missing values\n";
//cerr << " -j             treat as a descriptor file\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class ID_Activity
{
  private:
    const IWString _id;
    float _activity;

  public:
    ID_Activity(const IWString & i, float a) : _id(i), _activity(a) {}

    const IWString & id () const { return _id;}

    void set_activity (float a) { _activity = a;}
    float activity () const { return _activity;}
};

class Set_of_Numbers
{
  private:
    int _column;

    IWString _column_header;

    resizable_array_p<ID_Activity> _id_act;

    Accumulator<float> _acc;

  public:
    Set_of_Numbers();

    void set_column (int s) { _column = s;}
    void set_column_header (const IWString & s) { _column_header = s;}

    int n () const { return _id_act.number_elements();}

    int column() const { return _column;}

    const IWString & column_header() const { return _column_header;}

    int echo (std::ostream &) const;

    void extra (const IWString &, float);

    int ids_present (IW_STL_Hash_Set &) const;

    int build (const const_IWSubstring & buffer);

    int build_from_file (const char * fname);
    int build_from_file (iwstring_data_source &);

    void assign_ranks();

    float spearman_rank (const Set_of_Numbers & rhs) const;
};

Set_of_Numbers::Set_of_Numbers()
{
  _column = -1;

  return;
}

class Float_Comparator_Larger
{
  private:
  public:
    int operator() (float, float) const;
};

int
Float_Comparator_Larger::operator() (float i1, float i2) const
{
  if (i1 < i2)
    return -1;

  if (i1 > i2)
    return 1;

  return 0;
}

int
Set_of_Numbers::ids_present (IW_STL_Hash_Set & ids_present_in_expt) const
{
  int n = _id_act.number_elements();

  for (int i = 0; i < n; i++)
  {
    const IWString & id = _id_act[i]->id();

    ids_present_in_expt.insert(id);
  }

  return n;
}

int
Set_of_Numbers::build (const const_IWSubstring & buffer)
{
  IWString id(buffer);
  id.truncate_at_first(' ');

  const_IWSubstring token;

  if (! buffer.word(_column, token))
  {
    cerr << "Set_of_Numbers::build:cannot extract column " << (_column + 1) << endl;
    return 0;
  }

  float v;
  if (token.numeric_value(v))   // great
    ;
  else if (0 == _acc.n())    // first record non-numeric OK
  {
    _column_header = token;
    return 1;
  }
  else
  {
    cerr << "Set_of_Numbers::build:invalid numeric '" << token << "'\n";
    return 0;
  }

  extra(id, v);

  return 1;
}

void
Set_of_Numbers::extra (const IWString & s,
                       float a)
{
  ID_Activity * t = new ID_Activity(s, a);

  _acc.extra(a);

  _id_act.add(t);

  return;
}

int
Set_of_Numbers::build_from_file (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return build_from_file(input);
}

int
Set_of_Numbers::build_from_file (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! build (buffer))
    {
      cerr << "Invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_act.number_elements();
}

class ID_Activity_Comparator
{
  private:
  public:
    int operator () (const ID_Activity *, const ID_Activity *) const;
};

int
ID_Activity_Comparator::operator () (const ID_Activity * ida1,
                                     const ID_Activity * ida2) const
{
  float a1 = ida1->activity();
  float a2 = ida2->activity();

  if (a1 < a2)
    return -1;

  if (a1 > a2)
    return 1;

  return 0;
}

void
Set_of_Numbers::assign_ranks()
{

// Add something to the array to make processing easier

  ID_Activity * t = new ID_Activity("", _acc.maxval() + 3.14e+05);

  _id_act.add(t);

  ID_Activity_Comparator idac;

  _id_act.iwqsort (idac);

  int n = _id_act.number_elements();

  float aprev = _id_act[0]->activity();

  int items_in_current_group = 1;

  for (int i = 1; i < n; i++)
  {
    float ai = _id_act[i]->activity();

    if (ai == aprev)
    {
      items_in_current_group++;
      continue;
    }

    if (1 == items_in_current_group)
    {
      _id_act[i - 1]->set_activity(static_cast<float>(i - 1));
      aprev = ai;
      continue;
    }

    float group = (static_cast<float>(i - 1) + static_cast<float>(i - items_in_current_group)) * 0.5F;
    for (int j = i - items_in_current_group; j < i; j++)
    {
      _id_act[j]->set_activity(group);
    }

    aprev = ai;
    items_in_current_group = 1;
  }

  assert (1 == items_in_current_group);

  _id_act.pop();

  return;
}

float
Set_of_Numbers::spearman_rank (const Set_of_Numbers & rhs) const
{
  float sumx2 = 0.0;
  float sumy2 = 0.0;
  float sumxy = 0.0;
  float sumx = 0.0;
  float sumy = 0.0;

  ID_Activity_Hash rhs_activity;

  int nrhs = rhs._id_act.number_elements();

  for (int i = 0; i < nrhs; i++)
  {
    const ID_Activity * rhsi = rhs._id_act[i];

    const IWString & id = rhsi->id();
    float a = rhsi->activity();

    rhs_activity[id] = a;
  }

  int missing_identifiers = 0;

  int n = _id_act.number_elements();

  for (int i = 0; i < n; i++)
  {
    const ID_Activity * lhsi = _id_act[i];

    const IWString & id = lhsi->id();

    ID_Activity_Hash::const_iterator f = rhs_activity.find(id);

    if (f == rhs_activity.end())
    {
      if (ignore_missing_identifiers)
      {
        missing_identifiers++;
        continue;
      }
      cerr << "Set_of_Numbers::spearman_rank:no activity data for '" << id << "'\n";
      return 0;
    }

    float x = lhsi->activity();
    float y = (*f).second;

    sumx += x;
    sumy += y;

    sumx2 += x * x;
    sumy2 += y * y;

    sumxy += x * y;
  }

  n = n - missing_identifiers;

  float d1 = sqrt(n * sumx2 - sumx * sumx);
  if (0.0f == d1)
    return 0.0f;

  float d2 = sqrt(n * sumy2 - sumy * sumy);
  if (0.0f == d2)
    return 0.0f;

  float rc = (static_cast<float>(n) * sumxy - sumx * sumy) / 
              d1 / d2;

  return rc;
}

int
Set_of_Numbers::echo (std::ostream & os) const
{
  int n = _id_act.number_elements();

  os << "Set_of_Numbers:echo: " << n << endl;

  for (int i = 0; i < n; i++)
  {
    os << " i = " << i << ' ' << _id_act[i]->activity() << endl;
  }

  return 1;
}

static int
read_predicted_data_values_record (const const_IWSubstring & buffer,
                            const IW_STL_Hash_Set & ids_present_in_expt,
                            Set_of_Numbers * zdata,
                            int npred,
                            Set_of_Numbers & expt)
{
  if (ids_present_in_expt.size() > 0)   // already read from separate file
  {
    IWString id(buffer);
    id.truncate_at_first(' ');
    if (strip_leading_zeros_from_identifiers)
      id.remove_leading_chars('0');

    IW_STL_Hash_Set::const_iterator f = ids_present_in_expt.find(id);

    if (f == ids_present_in_expt.end())
    {
      cerr << "NO data for '" << id << "'\n";
      return 0;
    }
  }
  else if (! expt.build(buffer))
  {
    cerr << "Cannot read experimental value\n";
    return 0;
  }

  for (int i = 0; i < npred; i++)
  {
    if (! zdata[i].build(buffer))
    {
      cerr << "Cannot process\n";
      return 0;
    }
  }

  return 1;
}

static int
store_descriptor_names (const const_IWSubstring & buffer)
{
  descriptor_name = new IWString[buffer.nwords()];

  const_IWSubstring token;

  int i = 0;

  for (int ndx = 0; buffer.nextword(token, i); ndx++)
  {
    descriptor_name[ndx] = token;
  }

  return 1;
}

static int
read_predicted_data_values (iwstring_data_source & input,
                            const IW_STL_Hash_Set & ids_present_in_expt,
                            Set_of_Numbers * zdata,
                            int npred,
                            Set_of_Numbers & expt)
{
  const_IWSubstring buffer;
  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (0 == i && is_descriptor_file)
      store_descriptor_names(buffer);

    input.next_record(buffer);
  }

  while (input.next_record(buffer))
  {
    if (! read_predicted_data_values_record (buffer, ids_present_in_expt, zdata, npred, expt))
    {
      cerr << "Cannot read predicted data record '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_predicted_data_values (const char * fname,
                            const IW_STL_Hash_Set & ids_present_in_expt,
                            Set_of_Numbers * zdata,
                            int npred,
                            Set_of_Numbers & expt)
{  
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_predicted_data_values (input, ids_present_in_expt, zdata, npred, expt);
}

/*static int
read_id_activity_hash_record (const const_IWSubstring & buffer,
                              IW_STL_Hash_Map<IWString, float> & id_activity_hash,
                              int experimental_column)
{
  if (buffer.nwords() < 2)
  {
    cerr << "Not enough tokens '" << buffer << "'\n";
    return 0;
  }

  IWString id, act;

  if (experimental_column < 0)  // not set
  {
    if (! buffer.split (id, ' ', act))
    {
      cerr << "Cannot split buffer\n";
      return 0;
    }
  }
  else
  {
    buffer.word(0, id);

    if (! buffer.word(experimental_column, act))
    {
      cerr << "Cannot extract activity data '" << buffer << "'\n";
      return 0;
    }
  }

  if (strip_leading_zeros_from_identifiers)
    id.remove_leading_chars('0');

  act.remove_leading_chars(' ');

  float a;
  if (act.numeric_value(a))
    id_activity_hash[id] = a;
  else if (0 == id_activity_hash.size())   // non numeric header, OK
    ;
  else
  {
    cerr << "Illegal numeric '" << act << "'\n";
    return 0;
  }

  return 1;
}*/



static int
spearman_rank (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vE:e:p:zs:jy");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  int experimental_column = 1;

  if (cl.option_present ('e'))
  {
    if (! cl.value ('e', experimental_column) || experimental_column < 1)
    {
      cerr << "The experimental/measured value column (-e) must be a whole positive number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Experimental/Measured values in column " << experimental_column << endl;

    experimental_column--;
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', header_records_to_skip) || header_records_to_skip < 1)
    {
      cerr << "The number of records to skip (-s option) must be a whole positive number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will skip the first " << header_records_to_skip << " records\n";
  }

  if (cl.option_present('j'))
  {
    is_descriptor_file = 1;

    header_records_to_skip = 1;

    if (verbose)
      cerr << "Will treat as a descriptor file\n";
  }

  if (cl.option_present('y'))
  {
    ignore_missing_identifiers = 1;

    if (verbose)
      cerr << "Will ignore missing identifiers\n";
  }

  IW_STL_Hash_Map<IWString, float> id_activity_hash;

  IW_STL_Hash_Set ids_present_in_expt;

  Set_of_Numbers expt;
  expt.set_column (experimental_column);

  if (cl.option_present ('E'))
  {
    const char * e = cl.option_value ('E');

    if (cl.option_present ('z'))
    {
      strip_leading_zeros_from_identifiers = 1;

      if (verbose)
        cerr << "Will strip leading 0's from identifiers\n";
    }

    if (! expt.build_from_file (e))
    {
      cerr << "Cannot read experimental data from '" << e << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Read " << expt.n() << " activity values from '" << e << "'\n";

    if (! cl.option_present ('p'))
    {
      cerr << "Must specify one or more predicted columns via the -p option\n";
      usage (5);
    }

    expt.ids_present (ids_present_in_expt);

//  for (IW_STL_Hash_Map<IWString, float>::const_iterator i = id_activity_hash.begin(); i != id_activity_hash.end(); ++i)
//  {
//    cerr << "Activity for '" << (*i).first << "' is " << (*i).second << endl;
//  }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int npred = 1;

  if (cl.option_present('p'))
    npred = cl.option_count('p');

  Set_of_Numbers * zdata = new Set_of_Numbers[npred];

  if (cl.option_present('p'))
  {
    for (int i = 0; i < npred; i++)
    {
      int c;
      if (! cl.value('p', c, i) || c < 1)
      {
        cerr << "Prediced columns must be +ve whole numbers '" << c << "' invalid\n";
        return 3;
      }

      zdata[i].set_column(c - 1);
    }
  }
  else if (cl.option_present('E'))
    zdata[0].set_column(1);
  else
    zdata[1].set_column(2);

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, don't know hot to process multiple files\n";
    return 3;
  }

  if (! read_predicted_data_values (cl[0], ids_present_in_expt, zdata, npred, expt))
  {
    cerr << "Cannot read predicted values from '" << cl[0] << "'\n";
    return 3;
  }

  expt.assign_ranks();
//expt.echo(cerr);

  if (0 == expt.n())
  {
    cerr << "No experimental data\n";
    return 3;
  }

  for (int i = 0; i < npred; i++)
  {
    if (zdata[i].n() != expt.n())
    {
      cerr << "Data count mismatch, expt " << expt.n() << " vs expt " << zdata[i].n() << endl;
//    return  i + 1;
    }

    zdata[i].assign_ranks();
//  zdata[i].echo(cerr);
    float s = zdata[i].spearman_rank(expt);

    std::cout << "Rank with column " << (zdata[i].column() + 1) << ' ';
    if (zdata[i].column_header().length())
      std::cout << zdata[i].column_header() << ' ';
    std::cout << s << endl;
  }


  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = spearman_rank(argc, argv);

  return rc;
}
