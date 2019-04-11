/*
  Sort a file of TDT's based on the value of a given tag
*/

#include <stdlib.h>
#include <memory>

//#include <pair.h>

#include "cmdline_v2.h"
#include "iwstring_data_source.h"
#include "set_or_unset.h"
#include "iw_tdt.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "iwqsort.h"

static int verbose = 0;

/*
  What do we do with TDT's which don't have a dataitem with the right TAg
*/

static Set_or_Unset<float> missing_value_float;
static Set_or_Unset<IWString> missing_value_string;

static int numeric_sort = 1;

static int sort_by_number_of_records_in_each_tdt = 0;

static extending_resizable_array<int> records_per_tdt;

static int default_global_sort_order = 1;

static int default_global_column_for_sort_key = 0;

static IW_Regular_Expression sort_by_grepc;

/*
  We are often called upon to sort things that look like 'D100' of 'A3'
*/

static int strip_leading_non_numerics_from_strings = 0;

static int
numeric_after_stripping_leading_non_numerics(const_IWSubstring t,   // local copy
                                             float & f)
{
  while (t.length() && ! isdigit(t[0]))
  {
    t++;
  }

  if (0 == t.length())
    return 0;

  return t.numeric_value(f);
}

/*
  For each search criterion, we have a tag and a direction
*/

class TAG_and_Direction : public IWString
{
  private:
    int _which_one_to_retrieve;

    int _direction;

    int _column;

  public:
    TAG_and_Direction();

    int construct_from_command_line_token(const const_IWSubstring & t);

    int direction() const { return _direction;}
    void set_direction(int s) { _direction = s;}

    int which_one_to_retrieve() const { return _which_one_to_retrieve;}

    int build(const IW_TDT &, float &) const;
    int build(const IW_TDT &, IWString &) const;

    int do_comparison(float f1, float f2) const;
    int do_comparison(const IWString & s1, const IWString & s2) const;
};

TAG_and_Direction::TAG_and_Direction()
{
  _which_one_to_retrieve = 0;

  _direction = 1;

  _column = default_global_column_for_sort_key;

  return;
}

/*
*/

int
TAG_and_Direction::construct_from_command_line_token(const const_IWSubstring & t)
{
  const_IWSubstring myt(t);

  if (myt.starts_with('+'))
  {
    _direction = 1;
    myt.remove_leading_chars(1);
  }
  else if (myt.starts_with("-1,"))   // last instance
    ;
  else if (myt.starts_with('-'))
  {
    _direction = -1;
    myt.remove_leading_chars(1);
  }
  else
    _direction = default_global_sort_order;

  if (! myt.contains(','))
  {
    IWString::operator= (myt);

    return 1;
  }

  if (myt.starts_with(',') || myt.ends_with(','))
  {
    cerr << "TAG_and_Direction::construct_from_command_line_token:invalid specification '" << t << "'\n";
    return 0;
  }

  int ndx = myt.index(",col=");
  if (ndx > 0)
  {
    const_IWSubstring tmp(myt);
    tmp.remove_leading_chars(ndx + 5);
    if (! tmp.numeric_value(_column) || _column < 1)
    {
      cerr << "TAG_and_Direction::construct_from_command_line_token:invalid column '" << t << "'\n";
      return 0;
    }

    _column--;
    myt.truncate_at_last(',');

    if (! myt.contains(','))
    {
      IWString::operator=(myt);
      return 1;
    }
  }

// must be numeric qualifier followed by tag

  const_IWSubstring s1, s2;
  if (! myt.split(s1, ',', s2))
  {
    cerr << "TAG_and_Direction::construct_from_command_line_token:invalid directive '" << t << "'\n";
    return 0;
  }

  if (! s1.numeric_value(_which_one_to_retrieve))
  {
    cerr << "TAG_and_Direction::construct_from_command_line_token:invalid numeric '" << t << "'\n";
    return 0;
  }

  if (_which_one_to_retrieve > 0)
    ;
  else if (-1 == _which_one_to_retrieve)
    ;
  else
  {
    cerr << "TAG_and_Direction::construct_from_command_line_token:invalid index " << _which_one_to_retrieve << endl;
    return 0;
  }

  if (_which_one_to_retrieve > 0)
    _which_one_to_retrieve--;

  if (s2.starts_with('+'))
  {
    _direction = 1;
    s2.remove_leading_chars(1);
  }
  else if (s2.starts_with('-'))
  {
    _direction = -1;
    s2.remove_leading_chars(1);
  }

  IWString::operator=(s2);

  return 1;
}

int
TAG_and_Direction::build(const IW_TDT & tdt,
                          float & f) const
{
  const_IWSubstring tag_value;

  if (tdt.dataitem_value(*this, tag_value, _which_one_to_retrieve))
    ;
  else if (missing_value_float.value(f))
    return 1;
  else
  {
    cerr << "TAG_and_Direction::build:cannot extract " << _which_one_to_retrieve << " ocurrence of '" << (*this) << "' from tdt\n";
    return 0;
  }

  if (_column > 0)
    tag_value.remove_leading_words(_column);

  tag_value.truncate_at_first(' ');

  if (tag_value.numeric_value(f))
    return 1;

  if (strip_leading_non_numerics_from_strings && numeric_after_stripping_leading_non_numerics(tag_value, f))
    return 1;

  cerr << "TAG_and_Direction::build:invalid numeric '" << tag_value << "'\n";
  return 0;
}

int
TAG_and_Direction::build(const IW_TDT & tdt,
                          IWString & f) const
{
  if (tdt.dataitem_value(*this, f, _which_one_to_retrieve))
    ;
  else if (missing_value_string.value(f))
    return 1;
  else
  {
    cerr << "TAG_and_Direction::build:cannot extract " << _which_one_to_retrieve << " ocurrence of '" << (*this) << "' from tdt\n";
    return 0;
  }

  if (_which_one_to_retrieve > 0)
    f.remove_leading_words(_which_one_to_retrieve);
  else
    f.truncate_at_first(' ');

  return 1;
}

int
TAG_and_Direction::do_comparison(float f1, float f2) const
{
  if (f1 < f2)
    return - _direction;

  if (f1 > f2)
    return _direction;

  return 0;
}

int
TAG_and_Direction::do_comparison(const IWString & s1,
                                  const IWString & s2) const
{
  return s1.strcmp(s2);
}

static TAG_and_Direction * tag_and_direction = NULL;
static int ntags = 0;

/*
  Since we want to sort these items, we keep everything as pointers
*/

template <typename T>
class Offset_Value 
{
  protected:
    off_t _pos;

    T * _value;

  public:
    Offset_Value();
    ~Offset_Value();

    void set_offset(off_t o) { _pos = o;}

    int  set_single_value(T v);

    int build(IW_TDT &);

    const T * value() const { return _value;}
//  const T & zvalue() const { return _value;}
    off_t offset() const { return _pos;}
};

class Offset_Value_Float : public Offset_Value<float>
{
  private:
  public:
    int build(off_t, IW_TDT &, int &);
};

class Offset_Value_Int : public Offset_Value<int>
{
  private:
  public:
    int build(off_t, IW_TDT &, int &);
};

class Offset_Value_String : public Offset_Value<IWString>
{
  private:
  public:
    int build(off_t, IW_TDT &, int &);
};

int
Offset_Value_Float::build(off_t offset,
                           IW_TDT & tdt,
                           int & missing_values_found)
{
  _pos = offset;

  assert(NULL == _value);

  _value = new float[ntags];
  if (NULL == _value)
  {
    cerr << "Offset_Value_Float::build:cannot allocate " << ntags << " float values\n";
    return 0;
  }

  for (int i = 0; i < ntags; i++)
  {
    if (! tag_and_direction[i].build(tdt, _value[i]))
    {
      cerr << "Offset_Value_Float:build:invalid tdt " << tdt << endl;
      return 0;
    }
  }

  return 1;
}

int
Offset_Value_String::build(off_t offset,
                            IW_TDT & tdt,
                            int & missing_values_found)
{
  _pos = offset;

  _value = new IWString[ntags];
  if (NULL == _value)
  {
    cerr << "Offset_Value_Float::build:cannot allocate " << ntags << " float values\n";
    return 0;
  }

  for (int i = 0; i < ntags; i++)
  {
    if (! tag_and_direction[i].build(tdt, _value[i]))
    {
      cerr << "Offset_Value_Float:build:invalid tdt " << tdt << endl;
      return 0;
    }
  }

  return 1;
}

template class Offset_Value<float>;
template class Offset_Value<IWString>;

template <typename T>
Offset_Value<T>::Offset_Value()
{
  _value = NULL;

  _pos = 0;

  return;
}

template <typename T>
Offset_Value<T>::~Offset_Value()
{
  if (NULL != _value)
    delete [] _value;

  return;
}

template <typename T>
int
Offset_Value<T>::set_single_value(T v)
{
  assert (NULL == _value);

  _value = new T[1];

  _value[0] = v;

  return 1;
}

/*template <typename T>
class TDT_Comparitor_Base
{
  private:
};*/

class TDT_Comparitor_Float
{
  private:
  public:
    int operator () (const Offset_Value_Float &, const Offset_Value_Float &) const;
};

int
TDT_Comparitor_Float::operator ()(const Offset_Value_Float & ov1,
                                   const Offset_Value_Float & ov2) const
{
  const float * v1 = ov1.value();
  const float * v2 = ov2.value();

  for (int i = 0; i < ntags; i++)
  {
//  cerr << "i = " << i << " comparing " << v1[i] << " and " << v2[i] << endl;
    int tmp = tag_and_direction[i].do_comparison(v1[i], v2[i]);

    if (0 != tmp)
      return tmp;
  }

  return 0;
}

class TDT_Comparitor_String
{
  private:
  public:
    int operator()(const Offset_Value_String &, const Offset_Value_String &) const;
};

int
TDT_Comparitor_String::operator()(const Offset_Value_String & ov1,
                                   const Offset_Value_String & ov2) const
{
  const IWString * v1 = ov1.value();
  const IWString * v2 = ov2.value();

  for (int i = 0; i < ntags; i++)
  {
    int tmp = tag_and_direction[i].do_comparison(v1[i], v2[i]);

    if (0 != tmp)
      return tmp;
  }

  return 0;
}

class TDT_Comparitor_NRecs
{
  private:
  public:
    int operator()(const Offset_Value_Float &, const Offset_Value_Float &) const;
};

int
TDT_Comparitor_NRecs::operator()(const Offset_Value_Float & ov1,
                                   const Offset_Value_Float & ov2) const
{
  const float * v1 = ov1.value();
  const float * v2 = ov2.value();

  if (default_global_sort_order > 0)
  {
    if (v1[0] < v2[0])
      return -1;

    if (v1[0] > v2[0])
      return 1;
  }
  else
  {
    if (v1[0] < v2[0])
      return 1;

    if (v1[0] > v2[0])
      return -1;
  }

  return 0;
}

class TDT_Comparitor_int
{
  private:
  public:
    int operator()(const Offset_Value_Int &, const Offset_Value_Int &) const;
};

int
TDT_Comparitor_int::operator()(const Offset_Value_Int & ov1,
                                 const Offset_Value_Int & ov2) const
{
  const int * v1 = ov1.value();
  const int * v2 = ov2.value();

  if (v1[0] < v2[0])
    return 1;

  if (v1[0] > v2[0])
    return -1;

  return 0;
}

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Sorts a TDT file based on values of a tag\n";
  cerr << " -T <tag>       specify the sort field tag\n";
  cerr << " -T +tag        sort on tag in ascending order\n";
  cerr << " -T -tag        sort on tag in ascending order\n";
  cerr << " -T <tag>,col=n sort on column <n> of <tag>\n";
  cerr << " -T <n>,<tag>   which instance of <tag> to use\n";
  cerr << " -grepc <rx>    sort on number of matches to <rx> in each TDT\n";
  cerr << " -c <token>     sort key is column <token> in <tag>\n";
  cerr << " -r             reverse order\n";
  cerr << " -m <value>     value to be assigned to missing values\n";
  cerr << " -s             string comparisons rather than numeric\n";
  cerr << " -y             sort by number of records in each TDT\n";
//cerr << " -w <number>    which instance of <tag> to use (default 1)\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

/*static int
ov_comparitor_float (const void * ov1, const void * ov2)
{
  Offset_Value_Float * p1 = (Offset_Value_Float *) ov1;
  Offset_Value_Float * p2 = (Offset_Value_Float *) ov2;

  if (p1->zvalue () < p2->zvalue ())
    return -reverse_order;
  
  if (p1->zvalue () > p2->zvalue ())
    return reverse_order;

  return 0;
}

static int
ov_comparitor_string (const void * ov1, const void * ov2)
{
  Offset_Value_String * p1 = (Offset_Value_String *) ov1;
  Offset_Value_String * p2 = (Offset_Value_String *) ov2;

  if (p1->zvalue () < p2->zvalue ())
    return -reverse_order;
  
  if (p1->zvalue () > p2->zvalue ())
    return reverse_order;

  return 0;
}*/


template <typename T>
int
echo_the_tdts (iwstring_data_source & input,
               T * ov,
               int ntdts,
               int output_fd)
{
  IWString output_buffer;

  for (int i = 0; i < ntdts; i++)
  {
    const T & v = ov[i];

    off_t offset = v.offset();
    (void) input.seekg(offset);

    IW_TDT tdt;
    tdt.next(input);

    output_buffer << tdt;

    if (output_buffer.length() > 32768)
    {
      output_buffer.write(output_fd);
      output_buffer.resize_keep_storage(0);
    }
  }

  if (output_buffer.length())
    output_buffer.write(output_fd);

  return 1;
}

template int echo_the_tdts(iwstring_data_source &, Offset_Value_Float *, int, int);
template int echo_the_tdts(iwstring_data_source &, Offset_Value_String *, int, int);

template <typename T>
int
read_the_tdts (iwstring_data_source & input,
               T * ov,
               int & ntdts)
{
  ntdts = 0;
  off_t offset = input.tellg();

  int missing_values_found = 0;

  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (! ov[ntdts].build(offset, tdt, missing_values_found))
    {
      cerr << "Cannot process tdt " << ntdts << endl;
      cerr << tdt;
      return 0;
    }

    ntdts++;
    offset = input.tellg();
  }

  if (missing_values_found == ntdts)
  {
    cerr << "All TDT's missing\n";
    return 0;
  }

  if (verbose && missing_values_found)
    cerr << missing_values_found << " of " << ntdts << " TDT's are missing the tag\n";

  return ntdts;
}

template int read_the_tdts(iwstring_data_source &, Offset_Value_Float *, int &);
template int read_the_tdts(iwstring_data_source &, Offset_Value_String *, int &);

static int
do_grepc(const IW_TDT & tdt,
         IW_Regular_Expression & rx)
{
  int n = tdt.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    const_IWSubstring s;
    tdt.item(i, s);

    if (rx.matches(s))
      rc++;
  }

  return rc;
}

static int
read_the_tdts_grepc_version(iwstring_data_source & input,
                             Offset_Value_Int * ov,
                             int & ntdts)
{
  ntdts = 0;
  off_t offset = input.tellg();

  IW_TDT tdt;
  while (tdt.next(input))
  {
    int c = do_grepc(tdt, sort_by_grepc);
    ov[ntdts].set_single_value(c);
    ov[ntdts].set_offset(offset);

    if (verbose)
      records_per_tdt[c]++;

    ntdts++;
    offset = input.tellg();
  }

  if (verbose)
  {
    for (int i = 0; i < records_per_tdt.number_elements(); i++)
    {
      if (records_per_tdt[i])
        cerr << records_per_tdt[i] << " TDT's had " << i << " instances\n";
    }
  }

  return ntdts;
}

static int
read_the_tdts_nrecords_version (iwstring_data_source & input,
                                Offset_Value_Float * ov,
                                int & ntdts)
{
  ntdts = 0;
  off_t offset = input.tellg();

  IW_TDT tdt;
  while (tdt.next(input))
  {
    ov[ntdts].set_single_value(tdt.number_elements());
    ov[ntdts].set_offset(offset);

    if (verbose)
      records_per_tdt[tdt.number_elements()]++;

    ntdts++;
    offset = input.tellg();
  }

  if (verbose)
  {
    for (int i = 0; i < records_per_tdt.number_elements(); i++)
    {
      if (records_per_tdt[i])
        cerr << records_per_tdt[i] << " TDT's had " << i << " records\n";
    }
  }

  return ntdts;
}

static int
tdt_sort (iwstring_data_source & input,
          int output_fd)
{
  int ntdts = input.grep("^\\|$");

  if (verbose)
    cerr << "Input contains " << ntdts << " TDT's\n";

  if (0 == ntdts)
  {
    cerr << "No TDT's in the input\n";
    return 0;
  }

// I really should be able to combine these two blocks into a common template
// function, but the syntax was beyond me

  if (sort_by_number_of_records_in_each_tdt)
  {
    Offset_Value_Float * ov = new Offset_Value_Float[ntdts];
    if (NULL == ov)
    {
      cerr << "Very bad news, cannot allocate " << ntdts << " offset/value pairs\n";
      return 0;
    }

    if (! read_the_tdts_nrecords_version(input, ov, ntdts))
    {
      cerr << "Cannot read the TDT's\n";
      return 0;
    }

    TDT_Comparitor_NRecs tdtcnr;

    iwqsort(ov, ntdts, tdtcnr);

    return echo_the_tdts(input, ov, ntdts, output_fd);
  }

  if (sort_by_grepc.active())
  {
    Offset_Value_Int * ov = new Offset_Value_Int[ntdts];
    if (NULL == ov)
    {
      cerr << "Very bad news, cannot allocate " << ntdts << " offset/value pairs\n";
      return 0;
    }

    if (! read_the_tdts_grepc_version(input, ov, ntdts))
    {
      cerr << "Cannot read the TDT's\n";
      return 0;
    }

    TDT_Comparitor_int tdtcnr;

    iwqsort(ov, ntdts, tdtcnr);

    return echo_the_tdts(input, ov, ntdts, output_fd);
  }

  if (numeric_sort)
  {
    Offset_Value_Float * ov = new Offset_Value_Float[ntdts]; std::unique_ptr<Offset_Value_Float[]> free_ov(ov);
    if (NULL == ov)
    {
      cerr << "Very bad news, cannot allocate " << ntdts << " offset/value pairs\n";
      return 0;
    }

    int ntdts;
    if (! read_the_tdts(input, ov, ntdts) || 0 == ntdts)
    {
      cerr << "Cannot read the TDT's\n";
      return 0;
    }

    TDT_Comparitor_Float tdtcf;

    iwqsort(ov, ntdts, tdtcf);

    return echo_the_tdts(input, ov, ntdts, output_fd);
  }
  else
  {
    Offset_Value_String * ov = new Offset_Value_String[ntdts];
    if (NULL == ov)
    {
      cerr << "Very bad news, cannot allocate " << ntdts << " offset/value pairs\n";
      return 0;
    }

    int ntdts;
    if (! read_the_tdts(input, ov, ntdts) || 0 == ntdts)
    {
      cerr << "Cannot read the TDT's\n";
      return 0;
    }

    TDT_Comparitor_String tdtcs;

    iwqsort(ov, ntdts, tdtcs);

    return echo_the_tdts(input, ov, ntdts, output_fd);
  }

  return 0;    // should never come here
}

static int
tdt_sort(const char * fname,
          int output_fd)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tdt_sort(input, output_fd);
}

static int
tdt_sort(int argc, char ** argv)
{
//Command_Line cl(argc, argv, "vnT:m:rsc:w:zy");
  Command_Line_v2 cl(argc, argv, "-v-n-T=s-m=s-r-s-c=ipos-z-y-grepc=s");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('r'))    // must be done before -T option(s)
  {
    default_global_sort_order = -1;

    if (verbose)
      cerr << "Will sort in reverse order\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', default_global_column_for_sort_key) || default_global_column_for_sort_key < 1)
    {
      cerr << "INvalid identifier column (-c option)\n";
      usage(3);
    }

    if (verbose)
      cerr << "The sort key will be column " << default_global_column_for_sort_key << " in the tag\n";

    default_global_column_for_sort_key--;
  }

  if (cl.option_present('y') && cl.option_present('T'))
  {
    cerr << "Sorry, the -y and -T options are mutually exclusive, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    usage(4);
  }

  if (cl.option_present('y'))
  {
    sort_by_number_of_records_in_each_tdt = 1;
    if (verbose)
      cerr << "Will sort by number of records in each TDT\n";

    ntags = 1;
  }
  else if (cl.option_present("grepc"))
  {
    const_IWSubstring rx = cl.string_value("grepc");

    if (! sort_by_grepc.set_pattern(rx))
    {
      cerr << "Invalid sort by count rx '" << rx << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Will sort by number of matches to '" << sort_by_grepc.source() << "'\n";
  }
  else
  {
    ntags = cl.option_count('T');

    if (0 == ntags)
    {
      cerr << "Must specify sort field via the -T option\n";
      usage(3);
    }

    tag_and_direction = new TAG_and_Direction[ntags];
    if (NULL == tag_and_direction)
    {
      cerr << "Cannot allocate " << ntags << " tags\n";
      return 4;
    }
  }

  if (ntags)
  {
    int i = 0;
    const_IWSubstring t;

    while (cl.value('T', t, i))
    {
      if (! tag_and_direction[i].construct_from_command_line_token(t))
      {
        cerr << "Invalid sort criterion '" << t << "'\n";
        return i + 1;
      }

      i++;
    }
  }

  if (cl.option_present('z'))
  {
    strip_leading_non_numerics_from_strings = 1;
    if (verbose)
      cerr << "Will strip leading non-numeric chars from identifiers\n";
  }

/*if (cl.option_present('w'))
  {
    if (! cl.value ('w', which_one_to_retrieve) || which_one_to_retrieve < 1)
    {
      cerr << "The which tag item to sort (-w option) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will sort on the " << which_one_to_retrieve << " instance of '" << tag << "'\n";

    which_one_to_retrieve--;
  }*/

  if (cl.option_present('s'))
  {
    numeric_sort = 0;
    if (verbose)
      cerr << "Sort based on string comparisons\n";
  }

  if (cl.option_present('m'))
  {
    IWString m = cl.string_value('m');

    if (numeric_sort)
    {
      float tmp;
      if (! m.numeric_value(tmp))
      {
        cerr << "The -m (missing value) option must be followed by a valid number\n";
        usage(4);
      }

      missing_value_float.set(tmp);
      if (verbose)
        cerr << "Missing values will be set to " << tmp << endl;
    }
    else
    {
      missing_value_string.set(m);

      if (verbose)
        cerr << "Missing values will be set to '" << m << "'\n";
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Extra arguments ignored\n";
  }

  if (! tdt_sort(cl[0], 1))
    return 31;

  if (NULL != tag_and_direction)
    delete [] tag_and_direction;

  return 0;
}

int
main(int argc, char ** argv)
{
  int rc = tdt_sort(argc, argv);

  return rc;
}
