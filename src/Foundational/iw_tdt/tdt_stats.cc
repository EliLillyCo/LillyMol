/*
  Scans a TDT file and reports on the various items found there
*/

#include <stdlib.h>
#include <iostream>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "accumulator.h"
#include "iw_stl_hash_map.h"
#include "iw_stl_hash_set.h"
#include "iw_tdt.h"
#include "iwstring_data_source.h"

const char * prog_name = NULL;

static int verbose = 0;

static int tdts_read = 0;

static extending_resizable_array<int> dataitem_count;

static int check_for_numeric_values = 0;

static IW_STL_Hash_Set tags_to_check;
static IW_STL_Hash_Set tags_to_ignore;

static IW_STL_Hash_Map_int duplicate_tags_encountered;

class TDT_Tag
{
  private:
    IWString _name;
    int _times_found;
    int _numeric_values;
    Accumulator_Int<int> _size;

  public:
    TDT_Tag (const const_IWSubstring &);

    int extra (const const_IWSubstring &);

    int report (std::ostream &) const;
};

TDT_Tag::TDT_Tag (const const_IWSubstring & t) : _name(t)
{
  _times_found = 0;

  _numeric_values = 0;

  return;
}

int
TDT_Tag::extra (const const_IWSubstring & zdata)
{
  _times_found++;

  _size.extra(zdata.length());

  if (check_for_numeric_values && zdata.length() < 12)
  {
    double tmp;
    _numeric_values += zdata.numeric_value(tmp);
  }

  return 1;
}

int
TDT_Tag::report (std::ostream & os) const
{
  os << "Tag '" << _name << "' encountered " << _times_found << " times";

  if (_times_found)
  {
    if (_size.minval() == _size.maxval())
      os << " each " << _size.minval() << " bytes";
    else
    {
      os << " between " << _size.minval() << " and " << _size.maxval() << " bytes";
      if (_times_found > 1)
        os << " ave " << static_cast<float>(_size.average());
    }
  }

  if (check_for_numeric_values && _numeric_values > 0)
    os << ", " << _numeric_values << " times valid as a numeric";

  os << endl;

  return os.good();
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Scans a TDT file and determines characteristics of data stored\n";
  cerr << " -O <tag>       only check <tag>\n";
  cerr << " -X <tag>       do not check <tag>\n";
  cerr << " -n             try to interpret data as numeric\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

typedef std::unordered_map<IWString, TDT_Tag *, IWStringHash> TAGS;

static int
tdt_stats (const const_IWSubstring & buffer,
           TAGS & tags,
           IW_STL_Hash_Set & found_this_tdt,
           const int line_number)
{
  const_IWSubstring tag, data;

  if (! buffer.split(tag, '<', data))
  {
    cerr << "Cannot split record into before and after '<'\n";
    return 0;
  }

  if (0 == tag.length())
  {
    cerr << "Zero length tag, invalid\n";
    return 0;
  }

  if (! found_this_tdt.contains(tag))
    found_this_tdt.insert(tag);
  else
  {
    duplicate_tags_encountered[tag]++;
    if (verbose)
      cerr << "Duplicate tag '" << tag << "' near line " << line_number << endl;
  }

  if (data.ends_with('\n'))
    data.chop();

  if (! data.ends_with('>'))
  {
    cerr << "TDT data must end with '>'\n";
    return 0;
  }

  if (0 == data.length())    // how could this happen?
  {
    cerr << "Very messed up data!\n";
    return 0;
  }

  data.chop();         // get rid of the > character

  TAGS::iterator f = tags.find(tag);

  if (tags_to_ignore.size() && tags_to_ignore.contains(tag))
    return 1;
  else if (tags_to_check.size() && ! tags_to_check.contains(tag))
    return 1;

  if (f == tags.end())
  {
    TDT_Tag * t = new TDT_Tag(tag);
    tags[tag] = t;

    f = tags.find(tag);

    assert (f != tags.end());
  }

  TDT_Tag * t = (*f).second;

  t->extra(data);

  return 1;
}

static int
tdt_stats (const IW_TDT & tdt,
           TAGS & tags,
           const int line_number)
{
  int i = 0;
  const_IWSubstring d;

  int dataitem_count_here = 0;

  IW_STL_Hash_Set found_this_tdt;

  while (tdt.next_dataitem(d, i))
  {
    dataitem_count_here++;
    if (! tdt_stats(d, tags, found_this_tdt, line_number))
    {
      cerr << "INvalid TDT '" << d << "'\n";
      return 0;
    }
  }

  dataitem_count[dataitem_count_here]++;

  return 1;
}

static int
tdt_stats (iwstring_data_source & input,
           TAGS & tags)
{
  int l = input.lines_read();

  IW_TDT tdt;
  while (tdt.next(input))
  {
    tdts_read++;

    if (! tdt_stats(tdt, tags, l))
    {
      cerr << "Invalid TDT, at line " << l << endl;
      cerr << tdt;
      return 0;
    }

    l = input.lines_read();
  }

  return tags.size();
}

static int
tdt_stats (const char * fname,
           TAGS & tags)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tdt_stats(input, tags);
}

static int
process_tags (Command_Line & cl,
              char flag,
              IW_STL_Hash_Set & h)
{
  int i = 0;
  IWString s;

  while(cl.value(flag, s, i++))
  {
    if (s.ends_with('<'))
      s.chop();

    h.insert(s);
  }

  return h.size();
}

static int
tdt_stats (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vnO:X:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    check_for_numeric_values = 1;

    if (verbose)
      cerr << "Will try to deciper data as numeric\n";
  }

  if (cl.option_present('O') && cl.option_present('X'))
  {
    cerr << "Cannot use both -O and -X options\n";
    usage(2);
  }

  if (cl.option_present('O'))
  {
    if (! process_tags(cl, 'O', tags_to_check))
    {
      cerr << "Cannot process tags to check (-O)\n";
      return 2;
    }

    if (verbose)
      cerr << "Will ignore any of " << tags_to_ignore.size() << " tags\n";
  }
  else if (cl.option_present('X'))
  {
    if (! process_tags(cl, 'X', tags_to_ignore))
    {
      cerr << "Cannot process tags to ignore (-X)\n";
      return 2;
    }

    if (verbose)
      cerr << "Will ignore any of " << tags_to_ignore.size() << " tags\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  TAGS tags;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! tdt_stats(cl[i], tags))
    {
      rc = i + 1;
      break;
    }
  }

  std::cout << "Read " << tdts_read << " TDT's\n";

  for (int i = 0; i < dataitem_count.number_elements(); i++)
  {
    if (dataitem_count[i])
      cerr << dataitem_count[i] << " TDT's had " << i << " dataitems\n";
  }

  for (TAGS::const_iterator i = tags.begin(); i != tags.end(); i++)
  {
    const TDT_Tag * t = (*i).second;

    t->report(std::cout);

    delete t;
  }

//for (const auto i : duplicate_tags_encountered)
  for (auto i = duplicate_tags_encountered.begin(); i != duplicate_tags_encountered.end(); ++i)
  {
    cerr << i->second << " duplicate " << i->first << " values encountered\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tdt_stats(argc, argv);

  return rc;
}
