/*
  We have performed a nearest neighbour study and want to identify molecules
  that are active, but which have no close neighbours
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iw_tdt/iw_tdt.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IW_STL_Hash_Map_float activity;

static int items_read = 0;

static int items_written = 0;

static int items_removed_by_ignore_file = 0;

static float active_threshold = 0.0;

static float short_distance_cutoff = 0.0;

static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static Accumulator<float> acc_dist;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Identifies active molecules with no close neighbours\n";
  cerr << " -A <fname>     activity file\n";
  cerr << " -a <float>     threshold for the active class\n";
  cerr << " -I <fname>     ignore identifiers - likely items chosen as support vectors\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
read_activity_record (const const_IWSubstring & buffer,
                      IW_STL_Hash_Map_float & activity)
{
  IWString id;
  int i = 0;

  if (! buffer.nextword(id, i))
  {
    cerr << "Cannot extract first token of activity record\n";
    return 0;
  }

  const_IWSubstring token;

  if (! buffer.nextword (token, i))
  {
    cerr << "Cannot extract activity\n";
    return 0;
  }

  float a;

  if (token.numeric_value(a))
    ;
  else if (0 == activity.size())   // ignore errors in header record
    ;
  else
  {
    cerr << "INvalid numeric '" << token << "'\n";
    cerr << "Activity contains " << activity.size() << endl;
    return 0;
  }

  if (a < active_threshold)
    return 1;

  if (activity.find(id) != activity.end())
  {
    cerr << "Duplicate identifier '" << id << "'\n";
    return 0;
  }

  activity[id] = a;

  return 1;
}

static int
read_activity (iwstring_data_source & input,
               IW_STL_Hash_Map_float & activity)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (0 == buffer.length())
      continue;

    if (buffer.starts_with('#'))
      continue;

    if (read_activity_record (buffer, activity))
      ;
    else if (1 == input.lines_read())
      ;
    else
    {
      cerr << "Fatal error reading activity data '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }
  }

  return activity.size();
}

static int
read_activity (IW_STL_Hash_Map_float & activity,
               const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return read_activity (input, activity);
}

static int
read_items_to_ignore_record (const_IWSubstring & buffer,
                      IW_STL_Hash_Map_float & activity)
{
  if (0 == buffer.nwords())
    return 0;

  IWString id = buffer.word(0);

  IW_STL_Hash_Map_float::const_iterator f = activity.find(id);

  if (f == activity.end())
    return 1;

  activity.erase(id);

  items_removed_by_ignore_file++;

  return 1;
}

static int
read_items_to_ignore (iwstring_data_source & input,
                      IW_STL_Hash_Map_float & activity)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    read_items_to_ignore_record (buffer, activity);
  }

  return 1;
}

static int
read_items_to_ignore (const char * fname,
                      IW_STL_Hash_Map_float & activity)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return read_items_to_ignore (input, activity);
}

static int
nn_identify_active_outliers (const IW_TDT & tdt,
                             IWString_and_File_Descriptor & output)
{
  IWString id;

  if (! tdt.dataitem_value(identifier_tag, id))
  {
    cerr << "Cannot extract identifier '" << identifier_tag << "' from tdt\n";
    return 0;
  }

  float d;

  if (! tdt.dataitem_value(distance_tag, d) || d < 0.0)
  {
    cerr << "Cannot extract distance '" << distance_tag << "' from tdt\n";
    return 0;
  }

//cerr << "Item '" << id << " nndist " << d << ", compare " << short_distance_cutoff << ", size " << activity.size() << endl;

  if (d < short_distance_cutoff)  // too close
    return 1;

  IW_STL_Hash_Map_float::const_iterator f = activity.find(id);

  if (f == activity.end())
    return 1;

  output << id << ' ' << (*f).second << ' ' << d << '\n';

  acc_dist.extra(d);

  if (output.size() > 32768)
    output.write_whole_blocks_shift_unwritten();

  items_written++;

  return 1;
}

static int
nn_identify_active_outliers (iwstring_data_source & input,
                             IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    items_read ++;

    if (! nn_identify_active_outliers(tdt, output))
    {
      cerr << "Fatal error processing tdt\n";
      cerr << tdt;
      return 0;
    }
  }

  return 1;
}

static int
nn_identify_active_outliers (const char * fname,
                             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return nn_identify_active_outliers(input, output);
}


static int
nn_identify_active_outliers (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vt:a:A:I:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('A'))
  {
    cerr << "Must specify activity data via the -A option\n";
    usage(3);
  }

  if (! cl.option_present('a'))
  {
    cerr << "MUst specify active class threshold via the -t option\n";
    usage(2);
  }

  if (cl.option_present('a'))
  {
    if (! cl.value('a', active_threshold))
    {
      cerr << "Invalid active threshold (-a)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Molecules with activity above " << active_threshold << " considered active\n";
  }

  if (! cl.option_present('t'))
  {
    cerr << "Must specify the short-distance cutoff via the -t option\n";
    usage (3);
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', short_distance_cutoff) || short_distance_cutoff <= 0.0 || short_distance_cutoff > 1.0)
    {
      cerr << "The short distance cutoff value (-t) must be a valid distance\n";
      usage(5);
    }

    if (verbose)
      cerr << "MOlecules with nbrs within " << short_distance_cutoff << " are not considered singletons\n";
  }

  if (cl.option_present('A'))
  {
    const char * a = cl.option_value('A');

    if (! read_activity(activity, a))
    {
      cerr << "Cannot read activity file '" << a << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Read " << activity.size() << " active dataitems from '" << a << "'\n";
  }

  if (cl.option_present('I'))
  {
    const char * i = cl.option_value('I');

    if (! read_items_to_ignore(i, activity))
    {
      cerr << "Cannot read items to ignore from '" << i << "'\n";
      return 3;
    }

    if (verbose)
      cerr << items_removed_by_ignore_file << " items removed by -I file '" << i << "'\n";
  }

  if (0 == activity.size())
  {
    cerr << "At threshold " << active_threshold << " no items are active\n";
    return 6;
  }

#ifdef ECHO_ACTIVITIES
  for (IW_STL_Hash_Map_float::const_iterator i = activity.begin(); i != activity.end(); ++i)
  {
    cerr << (*i).first << " activity " << (*i).second << endl;
  }
#endif

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  output << "Name Activity Dist\n";

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! nn_identify_active_outliers(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << items_read << " items, wrote " << items_written << '\n';
    if (acc_dist.n() > 1)
      cerr << "Distances between " << acc_dist.minval() << " and " << acc_dist.maxval() << " ave " << acc_dist.average() << endl;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = nn_identify_active_outliers(argc, argv);

  return rc;
}
