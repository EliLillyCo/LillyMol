/*
  Scans a descriptor file, full of integer descriptors, and
  profiles the "bits" relative to an activity
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/accumulator/accumulator.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int strip_leading_zeros = 0;

static int number_descriptors = 0;
static Accumulator<double> * activity_set = nullptr;
static Accumulator<double> * activity_notset = nullptr;
static IWString * descriptor_name = nullptr;

static int ignore_identifiers_with_no_data = 0;
static int identifiers_with_no_data = 0;

static IWString missing_value('.');

static int skip_first_record_of_activity_file = 0;

static int output_precision = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Scans integer descriptors as fingerprints and compares with an activity\n";
  cerr << " -A <fname>  activity data is in a separate file\n";
  cerr << " -k          skip the first record in the activity file\n";
  cerr << " -z          strip leading zero's from identifiers\n";
  cerr << " -g          ignore identifiers in the descriptor file with no data\n";
//cerr << " -f          ignore floating point values in the descriptor file\n";
  cerr << " -p <prec>   output precision\n";
  cerr << " -v          verbose output\n";

  exit (rc);
}

static void
append_accumulator(IWString & output_buffer,
                   const Accumulator<double> & acc)
{
  output_buffer << acc.n() << ' ';

  if (output_precision > 0)
    output_buffer.append_number(acc.minval(), output_precision);
  else
   output_buffer << static_cast<float>(acc.minval());

  output_buffer << ' ';

  if (output_precision > 0)
    output_buffer.append_number(acc.maxval(), output_precision);
  else
   output_buffer << static_cast<float>(acc.maxval());

  output_buffer << ' ';

  if (output_precision > 0)
    output_buffer.append_number(acc.average_if_available_minval_if_not(), output_precision);
  else
   output_buffer << static_cast<float>(acc.average_if_available_minval_if_not());

  return;
}

static int
report_results (int output_fd)
{
  IWString output_buffer;

  output_buffer << "Name nset set_minval set_maxval set_average notset nset_minval nset_maxval nset_ave diff\n";

  for (int i = 0; i < number_descriptors; i++)
  {
    const Accumulator<double> & acci_set = activity_set[i];
    const Accumulator<double> & acci_nset = activity_notset[i];

    if (0 == acci_set.n())    // bit never set
      continue;

    if (0 == acci_nset.n())  // must be always set
      continue;

    output_buffer << descriptor_name[i] << ' ';
    append_accumulator(output_buffer, acci_set);
    output_buffer << ' ';
    append_accumulator(output_buffer, acci_nset);
    output_buffer << ' ';

    double diff = acci_set.average_if_available_minval_if_not() - acci_nset.average_if_available_minval_if_not();

    if (output_precision > 0)
      output_buffer.append_number(diff, output_precision);
    else
      output_buffer << static_cast<float> (diff);

    output_buffer << '\n';

    if (output_buffer.length() > 32768)
      output_buffer.write_whole_blocks_shift_unwritten(output_fd);
  }

  if (output_buffer.length())
    output_buffer.write(output_fd);

  return 1;
}

static int
parse_descriptor_file_record (const const_IWSubstring & buffer,
                              const IW_STL_Hash_Map_float & id_and_activity)
{
  if (buffer.nwords() != number_descriptors + 1)
  {
    cerr << "Invalid token count, expected " << (number_descriptors + 1) << ", bot " << buffer.nwords() << endl;
    return 0;
  }

  int i = 0;
  IWString id;

  buffer.nextword(id, i);

  IW_STL_Hash_Map_float::const_iterator f = id_and_activity.find(id);

  if (f == id_and_activity.end())
  {
    cerr << "Yipes, no activity data for '" << id << "'\n";
    identifiers_with_no_data++;
    if (ignore_identifiers_with_no_data)
      return 1;
  }

  const_IWSubstring token;
  int col = 0;
  while (buffer.nextword(token, i))
  {
    int b;
    if (! token.numeric_value(b))
    {
      if (token == missing_value)
      {
        col++;
        continue;
      }

      cerr << "Invalid numeric '" << token << "'\n";
      return 0;
    }

    if (b)
      activity_set[col].extra((*f).second);
    else
      activity_notset[col].extra((*f).second);

    col++;
  }

  return 1;
}

static int
descriptor_file_profile_activity (iwstring_data_source & input,
                                  const IW_STL_Hash_Map_float & id_and_activity,
                                  int output_fd)
{
  const_IWSubstring header;

  if (! input.next_record(header))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  number_descriptors = header.nwords() - 1;

  if (number_descriptors <= 0)
  {
    cerr << "Header record seems empty\n";
    return 0;
  }

  activity_set = new Accumulator<double>[number_descriptors];
  activity_notset = new Accumulator<double>[number_descriptors];
  descriptor_name = new IWString[number_descriptors];

  if (NULL == activity_set || NULL == activity_notset || NULL == descriptor_name)
  {
    cerr << "Cannot allocate " << number_descriptors << " accumulators\n";
    return 0;
  }

  int col = 0;
  int i = 0;
  const_IWSubstring token;

  header.nextword(token, i);    // skip first column

  while (header.nextword(token, i))
  {
    descriptor_name[col] = token;
    col++;
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! parse_descriptor_file_record (buffer, id_and_activity))
    {
      cerr << "Invalid record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
descriptor_file_profile_activity (const char * fname,
                                  const IW_STL_Hash_Map_float & id_and_activity,
                                  int output_fd)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptor_file_profile_activity (input, id_and_activity, output_fd);
}

static int
parse_id_and_activity_data(const const_IWSubstring & buffer,
                           IW_STL_Hash_Map_float & id_and_activity)
{
  if (0 == buffer.length() || '#' == buffer[0])
    return 1;

  if (buffer.nwords() < 2)
  {
    cerr << "parse_id_and_activity_data:must be at least two tokens\n";
    return 0;
  }

  IWString id;
  int i = 0;

  buffer.nextword(id, i);

  if (strip_leading_zeros)
    id.remove_leading_chars('0');

  const_IWSubstring token;
  buffer.nextword(token, i);

  float act;

  if (! token.numeric_value(act))
  {
    cerr << "Invalid activity value '" << token << "' id '" << id << "'\n";
    return 0;
  }

  id_and_activity[id] = act;

  return 1;
}

static int 
read_activity_data(iwstring_data_source & input,
                   IW_STL_Hash_Map_float & id_and_activity)
{
  const_IWSubstring buffer;

  if (skip_first_record_of_activity_file)
    input.next_record(buffer);

  while (input.next_record(buffer))
  {
    if (! parse_id_and_activity_data(buffer, id_and_activity))
    {
      cerr << "INvalid id and activity data '" << buffer << "'\n";
      return 0;
    }
  }

  return id_and_activity.size();
}

static int 
read_activity_data(const const_IWSubstring & fname,
                   IW_STL_Hash_Map_float & id_and_activity)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open activity file '" << fname << "'\n";
    return 0;
  }

  input.set_strip_trailing_blanks(1);
  input.set_translate_tabs(1);

  return read_activity_data(input, id_and_activity);
}


static int
descriptor_file_profile_activity (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:kzp:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present('A'))
  {
    cerr << "Must specity activity data via the -A option\n";
    usage(4);
  }

  if (cl.option_present('z'))
  {
    strip_leading_zeros = 1;

    if (verbose)
      cerr << "Will strip leading zero's from identifiers\n";
  }

  if (cl.option_present('k'))
  {
    skip_first_record_of_activity_file = 1;

    if (verbose)
      cerr << "Will skip the header record in the activity file\n";
  }

  IW_STL_Hash_Map_float id_and_activity;

  if (cl.option_present('A'))
  {
    const_IWSubstring a = cl.string_value('A');

    if (! read_activity_data(a, id_and_activity))
    {
      cerr << "Cannot read activity data from '" << a << "'\n";
      return 3;
    }

    if (verbose)
    {
      cerr << "Read " << id_and_activity.size() << " identifier/activity data pairs from '" << a << "'\n";

      Accumulator<double> acc;
      for (IW_STL_Hash_Map_float::const_iterator i = id_and_activity.begin(); i != id_and_activity.end() ; ++i)
      {
        acc.extra((*i).second);
      }

      cerr << "Activity values between " << static_cast<float> (acc.minval()) << " and " << static_cast<float> (acc.maxval()) << " ave " << static_cast<float> (acc.average()) << '\n';
    }
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', output_precision) || output_precision < 1)
    {
      cerr << "INvalid output precision (-p)\n";
      usage(4);
    }

    if (verbose)
      cerr << "output precision set to " << output_precision << endl;
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, cannot handle multiple input files\n";
    return 5;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! descriptor_file_profile_activity (cl[i], id_and_activity, 1))
    {
      rc = i + 1;
      break;
    }

    report_results(1);
  }

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_profile_activity (argc, argv);

  return rc;
}
