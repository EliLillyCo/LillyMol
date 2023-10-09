/*
  Randomly shuffle the contents of a file
*/

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"

using std::cerr;

const char * prog_name = NULL;

static int verbose = 0;

static int header_records_to_retain = 0;

static int lines_per_group = 1;

static std::unique_ptr<re2::RE2> iwcrex_end_of_group;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Must specify an input file\n";
  cerr << "Randomises the records in a file\n";
  cerr << " -h <number>  number of header records to NOT shuffle\n";
  cerr << " -l <number>  number of lines per item (default 1)\n";
  cerr << " -r <regexp>  regular expression for end of grouping\n";
  cerr << " -t <sdf,gfp> handle certain common file forms\n";
//cerr << " -s <seed>    specify random number seed\n";
  cerr << " -v           verbose output\n";

  exit(rc);
}

typedef std::pair<off_t, size_t> Start_Bytes;

static int
gather_offsets_regexp_defines_group_end(iwstring_data_source & input,
                                         std::unique_ptr<re2::RE2> & iwcrex,
                                         resizable_array_p<Start_Bytes> & sb)
{
  const_IWSubstring buffer;

  auto o = input.tellg();
  while (input.next_record(buffer))
  {
    while (! iwre2::RE2PartialMatch(buffer, *iwcrex))
    {
      if (! input.next_record(buffer))
      {
        cerr << "gather_offsets_regexp_defines_group:EOF\n";
        break;
      }
    }
    const auto nxt = input.tellg();

//  cerr << "Found '" << buffer << "' at " << nxt << '\n';

    sb.add(new Start_Bytes(o, nxt - o));

    o = nxt;
  }

  return sb.size();
}

static int
gather_offsets_many_records_per_group(iwstring_data_source & input,
                                       const int lines_per_group,
                                       resizable_array_p<Start_Bytes> & sb)
{
  const_IWSubstring buffer;

  auto o = input.tellg();

  while (input.next_record(buffer))
  {
    for (auto j = 0; j < lines_per_group; ++j)
    {
      if (! input.next_record(buffer))
      {
        cerr << "shuffle_file, premature eof reading " << lines_per_group << " records\n";
        break;
      }
    }

    const auto nxt = input.tellg();

    sb.add(new Start_Bytes(o, nxt - o));

    o = nxt;
  }

  return sb.size();
}

static int
gather_offsets_one_record_per_group(iwstring_data_source & input,
                                     resizable_array_p<Start_Bytes> & sb)
{
  const_IWSubstring buffer;

  auto o = input.tellg();

  while (input.next_record(buffer))
  {
    const auto nxt = input.tellg();

    sb.add(new Start_Bytes(o, nxt - o));

    o = nxt;
  }

  return sb.size();
}

static int
shuffle_file(iwstring_data_source & input,
              IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_retain; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Premature EOF reading header records\n";
      return 0;
    }

    output << buffer << "\n";

    output.write_if_buffer_holds_more_than(8192);
  }

  resizable_array_p<Start_Bytes> sb;

  sb.resize(10000);

  size_t o = input.tellg();

  if (iwcrex_end_of_group)
    gather_offsets_regexp_defines_group_end(input, iwcrex_end_of_group, sb);
  else if (1 == lines_per_group)
    gather_offsets_one_record_per_group(input, sb);
  else if (lines_per_group > 1)
    gather_offsets_many_records_per_group(input, lines_per_group, sb);
  else
  {
    cerr << "shuffle_file:not sure how to group records\n";
    return 0;
  }

  const auto n = sb.size();

  if (verbose)
    cerr << "Read " << n << " records\n";

  std::random_device rd;
  std::mt19937_64 g(rd());

  std::shuffle(sb.rawdata(), sb.rawdata() + n, g);

  for (int i = 0; i < static_cast<int>(n); ++i)
  {
    const Start_Bytes * s = sb[i];

    if (! input.seekg(s->first))
    {
      cerr << "shuffle_file:cannot seek to " << s->first << '\n';
      return 0;
    }

    size_t bytes_echod = 0;
    while (input.next_record(buffer))
    {
      output << buffer << '\n';

      output.write_if_buffer_holds_more_than(8192);

      bytes_echod += buffer.length() + 1;

      if (bytes_echod >= s->second)
        break;
    }
  }

  return 1;
}

static int
shuffle_file(const char * fname,
              IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return shuffle_file(input, output);
}


static int
shuffle_file(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vh:l:r:t:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  int selection_methods = 0;

  if (cl.option_present('h'))
  {
    if (! cl.value('h', header_records_to_retain) || header_records_to_retain < 1)
    {
      cerr << "The number of header records to retain (-h) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will retain " << header_records_to_retain << " records\n";

    selection_methods++;
  }

  if (cl.option_present('l'))
  {
    if (! cl.value('l', lines_per_group) || lines_per_group < 1)
    {
      cerr << "The number of lines per group (-l) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Items consist of " << lines_per_group << " lines\n";

    selection_methods++;
  }

  if (cl.option_present('r'))
  {
    const const_IWSubstring r = cl.option_value('r');

    if (! iwre2::RE2Reset(iwcrex_end_of_group, r))
    {
      cerr << "Cannot initialise regexp pattern '" << r << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Groups are terminated by record matching '" << iwcrex_end_of_group->pattern() << "'\n";

    selection_methods++;
  }

  if (cl.option_present('t'))
  {
    const_IWSubstring t = cl.string_value('t');
    cerr << "T = '" << t << "'\n";
    if ("sdf" == t || "mdl" == t)
    {
      const const_IWSubstring tmp("^\\$\\$\\$\\$");
      iwre2::RE2Reset(iwcrex_end_of_group, tmp);

      selection_methods++;
    }
    else if ("gfp" == t || "tdt" == t)
    {
      const const_IWSubstring tmp("^\\|");
      iwre2::RE2Reset(iwcrex_end_of_group, tmp);

      selection_methods++;
    }
    else
    {
      cerr << "Unrecognised file type '" << t << "'\n";
      usage(2);
    }
  }

  if (selection_methods > 1)
  {
    cerr << "Must specify one, and only one of  -l or -r\n";
    usage(1);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, do not know how to process more than one file at a time\n";
    return 2;
  }

  IWString_and_File_Descriptor output(1);

  if (! shuffle_file(cl[0], output))
  {
    cerr << "Randomisation of '" << cl[0] << "' failed\n";
    return 2;
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

  int rc = shuffle_file(argc, argv);

  return rc;
}
