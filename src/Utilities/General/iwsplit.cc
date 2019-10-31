/*
  splits files
*/

#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "cmdline_v2.h"
#include "misc.h"
#include "iwstring_data_source.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString stem("iwsplit");

static IWString suffix;

static int items_per_chunk = 0;

static int chunks_written = 0;

static IW_Regular_Expression rx;

static int number_to_create = 0;

static int append_to_existing_series = 0;

static int next_file_name_sequence_number = 0;

static int split_by_size = 0;

static int write_to_stdout = 0;

static int inherit_input_file_suffix = 0;

static int width = 0;

static char fill_char = '0';

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Splits a file into chunks based on regular expressions\n";

  cerr << " -n <number>    items per chunk\n";
  cerr << " -nc <number>   number of chunks to create - will count the file to\n";
  cerr << "                determine the required number of items per chunk\n";
  cerr << " -size          split into equally sized chunks - based on byte count\n";
  cerr << " -rx <rx>       the regular expression to use\n";
  cerr << "                SD  files need '^\\$'\n";
  cerr << "                TDT files need '^\\|'\n";
  cerr << " -s             chunks start with <rx> matched records rather than end (SDF and TDT)\n";
  cerr << " -g <col>       items are defined by identifiers in column <col>\n";
  cerr << " -g sep=x       column separator for the -g column\n";
  cerr << " -tdt           short-hand for TDT files (like .gfp)\n";
  cerr << " -sdf           short-hand for SDF files\n";
  cerr << " -stem <name>   create files with stem <name>\n";
  cerr << " -suffix <xx>   add <xx> to each file created\n";
  cerr << " -ssuffix       create chunks with same suffix as input\n";
  cerr << " -stop <number> create only <number> files\n";
  cerr << " -append        append to an existing set of split files\n";
  cerr << " -j <number>    repeat first <number> records from input in each file created\n";
  cerr << "                useful for splitting descriptor files\n";
  cerr << " -stdout        also write input stream to stdout\n";
  cerr << " -dd <fname>    do NOT do any splitting, just write dd commands to <fname>\n";
  cerr << " -w <n>         width for sequence number field\n";
  cerr << " -w char=x      character to use for filling to the -w width\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

/*
  There are several ways we determine the number of items in the file.
  Might just be a number of records, or a regular expression, or something
  based on grouping.
  All those methods need a common means of checking what they find
*/

static int
common_check_items_per_chunk (const int items_in_file,
                              const int number_chunks_reqested,
                              int & items_per_chunk)
{

  if (0 == items_in_file)
  {
    cerr << "No items in file!\n";
    return 0;
  }

  if (number_chunks_reqested > items_in_file)
  {
    cerr << "Requested " << number_chunks_reqested << " chunks, but only " << items_in_file << " items in input file\n";
    items_per_chunk = 1;

    return 1;
  }

  if (0 == items_in_file % number_chunks_reqested)
    items_per_chunk = items_in_file / number_chunks_reqested;
  else
  {
    items_per_chunk = items_in_file / number_chunks_reqested + 1;
  }

  if (verbose)
    cerr << "Input file contains " << items_in_file << " items, will put " << items_per_chunk << " items per chunk\n";

  return 1;
}

static int
determine_number_of_items_per_chunk (iwstring_data_source & input,
                                     int & items_in_file,
                                     int number_chunks_reqested,
                                     int & items_per_chunk)
{
  if (! rx.active())
    items_in_file = input.records_remaining();
  else
    items_in_file = input.grep(rx);

  return common_check_items_per_chunk(items_in_file, number_chunks_reqested, items_per_chunk);
}

static int
read_header_records (iwstring_data_source & input,
                     int records_to_read,
                     IWString & header)
{
  const_IWSubstring buffer;

  for (int i = 0; i < records_to_read; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "iwsplit:premature eof reading header records\n";
      return 0;
    }

    header << buffer << '\n';
  }

  return 1;
}

static void
determine_where_to_start (const IWString & stem,
                          const IWString & suffix,
                          int & next_file_name_sequence_number)
{
  while (1)
  {
    IWString tmp;
    tmp << stem << (next_file_name_sequence_number + 1) << suffix;

    if (! dash_s(tmp.null_terminated_chars()))
      return;

    next_file_name_sequence_number++;
  }
}

static int
new_item_found_in_column (const const_IWSubstring & buffer,
                          const int group_column,
                          char group_column_separator,
                          IWString & gid)
{
  const_IWSubstring token;
  int i = 0;

  for (auto col = 0; buffer.nextword_single_delimiter(token, i, group_column_separator); ++col)
  {
//  cerr << "In column " << col << " found '" << token << "'\n";

    if (col < group_column)
      continue;

//  cerr << "Compare '" << token << "' with '" << gid << "'\n";
    if (token == gid)    // not new
      return 0;

    gid = token;
    return 1;
  }

  cerr << "Did not find group column " << (group_column+1) << " in '" << buffer << "'\n";
  return 1;
}


static int
determine_number_of_items_per_chunk (iwstring_data_source & input,
                                     int & items_in_file,
                                     int number_chunks_reqested,
                                     int & items_per_chunk,
                                     const int group_column,
                                     const char group_column_separator)
{
  IWString gid;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (new_item_found_in_column(buffer, group_column, group_column_separator, gid))
      items_in_file++;
  }

  if (! input.seekg(0))
  {
    cerr << "Cannot seek back to beginning of file\n";
    return 0;
  }

  return common_check_items_per_chunk(items_in_file, number_chunks_reqested, items_per_chunk);
}

static int
iwsplit (int argc, char ** argv)
{
  Command_Line_v2 cl (argc, argv, "-v-suffix=s-ssuffix-stem=s-tdt-sdf-n=ipos-nc=ipos-rx=s-stop=ipos-append-j=ipos-suffix_len=ipos-gz-size-stdout-w=s-s-g=s-dd=s");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.good())
  {
    cerr << "Invalid command line option specification(s)\n";
    usage(5);
  }

  if (cl.option_present("suffix") && cl.option_present("ssuffix"))
  {
    cerr << "The -suffix and -ssuffix options are mutually incompatible\n";
    usage(3);
  }

  if (cl.option_present("suffix"))
  {
    cl.value("suffix", suffix);

    if (0 == suffix.length())
      ;
    else if ('.' != suffix[0])
      suffix.insert_before(0, '.');

    if (verbose)
      cerr << "Files created with suffix '" << suffix << "'\n";
  }

  if (cl.option_present("ssuffix"))
  {
    inherit_input_file_suffix = 1;

    if (verbose)
      cerr << "Will inherit suffix from (first) input file\n";
  }

  if (cl.option_present("stem"))
  {
    cl.value("stem", stem);

    if (verbose)
      cerr << "Files created with stem '" << stem << "'\n";
  }

  if (cl.option_present("sdf"))
  {
    rx.set_pattern("^\\$\\$\\$\\$");
    if (verbose)
      cerr << "Regular expression set for SDF input\n";
  }

  if (cl.option_present("tdt"))
  {
    rx.set_pattern("^\\|");
    if (verbose)
      cerr << "Regular expression set for TDT input\n";
  }

  int rx_starts_chunks = 0;

  if (cl.option_present("rx"))
  {
    const_IWSubstring s;
    cl.value("rx", s);

    if (! rx.set_pattern(s))
    {
      cerr << "Invalid record termination regular expression '" <<s << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Regular expression set to '" << rx.source() << "'\n";

    if (cl.option_present('s'))
    {
      rx_starts_chunks = 1;

      if (verbose)
        cerr << "Chunks start with records that match the regular expression\n";
    }
  }

  int group_column_separator = ' ';
  int group_column = -1;

  if (cl.option_present('g'))
  {
    if (cl.option_present('n') || cl.option_present("nc"))
    {
      cerr << "Sorry the -g option is not compatible with the -n or -nc option\n";
      usage(1);
    }

    const_IWSubstring g;
    for (auto i = 0; cl.value('g', g, i); ++i)
    {
      if (g.starts_with("sep="))
      {
        g.remove_leading_chars(4);
        if ("tab" == g)
          group_column_separator = '\t';
        else if ("comma" == g)
          group_column_separator = ',';
        else if (1 != g.length())
        {
          cerr << "The group by column separator must be a single character '" << g << "' invalid\n";
          return 1;
        }
        else
          group_column_separator = g[0];
      }
      else if (! g.numeric_value(group_column) || group_column < 1)
      {
        cerr << "The group column specifier must be a valid column number\n";
        usage(2);
      }
    }

    if (group_column < 0)
    {
      cerr << "Must specify the group identifier column via the -g option\n";
      usage(2);
    }

    if (verbose)
      cerr << "Chunks defined by groups of items in column " << group_column << endl;

    group_column--;
  }

  if (rx_starts_chunks && group_column >= 0)
  {
    cerr << "Cannot specify both the -rx and -g options\n";
    return 2;
  }

  if (cl.option_present("nc"))
  {
    if (cl.option_present("stop"))
    {
      cerr << "The -nc and -stop options are mutually exclusive\n";
      usage(11);
    }

    if (cl.option_present("n"))
    {
      cerr << "The -n and -nc options are mutually exclusive\n";
      usage(3);
    }
  }
  else if (cl.option_present("n"))
    ;
  else
  {
    cerr << "Must specify number of items per chunk via the -n option, or the\n";
    cerr << "number of chunks to create via the -nc option\n";
    usage(8);
  }

  int number_chunks_reqested = 0;

  if (cl.option_present("n"))
  {
    cl.value("n", items_per_chunk);
    if (verbose)
      cerr << "Will put " << items_per_chunk << " items in each chunk\n";
  }
  else if (cl.option_present("nc"))
  {
    cl.value("nc", number_chunks_reqested);
    if (verbose)
      cerr << "Will try to create " << number_chunks_reqested << " chunk files\n";
  }

  if (cl.option_present("size"))
  {
    split_by_size = 1;

    if (verbose)
      cerr << "Will split into equally sized chunks (byte count)\n";
  }

  if (cl.option_present("append"))
  {
    append_to_existing_series = 1;
    if (verbose)
      cerr << "Will append to any existing series\n";
  }

  if (cl.option_present("stop"))
  {
    cl.value("stop", number_to_create);
    if (verbose)
      cerr << "Will stop after creating " << number_to_create << " chunks\n";
  }

  if (cl.option_present("stdout"))
  {
    write_to_stdout = 1;

    if (verbose)
      cerr << "Will write to stdout\n";
  }

  int is_descriptor_file = 0;

  IWString header;

  if (cl.option_present("j"))
  {
    cl.value("j", is_descriptor_file);

    if (verbose)
      cerr << "First " << is_descriptor_file << " header records written to each chunk\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (split_by_size && cl.number_elements() > 1)
  {
    cerr << "Sorry, the -size option cannot be used with multiple input files\n";
    return 4;
  }

  if (split_by_size)
  {
    off_t input_size = dash_s(cl[0]);

    if (0 == input_size)
    {
      cerr << "Missing or empty input file '" << cl[0] << "'\n";
      return 3;
    }

    if (items_per_chunk)
    {
    }
    else if (number_chunks_reqested)
    {
      items_per_chunk = input_size / number_chunks_reqested;
      if (0 == items_per_chunk)   // very strange
        items_per_chunk = 1;

      if (verbose)
        cerr << "Will write " << items_per_chunk << " bytes per chunk\n";
    }
  }

  if (append_to_existing_series)
  {
    determine_where_to_start(stem, suffix, next_file_name_sequence_number);
  }

  if (inherit_input_file_suffix)
  {
    const_IWSubstring fname = cl[0];
    fname.remove_up_to_first('/');
    while(fname.remove_up_to_first('.'))
    {
    }

    suffix << '.' << fname;

    if (verbose)
      cerr << "Inherited '" << suffix << "' from first input file\n";
  }

  if (cl.option_present("w"))
  {
    int i = 0;
    const_IWSubstring w;
    while(cl.value("w", w, i++))
    {
      if (w.starts_with("char="))
      {
        w.remove_leading_chars(5);
        if (1 != w.length() || ' ' == w[0] || '\t' == w[0])
        {
          cerr << "The sequence width fill characters must be a single, non blank character, '" << w << "' invalid\n";
          usage(1);
        }

        fill_char = w[0];

        if (verbose)
          cerr << "Width fill char set to '" << fill_char << "'\n";
      }
      else if (! w.numeric_value(width) || width < 1)
      {
        cerr << "The sequence number width specification must be a whole +ve number, '" << w << "' invalid\n";
        usage(2);
      }
      else if (verbose)
        cerr << "Sequence number width specification " << width << "\n";
    }
  }

//int mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

  int items_written_this_chunk = 0;

  int lines_written = 0;
  int items_written = 0;

  IWString_and_File_Descriptor output;

  IWString_and_File_Descriptor mystdout(1);   // may or may not get used

  int items_in_file = 0;

  for (int i = 0; i < cl.number_elements(); i++)
  {
    const char * fname = cl[i];

    iwstring_data_source input(fname);
    if (! input.good())
    {
      cerr << "Cannot open '" << fname << "'\n";
      return 4;
    }

    if (split_by_size > 0)
      ;
    else if (number_chunks_reqested > 0)
    {
      bool success;

      if (group_column >= 0)
        success = determine_number_of_items_per_chunk(input, items_in_file, number_chunks_reqested, items_per_chunk, group_column, group_column_separator);
      else
        success = determine_number_of_items_per_chunk(input, items_in_file, number_chunks_reqested, items_per_chunk);

      if (! success)
      {
        cerr << "Cannot determine number of items per chunk\n";
        return 9;
      }

      items_in_file -= is_descriptor_file;
    }

    if (is_descriptor_file && 0 == header.length())
    {
      if (! read_header_records(input, is_descriptor_file, header))
      {
        cerr << "Cannot read " << is_descriptor_file << " header records from '" << fname << "'\n";
        return 5;
      }
    }

    IWString gid;    // used by group_column

    const_IWSubstring buffer;
    while (input.next_record(buffer))
    {
      if (write_to_stdout)
      {
        mystdout << buffer << '\n';
        mystdout.write_if_buffer_holds_more_than(32768);
      }

      int rx_starts_chunks_and_rx_matches_buffer = (rx_starts_chunks && rx.matches(buffer));

//    cerr << "From '" << buffer << "' match " << rx_starts_chunks_and_rx_matches_buffer << " this chunk " << items_written_this_chunk << endl;

//    if (rx.matches(buffer))
//      cerr << buffer << " : " << items_written_this_chunk << endl;

      bool need_new_file = false;

      if (group_column >= 0)
      {
        if (new_item_found_in_column(buffer, group_column, group_column_separator, gid))
          need_new_file = true;
      }
      else if (! output.is_open())
        need_new_file = true;
      else if (rx_starts_chunks_and_rx_matches_buffer && items_written_this_chunk >= items_per_chunk)
        need_new_file = true;

//#define DEBUG_IWSPLIT
#ifdef DEBUG_IWSPLIT
      if (need_new_file)
        cerr << "Just read " << buffer << "', new file needed, items_written_this_chunk " << items_written_this_chunk << ", per chunk " << items_per_chunk << ", gid '" << gid << "'\n";
#endif

      if (need_new_file)
      {
        if (output.is_open())
          output.close();

        next_file_name_sequence_number++;

        if (number_to_create > 0 && chunks_written >= number_to_create)
          break;

        IWString chunk_file_name;
        chunk_file_name << stem;
        if (gid.length() > 0)
        {
          chunk_file_name = stem;
          chunk_file_name << gid << suffix;
        }
        else if (width > 0)
        {
          int insert_at = chunk_file_name.length();
          int expected_length = chunk_file_name.length() + width;

          chunk_file_name << next_file_name_sequence_number;

          while (chunk_file_name.length() < expected_length)
          {
            chunk_file_name.insert_before(insert_at, fill_char);
          }
          chunk_file_name << suffix;
        }
        else
          chunk_file_name << next_file_name_sequence_number << suffix;

        if (! output.open(chunk_file_name.null_terminated_chars()))
        {
          cerr << "Cannot open '" << chunk_file_name << "'\n";
          return 3;
        }

        if (verbose > 1)
          cerr << "Opened '" << chunk_file_name << "'\n";

        lines_written += items_written_this_chunk;

        if (header.length() > 0)
          output << header;

        if (items_in_file > 0 && lines_written)    // user asked for a specific number of chunks, make sure we deliver that many
        {
          int still_to_write = items_in_file - lines_written;
#ifdef DEBUG_IWSPLIT
          cerr << "Starting " << chunk_file_name << ", started with " << items_in_file << " in file, have written " << lines_written << ", remaining " << (items_in_file - lines_written) << endl;
#endif
          int chunks_remaining = number_chunks_reqested - chunks_written;
          if (1 == chunks_remaining)   // write everything
            items_per_chunk = still_to_write;
          else if (chunks_remaining > 0)
          {
            if (0 == still_to_write % chunks_remaining)
              items_per_chunk = still_to_write / chunks_remaining;
            else
              items_per_chunk = still_to_write / chunks_remaining + 1;

#ifdef DEBUG_IWSPLIT
            cerr << "Remaining " << (items_per_chunk * chunks_remaining) << " compare " << still_to_write << " preliminary items_per_chunk " << items_per_chunk << endl;
#endif

            if (0 == items_per_chunk)
              items_per_chunk = 1;
            else if (items_per_chunk * chunks_remaining < still_to_write)
              items_per_chunk++;
          }

          if (verbose > 2)
            cerr << "still_to_write " << still_to_write << " in " << chunks_remaining << " chunks, items_per_chunk updated to " << items_per_chunk << endl;
        }

        items_written_this_chunk = 0;

        chunks_written++;
      }

      if (! output.is_open())
        cerr << "iwsplit: huh! output stream not open!\n";

      output << buffer << '\n';

      output.write_if_buffer_holds_more_than(32768);

      if (! rx.active())    // every record is an item boundary
        ;
      else if (rx_starts_chunks_and_rx_matches_buffer)  // definitely a boundary
        ;
      else if (! rx.matches(buffer))   // not an item boundary
        continue;

      if (split_by_size > 0)
        items_written_this_chunk += buffer.length() + 1;   // don't forget the newline
      else
        items_written_this_chunk++;

//    cerr << "items_written_this_chunk updated to " << items_written_this_chunk << endl;

      if (rx_starts_chunks)   // we close file before we process the next start record
        ;
      else if (group_column >= 0)
        ;
      else if (items_written_this_chunk >= items_per_chunk)
        output.close();            // force next file to open
    }
  }

  if (write_to_stdout)
    mystdout.flush();

  if (verbose)
    cerr << "Created " << chunks_written << " files\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = iwsplit(argc, argv);

  return rc;
}
