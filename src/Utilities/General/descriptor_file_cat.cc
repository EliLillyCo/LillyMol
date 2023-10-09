/*
  concatenates descriptor files
*/

#include <iostream>
#include <limits>

#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int check_tcount = 0;

static int columns_in_input = -1;

static int check_same_headers = 1;

static int check_unique_identifiers = 0;

static int header_records_to_skip = 1;

static int perform_excess_zero_trimming = 0;

static int ignore_empty_files = 0;

static int ignore_case_when_comparing_header_records = 0;

static IWString header;

static IW_STL_Hash_Set identifier_seen;

static char column_separator = ' ';

static IWString header_for_fname;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Concatenates descriptor files\n";
  cerr << " -c             ignore case when comparing descriptor names\n";
  cerr << " -t             do NOT check for consistent number of columns\n";
  cerr << " -i <char>      input column separator\n";
  cerr << " -H             suppress checking of same headers\n";
  cerr << " -h             ignore empty (zero size) files\n";
  cerr << " -u             check unique identifiers\n";
  cerr << " -s <number>    number header records to skip (default 1)\n";
  cerr << " -p <hdr>       append extra column with file name\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
identifier_seen_before(const const_IWSubstring & buffer)
{
  const_IWSubstring tmp(buffer);

  tmp.truncate_at_first(' ');

  IW_STL_Hash_Set::const_iterator f = identifier_seen.find(tmp);

  if (f != identifier_seen.end())
    return 1;

  identifier_seen.insert(tmp);

  return 0;
}

static int
descriptor_file_cat(iwstring_data_source & input,
                    const const_IWSubstring & buffer,
                    IWString_and_File_Descriptor & output)
{
  if (check_tcount)
  {
    int nw;
    if (' ' == column_separator)
      nw = buffer.nwords(column_separator);
    else
      nw = buffer.nwords_single_delimiter(column_separator);
    
    if (nw != columns_in_input)
    {
      cerr << "Token count mismatch, got " << buffer.nwords() << " expected " << columns_in_input << endl;
      return 0;
    }
  }

  if (! check_unique_identifiers)
    ;
  else if (identifier_seen_before(buffer))
    return 0;

#ifdef WHERE_IS_COMPRESS_DESCRIPTOR_FILE_RECORD
  if (perform_excess_zero_trimming) {
    IWString tmp;
    compress_descriptor_file_record(buffer, tmp);
    output << tmp;
  } else {
    output << buffer;
  }
#else
  output << buffer;
#endif

  if (header_for_fname.length())
    output << column_separator << input.fname();

  output << '\n';

  return 1;
}

static int
read_and_check_header_records(iwstring_data_source & input,
                               IWString_and_File_Descriptor & output)
{
  IWString buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return ignore_empty_files;
  }

  if (0 == header.length())
  {
    header = buffer;

    if (' ' == column_separator)
      columns_in_input = header.nwords(column_separator);
    else
      columns_in_input = header.nwords_single_delimiter(column_separator);

    if (columns_in_input < 2)
    {
      cerr << "Very strange, descriptor file contains " << columns_in_input << " tokens\n";
    }

    output << buffer;

    if (header_for_fname.length())
      output << column_separator << header_for_fname;
    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }
  else if (! check_same_headers)
    ;
  else if (buffer == header)   // great
    ;
  else 
  {
    int header_records_match = 0;

    if (ignore_case_when_comparing_header_records)
    {
      IWString t1(header);
      IWString t2(buffer);

      t1.to_lowercase();
      t2.to_lowercase();

      if(t1 == t2)
        header_records_match = 1;
      else
        cerr << "Header records do not match\n";
    }

    if (! header_records_match)
    {
      cerr << "Header record mismatch\n";
      cerr << "Prev '" << header << "'\n";
      cerr << "Now  '" << buffer << "'\n";
      return 0;
    }
  }

  for (int i = 1; i < header_records_to_skip; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot read header record " << i << endl;
      return 0;
    }
  }

  return 1;
}

static int
descriptor_file_cat (iwstring_data_source & input,
                     IWString_and_File_Descriptor & output)
{
  if (! read_and_check_header_records(input, output))
    return 0;

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! descriptor_file_cat(input, buffer, output))
    {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
descriptor_file_cat (const char * fname,
                     IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptor_file_cat(input, output);
}


static int
descriptor_file_cat (int argc, char ** argv)
{
  Command_Line_v2 cl(argc, argv, "-v-t-H-h-z-s=int-c-i=s-p=s");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('t'))
  {
    check_tcount = 1;

    if (verbose)
      cerr << "Will check token counts on each line\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', header_records_to_skip) || header_records_to_skip < 0)
    {
      cerr << "The number of header records to skip (-s) option, must be a whole +ve number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will skip " << header_records_to_skip << " header records\n";

    if (0 == header_records_to_skip)    // should use cat
      check_same_headers = 0;
  }

  if (cl.option_present('c'))
  {
    ignore_case_when_comparing_header_records = 1;

    if (verbose)
      cerr << "Will ignore case when comparing descriptor names\n";

    if (header_records_to_skip <= 0)
      header_records_to_skip = 1;
  }

  if (cl.option_present('p'))
  {
    header_for_fname = cl.string_value('p');

    if (verbose)
      cerr << "Will append the file name as an extra column\n";
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i))
    {
      cerr << "Column separator must be a character\n";
      return 2;
    }

    column_separator = i[0];
  }

  if (cl.option_present('z'))
  {
    perform_excess_zero_trimming = 1;

    if (verbose)
      cerr << "Will perform excess zero trimming\n";
  }

  if (cl.option_present('H'))
  {
    check_same_headers = 0;

    if (verbose)
      cerr << "Will NOT check for same headers!\n";
  }

  if (cl.option_present('h'))
  {
    ignore_empty_files = 1;

    if (verbose)
      cerr << "Will ignore empty files\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! descriptor_file_cat(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_cat(argc, argv);

  return rc;
}
