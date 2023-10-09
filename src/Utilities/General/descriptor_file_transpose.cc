/*
  transposes a descriptor file
*/

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;
const char * prog_name = NULL;

static int verbose = 0;

static char input_token_separator = ' ';
static char output_token_separator = ' ';

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
  cerr << "Transposes a descriptor file\n";
  cerr << " -l <nlines>    number of lines to process\n";
  cerr << " -t             input is tab separated\n";
  cerr << " -c             input is comma separated\n";
  cerr << " -i <sep>       input  token separator\n";
  cerr << " -o <sep>       output token separator\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
descriptor_file_transpose(const resizable_array_p<IWString> & zdata,
                          IWString_and_File_Descriptor & output)
{
  const unsigned int columns_in_input = zdata[0]->nwords(input_token_separator);

  if (verbose)
    cerr << "Input file contains " << columns_in_input << " columns\n";

  if (columns_in_input < 2)
  {
    cerr << "Huh, file contains " << columns_in_input << " columns, cannot process\n";
    return 0;
  }

  if (! output.make_room_for_extra_items(columns_in_input * 10))
  {
    cerr << "Problem sizing output buffer for " << columns_in_input << " columns\n";
    return 6;
  }

  const unsigned int number_records = zdata.size();

// We need an index into each row so we can get the next token

  int * row_ndx = new_int(number_records); std::unique_ptr<int[]> free_row_ndx(row_ndx);

  const_IWSubstring token;

  for (unsigned int i = 0; i < columns_in_input; i++)
  {
    for (unsigned int j = 0; j < number_records; j++)
    {
      if (! zdata[j]->nextword_single_delimiter(token, row_ndx[j], input_token_separator))
      {
        cerr << "Cannot extract column " << (i+1) << " line " << j << endl;
        return 0;
      }

      if (j > 0)
        output << output_token_separator;

      output << token;
    }

    output.add('\n');

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int
descriptor_file_transpose(const char * fname,
                          const int lines_to_process,
                          IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  resizable_array_p<IWString> zdata;
  zdata.resize(10000);

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    zdata.add(new IWString(buffer));

    if (zdata.number_elements() == lines_to_process)
      break;
  }

  const auto number_records = zdata.number_elements();

  if (verbose)
    cerr << "Read " << number_records << " records from '" << fname << "'\n";

  if (number_records < 2)
  {
    cerr << "Only " << number_records << " records in input, cannot continue\n";
    return 0;
  }

  return descriptor_file_transpose(zdata, output);
}

static int
descriptor_file_transpose(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vl:tci:o:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  int lines_to_process = std::numeric_limits<unsigned int>::max();

  if (cl.option_present('l'))
  {
    if (! cl.value('l', lines_to_process) || lines_to_process < 1)
    {
      cerr << "The number of lines to proces (-l) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only process the first " << lines_to_process << " lines in the file\n";
  }

  if (cl.option_present('t'))
  {
    input_token_separator = '\t';
    output_token_separator = '\t';
    if (verbose)
      cerr << "Input assumed to be tab separated\n";
  }

  if (cl.option_present('c'))
  {
    input_token_separator = ',';
    output_token_separator = ',';
    if (verbose)
      cerr << "Input assumed to be comma separated\n";
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i))
    {
      cerr << "Unrecognised input file separator '" << i << "'\n";
      return 1;
    }

    input_token_separator = i[0];
  }

  if (cl.option_present('o'))
  {
    IWString o = cl.string_value('o');
    if (! char_name_to_char(o))
    {
      cerr << "Unrecognised output file separator '" << o << "'\n";
      return 1;
    }

    output_token_separator = o[0];
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, can only handle one file at a time\n";
    return 4;
  }

  IWString_and_File_Descriptor output(1);

  if (! descriptor_file_transpose(cl[0], lines_to_process, output))
  {
    cerr << "Failed to transpose '" << cl[0] << "'\n";
    return 1;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_transpose (argc, argv);

  return rc;
}
