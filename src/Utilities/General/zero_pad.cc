/*
  Zero pad a given column
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int column = 0;

static int width = 0;

static char * prepadded = nullptr;

static char pad_char = '0';

static int header_records_to_skip = 0;

static int records_already_wide_enough = 0;

static int records_read = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Zero pad a given column\n";
  cerr << " -c <col>       which column to process, 1 by default\n";
  cerr << " -w <width>     zero pad to <width> width\n";
  cerr << " -p <char>      pad character (default '0')\n";
  cerr << " -s <number>    header records to skip\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
do_insert_padding (const const_IWSubstring & token,
                   IWString_and_File_Descriptor & output)
{
  int l = token.length();

  if (l < width)
    ;
  else if (l == width)
  {
    records_already_wide_enough++;
    output << token;
    return 1;
  }
  else
  {
    cerr << "Token '" << token << "' already wider than " << width << endl;
    output << token;
    return 1;
  }

  output.strncat(prepadded, width - l);

  output << token;

  return 1;
}

static int
zero_pad_record (const const_IWSubstring & buffer,
                 IWString_and_File_Descriptor & output)
{
  const_IWSubstring mybuffer(buffer);

  int col = 0;
  int i = 0;

  const_IWSubstring token;

  while (col < column)
  {
    if (! buffer.nextword(token, i))
    {
      cerr << "Cannot extract column " << (col + 1) << endl;
      return 0;
    }

    if (col > 0)
      output.add(' ');

    output << token;
    col++;
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "Cannot extract pad column " << (col + 1) << endl;
    return 0;
  }

  if (col > 0)
    output << ' ';

  do_insert_padding (token, output);

  while (buffer.nextword(token, i))
  {
    output << ' ';
    output << token;
  }

  output << '\n';

  return 1;
}

static int
zero_pad (iwstring_data_source & input,
          IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot fetch header record " << i << "'\n";
      return 0;
    }

    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);

    records_read++;
  }

  while (input.next_record(buffer))
  {
    records_read++;

    if (! zero_pad_record (buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
zero_pad (const char * fname,
          IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return zero_pad(input, output);
}


static int
zero_pad (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vc:w:p:s:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('c'))
  {
    if (! cl.value('c', column) || column < 1)
    {
      cerr << "The column must be a valid column number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will insert padding at column " << column << endl;

    column--;
  }

  if (cl.option_present('p'))
  {
    const_IWSubstring p = cl.string_value('p');
    if (1 != p.length())
    {
      cerr << "Sorry, the pad character (-p) must be a single character only\n";
      usage(3);
    }

    pad_char = p[0];

    if (verbose)
      cerr << "Pad character set to '" << pad_char << "'\n";
  }

  if (! cl.option_present('w'))
  {
    cerr << "Must specify width of column via the -w option\n";
    usage(3);
  }
  else
  {
    if (! cl.value('w', width) || width < 2)
    {
      cerr << "The width (-w) option must be a valid +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will pad to width " << width << endl;
  }

  prepadded = new char [width];

  memset(prepadded, pad_char, width);

  if (cl.option_present('s'))
  {
    if (! cl.value('s', header_records_to_skip) || header_records_to_skip < 0)
    {
      cerr << "The number of header records to skip (-s) must be a valid +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will skip the first " << header_records_to_skip << " header records from each file\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  output.resize(64000);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! zero_pad(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << records_read << " records read\n";

    if (records_already_wide_enough)
      cerr << records_already_wide_enough << " records already " << width << " or wider\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = zero_pad(argc, argv);

  return rc;
}
