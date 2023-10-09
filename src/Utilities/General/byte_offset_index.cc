/*
  Create an index of where records occur in a file.
*/

#include <fstream>
#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int initial_file_size_guess = 2000000;

static int number_chunks = 0;

static int items_per_chunk = 0;

static int write_diffs_as_bytes = 0;

static IWString index_file_name;

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
  cerr << "Makes a byte offset index of a file\n";
  cerr << " -s <records>   initial guess at number of records in file\n";
  cerr << " -n <chunks>    how many chunks\n";
  cerr << " -c <items>     how many items per chunk\n";
  cerr << " -S <fname>     name of the index file to produce\n";
  cerr << " -b             write the inter-position length as bytes rather than records";
  cerr << " -d             write a set of dd commands\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
byte_offset_index(const resizable_array<off_t> & offset,
                   const off_t file_size,
                   std::ostream & output)
{
  int nr = offset.number_elements();

  for (int i = 0; i < nr; i += items_per_chunk)
  {
    off_t o = offset[i];
    output << o << ' ';
    if (write_diffs_as_bytes)
    {
      if ((i + items_per_chunk) >= nr)
      {
        output << (file_size - o);
      }
      else
        output << (offset[i+items_per_chunk] - o);
    }
    else
      output << items_per_chunk;
      
    output << '\n';
  }

  return output.good();
}

static int
byte_offset_index(iwstring_data_source & input,
                   std::ostream & output)
{
  resizable_array<off_t> offset;

  offset.resize(2000000);

  offset.add(static_cast<off_t>(0));

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    const auto o = input.tellg();
    offset.add(o);
  }

  int nr = offset.number_elements();

  if (verbose)
    cerr << "Collected offset data on " << nr << " records\n";

  if (nr < 2)
  {
    cerr << "Sorry, only " << nr << " records in the file\n";
    return 0;
  }

// Compute items_per_chunk if needed

  if (0 == items_per_chunk)
  {
    if (number_chunks > nr)
    {
      cerr << "Sorry, input file contains " << nr << " records, so cannot partition into " << number_chunks << " chunks\n";
      return 0;
    }

    items_per_chunk = nr / number_chunks;

    if (verbose)
      cerr << "With " << nr << " records and " << number_chunks << " chunks, there are " << items_per_chunk << " items in each chunk\n";
  }

  return byte_offset_index(offset, input.tellg(), output);
}

static int
byte_offset_index(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << prog_name << " cannot open '" << fname << "'\n";
    return 0;
  }

  IWString name_of_index_file;
  if (index_file_name.length())
    name_of_index_file = index_file_name;
  else
    name_of_index_file << fname << ".ndx";

  std::ofstream index_file;

  index_file.open(name_of_index_file.null_terminated_chars(), std::ios::out);
  if (! index_file.good())
  {
    cerr << prog_name << " cannot open index file '" << name_of_index_file << "'\n";
    return 0;
  }

  return byte_offset_index(input, index_file);
}

static int
byte_offset_index(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:c:n:S:bd");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('s'))
  {
    if (! cl.value('s', initial_file_size_guess) || initial_file_size_guess < 2)
    {
      cerr << "The initial file size guess must be a positive whole number\n";
      usage(5);
    }

    if (verbose)
      cerr << "initial file record count guess " << initial_file_size_guess << endl;
  }

  if (cl.option_present('c') && cl.option_present('n'))
  {
    cerr << "Must specify just one of the -c or -n options, not both\n";
    usage(4);
  }
  else if (cl.option_present('c'))
  {
    if (! cl.value('c', items_per_chunk) || items_per_chunk < 1)
    {
      cerr << "The items per chunk specification (-c) must be a whole positive number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will have " << items_per_chunk << " items per chunk\n";
  }
  else if (cl.option_present('n'))
  {
    if (! cl.value('n', number_chunks) || number_chunks < 2)
    {
      cerr << "The number of chunks option (-n) must be a whole positive number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will divide into " << number_chunks << " chunks\n";
  }
  else
  {
    cerr << "Must specify one, and only one, of -c or -n options\n";
    usage(4);
  }

  if (cl.option_present('b'))
  {
    write_diffs_as_bytes = 1;

    if (verbose)
      cerr << "Intervals between offsets written as bytes\n";
  }

  if (cl.option_present('S'))
  {
    cl.value('S', index_file_name);
    if (verbose)
      cerr << "Index file name will be '" << index_file_name << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! byte_offset_index(cl[i]))
    {
      rc = i + 1;
      break;
    }
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

  int rc = byte_offset_index (argc, argv);

  return rc;
}
