/*
  Tester for mempory mapped file IO
*/

#include <stdlib.h>
#include <random>
#include <memory>

#include "cmdline.h"
#include "iwmmap.h"

const char * prog_name = nullptr;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "tester for memory mapped IWString I/O\n";
  cerr << " -R             use the RAM resident version\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

template <typename T>
int
test_iwmmap (T & input,
             IWString_and_File_Descriptor & output)
{
  int lines_read = 0;
  const_IWSubstring buffer;
  int bytes_read = 0;

  const int lines_in_file = input.records_remaining();

  cerr << "Input contains " << lines_in_file << " records\n";

  IWString * zdata = new IWString[lines_in_file]; std::unique_ptr<IWString[]> free_zdata(zdata);
  off_t * offsets = new off_t[lines_in_file + 1];std::unique_ptr<off_t[]> free_offsets(offsets);

  int next_report = 1000;

  offsets[lines_read] = 0;
  while (input.next_record(buffer))
  {
    output << buffer << '\n';
    zdata[lines_read] = buffer;
    lines_read++;
    bytes_read += buffer.length();
    output.flush();
    if (lines_read + input.records_remaining() != lines_in_file)
      cerr << "Mismatch on lines_read " << lines_read << " records_remaining " << input.records_remaining() << " and lines_in_file " << lines_in_file << endl;
//  cerr << lines_read << " wrote " << buffer.length() << " bytes\n";
    offsets[lines_read] = input.tellg();

    if (lines_read == next_report)
    {
      cerr << "read " << lines_read << " lines\n";
      next_report += 1000;
    }
  }


  if (lines_read != input.lines_read())
  {
    cerr << "test_iwmmap:line count mismatch, got " << lines_read << " object reports " << input.lines_read() << endl;
    return 0;
  }

  cerr << "bytes_read " << bytes_read << " lines_read " << lines_read << " tot " << (bytes_read + lines_read) << " cmp file size " << input.file_size() << endl;

  if (! input.seekg(0))
  {
    cerr << "Cannot seek back to beginning of file\n";
    return 0;
  }

  lines_read = 0;
  while (input.next_record(buffer))
  {
    if (buffer !=  zdata[lines_read])
    {
      cerr << "Data mismatch after seek, line " << lines_read << endl;
    }
    lines_read++;
  }

  std::random_device rd;
  std::uniform_int_distribution<int> u(0, lines_in_file-1);

  for (int i = 0; i < 100; ++i)
  {
    const int l = u(rd);

    if (! input.seekg(offsets[l]))
    {
      cerr << "Huh, cannot seek to " << offsets[l] << endl;
      continue;
    }

    if (! input.next_record(buffer))
    {
      cerr << "HUh, cannot fetch record after seeking to " << offsets[l] << endl;
      continue;
    }

    if (buffer != zdata[l])
    {
      cerr << "Data mismatch after random seek to " << offsets[l] << endl;
      continue;
    }

    input.push_record();

    const_IWSubstring b2;

    if (! input.next_record(b2))
    {
      cerr << "Huh, cannot fetch pushed record\n";
      continue;
    }

    if (b2 != buffer)
    {
      cerr << "Pushed record not the same as just fetched\n";
      cerr << " prev '" << buffer << "'\n";
      cerr << "  now '" << b2 << "'\n";
      continue;
    }
  }

  input.seekg(0);
  input.next_record(buffer);
  input.push_record();
  const_IWSubstring b2;
  if (! input.next_record(b2))
  {
    cerr << "Cannot retrieve first record (pushed)\n";
    return 0;
  }

  if (buffer != b2)
  {
    cerr << "First record buffered mismatch\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
test_iwmmap (const char * fname,
             IWString_and_File_Descriptor & output)
{
  T input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  cerr << "INput file '" << fname << "' has " << input.file_size() << " bytes\n";

  return test_iwmmap(input, output);
}


static int
test_iwmmap (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vR");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (cl.option_present('R'))
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! test_iwmmap<IWString_Data_Source_RAM>(cl[i], output))
      {
        cerr << "RAM file test " << cl[i] << " failed\n";
        rc = i + 1;
        break;
      }
    }
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! test_iwmmap<IWString_Data_Source_MMAP>(cl[i], output))
      {
        cerr << "MMAP file test " << cl[i] << " failed\n";
        rc = i + 1;
        break;
      }
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

  int rc = test_iwmmap(argc, argv);

  return rc;
}
