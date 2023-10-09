/*
  Tester for iwstring_data_source
*/

#include <stdlib.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <random>
#include <memory>
#include <random>

#include "Foundational/cmdline/cmdline.h"

#include "iwstring_data_source.h"

using std::cout;
using std::cerr;
using std::endl;

static int verbose= 0;
static int abort_on_error = 0;

char * prog_name = nullptr;

static IWString filter_pattern("^#");

int
test1 (const char * fname)
{
  std::cout << "Begin test1\n";

  iwstring_data_source file1 (fname, 4);
  if (! file1.ok ())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  cerr << file1.records_remaining () << " records in file\n";

  file1.set_filter_pattern (filter_pattern);

  cerr << "Filter pattern set\n";

  int rc = 1;
  IWString buffer;
  while (file1.next_record (buffer))
  {
    std::cout << buffer << endl;
    if (verbose)
      cerr << "Line of length " << buffer.length () << " characters\n";
    if (buffer.length () && '#' != buffer[0])
    {
      cerr << "Filter pattern failed\n";
      rc = 0;
    }
  }

  file1.debug_print(std::cout);

  return rc;
}

int
test2 (const char * fname)
{
  cout << "Begin test2\n";

  iwstring_data_source file2 (fname);
  if (! file2.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  IWString buffer;
  while (file2.next_record (buffer))
  {
    cout << "Record " << file2.lines_read () << " '" << buffer << "'\n";
  }

  file2.debug_print (cout);

  return 1;
}

int
test3 (const char * fname)
{
  cout << "Begin test3\n";

  iwstring_data_source file1 (fname, 4);
  if (! file1.ok ())
  {
    cerr << prog_name << ":test3 cannot open '" << fname << "'\n";
    return 0;
  }

  file1.set_ignore_pattern ("^#");
  file1.set_strip_leading_blanks ();

  IWString buffer;
  while (file1.next_record (buffer))
  {
    cout << buffer << endl;
    if (0 == buffer.length ())
      ;
    else if (isspace (buffer[0]))
    {
      cerr << "Strip leading blanks failed '" << buffer << "'\n";
      return 0;
    }
//  else if ('#' == buffer[0])
//  {
//    cerr << "ignore pattern failed\n";
//    return 1;
//  }
  }

  file1.debug_print (cout);

  return 1;
}

static int
test4 (const char * fname)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "test4 with '" << fname << "'\n";

  int nr = input.records_remaining ();

  if (0 == nr)
  {
    cerr << "Empty file '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Test4 file contains " << nr << " records\n";

// Add 1 to the size of the arrays in case unterminated record

  size_t * offset = new size_t[nr + 1];
  IWString * line = new IWString[nr + 1];

  assert (nullptr != offset && nullptr != line);

  size_t o = input.tellg ();
  cerr << "Initial offset " << o << endl;
  int ndx = 0;

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    line[ndx] = buffer;
    offset[ndx] = o;

//#define DEBUG_TEST4
#ifdef DEBUG_TEST4
    cerr << "Record " << ndx << " at " << o << endl;
#endif

    o = input.tellg ();

#ifdef DEBUG_TEST4
    cerr << "Next record starts at " << o << endl;
#endif

    ndx++;

    int r = input.records_remaining ();
    if (r != nr - ndx)
    {
      cerr << "Test4, " << r << " records remaining but on record " << ndx << endl;
      return 0;
    }
  }

#ifdef DEBUG_TEST4
  cerr << " nr " << nr << " ndx " << ndx << endl;
#endif
  if (nr != ndx)
  {
    cerr << "Bad news, read " << ndx << " records, but expected " << nr << endl;
    return 0;
  }

  for (int i = 0; i < nr; i++)
  {
    if (! input.seekg (offset[i]))
    {
      cerr << "Gack, cannot seek to " << offset[i] << endl;
      return 0;
    }

    if (! input.next_record (buffer))
    {
      cerr << "Cannot fetch record " << i << " at " << offset[i] << endl;
      return 0;
    }

    if (buffer != line[i])
    {
      cerr << "Mismatch on line " << i << " got '" << buffer << "', expected '" << line[i] << "'\n";
      return 0;
    }
  }

  int random_tests = nr / 10;
  if (0 == random_tests)
    random_tests = nr;

  std::default_random_engine generator;
  std::uniform_int_distribution<int> whole_file_distribution(1, nr - 1);
  std::uniform_int_distribution<int> zero_to_three(0, 3);


  for (int i = 0; i < random_tests; i++)
  {
    int j = whole_file_distribution(generator);

    size_t o = offset[j];

    if (! input.seekg(o))
    {
      cerr << "Gack, cannot seek to " << o << endl;
      return 0;
    }

    if (! input.next_record(buffer))
    {
      cerr << "Cannot fetch record " << j << " at " << o << endl;
      return 0;
    }

    if (buffer != line[j])
    {
      cerr << "Mismatch on line " << i << " got '" << buffer << "', expected '" << line[j] << "'\n";
      return 0;
    }

//  Do some more reading so we test reading after a seek

    int kstop = zero_to_three(generator);
    for (int k = 0; k < kstop; k++)
    {
      (void) input.next_record (buffer);
    }
  }

  return 1;
}

static int test_echo_bytes(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "test_echo_bytes:cannot open '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring buffer;

  const int total_records = input.records_remaining();

  if (total_records < 4)
  {
    cerr << "test_echo_bytes:file too small\n";
    return 0;
  }

  const size_t file_size = input.file_size();
  char * b = new char[file_size]; std::unique_ptr<char[]> free_b(b);

  cerr << "File " << fname << " contains " << total_records << " records, " << file_size << " bytes\n";

  std::random_device rd;
  std::uniform_int_distribution<int> u(1,total_records - 1);

  for (int i = 0; i < 10; ++i)
  {
    input.seekg(0);

    const int x = u(rd);

    for (int j = 0; j < x; ++j)
    {
      if (! input.next_record(buffer))
      {
        cerr << "Yipes, cannot read record " << j << " from '" << fname << "', read first " << x << " iteration " << i << endl;
        return 0;
      }
    }

    const off_t o = input.tellg();

    const size_t bytes_read = input.copy_raw_bytes(b, file_size - o);

    IWString as_read;
    as_read.resize(bytes_read);

    while (input.next_record(buffer))
    {
      if (as_read.length() > 0)
        as_read << '\n';
      as_read << buffer;
    }

    as_read << '\n';
    if (0 != as_read.strncmp(b, bytes_read))
    {
      cerr << "Mismatch on data, expected " << bytes_read << " bytes, got " << as_read.length() << endl;
      return 0;
    }
  }

  return 1;
}

static int
test_push_record (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  const int nr = input.records_remaining();

  std::random_device rd;
  std::uniform_int_distribution<int> u(0, nr - 1);

  for (int i = 0; i < 10; ++i)
  {
    if (! input.seekg(0))
    {
      cerr << "test_push_record:cannot seek to start of file, iteration " << i << endl;
      return 0;
    }

    const int x = u(rd);

    IWString buffer;
    for (int j = 0; j < x; ++j)
    {
      if (! input.next_record(buffer))
      {
        cerr << "test_push_record:cannot read record " << j << " file contains " << nr << " iteration " << i << endl;
        return 0;
      }
    }

    input.push_record();

    if (input.records_remaining() != nr - x)
    {
      cerr << "test_push_record:read " << x << " of " << nr << " records, but have " << input.records_remaining() << " remaining records\n";
      return 0;
    }

    IWString s;
    if (! input.next_record(s))
    {
      cerr << "test_push_record:cannot retrieve pushed record\n";
      return 0;
    }

    if (s != buffer)
    {
      cerr << "test_push_record:pushed record mismatch\n";
      cerr << "pushed '" << buffer << "'\n";
      cerr << " now   '" << s << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
do_rx_test (const char * fname, const const_IWSubstring & rx)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  int nr = input.records_remaining ();

  if (verbose)
    cerr << "Input file '" << fname << "' has " << nr << " records\n";

  re2::StringPiece rxs(rx.data(), rx.length());
  RE2 regex(rxs);
  if (! regex.ok()) {
    cerr << "Could not compile '" << rx << "'\n";
    return 0;
  }

  int ng = input.grep(regex);

  if (verbose)
    cerr << "Found " << ng << " instances of '" << rx << "'\n";

  const_IWSubstring buffer;

  if (! input.set_filter_pattern (rx))
  {
    cerr << "Could not set filter pattern to '" << rx << "'\n";
    return 0;
  }

  if (verbose)
    input.debug_print (cerr);

  int nfound = 0;
  while (input.next_record (buffer))
  {
    nfound++;
    cerr << "Match : " << buffer << endl;
  }

  if (verbose)
    cerr << nfound << " of " << input.lines_read () << " lines match\n";

  if (nfound != ng)
  {
    cerr << "Yipes, grep found " << ng << " but retrieved " << nfound << " records\n";
    return 0;
  }
  assert (nfound == ng);

  return 1;
}

static void
usage (int rc = 0)
{
  cerr << prog_name << ": usage \n";

  exit (rc);
}

int
tdss (int argc, char **argv)
{
  Command_Line cl (argc, argv, "vadi:F:r:");

  verbose = cl.option_count ('v');

  cerr << "Size of off_t " << sizeof(off_t) << " and streampos " << sizeof (std::streampos) << endl;

  if (cl.option_present ('a'))
    abort_on_error = 1;

  if (cl.option_present('F'))
  {
    cl.value('F', filter_pattern);

    if (verbose)
      cerr << "Filter pattern set to '" << filter_pattern << "'\n";
  }

  if (0 == cl.number_elements ())
  {
    usage (3);
  }

  int failures = 0;

  if (! test1 (cl[0]))
  {
    cerr << "test1 failed\n";
    failures++;
  }

  if (! test2 (cl[0]))
  {
    cerr << "test2 failed\n";
    failures++;
  }

  if (! test3 (cl[0]))
  {
    cerr << "test3 failed\n";
    failures++;
  }

  if (! test4 (cl[0]))
  {
    cerr << "Test4 failed\n";
    failures++;
  }

  if (! test_echo_bytes(cl[0]))
  {
    cerr << "Test echo bytes failed\n";
    failures++;
  }

  if (! test_push_record(cl[0]))
  {
    cerr << "Test push record failed\n";
    failures++;
  }

  if (cl.option_present ('r'))
  {
    const_IWSubstring r;
    cl.value ('r', r);

    do_rx_test (cl[0], r);
  }

  if (failures)
    cout << failures << " tests failed\n";
  else
    cout << "All tests successful\n";

  return failures;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tdss (argc, argv);

  return rc;
}
