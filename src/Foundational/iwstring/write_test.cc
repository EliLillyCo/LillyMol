/*
  IO test
*/

#include <stdlib.h>
#include <fstream>

#include "cmdline.h"
#include "iwstring_data_source.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int flush_every = 32768;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Tester for I/O, specify one or more files on command line (like cat)\n";
  cerr << " -S <fname>     name of file to create, otherwise uses stdout\n";
  cerr << " -f <bytes>     flush the IW output object every <bytes> types\n";
  cerr << " -g             use iostream\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
write_test(iwstring_data_source & input,
           ostream & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';
  }

  return 1;
}

static int
write_test_stdlib(const char * fname,
            ostream & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return write_test(input, output);
}

static int
write_test(iwstring_data_source & input,
           IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(flush_every);
  }

  return 1;
}

static int
write_test_iwstring(const char * fname,
                    IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return write_test(input, output);
}

static int
write_test_iwstring(Command_Line & cl,
                    IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < cl.number_elements(); i++)
  {
    write_test_iwstring(cl[i], output);
  }

  return 1;
}

static int
write_test_stdlib(Command_Line & cl,
                  ostream & output)
{
  for (int i = 0; i < cl.number_elements(); i++)
  {
    write_test_stdlib(cl[i], output);
  }

  return 1;
}

static int
write_test_iwstring(Command_Line & cl)
{
  if (! cl.option_present('S'))
  {
    IWString_and_File_Descriptor output(1);
    return write_test_iwstring(cl, output);
  }

  IWString_and_File_Descriptor output;

  const char * fname = cl.option_value('S');

  if (! output.open(fname))
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return write_test_iwstring(cl, output);
}

static int
write_test_stdlib(Command_Line & cl)
{
  if (! cl.option_present('S'))
    return write_test_stdlib(cl, cout);

  const char * fname = cl.option_value('S');

  ofstream output;

  output.open(fname, ios::out);

  if (! output.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return write_test_stdlib(cl, output);
}

static int
write_test(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vb:gS:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('b'))
  {
    if (! cl.value('b', flush_every) || flush_every < 1)
    {
      cerr << "The flush every option (-b) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
     cerr << "Will flush buffers when containing " << flush_every << " bytes\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('g'))
    write_test_stdlib(cl);
  else
    write_test_iwstring(cl);

  if (verbose)
  {
  }

  return 0;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = write_test(argc, argv);

  return rc;
}
