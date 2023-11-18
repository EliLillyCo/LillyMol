/*
  We want to convert a distance matrix to a file of distances
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int square_file_has_header = 0;

static int lower_triangular_form = 0;

/*
  Do they want fast I/O or not
*/

static int fast_io = 0;

static float maxval = static_cast<float>(0.0);

static IWString * string_result = NULL;

static int
establish_default_fast_io ()
{
  string_result = new IWString[101];

  string_result[0] = " 0";
  for (int i = 1; i < 100; i++)
  {
    IWString & s = string_result[i];

    s.resize(4);

    s = " .";
    if (i < 10)
      s << '0';

    s << i;
  }

  string_result[100] = " 1.0";

  return 1;
}

static int
establish_fast_io (int fast_io,
                   float maxval)
{
  string_result = new IWString[fast_io + 1];

  double dx = maxval / static_cast<double>(fast_io);

  assert (dx > 0.0);

  string_result[0] = " 0";

  for (int i = 1; i <= fast_io; i++)
  {
    IWString & s = string_result[i];
    s.resize(6);
    s << ' ';
    s.append_number(i * dx, 3);

//  cerr << " " << (i * dx) << " " << s << endl;
  }

  return 1;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Creates a file of distances from a distance matrix\n";
  cerr << " -k <form>      what kind of output\n";
  cerr << " -k 3           3 column format 'id1 id2 dist'\n";
  cerr << " -k sq          square data\n";
  cerr << " -k sqh         square data with a header record\n";
  cerr << " -k lt          lower triangular form\n";
  cerr << " -f def         fast I/O - distances between 0 and 1, 2 digits precision\n";
  cerr << " -f n           fast I/O - N buckets between min dist and max dist\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static void
append_distance (IWString & buffer,
                 float f)
{
  if (0 == fast_io)
  {
    buffer << ' ' << f;
  }
  else if (1 == fast_io)
  {
    int i = static_cast<int>(100.0 * f + 0.499);
    buffer << string_result[i];
  }
  else
  {
    int i = static_cast<int>(f * fast_io / maxval + 0.4999);
    buffer << string_result[i];
  }

  return;
}

static int
write_header (const IWDistanceMatrixBase<float> & dm,
              std::ostream & output)
{
  int n = dm.number_molecules();

  output << "Name";

  for (int i = 0; i < n; i++)
  {
    output << ' ' << dm.id(i);
  }

  output << endl;

  return output.good();
}

static int
write_as_square (const IWDistanceMatrixBase<float> & dm,
                 std::ostream & output)
{
  if (square_file_has_header)
  {
    if (! write_header(dm, output))
      return 0;
  }

  int n = dm.number_molecules();

  IWString output_buffer;
  output_buffer.resize(8 * n);   // 8 is arbitrary

  for (int i = 0; i < n; i++)
  {
    output_buffer.resize_keep_storage(0);

    output_buffer << dm.id(i);

    for (int j = 0; j < i; j++)
    {
//    cerr << "Distance between " << i << " and " << j << " is " << dm.zvalue(i, j) << endl;
      append_distance(output_buffer, dm.zvalue(i, j));
    }

    if (! lower_triangular_form)
    {
      output_buffer << " 0";               // distance from I to I

      for (int j = i + 1; j < n; j++)
      {
        append_distance(output_buffer, dm.zvalue(i, j));
      }
    }

    output_buffer << "\n";

    output << output_buffer;
  }

  return output.good();
}

static int
write_as_three_column_data (const IWDistanceMatrixBase<float> & dm,
                            std::ostream & output)
{
  int n = dm.number_molecules();

  IWString output_buffer;
  output_buffer.resize(80);

  for (int i = 0; i < n; i++)
  {
    const IWString & idi = dm.id(i);

    for (int j = i + 1; j < n; j++)
    {
      output_buffer.resize_keep_storage(0);

      output_buffer << idi << ' ' << dm.id(j);
      append_distance(output_buffer, dm.zvalue(i, j));
      output_buffer << '\n';

      output << output_buffer;
    }
  }

  return output.good();
}

static int
distance_matrix_to_distances (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vk:f:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('k'))
  {
    cerr << "Must specify output type via the -k option\n";
    usage(1);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
    cerr << "Extra arguments ignored\n";

  IWDistanceMatrixBase<float> dm;

  if (! dm.do_read(cl[0]))
  {
    cerr << "Cannot read distance matrix file '" << cl[0] << "'\n";
    return 4;
  }

  if (verbose)
    cerr << "Read distance matrix for " << dm.number_molecules() << " molecules\n";

  if (cl.option_present('f'))
  {
    maxval = dm.maxval();

    const_IWSubstring f = cl.string_value('f');

    if ("def" == f)
    {
      fast_io = 1;
      maxval = static_cast<float>(1.0);
      establish_default_fast_io();
    }
    else if (! f.numeric_value(fast_io) || fast_io < 2)
    {
      cerr << "The number of buckets for fast I/O must be a whole positive number\n";
      usage(4);
    }
    else if (! establish_fast_io(fast_io, maxval))
    {
      cerr << "Fatal error establishing fast I/O\n";
      return 7;
    }

    if (verbose)
      cerr << "Optimised for fast I/O\n";
  }

  const_IWSubstring k = cl.string_value('k');

  int rc = 0;

  if ('3' == k)
    rc = write_as_three_column_data(dm, std::cout);
  else if ("sq" == k)
    rc = write_as_square(dm, std::cout);
  else if ("sqh" == k)
  {
    square_file_has_header = 1;
    rc = write_as_square(dm, std::cout);
  }
  else if ("lt" == k)
  {
    lower_triangular_form = 1;
    rc = write_as_square(dm, std::cout);
  }
  else
  {
    cerr << "Unrecognised -k qualifier '" << k << "\n";
    usage(5);
  }

  if (0 == rc)
  {
    cerr << "Error writing distances\n";
    return 7;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_to_distances(argc, argv);

  return rc;
}
