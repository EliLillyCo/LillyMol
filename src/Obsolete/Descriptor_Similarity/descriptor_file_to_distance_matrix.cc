/*
  Makes a distance matrix from a descriptor file
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#define IWARAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

#define IWDESCRIPTOR_IMPLEMENTATION 1
#define SET_OF_DESCRIPTORS_IMPLEMENTATION 1

#include "iwdescriptor.h"
#include "set_of_descriptors.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Creates a distance matrix from a descriptor file\n";
  cerr << " -S <fname>       name of the distance matrix file to create\n";
  cerr << " -H               write three column format to stdout\n";
  display_distance_metrics (cerr, 'X');
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
descriptors_to_distance_matrix (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vX:S:H");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int distance_computation = DISTANCE_CARTESIAN;

  if (cl.option_present ('X'))
  {
    const_IWSubstring x = cl.string_value ('X');

    distance_computation = determine_distance_type (x, verbose);

    if (0 == distance_computation)
    {
      cerr << "Invalid distance metric '" << x << "'\n";
      usage (4);
    }
  }

  int write_to_stdout = 0;
  int create_distance_matrix = 0;

  if (cl.option_present ('H'))
  {
    write_to_stdout = 1;
  }

  if (cl.option_present ('S'))
  {
    create_distance_matrix = 1;
  }

  if (! write_to_stdout && ! create_distance_matrix)
  {
    cerr << "Must specify the name of the distance matrix via the -S or -H option\n";
    usage (4);
  }

  if (cl.number_elements () > 1)
    cerr << "Extra arguments ignored\n";

  typedef IWDescriptors<float, float> Descriptor;
  Set_of_Descriptors<Descriptor> descriptors;

  IWString header;
  if (! descriptors.build (cl[0], header))
  {
    cerr << "Cannot read descriptors from '" << cl[0] << "'\n";
    return 3;
  }

  int tokens_in_header = header.nwords ();

  if (verbose)
    cerr << "Header contains " << tokens_in_header << " tokens\n";

  int number_descriptors = descriptors.number_elements ();

  if (verbose)
    cerr << "Read " << number_descriptors << " records from '" << cl[0] << "'\n";

  IWDistanceMatrixBase<float> dm;

  if (create_distance_matrix)
  {
    if (! dm.resize (number_descriptors))
    {
      cerr << "Memory failure, cannot resize distance matrix for " << number_descriptors << " molecules\n";
      return 3;
    }

    for (int i = 0; i < number_descriptors; i++)
    {
      const Descriptor & di = descriptors[i];

      dm.set_id (i, di.id ());
    }
  }

// Our cross reference vector is simple

  int * xref = new int[tokens_in_header];

  assert (NULL != xref);

  std::unique_ptr<int[]> free_xref(xref);

  for (int i = 0; i < tokens_in_header; i++)
  {
    xref[i] = i;
  }

  Accumulator<double> acc;

  for (int i = 0; i < number_descriptors; i++)
  {
    const Descriptor & di = descriptors[i];

    for (int j = i + 1; j < number_descriptors; j++)
    {
      const Descriptor & dj = descriptors[j];

      float d = di.distance (dj, distance_computation, tokens_in_header, xref);

      if (write_to_stdout)
        std::cout << di.id () << ' ' << dj.id () << ' ' << d << endl;

      if (create_distance_matrix)
        dm.set (i, j, d);

      if (verbose)
        acc.extra (d);
    }
  }

  if (verbose && acc.n ())
  {
    cerr << "Computed " << acc.n () << " distances between " << acc.minval ();
    if (acc.n () > 1)
      cerr << " and " << acc.maxval () << " ave " << acc.average_if_available_minval_if_not ();
    cerr << endl;
  }

  if (cl.option_present ('S'))
  {
    const_IWSubstring s = cl.string_value ('S');

    if (! dm.do_write (s))
    {
      cerr << "Bad news, could not write distance matrix to '" << s << "'\n";
      return 3;
    }
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptors_to_distance_matrix (argc, argv);

  return rc;
}
