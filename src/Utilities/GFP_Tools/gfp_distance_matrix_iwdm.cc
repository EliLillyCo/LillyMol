/*
  Computes the distance matrix for a set of fingerprints
*/

#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <memory>
#include <limits>
using namespace std;

#include "cmdline.h"
#include "IWDistanceMatrixBase.h"

#include "iwstring_data_source.h"
#include "iw_tdt.h"
#include "iw_tdt_filter.h"
#include "accumulator.h"
#include "misc.h"

#include "gfp.h"
#include "tversky.h"

static int verbose = 0;

static IW_General_Fingerprint * pool = NULL;

static int pool_size = 0;

static int nfingerprints = 0;

static Accumulator<similarity_type_t> stats;

static IW_TDT_Filter filter;

static similarity_type_t zero_distance_value = 0.0;

static int output_distances = 1;

static Tversky tversky;

static int report = 0;

static int significant_digits = 0;

static float significant_digits_factor = static_cast<float>(0.0);

static int first_token_of_name_only = 1;

static int equal_weight_tanimoto = 0;

static similarity_type_t singleton_threshold = std::numeric_limits<similarity_type_t>::max();

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes the distance matrix for a pool of fingerprints\n";
  cerr << " -m               output as similarities rather than distances\n";
  cerr << " -b               store as bytes rather than floats\n";
  cerr << " -s <size>        specify the number of molecules\n";
  cerr << " -r <number>      report progress every <number> fingerprints\n";
  cerr << " -S <fname>       name of distance matrix file to create\n";
  cerr << " -F -P -W -Q      standard gfp options, enter '-F help' for info\n";
  cerr << " -V ...           standard Tversky options, enter '-V help' for info\n";
  cerr << " -O <option>      TDT filter options, enter '-O help' for info\n";
  cerr << " -d <number>      round distances to <number> significant digits\n";
  cerr << " -q               use equal weight tanimoto function on composite fingerprints\n";
  cerr << " -T <dist>        discard singletons, things with no nbrs within <dist>\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
allocate_pool()
{
  assert(pool_size > 0 && NULL == pool);

  pool = new IW_General_Fingerprint[pool_size];

  if (verbose)
    cerr << "Pool sized for " << pool_size << " molecules\n";

  assert (NULL != pool);

  return 1;
}

static similarity_type_t
round_to_significant_digits (int significant_digits,
                             similarity_type_t d)
{
  int tmp = static_cast<int>(d * significant_digits_factor + static_cast<float>(0.5));

  return static_cast<float>(tmp) / significant_digits_factor;
}

template <typename T>
int
do_output (T & dmc,
           IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < nfingerprints; i++)
  {
    const_IWSubstring n(pool[i].id());

    if (first_token_of_name_only && n.contains(' '))
      n.truncate_at_first(' ');

    dmc.set_id(i, n);
  }

  return dmc.do_write(output);
}

template <typename T>
int
distance_matrix2 (T & dmc,
                  IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < nfingerprints; i++)
  {
    IW_General_Fingerprint & fpi = pool[i];

    int jstart;    // need to be careful if we have a possibly asymmetric distance measure
    if (tversky.active() && ! tversky.optimistic_mode())
      jstart = 0;
    else
      jstart = i + 1;

    for (int j = jstart; j < nfingerprints; j++)
    {
      similarity_type_t result;

      if (tversky.active())
        result = fpi.tversky(pool[j], tversky);
      else if (equal_weight_tanimoto)
        result = fpi.equal_weight_tanimoto(pool[j]);
      else
        result = fpi.tanimoto(&(pool[j]));

//    cerr << "Between " << i << " and " << j << " distance " << result << endl;

      if (output_distances)
        result = static_cast<similarity_type_t>(1.0) - result;

      if (verbose)
        stats.extra(result);

      if (significant_digits)
        result = round_to_significant_digits(significant_digits, result);

      dmc.set(i, j, result);
    }

    if (0 != report && i > 0 && 0 == i % report)
      cerr << "processed " << i << " fingerprints\n";
  }

  return do_output(dmc, output);
}

static int
build_pool (IW_TDT & tdt)
{
  IW_General_Fingerprint & fp = pool[nfingerprints];
  int fatal;

  if (! fp.construct_from_tdt(tdt, fatal))
  {
    if (! fatal)
      return 1;

    return 0;
  }

  nfingerprints++;

  return 1;
}

static int
build_pool (iwstring_data_source & input)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (filter.active() && ! filter.matches(tdt))
      continue;

    if (! build_pool(tdt))
      return 0;

    if (nfingerprints >= pool_size)
    {
      if (verbose)
        cerr << "Pool is full, " << nfingerprints << endl;
      return 1;
    }
  }

  return 1;
}

static int
build_pool (const char * fname)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size)
  {
    pool_size = input.count_records_starting_with("PCN<");

    if (0 == pool_size)
    {
      cerr << "zero occurrences of '" << "PCN<" << "' in input\n";
      return 0;
    }

    if (! allocate_pool())
      return 0;
  }

  return build_pool(input);
}

template <typename T>
int
distance_matrix (T & iwdm,
                 IWString_and_File_Descriptor & output)
{
  cerr << "Resizing for " << nfingerprints << " fingerprints\n";
  if (! iwdm.resize(nfingerprints))
  {
    cerr << "Memory failure, cannot allocate distance matrix for " << nfingerprints << " fingerprints\n";
    return 8;
  }

  return distance_matrix2(iwdm, output);
}

//#ifdef __GNUG__
template int distance_matrix (IWDistanceMatrixMasquerading_as_Byte<float> &, IWString_and_File_Descriptor &);
template int distance_matrix(IWDistanceMatrixBase<float> &, IWString_and_File_Descriptor &);
template int distance_matrix2(IWDistanceMatrixMasquerading_as_Byte<float> &, IWString_and_File_Descriptor &);
template int distance_matrix2(IWDistanceMatrixBase<float> &, IWString_and_File_Descriptor &);
template int do_output(IWDistanceMatrixMasquerading_as_Byte<float> &, IWString_and_File_Descriptor &);
template int do_output(IWDistanceMatrixBase<float> &, IWString_and_File_Descriptor &);
//#endif

static int
distance_matrix_byte (IWString_and_File_Descriptor & output)
{
  IWDistanceMatrixMasquerading_as_Byte<float>iwdm;

  iwdm.set_range(static_cast<float>(0.0), static_cast<float>(1.0));

  return distance_matrix(iwdm, output);
}

static int
distance_matrix_float (IWString_and_File_Descriptor & output)
{
  IWDistanceMatrixBase<float>iwdm;

  return distance_matrix(iwdm, output);
}

static int
distance_matrix (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:mV:F:P:W:Q:O:r:S:bd:q");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! initialise_fingerprints(cl, verbose))
  {
    cerr << "Cannot initialise fingerprint specifications\n";
    usage(23);
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', report) || report < 1)
    {
      cerr << "Invalid report value (-r option)\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will report progress every " << report << " steps\n";
  }

  if (cl.option_present('O'))
  {
    const_IWSubstring o = cl.string_value('O');

    if ("help" == o)
    {
      display_tdt_filter_syntax(cerr);
      return 31;
    }

    if (! filter.build_from_string(o))
    {
      cerr << "Cannot construct tdt filter from '" << o << "'\n";
      return 14;
    }

    if (verbose)
      cerr << "Filter built from '" << o << "'\n";
  }

  if (cl.option_present('m'))
  {
    output_distances = 0;
    zero_distance_value = 1.0;

    if (verbose)
      cerr << "Will output similarity values rather than distances\n";
  }

  if (cl.option_present('d'))
  {
    if (! cl.value('d', significant_digits) || significant_digits < 1)
    {
      cerr << "Invalid significant digits directive (-d), must be whole number > 0\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will round distances to " << significant_digits << " significant digits only\n";

    significant_digits_factor = 1.0;
    for (int i = 0; i < significant_digits; i++)
    {
      significant_digits_factor *= 10.0;
    }
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      cerr << "The -s option must be followed by a whole number > 2\n";
      usage(5);
    }

    if (! allocate_pool())
      return 12;
  }

  if (cl.option_present('V'))
  {
    if (! tversky.parse_command_line(cl, 'V', verbose))
    {
      cerr << "Cannot get Tversky specifications\n";
      usage(18);
    }

    cerr << "Tversky parameters specified, these are asymetric, do you know what you are doing?\n";
  }

  if (cl.option_present('q'))
  {
    equal_weight_tanimoto = 1;

    if (verbose)
      cerr << "Will use equal weight Tanimoto measure\n";
  }

  if (! cl.option_present('S'))
  {
    cerr << "Must specify name of distance matrix file to create via the -S option\n";
    usage(5);
  }
  
  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (1 != cl.number_elements())
  {
    cerr << "Takes only one command line argument\n";
    usage(3);
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! build_pool(cl[i]))
      return i + 1;
  }

  if (verbose)
    cerr << "Read " << nfingerprints << " fingerprints\n";

  if (0 == nfingerprints)
  {
    cerr << "No fingerprints read\n";
    return 32;
  }

  IWString s = cl.string_value('S');

  IWString_and_File_Descriptor output;
  if (! output.open(s.null_terminated_chars()))
  {
    cerr << "Cannot open '" << s << "'\n";
    return 3;
  }

  int rc;
  if (cl.option_present('b'))
  {
    if (verbose)
      cerr << "Will store bytes masquerading as floats\n";

    rc = distance_matrix_byte(output);
  }
  else
  {
    if (verbose)
      cerr << "Values stored as float values\n";

    rc = distance_matrix_float(output);
  }

  delete [] pool;

  output.close();

  if (verbose)
  {
    cerr << "Computed " << stats.n() << " values, between " << stats.minval() << " and " << stats.maxval() << endl;
    if (stats.n() > 1)
      cerr << "Average " << stats.average() << " variance " << stats.variance() << endl;
  }

  if (0 == rc)
    return 8;

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = distance_matrix(argc, argv);

  return rc;
}
