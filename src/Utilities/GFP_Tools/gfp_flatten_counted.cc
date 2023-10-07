/*
  Flattens counted fingerprints to non counted forms
*/

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "sparsefp.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static IWString * tag = NULL;
static int ntags = 0;

//static resizable_array_p<IWString> detect_tags;

static int tdts_read = 0;

static int fingerprints_processed = 0;

static Accumulator_Int<int> non_zero_bits;
static Accumulator_Int<int> acc_nset;
static Accumulator_Int<int> acc_bits_chopped;

static int max_count = 1;

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
// clang-format off
  cerr << "Removes all count information from a gfp file with counted fingerprints\n";
  cerr << " -F <tag>       tag(s) to process, omit to process all\n";
  cerr << " -m <max>       truncate counts at <max> (default 1)\n";
  cerr << " -v             verbose output\n";
// clang-format on

  exit(rc);
}

static int
gfp_flatten_counted_record(const const_IWSubstring & buffer,
                           IWString_and_File_Descriptor & output)
{
  bool process = false;

  if (0 == ntags && buffer.starts_with("NC"))
    process = true;
  else
  {
    for (auto i = 0; i < ntags; ++i)
    {
      if (! buffer.starts_with(tag[i]))
        continue;

      process = true;
      break;
    }
  }

//cerr << process << " from " << buffer << endl;

  if (! process)
  {
    output << buffer << '\n';
    if ('|' == buffer)
      tdts_read++;

    output.write_if_buffer_holds_more_than(4096);

    return 1;
  }

  Sparse_Fingerprint sfp;

  if (! sfp.construct_from_tdt_record(buffer))
  {
    cerr << "Cannot parse input record as fingerprint\n";
    cerr << buffer << '\n';
    return 0;
  }

  acc_bits_chopped.extra(sfp.truncate_counts_at(max_count));

  acc_nset.extra(sfp.nbits());

  fingerprints_processed++;

  for (auto i = 0; i < buffer.length(); ++i)
  {
    output << buffer[i];
    if ('<' == buffer[i])
      break;
  }

  sfp.append_daylight_ascii_form_with_counts_encoded(output);

  output << ">\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
gfp_flatten_counted(iwstring_data_source & input,
                    IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! gfp_flatten_counted_record (buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }
  }

  return 1;
}

static int
gfp_flatten_counted(const char * fname,
     IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_flatten_counted(input, output);
}


static int
gfp_flatten_counted(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vF:m:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_count('F'))
  {
    ntags = cl.option_count('F');
    tag = new IWString[ntags];

    for (auto i = 0; cl.value('F', tag[i], i); i++)
    {
      if (! tag[i].ends_with('<'))
        tag[i] << '<';
    }
  }
  else
  {
    ntags = 1;
    tag = new IWString[ntags];
    tag[0] = "NC";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', max_count) || max_count < 1)
    {
      cerr << "The max feature count option (-m) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will truncate counts at " << max_count << endl;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! gfp_flatten_counted(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << tdts_read << " TDT's, processed " << fingerprints_processed << " fingerprints\n";
    cerr << acc_nset.n() << " fingerprints had between " << acc_nset.minval() << " and " << acc_nset.maxval() << " bits set";
    if (acc_nset.n() > 0)
      cerr << ", ave " << static_cast<float>(acc_nset.average());
    cerr << endl;
    cerr << "Truncated between " << acc_bits_chopped.minval() << " and " << acc_bits_chopped.maxval() << " bits";
    if (acc_bits_chopped.n() > 0)
      cerr << ", ave " << static_cast<float>(acc_bits_chopped.average());
    cerr << endl;
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_flatten_counted(argc, argv);

  return rc;
}
