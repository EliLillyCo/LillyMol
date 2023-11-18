/*
  Produces vector differences from sparse fingerprints
*/

#include <stdlib.h>
#include <iostream>
#include <limits>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;
using std::numeric_limits;

const char * prog_name = nullptr;

static int verbose = 0;

static float upper_distance_threshold = std::numeric_limits<float>::max();

static float minimum_activity_needed_for_pair = static_cast<float>(0.0);

static IWString identifier_tag("PCN<");

static IWString difference_fingerprint_tag("NCDIFF<");

static Tversky tversky;

static int pairs_written = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Produces vector differences between sparse fingerprints\n";
  cerr << " -A <fname>     activity data - must include header record\n";
  cerr << " -a <activity>  minimum activity needed by either member of pair\n";
  cerr << " -F ...         standard fingerprint options\n";
  cerr << " -V ...         Tversky specifications\n";
  cerr << " -T <dist>      upper distance threshold\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Fingerprint_and_Activity : public IW_General_Fingerprint
{
  private:
    float _activity;

  public:
    Fingerprint_and_Activity();

    float activity() const { return _activity;}
    void set_activity (float s) { _activity = s;}
};

static Fingerprint_and_Activity * pool = nullptr;
static int pool_size = 0;

Fingerprint_and_Activity::Fingerprint_and_Activity()
{
  _activity = static_cast<float>(0.0);

  return;
}

static int
parse_activity_file_record (const const_IWSubstring & buffer,
                            IW_STL_Hash_Map_float & id_to_activity)
{
  int i = 0;
  IWString id;

  if (! buffer.nextword(id, i))
  {
    cerr << "No identifier in activity record\n";
    return 0;
  }

  const_IWSubstring token;

  if (! buffer.nextword (token, i))
  {
    cerr << "Truncated activity file record\n";
    return 0;
  }

  double a;

  if (! token.numeric_value(a))
  {
    cerr << "Invalid activity value '" << token << "'\n";
    return 0;
  }

  id_to_activity[id] = static_cast<float>(a);

  return 1;
}

static int
assign_activity_to_each_vector (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Empty activity file\n";
    return 0;
  }

  IW_STL_Hash_Map_float id_to_activity;

  while (input.next_record (buffer))
  {
    if (! parse_activity_file_record (buffer, id_to_activity))
    {
      cerr << "Invalid activity data '" << buffer << "'\n";
      return 0;
    }
  }

  for (int i = 0; i < pool_size; i++)
  {
    Fingerprint_and_Activity & fpi = pool[i];

    const IWString & id = fpi.id();

    IW_STL_Hash_Map_float::const_iterator f = id_to_activity.find(id);

    if (f == id_to_activity.end())
    {
      cerr << "NO activity data for '" << id << "'\n";
      return 0;
    }

    fpi.set_activity((*f).second);
  }

  return pool_size;
}

class Add_128
{
  private:
  public:
    int operator() (int &) const;
};

int
Add_128::operator() (int & c) const
{
  c += 128;

  return c;
}

static Add_128 add_128;

static int
assign_activity_to_each_vector (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open activity file '" << fname << "'\n";
    return 0;
  }

  return assign_activity_to_each_vector (input);
}

//#define DEBUG_VDIFF

static int
do_output (const Fingerprint_and_Activity & fp1,
           const Fingerprint_and_Activity & fp2,
           float d,
           IWString_and_File_Descriptor & output)
{
  const Sparse_Fingerprint & sfp1 = fp1.sparse_fingerprint(0);
  const Sparse_Fingerprint & sfp2 = fp2.sparse_fingerprint(0);

  Sparse_Fingerprint diff;

  diff.vector_difference(sfp1, sfp2);

  diff.each_count(add_128);

  output << identifier_tag << fp1.id() << ' ' << fp1.activity() << ' ' << fp2.id() << ' ' << fp2.activity() << " dist " << d << " diff " << (fp1.activity() - fp2.activity()) << ">\n";

#ifdef DEBUG_VDIFF
  sfpc.debug_print(cerr);
#endif

  output << difference_fingerprint_tag;

  diff.append_daylight_ascii_form_with_counts_encoded(output);

  output << ">\n";

  output << "|\n";

  return 1;
}

static float
compute_the_distance (Fingerprint_and_Activity & fp1,
                      Fingerprint_and_Activity & fp2)
{
  if (tversky.active())
    return fp1.tversky (fp2, tversky);

  return fp1.distance(fp2);
}

static int
gfp_vector_differences (IWString_and_File_Descriptor & output)
{
  assert (pool_size > 1);

  for (int i = 0; i < pool_size; i++)
  {
    Fingerprint_and_Activity & fpi = pool[i];

    for (int j = 0; j < pool_size; j++)
    {
      if (j == i)
        continue;

      Fingerprint_and_Activity & fpj = pool[j];

      if (fpi.activity() < minimum_activity_needed_for_pair &&
          fpj.activity() < minimum_activity_needed_for_pair)
        continue;

      float d = compute_the_distance(fpi, fpj);

      if (d > upper_distance_threshold)
        continue;

      do_output(fpi, fpj, d, output);

      output.write_if_buffer_holds_more_than(32768);

      pairs_written++;
    }
  }

  return 1;
}

static int
build_pool (iwstring_data_source & input)
{
  int items_in_pool = 0;

  IW_TDT tdt;
  while (tdt.next(input))
  {
    int fatal;
    if (! pool[items_in_pool].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Fatal error processing '" << tdt << "'\n";
        return 0;
      }
    }

    items_in_pool++;

    if (items_in_pool >= pool_size)
      return items_in_pool;
  }

  pool_size = items_in_pool;

  return items_in_pool;
}


static int
build_pool (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size)
  {
    IWString tmp;
    tmp << '^' << identifier_tag;

    std::unique_ptr<re2::RE2> pcn;
    iwre2::RE2Reset(pcn, tmp);
    pool_size = input.grep(*pcn);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    pool = new Fingerprint_and_Activity[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }
  }

  return build_pool(input);
}

static void
truncate_counts(Fingerprint_and_Activity * pool,
                int pool_size)
{
  for (int i = 0; i < pool_size; i++)
  {
    Sparse_Fingerprint & sfp = pool[i].sparse_fingerprint(0);

    sfp.truncate_counts_at(128);
  }

  return;
}

static int
gfp_vector_differences (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:Q:V:T:A:a:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', upper_distance_threshold) || upper_distance_threshold < 0.0 || upper_distance_threshold > 1.0)
    {
      cerr << "The -T option must be followed by a valid similarity value\n";
      usage (12);
    }

    if (verbose)
      cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
  }

  if (cl.option_present ('V'))
  {
    if (! tversky.parse_command_line (cl, 'V', verbose))
    {
      cerr << "Cannot initialise Tversky parameters\n";
      return 8;
    }
  }

  if (! cl.option_present('A'))
  {
    cerr << "Must specify activity file via the -A option\n";
    usage(4);
  }

  if (cl.option_present('a'))
  {
    if (! cl.value('a', minimum_activity_needed_for_pair))
    {
      cerr << "Invalid minimum activity value (-a)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will only consider pairs where at least one has activity " << minimum_activity_needed_for_pair << endl;
  }
  else
    minimum_activity_needed_for_pair = -0.99 *  std::numeric_limits<float>::max();

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage (17);
    }
  }
  else if (! initialise_fingerprints (cl[0], verbose))
  {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Does not handle multiple files\n";
    return 0;
  }

  if (! build_pool (cl[0]))
  {
    cerr << "Cannot read fingerprints '" << cl[0] << "'\n";
    return 3;
  }

  if (verbose)
    cerr << "Read " << pool_size << " fingerprints from '" << cl[0] << "'\n";

  truncate_counts(pool, pool_size);

  if (number_fingerprints() || pool[0].molecular_properties_integer().active())
  {
    cerr << "Sorry, cannot have fixed width fingerprints or properties\n";
    return 5;
  }

  if (1 != number_sparse_fingerprints())
  {
    cerr << "Sorry, only works with 1 sparse fingerprint\n";
    return 6;
  }

  if (cl.option_present('A'))
  {
    const_IWSubstring a = cl.string_value('A');

    if (! assign_activity_to_each_vector (a))
    {
      cerr << "Problems with experimental data '" << a << "'\n";
      return 4;
    }
  }

  IWString_and_File_Descriptor output(1);

  gfp_vector_differences (output);

  output.flush();

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_vector_differences(argc, argv);

  return rc;
}
