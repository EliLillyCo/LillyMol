#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "maccskeys_fn5.h"

using std::cerr;

static int number_maccs_keys = 192;

// Always strip to largest fragment, the library uses several
// bonds_between calls.
static int strip_to_largest_fragment = 1;

static int revert_all_directional_bonds_to_non_directional = 0;

static int remove_all_chiral_centres = 0;

/*
  Another possibility for output is in TDT form for fpsim.
*/

static int nbits = 256;    // the default fingerprint width

/*
  When writing a TDT, we need to know the dataitem
*/

static IWString fingerprint_tag;

/*
  Perhaps the number of bits set would be an interesting descriptor.
*/

static int append_nset = 0;

static IWDigits iwdigits;

#ifdef NO_LONGER_USED
/*
  Aug 2005. We may want to delocalise H atoms on aromatic Nitrogen atoms.
*/

static int aromatic_nitrogens_do_not_have_hydrogens = 0;
#endif

/*
  MAY 2000. Introduce the concept of multi-stage fingerprints. For now, we
  use just two levels
*/

static IWString level2_fingerprint_tag;

static IWString non_colliding_fingerprint_tag;

static IWString fixed_size_counted_fingerprint_tag;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

/*
  We can also work as a tester, generating random smiles and testing
  for invariance of the keys
*/

static int ntest = 0;

static int report_test_failures = 1;

static Report_Progress report_progress;

/*
  By default, we stop on a test error, but we can optionally continue
*/

static int keep_going_after_test_failure = 0;

/*
  If we are keeping going after test failures, we should keep track of
  the number of errors
*/

static int molecules_with_test_failures = 0;

static int test_failures = 0;

static Chemical_Standardisation chemical_standardisation;

static int flush_after_each_molecule = 0;

/*
  If we are reporting statistics, we need some accumulators
*/

class int_accumulator : public Accumulator_Int<int>
{
  private:
    int _number_zero;

  public:
    int_accumulator();

    int extra (int);

    int number_zero() const { return _number_zero;}

    double average_number_of_hits() const;
};

int_accumulator::int_accumulator()
{
  _number_zero = 0;
}

int
int_accumulator::extra (int i)
{
  if (0 == i)
    _number_zero++;
  else
    Accumulator_Int<int>::extra (i);

  return 1;
}

/*
  The default average() method returns the average number of times the
  bit was hit - including the influence from non-matches. For this method,
  we exclude the influence of non-matches
*/

double
int_accumulator::average_number_of_hits() const
{
  assert (_number_zero >= 0);

  if (0 == _number_zero)
    return Accumulator_Int<int>::average();

  int number_non_zero = Accumulator_Int<int>::n() - _number_zero;

  if (0 == number_non_zero)
    return 0.0;

  return Accumulator_Int<int>::sum() / static_cast<double>(number_non_zero);
}

static int_accumulator * accumulators = nullptr;

static extending_resizable_array<int> keys_hit;

/*
  If we are collecting statistics, perhaps the user doesn't even want
  output of the individual keys
*/

static int write_key_values = 1;

static int write_header = 1;

static MACCSKeys mk;     // thing that actually does all the work

static int
do_write_header (IWString_and_File_Descriptor & output)
{
  output << "Name";

  for (int i = 0; i < number_maccs_keys; i++)
  {
    output << " mk_mk" << i;
  }

  if (append_nset)
    output << " mknset";

  output << '\n';

  return 1;
}

static void
MaybeFlush(IWString_and_File_Descriptor& output) {
  if (flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(4096);
  }
}

static int
do_output (const Molecule & m,
           const int * keys,
           IWString_and_File_Descriptor & output)
{
//const IWString & mname = m.name();

  append_first_token_of_name(m.name(), output);

  int nset = 0;

  IWString & s = output;

  for (int i = 0; i < number_maccs_keys; i++)
  {
    if (keys[i] > 0)
      nset++;

    iwdigits.append_number(s, keys[i]);
  }

  if (append_nset)
    iwdigits.append_number(s, nset);

  output += '\n';

  return output.good();
}

static int
all_keys_the_same (const int * mk1, const int * mk2,
                   resizable_array<int> & problematic_keys)
{
  int rc = 1;
  for (int i = 0; i < number_maccs_keys; i++)
  {
    if (mk1[i] != mk2[i])
    {
      if (report_test_failures)
        cerr << "Fingerprint mismatch, bit " << i << " values " << mk1[i] << " vs " << mk2[i] << '\n';
      problematic_keys.add_if_not_already_present(i);
      rc = 0;
    }
  }

  return rc;
}

static int
count_aromatic_rings (Molecule & m)
{
  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    if (m.ringi(i)->is_aromatic())
      rc++;
  }

  return rc;
}

static int
do_comparison_between_key_values (Molecule & m1,
                                  const IWString & initial_smiles,
                                  Molecule & m2,
                                  const IWString random_smiles,
                                  const int * k1,
                                  const int * k2,
                                  int permutation,
                                  resizable_array<int> & problematic_keys)
{
  if (all_keys_the_same (k1, k2, problematic_keys))
    return 1;

  test_failures++;

  if (! report_test_failures)
    return 0;

  cerr << "Yipes, key mismatch on permutation " << permutation << '\n';
  cerr << "Molecule '" << m1.name() << "'\n";
  cerr << initial_smiles << '\n';
  cerr << random_smiles << '\n';

// If these two forms have different unique smiles, that may indicate
// problems with aromaticity/sssr determination

  const IWString & m1usmi = m1.unique_smiles();
  const IWString & m2usmi = m2.unique_smiles();
  if (m1usmi == m2usmi)
    cerr << "Unique_smiles match\n";
  else
  {
    cerr << "Unique smiles mismatch\n";
    cerr << m1usmi << '\n';
    cerr << m2usmi << '\n';
    return 1;
  }

  int nr = m1.nrings();
  if (nr != m2.nrings())
  {
    cerr << "Ring count mismatch!, " << nr << " vs " << m2.nrings() << '\n';
    return 1;
  }

  int a1 = count_aromatic_rings(m1);
  int a2 = count_aromatic_rings(m2);

  if (a1 != a2)
  {
    cerr << "Aromatic ring count differs, " << a1 << " vs " << a2 << " rings\n";
    return 1;
  }

  cerr << "Both forms have " << a1 << " aromatic rings\n";

  return 0;
}

static int
test_maccskeys (Molecule & m,
                const int * keys)
{
  int tkeys[number_maccs_keys + 1];

  IWString initial_smiles = m.smiles();

  resizable_array<int> problematic_keys;

  int failures_this_molecule = 0;

  for (int i = 0; i < ntest; i++)
  {
    IWString random_smiles = m.random_smiles();
    Molecule mtmp;
    if (! mtmp.build_from_smiles (random_smiles))
    {
      cerr << "Yipes, cannot build from smiles '" << random_smiles << "'\n";
      return 0;
    }

    assert (mtmp.natoms() == m.natoms());
    assert (mtmp.nrings() == m.nrings());

    if (! mk (mtmp, tkeys))
    {
      cerr << "Yipes, cannot form maccs keys for variant " << i << '\n';
      cerr << "Smiles '" << random_smiles << "'\n";
      return keep_going_after_test_failure;
    }

    if (do_comparison_between_key_values (m, initial_smiles, mtmp, random_smiles, keys, tkeys, i, problematic_keys))
      continue;

    failures_this_molecule++;

    if (! keep_going_after_test_failure)
    {
      molecules_with_test_failures++;
      return 0;
    }
    break;
  }

// Need to check Kekule forms. Do it randomly

  if (m.contains_aromatic_atoms())
  {
    for (int i = 0; i < ntest; i++)
    {
      m.compute_aromaticity_if_needed();
      set_include_aromaticity_in_smiles (1);
      IWString random_smiles = m.random_smiles();
      set_include_aromaticity_in_smiles (0);

      Molecule mtmp;
      if (! mtmp.build_from_smiles (random_smiles))
      {
        cerr << "Cannot parse smiles '" << random_smiles << "', '" << m.name() << "'\n";
        continue;
      }

      if (! mk (mtmp, tkeys))
      {
        cerr << "Yipes, cannot form maccs keys for variant " << i << '\n';
        cerr << "Smiles '" << random_smiles << "'\n";
        return keep_going_after_test_failure;
      }

      if (do_comparison_between_key_values (m, initial_smiles, mtmp, random_smiles, keys, tkeys, i, problematic_keys))
        continue;

      failures_this_molecule++;
      if (! keep_going_after_test_failure)
      {
        molecules_with_test_failures++;
        return 0;
      }
      break;
    }
  }

  if (failures_this_molecule)
  {
    molecules_with_test_failures++;
    cerr << m.smiles() << ' ' << m.name() << " FAILED\n";
    for (int i = 0; i < problematic_keys.number_elements(); ++i)
    {
      cerr << " kfail " << problematic_keys[i] << '\n';
    }
  }

  if (report_progress())
    cerr << "Tested " << molecules_read << " molecules, " << molecules_with_test_failures << " molecules with failures\n";

  return 1;
}

static int
maccskeys (Molecule & m, 
           int * keys)
{
  int rc = mk (m, keys);

  if (rc && ntest)
    rc = test_maccskeys (m, keys);

  if (0 == rc)   // failure
    return 0;

  if (nullptr != accumulators)
  {
    int non_zero_keys = 0;
    for (int i = 0; i < number_maccs_keys; i++)
    {
      accumulators[i].extra (keys[i]);
      if (keys[i])
        non_zero_keys++;
    }

    keys_hit[non_zero_keys]++;
  }

  return rc;
}

static int
write_fixed_size_counted_fingerprint (const int * keys,
                                      IWString & output)
{
//int nset = count_non_zero_occurrences_in_array (keys + suppress_uninteresting_keys, nbits);

  output << fixed_size_counted_fingerprint_tag;
  (void) append_fixed_size_counted_fingerprint(keys, nbits, -1, -1, output);
  output << ">\n";

  return output.good();
}

static int
write_non_colliding_fingerprint (const int * keys,
                                 IWString & output)
{
  IWString dyascii;
  non_colliding_counted_fingerprint_daylight_representation(nbits, keys, dyascii);

  output << non_colliding_fingerprint_tag << dyascii << ">\n";
  
  return output.good();
}

/*
  Common code for writing any fingerprints that need to be written
*/

static int
write_fingerprints (int * keys,
                    IWString & output)
{
  if (fingerprint_tag.length())
  {
    IW_Bits_Base fp;

    (void) fp.construct_from_array_of_ints(keys, nbits);

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output << fingerprint_tag << tmp << ">\n";
  }

  if (level2_fingerprint_tag.length())
  {
    mk.set_level_2_fingerprint(keys);

    IW_Bits_Base fp(nbits);

    fp.construct_from_array_of_ints(keys, nbits);

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output << level2_fingerprint_tag << tmp << ">\n";
  }

  if (non_colliding_fingerprint_tag.length())
    write_non_colliding_fingerprint(keys, output);

  if (fixed_size_counted_fingerprint_tag.length())
    write_fixed_size_counted_fingerprint(keys, output);

  return output.good();
}

static void
preprocess (Molecule & m)
{
  m.remove_all(1);    // explicit hydrogens mess things up

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (strip_to_largest_fragment && m.number_fragments() > 1)
    m.reduce_to_largest_fragment();

  if (remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  if (revert_all_directional_bonds_to_non_directional)
    m.revert_all_directional_bonds_to_non_directional();

  // Ensure computed - to avoid incremental updates.
  m.recompute_distance_matrix();

  return;
}

static int
maccskeys (Molecule & m, 
           IWString_and_File_Descriptor & output)
{
  preprocess(m);

  int keys[256];     // make long enough for fingerprints however long they are

  set_vector(keys, 256, 0);

  int rc = maccskeys(m, keys);

  if (ntest)
    return rc;

  if (write_key_values)
  {
    if (! do_output(m, keys, output))
      return 0;
  }

  if (fingerprint_tag.length() || level2_fingerprint_tag.length() || non_colliding_fingerprint_tag.length() || fixed_size_counted_fingerprint_tag.length())
  {
    output << "$SMI<" << m.smiles() << ">\n";
    output << "PCN<" << m.name() << ">\n";
    write_fingerprints(keys, output);

    if (m.number_records_text_info())    // let's hope it is in TDT format!
      m.write_extra_text_info(output);

    output << "|\n";
  }

  return 1;
}

/*
  We function like Daylight's fingerprints utility
*/

static int
maccskeys_filter (iwstring_data_source & input,
                  int * keys,
                  IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with("$SMI<"))
    {
      buffer.remove_leading_chars(5);
      buffer.chop();

      Molecule m;
      if (! m.build_from_smiles(buffer))
      {
        cerr << "Yipes! cannot parse smiles '" << buffer << "'\n";
        return 0;
      }

      preprocess(m);

      if (! maccskeys(m, keys))
      {
        cerr << "Cannot compute keys for molecule at line " << input.lines_read() << '\n';
        cerr << buffer << '\n';
        return 0;
      }

      output << "$SMI<" << m.smiles() << ">\n";
      write_fingerprints(keys, output);
    }
    else if ('|' == buffer)
    {
      output << "|\n";

      set_vector(keys, nbits, 0);
    }
    else
    {
      output << buffer << '\n';
    }

    MaybeFlush(output);
  }

  return 1;
}

static int
maccskeys (data_source_and_type<Molecule> & input,
           IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while ((nullptr != (m = input.next_molecule())))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1)
      cerr << molecules_read << ' ' << m->name() << '\n';

    preprocess(*m);

    if (! maccskeys(*m, output))
      return 0;

    MaybeFlush(output);
  }

  return molecules_read;
}

static int
maccskeys_filter (const char * fname,
                  IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open stream file '" << fname << "'\n";
    return 0;
  }

  int * tmp = new_int(256);    // make it large - actually, no need to malloc it
  std::unique_ptr<int[]> free_tmp(tmp);

  return maccskeys_filter(input, tmp, output);
}

static int
maccskeys (const char * fname,
           FileType input_type,
           IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 2;
  }

  return maccskeys(input, output);
}

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
  cerr << "Usage: " << prog_name << " <options> <file1> <file2> ...\n";
  cerr << "  -t <number>    generate <number> random smiles for testing\n";
  cerr << "  -t rpt=<r>     report test progress every <r> molecules\n";
  cerr << "  -k             keep going after test failure\n";
  cerr << "  -S <fname>     gather and report statistics on key distributions\n";
  cerr << "  -n             no output of key values (only recognised with -s option)\n";
  cerr << "  -Y <fname>     write keys in Daylight form to <fname>\n";
  cerr << "  -J <tag>       tag for fingerprint data (default " << fingerprint_tag << ")\n";
  cerr << "  -J LEVEL2=XXX  use XXX as the level 2 tag (features with multiple occurrences)\n";
  cerr << "  -h             exclude hydrogens on aromatic nitrogen atoms\n";
  cerr << "  -b <nbits>     number of bits to produce (default 192)\n";
  cerr << "  -e             include the number of bits set as a descriptor\n";
  cerr << "  -x             remove cis trans bonds\n";
  cerr << "  -f             work as a filter\n";
  cerr << "  -i <type>      specify input type\n";
  cerr << "  -X ...         miscellaneous options, enter '-X help' for info\n";
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  cerr << "  -l             strip to largest fragment\n";
  cerr << "  -c             remove chiral centres\n";
  cerr << "  -E <el>        process element <el>\n";
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -v             verbose output\n";

  exit (rc);
}
static void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush         flush output after each molecule\n";

  ::exit(0);
}

static int
maccskeys (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "Y:b:nS:K:kt:J:fFvi:A:Clg:E:ehxcX:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(3);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      usage(6);
    }
  }

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot parse -E option\n";
      usage(18);
    }
  }

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose))
    {
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('b'))
  {
    if (! cl.value('b', number_maccs_keys) || number_maccs_keys <= 0)
    {
      cerr << "The number of maccs keys to produce (-b) must be a valid +ve number\n";
      usage(2);
    }

    if (! mk.set_nbits(number_maccs_keys))
    {
      cerr << "Invalid input for number MACCS keys " << number_maccs_keys << '\n';
      return 2;
    }

    if (verbose)
      cerr << "Will produce " << number_maccs_keys << " bits\n";
  }

  iwdigits.set_include_leading_space(1);
  iwdigits.initialise(255);

  if (0 == cl.number_elements())
  {
    cerr << "Must specify files to process on the command line\n";
    usage(5);
  }

  if (cl.option_present('t'))
  {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++))
    {
      if (t.starts_with("rpt="))
      {
        int rp;
        t.remove_leading_chars(4);
        if (! t.numeric_value(rp) || rp < 1)
        {
          cerr << "The report progress -t qualifier must be a whole +ve number\n";
          usage(4);
        }

        report_progress.set_report_every(rp);
      }
      else if ("norpt" == t)
      {
        report_test_failures = 0;
        if (verbose)
          cerr << "Will not report test failures\n";
      }
      else if (t.numeric_value(ntest) && ntest > 0)
      {
        if (verbose)
          cerr << "Will perform " << ntest << " random permutation tests on each molecule\n";
      }
      else
      {
        cerr << "Unrecognised -t qualifier '" << t << "'\n";
        usage(5);
      }
    }

    write_header = 0;
    set_input_aromatic_structures(1);
  }

  if (cl.option_present('k'))
  {
    if (! cl.option_present('t'))
    {
      cerr << "The -k option only makes sense with the -t option\n";
      usage(11);
    }

    keep_going_after_test_failure = 1;
    if (verbose)
      cerr << "Will continue execution even after a test failure\n";
  }

  IWString_and_File_Descriptor stream_for_statistics;

  if (cl.option_present('S'))
  {
    const char * s = cl.option_value('S');

    if (! stream_for_statistics.open(s))
    {
      cerr << "Cannot open stream for statistics '" << s << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Statistics on key frequency will be written to '" << s << "'\n";

    accumulators = new int_accumulator[number_maccs_keys + 1];
  }

  if (cl.option_present('n'))
  {
    if (0 == ntest && (! cl.option_present('J')) && 0 == stream_for_statistics.is_open())
    {
      cerr << "The -n option specified, but no other output\n";
      cerr << "Specify either -d or -J for some other output\n";
      usage(83);
    }

    write_key_values = 0;
    write_header = 0;
    if (verbose)
      cerr << "Writing keys to cout suppressed\n";
  }

  if (cl.option_present('h'))
  {
    mk.set_aromatic_nitrogens_do_not_have_hydrogens(1);
    if (verbose)
      cerr << "Will not consider Hydrogen atoms on aromatic Nitrogens\n";
  }

  if (cl.option_present('x'))
  {
    revert_all_directional_bonds_to_non_directional = 1;

    if (verbose)
      cerr << "Will remove directional bonds\n";
  }

  if (cl.option_present('c'))
  {
    remove_all_chiral_centres = 1;

    if (verbose)
      cerr << "Will remove all chiral centres\n";
  }

  if (cl.option_present('C'))
  {
    cerr << "The -C option is no longer recognised\n";
  }

  if (cl.option_present('F'))
  {
    cerr << "The -F option is no longer recognised\n";
  }

  if (cl.option_present('K'))
  {
    cerr << "The -K option is no longer recognised\n";
  }

  if (cl.option_present('e'))
  {
    append_nset = 1;

    if (verbose)
      cerr << "The number of bits set will be included as a descriptor\n";
  }

  if (cl.option_present('l'))
  {
    strip_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will strip input molecules to largest fragment\n";
  }


  if (cl.option_present('J'))
  {
    if (! cl.option_present('n'))
    {
      write_key_values = 0;
      write_header = 0;
      if (verbose)
        cerr << "Fingerprint output, array output suppressed\n";
    }

    int i = 0;
    const_IWSubstring j;
    while (cl.value('J', j, i++))
    {
      if (j.starts_with("LEVEL2="))
      {
        level2_fingerprint_tag = j;
        level2_fingerprint_tag.remove_leading_chars(7);

        if (verbose)
          cerr << "Level 2 fingerprint tag '" << level2_fingerprint_tag << "'\n";

        if (! level2_fingerprint_tag.ends_with('<'))
          level2_fingerprint_tag += '<';
      }
      else if (j.starts_with("NC="))
      {
        non_colliding_fingerprint_tag = j;
        non_colliding_fingerprint_tag.remove_up_to_first('=');

        if (verbose)
          cerr << "Non colliding fingerprint tag '" << non_colliding_fingerprint_tag << "'\n";

        if (! non_colliding_fingerprint_tag.ends_with('<'))
          non_colliding_fingerprint_tag += '<';
      }
      else if (j.starts_with("FC="))
      {
        fixed_size_counted_fingerprint_tag = j;
        fixed_size_counted_fingerprint_tag.remove_leading_chars(3);

        if (verbose)
          cerr << "Fixed size counted fingerprint tag " << fixed_size_counted_fingerprint_tag << "'\n";

        if (! fixed_size_counted_fingerprint_tag.ends_with('<'))
          fixed_size_counted_fingerprint_tag += '<';
      }
      else
      {
        fingerprint_tag = j;
        if (verbose)
          cerr << "Fingerprints written with dataitem '" << fingerprint_tag << "'\n";

        if (! fingerprint_tag.ends_with('<'))
          fingerprint_tag << '<';
      }
    }

    if (0 == level2_fingerprint_tag.length() && 0 == fingerprint_tag.length())
      ;
    else if (level2_fingerprint_tag == fingerprint_tag)
    {
      cerr << "Both default and level 1 fingerprint tags identical. Impossible\n";
      usage(11);
    }

    nbits = number_maccs_keys;
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  IWString_and_File_Descriptor output(1);

  FileType input_type = FILE_TYPE_INVALID;

  int rc = 0;

// If we are a filter, we take a different path

  if (cl.option_present('f'))
  {
    if (0 == fingerprint_tag.length() && 0 == non_colliding_fingerprint_tag.length())
    {
      fingerprint_tag = "FPMK<";
      cerr << "Fingerprint tag '" << fingerprint_tag << "'\n";
    }

    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! maccskeys_filter(cl[i], output))
      {
        cerr << "Error processing '" << cl[i] << "'\n";
        return 87;
      }
    }

    goto finished;
  }

  if (write_header)
  {
    if (! do_write_header(output))
    {
      return 8;
    }

  }

//int input_type;
  if (! process_input_type(cl, input_type))
  {
    cerr << "Cannot determine input type\n";
    usage(6);
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! maccskeys(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  finished:;

  if (verbose)
  {
    cerr << molecules_read << " molecules read\n";

    if (test_failures)
      cerr << molecules_with_test_failures << " molecules had " << test_failures << " test failures\n";
  }

  if (stream_for_statistics.is_open())
  {
    assert (nullptr != accumulators);

    stream_for_statistics << "Key statistics";
    if (0 == verbose)
      stream_for_statistics << " for " << molecules_read << " molecules";
    stream_for_statistics << '\n';

    for (int i = 0; i < number_maccs_keys; i++)
    {
      int_accumulator & a = accumulators[i];
      stream_for_statistics << "Key " << i << ' ' << a.n() << " molecules ";
      if (a.number_zero())
        stream_for_statistics << a.number_zero() << " misses";
      else
        stream_for_statistics << "min hits " << a.minval();

      stream_for_statistics << " max hits = " << a.maxval();
      if (a.n() > 1)
        stream_for_statistics << " ave = " << a.average();

      stream_for_statistics << '\n';
    }

    Accumulator_Int<int> keys_hit_per_molecule;
    for (int i = 0; i < keys_hit.number_elements(); i++)
    {
      if (keys_hit[i])
      {
        stream_for_statistics << keys_hit[i] << " molecules had " << i << " keys hit\n";
        keys_hit_per_molecule.extra(i, keys_hit[i]);
      }
    }

    if (keys_hit_per_molecule.n() > 1)
      stream_for_statistics << keys_hit_per_molecule.average() << " average bits hit per molecule\n";
  
// 2,   // key 6  1.916748672

    stream_for_statistics << "static int level2_threshold [] = {\n";
    for (auto i = 0; i < number_maccs_keys; ++i)
    {
      if (0 == accumulators[i].n())
        stream_for_statistics << "  0";
      else
        stream_for_statistics << "  " << static_cast<int>(accumulators[i].average() + 0.4999999);

      if (number_maccs_keys - 1 != i)
        stream_for_statistics << ',';

      stream_for_statistics << "     // key " << i << "   " << static_cast<float>(accumulators[i].average()) << "\n";
    }

    stream_for_statistics << "};\n";

    stream_for_statistics.close();
  }

  output.flush();

  if (nullptr != accumulators)
    delete [] accumulators;

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = maccskeys(argc, argv);

  return rc;
}
