//  Linear Fingerprint generation

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

const char * prog_name = nullptr;

static int verbose = 0;

static int function_as_filter = 0;

static int accumulate_statistics = 0;

static int molecules_read = 0;

static int tdts_read = 0;

static int nbits = 2048;

static int append_nset = 0;

/*
  We can write as 0's and 1's too
*/

static int write_as_numbers = 0;

static int write_non_zero_bits = 0;

/*
  We can collect statistics on which bits are hit
*/

static extending_resizable_array<int> bit_count;

/*
  Also on the number of bits hit per fingerprint
*/

static Accumulator_Int<int> bits_hit;

static int ntest = 0;

static int test_failures = 0;

static IWString fingerprint_tag("FPIW<");

static IWString descriptor_stem("chgfp_lfp");

static IWString level2_fingerprint_tag;

static IWString input_smiles_tag("$SMI<");

static int write_sparse_fingerprints = 0;

static int bits_with_overflowing_counts = 0;

/*
  The database can only hold 255 hits, so we can optionally cap
  the number of hits reported in a descriptor file
*/

static int max_hits = 0;

/*
  We can de-emphasise carbon atoms by specifying a minimum number of
  heteroatoms on the ends of the paths.
*/

static int min_heteroatoms_at_ends = 0;

static Element_Transformations element_transformations;

static Charge_Assigner charge_assigner;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int munge_to_single_bonds = 0;
static int munge_to_graph = 0;
static int munge_to_carbon = 0;
static int need_to_do_munging = 0;

static Atom_Typing_Specification atom_typing_specification;

#ifdef COUNT_TIMES_ATOM_IN_PATH
static int times_in_path_optimsation = 0;
#endif

static IWDigits iwdigits;
static IWDigits iwdigits_count;

static int flush_after_each_molecule = 0;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Produces hashed linear path fingerprints\n";

  cerr << prog_name << " <options> <input_file1> <input_file2> ... > newfile\n";
//cerr << "  -s min/max     min and max path length (default: 0/7) - obsolete\n";
  cerr << "  -r <min>       minimum path length (default 0)\n";
  cerr << "  -R <min>       maximum path length (default 7)\n";
  cerr << "  -k <number>    min number of heteroatoms at ends of paths: 0 (default) 1 or 2\n";
  cerr << "  -o             compute paths (heteroatom terminated) with CH2 groups omitted\n";
  cerr << "  -o             second -o option to CH2 suppressed paths for all paths\n";
  cerr << "  -T <ntest>     perform <ntest> random tests for each molecule\n";
  cerr << "  -B             accumulate statistics\n";
  cerr << "  -c <nbits>     creation size of fingerprint, bits (default: " << nbits << ")\n";
  cerr << "  -J <tag>       specify fingerprint tag (default '" << fingerprint_tag << "')\n";
  //cerr << "  -V             add Version TDT to output\n";
  cerr << "  -O <pathlen>   add formal charge bits to paths up to length <pathlen>\n";
  cerr << "  -e             include the number of bits set as a descriptor\n";
  cerr << "  -M ...         munge options, enter '-M help' for info\n";
  cerr << "  -f             function as filter\n";
  cerr << "  -S <tag>       input smiles tag (default " << input_smiles_tag << ")\n";
  cerr << "  -a             write as array (0,1). Repeat to write counts\n";
  cerr << "  -q             write non zero bit numbers, repeate to include counts\n";
  cerr << "  -p <prefix>    use <prefix> as the descriptor prefix\n";
  cerr << "  -P <atype>     specify atom typing\n";
  cerr << "  -x <nhits>     when producing descriptors, limit to <nhits> occurrences\n";
  cerr << "  -Y ...         standard fingerprint generation options, '-Y help' for info\n";
  cerr << "  -U ...         miscellaneous and obscure options, enter '-U help' for info\n";
#ifdef COUNT_TIMES_ATOM_IN_PATH
  cerr << "  -z <n>         perform <n> iterations of atom in path counter optimisation\n";
#endif
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -t E1=E2       element transformation options, enter '-t help' for info\n";
  (void) display_standard_charge_assigner_options (cerr, 'N');
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  cerr << "  -E ...         standard element options, enter '-E help' for info\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -v             verbose output\n";

  exit (rc);
}

static void
convert_to_single_bonds (Molecule & m)
{
  const int nb = m.nedges();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);
    if (! b->is_single_bond())
      continue;

    m.set_bond_type_between_atoms(b->a1(), b->a2(), SINGLE_BOND);
  }

  return;
}
static int
do_append_sparse_fingerprints (const IWMFingerprint & fp,
                               const IWString & tag,
                               IWString & output_buffer)
                               
{
  Sparse_Fingerprint_Creator sfp;

  const unsigned int n = iwmfingerprint_nbits();

  const int * v = fp.vector();

  for (unsigned int i = 0; i < n; i++)
  {
    if (0 == v[i])
      continue;

    if (v[i] <= 255)
    {
      sfp.hit_bit(i, v[i]);
      continue;
    }

//  use extra bits if we have overflowed 255

    bits_with_overflowing_counts++;

    int x = v[i];
    for (unsigned int b = i; ; b += nbits)
    {
      sfp.hit_bit(b, 255);

      x -= 255;

      if (x <= 255)
      {
        sfp.hit_bit(b + nbits, x);
        break;
      }
    }
  }

  if (verbose)
    bits_hit.extra(sfp.nbits());

  IWString tmp;

  sfp.daylight_ascii_form_with_counts_encoded(tmp);

  output_buffer << tag << tmp << ">\n";

  return 1;
}

static int
append_fingerprint_to_output_buffer (const IWMFingerprint & fp,
                                     const const_IWSubstring & tag,
                                     IWString & output_buffer)
{
  if (write_sparse_fingerprints)
    return do_append_sparse_fingerprints(fp, tag, output_buffer);

  output_buffer << tag;

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);

  output_buffer << tmp << ">\n";

  return output_buffer.length();
}


static int
do_munging(Molecule & m, IWString & output_buffer)
{
  if (munge_to_single_bonds)
  {
    IWMFingerprint fp;

    convert_to_single_bonds(m);

    fp.construct_fingerprint(m);

    append_fingerprint_to_output_buffer(fp, "SPMSB", output_buffer);
  }

  if (munge_to_graph)
  {
    m.change_to_graph_form();

    IWMFingerprint fp;

    fp.construct_fingerprint(m);

    append_fingerprint_to_output_buffer(fp, "FPMGR", output_buffer);
  }

  if (munge_to_carbon)
  {
    if (! munge_to_single_bonds)
      convert_to_single_bonds(m);

    int matoms = m.natoms();
    for (int i = 0; i < matoms; i++)
    {
      atomic_number_t z = m.atomic_number(i);
      if (6 != z)
        m.set_atomic_number(i, 6);
    }

    IWMFingerprint fp;
    fp.construct_fingerprint(m);

    append_fingerprint_to_output_buffer(fp, "FPMC", output_buffer);
  }

  return 1;
}

static int
do_write_non_zero_bits(const Molecule & m,
                       const IWMFingerprint & fp,
                       IWString & output_buffer)
{
  if (m.name().empty()) {
    output_buffer << '.';
  } else {
    append_first_token_of_name(m.name(), output_buffer);
  }

  int * b = const_cast<int *>(fp.vector());     // loss of const OK

  const int nb = fp.nbits();

  for (int i = 0; i < nb; ++i)
  {
    if (0 == b[i])
      continue;

    iwdigits.append_number(output_buffer, i);
    if (write_non_zero_bits > 1)
      iwdigits_count.append_number(output_buffer, b[i]);
  }

  return 1;
}

static int
do_output_numbers (const Molecule & m,
                   const IWMFingerprint & fp,
                   IWString & output_buffer)
{
  if (m.name().empty()) {
    output_buffer << '.';
  } else {
    append_first_token_of_name(m.name(), output_buffer);
  }

  output_buffer << ' ';

  int rc;
  if (1 == write_as_numbers)
    rc = fp.write_as_zero_and_one(output_buffer);
  else
    rc = fp.write_count(output_buffer);

  if (append_nset)
  {
    output_buffer << ' ' << fp.nset();
  }

  output_buffer << '\n';

  return rc;
}

static int
do_write_level2_fingerprint (IWMFingerprint & fp,
                   IWString & output_buffer)
{
  return 1;    // not implemented yet

  int * tmp = new int[nbits]; std::unique_ptr<int[]> free_tmp(tmp);

//fp.level2_fingerprint (tmp);

  IW_Bits_Base b;
  b.construct_from_array_of_ints(tmp, nbits);

  append_fingerprint_to_output_buffer(fp, level2_fingerprint_tag, output_buffer);

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
do_output (Molecule & m,
           IWMFingerprint & fp,
           IWString & output_buffer)
{
  if (write_as_numbers)
    return do_output_numbers(m, fp, output_buffer);

  if (write_non_zero_bits)
    return do_write_non_zero_bits(m, fp, output_buffer);

  output_buffer << "$SMI<" << m.smiles() << ">\n";

  append_fingerprint_to_output_buffer(fp, fingerprint_tag, output_buffer);

  if (level2_fingerprint_tag.length())
    do_write_level2_fingerprint(fp, output_buffer);

  output_buffer << "PCN<" << m.name() << ">\n";

  return 1;
}

static int
iwfp_tester (Molecule & m,
             const IWMFingerprint & initial_fingerprint)
{
///nt matoms = m.natoms();

  int rc = 1;

  for (int i = 0; i < ntest; i++)
  {
    IWString rsmiles = m.random_smiles();

    Molecule rmol;
    if (! rmol.build_from_smiles(rsmiles))
    {
      cerr << "Yipes, cannot parse smiles '" << rsmiles << "'\n";
      rc = 0;
      break;
    }

//  assert (matoms == rmol.natoms());

    IWMFingerprint fp;
    if (! fp.construct_fingerprint(rmol))
    {
      cerr << "Very strange, construct_fingerprint failed\n";
      rc = 0;
      break;
    }

    if (fp == initial_fingerprint)
      continue;
    
    int bits_the_same = 0;     // non zero bits only

    const int * expected_value = initial_fingerprint.vector();
    const int * myresult       = fp.vector();

    for (int j = 0; j < nbits; j++)
    {
      if (expected_value[j] == myresult[j])
      {
        if (myresult[j] > 0)
          bits_the_same++;
        continue;
      }

      if (rc)      // write header on first error
        cerr << "Bit value mismatch for random smiles " << i << " '" << rsmiles << "'\n";

      cerr << "Bit value difference, bit " << j << " expected " << expected_value[j] << " got " << myresult[j] << '\n';
      rc = 0;
    }

    if (0 == rc)
    {
      cerr << bits_the_same << " non zero bits the same\n";
      test_failures++;
      break;
    }
  }

  return rc;
}

static int
write_header (IWString & output_buffer,
              int write_as_numbers)
{
  if (0 == write_as_numbers)   // no header needed
    return 1;

  if (0 == verbose || 0 == write_as_numbers)    // don't say anything
    ;
  else if (1 == write_as_numbers)
    cerr << "Will write fingerprints as 0's and 1's\n";
  else if (2 == write_as_numbers)
    cerr << "Will write path counts\n";

  IWString header;

  for (int i = 0; i < nbits; i++)
  {
    header << ' ' << descriptor_stem << i;
  }

  if (append_nset) 
    header << ' ' << descriptor_stem << "nset";

  if (write_as_numbers)
    output_buffer << "Name" << header << '\n';

  return 1;
}

static void
do_accumulate_statistics (const IWMFingerprint & fp)
{
  const int * bvector = fp.vector();
  assert (nullptr != bvector);

  int bits_set_this_fingerprint = 0;

  for (int i = 0; i < nbits; i++)
  {
    if (bvector[i])
    {
      bit_count[i]++;
      bits_set_this_fingerprint++;
    }
  }

  bits_hit.extra(bits_set_this_fingerprint);

  if (0 == bits_set_this_fingerprint)
    cerr << "Yipes, sets no bits\n";

  return;
}

static int
iwfp (Molecule & m, IWMFingerprint & fp)
{
  if (min_heteroatoms_at_ends)
    fp.set_min_heteroatoms_at_path_ends(min_heteroatoms_at_ends);

  if (atom_typing_specification.active())
    fp.construct_fingerprint(m, atom_typing_specification, nullptr);
  else
    fp.construct_fingerprint(m);

  return 1;
}

static int
iwfp (Molecule & m,
      IWString_and_File_Descriptor & output)
{
  IWMFingerprint fp;

#ifdef COUNT_TIMES_ATOM_IN_PATH
  if (times_in_path_optimsation)
    fp.set_cycles_adjust_times_used(times_in_path_optimsation);
#endif

  if (! iwfp(m, fp))
    return 0;

  if (ntest)
    return iwfp_tester(m, fp);

  if (accumulate_statistics)
    do_accumulate_statistics(fp);

  if (max_hits > 0)
    fp.truncate_to_max_hits(max_hits);

  if (! do_output(m, fp, output))
    return 0;

  if (write_as_numbers)
    return 1;

  if (need_to_do_munging)
    do_munging(m, output);

  output << "|\n";

  MaybeFlush(output);

  return 1;
}

static int
iwfp (const const_IWSubstring & smi,
      IWString_and_File_Descriptor & output)
{
  Molecule m;
  if (! m.build_from_smiles(smi))
  {
    cerr << "Cannot parse smiles '" << smi << "'\n";
    return 0;
  }

  IWMFingerprint fp;
  if (! iwfp(m, fp))
    return 0;

  append_fingerprint_to_output_buffer(fp, fingerprint_tag, output);

  if (need_to_do_munging)
    do_munging(m, output);

  if (ntest)
  {
    if (! iwfp_tester(m, fp))
      return 0;
  }

  return 1;
}

static int
iwfp (iwstring_data_source & input,
      IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if ("|" == buffer)
    {
      tdts_read++;
      continue;
    }

    if (! buffer.starts_with(input_smiles_tag))
      continue;

    buffer.remove_leading_chars(input_smiles_tag.length());
    buffer.chop();

    if (! iwfp(buffer, output))
    {
      cerr << "Fatal error processing smiles '" << buffer << "', line " << input.lines_read() << '\n';
      return 0;
    }
  }

  return 1;
}

static int
iwfp (const char * fname,
      IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open filter source '" << fname << "'\n";
    return 0;
  }

  return iwfp(input, output);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (element_transformations.active())
    (void) element_transformations.process(m);

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (charge_assigner.active())
    (void) charge_assigner.process(m);

  return;
}

static int
iwfp (data_source_and_type<Molecule> & input,
      IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    if (! iwfp(*m, output))
      return 0;

    MaybeFlush(output);
  }

  return 1;
}

static int
iwfp (const char * fname, FileType input_type, 
      IWString_and_File_Descriptor & output)
{
  assert (nullptr != fname);

  if (function_as_filter)
    return iwfp(fname, output);

  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return iwfp(input, output);
}

static void
display_standard_munging_options (std::ostream & os, char mflag)
{
  os << "  -" << mflag << " single       munge to all single bonds\n";
  os << "  -" << mflag << " graph        munge to graph form\n";
  os << "  -" << mflag << " carbon       munge to all carbon\n";

  exit(1);
}

static void
DisplayDashUOptions(std::ostream& output) {
  output << " -U flush             flush output after each molecule\n";

  ::exit(0);
}

static int
iwfp (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vVc:E:aA:k:t:J:s:i:T:t:BfoN:O:g:lS:d:eM:x:Y:p:P:r:R:qz:U:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');


  if (! process_elements(cl))
    usage(2);

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose > 1))
    {
      cerr << "Cannot process aromaticity options (-A)\n";
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('t'))
  {
    if (! element_transformations.construct_from_command_line(cl, verbose, 't'))
    {
      cerr << "Cannot parse -t option\n";
      usage(11);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('N'))
  {
    if (! charge_assigner.construct_from_command_line(cl, verbose, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage(33);
    }
  }

  if (cl.option_present('N') && ! cl.option_present('O') && ! cl.option_present('y'))
  {
    set_include_formal_charge_bits(0);
    cerr << "Zero path length formal charge bits applied\n";
  }

  if (cl.option_present('O'))
  {
    int oo;
    if (! cl.value('O', oo) || oo < 0)
    {
      cerr << "The -O option must be followed by a whole non-negative number\n";
      usage(41);
    }

    set_include_formal_charge_bits(oo);
    if (verbose)
      cerr << "Formal charge bits set for paths up to length " << oo << '\n';
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', nbits) || nbits < 8)
    {
      cerr << "The -c option must be followed by a whole positive number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will produce fingerprints of " << nbits << " bits\n";
  }

  if (cl.option_present('p'))
  {
    cl.value('p', descriptor_stem);

    if (verbose)
      cerr << "Descriptors created with prefix '" << descriptor_stem << "'\n";
  }

  set_iwmfingerprint_nbits(nbits);

  if (cl.option_present('Y'))
  {
    if (! parse_misc_fingerprint_options(cl, 'Y', verbose))
      return 5;
  }

  if (cl.option_present('f'))
  {
    function_as_filter = 1;
    if (verbose)
      cerr << "Will function as a filter\n";

    if (cl.option_present('i'))
      cerr << "The -i option is ignored when functioning as a filter\n";
  }

  if (cl.option_present('S'))
  {
    if (! cl.option_present('f'))
    {
      cerr << "The -S option (input smiles tag) only makes sense with the -f (work as filter) option\n";
      usage(32);
    }

    cl.value('S', input_smiles_tag);
    if (verbose)
      cerr << "Input smiles tag '" << input_smiles_tag << "'\n";

    if (! input_smiles_tag.ends_with('<'))
      input_smiles_tag.add('<');
  }

  if (cl.option_present('o'))
  {
    int tmp = cl.option_count('o');

    set_omit_ch2(tmp);
    if (0 == verbose)
      ;
    else if (1 == tmp)
      cerr << "Will set bits for heteroatom terminated paths with omitted CH2 groups\n";
    else
      cerr << "Will set bits for all paths with omitted CH2 groups\n";
  }

#ifdef COUNT_TIMES_ATOM_IN_PATH
  if (cl.option_present('z'))
  {
    if (! cl.value('z', times_in_path_optimsation) || times_in_path_optimsation < 1)
    {
      cerr << "The number of times in fingerprints optimisation steps performed (-z) must be a whole +ve number\n";
      usage(1);
    }
  }
#endif

  if (cl.option_present('B'))
  {
    accumulate_statistics = 1;
    if (verbose)
      cerr << "Will accumulate bit statistics\n";

    bit_count.resize(nbits);
  }

  if (cl.option_present('k'))
  {
    if (! cl.value('k', min_heteroatoms_at_ends) || min_heteroatoms_at_ends < 0 || min_heteroatoms_at_ends > 2)
    {
      cerr << "The -k option must be followed by a whole number between 0 and 2\n";
      usage(3);
    }

    if (verbose)
      cerr << "Paths must have at least " << min_heteroatoms_at_ends << " heteroatoms at their ends\n";
  }

  if (cl.option_present('M'))
  {
    const_IWSubstring m;
    int i = 0;
    while (cl.value('M', m, i++))
    {
      if ("single" == m)
      {
        munge_to_single_bonds = 1;
        need_to_do_munging++;

        if (verbose)
          cerr << "Will munge to single bonds\n";
      }
      else if ("graph" == m)
      {
        munge_to_graph = 1;
        need_to_do_munging++;

        if (verbose)
          cerr << "Will munge to graph form\n";
      }
      else if ("carbon" == m)
      {
        munge_to_carbon = 1;
        need_to_do_munging++;

        if (verbose)
          cerr << "Will munge to carbon form\n";
      }
      else if ("help" == m)
      {
        display_standard_munging_options(cerr, 'M');
      }
      else
      {
        cerr << "Unrecognised munge specifier '" << m << "'\n";
        display_standard_munging_options(cerr, 'M');
      }
    }
  }

  if (cl.option_present('P'))
  {
    const const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p))
    {
      cerr << "Invalid atom type specification '" << p << "'\n";
      return 3;
    }
  }

  if (cl.option_present('U')) {
    const_IWSubstring u;
    for (int i = 0; cl.value('U', u, i); ++i) {
      if (u == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule processed\n";
        }
      } else if (u == "help") {
        DisplayDashUOptions(cerr);
      } else {
        cerr << "Unrecognised -U qualifier '" << u << "'\n";
        DisplayDashUOptions(cerr);
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_filter)     // no need to check input type
    ;
  else if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('T'))
  {
    if (! cl.value('T', ntest) || ntest < 1)
    {
      cerr << "The -T option must be followed by a whole positive number\n";
      usage(8);
    }

    if (verbose)
      cerr << "Will test " << ntest << " random smiles for each molecule\n";
  }

  if (cl.option_present('J'))
  {
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
      }
      else
      {
        fingerprint_tag = j;
        if (verbose)
          cerr << "Fingerprints written with dataitem '" << fingerprint_tag << "'\n";
      }
    }

    if (! fingerprint_tag.ends_with('<'))
      fingerprint_tag << '<';

    if (0 == level2_fingerprint_tag.length())
      ;
    else if (level2_fingerprint_tag.ends_with('<'))
      ;
    else
      level2_fingerprint_tag << '<';

    if (level2_fingerprint_tag == fingerprint_tag)
    {
      cerr << "Both default and level 1 fingerprint tags identical. Impossible\n";
      usage(11);
    }

    if (fingerprint_tag.starts_with("NC"))
    {
      set_form_bit_vector_during_fingerprint_creation(0);
      write_sparse_fingerprints = 1;
      if (verbose)
        cerr << "Will write fingerprints as non-colliding form\n";
    }
  }

  if (cl.option_present('s'))
  {
    if (! parse_path_length_options(cl, 's', verbose))
    {
      cerr << "Cannot parse path length options, -s\n";
      usage(7);
    }
  }

  if (cl.option_present('r'))
  {
    int p;
    if (! cl.value('r', p) || p < 0)
    {
      cerr << "The minimum path length option (-r) must be a whole +ve number\n";
      usage(3);
    }

    set_min_path_length(p);
    if (verbose)
      cerr << "Min path length set to " << p << '\n';
  }

  if (cl.option_present('R'))
  {
    int p;
    if (! cl.value('R', p) || p < 0)
    {
      cerr << "The maximum path length option (-R) must be a whole +ve number\n";
      usage(3);
    }

    set_max_path_length(p);
    if (verbose)
      cerr << "Max path length set to " << p << '\n';
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', max_hits) || max_hits < 1)
    {
      cerr << "The maximum number of hits per bit option (-x) must be a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Max hits per bit set to " << max_hits << '\n';
  }

  if (write_sparse_fingerprints && iwmfingerprint_nbits() <= 2048)
    cerr << "Warning, creating sparse fingerprints, but only " << iwmfingerprint_nbits() << " bits\n";

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  write_as_numbers = cl.option_count('a');

  write_non_zero_bits = cl.option_count('q');

  if (write_as_numbers && write_non_zero_bits)
  {
    cerr << "The -a and -q options are mutually exclusive\n";
    usage(1);
  }

  if (cl.option_present('e')) {
    if (write_non_zero_bits == 0) {
      write_as_numbers = 1;
    }
    append_nset = 1;

    if (verbose)
      cerr << "Will include the number of bits set as a descriptor\n";
  }

  if (write_non_zero_bits)
  {
    iwdigits.set_include_leading_space(1);
    iwdigits.initialise(nbits);

    if (write_non_zero_bits > 1)
    {
      iwdigits_count.set_leading_string(":");
      iwdigits_count.initialise(256);
    }
  }

  IWString_and_File_Descriptor output(1);
  write_header(output, write_as_numbers);

  if (cl.option_present('V'))
  {
    if (verbose)
      cerr << "A version tree will be added\n";

    output << "IWFP<Version 2.0 compiled " << __DATE__ << ">\n";
    output << "|\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! iwfp(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Produced fingerprints for ";
    if (molecules_read)
      cerr << molecules_read;
    else if (tdts_read)
      cerr << tdts_read;
    else
      cerr << 0;
    cerr << " molecules\n";

    if (bits_with_overflowing_counts)
      cerr << bits_with_overflowing_counts << " bits overflowed 255 and generated overflow bits\n";
  }

  if (accumulate_statistics)
  {
    int zero_hits = 0;
    int always_hit = 0;

    for (int i = 0; i < nbits; i++)
    {
      if (0 == bit_count[i])
        zero_hits++;
      else if (molecules_read == bit_count[i])
        always_hit++;

      if (bit_count[i])
        cerr << "Bit " << i << " hit " << bit_count[i] << " times\n";
    }

    if (zero_hits)
      cerr << zero_hits << " bits never hit\n";

    if (always_hit)
      cerr << always_hit << " bits always hit\n";

    cerr << "Bits hit between " << bits_hit.minval() << " and " << bits_hit.maxval();
    if (bits_hit.n())
      cerr << " average " << bits_hit.average();
    cerr << '\n';
  }

#ifdef CHECK_COLLISIONS
  cerr << "Found " << collisions << " collisions\n";
#endif

  if (test_failures)
    cerr << "Found " << test_failures << " test failures\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = iwfp(argc, argv);

  return rc;
}
