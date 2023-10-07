// Compute linear path fingerprints.

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/linear_fingerprint.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"

#include "Molecule_Tools/fingerprint_writer.h"

using std::cerr;

const char * prog_name = nullptr;

int verbose = 0;

int molecules_read = 0;

Chemical_Standardisation chemical_standardisation;

Atom_Typing_Specification atom_typing;

int reduce_to_largest_fragment = 0;

static fingerprint_writer::FingerprintWriter fp_writer;

bool function_as_tdt_filter = false;

bool check_for_collisions = false;

// By default, non colliding fingerprints are generated.
// If set, fixed width fingerprints are produced.
int fixed_width = 0;

const IWString smiles_tag("$SMI<");
const IWString identifier_tag("PCN<");

// Potentially interesting statistics on molecules processed.

bool gather_molecule_statistics = false;
Accumulator_Int<uint> atom_count;
Accumulator_Int<uint> longest_path_acc;
extending_resizable_array<int> longest_path;

Accumulator_Int<uint> bits_set;

using linear_fingerprint::LinearFingerprintGenerator;

void
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
  cerr << "Computes linear path fingerprints\n";
  cerr << "  -r <rad>      minimum path length (def 0)\n";
  cerr << "  -R <rad>      maximum path length (def 7)\n";
  cerr << "  -P ...        atom type specification\n";
  cerr << "  -J <tag>      tag for fingerprints\n";
  cerr << "  -f            function as a TDT filter\n";
  cerr << "  -X <fname>    look for bits in <fname> and provide explanations\n";
  cerr << "  -B <fname>    write all bits found to <fname>\n";
  cerr << "  -y            check for bit collisions\n";
  cerr << "  -s            gather statistics on molecules processed\n";
  cerr << "  -c            produce isotopically labelled smiles with coverage\n";
  cerr << "  -x            allow linear paths can cross\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

void
Preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

void
GatherStatistics(Molecule & m) {
  atom_count.extra(m.natoms());

  const int lp = m.longest_path();

  longest_path_acc.extra(lp);
  longest_path[lp]++;
}

int
LinearFingerprint(Molecule & m,
                  const uint64_t* atype,
                  LinearFingerprintGenerator& linear_fp_gen,
                  IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfc;

  // A single atom molecule would produce zero bits, or two atoms in different fragments...
  if (! linear_fp_gen.Fingerprint(m, nullptr, atype, sfc)) {
    cerr << "LinearFingerprintGenerator:fingerprinting failed " << m.smiles() << " ignored\n";
  }

  if (gather_molecule_statistics)
    GatherStatistics(m);

  if (verbose)
    bits_set.extra(sfc.nbits());

  if (function_as_tdt_filter) {
  } else if (fp_writer.FingerprintTag().length() > 0) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  fp_writer.WriteFingerprint(m.name(), sfc, output);

  if (function_as_tdt_filter) {
  } else if (fp_writer.FingerprintTag().length() > 0) {
    output << "|\n";
  }

  return 1;
}

int
LinearFingerprint(Molecule & m,
                  LinearFingerprintGenerator& linear_fp_gen,
                  IWString_and_File_Descriptor & output)
{
  uint64_t * atype = new uint64_t[m.natoms()]; std::unique_ptr<uint64_t[]> free_atype(atype);

  if (! atom_typing.assign_atom_types(m, atype)) {
    cerr << "LinearFingerprint::cannot assign atom types " << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  return LinearFingerprint(m, atype, linear_fp_gen, output);
}

int
LinearFingerprint(data_source_and_type<Molecule> & input,
                LinearFingerprintGenerator& linear_fp_gen,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    Preprocess(*m);

    if (! LinearFingerprint(*m, linear_fp_gen, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
LinearFingerprintPipe(const_IWSubstring line,   // Note local copy
                LinearFingerprintGenerator& linear_fp_gen,
                IWString_and_File_Descriptor& output) {
  assert(line.ends_with('>'));
  assert(line.starts_with(smiles_tag));

  line.remove_leading_chars(smiles_tag.length());
  line.chop();

  Molecule m;

  if (! m.build_from_smiles(line)) {
    cerr << "LinearFingerprintPipe:invalid smiles '" << line << "'\n";
    return 0;
  }

  Preprocess(m);

  return LinearFingerprint(m, linear_fp_gen, output);
}

int
LinearFingerprintPipe(iwstring_data_source& input,
                LinearFingerprintGenerator& linear_fp_gen,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    output << line << '\n';

    if (! line.starts_with(smiles_tag))
      continue;

    if (! LinearFingerprintPipe(line, linear_fp_gen, output)) {
      cerr << "LinearFingerprintPipe:invalid input " << line << "' line " << input.lines_read() << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
LinearFingerprint(const char * fname, FileType input_type, 
                LinearFingerprintGenerator& linear_fp_gen,
                IWString_and_File_Descriptor & output)
{
  assert(nullptr != fname);

  if (function_as_tdt_filter) {
    iwstring_data_source input(fname);
    if (! input.good()) {
      cerr << "LinearFingerprint::cannot open filter " << fname << '\n';
      return 0;
    }

    return LinearFingerprintPipe(input, linear_fp_gen, output);
  }

  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return LinearFingerprint(input, linear_fp_gen, output);
}

int
LinearFingerprint(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "E:A:K:lg:i:J:P:vfr:R:ysB:cx");

  if (cl.unrecognised_options_encountered())
    usage(1);

  verbose = cl.option_count('v');
  
  (void) process_elements(cl);

  if (! process_standard_smiles_options(cl, verbose))
  {
    usage(7);
  }

  set_global_aromaticity_type(Daylight);

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    usage(8);
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  int min_length = 0;
  if (cl.option_present('r'))
  {
    if (cl.value('r', min_length) || min_length < 1) {
      cerr << "The min atom pair radius (-r) must be > 0\n";
      return 1;
    }

    if (verbose)
      cerr << "Will only fingerprint pairs at least " << min_length << " bonds apart\n";
  }

  int max_length = 7;
  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_length) || max_length < min_length) {
      cerr << "Max radius (-R) must be larger than min_length " << min_length << '\n';
      return 1;
    }

    if (verbose)
      cerr << "Will only fingerprint pairs separated by " << max_length << " bonds or less\n";
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! atom_typing.build(p)) {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 1;
    }
  } else {
    atom_typing.build("UST:Y");
    if (verbose)
      cerr << "Default compressed atomic number atom typing\n";
  }

  if (cl.option_present('f')) {
    function_as_tdt_filter = 1;

    if (verbose)
      cerr << "Will function as a TDT filter\n";
  }

  if (! cl.option_present('J')) {
    cerr << "Must specify tag via the -J option\n";
    usage(1);
  }

  if (cl.option_present('J')) {
    if (! fp_writer.Initialise(cl, 'J', verbose)) {
      cerr << "Cannot initialise fingerprint writer (-J)\n";
      return 1;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_tdt_filter) {
    if (cl.number_elements() > 1) {
      cerr << "When working as a filter, only one argument possible\n";
      usage(1);
    }
  } else if (! cl.option_present('i')) {
    if (! all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot discern input type(s)\n";
      return 3;
    }
  } else if (! process_input_type(cl, input_type)) {
    cerr << prog_name << ": cannot discern input type\n";
    usage(1);
  }

  if (cl.option_present('s')) {
    gather_molecule_statistics = true;
    if (verbose) {
      cerr << "Will gather statistics on molecules processed\n";
    }
  }

  LinearFingerprintGenerator linear_fp_gen;
  linear_fp_gen.set_min_length(min_length);
  linear_fp_gen.set_max_length(max_length);
  if (cl.option_present('y')) {
    linear_fp_gen.set_check_for_collisions(true);
    check_for_collisions = true;
    if (verbose)
      cerr << "WIll check for atom pair collisions\n";
  }

  if (cl.option_present('c')) {
    linear_fp_gen.set_check_coverage(true);
    if (verbose)
      cerr << "Will produce isotopically labelled molecules with coverage\n";
  }

  if (cl.option_present('x')) {
    linear_fp_gen.set_paths_can_cross(true);
    if (verbose)
      cerr << "Paths can cross\n";
  }

  if (cl.option_present('B')) {
    const char * fname = cl.option_value('B');
    if (!linear_fp_gen.OpenStreamForBitMeanings(fname)) {
      cerr << "Cannot open stream for labelled paths '" << fname << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Bit meanings written to '" << fname << "'\n";
  }

  if (cl.empty())
  {
    cerr << prog_name << ": no inputs\n";
    usage(1);
  }

  IWString_and_File_Descriptor output(1);

  fp_writer.WriteHeaderIfNeeded(output);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); ++i) {
    if (! LinearFingerprint(cl[i], input_type, linear_fp_gen, output)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose || check_for_collisions)
  {
    cerr << molecules_read << " molecules read\n";
    if (check_for_collisions)
      linear_fp_gen.ReportCollisions(cerr);
    cerr << "Molecules set btw " << bits_set.minval() << " and " << bits_set.maxval() << " bits, ave " << bits_set.average() << "\n";
  }

  output.flush();

  if (gather_molecule_statistics) {
    cerr << "Molecules had btw " << atom_count.minval() << " and " << atom_count.maxval() << " atoms, ave " << atom_count.average() << "\n";
    cerr << "Longest paths btw " << longest_path_acc.minval() << " and " << longest_path_acc.maxval() << " bonds, ave " << longest_path_acc.average() << "\n";
    for (int i = 0; i < longest_path.number_elements(); ++i) {
      if (longest_path[i])
        cerr << longest_path[i] << " molecules had a longest path of " << i << " bonds\n";
    }
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = LinearFingerprint(argc, argv);

  return rc;
}
