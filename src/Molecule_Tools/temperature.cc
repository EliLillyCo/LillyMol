// Creates a temperature fingerprint based on the atom count

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/mpr.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static IWString temperature_tag;
static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int molecules_with_too_many_atoms = 0;

static Accumulator_Int<int> atom_statistics;

static Accumulator_Int<int>* accumulators = nullptr;

static int allow_multi_fragment_molecules = 0;

static int write_array_of_properties = 0;

static int function_as_filter = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

// we can optionally create a Sparse Fingerprint with certain properties.
// Which properties are computed will be stored as bits in this mask.
static unsigned char abbreviated_fingerprints = 0;

// Bits that govern what properties are computed.

#define ABBREV_NATOMS 1
#define ABBREV_LARGEST_RING_SIZE 2
#define ABBREV_NRINGS 4
#define ABBREV_RING_ATOMS 8
#define ABBREV_AROAMTIC_ATOMS 16
#define ABBREV_FUSED_RING_ATOMS 32
#define ABBREV_HETEROATOM_COUNT 64
#define ABBREV_UNSATURATION_COUNT 128


static void
DisplayDashYOptions(std::ostream& output)
{
  output << " -Y flush     flush output after every molecule\n";
  ::exit(1);
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
  // clang-format off

  cerr << "Computes molecular properties as a fingerprint - the -MPR option in gfp_make\n";
  cerr << "  -J <dataitem>  specify dataitem tag for properties bits (default "
       << temperature_tag << ")\n";
  cerr << "  -B ...         abbreviated fingerprint specification, enter '-B help'\n";
  cerr << "  -d <num>       divisor for atom count related properties\n";
  cerr << "  -r <num>       bit replicates for abbreviated properties\n";
  cerr << "  -a             output is a descriptor file\n";
  cerr << "  -m             allow multi fragment molecules\n";
  cerr << "  -u             count aromatic atoms as unsaturated\n";
  cerr << "  -f             work as a filter\n";
  cerr << "  -O ...         complete control on what gets computed, enter '-O help' for info\n";
  cerr << "  -i <type>      specify input file type (except with -f)\n";
  (void)display_standard_aromaticity_options(cerr);
  (void)display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -E <X>         create element 'X', 'autocreate' for auto creation\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

struct Options {
 public:
  // If dealing with large molecules we can divide atom counts by a factor.
  int atom_count_divide = 1;
  int bit_replicates = 8;
  int flush_after_each_molecule = 0;

  // New style interface.
  mpr::MolecularPropertiesGenerator gen;

 public:
  int Initialise(Command_Line& cl);

  template <typename T> int SetupOutVector(std::unique_ptr<T[]>& properties) const;
};

template <typename T>
int
Options::SetupOutVector(std::unique_ptr<T[]>& properties) const {
  int rc;
  if (gen.number_features() > 0) {
    rc = gen.number_features();
  } else {
    rc = NPROPERTIES;
  }

  properties.reset(new T[rc]);

  return rc;
}

int
Options::Initialise(Command_Line& cl)
{
  const int verbose = cl.option_present('v');

  if (cl.option_present('r')) {
    if (!cl.value('r', bit_replicates) || bit_replicates < 1) {
      cerr << "The number of bit replicates must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will create " << bit_replicates << " bit replicates\n";
    }
  }

  if (cl.option_present('d')) {
    if (!cl.value('d', atom_count_divide) || atom_count_divide < 1) {
      cerr << "The atom count divisor (-d) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Atom count related property divisor set to " << atom_count_divide << '\n';
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('O')) {
    if (! gen.Initialise(cl, 'O')) {
      cerr << "Cannot initialise -O feature generation\n";
      return 0;
    }
  }

  return 1;
}

static int
write_descriptor_file_header(IWString_and_File_Descriptor& output)
{
  output << "ID natoms lgrsize nrings ring_atoms aroma fused_atoms htroa unsat\n";

  return 1;
}

static int
do_write_fingerprint(Molecule& m, const unsigned char* properties,
                     const int nproperties,
                     IWString_and_File_Descriptor& output)
{
  IW_Bits_Base fp;
  fp.construct_from_array_of_bits(properties, nproperties * IW_BITS_PER_BYTE);

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);
  output << temperature_tag << tmp << ">\n";

  return 1;
}

static void
set_bit_in_each_creator(Sparse_Fingerprint_Creator& sfc, int n, int b, int c)
{
  for (int i = 0; i < n; i++) {
    sfc.hit_bit(i * n + b, c);
  }

  return;
}

static int
do_write_abbreviated_properties(const Molecule& m, const unsigned char* properties,
                                const Options& options,
                                IWString_and_File_Descriptor& output)
{
  Sparse_Fingerprint_Creator sfc;

  // cerr << "To write " << static_cast<int>(abbreviated_fingerprints) << '\n';

  if (ABBREV_NATOMS & abbreviated_fingerprints &&
      properties[0] >= options.atom_count_divide) {
    int c = static_cast<int>(properties[0]) / options.atom_count_divide;

    set_bit_in_each_creator(sfc, options.bit_replicates, 0, c);
  }

  if (ABBREV_LARGEST_RING_SIZE & abbreviated_fingerprints && properties[1] > 0) {
    int c = static_cast<int>(properties[1]);
    set_bit_in_each_creator(sfc, options.bit_replicates, 1, c);
  }

  if (ABBREV_NRINGS & abbreviated_fingerprints && properties[2] > 0) {
    int c = static_cast<int>(properties[2]);
    set_bit_in_each_creator(sfc, options.bit_replicates, 2, c);
  }

  if (ABBREV_RING_ATOMS & abbreviated_fingerprints && properties[3] > 0) {
    int c = static_cast<int>(properties[3]);
    set_bit_in_each_creator(sfc, options.bit_replicates, 3, c);
  }

  if (ABBREV_AROAMTIC_ATOMS & abbreviated_fingerprints &&
      properties[4] >= options.atom_count_divide) {
    int c = static_cast<int>(properties[4]) / options.atom_count_divide;
    set_bit_in_each_creator(sfc, options.bit_replicates, 4, c);
  }

  if (ABBREV_FUSED_RING_ATOMS & abbreviated_fingerprints && properties[5] > 0) {
    int c = static_cast<int>(properties[5]);
    set_bit_in_each_creator(sfc, options.bit_replicates, 5, c);
  }

  if (ABBREV_HETEROATOM_COUNT & abbreviated_fingerprints &&
      properties[6] >= options.atom_count_divide) {
    int c = static_cast<int>(properties[6] / options.atom_count_divide);
    set_bit_in_each_creator(sfc, options.bit_replicates, 6, c);
  }

  if (ABBREV_UNSATURATION_COUNT & abbreviated_fingerprints && properties[7] > 0) {
    int c = static_cast<int>(properties[7]);  // don't bother dividing
    set_bit_in_each_creator(sfc, options.bit_replicates, 7, c);
  }

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(temperature_tag, tmp);

  output << tmp << "\n";

  return 1;
}

static int
do_write_array_of_properties(const Molecule& m, const unsigned char* properties,
                             const int nproperties,
                             IWString_and_File_Descriptor& output)
{
  append_first_token_of_name(m.name(), output);

  static constexpr char kSep = ' ';
  for (int i = 0; i < nproperties; i++) {
    output << kSep << static_cast<int>(properties[i]);
  }

  output << '\n';

  return 1;
}

static int
preprocess(Molecule& m)
{
  if (allow_multi_fragment_molecules) {
    ;
  } else if (reduce_to_largest_fragment && m.number_fragments() > 1) {
    m.reduce_to_largest_organic_fragment();
  } else if (m.number_fragments() > 1) {
    cerr << "Error, molecule has " << m.number_fragments() << " components\n";
    return 0;
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return 1;
}

/*
  May 2000. Found differences on different platforms because of differing
  SSSR determinations. Therefore change things to use all the rings - doesn't
  matter, all will be treated the same
*/

static Molecular_Properties_Generator mpg;

static int
temperature(Molecule& m, Options& options, IWString_and_File_Descriptor& output)
{
  static int nproperties = 0;
  static std::unique_ptr<uint8_t[]> properties;

  if (! properties) {
    nproperties = options.SetupOutVector(properties);
    if (verbose) {
      accumulators = new Accumulator_Int<int>[nproperties];
    }
  }
  
  const int matoms = m.natoms();

  atom_statistics.extra(matoms);

  if (matoms > std::numeric_limits<unsigned char>::max()) {
    molecules_with_too_many_atoms++;
    if (verbose) {
      cerr << "Molecule with too many atoms " << matoms << ' ' << m.name() << '\n';
    }
  }

  int rc;
  if (options.gen.number_features()) {
    rc = options.gen.GenerateMolecularProperties(m, properties.get());
  } else {
    rc = mpg(m, properties.get());
  }
  if (rc == 0) {
    cerr << "molecular_properties_generation failed '" << m.name() << '\n';
    std::fill_n(properties.get(), nproperties, 0);
    return 0;
  }

  if (nullptr != accumulators) {
    for (int i = 0; i < nproperties; i++) {
      accumulators[i].extra(int(properties[i]));
    }
  }

  if (write_array_of_properties) {
    return do_write_array_of_properties(m, properties.get(), nproperties, output);
  } else if (0 != abbreviated_fingerprints) {
    return do_write_abbreviated_properties(m, properties.get(), options, output);
  } else if (temperature_tag.length() > 0) {
    return do_write_fingerprint(m, properties.get(), nproperties, output);
  }

  cerr << "What kind of output am I supposed to be doing???\n";

  return 0;
}

static int
temperature(iwstring_data_source& input, Options& options,
            IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!buffer.starts_with(smiles_tag)) {
      output << buffer << '\n';

      output.write_if_buffer_holds_more_than(32768);

      continue;
    }

    const_IWSubstring smi = buffer.substr(5);
    smi.chop();

    Molecule m;
    if (!m.build_from_smiles(smi)) {
      cerr << "Yipes, cannot parse smiles, line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }

    preprocess(m);

    output << smiles_tag << m.smiles() << ">\n";

    if (!temperature(m, options, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
temperature(data_source_and_type<Molecule>& input, Options& options,
            IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (!preprocess(*m)) {
      return 0;
    }

    if (!write_array_of_properties) {
      output << smiles_tag << m->smiles() << ">\n";
      output << identifier_tag << m->name() << ">\n";
    }

    int rc = temperature(*m, options, output);

    if (!write_array_of_properties) {
      output << "|\n";
    }

    if (options.flush_after_each_molecule) {
      output.flush();
    } else {
      output.write_if_buffer_holds_more_than(32768);
    }

    if (0 == rc) {
      return 0;
    }
  }

  // GH: Fix bug# 19048
  if (input.stopped_because_of_error()) {
    // Return fail due to error
    return 0;
  } else {
    // Return success
    return 1;
  }
}

static int
temperature(const char* fname, FileType input_type, Options& options,
            IWString_and_File_Descriptor& output)
{
  if (function_as_filter) {
    iwstring_data_source input(fname);
    if (!input.ok()) {
      cerr << "Cannot open '" << fname << "' for input\n";
      return 0;
    }

    return temperature(input, options, output);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    // GH: Fix bug# 19048
    return 0;  // Original: return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return temperature(input, options, output);
}

static void
display_dash_B_qualifiers(std::ostream& os)
{
  os << " -B n     natoms\n";
  os << " -B l     largest ring size\n";
  os << " -B r     nrings\n";
  os << " -B g     ring atoms\n";
  os << " -B a     aromatic atoms\n";
  os << " -B f     fused ring atoms\n";
  os << " -B h     heteroatom count\n";
  os << " -B u     unsaturated atoms\n";

  exit(1);
}

static int
temperature(int argc, char** argv)
{
  // GH: Fix bug# 19048
  // Added for checking the error code
  int error_code = 1;
  Command_Line cl(argc, argv, "J:vA:mauE:sfi:lg:B:d:r:Y:O:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl, verbose)) {
    cerr << "Cannot parse -E option\n";
    usage(7);
  }

  if (!process_standard_aromaticity_options(cl)) {
    cerr << "Cannot parse aromaticity options\n";
    usage(8);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      usage(6);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will strip to largest fragments\n";
    }
  }

  if (cl.option_present('B')) {
    IWString b = cl.string_value('B');
    b.to_lowercase();

    if ("help" == b) {
      display_dash_B_qualifiers(cerr);
    }

    for (int i = 0; i < b.length(); i++) {
      char bi = b[i];

      if ('u' == bi) {
        abbreviated_fingerprints |= ABBREV_UNSATURATION_COUNT;
      } else if ('h' == bi) {
        abbreviated_fingerprints |= ABBREV_HETEROATOM_COUNT;
      } else if ('f' == bi) {
        abbreviated_fingerprints |= ABBREV_FUSED_RING_ATOMS;
      } else if ('n' == bi) {
        abbreviated_fingerprints |= ABBREV_NATOMS;
      } else if ('r' == bi) {
        abbreviated_fingerprints |= ABBREV_NRINGS;
      } else if ('g' == bi) {
        abbreviated_fingerprints |= ABBREV_RING_ATOMS;
      } else if ('l' == bi) {
        abbreviated_fingerprints |= ABBREV_LARGEST_RING_SIZE;
      } else if ('a' == bi) {
        abbreviated_fingerprints |= ABBREV_AROAMTIC_ATOMS;
      } else {
        cerr << "Unrecognised -B qualifier '" << b << "'\n";
        cerr << "                           ";
        for (int j = 0; j < i; j++) {
          cerr << ' ';
        }
        cerr << "^\n";
        display_dash_B_qualifiers(cerr);
      }
    }

#ifdef NO_LONGER_NEEDED_QWEQWE
    if (cl.option_present('d')) {
      if (!cl.value('d', atom_count_divide) || atom_count_divide < 1) {
        cerr << "The atom count divisor (-d) must be a whole +ve number\n";
        usage(2);
      }

      if (verbose) {
        cerr << "Atom count related property divisor set to " << atom_count_divide
             << '\n';
      }
    }

    if (cl.option_present('r')) {
      if (!cl.value('r', bit_replicates) || bit_replicates < 1) {
        cerr << "The number of bit replicates must be a whole +ve number\n";
        usage(2);
      }

      if (verbose) {
        cerr << "Will create " << bit_replicates << " bit replicates\n";
      }
    }
#endif
  }

  if (cl.option_present('J')) {
    cl.value('J', temperature_tag);

    if (verbose) {
      cerr << "Temperature fingerprints stored as '" << temperature_tag << "' dataitem\n";
    }

    if (!temperature_tag.ends_with('<')) {
      temperature_tag << '<';
    }
  }

  if (cl.option_present('m')) {
    allow_multi_fragment_molecules = 1;
    if (verbose) {
      cerr << "Molecules with multiple fragments can be processed\n";
    }
  }

  if (!cl.option_present('a') && !cl.option_present('J')) {
    cerr << "To get output must specify one of -a or -J\n";
    usage(4);
  }

  if (cl.option_present('a') && cl.option_present('J')) {
    cerr << "Only one of -a and -J is allowed\n";
    usage(4);
  }

  if (cl.option_present('a')) {
    write_array_of_properties = 1;
    if (verbose) {
      cerr << "Molecular properties written as descriptor file\n";
    }
  }

  if (cl.option_present('i') && cl.option_present('f')) {
    cerr << "The -i (input type) and -f (filter) options are mutually exclusive\n";
    usage(17);
  }

  if (cl.option_present('f')) {
    function_as_filter = 1;
    if (verbose) {
      cerr << "Will work as a filter\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('f')) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('u')) {
    mpg.set_unsaturation_includes_aromatic(1);
    if (verbose) {
      cerr << "Unsaturated aromatic atoms included as unsaturated\n";
    }
  }

  if (cl.option_present('s')) {
    accumulators = new Accumulator_Int<int>[NPROPERTIES];
    if (verbose) {
      cerr << "Will gather statistics on themolecules\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(99);
  }

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  if (write_array_of_properties) {
    write_descriptor_file_header(output);
  }

  for (int i = 0; i < cl.number_elements(); i++) {
    // GH: Fix bug# 19048
    // The return value for temperature function will be checked
    if (!temperature(cl[i], input_type, options, output)) {
      error_code = 0;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Produced molecular_property fingerprints for " << atom_statistics.n()
         << " molecules\n";
    if (molecules_with_too_many_atoms) {
      cerr << molecules_with_too_many_atoms << " molecules had more than "
           << std::numeric_limits<unsigned char>::max() << " atoms\n";
    }

    cerr << "Molecules between " << atom_statistics.minval() << " and "
         << atom_statistics.maxval() << " atoms\n";
    if (atom_statistics.n() > 1) {
      cerr << "ave " << atom_statistics.average() << " variance "
           << atom_statistics.variance();
    }
    cerr << '\n';
  }

  if (nullptr != accumulators) {
    for (int i = 0; i < NPROPERTIES; i++) {
      Accumulator_Int<int>& acc = accumulators[i];
      cerr << "Property " << i << " between " << acc.minval() << " and " << acc.maxval();
      if (acc.n() > 1) {
        cerr << " ave " << acc.average() << " variance " << acc.variance();
      }
      cerr << '\n';
    }

    delete[] accumulators;
  }

  // GH: Fix bug# 19048
  // Return the error code
  return error_code;  // original:  return 0
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = temperature(argc, argv);
  // GH: Fix bug# 19048
  // Flip the error return based on Linux convention: 0 is success, otherwise fails
  rc = (0 == rc) ? 1 : 0;
  return rc;
}
