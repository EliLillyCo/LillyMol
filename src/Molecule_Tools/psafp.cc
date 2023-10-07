// Conputes PSA fingerprints

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/nvrtspsa.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int function_as_tdt_filter = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fp_tag("NCPSA<");

static int bit_replicates = 9;

static int flush_after_each_molecule = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off

  cerr << "Produces Novartis Polar Surface Area fingerprint\n";
  cerr << "  -J <tag>      fingerprint tag to use\n";
  cerr << "  -p <n>        number of bit replicates in fingerprints\n";
  cerr << "  -f            function as tdt filter\n";
  cerr << "  -q            quiet mode, do not display unclassified atom messages\n";
  cerr << "  -z            assign zero contribution to unclassified atoms\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

/*
  Handles both the case of molecule on the command line and filter
*/

static int
psafp(Molecule& m, IWString_and_File_Descriptor& output) {
  molecules_read++;

  preprocess(m);

  if (!function_as_tdt_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  m.reduce_to_largest_fragment();  // after we have written the smiles

  double psa = novartis_polar_surface_area(m);

  int int_psa = static_cast<int>(psa / 10.0 + 0.49999);
  if (int_psa <= 0) {
    int_psa = 1;
  }

  Sparse_Fingerprint_Creator sfc;

  for (int i = 0; i < bit_replicates; i++) {
    sfc.hit_bit(i, int_psa);
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(fp_tag, tmp);

  output << tmp << '\n';

  if (!function_as_tdt_filter) {
    output << "|\n";
  }

  if (flush_after_each_molecule) {
    output.flush();
  }

  return output.good();
}

static int
psafp(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (!psafp(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
psafp(const char* fname, FileType input_type, IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return psafp(input, output);
}

static int
psafp_molecule(const_IWSubstring buffer,  // local copy
               IWString_and_File_Descriptor& output) {
  buffer.remove_leading_chars(smiles_tag.length());
  assert(buffer.ends_with('>'));

  buffer.chop(1);

  Molecule m;

  if (!m.build_from_smiles(buffer)) {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return psafp(m, output);
}

static int
psafp(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(32768);

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!psafp_molecule(buffer, output)) {
      return 0;
    }
  }

  return 1;
}

static int
psafp(const char* input_fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(input_fname);

  if (!input.good()) {
    cerr << "psafp:cannot open '" << input_fname << "'\n";
    return 0;
  }

  return psafp(input, output);
}

static void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush       flush output after every molecule\n";

  ::exit(0);
}

static int
psafp(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lfJ:p:zX:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!verbose) {
    set_display_psa_unclassified_atom_mesages(0);
  }

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('f')) {
    function_as_tdt_filter = 1;

    if (verbose) {
      cerr << "Will work as a TDT filter\n";
    }
  }

  if (cl.option_present('z')) {
    set_return_zero_for_unclassified_atoms(1);

    if (verbose) {
      cerr << "Unclassified atoms assigned zero value in psa computation\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', fp_tag);

    if (!fp_tag.ends_with('<')) {
      fp_tag << '<';
    }

    if (verbose) {
      cerr << "PSA written with tag '" << fp_tag << "'\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', bit_replicates) || bit_replicates <= 0) {
      cerr << "The bit replicates value (-p) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Fingerprint bit replicates set to " << bit_replicates << endl;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (function_as_tdt_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (function_as_tdt_filter) {
    if (!psafp(cl[0], output)) {
      rc = 3;
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!psafp(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = psafp(argc, argv);

  return rc;
}
