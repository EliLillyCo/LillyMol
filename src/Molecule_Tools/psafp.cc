// Computes PSA fingerprints

#include <iostream>
#include <limits>
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

const char* prog_name = nullptr;

static int verbose = 0;

static uint64_t molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int function_as_tdt_filter = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fp_tag("NCPSA<");

static int bit_replicates = 9;

extending_resizable_array<int> acc_psa;

static int flush_after_each_molecule = 0;

static int generate_descriptor_file = 0;
static float min_psa = 0.0f;
static float max_psa = std::numeric_limits<float>::max();
static uint64_t rejected_for_out_of_range = 0;
static int min_or_max_specified = 0;
static char output_separator = ' ';

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

  cerr << R"(Produces Novartis Polar Surface Area fingerprint\n";
  -J <tag>      fingerprint tag to use, default 'NCPSA'
  -p <n>        number of bit replicates in fingerprints.
  -f            function as tdt filter.
  -q            quiet mode, do not display unclassified atom messages.
  -z            assign zero contribution to unclassified atoms.
  -T ...        tabular descriptor output, enter '-T help' for info. Can do filtering.
  -X ...        miscellaneous other options, enter '-X help' for info.
  -l            reduce to largest fragment.
  -i <type>     input specification.
  -g ...        chemical standardisation options.
  -E ...        standard element specifications.
  -A ...        standard aromaticity specifications.
  -v            verbose output.
)";
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

static int 
GenerateDescriptors(Molecule& m, float psa,
        IWString_and_File_Descriptor& output) {
  if (! min_or_max_specified) {
  } else if (psa < min_psa) {
    ++rejected_for_out_of_range;
    return 1;
  } else if (psa > max_psa) {
    ++rejected_for_out_of_range;
    return 1;
  }

  output << m.name() << output_separator << psa << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
psafp(Molecule& m, IWString_and_File_Descriptor& output) {
  molecules_read++;

  preprocess(m);

  if (generate_descriptor_file) {
    // do nothing yet.
  } else if (!function_as_tdt_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  m.reduce_to_largest_fragment();  // after we have written the smiles

  double psa = novartis_polar_surface_area(m);

  if (verbose) {
    ++acc_psa[static_cast<int>(psa + 0.49999)];
  }

  int int_psa = static_cast<int>(psa / 10.0 + 0.49999);
  if (int_psa <= 0) {
    int_psa = 1;
  }

  if (generate_descriptor_file) {
    return GenerateDescriptors(m, psa, output);
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

static void
DisplayDashTOptions(std::ostream& output) {
  output << R"(Options for generating a descriptor file.
-T min=<n>      minimum value for TPSA
-T max=<n>      maximum value for TPSA
-T sep=<char>   output separator, default ' '.
-T .            default, no filtering, all molecules written.
)";

  ::exit(0);
}

static int
psafp(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lfJ:p:zX:T:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!verbose) {
    nvrtspsa::set_display_psa_unclassified_atom_mesages(0);
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
    nvrtspsa::set_return_zero_for_unclassified_atoms(1);

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
      cerr << "Fingerprint bit replicates set to " << bit_replicates << '\n';
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

  if (cl.option_present('T')) {
    IWString t;
    for (int i = 0; cl.value('T', t, i); ++i) {
      if (t == '.') {
      } else if (t.starts_with("min=")) {
        t.remove_leading_chars(4);
        if (! t.numeric_value(min_psa)) {
          cerr << "Invalid min psa specification '" << t << "'\n";
          return 1;
        }
        if (verbose) {
          cerr << "Will not write molecules with TPSA values below " << min_psa << '\n';
        }
        min_or_max_specified = 1;
      } else if (t.starts_with("max=")) {
        t.remove_leading_chars(4);
        if (! t.numeric_value(max_psa)) {
          cerr << "Invalid max psa specification '" << t << "'\n";
          return 1;
        }
        if (verbose) {
          cerr << "Will not write molecules with TPSA values above " << max_psa << '\n';
        }
        min_or_max_specified = 1;
      } else if (t.starts_with("sep=")) {
        t.remove_leading_chars(4);
        if (! char_name_to_char(t)) {
          cerr << "Invalid descriptor file output separator '" << t << "'\n";
          return 0;
        }
        output_separator = t[0];
        if (verbose) {
          cerr << "Descriptor file output separator " << output_separator << '\n';
        }
      } else if (t == "help") {
        DisplayDashTOptions(cerr);
      } else {
        cerr << "Unrecognised -T qualifier '" << t << "'\n";
        DisplayDashTOptions(cerr);
      }
    }

    if (min_psa > max_psa) {
      cerr << "Invalid min " << min_psa << " max " << max_psa << " TPSA specification\n";
      return 1;
    }

    generate_descriptor_file = 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  if (generate_descriptor_file) {
    output << "ID" << output_separator << "nvrtspsa\n";
  }

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
    uint64_t n = 0;
    uint64_t sum = 0;
    int minv = -1;
    int maxv = 0;
    for (int i = 0; i < acc_psa.number_elements(); ++i) {
      if (acc_psa[i] == 0) {
        continue;
      }

      cerr << acc_psa[i] << " molecules with value " << i << '\n';

      if (minv < 0) {
        minv = i;
      }
      if (i > maxv) {
        maxv = i;
      }

      if (acc_psa[i] > 0) {
        n += acc_psa[i];
        sum += i * acc_psa[i];
      }
    }

    float ave = iwmisc::Fraction<double>(sum, n);
    cerr << "Values btw " << minv << " and " << maxv << " ave " << ave << '\n';

    if (generate_descriptor_file && min_or_max_specified) {
      cerr << rejected_for_out_of_range << " molecules rejected for out of range " <<
              min_psa << ' ' << max_psa << '\n';
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = psafp(argc, argv);

  return rc;
}
