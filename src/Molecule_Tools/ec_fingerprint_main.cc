// Compute ec fingerprints.

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/ec_fingerprint.h"
#include "Molecule_Tools/fingerprint_writer.h"

using std::cerr;

const char* prog_name = nullptr;

namespace ec_fingerprint {
int verbose = 0;

int molecules_read = 0;

Chemical_Standardisation chemical_standardisation;

Atom_Typing_Specification atom_typing;

int reduce_to_largest_fragment = 0;

JobParameters job_parameters;

const IWString smiles_tag("$SMI<");
const IWString identifier_tag("PCN<");

// Potentially interesting statistics on molecules processed.

bool gather_molecule_statistics = false;
Accumulator_Int<uint> atom_count;
Accumulator_Int<uint> longest_path_acc;
extending_resizable_array<int> longest_path;

void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Computes EC circular fingerprints\n";
  cerr << "  -r <rad>      minimum shell radius (def 0)\n";
  cerr << "  -R <rad>      maximum shell radius (def 3)\n";
  cerr << "  -P ...        atom type specification\n";
  cerr << "  -J <tag>      tag for fingerprints\n";
  cerr << "  -f            function as a TDT filter\n";
  cerr << "  -D ...        specify what operation is to be performed\n";
  cerr << "  -D file=<fname>  For those actions that require a file.\n";
  cerr << "  -D fp=<Tag>      Generate fingerprints with a given tag\n";
  cerr << "  -D explain       Look for bits in <fname> and provide explanations\n";
  cerr << "  -D collision     Look for bit collisions and write\n";
  cerr << "  -D pgen          Gather precedent data and write\n";
  cerr << "  -D puse          Read previously generated precedent data and use\n";
  cerr << "  -D writeall      Write all bits generated to <fname>\n";
  cerr << "  -D filter        Write smiles of molecules that do NOT contain bits in <fname>\n";
  cerr << "  -D args=...      Pass arbitrary text directives to the class doing the work\n";
  cerr << "  -n            suppress fingerprint output (useful with the some -D options)\n";
  cerr << "  -s            gather statistics on molecules processed\n";
  cerr << "  -p            precise fingerprints - includes attachment point in inner shell\n";
  cerr << "  -m            bit formation is multiplicative. shells not differentiated\n";
  cerr << "  -q <query>    query specification of atoms to be processed\n";
  cerr << "  -w <nbits>    generate fixed width binary fingerprints, `nbits` bits\n";
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
Preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

void
GatherStatistics(Molecule& m) {
  atom_count.extra(m.natoms());

  const int lp = m.longest_path();

  longest_path_acc.extra(lp);
  longest_path[lp]++;
}

template <typename T>
int
EcFingerprint(Molecule& m, const atom_type_t* atype, ECFingerprint& ec_fp_gen, T& op,
              IWString_and_File_Descriptor& output) {
  if (gather_molecule_statistics) {
    GatherStatistics(m);
  }

  // A single atom molecule would produce zero bits, or two atoms in different
  // fragments...
  if (!ec_fp_gen.Fingerprint(m, nullptr, atype, op)) {
    cerr << "EcFingerprint:fingerprinting failed " << m.smiles() << " ignored\n";
    return 1;
  }

  op.FingerprintingComplete(m);

  if (!job_parameters.produce_output) {
    return 1;
  }

  op.DoAnyOutput(m, job_parameters, output);

  if (job_parameters.fp_writer.IsWritingDescriptors()) {
  } else if (job_parameters.function_as_tdt_filter) {
  } else {
    output << "|\n";
  }

  return 1;
}

template <typename T>
int
EcFingerprint(Molecule& m, ECFingerprint& ec_fp_gen, T& op,
              IWString_and_File_Descriptor& output) {
  atom_type_t* atype = new atom_type_t[m.natoms()];
  std::unique_ptr<atom_type_t[]> free_atype(atype);

  m.compute_aromaticity();

  if (!atom_typing.assign_atom_types(m, atype)) {
    cerr << "EcFingerprint::cannot assign atom types " << m.smiles() << ' ' << m.name()
         << '\n';
    return 0;
  }

  return EcFingerprint(m, atype, ec_fp_gen, op, output);
}

template <typename T>
int
EcFingerprint(data_source_and_type<Molecule>& input, ECFingerprint& ec_fp_gen, T& op,
              IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    Preprocess(*m);

    if (!EcFingerprint(*m, ec_fp_gen, op, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

template <typename T>
int
EcFingerprintPipe(const_IWSubstring line,  // Note local copy
                  ECFingerprint& ec_fp_gen, T& op, IWString_and_File_Descriptor& output) {
  assert(line.ends_with('>'));
  assert(line.starts_with(smiles_tag));

  line.remove_leading_chars(smiles_tag.length());
  line.chop();

  Molecule m;

  if (!m.build_from_smiles(line)) {
    cerr << "EcFingerprintPipe:invalid smiles '" << line << "'\n";
    return 0;
  }

  Preprocess(m);

  return EcFingerprint(m, ec_fp_gen, op, output);
}

template <typename T>
int
EcFingerprintPipe(iwstring_data_source& input, ECFingerprint& ec_fp_gen, T& op,
                  IWString_and_File_Descriptor& output) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    output << line << '\n';

    if (!line.starts_with(smiles_tag)) {
      continue;
    }

    if (!EcFingerprintPipe(line, ec_fp_gen, op, output)) {
      cerr << "EcFingerprintPipe:invalid input " << line << "' line "
           << input.lines_read() << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

template <typename T>
int
EcFingerprint(const char* fname, FileType input_type, ECFingerprint& ec_fp_gen, T& op,
              IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  if (job_parameters.function_as_tdt_filter) {
    iwstring_data_source input(fname);
    if (!input.good()) {
      cerr << "EcFingerprint::cannot open filter " << fname << '\n';
      return 0;
    }

    return EcFingerprintPipe(input, ec_fp_gen, op, output);
  }

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return EcFingerprint(input, ec_fp_gen, op, output);
}

template <typename T>
int
DoEcFingerprint(const Command_Line& cl, const FileType input_type,
                ECFingerprint& ec_fp_gen, T& op, IWString_and_File_Descriptor& output) {
  for (int i = 0; i < cl.number_elements(); ++i) {
    if (!EcFingerprint(cl[i], input_type, ec_fp_gen, op, output)) {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return 0;
    }
  }

  return output.good();
}

std::tuple<IWString, IWString>
DirectiveAndMaybeFile(Command_Line& cl, char flag) {
  std::tuple<IWString, IWString> to_be_returned;

  const_IWSubstring d;
  for (int i = 0; cl.value(flag, d, i); ++i) {
    if (d.starts_with("file=")) {
      d.remove_leading_chars(5);
      std::get<1>(to_be_returned) = d;
    } else {
      std::get<0>(to_be_returned) = d;
    }
  }

  return to_be_returned;
}

// Examine the -flag options on the command line, and
// identify two special directives, file= and args=
// Return a tuple of
//  0: directive
//  1: any file= value
//  2: any args= value
std::tuple<IWString, IWString, IWString>
DirectiveAnd(Command_Line& cl, char flag) {
  std::tuple<IWString, IWString, IWString> to_be_returned;

  const_IWSubstring d;
  for (int i = 0; cl.value(flag, d, i); ++i) {
    if (d.starts_with("file=")) {
      d.remove_leading_chars(5);
      std::get<1>(to_be_returned) = d;
    } else if (d.starts_with("args=")) {
      d.remove_leading_chars(5);
      std::get<2>(to_be_returned) = d;
    } else {
      std::get<0>(to_be_returned) = d;
    }
  }

  return to_be_returned;
}

int
EcFingerprint(int argc, char** argv) {
  Command_Line cl(argc, argv, "A:E:K:lg:i:J:P:vfr:R:ysB:cmpM:nD:q:w:");

  if (cl.unrecognised_options_encountered()) {
    usage(1);
  }

  verbose = cl.option_count('v');

  (void)process_elements(cl);

  if (!process_standard_smiles_options(cl, verbose)) {
    usage(7);
  }

  set_global_aromaticity_type(Daylight);

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(8);
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 1;
    }
  }

  job_parameters.smiles_tag = smiles_tag;
  job_parameters.identifier_tag = identifier_tag;

  int min_radius = 0;
  if (cl.option_present('r')) {
    if (!cl.value('r', min_radius) || min_radius < 0) {
      cerr << "The min atom pair radius (-r) must be >= 0\n";
      return 1;
    }

    if (verbose) {
      cerr << "Will only fingerprint pairs at least " << min_radius << " bonds apart\n";
    }
  }

  int max_radius = 3;
  if (cl.option_present('R')) {
    if (!cl.value('R', max_radius) || max_radius < min_radius) {
      cerr << "Max radius (-R) must be larger than min_radius " << min_radius << '\n';
      return 1;
    }

    if (verbose) {
      cerr << "Will only fingerprints to radius " << max_radius << '\n';
    }
  }

  if (cl.option_present('P')) {
    cl.value('P', job_parameters.atom_type_string);
    if (!atom_typing.build(job_parameters.atom_type_string)) {
      cerr << "Invalid atom typing specification '" << job_parameters.atom_type_string
           << "'\n";
      return 1;
    }
  } else {
    job_parameters.atom_type_string = "UST:Y";
    atom_typing.build(job_parameters.atom_type_string);
    if (verbose) {
      cerr << "Default compressed atomic number atom typing\n";
    }
  }

  if (cl.option_present('f')) {
    job_parameters.function_as_tdt_filter = 1;

    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }
  }

  if (cl.option_present('w')) {
    if (!cl.value('w', job_parameters.write_fixed_width_fingerprint) ||
        job_parameters.write_fixed_width_fingerprint < 8) {
      cerr << "Invalid fixed width fingerprint directive (-w)\n";
      return 1;
    }
    if (verbose) {
      cerr << "Will generate fixed width fingerprints "
           << job_parameters.write_fixed_width_fingerprint << " bits\n";
    }
  }

  // if (! cl.option_present('D'))
  //{
  //   cerr << "Must specify what operation to perform via the -D option\n";
  //   usage(1);
  // }

  FileType input_type = FILE_TYPE_INVALID;
  if (job_parameters.function_as_tdt_filter) {
    if (cl.number_elements() > 1) {
      cerr << "When working as a filter, only one argument possible\n";
      usage(1);
    }
  } else if (!cl.option_present('i')) {
    if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot discern input type(s)\n";
      return 3;
    }
  } else if (!process_input_type(cl, input_type)) {
    cerr << prog_name << ": cannot discern input type\n";
    usage(1);
  }

  if (cl.option_present('s')) {
    gather_molecule_statistics = true;
    if (verbose) {
      cerr << "Will gather statistics on molecules processed\n";
    }
  }

  if (cl.option_present('n')) {
    job_parameters.produce_output = false;

    if (verbose) {
      cerr << "Will produce no output\n";
    }
  }

  if (job_parameters.produce_output && !cl.option_present('J')) {
    cerr << "Must specify tag via the -J option\n";
    usage(1);
  }

  if (cl.option_present('J')) {
    if (!job_parameters.fp_writer.Initialise(cl, 'J', verbose)) {
      cerr << "Cannot initiwlise fingerprint_writer\n";
      return 1;
    }
  }

  ECFingerprint ec_fp_gen;
  ec_fp_gen.set_min_radius(min_radius);
  ec_fp_gen.set_max_radius(max_radius);

  if (cl.option_present('p')) {
    ec_fp_gen.set_precise(true);
    if (verbose) {
      cerr << "Will generate precise fingerprints\n";
    }
  }

  if (cl.option_present('m')) {
    ec_fp_gen.set_additive(false);

    if (verbose) {
      cerr << "Will use multiplicative bit formation\n";
    }
  }

  if (cl.empty()) {
    cerr << prog_name << ": no inputs\n";
    usage(1);
  }

  IWString_and_File_Descriptor output(1);

  job_parameters.fp_writer.WriteHeaderIfNeeded(output);

  auto [directive, fname, args] = DirectiveAnd(cl, 'D');

  if (!cl.option_present('D') && cl.option_present('J')) {
    directive = "fp";
  }

  if (directive.empty()) {
    cerr << "The -D option must contain a valid action\n";
    usage(1);
  }

  int rc = 0;

  if ("fp" == directive) {
    ProduceFingerprint produce_fingerprint;
    rc = DoEcFingerprint(cl, input_type, ec_fp_gen, produce_fingerprint, output);
  } else if (directive.starts_with("cov")) {
    AtomMapCoverage coverage;
    rc = DoEcFingerprint(cl, input_type, ec_fp_gen, coverage, output);
  } else  // All others require a file
  {
    if (fname.empty()) {
      cerr << "Requested profile generation but no file= specified\n";
      usage(1);
    }

    if ("pgen" == directive) {
      ECBuildPrecedent precedent;
      rc = DoEcFingerprint(cl, input_type, ec_fp_gen, precedent, output);
      if (rc) {
        rc = precedent.WritePrecedentData(' ', job_parameters, fname);
        precedent.Report(std::cerr);
      }
    } else if ("puse" == directive) {
      ECUsePrecedent precedent(max_radius);
      if (!precedent.ReadPrecedentData(fname)) {
        cerr << "Cannot read precedent data from '" << fname << "'\n";
        return 1;
      }
      rc = DoEcFingerprint(cl, input_type, ec_fp_gen, precedent, output);
      precedent.Report(cerr);
    } else if (directive.starts_with("coll")) {
      ECCheckCollisions collisions;
      if (!collisions.Open(fname)) {
        return 1;
      }
      rc = DoEcFingerprint(cl, input_type, ec_fp_gen, collisions, output);
    } else if ("writeall" == directive) {
      WriteAllBits writebits;
      if (!writebits.Open(fname)) {
        return 1;
      }
      rc = DoEcFingerprint(cl, input_type, ec_fp_gen, writebits, output);
    } else if ("explain" == directive) {
      ECBitMeanings bit_meanings;
      if (!bit_meanings.ReadBitsToFind(fname)) {
        return 1;
      }
      rc = DoEcFingerprint(cl, input_type, ec_fp_gen, bit_meanings, output);
    } else if (directive == "filter") {
      ECFilterByBits filter_by_bits;
      if (!filter_by_bits.ReadBitsToFilter(fname)) {
        cerr << "Cannot read bits to filter from '" << fname << "'\n";
        return 1;
      }
      rc = DoEcFingerprint(cl, input_type, ec_fp_gen, filter_by_bits, output);
    } else {
      cerr << "Unrecognised -D qualifier '" << directive << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << molecules_read << " molecules read\n";
  }

  if (gather_molecule_statistics) {
    cerr << "Molecules had btw " << atom_count.minval() << " and " << atom_count.maxval()
         << " atoms, ave " << atom_count.average() << "\n";
    cerr << "Longest paths btw " << longest_path_acc.minval() << " and "
         << longest_path_acc.maxval() << " bonds, ave " << longest_path_acc.average()
         << "\n";
    for (int i = 0; i < longest_path.number_elements(); ++i) {
      if (longest_path[i]) {
        cerr << longest_path[i] << " molecules had a longest path of " << i << " bonds\n";
      }
    }
  }

  return !rc;
}

}  // namespace ec_fingerprint

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = ec_fingerprint::EcFingerprint(argc, argv);

  return rc;
}
