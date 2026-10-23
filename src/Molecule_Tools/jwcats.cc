/*
 * Compute the CATS descriptors according to the following paper
 * Angew. Chem. Int. Ed. 1999, 38, No. 19 by Gisbert Schneider, Werner Neidhart, Thomas
 * Giller, and Gerard Schmid
 */

#include <stdlib.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/allowed_elements.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/jwcats_lib.h"

namespace jwcats_driver {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off

  cerr << "Calculate CAT Search descriptors for molecules\n";

  cerr << "  -p             do not output values related to hydrophobic-hydrophobic pair\n";
  cerr << "  -m <n>         goes up to length of bonds(default 10)\n";
  cerr << "  -z <n>         minimum bond separation (features closer than <n> bonds are ignored\n";
  cerr << "  -q             use queries to define lipophilic agroups\n";
  cerr << "  -s             scaling type\n";
  cerr << "                 0 => no scaling\n";
  cerr << "                 1 => normalized by the number of heavy atom (default)\n";
  cerr << "                 2 => normalized by the sum of the number of two pharmacophore types\n";
  cerr << "                 3 => normalized by the atoms / (sum of the number of two pharmacophore types)\n";
  cerr << "  -J <tag>       create fingerprints\n";
  cerr << "  -H ...         donor acceptor assignment, enter '-H help' for info\n";
  cerr << "  -N ...         charge assigner, enter '-N help' for info\n";
  cerr << "  -h             make implicit hydrogens explicit\n";
  cerr << "  -f             function as a tdt filter\n";
  cerr << "  -o <sep>       output separator for descriptor output, recognised names like tab and comma\n";
  cerr << "  -X ...         more options\n";
  cerr << "  -i <type>      input type\n";
  cerr << "  -A ...         aromaticity options\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush      flush output after each molecule\n";

  ::exit(0);
}

void
ReportAccumulator(const Accumulator_Int<unsigned int>& acc, const char* atype,
                  std::ostream& os) {
  if (acc.n() == 0) {
    os << "No molecules had " << atype << " features\n";
    return;
  }

  os << "molecules had between " << acc.minval() << " and " << acc.maxval() << " "
     << atype << " ave " << static_cast<float>(acc.average_if_available_minval_if_not())
     << '\n';
}

class Options {
 private:
  Chemical_Standardisation _chemical_standardisation;
  Allowed_Elements _allowed_elements;

  int _verbose = 0;
  uint64_t _molecules_read = 0;
  int _reduce_to_largest_fragment = 0;
  uint64_t _molecules_containing_non_allowed_elements = 0;
  uint64_t _number_of_error = 0;

  IWString _fingerprint_tag;
  IWString _smiles_tag = "$SMI<";
  IWString _identifier_tag = "PCN<";

  int _function_as_gfp_filter = 0;
  int _flush_after_every_molecule = 0;
  char _output_separator = ' ';

  Accumulator_Int<unsigned int> _acceptor_acc;
  Accumulator_Int<unsigned int> _donor_acc;
  Accumulator_Int<unsigned int> _negative_acc;
  Accumulator_Int<unsigned int> _positive_acc;
  Accumulator_Int<unsigned int> _hydrophobe_acc;
  Accumulator_Int<unsigned int> _features_acc;

  Fraction_as_String _fraction_as_string;

  jwcats::JWCats _jwcats_calculator;

  int Preprocess(Molecule& m);
  int DoFingerprintOutput(Molecule& m, const std::vector<double>& scaled_counts,
                          IWString_and_File_Descriptor& output);
  int OutputResultHeader(IWString_and_File_Descriptor& output) const;
  int HandleMissingChargeDataFingerprint(Molecule& m, int array_size,
                                         IWString_and_File_Descriptor& output);
  int HandleMissingChargeData(Molecule& m, int array_size,
                              IWString_and_File_Descriptor& output);
  int Process(Molecule& m, IWString_and_File_Descriptor& output);
  int Process(data_source_and_type<Molecule>& input,
              IWString_and_File_Descriptor& output);
  int ProcessSmilesRecord(const const_IWSubstring& buffer,
                          IWString_and_File_Descriptor& output);
  int Process(iwstring_data_source& input, IWString_and_File_Descriptor& output);
  int Process(const char* fname, IWString_and_File_Descriptor& output);
  int Process(const char* fname, FileType input_type,
              IWString_and_File_Descriptor& output);
  int Report(std::ostream& output, std::ofstream& logfile, time_t current_time) const;

 public:
  int Main(int argc, char** argv);
};

int
Options::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  // Just assume this is getting rid of explicit hydrogens.
  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  } else {
    m.remove_all(1);
  }

  return 1;
}

int
Options::DoFingerprintOutput(Molecule& m, const std::vector<double>& scaled_counts,
                             IWString_and_File_Descriptor& output) {
  Sparse_Fingerprint_Creator sfc;

  const std::vector<int>& write_array_value = _jwcats_calculator.write_array_value();
  const int array_size = _jwcats_calculator.array_size();
  for (int i = 0; i < array_size; ++i) {
    if (scaled_counts[i] == 0.0) {
      continue;
    }

    if (!write_array_value[i]) {
      continue;
    }

    const int count = static_cast<int>(scaled_counts[i] + 0.01);
    sfc.hit_bit(i, count);
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tmp);

  output << _fingerprint_tag << tmp << ">\n";

  if (!_function_as_gfp_filter) {
    output << "|\n";
  }

  return 1;
}

int
Options::OutputResultHeader(IWString_and_File_Descriptor& output) const {
  for (const std::string& feature_name : _jwcats_calculator.FeatureNames()) {
    output << _output_separator << feature_name;
  }

  output << '\n';

  return 1;
}

int
Options::HandleMissingChargeDataFingerprint(Molecule& m, int array_size,
                                            IWString_and_File_Descriptor& output) {
  output << _fingerprint_tag << ">\n";

  ++_molecules_containing_non_allowed_elements;

  return 1;
}

int
Options::HandleMissingChargeData(Molecule& m, int array_size,
                                 IWString_and_File_Descriptor& output) {
  if (_fingerprint_tag.length() > 0) {
    return HandleMissingChargeDataFingerprint(m, array_size, output);
  }

  const_IWSubstring mname = m.name();
  mname.truncate_at_first(' ');

  output << mname;

  const std::vector<int>& write_array_value = _jwcats_calculator.write_array_value();
  for (int i = 0; i < array_size; ++i) {
    if (write_array_value[i]) {
      output << _output_separator << '.';
    }
  }

  output << '\n';

  return 1;
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  if (_function_as_gfp_filter) {
    ;
  } else if (_fingerprint_tag.length()) {  // write smiles before molecule gets changed
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  Preprocess(m);

  if (_allowed_elements.contains_non_allowed_atoms(m)) {
    ++_molecules_containing_non_allowed_elements;

    if (_fingerprint_tag.empty()) {  // descriptors, ignore molecule
      return 1;
    }

    output << _fingerprint_tag << ">\n";

    if (!_function_as_gfp_filter) {
      output << "|\n";
    }

    return 1;
  }

  jwcats::Result result;
  const jwcats::ComputeStatus status = _jwcats_calculator.Compute(m, result);
  if (status == jwcats::ComputeStatus::kMissingChargeData) {
    if (!_function_as_gfp_filter) {
      cerr << m.name() << " ERROR in calculation of partial charge, '" << m.name()
           << "'\n";
    }
    return HandleMissingChargeData(m, _jwcats_calculator.array_size(), output);
  }

  if (status != jwcats::ComputeStatus::kOk) {
    ++_number_of_error;
    return 0;
  }

  if (_verbose) {
    _acceptor_acc.extra(result.property_count[0]);
    _donor_acc.extra(result.property_count[1]);
    _positive_acc.extra(result.property_count[2]);
    _negative_acc.extra(result.property_count[3]);
    _hydrophobe_acc.extra(result.property_count[4]);

    unsigned int features_this_molecule = 0;
    for (int count : result.property_count) {
      features_this_molecule += count;
    }
    _features_acc.extra(features_this_molecule);
  }

  if (_fingerprint_tag.length() > 0) {
    DoFingerprintOutput(m, result.scaled_counts, output);

    if (_flush_after_every_molecule) {
      output.flush();
    }

    return 1;
  }

  append_first_token_of_name(m.name(), output);

  const std::vector<int>& write_array_value = _jwcats_calculator.write_array_value();
  for (int i = 0; i < _jwcats_calculator.array_size(); ++i) {
    if (!write_array_value[i]) {
      continue;
    }

    if (result.scaled_counts[i] == 0.0) {
      output << _output_separator << '0';
    } else {
      _fraction_as_string.append_number(output, result.scaled_counts[i]);
//    output << _output_separator << static_cast<float>(result.scaled_counts[i]);
    }
  }

  output << '\n';

  return 1;
}

int
Options::Process(data_source_and_type<Molecule>& input,
                 IWString_and_File_Descriptor& output) {
  Molecule* m;

  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (_verbose > 1) {
      cerr << _molecules_read << " processing '" << m->name() << "'\n";
    }

    if (!Process(*m, output)) {
      cerr << "ERROR in computing descriptors for " << m->name() << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
Options::ProcessSmilesRecord(const const_IWSubstring& buffer,
                             IWString_and_File_Descriptor& output) {
  const_IWSubstring mybuffer(buffer);
  mybuffer.remove_leading_chars(_smiles_tag.length());
  mybuffer.chop();

  Molecule m;

  if (!m.build_from_smiles(mybuffer)) {
    cerr << "Invalid smiles '" << mybuffer << "'\n";
    return 0;
  }

  return Process(m, output);
}

int
Options::Process(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(32768);

    if (!buffer.starts_with(_smiles_tag)) {
      continue;
    }

    if (!ProcessSmilesRecord(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Process(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, output);
}

int
Options::Process(const char* fname, FileType input_type,
                 IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose > 2) {
    input.set_verbose(1);
  }

  return Process(input, output);
}

int
Options::Report(std::ostream& output, std::ofstream& logfile, time_t current_time) const {
  time(&current_time);
  logfile << "\n\ncalculation ended at " << ctime(&current_time) << '\n';
  logfile << "Total Molecules read: " << _molecules_read
          << "\tTotal Error: " << _number_of_error << '\n';
  if (_molecules_read > 0) {
    logfile << "Error Rate: "
            << static_cast<double>(_number_of_error) /
                   static_cast<double>(_molecules_read)
            << '\n';
  }

  output << "Read " << _molecules_read << " molecules\n";

  output << "Processed " << _molecules_read << " molecules\n";
  ReportAccumulator(_acceptor_acc, "acceptor", output);
  ReportAccumulator(_donor_acc, "donor", output);
  ReportAccumulator(_positive_acc, "positive", output);
  ReportAccumulator(_negative_acc, "negative", output);
  ReportAccumulator(_hydrophobe_acc, "hydrophobe", output);
  ReportAccumulator(_features_acc, "features", output);

  if (_molecules_containing_non_allowed_elements) {
    output << _molecules_containing_non_allowed_elements
           << " molecules with non-allowed atom types\n";
  }

  return 1;
}

int
Options::Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:pi:H:q:N:F:m:s:g:J:K:flhz:X:o:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  _verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    Usage(2);
  }

  if (!process_standard_aromaticity_options(cl, _verbose)) {
    cerr << "Cannot process aromaticity options (-A)\n";
    Usage(5);
  }

  if (cl.option_count('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      Usage(4);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;

    if (_verbose) {
      cerr << "Will recuce to largest fragment\n";
    }
  }

  if (cl.option_present('p')) {
    _jwcats_calculator.SetIncludeHydrophobicPairs(0);

    if (_verbose) {
      cerr << "Will xxx hydrophobe-hydrophobe pairs\n";
    }
  }

  if (cl.option_present('s')) {
    int n = cl.option_count('s');

    int scaling_type = _jwcats_calculator.scaling_type();
    if (!cl.value('s', scaling_type, n - 1)) {
      cerr << "Invalid value for scaling type (-s) has to be a number\n";
      Usage(2);
    }

    if ((scaling_type < 0) || (scaling_type > 3)) {
      cerr << "The value for scaling type is not valid, should be 0, 1 or 2\n";
      cerr << "1 (normalized by number of heavy atoms) is used instead\n";
      scaling_type = 1;
    }
    _jwcats_calculator.SetScalingType(scaling_type);
  }

  if (cl.option_present('h')) {
    _jwcats_calculator.SetMakeImplicitHydrogensExplicit(1);

    if (_verbose) {
      cerr << "Will make implicit Hydrogens explicit\n";
    }
  }

  if (cl.option_present('z')) {
    int min_bond_separation = 0;
    if (!cl.value('z', min_bond_separation) || min_bond_separation < 1) {
      cerr << "Minimum bond separation (-z) must be a whole +ve number\n";
      Usage(2);
    }

    _jwcats_calculator.SetMinBondSeparation(min_bond_separation);

    if (_verbose) {
      cerr << "Features closer than " << min_bond_separation << " bonds apart ignored\n";
    }
  }

  if (cl.option_present('m')) {
    int max_bond_separation = _jwcats_calculator.max_bond_separation();
    if (!cl.value('m', max_bond_separation) || max_bond_separation < 1) {
      cerr << "Maximum bond separation (-m) must be a whole +ve number\n";
      Usage(2);
    }

    _jwcats_calculator.SetMaxBondSeparation(max_bond_separation);

    if (_verbose) {
      cerr << "Will perceive features up to " << max_bond_separation << " bonds apart\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _fingerprint_tag);

    if (_verbose) {
      cerr << "Will produce fingerprints with tag '" << _fingerprint_tag << "'\n";
    }

    if (!_fingerprint_tag.ends_with('<')) {
      _fingerprint_tag.add('<');
    }

    _jwcats_calculator.SetScalingType(0);

    if (cl.option_present('f')) {
      _function_as_gfp_filter = 1;

      if (_verbose) {
        cerr << "Will function as a filter\n";
      }
    }
  }

  if (cl.option_present('q')) {
    _jwcats_calculator.queries().resize(cl.option_count('q') + 100);
    if (!process_queries(cl, _jwcats_calculator.queries(), _verbose)) {
      cerr << "Cannot process queries from -q option(s)\n";
      return 6;
    }
    _jwcats_calculator.SetUseQueriesToDetermineHydrophobicity(1);
  }

  int nq = _jwcats_calculator.queries().number_elements();
  for (int i = 0; i < nq; ++i) {
    _jwcats_calculator.queries()[i]->set_find_unique_embeddings_only(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (_function_as_gfp_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('H')) {
    if (!_jwcats_calculator.donor_acceptor_assigner().construct_from_command_line(
            cl, 'H', _verbose)) {
      cerr << "Cannot initialise donor/acceptor assignment object\n";
      Usage(4);
    }
  }

  if (cl.option_present('N')) {
    if (!_jwcats_calculator.charge_assigner().construct_from_command_line(cl, 0, 'N')) {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      Usage(1);
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        _flush_after_every_molecule = 1;
        if (_verbose) {
          cerr << "Will flush after every molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  if (cl.option_present('o')) {
    IWString o;
    cl.value('o', o);
    if (!char_name_to_char(o)) {
      cerr << "Unrecognised output separator (-o) '" << o << "'\n";
      Usage(2);
    }
    _output_separator = o[0];

    if (_verbose) {
      cerr << "Descriptor output separator set to '" << o << "'\n";
    }
  }

  // After we have set output separator.
  constexpr int kOutputPrecision = 3;
  if (!_fraction_as_string.initialise(0.0f, 1.0f, kOutputPrecision)) {
    cerr << "Cannot initialise fractional descriptor output cache\n";
    return 0;
  }
  _fraction_as_string.set_leading_string(_output_separator);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  if (!_jwcats_calculator.Initialise()) {
    cerr << "Cannot initialise JWCats calculation\n";
    return 1;
  }

  set_default_iwstring_float_concatenation_precision(3);

  std::ofstream logfile;
  if (_verbose) {
    logfile.open("jwcatsearch_descriptor.log", std::ios::out);
  }

  if (_verbose && !logfile) {
    cerr << "jwcatsearch_descriptor.log file cannot be opened\n";
    return 0;
  }

  time_t current_time = 0;

  if (_verbose) {
    logfile << "This file collect all the info about error during the calculation\n";
    logfile << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
    time(&current_time);
    logfile << "calculation started at " << ctime(&current_time);
    logfile << "The complete Command was:" << '\n';
    for (int kk = 0; kk < argc; ++kk) {
      logfile << argv[kk] << " ";
    }
    logfile << '\n';
  }

  IWString_and_File_Descriptor output(1);

  if (_fingerprint_tag.empty()) {
    output << "Name";
    OutputResultHeader(output);
  }

#ifdef DEBUGGING_THING_TO_DETERMINE_PAIRS
  for (int i = 0; i < 5; i++) {
    for (int j = i; j < 5; j++) {
      int k = jwcats::PairNumber(i, j);
      cerr << " i = " << i << " j = " << j << " pair " << k << '\n';
    }
  }
#endif

  int rc = 0;

  if (_function_as_gfp_filter) {
    if (cl.number_elements() > 1) {
      cerr << "Fingerprint filter cannot have multiple arguments\n";
      exit(3);
    }

    rc = Process(cl[0], output) ? 0 : 1;
  } else {
    for (int i = 0; i < cl.number_elements(); ++i) {
      if (!Process(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (_verbose) {
    Report(cerr, logfile, current_time);
  }

  return rc;
}

}  // namespace jwcats_driver

int
main(int argc, char** argv) {
  jwcats_driver::Options options;
  return options.Main(argc, argv);
}
