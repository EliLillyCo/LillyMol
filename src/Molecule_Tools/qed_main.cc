// Command line interface to QED calculations.

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/qed.h"

namespace qed {

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
  cerr << R"(Computes QED
  -a            write all components of the QED score as descriptors
  -o <sep>      output separator (default ' ')
  -f            ignore molecules that fail the calculation - default is to stop processing.
  -M ...        miscellaneous options, enter '-M help' for info
  -v            verbose output\n";
)";

  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    Chemical_Standardisation _chemical_standardisation;

    int _write_descriptors;

    int _ignore_failed_calculations;

    // Sometimes it is useful to be able to prepend a suffix to the
    // default descriptor names.
    IWString _descriptor_prefix;

    // by default, we write space separated ouput.
    IWString _output_separator;

    qed::Qed _qed;

    Accumulator<double> _acc_amw;
    Accumulator<double> _acc_alogp;
    extending_resizable_array<uint32_t> _acc_hba;
    extending_resizable_array<uint32_t> _acc_hbd;
    Accumulator<double> _acc_psa;
    extending_resizable_array<uint32_t> _acc_rotb;
    extending_resizable_array<uint32_t> _acc_arom;
    extending_resizable_array<uint32_t> _acc_alerts;

    uint32_t _failed_calculations;

    int _molecules_read = 0;

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    int Preprocess(Molecule& m);

    int WriteHeader(IWString_and_File_Descriptor& output) const;

    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _molecules_read = 0;

  _write_descriptors = 0;
  _failed_calculations = 0;

  _output_separator = " ";
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('a')) {
    _write_descriptors = 1;
    if (_verbose) {
      cerr << "Will write all QED components\n";
    }
  }

  if (cl.option_present('f')) {
    _ignore_failed_calculations = 1;

    if (_verbose) {
      cerr << "Will ignore failed calculations\n";
    }
  }

  if (cl.option_present('o')) {
    cl.value('o', _output_separator);
    if (! char_name_to_char(_output_separator)) {
      cerr << "Invalid output separator specification (-o) '" << _output_separator << "'\n";
      return 0;
    }
  }

  if (cl.option_present('M')) {
    IWString m;
    for (int i = 0; cl.value('M', m, i); ++i) {
      if (m.starts_with("PREFIX=")) {
        m.remove_leading_chars(7);
        _descriptor_prefix = m;
      } else if (m == "help") {
      } else {
      }
    }
  }

  if (! cl.option_present('Q')) {
    cerr << "Must specify the QED alerts file via the -Q option\n";
    return 0;
  }

  if (! _qed.Initialise(cl, 'Q')) {
    cerr << "Cannot initialise QED configuration (-Q)\n";
    return 0;
  }

  return 1;
}

int
Options::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "Id";
  if (_write_descriptors) {
    output << _output_separator << _descriptor_prefix << "amw";
    output << _output_separator << _descriptor_prefix << "alogp";
    output << _output_separator << _descriptor_prefix << "hba";
    output << _output_separator << _descriptor_prefix << "hbd";
    output << _output_separator << _descriptor_prefix << "psa";
    output << _output_separator << _descriptor_prefix << "rotb";
    output << _output_separator << _descriptor_prefix << "arom";
    output << _output_separator << _descriptor_prefix << "alerts";
  }

  output << _output_separator << "QED\n";

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _failed_calculations << " failed\n";

  if (_molecules_read == 0) {
    return 1;
  }

  if (_acc_amw.empty()) {
    return 1;
  }

  output << "AMW  btw " << _acc_amw.minval() << " and " << _acc_amw.maxval() << 
            " ave " << static_cast<float>(_acc_amw.average()) << '\n';
  output << "ALOGP btw " << _acc_alogp.minval() << " and " << _acc_alogp.maxval() << 
            " ave " << static_cast<float>(_acc_alogp.average()) << '\n';
  output << "PSA   btw " << _acc_psa.minval() << " and " << _acc_psa.maxval() << 
            " ave " << static_cast<float>(_acc_psa.average()) << '\n';

  for (int i = 0; i < _acc_hba.number_elements(); ++i) {
    if (_acc_hba[i] == 0) {
      continue;
    }
    output << _acc_hba[i] << " molecules had " << i << " Hydrogen Bond Acceptors\n";
  }

  for (int i = 0; i < _acc_hbd.number_elements(); ++i) {
    if (_acc_hbd[i] == 0) {
      continue;
    }
    output << _acc_hbd[i] << " molecules had " << i << " Hydrogen Bond Donors\n";
  }

  for (int i = 0; i < _acc_rotb.number_elements(); ++i) {
    if (_acc_rotb[i] == 0) {
      continue;
    }
    output << _acc_rotb[i] << " molecules had " << i << " Rotatable Bonds\n";
  }

  for (int i = 0; i < _acc_arom.number_elements(); ++i) {
    if (_acc_arom[i] == 0) {
      continue;
    }
    output << _acc_arom[i] << " molecules had " << i << " Aromatic Rings\n";
  }

  for (int i = 0; i < _acc_alerts.number_elements(); ++i) {
    if (_acc_alerts[i] == 0) {
      continue;
    }
    output << _acc_alerts[i] << " molecules had " << i << " Alerts\n";
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  struct QEDProperties properties;
  if (! _qed.CalculateProperties(m, properties)) {
    ++_failed_calculations;
    if (_ignore_failed_calculations) {
      return 1;
    } else {
      return 0;
    }
  }

  output << m.name();

  if (_write_descriptors) {
    output << _output_separator <<  properties.amw;
    output << _output_separator <<  properties.alogp;
    output << _output_separator <<  properties.hba;
    output << _output_separator <<  properties.hbd;
    output << _output_separator <<  properties.psa;
    output << _output_separator <<  properties.rotb;
    output << _output_separator <<  properties.arom;
    output << _output_separator <<  properties.alerts;
  }

  if (_verbose) {
    _acc_amw.extra(properties.amw);
    _acc_alogp.extra(properties.alogp);
    ++_acc_hba[properties.hba];
    ++_acc_hbd[properties.hbd];
    _acc_psa.extra(properties.psa);
    ++_acc_rotb[properties.rotb];
    ++_acc_arom[properties.arom];
    ++_acc_alerts[properties.alerts];
  }

  float result = _qed.ComputeQed(properties);
  output << _output_separator << result << '\n';

  output.write_if_buffer_holds_more_than(4092);

  return 1;
}

int
QEDMain(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
QEDMain(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! QEDMain(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
QEDMain(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "QEDMain:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return QEDMain(options, input, output);
}

int
QEDMain(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:afo:M:L:Q:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  options.WriteHeader(output);

  for (const char * fname : cl) {
    if (! QEDMain(options, fname, input_type, output)) {
      cerr << "QEDMain::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace qed

int
main(int argc, char ** argv) {

  int rc = qed::QEDMain(argc, argv);

  return rc;
}
