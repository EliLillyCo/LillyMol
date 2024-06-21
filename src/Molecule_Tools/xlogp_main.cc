// Interface to xlogp

#include <cctype>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/xlogp.h"

namespace xlogp {

using std::cerr;

// By convention the Usage function tells how to use the tool.
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
  cerr << R"(Implementation of xlogp, Wang, Fu, Lai\n";
 -X <fname>  XLogP::XlogpParameters text proto with updated fragment contribution values.
 -3          use version 3 xlogp computation
 -M          use simplistic equation from Mannhold, Poda, Ostermann and Tetek
             logp = 1.46 + 0.11*carbon + - 0.11*heteratoms
 -J ...      fingerprint output, enter '-J help' for info
 -f          function as a TDT filter
 -y          display atom assignments - helpful for debugging
 -p <n>      number of bit replicates when generating fingerprint
 -S <fname>  write labelled smiles to <fname>
 -Y ...      other options, enter '-Y help' for info
 -v          verbose output
)";
// clang-format on

  ::exit(rc);
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    int _version3 = 0;

    int _mannhold = 0;

    int _molecules_read = 0;
    int _successful_calculations = 0;
    int _failed_calculations = 0;

    int _produce_descriptor_file = 1;

    IWString _output_separator = ' ';

    IWString_and_File_Descriptor _stream_for_labelled_molecules;

    IW_STL_Hash_Map_int _unclassified, _failed;

    IWString_and_File_Descriptor _stream_for_failed_calculations;

    IWString _smiles_tag;
    IWString _identifier_tag;
    IWString _xlogp_tag;
    int _bit_replicates;

    int _function_as_tdt_filter;

    int _flush_after_every_molecule;

    Accumulator<double> _acc;
    extending_resizable_array<int> _bucket;

    // Private functions
    int ProcessSuccessfulCalculation(Molecule& m, const int* status,
                double result,
                IWString_and_File_Descriptor& output);
    int ProcessFailedCalculaton(Molecule& m, const int* status);

    int PerformAnyFirstMoleculeRelated(IWString_and_File_Descriptor& output);

    void UpdateBucket(double result);
    int WriteFingerprint(Molecule& m,
                         const Sparse_Fingerprint_Creator& sfc,
                         IWString_and_File_Descriptor& output);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    const IWString& smiles_tag() const {
      return _smiles_tag;
    }

    int function_as_tdt_filter() const {
      return _function_as_tdt_filter;
    }
    int flush_after_every_molecule() const {
      return _flush_after_every_molecule;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 1;
  _remove_chirality = 0;
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
  _xlogp_tag = "NCXLOGP<";
  _function_as_tdt_filter = 0;
  _bit_replicates = 9;  // same default as clogp
  _flush_after_every_molecule = 0;
}

void
DisplayDashJOptions(std::ostream& output) {
  output << " -J <tag>          tag for fingerprints\n";
  output << " -J rep=<n>        number of bit replicates\n";
  ::exit(0);
}

void
DisplayDashYOptions(std::ostream& output) {
  output << " -Y flush          flush output after each molecule\n";
  ::exit(0);
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', _xlogp_tag);

    _xlogp_tag.EnsureEndsWith('<');

    _produce_descriptor_file = 0;
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _bit_replicates) || _bit_replicates < 1) {
      cerr << "The bit replicates (-p) option must be a while +ve integer\n";
      return 0;
    }

    _produce_descriptor_file = 0;
  }

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    _produce_descriptor_file = 0;
  }

  if (cl.option_present('3')) {
    _version3 = 1;
    if (_verbose) {
      cerr << "Will compute version 3 xlogp\n";
    }
  }

  if (cl.option_present('M')) {
    _mannhold = 1;
    if (_verbose) {
      cerr << "Will use the simplistic version of Mannhold et al\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "help") {
        DisplayDashYOptions(cerr);
      } else if (y == "flush") {
        _flush_after_every_molecule = 1;
        if (_verbose) {
          cerr << "Will flush after each molecule read\n";
        }
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    if (! fname.ends_with(".smi")) {
      fname << ".smi";
    }
    if (! _stream_for_failed_calculations.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for failed calculations\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Failed calculations written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! fname.ends_with(".smi")) {
      fname << ".smi";
    }
    if (! _stream_for_labelled_molecules.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for labelled smiles\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Labelled smiles written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
IstopLastZero(const extending_resizable_array<int>& values) {
  int rc = 0;
  for (int i = 0; i < values.number_elements(); ++i) {
    if (values[i] > 0) {
      rc = i + 1;
    }
  }

  return rc;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  output << _acc.n() << " successful " << _failed_calculations << " failed\n";
  output << "Values btw " << _acc.minval() << " and " << _acc.maxval() <<
            " mean " << static_cast<float>(_acc.average()) << '\n';

  const int bstop = IstopLastZero(_bucket);
  uint32_t tot = 0;
  for (int i = 0; i < bstop; ++i) {
    tot += _bucket[i];
  }

  constexpr char kTab = '\t';
  if (_version3)  {
    output << "xlogp3";
  } else {
    output << "xlogp2";
  }
  output << kTab << 'N' << kTab << "Fraction\n";
  for (int i = 0; i < bstop; ++i) {
    output << (-5 + i) << kTab << _bucket[i] << kTab <<
           static_cast<float>(_bucket[i]) / static_cast<float>(tot) << '\n';
  }

  for (const auto& [smt, count] : _failed) {
    output << count << " failed " << smt << '\n';
  }

  for (const auto& [smt, count] : _unclassified) {
    output << count << " unclassified " << smt << '\n';
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  m.remove_all(1);  // Always

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

void
AppendSpaceSuppressedName(const IWString& name,
                          IWString_and_File_Descriptor& output) {
  for (char c : name) {
    if (std::isspace(c)) {
      output << '_';
    } else {
      output << c;
    }
  }
}

// copied from clogp2descriptors_biobyte
int
convert_computed_to_positive_int(float f) {
  int rc = static_cast<int>(f + 5.4999F);

  if (rc <= 0) {
    return 1;
  } else {
    return rc;
  }
}

int
Options::WriteFingerprint(Molecule& m,
                          const Sparse_Fingerprint_Creator& sfc,
                          IWString_and_File_Descriptor& output) {
  if (! _function_as_tdt_filter) {
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(_xlogp_tag, tmp);

  output << tmp << '\n';

  if (!_function_as_tdt_filter) {
    output << "|\n";
  }

  return 1;
}

// Increment a position in _bucket based on `result`.
// Note that we arbitrarily truncate the range to [-5,10].
void
Options::UpdateBucket(double result) {
  int bucket;
  if (result < -5.0) {
    bucket = 0;
  } else if (result >= 10.0) {
    bucket = 15;
  } else {
    bucket = static_cast<int>(result + 5.0 + 0.49999);
  }

  ++_bucket[bucket];
}

int
Options::ProcessSuccessfulCalculation(Molecule& m,
                const int* status,
                double result,
                IWString_and_File_Descriptor& output) {
  ++_successful_calculations;
  _acc.extra(result);

  UpdateBucket(result);

  if (_stream_for_labelled_molecules.active()) {
    m.set_isotopes(status);
    _stream_for_labelled_molecules << m.smiles() << ' ' << m.name() << '\n';
    _stream_for_labelled_molecules.write_if_buffer_holds_more_than(4096);
  }

  if (_produce_descriptor_file) {
    AppendSpaceSuppressedName(m.name(), output);
    output << _output_separator << static_cast<float>(result) << '\n';
    output.write_if_buffer_holds_more_than(32768);
    return 1;
  }

  int int_logp = convert_computed_to_positive_int(result);

  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < _bit_replicates; ++i) {
    sfc.hit_bit(i, int_logp);
  }

  return WriteFingerprint(m, sfc, output);
}

int
Options::ProcessFailedCalculaton(Molecule& m,
                        const int* status) {
  ++_failed_calculations;

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i] > 0) {
      continue;
    }
    if (status[i] < 0) {
      ++_failed[m.smarts_equivalent_for_atom(i)];
      cerr << "NEgative status\n";
    } else if (status[i] == 0) {
      ++_unclassified[m.smarts_equivalent_for_atom(i)];
      cerr << "Status 0 " << m.smiles() << ' ' << m.name() << ' ' << m.smarts_equivalent_for_atom(i) << '\n';
    }
  }

  if (_stream_for_failed_calculations.is_open()) {
    for (int i = 0; i < matoms; ++i) {
      if (status[i] > 0) {
        continue;
      }
      m.set_isotope(i, 1);
    }
    _stream_for_failed_calculations << m.aromatic_smiles() << _output_separator << m.name() << '\n';
    _stream_for_failed_calculations.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
Options::PerformAnyFirstMoleculeRelated(IWString_and_File_Descriptor& output) {
  if (_produce_descriptor_file) {
    output << "Id" << _output_separator;
    if (_version3) {
      output << "xlogp3\n";
    } else if (_mannhold) {
      output << "mannhold\n";
    } else {
      output << "xlogp2\n";
    }
  }

  return 1;
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  if (_molecules_read == 0) {
    PerformAnyFirstMoleculeRelated(output);
  }

  ++_molecules_read;

  std::unique_ptr<int[]> status(new_int(m.natoms()));
  std::optional<double> x;
  if (_mannhold) {
    const int ncarbon = m.natoms(6);
    const int heteroatoms = m.natoms() - ncarbon;
    x = 1.46 + 0.11 * ncarbon - 0.11 * heteroatoms;
  } else {
    x = XLogP(m, status.get());
  }

  if (! x) {
    return ProcessFailedCalculaton(m, status.get());
  }

  return ProcessSuccessfulCalculation(m, status.get(), *x, output);
}

int
XlogPCalculation(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
XlogPCalculation(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! XlogPCalculation(options, *m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4092);
  }

  return 1;
}

int
XlogPCalculationRecord(Options& options,
                       const_IWSubstring buffer,  // note local copy
                       IWString_and_File_Descriptor& output) {
  if (! buffer.ends_with('>')) {
    return 0;
  }
  buffer.chop();
  buffer.remove_leading_chars(options.smiles_tag().length());

  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "XlogPCalculationRecord:cannot parse '" << buffer << "'\n";
    return 0;
  }

  return XlogPCalculation(options, m, output);
}

int
XlogPCalculation(Options& options,
                 iwstring_data_source& input,
                 IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);

    if (! buffer.starts_with(options.smiles_tag())) {
      continue;
    }
    if (! XlogPCalculationRecord(options, buffer, output)) {
      return 0;
    }

    if (options.flush_after_every_molecule()) {
      output.flush();
    }
  }

  return 1;
}

int
XlogPCalculation(Options& options,
                 const char* fname,
                 IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "XlogPCalculation:cannot open '" << fname << "'\n";
    return 0;
  }

  return XlogPCalculation(options, input, output);
}

int
XlogPCalculation(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "XlogPCalculation:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return XlogPCalculation(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:fJ:U:X:yY:p:3S:MC:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  if (cl.option_present('X')) {
    IWString fname = cl.string_value('X');
    if (! ReadNewFragmentParameters(fname)) {
      cerr << "Cannot process new fragment parameters (-X) '" << fname << "'\n";
      return 1;
    }
  }

  if (cl.option_present('y')) {
    SetDisplayAtomAssignments(1);
    if (verbose) {
      cerr << "Will write atom assigments\n";
    }
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  int function_as_tdt_filter = 0;
  if (cl.option_present('f')) {
    function_as_tdt_filter = 1;
    if (verbose) {
      cerr << "Will work as a TDT filter\n";
    }
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

  if (cl.option_present('C')) {
    const_IWSubstring c;
    for (int i = 0; cl.value('C', c, i); ++i) {
      TurnOffNitroxide();
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (function_as_tdt_filter) {
    for (const char* fname : cl) {
      if (! XlogPCalculation(options, fname, output)) {
        cerr << "XlogPCalculation::fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  } else {
    for (const char * fname : cl) {
      if (! XlogPCalculation(options, fname, input_type, output)) {
        cerr << "XlogPCalculation::fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  if (verbose) {
    output.flush();
    options.Report(cerr);
  }

  return 0;
}

}  // namespace xlogp

int
main(int argc, char ** argv) {

  int rc = xlogp::Main(argc, argv);

  return rc;
}
