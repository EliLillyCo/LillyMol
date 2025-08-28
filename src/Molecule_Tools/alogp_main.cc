// Computes alogp

#include <cctype>
#include <cstdio>
#include <iostream>
#include <limits>
#include <memory>

#include "absl/container/flat_hash_map.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/alogp.h"

namespace alogp {

using std::cerr;

// If we are being called from alogp_optimise we read a set of
// known values and rather than write out a prediction for each
// input molecule, we compute an RMS from the predicted values
// to the experimental values we store here.
class ObservedValues {
  private:
    // A mapping from id to observed value - as read from the file.
    absl::flat_hash_map<IWString, float> _id_to_obs;

    // What should we do if we encounter missing values. These are
    // molecules where we can do a computation, but we do not have
    // an experimental value for them.
    int _ignore_missing_values;

    // When computin gan RMS we need a means of accumulating the
    // squared differences.
    Accumulator<double> _acc;

  // Private functions.
    int BuildRecord(const const_IWSubstring& buffer);

  public:
    ObservedValues();

    int Initialise(Command_Line& cl, char opt);
    int Build(IWString& fname);
    int Build(iwstring_data_source& input);

    size_t size() const { 
      return _id_to_obs.size();
    }

    // Caller has made a prediction for `id`.
    // If we have `id` in _id_to_obs, add `value` to the
    // RMS accumulator.
    // If we do not have data for `id` return 0.
    int Extra(const IWString& id, float value);

    float RMS () const {
      return std::sqrt(_acc.average());
    }
};

ObservedValues::ObservedValues() {
  _ignore_missing_values = 0;
}

int
ObservedValues::Initialise(Command_Line& cl, char flag) {
  IWString fname;

  const int verbose = cl.option_present('v');

  IWString o;
  for (int i = 0; cl.value(flag, o, i); ++i) {
    if (o.starts_with("ignore")) {
      _ignore_missing_values = 1;
      if (verbose) {
        cerr << "Will ignore predictions without an experimental value\n";
      }
    } else {
      fname = o;
    }
  }

  if (fname.empty()) {
    cerr << "ObservedValues::Initialise:no file name specified\n";
    return 0;
  }

  if (verbose) {
    cerr << "Reading observed data values from '" << fname << "'\n";
  }

  return Build(fname);
}

int
ObservedValues::Build(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ObservedValues::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
ObservedValues::Build(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (input.lines_read() == 1) {
      continue;
    }
    if (! BuildRecord(buffer)) {
      cerr << "ObservedValues::Build:invalid input '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_to_obs.size();
}

int
ObservedValues::BuildRecord(const const_IWSubstring& buffer) {
  IWString id;
  const_IWSubstring token;

  int i = 0;
  if (! buffer.nextword(id, i) || id.empty() ||
      ! buffer.nextword(token, i) || token.empty()) {
    cerr << "ObservedValues::BuildRecord:invalid input\n";
    return 0;
  }

  float v;
  if (! token.numeric_value(v)) {
    cerr << "ObservedValues::BuildRecord:invalid numeric '" << token << "'\n";
    return 0;
  }

  auto iter = _id_to_obs.find(id);
  if (iter != _id_to_obs.end()) {
    cerr << "ObservedValues::BuildRecord:duplicte data for '" << id << "' ignored\n";
    return 0;
  }

  _id_to_obs[id] = v;

  return 1;
}

// Accumulate the square of the difference between `value` and _id_to_obs[id].
int
ObservedValues::Extra(const IWString& id, float value) {
  auto f = _id_to_obs.find(id);

  if (f != _id_to_obs.end()) {
    // Great, found an experimental value.
  } else if (_ignore_missing_values) {
    return 1;
  } else {
    cerr << "ObservedValues::Extra:no data for '" << id << "'\n";
    return 0;
  }

  float diff = f->second - value;
  _acc.extra(diff * diff);

  return 1;
}

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
  cerr << R"(Implementation of alogp, Wildman and Crippen\n";
 -J ...      fingerprint output, enter '-J help' for info.
 -f          function as a TDT filter.
 -p <n>      number of bit replicates when generating fingerprints.
 -Y ...      other ALogP options, enter '-Y help' for info
 -d          include logD in calculations - must also provide a charge assigner via the -N option.
 -N ...      charge assigner options, enter '-N help'.
 -C <fname>  an alogp::AlogpParameters textproto containing atomic contributions.
 -O <fname>  being used from alogp_optimise.jl. File with known values and only overall RMS is written to stdout.
 -m <min>    do NOT write molecules with computed values below <min>.
 -M <max>    do NOT write molecules with computed values above <max>.
 -l          reduce to largest fragment.
 -g ...      chemical standardisation options.
 -U <fname>  write molecules with failed calculations to <fname>.
 -o ...      options controlling the output - enter '-o help' for info.
 -v          verbose output.
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

    ALogP _alogp;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    uint64_t _molecules_read = 0;
    uint64_t _successful_calculations = 0;
    uint64_t _failed_calculations = 0;

    int _produce_descriptor_file = 1;

    IWString _output_separator = ' ';
    IWString _missing_value = '.';

    IWString_and_File_Descriptor _stream_for_labelled_molecules;

    IW_STL_Hash_Map_int _unclassified, _failed;

    IWString_and_File_Descriptor _stream_for_failed_calculations;

    // Sometimes a value must be produced, even if the calculation has
    // failed.
    double _default_failed_value = -10.0;

    IWString _smiles_tag;
    IWString _identifier_tag;
    IWString _alogp_tag;
    int _bit_replicates;

    // We can also write as a fixe width fingerprint where we set
    // bits according to the value of alogp.
    // Divide the range into intervals. For a given alogp value,
    // figure out the interval. Then set all bits up to that interval.
    // // Think of this as similar to a thermometer
    int _write_as_fixed_width_fingerprint;
    // We can also write all those 0's and 1's as descriptors.
    int _write_as_fixed_width_fingerprint_as_descriptors;

    int _function_as_tdt_filter;

    int _flush_after_every_molecule;

    Accumulator<double> _acc;
    extending_resizable_array<int> _bucket;

    Charge_Assigner _charge_assigner;
    int _rerun_charge_assigner_queries;

    // Computing logD is optional.
    int _compute_logd;

    // If we are being called from alogp_optimise we just accumulate
    // squared differences in this object.
    std::unique_ptr<ObservedValues> _obs;

    // We can also filter outinput molecules only writing those within specified bounds.
    float _min_value;
    float _max_value;
    uint64_t _molecules_below_min;
    uint64_t _molecules_above_max;

    // Private functions
    int ProcessSuccessfulCalculation(Molecule& m,
                double alogp,
                IWString_and_File_Descriptor& output);
    int ProcessFailedCalculaton(Molecule& m, IWString_and_File_Descriptor& output);

    int MaybeWriteSmiles(Molecule& m,
                         IWString_and_File_Descriptor& output) const;
    int PerformAnyFirstMoleculeRelated(IWString_and_File_Descriptor& output);

    std::optional<double> LogD(Molecule& m, double logp);
    double RerunChargeAssigner(Molecule& m, formal_charge_t* fc);

    void UpdateBucket(double result);
    int WriteFingerprint(Molecule& m,
                         const Sparse_Fingerprint_Creator& sfc,
                         IWString_and_File_Descriptor& output);
    int WriteFixedWidthFingerprint(Molecule& m, double result, 
                        IWString_and_File_Descriptor& output);
    int WriteFixedWidthDescriptors(const IWString& mname, int bucket,
                                IWString_and_File_Descriptor& output) const;

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
  
    // If we are being called from alogp_optimise the final task is to
    // write our accumlated RMS.
    int WriteRMSIfRequested(IWString_and_File_Descriptor& output) const;

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 1;
  _remove_chirality = 0;
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
  _alogp_tag = "NCALOGP<";
  _function_as_tdt_filter = 0;
  _bit_replicates = 9;  // same default as clogp
  _flush_after_every_molecule = 0;
  _write_as_fixed_width_fingerprint = 0;
  _write_as_fixed_width_fingerprint_as_descriptors = 0;

  _compute_logd = 0;

  _rerun_charge_assigner_queries = 0;

  _min_value = -std::numeric_limits<float>::max();
  _max_value = std::numeric_limits<float>::max();

  _molecules_below_min = 0;
  _molecules_above_max = 0;
}

void
DisplayDashJOptions(std::ostream& output) {
  output << R"( -J <tag>          tag for fingerprints
 -J fixed=<nn>     generate 'thermometer-like' fixed width fingerprints of width <nn> bits.
                   The number of bits set is a function of the alogp value.
                   Low values set hardly any of the <nn> bits high values set all of them.
)";

  ::exit(0);
}

void
DisplayDashYOptions(std::ostream& output) {
  output << " -Y flush          flush output after each molecule\n";
  output << " -Y label          label atoms with atom type assigned\n";
  output << " -Y quiet          suppress warning messages about unassigned atoms\n";
  output << " -Y alcacid        use the alcohol atom type for certain acids - rdkit compat\n";
  output << " -Y RDKIT.N+       treat N+ atoms same as RDKit\n";
  output << " -Y RDKIT.HP       treat phosphoric acids same as RDKit\n";
  output << " -Y ZWIT           explicit treatment of Zwitterions\n";
  output << " -Y okx            do NOT fail on encountering an unclassified atom. Values set to 0.0\n";

  ::exit(0);
}

void
DisplayDashoOptions(std::ostream& output) {
  output << " -o osep=<char>    output separator, default space: -o sep=tab, -o sep= ...\n";
  output << " -o missing=<char> token written for missing values, default '.'\n";

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
    IWString j;
    int got_tag = 0;
    for (int i = 0; cl.value('J', j, i); ++i) {
      if (j.starts_with("fixed=")) {
        j.remove_leading_chars(6);
        if (! j.numeric_value(_write_as_fixed_width_fingerprint) ||
                _write_as_fixed_width_fingerprint < 16) {
          cerr << "The width of a fixed width fingerprint must be at least 16 bits\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will write as fixed width fingerprints with " << 
                  _write_as_fixed_width_fingerprint << " bits\n";
        }
      } else if (j.starts_with("desc=")) {
        j.remove_leading_chars(5);
        if (! j.numeric_value(_write_as_fixed_width_fingerprint) ||
                _write_as_fixed_width_fingerprint < 16) {
          cerr << "The width of a fixed width fingerprint must be at least 16 bits\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will write fixed width fp as binary descriptors with " << 
                  _write_as_fixed_width_fingerprint << " columns\n";
        }
        _write_as_fixed_width_fingerprint_as_descriptors = 1;
      } else if (j == "help") {
        DisplayDashJOptions(cerr);
      } else {
        got_tag = 1;
        _alogp_tag = j;
      }
    }

    if (_write_as_fixed_width_fingerprint && ! got_tag) {
      _alogp_tag = "FPALP<";
    } else if (_write_as_fixed_width_fingerprint && ! _alogp_tag.starts_with("FP")) {
      cerr << "Fixed length fingerprints must start with 'FP', '" << _alogp_tag << "' invalid\n";
      return 0;
    }

    _alogp_tag.EnsureEndsWith('<');

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

  if (cl.option_present('N')) {
    if (! _charge_assigner.construct_from_command_line(cl, _verbose, 'N')) {
      cerr << "Options::Initialise:cannot initialise charge assigner (-N)\n";
      return 0;
    }
  }

  // Read this here before we process any command line options.
  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    if (! _alogp.ReadConfiguration(fname)) {
      cerr << "Options::Initialise:cannot read alogp parameters '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Read configuration from '" << fname << "'\n";
    }
  }

  if (cl.option_present('d')) {
    if (! _charge_assigner.active()) {
      cerr << "Options::initialise:logD requested (-d) but no charge assigner (-N)\n";
      return 0;
    }

    _compute_logd = 1;
    if (_verbose) {
      cerr << "Will also include logD\n";
    }
  }

  if (cl.option_present('o')) {
    static constexpr bool kMessageIfUnrecognised = false;

    IWString o;
    for (int i = 0; cl.value('o', o, i); ++i) {
      if (o.starts_with("osep=")) {
        o.remove_leading_chars(5);
        char_name_to_char(o, kMessageIfUnrecognised);
        _output_separator = o;
        if (_verbose) {
          cerr << "Output separator '" << _output_separator << "'\n";
        }
      } else if (o.starts_with("missing=")) {
        o.remove_leading_chars(8);
        _missing_value = o;
        char_name_to_char(_missing_value, kMessageIfUnrecognised);
        if (_verbose) {
          cerr << "Missing value set to '" << _missing_value << '\n';
        }
      } else if (o == "help") {
        DisplayDashoOptions(cerr);
      } else {
        cerr << "Unrecognised -o directive '" << o << "'\n";
        return 0;
      }
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
      } else if (y == "label") {
        _alogp.set_label_with_atom_type(1);
        if (_verbose) {
          cerr << "Will label molecules with atom type\n";
        }
      } else if (y == "quiet") {
        _alogp.set_display_error_messages(0);
        if (_verbose) {
          cerr << "Will NOT display unclassified atom error messages\n";
        }
      } else if (y == "alcacid") {
        _alogp.set_use_alcohol_for_acid(1);
        if (_verbose) {
          cerr << "Will use the alcohol H atom type for acids\n";
        }
      } else if (y == "RDKIT.N+") {
        _alogp.set_rdkit_charged_nitrogen(1);
        if (_verbose) {
          cerr << "Will use the RDKit H count calculation on charged nitrogen atoms\n";
        }
      } else if (y == "RDKIT.HP") {
        _alogp.set_rdkit_phoshoric_acid_hydrogen(1);
        if (_verbose) {
          cerr << "Will treat OH groups on phosphoric acids like RDKIt\n";
        }
      } else if (y == "ZWIT") {
        _alogp.set_apply_zwitterion_correction(1);
      } else if (y == "okx") {
        _alogp.set_fail_if_unclassified_atom(0);
        if (_verbose) {
          cerr << "Will ignore unclassified atoms\n";
        }
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('O')) {
    _obs = std::make_unique<ObservedValues>();
    if (! _obs->Initialise(cl, 'O')) {
      cerr << "Cannot initialised observed values (-O)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Read " << _obs->size() << " experimental values\n";
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_value)) {
      cerr << "Options::Initialise:invalid min value (-m)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will not write molecules with clogp values below " << _min_value << '\n';
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_value)) {
      cerr << "Options::Initialise:invalid max value (-M)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will not write molecules with clogp values above " << _max_value << '\n';
    }
  }

  if (_min_value > _max_value) {
    cerr << "Options::Initialise:inconsistent min " << _min_value << " max " << _max_value << " values\n";
    return 0;
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
    fname.EnsureEndsWith(".smi");
    if (! _stream_for_labelled_molecules.open(fname)) {
      cerr << "Options::Initialise:cannot open stream for labelled molecules '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Labelled molecules written to '" << fname << "'\n";
    }
    _alogp.set_label_with_atom_type(1);
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
  if (_acc.n() > 0) {
  output << "clogp values btw " << _acc.minval() << " and " << _acc.maxval() <<
            " mean " << static_cast<float>(_acc.average()) << '\n';
  }

  const int bstop = IstopLastZero(_bucket);
  uint32_t tot = 0;
  for (int i = 0; i < bstop; ++i) {
    tot += _bucket[i];
  }

  constexpr char kTab = '\t';
  output << "alogp" << kTab << 'N' << kTab << "Fraction\n";
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

  if (_molecules_below_min) {
    cerr << _molecules_below_min << " below " << _min_value << " discarded\n";
  }
  if (_molecules_above_max) {
    cerr << _molecules_above_max << " above " << _max_value << " discarded\n";
  }

  return 1;
}

int
Options::WriteRMSIfRequested(IWString_and_File_Descriptor& output) const {
  if (! _obs) {
    return 0;
  }

  output << _obs->RMS();  // no newline is deliberate.

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
AppendFirstTokenOfName(const IWString& name,
                          IWString_and_File_Descriptor& output) {
  for (char c : name) {
    if (std::isspace(c)) {
      return;
    } else {
      output << c;
    }
  }
}

#ifdef BUCKET_WIDTH_ONE
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
#endif

#define EXPERIMENT_WITH_NARROWER_BUCKETS
#ifdef EXPERIMENT_WITH_NARROWER_BUCKETS
// Divide into half log unit buckets.
int
convert_computed_to_positive_int(float f) {
  if (f <= -5.0f) {
    return 1;
  }

  if (f >= 10.0f) {
    return 20;
  }

  // convert to range 0,10
  f += 5.0f;

  return static_cast<int>(f + f + 0.4999f);

  int rc = static_cast<int>(f + 5.4999F);

  if (rc <= 0) {
    return 1;
  } else {
    return rc;
  }
}
#endif

int
Options::MaybeWriteSmiles(Molecule& m,
                         IWString_and_File_Descriptor& output) const {
  if (_obs) {
    return 1;
  }

  if (_function_as_tdt_filter) {
    return 1;
  }

  if (_write_as_fixed_width_fingerprint_as_descriptors) {
    return 1;
  }

  output << _smiles_tag << m.smiles() << ">\n";
  output << _identifier_tag << m.name() << ">\n";

  return 1;
}

int
Options::WriteFixedWidthFingerprint(Molecule& m, double result, 
                        IWString_and_File_Descriptor& output) {
  MaybeWriteSmiles(m, output);

  static constexpr double kMinval = -1.0;
  static constexpr double kMaxval = 7.0;

  int bucket;
  if (result <= kMinval) {
    bucket = 0;
  } else if (result >= kMaxval) {
    bucket = 10;
  } else {
    bucket = static_cast<int>((result - kMinval) / (kMaxval - kMinval) * _write_as_fixed_width_fingerprint);
  }

  assert(bucket < _write_as_fixed_width_fingerprint);

#ifdef DEBUG_FIXED_WIDTH_FINGERPRINT
  cerr << result << " in bucket " << bucket << '\n';
#endif

  if (_write_as_fixed_width_fingerprint_as_descriptors) {
    return WriteFixedWidthDescriptors(m.name(), bucket, output);
  }

  IW_Bits_Base bits(_write_as_fixed_width_fingerprint);
  for (int i = 0; i < bucket; ++i) {
    bits.set_bit(i, 1);
  }

  IWString fp;
  bits.daylight_ascii_representation_including_nset_info(fp);

  output << _alogp_tag << fp << ">\n";

  if (! _function_as_tdt_filter) {
    output << "|\n";
  }

  return 1;
}

int
Options::WriteFixedWidthDescriptors(const IWString& mname, int bucket,
                                IWString_and_File_Descriptor& output) const {
  output << mname;

  for (int i = 0; i <= bucket; ++i) {
    output << _output_separator << '1';
  }
  for (int i = bucket + 1; i < _write_as_fixed_width_fingerprint; ++i) {
    output << _output_separator <<  '0';
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::WriteFingerprint(Molecule& m,
                          const Sparse_Fingerprint_Creator& sfc,
                          IWString_and_File_Descriptor& output) {
  MaybeWriteSmiles(m, output);

  if (_write_as_fixed_width_fingerprint) {
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(_alogp_tag, tmp);

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
                double alogp,
                IWString_and_File_Descriptor& output) {
  ++_successful_calculations;
  _acc.extra(alogp);

  if (alogp < _min_value) {
    ++_molecules_below_min;
    return 1;
  } else if (alogp > _max_value) {
    ++_molecules_above_max;
    return 1;
  }

  UpdateBucket(alogp);

  std::optional<double> logd;
  if (_compute_logd) {
    logd = LogD(m, alogp);
  }

  if (_obs) {
    return _obs->Extra(m.name(), alogp);
  }

  if (_produce_descriptor_file) {
    AppendFirstTokenOfName(m.name(), output);
    output << _output_separator;
    static char buffer[32];
    sprintf(buffer, "%.3f", static_cast<float>(alogp));
    output << buffer;

    if (! _compute_logd) {
    } else if (logd) {
      sprintf(buffer, "%.3f", static_cast<float>(*logd));
      output << _output_separator << buffer;
    } else {
      output << _output_separator << _missing_value;
    }
    output << '\n';

    output.write_if_buffer_holds_more_than(32768);

    return 1;
  }

  if (_write_as_fixed_width_fingerprint) {
    return WriteFixedWidthFingerprint(m, alogp, output);
  }

  const int int_logp = convert_computed_to_positive_int(alogp);
  
  Sparse_Fingerprint_Creator sfc;

  for (int i = 0; i < _bit_replicates; ++i) {
    sfc.hit_bit(i, int_logp);
  }

  if (logd) {
    const int int_logd = convert_computed_to_positive_int(*logd);
    for (int i = 0; i < _bit_replicates; ++i) {
      sfc.hit_bit(_bit_replicates + i, int_logd);
    }
  }

  return WriteFingerprint(m, sfc, output);
}

// Return true if the molecule has a positive charge and no negative charges.
int
IsQuat(const Molecule& m) {
  const int matoms = m.natoms();

  int nneg = 0;
  int got_quat = 0;
  for (int i = 0; i < matoms; ++i) {
    formal_charge_t fc = m.formal_charge(i);
    if (fc == 0) [[ likely ]] {
      continue;
    }
    if (fc < 0) {
      ++nneg;
      continue;
    }

    if (m.atomic_number(i) == 7 && m.ncon(i) == 4) {
      got_quat = 1;
    }
  }

  return got_quat && nneg == 0;
}

std::optional<double>
Options::LogD(Molecule& m, double logp) {

  static constexpr double kNeutralUncharged = 0.34;
  static constexpr double kQuat = -0.29;
  static constexpr double kMultiMinus = 3.09;
  static constexpr double kPositive3 = 3.38;
  static constexpr double kPositive2 = 2.14;
  static constexpr double kZwit1Minus = 0.71;
  static constexpr double kZwit2Minus = 2.10;

  // cerr << "computing logD, quat " << IsQuat(m) << '\n';
  // Assume that anything with a positive charge coming in makes this a quat.
  if (IsQuat(m)) {
    return logp - kQuat;
  }

  std::vector<ChargeAndQuery> charges;
  if (! _charge_assigner.Process(m, charges)) {
    return logp - kNeutralUncharged;
  }

  int npos = 0;
  int nneg = 0;
  for (const ChargeAndQuery& afq : charges) {
    if (_verbose > 2) [[ unlikely ]] {
      cerr << "Got match " << afq << '\n';
    }

    if (afq.formal_charge < 0) {
      ++nneg;
    } else {
      ++npos;
    }
  }

  cerr << "Find nneg " << nneg << " npos " << npos << " logp " << logp << '\n';

  if (0 == npos && 0 == nneg) [[unlikely]] {  // really should have been picked up before
    return logp - kNeutralUncharged;
  }

  if (nneg > 1 && 0 == npos) {
    return logp - kMultiMinus;
  }

  if (npos >= 3 && 0 == nneg) {
    return logp - kPositive3;
  }

  if (2 == npos && 0 == nneg) {
    return logp - kPositive2;
  }

  if (npos > 0 && 1 == nneg) {
    return logp - kZwit1Minus;
  }

  if (npos > 0 && nneg > 1) {
    return logp - kZwit2Minus;
  }

  if (npos > 0 && 0 == nneg) {
    ;
  } else if (0 == npos && nneg > 0) {
    ;
  } else {
    cerr << "options::LogD:unusual charge state " << npos << "+ and " << nneg << "-\n";
    return std::nullopt;
  }

  // Doesn't quality as any of the special cases. Compute the offset from
  // the queries that match

  if (_rerun_charge_assigner_queries) {
    std::unique_ptr<formal_charge_t[]> formal_charge = std::make_unique<formal_charge_t[]>(m.natoms());
    double offset = RerunChargeAssigner(m, formal_charge.get());
    // cerr << "Offset from rerun " << offset << '\n';
    return logp + offset;
  }

  double offset = 0.0;
  for (const ChargeAndQuery& afq : charges) {
    double d;
    if (!_charge_assigner[afq.query_number]->numeric_value(d)) {
      continue;
    }

    offset += d;

    if (_verbose > 2) [[ unlikely ]] {
      cerr << "After " << afq << " d " << d << " total offset " << offset << '\n';
    }
  }

  cerr << "Offset " << offset << " returning " << (logp + offset) << '\n';

  return logp + offset;
}

double
Options::RerunChargeAssigner(Molecule& m, formal_charge_t* fc) {
  int matoms = m.natoms();

  std::fill_n(fc, matoms, 0);

  int n = _charge_assigner.number_elements();

  double rc = 0.0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < n; i++) {
    Substructure_Results sresults;

    Substructure_Hit_Statistics* q = _charge_assigner[i];

    int nhits = q->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    double logd_offset;
    if (!q->numeric_value(logd_offset)) {  // hmmm, no logd offset!!
      continue;
    }

    //  cerr << nhits << " to query '" << q->comment() << "' offset " << logd_offset << '\n';

    for (int j = 0; j < nhits; j++) {
      const Query_Atoms_Matched* qam = sresults.query_atoms_matching(j);

      const Set_of_Atoms* e = sresults.embedding(j);

      for (int k = 0; k < qam->number_elements(); k++) {
        atom_number_t l = e->item(k);

        if (0 != fc[l]) {  // already hit by this or something else
          continue;
        }

        const Substructure_Atom* a = qam->item(k);

        double charge_specification;
        if (!a->numeric_value(charge_specification)) {  // no charge placed here
          continue;
        }

        rc += logd_offset;
        if (charge_specification < 0.0) {
          fc[l] = -1;
        } else {
          fc[l] = 1;
        }
      }
    }
  }

  return rc;
}

int
Options::ProcessFailedCalculaton(Molecule& m,
                IWString_and_File_Descriptor& output) {
  ++_failed_calculations;

  if (_stream_for_failed_calculations.is_open()) {
    _stream_for_failed_calculations << m.aromatic_smiles() << _output_separator << m.name() << '\n';
    _stream_for_failed_calculations.write_if_buffer_holds_more_than(32768);
  }

  if (! _function_as_tdt_filter) {
  } else if (_write_as_fixed_width_fingerprint) {
    WriteFixedWidthFingerprint(m, _default_failed_value, output);
  } else {
    Sparse_Fingerprint_Creator sfc;

    IWString tmp;
    sfc.daylight_ascii_form_with_counts_encoded(_alogp_tag, tmp);

    output << tmp << '\n';
  }

  return 1;
}

int
Options::PerformAnyFirstMoleculeRelated(IWString_and_File_Descriptor& output) {
  if (_obs) {
    return 1;
  }

  if (_produce_descriptor_file) {
    output << "Id" << _output_separator << "alogp";
    if (_compute_logd) {
      output << _output_separator << "alogd";
    }
    output << '\n';
  }

  if (_write_as_fixed_width_fingerprint_as_descriptors) {
    output << "Id";
    for (int i = 0; i < _write_as_fixed_width_fingerprint; ++i) {
      output << _output_separator << "ALP" << i;
    }
    output << '\n';
  }

  return 1;
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  if (_molecules_read == 0) {
    PerformAnyFirstMoleculeRelated(output);
  }

  ++_molecules_read;

  std::optional<double> x = _alogp.LogP(m);

  if (! x) {
    return ProcessFailedCalculaton(m, output);
  }

  if (_stream_for_labelled_molecules.is_open()) {
    _stream_for_labelled_molecules << m.smiles() << ' ' << m.name() << ' ' << *x << '\n';
    _stream_for_labelled_molecules.write_if_buffer_holds_more_than(4096);
  }

  return ProcessSuccessfulCalculation(m, *x, output);
}

int
AlogPCalculation(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
AlogPCalculation(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! AlogPCalculation(options, *m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
AlogPCalculationRecord(Options& options,
                       const_IWSubstring buffer,  // note local copy
                       IWString_and_File_Descriptor& output) {
  if (! buffer.ends_with('>')) {
    return 0;
  }
  buffer.chop();
  buffer.remove_leading_chars(options.smiles_tag().length());

  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "AlogPCalculationRecord:cannot parse '" << buffer << "'\n";
    return 0;
  }

  return AlogPCalculation(options, m, output);
}

int
AlogPCalculation(Options& options,
                 iwstring_data_source& input,
                 IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);

    if (! buffer.starts_with(options.smiles_tag())) {
      continue;
    }
    if (! AlogPCalculationRecord(options, buffer, output)) {
      return 0;
    }

    if (options.flush_after_every_molecule()) {
      output.flush();
    }
  }

  return 1;
}

int
AlogPCalculation(Options& options,
                 const char* fname,
                 IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "AlogPCalculation:cannot open '" << fname << "'\n";
    return 0;
  }

  return AlogPCalculation(options, input, output);
}

int
AlogPCalculation(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "AlogPCalculation:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return AlogPCalculation(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:fJ:U:Y:p:N:do:C:S:O:m:M:");

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

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (function_as_tdt_filter) {
    for (const char* fname : cl) {
      if (! AlogPCalculation(options, fname, output)) {
        cerr << "AlogPCalculation::fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  } else {
    for (const char * fname : cl) {
      if (! AlogPCalculation(options, fname, input_type, output)) {
        cerr << "AlogPCalculation::fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  options.WriteRMSIfRequested(output);

  if (verbose) {
    output.flush();
    options.Report(cerr);
  }

  return 0;
}

}  // namespace alogp

int
main(int argc, char ** argv) {

  int rc = alogp::Main(argc, argv);

  return rc;
}
