// Labute Approximate Surface Area descriptors

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/e_state_computation_procedure.h"

namespace labute {

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
  cerr << R"(Computes Approximate Surface Area (ASA) descriptos from Labute, in MOE.
These are surface area weighted atomic properties.
 -Q <charge>    kind of partial charges, enter '-Q help' for info. Default is -Q Gasteiger.
 -X .           atomic properties are E-State values - takes an argument that is currently ignored.
 -J <tag>       generate fingerprints - must start with NC.
 -f             function as a TDT filter when generating fingerprints.
 -d             optionally generate more descriptors that involve computation of the distance matrix.
 -q             ignore failed partial charge calculations - default is to stop.
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

enum class PartialChargeType : int {
  kUnspecified = 0,
  kAbraham = 1,
  kGasteiger = 2,
  kHuckel = 3,
  kGasteigerHuckel = 4,
  kDelRe = 5,
  kPullman = 6,

  // This enum is used as an index into the _buckets array, so we also
  // need to accommodate E-State calculations;
  kEState = 7

  // Keep this in sync with the _buckets array.
  // Original plans were for computing multiple things at once, but
  // that was not implemented.
};

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

    PartialChargeType _partial_charge;

    // For each partial charge type, a vector of cut points.

    std::vector<float> _buckets[8];

    // This array is used to accumulate the contributions witin each range.
    // Renders the class not thread safe.
    float* _bucket_count;
    uint32_t _number_buckets;

    uint64_t _molecules_read = 0;

    // By default, we compute partial charge based descriptors, but we can
    // also do estate based descriptors.
    // SOmewhat of a dilemma. Observe that some RDKit Estate based descriptors
    // are good in tree based models. But the implementation of estate in RDKit 
    // is wrong - it only considers bonded atoms, an obvious bug.
    // But we have what looks like a more correct implementation of e state.
    // But will it be as effective?
    int _compute_estate_values;

    int _ignore_partial_charge_calcucation_failed;

    IWString _tag;
    int _function_as_tdt_filter;
    IWString _smiles_tag;
    IWString _identifier_tag;

    // When writing a fingerprint, we can control the relative
    // weight assigned to various features. These should be settable
    // parameters.
    float _minmax_scale = 20.0f;
    float _variance_scale = 20.0f;
    float _bucket_scale = 10.0f;

    // ASA values are large, so scale down so it does not dominate.
    // 100000 values between 19.257 and 253.343, tot 1.29403e+07 ave 129.403 std 34.4187
    float _asa_scale = 0.1;

    // We can optionally compute the bond separation between the
    // smallest and largest charged atoms.
    int _do_distance_matrix;

    IWString _descriptor_prefix;
    IWString _output_separator;
    
    // Private functions.

    int Labute(Molecule& m, const float* vi, IWString_and_File_Descriptor& output);

    int Abraham(Molecule& m, const float* vi, IWString_and_File_Descriptor& output);
    int Huckel(Molecule& m, const float* vi, IWString_and_File_Descriptor& output);
    int Gasteiger(Molecule& m, const float* vi, IWString_and_File_Descriptor& output);
    int GasteigerHuckel(Molecule& m, const float* vi, IWString_and_File_Descriptor& output);

    int EStateAtomicProperties(Molecule& m,
                                const float* vi,
                                IWString_and_File_Descriptor& output);

    template <typename T>
    int CommonOutput(Molecule& m, const float* vi, const T* atomic, const std::vector<float>& bucket, IWString_and_File_Descriptor& output);
    int FingerprintOutput(Molecule& m,
                           const float asa,
                           const Accumulator<double>& acc_charge,
                           const float total_abs,
                           int max_min_distance,
                           IWString_and_File_Descriptor& output);

  public:
    Options();
    ~Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    const IWString& smiles_tag() const {
      return _smiles_tag;
    }

    int verbose() const {
      return _verbose;
    }

    int WriteHeader(IWString_and_File_Descriptor& output) const;

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _ignore_partial_charge_calcucation_failed = 0;
  _bucket_count = nullptr;
  _number_buckets = 0;

  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";

  _function_as_tdt_filter = 0;

  _do_distance_matrix = 0;

  _compute_estate_values = 0;

  _descriptor_prefix = "asa_";
  _output_separator = ' ';

  // Gasteiger
  _buckets[static_cast<int>(PartialChargeType::kAbraham)] = { 1.80f, 3.0, 5.00f, 7.00f };
  _buckets[static_cast<int>(PartialChargeType::kGasteiger)] = { -0.30f, -0.25f, 0.08f, -0.01f, 0.05f, 0.10f, 0.15f, 0.25f, 0.5f };
  _buckets[static_cast<int>(PartialChargeType::kHuckel)] = { -0.50f, -0.25f, -0.15, 0.2f, 0.4f };
  _buckets[static_cast<int>(PartialChargeType::kGasteigerHuckel)] = { -0.50f, -0.30, -0.20f, -0.10f, 0.01f, 0.02f, 0.3f, 0.4f };
  _buckets[static_cast<int>(PartialChargeType::kEState)] = { -10.0, -2.0f, 3.0f, 7.0f, 10.0f, 15.0f, 20.0f};

  uint32_t max_size = 0;
  for (const auto& b : _buckets) {
    if (b.size() > max_size) {
      max_size = b.size();
    }
  }

  _bucket_count = new float[max_size + 1];
  _number_buckets = max_size;
}

Options::~Options() {
  if (_bucket_count != nullptr) {
    delete [] _bucket_count;
  }
}

void
DisplayPartialChargeSpecifiers(std::ostream& output) {
  output << "The following partial charge specifications are supported\n";
  output << " -Q Abraham\n";
  output << " -Q Gasteiger\n";
  output << " -Q Huckel\n";
  output << " -Q GasteigerHuckel\n";
//output << " -Q DelRe\n";
//output << " -Q Pullman\n";
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

  if (cl.option_present('Q') && cl.option_present('X')) {
    cerr << "Sorry, only computes either partial charge (-Q) or E-State descriptors (-X)\n";
    return 0;
  }

  if (cl.option_present('X')) {
    _compute_estate_values = 1;
    if (_verbose) {
      cerr << "Atomic values are E-State values\n";
    }
  } else if (cl.option_present('Q')) {
    const const_IWSubstring q = cl.string_value('Q');
    if (q == "Abraham" || q == "A") {
      _partial_charge = PartialChargeType::kAbraham;
    } else if (q == "Gasteiger" || q == "G") {
      _partial_charge = PartialChargeType::kGasteiger;
    } else if (q == "Huckel" || q == "H") {
      _partial_charge = PartialChargeType::kHuckel;
    } else if (q == "GasteigerHuckel" || q == "GH") {
      _partial_charge = PartialChargeType::kGasteigerHuckel;
//  } else if (q == "DelRe") {
//    _partial_charge = PartialChargeType::kDelRe;
//  } else if (q == "Pullman") {
//    _partial_charge = PartialChargeType::kPullman;
    } else if (q == "help") {
      DisplayPartialChargeSpecifiers(cerr);
      return 1;
    } else {
      cerr << "Unrecognised partial charge specification '" << q << "'\n";
      DisplayPartialChargeSpecifiers(cerr);
      return 0;
    }
  } else {
    _partial_charge = PartialChargeType::kGasteiger;
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    _tag.EnsureEndsWith('<');
    if (! _tag.starts_with("NC")) {
      cerr << "The tag must start with 'NC', '" << _tag << "' invalid\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Fingerprint output " << _tag << "\n";
    }
  }

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    if (_verbose) {
      cerr << "Will work as a tdt filter\n";
    }
  }

  if (cl.option_present('q')) {
    _ignore_partial_charge_calcucation_failed = 1;
    if (_verbose) {
      cerr << "Will ignore molecules with failed partial charge calculations\n";
    }
  }

  if (cl.option_present('d')) {
    _do_distance_matrix = 1;
    if (_verbose) {
      cerr << "Will compute distance matrix features\n";
    }
  }

  return 1;
}

int
Options::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "ID";
  output << _output_separator << _descriptor_prefix << "ASA";
  output << _output_separator << _descriptor_prefix << "total_abs";
  output << _output_separator << _descriptor_prefix << "var";
  output << _output_separator << _descriptor_prefix << "minval";
  output << _output_separator << _descriptor_prefix << "maxval";
  if (_do_distance_matrix) {
    output << _output_separator << _descriptor_prefix << "max_to_min";
  }

  for (uint32_t i = 0; i < _number_buckets; ++i) {
    output << _output_separator << _descriptor_prefix << "vsa" << i;
  }

  output << '\n';

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  // Other information about what has happened.

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

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

// `buckets` is an array of break points. Return the index of the
// first item for which `value` is greater
int
FindBucket(float value, const std::vector<float>& buckets) {
  // cerr << "Looking for " << value << " in " << buckets.size() << " items\n";
  if (value <= buckets[0]) [[unlikely]] {
    return 0;
  }

  if (value >= buckets.back()) [[unlikely]] {
    return buckets.size() - 1;
  }

  return std::upper_bound(buckets.begin(), buckets.end(), value) - buckets.begin();
}

// Return the number of bonds between the least partial charged atom
// the the largest partial charged atom.
template <typename T>
int
DistanceBetweenMinMax(Molecule& m, const T* atomic) {
  m.recompute_distance_matrix();

  const int matoms = m.natoms();

  T min_charge = atomic[0];
  atom_number_t amin = 0;
  T max_charge = atomic[0];
  atom_number_t amax = 0;

  for (int i = 1; i < matoms; ++i) {
    const T q = atomic[i];
    if (q < min_charge) {
      min_charge = q;
      amin = i;
    } else if (q > max_charge) {
      max_charge = q;
      amax = i;
    }
  }

  return m.bonds_between(amin, amax);
}

template <typename T>
int
Options::CommonOutput(Molecule& m, const float* vi,
                      const T* atomic,
                      const std::vector<float>& bucket,
                      IWString_and_File_Descriptor& output) {
  // cerr << "Filling into " << bucket.size() << " ranges\n";
  std::fill_n(_bucket_count, _number_buckets, 0.0f);

  const int matoms = m.natoms();

// #define ECHO_PARTIAL_CHARGES
#ifdef ECHO_PARTIAL_CHARGES
  for (int i = 0; i < matoms; ++i) {
    cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " " << atomic[i] << '\n';
  }
  return 1;
#endif

  Accumulator<double> acc_charge;
  float total_abs = 0.0f;

  for (int i = 0; i < matoms; ++i) {
    T q = atomic[i];
    acc_charge.extra(q);
    if (q >= 0.0f) {
      total_abs += q;
    } else {
      total_abs -= q;
    }

    int b = FindBucket(q, bucket);
    _bucket_count[b] += q * vi[i];
  }

  int max_min_distance = 0;
  if (_do_distance_matrix) {
    max_min_distance = DistanceBetweenMinMax(m, atomic);
  }

  const float asa = std::accumulate(vi, vi + matoms, 0.0f);

  if (_tag.length()) {
    return FingerprintOutput(m, asa, acc_charge, total_abs, max_min_distance, output);
  }

  // First token of name. Perhaps this should be optional.
  const IWString& id = m.name();
  for (char c : id) {
    if (c == ' ') {
      break;
    }
    output << c;
  }

  output << _output_separator << asa;
  output << _output_separator << total_abs;
  output << _output_separator << static_cast<float>(acc_charge.variance());
  output << _output_separator << static_cast<float>(acc_charge.minval());
  output << _output_separator << static_cast<float>(acc_charge.maxval());
  if (_do_distance_matrix) {
    output << _output_separator << max_min_distance;
  }

  for (uint32_t i = 0; i < _number_buckets; ++i) {
    output << _output_separator << _bucket_count[i];
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::FingerprintOutput(Molecule& m,
                           const float asa,
                           const Accumulator<double>& acc_charge,
                           const float total_abs,
                           int max_min_distance,
                           IWString_and_File_Descriptor& output) {
  if (! _function_as_tdt_filter) {
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  static constexpr int kTotalAbs = 0;
  static constexpr int kVariance = 1;
  static constexpr int kMinval = 2;
  static constexpr int kMaxval = 3;
  static constexpr int kMaxToMin = 4;
  static constexpr int kLabute = 5;
  static constexpr int kBStart = 10;

  Sparse_Fingerprint_Creator sfc;
  sfc.hit_bit(kTotalAbs, static_cast<int>(total_abs));
  sfc.hit_bit(kVariance, static_cast<int>(acc_charge.variance() * _variance_scale));
  if (acc_charge.minval() >= 0.0) {
    sfc.hit_bit(kMinval, static_cast<int>(acc_charge.minval() * _minmax_scale));
  } else {
    sfc.hit_bit(kMinval, static_cast<int>(- acc_charge.minval() * _minmax_scale));
  }

  if (acc_charge.maxval() >= 0.0) {
    sfc.hit_bit(kMaxval, static_cast<int>(acc_charge.maxval() * _minmax_scale));
  } else {
    sfc.hit_bit(kMaxval, static_cast<int>(- acc_charge.maxval() * _minmax_scale));
  }

  if (_do_distance_matrix) {
    sfc.hit_bit(kMaxToMin, max_min_distance);
  }

  sfc.hit_bit(kLabute, static_cast<int>(asa * _asa_scale));

  for (uint32_t i = 0; i < _number_buckets; ++i) {
    sfc.hit_bit(kBStart + i, _bucket_count[i] * _bucket_scale);
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(_tag, tmp);
  output << tmp << '\n';

  if (! _function_as_tdt_filter) {
    output << "|\n";
  }

  return 1;
}

std::unique_ptr<float[]> GetRadii(Molecule& m) {
  const int matoms = m.natoms();

  std::unique_ptr<float[]> result = std::make_unique<float[]>(matoms);

  for (int i = 0; i < matoms; ++i) {
    atomic_number_t z = m.atomic_number(i);
    if (z == 6) {
      result[i] = 1.7;
    } else if (z == 7) {
      result[i] = 1.6;
    } else if (z == 8) {
      result[i] = 1.55;
    } else if (z == 9) {
      result[i] = 1.5;
    } else if (z == 15) {
      result[i] = 1.95;
    } else if (z == 16) {
      result[i] = 1.8;
    } else if (z == 17) {
      result[i] = 1.8;
    } else if (z == 35) {
      result[i] = 1.9;
    } else if (z == 53) {
      result[i] = 2.1;
    } else if (z == 1) {
      result[i] = 0.33;
    } else {
      result[i] = 2.0;
    }
  }

  return result;
}

// Return the bond scaling factor associated with the bond type of `b`.
float
BondScalingFactor(const Bond& b) {
  if (b.is_aromatic()) {
    return 0.1f;
  }

  if (b.is_single_bond()) {
    return 0.0f;
  }

  if (b.is_double_bond()) {
    return 0.2f;
  }

  if (b.is_triple_bond()) {
    return 0.3f;
  }

  cerr << "Should never come here\n";
  return 0.0f;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  std::unique_ptr<float[]> radius = GetRadii(m);

  const int matoms = m.natoms();

  std::unique_ptr<float[]> vi = std::make_unique<float[]>(matoms);

  m.compute_aromaticity_if_needed();

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    float ri = radius[a1];
    float rj = radius[a2];
    float bij = ri + rj;

    bij -= BondScalingFactor(*b);

    float dij = std::min(std::max(std::abs(ri - rj), bij), ri + rj);

    vi[a1] += ri * rj - (rj - dij) * (rj - dij) / dij;
    vi[a2] += rj * ri - (ri - dij) * (ri - dij) / dij;
  }

  return Labute(m, vi.get(), output);
}

int
Options::Abraham(Molecule& m, const float* vi, IWString_and_File_Descriptor& output) {
  if (! m.compute_Abraham_partial_charges()) {
    cerr << "Options::Abraham:compute_Abraham_partial_charges failed\n";
    return _ignore_partial_charge_calcucation_failed;
  }

  return CommonOutput(m, m.partial_charges().rawdata(), vi, _buckets[static_cast<int>(PartialChargeType::kAbraham)], output);
}

int
Options::Huckel(Molecule& m, const float* vi, IWString_and_File_Descriptor& output) {
  if (! m.compute_Huckel_partial_charges()) {
    cerr << "Options::Gasteiger:compute_Huckel_partial_charges failed\n";
    return _ignore_partial_charge_calcucation_failed;
  }

  return CommonOutput(m, vi, m.partial_charges().rawdata(), _buckets[static_cast<int>(PartialChargeType::kHuckel)], output);
}
int
Options::Gasteiger(Molecule& m, const float* vi, IWString_and_File_Descriptor& output) {
  if (! m.compute_Gasteiger_partial_charges()) {
    cerr << "Options::Gasteiger:compute_Gasteiger_partial_charges failed\n";
    return _ignore_partial_charge_calcucation_failed;
  }

  return CommonOutput(m, vi, m.partial_charges().rawdata(), _buckets[static_cast<int>(PartialChargeType::kGasteiger)], output);
}

int
Options::GasteigerHuckel(Molecule& m, const float* vi, IWString_and_File_Descriptor& output) {
  if (! m.compute_Gasteiger_Huckel_partial_charges()) {
    cerr << "Options::Gasteiger:compute_Gasteiger_Huckel_partial_charges failed\n";
    return _ignore_partial_charge_calcucation_failed;
  }

  return CommonOutput(m, vi, m.partial_charges().rawdata(), _buckets[static_cast<int>(PartialChargeType::kGasteigerHuckel)], output);
}

int
Options::Labute(Molecule& m,
                const float* vi,
                IWString_and_File_Descriptor& output) {
  int rc = 1;

  if (_compute_estate_values) {
    return EStateAtomicProperties(m, vi, output);
  } else if (_partial_charge == PartialChargeType::kUnspecified) {
    cerr << "No partial charge specified\n";
    return 0;
  } else if (_partial_charge == PartialChargeType::kAbraham) {
    return Abraham(m, vi, output);
  } else if (_partial_charge == PartialChargeType::kGasteiger) {
    return Gasteiger(m, vi, output);
  } else if (_partial_charge == PartialChargeType::kHuckel) {
    return Huckel(m, vi, output);
  } else if (_partial_charge == PartialChargeType::kGasteigerHuckel) {
    return GasteigerHuckel(m, vi, output);
  } else {
    cerr << "Options::Labute:unrecognised partial charge type\n";
    return 0;
  }

  if (rc == 0) {
    cerr << "Options::Labute:partial charge calculation failed\n";
    return 0;
  }

  return 1;
}

int
Options::EStateAtomicProperties(Molecule& m,
                                const float* vi,
                                IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  std::unique_ptr<double[]> e_state = std::make_unique<double[]>(matoms);

  lillymol_estate::determine_atom_e_state_index(m, e_state.get());
  return CommonOutput(m, vi, e_state.get(), _buckets[static_cast<int>(PartialChargeType::kEState)], output);
}


int
LabuteApproximateSurfaceArea(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
LabuteApproximateSurfaceArea(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! LabuteApproximateSurfaceArea(options, *m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4092);
  }

  return 1;
}

int
LabuteApproximateSurfaceArea(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "LabuteApproximateSurfaceArea:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return LabuteApproximateSurfaceArea(options, input, output);
}

int
LabuteApproximateSurfaceAreaRecord(Options& options, const_IWSubstring buffer,  // note local copy
                IWString_and_File_Descriptor& output) {
  buffer.remove_leading_chars(options.smiles_tag().length());
  buffer.chop();

  Molecule m;
  if (! m.build_from_smiles(buffer)) [[ unlikely ]] {
    cerr << "LabuteApproximateSurfaceAreaRecord:invalid smiles\n";
    cerr << buffer << '\n';
    return 0;
  }

  return LabuteApproximateSurfaceArea(options, m, output);
}

int
LabuteApproximateSurfaceArea(Options& options, iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
    if (! buffer.starts_with(options.smiles_tag())) {
      continue;
    }
    if (! LabuteApproximateSurfaceAreaRecord(options, buffer, output)) {
      cerr << "LabuteApproximateSurfaceArea:error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
LabuteApproximateSurfaceArea(Options& options, const char* fname,
                IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "LabuteApproximateSurfaceArea:cannot open '" << fname << "'\n";
    return 0;
  }

  return LabuteApproximateSurfaceArea(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:Q:J:fqdX:");

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

  if (cl.option_present('f')) {
    if (cl.size() > 1) {
      cerr << "TDT pipeline input only works with one input file\n";
      Usage(1);
    }
  }

  if (cl.option_present('f')) {
    if (! LabuteApproximateSurfaceArea(options, cl[0], output)) {
      cerr << "LabuteApproximateSurfaceArea::error on tdt input\n";
      return 1;
    }
  } else {
    options.WriteHeader(output);
    for (const char * fname : cl) {
      if (! LabuteApproximateSurfaceArea(options, fname, input_type, output)) {
        cerr << "LabuteApproximateSurfaceArea::fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace labute

int
main(int argc, char ** argv) {

  int rc = labute::Main(argc, argv);

  return rc;
}
