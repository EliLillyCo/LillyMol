#include <stdlib.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "maccskeys_fn5.h"

namespace maccskeys_fn5_main {

using std::cerr;

class IntAccumulator : public Accumulator_Int<int> {
 private:
  int _number_zero = 0;

 public:
  int extra(int value);

  int
  number_zero() const {
    return _number_zero;
  }

  double average_number_of_hits() const;
};

int
IntAccumulator::extra(int value) {
  if (value == 0) {
    ++_number_zero;
  } else {
    Accumulator_Int<int>::extra(value);
  }

  return 1;
}

// The default average() method returns the average number of times the bit was
// hit - including the influence from non-matches. For this method, we exclude
// the influence of non-matches.
double
IntAccumulator::average_number_of_hits() const {
  assert(_number_zero >= 0);

  if (_number_zero == 0) {
    return Accumulator_Int<int>::average();
  }

  const int number_non_zero = Accumulator_Int<int>::n() - _number_zero;
  if (number_non_zero == 0) {
    return 0.0;
  }

  return Accumulator_Int<int>::sum() / static_cast<double>(number_non_zero);
}

void
DisplayDashXOptions(std::ostream& output) {
  output << R"(
 -X flush         flush output after each molecule
)";

  ::exit(0);
}

class Options {
 private:
  const char* _prog_name = nullptr;

  int _number_maccs_keys = 192;

  // Always strip to largest fragment, the library uses several bonds_between calls.
  bool _strip_to_largest_fragment = true;
  bool _revert_all_directional_bonds_to_non_directional = false;
  bool _remove_all_chiral_centres = false;

  int _nbits = 256;

  IWString _fingerprint_tag;
  IWString _level2_fingerprint_tag;
  IWString _non_colliding_fingerprint_tag;
  IWString _fixed_size_counted_fingerprint_tag;

  bool _append_nset = false;
  IWDigits _iwdigits;

  int _verbose = 0;
  uint64_t _molecules_read = 0;

  int _ntest = 0;
  bool _report_test_failures = true;
  Report_Progress _report_progress;
  bool _keep_going_after_test_failure = false;

  uint64_t _molecules_with_test_failures = 0;
  uint64_t _test_failures = 0;

  Chemical_Standardisation _chemical_standardisation;

  bool _flush_after_each_molecule = false;

  std::unique_ptr<IntAccumulator[]> _accumulators;
  extending_resizable_array<int> _keys_hit;

  bool _write_key_values = true;
  bool _write_header = true;

  MACCSKeys _mk;

  void Usage(int rc) const;
  bool WriteHeader(IWString_and_File_Descriptor& output) const;
  void MaybeFlush(IWString_and_File_Descriptor& output) const;
  bool DoOutput(const Molecule& m, const int* keys, IWString_and_File_Descriptor& output);
  bool AllKeysTheSame(const int* mk1, const int* mk2,
                      resizable_array<int>& problematic_keys) const;
  int CountAromaticRings(Molecule& m) const;
  bool DoComparisonBetweenKeyValues(Molecule& m1, const IWString& initial_smiles,
                                    Molecule& m2, const IWString& random_smiles,
                                    const int* k1, const int* k2, int permutation,
                                    resizable_array<int>& problematic_keys);
  bool TestMaccsKeys(Molecule& m, const int* keys);
  bool ComputeKeys(Molecule& m, int* keys);
  bool WriteFixedSizeCountedFingerprint(const int* keys, IWString& output) const;
  bool WriteNonCollidingFingerprint(const int* keys, IWString& output) const;
  bool WriteFingerprints(int* keys, IWString& output);
  void Preprocess(Molecule& m);
  bool Process(Molecule& m, IWString_and_File_Descriptor& output);
  bool ProcessFilter(iwstring_data_source& input, int* keys,
                     IWString_and_File_Descriptor& output);
  bool Process(data_source_and_type<Molecule>& input,
               IWString_and_File_Descriptor& output);
  bool ProcessFilter(const char* fname, IWString_and_File_Descriptor& output);
  bool Process(const char* fname, FileType input_type,
               IWString_and_File_Descriptor& output);
  void Report(std::ostream& output) const;
  void WriteStatistics(IWString_and_File_Descriptor& stream_for_statistics) const;

 public:
  int Main(int argc, char** argv);
};

void
Options::Usage(int rc) const {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << "Usage: " << _prog_name << " <options> <file1> <file2> ...\n";
  cerr << R"(MACCSKeys implementation. Started as a faithful implementation, but 
some key definitions changed and new keys added.
  -t <number>    generate <number> random smiles for testing
  -t rpt=<r>     report test progress every <r> molecules
  -k             keep going after test failure
  -S <fname>     gather and report statistics on key distributions
  -n             no output of key values (only recognised with -S or -J)
  -Y <fname>     write keys in Daylight form to <fname>
  -J <tag>       tag for fingerprint data
  -J LEVEL2=XXX  use XXX as the level 2 tag (features with multiple occurrences)
  -h             exclude hydrogens on aromatic nitrogen atoms
  -b <nbits>     number of bits to produce (default 192)
  -e             include the number of bits set as a descriptor
  -x             remove cis trans bonds
  -f             work as a filter
  -i <type>      specify input type
  -X ...         miscellaneous options, enter '-X help' for info
  -l             strip to the largest fragment
  -c             remove chiral centres
  -E <el>        process element <el>
  -g ...         chemcial standardisation options, enter '-g help' for info.
  -v             verbose output
)";

  exit(rc);
}

bool
Options::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "Name";

  for (int i = 0; i < _number_maccs_keys; ++i) {
    output << " mk_mk" << i;
  }

  if (_append_nset) {
    output << " mknset";
  }

  output << '\n';

  return output.good();
}

void
Options::MaybeFlush(IWString_and_File_Descriptor& output) const {
  if (_flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(4096);
  }
}

bool
Options::DoOutput(const Molecule& m, const int* keys,
                  IWString_and_File_Descriptor& output) {
  append_first_token_of_name(m.name(), output);

  int nset = 0;

  IWString& s = output;

  for (int i = 0; i < _number_maccs_keys; ++i) {
    if (keys[i] > 0) {
      ++nset;
    }

    _iwdigits.append_number(s, keys[i]);
  }

  if (_append_nset) {
    _iwdigits.append_number(s, nset);
  }

  output += '\n';

  return output.good();
}

bool
Options::AllKeysTheSame(const int* mk1, const int* mk2,
                        resizable_array<int>& problematic_keys) const {
  bool result = true;
  for (int i = 0; i < _number_maccs_keys; ++i) {
    if (mk1[i] != mk2[i]) {
      if (_report_test_failures) {
        cerr << "Fingerprint mismatch, bit " << i << " values " << mk1[i] << " vs "
             << mk2[i] << '\n';
      }
      problematic_keys.add_if_not_already_present(i);
      result = false;
    }
  }

  return result;
}

int
Options::CountAromaticRings(Molecule& m) const {
  const int nr = m.nrings();

  int result = 0;

  for (int i = 0; i < nr; ++i) {
    if (m.ringi(i)->is_aromatic()) {
      ++result;
    }
  }

  return result;
}

bool
Options::DoComparisonBetweenKeyValues(Molecule& m1, const IWString& initial_smiles,
                                      Molecule& m2, const IWString& random_smiles,
                                      const int* k1, const int* k2, int permutation,
                                      resizable_array<int>& problematic_keys) {
  if (AllKeysTheSame(k1, k2, problematic_keys)) {
    return true;
  }

  ++_test_failures;

  if (!_report_test_failures) {
    return false;
  }

  cerr << "Yipes, key mismatch on permutation " << permutation << '\n';
  cerr << "Molecule '" << m1.name() << "'\n";
  cerr << initial_smiles << '\n';
  cerr << random_smiles << '\n';

  // If these two forms have different unique smiles, that may indicate
  // problems with aromaticity/sssr determination.
  const IWString& m1usmi = m1.unique_smiles();
  const IWString& m2usmi = m2.unique_smiles();
  if (m1usmi == m2usmi) {
    cerr << "Unique_smiles match\n";
  } else {
    cerr << "Unique smiles mismatch\n";
    cerr << m1usmi << '\n';
    cerr << m2usmi << '\n';
    return true;
  }

  const int nr = m1.nrings();
  if (nr != m2.nrings()) {
    cerr << "Ring count mismatch!, " << nr << " vs " << m2.nrings() << '\n';
    return true;
  }

  const int a1 = CountAromaticRings(m1);
  const int a2 = CountAromaticRings(m2);

  if (a1 != a2) {
    cerr << "Aromatic ring count differs, " << a1 << " vs " << a2 << " rings\n";
    return true;
  }

  cerr << "Both forms have " << a1 << " aromatic rings\n";

  return false;
}

bool
Options::TestMaccsKeys(Molecule& m, const int* keys) {
  int tkeys[NUMBER_MACCS_KEYS + 1];

  IWString initial_smiles = m.smiles();

  resizable_array<int> problematic_keys;

  int failures_this_molecule = 0;

  for (int i = 0; i < _ntest; ++i) {
    IWString random_smiles = m.random_smiles();
    Molecule mtmp;
    if (!mtmp.build_from_smiles(random_smiles)) {
      cerr << "Yipes, cannot build from smiles '" << random_smiles << "'\n";
      return false;
    }

    assert(mtmp.natoms() == m.natoms());
    assert(mtmp.nrings() == m.nrings());

    if (!_mk(mtmp, tkeys)) {
      cerr << "Yipes, cannot form maccs keys for variant " << i << '\n';
      cerr << "Smiles '" << random_smiles << "'\n";
      return _keep_going_after_test_failure;
    }

    if (DoComparisonBetweenKeyValues(m, initial_smiles, mtmp, random_smiles, keys, tkeys,
                                     i, problematic_keys)) {
      continue;
    }

    ++failures_this_molecule;

    if (!_keep_going_after_test_failure) {
      ++_molecules_with_test_failures;
      return false;
    }
    break;
  }

  // Need to check Kekule forms. Do it randomly.
  if (m.contains_aromatic_atoms()) {
    for (int i = 0; i < _ntest; ++i) {
      m.compute_aromaticity_if_needed();
      set_include_aromaticity_in_smiles(1);
      IWString random_smiles = m.random_smiles();
      set_include_aromaticity_in_smiles(0);

      Molecule mtmp;
      if (!mtmp.build_from_smiles(random_smiles)) {
        cerr << "Cannot parse smiles '" << random_smiles << "', '" << m.name() << "'\n";
        continue;
      }

      if (!_mk(mtmp, tkeys)) {
        cerr << "Yipes, cannot form maccs keys for variant " << i << '\n';
        cerr << "Smiles '" << random_smiles << "'\n";
        return _keep_going_after_test_failure;
      }

      if (DoComparisonBetweenKeyValues(m, initial_smiles, mtmp, random_smiles, keys,
                                       tkeys, i, problematic_keys)) {
        continue;
      }

      ++failures_this_molecule;
      if (!_keep_going_after_test_failure) {
        ++_molecules_with_test_failures;
        return false;
      }
      break;
    }
  }

  if (failures_this_molecule) {
    ++_molecules_with_test_failures;
    cerr << m.smiles() << ' ' << m.name() << " FAILED\n";
    for (int i = 0; i < problematic_keys.number_elements(); ++i) {
      cerr << " kfail " << problematic_keys[i] << '\n';
    }
  }

  if (_report_progress()) {
    cerr << "Tested " << _molecules_read << " molecules, "
         << _molecules_with_test_failures << " molecules with failures\n";
  }

  return true;
}

bool
Options::ComputeKeys(Molecule& m, int* keys) {
  bool result = _mk(m, keys);

  if (result && _ntest) {
    result = TestMaccsKeys(m, keys);
  }

  if (!result) {
    return false;
  }

  if (_accumulators) {
    int non_zero_keys = 0;
    for (int i = 0; i < _number_maccs_keys; ++i) {
      _accumulators[i].extra(keys[i]);
      if (keys[i]) {
        ++non_zero_keys;
      }
    }

    _keys_hit[non_zero_keys]++;
  }

  return true;
}

bool
Options::WriteFixedSizeCountedFingerprint(const int* keys, IWString& output) const {
  output << _fixed_size_counted_fingerprint_tag;
  (void)append_fixed_size_counted_fingerprint(keys, _nbits, -1, -1, output);
  output << ">\n";

  return output.good();
}

bool
Options::WriteNonCollidingFingerprint(const int* keys, IWString& output) const {
  IWString dyascii;
  non_colliding_counted_fingerprint_daylight_representation(_nbits, keys, dyascii);

  output << _non_colliding_fingerprint_tag << dyascii << ">\n";

  return output.good();
}

bool
Options::WriteFingerprints(int* keys, IWString& output) {
  if (_fingerprint_tag.length()) {
    IW_Bits_Base fp;

    (void)fp.construct_from_array_of_ints(keys, _nbits);

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output << _fingerprint_tag << tmp << ">\n";
  }

  if (_level2_fingerprint_tag.length()) {
    _mk.set_level_2_fingerprint(keys);

    IW_Bits_Base fp(_nbits);

    fp.construct_from_array_of_ints(keys, _nbits);

    IWString tmp;
    fp.daylight_ascii_representation_including_nset_info(tmp);
    output << _level2_fingerprint_tag << tmp << ">\n";
  }

  if (_non_colliding_fingerprint_tag.length()) {
    WriteNonCollidingFingerprint(keys, output);
  }

  if (_fixed_size_counted_fingerprint_tag.length()) {
    WriteFixedSizeCountedFingerprint(keys, output);
  }

  return output.good();
}

void
Options::Preprocess(Molecule& m) {
  m.remove_all(1);  // explicit hydrogens mess things up

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_strip_to_largest_fragment && m.number_fragments() > 1) {
    m.reduce_to_largest_fragment();
  }

  if (_remove_all_chiral_centres) {
    m.remove_all_chiral_centres();
  }

  if (_revert_all_directional_bonds_to_non_directional) {
    m.revert_all_directional_bonds_to_non_directional();
  }

  // Ensure computed - to avoid incremental updates.
  m.recompute_distance_matrix();
}

bool
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  Preprocess(m);

  int keys[256];  // make long enough for fingerprints however long they are

  std::fill_n(keys, 256, 0);

  const bool result = ComputeKeys(m, keys);

  if (_ntest) {
    return result;
  }

  if (!result) {
    return false;
  }

  if (_write_key_values) {
    if (!DoOutput(m, keys, output)) {
      return false;
    }
  }

  if (_fingerprint_tag.length() || _level2_fingerprint_tag.length() ||
      _non_colliding_fingerprint_tag.length() ||
      _fixed_size_counted_fingerprint_tag.length()) {
    output << "$SMI<" << m.smiles() << ">\n";
    output << "PCN<" << m.name() << ">\n";
    WriteFingerprints(keys, output);

    if (m.number_records_text_info()) {  // let's hope it is in TDT format!
      m.write_extra_text_info(output);
    }

    output << "|\n";
  }

  return true;
}

bool
Options::ProcessFilter(iwstring_data_source& input, int* keys,
                       IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with("$SMI<")) {
      buffer.remove_leading_chars(5);
      buffer.chop();

      Molecule m;
      if (!m.build_from_smiles(buffer)) {
        cerr << "Yipes! cannot parse smiles '" << buffer << "'\n";
        return false;
      }

      ++_molecules_read;
      Preprocess(m);

      if (!ComputeKeys(m, keys)) {
        cerr << "Cannot compute keys for molecule at line " << input.lines_read() << '\n';
        cerr << buffer << '\n';
        return false;
      }

      output << "$SMI<" << m.smiles() << ">\n";
      WriteFingerprints(keys, output);
    } else if (buffer == '|') {
      output << "|\n";

      std::fill_n(keys, _nbits, 0);
    } else {
      output << buffer << '\n';
    }

    MaybeFlush(output);
  }

  return true;
}

bool
Options::Process(data_source_and_type<Molecule>& input,
                 IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    ++_molecules_read;

    if (_verbose > 1) {
      cerr << _molecules_read << ' ' << m->name() << '\n';
    }

    if (!Process(*m, output)) {
      return false;
    }

    MaybeFlush(output);
  }

  return true;
}

bool
Options::ProcessFilter(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open stream file '" << fname << "'\n";
    return false;
  }

  int tmp[256];
  std::fill_n(tmp, 256, 0);

  return ProcessFilter(input, tmp, output);
}

bool
Options::Process(const char* fname, FileType input_type,
                 IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open input '" << fname << "'\n";
    return false;
  }

  return Process(input, output);
}

void
Options::Report(std::ostream& output) const {
  output << _molecules_read << " molecules read\n";

  if (_test_failures) {
    output << _molecules_with_test_failures << " molecules had " << _test_failures
           << " test failures\n";
  }
}

void
Options::WriteStatistics(IWString_and_File_Descriptor& stream_for_statistics) const {
  assert(_accumulators != nullptr);

  stream_for_statistics << "Key statistics";
  if (_verbose == 0) {
    stream_for_statistics << " for " << _molecules_read << " molecules";
  }
  stream_for_statistics << '\n';

  for (int i = 0; i < _number_maccs_keys; ++i) {
    const IntAccumulator& a = _accumulators[i];
    stream_for_statistics << "Key " << i << ' ' << a.n() << " molecules ";
    if (a.number_zero()) {
      stream_for_statistics << a.number_zero() << " misses";
    } else {
      stream_for_statistics << "min hits " << a.minval();
    }

    stream_for_statistics << " max hits = " << a.maxval();
    if (a.n() > 1) {
      stream_for_statistics << " ave = " << a.average();
    }

    stream_for_statistics << '\n';
  }

  Accumulator_Int<int> keys_hit_per_molecule;
  for (int i = 0; i < _keys_hit.number_elements(); ++i) {
    if (_keys_hit[i]) {
      stream_for_statistics << _keys_hit[i] << " molecules had " << i << " keys hit\n";
      keys_hit_per_molecule.extra(i, _keys_hit[i]);
    }
  }

  if (keys_hit_per_molecule.n() > 1) {
    stream_for_statistics << keys_hit_per_molecule.average()
                          << " average bits hit per molecule\n";
  }

  stream_for_statistics << "static int level2_threshold [] = {\n";
  for (auto i = 0; i < _number_maccs_keys; ++i) {
    if (_accumulators[i].n() == 0) {
      stream_for_statistics << "  0";
    } else {
      stream_for_statistics << "  "
                            << static_cast<int>(_accumulators[i].average() + 0.4999999);
    }

    if (_number_maccs_keys - 1 != i) {
      stream_for_statistics << ',';
    }

    stream_for_statistics << "     // key " << i << "   ";
    if (_accumulators[i].n() == 0) {
      stream_for_statistics << ".\n";
    } else {
       stream_for_statistics << static_cast<float>(_accumulators[i].average()) << "\n";
    }
  }

  stream_for_statistics << "};\n";
}

int
Options::Main(int argc, char** argv) {
  _prog_name = argv[0];

  Command_Line cl(argc, argv, "Y:b:nS:K:kt:J:fFvi:A:Clg:E:ehxcX:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(3);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, _verbose, 'E')) {
      cerr << "Cannot parse -E option\n";
      Usage(18);
    }
  }

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, _verbose)) {
      Usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', _number_maccs_keys) || _number_maccs_keys <= 0) {
      cerr << "The number of maccs keys to produce (-b) must be a valid +ve number\n";
      Usage(2);
    }

    if (!_mk.set_nbits(_number_maccs_keys)) {
      cerr << "Invalid input for number MACCS keys " << _number_maccs_keys << '\n';
      return 2;
    }

    if (_verbose) {
      cerr << "Will produce " << _number_maccs_keys << " bits\n";
    }
  }

  _iwdigits.set_include_leading_space(1);
  _iwdigits.initialise(255);

  if (cl.empty()) {
    cerr << "Must specify files to process on the command line\n";
    Usage(5);
  }

  if (cl.option_present('t')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++)) {
      if (t.starts_with("rpt=")) {
        int rp;
        t.remove_leading_chars(4);
        if (!t.numeric_value(rp) || rp < 1) {
          cerr << "The report progress -t qualifier must be a whole +ve number\n";
          Usage(4);
        }

        _report_progress.set_report_every(rp);
      } else if ("norpt" == t) {
        _report_test_failures = false;
        if (_verbose) {
          cerr << "Will not report test failures\n";
        }
      } else if (t.numeric_value(_ntest) && _ntest > 0) {
        if (_verbose) {
          cerr << "Will perform " << _ntest
               << " random permutation tests on each molecule\n";
        }
      } else {
        cerr << "Unrecognised -t qualifier '" << t << "'\n";
        Usage(5);
      }
    }

    _write_header = false;
    set_input_aromatic_structures(1);
  }

  if (cl.option_present('k')) {
    if (!cl.option_present('t')) {
      cerr << "The -k option only makes sense with the -t option\n";
      Usage(11);
    }

    _keep_going_after_test_failure = true;
    if (_verbose) {
      cerr << "Will continue execution even after a test failure\n";
    }
  }

  IWString_and_File_Descriptor stream_for_statistics;

  if (cl.option_present('S')) {
    const char* s = cl.option_value('S');

    if (!stream_for_statistics.open(s)) {
      cerr << "Cannot open stream for statistics '" << s << "'\n";
      return 5;
    }

    if (_verbose) {
      cerr << "Statistics on key frequency will be written to '" << s << "'\n";
    }

    _accumulators = std::make_unique<IntAccumulator[]>(_number_maccs_keys + 1);
  }

  if (cl.option_present('n')) {
    if (_ntest == 0 && !cl.option_present('J') && !stream_for_statistics.is_open()) {
      cerr << "The -n option specified, but no other output\n";
      cerr << "Specify either -S or -J for some other output\n";
      Usage(83);
    }

    _write_key_values = false;
    _write_header = false;
    if (_verbose) {
      cerr << "Writing keys to cout suppressed\n";
    }
  }

  if (cl.option_present('h')) {
    _mk.set_aromatic_nitrogens_do_not_have_hydrogens(1);
    if (_verbose) {
      cerr << "Will not consider Hydrogen atoms on aromatic Nitrogens\n";
    }
  }

  if (cl.option_present('x')) {
    _revert_all_directional_bonds_to_non_directional = true;

    if (_verbose) {
      cerr << "Will remove directional bonds\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_all_chiral_centres = true;

    if (_verbose) {
      cerr << "Will remove all chiral centres\n";
    }
  }

  if (cl.option_present('C')) {
    cerr << "The -C option is no longer recognised\n";
  }

  if (cl.option_present('F')) {
    cerr << "The -F option is no longer recognised\n";
  }

  if (cl.option_present('K')) {
    cerr << "The -K option is no longer recognised\n";
  }

  if (cl.option_present('e')) {
    _append_nset = true;

    if (_verbose) {
      cerr << "The number of bits set will be included as a descriptor\n";
    }
  }

  if (cl.option_present('l')) {
    _strip_to_largest_fragment = true;
    if (_verbose) {
      cerr << "Will strip input molecules to largest fragment\n";
    }
  }

  if (cl.option_present('J')) {
    if (!cl.option_present('n')) {
      _write_key_values = false;
      _write_header = false;
      if (_verbose) {
        cerr << "Fingerprint output, array output suppressed\n";
      }
    }

    int i = 0;
    const_IWSubstring j;
    while (cl.value('J', j, i++)) {
      if (j.starts_with("LEVEL2=")) {
        _level2_fingerprint_tag = j;
        _level2_fingerprint_tag.remove_leading_chars(7);

        if (_verbose) {
          cerr << "Level 2 fingerprint tag '" << _level2_fingerprint_tag << "'\n";
        }

        _level2_fingerprint_tag.EnsureEndsWith('<');
      } else if (j.starts_with("NC=")) {
        _non_colliding_fingerprint_tag = j;
        _non_colliding_fingerprint_tag.remove_up_to_first('=');

        if (_verbose) {
          cerr << "Non colliding fingerprint tag '" << _non_colliding_fingerprint_tag
               << "'\n";
        }

        _non_colliding_fingerprint_tag.EnsureEndsWith('<');
      } else if (j.starts_with("FC=")) {
        _fixed_size_counted_fingerprint_tag = j;
        _fixed_size_counted_fingerprint_tag.remove_leading_chars(3);

        if (_verbose) {
          cerr << "Fixed size counted fingerprint tag "
               << _fixed_size_counted_fingerprint_tag << "'\n";
        }

        _fixed_size_counted_fingerprint_tag.EnsureEndsWith('<');
      } else {
        _fingerprint_tag = j;
        if (_verbose) {
          cerr << "Fingerprints written with dataitem '" << _fingerprint_tag << "'\n";
        }

        _fingerprint_tag.EnsureEndsWith('<');
      }
    }

    if (_level2_fingerprint_tag.empty() && _fingerprint_tag.empty()) {
      ;
    } else if (_level2_fingerprint_tag == _fingerprint_tag) {
      cerr << "Both default and level 1 fingerprint tags identical. Impossible\n";
      Usage(11);
    }

    _nbits = _number_maccs_keys;
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        _flush_after_each_molecule = true;
        if (_verbose) {
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

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  // If we are a filter, we take a different path.
  if (cl.option_present('f')) {
    if (_fingerprint_tag.empty() && _non_colliding_fingerprint_tag.empty()) {
      _fingerprint_tag = "FPMK<";
      cerr << "Fingerprint tag '" << _fingerprint_tag << "'\n";
    }

    for (int i = 0; i < cl.number_elements(); ++i) {
      if (!ProcessFilter(cl[i], output)) {
        cerr << "Error processing '" << cl[i] << "'\n";
        return 87;
      }
    }
  } else {
    if (_write_header) {
      if (!WriteHeader(output)) {
        return 8;
      }
    }

    FileType input_type = FILE_TYPE_INVALID;
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }

    for (int i = 0; i < cl.number_elements(); ++i) {
      if (!Process(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  if (_verbose) {
    Report(cerr);
  }

  if (stream_for_statistics.is_open()) {
    WriteStatistics(stream_for_statistics);
    stream_for_statistics.close();
  }

  output.flush();

  return rc;
}

}  // namespace maccskeys_fn5_main

int
main(int argc, char** argv) {
  maccskeys_fn5_main::Options options;
  return options.Main(argc, argv);
}
