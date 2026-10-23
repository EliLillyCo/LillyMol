// Twenty-third-pass target command-line driver shape for iwdescr after the library extraction.
// This file intentionally owns all I/O, filtering, output modes, preprocessing,
// and the command-line testing harness. IWDescr owns only descriptor calculation.

#include <cstring>
#include <iostream>
#include <memory>
#include <ranges>
#include <span>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwmisc/set_or_unset.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/iwdescr_lib.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/xlogp.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/iwdescr.pb.h"
#else
#include "iwdescr.pb.h"
#endif

namespace {

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
  cerr << R"(Computes molecular descriptors.
Use the shell wrapper iwdescr.sh to establish useful defaults.
iwdescr.sh file.smi > file.w
 -O ...                 specify feature sets to compute. Enter '-O help' for info.
 -F ...                 filter output, write smiles of molecules that pass. Enter '-F help' for info.
 -G ...                 fingerprint output, enter '-F help' for info.
 -B ...                 less used and obscure options, enter '-B help' for info.
 -H ...                 hbond  assigner specification - done in iwdescr.sh
 -T ...                 tests computed values across different random smiles, enter '-T help' for info.
 -N ...                 charge assigner specification - done in iwdescr.sh
 -u <char>              value for uncomputed values.
 -S                     the TPSA computation is done for maximum compatibility with RDKit.
 -l                     reduce input molecules to largest fragment - recommended.
 -i <type>              input type, enter '-i help' for info.
 -A ...                 standard aromaticity options, enter '-A help' for info.
 -E ...                 standard element     options, enter '-A help' for info.
 -v                     verbose output
)";
  // clang-format on

  ::exit(rc);
};

class IWDescrMainOptions {
 private:
  int _verbose = 0;

  int _reduce_to_largest_fragment = 0;

  Chemical_Standardisation _chemical_standardisation;

  // Existing command-line test harness. Keep this in the executable rather
  // than in IWDescr so the embeddable descriptor calculator stays small.
  int _ntest = 0;
  uint32_t _molecules_failing_tests = 0;
  int _keep_going_after_test_failure = 0;

  // If -T WRITE=<fname> is specified.
  IWString_and_File_Descriptor _stream_for_test_failures;

  // for each descriptor, the number of test failures.
  std::unique_ptr<uint32_t[]> _descriptor_failure;

  // Expected descriptor values used only by the command-line test harness.
  // Main owns allocation/lifetime and compares these against IWDescr::Process
  // results. IWDescr should know nothing about saved expected values.
  Set_or_Unset<float>* _saved_result = nullptr;


  // Filtering is based on computed descriptor values, but it is not inherent
  // to the descriptor calculator. The main program can compute descriptors,
  // apply filters, and decide whether/how to write output.
  // From the -F option
  std::vector<Descriptor_Filter> _descriptor_filter;
  uint32_t _rejected_by_filters = 0;

  // Output formatting and output mode controls belong to the caller.
  int _output_precision = 3;
  IWString _undefined_value = ".";
  int _include_smiles_as_descriptor = 0;

  // The namexref= option allows converting feature names.
  IW_STL_Hash_Map_String _name_translation;

  // Descriptor file pipeline controls govern alternate input/output plumbing.
  int _read_descriptor_file_pipeline = 0;
  int _write_descriptor_file_pipeline = 0;

  // Descriptor ranges govern optional output scaling/post-processing. IWDescr
  // computes raw descriptor values; main owns any scaled rendering.
  IWString _descriptor_range_proto_file_name;

  // TDT/input-output tags and related filter mode belong to the command-line
  // program, not IWDescr.
  IWString _tag;
  IWString _smiles_tag = "$SMI<";
  IWString _identifier_tag = "PCN<";
  int _work_as_tdt_filter = 0;

  IWDigits _iwdigits;
  Fraction_as_String _fraction_as_string;

  // The old alarm feature has not been used for many years. Prefer omitting it
  // from the IWDescr API. If retained, it should be owned here.
  unsigned int _alarm_time = 0;
  uint32_t _molecules_skipped_by_timer = 0;

  int _flush_after_each_molecule = 0;

  FileType _input_type = FILE_TYPE_INVALID;
  uint64_t _molecules_read = 0;

  // Will be prepended to the names of descriptors in tabular output.
  IWString _prefix = "w_";

  // Integer values will be written a "1." rather than "1".
  int _always_decimal_output = 0;

  // At the end of a run, report the statistics of each feature.
  int _report_descriptor_statistics = 0;

  Report_Progress _report_progress;

  char _output_separator = ' ';

  // private functions
  int WriteResults(Molecule& m, const float* results,
             const IWDescr& iwdescr, IWString_and_File_Descriptor& output) const;
  int InitialiseFingerprints(Command_Line& cl, IWDescr& iwdescr);
  int ShowPossibleMatches(const IWString& dname,
                    const std::span<Descriptor>& descriptors);
  int WriteFingerprint(Molecule& m, IWDescr& iwdescr,
                  const float* result,
                  IWString_and_File_Descriptor& output) const;
  bool OutputAsWholeNumber(float v, IWString& output) const;
  int Preprocess(Molecule& m);
  int PerformTests(Molecule& m, IWDescr& iwdescr);
  void MaybeReportProgress();
  int results_are_same(const std::optional<float>* saved_result,
                      const float* result,
                      IWDescr& iwdescr,
                      const IWString& mname,
                      resizable_array<int>& failed_descriptors) const;

  int WriteMaybeTranslatedName(const IWString& name,
                        IWString_and_File_Descriptor& output) const;

  // Functions called by InitialiseRangesAndFilters.

  // Once the IWDescr has been instantiated, parses the -F option.
  int ReadFilterSpecifications(Command_Line& cl, IWDescr& iwdescr);
  int NameToIndex(const IWString& name, const std::span<Descriptor>& descriptors) const;

  int ReadNameTranslation(const const_IWSubstring& fname);

  // This needs to be called after Initialise and when an IWDescr object has been
  // instantiated. Parses the -B ranges=<fname> directive.
  // <fname> is stored in _descriptor_range_proto_file_name.
  // If no ranges have been specified, returns true.
  int ReadDescriptorRanges(const IWDescr& iwdescr);

  int TestIwdescriptors(const IWString* rsmi, IWDescr& iwdescr,
                        const IWString& mname);
  int ApplyFilters(Molecule& m, const float* results,
                        IWDescr& iwdescr,
                        const IWString& starting_smiles,
                        IWString_and_File_Descriptor& output);

 public:
  IWDescrMainOptions();

  int Initialise(Command_Line& cl);

  int verbose() const {
    return _verbose;
  }

  const IWString& prefix() const {
    return _prefix;
  }

  FileType input_type() const {
    return _input_type;
  }

  int read_descriptor_file_pipeline() const {
    return _read_descriptor_file_pipeline;
  }
  int write_descriptor_file_pipeline() const {
    return _write_descriptor_file_pipeline;
  }

  char output_separator() const {
    return _output_separator;
  }

  int ntest() const {
    return _ntest;
  }

  // There are certain options that can only be applied once an IWDescr object
  // has been instantiated. Call this 
  int InitialiseRangesAndFilters(Command_Line& cl, IWDescr& iwdescr);

  int WriteHeader(const IWDescr& iwdescr,
            bool include_id,
            IWString_and_File_Descriptor& output) const;
  int Process(Molecule& m, IWDescr& iwdescr, IWString_and_File_Descriptor& output);

  void MaybeFlush(IWString_and_File_Descriptor& output) const;

  int Report(const IWDescr& iwdescr, std::ostream& output) const;
  int ReportTestingResults(const IWDescr& iwdescr, std::ostream& output) const;
};

IWDescrMainOptions::IWDescrMainOptions() {
  _iwdigits.initialise(1024);
}

void
DisplayDashBOptions(int rc) {
  cerr << R"(The following -B qualifiers are recognised.
 -B quiet               turn off unclassified atom warnings from TPSA computation.
 -B ranges=<fname>      read fingerprint ranges from <fname>. Input is textproto.
 -B sep=<char>          output separator, default space. -B sep=, or -B sep=tab, or...
 -B namexref=<fname>    w.Feature textproto of name translations. See docs for info.
 -B prefix=<s>          descriptor name prefix, 'w_' by default. Use `-B prefix=none` for no prefix.
 -B dstats              report individual statistics for each descriptor.
 -B flush               flush output after each molecule processed.
 -B float               integer values written as floats, "10" becomes "10."
 -B wpipe               include smiles as first column of output.
 -B rpipe               input is a descriptor file with a leading smiles column.
)";

  ::exit(rc);
}

static void
DisplayFilterOptions(std::ostream& output, int rc) {
  output << R"(
Use the -F option to filter molecules. The syntax uses fortran like comparison operators.
-F w_natoms.lt.50
only writes molecules that have fewer than 50 heavy atoms. Note that when using property filters
output is a smiles file.

Multiple filters are compiled as and conditions,
iwdescr ... -F w_natoms.gt.12 -F w_natoms.lt.50 -F w_nrings.gt.1 -F w_nrings.lt.5 file.smi > filtered.smi
would only write molecules that satisfy all conditions.
)";

  ::exit(rc);
}

static void
DisplayDashTOptions(std::ostream& output) {
  output << R"(The -T option performs tests.
For each molecule read, generate a number of random smiles and check that the computed values
are the same across all random variants.
 -T kg              if a molecule fails, do NOT terminate the compuation, keep going.
 -T rpt=<n>         report progress every <n> molecules tested.
 -T <n>             for each molecule, generate <n> random smiles variants.
 -T WRITE=<fname>   write molecules that fail to <fname>.
)";

  exit(0);
}


int
IWDescrMainOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('u')) {
    int i = 0;
    IWString u;
    while (cl.value('u', u, i++)) {
      _undefined_value = u;
    }
    if (_undefined_value == "none") {
      _undefined_value = "";
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose, 'g')) {
      return 0;
    }
  }

  // -B prefix= is descriptor metadata and is parsed by IWDescr. Main owns
  // output-specific -B options.
  if (cl.option_present('B')) {
    IWString b;
    for (int i = 0; cl.value('B', b, i); ++i) {
      if (b.starts_with("ranges=")) {
        b.remove_leading_chars(7);
        _descriptor_range_proto_file_name = b;
      } else if (b.starts_with("sep=")) {
        b.remove_leading_chars(4);
        IWString tmp(b);
        if (! char_name_to_char(tmp)) {
          cerr << "Unrecognised character specification '" << b << "'\n";
          return 0;
        }
        _output_separator = tmp[0];
      } else if (b == "flush") {
        _flush_after_each_molecule = 1;
      } else if (b == "rpipe") {
        _read_descriptor_file_pipeline = 1;
      } else if (b == "wpipe") {
        _write_descriptor_file_pipeline = 1;
      } else if (b == "quiet") {
        // intercepted in the library.
      } else if (b.starts_with("prefix=")) {
        b.remove_leading_chars(7);
        if (b == "none") {
          _prefix.resize(0);
        } else {
          _prefix = b;
        }
        if (_verbose) {
          cerr << "Descriptors generated with prefix '" << _prefix << "'\n";
        }
      } else if (b.starts_with("namexref=")) {
        b.remove_leading_chars(9);
        if (! ReadNameTranslation(b)) {
          cerr << "Cannot read name translation data from '" << b << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << _name_translation.size() << " name translations defined\n";
        }
      } else if (b.starts_with("mxsfsdf=") ||
                 b.starts_with("namexref=")) {
        // Handled by IWDescr or deliberately deferred in this driver.
        continue;
      } else if (b == "dstats") {
        _report_descriptor_statistics = 1;
        if (_verbose) {
          cerr << "Will report individual descriptor statistics\n";
        }
      } else if (b == "float") {
        _iwdigits.append_to_each_stored_string(".");
        if (_verbose) {
          cerr << "Whole numbers written as floats\n";
        }
      } else if (b == "help") {
        DisplayDashBOptions(0);
      } else {
        std::cerr << "Unrecognised -B qualifier '" << b << "'\n";
        DisplayDashBOptions(1);
        return 0;
      }
    }
  }

  if (cl.option_present('s')) {
    _include_smiles_as_descriptor = 1;
  }

  // Ran into problems trying to do test with a charge assigner present

  if (cl.option_present('T')) {
    if (cl.option_present('N')) {
      cerr << "Test mode does not work with a charge assigner\n";
      return 0;
    }

    const_IWSubstring t;
    for (int i = 0; cl.value('T', t, i); ++i) {
      if ("kg" == t) {
        _keep_going_after_test_failure = 1;
      } else if (t.starts_with("rpt=")) {
        t.remove_leading_chars(4);
        uint32_t rpt;
        if (! t.numeric_value(rpt)) {
          cerr << "Invalid report progress specification '" << t << "'\n";
          return 0;
        }
        _report_progress.set_report_every(rpt);
        if (_verbose) {
          cerr << "Will report progress every " << rpt << " molecules processed\n";
        }
      } else if (t.starts_with("WRITE=")) {
        t.remove_leading_chars(6);
        IWString fname(t);
        fname.EnsureEndsWith(".smi");
        if (! _stream_for_test_failures.open(fname)) {
          cerr << "Cannot open stream for test failures '" << fname << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Molecules failing tests written to '" << fname << "'\n";
        }
      } else if (t == "help") {
        DisplayDashTOptions(cerr);
      } else if (!t.numeric_value(_ntest) || _ntest < 1) {
        cerr << "The test option (-T) must be followed by a whole positive number '" <<
                t << "' invalid\n";
        return 0;
      }
    }

    if (0 == _ntest) {
      cerr << "Must specify the number of test to perform with the -T option\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will perform tests on " << _ntest << " random smiles permutations\n";
    }
  }

  if (_read_descriptor_file_pipeline) {
    ;  // No Molecule input, reading text.
  } else if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      std::cerr << "cannot discern input type\n";
      return 0;
    }
  } else if (cl.size() == 1 && std::strcmp(cl[0], "-") == 0) {
    _input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    std::cerr << "Cannot auto-detect file types\n";
    return 0;
  }

  if (!_fraction_as_string.initialise(0.0f, 1.0f, _output_precision)) {
    cerr << "Cannot initialise fractional descriptor output cache\n";
    return 0;
  }

  return 1;
}

int
IWDescrMainOptions::InitialiseRangesAndFilters(Command_Line& cl, IWDescr& iwdescr) {
  if (! ReadFilterSpecifications(cl, iwdescr)) {
    cerr << "Cannot initialise filter specifications\n";
    return 0;
  }

  if (! ReadDescriptorRanges(iwdescr)) {
    cerr << "Cannot initialise ranges\n";
    return 0;
  }

  if (! InitialiseFingerprints(cl, iwdescr)) {
    cerr << "Cannot initialise fingerprints\n";
    return 0;
  }

  if (_ntest > 0) {
    const int n = iwdescr.number_descriptors();
    _descriptor_failure = std::make_unique<uint32_t[]>(n);
  }

  return 1;
}

static void
DisplayFingerprintOptions(std::ostream& output) {
  // clang-format off
  output << R"(The following fingerprint directives are recognised.
 -G FILTER         work as a TDT filter, used internally by gfp_make.
 -G R=<n>          the number of buckets used in discretising the values.
 -G ALL            use all descriptors to generate the fingerprint.
 -G BEST           from calibration runs, certain descriptors have been designated.
                   as the 'best'. Use these designated to generate the fingerprint.
 -G d1,d2,d3...    specify individual descriptors to be fingerprinted.
 <n>               generate <n> replicates for each bit.

 If any descriptor name is followed by a ':n', that feature will include 'n'.
 replicates of that bit in the output. By default, all features get the same number.
 -G natoms:10,nrings:5
 -G natoms -G rotbond: 2
 The order in which options are specified is important.
 natoms and rotbond get 10 replicates, alogp gets 20.
 -G 10 -G natoms,rotbond,alogp:20
 natoms gets 10 replicates
 -G natoms:10 -G 20     # '-G 20' has no effect.
 natoms gets 10, rotbond 20 and nrings gets 3 replicates.
 -G 10 -G natoms -G 20 -G rotbond -G nrings:3
)";
}

static int
parse_replicates_specification(const_IWSubstring& dname, int& replicates) {
  if (!dname.contains(':')) {
    return 1;
  }

  const_IWSubstring tmp;
  if (!dname.split_into_directive_and_value(tmp, ':', replicates)) {
    return 0;
  }

  dname.truncate_at_first(':');

  return 1;
}

// Someone has made an attempt to look for `dname` in `descriptors`
// but no match was found. Try to explain what might be going on.
// Always return 0;
int
IWDescrMainOptions::ShowPossibleMatches(const IWString& dname,
                    const std::span<Descriptor>& descriptors) {
  IWString mydname(dname);
  if (_prefix.empty()) {
  } else if (dname.starts_with(_prefix)) {
    mydname.remove_leading_chars(_prefix.length());
  }

  cerr << "Looking for '" << mydname << '\n';
  int inactive = 0;
#ifdef LILLYMOL_HAS_CPP23
  for (auto [i, d] : std::views::enumerate(descriptors)) {
#else
  for (uint32_t i = 0; i < descriptors.size(); ++i) {
    const Descriptor& d = descriptors[i];
#endif
    const IWString& x = d.descriptor_name();
    if (! d.active()) {
      cerr << "Feature " << i << " inactive\n";
      ++inactive;
      continue;
    }

    cerr << x << '\n';

    if (x.contains(mydname)) {
      cerr << "Descriptor " << i << " name '" << x << "' contains '" << mydname <<
              "' active? " << d.active() << '\n';
    }
  }

  cerr << "No fingerprint match to '" << dname << "'\n";
  cerr << "Found " << inactive << " inactive descriptors\n";
  cerr << "If you have specified something like '-O none -G mxdst' that will fail.\n";
  cerr << "The first directive, '-O none' turns off all expensive desctiptors, including distance matrix\n";
  cerr << "Requesting a fingerprint based on a descriptor not being generated is impossible\n";

  return 0;
}

int
IWDescrMainOptions::InitialiseFingerprints(Command_Line& cl, IWDescr& iwdescr) {
  if (! cl.option_present('G')) {
    return 1;
  }

  std::span<Descriptor> descriptors = iwdescr.Descriptors();

  int replicates = 1;
  int resolution = 10;

  const_IWSubstring s;
  // Replicates and resolution must be discerned before anything else.
  for (int i = 0; cl.value('G', s, i); ++i) {

    if (s.starts_with("R=")) {
      s.remove_leading_chars(2);
      if (!s.numeric_value(resolution) || resolution < 2) {
        cerr << "The resolution on fingerprints (R=) must be a whole +ve number\n";
        return 1;
      }
      continue;
    }

    int tmp;
    if (s.numeric_value(tmp) && tmp > 0) {
      replicates = tmp;
      continue;
    }
  }

  for (int i = 0; cl.value('G', s, i); ++i) {
    if (s == "help") {
      DisplayFingerprintOptions(cerr);
      return 0;
    }

    if ("FILTER" == s) {
      _work_as_tdt_filter = 1;
      continue;
    }

    if (s.starts_with("TAG=")) {
      _tag = s;
      _tag.remove_leading_chars(4);
      _tag.EnsureEndsWith('<');
      continue;
    }

    if (s.starts_with("ALL")) {
      int tmpr = replicates;
      if (!parse_replicates_specification(s, tmpr)) {
        cerr << "Invalid ALL,replicates specification '" << s << "'\n";
        return 2;
      }

      for (Descriptor& d : descriptors) {
        d.set_produce_fingerprint(tmpr);
      }
      continue;
    }

    if (s.starts_with("BEST")) {  // leave open possibility of adding a number later
      int tmpr = replicates;
      if (!parse_replicates_specification(s, tmpr)) {
        cerr << "Invalid BEST,replicates specification '" << s << "'\n";
        return 2;
      }

      for (Descriptor& d : descriptors) {
        if (d.best_fingerprint()) {
          d.set_produce_fingerprint(tmpr);
        }
      }
      continue;
    }

    const_IWSubstring dname;
    for (int j = 0; s.nextword(dname, j, ',');) {
      int tmpr = replicates;
      if (!parse_replicates_specification(dname, tmpr)) {
        cerr << "Invalid descriptor:replicate specification '" << dname << "'\n";
        return 0;
      }

      int k = DescriptorNumber(_prefix, dname, descriptors);

      if (k < 0) {
        return ShowPossibleMatches(dname, descriptors);
      }

//    cerr << "descriptor " << descriptors[k].descriptor_name() << " getting " << tmpr
//         << " replicates " << replicates << '\n';
      descriptors[k].set_produce_fingerprint(tmpr);
    }
  }

  if (_tag.empty()) {
    _tag = "NCIWD<";
  }

  FillDescriptorExtremeties(descriptors, resolution);

  return 1;
}

int
IWDescrMainOptions::ReadFilterSpecifications(Command_Line& cl, IWDescr& iwdescr) {
  const int n = cl.option_count('F');
  if (n == 0) {
    return 1;
  }

  std::span<const Descriptor> descriptor_data(iwdescr.descriptor_data(),
                iwdescr.number_descriptors());

  _descriptor_filter.resize(n);
  for (int i = 0; i < n; ++i) {
    const_IWSubstring f = cl.string_value('F', i);
    if (f == "help") {
      DisplayFilterOptions(cerr, 0);
      continue;
    }

    if (! _descriptor_filter[i].build(_prefix, f, descriptor_data)) {
      cerr << "Cannot initialise descriptor filter '" << f << "'\n";
      return 0;
    }
  }

  if (_verbose) {
    cerr << "Defined " << _descriptor_filter.size() << " descriptor filters\n";
  }

  return 1;
}

int
IWDescrMainOptions::ReadNameTranslation(const const_IWSubstring& fname) {
  IWString tmp(fname);
  std::optional<w::Features> maybe_proto = iwmisc::ReadTextProto<w::Features>(tmp);
  if (!maybe_proto) {
    cerr << "ReadNameTranslation:cannot read '" << fname << "'\n";
    return 0;
  }

  for (const auto& feature : (*maybe_proto).feature()) {
    const IWString old_name = feature.computed_name();
    const IWString new_name = feature.name();
    _name_translation[old_name] = new_name;
  }

  return _name_translation.size();
}

int
IWDescrMainOptions::ReadDescriptorRanges(const IWDescr& iwdescr) {
  if (_descriptor_range_proto_file_name.empty()) {
    return 1;
  }

  std::optional<w::Ranges> maybe_proto =
    iwmisc::ReadTextProto<w::Ranges>(_descriptor_range_proto_file_name);
  if (!maybe_proto) {
    cerr << "ReadDescriptorRanges:cannot read '" << _descriptor_range_proto_file_name << "'\n";
    return 0;
  }

  const std::span<Descriptor> descriptors = iwdescr.Descriptors();

  // Ignore missing descriptors, those may not be turned on. Too complicated
  // otherwise. But what this means is that truly bad input will be silently
  // ignored. What we need is a hash containing all known descriptor names,
  // but that does not exist.
  int rc = 0;
  for (const auto& range : (*maybe_proto).range()) {
    int ndx = DescriptorNumber(_prefix, range.name(), descriptors);
    if (ndx < 0) {
      cerr << "ReadDescriptorRanges:no match for '" << range.name() << "'\n";
      continue;
    }

    descriptors[ndx].set_range(range.min(), range.max());
    ++rc;
  }

  return rc;
}


// We have our own special purpose preprocessing
int
IWDescrMainOptions::Preprocess(Molecule& m) {
  // we want to be able to compute a result for Hydrogen "molecule"
  if (1 == m.natoms() && 1 == m.atomic_number(0)) [[unlikely]] { 
    return 1;
  }

  // Revert any hydrogen isotopes. Various other functions may be
  // reluctant to remove them otherwise.
  // Neutralize any formal charge on a Hydrogen atom.
  // https://dot-jira.lilly.com/browse/GC3TK-652?jql=assignee%20in%20(RX87690)

  int explicit_hydrogen_found = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];

    if (1 != a.atomic_number()) {
      continue;
    }

    ++explicit_hydrogen_found;

    if (0 != a.isotope()) {
      m.set_isotope(i, 0);
    }

    if (a.formal_charge()) {
      m.set_formal_charge(i, 0);
    }
  }

  if (_chemical_standardisation.active()) {
    (void)_chemical_standardisation.process(m);
  }

  if (m.number_fragments() > 1) {
    if (0 == _reduce_to_largest_fragment) {
      cerr << "Fatal, '" << m.name() << " has " << m.number_fragments()
           << " components\n";
      // iwabort();
    }

    (void)m.reduce_to_largest_fragment_carefully();
  }

  // The charge assigner calls move_hydrogens_to_end_of_connection_table,
  // which messes up cached atomic properties. Make that call now to
  // hopefully make the subsequent call a no-op.
  if (explicit_hydrogen_found) {
    m.MoveToEndOfConnectionTable(1);
  }

  const IWString& mname = m.name();

  if (mname.empty()) {
    IWString newname;
    newname << "IWD" << _molecules_read;

    cerr << "No name, set to '" << newname << "'\n";

    m.set_name(newname);
  } else if (_descriptor_filter.size() > 0)
    ;
  else {
    IWString newname(mname);
    newname.truncate_at_first(' ');

    m.set_name(newname);
  }

  return 1;
}

void
IWDescrMainOptions::MaybeFlush(IWString_and_File_Descriptor& output) const {
  if (_flush_after_each_molecule) {
    output.flush();
  }

  output.write_if_buffer_holds_more_than(4096);
}

// Two test values have been reported as different. Are they close enough?

int
values_are_equivalent(const Set_or_Unset<float>& s1, const Set_or_Unset<float>& s2) {
  float v1, v2;

  if (!s1.value(v1) || !s2.value(v2)) {  // they should both be set
    return 0;
  }

  if (std::isnan(v1) && std::isnan(v2)) {
    return 1;
  }

  if (fabs(v1 - v2) < 1.0e-05) {
    return 1;
  }

  return 0;
}

// The indices of failed desctiptors go in `failed_descriptors`.
int
IWDescrMainOptions::results_are_same(const std::optional<float>* saved_result,
                      const float* result,
                      IWDescr& iwdescr,
                      const IWString& mname,
                      resizable_array<int>& failed_descriptors) const {
  failed_descriptors.resize_keep_storage(0);

  const int n = iwdescr.number_descriptors();

  for (int i = 0; i < n; i++) {
    if (! saved_result[i]) {  // not active
      continue;
    }

    if (saved_result[i] == result[i]) {
      continue;
    }

    if (values_are_equivalent(*saved_result[i], result[i])) {
      continue;
    }

    if (failed_descriptors.empty()) {  // one of these warnings per molecule.
      cerr << "Mismatch '" << mname << "'\n";
    }

    ++_descriptor_failure[i];

    const IWString& descriptor_name = iwdescr.descriptor_name(i);
    cerr << "Result mismatch, descriptor " << i << " '" << descriptor_name << "'\n";
    cerr << *saved_result[i] << " vs " << result[i] << " diff " <<
                (*saved_result[i] - result[i]) << '\n';
    failed_descriptors << i;
  }

  // Return true if there are no failed comparisons.
  return failed_descriptors.empty();
}

void
IWDescrMainOptions::MaybeReportProgress() {
  if (_report_progress()) {
    cerr << "Tested " << _molecules_read << " molecules, " <<
         _molecules_failing_tests << " failed testing\n";
  }
}

int
IWDescrMainOptions::TestIwdescriptors(const IWString* rsmi, IWDescr& iwdescr,
                const IWString& mname) {
  MaybeReportProgress();

  const int n = iwdescr.number_descriptors();

  // Save the results
  std::unique_ptr<std::optional<float>[]> saved_result =
                std::make_unique<std::optional<float>[]>(n);

  static std::unique_ptr<float[]> result = std::make_unique<float[]>(n);

  IWString failure_reasons;
  IWString failed_smiles;
  resizable_array<int> failed_descriptors;

  for (int i = 0; i < _ntest; i++) {
    Molecule tmp;
    if (!tmp.build_from_smiles(rsmi[i])) {
      cerr << "Yipes, cannot parse random smiles '" << rsmi[i] << "'\n";
      return 0;
    }

    if (! iwdescr.Process(tmp, result.get())) [[unlikely]] {
      cerr << "test_iwdescriptors::computation failed\n";
      cerr << rsmi[i] << ' ' << mname << '\n';
      return 0;
    }

    if (i == 0) {  // Store the first results.
      for (int j = 0; j < n; ++j) {
        if (iwdescr.descriptor_active(j)) {
          saved_result[j] = result[j];
        }
      }
      continue;
    }

    if (results_are_same(saved_result.get(), result.get(), iwdescr, mname, failed_descriptors)) {
      continue;
    }

    cerr << "Yipes, one or more result mismatches\n";

    for (int r : failed_descriptors) {
      failure_reasons.append_with_spacer(iwdescr.descriptor_name(r), ' ');
    }

    cerr << rsmi[i] << ' ' << mname << " rsmi test failure\n";

    if (failed_smiles.empty()) {
      failed_smiles = rsmi[i];
    }
  }

  // If no failures, done.
  if (failed_smiles.empty()) {
    return 1;
  }

  ++_molecules_failing_tests;

  if (! _stream_for_test_failures.is_open()) {
    return 0;
  }

  _stream_for_test_failures << failed_smiles << ' ' << mname << ' ' << failure_reasons << '\n';
  _stream_for_test_failures.flush();

  return 0;
}

//  Doing tests is complicated by the fact that the underlying functions change
//  the molecule. Generate all the unique smiles variants ahead of time so
//  all tests start with the unchanged molecule.
int
IWDescrMainOptions::PerformTests(Molecule& m, IWDescr& iwdescr) {
  std::unique_ptr<IWString[]> rsmi = std::make_unique<IWString[]>(_ntest);
  // cerr << "Parent smiles " << m.smiles() << '\n';
  for (int i = 0; i < _ntest; i++) {
    rsmi[i] = m.random_smiles();
    // cerr << rsmi[i] << ' ' << m.name() << '\n';
  }

  if (TestIwdescriptors(rsmi.get(), iwdescr, m.name())) {
    return 1;
  }

  cerr << "One or more test failures '" << m.name() << "' molecule "
       << _molecules_read << '\n';
  cerr << m.smiles() << '\n';

  if (_keep_going_after_test_failure) {
    return 1;
  }

  return 0;
}

int
IWDescrMainOptions::Report(const IWDescr& iwdescr, std::ostream& output) const {
  output << _molecules_read << " molecules read\n";
  if (_ntest > 0) {
    output << "Ran " << _ntest << " tests per molecule " <<
          _molecules_failing_tests << " molecules failed\n";
  }

  for (const Descriptor_Filter& f : _descriptor_filter) {
    f.Report(output);
  }

  if (_report_descriptor_statistics) {
    iwdescr.ReportDescriptorStatistics(output);
  }

  return 1;
}

int
IWDescrMainOptions::ReportTestingResults(const IWDescr& iwdescr,
                std::ostream& output) const {
  output << "Tested " << _molecules_read << " countered " << _molecules_failing_tests <<
            " molecules failing tests\n";
  const std::span<Descriptor> descriptors = iwdescr.Descriptors();

  int failing_descriptors = 0;
  for (uint32_t i = 0; i < descriptors.size(); ++i) {
    if (_descriptor_failure[i] == 0) {
      continue;
    }

    ++failing_descriptors;
    output << descriptors[i].name() << ' ' << _descriptor_failure[i] << " failures\n";
  }

  output << failing_descriptors << " descriptors with failed computations\n";

  return 1;
}

int
IWDescrMainOptions::WriteHeader(const IWDescr& iwdescr,
            bool include_id,
            IWString_and_File_Descriptor& output) const {
  if (_ntest > 0) {
    return 1;
  }

  if (_tag.length() > 0) {  // fingerprints, no header
    return 1;
  } 

  if (_descriptor_filter.size() > 0 || _work_as_tdt_filter) {
    return 1;
  }

  if (_read_descriptor_file_pipeline) {
  } else if (_write_descriptor_file_pipeline) {
    output << "smiles";
    if (include_id) {
      output << _output_separator << "Id";
    }
  } else {
    output << "Name";
  }

  if (_include_smiles_as_descriptor) {
    output << _output_separator << "smiles";
  }

  const int n = iwdescr.number_descriptors();
  for (int i = 0; i < n; ++i) {
    if (! iwdescr.descriptor_active(i)) {
      continue;
    }

    output << _output_separator;
    if (_name_translation.empty()) {
      output << _prefix << iwdescr.descriptor_name(i);
    } else {
      WriteMaybeTranslatedName(iwdescr.descriptor_name(i), output);
    }
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  if (_verbose) {
    std::cerr << "Output will contain " << n << " descriptors\n";
  }

  return output.good();
}

// this is complicated by the need to accommodate different forms
// of names specified.
// for example we want to accommodate both 'natoms' and "${prefix}natoms"
int
IWDescrMainOptions::WriteMaybeTranslatedName(const IWString& name,
                        IWString_and_File_Descriptor& output) const {
  if (_prefix.size() > 0) {
    output << _prefix;
  }

  auto iter = _name_translation.find(name);
  if (iter != _name_translation.end()) {
    output << iter->second;
    return 1;
  }

  // If we have a prefix, we can look for "${prefix}name"
  if (_prefix.empty()) {
    output << name;
  }

  IWString tmp(_prefix);
  tmp << name;
  iter = _name_translation.find(tmp);
  if (iter == _name_translation.end()) {
    output << name;
  } else {
    output << iter->second;
  }

  return 1;
}

int
IWDescrMainOptions::ApplyFilters(Molecule& m, const float* results,
                        IWDescr& iwdescr,
                        const IWString& starting_smiles,
                        IWString_and_File_Descriptor& output) {
  for (Descriptor_Filter& d : _descriptor_filter) {
    if (d.satisfied(iwdescr.descriptor_data())) {
      continue;
    }

    ++_rejected_by_filters;
    return 1;
  }

  output << starting_smiles << ' ' << m.name() << '\n';

  MaybeFlush(output);

  return 1;
}

int
IWDescrMainOptions::WriteFingerprint(Molecule& m, IWDescr& iwdescr,
                  const float* result,
                  IWString_and_File_Descriptor& output) const {
  if (!_work_as_tdt_filter) {  // reading a smiles file
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator sfc;

  const Descriptor* descriptor = iwdescr.descriptor_data();

  int bstart = 0;

  const int n = iwdescr.number_descriptors();
  for (int i = 0; i < n; i++) {
    if (! iwdescr.descriptor_active(i)) {
      continue;
    }

    if (!descriptor[i].produce_fingerprint()) {
      continue;
    }

    descriptor[i].produce_fingerprint(bstart, sfc);

    bstart += descriptor[i].bit_replicates();
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(_tag, tmp);
  output << tmp << '\n';

  if (!_work_as_tdt_filter) {
    output << "|\n";
  }

  MaybeFlush(output);

  return 1;
}

// If `v` is close to being a whole number, write as an integer to `output`.
bool
IWDescrMainOptions::OutputAsWholeNumber(float v, IWString& output) const {
  if (fabs(static_cast<int>(v) - v) > 1.0e-05) {
    return false;
  }

  _iwdigits.append_number(output, static_cast<int>(v));

  return true;
}

int
IWDescrMainOptions::WriteResults(Molecule& m, const float* results,
             const IWDescr& iwdescr, IWString_and_File_Descriptor& output) const {
  if (_read_descriptor_file_pipeline && _write_descriptor_file_pipeline) {
    const_IWSubstring tmp(m.name());
    tmp.remove_leading_words(1);
    output << tmp;   // includes all previously calculated descriptors.
  } else if (_write_descriptor_file_pipeline) {
    output << m.smiles() << _output_separator;
    append_first_token_of_name(m.name(), output);
  } else {
    append_first_token_of_name(m.name(), output);
  }

  if (_include_smiles_as_descriptor) {
    output << _output_separator << m.smiles();
  }

  const int n = iwdescr.number_descriptors();
  for (int i = 0; i < n; ++i) {
    if (! iwdescr.descriptor_active(i)) [[unlikely]] {
      continue;
    }

    output << _output_separator;

    const float v = results[i];

    if (std::isnan(v)) {
      output << _undefined_value;
    } else if (OutputAsWholeNumber(v, output)) {
      // good
    } else if (v > 100.0f) {
      output.append_number(v, 5);
    } else {
      _fraction_as_string.append_number(output, v);
    }
  }

  output << '\n';

  if (_report_descriptor_statistics) {
    for (Descriptor& d : iwdescr.Descriptors()) {
      if (d.active()) {
        d.update_statistics(_verbose);
      }
    }
  }

  return output.good();
}

int
IWDescrMainOptions::Process(Molecule& m,
                            IWDescr& iwdescr,
                            IWString_and_File_Descriptor& output) {
  ++_molecules_read;
  if (_ntest > 0) {
    Preprocess(m);
    return PerformTests(m, iwdescr);
  }

  std::unique_ptr<IWString> starting_smiles;
  if (_descriptor_filter.size() > 0) {
    starting_smiles = std::make_unique<IWString>(m.smiles());
  }

  // Or should this have been called before starting_smiles.
  // Here seems better.
  Preprocess(m);

  static std::unique_ptr<float[]> results =
    std::make_unique<float[]>(iwdescr.number_descriptors());

  if (! iwdescr.Process(m, results.get())) {
    cerr << m.smiles() << ' ' << m.name() << " computation failed\n";
    return 0;
  }

  if (_descriptor_filter.size() > 0) {
    return ApplyFilters(m, results.get(), iwdescr, *starting_smiles, output);
  }

  if (! _tag.empty()) {
    return WriteFingerprint(m, iwdescr, results.get(), output);
  }

  return WriteResults(m, results.get(), iwdescr, output);
}

int
ProcessMolecules(data_source_and_type<Molecule>& input,
                 IWDescr& iwdescr,
                 IWDescrMainOptions& options,
                 IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Process(*m, iwdescr, output)) {
      return 0;
    }

    options.MaybeFlush(output);
  }

  return 1;
}

// Note that we might change `buffer`.
static int
iwdescriptors_rpipe_line(const_IWSubstring& buffer, IWDescr& iwdescr,
                         IWDescrMainOptions& options,
                         IWString_and_File_Descriptor& output) {
  const_IWSubstring smiles;
  int i = 0;
  if (! buffer.nextword(smiles, i) || smiles.empty()) {
    cerr << "No first token\n";
    return 0;
  }

  Molecule m;
  // smiles is first token on line and all other tokens become the name.
  if (!m.build_from_smiles(smiles)) {
    cerr << "iwdescriptors_rpipe_line:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  const_IWSubstring id;
  buffer.nextword(id, i);
  m.set_name(id);
  
  if (options.write_descriptor_file_pipeline()) {
    output << buffer;
  } else {
    buffer.remove_leading_words(1);
    output << buffer;
  }
  output << options.output_separator();

  return options.Process(m, iwdescr, output);
}

// We are reading pipelined input.
// ID natoms nrings...
// SMILES ID 1 2 3
int
ProcessPipelinedDescriptorFile(iwstring_data_source& input,
            IWDescr& iwdescr,
            IWDescrMainOptions& options,
            IWString_and_File_Descriptor& output) {
  assert(options.read_descriptor_file_pipeline());

  static bool first_call = true;
  const_IWSubstring buffer;

  if (first_call) {
    // Process the header record.
    if (!input.next_record(buffer)) {
      cerr << "ProcessPipelinedDescriptorFile:cannot read header\n";
      return 0;
    }

    if (options.write_descriptor_file_pipeline()) {
      output << buffer << options.output_separator();
    }
    options.WriteHeader(iwdescr, false, output);

    first_call = false;
  }

  while (input.next_record(buffer)) {
    if (!iwdescriptors_rpipe_line(buffer, iwdescr, options, output)) {
      cerr << "Fatal error\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
ProcessPipelinedDescriptorFile(const char* fname,
            IWDescr& iwdescr,
            IWDescrMainOptions& options,
            IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "iwdescriptors_rpipe:cannot open '" << fname << "'\n";
    return 0;
  }

  return ProcessPipelinedDescriptorFile(input, iwdescr, options, output);
}

int
ProcessFile(const char* fname,
            IWDescr& iwdescr,
            IWDescrMainOptions& options,
            IWString_and_File_Descriptor& output) {
  if (options.read_descriptor_file_pipeline()) {
    return ProcessPipelinedDescriptorFile(fname, iwdescr, options, output);
  }

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.ok()) {
    std::cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ProcessMolecules(input, iwdescr, options, output);
}

int
iwdescr_main(int argc, char** argv) {
  Command_Line cl(argc, argv, "g:N:u:A:fq:E:vi:lH:b:T:O:F:a:G:s:SB:d");
  if (cl.unrecognised_options_encountered()) {
    std::cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  // Global options stay here: aromaticity, elements, input/output setup.

  if (!cl.option_present('A')) {
    set_global_aromaticity_type(Daylight);
  } else if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(1);
  }

  if (cl.option_present('E')) {
    process_elements(cl);
  }

  IWDescrMainOptions main_options;
  if (! main_options.Initialise(cl)) {
    Usage(1);
  }

  IWDescr iwdescr;

  if (! iwdescr.Initialise(cl)) {
    std::cerr << "Cannot initialise descriptor computation\n";
    return 1;
  }

  if (! main_options.InitialiseRangesAndFilters(cl, iwdescr)) {
    return 1;
  }

  int output_precision = 4;
  set_default_iwstring_float_concatenation_precision(output_precision);

  if (cl.empty()) {
    std::cerr << "no files specified\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (!main_options.read_descriptor_file_pipeline()) {
    if (! main_options.WriteHeader(iwdescr, true, output)) {
      return 1;
    }
  }

  for (const char* fname : cl) {
    if (! ProcessFile(fname, iwdescr, main_options, output)) {
      std::cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (cl.option_present('T')) {
    main_options.ReportTestingResults(iwdescr, cerr);
  }

  if (verbose) {
    main_options.Report(iwdescr, cerr);
  }

  return 0;
}

}  // namespace

int
main(int argc, char** argv) {
  return iwdescr_main(argc, argv);
}
