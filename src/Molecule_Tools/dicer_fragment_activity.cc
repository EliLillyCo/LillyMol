// Read output from dicer and associate the fragments with activity values.

#include <iostream>
#include <optional>
#include <string>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Tools/dicer_fragments.pb.h"
#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"

namespace dicer_fragment_activity {

using iw_tf_data_record::TFDataReader;
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
  cerr << R"(Associates dicer fragments with activity.
First run dicer, with either the '-B textproto' or '-B serialized_proto' options. This
will contain the diced molecules.
Using serialized_proto is dramatically faster.
This tool reads the fragment data and an activity file and for each fragment, reports the
average activity of the molecules containing the fragment.

A typical use might be

dicer -k 3 -B nbamide -B brcb -I 1 -B serialized_proto -S diced.dat -v file.smi
dicer_fragment_activity -A file.activity -tfdata -support 10 diced.dat > file.txt

The resulting tabular file can be imported into data visualisation environments.

The output might look like

Frag Parent Natoms N fraction Minval Ave Maxval Diff
smiles1 name1 11 59 0.0004928207 0 28.32882 87.1 -12.68904
smiles2 name1 13 22 0.0001837637 0 14.43561 49.6 -26.58226
smiles3 name1 6 49 0.0004092918 0 65.80041 100 24.78254
smiles4 name1 16 25 0.0002088223 0 47.8084 96.6 6.790536

Frag: the unique smiles of the fragment
Parent: the name of the first molecule containing this fragment
Natoms: number of atoms in the fragment
N: number of molecules containing this fragment
Fraction: N expressed as a fraction
Minval: minimum activity of molecules containing this fragment
Ave: average activity of molecules containing this fragment
Maxval: max activity of molecules containing this fragment
Diff: the difference between Ave and the average for all molecules.

The most interesting fragments are those with high values, positive or negative
for Diff. Hopefully also with a large N value, but that will seldom be the case.

 -A <fname>             activity file (id and activity, with a header)
 -tfdata                dicer output generated with '-B serialized_proto'
 -support <n>           discard fragments unles they have at least <n> examples
 -v                     verbose output.
)";

  // clang-format on

  ::exit(rc);
}

// Information about each fragment encountered.
// We do not record the smiles since this will be kept in a hash
// and the key will be the smiles.
class Fragment {
 private:
  // An accumulation of the activity values associated with this fragment.
  Accumulator<double> _activity;

  // The molecule in which this was first observed.
  std::string _parent;

  // Number of atoms in fragment - copied from proto.
  uint32_t _natoms;

 public:
  void
  set_parent(const std::string& s) {
    _parent = s;
  }

  void
  set_natoms(uint32_t s) {
    _natoms = s;
  }

  uint32_t
  times_encountered() const {
    return _activity.n();
  }

  void
  Extra(double a) {
    _activity.extra(a);
  }

  int WriteResult(const std::string& smiles, uint32_t nmolecules, double mean_activity,
                  IWString_and_File_Descriptor& output) const;
};

int
Fragment::WriteResult(const std::string& smiles, uint32_t nmolecules,
                      double mean_activity, IWString_and_File_Descriptor& output) const {
  static constexpr char kSep = ' ';

  float ave = _activity.average();

  output << smiles << kSep << _parent << kSep << _natoms << kSep << _activity.n() << kSep
         << iwmisc::Fraction<float>(_activity.n(), nmolecules) << kSep
         << _activity.minval() << kSep << ave << kSep << _activity.maxval() << kSep
         << static_cast<float>(ave - mean_activity) << '\n';

  return 1;
}

class Options {
 private:
  int _verbose;

  int _input_is_textproto = 1;

  // Only report fragments that occur this many times or more.
  uint32_t _support = 0;
  uint32_t _suppressed_by_support_requirement = 0;

  char _input_separator;

  // this is derived from the activity values read. We do not
  // determine it ourself, it is passed via set_mean_activity.
  double _mean_activity;

 public:
  Options();

  int Initialise(Command_Line_v2& cl);

  int ReadActivity(const char* fname, absl::flat_hash_map<std::string, double>& activity);
  int ReadActivity(iwstring_data_source& input,
                   absl::flat_hash_map<std::string, double>& activity);
  int ReadACtivityRecord(const const_IWSubstring& buffer,
                         absl::flat_hash_map<std::string, double>& activity);

  int
  input_is_textproto() const {
    return _input_is_textproto;
  }

  void
  set_mean_activity(double s) {
    _mean_activity = s;
  }

  int WriteResult(const std::string& smiles, uint32_t nmolecules,
                  const Fragment& fragment, IWString_and_File_Descriptor& output);

  int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _input_is_textproto = 1;
  _input_separator = ' ';
  _support = 0;
  _suppressed_by_support_requirement = 0;
  _mean_activity = 0.0;
}

int
Options::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present("tfdata")) {
    _input_is_textproto = 0;
    if (_verbose) {
      cerr << "Will read input as TFDataRecord\n";
    }
  }

  if (cl.option_present("support")) {
    cl.value("support", _support);
    if (_verbose) {
      cerr << "Will only write fragments that are found in " << _support
           << " or more molecules\n";
    }
  }

  return 1;
}

int
Options::ReadActivity(const char* fname,
                      absl::flat_hash_map<std::string, double>& activity) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Options:ReadActivityAsTextProto:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadActivity(input, activity);
}

int
Options::ReadActivity(iwstring_data_source& input,
                      absl::flat_hash_map<std::string, double>& activity) {
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Options::ReadACtivity:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (!ReadACtivityRecord(buffer, activity)) {
      cerr << "Options::ReadACtivity:invalid input\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return activity.size();
}

int
Options::ReadACtivityRecord(const const_IWSubstring& buffer,
                            absl::flat_hash_map<std::string, double>& activity) {
  const_IWSubstring token;
  int i = 0;
  if (!buffer.nextword(token, i, _input_separator)) {
    cerr << "Options::ReadACtivityRecord:cannot get id\n";
    return 0;
  }

  std::string id(token.AsString());

  if (!buffer.nextword(token, i, _input_separator)) {
    cerr << "Options::ReadACtivityRecord:cannot get activity\n";
    return 0;
  }

  double a;
  if (!token.numeric_value(a)) {
    cerr << "Options::ReadACtivityRecord:invalid activity value '" << token << "'\n";
    return 0;
  }

  activity[id] = a;

  return 1;
}

// for a given fragment, while smiles is `smiles`, write the results.
// In order to know prevalence of a fragment, the total number of
// molecules is needed, `nmolecules`.
int
Options::WriteResult(const std::string& smiles, uint32_t nmolecules,
                     const Fragment& fragment, IWString_and_File_Descriptor& output) {
  if (_support > 0 && fragment.times_encountered() < _support) {
    ++_suppressed_by_support_requirement;
    return 1;
  }

  return fragment.WriteResult(smiles, nmolecules, _mean_activity, output);
}

int
Options::Report(std::ostream& output) const {
  if (_support > 0) {
    output << _suppressed_by_support_requirement << " fragments suppressed by "
           << _support << " support requirement\n";
  }

  return 1;
}

int
AddToFragments(const dicer_data::DicedMolecule& parent,
               const dicer_data::DicerFragment& fragment, double activity,
               absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  auto iter = fragment_data.find(fragment.smi());
  if (iter == fragment_data.end()) {
    Fragment f;
    f.set_parent(parent.name());
    f.set_natoms(fragment.nat());
    f.Extra(activity);
    fragment_data.emplace(fragment.smi(), std::move(f));
  } else {
    iter->second.Extra(activity);
  }

  return 1;
}

int
ReadDicerData(const dicer_data::DicedMolecule& proto,
              const absl::flat_hash_map<std::string, double>& activity,
              absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  const auto iter = activity.find(proto.name());
  if (iter == activity.end()) {
    cerr << "ReadDicerData:no activity for '" << proto.name() << "'\n";
    return 0;
  }

  for (const auto& fragment : proto.fragment()) {
    AddToFragments(proto, fragment, iter->second, fragment_data);
  }

  return 1;
}

int
ReadDicerDataTfDataRecord(TFDataReader& reader, Options& options,
                          const absl::flat_hash_map<std::string, double>& activity,
                          absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  while (1) {
    std::optional<dicer_data::DicedMolecule> maybe_proto =
        reader.ReadProto<dicer_data::DicedMolecule>();

    if (!maybe_proto) {
      return 1;
    }

    if (!ReadDicerData(*maybe_proto, activity, fragment_data)) {
      cerr << "ReadDicerDataTfDataRecord:invalid proto\n";
      cerr << maybe_proto->ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
ReadDicerDataTfDataRecord(const char* fname, Options& options,
                          const absl::flat_hash_map<std::string, double>& activity,
                          absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  TFDataReader reader(fname);
  if (!reader.good()) {
    cerr << "ReadDicerDataTfDataRecord:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadDicerDataTfDataRecord(reader, options, activity, fragment_data);
}

int
ReadDicerDataTextprotoRecord(const const_IWSubstring& buffer,
                             const absl::flat_hash_map<std::string, double>& activity,
                             absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  dicer_data::DicedMolecule proto;
  const absl::string_view tmp(buffer.data(), buffer.length());
  if (!google::protobuf::TextFormat::ParseFromString(tmp, &proto)) {
    cerr << "ReadDicerDataTextprotoRecord:cannot parse text proto " << buffer << '\n';
    return 0;
  }

  return ReadDicerData(proto, activity, fragment_data);
}

int
ReadDicerDataTextproto(iwstring_data_source& input, Options& options,
                       const absl::flat_hash_map<std::string, double>& activity,
                       absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!ReadDicerDataTextprotoRecord(buffer, activity, fragment_data)) {
      cerr << "ReadDicerDataTextproto:invalid textproto record\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
ReadDicerDataTextproto(const char* fname, Options& options,
                       const absl::flat_hash_map<std::string, double>& activity,
                       absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "ReadDicerDataTextproto:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadDicerDataTextproto(input, options, activity, fragment_data);
}

int
ReadDicerData(const char* fname, Options& options,
              const absl::flat_hash_map<std::string, double>& activity,
              absl::flat_hash_map<std::string, Fragment>& fragment_data) {
  if (options.input_is_textproto()) {
    return ReadDicerDataTextproto(fname, options, activity, fragment_data);
  }

  return ReadDicerDataTfDataRecord(fname, options, activity, fragment_data);
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-A=sfile-tfdata-support=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!cl.option_present('A')) {
    cerr << "Must specify activity file with the -A option\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  absl::flat_hash_map<std::string, double> activity;
  double mean_activity = 0.0;

  if (cl.option_present('A')) {
    IWString fname = cl.string_value('A');
    if (!options.ReadActivity(fname.null_terminated_chars(), activity)) {
      cerr << "Cannot read activity '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Read " << activity.size() << " id->activity values from '" << fname
           << "'\n";
    }

    Accumulator<double> acc;
    for (const auto& [_, v] : activity) {
      acc.extra(v);
    }
    mean_activity = acc.average();
    if (verbose) {
      cerr << "Values btw " << acc.minval() << " and " << acc.maxval() << " ave "
           << acc.average() << '\n';
    }
  }

  options.set_mean_activity(mean_activity);

  // Mapping from unique smiles to data about the fragment.
  absl::flat_hash_map<std::string, Fragment> fragment_data;

  for (const char* fname : cl) {
    if (!ReadDicerData(fname, options, activity, fragment_data)) {
      cerr << "Error reading dicer data '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read information on " << fragment_data.size() << " dicer fragments\n";
  }

  IWString_and_File_Descriptor output(1);
  output.reserve(8192);

  static constexpr char kSep = ' ';
  output << "Frag" << kSep << "Parent" << kSep << "Natoms" << kSep << "N" << kSep
         << "Fraction" << kSep << "Minval" << kSep << "Ave" << kSep << "Maxval" << kSep
         << "Diff" << '\n';

  const uint32_t nmolecules = activity.size();
  for (const auto& [smiles, fragment] : fragment_data) {
    options.WriteResult(smiles, nmolecules, fragment, output);

    output.write_if_buffer_holds_more_than(4096);
  }
  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_fragment_activity

int
main(int argc, char** argv) {
  int rc = dicer_fragment_activity::Main(argc, argv);

  return rc;
}
