// Tool to aid duplicate detection.
// Reads molecules and groups by molecular formula.
// Molecules with the same molecular formula are written to the
// same file.

#include <iostream>
#include <limits>

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace formula_groupings {

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
  cerr << R"(Subdivides potentially larger numbers of molecules into files contianing molecules
with the same molecular formulae. This can facilitate uniqueness determinations or other groupings.
 -F <fmla>      The kind of formula used.
                MF (default) regular molecular formula
                AROMF aromatic distinguishing molecular formula - more precise, will create more files.
                Aromatic distinguishing formula is much more expensive to compute, and will create many more
                files. Chembl with default formula takes 18 seconds to process Chembl, generating 322k files.
                Using AROMF it takes 50 seconds and creates 1.502M files.
 -S <fname>  name stem for output files created - note there can be a LOT of files created.
 -v          verbose output
)";
// clang-format on

  ::exit(rc);
}

class Formula {
  private:
    // For each molecule in this formula group, the smiles and id.
    resizable_array_p<IWString> _line;

    // The sum of all the items in `_lines`.
    uint64_t _size;

    // The file to which these molecules will be written.
    IWString _fname;

    // If this is open, we can write molecules directly to the output stream,
    // but if not, we accumulate them in _line - above.
    std::unique_ptr<IWString_and_File_Descriptor> _output;

    uint64_t _molecules_written;

  public:
    Formula();

    void set_file_name(const IWString& s) {
      _fname = s;
    }

    uint64_t size() const {
      return _size;
    }

    uint64_t molecules_written() const {
      return _molecules_written;
    }

    bool IsOpen() const {
      return _output != nullptr;
    }

    int Extra(IWString&& s);

    uint64_t Write();

    int Close();
};

Formula::Formula() {
  _size = 0;
  _molecules_written = 0;
}

// If the file is open, write immediately.
// Otherwise cache.
int
Formula::Extra(IWString&& s) {
  if (_output) {
    *_output << s << '\n';
    _output->write_if_buffer_holds_more_than(4096);
    ++_molecules_written;
    return 1;
  }

  _line << new IWString(s);
  _size += _line.size();

  return 1;
}

uint64_t
Formula::Write() {
  assert (! _output);

  _output = std::make_unique<IWString_and_File_Descriptor>();

  IWString fname;
  if (_molecules_written >= 0) {
    fname << ">>";
  }
  fname << _fname;

  if (! _output->open(fname.null_terminated_chars())) {
    cerr << "Formula::Extra:cannot open " << fname << "'\n";
    return 0;
  }

  uint64_t rc = 0;
  for (const IWString* s : _line) {
    *_output << *s << '\n';
    _output->write_if_buffer_holds_more_than(4096);
    ++rc;
  }

  _molecules_written += _line.size();

  _line.resize_keep_storage(0);
  _size = 0;

  return rc;
}

int
Formula::Close() {
  Write();

  _output->close();

  _output.reset(nullptr);

  return 1;
}

enum class MFType {
  kUnspecified,
  kMf,
  kAromMf
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

    MFType _mftype;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    absl::flat_hash_map<IWString, Formula> _formula;

    uint64_t _size;
    uint64_t _write_threshold;

    IWString _name_stem;

    uint64_t _files_open;
    uint32_t _max_files_open;

  // Private functions

    int MaybeWrite();

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

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

    int Flush();

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _mftype = formula_groupings::MFType::kMf;
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;

  _size = 0;
  // 4GB RAM
  _write_threshold = std::numeric_limits<uint32_t>::max();

  _files_open = 0;
  _max_files_open = 512;
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

  if (cl.option_present('F')) {
    IWString f = cl.string_value('F');
    if (f == "MF") {
      _mftype = formula_groupings::MFType::kMf;
      if (_verbose) {
        cerr << "Formula is default molecular formula\n";
      }
    } else if (f == "AROMF") {
      _mftype = formula_groupings::MFType::kAromMf;
      if (_verbose) {
        cerr << "Formula is aromatic distinguishing type\n";
      }
    } else {
      cerr << "Unrecognised -F qualifier - must be 'MF' or 'AROMF' only\n";
      return 0;
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', _name_stem);
    if (_verbose) {
      cerr << "Files created with stem '" << _name_stem << '\n';
    }
  } else {
    _name_stem = "formula_group_";
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Find " << _formula.size() << " different formula\n";
  uint64_t molecules_written = 0;
  for (const auto& [key, value] : _formula) {
    molecules_written += value.molecules_written();
  }
  cerr << "Wrote " << molecules_written << " molecules\n";

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
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  IWString formula;
  switch (_mftype) {
    case formula_groupings::MFType::kMf:
      formula = m.molecular_formula();
      break;
    case formula_groupings::MFType::kAromMf:
      m.formula_distinguishing_aromatic(formula);
      break;
    case formula_groupings::MFType::kUnspecified:
      cerr << "Options::Process:formula type not set?\n";
      return 0;
    default:
      cerr << "Options::Process:formula type not set?\n";
      return 0;
  };

  IWString to_store;
  to_store << m.smiles() << ' ' << m.name();

  auto iter = _formula.find(formula);
  if (iter != _formula.end()) {
    iter->second.Extra(std::move(to_store));
    _size += to_store.size();
    return MaybeWrite();
  }

  auto iter2 = _formula.emplace(formula, Formula());

  IWString fname = _name_stem;
  fname << formula << ".smi";
  std::get<0>(iter2)->second.set_file_name(fname);

  std::get<0>(iter2)->second.Extra(std::move(to_store));

  if (_verbose > 1) {
    cerr << "Options::Process:encountered " << formula << '\n';
  }

  return MaybeWrite();
}

int
Options::MaybeWrite() {
  if (_size < _write_threshold) {
    return 1;
  }

  // We need to flush something. Identify the formula with the largest amount stored.

  uint32_t max_size = 0;
  IWString formula_with_max_size;
  for (const auto& [key, value] : _formula) {
    if (value.size() > max_size) {
      max_size = value.size();
      formula_with_max_size = key;
    }
  }

  // todo:ianwatson  finish this sometime...

  return 1;
}

int
Options::Flush() {
  // First write any formula that have an open file.
  for (auto& [_, value] : _formula) {
    if (value.IsOpen()) {
      value.Close();
    }
  }

  for (auto& [_, value] : _formula) {
    if (value.IsOpen()) {
      continue;
    }

    value.Close();
  }

  return 1;
}

int
ApplicationName(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
ApplicationName(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! ApplicationName(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ApplicationName(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:S:F:");

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

  for (const char * fname : cl) {
    if (! ApplicationName(options, fname, input_type, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.Flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace formula_groupings

int
main(int argc, char ** argv) {

  int rc = formula_groupings::Main(argc, argv);

  return rc;
}
