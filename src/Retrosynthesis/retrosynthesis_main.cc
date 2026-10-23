#include <cstdint>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/molecule.h"

#include "Retrosynthesis/retrosynthesis.h"

namespace retrosynthesis {

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
  cerr << R"(Retrosynthetic tool using reverse reactions.
 -C <fname>             Retrosynthesis::RetrosynthesisData textproto with reaction and reagent information.
 -v                     verbose output.
)";
  // clang-format on
  ::exit(rc);
}

class Options {
  private:
    int _verbose;
    uint64_t _molecules_read;

    Preprocessor _preprocess;

    // The thing that does all the work.
    Retrosynthesis _retrosynthesis;

    // Count of molecules that never match any of the top level queries.
    uint64_t _no_matches;

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _molecules_read = 0;
  _no_matches = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! _preprocess.Initialise(cl)) {
    cerr << "Retrosynthesis::Initialise:cannot initialise preprocessing\n";
    return 0;
  }

  if (! cl.option_present('C')) {
    cerr << "Options::Initialise:must specify data proto via the -C option\n";
    return 0;
  }

  if (cl.option_present('C')) {
    if (! _retrosynthesis.Initialise(cl, _preprocess)) {
      cerr << "Cannot initialise retrosynthesis data\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {

  ++_molecules_read;

  _preprocess.Process(m);

  MoleculeAndFragmentation result;
  result.set_parent_molecule(m);
  cerr << "Calling _retrosynthesis\n";
  if (! _retrosynthesis.Process(m, result)) {
    return 1;
  }

  cerr << "REsult is " << result << '\n';
  result.DebugPrint("", cerr);
  if (result.empty()) {
    ++_no_matches;
    return 1;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _no_matches << " molecules did not match any query\n";

  return 1;
}

int
RunRetrosynthesis(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
RunRetrosynthesis(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! RunRetrosynthesis(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
RunRetrosynthesis(Options& options,
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

  return RunRetrosynthesis(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:C:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(1);
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
    if (! RunRetrosynthesis(options, fname, input_type, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace retrosynthesis


int
main(int argc, char ** argv) {

  int rc = retrosynthesis::Main(argc, argv);

  return rc;
}
