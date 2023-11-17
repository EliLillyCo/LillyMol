// Main for minor_changes.

#include <stdlib.h>
#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/istream_and_type.h"

#include "Molecule_Tools/minor_changes.h"

namespace minor_changes_main {

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
  cerr << "Makes a potentially large number of small changes to a molecule.\n";
  cerr << " -C <fname>    textproto file with MinorChangesData options\n";
  cerr << " -F <fname>    one or more DicerFragment textproto files (from get_substituents)\n";
  cerr << " -B <fname>    one or more DicerFragment textproto files with bivalent (two isotopes) fragments\n";
  cerr << " -P <atype>    atom typing specification -needed if atom types in the get_substituents data\n";
  cerr << " -s <smarts>   specify atoms that can change\n";
  cerr << " -q <query>    specify atoms that can change\n";
  cerr << " -x            remove isotopes from product molecules\n";
  cerr << " -u <support>  support level for inclusion from fragment libraries\n";
  cerr << " -M <max>      maximum number of products per starting molecule\n";
  cerr << " -c            remove chirality\n";
  cerr << " -l            strip to largest fragment\n";
  cerr << " -g ...        chemical standardisation\n";
  cerr << " -v            verbose output\n";
// clang-format on

  ::exit(rc);
}

class LocalOptions {
  private:
    int _verbose;

    int _write_parent;

    int _remove_chirality;

    int _reduce_to_largest_fragment;

    Chemical_Standardisation _chemical_standardisation;

    IWDigits _digits;

    // Molecules are produces with isotopes. We can remove them.
    int _remove_isotopes;

  public:
    LocalOptions();

    int Initialise(Command_Line& cl);

    int write_parent() const {
      return _write_parent;
    }

    int remove_isotopes() const {
      return _remove_isotopes;
    }

    int Preprocess(Molecule& m);

    int AppendNumber(int number, IWString& destination) {
      return _digits.append_number(destination, number);
    }
};

LocalOptions::LocalOptions() : _digits(1000){
  _verbose = 0;
  _write_parent = 0;
  _remove_chirality = 0;
  _reduce_to_largest_fragment = 0;

  _digits.set_leading_string('.');
  _digits.append_to_each_stored_string('\n');

  _remove_isotopes = 0;
}

int
LocalOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('p')) {
    _write_parent = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will strip molecules to the largest fragment\n";
    }
  }

  if (cl.option_present('x')) {
    _remove_isotopes = 1;
    if (_verbose) {
      cerr << "Will remove isotopes from produce molecules\n";
    }
  }

  return 1;
}

int
LocalOptions::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
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
MinorChanges(LocalOptions& local_options,
                minor_changes::Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  if (! local_options.Preprocess(m)) {
    return 1;  // ignore.
  }

  resizable_array_p<Molecule> results;
  if (options.Process(m, results) < 0) {
    return 0;
  }

  if (results.empty()) {
    return 1;
  }

  static constexpr char kSep = ' ';

  if (local_options.write_parent()) {
    output << m.smiles() << kSep << m.name() << kSep << results.size() << '\n';
  }

  const int nvariants = results.number_elements();

  // Avoid reconstructing parts of the output.
  IWString fixed;
  fixed << kSep << m.name();

  for (int i = 0; i < nvariants; ++i) {
    if (local_options.remove_isotopes()) {
      results[i]->transform_to_non_isotopic_form();
    } 

    output << results[i]->smiles() << fixed;

    local_options.AppendNumber(i, output);
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
MinorChanges(LocalOptions& local_options,
                minor_changes::Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! MinorChanges(local_options, options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
MinorChanges(LocalOptions& local_options,
                minor_changes::Options& options,
                const char * fname,
                FileType input_type,
                IWString_and_File_Descriptor& output) {

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "MinorChanges:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return MinorChanges(local_options, options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:P:F:pC:B:xu:M:s:q:");

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

  LocalOptions local_options;
  if (! local_options.Initialise(cl)) {
    cerr << "Cannot initialise local options\n";
    return 1;
  }

  minor_changes::Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise minor_changes options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
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

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! MinorChanges(local_options, options, fname, input_type, output)) {
      cerr << "MinorChanges::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace minor_changes_main

int
main(int argc, char ** argv) {

  int rc = minor_changes_main::Main(argc, argv);

  return rc;
}
