#include <stdlib.h>

#include <iostream>

#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/scaffolds.h"

#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools/scaffolds.pb.h"

namespace scaffolds_main {

using std::cerr;
using scaffolds::ScaffoldFinder;

class LocalOptions {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    int _non_ring_molecules_skipped;

    Chemical_Standardisation _chemical_standardisation;

    int _write_parent;

    // Output can be either as smiles or proto.
    int _write_smiles;

    // We can also accumulate all the scaffolds encountered.
    int _accumulate_all_scaffolds;
    IWString_and_File_Descriptor _stream_for_all_scaffolds;
    absl::flat_hash_map<std::string, dicer_data::DicerFragment> _all_scaffolds;

  // private functions
 
    int Accumulate(Molecule& parent, const scaffolds::ScaffoldData& result);
    int WriteSmiles(Molecule& parent,
                    const scaffolds::ScaffoldData& result,
                    IWString_and_File_Descriptor& output) const;
    int WriteAccumulated(const std::string& usmi,
                const dicer_data::DicerFragment& proto);

  public:
    LocalOptions();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    void AnotherNonRingMolecule() {
      ++_non_ring_molecules_skipped;
    }

    // If accumulation has not been abled, this does nothing.
    int WriteAccumulated();

    int Write(Molecule& parent,
              const scaffolds::ScaffoldData& result,
              IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

LocalOptions::LocalOptions() {
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _non_ring_molecules_skipped = 0;
  _write_parent = 0;
  _write_smiles = 1;
  _accumulate_all_scaffolds = 0;
}

int
LocalOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('p')) {
    _write_parent = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('y')) {
    _write_smiles = 0;
    if (_verbose) {
      cerr << "Will write as proto form\n";
    }
  }

  if (cl.option_present('W')) {
    IWString fname = cl.string_value('W');
    fname.EnsureEndsWith(".textproto");
    if (! _stream_for_all_scaffolds.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for all scaffold subsets '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Proto data for all scaffold variants written to '" << fname << "'\n";
    }
    _accumulate_all_scaffolds = 1;
  }

  return 1;
}

int
LocalOptions::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (m.empty()) {
    return 0;
  }

  return 1;
}

int
LocalOptions::Write(Molecule& parent,
                    const scaffolds::ScaffoldData& result,
                    IWString_and_File_Descriptor& output) {
  if (_accumulate_all_scaffolds) {
    Accumulate(parent, result);
  }

  static constexpr char kSep = ' ';

  if (_write_smiles) {
    return WriteSmiles(parent, result, output);
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  if (_write_parent) {
    output << parent.smiles() << kSep << parent.name() << kSep << "parent\n";
  }

  std::string buffer;

  printer.PrintToString(result, &buffer);
  output << buffer << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
LocalOptions::WriteSmiles(Molecule& parent,
                    const scaffolds::ScaffoldData& result,
                    IWString_and_File_Descriptor& output) const {
  static constexpr char kSep = ' ';
  if (_write_parent) {
    output << parent.smiles() << kSep << parent.name() << '\n';
  }

  int ndx = 0;
  for (const auto& subset : result.subset()) {
    output << subset.smi() << kSep << parent.name() << '.' << ndx << kSep << subset.ring_sys() << '\n';
    ++ndx;
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
LocalOptions::Accumulate(Molecule& parent,
                         const scaffolds::ScaffoldData& result) {
  for (const auto& subset : result.subset()) {
    auto iter = _all_scaffolds.find(subset.smi());
    if (iter == _all_scaffolds.end()) {
      const std::string& smi = subset.smi();

      dicer_data::DicerFragment proto;
      proto.set_smi(smi);
      proto.set_par(parent.name().data(), parent.name().size());
      proto.set_n(1);
      _all_scaffolds[smi] = std::move(proto);
    } else {
      const auto n = iter->second.n();
      iter->second.set_n(n + 1);
    }
  }

  return 1;
}

// Not const because we write to _stream_for_all_scaffolds
int
LocalOptions::WriteAccumulated() {
  for (const auto& [usmi, proto] : _all_scaffolds) {
    WriteAccumulated(usmi, proto);
    _stream_for_all_scaffolds.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
LocalOptions::WriteAccumulated(const std::string& usmi,
                const dicer_data::DicerFragment& proto) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;

  printer.PrintToString(proto, &buffer);
  _stream_for_all_scaffolds << buffer << '\n';

  return 1;
}

void
Usage(int rc) {
  cerr << R"(Generates all scaffold combinations.
Generates all combinations of the ring systems in a molecule. For example of a molecule consists of
R1-R2-R3 the output will be {R1 R2 R3 R1-R2 R2-R3}.
 -G <fname>             config containing scaffolds::ScaffoldsOptions text proto
 -W <fname>             write dicer_data::DicerFragment data on all scaffolds
 -y                     output as scaffolds::ScaffoldData textproto
 -v                     verbose output
  )";

  ::exit(rc);
}
  
int
Scaffolds(LocalOptions& local_options,
          ScaffoldFinder& make_scaffolds,
          Molecule& m,
          IWString_and_File_Descriptor& output) {
  if (m.nrings() == 0) {
    local_options.AnotherNonRingMolecule();
    return 1;
  }

  scaffolds::ScaffoldData result;
  result.set_smi(m.smiles().AsString());
  result.set_par(m.name().AsString());

  if (! make_scaffolds.MakeScaffolds(m, result)) {
    cerr << "MakeScaffolds failed " << m.name() << '\n';
    return 1;
  }

  return local_options.Write(m, result, output);
}

int
Scaffolds(LocalOptions& local_options,
          ScaffoldFinder& make_scaffolds,
          data_source_and_type<Molecule>& input,
          IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! local_options.Preprocess(*m)) {
      return 0;
    }

    if (! Scaffolds(local_options, make_scaffolds, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Scaffolds(LocalOptions& local_options,
          ScaffoldFinder& make_scaffolds,
          const char * fname,
          FileType input_type,
          IWString_and_File_Descriptor& output) {

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return Scaffolds(local_options, make_scaffolds, input, output);
}

int
Scaffolds(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:E:g:clpG:yW:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    return 1;
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  LocalOptions local_options;
  if (! local_options.Initialise(cl)) {
    cerr << "Cannot initialise local options\n";
    return 1;
  }

  scaffolds::ScaffoldFinder make_scaffolds;

  if (! make_scaffolds.Initialise(cl, 'G')) {
    cerr << "Cannot initialise calculation\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  }

  if (FILE_TYPE_INVALID != input_type) {  // great, explicitly specified
  } else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  // reading from a pipe, assume smiles input
    if (verbose)
      cerr << "Assuming smiles input from pipe read\n";
    input_type = FILE_TYPE_SMI;
  } else if (all_files_recognised_by_suffix(cl)) {
  } else {
    cerr << "Cannot discern file types from names\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname: cl) {
    const const_IWSubstring as_string(fname);

    if (verbose) {
      cerr << "Processing '" << fname << "'\n";
    }

    if (! Scaffolds(local_options, make_scaffolds, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  // This is a no-op unless it has been activated.
  local_options.WriteAccumulated();

  if (verbose) {
    make_scaffolds.Report(cerr);
  }

  return 0;
}

}  // namespace scaffolds_main

int
main(int argc, char **argv) {
  int rc = scaffolds_main::Scaffolds(argc, argv);

  return rc;
}
