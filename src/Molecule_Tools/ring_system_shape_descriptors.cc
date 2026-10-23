#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"

#include "Molecule_Tools/ring_system_shape.h"

namespace ring_system_shape_descriptors {

using std::cerr;

using molecule_processing::MoleculePreprocessing;
using ring_system_shape::AnalyseRingSystemShape;
using ring_system_shape::NonRingBranchPointCount;
using ring_system_shape::RingAtomBranchPointCount;
using ring_system_shape::RingSystemShape;
using ring_system_shape::RingSystemShapeClass;

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
  cerr << R"(Computes simple ring-system rod-like shape descriptors.

For each molecule, ring systems are labelled and each ring system is classified
by the number and placement of meaningful exit points.

Output columns are:
  Name nringsys nringsys_terminal nringsys_applicable nringsys_rodlike nringsys_not_rodlike nringsys_multisub nringsys_invalid non_ring_branch_count ring_atom_branch_count rodlike_molecule rod_deficit_min rod_deficit_max rod_deficit_ave

Options:
 -x           ignore terminal single-atom ring substituents when larger substituents exist
 -I <iso>     only consider exit vectors from ring atoms with isotope <iso>
 -R <fname>   write rod-like molecules as smiles and identifier
 -N <fname>   write non-rod-like molecules as smiles and identifier
 -s           write a smiles file rather than a descriptor file
 -o <sep>     output separator, recognised names include tab and space
 -i <type>    input type
 -g ...       chemical standardisation
 -l           reduce to largest fragment
 -c           remove chirality
 -A ...       aromaticity options
 -E ...       element options
 -v           verbose output
)";
  // clang-format on

  ::exit(rc);
}

struct RingSystemShapeDescriptors {
  int number_ring_systems = 0;
  int terminal = 0;
  int applicable = 0;
  int rod_like = 0;
  int not_rod_like = 0;
  int multi_substituted = 0;
  int invalid = 0;
  int non_ring_branch_count = 0;
  int ring_atom_branch_count = 0;

  int rod_deficit_count = 0;
  int rod_deficit_min = 0;
  int rod_deficit_max = 0;
  int rod_deficit_sum = 0;

  void
  AddRodDeficit(int rod_deficit) {
    if (rod_deficit < 0) {
      return;
    }

    if (rod_deficit_count == 0) {
      rod_deficit_min = rod_deficit;
      rod_deficit_max = rod_deficit;
    } else {
      if (rod_deficit < rod_deficit_min) {
        rod_deficit_min = rod_deficit;
      }
      if (rod_deficit > rod_deficit_max) {
        rod_deficit_max = rod_deficit;
      }
    }

    ++rod_deficit_count;
    rod_deficit_sum += rod_deficit;
  }

  float
  AverageRodDeficit() const {
    if (rod_deficit_count == 0) {
      return 0.0f;
    }

    return iwmisc::Fraction<float>(rod_deficit_sum, rod_deficit_count);
  }

  int
  RodLikeMolecule() const {
    if (number_ring_systems == 0) {
      return 0;
    }

    return not_rod_like == 0 && multi_substituted == 0 && invalid == 0 &&
           non_ring_branch_count == 0 && ring_atom_branch_count == 0;
  }
};

class Options {
 private:
  int _verbose = 0;

  MoleculePreprocessing _preprocessing;

  char _output_separator = ' ';

  int _remove_terminal_single_atom_substituents = 0;

  isotope_t _only_process_isotope = 0;

  IWString_and_File_Descriptor _stream_for_rodlike;
  IWString_and_File_Descriptor _stream_for_non_rodlike;

  // By default we write a descriptor file, but we can also write a smiles file.
  bool _write_smiles = false;

  uint64_t _molecules_read = 0;
  uint64_t _atoms_removed = 0;
  uint64_t _rodlike_molecules_written = 0;
  uint64_t _non_rodlike_molecules_written = 0;

  uint64_t _total_ring_systems = 0;
  uint64_t _rodlike_molecules = 0;
  uint64_t _molecules_with_non_ring_branching = 0;
  uint64_t _molecules_with_ring_atom_branching = 0;

 public:
  Options();

  int Initialise(Command_Line& cl);

  int
  verbose() const {
    return _verbose;
  }

  int Preprocess(Molecule& m);

  int RemoveTerminalSingleAtomSubstituents(Molecule& m);

  int WriteIfRequested(const Molecule& m, const IWString& smiles,
                       const RingSystemShapeDescriptors& descriptors);

  bool write_smiles() const {
    return _write_smiles;
  }
  int WriteHeader(IWString_and_File_Descriptor& output) const;

  int Process(Molecule& m, IWString_and_File_Descriptor& output);

  void FlushOutputStreams();

  int Report(std::ostream& output) const;
};

Options::Options() {
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (!_preprocessing.Initialise(cl)) {
    cerr << "Options::Initialise:cannot initialise preprocessing\n";
    return 0;
  }

  if (cl.option_present('o')) {
    IWString sep = cl.string_value('o');
    if (!char_name_to_char(sep)) {
      cerr << "Invalid output separator '" << sep << "'\n";
      return 0;
    }
    _output_separator = sep[0];
  }

  if (cl.option_present('x')) {
    _remove_terminal_single_atom_substituents = 1;
    if (_verbose) {
      cerr << "Terminal single atom substituents removed\n";
    }
  }

  if (cl.option_present('I')) {
    int isotope;
    if (!cl.value('I', isotope) || isotope <= 0) {
      cerr << "The -I option requires a positive isotope value\n";
      return 0;
    }
    _only_process_isotope = isotope;
    if (_verbose) {
      cerr << "Only considering ring atoms with isotope " << _only_process_isotope
           << " as exit-vector atoms\n";
    }
  }

  if (cl.option_present('R')) {
    IWString fname = cl.string_value('R');
    fname.EnsureEndsWith(".smi");
    if (!_stream_for_rodlike.open(fname.null_terminated_chars())) {
      cerr << "Options::Initialise:cannot open rod-like output '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Rod-like molecules written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('N')) {
    IWString fname = cl.string_value('N');
    fname.EnsureEndsWith(".smi");
    if (!_stream_for_non_rodlike.open(fname.null_terminated_chars())) {
      cerr << "Options::Initialise:cannot open non-rod-like output '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "non-rod-like molecules written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('s')) {
    _write_smiles = true;
    if (_write_smiles) {
      cerr << "Will write smiles\n";
    }
  }

  return 1;
}

int
Options::RemoveTerminalSingleAtomSubstituents(Molecule& m) {
  const int matoms = m.natoms();
  if (matoms == 0) {
    return 1;
  }

  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(matoms);
  const int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());
  if (number_ring_systems == 0) {
    return 1;
  }

  std::vector<Set_of_Atoms> single_atom_substituents(number_ring_systems + 1);
  std::vector<int> has_larger_attachment(number_ring_systems + 1, 0);

  for (int i = 0; i < matoms; ++i) {
    const int ring_system_id = ring_system_membership[i];
    if (ring_system_id <= 0) {
      continue;
    }

    const Atom& atom = m[i];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(i);
      if (ring_system_membership[other] == ring_system_id) {
        continue;
      }

      if (bond->is_single_bond() && m.ncon(other) == 1) {
        single_atom_substituents[ring_system_id].add_if_not_already_present(other);
      } else if (bond->is_single_bond() || m.ncon(other) > 1) {
        has_larger_attachment[ring_system_id] = 1;
      }
    }
  }

  Set_of_Atoms atoms_to_remove;
  for (int ring_system_id = 1; ring_system_id <= number_ring_systems; ++ring_system_id) {
    if (!has_larger_attachment[ring_system_id]) {
      continue;
    }

    for (atom_number_t atom : single_atom_substituents[ring_system_id]) {
      atoms_to_remove.add_if_not_already_present(atom);
    }
  }

  if (atoms_to_remove.empty()) {
    return 1;
  }

  const int atoms_removed = m.remove_atoms(atoms_to_remove);
  _atoms_removed += atoms_removed;

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_preprocessing.active()) {
    _preprocessing.Process(m);
  }

  return 1;
}

void
Options::FlushOutputStreams() {
  if (_stream_for_rodlike.is_open()) {
    _stream_for_rodlike.flush();
  }
  if (_stream_for_non_rodlike.is_open()) {
    _stream_for_non_rodlike.flush();
  }
}

int
Options::Report(std::ostream& output) const {
  output << "Processed " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  output << _total_ring_systems << " ring systems, "
         << iwmisc::Fraction<float>(_total_ring_systems, _molecules_read)
         << " per molecule\n";
  output << _rodlike_molecules << " molecules classified as rod-like, fraction "
         << iwmisc::Fraction<float>(_rodlike_molecules, _molecules_read) << '\n';
  output << _molecules_with_non_ring_branching
         << " molecules with non-ring branch points, fraction "
         << iwmisc::Fraction<float>(_molecules_with_non_ring_branching, _molecules_read)
         << '\n';
  output << _molecules_with_ring_atom_branching
         << " molecules with ring-atom branch points, fraction "
         << iwmisc::Fraction<float>(_molecules_with_ring_atom_branching, _molecules_read)
         << '\n';

  if (_rodlike_molecules_written) {
    output << "Wrote " << _rodlike_molecules_written << " rod-like molecules\n";
  }
  if (_non_rodlike_molecules_written) {
    output << "Wrote " << _non_rodlike_molecules_written << " non-rod-like molecules\n";
  }

  if (_remove_terminal_single_atom_substituents) {
    output << "Removed " << _atoms_removed << " terminal single atom substituents\n";
  }
  return 1;
}

RingSystemShapeDescriptors
ComputeRingSystemShapeDescriptors(Molecule& m, int non_ring_branch_count,
                                  int ring_atom_branch_count,
                                  isotope_t only_process_isotope) {
  RingSystemShapeDescriptors result;
  result.non_ring_branch_count = non_ring_branch_count;
  result.ring_atom_branch_count = ring_atom_branch_count;

  const int matoms = m.natoms();
  if (matoms == 0) {
    return result;
  }

  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(matoms);
  result.number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());

  if (result.number_ring_systems == 0) {
    return result;
  }

  for (int ring_system_id = 1; ring_system_id <= result.number_ring_systems;
       ++ring_system_id) {
    RingSystemShape ring_system_shape;
    if (!AnalyseRingSystemShape(m, ring_system_membership.get(), ring_system_id,
                                only_process_isotope, ring_system_shape)) {
      ++result.invalid;
      continue;
    }

    switch (ring_system_shape.shape_class) {
      case RingSystemShapeClass::kNotApplicable:
        ++result.terminal;
        break;
      case RingSystemShapeClass::kRodLike:
        ++result.applicable;
        ++result.rod_like;
        result.AddRodDeficit(ring_system_shape.rod_deficit);
        break;
      case RingSystemShapeClass::kNotRodLike:
        ++result.applicable;
        ++result.not_rod_like;
        result.AddRodDeficit(ring_system_shape.rod_deficit);
        break;
      case RingSystemShapeClass::kMultiSubstituted:
        ++result.multi_substituted;
        break;
      case RingSystemShapeClass::kInvalid:
        ++result.invalid;
        break;
    }
  }

  return result;
}

int
Options::WriteIfRequested(const Molecule& m, const IWString& smiles,
                          const RingSystemShapeDescriptors& descriptors) {
  const int is_rodlike =
      descriptors.RodLikeMolecule() && descriptors.non_ring_branch_count == 0;

  if (is_rodlike && _stream_for_rodlike.is_open()) {
    _stream_for_rodlike << smiles << ' ' << m.name() << '\n';
    _stream_for_rodlike.write_if_buffer_holds_more_than(4096);
    ++_rodlike_molecules_written;
  }

  if (!is_rodlike && _stream_for_non_rodlike.is_open()) {
    _stream_for_non_rodlike << smiles << ' ' << m.name() << '\n';
    _stream_for_non_rodlike.write_if_buffer_holds_more_than(4096);
    ++_non_rodlike_molecules_written;
  }

  return 1;
}

int
Options::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "Name" << _output_separator << "nringsys" << _output_separator
         << "nringsys_terminal" << _output_separator << "nringsys_applicable"
         << _output_separator << "nringsys_rodlike" << _output_separator
         << "nringsys_not_rodlike" << _output_separator << "nringsys_multisub"
         << _output_separator << "nringsys_invalid" << _output_separator
         << "non_ring_branch_count" << _output_separator << "ring_atom_branch_count"
         << _output_separator << "rodlike_molecule" << _output_separator
         << "rod_deficit_min" << _output_separator << "rod_deficit_max"
         << _output_separator << "rod_deficit_ave\n";
  return 1;
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  IWString smiles;
  if (_stream_for_rodlike.is_open() || _stream_for_non_rodlike.is_open()) {
    smiles = m.smiles();
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(matoms);
  m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());

  const int non_ring_branch_count =
      NonRingBranchPointCount(m, ring_system_membership.get());
  const int ring_atom_branch_count =
      RingAtomBranchPointCount(m, ring_system_membership.get(),
                               _only_process_isotope);

  if (_remove_terminal_single_atom_substituents) {
    RemoveTerminalSingleAtomSubstituents(m);
  }

  const RingSystemShapeDescriptors descriptors =
      ComputeRingSystemShapeDescriptors(m, non_ring_branch_count, ring_atom_branch_count,
                                        _only_process_isotope);

  _total_ring_systems += descriptors.number_ring_systems;
  if (descriptors.RodLikeMolecule()) {
    ++_rodlike_molecules;
  }
  if (descriptors.non_ring_branch_count > 0) {
    ++_molecules_with_non_ring_branching;
  }
  if (descriptors.ring_atom_branch_count > 0) {
    ++_molecules_with_ring_atom_branching;
  }

  WriteIfRequested(m, smiles, descriptors);

  if (_write_smiles) {
    output << m.smiles() << _output_separator;
  }

  output << m.name() << _output_separator << descriptors.number_ring_systems
         << _output_separator << descriptors.terminal << _output_separator
         << descriptors.applicable << _output_separator << descriptors.rod_like
         << _output_separator << descriptors.not_rod_like << _output_separator
         << descriptors.multi_substituted << _output_separator << descriptors.invalid
         << _output_separator << descriptors.non_ring_branch_count << _output_separator
         << descriptors.ring_atom_branch_count << _output_separator
         << descriptors.RodLikeMolecule() << _output_separator
         << descriptors.rod_deficit_min << _output_separator
         << descriptors.rod_deficit_max << _output_separator
         << descriptors.AverageRodDeficit() << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
RingSystemShapeDescriptors(Options& options, Molecule& m,
                           IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
RingSystemShapeDescriptors(Options& options, data_source_and_type<Molecule>& input,
                           IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!options.Preprocess(*m)) {
      continue;
    }

    if (!RingSystemShapeDescriptors(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
RingSystemShapeDescriptors(Options& options, const char* fname, FileType input_type,
                           IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "RingSystemShapeDescriptors:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return RingSystemShapeDescriptors(options, input, output);
}

int
RingSystemShapeDescriptors(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:o:xR:N:sI:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(1);
  }
  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (cl.number_elements() == 1 && strcmp(cl[0], "-") == 0) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  if (! options.write_smiles()) {
    options.WriteHeader(output);
  }

  for (const char* fname : cl) {
    if (!RingSystemShapeDescriptors(options, fname, input_type, output)) {
      cerr << "RingSystemShapeDescriptors::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();
  options.FlushOutputStreams();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace ring_system_shape_descriptors

int
main(int argc, char** argv) {
  return ring_system_shape_descriptors::RingSystemShapeDescriptors(argc, argv);
}
