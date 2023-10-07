// Generate 3D atom distance fingerprint.

#include <iostream>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

#include "absl/algorithm/container.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/atomic_distance_fingerprint.pb.h"

namespace three_dimensional_fp {

using std::cerr;
using iw_tf_data_record::TFDataWriter;

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
  cerr << "Generates raw interatomic distances\n";
  cerr << "For all pairs of atoms, form a bit that depends on the atom type\n";
  cerr << "For all atom type pairs, accumulate the distances and write to\n";
  cerr << "three_dimensional_fp::Fingerprint serialised proto form\n";
  cerr << " -T <atype>     atom typing specification (default UST:AY)\n";
  cerr << " -c             remove chirality\n";
  cerr << " -X             remove explicit Hydrogens\n";
  cerr << " -S <fname>     write results to <fname>\n";
  cerr << " -v             verbose output\n";
// clang-format on

  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    int _remove_explicit_hydrogens = 0;

    Chemical_Standardisation _chemical_standardisation;

    Atom_Typing_Specification _atom_typing;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    // We can exclude distances based on a range.
    float _min_distance;
    float _max_distance;

    // Any distance beyond this gets truncated.
    float _distance_range;

    int _molecules_read = 0;

    extending_resizable_array<uint32_t> _number_pairs;
    uint32_t _molecules_written;

    Accumulator<double> _acc_dist;

  // private functions

    std::optional<uint32_t> DistanceToBucket(float d) const;
    bool OkDistance(float d) const;
    int WriteProto(Molecule& m, const std::unordered_map<uint32_t, std::vector<float>>& pair_to_distances,
                TFDataWriter& output);

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

    int Process(Molecule& mol, TFDataWriter& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;

  _min_distance = 2.0;
  _max_distance = std::numeric_limits<float>::max();

  _distance_range = 30.0;

  _remove_explicit_hydrogens = 0;

  _molecules_written = 0;
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

  if (cl.option_present('P')) {
    const IWString p = cl.string_value('P');
    if (! _atom_typing.build(p)) {
      cerr << "Options::Initialise:cannot process atom typing '" << p << "'\n";
      return 0;
    }
  } else {
    _atom_typing.build("UST:AY");
    if (_verbose) {
      cerr << "Default atom typeing 'UST:AY'\n";
    }
  }

  if (cl.option_present('X')) {
    _remove_explicit_hydrogens = 1;
    if (_verbose) {
      cerr << "Explicit Hydrogens will be removed\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules, wrote " << _molecules_written << " fingerprints\n";
  for (int i = 0; i < _number_pairs.number_elements(); ++i) {
    if (_number_pairs[i]) {
      output << _number_pairs[i] << " molecules generated " << i << " distance pairs\n";
    }
  }
  output << "Distances btw " << _acc_dist.minval() << " and " << _acc_dist.maxval() << " mean " << _acc_dist.average() << '\n';

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
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  if (_remove_explicit_hydrogens) {
    m.remove_explicit_hydrogens();
  }

  return 1;
}

uint32_t
BitNumber(uint32_t t1, uint32_t t2) {
  static constexpr uint32_t kMagic = 5009;
  if (t1 < t2) {
    return t1 * kMagic + t2;
  }
  return t2 * kMagic + t1;
}

std::optional<uint32_t>
Options::DistanceToBucket(float d) const {
  if (d < _min_distance) {
    return std::nullopt;
  }

  if (d > _max_distance) {
    return std::nullopt;
  }

  return static_cast<uint32_t>(d / _distance_range * 255 + 0.4999);
}

bool
Options::OkDistance(float d) const {
  if (d < _min_distance) {
    return false;
  }

  if (d > _max_distance) {
    return false;
  }

  return true;
}

void
SortDistances(std::unordered_map<uint32_t, std::vector<float>>& pair_to_distances) {
  for (auto& [_, dists] : pair_to_distances) {
    if (dists.size() == 1) {
      continue;
    }
    if (dists.size() == 2) {
      if (dists[0] > dists[1]) {
        std::swap(dists[0], dists[1]);
      }
      continue;
    }

    absl::c_sort(dists, [](float d1, float d2) {
      return d1 < d2;
    });
  }

#define CHECK_SORT
#ifdef CHECK_SORT
  for (auto& [bit, dists] : pair_to_distances) {
    if (dists.empty()) {
      cerr << "Empty set of distances for bit " << bit << '\n';
      continue;
    }
    if (dists.size() == 1) {
      continue;
    }
    for (uint32_t i = 1; i < dists.size(); ++i) {
      if (dists[i] < dists[i-1]) {
        cerr << "Out of order sort, bit " << bit << " i " << i << " prev " << dists[i-1] << " cmp " << dists[i] << '\n';
      }
    }
  }
#endif
} 

int
Options::Process(Molecule& m,
                 TFDataWriter& output) {
  ++_molecules_read;

  const int matoms = m.natoms();

  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);

  if (! _atom_typing.assign_atom_types(m, atype.get())) {
    cerr << "Options::Process:cannot assign atom types " << m.name() << "'\n";
    return 0;
  }

  std::unordered_map<uint32_t, std::vector<float>> pair_to_distances;

  for (int i = 0; i < matoms; ++i) {
    uint32_t atype_i = atype[i];
    for (int j = i + 1; j < matoms; ++j) {
      const float d = m.distance_between_atoms(i, j);
      if (! OkDistance(d)) {
        continue;
      }

      _acc_dist.extra(d);

      uint32_t bit = BitNumber(atype_i, atype[j]);
      auto iter = pair_to_distances.find(bit);
      if (iter == pair_to_distances.end()) {
        std::vector<float> dists;
        dists.push_back(d);
        pair_to_distances.emplace(std::make_pair(bit, std::move(dists)));
      } else {
        iter->second.push_back(d);
      }
    }
  }

  if (_verbose > 1) {
    cerr << m.name() << " generated " << pair_to_distances.size() << " pairs\n";
  }

  SortDistances(pair_to_distances);

  return WriteProto(m, pair_to_distances, output);
}

int
Options::WriteProto(Molecule& m,
                const std::unordered_map<uint32_t, std::vector<float>>& pair_to_distances,
                TFDataWriter& output) {
  ++_number_pairs[pair_to_distances.size()];
  ++_molecules_written;

  atomic_distance_fingerprint::Fingerprint proto;
  const IWString& mname = m.name();
  proto.set_name(mname.data(), mname.length());
  const IWString& smiles = m.smiles();
  proto.set_smiles(smiles.data(), smiles.length());

  for (const auto& [bit, dists] : pair_to_distances) {
    atomic_distance_fingerprint::Distances distances;
    for (uint32_t d : dists) {
      distances.add_distance(d);
    }
    (*proto.mutable_distances())[bit] = distances;
  }

  return output.WriteSerializedProto<atomic_distance_fingerprint::Fingerprint>(proto);
}

int
AtomicDistanceFingerprint(Options& options,
                Molecule& m,
                TFDataWriter& output) {
  return options.Process(m, output);
}

int
AtomicDistanceFingerprint(Options& options,
                data_source_and_type<Molecule>& input,
                TFDataWriter& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! AtomicDistanceFingerprint(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
AtomicDistanceFingerprint(Options& options,
             const char * fname,
             FileType input_type,
             TFDataWriter& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "AtomicDistanceFingerprint:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return AtomicDistanceFingerprint(options, input, output);
}

int
AtomicDistanceFingerprint(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:P:S:X");

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
    input_type = FILE_TYPE_SDF;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name via the -S option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  TFDataWriter output;
  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! output.Open(fname)) {
      cerr << "AtomicDistanceFingerprint::cannot open '" << fname << "'\n";
      return 1;
    }
    if (verbose) {
      cerr << "output to '" << fname << "'\n";
    }
  }

  for (const char * fname : cl) {
    if (! AtomicDistanceFingerprint(options, fname, input_type, output)) {
      cerr << "AtomicDistanceFingerprint::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace three_dimensional_fp

int
main(int argc, char ** argv) {

  int rc = three_dimensional_fp::AtomicDistanceFingerprint(argc, argv);

  return rc;
}
