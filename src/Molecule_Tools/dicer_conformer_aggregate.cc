#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/SVD"
#include "Molecule_Tools/dicer_conformers.pb.h"
#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_conformer_aggregate {

using dicer_conformers::AtomRole;
using dicer_conformers::ATTACHMENT;
using dicer_conformers::FRAGMENT_HEAVY;
using iw_tf_data_record::TFDataReader;
using iw_tf_data_record::TFDataWriter;
using std::cerr;

struct Options {
  int verbose = 0;
  int min_attachments = 2;
  int max_attachments = 2;
  int max_embeddings = 1000;
  double rmsd_tolerance = 0.25;
  double attachment_tolerance = 0.50;
  bool input_is_tfdata = false;
  bool output_is_tfdata = false;
  IWString output_name;

  uint64_t molecules_read = 0;
  uint64_t fragments_read = 0;
  uint64_t fragments_selected = 0;
  uint64_t duplicates = 0;
  uint64_t rejected_degenerate_alignment = 0;
  uint64_t rejected_too_many_embeddings = 0;
};

void
Usage(int rc) {
  cerr << R"(Aggregate and deduplicate 3D fragments produced by dicer -I geom.
  -i textproto|tfdata   input format (default textproto)
  -o textproto|tfdata   output format (default textproto)
  -S <file>             output file (stdout for textproto by default)
  -m <n>                minimum attachment atoms to retain (default 2)
  -M <n>                maximum attachment atoms to retain (default 2)
  -r <distance>         fragment heavy-atom RMSD tolerance (default 0.25)
  -a <distance>         maximum attachment displacement (default 0.50)
  -e <n>                maximum symmetry embeddings per fragment (default 1000)
  -v                    verbose output

All explicit Hydrogen atoms, including external Hydrogen attachment records,
are discarded. At least three non-collinear fragment heavy atoms are required.
)";
  exit(rc);
}

bond_type_t
LillyMolBondType(dicer_data::AttachmentGeometry::BondType btype) {
  switch (btype) {
    case dicer_data::AttachmentGeometry::BOND_SINGLE:
      return SINGLE_BOND;
    case dicer_data::AttachmentGeometry::BOND_DOUBLE:
      return DOUBLE_BOND;
    case dicer_data::AttachmentGeometry::BOND_TRIPLE:
      return TRIPLE_BOND;
    case dicer_data::AttachmentGeometry::BOND_AROMATIC:
      return AROMATIC_BOND;
    default:
      return INVALID_BOND_TYPE;
  }
}

struct Candidate {
  Molecule molecule;
  IWString fragment_usmi;
  IWString key;
  std::vector<AtomRole> role;
  std::vector<float> coordinates;
  uint64_t occurrences = 1;
};

// Kept as a separate function because the precise aggregation identity is
// expected to evolve as attachment typing is exercised on real datasets.
IWString
AggregationKey(Molecule& augmented, const std::vector<AtomRole>& role) {
  const IWString& usmi = augmented.unique_smiles();
  const auto& order = augmented.atom_order_in_smiles();
  IWString result(usmi);
  result << '|';
  for (int i = 0; i < order.number_elements(); ++i) {
    const atom_number_t atom = order[i];
    result << static_cast<int>(role[atom]);
    if (role[atom] == ATTACHMENT) {
      result << ':' << augmented.atomic_number(atom) << ':' << augmented.isotope(atom);
    }
    result << ',';
  }
  return result;
}

std::optional<Candidate>
BuildCandidate(const dicer_data::DicerFragment& fragment) {
  if (fragment.smi().empty() || fragment.attachment().empty()) {
    return std::nullopt;
  }

  Molecule mol;
  if (!mol.build_from_smiles(fragment.smi())) {
    cerr << "BuildCandidate:cannot parse fragment smiles '" << fragment.smi() << "'\n";
    return std::nullopt;
  }

  const int initial_atoms = mol.natoms();
  std::vector<int> old_to_new(initial_atoms, -1);
  int next_atom = 0;
  for (int i = 0; i < initial_atoms; ++i) {
    if (mol.atomic_number(i) != 1) {
      old_to_new[i] = next_atom++;
    }
  }
  mol.RemoveAllHydrogenAtoms();

  Candidate result;
  result.molecule = std::move(mol);
  result.role.resize(result.molecule.natoms(), FRAGMENT_HEAVY);
  result.occurrences = fragment.has_n() ? fragment.n() : 1;

  for (const dicer_data::AttachmentGeometry& attachment : fragment.attachment()) {
    if (attachment.atom() >= old_to_new.size() || old_to_new[attachment.atom()] < 0) {
      continue;
    }
    Molecule external;
    if (!external.build_from_smiles(attachment.ext()) || external.natoms() != 1) {
      cerr << "BuildCandidate:invalid external atom smiles '" << attachment.ext()
           << "'\n";
      return std::nullopt;
    }
    if (external.atomic_number(0) == 1) {
      continue;
    }
    const bond_type_t btype = LillyMolBondType(attachment.btype());
    if (btype == INVALID_BOND_TYPE) {
      cerr << "BuildCandidate:invalid attachment bond type\n";
      return std::nullopt;
    }

    const atom_number_t added = result.molecule.natoms();
    result.molecule.add(new Atom(external.atomi(0)));
    result.molecule.add_bond(old_to_new[attachment.atom()], added, btype);
    result.role.push_back(ATTACHMENT);
  }

  // Fragment smiles retain an implicit-H count appropriate to the severed
  // fragment. Recompute after restoring the external attachment bonds, and
  // also erase any explicit-H bookkeeping inherited from the input.
  for (int i = 0; i < result.molecule.natoms(); ++i) {
    result.molecule.unset_all_implicit_hydrogen_information(i);
  }

  Molecule fragment_only(result.molecule);
  while (fragment_only.natoms() > next_atom) {
    fragment_only.remove_atom(fragment_only.natoms() - 1);
  }
  result.fragment_usmi = fragment_only.unique_smiles();
  result.key = AggregationKey(result.molecule, result.role);

  result.coordinates.reserve(3 * result.molecule.natoms());
  for (int i = 0; i < result.molecule.natoms(); ++i) {
    result.coordinates.push_back(result.molecule.x(i));
    result.coordinates.push_back(result.molecule.y(i));
    result.coordinates.push_back(result.molecule.z(i));
  }
  return result;
}

class ConformerSet {
 private:
  Molecule _topology;
  IWString _fragment_usmi;
  IWString _key;
  IWString _first_parent;
  std::vector<AtomRole> _role;
  resizable_array_p<Set_of_Atoms> _embedding;
  std::vector<std::vector<float>> _conformer;
  uint64_t _occurrences = 0;
  uint64_t _geometries_seen = 0;

  bool RolePreserving(const Set_of_Atoms& embedding) const;
  bool BuildEmbeddings(int max_embeddings);
  bool Compare(const std::vector<float>& stored, const std::vector<float>& incoming,
               const Set_of_Atoms& embedding, double rmsd_tolerance,
               double attachment_tolerance, Eigen::Affine3d& transform,
               double& rmsd) const;

 public:
  bool Initialise(Candidate&& candidate, const std::string& parent, int max_embeddings);
  // Returns 1 if retained, 0 if a duplicate.
  int Add(Candidate&& candidate, double rmsd_tolerance, double attachment_tolerance);
  void ToProto(dicer_conformers::FragmentConformerSet& proto) const;
  int attachment_count() const;
};

bool
ConformerSet::RolePreserving(const Set_of_Atoms& embedding) const {
  if (embedding.size() != _role.size()) {
    return false;
  }
  for (int i = 0; i < embedding.number_elements(); ++i) {
    if (_role[i] != _role[embedding[i]]) {
      return false;
    }
  }
  return true;
}

// Determine a rigid transform such that target ~= transform * source. Scaling
// is deliberately excluded: differences in bond lengths are conformational
// differences, not something alignment should erase.
bool
RigidAlignment(const std::vector<Eigen::Vector3d>& source,
               const std::vector<Eigen::Vector3d>& target, Eigen::Affine3d& transform,
               double& rmsd) {
  if (source.size() != target.size() || source.size() < 3) {
    return false;
  }

  Eigen::Vector3d source_centroid = Eigen::Vector3d::Zero();
  Eigen::Vector3d target_centroid = Eigen::Vector3d::Zero();
  for (int i = 0; i < static_cast<int>(source.size()); ++i) {
    source_centroid += source[i];
    target_centroid += target[i];
  }
  source_centroid /= source.size();
  target_centroid /= target.size();

  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  for (int i = 0; i < static_cast<int>(source.size()); ++i) {
    covariance += (source[i] - source_centroid) * (target[i] - target_centroid).transpose();
  }
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(covariance,
                                        Eigen::ComputeFullU | Eigen::ComputeFullV);
  // Rank one means that rotation about the common line is undefined.
  if (svd.singularValues()[1] < 1.0e-8) {
    return false;
  }

  Eigen::Matrix3d correction = Eigen::Matrix3d::Identity();
  if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < 0.0) {
    correction(2, 2) = -1.0;
  }
  const Eigen::Matrix3d rotation = svd.matrixV() * correction * svd.matrixU().transpose();
  transform = Eigen::Affine3d::Identity();
  transform.linear() = rotation;
  transform.translation() = target_centroid - rotation * source_centroid;

  rmsd = 0.0;
  for (int i = 0; i < static_cast<int>(source.size()); ++i) {
    const double distance = (target[i] - transform * source[i]).norm();
    rmsd += distance * distance;
  }
  rmsd = std::sqrt(rmsd / source.size());
  return true;
}

bool
ConformerSet::BuildEmbeddings(int max_embeddings) {
  Molecule_to_Query_Specifications mqs;
  mqs.set_make_embedding(1);
  Substructure_Query query;
  if (!query.create_from_molecule(_topology, mqs)) {
    return false;
  }
  query.set_max_matches_to_find(max_embeddings + 1);
  Molecule_to_Match target(&_topology);
  Substructure_Results results;
  const int nhits = query.substructure_search(target, results);
  if (nhits > max_embeddings) {
    return false;
  }

  resizable_array_p<Set_of_Atoms> tmp;
  tmp.reserve(results.number_embeddings());
  for (uint32_t i = 0; i < results.number_embeddings(); ++i) {
    const Set_of_Atoms* embedding = results.embedding(i);
    if (RolePreserving(*embedding)) {
      tmp.add(new Set_of_Atoms(*embedding));
    }
  }
  _embedding.transfer_in(tmp);
  return !_embedding.empty();
}

int
ConformerSet::attachment_count() const {
  return std::count(_role.begin(), _role.end(), ATTACHMENT);
}

bool
ConformerSet::Initialise(Candidate&& candidate, const std::string& parent,
                         int max_embeddings) {
  _topology = std::move(candidate.molecule);
  _fragment_usmi = std::move(candidate.fragment_usmi);
  _key = std::move(candidate.key);
  _first_parent = parent;
  _role = std::move(candidate.role);
  _occurrences = candidate.occurrences;
  _geometries_seen = 1;
  _conformer.push_back(std::move(candidate.coordinates));
  return BuildEmbeddings(max_embeddings);
}

bool
ConformerSet::Compare(const std::vector<float>& stored,
                      const std::vector<float>& incoming, const Set_of_Atoms& embedding,
                      double rmsd_tolerance, double attachment_tolerance,
                      Eigen::Affine3d& transform, double& rmsd) const {
  transform = Eigen::Affine3d::Identity();
  rmsd = std::numeric_limits<double>::max();
  std::vector<Eigen::Vector3d> target;
  std::vector<Eigen::Vector3d> source;
  for (int i = 0; i < static_cast<int>(_role.size()); ++i) {
    if (_role[i] != FRAGMENT_HEAVY) {
      continue;
    }
    target.emplace_back(stored[3 * i], stored[3 * i + 1], stored[3 * i + 2]);
    const int j = embedding[i];
    source.emplace_back(incoming[3 * j], incoming[3 * j + 1], incoming[3 * j + 2]);
  }
  if (!RigidAlignment(source, target, transform, rmsd)) {
    return false;
  }
  if (rmsd > rmsd_tolerance) {
    return false;
  }

  double max_attachment_displacement = 0.0;
  for (int i = 0; i < static_cast<int>(_role.size()); ++i) {
    if (_role[i] != ATTACHMENT) {
      continue;
    }
    const int j = embedding[i];
    Eigen::Vector3d xyz(incoming[3 * j], incoming[3 * j + 1], incoming[3 * j + 2]);
    xyz = transform * xyz;
    const Eigen::Vector3d reference(stored[3 * i], stored[3 * i + 1], stored[3 * i + 2]);
    max_attachment_displacement =
        std::max(max_attachment_displacement, (reference - xyz).norm());
  }
  return max_attachment_displacement <= attachment_tolerance;
}

int
ConformerSet::Add(Candidate&& candidate, double rmsd_tolerance,
                  double attachment_tolerance) {
  _occurrences += candidate.occurrences;
  ++_geometries_seen;

  Eigen::Affine3d best_transform = Eigen::Affine3d::Identity();
  double best_rmsd = std::numeric_limits<double>::max();
  for (const std::vector<float>& stored : _conformer) {
    for (const Set_of_Atoms* embedding : _embedding) {
      Eigen::Affine3d transform;
      double rmsd;
      if (Compare(stored, candidate.coordinates, *embedding, rmsd_tolerance,
                  attachment_tolerance, transform, rmsd)) {
        return 0;
      }
      if (rmsd < best_rmsd) {
        best_rmsd = rmsd;
        best_transform = transform;
      }
    }
  }

  // Store all conformers in the coordinate frame of the first exemplar.
  for (int i = 0; i < static_cast<int>(_role.size()); ++i) {
    Eigen::Vector3d xyz(candidate.coordinates[3 * i], candidate.coordinates[3 * i + 1],
                        candidate.coordinates[3 * i + 2]);
    xyz = best_transform * xyz;
    candidate.coordinates[3 * i] = xyz[0];
    candidate.coordinates[3 * i + 1] = xyz[1];
    candidate.coordinates[3 * i + 2] = xyz[2];
  }
  _conformer.push_back(std::move(candidate.coordinates));
  return 1;
}

void
ConformerSet::ToProto(dicer_conformers::FragmentConformerSet& proto) const {
  proto.set_usmi(_fragment_usmi.data(), _fragment_usmi.length());
  Molecule topology(_topology);
  const IWString& smiles = topology.unique_smiles();
  proto.set_smiles(smiles.data(), smiles.length());
  proto.set_first_parent(_first_parent.data(), _first_parent.length());
  proto.set_occurrences(_occurrences);
  proto.set_geometries_seen(_geometries_seen);
  for (AtomRole role : _role) {
    proto.add_atom_role(role);
  }
  for (const atom_number_t atom : topology.atom_order_in_smiles()) {
    proto.add_smiles_atom_order(atom);
  }
  for (const std::vector<float>& coordinates : _conformer) {
    auto* conformer = proto.add_conformer();
    for (float value : coordinates) {
      conformer->add_xyz(value);
    }
  }
}

using Groups = std::unordered_map<std::string, std::unique_ptr<ConformerSet>>;

bool
Process(const dicer_data::DicedMolecule& proto, Options& options, Groups& groups) {
  ++options.molecules_read;
  for (const dicer_data::DicerFragment& fragment : proto.fragment()) {
    ++options.fragments_read;
    std::optional<Candidate> candidate = BuildCandidate(fragment);
    if (!candidate) {
      continue;
    }
    const int attachments =
        std::count(candidate->role.begin(), candidate->role.end(), ATTACHMENT);
    if (attachments < options.min_attachments || attachments > options.max_attachments) {
      continue;
    }
    const int heavy =
        std::count(candidate->role.begin(), candidate->role.end(), FRAGMENT_HEAVY);
    std::vector<Eigen::Vector3d> alignment_atoms;
    for (int i = 0; i < static_cast<int>(candidate->role.size()); ++i) {
      if (candidate->role[i] == FRAGMENT_HEAVY) {
        alignment_atoms.emplace_back(candidate->coordinates[3 * i],
                                     candidate->coordinates[3 * i + 1],
                                     candidate->coordinates[3 * i + 2]);
      }
    }
    Eigen::Affine3d unused_transform;
    double unused_rmsd;
    if (heavy < 3 || !RigidAlignment(alignment_atoms, alignment_atoms, unused_transform,
                                     unused_rmsd)) {
      ++options.rejected_degenerate_alignment;
      continue;
    }
    ++options.fragments_selected;

    std::string key(candidate->key.data(), candidate->key.length());
    const auto iter = groups.find(key);
    if (iter == groups.end()) {
      auto group = std::make_unique<ConformerSet>();
      if (!group->Initialise(std::move(*candidate), proto.name(),
                             options.max_embeddings)) {
        ++options.rejected_too_many_embeddings;
        continue;
      }
      groups.emplace(std::move(key), std::move(group));
    } else if (!iter->second->Add(std::move(*candidate), options.rmsd_tolerance,
                                  options.attachment_tolerance)) {
      ++options.duplicates;
    }
  }
  return true;
}

bool
ReadTextproto(const char* fname, Options& options, Groups& groups) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return false;
  }
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    dicer_data::DicedMolecule proto;
    if (!google::protobuf::TextFormat::ParseFromString(
            std::string(buffer.data(), buffer.length()), &proto) ||
        !Process(proto, options, groups)) {
      cerr << "Invalid DicedMolecule textproto in '" << fname << "'\n";
      return false;
    }
  }
  return true;
}

bool
ReadTfData(const char* fname, Options& options, Groups& groups) {
  TFDataReader input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return false;
  }
  while (true) {
    std::optional<dicer_data::DicedMolecule> proto =
        input.ReadProto<dicer_data::DicedMolecule>();
    if (!proto) {
      return true;
    }
    if (!Process(*proto, options, groups)) {
      return false;
    }
  }
}

bool
WriteText(const Groups& groups, IWString_and_File_Descriptor& output) {
  std::vector<std::string> keys;
  keys.reserve(groups.size());
  for (const auto& [key, unused] : groups) {
    keys.push_back(key);
  }
  std::sort(keys.begin(), keys.end());
  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);
  for (const std::string& key : keys) {
    dicer_conformers::FragmentConformerSet proto;
    groups.at(key)->ToProto(proto);
    std::string buffer;
    if (!printer.PrintToString(proto, &buffer)) {
      return false;
    }
    if (!buffer.empty() && buffer.back() == ' ') {
      buffer.pop_back();
    }
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }
  return output.good();
}

bool
WriteTfData(const Groups& groups, IWString fname) {
  TFDataWriter output;
  if (!output.Open(fname.null_terminated_chars())) {
    cerr << "Cannot open output '" << fname << "'\n";
    return false;
  }
  std::vector<std::string> keys;
  for (const auto& [key, unused] : groups) {
    keys.push_back(key);
  }
  std::sort(keys.begin(), keys.end());
  for (const std::string& key : keys) {
    dicer_conformers::FragmentConformerSet proto;
    groups.at(key)->ToProto(proto);
    if (!output.WriteSerializedProto(proto)) {
      return false;
    }
  }
  return true;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vS:i:o:m:M:r:a:e:");
  if (cl.unrecognised_options_encountered() || cl.empty()) {
    Usage(1);
  }
  Options options;
  options.verbose = cl.option_count('v');
  if (cl.option_present('m') &&
      (!cl.value('m', options.min_attachments) || options.min_attachments < 0)) {
    cerr << "Invalid min_attachments (-m)\n";
    Usage(1);
  }
  if (cl.option_present('M') &&
      (!cl.value('M', options.max_attachments) || options.max_attachments < 0)) {
    cerr << "Invalid max_attachments (-M)\n";
    Usage(1);
  }
  if (options.max_attachments < options.min_attachments) {
    cerr << "Invalid attachment range\n";
    return 1;
  }
  if (cl.option_present('r') &&
      (!cl.value('r', options.rmsd_tolerance) || options.rmsd_tolerance < 0.0)) {
    cerr << "Invalid rmsd_tolerance (-r)\n";
    Usage(1);
  }
  if (cl.option_present('a') && (!cl.value('a', options.attachment_tolerance) ||
                                 options.attachment_tolerance < 0.0)) {
    cerr << "Invalid attachment_tolerance (-a)\n";
    Usage(1);
  }
  if (cl.option_present('e') &&
      (!cl.value('e', options.max_embeddings) || options.max_embeddings < 1)) {
    cerr << "Invalid max_embeddings (-e)\n";
    Usage(1);
  }
  if (cl.option_present('S')) {
    cl.value('S', options.output_name);
  }
  if (cl.option_present('i')) {
    const_IWSubstring input = cl.string_value('i');
    if (input == "tfdata") {
      options.input_is_tfdata = true;
    } else if (input != "textproto") {
      Usage(1);
    }
  }
  if (cl.option_present('o')) {
    const_IWSubstring output = cl.string_value('o');
    if (output == "tfdata") {
      options.output_is_tfdata = true;
    } else if (output != "textproto") {
      Usage(1);
    }
  }
  if (options.output_is_tfdata && options.output_name.empty()) {
    cerr << "TFDataRecord output requires -S <file>\n";
    return 1;
  }

  Groups groups;
  for (const char* fname : cl) {
    const bool ok = options.input_is_tfdata ? ReadTfData(fname, options, groups)
                                            : ReadTextproto(fname, options, groups);
    if (!ok) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  bool written;
  if (options.output_is_tfdata) {
    written = WriteTfData(groups, options.output_name);
  } else {
    std::unique_ptr<IWString_and_File_Descriptor> output;
    if (options.output_name.empty()) {
      output = std::make_unique<IWString_and_File_Descriptor>(1);
    } else {
      output = std::make_unique<IWString_and_File_Descriptor>();
      if (!output->open(options.output_name.null_terminated_chars())) {
        cerr << "Cannot open output '" << options.output_name << "'\n";
        return 1;
      }
    }
    written = WriteText(groups, *output);
    output->flush();
  }
  if (!written) {
    return 1;
  }

  if (options.verbose) {
    cerr << "Read " << options.molecules_read << " molecules and "
         << options.fragments_read << " fragments; selected "
         << options.fragments_selected << ", retained " << groups.size()
         << " fragment types, discarded " << options.duplicates
         << " duplicate geometries\n";
    cerr << options.rejected_degenerate_alignment
         << " fragments rejected because heavy atoms do not define a 3D alignment\n";
    cerr << options.rejected_too_many_embeddings
         << " fragment types rejected for excessive or invalid symmetry embeddings\n";
  }
  return 0;
}

}  // namespace dicer_conformer_aggregate

int
main(int argc, char** argv) {
  return dicer_conformer_aggregate::Main(argc, argv);
}
