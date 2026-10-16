// Tool-style bindings for the nanobind LillyMol pilot.

#include "nanobind/lillymol_nb_internal.h"

#include <cmath>
#include <memory>
#include <vector>

#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>

#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Tools/dicer_api.h"
#include "Molecule_Tools/medchemwizard_lib.h"
#include "Molecule_Tools/ring_replacement_lib.h"
#include "Molecule_Tools/unique_molecules_api.h"
#include "Utilities/GFP_Tools/gfp_context.h"
#include "Utilities/GFP_Tools/truncated_distance_matrix.h"

namespace lillymol_nb {
namespace {

void
CheckDistanceMatrixIndex(const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
                         uint32_t ndx) {
  if (ndx >= dm.size()) {
    throw std::out_of_range("TruncatedDistanceMatrix index out of range");
  }
}

std::string
TruncatedDistanceMatrixName(const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
                            uint32_t ndx) {
  CheckDistanceMatrixIndex(dm, ndx);
  return dm.Name(ndx);
}

std::optional<float>
TruncatedDistanceMatrixDistance(const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
                                uint32_t i, uint32_t j) {
  CheckDistanceMatrixIndex(dm, i);
  CheckDistanceMatrixIndex(dm, j);
  return dm.Distance(i, j);
}

float
TruncatedDistanceMatrixDistanceOrDefault(
    const truncated_distance_matrix::TruncatedDistanceMatrix& dm, uint32_t i, uint32_t j) {
  CheckDistanceMatrixIndex(dm, i);
  CheckDistanceMatrixIndex(dm, j);
  return dm.DistanceOrDefault(i, j);
}

std::vector<float>
TruncatedDistanceMatrixDistancesOrDefault(
    const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
    const std::vector<uint32_t>& i, const std::vector<uint32_t>& j) {
  if (i.size() != j.size()) {
    throw std::invalid_argument("index arrays must have the same length");
  }
  for (uint32_t ndx : i) {
    CheckDistanceMatrixIndex(dm, ndx);
  }
  for (uint32_t ndx : j) {
    CheckDistanceMatrixIndex(dm, ndx);
  }
  return dm.DistancesOrDefault(i, j);
}

void
SetTruncatedDistanceMatrixDefaultDistance(
    truncated_distance_matrix::TruncatedDistanceMatrix& dm, float distance) {
  if (!dm.SetDefaultDistance(distance)) {
    throw std::invalid_argument("default distance cannot be below max stored distance");
  }
}

void
SetTruncatedDistanceMatrixDefaultDistanceByte(
    truncated_distance_matrix::TruncatedDistanceMatrix& dm, uint8_t distance) {
  if (!dm.SetDefaultDistanceByte(distance)) {
    throw std::invalid_argument("default distance cannot be below max stored distance");
  }
}

std::unordered_map<std::string, uint32_t>
DiceMolecule(dicer_api::Dicer& dicer, Molecule& mol) {
  std::unordered_map<std::string, uint32_t> result;
  dicer.Dice(mol, result);
  return result;
}

std::vector<Molecule>
MedchemWizardProcess(medchemwizard::MedchemWizard& wizard, const Molecule& mol) {
  Molecule mcopy(mol);
  resizable_array_p<Molecule> products;
  if (!wizard.ProcessToArray(mcopy, products)) {
    throw std::runtime_error("MedchemWizard product generation failed");
  }

  std::vector<Molecule> result;
  result.reserve(products.number_elements());
  for (Molecule* product : products) {
    result.emplace_back(*product);
  }

  return result;
}

std::unordered_map<std::string, uint64_t>
MedchemWizardStatsToMap(const medchemwizard::Stats& stats) {
  return {
      {"molecules_read", stats.molecules_read},
      {"molecules_produced", stats.molecules_produced},
      {"truncated_to_max_hits", stats.truncated_to_max_hits},
      {"duplicate_molecules_suppressed", stats.duplicate_molecules_suppressed},
      {"bad_valences_discarded", stats.bad_valences_discarded},
      {"discarded_for_too_many_atoms", stats.discarded_for_too_many_atoms},
      {"discarded_for_too_few_atoms", stats.discarded_for_too_few_atoms},
      {"molecules_not_matching_do_not_change_queries",
       stats.molecules_not_matching_do_not_change_queries},
      {"embeddings_rejected_for_changing_protected_atoms",
       stats.embeddings_rejected_for_changing_protected_atoms},
  };
}

struct GFPFactory {};

std::string
GFPTagToString(const IWString& tag) {
  return tag.AsString();
}

void
CheckGFPListIndex(const gfp_context::GFPList& gfp, int ndx, const char* argname) {
  if (ndx < 0 || ndx >= gfp.size()) {
    throw std::out_of_range(std::string(argname) + " index out of range");
  }
}

void
CheckGFPListMetadata(const gfp_context::GFPList& gfp) {
  if (!gfp.metadata_stored()) {
    throw std::runtime_error("GFPList does not store smiles/id metadata");
  }
}

void
CheckCompatibleFingerprint(const gfp_context::GFPList& gfp,
                           const gfp_context::GFPFingerprint& fingerprint) {
  if (fingerprint.context_hash() != gfp.context().context_hash()) {
    throw std::invalid_argument(
        "GFPFingerprint was generated with an incompatible GFPContext");
  }
}

std::vector<std::string>
GFPGeneratorSpecComponents(const gfp_context::GFPGeneratorSpec& spec) {
  std::vector<std::string> result;
  for (const gfp_context::Component& component : spec.Components()) {
    result.push_back(component.tag.AsString());
  }
  return result;
}

std::shared_ptr<gfp_context::GFPContext>
StandardGFPContext(bool preprocess) {
  auto result = std::make_shared<gfp_context::GFPContext>();
  if (!result->BuildStandard(preprocess)) {
    throw std::runtime_error("Cannot initialise standard GFP context");
  }
  return result;
}

std::shared_ptr<gfp_context::GFPContext>
GFPContextFromSpecs(const std::vector<gfp_context::GFPGeneratorSpec>& specs,
                    bool preprocess) {
  auto result = std::make_shared<gfp_context::GFPContext>();
  if (!result->BuildFromSpecs(specs, preprocess)) {
    throw std::runtime_error("Cannot initialise GFP context from specs");
  }
  return result;
}

std::unique_ptr<gfp_context::GFPFingerprint>
GFPContextFingerprint(gfp_context::GFPContext& context, Molecule& mol) {
  auto result = std::make_unique<gfp_context::GFPFingerprint>();
  if (!context.Fingerprint(mol, *result)) {
    throw std::runtime_error("Cannot generate GFP fingerprint");
  }
  return result;
}

void
GFPContextSetWeight(gfp_context::GFPContext& context, const std::string& tag,
                    float weight) {
  if (!context.SetWeight(IWString(tag), weight)) {
    throw std::runtime_error("Cannot set GFP weight for tag '" + tag + "'");
  }
}

void
GFPContextUseOnly(gfp_context::GFPContext& context,
                  const std::vector<std::string>& tags) {
  std::vector<IWString> iwtags;
  iwtags.reserve(tags.size());
  for (const std::string& tag : tags) {
    iwtags.emplace_back(tag);
  }
  if (!context.UseOnly(iwtags)) {
    throw std::runtime_error("Cannot restrict GFP components");
  }
}

std::shared_ptr<gfp_context::GFPList>
StandardGFPList(bool preprocess) {
  auto result = gfp_context::GFPList::Standard(preprocess);
  if (result == nullptr) {
    throw std::runtime_error("Cannot initialise standard GFP list");
  }
  return result;
}

std::shared_ptr<gfp_context::GFPList>
StandardGFPListFromMolecules(const std::vector<Molecule*>& molecules, bool preprocess,
                             bool store_metadata) {
  auto result = gfp_context::GFPList::StandardFromMolecules(
      molecules, preprocess, store_metadata);
  if (result == nullptr) {
    throw std::runtime_error("Cannot build standard GFP list from molecules");
  }
  return result;
}

std::shared_ptr<gfp_context::GFPList>
GFPListFromFile(const std::string& fname, int size_hint) {
  auto result = std::make_shared<gfp_context::GFPList>();
  if (!result->ReadFile(fname.c_str(), size_hint)) {
    throw std::runtime_error("Cannot read GFP file '" + fname + "'");
  }
  return result;
}

void
GFPListReadFile(gfp_context::GFPList& gfp, const std::string& fname, int size_hint) {
  if (!gfp.ReadFile(fname.c_str(), size_hint)) {
    throw std::runtime_error("Cannot read GFP file '" + fname + "'");
  }
}

std::string
GFPListSmiles(const gfp_context::GFPList& gfp, int ndx) {
  CheckGFPListIndex(gfp, ndx, "i");
  CheckGFPListMetadata(gfp);
  return GFPTagToString(gfp.smiles(ndx));
}

std::string
GFPListId(const gfp_context::GFPList& gfp, int ndx) {
  CheckGFPListIndex(gfp, ndx, "i");
  CheckGFPListMetadata(gfp);
  return GFPTagToString(gfp.id(ndx));
}

void
GFPListAdd(gfp_context::GFPList& gfp, Molecule& mol) {
  if (!gfp.Add(mol)) {
    throw std::runtime_error("Cannot add molecule to GFPList");
  }
}

void
GFPListAddMolecules(gfp_context::GFPList& gfp, const std::vector<Molecule*>& molecules,
                    bool store_metadata) {
  if (!gfp.AddMolecules(molecules, store_metadata)) {
    throw std::runtime_error("Cannot add molecules to GFPList");
  }
}

void
GFPListAddSmiles(gfp_context::GFPList& gfp, const std::vector<std::string>& smiles) {
  std::vector<std::unique_ptr<Molecule>> molecules;
  molecules.reserve(smiles.size());
  std::vector<Molecule*> molecule_ptrs;
  molecule_ptrs.reserve(smiles.size());

  for (size_t i = 0; i < smiles.size(); ++i) {
    auto molecule = std::make_unique<Molecule>();
    if (!molecule->build_from_smiles(smiles[i])) {
      throw std::runtime_error("Cannot parse smiles at index " + std::to_string(i) +
                               ": '" + smiles[i] + "'");
    }
    molecule_ptrs.push_back(molecule.get());
    molecules.push_back(std::move(molecule));
  }

  if (!gfp.AddMolecules(molecule_ptrs, false)) {
    throw std::runtime_error("Cannot add smiles to GFPList");
  }
}

float
GFPListDistanceIndices(const gfp_context::GFPList& gfp, int i, int j) {
  CheckGFPListIndex(gfp, i, "i");
  CheckGFPListIndex(gfp, j, "j");
  return gfp.Distance(i, j);
}

float
GFPListDistanceFingerprint(const gfp_context::GFPList& gfp,
                           const gfp_context::GFPFingerprint& fingerprint, int j) {
  CheckCompatibleFingerprint(gfp, fingerprint);
  CheckGFPListIndex(gfp, j, "j");
  const float result = gfp.Distance(fingerprint, j);
  if (!std::isfinite(result)) {
    throw std::runtime_error("GFP distance calculation failed");
  }
  return result;
}

std::vector<gfp_context::NearestNeighbour>
GFPListNearestNeighboursIndex(const gfp_context::GFPList& gfp, int query, int k) {
  CheckGFPListIndex(gfp, query, "query");
  return gfp.NearestNeighbours(query, k);
}

std::vector<gfp_context::NearestNeighbour>
GFPListNearestNeighboursFingerprint(const gfp_context::GFPList& gfp,
                                    const gfp_context::GFPFingerprint& query, int k) {
  CheckCompatibleFingerprint(gfp, query);
  return gfp.NearestNeighbours(query, k);
}

std::vector<gfp_context::NearestNeighbour>
GFPListNearestNeighboursWithinDistanceIndex(const gfp_context::GFPList& gfp,
                                            int query, float max_distance) {
  CheckGFPListIndex(gfp, query, "query");
  if (max_distance < 0.0f) {
    throw std::invalid_argument("max_distance must be non-negative");
  }
  return gfp.NearestNeighboursWithinDistance(query, max_distance);
}

std::vector<gfp_context::NearestNeighbour>
GFPListNearestNeighboursWithinDistanceFingerprint(
    const gfp_context::GFPList& gfp, const gfp_context::GFPFingerprint& query,
    float max_distance) {
  CheckCompatibleFingerprint(gfp, query);
  if (max_distance < 0.0f) {
    throw std::invalid_argument("max_distance must be non-negative");
  }
  return gfp.NearestNeighboursWithinDistance(query, max_distance);
}

void
GFPListSetWeight(gfp_context::GFPList& gfp, const std::string& tag, float weight) {
  if (!gfp.mutable_context().SetWeight(IWString(tag), weight)) {
    throw std::runtime_error("Cannot set GFP weight for tag '" + tag + "'");
  }
}

void
GFPListUseOnly(gfp_context::GFPList& gfp, const std::vector<std::string>& tags) {
  std::vector<IWString> iwtags;
  iwtags.reserve(tags.size());
  for (const std::string& tag : tags) {
    iwtags.emplace_back(tag);
  }
  if (!gfp.mutable_context().UseOnly(iwtags)) {
    throw std::runtime_error("Cannot restrict GFP components");
  }
}

std::string
NearestNeighbourRepr(const gfp_context::NearestNeighbour& neighbour) {
  return "GFPNearestNeighbour(index=" + std::to_string(neighbour.index) +
         ", distance=" + std::to_string(neighbour.distance) + ")";
}

}  // namespace

void
BindTools(nb::module_& m) {
  nb::enum_<truncated_distance_matrix::Storage>(m, "TruncatedDistanceMatrixStorage")
      .value("ROW_SPARSE", truncated_distance_matrix::Storage::kRowSparse)
      .value("ROW_HASH", truncated_distance_matrix::Storage::kRowHash);

  nb::enum_<truncated_distance_matrix::ProtoType>(m, "TruncatedDistanceMatrixProto")
      .value("NEARNEIGHBOURS", truncated_distance_matrix::ProtoType::kNearNeighbours)
      .value("NEARNEIGHBOURS_INDICES",
             truncated_distance_matrix::ProtoType::kNearNeighboursIndices);

  nb::class_<truncated_distance_matrix::TruncatedDistanceMatrix>(
      m, "TruncatedDistanceMatrix")
      .def("__init__",
           [](truncated_distance_matrix::TruncatedDistanceMatrix* dm,
              const std::string& fname, truncated_distance_matrix::Storage storage,
              truncated_distance_matrix::ProtoType proto_type) {
             new (dm) truncated_distance_matrix::TruncatedDistanceMatrix();
             if (!dm->Build(fname, storage, proto_type)) {
               throw std::runtime_error("Cannot build TruncatedDistanceMatrix from " + fname);
             }
           },
           nb::arg("fname"),
           nb::arg("storage") = truncated_distance_matrix::Storage::kRowSparse,
           nb::arg("proto_type") = truncated_distance_matrix::ProtoType::kNearNeighbours)
      .def("size", &truncated_distance_matrix::TruncatedDistanceMatrix::size)
      .def("__len__", &truncated_distance_matrix::TruncatedDistanceMatrix::size)
      .def("number_distances",
           &truncated_distance_matrix::TruncatedDistanceMatrix::number_distances)
      .def("duplicate_distances_differing_by_one",
           &truncated_distance_matrix::TruncatedDistanceMatrix::
               duplicate_distances_differing_by_one)
      .def("max_stored_distance",
           &truncated_distance_matrix::TruncatedDistanceMatrix::MaxStoredDistance)
      .def("default_distance",
           &truncated_distance_matrix::TruncatedDistanceMatrix::DefaultDistance)
      .def("max_stored_distance_byte",
           &truncated_distance_matrix::TruncatedDistanceMatrix::MaxStoredDistanceByte)
      .def("default_distance_byte",
           &truncated_distance_matrix::TruncatedDistanceMatrix::DefaultDistanceByte)
      .def("set_default_distance", &SetTruncatedDistanceMatrixDefaultDistance,
           nb::arg("distance"))
      .def("set_default_distance_byte", &SetTruncatedDistanceMatrixDefaultDistanceByte,
           nb::arg("distance"))
      .def("index",
           [](const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
              const std::string& name) { return dm.Index(name); },
           nb::arg("name"))
      .def("name", &TruncatedDistanceMatrixName, nb::arg("index"))
      .def("distance", &TruncatedDistanceMatrixDistance, nb::arg("i"), nb::arg("j"))
      .def("distance_or_default", &TruncatedDistanceMatrixDistanceOrDefault,
           nb::arg("i"), nb::arg("j"))
      .def("distances_or_default", &TruncatedDistanceMatrixDistancesOrDefault,
           nb::arg("i"), nb::arg("j"));

  nb::class_<gfp_context::GFPFingerprint>(m, "GFPFingerprint")
      .def("context_hash", &gfp_context::GFPFingerprint::context_hash);

  nb::class_<gfp_context::GFPGeneratorSpec>(m, "GFPGeneratorSpec")
      .def("components", &GFPGeneratorSpecComponents)
      .def("__repr__", &gfp_context::GFPGeneratorSpec::Repr);

  nb::class_<GFPFactory>(m, "GFP")
      .def_static("mpr", &gfp_context::GFPGeneratorSpec::MolecularProperties)
      .def_static("iw", &gfp_context::GFPGeneratorSpec::IWMFingerprint)
      .def_static("maccs", &gfp_context::GFPGeneratorSpec::MACCSKeys,
                  nb::arg("level2") = true)
      .def_static("formula", &gfp_context::GFPGeneratorSpec::FormulaFingerprint)
      .def_static("cats",
           [](int max_path_length, bool include_hydrophobic_pairs) {
             if (max_path_length < 1) {
               throw std::invalid_argument("max_path_length must be positive");
             }
             return gfp_context::GFPGeneratorSpec::CATS(max_path_length,
                                                        include_hydrophobic_pairs);
           },
           nb::arg("max_path_length") = 10,
           nb::arg("include_hydrophobic_pairs") = true)
      .def_static("alogp",
           [](int replicates) {
             if (replicates <= 0) {
               throw std::invalid_argument("replicates must be positive");
             }
             return gfp_context::GFPGeneratorSpec::ALogP(replicates);
           },
           nb::arg("replicates") = 9)
      .def_static("xlogp",
           [](int replicates) {
             if (replicates <= 0) {
               throw std::invalid_argument("replicates must be positive");
             }
             return gfp_context::GFPGeneratorSpec::XLogP(replicates);
           },
           nb::arg("replicates") = 9)
      .def_static("tpsa",
           [](int replicates) {
             if (replicates <= 0) {
               throw std::invalid_argument("replicates must be positive");
             }
             return gfp_context::GFPGeneratorSpec::TPSA(replicates);
           },
           nb::arg("replicates") = 9)
      .def_static("atom_pair",
           [](int min_separation, int max_separation, const std::string& atom_type,
              bool include_out_of_range) {
             if (min_separation < 0) {
               throw std::invalid_argument("min_separation must be non-negative");
             }
             if (max_separation < min_separation) {
               throw std::invalid_argument("max_separation must be >= min_separation");
             }
             if (atom_type.empty()) {
               throw std::invalid_argument("atom_type must be non-empty");
             }
             return gfp_context::GFPGeneratorSpec::AtomPair(
                 min_separation, max_separation, IWString(atom_type),
                 include_out_of_range);
           },
           nb::arg("min_separation") = 1, nb::arg("max_separation") = 10,
           nb::arg("atom_type") = "UST:Y", nb::arg("include_out_of_range") = false)
      .def_static("ec",
           [](int radius, const std::string& atom_type) {
             if (radius < 0) {
               throw std::invalid_argument("radius must be non-negative");
             }
             if (atom_type.empty()) {
               throw std::invalid_argument("atom_type must be non-empty");
             }
             return gfp_context::GFPGeneratorSpec::ECFingerprint(radius,
                                                                 IWString(atom_type));
           },
           nb::arg("radius") = 3, nb::arg("atom_type") = "UST:Z")
      .def_static("ring_substitution", &gfp_context::GFPGeneratorSpec::RingSubstitution)
      .def_static("spinach", &gfp_context::GFPGeneratorSpec::SpinachFingerprint,
                  nb::arg("label_join_points") = false)
      .def_static("scaffold", &gfp_context::GFPGeneratorSpec::ScaffoldFingerprint,
                  nb::arg("label_join_points") = false)
      .def_static("substructure",
           [](const std::string& smarts, int radius, const std::string& atom_type,
              const std::string& no_match) {
             if (smarts.empty()) {
               throw std::invalid_argument("smarts must be non-empty");
             }
             if (radius < 0) {
               throw std::invalid_argument("radius must be non-negative");
             }
             if (atom_type.empty()) {
               throw std::invalid_argument("atom_type must be non-empty");
             }
             bool no_match_is_empty;
             if (no_match == "empty") {
               no_match_is_empty = true;
             } else if (no_match == "error") {
               no_match_is_empty = false;
             } else {
               throw std::invalid_argument("no_match must be 'empty' or 'error'");
             }
             return gfp_context::GFPGeneratorSpec::SubstructureFingerprint(
                 IWString(smarts), radius, IWString(atom_type), no_match_is_empty);
           },
           nb::arg("smarts"), nb::arg("radius") = 0,
           nb::arg("atom_type") = "UST:ARY", nb::arg("no_match") = "empty");

  nb::class_<gfp_context::GFPContext>(m, "GFPContext")
      .def(nb::init<>())
      .def_static("standard", &StandardGFPContext, nb::arg("preprocess") = true,
                  "Create a context that generates the standard LillyMol GFP fingerprint")
      .def_static("from_specs", &GFPContextFromSpecs, nb::arg("specs"),
                  nb::arg("preprocess") = true,
                  "Create a context from GFP generator specifications")
      .def("tags", &gfp_context::GFPContext::Tags)
      .def("can_generate_fingerprints",
           &gfp_context::GFPContext::can_generate_fingerprints)
      .def("fingerprint", &GFPContextFingerprint, nb::arg("mol"))
      .def("distance", &gfp_context::GFPContext::Distance,
           nb::arg("lhs"), nb::arg("rhs"))
      .def("set_weight", &GFPContextSetWeight, nb::arg("tag"), nb::arg("weight"))
      .def("use_only", &GFPContextUseOnly, nb::arg("tags"))
      .def("use_all", &gfp_context::GFPContext::UseAll);

  nb::class_<gfp_context::NearestNeighbour>(m, "GFPNearestNeighbour")
      .def_ro("index", &gfp_context::NearestNeighbour::index)
      .def_ro("distance", &gfp_context::NearestNeighbour::distance)
      .def("__repr__", &NearestNeighbourRepr);

  nb::class_<gfp_context::GFPList>(m, "GFPList")
      .def(nb::init<>())
      .def(nb::init<std::shared_ptr<gfp_context::GFPContext>>(), nb::arg("context"))
      .def_static("standard", &StandardGFPList, nb::arg("preprocess") = true,
                  "Create an empty GFPList that generates standard LillyMol GFP fingerprints")
      .def_static("standard_from_molecules", &StandardGFPListFromMolecules,
                  nb::arg("molecules"), nb::arg("preprocess") = true,
                  nb::arg("store_metadata") = false,
                  "Build a standard GFPList from molecules")
      .def_static("from_file", &GFPListFromFile, nb::arg("fname"),
                  nb::arg("size_hint") = 0, "Read a GFP/TDT fingerprint file")
      .def("read_file", &GFPListReadFile, nb::arg("fname"),
           nb::arg("size_hint") = 0,
           "Read a GFP/TDT fingerprint file into this object")
      .def("__len__", &gfp_context::GFPList::size)
      .def("size", &gfp_context::GFPList::size)
      .def("metadata_stored", &gfp_context::GFPList::metadata_stored)
      .def("tags", [](const gfp_context::GFPList& gfp) {
        return gfp.context().Tags();
      })
      .def("smiles", &GFPListSmiles, nb::arg("i"))
      .def("id", &GFPListId, nb::arg("i"))
      .def("add", &GFPListAdd, nb::arg("mol"))
      .def("add_molecules", &GFPListAddMolecules, nb::arg("molecules"),
           nb::arg("store_metadata") = false)
      .def("add_smiles", &GFPListAddSmiles, nb::arg("smiles"),
           "Add fingerprints generated from a list of SMILES strings without storing metadata")
      .def("distance", &GFPListDistanceIndices, nb::arg("i"), nb::arg("j"))
      .def("distance", &GFPListDistanceFingerprint, nb::arg("fp"), nb::arg("j"))
      .def("nearest_neighbours", &GFPListNearestNeighboursIndex,
           nb::arg("query"), nb::arg("k"))
      .def("nearest_neighbours", &GFPListNearestNeighboursFingerprint,
           nb::arg("query"), nb::arg("k"))
      .def("nearest_neighbours_within_distance",
           &GFPListNearestNeighboursWithinDistanceIndex,
           nb::arg("query"), nb::arg("max_distance"))
      .def("nearest_neighbours_within_distance",
           &GFPListNearestNeighboursWithinDistanceFingerprint,
           nb::arg("query"), nb::arg("max_distance"))
      .def("set_weight", &GFPListSetWeight, nb::arg("tag"), nb::arg("weight"))
      .def("use_only", &GFPListUseOnly, nb::arg("tags"))
      .def("use_all", [](gfp_context::GFPList& gfp) {
        gfp.mutable_context().UseAll();
      });

  nb::class_<Mol2Graph>(m, "Mol2Graph")
      .def(nb::init<>())
      .def("set_exclude_triple_bonds_from_graph_reduction",
           &Mol2Graph::set_exclude_triple_bonds_from_graph_reduction,
           nb::arg("value"))
      .def("set_revert_all_directional_bonds_to_non_directional",
           &Mol2Graph::set_revert_all_directional_bonds_to_non_directional,
           nb::arg("value"))
      .def("set_preserve_cc_double_bonds_no_heteroatoms",
           &Mol2Graph::set_preserve_cc_double_bonds_no_heteroatoms,
           nb::arg("value"))
      .def("set_preserve_cc_double_bonds_saturated",
           &Mol2Graph::set_preserve_cc_double_bonds_saturated,
           nb::arg("value"))
      .def("set_append_molecular_formula", &Mol2Graph::set_append_molecular_formula,
           nb::arg("value"))
      .def("set_aromatic_distinguishing_formula",
           &Mol2Graph::set_aromatic_distinguishing_formula, nb::arg("value"))
      .def("set_remove_chiral_centres", &Mol2Graph::set_remove_chiral_centres,
           nb::arg("value"))
      .def("turn_on_most_useful_options", &Mol2Graph::TurnOnMostUsefulOptions)
      .def("set_active", &Mol2Graph::set_active, nb::arg("value"))
      .def("active", &Mol2Graph::active);

  nb::class_<dicer_api::Dicer>(m, "Dicer")
      .def(nb::init<>())
      .def("set_max_bonds_to_break", &dicer_api::Dicer::set_max_bonds_to_break,
           nb::arg("value"), "Set maximum number of bonds to break")
      .def("set_min_fragment_size", &dicer_api::Dicer::set_min_fragment_size,
           nb::arg("value"), "Set minimum fragment size")
      .def("set_max_fragment_size", &dicer_api::Dicer::set_max_fragment_size,
           nb::arg("value"), "Set maximum fragment size")
      .def("set_break_cc_bonds", &dicer_api::Dicer::set_break_cc_bonds,
           nb::arg("value"), "Control whether C-C bonds are broken")
      .def("set_break_cc_bonds_at_highly_connected",
           &dicer_api::Dicer::set_break_cc_bonds_at_highly_connected,
           nb::arg("value"), "Control whether C-[CD>2] bonds are broken")
      .def("set_break_amide_bonds", &dicer_api::Dicer::set_break_amide_bonds,
           nb::arg("value"), "Control whether amide bonds are broken")
      .def("set_label_join_points", &dicer_api::Dicer::set_label_join_points,
           nb::arg("value"), "Set isotope label for join points")
      .def("set_increment_isotope_for_join_points",
           &dicer_api::Dicer::set_increment_isotope_for_join_points,
           nb::arg("value"), "Increment isotopes at join points")
      .def("set_accumulate_global_fragment_count",
           &dicer_api::Dicer::set_accumulate_global_fragment_count,
           nb::arg("value"), "Accumulate fragments across dice() calls")
      .def("get_global_fragment_count", &dicer_api::Dicer::global_fragment_count,
           "Return accumulated global fragment counts")
      .def("set_perceive_symmetry_equivalent_matches",
           &dicer_api::Dicer::set_perceive_symmetry_equivalent_matches,
           nb::arg("value"), "Control symmetry-equivalent query matches")
      .def("set_determine_fragment_counts",
           &dicer_api::Dicer::set_determine_fragment_counts, nb::arg("value"),
           "Control whether per-molecule fragment counts are determined")
      .def("set_work_like_recap", &dicer_api::Dicer::set_work_like_recap,
           nb::arg("value"), "Work like Recap without recursion")
      .def("set_atom_type", &dicer_api::Dicer::set_atom_type, nb::arg("specification"),
           "Set atom typing used by dicer")
      .def("add_bond_break_smarts", &dicer_api::Dicer::AddBreakBondSmarts,
           nb::arg("smarts"), "Add bond-breaking SMARTS")
      .def("add_bond_break_query", &dicer_api::Dicer::AddBreakBondQuery,
           nb::arg("fname"), "Add bond-breaking query from a file directive")
      .def("add_fragment_requirement_smarts",
           &dicer_api::Dicer::AddFragmentRequirementSmarts, nb::arg("smarts"),
           "Require generated fragments to match SMARTS")
      .def("add_fragment_disqualifier_smarts",
           &dicer_api::Dicer::AddFragmentDisqualifierSmarts, nb::arg("smarts"),
           "Discard generated fragments matching SMARTS")
      .def("dice", &DiceMolecule, nb::arg("mol"),
           "Dice a molecule and return fragment unique smiles mapped to counts");

  nb::class_<medchemwizard::MedchemWizard>(m, "MedchemWizard")
      .def(nb::init<>())
      .def("initialise_from_environment",
           [](medchemwizard::MedchemWizard& wizard) {
             if (!wizard.InitialiseFromEnvironment()) {
               throw std::runtime_error(
                   "Cannot initialise MedchemWizard from "
                   "LILLYMOL_HOME/data/MedchemWizard/REACTIONS");
             }
           },
           "Initialise reactions from LILLYMOL_HOME/data/MedchemWizard/REACTIONS")
      .def("read_reactions",
           [](medchemwizard::MedchemWizard& wizard, const std::string& fname) {
             if (!wizard.ReadReactions(fname.c_str())) {
               throw std::runtime_error("Cannot read MedchemWizard reactions from '" +
                                        fname + "'");
             }
           },
           nb::arg("fname"), "Read reactions from a MedchemWizard REACTIONS file")
      .def("read_reaction",
           [](medchemwizard::MedchemWizard& wizard, const std::string& fname,
              bool textproto) {
             if (!wizard.ReadReactionFile(fname.c_str(), textproto)) {
               throw std::runtime_error("Cannot read MedchemWizard reaction from '" +
                                        fname + "'");
             }
           },
           nb::arg("fname"), nb::arg("textproto") = true,
           "Append a reaction file; resolve resources beside it. Use textproto=False for MSI.")
      .def("number_reactions", &medchemwizard::MedchemWizard::number_reactions)
      .def("process", &MedchemWizardProcess, nb::arg("mol"),
           "Generate MedchemWizard products for a molecule")
      .def("add_do_not_change_smarts",
           [](medchemwizard::MedchemWizard& wizard, const std::string& smarts) {
             auto query = std::make_unique<Substructure_Query>();
             if (!query->create_from_smarts(IWString(smarts))) {
               throw std::invalid_argument("Invalid SMARTS '" + smarts + "'");
             }
             wizard.do_not_change_queries() << query.release();
           },
           nb::arg("smarts"), "Protect atoms matching a SMARTS from reaction changes")
      .def("add_do_not_change_query",
           [](medchemwizard::MedchemWizard& wizard, const std::string& fname) {
             auto query = std::make_unique<Substructure_Query>();
             if (!query->read(fname.c_str())) {
               throw std::runtime_error("Cannot read do-not-change query '" + fname +
                                        "'");
             }
             wizard.do_not_change_queries() << query.release();
           },
           nb::arg("fname"), "Protect atoms matching a substructure query file")
      .def("set_ignore_do_not_change_queries_not_matching",
           [](medchemwizard::MedchemWizard& wizard, bool value) {
             wizard.options().ignore_do_not_change_queries_not_matching = value;
           },
           nb::arg("value"))
      .def("set_max_depth",
           [](medchemwizard::MedchemWizard& wizard, int value) {
             if (value < 0) {
               throw std::invalid_argument("max_depth must be non-negative");
             }
             wizard.options().max_depth = value;
           },
           nb::arg("value"))
      .def("set_max_hits",
           [](medchemwizard::MedchemWizard& wizard, int value) {
             if (value < 1) {
               throw std::invalid_argument("max_hits must be positive");
             }
             wizard.options().max_hits = value;
           },
           nb::arg("value"))
      .def("set_max_atoms",
           [](medchemwizard::MedchemWizard& wizard, int value) {
             if (value < 1) {
               throw std::invalid_argument("max_atoms must be positive");
             }
             wizard.options().max_atoms = value;
           },
           nb::arg("value"))
      .def("set_min_atoms",
           [](medchemwizard::MedchemWizard& wizard, int value) {
             if (value < 1) {
               throw std::invalid_argument("min_atoms must be positive");
             }
             wizard.options().min_atoms = value;
           },
           nb::arg("value"))
      .def("set_unique_within_molecule",
           [](medchemwizard::MedchemWizard& wizard, bool value) {
             wizard.options().unique_within_molecule = value;
           },
           nb::arg("value"))
      .def("set_unique_across_all_molecules",
           [](medchemwizard::MedchemWizard& wizard, bool value) {
             wizard.options().unique_across_all_molecules = value;
           },
           nb::arg("value"))
      .def("set_append_names",
           [](medchemwizard::MedchemWizard& wizard, bool value) {
             wizard.options().append_names = value;
           },
           nb::arg("value"))
      .def("set_name_separator",
           [](medchemwizard::MedchemWizard& wizard, const std::string& sep) {
             wizard.options().sep = sep;
           },
           nb::arg("sep"))
      .def("set_discard_bad_valences",
           [](medchemwizard::MedchemWizard& wizard, bool value) {
             wizard.options().discard_bad_valences = value;
           },
           nb::arg("value"))
      .def("stats",
           [](const medchemwizard::MedchemWizard& wizard) {
             return MedchemWizardStatsToMap(wizard.stats());
           },
           "Return MedchemWizard processing counters");

  nb::class_<ring_replacement::RingReplacement>(m, "RingReplacement")
      .def(nb::init<>())
      .def("set_ring_atom_smarts",
           &ring_replacement::RingReplacement::set_ring_atom_smarts,
           nb::arg("smarts"),
           "Set SMARTS for an atom in ring systems to be replaced")
      .def("set_unique_molecules_only",
           &ring_replacement::RingReplacement::set_unique_molecules_only,
           nb::arg("value"), "Control duplicate product suppression")
      .def("clear_unique_molecule_cache",
           &ring_replacement::RingReplacement::clear_unique_molecule_cache,
           "Clear products retained for duplicate suppression")
      .def("set_min_support_requirement",
           &ring_replacement::RingReplacement::set_min_support_requirement,
           nb::arg("value"),
           "Set minimum number of examples needed for a replacement ring")
      .def("set_max_formula_difference",
           &ring_replacement::RingReplacement::set_max_formula_difference,
           nb::arg("value"),
           "Set maximum formula difference between removed and replacement rings")
      .def("set_remove_isotopes",
           &ring_replacement::RingReplacement::set_remove_isotopes,
           nb::arg("value"), "Control whether isotopic labels are removed from products")
      .def("set_preserve_connectivity",
           &ring_replacement::RingReplacement::set_preserve_connectivity,
           nb::arg("value"),
           "Control whether replacement rings must preserve connectivity; set before "
           "reading replacement rings")
      .def("read_replacement_rings",
           [](ring_replacement::RingReplacement& replacement,
              const std::string& fname) -> uint32_t {
             return replacement.ReadReplacementRings(fname);
           },
           nb::arg("fname"), "Read replacement rings from a textproto file")
      .def("number_replacement_rings",
           &ring_replacement::RingReplacement::number_replacement_rings,
           "Return the number of replacement rings read")
      .def("process",
           [](ring_replacement::RingReplacement& replacement,
              Molecule& mol) -> std::vector<Molecule> {
             return replacement.Process(mol);
           },
           nb::arg("mol"), "Replace matching ring systems");

  nb::class_<unique_molecules::UniqueMolecules>(m, "UniqueMolecules")
      .def(nb::init<>())
      .def("set_include_chiral_info",
           &unique_molecules::UniqueMolecules::set_include_chiral_info,
           nb::arg("value"), "Control whether chirality is considered")
      .def("set_include_cis_trans_bonding_info",
           &unique_molecules::UniqueMolecules::set_include_cis_trans_bonding_info,
           nb::arg("value"), "Control whether cis/trans bonds are considered")
      .def("set_strip_to_largest_fragment",
           &unique_molecules::UniqueMolecules::set_strip_to_largest_fragment,
           nb::arg("value"), "Strip to largest fragment before comparison")
      .def("set_consider_isotopes",
           &unique_molecules::UniqueMolecules::set_consider_isotopes,
           nb::arg("value"), "Control whether isotopes are considered")
      .def("set_constant_isotope",
           &unique_molecules::UniqueMolecules::set_constant_isotope,
           nb::arg("value"), "Convert all isotope labels to a constant before comparison")
      .def("set_standardize_molecules",
           &unique_molecules::UniqueMolecules::set_standardize_molecules,
           nb::arg("value"), "Control whether molecules are standardised")
      .def("element_transformations",
           [](unique_molecules::UniqueMolecules& unique) -> Element_Transformations& {
             return unique.element_transformations();
           },
           nb::rv_policy::reference_internal, "Element transformations")
      .def("graph_specifications",
           [](unique_molecules::UniqueMolecules& unique) -> Mol2Graph& {
             return unique.graph_specifications();
           },
           nb::rv_policy::reference_internal, "Tautomer matching graph specifications")
      .def("add_to_hash", &unique_molecules::UniqueMolecules::AddToHashes,
           nb::arg("mol"), "Add molecule to internal uniqueness hashes")
      .def("is_unique", &unique_molecules::UniqueMolecules::IsUnique, nb::arg("mol"),
           "Return True if the molecule has not been seen before")
      .def("report", [](const unique_molecules::UniqueMolecules& unique) {
        return unique.Report(std::cerr);
      }, "Write a report to stderr");
}

}  // namespace lillymol_nb
