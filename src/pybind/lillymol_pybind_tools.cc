// Bindings for selected tools

#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/dicer_api.h"
#include "Molecule_Tools/iwdescr_lib.h"
#include "Molecule_Tools/jwcats_lib.h"
#include "Molecule_Tools/medchemwizard_lib.h"
#include "Molecule_Tools/mformula.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/qed.h"
#include "Molecule_Tools/ring_replacement_lib.h"
#include "Molecule_Tools/unique_molecules_api.h"
#include "Molecule_Tools/xlogp.h"
#include "Molecule_Tools_Bdb/iwecfp_database_lookup_lib.h"
#include "Molecule_Tools_Bdb/selimsteg.h"
#include "Molecule_Tools_Bdb/structure_database.h"
#include "Utilities/GFP_Tools/gfp_context.h"
#include "Utilities/GFP_Tools/truncated_distance_matrix.h"
#include "pybind11/stl.h"

namespace py = pybind11;

namespace {

struct GFPFactory {};

std::string
IWStringToStdString(const IWString& s) {
  return std::string(s.rawchars(), s.length());
}

void
CheckGFPListIndex(const gfp_context::GFPList& gfp, int i, const char* argname) {
  if (i < 0 || i >= gfp.size()) {
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
                           const gfp_context::GFPFingerprint& fp) {
  if (fp.context_hash() != gfp.context().context_hash()) {
    throw std::invalid_argument(
        "GFPFingerprint was generated with an incompatible GFPContext");
  }
}

// The counts live in Molecule_Lib. This used to be a local loop that capped
// the donor contribution of a nitrogen at two and scored every oxygen as one,
// which agreed with neither Lipinski nor the other copies in the tree.
std::tuple<int, int>
LipinskiHbaHbd(Molecule& m) {
  int hba = 0;
  int hbd = 0;
  m.LipinskiHbaHbd(hba, hbd);

  return std::make_tuple(hba, hbd);
}

class MolecularDescriptors {
  private:
    IWDescr _iwdescr;

    py::list FeatureNamesAsPyList() const {
      py::list result;
      const int ndescr = _iwdescr.number_descriptors();
      for (int i = 0; i < ndescr; ++i) {
        result.append(_iwdescr.descriptor_name(i).AsString());
      }
      return result;
    }

  public:
    MolecularDescriptors() {
      if (!_iwdescr.InitialiseAll()) {
        throw std::runtime_error(
            "Cannot initialise MolecularDescriptors; ensure LILLYMOL_HOME is defined "
            "and contains the standard charge and donor/acceptor queries");
      }
    }

    std::vector<std::string> Names() const {
      std::vector<std::string> result;
      const int ndescr = _iwdescr.number_descriptors();
      result.reserve(ndescr);
      for (int i = 0; i < ndescr; ++i) {
        result.push_back(_iwdescr.descriptor_name(i).AsString());
      }
      return result;
    }

    py::array_t<float> ComputeArray(Molecule& mol) {
      py::array_t<float> result(_iwdescr.number_descriptors());
      if (!_iwdescr.Process(mol, result.mutable_data())) {
        throw std::runtime_error("MolecularDescriptors calculation failed");
      }
      return result;
    }

    py::dict Compute(Molecule& mol) {
      py::array_t<float> values = ComputeArray(mol);
      const float* ptr = values.data();
      py::dict result;
      const int ndescr = _iwdescr.number_descriptors();
      for (int i = 0; i < ndescr; ++i) {
        result[py::str(_iwdescr.descriptor_name(i).AsString())] = ptr[i];
      }
      return result;
    }

    py::object ComputeList(std::vector<Molecule*>& mols, bool as_dataframe) {
      const int nmols = static_cast<int>(mols.size());
      const int ndescr = _iwdescr.number_descriptors();

      py::array_t<float> result({nmols, ndescr});
      float* ptr = result.mutable_data();

      for (int i = 0; i < nmols; ++i) {
        if (!_iwdescr.Process(*mols[i], ptr + i * ndescr)) {
          throw std::runtime_error(
              std::string("MolecularDescriptors calculation failed for molecule ") +
              mols[i]->name().AsString());
        }
      }

      if (!as_dataframe) {
        return py::object(result);
      }

      py::list index;
      for (int i = 0; i < nmols; ++i) {
        index.append(mols[i]->name().AsString());
      }

      py::module_ pd = py::module_::import("pandas");
      return pd.attr("DataFrame")(result, py::arg("columns") = FeatureNamesAsPyList(),
                                   py::arg("index") = index);
    }
};

std::unique_ptr<jwcats::JWCats>
MakeJWCats(bool initialise_default_assigners) {
  auto result = std::make_unique<jwcats::JWCats>();

  if (initialise_default_assigners) {
    if (!result->charge_assigner().BuildFromDefaultEnvs()) {
      throw std::runtime_error(
          "Cannot initialise JWCats charge assigner; ensure LILLYMOL_HOME is defined");
    }
    IWString default_dir("DEF");
    if (!result->donor_acceptor_assigner().BuildFromDir(default_dir, 0)) {
      throw std::runtime_error(
          "Cannot initialise JWCats donor/acceptor assigner; ensure LILLYMOL_HOME is "
          "defined");
    }
  }

  if (!result->Initialise()) {
    throw std::runtime_error("Cannot initialise JWCats");
  }

  return result;
}

py::array_t<double>
JWCatsResultToArray(const jwcats::JWCats& jwcats, const jwcats::Result& result) {
  const std::vector<int>& write_array_value = jwcats.write_array_value();
  int nfeatures = 0;
  for (int value : write_array_value) {
    if (value) {
      ++nfeatures;
    }
  }

  py::array_t<double> array(nfeatures);
  double* ptr = array.mutable_data();

  int ndx = 0;
  for (int i = 0; i < static_cast<int>(write_array_value.size()); ++i) {
    if (write_array_value[i]) {
      ptr[ndx] = result.scaled_counts[i];
      ++ndx;
    }
  }

  return array;
}

void
ThrowForJWCatsStatus(jwcats::ComputeStatus status) {
  switch (status) {
    case jwcats::ComputeStatus::kOk:
      return;
    case jwcats::ComputeStatus::kMissingChargeData:
      throw std::runtime_error("JWCats calculation failed: missing charge data");
    case jwcats::ComputeStatus::kNotInitialised:
      throw std::runtime_error("JWCats calculation failed: object is not initialised");
    case jwcats::ComputeStatus::kError:
      throw std::runtime_error("JWCats calculation failed");
  }

  throw std::runtime_error("JWCats calculation failed: unknown status");
}

py::dict
MedchemWizardStatsToDict(const medchemwizard::Stats& stats) {
  py::dict result;
  result["molecules_read"] = stats.molecules_read;
  result["molecules_produced"] = stats.molecules_produced;
  result["truncated_to_max_hits"] = stats.truncated_to_max_hits;
  result["duplicate_molecules_suppressed"] = stats.duplicate_molecules_suppressed;
  result["bad_valences_discarded"] = stats.bad_valences_discarded;
  result["discarded_for_too_many_atoms"] = stats.discarded_for_too_many_atoms;
  result["discarded_for_too_few_atoms"] = stats.discarded_for_too_few_atoms;
  result["molecules_not_matching_do_not_change_queries"] =
      stats.molecules_not_matching_do_not_change_queries;
  result["embeddings_rejected_for_changing_protected_atoms"] =
      stats.embeddings_rejected_for_changing_protected_atoms;

  return result;
}

}  // namespace

PYBIND11_MODULE(lillymol_tools, m) {
  using gfp_context::GFPContext;
  using gfp_context::GFPFingerprint;
  using gfp_context::GFPGeneratorSpec;
  using gfp_context::GFPList;
  using gfp_context::NearestNeighbour;
  using unique_molecules::UniqueMolecules;
  using truncated_distance_matrix::ProtoType;
  using truncated_distance_matrix::Storage;
  using truncated_distance_matrix::TruncatedDistanceMatrix;

  py::enum_<Storage>(m, "TruncatedDistanceMatrixStorage")
      .value("ROW_SPARSE", Storage::kRowSparse)
      .value("ROW_HASH", Storage::kRowHash)
      .export_values();

  py::enum_<ProtoType>(m, "TruncatedDistanceMatrixProto")
      .value("NEARNEIGHBOURS", ProtoType::kNearNeighbours)
      .value("NEARNEIGHBOURS_INDICES", ProtoType::kNearNeighboursIndices)
      .export_values();

  py::class_<TruncatedDistanceMatrix>(m, "TruncatedDistanceMatrix")
      .def(py::init([](const std::string& fname, Storage storage, ProtoType proto_type) {
             auto result = std::make_unique<TruncatedDistanceMatrix>();
             if (!result->Build(fname, storage, proto_type)) {
               throw std::runtime_error("Cannot build TruncatedDistanceMatrix from " + fname);
             }
             return result;
           }),
           py::arg("fname"), py::arg("storage") = Storage::kRowSparse,
           py::arg("proto_type") = ProtoType::kNearNeighbours)
      .def("size", &TruncatedDistanceMatrix::size)
      .def("__len__", &TruncatedDistanceMatrix::size)
      .def("number_distances", &TruncatedDistanceMatrix::number_distances)
      .def("duplicate_distances_differing_by_one",
           &TruncatedDistanceMatrix::duplicate_distances_differing_by_one)
      .def("max_stored_distance", &TruncatedDistanceMatrix::MaxStoredDistance)
      .def("default_distance", &TruncatedDistanceMatrix::DefaultDistance)
      .def("max_stored_distance_byte", &TruncatedDistanceMatrix::MaxStoredDistanceByte)
      .def("default_distance_byte", &TruncatedDistanceMatrix::DefaultDistanceByte)
      .def("set_default_distance",
           [](TruncatedDistanceMatrix& dm, float d) {
             if (!dm.SetDefaultDistance(d)) {
               throw std::invalid_argument("default distance cannot be below max stored distance");
             }
           },
           py::arg("distance"))
      .def("set_default_distance_byte",
           [](TruncatedDistanceMatrix& dm, uint8_t d) {
             if (!dm.SetDefaultDistanceByte(d)) {
               throw std::invalid_argument("default distance cannot be below max stored distance");
             }
           },
           py::arg("distance"))
      .def("index", &TruncatedDistanceMatrix::Index, py::arg("name"))
      .def("name",
           [](const TruncatedDistanceMatrix& dm, uint32_t ndx) -> std::string {
             if (ndx >= dm.size()) {
               throw py::index_error("index out of range");
             }
             return dm.Name(ndx);
           },
           py::arg("index"))
      .def("distance",
           [](const TruncatedDistanceMatrix& dm, uint32_t i, uint32_t j) {
             if (i >= dm.size() || j >= dm.size()) {
               throw py::index_error("index out of range");
             }
             return dm.Distance(i, j);
           },
           py::arg("i"), py::arg("j"))
      .def("distance_or_default",
           [](const TruncatedDistanceMatrix& dm, uint32_t i, uint32_t j) {
             if (i >= dm.size() || j >= dm.size()) {
               throw py::index_error("index out of range");
             }
             return dm.DistanceOrDefault(i, j);
           },
           py::arg("i"), py::arg("j"))
      .def("distances_or_default",
           [](const TruncatedDistanceMatrix& dm, const std::vector<uint32_t>& i,
              const std::vector<uint32_t>& j) {
             if (i.size() != j.size()) {
               throw std::invalid_argument("index arrays must have the same length");
             }
             for (uint32_t ndx : i) {
               if (ndx >= dm.size()) {
                 throw py::index_error("index out of range");
               }
             }
             for (uint32_t ndx : j) {
               if (ndx >= dm.size()) {
                 throw py::index_error("index out of range");
               }
             }
             return dm.DistancesOrDefault(i, j);
           },
           py::arg("i"), py::arg("j"));

  py::class_<GFPFingerprint>(m, "GFPFingerprint")
      .def("context_hash", &GFPFingerprint::context_hash);

  py::class_<GFPGeneratorSpec>(m, "GFPGeneratorSpec")
      .def("components",
           [](const GFPGeneratorSpec& spec) {
             std::vector<std::string> result;
             for (const gfp_context::Component& component : spec.Components()) {
               result.emplace_back(component.tag.rawchars(), component.tag.length());
             }
             return result;
           })
      .def("__repr__", &GFPGeneratorSpec::Repr);

  py::class_<GFPFactory>(m, "GFP")
      .def_static("mpr", &GFPGeneratorSpec::MolecularProperties)
      .def_static("iw", &GFPGeneratorSpec::IWMFingerprint)
      .def_static("maccs", &GFPGeneratorSpec::MACCSKeys, py::arg("level2") = true)
      .def_static("formula", &GFPGeneratorSpec::FormulaFingerprint)
      .def_static(
          "cats",
          [](int max_path_length, bool include_hydrophobic_pairs) {
            if (max_path_length < 1) {
              throw std::invalid_argument("max_path_length must be positive");
            }
            return GFPGeneratorSpec::CATS(max_path_length, include_hydrophobic_pairs);
          },
          py::arg("max_path_length") = 10, py::arg("include_hydrophobic_pairs") = true)
      .def_static(
          "alogp",
          [](int replicates) {
            if (replicates <= 0) {
              throw std::invalid_argument("replicates must be positive");
            }
            return GFPGeneratorSpec::ALogP(replicates);
          },
          py::arg("replicates") = 9)
      .def_static(
          "xlogp",
          [](int replicates) {
            if (replicates <= 0) {
              throw std::invalid_argument("replicates must be positive");
            }
            return GFPGeneratorSpec::XLogP(replicates);
          },
          py::arg("replicates") = 9)
      .def_static(
          "tpsa",
          [](int replicates) {
            if (replicates <= 0) {
              throw std::invalid_argument("replicates must be positive");
            }
            return GFPGeneratorSpec::TPSA(replicates);
          },
          py::arg("replicates") = 9)
      .def_static(
          "atom_pair",
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
            return GFPGeneratorSpec::AtomPair(min_separation, max_separation,
                                              IWString(atom_type), include_out_of_range);
          },
          py::arg("min_separation") = 1, py::arg("max_separation") = 10,
          py::arg("atom_type") = "UST:Y", py::arg("include_out_of_range") = false)
      .def_static(
          "ec",
          [](int radius, const std::string& atom_type) {
            if (radius < 0) {
              throw std::invalid_argument("radius must be non-negative");
            }
            if (atom_type.empty()) {
              throw std::invalid_argument("atom_type must be non-empty");
            }
            return GFPGeneratorSpec::ECFingerprint(radius, IWString(atom_type));
          },
          py::arg("radius") = 3, py::arg("atom_type") = "UST:Z")
      .def_static("ring_substitution", &GFPGeneratorSpec::RingSubstitution)
      .def_static("spinach", &GFPGeneratorSpec::SpinachFingerprint,
                  py::arg("label_join_points") = false)
      .def_static("scaffold", &GFPGeneratorSpec::ScaffoldFingerprint,
                  py::arg("label_join_points") = false)
      .def_static(
          "substructure",
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
            return GFPGeneratorSpec::SubstructureFingerprint(
                IWString(smarts), radius, IWString(atom_type), no_match_is_empty);
          },
          py::arg("smarts"), py::arg("radius") = 0,
          py::arg("atom_type") = "UST:ARY", py::arg("no_match") = "empty");

  py::class_<GFPContext, std::shared_ptr<GFPContext>>(m, "GFPContext")
      .def(py::init<>())
      .def_static(
          "standard",
          [](bool preprocess) {
            auto result = std::make_shared<GFPContext>();
            if (!result->BuildStandard(preprocess)) {
              throw std::runtime_error("Cannot initialise standard GFP context");
            }
            return result;
          },
          py::arg("preprocess") = true,
          "Create a context that generates the standard LillyMol GFP fingerprint")
      .def_static(
          "from_specs",
          [](const std::vector<GFPGeneratorSpec>& specs, bool preprocess) {
            auto result = std::make_shared<GFPContext>();
            if (!result->BuildFromSpecs(specs, preprocess)) {
              throw std::runtime_error("Cannot initialise GFP context from specs");
            }
            return result;
          },
          py::arg("specs"), py::arg("preprocess") = true,
          "Create a context from GFP generator specifications")
      .def("tags", &GFPContext::Tags)
      .def("can_generate_fingerprints", &GFPContext::can_generate_fingerprints)
      .def(
          "fingerprint",
          [](GFPContext& context, Molecule& mol) {
            auto result = std::make_unique<GFPFingerprint>();
            if (!context.Fingerprint(mol, *result)) {
              throw std::runtime_error("Cannot generate GFP fingerprint");
            }
            return result;
          },
          py::arg("mol"))
      .def("distance", &GFPContext::Distance, py::arg("lhs"), py::arg("rhs"))
      .def(
          "set_weight",
          [](GFPContext& context, const std::string& tag, float weight) {
            IWString iwtag(tag);
            if (!context.SetWeight(iwtag, weight)) {
              throw std::runtime_error("Cannot set GFP weight for tag '" + tag + "'");
            }
          },
          py::arg("tag"), py::arg("weight"))
      .def(
          "use_only",
          [](GFPContext& context, const std::vector<std::string>& tags) {
            std::vector<IWString> iwtags;
            iwtags.reserve(tags.size());
            for (const std::string& tag : tags) {
              iwtags.emplace_back(tag);
            }
            if (!context.UseOnly(iwtags)) {
              throw std::runtime_error("Cannot restrict GFP components");
            }
          },
          py::arg("tags"))
      .def("use_all", &GFPContext::UseAll);

  py::class_<NearestNeighbour>(m, "GFPNearestNeighbour")
      .def_readonly("index", &NearestNeighbour::index)
      .def_readonly("distance", &NearestNeighbour::distance)
      .def("__repr__", [](const NearestNeighbour& n) {
        return "GFPNearestNeighbour(index=" + std::to_string(n.index) +
               ", distance=" + std::to_string(n.distance) + ")";
      });

  py::class_<GFPList, std::shared_ptr<GFPList>>(m, "GFPList")
      .def(py::init<>())
      .def(py::init<std::shared_ptr<GFPContext>>(), py::arg("context"))
      .def_static(
          "standard",
          [](bool preprocess) {
            auto result = GFPList::Standard(preprocess);
            if (result == nullptr) {
              throw std::runtime_error("Cannot initialise standard GFP list");
            }
            return result;
          },
          py::arg("preprocess") = true,
          "Create an empty GFPList that generates standard LillyMol GFP fingerprints")
      .def_static(
          "standard_from_molecules",
          [](const std::vector<Molecule*>& molecules, bool preprocess,
             bool store_metadata) {
            auto result =
                GFPList::StandardFromMolecules(molecules, preprocess, store_metadata);
            if (result == nullptr) {
              throw std::runtime_error("Cannot build standard GFP list from molecules");
            }
            return result;
          },
          py::arg("molecules"), py::arg("preprocess") = true,
          py::arg("store_metadata") = false,
          "Build a standard GFPList from a Python sequence of Molecule objects")
      .def_static(
          "from_file",
          [](const std::string& fname, int size_hint) {
            auto result = std::make_shared<GFPList>();
            if (!result->ReadFile(fname.c_str(), size_hint)) {
              throw std::runtime_error("Cannot read GFP file '" + fname + "'");
            }
            return result;
          },
          py::arg("fname"), py::arg("size_hint") = 0, "Read a GFP/TDT fingerprint file")
      .def(
          "read_file",
          [](GFPList& gfp, const std::string& fname, int size_hint) {
            if (!gfp.ReadFile(fname.c_str(), size_hint)) {
              throw std::runtime_error("Cannot read GFP file '" + fname + "'");
            }
          },
          py::arg("fname"), py::arg("size_hint") = 0,
          "Read a GFP/TDT fingerprint file into this object")
      .def("__len__", &GFPList::size)
      .def("size", &GFPList::size)
      .def("metadata_stored", &GFPList::metadata_stored)
      .def("tags", [](const GFPList& gfp) { return gfp.context().Tags(); })
      .def(
          "smiles",
          [](const GFPList& gfp, int i) {
            CheckGFPListIndex(gfp, i, "i");
            CheckGFPListMetadata(gfp);
            return IWStringToStdString(gfp.smiles(i));
          },
          py::arg("i"))
      .def(
          "id",
          [](const GFPList& gfp, int i) {
            CheckGFPListIndex(gfp, i, "i");
            CheckGFPListMetadata(gfp);
            return IWStringToStdString(gfp.id(i));
          },
          py::arg("i"))
      .def(
          "add",
          [](GFPList& gfp, Molecule& mol) {
            if (!gfp.Add(mol)) {
              throw std::runtime_error("Cannot add molecule to GFPList");
            }
          },
          py::arg("mol"))
      .def(
          "add_molecules",
          [](GFPList& gfp, const std::vector<Molecule*>& molecules, bool store_metadata) {
            if (!gfp.AddMolecules(molecules, store_metadata)) {
              throw std::runtime_error("Cannot add molecules to GFPList");
            }
          },
          py::arg("molecules"), py::arg("store_metadata") = false)
      .def(
          "distance",
          [](const GFPList& gfp, int i, int j) {
            CheckGFPListIndex(gfp, i, "i");
            CheckGFPListIndex(gfp, j, "j");
            return gfp.Distance(i, j);
          },
          py::arg("i"), py::arg("j"))
      .def(
          "distance",
          [](const GFPList& gfp, const GFPFingerprint& fp, int j) {
            CheckCompatibleFingerprint(gfp, fp);
            CheckGFPListIndex(gfp, j, "j");
            const float result = gfp.Distance(fp, j);
            if (!std::isfinite(result)) {
              throw std::runtime_error("GFP distance calculation failed");
            }
            return result;
          },
          py::arg("fp"), py::arg("j"))
      .def(
          "nearest_neighbours",
          [](const GFPList& gfp, int query, int k) {
            CheckGFPListIndex(gfp, query, "query");
            return gfp.NearestNeighbours(query, k);
          },
          py::arg("query"), py::arg("k"))
      .def(
          "nearest_neighbours",
          [](const GFPList& gfp, const GFPFingerprint& query, int k) {
            CheckCompatibleFingerprint(gfp, query);
            return gfp.NearestNeighbours(query, k);
          },
          py::arg("query"), py::arg("k"))
      .def(
          "nearest_neighbours_within_distance",
          [](const GFPList& gfp, int query, float max_distance) {
            CheckGFPListIndex(gfp, query, "query");
            if (max_distance < 0.0f) {
              throw std::invalid_argument("max_distance must be non-negative");
            }
            return gfp.NearestNeighboursWithinDistance(query, max_distance);
          },
          py::arg("query"), py::arg("max_distance"))
      .def(
          "nearest_neighbours_within_distance",
          [](const GFPList& gfp, const GFPFingerprint& query, float max_distance) {
            CheckCompatibleFingerprint(gfp, query);
            if (max_distance < 0.0f) {
              throw std::invalid_argument("max_distance must be non-negative");
            }
            return gfp.NearestNeighboursWithinDistance(query, max_distance);
          },
          py::arg("query"), py::arg("max_distance"))
      .def(
          "set_weight",
          [](GFPList& gfp, const std::string& tag, float weight) {
            IWString iwtag(tag);
            if (!gfp.mutable_context().SetWeight(iwtag, weight)) {
              throw std::runtime_error("Cannot set GFP weight for tag '" + tag + "'");
            }
          },
          py::arg("tag"), py::arg("weight"))
      .def(
          "use_only",
          [](GFPList& gfp, const std::vector<std::string>& tags) {
            std::vector<IWString> iwtags;
            iwtags.reserve(tags.size());
            for (const std::string& tag : tags) {
              iwtags.emplace_back(tag);
            }
            if (!gfp.mutable_context().UseOnly(iwtags)) {
              throw std::runtime_error("Cannot restrict GFP components");
            }
          },
          py::arg("tags"))
      .def("use_all", [](GFPList& gfp) { gfp.mutable_context().UseAll(); });

  // This is a sub-optimal implementation. While functional, it is not efficient.
  // Much to my surprise I found that storing smiles in a python set() was much
  // faster than storing those same strings in the C++ map used in the current
  // implementation. TODO:ianwatson understand what is going on.
  // For now this is quite usable, just not efficient.
  py::class_<nvrtspsa::NovartisPolarSurfaceArea>(m, "TPSA")
      .def("set_rdkit_compatibility",
           [](nvrtspsa::NovartisPolarSurfaceArea& tpsa) { tpsa.SetRDKitCompatibility(); },
           "Set the options that make this agree with RDKit's TPSA - charge separated "
           "form, and zero contribution for sulphur and phosphorus.")
      .def(py::init<>())
      .def(
          "compute",
          [](nvrtspsa::NovartisPolarSurfaceArea& tpsa, Molecule& mol)
              -> std::optional<double> { return tpsa.PolarSurfaceArea(mol); },
          py::arg("mol"),
          "Compute topological polar surface area for a molecule, returning None for "
          "empty molecules")
      .def(
          "set_display_psa_unclassified_atom_messages",
          &nvrtspsa::NovartisPolarSurfaceArea::set_display_psa_unclassified_atom_messages,
          py::arg("value"), "Control warnings for unclassified atoms")
      .def("set_return_zero_for_unclassified_atoms",
           &nvrtspsa::NovartisPolarSurfaceArea::set_return_zero_for_unclassified_atoms,
           py::arg("value"), "Return zero for unclassified atoms")
      .def("set_non_zero_contribution_for_SD2",
           &nvrtspsa::NovartisPolarSurfaceArea::set_non_zero_contribution_for_SD2,
           py::arg("value"), "Control the SD2 sulphur contribution")
      .def("set_zero_for_all_sulphur_atoms",
           &nvrtspsa::NovartisPolarSurfaceArea::set_zero_for_all_sulphur_atoms,
           py::arg("value"), "Set all sulphur atom TPSA contributions to zero")
      .def("set_zero_for_all_phosphorus_atoms",
           &nvrtspsa::NovartisPolarSurfaceArea::set_zero_for_all_phosphorus_atoms,
           py::arg("value"), "Set all phosphorus atom TPSA contributions to zero")
      .def(
          "set_convert_to_charge_separated",
          &nvrtspsa::NovartisPolarSurfaceArea::set_convert_to_charge_separated,
          py::arg("value"), "Convert to charge-separated forms before computing TPSA");

  m.def(
      "HbaHbd", [](Molecule& mol) { return LipinskiHbaHbd(mol); }, py::arg("mol"),
      "Return Lipinski-style hydrogen-bond acceptor and donor counts as (hba, hbd)");

  m.def(
      "molecular_weight",
      [](const Molecule& mol) { return lillymol::MolecularWeightIsotopesAsLabels(mol); },
      py::arg("mol"),
      "Return molecular weight, treating isotope labels as labels rather than masses");

  m.def(
      "alogp",
      [](Molecule& mol) -> std::optional<float> {
        alogp::ALogP calc;
        calc.set_use_alcohol_for_acid(1);
        calc.set_rdkit_phoshoric_acid_hydrogen(1);
        std::optional<double> result = calc.LogP(mol);
        if (!result) {
          return std::nullopt;
        }
        return static_cast<float>(*result);
      },
      py::arg("mol"), "Compute ALogP, returning None if the calculation fails");

  m.def(
      "xlogp",
      [](Molecule& mol) -> std::optional<float> {
        xlogp::XLogPCalc calc;
        std::optional<double> result = calc.LogP(mol);
        if (!result) {
          return std::nullopt;
        }
        return static_cast<float>(*result);
      },
      py::arg("mol"), "Compute XLogP, returning None if the calculation fails");

  m.def(
      "tpsa",
      [](Molecule& mol) -> std::optional<double> {
        nvrtspsa::NovartisPolarSurfaceArea calc;
        return calc.PolarSurfaceArea(mol);
      },
      py::arg("mol"), "Compute topological polar surface area");

  m.def(
      "qed_score",
      [](const Molecule& mol) -> std::optional<float> {
        qed::Qed calc;
        if (!calc.InitialiseFromEnvironment()) {
          throw std::runtime_error(
              "Cannot initialise QED; ensure LILLYMOL_HOME points to a LillyMol tree");
        }
        Molecule mcopy(mol);
        return calc.qed(mcopy);
      },
      py::arg("mol"), "Compute QED, returning None if the calculation fails");

  py::class_<qed::Qed>(m, "QED")
      .def(py::init([](bool initialise_from_environment) {
             auto result = std::make_unique<qed::Qed>();
             if (initialise_from_environment && ! result->InitialiseFromEnvironment()) {
               throw std::runtime_error(
                   "Cannot initialise QED; ensure LILLYMOL_HOME points to a LillyMol tree");
             }
             return result;
           }),
           py::arg("initialise_from_environment") = true)
      .def(
          "initialise_from_environment",
          [](qed::Qed& qed) { return qed.InitialiseFromEnvironment(); },
          "Initialise QED from LILLYMOL_HOME/data/queries/QED")
      .def(
          "initialise_from_directory",
          [](qed::Qed& qed, const std::string& dirname) {
            IWString d(dirname);
            return qed.InitialiseFromDirectory(d);
          },
          py::arg("dirname"),
          "Initialise QED from a directory containing the QED smarts files")
      .def(
          "qed",
          [](qed::Qed& qed, const Molecule& mol) -> std::optional<float> {
            Molecule mcopy(mol);
            return qed.qed(mcopy);
          },
          py::arg("mol"), "Compute QED, returning None if the calculation fails");

  py::class_<medchemwizard::MedchemWizard>(m, "MedchemWizard")
      .def(py::init<>())
      .def(
          "initialise_from_environment",
          [](medchemwizard::MedchemWizard& wizard) {
            if (!wizard.InitialiseFromEnvironment()) {
              throw std::runtime_error(
                  "Cannot initialise MedchemWizard from "
                  "LILLYMOL_HOME/data/MedchemWizard/REACTIONS");
            }
          },
          "Initialise reactions from LILLYMOL_HOME/data/MedchemWizard/REACTIONS")
      .def(
          "read_reactions",
          [](medchemwizard::MedchemWizard& wizard, const std::string& fname) {
            if (!wizard.ReadReactions(fname.c_str())) {
              throw std::runtime_error("Cannot read MedchemWizard reactions from '" +
                                       fname + "'");
            }
          },
          py::arg("fname"), "Read reactions from a MedchemWizard REACTIONS file")
      .def("number_reactions", &medchemwizard::MedchemWizard::number_reactions)
      .def(
          "process",
          [](medchemwizard::MedchemWizard& wizard, const Molecule& mol) {
            Molecule mcopy(mol);
            resizable_array_p<Molecule> products;
            if (!wizard.ProcessToArray(mcopy, products)) {
              throw std::runtime_error("MedchemWizard product generation failed");
            }

            std::vector<std::shared_ptr<Molecule>> result;
            result.reserve(products.number_elements());
            for (Molecule* product : products) {
              result.emplace_back(std::shared_ptr<Molecule>(product));
            }
            products.resize_no_delete(0);
            return result;
          },
          py::arg("mol"), "Generate MedchemWizard products for a molecule")
      .def(
          "add_do_not_change_smarts",
          [](medchemwizard::MedchemWizard& wizard, const std::string& smarts) {
            auto query = std::make_unique<Substructure_Query>();
            if (!query->create_from_smarts(IWString(smarts))) {
              throw std::invalid_argument("Invalid SMARTS '" + smarts + "'");
            }
            wizard.do_not_change_queries() << query.release();
          },
          py::arg("smarts"), "Protect atoms matching a SMARTS from reaction changes")
      .def(
          "add_do_not_change_query",
          [](medchemwizard::MedchemWizard& wizard, const std::string& fname) {
            auto query = std::make_unique<Substructure_Query>();
            if (!query->read(fname.c_str())) {
              throw std::runtime_error("Cannot read do-not-change query '" + fname +
                                       "'");
            }
            wizard.do_not_change_queries() << query.release();
          },
          py::arg("fname"), "Protect atoms matching a substructure query file")
      .def(
          "set_ignore_do_not_change_queries_not_matching",
          [](medchemwizard::MedchemWizard& wizard, bool value) {
            wizard.options().ignore_do_not_change_queries_not_matching = value;
          },
          py::arg("value"))
      .def(
          "set_max_depth",
          [](medchemwizard::MedchemWizard& wizard, int value) {
            if (value < 0) {
              throw std::invalid_argument("max_depth must be non-negative");
            }
            wizard.options().max_depth = value;
          },
          py::arg("value"))
      .def(
          "set_max_hits",
          [](medchemwizard::MedchemWizard& wizard, int value) {
            if (value < 1) {
              throw std::invalid_argument("max_hits must be positive");
            }
            wizard.options().max_hits = value;
          },
          py::arg("value"))
      .def(
          "set_max_atoms",
          [](medchemwizard::MedchemWizard& wizard, int value) {
            if (value < 1) {
              throw std::invalid_argument("max_atoms must be positive");
            }
            wizard.options().max_atoms = value;
          },
          py::arg("value"))
      .def(
          "set_min_atoms",
          [](medchemwizard::MedchemWizard& wizard, int value) {
            if (value < 1) {
              throw std::invalid_argument("min_atoms must be positive");
            }
            wizard.options().min_atoms = value;
          },
          py::arg("value"))
      .def(
          "set_unique_within_molecule",
          [](medchemwizard::MedchemWizard& wizard, bool value) {
            wizard.options().unique_within_molecule = value;
          },
          py::arg("value"))
      .def(
          "set_unique_across_all_molecules",
          [](medchemwizard::MedchemWizard& wizard, bool value) {
            wizard.options().unique_across_all_molecules = value;
          },
          py::arg("value"))
      .def(
          "set_append_names",
          [](medchemwizard::MedchemWizard& wizard, bool value) {
            wizard.options().append_names = value;
          },
          py::arg("value"))
      .def(
          "set_name_separator",
          [](medchemwizard::MedchemWizard& wizard, const std::string& sep) {
            wizard.options().sep = sep;
          },
          py::arg("sep"))
      .def(
          "set_discard_bad_valences",
          [](medchemwizard::MedchemWizard& wizard, bool value) {
            wizard.options().discard_bad_valences = value;
          },
          py::arg("value"))
      .def(
          "stats",
          [](const medchemwizard::MedchemWizard& wizard) {
            return MedchemWizardStatsToDict(wizard.stats());
          },
          "Return MedchemWizard processing counters");

  py::class_<unique_molecules::UniqueMolecules>(m, "UniqueMolecules")
      .def(py::init<>())
      .def("set_include_chiral_info", &UniqueMolecules::set_include_chiral_info,
           "Control whether chirality is considered")
      .def("set_include_cis_trans_bonding_info",
           &UniqueMolecules::set_include_cis_trans_bonding_info,
           "Control whether cis trans bonds are considered")
      .def("set_strip_to_largest_fragment",
           &UniqueMolecules::set_strip_to_largest_fragment, "Strip to largest fragment")
      .def("set_consider_isotopes", &UniqueMolecules::set_consider_isotopes,
           "Control whether isotopes are cosidered")
      .def("set_constant_isotope", &UniqueMolecules::set_constant_isotope,
           "Convert all isotopes to a constant")
      .def("set_standardize_molecules", &UniqueMolecules::set_standardize_molecules,
           "Control whether or not molecules are standardised")
      .def(
          "element_transformations",
          [](UniqueMolecules& u) -> Element_Transformations& {
            return u.element_transformations();
          },
          py::return_value_policy::reference, "Element transformations")
  // Element Transformations temporarily not working - need to figure out duplicate
  // Element ptr issue
#ifdef ADD_ELEMENT_TRANSFORMATION_NOT_YET_WORKING
      .def(
          "add_element_transformation",
          [](UniqueMolecules& u, const std::string& trans) -> bool {
            const IWString s(trans);
            return u.element_transformations().Add(s);
          },
          "Add one element transformation")
#endif
      .def(
          "graph_specifications",
          [](UniqueMolecules& u) -> Mol2Graph& { return u.graph_specifications(); },
          py::return_value_policy::reference, "Tautomer matching specifications")
      .def("add_to_hash", &UniqueMolecules::AddToHashes,
           "Add molecule to internal structures")
      .def("is_unique", &UniqueMolecules::IsUnique, "True if the molecule is unique")
      .def(
          "report", [](const UniqueMolecules& u) { u.Report(std::cerr); }, "Report");

  // #ifdef THIS_WILL_REQUIRE_BERKELEY_DB_DYNAMIC_LIBRARIES
  py::class_<iwecfp_database_lookup::SP_Set_of_Databases>(m,
                                                          "SyntheticPrecedentDatabases")
      .def(py::init<>())
      .def(
          "add_database",
          [](iwecfp_database_lookup::SP_Set_of_Databases& databases,
             const std::string& dbname) -> bool { return databases.AddDatabase(dbname); },
          "Add an existing synethetic precedent database")
      .def(
          "set_max_radius",
          [](iwecfp_database_lookup::SP_Set_of_Databases& databases,
             int max_radius) -> bool { return databases.set_max_radius(max_radius); },
          "Set max radius for shells")
      .def(
          "per_shell_data",
          [](iwecfp_database_lookup::SP_Set_of_Databases& databases,
             Molecule& m) -> std::vector<int> {
            std::vector<int> result(3);

            // Need to do something better here, but does it ever fail?
            if (!databases.PerShellData(m, result)) {
              std::cerr << "Failed\n";
            }
            return result;
          },
          "Report lowest bit prevalence at each radius")
      .def(
          "slurp",
          [](iwecfp_database_lookup::SP_Set_of_Databases& databases,
             uint32_t min_examples) -> bool { return databases.slurp(min_examples); },
          "Slurp database entries with min_examples or more examples to memory")
      .def("__repr__", [](const iwecfp_database_lookup::Set_of_Databases& databases) {
        IWString result;
        result << "Synethetic Precedent database with " << databases.number_databases()
               << " databases";
        return result.AsString();
      });
  // #endif

  py::class_<selimsteg::Selimsteg>(m, "Selimsteg")
      .def(py::init<>())
      .def("open_database", &selimsteg::Selimsteg::OpenDatabase,
           "Open a BerkeleyDB datbase with Id->smiles mappints")
      .def("get_smiles", &selimsteg::Selimsteg::Lookup,
           "Fetch the smiles for an identifier")
      .def("get_molecule", &selimsteg::Selimsteg::GetMolecule,
           "Fetch a Molecule for an identifier")
      .def("get_molecules", &selimsteg::Selimsteg::GetMolecules,
           "Fetch a list of Molecules for list of identifiers")
      .def("close", &selimsteg::Selimsteg::Close, "Close the database. Idempotent")
      // A database holds a file handle, so give it the same shape as Reader - use
      // it in a with block, or call close(). Relying on the object being
      // collected works but is not something a caller can see or control.
      // Returns the python object, not a Selimsteg&. Returning a reference lets
      // the default return value policy copy it, and then the with block operates
      // on a copy while __exit__ closes the original.
      .def("__enter__", [](py::object self) { return self; })
      .def("__exit__",
           [](selimsteg::Selimsteg& s, py::object, py::object, py::object)->bool {
             s.Close();
             return false;   // do not suppress an exception from the block
           });

  py::enum_<structure_database::Lookup>(m, "LookupParams", py::arithmetic())
      .value("EXACT", structure_database::kExact)
      .value("STRIP", structure_database::kStrip)
      .value("NOCHIRAL", structure_database::kNoChiral)
      .value("GRAPH", structure_database::kGraph)
      .value("NOSTD", structure_database::kNoStandardise)
      .export_values();
  ;
  m.def("value", [](structure_database::Lookup e) { return static_cast<int>(e); });

  py::class_<structure_database::StructureDatabase>(m, "StructureDatabase")
      .def(py::init<>())
      .def(
          "open_read",
          [](structure_database::StructureDatabase& db,
             const std::string& dbname) -> bool {
            IWString tmp(dbname);
            return db.OpenForReading(tmp);
          },
          "Open a structure database for reading")
      .def(
          "lookup",
          [](structure_database::StructureDatabase& db,
             Molecule& m) -> std::optional<std::string> {
            IWString tmp;
            const uint32_t mask = structure_database::Lookup::kExact;
            int rc = db.Lookup(m, mask, tmp);
            if (rc == 0) {
              return std::nullopt;
            }
            std::string result = tmp.AsString();
            return result;
          },
          "")
      .def(
          "lookup",
          [](structure_database::StructureDatabase& db, Molecule& m,
             const uint32_t params) -> std::optional<std::string> {
            IWString tmp;
            int rc = db.Lookup(m, params, tmp);
            if (rc == 0) {
              return std::nullopt;
            }
            std::string result = tmp.AsString();
            return result;
          },
          "Lookup molecule in database, returning id's of equivalent molecules");
  py::class_<dicer_api::Dicer>(m, "Dicer")
      .def(py::init<>())
      .def("set_max_bonds_to_break", &dicer_api::Dicer::set_max_bonds_to_break,
           "set max bonds to break")
      .def("set_min_fragment_size", &dicer_api::Dicer::set_min_fragment_size,
           "set min fragment size")
      .def("set_max_fragment_size", &dicer_api::Dicer::set_max_fragment_size,
           "set max fragment size")
      .def("set_break_cc_bonds", &dicer_api::Dicer::set_break_cc_bonds,
           "set True of C-C bonds are to be broken")
      .def("set_break_cc_bonds_at_highly_connected",
           &dicer_api::Dicer::set_break_cc_bonds_at_highly_connected,
           "C-[CD>2] bonds are to be broken")
      .def("set_label_join_points", &dicer_api::Dicer::set_label_join_points,
           "Isotope for join points")
      .def("set_increment_isotope_for_join_points",
           &dicer_api::Dicer::set_increment_isotope_for_join_points,
           "Isotopes for join points added to any existing isotopic value")
      .def("set_accumulate_global_fragment_count",
           &dicer_api::Dicer::set_accumulate_global_fragment_count,
           "Set to True to accumulate fragments")
      .def("get_global_fragment_count", &dicer_api::Dicer::global_fragment_count,
           "Fetch the global fragment count")
      .def("set_perceive_symmetry_equivalent_matches",
           &dicer_api::Dicer::set_perceive_symmetry_equivalent_matches,
           "Set 0 to NOT perceive symmetry equivalent matches")
      .def("set_determine_fragment_counts",
           &dicer_api::Dicer::set_determine_fragment_counts,
           "If True determine the number of times a fragment appears in each starting "
           "molecule")
      .def("set_work_like_recap", &dicer_api::Dicer::set_work_like_recap,
           "Work like Recap - no recursion")
      .def(
          "add_bond_break_smarts",
          [](dicer_api::Dicer& dicer, const std::string& smarts) {
            return dicer.AddBreakBondSmarts(smarts);
          },
          "Add bond breaking smarts")
      .def(
          "add_bond_break_query",
          [](dicer_api::Dicer& dicer, const std::string& fname) {
            return dicer.AddBreakBondQuery(fname);
          },
          "Add bond breaking query from textproto file")
      .def(
          "dice",
          [](dicer_api::Dicer& dicer,
             Molecule& m) -> std::unordered_map<std::string, uint32_t> {
            std::unordered_map<std::string, uint32_t> result;
            dicer.Dice(m, result);
            return result;
          },
          "dice the molecule and return fragments and counts");

  py::class_<ring_replacement::RingReplacement>(m, "RingReplacement")
      .def(py::init<>())
      .def(
          "set_ring_atom_smarts",
          [](ring_replacement::RingReplacement& rp, const std::string& smarts) -> bool {
            return rp.set_ring_atom_smarts(smarts);
          },
          "Smarts for an atom in ring/ring systems to be removed")
      .def("set_unique_molecules_only",
           &ring_replacement::RingReplacement::set_unique_molecules_only,
           "only unique products will be generated")
      .def("clear_unique_molecule_cache",
           &ring_replacement::RingReplacement::clear_unique_molecule_cache,
           "Clear products retained for duplicate suppression")
      .def("set_min_support_requirement",
           &ring_replacement::RingReplacement::set_min_support_requirement,
           "min number of examples needed for a replacement ring to be used")
      .def("set_max_formula_difference",
           &ring_replacement::RingReplacement::set_max_formula_difference,
           "specify a maximum difference between the formula of the starting ring and "
           "the replacement")
      .def("set_remove_isotopes", &ring_replacement::RingReplacement::set_remove_isotopes,
           "isotopic labels are removed from products")
      .def(
          "read_replacement_rings",
          [](ring_replacement::RingReplacement& rp, const std::string& fname)
              -> uint32_t { return rp.ReadReplacementRings(fname); },
          "Add a set of replacement rings, "
          "/path/to/lillymol/data/ring_replacement/rings_6a.smi")
      .def("number_replacement_rings",
           &ring_replacement::RingReplacement::number_replacement_rings,
           "Return the number of replacement rings read by read_replacement_rings")
      .def(
          "process",
          [](ring_replacement::RingReplacement& rp,
             Molecule& m) -> std::vector<Molecule> { return rp.Process(m); },
          "replace rings");

  py::class_<jwcats::JWCats>(m, "JWCats")
      .def(py::init([](bool initialise_default_assigners) {
             return MakeJWCats(initialise_default_assigners);
           }),
           py::arg("initialise_default_assigners") = true,
           "Create a JWCats descriptor calculator. By default the standard charge and "
           "donor/acceptor assigners are initialised from LILLYMOL_HOME.")
      .def(
          "initialise", [](jwcats::JWCats& jwcats) { return jwcats.Initialise(); },
          "Initialise after changing settings")
      .def(
          "build_default_assigners",
          [](jwcats::JWCats& jwcats) {
            IWString default_dir("DEF");
            return jwcats.charge_assigner().BuildFromDefaultEnvs() &&
                   jwcats.donor_acceptor_assigner().BuildFromDir(default_dir, 0);
          },
          "Initialise charge and donor/acceptor assigners from LILLYMOL_HOME")
      .def("feature_names", &jwcats::JWCats::FeatureNames,
           "Return descriptor names in the same order as process() values")
      .def(
          "process",
          [](jwcats::JWCats& jwcats, Molecule& mol) {
            jwcats::Result result;
            const jwcats::ComputeStatus status = jwcats.Compute(mol, result);
            ThrowForJWCatsStatus(status);
            return JWCatsResultToArray(jwcats, result);
          },
          py::arg("mol"),
          "Compute JWCats descriptors for one molecule as a float64 NumPy array. The "
          "molecule may be modified by descriptor calculation.")
      .def(
          "process_list",
          [](jwcats::JWCats& jwcats, std::vector<Molecule*>& mols) {
            const std::vector<std::string> names = jwcats.FeatureNames();
            const int nmols = static_cast<int>(mols.size());
            const int nfeatures = static_cast<int>(names.size());

            py::array_t<double> array({nmols, nfeatures});
            double* ptr = array.mutable_data();

            for (int i = 0; i < nmols; ++i) {
              jwcats::Result result;
              const jwcats::ComputeStatus status = jwcats.Compute(*mols[i], result);
              ThrowForJWCatsStatus(status);

              const std::vector<int>& write_array_value = jwcats.write_array_value();
              int col = 0;
              for (int j = 0; j < static_cast<int>(write_array_value.size()); ++j) {
                if (write_array_value[j]) {
                  ptr[i * nfeatures + col] = result.scaled_counts[j];
                  ++col;
                }
              }
            }

            return array;
          },
          py::arg("mols"),
          "Compute JWCats descriptors for a list of molecules as a float64 NumPy array")
      .def("set_include_hydrophobic_pairs", &jwcats::JWCats::SetIncludeHydrophobicPairs,
           py::arg("value"),
           "Control whether hydrophobe-hydrophobe pair columns are emitted")
      .def("set_min_bond_separation", &jwcats::JWCats::SetMinBondSeparation,
           py::arg("value"), "Set the minimum bond separation")
      .def("set_max_bond_separation", &jwcats::JWCats::SetMaxBondSeparation,
           py::arg("value"), "Set the maximum bond separation")
      .def("set_scaling_type", &jwcats::JWCats::SetScalingType, py::arg("value"),
           "Set scaling type: 0 none, 1 heavy atoms, 2 feature counts, 3 heavy atoms / "
           "feature counts")
      .def("set_make_implicit_hydrogens_explicit",
           &jwcats::JWCats::SetMakeImplicitHydrogensExplicit, py::arg("value"),
           "Make implicit hydrogens explicit during calculation")
      .def("initialised", &jwcats::JWCats::initialised,
           "Return whether the object has been initialised");

  py::class_<MolecularDescriptors>(m, "MolecularDescriptors")
      .def(py::init<>())
      .def("names", &MolecularDescriptors::Names,
           "Return descriptor names in the same order as compute_array() values")
      .def("feature_names", &MolecularDescriptors::Names,
           "Alias for names()")
      .def("compute_array", &MolecularDescriptors::ComputeArray, py::arg("mol"),
           "Compute all descriptors for one molecule as a float32 NumPy array")
      .def("compute", &MolecularDescriptors::Compute, py::arg("mol"),
           "Compute all descriptors for one molecule as a dict keyed by descriptor name")
      .def("compute_list", &MolecularDescriptors::ComputeList,
           py::arg("mols"), py::arg("as_dataframe") = false,
           "Compute all descriptors for a list of molecules as an array or DataFrame");

  py::class_<IWDescr>(m, "IWDescr")
      .def(py::init([]() {
        auto result = std::make_unique<IWDescr>();
        if (!result->InitialiseAll()) {
          throw std::runtime_error(
              "Cannot initialise IWDescr; ensure LILLYMOL_HOME is defined and contains "
              "the standard charge and donor/acceptor queries");
        }
        return result;
      }))
      .def(
          "feature_names",
          [](const IWDescr& iwdescr) {
            std::vector<std::string> result;
            result.reserve(iwdescr.number_descriptors());
            for (int i = 0; i < iwdescr.number_descriptors(); ++i) {
              result.push_back(iwdescr.descriptor_name(i).AsString());
            }
            return result;
          },
          "Return descriptor names in the same order as process() values")
      .def(
          "process",
          [](IWDescr& iwdescr, Molecule& mol) {
            py::array_t<float> result(iwdescr.number_descriptors());
            if (!iwdescr.Process(mol, result.mutable_data())) {
              throw std::runtime_error("IWDescr calculation failed");
            }
            return result;
          },
          "Compute all descriptors for one molecule as a float32 NumPy array. "
          "The molecule may be modified by descriptor calculation")
      .def(
          "process_list",
          [](IWDescr& iwdescr, std::vector<Molecule*>& mols, bool as_dataframe) {
            const int nmols = static_cast<int>(mols.size());
            const int ndescr = iwdescr.number_descriptors();

            py::array_t<float> result({nmols, ndescr});
            float* ptr = result.mutable_data();

            for (int i = 0; i < nmols; ++i) {
              if (!iwdescr.Process(*mols[i], ptr + i * ndescr)) {
                throw std::runtime_error(
                    std::string("IWDescr calculation failed for molecule ") +
                    mols[i]->name().AsString());
              }
            }

            if (!as_dataframe) {
              return py::object(result);
            }

            // Build column names from descriptor names
            py::list columns;
            for (int i = 0; i < ndescr; ++i) {
              columns.append(iwdescr.descriptor_name(i).AsString());
            }

            // Build row index from molecule names
            py::list index;
            for (int i = 0; i < nmols; ++i) {
              index.append(mols[i]->name().AsString());
            }

            py::module_ pd = py::module_::import("pandas");
            py::object df = pd.attr("DataFrame")(result, py::arg("columns") = columns,
                                                 py::arg("index") = index);
            return df;
          },
          py::arg("mols"), py::arg("as_dataframe") = false,
          "Compute all descriptors for a list of molecules. "
          "Returns a float32 NumPy array of shape (n_molecules, n_descriptors) by "
          "default, "
          "or a pandas DataFrame (with molecule names as index) if as_dataframe=True.")

      ;

  py::class_<mformula::MFormula>(m, "MFormula")
      .def(py::init<>())
      .def(
          "build", [](mformula::MFormula& mf, Molecule& m) -> int { return mf.Build(m); },
          "Build the molecular for`m`, returns number of atoms in `m`")
      .def(
          "set_consider_aromatic",
          [](mformula::MFormula mf, bool s) { mf.set_consider_aromatic(s); },
          "Controls whether aromatic atoms are distinguished - true by default");
}
