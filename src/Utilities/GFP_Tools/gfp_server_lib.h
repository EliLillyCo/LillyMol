#ifndef UTILITIES_GFP_TOOLS_GFP_SERVER_LIB_H_
#define UTILITIES_GFP_TOOLS_GFP_SERVER_LIB_H_

#include <cstdint>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/standardise.h"
#include "Molecule_Tools/maccskeys_fn5.h"
#include "Molecule_Tools/mpr.h"
#include "Utilities/GFP_Tools/gfp_standard.h"
#include "Utilities/GFP_Tools/nn_request.pb.h"

class iwstring_data_source;
class Molecule;
class IW_TDT;
class IW_General_Fingerprint;

namespace gfp_server {

// Transport-neutral near-neighbour search engine.
//
// The class owns a read-only in-memory fingerprint pool after Build(). Search()
// converts an incoming SMILES to the same standard fingerprint representation,
// searches the pool, and returns an NnReply proto. Networking, HTTP, protobuf
// JSON conversion, and process lifecycle are intentionally outside this class.
class GfpNearNeighbourServer {
 public:
  GfpNearNeighbourServer();

  // `pool_size_hint` is optional. If not positive, the input is scanned to count
  // fingerprint records before loading. If `store_smiles` is true, any $SMI<...>
  // value present in the input TDT is copied to neighbour result records.
  int Build(const char* fname, int pool_size_hint = 0, bool store_smiles = true);
  int Build(iwstring_data_source& input, int pool_size_hint = 0, bool store_smiles = true);

  int pool_size() const {
    return static_cast<int>(_pool.size());
  }

  // Search one request or a batch. Failures are represented in the returned
  // proto status/message rather than by returning std::optional or throwing.
  NnReply Search(const NnRequest& request);
  NnBatchReply Search(const NnBatchRequest& request);

  // Convenience API for pybind and other transports that want to pass serialized
  // protobuf bytes rather than C++ protobuf objects across the boundary.
  std::string SearchSerialized(std::string_view serialized_request);
  std::string SearchBatchSerialized(std::string_view serialized_request);

 private:
  struct NeighbourDistance {
    float distance;
    int index;
  };

  struct NeighbourSimilarity {
    float similarity;
    int index;
  };

  int BuildGfpStandard(GFP_Standard& sfp, IW_General_Fingerprint& gfp) const;
  int BuildGfpStandard(GFP_Standard& sfp, IWString& id, IW_TDT& tdt) const;

  void Preprocess(Molecule& m);
  int SmilesToGfp(const std::string& smiles, GFP_Standard& gfp);

  void FillResultProto(const std::vector<NeighbourDistance>& did, nnbr::NearNeighbours& result) const;
  int TanimotoSingleNbr(const GFP_Standard& gfp, nnbr::NearNeighbours& result) const;
  int TanimotoWithinDistance(const GFP_Standard& gfp, float cutoff, nnbr::NearNeighbours& result) const;
  int Tanimoto(const GFP_Standard& gfp, int nbrs, nnbr::NearNeighbours& result) const;

  std::vector<GFP_Standard> _pool;
  std::vector<std::string> _id;
  std::vector<std::string> _smiles;

  Molecular_Properties_Generator _mpr;
  MACCSKeys _mk;
  Chemical_Standardisation _chemical_standardisation;

  // SmilesToGfp uses mutable helper objects above. Keep that mutation isolated
  // so multiple HTTP worker threads cannot race during query fingerprinting.
  std::mutex _fingerprint_mutex;
};

}  // namespace gfp_server

#endif  // UTILITIES_GFP_TOOLS_GFP_SERVER_LIB_H_
