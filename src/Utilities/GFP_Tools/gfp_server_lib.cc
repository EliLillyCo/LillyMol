#include "Utilities/GFP_Tools/gfp_server_lib.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <queue>
#include <utility>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwbits/iwbits.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"

namespace gfp_server {
namespace {

IWString smiles_tag("$SMI<");

}  // namespace

GfpNearNeighbourServer::GfpNearNeighbourServer() {
  _chemical_standardisation.activate_all();
}

int
GfpNearNeighbourServer::Build(const char* fname, int pool_size_hint, bool store_smiles) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    return 0;
  }

  return Build(input, pool_size_hint, store_smiles);
}

int
GfpNearNeighbourServer::Build(iwstring_data_source& input, int pool_size_hint, bool store_smiles) {
  int pool_size = pool_size_hint;
  if (pool_size <= 0) {
    pool_size = input.count_records_starting_with("|");
  }

  if (pool_size <= 0) {
    return 0;
  }

  _pool.clear();
  _id.clear();
  _smiles.clear();

  _pool.reserve(pool_size);
  _id.reserve(pool_size);
  if (store_smiles) {
    _smiles.reserve(pool_size);
  }

  for (int i = 0; i < pool_size; ++i) {
    IW_TDT tdt;
    if (! tdt.next(input)) {
      break;
    }

    GFP_Standard fp;
    IWString id;
    if (! BuildGfpStandard(fp, id, tdt)) {
      return 0;
    }

    _pool.push_back(fp);
    _id.push_back(id.AsString());

    if (store_smiles) {
      IWString smiles;
      tdt.dataitem_value(smiles_tag, smiles);
      _smiles.push_back(smiles.AsString());
    }
  }

  return ! _pool.empty();
}

int
GfpNearNeighbourServer::BuildGfpStandard(GFP_Standard& sfp,
                                         IW_General_Fingerprint& gfp) const {
  if (! gfp.molecular_properties_integer().active()) {
    return 0;
  }
  if (gfp.nfingerprints() < 3) {
    return 0;
  }

  sfp.build_molecular_properties(gfp.molecular_properties_integer());
  sfp.build_iw(gfp[0]);
  sfp.build_mk(gfp[1]);
  sfp.build_mk2(gfp[2]);

  return 1;
}

int
GfpNearNeighbourServer::BuildGfpStandard(GFP_Standard& sfp,
                                         IWString& id,
                                         IW_TDT& tdt) const {
  int fatal;
  IW_General_Fingerprint gfp;

  if (! gfp.construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  id = gfp.id();

  return BuildGfpStandard(sfp, gfp);
}

void
GfpNearNeighbourServer::Preprocess(Molecule& m) {
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  _chemical_standardisation.process(m);
}

int
GfpNearNeighbourServer::SmilesToGfp(const std::string& smiles, GFP_Standard& gfp) {
  Molecule m;
  if (! m.build_from_smiles(smiles)) {
    return 0;
  }

  std::lock_guard<std::mutex> lock(_fingerprint_mutex);

  Preprocess(m);

  IWMFingerprint iwfp;
  iwfp.construct_fingerprint(m);
  gfp.build_iwfp(iwfp.bits(), iwfp.nset());

  int tmp[2048];
  std::fill_n(tmp, 2048, 0);

  _mk(m, tmp);
  gfp.build_mk(tmp, _mk.nbits());

  _mk.set_level_2_fingerprint(tmp);
  gfp.build_mk2(tmp, _mk.nbits());

  int molecular_properties[8];
  _mpr(m, molecular_properties);
  gfp.build_molecular_properties(molecular_properties, 8);

  return 1;
}

void
GfpNearNeighbourServer::FillResultProto(const std::vector<NeighbourDistance>& did,
                                        nnbr::NearNeighbours& result) const {
  result.mutable_nbr()->Reserve(did.size());

  const bool have_smiles = _smiles.size() == _pool.size();
  for (const NeighbourDistance& item : did) {
    nnbr::Nbr* nbr = result.add_nbr();
    nbr->set_id(_id[item.index]);
    nbr->set_dist(item.distance);
    if (have_smiles) {
      nbr->set_smi(_smiles[item.index]);
    }
  }
}

int
GfpNearNeighbourServer::TanimotoSingleNbr(const GFP_Standard& gfp,
                                          nnbr::NearNeighbours& result) const {
  if (_pool.empty()) {
    return 0;
  }

  float best_similarity = _pool[0].tanimoto(gfp);
  int best_index = 0;

  for (int i = 1; i < static_cast<int>(_pool.size()); ++i) {
    const float similarity = _pool[i].tanimoto(gfp);
    if (similarity > best_similarity) {
      best_similarity = similarity;
      best_index = i;
    }
  }

  std::vector<NeighbourDistance> did;
  did.push_back(NeighbourDistance{1.0f - best_similarity, best_index});
  FillResultProto(did, result);

  return 1;
}

int
GfpNearNeighbourServer::TanimotoWithinDistance(const GFP_Standard& gfp,
                                               float cutoff,
                                               nnbr::NearNeighbours& result) const {
  cutoff = 1.0f - cutoff;  // convert distance to similarity.

  std::vector<NeighbourDistance> did;
  did.reserve(_pool.size());
  for (int i = 0; i < static_cast<int>(_pool.size()); ++i) {
    const float similarity = _pool[i].tanimoto(gfp);
    if (similarity >= cutoff) {
      did.push_back(NeighbourDistance{1.0f - similarity, i});
    }
  }

  std::sort(did.begin(), did.end(), [](const NeighbourDistance& lhs, const NeighbourDistance& rhs) {
    if (lhs.distance != rhs.distance) {
      return lhs.distance < rhs.distance;
    }
    return lhs.index < rhs.index;
  });

  FillResultProto(did, result);
  return 1;
}

int
GfpNearNeighbourServer::Tanimoto(const GFP_Standard& gfp,
                                 int nbrs,
                                 nnbr::NearNeighbours& result) const {
  if (nbrs <= 0 || _pool.empty()) {
    return 0;
  }

  if (nbrs > static_cast<int>(_pool.size())) {
    nbrs = _pool.size();
  }

  std::vector<NeighbourSimilarity> seed;
  seed.reserve(nbrs);
  for (int i = 0; i < nbrs; ++i) {
    seed.push_back(NeighbourSimilarity{gfp.tanimoto(_pool[i]), i});
  }

  auto compare_by_similarity = [](const NeighbourSimilarity& lhs,
                                  const NeighbourSimilarity& rhs) {
    return lhs.similarity > rhs.similarity;
  };
  std::priority_queue<NeighbourSimilarity,
                      std::vector<NeighbourSimilarity>,
                      decltype(compare_by_similarity)>
      pq(compare_by_similarity, std::move(seed));

  for (int i = nbrs; i < static_cast<int>(_pool.size()); ++i) {
    const float similarity = gfp.tanimoto(_pool[i]);
    const NeighbourSimilarity& worst_kept = pq.top();
    if (similarity < worst_kept.similarity) {
      continue;
    }

    pq.pop();
    pq.push(NeighbourSimilarity{similarity, i});
  }

  std::vector<NeighbourDistance> did(nbrs);
  for (int ndx = nbrs - 1; ! pq.empty(); --ndx) {
    const NeighbourSimilarity item = pq.top();
    did[ndx].distance = 1.0f - item.similarity;
    did[ndx].index = item.index;
    pq.pop();
  }

  FillResultProto(did, result);
  return 1;
}

NnReply
GfpNearNeighbourServer::Search(const NnRequest& request) {
  NnReply reply;

  if (_pool.empty()) {
    reply.set_status(Status::SERVER_ERROR);
    reply.set_message("fingerprint pool is empty");
    return reply;
  }

  if (request.smiles().empty()) {
    reply.set_status(Status::NO_SMILES);
    reply.set_message("empty smiles");
    return reply;
  }

  GFP_Standard gfp;
  if (! SmilesToGfp(request.smiles(), gfp)) {
    reply.set_status(Status::BAD_SMILES);
    reply.set_message("invalid smiles");
    return reply;
  }

  if (! request.has_nbrs() && ! request.has_distance()) {
    reply.set_status(Status::NO_DIRECTIVE);
    reply.set_message("request must specify nbrs or distance");
    return reply;
  }

  nnbr::NearNeighbours* result = reply.mutable_result();
  if (! request.id().empty()) {
    result->set_name(request.id());
  }
  result->set_smiles(request.smiles());

  if (request.has_distance()) {
    if (request.distance() < 0.0f) {
      reply.set_status(Status::NO_DIRECTIVE);
      reply.set_message("distance must be non-negative");
      return reply;
    }
    TanimotoWithinDistance(gfp, request.distance(), *result);
    reply.set_status(Status::OK);
    return reply;
  }

  if (request.nbrs() == 0) {
    reply.set_status(Status::NO_DIRECTIVE);
    reply.set_message("nbrs must be positive");
    return reply;
  }

  if (request.nbrs() == 1) {
    TanimotoSingleNbr(gfp, *result);
  } else {
    Tanimoto(gfp, request.nbrs(), *result);
  }

  reply.set_status(Status::OK);
  return reply;
}


std::string
GfpNearNeighbourServer::SearchSerialized(std::string_view serialized_request) {
  NnRequest request;
  NnReply reply;

  if (! request.ParseFromArray(serialized_request.data(), serialized_request.size())) {
    reply.set_status(Status::SERVER_ERROR);
    reply.set_message("cannot parse NnRequest proto");
  } else {
    reply = Search(request);
  }

  std::string result;
  if (! reply.SerializeToString(&result)) {
    return std::string();
  }
  return result;
}

std::string
GfpNearNeighbourServer::SearchBatchSerialized(std::string_view serialized_request) {
  NnBatchRequest request;
  NnBatchReply reply;

  if (! request.ParseFromArray(serialized_request.data(), serialized_request.size())) {
    NnReply* item = reply.add_reply();
    item->set_status(Status::SERVER_ERROR);
    item->set_message("cannot parse NnBatchRequest proto");
  } else {
    reply = Search(request);
  }

  std::string result;
  if (! reply.SerializeToString(&result)) {
    return std::string();
  }
  return result;
}

NnBatchReply
GfpNearNeighbourServer::Search(const NnBatchRequest& request) {
  NnBatchReply result;
  result.mutable_reply()->Reserve(request.request_size());

  for (const NnRequest& req : request.request()) {
    *result.add_reply() = Search(req);
  }

  return result;
}

}  // namespace gfp_server
