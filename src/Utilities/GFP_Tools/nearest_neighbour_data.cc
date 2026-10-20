#include "nearest_neighbour_data.h"

#include <cmath>
#include <iostream>
#include <utility>

#include "Foundational/data_source/tfdatarecord.h"

namespace nearest_neighbour_data {

bool
NearestNeighbourData::SetNeighboursToRetain(int count) {
  if (_loaded || count < 1) {
    std::cerr << "NearestNeighbourData: retention count must be positive and set before "
                 "reading data\n";
    return false;
  }
  _neighbours_to_retain = count;
  return true;
}

bool
NearestNeighbourData::Read(const char* fname) {
  iw_tf_data_record::TFDataReader input(fname);
  if (!input.good()) {
    std::cerr << "NearestNeighbourData: cannot open '" << fname << "'\n";
    return false;
  }
  std::vector<nnbr::NearNeighbours> neighbours;
  std::unordered_map<std::string, size_t> id_to_index;
  while (true) {
    const auto record = input.Next();
    if (!record) {
      if (!input.good() || !input.eof()) {
        std::cerr << "NearestNeighbourData: error reading '" << fname << "'\n";
        return false;
      }
      break;
    }
    nnbr::NearNeighbours source;
    if (!source.ParseFromArray(record->data(), record->length())) {
      std::cerr << "NearestNeighbourData: invalid proto in '" << fname << "'\n";
      return false;
    }
    if (!source.has_name() || source.name().empty()) {
      std::cerr << "NearestNeighbourData: missing query identifier\n";
      return false;
    }
    if (!id_to_index.emplace(source.name(), neighbours.size()).second) {
      std::cerr << "NearestNeighbourData: duplicate query identifier '" << source.name()
                << "'\n";
      return false;
    }

    // Truncating a repeated proto field can retain the discarded allocations.
    // Build a fresh proto so stored memory reflects only retained neighbours.
    nnbr::NearNeighbours retained;
    retained.set_name(source.name());
    if (source.has_smiles()) {
      retained.set_smiles(source.smiles());
    }
    if (source.has_cluster()) {
      retained.set_cluster(source.cluster());
    }
    if (source.has_csize()) {
      retained.set_csize(source.csize());
    }
    for (const nnbr::Nbr& nbr : source.nbr()) {
      if (nbr.id() == source.name()) {
        continue;
      }
      if (!std::isfinite(nbr.dist()) || nbr.dist() < 0.0f) {
        std::cerr << "NearestNeighbourData: invalid distance for '" << source.name()
                  << "'\n";
        return false;
      }
      *retained.add_nbr() = nbr;
      if (retained.nbr_size() == _neighbours_to_retain) {
        break;
      }
    }
    neighbours.push_back(std::move(retained));
  }
  _neighbours.swap(neighbours);
  _id_to_index.swap(id_to_index);
  _loaded = true;
  return true;
}

const nnbr::NearNeighbours*
NearestNeighbourData::Lookup(const std::string& id) const {
  const auto iter = _id_to_index.find(id);
  if (iter == _id_to_index.end()) {
    return nullptr;
  }
  return &_neighbours[iter->second];
}

std::optional<float>
NearestNeighbourData::NearestDistance(const std::string& id) const {
  const nnbr::NearNeighbours* proto = Lookup(id);
  if (!proto || proto->nbr().empty()) {
    return std::nullopt;
  }
  return proto->nbr(0).dist();
}
}  // namespace nearest_neighbour_data
