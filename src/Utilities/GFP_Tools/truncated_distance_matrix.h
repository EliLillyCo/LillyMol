#ifndef UTILITIES_GFP_TOOLS_TRUNCATED_DISTANCE_MATRIX_H
#define UTILITIES_GFP_TOOLS_TRUNCATED_DISTANCE_MATRIX_H

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"

namespace truncated_distance_matrix {

enum class Storage {
  kRowSparse = 0,
  kRowHash = 1,
};

enum class ProtoType {
  kNearNeighbours = 0,
  kNearNeighboursIndices = 1,
};

class TruncatedDistanceMatrix {
 private:
  Storage _storage;

  std::vector<std::string> _names;
  absl::flat_hash_map<std::string, uint32_t> _id_to_ndx;

  std::vector<uint64_t> _row_start;
  std::vector<uint32_t> _col;
  std::vector<uint8_t> _dist;

  std::vector<absl::flat_hash_map<uint32_t, uint8_t>> _row_hash;

  uint8_t _max_stored_distance;
  uint8_t _default_distance;
  float _max_stored_distance_float;
  float _default_distance_float;
  uint64_t _number_distances;
  uint64_t _duplicate_distances_differing_by_one;

  void Clear();

  std::optional<uint8_t> StoredDistance(uint32_t i, uint32_t j) const;

 public:
  TruncatedDistanceMatrix();

  int Build(const char* fname, Storage storage = Storage::kRowSparse,
            ProtoType proto_type = ProtoType::kNearNeighbours);
  int Build(const std::string& fname, Storage storage = Storage::kRowSparse,
            ProtoType proto_type = ProtoType::kNearNeighbours) {
    return Build(fname.c_str(), storage, proto_type);
  }

  uint32_t size() const {
    return _names.size();
  }

  uint64_t number_distances() const {
    return _number_distances;
  }

  uint64_t duplicate_distances_differing_by_one() const {
    return _duplicate_distances_differing_by_one;
  }

  Storage storage() const {
    return _storage;
  }

  std::optional<uint32_t> Index(std::string_view name) const;
  const std::string& Name(uint32_t ndx) const;

  std::optional<float> Distance(uint32_t i, uint32_t j) const;
  float DistanceOrDefault(uint32_t i, uint32_t j) const;

  std::vector<float> DistancesOrDefault(const std::vector<uint32_t>& i,
                                        const std::vector<uint32_t>& j) const;

  float MaxStoredDistance() const {
    return _max_stored_distance_float;
  }
  float DefaultDistance() const {
    return _default_distance_float;
  }
  uint8_t MaxStoredDistanceByte() const {
    return _max_stored_distance;
  }
  uint8_t DefaultDistanceByte() const {
    return _default_distance;
  }

  int SetDefaultDistance(float d);
  int SetDefaultDistanceByte(uint8_t d);
};

uint8_t DistanceToByte(float d);
float ByteToDistance(uint8_t d);

}  // namespace truncated_distance_matrix

#endif  // UTILITIES_GFP_TOOLS_TRUNCATED_DISTANCE_MATRIX_H
