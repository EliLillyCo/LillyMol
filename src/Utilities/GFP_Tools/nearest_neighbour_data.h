#ifndef UTILITIES_GFP_TOOLS_NEAREST_NEIGHBOUR_DATA_H_
#define UTILITIES_GFP_TOOLS_NEAREST_NEIGHBOUR_DATA_H_

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "Utilities/GFP_Tools/nearneighbours.pb.h"

namespace nearest_neighbour_data {

// An indexed, compact collection of sorted neighbour lists from TFDataRecord.
// IDs are exact strings. Input ordering is trusted; self neighbours are omitted.
class NearestNeighbourData {
 private:
  int _neighbours_to_retain = 1;
  bool _loaded = false;
  std::vector<nnbr::NearNeighbours> _neighbours;
  std::unordered_map<std::string, size_t> _id_to_index;

 public:
  // Positive counts only. Once Read succeeds, changing the count is an error,
  // including for an empty input file. Failure leaves the count unchanged.
  bool SetNeighboursToRetain(int count);

  int
  neighbours_to_retain() const {
    return _neighbours_to_retain;
  }

  size_t
  size() const {
    return _neighbours.size();
  }

  // Replace the collection on success. A failed read preserves existing data.
  // Reject missing/empty query IDs, duplicate query IDs and invalid retained
  // distances. Lists are assumed sorted; no sorting or order checks are done.
  bool Read(const char* fname);

  // Null if absent. An existing record may have no retained neighbours.
  // Pointers remain valid until the next successful Read or object destruction.
  const nnbr::NearNeighbours* Lookup(const std::string& id) const;
  std::optional<float> NearestDistance(const std::string& id) const;
};

}  // namespace nearest_neighbour_data
#endif  // UTILITIES_GFP_TOOLS_NEAREST_NEIGHBOUR_DATA_H_
