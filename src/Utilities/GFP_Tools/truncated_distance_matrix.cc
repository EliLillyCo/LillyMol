#include "Utilities/GFP_Tools/truncated_distance_matrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <utility>

#include "Foundational/data_source/tfdatarecord.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#else
#include "nearneighbours.pb.h"
#endif

namespace truncated_distance_matrix {

using std::cerr;

namespace {

using Row = std::vector<std::pair<uint32_t, uint8_t>>;

uint32_t
LowIndex(uint32_t i, uint32_t j) {
  return i < j ? i : j;
}

uint32_t
HighIndex(uint32_t i, uint32_t j) {
  return i < j ? j : i;
}

int
CheckIndex(uint32_t ndx, uint32_t n, const char* caller) {
  if (ndx < n) {
    return 1;
  }

  cerr << caller << ":index " << ndx << " outside valid range " << n << '\n';
  return 0;
}

template <typename Proto>
int
GatherIdentifiers(const char* fname, std::vector<std::string>& names,
                  absl::flat_hash_map<std::string, uint32_t>& id_to_ndx) {
  iw_tf_data_record::TFDataReader reader;
  if (! reader.Open(fname)) {
    cerr << "TruncatedDistanceMatrix::GatherIdentifiers:cannot open '" << fname << "'\n";
    return 0;
  }

  uint32_t ndx = 0;
  for (; ; ++ndx) {
    std::optional<Proto> proto = reader.ReadProto<Proto>();
    if (! proto) {
      break;
    }

    if (! proto->has_name() || proto->name().empty()) {
      cerr << "TruncatedDistanceMatrix::GatherIdentifiers:record " << ndx
           << " missing name\n";
      return 0;
    }

    const auto [_, inserted] = id_to_ndx.emplace(proto->name(), ndx);
    if (! inserted) {
      cerr << "TruncatedDistanceMatrix::GatherIdentifiers:duplicate name '"
           << proto->name() << "'\n";
      return 0;
    }
    names.emplace_back(proto->name());
  }

  return names.size();
}

int
AddDistance(uint32_t i, uint32_t j, uint8_t dist, std::vector<Row>& rows) {
  if (i == j) {
    if (dist != 0) {
      cerr << "TruncatedDistanceMatrix::AddDistance:non-zero self distance " << i
           << " dist " << static_cast<int>(dist) << '\n';
      return 0;
    }
    return 1;
  }

  const uint32_t row = LowIndex(i, j);
  const uint32_t col = HighIndex(i, j);
  rows[row].emplace_back(col, dist);

  return 1;
}

template <typename Proto>
int
ReadNameBasedNeighbours(const char* fname,
                        const absl::flat_hash_map<std::string, uint32_t>& id_to_ndx,
                        std::vector<Row>& rows, uint8_t& max_stored_distance) {
  iw_tf_data_record::TFDataReader reader;
  if (! reader.Open(fname)) {
    cerr << "TruncatedDistanceMatrix::ReadNameBasedNeighbours:cannot open '" << fname
         << "'\n";
    return 0;
  }

  uint32_t ndx = 0;
  for (; ; ++ndx) {
    std::optional<Proto> proto = reader.ReadProto<Proto>();
    if (! proto) {
      break;
    }

    if (! proto->has_name()) {
      cerr << "TruncatedDistanceMatrix::ReadNameBasedNeighbours:record " << ndx
           << " missing name\n";
      return 0;
    }
    const auto self = id_to_ndx.find(proto->name());
    if (self == id_to_ndx.end()) {
      cerr << "TruncatedDistanceMatrix::ReadNameBasedNeighbours:no index for '"
           << proto->name() << "'\n";
      return 0;
    }

    for (const nnbr::Nbr& nbr : proto->nbr()) {
      const auto iter = id_to_ndx.find(nbr.id());
      if (iter == id_to_ndx.end()) {
        cerr << "TruncatedDistanceMatrix::ReadNameBasedNeighbours:no index for neighbour '"
             << nbr.id() << "'\n";
        return 0;
      }

      const uint8_t d = DistanceToByte(nbr.dist());
      max_stored_distance = std::max(max_stored_distance, d);
      if (! AddDistance(self->second, iter->second, d, rows)) {
        return 0;
      }
    }
  }

  return 1;
}

template <typename Proto>
int
ReadIndexedNeighbours(const char* fname, uint32_t n, std::vector<Row>& rows,
                      uint8_t& max_stored_distance) {
  iw_tf_data_record::TFDataReader reader;
  if (! reader.Open(fname)) {
    cerr << "TruncatedDistanceMatrix::ReadIndexedNeighbours:cannot open '" << fname
         << "'\n";
    return 0;
  }

  uint32_t ndx = 0;
  for (; ; ++ndx) {
    std::optional<Proto> proto = reader.ReadProto<Proto>();
    if (! proto) {
      break;
    }
    if (! CheckIndex(ndx, n, "TruncatedDistanceMatrix::ReadIndexedNeighbours")) {
      return 0;
    }

    for (const nnbr::NbrNdx& nbr : proto->nbr()) {
      if (! CheckIndex(nbr.id(), n, "TruncatedDistanceMatrix::ReadIndexedNeighbours")) {
        return 0;
      }

      const uint8_t d = DistanceToByte(nbr.dist());
      max_stored_distance = std::max(max_stored_distance, d);
      if (! AddDistance(ndx, nbr.id(), d, rows)) {
        return 0;
      }
    }
  }

  return ndx == n;
}

int
SortAndCheckRows(std::vector<Row>& rows, uint64_t& number_distances,
                 uint64_t& duplicate_distances_differing_by_one) {
  number_distances = 0;
  for (uint32_t i = 0; i < rows.size(); ++i) {
    Row& row = rows[i];
    if (row.empty()) {
      continue;
    }

    std::sort(row.begin(), row.end(),
              [](const auto& p1, const auto& p2) { return p1.first < p2.first; });

    uint32_t write_ndx = 0;
    for (uint32_t read_ndx = 0; read_ndx < row.size(); ++read_ndx) {
      if (write_ndx > 0 && row[write_ndx - 1].first == row[read_ndx].first) {
        const uint8_t previous = row[write_ndx - 1].second;
        const uint8_t current = row[read_ndx].second;
        if (previous == current) [[likely]] {
          continue;
        }
        const int delta = std::abs(static_cast<int>(previous) - static_cast<int>(current));
        if (delta == 1) {
          row[write_ndx - 1].second = std::min(previous, current);
          ++duplicate_distances_differing_by_one;
          continue;
        }

        cerr << "TruncatedDistanceMatrix::SortAndCheckRows:conflicting distances for "
             << i << ',' << row[read_ndx].first << " "
             << static_cast<int>(previous) << " vs "
             << static_cast<int>(current) << '\n';
        return 0;
      }

      if (write_ndx != read_ndx) {
        row[write_ndx] = row[read_ndx];
      }
      ++write_ndx;
    }
    row.resize(write_ndx);
    number_distances += row.size();
  }

  return 1;
}

}  // namespace

TruncatedDistanceMatrix::TruncatedDistanceMatrix() {
  Clear();
}

void
TruncatedDistanceMatrix::Clear() {
  _storage = Storage::kRowSparse;
  _names.clear();
  _id_to_ndx.clear();
  _row_start.clear();
  _col.clear();
  _dist.clear();
  _row_hash.clear();
  _max_stored_distance = 0;
  _default_distance = 0;
  _max_stored_distance_float = 0.0f;
  _default_distance_float = 0.0f;
  _number_distances = 0;
  _duplicate_distances_differing_by_one = 0;
}

uint8_t
DistanceToByte(float d) {
  if (d <= 0.0f) {
    return 0;
  }
  if (d >= 1.0f) {
    return std::numeric_limits<uint8_t>::max();
  }

  return static_cast<uint8_t>(std::lround(d * 255.0f));
}

float
ByteToDistance(uint8_t d) {
  return static_cast<float>(d) / 255.0f;
}

int
TruncatedDistanceMatrix::Build(const char* fname, Storage storage, ProtoType proto_type) {
  Clear();
  _storage = storage;

  int n = 0;
  if (proto_type == ProtoType::kNearNeighbours) {
    n = GatherIdentifiers<nnbr::NearNeighbours>(fname, _names, _id_to_ndx);
  } else {
    n = GatherIdentifiers<nnbr::NearNeighboursIndices>(fname, _names, _id_to_ndx);
  }

  if (n <= 0) {
    cerr << "TruncatedDistanceMatrix::Build:no identifiers read from '" << fname << "'\n";
    Clear();
    return 0;
  }

  std::vector<Row> rows(n);

  int ok;
  if (proto_type == ProtoType::kNearNeighbours) {
    ok = ReadNameBasedNeighbours<nnbr::NearNeighbours>(fname, _id_to_ndx, rows,
                                                       _max_stored_distance);
  } else {
    ok = ReadIndexedNeighbours<nnbr::NearNeighboursIndices>(fname, n, rows,
                                                            _max_stored_distance);
  }
  if (! ok) {
    Clear();
    return 0;
  }

  if (! SortAndCheckRows(rows, _number_distances,
                         _duplicate_distances_differing_by_one)) {
    Clear();
    return 0;
  }

  if (_duplicate_distances_differing_by_one > 0) {
    cerr << "TruncatedDistanceMatrix::Build:" << _duplicate_distances_differing_by_one
         << " duplicate symmetric distances differed by 1 byte; stored smaller value\n";
  }

  _default_distance = _max_stored_distance;
  _max_stored_distance_float = ByteToDistance(_max_stored_distance);
  _default_distance_float = _max_stored_distance_float;

  if (storage == Storage::kRowSparse) {
    _row_start.resize(static_cast<size_t>(n) + 1);
    _col.reserve(_number_distances);
    _dist.reserve(_number_distances);

    uint64_t offset = 0;
    for (uint32_t i = 0; i < static_cast<uint32_t>(n); ++i) {
      _row_start[i] = offset;
      for (const auto& [col, dist] : rows[i]) {
        _col.push_back(col);
        _dist.push_back(dist);
        ++offset;
      }
    }
    _row_start[n] = offset;
  } else {
    _row_hash.resize(n);
    for (uint32_t i = 0; i < static_cast<uint32_t>(n); ++i) {
      _row_hash[i].reserve(rows[i].size());
      for (const auto& [col, dist] : rows[i]) {
        _row_hash[i].emplace(col, dist);
      }
    }
  }

  return 1;
}

std::optional<uint32_t>
TruncatedDistanceMatrix::Index(std::string_view name) const {
  const auto iter = _id_to_ndx.find(name);
  if (iter == _id_to_ndx.end()) {
    return std::nullopt;
  }

  return iter->second;
}

const std::string&
TruncatedDistanceMatrix::Name(uint32_t ndx) const {
  return _names[ndx];
}

std::optional<uint8_t>
TruncatedDistanceMatrix::StoredDistance(uint32_t i, uint32_t j) const {
  if (i >= _names.size() || j >= _names.size()) [[unlikely]] {
    return std::nullopt;
  }
  if (i == j) [[unlikely]] {
    return 0;
  }

  const uint32_t row = LowIndex(i, j);
  const uint32_t col = HighIndex(i, j);

  if (_storage == Storage::kRowHash) {
    const auto iter = _row_hash[row].find(col);
    if (iter == _row_hash[row].end()) [[likely]] {
      return std::nullopt;
    }
    return iter->second;
  }

  const uint64_t istart = _row_start[row];
  const uint64_t istop = _row_start[row + 1];
  const auto b = _col.begin() + istart;
  const auto e = _col.begin() + istop;
  const auto iter = std::lower_bound(b, e, col);
  if (iter == e || *iter != col) [[likely]] {
    return std::nullopt;
  }

  return _dist[iter - _col.begin()];
}

std::optional<float>
TruncatedDistanceMatrix::Distance(uint32_t i, uint32_t j) const {
  std::optional<uint8_t> d = StoredDistance(i, j);
  if (! d) [[likely]] {
    return std::nullopt;
  }

  return ByteToDistance(*d);
}

float
TruncatedDistanceMatrix::DistanceOrDefault(uint32_t i, uint32_t j) const {
  std::optional<uint8_t> d = StoredDistance(i, j);
  if (! d) [[likely]] {
    return _default_distance_float;
  }

  return ByteToDistance(*d);
}

std::vector<float>
TruncatedDistanceMatrix::DistancesOrDefault(const std::vector<uint32_t>& i,
                                            const std::vector<uint32_t>& j) const {
  std::vector<float> result;
  if (i.size() != j.size()) [[unlikely]] {
    return result;
  }

  result.reserve(i.size());
  for (uint32_t ndx = 0; ndx < i.size(); ++ndx) {
    result.push_back(DistanceOrDefault(i[ndx], j[ndx]));
  }

  return result;
}

int
TruncatedDistanceMatrix::SetDefaultDistance(float d) {
  return SetDefaultDistanceByte(DistanceToByte(d));
}

int
TruncatedDistanceMatrix::SetDefaultDistanceByte(uint8_t d) {
  if (d < _max_stored_distance) [[unlikely]] {
    cerr << "TruncatedDistanceMatrix::SetDefaultDistanceByte:cannot set default "
         << static_cast<int>(d) << " below max stored distance "
         << static_cast<int>(_max_stored_distance) << '\n';
    return 0;
  }

  _default_distance = d;
  _default_distance_float = ByteToDistance(d);
  return 1;
}

}  // namespace truncated_distance_matrix
