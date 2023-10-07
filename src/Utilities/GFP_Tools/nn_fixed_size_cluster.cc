// Implement fixed size cluter idea from Howard

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

namespace fixed_size_cluster {

using std::cerr;

IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");
IWString distance_tag("DIST<");

void
WriteSmiles(const IWString& smiles,
            IWString_and_File_Descriptor& output) {
  output << smiles_tag << smiles << ">\n";
}
void
WriteIdentifier(const IWString& id,
            IWString_and_File_Descriptor& output) {
  output << identifier_tag << id << ">\n";
}
void
WriteDistance(float d,
              IWString_and_File_Descriptor& output) {
  output << distance_tag << d << ">\n";
}

void
WriteSmilesIdDistance(const IWString& smiles,
                      const IWString& id,
                      float dist,
                      IWString_and_File_Descriptor& output) {
  WriteSmiles(smiles, output);
  WriteIdentifier(id, output);
  WriteDistance(dist, output);
  output.write_if_buffer_holds_more_than(4096);
}

int
DataItemValue(const_IWSubstring buffer,  // our own copy, we change it.
              const IWString& tag,
              IWString& result) {
  buffer.remove_leading_chars(tag.length());
  buffer.chop();
  result = buffer;
  if (result.empty()) {
    cerr << "No data\n";
    return 0;
  }

  return 1;
}

std::optional<float>
DataItemValue(const_IWSubstring buffer,
              const IWString& tag) {
  buffer.remove_leading_chars(tag.length());
  buffer.chop();
  float result;
  if (! buffer.numeric_value(result) || result < 0.0 || result > 1.0) {
    return std::nullopt;
  }

  return result;
}

class Neighbours {
  private:
    IWString _id;
    IWString _smiles;
    resizable_array<uint32_t> _nbr;
    resizable_array<float> _distance;

    float _ave_distance;

  public:
    Neighbours();

    int Build(iwstring_data_source& input,
              const IW_STL_Hash_Map<IWString, uint32_t>& id_to_ndx);

    float ave_distance() const {
      return _ave_distance;
    }

    const IWString& smiles() const {
      return _smiles;
    }

    const IWString& identifier() const {
      return _id;
    }

    const IWString* ptr_smiles() const {
      return &_smiles;
    }
    const IWString* ptr_identifier() const {
      return &_id;
    }

    float max_distance() const {
      return _distance.back();
    }

    // For each neighbour and distance, populate `distance_matrix`.
    int PopulateDistanceMatrix(const IW_STL_Hash_Map<IWString, uint32_t>& id_to_ndx, 
                               std::unordered_map<uint32_t, float>& distance_matrix) const;

    // A new cluster is being formed, with size `max_cluster_size`.
    // Add as any of our neighbours as we can to that cluster, updating
    // `already_assigned` as we go.
    int AddClusterMembers(uint32_t * already_assigned,
                int cluster_number,
                resizable_array<uint32_t>& cluster,
                uint32_t max_cluster_size,
                const IWString** smiles,
                const IWString** identifier,
                IWString_and_File_Descriptor& output) const;
};

Neighbours::Neighbours() {
  _ave_distance = -1.0f;
}

int
Neighbours::Build(iwstring_data_source& input,
                  const IW_STL_Hash_Map<IWString, uint32_t>& id_to_ndx) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Neighbours::Build:empty\n";
    return 0;
  }
  if (! DataItemValue(buffer, smiles_tag, _smiles)) {
    cerr << "Neighbours::Build:bad smiles '" << buffer << "'\n";
    return 0;
  }
  if (! input.next_record(buffer)) {
    cerr << "Neighbours::Build:empty\n";
    return 0;
  }
  if (! DataItemValue(buffer, identifier_tag, _id)) {
    cerr << "Neighbours::Build:bad identifier '" << buffer << "'\n";
    return 0;
  }

  IWString id;  // scope here for efficiency only.
  int current_id_ndx = -1;
  while (input.next_record(buffer)) {
    if (buffer.starts_with(smiles_tag)) {
      continue;
    } else if (buffer == '|') {
      break;
    } else if (buffer.starts_with(identifier_tag)) {
      DataItemValue(buffer, identifier_tag, id);
      const auto f = id_to_ndx.find(id);
      if (f == id_to_ndx.end()) {
        cerr << "Neighbours::Build:unrecognised identifier '" << id << "'\n";
        return 0;
      }
      current_id_ndx = f->second;
    } else if (buffer.starts_with(distance_tag)) {
      std::optional<float> maybe_distance = DataItemValue(buffer, distance_tag);
      if (! maybe_distance) {
        cerr << "Neighbours::Build:invalid distance '" << buffer << "'\n";
        return 0;
      }
      if (current_id_ndx < 0) {
        cerr << "Neighbours::Build:no PCN with distance\n";
        return 0;
      }
      _nbr << static_cast<uint32_t>(current_id_ndx);
      _distance << *maybe_distance;
      current_id_ndx = -1;
    }
  }

  float total = 0.0;
  for (float d : _distance) {
    total += d;
  }

  _ave_distance = total / static_cast<float>(_distance.size());

  return 1;
}


int
Neighbours::PopulateDistanceMatrix(const IW_STL_Hash_Map<IWString, uint32_t>& id_to_ndx, 
                               std::unordered_map<uint32_t, float>& distance_matrix) const {
  const auto f = id_to_ndx.find(_id);
  if (f == id_to_ndx.end()) {
    cerr << "Neighbours::PopulateDistanceMatrix:huh no hash match '" << _id << "'\n";
    return 0;
  }

  const uint32_t target_ndx = f->second;

  assert(_distance.size() == _nbr.size());

  const uint32_t n = id_to_ndx.size();

  const uint32_t number_neighbours = _nbr.size();
  for (uint32_t i = 0; i < number_neighbours; ++i) {
    uint32_t nbr_ndx = _nbr[i];
    uint32_t ndx = target_ndx * n + nbr_ndx;
    distance_matrix.try_emplace(ndx, _distance[i]);
    ndx = nbr_ndx * n + target_ndx;
    distance_matrix.try_emplace(ndx, _distance[i]);
  }

  return 1;
}

int
Neighbours::AddClusterMembers(uint32_t * already_assigned,
                int cluster_number,
                resizable_array<uint32_t>& cluster,
                uint32_t max_cluster_size,
                const IWString** smiles,
                const IWString** identifier,
                IWString_and_File_Descriptor& output) const {
  const uint32_t n = _nbr.size();

  int rc = 0;
  for (uint32_t i = 0; i < n && cluster.size() < max_cluster_size; ++i) {
    uint32_t j = _nbr[i];
    if (already_assigned[j]) {
      continue;
    }
    already_assigned[j] = cluster_number;
    cluster << j;
    WriteSmilesIdDistance(*smiles[j], *identifier[j], _distance[i], output);
    ++rc;
  }

  return rc;
}

class FixedSizeClusterImpl {
  private:
    int _verbose;

    uint32_t _cluster_size;

    Neighbours* _neighbours;

    // The size of the _neighbours array.
    uint32_t _nseeds;
  
    // A mapping from (i * _nseeds + j) to a distance.
    std::unordered_map<uint32_t, float> _distance_matrix;

    float _max_distance;

    // when finishing clusters, we need an array of prospective
    // ids and distances that we can sort.
    struct IDDistance {
      uint32_t id;
      float dist;
    };

    IDDistance* _sortable;

    extending_resizable_array<int> _sizes_formed;

  // private functions
    int CompleteCluster(uint32_t seed,
                uint32_t* already_assigned,
                int cluster_number,
                resizable_array<uint32_t>& cluster,
                const IWString** smiles,
                const IWString** identifier,
                IWString_and_File_Descriptor& output);

    float DistanceBetween(uint32_t i, uint32_t j) const;

  public:
    FixedSizeClusterImpl();
    ~FixedSizeClusterImpl();

    int Initialise(Command_Line& cl);

    int IngestData(const char* fname, IWString_and_File_Descriptor& output);
    int IngestData(iwstring_data_source& input, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

FixedSizeClusterImpl::FixedSizeClusterImpl() {
  _verbose = 0;
  _cluster_size = 0;
  _max_distance = 0.0f;
  _neighbours = nullptr;
  _sortable = nullptr;
}

FixedSizeClusterImpl::~FixedSizeClusterImpl() {
  if (_neighbours != nullptr) {
    delete [] _neighbours;
  }
  if (_sortable != nullptr) {
    delete [] _sortable;
  }
}

int
FixedSizeClusterImpl::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('s')) {
    cerr << "FixedSizeClusterImpl::Initialise:must specify cluster size via the -s option\n";
    return 0;
  }

  if (! cl.value('s', _cluster_size) || _cluster_size < 1) {
    cerr << "FixedSizeClusterImpl::Initialise:invalid cluster size (-s)\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Will generate clusters of size " << _cluster_size << '\n';
  }

  return 1;
}

int
FixedSizeClusterImpl::IngestData(const char * fname,
                        IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "FixedSizeClusterImpl::IngestData:cannot open '" << fname << "'\n";
    return 0;
  }

  return IngestData(input, output);
}

int
FixedSizeClusterImpl::IngestData(iwstring_data_source& input,
                        IWString_and_File_Descriptor& output) {
  // for each ID a mapping to a number.
  IW_STL_Hash_Map<IWString, uint32_t> id_to_ndx;
  // First pass through the data is to fill id_to_ndx
  IWString id;
  const_IWSubstring buffer;
  bool next_is_seed = true;
  while (input.next_record(buffer)) {
    if (buffer == '|') {
      next_is_seed = true;
      continue;
    }

    if (next_is_seed && buffer.starts_with(identifier_tag)) {
      DataItemValue(buffer, identifier_tag, id);
      auto f = id_to_ndx.find(id);
      if (f != id_to_ndx.end()) {
        cerr << "FixedSizeClusterImpl::IngestData:Duplicate target id '" << id << "'\n";
        return 0;
      }
      auto s = id_to_ndx.size();
      id_to_ndx.emplace(std::make_pair(id, s));
      next_is_seed = false;
    }
  }

  _nseeds = id_to_ndx.size();
  if (_verbose) {
    cerr << input.lines_read() << " lines contain " << _nseeds << " identifiers\n";
  }

  if (! input.seekg(0)) {
    cerr << "FixedSizeClusterImpl::IngestData:cannot seek back to start of file\n";
    return 0;
  }

  _neighbours = new Neighbours[_nseeds];
  _sortable = new IDDistance[_nseeds];

  for (uint32_t i = 0; i < _nseeds; ++i) {
    if (! _neighbours[i].Build(input, id_to_ndx)) {
      cerr << "FixedSizeClusterImpl::IngestData:cannot build item " << i << '\n';
      return 0;
    }
    //cerr << "i " << i << " dist " <<  _neighbours[i].ave_distance() << '\n';
    if (_neighbours[i].max_distance() > _max_distance) {
      _max_distance = _neighbours[i].max_distance();
    }
  }

  if (_verbose) {
    cerr << "Max distance " << _max_distance << '\n';
  }

  for (uint32_t i = 0; i < _nseeds; ++i) {
    _neighbours[i].PopulateDistanceMatrix(id_to_ndx, _distance_matrix);
  }

  std::unique_ptr<int[]> sorted = std::make_unique<int[]>(_nseeds);
  std::iota(sorted.get(), sorted.get() + _nseeds, 0);
  std::sort(sorted.get(), sorted.get() + _nseeds,
                [this](int i1, int i2) {
                  return _neighbours[i1].ave_distance() > _neighbours[i2].ave_distance();
                }
  );

  // We need an array of smiles and identifiers
  std::unique_ptr<const IWString*[]> smiles = std::make_unique<const IWString*[]>(_nseeds);
  std::unique_ptr<const IWString*[]> identifier(new const IWString*[_nseeds]);
  for (uint32_t i = 0; i < _nseeds; ++i) {
    smiles[i] = _neighbours[i].ptr_smiles();
    identifier[i] = _neighbours[i].ptr_identifier();
  }

  std::unique_ptr<uint32_t[]> already_assigned(new_unsigned_int(_nseeds));
  uint32_t clusters_to_form = _nseeds / _cluster_size;
  if (clusters_to_form * _cluster_size < _nseeds) {
    ++clusters_to_form;
  }
  if (_verbose) {
    cerr << "Clustering " << _nseeds << " into clusters of size " << _cluster_size
         << " will produce " << clusters_to_form << " clusters\n";
  }

  int cluster_number = 1;
  for (uint32_t i = 0; i < _nseeds; ++i, ++cluster_number) {
    if (already_assigned[i]) {
      continue;
    }

    already_assigned[i] = cluster_number;
    WriteSmiles(_neighbours[i].smiles(), output);
    WriteIdentifier(_neighbours[i].identifier(), output);

    resizable_array<uint32_t> cluster;
    cluster.resize(_cluster_size);
    cluster << i;
    _neighbours[i].AddClusterMembers(already_assigned.get(), cluster_number,
                cluster, _cluster_size,
                smiles.get(), identifier.get(),
                output);
    if (cluster.size() < _cluster_size) {
      CompleteCluster(i, already_assigned.get(), cluster_number,
                     cluster, smiles.get(), identifier.get(), output);
    }
    if (_verbose > 1) {
      cerr <<"Formed cluster of size " << cluster.size() << '\n';
    }
    _sizes_formed[cluster.size()]++;
    output << "|\n";
  }

  return 1;
}

float
FixedSizeClusterImpl::DistanceBetween(uint32_t i, uint32_t j) const {
  const auto f = _distance_matrix.find(i * _nseeds + j);
  if (f == _distance_matrix.end()) {
    return _max_distance;
  }

  return f->second;
}

// Cluster `cluster_number`, which is centered on `seed` has been partially
// formed and needs to be completed.
int
FixedSizeClusterImpl::CompleteCluster(uint32_t seed,
                uint32_t* already_assigned,
                int cluster_number,
                resizable_array<uint32_t>& cluster,
                const IWString** smiles,
                const IWString** identifier,
                IWString_and_File_Descriptor& output) {

  uint32_t candidates = 0;
  for (uint32_t i = 0; i < _nseeds; ++i) {
    if (already_assigned[i]) {
      continue;
    }
    _sortable[candidates].id = i;
    _sortable[candidates].dist = DistanceBetween(seed, i);
    ++candidates;
  }

  std::sort(_sortable, _sortable + candidates,
                [](const IDDistance& d1, IDDistance& d2) {
                  return d1.dist < d2.dist;
                }
  );

  // cerr << "Seed " << seed << " starts with " << cluster.size() << " members, have " << candidates << " candidates\n";

  uint32_t nprocess = _cluster_size - cluster.size();
  if (nprocess > candidates) {
    nprocess = candidates;
  }

  for (uint32_t i = 0; i < nprocess; ++i) {
    const uint32_t j = _sortable[i].id;
    already_assigned[j] = cluster_number;
    cluster << j;

    WriteSmilesIdDistance(*smiles[j], *identifier[j], _sortable[i].dist, output);
  }

  return 1;
}

int
FixedSizeClusterImpl::Report(std::ostream& output) const {
  output << "Clustered " << _nseeds << " seeds\n";
  for (int i = 0; i < _sizes_formed.number_elements(); ++i) {
    if (_sizes_formed[i] > 0) {
      output << _sizes_formed[i] << " clusters of size " << i << '\n';
    }
  }
  return output.good();
}

void
Usage(int rc) {
  cerr << "Generates fixed size clusters from a nearest neighbors file\n";
  cerr << " -s <size>     size of the fixed size clusters\n";
  cerr << " -v            verbose output\n";

  ::exit(rc);
}

int FixedSizeCluster(int argc, char** argv) {
  Command_Line cl(argc, argv, "vs:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered << '\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  FixedSizeClusterImpl fsc;
  if (! fsc.Initialise(cl)) {
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.size() > 1) {
    cerr << "Only takes one argument\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  if (! fsc.IngestData(cl[0], output)) {
    cerr << "Cannot initialise '" << cl[0] << "'\n";
    return 1;
  }

  if (verbose) {
    fsc.Report(cerr);
  }

  return 0;
}

}  // namespace fixed_size_cluster
int
main (int argc, char ** argv) {
  int rc = fixed_size_cluster::FixedSizeCluster(argc, argv);

  return rc;
}
