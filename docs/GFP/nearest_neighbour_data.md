# NearestNeighbourData

`nearest_neighbour_data::NearestNeighbourData` is a reusable C++ lookup over
TFDataRecord files of serialized `nnbr::NearNeighbours` protos. It is independent
of fingerprint generation and prediction analysis.

```cpp
nearest_neighbour_data::NearestNeighbourData neighbours;
if (!neighbours.SetNeighboursToRetain(1) || !neighbours.Read("dataset.tfrecord")) {
  // Configuration or input error; methods write diagnostics and return false.
}
const nnbr::NearNeighbours* proto = neighbours.Lookup("mol1");
std::optional<float> distance = neighbours.NearestDistance("mol1");
```

Include `Utilities/GFP_Tools/nearest_neighbour_data.h`; depend on the public
Bazel target `//Utilities/GFP_Tools:nearest_neighbour_data`.

The default retention count is one. Set a positive count before reading. Calling
`SetNeighboursToRetain` after any successful `Read` is an error, even if the file
was empty or the requested count is unchanged. It returns false without changing
the configuration. A new instance is needed for a different retention count.

A successful `Read` replaces the collection and builds an ID-to-index map. Failed
reads leave the previous collection intact. Query IDs must be present, nonempty,
and unique. IDs are exact strings; no tokenisation or normalisation is done.
Neighbour ordering is trusted. Copy the first configured number of non-self
neighbours, preserving their IDs, SMILES, and distances. Reject negative and
non-finite retained distances; distances are not restricted to a maximum of one.

Each stored proto is built afresh, copying metadata and only retained neighbours.
This avoids keeping allocations for discarded repeated-field entries. Reading
still temporarily requires one complete source proto and its serialized record.
Replacing a loaded collection also temporarily retains the previous collection
until the new file has been successfully read.

`Lookup` returns a const pointer without copying the proto, or null for an absent
ID. An existing record can have an empty neighbour list. The pointer remains
valid until the next successful `Read` or object destruction. `NearestDistance`
returns the first retained distance, or `std::nullopt` for an absent ID or empty
list. Zero distances to distinct IDs are valid. Data must not be reloaded while
other callers use pointers into the collection.

The reader accepts `NearNeighbours`, not `NearNeighboursIndices`. It does not
infer which dataset or training split the neighbours belong to. See
[mispredicted](../Utilities/mispredicted.md) for its first consumer.
