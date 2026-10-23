# Truncated Distance Matrix

`TruncatedDistanceMatrix` is a C++ and Python API for fast lookup of
precomputed pair-wise GFP distances when the full dense distance matrix is too
large to store. It is intended for large collections where only short distances
are useful and all omitted pairs can be treated as a constant large distance.

The object reads TFDataRecord files containing serialized near-neighbour protos,
builds an identifier table, and stores only the distances present in the file.
A query for a stored pair returns that stored distance. A query for an omitted
pair can either return `std::nullopt`/`None`, or return the object's default
distance. The default distance will default to the longest distance encountered
in the input file, but can be set to any larger value.

## Generating Input

Generate fingerprints, then write near-neighbour TFDataRecord output with a
distance threshold. The threshold form is important because it produces a
truncated symmetric distance matrix: all neighbours within the threshold are
written, and all omitted pairs are interpreted as far away.

```shell
gfp_make.sh file.smi > file.gfp
gfp_nearneighbours_single_file_tbb -T 0.4 -h 8 -S file.nn.tfdata file.gfp
```

The file must contain one record for every molecule whose identifier should be
known to the matrix, including molecules with no neighbours. Do not use `-z` for
files that will be loaded as a complete truncated distance matrix, because `-z`
omits no-neighbour records and the loader cannot infer those missing
identifiers.

By default, `-S` output from `gfp_nearneighbours_single_file_tbb` writes
`nnbr::NearNeighbours` records whose neighbours are identified by name. Add the
`indexed` keyword to write `nnbr::NearNeighboursIndices` records instead:

```shell
gfp_nearneighbours_single_file_tbb -T 0.4 -h 8 -S indexed -S file.nn.tfdata file.gfp
```

Indexed records are the preferred form for new truncated distance matrix files.
They are smaller and faster to load because neighbour identifiers are stored as
input row numbers rather than strings. The record `name` field is still written
for every molecule so the matrix can support identifier lookup.

## Distance Semantics

Distances are stored as bytes using

```text
byte = round(clamp(distance, 0.0, 1.0) * 255)
distance = byte / 255.0
```

The matrix is symmetric. Internally, each pair is stored once with canonical
ordering. If both directions are present in the input, duplicate distances must
agree. Duplicate distances that differ by one byte are accepted, counted, and
resolved by storing the smaller byte value. Larger disagreements are rejected as
malformed input.

The default distance is initialized to the largest stored distance observed in
the file. It can be raised by the caller, but cannot be lowered below the largest
stored distance. The diagonal distance is always zero.

## Storage Modes

Two storage modes are available.

`row_sparse` stores row offsets, sorted column indices, and byte distances. It is
the default and usually the best first choice for large data because it has low
and predictable memory overhead. Lookup performs a binary search within one
row.

`row_hash` stores one hash map per row. It uses more memory, but can be faster
for random lookup workloads where memory is available.

The public API is the same for both modes, so callers can switch storage mode
without changing lookup code.

## C++ API

The C++ class is `truncated_distance_matrix::TruncatedDistanceMatrix` in
`Utilities/GFP_Tools/truncated_distance_matrix.h`.

```cpp
#include "Utilities/GFP_Tools/truncated_distance_matrix.h"

using truncated_distance_matrix::ProtoType;
using truncated_distance_matrix::Storage;
using truncated_distance_matrix::TruncatedDistanceMatrix;

TruncatedDistanceMatrix dm;
if (! dm.Build("file.nn.tfdata", Storage::kRowSparse, ProtoType::kNearNeighbours)) {
  // handle malformed or unreadable input
}

std::optional<uint32_t> maybe_index = dm.Index("CHEMBL123");
std::optional<float> maybe_distance = dm.Distance(10, 17);
float distance = dm.DistanceOrDefault(10, 17);

dm.SetDefaultDistance(1.0f);
```

Important methods include:

| Method | Description |
| --- | --- |
| `Build(fname, storage, proto_type)` | Read a TFDataRecord near-neighbour file |
| `size()` | Number of molecule records known to the matrix |
| `number_distances()` | Number of stored pair distances after symmetric duplicate merging |
| `Index(name)` | Return the row index for an identifier, or `std::nullopt` |
| `Name(i)` | Return the identifier for row `i` |
| `Distance(i, j)` | Return a stored distance, zero for diagonal, or `std::nullopt` |
| `DistanceOrDefault(i, j)` | Return a stored distance, zero for diagonal, or the default distance |
| `DistancesOrDefault(i, j)` | Batch lookup for matching vectors of row indices |
| `MaxStoredDistance()` | Largest stored distance after byte discretisation |
| `DefaultDistance()` | Current default distance for omitted pairs |
| `SetDefaultDistance(d)` | Raise or set the default distance; fails if below max stored distance |
| `duplicate_distances_differing_by_one()` | Count of symmetric duplicate pairs resolved by one-byte tolerance |

## Python API

The Python binding is in `lillymol_tools`.

```python
from lillymol_tools import (
    TruncatedDistanceMatrix,
    TruncatedDistanceMatrixStorage,
    TruncatedDistanceMatrixProto,
)

dm = TruncatedDistanceMatrix(
    "file.nn.tfdata",
    storage=TruncatedDistanceMatrixStorage.ROW_SPARSE,
    proto_type=TruncatedDistanceMatrixProto.NEARNEIGHBOURS,
)

indexed_dm = TruncatedDistanceMatrix(
    "file.indexed.nn.tfdata",
    storage=TruncatedDistanceMatrixStorage.ROW_SPARSE,
    proto_type=TruncatedDistanceMatrixProto.NEARNEIGHBOURS_INDICES,
)

idx = dm.index("CHEMBL123")
print(dm.name(17))
print(dm.distance(10, 17))            # float or None
print(dm.distance_or_default(10, 17))
print(dm.distances_or_default([0, 1], [2, 3]))

dm.set_default_distance(1.0)
```

Invalid row indices raise `IndexError`. Unknown names return `None`. Build
failures, including duplicate identifiers or conflicting symmetric distances,
raise `RuntimeError` in Python.

For training loops, prefer `distances_or_default` over calling
`distance_or_default` once per pair from Python. The batch method still returns a
Python list today, but it avoids some per-call overhead and gives a natural place
for future NumPy support.
