#ifndef MOLECULE_TOOLS_BDB_RING_SYSTEM_LOOKUP_H_
#define MOLECULE_TOOLS_BDB_RING_SYSTEM_LOOKUP_H_

#include <optional>
#include <vector>

#include "ring_database_lookup.h"
#include "ring_extraction.h"
#include "ring_preprocessing.h"

namespace smi2rings_bdb {
struct RingSystemMatch {
  int _system_id = -1;
  int _ring_count = 0;
  bool _found = false;
  std::optional<RingKeyType> _representation;
  IWString _key;
  RingLookupResult _database;
};

// Reusable molecule-level lookup, independent of CLI output and filtering.
// Configure once, open databases once, then process molecules repeatedly.
class RingSystemLookup {
 private:
  RingExtractor _extractor;
  RingPreprocessor _preprocessor;
  RingDatabaseReader _databases;

 public:
  RingSystemLookup();

  RingExtractor&
  extractor() {
    return _extractor;
  }

  RingPreprocessor&
  preprocessor() {
    return _preprocessor;
  }

  bool
  Open(const std::vector<std::string>& filenames) {
    return _databases.Open(filenames);
  }

  // Preserve the input with one working copy per call. No per-database or
  // per-representation full-molecule copies. Failure clears the result vector.
  bool Lookup(const Molecule& m, std::vector<RingSystemMatch>& results);
  // The caller owns a disposable working molecule; preprocessing and extraction
  // may modify its isotopes, directional bonds, hydrogen counts, and charges.
  bool LookupInPlace(Molecule& m, std::vector<RingSystemMatch>& results);
};
}  // namespace smi2rings_bdb
#endif  // MOLECULE_TOOLS_BDB_RING_SYSTEM_LOOKUP_H_
