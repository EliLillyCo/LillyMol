#ifndef MOLECULE_TOOLS_BDB_RING_DATABASE_LOOKUP_H_
#define MOLECULE_TOOLS_BDB_RING_DATABASE_LOOKUP_H_

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "Foundational/iwstring/iwstring.h"

class Db;

namespace smi2rings_bdb {

enum class LookupStatus { kFound, kNotFound, kError };

struct RingLookupResult {
  uint32_t _count = 0;
  std::string _exemplar;
  uint32_t _databases_matched = 0;
};

// Lookup existing unique-SMILES keys. Counts are the maximum across databases,
// not their sum. For equal maxima, the last database's exemplar is retained.
// On a miss or error the result is cleared. No ring extraction or output occurs.
LookupStatus LookupRingInDatabases(std::span<Db* const> databases, const IWString& key,
                                   RingLookupResult& result);

// Own read-only BerkeleyDB handles for repeated C++ lookups. Each instance has
// independent databases; there is no shared mutable lookup state. The existing
// on-disk key and Smi2Rings::Ring textproto formats are unchanged.
class RingDatabaseReader {
 private:
  class Impl;
  std::unique_ptr<Impl> _impl;

 public:
  RingDatabaseReader();
  ~RingDatabaseReader();
  RingDatabaseReader(const RingDatabaseReader&) = delete;
  RingDatabaseReader& operator=(const RingDatabaseReader&) = delete;

  // Replace databases only after every file has opened successfully. Empty
  // filename lists are errors. Files are opened read-only and never created.
  bool Open(const std::vector<std::string>& filenames);
  bool is_open() const;
  LookupStatus Lookup(const IWString& key, RingLookupResult& result) const;
};

}  // namespace smi2rings_bdb
#endif  // MOLECULE_TOOLS_BDB_RING_DATABASE_LOOKUP_H_
