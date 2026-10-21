#include "ring_database_lookup.h"

#include <iostream>
#include <utility>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwaray/iwaray.h"

#include "db_cxx.h"

#include "Molecule_Tools_Bdb/smi2rings.pb.h"

namespace smi2rings_bdb {

LookupStatus
LookupRingInDatabases(std::span<Db* const> databases, const IWString& key,
                      RingLookupResult& result) {
  result = RingLookupResult{};
  if (databases.empty()) {
    std::cerr << "RingDatabaseReader: no databases open\n";
    return LookupStatus::kError;
  }
  RingLookupResult found;
  for (Db* database : databases) {
    Dbt dkey(const_cast<char*>(key.rawchars()), key.length());
    Dbt data;
    const int rc = database->get(nullptr, &dkey, &data, 0);
    if (rc == DB_NOTFOUND) {
      continue;
    }
    if (rc != 0) {
      database->err(rc, "RingDatabaseReader: lookup failed");
      return LookupStatus::kError;
    }
    google::protobuf::io::ArrayInputStream input(data.get_data(), data.get_size());
    Smi2Rings::Ring ring;
    if (!google::protobuf::TextFormat::Parse(&input, &ring)) {
      std::cerr << "RingDatabaseReader: invalid ring record for '" << key << "'\n";
      return LookupStatus::kError;
    }
    ++found._databases_matched;
    if (ring.n() >= found._count) {
      found._count = ring.n();
      found._exemplar = ring.ex();
    }
  }
  if (found._databases_matched == 0) {
    return LookupStatus::kNotFound;
  }
  result = std::move(found);
  return LookupStatus::kFound;
}

class RingDatabaseReader::Impl {
 public:
  resizable_array_p<Db> _databases;

  ~Impl() {
    for (Db* database : _databases) {
      database->close(0);
    }
  }
};

RingDatabaseReader::RingDatabaseReader() : _impl(std::make_unique<Impl>()) {
}

RingDatabaseReader::~RingDatabaseReader() = default;

bool
RingDatabaseReader::Open(const std::vector<std::string>& filenames) {
  if (filenames.empty()) {
    std::cerr << "RingDatabaseReader: no database filenames supplied\n";
    return false;
  }
  auto next = std::make_unique<Impl>();
  for (const std::string& filename : filenames) {
    auto database = std::make_unique<Db>(nullptr, DB_CXX_NO_EXCEPTIONS);
    if (const int rc =
            database->open(nullptr, filename.c_str(), nullptr, DB_UNKNOWN, DB_RDONLY, 0);
        rc != 0) {
      database->err(rc, "RingDatabaseReader: cannot open '%s'", filename.c_str());
      return false;
    }
    next->_databases.add(database.release());
  }
  _impl = std::move(next);
  return true;
}

LookupStatus
RingDatabaseReader::Lookup(const IWString& key, RingLookupResult& result) const {
  return LookupRingInDatabases(
      std::span<Db* const>(_impl->_databases.rawdata(), _impl->_databases.size()), key,
      result);
}
}  // namespace smi2rings_bdb
