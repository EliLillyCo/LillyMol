#include "ring_database_builder.h"

#include <iostream>

#include "google/protobuf/text_format.h"

#include "Foundational/iwstring/absl_hash.h"

#include "db_cxx.h"

#include "Molecule_Tools_Bdb/smi2rings.pb.h"

namespace smi2rings_bdb {
class RingDatabaseBuilder::Impl {
 public:
  std::unique_ptr<Db> _owned_database;
  Db* _database = nullptr;
  std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>> _cache;
  uint64_t _new_ring_systems = 0;
  bool _finished = false;

  ~Impl() {
    if (_owned_database) {
      _owned_database->close(0);
    }
  }

  bool
  Put(const IWString& key, const Smi2Rings::Ring& proto) {
    google::protobuf::TextFormat::Printer printer;
    printer.SetSingleLineMode(true);
    std::string text;
    if (!printer.PrintToString(proto, &text)) {
      return false;
    }
    Dbt k(const_cast<char*>(key.rawchars()), key.length());
    Dbt value(text.data(), text.size());
    if (const int rc = _database->put(nullptr, &k, &value, 0); rc != 0) {
      _database->err(rc, "RingDatabaseBuilder: cannot store key");
      return false;
    }
    return true;
  }
};

RingDatabaseBuilder::RingDatabaseBuilder() : _impl(std::make_unique<Impl>()) {
}

RingDatabaseBuilder::RingDatabaseBuilder(Db& database,
                                         const RingDatabaseBuildOptions& options)
    : RingDatabaseBuilder() {
  _impl->_database = &database;
  if (options._cache) {
    _impl->_cache = std::make_unique<absl::flat_hash_map<IWString, Smi2Rings::Ring>>();
  }
}

RingDatabaseBuilder::~RingDatabaseBuilder() = default;

bool
RingDatabaseBuilder::Open(const char* filename, const RingDatabaseBuildOptions& options) {
  if (_impl->_database) {
    std::cerr << "RingDatabaseBuilder: database already open\n";
    return false;
  }
  auto db = std::make_unique<Db>(nullptr, DB_CXX_NO_EXCEPTIONS);
  if (const int rc = db->open(nullptr, filename, nullptr, DB_BTREE, DB_CREATE, 0644);
      rc != 0) {
    db->err(rc, "RingDatabaseBuilder: cannot open '%s'", filename);
    return false;
  }
  _impl->_database = db.get();
  _impl->_owned_database = std::move(db);
  if (options._cache) {
    _impl->_cache = std::make_unique<absl::flat_hash_map<IWString, Smi2Rings::Ring>>();
  }
  return true;
}

bool
RingDatabaseBuilder::Store(const IWString& key, const IWString& exemplar) {
  if (!_impl->_database || _impl->_finished) {
    std::cerr << "RingDatabaseBuilder: no writable database or build already finished\n";
    return false;
  }
  if (_impl->_cache) {
    auto [iter, inserted] = _impl->_cache->try_emplace(key);
    if (inserted) {
      iter->second.set_ex(exemplar.AsString());
      ++_impl->_new_ring_systems;
    }
    iter->second.set_n(iter->second.n() + 1);
    return true;
  }
  Dbt k(const_cast<char*>(key.rawchars()), key.length());
  Dbt value;
  Smi2Rings::Ring proto;
  if (const int rc = _impl->_database->get(nullptr, &k, &value, 0); rc == DB_NOTFOUND) {
    proto.set_ex(exemplar.AsString());
    ++_impl->_new_ring_systems;
  } else if (rc != 0) {
    _impl->_database->err(rc, "RingDatabaseBuilder: cannot read key");
    return false;
  } else if (!google::protobuf::TextFormat::ParseFromString(
                 std::string(static_cast<const char*>(value.get_data()),
                             value.get_size()),
                 &proto)) {
    std::cerr << "RingDatabaseBuilder: invalid stored ring for '" << key << "'\n";
    return false;
  }
  proto.set_n(proto.n() + 1);
  return _impl->Put(key, proto);
}

bool
RingDatabaseBuilder::Finish() {
  if (!_impl->_database) {
    return false;
  }
  if (_impl->_finished) {
    return true;
  }
  if (_impl->_cache) {
    // Preserve the CLI's existing cached-load semantics: cached counts replace
    // existing values for matching keys, rather than adding to those counts.
    for (const auto& [key, proto] : *_impl->_cache) {
      if (!_impl->Put(key, proto)) {
        return false;
      }
    }
  }
  if (const int rc = _impl->_database->sync(0); rc != 0) {
    _impl->_database->err(rc, "RingDatabaseBuilder: cannot sync");
    return false;
  }
  _impl->_finished = true;
  _impl->_cache.reset();
  return true;
}

uint64_t
RingDatabaseBuilder::new_ring_systems() const {
  return _impl->_new_ring_systems;
}

size_t
RingDatabaseBuilder::pending_ring_systems() const {
  return _impl->_cache ? _impl->_cache->size() : 0;
}
}  // namespace smi2rings_bdb
