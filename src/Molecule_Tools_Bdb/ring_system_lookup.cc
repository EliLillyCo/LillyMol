#include "ring_system_lookup.h"

#include <iostream>
#include <utility>

namespace smi2rings_bdb {
RingSystemLookup::RingSystemLookup() {
  // Match the CLI's default when no representation options are specified.
  _extractor.options().set_label_attachment_points(1);
}

bool
RingSystemLookup::Lookup(const Molecule& m, std::vector<RingSystemMatch>& results) {
  Molecule working(m);
  return LookupInPlace(working, results);
}

bool
RingSystemLookup::LookupInPlace(Molecule& m, std::vector<RingSystemMatch>& results) {
  results.clear();
  if (!_databases.is_open()) {
    std::cerr << "RingSystemLookup: no databases open\n";
    return false;
  }
  if (_extractor.KeyTypes()._size == 0) {
    std::cerr << "RingSystemLookup: no key representations enabled\n";
    return false;
  }
  _preprocessor.Process(m);
  if (m.nrings() == 0) {
    return true;
  }
  RingExtractionData data(m, _extractor.options());
  if (!_extractor.Prepare(m, data)) {
    return false;
  }
  const int n = _extractor.IdentifySystems(m, data);
  for (int uid = 0; uid < n; ++uid) {
    if (!_extractor.PrepareSystem(m, data, uid)) {
      continue;
    }
    RingSystemMatch match;
    match._system_id = uid;
    match._ring_count = data._ring_system_size[uid];
    for (RingKeyType type : _extractor.KeyTypes()) {
      if (!_extractor.CreateKey(m, data, uid, type, match._key)) {
        results.clear();
        return false;
      }
      const LookupStatus status = _databases.Lookup(match._key, match._database);
      if (status == LookupStatus::kError) {
        results.clear();
        return false;
      }
      if (status == LookupStatus::kFound) {
        match._found = true;
        match._representation = type;
        break;
      }
      if (type != RingKeyType::kBare) {
        m.transform_to_non_isotopic_form();
      }
    }
    results.push_back(std::move(match));
  }
  return true;
}
}  // namespace smi2rings_bdb
