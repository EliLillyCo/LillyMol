#ifndef UTILITIES_GENERAL_PROTO_COLLATOR_H_
#define UTILITIES_GENERAL_PROTO_COLLATOR_H_

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"

namespace proto_collator {

// A typed hash collator for protobuf-like records.
//
// The Traits type supplies all proto-specific behaviour. It must define:
//
//   using Proto = ...;
//   using Value = ...;
//
//   bool IsValid(const Proto& proto, std::ostream& err) const;
//   Key(const Proto& proto) const;  // Returns something convertible to std::string_view.
//   Value MakeValue(Proto&& proto, const std::string& key) const;
//   void Merge(Value& destination, Proto&& incoming) const;
//   int Write(const std::string& key, const Value& value,
//             IWString_and_File_Descriptor& output) const;
//
// If WriteSorted is used, Traits must also define:
//
//   bool SortPrecedes(const Value& lhs, const Value& rhs) const;
//
// This deliberately avoids protobuf reflection. The compiler emits direct field
// access through Traits, while this class owns the common hash/merge/write loop.
template <typename Traits>
class ProtoCollator {
 private:
  using Proto = typename Traits::Proto;
  using Value = typename Traits::Value;

  Traits _traits;
  absl::flat_hash_map<std::string, Value> _items;

  uint64_t _items_read = 0;
  uint64_t _duplicates = 0;
  uint64_t _errors = 0;

 public:
  ProtoCollator() = default;
  explicit ProtoCollator(Traits traits) : _traits(std::move(traits)) {}

  const Traits& traits() const {
    return _traits;
  }
  Traits& mutable_traits() {
    return _traits;
  }

  uint64_t items_read() const {
    return _items_read;
  }
  uint64_t duplicates() const {
    return _duplicates;
  }
  uint64_t errors() const {
    return _errors;
  }
  size_t size() const {
    return _items.size();
  }
  bool empty() const {
    return _items.empty();
  }
  void reserve(size_t n) {
    _items.reserve(n);
  }
  void clear() {
    _items.clear();
    _items_read = 0;
    _duplicates = 0;
    _errors = 0;
  }

  const absl::flat_hash_map<std::string, Value>& items() const {
    return _items;
  }
  absl::flat_hash_map<std::string, Value>& mutable_items() {
    return _items;
  }

  int Add(Proto&& proto, std::ostream& err) {
    ++_items_read;

    if (! _traits.IsValid(proto, err)) {
      ++_errors;
      return 0;
    }

    const auto& key_ref = _traits.Key(proto);
    if (key_ref.empty()) {
      err << "ProtoCollator::Add:empty key\n";
      ++_errors;
      return 0;
    }

    auto iter = _items.find(key_ref);
    if (iter != _items.end()) {
      _traits.Merge(iter->second, std::move(proto));
      ++_duplicates;
      return 1;
    }

    // The map must own its key. Copy it before moving proto, since key_ref may
    // refer to storage inside proto.
    std::string key(key_ref);
    _items.emplace(key, _traits.MakeValue(std::move(proto), key));
    return 1;
  }

  int Add(const Proto& proto, std::ostream& err) {
    Proto copy(proto);
    return Add(std::move(copy), err);
  }

  int AccumulateTextProtoRecord(const const_IWSubstring& buffer, std::ostream& err) {
    google::protobuf::io::ArrayInputStream zero_copy_array(buffer.data(), buffer.nchars());

    Proto proto;
    if (! google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
      err << "ProtoCollator::AccumulateTextProtoRecord:cannot parse '" << buffer << "'\n";
      ++_errors;
      return 0;
    }

    return Add(std::move(proto), err);
  }

  int AccumulateTextProtoRecords(iwstring_data_source& input, std::ostream& err) {
    const_IWSubstring buffer;
    while (input.next_record(buffer)) {
      if (! AccumulateTextProtoRecord(buffer, err)) {
        err << "ProtoCollator::AccumulateTextProtoRecords:cannot process line " <<
               input.lines_read() << '\n';
        return 0;
      }
    }

    return 1;
  }

  int Write(IWString_and_File_Descriptor& output) const {
    for (const auto& [key, value] : _items) {
      if (! _traits.Write(key, value, output)) {
        return 0;
      }
    }

    output.flush();
    return 1;
  }

  int WriteSorted(IWString_and_File_Descriptor& output) const {
    std::vector<std::pair<const std::string*, const Value*>> values;
    values.reserve(_items.size());
    for (const auto& [key, value] : _items) {
      values.emplace_back(&key, &value);
    }

    std::sort(values.begin(), values.end(),
              [this](const auto& lhs, const auto& rhs) {
                return _traits.SortPrecedes(*lhs.second, *rhs.second);
              });

    for (const auto& [key, value] : values) {
      if (! _traits.Write(*key, *value, output)) {
        return 0;
      }
    }

    output.flush();
    return 1;
  }
};

}  // namespace proto_collator

#endif  // UTILITIES_GENERAL_PROTO_COLLATOR_H_
