#ifndef FOUNDATIONAL_IWMISC_TOML_SUPPORT_H_
#define FOUNDATIONAL_IWMISC_TOML_SUPPORT_H_

#include <iostream>
#include <optional>
#include <string>
#include <sstream>

#define TOML_EXCEPTIONS 0 // only necessary if you've left them enabled in your compiler
#define TOML_ENABLE_FORMATTERS 1

#include "absl/status/status.h"

#include "toml.hpp"
#include "google/protobuf/util/json_util.h"

#include "Foundational/iwstring/iwstring.h"

namespace iwmisc {


template <typename Proto>
std::optional<Proto>
TomlToProto(toml::table& tbl) {
  std::stringstream ss;
  ss << toml::json_formatter{tbl};
  const std::string as_json = ss.str();

  Proto proto;

  google::protobuf::util::JsonParseOptions options;
  auto status = google::protobuf::util::JsonStringToMessage(as_json, &proto, options);

  // cannot figure out how to do error checking properly. I never really understood Status
  // while I was at Google, still don't.
  if (status != absl::OkStatus()) {
    std::cerr << "TomlToProto:cannot parse JSON\n";
    std::cerr << as_json << '\n';
    return std::nullopt;
  }

  return proto;
}

template <typename Proto>
std::optional<Proto>
ReadTomlProto(const IWString& fname) {
  const std::string string_fname(fname.data(), fname.length());

  auto as_toml = toml::parse_file(string_fname);
  if (! as_toml) {
    std::cerr << "ReadTomlProto:parse_file failed " << string_fname << '\n';
    return std::nullopt;
  }

  return TomlToProto<Proto>(std::move(as_toml));
}

template <typename Proto>
std::optional<Proto>
ParseFromToml(const std::string& toml_string) {
  auto as_toml = toml::parse(toml_string);
  if (! as_toml) {
    std::cerr << "ParseFromToml:Cannot parse as TOML\n";
    std::cerr << toml_string << '\n';
    return std::nullopt;
  }

  return TomlToProto<Proto>(std::move(as_toml));
}

}

#endif // FOUNDATIONAL_IWMISC_TOML_SUPPORT_H_
