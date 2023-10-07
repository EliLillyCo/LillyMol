#include <filesystem>
#include <optional>
#include <string>

#include "Foundational/iwstring/iwstring.h"

namespace iwmisc {

namespace fs = std::filesystem;

std::optional<IWString>
FileOrPath(const IWString& outer_file_name,
           const IWString& inner_file_name) {
  const std::string possible_file(inner_file_name.data(), inner_file_name.length());

  if (fs::exists(fs::path(possible_file))) {
    // std::cerr << "File name " << possible_file << " alread exists, good\n";
    return inner_file_name;
  }

  // Try query_file_name in same directory as proto_file_name.

  // cerr << "outer_file_name '" << outer_file_name << "'\n";
  fs::path new_name = fs::path(std::string(outer_file_name.data(), outer_file_name.length())).replace_filename(possible_file);
  // cerr << "With directory '" << new_name.string() << "'\n";

  if (fs::exists(new_name)) {
    return new_name.string();
  }

  return std::nullopt;
}

}  // namespace iwmisc
