// Do not edit, overwritten by pre-commit githook.
#include <string>

#include "compile_time.h"

namespace compile_time {
const std::string&
CompileDate() {
  static const std::string compiled = "2022-Sep-07";
  return compiled;
}

const std::string&
GitHash() {
  static const std::string git_hash = "f74e3df";
  return git_hash;
}

}  // namespace compile_time
