
#include <iostream>
#include <string>

namespace compile_time {
// A fixed date string is placed in compile_time.cc and returned.
const std::string& CompileDate();

const std::string& GitHash();

}  // namespace compile_time
