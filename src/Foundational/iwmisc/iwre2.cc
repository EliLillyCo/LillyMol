#include <string_view>
#include "absl/strings/string_view.h"
#include "Foundational/iwmisc/iwre2.h"

namespace iwre2 {
bool RE2FullMatch(const IWString& s, RE2& rx) {
  absl::string_view tmp(s.data(), s.length());
  return RE2::FullMatch(tmp, rx);
}
bool RE2FullMatch(const const_IWSubstring& s, RE2& rx) {
  absl::string_view tmp(s.data(), s.length());
  return RE2::FullMatch(tmp, rx);
}
bool RE2PartialMatch(const IWString& s, RE2& rx) {
  absl::string_view tmp(s.data(), s.length());
  return RE2::PartialMatch(tmp, rx);
}
bool RE2PartialMatch(const const_IWSubstring& s, RE2& rx) {
  absl::string_view tmp(s.data(), s.length());
  return RE2::PartialMatch(tmp, rx);
}

bool RE2Reset(std::unique_ptr<RE2>& rx, const IWString& s) {
  absl::string_view tmp(s.data(), s.length());
  rx.reset(new RE2(tmp));
  return rx->ok();
}
bool RE2Reset(std::unique_ptr<RE2>& rx, const const_IWSubstring& s) {
  absl::string_view tmp(s.data(), s.length());
  rx.reset(new RE2(tmp));
  return rx->ok();
}

}  // namespace iwre2
