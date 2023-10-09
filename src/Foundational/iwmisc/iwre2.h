#ifndef FOUNDATIONAL_IWMISC_IW_RE2_H
#define FOUNDATIONAL_IWMISC_IW_RE2_H

#include <memory>

#include "re2/re2.h"

#include "Foundational/iwstring/iwstring.h"

namespace iwre2 {
bool RE2FullMatch(const IWString& s, RE2& rx);
bool RE2FullMatch(const const_IWSubstring& s, RE2& rx);

bool RE2PartialMatch(const IWString& s, RE2& rx);
bool RE2PartialMatch(const const_IWSubstring& s, RE2& rx);
// Reset a unique_ptr.
bool RE2Reset(std::unique_ptr<RE2>& rx, const IWString& s);
bool RE2Reset(std::unique_ptr<RE2>& rx, const const_IWSubstring& s);
}
#endif  // FOUNDATIONAL_IWMISC_IW_RE2_H
