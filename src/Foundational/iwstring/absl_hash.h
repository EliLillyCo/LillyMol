#ifndef FOUNDATIONAL_IWMISC_ABSL_HASH_H
#define FOUNDATIONAL_IWMISC_ABSL_HASH_H

#include <memory>

#include "absl/hash/hash.h"

#include "Foundational/iwstring/iwstring.h"


template <typename H>
H AbslHashValue(H state, const IWString& s) {
  if (s.empty()) {
    return H::combine(std::move(state), 709902);   // arbitrary number.
  }

  absl::string_view sv(s.data(), s.length());
  return H::combine(std::move(state), sv);
}

#endif  // FOUNDATIONAL_IWMISC_ABSL_HASH_H
