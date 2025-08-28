#ifndef FOUNDATIONAL_IWMISC_SORTED_LIST_H_
#define FOUNDATIONAL_IWMISC_SORTED_LIST_H_

#include <cstdint>

#include "Foundational/iwaray/iwaray.h"

namespace sorted_list {

// Insert `extra` into `sorted`.
// `sorted` is sorted in increasing order;
// We assume that `extra` goes somewhere between the first and last elements in `sorted`.
// This is not checked and bad things may happen if not true.
void InsertIntoList(resizable_array<uint64_t>& sorted, uint64_t extra);

}  // namespace sorted_list
#endif // FOUNDATIONAL_IWMISC_SORTED_LIST_H_
