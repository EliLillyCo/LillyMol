#include <iostream>

#include "sorted_list.h"

namespace sorted_list {

using std::cerr;

uint32_t search_threshold = 6;

void
InsertIntoList(resizable_array<uint64_t>& sorted, uint64_t extra) {
  assert(extra > sorted.front());
  assert(extra < sorted.back());

  // If the list is "short" do a linear search.
  if (sorted.size() < search_threshold) {
    const int n = sorted.number_elements();
    for (int i = 0; i < n; ++i) {
      if (extra <= sorted[i]) {
        // std::cerr << "Find " << extra << " at " << i << " cmp " << sorted[i] << " size " << sorted.size() << '\n';
        sorted.insert_before(i, extra);
        return;
      }
    }
    cerr << "Could not insert into sorted list size " << sorted.size() << '\n';
    return;
  }

  // binary search to find insertion point.
  // We rely on the fact that there will never be an equality comparison.
  int left = 0;
  int right = sorted.number_elements() - 1;
  int middle = (left + right) / 2;

  const uint64_t* s = sorted.rawdata();

  // cerr << "INserting " << extra << " into list of size " << sorted.size() << '\n';

  while (middle > left) {
    const uint64_t m = s[middle];

    // cerr << "  cmp " << extra << " vs " << m << '\n';
    if (extra < m) {
      right = middle;
      if (left + 1 == middle) {
        sorted.insert_before(right, extra);
        break;
      }
    } else if (extra > m) {
      left = middle;
      if (middle + 1 == right) {
        sorted.insert_after(middle, extra);
        break;
      }
    }

    middle = (left + right) / 2;

    if (middle == left) {
      break;
    }

    // cerr << "Update " << left << " " << middle << " " << right << '\n';
  }

#define CHECK_SORTED_
#ifdef CHECK_SORTED_
  for (int i = 1; i < sorted.number_elements(); ++i) {
    if (sorted[i] < sorted[i-1]) {
      cerr << "Out of order i " << i << " prev " << sorted[i-1] << " vs " << sorted[i] << '\n';
    }
  }
#endif

}

}  // namespace sorted_list
