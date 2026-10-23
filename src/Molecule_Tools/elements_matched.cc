#include <algorithm>
#include <iostream>

#include "elements_matched.h"

namespace elements_matched {

using std::cerr;

ElementsMatched::ElementsMatched() {
  std::fill_n(_atomic_number, HIGHEST_ATOMIC_NUMBER + 1, 0);
  _non_periodic_table = 0;
}

void
ElementsMatched::Reset() {
  std::fill_n(_atomic_number, HIGHEST_ATOMIC_NUMBER + 1, 0);
  _non_periodic_table = 0;
}

void
ElementsMatched::Extra(const Element* e) {
  return Extra(e->atomic_number());
}

void
ElementsMatched::Extra(atomic_number_t z) {
  if (z < 0) {
    ++_non_periodic_table;
  } else if (z <= (HIGHEST_ATOMIC_NUMBER + 1)) {
    ++_atomic_number[z];
  } else {
    ++_non_periodic_table;
  }
}

int
ElementsMatched::Report(std::ostream& output) const {
  for (int i = 0; i <= HIGHEST_ATOMIC_NUMBER; ++i) {
    if (_atomic_number[i] == 0) {
      continue;
    }

    const Element* e = get_element_from_atomic_number(i);
    output << e->symbol() << ' ' << _atomic_number[i] << '\n';
  }

  if (_non_periodic_table > 0) {
    output << "NPT " << _non_periodic_table << '\n';
  }

  return output.good();
}

}  // namespace elements_matched
