#ifndef MOLECULE_TOOLS_ELEMENTS_MATCHED_H
#define MOLECULE_TOOLS_ELEMENTS_MATCHED_H
// When filtering for non-organic elements, tools need
// a means of keeping track of just which elements have
// been discarded.

#include <cstdint>

#include "Molecule_Lib/element.h"

namespace elements_matched {

class ElementsMatched {
  private:
    // For each atomic number, the number of instances.
    uint64_t _atomic_number[HIGHEST_ATOMIC_NUMBER + 1];

    // For non periodic table elements.
    uint64_t _non_periodic_table;

  public:
    ElementsMatched();

    void Reset();

    void Extra(const Element* e);
    void Extra(const atomic_number_t z);

    const uint64_t* atomic_number() const {
      return _atomic_number;
    }

    uint64_t non_periodic_table() const {
      return _non_periodic_table;
    }

    int Report(std::ostream& output) const;
};

}  // namespace elements_matched
#endif // MOLECULE_TOOLS_ELEMENTS_MATCHED_H
