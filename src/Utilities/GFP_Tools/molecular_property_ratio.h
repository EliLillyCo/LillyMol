#ifndef UTILITIES_GFP_TOOLS_MOLECULAR_PROPERTY_RATIO_H_
#define UTILITIES_GFP_TOOLS_MOLECULAR_PROPERTY_RATIO_H_

#include <cstdint>

#include "Foundational/iwbits/iwbits.h"

namespace gfp_internal {

class MolecularPropertyRatioTable {
 private:
  similarity_type_t _ratio[256 * 256];

 public:
  MolecularPropertyRatioTable();

  const similarity_type_t* data() const {
    return _ratio;
  }

  similarity_type_t ratio(uint8_t lhs, uint8_t rhs) const {
    return _ratio[static_cast<int>(lhs) * 256 + static_cast<int>(rhs)];
  }

  similarity_type_t unnormalised_similarity(const uint8_t* lhs, const uint8_t* rhs,
                                            int nproperties) const;

  similarity_type_t similarity(const uint8_t* lhs, const uint8_t* rhs,
                               int nproperties) const;
};

const MolecularPropertyRatioTable& molecular_property_ratio_table();

}  // namespace gfp_internal

#endif  // UTILITIES_GFP_TOOLS_MOLECULAR_PROPERTY_RATIO_H_
