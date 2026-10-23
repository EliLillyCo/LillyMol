#include "Utilities/GFP_Tools/molecular_property_ratio.h"

namespace gfp_internal {

MolecularPropertyRatioTable::MolecularPropertyRatioTable() {
  for (int i = 0; i < 256; ++i) {
    const similarity_type_t si = static_cast<similarity_type_t>(i);
    for (int j = 0; j < 256; ++j) {
      if (i == j) {
        _ratio[i * 256 + j] = 1.0f;
      } else if (i == 0 || j == 0) {
        _ratio[i * 256 + j] = 0.5f;
      } else {
        const similarity_type_t sj = static_cast<similarity_type_t>(j);
        if (i < j) {
          _ratio[i * 256 + j] = si / sj;
        } else {
          _ratio[i * 256 + j] = sj / si;
        }
      }
    }
  }
}

similarity_type_t
MolecularPropertyRatioTable::unnormalised_similarity(const uint8_t* lhs,
                                                     const uint8_t* rhs,
                                                     int nproperties) const {
  similarity_type_t result = 0.0f;
  for (int i = 0; i < nproperties; ++i) {
    result += ratio(lhs[i], rhs[i]);
  }

  return result;
}

similarity_type_t
MolecularPropertyRatioTable::similarity(const uint8_t* lhs, const uint8_t* rhs,
                                        int nproperties) const {
  if (nproperties == 0) {
    return 0.0f;
  }

  return unnormalised_similarity(lhs, rhs, nproperties) /
         static_cast<similarity_type_t>(nproperties);
}

const MolecularPropertyRatioTable&
molecular_property_ratio_table() {
  static const MolecularPropertyRatioTable result;
  return result;
}

}  // namespace gfp_internal
