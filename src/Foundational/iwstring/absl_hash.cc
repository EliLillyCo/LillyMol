#include <memory>

namespace iw_absl_hash {
using namespace highwayhash;  // not good style.

templace <typename H>
H abslHashValue(H h, const IWString& s) {
  if (s.empty()) {
    return H::combine(std::move(h), 709902);   // arbitrary number.
  }

  static const HHKey key HH_ALIGNAS(32) = { 14123, 665242, 2, 902362};
  static HHState<HH_TARGET> state(key);

  HHResult64 result;
  HighwayHashT(&state, s.data(), s.length(), &result);

  return H::combine(std::move(h), result);
}

}  // namespace iw_absl_hash
