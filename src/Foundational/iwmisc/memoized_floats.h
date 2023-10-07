// When writing out coordinates, a lot of time can be spent doing floating
// point conversions. That process can be sped up by storing precomputed
// string representations

#include "Foundational/iwmisc/iwdigits.h"

namespace memoized_floats {
class MemoizedFloats {
  public:
    MemoizedFloats();

    int Build(int max_int, int float_digits);

    IWString Representation(float value) const;

  private:
    IWDigits _positive_ints;
    IWDigits _negative_ints;
    IWDigits _fractions;
    int _fraction_max_value;
};

}  // namespace memoized_floats

