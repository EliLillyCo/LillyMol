#ifndef MOLECULE_TOOLS_HIGHEST_RING_NUMBER_H
#define MOLECULE_TOOLS_HIGHEST_RING_NUMBER_H

#include <optional>

#include "Foundational/iwstring/iwstring.h"

namespace lillymol {

// Examine `smiles` and return the highest ring number used in the smiles.
// If there are no rings, 0 will be returned.
// If there are too many rings, the % sign is encountered, nullopt is returned.
// There is no smiles interpretation.
std::optional<int> HighestRingNumber(const IWString& smiles);

// Given a SAFE smiles, identify the unbalanced ring closures.
// by convention, the unbalanced rings that have been added are 
// all 2 digits and will be preceded by a % sign. This makes
// parsing easy.
int UnbalancedRingNumbers(const IWString& smiles, resizable_array<int>& ring_numbers);

// Given a smiles with isotopic labels, convert those to ring openings.
// [nnnC]C becomes [nnnC]%rC where `r` is derived from `ring_number`
// Note that `new_smiles` will be appended to.
int IsotopeToRingOpening(const IWString& smiles, int& ring_number, IWString& new_smiles);

// when assigning ring numbers in SAFE smiles the default behaviour is to just use
// a sequential ring number for each join. This class keeps track of which rings
// have been opened and closed and can re-use ring numbers.
class RingNumberControl {
  private:
    const int _max_rings;

    // For each ring number, is it allocated or free.
    int* _issued;
    int _next_ring_number;

  public:
    RingNumberControl(int lowest_ring_number, int max_rings);
    ~RingNumberControl();

    int GetRing();
    void OkToReuse(int ring_number);
};

}  // namespace lillymol

#endif
