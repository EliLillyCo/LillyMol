#ifndef MOLECULE_LIB_MOLECULAR_FORMULA_H_
#define MOLECULE_LIB_MOLECULAR_FORMULA_H_

#include <algorithm>
#include <iostream>

#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/molecule.h"

namespace molecular_formula {

// Note that these constants are different from any used by LillyMol.
constexpr uint32_t kCarbon = 0;
constexpr uint32_t kNitrogen = 1;
constexpr uint32_t kOxygen = 2;
constexpr uint32_t kFluorine = 3;
constexpr uint32_t kPhosphorus = 4;
constexpr uint32_t kSulphur = 5;
constexpr uint32_t kChlorine = 6;
constexpr uint32_t kBromine = 7;
constexpr uint32_t kIodine = 8;
constexpr uint32_t kBoron = 9;
constexpr uint32_t kOther = 10;

// The size of the _count array.
inline constexpr int kNTypes = kOther + 1;


// We have need for something to hold information about a molecular
// formula, especially in the case where we start with a smiles and
// do not necessarily instantiate a Molecule object.
// This is an efficiency play, so we make no allowance for all elements
// just the organic subset, with an indicator of presence of other types.
// But note that for atoms inside square brackets, any one or two letter
// combination is accepted, so [X] and [235Yz] are good. This does
// not know anything about the periodic table.
// Boron is not recognised, but would be simple to add.
// Note that since smiles interpretation is not performed, "c" is
// a perfectly valid input. Again, the assumption is that valid smiles
// will be passed to this object.
template <typename T>
class MolecularFormula {
  private:
    // For each kind of element, the number of instances.
    // T can be any integral type. Not sure what will work best.
    // TODO:ianwatson investigate different possibilities. If
    // we have something short, should we compare as uint64_t
    // instead of the hash?
    T _count[kNTypes];

    // The highway hash of the counts.
    uint64_t _hash;

    // Is the hash valid or not?
    int _hash_valid;

    // In some circumstances, computing the hash may not make sense -
    // the item will be used once.
    int _enable_hash_computation;

    IWDigits _digits;

  // Private functions.

    int IdentifyElementSQB(const const_IWSubstring& smiles, int& i);

  public:
    MolecularFormula();

    int Build(const const_IWSubstring& smiles);
    int Build(const Molecule& m);

    // This method is not const because if the hash has not been computed
    // it will be computed.
    int operator==(MolecularFormula& rhs);

    // Kind of dangerous function,
    const T* mf() const {
      return _count;
    }

    // Given _count, create a molecular formula in `s`. Existing
    // contents are overwritten.
    int MakeFormula(IWString& s) const;
    // A variant that appends instead.
    int AppendFormula(IWString& s) const;
};

};

#endif  // MOLECULE_LIB_MOLECULAR_FORMULA_H_
