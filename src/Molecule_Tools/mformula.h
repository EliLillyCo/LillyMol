#ifndef MOLECULE_TOOLS_MFORMULA_H_
#define MOLECULE_TOOLS_MFORMULA_H_

#include <cstdint>

// We have several situations where we want to restrict to minor changes
// to a molecle. One way of doing that is to restrict changes to the 
// molecular formula. This object holds a molecular formula for a molecule
// and can report a difference between two formulae.

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/set_of_atoms.h"

namespace mformula {

// This must be kept in sync with the constants defined in mformula.cc
inline constexpr int kMFOther = 17;

class MFormula {
  private:
    int _count[kMFOther + 1];

    // Some tools may need lazy evaluation.
    int _initialised;

  // Private functions
    void ZeroCountArray();
    int Build(Molecule& m, atom_number_t i);

  public:
    MFormula();

    int Build(Molecule& m);

    // Build only for the atoms in `embedding`.
    int Build(Molecule& m, const Set_of_Atoms& embedding);

    int initialised() const {
      return _initialised;
    }

    // The absolute difference between individual types.
    uint32_t Diff(const MFormula& rhs) const;

    int ToSparseFingerprint(IWString& destination) const;
    int ToFixedCountedFingerprint(IWString& destination) const;

};

}  // namespace mformula

#endif // MOLECULE_TOOLS_MFORMULA_H_
