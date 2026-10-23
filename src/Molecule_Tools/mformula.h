#ifndef MOLECULE_TOOLS_MFORMULA_H_
#define MOLECULE_TOOLS_MFORMULA_H_

#include <cstdint>

// We have several situations where we want to restrict to minor changes
// to a molecle. One way of doing that is to restrict changes to the 
// molecular formula. This object holds a molecular formula for a molecule
// and can report a difference between two formulae.

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/set_of_atoms.h"

class Sparse_Fingerprint_Creator;

namespace mformula {

// This must be kept in sync with the constants defined in mformula.cc
inline constexpr int kMFOther = 17;

class MFormula {
  private:
    int _count[kMFOther + 1];

    // We have the ability to separately consider aromatic forms.
    // by default, this will be true.
    // Note that ring perception is also turned off.
    int _consider_aromatic;

    // By default all elements are identical in weight. But in
    // reality, we would like to have more emphasis on rarer
    // elements. We want the difference between 20 and 21 carbon
    // atoms to be much less than the difference between zero sulphur
    // atoms and one sulphur atom.
    // So we introduce a log scale. THis is designed to keep all values
    // distinct up to 60 which is more then enough. From Claude
    // round(200 * ln(x + 1)/ln(41))
    //  call these parameters s (200) and d(41)
    // These result in transformations like
    //  1   37
    //  2   59
    //  3   75
    //  4   87
    //  5   96
    //  6  105
    //  7  112
    //  8  118
    //  9  124
    // 10  129
    // With the first collision at 60. Lower s and these numbers compress
    // and the first collision will occur sooner.

    float _log_s = 0.0f;
    float _log_d = 0.0f;

    // Some tools may need lazy evaluation.
    int _initialised;

    // Sum of _count.
    int _natoms;

  // Private functions
    void ZeroCountArray();
    int Build(Molecule& m, atom_number_t i);

    void ApplyLogCountScaling(Sparse_Fingerprint_Creator& sfc) const;

  public:
    MFormula();

    // This should be called before any calls to Build.
    void
    set_consider_aromatic(int s) {
      _consider_aromatic = s;
    }

    // See formula above.
    void
    set_log_scaling_factors(float s, float d) {
      _log_s = s;
      _log_d = d;
    }

    // Build from the text in `smiles`, no Molecule interpretation.
    // Note that something like [Np] would be counted as a Nitrogen atom.
    // Parsing is for speed, and it assumes well behaved input - organic only.
    // Nothing aromatic is perceived.
    int Build(const IWString& smiles);

    int Build(Molecule& m);

    // Build only for the atoms in `embedding`.
    int Build(Molecule& m, const Set_of_Atoms& embedding);

    int initialised() const {
      return _initialised;
    }

    // Only organic heavy atoms are counted. Note that this is only valid if
    // the formula has been built from a string. Change if ever needed.
    int natoms() const {
      return _natoms;
    }

    // The absolute difference between individual types.
    uint32_t Diff(const MFormula& rhs) const;

    // Returns true if `this` is a subset of `rhs`.
    bool IsSubset(const MFormula& rhs) const;
    // This version only checks the element counts, not heteroatoms and rings...
    bool IsElementCountSubset(const MFormula& rhs) const;

    int ToSparseFingerprint(IWString& destination) const;

    // Fill `destination` with the fixed counted fingerprint representation.
    // `n` must be at least kMFOther + 1. Counts are clamped to [0, 255].
    int ToFixedCountedFingerprint(int* destination, int n) const;

    // These are mostly used by tests.
    // Note that they do not consider aromatic variants.
    int Carbon() const;
    int Nitrogen() const;
    int Oxygen() const;
    int Fluorine() const;
    int Phosphorus() const;
    int Sulphur() const;
    int Chlorine() const;
    int Bromine() const;
    int Iodine() const;

};

}  // namespace mformula

#endif // MOLECULE_TOOLS_MFORMULA_H_
