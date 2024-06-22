#ifndef MOLECULE_LIB_RING_DATA_H_
#define MOLECULE_LIB_RING_DATA_H_

#include "molecule.h"
#include "set_of_atoms.h"

namespace lillymol {
// When making changes to a molecule, rings may get destroyed. So it can be
// convenient to store all the ring information before processing starts.

// Depending on what features are requested, we can store other information about
// the rings.
enum RIProperty {
  kNone = 0,
  kAromatic = 1,
  kFused = 2
};

class RingInformation {
  private:
    int _nrings;
    Set_of_Atoms*  _ring;

    // Optionally, we can also store variout other properties of the rings.
    // Will be stored if requested in the call to GatherRings.
    int* _aromatic;
    // Store fused_ring_neighbours() for each ring.
    int* _fused;

  public:
    RingInformation();
    ~RingInformation();

    int GatherRings(Molecule& m, RIProperty properties = RIProperty::kNone);

    int nrings() const {
      return _nrings;
    }

    const Set_of_Atoms& operator[](int ndx) const {
      assert(ndx >= 0 && ndx < _nrings);
      return _ring[ndx];
    }

    // Beware, no checking for null dereference.
    int aromatic(int ndx) const {
      assert(ndx >= 0 && ndx < _nrings);
      return _aromatic[ndx];
    }
    int fused(int ndx) const {
      assert(ndx >= 0 && ndx < _nrings);
      return _fused[ndx];
    }
};

}  // namespace lillymol

#endif // MOLECULE_LIB_RING_DATA_H_
