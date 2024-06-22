#include "molecule.h"
#include "path.h"
#include "set_of_atoms.h"

#include "ring_data.h"

namespace lillymol {
RingInformation::RingInformation() {
  _nrings = 0;
  _ring = nullptr;
  _aromatic = nullptr;
  _fused = nullptr;
}

RingInformation::~RingInformation() {
  delete[] _ring;
  delete[] _aromatic;
  delete[] _fused;
}

int
RingInformation::GatherRings(Molecule& m, RIProperty properties) {
  if (_ring) {
    delete [] _ring;
    _ring = nullptr;
  }
  if (_aromatic) {
    delete [] _aromatic;
    _aromatic = nullptr;
  }
  if (_fused) {
    delete [] _fused;
    _fused = nullptr;
  }

  _nrings = m.nrings();

  if (_nrings == 0) {
    return 1;
  }

  _ring = new Set_of_Atoms[_nrings];
  for (int i = 0; i < _nrings; ++i) {
    _ring[i] = *m.ringi(i);
  }

  if (properties & RIProperty::kAromatic) {
    _aromatic = new int[_nrings];
    for (int i = 0; i < _nrings; ++i) {
      if (m.ringi(i)->is_aromatic()) {
        _aromatic[i] = 1;
      } else {
        _aromatic[i] = 0;
      }
    }
  }

  if (properties & RIProperty::kFused) {
    _fused = new int[_nrings];
    for (int i = 0; i < _nrings; ++i) {
      _fused[i] = m.ringi(i)->fused_ring_neighbours();
    } 
  }

  return _nrings;
}

}  // namespace lillymol
