#ifndef MOLECULE_TOOLS_BDB_RING_PREPROCESSING_H_
#define MOLECULE_TOOLS_BDB_RING_PREPROCESSING_H_

#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

namespace smi2rings_bdb {
struct RingPreprocessingOptions {
  bool _largest_fragment = false;
  bool _remove_isotopes = false;
  bool _remove_chirality = false;
  // Preserve smi2rings_bdb's existing directional-bond handling.
  bool _strip_directional_bonds = true;

  void
  set_largest_fragment(bool value) {
    _largest_fragment = value;
  }

  void
  set_remove_isotopes(bool value) {
    _remove_isotopes = value;
  }

  void
  set_remove_chirality(bool value) {
    _remove_chirality = value;
  }

  void
  set_strip_directional_bonds(bool value) {
    _strip_directional_bonds = value;
  }
};

// Configured once and reused. Does not change process-wide configuration.
// Standardisation is initially inactive; activate a documented LillyMol
// directive ("all" is the profile for a future standardize=True binding).
class RingPreprocessor {
 private:
  RingPreprocessingOptions _options;
  Chemical_Standardisation _standardisation;
  Element_Transformations _element_transformations;

 public:
  RingPreprocessingOptions&
  options() {
    return _options;
  }

  const RingPreprocessingOptions&
  options() const {
    return _options;
  }

  Chemical_Standardisation&
  standardisation() {
    return _standardisation;
  }

  Element_Transformations&
  element_transformations() {
    return _element_transformations;
  }

  int SetStandardisation(const IWString& directive);
  void Process(Molecule& m);
};
}  // namespace smi2rings_bdb
#endif  // MOLECULE_TOOLS_BDB_RING_PREPROCESSING_H_
