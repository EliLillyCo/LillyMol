#include "ring_preprocessing.h"

namespace smi2rings_bdb {
int
RingPreprocessor::SetStandardisation(const IWString& directive) {
  return _standardisation.Activate(directive, 0);
}

void
RingPreprocessor::Process(Molecule& m) {
  if (_options._largest_fragment) {
    m.reduce_to_largest_fragment();
  }
  if (_element_transformations.active()) {
    _element_transformations.process(&m);
  }
  if (_standardisation.active()) {
    _standardisation.process(m);
  }
  // Remove input isotopes before atom typing; later attachment labels are kept.
  if (_options._remove_isotopes) {
    m.transform_to_non_isotopic_form();
  }
  if (_options._strip_directional_bonds) {
    m.revert_all_directional_bonds_to_non_directional();
  }
  if (_options._remove_chirality) {
    m.remove_all_chiral_centres();
  }
}
}  // namespace smi2rings_bdb
