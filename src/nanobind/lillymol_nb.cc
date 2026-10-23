#include <nanobind/nanobind.h>

#include "nanobind/lillymol_nb_internal.h"

namespace nb = nanobind;

NB_MODULE(lillymol, m) {
  lillymol_nb::BindIo(m);
  lillymol_nb::BindAtomBond(m);
  lillymol_nb::BindChirality(m);
  lillymol_nb::BindSetOfAtomsAndRing(m);
  lillymol_nb::BindSubstructure(m);
  lillymol_nb::BindTSubstructure(m);
  lillymol_nb::BindMolecule(m);
  lillymol_nb::BindDescriptors(m);
  lillymol_nb::BindStandardise(m);
  lillymol_nb::BindFingerprint(m);
  lillymol_nb::BindReaction(m);
  lillymol_nb::BindTools(m);
}
