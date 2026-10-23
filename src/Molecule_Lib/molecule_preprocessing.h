#ifndef MOLECULE_LIB_MOLECULE_PREPROCESSING_H_
#define MOLECULE_LIB_MOLECULE_PREPROCESSING_H_

class Molecule;
class Command_Line;

#include "Molecule_Lib/standardise.h"

namespace molecule_processing {

class MoleculePreprocessing {
 public:
  MoleculePreprocessing() = default;

  bool Initialise(Command_Line& cl);

  // non-const because _chemical_standardisation is not thread safe.
  // Returns > 0 if any changes are made.
  int Process(Molecule& m);

  bool active() const {
    return _active;
  }

  void set_reduce_to_largest_fragment(bool s) {
    _reduce_to_largest_fragment = s;
    if (s) {
      _active = true;
    }
  }

  void set_remove_chirality(bool s) {
    _remove_chirality = s;
    if (s) {
      _active = true;
    }
  }

  void set_remove_cis_trans_bonds(bool s) {
    _remove_cis_trans_bonds = s;
    if (s) {
      _active = true;
    }
  }

  void set_remove_isotopes(bool s) {
    _remove_isotopes = s;
    if (s) {
      _active = true;
    }
  }

  // Queries for current states.
  bool reduce_to_largest_fragment() const {
    return _reduce_to_largest_fragment;
  }
  bool remove_isotopes() const {
    return _remove_isotopes;
  }
  bool remove_chirality() const {
    return _remove_chirality;
  }

 private:
  bool _active = false;

  bool _reduce_to_largest_fragment = false;
  bool _remove_chirality = false;
  bool _remove_cis_trans_bonds = false;
  bool _remove_isotopes = false;

  Chemical_Standardisation _chemical_standardisation;
};

}  // namespace molecule_processing

#endif  // MOLECULE_LIB_MOLECULE_PREPROCESSING_H_
