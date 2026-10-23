#include <iostream>

#include "Molecule_Lib/molecule_preprocessing.h"

#include "Foundational/cmdline/cmdline.h"
#include "Molecule_Lib/molecule.h"

namespace molecule_processing {

using std::cerr;

bool
MoleculePreprocessing::Initialise(Command_Line& cl) {
  const int verbose = cl.option_present('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = true;
    _active = true;
    if (verbose) {
      cerr << "Will reduce to the largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = true;
    _active = true;
    if (verbose) {
      cerr << "Chirality will be removed\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      return false;
    }

    _active = true;
  }

  return true;
}

int
MoleculePreprocessing::Process(Molecule& m) {
  int rc = 0;
  if (_reduce_to_largest_fragment && m.number_fragments() > 1) {
    rc += m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality && m.chiral_centres()) {
    rc += m.remove_all_chiral_centres();
  }

  if (_remove_cis_trans_bonds) {
    rc += m.revert_all_directional_bonds_to_non_directional();
  }

  if (_remove_isotopes) {
    rc += m.transform_to_non_isotopic_form();
  }

  if (_chemical_standardisation.active()) {
    rc += _chemical_standardisation.process(m);
  }

  return rc;
}

}  // namespace molecule_processing
