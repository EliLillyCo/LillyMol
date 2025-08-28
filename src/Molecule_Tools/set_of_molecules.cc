#include "Molecule_Tools/set_of_molecules.h"

namespace set_of_molecules {

using std::cerr;

int
SetOfMolecules::Build(data_source_and_type<Molecule>& input) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    *this << m;
  }

  return this->number_elements();
}

int
SetOfMolecules::Build(IWString& fname) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "SetOfMolecules::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

}  // namespace set_of_molecules
