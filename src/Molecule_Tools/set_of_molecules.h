#ifndef MOLECULE_TOOLS_SET_OF_MOLECULES_H_
#define MOLECULE_TOOLS_SET_OF_MOLECULES_H_

#include <iostream>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmtypes.h"

namespace set_of_molecules {

template <typename M>
class SetOfMolecules : public resizable_array_p<M> {
  private:

  // private functions.
    int Build(data_source_and_type<M>& input);

  public:
    int Build(IWString& fname);
};

template <typename M>
int
SetOfMolecules<M>::Build(data_source_and_type<M>& input) {
  M* m;
  while ((m = input.next_molecule()) != nullptr) {
    *this << m;
  }

  return this->number_elements();
}

template <typename M>
int
SetOfMolecules<M>::Build(IWString& fname) {
  data_source_and_type<M> input(FILE_TYPE_SMI, fname);
  if (! input.good()) {
    std::cerr << "SetOfMolecules::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

}  // namespace set_of_molecules

#endif  // MOLECULE_TOOLS_SET_OF_MOLECULES_H_
