
#ifndef MOLECULE_LIB_CREATE_COMPONENTS_H_
#define MOLECULE_LIB_CREATE_COMPONENTS_H_

#include <memory>

#include "Foundational/iwaray/iwaray.h"

#include "molecule.h"

template <typename T>
int
Molecule::create_components(resizable_array_p<T>& components) {
  assert(ok());

  if (0 == _number_elements) {
    return 0;
  }

  if (_number_elements <= 1) {
    std::cerr << "Molecule::create_components: molecule one or fewer atoms\n";
    return 0;
  }

  const int nf = number_fragments();
  if (1 == nf) {
    std::cerr << "Molecule::create_components: molecule contains only one component\n";
    return 0;
  }

  // We need two arrays, one for fragment membership and the other for
  // the cross reference - that is not used here.
  std::unique_ptr<int[]> tmp =
      std::make_unique<int[]>(_number_elements + _number_elements);

  int* xref_not_used = tmp.get() + _number_elements;

  fragment_membership(tmp.get());

  assign_bond_numbers_to_bonds_if_needed();

  components.resize(nf);

  for (int i = 0; i < nf; i++) {
    //  std::cerr << "Creating component from fragment " << i << '\n';

    T* m = new T;
    create_subset(*m, tmp.get(), i, xref_not_used);
    components.add(m);
  }

  assert(nf == components.number_elements());

  return nf;
}

#endif // MOLECULE_LIB_CREATE_COMPONENTS_H_
