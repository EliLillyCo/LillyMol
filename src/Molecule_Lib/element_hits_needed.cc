#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "misc2.h"
#include "substructure.h"
#include "target.h"

using std::cerr;

Elements_Needed::Elements_Needed()
{
  _operator = Op::NOT_SET;

  return;
}

int
Elements_Needed::ok() const
{
  return 1;
}

int
Elements_Needed::debug_print(std::ostream & os, const IWString & ind) const
{
  os << ind << "Elements_Needed";
  for (const Element* e : _elements) {
    os << ' ' << e->symbol();
  }
  os << " hits_needed ";
  return _hits_needed.debug_print(os);
}

int
Elements_Needed::add_atomic_number(int atomic_number) {
  const Element* e = get_element_from_atomic_number(atomic_number);
  if (e == nullptr) {
    cerr << "Elements_Needed::add_atomic_number:unrecognised atomic number " << atomic_number << '\n';
    return 0;
  }

  _elements.add(e);
  return 1;
}

int
Elements_Needed::matches(Query_Atoms_Matched & qam) const
{
  if (_elements.empty()) {
    cerr << "Elements_Needed::matches: no elements\n";
    iwabort();
    return 0;
  }

  const int number_elements = _elements.number_elements();
  if (number_elements == 1) {
    return _matches_single_element(qam);
  }

  // How often is each element found in the matched atoms.
  int * found = new_int(number_elements); std::unique_ptr<int[]> free_found(found);

  int nhits = 0;
  for (const Substructure_Atom * a : qam) {
    if (! a->include_in_embedding())
      continue;

    const int ndx = _elements.index(a->current_hold_atom()->element());
    if (ndx < 0) {
      continue;
    }
    ++nhits;
    ++found[ndx];
  }

  return _determine_match(found, nhits);
}

// For each of the elements we are looking for, the `found` array
// contains the number of times that element has been encountered.
// `nhits` is the number of matched atoms found.
// Depending on _operator, return whether or not this is a match.
int
Elements_Needed::_determine_match(const int * found, const int nhits) const {
  const int number_elements = _elements.number_elements();

  const int elements_matched = count_non_zero_occurrences_in_array(found, number_elements);

  // If and, all elements must have been matched.
  if (_operator == Op::AND) {
    if (elements_matched != number_elements) {  // Each element must be matched.
      return 0;
    }
    return _hits_needed.matches(nhits);
  } 

  if (_operator == Op::OR) {
    if (elements_matched == 0) {  // Nothing matched
      return 0;
    }
    return _hits_needed.matches(nhits);
  } 

  if (_operator == Op::XOR) {
    if (elements_matched != 1) {
      return 0;
    }
    return _hits_needed.matches(nhits);
  }

  cerr << "Elements_Needed::_determine_match:operator not set\n";
  return 0;
}

int
Elements_Needed::_matches_single_element(const Query_Atoms_Matched& qam) const {
  const Element* e0 = _elements[0];

  int nhits = 0;
  for (const Substructure_Atom * a: qam) {
    if (! a->include_in_embedding()) {
      continue;
    }
    if (a->current_hold_atom()->element() == e0) {
      ++nhits;
    }
  }

  return _hits_needed.matches(nhits);
}

int
Elements_Needed::matches(Molecule_to_Match & target_molecule) const
{
  if (_elements.number_elements() == 1) {
    return _matches_single_element(target_molecule);
  }

  const int number_elements = _elements.number_elements();
  int * found = new_int(number_elements); std::unique_ptr<int[]> free_found(found);
  int nhits = 0;
  const Molecule* m = target_molecule.molecule();
  for (const Atom * a : *m) {
    const int ndx = _elements.index(a->element());
    if (ndx >= 0) {
      found[ndx]++;
      ++nhits;
    }
  }

  return _determine_match(found, nhits);
}

int
Elements_Needed::_matches_single_element(Molecule_to_Match & target_molecule) const
{
  const Molecule& m = *target_molecule.molecule();
  const Element* e0 = _elements[0];
  int nhits = 0;
  for (const Atom* a : m) {
    if (a->element() == e0) {
      ++nhits;
    }
  }
  return _hits_needed.matches(nhits);
}
