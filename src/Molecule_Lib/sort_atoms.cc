// Supporting functions to enable sorting the Atoms in a Molecule

#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "molecule.h"

using std::cerr;
using std::endl;

/*
  Helper class for sorting atoms.
  In my implementations so far, I haven't needed the pointer to the Atom
*/

template <typename T>
class Atom_and_Info {
 private:
  const Atom* _atom;

  atom_number_t _initial_atom_number;

  //  the property on which we sort

  T _property;

 public:
  Atom_and_Info();

  void initialise(const Atom*, atom_number_t, T);

  const Atom*
  atom() const {
    return _atom;
  }

  atom_number_t
  initial_atom_number() const {
    return _initial_atom_number;
  }

  T
  property() const {
    return _property;
  }
};

template <typename T>
Atom_and_Info<T>::Atom_and_Info() {
  _atom = nullptr;

  _initial_atom_number = INVALID_ATOM_NUMBER;

  _property = {};

  return;
}

template <typename T>
void
Atom_and_Info<T>::initialise(const Atom* a, atom_number_t b, T prop) {
  _atom = a;

  _initial_atom_number = b;

  _property = prop;

  return;
}

template class Atom_and_Info<int>;

template <typename T>
int
atom_and_info_comparitor_int_larger(const void* aip1, const void* aip2) {
  const Atom_and_Info<T>* ai1 = (const Atom_and_Info<T>*)aip1;
  const Atom_and_Info<T>* ai2 = (const Atom_and_Info<T>*)aip2;

  int p1 = ai1->property();
  int p2 = ai2->property();

  if (p1 > p2) {
    return -1;
  }

  if (p1 < p2) {
    return 1;
  }

  // Keep the sort stable wrt initial atom number

  if (ai1->initial_atom_number() > ai2->initial_atom_number()) {
    return -1;
  }

  if (ai1->initial_atom_number() < ai2->initial_atom_number()) {
    return 1;
  }

  // Should never come here

  return 0;
}

template <typename T>
int
atom_and_info_comparitor_int_smaller(const void* aip1, const void* aip2) {
  return -atom_and_info_comparitor_int_larger<T>(aip1, aip2);
}

template <typename T>
int
Molecule::sort(const T* criterion, int direction) {
  assert(ok());
  assert(0 != direction);

  Atom_and_Info<T>* ai = new Atom_and_Info<T>[_number_elements];
  std::unique_ptr<Atom_and_Info<T>[]> free_ai(ai);

  for (int i = 0; i < _number_elements; i++) {
    Atom_and_Info<T>& aii = ai[i];
    aii.initialise(_things[i], i, criterion[i]);
  }

  if (direction > 0) {
    ::qsort(ai, _number_elements, sizeof(Atom_and_Info<T>),
            atom_and_info_comparitor_int_larger<T>);
  } else if (direction < 0) {
    ::qsort(ai, _number_elements, sizeof(Atom_and_Info<T>),
            atom_and_info_comparitor_int_smaller<T>);
  }

  // We need to now tell each atom what its new atom number is

  int* new_number = new int[_number_elements];
  std::unique_ptr<int[]> free_new_number(new_number);

  for (int i = 0; i < _number_elements; i++) {
    const Atom_and_Info<T>& aii = ai[i];

    new_number[aii.initial_atom_number()] = i;
  }

  return renumber_atoms(new_number);
}

template int Molecule::sort(const int*, int);
template int Molecule::sort(const int64_t*, int);

/*
  We do this a very inefficient way.
  We already have the function swap_atoms. Rather than writing something
  that works on the whole molecule at once, we just repeatedly call swap_atoms.
*/

int
Molecule::renumber_atoms(const int* new_number) {
  assert(ok());

  for (int i = 0; i < _number_elements; i++) {
    int j = new_number[i];

    if (j < 0 || j >= _number_elements) {
      cerr << "Molecule::renumber_atoms: i = " << i << " new number " << new_number[i]
           << " is invalid. Natoms = " << _number_elements << endl;
      abort();
      return 0;
    }
  }

  int* complete = new_int(_number_elements);
  std::unique_ptr<int[]> free_complete(complete);

  return _renumber_atoms(new_number, complete);
}

/*
  This code is kind of dense, but it does work.

  We focus our attention at the position P. Every time we find an atom
  there that doesn't belong at P, we move it to where it should be,
  mark that other location as complete.
*/

int
Molecule::_renumber_atoms(const int* new_number, int* complete) {
  for (int i = 0; i < _number_elements; i++) {
    if (complete[i]) {
      continue;
    }

    int atom_at_i = i;

    while (new_number[atom_at_i] != i) {
      int a2 = new_number[atom_at_i];

      swap_atoms(i, a2);  // atom A2 is now at position I

      complete[a2] = 1;  // A2 went where it should be

      atom_at_i = a2;  // we must moved atom A2 to our position

      //    cerr << "Atom " << a2 << " now at position " << i << " its number " <<
      //    new_number[a2] << endl;
    }
  }

  return 1;
}
