#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "misc2.h"
#include "substructure.h"
#include "target.h"

namespace down_the_bond {

using std::cerr;

constexpr char open_brace = '{';
constexpr char close_brace = '}';

DownTheBond::DownTheBond() {
  _a1 = -1;
  _a2 = -1;
}

DownTheBond::DownTheBond(int a1) {
  _a1 = a1;
  _a2 = -1;
}

// Various directives can be placed here.
//   a  total number of atoms seen down the bond.
//   u  total number of unmatched atoms seen downthe bond.
//   d  longest distance of any atom from the _a2.
//      smarts for things that must NOT be present??
// For now, only the 'a' directive is supported
// a2, a>3, a<4, a[4-9]
// Long term other directives would be supported via logical operators
// a9;d<4
int
DownTheBond::Build(const const_IWSubstring& buffer) {
  // Temporary, pending more functionality being added.
  if (buffer[0] != 'a') {
    cerr << "DownTheBond::Build:not an atom directive '" << buffer << "'\n";
    return 0;
  }
  const char * s = buffer.data();
  if (! substructure_spec::SmartsNumericQualifier(s + 1, buffer.length() - 1, _natoms)) {
    cerr << "DownTheBond::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
DownTheBond::Matches(Molecule_to_Match& target,
                      Query_Atoms_Matched & matched_atoms,
                      int * visited) {
  const Molecule& m = *target.molecule();

  const atom_number_t a1 = matched_atoms[_a1]->current_hold_atom()->atom_number();
  const atom_number_t a2 = matched_atoms[_a2]->current_hold_atom()->atom_number();
  // cerr << "atoms " << a1 << " and " << a2 << '\n';

  if (! m.ok_atom_number(a1) || ! m.ok_atom_number(a2)) {
    cerr << "DownTheBond::Matches:invalid atom number " << a1 << " or " << a2 << '\n';
    return 0;
  }

  resizable_array<atom_number_t> atom_stack;
  const Atom& atom2 = m.atom(a2);
  for (const Bond* b : atom2) {
    const atom_number_t j = b->other(a2);
    if (j == a1) {
      continue;
    }
    atom_stack << j;
  }

  if (atom_stack.empty()) {
    return _natoms.matches(1);
  }

  std::fill_n(visited, m.natoms(), 0);

  visited[a1] = 1;
  visited[a2] = 1;
  int number_visited = 1;
  while (! atom_stack.empty()) {
    const atom_number_t i = atom_stack.pop();
    if (visited[i]) {
      continue;
    }
    visited[i] = 1;
    ++number_visited;
    const Atom& atom = m.atom(i);
    for (const Bond* bond : atom) {
      atom_number_t j = bond->other(i);
      if (j == a1) {  // Must be part of a loop, must fail.
        return 0;
      }
      if (visited[j]) {
        continue;
      }
      atom_stack << j;
    }
  }

  //cerr << "visited " << number_visited << " atom, matches " << _natoms.matches(number_visited) << '\n';

  return _natoms.matches(number_visited);
}

}  // namespace down_the_bond

using down_the_bond::DownTheBond;

int
Single_Substructure_Query::_down_the_bond_satisfied(Molecule_to_Match& target,
                Query_Atoms_Matched& matched_atoms) const {
  std::unique_ptr<int[]> visited = std::make_unique<int[]>(target.natoms());

  for (DownTheBond * dtb : _down_the_bond) {
    if (! dtb->Matches(target, matched_atoms, visited.get())) {
      return 0;
    }
  }

  return 1;
}
