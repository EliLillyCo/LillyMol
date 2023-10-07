// Read and write xyz files

#include <iostream>
#include <memory>

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"

using std::cerr;

int
Molecule::read_molecule_xyz_ds(iwstring_data_source & input) {
  IWString buffer;
  if (! input.next_record(buffer)) {
    cerr << "Molecule::read_molecule_xyz_ds:eof\n";
    return 0;
  }

  int n;
  if (! buffer.numeric_value(n) || n < 0) {
    cerr << "Molecule::read_molecule_xyz_ds:invalid natoms '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(buffer)) {
    cerr << "Molecule::read_molecule_xyz_ds:cannot read second record\n";
    return 0;
  }

  buffer.strip_leading_blanks();
  buffer.strip_trailing_blanks();

  _molecule_name = buffer;

  resize(n);

  static constexpr char kTab = '\t';
  static constexpr char kSpace = ' ';

  // TO make things easy, we translate all inputs to space seprated,
  // but we do not know if we have tabs in the input. Currently comma
  // separated tokens not supported. Easy to do if needed.
  int need_to_translate_tabs = 0;

  for (int i = 0; i < n; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "Molecule::read_molecule_xyz_ds:cannot read record " << i << '\n';
      return 0;
    }
    if (i == 0) {
      need_to_translate_tabs = buffer.contains(kTab);
    } 
    if (need_to_translate_tabs) {
      buffer.gsub(kTab, kSpace);
    }

    if (! AddXYZAtom(buffer)) {
      cerr << "Molecule::read_molecule_xyz_ds:invalid record '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Molecule::AddXYZAtom(const const_IWSubstring& buffer) {
  int i = 0;
  const_IWSubstring token;
  if (! buffer.nextword(token, i)) {
    cerr << "Molecule::AddXYZAtom:empty\n";
    return 0;
  }

  const Element* e = get_element_from_symbol_no_case_conversion(token);
  if (e == nullptr) {
    cerr << "Molecule:AddXYZAtom:unrecognised element '" << token << "'\n";
    return 0;
  }

  static constexpr int kPartialMolecule = 1;

  std::unique_ptr<Atom> atom = std::make_unique<Atom>(e);

  for (int j = 0; j < 3; ++j) {
    if (! buffer.nextword(token, i)) {
      cerr << "Molecule::AddXYZAtom:truncated\n";
      return 0;
    }

    float value;
    if (! token.numeric_value(value)) {
      cerr << "Molecule:AddXYZAtom:invalid numeric\n";
      return 0;
    }
    if (j == 0) {
      atom->x() = value;
    } else if (j == 1) {
      atom->y() = value;
    } else if (j == 2) {
      atom->z() = value;
    }
  }

  add(atom.release(), kPartialMolecule);

  return 1;
}

int
Molecule::write_molecule_xyz(std::ostream& output) const {
  output << _number_elements << '\n';
  output << _molecule_name << '\n';

  constexpr char kSep = ' ';

  for (int i = 0; i < _number_elements; ++i) {
    const Atom* a = _things[i];

    output << a->atomic_symbol() << kSep << a->x() <<
              kSep << a->y()  << kSep << a->z() << '\n';
  }

  return output.good();
}
