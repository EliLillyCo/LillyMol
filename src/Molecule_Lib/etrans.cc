#include "etrans.h"

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "element.h"
#include "molecule.h"
#include "target.h"

using std::cerr;

int
display_standard_etrans_options(std::ostream& os, char cflag) {
  // clang-format off
  os << "  -" << cflag << " E1=E2       specify elemental transformation, all E1 become E2\n";
  os << "               Special elements 'organic', 'nonorganic' and 'nonperiodic' are recognised\n";
  os << "\n";
  os << "               append '.valence' and the valence of the resulting molecule will be checked.\n";
  os << "               If the valence is not valid, the element is reverted\n";
  // clang-format on

  return os.good();
}

void
Element_Transformation::_default_values() {
  _to = nullptr;

  _molecules_processed = 0;
  _molecules_changed = 0;
  _atoms_changed = 0;

  _reverse_bad_valence = 0;

  _transform_every_atom_type = 0;

  _isotope = 0;

  return;
}

Element_Transformation::Element_Transformation() {
  _default_values();
}

int
Element_Transformation::ok() const {
  if (_from.element() && nullptr == _to) {
    return 0;
  }

  return 1;
}

int
Element_Transformation::debug_print(std::ostream& os) const {
  assert(ok());

  if (_transform_every_atom_type) {
    os << "Will change all atoms";
  } else {
    os << "Will change atoms of type " << _from;
  }

  if (_to) {
    os << " to atoms of type " << _to->symbol();
    if (_isotope) {
      os << ", isotope " << _isotope;
    }
    os << '\n';
  }

  if (_molecules_processed) {
    os << "Have processed " << _molecules_processed << " molecules\n";
  }

  if (_molecules_processed) {
    os << "Changed " << _atoms_changed << " atoms in " << _molecules_changed
       << " molecules\n";
  }

  return 1;
}

int
Element_Transformation::build(const IWString& value) {
  IWString esource(value);

  if (esource.ends_with(",valence")) {
    _reverse_bad_valence = 1;
    esource.chop(8);
  }

  int i = esource.index('=');
  if (i <= 0 || i == esource.length() - 1) {
    cerr << "String for Element Transformation must be 'E1=E2'\n";
    return 0;
  }

  IWString ee = esource.from_to(0, i - 1);

  if (ee == "all") {
    _transform_every_atom_type = 1;
  } else if (!_from.construct_from_string(ee)) {
    cerr << "Cannot determine E1 from '" << esource << "': '" << ee << "'\n";
    return 0;
  }

  ee = substr(esource, i + 1);

  if (isalnum(ee[0])) {  // great, elements start with letters or numbers(isotopes)
    ;
  } else if ('*' == ee) {
    ;
  } else {
    cerr << "Element_Transformation::build: elements must start with a letter or number '"
         << ee << "'\n";
    return 0;
  }

  _to = get_element_from_symbol_no_case_conversion(ee);

  if (nullptr == _to) {
    _to = create_element_with_symbol(ee);
  }

  return 1;
}

int
Element_Transformation::process(Molecule& m) {
  assert(ok());
  assert(m.ok());

  _molecules_processed++;

  int rc = 0;

  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    const Element* e = m.elementi(i);

    if (e == _to) {  // it is already what we are changing to (grammar?)
      if (_isotope) {
        m.set_isotope(i, _isotope);
        rc++;
      }
      continue;
    }

    if (_transform_every_atom_type || _from.matches(e)) {
      rc += MaybeChangeElement(m, i);
    }
  }

  if (rc) {
    _molecules_changed++;
    _atoms_changed += rc;
  }

  return rc;
}

int
Element_Transformation::MaybeChangeElement(Molecule& m, atom_number_t zatom) {
  const Element* from = m.elementi(zatom);
  m.set_element(zatom, _to);
  if (_reverse_bad_valence && ! m.valence_ok(zatom)) {
    m.set_element(zatom, from);
    return 0;
  }

  if (_isotope) {
    m.set_isotope(zatom, _isotope);
  }

  return 1;
}

int
Element_Transformation::process(Molecule_to_Match& m) {
  assert(ok());
  assert(m.ok());

  _molecules_processed++;

  int rc = 0;

  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    Target_Atom& a = m[i];

    const Element* e = a.element();

    if (e == _to)  // it is already what we are changing to (grammar?)
    {
      if (_isotope) {
        a.set_isotope(_isotope);
        rc++;
      }
      continue;
    }

    if (_transform_every_atom_type || _from.matches(e)) {
      a.set_element(_to);
      if (_isotope) {
        a.set_isotope(_isotope);
      }
      rc++;
    }
  }

  if (rc) {
    _molecules_changed++;
    _atoms_changed += rc;
  }

  return rc;
}

int
process_element_transformations(Command_Line& cl,
                                Element_Transformations& element_transformations,
                                int verbose, char eflag) {
  IWString c;
  int i = 0;
  int rc = 0;
  while (cl.value(eflag, c, i++)) {
    if (c.starts_with("HALOGEN")) {
      Element_Transformation* tmp = new Element_Transformation;
      tmp->build("I=Cl");
      element_transformations.add(tmp);
      tmp = new Element_Transformation;
      tmp->build("Br=Cl");
      element_transformations.add(tmp);
      continue;
    }

    Element_Transformation* tmp = new Element_Transformation;
    if (!tmp->build(c)) {
      cerr << "Cannot process option '" << eflag << "' number " << i << " '" << c
           << "'\n";
      return 0;
    }

    if (verbose) {
      tmp->debug_print(cerr);
    }

    element_transformations.add(tmp);
    rc++;
  }

  return rc + 1;
}

int
Element_Transformations::construct_from_command_line(Command_Line& cl, int verbose,
                                                     char eflag) {
  IWString c;
  for (int i = 0; cl.value(eflag, c, i); ++i) {;
    if ("help" == c) {
      display_standard_etrans_options(cerr, eflag);
      exit(2);
    }

    std::unique_ptr<Element_Transformation> tmp =
        std::make_unique<Element_Transformation>();
    if (!tmp->build(c)) {
      cerr << "Cannot process option '" << eflag << "' number " << i << " '" << c
           << "'\n";
      return 0;
    }

    if (verbose) {
      tmp->debug_print(cerr);
    }

    add(tmp.release());
  }

  return _number_elements;
}

int
Element_Transformations::Add(const IWString& directive) {
  std::unique_ptr<Element_Transformation> tmp =
      std::make_unique<Element_Transformation>();
  if (!tmp->build(directive)) {
    cerr << "Element_Transformations::Add:cannot interpret '" << directive << "'\n";
    return 0;
  }

  add(tmp.release());

  return 1;
}

int
Element_Transformations::ok() const {
  for (int i = 0; i < _number_elements; i++) {
    if (!_things[i]->ok()) {
      return 0;
    }
  }

  return 1;
}

int
Element_Transformations::debug_print(std::ostream& os) const {
  assert(os.good());

  os << "Info on " << _number_elements << " element transformation(s)\n";
  int rc = 0;
  for (int i = 0; i < _number_elements; i++) {
    rc += _things[i]->debug_print(os);
  }

  return rc;
}

int
Element_Transformations::process(Molecule* m) {
  assert(nullptr != m);

  return process(*m);
}

int
Element_Transformations::process(Molecule& m) {
  int rc = 0;
  for (int i = 0; i < _number_elements; i++) {
    rc += _things[i]->process(m);
  }

  return rc;
}

int
Element_Transformations::process(Molecule_to_Match& m) {
  int rc = 0;
  for (int i = 0; i < _number_elements; i++) {
    rc += _things[i]->process(m);
  }

  return rc;
}
