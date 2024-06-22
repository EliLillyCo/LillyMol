#include <stdlib.h>
#include <algorithm>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "molecule.h"
#include "allowed_elements.h"

using std::cerr;

void
Allowed_Elements::_default_values()
{
  set_vector(_allowed_element, HIGHEST_ATOMIC_NUMBER + 1, 0);

  _allowed_element[1] = 1;
  _allowed_element[6] = 1;
  _allowed_element[7] = 1;
  _allowed_element[8] = 1;
  _allowed_element[9] = 1;
  _allowed_element[15] = 1;
  _allowed_element[16] = 1;
  _allowed_element[17] = 1;
  _allowed_element[35] = 1;
  _allowed_element[53] = 1;
  _allowed_element[3]  = 1;     // Li
  _allowed_element[11] = 1;    // Na
  _allowed_element[12] = 1;    // Mg
  _allowed_element[19] = 1;    // K
  _allowed_element[20] = 1;    // Ca

  return;
}

Allowed_Elements::Allowed_Elements()
{
  _default_values();
}

void
Allowed_Elements::reset_to_defaults()
{
  _default_values();

  return;
}

int
Allowed_Elements::build_from_command_line(Command_Line & cl,
                                          char flag,
                                          int verbose)
{
  int i = 0;
  IWString e;

  while (cl.value(flag, e, i++))
  {
    const Element * o = get_element_from_symbol_no_case_conversion(e);
    if (nullptr == o)
    {
      cerr << "Sorry, non periodic table element '" << e << "', cannot be OK\n";
      return 0;
    }

    atomic_number_t z = o->atomic_number();

//  cerr << "Processing element '" << e << "', z = " << z << '\n';

    assert(z >= 0 && z <= HIGHEST_ATOMIC_NUMBER);

    _allowed_element[z] = 1;

    if (verbose)
      cerr << "Element " << e << " atomic number " << z << " allowed\n";
  }

  for (int i = 1; i <= HIGHEST_ATOMIC_NUMBER; i++)   // should we skip over Hydrogen or not?
  {
    if (_allowed_element[i])
      continue;
  
    const Element * e = get_element_from_atomic_number(i);
    if (! e->organic())
      continue;

    _allowed_element[i] = 1;
    if (verbose)
      cerr << "Organic element " << e->symbol() << " atomic number " << i << " allowed\n";
  }

  return 1;
}

int
Allowed_Elements::contains_non_allowed_atoms(const Molecule & m) const
{
  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    const Element * e = m.elementi(i);

    if (! e->is_in_periodic_table())
      return 1;

    if (! _allowed_element[e->atomic_number()])
      return 1;
  }

  return 0;   // all these elements are OK
}

void
Allowed_Elements::set_allow(atomic_number_t z, int s)
{
  _allowed_element[z] = s;

  return;
}

void
Allowed_Elements::exclude_metals() {
  for (int i = 0; i <= HIGHEST_ATOMIC_NUMBER; ++i) {
    const Element* e = get_element_from_atomic_number(i);
    if (e->is_metal()) {
      _allowed_element[i] = 0;
    }
  }
}

void
Allowed_Elements::Clear() {
  std::fill_n(_allowed_element, HIGHEST_ATOMIC_NUMBER + 1, 0);
}
