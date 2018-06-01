#include <stdlib.h>

#include "assert.h"

#include "cmdline.h"

#include "molecule.h"
#include "rmele.h"

void
Element_to_Remove::_default_values()
{
  _molecules_examined = 0;
  _molecules_changed = 0;
  _atoms_removed = 0;

  _maxcon_to_remove = -1;
  _add_bond_after_two_connected_removals = 0;

  _number_to_remove_per_molecule = 0;
}

Element_to_Remove::Element_to_Remove (const Element * e) : 
             Element_Matcher(e)
{
  _default_values();
}

Element_to_Remove::Element_to_Remove (atomic_number_t z) :
             Element_Matcher(z)
{
  _default_values();
}

Element_to_Remove::Element_to_Remove (const char * s) : 
             Element_Matcher(s)
{
  _default_values();
}

Element_to_Remove::Element_to_Remove (const IWString & s) : 
             Element_Matcher(s)
{
  _default_values();
}

Element_to_Remove::~Element_to_Remove()
{
  if (-14 == _molecules_examined)
    cerr << "Deleting already deleted Element_to_Remove\n";

  _molecules_examined = -14;

  return;
}

int
Element_to_Remove::ok() const
{
  if (_molecules_examined < 0)
    return 0;

  if (! Element_Matcher::ok())
    return 0;

  if (_molecules_changed < 0)
    return 0;

  if (_atoms_removed < 0)
    return 0;

  if (_atoms_removed > _molecules_changed)
    return 0;

  return 1;
}

int
Element_to_Remove::debug_print (std::ostream & os) const
{
  assert (os.good());

  (void) Element_Matcher::debug_print(os);
  os << _molecules_examined << " molecules examined, ";
  os << _molecules_changed << " molecules changed";
  if (_molecules_changed)
    os << ", " << _atoms_removed << " atoms removed";
  os << endl;

  return 1;
}

int
Element_to_Remove::report (std::ostream & os) const
{
  return debug_print(os);
}

int
Element_to_Remove::reset_counters()
{
  _default_values();

  return 1;
}

int
Element_to_Remove::_process (Molecule & m)
{
  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    const Element * e = m.elementi(i);

    int iso = m.isotope(i);

    if (! Element_Matcher::matches(e, iso))
      continue;

    int icon = m.ncon(i);
    if (_maxcon_to_remove < 0 ||
        (_maxcon_to_remove >= 0 && icon <= _maxcon_to_remove))
    {
      atom_number_t a1, a2;
      a1 = a2 = INVALID_ATOM_NUMBER;
      if (2 == icon && _add_bond_after_two_connected_removals)
      {
        a1 = m.other(i, 0);
        a2 = m.other(i, 1);
        if (a1 > i)
          a1--;
        if (a2 > i)
          a2--;
      }
      m.remove_atom(i);
      if (INVALID_ATOM_NUMBER != a1)
        m.add_bond(a1, a2, SINGLE_BOND);

      i--;
      matoms--;
      rc++;

      if (_number_to_remove_per_molecule && rc >= _number_to_remove_per_molecule)
        return rc;
    }
  }

  return rc;
}

/*
  Remove only those atoms for which ID == PROCESS_THESE[I]
*/

int
Element_to_Remove::_process (Molecule & m,
                             const int * process_these,
                             int id)
{
//cerr << "Element_to_Remove::_process:processing molecule with " << m.natoms() << " atoms\n";

  int rc = 0;
  for (int i = m.natoms() - 1; i >= 0; i--)    // scan backwards to preserve numbering
  {
    if (id != process_these[i])
      continue;

    int iso = m.isotope(i);

    const Element * e = m.elementi(i);

    if (! Element_Matcher::matches(e, iso))
      continue;

    cerr << "Element matcher matches atom " << i << endl;
    int icon = m.ncon(i);
    if (_maxcon_to_remove < 0 ||
        (_maxcon_to_remove >= 0 && icon <= _maxcon_to_remove))
    {
      atom_number_t a1, a2;
      a1 = a2 = INVALID_ATOM_NUMBER;
      if (2 == icon && _add_bond_after_two_connected_removals)
      {
        a1 = m.other(i, 0);
        a2 = m.other(i, 1);
        if (a1 > i)
          a1--;
        if (a2 > i)
          a2--;
      }

      m.remove_atom(i);
      if (INVALID_ATOM_NUMBER != a1)
        m.add_bond(a1, a2, SINGLE_BOND);

      rc++;

      if (_number_to_remove_per_molecule && rc >= _number_to_remove_per_molecule)
        return rc;
    }
  }

  return rc;
}

int
Element_to_Remove::process (Molecule & m)
{
  _molecules_examined++;

  int nr = _process(m);

  if (0 == nr)
    return 0;

  _molecules_changed++;
  _atoms_removed += nr;

  return nr;
}

int
Element_to_Remove::process (Molecule & m,
                            const int * process_these,
                            int id)
{
  _molecules_examined++;

  int nr = _process(m, process_these, id);

  if (0 == nr)
    return 0;

  _molecules_changed++;
  _atoms_removed += nr;

  return nr;
}

Elements_to_Remove::Elements_to_Remove()
{
  _remove_all_non_natural_elements = 0;
  _remove_all_isotopes = 0;

  return;
}

int
Elements_to_Remove::construct_from_command_line (Command_Line & cl,
                            int verbose,
                            char option)
{
  int i = 0;
  IWString ele;
  while (cl.value(option, ele, i++))
  {
/*  if ('*' == ele)
    {
      set_remove_all_non_natural_elements(1);
      continue;
    }*/

    if ("help" == ele)
    {
      cerr << "Removes elements, specify atomic symbol or atomic number\n";
      display_element_matcher_syntax(cerr);
      exit(0);
    }

    Element_to_Remove * tmp = new Element_to_Remove(ele);

    add(tmp);

    if (verbose)
      cerr << "Elements_to_Remove:: will remove atoms of type '" << ele << "'\n";
  }

  return _number_elements;
}

int
Elements_to_Remove::process (Molecule & m)
{
  int rc = 0;

#ifdef DEBUG_ELEMENTS_TO_REMOVE_PROCESS
  cerr << "Elements_to_Remove::process:has " << _number_elements << " items to check\n";
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i]->process(m);
  }

  if (_remove_all_non_natural_elements)
    rc += m.remove_all_non_natural_elements();

  return rc;
}

/*
  Sometimes the caller will want to process only certain atoms,
  those for which PROCESS_THESE[I] == ID

  Hmmm, this is more complex that it looks, because of
  different bonding, and possible interactions between
  atoms. Not implemented for now....
*/

/*int
Elements_to_Remove::process (Molecule & m,
                             const int * process_these,
                             const int id)
{
  resizable_array<atom_number_t> remove_these;

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i]->process(m, process_these, id, remove_these);
  }

  if (_remove_all_non_natural_elements)
    rc += m.remove_all_non_natural_elements();

  return rc;
}*/

int
Elements_to_Remove::report (std::ostream & os) const
{
  for (int i = 0; i < _number_elements; i++)
    _things[i]->report(os);

  return 1;
}

int
Elements_to_Remove::reset_counters()
{
  for (int i = 0; i < _number_elements; i++)
    _things[i]->reset_counters();

  return 1;
}
