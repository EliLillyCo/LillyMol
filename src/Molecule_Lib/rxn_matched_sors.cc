#include <stdlib.h>

#include "iwreaction.h"

Matched_Atom_in_Component::Matched_Atom_in_Component ()
{
  _component = MAIC_CURRENT_COMPONENT;
  _matched_atom = -1;

  return;
}

/*
  We are parsing 'Rnn.'
  Extract the nn part
*/

/*static int
fetch_number (const_IWSubstring & token,
              int & zresult)
{
  int rc = 0;     // the number of characters consumed

  while (token.length ())
  {
    char c = token[0];
    if ('.' == c)
    {
      token++;
      return rc;
    }

    if (c < '0' || c > '9')    // huh
      return 0;

    rc++;
    zresult = zresult * 10 + c - '0';
    token++;
  }

  return 0;
}*/

int
Matched_Atom_in_Component::construct(const const_IWSubstring & token,
                                     int default_component)
{
//const_IWSubstring mytoken = token;    // a local copy

  _component = -1;      // scaffold is the default - so (1 2) gets recognised as 1 in the scaffold

  if (token.contains('.'))
  {
    const_IWSubstring s1, s2;
    token.split(s1, '.', s2);
    if (! s1.numeric_value(_component) || _component < 0)
    {
      cerr << "Matched_Atom_in_Component::construct:invalid component '" << token << "'\n";
      return 0;
    }

    if (! s2.numeric_value(_matched_atom) || _matched_atom < 0)
    {
      cerr << "Matched_Atom_in_Component::construct:invalid atom '" << token << "'\n";
      return 0;
    }

    return 1;
  }

  if (! token.numeric_value(_matched_atom) || _matched_atom < 0)
  {
    cerr << "Matched_Atom_in_Component::construct:invalid atom '" << token << "'\n";
    return 0;
  }

  _component = default_component;

  return 1;
}

int
Matched_Atom_in_Component::debug_print(std::ostream & os) const
{
  os << "Matched_Atom_in_Component:component " << _component << " atom " << _matched_atom << endl;

  return 1;
}

int
Matched_Atom_in_Component::adjust_matched_atoms_in_component(const extending_resizable_array<int> & xref)
{
  if (_component < 0)   // we must be the scaffold
    return 1;

  int newc = xref[_component];

//cerr << "Changing from " << _component << " to " << newc << ", atom " << _matched_atom << endl;

  assert (newc >= -1);    // -1 is special meaning for scaffold

  _component = newc;

  return 1;
}

/*
  This gets pretty goofy with the wierd number assigned to the scaffold. Redo this sometime!
*/

std::ostream &
operator << (std::ostream & os, const Matched_Atom_in_Component & masos)
{
  if (masos.in_scaffold())
    os << "0.";
  else if (masos.in_component() >= 0)
    os << (masos.in_component() + 1) << '.';
  else if (masos.in_current_component())
    ;
  else
    os << masos.in_component() << " huh!";

  os << masos.matched_atom();

  return os;
}
