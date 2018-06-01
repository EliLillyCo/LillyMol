#include <stdlib.h>

#include "substructure.h"

/*
  A link record from ISIS looks like
M  LIN  1   4   6   3   5
                    ^   ^
                    connected to
                ^
                repeats

            ^
            atom number
*/

int
ISIS_Link_Atom::construct_from_M_ISIS_record (const IWString & buffer)
{
  assert (buffer.starts_with ("M  LIN"));

  if (buffer.nwords () < 6)
  {
    cerr << "Link_Atom::construct_from_M_ISIS_record:must have at least 6 words\n";
    cerr << buffer << endl;
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword (token, i);    // M
  buffer.nextword (token, i);    // LIN
  buffer.nextword (token, i);    // sqeuential number, not of interest

  buffer.nextword (token, i);   // atom number

  if (! token.numeric_value (_a) || _a < 1)
  {
    cerr << "ISIS_Link_Atom::construct_from_M_ISIS_record:invalid central atom\n";
    return 0;
  }

  _a--;

  buffer.nextword (token, i);  //   range

  int range;
  if (token.numeric_value (range) && range > 0)
  {
    _d.set_max (range);
  }
  else
  {
    _d.reset ();

    if (_d.initialise (token))
    {
      for (int j = 0; j < _d.number_elements (); j++)    // convert from atom count to bond count
      {
        _d[j] = _d[j] + 1;
      }
    }
    else
    {
      cerr << "ISIS_Link_Atom::construct_from_M_ISIS_record:invalid range '" << token << "'\n";
      return 0;
    }
  }

  assert (_d.ok ());

  buffer.nextword (token, i);
  if (! token.numeric_value (_a1) || _a1 < 1)
  {
    cerr << "ISIS_Link_Atom::construct_from_M_ISIS_record:invalid a1 '" << token << "'\n";
    return 0;
  }

  _a1--;

  buffer.nextword (token, i);
  if (! token.numeric_value (_a2) || _a2 < 1)
  {
    cerr << "ISIS_Link_Atom::construct_from_M_ISIS_record:invalid a2 '" << token << "'\n";
    return 0;
  }

  _a2--;

  return 1;
}

void
ISIS_Link_Atom::_default_values ()
{
  _d.set_min (2);    // must be at least 1 atom in the link

  return;
}

ISIS_Link_Atom::ISIS_Link_Atom ()
{
  _default_values ();
  
  _a = -1;

  return;
}

/*
  Note, NO call to _default_values, it messes up things in _d, because it is
  done after the constructor initialisation
*/

ISIS_Link_Atom::ISIS_Link_Atom (const ISIS_Link_Atom & rhs) : Link_Atom (rhs)
{
  _a = rhs._a;

  return;
}

int
ISIS_Link_Atom::_adjust_atom_number (const int * xref,
                                     atom_number_t & a)
{
  if (INVALID_ATOM_NUMBER == a)     // strange,
    cerr << "ISIS_Link_Atom::_adjust_atom_number:atom not set!\n";
  else if (xref[a] == a)       // no change
    ;
  else if (xref[a] < 0)
  {
    cerr << "ISIS_Link_Atom::adjust_atom_numbers:atom " << a << " being eliminated!\n";
    return 0;
  }
  else
    a = xref[a];

  return 1;
}

int
ISIS_Link_Atom::adjust_atom_numbers(int const * xref)
{
  if (! _adjust_atom_number (xref, _a1))
    return 0;

  return _adjust_atom_number (xref, _a2);
}

int
ISIS_Link_Atom::debug_print (std::ostream & os) const
{
  os << "ISIS_Link_Atom between " << _a1 << " through " << _a << " to " << _a2 << " range " << 

  _d.write_compact_representation (os);

  return os.good ();
}

std::ostream &
operator<< (std::ostream & os, const ISIS_Link_Atom & l)
{
  l.debug_print (os);
  return os;
}

int
ISIS_Link_Atom::write_M_LIN (std::ostream & output) const
{
  output << "M  LIN";

  int m;
  _d.max (m);

  output << "  " << (_a + 1) << ' ' << m << ' ' << (_a1 + 1) << ' ' << (_a2 + 1) << endl;

  output << endl;

  return output.good ();
}

/*
  
*/

int
ISIS_Link_Atom::swap_atoms (atom_number_t n1, atom_number_t n2)
{
  if (_a == n1)
    _a = n2;

  return Link_Atom::swap_atoms (n1, n2);
}
