#include <stdlib.h>
#include <iomanip>
using namespace std;

#include "misc.h"
#include "mdl_molecule.h"

ISIS_Atom_List::ISIS_Atom_List ()
{
  return;
}

bool
ISIS_Atom_List::operator== (const ISIS_Atom_List & rhs) const
{
  if (_normal_list != rhs._normal_list)
    return false;

  int n = _element.number_elements ();

  if (n != rhs._element.number_elements ())
    return false;

  for (int i = 0; i < n; i++)
  {
    const Element * e = _element[i];

    if (! rhs._element.contains (e))
      return false;
  }
  return true;
}

static int
fetch_int_value (const IWString & buffer,
                 int zfrom,
                 int zto,
                 int & zresult)
{
  const_IWSubstring token;

  buffer.from_to (zfrom, zto, token);

//cerr << "Token " << zfrom << " to " << zto << " in '" << buffer << "' is '" << token << "'\n";

  token.strip_leading_blanks ();

  return token.numeric_value (zresult);
}

int
ISIS_Atom_List::create_from_ALS_record (const IWString & buffer)
{
  assert (buffer.starts_with ("M  ALS "));

//if (buffer.length () < 24)
  if (buffer.length () < 20)   // ISIS Draw will create atom lists with just one member!
  {
    cerr << "ISIS_Atom_List::create_from_ALS_record: must be at least 24 chars\n";
    return 0;
  }

// M  ALS   1  2 F O   S   

  atom_number_t not_used_here;

  if (! fetch_int_value (buffer, 7, 9, not_used_here) || not_used_here < 1)
  {
    cerr << "ISIS_Atom_List::create_from_ALS_record: invalid atom number designation\n";
    return 0;
  }

  int n;

  if (! fetch_int_value (buffer, 10, 12, n) || n < 1)
  {
    cerr << "ISIS_Atom_List::create_from_ALS_record: invalid number of entries\n";
    return 0;
  }

  char tf = buffer[14];

  if ('F' == tf)
    _normal_list = 1;
  else if ('T' == tf)
    _normal_list = 0;
  else
  {
    cerr << "ISIS_Atom_List::create_from_ALS_record: must be T or F '" << buffer[14] << "'\n";
    return 0;
  }

  _element.resize (n);

  for (int i = 0; i < n; i++)
  {
    const_IWSubstring token;
    buffer.from_to (16 + i * 4, 16 + i * 4 + 3, token);

    token.strip_leading_blanks ();
    token.strip_trailing_blanks ();

    if (0 == token.length ())
    {
      cerr << "ISIS_Atom_List::create_from_ALS_record: empty token '" << buffer << "'\n";
      cerr << " n = " << n << " i = " << i << " from " << (16 + i * 4) << " to " << (16 + i * 4 + 3) << endl;
      return 0;
    }

    const Element * e = get_element_from_symbol_no_case_conversion(token);

    if (NULL == e)   // should not happen
    {
      cerr << "ISIS_RXN_FILE_Molecule::create_from_ALS_record: invalid element '" << token << "'\n";
      return 0;
    }

    _element.add (e);
  }

  return 1;
}

int
ISIS_Atom_List::debug_print (ostream & os) const
{
  os << "Atom list: normal " << _normal_list << endl;

  for (int i = 0; i < _element.number_elements (); i++)
  {
    const Element * e = _element[i];

    os << ' ' << e->atomic_number();
  }

  os << '\n';

  return os.good ();
}

int
ISIS_Atom_List::write_M_ALS (atom_number_t zatom,
                             ostream & output) const
{
  if (0 == _element.number_elements())
    return 1;

  output << "M  ALS" << setw (4) << (zatom + 1) << setw (3) << _element.number_elements ();

  if (_normal_list)
    output << " F";
  else
    output << " T";

  output << ' ';

  for (int i = 0; i < _element.number_elements (); i++)
  {
    const Element * e = _element[i];

    assert (NULL != e);

    output << e->symbol ();
    for (int j = e->symbol ().length (); j < 4; j++)
    {
      output << ' ';
    }
  }

  output << '\n';

  return output.good ();
}

/*
  'A' means any atom
*/

int
ISIS_Atom_List::initialise_from_mdl_A_symbol ()
{
  _normal_list = 1;

  initialise_from_mdl_Q_symbol();

  _element.add(get_element_from_atomic_number(6));    // C

  return 1;
}

/*
  Marvin adds AH, which is any atom, including Hydrogen
*/

int
ISIS_Atom_List::initialise_from_mdl_AH_symbol ()
{
  _normal_list = 1;

  initialise_from_mdl_Q_symbol();

  _element.add(get_element_from_atomic_number(6));    // C
  _element.add(get_element_from_atomic_number(1));    // H

  return 1;
}

/*
  Really don't like this, I should do something better with the Q symbol. Fix sometime...
*/

int
ISIS_Atom_List::initialise_from_mdl_Q_symbol ()
{
  _normal_list = 1;

  _element.add(get_element_from_atomic_number(7));    // N
  _element.add(get_element_from_atomic_number(8));    // O
  _element.add(get_element_from_atomic_number(9));    // F
  _element.add(get_element_from_atomic_number(15));   // P
  _element.add(get_element_from_atomic_number(16));   // S
  _element.add(get_element_from_atomic_number(17));   // Cl
  _element.add(get_element_from_atomic_number(35));   // Br
  _element.add(get_element_from_atomic_number(53));   // I

//cerr << "Initialised atom list, contains " << _element.number_elements() << " elements\n";

  return 1;
}

int
ISIS_Atom_List::initialise_from_mdl_QH_symbol ()
{
  initialise_from_mdl_Q_symbol();

  _element.add(get_element_from_atomic_number(1));

  return 1;
}

int
ISIS_Atom_List::convert_not_atom_lists_to_organic_lists ()
{
  int n = _element.number_elements();

  if (0 == n)
    return 0;

  if (_normal_list)
    return 0;

  resizable_array<atomic_number_t> current_set;

  for (int i = 0; i < n; i++)
  {
    current_set.add(_element[i]->atomic_number());
  }

  _element.resize_keep_storage(0);

  if (! current_set.contains(6))
    _element.add(get_element_from_atomic_number(6));
  if (! current_set.contains(7))
    _element.add(get_element_from_atomic_number(7));
  if (! current_set.contains(8))
    _element.add(get_element_from_atomic_number(8));
  if (! current_set.contains(9))
    _element.add(get_element_from_atomic_number(9));
  if (! current_set.contains(14))
    _element.add(get_element_from_atomic_number(14));
  if (! current_set.contains(15))
    _element.add(get_element_from_atomic_number(15));
  if (! current_set.contains(16))
    _element.add(get_element_from_atomic_number(16));
  if (! current_set.contains(17))
    _element.add(get_element_from_atomic_number(17));
  if (! current_set.contains(35))
    _element.add(get_element_from_atomic_number(35));
  if (! current_set.contains(53))
    _element.add(get_element_from_atomic_number(53));

  _normal_list = 1;

  return 1;
}

// M  V30 4 NOT[O,S,F,Cl,Br,I] 0.65 -1.80417 0 0

int
ISIS_Atom_List::initialise_atom_list_from_symbol (const const_IWSubstring & s)
{
  const_IWSubstring mys(s);

  if (mys.starts_with("NOT["))
  {
    mys.remove_leading_chars(4);
    _normal_list = 0;
  }
  else if (mys.starts_with('['))
  {
    _normal_list = 1;
    mys.remove_leading_chars(1);
  }

  if (! mys.ends_with(']'))
  {
    cerr << "MDL_Molecule::_convert_symbol_to_element:not sure what to do with '" << s << "'\n";
    return 0;
  }

  mys.chop();

  if (mys.length() > 1)
    ;
  else if ('A' == mys)
    return initialise_from_mdl_A_symbol();
  else if ('Q' == mys)
    return initialise_from_mdl_Q_symbol();

  int i = 0;
  const_IWSubstring token;
  while (mys.nextword(token, i, ','))
  {
    const Element * e = get_element_from_symbol_no_case_conversion(token);
    if (NULL == e)
    {
      cerr << "ISIS_Atom_List::initialise_atom_list_from_symbol:unrecognised element '" << token << "'\n";
      return 0;
    }

    _element.add(e);
  }

//cerr << "Atom list contains " << _element.number_elements() << " elements, list type " << _normal_list << endl;

  return _element.number_elements();
}

int
ISIS_Atom_List::invert_to_normal_list_based_on_organics ()
{
  int include_in_list[54];    // Iodine, 53 is the last organic

  set_vector(include_in_list, 54, 0);

  include_in_list[6] = 1;
  include_in_list[7] = 1;
  include_in_list[8] = 1;
  include_in_list[9] = 1;
  include_in_list[15] = 1;
  include_in_list[16] = 1;
  include_in_list[17] = 1;
  include_in_list[35] = 1;
  include_in_list[53] = 1;

  for (int i = 0; i < _element.number_elements(); i++)
  {
    atomic_number_t z = _element[i]->atomic_number();

    if (z < 0)
      ;
    else if (z > 53)
      ;
    else
      include_in_list[z] = 0;
  }

  _element.resize_keep_storage(0);

  for (int i = 0; i < 54; i++)
  {
    if (include_in_list[i])
      _element.add(get_element_from_atomic_number(i));
  }

  _normal_list = 1;

  return  _element.number_elements();
}
