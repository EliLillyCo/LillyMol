#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/data_source/string_data_source.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "atom_alias.h"

using std::cerr;
using std::endl;

Atom_Alias::Atom_Alias ()
{
  _atom = INVALID_ATOM_NUMBER;

  return;
}

Atom_Alias::Atom_Alias (const Atom_Alias & rhs)
{
  _copy (rhs);

  return;
}

Atom_Alias &
Atom_Alias::operator= (const Atom_Alias & rhs)
{
  _copy (rhs);

  return *this;
}

void
Atom_Alias::_copy (const Atom_Alias & rhs)
{
  _atom = rhs._atom;
  _alias = rhs._alias;

  return;
}

template <typename T>
int
Atom_Alias::build (const const_IWSubstring & buffer,
                   T & input)
{
  assert (buffer.starts_with ("A  "));

  const_IWSubstring tmp (buffer);

  tmp.remove_leading_chars (3);

  tmp.strip_leading_blanks ();

  if (! tmp.numeric_value (_atom) || _atom < 1)
  {
    cerr << "Atom_Alias::build: invalid atom number specification '" << buffer << "'\n";
    return 0;
  }

  _atom--;

  if (! input.next_record (_alias))
  {
    cerr << "Atom_Alias::build:premature EOF\n";
    return 0;
  }

  return 1;
}

template class resizable_array_p<Atom_Alias>;
template class resizable_array_base<Atom_Alias *>;

template int Atom_Alias::build<iwstring_data_source>(const_IWSubstring const&, iwstring_data_source&);
template int Atom_Alias::build<String_Data_Source>(const_IWSubstring const&, String_Data_Source&);
