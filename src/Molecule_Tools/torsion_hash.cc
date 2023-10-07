#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "torsion_hash.h"

using std::cerr;
using std::endl;

int
Torsion_Hash::add (const IWString & s)
{
  IW_STL_Hash_Map_int::const_iterator f = IW_STL_Hash_Map_int::find (s);

  if (f != end ())     // already present
    return 0;

  int n = size ();

  insert (IW_STL_Hash_Map_int::value_type (s, n));

  return n + 1;    // something definitely non-zero
}

/*
  Return a unique identifier for a string. Stores the string if not present
*/

int
Torsion_Hash::unique_identifier (const IWString & s)
{
  IW_STL_Hash_Map_int::const_iterator f = find (s);

  int rc;

  if (f != end ())     // already present, fetch the existing value
  {
    rc = (*f).second;
  }
  else
  {
    rc = size ();
    add (s);
  }

  return rc;
}

/*
  Check to see whether or not a string is in the archive. If so, set the
  value for ID
*/

int
Torsion_Hash::fetch_unique_identifier (const IWString & s, int & id) const
{
  IW_STL_Hash_Map_int::const_iterator f = find (s);

  if (end () == f)
    return 0;

  id = (*f).second;

  return 1;
}

int
Torsion_Hash::write (std::ostream & output) const
{
  output << "Stored = " << size () << endl;

  for (IW_STL_Hash_Map_int::const_iterator i = begin (); i != end (); ++i)
  {
    output << (*i).first << ' ' << (*i).second << endl;
  }

  return output.good ();
}

int
Torsion_Hash::do_write (IWString_and_File_Descriptor & output) const
{
  output << "Stored = " << static_cast<unsigned int>(size()) << '\n';

  for (IW_STL_Hash_Map_int::const_iterator i = begin (); i != end (); ++i)
  {
    output << (*i).first << ' ' << (*i).second << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good ();
}
