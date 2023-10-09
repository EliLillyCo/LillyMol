#include <stdlib.h>
#include <iostream>

#include "mdl.h"
#include "mdl_atom_record.h"

using std::cerr;
using std::endl;

MDL_Atom_Record::MDL_Atom_Record()
{
  _default_values();

  return;
}

void
MDL_Atom_Record::_default_values()
{
  _msdiff = 0;
  _chg = 0;
  _astere = 0;
  _hcount = 0;
  _stereo_care = 0;
  _valence = 0;
  _h0designator = 0;
  _atom_map = 0;
  _inversion = 0;
  _exact_change = 0;

  return;
}

/*
  In the interests of efficiency, we do not zero out the existing values when we start
  reading. In general, if there are the same number of columns from one record to the
  next, this will be fine. But if someone tries to feed first 16 columns, then 6 columns,
  that will create problems!
*/

int
MDL_Atom_Record::build(const const_IWSubstring & buffer)
{
  if (buffer.length() < 34)
  {
    cerr << "MDL_Atom_Record::build:record must have at least 34 characters '" << buffer << "'\n";
    return 0;
  }

  int nw = buffer.nwords();

  if (nw < 4)
  {
    cerr << "MDL_Atom_Record::build:record must have at least 4 tokens '" << buffer << "'\n";
    return 0;
  }

  _astere = 0;    // since that column is optional

  const_IWSubstring token;

  buffer.from_to(0, 9, token);
  token.strip_leading_blanks();
  token.strip_trailing_blanks();

  if (! token.numeric_value(_x))
  {
    cerr << "Bad X value '" << token << "'\n";
    return 0;
  }

  (void) buffer.from_to(10, 19, token);
  token.strip_leading_blanks();
  token.strip_trailing_blanks();

  if (! token.numeric_value(_y))
  {
    cerr << "Bad Y value '" << token << "'\n";
    return 0;
  }

  (void) buffer.from_to(20, 29, token);
  token.strip_leading_blanks();
  token.strip_trailing_blanks();

  if (! token.numeric_value(_z))
  {
    cerr << "Bad Z value '" << token << "'\n";
    return 0;
  }

  int i = 30;
  (void) buffer.nextword(_atomic_symbol, i);

  if (! buffer.nextword(token, i))    // very unusual, but OK
  {
    _msdiff = 0;
    _chg = 0;
    _astere = 0;
    return 1;
  }

  if (! token.numeric_value(_msdiff))
  {
    cerr << "Bad msdiff value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))
    return 0;

  if (! token.numeric_value(_chg))
  {
    cerr << "Bad chg value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))
    return 1;

  if (! token.numeric_value(_astere))
  {
    cerr << "Bad astere value '" << token << "'\n";
    return 0;
  }

// At this stage, we are done with all the essential information

  if (! buffer.nextword(token, i))
    return 1;

  if (! token.numeric_value(_hcount))
  {
    cerr << "Bad hcount value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))    // if we got this far, we should have all the columns
    return 1;    // argh, I get sdf files with all kinds of variability, just ignore...

  if (! token.numeric_value(_stereo_care))
  {
    cerr << "Bad stereo care value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))    // if we got this far, we should have all the columns
  {
    _valence = 0;
    return 1;
  }

  if (! token.numeric_value(_valence))
  {
    cerr << "Bad valence value '" << token << "'\n";
    return 0;
  }

// At this stage, we need to make some decisions. If nw is 16, then all tokens are
// present. If nw is 13, then h0, and two following fields are absent
// I think I need to change this to a column specific form...

  if (nw >= 15)
  {
    if (! buffer.nextword(token, i))    // if we got this far, we should have all the columns
    {
//    cerr << "MDL_Atom_Record::build:not enough tokens, ignored\n";
      return 1;
    }

    if (! token.numeric_value(_h0designator))
    {
      cerr << "Bad H0 value '" << token << "'\n";
      return 0;
    }

    (void) buffer.nextword(token, i);    // skip over unused column
    (void) buffer.nextword(token, i);
//  cerr << "Last unused column '" << token << "'\n";
  }
  else
    _h0designator = 0;              // should not be necessary

  if (buffer.length() < 60)
    return 1;

  if (isdigit(buffer[60]))
    buffer.from_to(60, 62, token);
  else
    buffer.nextword(token, i);

// Beware three digit atom maps
//    1.5000    7.4952    0.0000 C   0  0  0  0  0  0  0  0  2126  0  0
//    2.2500    8.7942    0.0000 N   0  0  0  0  0  0  0  0  2127  0  0

// not properly handling this...
//   13.1012  -23.8165    0.0000 C   0  0  0  0  0  0  0  0 10101  0  0

  if (! token.numeric_value(_atom_map))
  {
    cerr << "Bad atom map value '" << token << "'\n";
    return 0;
  }

/*
  Beware atom maps above 99

   35.0815   -0.0000    0.0000 C   0  0  0  0  0  0  0  0 11100  0  0
         1         2         3         4         5         6         7         8         9
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
*/

  if (! buffer.nextword(token, i))    // last 3 columns must be present if we got to here
  {
    cerr << "MDL_Atom_Record::build:not enough tokens, ignored\n";
    cerr << buffer << endl;
    return 1;
  }

  if (! token.numeric_value(_inversion))
  {
    cerr << "Bad inversion value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))    // last 3 columns must be present if we got to here
  {
    cerr << "MDL_Atom_Record::build:not enough tokens, ignored\n";
    cerr << buffer << endl;
    return 1;
  }

  if (! token.numeric_value(_exact_change))
  {
    cerr << "Bad exact change value '" << token << "'\n";
    return 0;
  }

  return 1;
}

Atom *
MDL_Atom_Record::create_atom() const
{
  assert (_atomic_symbol.length() > 0);

  int is_radical = 0;

  int formal_charge = convert_from_mdl_charge(_chg);
  if (MDL_RADICAL == formal_charge)
  {
    is_radical = 1;
    formal_charge = 0;
  }
  

  const MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  Atom * a = mdlfos->create_mdl_atom(_atomic_symbol, _msdiff, formal_charge, is_radical);

  if (nullptr == a)
    return nullptr;

  a->setxyz(_x, _y, _z);

  return a;
}

MDL_Bond_Record::MDL_Bond_Record()
{
  _bond_type_read_in = 0;
//_directionality = 0;
  _bond_stereo = 0;
  _bond_topology = 0;
  _reacting_center_status = 0;

  return;
}

/*
  This is complicated by the fact that a bond record might look like
 59100  1  0  0  0  0
  where the atom numbers have merged.
*/

int
MDL_Bond_Record::build(const const_IWSubstring & buffer,
                        int natoms,
                        const MDL_File_Supporting_Material & mdlfos)
{
  if (buffer.length() < 9)
  {
    cerr << "MDL_Bond_Record::build:invalid record, too short '" << buffer << "'\n";
    return 0;
  }

  int nw = buffer.nwords();
  if (isdigit(buffer[3]))   // atoms have run together
    nw++;

  if (nw < 3)
  {
    cerr << "MDL_Bond_Record::build:invalid record, too few tokens '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring token;
  buffer.from_to(0, 2, token);
  token.strip_leading_blanks();

  if (! token.numeric_value(_a1) || _a1 < 1 || _a1 > natoms)
  {
    cerr << "MDL_Bond_Record::build:invalid a1 '" << buffer << "'\n";
    return 0;
  }

  _a1--;

  buffer.from_to(3, 5, token);
  token.strip_leading_blanks();

  if (! token.numeric_value(_a2) || _a2 < 1 || _a2 > natoms)   // check self bonds elsewhere
  {
    cerr << "MDL_Bond_Record::build:invalid a1 '" << buffer << "'\n";
    return 0;
  }

  _a2--;

  if (_a1 != _a2)   // great, the most common case
    ;
  else if (mdlfos.ignore_self_bonds())
    cerr << "MDL_Bond_Record::build:ignoring self bond, atom " << _a1 << endl;
  else
  {
    cerr << "MDL_Bond_Record::build:self bond, atom " << _a1 << endl;
    return 0;
  }
  
  buffer.from_to(7, 8, token);
  token.strip_leading_blanks();

  if (! token.numeric_value(_bond_type_read_in) || _bond_type_read_in < 1 || _bond_type_read_in > 8)
  {
    cerr << "MDL_Bond_Record::build:invalid bond '" << buffer << "'\n";
    return 0;
  }

  if (3 == nw || buffer.length() < 12)
    return 1;

  int i = 9;   // position the iterator just after the bond type

  if (! buffer.nextword(token, i))
    return 1;

  if (! token.numeric_value(_bond_stereo) || _bond_stereo < 0 || _bond_stereo > 6)
  {
    cerr << "MDL_Bond_Record::build:invalid bond stereo '" << buffer << "'\n";
    return 0;
  }

  if (4 == nw)
    return 1;

  if (7 == nw)     // column not used. Skip if present
    buffer.nextword(token, i);

  buffer.nextword(token, i);

  if (! token.numeric_value(_bond_topology) || _bond_topology < 0 || _bond_topology > 13)
  {
    cerr << "MDL_Bond_Record::build:invalid bond topology '" << buffer << "'\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (! token.numeric_value(_reacting_center_status) || _reacting_center_status < 0 || _reacting_center_status > 13)
  {
    cerr << "MDL_Bond_Record::build:invalid reacting centre status '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

/*int
MDL_Bond_Record::convert_mdl_bond_type_read_in_to_query (bond_type_t & for_query,
                                        bond_type_t & for_building_a_molecule) const
{
  for_building_a_molecule = SINGLE_BOND;

  if (1 == _bond_type_read_in)
    for_query = SINGLE_BOND;
  else if (2 == _bond_type_read_in)
  {
    for_query = DOUBLE_BOND;
    for_building_a_molecule = DOUBLE_BOND;
  }
  else if (3 == _bond_type_read_in)
  {
    for_query = TRIPLE_BOND;
    for_building_a_molecule = TRIPLE_BOND;
  }
  else if (4 == _bond_type_read_in)
    for_query = AROMATIC_BOND;
  else if (5 == _bond_type_read_in)
    for_query = SINGLE_BOND | DOUBLE_BOND;
  else if (6 == _bond_type_read_in)
    for_query = SINGLE_BOND | AROMATIC_BOND;
  else if (7 == _bond_type_read_in)
    for_query = DOUBLE_BOND | AROMATIC_BOND;
  else if (8 == _bond_type_read_in)
    for_query = SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND;
  else
  {
    cerr << "MDL_Bond_Record::convert_bond_type_read_in_to_query: unrecognised bond type " << _bond_type_read_in << endl;
    return 0;
  }

  return 1;
}*/

int
MDL_Bond_Record::bond_type_for_molecule (bond_type_t & for_building_a_molecule) const
{
  for_building_a_molecule = SINGLE_BOND;

  if (1 == _bond_type_read_in)
    return 1;
  else if (2 == _bond_type_read_in)
  {
    for_building_a_molecule = DOUBLE_BOND;
    return 1;
  }
  else if (3 == _bond_type_read_in)
  {
    for_building_a_molecule = TRIPLE_BOND;
    return 1;
  }

  return 1;
}

int
MDL_Bond_Record::bond_type_for_query(bond_type_t & for_query) const
{
  if (1 == _bond_type_read_in)
    for_query = SINGLE_BOND;
  else if (2 == _bond_type_read_in)
    for_query = DOUBLE_BOND;
  else if (3 == _bond_type_read_in)
    for_query = TRIPLE_BOND;
  else if (4 == _bond_type_read_in)
    for_query = AROMATIC_BOND;
  else if (5 == _bond_type_read_in)
    for_query = SINGLE_BOND | DOUBLE_BOND;
  else if (6 == _bond_type_read_in)
    for_query = SINGLE_BOND | AROMATIC_BOND;
  else if (7 == _bond_type_read_in)
    for_query = DOUBLE_BOND | AROMATIC_BOND;
  else if (8 == _bond_type_read_in)
    for_query = SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND;
  else
  {
    cerr << "MDL_Bond_Record::bond_type_for_query: unrecognised bond type " << _bond_type_read_in << endl;
    return 0;
  }

  return 1;
}

