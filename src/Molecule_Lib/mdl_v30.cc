#include <stdlib.h>


#define COMPILING_MDL_CC

#include "mdl.h"
#include "molecule.h"
#include "aromatic.h"
#include "chiral_centre.h"
#include "iwstring_data_source.h"
#include "string_data_source.h"
#include "readmdl.h"

int
Molecule::_fill_empty_molecule_with_null_atoms (int na)
{
  assert (0 == _number_elements);

  resize (na);

// In general, it is a very bad idea to put NULL pointers into a resizable_array_p,
// but we do it as a check to make sure all atoms are assigned. It will be a
// fatal error if a NULL is found anyway

  for (int i = 0; i < na; i++)
  {
    resizable_array_p<Atom>::add (NULL);
  }

  return na;
}

/*
  V30 records may have continuation lines, indicated by a '-' at the end
  of the record
*/

template <typename T>
int
read_next_v30_record(T & input,
                     IWString & buffer)
{
  if (! input.next_record (buffer))
    return 0;

  while (buffer.ends_with (" _"))
  {
    const_IWSubstring cline;

    if (! input.next_record (cline))
    {
      cerr << "Premature eof in V30 record after continuation\n";
      return 0;
    }

    if (! cline.starts_with("M  V30"))
    {
      cerr << "V30 continuation line doesn't start with 'M  V30', line " << input.lines_read() << endl;
      cerr << buffer;
      return  0;
    }

    buffer.chop(1);    // remove the dash

    cline.remove_leading_chars(7);    // remove "M  V30 "

    buffer += cline;
  }

  return 1;
}

//template int read_next_v30_record(String_Data_Source &, IWString &);
template int read_next_v30_record<iwstring_data_source>(iwstring_data_source&, IWString&);
template int Molecule::_read_mdl_atom_connection_table_v30<iwstring_data_source>(iwstring_data_source&, int&, MDL_File_Supporting_Material&);
template int Molecule::_read_mdl_atom_connection_table_v30<String_Data_Source>(String_Data_Source&, int&, MDL_File_Supporting_Material&);
template int Molecule::_read_v30_bond_list<iwstring_data_source>(iwstring_data_source&, int, int*, int*);
template int Molecule::_read_v30_bond_list<String_Data_Source>(String_Data_Source&, int, int*, int*);
template int read_next_v30_record<String_Data_Source>(String_Data_Source&, IWString&);

template <typename T>
int
Molecule::_read_mdl_atom_connection_table_v30(T & input,
                                              int & nb, 
                                              MDL_File_Supporting_Material & mdlfos)
{
  input.set_strip_trailing_blanks(1);

  IWString buffer;
  if (! read_next_v30_record(input, buffer))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: cannot read first V30 file record\n";
    return 0;
  }

  if ("M  V30 BEGIN CTAB" != buffer)
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: possibly bad BEGIN CTAB record\n";
    cerr << buffer << endl;
    cerr << "Ignore possible error, continuing....\n";
  }

  if (! read_next_v30_record(input, buffer))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: cannot read V30 counts record\n";
    return 0;
  }

  if (! buffer.starts_with("M  V30 COUNTS "))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: this should be the counts record, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  const_IWSubstring token;
  if (! buffer.word(3, token))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: cannot extract atom count from count record\n";
    cerr << buffer << endl;
    return 0;
  }

  int na;
  if (! token.numeric_value(na) || na < 0)
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: invalid NA value '" << token << "' line " << input.lines_read() << endl;
    return 0;
  }

  if (_elements_allocated < na)
    _fill_empty_molecule_with_null_atoms(na);

  if (! buffer.word(4, token))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: cannot extract bond count from count record\n";
    cerr << buffer << endl;
    return 0;
  }

  if (! token.numeric_value(nb) || nb < 0)
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: invalid NB value '" << token << "' line " << input.lines_read() << endl;
    return 0;
  }

  if (! read_next_v30_record(input, buffer))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: premature EOF\n";
    return 0;
  }

// skip down (if necessary) until we get to the start of the atom table

  while ("M  V30 BEGIN ATOM" != buffer)
  {
    if (! read_next_v30_record(input, buffer))
    {
      cerr << "Molecule::_read_mdl_atom_connection_table_v30: premature EOF\n";
      return 0;
    }
  }

  for (int i = 0; i < na; i++)
  {
    if (! read_next_v30_record(input, buffer))
    {
      cerr << "Molecule::_read_mdl_atom_connection_table_v30: premature eof\n";
      return 0;
    }

    if (! _parse_v30_atom_record(buffer, 1, mdlfos))
    {
      cerr << "Molecule::_read_mdl_atom_connection_table_v30: invalid atom record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

// There should be 'END ATOM' at the end of the atom set

  if (! read_next_v30_record(input, buffer))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: premature eof\n";
    return 0;
  }

  if ("M  V30 END ATOM" != buffer)
  {
    cerr << "Molecule::_read_mdl_atom_connection_table_v30: hmmm, atom set not followed by 'END ATOM'\n";
    cerr << buffer << endl;
    return 0;
  }

// We need to make sure that all atoms have been assigned

  assert (_number_elements == na);

  for (int i = 0; i < _number_elements; i++)
  {
    if (NULL == _things[i])
    {
      cerr << "Molecule::_read_mdl_atom_connection_table_v30: atom " << i << " not assigned!!!\n";
      return 0;
    }
  }

  return 1;
}

/*
  An atom record looks like

  M  V30 index type x y z CHG= RAD= MASS= CFG=

  M  V30 4 NOT[O,S,F,Cl,Br,I] 0.65 -1.80417 0 0
*/

int
Molecule::_parse_v30_atom_record(const IWString & buffer,
                                 int convert_symbol_to_element,
                                 MDL_File_Supporting_Material & mdlfos)
{
// The mandatory parts of the record

  int ndx;
  const_IWSubstring zsymbol;
  float x, y, z;

// These are all optional

  int chg = 0;
  int mass = 0;
  int rad = 0;
  int cfg = 0;

  int ntokens = 0;

  int i = 7;    // start fetching tokens at this position

  const_IWSubstring token;

  while (buffer.nextword(token, i))
  {
    ntokens++;

    if (1 == ntokens)     // ndx value
    {
      if (! token.numeric_value(ndx) || ndx < 0 || ndx > _number_elements)
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid atom index '" << token << "', max " << _number_elements << endl;
        return 0;
      }

      ndx--;    // convert to array offset
    }
    else if (2 == ntokens)   // atom type
    {
//    cerr << "Atomic symbol is '" << token << "'\n";
      zsymbol = token;
/*    if (zsymbol.length() <= 2)
        ;
      else if (atomic_symbols_can_have_arbitrary_length())
        ;
      else
      {
        cerr << "Molecule::_parse_v30_atom_record: atomic symbol too long '" << zsymbol << "'\n";
        return 0;
      }*/
    }
    else if (3 == ntokens)   // x coordinate
    {
      if (! token.numeric_value(x))
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid X coordinate\n";
        return 0;
      }
    }
    else if (4 == ntokens)   // y coordinate
    {
      if (! token.numeric_value(y))
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid Y coordinate\n";
        return 0;
      }
    }
    else if (5 == ntokens)   // z coordinate
    {
      if (! token.numeric_value(z))
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid Z coordinate\n";
        return 0;
      }
    }
    else if (token.starts_with("CHG="))
    {
      token.remove_leading_chars(4);
      if (! token.numeric_value(chg))
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid charge '" << token << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("MASS="))
    {
      token.remove_leading_chars(5);
      if (! token.numeric_value(mass) || mass <= 0)
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid mass value '" << token << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("RAD="))
    {
      token.remove_leading_chars(4);
      if (! token.numeric_value(rad) || rad < 0 || rad > 3)
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid radical value '" << token << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("CFG="))
    {
      token.remove_leading_chars(4);
      if (! token.numeric_value(cfg) || cfg < 0 || cfg > 3)
      {
        cerr << "Molecule::_parse_v30_atom_record: invalid cfg value '" << token << "'\n";
        return 0;
      }
    }
  }

// At a minimum we must have 'index type x y z'

  if (ntokens < 5)
  {
    cerr << "Molecule::_parse_v30_atom_record: only " << ntokens << " tokens on atom record\n";
    return 0;
  }

  if (NULL != _things[ndx])
  {
    cerr << "Molecule::_parse_v30_atom_record: atom " << ndx << " multiply defined\n";
    return 0;
  }

  if (convert_symbol_to_element)
  {
    _things[ndx] = mdlfos.create_mdl_atom(zsymbol, mass, chg, rad);
    if (NULL == _things[ndx])
    {
      cerr << "Molecule::_parse_v30_atom_record: cannot create atom\n";
      return 0;
    }
  }
  else
    _things[ndx] = mdlfos.create_mdl_atom("*", mass, chg, rad);   // just to get something there

  _things[ndx]->setxyz(x, y, z);

  if (0 == cfg)
    ;
  else if (! _mdl_atom_is_chiral_centre(ndx, cfg, mdlfos))
  {
    cerr << "Molecule::_parse_v30_atom_record: invalid chirality " << cfg << endl;
    cerr << buffer << endl;
    return 0;
  }

  return 1;
}

template <typename T>
int
Molecule::_read_v30_bond_list (T & input,
                               int nb,
                               int * aromatic_atom,
                               int * aromatic_bond)
{
  if (0 == nb)
    return 1;

  if (_bond_list.elements_allocated() < nb)
    _bond_list.resize(nb);

  IWString buffer;
  if (! read_next_v30_record(input, buffer))
  {
    cerr << "Molecule::_read_v30_bond_list: cannot read first record\n";
    return 0;
  }

  while ("M  V30 BEGIN BOND" != buffer)
  {
    if (! read_next_v30_record(input, buffer))
    {
      cerr << "Molecule::_read_v30_bond_list: premature EOF\n";
      return 0;
    }
  }

  for (int i = 0; i < nb; i++)
  {
    if (! read_next_v30_record(input, buffer))
    {
      cerr << "Molecule::_read_v30_bond_list: premature EOF\n";
      return 0;
    }

    if (! _parse_v30_bond_record(buffer, aromatic_atom, aromatic_bond[i]))
    {
      cerr << "Molecule::_read_v30_bond_list: invalid bond record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    if (AROMATIC_BOND == _bond_list.last_item()->btype())
      aromatic_bond[i] = 1;
  }

  if (! read_next_v30_record(input, buffer))
  {
    if ("M  V30 END BOND" != buffer)
    {
      cerr << "Molecule::_read_v30_bond_list: bond list improperly terminated\n";
      cerr << buffer << endl;
      return 0;
    }
  }

  check_bonding();

  return 1;
}

int
Molecule::_parse_v30_bond_record (const const_IWSubstring & buffer,
                                  int * aromatic_atom,
                                  int & aromatic_bond,
                                  int reading_query_file)
{
  int ntokens = 0;

  int btype, a1, a2;      // mandatory

  int cfg = 0;            // optional

  int i = 7;     // start tokenising after the "M  V30 ";
  const_IWSubstring token;
  while (buffer.nextword(token, i))
  {
    ntokens++;

    if (2 == ntokens)    // bond type
    {
      if (! token.numeric_value(btype) || btype < 1)
      {
        cerr << "Molecule::_parse_v30_bond_record: invalid bond type '" << token << "'\n";
        return 0;
      }
    }
    else if (3 == ntokens)  // A1
    {
      if (! token.numeric_value(a1) || a1 < 0 || a1 > _number_elements)
      {
        cerr << "Molecule::_parse_v30_bond_record: invalid A1 specifier '" << token << "'\n";
        return 0;
      }
      a1--;
    }
    else if (4 == ntokens)
    {
      if (! token.numeric_value(a2) || a2 < 0 || a2 > _number_elements)
      {
        cerr << "Molecule::_parse_v30_bond_record: invalid A2 specifier '" << token << "'\n";
        return 0;
      }
      a2--;
    }
    else if (token.starts_with("CFG="))
    {
      token.remove_leading_chars(4);
      if (! token.numeric_value(cfg) || cfg < 0 || cfg > 3)
      {
        cerr << "Molecule::_parse_v30_bond_record: invalid cfg value '" << token << "'\n";
        return 0;
      }
    }
  }

  if (ntokens < 4)
  {
    cerr << "Molecule::_parse_v30_bond_record: too few tokens on record " << ntokens << endl;
    return 0;
  }

// Bit of a tough decision about whether to put the query combinations in
// here but it means we don't have to duplicate the code

  if (4 == btype)
  {
    if (reading_query_file)
      btype = AROMATIC_BOND;
    else
    {
      if (! input_aromatic_structures())
      {
        cerr << "Molecule::_parse_v30_bond_record:aromatic input not enabled\n";
        return 0;
      }
  
      aromatic_atom[a1] = 1;
      aromatic_atom[a2] = 1;
      btype = SINGLE_BOND;
    }
  }
  else if (3 == btype)
    btype = TRIPLE_BOND;
  else if (2 == btype)
    btype = DOUBLE_BOND;
  else if (1 == btype)
    btype = SINGLE_BOND;
  else if (5 == btype)
    btype = SINGLE_BOND | DOUBLE_BOND;
  else if (6 == btype)
    btype = SINGLE_BOND | AROMATIC_BOND;
  else if (7 == btype)
    btype = DOUBLE_BOND | AROMATIC_BOND;
  else if (8 == btype)
    btype = SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND;
  else
  {
    cerr << "Molecule::_read_v30_bond_list:unrecognised bond type " << btype << endl;
    return 0;
  }

  if (! add_bond(a1, a2, btype, 1))     // 1 means partially built molecule
  {
    if (are_bonded(a1, a2) && static_cast<bond_type_t>(btype) == bond_between_atoms(a1, a2)->btype())
      cerr << "molecule::_parse_v30_bond_record:ignoring duplicate bond specification\n";
    else
    {
      cerr << "Molecule::_parse_v30_bond_record: cannot add bond between " << a1 << " and " << a2 << endl;
      return 0;
    }
  }

  if (cfg)     // for backward compatibility, we need to convert badk to V2 numbers
  {
    if (1 == cfg)
      ;
    else if (2 == cfg)
      cfg = 4;
    else if (3 == cfg)
      cfg = 6;

    _mdl_set_bond_directionality(a1, a2, cfg);
  }

  return 1;
}

/*
  Complicated by the fact that input lines are a maximum of 80 characters long
*/

int
write_v30_record (IWString & buffer,
                  std::ostream & output)
{
  if (buffer.length() <= 80)
  {
    output << buffer << '\n';

    return output.good();
  }

  int max_length = 80;
  while (buffer.length() > max_length)
  {
    buffer.insert("-\nM  V30 ", max_length - 1);
    max_length += 80;
  }

  output << buffer << '\n';

  return output.good();
}
template int Molecule::write_molecule_mdl_v30<std::ostream>(std::ostream&, IWString const&, int) const;

int
Molecule::_convert_sgroups_to_elements(const resizable_array_p<IWString> & sgroup)
{
  for (int i = 0; i < sgroup.number_elements(); ++i)
  {
    if (! _convert_sgroup_to_elements(*sgroup[i]))
    {
      cerr << "Molecule::_convert_sgroups_to_elements:cannot process SGROUP\n";
      cerr << *sgroup[i] << '\n';
      return 0;
    }
  }

  return 1;
}

static int
get_atoms_in_set(const IWString & buffer,
                 int & i,
                 const const_IWSubstring & atoms_token,      // ATOMS=(1
                 Set_of_Atoms & atoms,
                 const int atoms_in_molecule)
{
  assert (atoms_token.starts_with("ATOMS=("));

  const_IWSubstring token(atoms_token);
  token.remove_leading_chars(7);

  atom_number_t n;
  if (! token.numeric_value(n) || n < 1)
  {
    cerr << "get_atoms_in_set::invalid atoms in set '" << atoms_token << "'\n";
    return 0;
  }

  atoms.resize(n);

  int got_last_token = 0;

  while (buffer.nextword(token, i))
  {
    if (token.ends_with(')'))
    {
      got_last_token = 1;
      token.chop();
    }

    atom_number_t x;
    if (! token.numeric_value(x) || x < 1 || x > atoms_in_molecule)
    {
      cerr << "get_atoms_in_set:invalid atom number '" << token << "'\n";
      return 0;
    }

    atoms.add(x-1);

    if (got_last_token)
      break;
  }

  if (! got_last_token)
    cerr << "get_atoms_in_set::end of group not found\n";

  if (n != atoms.number_elements())
    cerr << "get_atoms_in_set:count mismatch, expected " << n << " atoms, got " << atoms.number_elements() << endl;

//cerr << " atoms " << atoms << endl;

  return atoms.number_elements();
}

/*
  Will fail if there are multiple spaces in the label. For example LABEL="A   B"
*/

static int
fetch_rest_of_quoted_label(const IWString & sgroup,
                    int & i,
                    IWString & label)
{
  assert (label.starts_with('"'));

  label.remove_leading_chars(1);

  const_IWSubstring token;

  while (sgroup.nextword(token, i))
  {
    if (! token.ends_with('"'))
    {
      label.append_with_spacer(token);
      continue;
    }

    token.chop(1);
    label.append_with_spacer(token);
    return 1;
  }

  cerr << "fetch_rest_of_quoted_label:did not get end of quoted label\n";
  return 0;
}

int
Molecule::_convert_sgroup_to_elements(const IWString & sgroup)
{
  const_IWSubstring token;

  int i = 0;
  if (! sgroup.nextword(token, i) || ! sgroup.nextword(token, i))
    return 0;

  if ("SUP" != token)
    return 1;

  if (! sgroup.nextword(token, i))    // skip over external index
    return 0;

  IWString label;
  Set_of_Atoms atoms;

  while (sgroup.nextword(token, i))
  {
//  cerr << "Examining token '" << token << "'\n";

    if (token.starts_with("ATOMS="))
    {
      if (! get_atoms_in_set(sgroup, i, token, atoms, _number_elements))
      {
        cerr << "Molecule::_convert_sgroup_to_elements:cannot fetch atoms\n";
        return 0;
      }
    }
    else if (token.starts_with("LABEL="))
    {
      token.remove_leading_chars(6);
      label = token;
      if (label.starts_with('"'))
      {
        fetch_rest_of_quoted_label(sgroup, i, label);
      }
    }
  }

  if (0 == label.length())     // not all SGROUP blocks have a label. Might be just for display
  {
//  cerr << "Molecule::_convert_sgroup_to_elements:no LABEL\n";
//  cerr << sgroup << endl;
//  return 0;
    return 1;
  }

  if (0 == atoms.number_elements())
  {
    cerr << "Molecule::_convert_sgroup_to_elements:no ATOMS\n";
    cerr << sgroup << endl;
    return 0;
  }

  if (1 != atoms.number_elements())    // cannot handle these
    return 1;

  if (isdigit(label[0]))
    return 1;

  const Element * e = get_element_from_symbol_no_case_conversion(label);
  if (NULL == e)
  {
    if (! auto_create_new_elements())
    {
      cerr << "Molecule::_convert_sgroup_to_elements:cannot create element for '" << label << "'\n";
      return 0;
    }

    e = create_element_with_symbol(label);
    if (NULL == e)
    {
      cerr << "Molecule::_convert_sgroup_to_elements:cannot create element '" << label << "'\n";
      return 0;
    }
  }

  _things[atoms[0]]->set_element(e);

  return 1;
}
