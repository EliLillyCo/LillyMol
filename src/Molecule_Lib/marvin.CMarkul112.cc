#include <stdlib.h>

#include "Markup.h"

#include "iwstring_data_source.h"

#define COMPILING_MARVIN_CC

#include "molecule.h"

#include "marvin.h"

static const Marvin_Structure_Information * msi = nullptr;

void 
set_marvin_structure_information_for_writing (const Marvin_Structure_Information * s)
{
  msi = s;
}

void
Marvin_Structure_Information::reset ()
{
  _atom_colour.resize_keep_storage(0);
  _bond_colour.resize_keep_storage(0);

  reset_atoms_and_bonds();

  return;
}

void
Marvin_Structure_Information::reset_atoms_and_bonds ()
{
  _atom_number_colour_index.resize(0);
  _bond_number_colour_index.resize(0);

  return;
}

static int
do_add_looking_for_N_and_hash (const IWString & s,
                               resizable_array_p<IWString> & z)
{
  if ('N' == s)
  {
    z.add(new IWString('N'));
    return 1;
  }

  if (s.starts_with('#'))
  {
    IWString * t = new IWString (s);
    z.add(t);
  }

  IWString * t = new IWString;
  t->add('#');
  (*t) << s;
  
  z.add(t);

  return 1;
}

int
Marvin_Structure_Information::add_atom_colour (const const_IWSubstring & s)
{
  return do_add_looking_for_N_and_hash(s, _atom_colour);
}

int
Marvin_Structure_Information::add_bond_colour (const const_IWSubstring & s)
{
  return do_add_looking_for_N_and_hash(s, _bond_colour);
}

int
Marvin_Structure_Information::colour_index_for_atom (atom_number_t a) const
{
  hash_map<atom_number_t, unsigned int>::const_iterator f = _atom_number_colour_index.find(a);

  if (f == _atom_number_colour_index.end())
    return 0;

  return (*f).second;
}

int
Marvin_Structure_Information::colour_index_for_bond (int bond_number) const
{
  hash_map<atom_number_t, unsigned int>::const_iterator f = _bond_number_colour_index.find(bond_number);

  if (f == _atom_number_colour_index.end())
    return 0;

  return (*f).second;
}

int
Marvin_Structure_Information::write_atom_and_bond_colours (ostream & os) const
{
  if (_atom_colour.number_elements())
  {
    os << "        atomSetRGB=\"";
    for (int i = 0; i < _atom_colour.number_elements(); i++)
    {
      if (i > 0)
        os << ',';

      os << i << ":" << *(_atom_colour[i]);
    }
    os << "\"\n";
  }

  if (_bond_colour.number_elements())
  {
    os << "          bondSetRGB=\"";
    for (int i = 0; i < _bond_colour.number_elements(); i++)
    {
      if (i > 0)
        os << ',';

      os << i << ":" << *(_bond_colour[i]);
    }
    os << "\"\n";
  }

  return 1;
}

int
Molecule::write_molecule_mrv (ostream & os)
{
  os << "<MDocument";

  if (nullptr != msi)
    msi->write_atom_and_bond_colours (os);

  os << ">\n";

  int rc = _write_molecule_mrv(os);

  os << "</MDocument>\n";

  msi = nullptr;

  return rc;
}

int
Molecule::_write_molecule_mrv (ostream & os) const
{
  os << "  <MChemicalStruct>\n";
  os << "    <molecule title=\"" << _molecule_name << "\">\n";
  os << "      <atomArray\n";
  _write_atoms_mrv(os);
  os << "      />\n";
  os << "      <bondArray>\n";
  _write_bonds_mrv(os);
  os << "      </bondArray>\n";

  os << "    </molecule>\n";
  os << "  </MChemicalStruct>\n";

  return os.good();
}

int
Molecule::_write_atoms_mrv (ostream & os) const
{
  os << "          atomID=\"";
  for (int i = 0; i < _number_elements; i++)
  {
    if (i > 0)
      os << ' ';
    os << 'a' << (i + 1);
  }
  os << "\"\n";
  os << "          elementType=\"";
  for (int i = 0; i < _number_elements; i++)
  {
    if (i > 0)
      os << ' ';
    os << _things[i]->atomic_symbol();
  }
  os << "\"\n";

  if (maximum_isotope() > 0)
  {
    os << "          isotope=\"";
    for (int i = 0; i < _number_elements; i++)
    {
      if (i > 0)
        os << ' ';
      os << _things[i]->isotope();
    }
    os << "\"\n";
  }

  if (has_formal_charges())
  {
    os << "          formalCharge=\"";
    for (int i = 0; i < _number_elements; i++)
    {
      if (i > 0)
        os << ' ';
      os << _things[i]->formal_charge();
    }
    os << "\"\n";
  }

  if (nullptr == msi)   // no atom colours present
    ;
  else if (msi->atom_colour_specifications_present())
  {
    os << "          mrvSetSeq=\"";
    for (int i = 0; i < _number_elements; i++)
    {
      if (i > 0)
        os << ' ';
      if (nullptr == msi)
        os << '0';
      else
        os << msi->colour_index_for_atom(i);
    }
    os << "\"\n";
  }

  os << "          x2=\"";
  for (int i = 0; i < _number_elements; i++)
  {
    if (i > 0)
      os << ' ';
    os << _things[i]->x();
  }
  os << "\"\n";
  os << "          y2=\"";
  for (int i = 0; i < _number_elements; i++)
  {
    if (i > 0)
      os << ' ';
    os << _things[i]->y();
  }
  os << "\"\n";

  return 1;
}

int
Molecule::_write_bonds_mrv (ostream & os) const
{
  int n = _bond_list.number_elements();

  for (int i = 0; i < n; i++)
  {
    const Bond * b = _bond_list[i];
    os << "        <bond atomRefs2=\"a" << (b->a1() + 1) << " a" << (b->a2()+1) << "\" order=\"";
    if (b->is_single_bond())
      os << '1';
    else if (b->is_double_bond())
      os << '2';
    else if (b->is_triple_bond())
      os << '3';
    else if (b->is_aromatic())
      os << 'A';
    else
      os << '5';         // HUH!!

    os << "\" ";
    if (nullptr == msi)
    {
      os << "/>\n";
      continue;
    }

    int c = msi->colour_index_for_bond (i);
    if (c > 0)
      os << "mrvSetSeq=\"" << c << "\" ";
    os << " />\n";
  }

  return 1;
}

int
Molecule::read_molecule_mrv_ds (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  int found_mdocument = 0;
  int lines_read = 0;
  while (input.next_record(buffer))
  {
    lines_read++;

    if (! buffer.starts_with("<MDocument"))
      continue;

    found_mdocument = 1;
    break;
  }

  if (found_mdocument)    // great
    ;
  else if (0 == lines_read)   // normal eof
    return 0;
  else   // bad
  {
    cerr << "Molecule::read_molecule_mrv_ds:no MDocument record found\n";
    return 0;
  }

  IWString xml;

  xml << buffer;

  found_mdocument = 0;    // the end this time
  while (input.next_record(buffer))
  {
    xml << buffer;

//  cerr << "Check '" << buffer << "'\n";
    if (! buffer.starts_with("</MDocument"))
      continue;

    found_mdocument = 1;
    break;
  }

  if (! found_mdocument)
  {
    cerr << "Molecule::read_molecule_mrv_ds:no termination MDocument\n";
    return 0;
  }

  CMarkup cml;

  if (! cml.SetDoc(xml.null_terminated_chars()))
  {
    cerr << "Molecule::read_molecule_mrv_ds:cannot initialise xml document\n";
    return 0;
  }

  if (! cml.IsWellFormed ())
  {
    cerr << "Molecule::read_molecule_mrv_ds:malformed xml\n";
    return 0;
  }

//cerr << "Successfully initialised xml\n";

  if (! cml.FindElem("MDocument"))   // of course it is there!
  {
    cerr << "Molecule::read_molecule_mrv_ds:cannot seek to start of xml\n";
    return 0;
  }

  if (! cml.IntoElem())
  {
    cerr << "Molecule::read_molecule_mrv_ds:cannot descend into MDocument\n";
    return 0;
  }

  while (cml.FindElem())
  {
    std::string t = cml.GetTagName();

    if ("MChemicalStruct" != t)
      continue;

    cml.IntoElem();

    return read_molecule_mrv(cml);
  }

  cerr << "Molecule::read_molecule_mrv:no MChemicalStruct attribute\n";
  return 0;
}

int
Molecule::read_molecule_mrv (CMarkup & cml)
{
  while (cml.FindElem())
  {
    std::string t = cml.GetTagName();

    if ("molecule" != t)
      continue;

    return _read_molecule_mrv(cml);
  }

  cerr << "molecule::read_molecule_mrv:no 'molecule' tag\n";
  return 0;
}

int
Molecule::_read_molecule_mrv (CMarkup & cml)
{
  for (int i = 0; ; i++)
  {
    std::string n = cml.GetAttribName(i);
    if (0 == n.size())
      break;

    if ("title" != n)
      continue;

    std::string v = cml.GetAttrib(n);

    _molecule_name = v;
//  cerr << "Name " << i << " is '" << n << "', value '" << v << "'\n";
  }

  cml.IntoElem();    // into molecule object

  while (cml.FindElem())
  {
    std::string t = cml.GetTagName();

    if ("atomArray" == t)
    {
      if (! _read_atom_array_mrv(cml))
      {
        cerr << "Molecule::_read_molecule_mrv:cannot read atom array\n";
        return 0;
      }
    }
    else if ("bondArray" == t)
    {
      cml.IntoElem();
      if (! _read_bond_array_mrv(cml))
      {
        cerr << "Molecule::_read_molecule_mrv:cannot read bond array\n";
        return 0;
      }
      cml.OutOfElem();
    }
    else
      cerr << "Molecule::_read_molecule_mrv:unrecognised tag '" << t << "'\n";
  }

  cml.OutOfElem();

  return 1;
}

static int
zero_length_or_this_many_tokens (const IWString & s,
                                 int nw)
{
  if (0 == s.length())
    return 1;

  if (nw == s.nwords())
    return 1;

  cerr << "Molecule::_read_atom_array_mrv:inconsistent atom list lengths\n";
  cerr << "'" << s << "'\n";

  return 0;
}

/*
  Could not get a template for this to work, so we have two versions...
*/

template <typename T>
int
convert_to_numeric (const IWString & s,
                    resizable_array<T> & v,
                    T minval)
{
  if (0 == s.length())
    return 1;

  int i = 0;
  const_IWSubstring token;

  while (s.nextword(token, i))
  {
    if ('0' == token)    //  a common case
    {
      v.add(0);
      continue;
    }

    T t;
    if (! token.numeric_value(t))
    {
      cerr << "Invalid numeric '" << t << "'\n";
      return 0;
    }

    if (t < minval)
    {
      cerr << "Value out of range " << t << " compare " << minval << endl;
      return 0;
    }

    v.add(t);
  }
  
  return v.number_elements();
}

int
Molecule::_read_atom_array_mrv (CMarkup & cml)
{
  IWString elementType, x2,y2,x3,y3,z3,isotope,formalCharge;

  for (int i = 0; ; i++)
  {
    std::string d = cml.GetAttribName(i);

    if (0 == d.size())
      break;

    std::string v = cml.GetAttrib(d);

//  cerr << "Atom data '" << d << "' is '" << v << "'\n";

    if ("elementType" == d)
      elementType = v;
    else if ("isotope" == d)
      isotope = v;
    else if ("formalCharge" == d)
      formalCharge = v;
    else if ("x2" == d)
      x2 = v;
    else if ("y2" == d)
      y2 = v;
    else if ("x3" == d)
      x3 = v;
    else if ("y3" == d)
      y3 = v;
    else if ("z3" == d)
      z3 = v;
  }

  if (0 == elementType.length())
  {
    cerr << "Molecule::_read_atom_array_mrv:no elements specified\n";
    return 0;
  }

  int dimensionality = 0;

  if (x2.length() > 0 && y2.length() > 0 && 
      0 == x3.length() && 0 == y3.length() && 0 == z3.length())
    dimensionality = 2;
  else if (0 == x2.length() && 0 == y2.length() &&
           x3.length() && y3.length() && z3.length())
    dimensionality = 3;
  else if (0 == x2.length() && 0 == y2.length() &&
           0 == x3.length() && 0 == y3.length() && 0 == z3.length())
    ;      // zero dimensionality
  else
  {
    cerr << "Molecule::_read_atom_array_mrv:inconsistent geometry specifications\n";
    return 0;
  }

  int natoms = elementType.nwords();
  cerr << "Molecule contains " << natoms << " atoms\n";
  if (! zero_length_or_this_many_tokens(x2, natoms))
    return 0;
  if (! zero_length_or_this_many_tokens(y2, natoms))
    return 0;
  if (! zero_length_or_this_many_tokens(x3, natoms))
    return 0;
  if (! zero_length_or_this_many_tokens(y3, natoms))
    return 0;
  if (! zero_length_or_this_many_tokens(z3, natoms))
    return 0;
  if (! zero_length_or_this_many_tokens(isotope, natoms))
    return 0;
  if (! zero_length_or_this_many_tokens(formalCharge, natoms))
    return 0;

  int elementType_i = 0;
  int x2_i = 0;
  int x3_i = 0;
  int y2_i = 0;
  int y3_i = 0;
  int z3_i = 0;

  int notused = 0;    // no isotopes in atomic symbols here

  resizable_array<int> isotopes;
  if (! convert_to_numeric(isotope, isotopes, static_cast<int>(0)))
  {
    cerr << "Molecule::_read_atom_array_mrv:invalid isotopic specification\n";
    return 0;
  }

  resizable_array<formal_charge_t> formal_charges;
  if (! convert_to_numeric(formalCharge, formal_charges, -9999))
  {
    cerr << "Molecule::_read_atom_array_mrv:invalid formal charge\n";
    return 0;
  }

  resizable_array<coord_t> x, y, z;
  if (2 == dimensionality)
  {
    if (! convert_to_numeric(x2, x, - numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid x coordinate\n";
      return 0;
    }

    if (! convert_to_numeric(y2, y, - numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid y coordinate\n";
      return 0;
    }
  }
  else if (3 == dimensionality)
  {
    if (! convert_to_numeric(x3, x, - numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid x coordinate\n";
      return 0;
    }

    if (! convert_to_numeric(y3, y, - numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid y coordinate\n";
      return 0;
    }

    if (! convert_to_numeric(z3, z, - numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid z coordinate\n";
      return 0;
    }
  }

  const_IWSubstring elementType_token, x2_token, y2_token, x3_token, y3_token, z3_token, isotope_token, formalCharge_token;
  const_IWSubstring token;

  int ndx = 0;
  while (elementType.nextword(elementType_token, elementType_i))
  {
    const Element * e = get_element_from_symbol(elementType_token, notused);

    if (nullptr == e)
    {
      cerr << "Molecule::_read_atom_array_mrv:cannot create element from '" << elementType_token << "'\n";
      return 0;
    }
    
    Atom * a = new Atom(e);

    if (formal_charges.empty())
      ;
    else if (0 != formal_charges[ndx])
      a->set_formal_charge(formal_charges[ndx]);

    if (isotopes.empty())
      ;
    else if (0 != isotopes[ndx])
      a->set_isotope(isotopes[ndx]);

    if (2 == dimensionality)
    {
      a->setxyz(x[ndx], y[ndx], static_cast<coord_t>(0.0));
    }
    else if (3 == dimensionality)
    {
      a->setxyz(x[ndx], y[ndx], z[ndx]);
    }
      
    add(a);
    ndx++;
  }

  return 1;
}

int
Molecule::_read_bond_array_mrv (CMarkup & cml)
{
  while (cml.FindElem())
  {
    std::string t = cml.GetTagName();

//  cerr << "Bond tag is '" << t << "'\n";
    if ("bond" != t)
      continue;

    IWString atomRefs2, order;

    for (int i = 0; ; i++)
    {
      std::string n = cml.GetAttribName(i);

      if (0 == n.length())
        break;

      std::string d = cml.GetAttrib(n);

//    cerr << "Bond data for '" << n << " is '" << d << "'\n";

      if ("atomRefs2" == n)
        atomRefs2 = d;
      else if ("order" == n)
        order = d;
    }

    if (0 == atomRefs2.length() || 0 == order.length())
    {
      cerr << "Molecule::_read_bond_array_mrv:invalid bond\n";
      return 0;
    }

    if (2 != atomRefs2.nwords())
    {
      cerr << "Molecule::_read_bond_array_mrv:bond must have two atoms\n";
      return 0;
    }

    const_IWSubstring sa1, sa2;
    atomRefs2.split(sa1, ' ', sa2);

    if (! sa1.starts_with('a') || ! sa2.starts_with('a'))
    {
      cerr << "Molecule::_read_bond_array_mrv:bond specifications must start with 'a' '" << atomRefs2 << "'\n";
      return 0;
    }

    sa1++;
    sa2++;

    atom_number_t a1, a2;

    if (! sa1.numeric_value(a1) || a1 < 1 || a1 > _number_elements)
    {
      cerr << "Molecule::_read_bond_array_mrv:invalid atom in bond '" << atomRefs2 << "'\n";
      return 0;
    }
    if (! sa2.numeric_value(a2) || a2 < 1 || a2 > _number_elements || a1 == a2)
    {
      cerr << "Molecule::_read_bond_array_mrv:invalid atom in bond '" << atomRefs2 << "'\n";
      return 0;
    }

    a1--;
    a2--;

    bond_type_t bt;

    if ('1' == order)
      bt = SINGLE_BOND;
    else if ('2' == order)
      bt = DOUBLE_BOND;
    else if ('3' == order)
      bt = TRIPLE_BOND;
    else if ('A' == order)
    {
      cerr << "Molecule::_read_bond_array_mrv:aromatic input not build, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
      return 0;
    }

    add_bond(a1, a2, bt, 1);
  }

  return 1;
}
