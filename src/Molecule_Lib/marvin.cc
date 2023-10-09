#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/xmlParser/xmlParser.h"

#define COMPILING_MARVIN_CC

#include "molecule.h"
#include "marvin.h"

using std::cerr;
using std::endl;

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
  _atom_number_colour_index.clear();
  _bond_number_colour_index.clear();

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
  IW_Hash_Map<atom_number_t, unsigned int>::const_iterator f = _atom_number_colour_index.find(a);

  if (f == _atom_number_colour_index.end())
    return 0;

  return (*f).second;
}

int
Marvin_Structure_Information::colour_index_for_bond (int bond_number) const
{
  IW_Hash_Map<atom_number_t, unsigned int>::const_iterator f = _bond_number_colour_index.find(bond_number);

  if (f == _atom_number_colour_index.end())
    return 0;

  return (*f).second;
}

int
Marvin_Structure_Information::write_atom_and_bond_colours (std::ostream & os) const
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
Molecule::write_molecule_mrv (std::ostream & os)
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
Molecule::_write_molecule_mrv (std::ostream & os) const
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
Molecule::_write_atoms_mrv (std::ostream & os) const
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
Molecule::_write_bonds_mrv (std::ostream & os) const
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
    if (nullptr != msi)
    {
      int c = msi->colour_index_for_bond (i);
      if (c > 0)
        os << "mrvSetSeq=\"" << c << "\" ";
    }

    if (! b->is_wedge_definitive())
    {
      os << " />\n";
      continue;
    }

    os << ">\n";
    if (b->is_wedge_up())
      os << "          <bondStereo>H</bondStereo>\n";
    else if (b->is_wedge_down())
      os << "          <bondStereo>W</bondStereo>\n";

    os << "        </bond>\n";
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

    if ("</cml>" == buffer)
    {
      input.next_record(buffer);   // force eof
      return 0;
    }

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

  XMLResults xe;
  XMLNode xMainNode = XMLNode::parseString(xml.null_terminated_chars(), "MDocument", &xe);

  if (0 != xe.error)
  {
    cerr << "Molecule::read_molecule_mrv_ds:error in xml data, line " << xe.nLine << endl;
    return 0;
  }

#ifdef DEBUG_MARVIN_READ_MOLECULE
  cerr << "Successfully initialised xml, name '" << xMainNode.getName() << "'\n";
#endif

  XMLNode xNode=xMainNode.getChildNode("MChemicalStruct");

  if (nullptr == xNode.getName())
  {
    cerr << "molecule::read_molecule_mrv_ds:could not get 'MChemicalStruct' attribute\n";
    return 0;
  }

#ifdef DEBUG_MARVIN_READ_MOLECULE
  cerr << "Got '" << xNode.getName() << "'\n";
#endif

  return read_molecule_mrv_mchemical(xNode);
}

int
Molecule::read_molecule_mrv_mchemical (XMLNode & mchemicalstruct)
{
  XMLNode molecule = mchemicalstruct.getChildNode("molecule");

  if (nullptr == molecule.getName())
  {
    cerr << "molecule::read_molecule_mrv:no 'molecule' tag in xml\n";
    return 0;
  }

  return read_molecule_mrv_molecule(molecule);
}

int
Molecule::read_molecule_mrv_molecule (XMLNode & xml)
{
  XMLCSTR mname = xml.getAttribute("title", 0);

  if (nullptr != mname)
    _molecule_name = mname;

//cerr << "Molecule name set to '" << _molecule_name << "'\n";

  XMLNode atomArray = xml.getChildNode("atomArray");

  if (nullptr == atomArray.getName())
  {
    cerr << "Molecule::_read_molecule_mrv:no atomArray in xml\n";
    return 0;
  }

  if (! _read_atom_array_mrv(atomArray))
  {
    cerr << "Molecule::_read_molecule_mrv:cannot read atom array\n";
    return 0;
  }

  if (0 == _number_elements)   // cannot be anything else, strange...
  {
    cerr << "Molecule::_read_bond_array_mrv:no atoms\n";
    return 1;
  }

// There must be at least natoms bonds present, safe upper estimate

  int * aromatic_bonds = new_int(_number_elements + _number_elements); std::unique_ptr<int[]> free_aromatic_bonds(aromatic_bonds);

  XMLNode bondarray = xml.getChildNode("bondArray");
  if (nullptr == bondarray.getName())
    ;
  else if (! _read_bond_array_mrv(bondarray, aromatic_bonds))
  {
    cerr << "Molecule::_read_molecule_mrv:cannot read bond array\n";
    return 0;
  }

  if (locate_item_in_array(1, _bond_list.number_elements(), aromatic_bonds) < 0)   // no aromatic bonds
    ;
  else if (! _final_processing_of_aromatic_mdl_input(nullptr, aromatic_bonds))
  {
    cerr << "Molecule::_read_molecule_mrv:cannot resolve aromatic bonds\n";
    return 0;
  }

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

/*
  Atoms can appear as individual items or as an array
*/

int
Molecule::_read_atom_array_mrv(XMLNode & xml)
{
  XMLCSTR xml_atomID = xml.getAttribute("atomID", 0);

  if (nullptr == xml_atomID)
    return _read_atom_array_mrv_individual_attributes(xml);

  XMLCSTR xml_elementType = xml.getAttribute("elementType");
  if (nullptr == xml_elementType)
  {
    cerr << "Molecule::_read_atom_array_mrv:no elementType\n";
    return 0;
  }

  IWString elementType = xml_elementType;

  IWString x2 = xml.getAttribute("x2", 0);
  IWString y2 = xml.getAttribute("y2", 0);

  IWString x3 = xml.getAttribute("x3", 0);
  IWString y3 = xml.getAttribute("y3", 0);
  IWString z3 = xml.getAttribute("z3", 0);

  IWString isotope = xml.getAttribute("isotope");
  IWString formalCharge = xml.getAttribute("formalCharge");
  IWString mrvMap = xml.getAttribute("mrvMap");

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
//cerr << "Molecule contains " << natoms << " atoms\n";
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
//int x2_i = 0;
//int x3_i = 0;
//int y2_i = 0;
//int y3_i = 0;
//int z3_i = 0;

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
    if (! convert_to_numeric(x2, x, - std::numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid x coordinate\n";
      return 0;
    }

    if (! convert_to_numeric(y2, y, - std::numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid y coordinate\n";
      return 0;
    }
  }
  else if (3 == dimensionality)
  {
    if (! convert_to_numeric(x3, x, -std::numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid x coordinate\n";
      return 0;
    }

    if (! convert_to_numeric(y3, y, -std::numeric_limits<coord_t>::max()))
    {
      cerr << "Molecule::_read_atom_array_mrv:invalid y coordinate\n";
      return 0;
    }

    if (! convert_to_numeric(z3, z, -std::numeric_limits<coord_t>::max()))
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
    const Element * e = get_element_from_symbol_no_case_conversion(elementType_token);

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

template <typename T>
int
convert_to_numeric (XMLCSTR & s,
                    T & v)
{
  IWString tmp(s);

  return tmp.numeric_value(v);
}

int
Molecule::_read_atom_array_mrv_individual_attributes(const XMLNode & xml)
{
  for (int i = 0; ; i++)
  {
    XMLNode xml_a = xml.getChildNode(i);
    if (nullptr == xml_a.getName())
      break;

//  cerr << "Processing atom node '" << xml_a.getName() << "'\n";

    if (0 != strncmp("atom", xml_a.getName(), 4))
      continue;

    XMLCSTR elementType = xml_a.getAttribute("elementType");

//  cerr << "Element is '" << elementType << "'\n";

    if (nullptr == elementType)
    {
      cerr << "Molecule::_read_atom_array_mrv_individual_attributes:no elementType attribute\n";
      return 0;
    }

    XMLCSTR x2 = xml_a.getAttribute("x2");
    XMLCSTR y2 = xml_a.getAttribute("y2");
    XMLCSTR x3 = xml_a.getAttribute("x3");
    XMLCSTR y3 = xml_a.getAttribute("y3");
    XMLCSTR z3 = xml_a.getAttribute("z3");

    XMLCSTR isotope = xml_a.getAttribute("isotope");
    XMLCSTR formalCharge = xml_a.getAttribute("formalCharge");

    const Element * e = get_element_from_symbol_no_case_conversion(elementType);

    if (nullptr != e)
      ;
    else if (! auto_create_new_elements())
    {
      cerr << "Molecule::_read_atom_array_mrv_individual_attributes:cannot create element from '" << elementType << "'\n";
      return 0;
    }
    else
    {
      e = create_element_with_symbol(elementType);
      if (nullptr == e)
      {
        cerr << "Molecule::read_molecule_mrv_molecule:cannot create element '" << elementType << "'\n";
        return 0;
      }
    }

    Atom * a = new Atom(e);

    if (nullptr != formalCharge)
    {
      formal_charge_t q;
      if (! convert_to_numeric(formalCharge, q) || ! reasonable_formal_charge_value(q))
      {
        cerr << "Molecule::_read_atom_array_mrv_individual_attributes:invalid charge '" << formalCharge << "'\n";
        return 0;
      }

      a->set_formal_charge(q);
    }

    if (nullptr != isotope)
    {
      int iso;
      if (! convert_to_numeric(isotope, iso) || iso < 0)
      {
        cerr << "Molecule::_read_atom_array_mrv_individual_attributes:invalid isotopic specification '" << isotope << "'\n";
        return 0;
      }

      a->set_isotope(iso);
    }

    if (nullptr != x2)    // 2D coordinates specified
    {
      coord_t x;
      if (! convert_to_numeric(x2, x))
      {
        cerr << "Molecule::_read_atom_array_mrv_individual_attributes:invalid x2 '" << x2 << "'\n";
        return 0;
      }

      coord_t y;

      if (nullptr == y2)              // silently ignore this!!?
        y = static_cast<coord_t>(0.0);
      else if (! convert_to_numeric(y2, y))
      {
        cerr << "Molecule::_read_atom_array_mrv_individual_attributes:invalid y2 '" << y2 << "'\n";
        return 0;
      }

      a->setxyz(x, y, static_cast<coord_t>(0.0));
    }
    else if (nullptr != x3)
    {
      if (nullptr == y3 || nullptr == z3)
      {
        cerr << "Molecule::_read_atom_array_mrv_individual_attributes:incomplete 3d specification\n";
        return 0;
      }
      coord_t x, y, z;

      if (! convert_to_numeric(x3, x) ||
          ! convert_to_numeric(y3, y) ||
          ! convert_to_numeric(z3, z))
      {
        cerr << "Molecule::_read_atom_array_mrv_individual_attributes:invalid 3d specification '" << x3 << "', '" << y3 << "', '" << z3 << "'\n";
        return 0;
      }

      a->setxyz(x, y, z);
    }

    add(a);
  }

  return 1;
}

int
Molecule::_read_bond_array_mrv (XMLNode & xml,
                                int * aromatic_bonds)
{
  int wedge_bonds_encountered = 0;

  for (int i = 0; ; i++)
  {
    XMLNode xml_b = xml.getChildNode(i);
    if (nullptr == xml_b.getName())
      break;

#ifdef DEBUG_MARVIN_BOND_LIST
    cerr << "Examining child '" << xml_b.getName() << "'\n";
#endif

    if (0 != strncmp("bond", xml_b.getName(), 4))
      continue;

    IWString atomRefs2 = xml_b.getAttribute("atomRefs2", 0);
    IWString order = xml_b.getAttribute("order", 0);

    if (0 == atomRefs2.length() || 0 == order.length())
    {
      cerr << "Molecule::_read_bond_array_mrv:invalid bond\n";
      return 0;
    }

#ifdef DEBUG_MARVIN_BOND_LIST
    cerr << "atomRefs2 '" << atomRefs2 << "', order " << order << "'\n";
#endif

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
      int nb = _bond_list.number_elements();
      aromatic_bonds[nb] = 1;
      bt = SINGLE_BOND;
    }
    else
    {
      cerr << "Molecule::_read_bond_array_mrv:what kind of bond is '" << order << "', set to single\n";
      bt = SINGLE_BOND;
    }

    add_bond(a1, a2, bt, 1);

    for (int j = 0; ; j++)
    {
      XMLNode xml_bs = xml_b.getChildNode(j);
      if (nullptr == xml_bs.getName())
        break;

//    cerr << "Processing '" << xml_bs.getName() << "'\n";
      if (0 != strncmp("bondStereo", xml_bs.getName(), 10))
        continue;

      Bond * b = _bond_list.last_item();

      XMLCSTR att = xml_bs.getText();
      if (0 == strncmp(att, "H", 1))
      {
        b->set_wedge_up();
        wedge_bonds_encountered++;
      }
      else if (0 == strncmp(att, "W", 1))
      {
        b->set_wedge_down();
        wedge_bonds_encountered++;
      }
      else
      {
        cerr << "Molecule::_read_bond_array_mrv:unrecognised wedge bond attribute '" << att << "', ignored\n";
      }
    }
  }

  if (wedge_bonds_encountered)
    discern_chirality_from_wedge_bonds();

  return 1;
}
