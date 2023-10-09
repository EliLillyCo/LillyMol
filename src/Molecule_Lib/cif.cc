/*
  Read .cif files
*/

#include <algorithm>
#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

class CIF_Loop
{
  private:
    IWString _stem;

    resizable_array_p<IWString> _title;
    resizable_array_p<IWString> _zdata;

  public:
    CIF_Loop ();
    ~CIF_Loop ();

    const IWString & stem () const { return _stem;}

    int build (iwstring_data_source &);

    unsigned int dataitems () const { return _zdata.size();}
    unsigned int number_title_records () const { return _title.size();}

    int string_value (const char * key, IWString & v, int ndx) const;

    const IWString & title       (int ndx) const { return *(_title[ndx]);}
    const IWString & data_record (int ndx) const { return *(_zdata[ndx]);}
};

#define COMPILING_MOLECULE_CIF

#include "molecule.h"

CIF_Loop::CIF_Loop ()
{
  return;
}

CIF_Loop::~CIF_Loop ()
{
  return;
}

int
CIF_Loop::build (iwstring_data_source & input)
{
  _stem.resize_keep_storage(0);
  _title.resize_keep_storage(0);
  _zdata.resize_keep_storage(0);

  if (! input.skip_past("loop_"))
    return 0;

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
#ifdef DEBUG_CIF_LOOP_BUILD
    cerr << "CIF_Loop::build:examining '" << buffer << "' at line " << input.lines_read() << endl;
#endif

    if (! buffer.starts_with("_chem_comp"))
    {
      input.push_record();
      break;
    }

    if (0 == _stem.length())
    {
      _stem = buffer;
      _stem.truncate_at_first('.');
      _stem.remove_leading_chars(10);
      if (_stem.starts_with('_'))
        _stem.remove_leading_chars(1);
    }

    buffer.remove_up_to_first('.');
    _title.add(new IWString(buffer));
  }

#ifdef DEBUG_CIF_LOOP_BUILD
  cerr << "Reading data\n";
#endif

// now read the associated data records

  while (input.next_record(buffer))
  {
    if (buffer.starts_with("loop_") || buffer.starts_with("data_comp"))
    {
      input.push_record();
      break;
    }

//  cerr << "Adding '" << buffer << "' to zdata\n";
    _zdata.add(new IWString(buffer));
  }

  if (_zdata.empty())
  {
    cerr << "CIF_Loop::build:no data\n";
    return 0;
  }

#ifdef DEBUG_CIF_LOOP_BUILD
  cerr << "Found " << _zdata.number_elements() << " dataitems\n";
#endif

  return _zdata.number_elements();
}

int
CIF_Loop::string_value (const char * key,
                        IWString & v,
                        const int ndx) const
{
  const_IWSubstring mykey(key);

  if (! _zdata.ok_index(ndx))
  {
    cerr << "CIF_Loop::string_value:request for index " << ndx << " but only " << _zdata.number_elements() << " dataitems available\n";
    return 0;
  }

  const IWString & buffer = *(_zdata[ndx]);

#ifdef DEBUG_STRING_VALUE
  cerr << "Looking for '" << mykey << "' in";
  for (auto f = _title.cbegin(); f != _title.cend(); ++f)
  {
    cerr << ' ' << *(f);
  }
  cerr << endl;
#endif

#ifdef IW_HAVE_LAMBDA
  auto f = std::find_if(_title.cbegin(), _title.cend(), [&mykey] (const IWString * s1) {return *s1 == mykey;});
  const IWString * f = std::find_if(_title.cbegin(), _title.cend(), [&mykey] (const IWString * s1) {return *s1 == mykey;});

  if (f == _title.cend())
  {
    cerr << "CIF_Loop::string_value:no value for '" << key << "' available\n";
    return 0;
  }

  auto col = f - _title.cbegin();
#else
  int col = -1;

  for (unsigned int i = 0; i < _title.size(); ++i)
  {
    if (*(_title[i]) != mykey)
      continue;

    col = i;
    break;
  }

  if (col < 0)
  {
    cerr << "CIF_Loop::string_value:no value for '" << key << "' available\n";
    return 0;
  }
#endif

  if (! buffer.word(col, v, ' '))
  {
    cerr << "CIF_Loop::string_value:cannot extract column " << col << " from '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
Molecule::read_molecule_cif_ds (iwstring_data_source & input)
{
  input.set_strip_trailing_blanks(1);

  const_IWSubstring buffer;

  bool found_comp_list = false;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (buffer.starts_with("data_comp_list"))
    {
      found_comp_list = true;
      break;
    }
  }

  if (! found_comp_list)
  {
    if (input.eof())
    {
      cerr << "read mol cif eof\n";
      return 0;
    }

    cerr << "Molecule::read_molecule_cif_ds:did not find 'data_comp_list'\n";
    return 0;
  }

  CIF_Loop cifloop;

  if (! cifloop.build(input))
  {
    cerr << "Molecule::read_molecule_cif_ds:could not parse loop structure\n";
    return 0;
  }

  if (0 != cifloop.stem().length())
  {
    cerr << "Molecule::read_molecule_cif_ds:unexpected stem in loop '" << cifloop.stem() << "'\n";
    return 0;
  }

  if (! cifloop.string_value("id", _molecule_name, 0))
  {
    cerr << "Molecule::read_molecule_cif_ds:yipes, cannot extract molecule id\n";
    return 0;
  }

  IWString three_letter_code;
  if (! cifloop.string_value("three_letter_code", three_letter_code, 0))
  {
    cerr << "Molecule::read_molecule_cif_ds:cannot extract three_letter_code\n";
    return 0;
  }

  cerr << "Read '" << _molecule_name << "', three_letter_code '" << three_letter_code << "'\n";

  if (! input.next_record(buffer))
  {
    cerr << "Molecule::read_molecule_cif_ds:premature eof\n";
    return 0;
  }

  IWString should_be;
  should_be << "data_comp_" << three_letter_code;

  if (should_be != buffer)
    cerr << "Molecule::read_molecule_cif_ds:unexpected input record, expected '" << should_be << "', got '" << buffer << "', continuing, but watch for problems\n";

  if (! cifloop.build(input))
  {
    cerr << "Molecule::read_molecule_cif_ds:cannot read atom block\n";
    return 0;
  }

  if ("atom" != cifloop.stem())
  {
    cerr << "Molecule::read_molecule_cif_ds:block not 'atom', got '" << cifloop.stem() << "'\n";
    return 0;
  }

  const auto na = cifloop.dataitems();

  if (0 == na)
  {
    cerr << "Molecule::read_molecule_cif_ds:no atoms\n";
    return 0;
  }

  resize(na);

  int nh = cifloop.number_title_records();

  int comp_id_col = -1;
  int atom_id_col = -1;
  int type_symbol_col = -1;
  int x_col = -1;
  int y_col = -1;
  int z_col = -1;

  for (auto i = 0; i < nh; ++i)
  {
    const IWString & s = cifloop.title(i);
    if ("comp_id" == s)
      comp_id_col = i;
    else if ("atom_id" == s)
      atom_id_col = i;
    else if ("type_symbol" == s)
      type_symbol_col = i;
    else if ("x" == s)
      x_col = i;
    else if ("y" == s)
      y_col = i;
    else if ("z" == s)
      z_col = i;
  }

  if (comp_id_col < 0 || atom_id_col < 0 || type_symbol_col < 0 || x_col < 0 || y_col < 0 || z_col < 0)
  {
    cerr << "Molecule::read_molecule_cif_ds:could not identify one or more atom title columns\n";
    return 0;
  }

  IW_STL_Hash_Map_int atom_name_to_atom_number;

  for (unsigned int i = 0; i < na; ++i)
  {
    const IWString & s = cifloop.data_record(i);

    IWString comp_id, atom_id, type_symbol;
    IWString xs, ys, zs;

    const_IWSubstring token;
    for (auto col = 0, j = 0; s.nextword(token, j); ++col)
    {
      if (comp_id_col == col)
        comp_id = token;
      else if (atom_id_col == col)
        atom_id = token;
      else if (type_symbol_col == col)
        type_symbol = token;
      else if (x_col == col)
        xs = token;
      else if (y_col == col)
        ys = token;
      else if (z_col == col)
        zs = token;
    }

    coord_t x,y,z;
    if (! xs.numeric_value(x) || ! ys.numeric_value(y) || ! zs.numeric_value(z))
    {
      cerr << "Molecule::read_molecule_cif_ds:invalid coordinates '" << s << "'\n";
      return 0;
    }

    const Element * e = get_element_from_symbol_no_case_conversion(type_symbol);
    if (nullptr == e)
    {
      cerr << "Molecule::read_molecule_cif_ds:invalid element specification '" << s << "'\n";
      return 0;
    }

    atom_name_to_atom_number[atom_id] = _number_elements;

    Atom * a = new Atom(e);
    a->setxyz(x,y,z);

    add(a, 1);   // 1 means partial molecule
  }

  while (cifloop.build(input))
  {
    cerr << "Examining stem '" << cifloop.stem() << "'\n";

    if ("bond" == cifloop.stem())
    {
      if (! _cif_bond_list(cifloop, atom_name_to_atom_number))
      {
        cerr << "molecule::read_molecule_cif_ds:invalid bonding\n";
        return 0;
      }
    }
    else if ("chir" == cifloop.stem())
    {
      if (! _cif_chirality(cifloop, atom_name_to_atom_number))
      {
        cerr << "molecule::read_molecule_cif_ds:invalid chirality\n";
        return 0;
      }
    }
  }

  return _number_elements;
}

int
Molecule::_cif_bond_list(const CIF_Loop & cifloop, const IW_STL_Hash_Map_int & atom_name_to_atom_number)
{
  int nt = cifloop.number_title_records();

  int comp_id_col = -1;
  int atom_id_1_col = -1;
  int atom_id_2_col = -1;
  int type_col = -1;

  for (auto i = 0; i < nt; ++i)
  {
    const IWString & s = cifloop.title(i);
    if ("comp_id" == s)
      comp_id_col = i;
    else if ("atom_id_1" == s)
      atom_id_1_col = i;
    else if ("atom_id_2" == s)
      atom_id_2_col = i;
    else if ("type" == s)
      type_col = i;
  }

  if (atom_id_1_col < 0 || atom_id_2_col < 0 || type_col < 0)
  {
    cerr << "Molecule::_cif_bond_list:cannot find bond specification header records\n";
    return 0;
  }

  const auto nb = cifloop.dataitems();

  _bond_list.resize(nb);

  for (unsigned int i = 0; i < nb; ++i)
  {
    const IWString & s = cifloop.data_record(i);

    const_IWSubstring token;

    IWString sa1, sa2, sbt;

    for (int col = 0, j = 0; s.nextword(token, j); ++col)
    {
      if (atom_id_1_col == col)
        sa1 = token;
      else if (atom_id_2_col == col)
        sa2 = token;
      else if (type_col == col)
        sbt =  token;
    }

    (void) comp_id_col;  // quiet unused variable warning.

    auto f = atom_name_to_atom_number.find(sa1);
    if (f == atom_name_to_atom_number.end())
    {
      cerr << "Molecule::_cif_bond_list:unrecognised atoms in bond '" << s << "'\n";
      return 0;
    }

    atom_number_t a1 = f->second;

    f = atom_name_to_atom_number.find(sa2);
    if (f == atom_name_to_atom_number.end())
    {
      cerr << "Molecule::_cif_bond_list:unrecognised atoms in bond '" << s << "'\n";
      return 0;
    }

    atom_number_t a2 = f->second;

    bond_type_t bt;

    if ("single" == sbt)
      bt = SINGLE_BOND;
    else if ("double" == sbt)
      bt = DOUBLE_BOND;
    else if ("triple" == sbt)
      bt = TRIPLE_BOND;
    else
    {
      cerr << "Molecule::_cif_bond_list:unrecognised bond type '" << s << "'\n";
      return 0;
    }

    if (! add_bond(a1, a2, bt, 1))
    {
      cerr << "Molecule::_cif_bond_list:cannot add bond\n";
      return 0;
    }
  }

  return 1;
}

int
Molecule::_cif_chirality (const CIF_Loop & cifloop, const IW_STL_Hash_Map_int & atom_name_to_atom_number)
{
  cerr << "Warning, cif chirality not implemented\n";
  return 1;
}

int
Molecule::write_molecule_cif (std::ostream & os)
{
  cerr << "Molecule::write_molecule_cif:not implemented\n";
  return 0;
}
