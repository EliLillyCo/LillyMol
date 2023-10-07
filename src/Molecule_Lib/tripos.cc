#include <iostream>
#include <iomanip>
#include <memory>
using std::cerr;
using std::endl;

#define COMPILING_TRIPOS_CC

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"

#include "tripos.h"

#include "aromatic.h"
#include "charge_calculation.h"
#include "misc2.h"
#include "molecule.h"
#include "rwmolecule.h"

static int assign_default_formal_charges = 0;

void 
set_mol2_assign_default_formal_charges(int a)
{
  assign_default_formal_charges = a;

  return;
}

static int sybyl_write_assigned_atom_types = 0;

void
set_mol2_write_assigned_atom_types(int s)
{
  sybyl_write_assigned_atom_types = s;
}

static int mol2_write_formal_charge_as_partial_charge = 0;

void
set_mol2_write_formal_charge_as_partial_charge(int s)
{
  mol2_write_formal_charge_as_partial_charge = s;
}

static int mol2_read_charge_column_contains_formal_charges = 0;

void 
set_mol2_read_charge_column_contains_formal_charges(int s)
{
  mol2_read_charge_column_contains_formal_charges = s;
}

/*
  We don't want to burden the Molecule object with tripos residue information, so we can optionally place it 
  in the user specified void ptr
*/

static int place_mol2_residue_information_in_user_specified_void_ptr = 0;

void
set_place_mol2_residue_information_in_user_specified_void_ptr(int s)
{
  place_mol2_residue_information_in_user_specified_void_ptr = s;
}

Tripos_Residue_Information::Tripos_Residue_Information(int na)
{
  _number_residues = na;         // way too many, but avoid resizable arrays
  _residue_name = new IWString[na];         // way too many, but avoid resizable arrays

  _atom_number_to_residue_number = new_int(na, -1);

  return;
}

int
Tripos_Residue_Information::update_residue_information(int anum, int rnum, const IWString & rname)
{
  assert(anum >= 0 && anum < _number_residues);

  _residue_name[anum] = rname;
  _atom_number_to_residue_number[anum] = rnum;

  return 1;
}

int
Tripos_Residue_Information::remove_all_atoms_except(Molecule & m,
                                                     const IWString & rname) const
{
  const auto matoms = m.natoms();

  int * to_revove = new int[matoms]; std::unique_ptr<int[]> free_to_revove(to_revove);

  int rc = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (rname == _residue_name[i])
    {
      to_revove[i] = 0;
      rc++;
    }
    else
      to_revove[i] = 1;
  }

  if (0 == rc)
  {
    cerr << "Tripos_Residue_Information::remove_all_atoms_except:no atoms with residue name '" << rname << "'\n";
    return 0;
  }

  return m.remove_atoms(to_revove);
}

int
Molecule::read_molecule_mol2_ds(iwstring_data_source & input)
{
//input.set_skip_blank_lines(1);    will fail if the name is blank
  input.set_ignore_pattern("^#");

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (0 == buffer.length())
      continue;

    input.push_record();
    break;
  }

  resize(0);

  int rc = _read_molecule_mol2_ds(input);

// always skip to the next molecule

  if (! input.eof())
  {
    if (skip_to_string(input, "@<TRIPOS>MOLECULE", 0)) // (0 != rc)))   // only warn of skipping if error
      input.push_record();    // leave it there for reading next time
  }

  if (assign_default_formal_charges)
    _mol2_assign_default_formal_charges();

  if (nullptr == _atom_type)    // strange
    ;
  else if (_atom_type->number_elements() > 0 && 0 == _atom_type->ztype().length())
    _atom_type->set_type(ATOM_TYPE_SYBYL);

  int d = discern_chirality_from_3d_coordinates();

  if (d && 3 == highest_coordinate_dimensionality())
  {
    if (1 == d)
      (void) discern_chirality_from_3d_structure();
    else if (_number_elements <= d)
      (void) discern_chirality_from_3d_structure();
    else
      cerr << "Molecule::read_molecule_mol2_ds:skipped d@3d for too many atoms '" << name() << "' " << _number_elements << '\n';
  }

  return rc;
}

/*
  March 2004. This function looks very weak.
*/

int
Molecule::_mol2_assign_default_formal_charges()
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    atomic_number_t z = a->atomic_number();

    if (6 == z)
      continue;

    if (1 == z)
      continue;

    if (7 == z)
      ;
    else if (15 == z)
      ;
    else if (16 == z && 3 == a->nbonds())    // took this out, didn't help things
    {
//    a->set_formal_charge(1);
//    rc++;
      continue;
    }
    else if (53 == z && 2 == a->ncon() && 0 == explicit_hydrogens(i))   // -I-   becomes -[I+]-
    {
      a->set_formal_charge(1);
      rc++;
      continue;
    }
    else
      continue;

    if (3 == a->ncon() && 4 == a->nbonds() && ! is_ring_atom(i))
      ;
    else if (4 != a->ncon())
      continue;
    else if (4 != a->nbonds())
      continue;

    a->set_formal_charge(1);

    rc++;
  }

  return rc;
}

int
Molecule::_place_formal_charges_on_quat_n_from_mol2()
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (7 != a->atomic_number())
      continue;

    if (4 != a->ncon())
      continue;

    if (4 == a->nbonds())
    {
      a->set_formal_charge(1);
      rc++;
    }
  }

  return rc;
}


static int
parse_atom_and_bond_record (const_IWSubstring & buffer,
                            int & na,
                            int & nb)
{
  const_IWSubstring token;
  int i = 0;

  if (! buffer.nextword(token, i))
    return 0;

  if (! token.numeric_value(na) || na < 0)
    return 0;

  if (! buffer.nextword(token, i))
    return 0;

  if (! token.numeric_value(nb) || nb < 0)
    return 0;

  return 1;
}

int
Molecule::_read_molecule_mol2_ds (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  EXTRA_STRING_RECORD(input, buffer, "read mol mol2 eof");

  if (! buffer.starts_with("@<TRIPOS>MOLECULE"))
  {
    cerr << "Molecule::_read_molecule_mol2_ds: invalid starting record, line " << input.lines_read() << endl;
    cerr << buffer << endl;

    return 0;
  }

  EXTRA_STRING_RECORD(input, buffer, "cannot read mol2 name record");

  set_name(buffer);

  EXTRA_STRING_RECORD(input, buffer, "cannot read mol2 atom and bond record");

  int na, nb;

  if (! parse_atom_and_bond_record(buffer, na, nb))
  {
    cerr << "Invalid mol2 atom and bond record '" << buffer << "'\n";

    return 0;
  }

  EXTRA_STRING_RECORD(input, buffer, "cannot read mol2 type record");

  IWString charge_type;

  EXTRA_STRING_RECORD(input, charge_type, "cannot read mol2 charge type record");

// We skip the records between here and '@<TRIPOS>ATOM'

  while (1)
  {
    EXTRA_STRING_RECORD(input, buffer, "didn't reach atom table");

    if (buffer.starts_with("@<TRIPOS>ATOM"))
      break;
  }

  if (0 == na && 0 == nb)
  {
    return 1;
  }

  if (! resize(na) || ! _bond_list.resize(nb))
  {
    cerr << "Molecule::_read_molecule_mol2_ds: could not resize for " << na << " atoms and " << nb << " bonds\n";
    return 0;
  }

  int * aromatic_atom  = new_int(na); std::unique_ptr<int[]> free_aromatic_atom(aromatic_atom);

  int * aromatic_bonds;

  if (nb > 0)
    aromatic_bonds = new_int(nb);
  else
    aromatic_bonds = nullptr;

  std::unique_ptr<int[]> free_aromatic_bonds(aromatic_bonds);

  Tripos_Residue_Information * tri = nullptr;

  if (place_mol2_residue_information_in_user_specified_void_ptr)
  {
    tri = new Tripos_Residue_Information(na);
    _user_specified_void_ptr = reinterpret_cast<void *>(tri);
  }

  if (! _read_molecule_mol2_ds(input, na, nb, aromatic_atom, aromatic_bonds, tri))
    return 0;

  if (unconnect_covalently_bonded_non_organics_on_read())
    _do_unconnect_covalently_bonded_non_organics();

  process_delocalised_carbonyl_bonds(aromatic_atom, aromatic_bonds);   // may turn off aromaticity of some of the atoms

//iw_write_array(aromatic_atom, _number_elements, "aromatic atom", cerr);
//iw_write_array(aromatic_bonds, nedges(), "aromatic bonds", cerr);

// If something like 'CN1[N+](C)=NC(C)=N1 p18' is written with Pearlman aromaticity, there is no
// aromaticity, but it is a valence error

  if (put_formal_charges_on_neutral_ND3v4())
    _do_put_formal_charges_on_neutral_ND3v4();

  if (assign_default_formal_charges)
    _place_formal_charges_on_quat_n_from_mol2();

  if (locate_item_in_array(1, na, aromatic_atom) < 0 &&
      locate_item_in_array(1, nb, aromatic_bonds) < 0)    // no aromatic atoms or bonds to worry about
    ;
  else if (find_kekule_form(aromatic_atom, aromatic_bonds))   // great
    ;
  else
  {
    cerr << "Molecule::_read_molecule_mol2_ds: cannot find valid Kekule form\n";
    if (! allow_input_without_valid_kekule_form())
      return 0;
    _molecule_name += " (invalid KEKULE form)";
  }

  if (nullptr != _charges)
    _charges->set_type(charge_type);

#ifdef USE_IWMALLOC
  iwmalloc_check_malloc_magic(aromatic_atom);
  iwmalloc_check_malloc_magic(aromatic_bonds);
#endif

  return 1;
}

int
Molecule::_read_molecule_mol2_ds (iwstring_data_source & input,
                                  int na, int nb,
                                  int * aromatic_atom,
                                  int * aromatic_bond,
                                  Tripos_Residue_Information * tri)
{
  const_IWSubstring buffer;

  for (int i = 0; i < na; i++)
  {
    EXTRA_STRING_RECORD(input, buffer, "eof reading atom table");

    if (! _parse_tripos_atom_record(buffer, i, aromatic_atom[i], tri))
    {
      cerr << "Molecule::_read_molecule_mol2_ds: invalid atom record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  if (0 == nb)   // kind of strange...
    return 1;

  EXTRA_STRING_RECORD(input, buffer, "eof before bonds");

  if (! buffer.starts_with("@<TRIPOS>BOND"))
  {
    cerr << "Molecule::_read_molecule_mol2_ds: should be start of bonds, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  for (int i = 0; i < nb; i++)
  {
    EXTRA_STRING_RECORD(input, buffer, "eof during bond list");

    if (! _parse_tripos_bond_record(buffer, aromatic_atom, i, aromatic_bond))
    {
      cerr << "Molecule::_read_molecule_mol2_ds: invalid bond record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

int
Molecule::_parse_tripos_bond_record (const const_IWSubstring & buffer,
                                     int * aromatic_atoms,
                                     int bond_number,
                                     int * aromatic_bond)
{
  if (buffer.nwords() < 4)
    return 0;

  int i = 0;
  const_IWSubstring token;

  (void) buffer.nextword(token, i);    // skip over the index

  (void) buffer.nextword(token, i);

  atom_number_t a1;

  if (! token.numeric_value(a1) || a1 < 1 || a1 > _number_elements)
  {
    cerr << "MOlecule:_parse_tripos_bond_record: invalid atom 1 '" << token << "'\n";
    return 0;
  }

  a1--;

  (void) buffer.nextword(token, i);

  atom_number_t a2;

  if (! token.numeric_value(a2) || a2 < 1 || a2 > _number_elements)
  {
    cerr << "Molecule:_parse_tripos_bond_record: invalid atom 2 '" << token << "'\n";
    return 0;
  }

  if (a2 == a1 + 1)
  {
    cerr << "Molecule::_parse_tripos_bond_record: self bond found " << buffer << endl;
    if (ignore_self_bonds())
    {
      cerr << "Ignored\n";
      return 1;
    }

    return 0;
  }

  a2--;

  (void) buffer.nextword(token, i);     // the bond type

  bond_type_t bt = INVALID_BOND_TYPE;

  if ('1' == token)
    bt = SINGLE_BOND;
  else if ('2' == token)
    bt = DOUBLE_BOND;
  else if ('3' == token)
    bt = TRIPLE_BOND;
  else if ("ar" == token)
  {
    aromatic_atoms[a1] = 1;
    aromatic_atoms[a2] = 1;
    bt = SINGLE_BOND;
    aromatic_bond[bond_number] = 1;
  }
  else if ("am" == token)    // as near as I can gather, this means amide bond
  {
    bt = SINGLE_BOND;
  }
  else
  {
    cerr << "Molecule::_parse_tripos_bond_record: unrecognised bond type '" << token << "'\n";
    return 0;
  }

  return add_bond(a1, a2, bt, 1);
}

static int
convert_float_to_formal_charge (charge_t q,
                                formal_charge_t & f)
{
  f = static_cast<int>(q);

  if (fabs(static_cast<charge_t>(f) - q) < 1.0e-05)    // great, very close
    return 1;

  if (q > 0.0)     // maybe it was 0.999
    return convert_float_to_formal_charge(q + 0.001, f);
  else     // maybe it was -0.999
    return convert_float_to_formal_charge(q - 0.001, f);
}

/*
  A typical atom record looks like

  29 H29         0.4827   -2.5548    1.6833 H         1 <1>        0.0481
*/

int
Molecule::_parse_tripos_atom_record (const const_IWSubstring & buffer,
                                     atom_number_t zatom,      // which atom number are we processing
                                     int & arom,
                                     Tripos_Residue_Information * tri)
{
  int nw = buffer.nwords();

  if (nw < 6)
    return 0;

  int i = 0;
  const_IWSubstring token;

  (void) buffer.nextword(token, i);     // skip over the sequence number

  (void) buffer.nextword(token, i);     // H29. ignore this, use the 5th token

  coord_t x, y, z;

  (void) buffer.nextword(token, i);     // X coordinate

  if (! token.numeric_value(x))
  {
    cerr << "Molecule::_parse_tripos_atom_record: bad X coordinate\n";
    return 0;
  }

  (void) buffer.nextword(token, i);     // Y coordinate

  if (! token.numeric_value(y))
  {
    cerr << "Molecule::_parse_tripos_atom_record: bad Y coordinate\n";
    return 0;
  }

  (void) buffer.nextword(token, i);     // Y coordinate

  if (! token.numeric_value(z))
  {
    cerr << "Molecule::_parse_tripos_atom_record: bad Z coordinate\n";
    return 0;
  }

  (void) buffer.nextword(token, i);     // atom type

  if (token.ends_with(".ar"))
    arom = 1;
  else   
    arom = 0;

//cerr << "Type '" << token << " arom " << arom << endl;

  int atype = _tripos_atom_type_from_string(zatom, token);

// corina puts the type in the symbol

  (void) token.truncate_at_first('.');     // may do nothing

  const Element * e = get_element_from_symbol_no_case_conversion(token);
  if (nullptr == e)
  {
    if (! auto_create_new_elements())
    {
      cerr << "Molecue::_parse_tripos_atom_record: unrecosnised element '" << token << "'\n";
      return 0;
    }

    e = create_element_with_symbol(token);
    assert (nullptr != e);
  }

  Atom * a = new Atom(e);
  assert (nullptr != a);

  a->setxyz(x, y, z);

  add(a);

  set_atom_type(zatom, atype);

// If there are subsequent columns, we may have charges and/or residue information

  if (! buffer.nextword(token, i))    // we are done. If we are keeping going, token contains thr residue number
    return 1;

  int residue_number;

  if (token.starts_with('*'))
    residue_number = -1;
  else if (! token.numeric_value(residue_number) || residue_number < 1)
  {
    cerr << "Molecue::_parse_tripos_atom_record:invalid residue number '" << buffer << "'\n";
    return 0;
  }

  IWString residue_name;

  if (! buffer.nextword(residue_name, i))
  {
    cerr << "Molecue::_parse_tripos_atom_record:no residue name '" << buffer << "'\n";
    return 0;
  }

  if (nullptr != tri)
    tri->update_residue_information(_number_elements - 1, residue_number, residue_name);

  if (mol2_read_charge_column_contains_formal_charges)
    ;
  else if (nullptr == _charges)
  {
    allocate_charges();

    assert (nullptr != _charges);
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "Molecue::_parse_tripos_atom_record:no partial charge '" << buffer << "'\n";
    return 0;
  }

  float q;
  if (! token.numeric_value(q))
  {
    cerr << "Molecule::_parse_tripos_atom_record: invalid partial charge value '" << buffer << "'\n";
    return 0;
  }

  if (mol2_read_charge_column_contains_formal_charges)
  {
    formal_charge_t f;
    if (! convert_float_to_formal_charge(q, f))
    {
      cerr << "Molecule::_parse_tripos_atom_record: improper formal charge value '" << buffer << "'\n";
      return 0;
    }

    if (! reasonable_formal_charge_value(f))
    {
      cerr << "Molecule::_parse_tripos_atom_record:unreasonable formal charge value '" << buffer << "'\n";
      return 0;
    }

    _things[_number_elements - 1]->set_formal_charge(f);
  }
  else if (! reasonable_atomic_partial_charge_value(q))
  {
    cerr << "Molecule::_parse_tripos_atom_record: unreasonable partial charge value '" << buffer << "'\n";
    return 0;
  }
  else
    _charges->seti(_number_elements - 1, q);      // charge on last atom added

  return 1;
}

int
Molecule::_tripos_atom_type_from_string (atom_number_t zatom,
                             const const_IWSubstring & token)
{
  int rc = UNDEFINED_TRIPOS_ATOM_TYPE;

  if ('H' == token)
    rc = ATOM_TYPE_H;
  else if (token.starts_with("C."))
  {
    if ("C.ar" == token || "C.AR" == token)
      rc = ATOM_TYPE_CAR;
    else if ("C.3" == token)
      rc = ATOM_TYPE_C3;
    else if ("C.2" == token)
      rc = ATOM_TYPE_C2;
    else if ("C.1" == token)
      rc =  ATOM_TYPE_C1;
    else if ("C.cat" == token || "C.CAT" == token)
      rc =  ATOM_TYPE_CCAT;
  }
  else if (token.starts_with("N."))
  {
    if ("N.4" == token)
      rc =  ATOM_TYPE_N4;
    else if ("N.3" == token)
      rc =  ATOM_TYPE_N3;
    else if ("N.2" == token)
      rc =  ATOM_TYPE_N2;
    else if ("N.1" == token)
      rc =  ATOM_TYPE_N1;
    else if ("N.ar" == token || "N.AR" == token)
      rc =  ATOM_TYPE_NAR;
    else if ("N.pl3" == token || "N.PL3" == token)
      rc =  ATOM_TYPE_NPL3;
    else if ("N.am" == token || "N.AM" == token)
      rc =  ATOM_TYPE_NAM;
  }
  else if (token.starts_with("O."))
  {
    if ("O.3" == token)
      rc =  ATOM_TYPE_O3;
    else if ("O.2" == token)
      rc =  ATOM_TYPE_O2;
    else if ("O.co2" == token || "O.CO2" == token)
      rc =  ATOM_TYPE_OCO2;
  }
  else if (token.starts_with("S."))
  {
    if ("S.2" == token)
      rc =  ATOM_TYPE_S2;
    else if ("S.3" == token)
      rc =  ATOM_TYPE_S3;
    else if ("S.o" == token || "S.O" == token)
      rc =  ATOM_TYPE_SO;
    else if ("S.o2" == token || "S.O2" == token)
      rc =  ATOM_TYPE_SO2;
  }
  else if ("P.3" == token)
    rc =  ATOM_TYPE_P3;
  else if ("Br" == token)
    rc =  ATOM_TYPE_BR;
  else if ("Cl" == token || "CL" == token)
    rc =  ATOM_TYPE_CL;
  else if ("F" == token)
    rc =  ATOM_TYPE_F;
  else if ("I" == token)
    rc =  ATOM_TYPE_I;
  else
  {
    cerr << "Molecule::_assign_tripos_atom_type_from_string:unrecognised type '" << token << "'\n";
  }

  return rc;
}

/*
  Common code for writing the atomic symbol in 5 spaces
*/

static int
sybyl_append_atomic_symbol (std::ostream & os, const IWString & asymbol)
{
  os << asymbol;
  if (1 == asymbol.length())
    os << "    ";
  else
    os << "   ";

  return os.good();
}

int
append_sybyl_atom_type (std::ostream & os, int atype, const IWString & asymbol)
{
  switch (atype)
  {
    case ATOM_TYPE_H:
      os << "H    ";
      return os.good();

    case ATOM_TYPE_C3:
      os << "C.3  ";
      return os.good();

    case ATOM_TYPE_C2:
      os << "C.2  ";
      return os.good();

    case ATOM_TYPE_CAR:
      os << "C.ar ";
      return os.good();

    case ATOM_TYPE_C1:
      os << "C.1  ";
      return os.good();

    case ATOM_TYPE_N4:
      os << "N.4  ";
      return os.good();

    case ATOM_TYPE_N3:
      os << "N.3  ";
      return os.good();

    case ATOM_TYPE_N2:
      os << "N.2  ";
      return os.good();

    case ATOM_TYPE_N1:
      os << "N.1  ";
      return os.good();

    case ATOM_TYPE_NAR:
      os << "N.ar ";
      return os.good();

    case ATOM_TYPE_NPL3:
      os << "N.pl3";
      return os.good();

    case ATOM_TYPE_NAM:
      os << "N.am ";
      return os.good();

    case ATOM_TYPE_O3:
      os << "O.3  ";
      return os.good();

    case ATOM_TYPE_O2:
      os << "O.2  ";
      return os.good();

    case ATOM_TYPE_S3:
      os << "S.3  ";
      return os.good();

    case ATOM_TYPE_S2:
      os << "S.2  ";
      return os.good();

    case ATOM_TYPE_P3:
      os << "P.3  ";
      return os.good();

    case ATOM_TYPE_BR:
      os << "Br   ";
      return os.good();

    case ATOM_TYPE_CL:
      os << "Cl   ";
      return os.good();

    case ATOM_TYPE_F :
      os << "F    ";
      return os.good();

    case ATOM_TYPE_I:
      os << "I    ";
      return os.good();

    case ATOM_TYPE_SO:
      os << "S.o  ";
      return os.good();

    case ATOM_TYPE_SO2:
      os << "S.o2 ";
      return os.good();

    case ATOM_TYPE_OCO2:
      os << "O.co2 ";
      return os.good();

    case ATOM_TYPE_CCAT:
      os << "C.cat ";
      return os.good();
    default:
    {
      sybyl_append_atomic_symbol(os, asymbol);

      return 0;
    }
  }

  assert (nullptr == "should not come to here");

  return 1;
}

static int
is_amide (const Molecule & m, const Bond * b, const int * atype)
{
  if (! b->is_single_bond())
    return 0;

  atom_number_t n = b->a1();
  atom_number_t a2 = b->a2();

// Set things up so that N is the nitrogen and A2 is the carbon

  if (ATOM_TYPE_NAM == atype[n] && ATOM_TYPE_C2 == atype[a2])
    ;
  else if (ATOM_TYPE_NAM == atype[a2] && ATOM_TYPE_C2 == atype[n])
    a2 = n;
  else
    return 0;

  const Atom * c = m.atomi(a2);

  assert (6 == c->atomic_number());

  if (3 != c->ncon())
    return 0;

  for (int i = 0; i < 3; i++)
  {
    const Bond * cbi = c->item(i);

    if (! cbi->is_double_bond())
      continue;

    atom_number_t o = cbi->other(a2);

    if (8 == m.atomic_number(o))
      return 1;
  }

  return 0;
}

int
Molecule::write_molecule_mol2 (std::ostream & os)
{
  compute_aromaticity();

  int * atype = new int[_number_elements]; std::unique_ptr<int[]> free_atype(atype);

  if (! find_simplified_sybyl_atom_type_sybyl_style(atype))
  {
    cerr << "Molecule::write_molecule_mol2: warning, some atoms unclassified\n";
  }

//for (int i = 0; i < _number_elements; i++)
//{
//  cerr << "Atom " << i << " '" << smarts_equivalent_for_atom (i) << "' type " << atype[i] << endl;
//}

  return _write_molecule_mol2(os, atype);
}

int
Molecule::_write_molecule_mol2 (std::ostream & os, const int * atype)
{
  if (sybyl_write_assigned_atom_types && nullptr == _atom_type)
  {
    cerr << "MOLECULE::_write_molecule_mol2:no atom types available\n";
    return 0;
  }

  os << "#\n";
  os << "# Written by tripos.cc compiled " << __DATE__ << ' ' << __TIME__ << endl;
  os << "#\n";
  os << "@<TRIPOS>MOLECULE\n";

  os << _molecule_name << '\n';

  int nb = _bond_list.number_elements();

  os <<  std::setw(5) << _number_elements << std::setw(6) << nb << "      1" << '\n';
  os << "SMALL\n";

  if (mol2_write_formal_charge_as_partial_charge)
    os << "USER_CHARGES\n";
  else if (nullptr == _charges)
    os << "NO_CHARGES\n";
  else if (_charges->ztype().length())
    os << _charges->ztype() << '\n';
  else
    os << "UNKNOWN_CHARGE_TYPE\n";

  os << '\n';
  os << '\n';

  os << "@<TRIPOS>ATOM\n";

  for (int i = 0; i < _number_elements; i++)
  {
    os << std::setw(7) << (i + 1) << ' ';      // sequence number

    const Atom * a = _things[i];

    os << a->atomic_symbol() << (i + 1);

    int strlen_atomic_symbol = a->atomic_symbol().length();

    if (1 == strlen_atomic_symbol)
      os << ' ';

    if (i < 9)
      os << ' ';
    if (i < 99)
      os << ' ';

    os << "   ";

    a->write_coordinates (os);

    if (sybyl_write_assigned_atom_types)
      os << 'T' << std::left << std::setw(5) << _atom_type->item(i);
    else if (! append_sybyl_atom_type(os, atype[i], a->atomic_symbol()))
    {
      cerr << "molecule::_write_molecule_mol2: unrecognised type " << atype[i] << endl;
      cerr << "Atom " << i << " type " << a->atomic_symbol() << ", " << a->ncon() << " connections\n";
    }

    if (mol2_write_formal_charge_as_partial_charge || nullptr != _charges)
    {
      os << "     1 <1>       ";

      charge_t q;

      if (mol2_write_formal_charge_as_partial_charge)
        q = static_cast<charge_t>(_things[i]->formal_charge());
      else
        q = _charges->item(i);

#if defined(__GNUG__) || defined(IW_INTEL_COMPILER)
      std::ios::fmtflags old_flags = os.flags(std::ios::fixed);
#else
      long old_flags = os.flags (std::ios::fixed);
#endif
      int old_precision  = os.precision (4);

      os << std::setw(7) << q;

      os.precision(old_precision);
      os.flags(old_flags);
    }

    os << '\n';
  }

  os << "@<TRIPOS>BOND\n";

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    os << std::setw(6) << (i + 1);
      
    os << std::setw(5) << (b->a1() + 1) << std::setw(5) << (b->a2() + 1) << ' ';

    if (b->is_aromatic())
      os << "ar";
    else if (is_amide(*this, b, atype))
    {
      os << "am";
    }
    else if (b->is_single_bond())
      os << '1';
    else if (b->is_double_bond())
      os << '2';
    else if (b->is_triple_bond())
      os << '3';
    else
    {
      cerr << "Molecule::write_molecule_mol2: what kind of bond is this " << b->a1() << " to " << b->a2() << endl;
      os << '?';
    }

    os << '\n';
  }

  os << "@<TRIPOS>SUBSTRUCTURE\n";
  os << "1 ****        1 TEMP              0 ****  ****    0 ROOT\n";


  return os.good();
}

int
Molecule::_do_put_formal_charges_on_neutral_ND3v4()
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];

    if (7 != a->atomic_number())
      continue;

    if (3 != a->ncon())
      continue;

    if (0 != a->formal_charge())
      continue;

    if (4 != a->nbonds())
      continue;

    if (is_ring_atom(i) && _doubly_bonded_to_oxygen(i))
      continue;

    a->set_formal_charge(1);
    rc++;
  }

  if (rc)
    _set_modified();

  return rc;
}

int
Molecule::_doubly_bonded_to_oxygen (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon == a->nbonds())
    return 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(zatom);

    if (8 == _things[o]->atomic_number())
      return 1;
  }

  return 0;
}
