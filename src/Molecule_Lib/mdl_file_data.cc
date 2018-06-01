#include <stdlib.h>

#include "misc.h"
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "mdl_file_data.h"

void
MDL_Atom_Data::_default_values ()
{
  _atom_number = INVALID_ATOM_NUMBER;

  _hcount = 0;
  _h0designator = 0;
  _unsaturated = 0;
  _substitution = 0;
  _ring_bond = 0;
  _valence = 0;

  _atom_map = -1;
  _inversion = 0;
  _exact_change = 0;

// and the settable quantities

  _min_ncon = 0;
  _max_ncon = 0;
  _min_hcount = 0;
  _explicit_hydrogen_atoms_removed = 0;

// We initialise aromatic with a negative value so we can tell
// the difference between 0, which means not aromatic, and 1 which
// means known aromatic. The negative value means not computed

  _aromatic = -1;

  _atom = NULL;

  _nbonds = -1;

//  don't forget the copy constructor (below) if you add things here

  return;
}

MDL_Atom_Data::MDL_Atom_Data()
{
  _default_values();

  return;
}

void
MDL_Atom_Data::_do_copy (const MDL_Atom_Data & rhs)
{
  _atom_number = rhs._atom_number;

  _hcount = rhs._hcount;;
  _h0designator = rhs._h0designator;;
  _unsaturated = rhs._unsaturated;;
  _substitution = rhs._substitution;;
  _ring_bond = rhs._ring_bond;;
  _valence = rhs._valence;;

  _atom_map = rhs._atom_map;;
  _inversion = rhs._inversion;;
  _exact_change = rhs._exact_change;;

  _min_ncon = rhs._min_ncon;;
  _max_ncon = rhs._max_ncon;;
  _min_hcount = rhs._min_hcount;;
  _explicit_hydrogen_atoms_removed = rhs._explicit_hydrogen_atoms_removed;;

  _aromatic = rhs._aromatic;;

  _atom = rhs._atom;;

  _atom_list = rhs._atom_list;

  _alias = rhs._alias;

  _initial_atomic_symbol = rhs._initial_atomic_symbol;
  
  _nbonds = rhs._nbonds;

  return;
}

MDL_Atom_Data::MDL_Atom_Data (const MDL_Atom_Data & rhs)
{
  _do_copy(rhs);
}

MDL_Atom_Data &
MDL_Atom_Data::operator= (const MDL_Atom_Data & rhs)
{
  _do_copy(rhs);

  return *this;
}

int
MDL_Atom_Data::extract_info_from_mdl_file_record(const MDL_Atom_Record & mar)
{
  _hcount = mar.hcount();
  _h0designator = mar.h0designator();
  _valence = mar.valence();

  _atom_map = mar.atom_map();
  _inversion = mar.inversion();
  _exact_change = mar.exact_change();

  _initial_atomic_symbol = mar.atomic_symbol();

  return 1;
}

MDL_Bond_Data::MDL_Bond_Data()
{
  _bond_type_read_in = 0;
  _bond_topology = 0;
  _btype = INVALID_BOND_TYPE;

  _reacting_centre_status = 0;    // not used anywhere

  return;
}

MDL_Bond_Data::MDL_Bond_Data (const MDL_Bond_Data & rhs)
{
  _bond_type_read_in = rhs._bond_type_read_in;
  _bond_topology = rhs._bond_topology;
  _btype = rhs._btype;

  _reacting_centre_status = rhs._reacting_centre_status;

  return;
}

int
MDL_Bond_Data::extract_info_from_mdl_file_record(const MDL_Bond_Record & mbr)
{
  _bond_type_read_in = mbr.bond_type_read_in();

  _bond_topology = mbr.bond_topology();

  if (! mbr.bond_type_for_query(_btype))
  {
    cerr << "MDL_Bond_Data::extract_info_from_mdl_file_record:invalid query bond type\n";
    return 0;
  }

//cerr << "MDL_Bond_Data::extract_info_from_mdl_file_record:got " << _btype << endl;

  _reacting_centre_status = mbr.reacting_center_status();

  return 1;
}

MDL_File_Data::MDL_File_Data()
{
#ifdef USE_IWMALLOC
  cerr << "Checking MDL_File_Data\n";
  iwmalloc_check_all_malloced(stderr);
#endif

  return;
}

int
MDL_File_Data::_do_copy (const MDL_File_Data & rhs)
{
  int na = rhs._mdl_atom.number_elements();

  if (0 == na)
  {
    cerr << "MDL_File_Data::_do_copy:empty rhs\n";
    return 0;
  }

  int nb = rhs._mdl_bond.number_elements();

  for (int i = 0; i < na; i++)
  {
    MDL_Atom_Data * t = new MDL_Atom_Data(*(rhs._mdl_atom[i]));
    _mdl_atom.add(t);
  }

  for (int i = 0; i < nb; i++)
  {
    MDL_Bond_Data * t = new MDL_Bond_Data(*(rhs._mdl_bond[i]));
    _mdl_bond.add(t);
  }

  int nl = _link_atom.number_elements();

  for (int i = 0; i < nl; i++)
  {
    Link_Atom * t = new Link_Atom(*(rhs._link_atom[i]));
    _link_atom.add(t);
  }

  return 1;
}

MDL_File_Data::MDL_File_Data(const MDL_File_Data & rhs)
{
  _do_copy(rhs);

  return;
}

MDL_File_Data::~MDL_File_Data()
{
  return;
}

MDL_File_Data &
MDL_File_Data::operator= (MDL_File_Data && rhs)
{
  _third_line_of_input_sdf_file = std::move(rhs._third_line_of_input_sdf_file);
  _mdl_atom = std::move(rhs._mdl_atom);
  _mdl_bond = std::move(rhs._mdl_bond);
  _link_atom = std::move(rhs._link_atom);

  return *this;
}

int
MDL_File_Data::allocate_arrays(int na, int nb)
{
  _mdl_atom.resize_keep_storage(na);
  _mdl_bond.resize_keep_storage(nb);

  _mdl_atom.resize(na);
  for (int i = 0; i < na; i++)
  {
    _mdl_atom.add(new MDL_Atom_Data);
  }

  if (nb > 0)
  {
    _mdl_bond.resize(nb);

    for (int i = 0; i < nb; i++)
    {
      _mdl_bond.add(new MDL_Bond_Data);
    }
  }

  return 1;
}

int
MDL_Bond_Data::is_or_aromatic() const
{
  if (AROMATIC_BOND == _btype)
    return 0;

  if (_btype | AROMATIC_BOND)
    return 1;

  return 0;
}

int
MDL_File_Data::or_aromatic_bonds_present() const
{
  int rc = 0;

  int nb = _mdl_bond.number_elements();

  for (int i = 0; i < nb; i++)
  {
    const MDL_Bond_Data * b = _mdl_bond[i];

    if (b->is_or_aromatic())
      rc++;
  }

  return rc;
}

int
MDL_File_Data::build (Molecule & m)
{
  _mdl_atom.resize_keep_storage(0);
  _mdl_bond.resize_keep_storage(0);

  int matoms = m.natoms();
  int nbonds = m.nedges();

  _mdl_atom.resize(matoms);
  _mdl_bond.resize(nbonds);

  for (int i = 0; i < matoms; i++)
  {
    MDL_Atom_Data * madi = new MDL_Atom_Data;

    if (m.is_aromatic(i))
      madi->set_aromatic(1);
    else
      madi->set_aromatic(0);

    _mdl_atom.add(madi);
  }

  for (int i = 0; i < nbonds; i++)
  {
    MDL_Bond_Data * mdlb = new MDL_Bond_Data;

    const Bond * b = m.bondi(i);

    mdlb->set_btype(BOND_TYPE_ONLY(b->btype()));

    _mdl_bond.add(mdlb);
  }

  return 1;
}

int
MDL_Atom_Data::build_atom_list(const const_IWSubstring & buffer)
{
  return _atom_list.create_from_ALS_record(buffer);
}

int
MDL_Atom_Data::initialise_atom_list_from_symbol(const const_IWSubstring & token)
{
  _initial_atomic_symbol = "L";

  return _atom_list.initialise_atom_list_from_symbol(token);
}

int
MDL_Atom_Data::convert_a_or_q_atoms_to_atom_list (const IWString & s)
{
//cerr << "MDL_Atom_Data::convert_a_or_q_atoms_to_atom_list:string is '" << s << "'\n";

  if ('A' == s)
    _atom_list.initialise_from_mdl_A_symbol();
  else if ('Q' == s)
    _atom_list.initialise_from_mdl_Q_symbol();
  else if ("AH" == s)
    _atom_list.initialise_from_mdl_AH_symbol();
  else if ("QH" == s)
    _atom_list.initialise_from_mdl_QH_symbol();

  return 1;
}

int
MDL_Atom_Data::convert_not_atom_lists_to_organic_lists ()
{
  if (_atom_list.active() && ! _atom_list.normal_list())
    return _atom_list.convert_not_atom_lists_to_organic_lists();

  return 0;
}

int
MDL_Atom_Data::write_M_ALS(atom_number_t zatom,
                           std::ostream & output) const
{
  return _atom_list.write_M_ALS(zatom, output);
}

int
MDL_Atom_Data::swap_atoms (atom_number_t a1, atom_number_t a2)
{
  if (_atom_number == a1)
    _atom_number = a2;
  else if (_atom_number == a2)
    _atom_number = a1;

  return 1;
}

int
MDL_Atom_Data::connected_atom_is_being_removed (atomic_number_t zremove)
{
  if (_substitution > 0)
    _substitution--;

  if (_min_ncon > 0)
    _min_ncon--;

// Problems with _min_hcount if removing a Hydrogen. Too hard to resolve

  if (_nbonds > 0)
    _nbonds--;     // what if we removed a doubly bonded atom??

  return 1;
}

int
MDL_File_Data::swap_atoms (atom_number_t a1, atom_number_t a2)
{
  int n = _mdl_atom.number_elements();

  if (0 == n)
  {
    cerr << "MDL_File_Data::swap_atoms:no atoms!\n";
    return 0;
  }

  if (a1 < 0 || a2 < 0 || a1 == a2 || a1 >= n || a2 >= n)
  {
    cerr << "MDL_File_Data::swap_atoms:invalid atom number(s) " << a1 << " and " << a2 << " have " << n << " atoms\n";
    return 0;
  }

  for (int i = 0; i < n; i++)
  {
    _mdl_atom[i]->swap_atoms(a1, a2);
  }

  _mdl_atom.swap_elements(a1, a2);

  n = _mdl_bond.number_elements();

  for (int i = 0; i < n; i++)
  {
    _mdl_bond[i]->swap_atoms(a1, a2);
  }

  n = _link_atom.number_elements();

  for (int i = 0; i < n; i++)
  {
    _link_atom[i]->swap_atoms(a1, a2);
  }

  return 0;
}

void
MDL_File_Data::discard_atom_map ()
{
  const auto n = _mdl_atom.number_elements();

  for (int i = 0; i < n; ++i)
  {
    _mdl_atom[i]->set_atom_map(0);
  }

  return;
}
const MDL_Bond_Data *
MDL_File_Data::mdl_bond_data(int i) const
{
  if (! _mdl_bond.ok_index(i))
    return nullptr;
  else
    return _mdl_bond[i];
}

MDL_Bond_Data *
MDL_File_Data::mdl_bond_data(int i)
{ 
  if (! _mdl_bond.ok_index(i))
    return nullptr;
  else
   return _mdl_bond[i];
}

const MDL_Bond_Data *
MDL_File_Data::mdl_bond(int i) const
{
  if (! _mdl_bond.ok_index(i))
    return nullptr;
  else
    return _mdl_bond[i];
}

MDL_Bond_Data *
MDL_File_Data::mdl_bond(int i)
{
  if (! _mdl_bond.ok_index(i))
    return nullptr;
  else
   return _mdl_bond[i];
}

#ifdef IMPLEMENT_THIS_SOMETIME
int
MDL_File_Data::add (const Element * e)
{
  MDL_Atom_Data * t = new MDL_Atom_Data();

  t->set_atom_number(_mdl_atom.number_elements());

  t->set_symbol

  _mdl_atom.add(t);

  return 1;
}
#endif
