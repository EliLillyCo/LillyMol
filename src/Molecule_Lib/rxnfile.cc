#include <memory>
#include <utility>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <exception>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "iw_stl_hash_map.h"

using std::cerr;
using std::endl;

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "misc.h"
#include "iwqsort.h"

#define COMPILING_MDL_CC
#define COMPILING_RXN_FILE

#include "atom.h"
#include "mdl.h"
#include "smiles.h"
#include "chiral_centre.h"
#include "rwmolecule.h"
#include "rxn_file.h"
#include "molecule_to_query.h"
#include "path.h"
#include "atom_typing.h"
#include "msi_object.h"


static int warn_no_mapped_atoms = 1;

void
set_warn_no_mapped_atoms(const int s)
{
  warn_no_mapped_atoms = s;
}

Changing_Atom_Conditions::Changing_Atom_Conditions()
{
  _is_changing_if_different_neighbours = 0;

  _ignore_lost_atom_if_isolated = 0;

  _only_consider_largest_reagent_fragment = 0;

  _include_changing_bonds_in_changing_atom_count = 1;

  _discern_changing_atoms_only_in_first_fragment = 0;

  _consider_aromatic_bonds = 0;

  return;
}

RXN_File::RXN_File ()
{
  _nr = 0;
  _reagent = NULL;

  _np = 0;
  _product = NULL;

  _na = 0;
  _agent = NULL;

  _remove_product_fragments = 0;

  _reagent_locator = NULL;
  _product_locator = NULL;

  _atom_map_to_reagent_atom_number = NULL;
  _atom_map_to_product_atom_number = NULL;

  _initial_bond_type = NULL;
  _final_bond_type = NULL;

  _remove_unmapped_atoms_that_disappear = 0;

  _aromatic_bonds_lose_kekule_identity = 0;

  _preserve_kekule_forms = 0;

  _swap_atoms_to_put_rare_atoms_first = 0;

  _do_automatic_atom_mapping = 1;

  _remove_explicit_hydrogens = 0;

  _involved_in_square_bond = NULL;

  _unconnect_unmapped_atoms_that_exceed_product_valence = 0;

  _convert_A_to_C_for_aromaticity = 0;
    
  _convert_atom_aliases_to_isotopes = 0;

  _interpret_atom_alias_as_smarts = 1;

  _interpret_carbon_on_rhs_of_list_as_no_change = 1;

  _auto_fix_orphans = 0;

//_is_changing_if_different_neighbours = 0;

  _mol2qry_isotope_special_meaning = 0;

  _mark_atoms_changed_when_kekule_form_of_bond_changes = 1;  // by default, consider these atoms changed

  _queryOutStream = NULL;
  
  return;
}

RXN_File::~RXN_File ()
{
  if (NULL != _reagent)
    delete [] _reagent;

  if (NULL != _product)
    delete [] _product;

  if (NULL != _agent)
    delete [] _agent;

  if (NULL != _reagent_locator)
    delete [] _reagent_locator;

  if (NULL != _product_locator)
    delete [] _product_locator;

  if (NULL != _atom_map_to_reagent_atom_number)
    delete [] _atom_map_to_reagent_atom_number;
  if (NULL != _atom_map_to_product_atom_number)
    delete [] _atom_map_to_product_atom_number;

  if (NULL != _initial_bond_type)
    delete [] _initial_bond_type;

  if (NULL != _final_bond_type)
    delete [] _final_bond_type;

  if (NULL != _involved_in_square_bond)
    delete [] _involved_in_square_bond;
 

  return;
}


int
RXN_File::debug_print (std::ostream & os) const
{
  os << "ISIS REACTION with " << _nr << " reagents and " << _np << " products\n";
 
  os << "reagent molecules contain";
  for (int i = 0; i < _nr; ++i)
  {
    os << ' ' << _reagent[i].natoms();
  }
  os << " atoms\n";

  for (int i = 0; i < _nr; i++)
  {
    Molecule & m = _reagent[i];

    if (i > 0)
      os << " + ";
    os << m.smiles();
  }

  os << " -> ";

  for (int i = 0; i < _np; i++)
  {
    Molecule & m = _product[i];

    if (i > 0)
      os << " + ";

    os << m.smiles();
  }

  os << endl;

  for (int i = 0; i < _nr; ++i)
  {
    const ISIS_RXN_FILE_Molecule & r = _reagent[i];

    for (int j = 0; j < r.natoms(); ++j)
    {
      const int mapj = r.atom_map()[j];

      os << " reagent " << i << " atom " << j << " " << r.smarts_equivalent_for_atom(j) << " map " << mapj;
      if (mapj > 0 && NULL != _product_locator)
        os << " in product " << _product_locator[mapj];
      os << "\n";
    }
  }

  return os.good();
}

int
RXN_File::print_atom_map_into (std::ostream & output) const
{
  output << "RXN_File::print_atom_map_into:reaction with " << _nr << " reagents\n";
  for (int i = 0; i < _nr; ++i)
  {
    const int matoms = _reagent[i].natoms();

    for (int j = 0; j < matoms; ++j)
    {
      const int amap = _reagent[i].atom_map(j);

      output << " reagent " << i << " atom " << j << " map " << amap << endl;
      if (_reagent_locator[amap] != i)
        output << "     mismatch in reagent locator, " << _reagent_locator[amap] << endl;
        
    }
  }

  output << "And " << _np << " products\n";

  for (int i = 0; i < _np; ++i)
  {
    const int matoms = _product[i].natoms();

    for (int j = 0; j < matoms; ++j)
    {
      const int amap = _product[i].atom_map(j);
      if (0 == amap)
        continue;

      output << " product " << i << " atom " << j << " map " << amap;
      const int r = _reagent_locator[amap];
      if (r < 0)
        output << " not in reagents\n";
      else
        output << " in reagent " << r << endl;
    }
  }

  return 1;
} 

static int
write_atom_mapped_smiles (const ISIS_RXN_FILE_Molecule & m,
                          std::ostream & os)
{
  Molecule mcopy(m);

  const int matoms = m.natoms();

  const int * amap = m.atom_map();

  if (NULL == amap)
  {
    os << mcopy.smiles();
    return 1;
  }

  for (int i = 0; i < matoms; ++i)
  {
    mcopy.set_isotope(i, amap[i]);
  }

  os << mcopy.smiles();

  return 1;
}

void
RXN_File::set_aromatic_bonds_lose_kekule_identity (int s)
{
  _aromatic_bonds_lose_kekule_identity = s;

  for (int i = 0; i < _nr; ++i)
  {
    _reagent[i].set_aromatic_bonds_lose_kekule_identity(s);
  }

  for (int i = 0; i < _np; ++i)
  {
    _product[i].set_aromatic_bonds_lose_kekule_identity(s);
  }

  return;
}

void
RXN_File::set_preserve_kekule_forms (int s)
{
  _preserve_kekule_forms = s;

  for (int i = 0; i < _nr; ++i)
  {
    _reagent[i].set_preserve_kekule_forms(s);
  }

  for (int i = 0; i < _np; ++i)
  {
    _product[i].set_preserve_kekule_forms(s);
  }

  return;
}

int
RXN_File::write_mapped_reaction (std::ostream & os) const
{
  os << "ISIS REACTION with " << _nr << " reagents and " << _np << " products\n";

  for (int i = 0; i < _nr; i++)
  {
    ISIS_RXN_FILE_Molecule & m = _reagent[i];

    if (i > 0)
      os << " + ";
    write_atom_mapped_smiles(m, os);
  }

  os << " -> ";

  for (int i = 0; i < _np; i++)
  {
    ISIS_RXN_FILE_Molecule & m = _product[i];

    if (i > 0)
      os << " + ";

    write_atom_mapped_smiles(m, os);
  }

  os << endl;

  return os.good();
}

static int
count_atoms (const ISIS_RXN_FILE_Molecule * m,
             int n)
{
  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    rc += m[i].natoms();
  }

  return rc;
}

int
RXN_File::_number_reagent_atoms () const
{
  return count_atoms(_reagent, _nr);
}

int
RXN_File::_number_product_atoms () const
{
  return count_atoms(_product, _np);
}

ISIS_RXN_FILE_Molecule::ISIS_RXN_FILE_Molecule ()
{
  _atom_map = NULL;

  _computed_implicit_hydrogen_count = NULL;
  _explicit_hydrogen_count = NULL;

  _btype = NULL;

  _not_atom_lists_present = 0;

  _explicit_hydrogen_atom_removed = NULL;

  _connections_lost = NULL;

  _convert_A_to_C_for_aromaticity = 0;

  _aromatic_bonds_lose_kekule_identity = 0;

  _preserve_kekule_forms = 0;

  _nbonds = NULL;

  _use_as_representative = NULL;

  _swap_atoms_to_put_rare_atoms_first = 0;

  _molecule_is_sole_reagent = 0;

  _remove_explicit_hydrogens = 0;

  _structure_value = 0;

  _is_duplicate_fragment = 0;

  _interpret_atom_alias_as_smarts = 1;

  _atype = nullptr;

  _aromatic = nullptr;

  return;
}

ISIS_RXN_FILE_Molecule::~ISIS_RXN_FILE_Molecule ()
{
  DELETE_IF_NOT_NULL_ARRAY(_atom_map);
  DELETE_IF_NOT_NULL_ARRAY(_computed_implicit_hydrogen_count);
  DELETE_IF_NOT_NULL_ARRAY(_explicit_hydrogen_atom_removed);
  DELETE_IF_NOT_NULL_ARRAY(_connections_lost);
  DELETE_IF_NOT_NULL_ARRAY(_btype);
  DELETE_IF_NOT_NULL_ARRAY(_nbonds);
  DELETE_IF_NOT_NULL_ARRAY(_use_as_representative);
  DELETE_IF_NOT_NULL_ARRAY(_explicit_hydrogen_count);
  DELETE_IF_NOT_NULL_ARRAY(_atype);

  return;
}

int
ISIS_RXN_FILE_Molecule::allocate_arrays (int na, int nb)
{
  assert(na > 0);

  if (! MDL_Molecule::allocate_arrays(na, nb))
    return 0;

  _atom_map = new_int(na);

  if (0 == nb)
    return 1;

  return 1;
}

int
ISIS_RXN_FILE_Molecule::highest_atom_map_number () const
{
  int matoms = Molecule::natoms();

  if (0 == matoms)
  {
    cerr << "ISIS_RXN_FILE_Molecule:highest_atom_map_number: empty molecule\n";
    return 0;
  }

  int rc = _atom_map[0];
  for (int i = 1; i < matoms; i++)
  {
    if (_atom_map[i] > rc)
      rc = _atom_map[i];
  }

  return rc;
}

template <typename T>
int
resize_array (T * & v,
              const int initial_size,
              const int nextra,
              const T filler)
{
  if (NULL == v)
    return 0;

  T * newv = new T[initial_size + nextra];
  std::copy_n(v, initial_size, newv);
  std::fill_n(newv + initial_size, nextra, filler);

  delete [] v;

  v = newv;

  return 1;
}

int
ISIS_RXN_FILE_Molecule::add (Element const* e)
{
  int initial_matoms = natoms();

  if (! MDL_Molecule::add(e))
    return 0;

  const int matoms = natoms();

  if (NULL == _atom_map)
    _atom_map = new_int(matoms, -1);
  else
    resize_array(_atom_map, initial_matoms, 1, -1);

  if (NULL == _explicit_hydrogen_count)
    _explicit_hydrogen_count = new_int(matoms, 0);
  else
    resize_array(_explicit_hydrogen_count, initial_matoms, 1, 0);

  if (NULL == _computed_implicit_hydrogen_count)
    _computed_implicit_hydrogen_count = new_int(matoms, -1);
  else
    resize_array(_computed_implicit_hydrogen_count, initial_matoms, 1, 0);

  if (NULL == _explicit_hydrogen_atom_removed)
    _explicit_hydrogen_atom_removed = new_int(matoms, -1);
  else
    resize_array(_explicit_hydrogen_atom_removed, initial_matoms, 1, 0);

  if (NULL == _connections_lost)
    _connections_lost = new_int(matoms, -1);
  else
    resize_array(_connections_lost, initial_matoms, 1, 0);

  if (NULL == _btype)
  {
    _btype = new bond_type_t[matoms];
    std::fill_n(_btype, matoms, static_cast<bond_type_t>(INVALID_BOND_TYPE));
  }
  else
    resize_array(_btype, initial_matoms, 1, static_cast<bond_type_t>(INVALID_BOND_TYPE));

  _unique_smiles.resize_keep_storage(0);

  return 1;
}

int
ISIS_RXN_FILE_Molecule::remove_atoms (const int * to_remove)
{
  if (NULL != _atom_map)   // push current atom map back to the atoms in order to save it
  {
    const int matoms = natoms();

    for (int i = 0; i < matoms; ++i)
    {
      _mdl_atom[i]->set_atom_map(_atom_map[i]);
    }
  }

#ifdef DEBUG_REMOVE_ATOMS_QQQQ
  cerr << "ISIS_RXN_FILE_Molecule::remove_atoms:begin with " << natoms() << " atoms\n";
  for (int i = 0; i < natoms(); i++)
  {
    cerr << "      atom " << i << " map " << _atom_map[i] << endl;
  }
  cerr << "REMOVED line " << __LINE__ << endl;
#endif

  const int rc = MDL_Molecule::remove_atoms(to_remove);

  if (NULL == _atom_map)
    return rc;

  const int matoms = natoms();

  for (int i = 0; i < matoms; i++)
  {
    _atom_map[i] = _mdl_atom[i]->atom_map();
    set_atom_map_number(i, _atom_map[i]);
  }

  return rc;
}

/*
  HUH, turns out this was a bad idea, I cannot determine this just by looking at an
  individual component, just return 1
*/

int
ISIS_RXN_FILE_Molecule::has_inter_fragment_changes () const
{
  return 1;
}

/*
  When assembling a reaction, we need to know where each atom is.
*/

int
ISIS_RXN_FILE_Molecule::identify_which_atoms_we_have(int * locator_array,
                                                     const int mark,
                                                     int * atom_map_to_atom_number) const
{
  const int matoms = Molecule::natoms();

  if (0 == matoms)
  {
    cerr << "ISIS_RXN_FILE_Molecule:identify_which_atoms_we_have: empty molecule\n";
    return 0;
  }

  for (int i = 0; i < matoms; i++)
  {
    int m = _atom_map[i];
    if (m <= 0)
      continue;
    locator_array[m] = mark;
    atom_map_to_atom_number[m] = i;
  }

  return matoms;
}

int
ISIS_RXN_FILE_Molecule::fill_bonds_between_mapped_atoms_array(bond_type_t * barray,
                                                              const int dim)
{
//cerr << "ISIS_RXN_FILE_Molecule::fill_bonds_between_mapped_atoms_array:_preserve_kekule_forms " << _preserve_kekule_forms << " _aromatic_bonds_lose_kekule_identity " << _aromatic_bonds_lose_kekule_identity << endl;

  compute_aromaticity_if_needed();

  const int matoms = Molecule::natoms();

  for (int i = 0; i < matoms; i++)
  {
    int mi = _atom_map[i];

    if (0 == mi)
      continue;

    const Atom * a = Molecule::atomi(i);

    for (int j = 0; j < a->ncon(); j++)
    {
      const Bond * b = a->item(j);

      atom_number_t k = b->other(i);

      int mk = _atom_map[k];

      if (0 == mk)
        continue;

      bond_type_t bt;
      if (_preserve_kekule_forms)
        bt = BOND_TYPE_ONLY(b->btype());
      else if (_aromatic_bonds_lose_kekule_identity && b->is_aromatic())
        bt = (AROMATIC_BOND | SINGLE_BOND);
      else
        bt = BOND_TYPE_ONLY(b->btype());

      if (0 == bt)
      {
        cerr << "RXN_File::fill_bonds_between_mapped_atoms_array:unrecognised bond type " << b->btype() << " atoms " << b->a1() << " and " << b->a2() << endl;
        return 0;
      }

//    cerr << "Bond btw mapped atom " << mi << " (" << i << ") and mapped atom " << mk << " (" << k << ") " << bt << endl;

      barray[mi * dim + mk] = bt;
      barray[mk * dim + mi] = bt;
    }
  }

  return 1;
}

/*
  We need to convert between a number in the atom map and the number of a query atom
*/

atom_number_t
ISIS_RXN_FILE_Molecule::which_is_mapped_atom (int m) const
{
  assert (NULL != _atom_map);

  int matoms = Molecule::natoms();
//cerr << "Searching " << matoms << " atoms for map " << m << endl;

  for (int i = 0; i < matoms; i++)
  {
//  cerr << " i = " << i << " check " << _atom_map[i] << " match? " << (_atom_map[i] == m) << endl;
    if (_atom_map[i] == m)
      return i;
  }

  return -1;
}

/*
  We are trying to extend an atom mapping. Does mapped atom M have a single unmapped
  neighbour?

  Ran into issues with Kekule forms. Modify this so that if we are at a 3 connected
  aromatic, and if there are two symmetric equivalent neighbours, go with the one
  with a single bond
*/

//#define DEBUG_IDENTIFY_UNMAPPED_NEIGHBOUR

int
ISIS_RXN_FILE_Molecule::identify_unmapped_neighbour (int m,
                                                     atom_number_t & n,
                                                     int & single_bond_preferentially_chosen)
{
  n = INVALID_ATOM_NUMBER;

  atom_number_t a = which_is_mapped_atom(m);

  if (INVALID_ATOM_NUMBER == a)
  {
    cerr << "ISIS_RXN_FILE_Molecule::identify_unmapped_neighbour:gack, we don't have mapped atom " << m << endl;
    return 0;
  }

  single_bond_preferentially_chosen = 0;

  const Atom * aa = Molecule::atomi(a);

  int is_three_connected_aromatic = ((3 == aa->ncon()) && is_aromatic(a));

#ifdef DEBUG_IDENTIFY_UNMAPPED_NEIGHBOUR
  cerr << "Mapped atom " << m << " is atom " << a << ", z = " << aa->atomic_number() << " " << aa->ncon() << " connections\n";
#endif

  atom_number_t unmapped_singly_bonded_neighbour = INVALID_ATOM_NUMBER;

  int mapped_connections_found = 0;

  for (int i = 0; i < aa->ncon(); i++)
  {
    const Bond * b = aa->item(i);

    atom_number_t j = b->other(a);

    if (_atom_map[j] > 0)    // is mapped
    {
      mapped_connections_found++;
      continue;
    }

    if (! _use_as_representative[j])
    {
      if (is_three_connected_aromatic && INVALID_ATOM_NUMBER == unmapped_singly_bonded_neighbour && b->is_single_bond() && is_aromatic(j))
        unmapped_singly_bonded_neighbour = j;
      continue;
    }

    if (INVALID_ATOM_NUMBER != n)    // have found > 1 unmapped neighbours
      return 0;

    n = j;
  }

  if (INVALID_ATOM_NUMBER == n)
    return 0;

  if (INVALID_ATOM_NUMBER != unmapped_singly_bonded_neighbour && 1 == mapped_connections_found)
  {
    single_bond_preferentially_chosen = 1;
    n = unmapped_singly_bonded_neighbour;
  }

  return 1;
}

int
ISIS_RXN_FILE_Molecule::identify_unmapped_neighbours (int m,
                                                      Set_of_Atoms & nbr)
{
  atom_number_t a = which_is_mapped_atom(m);

  if (INVALID_ATOM_NUMBER == a)
  {
    cerr << "ISIS_RXN_FILE_Molecule::identify_unmapped_neighbours:gack, we don't have mapped atom " << m << endl;
    abort();
    return 0;
  }

  if (0 == _use_as_representative[a])
    return 0;

  const Atom * aa = Molecule::atomi(a);

  const int is_three_connected_aromatic = ((3 == aa->ncon()) && is_aromatic(a));

#ifdef DEBUG_IDENTIFY_UNMAPPED_NEIGHBOUR
  cerr << "Mapped atom " << m << " is atom " << a << ", z = " << aa->atomic_number() << " " << aa->ncon() << " connections\n";
#endif

  atom_number_t unmapped_singly_bonded_neighbour = INVALID_ATOM_NUMBER;
  atom_number_t unmapped_doubly_bonded_neighbour = INVALID_ATOM_NUMBER;

  int mapped_connections_found = 0;

  for (int i = 0; i < aa->ncon(); i++)
  {
    const Bond * b = aa->item(i);

    atom_number_t j = b->other(a);

    if (_atom_map[j] > 0)    // is mapped
    {
      mapped_connections_found++;
      continue;
    }

    if (! _use_as_representative[j])
    {
      if (is_three_connected_aromatic && INVALID_ATOM_NUMBER == unmapped_singly_bonded_neighbour && b->is_single_bond() && is_aromatic(j))
        unmapped_singly_bonded_neighbour = j;
      continue;
    }
    else if (is_three_connected_aromatic && b->is_double_bond())
      unmapped_doubly_bonded_neighbour = j;
    else
      nbr.add(j);
  }

  if (INVALID_ATOM_NUMBER != unmapped_singly_bonded_neighbour)
    nbr.add (unmapped_singly_bonded_neighbour);
  else if (INVALID_ATOM_NUMBER != unmapped_doubly_bonded_neighbour)
    nbr.add(unmapped_doubly_bonded_neighbour);

  return nbr.number_elements();
}

int
ISIS_RXN_FILE_Molecule::identify_singly_bonded_mapped_neighbours (atom_number_t zatom,
                                                    resizable_array<int> & neighbour) const
{
  const Atom * a = atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(zatom);

    if (_atom_map[j] > 0)
      neighbour.add(_atom_map[j]);
  }

  return neighbour.number_elements();
}


int
ISIS_RXN_FILE_Molecule::identify_mapped_neighbours (atom_number_t zatom,
                                                    resizable_array<int> & neighbour) const
{
  const Atom * a = atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (_atom_map[j] > 0)
      neighbour.add(_atom_map[j]);
  }

  return neighbour.number_elements();
}

/*
  We want to skip over the old style atom lists
*/

static IW_Regular_Expression atom_list_rx("[ 0-9][ 0-9][1-9] [F,T]    [2-5] ");

int
ISIS_RXN_FILE_Molecule::do_read(iwstring_data_source & input)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "ISIS_RXN_FILE_Molecule::do_read: cannot read first record\n";
    return 0;
  }

  if (! buffer.starts_with("$MOL"))
  {
    cerr << "ISIS_RXN_FILE_Molecule::do_read: first record must be '$MOL', '" << buffer << " is not\n";
    return 0;
  }

  if (! MDL_Molecule::read_molecule_mdl_ds(input, 1))   // extra arg means return if M  END encountered
  {
    cerr << "ISIS_RXN_FILE_Molecule::do_read:cannot read underlying molecule info\n";
    return 0;
  }

  int contains_explicit_hydrogen = MDL_Molecule::natoms(1);

  int matoms = Molecule::natoms();

  int A_atoms_present = 0;

// If element A has not been created, then there are no A elements

  const Element * notused = get_element_from_symbol_no_case_conversion("A");

  if (NULL != notused)    // Element A has been created, there may be some in our molecule
  {
    for (int i = 0; i < matoms; i++)
    {
      if ("A" == atomic_symbol(i))
        A_atoms_present++;
    }
  }

  int atom_map_present = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (_mdl_atom[i]->atom_map() > 0)
      atom_map_present++;
  }

  if (0 == atom_map_present)
  {
    if (warn_no_mapped_atoms)
      cerr << "ISIS_RXN_FILE_Molecule::do_read: warning, no mapped atoms!\n";
  }

//#define ECHO_CHIRAL_CENTRES
#ifdef ECHO_CHIRAL_CENTRES
  int nc = chiral_centres();

  cerr << "Discerned " << nc << " chiral centres\n";
  for (int i = 0; i < nc; i++)
  {
    const Chiral_Centre * c = chiral_centre_in_molecule_not_indexed_by_atom_number(i);
    c->debug_print(cerr);
  }
#endif

// Look for NOT atom lists

  _not_atom_lists_present = MDL_Molecule::not_atom_lists_present();

  if (_swap_atoms_to_put_rare_atoms_first)
    _do_swap_atoms_to_put_rare_atoms_first();

  if (! _check_non_periodic_table_elements_and_atom_lists())
  {
    cerr << "ISIS_RXN_FILE_Molecule::do_read: non periodic table elements found\n";
    return 0;
  }

  _computed_implicit_hydrogen_count = new int[matoms];
  _explicit_hydrogen_count = new int[matoms];

  _nbonds = new int[matoms];
  for (int i = 0; i < matoms; i++)
  {
    _nbonds[i] = nbonds(i);
  }

// If something has a chiral centre, we may be able to infer something about the bonding
// Turns out this is wrong. The unspecified connection in the product may be an unspecified
// connection in the reagents

#ifdef EXAMINE_CHIRAL_CENTRES_IN_PRODUCTS
  int n = Molecule::chiral_centres();
  for (int i = 0; i < n; i++)
  {
    const Chiral_Centre * c = Molecule::chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    atom_number_t a = c->a();

    if (4 == _nbonds[a])
      continue;

    _nbonds[a] = 4;
  }
#endif

  if (_remove_explicit_hydrogens && contains_explicit_hydrogen)
  {
    MDL_Molecule::remove_explicit_hydrogens(1);
    matoms = Molecule::natoms();
  }

#ifdef DEBUG_TRANSLATE_NP_ELEMENTS
  cerr << "_convert_A_to_C_for_aromaticity " << _convert_A_to_C_for_aromaticity << " A_atoms_present " << A_atoms_present << endl;
#endif

  if (_convert_A_to_C_for_aromaticity && A_atoms_present)
  {
    for (int i = 0; i < matoms; i++)
    {
      const Atom * a = atomi(i);

      if (a->element()->is_in_periodic_table())
        continue;

      if (is_non_ring_atom(i))    // we only change atoms in possibly aromatic rings
        continue;

      if ('A' == a->atomic_symbol() || '*' == a->atomic_symbol())
      {
        Molecule::set_atomic_number(i, 6);

#ifdef DEBUG_TRANSLATE_NP_ELEMENTS
        cerr << "Atom " << i << " set to carbon\n";
#endif
      }
    }
#ifdef DEBUG_TRANSLATE_NP_ELEMENTS
    cerr << "After switching A atoms to C\n";
    for (int i = 0; i < natoms(); i++)
    {
      cerr << "Atom " << i << " type " << atomic_symbol(i) << " aromatic " << is_aromatic(i) << endl;
    }
#endif
  }

  Molecule::compute_aromaticity ();    // may fail if non periodic table elements present

  if (MDL_File_Data::or_aromatic_bonds_present())
    _look_for_rings_that_are_supposed_to_be_aromatic();

  _compute_implicit_hydrogens();

  _unique_smiles = Molecule::unique_smiles();

  if (_aromatic_bonds_lose_kekule_identity)
    _change_explicit_kekule_forms_to_aromatic();

  if (NULL != _atom_map)
    delete [] _atom_map;

  _atom_map = new int[matoms];

  for (int i = 0; i < matoms; i++)
  {
    _atom_map[i] = _mdl_atom[i]->atom_map();
  }

  return 1;
}

void
ISIS_RXN_FILE_Molecule::recompute_unique_smiles ()
{
  _unique_smiles = Molecule::unique_smiles();

  return;
}

/*
  We have read in one Kekule form, but now need to change the _btype array to 
  aromatic rather than single/double
*/

int
ISIS_RXN_FILE_Molecule::_change_explicit_kekule_forms_to_aromatic ()
{
  int nr = Molecule::nrings();

  cerr << "Checking " << nr << " rings for aromatic character\n";

  if (0 == nr)
    return 1;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

//  cerr << "Aromaticity for ring " << (*ri) << endl;

    if (ri->is_aromatic())
      _change_explicit_kekule_forms_to_aromatic(*ri);
  }

  return 1;
}

int
ISIS_RXN_FILE_Molecule::_change_explicit_kekule_forms_to_aromatic (const Ring & r)
{
  for (Ring_Bond_Iterator i(r); i != r.end(); i++)
  {
    atom_number_t a1 = i.a1();
    atom_number_t a2 = i.a2();

    int j = Molecule::which_bond(a1, a2);

    cerr << "Atoms " << a1 << " and " << a2 << " are bond " << j << endl;

    MDL_Bond_Data * mdlb = mdl_bond(j);

    mdlb->set_btype(AROMATIC_BOND);
  }

  return 1;
}

/*
  Our molecule has some aromatic bonds.  first we change all the
  possibly double bonds to double bonds and see if they all become
  aromatic.  Not sure why someone would do this, but I've seen it.
*/

int
ISIS_RXN_FILE_Molecule::_look_for_rings_that_are_supposed_to_be_aromatic ()
{
  int matoms = Molecule::natoms();

  int * aromatic = new_int(matoms); std::unique_ptr<int[]> free_aromatic(aromatic);

  int rc = __look_for_rings_that_are_supposed_to_be_aromatic(aromatic);

  for (int i = 0; i < matoms; i++)
  {
    MDL_Atom_Data * madi = _mdl_atom[i];

    if (1 == aromatic[i])
      madi->set_aromatic(AROMATIC);
//  else
//    madi->set_aromatic(NOT_AROMATIC);
  }

// Make sure all the atoms at the end of possibly aromatic bonds are
// marked aromatic or undetermined

  int nb = Molecule::nedges();

  for (int i = 0; i < nb; i++)
  {
    MDL_Bond_Data * mdlb = _mdl_bond[i];

    if ((DOUBLE_BOND | AROMATIC_BOND) != mdlb->btype())
      continue;

    const Bond * b = Molecule::bondi(i);

    MDL_Bond_Data * bb = _mdl_bond[b->a1()];

    if (NOT_AROMATIC == bb->btype())
      bb->set_btype(AROMATICITY_NOT_DETERMINED);

    bb = _mdl_bond[b->a2()];

    if (NOT_AROMATIC == bb->btype())
      bb->set_btype(AROMATICITY_NOT_DETERMINED);
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::__look_for_rings_that_are_supposed_to_be_aromatic (int * aromatic)
{
  Molecule::aromaticity(aromatic);   // preserve our initial value

  int nb = Molecule::nedges();

  set_display_abnormal_valence_messages(0);

// First try setting all to double bonds

  int bonds_changed = 0;

  for (int i = 0; i < nb; i++)
  {
    const MDL_Bond_Data * mdlb = mdl_bond(i);

    if ((DOUBLE_BOND | AROMATIC_BOND) != mdlb->btype())
      continue;
     
    const Bond * b = Molecule::bondi(i);

    cerr << "Bond between atoms " << b->a1() << " and " << b->a2() << " going to double\n";

    Molecule::set_bond_type_between_atoms(b->a1(), b->a2(), DOUBLE_BOND);

    if (Molecule::valence_ok(b->a1()) && Molecule::valence_ok(b->a2()))
      bonds_changed++;
    else
      Molecule::set_bond_type_between_atoms(b->a1(), b->a2(), SINGLE_BOND);
  }

  set_display_abnormal_valence_messages(1);

  if (0 == bonds_changed)   // strange
    return 1;

  Molecule::compute_aromaticity();

  int aromaticity_gained = 0;
  int aromaticity_lost = 0;

  int matoms = Molecule::natoms();

  for (int i = 0; i < matoms; i++)
  {
    int a = Molecule::is_aromatic(i);

    if (aromatic[i] && a)    // no change, both aromatic
      ;
    else if (0 == aromatic[i] && 0 == a)   // no change
      ;
    else if (aromatic[i] && 0 == a)
      aromaticity_lost++;
    else
      aromaticity_gained++;
  }

//cerr << "Gained " << aromaticity_gained << " Lost " << aromaticity_lost << endl;
  if (aromaticity_lost >= aromaticity_gained)
  {
    _back_to_single_bonds(DOUBLE_BOND | AROMATIC_BOND);
    return 1;
  }

// Excellent, we have gained aromaticity. We should actually check that each of our DOUBLE or AROMATIC bonds have become aromatic, but we don't do that...

  Molecule::aromaticity(aromatic);

  return 1;
}

int
ISIS_RXN_FILE_Molecule::_back_to_single_bonds (bond_type_t bt)
{
  int nb = Molecule::nedges();

  int rc = 0;
  for (int i = 0; i < nb; i++)
  {
    MDL_Bond_Data * mdlb = mdl_bond(i);

    if (bt != mdlb->btype())
      continue;

    const Bond * b = Molecule::bondi(i);

    Molecule::set_bond_type_between_atoms(b->a1(), b->a2(), SINGLE_BOND);
    rc++;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::_identify_explicit_hydrogens_to_be_removed (Set_of_Atoms & to_be_removed,
                                                   int * xref)
{
  int na = Molecule::natoms();

  _explicit_hydrogen_atom_removed = new_int(na);

// First thing is to discern any query info about the atoms that will lose Hydrogens

  int delta = 0;

  for (int i = 0; i < na; i++)
  {
    const Atom * a = Molecule::atomi(i);

    if (1 == a->atomic_number() && 1 == a->ncon())
    {
      atom_number_t o = a->other(i, 0);

      assert(0 >= 0 && o < na);

      _explicit_hydrogen_atom_removed[o]++;
      xref[i] = -1;

      delta--;

      to_be_removed.add(i);
    }
    else
    {
      xref[i] += delta;
    }
  }

  return (0 != delta);    // 0 == delta means nothing happened
}

int
ISIS_RXN_FILE_Molecule::is_an_atom_list (atom_number_t a) const
{
  const MDL_Atom_Data * madi = _mdl_atom[a];

  const ISIS_Atom_List & l = madi->atom_list();

  return l.active();
}

//#define DEBUG_MAPPED_ATOMS_ARE_BONDED

bond_type_t
ISIS_RXN_FILE_Molecule::mapped_atoms_are_bonded (int m1, int m2) const
{
  assert(m1 > 0);
  assert(m2 > 0);

  int matoms = Molecule::natoms();

  atom_number_t a1 = INVALID_ATOM_NUMBER;
  atom_number_t a2 = INVALID_ATOM_NUMBER;

#ifdef DEBUG_MAPPED_ATOMS_ARE_BONDED
  cerr << "Searching " << matoms << " atoms for mapped atoms " << m1 << " and " << m2 << endl;
#endif

  for (int i = 0; i < matoms; i++)
  {
    if (m1 == _atom_map[i])
    {
      assert(INVALID_ATOM_NUMBER == a1);
      a1 = i;
      if (INVALID_ATOM_NUMBER != a2)
        break;
    }
    if (m2 == _atom_map[i])
    {
      assert(INVALID_ATOM_NUMBER == a2);
      a2 = i;
      if (INVALID_ATOM_NUMBER != a1)
        break;
    }
  }

#ifdef DEBUG_MAPPED_ATOMS_ARE_BONDED
  cerr << "Atom numbers " << a1 << " and " << a2 << endl;
#endif

  if (INVALID_ATOM_NUMBER == a1 && INVALID_ATOM_NUMBER == a2)     // must be in different fragments
    return 0;

  if (INVALID_ATOM_NUMBER == a1 || INVALID_ATOM_NUMBER == a2)     // must be in different fragments
  {
    cerr << "ISIS_RXN_FILE_Molecule::mapped_atoms_are_bonded: atom map error " << m1 << " and " << m2 << endl;
    return INVALID_BOND_TYPE;
  }
    
#ifdef DEBUG_MAPPED_ATOMS_ARE_BONDED
  cerr << "Atoms " << a1 << " and " << a2 << " bonded " <<  Molecule::are_bonded(a1, a2) << endl;
#endif

  if (! Molecule::are_bonded(a1, a2))
    return INVALID_BOND_TYPE;

  const Bond * b = Molecule::bond_between_atoms(a1, a2);

  assert (NULL != b);

  if (b->is_single_bond())
    return SINGLE_BOND;
  if (b->is_double_bond())
    return DOUBLE_BOND;
  if (b->is_triple_bond())
    return TRIPLE_BOND;

  cerr << "ISIS_RXN_FILE_Molecule::mapped_atoms_are_bonded: what kind of bond between " << a1 << " and " << a2 << endl;

  return INVALID_BOND_TYPE;
}

int
ISIS_RXN_FILE_Molecule::_grow_symmetry_class(const atom_number_t zatom,
                                             int * symmetry_equivalent_atoms)
{
#ifdef DEBUG_GROW_SYMMETRY_CLASS
  cerr << "Growing symmetry class to atom " << zatom << endl;
#endif

  _use_as_representative[zatom] = 1;
  symmetry_equivalent_atoms[zatom] = 0;

  const int * symm = symmetry_classes();

  const int s = symmetry_class(zatom);

  const int matoms = natoms();

  const Atom * a = atomi(zatom);

  for (int i = 0; i < matoms; i++)     // mark all symmetry equivalent atoms OFF
  {
    if (i == zatom)
      continue;

    if (s != symm[i])
      continue;

//  if (a->is_bonded_to(i))      Mar 2017. Comment out, not sure why this was here. Breaks alkene_oxidation reaction
//    continue;

    _use_as_representative[i] = 0;
    symmetry_equivalent_atoms[i] = 0;
  }

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (symmetry_equivalent_atoms[j] > 1)
      _grow_symmetry_class(j, symmetry_equivalent_atoms);
  }

  return 1;
}

/*
  Major complication. We need to make sure that all the _use_as_representative atoms are together.
*/

int
ISIS_RXN_FILE_Molecule::identify_symmetry_classes()
{
  int matoms = natoms();

  if (0 == matoms)
    return 0;

  if (NULL == _use_as_representative)
    _use_as_representative = new_int(matoms, 1);
  else
    std::fill_n(_use_as_representative, matoms, 1);

  if (number_symmetry_classes() == matoms)     // every atom in its own class
    return 1;

  int * symmetry_equivalent_atoms = new_int(matoms, 1); std::unique_ptr<int[]> free_seq(symmetry_equivalent_atoms);

  const int * symm = symmetry_classes();

  for (int i = 0; i < matoms; i++)
  {
    for (int j = (i+1); j < matoms; j++)
    {

      if (symm[j] == symm[i])
      {
        symmetry_equivalent_atoms[i]++;
        symmetry_equivalent_atoms[j]++;
      }
    }
  }

//#define DEBUG_IDENTIFY_SYMMETRY_CLASSES
#ifdef DEBUG_IDENTIFY_SYMMETRY_CLASSES
  cerr << "Molecule with " << matoms << " determining symmetry\n";
  for (int i = 0; i < matoms; i++)
  {
    cerr << "Atom " << i << " " << smarts_equivalent_for_atom(i) << " " << symmetry_equivalent_atoms[i] << " symmetry equivalent atoms\n";
  }
#endif

  for (int i = 0; i < matoms; i++)
  {
    if (symmetry_equivalent_atoms[i] > 1)    // include all attached atoms
      _grow_symmetry_class(i, symmetry_equivalent_atoms);
  }

//#define DEBUG_IDENTIFY_SYMMETRY_CLASSES
#ifdef DEBUG_IDENTIFY_SYMMETRY_CLASSES
  cerr << "identify_symmetry_classes result\n";
  for (int i = 0; i < matoms; i++)
  {
    cerr << "Atom " << i << " " << smarts_equivalent_for_atom(i) << " representative " << _use_as_representative[i] << endl;
  }
#endif

  return 1;
}

/*int
ISIS_RXN_FILE_Molecule::break_symmetry (atom_number_t a)
{
  int matoms = Molecule::natoms ();

  assert(a >= 0 && a < matoms);

  const int * symm = symmetry_classes();

  int s = symm[a];

  for (int i = 0; i < matoms; i++)
  {
    if (symm[i] == s)
      _use_as_representative[i] = 1;
  }

  return 1;
}*/

void
RXN_File::set_molecule_is_sole_reagent (int r, int s)
{
  assert (r >= 0);

  _molecule_is_sole_reagent[r] = s;

  return;
}

static int
assign_identical_fragments (ISIS_RXN_FILE_Molecule * m, int n)
{
  for (int i = 0; i < n; i++)
  {
    const IWString iusmi = m[i].unique_smiles();

    for (int j = i + 1; j < n; j++)
    {
      if (m[j].is_duplicate_fragment())
        continue;

      if (iusmi == m[j].unique_smiles())
      {
        m[i].set_structure_value(i + 1);
        m[j].set_structure_value(i + 1);
        m[j].set_is_duplicate_fragment(1);
      }
    }
  }

  return 1;
}

int
RXN_File::do_read (iwstring_data_source & input)
{
  input.set_dos(1);    // just in case

  const_IWSubstring buffer;

  if (! input.next_record(buffer))   // normal EOF in a multi-reaction file
  {
    return 0;
  }

  if ("$RXN V3000" == buffer)
  {
    cerr << "RXN_File::do_read:sorry, don't know how to read V3000 reaction files\n";
    cerr << "Try removing highlighting and other advanced features\n";
    return 0;

    return _do_read_v3000(input);
  }

  if (! buffer.starts_with("$RXN"))    // may be DOS 
  {
    cerr << "RXN_File::do_read:invalid first record of RXN file '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(_comment))
  {
    cerr << "RXN_File::do_read: cannot read title\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::do_read: cannot read ISIS record\n";
    return 0;
  }

  if (buffer.contains(" ISIS "))
    ;
  else if (buffer.contains(" Marvin"))
    ;
  else if (buffer.contains(" SMMXDraw "))
    ;
  else if (buffer.contains(" Accelrys "))
    ;
  else if (buffer.contains(" ACCLDraw"))
    ;
  else if (buffer.contains("-NextMove-"))
    ;
  else
    cerr << "RXN_File::do_read: ISIS record possibly invalid '" << buffer << "'\n";
   
  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::do_read: premature eof\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::do_read: premature eof\n";
    return 0;
  }

// Next record should be the number of reagents and products

//cerr << "Discerning numbers from '" << buffer << "'\n";

  _na = 0;

  if (2 == buffer.nwords())
  {
    if (2 != int3d(buffer, _nr, _np))
    {
      cerr << "RXN_File::do_read: invalid Reagents and Products record '" << buffer << "'\n";
      return 0;
    }
  }
  else if (3 == buffer.nwords())
  {
    if (3 != int3d(buffer, _nr, _np, &_na))
    {
      cerr << "RXN_File::do_read: invalid Reagents and Products record '" << buffer << "'\n";
      return 0;
    }
  }

//cerr << "Contains " << _nr << " reagents and " << _np << " products " << _na << " agents\n";

  _reagent = new ISIS_RXN_FILE_Molecule[_nr];
  _product = new ISIS_RXN_FILE_Molecule[_np];
  if (_na > 0)
    _agent = new ISIS_RXN_FILE_Molecule[_na];

  for (int i = 0; i < _nr; i++)
  {
    if (_aromatic_bonds_lose_kekule_identity)
      _reagent[i].set_aromatic_bonds_lose_kekule_identity(1);

    if (_swap_atoms_to_put_rare_atoms_first)
      _reagent[i].set_swap_atoms_to_put_rare_atoms_first(1);

    if (_molecule_is_sole_reagent[i])
      _reagent[i].set_molecule_is_sole_reagent(1);

    if (_remove_explicit_hydrogens)
      _reagent[i].set_remove_explicit_hydrogens(1);

    if (_convert_A_to_C_for_aromaticity)
      _reagent[i].set_convert_A_to_C_for_aromaticity(1);

    if (! _reagent[i].do_read(input))
    {
      cerr << "RXN_File::do_read: cannot read reagent " << i << ", " << input.lines_read() << " lines read\n";
      return 0;
    }

    if (_convert_atom_aliases_to_isotopes)
      _reagent[i].transfer_atom_alias_to_isotope();
  }

  for (int i = 0; i < _np; i++)
  {
    if (_aromatic_bonds_lose_kekule_identity)
      _product[i].set_aromatic_bonds_lose_kekule_identity(1);

    if (_remove_explicit_hydrogens)
      _product[i].set_remove_explicit_hydrogens(1);

    if (_convert_A_to_C_for_aromaticity)
      _product[i].set_convert_A_to_C_for_aromaticity(1);

    if (! _product[i].do_read(input))
    {
      cerr << "RXN_File::do_read: cannot read product " << i << ", " << input.lines_read() << " lines read\n";
      return 0;
    }

    if (_convert_atom_aliases_to_isotopes)
      _product[i].transfer_atom_alias_to_isotope();
//  cerr << "Product smiles " << _product[i].smiles() << endl;
  }

  if (_na > 0)
  {
    for (int i = 0; i < _na; ++i)
    {
      if (! _agent[i].do_read(input))
      {
        cerr << "RXN_File::do_read:cannot read agent " << i << ", " << input.lines_read() << " lines read\n";
        return 0;
      }
    }
  }

  if (_nr > 1)
    assign_identical_fragments(_reagent, _nr);

  if (_np > 1)
    assign_identical_fragments(_product, _np);

  return 1;
}

int
RXN_File::do_read (const const_IWSubstring & f)
{
  iwstring_data_source input(f);

  if (! input.ok())
  {
    cerr << "RXN_File::do_read: cannot open '" << f << "'\n";
    return 0;
  }

  _fname = f;

  if (do_read(input))
    return 1;

  cerr << "RXN_File::do_read: fatal error reading '" << f << "'\n";
  return 0;
}

/*
$RXN V3000

      ISIS     052720051317

M  V30 COUNTS 4 4
M  V30 BEGIN REACTANT
M  V30 BEGIN CTAB
M  V30 COUNTS 18 18 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.7431 -0.312461 0 0
M  V30 2 C -5.74425 -1.13983 0 0
M  V30 3 C -5.02944 -1.5527 0 0
M  V30 4 C -4.313 -1.13935 0 0
M  V30 5 C -4.31586 -0.308828 0 0
M  V30 6 C -5.03125 0.100306 0 0
M  V30 7 N -3.59788 -1.55073 0 0 CFG=3
M  V30 8 C -3.60294 0.106361 0 0
M  V30 9 O -2.88692 -0.303444 0 0
M  V30 10 O -3.60605 0.931355 0 0
M  V30 11 C -2.88498 -1.13552 0 0
M  V30 12 C -3.59567 -2.37573 0 0
M  V30 13 O -4.30903 -2.79014 0 0
M  V30 14 O -2.8801 -2.78631 0 0
M  V30 15 C -2.87788 -3.61131 0 0 CFG=3
M  V30 16 C -2.16231 -4.02189 0 0
M  V30 17 C -3.59125 -4.02572 0 0
M  V30 18 C -2.15442 -3.21388 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 15 17
M  V30 2 1 15 18
M  V30 3 1 5 8
M  V30 4 1 2 3
M  V30 5 1 8 9
M  V30 6 2 8 10
M  V30 7 2 3 4
M  V30 8 1 7 11
M  V30 9 1 4 5
M  V30 10 1 7 12
M  V30 11 2 12 13
M  V30 12 2 5 6
M  V30 13 1 12 14
M  V30 14 1 6 1
M  V30 15 1 14 15
M  V30 16 2 1 2
M  V30 17 1 15 16
M  V30 18 1 4 7
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(18 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18) -
M  V30 BONDS=(18 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18)
M  V30 END COLLECTION
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.2375 -0.808333 0 0
M  V30 2 N 0.5875 -0.808333 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(2 1 2) BONDS=(1 1)
M  V30 END COLLECTION
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.6375 -0.8125 0 0
M  V30 2 C 3.4625 -0.8125 0 0
M  V30 3 O 3.875 -0.0980291 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(3 1 2 3) BONDS=(2 1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 8.03333 -0.795833 0 0 CHG=1
M  V30 2 C 8.85833 -0.795833 0 0 CHG=-1
M  V30 3 C 7.20833 -0.795833 0 0
M  V30 4 C 6.79707 -0.0795405 0 0
M  V30 5 C 5.97567 -0.0775696 0 0
M  V30 6 C 5.55931 -0.790207 0 0
M  V30 7 C 5.97057 -1.5065 0 0
M  V30 8 C 6.7982 -1.51016 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 3 1 2
M  V30 2 1 1 3
M  V30 3 2 3 4
M  V30 4 1 3 8
M  V30 5 1 4 5
M  V30 6 1 5 6
M  V30 7 1 6 7
M  V30 8 1 7 8
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(3 1 2 3) BONDS=(8 1 2 3 4 5 6 7 8)
M  V30 END COLLECTION
M  V30 END CTAB
M  V30 END REACTANT
M  V30 BEGIN PRODUCT
M  V30 BEGIN CTAB
M  V30 COUNTS 16 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 11.7444 -0.2291 0 0
M  V30 2 C 11.7433 -1.0565 0 0
M  V30 3 C 12.4581 -1.4694 0 0
M  V30 4 C 12.4563 0.1836 0 0
M  V30 5 C 13.1765 -1.0536 0 0
M  V30 6 C 13.1717 -0.2255 0 0
M  V30 7 C 13.8188 0.3004 0 0
M  V30 8 N 13.8327 -1.5681 0 0 CFG=3
M  V30 9 N 14.6309 0.1261 0 0 CFG=3
M  V30 10 C 14.6434 -1.375 0 0
M  V30 11 C 14.9957 -0.6203 0 0 CFG=3
M  V30 12 O 13.6263 1.1026 0 0
M  V30 13 O 15.1647 -2.0145 0 0
M  V30 14 C 15.8207 -0.6129 0 0
M  V30 15 C 13.8354 -2.40733 0 0
M  V30 16 C 15.1392 0.775913 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 6 4
M  V30 2 1 7 9
M  V30 3 1 4 1
M  V30 4 1 8 10
M  V30 5 1 9 11
M  V30 6 1 10 11
M  V30 7 1 2 3
M  V30 8 2 7 12
M  V30 9 1 5 6
M  V30 10 2 10 13
M  V30 11 2 3 5
M  V30 12 1 11 14
M  V30 13 1 6 7
M  V30 14 2 1 2
M  V30 15 1 8 15
M  V30 16 1 5 8
M  V30 17 1 9 16
M  V30 END BOND
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 13.396 -3.80073 0 0
M  V30 2 O 12.6826 -4.21514 0 0
M  V30 3 O 14.1116 -4.21131 0 0
M  V30 4 C 14.1138 -5.03631 0 0 CFG=3
M  V30 5 C 14.8294 -5.44689 0 0
M  V30 6 C 13.4004 -5.45072 0 0
M  V30 7 C 14.9331 -4.74722 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 4 6
M  V30 2 1 4 7
M  V30 3 2 1 2
M  V30 4 1 1 3
M  V30 5 1 3 4
M  V30 6 1 4 5
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(7 1 2 3 4 5 6 7) BONDS=(6 1 2 3 4 5 6)
M  V30 END COLLECTION
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 16.8042 -4.4 0 0
M  V30 END ATOM
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 15.3833 -7.1875 0 0 CHG=1
M  V30 2 C 14.5583 -7.1875 0 0
M  V30 3 C 14.1471 -6.47121 0 0
M  V30 4 C 13.3257 -6.46924 0 0
M  V30 5 C 12.9093 -7.18187 0 0
M  V30 6 C 13.3206 -7.89817 0 0
M  V30 7 C 14.1482 -7.90182 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 2 7
M  V30 4 1 3 4
M  V30 5 1 4 5
M  V30 6 1 5 6
M  V30 7 1 6 7
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(2 1 2) BONDS=(7 1 2 3 4 5 6 7)
M  V30 END COLLECTION
M  V30 END CTAB
M  V30 END PRODUCT
M  END
*/

int
RXN_File::_do_read_v3000(iwstring_data_source & input)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::_do_read_v3000:second record missing\n";
    return 0;
  }

  if (buffer.length() > 0)
    cerr << "RXN_File::_do_read_v3000:second record not blank '" << buffer << "'\n";

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::_do_read_v3000:ISIS record not present\n";
    return 0;
  }

  if (! buffer.contains(" ISIS "))
  {
    cerr << "RXN_File::_do_read_v3000:record must contains ISIS '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::_do_read_v3000:second blank line missing\n";
    return 0;
  }

  if (buffer.length())
    cerr << "RXN_File::_do_read_v3000:not blank '" << buffer << "'\n";

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::_do_read_v3000:counts line missing\n";
    return 0;
  }

// M  V30 COUNTS 4 4

  if (! buffer.starts_with("M  V30 COUNTS"))
  {
    cerr << "RXN_File::_do_read_v3000:must be counts '" << buffer << "'\n";
    return 0;
  }


  buffer.remove_leading_chars(static_cast<int>(::strlen("M  V30 COUNTS")));

  if (2 != int3d(buffer, _nr, _np))
  {
    cerr << "RXN_File::_do_read_v3000:invalid Reagents and Products record '" << buffer << "'\n";
    return 0;
  }

//cerr << "Contains " << _nr << " reagents and " << _np << " products\n";

  _reagent = new ISIS_RXN_FILE_Molecule[_nr];
  _product = new ISIS_RXN_FILE_Molecule[_np];

  if (! input.next_record(buffer))
  {
    cerr << "RXN_File::_do_read_v3000:reactant record not present\n";
    return 0;
  }

  if ("M  V30 BEGIN REACTANT" != buffer)
  {
    cerr << "RXN_File::_do_read_v3000:not reactant '" << buffer << "'\n";
    return 0;
  }

  return 0;
}

/*
  We need to encode a reagent number and atom number in a single int.
  We use base 200000 arithmetic, which will of course break if we ever
  get a reagent with more than 200000 atoms
*/

#define REAGENT_ATOM_NUMBER_BASE 200000


template <typename C> int
identify_unmapped_atoms (ISIS_RXN_FILE_Molecule * m,
                         int n,
                         std::unordered_map<int, int> & t,
                         C & c)
{
//cerr << "identify_unmapped_atoms, processing " << n << " molecules\n";

  int rc = 0;     // the number of unmapped atoms

  for (int i = 0; i < n; i++)
  {
    ISIS_RXN_FILE_Molecule & mi = m[i];

    if (mi.is_duplicate_fragment())
      continue;

    const int matoms = mi.natoms();

    for (int j = 0; j < matoms; j++)
    {
//    cerr << " rgnt " << 0 << " atom " << j << " map " << mi.atom_map(j) << " symmetry " << mi.use_atom_as_symmetry_representative(j) << ' ' << mi.smarts_equivalent_for_atom(j) << endl;
      if (mi.atom_map(j) > 0)
        continue;

      if (mi.is_an_atom_list(j))
        continue;

      if (! mi.use_atom_as_symmetry_representative(j))
        continue;

      const auto h = c(mi, j);

//    cerr << i << " j = " << j << " h = " << h << " atom type " << mi.smarts_equivalent_for_atom(j) << endl;

      std::unordered_map<int, int>::iterator f = t.find(h);

      if (f == t.end())      // first occurrence
        t[h] = i * REAGENT_ATOM_NUMBER_BASE + j;
      else if ((*f).second >= 0)     // second or greater occurrence, we are only interested in finding elements with single instances
        (*f).second = -1;

      rc++;
    }
  }

//cerr << "Found " << rc << " unmapped atoms\n";
  return rc;
}

class Just_Atomic_Number
{
  private:
  public:
    int operator () (const ISIS_RXN_FILE_Molecule &, atom_number_t) const;
};

int
Just_Atomic_Number::operator () (const ISIS_RXN_FILE_Molecule & m, atom_number_t a) const
{
  const Element * e = m.elementi(a);

  return e->atomic_symbol_hash_value();
}

class Atomic_Number_and_Aromaticity
{
  private:
  public:
    int operator () (ISIS_RXN_FILE_Molecule &, atom_number_t) const;
};

int
Atomic_Number_and_Aromaticity::operator () (ISIS_RXN_FILE_Molecule & m, atom_number_t a) const
{
  const Element * e = m.elementi(a);

  int rc = e->atomic_symbol_hash_value() * 2;

  if (m.is_aromatic(a))
    rc += 1;

  return rc;
}

class Atomic_Number_and_Symmetry
{
  private:
  public:
    int operator() (ISIS_RXN_FILE_Molecule & m, atom_number_t) const;
};

int
Atomic_Number_and_Symmetry::operator() (ISIS_RXN_FILE_Molecule & m, const atom_number_t zatom) const
{
  const Element * e = m.elementi(zatom);

  int rc = e->atomic_symbol_hash_value() * 2;

  Set_of_Atoms symm;

  m.symmetry_equivalents(zatom, symm);

  const int ns = symm.number_elements();
//cerr << " atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << " has " << ns << " symmetry equivalents " << symm << endl;
  if (0 == ns)
    return rc;

  for (int i = 0; i < ns; ++i)
  {
    const atom_number_t j = symm[i];

    if (zatom < j)
      return rc + i + 1;
  }

  return rc + ns + 1;
}

class Connections_and_Bonds
{
  private:
  public:
    int operator () (ISIS_RXN_FILE_Molecule &, atom_number_t) const;
};

/*
  We generate arbitrary numbers, based on the bonds and connected atoms. We need
  a means of assigning those to a unique number - that's what the hash_map does
*/

class Neighbours_and_Connections
{
  private:
    int _natoms;

    int * _multiplier;

    int * _values_found;

    IW_Hash_Map<int, int> _value_to_unique_id;

  public:
    Neighbours_and_Connections ();
    ~Neighbours_and_Connections ();

    int set_natoms (int n);
    int set_max_connectivity(int s);

    int operator () (const ISIS_RXN_FILE_Molecule &, atom_number_t);
};

Neighbours_and_Connections::Neighbours_and_Connections ()
{
  _natoms = 0;
  _multiplier = NULL;

  return;
}

Neighbours_and_Connections::~Neighbours_and_Connections ()
{
  if (NULL != _multiplier)
    delete [] _multiplier;

  if (NULL != _values_found)
    delete [] _values_found;

  return;
}

int
Neighbours_and_Connections::set_natoms (int n)
{
  assert (0 == _natoms);
  assert (n > 0);

  _multiplier = new int[n];

  _values_found = NULL;

  _multiplier[0] = 1;

  for (int i = 1; i < n; i++)
  {
    _multiplier[i] = _multiplier[i - 1] * n;
  }

  _natoms = n;

  return 1;
}

int
Neighbours_and_Connections::set_max_connectivity(const int s)
{
  assert (NULL == _values_found);

  _values_found = new int[s];

  return 1;
}

/*
*/

int
Neighbours_and_Connections::operator () (const ISIS_RXN_FILE_Molecule & m,
                                         atom_number_t zatom)
{
  assert (NULL != _values_found);

//cerr << "Neighbours_and_Connections::operator: atom " << zatom << " hash contains " << _value_to_unique_id.size() << endl;

  const Atom * a = m.atomi(zatom);

  const int ncon = a->ncon();

  for (int i = 0; i < ncon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    const Atom * aj = m.atomi(j);

    int h = aj->element()->atomic_symbol_hash_value();

    int t;
    if (b->is_aromatic())
      t = 4 * h;
    else
      t = b->number_of_bonds() * h;

    int v;

    IW_Hash_Map<int, int>::const_iterator f = _value_to_unique_id.find(t);
    if (f == _value_to_unique_id.end())
    {
      v = static_cast<int>(_value_to_unique_id.size());
      _value_to_unique_id[t] = v;
    }
    else
      v = (*f).second;

    assert (v >= 0 && v < _natoms);

    _values_found[i] = v;
  }

  assert (_natoms > 0);

  if (ncon > 1)
    std::sort(_values_found, _values_found + ncon);

// We have a sequence of up to 4 numbers, in the range 1-natoms

  int rc = 0;

  for (int i = 0; i < ncon; i++)
  {
    const int v = _values_found[i];

    rc += _multiplier[v];
  }

  return rc;
}

int
Connections_and_Bonds::operator () (ISIS_RXN_FILE_Molecule & m, atom_number_t zatom) const
{
  const Atom * a = m.atomi(zatom);

  int z = a->element()->atomic_symbol_hash_value();

  int ncon = a->ncon();

  int nbonds = a->nbonds();

  int rc = 2 * (500 * z + 6 * ncon + nbonds);    // could conceivably overflow

  if (m.is_aromatic(zatom))
    rc += 1;

  return rc;
}

static int
bond_constant (const Bond & b)
{
  if (b.is_aromatic())
    return 4;
  if (b.is_single_bond())
    return 1;
  if (b.is_double_bond())
    return 2;

  return 3;
}

static int
ncon2 (const ISIS_RXN_FILE_Molecule & m, 
       const atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int rc = 0;

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const auto j = a->other(zatom, i);

    rc += 500 * m.ncon(j) + m.atom_map(j);
  }

//cerr << " ncon2 returning " << rc << endl;
  return rc;
}

class Shell_Based
{
  private:
    int _radius;
  public:
    Shell_Based();

    void set_radius (int s) { _radius = s;}

    int operator () (ISIS_RXN_FILE_Molecule &, atom_number_t) const;
};

Shell_Based::Shell_Based()
{
  _radius = 1;
}

int
Shell_Based::operator() (ISIS_RXN_FILE_Molecule & m, const atom_number_t zatom) const
{
  m.compute_aromaticity_if_needed();

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  const int * amap = m.atom_map();

  uint64_t rc = a->atomic_number();

  if (m.is_aromatic(zatom))
    rc = rc * 81 + 30;

  rc += 7 * acon + 2 * a->nbonds() + m.formal_charge(zatom);

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    const atom_number_t j = b->other(zatom);

    const int jtype = 100002 * m.atomic_number(j) + 10000 * amap[j] * 2000 * m.ncon(j);

    rc += bond_constant(*b) * jtype * ncon2(m, j);
  }

  return static_cast<int>(rc);
}

/*
  We need to decode a number encoded with base REAGENT_ATOM_NUMBER_BASE
*/

static void
identify_reagent_and_atom_number (int t,
                                  int & m,
                                  int & a)
{
  m = t / REAGENT_ATOM_NUMBER_BASE;

  a = t % REAGENT_ATOM_NUMBER_BASE;

  return;
}

//#define DEBUG_ESTABLISH_ATOM_MAPPING

/*
  We have an unmapped, or poorly mapped reaction. See if we can
  establish some atom mapping
*/

int
RXN_File::_establish_atom_mapping (int & highest_atom_map_number)
{
// First check elements. If there is only one occurrence of an unmapped element in the products, we can make a match

  std::unordered_map<int, int> unmapped_occurrences_in_reagents;
  std::unordered_map<int, int> unmapped_occurrences_in_products;

  Just_Atomic_Number just_atomic_number;    // function object

  _unmapped_atoms_in_reagents = identify_unmapped_atoms(_reagent, _nr, unmapped_occurrences_in_reagents, just_atomic_number);
  _unmapped_atoms_in_products = identify_unmapped_atoms(_product, _np, unmapped_occurrences_in_products, just_atomic_number);

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "unmapped_atoms_in_reagents just atomic number " << _unmapped_atoms_in_reagents << " unmapped_atoms_in_products " << _unmapped_atoms_in_products << endl;
#endif

  if (0 == _unmapped_atoms_in_products)    // did we remove all unmapped atoms in the reagents
    return 0;

  if (_establish_atom_mapping(unmapped_occurrences_in_reagents, unmapped_occurrences_in_products, highest_atom_map_number))
    return 1;

  if (0 == _unmapped_atoms_in_reagents)
    return 1;

  unmapped_occurrences_in_reagents.clear();
  unmapped_occurrences_in_products.clear();

  Atomic_Number_and_Aromaticity atomic_number_and_aromaticity;

  _unmapped_atoms_in_reagents = identify_unmapped_atoms(_reagent, _nr, unmapped_occurrences_in_reagents, atomic_number_and_aromaticity);
  _unmapped_atoms_in_products = identify_unmapped_atoms(_product, _np, unmapped_occurrences_in_products, atomic_number_and_aromaticity);

  if (_establish_atom_mapping(unmapped_occurrences_in_reagents, unmapped_occurrences_in_products, highest_atom_map_number))
    return 1;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "unmapped_atoms_in_reagents atomic number and aromaticity " << _unmapped_atoms_in_reagents << " unmapped_atoms_in_products " << _unmapped_atoms_in_products << endl;
#endif

  if (0 == _unmapped_atoms_in_reagents)
    return 1;

// Let's try to break symmetry early

  unmapped_occurrences_in_reagents.clear();
  unmapped_occurrences_in_products.clear();

#ifdef DOES_SYMMETRY_HERE_REALLY_HELP
  Atomic_Number_and_Symmetry atomic_number_and_symmetry;

  _unmapped_atoms_in_reagents = identify_unmapped_atoms(_reagent, _nr, unmapped_occurrences_in_reagents, atomic_number_and_symmetry);
  _unmapped_atoms_in_products = identify_unmapped_atoms(_product, _np, unmapped_occurrences_in_products, atomic_number_and_symmetry);

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "unmapped_atoms_in_reagents atomic number and symmetry " << _unmapped_atoms_in_reagents << " unmapped_atoms_in_products " << _unmapped_atoms_in_products << endl;
#endif

  if (0 == _unmapped_atoms_in_reagents)
    return 1;
#endif

  unmapped_occurrences_in_reagents.clear();
  unmapped_occurrences_in_products.clear();

  Connections_and_Bonds connections_and_bonds;

  _unmapped_atoms_in_reagents = identify_unmapped_atoms(_reagent, _nr, unmapped_occurrences_in_reagents, connections_and_bonds);
  _unmapped_atoms_in_products = identify_unmapped_atoms(_product, _np, unmapped_occurrences_in_products, connections_and_bonds);

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "unmapped_atoms_in_reagents connections and bonds " << _unmapped_atoms_in_reagents << " unmapped_atoms_in_products " << _unmapped_atoms_in_products << endl;
#endif

  if (0 == _unmapped_atoms_in_reagents)
    return 1;

  if (_establish_atom_mapping(unmapped_occurrences_in_reagents, unmapped_occurrences_in_products, highest_atom_map_number))
    return 1;

  if (0 == _unmapped_atoms_in_reagents)
    return 1;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "RXN_File::_establish_atom_mapping: " << _unmapped_atoms_in_reagents << " unmapped reagents, " << _unmapped_atoms_in_products << " unmapped products\n";
#endif

  Neighbours_and_Connections neighbours_and_connections;

//const int reagent_atoms = _number_reagent_atoms();
//const int product_atoms = _number_product_atoms();

  neighbours_and_connections.set_natoms(_number_reagent_atoms() + _number_product_atoms());
  neighbours_and_connections.set_max_connectivity(20);    // surely that's enough

  unmapped_occurrences_in_reagents.clear();
  unmapped_occurrences_in_products.clear();

  _unmapped_atoms_in_reagents = identify_unmapped_atoms(_reagent, _nr, unmapped_occurrences_in_reagents, neighbours_and_connections);
  _unmapped_atoms_in_products = identify_unmapped_atoms(_product, _np, unmapped_occurrences_in_products, neighbours_and_connections);

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "unmapped_atoms_in_reagents neighbours and connections " << _unmapped_atoms_in_reagents << " unmapped_atoms_in_products " << _unmapped_atoms_in_products << endl;
#endif

  if (0 == _unmapped_atoms_in_reagents)
    return 1;

  if (_establish_atom_mapping(unmapped_occurrences_in_reagents, unmapped_occurrences_in_products, highest_atom_map_number))
    return 1;

  if (0 == _unmapped_atoms_in_reagents)
    return 1;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "RXN_File::_establish_atom_mapping: " << _unmapped_atoms_in_reagents << " unmapped reagents, " << _unmapped_atoms_in_products << " unmapped products\n";
#endif

  Shell_Based sb;
//cerr << "Begin shell based expansions\n";

  for (int i = 0; i < 6; ++i)
  {
    unmapped_occurrences_in_reagents.clear();
    unmapped_occurrences_in_products.clear();

    _unmapped_atoms_in_reagents = identify_unmapped_atoms(_reagent, _nr, unmapped_occurrences_in_reagents, sb);
    _unmapped_atoms_in_products = identify_unmapped_atoms(_product, _np, unmapped_occurrences_in_products, sb);

    if (0 == _unmapped_atoms_in_reagents)
      return 1;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
    cerr << "RXN_File::_establish_atom_mapping: shell based " << i << ' ' << _unmapped_atoms_in_reagents << " unmapped reagents, " << _unmapped_atoms_in_products << " unmapped product atoms, trying shell based\n";
#endif

    if (_establish_atom_mapping(unmapped_occurrences_in_reagents, unmapped_occurrences_in_products, highest_atom_map_number))
      return 1;

    if (0 == _unmapped_atoms_in_reagents)
      return 1;
  }

  return 1;
}

//#define DEBUG_SET_ATOM_MAP

int
RXN_File::_set_atom_map (int & highest_atom_map_number,
                         int r,
                         atom_number_t ar,
                         int p,
                         atom_number_t ap)
{
  highest_atom_map_number++;

  assert(_reagent[r].atom_map(ar) <= 0);
  assert(_product[p].atom_map(ap) <= 0);

  _reagent[r].set_atom_map(ar, highest_atom_map_number);
  _product[p].set_atom_map(ap, highest_atom_map_number);

#ifdef DEBUG_SET_ATOM_MAP
  cerr << "Assign map " << highest_atom_map_number << " for reagent " << r << " atom " << ar << " and product " << p << " atom " << ap << endl;
#endif

  _reagent_locator[highest_atom_map_number] = r;
  _product_locator[highest_atom_map_number] = p;

  return 1;
}

static int
break_symmetry (ISIS_RXN_FILE_Molecule * m,     // reagents or products
                const int n,                          // number reagents or products
                const int f,
                const atom_number_t a)
{
  int s = m[f].symmetry_class(a);
//cerr << "break_symmetry, class " << s << " from fragment " << f << ", atom " << a << endl;

  for (int i = 0; i < n; i++)
  {
    if (i == f)
      ;
    else if (m[f].structure_value() > 0 && m[i].structure_value() == m[f].structure_value())
      ;
    else
      continue;

    m[i].set_is_duplicate_fragment(0);

    int matoms = m[i].natoms();

    for (int j = 0; j < matoms; j++)
    {
      if (s == m[i].symmetry_class(j))
        m[i].set_use_as_symmetry_representative(j, 0);
    }
  }

  return 1;
}

static int
identify_symmetry_equivalents(ISIS_RXN_FILE_Molecule * m,
                              const int n,
                              resizable_array_p<Atom_in_Fragment> & af,
                              const int f,
                              const atom_number_t a)
{
  ISIS_RXN_FILE_Molecule & mf = m[f];

  const int s = mf.symmetry_class(a);

  for (int i = 0; i < n; i++)
  {
    ISIS_RXN_FILE_Molecule & mi = m[i];

    if (i == f)      // molecule F, definitely do it
      ;
    else if (mf.structure_value() > 0 && mi.structure_value() == mf.structure_value())   // and all that are identical
      ;
    else
      continue;

    const int matoms = mi.natoms();

    for (int j = 0; j < matoms; j++)
    {
      if (mi.atom_map(j) > 0)
        continue;

      if (s == mi.symmetry_class(j))
      {
        Atom_in_Fragment * tmp = new Atom_in_Fragment(j, i);
        af.add(tmp);
      }
    }
  }

  return 1;
}

//#define DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
static int
write_mapped_neighbours (const char * s,
                         ISIS_RXN_FILE_Molecule * m,
                         int f,
                         atom_number_t a,
                         const resizable_array<int> & nbr,
                         std::ostream & os)
{
  os << s << " fragment " << f << " atom " << a << ' ';

  ISIS_RXN_FILE_Molecule & mf = m[f];

  os << mf.smarts_equivalent_for_atom(a);

  os << " has " << nbr.number_elements() << " mapped neighbours";

  for (int i = 0; i < nbr.number_elements(); i++)
  {
    os << ' ' << nbr[i];
  }

   os << '\n';

  return os.good();
}
#endif

static int
atoms_in_common (const resizable_array<int> & m1,
                 const resizable_array<int> & m2)
{
  int n = m1.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    if (m2.contains(m1[i]))
      rc++;
  }

  return rc;
}

int
RXN_File::_map_symmetry_equivalent_atoms(int & highest_atom_map_number,
                                         int r,
                                         atom_number_t ar,
                                         int p,
                                         atom_number_t ap)
{
#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
  cerr << "SYMM:Atom " << ar << " in reagent " << r << " " << _reagent[r].smarts_equivalent_for_atom(ar) << " mapped with atom " << ap << " in product molecule " << p << " " << _product[p].smarts_equivalent_for_atom(ap) << ", map " << (highest_atom_map_number + 1) << endl;
#endif

  resizable_array_p<Atom_in_Fragment> s1, s2;

  identify_symmetry_equivalents(_reagent, _nr, s1, r, ar);
  identify_symmetry_equivalents(_product, _np, s2, p, ap);

  const int n1 = s1.number_elements();
  const int n2 = s2.number_elements();

  assert(n1 > 0);
  assert(n2 > 0);

  if (1 == n1 && 1 == n2)    // great, no symmetry equivalents, just the atoms themselves
  {
    _set_atom_map(highest_atom_map_number, r, ar, p, ap);
    return 1;
  }

#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
  cerr << s1.number_elements() << " and " << s2.number_elements() << " symmetry equivalents\n";
  cerr << "atom " << ar << " reagent " << r;
  for (int i = 0; i < n1; ++i)
  {
    cerr << " atom " << s1[i]->atom() << " frag " << s1[i]->fragment();
  }
  cerr << endl;
  cerr << "atom " << ap << " product " << p;
  for (int i = 0; i < n2; ++i)
  {
    cerr << " atom " << s2[i]->atom() << " frag " << s2[i]->fragment();
  }
  cerr << endl;
#endif

  int need_to_break_product_symmetry = 0;
  int need_to_break_reagent_symmetry = 0;

  if (n1 != n2)    // different symmetry in product, or symmetry destroyed
  {
//  cerr << "Break symmetry for reagent " << r << " atom " << ar << " " << rm.atomic_symbol(ar) << endl;
//  cerr << "Break symmetry for product " << p << " atom " << ap << " " << pm.atomic_symbol(ap) << endl;
    need_to_break_reagent_symmetry = 1;
    need_to_break_product_symmetry = 1;
  }

// We need to match these up. The items in S1 don't necessarily match with the corresponding items in S2
// We need to identify the best pair-wise matches

  int * pair_wise_score = new_int(n1 * n1 + n2); std::unique_ptr<int[]> free_pws(pair_wise_score);

  for (int i = 0; i < n1; i++)
  {
    const Atom_in_Fragment * airf = s1[i];        // atom in reagent fragment

    int rf = airf->fragment();           // reagent fragment

    atom_number_t ra = airf->atom();           // reagent atom

    resizable_array<int> reagent_mapped_atoms_attached;
    _reagent[rf].identify_mapped_neighbours(ra, reagent_mapped_atoms_attached);

#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
    write_mapped_neighbours("Reagent", _reagent, rf, ra, reagent_mapped_atoms_attached, cerr);
#endif

    for (int j = 0; j < n2; j++)       // shouldn't bother doing this of reagent_mapped_atoms_attached is empty
    {
      Atom_in_Fragment * aipf = s2[j];        // atom in product fragment

      atom_number_t pa = aipf->atom();

//    cerr << "j = " << j << " atom " << pa << " in product fragment " << aipf->fragment() << endl;
      
      if (INVALID_ATOM_NUMBER == pa)    // already matched
        continue;

      int pf = aipf->fragment();

      resizable_array<int> product_mapped_atoms_attached;

      _product[pf].identify_mapped_neighbours(pa, product_mapped_atoms_attached);

#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
      write_mapped_neighbours("Product", _product, pf, pa, product_mapped_atoms_attached, cerr);
#endif

      int c = atoms_in_common(reagent_mapped_atoms_attached, product_mapped_atoms_attached);

      pair_wise_score[i * n1 + j] = c;

#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
      cerr << "reagent atom " << ra << " product atom " << pa << ",  " << c << " atoms in common\n";
#endif
    }
  }

// Now identify the best pair-wise matches

#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
  cerr << "Looking for best pair-wise matches, n1 = " << n1 << " n2 = " << n2 << endl;
#endif

  while (1)
  {
    int best_score = -1;
    int best_i = -1;
    int best_j = -1;

    for (int i = 0; i < n1; i++)
    {
      if (INVALID_ATOM_NUMBER == s1[i]->atom())    // already paired
        continue;

      for (int j = 0; j < n2; j++)
      {
        if (INVALID_ATOM_NUMBER == s2[j]->atom())    // already paired
          continue;

//      cerr << "j = " << j << " product atom " << s1[i]->atom() << " and reagent atom " << s2[j]->atom() << " score " << pair_wise_score[i * n1 + j] << endl;

        if (pair_wise_score[i * n1 + j] > best_score)
        {
          best_score = pair_wise_score[i * n1 + j];
          best_i = i;
          best_j = j;
        }
      }
    }

#ifdef DEBUG_MAP_SYMMETRY_EQUIVALENT_ATOMS
    cerr << "Highest score " << best_score << " best_i " << best_i << " best_j " << best_j << endl;
    if (best_i >= 0 && best_j >= 0)
      cerr << " atoms " << s1[best_i]->atom() << " and " << s2[best_j]->atom() << " atom map to assign " << highest_atom_map_number << endl;
#endif

    if (best_score < 0)
      break;

    _set_atom_map(highest_atom_map_number, s1[best_i]->fragment(), s1[best_i]->atom(), 
                                           s2[best_j]->fragment(), s2[best_j]->atom());

    s1[best_i]->set_atom(INVALID_ATOM_NUMBER);
    s2[best_j]->set_atom(INVALID_ATOM_NUMBER);

    pair_wise_score[best_i * n1 + best_j] = -2;
  }

  if (need_to_break_reagent_symmetry)
    break_symmetry(_reagent, _nr, r, ar);
    
  if (need_to_break_product_symmetry)
    break_symmetry(_product, _np, p, ap);

  return 1;
}

//#define DEBUG_ESTABLISH_ATOM_MAPPING

/*
  We have used identify_unmapped_atoms to identify possibly common
  atoms in reagents and products.  Map together any unique matches.
  Take symmetry into account
*/

int
RXN_File::_establish_atom_mapping(const IW_Hash_Map<int, int> & unmapped_occurrences_in_reagents,
                                  const IW_Hash_Map<int, int> & unmapped_occurrences_in_products,
                                  int & highest_atom_map_number)
{
  int rc = 0;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
  cerr << "Hashes have " << unmapped_occurrences_in_reagents.size() << " and " << unmapped_occurrences_in_products.size() << " members\n";
  for (auto f : unmapped_occurrences_in_reagents)
  {
    cerr <<  " reagent " << f.first << ' ' << f.second << endl;
  }
  for (auto f : unmapped_occurrences_in_products)
  {
    cerr <<  " product " << f.first << ' ' << f.second << endl;
  }
#endif

  for (IW_Hash_Map<int, int>::const_iterator fr = unmapped_occurrences_in_reagents.begin(); fr != unmapped_occurrences_in_reagents.end(); ++fr)
  {
//  cerr << "Reagent key " << (*fr).first << " value " << (*fr).second << endl;;
    if ((*fr).second < 0)    // more than one occurrence
      continue;

    int h = (*fr).first;

    IW_Hash_Map<int, int>::const_iterator fp = unmapped_occurrences_in_products.find(h);
    if (fp == unmapped_occurrences_in_products.end())    // appeared in reagents, but not products
      continue;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
    cerr << "In second hash " << (*fp).second << endl;
#endif

    if ((*fp).second < 0)    // more than one occurrence in the products
      continue;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
    IWString tmp;
    symbol_for_atomic_symbol_hash_value(h, tmp);
    cerr << "Possible atom map for element hash " << h << " element '" << tmp << "'\n";
#endif

    int mr, ar;
    identify_reagent_and_atom_number((*fr).second, mr, ar);
    assert(mr >= 0 && mr < _nr);

    if (_reagent[mr].atom_map(ar) > 0)
      continue;

    int mp, ap;
    identify_reagent_and_atom_number((*fp).second, mp, ap);
    assert (mp >= 0 && mp < _np);

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
    cerr << "Current map for atom " << ap << " in product " << mp << " " << _product[mp].atom_map(ap) << endl;
#endif

    if (_product[mp].atom_map(ap) > 0)
      continue;

#ifdef DEBUG_ESTABLISH_ATOM_MAPPING
    cerr << "Atom " << ar << " " << _reagent[mr].smarts_equivalent_for_atom(ar) << " in reagent " << mr << " mapped to ";
    cerr << "atom " << ap << " " << _product[mp].smarts_equivalent_for_atom(ap) << " in product " << mp << ", map number " << highest_atom_map_number << " hash " << (*fp).first << endl;
#endif

//  Now we must consider any symmetry equivalent atoms

    _map_symmetry_equivalent_atoms(highest_atom_map_number, mr, ar, mp, ap);

    rc++;
  }

  return rc;
}

/*
  In the absence of serious re-arrangements, we may be able to do more
  aggressive atom mapping.

  Beware!
*/

int
RXN_File::_aggressively_map_unmapped_atoms (int highest_atom_map_number)
{
  return 1;
}

/*
  We report a failure if we find a fragment with 0 mapped atoms
*/

static int
molecules_with_at_least_one_mapped_atom(const ISIS_RXN_FILE_Molecule * m,
                                        const int n)
{
  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    const ISIS_RXN_FILE_Molecule & mi = m[i];

    int tmp = count_non_zero_occurrences_in_array(mi.atom_map(), mi.natoms());
    if (tmp)
      rc++;
  }
  
  return rc;
}

static int
molecules_with_at_least_one_inter_fragment_change(const ISIS_RXN_FILE_Molecule * m,
                                                  int n)
{
  int rc = 0;

  for (int i = 1; i < n; i++)
  {
    rc += m[i].has_inter_fragment_changes();
  }

  return rc;
}
static int
report_mapping_result(ISIS_RXN_FILE_Molecule * m,
                      const int n,
                      std::ostream & output)
{
  int atoms_in_set = 0;
  int atoms_mapped = 0;

  for (int i = 0; i < n; i++)
  {
    const ISIS_RXN_FILE_Molecule & mi = m[i];

    const int matoms = mi.natoms();

    atoms_in_set += matoms;

    const int atoms_mapped_this_molecule = count_non_zero_occurrences_in_array(mi.atom_map(), matoms);

//  output << "Component " << i << " " << matoms << " atoms, " << atoms_mapped_this_molecule << " atoms mapped\n";

    atoms_mapped += atoms_mapped_this_molecule;
  }

  output << n << " molecules, " << atoms_in_set << " atoms, " << atoms_mapped << " atoms mapped\n";

  return output.good();
}

static int
look_very_similar(ISIS_RXN_FILE_Molecule & m1,
                  const atom_number_t a1,
                  ISIS_RXN_FILE_Molecule & m2,
                  const atom_number_t a2)
{
  const Element * re = m1.elementi(a1);
  const Element * pe = m2.elementi(a2);

  if (re->atomic_symbol_hash_value() != pe->atomic_symbol_hash_value())
    return 0;

  if (m1.ncon(a1) != m2.ncon(a2))
    return 0;

  const int arom1 = m1.is_aromatic(a1);
  const int arom2 = m2.is_aromatic(a2);

  if (arom1 != arom2)
    return 0;

  if (m1.is_ring_atom(a1) != m2.is_ring_atom(a2))
    return 0;

  if (m1.nbonds(a1) != m2.nbonds(a2))
    return 0;

  return 1;
}

//#define DEBUG__MAP_UNMAPPED_ATOMS

/*
  Let's see if we can extend the atom mapping. Can we identify mapped atoms
  that have 1 unmapped connection in both the reagent and the product.
*/

int
RXN_File::__map_unmapped_atoms(int & highest_atom_map_number)
{
  for (int i = 1; i <= highest_atom_map_number; i++)
  {
    int ri = _reagent_locator[i];
    if (ri < 0)
      continue;

    ISIS_RXN_FILE_Molecule & r = _reagent[ri];

#ifdef DEBUG__MAP_UNMAPPED_ATOMS
    cerr << "Mapped atom " << i << " should be in reagent " << ri << endl;
#endif

    atom_number_t unmapped_neighbour_reagent;
    int single_bond_preferentially_chosen_reagent;

    if (! r.identify_unmapped_neighbour(i, unmapped_neighbour_reagent, single_bond_preferentially_chosen_reagent))
      continue;

    int pi = _product_locator[i];
    if (pi < 0)
    {
      cerr << "RXN_File::__map_unmapped_atoms:strange, mapped atom " << i << " in reagent " << ri << " not in any product " << _comment << " ignored\n";
      continue;
    }

#ifdef DEBUG__MAP_UNMAPPED_ATOMS
    cerr << "Mapped atom " << i << " should be in product " << pi << endl;
#endif

    ISIS_RXN_FILE_Molecule & p = _product[pi];

    atom_number_t unmapped_neighbour_product;
    int single_bond_preferentially_chosen_product;

    if (! p.identify_unmapped_neighbour(i, unmapped_neighbour_product, single_bond_preferentially_chosen_product))
      continue;

//  Compare via element symbol hash values so we can do things like R1

    const Element * re = r.elementi(unmapped_neighbour_reagent);
    const Element * pe = p.elementi(unmapped_neighbour_product);

    if (re->atomic_symbol_hash_value() != pe->atomic_symbol_hash_value())
      continue;

//  cerr << "From unmapped neighbours to mapped atom " << i << endl;

    if (single_bond_preferentially_chosen_reagent && single_bond_preferentially_chosen_product)
      _set_atom_map(highest_atom_map_number, ri, unmapped_neighbour_reagent, pi, unmapped_neighbour_product);
    else
      _map_symmetry_equivalent_atoms(highest_atom_map_number, ri, unmapped_neighbour_reagent, pi, unmapped_neighbour_product);

    return highest_atom_map_number;
  }

//cerr << "Looking for mapped atoms with multiple unmapped neighbours\n";

// We looked for cases where there was a single unmapped atom hanging off a mapped atom. Now the more
// complex case where there might be more than 1 unmapped atom hanging off

  for (int i = 1; i <= highest_atom_map_number; i++)
  {
    int ri = _reagent_locator[i];
    if (ri < 0)
      continue;

    ISIS_RXN_FILE_Molecule & r = _reagent[ri];

//  cerr << "Mapped atom " << i << " should be in reagent " << ri << endl;

    Set_of_Atoms unmapped_neighbour_reagent;

    if (! r.identify_unmapped_neighbours(i, unmapped_neighbour_reagent))
      continue;

    int pi = _product_locator[i];
    if (pi < 0)
    {
      cerr << "RXN_File::__map_unmapped_atoms:strange, mapped atom " << i << " in reagent " << ri << " not in products " << _comment << " ignored\n";
      continue;
    }

//  cerr << "Mapped atom " << i << " should be in product " << pi << endl;

    ISIS_RXN_FILE_Molecule & p = _product[pi];

    Set_of_Atoms unmapped_neighbour_product;
    if (! p.identify_unmapped_neighbours(i, unmapped_neighbour_product))
      continue;

//  cerr << "Adjacent to mapped atom " << i << " reagents " << unmapped_neighbour_reagent.number_elements() << " neighbours, " << unmapped_neighbour_product.number_elements() << " around product\n";

//  do all pair-wise comparisons to see if we get a match. Make sure there is just one match to each item in unmapped_neighbour_reagent

    for (int j = 0; j < unmapped_neighbour_reagent.number_elements(); j++)
    {
      atom_number_t aj = unmapped_neighbour_reagent[j];

      atom_number_t matching_product_atom = INVALID_ATOM_NUMBER;
      int item_to_remove = -1;

      for (int k = 0; k < unmapped_neighbour_product.number_elements(); k++)
      {
        atom_number_t ak = unmapped_neighbour_product[k];

//      cerr << "Check similarity between atom " << aj << " in reagent " << ri << " '" << _reagent[ri].smarts_equivalent_for_atom(aj) << " and atom " << unmapped_neighbour_product[0] << " in product " << pi << " " << _product[pi].smarts_equivalent_for_atom(ak) << endl;

        if (! look_very_similar(_reagent[ri], aj, _product[pi], ak))
          continue;

        if (INVALID_ATOM_NUMBER != matching_product_atom)    // already got a match, can't handle multiple matches
        {
          matching_product_atom = INVALID_ATOM_NUMBER;
          break;
        }
        
        matching_product_atom = ak;
        item_to_remove = k;
      }

      if (INVALID_ATOM_NUMBER != matching_product_atom)
      {
        _map_symmetry_equivalent_atoms(highest_atom_map_number, ri, aj, pi, matching_product_atom);
        return highest_atom_map_number;
        unmapped_neighbour_product.remove_item(item_to_remove);
      }
    }
  }

  return highest_atom_map_number;
}

/*
  Try to extend the mapping iteratively until no change
*/

int
RXN_File::_map_unmapped_atoms(int & highest_atom_map_number)
{
  _map_any_atom_lists_on_either_side(highest_atom_map_number);

  while (1)
  {
    const int h = highest_atom_map_number;

    __map_unmapped_atoms(highest_atom_map_number);

// We may have the case of [mapped]-[unmapped]-[mapped]

    _map_unmapped_surrounded_by_mapped(highest_atom_map_number);

    _establish_atom_mapping(highest_atom_map_number);

    if (h == highest_atom_map_number)    // no change, we are done
      return 1;
  }
}

/*
  Stricly speaking, this is incorrect. We should check throughout the reagents
  that there is only one combination of a particular element and attached mapped
  atoms. Let's hope that's too rare to worry about
*/

int
RXN_File::_map_unmapped_surrounded_by_mapped(int & highest_atom_map_number)
{
  for (int i = 0; i < _nr; i++)
  {
    ISIS_RXN_FILE_Molecule & ri = _reagent[i];

    int matoms = ri.natoms();

    for (int j = 0; j < matoms; j++)
    {
      if (ri.atom_map(j) > 0)
        continue;

//    if (! ri.use_atom_as_symmetry_representative(j))    // maybe the reaction destroys symmetry
//      continue;

      resizable_array<int> mapped_neighbours;

      if (0 == ri.identify_mapped_neighbours(j, mapped_neighbours))
        continue;

      int h = ri.elementi(j)->atomic_symbol_hash_value();

      int pm;
      atom_number_t pa;

      if (! _identify_product_with_same_mapped_neighbours(h, mapped_neighbours, pm, pa))
        continue;

//    cerr << "From being surrounded\n";

      _map_symmetry_equivalent_atoms(highest_atom_map_number, i, j, pm, pa);
    }
  }

  return 1;
}

static int
lists_the_same (const resizable_array<int> & l1,
                const resizable_array<int> & l2)
{
  int n = l1.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (! l2.contains(l1[i]))
      return 0;
  }

  return 1;
}

int
RXN_File::_identify_product_with_same_mapped_neighbours (int h,
                                                         const resizable_array<int> & mapped_neighbours,
                                                         int & pm,
                                                         atom_number_t & pa) const
{
#ifdef DEBUG_PRODUCT_WITH_SAME_NEIGHBOURS
  cerr << "Looking for product with h = " << h << " nbrs";
  for (int i = 0; i < mapped_neighbours.number_elements(); i++)
  {
    cerr << ' ' << mapped_neighbours[i];
  }
  cerr << endl;
#endif

  pm = -1;
  pa = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _np; i++)
  {
    const ISIS_RXN_FILE_Molecule & m = _product[i];

    int matoms = m.natoms();

    for (int j = 0; j < matoms; j++)
    {
      if (m.atom_map(j) > 0)    // we are looking for unmapped atoms
        continue;

      if (! m.use_atom_as_symmetry_representative(j))
        continue;

      if (m.elementi(j)->atomic_symbol_hash_value() != h)
        continue;

      resizable_array<int> mn;

      if (m.identify_mapped_neighbours(j, mn) != mapped_neighbours.number_elements())
        continue;

      if (! lists_the_same(mapped_neighbours, mn))
        continue;

      if (INVALID_ATOM_NUMBER != pa)
        return 0;

      pm = i;
      pa = j;
    }
  }

  return INVALID_ATOM_NUMBER != pa;
}

int
RXN_File::_map_any_atom_lists_on_either_side (int & highest_atom_map_number)
{
  int rc = 0;

//cerr << "Looking to map corresponding atom lists\n";

  for (int i = 0; i < _nr; i++)
  {
    ISIS_RXN_FILE_Molecule & ri = _reagent[i];

    int matoms = ri.natoms();

    for (int j = 0; j < matoms; j++)
    {
      if (ri.atom_map(j) > 0)    // already mapped, we are looking for unmapped lists
        continue;

      const ISIS_Atom_List * li = ri.atom_list_for_atom(j);
      if (NULL == li)
        continue;

//    cerr << "Got list on atom " << j << " of reagent " << i << endl;

      int pf;
      atom_number_t pa;

      if (_identify_corresponding_list_in_products(*li, pf, pa))
      {
        _set_atom_map(highest_atom_map_number, i, j, pf, pa);
        rc++;
      }
    }
  }

  return rc;
}

/*int
RXN_File::_map_corresponding_atom_list (int highest_atom_map_number,
                                        const ISIS_Atom_List & als)
{
  int rcp = -1;
  atom_number_t rca = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _np; i++)
  {
    const ISIS_RXN_FILE_Molecule & pi = _product[i];

    int matoms = pi.natoms ();

    for (int j = 0; j < matoms; j++)
    {
      if (pi.atom_map (j) >= 1)    // great, already mapped
        continue;

      const ISIS_Atom_List * lj = pi.atom_list_for_atom (j);
      if (NULL == lj)
        continue;

      if (*lj != als)    // atom lists are different
        continue;

      if (rcp >= 0)     // more than one matching atom list
        return 0;

      rcp = i;
      rca = j;
    }
  }

  if (rcp < 0)    // never found a match
    return 0;

  _product[rcp].set_atom_map (rca, highest_atom_map_number + 1);
  _product_locator[highest_atom_map_number + 1] = rcp;

  return 1;
}*/

int
RXN_File::_identify_corresponding_list_in_products (const ISIS_Atom_List & als,
                                                    int & pf,
                                                    atom_number_t & pa) const
{
  pf = -1;
  pa = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _np; i++)
  {
    const ISIS_RXN_FILE_Molecule & pi = _product[i];

    int matoms = pi.natoms();

    for (int j = 0; j < matoms; j++)
    {
      if (pi.atom_map(j) >= 1)    // great, already mapped
        continue;

      const ISIS_Atom_List * lj = pi.atom_list_for_atom(j);
      if (NULL == lj)
        continue;

      if (*lj != als)    // atom lists are different
        continue;

      if (pf >= 0)     // more than one matching atom list
        return 0;

      pf = i;
      pa = j;
    }
  }

  if (pf < 0)    // never found a match
    return 0;

  return 1;

  return 1;
}

void
RXN_File::_fill_reagent_and_product_locator_arrays ()
{
  for (int i = 0; i < _nr; i++)
  {
    _reagent[i].identify_which_atoms_we_have(_reagent_locator, i, _atom_map_to_reagent_atom_number);
  }

  for (int i = 0; i < _np; i++)
  {
    _product[i].identify_which_atoms_we_have(_product_locator, i, _atom_map_to_product_atom_number);
  }

#ifdef DEBUG_FILL_REAGENT_AND_PRODUCT_LOCATOR_ARRAYS
  const int h = highest_atom_map_number();
  cerr << "RXN_File::_fill_reagent_and_product_locator_arrays h = " << h << endl;
  for (int i = 0; i <= h; ++i)
  {
    cerr << " map " << i << " reagent " << _reagent_locator[i] << " product " << _product_locator[i] << endl;
  }
#endif

  return;
}

int
RXN_File::highest_atom_map_number () const
{
  int rc = _reagent[0].highest_atom_map_number();

  for (int i = 1; i < _nr; ++i)
  {
    const auto t = _reagent[i].highest_atom_map_number();
    if (t > rc)
      rc = t;
  }

  for (int i = 0; i < _np; ++i)   // we only need to check products if mapping is asymmetric
  {
    const auto t = _product[i].highest_atom_map_number();
    if (t > rc)
      rc = t;
  }

  return rc;
}

/*
  We need a suitable dimension for various arrays
*/

#ifdef IMPLEMENT_IF_NEEDED
int
RXN_File:_max_atoms_in_any_reagent() const
{
  int rc = _reagent[0].natoms();
  
  for (int i = 1; i < _nr; ++i)
  {
    if (_reagent[i].natoms() > rc)
      rc = _reagent[i].natoms();
  }

  return rc;
}
#endif

int
RXN_File::_highest_atom_map_number_or_atoms_in_reagents () const
{
  int rc = _reagent[0].highest_atom_map_number();

  int atoms_in_reagents = 0;

  for (int i = 1; i < _nr; ++i)
  {
    const auto t = _reagent[i].highest_atom_map_number();
    if (t > rc)
      rc = t;

    atoms_in_reagents += _reagent[i].natoms();
  }

  for (int i = 0; i < _np; ++i)   // we only need to check products if mapping is asymmetric
  {
    const auto t = _product[i].highest_atom_map_number();
    if (t > rc)
      rc = t;
  }

  if (rc > 0)
    return rc;

  return atoms_in_reagents;
}

static void
transfer_atom_numbers(ISIS_RXN_FILE_Molecule * reagent_or_product,
                      const int n)
{
  for (int i = 0; i < n; ++i)
  {
    const int matoms = reagent_or_product[i].natoms();

    for (int j = 0; j < matoms; ++j)
    {
      reagent_or_product[i].set_atom_map_number(j, reagent_or_product[i].atom_map(j));
    }
  }

  return;
}

int
RXN_File::_setup_reagent_product_locator_arrays()
{
  DELETE_IF_NOT_NULL_ARRAY(_reagent_locator);
  DELETE_IF_NOT_NULL_ARRAY(_product_locator);

  if (0 == _nr && 0 == _np)
  {
    cerr << "RXN_File:_setup_reagent_product_locator_arrays: empty reaction\n";
    return 0;
  }

  int atoms_in_reagents = 0;
  int highest_atom_map = 0;   // we only need to check reagents, so don't call the member function

  for (int i = 0; i < _nr; i++)
  {
    atoms_in_reagents += _reagent[i].natoms();

    const auto h = _reagent[i].highest_atom_map_number();
    if (h > highest_atom_map)
      highest_atom_map = h;
  }

  for (int i = 0; i < _np; ++i)
  {
    const auto x = _product[i].highest_atom_map_number();
    if (x > highest_atom_map)
      highest_atom_map = x;
  }

//cerr << "RXN_File::_setup_reagent_product_locator_arrays: highest atom map number " << highest_atom_map << endl;

  if (0 == highest_atom_map)
    highest_atom_map = atoms_in_reagents;    // one for each atom, should be fine

// When we map unmapped atoms, we want to make sure we have room in the _locator arrays, so allocate them for the maximum size

  int mxdim = highest_atom_map + 1 + atoms_in_reagents;

  _reagent_locator = new_int(mxdim, -1);
  _product_locator = new_int(mxdim, -1);

  _atom_map_to_reagent_atom_number = new_int(mxdim, -1);
  _atom_map_to_product_atom_number = new_int(mxdim, -1);

  _fill_reagent_and_product_locator_arrays();

  return 1;
}

//#define DEBUG_PREPARE_FOR_REACTION_CONSTRUCTION

int
RXN_File::prepare_for_reaction_construction()
{
  _setup_reagent_product_locator_arrays();

//cerr << "ORPHANS: _do_automatic_atom_mapping " << _do_automatic_atom_mapping << " _auto_fix_orphans " << _auto_fix_orphans << " orphans " << _orphan_atoms.natoms() << endl;

  if (_orphan_atoms.natoms() > 0)    // has already been run
    ;
  else if (! _do_automatic_atom_mapping)
  {
    if (_auto_fix_orphans && ! check_for_widows_and_orphans())   // hard to know if we should let orphans be discovered, just so the user knows they exist
    {
      cerr << "RXN_File::create_reaction:failed widow/orphan check\n";
      return 0;
    }
  }

  for (int i = 0; i < _nr; i++)
  {
    _reagent[i].identify_symmetry_classes();
  }

  for (int i = 0; i < _np; i++)
  {
    _product[i].identify_symmetry_classes();
//  for (int j = 0; j < _product[i].natoms(); ++j)
//  {
//    cerr << " product " << i << " atom " << j << " symm " << _product[i].use_atom_as_symmetry_representative(j) << " type " << _product[i].smarts_equivalent_for_atom(j) << endl;
//  }
  }

#ifdef DEBUG_PREPARE_FOR_REACTION_CONSTRUCTION
  cerr << "RXN_File::prepare_for_reaction_construction:mapping? " << _do_automatic_atom_mapping << endl;
#endif

  if (_do_automatic_atom_mapping)
  {
    int highest_atom_map = highest_atom_map_number();

    int h = highest_atom_map;

    while (_do_automatic_atom_mapping)
    {
      _establish_atom_mapping(highest_atom_map);

//    cerr << "After mapping unmapped atoms, highest atom map number " << highest_atom_map << endl;

      _map_unmapped_atoms(highest_atom_map);

      if (highest_atom_map == h)
        break;

      h = highest_atom_map;
    }

    report_mapping_result(_reagent, _nr, cerr);
    report_mapping_result(_product, _np, cerr);
  
    transfer_atom_numbers(_reagent, _nr);
    transfer_atom_numbers(_product, _np);

    _fill_reagent_and_product_locator_arrays();
  }

#ifdef DEBUG_PREPARE_FOR_REACTION_CONSTRUCTION
  cerr << "RXN_File::prepare_for_reaction_construction final\n";
  Reaction_Smiles_Options rxnsmo;
  write_rxn_smiles(rxnsmo, cerr);

  for (int i = 0; i < _nr; ++i)
  {
    const ISIS_RXN_FILE_Molecule & r = _reagent[i];

    for (int j = 0; j < r.natoms(); ++j)
    {
      cerr << " reagent " << i << " atom " << j << ' ' << r.smarts_equivalent_for_atom(j) << " map " << r.atom_map()[j] << endl;
    }
  }
  debug_print(cerr);
#endif

  return 1;
}

/*
*/

int
RXN_File::create_reaction(IWReaction & rxn,
                          const RXN_File_Create_Reaction_Options & rxnfcro,
                          const int * include_these_atom_map_numbers)
{
  Molecule_to_Query_Specifications mqs;
  mqs.set_make_embedding(1);

  return create_reaction(rxn, rxnfcro, mqs, include_these_atom_map_numbers);
}

int
RXN_File::create_reaction(IWReaction & rxn,
                          const RXN_File_Create_Reaction_Options & rxnfcro,
                          Molecule_to_Query_Specifications & mqs,
                          const int * include_these_atom_map_numbers)
{
//if (NULL == _reagent_locator || NULL == _product_locator)
  {
    if (! prepare_for_reaction_construction())
      return 0;
  }

  if (0 == _nr && 0 == _np)
  {
    cerr << "RXN_File::create_reaction: empty reaction\n";
    return 0;
  }

  const int highest_atom_map = highest_atom_map_number();
  if (highest_atom_map <= 0)
  {
    cerr << "RXN_File::create_reaction:highest atom map not set, cannot continue\n";
    return 0;
  }

  _btype_dim = highest_atom_map + 1;

//cerr << "Dimensionality of bond type arrays " << _btype_dim << endl;

  int t = _btype_dim * _btype_dim;

  DELETE_IF_NOT_NULL_ARRAY(_initial_bond_type);
  DELETE_IF_NOT_NULL_ARRAY(_final_bond_type);

  _initial_bond_type = new bond_type_t[t];
  _final_bond_type   = new bond_type_t[t];

  set_vector(_initial_bond_type, t, static_cast<bond_type_t>(INVALID_BOND_TYPE));
  set_vector(_final_bond_type, t, static_cast<bond_type_t>(INVALID_BOND_TYPE));

//cerr << "RXN_File::create_reaction:input has " << _nr << " reagents and " << _np << " products\n";

  for (int i = 0; i < _nr; i++)
  {
    if (! _reagent[i].fill_bonds_between_mapped_atoms_array(_initial_bond_type, _btype_dim))
    {
      cerr << "RXN_File::create_reaction:fill_bonds_between_mapped_atoms_array failed for reagent " << i << " " << _comment << endl;
      return 0;
    }
  }

  for (int i = 0; i < _np; i++)
  {
    _product[i].fill_bonds_between_mapped_atoms_array(_final_bond_type, _btype_dim);
  }

  _orphan_atoms.fill_bonds_between_mapped_atoms_array(_final_bond_type, _btype_dim);

#ifdef ECHO_BONDING
  for (int i = 0; i < _btype_dim; i++)
  {
    for (int j = 0; j < _btype_dim; j++)
    {
      cerr << " i = " << i << " j = " << j << " ini " << _initial_bond_type[i * _btype_dim + j] << " final " << _final_bond_type[i + _btype_dim + j] << endl;

      assert (_initial_bond_type[_btype_dim * i + j] == _initial_bond_type[_btype_dim * j + i]);
      assert (_final_bond_type[_btype_dim * i + j] == _final_bond_type[_btype_dim * j + i]);
    }
  }
#endif

  int reagents_mapped = molecules_with_at_least_one_mapped_atom(_reagent, _nr);

  if (_orphan_atoms.natoms() > 0)
    reagents_mapped++;

  if (reagents_mapped < _nr)
  {
    cerr << "RXN_File::create_reaction: only " << reagents_mapped << " of " << _nr << " reagents contain mapped atoms, cannot continue\n";
    return 0;
  }
   
  if (_nr > 1)     // wrong if someone is only doing changes in multiple scaffold embeddings
  {
    int changing_reagents = molecules_with_at_least_one_inter_fragment_change(_reagent, _nr);

    if (changing_reagents != (_nr - 1))
    {
      cerr << "RXN_File::create_reaction: only " << reagents_mapped << " of " << _nr << " reagents undergo inter-particle changes, cannot continue\n";
      return 0;
    }
  }

  if (! _create_reaction(rxn, rxnfcro, mqs, highest_atom_map, include_these_atom_map_numbers))
    return 0;

  if (_fname_for_echo.length())
  {
    std::ofstream echo_rxn(_fname_for_echo.null_terminated_chars(), std::ios::out);
    if (! echo_rxn.good())
      cerr << "RXN_File::create_reaction:cannot open echo file '" << _fname_for_echo << "'\n";
    else
      do_write(echo_rxn);
  }

  return 1;
}

/*

*/

int
ISIS_RXN_FILE_Molecule::create_query (Reaction_Site & r,
                                      const int * include_these_atoms,
																			std::ofstream *queryOutStream)
{
  Molecule_to_Query_Specifications mqs;

  mqs.set_make_embedding(1);
  if (! _interpret_atom_alias_as_smarts)
    mqs.set_interpret_atom_alias_as_smarts(0);

  mqs.set_set_element_hits_needed_during_molecule_to_query(0);    // arbitrary choice in Jun 2016. We are doing a lot of small radius things, and the element count is just slowing things down...

  return create_query(r, include_these_atoms, mqs, queryOutStream);
}

int
ISIS_RXN_FILE_Molecule::create_query (Reaction_Site & r,
                                      const int * include_these_atoms,
                                      Molecule_to_Query_Specifications & mqs,
																			std::ofstream *queryOutStream)
{
//cerr << "ISIS_RXN_FILE_Molecule::create_query: input " << smiles() << endl;

  if (! r.create_from_molecule(*this, mqs, include_these_atoms))
  {
    cerr << "ISIS_RXN_FILE_Molecule::create_from_molecule: cannot create query from '" << smiles() << "'\n";
    return 0;
  }
  r.set_comment(this->name());
 

  if (0 == r[0]->root_atoms())
  {
    cerr << "ISIS_RXN_FILE_Molecule::create_from_molecule:no root atoms, impossible\n";
    return 0;
  }
  
  
  
  if (queryOutStream != NULL && queryOutStream->is_open())
  	r.write_msi (*queryOutStream);
  	
//  std::ofstream queryOut ("queryout.msi", std::ofstream::out);
//  r.write_msi (queryOut);
//  queryOut.close();
//  
//  iwstring_data_source input("queryout.msi");
//
//  if (! input.good())
//  {
//    cerr << "Cannot open reaction 'queryout.msi'\n";
//    return 0;
//  }
//
//  msi_object msi;
//  msi.set_display_no_data_error_message(0);
//  Sidechain_Match_Conditions  smc;
//  
//  int thisIndex = 1;
//  while(msi.read(input))
//  {
//  	Single_Substructure_Query thisQuery;
//  	thisQuery.construct_from_msi_object(msi);
//  	
//  	std::ostringstream thisFilename;
//
//  	thisFilename << "queryout" << thisIndex << ".msi";
//
//  	std::ofstream queryOut (thisFilename.str().c_str(),  std::ofstream::out);
//
//  	thisQuery.write_msi (queryOut);
//  	queryOut.close();
//
//		thisIndex++;
//  }
//  
  

  return 1;
}

int
RXN_File::_create_query (Reaction_Site & r,
                         ISIS_RXN_FILE_Molecule & m,
                         Molecule_to_Query_Specifications & mqs,
                         const int * include_these_atoms)
{
	m.set_name(name());
  return m.create_query(r, include_these_atoms, mqs, _queryOutStream);
}

static int
first_positive_item (const int * a, int n)
{
  for (int i = 0; i < n; i++)
  {
    if (a[i] > 0)
      return a[i];
  }

  return -1;
}

//#define DEBUG_CREATE_REACTION

int
RXN_File::_create_reaction(IWReaction & rxn,
                           const RXN_File_Create_Reaction_Options & rxnfcro,
                           Molecule_to_Query_Specifications & mqs,
                           int highest_atom_map,
                           const int * include_these_atom_map_numbers)
{
#ifdef DEBUG_CREATE_REACTION
  cerr << "RXN_File::_create_reaction\n";
  write_mapped_reaction(cerr);

  for (int i = 1; i <= highest_atom_map; ++i)
  {
    cerr << i << " include_these_atom_map_numbers " << include_these_atom_map_numbers[i] << endl;
  }
#endif

  assert (_nr > 0);

//if (_nr > 1)
//  cerr << "RXN_File::_create_reaction:warning reaction " << _comment << " contains " << _nr << " reagents, only first one processed. Contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";

  rxn.set_comment(_comment);

  int x = _nr - 1;

  if (_nr > 1 && rxnfcro.only_create_query_from_first_reagent())   // second sidechain will contain all fixed reagents
    x = 1;

  if (_orphan_atoms.natoms() > 0)
    x++;

#ifdef DEBUG_CREATE_REACTION
  cerr << "NR " << _nr << " _orphan_atoms " << _orphan_atoms.natoms() << " sidechains in rxn " << x << endl;
#endif

// Note that this can create lots of problems if there are important mapped atoms in the other sidechains...

  if (! rxn.set_number_sidechains(x))
  {
    cerr << "RXN_File::_create_reaction: cannot resize reaction for " << x << " reagents\n";
    return 0;
  }

//int is_subset = 0;
//if (NULL == include_these_atom_map_numbers)
//  ;
//else if (std::count_if(include_these_atom_map_numbers, include_these_atom_map_numbers + highest_atom_map + 1, [] (const int i) { return i > 0;}))
//  is_subset = 1;

// If we are dealing with a subset, we need to convert this subset specified by
// atom map numbers, to atom numbers in each of the reagents.

  Reaction_Subset subset(_nr);

  if (! subset.build(_reagent, _reagent_locator, include_these_atom_map_numbers, highest_atom_map))
  {
    cerr << "RXN_File::_create_reaction:cannot create subset\n";
    return 0;
  }

//cerr << "LINE " << __LINE__ << endl;
//subset.debug_print(cerr);
  
#ifdef DEBUG_CREATE_REACTION
  cerr << "Set " << x << " sidechains in RXN\n";
#endif

  _identify_square_bonding_changes(highest_atom_map);

  resizable_array_p<Bond> connections_involved_in_substitutions;    // both makes and breaks
  resizable_array<Replace_Atom *> atoms_to_replace;

//_identify_atom_substitutions(highest_atom_map, atoms_to_replace, connections_involved_in_substitutions);

  resizable_array_p<Bond> bonds_to_be_broken;

  _identify_bonds_to_be_broken(highest_atom_map, bonds_to_be_broken, connections_involved_in_substitutions, include_these_atom_map_numbers);

  resizable_array_p<Bond> bonds_to_be_made;

  _identify_bonds_to_be_made(highest_atom_map, bonds_to_be_made, connections_involved_in_substitutions, include_these_atom_map_numbers);

#ifdef DEBUG_CREATE_REACTION
  for (int i = 0; i < bonds_to_be_made.number_elements(); ++i)
  {
    cerr << " bond 2 be made " << *(bonds_to_be_made[i]) << endl;
  }
#endif

  resizable_array_p<Reaction_Place_Isotope> isotopes;

  _identify_isotopes_to_be_placed(highest_atom_map, isotopes);

//cerr << "Identified " << bonds_to_be_made.number_elements() << " bonds to be made\n";

  if (! _interpret_atom_alias_as_smarts)
  {
    for (int i = 0; i < _nr; ++i)
    {
      _reagent[i].set_interpret_atom_alias_as_smarts(0);
    }
  }

// Now this gets ugly. We need to create an atom map for each reagent fragment.

  if (! _create_query(rxn, _reagent[0], mqs, subset.include_atom(0)))
  {
    cerr << "RXN_File::create_reaction:cannot create query from reagent 0\n";
    return 0;
  }
  

  _identify_kekule_forms_to_be_toggled(include_these_atom_map_numbers);

  _reagent[0].add_toggle_kekule_forms(rxn, subset, 0, rxnfcro);

  rxn.set_find_kekule_forms_for_bad_valence(1);

//cerr << " LINE " << __LINE__ << " nr " << _nr << " orphan " << _orphan_atoms.natoms() << endl;
  if (1 == _nr)    // nothing to worry about
    ;
  else if (rxnfcro.only_create_query_from_first_reagent())
  {
    if (2 == _nr)
      rxn.sidechain(0)->set_single_reagent(_reagent[1]);
    else
    {
      Molecule single_reagent;
      for (int i = 1; i < _nr; ++i)
      {
        const Molecule & ri = _reagent[i];
        single_reagent.add_molecule(&ri);
      }
      rxn.sidechain(0)->set_single_reagent(single_reagent);
    }
  }
  else
  {
    for (int i = 1; i < _nr; i++)
    {
      Sidechain_Reaction_Site * r = rxn.sidechain(i - 1);
  
      if (_reagent[i].single_reagent_only())
        r->set_single_reagent(_reagent[i]);
      else
      {
        if (! _create_query(*r, _reagent[i], mqs, subset.include_atom(i)))
        {
          cerr << "Cannot create query for reagent " << i << " " << _reagent[i].smiles() << endl;
          return 0;
        }
      }

      _reagent[i].add_toggle_kekule_forms(*r, subset, i, rxnfcro);
    }
  }

  if (_orphan_atoms.natoms() > 0)   // do not create a query for these, but add it as a reagent
  {
    Sidechain_Reaction_Site * r = rxn.sidechain(x-1);
    r->set_single_reagent(_orphan_atoms);
  }

// Now that the queries are in place, make them unique matches only

// Complication with atom replacements. Some naturally belong to the scaffold, but convention
// dictates that they belong to a sidechain

//cerr << "Found " << atoms_to_replace.number_elements() << " replace atoms\n";
  for (int i = 0; i < atoms_to_replace.number_elements(); i++)
  {
    Replace_Atom * r = atoms_to_replace[i];

    const int c1 = r->a1().in_component();
    const int c2 = r->a2().in_component();

//  We need to convert things to the numbering used in the reaction object

    if (0 == c1)
      r->a1().set_in_scaffold();
    else
      r->a1().set_in_component(c1 - 1);

    if (0 == c2)
      r->a2().set_in_scaffold();
    else
      r->a2().set_in_component(c2 - 1);

    if (0 == c1 && 0 == c2)     // intra-scaffold replacement
      rxn.add_replace_atom(r);
    else if (0 == c1)
      rxn.sidechain(c2 - 1)->add_replace_atom(r);
    else if (0 == c2)
      rxn.sidechain(c1 - 1)->add_replace_atom(r);
    else
      rxn.sidechain(c2 - 1)->add_replace_atom(r);
  }

  int nb = bonds_to_be_broken.number_elements();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = bonds_to_be_broken[i];

    const atom_number_t m1 = b->a1();
    const atom_number_t m2 = b->a2();

//  cerr << "Bond between mapped atoms " << m1 << " and " << m2 << " being broken\n";

    int r = _reagent_locator[m1];

    int q1 = _reagent[r].which_is_mapped_atom(m1);
    int q2 = _reagent[r].which_is_mapped_atom(m2);

//  cerr << " which corresponds to mapped atoms " << q1 << ' ' << _reagent[r].smarts_equivalent_for_atom(q1) << " and " << q2 << ' ' << _reagent[r].smarts_equivalent_for_atom(q2) << endl;
//  subset.debug_print(cerr);

    assert (q1 >= 0 && q2 >= 0);

    if (r == 0 || !rxnfcro.only_create_query_from_first_reagent())
    {
			q1 = subset.atom_number_in_subset(q1, r);
			q2 = subset.atom_number_in_subset(q2, r);
    }

    assert (q1 >= 0 && q2 >= 0);

    int rc;

    if (0 == r)
    {
      rc = rxn.add_bond_to_be_broken(q1, q2);
    }
    else if (r <= rxn.number_sidechains())
    {
      Sidechain_Reaction_Site * s = rxn.sidechain(r - 1);

      rc = s->add_bond_to_be_broken(q1, q2);
    }

    if (0 == rc)
    {
      cerr << "RXN_File::_create_reaction: cannot add bond " << (*b) << " in fragment " << r << endl;
      return 0;
    }
  }

  nb = bonds_to_be_made.number_elements();

#ifdef DEBUG_CREATE_REACTION
  cerr << "RXN_File::_create_reaction:processing " << nb << " bonds to be made\n";
#endif

  for (int i = 0; i < nb; i++)   // first process bonds within the same fragment
  {
    const Bond * b = bonds_to_be_made[i];

    int m1 = b->a1();
    int m2 = b->a2();

    int r1 = _reagent_locator[m1];
    int r2 = _reagent_locator[m2];
    
    if (r1 != r2)    // this loop processes bonds within the same fragment
      continue;

    if (r1 < 0 || r2 < 0)    // orphan
      continue;

//#define DEBUG_BONDS_TO_BE_MADE
#ifdef DEBUG_BONDS_TO_BE_MADE
    cerr << "Making ";
    if (b->is_single_bond())
      cerr << "single";
    else if (b->is_double_bond())
      cerr << "double";
    else if (b->is_triple_bond())
      cerr << "triple";
    cerr << " bond between mapped atoms " << m1 << " and " << m2 << " in reagent " << r1 << endl;
#endif

    int q1 = _reagent[r1].which_is_mapped_atom(m1);
    int q2 = _reagent[r1].which_is_mapped_atom(m2);    	

    assert (q1 >= 0 && q2 >= 0);

#ifdef DEBUG_BONDS_TO_BE_MADE
    cerr << "which are atoms " << q1 << " and " << q2 << ' ' << _reagent[r1].smarts_equivalent_for_atom(q1) << ' ' << _reagent[r1].smarts_equivalent_for_atom(q2) << endl;
#endif

    if (r2 == 0 || !rxnfcro.only_create_query_from_first_reagent())
    {
	    q1 = subset.atom_number_in_subset(q1, r1);
	    q2 = subset.atom_number_in_subset(q2, r2);
	  }
	  
    if (q1 < 0 || q2 < 0) 
    {
      cerr << "Matched atom " << _reagent[r1].which_is_mapped_atom(m1) << " in subset " << q1 << " or " << _reagent[r1].which_is_mapped_atom(m2) << " in subset " << q2 << " invalid\n";
      Molecule mcopy(_reagent[r1]);
      write_isotopically_labelled_smiles(mcopy, cerr);
      subset.debug_print(cerr);
      abort();
    }

#ifdef DEBUG_BONDS_TO_BE_MADE
    cerr << "In subset " << q1 << " and " << q2 << endl;
#endif

    assert (q1 >= 0 && q2 >= 0);

    int rc;

    if (0 == r1)
      rc = rxn.add_bond_to_be_made(q1, q2, b->btype());
    else if ((r1 - 1) >= rxn.number_sidechains())
    {
      cerr << "RXN_File::_create_reaction:reagent " << r1 << " out of range, ignored\n";   // perhaps only creating reaction from first reagent, but other reagents have changing atoms.
      continue;
    }
    else
      rc = rxn.sidechain(r1-1)->add_bond_to_be_made(q1, q2, b->btype());

    if (0 == rc)
    {
      cerr << "RXN_File::_create_reaction: cannot add breaking bond " << (*b) << " in fragment " << r1 << endl;
      return 0;
    }
  }

// Another loop for inter particle bonds

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = bonds_to_be_made[i];

    int m1 = b->a1();
    int m2 = b->a2();

    int r1 = _reagent_locator[m1];
    int r2 = _reagent_locator[m2];

    if (r1 == r2)    // this loop processes bonds between fragments
      continue;

    if (r1 < 0 || r2 < 0)
      continue;

    if (r1 > r2)
    {
      std::swap(r1,r2);
      std::swap(m1,m2);
    }

//  cerr << "Inter particle bond between mapped atoms " << m1 << " (" << r1 << ") and " << m2 << " (" << r2 << ") being made\n";

    const ISIS_RXN_FILE_Molecule & mol1 = _reagent[r1];
    int q1 = mol1.which_is_mapped_atom(m1);

    const ISIS_RXN_FILE_Molecule & mol2 = _reagent[r2];
    int q2 = mol2.which_is_mapped_atom(m2);

//  cerr << " atoms " << q1 << " and " << q2 << endl;
    if (r1 == 0 || ! rxnfcro.only_create_query_from_first_reagent())
    {
    	q1 = subset.atom_number_in_subset(q1, r1);
    }
    if (r2 == 0 || ! rxnfcro.only_create_query_from_first_reagent())
    { 	
    	q2 = subset.atom_number_in_subset(q2, r2);
    }

    if (q1 < 0 || q2 < 0)
      subset.debug_print(cerr);

    assert (q1 >= 0 && q2 >= 0);
//  cerr << "Atoms " << q1 << " and " << q2 << endl;

    if (r2 > rxn.number_sidechains())
      continue;

    Sidechain_Reaction_Site * s = rxn.sidechain(r2 - 1);

    if (! s->add_inter_particle_bond(r1 - 1, q1, q2, b->btype()))   // In the reaction object the scaffold is designated -1
    {
      cerr << "RXN_File::_create_reaction: cannot add inter particle bond " << (*b) << " between fragments " << r1 << " and " << r2 << endl;
      return 0;
    }
  }

// Deal with orphans

//cerr << "Orphan is " << _orphan_atoms.smiles() << endl;

  if (_orphan_atoms.natoms())
  {
    Sidechain_Reaction_Site * s = rxn.sidechain(_nr - 1);   // the extra sidechain for the orphans

    for (int i = 0; i < nb; ++i)
    {
      const Bond * b = bonds_to_be_made[i];

      int m1 = b->a1();
      int m2 = b->a2();

      int r1 = _reagent_locator[m1];
      int r2 = _reagent_locator[m2];

      if (r1 == r2)    // orphans are between reagents
        continue;

      if (r1 >= 0 && r2 >= 0)   // non-orphan case
        continue;

      if (r1 < 0 && r2 < 0)    // should not happen
      {
        cerr << "RXN_File::_create_reaction:bond to be made " << m1 << " (" << r1 << ") to " << m2 << " (" << r2 << ") invalid reagents\n";
        return 0;
      }

//    make sure m2 is the orphan

      if (r1 < 0)
      {
        std::swap(r1,r2);
        std::swap(m1,m2);
      }

      r2 = _product_locator[m2];

//    cerr << "Inter particle orphan bond between mapped atoms " << m1 << " (" << r1 << ") and " << m2 << " (" << r2 << ") being made\n";

//    So now we have M1 as an atom on the LHS and M2 is an orphan on the RHS

      const ISIS_RXN_FILE_Molecule & mol1 = _reagent[r1];
      int q1 = mol1.which_is_mapped_atom(m1);

      int q2 = _orphan_atoms.which_is_mapped_atom(m2);

      if (r1 == 0 || ! rxnfcro.only_create_query_from_first_reagent())
      {
      	q1 = subset.atom_number_in_subset(q1, r1);    // note that we do NOT translate q2
			}
//    cerr << "Final addition to " << (r1-1) << " btw " << q1 << " and " << q2 << endl;
      s->add_inter_particle_bond(r1 - 1, q1, q2, b->btype());
    }
  }

  if (_remove_product_fragments && _np > 1)
  {
    for (int i = 1; i < _np; i++)
    {
      const ISIS_RXN_FILE_Molecule & mp = _product[i];

      const int * amap = mp.atom_map();

      int mpd1 = first_positive_item(amap, mp.natoms());
      if (mpd1 <= 0)    // strange, all atoms mapped to 0 in that fragment
        continue;

      int r = _reagent_locator[mpd1];     // which reagent has mapped atom MPD1

      const ISIS_RXN_FILE_Molecule & mr = _reagent[r];

      int mpd2 = mr.which_is_mapped_atom(mpd1);    // which atom number is it in the query

      if (r == 0 || ! rxnfcro.only_create_query_from_first_reagent())
      {
				mpd2 = subset.atom_number_in_subset(mpd2, r);
			}
      assert (mpd2 >= 0);

//    cerr << "Reagent " << r << " will remove the fragment containing matched atom " << mpd2 << endl;

      int rc;
      if (0 == r)
      {
        rc = rxn.add_fragment_to_be_removed(mpd2);
      }
      else
      {
        Sidechain_Reaction_Site * s = rxn.sidechain(r - 1);

        rc = s->add_fragment_to_be_removed(mpd2);
      }

      if (0 == rc)
        return 0;
    }
  }

// Do any of our reagents have stereo centres to be inverted

  for (int i = 0; i < _nr; i++)
  {
    const ISIS_RXN_FILE_Molecule & m = _reagent[i];

    if (0 == i)
      m.add_chiral_centres_to_be_inverted(rxn);
    else if (i <= rxn.number_sidechains())
    {
      Sidechain_Reaction_Site * s = rxn.sidechain(i - 1);

//    need to adjust for subsetting

      m.add_chiral_centres_to_be_inverted(*s);
    }
  }

// Are there any stereo centres formed in the products

  _look_for_stereo_centres_made(rxn);

// Are there any element transformations

  for (int i = 1; i <= highest_atom_map; i++)
  {
    int r = _reagent_locator[i];
    int p = _product_locator[i];

    if (r < 0 || p < 0)
      continue;

    int ra = _reagent[r].which_is_mapped_atom(i);
    int pa = _product[p].which_is_mapped_atom(i);

    const Element * re = _reagent[r].elementi(ra);
    const Element * pe = _product[p].elementi(pa);

    if (re->atomic_number() == pe->atomic_number())   // no change
      continue;

    if (! pe->is_in_periodic_table())    // not sure what to do with these
      continue;

    if (r == 0 || ! rxnfcro.only_create_query_from_first_reagent())
    {
			ra = subset.atom_number_in_subset(ra, r);
		}

//  this is kind of difficult, what to do if the elements are different.
//  Many times, if there was a list on the LHS, they will just draw a carbon in the RHS.

    const MDL_Atom_Data * rad = _reagent[r].mdl_atom_data(ra);

//  cerr << "Initial atomis symbol() " << rad->initial_atomic_symbol() << endl;

    if ((_reagent[r].is_an_atom_list(ra) || 'Q' == rad->initial_atomic_symbol() || 'A' == rad->initial_atomic_symbol()) &&
         6 == pe->atomic_number() && _interpret_carbon_on_rhs_of_list_as_no_change)
      continue;

//  cerr << "Different atomic numbers " << re->atomic_number() << " and " << pe->atomic_number() << " ipt " << re->is_in_periodic_table() << " " << pe->is_in_periodic_table() << endl;
//  cerr << _reagent[r].smiles() << " and " << _product[p].smiles() << endl;
     
    Reaction_Change_Element * tmp = new Reaction_Change_Element;
    tmp->set_atom(ra);
    tmp->set_element(pe);

    if (0 == r)
    {
      rxn.add_element_to_be_changed(tmp);
    }
    else
    {
      Sidechain_Reaction_Site * s = rxn.sidechain(r - 1);
      s->add_element_to_be_changed(tmp);
    }
  }

// Any formal charges to be set

  for (int i = 1; i <= highest_atom_map; i++)
  {
    const int r = _reagent_locator[i];
    const int p = _product_locator[i];

    if (r < 0 || p < 0)
      continue;

    int ra = _reagent[r].which_is_mapped_atom(i);
    int pa = _product[p].which_is_mapped_atom(i);
    assert(ra >= 0 && pa >= 0);

    const int fcr = _reagent[r].formal_charge(ra);
    const int fcp = _product[p].formal_charge(pa);

    if (fcr == fcp)    // no change
      continue;

	  if (r == 0 || ! rxnfcro.only_create_query_from_first_reagent())
    {
			ra = subset.atom_number_in_subset(ra, r);    // re-use the variable RA
		}
		
    if (ra < 0)
    {
      cerr << "formal_charge:atom map " << i << " in reagent " << r << " atom " << ra << " product " << p << " atom " << pa << " charges " << fcr << ' ' << fcp << " not in subset\n";
      return 0;
    }

    if (0 == r)
      rxn.add_formal_charge_to_assign(ra, fcp);
    else
    {
      Sidechain_Reaction_Site * s = rxn.sidechain(r - 1);
      s->add_formal_charge_to_assign(ra, fcp);
    }
  }

// How about isotopes

  for (int i = 0; i < isotopes.number_elements(); i++)
  {
    const Reaction_Place_Isotope * rpi = isotopes[i];

    int m = rpi->atom();

    int r = _reagent_locator[m];

    int a = _reagent[r].which_is_mapped_atom(m);

	  if (r == 0 || ! rxnfcro.only_create_query_from_first_reagent())
    {
			a = subset.atom_number_in_subset(a, r);
		}
		
    assert (a >= 0);

    if (0 == r)
      rxn.add_isotope_to_be_placed(a, rpi->isotope());
    else
    {
      Sidechain_Reaction_Site * s = rxn.sidechain(r - 1);
      s->add_isotope_to_be_placed(a, rpi->isotope());
    }
  }

// Are there any elements that disappear

  if (_remove_unmapped_atoms_that_disappear)
    _look_for_unmapped_atoms_that_disappear(highest_atom_map, rxn, subset,rxnfcro);

  if (_unconnect_unmapped_atoms_that_exceed_product_valence)
    _look_for_unmapped_atoms_that_exceed_product_valence(highest_atom_map, rxn);

//cerr << "LINE " << __LINE__ << endl;
//rxn.debug_print(cerr);

  if (! rxn.check_internal_consistency())
  {
    cerr << "RXN_File::_create_reaction: reaction internal consistency failed\n";
    rxn.debug_print(cerr);

    return 0;
  }

  return 1;
}

static int
gather_unmapped_elements (int n,
                          const ISIS_RXN_FILE_Molecule * molecules,
                          extending_resizable_array<int> & count)
{
  if (0 == count.elements_allocated())
    count.extend(1324);     // a likely value of the largest symbol hash

  for (int i = 0; i < n; i++)
  {
    const ISIS_RXN_FILE_Molecule & m = molecules[i];

    int matoms = m.natoms();

    const int * atom_map = m.atom_map();

    for (int j = 0; j < matoms; j++)
    {
      if (atom_map[j] >= 1)
        continue;

      if (m.is_an_atom_list(j))
        continue;

      const Element * ej = m.elementi(j);

      int h = ej->atomic_symbol_hash_value();

//    cerr << "Atom " << j << " type " << ej->symbol() << " ashv " << h << endl;

      count[h]++;
    }
  }

  return 1;
}

/*
  Do we have any unmapped atoms that disappear - more present in the reagents than present
  in the products
*/

int
RXN_File::_look_for_unmapped_atoms_that_disappear (int & highest_atom_map, 
                                                   IWReaction & rxn,
                                                   const Reaction_Subset & subset,
                                                   const RXN_File_Create_Reaction_Options & rxnfcro)
{
  extending_resizable_array<int> in_reagents;

  gather_unmapped_elements(_nr, _reagent, in_reagents);

  extending_resizable_array<int> in_products;

  gather_unmapped_elements(_np, _product, in_products);
  
// For simplicity, size the two arrays the same

  if (in_reagents.number_elements() > in_products.number_elements())
    in_products.extend(in_reagents.number_elements());
  else if (in_products.number_elements() > in_reagents.number_elements())
    in_reagents.extend(in_products.number_elements());

  int n = in_reagents.number_elements();

  for (int i = 0; i < n; i++)
  {
    if (in_reagents[i] == in_products[i])    // may be 0 occurrences
      continue;

    if (in_reagents[i] < in_products[i])
    {
      cerr << "RXN_File::_look_for_unmapped_atoms_that_disappear: strange, element with hash " << i << " created from nothing\n";
      continue;
    }

//  cerr << "Element hash " << i << " in reagents " << in_reagents[i] << " but " << in_products[i] << " in the products\n";

    if (0 == in_products[i])
      _add_atom_removals(rxn, i);
  }

// RXN files from NextMove canhave mapped atoms that do not appear on the RHS

  for (int i = 1; i <= highest_atom_map; ++i)
  {
    if (_reagent_locator[i] >= 0 && _product_locator[i] < 0)
    {
      const int r = _reagent_locator[i];
      int a = _reagent[r].which_is_mapped_atom(i);
      assert(a >= 0);

      if (r == 0 || ! rxnfcro.only_create_query_from_first_reagent())
      {
      	a = subset.atom_number_in_subset(a, r);
      }

      if (0 == r)
        _add_atom_removal(rxn, a);
      else
      {
        Sidechain_Reaction_Site * s = rxn.sidechain(r - 1);
        _add_atom_removal(*s, a);
      }
    }
  }

  return 1;
}

int
RXN_File::_look_for_unmapped_atoms_that_exceed_product_valence (int highest_atom_map_number,
                                IWReaction & rxn) const
{
  for (int i = 1; i <= highest_atom_map_number; i++)
  {
    int ri = _reagent_locator[i];
    if (ri < 0)
      continue;

    ISIS_RXN_FILE_Molecule & r = _reagent[ri];

    atom_number_t unmapped_neighbour_reagent;
    int notused;
    if (! r.identify_unmapped_neighbour(i, unmapped_neighbour_reagent, notused))
      continue;

    int pi = _product_locator[i];

    ISIS_RXN_FILE_Molecule & p = _product[pi];

    atom_number_t unmapped_neighbour_product;
    if (p.identify_unmapped_neighbour(i, unmapped_neighbour_product, notused))   // has an unmapped neighbour, can't do anything
      continue;

    atom_number_t a1 = r.which_is_mapped_atom(i);

    if (0 == ri)
      rxn.add_bond_to_be_broken (a1, unmapped_neighbour_reagent);
    else
    {
      Sidechain_Reaction_Site * s = rxn.sidechain (ri - 1);

      s->add_bond_to_be_broken(a1, unmapped_neighbour_reagent);
    }
  }

  return 1;
}

/*
  We know that all elements of type H disappear. 
*/

int
RXN_File::_add_atom_removals (IWReaction & rxn,
                              int h) const
{
  for (int i = 0; i < _nr; i++)
  {
    const ISIS_RXN_FILE_Molecule & m = _reagent[i];

    int matoms = m.natoms();

    const int * atom_map = m.atom_map();

    for (int j = 0; j < matoms; j++)
    {
      if (atom_map[j] >= 1)
        continue;

      if (m.is_an_atom_list(j))
        continue;

      const Element * e = m.elementi(j);

      if (h != e->atomic_symbol_hash_value())
        continue;

//    cerr << "Atom " << j << " in reagent " << i << " will be removed, map " << atom_map[j] << endl;

//    Need to remove the atom.

      if (0 == i)
        _add_atom_removal(rxn, j);
      else
      {
        Sidechain_Reaction_Site * s = rxn.sidechain(i - 1);
        _add_atom_removal(*s, j);
      }
    }
  }

  return 1;
}

int
RXN_File::_add_atom_removal (Reaction_Site & r,
                             int a) const
{
  Substructure_Atom * q = r.query_atom_with_initial_atom_number(a);

  if (NULL == q)
  {
    cerr << "RXN_File::_add_atom_removal:yipes, no query atom corresponds to atom number " << a << endl;
    return 0;
  }

//cerr << "Initial atom number " << a << " query_unique_id " << q->unique_id() << endl;

// 23 Sept 2010. Not sure what is going on here. Things were breaking when I used
// the unique_id, but worked OK when I jused the atom number. But that should not
// work. Need to dig into this some more.

//r.add_atom_to_be_removed(q->unique_id());

  r.add_atom_to_be_removed(a);

  return 1;
}

/*
  We need to check whether or not a bond between two atoms is in a pre-existing list of bonds
*/

static int
bond_already_present (const resizable_array_p<Bond> & bonds,
                      int m1,
                      int m2)
{
  assert(m1 != m2);

  int nb = bonds.number_elements();
  for (int i = 0; i < nb; i++)
  {
    if (bonds[i]->involves(m1, m2))
      return 1;
  }

  return 0;
}

int
RXN_File::_mapped_atoms_bonded_in_products (int m1, int m2) const
{
  int p1 = _product_locator[m1];
  int p2 = _product_locator[m2];

  if (p1 != p2)   // mapped atoms appear in different product fragments
    return 0;

  return _product[p1].mapped_atoms_are_bonded(m1, m2);
}

/*
*/

int
ISIS_RXN_FILE_Molecule::identify_connected_mapped_atoms (atom_number_t zatom,
                                                         resizable_array<int> & connected_to) const
{
  const Atom * a = Molecule::atomi(zatom);

  for (int i = 0; i < a->ncon(); i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (_atom_map[j] < 1)
      continue;

    connected_to.add(_atom_map[j]);    // we return map numbers
  }

  return connected_to.number_elements();
}

/*
  We are looking for atom replacements. For mapped atom M, we have the set of
  mapped atoms initially connected, and the set of mapped atoms ultimately connected.
  All the connections should be single bonds
*/

int
RXN_File::_looks_like_atom_replacement (int m,
                                        resizable_array<int> & initial_connections,
                                        resizable_array<int> & final_connections) const
{
  for (int i = initial_connections.number_elements() - 1; i >= 0; i--)  // first remove all fully preserved connections
  {
    int j = initial_connections[i];

    if (! final_connections.contains(j))
      continue;

    assert (SINGLE_BOND == _initial_bond_type[_btype_dim * m + j] && SINGLE_BOND == _final_bond_type[_btype_dim * m + j]);

//  At this stage, we know that the connection is preserved during the reaction

    initial_connections.remove_item(i);
    final_connections.remove_all(j);
  }

  if (0 == initial_connections.number_elements() || 0 == final_connections.number_elements())
    return 0;

// If there is just one connection changed, then that looks like a replacement

  if (1 != initial_connections.number_elements())
    return 0;

  if (1 != final_connections.number_elements())
    return 0;

  return 1;
}

//#define DEBUG_IDENTIFY_ATOM_SUBSTITUTIONS

/*
  What is a substitution?

  A substitution is a reagent atom that has one, and only one, singly bonded substituent
  that is lost and replaced by another singly bonded atom
*/

int
RXN_File::_identify_atom_substitutions (int highest_atom_map,
                                         resizable_array<Replace_Atom *> & atom_substitutions,
                                         resizable_array_p<Bond> & connections_involved_in_substitutions) const
{
  int ham1 = highest_atom_map + 1;

#ifdef DEBUG_IDENTIFY_ATOM_SUBSTITUTIONS
  cerr << "RXN_File::_identify_atom_substitutions: highest atom map " << ham1 << endl;
#endif

  for (int i = 0; i < _nr; i++)    // loop over all reagents
  {
    const ISIS_RXN_FILE_Molecule & ri = _reagent[i];

    const int * amap = ri.atom_map();

    int matoms = ri.natoms();
    for (int j = 0; j < matoms; j++)     // loop over every atom in each reagent
    {
      if (amap[j] < 1)
        continue;

      if (6 != ri.atomic_number(j))     // We only do replacements adjacent to Carbon atoms
        continue;

      if (ri.nbonds(j) > ri.ncon(j))   // unsaturated, can't have a chiral centre, no atom replacements adjacent to this atom
        continue;

      resizable_array<int> initially_connected_to;
      int icon = ri.identify_connected_mapped_atoms(j, initially_connected_to);

      if (0 == icon)     // no possible substitutions here
        continue;

      int pj = _product_locator[amap[j]];     // in which product is mapped atom J found

      assert(pj >= 0);

      const ISIS_RXN_FILE_Molecule & p = _product[pj];

      atom_number_t k = p.which_is_mapped_atom(amap[j]);

      if (p.nbonds(k) > p.ncon(k))    // unsaturated, can't have a chiral centre, no atom replacements possible
        continue;

      resizable_array<int> ultimately_connected_to;
      int fcon = p.identify_connected_mapped_atoms(k, ultimately_connected_to);

      if (0 == fcon)     // atom probably being eliminated, no substitutions possible
        continue;

#ifdef DEBUG_IDENTIFY_ATOM_SUBSTITUTIONS
      cerr << "atom " << j << " (mapped " << amap[j] << ") has " << initially_connected_to.number_elements() << " initial connections and " << ultimately_connected_to.number_elements() << " final connections\n";
#endif

      if (! _looks_like_atom_replacement(amap[j], initially_connected_to, ultimately_connected_to))
        continue;

      int m1 = initially_connected_to[0];
      int m2 = ultimately_connected_to[0];

      if (_involved_in_square_bond[m1 * ham1 + m2])
        continue;

      if (_mapped_atoms_bonded_in_products(m1, m2))    // Beckmann rearrangement problem
        continue;

#ifdef DEBUG_IDENTIFY_ATOM_SUBSTITUTIONS
      cerr << "From atom " << amap[j] << " we discern replace atom " << m1 << " in reagent " << i << " with atom " << m2 << " in reagent " << _reagent_locator[m2] << endl;
#endif

      Replace_Atom * r = new Replace_Atom;

      Matched_Atom_in_Component & a1 = r->a1();
      a1.set_in_component(i);
      a1.set_matched_atom(ri.which_is_mapped_atom(m1));

#ifdef DEBUG_IDENTIFY_ATOM_SUBSTITUTIONS
      cerr << "First part, reagent " << i << " matched atom " << ri.which_is_mapped_atom(m1) << endl;
#endif

      int rm2 = _reagent_locator[m2];

      const ISIS_RXN_FILE_Molecule & mrm2 = _reagent[rm2];
      atom_number_t l = mrm2.which_is_mapped_atom(m2);

#ifdef DEBUG_IDENTIFY_ATOM_SUBSTITUTIONS
      cerr << "Second part, reagent " << rm2 << " matched atom " << l << endl;
#endif

      Matched_Atom_in_Component & a2 = r->a2();
      a2.set_in_component(rm2);
      a2.set_matched_atom(l);

      atom_substitutions.add(r);

      Bond * b = new Bond(amap[j], m1, SINGLE_BOND);
      connections_involved_in_substitutions.add(b);
      b = new Bond(amap[j], m2, SINGLE_BOND);
      connections_involved_in_substitutions.add(b);
    }
  }

  return 1;
}

/*
  Lots of problems with the reaction


$RXN

      ISIS     102120021143

  2  2
$MOL

  -ISIS-  10210211432D

  4  3  0     0  0              1 V2000
    3.4030    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  1  0  0
    2.8660    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    2.0000    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  5  0  0
    2.8660   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END
$MOL

  -ISIS-  10210211432D

  2  1  0     0  0              1 V2000
    1.2500    0.0000    0.0000 Li  0  0  0  0  0  0  0  0  0  2  0  0
    2.2500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0
  1  2  1  0  0  0  0
M  END
$MOL

  -ISIS-  10210211432D

  4  3  0     0  0              1 V2000
    5.9030    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  1  0  0
    5.3660    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    4.5000    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  5  0  0
    5.3660   -0.7500    0.0000 Li  0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END
$MOL

  -ISIS-  10210211432D

  2  1  0     0  0              1 V2000
    7.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0
    8.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0
  1  2  1  0  0  0  0
M  END

  If you just process it, you get two atom substitutions, which turns out to
  be ghastly. It looks as if we need to identify the case

     1    3            1----3
     |    |
     |    |        ->
     |    |
     2    4            2----4

   Note that this ghastly procedure discovers every square bond pattern four times
*/

int
RXN_File::_identify_square_bonding_changes (int highest_atom_map)
{
  int ham1 = highest_atom_map + 1;

  DELETE_IF_NOT_NULL_ARRAY(_involved_in_square_bond);
  _involved_in_square_bond = new_int(ham1 * ham1);

  for (int i = 0; i <= highest_atom_map; i++)
  {
    int ri = _reagent_locator[i];
    if (ri < 0)
      continue;

    resizable_array<int> mapped_atoms_connected_to_i;

    const ISIS_RXN_FILE_Molecule & r = _reagent[ri];

    atom_number_t k = r.which_is_mapped_atom(i);

    if (k < 0)    // may have been removed
      continue;

    resizable_array<int> initially_connected_to_i;

    _reagent[ri].identify_singly_bonded_mapped_neighbours(k, initially_connected_to_i);

    if (0 == initially_connected_to_i.number_elements())
      continue;

    int p = _product_locator[i];   
    if (p < 0)    // mapped on LHS, but lost during reaction
      continue;

    const ISIS_RXN_FILE_Molecule & pm = _product[p];

    k = pm.which_is_mapped_atom(i);

    assert (INVALID_ATOM_NUMBER != k);

    resizable_array<int> ultimately_connected_to;

    pm.identify_singly_bonded_mapped_neighbours(k, ultimately_connected_to);

    if (0 == ultimately_connected_to.number_elements())
      continue;

//  At this stage, we have triangles.

    for (int j = 0; j < initially_connected_to_i.number_elements(); j++)
    {
      int k = initially_connected_to_i[j];

      for (int l = 0; l < ultimately_connected_to.number_elements(); l++)
      {
        int u = ultimately_connected_to[l];

//#define DEBUG_IDENTIFY_SQUARE_BONDS
#ifdef DEBUG_IDENTIFY_SQUARE_BONDS
        cerr << "Is there a square involving " << k << ' ' << i << " and " << u << endl;
#endif

        int m = _identify_square(k, i, u);
        if (m < 0)
          continue;

        _involved_in_square_bond[i * ham1 + k] = _involved_in_square_bond[k * ham1 + i] = 1;
        _involved_in_square_bond[i * ham1 + u] = _involved_in_square_bond[u * ham1 + i] = 1;
        _involved_in_square_bond[m * ham1 + u] = _involved_in_square_bond[u * ham1 + m] = 1;
        _involved_in_square_bond[m * ham1 + k] = _involved_in_square_bond[k * ham1 + m] = 1;

//      Note we do the "diagonal" bonds too

        _involved_in_square_bond[i * ham1 + m] = _involved_in_square_bond[m * ham1 + i] = 1;
        _involved_in_square_bond[u * ham1 + k] = _involved_in_square_bond[k * ham1 + u] = 1;

#ifdef DEBUG_IDENTIFY_SQUARE_BONDS
        cerr << "Mapped " << i << ' ' << k << ' ' << u << ' ' << m << " involved in square bond\n";
#endif
      }
    }
  }

  return 1;
}

/*
  M1-M2-M3

  M1-M2 in the reagents, M2-M3 in the products.

  can we find M4 that is bonded to both M1 in the reagents and M3 in the products
*/

int
RXN_File::_identify_square (int m1, int m2, int m3) const
{
  int r = _reagent_locator[m3];

  if (r < 0)
    return 0;

  const ISIS_RXN_FILE_Molecule & rm = _reagent[r];

  atom_number_t i = rm.which_is_mapped_atom(m3);

  assert (i >= 0);

  resizable_array<int> initially_connected_m3;

  rm.identify_singly_bonded_mapped_neighbours(i, initially_connected_m3);

  int p = _product_locator[m1];

  if (p < 0)
  {
//  cerr << "RXN_File::_identify_square:no product contains mapped atom " << m1 << endl;
    return 0;
  }

  const ISIS_RXN_FILE_Molecule & pm = _product[p];

  i = pm.which_is_mapped_atom(m1);

  assert (i >= 0);

  resizable_array<int> ultimately_connected_m1;

  pm.identify_singly_bonded_mapped_neighbours(i, ultimately_connected_m1);

  for (int i = 0; i < initially_connected_m3.number_elements(); i++)
  {
    int j = initially_connected_m3[i];

#ifdef DEBUG_IDENTIFY_SQUARE_BONDS
    cerr << j << " is initialially connected to " << m3 << endl;
#endif

    if (j == m2)
      continue;

    if (ultimately_connected_m1.contains(j))
      return j;
  }

  return -1;
}

//#define DEBUG_IDENTIFY_BONDS_TO_BE_BROKEN

int
RXN_File::_identify_bonds_to_be_broken (int highest_atom_map,
                                        resizable_array_p<Bond> & bonds_to_be_broken,
                                        const resizable_array_p<Bond> & connections_involved_in_substitutions,
                                        const int * include_these_atom_map_numbers) const
{
  for (int i = 1; i <= highest_atom_map; i++)
  {
    if (NULL == include_these_atom_map_numbers)
      ;
    else if (0 == include_these_atom_map_numbers[i])
      continue;

    int ri = _reagent_locator[i];
    if (ri < 0)
      continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_BROKEN
    cerr << "Mapped atom " << i << " in reagent " << _reagent_locator[i] << " and product " << _product_locator[i] << endl;
#endif

    for (int j = i + 1; j <= highest_atom_map; j++)
    {
      if (NULL == include_these_atom_map_numbers)
        ;
      else if (0 == include_these_atom_map_numbers[j])
        continue;

      bond_type_t initial_bond_type = _initial_bond_type[_btype_dim * i + j];

      if (INVALID_BOND_TYPE == initial_bond_type)   // not bonded
        continue;

      bond_type_t final_bond_type   = _final_bond_type[_btype_dim * i + j];

      if (initial_bond_type == final_bond_type)
        continue;

      if (INVALID_BOND_TYPE != final_bond_type)   // still bonded, just the type of bond has changed. No need to remove the bond
        continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_BROKEN
      cerr << " mapped atom " << j << " in reagent " << _reagent_locator[j] << " and product " << _product_locator[j] << endl;
      cerr << " initially bonded. Products " << _product_locator[i] << " and " << _product_locator[j] << endl;
#endif

      if (bond_already_present(connections_involved_in_substitutions, i, j))
        continue;

      Bond * b = new Bond(i, j, SINGLE_BOND);
      bonds_to_be_broken.add(b);
    }
  }

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_BROKEN
  cerr << "Found " << bonds_to_be_broken.number_elements() << " bonds to be broken\n";
#endif

  return 1;
}

//#define DEBUG_IDENTIFY_BONDS_TO_BE_MADE

/*
  Bonds to be made include the case where there are two atoms that remain bonded, but
  with a different bond type
*/

int
RXN_File::_identify_bonds_to_be_made (int highest_atom_map,
                                      resizable_array_p<Bond> & bonds_to_be_made,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions,
                                      const int * include_these_atom_map_numbers)
{
#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_MADE
  cerr << "RXN_File::_identify_bonds_to_be_made:highest_atom_map " << highest_atom_map << endl;
#endif

  resizable_array<int> orphans;

  for (int i = 1; i <= highest_atom_map; i++)
  {
    if (NULL == include_these_atom_map_numbers)
      ;
    else if (! include_these_atom_map_numbers[i])
      continue;

    const int ri = _reagent_locator[i];
    const int pi = _product_locator[i];

    if (ri >= 0)    // great, appears in a reagent
      ;
    else if (pi >= 0)    // only appears in products, is orphan
    {
      orphans.add(i);
      continue;
    }
    else          // not sure what is happening here, don't think this will ever happen
      continue;

    if (pi < 0)            // mapped atom disappears - but then it would not be mapped, should not happen, just being safe
      continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_MADE
    atom_number_t t = _reagent[ri].which_is_mapped_atom(i);
    cerr << "Mapped atom " << i << " (" << _reagent[ri].atomic_symbol(t) << ") in reagent " << _reagent_locator[i] << " and product " << _product_locator[i] << endl;
#endif

    _identify_bonding_changes_involving_matched_atom(i, _reagent[ri], _product[pi], highest_atom_map,
                                                     bonds_to_be_made, connections_involved_in_substitutions, include_these_atom_map_numbers);
  }

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_MADE
  cerr << "Found " << bonds_to_be_made.number_elements() << " bonds to be made, orphans? " << orphans.number_elements() << endl;
#endif

  if (_orphan_atoms.natoms())
    _identify_bonds_to_be_made_involving_orphan(bonds_to_be_made, connections_involved_in_substitutions);
    
  if (0 == orphans.number_elements())
    return 1;

  return _identify_bonds_to_be_made_involving_orphan(orphans, bonds_to_be_made, connections_involved_in_substitutions);
}

int
RXN_File::_identify_bonding_changes_involving_matched_atom(const int mstart,
                                      ISIS_RXN_FILE_Molecule & reagent,
                                      ISIS_RXN_FILE_Molecule & product,
                                      const int highest_atom_map,
                                      resizable_array_p<Bond> & bonds_to_be_made,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions,
                                      const int * include_these_atom_map_numbers) const
{
  for (int j = mstart + 1; j <= highest_atom_map; j++)   // check all the rest of the matched atoms
  {
    int rj = _reagent_locator[j];
    if (rj < 0)
      continue;

    const bond_type_t final_bond_type = _final_bond_type[_btype_dim * mstart + j];
      
    if (INVALID_BOND_TYPE == final_bond_type)   // not bonded when done
      continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_MADE
    cerr << " mapped atom " << j << " in reagent " << _reagent_locator[j] << " and product " << _product_locator[j] << ". Initial " << _initial_bond_type[_btype_dim * mstart + j] << " final " << final_bond_type << endl;
#endif

    const bond_type_t initial_bond_type = _initial_bond_type[_btype_dim * mstart + j];

    if (initial_bond_type == final_bond_type)     // no change
      continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_MADE
    cerr << "Making bond between mapped atoms " << mstart << " and " << j << " type " << final_bond_type << endl;
#endif

//  Maybe we can also discern some info about the required bonding in the reagents

    int final_nbonds = product.nbonds_for_mapped_atom(mstart);

    reagent.discern_initial_conditions(mstart, final_nbonds);

    final_nbonds = product.nbonds_for_mapped_atom(j);

    _reagent[rj].discern_initial_conditions(j, final_nbonds);

    if (bond_already_present(connections_involved_in_substitutions, mstart, j))
      continue;

    Bond * b = new Bond(mstart, j, final_bond_type);
    bonds_to_be_made.add(b);
  }

  return 1;
}

//#define DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN

int
RXN_File::_identify_bonds_to_be_made_involving_orphan(const resizable_array<int> & orphan,
                                      resizable_array_p<Bond> & bonds_to_be_made,
                                      const resizable_array_p<Bond> & connections_involved_in_substitutions)
{
#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
  for (int i = 0; i < _np; ++i)
  {
    cerr << "product " << i << ' ' << _product[i].smiles() << endl;
  }

  for (int i = 1; i <= 13; ++i)
  {
    cerr << "_reagent_locator " << _reagent_locator[i] << " _product_locator " << _product_locator[i] << endl;
  }
#endif

  for (int i = 0; i < orphan.number_elements(); ++i)
  {
    const int oi = orphan[i];
    _identify_bonds_to_be_made_involving_orphan_atom(oi, bonds_to_be_made, connections_involved_in_substitutions);
    continue;

//  once I am convinced everything is working, get rid of all the rest of this code

    int pi = _product_locator[oi];

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
    cerr << "Orphan atom " << oi << " in product " << pi << endl;
#endif

    const atom_number_t a = _product[pi].which_is_mapped_atom(oi);
    assert (a >= 0);

    const Atom * o = _product[pi].atomi(a);

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
    cerr << "Bonding change involving orphan, mapped atom " << i << " ncon " << o->ncon() << endl;
#endif

    for (int j = 0; j < o->ncon(); ++j)
    {
      const Bond * b = o->item(j);

      const atom_number_t k = b->other(a);

      const int mk = _product[pi].atom_map(k);

      const int rk = _reagent_locator[mk];

      if (rk < 0)     // another orphan
        continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
      cerr << "  bonded to " << k << " which is map " << mk << " in reagent " << rk << endl;
#endif

      Bond * nb = new Bond(oi, mk, (int)b->btype());
      bonds_to_be_made.add(nb);

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
      cerr << "Added bond btw mapped atoms " << oi << " and " << mk << endl;
#endif
    }
  }

  return 1;
}

int
RXN_File::_identify_bonds_to_be_made_involving_orphan(resizable_array_p<Bond> & bonds_to_be_made,
                                                      const resizable_array_p<Bond> & connections_involved_in_substitutions) const
{
//cerr << "ORPHA has " << _orphan_atoms.natoms() << " atoms\n";

  for (int i = 0; i < _orphan_atoms.natoms(); ++i)
  {
    const int oi = _orphan_atoms.atom_map(i);
//  cerr << "atom map " << oi << " is in orphan\n";

    _identify_bonds_to_be_made_involving_orphan_atom(oi, bonds_to_be_made, connections_involved_in_substitutions);
  }

//cerr << "LINE " << __LINE__ << endl;

  return 1;
}

/*
  Common code for dealing with an orphan atom, regardless of how it is perceived
*/

int
RXN_File::_identify_bonds_to_be_made_involving_orphan_atom(const int oi,     // atom map number of the orphan atom
                                                      resizable_array_p<Bond> & bonds_to_be_made,
                                                      const resizable_array_p<Bond> & connections_involved_in_substitutions) const
{
  int pi = _product_locator[oi];
  if (pi < 0)     // likely no atom map number
    return 1;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
  cerr << "Orphan atom " << oi << " in product " << pi << endl;
#endif

  const atom_number_t a = _product[pi].which_is_mapped_atom(oi);
  assert (a >= 0);

  const Atom * o = _product[pi].atomi(a);

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
  cerr << "Bonding change involving orphan, atom map number " << oi << " ncon " << o->ncon() << endl;
#endif

  for (int j = 0; j < o->ncon(); ++j)
  {
    const Bond * b = o->item(j);

    const atom_number_t k = b->other(a);

    const int mk = _product[pi].atom_map(k);

    const int rk = _reagent_locator[mk];

    if (rk < 0)     // another orphan
      continue;

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
    cerr << "  bonded to " << k << " which is map " << mk << " in reagent " << rk << endl;
#endif

    Bond * nb = new Bond(oi, mk, (int)b->btype());
    bonds_to_be_made.add(nb);

#ifdef DEBUG_IDENTIFY_BONDS_TO_BE_ADE_ORPHAN
    cerr << "Added bond btw mapped atoms " << oi << " and " << mk << endl;
#endif
  }

  return 1;
}

int
RXN_File::_identify_isotopes_to_be_placed (int highest_atom_map, 
                                           resizable_array_p<Reaction_Place_Isotope> & isotopes) const
{
  for (int i = 1; i <= highest_atom_map; i++)
  {
    int ra = _reagent_locator[i];
    if (ra <= 0)
      continue;

    int rb = _reagent[ra].which_is_mapped_atom(i);
    if (rb < 0)
      continue;

    int ri = _reagent[ra].isotope(rb);

    int pa = _product_locator[i];
    if (pa < 0)
      continue;
    int pb = _product[pa].which_is_mapped_atom(i);
    int pi = _product[pa].isotope(pb);

    if (ri == pi)
      continue;

    Reaction_Place_Isotope * tmp = new Reaction_Place_Isotope;
    tmp->set_atom(i);
    tmp->set_isotope(pi);

    isotopes.add(tmp);
  }

  return 1;
}

/*
  Atom A is involved in a bond addition - either changing an existing bond
  or adding a bond.
  We look at the difference in the number of bonds - not connectivity
*/

//#define DEBUG_DISCERN_INITIAL_CONDITIONS

int
ISIS_RXN_FILE_Molecule::discern_initial_conditions (int m,
                                                    int final_nbonds)
{
  atom_number_t a = which_is_mapped_atom(m);

  if (a < 0)
  {
    cerr << "ISIS_RXN_FILE_Molecule::discern_initial_conditions:no mapped atom " << m << endl;
    return 0;
  }

  MDL_Atom_Data * mad = _mdl_atom[a];

#ifdef DEBUG_DISCERN_INITIAL_CONDITIONS
  cerr << "discern_initial_conditions:mapped atom " << m << " is atom " << a << ", hcount " << mad->hcount() << endl;
#endif

  if (0 != mad->hcount())    // already set by user
    return 1;

  int initial_nbonds = nbonds(a);

  if (NULL != _connections_lost)
  {
//  cerr << "Atom " << a << " lost " << _connections_lost[a] << " connections\n";
    initial_nbonds += _connections_lost[a];   // assumes single bonds
  }

  int bonding_increase = final_nbonds - initial_nbonds;

#ifdef DEBUG_DISCERN_INITIAL_CONDITIONS
  cerr << "Mapped atom " << m << " atom " << a << " bonding_increase " << bonding_increase << endl;
#endif

//int min_hount;

  if (0 == bonding_increase)
    mad->set_min_hcount(_explicit_hydrogen_count[a]);
  else if (0 == _computed_implicit_hydrogen_count[a])   // it has no open valence
    mad->set_min_hcount(_explicit_hydrogen_count[a]);
  else if (_aromatic_bonds_attached(a))    // information is lost
    ;
  else if (1 == bonding_increase)
    mad->set_min_hcount(1);
  else if (2 == bonding_increase)
    mad->set_min_hcount(2);
  else if (3 == bonding_increase)
    mad->set_min_hcount(3);

#ifdef DEBUG_DISCERN_INITIAL_CONDITIONS
  cerr << "Mapped atom " << m << " (" << atomi(a)->atomic_symbol() << ") our atom " << a << " initial_nbonds " << initial_nbonds << " final_nbonds " << final_nbonds << " min_hcount " << mad->min_hcount() << endl;
#endif

  return 1;
}

int
ISIS_RXN_FILE_Molecule::_aromatic_bonds_attached (atom_number_t zatom) const
{
  const Atom * a = atomi(zatom);

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    int j = which_bond(b->a1(), b->a2());

    const MDL_Bond_Data * mdlb = mdl_bond(j);

    if (mdlb->btype() & AROMATIC_BOND)
      rc++;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::nbonds_for_mapped_atom (int m)
{
  atom_number_t zatom = which_is_mapped_atom(m);

  if (zatom < 0)
  {
    cerr << "ISIS_RXN_FILE_Molecule::nbonds_for_mapped_atom:no map for " << m << endl;
    for (int i = 0; i < Molecule::natoms(); i++)
    {
      cerr << " atom " << i << " map " << _atom_map[i] << endl;
    }
    abort();
  }

  return _nbonds[zatom];
}

/*
  Make sure there is an atom list present for every 'L' element
*/

int
ISIS_RXN_FILE_Molecule::_check_non_periodic_table_elements_and_atom_lists ()
{
  int matoms = Molecule::natoms();

  int rc = 1;

  for (int i = 0; i < matoms; i++)
  {
    const Element * e = Molecule::elementi(i);

    if (e->is_in_periodic_table())
      continue;

    const IWString & s = e->symbol();

    if ('A' == s)    // we know what to do with these
      continue;
    if ('Q' == s)    // we know what to do with these
      continue;

    if ('R' == s)    // just assume they wanted A
    {
      const Element * r = get_element_from_symbol_no_case_conversion("R", 1);
      Molecule::set_element(i, r);
      continue;
    }

    if ('L' != s)
    {
      cerr << "ISIS_RXN_FILE_Molecule::_check_non_periodic_table_elements_and_atom_lists:warning, non periodic table element encountered '" << s << "'\n";
      continue;
    }

    const MDL_Atom_Data * mad = _mdl_atom[i];

    const ISIS_Atom_List & ali = mad->atom_list();

    if (ali.active())    // great, L symbol is indeed an atom list
      continue;

    cerr << "ISIS_RXN_FILE_Molecule::_check_non_periodic_table_elements_and_atom_lists\n";
    cerr << "Atom " << i << " is 'L' type, but not an atom list\n";
    rc = 0;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::add_chiral_centres_to_be_inverted (Reaction_Site & r) const
{
  int matoms = Molecule::natoms();
  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * madi = _mdl_atom[i];

    if (0 == madi->inversion())
      continue;

    r.add_chiral_centre_to_be_inverted(i);

    Substructure_Atom * a = r.query_atom_with_initial_atom_number(i);
    a->set_chirality(1);
  }

  return 1;
}

/*static void
set_matched_atom (Stereo_Centre_Component & s,
                  int z)
{
  if (z >= 0)
    s.set_matched_atom(z);
  else
    s.set_implicit_hydrogen(1);

  return;
}*/

/*
  Strange function. We are building a chiral centre specification. 
  We need to specify the atom number of one of the appendages to the centre atom.

  CENTRE_ATOM_REAGENT is the reagent number of the centre atom.
  OUR_REAGENT is the reagent number of the atom we are interested in.

  If we are in the same reagent, we just use the atom number A.

  If we are from a different reagent, we need to ask that reagent which query
  atom number corresponds to the mapped atom by which we are known
*/

void
RXN_File::_specify_stereo_centre_component (Stereo_Centre_Component & s,
                                            int centre_atom_reagent,
                                            int our_reagent,
                                            int mapped) const
{
  if (mapped < 0)     // must be a lone pair or implicit hydrogen
  {
    s.set_implicit_hydrogen(1);
    return;
  }

  assert (our_reagent >= 0 && our_reagent < _nr);
 
// Need to get the atom number in that reagent

  const ISIS_RXN_FILE_Molecule & r = _reagent[our_reagent];

  int m = r.which_is_mapped_atom(mapped);

#ifdef DEBUG_SPECIFY_STEREO_CENTRE_COMPONENT
  cerr << "Mapped atom " << mapped << " component " << our_reagent << " m = " << m << endl;
#endif

  assert (m >= 0);

  if (0 == our_reagent)
    s.set_in_scaffold();
  else
    s.set_in_component(our_reagent - 1);

#ifdef DEBUG_SPECIFY_STEREO_CENTRE_COMPONENT
  cerr << "In component " << our_reagent << " mapped atom " << mapped << " is atom " << m << endl;
#endif

  s.set_matched_atom(m);

  return;
}

/*
  We have an atom that is part of a stereo centre specification. But it may be a lone pair
  or implicit hydrogen
*/

static int
find_reagent (const int * reagent_locator,
              int m)
{
  if (m < 0)
    return m;

  return reagent_locator[m];
}

static int
convert_to_mapped_atom (const ISIS_RXN_FILE_Molecule & m,
                        atom_number_t a,
                        int & ma)
{
  if (a < 0)     // implicit hydrogen or lone pair
  {
    ma = a;
    return 1;
  }

  ma = m.atom_map(a);

  if (ma < 1)
    return 0;

  return 1;
}

static int
collect_mapped_neighbours(const ISIS_RXN_FILE_Molecule & m,
                          const Chiral_Centre * c,
                          int & mtf,
                          int & mtb,
                          int & mld,
                          int & mrd)
{
//c->debug_print(cerr);
  atom_number_t a = c->top_front();
  if (! convert_to_mapped_atom(m, a, mtf))
    return 0;

  a = c->top_back();
  if (! convert_to_mapped_atom(m, a, mtb))
    return 0;

  a = c->left_down();
  if (! convert_to_mapped_atom(m, a, mld))
    return 0;

  a = c->right_down();
  if (! convert_to_mapped_atom(m, a, mrd))
    return 0;

  return 1;
}

//#define DEBUG_LOOK_FOR_STEREO_CENTRES_MADE

int
RXN_File::_look_for_stereo_centres_made (IWReaction & rxn)
{
  int return_code = 0;

  for (int i = 0; i < _np; i++)
  {
    const ISIS_RXN_FILE_Molecule & m = _product[i];

#ifdef DEBUG_LOOK_FOR_STEREO_CENTRES_MADE
    cerr << "Product " << i << " contains " << m.chiral_centres() << " chiral centres\n";
#endif
    for (int j = 0; j < m.chiral_centres(); j++)
    {
      const Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(j);

      int mc = m.atom_map(c->a());
//    cerr << "j = " << j << " mapped atom is " << mc << endl;
      if (mc <= 0)                     // not a  mapped atom, can't process these
        continue;

      int mtf, mtb, mld, mrd;
      if (! collect_mapped_neighbours(m, c, mtf, mtb, mld, mrd))
        continue;

#ifdef DEBUG_LOOK_FOR_STEREO_CENTRES_MADE
      cerr << "Mapped neighbours " << mc << " tf " << mtf << " tb " << mtb << " ld " << mld << " rd " << mrd << endl;
#endif

      int rc  = find_reagent(_reagent_locator, mc);

      if (rc < 0)    // chiral centre on atom not on LHS,
        continue;

      int rtf = find_reagent(_reagent_locator, mtf);     // in which reagent is mapped atom MTF
      if (mtf >= 0 && rtf < 0 )
        continue;   
      int rtb = find_reagent(_reagent_locator, mtb);
      if (mtb >= 0 &&  rtb < 0)
        continue;
      int rld = find_reagent(_reagent_locator, mld);
      if (mld >= 0 && rld < 0)
        continue;
      int rrd = find_reagent(_reagent_locator, mrd);
      if (mrd >= 0 && rrd < 0)
        continue;

#ifdef DEBUG_LOOK_FOR_STEREO_CENTRES_MADE
      cerr << "In components " << rc << " tf " << rtf << " tb " << rtb << " ld " << rld << " rd " << rrd << endl;
#endif

      Reaction_Stereo_Centre * r = new Reaction_Stereo_Centre;

      Stereo_Centre_Component & sc = r->centre();
      _specify_stereo_centre_component(sc, rc, rc, mc);
      Stereo_Centre_Component & stf = r->top_front();
      _specify_stereo_centre_component(stf, rc, rtf, mtf);
      Stereo_Centre_Component & stb = r->top_back();
      _specify_stereo_centre_component(stb, rc, rtb, mtb);
      Stereo_Centre_Component & sld = r->left_down();
      _specify_stereo_centre_component(sld, rc, rld, mld);
      Stereo_Centre_Component & srd = r->right_down();
      _specify_stereo_centre_component(srd, rc, rrd, mrd);

#ifdef DEBUG_LOOK_FOR_STEREO_CENTRES_MADE
      int pc = _reagent_locator[mc];

      cerr << "Stereo centre for mapped atom " << mc << " from reagent " << pc << endl;

      r->debug_print(cerr);
#endif

      rxn.add_reaction_stereo_centre(r);

      return_code++;
    }
  }

#ifdef DEBUG_LOOK_FOR_STEREO_CENTRES_MADE
  cerr << "Placed " << return_code << " stereo centres\n";
#endif

  return return_code;
}

int
write_isis_reaction_file_header (const IWString & comment,
                                 const int nr,
                                 const int np,
                                 std::ostream & output)
{
  output << "$RXN\n";
  output << comment << '\n';
  output << "      ISIS     102120021143\n";
  output << '\n';
  output << std::setw(3) << nr << std::setw(3) << np << '\n';

  return output.good();
}

int
RXN_File::do_write (std::ostream & output) const
{
  write_isis_reaction_file_header(_comment, _nr, _np, output);

  for (int i = 0; i < _nr; i++)
  {
    output << "$MOL\n";
    _reagent[i].do_write(output);
  }

  for (int i = 0; i < _np; i++)
  {
    output << "$MOL\n";
    _product[i].do_write(output);
  }

  return output.good();
}

int
RXN_File::do_write (const char * fname) const
{
  std::ofstream output(fname, std::ios::out);

  if (! output.good())
  {
    cerr << "RXN_File::do_write:cannot open '" << fname << "'\n";
    return 0;
  }

  return do_write(output);
}

static int
write_m_values (int matoms,
                int n,
                int max_per_record,
                const char * prefix,
                const int * values,
                std::ostream & output)
{
  assert (3 == ::strlen(prefix));

  int j = 0;    // index into VALUES array

  while (n > 0)
  {
    int items_this_record;

    if (n > max_per_record)
      items_this_record = max_per_record;
    else
     items_this_record = n;

    output << "M  " << prefix << std::setw(3) << items_this_record;

    while (j < matoms)
    {
      int s = values[j];
      j++;
      if (0 == s)
        continue;

      output << std::setw(4) << j << std::setw(4) << s;
    }

    output << '\n';

    n -= items_this_record;
  }

  return output.good();
}

int
ISIS_RXN_FILE_Molecule::_write_m_sub_records (int n,
                                        std::ostream & output) const
{
  int matoms = Molecule::natoms();

  int * substitution = new int[matoms]; std::unique_ptr<int[]> free_substitution(substitution);

  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * mad = MDL_Molecule::mdl_atom_data(i);

    substitution[i] = mad->substitution();
  }

  return write_m_values(matoms, n, 8, "SUB", substitution, output);
}


int
ISIS_RXN_FILE_Molecule::_write_m_uns_records (int n,
                                              std::ostream & output) const
{
  int matoms = Molecule::natoms();

  int * unsaturated = new int[matoms]; std::unique_ptr<int[]> free_unsaturated(unsaturated);

  for (int i = 0; i < matoms; i++)
  {
    unsaturated[i] = _mdl_atom[i]->unsaturated();
  }

  return write_m_values(matoms, n, 8, "UNS", unsaturated, output);
}

int
ISIS_RXN_FILE_Molecule::_write_m_rbc_records (int n,
                                              std::ostream & output) const
{
  int matoms = Molecule::natoms();

  int * rbc = new int[matoms]; std::unique_ptr<int[]> free_rbc(rbc);

  for (int i = 0; i < matoms; i++)
  {
    rbc[i] = _mdl_atom[i]->ring_bond();
  }

  return write_m_values(matoms, n, 8, "RBC", rbc, output);
}

static int
btype_back_to_mdl_form (bond_type_t b)
{
  if (SINGLE_BOND == b)
    return 1;

  if (DOUBLE_BOND == b)
    return 2;

  if (TRIPLE_BOND == b)
    return 3;

  if (AROMATIC_BOND == b)
    return 4;

  if ((SINGLE_BOND | DOUBLE_BOND) == b)
    return 5;

  if ((AROMATIC_BOND | SINGLE_BOND) == b)
    return 6;

  if ((AROMATIC_BOND | DOUBLE_BOND) == b)
    return 7;

  if ((SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND) == b)
    return 8;

  cerr << "btype_back_to_mdl_form: warning, unimplemented feature " << b << endl;

  return 5;
}

int
ISIS_RXN_FILE_Molecule::_number_atoms_with_unsaturation_directives() const
{
  int matoms = Molecule::natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * madi = _mdl_atom[i];

    if (0 != madi->unsaturated())
      rc++;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::_number_atoms_with_substitution_directives() const
{
  int matoms = Molecule::natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * madi = _mdl_atom[i];

    if (0 != madi->substitution())
      rc++;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::_number_atoms_with_ring_bond_directives() const
{
  int matoms = Molecule::natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * madi = _mdl_atom[i];

    if (0 != madi->ring_bond())
      rc++;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::do_write (std::ostream & output) const
{
  MDL_File_Supporting_Material * mdlsm = global_default_MDL_File_Supporting_Material();
  set_write_isis_standard(1);

//output << "$MOL\n";
  output << molecule_name() << '\n';
  output << "  -ISIS-  10210211432D\n";
  output << '\n';

  int nfc = number_formally_charged_atoms();
  int iat = number_isotopic_atoms();

  _mdl_write_atoms_and_bonds_record(output, nfc, iat, *mdlsm);

  const int matoms = Molecule::natoms();

  int write_madi_atom_map = 1;
  if (NULL != _atom_map)
    write_madi_atom_map = ! std::find_if(_atom_map, _atom_map + matoms, [] (int s) { return s > 0;});

  IWString buffer;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = atomi(i);

    a->write_coordinates(output);

    _write_mdl_atom_record_element_charge_and_chirality(i, buffer, *mdlsm);

    output << buffer;

    const MDL_Atom_Data * madi = _mdl_atom[i];

    output << std::setw(3) << madi->hcount() << "  0  0";

    output << std::setw(3) << madi->h0designator() << "  0  0";

    if (write_madi_atom_map)
      output << std::setw(3) << madi->atom_map();
    else
      output << std::setw(3) << _atom_map[i];

    output << std::setw(3) << madi->inversion();

    output << "  0";

    output << '\n';
  }

  int nb = Molecule::nedges();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = Molecule::bondi(i);

    output << std::setw(3) << (b->a1() + 1) << std::setw(3) << (b->a2() + 1);

    const MDL_Bond_Data * mdlb = mdl_bond(i);

    output << std::setw(3) << btype_back_to_mdl_form(mdlb->btype());

    output << "  0";    // bond stereochemistry
    output << "  0";    // not used by MDL

    const MDL_Bond_Data * mdlbi = _mdl_bond[i];

    output << std::setw(3) << mdlbi->bond_topology();
    output << std::setw(3) << mdlbi->reacting_centre_status();
    output << '\n';
  }

  int n = number_isotopic_atoms();
  if (n)
    _write_m_iso_records(output, n);

  n = number_formally_charged_atoms();
  if (n)
    _write_m_chg_records(output, n);

  n = _number_atoms_with_substitution_directives();
  if (n)
    _write_m_sub_records(n, output);

  n = _number_atoms_with_unsaturation_directives();
  if (n)
    _write_m_uns_records(n, output);

  n = _number_atoms_with_ring_bond_directives();
  if (n)
    _write_m_rbc_records(n, output);

  for (int i = 0; i < matoms; i++)
  {
    _mdl_atom[i]->write_M_ALS(i, output);
  }

  for (int i = 0; i < _link_atom.number_elements(); i++)
  {
    _link_atom[i]->write_M_LIN(0, output);   // we have removed the L atoms, and don't bother putting them back
  }

  output << "M  END\n";

  return output.good();
}

/*
  We have found an aromatic ring that may need its Kekule form toggled. We need
  to identify all aromatic rings that are fused
*/

static int
mark_all_fused_rings_done (Molecule & m, 
                           int * ring_already_done,
                           const Ring & r)
{
  ring_already_done[r.ring_number()] = 1;

  int n = r.fused_ring_neighbours();

  for (int i = 0; i < n; i++)
  {
    const Ring * ri = r.fused_neighbour(i);

    int rn = ri->ring_number();

    if (! ri->is_aromatic())
      continue;

    if (ring_already_done[rn])
      continue;

    mark_all_fused_rings_done(m, ring_already_done, *ri);
  }

  return 1;
}

static int
contains_adjacent_mapped_atoms (const ISIS_RXN_FILE_Molecule & m,
                                const Ring & r,
                                atom_number_t & a1,
                                atom_number_t & a2)
{
  for (Ring_Bond_Iterator i(r); i != r.end(); i++)
  {
    const atom_number_t ra1 = i.a1();
    const atom_number_t ra2 = i.a2();

    if (m.atom_map(ra1) < 1 || m.atom_map(ra2) < 1)
      continue;

    a1 = ra1;
    a2 = ra2;

    return 1;
  }

  return 0;
}

int
RXN_File::_reagent_ring_contains_adjacent_mapped_atoms (const ISIS_RXN_FILE_Molecule & m,
                                const Ring & r,
                                atom_number_t & a1,
                                atom_number_t & a2,
                                const int * include_these_atom_map_numbers) const
{
  for (Ring_Bond_Iterator i(r); i != r.end(); i++)
  {
    const atom_number_t ra1 = i.a1();
    const atom_number_t ra2 = i.a2();

    const int amap1 = m.atom_map(ra1);
    if (amap1 < 1)
      continue;

    const int amap2 = m.atom_map(ra2);
    if (amap2 < 1)
      continue;

    if (NULL == include_these_atom_map_numbers)
      ;
    else if (! include_these_atom_map_numbers[amap1] || ! include_these_atom_map_numbers[amap2])
      continue;

    a1 = ra1;
    a2 = ra2;

    return 1;
  }

  return 0;
}

int
RXN_File::_max_rings_in_any_reagent() const
{
  if (0 == _nr)
    return 0;

  int rc = _reagent[0].nrings();

  for (int i = 1; i < _nr; ++i)
  {
    if (_reagent[i].nrings() > rc)
      rc = _reagent[i].nrings();
  }

  return rc;
}

//#define DEBUG_IDENTIFY_KEKULE_FORMS_TO_BE_TOGGLED

/*
  This is not correct.
  We cannot assume that all rings in a fused system can be covered by just one toggling, we need
  to look at each fused aromatic system separately. Fix sometime!
*/

int
RXN_File::_identify_kekule_forms_to_be_toggled (const int * include_these_atom_map_numbers)
{
  const int maxrings = _max_rings_in_any_reagent();

#ifdef DEBUG_IDENTIFY_KEKULE_FORMS_TO_BE_TOGGLED
  cerr << "_identify_kekule_forms_to_be_toggled:maxrings " << maxrings << endl;
#endif

  if (0 == maxrings)
    return 1;

  int * ring_already_done = new int[maxrings + 1]; std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int rc = 0;

  for (int i = 0; i < _nr; i++)
  {
    ISIS_RXN_FILE_Molecule & ri = _reagent[i];

    const int nr = ri.nrings();

    if (0 == nr)    // no rings in the reagent, nothing to do
      continue;

    ri.compute_aromaticity_if_needed();    // probably not needed...

    std::fill_n(ring_already_done, nr + 1, 0);

    for (int j = 0; j < nr; j++)
    {
      if (ring_already_done[j])
        continue;

      const Ring * r = ri.ringi(j);

      if (6 != r->number_elements())    // only 6 membered rings can have differing Kekule forms
        continue;

      if (! r->is_aromatic())
        continue;

      atom_number_t a1, a2;
      if (! _reagent_ring_contains_adjacent_mapped_atoms(ri, *r, a1, a2, include_these_atom_map_numbers))
        continue;

#ifdef DEBUG_IDENTIFY_KEKULE_FORMS_TO_BE_TOGGLED
      cerr << "What about ring " << (*r) << endl;
      cerr << " adjacent mapped atoms " << ri.atom_map(a1) << " and " << ri.atom_map(a2) << endl;
#endif

      if (0 == _changes_in_ring_membership(ri, *r))
        continue;

      bond_type_t bt = ri.btype_between_atoms(a1, a2);

#ifdef DEBUG_IDENTIFY_KEKULE_FORMS_TO_BE_TOGGLED
      cerr << "Adding toggle Kekule form beetween atom " << a1 << " (map " << ri.atom_map(a1) << ") and " << a2 << "(map " << ri.atom_map(a2) << ") type " << bt << endl;
#endif

      ri.add_toggle_kekule_form(a1, a2, bt);

      ring_already_done[j] = 1;

      if (r->is_fused())
        mark_all_fused_rings_done(ri, ring_already_done, *r);

      rc++;
    }
  }

  return rc;
}

int
RXN_File::_changes_in_ring_membership (ISIS_RXN_FILE_Molecule & mfrom,
                                       const Ring & r) const
{
  int n = r.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    int m = mfrom.atom_map(j);

    if (m <= 0)
      continue;

    int p = _product_locator[m];

    if (p < 0)     // mapped atom does not appear on RHS
      continue;

    ISIS_RXN_FILE_Molecule & pm = _product[p];

    atom_number_t k = pm.which_is_mapped_atom(m);

    assert (k >= 0 && k < pm.natoms());

    if (pm.ring_bond_count(k) == mfrom.ring_bond_count(j))    // no change in ring membership
        continue;

    rc++;
  }

  return rc;
}

int
ISIS_RXN_FILE_Molecule::add_toggle_kekule_form (atom_number_t a1,
                                                atom_number_t a2,
                                                bond_type_t bt)
{
  bt = BOND_TYPE_ONLY(bt);

  Bond * b = new Bond(a1, a2, bt);

  _toggle_kekule_form.add(b);

#ifdef DEBUG_IDENTIFY_KEKULE_FORMS_TO_BE_TOGGLED
  cerr << "Toggle Kekule form between atoms " << a1 << " and " << a2 << " type " << bt << endl;
#endif

  return 1;
}

int
ISIS_RXN_FILE_Molecule::add_toggle_kekule_forms (Reaction_Site & rxn,
                                                 const Reaction_Subset & subset,
                                                 const int ndx,
                                                 const RXN_File_Create_Reaction_Options & rxnfcro)
                                                  const
{
  int nt = _toggle_kekule_form.number_elements();

//cerr << "ISIS_RXN_FILE_Molecule::add_toggle_kekule_forms:adding " << nt << " Kekule togglings\n";

  for (int i = 0; i < nt; i++)
  {
    const Bond * b = _toggle_kekule_form[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
  
	  if (ndx == 0 || ! rxnfcro.only_create_query_from_first_reagent())
    {
	    a1 = subset.atom_number_in_subset(a1, ndx);
	    a2 = subset.atom_number_in_subset(a2, ndx);
	  }
	  
	  const atom_number_t s1 = a1;
	  const atom_number_t s2 = a2;


    if (s1 < 0 || s2 < 0)
    {
      cerr << "ISIS_RXN_FILE_Molecule::add_toggle_kekule_forms:atom " << a1 << " (" << s1 << ") and/or atom " << a2 << " (" << s2 << ") not found in subset, should not happen\n";
      abort();
      return 0;
    }

//  cerr << "Toggle kekule form, " << b->a1() << " and " << b->a2() << " type " << b->btype() << endl;

    rxn.add_toggle_kekule_form(s1, s2, b->btype());
  }

  return 1;
}

int
ISIS_RXN_FILE_Molecule::_compute_implicit_hydrogens ()
{
  int matoms = natoms();

//cerr << "Computing implicit Hydrogens for " << matoms << " atoms\n";

  if (NULL == _computed_implicit_hydrogen_count)
    _computed_implicit_hydrogen_count = new int[matoms];

  if (NULL == _explicit_hydrogen_count) 
    _explicit_hydrogen_count = new int[matoms];

  for (int i = 0; i < matoms; i++)
  {
    _explicit_hydrogen_count[i] = explicit_hydrogens(i);

    _computed_implicit_hydrogen_count[i] = implicit_hydrogens(i);
  }

  int nb = nedges();

  for (int i = 0; i < nb; i++)
  {
    const MDL_Bond_Data * mdlb = mdl_bond(i);

    if ((DOUBLE_BOND | AROMATIC_BOND) == mdlb->btype())
      ;
    else if (AROMATIC_BOND == mdlb->btype())
      ;
    else
      continue;

    const Bond * b = bondi(i);

    _adjust_computed_implicit_hydrogens_for_unstautration(b->a1());
    _adjust_computed_implicit_hydrogens_for_unstautration(b->a2());
  }

  return 1;
}

/*
  What if we have a 2 connected atom where one of the bonds is double or aromatic.
  We don't know anything about the required hcount, since the bonds drawn may be
  part of the aromatic ring or not part of the ring
*/

int
ISIS_RXN_FILE_Molecule::_adjust_computed_implicit_hydrogens_for_unstautration (atom_number_t zatom)
{
  if (2 == ncon(zatom))
    _computed_implicit_hydrogen_count[zatom] = 0;
  else if (_computed_implicit_hydrogen_count[zatom] > 0)
  {
    _computed_implicit_hydrogen_count[zatom]--;
//  cerr << "_computed_implicit_hydrogen_count for atom " << zatom << " set to " << _computed_implicit_hydrogen_count[zatom] << endl;
  }
 
  return 1;
}

static int
display_standard_rxn_file_options (std::ostream & os,
                                   char flag)
{
  os << "  -" << flag << " rmfrag    remove small fragments from products\n";
  os << "  -" << flag << " rmunmp    remove unmapped reagent atoms that don't appear in the products\n";
  os << "  -" << flag << " aromb     aromatic bonds lose their Kekule forms and become pure aromatic\n";
  os << "  -" << flag << " hswap     rearrange the atoms to put heteroatoms first - for efficiency\n";
  os << "  -" << flag << " noautoam  suppress automatic atom mapping\n";
  os << "  -" << flag << " SR<nn>    reagent <nn> is a single reagent, not a substructure query\n";
  os << "  -" << flag << " rmxh      remove explicit Hydrogen atoms in the query\n";
  os << "  -" << flag << " ucumrapve unconnect unmapped reagent atoms that would violate a\n";
  os << "                    valence in the product\n";
  os << "  -" << flag << " alias2iso convert atom aliases to isotopes\n";
  os << "  -" << flag << " cmte      a carbon atom on the RHS of an atom list means transform the element\n";
  os << "  -" << flag << "aromAC     detect aromaticity by temporarily changing A atoms to C\n";

  return os.good();
}

int
parse_isis_rxn_file_options (Command_Line & cl,
                             char flag,
                             RXN_File & ISIS_rxn)
{
  int i = 0;
  const_IWSubstring d;
  const_IWSubstring fname;

  while (cl.value(flag, d, i++))
  {
    if ("rmfrag" == d)
    {
      ISIS_rxn.set_remove_product_fragments(1);
    }
    else if ("rmunmp" == d)
    {
      ISIS_rxn.set_remove_unmapped_atoms_that_disappear(1);
    }
    else if ("aromb" == d)
    {
      ISIS_rxn.set_aromatic_bonds_lose_kekule_identity(1);
    }
    else if ("hswap" == d)
    {
      ISIS_rxn.set_swap_atoms_to_put_rare_atoms_first(1);
    }
    else if ("noautoam" == d)
    {
      ISIS_rxn.set_do_automatic_atom_mapping(0);
    }
    else if ("orphan" == d)
    {
      ISIS_rxn.set_auto_fix_orphans(1);
    }
    else if (d.starts_with("echo="))
    {
      d.remove_leading_chars(5);
      ISIS_rxn.set_fname_for_echo(d);
    }
    else if (d.starts_with("SR"))
    {
      d.remove_leading_chars(2);
      int i;
      if (! d.numeric_value(i) || i < 1)
      {
        cerr << "Invalid reagent number '" << d << "'\n";
        return 0;
      }

      i--;

      ISIS_rxn.set_molecule_is_sole_reagent(i, 1);
    }
    else if ("rmxh" == d)
    {
      ISIS_rxn.set_remove_explicit_hydrogens(1);
    }
    else if ("ucumrapve" == d)
    {
      ISIS_rxn.set_unconnect_unmapped_atoms_that_exceed_product_valence(1);
    }
    else if ("aromAC" == d)
    {
      ISIS_rxn.set_convert_A_to_C_for_aromaticity(1);
    }
    else if ("cmte" == d)
    {
      ISIS_rxn.set_interpret_carbon_on_rhs_of_list_as_no_change(0);
    }
    else if ("alias2iso" == d)
    {
      ISIS_rxn.set_convert_atom_aliases_to_isotopes(1);
    }
    else if ("help" == d)
    {
      display_standard_rxn_file_options(cerr, flag);
      exit(0);
    }
    else if (0 == fname.length())
      fname = d;
    else
    {
      cerr << "Unrecognised -D qualifier '" << d << "'\n";
      return 0;
    }
  }

  if (0 == fname.length())
  {
    cerr << "Must specify the ISIS reaction file to use via the -" << flag << " option\n";
    return 0;
  }

  if (! ISIS_rxn.do_read(fname))
  {
    cerr << "Cannot read ISIS reaction '" << fname << "'\n";
    return 0;
  }

  return 1;
}

template <typename T>
void
swap_array_elements (T * v,
                     int i1,
                     int i2)
{
  T tmp = v[i1];
  v[i1] = v[i2];
  v[i2] = tmp;

  return;
}

int
ISIS_RXN_FILE_Molecule::_do_swap_atoms_to_put_rare_atoms_first ()
{
  int matoms = Molecule::natoms();

  int first_heteroatom = -1;
  int first_atom_list = -1;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = Molecule::atomic_number(i);

    if (1 == z)    // Yuck, an explicit Hydrogen
      continue;

    if (is_an_atom_list(i))
    {
      first_atom_list = i;
      continue;
    }

    if (6 == z)
      continue;

    first_heteroatom = i;
    break;
  }

  if (0 == first_heteroatom)    // great, already got a heteroatom first
    return 1;

  if (first_heteroatom < 0 && first_atom_list < 0)   // all carbons and no atom lists!
    return 0;

  if (first_heteroatom < 0)
    first_heteroatom = first_atom_list;
  else
    first_atom_list = -1;

  MDL_Molecule::swap_atoms(0, first_heteroatom);

  swap_array_elements(_atom_map, 0, first_heteroatom);

  if (NULL != _computed_implicit_hydrogen_count)
    swap_array_elements(_computed_implicit_hydrogen_count, 0, first_heteroatom);

  if (NULL != _explicit_hydrogen_count)
    swap_array_elements(_explicit_hydrogen_count, 0, first_heteroatom);

  if (NULL != _explicit_hydrogen_atom_removed)
    swap_array_elements(_explicit_hydrogen_atom_removed, 0, first_heteroatom);

  if (NULL != _connections_lost)
    swap_array_elements(_connections_lost, 0, first_heteroatom);

  if (NULL != _nbonds)
    swap_array_elements(_nbonds, 0, first_heteroatom);

  for (int i = 0; i < _link_atom.number_elements(); i++)
  {
    _link_atom[i]->swap_atoms(first_heteroatom, 0);
  }

  return 1;
}

int
ISIS_RXN_FILE_Molecule::transfer_atom_alias_to_isotope()
{
  int n = natoms();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    const MDL_Atom_Data * madi = mdl_atom(i);

    const IWString alias = madi->alias();

    if (0 == alias.length())
      continue;

    int iso;
    if (! alias.numeric_value(iso) || iso < 0)
    {
      cerr << "ISIS_RXN_FILE_Molecule::transfer_atom_alias_to_isotope:invalid numeric alias '" << alias << "', ignored\n";
      continue;
    }

    set_isotope(i, iso);

//  cerr << "ISIS_RXN_FILE_Molecule::atom " << i << " set to isotope " << iso << endl;

    rc++;
  }

  return rc;
}

void
RXN_File::discard_atom_map ()
{
  for (int i = 0; i < _nr; ++i)
  {
    _reagent[i].discard_atom_map();
  }
  for (int i = 0; i < _np; ++i)
  {
    _product[i].discard_atom_map();
  }

  return;
}

void
ISIS_RXN_FILE_Molecule::discard_atom_map()
{
  MDL_Molecule::discard_atom_map();

  if (NULL == _atom_map)
    return;

  std::fill_n(_atom_map, natoms(), 0);

  return;
}

int
ISIS_RXN_FILE_Molecule::assign_atom_types (Atom_Typing_Specification & ats)
{
  const int matoms = natoms();

  if (0 == matoms)
    return 1;

  if (nullptr == _atype)
    _atype = new_int(matoms);

  return ats.assign_atom_types(*this, _atype);
}

/*int
RXN_File::transfer_atom_aliases_to_isotopes ()
{
  int rc = 0;

  cerr << "RXN_File::transfer_atom_aliases_to_isotopes: " << _nr << " reagents\n";

  for (int i = 0; i < _nr; i++)
  {
    rc += _reagent[i].transfer_atom_alias_to_isotope();
  }

  for (int i = 0; i < _np; i++)
  {
    rc += _product[i].transfer_atom_alias_to_isotope();
  }

  return rc;
}*/

/*
  We need to look for mapped atoms that appear in product molecules, but which are
  not present in reagents
*/

//#define DEBUG_CHECK_FOR_WIDOWS_AND_ORPHANS

int
RXN_File::check_for_widows_and_orphans()
{
#ifdef DEBUG_CHECK_FOR_WIDOWS_AND_ORPHANS
  cerr << "RXN_File::_check_for_widows_and_orphans:checking " << _np << " product molecules\n";
#endif

  if (NULL == _reagent_locator || NULL == _product_locator)
    _setup_reagent_product_locator_arrays();

  for (int i = 0; i < _np; ++i)
  {
    const int matoms = _product[i].natoms();

#ifdef DEBUG_CHECK_FOR_WIDOWS_AND_ORPHANS
    cerr << "RXN_File::_check_for_widows_and_orphans:product " << i << " has " << matoms << " atoms\n";
#endif

    Set_of_Atoms orphans;

    for (int j = 0; j < matoms; ++j)
    {
      const int a = _product[i].atom_map(j);

      const int r = _reagent_locator[a];

#ifdef DEBUG_CHECK_FOR_WIDOWS_AND_ORPHANS
      cerr << "RXN_File::_check_for_widows_and_orphans:atom " << j << " map " << a << " in reagent " << r << endl;
#endif

      if (r < 0)
      {
#ifdef DEBUG_CHECK_FOR_WIDOWS_AND_ORPHANS
        cerr << "RXN_File::_check_for_widows_and_orphans:found orphan, atom " << j << " (map " << a << ") in product " << i << " type " << _product[i].smarts_equivalent_for_atom(j) << "\n";
#endif

        orphans.add(j);
      }
    }

#ifdef DEBUG_CHECK_FOR_WIDOWS_AND_ORPHANS
    cerr << "Found " << orphans.number_elements() << " orphan atoms\n";
#endif
    if (0 == orphans.number_elements())
      continue;

    if (orphans.number_elements() == matoms)    // every atom in this reagent is not on LHS, does not matter
      continue;

//  cerr << "_auto_fix_orphans " << _auto_fix_orphans << endl;

    if (! _auto_fix_orphans)
      return 0;

    if (! _fix_orphan_condition(orphans, i))
      return 0;
  }

  return 1;
}

//#define DEBUG_FIX_ORPHAN_CONDITION

int
RXN_File::_fix_orphan_condition (const Set_of_Atoms & orphans,
                                 const int p)
{
#ifdef DEBUG_FIX_ORPHAN_CONDITION
  cerr << "_fix_orphan_condition processing " << orphans << endl;
#endif

  const int norphan = orphans.number_elements();

  const int initial_orphan_natoms = _orphan_atoms.natoms();

  for (int i = 0; i < norphan; ++i)
  {
    const atom_number_t o = orphans[i];

    const Element * e = _product[p].elementi(o);
    _orphan_atoms.add(e);
    const auto q = _product[p].formal_charge(o);
    if (0 != q)
      _orphan_atoms.set_formal_charge(_orphan_atoms.natoms() - 1, q);
#ifdef DEBUG_FIX_ORPHAN_CONDITION
    cerr << "setting atom map for atom " << (_orphan_atoms.natoms() - 1) << " to " << _product[p].atom_map(o) << " product " << p << " orphan atom number " << o << endl;
#endif

    _orphan_atoms.set_atom_map(_orphan_atoms.natoms() - 1, _product[p].atom_map(o));
    _orphan_atoms.set_atom_map_number(_orphan_atoms.natoms() - 1, _product[p].atom_map(o));
  }

#ifdef DEBUG_FIX_ORPHAN_CONDITION
  cerr << "After adding orphan atom(s) " << _orphan_atoms.smiles() << endl;
#endif

  if (1 == norphan)    // no bonding to worry about
    return 1;

  for (int i = 0; i < norphan; ++i)
  {
    const atom_number_t oi = orphans[i];

    for (int j = i + 1; j < norphan; ++j)
    {
      const atom_number_t oj = orphans[j];

      if (! _product[p].are_bonded(oi, oj))
        continue;

      const Bond * b = _product[p].bond_between_atoms(oi, oj);
      _orphan_atoms.add_bond(initial_orphan_natoms + i, initial_orphan_natoms + j, (int)b->btype(), (int)b->btype(), 0);
    }
  }

  return 1;
}

static int 
bond_constant (const Bond * bondi)
{
  if (bondi->is_aromatic())
    return 11;
  if (bondi->is_triple_bond())
    return 7;
  if (bondi->is_double_bond())
    return 5;
  
  return 3;
}

static int
smallest_ring_containing (Molecule & m,
                          const atom_number_t zatom)
{
  const int nr = m.nrings();
  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (ri->contains(zatom))
      return ri->number_elements();
  }

  cerr << "smallest_ring_containing:ring atom not found in any ring!\n";

  return 0;    // should not happen
}

static int64_t
connected_atoms (Molecule & m,
                 const atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int64_t sum = 0;

  for (int i = 0; i < a->ncon(); ++i)
  {
    const Bond * b = a->item(i);
    const atom_number_t j = b->other(zatom);

    const int r = m.is_ring_atom(j);

    int smr;
    if (r)
      smr = smallest_ring_containing(m, j);
    else
      smr = 0;

    const Atom * aj = m.atomi(j);

//  initial attempt is too fine grained. Would catch atoms adjacent
//  const int t = 10000003 * aj->atomic_number() + 1000002 * aj->ncon() + 100007 * m.nbonds(j) +
//          40002 * r + 15007 * m.is_aromatic(j) + 2001 * m.hcount(j) + 581 * smr +
//          71 * aj->isotope() + bond_constant(b);

    const int t = 10000003 * aj->atomic_number() + 40002 * r + 15007 * m.is_aromatic(j) + 581 * smr +
            71 * aj->isotope() + bond_constant(b);

    sum += t;
  }

  return sum;
}

//#define DEBUG_ATOMS_THE_SAME

static int
atoms_the_same (ISIS_RXN_FILE_Molecule & m1,
                const atom_number_t x1,
                ISIS_RXN_FILE_Molecule & m2,
                const atom_number_t x2)
{
  const Atom * a1 = m1.atomi(x1);
  const Atom * a2 = m2.atomi(x2);

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare atomic numbers " << a1->atomic_number() << ' ' << a2->atomic_number() << endl;
#endif

  if (a1->atomic_number() != a2->atomic_number())
    return 0;

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare ncon " << a1->ncon() << ' ' << a2->ncon() << endl;
#endif

  if (a1->ncon() != a2->ncon())
    return 0;

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare nbonds " << m1.nbonds(x1) << ' ' << m2.nbonds(x2) << endl;
#endif

  if (m1.nbonds(x1) != m2.nbonds(x2))
    return 0;

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare hcount " << m1.hcount(x1) << ' ' << m2.hcount(x2) << endl;
#endif

  if (m1.hcount(x1) != m2.hcount(x2))
    return 0;

  const int r1 = m1.is_ring_atom(x1);
  const int r2 = m2.is_ring_atom(x2);

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare ring atom " << r1 << ' ' << r2 << endl;
#endif

  if (r1 != r2)
    return 0;

  if (r1)
  {
    const int a1 = m1.is_aromatic(x1);
    const int a2 = m2.is_aromatic(x2);

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare aromaticity " << a1 << ' ' << a2 << endl;
#endif

    if (a1 != a2)
      return 0;
  }

  if (a1->isotope() != a2->isotope())
    return 0;

  const auto t1 = connected_atoms(m1, x1);
  const auto t2 = connected_atoms(m2, x2);

#ifdef DEBUG_ATOMS_THE_SAME
  cerr << "Compare connected atoms " << t1 << ' ' << t2 << endl;
#endif

  if (t1 != t2)
    return 0;

  return 1;
}

int
RXN_File::identify_atoms_changing_reagent(const int rgnt,
                                          int * changed,
                                          const Changing_Atom_Conditions & cac)
{
//if (_nr > 1)
//  cerr << "RXN_File::identify_atoms_changing:warning, only examining first of " << _nr << " reagents\n";

  ISIS_RXN_FILE_Molecule & r = _reagent[rgnt];

  const int matoms = r.natoms();

  std::fill_n(changed, matoms, 0);

  const int * amap = r.atom_map();

  if (NULL == amap)
  {
    cerr << "RXN_File::identify_atoms_changing:no atom map data available\n";
    return 0;
  }

  if (NULL == _reagent_locator || NULL == _product_locator)
  {
    cerr << "RXN_File::identify_atoms_changing:reagent/product locator arrays not set\n";
    return 0;
  }

  int only_consider_fragment = -1;

  if (cac.only_consider_largest_reagent_fragment() && r.number_fragments() > 1)
    only_consider_fragment = r.largest_fragment();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == amap[i])
      continue;

    if (only_consider_fragment >= 0 && r.fragment_membership(i) != only_consider_fragment)
      continue;

    if (cac.ignore_lost_atom_if_isolated() && 0 == r.ncon(i))
      continue;

    const auto p = _product_locator[amap[i]];

    if (p < 0)   // atom is lost, so definitely changes
    {
      changed[i] = 1;
      rc++;
      continue;
    }

    const atom_number_t j = _product[p].which_is_mapped_atom(amap[i]);

    if (j < 0)    // should not happen
      continue;

#ifdef DEBUG_ATOMS_THE_SAME
    cerr << "Comparing atoms: Reagent atom number " << i << " " << r.smarts_equivalent_for_atom(i) << " atom RHS " << j << ' ' << _product[p].smarts_equivalent_for_atom(j) << " MAP " << amap[i] << endl;
#endif

    if (! atoms_the_same(r, i, _product[p], j))
    {
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "DIFFERENT\n";
#endif
      changed[i] = 1;
      rc++;
    }
    else if (cac.is_changing_if_different_neighbours() && _different_mapped_atoms_attached(r, i, amap[i]))
    {
      changed[i] = 1;
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "DIFFERENT - different atoms mapped, atom " << i << endl;
#endif
      rc++;
    }
  }

#ifdef DEBUG_ATOMS_THE_SAME
  for (int i = 0; i < matoms; ++i)
  {
    if (changed[i])
      cerr << " atom " << i << " (map" << r.atom_map(i) << ") changing\n";
  }

  cerr << "rc " << rc << endl;
#endif
  
  return rc + _identify_atoms_changing_by_bond_reagent(r, changed, cac);
}

int
RXN_File::identify_atoms_changing_product (const int prdct,
                                           int * changed,
                                           const Changing_Atom_Conditions & cac)
{
  ISIS_RXN_FILE_Molecule & p = _product[prdct];

  const int matoms = p.natoms();

  std::fill_n(changed, matoms, 0);

  const int * amap = p.atom_map();

  if (NULL == amap)
  {
    cerr << "RXN_File::identify_atoms_changing_product:no atom map data available\n";
    return 0;
  }

  if (NULL == _reagent_locator || NULL == _product_locator)
  {
    cerr << "RXN_File::identify_atoms_changing_product:reagent/product locator arrays not set\n";
    return 0;
  }

  int only_consider_fragment = -1;

  if (cac.only_consider_largest_reagent_fragment() && p.number_fragments() > 1)
    only_consider_fragment = p.largest_fragment();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == amap[i])
      continue;
 
    if (only_consider_fragment >= 0 && p.fragment_membership(i) != only_consider_fragment)
      continue;

    const auto r = _reagent_locator[amap[i]];

    if (r < 0)   // not present in reagents, we are an atom formed from nothing
    {
      changed[i] = 1;
      rc++;
      continue;
    }

    const atom_number_t j = _reagent[r].which_is_mapped_atom(amap[i]);

    if (j < 0)    // should not happen
      continue;

#ifdef DEBUG_ATOMS_THE_SAME
    cerr << "Comparing atoms mapped " << i << " atom " << i << " in RHS " << p.smarts_equivalent_for_atom(i) << " atom " << j << ' ' << _reagent[r].smarts_equivalent_for_atom(j) << endl;
#endif

    if (! atoms_the_same(p, i, _reagent[r], j))
    {
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "DIFFERENT\n";
#endif
      changed[i] = 1;
      rc++;
    }
  }

#ifdef DEBUG_ATOMS_THE_SAME
  for (int i = 0; i < matoms; ++i)
  {
    if (changed[i])
      cerr << " atom " << i << " (map " << p.atom_map(i) << ") changing\n";
  }

  cerr << " rc " << rc << endl;
#endif
  
  return 1;
}

/*
  Identify aotms that are aromatic, and where the different bonds might
  be just due to a different Kekule form
*/

int
RXN_File::_identify_just_changed_kekule_form(ISIS_RXN_FILE_Molecule & r,
                                             const int * changed,
                                             int * just_changed_kekule_form) const
{
  const int matoms = r.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (changed[i])
      continue;

    if (! r.is_aromatic(i))
      continue;

    just_changed_kekule_form[i] = _identify_just_changed_kekule_form(r, i);
  }

  return 1;
}

int
RXN_File::_identify_just_changed_kekule_form(ISIS_RXN_FILE_Molecule & r,
                                             const atom_number_t zatom) const
{
  const int * amap = r.atom_map();

  const int zmap = amap[zatom];

  const int pz = _product_locator[zmap];

  if (pz < 0)
    return 0;

  const atom_number_t a1 = _product[pz].which_is_mapped_atom(zmap);

  const Atom * a = r.atomi(zatom);

  const int acon = a->ncon();

  if (acon < 2 || acon > 3)    // hard to imagine either of these happening
    return 0;

  int rc = 0;

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    const atom_number_t j = b->other(zatom);

    const int amapj = amap[j];

    const int p2 = _product_locator[amapj];

    if (p2 != pz)    // atoms not in same product molecule, much larger change
      continue;

    const atom_number_t a2 = _product[p2].which_is_mapped_atom(amapj);

    const Bond * b2 = _product[p2].bond_between_atoms(a1, a2);

    if (b->btype() == b2->btype())    // no change in bond type
      continue;

    if (b->is_aromatic() && b2->is_aromatic())   // both aromatic, but different types, possible Kekule switch
      rc = 1;
    else
      return 0;   // non aromatic bond is different. This is more than a Kekule change
  }

  return rc;
}

//#define DEBUG_IDENTIFY_ATOMS_CHANGING_BY_BOND_REAGENT

int
RXN_File::_identify_atoms_changing_by_bond_reagent(ISIS_RXN_FILE_Molecule & r,
                                                   int * changed,
                                                   const Changing_Atom_Conditions & cac) const
{
#ifdef DEBUG_IDENTIFY_ATOMS_CHANGING_BY_BOND_REAGENT
  Molecule x(r);
  x.set_isotopes(changed);
  cerr << "RXN_File::_identify_atoms_changing_by_bond_reagent:initially " << x.smiles() << endl;
#endif

  if (cac.consider_aromatic_bonds())
  {
    r.compute_aromaticity_if_needed();
    for (int i = 0; i < _np; ++i)
    {
      _product[i].compute_aromaticity_if_needed();
    }
  }

  const int * amap = r.atom_map();

  const int nb = r.nedges();

  int rc = 0;

  int * just_changed_kekule_form = new_int(r.natoms()); std::unique_ptr<int[]> free_just_changed_kekule_form(just_changed_kekule_form);

//_identify_atoms_in_rings_separated_from_changing_atoms(r, changed, just_changed_kekule_form);

  _identify_just_changed_kekule_form(r, changed, just_changed_kekule_form);

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = r.bondi(i);

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    if (changed[a1] && changed[a2])
      continue;

// the following code, when not commented out, causes the two atoms to NOT be marked as changing.
// Down the line, this causes a core dump in retrosynthetic_quick because the bond between the atoms IS changing.

    if (!_mark_atoms_changed_when_kekule_form_of_bond_changes)
    {
      if (just_changed_kekule_form[a1] && just_changed_kekule_form[a2])
        continue;      
    }

    const int amap1 = amap[a1];
    const int amap2 = amap[a2];

    const int p1 = _product_locator[amap1];
    const int p2 = _product_locator[amap2];

    if (p1 < 0 || p2 < 0)
      continue;

#ifdef DEBUG_IDENTIFY_ATOMS_CHANGING_BY_BOND_REAGENT
    cerr << "amap1 " << amap1 << " amap2 " << amap2 << " in product " << p1 << " and " << p2 << " bond " << *b << endl;
#endif

    if (p1 != p2)    // sides of the bond end up in different fragments
      continue;

    const atom_number_t pa1 = _product[p1].which_is_mapped_atom(amap1);
    const atom_number_t pa2 = _product[p2].which_is_mapped_atom(amap2);

    const Bond * b2 = _product[p1].bond_between_atoms_if_present(pa1, pa2);
    if (NULL == b2)    // same product, but different fragment, definitely changed
    {
      if (0 == changed[a1])
      {
        changed[a1] = 2;
        rc++;
      }
      if (0 == changed[a2])
      {
        changed[a2] = 2;
        rc++;
      }
      continue;
    }

//  cerr << "atom numbers " << pa1 << " and " << pa2 << *b2 << endl;

//  if (cac.consider_aromatic_bonds() && b->is_aromatic() && b2->is_aromatic())
//    continue;

    if (b->is_aromatic() != b2->is_aromatic())    // if aromaticity changed, then definitely different
      ;
    else if (BOND_TYPE_ONLY(b->btype()) == BOND_TYPE_ONLY(b2->btype()))
      continue;

#ifdef DEBUG_IDENTIFY_ATOMS_CHANGING_BY_BOND_REAGENT
    cerr << " atom maps " << amap1 << " and " << amap2 << " marked as changing\n";
#endif

    if (0 == changed[a1])
    {
      changed[a1] = 2;     // note that atoms changed by their bond types changing get a different flag
      rc++;
    }
    if (0 == changed[a2])
    {
      changed[a2] = 2;
      rc++;
    }
  }

#ifdef DEBUG_IDENTIFY_ATOMS_CHANGING_BY_BOND_REAGENT
  x.set_isotopes(changed);
  cerr << "RXN_File::_identify_atoms_changing_by_bond_reagent:after " << x.smiles() << endl;
#endif

  return rc;
}

int
RXN_File::identify_atoms_changing_reagent (const int rgnt,
                                           Atom_Typing_Specification & ats,
                                           int * changed,
                                           const Changing_Atom_Conditions & cac)
{
//cerr << "LINE " << __LINE__ << " " << cac.is_changing_if_different_neighbours() << endl;
  ISIS_RXN_FILE_Molecule & r = _reagent[rgnt];

  const int matoms = r.natoms();

  std::fill_n(changed, matoms, 0);

  const int * amap = r.atom_map();

  if (NULL == amap)
  {
    cerr << "RXN_File::identify_atoms_changing_reagent:no atom map data available\n";
    return 0;
  }

  if (NULL == _reagent_locator || NULL == _product_locator)
  {
    cerr << "RXN_File::identify_atoms_changing_reagent:reagent/product locator arrays not set\n";
    return 0;
  }

  for (int i = 0; i < _nr; ++i)
  {
    _reagent[i].assign_atom_types(ats);
  }

  for (int i = 0; i < _np; ++i)
  {
    _product[i].assign_atom_types(ats);
  }

//cerr << " nf " << r.number_fragments() << " ignore " << cac.only_consider_largest_reagent_fragment() << endl;
  int only_consider_fragment = -1;
  if (cac.only_consider_largest_reagent_fragment() && r.number_fragments() > 1)
    only_consider_fragment = r.largest_fragment();

//cerr << "Largest frag " << only_consider_fragment << endl;

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
//  cerr << "Atom " << i << " mapped " << amap[i] << endl;
    if (0 == amap[i])
      continue;

//  cerr << " atom " << i << " type " << r.smarts_equivalent_for_atom(i) << " in frgament " << r.fragment_membership(i) << endl;
    if (only_consider_fragment >= 0 && r.fragment_membership(i) != only_consider_fragment)
      continue;

    const auto p = _product_locator[amap[i]];

    if (p < 0)   // not present in products, atom is disappearing
    {
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "Reagent atom " << i << ' ' << r.smarts_equivalent_for_atom(i) << " disappears\n";
#endif
      if (0 == r.ncon(i) && cac.ignore_lost_atom_if_isolated())    // we ignore disappearing atoms
        ;
      else
      {
        changed[i] = 1;
        rc++;
      }
      continue;
    }

    const atom_number_t j = _product[p].which_is_mapped_atom(amap[i]);

    if (j < 0)    // should not happen
      continue;

#ifdef DEBUG_ATOMS_THE_SAME
    cerr << "Comparing atoms: Reagent atom number " << i << " " << r.smarts_equivalent_for_atom(i) << " atom RHS " << j << ' ' << _product[p].smarts_equivalent_for_atom(j) << " MAP " << amap[i] << endl;
#endif

    if (_product[p].atom_type()[j] != r.atom_type()[i])
    {
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "DIFFERENT\n";
#endif
      changed[i] = 1;
      rc++;
    }
    else if (cac.is_changing_if_different_neighbours() && _different_mapped_atoms_attached(r, i, amap[i]))
    {
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "DIFFERENT:atom " << i << " (map " << amap[i] << ") has different mapped atoms attached\n";
#endif
      changed[i] = 1;
      rc++;
    }
  }

#ifdef DEBUG_ATOMS_THE_SAME
  Molecule mcopy(r);
  mcopy.set_isotopes(changed);
  cerr << "identify_atoms_changing_reagent:changing atoms " << mcopy.smiles() << endl;
  for (int i = 0; i < matoms; ++i)
  {
    if (changed[i])
      cerr << " atom " << i << " (map " << r.atom_map(i) << ") changing\n";
  }

  cerr << " rc " << rc << endl;
#endif

  if (r.nrings())
    _identify_changing_aromatic_systems(rgnt, changed, cac);
  
  if (cac.include_changing_bonds_in_changing_atom_count())
    rc += _identify_atoms_changing_by_bond_reagent(r, changed, cac);

  return rc;
}

//#define DEBUG_IDENTIFY_CHANGING_AROMATIC_SYSTEMS

static void
update_changed_array(const Ring & r,
                     int * changed,
                     const int flag)
{
  const int n = r.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const int j = r[i];

    if (0 == changed[j])
      changed[j] = flag;
  }

  return;
}

void
RXN_File::_update_aromatic_ring_system(ISIS_RXN_FILE_Molecule & m, 
                            const Ring & r,
                            int * changed,
                            const int flag,
                            int * ring_already_done) const
{
#ifdef DEBUG_IDENTIFY_CHANGING_AROMATIC_SYSTEMS
  cerr << "update_changed_array:ring has " << r.fused_ring_neighbours() << " fused neighbours\n";
#endif

  for (int i = 0; i < r.fused_ring_neighbours(); ++i)
  {
    const Ring * n = r.fused_neighbour(i);

    const int nrn = n->ring_number();

    if (ring_already_done[nrn])
      continue;

    ring_already_done[nrn] = 1;

    if (! n->is_aromatic())
      continue;

    if (! _all_bonds_unchanged(m, r))
      update_changed_array(*n, changed, flag);
  }

  return;
}

int
RXN_File::_identify_changing_aromatic_systems(const int reagent_number,
                                              int * changed,
                                              const Changing_Atom_Conditions & cac) const
{
  ISIS_RXN_FILE_Molecule & r = _reagent[reagent_number];

  const int nr = r.nrings();

  if (0 == nr)
    return 1;

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  r.compute_aromaticity_if_needed();

  for (int i = 0; i < nr; ++i)
  {
    if (ring_already_done[i])
      continue;

    ring_already_done[i] = 1;

    const Ring * ri = r.ringi(i);

    if (! ri->is_aromatic())
      continue;

#ifdef DEBUG_IDENTIFY_CHANGING_AROMATIC_SYSTEMS
    cerr << "aromatic ring " << *ri << " hits " << ri->number_members_non_zero_in_array(changed) << " changed atoms\n";
#endif

    if (ri->number_members_non_zero_in_array(changed) < 2)
      continue;

#ifdef DEBUG_IDENTIFY_CHANGING_AROMATIC_SYSTEMS
    cerr << "_all_bonds_unchanged " << _all_bonds_unchanged(r, *ri) << endl;
#endif

    if (_all_bonds_unchanged(r, *ri))   // nothing to worry about
      continue;

    update_changed_array(*ri, changed, 2);

    if (0 == ri->fused_ring_neighbours())
      continue;

    _update_aromatic_ring_system(r, *ri, changed, 2, ring_already_done);
  }

  return 1;
}

/*
  We are looking at aromatic bonds and want to figure out whether
  or not any of the bonds in a ring are changed in a reaction
*/

int
RXN_File::_all_bonds_unchanged(const ISIS_RXN_FILE_Molecule & m,
                               const Ring & r) const
{
  const int * amap = m.atom_map();

  for (Ring_Bond_Iterator i(r); i != r.end(); i++)
  {
    const atom_number_t ra1 = i.a1();
    const atom_number_t ra2 = i.a2();

    const int amap1 = amap[ra1];
    if (amap1 < 1)
      continue;

    const int amap2 = amap[ra2];
    if (amap2 < 1)
      continue;

    const int p1 = _product_locator[amap1];
    const int p2 = _product_locator[amap2];

    if (p1 < 0 || p1 != p2)
      return 0;

    const atom_number_t pa1 = _product[p1].which_is_mapped_atom(amap1);
    const atom_number_t pa2 = _product[p1].which_is_mapped_atom(amap2);

    const Bond * b1 = m.bond_between_atoms_if_present(ra1, ra2);
    if (NULL == b1)    // aromatic bond being destroyed
      continue;

    const Bond * b2 = _product[p1].bond_between_atoms_if_present(pa1, pa2);
    if (NULL == b2)
      continue;

    if (b1->btype() != b2->btype())
      return 0;     // some bonds are changed
  }

  return 1;    // all bonds unchanged
}


/*
  We have a ring that involves one or more changing atoms, CHANGED==1
  If the bonding is changed, is it just a different Kekule form?

  Oct 2016. Found this was problematic with

  [CH:3]1=[C:8]([CH3:10])[1C:7]2=[1C:6]([CH:5]=[C:4]1[Br:9])[1NH:12][1C:1](=[1N:11]2)[CH3:2]>>[O:13]([C:1]([CH3:2])([O:14][CH2:17][CH3:20])[O:15][CH2:18][CH3:21])[CH2:16][CH3:19]+[CH:3]1=[C:8]([CH3:10])[C:7]([NH2:11])=[C:6]([NH2:12])[CH:5]=[C:4]1[Br:9] 01K-E00124-063'

  Within the ring, it is just changes in Kekule form, but that change might be important for the reaction.

  So, until I can get a clearer view on this, I am just making this function a no-op
*/

#ifdef OLD_VERSION_NOW_I_HAVE_SOMETHING_BETTER
int
RXN_File::_identify_atoms_in_rings_separated_from_changing_atoms(ISIS_RXN_FILE_Molecule & reagent,
                                                                 const Ring & r,
                                                                 const int * changed,
                                                                 int * just_changed_kekule_form) const
{
  assert (r.is_aromatic());

  return 1;

  const int * amap = reagent.atom_map();

  int single_bond_count1 = 0;
  int double_bond_count1 = 0;
  int single_bond_count2 = 0;
  int double_bond_count2 = 0;

  int differences_found = 0;
  Set_of_Atoms differences_at_changing_atoms;

  int different_exocyclic_bonding = 0;

  const int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; ++i)
  {
    const atom_number_t j = r[i];

    const Atom * aj = reagent.atomi(j);

    const int jcon = aj->ncon();

    const int p = _product_locator[amap[j]];

    const atom_number_t pj = _product[p].which_is_mapped_atom(amap[j]);

    for (int k = 0; k < jcon; ++k)
    {
      const Bond * b = aj->item(k);

      const atom_number_t j2 = b->other(j);

      const int p2 = _product_locator[amap[j2]];

      if (p2 != p)     // atom in different product, cannot process here
        continue;

      const atom_number_t j2p = _product[p2].which_is_mapped_atom(amap[j2]);

      const Bond * b2 = _product[p2].bond_between_atoms(pj, j2p);
    
      if (BOND_TYPE_ONLY(b->btype()) == BOND_TYPE_ONLY(b2->btype()))
        continue;

      if (! r.contains(j2))   // exocyclic bond is different, tis is not just a Kekule thing
        return 1;
    }
  }

  r.set_vector(just_changed_kekule_form, 1);

  return 1;
}
#endif

#ifdef THIS_IS_NOT_IN_WORKING_CONDITION_____
int
RXN_File::_identify_atoms_in_rings_separated_from_changing_atoms(ISIS_RXN_FILE_Molecule & reagent,
                                                                 const Ring & r,
                                                                 const int * changed,
                                                                 int * just_changed_kekule_form) const
{
  no this is not right. We need to look at the actual bonding on both LHS and RHS. Cannot just look from LHS.

  assert (r.is_aromatic());

  const int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; ++i)
  {
    const atom_number_t j = r[i];

    const Atom * aj = reagent.atomi(j);

    const int jcon = aj->ncon();

    const int p = _product_locator[amap[j]];

    const atom_number_t pj = _product[p].which_is_mapped_atom(amap[j]);

    if (! _product[p].is_aromatic(pj))   // aromatic atom no longer aromatic, definitely not a Kekule thing
      return 0;

    for (int k = 0; k < jcon; ++k)     // examine connected atoms in LHS and RHS
    {
      const Bond * b = aj->item(k);

      const atom_number_t j2 = b->other(j);

      const int p2 = _product_locator[amap[j2]];

      if (p2 != p)    // atom in a different product. Presumably exocyclic stripped off
      {
        if (! b->is_single_bond())     // but currently connected via a double bond, definitely not Kekule form
          return 1;

        continue;     // single bond being broken, still OK to be Kekule variant
      }

      const atom_number_t j2p = _product[p2].which_is_mapped_atom(amap[j2]);

      const Bond * b2 = _product[p2].bond_between_atoms(pj, j2p);
    
      if (BOND_TYPE_ONLY(b->btype()) == BOND_TYPE_ONLY(b2->btype()))   // no change, great
        continue;

      if (! r.contains(j2))   // exocyclic bond is different. Cannot be Kekule variant
        return 1;
    }
  }

  r.set_vector(just_changed_kekule_form, 1);

  return 1;
}

#endif

int
ISIS_RXN_FILE_Molecule::gather_mapped_neighbours (const atom_number_t zatom,
                                                  resizable_array<int> & nbrs) const
{
  const Atom * a = atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);

    if (0 == _atom_map[j])
      continue;

    nbrs.add(_atom_map[j]);
  }

  return nbrs.number_elements();
}

int
ISIS_RXN_FILE_Molecule::number_mapped_atoms() const
{
  const int matoms = Molecule::natoms();

  if (NULL == _atom_map)
  {
    cerr << "ISIS_RXN_FILE_Molecule::number_mapped_atoms:no atom map\n";
    return 0;
  }

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (_atom_map[i] > 0)
      rc++;
  }

  return rc;
}

template <typename T>
void
transfer_array(T * & lhs, T * & rhs)
{
  if (NULL != lhs)
    delete [] lhs;

  lhs = rhs;
  rhs = NULL;

  return;
}

ISIS_RXN_FILE_Molecule &
ISIS_RXN_FILE_Molecule::operator= (ISIS_RXN_FILE_Molecule && rhs)
{
//cerr << "ISIS_RXN_FILE_Molecule::operator= move\n";

  transfer_array(_atom_map, rhs._atom_map);
  transfer_array(_explicit_hydrogen_count, rhs._explicit_hydrogen_count);
  transfer_array(_computed_implicit_hydrogen_count, rhs._computed_implicit_hydrogen_count);
  transfer_array(_explicit_hydrogen_atom_removed, rhs._explicit_hydrogen_atom_removed);
  transfer_array(_connections_lost, rhs._connections_lost);
  transfer_array(_btype, rhs._btype);
  transfer_array(_aromatic, rhs._aromatic);
  transfer_array(_nbonds, rhs._nbonds);
  transfer_array(_use_as_representative, rhs._use_as_representative);
  transfer_array(_atype, rhs._atype);

  _not_atom_lists_present = rhs._not_atom_lists_present;
  _convert_A_to_C_for_aromaticity = rhs._convert_A_to_C_for_aromaticity;
  _aromatic_bonds_lose_kekule_identity = rhs._aromatic_bonds_lose_kekule_identity;
  _swap_atoms_to_put_rare_atoms_first = rhs._swap_atoms_to_put_rare_atoms_first;
  _molecule_is_sole_reagent = rhs._molecule_is_sole_reagent;
  _remove_explicit_hydrogens = rhs._remove_explicit_hydrogens;
  _structure_value = rhs._structure_value;
  _is_duplicate_fragment = rhs._is_duplicate_fragment;
  _interpret_atom_alias_as_smarts = rhs._interpret_atom_alias_as_smarts;

  _unique_smiles = std::move(rhs._unique_smiles);
//_toggle_kekule_form = std::move(rhs._toggle_kekule_form);

  Molecule & mlhs = *this;
  Molecule & mrhs = rhs;
  mlhs = std::move(mrhs);

  MDL_File_Data & mdlfdlhs = *this;
  MDL_File_Data & mdlfdrhs = rhs;
  mdlfdlhs = std::move(mdlfdrhs);

  return *this;
}

//#define DEBUG_DIFFERENT_MAPPED_ATOMS_ATTACHED

/*
  We look at the mapped atoms attached to see if an atom in the reagents should be
  considered a changing atom
*/

int
RXN_File::_different_mapped_atoms_attached(const ISIS_RXN_FILE_Molecule & r,
                                           const atom_number_t zatom,
                                           const int zmap) const
{
  const auto p = _product_locator[zmap];

  const atom_number_t j = _product[p].which_is_mapped_atom(zmap);

#ifdef DEBUG_DIFFERENT_MAPPED_ATOMS_ATTACHED
  cerr << " atom " << zatom << " " << r.smarts_equivalent_for_atom(zatom) << " in reagents compare to " << j << ' ' << _product[p].smarts_equivalent_for_atom(j) << " in products\n";
#endif

  if (r.ncon(zatom) != _product[p].ncon(j))
    return 1;

  resizable_array<int> rmap, pmap;

  r.gather_mapped_neighbours(zatom, rmap);

  _product[p].gather_mapped_neighbours(j, pmap);

  const int nrmap = rmap.number_elements();

  if (nrmap != pmap.number_elements())
    return 1;

  if (0 == nrmap)   // no other mapped atoms attached on either side
    return 0;

#ifdef DEBUG_DIFFERENT_MAPPED_ATOMS_ATTACHED
  cerr << "From atom " << zatom << " map " << zmap << " which is atom " << j << " in product\n";
  cerr << r.smarts_equivalent_for_atom(zatom) << " and " << _product[p].smarts_equivalent_for_atom(j) << endl;
  cerr << "R";
  for (int i = 0; i < nrmap; ++i)
  {
    cerr << ' ' << rmap[i];
  }
  cerr << endl;
  cerr << "P";
  for (int i = 0; i < nrmap; ++i)
  {
    cerr << ' ' << pmap[i];
  }
  cerr << endl;
#endif

  if (1 == nrmap)
    return rmap[0] != pmap[0];

  if (2 == nrmap)
  {
    if (rmap[0] == pmap[0] && rmap[1] == pmap[1])
      return 0;
    if (rmap[0] == pmap[1] && rmap[1] == pmap[0])
      return 0;

    return 1;
  }

  Int_Comparator_Larger icl;
  rmap.iwqsort(icl);
  pmap.iwqsort(icl);

  for (int i = 0; i <nrmap; ++i)
  {
    if (rmap[i] != pmap[i])
      return 1;
  }

  return 0;
}

int
RXN_File::identify_atoms_changing_product (const int prdct,
                                           Atom_Typing_Specification & ats,
                                           int * changed,
                                           const Changing_Atom_Conditions & cac)
{
  ISIS_RXN_FILE_Molecule & p = _product[prdct];

  const int matoms = p.natoms();

  std::fill_n(changed, matoms, 0);

  const int * amap = p.atom_map();

  if (NULL == amap)
  {
    cerr << "RXN_File::identify_atoms_changing_product:no atom map data available\n";
    return 0;
  }

  if (NULL == _reagent_locator || NULL == _product_locator)
  {
    cerr << "RXN_File::identify_atoms_changing_product:reagent/product locator arrays not set\n";
    return 0;
  }

  for (int i = 0; i < _nr; ++i)
  {
    _reagent[i].assign_atom_types(ats);
  }

  for (int i = 0; i < _np; ++i)
  {
    _product[i].assign_atom_types(ats);
  }

  int only_consider_fragment = -1;

  if (cac.only_consider_largest_reagent_fragment() && p.number_fragments() > 1)
    only_consider_fragment = p.largest_fragment();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == amap[i])
      continue;

    if (only_consider_fragment >= 0 && p.fragment_membership(i) != only_consider_fragment)
      continue;

    const auto r = _reagent_locator[amap[i]];

    if (r < 0)   // not present in reagents, we are an atom formed from nothing
    {
      changed[i] = 1;
      rc++;
      continue;
    }

    const atom_number_t j = _reagent[r].which_is_mapped_atom(amap[i]);

    if (j < 0)    // should not happen
      continue;

#ifdef DEBUG_ATOMS_THE_SAME
    cerr << "Comparing atoms mapped " << i << " atom " << i << " in RHS " << p.smarts_equivalent_for_atom(i) << " atom " << j << ' ' << _reagent[r].smarts_equivalent_for_atom(j) << endl;
#endif

    if (_reagent[r].atom_type()[j] != p.atom_type()[i])
    {
#ifdef DEBUG_ATOMS_THE_SAME
      cerr << "DIFFERENT\n";
#endif
      changed[i] = 1;
      rc++;
    }
  }
  
  return rc;
}


#ifdef INTERESTING_IDEA_BUT_NOT_NEEDED
/*
  Does not work....

  We are looking for fragments that are really agents - they do not participate in the reaction.
  We look for their atom maps in the reagent, or product locator arrays. If all atoms are absent, then
  that fragment is an agent.

  Does not really need to be a member function
*/

int
RXN_File::_all_atom_maps_only_here (const ISIS_RXN_FILE_Molecule & m,
                                    const int * locator) const
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    const int a = m.atom_map(i);
    if (a < 1)
      continue;

    if (locator[a] >= 0)    // atom map appears on other side of reaction, not an agent
      return 0;
  }

  return 1;
}

/*
  Are there fragments that appear unchaged?
  When we are dealing with forward reactions, these will appear on the LHS, but when dealing with
  reversed reactions, they will appear as products.
  All their atom map numbers will appear in just that one place
*/

int
RXN_File::convert_unchanging_fragments_to_agents ()
{
  if (NULL == _atom_map)
  {
    cerr << "RXN_File::convert_unchanging_fragments_to_agents:no atom map present\n";
    return 0;
  }
  
  resizable_array<int> reagents_to_agents;

  for (int i = 0; i < _nr; ++i)
  {
    if (_all_atom_maps_only_here(_reagent[i], _product_locator))
      reagents_to_agents.add(i);
  }

  resizable_array<int> products_to_agents;

  for (int i = 0; i < _np; ++i)
  {
    if (_all_atom_maps_only_here(_product[i], _reagent_locator))
      products_to_agents.add(i);
  }

  const int rc = reagents_to_agents.number_elements() + products_to_agents.number_elements();

  if (0 == rc)
    return 0;

  if (_na > 0)
  {
    cerr << "RXN_File::convert_unchanging_fragments_to_agents:agents already present, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    return 0;
  }

  _na = rc;
  _agent = new ISIS_RXN_FILE_Molecule[_na];

  int ndx = 0;

  return rc;
}
#endif

Reaction_Subset::Reaction_Subset(const int nr) : _nr(nr)
{
  _xref = nullptr;

  _reagent_offset = new int[_nr];

  _include_atom = nullptr;

  _total_atoms = 0;

  return;
}

Reaction_Subset::~Reaction_Subset()
{
  if (nullptr != _xref)
    delete [] _xref;

  if (nullptr != _reagent_offset)
    delete [] _reagent_offset;

}

//#define DEBUG_REACTION_SUBSET_BUILD

int
Reaction_Subset::build(const ISIS_RXN_FILE_Molecule * reagent,
                       const int * reagent_locator,
                       const int * include_these_atom_map_numbers,
                       const int highest_atom_map)
{
  assert(nullptr != _reagent_offset);
#ifdef DEBUG_REACTION_SUBSET_BUILD
  cerr << "Reaction_Subset::build:building for reaction with " << _nr << " reagents\n";
#endif

  if (NULL == include_these_atom_map_numbers)    // nothing to do...
    return 1;

  _reagent_offset[0] = 0;

  _total_atoms = reagent[0].natoms();

  for (int i = 1; i < _nr; ++i)
  {
    _reagent_offset[i] = _reagent_offset[i-1] + reagent[i-1].natoms();

    _total_atoms += reagent[i].natoms();

    if (NULL == reagent[i].atom_map())
      cerr << "Reaction_Subset::build:reagent " << i << " has no atom map\n";
  }

#ifdef DEBUG_REACTION_SUBSET_BUILD
  cerr << "Total atoms " << _total_atoms << endl;
#endif

  _xref = new_int(_total_atoms + _total_atoms, -1);    // by default, nothing is included

  _include_atom = _xref + _total_atoms;    // avoid allocating two arrays

  int atoms_in_subset = 0;

  for (int i = 1; i <= highest_atom_map; ++i)
  {
#ifdef DEBUG_REACTION_SUBSET_BUILD
    cerr << "Atom map " << i << " inc " << include_these_atom_map_numbers[i] << endl;
#endif

    if (0 == include_these_atom_map_numbers[i])
      continue;

    const int r = reagent_locator[i];
    if (r < 0)    // may have been removed
      continue;

    const int a = reagent[r].which_is_mapped_atom(i);

#ifdef DEBUG_REACTION_SUBSET_BUILD
    cerr << " atom map " << i << " found in reagent " << r << " atom " << a << endl;
    cerr << " atom map " << i << " found in reagent " << r << " offset " << _reagent_offset[r] << endl;
#endif

    if (a < 0)    // might have been removed during fragment cleaning, so ignore. Might not be a safe thing to do...
    {
      cerr << "Reaction_Subset::build:no reagent atom is mapped atom " << i << endl;
      continue;
    }

#ifdef DEBUG_REACTION_SUBSET_BUILD
    cerr << " is atom number " << a << ' ' << reagent[r].smarts_equivalent_for_atom(a) << endl;
#endif

    _xref[_reagent_offset[r] + a] = 1;   // just mark with a positive value to show included. Will convert to xref later

    atoms_in_subset++;
  }

//cerr << atoms_in_subset << " atoms in subset\n";

  if (0 == atoms_in_subset)
  {
    cerr << "Reaction_Subset::build:no atoms in subset, highest atom map " << highest_atom_map << endl;
    return 0;
  }

  for (int i = 0; i < _nr; ++i)    // convert to cross reference values
  {
    int ndx = 0;
    const int matoms = reagent[i].natoms();
    int * xi = _xref + _reagent_offset[i];
    int * inc = _include_atom + _reagent_offset[i];

    for (int j = 0; j < matoms; ++j)
    {
      if (xi[j] < 0)
        inc[j] = 0;
      else
      {
        xi[j] = ndx;    // the new atom number for atom J in the subset
        ndx++;
        inc[j] = 1;
      }
    }
  }

#ifdef DEBUG_REACTION_SUBSET_BUILD
  debug_print(cerr);
#endif

  return 1;
}

atom_number_t
Reaction_Subset::atom_number_in_subset (const atom_number_t a, const int r) const
{
  if (nullptr == _xref)
    return a;

  const int * x = _xref + _reagent_offset[r];

  return x[a];
}

int
Reaction_Subset::debug_print(std::ostream & output) const
{
  output << "Reaction_Subset with " << _nr << " reagents\n";
  for (int i = 0; i < _nr; ++i)
  {
    const int * x = _xref + _reagent_offset[i];
    const int * inc = _include_atom + _reagent_offset[i];

    int jstart = _reagent_offset[i];
    int jstop;
    if (i == (_nr-1))
      jstop = _total_atoms;
    else
      jstop = _reagent_offset[i+1];

    output << " reagent " << i << " has " << (jstop - jstart) << " atoms\n";
    for (int j = jstart; j < jstop; ++j)
    {
      output << "   " << j << ' ' << x[j] << " inc " << inc[j] << '\n';
    }
  }

  return 1;
}

void
Reaction_Subset::convert_cross_reference_to_boolean_include()
{
  for (int i = 0; i < _nr; ++i)
  {
    int jstart = _reagent_offset[i];
    int jstop;
    if (i == (_nr-1))
      jstop = _total_atoms;
    else
      jstop = _reagent_offset[i+1];

    for (int j = jstart; j < jstop; ++j)
    {
      if (_xref[j] >= 0)
        _xref[j] = 1;
      else
        _xref[j] = 0;
    }
  }

  return;
}

//#define DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING

int
RXN_File::remove_fragments_not_participating ()
{
#ifdef DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING
  cerr << "RXN_File::remove_fragments_not_participating, _nr " << _nr << endl;
#endif

  if (NULL == _reagent_locator)
    prepare_for_reaction_construction();

  int max_atoms_in_any_reagent = _reagent[0].natoms();
  int max_frag = _reagent[0].number_fragments();

  for (int i = 1; i < _nr; ++i)
  {
    if (_reagent[i].natoms() > max_atoms_in_any_reagent)
      max_atoms_in_any_reagent = _reagent[i].natoms();
    if (_reagent[i].number_fragments() > max_frag)
      max_frag = _reagent[i].number_fragments();
  }

#ifdef DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING
  cerr << "Max frag " << max_frag << endl;
#endif

//  if (1 == max_frag)    // all reagents have one fragment, cannot remove anything
//    return 0;

  int * atoms_appearing_in_products = new int[max_frag + max_frag + max_atoms_in_any_reagent]; std::unique_ptr<int[]> free_atoms_appearing_in_products(atoms_appearing_in_products);
  int * to_remove = atoms_appearing_in_products + max_frag + max_frag;

  int rc = 0;

  for (int i = 0; i < _nr; ++i)
  {
    const int nf = _reagent[i].number_fragments();

//    if (1 == nf)
//      continue;

    std::fill_n(atoms_appearing_in_products, max_frag + max_frag, 0);
    int * atoms_disappearing = atoms_appearing_in_products + max_frag;

    const int matoms = _reagent[i].natoms();

#ifdef DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING
    cerr << "remove_fragments_not_participating:reagent " << i << " has " << matoms << " atoms\n";
#endif

    for (int j = 0; j < matoms; ++j)
    {
      const int amap = _reagent[i].atom_map(j);
      if (amap < 0)
        continue;

      const int f = _reagent[i].fragment_membership(j);

//    cerr << "RXN_File::remove_fragments_not_participating:reagent " << i << " atom " << j << " map " << amap << " in product " << _product_locator[amap] << endl;

      if (_product_locator[amap] >= 0)
        atoms_appearing_in_products[f]++;
      else
        atoms_disappearing[f]++;
    }

    if (nf == count_non_zero_occurrences_in_array(atoms_appearing_in_products, nf))  // every fragment has one or more atoms appearing in the products
      continue;

    std::fill_n(to_remove, matoms, 0);

    for (int f = 0; f < nf; ++f)
    {
      if (atoms_appearing_in_products[f])    // must keep this fragment
        continue;

      if (0 == atoms_disappearing[f])     // no atoms appearing in the product, but nothing lost either, who cares - seems unlikely
        continue;

//    At this stage, we have a fragment with no atoms on the RHS, but mapped atoms that disappear

#ifdef DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING
      cerr << " removing fragment " << f << " from reagent " << i << " atoms_disappearing " << atoms_disappearing[f] << endl;
#endif

      for (int j = 0; j < matoms; ++j)
      {
        if (f == _reagent[i].fragment_membership(j))
          to_remove[j] = 1;
      }

#ifdef DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING
      cerr << "Reagent " << i << " contains " << _reagent[i].natoms() << " atoms before removal\n";
      for (int q = 0; q < matoms; ++q)
      {
        cerr << q << " to_remove " << to_remove[q];
        if (to_remove[q])
          cerr << ' ' << _reagent[i].smarts_equivalent_for_atom(q);
        cerr << endl;
      }
#endif
    }

    _reagent[i].remove_atoms(to_remove);
    _reagent[i].recompute_unique_smiles();
    rc++;
  }

#ifdef DEBUG_REMOVE_FRAGMENTS_NOT_PARTICIPATING
  for (int i = 0; i < _nr; ++i)
  {
    cerr << "Reagent " << i << " " << _reagent[i].unique_smiles() << " " << _reagent[i].natoms() << " atoms\n";
  }
#endif

  int ndx = 0;

  for (int i = 0; i < _nr; ++i)
  {
    if (0 == _reagent[i].natoms())
      continue;

    if (ndx != i)
      _reagent[ndx] = std::move(_reagent[i]);
    ndx++;
  }

  _nr = ndx;

  if (NULL != _reagent_locator)
    _fill_reagent_and_product_locator_arrays();

  return rc;
}

/*
  We can have reagents where all the mapped atoms do not appear on the RHS.
  Eliminate them.
*/

int
RXN_File::eliminate_reagents_not_participating()
{
  cerr << "RXN_File::eliminate_reagents_not_participating, nr " << _nr << endl;

  if (1 == _nr)   // nothing to do
    return 0;

  if (NULL == _reagent_locator)
    prepare_for_reaction_construction();

  int * remove_reagent = new_int(_nr); std::unique_ptr<int[]> free_remove_reagent(remove_reagent);
  int nremove = 0;
  int mapped_atoms = 0;

  for (int i = 0; i < _nr; ++i)
  {
    const int matoms = _reagent[i].natoms();

    int atoms_in_rhs = 0;

    for (int j = 0; j < matoms; ++j)
    {
      const int amap = _reagent[i].atom_map(j);

      if (amap <= 0)
        continue;

      mapped_atoms++;

      const int p = _product_locator[amap];
      if (p < 0)
        continue;

      atoms_in_rhs = 1;
      break;
    }

    //cerr << "Reagent " << i << " with " << matoms << " atoms has " << atoms_in_rhs << " atoms appearing as products\n";

    if (0 == mapped_atoms)
      continue;

    if (0 == atoms_in_rhs)
    {
      remove_reagent[i] = 1;
      nremove++;
    }
  }

  if (0 == nremove)
    return 0;

  if (_nr == nremove)
  {
    cerr << "RXN_File::eliminate_reagents_not_participating:all reagents slated for removal, cannot do\n";
    return 0;
  }

  int ndx = 0;
  for (int i = 0; i < _nr; ++i)
  {
    if (remove_reagent[i])
      continue;

//  if (ndx != i)
//    cerr << "Moving " << i << " to " << ndx << endl;
    if (ndx != i)
      _reagent[ndx] = std::move(_reagent[i]);

//  cerr << "new molecule smiles " << _reagent[ndx].smiles() << endl;
    ndx++;
  }

  _nr = ndx;

  _reestablish_reagent_locator_array();

  return 1;
}

/*
  This is specifically designed for reversed reactions.
  Look through the products for fragments that do not appear as reagents.
  Essentially we are reversing the process of adding that we probably did earlier...
*/

int
RXN_File::remove_non_participating_fragments()
{
  if (NULL == _reagent_locator)
    prepare_for_reaction_construction();

  int rc = 0;

  for (int i = 0; i < _np; ++i)
  {
    auto & pi = _product[i];

    const int nf = pi.number_fragments();
//  cerr << "Product " << i << " has " << nf << " fragments\n";

    if (1 == nf)
      continue;

    if (_remove_non_participating_fragments(pi))
      rc++;
  }

  return rc;
}

//#define DEBUG_REMOVE_NON_PARTICIPATING_FRAGMENTS

int
RXN_File::_remove_non_participating_fragments(ISIS_RXN_FILE_Molecule & p)   // one of the products
{
  const int nf = p.number_fragments();

  const int matoms = p.natoms();

#ifdef DEBUG_REMOVE_NON_PARTICIPATING_FRAGMENTS
  cerr << "RXN_File::_remove_non_participating_fragments:molecule with " << matoms << " has " << nf << " fragments\n";
#endif

  int * remove = new_int(matoms + matoms); std::unique_ptr<int[]> free_remove(remove);
  int * fragment_membership = remove + matoms;
  p.fragment_membership(fragment_membership);

  int rc = 0;

  for (int f = 0; f < nf; ++f)
  {
    int all_atoms_absent = 1;

    for (int i = 0; i < matoms; ++i)
    {
      if (fragment_membership[i] != f)
        continue;

      const int amap = p.atom_map(i);
#ifdef DEBUG_REMOVE_NON_PARTICIPATING_FRAGMENTS
      cerr << "Fragment " << f << " atom " << i << ' ' << p.smarts_equivalent_for_atom(i) << " map " << amap << " in reagent " << _reagent_locator[amap] << endl;
#endif

      if (_reagent_locator[amap] >= 0)
      {
        all_atoms_absent = 0;    // atom is in a reagent, we cnanot remove this fragment
        break;
      }
    }

    if (0 == all_atoms_absent)
      continue;

#ifdef DEBUG_REMOVE_NON_PARTICIPATING_FRAGMENTS
    cerr << "Will be removing all the atoms in fragment " << f << endl;
#endif

    for (int i = 0; i < matoms; ++i)
    {
      if (fragment_membership[i] == f)
        remove[i] = 1;
    }
    rc++;
  }

  if (0 == rc)
    return 0;

//do_write("B4.rxn");
  p.remove_atoms(remove);

  _reestablish_locator_array(_product, _np, _product_locator);
//do_write("AFTER.rxn");

  return rc;
}




int
RXN_File::_reestablish_reagent_locator_array()
{
	if (NULL == _reagent_locator)
    prepare_for_reaction_construction();

#ifdef DEBUG_REESTABLISH_REAGENT_LOCATOR_ARRAY
  cerr << "RXN_File::_reestablish_reagent_locator_array:have " << _nr << " reagents\n";
#endif

  const int h = highest_atom_map_number();

  std::fill_n(_reagent_locator, h+1, -1);

  for (int i = 0; i < _nr; ++i)
  {
    const int matoms = _reagent[i].natoms();

#ifdef DEBUG_REESTABLISH_REAGENT_LOCATOR_ARRAY
    cerr << "RXN_File::_reestablish_reagent_locator_array: processing " << _reagent[i].smiles() << endl;
#endif

    for (int j = 0; j < matoms; ++j)
    {
      const int amap = _reagent[i].atom_map(j);

      if (0 == amap)
        continue;

      _reagent_locator[amap] = i;

#ifdef DEBUG_REESTABLISH_REAGENT_LOCATOR_ARRAY
      cerr << "RXN_File::_reestablish_reagent_locator_array: reagent " << i << " atom " << j << " map " << amap << endl;
#endif
    }
  }

  return 1;
}

int
RXN_File::_reestablish_locator_array(ISIS_RXN_FILE_Molecule * x, const int n,
                                     int * locator)
{
  assert (NULL != locator);

  for (int i = 0; i < n; ++i)
  {
    const int matoms = x[i].natoms();

#ifdef DEBUG_REESTABLISH_REAGENT_LOCATOR_ARRAY
    cerr << "RXN_File::_reestablish_locator_array: processing " << x[i].smiles() << endl;
#endif

    for (int j = 0; j < matoms; ++j)
    {
      const int amap = x[i].atom_map(j);

      if (0 == amap)
        continue;

      locator[amap] = i;

#ifdef DEBUG_REESTABLISH_REAGENT_LOCATOR_ARRAY
      cerr << "RXN_File::_reestablish_locator_array: r/p " << i << " atom " << j << " map " << amap << endl;
#endif
    }
  }

  return 1;
}

/*
  Note that we do not check chirality here
*/

int
RXN_File::all_reagents_the_same ()
{
  if (_nr <= 1)
    return 0;

  const int matoms = _reagent[0].natoms();
  for (int i = 1; i < _nr; ++i)
  {
    if (matoms != _reagent[i].natoms())
      return 0;
  }

// At this stage, they all have the same atom count

//cerr << __LINE__ << endl;

  Int_Comparator_Larger icl;

  for (int i = 0; i < matoms; ++i)
  {
    const atomic_number_t z0 = _reagent[0].atomic_number(i);
    const int amap0 = _reagent[0].atom_map(i);

    resizable_array<int> i_mapped_nbrs;
    _reagent[0].identify_connected_mapped_atoms(i, i_mapped_nbrs);
    i_mapped_nbrs.iwqsort(icl);

//  cerr << "Atom " << i << " Z " << z0 << " map " << amap0 << endl;

    for (int j = 1; j < _nr; ++j)
    {
      const int a2 = _reagent[j].which_is_mapped_atom(amap0);
      if (a2 < 0)
        return 0;

      if (z0 != _reagent[j].atomic_number(a2))
        return 0;

      resizable_array<int> a2_mapped_nbrs;
      _reagent[j].identify_connected_mapped_atoms(a2, a2_mapped_nbrs);
      a2_mapped_nbrs.iwqsort(icl);

      if (! (i_mapped_nbrs == a2_mapped_nbrs))    // there is no operator != defined
        return 0;
    }
  }

// they are the same, forget about the dups

  _nr = 1;

  if (NULL == _reagent_locator)
    return 1;

  return _reestablish_reagent_locator_array();
}

static int
atom_map_numbers_the_same(const ISIS_RXN_FILE_Molecule & r1,
                          const ISIS_RXN_FILE_Molecule & r2)
{
  extending_resizable_array<int> map1;

  const int matoms = r1.natoms();

  map1.resize(2 * matoms);

  assert (matoms == r2.natoms());

  for (int i = 0; i < matoms; ++i)
  {
    const int amap = r1.atom_map(i);

    map1[amap] = 1;
  }

  for (int i = 0; i < matoms; ++i)
  {
    const int amap = r2.atom_map(i);

    if (0 == map1[amap])
      return 0;
  }

  return 1;
}

//reduce the product to the largest fragment in any component

int
RXN_File::reduce_to_largest_product()
{ 
	return reduce_to_largest_component(_product,_np);
}       
int

RXN_File::reduce_to_largest_reactant()
{ 
	return reduce_to_largest_component(_reagent,_nr);
}       
	
int
RXN_File::reduce_to_largest_component(ISIS_RXN_FILE_Molecule *component, int &n)
{        
	 
  int largestFragmentComponentIndex, largestFragmentFragmentindex, largestFragmentSize = 0;
  
  for (int i = 0; i < n; ++i)
  {
    auto & pi = component[i];

    const int nf = pi.number_fragments();
    
    for (int j=0;  j != nf ; ++j)
    {
    	int fragAtoms = pi.atoms_in_fragment (j);
    	if (fragAtoms > largestFragmentSize)
    	{
    		largestFragmentSize = fragAtoms;
    		largestFragmentComponentIndex = i;
    		largestFragmentFragmentindex = j;
    	}
    }
  }
    
  if (largestFragmentSize == 0)
  	return 0;  // nothing to remove
    	
   if (largestFragmentComponentIndex != 0)
   	 component[0] = std::move(component[largestFragmentComponentIndex]); 	
   n = 1;
	
	// now we have just the component that has the largest fragment, so just trim it to its largest fragment!
	
	component[0].reduce_to_largest_fragment();
	
  _reestablish_reagent_locator_array();
	
  return 1;
}

/*
  We may have duplicate reagents, but the atom maps are scrambled.
  The unique smiles must match, AND they must have the same mapped atoms
*/

int
RXN_File::remove_duplicate_reagents_atom_maps_scrambled ()
{
  int rc = 0;

// We exclude chirality in these comparisons, but it is reasonable to ask if that is correct or not....

  set_include_chiral_info_in_smiles(0);
  for (int i = 0; i < _nr; ++i)
  {
    _reagent[i].unique_smiles();
  }
  set_include_chiral_info_in_smiles(1);

  for (int i = 0; i < _nr; ++i)
  {
    const IWString & iusmi = _reagent[i].unique_smiles();

    const int matomsi = _reagent[i].natoms();

//  cerr << "Reagent " << i << " usmi " << iusmi << endl;

    for (int j = (i+1); j < _nr; ++j)
    {
//    cerr << "  reagent " << j << " usmi " << _reagent[j].unique_smiles() << endl;

      if (_reagent[j].natoms() != matomsi)
        continue;

      if (_reagent[j].unique_smiles() != iusmi)   // definitely not the same
        continue;

      if (! atom_map_numbers_the_same(_reagent[i], _reagent[j]))
        continue;

//    cerr << "in " << name() << " duplicate reagents i = " << i << " j = " << j << endl;

      for (int k = j; k < (_nr-1); ++k)     // reagent J is a duplicate of I. Shift everything down
      {
        _reagent[k] = std::move(_reagent[k+1]);
      }
      _nr--;
      j--;
      rc++;
    }
  }

  if (0 == rc)
    return 0;

  return _reestablish_reagent_locator_array();
}

int
RXN_File::_remove_unmapped_components (ISIS_RXN_FILE_Molecule * component, int & n)
{
  int rc = 0;
  
  for (int i = 0; i < n; ++i)
  {
  	const int matoms = component[i].natoms();
		int foundAMapping = 0;


	  for (int j = 0; j < matoms; ++j)
	  {
	    if (component[i].atom_map()[j] > 0)
	    {   	
	    	foundAMapping = 1;
	    	break;
	    }
	  }
	  
	  if (foundAMapping)
	        continue;


    for (int k = i; k < (n-1); ++k)     // reagent J is a duplicate of I. Shift everything down
    {
      component[k] = std::move(component[k+1]);
    }
    n--;
    i--;
    rc++;
  }

  if (0 == rc)
    return 0;

  return _reestablish_reagent_locator_array();
}

int
RXN_File::remove_unmapped_components ()
{
	int rc = 1;
	if (!_remove_unmapped_components(_reagent, _nr))
  	rc = 0;
	if (!_remove_unmapped_components(_agent, _na))
  	rc = 0;
	if (!_remove_unmapped_components(_product, _np))
  	rc = 0;

  return rc;
}
