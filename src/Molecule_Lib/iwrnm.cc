#include <iostream>
#include "Foundational/iwmisc/misc.h"
#include "smiles.h"
#include "iwrnm.h"
#include "misc2.h"

using std::cerr;
using std::endl;

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

static int FileScope_include_directionality_in_ring_closure_bonds = 0;

int
include_directionality_in_ring_closure_bonds()
{
  return FileScope_include_directionality_in_ring_closure_bonds;
}

void
set_include_directionality_in_ring_closure_bonds (int s)
{
  FileScope_include_directionality_in_ring_closure_bonds = s;
}

void
Ring_Number_Manager::_default_values()
{
  _nr = 0;
  _ring_id = nullptr;
  _bond = nullptr;
//_bt = nullptr;
  _from_atom = nullptr;
  _include_aromaticity_in_smiles = get_include_aromaticity_in_smiles();
  _include_cis_trans_in_smiles = include_cis_trans_in_smiles();

  return;
}

Ring_Number_Manager::Ring_Number_Manager()
{
  _default_values();
}

Ring_Number_Manager::Ring_Number_Manager (int nrings)
{
  _default_values();

  if (0 == nrings)
  {
    return;
  }

  activate(nrings + 1);      // we don't use ring number 0, so allocate space for one extra

  return;
}

#define UNUSED_RING_NUMBER -1
#define NEVER_USE_THIS_RING_AGAIN -2

int
Ring_Number_Manager::ok() const
{
  if (_nr < 0)
    return 0;

  if (_nr > 10000)
    cerr << "Ring_Number_Manager::ok:improbably large nrings value " << _nr << endl;

  if (nullptr != _ring_id)
  {
    for (int i = 1; i <= _nr; i++)
    {
      if (UNUSED_RING_NUMBER == _ring_id[i])
        continue;

      if (NEVER_USE_THIS_RING_AGAIN == _ring_id[i])
        continue;

      if (_ring_id[i] < 0)
      {
        cerr << "Ring_Number_Manager::ok:ring id " << i << " is " << _ring_id[i] << endl;
        return 0;
      }

      if (_ring_id[i] > 100000)
        cerr << "Ring_Number_Manager::ok:improbably large ring number " << _ring_id[i] << endl;
    }
  }

  if (nullptr != _from_atom)
  {
    for (int i = 1; i <= _nr; i++)
    {
      if (_ring_id[i] < 0)
        continue;

      if (_from_atom[i] < 0)
      {
        cerr << "Ring_Number_Manager::ok:from atom " << i << " is " << _from_atom[i] << endl;
        return 0;
      }
    }
  }

  if (nullptr != _bond)
  {
    for (int i = 1; i < _nr; i++)
    {
      if (_ring_id[i] < 0)
        continue;

      if (! _bond[i]->ok())
      {
        cerr << "Ring_Number_Manager::ok:bond " << i << " is bad\n";
        return 0;
      }
    }
  }

  return 1;
}

int
Ring_Number_Manager::activate (int nrings)
{
  assert (nrings > 0);
  assert (nullptr == _ring_id);

  _nr = nrings;
  _ring_id = new_int(_nr + 1, UNUSED_RING_NUMBER);
  _bond    = new const Bond * [_nr + 1];
  _from_atom = new atom_number_t[_nr + 1];

  _ring_id[0] = -999;    // ring 0 is never used for smiles

  return 1;
}

Ring_Number_Manager::~Ring_Number_Manager()
{
  if (nullptr != _ring_id)
  {
    delete [] _ring_id;
    delete [] _bond;
    delete []_from_atom;
  }
}

int
Ring_Number_Manager::debug_print (std::ostream & os) const
{
  os << "Info on Ring_Number_Manager for " << _nr << " rings\n";

  for (int i = 1; i <= _nr; i++)
  {
    if (UNUSED_RING_NUMBER == _ring_id[i])
    {
      os << "ring " << i << " not active\n";
      continue;
    }

    if (NEVER_USE_THIS_RING_AGAIN == _ring_id[i])
    {
      os << "ring " << i << " never to be re-used\n";
      continue;
    }

    os << "Ring " << i << " atom " << _ring_id[i] << " bt = " <<
          _bond[i]->btype() << " other = " << _from_atom[i] << endl;
  }

  return os.good();
}

/*static void
append_bond_type (IWString & smiles, const Bond * b)
{
  if (b->is_single_bond())
    ;
  else if (b->is_double_bond())
    smiles += DOUBLE_BOND_SYMBOL;
  else if (b->is_triple_bond())
    smiles += TRIPLE_BOND_SYMBOL;
  else if (b->is_aromatic())
    smiles += AROMATIC_BOND_SYMBOL;
  else if (IS_NOT_A_BOND (b->btype()))
    smiles += '.';

//  all other bond types are silently ignored

  return;
}*/

/*
  Jul 2003. John Lampe had a case where we may not want to use low ring numbers
  Allow an offset
*/

static int ring_number_offset = 0;

void
set_smiles_ring_number_offset (int s)
{
  ring_number_offset = s;

  return;
}

static char single_digit[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

//#define DEBUG_APPEND_RING_CLOSURE_DIGITS

void
Ring_Number_Manager::_append_ring_closure_digits (IWString & smiles,
                            int ring_closure_number,
                            const Bond * b,
                            atom_number_t ato) const
{
#ifdef DEBUG_APPEND_RING_CLOSURE_DIGITS
  if (nullptr == b)
    cerr << "No ring closure bond\n";
  else if (b->is_aromatic())
    cerr << "closed with an aromatic bond\n";
  else if (b->is_single_bond())
    cerr << "closed with a single bond\n";
  else if (! b->is_single_bond())
    cerr << "Non single bond found\n";
#endif

  ring_closure_number += ring_number_offset;

/*if (nullptr == b)
    ;
  else if (b->is_aromatic() && include_aromaticity_in_smiles)
    ;
  else if (b->is_permanent_aromatic())   // never write these
    ;
  else if (! b->is_single_bond())
    append_bond_type (smiles, b);*/

// Oct 2009. Do not include directional bonds when closing rings - corina hates it
//  C1(=C/N(C)C)/C(=O)NN=C1 PBCHM737863           OK
//  C1(=C/N(C)C)/C(=O)NN=C/1 PBCHM737863          Corina complains incomplete specification and ignores
//  C1(=C/N(C)C)/C(=O)NN=C\1 PBCHM737863          Corina complains inconsistent and fails

  if (nullptr == b)   // is a ring opening
    ;
  else if (! b->is_single_bond())
    b->append_bond_type(smiles, ato, _include_aromaticity_in_smiles);
  else if (b->is_aromatic())
    b->append_bond_type(smiles, ato, _include_aromaticity_in_smiles);
  else if (! b->is_directional())
    b->append_bond_type(smiles, ato, _include_aromaticity_in_smiles);
  else if (! _include_cis_trans_in_smiles)
    b->append_bond_type(smiles, ato, _include_aromaticity_in_smiles);
  else if (FileScope_include_directionality_in_ring_closure_bonds)
    b->append_bond_type(smiles, ato, _include_aromaticity_in_smiles);
  else    // is directional, but don't want it to be marked directional
  {
    set_include_cis_trans_in_smiles(0);
    b->append_bond_type(smiles, ato, _include_aromaticity_in_smiles);
    set_include_cis_trans_in_smiles(1);
  }

  if (ring_closure_number < 10)
    smiles += single_digit[ring_closure_number];
  else if (ring_closure_number < 100)
  {
    smiles += '%';
    smiles.append_number(ring_closure_number);
  }
  else
  {
    smiles << "%(";
    smiles.append_number(ring_closure_number);
    smiles << ')';
  }

  return;
}

//#define DEBUG_STORE_RING

/*
  A ring opening digit is being created from atom AFROM to atom ATO

*/


int
Ring_Number_Manager::_generate_ring_opening_chars (IWString & ring_opening_chars,
                                                  const resizable_array<const Bond *> & ring_opening_bonds,
                                                  atom_number_t ato)
{
  const int nro = ring_opening_bonds.number_elements();

  ring_opening_chars.resize(2 * nro);

  for (int i = 0; i < nro; ++i)
  {
    const Bond * b = ring_opening_bonds[i];

    const atom_number_t afrom = b->other(ato);

    int free_ring = locate_item_in_array(UNUSED_RING_NUMBER, _nr + 1, _ring_id);
    assert (free_ring > 0);

    _ring_id[free_ring] = afrom;
    _from_atom[free_ring] = ato;
    _bond[free_ring] = b;

    _append_ring_closure_digits(ring_opening_chars, free_ring, nullptr, ato);
  }

  return 1;
}

/*
  Handle the Concord required ring opening/closing convention as appropriate
*/

int
Ring_Number_Manager::_process_ring (IWString & smiles,
                                    int ring_number,
                                    atom_number_t ato)
{
  assert (ring_number >= 0 && ring_number < _nr);

  _append_ring_closure_digits(smiles, ring_number, _bond[ring_number], ato);

  if (smiles_reuse_ring_closure_numbers())
    _ring_id[ring_number] = UNUSED_RING_NUMBER;
  else
    _ring_id[ring_number] = NEVER_USE_THIS_RING_AGAIN;

  return 1;
}

int
Ring_Number_Manager::_place_ring_closure (IWString & smiles,
                   atom_number_t a,
                   atom_number_t afrom)
{
  for (int i = 0; i < _nr; i++)
  {
    if (a == _ring_id[i] && afrom == _from_atom[i])
    {
      _process_ring(smiles, i, a);
      return 1;
    }
  }

  cerr << "Ring_Number_Manager::_place_ring_closure: no ring closure from " << a << " to " << afrom << endl;
  debug_print(cerr);
  iwabort();

  return 0;
}

/*
  Atom A is a chiral centre, and some ring closures end on it.
  We must place the ring closure numbers (bonds) in the correct chiral
  order. The calling routine must have placed RING_CLOSURES_FOUND in
  the correct order for the smiles
*/

int
Ring_Number_Manager::_append_ring_closures_for_chiral_atom (IWString & smiles,
                           atom_number_t a,
                           const resizable_array<atom_number_t> & ring_closures_found)
{
//if (2 != ring_closures_found.number_elements())
//{
//  cerr << "Ring_Number_Manager::_append_ring_closures_for_chiral_atom: problem at atom " << a << endl;
//  cerr << "Found " << ring_closures_found.number_elements() << " ring closures:";
//  for (int i = 0; i < ring_closures_found.number_elements(); i++)
//    cerr << ' ' << ring_closures_found[i];
//  cerr << endl;
//}

  int nr = ring_closures_found.number_elements();
  for (int i = 0; i < nr; i++)
  {
    _place_ring_closure(smiles, a, ring_closures_found[i]);
  }

  return 1;
}

static int
locate_item_in_array (int needle, int size_of_haystack, const int * haystack, int istart)
{
  for (int i = istart; i < size_of_haystack; i++)
  {
    if (needle == haystack[i])
      return i;
  }
  
  return -1;
};

//#define DEBUG_APPEND_RING_CLOSURES_FOR_ATOM

/*
   Our implementation of chiral atoms requires ring closures to be added
   before ring openings. BUT, we really don't want the same ring number
   appearing twice on one atom (once as a ring closure, and then as a
   new ring opening). So, we process the
   ring openings first, but put the results into a temporary

   Handle chiral atoms specially
*/

int
Ring_Number_Manager::append_ring_closing_and_opening_digits (IWString & smiles,
                   atom_number_t zatom,
                   const resizable_array<atom_number_t> & ring_closures_found,
                   const resizable_array<const Bond *> & ring_opening_bonds,
                   const Chiral_Centre * c)
{
  assert (ok());

  if (0 == _nr)
    return 1;

#ifdef DEBUG_APPEND_RING_CLOSURES_FOR_ATOM
  cerr << "Adding ring closure digits for atom " << zatom << endl;
  for (int i = 0; i < _nr; i++)
  {
    cerr << "_ring_id[" << i << "] = " << _ring_id[i] << endl;
  }
#endif

  IWString ring_opening_string;

  _generate_ring_opening_chars(ring_opening_string, ring_opening_bonds, zatom);

  if (nullptr != c && ring_closures_found.number_elements() > 1)
    _append_ring_closures_for_chiral_atom(smiles, zatom, ring_closures_found);
  else
  {
    int i = 0;
    while ((i = locate_item_in_array(zatom, _nr + 1, _ring_id, i)) > 0)
    {
      _process_ring(smiles, i, zatom);
    }
  }

  if (ring_opening_string.length())
    smiles += ring_opening_string;

  return 1;
}
