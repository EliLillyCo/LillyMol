#include <stdlib.h>
#include <iostream>
#include <memory>
#include <random>
#include <algorithm>

#include "assert.h"

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwbits/iwbits.h"

// Define this to get the molecule private functions defined here

#define COMPILING_MOLECULER_CC
#define COMPILING_SMILES_CC
#define COMPILING_CTB

#include "molecule.h"
#include "path.h"
#include "smiles.h"
#include "aromatic.h"
#include "misc2.h"
#include "iwrnm.h"
#include "iwrcb.h"

using std::cerr;

// Values for _smiles_order_type

#define INVALID_SMILES_ORDER_TYPE 0
#define UNIQUE_SMILES_ORDER_TYPE 2
#define RANDOM_SMILES_ORDER_TYPE 3
#define DEFAULT_SMILES_ORDER_TYPE 4
#define SUBSET_SMILES_ORDER_TYPE 5
#define USER_SPECIFIED_SMILES_ORDER 6

enum class SmilesOrderType {
  kInvalid,
  kUnique,
  kRandom,
  kDefault,
  kSubset,
  kUserSpecified
};

static unsigned int random_smiles_default_seed = 3172776704;

/*
  We have various classes that guide the atom ordering when building a smiles
  Each class must have two methods, next_starting_atom, and next_atom with the
  signatures below.
  Obviously I could have done this with virtual class heirarchy.
*/

class Atom_Chooser_Default
{
  private:
  public:
    Atom_Chooser_Default() {}

    int next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom) const;
    int next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom) const;
//  int next_atom(const Molecule & m, const atom_number_t zatom, Set_of_Atoms & unprocessed_connections, const int * zorder, atom_number_t & b) const;
};

int
Atom_Chooser_Default::next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom) const
{
  const int matoms = m.natoms();

  const int include_chiral_info = include_chiral_info_in_smiles();

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  for (int i = 0; i < matoms; ++i)
  {
    if (zorder[i] >= 0)   // already in smiles
      continue;

    if (include_atom != nullptr && 0 == include_atom[i])
      continue;

    const Chiral_Centre * c = m.chiral_centre_at_atom(i);

    if (c && include_chiral_info)
    {
      if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found)
        atom_to_return_if_nothing_else_found = i;
    }
    else
    {
      a = i;
      return 1;
    }
  }

  if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found)
    return 0;

  a = atom_to_return_if_nothing_else_found;

  return 1;
}

#ifdef NEW_VERSION_MAYBE
int
Atom_Chooser_Default::next_atom(const Molecule & m,
                                const atom_number_t zatom,
                                Set_of_Atoms & unprocessed_connections,
                                const int * zorder,
                                atom_number_t & b) const
{
  const int n = unprocessed_connections.number_elements();
  if (0 == n)
    return 0;

  if (1 == n)
  {
    b = unprocessed_connections[0];
    unprocessed_connections.resize(0);
    return 1;
  }

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = unprocessed_connections[i];

    if (zorder[j] >= 0)   // we could remove it
      continue;

    if (1 == m.ncon(j))
    {
      b = j;
      return 1;
    }
  }
}
#endif

int
Atom_Chooser_Default::next_atom (const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & next_atom, const int * include_atom) const
{
  next_atom = INVALID_ATOM_NUMBER;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();
  for (int i = 0; i < acon; ++i)
  {
    atom_number_t j;
    bond_type_t   bt;
    a->other_and_type(zatom, i, j, bt);

    if (zorder[j] >= 0)    // atom already classified
      continue;

    if (include_atom != nullptr && 0 == include_atom[j])
      continue;

    if (1 == m.ncon(j))
    {
      next_atom = j;
      return 1;
    }

    if (! IS_SINGLE_BOND(bt))
    {
      next_atom = j;
      return 1;
    }

    if (INVALID_ATOM_NUMBER == next_atom)
      next_atom = j;
  }

  return INVALID_ATOM_NUMBER != next_atom;
}


class Atom_Chooser_Lowest_Rank : public Atom_Chooser_Default
{
  private:
    const int * _rank;

  public:
    Atom_Chooser_Lowest_Rank(const int * r) : _rank(r) {}

    int next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom) const;
    int next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom) const;
};

int
Atom_Chooser_Lowest_Rank::next_starting_atom (const Molecule & m,
                                              const int * zorder,
                                              atom_number_t & a,
                                              const int * include_atom) const
{
  const int matoms = m.natoms();

  int lowest_rank = 0;
  a = INVALID_ATOM_NUMBER;

  const int include_chiral_info = include_chiral_info_in_smiles();

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;
  int rank_last_resort_atom = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (zorder[i] >= 0)   // already in smiles
      continue;

    if (include_atom != nullptr && 0 == include_atom[i])
      continue;

    const Chiral_Centre * c = m.chiral_centre_at_atom(i);

    if (c && include_chiral_info)
    {
      if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found || _rank[i] < rank_last_resort_atom)
      {
        atom_to_return_if_nothing_else_found = i;
        rank_last_resort_atom = _rank[i];
      }
    }
    else if (INVALID_ATOM_NUMBER == a || _rank[i] < lowest_rank)
    {
      lowest_rank = _rank[i];
      a = i;
    }
  }

  if (INVALID_ATOM_NUMBER != a)
    return 1;

  if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found)
    return 0;

  a = atom_to_return_if_nothing_else_found;

  return 1;
}

int
Atom_Chooser_Lowest_Rank::next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom) const
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  int lowest_score = 0;
  atom_number_t atom_with_lowest_score = INVALID_ATOM_NUMBER;

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);
    if (zorder[j] >= 0)
      continue;

    if (include_atom != nullptr && 0 == include_atom[j])
      continue;

    if (INVALID_ATOM_NUMBER == atom_with_lowest_score || _rank[j] < lowest_score)
    {
      atom_with_lowest_score = j;
      lowest_score = _rank[j];
    }
  }

  if (INVALID_ATOM_NUMBER == atom_with_lowest_score)
    return 0;

  b = atom_with_lowest_score;

  return 1;
}

class Atom_Chooser_Random
{
  private:
    std::default_random_engine _rng;

  public:
    Atom_Chooser_Random(const unsigned int s = 0);

    int next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom);
    int next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom);
};

Atom_Chooser_Random::Atom_Chooser_Random(const unsigned int s)
{
  if (0 == s)
  {
    std::random_device rd;
    _rng.seed(rd());
  }
  else
    _rng.seed(random_smiles_default_seed);

  return;
}

int
Atom_Chooser_Random::next_starting_atom (const Molecule & m,
                                         const int * zorder,
                                         atom_number_t & a,
                                         const int * include_atom)
{
  const int matoms = m.natoms();

  std::uniform_int_distribution<int> u(0, matoms - 1);

  int ndx = u(_rng);
//cerr << " ndx " << ndx << '\n';

  const int include_chiral_info = include_chiral_info_in_smiles();

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  for (int i = 0; i < matoms; ++i, ++ndx)
  {
    if (ndx == matoms)
      ndx = 0;

    if (zorder[ndx] >= 0)
      continue;

    if (include_atom != nullptr && 0 == include_atom[ndx])
      continue;

    const Chiral_Centre * c= m.chiral_centre_at_atom(ndx);

    if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found && include_chiral_info && c)
      atom_to_return_if_nothing_else_found = ndx;
    else if (include_chiral_info && c)
      ;
    else
    {
      a = ndx;
      return 1;
    }
  }

  if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found)
    return 0;

  a = atom_to_return_if_nothing_else_found;

  return 1;
}

int
Atom_Chooser_Random::next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & next_atom, const int * include_atom)
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  next_atom = INVALID_ATOM_NUMBER;
  atom_number_t multiple_bond = INVALID_ATOM_NUMBER;

  std::bernoulli_distribution u(0.5);

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t   bt;
    a->other_and_type(zatom, i, j, bt);

    if (zorder[j] >= 0)
      continue;

    if (include_atom != nullptr && 0 == include_atom[j])
      continue;

    if (! IS_SINGLE_BOND(bt))
    {
      if (INVALID_ATOM_NUMBER == multiple_bond)
        multiple_bond = j;
      else if (u(_rng))
        multiple_bond = j;
    }
    else if (INVALID_ATOM_NUMBER == next_atom)
      next_atom = j;
    else if (u(_rng))
      next_atom = j;
  }

  if (INVALID_ATOM_NUMBER != multiple_bond)
  {
    next_atom = multiple_bond;
    return 1;
  }

  return INVALID_ATOM_NUMBER != next_atom;
}

class Atom_Chooser_Unique
{
  private:
    const int * _canonical_rank;

  public:
    Atom_Chooser_Unique(const int * r) : _canonical_rank(r) {}

    int next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom) const;
    int next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom) const;
};

int
Atom_Chooser_Unique::next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & first_atom, const int * include_atom) const
{
  int include_chiral_info = include_chiral_info_in_smiles();

  const int matoms = m.natoms();

  int min_ncon = matoms;    // greater than all ncon() values

  first_atom = INVALID_ATOM_NUMBER;
  int rank_first_atom = 0;

  atom_number_t chiral_atom_if_we_need_to = INVALID_ATOM_NUMBER;   // in case everything chiral [P@]12[P@]3[P@@]1[P@@]23
  int rank_chiral_atom = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (zorder[i] >= 0)    // already processed
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && m.chiral_centre_at_atom(i))
    {
      if (chiral_atom_if_we_need_to < 0 || _canonical_rank[i] > rank_chiral_atom)
      {
        chiral_atom_if_we_need_to = i;
        rank_chiral_atom = _canonical_rank[i];
      }
      continue;
    }

    const Atom * a = m.atomi(i);

    const int ncon = a->ncon();

    if (ncon < min_ncon)
    {
      min_ncon = ncon;
      first_atom = i;
      rank_first_atom = _canonical_rank[i];
    }
    else if (ncon == min_ncon && _canonical_rank[i] > rank_first_atom)
    {
      first_atom = i;
      rank_first_atom = _canonical_rank[i];
    }
  }

  if (INVALID_ATOM_NUMBER != first_atom)
    return 1;

  if (chiral_atom_if_we_need_to >= 0)
  {
    first_atom = chiral_atom_if_we_need_to;
    return 1;
  }

  return 0;
}

//#define DEBUG_UNIQUE_SMILES_ORDERING



int
Atom_Chooser_Unique::next_atom(const Molecule & m, const int * zorder,
                               const atom_number_t current_atom, atom_number_t & b,
                               const int * include_atom) const
{
// Variables for bond order decisions

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int hbc_save = 0;                 // initialised to shut gcc up
  int highest_bond_count = 0;     // definitely initialised to 0

  const Atom * a = m.atomi(current_atom);

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(current_atom);

    if (zorder[j] >= 0)
      continue;

    if (include_atom != nullptr && 0 == include_atom[j])
      continue;

    int bcount;

//  Note that the counting here is different from in the function number_of_bonds().
//  In that function, single and double bonds are counted first, here we want to catch
//  aromatic bonds first as contributing one to the bond count,
//  modulo include_bond_aromaticity_in_smiles

    if (b->is_aromatic()) {
      if (! include_bond_aromaticity_in_smiles())
        bcount = b->number_of_bonds();
      else
        bcount = 1;
    } else {
      bcount = b->number_of_bonds();
    }

#ifdef DEBUG_UNIQUE_SMILES_ORDERING
    cerr << "Choosing next unique atom, " << j << " (bcount = " << bcount << ", rj = " << canonical_rank(j) << ")\n";
#endif

    if (bcount < highest_bond_count)   // definitely not
      continue;

    const int rj = _canonical_rank[j];

    if (bcount > highest_bond_count)
    {
      highest_bond_count = bcount;
      zdefault = j;
      hbc_save = rj;
    }
    else if (rj > hbc_save)    // bcount == highest_bond_count
    {
#ifdef DEBUG_UNIQUE_SMILES_ORDERING
      cerr << "Bonds equal, rj = " << rj << " hbc_save = " << hbc_save << '\n';
#endif

      zdefault = j;
      hbc_save = rj;
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    b = zdefault;

#ifdef DEBUG_UNIQUE_SMILES_ORDERING
    cerr << "Next unique atom is atom " << b << '\n';
#endif

    return 1;
  }

  return 0;
}

class Atom_Chooser_Specific_Atom : public Atom_Chooser_Default
{
  private:
    Set_of_Atoms _start_atom;

// private functions

    int _first_unselected (const int matoms, const int * zorder, atom_number_t & a, const int * include_atom) const;

  public:
    Atom_Chooser_Specific_Atom(const atom_number_t a) { _start_atom.add(a);}
    Atom_Chooser_Specific_Atom(const Atom_Chooser_Specific_Atom & s) : _start_atom(s._start_atom) {}

    int next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom);
    int next_atom(const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom) const;
};

int
Atom_Chooser_Specific_Atom::_first_unselected (const int matoms,
                                               const int * zorder,
                                               atom_number_t & a,
                                               const int * include_atom) const
{
  for (int i = 0; i < matoms; ++i)
  {
    if (zorder[i] >= 0)
      continue;

    if (include_atom != nullptr && 0 == include_atom[i])
      continue;

    a = i;
    return 1;
  }

  return 0;
}

int
Atom_Chooser_Specific_Atom::next_starting_atom(const Molecule & m, const int * zorder, atom_number_t & a, const int * include_atom)
{
  const int n = _start_atom.number_elements();

  if (0 == n)
    return _first_unselected(m.natoms(), zorder, a, include_atom);

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  const int include_chiral_info = include_chiral_info_in_smiles();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = _start_atom[i];
    if (zorder[j] >= 0)     // at this stage we could remove j from the array if we wanted to
      continue;

    if (include_atom != nullptr && 0 == include_atom[j])
    {
      cerr << "Atom_Chooser_Specific_Atom:warning, selected start atom not included " << m.name() << '\n';
      continue;
    }

    if (include_chiral_info && m.chiral_centre_at_atom(j))
      atom_to_return_if_nothing_else_found = j;
    else
    {
      a = j;
      _start_atom.remove_item(i);
      return 1;
    }
  }

  if (INVALID_ATOM_NUMBER != atom_to_return_if_nothing_else_found)
  {
    a = atom_to_return_if_nothing_else_found;
    return 1;
  }
  
// None of the specified atoms available, choose first

  return _first_unselected(m.natoms(), zorder, a, include_atom);
}

int
Atom_Chooser_Specific_Atom::next_atom (const Molecule & m, const int * zorder, const atom_number_t zatom, atom_number_t & b, const int * include_atom) const
{
  return this->Atom_Chooser_Default::next_atom(m, zorder, zatom, b, include_atom);
}

/*
  I need an object to describe how to choose the first atom in a smiles
*/

#define SMILES_FIRST_ATOM_DEFAULT 1
#define SMILES_FIRST_ATOM_UNIQUE 4

/*
  Note that if 
*/

class Smiles_First_Atom
{
  private:
    int _how_to_choose;    // one of the defined values above

    atom_number_t _a;

  public:
    Smiles_First_Atom();

    int smdeflt() const { return SMILES_FIRST_ATOM_DEFAULT == _how_to_choose;}
    int unique () const { return SMILES_FIRST_ATOM_UNIQUE  == _how_to_choose;}

    void set_build_type (int s) { _how_to_choose = s;}

    void set_atom_number (atom_number_t a) { _a = a;}

    int atom_specified (atom_number_t &) const;
    int atom_specified() const { return INVALID_ATOM_NUMBER != _a;}
    void unset_first_atom() { _a = INVALID_ATOM_NUMBER; _how_to_choose = SMILES_FIRST_ATOM_DEFAULT;}
};

Smiles_First_Atom::Smiles_First_Atom()
{
  _how_to_choose = SMILES_FIRST_ATOM_DEFAULT;

  _a = INVALID_ATOM_NUMBER;

  return;
}

int
Smiles_First_Atom::atom_specified (atom_number_t & first_atom) const
{
  if (INVALID_ATOM_NUMBER == _a)
    return 0;

  first_atom = _a;

  return 1;
}

//#define DEBUG_SMILES_CHOOSE_FIRST_ATOM
#ifdef DEBUG_SMILES_CHOOSE_FIRST_ATOM
static ostream &
operator << (ostream & os, const Smiles_First_Atom & smfa)
{
  os << "SMFA: ";
  if (smfa.atom_specified())
  {
    atom_number_t a;
    (void) smfa.atom_specified(a);
    os << "atom " << a;
  }
  else if (smfa.smdeflt())
    os << "default";
  else if (smfa.nice())
    os << "nice";
  else if (smfa.unique())
    os << "unique";
  else if (smfa.random())
    os << "random";
  else 
    os << " HUH";

  return os;
}
#endif

static int _write_aromatic_bonds_as_colons = 0;

void
set_write_smiles_aromatic_bonds_as_colons(int s)
{
  _write_aromatic_bonds_as_colons = s;
}

int
write_smiles_aromatic_bonds_as_colons()
{
  return _write_aromatic_bonds_as_colons;
}

void
Smiles_Information::_default_values()
{
  _smiles_order_type = INVALID_SMILES_ORDER_TYPE;

  _smiles_order = nullptr;

  _smiles_is_smarts = 0;

  _create_smarts_embedding = nullptr;

  _user_specified_atomic_smarts = nullptr;

  return;
}

Smiles_Information::Smiles_Information(int natoms) : _natoms(natoms)
{
  _default_values();

  return;
};

Smiles_Information::Smiles_Information() : _natoms(-1)
{
  _default_values();

  return;
}

Smiles_Information::~Smiles_Information()
{
  if (nullptr != _smiles_order)
    delete [] _smiles_order;

  if (nullptr != _create_smarts_embedding)
    delete [] _create_smarts_embedding;

  if (nullptr != _user_specified_atomic_smarts)
    delete [] _user_specified_atomic_smarts;

  return;
}

int
Smiles_Information::debug_print(std::ostream & os) const
{
  if (nullptr != _smiles_order)
  {
    if (UNIQUE_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Unique smiles order computed\n";
    else if (RANDOM_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Random smiles order computed\n";
    else if (DEFAULT_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Default smiles order computed\n";
    else if (SUBSET_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Subset smiles order computed\n";
    else
      os << "Hmmm, smiles order type is " << _smiles_order_type << '\n';
  }

  _ring_closure_bonds.write_bonds(os);

  if (_smiles_start_atom.number_elements() > 0)
  {
    for (int i = 0; i < _smiles_start_atom.number_elements(); i++)
    {
      os << "smiles in fragment " << i << " starts with atom " << _smiles_start_atom[i] << '\n';
    }
  }

  if (_smiles.length())
    os << "Smiles is '" << _smiles << "'\n";

  return os.good();
}

void
Smiles_Information::make_empty()
{
  _smiles = EMPTY_MOLECULE_SMILES;

  if (_smiles_order != nullptr)
  {
    delete [] _smiles_order;
    _smiles_order = nullptr;
  }

  return;
}

int
Smiles_Information::prepare_to_build_ordering(int matoms)
{
  _smiles_order_type = INVALID_SMILES_ORDER_TYPE;

  if (_smiles_order == nullptr)
    _smiles_order = new_int(matoms, -1);
  else
    set_vector(_smiles_order, matoms, -1);

  _natoms = matoms;

  if (! _ring_closure_bonds.activate(matoms))
    return 0;

  _smiles_start_atom.resize_keep_storage(0);

  return 1;
}


int
Smiles_Information::prepare_to_build_smiles(int matoms)
{
  _smiles.resize_keep_storage(0);

  if (_smiles.elements_allocated() < 3 * matoms)
    _smiles.resize(3 * matoms);

  if (_atom_order_in_smiles.number_elements() > 0)
    _atom_order_in_smiles.resize_keep_storage(0);
  else
    _atom_order_in_smiles.resize(matoms);

  _natoms = matoms;

  return 1;
}

void
Smiles_Information::invalidate()
{
  _natoms = -1;               // added IAW Aug 2014

  _smiles_order_type = INVALID_SMILES_ORDER_TYPE;

  _smiles.resize_keep_storage(0);

  if (_smiles_order != nullptr)
  {
    delete [] _smiles_order;
    _smiles_order = nullptr;
  }

  return;
}

int
Smiles_Information::create_smarts_embedding(atom_number_t zatom) const
{
  if (nullptr == _create_smarts_embedding)  // should be a fatal error
    return 0;

  return _create_smarts_embedding[zatom];
}

int
Smiles_Information::set_create_smarts_embedding(int s)
{
  assert(_natoms > 0);

  if (nullptr == _create_smarts_embedding)
    _create_smarts_embedding = new_int(_natoms, s);
  else
    set_vector(_create_smarts_embedding, _natoms, s);

  return 1;
}

int
Smiles_Information::set_create_smarts_embedding(atom_number_t zatom,
                                                int s)
{
  assert(_natoms > 0 && zatom >= 0 && zatom < _natoms);

  if (nullptr == _create_smarts_embedding)
    _create_smarts_embedding = new_int(_natoms);

  _create_smarts_embedding[zatom] = s;

  return 1;
}

/*
*/

int
Molecule::_smiles_choose_first_atom(const int * zorder,
                                    Smiles_First_Atom & smfa,
                                    atom_number_t & first_atom,
                                    const int * include_atom)
{
#ifdef DEBUG_SMILES_CHOOSE_FIRST_ATOM
  cerr << "INto _smiles_choose_first_atom " << smfa << '\n';
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " order " << zorder[i] << '\n';
  }
#endif

  if (smfa.smdeflt())
    return _smiles_choose_first_atom(zorder, first_atom, include_atom);

  if (smfa.unique())
  {
    if (include_atom == nullptr)
      return _smiles_choose_unique_first_atom(zorder, first_atom);
    else
      return _smiles_choose_unique_first_atom(zorder, first_atom, include_atom);
  }

  cerr << "Molecule::_smiles_choose_first_atom: no method specified\n";
  iwabort();

  return 0;
}

/*
  Note that this is incorrect, a chiral atom can start a smiles, but
  not sure if that works or not...

  Mar 2004. Ran into problems with a subset. Every member of the molecule subset
  was a chiral centre. Therefore save a possible match in ATOM_TO_RETURN_IF_NOTHING_ELSE_FOUND
*/

int
Molecule::_smiles_choose_first_atom(const int * zorder,
                                    atom_number_t & first_atom,
                                    const int * include_atom)
{
  int include_chiral_info = include_chiral_info_in_smiles();

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already done
      continue;

    if (include_atom != nullptr && 0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

#ifdef DEBUG_SMILES_CHOOSE_FIRST_ATOM
    cerr << "Can the smiles start with atom " << i << '\n';
    cerr << "include_chiral_info " << include_chiral_info << '\n';
    cerr << "ncon " << _things[i]->ncon() << '\n';
    cerr << "chiral " << chiral_centre_at_atom(i) << '\n';
    if (include_chiral_info && _things[i]->ncon() > 2 && chiral_centre_at_atom(i))
      cerr << "Nope, that looks chiral\n";
#endif

    if (include_chiral_info && _things[i]->ncon() > 2 && chiral_centre_at_atom(i))
    {
      atom_to_return_if_nothing_else_found = i;
      continue;
    }

    first_atom = i;
    return 1;
  }

  if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found)
    return 0;

  first_atom = atom_to_return_if_nothing_else_found;

  return 1;
}


/*
  External entry point for setting the random number used for
  random smiles generation
*/

void
set_smiles_random_number_seed (random_number_seed_t seed)
{
  random_smiles_default_seed = seed;
}

random_number_seed_t
set_smiles_random_number_seed_random()
{
  std::random_device rd;

  random_smiles_default_seed = rd();

  return random_smiles_default_seed;
}

/*
  As you can see, this is really not truly random, in that the
  first atoms will be favoured.
*/

/*int
Molecule::_smiles_choose_random_first_atom(const int * zorder,
                    atom_number_t & first_atom,
                    const int * include_atom)
{
  int include_chiral_info = include_chiral_info_in_smiles();

  int istart = smiles_random_number_stream.intbtwij(0, _number_elements);

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  for (int i = istart; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already done
      continue;

    if (include_atom != nullptr && 0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
    {
      atom_to_return_if_nothing_else_found = i;
      continue;
    }

    first_atom = i;
    return 1;
  }

// If we come to here, we did not find a suitable atom in the
// range [istart..matoms). How about (0..istart)

  for (int i = 0; i < istart; i++)
  {
//  cerr << "Checking random start atom " << i << " zorder = " << zorder[i] << '\n';
    if (zorder[i] >= 0)    // already processed
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info_in_smiles() && chiral_centre_at_atom(i))
    {
      atom_to_return_if_nothing_else_found = i;
      continue;
    }

    first_atom = i;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != atom_to_return_if_nothing_else_found)
  {
    first_atom = atom_to_return_if_nothing_else_found;
    return 1;
  }

//cerr << "No random start atom available\n";

// Nothing possible

  return 0;
}*/

/*
  For unique smiles, we choose the highest ranked, most lowly connected atom
*/

int
Molecule::_smiles_choose_unique_first_atom(const int * zorder,
                                           atom_number_t & first_atom)
{
  const int * canonical_rank = _symmetry_class_and_canonical_rank.canonical_rank();

  assert(nullptr != canonical_rank);

  int include_chiral_info = include_chiral_info_in_smiles();

  int min_ncon = _number_elements;    // greater than all ncon() values

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int rsave = 0;

  atom_number_t chiral_atom_if_we_need_to = INVALID_ATOM_NUMBER;   // in case everything chiral [P@]12[P@]3[P@@]1[P@@]23

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already processed
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
    {
      if (chiral_atom_if_we_need_to < 0)
        chiral_atom_if_we_need_to = i;
      continue;
    }

    const Atom * a = _things[i];

    int ncon = a->ncon();

    if (ncon < min_ncon)
    {
      min_ncon = ncon;
      zdefault = i;
      rsave = canonical_rank[i];
    }
    else if (ncon == min_ncon && canonical_rank[i] > rsave)
    {
      zdefault = i;
      rsave = canonical_rank[i];
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    first_atom = zdefault;
    return 1;
  }

  if (chiral_atom_if_we_need_to >= 0)
  {
    first_atom = chiral_atom_if_we_need_to;
    return 1;
  }

  return 0;
}

/*
  When dealing with a subset, we need to take into account only the number
  of atoms connected in the subset
*/

int
Molecule::_smiles_choose_unique_first_atom(const int * zorder,
                                           atom_number_t & first_atom,
                                           const int * include_atom)
{
  const int * canonical_rank = _symmetry_class_and_canonical_rank.canonical_rank();

  assert(nullptr != canonical_rank);
  assert(include_atom != nullptr);

  int include_chiral_info = include_chiral_info_in_smiles();

  int min_ncon = _number_elements;    // greater than all ncon() values

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int rsave = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already processed
      continue;

    if (0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
      continue;

    const Atom * a = _things[i];

    int ncon = a->ncon(i, include_atom);

    if (ncon < min_ncon)
    {
      min_ncon = ncon;
      zdefault = i;
      rsave = canonical_rank[i];
    }
    else if (ncon == min_ncon && canonical_rank[i] > rsave)
    {
      zdefault = i;
      rsave = canonical_rank[i];
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    first_atom = zdefault;
    return 1;
  }

  return 0;
}


/*
  The default (fast) behaviour is to find the first singly connected
  and follow that. Otherwise follow the first multiple bond
*/

int
Molecule::_smiles_choose_next_atom(const int * zorder,
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom)
{
  atom_number_t zdefault = INVALID_ATOM_NUMBER;

  const Atom * a = _things[current_atom];

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t   bt;
    a->other_and_type(current_atom, i, j, bt);

    if (zorder[j] >= 0)    // atom already classified
      continue;

    if (include_atom != nullptr && 0 == include_atom[j])
      continue;

    if (1 == _things[j]->ncon())
    {
      next_atom = j;
      return 1;
    }

    if (! IS_SINGLE_BOND(bt))
    {
      next_atom = j;
      return 1;
    }
    else if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = j;
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;
    return 1;
  }

  return 0;
}

/*
  The rule for unique order is connection with the highest bond order
  first, with ties broken by unique ordering.
*/


//#define DEBUG_BUILD_SMILES_ORDERING

/*
  Look for ring closure bonds attached to the atom, then continue
  the search.

  Turns out we can get slightly quicker run times if we can pass N by value rather than
  by reference. But if it is a complex object, then that fails.
*/

#ifdef NEW_VEWRSIONASDLKASD

template <typename N>
int
Molecule::_build_smiles_ordering_fctr(N identify_next_atom,
                                      const atom_number_t previous_atom,
                                      const atom_number_t zatom,
                                      int & icounter,
                                      const int * include_atom,
                                      Smiles_Information & smi_info)
{
  int * zorder = smi_info.smiles_order();

  assert(zatom >= 0 && zatom < _number_elements && zorder[zatom] < 0);

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "_build_smiles_ordering continues with atom " << zatom << '\n';
#endif

  zorder[zatom] = icounter;

  icounter++;

  const Atom * a = _things[zatom];

  const int acon = a->ncon();

  if (INVALID_ATOM_NUMBER != previous_atom && 1 == acon)  // got to a terminal atom
    return 0;

  Set_of_Atoms unprocessed_connections;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);;

    if (previous_atom == j)      // the atom from which we came
      continue;

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "  Atom " << j << " is connected";
  if (zorder[j] >= 0)
    cerr << ". Ring closure detected";
  cerr << '\n';
#endif

    if (zorder[j] >= 0)     // we have found a ring closure
      smi_info.add_ring_closure_bond(zatom, j);
    else if (include_atom != nullptr && 0 == include_atom[j])
      ;
    else 
      unprocessed_connections.add(j);
  }

  if (unprocessed_connections.empty())
    return 0;

  if (1 == unprocessed_connections.number_elements())
    return 1 + _build_smiles_ordering_fctr(identify_next_atom, zatom, unprocessed_connections[0], icounter, include_atom, smi_info);

  atom_number_t b;
  int rc = 0;
  while (identify_next_atom.next_atom(*this, unprocessed_connections, zorder, zatom, b))
  {
    rc += _build_smiles_ordering_fctr(identify_next_atom, zatom, b, icounter, include_atom, smi_info);
  }

  return rc;
}
#endif

template <typename N>
int
Molecule::_build_smiles_ordering_fctr(N identify_next_atom,
                                      atom_number_t previous_atom,
                                      atom_number_t zatom,
                                      int & icounter,
                                      const int * include_atom,
                                      Smiles_Information & smi_info)
{
  int * zorder = smi_info.smiles_order();

  assert(zatom >= 0 && zatom < _number_elements && zorder[zatom] < 0);

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "_build_smiles_ordering continues with atom " << zatom << '\n';
#endif

  while (1)
  {
    zorder[zatom] = icounter;

    icounter++;

    const Atom * a = _things[zatom];

    const int acon = a->ncon();

    if (INVALID_ATOM_NUMBER != previous_atom && 1 == acon)  // got to a terminal atom
      return 0;

    atom_number_t next_unprocessed_connection = INVALID_ATOM_NUMBER;

    for (int i = 0; i < acon; i++)
    {
      const atom_number_t j = a->other(zatom, i);;

      if (previous_atom == j)      // the atom from which we came
        continue;

#ifdef DEBUG_BUILD_SMILES_ORDERING
    cerr << "  Atom " << j << " is connected";
    if (zorder[j] >= 0)
      cerr << ". Ring closure detected";
    cerr << '\n';
#endif

      if (zorder[j] >= 0)     // we have found a ring closure
        smi_info.add_ring_closure_bond(zatom, j);
      else if (include_atom != nullptr && 0 == include_atom[j])
        ;
      else if (INVALID_ATOM_NUMBER == next_unprocessed_connection)
        next_unprocessed_connection = j;
      else
        next_unprocessed_connection = _number_elements;    // special value that cannot be an atom number
    }

    if (INVALID_ATOM_NUMBER == next_unprocessed_connection)
      return 0;

    if (_number_elements == next_unprocessed_connection)   // remember special value used above. More than one path from her
      break;

    previous_atom = zatom;
    zatom = next_unprocessed_connection;
  }

//if (_number_elements != next_unprocessed_connection)     // remember special value used above
//  return 1 + _build_smiles_ordering_fctr(identify_next_atom, zatom, next_unprocessed_connection, icounter, include_atom, smi_info);

  atom_number_t b;
  int rc = 0;
  while (identify_next_atom.next_atom(*this, zorder, zatom, b, include_atom))
  {
    rc += _build_smiles_ordering_fctr(identify_next_atom, zatom, b, icounter, include_atom, smi_info);
  }

  return rc;
}
 
/*
  All smiles building comes through here.
*/

template <typename N>
int
Molecule::_build_smiles_ordering_fctr(N next_atom_selector,
                                      const int * include_atom,
                                      Smiles_Information & smi_info)
{
#ifdef DEBUG_BUILD_SMILES_ORDERING
  if (! _fragment_information.contains_valid_data()) {
    cerr << "No valid fragment data\n";
    cerr << "nfrag " << number_fragments() << '\n';
  } else {
    cerr << "nfrag " << number_fragments() << '\n';
  }
#endif
  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

//cerr << "Molecule::_build_smiles_ordering:there are " << _fragment_information.number_fragments() << " fragments in full molecule with " << _number_elements << " atoms\n";

  smi_info.prepare_to_build_ordering(_number_elements);

  int zcounter = 0;    // allocated in increasing order
  atom_number_t a;     // start atom within each fragment

  int frag = 0;

  while (next_atom_selector.next_starting_atom(*this, smi_info.smiles_order(), a, include_atom))
  {
#ifdef DEBUG_BUILD_SMILES_ORDERING
    cerr << "Starting fragment  with atom " << a << " " << const_smarts_equivalent_for_atom(a) << '\n';
#endif

    smi_info.add_start_atom(a);

    _build_smiles_ordering_fctr(next_atom_selector, INVALID_ATOM_NUMBER, a, zcounter, include_atom, smi_info);

    frag++;
  }

// If we built a subset, we may have fewer or more fragments detected

  if (0 == frag)
  {
    cerr << "Molecule::_build_smiles_ordering:no atoms selected!\n";
    return 0;
  }

  assert(include_atom == nullptr ? frag == _fragment_information.number_fragments() : frag > 0);

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "Smiles order array constructed, " << frag << " fragments\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i; // " order is " << zorder[i];
    if (_fragment_information.number_fragments() > 1)
      cerr << " fragment " << _fragment_information.fragment_membership(i);
    cerr << '\n';
  }

//if (_ring_closure_bonds.number_elements())
//{
//  for (int i = 0; i < _ring_closure_bonds.number_elements(); i++) 
//  {
//    const Bond * b = _ring_closure_bonds[i];
//    cerr << "Ring closure bond " << i << " " << *b << '\n';
//  }
//}
#endif

  return 1;
}


/*
  In a subset, is a given bond still in a ring
*/

#ifdef SEEMS_NOT_BEING_USED
int
Molecule::_ring_bond_in_subset(const int * include_atom,
                               atom_number_t a1,
                               atom_number_t a2)
{
  int nr = nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    if (! ri->contains_bond(a1, a2))
      continue;

    if (ri->all_members_set_in_array(include_atom, 1))
      return 1;
  }

  return 0;
}
#endif


/*
  We are processing the connections to ANCHOR, and need to build
  a sorted array of the adjacent bonds which need to be processed.
  We do an insertion sort, based on the ZORDER value of the atom
  at the other end of the bond
*/

static void
insert_bond (atom_number_t anchor,
             const int * zorder,
             resizable_array<const Bond *> & bonds,
             const Bond * b)
{
  int nb = bonds.number_elements();
  atom_number_t zatom = b->other (anchor);

  for (int i = 0; i < nb; i++)
  {
    atom_number_t ai = bonds[i]->other (anchor);
    if (zorder[zatom] < zorder[ai])
    {
      bonds.insert_before (i, b);
      return;
    }
  }

// If we come out here, either the list is empty, or the new atom at
// the end of the bond is ordered higher than all others

  bonds.add(b);
  return;
}

/*
  We need to put the ring closure bonds into canonical order
*/

static void
sort_ring_closures_found(resizable_array<atom_number_t> & ring_closures_found,
                         const int * zorder)
{
  if (2 == ring_closures_found.number_elements())
  {
    atom_number_t a0 = ring_closures_found[0];
    atom_number_t a1 = ring_closures_found[1];

    if (zorder[a0] > zorder[a1])
      ring_closures_found.swap_elements(0, 1);
  }
  else if (3 == ring_closures_found.number_elements())   // poor-man's sort
  {
    atom_number_t a0 = ring_closures_found[0];
    atom_number_t a1 = ring_closures_found[1];
    atom_number_t a2 = ring_closures_found[2];

    if (zorder[a0] > zorder[a1])
    {
      std::swap(a0, a1);
      ring_closures_found.swap_elements(0, 1);
    }
    if (zorder[a1] > zorder[a2])
    {
      std::swap(a1, a2);
      ring_closures_found.swap_elements(1, 2);
    }
    if (zorder[a0] > zorder[a1])
      ring_closures_found.swap_elements(0, 1);
  }

  return;
}

//#define DEBUG_SMILES_FORMATION

/*
  Ran into problems with deciding which atoms should be in the smiles
  and which should be omitted. For now, if it is in the CT, it will
  be processed. 

  Very strange stuff when dealing with directional bonds.
  Consider

     2                                     *
      \                                    *
       1==3                                *
      /    \                               *
     0      4

  And consider starting with atom 2. The smiles is 2\1(\0)=3\4
  BUT, the 0-1 bond will be an / bond.
  Similarly, when starting with atom 0, the smiles is 0/1(/2)=3\4
  BUT, the 1-2 bond will be a \ bond
  When we detect that, we must react appropriately

  The implementation is deliberately not recursive. Recursive stops processing larger molecules
*/

int
Molecule::_construct_smiles_for_fragment(Smiles_Formation_Info & sfi,
                                         Smiles_Information & smi_info)
{
  int * already_done = sfi.already_done();

  IWString & smiles = smi_info.smiles();

  atom_number_t previous_atom = sfi.previous_atom();
  atom_number_t zatom = sfi.zatom();

  const int * include_atom = sfi.include_atom();

#ifdef DOES_NOT_HELP
  if (0 == _things[zatom]->ncon() && include_atom != nullptr && include_atom[zatom])     // a counterion
  {
    _process_atom_for_smiles(sfi, smiles);
    if (include_atom != nullptr)
      include_atom[zatom] = 1;
    return 1;
  }
#endif

  resizable_array<atom_number_t> prev, current;    // stack of various things to avoid recursion
  IWString open_paren, clse_paren, bond_char;

//cerr << "_construct_smiles_for_fragment previous_atom " << previous_atom << '\n';
  prev.add(previous_atom);
  current.add(zatom);
  open_paren.add(' ');
  clse_paren.add(' ');
  bond_char.add(' ');

  // How cis trans bonds and aromaticity are handled.

  const int inc_ctb = include_cis_trans_in_smiles();    // do the call once for efficiency

  unsigned int inc_arom;    // aromatic bonds or not
  if (! sfi.write_smiles())
    inc_arom = 1;
  else
    inc_arom = include_bond_aromaticity_in_smiles();

  if (write_single_bonds_in_smiles())
    inc_arom |= 2;

  const int * zorder = smi_info.smiles_order();

  resizable_array<const Bond *>  ring_opening_bonds;
  resizable_array<atom_number_t> ring_closures_found;

  resizable_array<const Bond *> process_these_bonds;

  while (current.number_elements() > 0)
  {
#ifdef DEBUG_SMILES_FORMATION
    cerr << "stack";
    for (int i = 0; i < current.number_elements(); ++i)
    {
      if (i > 0)
        cerr << '|';

      cerr << ' ' << current[i] << ' ' << prev[i] << ' ' << open_paren[i] << ' ' << bond_char[i];
    }
    cerr  << '\n';
#endif

//  if (current.number_elements() > 1000)
//    cerr << "stack " << current.number_elements() << '\n';
//  cerr << smiles.length() << '\n';

//  cerr << current.number_elements() << ' ' << prev.number_elements() << ' ' << bond_char.number_elements() << ' ' << open_paren.number_elements() << '\n';
    zatom = current.pop();
    previous_atom = prev.pop();
    const char bsymb  = bond_char.pop();
    const char oparen = open_paren.pop();

    if (' ' != oparen)
      clse_paren += ')';

    sfi.set_prev_and_zatom(previous_atom, zatom);

    if (' ' != oparen)
      smiles += oparen;
    if (' ' != bsymb)
      smiles += bsymb;

    already_done[zatom] = 1;

    smi_info.add_atom(zatom);

    const Atom * a = _things[zatom];

    int acon = a->ncon();

//  If there are no other connections from here, no need to go looking
//  for ring openings or prioritising connections

    if (acon > 1)    // most likely case
      ;
    else if ((1 == acon && INVALID_ATOM_NUMBER != previous_atom) || (0 == acon))
    {
      _process_atom_for_smiles(sfi, smiles);
      if (clse_paren.number_elements() && ' ' != clse_paren.pop())
        smiles += ')';
      continue;
    }

//  Build a sorted list of the connections to be processed, identifying
//  any ring openings or closings.
//  Note that we duplicate a lot of code here to avoid repeatedly testing INCLUDE_ATOM

    ring_opening_bonds.resize_keep_storage(0);
    ring_closures_found.resize_keep_storage(0);
    process_these_bonds.resize_keep_storage(0);

    if (include_atom == nullptr)
    {
      for (int i = 0; i < acon; i++)
      {
        const Bond * b = a->item(i);

        atom_number_t j = b->other(zatom);

        if (previous_atom == j)
          continue;

        if (already_done[j])              // closing a ring
          ring_closures_found.add(j);
        else if (smi_info.contains_ring_closure_bond(zatom, j))
          insert_bond(zatom, zorder, ring_opening_bonds, b);
        else                                      // just a regular connection
          insert_bond(zatom, zorder, process_these_bonds, b);
      }
    }
    else
    {
      for (int i = 0; i < acon; i++)
      {
        const Bond * b = a->item(i);
  
        atom_number_t j = b->other(zatom);
  
        if (previous_atom == j)
          continue;

        if (0 == include_atom[j])
          continue;

        if (already_done[j])              // closing a ring
          ring_closures_found.add(j);
        else if (smi_info.contains_ring_closure_bond(zatom, j))
          insert_bond(zatom, zorder, ring_opening_bonds, b);
        else                                      // just a regular connection
          insert_bond(zatom, zorder, process_these_bonds, b);
      }
    }

#ifdef DEBUG_INSERT_BOND
    int npb = process_these_bonds.number_elements();
    for (int i = 0; i < npb; i++)
    {
      const Bond * b = process_these_bonds[i];
      cerr << " i = " << i << " bond to " << b->other(zatom) << " order " << zorder[b->other(zatom)] << '\n';
    }
#endif

    const Chiral_Centre * c = nullptr;
    if (include_chiral_info_in_smiles())
      c = chiral_centre_at_atom(zatom);     // will be nullptr if atom A is not a chiral centre

    if (nullptr != c && ring_closures_found.number_elements() > 1)
      sort_ring_closures_found(ring_closures_found, zorder);

//  Now that we have determined any ring openings, we can append the smiles symbol.
//  We must wait until ring openings are determined for chiral atoms

    (void) _process_atom_for_smiles(sfi, zorder, ring_opening_bonds, ring_closures_found, c, smiles);

    Ring_Number_Manager & rnm = sfi.rnm();
    assert(rnm.ok());

    rnm.append_ring_closing_and_opening_digits(smiles, zatom, ring_closures_found, ring_opening_bonds, c);

#ifdef DEBUG_SMILES_FORMATION
    cerr << "After atom " << zatom << " smiles is now '" << smiles << "'\n";
    if (ring_closures_found.number_elements())
      cerr << "Found " << ring_closures_found.number_elements() << " ring closures\n";
#endif

//  Handle the case of no further connections
//  If there are no connections, there should not be any ring openings

    acon = process_these_bonds.number_elements();
    if (0 == acon) {
      if (ring_opening_bonds.number_elements())
        cerr << "No connections, but " << ring_opening_bonds.number_elements() << " ring openings\n";
      assert( (include_atom == nullptr) ? (ring_opening_bonds.empty()) : 1);
      if (clse_paren.number_elements() && ' ' != clse_paren.pop())
        smiles += ')';
      continue;
    }

    for (int i = acon - 1; i >= 0; --i) {
      const Bond * b = process_these_bonds[i];

      const atom_number_t j = b->other(zatom);

      if (inc_ctb && b->is_directional())
        _process_directional_bond_for_smiles(bond_char, b, j);
      else
        b->append_bond_type_space_for_nothing(bond_char, j, inc_arom);

      if (i != (acon - 1))
        open_paren.add('(');
      else
        open_paren.add(' ');

      prev.add(zatom);
      current.add(j);
    }
  }

  while (clse_paren.number_elements()) {
    if (' ' != clse_paren.pop())
      smiles += ')';
  }

  return 1;
}

/*
  A smiles for a fragment will begin with atom A
  We need to identify the bond down which the smiles will start
*/

const Bond *
Molecule::_identify_first_smiles_bond(atom_number_t zatom,
                                      const int * zorder)
{
  const Atom * a = _things[zatom];
  int acon = a->ncon();

  int smallest_zorder = -1;
  const Bond * rc = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    atom_number_t j = b->other(zatom);

    if (smallest_zorder < 0)
    {
      smallest_zorder = zorder[j];
      rc = b;
    }
    else if (zorder[j] < smallest_zorder)
    {
      smallest_zorder = zorder[j];
      rc = b;
    }
  }

  return rc;
}

int
Molecule::_construct_smiles(const Fragment_Information & frag_info,
                            Smiles_Information & smi_info,
                            const int * include_atom)
{
  int * already_done = new_int(_number_elements); std::unique_ptr<int[]> free_already_done(already_done);

  smi_info.prepare_to_build_smiles(_number_elements);

  const int * zorder = smi_info.smiles_order();

#ifdef DEBUG_SMILES_FORMATION
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " (" << atomic_symbol(i) << ") order is " << zorder[i];
    if (include_atom != nullptr)
      cerr << " include_atom " << include_atom[i];
    cerr << '\n';
  }
#endif

  const Set_of_Atoms & smiles_start_atom = smi_info.smiles_start_atom();

  int n = smiles_start_atom.number_elements();

  assert(n == frag_info.number_fragments());   // must be a smiles start atom in each fragment

  int rc = 0;

  int need_dot = 0;

#ifdef DEBUG_SMILES_FORMATION
  cerr << "Generating smiles for molecule with " << frag_info.number_fragments() << " fragments\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " atom " << i << " '" << smarts_equivalent_for_atom(i) << " in fragment " << frag_info.fragment_membership(i);
    if (include_atom != nullptr)
      cerr << " include " << include_atom[i];
    cerr << '\n';
  }
#endif

  IWString & smiles = smi_info.smiles();

  for (int i = 0; i < n; i++)
  {
    atom_number_t astart = smiles_start_atom[i];

    int f = frag_info.fragment_membership(astart);

    if (include_atom == nullptr)    // no need to check anything
      ;
    else if (include_atom[astart])    // great, that atom is being processed
      ;
    else             // see if we can find something in fragment F
    {
      astart = _choose_highest_canonical_order_in_fragment(f, zorder, include_atom);
      if (INVALID_ATOM_NUMBER == astart)
        continue;
    }

    if (need_dot)
    	smiles += getUseThisCharAsDotInSmiles();
      //smiles += '.';

    int nr = frag_info.rings_in_fragment(f);

#ifdef DEBUG_SMILES_FORMATION
    cerr << "Fragment " << f << " contains " << frag_info.bonds_in_fragment(f) << " bonds and " << frag_info.atoms_in_fragment(f) << " atoms, nr = " << nr << ", start atom " << astart << '\n';
#endif

    Smiles_Formation_Info sfi(_number_elements, nr);
//  cerr << "Doing smarts? " << smi_info.smiles_is_smarts() << '\n';

    if (nullptr != smi_info.user_specified_atomic_smarts())
      sfi.set_user_specified_atomic_smarts(smi_info.user_specified_atomic_smarts());

    if (smi_info.smiles_is_smarts())
    {
      sfi.set_make_smarts_embedding(smi_info.create_smarts_embedding());
//    cerr << "Adding user specified atomic smarts " << smi_info.user_specified_atomic_smarts() << '\n';
    }

    sfi.set_zatom(astart);
    sfi.set_already_done(already_done);
    sfi.set_include_atom(include_atom);

    if (write_smiles_with_smarts_atoms())
      sfi.set_write_smiles(0);

    rc += _construct_smiles_for_fragment(sfi, smi_info);

    need_dot = 1;
  }

  return rc;
}

/*
  The array _smiles_order can be ordered in a number of ways.
  _smiles_order_type must be updated to correspond with the
  ordering type
*/

const IWString &
Molecule::smiles()
{
  assert(ok());

  if (! _smiles_information.contains_smiles())   // need to recompute
    ;
  else if (_smiles_information.smiles_is_smarts())
  {
    _smiles_information.make_empty();
    _smiles_information.set_smiles_is_smarts(0);
  }
  else
    return _smiles_information.smiles();

  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  (void) number_fragments();

  if (! _smiles_information.contains_valid_ordering())
  {
    Atom_Chooser_Default acd;
    if (! _build_smiles_ordering_fctr(acd, nullptr, _smiles_information))
    {
      cerr << "Molecule::smiles: cannot construct ordering\n";
      _smiles_information.set_error();
      return _smiles_information.smiles();
    }

    _smiles_information.set_smiles_order_type(DEFAULT_SMILES_ORDER_TYPE);
  }

  _smiles_information.smiles().resize(4 * _number_elements);

  _construct_smiles(_fragment_information, _smiles_information, nullptr);

  return _smiles_information.smiles();
}


const IWString &
Molecule::smiles(Smiles_Information & smi_info,
                 const int * include_atom)
{
  if (0 == _number_elements) {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  Atom_Chooser_Default acd;
  if (! _build_smiles_ordering_fctr(acd, include_atom, smi_info))
  {
    cerr << "Molecule::smiles:cannot build subset smiles info\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  Fragment_Information frag_info;

  if (! compute_fragment_information(frag_info, include_atom))
  {
    cerr << "Molecule::smiles:cannot compute fragment info of subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  smi_info.set_smiles_order_type(SUBSET_SMILES_ORDER_TYPE);

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}

// There is quite a bit of code duplication with smiles().
// Maybe sometime rationalise those...
int
Molecule::smiles(Smiles_Information& smi_info) {
  if (_number_elements == 0) {
    smi_info.make_empty();
    return 1;
  }

  (void) number_fragments();

  Atom_Chooser_Default acd;
  if (! _build_smiles_ordering_fctr(acd, nullptr, smi_info)) {
    cerr << "Molecule::smiles: cannot construct ordering\n";
    smi_info.set_error();
    return 0;
  }

  smi_info.set_smiles_order_type(DEFAULT_SMILES_ORDER_TYPE);

  smi_info.smiles().resize(4 * _number_elements);

  return _construct_smiles(_fragment_information, smi_info, nullptr);
}

const IWString &
Molecule::random_smiles()
{
  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  invalidate_smiles();

  Atom_Chooser_Random acr(set_smiles_random_number_seed_random());

  (void) _build_smiles_ordering_fctr(acr, nullptr, _smiles_information);

  _smiles_information.set_smiles_order_type(RANDOM_SMILES_ORDER_TYPE);

  _construct_smiles(_fragment_information, _smiles_information, nullptr);

  return _smiles_information.smiles();
}


/*
  For helping people with text based programmes, we have the ability to
  start a smiles with any atom. This is a kludge, because this doesn't
  really fit well with how the smiles are built.
*/

const IWString &
Molecule::smiles_starting_with_atom(atom_number_t astart)
{
  return smiles_starting_with_atom(astart, _smiles_information, nullptr);
}

const IWString &
Molecule::smiles_starting_with_atom(atom_number_t astart,
                                    Smiles_Information & smi_info,
                                    const int * include_atom)
{
  if (0 == _number_elements)
  {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  assert(include_atom == nullptr ? 1 : 0 != include_atom[astart]);

  smi_info.invalidate();

  (void) number_fragments();

  Atom_Chooser_Specific_Atom acsa(astart);

  (void) _build_smiles_ordering_fctr(acsa, include_atom, smi_info);

#ifdef NOT_SURE_WHY_I_HAD_DONE_THIS
  Fragment_Information frag_info;
  if (! compute_fragment_information(frag_info, include_atom))
  {
    cerr << "Molecule::smiles_starting_with_atom:cannot find fragment info for subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }
#endif

  smi_info.set_smiles_order_type(RANDOM_SMILES_ORDER_TYPE);

  _construct_smiles(_fragment_information, smi_info, include_atom);

  return smi_info.smiles();
}

const IWString&
Molecule::aromatic_smiles() {
  auto asave = get_include_aromaticity_in_smiles();
  set_include_aromaticity_in_smiles(1);
  invalidate_smiles();
  smiles();
  set_include_aromaticity_in_smiles(asave);
  return smiles();
}

/*
  Someone may need to know the order of the atoms in the smiles
*/

int
Molecule::smiles_atom_order(int * s)
{
  (void) smiles();     // will force construction of the array(s)

  const int * so = _smiles_information.smiles_order();

  copy_vector(s, so, _number_elements);

  return 1;
}

//#define DEBUG_UNIQUE_SMILES

const IWString &
Molecule::_unique_smiles(const Fragment_Information & frag_info,
                         Smiles_Information & smi_info,
                         Symmetry_Class_and_Canonical_Rank & sccr,
                         const int * include_atom)
{
//cerr << "Allocated? " << sccr.arrays_allocated() << '\n';

  if (0 == _number_elements) {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  smi_info.prepare_to_build_ordering(_number_elements);

  compute_aromaticity_if_needed();

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Computing unique smiles for " << _number_elements << " atoms\n";
  if (! sccr.arrays_allocated())
    cerr << "Computing canonical rank\n";
  else
    cerr << "Canonical rank already computed\n";
#endif

  if (! sccr.arrays_allocated()) {
    compute_canonical_ranking(sccr, include_atom);
  }

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Canonical rank computed\n";
  const int * c = sccr.canonical_rank();
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " type " << _things[i]->atomic_symbol();
    if (include_atom != nullptr)
      cerr << " include " << include_atom[i];
    cerr << " rank " << c[i] << '\n';
  }

  cerr << "Order type " << smi_info.smiles_order_type() << '\n';
#endif

  if (include_cis_trans_in_smiles()) {
    _adjust_cis_trans_bonds_to_canonical_form(sccr.canonical_rank());
  }

  assert(nullptr != _aromaticity);    // aromaticity computed in compute_canonical_ranking

  Smiles_First_Atom smfa;

  smfa.set_build_type(SMILES_FIRST_ATOM_UNIQUE);

  Atom_Chooser_Unique acn(canonical_ranks());
  (void) _build_smiles_ordering_fctr(acn, include_atom, smi_info);

  if (include_atom == nullptr)
    smi_info.set_smiles_order_type(UNIQUE_SMILES_ORDER_TYPE);
  else
    smi_info.set_smiles_order_type(SUBSET_SMILES_ORDER_TYPE);

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Smiles unique order is\n";
  const int * s = smi_info.smiles_order();

  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " i = " << i << " order = " << s[i] << '\n';
  }
#endif

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}

#ifdef VERSION_USING_FUNCTION_POINTER
const IWString &
Molecule::_unique_smiles(const Fragment_Information & frag_info,
                         Smiles_Information & smi_info,
                         Symmetry_Class_and_Canonical_Rank & sccr,
                         const int * include_atom)
{
//cerr << "Allocated? " << sccr.arrays_allocated() << '\n';

  if (0 == _number_elements)
  {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  smi_info.prepare_to_build_ordering(_number_elements);

  compute_aromaticity_if_needed();

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Computing unique smiles for " << _number_elements << " atoms\n";
  if (! sccr.arrays_allocated())
    cerr << "Computing canonical rank\n";
  else
    cerr << "Canonical rank already computed\n";
#endif

  if (! sccr.arrays_allocated()) {
    compute_canonical_ranking(sccr, include_atom);
  }

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Canonical rank computed\n";
  const int * c = sccr.canonical_rank();
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " type " << _things[i]->atomic_symbol();
    if (include_atom != nullptr)
      cerr << " include " << include_atom[i];
    cerr << " rank " << c[i] << '\n';
  }

  cerr << "Order type " << smi_info.smiles_order_type() << '\n';
#endif

  if (include_cis_trans_in_smiles())
    _adjust_cis_trans_bonds_to_canonical_form(sccr.canonical_rank());

  assert(nullptr != _aromaticity);    // aromaticity computed in compute_canonical_ranking

  Smiles_First_Atom smfa;

  smfa.set_build_type(SMILES_FIRST_ATOM_UNIQUE);

  (void) _build_smiles_ordering(smfa,
                                 &Molecule::_smiles_choose_unique_next_atom,
                                 include_atom,
                                 smi_info);

  if (include_atom != nullptr)
    smi_info.set_smiles_order_type(UNIQUE_SMILES_ORDER_TYPE);
  else
    smi_info.set_smiles_order_type(SUBSET_SMILES_ORDER_TYPE);

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Smiles unique order is\n";
  const int * s = smi_info.smiles_order();

  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " i = " << i << " order = " << s[i] << '\n';
  }
#endif

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}
#endif

/*
  Be careful using unique_smiles() and non_aromatic_unique_smiles()
  If you call unique_smiles() and then non_aromatic_unique_smiles()
  you will likely get the same smiles back. To force a recomputation,
  call invalidate_smiles() in between...

  Jul 2000. Introduce a default aromaticity that is used for all
  unique smiles
*/

static int default_unique_smiles_aromaticity = Daylight;

int
set_default_unique_smiles_aromaticity(int a)
{
  default_unique_smiles_aromaticity = a;

  return 1;
}

/*
  We have some common tasks that need to happen when doing unique smiles determinations.
  Some things need to be set and then restored.
*/

class Hold_and_Restore_Global_Settings
{
  private:
    int _aromsave;     // save the global aromaticity definition

    int _incaromsave;  // save the include aromaticity in smiles definition

    int _include_bond_aromaticity_save;

    int _inc_atom_map_save;

  public:
    Hold_and_Restore_Global_Settings (int aromatic_atom, int aromatic_bond);
    ~Hold_and_Restore_Global_Settings();

    int aromaticity_changed() const { return _aromsave != default_unique_smiles_aromaticity;}
};

Hold_and_Restore_Global_Settings::Hold_and_Restore_Global_Settings(int aromatic_atom, int aromatic_bond)
{
  _aromsave = global_aromaticity_type();

  set_global_aromaticity_type(default_unique_smiles_aromaticity);

  _incaromsave = get_include_aromaticity_in_smiles();

  _include_bond_aromaticity_save = include_bond_aromaticity_in_smiles();

  set_include_aromaticity_in_smiles(aromatic_atom, aromatic_bond);

  _inc_atom_map_save = include_atom_map_with_smiles();

  set_include_atom_map_with_smiles(0);

  return;
}

Hold_and_Restore_Global_Settings::~Hold_and_Restore_Global_Settings()
{
  set_global_aromaticity_type(_aromsave);

  set_include_aromaticity_in_smiles(_incaromsave, _include_bond_aromaticity_save);

  set_include_atom_map_with_smiles(_inc_atom_map_save);

  return;
}

const IWString &
Molecule::unique_smiles()
{
  if (UNIQUE_SMILES_ORDER_TYPE == _smiles_information.smiles_order_type())
    return _smiles_information.smiles();

  Hold_and_Restore_Global_Settings hrgs(1, 1);

  if (hrgs.aromaticity_changed())    // we may have Pearlman aromaticity, but need Daylight for unique smiles
    compute_aromaticity();
  else 
    compute_aromaticity_if_needed();

  _smiles_information.set_smiles_is_smarts(0);

  return _unique_smiles(_fragment_information, _smiles_information, _symmetry_class_and_canonical_rank, nullptr);
}

const IWString &
Molecule::unique_smiles(Smiles_Information & smi_info,
                        const int * include_atom)
{
  assert(include_atom != nullptr);

  if (0 == _number_elements) {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  _smiles_information.set_smiles_is_smarts(0);

//_fragment_information.debug_print(cerr);

  Hold_and_Restore_Global_Settings hrgs(1, 1);

  if (hrgs.aromaticity_changed()) {
    compute_aromaticity();
  }

  Fragment_Information frag_info;

  if (! compute_fragment_information(frag_info, include_atom)) {
    cerr << "Molecule::unique_smiles:cannot compute fragment info for subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  Symmetry_Class_and_Canonical_Rank sccr;

  if (! sccr.allocate_arrays(_number_elements)){
    return smi_info.set_error();
  }

  compute_canonical_ranking(sccr, include_atom);

#ifdef DEBUG_SMILES_FORMATION
  cerr << "Canonical ranking computed (external Smiles_Information)\n";
  for (int i = 0; i < _number_elements; ++i) {
    cerr << " atom " << i << " " << _things[i]->atomic_symbol() << " rank " << sccr.canonical_rank()[i] << '\n';
  }
  cerr << "Order type " << smi_info.smiles_order_type() << '\n';
#endif

  if (include_cis_trans_in_smiles()) {
    _adjust_cis_trans_bonds_to_canonical_form(sccr.canonical_rank());
  }

// _smiles_choose_unique_*_atom need to have the molecule's canonical order fixed. Make
// a copy of any existing data in _symmetry_class_and_canonical_rank and store the
// values from sccr into _symmetry_class_and_canonical_rank

  Symmetry_Class_and_Canonical_Rank sccr_save;
  sccr_save.store_values_from(_symmetry_class_and_canonical_rank, _number_elements);

  _symmetry_class_and_canonical_rank.store_values_from(sccr, _number_elements);

  _unique_smiles(frag_info, smi_info, sccr, include_atom);

  _symmetry_class_and_canonical_rank.store_values_from(sccr_save, _number_elements);

  return smi_info.smiles();
}

const IWString &
Molecule::non_aromatic_unique_smiles()
{
  if (0 == _number_elements) {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  _smiles_information.set_smiles_is_smarts(0);

  if (UNIQUE_SMILES_ORDER_TYPE != _smiles_information.smiles_order_type()) {
    invalidate_smiles();
  }

  Hold_and_Restore_Global_Settings hrgs(0, 0);   // kekule atoms, kekule bonds

  if (hrgs.aromaticity_changed())
    compute_aromaticity();

  return _unique_smiles(_fragment_information, _smiles_information, _symmetry_class_and_canonical_rank, nullptr);
}

const IWString&
Molecule::UniqueKekuleSmiles() {
  if (! contains_aromatic_atoms()) {
    return unique_smiles();
  }

  IWString usmi = unique_smiles();  // Our own copy.
  if (! build_from_smiles(usmi)) {
    cerr << "Molecule::UniqueKekuleSmiles:cannot interpret " << usmi << '\n';
    return _smiles_information.set_error();
  }

  Hold_and_Restore_Global_Settings hrgs(1, 0);  // aromatic atoms, kekule bonds.
  _smiles_information.set_smiles_is_smarts(0);
  compute_aromaticity();  // Always recompute.

  return _unique_smiles(_fragment_information, _smiles_information, _symmetry_class_and_canonical_rank, nullptr);
}

/*
  In processing a fused system, it has been discerned that some rings
  which had previously been assigned separate fused system identifiers
  are in fact part of the same fused system.
  We examine rings in RINGS, starting with ring RSTART.
  The common fused system identifier is FUSED_SYSTEM_IDENTIFIER.
  The fused system identifiers which need to be changed are in FUSED_SYS_IDS_TO_BE_CHANGES
*/

int
Molecule::_merge_fused_system_identifiers(resizable_array<Ring *> & rings,
                                          int rstart,
                                          int fused_system_identifier,
                                          resizable_array<int> & fused_sys_ids_to_be_changed)
{
  int rc = 0;
  int nr = rings.number_elements();
  for (int i = rstart; i < nr; i++)
  {
    Ring * r = rings[i];
    if (! r->is_fused())
      continue;

    int rfsysid = r->fused_system_identifier();
    if (rfsysid == fused_system_identifier)
      continue;

    if (fused_sys_ids_to_be_changed.contains(rfsysid))
    {
      r->set_fused_system_identifier(fused_system_identifier);
      r->set_is_fused(1);
      rc++;
    }
  }

  return rc;
}

int
Molecule::_find_raw_rings(const atom_number_t previous_atom,
                          const atom_number_t current_atom,
                          resizable_array<Ring *> & rings,
                          resizable_array<atom_number_t> & active_rings,
                          int * already_done)
{
  assert(0 == already_done[current_atom]);
  assert(rings.number_elements() == active_rings.number_elements());

  already_done[current_atom] = 1;
//cerr << "GBORETN: processing atom " << current_atom << ", ncon " << _things[current_atom]->ncon() << '\n';

  const Atom * c = _things[current_atom];

  int acon = c->ncon();
  if (0 == acon)
    return rings.number_elements();

  Set_of_Atoms to_process;
  to_process.resize(acon);

  for (int i = 0; i < acon; i++)
  {
    const atom_number_t j = c->other(current_atom, i);
    if (previous_atom == j)
      continue;

    if (already_done[j])
    {
//    cerr << "From atom " << current_atom << " found new ring to atom " << j << '\n';
      Ring * tmp = new Ring;
      tmp->resize(6);
      tmp->add(j);
      tmp->add(current_atom);
      rings.add(tmp);
      tmp->set_fragment_membership(_fragment_information.fragment_membership(current_atom));

      active_rings.add(j);
    }
    else
      to_process.add(j);
  }

// Recursively call this function for each bond attached to CURRENT_ATOM

  acon = to_process.number_elements();

  for (int i = 0; i < acon; i++)
  {
    const atom_number_t j = to_process[i]; // c->other(current_atom, i);
    if (already_done[j])
      continue;

    const int rstart = rings.number_elements();
    (void) _find_raw_rings(current_atom, j, rings, active_rings, already_done);

    const int nrings_now = rings.number_elements();
//  cerr << "Looking from " << current_atom << " to " << j << " now have " << nrings_now << " rings, compare " << rstart << '\n';
    if (rstart == nrings_now)   // no new rings found down this bond
      continue;

//  Count the number of new rings down this bond, and try to determine
//  any existing fusion specifications. Note that we may find different
//  fused sys identifiers in the list of new rings. These numbers must
//  all be consolidated, as they now are known to belong to the same
//  fused system,

    int number_new_rings = 0;
    int fused_system_identifier = -99;

//  cerr << "Current = " << current_atom << " to " << j << 
//          " now " << rings.number_elements() << " rings\n";

    resizable_array<int> fused_sys_ids_to_be_changed;

    for (int k = rstart; k < nrings_now; k++)
    {
      if (INVALID_ATOM_NUMBER == active_rings[k])    // ring has been processed
        continue;

      number_new_rings++;
      Ring * r = rings[k];

      if (! r->is_fused())    // not interested in isolated rings
        continue;

      int rsysid = r->fused_system_identifier();    

      if (fused_system_identifier < 0)        // first fused ring found here
        fused_system_identifier = rsysid;
      else if (rsysid == fused_system_identifier)    // already has same id, no change
        ;
      else     // different systems need to be merged
        fused_sys_ids_to_be_changed.add(rsysid);

//    if (r->is_fused())
//      cerr << "FBLOGD current = " << current_atom << " con = " << i << " atom " << j << 
//              " k = " << k << " found fused sys identifier " << r->fused_system_identifier() << 
//              " to " << active_rings[k] << '\n';
    }

    if (0 == number_new_rings)
      continue;

    if (-99 == fused_system_identifier)
      fused_system_identifier = current_atom;
    else
      _merge_fused_system_identifiers(rings, rstart, fused_system_identifier, fused_sys_ids_to_be_changed);

//  cerr << "Will assign fused system identifier " << fused_system_identifier << '\n';

    for (int k = rstart; k < nrings_now; k++)
    {
      if (INVALID_ATOM_NUMBER == active_rings[k])
        continue;

      Ring * r = rings[k];
      if (number_new_rings > 1)
      {
        r->set_fused_system_identifier(fused_system_identifier);
        r->set_is_fused(1);
      }

      if (active_rings[k] == current_atom)    // ring is complete, terminate it
        active_rings[k] = INVALID_ATOM_NUMBER;
      else
        r->add(current_atom);
    }
  }

  assert(rings.number_elements() == active_rings.number_elements());

  return rings.number_elements();
}

//#define DEBUG_FIND_RAW_RINGS_FOR_FRAGMENT

/*
  Process the raw rings for atoms with _fragment_membership == id

  We need to be careful with systems with spiro fusions to ring systems,

  for example. Since the spiro rings are not fused, we can correclty set
  their final ring membership.
  By default, we set the ring membership of fused systems to IW_RING_MEMBERSHIP_IS_A_RING_ATOM,
  but that would overwrite the values correctly found for the spiro rings. So, if we have
  the case of spiro rings joined to a ring system, force a complete SSSR
*/

int
Molecule::_find_raw_rings_for_fragment(int id, int * already_done)
{
  if (0 == nrings())
    return 1;

  if (id < 0 || id >= _fragment_information.number_fragments())
  {
    cerr << "Molecule::_find_raw_rings_for_fragment:finding rings in fragment " << id << " but only " << _fragment_information.number_fragments() << " fragments\n";
    debug_print(cerr);
    assert(NULL == "This is very bad");
  }

  if (nullptr == _ring_membership)
    _initialise_ring_membership();

// Initialise all these atoms as 0 ring membership

  int atoms_being_processed = 0;

  atom_number_t start_atom = INVALID_ATOM_NUMBER;   // first atom in the fragment

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    if (id == fragment_membership[i])
    {
      _ring_membership[i] = 0;
      atoms_being_processed++;
      if (INVALID_ATOM_NUMBER == start_atom)
        start_atom = i;
    }
  }

  int nr = _fragment_information.rings_in_fragment(id);

  assert(nr >= 0);

  if (0 == nr)      // no rings in this fragment
    return 1;

  assert(atoms_being_processed > 2);

  resizable_array<Ring *> rings;

  resizable_array<atom_number_t> active_rings;
  active_rings.resize(nrings());

  _find_raw_rings(INVALID_ATOM_NUMBER, start_atom, rings, active_rings, already_done);

  nr = rings.number_elements();

  assert(nr > 0);

#ifdef DEBUG_FIND_RAW_RINGS_FOR_FRAGMENT
  cerr << "Found " << nr << " rings for fragment " << id << '\n';
  for (int i = 0; i < nr; i++)
  {
    cerr << "Ring " << i << " ";
    const Ring * ri = rings[i];
    cerr << (*ri) << '\n';
    if (! ok_ring(ri))
    {
      cerr << "Very bad news, not a valid ring\n";
      iwabort();
    }
  }
#endif

// Update ring membership with the details

  int number_fused_rings = 0;

  for (int i = 0; i < nr; i++)
  {
    Ring * ri = rings[i];

    if (ri->is_fused())
      number_fused_rings++;
    else
    {
      ri->increment_vector(_ring_membership, 1);
      _add_ring_to_sssr(ri);
    }
  }

  if (0 == number_fused_rings)
    return _assign_fsid_values_to_isolated_rings();

  assert(nr > 1);     // can't be just one fused ring


// Accumulate the fused_system_identifiers of the fused systems that need to be processed with their attached spiro rings

  resizable_array<int> spiro_between_isolated_and_fused;

  if (number_fused_rings < nr)    // must be both fused and isolated rings present - check for spiro fusions between them
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = rings[i];

      if (! ri->is_fused())
        continue;

      if (ri->fused_ring_check_for_spiro_fusion(_ring_membership))
        spiro_between_isolated_and_fused.add_if_not_already_present(ri->fused_system_identifier());
    }
  }

#ifdef DEBUG_FIND_RAW_RINGS_FOR_FRAGMENT
  cerr << "After examining rings, spiro between isolated and fused = " << spiro_between_isolated_and_fused.number_elements() << ", nr = " << nr << '\n';
#endif

// Fused rings that are not bonded to a spiro ring to go raw rings

  for (int i = 0; i < nr; i++)
  {
    Ring * ri = rings[i];

    if (! ri->is_fused())
      continue;

    int fsid = ri->fused_system_identifier();
    if (spiro_between_isolated_and_fused.contains(fsid))
      continue;

    ri->set_vector(_ring_membership, kRingMembershipIsRingAtom);

    _raw_rings.add(ri);
//  _experimental_raw_rings.add(ri);
  }

  if (spiro_between_isolated_and_fused.number_elements())
  {
//  cerr << "Spiro fusion to fused system\n";
    for (int i = 0; i < spiro_between_isolated_and_fused.number_elements(); i++)
    {
      int fsid = spiro_between_isolated_and_fused[i];

      _handle_spiro_between_isolated_and_fused(rings, fsid, already_done);
    }

    for (int i = 0; i < nr; i++)
    {
      Ring * ri = rings[i];

      if (! ri->is_fused())
        continue;

      int fsid = ri->fused_system_identifier();
      if (spiro_between_isolated_and_fused.contains(fsid))
        delete ri;
    }
  }

  return nr;
}

// Return the first non zero index in `values` and set that
// index to 1.
int
Unused(extending_resizable_array<int>& values) {
  for (int i = 0; i < values.number_elements(); ++i) {
    if (values[i] == 0) {
      values[i] = 1;
      return i;
    }
  }

  // Nothing available, append.

  int result = values.size();
  values[result] = 1;
  return result;
}

// Mar 2015. Assign non-negative fsid values to isolated rings
// Oct 2022. Make sure this works for multiple fragments. 

int
Molecule::_assign_fsid_values_to_isolated_rings()
{
  // cerr << "_assign_fsid_values_to_isolated_rings: have " << _sssr_rings.size() << " rings\n";

  for (Ring* ri : _sssr_rings) {
    if (ri->fused_system_identifier() >= 0) {
      continue;
    }
    const int fsid = _unused_fused_system_identifier();
    ri->set_fused_system_identifier(fsid);
    ri->set_is_fused(0);
  }

  return _sssr_rings.size();
}

/*
  A new isolated ring has been found.  Update the _ring_membership,
  but only values greater than 0 - the others are undetermined yet.

  Sept 97. Previously this set ring membership to 1 only if the
  existing value of _ring_membership was zero. This failed in
  the case of spiro fused rings, so now we check for >= 0 values
  and always increment.

  also tell the bonds about the new ring
*/

#ifdef IS_THIS_BEING_CALLED_NO
int
Molecule::_update_ring_membership(const Ring * r)
{
  int ring_size = r->number_elements();
  atom_number_t prev_atom = r->last_item();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r->item(i);

    cerr << "Molecule::_update_ring_membership: atom " << j << " rm = " << _ring_membership[j] << '\n';

    if (kRingMembershipIsRingAtom  == _ring_membership[j])   // probably part of a fused system - spiro fusion
      _ring_membership[j] = 1;
    else if (_ring_membership[j] >= 0)
      _ring_membership[j]++;

    cerr << "Molecule::_update_ring_membership: atom " << j << " rm = " << _ring_membership[j] << '\n';
    const Atom * aj = _things[j];

    Bond * b = const_cast<Bond *>(aj->bond_to_atom(prev_atom));   // loss of const OK
    assert(b);
    if (b == nullptr)
    {
    	return 0;
    }

    b->set_nrings(1);

    prev_atom = j;
  }

  return 1;
}
#endif

int
Molecule::_find_raw_rings(int * already_done)
{
  int nf = _fragment_information.number_fragments();

//cerr << "Finding raw rings for " << nf << " fragments\n";

  for (int i = 0; i < nf; i++)
  {
//  cerr << "Finding raw rings for fragment " << i << '\n';
    if (! _find_raw_rings_for_fragment(i, already_done))
    {
      cerr << "Molecule::_find_raw_rings: Bad news, cannot get raw rings for fragment " << i << '\n';
      debug_print(cerr);
      iwabort();
    }
  }

  return _raw_rings.number_elements();
}

int
Molecule::_find_raw_rings()
{
  assert(_raw_rings.empty());

  (void) number_fragments();

  if (0 == nrings())
    return 1;

  int * tmp = new_int(_number_elements); std::unique_ptr<int[]> free_tmp(tmp);

  return _find_raw_rings(tmp);
}

int
smiles_error_message (const char * smiles,
                      int length_of_smiles,
                      int characters_processed,
                      const char * message)
{
  if (! display_smiles_interpretation_error_messages())
    return 1;

  assert(message);

  cerr << message << '\n';

  int smiles_chars_to_print = characters_processed + 10;
  if (smiles_chars_to_print > length_of_smiles || smiles_chars_to_print > 80)
    smiles_chars_to_print = length_of_smiles;

  for (int i = 0; i < smiles_chars_to_print; i++)
  {
    cerr << smiles[i];
  }
  cerr << '\n';

//cerr << "                     ";
  for (int i = 0; i < characters_processed; i++)
    cerr << ' ';
  cerr << "^\n";

  return 1;
}

void
Molecule::_add_ring_to_sssr(Ring * ri)
{
  ri->set_ring_number(_sssr_rings.number_elements());
  _sssr_rings.add(ri);

// _experimental_sssr_rings.add(ri);

// Update the bond ring membership

  atom_number_t aprev = ri->last_item();

  int ring_size = ri->number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = ri->item(i);

//  _things[j]->in_another_ring();

    Bond * b = const_cast<Bond *>(_things[j]->bond_to_atom(aprev));
   
    b->in_another_ring();

    aprev = j;
  }

  return;
}

//#define DEBUG_HANDLE_SPIRO_BETWEEN_ISOLATED_AND_FUSED

/*
  Spiro fusions between isolated and fused systems are problematic. We force SSSR
  determination of all the fused rings in the fragment
*/

int
Molecule::_handle_spiro_between_isolated_and_fused(const resizable_array<Ring *> & rings,
                                                    int fsid,
                                                    int * atmp)
{
  set_vector(atmp, _number_elements, 0);

  int nr = rings.number_elements();

#ifdef DEBUG_HANDLE_SPIRO_BETWEEN_ISOLATED_AND_FUSED
  cerr << "Handling " << nr << " rings for possible spiro/fused systems, _sssr_rings.number_elements  = " << _sssr_rings.number_elements() << '\n';
#endif

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = rings[i];

    if (! ri->is_fused())
      continue;

    if (fsid != ri->fused_system_identifier())
      continue;

    ri->set_vector(atmp, 1);
  }

  return _pearlman_sssr(atmp, 1);
}

#ifdef NOT_USED
int
Molecule::_all_atoms_are_chain_atoms(const int * process_these_atoms)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (is_ring_atom(i))
      return 0;
  }

  return 1;
}
#endif

/*int
Molecule::_determine_ring_closure_bonds(const int * zorder,
                                        const int * include_atom)
{
  assert(include_atom != nullptr);

  int * already_done = new_int (_number_elements); std::unique_ptr<int[]> free_already_done (already_done);

  int n = _smiles_start_atom.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t astart = _smiles_start_atom[i];

    if (include_atom == nullptr)
      ;
    else if (include_atom[astart])
      ;
    else
    {
      astart = _choose_highest_canonical_order_in_fragment (i, zorder, include_atom);
      if (INVALID_ATOM_NUMBER == astart)
        continue;
    }

    _determine_ring_closure_bonds (INVALID_ATOM_NUMBER, astart, zorder, include_atom, already_done);
  }

  return 1;
}*/

/*
  We need to re-determine the ring closure bonds when dealing with a subset
  The canonical order is already known
*/

/*int
Molecule::_determine_ring_closure_bonds(atom_number_t aprev,
                                        atom_number_t zatom,
                                        const int * zorder,
                                        const int * include_atom,
                                        int * already_done)
{
  assert(include_atom != nullptr);    // subset only

  already_done[zatom] = 1;

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  resizable_array<const Bond *> process_these_bonds;

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item (i);

    atom_number_t j = b->other (zatom);

    if (aprev == j)
      continue;

    if (0 == include_atom[j])
      continue;

    if (already_done[j])
    {
      _ring_closure_bonds.add (zatom, j);
      rc++;
    }
    else
    {
      insert_bond (zatom, zorder, process_these_bonds, b);
    }
  }

  int n = process_these_bonds.number_elements();

  for (int i = 0; i < n; i++)
  {
    const Bond * b = process_these_bonds[i];

    atom_number_t j = b->other (zatom);

    if (already_done[j])
      continue;

    rc += _determine_ring_closure_bonds (zatom, j, zorder, include_atom, already_done);
  }

  return rc;
}*/

/*
  We are doing the smiles for a fragment, but there is a subset, and the molecule's default
  _smiles_start_atom was not part of the subset. Find another suitable starting point
*/

atom_number_t
Molecule::_choose_highest_canonical_order_in_fragment(int f,
                                                      const int * zorder,
                                                      const int * include_atom) const
{
  atom_number_t rc = INVALID_ATOM_NUMBER;
  int z;

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == include_atom[i])
      continue;

    if (f != fragment_membership[i])
      continue;

    if (INVALID_ATOM_NUMBER == rc || z > zorder[i])
    {
      rc = i;
      z = zorder[i];
    }
  }

  return rc;
}

int
Smiles_Information::allocate_user_specified_atomic_smarts()
{
  assert(nullptr == _user_specified_atomic_smarts);

  assert(_natoms > 0);

  _user_specified_atomic_smarts = new IWString[_natoms];

  assert(nullptr != _user_specified_atomic_smarts);

  return 1;
}

IWString &
Smiles_Information::user_specified_atomic_smarts(atom_number_t zatom)
{
  return _user_specified_atomic_smarts[zatom];
}

const IWString &
Smiles_Information::user_specified_atomic_smarts(atom_number_t zatom) const
{
  return _user_specified_atomic_smarts[zatom];
}

void
Smiles_Information::set_user_specified_atomic_smarts(atom_number_t zatom,
                                                const IWString & s)
{
  if (nullptr != _user_specified_atomic_smarts)
    ;
  else if (_natoms <= 0)
  {
    cerr << "Smiles_Information::set_user_specified_atomic_smarts:atom count unknown\n";
    return;
  }
  else
    _user_specified_atomic_smarts = new IWString[_natoms];

  _user_specified_atomic_smarts[zatom] = s;

  return;
}

int
Smiles_Information::contains_valid_ordering() const {
  return INVALID_SMILES_ORDER_TYPE != _smiles_order_type;
}

IWString
Molecule::isotopically_labelled_smiles()
{
  if (0 == _number_elements)
    return (".");

  int * isave = nullptr;

  for (int i = 0; i < _number_elements; i++)
  {
    int iso = _things[i]->isotope();

    if (0 == iso)
    {
      _things[i]->set_isotope(i);
      continue;
    }

    if (isave == nullptr)
      isave = new_int(_number_elements);

    isave[i] = iso;

    _things[i]->set_isotope(i);
  }

  Smiles_Information sminfo;

  (void) number_fragments();

  Atom_Chooser_Default acd;

  if (! _build_smiles_ordering_fctr(acd, nullptr, sminfo))
  {
    cerr << "Molecule::isotopically_labelled_smiles: cannot construct ordering\n";
    _smiles_information.set_error();
    return sminfo.smiles();
  }

  sminfo.set_smiles_order_type(DEFAULT_SMILES_ORDER_TYPE);

  _construct_smiles(_fragment_information, sminfo, nullptr);

  if (isave != nullptr)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->set_isotope(isave[i]);
    }

    delete [] isave;
  }
  else
  {
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->set_isotope(0);
    }
  }

  return (sminfo.smiles());
}

const IWString &
Molecule::smiles_using_order(const int * user_specified_rank)
{
  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  _smiles_information.invalidate();

  (void) number_fragments();

  _smiles_information.prepare_to_build_ordering(_number_elements);

  _smiles_information.set_smiles_order_type(USER_SPECIFIED_SMILES_ORDER);

  Atom_Chooser_Lowest_Rank aclr(user_specified_rank);

  int rc = _build_smiles_ordering_fctr(aclr, nullptr, _smiles_information);

  if (0 == rc)
    return _smiles_information.smiles();

  _smiles_information.set_smiles_order_type(USER_SPECIFIED_SMILES_ORDER);

  _construct_smiles(_fragment_information, _smiles_information, nullptr);

  return _smiles_information.smiles();
}

std::string
Molecule::Smiles() {
  const IWString& s = smiles();
  return std::string(s.data(), s.length());
}

std::string
Molecule::UniqueSmiles() {
  const IWString& s = unique_smiles();
  return std::string(s.data(), s.length());
}

std::string
Molecule::RandomSmiles() {
  const IWString& s = random_smiles();
  return std::string(s.data(), s.length());
}
