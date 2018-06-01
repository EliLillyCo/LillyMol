#include <iostream>
#include <iomanip>
#include <memory>
#include <limits>
#include <random>
#include <algorithm>

using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#include "primes.h"
#include "iwbits.h"
#include "misc.h"
#include "iwqsort.h"

#include "timsort.hpp"

#include "molecule.h"
#include "path.h"
#include "misc2.h"
#include "smiles.h"
#include "target.h"
#include "chiral_centre.h"

static int file_scope_include_isotopic_information_in_unique_smiles = 1;

void
set_include_isotopic_information_in_unique_smiles (int s)
{
  file_scope_include_isotopic_information_in_unique_smiles = s;

  return;
}

int
include_isotopic_information_in_unique_smiles () 
{
  return file_scope_include_isotopic_information_in_unique_smiles;
}

static int include_directional_bonding_information_in_unique_smiles = 1;

void
set_include_directional_bonding_information_in_unique_smiles (int s)
{
  include_directional_bonding_information_in_unique_smiles = s;

  return;
}

static int consider_implicit_hydrogens_in_unique_smiles = 1;

void
set_consider_implicit_hydrogens_in_unique_smiles(int s)
{
  consider_implicit_hydrogens_in_unique_smiles = s;
}

static int resolve_ties_by_geometry = 0;

void
set_resolve_unique_smiles_ties_by_geometry(int s)
{
  resolve_ties_by_geometry = s;
}

static int file_scope_consider_isotopes_as_zero_and_non_zero = 0;

void
set_consider_isotopes_as_zero_and_non_zero (int s)
{
  file_scope_consider_isotopes_as_zero_and_non_zero = s;

  if (0 == file_scope_consider_isotopes_as_zero_and_non_zero)   // turned off, no problems
    ;
  else if (file_scope_include_isotopic_information_in_unique_smiles)  // isotopes included, good
    ;
  else
    cerr << "set_consider_isotopes_as_zero_and_non_zero::warning, isotopic information not included in unique smiles\n";
}

int
consider_isotopes_as_zero_and_non_zero ()
{
  return file_scope_consider_isotopes_as_zero_and_non_zero;
}

static int file_scope_use_version_two_initial_rank_assignment = 0;

void
set_unique_determination_version(const int s)
{
  if (s > 1)
    file_scope_use_version_two_initial_rank_assignment = 1;
  else
    file_scope_use_version_two_initial_rank_assignment = 0;
}

static void
get_bonds (const Atom * a,
           int & aromatic_bonds,
           int & single_bonds,
           int & double_bonds,
           int & triple_bonds)
{
  aromatic_bonds = 0;
  single_bonds = 0;
  double_bonds = 0;
  triple_bonds = 0;

  const int acon = a->number_elements();

  Bond *const* rawb = a->rawdata();

  for (int i = 0; i < acon; i++)
  {
//  const Bond * b = a->item(i);
    const Bond * b = rawb[i];

    if (b->is_aromatic())
      aromatic_bonds++;
    else if (b->is_single_bond())
      single_bonds++;
    else if (b->is_double_bond())
      double_bonds++;
    else if (b->is_triple_bond())
      triple_bonds++;
    else if (IS_COORDINATION_BOND(b->btype()))     // clearly wrong, but I don't want to slow everything down for coordination bonds
      triple_bonds++;
    else
    {
      cerr << "What kind of bond is this!! " << (*b) << '\n';
      iwabort();
    }
  }

  return;
}

#define SEEMS_FASTER_WITH_INTEL
#ifdef SEEMS_FASTER_WITH_INTEL
static int
set_bond_score_in_ncon2 (const Atom * a,
                         Target_Atom * p)
{
  const int acon = a->ncon();

  Bond *const* rawb = a->rawdata();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = rawb[i];

    if (b->is_aromatic())
      rc += 1000;
    else if (b->is_single_bond())
      rc += 100;
    else if (b->is_double_bond())
      rc += 10;
    else if (b->is_triple_bond())
      rc++;
    else if (IS_COORDINATION_BOND(b->btype()))
      rc++;
  }
  
  p->set_ncon2(rc);

  return rc;
}
#endif

//#define SEEMS_FASTER_WITH_GCC
#ifdef SEEMS_FASTER_WITH_GCC
static int
set_bond_score_in_ncon2 (const Atom * a,
                         Target_Atom * p)
{
  int aromatic_bonds, single_bonds, double_bonds, triple_bonds;

  get_bonds(a, aromatic_bonds, single_bonds, double_bonds, triple_bonds);

  int rc = 1000 * aromatic_bonds + 100 * single_bonds + 10 * double_bonds + triple_bonds;

  p->set_ncon2(rc);

  return rc;
}
#endif

//#define DEBUG_TARGET_ATOM_COMPARITOR

/*
   The first part of the ranking is to order the atoms. We use a molecule target to
   do that.
*/

static int
target_atom_comparitor (Target_Atom * const * ppa1, Target_Atom * const * ppa2)
{
  Target_Atom *p1 = *ppa1;
  Target_Atom *p2 = *ppa2;

#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  assert (p1->ok());      // too expensive normally
  assert (p2->ok());
  atom_number_t t1 = p1->atom_number();
  atom_number_t t2 = p2->atom_number();
  cerr << "Comparing atom " << t1 << " (" << p1->m()->atomic_symbol(t1) << ", " << p1->ncon() << " connections) and " <<
                               t2 << " (" << p2->m()->atomic_symbol(t2) << ", " << p1->ncon() << " connections)\n";
//cerr << "Atomic numbers " << p1->m()->atomic_number(t1) << " and " << p2->m()->atomic_number(t2) << "\n";
#endif

  const Element * e1 = p1->element();
  const Element * e2 = p2->element();

// Lots of different possibilities here because elements may be NULL (excluded
// from the computation) or not from the periodic table (atomic number not valid)

  if (e1 == e2)    // covers case of both NULL
    ;
  else if (NULL == e1)
    return -1;
  else if (NULL == e2)
    return 1;
  else if (e1->atomic_number() < e2->atomic_number())
    return -1;
  else if (e1->atomic_number() > e2->atomic_number())
    return 1;
  else if (! e1->is_in_periodic_table() && ! e2->is_in_periodic_table())    // hash values really should be different
  {
    if (e1->atomic_symbol_hash_value() < e2->atomic_symbol_hash_value())
      return -1;
    else if (e1->atomic_symbol_hash_value() > e2->atomic_symbol_hash_value())
      return 1;
  }
  else if (! e1->is_in_periodic_table())
    return -1;
  else if (! e2->is_in_periodic_table())
    return 1;
  else if (e1->permanent_aromatic() == e2->permanent_aromatic())
    ;
  else if (e1->permanent_aromatic())
    return -1;
  else
    return 1;

  if (INVALID_ATOMIC_NUMBER == p1->atomic_number())    // these are atoms that have been excluded from consideration, they are always considered equivalent
    return 0;

  if (p1->ncon() < p2->ncon())
    return -1;
  else if (p1->ncon() > p2->ncon())
    return 1;

  if (p1->isotope() == p2->isotope())
    ;
  else if (! file_scope_include_isotopic_information_in_unique_smiles)
    ;
  else if (file_scope_consider_isotopes_as_zero_and_non_zero)
  {
    int i1 = p1->isotope() > 0;
    int i2 = p2->isotope() > 0;
    if (i1 < i2)
      return -1;
    else if (i1 > i2)
      return 1;
  }
  else if (p1->isotope() < p2->isotope())
    return -1;
  else
    return 1;

// Rings are strange because of the definition of the SSSR. For example,
// in cubane, all the atoms are equivalent initially, although some are
// in different numbers of rings.
// If either of the atoms is in 1 ring, then compare nrings, otherwise no

  if (p1->is_non_ring_atom() && p2->is_ring_atom())
    return -1;
  else if (p1->is_ring_atom() && p2->is_non_ring_atom())
    return 1;

  if (p1->formal_charge() == p2->formal_charge())
    ;
  else if (p1->formal_charge() < p2->formal_charge())
    return -1;
  else if (p1->formal_charge() > p2->formal_charge())
    return 1;

  Atom * a1 = const_cast<Atom *>(p1->atom());        // loss of const OK
  Atom * a2 = const_cast<Atom *>(p2->atom());        // loss of const OK

// Jan 2007. Restrict the case of ignoring implicit Hydrogens to heteroatoms
// only. Otherwise molecules have different unique smiles from their default

  if (! consider_implicit_hydrogens_in_unique_smiles)
    ;
  else if (a1->implicit_hydrogens() < a2->implicit_hydrogens())
    return -1;
  else if (a1->implicit_hydrogens() > a2->implicit_hydrogens())
    return 1;

// We store the bond scores in the target_atom's ncon2 attribute

  int nb1;

  if (! p1->ncon2_value_set())
    nb1 = set_bond_score_in_ncon2(a1, p1);
  else
    nb1 = p1->ncon2();
  
  int nb2;  
  if (! p2->ncon2_value_set())
    nb2 = set_bond_score_in_ncon2(a2, p2);
  else
    nb2 = p2->ncon2();

  if (nb1 < nb2)
    return -1;
  else if (nb1 > nb2)
    return 1;
    
#ifdef VERSION_CALLING_GET_BONDS_EACH_TIME
  int arom1, single1, double1, triple1;
  get_bonds(a1, arom1, single1, double1, triple1);

  int arom2, single2, double2, triple2;
  get_bonds(a2, arom2, single2, double2, triple2);

#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  cerr << "Not yet resolved, bonds are:";
  if (arom1 || arom2)
    cerr << " arom(" << arom1 << ',' << arom2 << ')';
  if (single1 | single2)
    cerr << " single (" << single1 << ',' << single2 << ')';
  if (double1 || double2)
    cerr << " double (" << double1 << ',' << double2 << ')';
  if (triple1 || triple2)
    cerr << " triple (" <<
          triple1 << ',' << triple2 << ')';
  cerr << '\n';
#endif

  if (arom1 < arom2)
    return -1;
  else if (arom1 > arom2)
    return 1;

  if (single1 < single2)
    return -1;
  else if (single1 > single2)
    return 1;
  
  if (double1 < double2)
    return -1;
  else if (double1 > double2)
    return 1;
  
  if (triple1 < triple2)
    return -1;
  else if (triple1 > triple2)
    return 1;
#endif
  
#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  cerr << "Atoms " << p1->atom_number() << " and " << p2->atom_number() << " not resolved yet, try rings\n";
#endif

// In the test for ring/non-ring above, we split the atoms if they were
// in a ring and out of a ring. At this stage, we cannot resolve chain
// atoms, so if one of them is in a chain, we know the other is as well.

  if (p1->is_non_ring_atom())      // only need to test one of them
    return 0;

// At this stage they are both in one or more rings. 

#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  cerr << "Resolving atoms " << p1->atom_number() << " and " << p2->atom_number() << " with rings\n";
#endif

  const List_of_Ring_Sizes * rs1 = p1->ring_sizes_including_non_sssr();
  const List_of_Ring_Sizes * rs2 = p2->ring_sizes_including_non_sssr();

// Need to be very careful here because of macrocycles. Because of the
// vagaries of ring perception, just check to make sure the smallest
// ring size of each is contained in the list of the other.

  int sr1 = rs1->item(0);
  int sr2 = rs2->item(0);

  if (sr1 < sr2)
  {
    if (! rs2->contains(sr1))
      return -1;

    if (! rs1->contains(sr2))
      return 1;
  }
  else if (sr1 > sr2)
  {
    if (! rs1->contains(sr2))
      return 1;

    if (! rs2->contains(sr1))
      return -1;
  }

#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  cerr << "chirlaity ? " << include_chiral_info_in_smiles() << " " << (NULL != p1->chiral_centre()) << " and " << (NULL != p2->chiral_centre()) << '\n';
#endif

  if (! include_chiral_info_in_smiles())
    ;
  else if (NULL == p1->chiral_centre() && NULL == p2->chiral_centre())
    ;
  else if (NULL != p1->chiral_centre() && NULL != p2->chiral_centre())
    ;
  else if (p1->chiral_centre())
    return -1;
  else
    return 1;

#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  cerr << "Atoms not resolved\n";
#endif

  return 0;
}

class Unique_Determination;


// We need to keep track of atoms and their associated rank.
// Also, each atom keeps track of its neighbours and the bonds
// via which they are attached

class Atom_and_Rank : public resizable_array<Atom_and_Rank *>
{
  friend
    class Unique_Determination;
  private:
    atom_number_t _a;
    const Atom *  _atom;
    unsigned int  _rank;

//  whenever it will fit into a 32 bit int, we compute the prime product, ala Weininger

    unsigned int _prime_product_of_neighbours;

//  The ranks of my neighbours

    int * _ranks_of_neighbours;

//  related to the bond types of the neighbours

    resizable_array<int> _bt;

//  Whenever we cannot compute the prime sum, we just store the neighbours ranks
//  and compare them. To facilitate that, we also store the sum

    int _sum_of_neighbour_ranks;

//  We keep track of when this atom has been identified as unique and is
//  awaiting processing.

    int _resolved_pending_processing;

//  During unique determination, the neighbours will change as they are resolved,
//  but the distance matrix computation requires knowledge of all atoms

    resizable_array<const Atom_and_Rank *> _connected_atoms;

    int _ncon;     // _ncon is just _connected_atoms.number_elements();

//  Chirality is a major complicating factor. We run the uniqueness
//  determination up to storing symmetry without chirality, then turn it on

    int _considering_chirality;

    const Chiral_Centre * _chiral_centre;

//  Based on the ranks of the neighbours, we can have a chirality score
//  that can differentiate stereoisomers

    int _chirality_score;

//  We may also get chirality indications from our chiral neighbours

    int _chiral_neighbours;

    resizable_array<int> _neighbour_chirality_score;

//  private functions

    void _default_values();

    int _index_of_atom(atom_number_t a) const;

    int _compare_by_chirality(const Atom_and_Rank * r2) const;
    int _compare_by_chirality_for_timsort(const Atom_and_Rank * r2) const;

    int _identify_two_unresolved_connections(const unsigned int * rank,
                                               atom_number_t & a1, 
                                               atom_number_t & i2) const;

  public:
    Atom_and_Rank(atom_number_t za, int zr, int nc);
    ~Atom_and_Rank();

    int ok () const;
    int debug_print (std::ostream &) const;

    int terse_print_neighbour_ranks (std::ostream &) const;

    int bt (int i) const { return _bt[i];}

    const Atom * atom () const { return _atom;}

    void set_chirality (const Chiral_Centre * c) { _chiral_centre = c;}
    const Chiral_Centre * chiral_centre () const { return _chiral_centre;}

    void set_consider_chirality (int s) { _considering_chirality = s;}

//  This function establishes the connectivity

    int is_connected_to (Atom_and_Rank *, const Bond *);

//  After we know about our connections, we can establish bond types

    void establish_neighbours ();

    int compare (const Atom_and_Rank *) const;
    int compare_for_timsort (const Atom_and_Rank *) const;
    int compare_less (const Atom_and_Rank *) const;

//  called whenever a neighboring atom has been classified, and is
//  removed from the computation

    int a_neighbour_has_been_classified (const Atom_and_Rank * c);

    unsigned int  rank () const { return _rank;}

    void set_rank (int r) { _rank = r;}

    int  resolved_pending_processing () const { return _resolved_pending_processing;}
    void set_resolved_pending_processing () { _resolved_pending_processing = 1;}

    int  set_rank (int r, int * unused_ranks);
    int  choose_an_unused_rank (int * unused_ranks, int n, int rank_delta);

    void set_atom_number (atom_number_t na) { _a = na;}
    atom_number_t atom_number () const { return _a;}

    int  ncon () const { return _connected_atoms.number_elements ();}

    void collect_neighbour_ranks ();

    int  is_connected_to (const Atom_and_Rank * z) const { return _connected_atoms.contains (z);}

//  When including chirality in a smiles, we need a means of getting a number
//  that reflects the chiral arrangement around an atom

    int  chiral_score (const unsigned int *) const;

    int  chirality_score () const { return _chirality_score;};
    void  compute_chirality_score (const unsigned int *);

    void propagate_chirality_influence (const unsigned int *) const;

//  int account_for_cis_trans_bonds(const Molecule &);
};

//#if defined(__GNUG__) || defined (__SUNPRO_CC)
template class resizable_array_p<Atom_and_Rank>;
template class resizable_array<Atom_and_Rank *>;
template class resizable_array<const Atom_and_Rank *>;
template class resizable_array_base<Atom_and_Rank *>;
template class resizable_array_base<const Atom_and_Rank *>;
//#endif

void
Atom_and_Rank::_default_values ()
{
  _ranks_of_neighbours = NULL;
  _sum_of_neighbour_ranks = 0;

  _chiral_centre = NULL;
  _chirality_score = 0;

  _prime_product_of_neighbours = 0;

  _resolved_pending_processing = 0;

  _ncon = 0;

  _chiral_neighbours = 0;

  _considering_chirality = 0;

  return;
}

Atom_and_Rank::Atom_and_Rank (atom_number_t za, int zr, int nc)
{
  _default_values();

  resize(nc);
  _bt.resize(nc);
  _connected_atoms.resize(nc);

  _a = za;
  _rank = zr;

  return;
}

Atom_and_Rank::~Atom_and_Rank ()
{
  delete [] _ranks_of_neighbours;
  _ranks_of_neighbours = NULL;

  return;
}

/*
  Scan the remaining connections and return the index of a given atom number
*/

int
Atom_and_Rank::_index_of_atom (atom_number_t a) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (a == _things[i]->atom_number())
      return i;
  }

  cerr << "Yipes, atom " << _a << " no remaining connection to " << a << '\n';
  abort();

  return -1;
}

int
Atom_and_Rank::set_rank (int new_rank, int * rank_in_use)
{
  rank_in_use[_rank]--;

  _rank = new_rank;
  rank_in_use[_rank]++;

  return 1;
}

/*
  The numbering used for the different bond types is arbitrary, all
  that is required is that the numbers be different for each type
*/

#define AR_AROMATIC_BOND 0
#define AR_SINGLE_BOND 1
#define AR_DOUBLE_BOND 2
#define AR_TRIPLE_BOND 3
#define AR_UP_BOND 4
#define AR_DOWN_BOND 5

int
Atom_and_Rank::choose_an_unused_rank (int * rank_in_use,
                                      int n,
                                      int rank_delta)
{
  if (1 == rank_in_use[_rank])    // it is already unique
    return 1;

  int new_rank = -1;
  for (int i = 0; i < n; i += rank_delta)
  {
    if (0 == rank_in_use[i])
    {
      new_rank = i;
      break;
    }
  }

  if (new_rank < 0)
  {
    cerr << "Yipes, no unused ranks\n";
    assert (NULL == "I don't know what to do");
    return 0;
  }

  rank_in_use[_rank]--;
  _rank = new_rank;
  rank_in_use[_rank] = 1;

  return 1;
}

/*
  Notify one atom that another atom is bonded to it.
*/

int
Atom_and_Rank::is_connected_to (Atom_and_Rank * n,
                                const Bond * b)
{
  assert (_number_elements == _bt.number_elements());

  add(n);

  if (b->is_aromatic())
    _bt.add(AR_AROMATIC_BOND);
  else if (b->is_single_bond())
  {
    if (b->is_directional() && include_directional_bonding_information_in_unique_smiles)
    {
      _bt.add(AR_DOWN_BOND);
/*    if (_a == b->a1())
      {
        if (b->is_directional_up())
          _bt.add(AR_UP_BOND);
        else
          _bt.add(AR_DOWN_BOND);
      }
      else
      {
        if (b->is_directional_up())
          _bt.add(AR_DOWN_BOND);
        else
          _bt.add(AR_UP_BOND);
      }*/
    }
    else
      _bt.add(AR_SINGLE_BOND);
  }
  else if (b->is_double_bond())
    _bt.add(AR_DOUBLE_BOND);
  else if (b->is_triple_bond())
    _bt.add(AR_TRIPLE_BOND);
  else
  {
    cerr << "What kind of bond is this!!! " << (*b) << '\n';
    iwabort();
  }

  _connected_atoms.add(n);
  _ncon = _connected_atoms.number_elements();

  if (! include_chiral_info_in_smiles())
    ;
  else if (NULL != n->chiral_centre())
    _chiral_neighbours++;

  return _number_elements;
}

//#define USING_QSORT

#define USING_TIMSORT 1

#if defined(USING_QSORT)
#elif defined(USING_TIMSORT)
#else
class Int_Comparitor_Larger_fo
{
  private:
  public:
    int operator() (const int & i1, const int & i2) const { if (i1 < i2) 
                                                              return -1;
                                                            else if (i1 > i2)
                                                              return 1;
                                                            else
                                                              return 0;
                                                           }
};

static Int_Comparitor_Larger_fo int_comparitor_larger_fo;
#endif

void
Atom_and_Rank::establish_neighbours ()
{
  assert (NULL == _ranks_of_neighbours);

  if (_number_elements)
    _ranks_of_neighbours = new int[_number_elements];
//  _ranks_of_neighbours = new_int(_number_elements);    // not necessary for it to be 0'd, aesthetics only

  return;
}

/*
  Single step in the Morgan-like algorithm
*/

void
Atom_and_Rank::collect_neighbour_ranks ()
{
  if (0 == _number_elements)
    return;

  _prime_product_of_neighbours = 0;

  unsigned int r0 = _bt[0] + _things[0]->rank();

  if (1 == _number_elements)
  {
    _prime_product_of_neighbours = primes[r0];
    return;
  }

  unsigned int r1 = _bt[1] + _things[1]->rank();

  if (r0 > r1)
    std::swap(r0, r1);

  if (2 == _number_elements)
  {
    if (r1 < IWNPRIMES)  // r1 is larger, both in range of primes array
      _prime_product_of_neighbours = primes[r0] * primes[r1];
    else
    {
      _sum_of_neighbour_ranks = r0 + r1;
      _ranks_of_neighbours[0] = r0;
      _ranks_of_neighbours[1] = r1;
    }

    return;
  }

  if (3 == _number_elements)
  {
    unsigned int r2 = _bt[2] + _things[2]->rank();

//  Ensure that r2 is the largest - remember, r0 and r1 are sorted

    if (r1 > r2)
      std::swap(r1, r2);

    assert (r0 <= r2 && r1 <= r2);

//  If this will fit in a 32 bit int, do it!
  
    if (r2 < 258)     // a heuristic, The 258'th prime is 1627 and 1619*1621*1627==2**32
    {
      _prime_product_of_neighbours = primes[r0] * primes[r1] * primes[r2];
      return;
    }

    if (r2 >= IWNPRIMES)   // only the first IWNPRIMES primes are in the header file
      ;
    else if (std::numeric_limits<unsigned int>::max() / primes[r2] < (primes[r0] * primes[r1]))
    {
      _prime_product_of_neighbours = primes[r0] * primes[r1] * primes[r2];
      return;
    }

//  Will not fit into 32 bit integer

    if (r0 < r1)
    {
      _ranks_of_neighbours[0] = r0;
      _ranks_of_neighbours[1] = r1;
    }
    else
    {
      _ranks_of_neighbours[0] = r1;
      _ranks_of_neighbours[1] = r0;
    }
    _ranks_of_neighbours[2] = r2;

    _sum_of_neighbour_ranks = r0 + r1 + r2;

    return;
  }

// When we have more then 3, we just record the elements

  _sum_of_neighbour_ranks = r0 + r1;
  if (r0 < r1)
  {
    _ranks_of_neighbours[0] = r0;
    _ranks_of_neighbours[1] = r1;
  }
  else
  {
    _ranks_of_neighbours[0] = r1;
    _ranks_of_neighbours[1] = r0;
  }

  for (int i = 2; i < _number_elements; i++)
  {
    _ranks_of_neighbours[i] = _bt[i] + _things[i]->rank();
    _sum_of_neighbour_ranks += _ranks_of_neighbours[i];
  }

#if defined(USING_QSORT)
  qsort(_ranks_of_neighbours, _number_elements, sizeof(int), (int (*) (const void *, const void *)) int_comparitor_larger);
#elif defined(USING_TIMSORT)
  gfx::timsort (_ranks_of_neighbours, _ranks_of_neighbours + _number_elements);
#else
  ::iwqsort (_ranks_of_neighbours, _number_elements, int_comparitor_larger_fo);
#endif

  return;
}

//#define DEBUG_A_NEIGHBOUR_HAS_BEEN_CLASSIFIED

/*
  One of the atoms attached to this one has been classified.
*/

int
Atom_and_Rank::a_neighbour_has_been_classified (const Atom_and_Rank * c)
{
#ifdef DEBUG_A_NEIGHBOUR_HAS_BEEN_CLASSIFIED
  cerr << "Neighbour " << c->atom_number() << " of atom " << _a << " has been classified\n";
  cerr << "I have " << _number_elements << " neighbours\n";
  for (int j = 0; j < _number_elements; j++)
  {
    const Atom_and_Rank * n = _things[j];
    cerr << "Atom " << n->atom_number() << " is a neighbour\n";
  }
#endif

  int i = index(const_cast<Atom_and_Rank *>(c));    // why do I need to cast away const, index takes a const object!!
  assert (i >= 0);

  remove_item(i);

  assert (_bt.ok_index(i));

  _bt.remove_item(i);

  return 1;
}

int
Atom_and_Rank::ok () const
{
  if (! resizable_array<Atom_and_Rank *>:: ok())
    return 0;

  if (_ncon < _number_elements)
    return 0;

  return 1;
}

int
Atom_and_Rank::debug_print (std::ostream & os) const
{
  os << "Atom " << _a << ", " << _number_elements << " neighbours, rank " << _rank << '\n';
  if (_number_elements < _ncon)
    os << _ncon << " connections within the molecule\n";

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom_and_Rank * n = _things[i];
    os << "Neighbour " << i << " is atom " << n->atom_number() << " btype " << _bt[i] << " rank " << n->rank() << endl;
  }

  if (_prime_product_of_neighbours)
    os << "Prime product of neighbours = " << _prime_product_of_neighbours << '\n';
  else if (_sum_of_neighbour_ranks)
  {
    os << "Sum of neighbours ranks = " << _sum_of_neighbour_ranks << '\n';
    for (int i = 0; i < _number_elements; i++)
    {
      os << " neighbour " << i << " rank = " << _ranks_of_neighbours[i] << '\n';
    }
  }

  return os.good();
}

int
Atom_and_Rank::terse_print_neighbour_ranks (std::ostream & os) const
{
  if (_prime_product_of_neighbours)
    os << "prime p = " << _prime_product_of_neighbours;
  else if (_sum_of_neighbour_ranks)
  {
    os << "neighbours";
    for (int i = 0; i < _number_elements; i++)
    {
      os << ' ' << _ranks_of_neighbours[i];
    }
  }

  return os.good();
}

#if defined(USING_QSORT)

static int
atom_and_rank_comparitor (Atom_and_Rank * const * pai1,
                          Atom_and_Rank * const * pai2)
{
  int r1 = (*pai1)->rank();
  int r2 = (*pai2)->rank();

  if (r1 < r2)
    return -1;
  else if (r1 == r2)
    return 0;
  else
    return 1;
}


#elif defined(USING_TIMSORT)

class Atom_and_Rank_Comparitor_Timsort
{
  private:
  public:
    int operator () (const Atom_and_Rank * pai1, const Atom_and_Rank * pai2) const
                                        {
                                          int r1 = pai1->rank();
                                          int r2 = pai2->rank();
                                          if (r1 < r2)
                                            return 1;
                                          return 0;
                                        }
};

static Atom_and_Rank_Comparitor_Timsort arc_timsort;

#else

class Atom_and_Rank_Comparitor
{
  private:
  public:
    int operator () (const Atom_and_Rank * pai1, const Atom_and_Rank * pai2) const
                                        {
                                          int r1 = pai1->rank();
                                          int r2 = pai2->rank();
                                          if (r1 < r2)
                                            return -1;
                                          else if (r1 > r2)
                                            return 1;
                                          else
                                            return 0;
                                        }
};

static Atom_and_Rank_Comparitor arc;
#endif


//#define DEBUG_COMPARE

int
Atom_and_Rank::compare (const Atom_and_Rank * r2) const
{
#ifdef DEBUG_COMPARE
  cerr << "Comparing " << _a << " with " << r2->_a << '\n';
  debug_print(cerr);
  cerr << "and\n";
  r2->debug_print(cerr);
#endif

// first check whether they are differentiated by their current rank

  if (_rank < r2->_rank)
    return -1;
  if (_rank > r2->_rank)
    return 1;

// They are not distinguished by their existing rank, check their neighbours

  if (_number_elements < r2->_number_elements)
    return -1;
  if (_number_elements > r2->_number_elements)
    return 1;

// They both have the same number of neighbours. Compare products, even
// if not computed

  if (_prime_product_of_neighbours < r2->_prime_product_of_neighbours)
    return -1;
  if (_prime_product_of_neighbours > r2->_prime_product_of_neighbours)
    return 1;

  if (0 == _prime_product_of_neighbours && 0 == r2->_prime_product_of_neighbours)   // neither one is computed
    ;
  else 
    return _compare_by_chirality(r2);

// prime products were not computed. Check the sums

  if (_sum_of_neighbour_ranks < r2->_sum_of_neighbour_ranks)
    return -1;
  if (_sum_of_neighbour_ranks > r2->_sum_of_neighbour_ranks)
    return 1;

  if (0 == _sum_of_neighbour_ranks)    // part of a subset, no neighbours computed
    return 0;

// At this stage, there is nothing left but an element by element comparison

  for (int i = 0; i < _number_elements; i++)
  {
    if (_ranks_of_neighbours[i] < r2->_ranks_of_neighbours[i])
      return -1;
    else if (_ranks_of_neighbours[i] > r2->_ranks_of_neighbours[i])
      return 1;
  }

#ifdef DEBUG_COMPARE
  cerr << "What about chirality " << _considering_chirality << '\n';
#endif

  return _compare_by_chirality(r2);
}

int
Atom_and_Rank::compare_for_timsort (const Atom_and_Rank * r2) const
{
#ifdef DEBUG_COMPARE
  cerr << "Comparing " << _a << " with " << r2->_a << '\n';
  debug_print(cerr);
  cerr << "and\n";
  r2->debug_print(cerr);
#endif

// first check whether they are differentiated by their current rank

  if (_rank < r2->_rank)
    return 1;
  if (_rank > r2->_rank)
    return 0;

// They are not distinguished by their existing rank, check their neighbours

  if (_number_elements < r2->_number_elements)
    return 1;
  if (_number_elements > r2->_number_elements)
    return 0;

// They both have the same number of neighbours. Compare products, even
// if not computed

  if (_prime_product_of_neighbours < r2->_prime_product_of_neighbours)
    return 1;
  if (_prime_product_of_neighbours > r2->_prime_product_of_neighbours)
    return 0;

  if (0 == _prime_product_of_neighbours && 0 == r2->_prime_product_of_neighbours)   // neither one is computed
    ;
  else 
    return _compare_by_chirality_for_timsort(r2);

// prime products were not computed. Check the sums

  if (_sum_of_neighbour_ranks < r2->_sum_of_neighbour_ranks)
    return 1;
  if (_sum_of_neighbour_ranks > r2->_sum_of_neighbour_ranks)
    return 0;

  if (0 == _sum_of_neighbour_ranks)    // part of a subset, no neighbours computed
    return 0;

// At this stage, there is nothing left but an element by element comparison

  for (int i = 0; i < _number_elements; i++)
  {
    if (_ranks_of_neighbours[i] < r2->_ranks_of_neighbours[i])
      return 1;
    else if (_ranks_of_neighbours[i] > r2->_ranks_of_neighbours[i])
      return 0;
  }

#ifdef DEBUG_COMPARE
  cerr << "What about chirality " << _considering_chirality << '\n';
#endif

  return _compare_by_chirality_for_timsort(r2);
}




//#define DEBUG_COMPARE_BY_CHIRAL_INFLUENCE

/*
  We've tried everything else to resolve two molecules, what about
  any chirality influences
*/

int
Atom_and_Rank::_compare_by_chirality (const Atom_and_Rank * r2) const
{
#ifdef DEBUG_COMPARE_BY_CHIRAL_INFLUENCE
  cerr << "Comparing chirality " << _considering_chirality << ", " << _chirality_score << " vs " << r2->_chirality_score << '\n';
#endif

  if (! _considering_chirality)
    return 0;

  if (_chirality_score == r2->_chirality_score)
    ;
  else if (_chirality_score < r2->_chirality_score)
    return -1;
  else
    return 1;

// Maybe resolved by configuration of neighbours

  int ncs = _neighbour_chirality_score.number_elements();

#ifdef DEBUG_COMPARE_BY_CHIRAL_INFLUENCE
  cerr << "Comparing " << ncs << " and " << r2->_neighbour_chirality_score.number_elements() << " neighbour chirality scores\n";
#endif

  if (ncs == r2->_neighbour_chirality_score.number_elements())
    ;
  else if (ncs > r2->_neighbour_chirality_score.number_elements())
    return -1;
  else
    return 1;

  if (0 == ncs)
    return 0;

  for (int i = 0; i < ncs; i++)
  {
    int n1 = _neighbour_chirality_score[i];
    int n2 = r2->_neighbour_chirality_score[i];

//  cerr << "   comparing neighbour chirality score " << n1 << " and " << n2 << '\n';

    if (n1 == n2)
      continue;

    if (n1 < n2)
      return -1;

    return 1;
  }

// If we come to here, we cannot resolve them.

  return 0;
}

int
Atom_and_Rank::_compare_by_chirality_for_timsort (const Atom_and_Rank * r2) const
{
#ifdef DEBUG_COMPARE_BY_CHIRAL_INFLUENCE
  cerr << "Comparing chirality " << _considering_chirality << ", " << _chirality_score << " vs " << r2->_chirality_score << '\n';
#endif

  if (! _considering_chirality)
    return 0;

  if (_chirality_score == r2->_chirality_score)
    ;
  else if (_chirality_score < r2->_chirality_score)
    return 1;
  else
    return 0;

// Maybe resolved by configuration of neighbours

  int ncs = _neighbour_chirality_score.number_elements();

#ifdef DEBUG_COMPARE_BY_CHIRAL_INFLUENCE
  cerr << "Comparing " << ncs << " and " << r2->_neighbour_chirality_score.number_elements() << " neighbour chirality scores\n";
#endif

  if (ncs == r2->_neighbour_chirality_score.number_elements())
    ;
  else if (ncs > r2->_neighbour_chirality_score.number_elements())
    return 1;
  else
    return 0;

  if (0 == ncs)
    return 0;

  for (int i = 0; i < ncs; i++)
  {
    int n1 = _neighbour_chirality_score[i];
    int n2 = r2->_neighbour_chirality_score[i];

//  cerr << "   comparing neighbour chirality score " << n1 << " and " << n2 << '\n';

    if (n1 == n2)
      continue;

    if (n1 < n2)
      return 1;

    return 0;
  }

// If we come to here, we cannot resolve them.

  return 0;
}

int
Atom_and_Rank::chiral_score (const unsigned int * rank) const
{
  assert (NULL != _chiral_centre);

  return _chiral_centre->orientation(rank);
}

/*
  We need to know the rank of an atom that is part of a chiral centre.
  Beware of lone pairs and implicit hydrogens
*/

static unsigned int
get_rank_for_connection (const unsigned int * rank, atom_number_t a)
{
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == a)
    return std::numeric_limits<unsigned int>::max(); 

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == a)
    return std::numeric_limits<unsigned int>::max() - 1;

  return  rank[a];
}

/*
  A chiral atom has been resolved. Does it have just two neighbours with
  the same rank. A1 and A2 are atom numbers
*/

int
Atom_and_Rank::_identify_two_unresolved_connections (const unsigned int * rank,
                                                     atom_number_t & a1, 
                                                     atom_number_t & a2) const
{
  unsigned int rank_top_front  = get_rank_for_connection(rank, _chiral_centre->top_front());
  unsigned int rank_top_back   = get_rank_for_connection(rank, _chiral_centre->top_back());
  unsigned int rank_left_down  = get_rank_for_connection(rank, _chiral_centre->left_down());
  unsigned int rank_right_down = get_rank_for_connection(rank, _chiral_centre->right_down());

  resizable_array<unsigned int> ranks;
  ranks.resize(4);

  ranks.add(rank_top_front);
  ranks.add_if_not_already_present(rank_top_back);
  ranks.add_if_not_already_present(rank_left_down);
  ranks.add_if_not_already_present(rank_right_down);

// We must have 3 distinct values

//cerr << "Found " << ranks.number_elements() << " distinct ranks\n";

  if (3 != ranks.number_elements())
    return 0;

// Now, identify the matched pair.

  a1 = INVALID_ATOM_NUMBER;     // assist checking later
  a2 = INVALID_ATOM_NUMBER;

  if (rank_top_front == rank_top_back)
  {
    a1 = _chiral_centre->top_front();
    a2 = _chiral_centre->top_back();
  }

  if (rank_top_front == rank_left_down)
  {
    a1 = _chiral_centre->top_front();
    a2 = _chiral_centre->left_down();
  }

  if (rank_top_front == rank_right_down)
  {
    a1 = _chiral_centre->top_front();
    a2 = _chiral_centre->right_down();
  }

  if (rank_top_back == rank_left_down)
  {
    a1 = _chiral_centre->top_back();
    a2 = _chiral_centre->left_down();
  }

  if (rank_top_back == rank_right_down)
  {
    a1 = _chiral_centre->top_back();
    a2 = _chiral_centre->right_down();
  }

  if (rank_left_down == rank_right_down)
  {
    a1 = _chiral_centre->left_down();
    a2 = _chiral_centre->right_down();
  }

  assert (a1 >= 0);
  assert (a2 >= 0);

  return 1;
}

//#define DEBUG_COMPUTE_CHIRALITY_SCORE

/*
  If the atom is non-chiral, its chirality score is 0.

  If chiral, call orientation . That returns -1, 0 or 1.
  In order to differentiate these numbers from 0, add 2
*/

void
Atom_and_Rank::compute_chirality_score (const unsigned int * rank)
{
  if (NULL == _chiral_centre)
    _chirality_score = 0;
  else
    _chirality_score = 2 + _chiral_centre->orientation(rank);

#ifdef DEBUG_COMPUTE_CHIRALITY_SCORE
  cerr << "Atom " << _a << " chirality score " << _chirality_score << " chiral neighbours? " << _chiral_neighbours << '\n';
#endif

  if (0 == _chiral_neighbours)
    return;

  _neighbour_chirality_score.resize_keep_storage(0);

  for (int i = 0; i < _ncon; i++)
  {
    const Atom_and_Rank * ari = _connected_atoms[i];

    const Chiral_Centre * c = ari->chiral_centre();

    if (NULL == c)
      continue;

    int s = c->influence(rank, _a);

#ifdef DEBUG_COMPUTE_CHIRALITY_SCORE
    cerr << "influence on atom " << _a << " from atom " << ari->atom_number() << " is " << s << '\n';
#endif

    if (0 == s)
      continue;

    _neighbour_chirality_score.insert_in_order(s);
  }

  return;
}

class Unique_Determination : public resizable_array_p<Atom_and_Rank>
{
  private:
    int _matoms;
    int _ma9;      // the size of the _rank_in_use array

    int _nactive;   // the number which have not yet been assigned a canonical rank

//  We assign canonical id's in decreasing order

    int _next_canonical_rank_to_assign;

//  We keep track of whether or not symmetry has been completely stored

    int _symmetry_stored;

//  Symmetry classes are assigned in increasing order

    int _next_symmetry_class_to_assign;

//  If chirality is included in smiles, it is helpful to keep track of the
//  number of chiral atoms remaining to be processed

    int _nchiral;

//  But we don't actually use chirality until symmetry has been perceived

    int _use_chirality;

//  to avoid repeated calls to include_chiral_info_in_smiles(), we store
//  the value here

    int _include_chiral_info_in_smiles;

//  the number of cis-trans bonds remaining

    int _cis_trans_bonds;

    Molecule * _m;    // should be const, but some functions needed are non const

//  For chirality scoring, we need an array of ranks

    unsigned int * _rank;

//  We need to keep track of the ranks which have been assigned
//  so we can get an unused one when we need it.

    int * _rank_in_use;

//  We use a different rank_delta if cis_trans bonds are present. This
//  is to preserve compatability with existing unique smiles. If that
//  isn't required, we could just use (AR_DOWN_BOND + 1)

    int _rank_delta;

//  The actual canonical rank

    int * _canonical_rank;

//  As we do the unique determination, we also perceive symmetry.

    int * _symmetry;

//  As we iterate, we need to keep track of the previous ranks of each atom
//  so we can determine when the ranks are stable.

    unsigned int * _old_rank;

//  We also need a cross reference between atom number and Atom_and_Rank

    Atom_and_Rank ** _atom_xref;

//  With larger molecules, we run into serious problems sorting a partially
//  sorted array, so we shuffle

    std::minstd_rand _rng;

//  private functions;
  
    int _allocate_atom_arrays (int, int);
    int _free_all_arrays ();

    int _initialise (Molecule &, const int *);

    void _assign_initial_ranks();
    void _assign_initial_ranks_v2();
    void _assign_initial_ranks (const int * include_atom);
    void _initialise_rank_in_use ();
    void _reassign_ranks ();
    void _expand_around_cis_trans_bonds();

//  Internal function to perform the analysis up to the point of full symmetry perception

    int _perceive_symmetry ();

//  Have the ranks changed from one iteration to the next

    int _ranks_changed ();

//  Does a comparison of stored ranks for the extended neighbour list of
//  a set of atoms.

    int _rank_values_resolved (resizable_array<atom_number_t> &, int);

//  Called only from _single_step_process_unique_atoms. 

    void _atom_is_unique (atom_number_t a);

//  when computing chirality influences, we need the _rank array filled

    void _fill_rank_array_for_chirality ();

//  Once we have perceived symmetry, we start to include the influence of chirality.

    void _turn_on_chirality_considerations_and_reassign_ranks ();

//  Scan through the atom_and_rank array and identify atoms with unique ranks

    int _get_indices_of_unique_atoms (resizable_array<int> & unique_atoms) const;

//  These functions are used when unique atoms are identified in the
//  sorted array of ranks

    int _single_step_process_unique_atoms ();

//  If atoms become disconnected, they are automatically done

    int _process_all_now_disconnected_atoms ();

//  Once the ranking remains unchanged from one iteration to the next,
//  we store the resulting symmetry information

    int _store_symmetry_info ();

//  Use this to break ties

    int _break_a_tie ();

//  If all else fails, we can break ties by geometry

    int _choose_tie_breaker_by_geometry() const;

    template <typename T> int _identify_extreme_value(T & c) const;

    int _choose_tie_breaker_atom ();

//  Fill the _chirality_score array

    void _compute_chirality_scores ();

//  Look at a sequence of equivalent atoms and analyse their chirality scores

    void _analyse_chirality_in_sequence (int sstart, int & next_starting_position, int * cs, int & chiral_atoms_in_sequence) const;
    void _identify_next_sequence (int sstart, int & next_starting_position, int & chiral_atoms_in_sequence) const;

//  An atom has been classified as unique. Move it to the inactive part of the array

    int _move_to_inactive (int);

    int _process_all_unique_atoms ();

    int _identify_unused_rank () const;
    int _identify_two_unused_ranks (int & r1, int & r2) const;
    int _identify_some_unused_ranks (int ranks_needed,
                      resizable_array<int> & ranks_identified) const;

//  When an atom is ranked, we change the rank values for those atoms connected

    int __adjust_rank_of_atoms_attached_to (Atom_and_Rank * r);
    int _adjust_rank_of_atoms_attached_to (Atom_and_Rank * r);

//  A single step of the expansion

    int _expand (int);

    int _canonical_order (int);

    int _index_of_atom (atom_number_t zatom) const;

    int _get_rank (atom_number_t a) const;
    int _index_if_active (atom_number_t a) const;
    void _adjust_initial_ranks_for_cis_trans_bonds ();
    void _adjust_initial_ranks_for_cis_trans_bonds (atom_number_t a1, atom_number_t a2);
    unsigned int _compute_cis_trans_rank (Atom_and_Rank * a, 
                       atom_number_t nw,
                       atom_number_t sw,
                       atom_number_t a2,
                       atom_number_t ne,
                       atom_number_t se);

    int _expand_around_cis_trans_bond (atom_number_t a1,
                                atom_number_t a2,
                                unsigned int * new_rank);

    int _identify_directionally_attached_atoms (const Atom_and_Rank * ar,
                                atom_number_t & nw,
                                atom_number_t & sw) const;
    int _identify_directionally_attached_bonds (const Atom_and_Rank * ar,
                                                const Bond * & b1,
                                                const Bond * & b2) const;
  public:
    Unique_Determination ();
    Unique_Determination (int);
    ~Unique_Determination ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int symmetry (int, int *) const;    // return previously computed values
    int symmetry (Molecule *, int *, const int *);

    int canonical_order (Molecule &, int *, const int *);
};

//#define DEBUG_UNIQUE_DETERMINATION

Unique_Determination::Unique_Determination ()
{
  _matoms = 0;
  _ma9 = 0;
  _nactive = 0;

  _rank = NULL;

  _rank_in_use = NULL;

  _symmetry = NULL;
  _canonical_rank = NULL;

  _old_rank = NULL;

  _atom_xref = NULL;

  _use_chirality = 0;

  _include_chiral_info_in_smiles = include_chiral_info_in_smiles();

// We leave gaps in the rank array so we can distinguish different
// ranks attached via different bond types. We need a different rank
// if cis-trans bonds are present.

  _rank_delta = AR_TRIPLE_BOND + 1;

  std::random_device rd;

  _rng.seed(3172776704);    // just used for shuffling, so proper randomisation not important

  return;
}

/*
  To cut down on the number of allocations, we allocate one large array and have
  the sub-arrays be part of that.
*/

int
Unique_Determination::_allocate_atom_arrays (int matoms,
                                             int size_ma9)
{
  assert (matoms > 0);

  _matoms = matoms;

  _ma9 = size_ma9;

  assert (sizeof(int) == sizeof(unsigned int));

  _rank = new unsigned int[_matoms + _ma9 + _matoms + _matoms + _matoms];

  _rank_in_use    = reinterpret_cast<int *>(_rank + _matoms);
  _canonical_rank = reinterpret_cast<int *>(_rank + _matoms + _ma9);
  _symmetry       = reinterpret_cast<int *>(_rank + _matoms + _ma9 + _matoms);
  _old_rank       =                         _rank + _matoms + _ma9 + _matoms + _matoms;
  
  std::fill_n(_canonical_rank, _matoms, -1);
  std::fill_n(_symmetry, _matoms, -1);

  _atom_xref = new Atom_and_Rank *[_matoms];

  return 1;
}

Unique_Determination::~Unique_Determination ()
{
  _free_all_arrays();

  return;
}

int
Unique_Determination::_free_all_arrays ()
{
  if (NULL != _rank)
  {
    delete [] _rank;
//  delete [] _rank_in_use;
//  delete []_canonical_rank;
//  delete []_symmetry;
//  delete []_old_rank;

    delete []_atom_xref;
  }

  return 1;
}

int
Unique_Determination::ok () const
{
  if (NULL == _m)
    return 0;

  if (! _m->ok())
    return 0;

  return 1;
}

int
Unique_Determination::debug_print (std::ostream & os) const
{
  os << "Details on Unique_Determination object, for " << _number_elements << " atoms, " << _nactive << " active\n";
  os << "use chirality " << _use_chirality << " include_chiral_info_in_smiles " << _include_chiral_info_in_smiles << endl;

  for (int i = 0; i < _nactive; i++)
  {
    const Atom_and_Rank * ar = _things[i];
    atom_number_t a = ar->atom_number();

    os << "I = " << std::setw(3) << i << " atom " << std::setw(3) << a << " (" << _m->atomic_symbol(a) << ") rank = " <<
           ar->rank() << " (" << _rank_in_use[ar->rank()] << ')';
    os << " ncon = " << ar->number_elements();
    if (_symmetry && _symmetry[a] >= 0)    // equals 0 is illegal
      os << " symmetry = " << _symmetry[ar->atom_number()];

    os << ' ';

    if (ar->chirality_score())
      os << " CHIRAL " << ar->chirality_score();

    os << ' ';
    ar->terse_print_neighbour_ranks(os);

    os << '\n';
  }

  return os.good();
}

//#define DEBUG_STORE_SYMMETRY

/*
  The ordering has remained constant from one iteration to the next. Store
  symmetry classes
*/

int
Unique_Determination::_store_symmetry_info ()
{
#ifdef DEBUG_STORE_SYMMETRY
  cerr << "Storing symmetry classes starting with " << _next_symmetry_class_to_assign << '\n';
  debug_print(cerr);
#endif

  unsigned int rprev = 0;

  for (int i = 0; i < _nactive; i++)
  {
    Atom_and_Rank * r = _things[i];

    if (r->rank() != rprev)
      _next_symmetry_class_to_assign++;

    atom_number_t a = r->atom_number();
    _symmetry[a] = _next_symmetry_class_to_assign;

    rprev = r->rank();
  }

  return 1;
}

void
Unique_Determination::_compute_chirality_scores ()
{
  _fill_rank_array_for_chirality();

  for (int i = 0; i < _nactive; i++)
  {
    _things[i]->compute_chirality_score(_rank);
  }

  return;
}

/*
  Analyse the atoms starting at SSTART (descending). For each atom
  with the same rank as SSTART, increment the CS array with its chirality score

  We also accumulate the number of chiral atoms we visit in this grouping
*/

void
Unique_Determination::_identify_next_sequence (int sstart,
                                               int & next_starting_position,
                                               int & chiral_atoms_in_sequence) const
{
  unsigned int zrank = _things[sstart]->rank();    // keep processing while rank is the same

  next_starting_position = sstart;

  chiral_atoms_in_sequence = 0;

  while (next_starting_position >= 0)
  {
    const Atom_and_Rank * ari = _things[next_starting_position];

    if (zrank != ari->rank())    // finished with the grouping of equivalent atoms
      return;

    if (ari->chiral_centre())
      chiral_atoms_in_sequence++;

    next_starting_position--;
  }

  return;
}

template <typename T>
int
Unique_Determination::_identify_extreme_value(T & c) const
{
  coord_t qmin = c(_things[0]->atom());
  coord_t qmax = qmin;
  int which_min = 0;
  int which_max = 0;

  for (int i = 1; i < _nactive; i++)
  {
    coord_t q = c(_things[i]->atom());

    if (q > qmax)
    {
      qmax = q;
      which_max = i;
    }
    else if (q == qmax)
      which_max = -1;
    else if (q < qmin)
    {
      qmin = q;
      which_min = i;
    }
    else if (q == qmin)
      which_min = -1;
  }

  if (which_max >= 0)
    return which_max;

  if (which_min >= 0)
    return which_min;

  return -1;
}

class Xfetcher
{
  private:
  public:
    coord_t operator() (const Atom * a) const { return a->x();}
};

class Yfetcher
{
  private:
  public:
    coord_t operator() (const Atom * a) const { return a->y();}
};

class Zfetcher
{
  private:
  public:
    coord_t operator() (const Atom * a) const { return a->z();}
};

template int Unique_Determination::_identify_extreme_value(Xfetcher &) const;
template int Unique_Determination::_identify_extreme_value(Yfetcher &) const;
template int Unique_Determination::_identify_extreme_value(Zfetcher &) const;

/*
  There are no other means of differentiating atoms. Choose one by geometry
*/

int
Unique_Determination::_choose_tie_breaker_by_geometry() const
{
  if (! resolve_ties_by_geometry)
    return _nactive - 1;

  Xfetcher xf;
  int i = _identify_extreme_value(xf);
  if (i >= 0)
    return i;

  Yfetcher yf;
  i = _identify_extreme_value(yf);
  if (i >= 0)
    return i;

  Zfetcher zf;
  i = _identify_extreme_value(zf);
  if (i >= 0)
    return i;

  return _nactive - 1;
}

//#define DEBUG_CHOOSE_TIE_BREAKER_ATOM

/*
  We need to break a tie. For all ranks, there are at least two atoms with
  the same rank.
  Our strategy is to first break the rank of a chiral atom
*/

int
Unique_Determination::_choose_tie_breaker_atom ()
{
#ifdef DEBUG_CHOOSE_TIE_BREAKER_ATOM
  cerr << "Choosing tie breaker atom, _nchiral = " << _nchiral << '\n';
#endif

  if (0 == _nchiral)      // no chiral atoms, just return the last atom on the list
    return _choose_tie_breaker_by_geometry();

  if (! _include_chiral_info_in_smiles)
    return _choose_tie_breaker_by_geometry();

  int chiral_atoms_still_active = 0;

  for (int i = 0; i < _nactive; i++)
  {
    if (_things[i]->chiral_centre())
    {
      chiral_atoms_still_active = 1;
      break;
    }
  }

  if (0 == chiral_atoms_still_active)
  {
#ifdef DEBUG_CHOOSE_TIE_BREAKER_ATOM
    cerr << "No chiral atoms remaining\n";
#endif

    return _choose_tie_breaker_by_geometry();
  }

// If we don't find something that can be resolved, we'll break the sequence
// with the smallest number of chiral atoms in it

  int shortest_sequence = _nactive;
  int start_of_shortest_sequence = _nactive - 1;

  int next_starting_position = _nactive - 1;

  while (next_starting_position > 0)
  {
    int sstart = next_starting_position;
    int chiral_atoms_in_sequence;

    _identify_next_sequence(sstart, next_starting_position, chiral_atoms_in_sequence);

//#define DEBUG_SEQUENCE_FOUND
#ifdef DEBUG_SEQUENCE_FOUND
    cerr << "Sequence contains " << chiral_atoms_in_sequence << " starting with " << sstart << ",";
    for (int j = 0; j < 4; j++)
    {
      cerr << ' ' << cs[j];
    }
    cerr << '\n';
#endif

    if (chiral_atoms_in_sequence > 0 && chiral_atoms_in_sequence < shortest_sequence)
    {
      shortest_sequence = chiral_atoms_in_sequence;
      start_of_shortest_sequence = sstart;
    }
  }

// If we come to here, we weren't able to resolve things by chirality. Let's
// break something in the shortest sequence

#ifdef DEBUG_CHOOSE_TIE_BREAKER_ATOM
  cerr << "Shortest chiral containing sequence starts at " << start_of_shortest_sequence << '\n';
#endif

  for (int i = start_of_shortest_sequence; i >= 0; i--)
  {
    if (_things[i]->chirality_score())
      return i;
  }

// should not come to here

#ifdef DEBUG_CHOOSE_TIE_BREAKER_ATOM
  cerr << "Hmmm, choosing last atom\n";
#endif

  return _nactive - 1;
}

//#define DEBUG_BREAK_A_TIE

int
Unique_Determination::_break_a_tie ()
{
  int t = _choose_tie_breaker_atom();

  Atom_and_Rank * r = _things[t];
  atom_number_t a = r->atom_number();

#ifdef DEBUG_BREAK_A_TIE
  cerr << "Atom " << a << " chosen as tie breaker\n";
#endif

  _adjust_rank_of_atoms_attached_to(r);

  _atom_is_unique(a);

  _move_to_inactive(t);

  _reassign_ranks();

#ifdef DEBUG_BREAK_A_TIE
  cerr << "After breaking tie\n";
  debug_print(cerr);
#endif

  return 1;
}

int
Unique_Determination::_identify_unused_rank () const
{
  for (int i = 0; i < _ma9; i += _rank_delta)
  {
    if (0 == _rank_in_use[i])
      return i;
  }

  cerr << "Yipes, no ranks available\n";
  debug_print(cerr);
  assert (NULL == "Cannot continue");

  return -1;
}

int
Unique_Determination::_identify_two_unused_ranks (int & r1, int & r2) const
{
  int rc = 0;
  for (int i = 0; i < _ma9; i += _rank_delta)
  {
    if (0 == _rank_in_use[i])
    {
      if (0 == rc)
      {
        r1 = i;
        rc++;
      }
      else
      {
        r2 = i;
        return 1;
      }
    }
  }

  cerr << "Yipes, cannot find two unused ranks\n";
  debug_print(cerr);

  return 0;
}

int
Unique_Determination::_identify_some_unused_ranks (int ranks_needed,
                      resizable_array<int> & ranks_identified) const
{
  ranks_identified.resize(ranks_needed);

  int rc = 0;
  for (int i = 0; i < _ma9; i += _rank_delta)
  {
    if (0 == _rank_in_use[i])
    {
      ranks_identified.add(i);
      rc++;
      if (rc == ranks_needed)
        return rc;
    }
  }

  cerr << "Yipes, cannot find " << ranks_needed << " unused ranks, got " << ranks_identified.number_elements() << endl;
  debug_print(cerr);
  iwabort();

  return 0;
}

//#define DEBUG_ADJUST_RANK_OF_ATOMS_ATTACHED

/*
  Atom r->atom_number () has been assigned a unique number.  We must
  now adjust all its neighbours to reflect the fact that it is now "different"

  We compute a new rank for each neighbour, based on an offset (basically
  a number which will be unique for each newly unique atom), the bond
  type and the current rank for the neighbour.

  We assign new ranks to these neighbours, but must make sure that we
  assign those ranks in a definite order.
*/

int
Unique_Determination::__adjust_rank_of_atoms_attached_to (Atom_and_Rank * r)
{

// scan through the neighbours

  int neighbours_of_r = r->number_elements();

#ifdef DEBUG_ADJUST_RANK_OF_ATOMS_ATTACHED
  cerr << "Adjusting ranks for " << neighbours_of_r << " atoms attached to atom " << r->atom_number() << '\n';
#endif

  if (0 == neighbours_of_r)     // atom has already been isolated
    return 1;

// treat the case of one neighbour specially

  if (1 == neighbours_of_r)
  {
    Atom_and_Rank * n = r->item(0);
    n->choose_an_unused_rank(_rank_in_use, _ma9, _rank_delta);

#ifdef DEBUG_ADJUST_RANK_OF_ATOMS_ATTACHED
    cerr << "One neighbour, atom " << n->atom_number() << " assigned rank " << n->rank() << '\n';
#endif

    return 1;
  }

// Two neighbours is probably pretty common too

  if (2 == neighbours_of_r)
  {
    Atom_and_Rank * n0 = r->item(0);
    Atom_and_Rank * n1 = r->item(1);

    int r0 = n0->rank() + r->bt(0);
    int r1 = n1->rank() + r->bt(1);

    if (r0 == r1)
    {
      unsigned int oldrank = n0->rank();
      if (2 != _rank_in_use[oldrank])
        n0->choose_an_unused_rank(_rank_in_use, _ma9, _rank_delta);

      n1->set_rank(n0->rank(), _rank_in_use);
    }
    else
    {
      int nr0, nr1;
      (void) _identify_two_unused_ranks(nr0, nr1);

      if (r0 < r1)
      {
        n0->set_rank(nr0, _rank_in_use);
        n1->set_rank(nr1, _rank_in_use);
      }
      else
      {
        n0->set_rank(nr1, _rank_in_use);
        n1->set_rank(nr0, _rank_in_use);
      }
    }

#ifdef DEBUG_ADJUST_RANK_OF_ATOMS_ATTACHED
    cerr << "Two neighbours, atom " << n0->atom_number() << " assigned " << n0->rank() <<
            " atom " << n1->atom_number() << " assigned " << n1->rank() << '\n';
#endif

    return 1;
  }

// Three or more neighbours is done the hard way, but remember, there will
// hardly ever be more than 4, so no quick sort below

  resizable_array<unsigned int> old_ranks;
  old_ranks.resize_keep_storage(neighbours_of_r);

  for (int i = 0; i < neighbours_of_r; i++)
  {
    Atom_and_Rank * ni = r->item(i);

//  cerr << " i = " << i << " bt " << r->bt(i) << endl;

    int orank = r->bt(i) + ni->rank();
    old_ranks.insert_in_order_if_not_already_present(orank);
  }

#ifdef DEBUG_ADJUST_RANK_OF_ATOMS_ATTACHED
  cerr << "Sorted old ranks of neighbours are";
  for (int i = 0; i < old_ranks.number_elements(); i++)
  {
    cerr << ' ' << old_ranks[i];
  }

  cerr << '\n';
#endif

// Now convert all the old ranks to new ones, in order

  int nr = old_ranks.number_elements();
  if (1 == nr)     // only one unique value - neighbours still unresolved
  {
    if (neighbours_of_r == _rank_in_use[old_ranks[0]])   // these are the only instances of a rank. Don't change them
      return 1;

    int new_rank = _identify_unused_rank();

    for (int i = 0; i < neighbours_of_r; i++)
    {
      Atom_and_Rank * n = r->item(i);
      n->set_rank(new_rank, _rank_in_use);
    }

    return 1;
  }

// All the others are done the hard way

  resizable_array<int> new_ranks;
  (void) _identify_some_unused_ranks(nr, new_ranks);

  for (int i = 0; i < nr; i++)
  {
    unsigned int orank = old_ranks[i];
    int new_rank = new_ranks[i];
    for (int j = 0; j < neighbours_of_r; j++)
    {
      Atom_and_Rank * n = r->item(j);
      if (orank == n->rank() + r->bt(j))
      {
        n->set_rank(new_rank, _rank_in_use);
      }
    }
  }

  return 1;
}

int
Unique_Determination::_adjust_rank_of_atoms_attached_to (Atom_and_Rank * r)
{
  int rc = __adjust_rank_of_atoms_attached_to(r);

// Let r's neighbours know that it has been classified

  int n = r->number_elements();

  for (int i = 0; i < n; i++)
  {
//  cerr << "Letting atom " << r->item(i)->atom_number() << " know about it\n";
    r->item(i)->a_neighbour_has_been_classified(r);
  }

#ifdef DEBUG_ADJUST_RANK_OF_ATOMS_ATTACHED
  cerr << "After making adjustments\n";
  debug_print(cerr);
#endif

  return rc;
}

//#define DEBUG_ATOM_IS_UNIQUE

/*
  Atom A has been determined to be unique.
*/

void
Unique_Determination::_atom_is_unique (atom_number_t a)
{
#ifdef DEBUG_ATOM_IS_UNIQUE
  cerr << "Atom " << a << " has been determined to be unique";
  if (_symmetry[a] < 0)
    cerr << ". Symmetry " << _next_symmetry_class_to_assign;
  cerr << '\n';
#endif

  assert (_canonical_rank[a] < 0);

  if (_symmetry[a] < 0)    // if it has not already been done
    _symmetry[a] = _next_symmetry_class_to_assign++;

  _canonical_rank[a] = _next_canonical_rank_to_assign--;

  return;
}

/*
  Scan through the array and identify the indices of atoms with
  unique ranks
*/

int
Unique_Determination::_get_indices_of_unique_atoms (resizable_array<int> & unique_atoms) const
{
  if (1 == _nactive)
  {
    unique_atoms.add(0);
    return 1;
  }

  unsigned int prev_rank = _things[0]->rank();
  int count = 1;     // the number of instances of a rank value

  for (int i = 1; i < _nactive; i++)
  {
    Atom_and_Rank * aii = _things[i];

    if (aii->rank() == prev_rank)    // same as the one before
    {
      count++;
      continue;
    }

    assert (aii->rank() > prev_rank);

//  Handle the cases where I is greater than I - 1

    if (1 == count)                   // there was only one of the previous rank
    {
      unique_atoms.add(i - 1);      // the previous atom is unique
      if (_nactive - 1 == i)               // the last atom in the list is also unique
        unique_atoms.add(i);
    }
    else if (_nactive - 1 == i)            // the last atom is different from the one before
      unique_atoms.add(i);
  
    count = 1;
    prev_rank = aii->rank();
  }

  return unique_atoms.number_elements();
}

//#define DEBUG_MOVE_TO_INACTIVE

/*
  Atom_and_Rank object number A has been classified as unique. We must now move
  it to the inactive section.
*/

int
Unique_Determination::_move_to_inactive (int a)
{
#ifdef DEBUG_MOVE_TO_INACTIVE
  cerr << "Item number " << a << " is to be made inactive\n";
#endif

  assert (_nactive);
  if (1 == _nactive)
  {
    _nactive = 0;
    return 1;
  }

  Atom_and_Rank * r = _things[a];

  for (int i = a; i < _nactive - 1; i++)
  {
    _things[i] = _things[i + 1];
  }

  _nactive--;
  _things[_nactive] = r;

#ifdef DEBUG_MOVE_TO_INACTIVE
  cerr << "After moving inactive item\n";
  debug_print(cerr);
#endif

  return 1;
}

/*
   We are doing something with chirality, so we need to have all the
   ranks assembed in the _rank array
*/

void
Unique_Determination::_fill_rank_array_for_chirality ()
{
  for (int i = 0; i < _matoms; i++)     // _canonical_rank will only be partly filled
  {
    _rank[i] = _canonical_rank[i];
//  cerr << "Assigned " << _rank[i] << " to atom " << i << '\n';
  }
  

// For the still active atoms, put in some number guaranteed not to
// collide with anything in the array

  for (int i = 0; i < _nactive; i++)
  {
    const Atom_and_Rank * ari = _things[i];

    atom_number_t a = ari->atom_number();
    assert (a >= 0 && a < _matoms);

    _rank[a] = _matoms + 1 + ari->rank();   // will be different from other values present
//  cerr << "active: assigned " << _rank[a] << " to atom " << a << '\n';
  }

  return;
}

//#define DEBUG_SINGLE_STEP_PROCESS_UNIQUE_ATOMS

/*
  Step through the array of ranked atoms - assumed to be sorted
*/

int
Unique_Determination::_single_step_process_unique_atoms ()
{
#ifdef DEBUG_SINGLE_STEP_PROCESS_UNIQUE_ATOMS
  cerr << "At beginning of _single_step_process_unique_atoms\n";
  debug_print(cerr);
#endif

  resizable_array<int> unique_atoms;    // we record their index in the array 
  unique_atoms.resize_keep_storage(_matoms);        // lets be optimistic

  int nu = _get_indices_of_unique_atoms(unique_atoms);

// If we got any unique atoms, we must adjust the ranks of their neighbours and
// ensure that the list is still sorted

#ifdef DEBUG_SINGLE_STEP_PROCESS_UNIQUE_ATOMS
  cerr << "Found " << nu << " unique atoms\n";
  for (int i = 0; i < nu; i++)
  {
    int j = unique_atoms[i];
    const Atom_and_Rank * r = _things[j];

    cerr << " i = " << i << " j = " << j << " atom " << r->atom_number() << " rank " << r->rank() << '\n';
  }
#endif

  if (0 == nu)
   return nu;

// If we are including chirality, we need to initialise the _rank array

  if (! _include_chiral_info_in_smiles)
    ;
  else if (_use_chirality)
    _fill_rank_array_for_chirality();

  for (int i = nu - 1; i >= 0; i--)    // we will be removing elements
  {
    int j = unique_atoms[i];

    Atom_and_Rank * r = _things[j];
    atom_number_t   a = r->atom_number();

    _adjust_rank_of_atoms_attached_to(r);

    _atom_is_unique(a);

    _move_to_inactive(j);
  }

// For now, re-sort, but later fix to just move the ones changed

  _reassign_ranks();     // does a sort and renumbers

#ifdef DEBUG_SINGLE_STEP_PROCESS_UNIQUE_ATOMS
  cerr << "At end of single_step_process_unique_atoms\n";
  debug_print(cerr);
#endif

  return nu;
}

//#define DEBUG_REASSIGN_RANKS

void
Unique_Determination::_reassign_ranks ()
{
#ifdef DEBUG_REASSIGN_RANKS
  cerr << "Reassigning ranks. Current ranks are\n";
  debug_print(cerr);
#endif

  _initialise_rank_in_use();

  if (! _include_chiral_info_in_smiles)
    ;
  else if (_use_chirality)
    _compute_chirality_scores();

#if defined(USING_QSORT)
  std::shuffle (_things, _things + _nactive, _rng);
  qsort(_things, _nactive, sizeof(Atom_and_Rank *), (int (*) (const void *, const void *)) atom_and_rank_comparitor);
#elif defined(USING_TIMSORT)
  gfx::timsort(_things, _things + _nactive, arc_timsort);
#else
  if (_nactive > 20)
    std::shuffle (_things, _things + _nactive, _rng);
  ::iwqsort (_things, _nactive, arc);
#endif

  unsigned int rprev = 0;

  int rank_to_assign = 0;

  for (int i = 0; i < _nactive; i++)
  {
    Atom_and_Rank * t = _things[i];
    if (t->rank() != rprev)
    {
      rank_to_assign += _rank_delta;
      rprev = t->rank();
    }

    t->set_rank(rank_to_assign);
    _rank_in_use[rank_to_assign]++;
  }

#ifdef DEBUG_REASSIGN_RANKS
  cerr << "Ranks for " << _nactive << " atoms reassigned\n";
  for (int i = 0; i < _nactive; i++)
  {
    const Atom_and_Rank * r = _things[i];
    cerr << "i = " << i << " atom " << r->atom_number() << " rank " << r->rank() << '\n';
  }
#endif

  return;
}

//#define DEBUG_PROCESS_ALL_NOW_DISCONNECTED_ATOMS

/*
  During processing atoms can have all their neighbours classified.
  When that happens, they too can be classified
*/

int
Unique_Determination::_process_all_now_disconnected_atoms ()
{
  resizable_array<int> unique_atoms;     // actually indices in the array
  unique_atoms.resize_keep_storage(_nactive);

  for (int i = 0; i < _nactive; i++)
  {
    Atom_and_Rank * r = _things[i];
    if (0 == r->number_elements())
      unique_atoms.add(i);
  }

  int nu = unique_atoms.number_elements();

#ifdef DEBUG_PROCESS_ALL_NOW_DISCONNECTED_ATOMS

  cerr << "Found " << nu << " now disconnected atoms\n";
  for (int i = 0; i < nu; i++)
  {
    int j = unique_atoms[i];
    const Atom_and_Rank * r = _things[j];
    cerr << "j = " << j << " atom " << r->atom_number() << " current rank " << r->rank() << '\n';
  }
#endif

  if (0 == nu)
    return 0;

// Ran into a problem with symmetries. Consider CF3. After processing the C
// atom, we have three disconnected atoms. These are passed here. But, if
// we call _atom_is_unique, they will be assigned an arbitrary symmetry

// Oct 2000. Is this correct? What if there are two separate CF3 groups in the
// molecule? Look into this sometime...

  unsigned int prev_rank;    // don't worry, no uninitialised access - see below
  int prev_sym = -1;
  for (int i = nu - 1; i >= 0; i--)
  {
    int j = unique_atoms[i];
    Atom_and_Rank * r = _things[j];
    atom_number_t   a = r->atom_number();

    int need_to_store_symmetry = 0;
    if (_symmetry[a] < 0 && (i < nu - 1) &&  r->rank() == prev_rank)
      need_to_store_symmetry = 1;

    if (need_to_store_symmetry)
      _symmetry[a] = prev_sym;

    _atom_is_unique(a);

    _move_to_inactive(j);

    prev_rank = r->rank();
    prev_sym = _symmetry[a];
  }

// No need to resort the array, as all we did was remove items

#ifdef DEBUG_PROCESS_ALL_NOW_DISCONNECTED_ATOMS
  cerr << "At end of process_all_now_disconnected atoms\n";
  debug_print(cerr);
#endif

  return nu;
}

int
Unique_Determination::_process_all_unique_atoms ()
{
  int rc = 0;      // the number of unique atoms we find
  int found_each_iteration = 0;
  while ((found_each_iteration = _single_step_process_unique_atoms()))
  {
    rc += found_each_iteration;
    if (0 == _nactive)
      return rc;
  }

  return rc;
}

/*
  At each step, we check to see if the ranks have changed. As we check, we store
  the old ranks
*/

int
Unique_Determination::_ranks_changed ()
{
  int rc = 0;
  for (int i = 0; i < _nactive; i++)
  {
    const Atom_and_Rank * r = _things[i];
    atom_number_t a = r->atom_number();

    if (_old_rank[a] != r->rank())
    {
      _old_rank[a] = r->rank();
      rc++;
    }
  }

  return rc;
}

//#undef USING_TIMSORT

#if defined(USING_QSORT)
static int
atom_and_neighbour_rank_comparitor (Atom_and_Rank * const * pai1,
                                    Atom_and_Rank * const * pai2)
{
  Atom_and_Rank * a1 = *pai1;
  Atom_and_Rank * a2 = *pai2;

  return a1->compare(a2);
}
#elif defined(USING_TIMSORT)
class Atom_and_Neighbour_Rank_Comparitor_Timsort
{
  private:
  public:
    int operator() (const Atom_and_Rank * pai1, const Atom_and_Rank * pai2) const { return (pai1)->compare_for_timsort(pai2);}
};

static Atom_and_Neighbour_Rank_Comparitor_Timsort atom_and_neighbour_rank_comparitor_timsort;

#else

class Atom_and_Neighbour_Rank_Comparitor
{
  private:
  public:
    int operator() (const Atom_and_Rank * pai1, const Atom_and_Rank * pai2) const { return (pai1)->compare(pai2);}
};

static Atom_and_Neighbour_Rank_Comparitor atom_and_neighbour_rank_comparitor;

#endif


//#define DEBUG_EXPAND

/*
  A single step of our process involves having each atom collect the ranks
  from its neighbours, and then sorting the resulting list. We then
  re-assign the ranks.

  This is called during each step of processing, but also when chirality
  is introduced. When chirality is introduced, it won't be necessary
  to call collect_neighbour_ranks, so that's why we have an argument
*/

int
Unique_Determination::_expand (int collect_neighbours)
{
#ifdef DEBUG_EXPAND
  cerr << "At start of expand\n";
  debug_print(cerr);
#endif

  if (! _include_chiral_info_in_smiles)
    ;
  else if (_use_chirality)
    _compute_chirality_scores();

  if (collect_neighbours)
  {
    for (int i = 0; i < _nactive; i++)
    {
      Atom_and_Rank * r = _things[i];
      (void) r->collect_neighbour_ranks();
    }
  }

#ifdef DEBUG_EXPAND
  cerr << "After collecting neighbour ranks\n";
  debug_print(cerr);
#endif

// we must now sort them according to their neighbour ranks

#if defined(USING_QSORT)
  std::shuffle (_things, _things + _nactive, _rng);
  qsort(_things, _nactive, sizeof(Atom_and_Rank *), (int (*) (const void *, const void *)) atom_and_neighbour_rank_comparitor);
#elif defined(USING_TIMSORT)
  gfx::timsort(_things, _things + _nactive, atom_and_neighbour_rank_comparitor_timsort);
#else
  if (_nactive > 10)
    std::shuffle (_things, _things + _nactive, _rng);
  ::iwqsort (_things, _nactive, atom_and_neighbour_rank_comparitor);
#endif

#ifdef DEBUG_EXPAND
  cerr << "After sorting based on neighbours\n";
  debug_print(cerr);
#endif

  _initialise_rank_in_use();

// and assign ranks according to the neighbours rank. The only trickery here
// is that we must assign the new rank for the previous atom after we do the
// comparison, because atoms are ordered primarily by their current rank value

  int next_rank_to_assign = 0;

  for (int i = 1; i < _nactive; i++)
  {
    Atom_and_Rank * rp = _things[i - 1];
    Atom_and_Rank * r  = _things[i];

#if defined(USING_TIMSORT)
    int tmp = rp->compare_for_timsort(r);
#else
    int tmp = r->compare(rp);
#endif

#ifdef DEBUG_EXPAND
    cerr << " i = " << i << " rank " << r->rank() << " rp rank " << rp->rank() << " tmp " << tmp << endl;
#endif

    rp->set_rank(next_rank_to_assign);   // change its rank after we invoke compare()
    _rank_in_use[next_rank_to_assign]++;

    assert (tmp >= 0);

    if (tmp > 0)
      next_rank_to_assign += _rank_delta;
  }

  _things[_nactive - 1]->set_rank(next_rank_to_assign);
  _rank_in_use[next_rank_to_assign]++;

#ifdef DEBUG_EXPAND
  cerr << "After assigning new ranks\n";
  debug_print(cerr);
#endif

  return 1;
}

//#define DEBUG_ASSIGN_INITIAL_RANKS

/*
  Both the canonical order and symmetry functions need to initialise the
  Unique_Determination object for a given molecule.
*/

int
Unique_Determination::_initialise (Molecule & m,
                                   const int * include_atom)
{
  _m = &m;
  _matoms = m.natoms();

// Handle the two special cases specifically

  if (0 == _matoms)
    return 1;


  if (include_directional_bonding_information_in_unique_smiles && m.cis_trans_bonds_present())
  {
    _rank_delta = AR_DOWN_BOND + 1;
    (void) _allocate_atom_arrays(_matoms, _matoms * 25);   
  }
  else
    (void) _allocate_atom_arrays(_matoms, _matoms * 9);   

  if (1 == _matoms)
  {
    _canonical_rank[0] = 1;
    _symmetry[0] = 1;
    return 1;
  }

  set_vector(_old_rank, _matoms, static_cast<unsigned int>(0));

  _next_canonical_rank_to_assign = _matoms;
  _next_symmetry_class_to_assign = 1;
  _symmetry_stored = 0;

  _nchiral = 0;

  _cis_trans_bonds = 0;

  if (NULL != include_atom)
    _assign_initial_ranks(include_atom);
  else if (file_scope_use_version_two_initial_rank_assignment)
    _assign_initial_ranks_v2();
  else
    _assign_initial_ranks();

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  for (int i = 0; i < _matoms; i++)
  {
    _atom_xref[i] = NULL;
  }
#endif

  return 1;
}

//#define DEBUG_SYMMETRY

int
Unique_Determination::_perceive_symmetry ()
{
  return 1;
}

int
Unique_Determination::symmetry (Molecule * m, int * symmetry_class,
                                const int * include_atom)
{
  assert (OK_MOLECULE (m));
  assert (symmetry_class);

  _initialise(*m, include_atom);

  if (_matoms < 2)
    return 1;

  int rc = _perceive_symmetry();

#ifdef DEBUG_SYMMETRY
  cerr << "After symmetry perception\n";
#endif

  for (int i = 0; i < _matoms; i++)
  {
    symmetry_class[i] = _symmetry[i];

#ifdef DEBUG_SYMMETRY
    cerr << "Atom " << i << " symmetry class " << _symmetry[i] << '\n';
#endif
  }

  return rc;
}

int
Unique_Determination::symmetry (int matoms, int * symmetry_class) const
{
  assert (ok());
  assert (matoms == _matoms);    // some slight assurance that the right molecule is calling us

  copy_vector(symmetry_class, _symmetry, _matoms);

  return 1;
}

//#define DEBUG_CANONICAL_ORDER

int
Unique_Determination::_canonical_order (int stop_when_symmetry_perceived)
{
  assert (_nactive);

  _symmetry_stored = 0;

  int iterations = 0;
  while (_nactive > 0)
  {
    iterations++;

#ifdef DEBUG_CANONICAL_ORDER
    cerr << "Unique_Determination::_canonical_order: iteration " << iterations << " nactive = " << _nactive << '\n';
#endif

    if (iterations > 1 && ! _ranks_changed())
    {
      if (_cis_trans_bonds > 0)
      {
        _expand_around_cis_trans_bonds();
        continue;
      }

      if (0 == _symmetry_stored)
      {
        _store_symmetry_info();
        _symmetry_stored = 1;
        if (stop_when_symmetry_perceived)
          break;

//      cerr << "ranks not changed, chirality " << _nchiral << '\n';

        if (_nchiral && _include_chiral_info_in_smiles)
        {
          _turn_on_chirality_considerations_and_reassign_ranks();
          continue;
        }
      }
      
      _break_a_tie();
    }

    if (_process_all_unique_atoms())
    {
      if (0 == _nactive)
        break;
    }

    if (iterations > 1 && _process_all_now_disconnected_atoms())
    {
      if (0 == _nactive)
        break;
    }

    _expand(1);     // 1 means collect neighbour ranks
    if (_cis_trans_bonds > 0)
      _expand_around_cis_trans_bonds();

#ifdef DEBUG_CANONICAL_ORDER
    cerr << "At completion of iteration\n";
    debug_print(cerr);
#endif
  }

  if (0 == _symmetry_stored)
    _store_symmetry_info();

#ifdef DEBUG_CANONICAL_ORDER
  cerr << "Canonical order determined\n";
  for (int i = 0; i < _matoms; i++)
  {
    const Atom_and_Rank * r = _things[i];

    atom_number_t a = r->atom_number();

    cerr << " i = " << i << " atom " << a << ' ' << _m->const_smarts_equivalent_for_atom(a) << " rank " << r->rank() << '\n';
  }
#endif

  return 1;
}

/*
  We have just got to the stage where symmetry is perceived. We now need to 
  start including chirality into considerations
*/

void
Unique_Determination::_turn_on_chirality_considerations_and_reassign_ranks()
{
//cerr << "Chirality now being considered\n";

  for (int i = 0; i < _nactive; i++)
  {
    _things[i]->set_consider_chirality(1);
  }

  _use_chirality = 1;

  _expand(0);    // 0 means don't collect neighbour ranks again

//cerr << "_turn_on_chirality_considerations_and_reassign_ranks\n";
//debug_print(cerr);

  return;
}

void
Unique_Determination::_expand_around_cis_trans_bonds()
{
  int matoms = _m->natoms();

  unsigned int * rank_delta = new unsigned int[matoms]; std::unique_ptr<unsigned int[]> free_rank_delta(rank_delta);

  set_vector(rank_delta, matoms, static_cast<unsigned int>(0));

  int ne = _m->nedges();

  for (int i = 0; i < ne && _cis_trans_bonds > 0; i++)
  {
    const Bond * b = _m->bondi(i);

    if (! b->is_double_bond())
      continue;

    if (! b->part_of_cis_trans_grouping())
      continue;

    _expand_around_cis_trans_bond (b->a1(), b->a2(), rank_delta);

 // _cis_trans_bonds--;
  }

  for (int i = 0; i < _nactive; i++)
  {
//  if (rank_delta[i] > 0)
//    cerr << "Item " << i << " new rank_delta " <<rank_delta[i] << endl;
    if (rank_delta[i] > 0)
      _things[i]->set_rank(_things[i]->rank() + rank_delta[i]);
  }

// Need to turn these off until I get the basic algorithm working properly...

  _cis_trans_bonds = 0;    // should already be zero, maybe not..

#ifdef DEBUG_EXPAND_AROUND_CIS_TRANS_BONDS
  cerr << "Expanded " << _cis_trans_bonds << " CT bonds, new ranks\n";
  for (int i = 0; i < _nactive; i++)
  {
    cerr << " atom " << i << " rank " << _rank[i] << endl;
  }
#endif

  _reassign_ranks();

#ifdef DEBUG_EXPAND_AROUND_CIS_TRANS_BONDS
  cerr << "After expanding around " << _cis_trans_bonds << " cis trans bonds, ranks are\n";

  for (int i = 0; i < _nactive; i++)
  {
    cerr << " atom " << i << " rank " << _rank[i] << endl;
  }
#endif

  return;
}

int
Unique_Determination::_expand_around_cis_trans_bond(atom_number_t a1,
                                atom_number_t a2,
                                unsigned int * rank_delta)
{
  int ndx1 = _index_if_active(a1);
  if (ndx1 < 0)
    return 0;

  int ndx2 = _index_if_active(a2);
  if (ndx2 < 0)
    return 0;

  atom_number_t nw, sw;

  if (! _identify_directionally_attached_atoms(_things[ndx1], nw, sw))
    return 0;

  atom_number_t ne, se;

  if (! _identify_directionally_attached_atoms(_things[ndx2], ne, se))
    return 0;

  rank_delta[ndx1] += _compute_cis_trans_rank(_things[a1], nw, sw, a2, ne, se);
  rank_delta[ndx2] += _compute_cis_trans_rank(_things[a2], ne, se, a1, nw, sw);

//cerr << "Centre atom deltas " << rank_delta[ndx1] << " and " << rank_delta[ndx2] << endl;

  const Bond * bnw, * bsw;
  if (! _identify_directionally_attached_bonds(_things[a1], bnw, bsw))
    return 0;

  const Bond * bne, * bse;
  if (! _identify_directionally_attached_bonds(_things[a2], bne, bse))
    return 0;

//int matoms = _m->natoms();

  int rnw;
//atom_number_t anw;
  if (NULL != bnw)
  {
    nw = bnw->other(a1);
    rnw = _get_rank(nw);
  }
  else
  {
    nw = INVALID_ATOM_NUMBER;
    rnw = -1;
  }

  int rsw;
//atom_number_t asw;
  if (NULL != bsw)
  {
    se = bsw->other(a1);
    rsw = _get_rank(se);
  }
  else
  {
    sw = INVALID_ATOM_NUMBER;
    rsw = -1;
  }

  int rne;
//atom_number_t ane;
  if (NULL != bne)
  {
    ne = bne->other(a1);
    rne = _get_rank(ne);
  }
  else
  {
    ne = INVALID_ATOM_NUMBER;
    rne = -1;
  }

  int rse;
//atom_number_t ase;
  if (NULL != bse)
  {
    se = bse->other(a1);
    rse = _get_rank(se);
  }
  else
  {
    se = INVALID_ATOM_NUMBER;
    rse = -1;
  }

//int ra1 = _get_rank(a1);
//int ra2 = _get_rank(a2);

  if (rnw >= 0 && rne >= 0)
  {
    rank_delta[a1] += 2 * (rnw + rne);
    rank_delta[a2] += 2 * (rnw + rne);
  }

  if (rsw >= 0 && rse >= 0)
  {
    rank_delta[a1] += 2 * (rsw + rse);
    rank_delta[a2] += 2 * (rsw + rse);
  }

  if (rnw >= 0 && rse >= 0)
  {
    rank_delta[a1] += 3 * (rnw + rse);
    rank_delta[a2] += 3 * (rnw + rse);
  }

  if (rne >= 0 && rsw >= 0)
  {
    rank_delta[a1] += 3 * (rne + rsw);
    rank_delta[a2] += 3 * (rne + rsw);
  }

  return 1;
}

int
Unique_Determination::_identify_directionally_attached_atoms (const Atom_and_Rank * ar,
                                atom_number_t & nw,
                                atom_number_t & sw) const
{
  nw = INVALID_ATOM_NUMBER;
  sw = INVALID_ATOM_NUMBER;

  atom_number_t zatom = ar->atom_number();

  const Atom * a = _m->atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    if (b->is_double_bond())
      continue;

    if (! b->is_directional())
      continue;

    atom_number_t j = b->other(zatom);

    int k = _index_if_active(j);
    if (k < 0)
      continue;

    if (b->is_directional_up())
    {
      if (zatom == b->a1())
        nw = j;
      else
        sw = j;
    }
    else if (b->is_directional_down())
    {
      if (zatom == b->a1())
        sw = j;
      else
        nw = j;
    }
  }

  if (INVALID_ATOM_NUMBER == nw && INVALID_ATOM_NUMBER == sw)
    return 0;

  return 1;
}

int
Unique_Determination::_identify_directionally_attached_bonds (const Atom_and_Rank * ar,
                                                const Bond * & bnw,
                                                const Bond * & bsw) const
{
  bnw = NULL;
  bsw = NULL;

  atom_number_t zatom = ar->atom_number();

  const Atom * a = _m->atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    if (b->is_double_bond())
      continue;

    if (! b->is_directional())
      continue;

    atom_number_t j = b->other(zatom);

    if (_index_if_active(j) < 0)    // seems unlikely
      continue;

    if (b->is_directional_up())
    {
      if (zatom == b->a1())
        bnw = b;
      else
        bsw = b;
    }
    else if (b->is_directional_down())
    {
      if (zatom == b->a1())
        bsw = b;
      else
        bnw = b;
    }
  }

  if (NULL != bnw)
    return 1;

  if (NULL != bsw)
    return 1;

  return 0;
}

/*int
Atom_and_Rank::account_for_cis_trans_bonds(const Molecule & m)
{
  if (_number_elements < 2 || _number_elements > 3)
    return 0;

  const Atom * a = m.atomi(_a);

  int acon = a->ncon();

  if (acon != _number_elements)   // cannot process otherwise
    return 0;

  if (a->nbonds() - 1 != acon)
    return 0;

  const Bond * b1 = NULL;
  const Bond * b2 = NULL;
  atom_number_t west = INVALID_ATOM_NUMBER;

  int b1rank = 0;
  int b2rank = 0;
  int west_rank = -1;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
    {
      if (! b->part_of_cis_trans_grouping())
        return 0;
      west = b->other(_a);
      west_rank = _ranks_of_neighbours[i];
    }
    else if (! b->is_directional())
      return 0;
    else if (NULL == b1)
    {
      b1 = b;
      b1rank = _ranks_of_neighbours[i];
    }
    else
    {
      b2 = b;
      b2rank = _ranks_of_neighbours[i];
    }
  }

  if (NULL == b1)   // how could that happen?
    return 0;

  atom_number_t ne = INVALID_ATOM_NUMBER;
  atom_number_t se = INVALID_ATOM_NUMBER;

  if (b1->is_directional_up())
  {
    if (_a == b1->a1())
      ne = b1->other(_a);
    else
      se = b1->other(_a);
  }
  else if (b1->is_directional_down())
  {
    if (_a == b1->a1())
      se = b1->other(_a);
    else
      ne = b1->other(_a);
  }

  if (NULL == b2)
    ;
  else if (b2->is_directional_up())
  {
    if (_a == b2->a1())
      ne = b2->other(_a);
    else
      se = b2->other(_a);
  }
  else if (b2->is_directional_down())
  {
    if (_a == b2->a1())
      se = b2->other(_a);
    else
      ne = b2->other(_a);
  }

  _rank += west_rank;
  if (se >= 0)
  {
    int sei = _index_of_atom(se);
    _rank += 100 * _ranks_of_neighbours[sei];
  }

  if (ne >= 0)
  {
    int nei = _index_of_atom(ne);
    _rank += 601 * _ranks_of_neighbours[nei];
  }

  return 1;
}*/

void
Unique_Determination::_initialise_rank_in_use ()
{
  set_vector(_rank_in_use, _ma9, 0);

  return;
}

#define COMPILE_IN_FASTER_VERSION
#ifdef COMPILE_IN_FASTER_VERSION

/*
  for the initial sort, we collect as much info about the atom as we can into one of these objects.
  We have three invariants, which are computed different ways, so as to lessen the probability
  of collisions
*/

class Atom_Info
{
  private:
    uint64_t _invariant[3];
    atom_number_t _atom_number;
    resizable_array<int> _list_of_ring_sizes;

  public:

    int initialise(Molecule & m, const atom_number_t zatom);
    int initialise(Molecule & m, const atom_number_t zatom, const Chiral_Centre * c);

    atom_number_t atom_number() const { return _atom_number;}

    int compare(const Atom_Info & rhs) const;

    template <typename O> int debug_print(O & output) const;
};

#ifdef AT_HAS_CONSTRUCTOR
Atom_Info::Atom_Info()
{
}
#endif

template <typename O>
int
Atom_Info::debug_print(O & output) const
{
  output << "Atom_Info::debug_print:atom " << _atom_number;
  for (int i = 0; i < 3; ++i)
  {
    output << ' ' << _invariant[i];
  }

  for (int i = 0; i < _list_of_ring_sizes.number_elements(); ++i)
  {
    output << ' ' << _list_of_ring_sizes[i];
  }

  output << '\n';

  return 1;
}

int
Atom_Info::compare(const Atom_Info & rhs) const
{
//#define DEBUG_ATOM_INFO_CMP
#ifdef DEBUG_ATOM_INFO_CMP
  cerr << "CMP  ";
  debug_print(cerr);
  cerr << "WITH ";
  rhs.debug_print(cerr);
#endif

  for (int i = 0; i < 3; ++i)
  {
    if (_invariant[i] == rhs._invariant[i])
      ;
    else if (_invariant[i] < rhs._invariant[i])
      return -1;
    else
      return 1;
  }

//cerr << "Atom_Info::compare:not resolve by invariants\n";

  const int n1 = _list_of_ring_sizes.number_elements();
  const int n2 = rhs._list_of_ring_sizes.number_elements();

  if (0 == n1 && 0 == n2)    // neither atom in a ring
    return 0;
  
  if (0 == n1)     // probably never happens because rbc is part of the invariant
    return -1;
  if (0 == n2)     // probably never happens because rbc is part of the invariant
    return 1;

#ifdef DEBUG_RING_SIZE
  cerr << "RS1";
  for (int i = 0; i< n1; ++i)
  {
    cerr << ' ' << _list_of_ring_sizes[i];
  }
  cerr << endl;
  cerr << "RS2";
  for (int i = 0; i< n2; ++i)
  {
    cerr << ' ' << rhs._list_of_ring_sizes[i];
  }
  cerr << endl;
#endif

//if (n1 < n2)
//  return -1;
//else if (n1 > n2)
//  return 1;

  const int sr1 = _list_of_ring_sizes[0];
  const int sr2 = rhs._list_of_ring_sizes[0];

  if (sr1 < sr2)
  {
    if (! rhs._list_of_ring_sizes.contains(sr1))
      return -1;

    if (! _list_of_ring_sizes.contains(sr2))
      return 1;
  }
  else if (sr1 > sr2)
  {
    if (! _list_of_ring_sizes.contains(sr2))
      return 1;

    if (! rhs._list_of_ring_sizes.contains(sr1))
      return -1;
  }

  return 0;

  if (! _list_of_ring_sizes.contains(rhs._list_of_ring_sizes[0]))
    return -1;
  if (! rhs._list_of_ring_sizes.contains(_list_of_ring_sizes[0]))
    return 1;
  return 0;


  int n = n1 <= n2 ? n1 : n2;

  for (int i = 0; i < n; ++i)
  {
    if (_list_of_ring_sizes[i] < rhs._list_of_ring_sizes[i])
      return -1;
    if (_list_of_ring_sizes[i] > rhs._list_of_ring_sizes[i])
      return  1;
  }

  return 0;
}

static int
bond_and_atom_score(Molecule & m,
                    const atom_number_t zatom,
                    int & rbc,
                    uint64_t & zsum)
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  rbc = 0;

  Bond *const* rawb = a->rawdata();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = rawb[i];

    const atom_number_t j = b->other(zatom);

    const Atom * aj = m.atomi(j);

    if (b->is_aromatic())
    {
      rc += 1000;
      zsum += 11 * aj->element()->atomic_symbol_hash_value() + 2 * aj->ncon();
      rbc++;
    }
    else if (b->is_single_bond())
    {
      rc += 100;
      zsum += 3 * aj->element()->atomic_symbol_hash_value() + 2 * aj->ncon();
      if (b->nrings())
        rbc++;
    }
    else if (b->is_double_bond())
    {
      rc += 10;
      zsum += 5 * aj->element()->atomic_symbol_hash_value() + 2 * aj->ncon();
      if (b->nrings())
        rbc++;
    }
    else if (b->is_triple_bond())
    {
      rc++;
      zsum += 7 * aj->element()->atomic_symbol_hash_value() + 2 * aj->ncon();
      if (b->nrings())
        rbc++;
    }
    else if (IS_COORDINATION_BOND(b->btype()))
    {
      rc++;
      zsum += 13 * aj->element()->atomic_symbol_hash_value() + 2 * aj->ncon();
      if (b->nrings())
        rbc++;
    }
  }
  
  return rc;
}

int
Atom_Info::initialise(Molecule & m,
                      const atom_number_t zatom,
                      const Chiral_Centre * c)
{
  _atom_number = zatom;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  int rbc = 0;
  uint64_t zsum = 0;

  const int bs = bond_and_atom_score(m, zatom, rbc, zsum);

  _invariant[0] = a->element()->atomic_symbol_hash_value() + 82828282 * zsum;

  const int arom = m.is_aromatic(zatom);

  uint64_t x1 = 10000000000 * acon + 288578282 * rbc                + 26012422 * a->formal_charge() - 1942 * bs;
  uint64_t x2 = 21928000000 * rbc  + 911971781 * a->formal_charge() + 10970988 * acon + 19204 * bs;

  x1 -= 9076 * arom;
  x2 += 58 * arom;

  if (0 == a->isotope())
    ;
  else if (include_isotopic_information_in_unique_smiles())
  {
    x1 += 571 * a->isotope();
    x2 += 82122 * a->isotope();
  }

  if (consider_implicit_hydrogens_in_unique_smiles)
  {
    const int ih = m.implicit_hydrogens(zatom);

    x1 += 1972 * ih;
    x2 -= 82221 * ih;
  
//  if (a->implicit_hydrogens_known())    bad idea
//    x2 -= 5;
  }

  if (NULL != c)
  {
    x1 -= 821 * bs;
    x2 += 741 * acon;
  }

  if (0 == rbc)
  {
    _invariant[1] = x1;
    _invariant[2] = x2;

     return 1;
  }

  const int nr = m.nrings();

  for (int i = 0; i < nr; ++i)
  {
    const Ring * r = m.ringi(i);
    if (r->contains(zatom))
      _list_of_ring_sizes.add_if_not_already_present(r->number_elements());
  }

  for (int i = 0; i < m.non_sssr_rings(); ++i)
  {
    const Ring * r = m.non_sssr_ring(i);

    if (r->contains(zatom))
      _list_of_ring_sizes.add_if_not_already_present(r->number_elements());
  }

  x1 -=  731 * _list_of_ring_sizes[0];
  x1 += 4222 * _list_of_ring_sizes[0];

  _invariant[1] = x1;
  _invariant[2] = x2;

  return 1;
}

int
Atom_Info::initialise(Molecule & m,
                      const atom_number_t zatom)
{
  _atom_number = zatom;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  int rbc = 0;
  uint64_t zsum = 0;

  const int bs = bond_and_atom_score(m, zatom, rbc, zsum);

  _invariant[0] = a->element()->atomic_symbol_hash_value() + 82828282 * zsum;

  const int arom = m.is_aromatic(zatom);

  uint64_t x1 = 10000000000 * acon + 288578282 * rbc                + 26012422 * a->formal_charge() - 1942 * bs;
  uint64_t x2 = 21928000000 * rbc  + 911971781 * a->formal_charge() + 10970988 * acon + 19204 * bs;

  x1 -= 9076 * arom;
  x2 += 58 * arom;

  if (0 == a->isotope())
    ;
  else if (include_isotopic_information_in_unique_smiles())
  {
    x1 += 571 * a->isotope();
    x2 += 82122 * a->isotope();
  }

  if (consider_implicit_hydrogens_in_unique_smiles)
  {
    const int ih = m.implicit_hydrogens(zatom);

    x1 += 1972 * ih;
    x2 -= 82221 * ih;
  
//  if (a->implicit_hydrogens_known())       bad idea, breaks things
//    x2 -= 5;
  }

  if (include_chiral_info_in_smiles())
  {
    const Chiral_Centre * c = m.chiral_centre_at_atom(zatom);
    if (NULL != c)
    {
      x1 -= 821 * bs;
      x2 += 741 * acon;
    }
  }

  if (0 == rbc)
  {
    _invariant[1] = x1;
    _invariant[2] = x2;

     return 1;
  }

  List_of_Ring_Sizes rs1, rs2;
  m.ring_sizes_for_atom(zatom, rs1);
  m.ring_sizes_for_non_sssr_rings(zatom, rs2);

  for (int i = 0; i < rs1.number_elements(); ++i)
  {
    _list_of_ring_sizes.add(rs1[i]);
  }

  x1 -= 731 * rs1[0];
  x2 += 4222 * rs1[0];

  for (int i = 0; i < rs2.number_elements(); ++i)
  {
    _list_of_ring_sizes.add(rs2[i]);
  }

  _invariant[1] = x1;
  _invariant[2] = x2;

  return 1;
}

//#define DEBUG_TARGET_ATOM_COMPARITOR

/*
   The first part of the ranking is to order the atoms. We use a molecule target to
   do that.
*/

static int
atom_info_comparitor(Atom_Info * const * ppa1, Atom_Info * const * ppa2)
{
  const Atom_Info *p1 = *ppa1;
  const Atom_Info *p2 = *ppa2;

  return p1->compare(*p2);
}

/*
  Use a molecule_to_Match object to order the atoms in the molecule.
  This is kind of sub-optimal in that we have to do extra
  comparisons after the sort to assign the ranks

  the SINGLE_ARRAY business is strange.

  Timing seems to vary with compiler, and usually individual allocations
  rather than a while array, are faster! Not sure why...
*/

void
Unique_Determination::_assign_initial_ranks_v2()
{
  assert (0 == _number_elements);

  this->resize(_matoms);

//#define SINGLE_ARRAY
#ifdef SINGLE_ARRAY
  Atom_Info * x = new Atom_Info[_matoms]; std::unique_ptr<Atom_Info[]> free_x;

  resizable_array<Atom_Info *> target;
#else
  resizable_array_p<Atom_Info> target;
#endif

  target.resize(_matoms);

  int molecule_contains_chirality = 0;
  _nchiral = 0;

  assert (_include_chiral_info_in_smiles == include_chiral_info_in_smiles());

  if (include_chiral_info_in_smiles() && (_nchiral = _m->chiral_centres()) > 0)
  {
    for (int i = 0; i < _matoms; i++)
    {
#ifdef SINGLE_ARRAY
      Atom_Info * a = x + i;
#else
      Atom_Info * a = new Atom_Info();
#endif

      const Chiral_Centre * c = _m->chiral_centre_at_atom(i);

      a->initialise(*_m, i, c);

      target.add(a);
    }
  }
  else
  {
    for (int i = 0; i < _matoms; i++)
    {
#ifdef SINGLE_ARRAY
      Atom_Info * a = x + i;
#else
      Atom_Info * a = new Atom_Info();
#endif

      a->initialise(*_m, i, NULL);

      target.add(a);
    }
  }

  molecule_contains_chirality = _nchiral;

  target.sort(atom_info_comparitor);

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After sorting targets\n";
  for (int i = 0; i < _matoms; i++)
  {
    const Atom_Info * t = target[i];
    const atom_number_t a = t->atom_number();
    cerr << "i = " << i << " is atom " << a << " (" << _m->smarts_equivalent_for_atom(a) << ")" << endl;
  }
#endif

  _initialise_rank_in_use();

  int rank_to_assign = 0;

  for (int i = 0; i < _matoms; i++)
  {
    int tmp;
    if (i > 0)    // only compare with previous if there is a previous
      tmp = atom_info_comparitor(&(target[i]), &(target[i - 1]));
    else
      tmp = 0;

    if (tmp < 0)
    {
      cerr << "Atoms out of order\n";
      target[i]->debug_print(cerr);
      target[i - 1]->debug_print(cerr);
    }
    assert (tmp >= 0);    // the array is supposed to be sorted

    if (tmp > 0)     // different from one before, increment rank
      rank_to_assign += _rank_delta;

    const atom_number_t a = target[i]->atom_number();

    Atom_and_Rank * r = new Atom_and_Rank(a, rank_to_assign, _m->ncon(a));

    if (molecule_contains_chirality)
    {
      r->set_chirality(_m->chiral_centre_at_atom(a));
      if (r->chiral_centre())
        molecule_contains_chirality--;
    }

    add(r);
    _atom_xref[a] = r;
    _rank_in_use[rank_to_assign]++;
  }

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  for (int i = 0; i < _matoms; i++)
  {
    if (NULL == _atom_xref[i])
    {
      cerr << "Yipes, the point for atom " << i << " is null\n";
    }
  }
#endif

// Now add connections to each atom

  for (int i = 0; i < _matoms; i++)
  {
    Atom_and_Rank * r = _things[i];
    atom_number_t a = r->atom_number();

    const Atom * ai = _m->atomi(a);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = ai->item(j);

      if (! include_directional_bonding_information_in_unique_smiles)
        ;
      else if (b->is_double_bond() && b->part_of_cis_trans_grouping())
        _cis_trans_bonds++;

      atom_number_t k = b->other(a);

      r->is_connected_to(_atom_xref[k], b);
    }

    _things[i]->establish_neighbours();
  }

  _nactive = _matoms;

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After establishing initial ranks\n";
  debug_print(cerr);
#endif

  return;
}

#endif

/*
  Use a molecule_to_Match object to order the atoms in the molecule.
  This is kind of sub-optimal in that we have to do extra
  comparisons after the sort to assign the ranks
*/

void
Unique_Determination::_assign_initial_ranks()
{
  assert (0 == _number_elements);

  resizable_array_p<Target_Atom> target;
  target.resize(_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    Target_Atom * a = new Target_Atom();

    a->initialise(_m, i, const_cast<Atom *>(_m->atomi(i)), NULL);

    target.add(a);
  }

  target.sort(target_atom_comparitor);

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After sorting targets\n";
  for (int i = 0; i < _matoms; i++)
  {
    const Target_Atom * t = target[i];
    atom_number_t a = t->atom_number();
    cerr << "i = " << i << " is atom " << a << " (" << _m->atomic_symbol(a) << ")" <<
            t->ncon() << " connections\n";
  }
#endif

  _initialise_rank_in_use();

  int rank_to_assign = 0;

  for (int i = 0; i < _matoms; i++)
  {
    int tmp;
    if (i > 0)    // only compare with previous if there is a previous
      tmp = target_atom_comparitor(&(target[i]), &(target[i - 1]));
    else
      tmp = 0;

    if (tmp < 0)
    {
      cerr << "Atoms out of order\n";
      target[i]->debug_print(cerr);
      target[i - 1]->debug_print(cerr);
    }
    assert (tmp >= 0);    // the array is supposed to be sorted

    if (tmp > 0)     // different from one before, increment rank
      rank_to_assign += _rank_delta;

    atom_number_t a = target[i]->atom_number();

    Atom_and_Rank * r = new Atom_and_Rank(a, rank_to_assign, target[i]->ncon());

    if (_include_chiral_info_in_smiles)
    {
      r->set_chirality(_m->chiral_centre_at_atom(a));
      if (r->chiral_centre())
        _nchiral++;
    }

    add(r);
    _atom_xref[a] = r;
    _rank_in_use[rank_to_assign]++;
  }

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  for (int i = 0; i < _matoms; i++)
  {
    if (NULL == _atom_xref[i])
    {
      cerr << "Yipes, the point for atom " << i << " is null\n";
    }
  }
#endif

// Now add connections to each atom

  for (int i = 0; i < _matoms; i++)
  {
    Atom_and_Rank * r = _things[i];
    atom_number_t a = r->atom_number();

    const Atom * ai = _m->atomi(a);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = ai->item(j);

      if (! include_directional_bonding_information_in_unique_smiles)
        ;
      else if (b->is_double_bond() && b->part_of_cis_trans_grouping())
        _cis_trans_bonds++;

      atom_number_t k = b->other(a);

      r->is_connected_to(_atom_xref[k], b);
    }

    _things[i]->establish_neighbours();
  }

  _nactive = _matoms;

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After establishing initial ranks\n";
  debug_print(cerr);
#endif

  return;
}

/*
  I have deliberately duplicated the code for assigning initial ranks with a subset
  in order to avoid slowing down the version that doesn't do subsetting
*/

void
Unique_Determination::_assign_initial_ranks (const int * include_atom)
{
  assert (0 == _number_elements);

  resizable_array_p<Target_Atom> target;
  target.resize(_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    Target_Atom * s = new Target_Atom();

    Atom * a = const_cast<Atom *>(_m->atomi(i));   // loss of const OK

    s->initialise(_m, i, a, NULL);

    if (! include_atom[i])
    {
      s->set_element(NULL);
      s->set_ncon(0);
    }
    else    // only find neighbours for molecules that are in the subset
    {
      int neighbours_included = 0;

      int acon = a->ncon();

      for (int j = 0; j < acon; j++)
      {
        atom_number_t k = a->other(i, j);

        if (include_atom[k])
          neighbours_included++;
      }

      if (neighbours_included != acon)
        s->set_ncon(neighbours_included);
    }

    target.add(s);
  }

  target.sort(target_atom_comparitor);

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After sorting targets\n";
  for (int i = 0; i < _matoms; i++)
  {
    const Target_Atom * t = target[i];
    atom_number_t a = t->atom_number();
    cerr << "i = " << i << " is atom " << a << " (" << _m->atomic_symbol(a) << ") " << t->ncon() << " connections" << " inc= " << include_atom[a] << '\n';
  }
#endif

  _initialise_rank_in_use();

//set_vector(_atom_xref, _m->natoms(), static_cast<Atom_and_Rank *>(NULL));

  int rank_to_assign = 0;

  for (int i = 0; i < _matoms; i++)
  {
    int tmp;
    if (i > 0)    // only compare with previous if there is a previous
      tmp = target_atom_comparitor(&(target[i]), &(target[i - 1]));
    else
      tmp = 0;

    if (tmp < 0)
    {
      cerr << "Atoms out of order\n";
      target[i]->debug_print(cerr);
      target[i - 1]->debug_print(cerr);
    }
    assert (tmp >= 0);    // the array is supposed to be sorted

//  if (tmp > 0 || INVALID_ATOMIC_NUMBER == target[i]->atomic_number())     // different from one before, increment rank
    if (tmp > 0 || 0 == include_atom[i])   // different from one before, increment rank
      rank_to_assign += _rank_delta;

    atom_number_t a = target[i]->atom_number();

    Atom_and_Rank * r = new Atom_and_Rank(a, rank_to_assign, target[i]->ncon());

    if (! _include_chiral_info_in_smiles)
      ;
    else if (target[i]->ncon() < 3)
      ;
    else
    {
      const Chiral_Centre * c = _m->chiral_centre_at_atom(a);
      if (NULL == c)
        ;
      else if (c->all_atoms_in_subset(include_atom, 1))
      {
        r->set_chirality(c);
        _nchiral++;
      }
    }

    add(r);
    _atom_xref[a] = r;
    _rank_in_use[rank_to_assign]++;
  }

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  for (int i = 0; i < _matoms; i++)
  {
    if (NULL == _atom_xref[i])
    {
      cerr << "Yipes, the point for atom " << i << " is null\n";
    }
  }
#endif

// Now add connections to each atom

  for (int i = 0; i < _matoms; i++)
  {
    Atom_and_Rank * r = _things[i];
    atom_number_t a = r->atom_number();

    if (! include_atom[a])
      continue;

    const Atom * ai = _m->atomi(a);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = ai->item(j);
      atom_number_t k = b->other(a);

      if (include_atom[k])
        r->is_connected_to(_atom_xref[k], b);
    }

    _things[i]->establish_neighbours();
  }

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After establishing neighbours\n";
  for (int i = 0; i < _matoms; i++)
  {
    const Atom_and_Rank * r = _things[i];

    atom_number_t a = r->atom_number();

    if (! include_atom[a])
      continue;

    cerr << "i = " << i << ' ';
    _things[i]->debug_print(cerr);
  }
#endif

  _nactive = _matoms;

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After establishing initial ranks\n";
  debug_print(cerr);
#endif

  return;
}

/*void
Unique_Determination::_adjust_initial_ranks_for_cis_trans_bonds ()
{
  return;
  int ne = _m->nedges();

#ifdef DEBUG_CIS_TRANS_STUFF
  cerr << "Checking " << ne << " bonds for cis-trans\n";
#endif

  for (int i = 0; i < ne; i++)
  {
    const Bond * b = _m->bondi(i);

    if (! b->part_of_cis_trans_grouping())
      continue;

    _adjust_initial_ranks_for_cis_trans_bonds(b->a1(), b->a2());
  }

  return;
}*/

int
Unique_Determination::_get_rank (atom_number_t a) const
{
  if (INVALID_ATOM_NUMBER == a)
    return -1;

  const Atom_and_Rank * ar = _atom_xref[a];

  return ar->rank();
}

int
Unique_Determination::_index_if_active (atom_number_t a) const
{
  for (int i = 0; i < _nactive; i++)
  {
    if (a == _things[i]->atom_number())
      return i;
  }

  return -1;
}

/*
   NW        NE
     \      /
      C == C
     /      \
   SW        SE

   In this function A is the LHS carbon atom
*/

unsigned int
Unique_Determination::_compute_cis_trans_rank (Atom_and_Rank * a, 
                       atom_number_t nw,
                       atom_number_t sw,
                       atom_number_t a2,
                       atom_number_t ne,
                       atom_number_t se)
{
  int rnw = _get_rank(nw);
  int rsw = _get_rank(sw);
  int ra2 = _get_rank(a2);
  int rne = _get_rank(ne);
  int rse = _get_rank(se);

//#define DEBUG_CIS_TRANS_STUFF
#ifdef DEBUG_CIS_TRANS_STUFF
  _m->debug_print(cerr);
  cerr << "A1 " << a->atom_number() << endl;
  cerr << "A2 " << a2 << endl;
  cerr << "NW " << nw << " rank " << rnw << endl;
  cerr << "SW " << sw << " rank " << rsw << endl;
  cerr << "NE " << ne << " rank " << rne << endl;
  cerr << "SE " << se << " rank " << rse << endl;
#endif

  unsigned int r = a->rank();

  r += ra2 * _m->natoms();
  if (rnw >= 0)
    r += 177 * rnw;
  if (rsw >= 0)
    r += 217 * rsw;
  if (rne >= 0)
    r += 303 * rne;
  if (rse >= 0)
    r += 419 * rse;

#ifdef DEBUG_CIS_TRANS_STUFF
  cerr << "Updated rank " << r << endl;
#endif

  return r;
}

/*
   NW   
     \   
      C ==
     / 
   SW 
*/

static int
identify_directionally_attached_atoms (const Atom * a,
                                       atom_number_t zatom,
                                       atom_number_t & nw,
                                       atom_number_t & sw)
{
  nw = INVALID_ATOMIC_NUMBER;
  sw = INVALID_ATOMIC_NUMBER;

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item (i);

    if (b->is_double_bond())
      continue;

    if (! b->is_directional())
      continue;

    if (b->is_directional_up())
    {
      if (b->a1() == zatom)
        nw = b->other(zatom);
      else
        sw = b->other(zatom);
    }
    else if (b->is_directional_down())
    {
      if (b->a1() == zatom)
        sw = b->other(zatom);
      else
        nw = b->other(zatom);
    }
  }

  if (nw == INVALID_ATOMIC_NUMBER && sw == INVALID_ATOMIC_NUMBER)
    return 0;

  return 1;
}


/*void
Unique_Determination::_adjust_initial_ranks_for_cis_trans_bonds (atom_number_t a1, atom_number_t a2)
{
  atom_number_t nw, sw;

  if (! identify_directionally_attached_atoms (_m->atomi(a1), a1, nw, sw))
    return;

  atom_number_t ne, se;

  if (! identify_directionally_attached_atoms (_m->atomi(a2), a2, ne, se))
    return;

  int ar1 = _index_of_atom(a1);
  int ar2 = _index_of_atom(a2);

  int newr1 = _compute_cis_trans_rank (_things[ar1], nw, sw, a2, ne, se);
  int newr2 = _compute_cis_trans_rank (_things[ar2], ne, se, a1, nw, sw);

  _things[ar1]->set_rank(newr1);
  _things[ar2]->set_rank(newr2);

  return;
}*/

int
Unique_Determination::_index_of_atom (atom_number_t zatom) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (zatom == _things[i]->atom_number())
      return i;
  }

  return -1;
}

int
Unique_Determination::canonical_order (Molecule & m,
                                       int * canonical_rank,
                                       const int * include_atom)
{
#ifdef DEBUG_UNIQUE_DETERMINATION
  cerr << "Begin canonical order computation\n";
#endif

  assert (NULL != canonical_rank);
  assert (m.ok());

  (void) _initialise(m, include_atom);

  if (_matoms < 2)
    return 1;

  int rc = _canonical_order(0);    // 0 means don't stop when symmetry is perceived

  copy_vector(canonical_rank, _canonical_rank, _matoms);

#ifdef DEBUG_UNIQUE_DETERMINATION
  cerr << "Canonical order is\n";
  for (int i = 0; i < _matoms; i++)
  {
    cerr << "Atom " << i << " (" << std::setw(2) << m.atomic_symbol(i) << " " <<
            m.ncon(i) << " connections, " << m.nbonds(i) << " bonds ) " << 
            canonical_rank[i] << " symmetry " << _symmetry[i];
    if (NULL != include_atom)
      cerr << " inc " << include_atom[i];
    cerr << '\n';
  }
#endif

  return rc;
}

int
Molecule::compute_canonical_ranking (Symmetry_Class_and_Canonical_Rank & sccr,
                                     const int * include_atom)
{
  assert (ok());

  if (0 == _number_elements)
    return 1;

  compute_aromaticity_if_needed();

  if (! sccr.arrays_allocated())
    sccr.allocate_arrays(_number_elements);

  Unique_Determination unqd;

  int rc = unqd.canonical_order(*this, sccr.canonical_rank(), include_atom);

  unqd.symmetry(_number_elements, sccr.symmetry_class());    // must do this after the canonical ranking

  return rc;
}

int
Molecule::compute_canonical_ranking ()
{
  return compute_canonical_ranking(_symmetry_class_and_canonical_rank, NULL);
}

int
Molecule::canonical_rank (atom_number_t a)
{
  assert(ok_atom_number(a));

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  return _symmetry_class_and_canonical_rank.canonical_rank(a);
}

int
Molecule::canonical_ranks (int * r)
{
  if (0 == _number_elements)
    return 0;

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  copy_vector(r, _symmetry_class_and_canonical_rank.canonical_rank(), _number_elements);

  return _number_elements;
}

const int *
Molecule::canonical_ranks ()
{
  if (0 == _number_elements)
    return NULL;

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  return _symmetry_class_and_canonical_rank.canonical_rank();
}

const int *
Molecule::symmetry_classes () 
{
  if (0 == _number_elements)
    return NULL;

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  assert (NULL != _symmetry_class_and_canonical_rank.symmetry_class());

  return _symmetry_class_and_canonical_rank.symmetry_class();
}

int
Molecule::symmetry_class (atom_number_t a)
{
  assert (ok_atom_number(a));

  (void) symmetry_classes();

  return _symmetry_class_and_canonical_rank.symmetry_class(a);
}

int
Molecule::symmetry_equivalents (atom_number_t a, Set_of_Atoms & sym)
{
  assert (ok_atom_number(a));

  const int * symmc = symmetry_classes();

  int s = symmc[a];

  sym.resize_keep_storage(0);

  for (int i = 0; i < _number_elements; i++)
  {
    if (s != symmc[i])
      continue;

    if (i == a)
      continue;

    sym.add(i);
  }

  return sym.number_elements();
}

int
Molecule::number_symmetry_classes ()
{
  const int * symmc = symmetry_classes();     // force symmetry perception

  const int array_size = _number_elements + 1;
  int * tmp = new_int(array_size); std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < _number_elements; i++)
  {
    tmp[symmc[i]]++;
  }

  int rc = 0;
  for (int i = 0; i < array_size; i++)
  {
    if (tmp[i])
      rc++;
  }

  return rc;
}

int
Molecule::bond_symmetry_class_small_memory (int * s)
{
  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();
 
  int nb = _bond_list.number_elements();

  const int * symm = _symmetry_class_and_canonical_rank.symmetry_class();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (symm[a1] < symm[a2])
      s[i] = symm[a2] * _number_elements + symm[a1];
    else
      s[i] = symm[a1] * _number_elements + symm[a2];

    s[i] = - s[i];    // convert to a negative number
  }

  int next_number_to_assign = -1;    // we increment it before using it

  for (int i = 0; i < nb; i++)
  {
    int change_from = s[i];

    if (change_from >= 0)    // already processed
      continue;

    next_number_to_assign++;

    s[i] = next_number_to_assign;
    for (int j = i + 1; j < nb; j++)
    {
      if (s[j] == change_from)
        s[j] = next_number_to_assign;
    }
  }

  return next_number_to_assign + 1;   // the number of different bond classes
}

int
Molecule::bond_symmetry_class_large_memory (int * s)
{
  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  int ns = number_symmetry_classes();

//cerr << "There are " << ns << " symmetry classes\n";

  const int * symm = _symmetry_class_and_canonical_rank.symmetry_class();

  int * tmp = new_int((ns + 1) * (ns + 1), -1); std::unique_ptr<int[]> free_tmp(tmp);

  int nb = _bond_list.number_elements();

  int next_number_to_assign = -1;     // incremented before use

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    int j;
    if (symm[a1] < symm[a2])
      j = symm[a2] * ns + symm[a1];
    else
      j = symm[a1] * ns + symm[a2];

    if (tmp[j] < 0)
    {
      next_number_to_assign++;
      tmp[j] = next_number_to_assign;
//    cerr << " j = " << j << " assigned " << tmp[j] << '\n';
    }

//  cerr << "Classes " << symm[a1] << " to " << symm[a2] << " j = " << j << " assigned " << tmp[j] << '\n';

    s[i] = tmp[j];
  }

  return next_number_to_assign + 1;
}

void
reset_unique_file_scope_variables()
{
  file_scope_include_isotopic_information_in_unique_smiles = 1;
  include_directional_bonding_information_in_unique_smiles = 1;
  consider_implicit_hydrogens_in_unique_smiles = 1;
  resolve_ties_by_geometry = 0;
  file_scope_consider_isotopes_as_zero_and_non_zero = 0;

  return;
}
