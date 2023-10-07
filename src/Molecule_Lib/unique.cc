#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION

#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/timsort.hpp"
#include "Foundational/iwqsort/iwqsort.h"

#include "chiral_centre.h"
#include "misc2.h"
#include "molecule.h"
#include "path.h"
#include "smiles.h"
#include "target.h"

using std::cerr;
using std::endl;

static int file_scope_include_isotopic_information_in_unique_smiles = 1;

void
set_include_isotopic_information_in_unique_smiles(int s)
{
  file_scope_include_isotopic_information_in_unique_smiles = s;

  return;
}

int
include_isotopic_information_in_unique_smiles() 
{
  return file_scope_include_isotopic_information_in_unique_smiles;
}

static int include_directional_bonding_information_in_unique_smiles = 1;

void
set_include_directional_bonding_information_in_unique_smiles(int s)
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
set_consider_isotopes_as_zero_and_non_zero(int s)
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
consider_isotopes_as_zero_and_non_zero()
{
  return file_scope_consider_isotopes_as_zero_and_non_zero;
}

// By default use faster canonicalization.
static int file_scope_use_legacy_initial_atom_ordering = 0;

void
set_unique_smiles_legacy_atom_ordering(int s) {
  file_scope_use_legacy_initial_atom_ordering = s;
}

// Very simplistic tool to keep track of which ranks have been assigned.
// Quite dangerous since we expose operator[] which means callers can basically
// do anything.
// When a caller requests an unused rank, only then do we zero a chunk of ranks
// for them.
class RankInUse {
  private:
    int * _rank_in_use;
    int _rank_delta;

    // The next unused rank to be given out.
    int _next_to_return;

    // Not clear that this is needed beyong the call to Initialise.
    // Retained for now, useful for checking.
    int _matoms;

  public:
    RankInUse();
    ~RankInUse();

    // The internal _rank_in_use array will be sized (matoms * expansion).
    // rank_delta is the offset between ranks given out by IdentifyUnusedRank.
    int Initialise(int matoms, int expansion, int rank_delta);

    int & operator[](int i) { return _rank_in_use[i];}
    int operator[](int i) const { return _rank_in_use[i];}

    // Used to restart rank assignment.
    void reset();

    // Return the index of an unused rank. The _rank_in_use array will be
    // zero'd from this newly issued rank up till the next rank that might
    // be issued.
    int IdentifyUnusedRank();
};

RankInUse::RankInUse() {
  _rank_in_use = nullptr;
  _rank_delta = 0;
  _next_to_return = 0;
  _matoms = 0;
}

RankInUse::~RankInUse() {
  if (_rank_in_use != nullptr) {
    delete [] _rank_in_use;
  }
}

int
RankInUse::Initialise(int matoms, int expansion, int rank_delta) {
  if (_rank_in_use != nullptr) {
    delete [] _rank_in_use;
  }
  _rank_in_use = new int[matoms * expansion];
  _next_to_return = 0;
  _rank_delta = rank_delta;
  _matoms = matoms;
  return 1;
}

void
RankInUse::reset() {
  _next_to_return = 0;
}

int RankInUse::IdentifyUnusedRank() {
#ifdef DEBUG_IDENTIFY_UNUSED_RANK
  cerr << "RankInUse::IdentifyUnusedRank:_next_to_return " << _next_to_return << " matoms " << _matoms << " delta " << _rank_delta << endl;
#endif
  std::fill_n(_rank_in_use + _next_to_return, _rank_delta, 0);
  int rc = _next_to_return;
  _next_to_return += _rank_delta;
  return rc;
}


std::tuple<int, int>
Ncon2Ahc(const Molecule& m, atom_number_t zatom) {
  const Atom * a = m.atomi(zatom);
  int ncon2 = 0;
  int attached_heteroatom_count = 0;
  for (const Bond * b : *a) {
    atom_number_t o = b->other(zatom);
    const Atom * other_atom = m.atomi(o);
    ncon2 += other_atom->ncon();
    if (other_atom->atomic_number() == 6)
      ;
    else if (other_atom->atomic_number() == 1)
      ;
    else
      attached_heteroatom_count++;
  }

  return {ncon2, attached_heteroatom_count};
}

// Return the tuple of properties, derived from the neighbors of `zatom`.
// 1. number of connections of connnected to the neighbors of zatom.
// 2. the sum of the atomic numbers of those neighbors.
// 3. the sum of the ring bond counts of the neighbors.
std::tuple<int, int, int>
SingleShellProperties(Molecule& m, atom_number_t zatom) {
  const Atom * a = m.atomi(zatom);
  int ncon2 = 0;
  int sum_atomic_numbers = 0;
  int rbc = 0;
  for (const Bond * b : *a) {
    atom_number_t o = b->other(zatom);
    const Atom * other_atom = m.atomi(o);
    ncon2 += other_atom->ncon();
    sum_atomic_numbers += other_atom->atomic_number();
    rbc += m.ring_bond_count(o);
  }

  return {ncon2, sum_atomic_numbers, rbc};
}

// 64 bits of data holding atomic properties.
// Note that this only works for elements that are in the periodic table,
// and for molecules that are in general like typical organic molecules.
struct AtomProperties {
  // Only if part the periodic table.
  int16_t _atomic_symbol_hash;
  // 4 bits each for ncon and nbonds
  uint8_t _ncon_nbonds;
  // Formal charge.
  int8_t _fc;
  // Bits set for is aromatic, is chiral, and 6 bits for implicit Hydrogens.
  uint8_t _arom_chiral_imph;
  // Ncon2
  uint8_t _ncon2;
  // ring bond count.
  uint8_t _ring_bond_count;
  // Sum of neighboring atomic numbers.
  uint8_t _attached_heteroatom_count;
};

// Ring information is stored in words, to facilitate quick OR operations
// to tell if a given ring size is in another list of ring sizes;
// 16 membered rings are the max that can be stored.
struct RingInfo {
  // The smallest ring containing the atom.
  uint16_t _smallest_ring_size;
  // A bit is set for each ring size.
  uint16_t _ring_sizes;
};

// A class to hold atomic properties used in uniqueness determinations.
// All properties are computed at creation.
// Values are packed into 64 bit ints that can be compared as one.
class AtomPropertiesForRanking {
  private:
    // This will be constructed as an AtomProperties and which gets copied.
    uint64_t _hash;
    uint32_t _ring_size_hash;
    atom_number_t _atom_number;

  public:
    // A constructor for atoms not being processed.
    AtomPropertiesForRanking();  
    // A constructor for atoms being ranked.
    AtomPropertiesForRanking(Molecule& m, atom_number_t zatom, const Chiral_Centre* chiral);
    void initialise(Molecule& m, atom_number_t zatom, const Chiral_Centre* chiral);

    int debug_print(std::ostream& output) const;

    uint64_t hash() const { return _hash;}
    uint32_t ring_size_hash() const { return _ring_size_hash;}
    atom_number_t atom_number() const { return _atom_number;}
};

AtomPropertiesForRanking::AtomPropertiesForRanking() {
  _hash = 0;
  _ring_size_hash = 0;
  _atom_number = INVALID_ATOM_NUMBER;
}

AtomPropertiesForRanking::AtomPropertiesForRanking(Molecule & m, atom_number_t zatom,
                        const Chiral_Centre* chiral)
{
  initialise(m, zatom, chiral);
}

void
AtomPropertiesForRanking::initialise(Molecule& m,
                        atom_number_t zatom,
                        const Chiral_Centre* chiral)
{
  _hash = 0;
  _ring_size_hash = 0;
  _atom_number = zatom;

  AtomProperties * aprop = reinterpret_cast<AtomProperties*>(&_hash);

  Atom * a = const_cast<Atom*>(m.atomi(zatom));  // loss of const OK.

  aprop->_atomic_symbol_hash = a->element()->atomic_symbol_hash_value();

  aprop->_ncon_nbonds = a->ncon() << 4;
  aprop->_ncon_nbonds |= a->nbonds();

//cerr << "ncon " << a->ncon() << " nbonds " << a->nbonds() << " result " << int(aprop->_ncon_nbonds) << " type " << m.smarts_equivalent_for_atom(zatom) << '\n';

  aprop->_fc = a->formal_charge();

  const int rbc = m.ring_bond_count(zatom);
  int aromatic = 0;
  if (rbc > 0 && m.is_aromatic(zatom))
    aromatic = 1;

  int implicit_h = 0;
  if (consider_implicit_hydrogens_in_unique_smiles)
    implicit_h = a->implicit_hydrogens();

  // _ring_arom_imph a uint8_t type.
  aprop->_arom_chiral_imph = implicit_h;
  if (aromatic) {
    aprop->_arom_chiral_imph |= (1 << 7);
  }
  if (chiral != nullptr) {
    aprop->_arom_chiral_imph |= (1 << 6);
  }

  auto [ncon2, sum_atomic_numbers, neighbour_rbc] = SingleShellProperties(m, zatom);
  aprop->_attached_heteroatom_count = sum_atomic_numbers;
  aprop->_ncon2 = ncon2;
  aprop->_ring_bond_count = 20 * rbc + neighbour_rbc;

  if (rbc == 0) {
    return;
  }

  RingInfo* ring_info = reinterpret_cast<RingInfo*>(&_ring_size_hash);

  int nr = m.nrings();
  for (int i = 0; i < nr; ++i) {
    const Ring * r = m.ringi(i);
    if (! r->contains(zatom)) {
      continue;
    }
    if (ring_info->_smallest_ring_size == 0) {
      ring_info->_smallest_ring_size = 1 << r->size();
    }
    ring_info->_ring_sizes |= (1 << r->size());
  }
  nr = m.non_sssr_rings();
  for (int i = 0; i < nr; ++i) {
    const Ring * r = m.non_sssr_ring(i);
    if (r->contains(zatom)) {
      ring_info->_ring_sizes |= (1 << r->size());
    }
  }

}

int
AtomPropertiesForRanking::debug_print(std::ostream & output) const {
  output << "AtomPropertiesComparitor for atom " << _atom_number << '\n';
  return output.good();
}

int
AtomPropertiesComparitor(AtomPropertiesForRanking* const * p1,
                         AtomPropertiesForRanking* const * p2) {
  if ((*p1)->hash() == (*p2)->hash()) {
    ;
  } else if ((*p1)->hash() < (*p2)->hash()) {
    return 1;
  } else {
    return -1;
  }

  // Atomic properties do not resolve, what about ring sizes.
  const auto rs1 = (*p1)->ring_size_hash();
  const auto rs2 = (*p2)->ring_size_hash();

  if (rs1 == rs2) {
    return 0;
  }
  const RingInfo* ri1 = reinterpret_cast<const RingInfo*>(&rs1);
  const RingInfo* ri2 = reinterpret_cast<const RingInfo*>(&rs2);

  const auto sr1 = ri1->_smallest_ring_size;
  const auto sr2 = ri2->_smallest_ring_size;

  if (sr1 < sr2)
  {
    if ((sr1 & ri2->_ring_sizes) == 0)
      return -1;

    if ((sr2 & ri1->_ring_sizes) == 0)
      return 1;
  }
  else if (sr1 > sr2)
  {
    if ((sr2 & ri1->_ring_sizes) == 0)
      return 1;

    if ((sr1 & ri2->_ring_sizes) == 0)
      return -1;
  }

  return 0;   // Should never come here.
}

// If we are going to use the AtomPropertiesForRanking struct to do
// uniqueness determinations, then we must impose some constraints on
// the molecule being processed, since that struct packs properties
// into narrow fields. FOr most reasonable organic structures, it
// should be fine.
bool
Molecule::_ok_for_fast_atom_comparisons()  {
  for (int i = 0; i < _number_elements; ++i) {
    const Atom * a = _things[i];
    if (a->ncon() > 16)
      return false;
    if (! a->element()->is_in_periodic_table())
      return false;
    if (a->formal_charge() == 0)
      ;
    else if (a->formal_charge() < -255)
      return false;
    else if (a->formal_charge() > 255)
      return false;
    if (a->isotope()) {
      return false;
    }
  }

  const int nr = nrings();
  if (nr == 0) {
    return true;
  }

  // Just check the largest ring for compatibility.
  if (ringi(nr - 1)->size() > 15)
    return false;

  return true;
}



#ifdef ONLY_USED_WITH_SOME_COMPILATION_OPTIONS
static void
get_bonds(const Atom * a,
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
#endif

#define SEEMS_FASTER_WITH_INTEL
#ifdef SEEMS_FASTER_WITH_INTEL
static int
set_bond_score_in_aryl(const Atom * a,
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
  
  p->set_aryl(rc);

  return rc;
}
#endif

//#define SEEMS_FASTER_WITH_GCC
#ifdef SEEMS_FASTER_WITH_GCC
static int
set_bond_score_in_aryl (const Atom * a,
                         Target_Atom * p)
{
  int aromatic_bonds, single_bonds, double_bonds, triple_bonds;

  get_bonds(a, aromatic_bonds, single_bonds, double_bonds, triple_bonds);

  int rc = 1000 * aromatic_bonds + 100 * single_bonds + 10 * double_bonds + triple_bonds;

  p->set_aryl(rc);

  return rc;
}
#endif

//#define DEBUG_TARGET_ATOM_COMPARITOR

/*
   The first part of the ranking is to order the atoms. We use a molecule target to
   do that.
*/

static int
target_atom_comparitor(Target_Atom * const * ppa1, Target_Atom * const * ppa2)
{
  Target_Atom *p1 = *ppa1;
  Target_Atom *p2 = *ppa2;

#ifdef DEBUG_TARGET_ATOM_COMPARITOR
  assert(p1->ok());      // too expensive normally
  assert(p2->ok());
  atom_number_t t1 = p1->atom_number();
  atom_number_t t2 = p2->atom_number();
  cerr << "Comparing atom " << t1 << " (" << p1->m()->atomic_symbol(t1) << ", " << p1->ncon() << " connections) and " <<
                               t2 << " (" << p2->m()->atomic_symbol(t2) << ", " << p1->ncon() << " connections)\n";
//cerr << "Atomic numbers " << p1->m()->atomic_number(t1) << " and " << p2->m()->atomic_number(t2) << "\n";
#endif

  const Element * e1 = p1->element();
  const Element * e2 = p2->element();

// Lots of different possibilities here because elements may be nullptr (excluded
// from the computation) or not from the periodic table (atomic number not valid)

  if (e1 == e2)    // covers case of both NULL
    ;
  else if (nullptr == e1)
    return -1;
  else if (nullptr == e2)
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

  if (kInvalidAtomicNumber == p1->atomic_number())    // these are atoms that have been excluded from consideration, they are always considered equivalent
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

// We store the bond scores in the target_atom's aryl attribute

  int nb1;

  if (! p1->aryl_value_set())
    nb1 = set_bond_score_in_aryl(a1, p1);
  else
    nb1 = p1->aryl();
  
  int nb2;  
  if (! p2->aryl_value_set())
    nb2 = set_bond_score_in_aryl(a2, p2);
  else
    nb2 = p2->aryl();

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
  cerr << "chirlaity ? " << include_chiral_info_in_smiles() << " " << (nullptr != p1->chiral_centre()) << " and " << (nullptr != p2->chiral_centre()) << '\n';
#endif

  if (! include_chiral_info_in_smiles())
    ;
  else if (nullptr == p1->chiral_centre() && nullptr == p2->chiral_centre())
    ;
  else if (nullptr != p1->chiral_centre() && nullptr != p2->chiral_centre())
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

    int ok() const;
    int debug_print(std::ostream &) const;

    int terse_print_neighbour_ranks(std::ostream &) const;

    int bt(int i) const { return _bt[i];}

    const Atom * atom() const { return _atom;}

    void set_chirality(const Chiral_Centre * c) { _chiral_centre = c;}
    const Chiral_Centre * chiral_centre() const { return _chiral_centre;}

    void set_consider_chirality(int s) { _considering_chirality = s;}

//  This function establishes the connectivity

    int is_connected_to(Atom_and_Rank *, const Bond *);

//  After we know about our connections, we can establish bond types

    void establish_neighbours();

    int compare(const Atom_and_Rank *) const;
    int compare_for_timsort(const Atom_and_Rank *) const;
    int compare_less(const Atom_and_Rank *) const;

//  called whenever a neighboring atom has been classified, and is
//  removed from the computation

    int a_neighbour_has_been_classified(const Atom_and_Rank * c);

    unsigned int  rank() const { return _rank;}

    void set_rank(int r) { _rank = r;}

    int  resolved_pending_processing() const { return _resolved_pending_processing;}
    void set_resolved_pending_processing() { _resolved_pending_processing = 1;}

    int  set_rank(int new_rank, RankInUse& rank_in_use);
    int  choose_an_unused_rank(RankInUse& rank_in_use);

    void set_atom_number(atom_number_t na) { _a = na;}
    atom_number_t atom_number() const { return _a;}

    int  ncon() const { return _connected_atoms.number_elements();}

    void collect_neighbour_ranks();

    int  is_connected_to(const Atom_and_Rank * z) const { return _connected_atoms.contains(z);}

//  When including chirality in a smiles, we need a means of getting a number
//  that reflects the chiral arrangement around an atom

    int  chiral_score(const unsigned int *) const;

    int  chirality_score() const { return _chirality_score;};
    void  compute_chirality_score(const unsigned int *);

    void propagate_chirality_influence(const unsigned int *) const;

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
Atom_and_Rank::_default_values()
{
  _ranks_of_neighbours = nullptr;
  _sum_of_neighbour_ranks = 0;

  _chiral_centre = nullptr;
  _chirality_score = 0;

  _resolved_pending_processing = 0;

  _ncon = 0;

  _chiral_neighbours = 0;

  _considering_chirality = 0;

  return;
}

Atom_and_Rank::Atom_and_Rank(atom_number_t za, int zr, int nc)
{
  _default_values();

  resize(nc);
  _bt.resize(nc);
  _connected_atoms.resize(nc);

  _a = za;
  _rank = zr;

  return;
}

Atom_and_Rank::~Atom_and_Rank()
{
  delete [] _ranks_of_neighbours;
  _ranks_of_neighbours = nullptr;

  return;
}

/*
  Scan the remaining connections and return the index of a given atom number
*/

int
Atom_and_Rank::_index_of_atom(atom_number_t a) const
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
Atom_and_Rank::set_rank(int new_rank, RankInUse& rank_in_use)
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
Atom_and_Rank::choose_an_unused_rank(RankInUse& rank_in_use)
{
  if (rank_in_use[_rank] == 1)    // it is already unique
    return 1;

  int new_rank = rank_in_use.IdentifyUnusedRank();

  rank_in_use[_rank]--;
  _rank = new_rank;
  rank_in_use[_rank] = 1;

  return _rank;
}

/*
  Notify one atom that another atom is bonded to it.
*/

int
Atom_and_Rank::is_connected_to(Atom_and_Rank * n,
                               const Bond * b)
{
  assert (_number_elements == _bt.number_elements());

  add(n);

  if (b->is_aromatic())
    _bt.add(AR_AROMATIC_BOND);
  else if (b->is_single_bond())
  {
    if (include_directional_bonding_information_in_unique_smiles && b->is_directional()) {
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
  else if (b->btype() == UNKNOWN_BOND_TYPE) {
    _bt.add(UNKNOWN_BOND_TYPE);
  } else {
    cerr << "What kind of bond is this!!! " << (*b) << '\n';
    debug_print(cerr);
    iwabort();
  }

  _connected_atoms.add(n);
  _ncon = _connected_atoms.number_elements();

  if (! include_chiral_info_in_smiles())
    ;
  else if (nullptr != n->chiral_centre())
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
Atom_and_Rank::establish_neighbours()
{
  assert (nullptr == _ranks_of_neighbours);

  if (_number_elements)
    _ranks_of_neighbours = new int[_number_elements];
//  _ranks_of_neighbours = new_int(_number_elements);    // not necessary for it to be 0'd, aesthetics only

  return;
}

/*
  Single step in the Morgan-like algorithm
*/
//#define DEBUG_COLLECT_NEIGHBOUR_RANKS

void
Atom_and_Rank::collect_neighbour_ranks()
{
  if (0 == _number_elements)
    return;

#ifdef DEBUG_COLLECT_NEIGHBOUR_RANKS
  cerr << "Atom_and_Rank::collect_neighbour_ranks: ranks";
  for (int i = 0; i < _number_elements; ++i) {
    cerr << ' ' << _things[i]->rank();
  }
  cerr << '\n';
#endif

  _sum_of_neighbour_ranks = 0;
  for (int i = 0; i < _number_elements; ++i) {
    _ranks_of_neighbours[i] = _bt[i] + _things[i]->rank();
    _sum_of_neighbour_ranks += _ranks_of_neighbours[i];
  }

  if (_number_elements == 1)
    return;

  if (_ranks_of_neighbours[0] > _ranks_of_neighbours[1])
    std::swap(_ranks_of_neighbours[0], _ranks_of_neighbours[1]);

  if (_number_elements == 2) {
    return;
  }

  // Some cheap pseudo sorting for the common cases.

  if (_ranks_of_neighbours[1] > _ranks_of_neighbours[2])
    std::swap(_ranks_of_neighbours[1], _ranks_of_neighbours[2]);
  if (_ranks_of_neighbours[0] > _ranks_of_neighbours[1])
    std::swap(_ranks_of_neighbours[0], _ranks_of_neighbours[1]);

  if (3 == _number_elements)
    return;

  // 4 is the last common case.
  if (_ranks_of_neighbours[2] > _ranks_of_neighbours[3])
    std::swap(_ranks_of_neighbours[2], _ranks_of_neighbours[3]);
  if (_ranks_of_neighbours[1] > _ranks_of_neighbours[2])
    std::swap(_ranks_of_neighbours[1], _ranks_of_neighbours[2]);
  if (_ranks_of_neighbours[0] > _ranks_of_neighbours[1])
    std::swap(_ranks_of_neighbours[0], _ranks_of_neighbours[1]);

#ifdef REMOVE_THIS
  for (int i = 0; i < 4; ++i) {
    cerr << "Check sort " << i << " _ranks_of_neighbours " << _ranks_of_neighbours[i] << endl;
    if (i > 0 && _ranks_of_neighbours[i] < _ranks_of_neighbours[i-1]) {
      cerr << "OUT OF ORDER\n";
    }
  }
#endif
  if (_number_elements == 4) {
    return;
  }

  // This should hardly ever happen.
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
Atom_and_Rank::ok() const
{
  if (! resizable_array<Atom_and_Rank *>::ok())
    return 0;

  if (_ncon < _number_elements)
    return 0;

  return 1;
}

int
Atom_and_Rank::debug_print(std::ostream & os) const
{
  os << "Atom " << _a << ", " << _number_elements << " neighbours, rank " << _rank << '\n';
  if (_number_elements < _ncon)
    os << _ncon << " connections within the molecule\n";

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom_and_Rank * n = _things[i];
    os << "Neighbour " << i << " is atom " << n->atom_number() << " btype " << _bt[i] << " rank " << n->rank() << endl;
  }

  if (_sum_of_neighbour_ranks)
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
Atom_and_Rank::terse_print_neighbour_ranks(std::ostream & os) const
{
  if (_sum_of_neighbour_ranks)
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
atom_and_rank_comparitor(Atom_and_Rank * const * pai1,
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
Atom_and_Rank::compare(const Atom_and_Rank * r2) const
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

  if (_sum_of_neighbour_ranks < r2->_sum_of_neighbour_ranks)
    return -1;
  if (_sum_of_neighbour_ranks > r2->_sum_of_neighbour_ranks)
    return 1;

  if (0 == _sum_of_neighbour_ranks)    // part of a subset, no neighbours computed
    return 0;

// Element by element comparison.

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
Atom_and_Rank::compare_for_timsort(const Atom_and_Rank * r2) const
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
Atom_and_Rank::_compare_by_chirality(const Atom_and_Rank * r2) const
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
Atom_and_Rank::_compare_by_chirality_for_timsort(const Atom_and_Rank * r2) const
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
Atom_and_Rank::chiral_score(const unsigned int * rank) const
{
  assert (nullptr != _chiral_centre);

  return _chiral_centre->orientation(rank);
}

/*
  We need to know the rank of an atom that is part of a chiral centre.
  Beware of lone pairs and implicit hydrogens
*/

static unsigned int
get_rank_for_connection(const unsigned int * rank, atom_number_t a)
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
Atom_and_Rank::_identify_two_unresolved_connections(const unsigned int * rank,
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
Atom_and_Rank::compute_chirality_score(const unsigned int * rank)
{
  if (nullptr == _chiral_centre)
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

    if (nullptr == c)
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

    RankInUse _rank_in_use;

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

    // If we are processing a molecule that does not fit into the compact
    // hash form, we use the legacy ordering.
    const bool _legacy_atom_ordering;

//  private functions;
  
    int _allocate_atom_arrays (int matoms, int rank_delta);
    int _free_all_arrays ();

    int _initialise (Molecule &, const int *);

    void _assign_initial_ranks();
    void _assign_initial_ranks_legacy();
    void _assign_initial_ranks_legacy (const int * include_atom);

    void _reassign_ranks ();
    void _expand_around_cis_trans_bonds();

//  Internal function to perform the analysis up to the point of full symmetry perception

    int _perceive_symmetry ();

//  Have the ranks changed from one iteration to the next

    int _ranks_changed();

//  Does a comparison of stored ranks for the extended neighbour list of
//  a set of atoms.

    int _rank_values_resolved(resizable_array<atom_number_t> &, int);

//  Called only from _single_step_process_unique_atoms. 

    void _atom_is_unique(atom_number_t a);

//  when computing chirality influences, we need the _rank array filled

    void _fill_rank_array_for_chirality();

//  Once we have perceived symmetry, we start to include the influence of chirality.

    void _turn_on_chirality_considerations_and_reassign_ranks();

//  Scan through the atom_and_rank array and identify atoms with unique ranks

    int _get_indices_of_unique_atoms(resizable_array<int> & unique_atoms) const;

//  These functions are used when unique atoms are identified in the
//  sorted array of ranks

    int _single_step_process_unique_atoms();

//  If atoms become disconnected, they are automatically done

    int _process_all_now_disconnected_atoms();

//  Once the ranking remains unchanged from one iteration to the next,
//  we store the resulting symmetry information

    int _store_symmetry_info();

//  Use this to break ties

    int _break_a_tie();

//  If all else fails, we can break ties by geometry

    int _choose_tie_breaker_by_geometry() const;

    template <typename T> int _identify_extreme_value(T & c) const;

    int _choose_tie_breaker_atom();

//  Fill the _chirality_score array

    void _compute_chirality_scores();

//  Look at a sequence of equivalent atoms and analyse their chirality scores

    void _analyse_chirality_in_sequence(int sstart, int & next_starting_position, int * cs, int & chiral_atoms_in_sequence) const;
    void _identify_next_sequence(int sstart, int & next_starting_position, int & chiral_atoms_in_sequence) const;

//  An atom has been classified as unique. Move it to the inactive part of the array

    int _move_to_inactive(int);

    int _process_all_unique_atoms();

    int _identify_unused_rank() const;
    int _identify_two_unused_ranks(int & r1, int & r2) const;
    int _identify_some_unused_ranks(int ranks_needed,
                      resizable_array<int> & ranks_identified) const;

//  When an atom is ranked, we change the rank values for those atoms connected

    int __adjust_rank_of_atoms_attached_to(Atom_and_Rank * r);
    int _adjust_rank_of_atoms_attached_to(Atom_and_Rank * r);

//  A single step of the expansion

    int _expand(int);

    int _canonical_order(int);

    int _index_of_atom(atom_number_t zatom) const;

    int _get_rank(atom_number_t a) const;
    int _index_if_active(atom_number_t a) const;
    void _adjust_initial_ranks_for_cis_trans_bonds();
    void _adjust_initial_ranks_for_cis_trans_bonds(atom_number_t a1, atom_number_t a2);
    unsigned int _compute_cis_trans_rank(Atom_and_Rank * a, 
                       atom_number_t nw,
                       atom_number_t sw,
                       atom_number_t a2,
                       atom_number_t ne,
                       atom_number_t se);

    int _expand_around_cis_trans_bond(atom_number_t a1,
                                atom_number_t a2,
                                unsigned int * new_rank);

    int _identify_directionally_attached_atoms(const Atom_and_Rank * ar,
                                atom_number_t & nw,
                                atom_number_t & sw) const;
    int _identify_directionally_attached_bonds(const Atom_and_Rank * ar,
                                                const Bond * & b1,
                                                const Bond * & b2) const;
    // Utility functions used in debugging.
    int _print_canonical_order(const Molecule & m,
                        const int * include_atom,
                        std::ostream& output) const;
    int _print_canonical_order_by_canonical_order(const Molecule & m,
                        const int * include_atom,
                        std::ostream& output) const;
  public:
    Unique_Determination(bool = false);
    ~Unique_Determination();

    int ok() const;
    int debug_print(std::ostream &) const;

    int symmetry(int, int *) const;    // return previously computed values
    int symmetry(Molecule *, int *, const int *);

    int canonical_order(Molecule &, int *, const int *);
};

//#define DEBUG_UNIQUE_DETERMINATION

Unique_Determination::Unique_Determination(bool legacy_atom_ordering) : _legacy_atom_ordering(legacy_atom_ordering)
{
  _matoms = 0;
  _nactive = 0;

  _rank = nullptr;

  _symmetry = nullptr;
  _canonical_rank = nullptr;

  _old_rank = nullptr;

  _atom_xref = nullptr;

  _use_chirality = 0;

  _include_chiral_info_in_smiles = include_chiral_info_in_smiles();

  std::random_device rd;

  _rng.seed(3172776704);    // just used for shuffling, so proper randomisation not important

  // Will be adjusted upwards if chirality is included.
  return;
}

/*
  To cut down on the number of allocations, we allocate one large array and have
  the sub-arrays be part of that.
*/

int
Unique_Determination::_allocate_atom_arrays(int matoms, int rank_delta)
{
  assert (matoms > 0);

  _matoms = matoms;

  assert (sizeof(int) == sizeof(unsigned int));

  _rank = new unsigned int[_matoms + _matoms + _matoms + _matoms];

  _canonical_rank = reinterpret_cast<int *>(_rank + _matoms);
  _symmetry       = reinterpret_cast<int *>(_rank + _matoms + _matoms);
  _old_rank       =                         _rank + _matoms + _matoms + _matoms;
  
  std::fill_n(_canonical_rank, _matoms, -1);
  std::fill_n(_symmetry, _matoms, -1);

  _atom_xref = new Atom_and_Rank *[_matoms];
  
  _rank_in_use.Initialise(matoms, 11, AR_TRIPLE_BOND + 1);

  return 1;
}

Unique_Determination::~Unique_Determination()
{
  _free_all_arrays();

  return;
}

int
Unique_Determination::_free_all_arrays()
{
  if (nullptr != _rank)
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
Unique_Determination::ok() const
{
  if (nullptr == _m)
    return 0;

  if (! _m->ok())
    return 0;

  return 1;
}

int
Unique_Determination::debug_print(std::ostream & os) const
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
Unique_Determination::_compute_chirality_scores()
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
Unique_Determination::_identify_next_sequence(int sstart,
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
Unique_Determination::_choose_tie_breaker_atom()
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
Unique_Determination::_break_a_tie()
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
Unique_Determination::__adjust_rank_of_atoms_attached_to(Atom_and_Rank * r)
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
    n->choose_an_unused_rank(_rank_in_use);

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
        n0->choose_an_unused_rank(_rank_in_use);

      n1->set_rank(n0->rank(), _rank_in_use);
    }
    else
    {
      int nr0 = _rank_in_use.IdentifyUnusedRank();
      int nr1 = _rank_in_use.IdentifyUnusedRank();

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

    const int new_rank = _rank_in_use.IdentifyUnusedRank();

    for (int i = 0; i < neighbours_of_r; i++)
    {
      Atom_and_Rank * n = r->item(i);
      n->set_rank(new_rank, _rank_in_use);
    }

    return 1;
  }

// All the others are done the hard way

  for (int i = 0; i < nr; i++)
  {
    unsigned int orank = old_ranks[i];
    int new_rank = _rank_in_use.IdentifyUnusedRank();
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
Unique_Determination::_adjust_rank_of_atoms_attached_to(Atom_and_Rank * r)
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
Unique_Determination::_atom_is_unique(atom_number_t a)
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
Unique_Determination::_get_indices_of_unique_atoms(resizable_array<int> & unique_atoms) const
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
Unique_Determination::_move_to_inactive(int a)
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
Unique_Determination::_fill_rank_array_for_chirality()
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
Unique_Determination::_single_step_process_unique_atoms()
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
Unique_Determination::_reassign_ranks()
{
#ifdef DEBUG_REASSIGN_RANKS
  cerr << "Reassigning ranks. Current ranks are\n";
  debug_print(cerr);
#endif

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

  _rank_in_use.reset();
  int rank_to_assign = _rank_in_use.IdentifyUnusedRank();

  for (int i = 0; i < _nactive; i++)
  {
    Atom_and_Rank * t = _things[i];
    if (t->rank() != rprev)
    {
      rank_to_assign = _rank_in_use.IdentifyUnusedRank();
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
Unique_Determination::_process_all_now_disconnected_atoms()
{
  resizable_array<int> unique_atoms;     // actually indices in the array
  unique_atoms.resize_keep_storage(_nactive);

  for (int i = 0; i < _nactive; i++)
  {
    Atom_and_Rank * r = _things[i];
    if (r->empty())
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

  unsigned int prev_rank = 0;   // Initialise to keep the compiler quiet.
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
Unique_Determination::_process_all_unique_atoms()
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
Unique_Determination::_ranks_changed()
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
atom_and_neighbour_rank_comparitor(Atom_and_Rank * const * pai1,
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
Unique_Determination::_expand(int collect_neighbours)
{
#ifdef DEBUG_EXPAND
  cerr << "At start of expand, collect_neighbours " << collect_neighbours << '\n';
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

// and assign ranks according to the neighbours rank. The only trickery here
// is that we must assign the new rank for the previous atom after we do the
// comparison, because atoms are ordered primarily by their current rank value

  _rank_in_use.reset();
  int next_rank_to_assign = _rank_in_use.IdentifyUnusedRank();

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
      next_rank_to_assign = _rank_in_use.IdentifyUnusedRank();
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
Unique_Determination::_initialise(Molecule & m,
                                  const int * include_atom)
{
  _m = &m;
  _matoms = m.natoms();

// Handle the two special cases specifically

  if (0 == _matoms)
    return 1;

  if (include_directional_bonding_information_in_unique_smiles && m.cis_trans_bonds_present()) {
    (void) _allocate_atom_arrays(_matoms, _matoms * 25);   
  } else {
    (void) _allocate_atom_arrays(_matoms, _matoms * 9);   
  }

  if (1 == _matoms) {
    _canonical_rank[0] = 1;
    _symmetry[0] = 1;
    return 1;
  }

  std::fill_n(_old_rank, _matoms, 0);

  _next_canonical_rank_to_assign = _matoms;
  _next_symmetry_class_to_assign = 1;
  _symmetry_stored = 0;

  _nchiral = 0;

  _cis_trans_bonds = 0;

  if (include_atom != nullptr) {
    if (_legacy_atom_ordering) {
      _assign_initial_ranks_legacy(include_atom);
    } else {
      // Would be a bunch of work to implement this, not sure if it is needed.
      cerr << "Unimplemented feature, see Ian\n";
      //_assign_initial_ranks(include_atom);
    }
  } else if (_legacy_atom_ordering) {
    _assign_initial_ranks_legacy();
  } else {
    _assign_initial_ranks();
  }

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  for (int i = 0; i < _matoms; i++)
  {
    _atom_xref[i] = nullptr;
  }
#endif

  return 1;
}

//#define DEBUG_SYMMETRY

int
Unique_Determination::_perceive_symmetry()
{
  return 1;
}

int
Unique_Determination::symmetry(Molecule * m, int * symmetry_class,
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
Unique_Determination::_canonical_order(int stop_when_symmetry_perceived)
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
#ifdef CIS_TRANS_DOES_NOT_WORK
      if (_cis_trans_bonds > 0)
      {
        _expand_around_cis_trans_bonds();
        continue;
      }
#endif

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
#ifdef CIS_TRANS_DOES_NOT_WORK
    if (_cis_trans_bonds > 0)
      _expand_around_cis_trans_bonds();
#endif

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

#ifdef CIS_TRANS_DOES_NOT_WORK
void
Unique_Determination::_expand_around_cis_trans_bonds()
{
  int matoms = _m->natoms();

  int ne = _m->nedges();

  for (int i = 0; i < ne && _cis_trans_bonds > 0; i++)
  {
    const Bond * b = _m->bondi(i);

    if (! b->is_double_bond())
      continue;

    if (! b->part_of_cis_trans_grouping())
      continue;

    _expand_around_cis_trans_bond(b->a1(), b->a2(), rank_delta);

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
  if (nullptr != bnw)
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
  if (nullptr != bsw)
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
  if (nullptr != bne)
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
  if (nullptr != bse)
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
Unique_Determination::_identify_directionally_attached_atoms(const Atom_and_Rank * ar,
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
Unique_Determination::_identify_directionally_attached_bonds(const Atom_and_Rank * ar,
                                                const Bond * & bnw,
                                                const Bond * & bsw) const
{
  bnw = nullptr;
  bsw = nullptr;

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

  if (nullptr != bnw)
    return 1;

  if (nullptr != bsw)
    return 1;

  return 0;
}
#endif

// Using SINGLE_ARRAY is much faster.
void
Unique_Determination::_assign_initial_ranks()
{
  assert (0 == _number_elements);

  this->resize(_matoms);

#define SINGLE_ARRAY
#ifdef SINGLE_ARRAY
  AtomPropertiesForRanking *  x = new AtomPropertiesForRanking[_matoms]; std::unique_ptr<AtomPropertiesForRanking[]> free_x(x);

  resizable_array<AtomPropertiesForRanking *> target;
#else
  resizable_array_p<AtomPropertiesForRanking> target;
#endif

  target.resize(_matoms);

  int molecule_contains_chirality = 0;
  _nchiral = 0;

  assert (_include_chiral_info_in_smiles == include_chiral_info_in_smiles());

  if (include_chiral_info_in_smiles() && (_nchiral = _m->chiral_centres()) > 0)
  {
    for (int i = 0; i < _matoms; i++)
    {
      const Chiral_Centre * c = _m->chiral_centre_at_atom(i);

#ifdef SINGLE_ARRAY
      AtomPropertiesForRanking * a = x + i;
      a->initialise(*_m, i, c);
#else
      AtomPropertiesForRanking * a = new AtomPropertiesForRanking(*_m, i, c);
#endif

      target.add(a);
    }
  }
  else
  {
    for (int i = 0; i < _matoms; i++)
    {
#ifdef SINGLE_ARRAY
      AtomPropertiesForRanking * a = x + i;
      a->initialise(*_m, i, nullptr);
#else
      AtomPropertiesForRanking * a = new AtomPropertiesForRanking(*_m, i, nullptr);
#endif

      target.add(a);
    }
  }

  molecule_contains_chirality = _nchiral;

  target.sort(AtomPropertiesComparitor);

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After sorting targets\n";
  for (int i = 0; i < _matoms; i++)
  {
    const AtomPropertiesForRanking * t = target[i];
    const atom_number_t a = t->atom_number();
    cerr << "i = " << i << " is atom " << a << " (" << _m->smarts_equivalent_for_atom(a) << ") hash " << t->hash() << '\n';
  }
#endif

  _rank_in_use.reset();
  int rank_to_assign = _rank_in_use.IdentifyUnusedRank();

//cerr << "Assigining ranks to " << _matoms << " atoms\n";
  for (int i = 0; i < _matoms; i++)
  {
    int tmp;
    if (i > 0)    // only compare with previous if there is a previous
      tmp = AtomPropertiesComparitor(&(target[i]), &(target[i - 1]));
    else
      tmp = 0;

//  cerr << "Comparison outcome " << tmp << endl;

    if (tmp < 0)
    {
      cerr << "Atoms out of order\n";
      target[i]->debug_print(cerr);
      target[i - 1]->debug_print(cerr);
    }
    assert (tmp >= 0);    // the array is supposed to be sorted

    if (tmp > 0)     // different from one before, change rank
      rank_to_assign  = _rank_in_use.IdentifyUnusedRank();

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
    if (nullptr == _atom_xref[i])
    {
      cerr << "Yipes, the point for atom " << i << " is nullptr\n";
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
  Molecule mcopy(*_m);
  write_atom_map_number_labelled_smiles(mcopy, false, cerr) << '\n';
  debug_print(cerr);
#endif

  return;
}


/*
  Use a molecule_to_Match object to order the atoms in the molecule.
  This is kind of sub-optimal in that we have to do extra
  comparisons after the sort to assign the ranks
  Obsolete, retained for comatibility only.
*/

void
Unique_Determination::_assign_initial_ranks_legacy()
{
  assert (0 == _number_elements);

  resizable_array_p<Target_Atom> target;
  target.resize(_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    Target_Atom * a = new Target_Atom();

    a->initialise(_m, i, const_cast<Atom *>(_m->atomi(i)), nullptr);

    target.add(a);
  }

  target.sort(target_atom_comparitor);

#ifdef DEBUG_ASSIGN_INITIAL_RANKS
  cerr << "After sorting targets\n";
  for (int i = 0; i < _matoms; i++)
  {
    const Target_Atom * t = target[i];
    atom_number_t a = t->atom_number();
    cerr << "i = " << i << " is atom " << a << " (" << _m->atomic_symbol(a) << ") " <<
            t->ncon() << " connections\n";
  }
#endif

  _rank_in_use.reset();
  int rank_to_assign = _rank_in_use.IdentifyUnusedRank();

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
      rank_to_assign = _rank_in_use.IdentifyUnusedRank();

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
    if (nullptr == _atom_xref[i])
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
  Molecule mcopy(*_m);
  write_atom_map_number_labelled_smiles(mcopy, false, cerr) << '\n';
  debug_print(cerr);
#endif

  return;
}

/*
  I have deliberately duplicated the code for assigning initial ranks with a subset
  in order to avoid slowing down the version that doesn't do subsetting
*/

void
Unique_Determination::_assign_initial_ranks_legacy(const int * include_atom)
{
  assert (0 == _number_elements);

  resizable_array_p<Target_Atom> target;
  target.resize(_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    Target_Atom * s = new Target_Atom();

    Atom * a = const_cast<Atom *>(_m->atomi(i));   // loss of const OK

    s->initialise(_m, i, a, nullptr);

    if (! include_atom[i])
    {
      s->set_element(nullptr);
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

  _rank_in_use.reset();

//set_vector(_atom_xref, _m->natoms(), static_cast<Atom_and_Rank *>(nullptr));

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

//  if (tmp > 0 || kInvalidAtomicNumber == target[i]->atomic_number())     // different from one before, increment rank
    if (tmp > 0 || 0 == include_atom[i])   // different from one before, increment rank
      rank_to_assign = _rank_in_use.IdentifyUnusedRank();

    atom_number_t a = target[i]->atom_number();

    Atom_and_Rank * r = new Atom_and_Rank(a, rank_to_assign, target[i]->ncon());

    if (! _include_chiral_info_in_smiles)
      ;
    else if (target[i]->ncon() < 3)
      ;
    else
    {
      const Chiral_Centre * c = _m->chiral_centre_at_atom(a);
      if (nullptr == c)
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
    if (nullptr == _atom_xref[i])
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
Unique_Determination::_get_rank(atom_number_t a) const
{
  if (INVALID_ATOM_NUMBER == a)
    return -1;

  const Atom_and_Rank * ar = _atom_xref[a];

  return ar->rank();
}

int
Unique_Determination::_index_if_active(atom_number_t a) const
{
  for (int i = 0; i < _nactive; i++)
  {
    if (a == _things[i]->atom_number())
      return i;
  }

  return -1;
}

/*
   NW        NE                         *
     \      /                           *
      C == C
     /      \                           *
   SW        SE                         *

   In this function A is the LHS carbon atom
*/

unsigned int
Unique_Determination::_compute_cis_trans_rank(Atom_and_Rank * a, 
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

#ifdef THIS_IS_NOT_USED
static int
identify_directionally_attached_atoms(const Atom * a,
                                      const atom_number_t zatom,
                                      const atom_number_t & nw,
                                      const atom_number_t & sw)
{
  nw = kInvalidAtomicNumber;
  sw = kInvalidAtomicNumber;

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

  if (nw == kInvalidAtomicNumber && sw == kInvalidAtomicNumber)
    return 0;

  return 1;
}
#endif


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
Unique_Determination::_index_of_atom(atom_number_t zatom) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (zatom == _things[i]->atom_number())
      return i;
  }

  return -1;
}

int
Unique_Determination::canonical_order(Molecule & m,
                                      int * canonical_rank,
                                      const int * include_atom)
{
#ifdef DEBUG_UNIQUE_DETERMINATION
  cerr << "Begin canonical order computation\n";
  if (include_atom != nullptr) {
    for (int i = 0; i < m.natoms(); ++i) {
      cerr << "include_atom[" << i <<"] " << include_atom[i] << '\n';
    }
  }
#endif
  assert (nullptr != canonical_rank);
  assert (m.ok());

  (void) _initialise(m, include_atom);

  if (_matoms < 2) {
    if (m.natoms() == 0) {
      return 1;
    }
    canonical_rank[0] = 1;
    return 1;
  }

  int rc = _canonical_order(0);    // 0 means don't stop when symmetry is perceived

  copy_vector(canonical_rank, _canonical_rank, _matoms);

#ifdef DEBUG_UNIQUE_DETERMINATION
  cerr << "Canonical order is\n";
  _print_canonical_order_by_canonical_order(m, include_atom, std::cerr);
#endif

  return rc;
}

int
Unique_Determination::_print_canonical_order(const Molecule & m,
                        const int * include_atom,
                        std::ostream& output) const {
  for (int i = 0; i < _matoms; ++i) {
    output << " atom " << i << " (" << std::setw(2) << m.atomic_symbol(i) << ' ' <<
                   m.ncon(i) << " connections, " << m.nbonds(i) << " bonds) " <<
                   _canonical_rank[i] << " symmetry " << _symmetry[i];
    if (nullptr != include_atom)
      cerr << " inc " << include_atom[i];
    cerr << '\n';
  }

  return output.good();
}

int
Unique_Determination::_print_canonical_order_by_canonical_order(const Molecule & m,
                        const int * include_atom,
                        std::ostream& output) const {
  std::unique_ptr<int[]> xref(std::make_unique<int[]>(_matoms));
  std::iota(xref.get(), xref.get() + _matoms, 0);
  std::sort(xref.get(), xref.get() + _matoms,
      [&](int i1, int i2) {
        return _canonical_rank[i1] < _canonical_rank[i2];
      });

  for (int i = 0; i < _matoms; ++i) {
    atom_number_t a = xref[i];
    output << i << " atom " << a << " (" << std::setw(2) << m.atomic_symbol(a) << ' ' <<
                   m.ncon(a) << " connections, " << m.nbonds(a) << " bonds) " <<
                   _canonical_rank[a] << " symmetry " << _symmetry[a];
    if (nullptr != include_atom)
      cerr << " inc " << include_atom[a];
    cerr << '\n';
  }

  return output.good();
}

int
Molecule::compute_canonical_ranking(Symmetry_Class_and_Canonical_Rank & sccr,
                                    const int * include_atom)
{
  assert (ok());

  if (0 == _number_elements)
    return 1;

  compute_aromaticity_if_needed();

  bool legacy_atom_ordering = file_scope_use_legacy_initial_atom_ordering;
  if (! _ok_for_fast_atom_comparisons())
    legacy_atom_ordering = true;

//cerr << "legacy_atom_ordering " << legacy_atom_ordering << endl;

  if (! sccr.arrays_allocated())
    sccr.allocate_arrays(_number_elements);

  Unique_Determination unqd(legacy_atom_ordering);

  int rc = unqd.canonical_order(*this, sccr.canonical_rank(), include_atom);

  unqd.symmetry(_number_elements, sccr.symmetry_class());    // must do this after the canonical ranking

  return rc;
}

int
Molecule::compute_canonical_ranking()
{
  return compute_canonical_ranking(_symmetry_class_and_canonical_rank, nullptr);
}

int
Molecule::canonical_rank(atom_number_t a)
{
  assert(ok_atom_number(a));

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  return _symmetry_class_and_canonical_rank.canonical_rank(a);
}

int
Molecule::canonical_ranks(int * r)
{
  if (0 == _number_elements)
    return 0;

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  copy_vector(r, _symmetry_class_and_canonical_rank.canonical_rank(), _number_elements);

  return _number_elements;
}

const int *
Molecule::canonical_ranks()
{
  if (0 == _number_elements)
    return nullptr;

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  return _symmetry_class_and_canonical_rank.canonical_rank();
}

const int *
Molecule::symmetry_classes() 
{
  if (0 == _number_elements)
    return nullptr;

  if (! _symmetry_class_and_canonical_rank.arrays_allocated())
    compute_canonical_ranking();

  assert (nullptr != _symmetry_class_and_canonical_rank.symmetry_class());

  return _symmetry_class_and_canonical_rank.symmetry_class();
}

int
Molecule::symmetry_class(atom_number_t a)
{
  assert (ok_atom_number(a));

  (void) symmetry_classes();

  return _symmetry_class_and_canonical_rank.symmetry_class(a);
}

int
Molecule::symmetry_equivalents(atom_number_t a, Set_of_Atoms & sym)
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
Molecule::number_symmetry_classes()
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
Molecule::bond_symmetry_class_small_memory(int * s)
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
Molecule::bond_symmetry_class_large_memory(int * s)
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
