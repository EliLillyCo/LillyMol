#include <algorithm>
#include <iostream>
#include <memory>

/*
  In order to get the private Molecule:: functions in this file,
  make sure we define this to the preprocessor before including molecule.h
*/

#define COMPILING_AROMATIC_CC
#define COMPILING_MOLECULER_CC

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "aromatic.h"
#include "charge_calculation.h"
#include "misc2.h"
#include "molecule.h"
#include "path.h"
#include "smiles.h"
#include "temp_detach_atoms.h"
#include "toggle_kekule_form.h"

using std::cerr;
using std::endl;

/*
  Enable dual nature aromaticity storage.
  That is, record both aliphatic and aromatic nature of some atoms
  to aid in substructure things.
*/

static int record_dual_nature_aromaticity = 0;

static int file_scope_display_no_kekule_form_message = 1;

void
set_display_no_kekule_form_message(int s)
{
  file_scope_display_no_kekule_form_message = s;
}

static int max_aromatic_ring_size = 8;

void
set_max_aromatic_ring_size(int s)
{
  assert(s >= 6);

  max_aromatic_ring_size = s;

  return;
}

static int min_aromatic_ring_size = 4;

void
set_min_aromatic_ring_size(int s) {
  min_aromatic_ring_size = s;
}

int
display_no_kekule_form_message()
{
  return file_scope_display_no_kekule_form_message;
}

static int perform_kekule_perception = 1;

void
set_perform_kekule_perception(int s)
{
  perform_kekule_perception = s;
}

/*
  Just allow these Pipeline Pilot smiles to come in, who cares
*/

static int allow_pipeline_pilot_aromaticity_on_input = 1;

void
set_allow_pipeline_pilot_aromaticity_on_input(int s)
{
  allow_pipeline_pilot_aromaticity_on_input = s;
}

static int strongly_fused_rings_can_be_aromatic = 1;

void
set_strongly_fused_rings_can_be_aromatic(int s)
{
  strongly_fused_rings_can_be_aromatic = s;
}

static int file_scope_allow_any_even_number_of_pi_electrons = 0;

void
set_allow_any_even_number_of_pi_electrons(int s)
{
  file_scope_allow_any_even_number_of_pi_electrons = s;
}

static int file_scope_aromatic_rings_must_contain_unsaturation = 1;

void
set_aromatic_rings_must_contain_unsaturation(const int s)
{
  file_scope_aromatic_rings_must_contain_unsaturation = s;
}

int
Molecule::_allocate_aromaticity()
{
  if (nullptr == _aromaticity)
  {
    if (0 == _number_elements)
      return 1;

    _aromaticity = (aromaticity_type_t *)new_int(_number_elements, AROMATICITY_NOT_DETERMINED);

    assert(nrings() >= 0);

    return 1;
  }

  return 0;
}

int
Molecule::aromaticity(atom_number_t a, aromaticity_type_t & result)
{
  assert(ok_atom_number(a));

  if (nullptr == _aromaticity)
    _compute_aromaticity();

  result = _aromaticity[a];

  if (AROMATICITY_NOT_DETERMINED == _aromaticity[a])    // must be in a strange ring
    return 0;

  return 1;
}

int
Molecule::is_aromatic(atom_number_t a)
{
  aromaticity_type_t arom;

  if (! aromaticity(a, arom))
    return 0;

  // Should use this, but certain atom typing codes depend on the old behaviour.
  // return is_aromatic_atom(arom);  // correct, but breaks compatibility.

  // Retained for compatibility.
  return _AROM_BIT & arom;
}

bool
Molecule::IsAromatic(atom_number_t a) {
  aromaticity_type_t arom;

  if (! aromaticity(a, arom)) {
    return false;
  }

  return is_aromatic_atom(arom);
}

int
Molecule::is_aromatic_no_computation(atom_number_t zatom) const
{
  if (nullptr == _aromaticity)
    return -1;

  return is_aromatic_atom(_aromaticity[zatom]);
}

int
Molecule::is_permanent_aromatic(const atom_number_t & a) const
{
  assert(ok_atom_number(a));

  return _things[a]->permanent_aromatic();
}

int
Molecule::aromaticity(int * result)
{
  if (nullptr == _aromaticity)
    _compute_aromaticity();

  int rc = 1;    // unless there is a problem
  for (int i = 0; i < _number_elements; i++)
  {
    if (AROMATICITY_NOT_DETERMINED == _aromaticity[i])
    {
      rc = 0;
      result[i] = 0;
    }
    else if (is_aromatic_atom(_aromaticity[i]))
      result[i] = 1;
    else
      result[i] = 0;
  }

  return rc;
}

int
Molecule::compute_aromaticity_if_needed()
{
  if (nullptr != _aromaticity)
    return 0;

  return _compute_aromaticity();
}

int
Molecule::contains_aromatic_atoms()
{
  assert(ok());

  if (0 == _number_elements)
    return 0;

  if (nullptr == _aromaticity)
    _compute_aromaticity();

  for (int i = 0; i < _number_elements; i++)
  {
    if (is_aromatic_atom(_aromaticity[i]))
      return 1;
  }

  return 0;
}

int
Molecule::aromatic_atom_count()
{
  assert(ok());

  if (0 == _number_elements)
    return 0;

  if (nullptr == _aromaticity)
    _compute_aromaticity();

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (is_aromatic_atom(_aromaticity[i]))
      rc++;
  }

  return rc;
}

int
Molecule::aromatic_ring_count()
{
  if (0 == _number_elements)
    return 0;

  const int nr = nrings();
  if (0 == nr)
    return 0;

  if (nullptr == _aromaticity)
    _compute_aromaticity();

  int rc = 0;
  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = _sssr_rings[i];

    if (ri->is_aromatic())
      rc++;
  }

  return rc;
}

/*
  A quick test to see whether or not aromaticity has been computed.
  Note that this is not watertight, as we only check one atom. If
  ever we change to computing aromaticity as needed, we would need
  to change this.
*/

int
Molecule::aromaticity_computed() const
{
  assert(ok());

  if (nullptr == _aromaticity)
    return 0;

  return AROMATICITY_NOT_DETERMINED != _aromaticity[0];
}

int
Molecule::lone_pair_count(atom_number_t a, int & result)
{
  assert(ok_atom_number(a));

  return _things[a]->lone_pair_count(result);
}

int
Molecule::pi_electrons(atom_number_t zatom, int & pe)
{
  assert(ok_atom_number(zatom));

  return _things[zatom]->pi_electrons(pe);
}

/*
  this option is only used when inputting Daylight derived aromatic
  smiles.
  Oct 2022. Make default.
*/

static int _allow_two_electron_systems_to_be_aromatic = 1;

int
allow_two_electron_systems_to_be_aromatic()
{
  return _allow_two_electron_systems_to_be_aromatic;
}

/*
  Sept 97, ran into problems with Ludi producing aromatic bonds in
  chains, so we have an option to suppress 'aromatic' chain bonds
*/

static int x_convert_chain_aromatic_bonds = 0;

int
convert_chain_aromatic_bonds()
{
  return x_convert_chain_aromatic_bonds;
}

void
set_convert_chain_aromatic_bonds(int s)
{
  x_convert_chain_aromatic_bonds = s;
}

static int _aromatic_chain_bonds_are_ok = 0;

int
aromatic_chain_bonds_are_ok()
{
  return _aromatic_chain_bonds_are_ok;
}

void
set_aromatic_chain_bonds_are_ok(int s)
{
  _aromatic_chain_bonds_are_ok = s;
}

/*
  Do we allow things like c1cccc1 to be aromatic?
*/

static int _non_kekule_systems_ok_to_be_aromatic = 0;

void
set_non_kekule_systems_ok_to_be_aromatic(int s)
{
  _non_kekule_systems_ok_to_be_aromatic = s;
}

int
non_kekule_systems_ok_to_be_aromatic()
{
  return _non_kekule_systems_ok_to_be_aromatic;
}

static int global_aromaticity_determination_type = Daylight;

int
set_global_aromaticity_type(int na)
{
  if (na == global_aromaticity_determination_type)    // not being changed
    return 1;

  if (Simple_4n_plus_2 == na)
    global_aromaticity_determination_type = Simple_4n_plus_2;
  else if (Daylight == na)
    global_aromaticity_determination_type = Daylight;
  else if (Pearlman == na)
    global_aromaticity_determination_type = Pearlman;
  else if (WangFuLai == na)
    global_aromaticity_determination_type = WangFuLai;
  else if (Vijay_Gombar == na)
    global_aromaticity_determination_type = Vijay_Gombar;
  else if (Pipeline_Pilot == na)
    global_aromaticity_determination_type = Pipeline_Pilot;
  else if (PUBCHEM_AROMATICITY == na)
    global_aromaticity_determination_type = PUBCHEM_AROMATICITY;
  else if (EVERYTHING_HAS_A_PI_ELECTRON == na)
    global_aromaticity_determination_type = EVERYTHING_HAS_A_PI_ELECTRON;
  else if (na == ANY_EVEN_NUMBER_OF_PI_ELECTRONS) {
    global_aromaticity_determination_type = ANY_EVEN_NUMBER_OF_PI_ELECTRONS;
  } else {
    cerr << "set_global_aromaticity_type: unknown aromaticity type " << na << '\n';
    return 0;
  }

  return 1;
}

int
global_aromaticity_type()
{
  return global_aromaticity_determination_type;
}

/*
  Jun 2004. Ran into cases from Skelgen where some of the bonds
  in an aromatic ring were not marked as aromatic bonds. This
  became impossible because the atoms were not marked as C.ar,
  but as C.2

*/

static int all_bonds_in_aromatic_ring_must_be_aromatic = 1;

void
set_all_bonds_in_aromatic_ring_must_be_aromatic(int s)
{
  all_bonds_in_aromatic_ring_must_be_aromatic = s;
}

// Return true if `zatom` is set in `arom_data.in_all_pi_ring`.
int
Molecule::_in_all_pi_ring(AromData& arom_data, atom_number_t zatom)
{
  if (arom_data.in_all_pi_ring == nullptr) {
    _determine_in_all_pi_ring(arom_data);
  }

  return arom_data.in_all_pi_ring[zatom];
}

// `in_ring` is an atom in a ring.
// `maybe_n` is doubly bonded to `in_ring`, but outside its ring (but might be in a ring).
//  Does `x` look like an amide...
int
Molecule::_amide_like(atom_number_t in_ring,
                      atom_number_t maybe_n) const {
  const Atom * n = _things[maybe_n];
  if (n->atomic_number() != 7) {
    return 0;
  }

  if (n->ncon() != 2) {
    return 0;
  }

  for (const Bond * b : *n) {
    if (b->is_double_bond()) {  // Looping back to the ring..
      continue;
    }

    const atom_number_t c = b->other(maybe_n);
    // Amide or sulfonamide OK
    if (_things[c]->atomic_number() == 6)
      ;
    else if (_things[c]->atomic_number() == 16)
      ;
    else
      return 0;

    if (_things[c]->ncon() < 3) {
      return 0;
    }

    for (const Bond * b2 : *_things[c]) {
      if (! b2->is_double_bond()) {
        continue;
      }
      const atom_number_t o = b2->other(c);
      if (_things[o]->atomic_number() == 8 || _things[o]->atomic_number() == 16) {
        return 1;
      }
      // Found cases where there is [R]=N-C=C, why not...
      if (_things[o]->atomic_number() == 6) {
        return 1;
      }

      return 0;
    }
    return 0;
  }

  return 0;
}

// Identify rings in which every atoms has pi electrons. For every
// atom in such a ring, set the corresponding value in `arom_data.in_all_pi_ring`.
// Return the number of such rings.
int
Molecule::_determine_in_all_pi_ring(AromData& arom_data)
{
  arom_data.in_all_pi_ring = new_int(_number_elements);

  int rc = 0;
  const int nr = nrings();
  for (int i = 0; i < nr; ++i) {
    const Ring * r = ringi(i);
    const int ring_size = r->number_elements();
    int has_pi_electrons = 0;  // Count atoms with pi electrons.
    for (int j = 0; j < ring_size; ++j) {
      const atom_number_t k = r->item(j);
      if (arom_data.in_all_pi_ring[k]) {  // Already identified.
        ++has_pi_electrons;
        continue;
      }

      int pi;
      pi_electrons(k, pi);  // Do not check success/failure.

      // Assume non-organics contribute pi electrons....
      if (pi > 0 || ! _things[k]->element()->organic()) {
        ++has_pi_electrons;
      }
    }

    if (has_pi_electrons == r->number_elements()) {
      r->set_vector(arom_data.in_all_pi_ring, 1);
      ++rc;
    }
  }

  return rc;
}

/*
  Depends on ring membership having been called - so the bonds know their ring membership
*/

int
Molecule::_doubly_bonded_to_something_outside_ring(atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_single_bond())
      continue;

    if (0 == b->nrings())
      return 1;
  }

  return 0;
}

/*
  Need to differentiate two kinds of Phosphorus environments
    P1(C(=P(=C1C1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C)C)C1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C)C(C)(C)C PBCHM101753934

    O=[p]1[nH][pH][n]1 PBCHM22761484
    O=P1=NPN1

  In the first one, both double bonds are in the ring, whereas in the second one, there is a non ring double bond
*/

int
Molecule::_both_double_bonds_in_set(const Set_of_Atoms & p, const atom_number_t zatom) const
{
  const Atom * a = _things[zatom];
  const int acon = a->ncon();
  assert(3 == acon);

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    const atom_number_t o = b->other(zatom);

    if (8 != _things[o]->atomic_number())
      continue;

    if (0 == b->nrings())    // definitely not part of a ring system
      return 0;

    if (! p.contains(
            o))    // a little tricky here, perhaps it is in a ring, but not being processed right now, hard to know what is right
      return 0;
  }

  return 1;
}

//#define DEBUG_AROMATICITY_DETERMINATION

/*
  Aromaticity detection.
  Note that argument P is not necessarily a single ring, but can be all atoms in
  a fused ring system. Therefore we make no assumptions about adjacent entries
  being bonded, ring membership or anything like that.
  If we discern that this ring could not possibly be aromatic, then we
  set argument IMPOSSIBLE to 1
*/

int
Molecule::_determine_aromaticity(const Set_of_Atoms & p,
                                 int& impossible_aromatic,
                                 aromaticity_type_t & result,
                                 AromData& arom_data)
{
  int np = p.number_elements();
  assert(np > 2);

  impossible_aromatic = 0;

#ifdef DEBUG_AROMATICITY_DETERMINATION
  cerr << "Begin aromaticity determination for " << p << " arom " << arom_data.aromaticity_rule << '\n';
#endif

  int unsaturation = 0;

  int pe = 0;    // accumulator for pi electrons
  for (int i = 0; i < np; i++)
  {
    atom_number_t j = p[i];

    Atom * aj = _things[j];

#ifdef DEBUG_AROMATICITY_DETERMINATION
    cerr << "Atom " << j << " check > 3 connected ncon = " << aj->ncon() << " + "
         << aj->implicit_hydrogens() << " = " << aj->ncon() + aj->implicit_hydrogens() << '\n';
#endif

    if (aj->element()->organic()) {  // great
    } else if (aj->permanent_aromatic()) {
    } else {
      impossible_aromatic = 1;
      return 0;
    }

    if (! aj->valence_ok()) {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 0;
    }

    int jcon = aj->ncon();
    if (jcon + aj->implicit_hydrogens() > 3)
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;    // this is a successful determination.
    }

#ifdef DEBUG_AROMATICITY_DETERMINATION
    int foo;
    cerr << "Atom " << j << " (atomic number " << atomic_number(j) << ")";
    if (pi_electrons(j, foo))
      cerr << " has " << foo << " pi electrons\n";
    else
      cerr << " cannot determine pi electron count\n";
    cerr << " precomputed " << arom_data.pi_electrons[j] << '\n';
#endif

    //  There is some complexity associated with the pi electron count.
    //  If the value in pi_electron_count is valid, we want to use it.
    //  Otherwise, we must compute the value here.

    int tmp;
    if (arom_data.pi_electrons[j] >= 0)    // value already computed
      tmp = arom_data.pi_electrons[j];
    else if (! aj->pi_electrons(tmp))    // try to compute it
    {
      impossible_aromatic = 1;    // we don't know how!!
      return 0;
    }
    else
    {
      //    cerr << "Computing pi electron count for " << j << " type " << smarts_equivalent_for_atom(j) << " value " << tmp << '\n';
      if (2 == jcon && 6 == aj->atomic_number() &&
          4 == aj->nbonds())    // Oct 2001, the ring N1=C=N-C=C1 was aromatic
      {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }
      arom_data.pi_electrons[j] = tmp;
    }

    //  If it has 0 pi electrons, and all single bonds, then this ring is non aromatic
    //  Mostly this is to catch saturated carbon, but may also catch others...

#ifdef DEBUG_AROMATICITY_DETERMINATION
    cerr << "Atom " << j << " has " << tmp << " pi electrons, ncon = " << ncon(j)
         << " nbonds() = " << nbonds(j) << '\n';
#endif

    const int jbonds = aj->nbonds();
    if (0 == tmp && jcon == jbonds)    // no pi electrons and fully saturated
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    if (jbonds > jcon)
      unsaturation++;

    const atomic_number_t zj = aj->atomic_number();

    //  Jun 97, ran into problems with a 3 connected iodine in a ring. Reject them
    //  Jul 2007 extend to Cl and Br
    // Aug 2023. TODO:ianwatson allow rings like this to be aromatic.
    // S1(=NC(=NC2=CC(=C(C=C12)OC)OC)N1CCN(C(=O)C2=CC=CO2)CC1)(=O)C1=CC=CC=C1 CHEMBL1196444

    if (6 == zj || 7 == zj)    // the most common cases
      ;
    else if (53 == zj || 35 == zj || 17 == zj)
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }
#ifdef SD3V4_NOT_AROMATIC
    else if (
        16 == zj && 3 == jcon &&
        4 == jbonds)  // Jan 2017 S12=NN(O)NC1=CC(=NC2=C)C PBCHM57260347 - do not aromatise. But what about CC1=N(=O)C=CC=C1
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }
#endif
    // P1(C(=P(=C1C1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C)C)C1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C)C(C)(C)C PBCHM101753934
    else if ( 15 == zj && 3 == jcon && 5 == jbonds && _both_double_bonds_in_set( p, j)) {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    //  Sept 2007. Do not aromatise O1N=C2C3=C(N4N(=N3)=C3C(=N4)C=CC=C3)C=CC2=N1 PBCHM389865
    //  and [N+]12-[N-]N(C)N=C1C=CC=2 PBCHM10606729. but would be nice to have these aromatic...

    if (7 == zj && 3 == aj->ncon() && 0 == aj->formal_charge() && 5 == aj->nbonds() &&
        ! _doubly_bonded_to_something_outside_ring(j)) {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    //  According to Pearlman's rules, 3 connected N automatically disqualifies
    //  On 30 Aug 96 I changed this to only exclude 3 connected N when it had 3
    //  bonds. The problem came about for some molecules which had an exocyclic =O
    //  bond, which Concord considers aromatic

    if (Pearlman == arom_data.aromaticity_rule || WangFuLai == arom_data.aromaticity_rule)
    {
      if (7 == zj && 0 == aj->formal_charge() &&
          ((3 == jcon && 3 == jbonds) ||    // 30 aug added the condition '3 == nbonds(j)'
           (2 == jcon && 1 == hcount(j))))
      {
#ifdef DEBUG_AROMATICITY_DETERMINATION
        cerr << "Found 3 connected N in ring\n";
#endif

        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }
    }

    //  Jan 2001. The algorithm was aromatising the ring
    //  Don't aromatise any 5 connected Nitrogen with two double bonds in the ring

    // Aug 2002. Wow, this is hard. If you draw the ring in its charge separated form, then 
    // it is aromatic. But in the form drawn, Daylight doesn't consider it aromatic.
    // For now I'm leaving it non-aromatic..., but not happy about this.. It affects
    // a relatively small number of molecules

    // Don't forget about N2=CSC(=NC#N)N=2

    if (7 != zj)    // only interested in Nitrogen atoms here
      ;
    else if (0 != aj->formal_charge())
    {
      if (1 == aj->formal_charge() && 2 == aj->ncon() &&
          4 == aj->nbonds())    // N1C(=O)N=[N+]=CC1=O PBCHM56995681
      {
        result = NOT_AROMATIC;
        return 1;
      }
    }
    else if (jcon == jbonds)    // fully saturated, not interested
      ;
    else if (2 == jcon && 4 == jbonds)    // case of ring from N=C1N=N(C)=CS1
    {
      result = NOT_AROMATIC;
      return 1;
    }
    else if (3 == jcon && 5 == jbonds)    // if both bonds in ring, not aromatic
    {
      // Look for multiple doubly bonded atoms in the system.
      atom_number_t double_bond_1_in_system = INVALID_ATOM_NUMBER;
      atom_number_t double_bond_2_in_system = INVALID_ATOM_NUMBER;

      for (int k = 0; k < jcon; k++)
      {
        const Bond * b = aj->item(k);
        if (b->is_single_bond())
          continue;

        atom_number_t ak = b->other(j);

        // If there is a double bond outside the system, there cannot be two in the ring.
        if (! p.contains(ak))
          break;

        if (INVALID_ATOM_NUMBER == double_bond_1_in_system)
          double_bond_1_in_system = ak;
        else
          double_bond_2_in_system = ak;
      }

      if (INVALID_ATOM_NUMBER == double_bond_2_in_system)
        ;
      else if (in_same_ring(double_bond_1_in_system, double_bond_2_in_system))
      {
        result = NOT_AROMATIC;
        return 1;
      }
    }

    //  Sept 96. Ran into a problem with -N=S=N- in a ring, which got aromatised,
    //  so we create a special case to suppress that.
    //  Also found =S- as a problem too
    //  Oct 2000. [S+] can be aromatic
    //  Jul 2023. Make sure all of these are aromatic
    // S1=CC=CC=C1 p2.a
    // [SH]1=CC=CC=C1 p2.b
    // CS1=CC=CC=C1 p2.c

    if (16 == zj) {
      if (jcon == 3 && jbonds == 4 && aj->formal_charge() == 0) {
      } else if (jcon == 2 && jbonds == 3 && hcount(j) == 1) {
      } else if ((2 == jcon && 4 == jbonds) || (2 == jcon && 3 == jbonds && 0 == aj->formal_charge()))
      {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }

      // dec 2002, Sulphur never aromatic in Pearlman's world
      if (2 == jcon && 2 == jbonds &&
          (Pearlman == arom_data.aromaticity_rule || arom_data.aromaticity_rule == WangFuLai)) {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }

      // Jan 2006: MDDR242489. make it not aromatic
      if (3 == jcon && 4 == jbonds && _doubly_bonded_to_something_outside_ring(j)) {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }
    }

    //  Jun97, Pearlman disqualifies any ring Oxygen as aromatic - furan for example

    if (8 == zj && 
        (Pearlman == arom_data.aromaticity_rule || arom_data.aromaticity_rule == WangFuLai)) {
      //impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    //  If the atom looks like C=O, then that is OK to have no pi electrons

    if (EVERYTHING_HAS_A_PI_ELECTRON == arom_data.aromaticity_rule) {
      if (tmp > 0)
        ;
      else if (3 == jcon && 6 == aj->atomic_number() && 1 == doubly_bonded_oxygen_count(j))
        ;
      else
      {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }

      pe += tmp;
      continue;
    }

    if (PUBCHEM_AROMATICITY == arom_data.aromaticity_rule &&
        _doubly_bonded_to_something_outside_ring(j))
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    // Does this atom have a double bond to an atom not in P. We decrement
    // the pi electron count if there is a double bond to an atom of
    // different type. Also, check for triple bonds in the ring.

    if (jcon < jbonds)    // there is at least one multiple bond
    {
#ifdef DEBUG_AROMATICITY_DETERMINATION
      cerr << "Checking for multiple bonds outside the ring from atom " << j << '\n';
#endif

      for (int k = 0; k < jcon; k++)
      {
        const Bond * b = aj->item(k);

        if (b->is_single_bond())
          continue;

        if (b->is_triple_bond())
        {
          impossible_aromatic = 1;
          result = NOT_AROMATIC;
          return 1;
        }

        atom_number_t l = b->other(j);

        const Atom * al = _things[l];

        // Apr 2004. Conversations with Vijay. c=C does not contribute

        if (Vijay_Gombar == arom_data.aromaticity_rule && 6 == zj && 6 == al->atomic_number() &&
            ! p.contains(l)) {
          tmp--;
          continue;
        }

        //      The atom types must be different, and either
        //        atom l is non ring
        //        atom l is not in this ring, nor are atoms j and l in any ring together.

        //      Note an extra wrinkle here. In order to maximise our ability to read
        //      Daylight "aromatic" smiles, we only decrement the pi count if the
        //      ring atom is a carbon.

        if (zj == al->atomic_number())    // both atom types the same, ignore
          continue;

        if (p.contains(l))  // Is in the ring being determined.
          continue;

        // Atom J is doubly bonded to something outside the ring

//      cerr << "Consider atom " << j << ' ' << smarts_equivalent_for_atom(j) << " bonded to atom " << l << " " << smarts_equivalent_for_atom(l) << " in ring " << p.contains(l) << '\n';
        // Jan 2022: make sure C1Cn2cn[nH]c2=N1 is handled.
        // This is really hard, lots of things break, needs attention....
        if (1 == al->ncon() || (! _in_all_pi_ring(arom_data, l) && ! _amide_like(j, l)) ||
            is_non_ring_atom(l)) {
          if (16 == zj)    // mar 2004. The S atom in O=C1NS(=O)NC2=CC=CC=C12 PBCHM71359875 contributes electrons
            ;
          else
            tmp--;

          if (Pearlman == arom_data.aromaticity_rule || WangFuLai == arom_data.aromaticity_rule) {
            if (6 == zj &&
                (8 == al->atomic_number() || 16 == al->atomic_number()))    // carbonyl disqualifies
            {
              impossible_aromatic = 1;
              result = NOT_AROMATIC;
              return 1;
            }
          }
        }
      }
    }

    if (tmp < 0)    // happens when parsing very broken smiles
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    pe += tmp;

#ifdef DEBUG_AROMATICITY_DETERMINATION
    cerr << "After adjustments, pi count is " << tmp << " total = " << pe << " finished atom " << j << '\n';
#endif
  }

  if (0 == unsaturation &&
      file_scope_aromatic_rings_must_contain_unsaturation)    // P12P3P1[P-]P(P2)[P-]3 PBCHM15961351
  {
    result = NOT_AROMATIC;
    return 1;
  }

  if (EVERYTHING_HAS_A_PI_ELECTRON == arom_data.aromaticity_rule) {
    result = AROMATIC;
  } else if (Simple_4n_plus_2 == arom_data.aromaticity_rule && 
           allow_two_electron_systems_to_be_aromatic() && 2 == pe) {
    result = AROMATIC;
  } else if (pe > 2 && 2 == pe % 4) {
    result = AROMATIC;
  } else if (pe == 2 && allow_two_electron_systems_to_be_aromatic()) {
    result = AROMATIC;
  } else if (arom_data.aromaticity_rule == ANY_EVEN_NUMBER_OF_PI_ELECTRONS && (pe / 2) * 2 == pe) {
    result = AROMATIC;
  } else {
    result = NOT_AROMATIC;
  }

#ifdef DEBUG_AROMATICITY_DETERMINATION
  cerr << "Pi electron count " << pe << " arom is " << (AROMATIC == result) << '\n';
#endif

  return 1;
}

#ifdef TOO_MANY_ARGS
int
Molecule::_determine_aromaticity(const Set_of_Atoms & p, aromaticity_type_t & result,
                                 int & impossible_aromatic, int * pi_electron_count,
                                 int aromaticity_determination_type)
{
  int np = p.number_elements();
  assert(np > 2);

  impossible_aromatic = 0;

#ifdef DEBUG_AROMATICITY_DETERMINATION
  cerr << "Begin aromaticity determination for " << p << " arom " << aromaticity_determination_type << '\n';
#endif

  int unsaturation = 0;

  int pe = 0;    // accumulator for pi electrons
  for (int i = 0; i < np; i++)
  {
    atom_number_t j = p[i];

    Atom * aj = _things[j];

#ifdef DEBUG_AROMATICITY_DETERMINATION
    cerr << "Atom " << j << " check > 3 connected ncon = " << aj->ncon() << " + "
         << aj->implicit_hydrogens() << " = " << aj->ncon() + aj->implicit_hydrogens() << '\n';
#endif

    if (aj->element()->organic()) {  // great
    } else if (aj->permanent_aromatic()) {
    } else {
      impossible_aromatic = 1;
      return 0;
    }

    if (! aj->valence_ok()) {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 0;
    }

    int jcon = aj->ncon();
    if (jcon + aj->implicit_hydrogens() > 3)
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;    // this is a successful determination.
    }

#ifdef DEBUG_AROMATICITY_DETERMINATION
    int foo;
    cerr << "Atom " << j << " (atomic number " << atomic_number(j) << ")";
    if (pi_electrons(j, foo))
      cerr << " has " << foo << " pi electrons\n";
    else
      cerr << " cannot determine pi electron count\n";
    cerr << " precomputed " << pi_electron_count[j] << '\n';
#endif

    //  There is some complexity associated with the pi electron count.
    //  If the value in pi_electron_count is valid, we want to use it.
    //  Otherwise, we must compute the value here.

    int tmp;
    if (pi_electron_count[j] >= 0)    // value already computed
      tmp = pi_electron_count[j];
    else if (! aj->pi_electrons(tmp))    // try to compute it
    {
      impossible_aromatic = 1;    // we don't know how!!
      return 0;
    }
    else
    {
      //    cerr << "Computing pi electron count for " << j << " type " << smarts_equivalent_for_atom(j) << " value " << tmp << '\n';
      if (2 == jcon && 6 == aj->atomic_number() &&
          4 == aj->nbonds())    // Oct 2001, the ring N1=C=N-C=C1 was aromatic
      {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }
      pi_electron_count[j] = tmp;
    }

    //  If it has 0 pi electrons, and all single bonds, then this ring is non aromatic
    //  Mostly this is to catch saturated carbon, but may also catch others...

#ifdef DEBUG_AROMATICITY_DETERMINATION
    cerr << "Atom " << j << " has " << tmp << " pi electrons, ncon = " << ncon(j)
         << " nbonds() = " << nbonds(j) << '\n';
#endif

    const int jbonds = aj->nbonds();
    if (0 == tmp && jcon == jbonds)    // no pi electrons and fully saturated
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    if (jbonds > jcon)
      unsaturation++;

    const atomic_number_t zj = aj->atomic_number();

    //  Jun 97, ran into problems with a 3 connected iodine in a ring. Reject them
    //  Jul 2007 extend to Cl and Br

    if (6 == zj || 7 == zj)    // the most common cases
      ;
    else if (53 == zj || 35 == zj || 17 == zj)
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }
#ifdef SD3V4_NOT_AROMATIC
    else if (
        16 == zj && 3 == jcon &&
        4 == jbonds)  // Jan 2017 S12=NN(O)NC1=CC(=NC2=C)C PBCHM57260347 - do not aromatise. But what about CC1=N(=O)C=CC=C1
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }
#endif
    else if (
        15 == zj && 3 == jcon && 5 == jbonds &&
        _both_double_bonds_in_set(
            p,
            j))    // P1(C(=P(=C1C1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C)C)C1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C)C(C)(C)C PBCHM101753934
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    //  Sept 2007. Do not aromatise O1N=C2C3=C(N4N(=N3)=C3C(=N4)C=CC=C3)C=CC2=N1 PBCHM389865
    //  and [N+]12-[N-]N(C)N=C1C=CC=2 PBCHM10606729. but would be nice to have these aromatic...

    if (7 == zj && 3 == aj->ncon() && 0 == aj->formal_charge() && 5 == aj->nbonds() &&
        ! _doubly_bonded_to_something_outside_ring(j))
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    //  According to Pearlman's rules, 3 connected N automatically disqualifies
    //  On 30 Aug 96 I changed this to only exclude 3 connected N when it had 3
    //  bonds. The problem came about for some molecules which had an exocyclic =O
    //  bond, which Concord considers aromatic

    if (Pearlman == aromaticity_determination_type || WangFuLai == aromaticity_determination_type)
    {
      if (7 == zj && 0 == aj->formal_charge() &&
          ((3 == jcon && 3 == jbonds) ||    // 30 aug added the condition '3 == nbonds(j)'
           (2 == jcon && 1 == hcount(j))))
      {
#ifdef DEBUG_AROMATICITY_DETERMINATION
        cerr << "Found 3 connected N in ring\n";
#endif

        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }
    }

    //  Jan 2001. The algorithm was aromatising the ring
    //  Don't aromatise any 5 connected Nitrogen with two double bonds in the ring

    // Aug 2002. Wow, this is hard. If you draw the ring in its charge separated form, then it is aromatic. But in
    // the form drawn, Daylight doesn't consider it aromatic. For now I'm leaving it non-aromatic..., but not
    // happy about this.. It affects a relatively small number of molecules

    // Don't forget about N2=CSC(=NC#N)N=2

    if (7 != zj)    // only interested in Nitrogen atoms here
      ;
    else if (0 != aj->formal_charge())
    {
      if (1 == aj->formal_charge() && 2 == aj->ncon() &&
          4 == aj->nbonds())    // N1C(=O)N=[N+]=CC1=O PBCHM56995681
      {
        result = NOT_AROMATIC;
        return 1;
      }
    }
    else if (jcon == jbonds)    // fully saturated, not interested
      ;
    else if (2 == jcon && 4 == jbonds)    // case of ring from N=C1N=N(C)=CS1
    {
      result = NOT_AROMATIC;
      return 1;
    }
    else if (3 == jcon && 5 == jbonds)    // if both bonds in ring, not aromatic
    {
      // Look for multiple doubly bonded atoms in the system.
      atom_number_t double_bond_1_in_system = INVALID_ATOM_NUMBER;
      atom_number_t double_bond_2_in_system = INVALID_ATOM_NUMBER;

      for (int k = 0; k < jcon; k++)
      {
        const Bond * b = aj->item(k);
        if (b->is_single_bond())
          continue;

        atom_number_t ak = b->other(j);

        // If there is a double bond outside the system, there cannot be two in the ring.
        if (! p.contains(ak))
          break;

        if (INVALID_ATOM_NUMBER == double_bond_1_in_system)
          double_bond_1_in_system = ak;
        else
          double_bond_2_in_system = ak;
      }

      if (INVALID_ATOM_NUMBER == double_bond_2_in_system)
        ;
      else if (in_same_ring(double_bond_1_in_system, double_bond_2_in_system))
      {
        result = NOT_AROMATIC;
        return 1;
      }
    }

    //  Sept 96. Ran into a problem with -N=S=N- in a ring, which got aromatised,
    //  so we create a special case to suppress that.
    //  Also found =S- as a problem too
    //  Oct 2000. [S+] can be aromatic

    if (16 == zj)
    {
      if ((2 == jcon && 4 == jbonds) || (2 == jcon && 3 == jbonds && 0 == aj->formal_charge()))
      {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }

      // dec 2002, Sulphur never aromatic in Pearlman's world
      if (2 == jcon && 2 == jbonds && Pearlman == aromaticity_determination_type) {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }

      // Jan 2006: MDDR242489. make it not aromatic
      if (3 == jcon && 4 == jbonds && _doubly_bonded_to_something_outside_ring(j)) {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }
    }

    //  Jun97, Pearlman disqualifies any ring Oxygen as aromatic - furan for example

    if (8 == zj && Pearlman == aromaticity_determination_type)
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    //  If the atom looks like C=O, then that is OK to have no pi electrons

    if (EVERYTHING_HAS_A_PI_ELECTRON == aromaticity_determination_type)
    {
      if (tmp > 0)
        ;
      else if (3 == jcon && 6 == aj->atomic_number() && 1 == doubly_bonded_oxygen_count(j))
        ;
      else
      {
        impossible_aromatic = 1;
        result = NOT_AROMATIC;
        return 1;
      }

      pe += tmp;
      continue;
    }

    if (PUBCHEM_AROMATICITY == aromaticity_determination_type &&
        _doubly_bonded_to_something_outside_ring(j))
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    // Does this atom have a double bond to an atom not in P. We decrement
    // the pi electron count if there is a double bond to an atom of
    // different type. Also, check for triple bonds in the ring.

    if (jcon < jbonds)    // there is at least one multiple bond
    {
#ifdef DEBUG_AROMATICITY_DETERMINATION
      cerr << "Checking for multiple bonds outside the ring from atom " << j << '\n';
#endif

      for (int k = 0; k < jcon; k++)
      {
        const Bond * b = aj->item(k);

        if (b->is_single_bond())
          continue;

        if (b->is_triple_bond())
        {
          impossible_aromatic = 1;
          result = NOT_AROMATIC;
          return 1;
        }

        atom_number_t l = b->other(j);

        const Atom * al = _things[l];

        // Apr 2004. Conversations with Vijay. c=C does not contribute

        if (Vijay_Gombar == aromaticity_determination_type && 6 == zj && 6 == al->atomic_number() &&
            ! p.contains(l))
        {
          tmp--;
          continue;
        }

        //      The atom types must be different, and either
        //        atom l is non ring
        //        atom l is not in this ring, nor are atoms j and l in any ring together.

        //      Note an extra wrinkle here. In order to maximise our ability to read
        //      Daylight "aromatic" smiles, we only decrement the pi count if the
        //      ring atom is a carbon.

        if (zj == al->atomic_number())    // both atom types the same, ignore
          continue;

        if (p.contains(l))  // Is in the ring being determined.
          continue;

        //      Atom J is doubly bonded to something outside the ring

//      cerr << "Consider atom " << j << ' ' << smarts_equivalent_for_atom(j) << " bonded to atom " << l << " " << smarts_equivalent_for_atom(l) << " in ring " << p.contains(l) << '\n';
        // Jan 2022: make sure C1Cn2cn[nH]c2=N1 is handled.
        if (1 == al->ncon() || is_non_ring_atom(l) || (! p.contains(l)))
        {
          if (16 == zj)    // mar 2004. The S atom in O=C1NS(=O)NC2=CC=CC=C12 PBCHM71359875 contributes electrons
            ;
          else
            tmp--;

          if (Pearlman == aromaticity_determination_type ||
              WangFuLai == aromaticity_determination_type)
          {
            if (6 == zj &&
                (8 == al->atomic_number() || 16 == al->atomic_number()))    // carbonyl disqualifies
            {
              impossible_aromatic = 1;
              result = NOT_AROMATIC;
              return 1;
            }
          }
        }
      }
    }

    if (tmp < 0)    // happens when parsing very broken smiles
    {
      impossible_aromatic = 1;
      result = NOT_AROMATIC;
      return 1;
    }

    pe += tmp;

#ifdef DEBUG_AROMATICITY_DETERMINATION
    cerr << "After adjustments, pi count is " << tmp << " total = " << pe << " finished atom " << j << '\n';
#endif
  }

  if (0 == unsaturation &&
      file_scope_aromatic_rings_must_contain_unsaturation)    // P12P3P1[P-]P(P2)[P-]3 PBCHM15961351
  {
    result = NOT_AROMATIC;
    return 1;
  }

  if (EVERYTHING_HAS_A_PI_ELECTRON == aromaticity_determination_type) {
    result = AROMATIC;
  } else if (Simple_4n_plus_2 == aromaticity_determination_type && 
           allow_two_electron_systems_to_be_aromatic() && 2 == pe) {
    result = AROMATIC;
  } else if (pe > 2 && 2 == pe % 4) {
    result = AROMATIC;
  } else if (pe == 2 && allow_two_electron_systems_to_be_aromatic()) {
    result = AROMATIC;
  } else if (aromaticity_determination_type == ANY_EVEN_NUMBER_OF_PI_ELECTRONS && (pe / 2) * 2 == pe) {
    result = AROMATIC;
  } else {
    result = NOT_AROMATIC;
  }

#ifdef DEBUG_AROMATICITY_DETERMINATION
  cerr << "Pi electron count " << pe << " arom is " << (AROMATIC == result) << '\n';
#endif

  return 1;
}
#endif

/*
  When working with iwfrag, we often set chain atoms to be aromatic.
  We therefore need a means of suppressing the warning which occurs
  when setting a non-ring atom aromatic
*/

static int warn_aromatic_chain_atoms = 1;

void
set_warn_aromatic_chain_atoms(int i)
{
  warn_aromatic_chain_atoms = i;
}

/*
  May 2000. Michal Vieth had some molecules that needed a positive
  charge in order to become aromatic
*/

static int kekule_try_positive_nitrogen = 0;

void
set_kekule_try_positive_nitrogen(int s)
{
  kekule_try_positive_nitrogen = s;
}

static int kekule_try_positive_oxygen = 0;

void
set_kekule_try_positive_oxygen(int s)
{
  kekule_try_positive_oxygen = s;
}

namespace aromatic {
// March 2022. Have been too strict reading aromatic structures.
// Previously I only accepted a Kekule pattern if I thought it was
// aromatic. But that is not correct, I may disagree with the tool
// that wrote the smiles. So, instead, we can accept the first
// valid alternating bond pattern, and then make our own decisions
// about aromaticity.
int any_kekule_bonding_pattern_ok_for_aromatic_input = 1;

void
set_any_kekule_bonding_pattern_ok_for_aromatic_input(int s) {
  any_kekule_bonding_pattern_ok_for_aromatic_input = s;
}

}  // namespace aromatic

/*
  A ring has been found to be aromatic. Update all the bonds
  associated with this ring. We rely on the atoms in the
  Ring object being ordered.
*/

/*int
Molecule::_make_bonds_aromatic (const Ring * r)
{
  int rsize = r->number_elements();

  for (int i = 0; i < rsize; i++)
  {
    atom_number_t a1 = r->item(i);
    Atom * a = _things[a1];

    int acon = a->ncon();
    for (int j = 0; j < acon; j++)
    {
      Bond * b = a->item(j);

      if (b->is_aromatic())    // may have been done here
        continue;

      atom_number_t k = b->other(a1);

      if (r->contains(k))     // the (a1,k) bond is in the ring
        b->set_aromatic();
    }
  }

  return 1;
}*/

int
Molecule::set_aromaticity(atom_number_t zatom, aromaticity_type_t arom,
                          int invalidate_smiles_when_done)
{
  assert(ok_atom_number(zatom));
  assert(OK_ATOM_AROMATICITY(arom));

  if (AROMATIC == arom && is_non_ring_atom(zatom) && warn_aromatic_chain_atoms)
    cerr << "Molecule::set_aromaticity: warning atom " << zatom << " is not in a ring\n";

  if (nullptr == _aromaticity)
    _compute_aromaticity();

  if (arom == _aromaticity[zatom])    // no change
    ;
  else
  {
    _aromaticity[zatom] = arom;
    if (invalidate_smiles_when_done)
      invalidate_smiles();
  }

  return 1;
}

int
Molecule::set_permanent_aromatic(atom_number_t zatom, int s)
{
  _things[zatom]->set_permanent_aromatic(s);

  return 1;
}

int
Molecule::set_bond_permanent_aromatic(atom_number_t a1, atom_number_t a2) {
  for (const Bond* b : *_things[a1]) {
    atom_number_t o = b->other(a1);
    if (o != a2) {
      continue;
    }
    const_cast<Bond*>(b)->set_permanent_aromatic(1);
    return 1;
  }

  return 0;
}

/*
  This function was for Bob Coner, doing work on partial molecule.
  The optional parameter, invalidate_smiles_when_done, can be set to 0,
  which will avoid the call to invalidate_smiles().
  Note that this is specifically designed to put the molecule into
  a strange state - aromatic atoms not in rings for example.
*/

int
Molecule::set_aromaticity_two_atoms(atom_number_t a1, atom_number_t a2,
                                    const aromaticity_type_t arom, int invalidate_smiles_when_done)
{
  assert(are_bonded(a1, a2));

  if (AROMATIC == arom && (is_non_ring_atom(a1) || is_non_ring_atom(a2)) &&
      warn_aromatic_chain_atoms)
    cerr << "Molecule::set_aromaticity_two_atoms: warning atom " << a1 << " or " << a2
         << " is not in a ring\n";

  compute_aromaticity_if_needed();

  int need_to_invalidate_smiles = 0;    // if any changes made

  if (arom != _aromaticity[a1])
  {
    _aromaticity[a1] = arom;
    need_to_invalidate_smiles = 1;
  }

  if (arom != _aromaticity[a2])
  {
    _aromaticity[a2] = arom;
    need_to_invalidate_smiles = 1;
  }

  Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a1, a1));
  if (is_aromatic_atom(arom))
    b->set_aromatic();
  else
    b->set_non_aromatic();

  if (need_to_invalidate_smiles && invalidate_smiles_when_done)
    invalidate_smiles();

  return 1;
}

/*
  Update aromaticity is different from set aromaticity.
  set_aromaticity always does an overwrite.
  _update_aromaticity is for use when discerning aromaticity in a molecule.

  This is unfortunately complicated due the number of special cases.

  A value of AROMATIC overwrites whatever is already stored.
  A value of NOT_AROMATIC will not overwrite a value of AROMATIC - this is
  to accommodate the case of an aromatic ring (already determined) fused to
  a non aromatic ring (just discerned).
*/

//#define DEBUG_UPDATE_AROMATICITY

int
Molecule::_update_aromaticity(const Ring & zring, aromaticity_type_t arom)
{
  int na = zring.number_elements();

#ifdef DEBUG_UPDATE_AROMATICITY
  cerr << "update arom " << zring << " arom = " << arom << '\n';
#endif

  if (AROMATIC == arom)
  {
    if (record_dual_nature_aromaticity)
    {
      for (int i = 0; i < na; i++)
      {
        atom_number_t j = zring[i];
        if (AROMATICITY_NOT_DETERMINED == _aromaticity[j])
          _aromaticity[j] = AROMATIC;
        else
          add_aromatic(_aromaticity[j]);
      }
    }
    else
    {
      for (int i = 0; i < na; i++)
      {
        atom_number_t j = zring[i];
        _aromaticity[j] = AROMATIC;
      }
    }
  }
  else if (NOT_AROMATIC == arom)
  {
    if (record_dual_nature_aromaticity)
    {
      for (int i = 0; i < na; i++)
      {
        atom_number_t j = zring[i];
        if (AROMATICITY_NOT_DETERMINED == _aromaticity[j])
          _aromaticity[j] = NOT_AROMATIC;
        else
          add_aliphatic(_aromaticity[j]);
      }
    }
    else
    {
      for (int i = 0; i < na; i++)
      {
        atom_number_t j = zring[i];
        if (AROMATICITY_NOT_DETERMINED == _aromaticity[j])
          _aromaticity[j] = NOT_AROMATIC;
      }
    }
  }
  else
  {
    cerr << "Strange aromatic value " << arom << endl;
    assert(nullptr == "Huh?");
  }

#ifdef DEBUG_UPDATE_AROMATICITY
  for (int i = 0; i < na; i++)
  {
    atom_number_t j = zring[i];
    cerr << "After update, atom " << j << " arom = " << _aromaticity[j] << endl;
  }
#endif

  return 1;
}

int
Molecule::_update_aromaticity(int ring_number, aromaticity_type_t arom)
{
  assert(OK_ATOM_AROMATICITY(arom));

  if (nullptr == _aromaticity)
    _allocate_aromaticity();

  Ring * r = _sssr_rings[ring_number];
  r->set_aromaticity(arom);

  if (AROMATIC == arom)
    _make_bonds_aromatic(r);

  return _update_aromaticity(*r, arom);
}

//#define DEBUG_COMPUTE_AROMATICITY_FOR_RING

int
Molecule::__compute_aromaticity_for_ring(const Ring & p, aromaticity_type_t & arom,
                                         int & impossible_aromatic,
                                         int& unshared_pi_electrons, AromData& arom_data)
{
#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Determining aromaticity for ring " << p << endl;
#endif

  unshared_pi_electrons = 0;

  if (! _determine_aromaticity(p, impossible_aromatic, arom, arom_data))
    return 0;

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "_determine aromaticity succeeded, impossible = " << impossible_aromatic
       << " arom = " << arom << endl;
#endif

  if (AROMATIC == arom || impossible_aromatic)
    return 1;

  /*            O
          _   ||
           \ /  \
    arom    |    =C-
            |   /
          _/ \ /
              ||
              O

// The ring has been identified as not aromatic.  But, it might be
// combined with something else and become aromatic.  We implement a
// rule that if the unshared atoms of a ring contribute no pi
// electrons, then the ring is not arom.

// I also tried lots of other things here, but none seemed to work as
// well as the above.

// Some things tried:
//   If the extra bonds in the ring were all single bonds, reject it.
//   If there were more extra atoms contributing no pi electrons than atoms
//     contributing an electron.
//   Adjacent atoms each contributing no pi electrons.
*/

  int ring_size = p.number_elements();

#ifdef ONLY_SINGLE_BONDS
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = p[i];
    if (nrings(j) > 1)
      continue;

    int jcon = ncon(j);
    for (int k = 0; k < jcon; k++)
    {
      atom_number_t l;
      bond_type_t bt;
      other_and_type(j, k, l, bt);
      if (p->contains(l) && ! IS_SINGLE_BOND(bt))
        return 1;
    }
  }

  impossible_aromatic = 1;
  return 1;
#else

  // Sept 96. Implement a rule that if there is only one unshared atom
  // and it is an O, S or Se, then the ring is not aromatic.
  // Pearlman's rules make any such ring non-aromatic. Daylight makes
  // it non aromatic if there is one such group, but allows aromatic
  // for more than one occurrence
  // For example:   O1C2=CC=CC=C2OOC2=C1C=CC=C2 PBCHM20184165
  // or             C12=C3C4=C(C=CC=C4)OC1=CC=CC2=CC=C3 PBCHM67456

  // Jun 97, hmmm, seems that isn't correct either. Consider
  // ClC1=C(F)C=C2C(=C1)N1C(=C(C2=O)C(=O)O)SC2=CC=CC=C12 PBCHM13626885
  // which has a five membered ring with an S (or O) atom as the only
  // unshared atom. Daylight 4.51 makes that ring aromatic. Therefore
  // we make the rule for five membered rings only

  int sose_count = 0;    // the number of single ring O, S or Se
  int unshared_atom_count = 0;
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = p[i];

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
    cerr << "Atom " << j << " (" << _things[j]->atomic_symbol() << ") contributes " << arom_data.pi_electrons[j]
         << " pi electrons\n";
#endif

    if (arom_data.pi_electrons[j] < 0)    // perhaps incomplete determination. Skip rest.
      return 1;

    if (0 == arom_data.pi_electrons[j])    // has no pi electrons
      continue;

    if (2 == _things[j]->ncon() || 1 == nrings(j))    // an atom not shared by other rings.
    {
      unshared_atom_count++;
      atomic_number_t z = _things[j]->atomic_number();
      if (8 == z || 16 == z || 34 == z)
        sose_count++;
    }
  }
#endif

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Ring has " << unshared_atom_count << " unshared atoms, sose = " << sose_count << endl;
#endif

  if (1 == unshared_atom_count && 1 == sose_count)
  {
    if (5 == ring_size)
      ;
    else if (7 == ring_size)
      ;
    else if (ring_size < min_aromatic_ring_size)
    {
      impossible_aromatic = 1;
      return 1;
    }
    //  else if (3 == ring_size)     // O1C2=CN=C12 PBCHM20160574. No pubchem does not make these aromatic
    //    ;
    else
    {
      impossible_aromatic = 1;
      return 1;
    }
  }

  if (unshared_atom_count)    // because of the way they are counted above, these are unshared atoms contributing one or more pi electrons
    return 1;

  // None of the unshared atoms contributed any pi electrons.

  arom = NOT_AROMATIC;
  impossible_aromatic = 1;

  return 1;
}

#ifdef TOO_MANY_ARGS
int
Molecule::__compute_aromaticity_for_ring(const Ring & p, aromaticity_type_t & arom,
                                         int & impossible_aromatic, int * e_count,
                                         int & unshared_pi_electrons, int aromaticity_rules)
{
#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Determining aromaticity for ring " << p << endl;
#endif

  unshared_pi_electrons = 0;

  if (! _determine_aromaticity(p, arom, impossible_aromatic, e_count, aromaticity_rules))
    return 0;

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "_determine aromaticity succeeded, impossible = " << impossible_aromatic
       << " arom = " << arom << endl;
#endif

  if (AROMATIC == arom || impossible_aromatic)
    return 1;

  /*            O
          _   ||
           \ /  \
    arom    |    =C-
            |   /
          _/ \ /
              ||
              O

// The ring has been identified as not aromatic.  But, it might be
// combined with something else and become aromatic.  We implement a
// rule that if the unshared atoms of a ring contribute no pi
// electrons, then the ring is not arom.

// I also tried lots of other things here, but none seemed to work as
// well as the above.

// Some things tried:
//   If the extra bonds in the ring were all single bonds, reject it.
//   If there were more extra atoms contributing no pi electrons than atoms
//     contributing an electron.
//   Adjacent atoms each contributing no pi electrons.
*/

  int ring_size = p.number_elements();

#ifdef ONLY_SINGLE_BONDS
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = p[i];
    if (nrings(j) > 1)
      continue;

    int jcon = ncon(j);
    for (int k = 0; k < jcon; k++)
    {
      atom_number_t l;
      bond_type_t bt;
      other_and_type(j, k, l, bt);
      if (p->contains(l) && ! IS_SINGLE_BOND(bt))
        return 1;
    }
  }

  impossible_aromatic = 1;
  return 1;
#else

  // Sept 96. Implement a rule that if there is only one unshared atom
  // and it is an O, S or Se, then the ring is not aromatic.
  // Pearlman's rules make any such ring non-aromatic. Daylight makes
  // it non aromatic if there is one such group, but allows aromatic
  // for more than one occurrence
  // For example:   O1C2=CC=CC=C2OOC2=C1C=CC=C2 PBCHM20184165
  // or             C12=C3C4=C(C=CC=C4)OC1=CC=CC2=CC=C3 PBCHM67456

  // Jun 97, hmmm, seems that isn't correct either. Consider
  // ClC1=C(F)C=C2C(=C1)N1C(=C(C2=O)C(=O)O)SC2=CC=CC=C12 PBCHM13626885
  // which has a five membered ring with an S (or O) atom as the only
  // unshared atom. Daylight 4.51 makes that ring aromatic. Therefore
  // we make the rule for five membered rings only

  int sose_count = 0;    // the number of single ring O, S or Se
  int unshared_atom_count = 0;
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = p[i];

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
    cerr << "Atom " << j << " (" << _things[j]->atomic_symbol() << ") contributes " << e_count[j]
         << " pi electrons\n";
#endif

    if (e_count[j] < 0)    // perhaps incomplete determination. Skip rest.
      return 1;

    if (0 == e_count[j])    // has no pi electrons
      continue;

    if (2 == _things[j]->ncon() || 1 == nrings(j))    // an atom not shared by other rings.
    {
      unshared_atom_count++;
      atomic_number_t z = _things[j]->atomic_number();
      if (8 == z || 16 == z || 34 == z)
        sose_count++;
    }
  }
#endif

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Ring has " << unshared_atom_count << " unshared atoms, sose = " << sose_count << endl;
#endif

  if (1 == unshared_atom_count && 1 == sose_count)
  {
    if (5 == ring_size)
      ;
    else if (7 == ring_size)
      ;
    else if (ring_size < min_aromatic_ring_size)
    {
      impossible_aromatic = 1;
      return 1;
    }
    //  else if (3 == ring_size)     // O1C2=CN=C12 PBCHM20160574. No pubchem does not make these aromatic
    //    ;
    else
    {
      impossible_aromatic = 1;
      return 1;
    }
  }

  if (unshared_atom_count)    // because of the way they are counted above, these are unshared atoms contributing one or more pi electrons
    return 1;

  // None of the unshared atoms contributed any pi electrons.

  arom = NOT_AROMATIC;
  impossible_aromatic = 1;

  return 1;
}
#endif

//#define DEBUG_BUILD_FUSED_SYSTEM

/*
  Recursively build a fused system starting with ring R.
  For each ring included, set the corresponding element in RINGS_IN_FUSED_SYSTEM
*/

static int
build_fused_system(Set_of_Atoms & ring_system, const Ring * r, int system_identifier,
                   int * rings_in_fused_system, const int * ok_to_include,
                   const int * impossible_aromatic)
{
  rings_in_fused_system[r->ring_number()] = system_identifier;
  assert(ok_to_include[r->ring_number()]);

#ifdef DEBUG_BUILD_FUSED_SYSTEM
  cerr << "Build fused system " << system_identifier << " continues with ring " << r->ring_number()
       << ", " << r->fused_ring_neighbours() << " fused neighbours\n";
#endif

  int rc = 1;

  int neighbours = r->fused_ring_neighbours();
  for (int i = 0; i < neighbours; i++)
  {
    const Ring * n = r->fused_neighbour(i);

    int nrn = n->ring_number();
    if (system_identifier == rings_in_fused_system[nrn])    // already detected
      continue;

    if (! ok_to_include[nrn] || n->number_elements() < min_aromatic_ring_size)
      continue;

    //  cerr << "Ring " << (*r) << " and " << (*n) << " share " << r->compute_bonds_shared_with(*n) << " bonds\n";

    if (0 == r->strongly_fused_ring_neighbours())    // no need to check
      ;
    else if (r->compute_bonds_shared_with(*n) > 1)    // no, cannot aromatise ring R
      continue;

    rings_in_fused_system[nrn] = system_identifier;
    if (ring_system.add_non_duplicated_elements(*n))
      rc += build_fused_system(ring_system, n, system_identifier, rings_in_fused_system,
                               ok_to_include, impossible_aromatic);
  }

  return rc;
}

//#define DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS

/*
  Note there is a problem with non-uniqueness here. We should have a canonical means of choosing which
  ring gets added first.

  c12coc([n]1)c1ccccc12 PBCHM70404328
*/

int
Molecule::_determine_aromaticity_by_single_fusions_to_ring(AromData& arom_data,
    int * rings_in_fused_system,
    int * ok_to_include,
    const Ring * r)
{
#ifdef DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS
  cerr << "Trying to expand " << (*r) << endl;
#endif

  if (1 != r->fused_ring_neighbours())
    return 0;

  Set_of_Atoms ring_system;
  ring_system += *r;

  const int system_identifier = r->ring_number() + 1;    // just a number to use, could be anything

  assert(0 == ok_to_include[r->ring_number()]);

  ok_to_include[r->ring_number()] = 1;
  int system_size = build_fused_system(ring_system, r, system_identifier, rings_in_fused_system,
                                       ok_to_include, arom_data.impossible_aromatic);
  ok_to_include[r->ring_number()] = 0;

#ifdef DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS
  cerr << "Found " << system_size << " fused system size, id " << system_identifier << " :";
  const int nr = nrings();
  for (int i = 0; i < nr; ++i)
  {
    cerr << ' ' << rings_in_fused_system[i];
  }
  cerr << endl;
#endif

  if (1 == system_size)    // probably two aliphatic rings together
    return 0;

  aromaticity_type_t arom;
  int tmp;    // no interest in impossible aromatic at this stage

  if (! _determine_aromaticity(ring_system, arom, tmp, arom_data))
    return 0;

#ifdef DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS
  cerr << "Aromaticity found to be " << arom << " arom? " << (AROMATIC == arom) << endl;
#endif

  return (AROMATIC == arom);
}

#ifdef TOO_MANY_ARGS
int
Molecule::_determine_aromaticity_by_single_fusions_to_ring(
    int * pi_electron_count, int aromaticity_rules, int * rings_in_fused_system,
    int * ok_to_include, const int * impossible_aromatic, const Ring * r)
{
#ifdef DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS
  cerr << "Trying to expand " << (*r) << endl;
#endif

  if (1 != r->fused_ring_neighbours())
    return 0;

  Set_of_Atoms ring_system;
  ring_system += *r;

  const int system_identifier = r->ring_number() + 1;    // just a number to use, could be anything

  assert(0 == ok_to_include[r->ring_number()]);

  ok_to_include[r->ring_number()] = 1;
  int system_size = build_fused_system(ring_system, r, system_identifier, rings_in_fused_system,
                                       ok_to_include, impossible_aromatic);
  ok_to_include[r->ring_number()] = 0;

#ifdef DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS
  cerr << "Found " << system_size << " fused system size, id " << system_identifier << " :";
  const int nr = nrings();
  for (int i = 0; i < nr; ++i)
  {
    cerr << ' ' << rings_in_fused_system[i];
  }
  cerr << endl;
#endif

  if (1 == system_size)    // probably two aliphatic rings together
    return 0;

  aromaticity_type_t arom;
  int tmp;    // no interest in impossible aromatic at this stage

  if (! _determine_aromaticity(ring_system, arom, tmp, pi_electron_count, aromaticity_rules))
    return 0;

#ifdef DEBUG_DETERMINE_AROMATICITY_BY_SINGLE_FUSIONS
  cerr << "Aromaticity found to be " << arom << " arom? " << (AROMATIC == arom) << endl;
#endif

  return (AROMATIC == arom);
}
#endif

/*
  This is the last ditch attempt at finding aromaticity.
  We have some still undetermined rings in fused systems.

  Note there is a problem with non-uniqueness here. We should have a canonical means of choosing which
  ring gets added first.
*/

//#define DEBUG_AROM_BY_SINGLE_FUSIONS

int
Molecule::_determine_aromaticity_by_single_fusions(AromData& arom_data)
{
  int rc = 0;

  const int nr = _sssr_rings.number_elements();

  int * rings_in_fused_system = new int[nr];
  std::unique_ptr<int[]> free_rings_in_fused_system(rings_in_fused_system);

#ifdef DEBUG_AROM_BY_SINGLE_FUSIONS
  cerr << "Performing single fusion aromaticity determination\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "i = " << i << " ring already done = " << ring_already_done[i] << endl;
  }
#endif

  for (int i = 0; i < nr; i++)
  {
    if (arom_data.ring_already_done[i] || arom_data.impossible_aromatic[i])
      continue;

    const Ring * ri = _sssr_rings[i];

#ifdef DEBUG_AROM_BY_SINGLE_FUSIONS
    cerr << "Looking for single fusions to :" << ringi(i) << endl;
#endif

    std::fill_n(rings_in_fused_system, nr, 0);

    if (_determine_aromaticity_by_single_fusions_to_ring(arom_data,
                                                         rings_in_fused_system,
                                                         arom_data.ring_already_done, ri))
    {
      for (auto j = 0; j < nr; ++j)
      {
        if (rings_in_fused_system[j])
        {
          _update_aromaticity(j, AROMATIC);
          arom_data.ring_already_done[j] = 1;
          rc++;
        }
      }
//    _update_aromaticity (i, AROMATIC);
//    ring_already_done[i] = 1;
//    rc++;
#ifdef DEBUG_AROM_BY_SINGLE_FUSIONS
      cerr << *ri << " made aromatic by single ring fusion\n";
#endif
    }
  }

  return rc;
}

#ifdef TOO_MANY_ARGS
int
Molecule::_determine_aromaticity_by_single_fusions(int * ring_already_done, int * pi_electron_count,
                                                   const int * impossible_aromatic,
                                                   int aromaticity_rules)
{
  int rc = 0;

  const int nr = _sssr_rings.number_elements();

  int * rings_in_fused_system = new int[nr];
  std::unique_ptr<int[]> free_rings_in_fused_system(rings_in_fused_system);

#ifdef DEBUG_AROM_BY_SINGLE_FUSIONS
  cerr << "Performing single fusion aromaticity determination\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "i = " << i << " ring already done = " << ring_already_done[i] << endl;
  }
#endif

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i] || impossible_aromatic[i])
      continue;

    const Ring * ri = _sssr_rings[i];

#ifdef DEBUG_AROM_BY_SINGLE_FUSIONS
    cerr << "Looking for single fusions to :" << ringi(i) << endl;
#endif

    std::fill_n(rings_in_fused_system, nr, 0);

    if (_determine_aromaticity_by_single_fusions_to_ring(pi_electron_count, aromaticity_rules,
                                                         rings_in_fused_system, ring_already_done,
                                                         impossible_aromatic, ri))
    {
      for (auto j = 0; j < nr; ++j)
      {
        if (rings_in_fused_system[j])
        {
          _update_aromaticity(j, AROMATIC);
          ring_already_done[j] = 1;
          rc++;
        }
      }
//    _update_aromaticity (i, AROMATIC);
//    ring_already_done[i] = 1;
//    rc++;
#ifdef DEBUG_AROM_BY_SINGLE_FUSIONS
      cerr << *ri << " made aromatic by single ring fusion\n";
#endif
    }
  }

  return rc;
}
#endif

//#define DEBUG_COMPUTE_AROMATICITY_FOR_RING

int
Molecule::_compute_aromaticity_for_ring(int ring_number, aromaticity_type_t & arom,
                                        int & impossible_aromatic, AromData& arom_data,
                                        int & unshared_pi_electrons)
{
  const Ring * r = ringi(ring_number);

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Compute arom for ring " << (*r) << endl;
#endif

  if (r->number_elements() > max_aromatic_ring_size)
  {
    arom = NOT_AROMATIC;
    _update_aromaticity(ring_number, arom);
    return 1;
  }

  if (r->number_elements() < min_aromatic_ring_size)
  {
    arom = NOT_AROMATIC;
    _update_aromaticity(ring_number, arom);
    if (EVERYTHING_HAS_A_PI_ELECTRON != arom_data.aromaticity_rule)
      impossible_aromatic = 2;    // hard stop
    return 1;
  }

  if (2 == r->largest_number_of_bonds_shared_with_another_ring() &&
      r->number_elements() < 4)    // heuristic
  {
    arom = NOT_AROMATIC;
    impossible_aromatic = 1;    // may come back later
    return 1;
  }

  int rc = __compute_aromaticity_for_ring(*r, arom, impossible_aromatic, unshared_pi_electrons, arom_data);

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Compute_aromaticity_for_ring " << (*r) << endl;
  cerr << "Arom = " << arom << " impossible = " << impossible_aromatic << " rc " << rc << endl;
#endif

  if (rc)
    _update_aromaticity(ring_number, arom);

  return rc;
}
#ifdef TOO_MANY_ARGS
int
Molecule::_compute_aromaticity_for_ring(int ring_number, aromaticity_type_t & arom,
                                        int & impossible_aromatic, int * pi_electron_count,
                                        int & unshared_pi_electrons, int aromaticity_rules)
{
  const Ring * r = ringi(ring_number);

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Compute arom for ring " << (*r) << endl;
#endif

  if (r->number_elements() > max_aromatic_ring_size)
  {
    arom = NOT_AROMATIC;
    _update_aromaticity(ring_number, arom);
    return 1;
  }

  if (r->number_elements() < min_aromatic_ring_size)
  {
    arom = NOT_AROMATIC;
    _update_aromaticity(ring_number, arom);
    if (EVERYTHING_HAS_A_PI_ELECTRON != aromaticity_rules)
      impossible_aromatic = 2;    // hard stop
    return 1;
  }

  if (2 == r->largest_number_of_bonds_shared_with_another_ring() &&
      r->number_elements() < 4)    // heuristic
  {
    arom = NOT_AROMATIC;
    impossible_aromatic = 1;    // may come back later
    return 1;
  }

  int rc = __compute_aromaticity_for_ring(*r, arom, impossible_aromatic, pi_electron_count,
                                          unshared_pi_electrons, aromaticity_rules);

#ifdef DEBUG_COMPUTE_AROMATICITY_FOR_RING
  cerr << "Compute_aromaticity_for_ring " << (*r) << endl;
  cerr << "Arom = " << arom << " impossible = " << impossible_aromatic << " rc " << rc << endl;
#endif

  if (rc)
    _update_aromaticity(ring_number, arom);

  return rc;
}
#endif

//#define DEBUG_COMBINE_NON_AROM_RING

/*
  Ring ZRING has not been classified. Try to build a contiguous ring system
  from it.
*/

int
Molecule::_combine_non_arom_ring(int zring, int system_identifier, AromData& arom_data,
                                 const int * include_in_ring_systems,
                                 int * rings_in_system)
{
  assert(include_in_ring_systems[zring]);

#ifdef DEBUG_COMBINE_NON_AROM_RING
  cerr << "Combine non arom ring processing ring " << zring << endl;
#endif

  const Ring * r = ringi(zring);

  if (zring != r->ring_number())
  {
    cerr << "Yipes, ring " << zring << " is " << (*r);
  }

  assert(zring == r->ring_number());

  Ring ring_system;
  ring_system.resize(_number_elements);
  ring_system += *r;

  int system_size = build_fused_system(ring_system, r, system_identifier, rings_in_system,
                                       include_in_ring_systems, nullptr);

#ifdef DEBUG_COMBINE_NON_AROM_RING
  cerr << "System size from ring " << zring << " is " << system_size << endl;
  cerr << "Atoms " << ring_system << endl;
#endif

  if (1 == system_size)    // must be fused to a definitely aliphatic ring
    return 0;

  int nr = _sssr_rings.number_elements();

  int tmp;    // we don't worry about the impossible_aromatic return at this stage
  aromaticity_type_t arom;

  int rc = 0;
  if (_determine_aromaticity(ring_system, tmp, arom, arom_data) && arom == AROMATIC)
  {
    for (int i = 0; i < nr; i++)    // update all unprocessed rings in the system
    {
      if (system_identifier == rings_in_system[i])
      {
        rc++;
        _update_aromaticity(i, arom);
        arom_data.ring_already_done[i] = 1;
      }
    }
  }

  return rc;
}


#ifdef TOO_MANY_ARGS
int
Molecule::_combine_non_arom_ring(int zring, int system_identifier, int * ring_already_done,
                                 int * pi_electron_count, const int * include_in_ring_systems,
                                 const int aromaticity_rules, int * rings_in_system)
{
  assert(include_in_ring_systems[zring]);

#ifdef DEBUG_COMBINE_NON_AROM_RING
  cerr << "Combine non arom ring processing ring " << zring << endl;
#endif

  const Ring * r = ringi(zring);

  if (zring != r->ring_number())
  {
    cerr << "Yipes, ring " << zring << " is " << (*r);
  }

  assert(zring == r->ring_number());

  Ring ring_system;
  ring_system.resize(_number_elements);
  ring_system += *r;

  int system_size = build_fused_system(ring_system, r, system_identifier, rings_in_system,
                                       include_in_ring_systems, nullptr);

#ifdef DEBUG_COMBINE_NON_AROM_RING
  cerr << "System size from ring " << zring << " is " << system_size << endl;
  cerr << "Atoms " << ring_system << endl;
#endif

  if (1 == system_size)    // must be fused to a definitely aliphatic ring
    return 0;

  int nr = _sssr_rings.number_elements();

  aromaticity_type_t arom;
  int tmp;    // we don't worry about the impossible_aromatic return at this stage

  int rc = 0;
  if (_determine_aromaticity(ring_system, arom, tmp, pi_electron_count, aromaticity_rules) &&
      AROMATIC == arom)
  {
    for (int i = 0; i < nr; i++)    // update all unprocessed rings in the system
    {
      if (system_identifier == rings_in_system[i])
      {
        rc++;
        _update_aromaticity(i, arom);
        ring_already_done[i] = 1;
      }
    }
  }

  return rc;
}
#endif

/*
  Try to grow systems from unclassified rings to see if they can be
  made into aromatic systems. In this phase of processing, we grow
  the largest possible ring system starting with each unclassified
  ring
*/


#ifdef TOO_MANY_ARGS
int
Molecule::_determine_aromaticity_of_fused_systems(int * ring_already_done, int * pi_electron_count,
                                                  const int * include_in_ring_systems,
                                                  const int aromaticity_rules,
                                                  int & fused_system_identifier)
{
  int nr = _sssr_rings.number_elements();

  int * tmp = new_int(nr);
  std::unique_ptr<int[]> free_tmp(tmp);

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i] || ! include_in_ring_systems[i])
      continue;

    rc += _combine_non_arom_ring(i, fused_system_identifier++, ring_already_done, pi_electron_count,
                                 include_in_ring_systems, aromaticity_rules, tmp);
  }

  return rc;
}
#endif

int
Molecule::_determine_aromaticity_of_fused_systems(AromData& arom_data,
                                                  const int * include_in_ring_systems,
                                                  int & fused_system_identifier)
{
  int nr = _sssr_rings.number_elements();

  int * tmp = new_int(nr);
  std::unique_ptr<int[]> free_tmp(tmp);

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    if (arom_data.ring_already_done[i] || ! include_in_ring_systems[i])
      continue;

    rc += _combine_non_arom_ring(i, fused_system_identifier++, arom_data,
                                 include_in_ring_systems, tmp);
  }

  return rc;
}

/*
  When doing aromaticity determinations, we often need to determine
  whether or not a given ring has any unshared pi electrons.
*/

int
Molecule::_unshared_pi_electrons(const Ring * r, const int * pi_electrons)
{
  int rc = 0;

  int rs = r->number_elements();
  for (int i = 0; i < rs; i++)
  {
    atom_number_t j = r->item(i);

    if (0 == pi_electrons[j])
      continue;

    if (nrings(j) > 1)    // we are only looking for unshared electrons
      continue;

    Atom * a = _things[j];
    int acon = a->ncon();

    //  Is atom J multiply bonded to a heteroatom outside this ring system?
    //  We only do an approximate job here. The detailed determination will
    //  be done by _determine_aromaticity

    if (2 == acon || acon == a->nbonds() || 6 != a->atomic_number())
    {
      rc += pi_electrons[j];
      continue;
    }

    //  Is there a multiple bond to a heteroatom outside the system
    //  We should be able to combine the code from _determine_aromaticity...

    int multiple_bond_to_heteroatom_outside_ring = 0;
    for (int k = 0; k < acon; k++)
    {
      const Bond * b = a->item(k);

      if (b->is_single_bond())
        continue;

      atom_number_t l = b->other(j);
      if (a->atomic_number() == _things[l]->atomic_number())
        continue;

      if (1 == _things[l]->ncon() || is_non_ring_atom(l) || ! r->contains(l))
      {
        multiple_bond_to_heteroatom_outside_ring = 1;
        break;
      }
    }

    if (! multiple_bond_to_heteroatom_outside_ring)
      rc += pi_electrons[j];
  }

  //cerr << "Ring " << r << " has " << rc << " unshared pi electrons\n";
  return rc;
}

//#define DEBUG_INCLUDE_IN_RING_SYSTEMS

/*
  When building a ring system from an unprocessed ring, we need to
  decide which rings will be added to such a system.
  Impossible aromatic's and things with no unshared pi electrons are excluded
*/

int
Molecule::_determine_aromaticity_of_fused_systems(AromData& arom_data)
{
  int nr = _sssr_rings.number_elements();

  int * include_in_ring_systems = new_int(nr, 1); std::unique_ptr<int[]> free_include_in_ring_systems(include_in_ring_systems);

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
  cerr << "_determine_aromaticity_of_fused_systems:checking " << nr << " rings\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << " ring " << i << " done " << ring_already_done[i] << " include "
         << include_in_ring_systems[i] << " impossible " << impossible_aromatic[i] << endl;
  }
#endif

  int callabort = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = ringi(i);
    if (r->ring_number() != i)
    {
      cerr << "Molecule::_determine_aromaticity_of_fused_systems ring not numbered properly, nr = "
           << nr << endl;
      cerr << "Yipes, ring number mismatch, " << i << " vs " << r->ring_number() << endl;
      cerr << "ring " << (*r) << endl;
      callabort = 1;
    }
  }

  if (callabort)
    abort();

  resizable_array<int> rings_with_zero_unshared_pi_electrons;
  for (int i = 0; i < nr; i++)
  {
    if (arom_data.impossible_aromatic[i])
      include_in_ring_systems[i] = 0;
    else if (arom_data.ring_already_done[i])    // aromatic, yes, include it
      ;
    else if (0 == _unshared_pi_electrons(ringi(i), arom_data.pi_electrons))
    {
      include_in_ring_systems[i] = 0;
      rings_with_zero_unshared_pi_electrons.add(i);
    }

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
    cerr << "Ring " << i << ' ' << ringi(i) << " done = " << ring_already_done[i]
         << " include in systems = " << include_in_ring_systems[i] << endl;
#endif
  }

  // _determine_aromaticity_of_fused_systems needs a unique identifier for each system.
  int fused_system_identifier = _number_elements + 1;

  int rc = _determine_aromaticity_of_fused_systems(arom_data,
                                                   include_in_ring_systems,
                                                   fused_system_identifier);

  if (locate_item_in_array(0, nr, arom_data.ring_already_done) < 0)    // all rings done!
    return rc;

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
  cerr << "Rings remain unprocessed\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "ring " << i << " status " << arom_data.ring_already_done[i] << endl;
  }
  cerr << rings_with_zero_unshared_pi_electrons.number_elements()
       << " rings have zero pi electrons\n";
#endif

  // We will now be trying several means of maximising perceived aromaticity

  if (rings_with_zero_unshared_pi_electrons.empty())
    return _determine_aromaticity_by_single_fusions(arom_data);

  for (int i = 0; i < rings_with_zero_unshared_pi_electrons.number_elements(); i++)
  {
    include_in_ring_systems[rings_with_zero_unshared_pi_electrons[i]] = 1;
    //  cerr << "Ring " << rings_with_zero_unshared_pi_electrons[i] << " being included\n";
  }

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
  cerr << "After re-enabling rings with zero unshared pi electrons\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "i = " << i << " include = " << include_in_ring_systems[i] << endl;
  }
#endif

  rc = _determine_aromaticity_of_fused_systems(arom_data, 
                                               include_in_ring_systems,
                                               fused_system_identifier);

  if (locate_item_in_array(0, nr, arom_data.ring_already_done) < 0)
    return rc;

  // We still have one or more aliphatic rings.

  return _determine_aromaticity_by_single_fusions(arom_data);
}

#ifdef TOO_MANY_ARGS
int
Molecule::_determine_aromaticity_of_fused_systems(int * ring_already_done,
                                                  int * pi_electron_count,
                                                  const int * impossible_aromatic,
                                                  const int aromaticity_rules)
{
  int nr = _sssr_rings.number_elements();

  int * include_in_ring_systems = new_int(nr, 1); std::unique_ptr<int[]> free_include_in_ring_systems(include_in_ring_systems);

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
  cerr << "_determine_aromaticity_of_fused_systems:checking " << nr << " rings\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << " ring " << i << " done " << ring_already_done[i] << " include "
         << include_in_ring_systems[i] << " impossible " << impossible_aromatic[i] << endl;
  }
#endif

  int callabort = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = ringi(i);
    if (r->ring_number() != i)
    {
      cerr << "Molecule::_determine_aromaticity_of_fused_systems ring not numbered properly, nr = "
           << nr << endl;
      cerr << "Yipes, ring number mismatch, " << i << " vs " << r->ring_number() << endl;
      cerr << "ring " << (*r) << endl;
      callabort = 1;
    }
  }

  if (callabort)
    abort();

  resizable_array<int> rings_with_zero_unshared_pi_electrons;
  for (int i = 0; i < nr; i++)
  {
    if (impossible_aromatic[i])
      include_in_ring_systems[i] = 0;
    else if (ring_already_done[i])    // aromatic, yes, include it
      ;
    else if (0 == _unshared_pi_electrons(ringi(i), pi_electron_count))
    {
      include_in_ring_systems[i] = 0;
      rings_with_zero_unshared_pi_electrons.add(i);
    }

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
    cerr << "Ring " << i << ' ' << ringi(i) << " done = " << ring_already_done[i]
         << " include in systems = " << include_in_ring_systems[i] << endl;
#endif
  }

  // _determine_aromaticity_of_fused_systems needs a unique identifier for each system.
  int fused_system_identifier = _number_elements + 1;

  int rc = _determine_aromaticity_of_fused_systems(ring_already_done, pi_electron_count,
                                                   include_in_ring_systems, aromaticity_rules,
                                                   fused_system_identifier);

  if (locate_item_in_array(0, nr, ring_already_done) < 0)    // all rings done!
    return rc;

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
  cerr << "Rings remain unprocessed\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "ring " << i << " status " << ring_already_done[i] << endl;
  }
  cerr << rings_with_zero_unshared_pi_electrons.number_elements()
       << " rings have zero pi electrons\n";
#endif

  // We will now be trying several means of maximising perceived aromaticity

  if (rings_with_zero_unshared_pi_electrons.empty())
    return _determine_aromaticity_by_single_fusions(ring_already_done, pi_electron_count,
                                                    impossible_aromatic, aromaticity_rules);

  for (int i = 0; i < rings_with_zero_unshared_pi_electrons.number_elements(); i++)
  {
    include_in_ring_systems[rings_with_zero_unshared_pi_electrons[i]] = 1;
    //  cerr << "Ring " << rings_with_zero_unshared_pi_electrons[i] << " being included\n";
  }

#ifdef DEBUG_INCLUDE_IN_RING_SYSTEMS
  cerr << "After re-enabling rings with zero unshared pi electrons\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "i = " << i << " include = " << include_in_ring_systems[i] << endl;
  }
#endif

  rc = _determine_aromaticity_of_fused_systems(ring_already_done, pi_electron_count,
                                               include_in_ring_systems, aromaticity_rules,
                                               fused_system_identifier);

  if (locate_item_in_array(0, nr, ring_already_done) < 0)
    return rc;

  // We still have one or more aliphatic rings.

  return _determine_aromaticity_by_single_fusions(ring_already_done, pi_electron_count,
                                                  impossible_aromatic, aromaticity_rules);
}
#endif

/*
  Locate COUNT extra rings which can be fused to the ring numbers in RINGS
*/

/*int
identify_fused_rings (Molecule * m,
                      resizable_array<int> & rings,
                      int count,
                      const int * impossible_aromatic)
{
  int rc = 0;
  int rings_in_molecule = m->nrings();

#ifdef DEBUG_IDENTIFY_FUSED_RINGS
  cerr << "Looking for rings fused to " << rings << endl;
#endif

  while (rc < count)
  {
    int found_this_iteration = 0;
    for (int i = 0; i < rings.number_elements(); i++)
    {
      const Ring * ri = m->ringi (rings[i]);
      for (int j = 0; j < rings_in_molecule; j++)
      {
        if (impossible_aromatic[j])
          continue;

        if (rings.contains (j))
          continue;

        const Ring * rj = m->ringi (j);

#ifdef DEBUG_IDENTIFY_FUSED_RINGS
        cerr << "Let's see if ring " << j << " can join\n";
#endif
        if (rj->fused_system_identifier() != ri->fused_system_identifier())
          continue;

        if (ri->is_fused_to (rj))
        {
#ifdef DEBUG_IDENTIFY_FUSED_RINGS
          cerr << "Ring " << * rj << " fused to " << *ri << endl;
#endif
          rings.add (j);
          rc++;
          if (rc >= count)
            return rc;
          found_this_iteration++;
        }
      }
    }
    if (0 == found_this_iteration)
      return rc;
  }

  return rc;
}*/

//#define DEBUG_COMPUTE_AROMATICITY_WITH_FUSED_RINGS

/*
  We have one or more rings which have not been classified as aromatic,
  but which are fused to other rings. See if we can get aromaticity by
  combining them.
*/

/*int
Molecule::_compute_aromaticity_with_fused_rings (int * already_done,
                                const int * impossible_aromatic,
                                int * pi_electron_count,
                                int aromaticity_rules)
{
  int nr = nrings();

#ifdef DEBUG_COMPUTE_AROMATICITY_WITH_FUSED_RINGS
  cerr << "Compute_aromaticity_with_fused_rings\n";
  for (int i = 0; i < nr; i++)
  {
    cerr << "Ring " << i << " already_done = " << already_done[i];
    if (impossible_aromatic[i])
      cerr << " impossible aromatic";
    cerr << endl;
  }
#endif

// The for "count =" loop successively adds rings to ring i, until we get an aromatic combo.
// We should probably take the largest fused set, and then successively remove rings, and
// take the largest aromatic combination(s)
// Will do that some day...

  for (int i = 0; i < nr; i++)
  {
    if (already_done[i])
      continue;

    for (int count = 1; 0 == already_done[i]; count++)
    {
      resizable_array<int> rings;
      rings.add (i);

      if (count > identify_fused_rings (this, rings, count, impossible_aromatic))   // did not find that many rings to fuse.
      {
        already_done[i] = 1;    // give up!
        break;
      }

      Set_of_Atoms super_ring;
      for (int j = 0; j < rings.number_elements(); j++)
      {
        const Ring * rj = ringi (rings[j]);
        super_ring.add_non_duplicated_elements (*rj);
      }

      aromaticity_type_t arom;
      int impossible_aromatic;

      if (_determine_aromaticity (super_ring, arom, impossible_aromatic,
               pi_electron_count, aromaticity_rules) &&
           AROMATIC == arom)
      {
        for (int j = 0; j < rings.number_elements(); j++)
        {
          int k = rings[j];
          _update_aromaticity (k, arom);
          already_done[k] = 1;
        }
      }
    }
  }

  return 1;
}*/

/*
  Main entry point for aromaticity determinations
*/

//#define DEBUG_COMPUTE_AROMATICITY

int
Molecule::_compute_aromaticity(AromData& arom_data)
{
  // Force SSSR determination - or maybe ESSSR

  (void)ring_membership();

  int nr = _sssr_rings.number_elements();

  if (doNotComputeAromaticity)
    return 1;

  int undetermined_rings = nr;

#ifdef DEBUG_COMPUTE_AROMATICITY
  cerr << "Finding aromaticity for " << nr << " rings\n";
#endif

  int rc = 1;    // by default we assume we can do everything

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    aromaticity_type_t arom = AROMATICITY_NOT_DETERMINED;

    if (! strongly_fused_rings_can_be_aromatic &&
        ri->strongly_fused_ring_neighbours())    // no ring strongly fused ring can be aromatic
    {
      arom_data.ring_already_done[i] = 1;
      undetermined_rings--;
      arom_data.impossible_aromatic[i] = 1;
    }
    else if (! _compute_aromaticity_for_ring(i, arom, arom_data.impossible_aromatic[i], arom_data,
                                             arom_data.unshared_pi_electrons[i]))
    {
      arom_data.ring_already_done[i] = 1;    // unable to discern is the same as non-aromatic at this level
      undetermined_rings--;
      rc = 0;    // we must return 0 to indicate a failure (but maybe partial success)
    }
    else if (arom_data.impossible_aromatic[i] || AROMATIC == arom)
    {
//    cerr << "Ring " << i << " atoms " << *ri << " either aromatic " << arom << " or impossible_aromatic " << impossible_aromatic[i] << '\n';
      arom_data.ring_already_done[i] = 1;
      undetermined_rings--;
    }
    else if (ri->number_elements() > max_aromatic_ring_size)
    {
      arom_data.ring_already_done[i] = 1;
      undetermined_rings--;
      arom_data.impossible_aromatic[i] = 1;
    }
    else if (ri->number_elements() < min_aromatic_ring_size)
    {
      arom_data.ring_already_done[i] = 1;
      undetermined_rings--;
      arom_data.impossible_aromatic[i] = 1;
    }
    else if (ri->is_fused())    // let's wait and see if it can be combined into an aromatic system
      ;
    else if (NOT_AROMATIC == arom)    // stand-alone non-aromatic
    {
      arom_data.ring_already_done[i] = 1;
      undetermined_rings--;
    }
    else
    {
      cerr << "Molecule::_compute_aromaticity:Yipes! arom = " << arom
           << " impossible = " << arom_data.impossible_aromatic[i] << endl;
      iwabort();
    }
  }

#ifdef DEBUG_COMPUTE_AROMATICITY
  cerr << "There are " << undetermined_rings << " rings remaining\n";
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = ringi(i);
    cerr << "Ring " << i << ' ' << (*r) << ' ';
    if (r->is_aromatic())
      cerr << "aromatic ";
    if (r->is_non_aromatic())
      cerr << "non aromatic ";
    if (impossible_aromatic[i])
      cerr << "impossible ";
    if (r->undetermined_aromaticity())
      cerr << "undetermined ";
    if (ring_already_done[i])
      cerr << "COMPLETE";
    cerr << endl;
  }
#endif

  if (0 == undetermined_rings)
    return rc;

  // Now this gets ugly. There are one or more fused rings which were classified as
  // non aromatic by themselves. First look at entire fused systems and see if
  // any of them can be classified as aromatic in their entirety

  rc = _determine_aromaticity_of_fused_systems(arom_data);

  return rc;
}

#ifdef TOO_MANY_ARGS
int
Molecule::_compute_aromaticity(int * already_done, int * impossible_aromatic,
                               int * pi_electron_count, int * unshared_pi_electrons,
                               int aromaticity_rules)
{
  // Force SSSR determination - or maybe ESSSR

  (void)ring_membership();

  int nr = _sssr_rings.number_elements();

  if (doNotComputeAromaticity)
    return 1;

  int undetermined_rings = nr;

#ifdef DEBUG_COMPUTE_AROMATICITY
  cerr << "Finding aromaticity for " << nr << " rings\n";
#endif

  int rc = 1;    // by default we assume we can do everything

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    aromaticity_type_t arom = AROMATICITY_NOT_DETERMINED;

    if (! strongly_fused_rings_can_be_aromatic &&
        ri->strongly_fused_ring_neighbours())    // no ring strongly fused ring can be aromatic
    {
      already_done[i] = 1;
      undetermined_rings--;
      impossible_aromatic[i] = 1;
    }
    else if (! _compute_aromaticity_for_ring(i, arom, impossible_aromatic[i], pi_electron_count,
                                             unshared_pi_electrons[i], aromaticity_rules))
    {
      already_done[i] = 1;    // unable to discern is the same as non-aromatic at this level
      undetermined_rings--;
      rc = 0;    // we must return 0 to indicate a failure (but maybe partial success)
    }
    else if (impossible_aromatic[i] || AROMATIC == arom)
    {
//    cerr << "Ring " << i << " atoms " << *ri << " either aromatic " << arom << " or impossible_aromatic " << impossible_aromatic[i] << '\n';
      already_done[i] = 1;
      undetermined_rings--;
    }
    else if (ri->number_elements() > max_aromatic_ring_size)
    {
      already_done[i] = 1;
      undetermined_rings--;
      impossible_aromatic[i] = 1;
    }
    else if (ri->number_elements() < min_aromatic_ring_size)
    {
      already_done[i] = 1;
      undetermined_rings--;
      impossible_aromatic[i] = 1;
    }
    else if (ri->is_fused())    // let's wait and see if it can be combined into an aromatic system
      ;
    else if (NOT_AROMATIC == arom)    // stand-alone non-aromatic
    {
      already_done[i] = 1;
      undetermined_rings--;
    }
    else
    {
      cerr << "Molecule::_compute_aromaticity:Yipes! arom = " << arom
           << " impossible = " << impossible_aromatic[i] << endl;
      iwabort();
    }
  }

#ifdef DEBUG_COMPUTE_AROMATICITY
  cerr << "There are " << undetermined_rings << " rings remaining\n";
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = ringi(i);
    cerr << "Ring " << i << ' ' << (*r) << ' ';
    if (r->is_aromatic())
      cerr << "aromatic ";
    if (r->is_non_aromatic())
      cerr << "non aromatic ";
    if (impossible_aromatic[i])
      cerr << "impossible ";
    if (r->undetermined_aromaticity())
      cerr << "undetermined ";
    if (already_done[i])
      cerr << "COMPLETE";
    cerr << endl;
  }
#endif

  if (0 == undetermined_rings)
    return rc;

  // Now this gets ugly. There are one or more fused rings which were classified as
  // non aromatic by themselves. First look at entire fused systems and see if
  // any of them can be classified as aromatic in their entirety

  rc = _determine_aromaticity_of_fused_systems(already_done, pi_electron_count, impossible_aromatic,
                                               aromaticity_rules);

  return rc;
}
#endif

/*
  Strategy for aromaticity is complex - and maybe wrong (yes)
  We first try each ring as a stand alone ring.
  If all non-fused rings are found to be aromatic, we are done.
  Function _compute_aromaticity sets the impossible flag for any
  ring which cannot be aromatic under any circumstance.
*/

int
Molecule::_compute_aromaticity()
{
  if (0 == _number_elements)
    return 1;

  (void)ring_membership();    // ensure rings computed

  if (nullptr == _aromaticity)
    _allocate_aromaticity();

  int nr = _sssr_rings.number_elements();

  if (nr > 0)
    ;
  else if (aromatic_chain_bonds_are_ok())
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (_things[i]->permanent_aromatic())
        _aromaticity[i] = AROMATIC;
      else
        _aromaticity[i] = NOT_AROMATIC;
    }

    return 1;
  }
  else
  {
    set_vector(_aromaticity, _number_elements, NOT_AROMATIC);

    return 1;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->permanent_aromatic())
      _aromaticity[i] = AROMATIC;
    else if (is_non_ring_atom(i))
      _aromaticity[i] = NOT_AROMATIC;
    else
      _aromaticity[i] = AROMATICITY_NOT_DETERMINED;
  }

#ifdef OLD_VERSION_ASDA
  // We need some temporary arrays, allocate from one

  const int needed = nr + nr + nr + _number_elements;
  int * already_done = new int[needed];
  std::unique_ptr<int[]> free_already_done(already_done);

  std::fill_n(already_done, nr + nr + nr, 0);

  // We need to know which rings cannot be aromatic

  int * impossible_aromatic = already_done + nr;
  int * unshared_pi_electrons = already_done + nr + nr;
  int * pi_electron_count = already_done + nr + nr + nr;

  std::fill_n(pi_electron_count, _number_elements, -99);

  int rc = _compute_aromaticity(already_done, impossible_aromatic, pi_electron_count,
                                unshared_pi_electrons, global_aromaticity_determination_type);
#endif
  AromData arom_data(_number_elements, nr, global_aromaticity_determination_type);
  return _compute_aromaticity(arom_data);
}

/*
  Not sure this is really necessary, why not just make _compute_aromaticity
  public?
*/

int
Molecule::compute_aromaticity()
{
  if (nullptr != _aromaticity)
  {
    delete[] _aromaticity;
    _aromaticity = nullptr;
  }

  for (int i = 0; i < _bond_list.number_elements(); i++)
  {
    _bond_list[i]->set_non_aromatic();
  }

  int nr = _sssr_rings.number_elements();
  for (int i = 0; i < nr; i++)
  {
    _sssr_rings[i]->set_aromaticity_to_not_determined();
  }

  return _compute_aromaticity();
}

int
Molecule::in_same_aromatic_ring(atom_number_t a1, atom_number_t a2)
{
  (void)ring_membership();

  int nr = _sssr_rings.number_elements();

  if (0 == nr)
    return 0;

  aromaticity_type_t arom;
  if (! aromaticity(a1, arom) || (! is_aromatic_atom(arom)))
    return 0;

  if (! aromaticity(a2, arom) || (! is_aromatic_atom(arom)))
    return 0;

  // At this stage, both atoms are aromatic. See if we can find an aromatic
  // ring which contains them both.

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);
    if (! ri->is_aromatic())
      continue;

    if (ri->contains(a1) && ri->contains(a2))
      return 1;
  }

  return 0;
}

/*
  Aug 2000.
  this is incorrect, and no longer needed

int
Molecule::add_aromaticity_to_bonds()
{
  if (0 == nrings())
    return 0;

  if (nullptr == _aromaticity)
    _compute_aromaticity();

  int nb = _bond_list.number_elements();
  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (is_aromatic_atom (_aromaticity[a1]) && is_aromatic_atom (_aromaticity[a2]) &&
        in_same_ring (a1, a2))
    {
      b->set_aromatic();
    }
  }

  return 1;
}*/

int
display_standard_aromaticity_options(std::ostream & os)
{
  os << "  -A <qualifier> Aromaticity, enter \"-A help\" for options\n";

  return 1;
}

int
display_all_aromaticity_options(std::ostream & os)
{
  os << "  -A D           use Daylight aromaticity definitions\n";
  os << "  -A P           use Pearlman aromaticity definitions\n";
  os << "  -A WFL         use Wang Fu Lai aromaticity definitions (modified Pearlman)\n";
  os << "  -A VJ          use Vijay Gombar definitions\n";
  os << "  -A ALLPI       aromatic ring if every atom has pi electrons!\n";
  os << "  -A EVEN        aromatic ring if any even number of pi electrons\n";
  os << "  -A M           enable mixed mode aromaticity (arom and aliph)\n";
  os << "  -A I           enable input of aromatic structures\n";
  os << "  -A O           output aromatic structure files\n";
  os << "  -A K           allow input of structures with no valid kekule form\n";
  os << "  -A 2           allow systems with just two pi electrons to be arom\n";
  os << "  -A C           allow input of delocalised carbonyl bonds (MDL only)\n";
  os << "  -A S           discard obviously non aromatic Kekule input\n";
  os << "  -A N           ignore aromatic atoms/bonds in chains\n";
  os << "  -A okchain     aromatic atoms/bonds in a chain are OK\n";
  os << "  -A tryn+       when an aromatic structure cannot be found, try protonating N atoms\n";
  os << "  -A anyK        when reading aromatic smiles, any Kekule form is accepted\n";
  os << "  -A usmidflt=X  set default aromaticity type for unique smiles\n";
  os << "  -A mnars=N     min size of an aromatic ring, def " << min_aromatic_ring_size << '\n';
  os << "  -A mxars=N     max size of an aromatic ring, def " << max_aromatic_ring_size << '\n';
  os << "  -A perm=a      element 'a' is permanent aromatic type\n";
  os << "  -A oknk        non Kekule rings like c1cccc1 can be aromatic\n";
  os << "  -A nokekule    do not attempt Kekule form determination\n";
  os << "  -A nabar       try to aromatise rings with just some bonds marked aromatic\n";
  os << "  -A arom:       aromatic bonds in smiles written as :\n";
  os << "  -A ipp         when reading smiles, allow aromatic forms from Pipeline Pilot\n";
  os << "  -A sfrna       strongly fused can never be aromatic\n";
  os << "  -A arallsat    a ring can be aromatic even if it contains all fully saturated atoms\n";
  os << "  -A nbh0        no square bracket on an atom means H0\n";

  return os.good();
}

static int _input_aromatic_structures = 1;    // Apr 2010, change default behaviour

void
set_input_aromatic_structures(int i)
{
  _input_aromatic_structures = i;
}

int
input_aromatic_structures()
{
  return _input_aromatic_structures;
}

static int _allow_input_without_valid_kekule_form = 0;

int
allow_input_without_valid_kekule_form()
{
  return _allow_input_without_valid_kekule_form;
}

void
set_allow_input_without_valid_kekule_form(int i)
{
  _allow_input_without_valid_kekule_form = i;

  return;
}

void
set_allow_two_electron_systems_to_be_aromatic(int i)
{
  _allow_two_electron_systems_to_be_aromatic = i;
}

static int _allow_delocalised_carbonyl_bonds = 0;

void
set_allow_delocalised_carbonyl_bonds(int i)
{
  _allow_delocalised_carbonyl_bonds = i;
}

int
allow_delocalised_carbonyl_bonds()
{
  return _allow_delocalised_carbonyl_bonds;
}

/*
  When reading a kekule structure I encountered cases where the input
  specified aromatic, but we would definitely not aromatise a particular
  ring - most commonly because of C=O bonds.
  With this flag, we can remove such atoms from kekule consideration
*/

static int _discard_non_aromatic_kekule_input = 0;

void
set_discard_non_aromatic_kekule_input(int d)
{
  _discard_non_aromatic_kekule_input = d;
}

int
discard_non_aromatic_kekule_input()
{
  return _discard_non_aromatic_kekule_input;
}

// Aug 2023. Realise that letting all Nitrogens vary can result
// in ambiguous, and undesirable Kekule determinations.
// C12=NC3=C(N=C1NN=C2)NN=C3 PBCHM4110374 when aromatized is
// [n]1c2c[n][nH]c2[n]c2[nH][n]cc12.
// But Chemaxon smiles is
// c1n[nH]c2nc3[nH]ncc3nc12
// where 'n' is to be interpreted as [n].
// Traditionally varying hcount was to enable reading smiles from
// sources that had not properly specified hcounts.
static int _h_unspecified_means_zero = 0;

void set_h_unspecified_means_zero(int s) {
  _h_unspecified_means_zero = s;
}

static const Element *
create_permanent_aromatic(const const_IWSubstring & c)
{
  const Element * rc = get_element_from_symbol_no_case_conversion(c);
  if (nullptr == rc)
    rc = create_element_with_symbol(c);

  assert(nullptr != rc);

  const_cast<Element *>(rc)->set_permanent_aromatic(1);

  return rc;
}

/*
  Small function to recognise one of the aromaticity types
*/

static int
string_to_aromaticity_type(const const_IWSubstring & a, IWString & atype)
{
  if ('D' == a)
  {
    atype = "Daylight";
    return Daylight;
  }

  if ('P' == a)
  {
    atype = "Pearlman";
    return Pearlman;
  }

  if ("WFL" == a)
  {
    atype = "Wang Fu Lai";
    return WangFuLai;
  }

  if ("VJ" == a)
  {
    atype = "Vijay_Gombar";
    return Vijay_Gombar;
  }

  if ("ALLPI" == a)
  {
    atype = "All PI";
    return EVERYTHING_HAS_A_PI_ELECTRON;
  }

  if ("EVEN" == a)
  {
    atype = "EVEN";
    file_scope_allow_any_even_number_of_pi_electrons = 1;
    return ANY_EVEN_NUMBER_OF_PI_ELECTRONS;
  }

  if ("PP" == a)
  {
    atype = "PipeLine_Pilot";
    return Pipeline_Pilot;
  }

  if ("PUBCHEM" == a)
  {
    atype = "Pubchem";
    return PUBCHEM_AROMATICITY;
  }

  return 0;
}

int
process_standard_aromaticity_options(Command_Line & cl, int verbose, char aflag)
{
  int i = 0;
  const_IWSubstring c;
  while (cl.value(aflag, c, i++))
  {
    IWString atype;

    int a = string_to_aromaticity_type(c, atype);

    if (0 != a)    // recognised as one of the basic types
    {
      global_aromaticity_determination_type = a;
      if (verbose)
        cerr << "Aromaticity done by " << atype << " rules\n";
    }
    else if ('O' == c)
    {
      set_include_aromaticity_in_smiles(1);
      if (verbose)
        cerr << "Aromaticity included in all structure files written\n";
    }
    else if ('M' == c)
    {
      record_dual_nature_aromaticity = 1;
      if (verbose)
        cerr << "Dual nature aromaticity enabled\n";
    }
    else if ('I' == c)
    {
      set_input_aromatic_structures(1);
      if (verbose)
        cerr << "Aromatic structures can be input\n";
    }
    else if ('K' == c)
    {
      set_input_aromatic_structures(1);
      set_allow_input_without_valid_kekule_form(1);
      if (verbose)
        cerr << "\"Aromatic\" structures with no valid kekule form can be input\n";
    }
    else if ('2' == c)
    {
      set_allow_two_electron_systems_to_be_aromatic(1);
      if (verbose)
        cerr << "\"Aromatic\" rings may have only two pi electrons\n";
    }
    else if ('C' == c)
    {
      set_input_aromatic_structures(1);
      set_allow_delocalised_carbonyl_bonds(1);
      if (verbose)
        cerr << "Delocalised carbonyl bonds can be input (MDL file only)\n";
    }
    else if ('S' == c)
    {
      set_input_aromatic_structures(1);
      set_discard_non_aromatic_kekule_input(1);
      if (verbose)
        cerr << "Non aromatic Kekule input will be discarded\n";
    }
    else if ('N' == c)
    {
      set_input_aromatic_structures(1);
      set_convert_chain_aromatic_bonds(1);
      if (verbose)
        cerr << "Aromatic bonds in chains will be ingored\n";
    }
    else if ("okchain" == c)
    {
      set_input_aromatic_structures(1);
      set_aromatic_chain_bonds_are_ok(1);
      if (verbose)
        cerr << "Aromatic chain atoms allowed!\n";
    }
    else if ("tryn+" == c)
    {
      set_input_aromatic_structures(1);
      kekule_try_positive_nitrogen = 1;

      if (verbose)
        cerr << "Will try putting a positive charge on otherwise failing rings\n";
    }
    else if ("tryo+" == c)
    {
      set_input_aromatic_structures(1);
      kekule_try_positive_oxygen = 1;

      if (verbose)
        cerr << "Will try putting a positive charge on otherwise failing rings (oxygen)\n";
    }
    else if (c == "anyK") {
      aromatic::set_any_kekule_bonding_pattern_ok_for_aromatic_input(1);
      if (verbose) {
        cerr << "First kekule form for aromatic smiles will be used\n";
      }
    }
    else if (c.starts_with("usmidflt="))
    {
      c.remove_leading_chars(9);

      a = string_to_aromaticity_type(c, atype);
      if (0 == a)
      {
        cerr << "Unrecognised default unique smiles aromaticity type '" << c << "'\n";
        return 0;
      }

      set_default_unique_smiles_aromaticity(a);

      if (verbose)
        cerr << "All unique smiles generated with " << atype << " aromaticity rules\n";
    }
    else if (c.starts_with("mnars=")) {
      c.remove_leading_chars(6);
      if (! c.numeric_value(min_aromatic_ring_size) || min_aromatic_ring_size < 3) {
        cerr << "Invalid min aromatic ring size directive '" << c << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Min aromatic ring size " << min_aromatic_ring_size << '\n';
      }
    }
    else if (c.starts_with("mxars=")) {
      c.remove_leading_chars(6);
      if (! c.numeric_value(max_aromatic_ring_size) || max_aromatic_ring_size < 3) {
        cerr << "Invalid max_aromatic ring size directive '" << c << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Max aromatic ring size " << max_aromatic_ring_size << '\n';
      }
    }
    else if (c.starts_with("perm="))
    {
      c.remove_leading_chars(5);

      create_permanent_aromatic(c);
    }
    else if ("oknk" == c)
    {
      set_input_aromatic_structures(1);
      set_non_kekule_systems_ok_to_be_aromatic(1);
      if (verbose)
        cerr << "Non Kekule ring systems can be aromatic\n";
    }
    else if ("nokekule" == c)
    {
      set_input_aromatic_structures(1);
      set_perform_kekule_perception(0);
    }
    else if ("nabar" == c)
    {
      set_all_bonds_in_aromatic_ring_must_be_aromatic(0);
      if (verbose)
        cerr << "Will try to aromatise rings with just some bonds aromatic\n";
    }
    else if ("arom:" == c)
    {
      set_write_smiles_aromatic_bonds_as_colons(1);
      if (verbose)
        cerr << "Aromatic bonds in smiles written as :\n";
    }
    else if ("ipp" == c)
    {
      allow_pipeline_pilot_aromaticity_on_input = 1;
      if (verbose)
        cerr << "When reading, will allow aromaticity types from Pipeline Pilot\n";
    }
    else if ("sfrna" == c)
    {
      strongly_fused_rings_can_be_aromatic = 0;
      if (verbose)
        cerr << "Strongly fused rings cannot be aromatic\n";
    }
    else if ("arallsat" == c)
    {
      file_scope_aromatic_rings_must_contain_unsaturation = 0;
      if (verbose)
        cerr << "Rings with no unsaturation can be aromatic\n";
    }
    else if (c == "nbh0") {
      set_h_unspecified_means_zero(1);
      if (verbose) {
        cerr << "Will interpret 'n' as '[n]' in aromatic smiles\n";
      }
    }
    else if ("help" == c)
    {
      display_all_aromaticity_options(cerr);
      exit(1);    // Note very different behaviour for help!
    }
    else
    {
      cerr << "Unrecognised aromaticity option '" << c << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  During Kekule finding we have lots of times we need various temporary arrays.
  Rather than very long argument lists, we package these things into an object.

  Note that the aromatic_atoms and aromatic_bonds arrays come from somewhere else

  Jan 2005. Ran into the case of 
  N12C(=NC=C1)C1=CC=CC=C1C=C2 PBCHM597765
  as an sd file, where the two 6 membered rings in the fused system had aromatic
  bonds between all the atoms, but the fused 5 membered ring had Kekule single
  and double bonds. The whole system cannot be made aromatic without considering
  the fused ring. For this reason, we introduce the concept of additional fused
  pi electrons.

  Maybe it would simply be better to do all rings in each fused system, but then
  I'd need to do some careful setting of the vary_bonds and such arrays.
*/

class Kekule_Temporary_Arrays
{
 private:
  int _natoms;
  int _nr;

  int _a_atoms_present;

  int * _aromatic_atoms;

  const int * _aromatic_bonds;

  int * _process_these_atoms;

  int * _process_these_rings;

  int * _vary_hcount;

  int * _vary_bonds;    // per atom

  int * _implicit_hydrogens_needed;

  int * _pi_electrons;

  //  Some variables are only set if the molecule is larg//  Some variables are only set if the molecule is large.

  int * _smallest_ring_size;    // for each atom, what is the size of the smallest ring containing it

  int _additional_fused_pi_electrons;

 public:
  Kekule_Temporary_Arrays(int matoms, int nr, int * a, const int * b);
  ~Kekule_Temporary_Arrays();

  void set_a_atoms_found(int s) { _a_atoms_present = s; }
  int a_atoms_found() const { return _a_atoms_present; }

  const int * aromatic_atoms() const { return _aromatic_atoms; }
  const int * aromatic_bonds() const { return _aromatic_bonds; }

  int * aromatic_atoms() { return _aromatic_atoms; }

  int * vary_bonds() { return _vary_bonds; }
  int * vary_hcount() { return _vary_hcount; }

  int * implicit_hydrogens_needed() { return _implicit_hydrogens_needed; }

  int * process_these_atoms() { return _process_these_atoms; }
  int * process_these_rings() { return _process_these_rings; }

  int * pi_electrons() { return _pi_electrons; }

  void set_additional_fused_pi_electrons(int s) { _additional_fused_pi_electrons = s; }
  int additional_fused_pi_electrons() const { return _additional_fused_pi_electrons; }

  int establish_smallest_ring_size(Molecule &);
};

Kekule_Temporary_Arrays::Kekule_Temporary_Arrays(int matoms, int nr, int * a, const int * b)
    : _aromatic_atoms(a), _aromatic_bonds(b)
{
  _natoms = matoms;
  _nr = nr;

  _a_atoms_present = 0;

  _process_these_atoms = new int[matoms];

  copy_vector(_process_these_atoms, _aromatic_atoms, matoms);

  _process_these_rings = new_int(nr, 1);

  _vary_hcount = new int[matoms];
  _vary_bonds = new_int(matoms, 1);

  _implicit_hydrogens_needed = new_int(matoms);

  _pi_electrons = new int[matoms];

  if (_nr > 6)
    _smallest_ring_size = new int[matoms];
  else
    _smallest_ring_size = nullptr;

  _additional_fused_pi_electrons = 0;

  return;
}

/*
  Destructor is complicated by the fact that we don't own the _aromatic_atoms
  array, but we may have made a copy
*/

Kekule_Temporary_Arrays::~Kekule_Temporary_Arrays()
{
  delete[] _process_these_atoms;
  delete[] _process_these_rings;

  delete[] _vary_hcount;
  delete[] _vary_bonds;

  delete[] _implicit_hydrogens_needed;

  delete[] _pi_electrons;

  if (nullptr != _smallest_ring_size)
    delete[] _smallest_ring_size;

  return;
}

int
Kekule_Temporary_Arrays::establish_smallest_ring_size(Molecule & m)
{
  if (nullptr == _smallest_ring_size)
    return 0;

  assert(_nr == m.nrings());

  // the non ring atoms will never be changed, OK
  set_vector(_smallest_ring_size, m.natoms(), m.natoms());

  for (auto i = 0; i < _nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    const int ring_size = ri->number_elements();

    for (int j = 0; j < ring_size; ++j)
    {
      const atom_number_t k = ri->item(j);

      if (ring_size < _smallest_ring_size[k])
        _smallest_ring_size[k] = ring_size;
    }
  }

  return 1;
}

/*
  We have determined that the bonds to a given atom cannot vary
*/

/*int
Molecule::_bonds_to_atom_cannot_vary (Kekule_Temporary_Arrays & kta,
                                      atom_number_t zatom)
{
  int * vary_bonds = kta.vary_bonds();

  vary_bonds[zatom] = 0;

  return 1;
}*/

/*
  Is there just one atom in the aromatic system continuing on from ZATOM,
  but not including APREV
*/

int
Molecule::_identify_continuation_atom(atom_number_t stop_atom, atom_number_t aprev,
                                      atom_number_t zatom, const int * aromatic_atom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  atom_number_t rc = INVALID_ATOM_NUMBER;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (aprev == j)
      continue;

    if (stop_atom == j)
      continue;

    if (aromatic_atom[j] != aromatic_atom[zatom])
      continue;

    if (INVALID_ATOM_NUMBER != rc)    // already found a match, more than one possibility
      return INVALID_ATOM_NUMBER;

    rc = j;
  }

  return rc;
}

//#define DEBUG_MUST_BE_SINGLE_BOND_BETWEEN

/*
  By looking at atoms in the ring, we have figured out that the
  bond between a pair of atoms in an aromatic system must be
  a single bond. Using that information, can we then infer
  some information about the next atom in the ring.

  Because of things like OC1=CC=CC=C1N1N2N1C1=C2C=C(NC(=O)C(=C)C)C=C1 PBCHM59858191
  we need to keep track of one atom on the ring to stop at
*/

int
Molecule::_must_be_single_bond_between(Kekule_Temporary_Arrays & kta, atom_number_t stop_atom,
                                       atom_number_t aprev, atom_number_t zatom)
{
  Atom * a = _things[zatom];

  if (a->ncon() < a->nbonds())
    return 1;

  const int * process_these_atoms = kta.process_these_atoms();

  atom_number_t continuation_atom =
      _identify_continuation_atom(stop_atom, aprev, zatom, process_these_atoms);

#ifdef DEBUG_MUST_BE_SINGLE_BOND_BETWEEN
  cerr << "Molecule::_must_be_single_bond_between:from " << aprev << " to " << zatom
       << " continue-> " << continuation_atom << endl;
#endif

  if (INVALID_ATOM_NUMBER == continuation_atom)
    return 1;

  // Can continuation_atom take a double bond?

  Atom * c = _things[continuation_atom];

  // already got a double bond, cannot have another one - what about n=O groups?
  if (c->ncon() < c->nbonds())
    return 1;

  if (0 == implicit_hydrogens(continuation_atom))
    return 1;

  // Are we dealing with a pyrrole type Nitrogen atom?

  int * vary_bonds = kta.vary_bonds();

  if (6 == a->atomic_number() && 0 == a->formal_charge())
    ;
  else if (7 == a->atomic_number())
  {
    if (3 == a->ncon() && 0 == a->formal_charge())
      ;
    else if (a->implicit_hydrogens_known() && 1 == a->implicit_hydrogens() &&
             0 == a->formal_charge())
      ;
    else
      return 1;

    // cerr << "Atom " << zatom << " no longer varying bonds\n";
    vary_bonds[zatom] = 0;
    return _must_be_single_bond_between(kta, stop_atom, zatom, continuation_atom);
  }
  else
    return 1;

  //cerr << "LINE " << __LINE__ << " setting db btw " << zatom << " " << smarts_equivalent_for_atom(zatom) << " and " << continuation_atom << " " << smarts_equivalent_for_atom(continuation_atom) << endl;

  a->set_bond_type_to_atom(continuation_atom, DOUBLE_BOND);
  _things[continuation_atom]->set_modified();

#ifdef DEBUG_MUST_BE_SINGLE_BOND_BETWEEN
  cerr << "Set double bond between " << zatom << " and " << continuation_atom << endl;
#endif

  if (7 == a->atomic_number() && 2 == a->ncon() && 1 == a->formal_charge() &&
      0 == a->implicit_hydrogens())
    ;
  else {
    // cerr << "Atom " << zatom << " turned off qq\n";
    vary_bonds[zatom] = 0;
  }

  if (7 == c->atomic_number() && 2 == c->ncon() && 1 == c->formal_charge() &&
      0 == c->implicit_hydrogens())
    ;
  else {
    vary_bonds[continuation_atom] = 0;
    // cerr << "Unset continuation_atom " << continuation_atom << '\n';
  }

  // We have now set PREV-ZATOM=CONTINUATION_ATOM
  // Can we go further around the ring

  atom_number_t c2 =
      _identify_continuation_atom(stop_atom, zatom, continuation_atom, process_these_atoms);

  if (INVALID_ATOM_NUMBER == c2)
    return 1;

  return _must_be_single_bond_between(kta, stop_atom, continuation_atom, c2);
}

int
Molecule::_is_pyrrole_type_nitrogen(atom_number_t zatom, atom_number_t & a1, atom_number_t & a2,
                                    atom_number_t & a3, const int * process_these_atoms) const
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;
  a3 = INVALID_ATOM_NUMBER;

  Atom * n = _things[zatom];

  if (7 != n->atomic_number())
    return 0;

  if (0 != n->formal_charge())
    return 0;

  if (n->ncon() < n->nbonds())
    return 0;

  if (3 == n->ncon() && 0 == n->formal_charge())
    ;
  else if (n->implicit_hydrogens_known() && 1 == n->implicit_hydrogens())
    ;
  else
    return 0;

  if (2 == n->ncon())
  {
    a1 = n->other(zatom, 0);
    a2 = n->other(zatom, 1);
    return 1;
  }

  for (int i = 0; i < 3; i++)
  {
    atom_number_t j = n->other(zatom, i);

    if (process_these_atoms[j] != process_these_atoms[zatom])
      continue;

    if (INVALID_ATOM_NUMBER == a1)
      a1 = j;
    else if (INVALID_ATOM_NUMBER == a2)    // 3 connections all in the system
      a2 = j;
    else
      a3 = j;
  }

  return INVALID_ATOM_NUMBER != a2;
}

int
Molecule::_is_furan_or_thiophene(atom_number_t zatom, atom_number_t & a1, atom_number_t & a2) const
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  Atom * o = _things[zatom];

  if (8 == o->atomic_number())
    ;
  else if (16 == o->atomic_number())
    ;
  else
    return 0;

  // Make sure our oplus test case works

  if (0 == o->formal_charge() && kekule_try_positive_nitrogen)
    return 0;

  if (2 != o->ncon())    // rare
    return 0;

  //cerr << "Atom " << o->atomic_symbol() << " has " << o->implicit_hydrogens() << " IH\n";

  if (-1 == o->formal_charge() && 1 == o->implicit_hydrogens())
    ;
  else if (16 == o->atomic_number() && 1 == o->formal_charge() && o->implicit_hydrogens_known() &&
           1 == o->implicit_hydrogens())
    ;
  else if (0 != o->formal_charge())
    return 0;

  a1 = o->other(zatom, 0);
  a2 = o->other(zatom, 1);

  return 1;
}

int
Molecule::_is_nitrogen_double_bond_to_something_outside_ring(atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  atomic_number_t z = a->atomic_number();

  if (6 == z)    // the most common case
    return 0;

  if (7 == z)
    ;
  else if (15 == z)
    ;
  else
    return 0;

  if (3 != a->ncon())
    return 0;

  if (4 != a->nbonds())
    return 0;

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other(zatom);

    if (! in_same_ring(zatom, j))
      return 1;
  }

  return 0;
}

int
Molecule::_is_aromatic_carbonyl(atom_number_t zatom, atom_number_t & a1, atom_number_t & a2) const
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  atom_number_t doubly_bonded_oxygen = INVALID_ATOM_NUMBER;

  const Atom * c = _things[zatom];

  if (6 != c->atomic_number())
    return 0;

  if (3 != c->ncon())
    return 0;

  if (4 != c->nbonds())
    return 0;

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = c->item(i);

    atom_number_t j = b->other(zatom);

    if (b->is_double_bond())
    {
      if (8 == _things[j]->atomic_number())
        doubly_bonded_oxygen = j;
    }
    else if (INVALID_ATOM_NUMBER == a1)
      a1 = j;
    else
      a2 = j;
  }

  return INVALID_ATOM_NUMBER != doubly_bonded_oxygen;
}

//#define DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS

/*
  We need to identify those bonds will be adjusted during the fit.
  Anything with a full valence is excluded, as is anything with
  a double bond (presumably out of the ring, although not necessarily).

  Note too there there is redundancy built into this.
  I will guess, that almost always the first test (checking for
  full valence) will also cover the second case - a double bond
  present. Being paranoid, we do both tests.
*/

int
Molecule::_do_obvious_bond_order_settings(Kekule_Temporary_Arrays & kta)
{
  const int * process_these_atoms = kta.process_these_atoms();

  int * vary_bonds = kta.vary_bonds();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    Atom * a = _things[i];

    const Element * e = a->element();

    if (0 == e->normal_valence())
    {
      cerr << "Molecule::_do_obvious_bond_order_settings: atom " << i << " (" << e->symbol()
           << ") has zero valence\n";
      return 0;
    }

#ifdef DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS
    cerr << "Atom " << i << " (" << atomic_symbol(i) << ") has " << a->nbonds()
         << " bonds, ih = " << a->implicit_hydrogens();
    if (a->implicit_hydrogens_known())
      cerr << ", H fixed";
    ;
    if (a->formal_charge())
      cerr << ", formal charge " << a->formal_charge();

    cerr << " VB? " << vary_bonds[i] << endl;
#endif

    atomic_number_t z = a->atomic_number();
    int icon = a->ncon();
    int ibonds = a->nbonds();

    //  cerr << "Atom " << i << ' ' << e->symbol() << " ncon " << icon << " ibonds " << ibonds << " charge " << a->formal_charge() << " ih = " << a->implicit_hydrogens() << endl;

    if (7 == z && 3 == icon && 4 == ibonds && 0 == a->formal_charge())
      ;
    else if (2 == a->ncon() && (8 == z || 16 == z))
    {
      if (0 == a->formal_charge() && 0 == kekule_try_positive_nitrogen)
      {
        vary_bonds[i] = 0;
        continue;
      }

      if (-1 == a->formal_charge())    // not sure this really exists..
      {
        vary_bonds[i] = 0;
        continue;
      }

      //    Beware S+. If we ask for ->implicit_hydrogens() we may get 1, but likely it has 0

      int ih;
      if (a->implicit_hydrogens_known())
        ih = a->implicit_hydrogens();
      else
        ih = 0;

      if (ibonds + ih == e->normal_valence() + a->formal_charge())
      {
        vary_bonds[i] = 0;
        continue;
      }
    }
    else if (a->atomic_number() == 7 && a->formal_charge() == -1 && a->ncon() == 2 && a->implicit_hydrogens() == 0) {
      // [nH]1ccc[n-]c1=O
    }
    else if (ibonds + a->implicit_hydrogens() == e->normal_valence() + a->formal_charge())
    {
#ifdef DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS
      cerr << "vary_bonds for atom " << i << " set to zero, ibonds " << ibonds << " ih "
           << a->implicit_hydrogens() << endl;
#endif
      // cerr << "Unset atom " << i << " here\n";
      vary_bonds[i] = 0;
      continue;
    }

    if (6 == z && -1 == a->formal_charge() && ((3 == a->ncon() || 1 == a->implicit_hydrogens())))
    {
      vary_bonds[i] = 0;
      continue;
    }

#ifdef DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS
    cerr << "Check " << icon << " connections\n";
#endif

    for (int j = 0; j < icon; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_triple_bond())
      {
        if (file_scope_display_no_kekule_form_message)
          cerr << "Molecule::_do_obvious_bond_order_settings: atom " << i << " has a triple bond\n";
        return 0;
      }

      //    To ensure max compatibility with Daylight, allow n(=O) to vary
      //    12 Sept 96. Ran into this case,  C=N(=C)C which had been aromatisised.
      //    Therefore we generalise the rule to any kind of 4 bonded Nitrogen.

      if (b->is_double_bond())
      {
        if (7 == a->atomic_number() && 3 == icon && 4 == ibonds && 0 == a->formal_charge())
          ;
        else if (15 == a->atomic_number() && 3 == icon && 4 == ibonds && 0 == a->formal_charge())
          ;
        else if (16 == a->atomic_number() && 4 == icon && 5 == ibonds &&
                 allow_pipeline_pilot_aromaticity_on_input && 0 == a->formal_charge())
          ;
        else if (16 == a->atomic_number() && 3 == icon && 4 == ibonds &&
                 allow_pipeline_pilot_aromaticity_on_input && 0 == a->formal_charge())
          ;
        else
          vary_bonds[i] = 0;
      }
    }
  }

#ifdef DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS
  cerr << "Midway through _do_obvious_bond_order_settings\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " i = " << i << " VB? " << vary_bonds[i] << endl;
  }
#endif

  // Oct 2004. Whenever you have an O=c group, we can get some information
  // about the bonding of the atoms ortho to that

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    atom_number_t a1, a2, a3;

    a3 = INVALID_ATOM_NUMBER;

    if (_is_furan_or_thiophene(i, a1, a2))
      ;
    else if (_is_pyrrole_type_nitrogen(i, a1, a2, a3, process_these_atoms))
      ;
    else if (_is_aromatic_carbonyl(i, a1, a2))
      ;
    else
      continue;

#ifdef DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS
    cerr << "Atom " << i << " is an aromatic carbonyl, a1=" << a1 << " a2=" << a2;
    if (INVALID_ATOM_NUMBER != a3)
      cerr << " a3 = " << a3;
    cerr << endl;
#endif

    _must_be_single_bond_between(kta, i, i, a1);
    _must_be_single_bond_between(kta, i, i, a2);
    if (INVALID_ATOM_NUMBER != a3)
      _must_be_single_bond_between(kta, i, i, a3);
  }

#ifdef DEBUG_DO_OBVIOUS_BOND_ORDER_SETTINGS
  cerr << "At end of _do_obvious_bond_order_settings\n";
  for (int i = 0; i < _number_elements; i++)
  {
    if (process_these_atoms[i])
      cerr << " atom " << i << " vary_bonds " << vary_bonds[i] << endl;
  }
#endif

  return 1;
}

// A system may contain both explicit aromatic atoms, and rings
// represented in their Kekule forms: c1cnc2n1C=CC3=C2C=CC=C3
// Count the number of pi electrons in rings that do not touch
// the aromatic atoms being processed. For that molecuile, there are
// 4 extra pi electrons.
int
Molecule::_count_pi_electrons_in_any_fused_but_kekule_form_rings(
    Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings)
{
  kta.set_additional_fused_pi_electrons(0);

  const int fsid = rings[0]->fused_system_identifier();

//cerr << "fsid " << fsid << " rngs " << rings.size() << '\n';
  const int nr = _sssr_rings.number_elements();

  if (nr == rings.number_elements())    // all rings were designated aromatic
    return 1;

  int extra_pi_electrons = 0;

  const int * process_these_atoms = kta.process_these_atoms();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    if (ri->fused_system_identifier() != fsid)
      continue;

    if (rings.contains(ri))
      continue;

    // maybe a ring system was divided into disjoint pieces
    if (! ri->any_members_set_in_array(process_these_atoms))
      continue;

    for (int j = 0; j < ri->number_elements(); j++)
    {
      atom_number_t k = ri->item(j);
      if (process_these_atoms[k])
        continue;

      int pi;
      (void)_things[k]->pi_electrons(pi);
      extra_pi_electrons += pi;
    }
  }

  kta.set_additional_fused_pi_electrons(extra_pi_electrons);

  return 1;
}

/*
  Somewhere in this molecule, there are some 'a' atoms. Are there any
  in this set of rings
*/

int
Molecule::_a_atoms_found_this_set_of_rings(const Kekule_Temporary_Arrays & kta,
                                           const resizable_array<const Ring *> & rings)
{
  int nr = rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = rings[i];

    int ring_size = ri->number_elements();

    for (int j = 0; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);

      const Atom * ak = _things[k];

      if (ak->element()->organic())    // the most common case
        continue;

      if ("a" == ak->element()->symbol())
        return 1;
    }
  }

  return 0;
}

AromData::AromData(int natoms, int nrings, int arom) : aromaticity_rule(arom) {
  ring_already_done = new_int(nrings + nrings + nrings);
  impossible_aromatic = ring_already_done + nrings;
  unshared_pi_electrons = ring_already_done + nrings + nrings;

  pi_electrons = new int[natoms];
  std::fill_n(pi_electrons, natoms, -99);

  in_all_pi_ring = nullptr;
}

AromData::~AromData() {
  delete [] ring_already_done;
  delete [] pi_electrons;

  if (in_all_pi_ring != nullptr) {
    delete [] in_all_pi_ring;
  }
}

/*
  We have a ring or ring system in which there are 'a' atoms. Make sure
  all atoms within the system are permanent_aromatic, and all bonds
  between member atoms are permanent_aromatic
*/

int
Molecule::_kekule_all_atoms_and_bonds_in_system_aromatic(
    Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings)
{
  //cerr << "Calling _kekule_all_atoms_and_bonds_in_system_aromatic\n";

  const int * aromatic_atoms = kta.aromatic_atoms();

  int nr = rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = rings[i];

    int ring_size = ri->number_elements();

    for (int j = 0; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);

      Atom * ak = _things[k];

      ak->set_permanent_aromatic(1);
      //    cerr << "Atom " << k << " becomes permanent aromatic\n";

      int kcon = ak->ncon();

      for (int l = 0; l < kcon; l++)
      {
        Bond * b = ak->item(l);

        atom_number_t m = b->other(k);

        // bond either outside the ring, or to a different aromatic system (biphenyl)
        if (aromatic_atoms[k] != aromatic_atoms[m])
          continue;

        if (b->is_aromatic())
          ;
        else if (! b->is_single_bond())
        {
          // cerr << "Molecule::_kekule_all_atoms_and_bonds_in_system_aromatic:bond not aromatic\n";
          // cerr << "'" << name() << "' atoms " << b->a1() << " to " << b->a2() << endl;
        }

        b->set_permanent_aromatic(1);
      }
    }
  }

  return 1;
}

//#define DEBUG_KEKULE

/*
  Make any obvious hcount adjustments. For those atoms which are
  not obvious, set VARY_HCOUNT to 1, which means that as aromatic
  systems are determined, the hcount will be allowed to vary.
*/

int
Molecule::_do_obvious_hcount_adjustments(Kekule_Temporary_Arrays & kta)
{
  const int * process_these_atoms = kta.process_these_atoms();
  int * vary_hcount = kta.vary_hcount();
  const int * implicit_hydrogens_needed = kta.implicit_hydrogens_needed();

#ifdef DEBUG_KEKULE
  cerr << "Doing obvious vary_hcount stuff\n";
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    Atom * a = _things[i];

    int icon = a->ncon();

    const atomic_number_t z = a->atomic_number();

    if (2 == icon)
      ;
    else if (16 == z && icon > 2 && allow_pipeline_pilot_aromaticity_on_input &&
             0 == a->formal_charge())
    {
      if (4 == icon)
      {
        vary_hcount[i] = 0;
        a->set_implicit_hydrogens(0, 1);
        continue;
      }
      else if (3 == icon && 4 == a->nbonds())
      {
        vary_hcount[i] = 0;
        a->set_implicit_hydrogens(1, 1);
        continue;
      }
    }
    else if (3 == icon)
      ;
    else if (15 == z && 4 == icon && allow_pipeline_pilot_aromaticity_on_input)
    {
      vary_hcount[i] = 0;
      a->set_implicit_hydrogens(0, 1);
      a->set_implicit_hydrogens_known(1);
    }
    else
    {
      if (file_scope_display_no_kekule_form_message)
        cerr << "Atom " << i << " has " << icon << " connections but is supposed to be aromatic\n";
      return 0;
    }

#ifdef DEBUG_KEKULE
    cerr << "Atom " << i << " (" << a->atomic_symbol() << "), ncon " << a->ncon() << " needs "
         << implicit_hydrogens_needed[i] << " ih, compute " << a->implicit_hydrogens() << endl;
#endif

    if (implicit_hydrogens_needed[i] >= 0) {
      vary_hcount[i] = 0;
      if (a->implicit_hydrogens() != implicit_hydrogens_needed[i])
        a->set_implicit_hydrogens(implicit_hydrogens_needed[i], 1);
      continue;
    }

    // don't check valences. Ran into problems with [n] in a ring. That isn't a valid N state just now
    if (a->implicit_hydrogens_known())
      ;
    else if (
        7 == z && 3 == a->ncon() && 0 == a->formal_charge() &&
        4 == a->nbonds())    // one of those wierd ND3v4 that hopefully will acquire a positive charge
      ;
    else if (! a->valence_ok())
    {
      if (display_abnormal_valence_messages())
        cerr << "Atom " << i << " (" << atomic_symbol(i) << ") ncon = " << icon
             << " has an abnormal valence\n";
      return 0;
    }

    vary_hcount[i] = 0;

    // known value from input, don't change
    if (a->implicit_hydrogens_known()) {
      continue;
    }

    int ih;
    if (! a->compute_implicit_hydrogens(ih))    // how could this happen
    {
      cerr << "Yipes, cannot compute implicit hydrogens for atom " << i << endl;
      continue;
    }

    if (0 == ih)    // don't mess with it
      continue;

    if (3 == icon)
      a->set_implicit_hydrogens(0);
    else if (6 == z)
      a->set_implicit_hydrogens(1);
    else if (z == 7) {
      if (_h_unspecified_means_zero) {
        a->set_implicit_hydrogens(ih);
        vary_hcount[i] = 0;
      } else if (ih == 1) {
        a->set_implicit_hydrogens(0);
        vary_hcount[i] = 1;
      }
    }
    else if (8 == z)
      a->set_implicit_hydrogens(0);    // especially to deal with O+
    // P can be aromatic - for now I'm not dealing with the case of [PH] as an aromatic
    else if (15 == z)
    {
      if (2 == a->ncon())
        a->set_implicit_hydrogens(0);
      vary_hcount[i] = 0;
    }
    else if (16 == z)
    {
      if (1 == a->formal_charge())
        vary_hcount[i] = 1;

      a->set_implicit_hydrogens(0);
    }
    else if (2 == ih)
    {
      vary_hcount[i] = 1;
      a->set_implicit_hydrogens(0);    // may change it later
    }
  }

#ifdef DEBUG_KEKULE
  cerr << "After making obvious hcount adjustments\n";
  for (int i = 0; i < _number_elements; i++)
  {
    if (process_these_atoms[i])
    {
      cerr << "Atom " << i << " (" << atomic_symbol(i) << ") ncon = " << ncon(i)
           << " hcount = " << implicit_hydrogens(i);
      if (formal_charge(i))
        cerr << " fc = " << formal_charge(i);

      cerr << " vary hcount = " << vary_hcount[i] << endl;
    }
  }
#endif

  return 1;
}

int
Molecule::__kekule_arom_test_rings(const int * process_these_atoms,
                                   const resizable_array<const Ring *> & rings,
                                   aromaticity_type_t * ring_aromaticity,
                                   aromaticity_type_t * atom_aromaticity, 
                                   AromData& arom_data)
{

#ifdef DEBUG_KEKULE
  cerr << "Finding kekule form for " << rings.number_elements() << " rings\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " arom " << atom_aromaticity[i] << endl;
  }
#endif

  // First compute the pi electrons for the atoms being processed

  int nr = rings.number_elements();

  // Even though _compute_aromaticity_for_ring may fail,

  for (int i = 0; i < nr; i++)
  {
    aromaticity_type_t arom;
    int impossible_aromatic;
    const Ring * r = rings[i];
    if (! __compute_aromaticity_for_ring(*r, arom, impossible_aromatic, arom_data.unshared_pi_electrons[i], arom_data))
      return 0;

    if (impossible_aromatic)
      continue;

    //  Even though this ring is non aromatic, perhaps all its atoms will be -
    //  requires every atom in R to be fused to an aromatic ring.

    ring_aromaticity[i] = arom;
    //  cerr << "Ring " << *r << " is " << arom << endl;
    r->set_vector(atom_aromaticity, arom);
  }

#ifdef DEBUG_KEKULE
  cerr << "Computed initial aromaticity for all " << nr << " rings\n";
  for (int i = 0; i < _number_elements; i++)
  {
    if (! process_these_atoms[i])
      continue;

    cerr << "Atom " << i;

    if (AROMATICITY_NOT_DETERMINED == atom_aromaticity[i])
      cerr << " not determined\n";
    else if (AROMATIC == atom_aromaticity[i])
      cerr << " aromatic\n";
    else if (NOT_AROMATIC == atom_aromaticity[i])
      cerr << " aliphatic\n";
    else
      cerr << " huh?? " << atom_aromaticity[i] << endl;
  }
#endif

  // Even if not all rings were aromatic, we may still have classified every atom
  // as aromatic - think of the 4 membered ring between two phenyl rings

  int unclassified_atoms = 0;
  int misclassified_atoms = 0;
  for (int i = 0; i < _number_elements && 0 == unclassified_atoms; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (AROMATICITY_NOT_DETERMINED == atom_aromaticity[i])
      unclassified_atoms++;
    else if (NOT_AROMATIC == atom_aromaticity[i])
      misclassified_atoms++;
  }

#ifdef DEBUG_KEKULE
  cerr << misclassified_atoms << " misclassified_atoms\n";
#endif

  if (misclassified_atoms)
    return 0;

  if (0 == unclassified_atoms)
    return 1;

#ifdef DEBUG_KEKULE
  cerr << "Will call fused rings determination\n";
#endif

  // We must now call _determine_aromaticity_of_fused_systems.
  // It needs a vector already_done, or length nrings(), which describes
  // which rings have already been classified.

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = rings[i];
    arom_data.ring_already_done[r->ring_number()] = 0;    // process ring J
  }

  int rc = _determine_aromaticity_of_fused_systems(arom_data);

  // _compute_aromaticity_with_fused_rings updates the aromaticity information
  // in the molecule. For now, we are leaving it in place, but watch out for
  // problems later

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = rings[i];
    if (! r->is_aromatic())
      return 0;
  }

  return rc;
}

#ifdef TOO_MANY_ARGS
int
Molecule::__kekule_arom_test_rings(const int * process_these_atoms,
                                   const resizable_array<const Ring *> & rings,
                                   aromaticity_type_t * ring_aromaticity,
                                   aromaticity_type_t * atom_aromaticity, int * pi_e,
                                   int * unshared_pi_electrons)
{

#ifdef DEBUG_KEKULE
  cerr << "Finding kekule form for " << rings.number_elements() << " rings\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " arom " << atom_aromaticity[i] << endl;
  }
#endif

  // First compute the pi electrons for the atoms being processed

  int nr = rings.number_elements();

  // Even though _compute_aromaticity_for_ring may fail,

  for (int i = 0; i < nr; i++)
  {
    aromaticity_type_t arom;
    int impossible_aromatic;
    const Ring * r = rings[i];
    if (! __compute_aromaticity_for_ring(*r, arom, impossible_aromatic, pi_e,
                                         unshared_pi_electrons[i], Simple_4n_plus_2))
      return 0;

    if (impossible_aromatic)
      continue;

    //  Even though this ring is non aromatic, perhaps all its atoms will be -
    //  requires every atom in R to be fused to an aromatic ring.

    ring_aromaticity[i] = arom;
    //  cerr << "Ring " << *r << " is " << arom << endl;
    r->set_vector(atom_aromaticity, arom);
  }

#ifdef DEBUG_KEKULE
  cerr << "Computed initial aromaticity for all " << nr << " rings\n";
  for (int i = 0; i < _number_elements; i++)
  {
    if (! process_these_atoms[i])
      continue;

    cerr << "Atom " << i;

    if (AROMATICITY_NOT_DETERMINED == atom_aromaticity[i])
      cerr << " not determined\n";
    else if (AROMATIC == atom_aromaticity[i])
      cerr << " aromatic\n";
    else if (NOT_AROMATIC == atom_aromaticity[i])
      cerr << " aliphatic\n";
    else
      cerr << " huh?? " << atom_aromaticity[i] << endl;
  }
#endif

  // Even if not all rings were aromatic, we may still have classified every atom
  // as aromatic - think of the 4 membered ring between two phenyl rings

  int unclassified_atoms = 0;
  int misclassified_atoms = 0;
  for (int i = 0; i < _number_elements && 0 == unclassified_atoms; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (AROMATICITY_NOT_DETERMINED == atom_aromaticity[i])
      unclassified_atoms++;
    else if (NOT_AROMATIC == atom_aromaticity[i])
      misclassified_atoms++;
  }

#ifdef DEBUG_KEKULE
  cerr << misclassified_atoms << " misclassified_atoms\n";
#endif

  if (misclassified_atoms)
    return 0;

  if (0 == unclassified_atoms)
    return 1;

#ifdef DEBUG_KEKULE
  cerr << "Will call fused rings determination\n";
#endif

  // We must now call _determine_aromaticity_of_fused_systems.
  // It needs a vector already_done, or length nrings(), which describes
  // which rings have already been classified.

  int rings_in_molecule = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = rings[i];
    arom_data.ring_already_done[r->ring_number()] = 0;    // process ring J
  }

//int rc = _determine_aromaticity_of_fused_systems(already_done, pi_e, impossible_aromatic_not_used,
//                                                 Simple_4n_plus_2);

  int rc = _determine_aromaticity_of_fused_systems(arom_data);

  // _compute_aromaticity_with_fused_rings updates the aromaticity information
  // in the molecule. For now, we are leaving it in place, but watch out for
  // problems later

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = rings[i];
    if (! r->is_aromatic())
      return 0;
  }

  return rc;
}
#endif

/*
  Looks as if there is some kind of fused system.
  We set the aromaticity determination type to simple 4n+2
  This function is just an interface to __kekule_arom_test_rings
*/

int
Molecule::_kekule_arom_test_rings(Kekule_Temporary_Arrays & kta,
                                  const resizable_array<const Ring *> & rings)
{
  const int * process_these_atoms = kta.process_these_atoms();

  int nr = rings.number_elements();

  aromaticity_type_t * ra = new aromaticity_type_t[nr];
  std::unique_ptr<aromaticity_type_t[]> free_ra(ra);

  set_vector(ra, nr, AROMATICITY_NOT_DETERMINED);

  aromaticity_type_t * aa = new aromaticity_type_t[_number_elements];
  std::unique_ptr<aromaticity_type_t[]> free_aa(aa);

  set_vector(aa, _number_elements, AROMATICITY_NOT_DETERMINED);

  int * pi_e = kta.pi_electrons();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    int tmp;
    if (! _things[i]->pi_electrons(tmp))
      return 0;

    pi_e[i] = tmp;
  }

  AromData arom_data(_number_elements, _sssr_rings.number_elements(), Simple_4n_plus_2);
  std::copy_n(pi_e, _number_elements, arom_data.pi_electrons);
  return __kekule_arom_test_rings(process_these_atoms, rings, ra, aa, arom_data);


#ifdef TOO_MANY_ARGS
  int * unshared_pi_electrons = new_int(nr);
  std::unique_ptr<int[]> free_unshared_pi_electrons(unshared_pi_electrons);

  return __kekule_arom_test_rings(process_these_atoms, rings, ra, aa, pi_e, unshared_pi_electrons);
#endif
}

#define KEKULE_READY_TO_PROCESS -1
#define KEKULE_PROCESSED 1

/*
  This code lifted from the code in the aromaticity determination.
  Should be combined.
*/

int
Molecule::_bonded_to_heteroatom_outside_system(atom_number_t zatom, const int * aromatic_system)
{
  Atom * a = _things[zatom];

  int acon = a->ncon();
  //cerr << "Atom " << zatom << " ncon " << acon << " ih " << implicit_hydrogens(zatom) << " bonds " << a->nbonds() << endl;
  if (acon + implicit_hydrogens(zatom) == a->nbonds())
    return 0;

#ifdef DEBUG_KEKULE
  cerr << "Checking for multiple bonds outside the ring\n";
#endif

  const atomic_number_t za = a->atomic_number();

  // For maximum compatibility with Weininger, this only applies to C atoms
  // in the ring. This is because I found lots of n(=O) systems in which
  // it was clear that he had allowed one electron in the 'n'

  // Mar 2004. but it looks like the S atoms contributes electrons O=C1NS(=O)NC2=CC=CC=C12 PBCHM71359875

  if (16 == za)
    ;
  else if (7 == za && 1 == a->formal_charge())
    ;
  else if (6 != za)    // do not consider heteroatoms
    return 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    if (b->is_single_bond())
      continue;

    atom_number_t j = b->other(zatom);

    //  The atom types must differ

    if (za == _things[j]->atomic_number())
      continue;

    //  Either
    //    atom j is non ring
    //    atom j is not in this ring, nor are atoms j and l in any ring together.

    // A terminal atom is definitely outside the ring, likely most common case.
    if (_things[j]->ncon() == 1) {
      return 1;
    }

    // Jan 2022. Make sure we handle C1Cn2cn[nH]c2=N1
    // if (0 == aromatic_system[j] && ! in_same_ring(zatom, j))
    //   return 1;
    // Should also check in all pi system, but AromData not avaiable. monitor...
#ifdef DEBUG_KEKULE
    cerr << "Checking _amide_like for atoms " << zatom << " and " << j  << " _amide_like " << _amide_like(zatom, j) << '\n';
#endif
    if (aromatic_system[j] == 0) {
      return 1;
    }
    if (_amide_like(zatom, j)) {
      return 0;
    }
  }

  return 0;
}

/*
  A branch has been exhausted. Can we identify another unprocessed
  atom. It must be bonded to one already processed. Be sure th check
  the error case of unprocessed atoms, but none connected to anything
  processed
*/

int
Molecule::_kekule_identify_next_atom(const int * process_these, atom_number_t & zatom) const
{
  int found = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (KEKULE_READY_TO_PROCESS != process_these[i])
      continue;

    found++;

    //  We have an unprocessed atom. Is it bonded to an already processed atom?

    const Atom * a = _things[i];

    int icon = a->ncon();
    for (int j = 0; j < icon; j++)
    {
      atom_number_t k = a->other(i, j);
      if (KEKULE_PROCESSED == process_these[k])
      {
        zatom = i;
        return 1;
      }
    }
  }

  // We did not find any suitable atoms.

  //assert (0 == found);   sep 06, not sure why this was here

  return 0;
}

/*
  One of the Kekule finders has found a system with PI_ELECTRONS pi electrons.
  Does it look OK?
*/

static int
looks_like_valid_kekule_system(const Kekule_Temporary_Arrays & kta, int pi_electron_count)
{
  //cerr << "Can " << pi_electron_count << " be aromatic?\n";

  if (allow_two_electron_systems_to_be_aromatic() && 2 == pi_electron_count)
    return 1;

  if (pi_electron_count < 6)
    return 0;

  if (2 == (pi_electron_count % 4))
    return 1;

  int extra_pi_electrons = kta.additional_fused_pi_electrons();

  if (extra_pi_electrons)
  {
    pi_electron_count += extra_pi_electrons;

    if (pi_electron_count < 6)
      return 0;

    if (2 == (pi_electron_count % 4))
      return 1;
  }

  if (file_scope_allow_any_even_number_of_pi_electrons && 0 == pi_electron_count % 2)
    return 1;

  return 0;
}

/*
  Atom A is in a possibly aromatic ring, determine the
  number of pi electrons it will contribute
*/

int
Molecule::_count_pi_electrons(atom_number_t zatom, const int * process_these_atoms, int & result)
{
  Atom * a = _things[zatom];

  if (6 == a->atomic_number())    // the most common case
  {
    if (0 == a->formal_charge())
      result = 1;
    else if (-1 == a->formal_charge())
    {
      if (3 == a->ncon() || a->implicit_hydrogens())
        result = 2;
      else
        result = 1;
    }
    else
      result = 1;
  }
  else if (! a->pi_electrons(result))
    return 0;

  //cerr << "Atom " << zatom << " computed to have " << result << " pi electrons\n";

  const int acon = a->ncon();

  // make special allowance for N=O groups
  // Probably should consolidate all the 4 connected Sulphur rules, monitor behaviour...

  if (2 == acon)    // cannot have any bonds outside the ring system
    ;
  else if (_bonded_to_heteroatom_outside_system(zatom, process_these_atoms))
  {
    //  cerr << "Atom " << zatom << " bonded to heteroatom outside ring\n";
    if (16 == a->atomic_number() && 4 == acon && 5 == a->nbonds() &&
        allow_pipeline_pilot_aromaticity_on_input && 0 == a->formal_charge())
      result = 1;
    else if (16 == a->atomic_number() && 4 == acon && 6 == a->nbonds() &&
             allow_pipeline_pilot_aromaticity_on_input && 0 == a->formal_charge())
      result = 1;
    else
      result = 0;

    //  cerr << "Pi electrons outside ring for S " << result << endl;
  }
  else if (7 == a->atomic_number() && 3 == acon && a->nbonds() >= 4 && 0 == a->formal_charge())
    result = 1;

  return 1;
}

/*
  During Kekule assignment, we try to avoid putting double bonds in small carbocycles
*/

int
Molecule::_bond_only_in_small_carbocycle(atom_number_t a0, atom_number_t a1)
{
  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    if (! ri->contains_bond(a0, a1))
      continue;

    if (ri->number_elements() > 5)    // in a larger ring
      return 0;

    if (count_heteroatoms(*ri) > 0)    // found in a non-carbocycle
      return 0;
  }

  return 1;
}

/*
  Someone is trying to choose which way to go. We always direct them into the
  smallest ring. We assume just two items in P, but don't enforce that...

  We try to avoid placing double bonds in small carbocycles

*/

int
Molecule::_sort_possible_continuations_by_ring_size(atom_number_t anchor, Set_of_Atoms & p)
{
  assert(p.number_elements() > 1);

  atom_number_t a0 = p[0];
  atom_number_t a1 = p[1];

  if (_bond_only_in_small_carbocycle(anchor, a0))
  {
    p.swap_elements(0, 1);
    return 1;
  }

  if (_bond_only_in_small_carbocycle(anchor, a1))
    return 1;

  resizable_array<int> ring_sizes;

  for (int i = 0; i < p.number_elements(); i++)
  {
    atom_number_t j = p[i];

    ring_sizes.add(ring_containing_atom(j)->number_elements());
  }

  if (ring_sizes[0] > ring_sizes[1])
    p.swap_elements(0, 1);

  return 1;
}

/*
  During kekule determination, the process_these array serves several purposes.
  On input, all atoms to be determined should have their PROCESS_THESE
  entry set to KEKULE_READY_TO_PROCESS.

*/

int
Molecule::__find_kekule_form_current_config(const resizable_array<const Ring *> & rings,
                                            Kekule_Temporary_Arrays & kta, atom_number_t zatom,
                                            int pi_electron_count)
{
  int * process_these_atoms = kta.process_these_atoms();

  assert(KEKULE_PROCESSED == process_these_atoms[zatom]);

  // Need to be a little careful. This atom may not yet be in final form,
  // so asking for pi electrons for a carbon, may be misleading

  int pi;
  if (! _count_pi_electrons(zatom, process_these_atoms, pi))
    return 0;

  pi_electron_count += pi;

  Atom * a = _things[zatom];

#ifdef DEBUG_KEKULE
  cerr << "current_config continues with atom " << zatom << " ncon = " << a->ncon()
       << " pi = " << pi << " pi count = " << pi_electron_count << endl;
#endif

  // Dec 2004. Heuristic to first follow things in small rings

  Set_of_Atoms possible_continuations;

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (KEKULE_READY_TO_PROCESS == process_these_atoms[j])
      possible_continuations.add(j);
  }

  if (possible_continuations.number_elements() > 1)
    _sort_possible_continuations_by_ring_size(zatom, possible_continuations);

  // Why is this a loop, it only processes the first one?

  for (int i = 0; i < possible_continuations.number_elements(); i++)
  {
    atom_number_t j = possible_continuations[i];

#ifdef DEBUG_KEKULE
    cerr << "From atom " << zatom << " continuing kekule search with atom " << j << endl;
#endif

    return _find_kekule_form(rings, kta, j, pi_electron_count);
  }

  // Could not find a continuation atom directly bonded to atom A,
  // look for one elsewhere in the molecule

  atom_number_t astart;
#ifdef DEBUG_KEKULE
  cerr << "Looking for the enxt start atom, astart = " << astart << '\n';
  if (_kekule_identify_next_atom(process_these_atoms, astart))
    cerr << "Atom " << astart << " identified as new start atom\n";
#endif

  if (_kekule_identify_next_atom(process_these_atoms, astart))
    return _find_kekule_form(rings, kta, astart, pi_electron_count);

  // Cound not find any unprocessed atoms. Is this a kekule form????

#ifdef DEBUG_KEKULE
  cerr << "Checking for possible aromatic " << pi_electron_count << endl;
#endif

  if (looks_like_valid_kekule_system(kta, pi_electron_count))
    return 1;

  if (aromatic::any_kekule_bonding_pattern_ok_for_aromatic_input) {
    return 1;
  }

  if (pi_electron_count < 6)
    return 0;

  // When I first implemented this, I rejected things which had an odd number
  // of electrons, but I quickly ran into a counter-example,
  // N12C(=NC3=CC=CC=C3C1=O)SCC2=O PBCHM11481368

  // Now things get interesting. We do not have a 4n + 2 system, but the
  // system may consist of aromatic rings.

  int rc = _kekule_arom_test_rings(kta, rings);

#ifdef DEBUG_KEKULE
  cerr << "Final determination " << rc << endl;
#endif

  return rc;
}

/*
  As we examine atoms in the molecule, we need to mark each atom as having been
  processed. As we back out, we must release that hold.
  This function centralises those actions;
*/

int
Molecule::_find_kekule_form_current_config(const resizable_array<const Ring *> & rings,
                                           Kekule_Temporary_Arrays & kta, atom_number_t a,
                                           int pi_electron_count)
{
  //cerr << "In _find_kekule_form_current_config\n";
  int * process_these_atoms = kta.process_these_atoms();

  assert(KEKULE_READY_TO_PROCESS == process_these_atoms[a]);

  process_these_atoms[a] = KEKULE_PROCESSED;

  int rc = __find_kekule_form_current_config(rings, kta, a, pi_electron_count);

  process_these_atoms[a] = KEKULE_READY_TO_PROCESS;

  //cerr << "_find_kekule_form_current atom " << a << " returning " << rc << endl;
  return rc;
}

/*
  Twiddle bonds and implicit hydrogens in order to find a kekule form
*/

int
Molecule::_find_kekule_form(const resizable_array<const Ring *> & rings,
                            Kekule_Temporary_Arrays & kta, atom_number_t zatom,
                            int pi_electron_count)
{
  int * process_these_atoms = kta.process_these_atoms();
  const int * implicit_hydrogens_needed = kta.implicit_hydrogens_needed();

  assert(KEKULE_READY_TO_PROCESS == process_these_atoms[zatom]);

  int * vary_bonds = kta.vary_bonds();

#ifdef DEBUG_KEKULE
  cerr << "Into _find_kekule_form, atom " << zatom << " vary_bonds " << vary_bonds[zatom] << endl;
#endif

  if (0 == vary_bonds[zatom])
    return _find_kekule_form_current_config(rings, kta, zatom, pi_electron_count);

  int * vary_hcount = kta.vary_hcount();

  Atom * a = _things[zatom];

#ifdef DEBUG_KEKULE
  cerr << "_find_kekule_form:continuing with atom " << zatom << " (" << a->atomic_symbol()
       << ") nbonds = " << nbonds(zatom) << " hcount = " << a->implicit_hydrogens()
       << " vary hcount = " << vary_hcount[zatom] << " vary bonds = " << vary_bonds[zatom] << endl;
#endif

  // first try no implicit hydrogens
  if (vary_hcount[zatom]) {
    a->set_implicit_hydrogens(0);
  }

  int bonds_found = 0;
  int unprocessed_neighbours = 0;

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    Bond * b = a->item(i);

    if (b->is_single_bond())
      bonds_found++;
    else if (b->is_double_bond())
      bonds_found += 2;

    atom_number_t j = b->other(zatom);

    if (KEKULE_READY_TO_PROCESS == process_these_atoms[j])    // other connection not yet processed
    {
      unprocessed_neighbours++;
      //    if (! b->is_single_bond())     // why was I doing this?
      //      b->set_bond_type (SINGLE_BOND);
    }

#ifdef DEBUG_KEKULE
    cerr << "  Atom " << j;
    if (b->is_single_bond())
      cerr << " single bond";
    else if (b->is_double_bond())
      cerr << " double bond";
    else
      cerr << " unrecognised bond";
    if (KEKULE_READY_TO_PROCESS == process_these_atoms[j] && vary_bonds[j])
      cerr << " varying";
    cerr << ", process? " << (KEKULE_READY_TO_PROCESS == process_these_atoms[j]) << endl;
#endif
  }

  if (unprocessed_neighbours)
    a->recompute_nbonds();

#ifdef DEBUG_KEKULE
  if (a->implicit_hydrogens())
    cerr << "  and " << a->implicit_hydrogens() << " implicit hydrogens\n";
#endif

  const Element * e = a->element();

  if (3 == acon)    // cannot be any implicit hydrogens on a 3 connected ring atom
    ;
  else if (bonds_found >= 4)
    ;
  else if (implicit_hydrogens_needed[zatom] >= 0)    // known from input
    bonds_found += implicit_hydrogens_needed[zatom];
  else if (a->implicit_hydrogens())    // computed/stored value likely to be wrong
    bonds_found++;

  int bonds_needed = e->normal_valence() + a->formal_charge();
  if (bonds_needed < a->nbonds() && e->number_alternate_valences())
  {
    bonds_needed = e->alternate_valence(0) + a->formal_charge();
    if (bonds_needed < a->nbonds() && e->number_alternate_valences() > 1)
      bonds_needed = e->alternate_valence(1) + a->formal_charge();
  }

#ifdef DEBUG_KEKULE
  cerr << "  Has " << bonds_found << " bonds, needs " << bonds_needed << " from "
       << unprocessed_neighbours << " connections\n";
  cerr << "Needs " << implicit_hydrogens_needed[zatom] << " implicit hydrogens\n";
#endif

  if (bonds_found <= bonds_needed)    // great
    ;
  else if (6 == bonds_found && 4 == bonds_needed && 16 == e->atomic_number() &&
           allow_pipeline_pilot_aromaticity_on_input)
    ;
  else if (bonds_found > bonds_needed)    // obviously wrong
    return 0;

  // Intercept c+ here

  if (5 == bonds_needed && 6 == e->atomic_number() && 1 == a->formal_charge() &&
      (2 == acon || 3 == acon))
    bonds_needed = 3;

  //cerr << "found " << bonds_found << " need " << bonds_needed << endl;

  // We only allow the addition of one extra bond over what might be there already

  //assert (bonds_found == bonds_needed || (bonds_found + 1 == bonds_needed));

  if (bonds_found == bonds_needed || (bonds_found + 1 == bonds_needed))
    ;
  else if (6 == bonds_found && 4 == bonds_needed && 16 == e->atomic_number() &&
           allow_pipeline_pilot_aromaticity_on_input)
    ;
  else
  {
    if (file_scope_display_no_kekule_form_message)
    {
      cerr << _molecule_name << " atom " << zatom << " (" << atomic_symbol(zatom) << " D"
           << a->ncon() << ") has " << bonds_found << " and needs " << bonds_needed << endl;
      cerr << "Atom has formal charge " << formal_charge(zatom) << " and "
           << implicit_hydrogens(zatom) << " implicit hydrogens\n";
      cerr << "Impossible\n";
    }
    return 0;
  }

  // If all our neighbours have been processed, can we add a Hydrogen

  if (0 == unprocessed_neighbours && bonds_found < bonds_needed)
  {
    if (0 == vary_hcount[zatom])    // not varying the H count here, oh well
      return 0;

    if (a->implicit_hydrogens())    // already has an H
      return 0;

    a->set_implicit_hydrogens(1);
    bonds_found++;
  }

  //cerr << "At line " << __LINE__ << " atom has " << a->implicit_hydrogens() << " implicit hydrogens\n";

  // If we already have a full complement of bonds, proceed with this configuration

  if (bonds_found >= bonds_needed)
    return _find_kekule_form_current_config(rings, kta, zatom, pi_electron_count);

  // At this stage, we are short a bond. That bond must be gained either by
  // placing a bond, or adding a hydrogen.

  // always try placing a bond first.

  // Dec 2004. Heuristic to always place the bond in the smallest ring

  Set_of_Atoms try_double_bond_to;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

#ifdef DEBUG_KEKULE
    cerr << "From atom " << zatom << " what about atom " << j << " vary_bonds = " << vary_bonds[j]
         << endl;
#endif

    if (KEKULE_READY_TO_PROCESS != process_these_atoms[j] || 0 == vary_bonds[j])
      continue;

    Atom * aj = _things[j];

    if (7 == aj->atomic_number() && 3 == aj->nbonds() &&
        0 == aj->formal_charge())    // never create [nD4v4]
      continue;

    //  Make sure we handle O=c1[nH]c(=O)[n][n+]c1 PBCHM56995681
    //  We may need to form =[N+]= within an aromatic ring. Strange...

    if (_is_nitrogen_double_bond_to_something_outside_ring(j))
      ;
    else if (7 == aj->atomic_number() && 2 == aj->ncon() && 1 == aj->formal_charge() &&
             3 == aj->nbonds() && 0 == aj->implicit_hydrogens())
      ;
    else if (aj->unsaturated())
      continue;

    try_double_bond_to.add(j);
  }

#ifdef DEBUG_KEKULE
  cerr << "Found " << try_double_bond_to.number_elements() << " possible places for double bonds\n";
#endif

  // If more than one possibility, shuffle to get smallest ring size first. Note that I only allow
  // for two possibilities

  if (try_double_bond_to.number_elements() > 1)
    _sort_possible_continuations_by_ring_size(zatom, try_double_bond_to);

  for (int i = 0; i < try_double_bond_to.number_elements(); i++)
  {
    atom_number_t j = try_double_bond_to[i];

    Bond * b = const_cast<Bond *>(a->bond_to_atom(j));

#ifdef DEBUG_KEKULE
    cerr << "From atom " << zatom << " what about atom " << j << " vary_bonds = " << vary_bonds[j]
         << endl;
#endif

    Atom * aj = _things[j];

    if (_is_nitrogen_double_bond_to_something_outside_ring(j))
      ;
    else if (7 == aj->atomic_number() && 2 == aj->ncon() && 1 == aj->formal_charge() &&
             3 == aj->nbonds() && 0 == aj->implicit_hydrogens())
      ;
    else if (16 == aj->atomic_number() && 4 == aj->ncon() &&
             allow_pipeline_pilot_aromaticity_on_input && 5 == aj->nbonds() &&
             0 == aj->formal_charge())
      ;
    else if (7 == aj->atomic_number() && 3 == aj->ncon() && 0 == aj->formal_charge() &&
             3 == aj->nbonds())
      continue;
    else if (aj->nbonds() > aj->ncon())    // likely already been set recursively
      continue;

    b->set_bond_type(DOUBLE_BOND);
    a->invalidate_computed_values_after_bond_change();
    aj->invalidate_computed_values_after_bond_change();

#ifdef DEBUG_KEKULE
    cerr << "Just set double bond between atoms " << zatom << " and " << j << endl;
#endif

    if (_find_kekule_form_current_config(rings, kta, zatom, pi_electron_count))
      return 1;

#ifdef DEBUG_KEKULE
    cerr << "Reset bond between " << zatom << " and " << j << " back to single\n";
#endif

    b->set_bond_type(SINGLE_BOND);
    a->invalidate_computed_values_after_bond_change();
    aj->invalidate_computed_values_after_bond_change();
  }

#ifdef DEBUG_KEKULE
  cerr << "Could not place double bond on atom " << zatom << endl;
#endif

  // If we don't have enough bonds, but we are allowed to vary the hcount, try adding
  // a hydrogen

  if (vary_hcount[zatom])
  {
#ifdef DEBUG_KEKULE
    cerr << "Adding a hydrogen to atom " << zatom << endl;
#endif

    a->set_implicit_hydrogens(1);
    if (_find_kekule_form_current_config(rings, kta, zatom, pi_electron_count))
      return 1;

#ifdef DEBUG_KEKULE
    cerr << "Removing implicit hydrogen from atom " << zatom << endl;
#endif

    a->set_implicit_hydrogens(0);
  }

#ifdef DEBUG_KEKULE
  cerr << "No luck with atom " << zatom << endl;
#endif

  return 0;
}

/*
  We have marked a set of atoms to be processed.

  Set the process_these array to KEKULE_READY_TO_PROCESS, and identify our
  starting atom.  Note that we choose an atom with variability (of some kind)
  if at all possible.  The reason for this is that I ran into a case where
  the first atom had no variability, so it was not possible to construct a kekule
  system from it - its bonding or hcount did not change!

*/

int
Molecule::_identify_kekule_search_starting_atom(Kekule_Temporary_Arrays & kta,
                                                atom_number_t & astart) const
{
  int * process_these_atoms = kta.process_these_atoms();
  const int * vary_hcount = kta.vary_hcount();
  const int * vary_bonds = kta.vary_bonds();

  astart = INVALID_ATOM_NUMBER;

  int system_size = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (process_these_atoms[i])
    {
      system_size++;
      process_these_atoms[i] = KEKULE_READY_TO_PROCESS;
      if (vary_hcount[i] || vary_bonds[i])
      {
        if (0 == _things[i]->formal_charge())    // always prefer a neutral starting atom
          astart = i;
        else if (INVALID_ATOM_NUMBER == astart)
          astart = i;
      }
      else if (INVALID_ATOM_NUMBER == astart)
        astart = i;
    }
  }

  assert(INVALID_ATOM_NUMBER != astart);

  return system_size;
}

/*
*/

int
Molecule::_find_kekule_form_ring_system(Kekule_Temporary_Arrays & kta,
                                        resizable_array<const Ring *> & rings)
{
  int * process_these_atoms = kta.process_these_atoms();

#ifdef DEBUG_KEKULE
  iw_write_array(process_these_atoms, _number_elements, "process_these_atoms", cerr);
#endif

  // If 'a' atoms are present, test to see if they are present in this set of rings

  //cerr << "'a' atoms " << kta.a_atoms_found() << endl;

  if (kta.a_atoms_found() && _a_atoms_found_this_set_of_rings(kta, rings))
    return _kekule_all_atoms_and_bonds_in_system_aromatic(kta, rings);

  if (! _do_obvious_hcount_adjustments(kta))
  {
    if (file_scope_display_no_kekule_form_message)
      cerr << "Molecule::_find_kekule_form: bad hcount detected\n";
    return 0;
  }

  if (! _do_obvious_bond_order_settings(kta))
  {
    if (file_scope_display_no_kekule_form_message)
      cerr << "Molecule::_find_kekule_form: bad bonding detected\n";
    return 0;
  }

  // The case of PBCHM597765 C1=CC=C2C(=C1)C=CN3C2=NC=C3, or c1cnc2n1C=CC3=C2C=CC=C3

  _count_pi_electrons_in_any_fused_but_kekule_form_rings(kta, rings);

#ifdef DEBUG_KEKULE
  cerr << "Start kekule search\n";
  for (int i = 0; i < _number_elements; i++)
  {
    if (! process_these_atoms[i])
      continue;

    cerr << "Atom " << i << " (" << atomic_symbol(i) << ") ";
    if (implicit_hydrogens(i))
      cerr << implicit_hydrogens(i) << " implicit_hydrogens ";
    if (kta.vary_hcount()[i])
      cerr << "vary hcount ";
    if (kta.vary_bonds()[i])
      cerr << "vary bonds ";
    cerr << endl;
  }
#endif

  atom_number_t astart;
  int system_size = _identify_kekule_search_starting_atom(kta, astart);

#ifdef DEBUG_KEKULE
  cerr << "Starting Kekule search with atom " << astart << " '"
       << smarts_equivalent_for_atom(astart) << "', system_size " << system_size << " atoms\n";
#endif

  int pi_electron_count = 0;

  int rc = _find_kekule_form(rings, kta, astart, pi_electron_count);

  // Feb 2022. Open question, should this be done if the calculation has failed?
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    Atom * a = _things[i];

    a->recompute_nbonds();
    if (! a->implicit_hydrogens_known())
      continue;

    //  If the computed implicit hydrogens are the same as the number specified, no need for the known flag

    int ih;
    a->compute_implicit_hydrogens(ih);
    if (ih == a->implicit_hydrogens())
      a->set_implicit_hydrogens_known(0);
  }

  if (rc)
    return rc;

  // Could not do the aromaticity of the whole system, can we do it in
  // subsets.  Actually, that's a little hard now, let's just do it one
  // ring at a time.  We record success only if all rings are
  // individually found to be aromatic.  Obviously this is wrong.  Fix
  // it if it ever becomes an issue.  Note too that we may leave a
  // partial result here...

  int nr = rings.number_elements();

  if (1 == nr)    // Just one ring to start with, cannot do subsets
  {
    if (file_scope_display_no_kekule_form_message)
      cerr << "Molecule::_find_kekule_form: no kekule form for " << system_size << " atoms\n";
    return _maybe_all_atoms_are_permanent_aromatic(kta, rings);
  }

  /*
  If the algorithm starts with the Sulphur containing ring, we fail, because we need
  the double bond to the Nitrogen outside the ring. Therefore, we do an iterative
  process

  But since we are not doing an interative process now, we turn off the iterative stuff 
  for efficiency
  */

#ifdef CODE_FOR_WHEN_WE_DO_PAIRWISE_FINDING
  int * ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  while (1)
  {
    int rings_found_aromatic_this_iteration = 0;

    for (int i = 0; i < nr; i++)
    {
      if (ring_already_done[i])
        continue;

      const Ring * ri = rings[i];

#ifdef DEBUG_KEKULE
      cerr << "Doing " << (*ri) << " separately\n";
#endif

      set_vector(process_these_atoms, _number_elements, 0);

      ri->set_vector(process_these_atoms, 1);

      resizable_array<const Ring *> tmp;
      tmp.add(ri);

      if (_find_kekule_form_ring_system(kta, tmp))    // no recursion because only one item in tmp
      {
        rings_found_aromatic_this_iteration++;
        ring_already_done[i] = 1;
        rc++;
      }
    }

    if (0 == rings_found_aromatic_this_iteration)
      break;
  }
#else
  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = rings[i];
    set_vector(process_these_atoms, _number_elements, 0);

    ri->set_vector(process_these_atoms, 1);

    resizable_array<const Ring *> tmp;
    tmp.add(ri);

    if (_find_kekule_form_ring_system(kta, tmp))    // no recursion because only one item in tmp
    {
      rc++;
    }
  }
#endif

#ifdef DEBUG_KEKULE
  cerr << "rc = " << rc << endl;
#endif

  if (rc == nr)
    return rc;

  if (file_scope_display_no_kekule_form_message)
    cerr << "Molecule::_find_kekule_form: no kekule form for " << system_size << " atoms\n";

  return _maybe_all_atoms_are_permanent_aromatic(kta, rings);
}

static void
unset_all_atoms_in_rings(const resizable_array<const Ring *> & rings, int * aromatic_atoms,
                         int * process_these_rings)
{
  for (int i = 0; i < rings.number_elements(); i++)
  {
    rings[i]->set_vector(aromatic_atoms, 0);

    int rn = rings[i]->ring_number();

    process_these_rings[rn] = 0;
  }

  return;
}

//#define DEBUG_FIND_NEXT_RING_OR_RING_SYSTEM

int
Molecule::_find_next_ring_or_ring_system(Kekule_Temporary_Arrays & kta, int flag,
                                         resizable_array<const Ring *> & rings)
{
  int * process_these_atoms = kta.process_these_atoms();

  set_vector(process_these_atoms, _number_elements, 0);

  const int * process_these_rings = kta.process_these_rings();

  int nr = _sssr_rings.number_elements();

#ifdef DEBUG_FIND_NEXT_RING_OR_RING_SYSTEM
  cerr << "In _find_next_ring_or_ring_system, flag " << flag << endl;
  iw_write_array(process_these_rings, nr, "process_these_rings", cerr);
#endif

  if (rings.number_elements())
    rings.resize_keep_storage(0);
  else
    rings.resize(nr);

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    if (0 == process_these_rings[i])
      continue;

    if (flag != process_these_rings[i])
      continue;

    const Ring * ri = ringi(i);

#ifdef DEBUG_FIND_NEXT_RING_OR_RING_SYSTEM
    cerr << "Doing ring " << i << " " << (*ri) << endl;
#endif

    ri->set_vector(process_these_atoms, 1);
    rings.add(ri);
    rc++;
  }

  return rc;
}

//#define DEBUG_FIND_KEKULE_FORM_MAIN

int
Molecule::_find_kekule_form(Kekule_Temporary_Arrays & kta)
{
  int failures = 0;
  resizable_array<const Ring *> rings;

  int flag = 2;
  while (_find_next_ring_or_ring_system(kta, flag, rings))
  {
    flag++;

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
    cerr << "Ring system contains " << rings.number_elements() << " rings\n";
#endif

    if (_find_kekule_form_ring_system(kta, rings))
      unset_all_atoms_in_rings(rings, kta.aromatic_atoms(),
                               kta.process_these_rings());    // mark atoms and rings as done
    else
      failures++;

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
    cerr << "After iteration, failures = " << failures << endl;
    for (int i = 0; i < _number_elements; i++)
    {
      cerr << "i = " << i << " aromatic = " << kta.aromatic_atoms()[i] << endl;
    }
#endif
  }

  return 0 == failures;
}

static const Element * aromatic_carbon = nullptr;
static const Element * aromatic_nitrogen = nullptr;
static const Element * aromatic_oxygen = nullptr;
static const Element * aromatic_phosphorus = nullptr;
static const Element * aromatic_sulphur = nullptr;
static const Element * aromatic_a = nullptr;

static void
my_fetch_or_create_permanent_aromatic_element(const char * s, const Element *& e,
                                              const Element * copy_element_data_from)
{
  if (nullptr != e)
    return;

  e = get_element_from_symbol_no_case_conversion(s);
  if (nullptr == e)
    e = create_element_with_symbol(s);

  assert(nullptr != e);

  if (nullptr != copy_element_data_from)
    const_cast<Element *>(e)->copy_element_data(copy_element_data_from);

  const_cast<Element *>(e)->set_permanent_aromatic(1);

  return;
}

int
Molecule::_convert_any_chain_aromatic_atoms_to_permanent_aromatic(const int * process_these_atoms)
{

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (nrings(i))
      continue;

    _convert_atom_to_permanent_aromatic(i);
  }

  return 1;
}

static void
create_permanent_aromatic_elements()
{
  my_fetch_or_create_permanent_aromatic_element("c", aromatic_carbon,
                                                get_element_from_atomic_number(6));
  my_fetch_or_create_permanent_aromatic_element("n", aromatic_nitrogen,
                                                get_element_from_atomic_number(7));
  my_fetch_or_create_permanent_aromatic_element("o", aromatic_oxygen,
                                                get_element_from_atomic_number(8));
  my_fetch_or_create_permanent_aromatic_element("p", aromatic_phosphorus,
                                                get_element_from_atomic_number(15));
  my_fetch_or_create_permanent_aromatic_element("s", aromatic_sulphur,
                                                get_element_from_atomic_number(16));
  my_fetch_or_create_permanent_aromatic_element("a", aromatic_a, nullptr);

  return;
}

int
Molecule::_convert_atom_to_permanent_aromatic(atom_number_t zatom)
{
  if (nullptr == aromatic_carbon)
    create_permanent_aromatic_elements();

  Atom * a = _things[zatom];

  atomic_number_t z = a->atomic_number();

  //  cerr << "Converting atom " << zatom << " to permanent aromatic\n";

  if (6 == z)
    a->set_element(aromatic_carbon);
  else if (7 == z)
    a->set_element(aromatic_nitrogen);
  else if (8 == z)
    a->set_element(aromatic_oxygen);
  else if (16 == z)
    a->set_element(aromatic_sulphur);
  else if (15 == z)
    a->set_element(aromatic_phosphorus);
  else if ("a" == a->element()->symbol())
    a->set_element(aromatic_a);
  else if (a->element()->permanent_aromatic())
    ;
  else
    cerr << "Molecule::_convert_atom_to_permanent_aromatic:what to do with atomic number " << z
         << endl;

  return 1;
}

/*
  Do some rudimentary tests on whether or not any of the atoms in
  `process_these_atoms` can be aromatic.
  As atoms are examined, update `aromatic_atoms_found`.
  Return 0 if any of the atoms in `process_these_atoms` cannot be aromatic.
*/

int
Molecule::_kekule_could_be_aromatic(int * process_these_atoms, int & aromatic_atoms_found)
{

  // To avoid recomputing ring membership, keep a copy of the initial ring membmership
  int * rm = new int[_number_elements]; std::unique_ptr<int[]> free_rm(rm);

  ring_membership(rm);    // force SSSR

  if (0 == nrings())
  {
    if (convert_chain_aromatic_bonds())
    {
      (void)_convert_chain_aromatic_bonds(process_these_atoms, rm);
      return 1;    // we are happy
    }
    else if (aromatic_chain_bonds_are_ok())
    {
      _convert_any_chain_aromatic_atoms_to_permanent_aromatic(process_these_atoms);
      set_vector(process_these_atoms, _number_elements, 0);
      return 1;    // that's good!
    }
    else
    {
      if (warn_aromatic_chain_atoms)
        cerr << "Molecule contains no rings, but has aromatic atoms. IMPOSSIBLE\n";
      return 0;    // this cannot be aromatic
    }
  }

  // Function _convert_chain_aromatic_bonds does all the chain atoms in
  // the molecule, so we need to keep track of whether or not we have
  // called it

  int called_convert_chain_aromatic_bonds = 0;

  aromatic_atoms_found = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (0 == rm[i])
    {
      if (convert_chain_aromatic_bonds())
        ;
      else if (aromatic_chain_bonds_are_ok())
      {
        process_these_atoms[i] = 0;
        _convert_atom_to_permanent_aromatic(i);
      }
      else
      {
        if (warn_aromatic_chain_atoms)
          cerr << "Molecule::_kekule_could_be_aromatic: atom " << i << " is not in a ring\n";
        return 0;
      }

      if (! called_convert_chain_aromatic_bonds) {
        (void)_convert_chain_aromatic_bonds(process_these_atoms, rm);
        called_convert_chain_aromatic_bonds = 1;
      }

      ++aromatic_atoms_found;
      continue;  // No further checking below...
    }

    Atom * a = _things[i];

    if (a->implicit_hydrogens_known() && a->implicit_hydrogens() > 1) {
      if (display_unusual_hcount_warning_messages()) {
        cerr << "Molecule::_kekule_cannot_be_aromatic: atom " << i << " has "
             << a->implicit_hydrogens() << " KNOWN implicit hydrogens\n";
      }
      return 0;
    }

    const int icon = a->ncon();

    // Determination of attached Hydrogen is complicated by the fact that the Kekule
    // bonds are not present at this point. Mostly ignore...

    if (icon > 3) {
      if (allow_pipeline_pilot_aromaticity_on_input && 4 == icon && 16 == a->atomic_number() &&
          1 == rm[i] && (a->nbonds() > 4) && 0 == a->formal_charge())
        ;
      else if (allow_pipeline_pilot_aromaticity_on_input && 4 == icon && 15 == a->atomic_number() &&
               4 == a->nbonds() && 0 == a->formal_charge())
        ;
      else {
        if (warn_aromatic_chain_atoms) {  // not exactly, but OK
          cerr << "Molecule::_kekule_could_be_aromatic: atom " << i << " has " << icon << " connections\n";
        }
        return 0;
      }
    }

    aromatic_atoms_found++;
  }

  //cerr << "Found " << aromatic_atoms_found << " aromatic_atoms_found\n";

  if (0 == aromatic_atoms_found)
    return 1;    // no aromatic atoms specified, no further processing

  // Apr 2001. Remove this test because of molecules like Cn1on1C which are aromatic
  // in the Daylight world

  //if (aromatic_atoms_found < 4)  // no 3 membered aromatic rings
  //{
  //  cerr << "Molecule::_kekule_could_be_aromatic: only " << aromatic_atoms_found << " aromatic atoms\n";
  //  return 0;
  //}

  return 1;    // this collection of atoms could all be aromatic
}

/*
  This ring will NOT be aromatic. Any bonds that are unshared with other
  rings cannot be aromatic
*/

int
Molecule::_make_all_unshared_bonds_non_aromatic(const Ring & r)
{
  for (Ring_Bond_Iterator i(r); i != r.zend(); i++)
  {
    atom_number_t a1 = i.a1();
    atom_number_t a2 = i.a2();

    Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));

    if (! b->is_permanent_aromatic())
      continue;

    if (b->nrings() > 1)
      continue;

    b->set_permanent_aromatic(0);
  }

  return 1;
}

/*
  Very similar functions, consolidate?
*/

int
Molecule::_make_bonds_aromatic(const Ring * zring)
{
  int rsize = zring->number_elements();

  const atom_number_t * r = zring->rawdata();

  for (int i = 0; i < rsize; i++)
  {
    atom_number_t a1 = r[i];

    atom_number_t a2;

    if (0 == i)
      a2 = r[rsize - 1];
    else
      a2 = r[i - 1];

    Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));

    b->set_aromatic();
  }

  return 1;
}

int
Molecule::_make_all_bonds_aromatic(const Ring & r)
{
  for (Ring_Bond_Iterator i(r); i != r.zend(); i++)
  {
    atom_number_t a1 = i.a1();
    atom_number_t a2 = i.a2();

    Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));

    if (! b->is_permanent_aromatic())
      b->set_permanent_aromatic(1);
  }

  return 1;
}

/*
  We have found a ring that has some, but not all, bonds aromatic. We
  are trying to figure out whether or not that ring should be converted
  into a wholly aromatic ring. If we find an aromatic bond that is not
  shared with another ring, then we will want to try to aromatise this
  whole ring
*/

int
Molecule::_count_aromatic_bonds_in_just_one_ring(const Ring & r)
{
  int rc = 0;
  for (Ring_Bond_Iterator i(r); i != r.zend(); i++)
  {
    atom_number_t a1 = i.a1();
    atom_number_t a2 = i.a2();

    const Bond * b = _things[a1]->bond_to_atom(a2);

    if (! b->is_permanent_aromatic())
      continue;

    if (1 == b->nrings())    // aromatic bond not shared
      rc++;
  }

  return rc;
}

int
Molecule::_all_bonds_aromatic(const Ring & r) const
{
  for (Ring_Bond_Iterator i(r); i != r.zend(); i++)
  {
    atom_number_t a1 = i.a1();
    atom_number_t a2 = i.a2();

    const Bond * b = _things[a1]->bond_to_atom(a2);

    //  cerr << "Bond between " << a1 << " and " << a2 << " aromatic? " << b->is_permanent_aromatic() << endl;

    if (! b->is_permanent_aromatic())
      return 0;
  }

  return 1;
}

/*
  Originally this was designed to just discard any ring that had aliphatic bonds in
  it, but when Skelgen came along, I found examples where the atoms were not aromatic,
  and some of the bonds in a ring were not aromatic
*/

int
Molecule::_kekule_check_rings_containing_aliphatic_bonds(Kekule_Temporary_Arrays & kta)
{
  int * process_these_rings = kta.process_these_rings();

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    if (0 == process_these_rings[i])
      continue;

    const Ring * ri = ringi(i);

    if (_all_bonds_aromatic(*ri))    // excellent, as it should be
      continue;

//  cerr << _molecule_name << " not all bonds aromatic " << (*ri) << endl;

    if (all_bonds_in_aromatic_ring_must_be_aromatic)    // cannot make aromatic
    {
      _make_all_unshared_bonds_non_aromatic(*ri);

      process_these_rings[i] = 0;
      continue;
    }

    int abj1r = _count_aromatic_bonds_in_just_one_ring(*ri);    // not shared with another ring

    if (0 == abj1r)
    {
      process_these_rings[i] = 0;
      continue;
    }

    //  cerr << "Ring " << (*ri) << " has " << abj1r << " unshared aromatic bonds\n";

    _make_all_bonds_aromatic(*ri);
  }

  return 1;
}

/*
  We have just read the molecule and try to decide as much about ultimate Hydrogen counts
  as possible.
  Note that explicit hydrogens have not yet been detached
*/

int
Molecule::_determine_hydrogen_counts(Kekule_Temporary_Arrays & kta)
{
  const int * aromatic_atoms = kta.aromatic_atoms();

  int * implicit_hydrogens_needed = kta.implicit_hydrogens_needed();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == aromatic_atoms[i])
      continue;

    Atom * a = _things[i];

    int acon = a->ncon();

    const atomic_number_t z = a->atomic_number();

    const int exph = explicit_hydrogens(i);

    //  cerr << "Checking aromatic atom " << i << " type " << a->atomic_symbol() << " exph " << exph << " acon " << acon << " Hknown " <<   a->implicit_hydrogens_known() << " ih " << a->implicit_hydrogens() << endl;

    if (exph)    // dangerous if someone has just some of the Hydrogen atoms explicit
    {
      implicit_hydrogens_needed[i] = exph;
      continue;
    }

    if (3 == acon)    // after we check for explicit hydrogens
    {
      implicit_hydrogens_needed[i] = 0;
      continue;
    }

    if (a->implicit_hydrogens_known())
    {
      if (7 == a->atomic_number())    // the most common case
        implicit_hydrogens_needed[i] = a->implicit_hydrogens();
      else if (0 == a->implicit_hydrogens() && 16 == a->atomic_number() && 1 == a->formal_charge())
      {
        a->set_implicit_hydrogens_known(0);
        a->set_implicit_hydrogens(0, 1);
        implicit_hydrogens_needed[i] = 0;
      }
      else
        implicit_hydrogens_needed[i] = a->implicit_hydrogens();
      continue;
    }

    if (6 == z)
    {
      if (0 == a->formal_charge())
        implicit_hydrogens_needed[i] = 1;
      else if (-1 == a->formal_charge() && 2 == acon)
        implicit_hydrogens_needed[i] = 1;
      else
        implicit_hydrogens_needed[i] = 0;
    }
    else if (7 == z)
    {
      if (_h_unspecified_means_zero) {
        implicit_hydrogens_needed[i] = 0;
      } else if (nullptr == _atom_type)
        implicit_hydrogens_needed[i] = -1;              // cannot say, pyrrole or pyridine
      else if (ATOM_TYPE_NAR == _atom_type->item(i))    // Tripos aromatic Nitrogen, no hydrogens
        implicit_hydrogens_needed[i] = 0;
      else if (ATOM_TYPE_NPL3 == _atom_type->item(i))    // Tripos designation for Pyrrole
        implicit_hydrogens_needed[i] = 1;
      else
        implicit_hydrogens_needed[i] = -1;    // cannot figure it out
    }
    else if (16 == z)
    {
      if (-1 == a->formal_charge())
        implicit_hydrogens_needed[i] = 1;
      else if (1 == a->formal_charge() && 2 == a->ncon())
        implicit_hydrogens_needed[i] = 0;
      else if (4 == a->ncon())
        implicit_hydrogens_needed[i] = 0;
    }
    else if (15 == z)
      implicit_hydrogens_needed[i] = -1;
    else if (NOT_AN_ELEMENT != z)    // some other member of the periodic table
      implicit_hydrogens_needed[i] = 0;
    else if ("a" == a->element()->symbol())    // always aromatic
      implicit_hydrogens_needed[i] = -1;
  }

  return 1;
}

/*
  Aug 97. Ran into problems reading

                      O
                      ||
                      ||
                     /  \
                   \/    \/
               arom |    | arom
                    |    |
                   / \  / \
                      \/ arom
                      ||
                      ||
                       O

  where ALL atoms in this ring were entered as aromatic. Our
  current rules for aromaticity specify that this middle ring
  will be non-aromatic.
  We look for such rings, and when we find one, unset those
  atoms from the aromatic_atoms array
*/

int
Molecule::_kekule_suppress_non_aromatic_rings(Kekule_Temporary_Arrays & kta)
{
  int * aromatic_atoms = kta.aromatic_atoms();

  int rc = 0;

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = ringi(i);
    if (! r->is_fused())
      continue;

    if (! r->all_members_set_in_array(aromatic_atoms, 1))
      continue;

    int atoms_with_unshared_pi_electrons = 0;

    const int ring_size = r->number_elements();
    for (int j = 0; j < ring_size; j++)
    {
      const atom_number_t k = r->item(j);

      if (nrings(k) > 1)    // whatever pi electrons this atom has are shared: This is wrong, C1Cn2cn[nH]c2=N1 TODO monitor...
        continue;

      Atom * ak = _things[k];

      const int kcon = ak->ncon();
      if (kcon == 2 || 6 != ak->atomic_number() || kcon == ak->nbonds())
      {
        atoms_with_unshared_pi_electrons++;
        continue;
      }

      //    At this stage, we know atom K has a double bond outside the ring. Is it to a heteroatom?

      int heteroatom_outside_ring = 0;
      for (int l = 0; l < kcon; l++) {
        const Bond * b = ak->item(l);
        if (! b->is_double_bond())
          continue;

        atom_number_t o = b->other(k);
        if (r->contains(o)) {
          continue;
        }
        if (6 != _things[o]->atomic_number()) {
          heteroatom_outside_ring = 1;
          break;
        }
      }

      if (0 == heteroatom_outside_ring)
        atoms_with_unshared_pi_electrons++;
    }

    if (0 == atoms_with_unshared_pi_electrons)
    {
      for (int j = 0; j < ring_size; j++)
      {
        atom_number_t k = r->item(j);
        if (aromatic_atoms[k] && 1 == nrings(k))
        {
          aromatic_atoms[k] = 0;
          //        cerr << "Exclude atom " << k << " as non aromatic\n";
          rc++;
        }
      }
    }
  }

  return rc;
}

static int
all_aromatic_atoms_in_first_contained_in_second(const Ring & r1, const Ring & r2,
                                                const int * aromatic_atoms)
{
  int n = r1.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r1[i];

    if (0 == aromatic_atoms[j])    // we detect all aromatic in R1 contained in R2
      continue;

    if (! r2.contains(j))
      return 0;
  }

  return 1;
}

//#define DEBUG_ONLY_PROCESS_RINGS_WITH_ALL_AROMATIC_ATOMS

int
Molecule::_only_process_rings_with_all_aromatic_atoms(Kekule_Temporary_Arrays & kta)
{
  int * process_these_rings = kta.process_these_rings();

  const int * aromatic_atoms = kta.aromatic_atoms();

#ifdef DEBUG_ONLY_PROCESS_RINGS_WITH_ALL_AROMATIC_ATOMS
  iw_write_array(aromatic_atoms, _number_elements,
                 "_only_process_rings_with_all_aromatic_atoms:aromatic_atoms", cerr);
#endif

  int nr = _sssr_rings.number_elements();

  // Ring numbers that have only some of their atoms in `aromatic_atoms`.
  resizable_array<int> problematic_rings;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    int n = ri->count_members_set_in_array(aromatic_atoms, 1);

    //  cerr << "Ring " << i << ' ' << (*ri) << " N = " << n << endl;

    if (n == ri->number_elements())  // All ring members are aromatic
      continue;

    if (0 == n) {  // No ring atoms are aromatic.
      process_these_rings[i] = 0;
      continue;
    }

    //  if (ri->all_members_set_in_array(aromatic_atoms, 1))
    //    continue;

#ifdef DEBUG_ONLY_PROCESS_RINGS_WITH_ALL_AROMATIC_ATOMS
    cerr << "Adding problematic ring " << (*ri) << ", " << n << " of " << ri->number_elements()
         << " atoms aromatic\n";
#endif

    problematic_rings.add(i);
  }

  if (problematic_rings.empty())
    return nr;

  // For each of the potentially problematic rings, is there another ring
  // that covers all the aromatic atoms
  // In some molecule there is a 4 membered ring fused to the aromatic. Three of
  // the atoms in the 4 membered ring are aromatic. But all those atoms are in another
  // aromatic ring

  for (int i = problematic_rings.number_elements() - 1; i >= 0; i--)
  {
    //  cerr << "Ring " << problematic_rings[i] << " is potentially problematic\n";

    const Ring * pri = ringi(problematic_rings[i]);

    for (int j = 0; j < nr; j++)
    {
      if (problematic_rings.contains(j))
        continue;

      const Ring * rj = ringi(j);

#ifdef DEBUG_ONLY_PROCESS_RINGS_WITH_ALL_AROMATIC_ATOMS
      cerr << "Examining ring " << j << " fused? " << pri->is_fused_to(rj) << endl;
#endif

      if (! pri->is_fused_to(rj))
        continue;

      if (all_aromatic_atoms_in_first_contained_in_second(*pri, *rj, aromatic_atoms))
      {
        process_these_rings[problematic_rings[i]] = 0;
        problematic_rings.remove_item(i);
        break;
      }
    }
  }

  if (problematic_rings.empty())
    return 1;

#ifdef DEBUG_ONLY_PROCESS_RINGS_WITH_ALL_AROMATIC_ATOMS
  cerr << "Have " << problematic_rings.number_elements() << " rings that are problematic\n";
#endif

  // Now this gets ugly. Can we change the SSSR set so that we get rings
  // in there that include all the aromatic atoms

  for (int i = 0; i < _non_sssr_rings.number_elements(); i++)
  {
    const Ring * ri = _non_sssr_rings[i];

    if (! ri->all_members_set_in_array(aromatic_atoms, 1))
      continue;

    //  Can we have RI change places with one of the problematic ones

    for (int j = 0; j < problematic_rings.number_elements(); j++)
    {
      int k = problematic_rings[j];

      const Ring * rk = _sssr_rings[k];

#ifdef DEBUG_ONLY_PROCESS_RINGS_WITH_ALL_AROMATIC_ATOMS
      cerr << "Can it replace " << (*rk) << " comm " << rk->compute_bonds_shared_with(*ri) << endl;
#endif

      if (0 == rk->compute_bonds_shared_with(*ri))
        continue;

      _transfer_from_non_sssr_to_sssr_ring_set(i, k);
      problematic_rings.remove_item(j);
    }
  }

  if (problematic_rings.empty())
    return 1;

  for (int i = 0; i < problematic_rings.number_elements(); i++)
  {
    int j = problematic_rings[i];

    //  const Ring * rj = _sssr_rings[j];

    //  cerr << "Molecule::_only_process_rings_with_all_aromatic_atoms:ring " << (*rj) << " not all atoms aromatic, ignored\n";
    process_these_rings[j] = 0;
  }

  return 0;
}

int
Molecule::_has_aromatic_bonds(const Atom * a) const
{
  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_permanent_aromatic())
      return 1;
  }

  return 0;
}

/*
  When looking for kekule forms, we need to identify certain rings which
  will not be aromatic. The specific case which lead to this is:

                 \____/
                  |  |
          arom    |  | arom
                  |__|
                 /    \

  The inner 4 membered ring is not aromatic, even though all its atoms
  may have been specified as being aromatic. We use the rule that if
  it has no unshared pi electrons, it is aliphatic

  Well, it's not quite so simple. When I apply the rule above, I
  immediately ran into problems with something like

                      O
                      ||
                      ||
                     /  \
                   \/    \/
               arom |    | arom
                    |    |
                   / \  / \
                      \/ arom
                      ||
                      ||
                       O

  A strict interpretation of the above rule would classify this ring
  as non aromatic, which it is under Pearlman's rules. But we cannot
  exclude the ring, because there is an atom (the C in the Carbonyl)
  which is only in this ring. So, there is a kludge. If all atoms
  are in multiple rings, and it is a 4 membered ring, classify as
  non aromatic.
*/

int
Molecule::_kekule_identify_non_arom_rings(Kekule_Temporary_Arrays & kta)
{
  const int * aromatic_atoms = kta.aromatic_atoms();
  const int * aromatic_bonds = kta.aromatic_bonds();
  int * process_these_rings = kta.process_these_rings();

  int nr = _sssr_rings.number_elements();

//#define DEBUG_KEKULE_IDENTIFY_NON_AROM_RINGS
#ifdef DEBUG_KEKULE_IDENTIFY_NON_AROM_RINGS
  cerr << "At start of _kekule_identify_non_arom_rings\n";
  iw_write_array(process_these_rings, nr, "process_these_rings", cerr);
#endif

  for (int i = 0; i < nr; i++)
  {
    if (0 == process_these_rings[i])
      continue;

    const Ring * r = ringi(i);

    //  cerr << "ring i = " << i << " fused? " << r->is_fused() << " all_set " << r->all_members_set_in_array(aromatic_atoms, 1) << endl;

    if (! r->is_fused())    // only fused rings have this problem
      continue;

    if (! r->all_members_set_in_array(aromatic_atoms, 1))
      continue;

#ifdef DEBUG_KEKULE_IDENTIFY_NON_AROM_RINGS
    cerr << "What about ring " << i << " " << (*r) << endl;
#endif

    int n = r->number_elements();
    int found_atom_in_one_ring = 0;
    int found_atom_in_one_ring_with_pi_electrons = 0;
    for (int j = 0; j < n; j++)
    {
      atom_number_t k = r->item(j);
      if (nrings(k) > 1)    // atom shared with another ring.
        continue;

      found_atom_in_one_ring++;

      //    Crude test for number of pi electrons

      if (nbonds(k) < 4)
        found_atom_in_one_ring_with_pi_electrons++;
    }

    if (found_atom_in_one_ring)
      continue;

    if (0 == found_atom_in_one_ring_with_pi_electrons && 4 == n)
      process_these_rings[i] = 0;

#ifdef DEBUG_KEKULE_IDENTIFY_NON_AROM_RINGS
    if (0 == found_atom_in_one_ring_with_pi_electrons)
      cerr << "Non arom ring found:" << i << " " << *r << endl;
#endif
  }

  // If we have aromatic bonds we may be able to exclude rings that have non-aromatic bonds

  if (nullptr == aromatic_bonds)
    return 1;

  for (int i = 0; i < nr; i++)
  {
    if (! process_these_rings[i])
      continue;

    const Ring * ri = ringi(i);

    int ring_size = ri->number_elements();

    for (int j = 0; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);

      const Atom * ak = _things[k];

      if (6 != ak->atomic_number())
        continue;

      if (ak->nbonds() > ak->ncon())    // possibly aromatic
        continue;

      if (_has_aromatic_bonds(ak))
        continue;

      process_these_rings[i] = 0;
      //    cerr << "Turning off ring " << i << endl;
      break;
    }
  }

#ifdef DEBUG_KEKULE_IDENTIFY_NON_AROM_RINGS
  cerr << "After _kekule_identify_non_arom_rings\n";
  iw_write_array(process_these_rings, nr, "process_these_rings", cerr);
#endif

  return 1;
}

//#define DEBUG_GROUP_RINGS_INTO_GROUPS_TO_BE_PROCESSED_TOGETHER

/*
  For each set of rings that will need to be processed together,
  assign a unique number.

  Because the numbers in process_these_rings initially are 0 or 1, we
  start our numbering with 2
*/

int
Molecule::_group_rings_into_groups_to_be_processed_together(Kekule_Temporary_Arrays & kta)
{
  int * process_these_rings = kta.process_these_rings();

  int nr = _sssr_rings.number_elements();

#ifdef DEBUG_GROUP_RINGS_INTO_GROUPS_TO_BE_PROCESSED_TOGETHER
  cerr << "_group_rings_into_groups_to_be_processed_together\n";
  iw_write_array(process_these_rings, nr, "process_these_rings", cerr);
#endif

  int flag = 2;    // start numbering here because initially rings being processed will have value 1

  for (int i = 0; i < nr; i++)
  {
    if (0 == process_these_rings[i])    // ring not being processed
      continue;

    if (process_these_rings[i] >= 2)    // already been put into a ring or ring system
      continue;

    const Ring * ri = ringi(i);

    process_these_rings[i] = flag;

    if (ri->is_fused())
      _grow_fused_system(*ri, flag, process_these_rings);

    flag++;
  }

#ifdef DEBUG_GROUP_RINGS_INTO_GROUPS_TO_BE_PROCESSED_TOGETHER
  cerr << "After grouping\n";
  iw_write_array(process_these_rings, nr, "process_these_rings", cerr);

  for (int i = 2; i < flag; i++)
  {
    cerr << "Ring system with flag " << i;
    for (int j = 0; j < nr; j++)
    {
      if (process_these_rings[j] == i)
        cerr << ' ' << j;
    }
    cerr << endl;
  }
#endif

  return 1;
}

int
Molecule::_grow_fused_system(const Ring & r, int flag, int * in_system)
{
  for (int i = 0; i < r.fused_ring_neighbours(); i++)
  {
    const Ring * rn = r.fused_neighbour(i);

    int ring_number = rn->ring_number();

    if (0 == in_system[ring_number])    // not being processed at all
      continue;

    if (flag == in_system[ring_number])    // already been found to be part of this ring system
      continue;

    in_system[ring_number] = flag;
    _grow_fused_system(*rn, flag, in_system);
  }

  return 1;
}

int
Molecule::find_kekule_form(int * aromatic_atoms, const int * aromatic_bonds)
{
  if (nullptr == aromatic_atoms)
    ;
  else if (! perform_kekule_perception)
    return _process_without_kekule_perception(aromatic_atoms, aromatic_bonds);
  else
  {
    if (nullptr == _aromaticity)
      _aromaticity = new aromaticity_type_t[_number_elements];

    for (int i = 0; i < _number_elements; i++)
    {
      if (aromatic_atoms[i])
        _aromaticity[i] = AROMATIC;
      else
        _aromaticity[i] = NOT_AROMATIC;
    }
  }

  (void)ring_membership();    // ensure SSSR computed

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  iw_write_array(aromatic_atoms, _number_elements, "find_kekule_form:aromatic_atoms", cerr);
#endif

  int nb = _bond_list.number_elements();

  // first assign any aromatic bonds. Try to transfer that information to the atom array

  bond_type_t single_and_permanent_aromatic = (SINGLE_BOND | PERMANENT_AROMATIC_BOND);

  if (nullptr != aromatic_bonds)
  {
    for (int i = 0; i < nb; i++)
    {
      if (! aromatic_bonds[i])
        continue;

      const Bond * b = _bond_list[i];
      if (0 == b->nrings())
        continue;

      _bond_list[i]->set_bond_type(single_and_permanent_aromatic);    // does not perturb ring membership

      const atom_number_t a1 = b->a1();
      const atom_number_t a2 = b->a2();

      aromatic_atoms[a1] = 1;    // might already be set
      aromatic_atoms[a2] = 1;
    }
  }

  // Look for explicit hydrogen atoms bonded to the aromatic atoms,
  // and non periodic table atoms in the aromatic list.
  bool explicit_hydrogens_bonded_to_aromatics = false;
  int a_atoms_found = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom* a = _things[i];
    const atomic_number_t z = a->atomic_number();
    if (z == 6 || z == 7 || z == 8 || z == 9) {
      continue;
    }

    if (aromatic_atoms[i] && ! a->element()->is_in_periodic_table()) {
      ++a_atoms_found;
      continue;
    }

    if (1 != a->atomic_number())
      continue;

    if (1 != a->ncon())
      continue;

    const atom_number_t c = a->other(i, 0);

    if (aromatic_atoms[c]) {
      explicit_hydrogens_bonded_to_aromatics = true;
    }
  }

//#define ECHO_ATOMS_FOR_CONSIDERATION
#ifdef ECHO_ATOMS_FOR_CONSIDERATION
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " '" << smarts_equivalent_for_atom(i) << "' aromatic? "
         << aromatic_atoms[i] << endl;
  }
#endif

  int aromatic_atoms_found = 0;
  if (! _kekule_could_be_aromatic(aromatic_atoms, aromatic_atoms_found)) {
    cerr << "Could not be aromatic\n";
    return 0;
  }

  //cerr << "Any 'a' atoms found " << a_atoms_found << endl;
  //cerr << "Got any aromatic atoms " << aromatic_atoms_found << endl;

  if (0 == aromatic_atoms_found)    // no aromatic atoms
    return 1;

  int nr = _sssr_rings.number_elements();
  if (0 == nr)    // user must have asked to ignore chain aromatics
    return 1;

  Kekule_Temporary_Arrays kta(_number_elements, nr, aromatic_atoms, aromatic_bonds);

  if (nr > 6)
    kta.establish_smallest_ring_size(*this);

  _determine_hydrogen_counts(kta);

  kta.set_a_atoms_found(a_atoms_found);

  /*if (! _do_obvious_hcount_adjustments(kta))
  {
    cerr << "Molecule::find_kekule_form: bad hcount detected\n";
    return 0;
  }*/

  // Detaching Hydrogens will force a re-determination of rings. We depend on ring numbers,
  // so detach Hydrogens before doing anything else

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  if (explicit_hydrogens_bonded_to_aromatics)
    cerr << "There are explicit Hydrogens bonded to aromatic atoms\n";
#endif

  Temp_Detach_Atoms tda;

  if (explicit_hydrogens_bonded_to_aromatics) {
    tda.detach_atoms(*this);
    ring_membership();
  }

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  _only_process_rings_with_all_aromatic_atoms(kta);

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  if (_discard_non_aromatic_kekule_input)
    _kekule_suppress_non_aromatic_rings(kta);

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  // Jul 2009. The problem is that we have allowed these rings to participate
  // when making something aromatic, so if we do not allow them now, we get
  // cases where we cannot decode our own unique smiles. See if leaving this
  // out causes problems - like double bonds put along the unshared bond!
  //_kekule_identify_non_arom_rings (kta);

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  if (nullptr != aromatic_bonds)
    _kekule_check_rings_containing_aliphatic_bonds(kta);

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  _group_rings_into_groups_to_be_processed_together(kta);

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  // If we are going to re-try with positive nitrogen or whatever,
  // we need to make a copy of the aromatic_atoms array

  int rc = _find_kekule_form(kta);

#ifdef DEBUG_FIND_KEKULE_FORM_MAIN
  cerr << "at " << __LINE__ << ", rc = " << rc << endl;
  iw_write_array(kta.process_these_rings(), nr, "process_these_rings", cerr);
#endif

  if (rc)
    ;
  else if (kekule_try_positive_nitrogen || kekule_try_positive_oxygen)
  {
    rc = _find_kekule_form_positive_nitrogen(kta);
  }

  if (tda.natoms())
    tda.reattach_atoms(*this);

  _set_modified();

  // Any permanent aromatic atoms for which a kekule form was found can be unset

  nb = _bond_list.number_elements();    // may have changed

  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (! b->is_permanent_aromatic())
      continue;

    if (b->is_single_bond() || b->is_double_bond())    // great, kekule form found
      b->set_permanent_aromatic(0);
  }

  return rc;
}

/*
  We failed to find a Kekule form for the molecule as read in. Can we
  put a positive charge on a Nitrogen atom in order to achieve aromaticity.

  We look for 3 connected Nitrogen atoms.

  Note that this function doesn't work if there are two such atoms in the
  molecule...

  AROMATIC_ATOMS contains the identity of the atoms that are known to
  be aromatic but haven't yet been associated with an aromatic ring
*/

//#define DEBUG_VARY_NPLUS

int
Molecule::_find_kekule_form_positive_nitrogen(Kekule_Temporary_Arrays & tka)
{
  int * aromatic_atoms = tka.aromatic_atoms();
  int * vary_bonds = tka.vary_bonds();

  int atoms_to_process = 0;

#ifdef DEBUG_VARY_NPLUS
  cerr << "Any positive Nitrogen candidates?\n";
#endif

  for (int i = 0; i < _number_elements; i++)
  {
#ifdef DEBUG_VARY_NPLUS
    cerr << "Atom " << i << " aromatic " << aromatic_atoms[i] << " vary bonds " << vary_bonds[i]
         << " vary hcount " << tka.vary_hcount()[i] << endl;
#endif

    if (aromatic_atoms[i])
      atoms_to_process++;
  }

  assert(atoms_to_process > 0);

#ifdef DEBUG_VARY_NPLUS
  cerr << "Let's see if we can place a positive charge\n";
#endif

  Set_of_Atoms possible_nplus;

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == aromatic_atoms[i])
      continue;

    Atom * a = _things[i];

    if (0 != a->formal_charge())
      continue;

    if (7 == a->atomic_number())
      ;
    else if (15 == a->atomic_number())
      ;
    else
      continue;

    if (3 != a->ncon())
      continue;

    if (4 == a->nbonds())    // things like *=N(=O)-*, which should have been dealt with previously
      continue;

    //  Make sure we first process those N atoms in multiple rings.

    if (possible_nplus.empty())
      possible_nplus.add(i);
    else if (1 == nrings(i))
      possible_nplus.add(i);
    else
      possible_nplus.insert_before(0, i);

    vary_bonds[i] = 1;

#ifdef DEBUG_VARY_NPLUS
    cerr << "Atom " << i << " is a possible positively charged " << a->atomic_symbol() << "\n";
#endif
  }

  // Now look for possible positively charged oxygens

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == aromatic_atoms[i])
      continue;

    Atom * a = _things[i];

    if (0 != a->formal_charge())
      continue;

    if (8 == a->atomic_number())
      ;
    else if (16 == a->atomic_number())
      ;
    else
      continue;

    //  if (2 != a->ncon())    not needed, otherwise cannot read Co1cccn1
    //      continue;

    //  Make sure we first process those N atoms in multiple rings.

    if (possible_nplus.empty())
      possible_nplus.add(i);
    else if (1 == nrings(i))
      possible_nplus.add(i);
    else
      possible_nplus.insert_before(0, i);

#ifdef DEBUG_VARY_NPLUS
    cerr << "Atom " << i << " is a possible positively charged Oxygen\n";
#endif
  }

  if (possible_nplus.empty())
    return 0;

  return _find_kekule_form_positive_nitrogen(possible_nplus, tka);
}

void
Molecule::_set_all_bonds_in_system_to_single(const int * process_these_atoms)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    const Atom * a = _things[i];

    int acon = a->ncon();

    if (acon == a->nbonds())
      continue;

    for (int j = 0; j < acon; j++)
    {
      Bond * b = a->item(j);

      if (b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (0 == process_these_atoms[k])
        continue;

      b->set_bond_type(SINGLE_BOND);
    }
  }

  return;
}

/*
  We have identified some Nitrogen atoms that could possibly have a
  positive charge put on them.
  Our strategy is to first try to make the rings aromatic by considering
  one ring at a time. If that doesn't work, then try combining rings
*/

int
Molecule::_find_kekule_form_positive_nitrogen(Set_of_Atoms & possible_nplus,
                                              Kekule_Temporary_Arrays & kta)
{
  int * process_these_rings = kta.process_these_rings();

#ifdef DEBUG_VARY_NPLUS
  for (int i = 0; i < nrings(); i++)
  {
    cerr << "Processing ring " << i << " " << process_these_rings[i] << endl;
  }
#endif

  resizable_array<const Ring *> rings;    // scope here for efficiency

  int rc = 1;

  // We need to process all ring groupings individually.

  // There may have been some partial results previously, so we save the existing bonds in case of failure

  bond_type_t * bsave = new bond_type_t[_bond_list.number_elements()];
  std::unique_ptr<bond_type_t[]> free_bsave(bsave);

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    if (0 == process_these_rings[i])
      continue;

    int flag = process_these_rings[i];

    _find_next_ring_or_ring_system(kta, flag, rings);

    _bond_list.copy_bond_types(bsave);

    _set_all_bonds_in_system_to_single(kta.process_these_atoms());

#ifdef DEBUG_VARY_NPLUS
    cerr << "Processing ring system with " << rings.number_elements() << " rings\n";
#endif

    Set_of_Atoms possible_nplus_in_system;

    _identify_possible_nplus_in_system(possible_nplus, possible_nplus_in_system,
                                       kta.process_these_atoms());

    if (_find_kekule_form_positive_nitrogen_ring_system(kta, rings, possible_nplus_in_system))
      unset_all_atoms_in_rings(rings, kta.aromatic_atoms(),
                               process_these_rings);    // mark atoms and rings as done
    else
    {
      _restore_existing_bond_types(bsave);
      rc = 0;
    }

    if (possible_nplus.empty())
      break;
  }

  return rc;
}

/*
  We have a set of fused rings and a set of atoms in those rings. By placing positive
  charges on the Nitrogens, try to find a Kekule form
*/

int
Molecule::_find_kekule_form_positive_nitrogen_ring_system(
    Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings,
    const Set_of_Atoms & possible_nplus_in_system)
{
  int * process_these_atoms = kta.process_these_atoms();

  for (int i = 0; i < _number_elements; i++)
  {
    if (process_these_atoms[i])
      process_these_atoms[i] = KEKULE_READY_TO_PROCESS;
  }

  int nn = possible_nplus_in_system.number_elements();

  for (int i = 1; i <= nn; i++)    // do the nitrogens 1 at a time, 2 at a time, 3 at a time ...
  {
    Set_of_Atoms tmp;

    if (_find_kekule_form_positive_nitrogen_ring_system(kta, rings, possible_nplus_in_system, tmp,
                                                        i))
      return 1;
  }

  return 0;
}

int
Molecule::_find_kekule_form_positive_nitrogen_ring_system(
    Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings,
    const Set_of_Atoms & possible_nplus_in_system, Set_of_Atoms & s, int maxdepth)
{
  int nn = possible_nplus_in_system.number_elements();

  for (int i = 0; i < nn; i++)
  {
    atom_number_t n = possible_nplus_in_system[i];

    if (s.contains(n))
      continue;

    s.add(n);

    if (_find_kekule_form_positive_nitrogen_ring_system_place_charges(kta, rings, s))
      return 1;

    if (s.number_elements() == maxdepth)    // at max depth, cannot continue
      ;
    else if (_find_kekule_form_positive_nitrogen_ring_system(kta, rings, possible_nplus_in_system,
                                                             s, maxdepth))
      return 1;

    s.chop();
  }

  return 0;
}

/*
  Mojor deficiency. Charges are applied one at a time. Should allow multiple charges
*/

int
Molecule::_find_kekule_form_positive_nitrogen_ring_system_place_charges(
    Kekule_Temporary_Arrays & kta, const resizable_array<const Ring *> & rings,
    const Set_of_Atoms & possible_nplus_in_system)
{
  int nn = possible_nplus_in_system.number_elements();
  int * vary_bonds = kta.vary_bonds();

  for (int i = 0; i < nn; i++)
  {
    atom_number_t n = possible_nplus_in_system[i];

    _things[n]->set_formal_charge(1);
    _things[n]->set_implicit_hydrogens(0);
    vary_bonds[n] = 1;    // maybe we should save the original value and restore it later

#ifdef DEBUG_VARY_NPLUS
    cerr << "Just put formal charge on atom " << n << endl;
#endif
  }

  int pi_electron_count = 0;

  int rc = _find_kekule_form(rings, kta, possible_nplus_in_system[0], pi_electron_count);

#ifdef DEBUG_VARY_NPLUS
  cerr << "Result from _find_kekule_form " << rc << endl;
#endif

  if (rc)
  {
    for (int i = 0; i < rings.number_elements(); i++)
    {
      const Ring * ri = rings[i];

      ri->set_vector(kta.aromatic_atoms(), 0);
    }
  }
  else
  {
    for (int i = 0; i < nn; i++)
    {
      atom_number_t n = possible_nplus_in_system[i];

      _things[n]->set_formal_charge(0);
      _things[n]->set_modified();
      vary_bonds[n] = 0;

#ifdef DEBUG_VARY_NPLUS
      cerr << "Removed formal charge from atom " << n << endl;
#endif
    }
  }

  return rc;
}

/*
  A molecule was input which had either no rings or some of its
  aromatic atoms were not in a ring. Convert those bonds

  RM is the same as RING_MEMBERSHIP
*/

int
Molecule::_convert_chain_aromatic_bonds(int * aromatic_atoms, const int * rm)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (! aromatic_atoms[i])
      continue;

    if (rm[i] > 0)    // in a ring, we only process chain atoms
      continue;

    const Atom * a = _things[i];
    int icon = a->ncon();

    Set_of_Atoms aromatic_atoms_found;

    for (int j = 0; j < icon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (aromatic_atoms[k])
        aromatic_atoms_found.add(k);
    }

    int na = aromatic_atoms_found.number_elements();

    if (0 == na)
      cerr << "Molecule::_convert_chain_aromatic_bonds: atom " << i
           << " was aromatic but no connections were??\n";
    else if (1 == na)
    {
      atom_number_t j = aromatic_atoms_found[0];

      if (j < i)
        continue;

      if (_things[i]->implicit_hydrogens() && _things[j]->implicit_hydrogens())
        set_bond_type_between_atoms(i, j, DOUBLE_BOND);
      else
        set_bond_type_between_atoms(i, j, SINGLE_BOND);    // that should already be what it is
      //      cerr << "Molecule::_convert_chain_aromatic_bonds: chain atoms " << i << " and " << j << " were aromatic\n";
      aromatic_atoms[i] = aromatic_atoms[j] = 0;
    }
    else
    {
      cerr << "Molecule::_convert_chain_aromatic_bonds:found " << aromatic_atoms_found
           << " aromatic atoms attached\n";
    }
  }

  return 1;
}

/*
  We want to figure out if an atom is in a ring that is surrounded by aliphatic rings. If so,
  that ring should be processed as a single ring. All we do is make sure that there is only
  one active ring that contains the atom
*/

const Ring *
Molecule::_single_active_ring_containing_atom(atom_number_t n, const int * process_these_rings)
{
  const Ring * rc = nullptr;

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    if (! process_these_rings[i])
      continue;

    const Ring * ri = ringi(i);

    if (! ri->contains(n))
      continue;

    if (nullptr != rc)
      return nullptr;

    rc = ri;
  }

  return rc;
}

/*
  Which of the atoms in POSSIBLE_NPLUS are set in PROCESS_THESE_ATOMS
*/

int
Molecule::_identify_possible_nplus_in_system(Set_of_Atoms & possible_nplus,
                                             Set_of_Atoms & possible_nplus_in_system,
                                             const int * process_these_atoms) const
{
  possible_nplus_in_system.resize_keep_storage(0);

  int nn = possible_nplus.number_elements();

  for (int i = nn - 1; i >= 0; i--)
  {
    atom_number_t n = possible_nplus[i];

    if (process_these_atoms[n])
    {
      possible_nplus.remove_item(i);
      possible_nplus_in_system.add(n);
    }
  }

  return possible_nplus_in_system.number_elements();
}

int
Molecule::_restore_existing_bond_types(const bond_type_t * bsave)
{
  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (BOND_TYPE_ONLY(b->btype()) == bsave[i])
      continue;

    b->set_bond_type(bsave[i]);
    _things[b->a1()]->recompute_nbonds();
    _things[b->a2()]->recompute_nbonds();
  }

  return 1;
}

/*
  First check the elements. If they are all permanent aromatic, then set the atoms
  to permanent aromatic
*/

int
Molecule::_maybe_all_atoms_are_permanent_aromatic(Kekule_Temporary_Arrays & kta,
                                                  const resizable_array<const Ring *> & rings)
{
  int nr = rings.number_elements();

  int all_aromatic = 1;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = rings[i];

    int ring_size = ri->number_elements();

    for (int j = 0; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);

      if (_things[k]->element()->permanent_aromatic())
        continue;

      all_aromatic = 0;
      break;
    }

    if (! all_aromatic)
      break;
  }

  if (all_aromatic)
    return _kekule_all_atoms_and_bonds_in_system_aromatic(kta, rings);

  if (non_kekule_systems_ok_to_be_aromatic())
    return _kekule_all_atoms_and_bonds_in_system_aromatic(kta, rings);

  return 0;
}

int
Molecule::unset_all_permanent_aromatic_atoms()
{
  assert(ok());

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_permanent_aromatic(0);
  }

  return 1;
}

int
Molecule::change_double_bonds_between_permanent_aromatic_to_single(int call_set_modified)
{
  int rc = 0;

  int nb = _bond_list.number_elements();

  Bond ** b = const_cast<Bond **>(_bond_list.rawdata());

  for (int i = 0; i < nb; i++)
  {
    Bond * bi = b[i];

    if (! bi->is_double_bond())
      continue;

    if (! _things[bi->a1()]->permanent_aromatic())
      continue;

    if (! _things[bi->a2()]->permanent_aromatic())
      continue;

    bi->set_bond_type(SINGLE_BOND | AROMATIC_BOND);
    bi->set_permanent_aromatic(1);
    rc++;
  }

  if (0 == rc)
    ;
  else if (call_set_modified)
    _set_modified();

  return rc;
}

//  Do not do Kekule perception, but set aromatic elements to
//  permanent aromatic forms.
int
Molecule::_process_without_kekule_perception(const int * aromatic_atoms, const int * aromatic_bonds)
{
  if (nullptr == aromatic_carbon)
    create_permanent_aromatic_elements();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == aromatic_atoms[i])
      continue;

    Atom * a = _things[i];

    const Element * e = a->element();

    if (e->permanent_aromatic())
      continue;

    if (6 == e->atomic_number())
      a->set_element(aromatic_carbon);
    else if (7 == e->atomic_number())
      a->set_element(aromatic_nitrogen);
    else if (8 == e->atomic_number())
      a->set_element(aromatic_oxygen);
    else if (16 == e->atomic_number())
      a->set_element(aromatic_sulphur);
  }

  // Should process aromatic bonds too...

  return 1;
}

/*int
Molecule::change_ring_bonds_between_permanent_aromatic_to_aromatic (int call_set_modified)
{
  if (0 == nrings())
    return 0;

  int rc = 0;

  ring_membership();

  int nb = _bond_list.number_elements();

  Bond ** b = const_cast<Bond **>(_bond_list.rawdata());

  for (int i = 0; i < nb; i++)
  {
    Bond * bi = b[i];

    if ( 0 == bi->nrings())
      continue;

    if ( !_things[bi->a1()]->permanent_aromatic() )
      continue;

    if (! _things[bi->a2()]->permanent_aromatic())
      continue;

    bond_type_t bnd_type= bi->btype();
    bi->set_bond_type( bnd_type | AROMATIC_BOND | PERMANENT_AROMATIC_BOND );
    bi->set_permanent_aromatic(1);
    rc++;
  }

  if (0 == rc)
    ;
  else if (call_set_modified)
    _set_modified();

  return rc;
}*/

/*
  helper function for change_ring_bonds_between_permanent_aromatic_to_aromatic

  We have two atoms at either ends of a bond. Are they part of a ring for
  which every member of the ring is marked permanent aromatic
*/

int
Molecule::_bond_in_ring_with_all_members_permanent_aromatic(atom_number_t a1, atom_number_t a2,
                                                            const int * is_permanent_aromatic) const
{
  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    //  if (ri->number_elements() <= 4)    uncomment for no 4 membered aromatics, but beware O=C1C(=O)C=C1 ???
    //    continue;

    if (! ri->contains_bond(a1, a2))
      continue;

    if (ri->all_members_non_zero_in_array(is_permanent_aromatic))
      return 1;
  }

  return 0;
}

int
Molecule::change_ring_bonds_between_permanent_aromatic_to_aromatic(int call_set_modified)
{
  if (0 == nrings())
    return 0;

  ring_membership();    // force sssr

  int * is_permanent_aromatic = new int[_number_elements];
  std::unique_ptr<int[]> free_is_permanent_aromatic(is_permanent_aromatic);

  int permanent_aromatic_count = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->permanent_aromatic())
    {
      is_permanent_aromatic[i] = 1;
      permanent_aromatic_count++;
    }
    else
      is_permanent_aromatic[i] = 0;
  }

  if (0 == permanent_aromatic_count)    // nothing to process
    return 0;

  int rc = 0;

  int nb = _bond_list.number_elements();

  Bond ** b = const_cast<Bond **>(_bond_list.rawdata());

  for (int i = 0; i < nb; i++)
  {
    Bond * bi = b[i];

    if (bi->is_permanent_aromatic())
      continue;

    if (! is_permanent_aromatic[bi->a1()])
      continue;

    if (! is_permanent_aromatic[bi->a2()])
      continue;

    if (0 == bi->nrings())
      continue;

    if (bi->is_aromatic())
      ;
    else if (_bond_in_ring_with_all_members_permanent_aromatic(bi->a1(), bi->a2(),
                                                               is_permanent_aromatic))
      continue;

    bond_type_t bnd_type = bi->btype();
    bi->set_bond_type(bnd_type | AROMATIC_BOND | PERMANENT_AROMATIC_BOND);
    bi->set_permanent_aromatic(1);
    rc++;
  }

  if (0 == rc)
    ;
  else if (call_set_modified)
    _set_modified();

  return rc;
}

/*
  This is related to Valhalla substructure searching where we want some
  kinds of aromatic bonds to match both aromatic forms and the underlying
  Kekule bonds. So, anything that is NOT a fixed Kekule form, we lose
  the Kekule forms.
  Note that we leave the molecule in an awkward state, and don't call
  set_modified

  To deal with fused rings, where one ring might be an alternating
  Kekule form, but the other fixed, most 6-5 systems, we first
  identify the bonds that will retain their characteritics
*/

//#define DEBUG_SET_BOND_TYPES_FOR_ISIS_AROMATICITY_MATCHING

int
Molecule::set_bond_types_for_isis_aromaticity_matching()
{
  compute_aromaticity_if_needed();

  int nr = nrings();

#ifdef DEBUG_SET_BOND_TYPES_FOR_ISIS_AROMATICITY_MATCHING
  cerr << "set_bond_types_for_isis_aromaticity_matching:molecule has " << nr << " rings\n";
#endif

  int rc = 0;

  int * keep_as_is = new_int(_number_elements * _number_elements);
  std::unique_ptr<int[]> free_keep_as_is(keep_as_is);

  // First mark bonds in the fixed kekule form rings as unchanging

  int alternating_kekule_form_rings_present = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    if (! ri->is_aromatic())
      continue;

    if (! is_fixed_kekule_form(*this, *ri))
    {
      alternating_kekule_form_rings_present++;
      continue;
    }

    for (Ring_Bond_Iterator j(*ri); j != ri->zend(); j++)
    {
      atom_number_t a1 = j.a1();
      atom_number_t a2 = j.a2();

      keep_as_is[a1 * _number_elements + a2] = 1;
      keep_as_is[a2 * _number_elements + a1] = 1;
    }
  }

  if (0 == alternating_kekule_form_rings_present)
    return 0;

  // Now set all the bonds in the alternating rings to aromatic

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    if (! ri->is_aromatic())
      continue;

    if (is_fixed_kekule_form(*this, *ri))
      continue;

#ifdef DEBUG_SET_BOND_TYPES_FOR_ISIS_AROMATICITY_MATCHING
    cerr << "Ring " << (*ri) << " is alternating aromatic\n";
#endif

    for (Ring_Bond_Iterator j(*ri); j != ri->zend(); j++)
    {
      atom_number_t a1 = j.a1();
      atom_number_t a2 = j.a2();

      if (keep_as_is[a1 * _number_elements + a2])
        continue;

      Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));

      Connection * c = b;

#ifdef DEBUG_SET_BOND_TYPES_FOR_ISIS_AROMATICITY_MATCHING
      cerr << "Set bond btw " << a1 << " and " << a2 << endl;
#endif

      c->set_bond_type(
          AROMATIC_BOND);    // use the underlying Connection method so we avoid an OK_BOND_TYPE check
    }

    rc++;
  }

  return rc;
}

void
reset_aromatic_file_scope_variables()
{
  record_dual_nature_aromaticity = 0;
  file_scope_display_no_kekule_form_message = 1;
  max_aromatic_ring_size = 8;
  min_aromatic_ring_size = 4;
  perform_kekule_perception = 1;
  allow_pipeline_pilot_aromaticity_on_input = 1;
  _allow_two_electron_systems_to_be_aromatic = 0;
  x_convert_chain_aromatic_bonds = 0;
  _aromatic_chain_bonds_are_ok = 0;
  _non_kekule_systems_ok_to_be_aromatic = 0;
  global_aromaticity_determination_type = Daylight;
  all_bonds_in_aromatic_ring_must_be_aromatic = 1;
  warn_aromatic_chain_atoms = 1;
  kekule_try_positive_nitrogen = 0;
  kekule_try_positive_oxygen = 0;
  _input_aromatic_structures = 1;
  _allow_input_without_valid_kekule_form = 0;
  _allow_delocalised_carbonyl_bonds = 0;
  _discard_non_aromatic_kekule_input = 0;

  aromatic_carbon = nullptr;
  aromatic_nitrogen = nullptr;
  aromatic_oxygen = nullptr;
  aromatic_phosphorus = nullptr;
  aromatic_sulphur = nullptr;
  aromatic_a = nullptr;

  return;
}

static int
identify_nitrogen_and_adjacent(const Molecule & m, const Set_of_Atoms & r, atom_number_t & n,
                               atom_number_t & adjacent, bond_type_t & bt)
{
  const int ring_size = r.size();

  for (int i = 0; i < ring_size; ++i)
  {
    const auto j = r[i];

    if (7 != m.atomic_number(j))
      continue;

    n = j;
    adjacent = r.next_after_wrap(i);

    bt = m.btype_between_atoms(j, adjacent);
    return 1;
  }

  return 0;
}

/*
  We do not attempt to generate all switched forms. We just cycle through
  the rings and process them one at a time. Obviously more combinations
  could be possible by making this recursive...
*/

int
Molecule::generate_switched_kekule_forms(resizable_array_p<Molecule> & variant)
{
  const auto nr = nrings();

  if (0 == nr)
    return 0;

  IW_STL_Hash_Set already_done;

  compute_aromaticity_if_needed();

  already_done.insert(smiles());    // note we do not use the unique smiles

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = ringi(i);

    if (! ri->is_aromatic())
      continue;

    atom_number_t n, adjacent;
    bond_type_t bt;
    if (! identify_nitrogen_and_adjacent(*this, *ri, n, adjacent, bt))
      continue;

    Toggle_Kekule_Form tkf;
    tkf.set_display_error_messages(0);

    if (IS_SINGLE_BOND(bt))
      tkf.add_bond(0, 1, DOUBLE_BOND);
    else
      tkf.add_bond(0, 1, SINGLE_BOND);

    Set_of_Atoms s;
    s.add(n);
    s.add(adjacent);

    Molecule * mcopy = new Molecule(*this);
    int changed = 0;
    if (! tkf.process(*mcopy, s, changed) || ! changed || already_done.contains(mcopy->smiles()))
    {
      delete mcopy;
      continue;
    }

    already_done.insert(mcopy->smiles());

    variant.add(mcopy);
  }

  return variant.number_elements();
}

/*
  We need to know if an atom is in rings that are likely to have alternating forms.
  Note that this is not a robust check,
*/

int
Molecule::all_rings_containing_atom_are_kekule(const atom_number_t zatom)
{
  compute_aromaticity_if_needed();

  const int nr = nrings();

  int found_in_ring = 0;

  //cerr << "Molecule contains " << nr << " rings\n";

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (! ri->contains(zatom))
      continue;

    found_in_ring++;

    if (! _has_kekule_forms(*ri))
      return 0;
  }

  if (0 == found_in_ring)    // should complain
    return 0;

  return 1;
}

int
Molecule::_has_kekule_forms(const Set_of_Atoms & r) const
{
  for (auto i : r)
  {
    //  cerr << "Checking atom " << i << ' ' << smarts_equivalent_for_atom(i) << endl;

    const Atom * a = _things[i];

    const auto z = a->atomic_number();

    if (8 == z || 16 == z)
      return 0;

    if (7 == z && 3 == a->ncon() && 0 == a->formal_charge())
      return 0;

    if (_doubly_bonded_to_something_outside_ring(i))
      return 0;
  }

  return 1;
}

void add_aromatic(aromaticity_type_t & arom) {
  arom |= _AROM_BIT;
}

void add_aliphatic(aromaticity_type_t & arom) {
  arom |= _ALIPH_BIT;
}

bool is_aromatic_atom(aromaticity_type_t arom) {
  return (_AROM_BIT & arom) != 0;
}

bool is_aliphatic_atom(aromaticity_type_t arom) {
  return (_ALIPH_BIT & arom) != 0; 
}
