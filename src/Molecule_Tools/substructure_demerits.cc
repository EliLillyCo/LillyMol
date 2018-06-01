#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "misc.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring_data_source.h"

#include "path.h"
#include "demerit.h"
#include "charge_assigner.h"

#include "qry_and_demerit.h"

namespace substructure_demerits
{
static int verbose = 0;

void
set_verbose(int s)
{
  verbose = s;
}

static int keep_going_after_rejection = 0;

void
set_keep_going_after_rejection(int s)
{
  keep_going_after_rejection = s;
}

static Charge_Assigner _charge_assigner;

int
initialise_charge_assigner(const char * s)
{
  return _charge_assigner.build(s);
}

static int too_many_rings_cutoff = 6;

void
set_substructure_demerits_too_many_rings(int s)
{
  too_many_rings_cutoff = s;

  return;
}

static int ring_size_too_large = 7;

void
set_substructure_demerits_ring_size_too_large(int s)
{
  ring_size_too_large = s;

  return;
}

static int all_demerits_same_numeric_value = -1;

static int
demerit_or_default_if_specified(int d)
{
  if (all_demerits_same_numeric_value >= 0)
    return all_demerits_same_numeric_value;

  return d;
}

void
set_all_numeric_demerit_values(int s)
{
  all_demerits_same_numeric_value = s;

  return;
}

/*
  We need to initialise the charge assigner from the command line, but this file
  never sees the command line. Therefore we can pass a reference to our charge
  assigner to the controlling programme.
  Pretty awful...
*/

extern Charge_Assigner &
charge_assigner()
{
  return _charge_assigner;
}

static int
identify_largest_fragment(int matoms,
                           int nf,
                           const int * fragment_membership)
{
  int largest_fragment = -1;
  int atoms_in_largest_fragment = -1;

  for (int i = 0; i < nf; i++)
  {
    int tmp = count_occurrences_of_item_in_array(i, matoms, fragment_membership);

    if (tmp > atoms_in_largest_fragment)
    {
      largest_fragment = i;
      atoms_in_largest_fragment = tmp;
    }
  }

  return largest_fragment;
}

static int apply_positive_charge_demerit = 1;
static int apply_negative_charge_demerit = 1;

static int
too_many_charges(Molecule & m, Demerit & demerit)
{
  if (0 == apply_positive_charge_demerit && 0 == apply_negative_charge_demerit)
    return 0;

  int matoms = m.natoms();

  formal_charge_t * f = new formal_charge_t[matoms]; std::unique_ptr<formal_charge_t[]> free_f(f);

  set_vector(f, matoms, static_cast<formal_charge_t>(0));

  _charge_assigner.set_apply_charges_to_molecule(0);

  if (0 == _charge_assigner.process(m, f))
    return 0;

  int * fragment_membership = new int[matoms]; std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  int nf = m.fragment_membership(fragment_membership);

  int largest_fragment;

  if (1 == nf)
    largest_fragment = 0;
  else
    largest_fragment = identify_largest_fragment(matoms, nf, fragment_membership);

  int nneg = 0;
  int npos = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (fragment_membership[i] != largest_fragment)
      continue;

    if (0 == f[i])
      continue;

    if (f[i] < 0)
      nneg++;
    else 
      npos++;
  }

//cerr << m.name() << " nneg " << nneg << " npos " << npos << endl;

  int d = demerit_or_default_if_specified(50);

  if (nneg > 1)
    demerit.extra((nneg - 1) * d, "negative");

  if (npos > 1)
    demerit.extra((npos - 1) * d, "positive");

  return demerit.rejected();
}

static int
compute_max_nrings(Molecule & m,
                    const Ring & r)
{
  int rc = 0;

  for (Ring_Bond_Iterator i(r); i != r.end(); ++i)
  {
    const Bond * b = m.bond_between_atoms(i.a1(), i.a2());

    int rm = b->nrings();

    if (rm > rc)
      rc = rm;
  }

  return rc;
}

static int complex_fused_rings_count = 0;

static int
complex_fused_rings(Molecule & m)
{
  int nr = m.nrings();

  if (nr < 3)
    return 0;

  (void) m.ring_membership();

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  iwmax<int> max_nrings_for_a_bond(0);
  iwmin<int> min_fused_neighbours(nr);

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    if (0 == ri->fused_ring_neighbours())   // isolated ring
      continue;

    min_fused_neighbours.extra(ri->fused_ring_neighbours());

    int rings_in_system = 1;

    int fsid = ri->fused_system_identifier();

    ring_already_done[i] = fsid;

    max_nrings_for_a_bond.extra(compute_max_nrings(m, *ri));

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi(j);

      if (rj->fused_system_identifier() != fsid)
        continue;

      min_fused_neighbours.extra(rj->fused_ring_neighbours());

      ring_already_done[j] = fsid;

      rings_in_system++;

      max_nrings_for_a_bond.extra(compute_max_nrings(m, *rj));
    }

    if (max_nrings_for_a_bond.maxval() >= 3)
      return 1;
    else if (rings_in_system > 3 && min_fused_neighbours.minval() >= 2)
      return 1;
  }

  return 0;
}

static int
complex_fused_rings(Molecule & m,
                     Demerit & demerit)
{
  if (complex_fused_rings(m))
  {
    demerit.reject("complexfusedrings");
    complex_fused_rings_count++;
    return 1;
  }

  return 0;
}

static int
identify_largest_saturated_carbon_section(Molecule & m,
                                           const atomic_number_t * z,
                                           const int * ncon,
                                           const int * attached_heteroatom_count,
                                           atom_number_t zatom,
                                           int flag,
                                           int * already_done)
{
  int rc = 1;

  already_done[zatom] = flag;

  const Atom * a = m.atomi(zatom);

  for (int i = 0; i < ncon[zatom]; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (already_done[j])
      continue;

    if (attached_heteroatom_count[j])
      continue;

    if (6 != z[j])
      continue;

    if (ncon[j] < m.nbonds(j))
      continue;

    rc += identify_largest_saturated_carbon_section(m, z, ncon, attached_heteroatom_count, j, flag, already_done);
  }

  return rc;
}

static int satcg_count = 0;
static int apply_satcg = 0;   // turn off until we get this debugged

static int
large_saturated_carbon_sections_including_rings(Molecule & m,
                             Demerit & demerit,
                             const atomic_number_t * z,
                             const int * ncon,
                             const int * ring_membership)
{
  int matoms = m.natoms();

  int * attached_heteroatom_count = new_int(matoms); std::unique_ptr<int[]> free_attached_heteroatom_count(attached_heteroatom_count);

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (6 != z[k])
        attached_heteroatom_count[i] = 1;
      else if (m.is_aromatic(k))
        attached_heteroatom_count[i] = 1;
      else if (m.multiple_bond_to_heteroatom(k))
        attached_heteroatom_count[i] = 1;

      if (attached_heteroatom_count[i])
        break;
    }
  }

  int * tmp = new_int(matoms); std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < matoms; i++)
  {
    if (tmp[i])   // already been discerned
      continue;

    if (attached_heteroatom_count[i])
      continue;

    if (6 != z[i])
      continue;

    if (ncon[i] < m.nbonds(i))
      continue;

    int flag = i + 1;    // a unique identifier for this grouping

    int carbon_atoms_in_group = identify_largest_saturated_carbon_section(m, z, ncon, attached_heteroatom_count, i, flag, tmp);

    if (carbon_atoms_in_group < 7)
      continue;

//  Must be at least 5 non-ring atoms in group

    int non_ring_atoms_in_grouping = 0;
    int nch2 = 0;
    int atoms_in_multiple_rings = 0;

    for (int j = 0; j < matoms; j++)
    {
      if (flag != tmp[j])
        continue;

      if (2 == ncon[j])
        nch2++;

      if (0 == ring_membership[j])
        non_ring_atoms_in_grouping++;
      else if (ring_membership[j] > 1)
        atoms_in_multiple_rings++;
    }


//#define DEBUG_LARGE_SATURATED_CARBON_SECTIONS_INCLUDING_RINGS
#ifdef DEBUG_LARGE_SATURATED_CARBON_SECTIONS_INCLUDING_RINGS
    cerr << "carbon_atoms_in_group " << carbon_atoms_in_group << endl;
    cerr << "ring_atoms_in_grouping " << (carbon_atoms_in_group - non_ring_atoms_in_grouping) <<  endl;
    cerr << "atoms_in_multiple_rings " << atoms_in_multiple_rings << endl;
    cerr << "nch2 " << nch2 << endl;
#endif

    if (nch2 < 5)
      continue;

    if (atoms_in_multiple_rings > 0)
      continue;

    int ring_atoms_in_grouping = carbon_atoms_in_group - non_ring_atoms_in_grouping;

    if (ring_atoms_in_grouping >= 14)   // don't want to hit steroids
      continue;

    if (non_ring_atoms_in_grouping >= 5)
    {
      demerit.reject("SATCG");
      satcg_count++;
      return 1;
    }
  }

  return 0;
}

static int
bonded_to_heteroatom_outside_ring(Molecule & m,
                                   const Ring & r,
                                   const Ring_Atom_Iterator & rai,
                                   const atomic_number_t * z)
{
  atom_number_t zatom = rai.current();

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (rai.is_next_or_previous(j))
      continue;

    if (6 != z[j])
      rc++;
  }
  
  return rc;
}

static int
is_rejectable_large_ring(Molecule & m,
                          const Ring & ri,
                          const atomic_number_t * z,
                          const int * ncon)
{
  int heteroatoms_encountered = 0;
  int fused_ring_atoms_encountered = 0;
  int heteroatom_outside_ring = 0;

  for (Ring_Atom_Iterator i(ri); i != ri.end(); i++)
  {
    atom_number_t j = i.current();

    if (6 != z[j])
//    heteroatoms_encountered++;    relaxed form, not used
      return 0;

    const Atom * aj = m.atomi(j);

    if (aj->unsaturated())    // what about -C(=*)- too rare to worry about
      return 0;

    if (2 == ncon[j])
      continue;

    if (m.nrings(j) > 1)
      fused_ring_atoms_encountered++;

    if (bonded_to_heteroatom_outside_ring(m, ri, i, z))
      heteroatom_outside_ring++;
  }

  if (heteroatoms_encountered > 1)
    return 0;

  if (ri.number_elements() - fused_ring_atoms_encountered < 6)   // at least 6 non-fused ring atoms
    return 0;

  if (heteroatom_outside_ring > 1)
    return 0;

  if (heteroatoms_encountered && heteroatom_outside_ring)
    return 0;

  return 1;
}

static int c7ring_count = 0;
static int apply_c7ring = 1;

static int
determine_c7_ring(Molecule & m,
                   Demerit & demerit,
                   const atomic_number_t * z,
                   const int * ncon)
{
  if (! apply_c7ring)
    return 0;

  int nr = m.nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    int ring_size = ri->number_elements();

    if (ring_size < ring_size_too_large)           // 7 by default
      continue;

    if (! is_rejectable_large_ring(m, *ri, z, ncon))
      continue;

    demerit.reject("LargeCarbocycle");
    c7ring_count++;
    return 1;
  }

  return 0;
}

/*
  Dealing with long carbon chains proved difficult with substructure
  search 

  Atom A is part of a carbon chain. Grow the chain
*/

static int
grow_chain(Molecule & m,
            atom_number_t zatom,
            const atomic_number_t * mz, const int * ncon,
            int * already_done)
{
  int rc = 0;
  while (1)
  {
    assert(0 == already_done[zatom]);

    already_done[zatom]++;
    rc++;

    if (1 == ncon[zatom])
      return rc;

    atom_number_t tmp = m.other(zatom, 0);
    if (already_done[tmp])
      tmp = m.other(zatom, 1);

    assert(0 == already_done[tmp]);

    if (6 != mz[tmp] || ncon[tmp] > 2 || m.is_ring_atom(tmp) || m.formal_charge(tmp))
      return rc;

    zatom = tmp;
  }
}

static int
determine_chain(Molecule & m,
                 atom_number_t zatom,
                 const atomic_number_t * mz, const int * ncon,
                 int * already_done)
{
//int matoms = m.natoms();

  already_done[zatom] = 1;
  int rc = 1;

  assert(2 == ncon[zatom]);

  const Atom * a = m.atomi(zatom);

  for (int i = 0; i < 2; i++)
  {
    atom_number_t a1 = a->other(zatom, i);
    assert(a1 >= 0 && a1 < m.natoms() && 0 == already_done[a1]);

    if (6 == mz[a1] && ncon[a1] <= 2 && m.is_non_ring_atom(a1) && 0 == m.formal_charge(a1))
    {
      rc += grow_chain(m, a1, mz, ncon, already_done);
    }
  }

  return rc;
}            

static int cx_chain_rejection_length = 7;

void
set_cx_chain_rejection_length(int s)
{
  cx_chain_rejection_length = s;
}

static int c7_count = 0;
static int apply_c7 = 1;
static int c6_count = 0;
static int apply_c6 = 1;
static int c5_count = 0;
static int apply_c5 = 1;
static int c4_count = 0;
static int apply_c4 = 1;

static int
long_carbon_chains(Molecule & m, Demerit & demerit,
                    const atomic_number_t * mz, const int * ncon)
{
  extending_resizable_array<int> long_chain;

  int matoms = m.natoms();

  int * already_done = new_int(matoms); std::unique_ptr<int[]> free_already_done(already_done);

  for (int i = 0; i < matoms; i++)
  {
    if (6 != mz[i])
      continue;

    if (2 != ncon[i])
      continue;

    if (already_done[i])
      continue;

    if (m.is_ring_atom(i))
      continue;

    const Atom * a = m.atomi(i);

    if (2 == a->nbonds() && 0 == a->formal_charge())
    {
      int path_length = determine_chain(m, i, mz, ncon, already_done);

      long_chain[path_length]++;
    }
  }

  int nl = long_chain.number_elements();

//#define ECHO_CHAIN_LENGTHS
#ifdef ECHO_CHAIN_LENGTHS
  for (int i = 0; i < nl; i++)
  {
    if (long_chain[i])
      cerr << "Found " << long_chain[i] << " chain sections of length " << i << endl;
  }
#endif

  if (nl < 4)
    return 0;

  if (apply_c4 && long_chain[cx_chain_rejection_length - 3])
  {
    c4_count++;
    demerit.extra(demerit_or_default_if_specified(10) * long_chain[4], "C4");
  }

  if (5 == nl)
    return 0;

  if (apply_c5 && long_chain[cx_chain_rejection_length - 2])
  {
    c5_count++;
    demerit.extra(demerit_or_default_if_specified(25) * long_chain[5], "C5");
  }

  if (6 == nl)
    return 0;

  if (apply_c6 && long_chain[cx_chain_rejection_length - 1])
  {
    c6_count++;
    demerit.extra(demerit_or_default_if_specified(50) * long_chain[6], "C6");
  }

  if (! apply_c7)
    return 0;

  int total_demerit = 0;

  for (int i = cx_chain_rejection_length; i < nl; i++)
  {
    if (0 == long_chain[i])
      continue;

    total_demerit += 100 + (i - cx_chain_rejection_length) * 50 * long_chain[i];
    c7_count++;
  }

  if (total_demerit)
    demerit.extra(total_demerit, "LongCChain");

  return 0;
}

/*
  The alkyl halide rule is complicated by the decision in Jul 97
  that CCl3 does not count as an alkyl halide
*/

static int alkyl_halides_count = 0;

#ifdef HALOGENS
static int too_many_bromines = 0;
static int too_many_chlorine_count = 0;
#endif

static int apply_alkyl_halides = 1;

static int
alkyl_halides(Molecule & m, Demerit & demerit,
               const atomic_number_t * mz, const int * ncon)
{
  if (! apply_alkyl_halides)
    return 0;

  resizable_array<atom_number_t> ccl3;    // carbons which are CCl3 centres

  int matoms = m.natoms();
  int nbromine = 0;
  int nchlorine = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (1 != ncon[i])
      continue;

    if (17 == mz[i])
      nchlorine++;
    else if (35 == mz[i])
      nbromine++;
    else if (53 == mz[i])
      ;
    else
      continue;

    const Atom * a = m.atomi(i);

    atom_number_t c = a->other(i, 0);
    if (6 != mz[c])
      continue;

    if (ccl3.contains(c))
      continue;

    if (ncon[c] < m.nbonds(c))
      continue;

//  We have a halogen attached to a saturated carbon.
//  Jul 97, do not consider CCl3 as an alkyl halide

    if (4 == ncon[c])
    {
      const Atom * ac = m.atomi(c);

      int nc = 0;    // count the chlorines attached to 'c'
      for (int j = 0; j < 4; j++)
      {
        atom_number_t k = ac->other(c, j);
        if (17 == mz[k])
          nc++;
      }

      if (nc >= 3)
      {
        ccl3.add(c);
        continue;      // not considered an alkyl halide
      }
    }

    demerit.reject("alkyl_halide");
    alkyl_halides_count++;
    return 1;
  }

  return 0;
}

/*
  Sulphur bonded to halogen, or another sulphur
*/

static int
is_phosphorus_acid (const Atom * p,
                    const atom_number_t zatom,
                    const atomic_number_t * mz, const int * ncon)
{
  int doubly_bonded_oxygen = 0;
  int singly_bonded_oxygen = 0;

  int pcon = p->ncon();

  for (int i = 0; i < pcon; i++)
  {
    const Bond * b = p->item(i);

    atom_number_t j = b->other(zatom);

    if (8 != mz[j])
      continue;

    if (b->is_double_bond() && 1 == ncon[j])
      doubly_bonded_oxygen++;
    else if (b->is_single_bond())   // note, no restriction on singly connected
      singly_bonded_oxygen++;
  }

  if (doubly_bonded_oxygen && singly_bonded_oxygen)
    return 1;
  else
    return 0;
}

/*
  Aug 2004. Differentiate between isolated cf3 groups
  and cf3 groups adjacent to cf2.
*/

static int molecules_with_cf2 = 0;
static int molecules_with_cf3 = 0;

static int
identify_cfx_section (const Molecule & m,
                      const atomic_number_t * z,
                      atom_number_t c,
                      int * counted,
                      int & cf2,
                      int & cf3)
{
  if (counted[c])
    return 0;

  counted[c] = 1;

  const Atom * a = m.atomi(c);

  int fluorine_connections = 0;

  Set_of_Atoms carbon_connections;

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(c, i);

    if (counted[j])
      continue;

    if (9 == z[j])
    {
      counted[j] = 1;
      fluorine_connections++;
    }
    else if (6 == z[j])
      carbon_connections.add(j);
  }

  if (2 == fluorine_connections)
    cf2++;
  else if (3 == fluorine_connections)
    cf3++;
  else
    return 0;

  int rc = 1;

  for (int i = 0; i < carbon_connections.number_elements(); i++)
  {
    atom_number_t c2 = carbon_connections[i];

    rc += identify_cfx_section(m, z, c2, counted, cf2, cf3);
  }

  return rc;
}

/*static int
fluorine (Molecule & m,
          Demerit & demerit,
          const atomic_number_t * z,
          const int * ncon)
{
  int cf2 = 0;
  int cf3 = 0;

  int matoms = m.natoms();

  int * counted = new_int (matoms); std::unique_ptr<int[]> free_counted (counted);

  for (int i = 0; i < matoms; i++)
  {
    if (counted[i])
      continue;

    if (9 != z[i])
      continue;

    if (0 == ncon[i])
      continue;

    atom_number_t c = m.other(i, 0);

    if (counted[c])
      continue;

    int tmpcf2 = 0;
    int tmpcf3 = 0;

    if (! identify_cfx_section (m, z, c, counted, tmpcf2, tmpcf3))
      continue;

//  If there was a cf2 in this grouping, then the whole thing is a cf2 grouping

    if (tmpcf2)
      cf2 += tmpcf2 + tmpcf3;
    else
      cf3++;
  }

  if (cf2)
  {
    demerit.extra (cf2 * 35, "CF2");
    molecules_with_cf2++;
  }

  if (cf3 < 2)
    return demerit.rejected();

  molecules_with_cf3++;

  if (2 == cf3)
    demerit.extra (10, "CF3");
  else if (3 == cf3)
    demerit.extra (60, "CF3");
  else
    demerit.reject("CF3");

  return demerit.rejected();
}*/


/*
  Note that this query is CPC, with no restrictions on other connections to the P

  Jun 97.
  Changed the rule for more than two phosphorus to not count phosphoric
  acids
*/

static int phosphorus_to_two_carbons_count = 0;
static int apply_phosphorus_to_two_carbons = 1;

static int phosphorus_to_nitrogen_count = 0;
static int apply_phosphorus_to_nitrogen = 1;

static int sulphur_phosphorus_bond_count = 0;
static int apply_sulphur_phosphorus_bond = 1;

static int more_than_two_phosporus_count = 0;
static int apply_more_than_two_phosporus = 1;

static int
phosphorus (Molecule & m, Demerit & demerit,
            const atomic_number_t * mz, const int * ncon)
{
  const int matoms = m.natoms();

  int phosphoric_acids = 0;
  int np = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (15 != mz[i])
      continue;

    const Atom * p = m.atomi(i);

    np++;

    if (ncon[i] > 2 && is_phosphorus_acid(p, i, mz, ncon))
      phosphoric_acids++;
    else if (np - phosphoric_acids > 1 && apply_more_than_two_phosporus)
    {
      demerit.reject("more_than_1_P");
      more_than_two_phosporus_count++;
      return 1;
    }

    int icon = ncon[i];

    int ncarbon = 0;
    for (int j = 0; j < icon; j++)
    {
      atomic_number_t z = mz[p->other(i, j)];

      if (6 == z)
      {
        ncarbon++;
        if (ncarbon > 1 && apply_phosphorus_to_two_carbons)
        {
          demerit.reject("CPC");
          phosphorus_to_two_carbons_count++;
          return 1;
        }
      }
      else if (7 == z && apply_phosphorus_to_nitrogen)
      {
        demerit.reject("PN");
        phosphorus_to_nitrogen_count++;
        return 1;
      }
      else if (16 == z && apply_sulphur_phosphorus_bond)        // prior to 21 Aug, this also specified a double bond
      {
        demerit.reject("P-S");
        sulphur_phosphorus_bond_count++;
        return 1;
      }
    }
  }

  return 0;
}

/*
  Molecules are rejected if they have two or more halogens attached to
  different aliphatic carbons.

  Hmmm, may 97, I think that the alkyl halide rejection supercedes this
*/

static int two_halogens_at_different_attach_points_count = 0;
static int apply_two_halogens_at_different_attach_points = 1;

static int
two_halogens_at_different_attach_points (Molecule & m, Demerit & demerit,
               const atomic_number_t * mz, const int * ncon)
{
  if (! apply_two_halogens_at_different_attach_points)
    return 0;

  int first_carbon_centre = INVALID_ATOM_NUMBER;

  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = mz[i];
    if ((17 == z || 35 == z || 53 == z) && 1 == ncon[i])
    {
      atom_number_t c = m.other(i, 0);
      atomic_number_t cz = mz[c];
      if (6 == cz && ncon[c] == m.nbonds(c))    // SP3 carbon only
      {
        if (INVALID_ATOM_NUMBER == first_carbon_centre)
          first_carbon_centre = c;
        else if (c != first_carbon_centre)
        {
          demerit.reject("halogens_to_multiple_aliphatic_C");
          two_halogens_at_different_attach_points_count++;
          return 1;
        }
      }
    }
  }

  return 0;
}

/*
  More than two [SD2]
*/

static int more_than_two_sulphur_count = 0;
static int apply_more_than_two_sulphur = 1;

static int
more_than_two_sulphur (Molecule & m, Demerit & demerit,
                       const atomic_number_t * mz, const int * ncon)
{
  if (! apply_more_than_two_sulphur)
    return 0;

  int nsulphur = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (16 == mz[i] && 2 == ncon[i] && 0 == m.nrings(i))
      nsulphur++;
  }

  if (nsulphur < 3)
    return 0;

  more_than_two_sulphur_count++;
  demerit.extra(demerit_or_default_if_specified(50) * (nsulphur - 1), "more_than_2_[SD2]");

  return 1;
}

/*
  Rule here is [SD2]-N
  And [SD2]=N
  And [SD3]=N

  In Apr 96 this was expanded to include
    [SD2R]=N
  In May 96 this expanded to include all S=N outside a ring
*/

static int sulphur_nitrogen_count = 0;
static int apply_sulphur_nitrogen = 1;

static int
sulphur_nitrogen (Molecule & m, Demerit & demerit,
               const atomic_number_t * mz, const int * ncon)
{
  if (! apply_sulphur_nitrogen)
    return 0;

  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (16 != mz[i] || ncon[i] < 2 || ncon[i] > 4 || m.formal_charge(i))
      continue;

    if (m.is_aromatic(i))
      continue;

//  We now have a 2 or 3 connected S

    const Atom * as = m.atomi(i);
    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = as->item(j);

      atom_number_t k = b->other(i);

      if (7 != mz[k])
        continue;

      if (m.is_ring_atom(k))
        continue;

      if (2 == ncon[i] && b->is_single_bond())
        demerit.reject("[SD2]-N");
      else if (2 == ncon[i] && b->is_double_bond())
        demerit.reject("[SD2]=N");
//    else if (3 == ncon[i] && b->is_double_bond())    // removed 25 sept 2010
//      demerit.reject("[SD3]=N");
      else if (4 == ncon[i] && 6 == as->nbonds())    // sulphone type OK
        continue;
      else
        demerit.reject("sulfinamide");

      sulphur_nitrogen_count++;
      return 1;
    }
  }

  return 0;
}


static int nitrogen_nitrogen_double_bond_count = 0;
static int apply_nitrogen_nitrogen_double_bond = 1;

static int nitrogen_nitrogen_triple_bond_count = 0;
static int apply_nitrogen_nitrogen_triple_bond = 1;

static int two_nitrogens_with_double_bonds_count = 0;
static int apply_two_nitrogens_with_double_bonds = 1;

static int
nitrogen_nitrogen (Molecule & m, Demerit & demerit,
               const atomic_number_t * mz, const int * ncon)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (7 != mz[i])
      continue;

    if (m.nrings(i))
      continue;

    const Atom * a = m.atomi(i);

    int icon = ncon[i];
    for (int j = 0; j < icon; j++)
    {
      const Bond * b = a->item(j);
      atom_number_t k = b->other(i);

      if (7 == mz[k] && 0 == m.nrings(k) && apply_nitrogen_nitrogen_triple_bond)     // two N's in a chain.
      {
        if (b->is_triple_bond())
        {
          demerit.reject("N#N");
          nitrogen_nitrogen_triple_bond_count++;
          return 1;
        }
        if (b->is_double_bond())
        {
          demerit.reject("N=N");
          nitrogen_nitrogen_double_bond_count++;
          return 1;
        }
        if (b->is_single_bond())
        {
          if (a->nbonds() > ncon[i] && m.nbonds(k) > ncon[k])
          {
            demerit.reject("=N-N=");
            two_nitrogens_with_double_bonds_count++;
            return 1;
          }
        }
      }
    }
  }

  return 0;
}

static int nitrogen_single_bond_nitrogen_count = 0;
static int apply_nitrogen_single_bond_nitrogen = 1;

/*
  We have an N-N bond, need to see if the second nitrogen is doubly bonded to a carbon
*/

static int
is_hydrazid (Molecule & m,
             const atom_number_t n,
             const atomic_number_t * mz,
             const int * ncon)
{
  const Atom * a = m.atomi(n);

  const int acon = ncon[n];

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())    // redudant
      continue;

    const atom_number_t j = b->other(n);

    const atomic_number_t zj = mz[j];

    if (7 == zj)   // looping back to the starting N atom
      continue;

    if (6 == zj || 16 == zj)
      ;
    else
      continue;

    if (m.is_aromatic(j))
      continue;

    const int jcon = ncon[j];

    if (6 == zj && 3 == jcon)
      ;
    else if (32 == zj && 4 == jcon)
      ;
    else
      continue;

    const Atom * aj = m.atomi(j);

    for (int x = 0; x < jcon; ++x)
    {
      const Bond * b = aj->item(x);

      if (! b->is_double_bond())
        continue;

      const atom_number_t o = b->other(j);

      if (8 == mz[o])
        return 1;
      else
        return 0;
    }
  }

  return 0;
}

static int
is_hydrazone (const Molecule & m,
              const atom_number_t n,
              const atomic_number_t * mz,
              const int * ncon)
{
  const Atom * a = m.atomi(n);

  if (0 != a->formal_charge())
    return 0;

  const int acon = ncon[n];

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    const atom_number_t c = b->other(n);

    if (6 == mz[c])
      return 1;
  }

  return 0;
}

static int
is_hydrazone_or_hydrazid (Molecule & m,
                          const atom_number_t n,
                          const atomic_number_t * mz,
                          const int * ncon)
{
  const Atom * a = m.atomi(n);

  const int acon = ncon[n];

  if (2 != acon)
    return 0;

  const int nb = a->nbonds();

  if (3 == nb)
    return is_hydrazone(m, n, mz, ncon);

  if (2 == nb)
    return is_hydrazid(m, n, mz, ncon);

  return 0;
}

// Make sure we do not hit Hydrazones, or Hydrazids

static int
nitrogen_single_bond_nitrogen (Molecule & m,
                        Demerit & demerit,
                        const atomic_number_t * mz,
                        const int * ncon)
{
  if (! apply_nitrogen_single_bond_nitrogen)
    return 0;

  const int matoms = m.natoms();
  int NNcount = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (7 != mz[i])
      continue;

    const Atom * a = m.atomi(i);

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (k < i)   // make sure we only count these once
        continue;

      if (7 != mz[k])
        continue;

      if (1 == ncon[i] && is_hydrazone_or_hydrazid(m, k, mz, ncon))
        continue;
      else if (1 == ncon[k] && is_hydrazone_or_hydrazid(m, i, mz, ncon))
        continue;

      if (! m.in_same_rings(i, k))
        NNcount++;
    }
  }

  if (NNcount)
  {
    demerit.extra(NNcount * demerit_or_default_if_specified(75), "N-N");   // jun97 change from 50 to 75 each
    nitrogen_single_bond_nitrogen_count++;
  }

  return 0;
}

static int apply_too_many_rings_demerit = 1;

static  int
too_many_rings (Molecule & m,
                Demerit & demerit)
{
  if (! apply_too_many_rings_demerit)
    return 0;

  int nr = m.nrings();

  if (nr > too_many_rings_cutoff)
    demerit.extra((nr - too_many_rings_cutoff) * demerit_or_default_if_specified(20), "too_many_rings");

  return demerit.rejected();
}

/*static int apply_aromatic_rings_demerit = 1;

static int
too_many_aromatic_rings (Molecule & m,
                         Demerit & demerit)
{
  if (! apply_aromatic_rings_demerit)
    return 0;

  int nr = m.nrings();

  if (nr < 5)
    return 0;

  m.compute_aromaticity_if_needed();

  int aromatic_rings = 0;

  for (int i = 0; i < m.nrings(); i++)
  {
    const Ring * ri = m.ringi (i);

    if (ri->is_aromatic())
      aromatic_rings++;
  }

  if (aromatic_rings < 5)
    return 0;

  demerit.extra (40 * (aromatic_rings - 4), "Aromatic Rings");

  return demerit.rejected();
}*/

static int
hard_coded_queries (Molecule & m,
                    Demerit & demerit,
                    const atomic_number_t * z,
                    const int * ncon)
{
  if (alkyl_halides(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

//if (fluorine(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
//  return 1;

  if (phosphorus(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (two_halogens_at_different_attach_points(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (nitrogen_nitrogen(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (nitrogen_single_bond_nitrogen(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (too_many_rings(m, demerit) && 0 == keep_going_after_rejection)
    return 1;

  if (more_than_two_sulphur(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (sulphur_nitrogen(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (too_many_charges(m, demerit) && 0 == keep_going_after_rejection)
    return 1;

  if (long_carbon_chains(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (determine_c7_ring(m, demerit, z, ncon) && 0 == keep_going_after_rejection)
    return 1;

  if (apply_satcg)
  {
    int * ring_membership = new int[m.natoms()]; std::unique_ptr<int[]> free_ring_membership(ring_membership);
    m.ring_membership(ring_membership);

    if (large_saturated_carbon_sections_including_rings(m, demerit, z, ncon, ring_membership) && 0 == keep_going_after_rejection)
      return 1;
  }

//if (complex_fused_rings(m, demerit))
//  return 1;

  return 0;
}

int
hard_coded_queries(Molecule & m, Demerit & demerit)
{
  int matoms = m.natoms();

  int * ncon = new int[matoms]; std::unique_ptr<int[]> free_ncon(ncon);
  (void) m.ncon(ncon);

  atomic_number_t * z = new atomic_number_t[matoms]; std::unique_ptr<atomic_number_t[]> free_z(z);
  (void) m.atomic_numbers(z);

  return hard_coded_queries(m, demerit, z, ncon);
}

int
hard_coded_queries_statistics(std::ostream & os)
{
  os << "alkyl_halides_count = " << alkyl_halides_count << endl;

#ifdef HALOGENS
  os << "too_many_bromines = " << too_many_bromines << endl;
  os << "too_many_chlorine_count = " << too_many_chlorine_count << endl;
#endif

  os << molecules_with_cf2 << " molecules containing CF2, " << molecules_with_cf3 << " with CF3 groups\n";

  os << "phosphorus_to_two_carbons_count = " << phosphorus_to_two_carbons_count << endl;

  os << "phosphorus_to_nitrogen_count = " << phosphorus_to_nitrogen_count << endl;

  os << "more_than_two_phosporus_count = " << more_than_two_phosporus_count << endl;

  os << "two_halogens_at_different_attach_points_count = " << two_halogens_at_different_attach_points_count << endl;

  os << "more_than_two_sulphur_count = " << more_than_two_sulphur_count << endl;

  os << "sulphur_nitrogen_count = " << sulphur_nitrogen_count << endl;

  os << "two_nitrogens_with_double_bonds_count = " << two_nitrogens_with_double_bonds_count << endl;

  os << "nitrogen_nitrogen_double_bond_count = " << nitrogen_nitrogen_double_bond_count << endl;

  os << "nitrogen_nitrogen_triple_bond_count = " << nitrogen_nitrogen_triple_bond_count << endl;

  os << "nitrogen_single_bond_nitrogen_count = " << nitrogen_single_bond_nitrogen_count << endl;

  os << "satcg_count = " << satcg_count << endl;

  os << "c7_count = " << c7_count << endl;

  os << "c6_count = " << c6_count << endl;

  os << "c5_count = " << c5_count << endl;

  os << "c4_count = " << c4_count << endl;

  os << "c7ring_count = " << c7ring_count << endl;

  os << "complexfusedrings " << complex_fused_rings_count << endl;

  return 1;
}

/*
  Dec 2001.
  We want the ability to turn on and turn off individual tests

  The first token on the line must be a directive we recognise.

  If that is followed by a 0, that test will not be done
*/

int
suppress_query(const const_IWSubstring & query_name)
{
  if ("c4" == query_name)
    apply_c4 = 0;
  else if ("c5" == query_name)
    apply_c5 = 0;
  else if ("c6" == query_name)
    apply_c6 = 0;
  else if ("c7" == query_name)
    apply_c7 = 0;
  else if ("satcg" == query_name)
    apply_satcg = 0;
  else if ("c7ring" == query_name)
    apply_c7ring = 0;
  else if (query_name.starts_with ("alkyl_halide"))
    apply_alkyl_halides = 0;
  else if ("phosphorus_to_two_carbon" == query_name)
    apply_phosphorus_to_two_carbons = 0;
  else if ("phosphorus_nitrogen" == query_name)
    apply_phosphorus_to_nitrogen = 0;
  else if ("phosphorus_sulphur" == query_name)
    apply_sulphur_phosphorus_bond = 0;
  else if ("too_many_phosphorus" == query_name)
    apply_more_than_two_phosporus = 0;
  else if ("halogens_at_different_attach_points" == query_name)
    apply_two_halogens_at_different_attach_points = 0;
  else if ("more_than_two_sulphur" == query_name)
    apply_more_than_two_sulphur = 0;
  else if ("sulphur_nitrogen" == query_name)
    apply_sulphur_nitrogen = 0;
#ifdef SCHIFF_BASE
  else if ("nd3_plus" == query_name)
    apply_nd3_plus = 0;
#endif
  else if ("nitrogen_nitrogen_double_bond" == query_name)
     apply_nitrogen_nitrogen_double_bond = 0;
  else if ("nitrogen_nitrogen_triple_bond" == query_name)
     apply_nitrogen_nitrogen_triple_bond = 0;
  else if ("two_nitrogens_with_double_bonds" == query_name)
    apply_two_nitrogens_with_double_bonds = 0;
  else if ("nitrogen_single_bond_nitrogen" == query_name)
    apply_nitrogen_single_bond_nitrogen = 0;
  else if ("positive_charge" == query_name)
    apply_positive_charge_demerit = 0;
  else if ("negative_charge" == query_name)
    apply_negative_charge_demerit = 0;
  else if ("rings" == query_name)
    apply_too_many_rings_demerit = 0;
  else
    return 0;

  return 1;
}

int
initialise_queries_to_do (const const_IWSubstring & buffer)
{
  const_IWSubstring t1;

  buffer.word (0, t1);

  if (buffer.nwords() < 2)     // all queries are on by default, so just one token on the line must mean do the query, therefore nothing to do here
    return 1;

  const_IWSubstring t2;
  buffer.word(1, t2);

  if ('1' == t2 || "yes" == t2)
    return 1;

  if ('0' == t2)
    ;
  else if ("no" == t2)
    ;
  else
  {
    cerr << "Invalid qualifier, '" << t2 << "' to '" << t1 << " directive\n";
    return 0;
  }

  if (suppress_query(t1))
  {
    cerr << "Demerit for '" << t1 << "' suppressed\n";
    return 1;
  }

  cerr << "Unrecognised directive '" << buffer << "'\n";
  return 0;
}

int
initialise_queries_to_do(iwstring_data_source & input)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    buffer.strip_trailing_blanks();

    if (0 == buffer.length())
      continue;

    if (! initialise_queries_to_do(buffer))
    {
      cerr << "initialise_queries_to_do: invalid or unrecognised directive on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

int
initialise_hard_coded_queries_to_do(IWString & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return initialise_queries_to_do(input);
}

void
set_only_apply_rejection_rules()
{
  suppress_query("c4");
  suppress_query("c5");
  suppress_query("c6");
  suppress_query("positive_charge");
  suppress_query("negative_charge");
  suppress_query("nitrogen_single_bond_nitrogen");
  suppress_query("rings");

  return;
}
}    // end of substructure_demerits namespace
