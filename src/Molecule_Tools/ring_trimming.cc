/*
  Exhaustively trim rings from ring systems
*/

#include <stdlib.h>
#include <memory>
#include <algorithm>
using namespace std;

#include "cmdline.h"
#include "iw_auto_array.h"
#include "iw_stl_hash_set.h"
#include "misc.h"

#include "istream_and_type.h"
#include "path.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int write_parent_molecule = 0;
static int write_ring_systems = 0;
static int write_molecular_scaffold = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static extending_resizable_array<int> fragment_count;

static int isotope_for_scaffold_join_pts = 0;

static int isotope_for_ring_fusion = 0;

static int max_rings_in_system = 4;

static int max_rings_to_remove = 5;

static int remove_chiral_centres = 0;

static int suppress_duplicate_structures = 0;

static int do_preference_ranking = 0;

static int allow_destruction_of_aromaticity = 0;

static int add_first_removed_atom_number = 0;

static int remove_strongly_fused_as_a_unit = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Exhaustively trim rings from ring systems, preserving aromaticity\n";
  cerr << "  -w parent     write parent molecule\n";
  cerr << "  -w rings      write isolated ring systems\n";
  cerr << "  -w scaffold   write the molecular scaffold\n";
  cerr << "  -j <iso>      isotope for where scaffold joined the rest of the molecule\n";
  cerr << "  -J <iso>      isotope for where ring joins are broken\n";
  cerr << "  -R <nrings>   maximum number of rings in a ring system (def 4)\n";
  cerr << "  -m <n>        maximum number of rings to remove from a ring system\n";
  cerr << "  -c            remove all chiral centres\n";
  cerr << "  -u            only unique structures from each input molecule\n";
  cerr << "  -p            form rings in preference order - also annotate output\n";
  cerr << "  -a            allow destruction of aromaticity\n";
  cerr << "  -f            remove strongly fused ring systems as a unit\n";
  cerr << "  -y            append atom number of first atom removed (Jibo)\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (remove_chiral_centres)
    m.remove_all_chiral_centres();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

#define ATOM_STATUS_RING_ATOM 1
#define ATOM_STATUS_CHAIN 2
#define ATOM_STATUS_DOUBLY_BONDED_TO_RING 4
#define ATOM_STATUS_DOUBLY_BONDED_TO_CHAIN 8
#define ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT 256

/*
  Processing exocyclic double bonds, and doubly bonded atoms attached to the scaffold create lots of
  problems for us. We set up a class that can hold the status of each atom. We then set the
  user_specified_atom_void_ptr of each atom to point to an instance of this class.
*/

class Doubly_Bonded_to_Singly_Connected_Status
{
  private:
    int _status;

  public:
    Doubly_Bonded_to_Singly_Connected_Status();
};

/*
  We need to keep track of the parent child relationships.
  Every fragment gets assigned a unique ID
*/

class Parents
{
  private:
    int _next_to_give;

  public:
    Parents ();

    int unique_id();
};

Parents::Parents ()
{
  _next_to_give = 0;

  return;
}

int
Parents::unique_id ()
{
  const auto rc = _next_to_give;
  _next_to_give++;
  return rc;
}

static int
do_output (Molecule & m,
           const int parent_id,
           const int my_id,
           const int rings_removed,
           const char preference,
           IW_STL_Hash_Set & unique_smiles,
           int first_atom_removed,
           const IWString & description,
           IWString_and_File_Descriptor & output)
{
  if (suppress_duplicate_structures && unique_smiles.contains(m.unique_smiles()))
    return 1;

  unique_smiles.insert(m.unique_smiles());

  output << m.smiles() << ' ' << m.name() << ' ' << rings_removed << " nr=" << m.nrings() << ' ' << parent_id << ' ' << my_id << ' ' << preference << ' ' << description;

  if (add_first_removed_atom_number)
    output << ' ' << first_atom_removed;

  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
heteroatom_count (const Molecule & m,
                  const Set_of_Atoms & s)
{
  const auto n = s.size();

  int rc = 0;

  for (auto i = 0; i < n; ++i)
  {
    const auto j = s[i];

    if (6 != m.atomic_number(j))
      rc++;
  }

  return rc;
}

static int
heteroatom_count (const Molecule & m,
                  const int * being_processed)
{
  const auto matoms = m.natoms();

  int rc = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (0 == being_processed[i])
      continue;

    if (6 != m.atomic_number(i))
      rc++;
  }

  return rc;
}

static int
size_of_fused_system (Molecule & m,
                      const int fsid)
{
  const auto nr = m.nrings();

  int rc = 0;

  for (auto i = 0; i < nr; ++i)
  {
    if (m.ringi(i)->fused_system_identifier() == fsid)
      rc++;
  }

  return rc;
}

static int
ring_preference (Molecule & m,
                 const int rn)
{
  const Ring & r = *(m.ringi(rn));

  int rc = 1000000 * r.number_elements();

  rc += 1000 * size_of_fused_system(m, r.fused_system_identifier());

  if (r.is_aromatic())
    rc += rc + 232;

  rc += 10 * heteroatom_count(m, r);

  rc += r.largest_number_of_bonds_shared_with_another_ring();

  return rc;
}

static void
ring_description (Molecule & m,
                  const int rn,
                  IWString & s)
{
  const Ring & r = *(m.ringi(rn));

  s << 'F' << r.fused_ring_neighbours();

  if (r.is_aromatic())
    s << 'a';
  else
    s << 'A';

  s << r.number_elements();

  s << ',' << heteroatom_count(m, r);

  s << ',' << r.largest_number_of_bonds_shared_with_another_ring();

  return;
}

static void
rank_rings (Molecule & m,
            std::pair<int, int> * rings)
{
  const auto nr = m.nrings();

  m.compute_aromaticity_if_needed();

  for (auto i = 0; i < nr; ++i)
  {
    rings[i].first = i;

    rings[i].second = ring_preference(m, i);
  }

  std::sort(rings, rings + nr, [] (const std::pair<int, int> & lhs, const std::pair<int, int> & rhs) { return lhs.second > rhs.second;});

#ifdef ECHO_SORTED_VALUES
  for (auto i = 0; i < nr; ++i)
  {
    cerr << " after sort i = " << i << " ring " << rings[i].first << " preference " << rings[i].second << endl;
  }
#endif
}

static int
first_non_zero (const int * v, const int n)
{
  for (auto i = 0; i < n; ++i)
  {
    if (0 != v[i])
      return i;
  }

  return -1;
}

static void
expand_two_connected (const Molecule & m,
                      const atom_number_t zatom,
                      int * ring_system,
                      int * to_remove)
{
  if (to_remove[zatom])
    return;

  const Atom * a = m.atomi(zatom);

  if (2 != a->ncon())
    return;

  to_remove[zatom] = 1;
  ring_system[zatom] = -1;

  expand_two_connected(m, a->other(zatom, 0), ring_system, to_remove);
  expand_two_connected(m, a->other(zatom, 1), ring_system, to_remove);

  return;
}

static int
next_set_of_atoms_from_strongly_fused_ring (Molecule & m,
                                            int & ndx,
                                            int * ring_system,
                                            const int lbl,
                                            int * to_remove,
                                            atom_number_t & first_ring_atom_removed)
{
  const auto matoms = m.natoms();

  std::fill_n(to_remove, matoms, 0);

  for (; ndx < matoms; ++ndx)
  {
    if (lbl != ring_system[ndx])
      continue;

    const auto a = m.atomi(ndx);

    if (2 != a->ncon())
      continue;

    to_remove[ndx] = 1;
    ring_system[ndx] = -1;
    expand_two_connected(m, a->other(ndx, 0), ring_system, to_remove);
    expand_two_connected(m, a->other(ndx, 1), ring_system, to_remove);

    first_ring_atom_removed = first_non_zero(to_remove, matoms);
    return 1;
  }

  return 0;
}

/*
  We are thining about removing an aromatic ring from the molecule.
  If that would leave an adjacent ring with single bonds at the
  join point, we need to convert that bond to a double bond
*/

static int
double_bond_to_retained_atoms (const Molecule & m,
                               const int * atoms_being_removed,
                               const atom_number_t a1,
                               const atom_number_t a2)
{
  const Atom * a = m.atomi(a1);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    const auto j = b->other(a1);

    if (! atoms_being_removed[j])   // we have found a double bond to an atom being retained
      return 1;

    if (j == a2)     // the bond across the two rings is already a double bond
      return 0;
  }

  return 0;
}

static int
gather_doubly_bonded_offshoots (const Molecule & m,
                                const atom_number_t zatom,
                                const int * ignore,
                                Set_of_Atoms & doubly_bonded_offshoots)
{
  const auto a = m.atomi(zatom);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    const auto j = b->other(zatom);

    if (ignore[j])
      continue;

    if (1 != m.ncon(j))
      continue;

    doubly_bonded_offshoots.add(j);
  }

  return doubly_bonded_offshoots.size();
}




static void
identify_exposed_chains (const Molecule & m,
                         const atom_number_t zatom,
                         int * to_remove)
{
  const Atom * a = m.atomi(zatom);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    const auto j = b->other(zatom);

    if (to_remove[j])    // already processed
      continue;

    const Atom * aj = m.atomi(j);

    const int v = * reinterpret_cast<const int *>(aj->user_specified_void_ptr());

    if (ATOM_STATUS_RING_ATOM & v)    // we got to a ring system, might be done
      continue;

    const auto jcon = aj->ncon();

    if (2 == jcon)               // 2 connected, can continue down here
    {
      to_remove[j] = 1;
      identify_exposed_chains(m, j, to_remove);
      continue;
    }

//  If J has just one connection, it will be removed

    if (1 == jcon)
    {
      to_remove[j] = 1;
      continue;
    }

    continue;    // I don't think the code below is now needed

//  If connections are to doubly bonded offshoots, we may be able to continue

    if (jcon == aj->nbonds())    // fully saturated and > 2 connections, cannot continue down here
      continue;

    if (0 == (ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT & v))    // atom is not part of a singly connected double bond
      continue;

    to_remove[j] = 1;

    Set_of_Atoms doubly_bonded_offshoots;
    gather_doubly_bonded_offshoots(m, j, to_remove, doubly_bonded_offshoots);

    doubly_bonded_offshoots.set_vector(to_remove, 1);

    identify_exposed_chains(m, j, to_remove);
  }
}

/*
  We have removed a ring. If there are any chains that have been exposed, we need to trim them
  down so we get to a scaffold structure

  This is complicated because we need to look for two cases

  (1) atoms that started as chain atoms, and now have 1 connection 
  (2) atoms that started as ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT and now have two connections (C=O adjacent to ring for example)
  (3) atoms that started as ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT and now have three connections (SO2 adjacent to ring for example)
*/

static void
identify_exposed_chains (const Molecule & m,
                         int * to_remove)
{
  const auto matoms = m.natoms();

  std::fill(to_remove, to_remove + m.natoms(), 0);

  for (auto i = 0; i < matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    const auto acon = a->ncon();

//  look for atoms that started as chain atoms, but now are singly connected

    const int v = * reinterpret_cast<const int *>(a->user_specified_void_ptr());

    if (ATOM_STATUS_CHAIN == v && 1 == acon)      // case (1)
    {
      to_remove[i] = 1;
      identify_exposed_chains(m, i, to_remove);
    }
    else if (2 == acon && (ATOM_STATUS_CHAIN | ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT) == v && 3 == a->nbonds())
    {
      to_remove[i] = 1;
      identify_exposed_chains(m, i, to_remove);
    }
    else if (3 == acon && (ATOM_STATUS_CHAIN | ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT) == v && 5 == a->nbonds())    // SO2, maybe phosphorus
    {
      to_remove[i] = 1;
      identify_exposed_chains(m, i, to_remove);
    }
  }

  return;
}

/*
  Atom zatom is part of a ring being removed. There is an exocyclic double bond,
  so it must be removed too
*/

static void
identify_doubly_bonded_non_ring_atom (Molecule & m, 
                                      const atom_number_t zatom,
                                      int * atoms_to_remove)
{
  const auto a = m.atomi(zatom);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    const auto j = b->other(zatom);

    if (1 == m.ncon(j))
      atoms_to_remove[j] = 1;
  }

  return;
}

/*
  Remove a specific ring. Note that this will NOT properly handle the case
  of a ring that has every atom also part of another ring. Those are too
  rare to worry about
*/

static int
remove_ring (Molecule & m, 
             const int * ring_membership,
             const Ring & ring_being_removed,
             int * to_be_removed,
             const int removed_ring_was_aromatic,
             atom_number_t & first_ring_atom_removed)
{
  const auto matoms = m.natoms();

  std::fill(to_be_removed, to_be_removed + matoms, 0);    // we will remove these atoms

  int atoms_removed = 0;

  for (Ring_Atom_Iterator i(ring_being_removed); i != ring_being_removed.end(); i++)
  {
    atom_number_t j = i.current();

    if (ring_membership[j] > 1)     // shared with other rings, cannot be removed
      continue;

    atoms_removed++;

    to_be_removed[j] = 1;

    const auto a = m.atomi(j);

    if (2 == a->ncon())    // nothing branching off here
      continue;

    const int v = *reinterpret_cast<const int *>(a->user_specified_void_ptr());

    if (ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT & v)
      identify_doubly_bonded_non_ring_atom(m, j, to_be_removed);
  }

  if (0 == atoms_removed)
  {
    if (verbose)
      cerr << "Ring removal: no unshared atoms, ring removal has failed '" << m.name() << "' " << m.smiles() << endl;
    return 0;
  }

  if (add_first_removed_atom_number)
    first_ring_atom_removed = first_non_zero(to_be_removed, matoms);

// Are there any bonds that should be converted to double bonds. 

  if (removed_ring_was_aromatic)
  {
    Set_of_Atoms db1, db2;

    for (Ring_Bond_Iterator i(ring_being_removed); i != ring_being_removed.end(); i++)
    {
      const auto a1 = i.a1();

      if (1 == to_be_removed[a1])    // atom is being removed
        continue;

      const auto a2 = i.a2();
      if (1 == to_be_removed[a2])    // atom is being removed
        continue;

      const Bond * b = m.bond_between_atoms(a1, a2);

      if (b->is_double_bond())
        continue;

//    But if there are already double bonds to retained atoms, we cannot place a double bond

     if (double_bond_to_retained_atoms(m, to_be_removed, a1, a2))
       continue;

     if (double_bond_to_retained_atoms(m, to_be_removed, a2, a1))
       continue;

      db1.add(a1);
      db2.add(a2);
    }

    for (auto i = 0; i < db1.number_elements(); ++i)
    {
      m.set_bond_type_between_atoms(db1[i], db2[i], DOUBLE_BOND);
    }
  }

  if (isotope_for_ring_fusion)
  {
    for (Ring_Atom_Iterator i(ring_being_removed); i != ring_being_removed.end(); i++)
    {
      const auto j = i.current();

      if (0 == to_be_removed[j])    // atom will be kept, but it is part of our ring
        m.set_isotope(j, isotope_for_ring_fusion);
    }
  }

#ifdef DEBUG_RINT_TRIMMING
  for (auto i = 0; i < m.natoms(); ++i)
  {
    if (to_be_removed[i])
      cerr << "Will remove atom " << i << " " << m.smarts_equivalent_for_atom(i) << endl;
  }
#endif

  m.remove_atoms(to_be_removed);

#ifdef DEBUG_RINT_TRIMMING
  cerr << "Removed the atoms " << m.smiles() << ", valence " << m.valence_ok() << endl;
#endif

  if (! m.valence_ok())   // the double bond we added might have created problems
    return 0;

  if (1 != m.number_fragments())    // we have removed something from the inside
    return 0;

  identify_exposed_chains(m, to_be_removed);

  m.remove_atoms(to_be_removed);

  return 1;
}

static int
identify_chain_atoms (Molecule & m,
                      const atom_number_t prev,
                      const atom_number_t zatom,
                      int * ring_membership)
{
  const auto a = m.atomi(zatom);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    const auto j = a->other(zatom, i);

    if (-1 == ring_membership[j])
      continue;
  }

  return 1;
}


static int
identify_strongly_fused_ring_group (Molecule & m, 
                                  const Ring & r,
                                  int * to_remove,
                                  int * ring_already_done)
{
  const auto ring_number = r.ring_number();

  ring_already_done[ring_number] = 1;

  r.set_vector(to_remove, 1);

  int n = r.fused_ring_neighbours();

  for (auto i = 0; i < n; ++i)
  {
    const auto nbr = r.fused_neighbour(i);

    const auto n = nbr->ring_number();

    if (ring_already_done[n])
      continue;

    if (r.compute_bonds_shared_with(*nbr) < 2)        // not strongly fused
      continue;

    identify_strongly_fused_ring_group(m, *nbr, to_remove, ring_already_done);
  }

  return 1;
}

//#define DEBUG_RINT_TRIMMING

static int
remove_strongly_fused_ring_group (Molecule & m, 
                                  const Ring & r,
                                  int * to_remove,
                                  int * ring_already_done,
                                  atom_number_t & first_ring_atom_removed)
{
  const auto matoms = m.natoms();

  std::fill(to_remove, to_remove + matoms, 0);

  auto rc = identify_strongly_fused_ring_group(m, r, to_remove, ring_already_done);

  if (add_first_removed_atom_number)
    first_ring_atom_removed = first_non_zero(to_remove, matoms);

  m.remove_atoms(to_remove);

#ifdef DEBUG_RINT_TRIMMING
    iw_write_array(to_remove, matoms, "Strongly fused", cerr);
    m.remove_atoms(to_remove);
    cerr << "After removal " << m.smiles() << endl;
#endif

  return rc;
}

static char
symbol_representing_preference (const std::pair<int, int> * rscore,
                                const int nr,
                                int ndx)
{
  if (0 == ndx)
  {
    if (1 == nr)
      return 'F';

    if (rscore[0].second == rscore[1].second)    // there is a tie
      return 'f';
    else
      return 'F';
  }

  if (rscore[ndx].second == rscore[0].second)   // tied for first
    return 'f';

  return '*';
}

#define MAGIC_NUMBER_INDICATING_STRONGLY_FUSED 23412341

/*
  We are processing a ring that is strongly fused. We now need to mark all the other rings
  that are in that same system as done
*/

static int
mark_all_other_rings_in_system (Molecule & m,
                                int istart,
                                int * ring_already_done)
{
  const auto nr = m.nrings();

  const auto fsid = m.ringi(istart)->fused_system_identifier();

  int rc = 1;

  for (auto i = 0; i < nr; ++i)
  {
    if (i == istart)
      continue;

    if (1 == ring_already_done[i])    // not strongly fused, was processed normally
      continue;

    if (fsid == m.ringi(i)->fused_system_identifier())
    {
      ring_already_done[i] = 1;
      rc++;
    }
  }

  return rc;
}

static int
ok_to_keep (Molecule & m,
            const int aromatic_ring_removed,   // aromaticity of the ring just removed
            const int parent_aromatic_ring_count)     // total aromatic ring count in parent
{
  if (0 == m.natoms() || 1 != m.number_fragments())    // definitely no good
    return 0;

  const auto arc = m.aromatic_ring_count();

  if (arc == parent_aromatic_ring_count)    // aromaticity not changed, we must have removed an aliphatic ring
    ;
  else if (aromatic_ring_removed && arc + 1 == parent_aromatic_ring_count)    // we removed an aromatic ring and now there are one fewer aromatic rings, good
    ;
  else if (allow_destruction_of_aromaticity) // we seem to have destroyed aromaticity
    ;
  else
    return 0;

  return 1;
}

static int
ring_membership_preserved (Molecule & m)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    const int v = *reinterpret_cast<const int *>(a->user_specified_void_ptr());

    if (0 == (ATOM_STATUS_RING_ATOM & v))    // atom was not a ring atom in the parent
      continue;

    if (! m.is_ring_atom(i))
      return 0;
  }

  return 1;
}

/*
  This is the main processing loop. It will be called recursively.

  First we try removing each ring in the molecule one at a time.
  Reject if there are now two fragments.
  Remove any chains thereby created.
  Recurse
*/

static int
ring_trimming_from_scaffold (Molecule & m,
                             int * to_remove,
                             Parents & parents,
                             const int parent_id,
                             IW_STL_Hash_Set & unique_smiles,
                             int rings_removed_so_far,
                             IWString_and_File_Descriptor & output)
{
  const auto matoms = m.natoms();

  m.compute_aromaticity_if_needed ();

  const auto nr = m.nrings();

  int * ring_already_done = new_int(nr); unique_ptr<int[]> free_ring_already_done(ring_already_done);

  const auto aromatic_ring_count = m.aromatic_ring_count();

  std::pair<int, int> * rscore = new std::pair<int, int>[nr]; unique_ptr<std::pair<int, int>[] > free_rscore(rscore);

  if (do_preference_ranking)
    rank_rings(m, rscore);
  else
  {
    for (auto i = 0; i < nr; ++i)
    {
      rscore[i].first = i;
      rscore[i].second = i;   // not really needed
    }
  }

  int strongly_fused_rings_present = 0;

  if (! remove_strongly_fused_as_a_unit)    // then mark them as "done" so they can be processed separately
  {
    for (auto i = 0; i < nr; ++i)
    {
      const auto ri = m.ringi(i);
      if (ri->strongly_fused_ring_neighbours())
      {
        ring_already_done[i] = MAGIC_NUMBER_INDICATING_STRONGLY_FUSED;          // special number to denote strongly fused rings
        strongly_fused_rings_present++;
      }
    }
  }

// Now we start removing rings one at a time

#ifdef DEBUG_RINT_TRIMMING
  cerr << rings_removed_so_far << ", molecule contains " << matoms << " atoms and " << nr << " rings\n";
#endif

  for (auto i = 0; i < nr; ++i)
  {
    const auto j = rscore[i].first;

    if (ring_already_done[j])
      continue;

    ring_already_done[j] = 1;

    const auto ri = m.ringi(j);

    Molecule mcopy(m);

#ifdef DEBUG_RINT_TRIMMING
    cerr << " i = " << i << ", j = " << j << " trying to remove ring " << *ri << " from " << m.smiles() << endl;
#endif

    atom_number_t first_ring_atom_removed;

    if (ri->strongly_fused_ring_neighbours())
      remove_strongly_fused_ring_group(mcopy, *ri, to_remove, ring_already_done, first_ring_atom_removed);
    else if (! remove_ring(mcopy, m.ring_membership(), *ri, to_remove, ri->is_aromatic(), first_ring_atom_removed))
      continue;

#ifdef DEBUG_RINT_TRIMMING
    cerr << "Ring " << j << " removed " << mcopy.smiles() << endl;
#endif

    identify_exposed_chains(mcopy, to_remove);

    mcopy.remove_atoms(to_remove);

#ifdef DEBUG_RINT_TRIMMING
    cerr << "After removing chains " << mcopy.smiles() << endl;
#endif

    if (! ok_to_keep(mcopy, ri->is_aromatic(), aromatic_ring_count))
    {
       cerr << "NOT OK TO KEEP\n";
      continue;
    }

    const auto my_id = parents.unique_id();

    IWString description;
    if (do_preference_ranking)
      ring_description(m, j, description);
    else
      description = '*';

    do_output(mcopy, parent_id, my_id, rings_removed_so_far + 1, symbol_representing_preference(rscore, nr, i), unique_smiles, first_ring_atom_removed, description, output);

    if (rings_removed_so_far + 1 < max_rings_to_remove && mcopy.nrings() > 1)
      ring_trimming_from_scaffold(mcopy, to_remove, parents, my_id, unique_smiles, rings_removed_so_far+1, output);
  }

// now do the strongly fused rings

  if (strongly_fused_rings_present)
  {
    int * ring_system = new int[matoms + matoms]; unique_ptr<int[]> free_ring_system(ring_system);   // we share this
    m.label_atoms_by_ring_system(ring_system);

    int * ring_membership = ring_system + matoms;
    m.ring_membership(ring_membership);

    for (auto i = 0; i < nr; ++i)
    {
      if (MAGIC_NUMBER_INDICATING_STRONGLY_FUSED != ring_already_done[i])
        continue;

      const Ring & ri = *(m.ringi(i));

      const int lbl = ring_system[ri[0]];   // label applied to this ring

      std::fill_n(to_remove, matoms, 0);

      int system_size = mark_all_other_rings_in_system(m, i, ring_already_done);

      int ndx = 0;
      atom_number_t first_ring_atom_removed;

//    cerr << "Ring " << ri << " part of strongly fused system\n";

      while (next_set_of_atoms_from_strongly_fused_ring(m, ndx, ring_system, lbl, to_remove, first_ring_atom_removed))
      {
        Molecule mcopy(m);
        const auto atoms_removed = mcopy.remove_atoms(to_remove);

        if (! ok_to_keep(mcopy, 0, aromatic_ring_count))
          continue;

        if (! ring_membership_preserved(mcopy))
          continue;

        const auto my_id = parents.unique_id();

        IWString description("S*A");
        description << atoms_removed << ',' << heteroatom_count(m, to_remove) << ",2";

        do_output(mcopy, parent_id, my_id, rings_removed_so_far + 1, '*', unique_smiles, first_ring_atom_removed, description, output);

        if (rings_removed_so_far + 1 < max_rings_to_remove && mcopy.nrings() > 1)
          ring_trimming_from_scaffold(mcopy, to_remove, parents, my_id, unique_smiles, rings_removed_so_far+1, output);
      }
    }
  }

  return 1;
}

static void
assign_atom_status (Molecule & m,
                    int * spinach,
                    int * atom_status)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    if (m.nrings(i))
    {
      atom_status[i] = ATOM_STATUS_RING_ATOM;
      continue;
    }

    if (0 == spinach[i])    // not in a ring, but not part of spinach. Is a chain
    {
      atom_status[i] = ATOM_STATUS_CHAIN;
      continue;
    }

//  Chain atom, currently in the spinach.
//  If singly connected and doubly bonded to a scaffold atom, then this becomes part of the scaffold

    const Atom * a = m.atomi(i);

    const auto acon = a->ncon();

    if (1 != acon)     // not singly connected
      continue;

    const Bond * b = a->item(0);

    if (! b->is_double_bond())
      continue;

    const auto j = b->other(i);

    if (spinach[j])     // atom I is not connected to the scaffold
      continue;

    spinach[i] = 0;

    if (m.nrings(j))
      atom_status[i] = ATOM_STATUS_DOUBLY_BONDED_TO_RING;
    else
      atom_status[i] = ATOM_STATUS_DOUBLY_BONDED_TO_CHAIN;
  }

// For all the atoms that are singly connected, mark their anchor points

  for (auto i = 0; i < matoms; ++i)
  {
    m.set_user_specified_atom_void_ptr(i, atom_status + i);

    if (ATOM_STATUS_DOUBLY_BONDED_TO_RING == atom_status[i] || ATOM_STATUS_DOUBLY_BONDED_TO_CHAIN == atom_status[i])
    {
      const auto j = m.other(i, 0);

      atom_status[j] |= ATOM_STATUS_ROOT_OF_DOUBLE_BOND_ATTACHMENT;
    }
  }

  return;
}

static void
place_isotopes (Molecule & m,
                const int * spinach)
{
  const auto matoms = m.natoms();

  Set_of_Atoms join_points;

  for (auto i = 0; i < matoms; ++i)
  {
    if (spinach[i])    // we only look at ring atoms
      continue;

    if (2 == m.ncon(i))   // nothing branching off here
      continue;

    if (1 == m.nrings(i))
      join_points.add(i);
  }

  if (join_points.size())
    m.set_isotope(join_points, isotope_for_scaffold_join_pts);

  return ;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
ring_trimming (Molecule & m,
               IWString_and_File_Descriptor & output)
{
  const auto nr = m.nrings();

  if (0 == nr)
    return 1;

  if (write_parent_molecule)
    output << m.smiles() << ' ' << m.name() << " PARENT\n";

  auto matoms = m.natoms();

  int * spinach = new_int(matoms + matoms); unique_ptr<int[]> free_spinach(spinach);

  m.identify_spinach(spinach);

  int * atom_status = spinach + matoms;

  assign_atom_status(m, spinach, atom_status);

  if (isotope_for_scaffold_join_pts)
    place_isotopes(m, spinach);

  m.remove_atoms(spinach);

  if (write_molecular_scaffold)
    output << m.smiles() << ' ' << m.name() << " SCAFFOLD " << m.nrings() << '\n';

#ifdef DEBUG_RINT_TRIMMING
  Molecule mcopy(m);
  for (auto i = 0; i < mcopy.natoms(); ++i)
  {
    mcopy.set_isotope(i, i);
  }
  cerr << mcopy.smiles() << " " << m.name() << " SCAFFOLD\n";
#endif

  IW_STL_Hash_Set unique_smiles;

  if (1 == nr)
  {
    IWString description;
    ring_description(m, 0, description);

    return do_output(m, 0, 1, 0, 'F', unique_smiles, 0, description, output);
  }

  Parents parents;

  return ring_trimming_from_scaffold (m, spinach, parents, 0, unique_smiles, 0, output);   // don't need spinach array any more
}

static int
ring_trimming (data_source_and_type<Molecule> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ring_trimming(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
ring_trimming (const char * fname, int input_type, 
                IWString_and_File_Descriptor & output)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return ring_trimming(input, output);
}
static int
ring_trimming (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lw:R:j:J:cm:upayf");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  set_copy_name_in_molecule_copy_constructor(1);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('c'))
  {
    remove_chiral_centres = 1;

    if (verbose)
      cerr << "Chirality removed from all input molecules\n";
  }

  set_copy_atom_based_user_specified_void_pointers_during_add_molecle(1);

  if (cl.option_present('u'))
  {
    suppress_duplicate_structures = 1;

    if (verbose)
      cerr << "Duplicat structures suppressed\n";
  }

  if (cl.option_present('p'))
  {
    do_preference_ranking = 1;

    if (verbose)
      cerr << "Will sort rings by preference\n";
  }

  if (cl.option_present('a'))
  {
    allow_destruction_of_aromaticity = 1;

    if (verbose)
      cerr << "Will create molecules where aromaticity has been destroyed\n";
  }

  if (cl.option_present('y'))
  {
    add_first_removed_atom_number = 1;

    if (verbose)
      cerr << "Will add the atom number of the first atom removed to the output\n";
  }

  if (cl.option_present('f'))
  {
    remove_strongly_fused_as_a_unit = 1;

    if (verbose)
      cerr << "Will remove all atoms in a strongly fused ring system\n";
  }

  if (cl.option_present('w'))
  {
    const_IWSubstring w;

    for (auto i = 0; cl.value('w', w, i); ++i)
    {
      if ("parent" == w)
        write_parent_molecule = 1;
      else if ("rings" == w)
        write_ring_systems = 1;
      else if (w.starts_with("scaf"))
        write_molecular_scaffold = 1;
      else
      {
        cerr << "Unrecognised -w qualifier '" << w << "'\n";
        usage(1);
      }
    }
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_rings_in_system) || max_rings_in_system < 1)
    {
      cerr << "The max rings in a ring system (-R) option must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will only process ring systems with " << max_rings_in_system << " or fewer rings\n";
  }

  if (cl.option_present('j'))
  {
    if (! cl.value('j', isotope_for_scaffold_join_pts) || isotope_for_scaffold_join_pts < 1)
    {
      cerr << "The isotopic label for scaffold attachment points (-j) must be a valid isotope\n";
      usage(2);
    }

    if (verbose)
      cerr << "Scaffold join points labelled with isotope " << isotope_for_scaffold_join_pts << endl;
  }

  if (cl.option_present('J'))
  {
    if (! cl.value('J', isotope_for_ring_fusion) || isotope_for_ring_fusion < 1)
    {
      cerr << "The isotopic label for ring fusion (-J) must be a valid isotope\n";
      usage(2);
    }

    if (verbose)
      cerr << "Broken ring fusion atoms labelled with isotope " << isotope_for_ring_fusion << endl;
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', max_rings_to_remove) || max_rings_to_remove < 1)
    {
      cerr << "The max number of rings to remove (-m) must be a whole +ve integer\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will remove a max of " << max_rings_to_remove << " rings from a ring system\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! ring_trimming(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    for (auto i = 0; i < fragment_count.number_elements(); ++i)
    {
      if (fragment_count[i] > 0)
        cerr << fragment_count[i] << " molecules produced " << i << " separate ring systems\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ring_trimming(argc, argv);

  return rc;
}
