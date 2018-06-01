#include <stdlib.h>
#include <memory>

#include "iwminmax.h"
#include "misc.h"
#include "substructure.h"
#include "target.h"
#include "path.h"
#include "misc2.h"

Substructure_Ring_System_Specification::Substructure_Ring_System_Specification()
{
  _need_per_atom_array = -1;
  
  return;
}

Substructure_Ring_System_Specification::~Substructure_Ring_System_Specification()
{
  if (-5 == _match_as_match_or_rejection)
    cerr << "Deleting already deleted Substructure_Ring_System_Specification\n";

  _match_as_match_or_rejection = -5;

  return;
}

int
Substructure_Ring_System_Specification::ok() const
{
  if (_match_as_match_or_rejection < 0)
    return 0;

  return 1;
}

int
Substructure_Ring_System_Specification::debug_print (std::ostream & os,
                                              const IWString & indentation) const
{
  IWString ind = indentation;
  ind += "  ";

  os << ind << "Details on Substructure_Ring_System_Specification";

  if (0 == _match_as_match_or_rejection)
    os << ", REJECTION";

  os << endl;

  if (_comment.length())
    os << ind << "  " << _comment << endl;

  os << ind << "  Hits needed " << _hits_needed << endl;
  os << ind << "  System size " << _rings_in_system << endl;
  os << ind << "  Ncon        " << _ncon << endl;
  os << ind << "  Fused       " << _degree_of_fusion << endl;
  os << ind << "  Heteroatoms " << _heteroatom_count << endl;
  os << ind << "  Aromatic    " << _aromatic_ring_count << endl;
  os << ind << "  Non arom    " << _non_aromatic_ring_count << endl;
  os << ind << "  Max Fused N " << _largest_number_of_bonds_shared_with_another_ring << endl;
  os << ind << "  Strongly Fs " << _strongly_fused_ring_neighbours << endl;
  os << ind << "  Nbr Spch Gp " << _number_spinach_groups << endl;
  os << ind << "  Nbr NSPCH Gp " << _number_non_spinach_groups << endl;
  os << ind << "  Atm in Spch " << _atoms_in_spinach_group << endl;
  os << ind << "  Len of Spch " << _length_of_spinach_group << endl;
  os << ind << "  Dist to Rng " << _distance_to_another_ring << endl;

  return os.good();
}

int
Substructure_Ring_System_Specification::terse_details (std::ostream & os,
                                      const IWString & indentation) const
{
  IWString ind = indentation;
  ind += "  ";

  os << ind << "Details on Substructure_Ring_System_Specification";

  if (0 == _match_as_match_or_rejection)
    os << ", REJECTION";

  os << endl;

  if (_hits_needed.is_set())
    os << ind << "  Hits needed " << _hits_needed << endl;
  if (_rings_in_system.is_set())
    os << ind << "  Ring size   " << _rings_in_system << endl;
  if (_ncon.is_set())
    os << ind << "  Ncon        " << _ncon << endl;
  if (_degree_of_fusion.is_set())
    os << ind << "  Fused       " << _degree_of_fusion << endl;
  if (_heteroatom_count.is_set())
    os << ind << "  Heteroatoms " << _heteroatom_count << endl;
  if (_aromatic_ring_count.is_set())
    os << ind << "  Aromatic    " << _aromatic_ring_count << endl;
  if (_non_aromatic_ring_count.is_set())
    os << ind << "  Non arom    " << _non_aromatic_ring_count << endl;
  if (_number_spinach_groups.is_set())
    os << ind << "  number spinach groups " << _number_spinach_groups << endl;
  if (_number_non_spinach_groups.is_set())
    os << ind << "  number non-spinach groups " << _number_non_spinach_groups << endl;
  if (_atoms_in_spinach_group.is_set())
    os << ind << "  atoms in spinach " << _atoms_in_spinach_group << endl;
  if (_length_of_spinach_group.is_set())
    os << ind << "  length of spinach " << _length_of_spinach_group << endl;
  if (_distance_to_another_ring.is_set())
    os << ind << "  dist to other ring " << _distance_to_another_ring << endl;

  return os.good();
}

/*
  ATOMS_IN_SYSTEM is an array over target.natoms().
  It will be 1 for each atom in the ring system
*/

int
Substructure_Ring_System_Specification::_check_ncon (const int * atoms_in_system,
                    Molecule_to_Match & target) const
{
  int atoms_in_molecule = target.natoms();

  int ncon = 0;
  int ahc = 0;     // attached heteroatoms(outside the ring system)

  for (int i = 0; i < atoms_in_molecule; i++)
  {
    if (0 == atoms_in_system[i])
      continue;

    const Atom * a = target[i].atom();
    int acon = a->ncon();
    if (2 == acon)    // can only be in 1 ring
      continue;

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);
      if (atoms_in_system[k])    // atom K is in the ring system
        continue;

      ncon++;
      if (! _attached_heteroatom_count.is_set())    // no need to determine heteroatoms
        continue;

      if (6 != target[k].atomic_number())
        ahc++;
    }
  }

  if (_ncon.is_set() && ! _ncon.matches(ncon))
    return 0;

  if (_attached_heteroatom_count.is_set() && ! _attached_heteroatom_count.matches(ahc))
    return 0;

  return 1;
}

static int
compute_atoms_with_pi_electrons (const atom_number_t * in_ring, 
                                 Molecule_to_Match & target)
{
  assert (NULL != in_ring);

  int rc = 0;

  int matoms = target.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (0 == in_ring[i])
      continue;

    Target_Atom & a = target[i];
  
    int acon = a.ncon();

    if (a.nbonds() > acon)     // unsaturation so we assume pi electrons
      rc++;
    else
    {
      Atom * a1 = const_cast<Atom *>(a.atom());

      int pi;
      if (a1->pi_electrons(pi) && pi > 0)
        rc++;
    }
  }

  cerr << "ring system has " << rc << " atoms with pi electrons\n";

  return rc;
}


int
Substructure_Ring_System_Specification::_check_heteroatoms (const int * atoms_in_system,
             Molecule_to_Match & target) const
{
  int atoms_in_molecule = target.natoms();

  int hac = 0;
  for (int i = 0; i < atoms_in_molecule; i++)
  {
    if (atoms_in_system[i] && 6 != target[i].atomic_number())
      hac++;
  }

  if (! _heteroatom_count.matches(hac))
    return 0;

  return 1;
}
              

//#define DEBUG_RING_SYS_MATCHES

int
Substructure_Ring_System_Specification::_matches (Molecule_to_Match & target,
               int * ring_already_done,
               atom_number_t * atoms_in_system)
{
  int nr = target.nrings();

#ifdef DEBUG_RING_SYS_MATCHES
  cerr << "Substructure_Ring_System_Specification::_matches: checking " << nr << " rings in target\n";
#endif

  int nhits = 0;

  int matoms = target.natoms();

// When _all_hits_in_same_fragment is set, we need to keep track of the number
// of hits in each fragment

  extending_resizable_array<int> hits_in_fragment;
  Molecule * m;
  if (_all_hits_in_same_fragment)
    m = target.molecule();
  else
    m = NULL;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = target.ringi(i);

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Starting with ring " << i << " size " << ri->number_elements() << " fused = " << ri->is_fused() << endl;
#endif

//  We have the start of a fused system, build it.
//  Note that if we find that a given ring size causes the system 
//  to be rejected, we must continue scanning all the rings in the
//  system so as to identify them

    int system_rejected = 0;
    int ring_sizes_matched = 0;
    if (_ring_sizes.matches(ri->number_elements()))
      ring_sizes_matched = 1;
    else
      system_rejected = 1;

    int sys_size = 1;
    ring_already_done[i] = i + 1;

    int arom     = ri->is_aromatic();
    int non_arom = ! arom;

    iwmax<int> max_fused(ri->fused_ring_neighbours());
    iwmax<int> largest_number_of_bonds_shared_with_another_ring(ri->largest_number_of_bonds_shared_with_another_ring());
    iwmax<int> strongly_fused_ring_neighbours(ri->strongly_fused_ring_neighbours());

    int rings_with_strongly_fused_neighbours = ri->strongly_fused_ring_neighbours();

    if (NULL != atoms_in_system)
    {
      set_vector(atoms_in_system, matoms, 0);
      ri->set_vector(atoms_in_system, 1);
    }

    if (ri->is_fused())
    {
      for (int j = i + 1; j < nr; j++)
      {
        if (ring_already_done[j])
          continue;

        const Ring * rj = target.ringi(j);
        if (rj->fused_system_identifier() != ri->fused_system_identifier())
          continue;
    
        sys_size++;
        ring_already_done[j] = i + 1;
  
        if (_ring_sizes.matches(rj->number_elements()))
          ring_sizes_matched++;
        else
          system_rejected = 1;
  
        if (NULL != atoms_in_system)
          rj->set_vector(atoms_in_system, 1);

        max_fused.extra(rj->fused_ring_neighbours());
        largest_number_of_bonds_shared_with_another_ring.extra(rj->largest_number_of_bonds_shared_with_another_ring());
        strongly_fused_ring_neighbours.extra(rj->strongly_fused_ring_neighbours());
        if (rj->strongly_fused_ring_neighbours())
          rings_with_strongly_fused_neighbours++;

        if (rj->is_aromatic())
          arom++;
        else
          non_arom++;
      }
    }

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Found system with " << sys_size << " rings, rejected = " << system_rejected << endl;
    cerr << "_rings_that_must_match_ring_sizes? " << _rings_that_must_match_ring_sizes.is_set() << " matched " << ring_sizes_matched << " rings\n";
#endif

    if (! _rings_that_must_match_ring_sizes.is_set())    // every ring must match
    {
      if (system_rejected)
        continue;
    }
    else if (! _rings_that_must_match_ring_sizes.matches(ring_sizes_matched))
      continue;

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Now checking _rings_in_system " << _rings_in_system.matches(sys_size) << endl;
#endif

    if (! _rings_in_system.matches(sys_size))
      continue;

    if (! _degree_of_fusion.matches(max_fused.maxval()))
      continue;

    if (! _aromatic_ring_count.matches(arom))
      continue;

    if (! _non_aromatic_ring_count.matches(non_arom))
      continue;

    if (! _largest_number_of_bonds_shared_with_another_ring.matches(largest_number_of_bonds_shared_with_another_ring.maxval()))
      continue;

    if (! _strongly_fused_ring_neighbours.matches(strongly_fused_ring_neighbours.maxval()))
      continue;

    if (! _strongly_fused_ring_count.matches(rings_with_strongly_fused_neighbours))
      continue;

//  Looking good, count heteroatoms in this system if needed. These
//  cases must be handled specially because they require identification
//  of every atom in the system

    if (_heteroatom_count.is_set())
    {
      if (! _check_heteroatoms(atoms_in_system, target))
        continue;
    }

    if (_ncon.is_set() || _attached_heteroatom_count.is_set())
    {
      if (! _check_ncon(atoms_in_system, target))
        continue;
    }

    if (_atoms_in_system.is_set())
    {
      int ais = count_occurrences_of_item_in_array(1, matoms, atoms_in_system);
      if (! _atoms_in_system.matches(ais))
        continue;
    }

    if (_atoms_with_pi_electrons.is_set())
    {
      int awpe = compute_atoms_with_pi_electrons(atoms_in_system, target);
      if (! _atoms_with_pi_electrons.matches(awpe))
        continue;
    }

#ifdef DEBUG_RING_SYS_MATCHES
    cerr << "Checking environment " << _environment_atom.number_elements() << endl;
#endif

    if (_environment_atom.number_elements())
    {
      if (! _environment_matches(target, atoms_in_system))
        continue;
    }

    if (_number_spinach_groups.is_set() || _atoms_in_spinach_group.is_set() || _length_of_spinach_group.is_set() || _number_non_spinach_groups.is_set())
    {
      if (! _spinach_matches(target, atoms_in_system))
        continue;
    }

    if (_distance_to_another_ring.is_set())
    {
      if (! _match_distance_to_another_ring(target, atoms_in_system))
        continue;
    }

//  We have a ring system which matches!

    nhits++;
    if (_all_hits_in_same_fragment)
    {
      atom_number_t j = ri->item(0);
      hits_in_fragment[m->fragment_membership(j)]++;
    }
  }

#ifdef DEBUG_RING_SYS_MATCHES
  cerr << "Found " << nhits << " matching fused systems, matches nhits = " << _hits_needed.matches(nhits) << endl;
#endif

  if (0 == nhits)
    return ! _match_as_match_or_rejection;

  if (_all_hits_in_same_fragment)
  {
    for (int i = 0; i < hits_in_fragment.number_elements(); i++)
    {
      if (_hits_needed.matches(hits_in_fragment[i]))
        return _match_as_match_or_rejection;
    }

    return ! _match_as_match_or_rejection;
  }

  if (_hits_needed.matches(nhits))
    return _match_as_match_or_rejection;

  return ! _match_as_match_or_rejection;
}

int
Substructure_Ring_System_Specification::matches (Molecule_to_Match & target)
{
  int nr = target.nrings();

  if (0 == nr)
    return ! _match_as_match_or_rejection;

  if (_aromatic_ring_count.is_set() || _non_aromatic_ring_count.is_set())
    target.molecule()->compute_aromaticity_if_needed();

  int * rtmp = new_int(nr); std::unique_ptr<int[]> free_rtmp(rtmp);

// If we will be examining the individual atoms in the system, we need an array for them

  if (_need_per_atom_array < 0)
  {
    if (_heteroatom_count.is_set() || _ncon.is_set() || _atoms_in_system.is_set() || _environment_atom.number_elements() || _number_spinach_groups.is_set() || _atoms_in_spinach_group.is_set() || _length_of_spinach_group.is_set() || _distance_to_another_ring.is_set())
      _need_per_atom_array = 1;
    else
      _need_per_atom_array = 0;
  }

  int * atmp;
  if (_need_per_atom_array)
    atmp = new atom_number_t[target.natoms()];
  else 
    atmp = NULL;

  int rc = _matches(target, rtmp, atmp);

  if (NULL != atmp)
    delete [] atmp;

  return rc;
}

int
Substructure_Ring_System_Specification::_spinach_matches (Molecule_to_Match & target,
                                                          const int * in_system) const
{
  int natoms = target.natoms();

  const Molecule * m = target.molecule();

  int number_spinach_groups = 0;

  int between_ring_connections = 0;

  int got_match_to_atoms_in_group = 0;
  int got_match_to_length_of_spinach = 0;

  for (int i = 0; i < natoms; i++)
  {
    if (! in_system[i])
      continue;

//  cerr << "System contains atom " << i << endl;

    const Atom * a = m->atomi(i);

    int acon = a->ncon();

    if (2 == acon)      // 2 connections, definitely no spinach
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item(j);

      atom_number_t k = b->other(i);

      if (in_system[k])
        continue;

      int spinach = target.is_spinach(k);
//    cerr << " ATom " << k << " spinach " << spinach << endl;

      if (spinach <= 0)
      {
        between_ring_connections++;
        continue;
      }

      if (b->is_double_bond() && 1 == target[k].ncon())   // considered part of the ring
        continue;

      number_spinach_groups++;

      if (_atoms_in_spinach_group.is_set() && _atoms_in_spinach_group.matches(spinach))
        got_match_to_atoms_in_group++;

      if (_length_of_spinach_group.is_set() && ! got_match_to_length_of_spinach)   // expensive so only do it once
        got_match_to_length_of_spinach = _check_length_of_spinach(*m, in_system, i, k);
    }
  }

//cerr << "Checking spinach constraints, spinach " << number_spinach_groups << " between ring " << between_ring_connections << endl;
  if (_number_spinach_groups.is_set() && ! _number_spinach_groups.matches(number_spinach_groups))
    return 0;

  if (_number_non_spinach_groups.is_set() && ! _number_non_spinach_groups.matches(between_ring_connections))
    return 0;

// Maybe we have the case where there

  if (! _atoms_in_spinach_group.is_set())   // nothing to check
    ;
  else if (got_match_to_atoms_in_group)
    ;
  else if (number_spinach_groups > 0)
    return 0;
  else if (_atoms_in_spinach_group.number_elements() > 0)
    return 0;
  else      // if the only specification is a max, and we never got anything, that's OK
  {
    int notused;
    if (_atoms_in_spinach_group.min(notused))   // if min was set, and we didn't test anything, that's a fail
      return 0;
  }

  if (_length_of_spinach_group.is_set() && ! got_match_to_length_of_spinach)
    return 0;

//cerr << "Spinach specification returning match\n";

  return 1;
}

static int 
max_length_of_spinach (const Molecule & m,
                       int * already_done,
                       atom_number_t zatom)
{
  already_done[zatom] = 1;

  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (already_done[j])
      continue;

    int tmp = max_length_of_spinach(m, already_done, j);

    if (tmp > rc)
      rc = tmp;
  }

  return rc + 1;
}

int
Substructure_Ring_System_Specification::_check_length_of_spinach (const Molecule & m,
                                                const int * in_system,
                                                atom_number_t atom_in_ring, 
                                                atom_number_t first_spinach_atom) const
{
  int * already_done = new_int(m.natoms());

  std::unique_ptr<int[]> free_already_done(already_done);

  already_done[atom_in_ring] = 1;

  int maxdist = max_length_of_spinach(m, already_done, first_spinach_atom);

  return _length_of_spinach_group.matches(maxdist);
}

static int
shortest_distance_to_another_ring (Molecule_to_Match & target,
                                   int * already_done,
                                   atom_number_t zatom)
{
#ifdef DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING
  cerr << "shortest_distance_to_another_ring atom " << zatom << endl;
#endif

  already_done[zatom] = 1;

  const Molecule * m = target.molecule();

  int rc = 0;

  const Atom * a = m->atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (already_done[j])
      continue;

    if (target[j].is_ring_atom())
    {
      rc = 1;
      continue;
    }

    if (1 == target[j].ncon())
      continue;

    if (target.is_spinach(j))
      continue;

    int tmp = shortest_distance_to_another_ring(target, already_done, j);

    if (0 == rc)
      rc = tmp;
    else if (tmp < rc)
      rc = tmp;
  }

  return 1 + rc;
}

//#define DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING

int
Substructure_Ring_System_Specification::_match_distance_to_another_ring (Molecule_to_Match & target,
                                                const int * in_ring_system) const
{
  const int matoms = target.natoms();

  int * already_done = new int[matoms]; std::unique_ptr<int[]> free_already_done(already_done);

  const Molecule * m = target.molecule();

  int rc = 0;

#ifdef DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING
  cerr << "_match_distance_to_another_ring:checking " << matoms << " atoms\n";
#endif

  for (int i = 0; i < matoms; i++)
  {
    if (! in_ring_system[i])
      continue;

    const Atom * a = m->atomi(i);

    int acon = a->ncon();

    if (2 == acon)     // no branches outside the ring
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item(j);

      const atom_number_t k = b->other(i);

      if (in_ring_system[k])
        continue;
      
      int d;
      if (0 == b->nrings() && target[k].nrings())
        d = 1;
      else
      {
        if (! target.is_between_rings(k))
          continue;

        set_vector(already_done, matoms, 0);

        already_done[i] = 1;

        d = shortest_distance_to_another_ring(target, already_done, k);
      }

#ifdef DEBUG_SHORTEST_DISTANCE_TO_ANOTHER_RING
      cerr << "From atom " << i << " shortest dist to another ring " << d << endl;
#endif

      const auto mm = _distance_to_another_ring.matches(d);

      if (0 == _distance_to_another_ring.number_elements())    // we have a min and/or max specification only
      {
        if (! mm)         // did not match, violated constraint, we are done
          return 0;

        rc = 1;    // we found a ring that did NOT violate the constraint
      }
      else if (mm)                  // found a specific distance that matches
        return 1;

      if (3 == acon)     // just one bond outside the ring, hopefully the most common case
        break;
    }
  }

  return rc;
}
