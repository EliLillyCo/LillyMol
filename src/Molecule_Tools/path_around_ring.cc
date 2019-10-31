#include <stdlib.h>
#include <memory>
using namespace std;

#include "misc.h"

#include "molecule.h"
#include "path.h"
#include "path_around_ring.h"

//#define DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM

int
path_around_edge_of_ring_system (Molecule & m,
                                 const int * process_these_atoms,
                                 int flag,
                                 Set_of_Atoms & s)
{
  int matoms = m.natoms();

  if (s.number_elements())
    s.resize_keep_storage(0);

  const int * ring_membership = m.ring_membership();

#ifdef DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM
  for (int i = 0; i < matoms; i++)
  {
    if (flag == process_these_atoms[i])
      cerr << "Atom " << i << " type " << m.atomic_symbol(i) << " in " << ring_membership[i] << " rings\n";
  }
#endif

// The start atom must be one in just one ring

  atom_number_t start_atom = INVALID_ATOM_NUMBER;

  for (int i = 0; i < matoms; i++)
  {
    if (flag != process_these_atoms[i])
      continue;

    if (INVALID_ATOM_NUMBER != start_atom)   // already found a starting point
      continue;

    if (1 != ring_membership[i])
      continue;

    start_atom = i;
  }

  if (INVALID_ATOM_NUMBER == start_atom)
  {
//  cerr << "path_around_edge_of_ring_system::no start atom in '" << m.name() << "'\n";
    return 0;
  }

  int * already_visited = new_int(matoms); std::unique_ptr<int[]> free_already_visited(already_visited);

  already_visited[start_atom] = 1;
  s.add(start_atom);

// If there are any non sssr rings that touch these atoms, don't continue

  int n = m.non_sssr_rings();

  for (int i = 0; i < n; i++)
  {
    const Ring * r = m.non_sssr_ring(i);
    if (r->any_members_set_in_array(process_these_atoms, flag))
      return 0;
  }

#ifdef DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM
  cerr << "First atom in traversal " << s[0] << endl;
#endif

  atom_number_t current_atom = s[0];

  atom_number_t first_atom = current_atom;

  while (s.number_elements() <= matoms + 1)
  {
    const Atom * a = m.atomi(current_atom);

    int acon = a->ncon();

#ifdef DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM
    cerr << "Current atom " << current_atom << ", ncon " << acon << endl;
#endif

    atom_number_t next_atom = INVALID_ATOM_NUMBER;

    for (int i = 0; i < acon; i++)
    {
      const Bond * b = a->item(i);

#ifdef DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM
      cerr << "From " << current_atom << " can we go to to atom " << b->other(current_atom) << " bond in " << b->nrings() << " rings , flag? " << (flag == process_these_atoms[b->other(current_atom)]) << endl;
#endif

      atom_number_t j = b->other(current_atom);

      if (flag != process_these_atoms[j])
        continue;

      if (b->nrings() > 2)   // either non-planar ring, or 3 rings intersecting
        continue;

      if (2 == b->nrings())    // fused ring, not going around the outside
        continue;

#ifdef DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM
      cerr << "Check first atom " << first_atom << " vs " << j << " have " << s.number_elements() << " items\n";
#endif
      if (first_atom == j && s.number_elements() > 2)    // we are done!
        return s.number_elements();

      if (already_visited[j])
        continue;

//    if (ring_membership[j] > 2)
//      return 0;

      if (INVALID_ATOM_NUMBER == next_atom)
        next_atom = j;
      else if (s.number_elements() > 1)   // first step will have two choices
        return 0;
    }

    if (INVALID_ATOM_NUMBER == next_atom)
      return 0;

#ifdef DEBUG_PATH_AROUND_EDGE_OF_RING_SYSTEM
      cerr << "path_around_edge_of_ring_system '" << m.name() << " continuing from " << current_atom << " to " << next_atom << endl; //" '" << m.smarts_equivalent_for_atom (current_atom) << "'\n";
#endif

    s.add(next_atom);
    already_visited[next_atom] = 1;
    current_atom = next_atom;
  }

  cerr << "path_around_edge_of_ring_system:no termination\n";
  cerr << s << endl;

  assert (NULL == "Should not come to here");

  return 0;
}
