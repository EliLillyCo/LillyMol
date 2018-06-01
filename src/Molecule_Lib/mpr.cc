#include <stdlib.h>
#include <limits>
#include <memory>


#include "molecule.h"
#include "path.h"

#include "mpr.h"

Molecular_Properties_Generator::Molecular_Properties_Generator ()
{
  _unsaturation_includes_aromatic = 0;

  return;
}

template <typename T>
int
Molecular_Properties_Generator::_molecular_properties_generation (Molecule & m,
                                 const int * ring_membership,
                                 T * properties) const
{
  int matoms = m.natoms();

  if (matoms > std::numeric_limits<T>::max())    // check for overflow. Note that if natoms overflows, so will many of the others. We don't check this!!! - too lazy
  {
    properties[0] = std::numeric_limits<T>::max();
  }
  else
    properties[0] = matoms;

  int nr = m.nrings();

  properties[1] = 0;
  properties[2] = nr;

  if (nr)
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * r = m.ringi(i);
      if (r->number_elements() > properties[1])
        properties[1] = r->number_elements();
    }
  }

  int ring_atoms = 0;
  int aromatic_atoms = 0;
  int fused_ring_atoms = 0;

  if (nr)    // the molecule has rings
  {
    assert (NULL != ring_membership);

    for (int i = 0; i < matoms; i++)
    {
      if (0 == ring_membership[i])
        continue;

      ring_atoms++;
      if (m.is_aromatic(i))
        aromatic_atoms++;

      if (nr > 1 && ring_membership[i] > 1)
        fused_ring_atoms++;
    }
  }
  
  properties[3] = ring_atoms;
  properties[4] = aromatic_atoms;
  properties[5] = fused_ring_atoms;

  int heteroatom_count = 0;
  int unsaturation = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (6 != a->atomic_number())
      heteroatom_count++;

    if (a->nbonds() > a->ncon())
    {
      if (_unsaturation_includes_aromatic)
        unsaturation++;
      else if (nr && ring_membership[i] && m.is_aromatic(i))
        ;
      else
        unsaturation++;
    }
  }

  properties[6] = heteroatom_count;
  properties[7] = unsaturation;

//#define DEBUG_PROPERTIES
#ifdef DEBUG_PROPERTIES
  cerr << "Molecule '" << m.name() << "' properties";
  for (int i = 0; i < NPROPERTIES; i++)
  {
    cerr << ' ' << int(properties[i]);
  }

  cerr << endl;
#endif

  return 1;
}


template <typename T>
int
Molecular_Properties_Generator::operator() (Molecule & m,
                                 T * properties) const
{
  int nr = m.nrings();

  if (0 == nr)
    return _molecular_properties_generation(m, NULL, properties);

  int * ring_membership = new int[m.natoms()]; std::unique_ptr<int[]> free_ring_membership(ring_membership);

  m.ring_membership_including_non_sssr_rings(ring_membership);

#ifdef ECHO_RING_MEMBERSHIP
  for (int i = 0; i < m.natoms(); i++)
  {
    if (ring_membership[i])
      cerr << "Atom " << i << " in " << ring_membership[i] << " rings, '" << m.smarts_equivalent_for_atom(i) << "'\n";
  }
#endif

  return _molecular_properties_generation(m, ring_membership, properties);
}

template int Molecular_Properties_Generator::operator()(Molecule &, unsigned char *) const;
template int Molecular_Properties_Generator::operator()(Molecule &, int *) const;
