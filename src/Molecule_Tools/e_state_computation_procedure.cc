/*
 * This program compute the MEDV-13 descriptors
 * Reference:
 * 1. JCICS, 2001, 41, 321-329, Liu et al
 * 2. JCICS, 1999, 39, 951-957, Liu et al 

 It looks as if this version assumes no explicit Hydrogen atoms.
 */

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <assert.h>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/molecule.h"
#include "e_state_computation_procedure.h"

static int
principal_quantum_number (const atomic_number_t atomic_number_i)
{
  // C N O F
  if ((5<atomic_number_i) && (10>atomic_number_i))
    return 2;
  // P S Cl
  else if ((14<atomic_number_i) && (18>atomic_number_i))
    return 3;
  // Br
  else if (35 == atomic_number_i)
    return 4;
  // I
  else if (53 == atomic_number_i)
    return 5;
  if (1==atomic_number_i)
    return 1;
  // everything else, unrecognized
  else return -1;
}

static int
principal_quantum_number (Molecule &m, const int i)
{
  int atomic_number_i = m.atomic_number(i);

  return principal_quantum_number(atomic_number_i);
}

static int
Kier_and_Hall_simple_delta_value (const Atom * a)
{
  return a->ncon();
}

static double
Kier_and_Hall_valence_delta_value (atomic_number_t atomic_number_i,
                                   const int hcount)
{
  int valance = 0;

  // C
  if (6== atomic_number_i) valance = 4;
  // F Cl Br I
  else if ((9==atomic_number_i) ||(17==atomic_number_i) ||(35==atomic_number_i) || (53==atomic_number_i))  valance = 7;
  // O S
  else if ((8==atomic_number_i) ||(16==atomic_number_i)) valance = 6;
  // N P
  else if ((7==atomic_number_i) ||(15==atomic_number_i)) valance = 5;
  // H
  else if (1== atomic_number_i) valance = 1;

  //  cerr<<"atom i="<<i<<"\t valance_delta ="<<(valance - m.hcount(i))<<endl;
  return valance - hcount;
}

double
value_of_Kier_and_Hall_atom_intrinsic_state (Molecule &m, const int i,
                                             const Atom * const * atoms,
                                             const atomic_number_t * z,
                                             const int * hcount)
{
  double delta = static_cast<double>(Kier_and_Hall_simple_delta_value(atoms[i]));

  if (0.0 == delta)
    return 0.0;

  double n = static_cast<double>(principal_quantum_number(z[i]));

  double delta_v = Kier_and_Hall_valence_delta_value (z[i], hcount[i]);

  return ((2.0/n)*(2.0/n) * delta_v + 1.0)/delta;
}
 
static double
value_of_Kier_and_Hall_relative_electronegativity (Molecule &m, const int i)
{
  if (1==m.atomic_number(i)) return -0.2;

  double n = (double) principal_quantum_number (m, i);
  double delta = (double) Kier_and_Hall_simple_delta_value (m.atomi (i));
  double delta_v = (double) Kier_and_Hall_valence_delta_value (m.atomic_number (i), m.hcount (i));

  return (delta_v - delta) / n /n;
  //  return 1;
}

static int
determine_Kier_and_Hall_atom_intrinsic_state (Molecule &m, double i_state [],
                                const Atom * const * atoms,
                                const atomic_number_t * z,
                                const int * hcount)
{
  int n_atoms = m.natoms();

  for (int i=0; i<n_atoms; i++)
    i_state[i] = value_of_Kier_and_Hall_atom_intrinsic_state(m, i, atoms, z, hcount);

//#define DEBUG_INTRINSIC_STATE
#ifdef DEBUG_INTRINSIC_STATE
  for (int i=0; i<n_atoms; i++)
    cerr<<"atom i="<<i<<"\ti_state ="<<i_state[i]<<endl;
#endif
  
  //  cerr<<endl<<endl;

  return 1;
}

static int 
determine_Kier_and_Hall_relative_electronegativity (Molecule &m, double relative_electronegativity [])
{
  int n_atoms = m.natoms ();

  for (int i=0; i<n_atoms; i++)
    relative_electronegativity[i] = value_of_Kier_and_Hall_relative_electronegativity (m, i);

  return 1;
}

/*
 * determine the E state index from intrinsic state
 */

int
determine_atom_e_state_index (Molecule &m,
                              double e_state_index [],
                              const Atom * const * atoms,
                              const atomic_number_t * z,
                              const int * hcount)
{
  int n_atoms = m.natoms();
  
  double * i_state = new double[n_atoms]; std::unique_ptr<double[]> free_i_state(i_state);

  (void) determine_Kier_and_Hall_atom_intrinsic_state (m, i_state, atoms, z, hcount);

  set_vector (e_state_index, n_atoms, 0.0);

  for (int i=0; i<n_atoms; i++)
  {
	e_state_index[i] = i_state[i];
//      cerr << "Atom " << i << " type " << m.atomic_symbol(i) << " e_state_index from istate " << e_state_index[i] << endl;
	for (int j=0; j<n_atoms; j++)
        {
	  if (j!=i)
	  {
	    int graph_distance = m.bonds_between (i, j) +1;

	    e_state_index[i] += (i_state[i] - i_state[j]) / (graph_distance * graph_distance);
//          cerr << "Atom " << i << " to " << j << ", i_state " << i_state[j] << " dist " << graph_distance << " value now " << e_state_index[i] << endl;
	  }
        }
      }

//#define DEBUG_ATOM_E_STATE_INDEX
#ifdef DEBUG_ATOM_E_STATE_INDEX
      for (int i = 0; i < n_atoms; i++)
      {
        cerr << "Atom " << i << " type " << m.atomic_symbol(i) << " e_state_index " << e_state_index[i] << endl;
      }
#endif

  return 1;
}

int
determine_atom_e_state_index (Molecule & m,
                              double e_state_index [])
{
  int matoms = m.natoms();

  atomic_number_t * z = new atomic_number_t[matoms]; std::unique_ptr<atomic_number_t[]> free_z(z);
  m.atomic_numbers(z);

  int * hcount = new int[matoms]; std::unique_ptr<int[]> free_hcount(hcount);

  for (int i = 0; i < matoms; i++)
  {
    hcount[i] = m.hcount(i);
  }

  const Atom * * atoms = new const Atom *[matoms];
  m.atoms(atoms);

  int rc = determine_atom_e_state_index(m, e_state_index, atoms, z, hcount);

  delete [] atoms;

  return rc;
}

int
determine_hydrogen_e_state_index (Molecule &m, double h_e_state_index [],
                                  const atomic_number_t * z)
{
  int n_atoms = m.natoms();
  double * relative_electronegativity = new double[n_atoms]; std::unique_ptr<double[]> free_relative_electronegativity(relative_electronegativity);
  (void) determine_Kier_and_Hall_relative_electronegativity (m, relative_electronegativity);

  //  for (int i=0; i<n_atoms; i++)
  //    cerr<<"atom i="<<i<<"\t ele="<<relative_electronegativity[i]<<endl;
  //  cerr<<endl<<endl;

  for (int i=0; i<n_atoms; i++)
  {
    if (1 == z[i])
      continue;

    if (0==m.hcount(i))
      h_e_state_index[i] = 0.0;
    else
    {
      h_e_state_index[i] = -0.2 - relative_electronegativity[i] *2.0;
      for (int j=0; j<n_atoms; j++)
        if ((j!=i) && (1 != z[j]))
          {
            int graph_distance = m.bonds_between (i, j) +1;

            h_e_state_index[i] -= - (-0.2 - relative_electronegativity[j]) / (graph_distance * graph_distance);
          }
      }
  }
  
  for (int i=0; i<n_atoms; i++)
  {
//  if (0.0 != h_e_state_index[i])
      h_e_state_index[i] = -h_e_state_index[i];
  }

  return 1;
}

int
determine_hydrogen_e_state_index (Molecule & m,
                                  double h_e_state_index [])
{
  const int n_atoms = m.natoms();

  atomic_number_t * z = new atomic_number_t[n_atoms]; std::unique_ptr<atomic_number_t[]> free_z(z);

  m.atomic_numbers(z);

  return determine_hydrogen_e_state_index(m, h_e_state_index, z);
}

int
number_of_bond_in_conjugated_system (Molecule &m, const int ith_atom)
{
  int conjugated_bond = 0;
  int nconi = m.ncon(ith_atom);
  for (int j=0; j<nconi; j++)
    {
      if (2 == (m.btype_to_connection(ith_atom, j) % 8))
        {
          atom_number_t jth_atom = m.other (ith_atom, j);
          int nconj = m.ncon(jth_atom);
            for (int k=0; k<nconj; k++)
              {
                atom_number_t kth_atom = m.other(jth_atom, k);
                if ((kth_atom != ith_atom) && (1 != m.atomic_number(kth_atom)) && (1==m.nbonds(kth_atom) - m.ncon(kth_atom)))
                  {
                    conjugated_bond ++;
                    break;
                  }
              }
        }
      else if (1 == (m.btype_to_connection(ith_atom, j) % 8))
        {
          atom_number_t jth_atom = m.other (ith_atom, j);
          if ((1 != m.atomic_number(jth_atom)) && (1 == m.nbonds(jth_atom) - m.ncon(jth_atom)))
            conjugated_bond ++;
        }
    }
  return conjugated_bond;
}

int
is_triple_bonded (Molecule &m, const int ith_atom)
{
  int nconi = m.ncon(ith_atom);
  for (int j=0; j<nconi; j++)
    if (3 == (m.btype_to_connection (ith_atom, j) % 8))
      return 1;
  return 0;
}

