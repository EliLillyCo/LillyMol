/*
 * This program compute the MEDV-13 descriptors
 * Reference:
 * 1. JCICS, 2001, 41, 321-329, Liu et al
 * 2. JCICS, 1999, 39, 951-957, Liu et al 

 This version deals with explicit hydrogens
 */

#include <stdlib.h>
#include <assert.h>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "e_state_computation_procedure_ori.h"

#include "Molecule_Lib/molecule.h"

static int
principal_quantum_number (Molecule &m, const int i)
{
  int atomic_number_i = m.atomic_number(i);

  // H
  if (1==atomic_number_i)
    return 1;
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
  // everything else, unrecognized
  else return -1;
}

static int
Kier_and_Hall_simple_delta_value (Molecule &m, const int i)
{
  return (m.ncon(i) - m.explicit_hydrogens(i));
}

static double
Kier_and_Hall_valence_delta_value (Molecule &m, const int i)
{
  int valance = 0;

  int atomic_number_i = m.atomic_number (i);
  // H
  if (1== atomic_number_i) valance = 1;
  // C
  if (6== atomic_number_i) valance = 4;
  // F Cl Br I
  else if ((9==atomic_number_i) ||(17==atomic_number_i) ||(35==atomic_number_i) || (53==atomic_number_i))  valance = 7;
  // O S
  else if ((8==atomic_number_i) ||(16==atomic_number_i)) valance = 6;
  // N P
  else if ((7==atomic_number_i) ||(15==atomic_number_i)) valance = 5;

  //  cerr<<"atom i="<<i<<"\t valance_delta ="<<(valance - m.hcount(i))<<endl;
  return valance - m.hcount(i);
}

double
value_of_Kier_and_Hall_atom_intrinsic_state (Molecule &m, const int i)
{
  if (1 == m.atomic_number (i)) return 0;
  
  double n = (double) principal_quantum_number (m, i);
  double delta = (double) Kier_and_Hall_simple_delta_value (m, i);
  double delta_v = (double) Kier_and_Hall_valence_delta_value (m, i);

  if (0.0==delta)
    return 0.0;
  else 
    return ((2.0/n)*(2.0/n)*delta_v +1.0)/delta;
}
 
static double
value_of_Kier_and_Hall_relative_electronegativity (Molecule &m, const int i)
{
  if (1==m.atomic_number(i)) return -0.2;

  double n = (double) principal_quantum_number (m, i);
  double delta = (double) Kier_and_Hall_simple_delta_value (m, i);
  double delta_v = (double) Kier_and_Hall_valence_delta_value (m, i);

  return (delta_v - delta) / n /n;
  //  return 1;
}

static int
determine_Kier_and_Hall_atom_intrinsic_state (Molecule &m, double i_state [])
{
  int n_atoms = m.natoms();

  for (int i=0; i<n_atoms; i++)
    i_state[i] = value_of_Kier_and_Hall_atom_intrinsic_state(m, i);

#ifdef DEBUG_DETERMINE_KIER_AND_HALL_ATOM_INTRINSIC_STATE
   cerr << "determine_Kier_and_Hall_atom_intrinsic_state on " << n_atoms << " atoms\n";
   for (int i=0; i<n_atoms; i++)
     cerr<<"atom i="<<i<<"\ti_state ="<<i_state[i]<<endl;
#endif

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
determine_atom_e_state_index (Molecule &m, double e_state_index [])
{
  int n_atoms = m.natoms();
  
  atomic_number_t * z = new atomic_number_t[n_atoms]; std::unique_ptr<atomic_number_t[]> free_z(z);

  m.atomic_numbers(z);

  return determine_atom_e_state_index(m, e_state_index, z);
}

int
determine_atom_e_state_index (Molecule & m, 
                              double e_state_index[], 
                              const atomic_number_t * z)
{
  const int n_atoms = m.natoms();

  double * i_state = new double[n_atoms]; std::unique_ptr<double[]> free_i_state(i_state);

  (void) determine_Kier_and_Hall_atom_intrinsic_state (m, i_state);

  return determine_atom_e_state_index(m, e_state_index, i_state, z);
}

int
determine_atom_e_state_index (Molecule &m, double e_state_index [],
                              const double * i_state,
                              const atomic_number_t * z)
{
  const int n_atoms = m.natoms();

  set_vector (e_state_index, n_atoms, 0.0);

  for (int i=0; i<n_atoms; i++)
  {
    if (1 == z[i])
      continue;

    e_state_index[i] = i_state[i];
    for (int j=0; j<n_atoms; j++)
    {
      if ((j!=i) && 1 != z[j])
      {
        int graph_distance = m.bonds_between (i, j) +1;

        e_state_index[i] += (i_state[i] - i_state[j]) / (graph_distance * graph_distance);
      }
    }
  }

  return 1;
}

int
determine_hydrogen_e_state_index (Molecule &m, double h_e_state_index [],
                                  const double * relative_electronegativity,
                                  const atomic_number_t * z)
{
  int n_atoms = m.natoms();

  for (int i=0; i<n_atoms; i++)
    if (0==m.hcount(i))
      h_e_state_index[i] = 0.0;
    else if (1 != z[i])
      {
        h_e_state_index[i] = -0.2 - relative_electronegativity[i] *2.0;
        for (int j=0; j<n_atoms; j++)
          if ((j!=i) && (1 != z[j]))
            {
              int graph_distance = m.bonds_between (i, j) +1;
              h_e_state_index[i] -=  - (-0.2 - relative_electronegativity[j]) / graph_distance / graph_distance;
            }
      }
  
  for (int i=0; i<n_atoms; i++)
  {
    h_e_state_index[i] = -h_e_state_index[i];
  }

  return 1;
}

int
determine_hydrogen_e_state_index (Molecule & m, 
                                  double h_e_state_index [],
                                  const atomic_number_t * z)
{
  const int n_atoms = m.natoms();

  double * relative_electronegativity = new double[n_atoms]; std::unique_ptr<double[]> free_relative_electronegativity(relative_electronegativity);

  (void) determine_Kier_and_Hall_relative_electronegativity (m, relative_electronegativity);

  return determine_hydrogen_e_state_index(m, h_e_state_index, relative_electronegativity, z);
}

int
number_of_bond_in_conjugated_system (Molecule &m, const int ith_atom)
{
  const Atom * a = m.atomi(ith_atom);
  const int nconi = a->ncon();

  int conjugated_bond = 0;
  for (int j=0; j<nconi; j++)
  {
    const Bond * b = a->item(j);
    const atom_number_t jth_atom = b->other (ith_atom);

    if (b->is_double_bond())
    {
      const Atom * aj = m.atomi(jth_atom);
      int nconj = aj->ncon();

      for (int k=0; k<nconj; k++)
      {
        atom_number_t kth_atom = aj->other(jth_atom, k);
//        if ((kth_atom != ith_atom) && (1 != m.atomic_number(kth_atom)) && (1==m.nbonds(kth_atom) - m.ncon(kth_atom)))
        if ((kth_atom != ith_atom) && ! m.saturated(kth_atom))
        {
          conjugated_bond ++;
          break;
        }
      }
    }
    else if (b->is_single_bond())
    {
//      if ((1 != m.atomic_number(jth_atom)) && (1 == m.nbonds(jth_atom) - m.ncon(jth_atom)))
      if ((1 != m.atomic_number(jth_atom)) && ! m.saturated(jth_atom))
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
    if (IS_TRIPLE_BOND (m.btype_to_connection (ith_atom, j)))
      return 1;

  return 0;
}

