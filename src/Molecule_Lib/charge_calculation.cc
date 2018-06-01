/*
 * Calculate partial charge for molecules
 * Method 1, Abraham's method
 * The last paper for partial charge is Journal of Computational Chemistry (1992) Vol 13 492-504
 * Method 2, Gasteiger's method
 * Reference Tetrahedron, 36, 3219-3228 (1980)
 * Method 3, Huckel's method
 * Reference Streitwieser, Molecular Orbital Theory for Organic Chemists, 1961
 *           Wilkinson, The Algebraic Eigenvalue Problem, 1967
 */

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <memory>

#include <iomanip>
#include <iostream>

using namespace std;

#include "misc.h"

#include "charge_calculation.h"
#include "path.h"
#include "molecule.h"
#include "aromatic.h"

//DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

// define the values for huckel charge calculation
static const double huckel_parameters [NUMBER_OF_ATOM_TYPE] = 
{ 0, 2.0, 0.0, 0, 0.0, 0.3, 2.0, 0.4, 1.0, 0.4, 1.6, 1.1, 2.0, 0.7, 2.0, 0.35, -1.0, 2.2, 2.4, 2.6, 2.0 , 0.35, 0.35, 0.7, -1};

// defube the number of electrons
static const double electron_parameters [NUMBER_OF_ATOM_TYPE] =
{ 0, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2,1,0};
// the following are defined for the partial charge calculation (Abraham's method)
// first defined are the constants

static const int JURS_NUMBER_OF_ATOMS =17;
static const int HYDROGEN = 1, CARBON = 6, NITROGEN =7, OXYGEN =8, FLURINE=9, CLORINE=17, BROMINE=35, IODINE=53, SULFUR=16, PHOSPHOR=15;
static const int H_ATOM=0, CSP3_ATOM=1, CSP2_ATOM=2, CSP1_ATOM=3, NSP3_t_ATOM=4, NSP2_ATOM=5, NSP1_ATOM=6, OSP3_ATOM=7, OSP2_ATOM=8, F_ATOM=9,CL_ATOM=10, BR_ATOM=11, I_ATOM=12, SSP3_ATOM=13, SSP2_ATOM=14, P_ATOM=15, NSP3_p_ATOM=16;
static const double constant_a1=40.4; // for bonds between second row atoms
static const double constant_a2=23.8; // for bonds involving Hydrogen
static const double constant_a3=21.0; // for all other bonds (all except above or below)
static const double constant_a4=30.8; // for bond between N & O in NO2 group and N-C(SP2) in NO2-C= and N-C=O in amide
static const double constant_a5=13.6; // for S-H bond in R-S-H  and O=S-H

static const double constant_b=198.4; // constant used in partial charge calculation for two-bond effect
static const double constant_c=3.0; // constant used in polarizability calculation
static const double constant_d=5.0; // constant used in partial charge calculation for the three-bond effect

static double MAX_CHARGE_DIFFERENCE = 1e-4; // constant used to exit the do while loop in function calculate_two_and_three_bond_partial_charge() 

static int MAX_LOOP = 200; // while balancing the charge, at most 200 loop is run

static const double electronegativity[JURS_NUMBER_OF_ATOMS]={7.17, 7.98, 8.79, 10.38, 11.54, 12.87, 15.68, 14.20, 17.07, 16.96, 11.85, 11.60, 11, 10.21, 10.21, 8.90, 12.32};
static const double polarizability[JURS_NUMBER_OF_ATOMS]={0.628, 1.056, 1.166, 1.137, 1.085, 0.954, 0.986, 0.912, 0.860, 0.697, 1.518,1.885,2.423, 1.748, 1.532, 1.549, 0.922};
static const double std_state_charge[JURS_NUMBER_OF_ATOMS]={0.034, -0.068, -0.126, -0.238, -0.264, -0.242, -0.191, -0.308, -0.225, -0.222, -0.184, -0.172, -0.144, -0.212, -0.212, -0.131, -0.322 };
// end of constants and arrays for the partial charge calculation

// The following are defined for the partial charge calculation (Gasteiger's method) sigma charge
static const double gasteiger_constant_a [NUMBER_OF_ATOM_TYPE] =                       {7.17, 7.98,  8.79, 8.79, 10.39,  0.00, 11.54, 12.87, 15.68, 12.87, 12.32, 12.32, 14.18, 17.07, 10.14, 10.88,  8.9, 10.08, 11.00, 14.66, 9.90, 10.14, 12.00, 17.07,  8.79};
static const double gasteiger_constant_b [NUMBER_OF_ATOM_TYPE] =                       {6.24, 9.18,  9.32, 9.32,  9.45, 11.86, 10.82, 11.15, 11.70, 11.15, 11.20, 11.20, 12.92, 13.79,  9.13,  9.49,  8.24, 8.47,  9.69, 13.85, 7.96,  9.13, 10.80, 13.79,  9.32};
static const double gasteiger_constant_c [NUMBER_OF_ATOM_TYPE] =                       {-0.56,1.88,  1.51, 1.51,  0.73, 11.86,  1.36,  0.85, -0.27,  0.85,  1.34,  1.34,  1.39,  0.47,  1.38,  1.33,  0.96, 1.16,  1.35,  2.31, 0.96,  1.38,  1.20,  0.47,  1.51};

// The following are defined for the partial charge calculation (Gasteiger's method)  pi charges
static const double gasteiger_pi_constant_a [NUMBER_OF_ATOM_TYPE] = {0,0,5.6,5.6,5.64,0,4.54,7.95,7.92,7.95,3.57,3.57,7.91,10.09,6.6,7.73,0,5.2,6.5,7.34,8.81,6.6,6.6,10.09,5.6};
static const double gasteiger_pi_constant_b [NUMBER_OF_ATOM_TYPE] = {0,0,8.93,8.93,8.4,11.86,11.86,9.73,9.775,9.73,10.17,10.17,14.76,11.73,10.32,8.16,0,9.68,9.69,13.86,8.81,10.32,10.32,11.73,8.93};
static const double gasteiger_pi_constant_c [NUMBER_OF_ATOM_TYPE] = {0,0,2.94,2.94,2.81,11.86,7.32,2.67,2.685,2.67,6.6,6.6,6.85,2.87,3.72,1.81,0,4.48,5.49,9.64,0,3.72,3.72,2.87,2.94};
// For everything except hydrogen, cation's electronegativity is the additon of constant a, b and c
// for hydrogen it is 20.02
static const double gasteiger_element_cation_electronegativity [NUMBER_OF_ATOM_TYPE] = {20.02,19.04,19.62,19.62, 20.57, 23.72, 23.72, 24.87, 27.11, 24.87, 24.86, 24.86, 28.49, 31.33, 20.65, 21.70, 18.1, 30.82, 22.04, 19.71,18.82,  20.65,24.00, 31.33, 19.62};
// end of constant for Gasteiger's partial charge 
/*
 * This funciton find out the a constant to be used in the partial_charge calculation (except special cases)
 * If one of the atom is hydrogen (atomic_number = 1) return a2
 * If both of them is second row atom (atomic_number < 11) return a1
 * For the rest of the situations, return a3
 */

double 
Molecule::__a_const(atomic_number_t ith_atom_atomic_number, atomic_number_t jth_atom_atomic_number)
{
  if((ith_atom_atomic_number == 1) || (jth_atom_atomic_number == 1)) 
    {
      return constant_a2;
    }
  else 
    {
      if ((ith_atom_atomic_number <11) && (jth_atom_atomic_number < 11)) 
        {
          return constant_a1;
        }
      else return constant_a3;
    }
}

/*
 * This funciton find out the a constant to be used in the partial_charge calculation (including special cases)
 * If one of the atom is hydrogen (atomic_number = 1) return a2
 * If both of them is second row atom (atomic_number < 11) return a1
 * Special conditions return a4 & a5  -- detail in comment below
 * For the rest of the situations, return a3
 */

double 
Molecule::_a_const (atom_number_t i, atom_number_t j)
{
  const Atom * ith_atom = _things[i];
  const Atom * jth_atom = _things[j];

  int ith_atom_atomic_number = ith_atom->atomic_number();
  int jth_atom_atomic_number = jth_atom->atomic_number();

  // ith_atom_atomic number is always <= jth_atom_atomic number 
  if (ith_atom_atomic_number>jth_atom_atomic_number)
    {
      const Atom * tmp_atom;

      int tmp=ith_atom_atomic_number;
      tmp_atom=ith_atom;
 
      ith_atom_atomic_number=jth_atom_atomic_number;
      ith_atom=jth_atom;
      
      jth_atom_atomic_number=tmp;
      jth_atom=tmp_atom;
    }

  // if it is SH bond, return constant a5;
  if ((jth_atom_atomic_number == SULFUR)&&(ith_atom_atomic_number==HYDROGEN)) 
    return constant_a5;

  // if it is NO bond in NO2 group, return constant a4
  if ((ith_atom_atomic_number==NITROGEN)&&(jth_atom_atomic_number==OXYGEN))
    if (ith_atom->nbonds()>3) return constant_a4;

  // if it is NC bond in NO2C return constant a4 
  if ((ith_atom_atomic_number==CARBON)&&(jth_atom_atomic_number==NITROGEN))
    if ((jth_atom->nbonds()>3)&&(ith_atom->nbonds() != ith_atom->ncon())) return constant_a4;

  // if it is NC bond in -NC(=O)- return constant a4
  if ((ith_atom_atomic_number==CARBON)&&(jth_atom_atomic_number==NITROGEN))
    if ((jth_atom->nbonds() - jth_atom->ncon()) ==1) return constant_a4;

  return __a_const (ith_atom_atomic_number, jth_atom_atomic_number);
}

/* 
 * Find the index so that program can locate the parameters to use for
 * the abrahms partial charge calculation In order to use it for the
 * gasteiger charge calculation, special attention is used to treat
 * special return from this funciton in gaister_row_number (Molecule
 * &m, atom_number_t ith_atom)
 */

int 
Molecule::_row_number(atom_number_t ith_atom)
{
  const Atom * atom = _things[ith_atom];
#if (DEBUG_SWITCH)
  cerr<<"atomic_number:   "<<atom->atomic_number()<<"  ncon:   "<<atom->ncon()<<"  nbond:   "<<atom->nbonds()<<endl;
#endif
  //
  switch (atom->atomic_number()) {
  case HYDROGEN: return H_ATOM;
  case CARBON:
    switch (atom->nbonds()-atom->ncon()) {
    case 2: 
#if (DEBUG_SWITCH) 
      cerr<<"SP Carbon"<<endl;
#endif
      return CSP1_ATOM;
      break;
    case 1:
#if (DEBUG_SWITCH) 
      cerr<<"SP2 Carbon"<<endl;
#endif
      return CSP2_ATOM;
      break;
    default: 
#if (DEBUG_SWITCH) 
      cerr<<"SP3 Carbon"<<endl;
#endif
      return CSP3_ATOM;
    }
    break;
  case NITROGEN:
    switch (atom->nbonds()-atom->ncon()) {
    case 2:
#if (DEBUG_SWITCH) 
      cerr<<"SP1 Nigrogen"<<endl;
#endif
      if (atom->ncon()==1) return NSP1_ATOM;   // N#C group
      else return NSP3_p_ATOM;  // NO2 group
      break;
    case 1:
#if (DEBUG_SWITCH) 
      cerr<<"SP2 Nigrogen"<<endl;
#endif
      return NSP2_ATOM;    // N=C group
      break;
    default:
      int planar_atom = 0;
      int tmp_connection = atom->ncon();
      for (int i=0; i<tmp_connection; i++) 
        {
          const Atom *atomj = _things[other(ith_atom, i)];
          if (atomj->nbonds()!= atomj->ncon()) 
            {
              planar_atom = 1;
              break;
            }
        }
#if (DEBUG_SWITCH) 
      if (planar_atom) cerr<<"Planar N Atom"<<endl;
      else cerr<<"Not planar N Atom"<<endl;
#endif
      if (planar_atom) return NSP3_p_ATOM;
      else return NSP3_t_ATOM;
    }
    break;
  case OXYGEN: 
    if (atom->ncon()==atom->nbonds()) return OSP3_ATOM; else return OSP2_ATOM;
    break;
  case FLURINE: return F_ATOM;
    break;
  case CLORINE: return CL_ATOM;
    break;
  case BROMINE: return BR_ATOM;
    break;
  case IODINE: return I_ATOM;
    break;
  case SULFUR:
#if (DEBUG_SWITCH) 
    cerr<<"Sulfur Atom"<<endl;
#endif
    if (atom->ncon()==atom->nbonds()) return SSP3_ATOM; else return SSP2_ATOM;
    break;
  case PHOSPHOR: 
#if (DEBUG_SWITCH) 
    cerr<<"Phosphor Atom"<<endl;
#endif
    return P_ATOM;
  default: return OUT_OF_RANGE;
  }
}

/*
 * This funciton find neighbors (one, two and three bond) for all atoms of molecules
 */

int
Molecule::_find_neighbors_for_molecule(arrays_holder neighbors[])
{

#if (DEBUG_SWITCH) 
  cerr<<"I am inside find_neighbors_for_molecule (array_holder) now"<<endl;
#endif
  int rc=1;
  int n_atom = natoms();

  // find the one_bond_neighbors for all atoms of molecule
  for (atom_number_t i=0;i<n_atom;i++)
    { 
      // resize the one_bond_neighbors for each atom of each molecule
      neighbors[i].one_bond_neighbors.resize(0);

      Set_of_Atoms temp_array;
      rc = connections(i, temp_array);

      //do a union -- in one_bond case, only one array
      neighbors[i].one_bond_neighbors += temp_array;

      //resize the array to release the memory
      temp_array.resize(0);
      
      // printout the one_bond neighbors for the molecule
#if (DEBUG_SWITCH)
      cerr<<"Inside find_neighbors_for_molecule() -- one_bond_neighbor"<<endl;
      cerr<<"all one-bond neighbors   "<<i<<endl;
      cerr<<"     Atom     number"<<endl;
      for (int kk=0;kk<neighbors[i].one_bond_neighbors.number_elements();kk++)
      {
        cerr<<"      "<<atomi(neighbors[i].one_bond_neighbors.item(kk))->atomic_symbol() <<"      "<<neighbors[i].one_bond_neighbors.item(kk)<<endl;
      }
      cerr<<endl<<endl<<endl;
#endif
    }

  // find two_bond neighbors for all atoms of the molecule
  for (atom_number_t i=0;i<n_atom;i++)
    {
      // The idea here is that the two_bond neighbors of an atom is the union of the
      // one_bond neighbors of its one_bond neighbors
 
      int number_of_one_bond_neighbors=neighbors[i].one_bond_neighbors.number_elements();

      // initialize the two_bond_neighbors for each atom of the molecule
      neighbors[i].two_bond_neighbors.resize(0);
      neighbors[i].two_bond_neighbors_predecessor.resize(0);
      for (int j=0;j<number_of_one_bond_neighbors;j++)
        {
          int number_of_neighbors = neighbors[neighbors[i].one_bond_neighbors.item(j)].one_bond_neighbors.number_elements();
          for (int k=0;k<number_of_neighbors;k++) 
            {
              // get rid of the original atom (ith atom cannot be a two_bond neighbor of ith atom)
              int temp = neighbors[neighbors[i].one_bond_neighbors.item(j)].one_bond_neighbors.item(k);
              if (i != temp)
                {
                  neighbors[i].two_bond_neighbors.add(temp);
                  neighbors[i].two_bond_neighbors_predecessor.add(neighbors[i].one_bond_neighbors.item(j));
                }
            }
        }


      // printout the two_bond neighbors for the molecule
#if (DEBUG_SWITCH)
      cerr<<"Inside find_neighbors_for_molecule() -- two_bond_neighbors"<<endl;
      cerr<<"all two-bond neighbors   "<<i<<endl;
      cerr<<"     Atom     number"<<endl;
      for (int kk=0;kk<neighbors[i].two_bond_neighbors.number_elements();kk++)
      {
        cerr<<"      "<<atomi(neighbors[i].two_bond_neighbors.item(kk))->atomic_symbol() <<"      "<<neighbors[i].two_bond_neighbors.item(kk)<<endl;
      }
      cerr<<endl<<endl<<endl;
#endif

    }


  // find three_bond neighbors for all atoms of the molecule
  for (atom_number_t i=0;i<n_atom;i++)
    {
      // The idea here is that the three_bond neighbors of an atom is the union of the
      // one_bond neighbors of its two_bond neighbors except going backwards.
      int number_of_two_bond_neighbors=neighbors[i].two_bond_neighbors.number_elements();
      
      // initialize the two_bond_neighbors for each atom of the molecule
      neighbors[i].three_bond_neighbors.resize(0);
      neighbors[i].three_bond_neighbors_predecessor.resize(0);

      for (int j=0;j<number_of_two_bond_neighbors;j++)
        {
          int number_of_neighbors = neighbors[neighbors[i].two_bond_neighbors.item(j)].one_bond_neighbors.number_elements();
          // temp1 keep track of each two_bond_neighbors' predecessor
          int temp1 = neighbors[i].two_bond_neighbors_predecessor.item(j);        
          for (int k=0;k<number_of_neighbors;k++)
            {
              // temp2 is one of the value of possible three_bond_neighbor
      int temp2 = neighbors[neighbors[i].two_bond_neighbors.item(j)].one_bond_neighbors.item(k);

              // when ith atom's possible three bond neighbor is not itself and is not its two bond neighbor's predecessor, it is a true three bond neighbor 
              if ((temp2 != i) && (temp2 != temp1)) 
                {
                  neighbors[i].three_bond_neighbors.add(temp2);
                  neighbors[i].three_bond_neighbors_predecessor.add(temp1);
                }
            }
        }

      // printout the three_bond neighbors for the molecule
#if (DEBUG_SWITCH)
      cerr<<"Inside find_neighbors_for_molecule() -- three_bond_neighbors"<<endl;
      cerr<<"all three-bond neighbors   "<<i<<endl;
      cerr<<"     Atom     number"<<endl;
      for (int kk=0;kk<neighbors[i].three_bond_neighbors.number_elements();kk++)
      {
        cerr<<"      "<<atomi(neighbors[i].three_bond_neighbors.item(kk))->atomic_symbol() <<"      "<<neighbors[i].three_bond_neighbors.item(kk)<<endl;
      }
      cerr<<endl<<endl<<endl;
#endif

    }
      
  return rc;
}

/*
 * This function caluclates all the partial charge of all atoms of molecules
 * resulting from the one-bond neighbours effect
 * The result is stored in charge[] and returned to user -- call by reference
 */

int 
Molecule::_calculate_one_bond_partial_charge(arrays_holder neighbors[], double charge[])
{
  int n_atom = natoms();

#if (DEBUG_SWITCH)
  cerr<<"inside calculate_one_bond_partial_charge"<<endl;
#endif

  for (int ith_atom=0;ith_atom<n_atom; ith_atom++)
    {
#if (DEBUG_SWITCH)
      cerr<<"For "<<ith_atom<<"th atom"<<endl;
#endif
      // int ith_atom_atomic_number= ith_atom_ptr->atomic_number();
      double ith_atom_electronegativity;
      int row_number_i = _row_number(ith_atom);

      // if row_number returned is invalid -- exit else get the Ei

      if (row_number_i == OUT_OF_RANGE) 
        return 0;
      else 
        ith_atom_electronegativity=electronegativity[row_number_i];

      charge[ith_atom]=0;

      int n_neighbor = neighbors[ith_atom].one_bond_neighbors.number_elements();

      for (int j=0;j<n_neighbor;j++)
        {
          int jth_atom = neighbors[ith_atom].one_bond_neighbors.item(j);

          double jth_atom_electronegativity;
          
          // if row_number returned is invalid -- exit else get the Ej
          int row_number_j=_row_number(jth_atom);
          if (row_number_j == OUT_OF_RANGE)
            {
              return 0;
            }
          else 
            {
              jth_atom_electronegativity=electronegativity[row_number_j];
            }
          
          // calculate the partial charge
          charge[ith_atom] += (jth_atom_electronegativity-ith_atom_electronegativity)/_a_const(ith_atom, jth_atom);
#if (DEBUG_SWITCH)
          cerr<<"After "<<j+1<<"th neighbors    charge="<<charge[ith_atom]<<endl;
#endif
        }
    }
  return 1;
}


/*
 * This function calculates all the partial charge resulted from two and three bond neighbors
 * It is calculated this way:
 * qi(beta) = constant1 * Polarizability_of_i
 * qi(gamma) = constant2 * Polarizability_of_i
 * qi(beta)+qi(gamma)=(constant1+constant2) * Polarizability-of-i
 * Polarizability_of_i depends on the total charge, do iteration with calculated value of constant
 * Charge is changed as follows:  in the case of both beta and gamma, equal amount (opposit sign) is transfered
 * to the corresponding one_bond_neighbors.  (need to trace the neighbors)  
 */
int
Molecule::_calculate_two_and_three_bond_partial_charge(arrays_holder neighbors[], double charge[])
{
  int n_atoms = natoms();
  
  double * current_tmp_charge  = new double [n_atoms]; std::unique_ptr<double[]> free_current_tmp_charge(current_tmp_charge);
  double * previous_tmp_charge = new double [n_atoms]; std::unique_ptr<double[]> free_previous_tmp_charge(previous_tmp_charge);
  
  return _calculate_two_and_three_bond_partial_charge(neighbors, charge, current_tmp_charge, previous_tmp_charge);
}

int 
Molecule::_calculate_two_and_three_bond_partial_charge(arrays_holder neighbors[], double charge[],
                                                         double current_tmp_charge[], double previous_tmp_charge[])
{
  int n_atoms = natoms();
  
  for (int i=0; i<n_atoms;i++) 
    current_tmp_charge [i] = 0.0;
  for (int i=0; i<n_atoms;i++) 
    previous_tmp_charge [i] = 0.0;
#if (DEBUG_SWITCH)
  cerr<<"inside calculate_two_and_three_bond_partial_charge"<<endl;
#endif
  
  int not_done = 0;
  for (int loop_number =0; loop_number<MAX_LOOP; loop_number++) 
    {
      // most of the result should converge within 20 loops.  Within 20 loops fast convergence is used
      // If after 20 loops, it is still not converged, this might due to the polarizability flip-flop
      // A slow convergence method is used -- Although it is slow, it will converge
      if (loop_number>20)
        for (int i=0; i<n_atoms; i++) 
          previous_tmp_charge[i] = (current_tmp_charge[i]+previous_tmp_charge[i])/2;
      else
        for (int i=0; i<n_atoms; i++) 
          previous_tmp_charge[i] = current_tmp_charge[i];

      // reset the current_charge
      for (int i=0; i<n_atoms; i++) 
        current_tmp_charge[i]=0.0;
      
      for (int i=0; i<n_atoms; i++)
        {
          
#if (DEBUG_SWITCH)
          cerr<<"For "<<i<<"th atom"<<endl;
#endif
          
          double ith_atom_std_polarizability;
          double ith_atom_std_state_charge;
          int row_number_i = _row_number(i);
          
          // if row_number returned is invalid -- exit else get the Ei
          if (row_number_i == OUT_OF_RANGE) 
            {
              return 0;
            }
          else 
            {
              ith_atom_std_polarizability=polarizability[row_number_i];
              ith_atom_std_state_charge = std_state_charge[row_number_i];
            }
          
          double ith_atom_polarizability = ith_atom_std_polarizability * (1+constant_c*(ith_atom_std_state_charge - (charge[i]+previous_tmp_charge[i])));
          
          if (ith_atom_polarizability>=0) 
            {
              //go through each two bond neighbors and do the calculation
              int n_neighbor1 = neighbors[i].two_bond_neighbors.number_elements();
              for (int j=0; j<n_neighbor1; j++)
                {
                  int jth_atom = neighbors[i].two_bond_neighbors.item(j);
                  double jth_atom_electronegativity;
                  
                  // if row_number returned is invalid -- exit else get the Ej
                  int row_number_j=_row_number(jth_atom);
                  if (row_number_j == OUT_OF_RANGE)
                    {
                      return 0;
                    }
                  else 
                    {
                      jth_atom_electronegativity=electronegativity[row_number_j];
                    }
                  double tmp_charge = (jth_atom_electronegativity-electronegativity[H_ATOM])*ith_atom_polarizability/constant_b;
                  
                  // move charge from current atom to its corresponding one_bond_neighbor (two_bond)
                  current_tmp_charge [i] += tmp_charge;
                  current_tmp_charge [neighbors[i].two_bond_neighbors_predecessor.item(j)] += -tmp_charge;
                }
              
              int n_neighbor2 = neighbors[i].three_bond_neighbors.number_elements();
              for (int j=0; j<n_neighbor2; j++)
                {
                  int jth_atom = neighbors[i].three_bond_neighbors.item(j);
                  double jth_atom_electronegativity;
                  
                  // if row_number returned is invalid -- exit else get the Ej
                  int row_number_j=_row_number(jth_atom);
                  if (row_number_j == OUT_OF_RANGE)
                    {
                      return 0;
                    }
                  else 
                    {
                      jth_atom_electronegativity=electronegativity[row_number_j];
                    }
                  double tmp_charge = (jth_atom_electronegativity-electronegativity[H_ATOM])*ith_atom_polarizability/constant_b/constant_d;
                  
                  // move charge from current atom to its corresponding one_bond_neighbor (three_bond)
                  current_tmp_charge[i] += tmp_charge;
                  current_tmp_charge [neighbors[i].three_bond_neighbors_predecessor.item(j)] += -tmp_charge;
                }
            }
#if (DEBUG_SWITCH)
          else cerr<<"probability of atom "<<i<<" <0\n";
#endif
        }
      // end of calculating constant

      not_done = 0;
      for (int i=0; i<n_atoms; i++) 
        {
          double charge_difference;
          if (current_tmp_charge[i]>previous_tmp_charge[i]) charge_difference=current_tmp_charge[i]-previous_tmp_charge[i];
          else charge_difference=previous_tmp_charge[i]-current_tmp_charge[i];
          if (charge_difference > MAX_CHARGE_DIFFERENCE) 
            {
              not_done=1;
              break;
            }
        }
#if (DEBUG_SWITCH)
      for (int i=0; i<n_atoms; i++) 
        {
          cerr<<"Atom "<<i<<"\t";
          cerr<<"current charge="<<current_tmp_charge[i]<<"\tprevious charge="<<previous_tmp_charge[i]<<"\n";
        }
#endif
      if (!not_done) break;
    }

  // if the for loop was ended due to the MAX_LOOP, not due to convergence do something
  if (not_done) 
    {
      for (int i=0; i<n_atoms; i++)
        current_tmp_charge[i]=(current_tmp_charge[i]+previous_tmp_charge[i])/2;
    }
  for (int i=0; i<n_atoms; i++) 
    {
      charge[i] +=current_tmp_charge[i];
    }
  
  return 1;
}

/*
 * This function find the conjugated_fragment for each atom
 * It is done recursively. For C1=C2-C3 C3 is recursively searched.
 * This function and the function below is used to calculate the pi effect
 */

int
Molecule::_find_conjugated_fragment_for_atom(int ith_atom_number, int fragment_array[], int &fragment_no, int is_double)
{
  if (fragment_array [ith_atom_number]) 
    return 0;

  const Atom * ith_atom = _things[ith_atom_number];

  // first atom of the fragment
  if (is_double ==2)
    {
      if (ith_atom->ncon() == ith_atom->nbonds()) return 0;
      else
        {
          fragment_array[ith_atom_number] = fragment_no;
          int tmp_connection = ith_atom->ncon();
          for (int j=0;j<tmp_connection;j++) 
            {
              int jth_atom_number = other(ith_atom_number, j);
              
              int is_this_bond_double = 0;
              if (ith_atom->btype_to_connection(j)>1) is_this_bond_double = 1;
              _find_conjugated_fragment_for_atom(jth_atom_number, fragment_array, fragment_no, is_this_bond_double);
            }
          return 1;
        }
    }
  // bond leading to this atom is single
  if (is_double == 0)
    {
      fragment_array[ith_atom_number] = fragment_no;
      if (ith_atom->ncon() !=ith_atom->nbonds())
        {
          int tmp_connection = ith_atom->ncon();
          for (int j=0;j<tmp_connection;j++) 
            {
              int jth_atom_number = other(ith_atom_number, j);
              int is_this_bond_double = 0;
              
              if (ith_atom->btype_to_connection(j)>1)
                is_this_bond_double = 1;
              _find_conjugated_fragment_for_atom(jth_atom_number, fragment_array, fragment_no, is_this_bond_double);
            }
        }
      else
        {
          int tmp_connection = ith_atom->ncon();
          for (int j=0; j<tmp_connection; j++)
            {
              int jth_atom_number = other (ith_atom_number, j);
              if (ncon(jth_atom_number) != nbonds(jth_atom_number))
                _find_conjugated_fragment_for_atom(jth_atom_number, fragment_array, fragment_no, 0);
            }
        }
    }
  if (is_double == 1)
    {
      fragment_array[ith_atom_number] = fragment_no;
      int tmp_connection = ith_atom->ncon();
      for (int j=0; j<tmp_connection; j++)
        {
          int is_this_bond_double = 0;
          
          int jth_atom_number = other(ith_atom_number, j);
          
          if (ith_atom->btype_to_connection(j)>1) is_this_bond_double = 1;
          _find_conjugated_fragment_for_atom (jth_atom_number, fragment_array, fragment_no, is_this_bond_double);
        }
    }
  return 0;
}


/*
 * This function find all of the conjugated system (or aromatic) for the molecule
 * This function and the above function is used to calculate the pi effect.
 */

int
Molecule :: _find_conjugated_fragment (int conjugated_fragment_array[])
{
  (void) compute_aromaticity_if_needed();
  int n_atoms=natoms();

  // If any conjugated fragments are found, label it with fragment number
  int fragment_i = 1;

  for (int i=0;i<n_atoms;i++) conjugated_fragment_array [i]=0;

  for (int i=0;i<n_atoms;i++) 
    if (_find_conjugated_fragment_for_atom(i, conjugated_fragment_array, fragment_i, 2)) fragment_i++;

  return fragment_i - 1;
}


/* 
 *This function calculates all the partial charge for all atoms of the molecule
 */

int
Molecule :: _calculate_abraham_partial_charge (double partial_charge [])
{
  int n_atom = natoms();

  arrays_holder *neighbors = new arrays_holder [n_atom];

  int rc = _calculate_abraham_partial_charge (partial_charge, neighbors);

  delete [] neighbors;

  return rc;
}

int
Molecule::_calculate_abraham_partial_charge (double partial_charge[], arrays_holder neighbors [])
{
  int rc;

  int n_atom = natoms();

  // first find all the neighbors for each atoms of the molecule, the results are in neighbors[n_atom]
  rc = _find_neighbors_for_molecule(neighbors);
  if (0==rc) return 0;
  
  // checking to see if neighbors is ok
#if (DEBUG_SWITCH) 
  for (int i=0;i<n_atom;i++)
    {
      // printout the three_bond neighbors for the molecule
      cerr<<"Inside calculate_partial_charge, just exited from find_neighbors_for_molecule()"<<endl;
      cerr<<"all three-bond neighbors   "<<i<<endl;
      cerr<<"     Atom     number"<<endl;
      for (int kk=0;kk<neighbors[i].three_bond_neighbors.number_elements();kk++)
        cerr<<"      "<<atomi(neighbors[i].three_bond_neighbors.item(kk))->atomic_symbol() <<"      "<<neighbors[i].three_bond_neighbors.item(kk)<<endl;
      cerr<<endl<<endl<<endl;
    }
#endif
 
  // calculate partial charge resulted from one-bond neighbors
  rc = _calculate_one_bond_partial_charge(neighbors, partial_charge);
  if (0==rc) return 0;

  // checking to see if one_bond_partial_charges are ok
#if (DEBUG_SWITCH) 
  double total_charge=0;
  cerr<<"Inside calculate_partial_charge, just exited from calculate_one_bond_partial_charge()"<<endl;

  for (int i=0;i<n_atom;i++)
    {
      total_charge+=partial_charge[i];
      // printout the partial charge resulted from one_bond neighbors of the molecule
      cerr<<partial_charge[i]<< "      ";
    }
  cerr<<endl<<"Total Charge is:  "<<total_charge<<endl;
#endif

  
  //calculate partial charge resulted from two-bond and three-bond neighbors
  rc = _calculate_two_and_three_bond_partial_charge(neighbors, partial_charge);
  if (0==rc) return 0;
  
  // checking to see if two_three_bond_partial_charges are ok
#if (DEBUG_SWITCH) 
  cerr<<"Inside calculate_partial_charge, just exited from calculate_two_and_three_bond_partial_charge()"<<endl;
  total_charge =0;
  for (int i=0;i<n_atom;i++)
    {
      // printout the final partial charge for each atom of the molecule
      total_charge+=partial_charge[i];
      cerr<<partial_charge[i]<< "      ";
    }
  cerr<<endl<<"Total Charge is:  "<<total_charge<<endl;
#endif

  //release memory

  for (int i=0;i<n_atom;i++)
    {
      neighbors[i].one_bond_neighbors.resize(0);
      neighbors[i].two_bond_neighbors.resize(0);
      neighbors[i].two_bond_neighbors_predecessor.resize(0);
      neighbors[i].three_bond_neighbors.resize(0);
      neighbors[i].three_bond_neighbors_predecessor.resize(0);
    }

  return 1;
}

int
Molecule :: compute_Abraham_partial_charges ()
{
  //  make_implicit_hydrogens_explicit();
  int n_atoms = natoms();
  double * partial_charge = new double[n_atoms];   // deleted several lines later  #####

  for (int i=0; i<n_atoms; i++) partial_charge[i]=0;
  int rc = _calculate_abraham_partial_charge (partial_charge);


  // assign _charge;
  for (int i=0; i<n_atoms; i++) set_charge(i, static_cast<charge_t> (partial_charge[i] ) );

  // assign _charge_type and array of _charge;
  _charges->set_type ("ABRAHAM");

  delete [] partial_charge;   // #####

  return rc;// successful
}
/*
 * _gasteiger_row_number () return the row number to access constants for gasteiger partial charge calculation
 */
/*
int 
Molecule :: _gasteiger_row_number (atom_number_t ith_atom)
{
  int return_value = _row_number (ith_atom);

  if ((return_value == SSP3_ATOM) || (return_value == SSP2_ATOM)) return GASTEIGER_SULFUR;
  if ((return_value == NSP3_p_ATOM) || (return_value == NSP3_t_ATOM))
    {
      if (ncon(ith_atom) == 4) return GASTEIGER_N4_ATOM;
      else return NSP3_t_ATOM;
    }
  return return_value;
}
*/

/*
 */
/*int
Molecule :: _gasteiger_partial_charge_procedure (double partial_charge[])
{
  int n_atoms = natoms();
  int * atom_type = new int [n_atoms];
  if (find_simplified_sybyl_atom_type_sybyl_style (atom_type))
    {
      _gasteiger_partial_charge_procedure (partial_charge, atom_type);
      delete [] atom_type;
    }
  else
    {
      delete [] atom_type;
      cerr<<"Unreconized Atom Type Encountered  -- "<<molecule_name()<<endl;
      return 0;
    }

#if (DEBUG_PRINT)
  for (int i=0; i<natoms(); i++)
    cerr<<"Atom "<<i<<"   Type ="<<atom_type[i]<<endl;
#endif

  return 1;
}
*/


/*
 * main procedure for the gasteiger charge calculation
 * in order for future expantion (gasteiger_huckle charge) the function take the initial value of
 * partial_chage array before passing on to the partial charge calculation
 */

int
Molecule :: _gasteiger_partial_charge_procedure (double partial_charge[], atom_type_t atom_type[])
{
  int n_atom = natoms();

  double * atom_elec = new double [n_atom];
  Set_of_Atoms * neighbors = new Set_of_Atoms [n_atom];

  //  resizable_array<int> * neighbors = new resizable_array<int> [n_atom];
  double * old_charge = new double [n_atom];
  double * const_a = new double [n_atom];
  double * const_b = new double [n_atom];
  double * const_c = new double [n_atom];
  double * cation_elec = new double [n_atom];

  int rc = _gasteiger_partial_charge_procedure  (partial_charge, atom_elec, neighbors, old_charge, const_a, const_b, const_c, cation_elec, atom_type);

  //  cerr<<"THERE "<<rc<<endl;
  //release memory
  delete [] atom_elec;
  delete [] neighbors;
  delete [] old_charge;
  delete [] const_a;
  delete [] const_b;
  delete [] const_c;
  delete [] cation_elec;

  return rc;
}

int
Molecule :: _gasteiger_partial_charge_procedure (double partial_charge[], double gasteiger_atom_electronegativity [],
                                                 //resizable_array<int> one_bond_neighbors[],
                                                 Set_of_Atoms one_bond_neighbors [],
                                                 double previous_partial_charge [],
                                                 double constant_a [],
                                                 double constant_b [],
                                                 double constant_c [],
                                                 double gasteiger_atom_cation_electronegativity [],
                                                 atom_type_t atom_type [])
{
  int rc;

  // assign constant for calculation of electronegativity

  int n_atom = natoms();

  for (int i=0; i<n_atom; i++)
    {
      int tmp= atom_type [i];
      if (tmp==OUT_OF_RANGE) return 0;  // if there is unrecognized atom, return
      constant_a[i]=gasteiger_constant_a[tmp];
      constant_b[i]=gasteiger_constant_b[tmp];
      constant_c[i]=gasteiger_constant_c[tmp];
      gasteiger_atom_cation_electronegativity [i] = gasteiger_element_cation_electronegativity[tmp];
    }

  // find all bonding partners  -- NECESSARY???

  for (int i=0; i<n_atom; i++)
    {
      one_bond_neighbors[i].resize(0);
      rc = connections (i, one_bond_neighbors[i]);
      //      cerr<<"i="<<i<<"   THERE1 "<<rc<<endl;
      if (0==rc) return 0;
    }

  // not_done is used to end the loop.
  int not_done=0;

  // damping factor (1/2)^n
  double factor = 1.0;
  // In both sybyl & concord, the loop_number is set to 6 
  for (int loop_number=0; loop_number < 6 /*MAX_LOOP*/; loop_number ++)
    {
      factor=factor * 0.5;

      // keep a copy of the previous partial charge

      for (int i=0; i<n_atom; i++)
        previous_partial_charge[i]=partial_charge[i];

#if (DEBUG_PRINT)
      cerr<<"loop number ="<<loop_number<<endl;
      for (int i=0; i<n_atom; i++)
        cerr<<"Atom "<<i<<"\tcharge ="<<partial_charge [i]<<endl;
      cerr<<endl<<endl;
#endif
      for (int i=0; i<n_atom; i++)
        {
          gasteiger_atom_electronegativity [i] = constant_a[i] + (constant_b[i] + constant_c[i] * partial_charge[i])* partial_charge[i];
          //      cout<<"atom i="<<i<<"  electronegativity ="<<gasteiger_atom_electronegativity [i]<<endl;
        }

      for (int i=0; i<n_atom; i++)
        {
          // do the calculation for each atom
          double cation_electronegativity_current;
          int number_of_neighbors = one_bond_neighbors[i].number_elements();

          for (int j=0; j<number_of_neighbors; j++)
            {
              if (gasteiger_atom_electronegativity[i]>gasteiger_atom_electronegativity[one_bond_neighbors[i].item(j)])
                cation_electronegativity_current = gasteiger_atom_cation_electronegativity [one_bond_neighbors[i].item(j)];
              else 
                cation_electronegativity_current = gasteiger_atom_cation_electronegativity [i];
              //              cout<<"Atom "<<i<<"  charge increase ="<<(gasteiger_atom_electronegativity[one_bond_neighbors[i].item(j)] - gasteiger_atom_electronegativity[i])/cation_electronegativity_current*factor<<endl;
              partial_charge[i] += (gasteiger_atom_electronegativity[one_bond_neighbors[i].item(j)] - gasteiger_atom_electronegativity[i])/cation_electronegativity_current*factor;
            }
        }

      not_done = 0;
      for (int i=0; i<n_atom; i++) 
        {
          double charge_difference;
          if (partial_charge[i]>previous_partial_charge[i]) charge_difference=partial_charge[i]-previous_partial_charge[i];
          else charge_difference=previous_partial_charge[i] - partial_charge[i];
          if (charge_difference > MAX_CHARGE_DIFFERENCE) 
            {
              not_done=1;
              break;
            }
        }

      if (!not_done) break;
    }
  /* no convergence is necessary here because in concord and sybyl, loop is set to 6
  // if after MAX_LOOP, the charge calculation has not converge
  if (not_done)
    {
      for (int i=0;i<n_atom; i++)
        partial_charge[i]=(partial_charge[i]+previous_partial_charge[i])/2;
      cerr <<_molecule_name<<"  --charge did not converged"<<endl;
    }
  */

  //release memory
  for (int i=0;i<n_atom;i++)
    one_bond_neighbors[i].resize(0);

  return 1;
}


// wrapper function 
int
Molecule::compute_Gasteiger_partial_charges ()
{
  /*  int rc;
  int n_atoms = natoms();
  int * atom_type = new int [n_atoms];

  if (find_simplified_sybyl_atom_type_sybyl_style (atom_type))
    {
      rc = compute_Gasteiger_partial_charges (atom_type);
      delete [] atom_type;
    }
  else
    {
      delete [] atom_type;
      cerr<<"Unreconized Atom Type Encountered  -- "<<molecule_name()<<endl;
      return 0;
    }
  */
  // new version

  int n_atoms = natoms();
  atom_type_t * sybyl_atom_type = new unsigned int [n_atoms];
  std::unique_ptr<atom_type_t[]> free_sybyl_atom_type(sybyl_atom_type);

  if (NULL == _atom_type)
    {
      allocate_atom_types ();
      assert ( NULL != _atom_type);
    }

  if ("SYBYL" == _atom_type->ztype())
    {
      for (int i=0; i<n_atoms; i++)
        sybyl_atom_type[i] = atom_type(i);
    }
  else 
  {
    if (0 == find_simplified_sybyl_atom_type_sybyl_style(sybyl_atom_type))
    {
      cerr<<"Unreconized Atom Type Encountered  -- "<<molecule_name()<<endl;
      return 0;
    }
    if (0 == _atom_type->ztype().length())
    {
      for (int i=0; i<n_atoms; i++)
      {
        set_atom_type (i, sybyl_atom_type[i]);
      }
      _atom_type->set_type ("SYBYL");
    }
  }
 
  int rc = compute_Gasteiger_partial_charges (sybyl_atom_type, 0);

#if (DEBUG_PRINT)
  for (int i=0; i<natoms(); i++)
    cerr<<"Atom "<<i<<"   Type ="<<atom_type[i]<<endl;
#endif
  return rc;
}

// function for jon to include the pi charge part of the gasteiger partial charge computation

int
Molecule:: compute_Gasteiger_partial_charges_with_pi_charges ()
{
  int n_atoms = natoms();
  atom_type_t * sybyl_atom_type = new unsigned int [n_atoms];
  std::unique_ptr<atom_type_t[]> free_sybyl_atom_type(sybyl_atom_type);

  if (NULL == _atom_type)
    {
      allocate_atom_types ();
      assert ( NULL != _atom_type);
    }

  if ("SYBYL" == _atom_type->ztype())
    {
      for (int i=0; i<n_atoms; i++)
        sybyl_atom_type [i] = atom_type (i);
    }
  else 
    {
      int rc = find_simplified_sybyl_atom_type_sybyl_style(sybyl_atom_type);
      if (0==rc)
        {
          cerr<<"Unreconized Atom Type Encountered  -- "<<molecule_name()<<endl;
          return 0;
        }
      if (!_atom_type->ztype().length())
        {
          for (int i=0; i<n_atoms; i++)
            {
              set_atom_type (i, sybyl_atom_type[i]);
            }
          _atom_type->set_type ("SYBYL");
        }
    }
 
  int rc = compute_Gasteiger_partial_charges (sybyl_atom_type, 1);

#if (DEBUG_PRINT)
  for (int i=0; i<natoms(); i++)
    cerr<<"Atom "<<i<<"   Type ="<<atom_type[i]<<endl;
#endif

  return rc;
}


void
Molecule :: _assign_initial_charge (double partial_charge[], atom_type_t atom_type[])
{
  int n_atoms = natoms();
  for (int i=0; i<n_atoms; i++)
    {
      // deal with the special case for carboxyl group here
      if (atom_type [i] != ATOM_TYPE_OCO2)
        partial_charge[i]=formal_charge (i);
      else
        partial_charge [i] = -0.5;
    }
}

int
Molecule::_pi_gasteiger_partial_charge_computation_procedure (double partial_charge [], atom_type_t atom_type[])
{
  return 1;
} 
int
Molecule:: compute_Gasteiger_partial_charges (atom_type_t atom_type [], int compute_pi_charge)
{
  int n_atoms = natoms();

  double * partial_charge = new double [n_atoms]; std::unique_ptr<double[]> free_partial_charge(partial_charge);

  _assign_initial_charge (partial_charge, atom_type);

  int rc = _gasteiger_partial_charge_procedure (partial_charge, atom_type);

  if (compute_pi_charge)
  {
    (void) _pi_gasteiger_partial_charge_computation_procedure (partial_charge, atom_type);
  }

  for (int i=0; i<n_atoms; i++)
  {
    set_charge(i, static_cast<charge_t> (partial_charge[i]) );
  }

  // assign _charge_type and array of _charge;
  _charges->set_type ("GASTEIGER");

  return rc;// successful
}


// for this procedure, the atom type setting is the same as concord 4.02
// there is another function that set the atom type as sybyl called find_simplified_sybyl_atom_type_sybyl_style
int Molecule :: find_simplified_sybyl_atom_type (atom_type_t atom_type [])
{
  int rc = 1;
  int n_atoms = natoms();
  for (int i=0; i<n_atoms; i++)
    {
      const Atom * ith_atom = atomi(i);
      //      cout<<"Atom Number"<<ith_atom->atomic_number()<<"  Atom i="<<i<<" Connection="<<ncon(i)<<" Bond="<<nbonds(i)<<endl;
      switch (ith_atom->atomic_number()) {
      case 1: 
        atom_type [i] = ATOM_TYPE_H;
        break;
      case 6:
        // deal with descrepency between Ian's Pealman definition and Concord's definition
        if (is_aromatic(i) && in_ring_of_given_size (i, 6)) atom_type [i] = ATOM_TYPE_CAR; 
        else if ((nbonds(i) - ncon(i)) == 0) atom_type [i] = ATOM_TYPE_C3;
        else if ((nbonds(i) - ncon(i)) == 1) atom_type [i] = ATOM_TYPE_C2;
        else if ((nbonds(i) - ncon(i)) == 2) atom_type [i] = ATOM_TYPE_C1; 
        else atom_type [i] = ATOM_TYPE_C3;
        break;
      case 7:
        if (is_aromatic(i) && in_ring_of_given_size (i, 6)) atom_type [i] = ATOM_TYPE_NAR;
        else if (ncon(i) == 4) atom_type [i] = ATOM_TYPE_N4;
        else if ((ncon(i) == 3) && (nbonds(i) == 5)) atom_type [i] = ATOM_TYPE_NPL3;
        else if ((ncon(i) == 2) && (nbonds(i) == 3)) atom_type [i] = ATOM_TYPE_N2;
        else if ((ncon(i) == 3) && (nbonds(i) == 4) && (atomi(i)->formal_charge() == 1)) atom_type [i] = ATOM_TYPE_N2;
        else if ((ncon(i) == 1) && (nbonds(i) == 3)) atom_type [i] = ATOM_TYPE_N1;
        else if (is_aromatic(i)) atom_type [i] = ATOM_TYPE_NPL3;
        else 
          {
            atom_type [i] = ATOM_TYPE_N3;
            if (ncon(i)==3)
              {
                int car_count =0;
                int hetero_atom_count =0;
                int c2_count = 0;
                if (is_ring_atom (i))
                  {
                    for (int j=0; j<3; j++)
                      {
                        int atomj = other(i, j);
                        if ((nbonds(atomj) > ncon(atomj)) && (!is_aromatic(atomj))) 
                          {
                            // special case N-S(=O)(=O)-C, this S does not count as c2 count
                            int no_special_case =1;
                            if ((atomic_number(atomj) == 16) && is_ring_atom (i))
                              {
                                int O2_count = 0;
                                int tmp_connection = ncon(atomj);
                                for (int jj=0; jj<tmp_connection; jj++)
                                  {
                                    int jjth_atom = other (atomj, jj);
                                    if ((atomic_number (jjth_atom) == 8) && (ncon(jjth_atom) == 1) && (nbonds(jjth_atom) == 2))
                                      O2_count ++;
                                  }
                                if ((2==O2_count) && (0==hcount(atomj))) no_special_case = 0;
                              }
                            if (no_special_case) c2_count ++;
                          }
                        else if (is_aromatic(atomj) && ((!in_ring_of_given_size(i, 6)) || (!in_same_ring (i, atomj)))) car_count ++;
                        
                        if ((atomic_number(atomj) != 6) && (atomic_number(atomj) != 1))
                          hetero_atom_count ++;
                      }
                    if (in_ring_of_given_size(i, 5)) c2_count +=car_count;
                  }
                else
                  for (int j=0; j<3; j++)
                    {
                      int atomj = other(i, j);
                      if (nbonds(atomj) > ncon(atomj)) c2_count ++;
                      if (is_aromatic(atomj)) car_count ++;
                      
                      if ((atomic_number(atomj) != 6) && (atomic_number(atomj) != 1))
                        hetero_atom_count ++;
                    }
                //              cerr<<"HERE Atom ="<<i<<"  c2 count ="<<c2_count<<"   car count ="<<car_count<<"   heter count ="<<hetero_atom_count<<endl;
                if ((c2_count>1) || (car_count>1) || ((1==car_count) && (0==hetero_atom_count))) atom_type [i] = ATOM_TYPE_NPL3;
              }
          }
        {
          // special case N-C=S  (NPL3)
          int tmp_connection = ncon(i);
          for (int j=0; j<tmp_connection; j++)
            {
              int atomj = other (i, j);
              if ((atomi(atomj)->atomic_number() !=6) || (ncon(atomj) != 3))  continue;
              else
                for (int k=0; k<3; k++)
                  if ((atomi(other (atomj, k))->atomic_number() == 16) && (btype_to_connection (atomj, k) == 2))
                    {
                      atom_type [i] = ATOM_TYPE_NPL3;
                      break;
                    }
            }
          
          // special case N-C=O  (Nam)
          for (int j=0; j<tmp_connection; j++)
            {
              int atomj = other (i, j);
              if ((atomi(atomj)->atomic_number() !=6) || (ncon(atomj) != 3))  continue;
              else
                for (int k=0; k<3; k++)
                  if ((atomi(other (atomj, k))->atomic_number() == 8) && (btype_to_connection (atomj, k) == 2))
                    {
                      atom_type [i] = ATOM_TYPE_NAM;
                      break;
                    }
            }
        }
        break;
      case 8:
        if ((ncon(i) == 1) && (nbonds(i) == 2)) atom_type [i] = ATOM_TYPE_O2;
        else atom_type [i] = ATOM_TYPE_O3;
        break;
      case 9:
        atom_type [i] = ATOM_TYPE_F;
        break;
      case 15:
        atom_type [i] = ATOM_TYPE_P3;
        break;
      case 16:
        if (ncon(i) == nbonds(i)) atom_type [i] = ATOM_TYPE_S3;
        else 
          {
            atom_type [i] = ATOM_TYPE_S2;
            int O2_count = 0;
            int tmp_connection = ncon(i);
            for (int j=0; j<tmp_connection; j++)
              {
                int jth_atom = other (i, j);
                if ((atomic_number (jth_atom) == 8) && (ncon(jth_atom) == 1) && (nbonds(jth_atom) == 2))
                  O2_count ++;
              }
            if (1==O2_count) atom_type [i] = ATOM_TYPE_SO;
            if (2==O2_count) atom_type [i] = ATOM_TYPE_SO2;
          }
        break;
      case 17:
        atom_type [i] = ATOM_TYPE_CL;
        break;
      case 35:
        atom_type [i] = ATOM_TYPE_BR;
        break;
      case 53:
        atom_type [i] = ATOM_TYPE_I;
        break;
      default:
        atom_type [i] = OUT_OF_RANGE;
        rc = 0;
      }
    }

  return rc;
}

/*
 *  Find all doubly bonded oxygens bound to the central atom
 */

static int
count_doubly_bonded_oxygens_attached (const Molecule & m,
                                      atom_number_t centre,
                                      atom_number_t & o1,
                                      atom_number_t & o2)
{
  int doubly_bonded_oxygens_found = 0;
  
  const Atom * acs = m.atomi (centre);
  
  for (int i = 0; i < acs->ncon (); i++)
    {
      const Bond * b = acs->item (i);
      if (! b->is_double_bond ())
        continue;
      
      atom_number_t j = b->other (centre);
      
      if (8 != m.atomic_number (j))
        continue;
      
      if (0 == doubly_bonded_oxygens_found)
        o1 = j;
      else if (1 == doubly_bonded_oxygens_found)
        o2 = j;
      else
        {
          //      cerr << "is_charged_acid: too many doubly bonded oxygens, ignoring atom " << centre << endl;
          return 1;
        }
      
      doubly_bonded_oxygens_found++;
    }
  
  //  cerr << "Atom " << centre << " has " << doubly_bonded_oxygens_found << " doubly bonded oxygens\n";
  
  return doubly_bonded_oxygens_found;
}

/*
 *  Atom O1 is a charged oxygen atom. Is it part of a carboxyllic or sulf* acid
 * If a carboxyllic acid, we set O2. If Sulphonic, we also set O3
 */

static int
is_charged_acid (const Molecule & m,
                 atom_number_t o1,
                 atom_type_t atom_type [])
{
  const Atom * a1 = m.atomi (o1);
  
  assert (8 == a1->atomic_number ());
  assert (-1 == a1->formal_charge ());
  
  atom_number_t centre = a1->other (o1, 0);    // 1 negatively charged oxygen has only one connection
  
  const Atom * acs = m.atomi (centre);
  
  if (acs->ncon () < 3)
    return 0;
  
  if (acs->ncon () == acs->nbonds ())    // central must be unsaturated
    return 0;
  
  if (6 == acs->atomic_number ())     // possible carboxyllic acid
    {
      if (3 != acs->ncon ())
        return 0;
    }
  else if (16 == acs->atomic_number ())    // possible sulph* acid
    {
      if (4 != acs->ncon ())
        return 0;
    }
  else
    return 0;
  
  atom_number_t o2, o3;
  int no = count_doubly_bonded_oxygens_attached (m, centre, o2, o3);
  
  //  cerr << "Atom " << centre << " has " << no << " doubly bonded oxygens\n";
  if (0 == no)
    return 0;
  
  atom_type[o1] = ATOM_TYPE_OCO2;
  
  // Carboxyllic or sulph* type acid?
  
  if (16 == acs->atomic_number ())
    atom_type[centre] = ATOM_TYPE_SO2;
  else
    atom_type[centre] = ATOM_TYPE_C2;
  
  atom_type[o2] = ATOM_TYPE_OCO2;
  if (no > 1)
    atom_type[o3] = ATOM_TYPE_OCO2;

  return 1;
}


/*
 * this is the function checking if this aromatic ring has substitute in its neighbors
 */

int number_neighbor_has_branched_aromatic_atoms (const Molecule &m, atom_number_t this_atom, atom_number_t N_atom)
{
  int number = 0;

  int n_con = m.ncon(this_atom);

  for (int i=0; i<n_con; i++)
    {
      atom_number_t atomi = m.other(this_atom, i);
      int n_not_h_atom = 0;

      if ((N_atom != atomi) && (1 != m.atomic_number(atomi)))
        {
          int n_con_atomi = m.ncon(atomi);
          for (int j=0; j<n_con_atomi; j++)
            if (1 == m.atomic_number(m.other(atomi, j)))
              continue;
            else if ((8 == m.atomic_number(m.other(atomi, j))) && (m.btype_to_connection (atomi, j) ==2))
              return 0;
            else n_not_h_atom ++;
        }
      if (n_not_h_atom > 2) number++;
    }
  return number;
}

void Molecule :: _find_number_of_unsaturated_carbon_of_nitrogen_neighbor (atom_number_t i, int & c2_count, int & car_count, 
                                                                atom_number_t & one_car_atom)
{
  c2_count = 0;
  car_count =0;
  one_car_atom = 0;

  int nconi = ncon (i);
  if (is_ring_atom (i))
    {
      for (int j=0; j<nconi; j++)
        {
          int atomj = other(i, j);
          //cerr<<"in ring 5 ="<<in_ring_of_given_size(atomj, 5)<<"    is part of fused ring ="<<is_part_of_fused_ring_system(atomj)<<endl;
          if (is_aromatic(atomj))
            {
              car_count ++;
              one_car_atom = atomj;
            }
          else if ((nbonds(atomj) > ncon(atomj)) && in_ring_of_given_size (atomj, 5) && is_part_of_fused_ring_system (atomj) && (!in_same_ring(i, atomj)))
            {
              car_count ++;
              one_car_atom = atomj;
            }

          else if ((nbonds(atomj) > ncon(atomj)) && in_ring_of_given_size (atomj, 5) && !is_part_of_fused_ring_system (atomj) && (!in_same_ring(i, atomj)))
            {
              car_count ++;
              one_car_atom = atomj;
            }

          else if ((nbonds(atomj) > ncon(atomj)) && (!is_aromatic(atomj))) 
            {
              c2_count++;
            }
          
          else if ((atomic_number(atomj) == 6) && is_aromatic(atomj)) 
            {
              car_count++;
              one_car_atom = atomj;
            }
          else if (is_aromatic(atomj) && ((!in_ring_of_given_size(i, 6)) || (!in_same_ring (i, atomj)))) 
            {
              car_count ++;
              one_car_atom = atomj;
            }
          
          //                    if ((atomic_number(atomj) != 6) && (atomic_number(atomj) != 1))
          //                      hetero_atom_count ++;
        }
      //                    if (in_ring_of_given_size(i, 5)) c2_count +=car_count;
    }
  else
    for (int j=0; j<nconi; j++)
      {
        int atomj = other(i, j);
        // cerr<<"in ring 5 ="<<in_ring_of_given_size(atomj, 5)<<"    is part of fused ring ="<<is_part_of_fused_ring_system(atomj)<<endl;
        if ((is_aromatic(atomj)) || ( (nbonds(atomj) > ncon(atomj)) && in_ring_of_given_size (atomj, 5)/* && is_part_of_fused_ring_system (atomj)*/))
          {
            car_count ++;
            one_car_atom = atomj;
          }
        else if (nbonds(atomj) > ncon(atomj)) c2_count ++;
        
        //                    if ((atomic_number(atomj) != 6) && (atomic_number(atomj) != 1))
        //                      hetero_atom_count ++;
      }
  //cerr<<"HERE Atom ="<<i<<"  c2 count ="<<c2_count<<"   car count ="<<car_count<<endl;
  //            if ((c2_count>1) || (car_count>1) || ((1==car_count) && (0==hetero_atom_count))) atom_type [i] = ATOM_TYPE_NPL3;
}

// set the atom type for MCS computation, much simplier than sybyl style
int Molecule :: find_mcs_atom_type_similar_to_sybyl (int atom_type [])
{
  int n_atoms = natoms();

  set_vector ((int *) atom_type, n_atoms, UNDEFINED_TRIPOS_ATOM_TYPE);

  for (int i=0; i<n_atoms; i++)
    {
      if (UNDEFINED_TRIPOS_ATOM_TYPE != atom_type[i])  // atom already done as part of a group
        continue;
      
      Atom * ith_atom = _things[i];

      int nconi = ith_atom->ncon ();
      int nbondsi = ith_atom->nbonds();

      switch (ith_atom->atomic_number()) {
      case 1: 
        atom_type [i] = ATOM_TYPE_H;
        break;
      case 6:
        // deal with descrepency between Ian's Pealman definition and Concord's definition
        if (is_aromatic(i) && in_ring_of_given_size (i, 6)) atom_type [i] = ATOM_TYPE_CAR; 
        else if ((nbondsi - nconi) == 0) 
          {
            if (1 ==formal_charge(i)) atom_type [i] = ATOM_TYPE_CCAT;
            else if (-1==formal_charge(i)) atom_type [i] = ATOM_TYPE_C2; 
            else atom_type [i] = ATOM_TYPE_C3;
          }
        else if ((nbondsi - nconi) == 1) atom_type [i] = ATOM_TYPE_C2;
        else if ((nbondsi - nconi) == 2) atom_type [i] = ATOM_TYPE_C1; 
        else atom_type [i] = ATOM_TYPE_C3;
        break;
      case 7:
        //      cerr<<"atom i="<<i<<"\tncon="<<nconi<<"\tnbonds="<<nbondsi<<endl;
        if (is_aromatic(i) && in_ring_of_given_size (i, 6)) atom_type [i] = ATOM_TYPE_NAR;
        else if ((nconi == 4) && (nbondsi == 5)) atom_type [i] = ATOM_TYPE_N3;
        else if (nconi == 4) atom_type [i] = ATOM_TYPE_N4;
        else if ((nconi == 2) && (nbondsi == 4)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 3) && (nbondsi == 5)) atom_type [i] = ATOM_TYPE_NPL3;
        else if ((nconi == 2) && (nbondsi == 3)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 2) && (nbondsi == 2) && (atomi(i)->formal_charge() ==-1)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 3) && (nbondsi == 4) && (atomi(i)->formal_charge() == 1)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 1) && (nbondsi == 3)) atom_type [i] = ATOM_TYPE_N1;
        else if (is_aromatic(i)) atom_type [i] = ATOM_TYPE_NPL3;
        else 
          {
            atom_type [i] = ATOM_TYPE_N3;
          }

        break;
        
      case 8:
        if ((nconi == 1) && (nbondsi == 2))
          {
            atom_type [i] = ATOM_TYPE_O2;
            int atomj = other(i, 0);
            
            // Deal with Oco2 type
            // in case of C(=O)O-
            if ((atomic_number(atomj) == 6) && (ncon(atomj) >2))
              for (int k = 0; k<ncon(atomj); k++)
                {
                  int atomk = other (atomj, k);
                  if ((atomic_number(atomk) == 8) && (formal_charge(atomk) == -1))
                    {
                      atom_type [i] = ATOM_TYPE_OCO2;
                      break;
                    }
                }
            // in case of S(=O)(=O)O-
            if ((atomic_number(atomj) == 16) && (ncon(atomj) >3))
              {
                int O_minus_count =0;
                int double_bond_O_count =0;
                for (int k = 0; k<ncon(atomj); k++)
                  {
                    int atomk = other (atomj, k);
                    if ((atomic_number(atomk) == 8) && (formal_charge(atomk) == -1))
                      O_minus_count ++;
                    else if ((atomic_number(atomk) ==8) && (ncon(atomk) == 1) && (nbonds(atomk) == 2))
                      double_bond_O_count ++;
                  }
                if ((1==O_minus_count) && (2== double_bond_O_count))
                  atom_type [i] = ATOM_TYPE_OCO2;
              }
          }
        else if ((nconi == 1) && (formal_charge(i) == -1))
          {
            atom_type [i] = ATOM_TYPE_O3;
            int atomj = other(i, 0);

            // Deal with Oco2 type
            // in case of C(=O)O-
            if ((atomic_number(atomj) == 6) && (ncon(atomj) >1))
              for (int k = 0; k<ncon(atomj); k++)
                {
                  int atomk = other (atomj, k);
                  if ((atomic_number(atomk) == 8) && (ncon(atomk) == 1) && (nbonds(atomk) == 2))
                    {
                      atom_type [i] = ATOM_TYPE_OCO2;
                      break;
                    }
                }
            // in case of S(=O)(=O)O-
            if ((atomic_number(atomj) == 16) && ncon(atomj) >1)
              {
                int double_bond_O_count =0;
                for (int k=0; k<ncon(atomj); k++)
                  {
                    int atomk = other (atomj, k);
                    if ((atomic_number(atomk) == 8) && (ncon(atomk) == 1) && (nbonds(atomk) == 2))
                      double_bond_O_count ++;
                  }
                if (2==double_bond_O_count) atom_type [i] = ATOM_TYPE_OCO2;
              }
          }
           else atom_type [i] = ATOM_TYPE_O3;
        break;
      case 9:
        atom_type [i] = ATOM_TYPE_F;
        break;
      case 15:
        atom_type [i] = ATOM_TYPE_P3;
        break;
      case 16:
        if (nconi == nbondsi) atom_type [i] = ATOM_TYPE_S3;
        else 
          {
            atom_type [i] = ATOM_TYPE_S2;

            int O2_S2_count = 0;
            int tmp_connection = nconi;
            for (int j=0; j<tmp_connection; j++)
              {
                int jth_atom = other (i, j);
                if (((atomic_number (jth_atom) == 8) || (atomic_number (jth_atom) == 16)) && (ncon(jth_atom) == 1) && (nbonds(jth_atom) == 2))
                  O2_S2_count ++;
              }
            if (1==O2_S2_count) atom_type [i] = ATOM_TYPE_SO;
            if (2==O2_S2_count) 
              {
                atom_type [i] = ATOM_TYPE_SO2;
              }
          }
        break;
      case 17:
        atom_type [i] = ATOM_TYPE_CL;
        break;
      case 35:
        atom_type [i] = ATOM_TYPE_BR;
        break;
      case 53:
        atom_type [i] = ATOM_TYPE_I;
        break;
      default:
        atom_type [i] = OUT_OF_RANGE;
        return 0;
      }
    }
  
  return 1;
}

// set the atom type as defined by sybyl66
// there is another function that set the atom type as concord4.02 called find_simplified_sybyl_atom_type
int Molecule :: find_simplified_sybyl_atom_type_sybyl_style (int atom_type [])
{
  return find_simplified_sybyl_atom_type_sybyl_style ((atom_type_t *) atom_type);
  //  int n_atoms = natoms();
  //  atom_type_t * sybyl_atom_type = new atom_type_t [n_atoms];
  
}

// set the atom type as defined by sybyl66
// there is another function that set the atom type as concord4.02 called find_simplified_sybyl_atom_type
int Molecule :: find_simplified_sybyl_atom_type_sybyl_style (atom_type_t atom_type [])
{
  int n_atoms = natoms();

  set_vector ((int *) atom_type, n_atoms, UNDEFINED_TRIPOS_ATOM_TYPE);

  for (int i=0; i<n_atoms; i++)
    {
      if (UNDEFINED_TRIPOS_ATOM_TYPE != atom_type[i])  // atom already done as part of a group
        continue;
      
      Atom * ith_atom = _things[i];

      int nconi = ith_atom->ncon ();
      int nbondsi = ith_atom->nbonds();

      switch (ith_atom->atomic_number()) {
      case 1: 
        atom_type [i] = ATOM_TYPE_H;
        break;
      case 6:
        // deal with descrepency between Ian's Pealman definition and Concord's definition
        if (is_aromatic(i) && in_ring_of_given_size (i, 6)) atom_type [i] = ATOM_TYPE_CAR; 
        else if ((nbondsi - nconi) == 0) 
          {
            if (1 ==formal_charge(i)) atom_type [i] = ATOM_TYPE_CCAT;
            else if (-1==formal_charge(i)) atom_type [i] = ATOM_TYPE_C2; 
            else atom_type [i] = ATOM_TYPE_C3;
          }
        else if ((nbondsi - nconi) == 1) atom_type [i] = ATOM_TYPE_C2;
        else if ((nbondsi - nconi) == 2) atom_type [i] = ATOM_TYPE_C1; 
        else atom_type [i] = ATOM_TYPE_C3;
        break;
      case 7:
        //      cerr<<"atom i="<<i<<"\tncon="<<nconi<<"\tnbonds="<<nbondsi<<endl;
        if (is_aromatic(i) && in_ring_of_given_size (i, 6)) atom_type [i] = ATOM_TYPE_NAR;
        else if ((nconi == 4) && (nbondsi == 5)) atom_type [i] = ATOM_TYPE_N3;
        else if (nconi == 4) atom_type [i] = ATOM_TYPE_N4;
        else if ((nconi == 2) && (nbondsi == 4)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 3) && (nbondsi == 5)) atom_type [i] = ATOM_TYPE_NPL3;
        else if ((nconi == 2) && (nbondsi == 3)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 2) && (nbondsi == 2) && (atomi(i)->formal_charge() ==-1)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 3) && (nbondsi == 4) && (atomi(i)->formal_charge() == 1)) atom_type [i] = ATOM_TYPE_N2;
        else if ((nconi == 1) && (nbondsi == 3)) atom_type [i] = ATOM_TYPE_N1;
        else if (is_aromatic(i)) atom_type [i] = ATOM_TYPE_NPL3;
        else 
          {

            atom_type [i] = ATOM_TYPE_N3;
            if (nconi==3)
              {
                int car_count =0;
                atom_number_t one_car_atom = 0;
                int c2_count = 0;

                _find_number_of_unsaturated_carbon_of_nitrogen_neighbor (i, c2_count, car_count, one_car_atom);

                if (c2_count >1) atom_type [i] = ATOM_TYPE_NPL3;
                else if (c2_count >0)
                  {
                    atom_type [i] = ATOM_TYPE_NPL3;
                  }
                else if (c2_count + car_count >1)
                  {
                    if (!in_ring_of_given_size(i, 7)) atom_type [i] = ATOM_TYPE_NPL3;
                    else
                      {
                        const Ring * this_ring = ring_containing_atom (i);
                        if (this_ring->number_elements() !=7) atom_type [i] = ATOM_TYPE_NPL3;
                        else
                          {
                            int n_c_double_bond_in_the_ring =0;
                            int number_of_double_bond_in_the_ring =0;
                            for (int ii=0; ii<7; ii++)
                              {
                                int atom1 = this_ring->item(ii);
                                int atom2 = this_ring->item((ii+1)%7);
                                if (btype_between_atoms (atom1, atom2)>1)
                                  {
                                    number_of_double_bond_in_the_ring++;
                                    if (((atomic_number(atom1) == 7) && (atomic_number(atom2) == 6)) ||  ((atomic_number(atom1) == 6) && (atomic_number(atom2) == 7)))
                                      n_c_double_bond_in_the_ring ++;
                                  }
                              }
                            if ((number_of_double_bond_in_the_ring ==2) || (n_c_double_bond_in_the_ring))
                              atom_type [i] = ATOM_TYPE_NPL3;
                          }
                      }
                    
                  } 
                
                else if (car_count == 1)
                  {
                    if (hcount(i)>1) atom_type [i] = ATOM_TYPE_NPL3;
                    else if ((hcount(i) == 1) && (number_neighbor_has_branched_aromatic_atoms(*this, one_car_atom, i)<2))
                      atom_type [i] = ATOM_TYPE_NPL3; 

                    else if (is_ring_atom (i) && (!is_part_of_fused_ring_system (i)) && (!in_same_ring(i, one_car_atom) || in_ring_of_given_size(i, 5)))
                      {
                        // there is problem in aromaticity
                        /*                      if (!in_same_ring(i, one_car_atom) && in_ring_of_given_size(one_car_atom, 6))
                          {
                            const Ring * this_ring = ring_containing_atom(one_car_atom);
                            if (this_ring->number_elements() == 6)
                              {
                                int all_carbon_in_ring = 1;
                                int double_bond_in_ring = 0;
                                int all_atoms_unsaturated = 1;
                                for (int j=0; j<6; j++)
                                  {
                                    if ((btype_between_atoms(this_ring->item(j), this_ring->item((j+1) %6)) % 8) > 1)
                                      double_bond_in_ring++;
                                    if (6!=atomic_number(this_ring->item(j)))
                                      all_carbon_in_ring = 0;
                                    if (ncon(this_ring->item(j)) == nbonds(this_ring->item(j)))
                                      all_atoms_unsaturated = 0;
                                  }
                                if ((2==double_bond_in_ring) && (all_carbon_in_ring) && (all_atoms_unsaturated))
                                  atom_type [i] = ATOM_TYPE_NPL3;
                              }
                              }*/
                        if ((ATOM_TYPE_NPL3 != atom_type [i]) && !in_same_ring(i, one_car_atom) && in_ring_of_given_size(one_car_atom, 5))
                          {
                            const Ring * this_ring = ring_containing_atom(one_car_atom);
                            
                            int not_exception = 1;
                            if (this_ring->number_elements() ==5)
                              {
                                int all_carbon_in_ring = 1;
                                int double_bond_in_ring = 0;
                                int single_bonded_item = -1;
                                for (int j=0; j<5; j++)
                                  {
                                    // cerr<<"bond type="<<btype_between_atoms(this_ring->item(j), this_ring->item((j+1) %5))<<endl;
                                    if (((btype_between_atoms(this_ring->item(j), this_ring->item((j+1) %5)) %8) == 1) && ((btype_between_atoms(this_ring->item((j+2) % 5), this_ring->item((j+1) %5)) %8 )== 1))
                                      single_bonded_item = (j+1) % 5;
                                    if ((btype_between_atoms(this_ring->item(j), this_ring->item((j+1) %5)) % 8) > 1)
                                      double_bond_in_ring++;
                                    if (6!=atomic_number(this_ring->item(j)))
                                      all_carbon_in_ring = 0;
                                  }
                                
                                if ((2==double_bond_in_ring) && (-1 != single_bonded_item)) 
                                  {
                                    atom_number_t atomj = this_ring->item(single_bonded_item);
                                    if (nbonds(atomj) != ncon(atomj)) atom_type [i] = ATOM_TYPE_NPL3;
                                    else not_exception = 0;
                                  }
                                //cerr<<"atomi="<<i<<"\t double bond="<<double_bond_in_ring<<"\tsigle="<<single_bonded_item<<"\tnot exception ="<<not_exception<<endl;
                              }
                            

                            if (ATOM_TYPE_NPL3 != atom_type [i])
                              {

                                if (6==atomic_number(one_car_atom) && not_exception)
                                  {
                                    int conn = ncon(one_car_atom);
                                    
                                    for (int j=0; j<conn; j++)
                                      {
                                        if ((6==atomic_number(other(one_car_atom,j))) && in_same_ring(one_car_atom, other (one_car_atom, j)) && ((btype_to_connection (one_car_atom, j) % 8) >1))
                                          atom_type [i] = ATOM_TYPE_NPL3;
                                      }
                                  }
                                
                                
                                
                              }                     
                          }
                        
                        if (ATOM_TYPE_NPL3 != atom_type [i])
                          {
                            if (( 3 == ncon(one_car_atom)) && (4 == nbonds(one_car_atom)) && (6== atomic_number(one_car_atom)) && in_ring_of_given_size (one_car_atom, 5))
                              {
                                int is_in_fused_ring = is_part_of_fused_ring_system (one_car_atom);
                                
                                int one_car_not_in_exception_ring_for_c7_above_ring = 1;
                                int one_car_in_exception_ring_for_c5_c6_ring = 0;
                                
                                const Ring * this_ring = ring_containing_atom(one_car_atom);
                                if (this_ring->number_elements() ==5)
                                  {
                                    int n_of_o_or_s = 0;
                                    int n_of_n =0;
                                    int double_bond_in_ring = 0;
                                    for (int j=0; j<5; j++)
                                      {
                                        if (7 == atomic_number(this_ring->item(j))) n_of_n ++; 
                                        if ((8 == atomic_number(this_ring->item (j))) || (16 == atomic_number(this_ring->item (j))))
                                          n_of_o_or_s ++;

                                        if ((btype_between_atoms(this_ring->item(j), this_ring->item((j+1) %5)) %8) > 1)
                                          double_bond_in_ring++;
                                      }
                                //cerr<<"o or s="<<n_of_o_or_s<<"    no of n="<<n_of_n<<endl;
                                if ((1==n_of_o_or_s) && (n_of_n ==2)) one_car_not_in_exception_ring_for_c7_above_ring = 0;
                                if ((n_of_o_or_s<2) && (n_of_n == 1) && (!is_in_fused_ring) && (double_bond_in_ring<2)) 
                                  one_car_in_exception_ring_for_c5_c6_ring = 1;
                                  }
                                
                                //  cerr<<"HERE atom="<<i<<"   ring size="<<ring_containing_atom(i)->number_elements()<<endl;
                                if (((ring_containing_atom (i)->number_elements() > 6) /*&& (one_car_not_in_exception_ring_for_c7_above_ring)*/) || (((ring_containing_atom (i)->number_elements() ==5) || (ring_containing_atom (i)->number_elements() ==6)) && (one_car_in_exception_ring_for_c5_c6_ring)))
                                  {
                                    for (int j=0; j<3; j++)
                                      {
                                        atom_number_t atomj = other (one_car_atom, j);
                                        if ((7== atomic_number (atomj)) && ((btype_to_connection(one_car_atom, j) %8)>1))
                                          {
                                            atom_type [i] = ATOM_TYPE_NPL3;
                                            break;
                                          }
                                      }                         
                                  }
                              }
                          }
                        
                        if (ATOM_TYPE_NPL3 != atom_type [i])
                          {
                            int npl3_type = 0;
                            const Ring * this_ring = ring_containing_atom (i);
                            
                            for (int ii=0; ii<3; ii++)
                              {
                                atom_number_t atomkk = other (i, ii);
                                if (this_ring->contains(atomkk))
                                  {
                                    if (((atomic_number(atomkk) !=7) && (atomic_number(atomkk) !=8)) || (ncon(atomkk) != nbonds(atomkk))) continue;
                                    else
                                      for (int kk=0; kk<ncon(atomkk); kk++)
                                        //                                if (this_ring->contains(other(atomkk, kk)))
                                        {
                                          atom_number_t atomjj = other (atomkk, kk);
                                          
                                          if (atomjj!=i) 
                                            for (int iii=0; iii<ncon(atomjj); iii++)
                                              if ((!this_ring->contains(other(atomjj, iii))) && ((btype_to_connection (atomjj, iii) %8)>1))
                                                // if ((this_ring->contains(atomjj)) && (ncon(atomjj)<nbonds(atomjj)))
                                                npl3_type = 1;
                                        }
                                  }
                                /*          else 
                                            {
                                            if (atomic_number(atomkk) == 1) continue;  
                                            else 
                                            }*/
                              }
                            //}
                            // else if (is_ring_atom (i) && (!is_part_of_fused_ring_system (i)) && (!in_same_ring(i, one_car_atom) || in_ring_of_given_size(i, 5)))
                            
                            int number_of_atoms_in_ring = this_ring->number_elements();
                            int number_of_unsaturated_atom =0;
                            int n_in_ring =0;
                            for (int ii=0; ii<number_of_atoms_in_ring; ii++)
                              {
                                atom_number_t atomkk = this_ring->item(ii);
                                if (ncon(atomkk) < nbonds(atomkk)) number_of_unsaturated_atom ++;
                                if (atomic_number(atomkk) == 7) n_in_ring ++;
                              }
                            if ((number_of_atoms_in_ring == 6) && (number_of_unsaturated_atom >1) && (n_in_ring>1)) npl3_type = 1;
                            else if ((number_of_atoms_in_ring == 5) && (number_of_unsaturated_atom >1) && (n_in_ring>1)) npl3_type = 1;
                            
                            if (npl3_type) atom_type [i] = ATOM_TYPE_NPL3;
                            else atom_type [i] = ATOM_TYPE_N3;
                          }
                      }
                    
                    else if (is_part_of_fused_ring_system(i) && (in_same_ring(i, one_car_atom)) && in_ring_of_given_size(i, 5))
                      {                 
                        int npl3_type = 0;
                        const Ring * this_ring = ring_containing_atom (i);
                        int number_of_atoms_in_ring = this_ring->number_elements();
                        int number_of_unsaturated_atom =0;
                        int n_in_ring =0;
                        for (int ii=0; ii<number_of_atoms_in_ring; ii++)
                          {
                            atom_number_t atomkk = this_ring->item(ii);
                            if (ncon(atomkk) < nbonds(atomkk)) number_of_unsaturated_atom ++;
                            if (atomic_number(atomkk) == 7) n_in_ring ++;
                          }

                        if ((number_of_atoms_in_ring == 5) && (number_of_unsaturated_atom >2) && (n_in_ring>1)) npl3_type = 1;
                        
                        if (npl3_type) atom_type [i] = ATOM_TYPE_NPL3;
                        else atom_type [i] = ATOM_TYPE_N3;
                        
                      }
                    else if ( (is_part_of_fused_ring_system (i)) && (!in_same_ring(i, one_car_atom)))
                      {
                        const Ring * this_ring = ring_containing_atom (i);
                        int number_of_atoms_in_ring = this_ring->number_elements();
                        int number_of_unsaturated_atom =0;
                        for (int ii=0; ii<number_of_atoms_in_ring; ii++)
                          {
                            atom_number_t atomkk = this_ring->item(ii);
                            if (ncon(atomkk) < nbonds(atomkk)) number_of_unsaturated_atom ++;
                          }
                        //                      else if ((number_of_atoms_in_ring == 5) && (number_of_unsaturated_atom >0)) np3_type = 1
                        
                        int n_in_ring = 0;
                        
                        for (int ii=0; ii<3; ii++)
                          if (this_ring->contains (other (i, ii)) && (atomic_number(other(i, ii)) ==7))
                            n_in_ring ++;
                        
                        if ((number_of_atoms_in_ring == 6) && (number_of_unsaturated_atom >1) && (n_in_ring>0)) atom_type [i] = ATOM_TYPE_NPL3;
                        else atom_type [i] = ATOM_TYPE_N3;
                      }

                    else if (!is_ring_atom (i))
                      {
                        if (/*(is_part_of_fused_ring_system (one_car_atom)) && */(in_ring_of_given_size (one_car_atom, 5)))
                          {
                            int has_C_N_double_bond = 0;
                         
                            if (( 3 == ncon(one_car_atom)) && (4 == nbonds(one_car_atom)) && (6== atomic_number(one_car_atom)))
                              {
                                for (int j=0; j<3; j++)
                                  {
                                    atom_number_t atomj = other (one_car_atom, j);

                                    if ((7== atomic_number (atomj)) && ((btype_to_connection(one_car_atom, j) % 8)>1))
                                      {
                                        has_C_N_double_bond ++;
                                        //                                      atom_type [i] = ATOM_TYPE_NPL3;
                                        break;
                                      }
                                  }
                              }
                            
                            if ((is_part_of_fused_ring_system (one_car_atom)))// && (ATOM_TYPE_NPL3 != atom_type [i]))
                              {
                                const Ring * this_ring = ring_containing_atom (one_car_atom);
                            
                                int number_of_atoms_in_ring = this_ring->number_elements();
                                int number_of_unsaturated_atom =0;
                                int number_of_nitrogen_atom =0;
                                for (int ii=0; ii<number_of_atoms_in_ring; ii++)
                                  {
                                    atom_number_t atomkk = this_ring->item(ii);
                                    if (ncon(atomkk) < nbonds(atomkk)) number_of_unsaturated_atom ++;
                                    if (7== atomic_number(atomkk)) number_of_nitrogen_atom ++;
                                  }

                                if (5 == number_of_unsaturated_atom)
                                  atom_type [i] = ATOM_TYPE_NPL3;
                                if (((1==number_of_nitrogen_atom) || (3==number_of_nitrogen_atom)) && (has_C_N_double_bond))
                                  atom_type [i] = ATOM_TYPE_NPL3;
                              }
                          }
                      
                        if (number_neighbor_has_branched_aromatic_atoms (*this, one_car_atom, i)==0)
                          atom_type [i] = ATOM_TYPE_NPL3;
                      }
                    
                    else if (is_part_of_fused_ring_system (i))
                      {
                        if ((number_neighbor_has_branched_aromatic_atoms (*this, one_car_atom, i)==0) && (number_neighbor_has_branched_aromatic_atoms (*this, i, one_car_atom)==0))
                          atom_type [i] = ATOM_TYPE_NPL3;
                        //else if (in_ring_of_given_size(one_car_atom))
                    // {

                    //}
                      }
                  }
              }
          }
        
        // take care of special cases for NPL3 and Nam types
        if ((atom_type [i] != ATOM_TYPE_N2) && (atom_type [i] != ATOM_TYPE_N4) && (atom_type [i] != ATOM_TYPE_N1) && (atom_type [i] != ATOM_TYPE_NAR) && (nbondsi <4)) 
          {
            // special case N-C=S  (NPL3)
            for (int j=0; j<nconi; j++)
              {
                atom_number_t atomj = other (i, j);
                if ((atomi(atomj)->atomic_number() !=6) || (ncon(atomj) != 3))  continue;
                else
                  for (int k=0; k<3; k++)
                    if ((atomi(other (atomj, k))->atomic_number() == 16) && (btype_to_connection (atomj, k) == 2))
                      {
                        atom_type [i] = ATOM_TYPE_NPL3;
                        break;
                      }
              }
            
            // special case N-C=O and N-C=S (Nam)
            for (int j=0; j<nconi; j++)
              {
                int atomj = other (i, j);
                if ((atomi(atomj)->atomic_number() !=6) || (ncon(atomj) != 3))  continue;
                else
                  for (int k=0; k<3; k++)
                    {
                      // N-C=O
                      atom_number_t notused1, notused2;
                      
                      if (1 == count_doubly_bonded_oxygens_attached (*this, atomj, notused1, notused2))
                        {
                          atom_type [i] = ATOM_TYPE_NAM;
                          break;
                        }
                      // N-C=S
                      if ((atomi(other (atomj, k))->atomic_number() == 16) && (btype_to_connection (atomj, k) == 2))
                        {
                          atom_type [i] = ATOM_TYPE_NAM;
                          break;
                        }

                    }
              }
            
          }
        break;
        
      case 8:
        if (-1 == ith_atom->formal_charge () && is_charged_acid (*this, i, atom_type))
          ;
        else if ((nconi == 1) && (nbondsi == 2))
          atom_type [i] = ATOM_TYPE_O2;
        else atom_type [i] = ATOM_TYPE_O3;
          break;

          /*

        if ((nconi == 1) && (nbondsi == 2))
          {
            atom_type [i] = ATOM_TYPE_O2;
            int atomj = other(i, 0);
            
            // Deal with Oco2 type
            // in case of C(=O)O-
            if ((atomic_number(atomj) == 6) && (ncon(atomj) >2))
              for (int k = 0; k<ncon(atomj); k++)
                {
                  int atomk = other (atomj, k);
                  if ((atomic_number(atomk) == 8) && (formal_charge(atomk) == -1))
                    {
                      atom_type [i] = ATOM_TYPE_OCO2;
                      break;
                    }
                }
            // in case of S(=O)(=O)O-
            if ((atomic_number(atomj) == 16) && (ncon(atomj) >3))
              {
                int O_minus_count =0;
                int double_bond_O_count =0;
                for (int k = 0; k<ncon(atomj); k++)
                  {
                    int atomk = other (atomj, k);
                    if ((atomic_number(atomk) == 8) && (formal_charge(atomk) == -1))
                      O_minus_count ++;
                    else if ((atomic_number(atomk) ==8) && (ncon(atomk) == 1) && (nbonds(atomk) == 2))
                      double_bond_O_count ++;
                  }
                if ((1==O_minus_count) && (2== double_bond_O_count))
                  atom_type [i] = ATOM_TYPE_OCO2;
              }
          }
        else if ((nconi == 1) && (formal_charge(i) == -1))
          {
            atom_type [i] = ATOM_TYPE_O3;
            int atomj = other(i, 0);

            // Deal with Oco2 type
            // in case of C(=O)O-
            if ((atomic_number(atomj) == 6) && (ncon(atomj) >1))
              for (int k = 0; k<ncon(atomj); k++)
                {
                  int atomk = other (atomj, k);
                  if ((atomic_number(atomk) == 8) && (ncon(atomk) == 1) && (nbonds(atomk) == 2))
                    {
                      atom_type [i] = ATOM_TYPE_OCO2;
                      break;
                    }
                }
            // in case of S(=O)(=O)O-
            if ((atomic_number(atomj) == 16) && ncon(atomj) >1)
              {
                int double_bond_O_count =0;
                for (int k=0; k<ncon(atomj); k++)
                  {
                    int atomk = other (atomj, k);
                    if ((atomic_number(atomk) == 8) && (ncon(atomk) == 1) && (nbonds(atomk) == 2))
                      double_bond_O_count ++;
                  }
                if (2==double_bond_O_count) atom_type [i] = ATOM_TYPE_OCO2;
              }
          }
          else atom_type [i] = ATOM_TYPE_O3;
          break;*/
          case 9:
        atom_type [i] = ATOM_TYPE_F;
        break;
      case 15:
        atom_type [i] = ATOM_TYPE_P3;
        break;
      case 16:
        if (nconi == nbondsi) atom_type [i] = ATOM_TYPE_S3;
        else 
          {
            atom_type [i] = ATOM_TYPE_S2;

            /* lifted from Ian's modification  -- need further modification
             
               int notused1, notused2;
               
            int O2_count = count_doubly_bonded_oxygens_attached (*this, i, notused1, notused2);

            if (1==O2_count) atom_type [i] = ATOM_TYPE_SO;
            if (2==O2_count) atom_type [i] = ATOM_TYPE_SO2;
          }
        break;
        
            */
            
            int O2_S2_count = 0;
            int tmp_connection = nconi;
            for (int j=0; j<tmp_connection; j++)
              {
                int jth_atom = other (i, j);
                if (((atomic_number (jth_atom) == 8) || (atomic_number (jth_atom) == 16)) && (ncon(jth_atom) == 1) && (nbonds(jth_atom) == 2))
                  O2_S2_count ++;
              }
            if (1==O2_S2_count) atom_type [i] = ATOM_TYPE_SO;
            if (2==O2_S2_count) 
              {
                atom_type [i] = ATOM_TYPE_SO2;
              }

            
          }
        break;
          case 17:
        atom_type [i] = ATOM_TYPE_CL;
        break;
      case 35:
        atom_type [i] = ATOM_TYPE_BR;
        break;
      case 53:
        atom_type [i] = ATOM_TYPE_I;
        break;
      default:
        atom_type [i] = OUT_OF_RANGE;
        return 0;
      }
    }

  // special cases: if N3 atom in the same ring with two car type atom and these two  
  // car atoms are in the same ring

  /*  Set_of_Atoms ring_n3_atoms;
  Set_of_Atoms car_atoms;
  for (int i=0; i<n_atoms; i++)
    if ((atom_type [i] == ATOM_TYPE_N3) && is_ring_atom (i)) ring_n3_atoms.add(i);
    else if (atom_type [i] == ATOM_TYPE_CAR) car_atoms.add(i);
  
  int n_ring_n3_atoms = ring_n3_atoms.number_elements();
  int n_car_atoms = car_atoms.number_elements();
  if ((n_ring_n3_atoms > 0) && (n_car_atoms >1))
    {
      for (int i=0; i<n_ring_n3_atoms; i++)
        for (int j=0; j<n_car_atoms; j++)
          for (int k=j+1; k<n_car_atoms; k++)
            {
              int n3_atom = ring_n3_atoms.item(i);
              int car_1_atom = car_atoms.item(j);
              int car_2_atom = car_atoms.item(k);
              if (in_same_ring (n3_atom, car_1_atom) && in_same_ring (n3_atom, car_2_atom) && in_same_aromatic_ring (car_1_atom, car_2_atom))
                atom_type [n3_atom] = ATOM_TYPE_NPL3;
            }
            }*/



  for (int i=0; i<n_atoms; i++)
    {
      // all the neighbor nitrogen atom of Nam type, if car =1 then it is NPL3 type 
      if (ATOM_TYPE_NAM == atom_type [i])
        {
          int nconi = ncon (i);
          for (int j=0; j<nconi; j++)
            {
              atom_number_t atomj = other (i, j);
            if (ATOM_TYPE_N3 == atom_type [atomj])
              {
                int c2_count =0;
                int car_count = 0;
                atom_number_t one_car_atom =0;
                _find_number_of_unsaturated_carbon_of_nitrogen_neighbor (atomj, c2_count, car_count, one_car_atom);
                if (1 == car_count) atom_type [atomj] = ATOM_TYPE_NPL3;
              }
          }
        }

      // N in S(=O)(=O)[N-]S(=O)(=O) is N2 
      if ((ATOM_TYPE_N2 == atom_type [i]) && (-1 == formal_charge(i)) && (2==ncon(i)))
        {
          if ((ATOM_TYPE_SO2 == atom_type [other(i, 0)]) && (ATOM_TYPE_SO2 == atom_type [other(i, 1)]))
            atom_type [i] = ATOM_TYPE_N3;
        }


    }

  return 1;
}

void Molecule ::
_setup_conversion_array_between_huckel_matrix_and_atomi (resizable_array<int> * column_number_to_atomi, 
                                                         int atomi_to_column_number[])
{
  int n_atoms = natoms();

  column_number_to_atomi->resize(n_atoms);

  for (int i=0; i<n_atoms; i++)
    atomi_to_column_number [i] = -1;

  for (int i=0; i<n_atoms; i++)
    if (atomi(i)->atomic_number() != 1)
      {
        column_number_to_atomi->add(i);
        atomi_to_column_number [i] = column_number_to_atomi->number_elements()-1;
      }

#if (DEBUG_PRINT)
  cerr<<"Number of columns in matrix ="<<column_number_to_atomi->number_elements()<<endl;
  for (int i=0; i<column_number_to_atomi->number_elements(); i++)
    cerr<<"Column number "<<i<<"\tAtom number "<<column_number_to_atomi->item(i)<<endl;

  cerr<<endl;

  cerr<<"Number of atoms ="<<n_atoms<<endl;
  for (int i=0; i<n_atoms; i++)
    if (atomi_to_column_number[i] != -1)
      cerr<<"Atom number "<<i<<"\tColumn number "<<atomi_to_column_number[i]<<endl;

  cerr<<endl;
#endif
}

int
Molecule :: _setup_huckel_element_per_bond (const Bond * bondi, resizable_array<int> * column_to_atomi, 
                                            int column_number[], atom_type_t atom_type[], double ** huckel_matrix)
{
  int ith_atom = bondi->a1();
  int jth_atom = bondi->a2();

  // for atom type, atomi is always smaller than atomj
  if (atom_type [jth_atom]<atom_type [ith_atom])
    {
      int atom_tmp = ith_atom;
      ith_atom = jth_atom;
      jth_atom = atom_tmp;
    }

  int ci = column_number [ith_atom];
  int cj = column_number [jth_atom];

  switch (atom_type [ith_atom]) {
  case ATOM_TYPE_H:
    break;

    // if one atom is C3
  case ATOM_TYPE_C3:
    if ((atom_type [jth_atom] == ATOM_TYPE_C2) || (atom_type [jth_atom] == ATOM_TYPE_CAR)) 
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.7;
    break;

    // if one atom is C2
  case ATOM_TYPE_C2:
    switch (atom_type [jth_atom]) {

      // C2 C2 bond
    case ATOM_TYPE_C2: 
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.0;
      break;

      // C2-Car bond
    case ATOM_TYPE_CAR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.0;
      break;

      // C2 C1 bond
    case ATOM_TYPE_C1:
      // C2-C1
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.97;
      // C2=C1
      else if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.18;
      break;

      // C2-N4 bond
    case ATOM_TYPE_N4:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // C2-N3 bond
    case ATOM_TYPE_N3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // C2 N2 bond
    case ATOM_TYPE_N2:
      // C2-N2
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      // C2=N2
      else if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.10;
      break;

      // C2-N1 bond
    case ATOM_TYPE_N1:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // C2-O3 bond
    case ATOM_TYPE_O3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // C2-O2 bond
    case ATOM_TYPE_O2:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 2.00;
      break;

      // C2-S3 bond
    case ATOM_TYPE_S3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // C2-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // C2-Br
    case ATOM_TYPE_BR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // C2-Cl
    case ATOM_TYPE_CL:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // C2-F
    case ATOM_TYPE_F:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // C2-I
    case ATOM_TYPE_I:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.20;
      break;

      // C2=S
    case ATOM_TYPE_S2:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.20;
      break;

      // C2 NPL3 bond
    case ATOM_TYPE_NPL3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;

      // according to sybyl, the above is fine.  But according to concord, the parameter for NPL3-C2 bond
      // in N-(C=S)- is 0.00 
      if ((6 == atomic_number(ith_atom)) && (4 == nbonds(ith_atom)) && (3 == ncon(ith_atom)))
        {
          int tmp_connection = ncon(ith_atom);
          for (int kk=0; kk<tmp_connection; kk++)
            {
              int tmp_atom = other (ith_atom, kk);
              if (jth_atom == tmp_atom) continue;
              else if ((16 == atomic_number (tmp_atom)) && (1 == ncon(tmp_atom)) && (2 == nbonds(tmp_atom)))
                huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.00;
            }
        }
      break;
      
      // C2 Nam bond
    case ATOM_TYPE_NAM:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // C2-ar -Oco2
    case ATOM_TYPE_OCO2:
      if (bondi->is_aromatic())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 2.00;
      break;


    }
    break;

    // one atom is Car (aromatic)
  case ATOM_TYPE_CAR:
    switch (atom_type [jth_atom]) {

      // Car-Car bond
    case ATOM_TYPE_CAR:
      if (bondi->is_aromatic()) huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.0;
      else if (bondi->is_single_bond()) huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.9;
      break;

      // Car-C1 bond
    case ATOM_TYPE_C1:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // Car-N4 bond
    case ATOM_TYPE_N4:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // Car-N3 bond
    case ATOM_TYPE_N3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // Car-N2 bond
    case ATOM_TYPE_N2:

      //      if (bondi->is_aromatic())
      //        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.10;
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.10;
      break;

      // Car-N1 bond
    case ATOM_TYPE_N1:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // Car-O3 bond
    case ATOM_TYPE_O3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

    case ATOM_TYPE_O2:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 2.00;
      // Car S3
    case ATOM_TYPE_S3:
      // Car-S3 bond
      if (bondi->is_aromatic())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;

      else if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // Car Nar bond
    case ATOM_TYPE_NAR:
      if (bondi->is_aromatic())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      //      if (bondi->is_single_bond())
      //        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.9;
      //      if (bondi->is_double_bond())
      //        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.10;
      break;

      // Car-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // Car-Br
    case ATOM_TYPE_BR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // Car-Cl
    case ATOM_TYPE_CL:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // Car-F
    case ATOM_TYPE_F:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // Car-I
    case ATOM_TYPE_I:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.20;
      break;

      // Car NPL3 bond
    case ATOM_TYPE_NPL3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;
      
      // Car Nam bond
    case ATOM_TYPE_NAM:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

    }
    break;

  case ATOM_TYPE_C1:
    switch (atom_type [jth_atom]) {

      // C1-C1 bond
    case ATOM_TYPE_C1:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.04;
      else if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.24;
      else if (bondi->is_triple_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.40;
      break;

      // C1-N3 bond
    case ATOM_TYPE_N3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // C1-N2 bond
    case ATOM_TYPE_N2:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      else if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.20;
      break;

      // C1-N1 bond
    case ATOM_TYPE_N1:
      if (bondi->is_triple_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.50;
      break;

      // C1-O3 bond
    case ATOM_TYPE_O3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // C1=O2 bond
    case ATOM_TYPE_O2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.30;
      break;

      // C1 S3
    case ATOM_TYPE_S3:
      // Car-S3 bond
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // C1-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // C1-Br
    case ATOM_TYPE_BR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // C1-Cl
    case ATOM_TYPE_CL:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // C1-F
    case ATOM_TYPE_F:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // C1-I
    case ATOM_TYPE_I:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.30;
      break;

      // C1 NPL3 bond
    case ATOM_TYPE_NPL3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;
      
      // C1 Nam bond
    case ATOM_TYPE_NAM:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;
    }

    break;

    //N4-N2
  case ATOM_TYPE_N4:
    if ((atom_type [jth_atom] == ATOM_TYPE_N2) && (bondi->is_single_bond()))
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
    break;

    //N3-N2 & N3-N1
  case ATOM_TYPE_N3:
    // N3-N2
    if ((atom_type [jth_atom] == ATOM_TYPE_N2) && (bondi->is_single_bond()))
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
    // N3-N1
    else if ((atom_type [jth_atom] == ATOM_TYPE_N1) && (bondi->is_single_bond()))
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
    break;

  case ATOM_TYPE_N2:
    switch (atom_type [jth_atom]) {

      // N2-N2 bond
    case ATOM_TYPE_N2:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      else if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.10;
      break;

      // N2-N1 bond
    case ATOM_TYPE_N1:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.10;
      break;

      // N2-O3 bond
    case ATOM_TYPE_O3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // N2=O2 bond
    case ATOM_TYPE_O2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.20;
      break;

      // N2 S3
    case ATOM_TYPE_S3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.30;
      break;

      // N2 Nar
    case ATOM_TYPE_NAR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // N2-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // N2=S2 bond
    case ATOM_TYPE_S2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // N2 NPL3 bond
    case ATOM_TYPE_NPL3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;
      
      // N2 Nam bond
    case ATOM_TYPE_NAM:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;
    }
    break;

  case ATOM_TYPE_N1:
    switch (atom_type [jth_atom]) {

      // N1-N1 N1=N1 N1#N1
    case ATOM_TYPE_N1:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.20;
      break;

      // N1-O3 bond
    case ATOM_TYPE_O3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // N1=O2 bond
    case ATOM_TYPE_O2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.20;
      break;

      // N1-S3
    case ATOM_TYPE_S3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // N1=S2 bond
    case ATOM_TYPE_S2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

    }
    break;

  case ATOM_TYPE_NAR:
    switch (atom_type [jth_atom]) {
      // Nar-O3  Nar-ar-O3
    case ATOM_TYPE_O3:
      if ((bondi->is_single_bond()) || (bondi->is_aromatic()))
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // Nar-S3
    case ATOM_TYPE_S3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // Nar-Npl3
    case ATOM_TYPE_NPL3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // Nar-Nam
    case ATOM_TYPE_NAM:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // Nar-Nar
    case ATOM_TYPE_NAR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.0001;
      break;
    }
    break;

  case ATOM_TYPE_NPL3:

    switch (atom_type [jth_atom]) {

      // Npl3-O3
    case ATOM_TYPE_O3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // Npl3-S3
    case ATOM_TYPE_S3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // Npl3-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // Npl3-Ccat
    case ATOM_TYPE_CCAT:
      if (bondi->is_aromatic())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.9;
      break;
    }

    break;
  case ATOM_TYPE_NAM:
    switch (atom_type [jth_atom]) {

      // Nam-O3
    case ATOM_TYPE_O3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // Nam-S3
    case ATOM_TYPE_S3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // Nam-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;
    }

    break;
  case ATOM_TYPE_O3:
    switch (atom_type [jth_atom]) {

      // O3-O3 bond
    case ATOM_TYPE_O3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // O3-S3
    case ATOM_TYPE_S3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // O3-P3
    case ATOM_TYPE_P3:
      if (bondi->is_single_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // O3-Br
    case ATOM_TYPE_BR:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.40;
      break;

      // O3-Cl
    case ATOM_TYPE_CL:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.50;
      break;

      // O3-F
    case ATOM_TYPE_F:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.60;
      break;

      // O3-I
    case ATOM_TYPE_I:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.30;
      break;
    }
    break;

  case ATOM_TYPE_O2:
    switch (atom_type [jth_atom]) {

      // O2=S3
    case ATOM_TYPE_S3:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.30;
      break;

      // O2=P3
    case ATOM_TYPE_P3:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.40;
      break;

      // O2=Br
    case ATOM_TYPE_BR:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.90;
      break;

      // O2=Cl
    case ATOM_TYPE_CL:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // O2=F
    case ATOM_TYPE_F:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.70;
      break;

      // O3=I
    case ATOM_TYPE_I:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;

      // O2-S2
    case ATOM_TYPE_S2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;
    }
    break;

  case ATOM_TYPE_S3:
    switch (atom_type [jth_atom]) {

      // S3-S3
    case ATOM_TYPE_S3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.30;
      break;

      // S3-S2
    case ATOM_TYPE_S2:
      if (bondi->is_double_bond())
        huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
      break;

      // Nam-P3
    case ATOM_TYPE_P3:
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.00;
      break;
    }
    break;

    // S2-P3
  case ATOM_TYPE_S2:
    if (atom_type [jth_atom] == ATOM_TYPE_P3)
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 0.80;
    break;

    // P3-ar-Oco2
  case ATOM_TYPE_P3:
    if ((atom_type [jth_atom] == ATOM_TYPE_OCO2) && bondi->is_aromatic())
      huckel_matrix [ci][cj] = huckel_matrix [cj][ci] = 1.4;
    break;
  }
  
  return 1;
}

void
Molecule :: _setup_huckel_diagonal_element (int column_i, resizable_array<int> * column_number_to_atomi, atom_type_t atom_type[], double ** huckel_matrix, double number_of_electron [])
{
  int atom_i = column_number_to_atomi->item(column_i);
  number_of_electron [column_i] = electron_parameters [atom_type [atom_i]];

  huckel_matrix [column_i][column_i] = huckel_parameters [atom_type [atom_i]];

  if (atom_type [atom_i] == ATOM_TYPE_C2)
    {
      int number_of_connection = ncon(atom_i);
      double huckel_parameters = 0;
      for (int j=0; j<number_of_connection; j++)
        {
          int atom_j = other(atom_i, j);
          switch (atom_type [atom_j]) {
            // c2=N2
            //    case ATOM_TYPE_N2:
            //      if (bond_between_atoms (atom_i, atom_j)->is_double_bond())
            //        if (huckel_parameters<0.1) huckel_parameters =0.1;
            //      break;
          case ATOM_TYPE_BR:
            if (huckel_parameters<0.2) huckel_parameters =0.2;
            break;
          case ATOM_TYPE_CL:
            if (huckel_parameters<0.2) huckel_parameters =0.2;
            break;
          case ATOM_TYPE_F:
            if (huckel_parameters<0.1) huckel_parameters =0.1;
            break;
          case ATOM_TYPE_I:
            if (huckel_parameters<0.1) huckel_parameters =0.1;
            break;
          default:
            break;
          }
        }
      huckel_matrix [column_i][column_i] = huckel_parameters;
    }
  else if (atom_type [atom_i] == ATOM_TYPE_CAR)
    {
      int number_of_connection = ncon(atom_i);
      double huckel_parameters = 0;
      for (int j=0; j<number_of_connection; j++)
        {
          int atom_j = other(atom_i, j);
          switch (atom_type [atom_j]) {
            // CAR-ar-NAR
          case ATOM_TYPE_NAR:
            if (bond_between_atoms (atom_i, atom_j)->is_aromatic())
              if (huckel_parameters<0.1) huckel_parameters =0.1;
            break;
          case ATOM_TYPE_BR:
            if (huckel_parameters<0.2) huckel_parameters =0.2;
            break;
          case ATOM_TYPE_CL:
            if (huckel_parameters<0.2) huckel_parameters =0.2;
            break;
          case ATOM_TYPE_F:
            if (huckel_parameters<0.1) huckel_parameters =0.1;
            break;
          case ATOM_TYPE_I:
            if (huckel_parameters<0.1) huckel_parameters =0.1;
            break;
          default:
            break;
          }
        }
      huckel_matrix [column_i][column_i] = huckel_parameters;
    }
  else if (atom_type [atom_i] == ATOM_TYPE_N3)
    {
      int number_of_connection = ncon(atom_i);
      int number_of_hydrogen_connected =0;
      int number_of_carbon3_connected =0;
      for (int j=0; j<number_of_connection; j++)
        {
          int atom_j = other(atom_i, j);
          if (atom_type [atom_j] == ATOM_TYPE_H) number_of_hydrogen_connected ++;
          else if (atom_type [atom_j] == ATOM_TYPE_C3) number_of_hydrogen_connected ++; 
        }
      if (2 == number_of_hydrogen_connected) huckel_matrix [column_i][column_i] = 1.8;
      else if ((1 == number_of_hydrogen_connected) && (1 == number_of_carbon3_connected)) huckel_matrix [column_i][column_i] = 1.7;
      else if (2 == number_of_carbon3_connected) huckel_matrix [column_i][column_i] = 1.6;
    }
  else if (atom_type [atom_i] == ATOM_TYPE_C3)
    {
      int number_of_hydrogen =0;
      int tmp_connection = ncon(atom_i);
      for (int j=0; j<tmp_connection; j++)
        if (atom_type [other (atom_i, j)] == ATOM_TYPE_H)
          number_of_hydrogen ++;
      if (0 == number_of_hydrogen) 
        {
          huckel_matrix[column_i][column_i] = 0;
          number_of_electron [column_i] = 0;
        }
    }
}

int 
Molecule ::  _setup_molecule_huckel_matrix (resizable_array<int> * column_number_to_atomi, 
                                   int atomi_to_column_number[], atom_type_t atom_type[], 
                                   double ** huckel_matrix, double number_of_electron [])
{
  int matrix_size = column_number_to_atomi->number_elements();

  // initialize the matrix
  for (int i=0; i<matrix_size; i++)
    for (int j=0; j<matrix_size; j++)
      huckel_matrix [i][j] = 0;

  for (int i=0; i<matrix_size; i++)
    _setup_huckel_diagonal_element (i, column_number_to_atomi, atom_type, huckel_matrix, number_of_electron);

  int rc;

  int n_bonds = nedges();

  for (int i=0; i<n_bonds; i++)
    {
      const Bond * ith_bond = bondi(i);
      rc = _setup_huckel_element_per_bond (ith_bond, column_number_to_atomi, atomi_to_column_number, atom_type, huckel_matrix);
    }

  return 1;
}

int
Molecule :: _householder_tri_diagonalization (int matrix_size, double ** transformation_matrix, double ** huckel_matrix, 
                                              double temp_array [], double huckel_tmp_vector1[], double huckel_tmp_vector2[])
{
  for (int i=0; i<matrix_size-2; i++)
    {
      double sum =0;
      for (int j=i+2; j<matrix_size; j++)
        sum += huckel_matrix[i][j] * huckel_matrix[i][j];

      //      if (sum == 0)
      if (sum < 1e-14)
        {
          for (int j=i+2; j<matrix_size; j++)
            huckel_matrix [i][j] = huckel_matrix [j][i] = 0;

          double tmp = - huckel_matrix [i][i+1];

          huckel_matrix [i][i+1] = huckel_matrix[i+1][i] = tmp;
          for (int j=0; j<matrix_size; j++)
            transformation_matrix [i][j] = 0;

          transformation_matrix [i][i+1] = 2;

          transformation_matrix[i][matrix_size] = 1;
          transformation_matrix[i][matrix_size] = 2;
          continue;
        }

      if (fabs(huckel_matrix [i][i+1]) < 1e-5) huckel_matrix [i][i+1] = huckel_matrix [i+1][i] = 0;
      sum += huckel_matrix[i][i+1] * huckel_matrix[i][i+1];
      
      double sum_sqrt = sqrt (sum);

      // keep the sign the same as a [i][i+1]
      if (huckel_matrix[i][i+1]<0) sum_sqrt = -sum_sqrt;

      double const_2_k_squared = sum + huckel_matrix [i][i+1]* sum_sqrt;

      // setup the tranformation vector
      for (int k=0; k<i+1; k++)
        temp_array [k] = 0;
      temp_array [i+1] = huckel_matrix[i][i+1] + sum_sqrt;
      for (int k=i+2; k<matrix_size; k++)
        temp_array [k] = huckel_matrix [i][k];
      
      // for the ease of computation, compute 1 vector first
      for (int j=0; j<matrix_size; j++)
        {
          huckel_tmp_vector1 [j] = 0;
          for (int k=0; k<matrix_size; k++)
            huckel_tmp_vector1 [j] += huckel_matrix [j][k] * temp_array [k];
          huckel_tmp_vector1 [j] = huckel_tmp_vector1[j] / const_2_k_squared;
        }

      double temp2 = 0;
      for (int j=0; j<matrix_size; j++)
        temp2 += temp_array [j] * huckel_tmp_vector1 [j];
      temp2 = temp2 / const_2_k_squared;

      for (int j=0; j<matrix_size; j++)
        huckel_tmp_vector2 [j] = huckel_tmp_vector1 [j] - temp_array [j] * temp2 /2;

#if (DEBUG_PRINT)
      cerr<<"temp vector  (U)"<<endl;
      for (int iii=0; iii<matrix_size; iii++)
        cerr<<temp_array [iii]<<"  ";
      cerr<<endl;
      cerr<<"transformation vector 1"<<endl;
      for (int iii=0; iii<matrix_size; iii++)
        cerr<<huckel_tmp_vector1 [iii]<<"  ";
      cerr<<endl;
      cerr<<"transformation vector 2"<<endl;
      for (int iii=0; iii<matrix_size; iii++)
        cerr<<huckel_tmp_vector2 [iii]<<"  ";
      cerr<<endl;
#endif

      for (int j=0; j<matrix_size; j++)
        // compute only needed portion of the "new" huckel matrix
        if (j<i) continue;
        else if (j == i)
          {
            huckel_matrix[j][j+1] = huckel_matrix [j+1][j] = - sum_sqrt;
            for (int k = j+2; k<matrix_size; k++)
              huckel_matrix [j][k] = huckel_matrix [k][j] = 0;
          }
        else
          for (int k=j; k<matrix_size; k++)
            huckel_matrix [k][j] = huckel_matrix [j][k] = huckel_matrix [j][k] - temp_array [j] * huckel_tmp_vector2 [k] - temp_array [k] * huckel_tmp_vector2 [j];
      
      for (int j=0; j<matrix_size; j++)
        transformation_matrix[i][j] = temp_array[j];

      transformation_matrix[i][matrix_size] = const_2_k_squared;
    }
  return 1;
}
 
int 
Molecule :: _householder_tri_diagonalization (int matrix_size, double ** transformation_matrix, double ** huckel_matrix)
{
  double * temp_array = new double [matrix_size];
  double * tmp_vector1 = new double [matrix_size];
  double * tmp_vector2 = new double [matrix_size];

  int rc = _householder_tri_diagonalization (matrix_size, transformation_matrix, huckel_matrix, temp_array, tmp_vector1, tmp_vector2);

  delete [] temp_array;
  delete [] tmp_vector1;
  delete [] tmp_vector2;

  return rc;
}

double Molecule :: _infinite_norm (int matrix_size, double ** huckel_matrix)
{
  double norm = 0;
  for (int i=0; i<matrix_size; i++)
    {
      double column_value = 0; 
      for (int j=0; j<matrix_size; j++)
        column_value += fabs (huckel_matrix [i][j]);
      if (column_value > norm) norm = column_value;
    }
  return norm;
}

int
Molecule :: _number_of_sign_agreement (int size_of_matrix, double median, double sub_alpha[], double sub_beta_squared[], 
                                       double principal_minor_serie [])
{
  int number_of_sign_agreement = 0;
  
  principal_minor_serie [0] = 1;
  int positive_sign = 1;
  principal_minor_serie [1] = sub_alpha[0] - median;
  
  if (principal_minor_serie [1]>0) number_of_sign_agreement ++;
  else positive_sign = 0;
  
  for (int i=2; i<size_of_matrix+1; i++)
    {
      principal_minor_serie [i] = (sub_alpha[i-1]-median)*principal_minor_serie[i-1] - sub_beta_squared [i-1] * principal_minor_serie [i-2];
      if ((principal_minor_serie [i]>0 && positive_sign) || 
          (principal_minor_serie [i]<0 && !positive_sign)) number_of_sign_agreement++;
      else positive_sign = !positive_sign;
    }
  
  return number_of_sign_agreement;
}


int
Molecule :: _number_of_sign_agreement (int size_of_matrix, double median, double sub_alpha[], double sub_beta_squared[])
{
  double * principal_minor_serie = new double [size_of_matrix + 1];
  int number_of_sign_agreement = _number_of_sign_agreement (size_of_matrix, median, sub_alpha, sub_beta_squared, principal_minor_serie);
  delete [] principal_minor_serie;
  
  return number_of_sign_agreement;
}


void
Molecule :: _compute_sub_eigen_value (int number_of_eigen_value, int starting, 
                                      int size_of_matrix, int number_of_eigen_value_above_region,
                                      double sub_alpha[], double sub_beta_squared[], 
                                      double min_eigen_value, double max_eigen_value, 
                                      double sub_eigen_value[])
{
  if ((number_of_eigen_value == 1) && ((max_eigen_value - min_eigen_value)<1e-11*fabs (max_eigen_value)))
    {
      sub_eigen_value[starting] = (min_eigen_value + max_eigen_value) /2;
    }
  else if (max_eigen_value - min_eigen_value < number_of_eigen_value *1e-11)
    {
      for (int i=0; i<number_of_eigen_value; i++)
        sub_eigen_value[starting + i] = min_eigen_value + i * 1e-11;
      return;
    }
  else
    {
      double median = (min_eigen_value + max_eigen_value) /2;

      int number_of_sign_agreement = _number_of_sign_agreement (size_of_matrix, median, sub_alpha, sub_beta_squared);

      if ( 0 == number_of_sign_agreement - number_of_eigen_value_above_region)
        _compute_sub_eigen_value (number_of_eigen_value, starting, size_of_matrix, number_of_eigen_value_above_region, sub_alpha, sub_beta_squared, min_eigen_value, median, sub_eigen_value);
      else if (number_of_eigen_value == number_of_sign_agreement - number_of_eigen_value_above_region)
        _compute_sub_eigen_value (number_of_eigen_value, starting, size_of_matrix, number_of_eigen_value_above_region, sub_alpha, sub_beta_squared, median, max_eigen_value, sub_eigen_value);
      else
        {
          int new_number_of_eigen_value = number_of_eigen_value - number_of_sign_agreement + number_of_eigen_value_above_region;

          _compute_sub_eigen_value (new_number_of_eigen_value, starting, size_of_matrix, number_of_sign_agreement, sub_alpha, sub_beta_squared, min_eigen_value, median, sub_eigen_value);

          int new_starting = starting + new_number_of_eigen_value;
          new_number_of_eigen_value = number_of_sign_agreement - number_of_eigen_value_above_region;

          _compute_sub_eigen_value (new_number_of_eigen_value, new_starting, size_of_matrix, number_of_eigen_value_above_region, sub_alpha, sub_beta_squared, median, max_eigen_value, sub_eigen_value);
        }


    }
}

// sort the eigen vector according to the desending order of corresponding eigen value

void Molecule :: _sort_eigenvalue_eigenvector (int matrix_size, double eigenvalue[], 
                                               double ** eigen_vector, double tmp_array[])
{
  for (int i=0; i<matrix_size; i++)
    {
      double i_value = eigenvalue [i];
      double max_value = eigenvalue [i];
      for (int j=0; j<matrix_size; j++)
        tmp_array [j] = eigen_vector [j][i]; 
      
      int max_position = i;
      for (int j=i+1; j<matrix_size; j++)
        if (eigenvalue[j] > max_value)
          {
            max_value = eigenvalue [j];
            max_position = j;
          }
      
      eigenvalue[i] = max_value;
      eigenvalue[max_position] = i_value;
      
      for (int j=0; j<matrix_size; j++)
        {
          eigen_vector [j][i] = eigen_vector [j][max_position];
          eigen_vector [j][max_position] = tmp_array [j];
        }
    }
}

void Molecule :: _sort_eigenvalue_eigenvector (int matrix_size, double eigenvalue[], double ** eigen_vector)
{
  double * tmp_array = new double [matrix_size]; std::unique_ptr<double[]> free_tmp_array(tmp_array);

  _sort_eigenvalue_eigenvector (matrix_size, eigenvalue, eigen_vector, tmp_array);

  return;
}


void Molecule :: _triangulation (int matrix_size, double ** matrix)
{

  for (int i=1; i<matrix_size; i++)
    {
      // step 1:doing switch if needed  HERE, might need change in the future
      if (fabs(matrix [i][i-1])> fabs(matrix [i-1][i-1]))
        {
          double tmp = matrix [i][i-1];
          matrix[i][i-1] = matrix[i-1][i-1];
          matrix[i-1][i-1] = tmp;

          tmp = matrix [i][i];
          matrix[i][i] = matrix [i-1][i];
          matrix[i-1][i] = tmp;
          if (i<matrix_size-1)
            {
              tmp = matrix [i][i+1];
              matrix [i][i+1] = matrix [i-1][i+1];
              matrix [i-1][i+1] = tmp;
            }
        }
      

#if (DEBUG_PRINT)
      cerr<<"After Triangulation step 1 -- Matrix"<<endl;
      for (int ii=0; ii<matrix_size; ii++)
        {
          for (int jj=0; jj<matrix_size; jj++)
            cerr<<matrix[ii][jj]<<"\t";
          cerr<<endl;
        }
#endif

      // step 2:



      if (fabs(matrix[i][i-1])<1e-4)
        {
          matrix[i][i-1] = 0;
          continue;
        }


      double tmp_m = matrix[i][i-1] / matrix [i-1][i-1];
      matrix [i][i-1] = 0;
      // step 3:

      matrix[i][i] = matrix[i][i] - tmp_m * matrix[i-1][i];

      if (i<matrix_size-1)
        {
          matrix[i][i+1] = matrix[i][i+1] - tmp_m * matrix [i-1][i+1];
          //if (fabs(matrix[i][i+1]<1e-7)) matrix[i][i+1] = matrix[i][i+1]/fabs(matrix[i][i+1]) * 1e-7;
        }
#if (DEBUG_PRINT)
      cerr<<"After Triangulation step 2 & 3 -- Matrix"<<endl;
      for (int i=0; i<matrix_size; i++)
        {
          for (int j=0; j<matrix_size; j++)
            cerr<<matrix[i][j]<<"\t";
          cerr<<endl;
        }
#endif

    }
}

void Molecule :: _normalize_vector (int vector_size, double vector [])
{
  double squared_sum = 0;
  double sum_sqrt = 0;
  for (int i=0; i<vector_size; i++)
    squared_sum += vector [i]*vector[i];

  if (0 == squared_sum) return;
  else 
    {
      sum_sqrt = sqrt (squared_sum);
      for (int i=0; i<vector_size; i++)
        vector [i] = vector [i]/sum_sqrt;
      return;
    }
}

void Molecule :: _inverse_iteration (int matrix_size, double ** matrix, double eigen_vector [], double tmp_solution [])
{
  // initialize the vector
  for (int i=0; i<matrix_size; i++)
    tmp_solution [i] = 1.0;

  for (int ii=0; ii<3; ii++)
    {
      // solve the equation to get the eigen vector
      if (fabs(matrix [matrix_size-1][matrix_size-1]) == 0 ) tmp_solution [matrix_size-1] = tmp_solution [matrix_size-1]/1e-8;
      else tmp_solution [matrix_size-1] = tmp_solution [matrix_size-1]/fabs(matrix [matrix_size-1][matrix_size-1]);
      
      tmp_solution [matrix_size-2] = (tmp_solution [matrix_size-2] - matrix [matrix_size-2][matrix_size-1] * tmp_solution [matrix_size -1])/matrix [matrix_size-2][matrix_size-2];
      
      for (int i=matrix_size-3; i>-1; i--)
        tmp_solution [i] =  (tmp_solution [i] - tmp_solution [i+1]*matrix [i][i+1] - tmp_solution [i+2]*matrix[i][i+2]) /matrix [i][i];

      _normalize_vector (matrix_size, tmp_solution);      
    }

  for (int i=0; i<matrix_size; i++)
    if (fabs(tmp_solution [i]) > 1e-12)
      eigen_vector [i] = tmp_solution [i];
    else
      eigen_vector [i] = 0;
  return;
}

void
Molecule :: _inverse_iteration (int matrix_size, double ** matrix, double eigen_vector [])
{
  double * tmp_solution = new double [matrix_size]; std::unique_ptr<double[]> free_tmp_solution(tmp_solution);

  _inverse_iteration (matrix_size, matrix, eigen_vector, tmp_solution);

  return;
}

void
Molecule :: _compute_sub_eigenvector (int matrix_size, double eigen_value[], double alpha[], 
                                      double beta [], double **sub_eigen_vector, double ** matrix)
{

  for (int eigen_no = 0; eigen_no<matrix_size; eigen_no++)
    {
      for (int i=0; i<matrix_size; i++)
        for (int j=0; j<matrix_size; j++)
          matrix [i][j] = 0;
      
      for (int i=0; i<matrix_size; i++)
        {
          matrix [i][i] = alpha [i] - eigen_value [eigen_no];;
          if (i>0) matrix [i-1][i] = matrix [i][i-1] = beta [i];
        }
      
#if (DEBUG_PRINT)
      cerr<<"Before Triangulation -- Matrix"<<endl;
      for (int i=0; i<matrix_size; i++)
        {
          for (int j=0; j<matrix_size; j++)
            cerr<<matrix[i][j]<<"  ";
          cerr<<endl;
        }
#endif

      _triangulation (matrix_size, matrix);

#if (DEBUG_PRINT)
      cerr<<"After Triangulation -- Matrix"<<endl;
      for (int i=0; i<matrix_size; i++)
        {
          for (int j=0; j<matrix_size; j++)
            cerr<<matrix[i][j]<<"  ";
          cerr<<endl;
        }
      cerr<<endl<<endl;
#endif

      // each row stores 1 eigen vector for 1 eigen value
      _inverse_iteration (matrix_size, matrix, sub_eigen_vector [eigen_no]);
    }
}

void Molecule :: _compute_sub_eigenvector (int matrix_size, double eigen_value[], double alpha[], 
                                           double beta [], double **sub_eigen_vector)
{
  double ** matrix = new double * [matrix_size];
  for (int i=0; i<matrix_size; i++)
    matrix [i] = new double [matrix_size];

  _compute_sub_eigenvector (matrix_size, eigen_value, alpha, beta, sub_eigen_vector, matrix);

  // housekeeping
  for (int i=0; i<matrix_size; i++)
    delete [] matrix [i];
  delete [] matrix;
}

void Molecule :: _compute_eigenvector_procedure 
(int matrix_size, double ** huckel_matrix, double ** eigen_vector, double eigenvalue[], int starting_point, 
 int ending_point, double sub_alpha[], double sub_beta_squared[], double sub_eigen_value[], 
 double sub_beta[], double ** sub_eigen_vector)
{
  int number_of_eigen_value = ending_point - starting_point;
  double max_eigen_value = _infinite_norm (matrix_size, huckel_matrix);

  for (int i=0; i<number_of_eigen_value; i++)
    {
      sub_alpha [i] = huckel_matrix[starting_point+i][starting_point+i];
      if (i>0) sub_beta_squared[i] = huckel_matrix[starting_point+i-1][starting_point+i] * huckel_matrix[starting_point+i-1][starting_point+i];
    }
  
  double min_eigen_value = - max_eigen_value;
  _compute_sub_eigen_value (number_of_eigen_value, 0, number_of_eigen_value, 0, sub_alpha, sub_beta_squared, min_eigen_value, max_eigen_value, sub_eigen_value);
  
  for (int i=1; i<number_of_eigen_value; i++)
    sub_beta [i] = huckel_matrix[starting_point+i-1][starting_point+i];
  
  // for i th eigenvalue, eigen vector is stored in ith row
  _compute_sub_eigenvector (number_of_eigen_value, sub_eigen_value, sub_alpha, sub_beta, sub_eigen_vector); 
  
  for (int i=0; i<number_of_eigen_value; i++)
    {
      eigenvalue [starting_point +i] = sub_eigen_value [i];
      for (int j=0; j<number_of_eigen_value; j++)
        eigen_vector[starting_point+i][starting_point+j] = sub_eigen_vector [i][j];
    }
} 
  
void Molecule :: _transform_eigen_vector (int matrix_size, double ** eigen_vector, 
                                          double **transformation_matrix, double **tmp_eigen_vector)
{
  for (int i=0; i<matrix_size; i++)
    {
      for (int j=0; j<matrix_size; j++)
        tmp_eigen_vector [i][j] = eigen_vector [i][j];
      for (int j=0; j<matrix_size-2; j++)
        {
          double tmp = 0;
          for (int k=0; k<matrix_size; k++)
            tmp += transformation_matrix[matrix_size-3-j][k]*tmp_eigen_vector [i][k];
          for (int k=0; k<matrix_size; k++)
            tmp_eigen_vector [i][k] = tmp_eigen_vector [i][k] - tmp / transformation_matrix [matrix_size-3-j][matrix_size] * transformation_matrix[matrix_size-3-j][k];
        }
    }

  // after the exchange, eigen vector of ith eigen value is stored @ ith column     
  for (int i=0; i<matrix_size; i++)
    for (int j=0; j<matrix_size; j++)
      eigen_vector [i][j] = tmp_eigen_vector [j][i];
}

void Molecule :: _compute_eigenvector (int matrix_size, double ** huckel_matrix, 
                                       double ** eigen_vector, double eigenvalue[],
                                       double ** transformation_matrix)
{
  int starting_point;
  int ending_point;

  starting_point = 0;
  ending_point = 1;

  while (starting_point < matrix_size)
    {
      while ((ending_point < matrix_size) && (huckel_matrix[ending_point-1][ending_point] !=0))
        ending_point ++;

      int number_of_eigen_value = ending_point - starting_point;

      if (1 == number_of_eigen_value)
        {
          eigenvalue [starting_point] = huckel_matrix [starting_point][starting_point];
          eigen_vector [starting_point][starting_point] = 1.0;

          starting_point = ending_point;
          ending_point = starting_point +1;
        }
      else 
        {
          double * sub_alpha = new double [number_of_eigen_value];
          double * sub_beta_squared = new double [number_of_eigen_value];
          double * sub_eigen_value = new double [number_of_eigen_value];
          double * sub_beta = new double [number_of_eigen_value];
          double ** sub_eigen_vector = new double * [number_of_eigen_value];

          for (int i=0; i<number_of_eigen_value; i++)
            sub_eigen_vector [i] = new double [number_of_eigen_value];
          
          _compute_eigenvector_procedure (matrix_size, huckel_matrix, eigen_vector, eigenvalue, starting_point, ending_point, sub_alpha, sub_beta_squared, sub_eigen_value, sub_beta, sub_eigen_vector);
          
          delete [] sub_alpha;
          delete [] sub_beta_squared;
          delete [] sub_eigen_value;
          delete [] sub_beta;
          for (int i=0; i<number_of_eigen_value; i++)
            delete [] sub_eigen_vector [i];
          
          delete [] sub_eigen_vector;

          starting_point = ending_point;
          ending_point = starting_point +1;
        }
    }
  
  double ** tmp_eigen_vector = new double * [matrix_size];
  for (int i=0; i<matrix_size; i++)
    tmp_eigen_vector [i] = new double [matrix_size];

  _transform_eigen_vector (matrix_size, eigen_vector, transformation_matrix, tmp_eigen_vector);

  for (int i=0; i<matrix_size; i++)
    delete [] tmp_eigen_vector [i];
  delete [] tmp_eigen_vector;

  _sort_eigenvalue_eigenvector (matrix_size, eigenvalue, eigen_vector);
  
#if (DEBUG_PRINT)
  cerr<<"After sorting"<<endl;
  cerr<<"eigen value"<<endl;
  for (int i=0; i<matrix_size; i++)
    cerr<<eigenvalue[i]<<"   ";
  cerr<<endl<<"eigen_vector"<<endl;
  for (int i=0; i<matrix_size; i++)
    {
      for (int j=0; j<matrix_size; j++)
        cerr<<eigen_vector [i][j]<<"   ";
      cerr<<endl;
    }
#endif

  return;
}

int
Molecule :: _assign_huckel_charge_to_atom (resizable_array<int> * huckel_matrix_column_number_to_atomi, 
                                         double number_of_electron [],
                                         double **eigenvector,
                                         double partial_charge [],
                                         double electron_density [])
{
  int matrix_size = huckel_matrix_column_number_to_atomi->number_elements();

  double total_electrons = 0;
  for (int i=0; i<matrix_size; i++)
    total_electrons +=number_of_electron [i];

  for (int i=0; i<matrix_size; i++)
    {
      electron_density[i] = 0;
      double tmp_electrons = total_electrons;
      for (int j=0; j<matrix_size; j++)
        {
          if (tmp_electrons >= 2)
            {
              electron_density [i] += 2 * eigenvector [i][j] * eigenvector [i][j];
              tmp_electrons = tmp_electrons -2;
              if (tmp_electrons <=0) break;
            }
          else
            {
              electron_density[i] += tmp_electrons *eigenvector[i][j] * eigenvector [i][j];
              break;
            }
        }
      partial_charge[huckel_matrix_column_number_to_atomi->item(i)] = number_of_electron [i] - electron_density [i];
    }
  
  return 1;
}

int
Molecule::_assign_huckel_charge_to_atom (resizable_array<int> * huckel_matrix_column_number_to_atomi,
                        double number_of_electron [],
                        double **eigenvector,
                        double partial_charge [])
{
  int matrix_size = huckel_matrix_column_number_to_atomi->number_elements();

  double * electron_density = new double [matrix_size]; std::unique_ptr<double[]> free_electron_density(electron_density);

  return _assign_huckel_charge_to_atom(huckel_matrix_column_number_to_atomi, number_of_electron, eigenvector, partial_charge, electron_density);
}

int Molecule :: find_eigen_value_for_matrix (int matrix_size,
                                             double ** matrix,
                                             double eigenvalue [],
                                             double ** eigenvector_matrix)
{
  double ** transformation_matrix = new double * [matrix_size];
  for (int i=0; i<matrix_size; i++)
    transformation_matrix[i] = new double [matrix_size+1];

  int rc = _householder_tri_diagonalization (matrix_size, transformation_matrix, matrix);

  for (int i=0; i<matrix_size-1; i++)
    if (fabs(matrix[i][i+1]) < 1e-7)
      matrix[i+1][i] = matrix[i][i+1] = 0.0;

  // initialize the eigenvector and eigen value
  for (int i=0; i<matrix_size; i++)
    {
      eigenvalue [i] =0.0;
      for (int j=0; j<matrix_size; j++)
        eigenvector_matrix [i][j] = 0.0;
    }
  _compute_eigenvector (matrix_size, matrix, eigenvector_matrix, eigenvalue, transformation_matrix);

  for (int i=0; i<matrix_size; i++)
    delete [] transformation_matrix [i];
  delete [] transformation_matrix;

  return rc;
}






// new method

void
Molecule :: _compute_sub_eigen_value_with_low_accuracy (int number_of_eigen_value,
                        int starting, int size_of_matrix, int number_of_eigen_value_above_region,
                         double sub_alpha[], double sub_beta_squared[], 
                         double min_eigen_value, double max_eigen_value,
                         double sub_eigen_value[])
{
  if ((number_of_eigen_value == 1) && ((max_eigen_value - min_eigen_value)<1e-7*fabs (max_eigen_value)))
    {
      sub_eigen_value[starting] = (min_eigen_value + max_eigen_value) /2.0;
    }
  else if (max_eigen_value - min_eigen_value < number_of_eigen_value *1e-7)
    {
      for (int i=0; i<number_of_eigen_value; i++)
        sub_eigen_value[starting + i] = min_eigen_value + i * 1e-7;
      return;
    }
  else
    {
      double median = (min_eigen_value + max_eigen_value) /2.0;

      int number_of_sign_agreement = _number_of_sign_agreement (size_of_matrix, median, sub_alpha, sub_beta_squared);

      if ( 0 == number_of_sign_agreement - number_of_eigen_value_above_region)
        _compute_sub_eigen_value_with_low_accuracy (number_of_eigen_value, starting, size_of_matrix, number_of_eigen_value_above_region, sub_alpha, sub_beta_squared, min_eigen_value, median, sub_eigen_value);
      else if (number_of_eigen_value == number_of_sign_agreement - number_of_eigen_value_above_region)
        _compute_sub_eigen_value_with_low_accuracy (number_of_eigen_value, starting, size_of_matrix, number_of_eigen_value_above_region, sub_alpha, sub_beta_squared, median, max_eigen_value, sub_eigen_value);
      else
        {
          int new_number_of_eigen_value = number_of_eigen_value - number_of_sign_agreement + number_of_eigen_value_above_region;

          _compute_sub_eigen_value_with_low_accuracy (new_number_of_eigen_value, starting, size_of_matrix, number_of_sign_agreement, sub_alpha, sub_beta_squared, min_eigen_value, median, sub_eigen_value);

          int new_starting = starting + new_number_of_eigen_value;
          new_number_of_eigen_value = number_of_sign_agreement - number_of_eigen_value_above_region;

          _compute_sub_eigen_value_with_low_accuracy (new_number_of_eigen_value, new_starting, size_of_matrix, number_of_eigen_value_above_region, sub_alpha, sub_beta_squared, median, max_eigen_value, sub_eigen_value);
        }
    }
}


// sort the eigen according to the desending order of corresponding eigen value
void
Molecule :: _sort_eigenvalue_without_eigenvector (int matrix_size,
                                                  double eigenvalue[])
{
  for (int i=0; i<matrix_size; i++)
    {
      double i_value = eigenvalue [i];
      double max_value = eigenvalue [i];
      
      int max_position = i;
      for (int j=i+1; j<matrix_size; j++)
        if (eigenvalue[j] > max_value)
          {
            max_value = eigenvalue [j];
            max_position = j;
          }
      
      eigenvalue[i] = max_value;
      eigenvalue[max_position] = i_value;
    }
}

void
Molecule :: _compute_eigen_value_without_eigen_vector_procedure (int matrix_size,
                                        double ** huckel_matrix,
                                        double eigenvalue[], 
                                        int starting_point,
                                        int ending_point,
                                        double sub_alpha[],
                                        double sub_beta_squared[],
                                        double sub_eigen_value[])
{
  int number_of_eigen_value = ending_point - starting_point;
  double max_eigen_value = _infinite_norm (matrix_size, huckel_matrix);
  
  for (int i=0; i<number_of_eigen_value; i++)
    {
      sub_alpha[i] = huckel_matrix[starting_point+i][starting_point+i];
      if (i>0) sub_beta_squared[i] = huckel_matrix[starting_point+i-1][starting_point+i] * huckel_matrix[starting_point+i-1][starting_point+i];
    }
  
  double min_eigen_value = - max_eigen_value;

  _compute_sub_eigen_value_with_low_accuracy (number_of_eigen_value, 0, number_of_eigen_value, 0, sub_alpha, sub_beta_squared, min_eigen_value, max_eigen_value, sub_eigen_value);
  
  for (int i=0; i<number_of_eigen_value; i++)
    eigenvalue [starting_point +i] = sub_eigen_value [i];
} 


void Molecule :: _compute_eigen_value_without_eigen_vector (int matrix_size, double ** huckel_matrix, 
                                                            double eigenvalue[], double ** transformation_matrix)
{
  int starting_point;
  int ending_point;
  
  starting_point = 0;
  ending_point = 1;
  
  while (starting_point < matrix_size)
    {
      while ((huckel_matrix [ending_point-1][ending_point] !=0) && (ending_point<matrix_size)) ending_point ++;
      
      int number_of_eigen_value = ending_point - starting_point;
      
      if (1 == number_of_eigen_value)
        {
          eigenvalue [starting_point] = huckel_matrix [starting_point][starting_point];
          
          starting_point = ending_point;
          ending_point = starting_point +1;
        }
      else 
        {
          double * sub_alpha = new double [number_of_eigen_value];
          double * sub_beta_squared = new double [number_of_eigen_value];
          double * sub_eigen_value = new double [number_of_eigen_value];
          
          _compute_eigen_value_without_eigen_vector_procedure (matrix_size, huckel_matrix, eigenvalue, starting_point, ending_point, sub_alpha, sub_beta_squared, sub_eigen_value);
          
          delete [] sub_alpha;
          delete [] sub_beta_squared;
          delete [] sub_eigen_value;

          starting_point = ending_point;
          ending_point = starting_point +1;
        }
    }
  
  _sort_eigenvalue_without_eigenvector (matrix_size, eigenvalue);
  
#if (DEBUG_PRINT)
  cerr<<"After sorting"<<endl;
  cerr<<"eigen value"<<endl;
  for (int i=0; i<matrix_size; i++)
    cerr<<eigenvalue[i]<<"   ";
  cerr<<endl<<"eigen_vector"<<endl;
  for (int i=0; i<matrix_size; i++)
    {
      for (int j=0; j<matrix_size; j++)
        cerr<<eigen_vector [i][j]<<"   ";
      cerr<<endl;
    }
#endif

  return;
}


int Molecule :: find_eigen_value_for_matrix_without_eigen_vector (int matrix_size, double ** matrix, double eigenvalue [])
{
  double ** transformation_matrix = new double * [matrix_size];
  for (int i=0; i<matrix_size; i++)
    transformation_matrix [i] = new double [matrix_size+1];
  
  int rc = _householder_tri_diagonalization (matrix_size, transformation_matrix, matrix);
  
  for (int i=0; i<matrix_size-1; i++)
    if (fabs(matrix[i][i+1]) < 1e-7)
      matrix[i+1][i] = matrix[i][i+1] = 0.0;
  
  // initialize the eigenvector and eigen value
  for (int i=0; i<matrix_size; i++)
    eigenvalue [i] =0.0;

  _compute_eigen_value_without_eigen_vector (matrix_size, matrix, eigenvalue, transformation_matrix);

  for (int i=0; i<matrix_size; i++)
    delete [] transformation_matrix [i];
  delete [] transformation_matrix;

  return rc;
}

// end of new method






int Molecule :: _compute_charge_for_fragment 
(double partial_charge[], int atomi_to_huckel_matrix_column_number[], double ** molecule_huckel_matrix, double number_of_electron [], resizable_array<int> * fragment_matrix_no_to_atomi, double ** huckel_matrix, double ** transformation_matrix, double ** eigenvector_matrix, double fragment_electron [], double eigen_value [])
{
  int matrix_size = fragment_matrix_no_to_atomi->number_elements();

  for (int i=0; i<matrix_size; i++)
    {
      int i_to_molecule = atomi_to_huckel_matrix_column_number [fragment_matrix_no_to_atomi->item(i)];

      fragment_electron [i] = number_of_electron [i_to_molecule];

      for (int j=i; j<matrix_size; j++)
        {
          int j_to_molecule = atomi_to_huckel_matrix_column_number [fragment_matrix_no_to_atomi->item(j)];
          huckel_matrix [i][j] = huckel_matrix [j][i] = molecule_huckel_matrix [i_to_molecule][j_to_molecule];
        }
    }
  int rc = _householder_tri_diagonalization (matrix_size, transformation_matrix, huckel_matrix);

  for (int i=0; i<matrix_size-1; i++)
    if (fabs(huckel_matrix[i][i+1]) < 1e-7)
      huckel_matrix[i+1][i] = huckel_matrix[i][i+1] = 0;

#if (DEBUG_PRINT)
  cerr<< "Huckel Matrix "<<endl;
  for (int i=0; i<matrix_size; i++)
    {
      for (int j=0; j<matrix_size;j++)
        cerr<<huckel_matrix [i][j]<<"    ";
      cerr<<endl;
    }

  cerr<<endl<<"Transformation Matrix "<<endl;
  for (int i=0; i<matrix_size-2; i++)
    {
      for (int j=0; j<matrix_size+1;j++)
        cerr<<transformation_matrix [i][j]<<"    ";
      cerr<<endl;
    }

#endif

  // initialize the eigenvector
  for (int i=0; i<matrix_size; i++)
    for (int j=0; j<matrix_size; j++)
      eigenvector_matrix [i][j] = 0;

  _compute_eigenvector (matrix_size, huckel_matrix, eigenvector_matrix, eigen_value, transformation_matrix);

#if (DEBUG_PRINT)
  cerr<< "Eigen Value "<<endl;
  for (int i=0; i<matrix_size; i++)
    cerr<<eigen_value [i]<<"    ";
  cerr<<endl;

  cerr<<endl<<"Eigen Vector "<<endl;
  for (int i=0; i<matrix_size; i++)
    {
      for (int j=0; j<matrix_size;j++)
        cerr<<eigenvector_matrix [i][j]<<"    ";
      cerr<<endl;
    }

#endif

  rc = _assign_huckel_charge_to_atom (fragment_matrix_no_to_atomi, fragment_electron, eigenvector_matrix, partial_charge);

  return rc;
}

int Molecule :: _compute_charge_for_fragment 
(double partial_charge[], int atomi_to_huckel_matrix_column_number[], double ** molecule_huckel_matrix, double number_of_electron [], resizable_array<int> * fragment_matrix_no_to_atomi)
{
  int matrix_size = fragment_matrix_no_to_atomi->number_elements();

  double ** huckel_matrix = new double * [matrix_size];
  double ** transformation_matrix = new double * [matrix_size];
  double ** eigenvector_matrix = new double * [matrix_size];

  for (int i=0; i<matrix_size; i++)
    {
      huckel_matrix [i] = new double [matrix_size];
      transformation_matrix [i] = new double [matrix_size+1];
      eigenvector_matrix [i] = new double [matrix_size];
    }

  double * fragment_electron = new double [matrix_size];
  double * eigen_value = new double [matrix_size];

  int rc = _compute_charge_for_fragment (partial_charge, atomi_to_huckel_matrix_column_number, molecule_huckel_matrix, number_of_electron, fragment_matrix_no_to_atomi, huckel_matrix, transformation_matrix, eigenvector_matrix, fragment_electron, eigen_value);

  for (int i=0; i<matrix_size; i++)
    {
      delete [] transformation_matrix [i];
      delete [] eigenvector_matrix [i];
      delete [] huckel_matrix [i];
    }

  delete [] huckel_matrix;
  delete [] transformation_matrix;
  delete [] eigenvector_matrix;

  delete [] fragment_electron;
  delete [] eigen_value;

  return rc;
}

int Molecule :: _compute_charge_for_fragment 
(double partial_charge[], int atomi_to_huckel_matrix_column_number[], int conjugated_fragment_array[], double ** molecule_huckel_matrix, double number_of_electron [], int fragment_no)
{
  int rc;
  int n_atoms = natoms();
  int number_in_fragment = 0;


  resizable_array<int> fragment_matrix_no_to_atomi;
  fragment_matrix_no_to_atomi.resize(n_atoms);

  for (int i=0; i<n_atoms; i++)
    if ((conjugated_fragment_array [i] == fragment_no) && (atomic_number(i) != 1))
      {
        number_in_fragment ++;
        fragment_matrix_no_to_atomi.add(i);
      }

  rc = _compute_charge_for_fragment (partial_charge, atomi_to_huckel_matrix_column_number, molecule_huckel_matrix, number_of_electron, &fragment_matrix_no_to_atomi);

  fragment_matrix_no_to_atomi.resize(0);
  return rc;
}

int Molecule :: _adjust_with_formal_charge (double partial_charge [], int fragment_array [], int fragment_no)
{
  int n_atoms = natoms();
  if (0==fragment_no)
    for (int i=0; i<n_atoms; i++)
        partial_charge [i] +=formal_charge (i);
  else 
    {
      double total_formal_charge =0;
      int number_of_atom_in_fragment =0;

      for (int i=0; i<n_atoms; i++)
        if ((1 != atomic_number(i)) && (fragment_no == fragment_array[i]))
          {
            total_formal_charge += formal_charge (i);
            number_of_atom_in_fragment ++;
          }
      if (total_formal_charge > 0)
        for (int i=0; i<n_atoms; i++)
          if ((1 != atomic_number(i)) && (fragment_no == fragment_array[i]))
            partial_charge [i] += total_formal_charge / number_of_atom_in_fragment;
    }

  return 1;
}

// this method deal with O=N, O=S, O=P group 
// to conform it to the sybyl definition for the huckel charge calculation
void Molecule :: _specific_modification_of_conjugated_fragment (int conjugated_fragment[], atom_type_t atom_type [])
{
  int n_atom = natoms();
  for (int i=0; i<n_atom; i++)
    if ((atom_type [i] == ATOM_TYPE_O2) && (ncon (i) == 1))
      { 
        int neighbor_type = atom_type [other(i, 0)];
        if ((neighbor_type == ATOM_TYPE_NPL3) || (neighbor_type == ATOM_TYPE_S2) || (neighbor_type == ATOM_TYPE_SO) || (neighbor_type == ATOM_TYPE_SO2))
          conjugated_fragment[i] = 0;
      }
}

int Molecule :: _calculate_Huckel_partial_charges 
(double partial_charge [], int atomi_to_huckel_matrix_column_number [], 
int conjugated_fragment_array [], atom_type_t atom_type [], resizable_array<int> * huckel_matrix_column_number_to_atomi, 
double ** huckel_matrix, double number_of_electron [])
{
  int rc;

  rc = _setup_molecule_huckel_matrix (huckel_matrix_column_number_to_atomi, atomi_to_huckel_matrix_column_number, atom_type, huckel_matrix, number_of_electron);

  if (rc ==0)
    {
      cerr<<"Error in setting up huckel matrix  -- "<<molecule_name()<<endl;
      return 0;
    }
#if (DEBUG_PRINT)
  int matrix_size = huckel_matrix_column_number_to_atomi->number_elements();
  cerr<< "Huckel matrix "<<endl;
  for (int i=0; i<matrix_size; i++)
    {
      for (int j=0; j<matrix_size;j++)
        cerr<<huckel_matrix [i][j]<<"    ";
      cerr<<endl;
    }
  cerr<<"Electron Array "<<endl;
  for (int i=0; i<matrix_size; i++)
    cerr<<number_of_electron [i]<<"    ";
  cerr<<endl;
#endif

  //Find all of the conjugated fragment
  int number_of_fragment = _find_conjugated_fragment ( conjugated_fragment_array);

#if (DEBUG_PRINT)
  for (int i=0; i<natoms(); i++)
    cerr<<"Atom "<<i<<"   Fragment="<<conjugated_fragment_array[i]<<endl;
#endif

  _specific_modification_of_conjugated_fragment (conjugated_fragment_array, atom_type);

  for (int i=1; i<=number_of_fragment; i++)
    rc =_compute_charge_for_fragment (partial_charge, atomi_to_huckel_matrix_column_number, conjugated_fragment_array, huckel_matrix, number_of_electron, i);
  
  // if molecule has formal charge, distribute among conjugated fragments
  if (has_formal_charges())
    for (int i=0; i<=number_of_fragment; i++)
      rc =_adjust_with_formal_charge (partial_charge, conjugated_fragment_array, i);
  return rc;
}

int Molecule :: _calculate_Huckel_partial_charges (double partial_charge [], atom_type_t atom_type [])
{
  int n_atoms = natoms();

  for (int i=0; i<n_atoms; i++)
    partial_charge [i] = 0;

  int * atomi_to_huckel_matrix_column_number = new int [n_atoms];
  int * conjugated_fragment_array = new int [n_atoms];

  resizable_array<int> huckel_matrix_column_number_to_atomi; 
  _setup_conversion_array_between_huckel_matrix_and_atomi (&huckel_matrix_column_number_to_atomi, atomi_to_huckel_matrix_column_number);

  int matrix_size = huckel_matrix_column_number_to_atomi.number_elements();

  double * number_of_electron = new double [matrix_size];
  
  double ** huckel_matrix = new double * [matrix_size];
  for (int i=0; i<matrix_size; i++)
    huckel_matrix [i] = new double [matrix_size];
  
  int rc= _calculate_Huckel_partial_charges(partial_charge, atomi_to_huckel_matrix_column_number, conjugated_fragment_array, atom_type, &huckel_matrix_column_number_to_atomi, huckel_matrix, number_of_electron);
  
  huckel_matrix_column_number_to_atomi.resize(0);

  for (int i=0; i<matrix_size; i++)
    delete [] huckel_matrix [i];

  delete [] huckel_matrix;
  delete [] number_of_electron;
  delete [] conjugated_fragment_array;
  delete [] atomi_to_huckel_matrix_column_number;
  return rc;
}

int
Molecule:: compute_Huckel_partial_charges ()
{
  int n_atoms = natoms();
  atom_type_t * sybyl_atom_type = new unsigned int [n_atoms]; std::unique_ptr<atom_type_t[]> free_sybyl_atom_type(sybyl_atom_type);

  if (NULL == _atom_type)
    {
      allocate_atom_types ();
      assert ( NULL != _atom_type);
    }

  if ("SYBYL" == _atom_type->ztype())
    {
      for (int i=0; i<n_atoms; i++)
        sybyl_atom_type [i] = atom_type (i);
    }
  else 
    {
      int rc = find_simplified_sybyl_atom_type_sybyl_style(sybyl_atom_type);
      if (0==rc)
        {
          cerr<<"Unreconized Atom Type Encountered  -- "<<molecule_name()<<endl;
          return 0;
        }
      if (!_atom_type->ztype().length())
        {
          for (int i=0; i<n_atoms; i++)
            {
              set_atom_type (i, sybyl_atom_type[i]);
            }
          _atom_type->set_type ("SYBYL");
        }
    }
 
  (void) compute_Huckel_partial_charges (sybyl_atom_type);

#if (DEBUG_PRINT)
  for (int i=0; i<natoms(); i++)
    cerr<<"Atom "<<i<<"   Type ="<<sybyl_atom_type[i]<<endl;
#endif
  return 1;
}

int
Molecule:: compute_Huckel_partial_charges (atom_type_t atom_type [])
{
  compute_aromaticity_if_needed();

  int n_atoms = natoms();

  double * partial_charge = new double [n_atoms]; std::unique_ptr<double[]> free_partial_charge(partial_charge);

  int rc = _calculate_Huckel_partial_charges (partial_charge, atom_type);

#if (DEBUG_PRING)
  for (int i=0; i<n_atoms; i++)
    cerr<<partial_charge [i]<<"   ";
  cerr<<endl;
#endif

  // assign _charge;

  for (int i=0; i<n_atoms; i++) 
  {
    set_charge(i, static_cast<charge_t> (partial_charge[i]) );
  }

  // assign _charge_type and array of _charge;
  _charges->set_type ("HUCKEL");

  return rc;
}

int
Molecule:: compute_Gasteiger_Huckel_partial_charges ()
{
  int rc;
  compute_aromaticity_if_needed();
  int n_atoms = natoms();
  
  atom_type_t * sybyl_atom_type = new atom_type_t [n_atoms]; std::unique_ptr<atom_type_t[]> free_sybyl_atom_type(sybyl_atom_type);
  
  if (NULL == _atom_type)
    {
      allocate_atom_types ();
      assert ( NULL != _atom_type);
    }
  
  if ("SYBYL" == _atom_type->ztype())
  {
    for (int i=0; i<n_atoms; i++)
    {
      sybyl_atom_type [i] = atom_type (i);
    }
  }
  else 
  {
    if (0 == find_simplified_sybyl_atom_type_sybyl_style(sybyl_atom_type))
    {
      cerr<<"Unreconized Atom Type Encountered  -- "<<molecule_name()<<endl;
      return 0;
    }

    if (0 == _atom_type->ztype().length())
    {
      for (int i=0; i<n_atoms; i++)
      {
        set_atom_type (i, sybyl_atom_type[i]);
      }
      _atom_type->set_type ("SYBYL");
    }
  }

  rc = compute_Gasteiger_Huckel_partial_charges (sybyl_atom_type);

#if (DEBUG_PRINT)
  for (int i=0; i<natoms(); i++)
    cerr<<"Atom "<<i<<"   Type ="<<sybyl_atom_type[i]<<endl;
#endif
  return rc;
}

int
Molecule:: compute_Gasteiger_Huckel_partial_charges (atom_type_t atom_type [])
{ 
  int n_atoms = natoms();

  compute_aromaticity_if_needed();

  double * partial_charge = new double [n_atoms]; std::unique_ptr<double[]> free_partial_charge(partial_charge);

  if (!_calculate_Huckel_partial_charges (partial_charge, atom_type))
    return 0;
      
  if (! _gasteiger_partial_charge_procedure (partial_charge, atom_type))   // failed
    return 0;     

  // assign _charge;

  for (int i=0; i<n_atoms; i++)
  {
    if ((partial_charge[i]<3.0) && (partial_charge[i]>-3.0)) 
      set_charge(i, static_cast<charge_t> (partial_charge[i]) );
    else 
      return 0;
  }

  _charges->set_type("GAST_HUCK");
  
  return 1;
}
