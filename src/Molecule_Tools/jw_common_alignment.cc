#include <iostream>

#include "jw_common_alignment.h"

using std::cerr;
using std::endl;

int
orient_molecule_to_rotational_axis (Molecule &m, double ** coordinates, double ** original_matrix, 
				    double temp_array [], double ** axis)
{
  int n_atoms = m.natoms();

  for (int i=0; i<n_atoms; i++)
    {
      coordinates [i][0] = m.x(i);
      coordinates [i][1] = m.y(i);
      coordinates [i][2] = m.z(i);
    }

  // initialize the matrices
  for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
	{
	  original_matrix [i][j] = 0;
	  axis [i][j] = 0;
	}
      temp_array [i] = 0;
    }

  // compute the center, stored in temp_array
  for (int i=0; i<n_atoms; i++) 
    for (int j=0; j<3; j++)
      temp_array [j] += coordinates [i][j];
  
  for (int i=0; i<3; i++)
    temp_array [i] = temp_array [i] / (double) n_atoms;
  
  m.translate_atoms (-temp_array [0], -temp_array [1], -temp_array [2]);
  
  for (int i=0; i<n_atoms; i++)
    {
      coordinates [i][0] = m.x(i);
      coordinates [i][1] = m.y(i);
      coordinates [i][2] = m.z(i);
    }
  
  for (int i=0; i<n_atoms; i++)
    {
      double tmp0 = coordinates [i][0];
      double tmp1 = coordinates [i][1];
      double tmp2 = coordinates [i][2];
      
      double tmp0_square = tmp0 * tmp0;
      double tmp1_square = tmp1 * tmp1;
      double tmp2_square = tmp2 * tmp2;

      original_matrix [0][0] += (tmp1_square + tmp2_square);
      original_matrix [1][1] += (tmp0_square + tmp2_square);
      original_matrix [2][2] += (tmp0_square + tmp1_square);

      original_matrix [0][1] -= tmp0 * tmp1;
      original_matrix [1][2] -= tmp1 * tmp2;
      original_matrix [0][2] -= tmp0 * tmp2;
    }

  /* old way of computation.  Something is not quite right
  // compute the original_matrix
  for (int i=0; i<n_atoms; i++)
    {
      // diagonal elemetns
      for (int j=0; j<3; j++)
	original_matrix [j][j] = coordinates [i][j] * coordinates [i][j];
      
      original_matrix [0][1] -= coordinates [i][0] * coordinates [i][1];
      original_matrix [1][2] -= coordinates [i][1] * coordinates [i][2];
      original_matrix [0][2] -= coordinates [i][0] * coordinates [i][2];
    }
  
  */

  original_matrix [1][0] = original_matrix [0][1];
  original_matrix [2][1] = original_matrix [1][2];
  original_matrix [2][0] = original_matrix [0][2];
  
  (void) m.find_eigen_value_for_matrix (3, original_matrix, temp_array, axis);
  
  for (int i=0; i<n_atoms; i++)
    {
      double tempx = 0;
      double tempy = 0;
      double tempz = 0;
      
      for (int j=0; j<3; j++)
	{
	  tempx += coordinates [i][j] * axis [j][0];
	  tempy += coordinates [i][j] * axis [j][1];
	  tempz += coordinates [i][j] * axis [j][2];
	}
      m.setxyz(i, -tempz, tempy, tempx);
      //m.setxyz(i, tempx, tempy, tempz);
    }

  return 1;
}
