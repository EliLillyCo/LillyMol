#include "molecule.h"
#include "jw_vector_array_operation.h"
#include "qry_wstats.h"
#include "accumulator.h"
#include "istream_and_type.h"
#include "comma_related_procedures.h"

/*
 * compute the momet of inertia using center and axis as frame
 */
void compute_moment_of_inertia_using_center_and_axis (Molecule &m, int n_points, double ** coordinate, double center [], double **axis, 
						      double field [], double field_moment_of_inertia [])
{
  double ** field_moment_of_inertia_tensor = new double * [3];
  double ** field_moment_of_inertia_axis = new double * [3];

  for (int i=0; i<3; i++)
  {
    field_moment_of_inertia_tensor[i] = new double [3];
    field_moment_of_inertia_axis[i] = new double [3];
  }

  // initialize the array
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      field_moment_of_inertia_tensor [i][j] = 0.0;
  
  // compute the moment of inertia tensor for field around center using axis as frame;
  for (int i=0; i<n_points; i++)
    {
      double tmp[3];
      double tmp_after_rotation[3];
      for (int j=0; j<3; j++)
	tmp[j] = coordinate[i][j] - center[j];
      
      for (int j=0; j<3; j++)
      {
        tmp_after_rotation[j] = 0.0;
        for (int k=0; k<3; k++)
          tmp_after_rotation[j] += tmp[k]* axis[k][j];
      }

      double tmp0_square = tmp_after_rotation[0] * tmp_after_rotation[0];
      double tmp1_square = tmp_after_rotation[1] * tmp_after_rotation[1];
      double tmp2_square = tmp_after_rotation[2] * tmp_after_rotation[2];
      
      double tmp_field = field [i];
      
      field_moment_of_inertia_tensor[0][0] += (tmp1_square + tmp2_square) * tmp_field;
      field_moment_of_inertia_tensor[1][1] += (tmp0_square + tmp2_square) * tmp_field;
      field_moment_of_inertia_tensor[2][2] += (tmp0_square + tmp1_square) * tmp_field;
      
      field_moment_of_inertia_tensor[0][1] -= tmp_after_rotation[0] * tmp_after_rotation[1] * tmp_field;
      field_moment_of_inertia_tensor[1][2] -= tmp_after_rotation[1] * tmp_after_rotation[2] * tmp_field;
      field_moment_of_inertia_tensor[0][2] -= tmp_after_rotation[0] * tmp_after_rotation[2] * tmp_field;
    }
  
  field_moment_of_inertia_tensor[1][0] = field_moment_of_inertia_tensor[0][1];
  field_moment_of_inertia_tensor[2][1] = field_moment_of_inertia_tensor[1][2];
  field_moment_of_inertia_tensor[2][0] = field_moment_of_inertia_tensor[0][2];

  (void) m.find_eigen_value_for_matrix (3, field_moment_of_inertia_tensor, field_moment_of_inertia, field_moment_of_inertia_axis);

  for (int i=0; i<3; i++)
    {
      delete [] field_moment_of_inertia_axis [i];
      delete [] field_moment_of_inertia_tensor [i];
    }
  
  delete [] field_moment_of_inertia_axis;
  delete [] field_moment_of_inertia_tensor;
}

/*
 * compute the moment of inertia
 */

void
compute_moment_of_inertia (Molecule &m,
                           int n_points,
                           const double * const * coordinate,
                           const double field_center [],
                           const double field [],
                           double ** moment_of_inertia_axis, 
                           double moment_of_inertia [])
{
  double ** moment_of_inertia_tensor = new double * [3];
  
  for (int i=0; i<3; i++)
    moment_of_inertia_tensor[i] = new double [3] ;
  
  // initialize the array
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      moment_of_inertia_tensor[i][j] = 0.0;
  
  // compute the moment of inertia tensor for field
  for (int i=0; i<n_points; i++)
    {
      double tmp0 = coordinate[i][0] - field_center[0];
      double tmp1 = coordinate[i][1] - field_center[1];
      double tmp2 = coordinate[i][2] - field_center[2];
      
      double tmp0_square = tmp0 * tmp0;
      double tmp1_square = tmp1 * tmp1;
      double tmp2_square = tmp2 * tmp2;
      
      double tmp_field = field [i];
      
      moment_of_inertia_tensor[0][0] += (tmp1_square + tmp2_square) * tmp_field;
      moment_of_inertia_tensor[1][1] += (tmp0_square + tmp2_square) * tmp_field;
      moment_of_inertia_tensor[2][2] += (tmp0_square + tmp1_square) * tmp_field;
      
      moment_of_inertia_tensor[0][1] -= tmp0 * tmp1 * tmp_field;
      moment_of_inertia_tensor[1][2] -= tmp1 * tmp2 * tmp_field;
      moment_of_inertia_tensor[0][2] -= tmp0 * tmp2 * tmp_field;
    }
  
  moment_of_inertia_tensor[1][0] = moment_of_inertia_tensor[0][1];
  moment_of_inertia_tensor[2][1] = moment_of_inertia_tensor[1][2];
  moment_of_inertia_tensor[2][0] = moment_of_inertia_tensor[0][2];

  (void) m.find_eigen_value_for_matrix (3, moment_of_inertia_tensor, moment_of_inertia, moment_of_inertia_axis);
  
  for (int i=0; i<3; i++)
    delete [] moment_of_inertia_tensor [i];
  
  delete [] moment_of_inertia_tensor;
}
 
/*
 * compute the field center, (the place where the first order moment is zero)
 */
void
compute_field_center (int n_points,
                      const double * const * coordinate,
                      double field_center [],
                      const double field_array [])
{
  // initialize the four array
  for (int i=0; i<3; i++)
    field_center[i] = 0.0;

  // compute the field_center 
  for (int i=0; i<n_points; i++)
    for (int j=0; j<3; j++)
      field_center[j] += coordinate[i][j] * field_array[i];
  
  for (int i=0; i<3; i++)
    field_center[i] = field_center[i] / n_points;
}


/*
 * The procedrue that is acutally do the working
 */
void comma2_related_descriptor_computation_procedure (Molecule &m, int n_points, double field [], double comma2_descriptor [], 
                                                      double ** coordinate,
                                                      double geometry_center [],
                                                      double ** geometry_moment_of_inertia_axis)
{
  double field_center[3];
  double field_center_to_geometry_center[3];

  double field_moment_of_inertia[3];
  double field_moment_of_inertia_geometry_center_frame[3];

  double ** field_moment_of_inertia_axis = new double * [3];

  for (int i=0; i<3; i++)
    field_moment_of_inertia_axis [i] = new double [3];

  compute_field_center(n_points, coordinate, field_center, field);
  
  for (int i=0; i<3; i++)
    field_center_to_geometry_center[i] = field_center[i] - geometry_center[i];

  compute_moment_of_inertia(m, n_points, coordinate, field_center, field, field_moment_of_inertia_axis, field_moment_of_inertia);

  compute_moment_of_inertia_using_center_and_axis(m, n_points, coordinate, geometry_center, geometry_moment_of_inertia_axis, field, field_moment_of_inertia_geometry_center_frame);

  //assign comma descriptor value
  comma2_descriptor[0] = vector_length(field_center_to_geometry_center);

  for (int i=0; i<3; i++)
    {
      comma2_descriptor[1+i] = field_moment_of_inertia[i];
      comma2_descriptor[4+i] = field_moment_of_inertia_geometry_center_frame[i];
    }
  
  for (int i=0; i<3; i++)
    delete [] field_moment_of_inertia_axis [i];
  
  delete [] field_moment_of_inertia_axis;

  return;
}

void output_ms_comma_result_header_to_buffer (IWString &buffer)
{
  for (int i=0; i<SET_OF_COMMA_DESCRIPTOR*7; i++)
    buffer<<"whim_JWMSCM"<<i<<" ";
}
