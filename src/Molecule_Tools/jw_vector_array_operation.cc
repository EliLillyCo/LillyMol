/* 
 * Funcitons and subroutines from fortran code (connolly.f) 
 * All of the arrays are 3-dimension
 */
#include "jw_vector_array_operation.h"

double
distance_square (const double a [3], const double b [3])
{
  double distance_squa =0.0;
  for (int i =0; i< 3; i++)
  {
    distance_squa += (a[i]-b[i]) * (a[i]-b[i]);
  }
  return distance_squa;
}

double
jw_distance (const double a [3], const double b [3])
{
  return sqrt(distance_square (a, b));
} 

/* 
 * this function is a computational inexpensive version of distance check
 * to see if two points are within distance R   true if <=R
 */
int
two_points_within_distance (const double i_coordinate[3], const double j_coordinate[3], double dist)
{
  double dx = fabs(i_coordinate[0] - j_coordinate[0]);
  if (dx > dist)
    return 0;

  double dy = fabs(i_coordinate[1] - j_coordinate[1]);
  if (dy > dist)
    return 0;

  double dz = fabs(i_coordinate[2] - j_coordinate[2]);
  if (dz > dist)
    return 0;

  return sqrt(dx * dx + dy * dy + dz * dz) < dist;
}

// true if distance <R 
int
two_points_less_than_distance (const double i_coordinate[3],
                               const double j_coordinate[3],
                               double dist)
{
  double dx = i_coordinate[0] - j_coordinate[0];
  if (dx > dist)
    return 0;
  else if (dx < 0.0)
  {
    dx = -dx;
    if (dx > dist)
      return 0;
  }

  double dy = i_coordinate[1] - j_coordinate[1];
  if (dy > dist)
    return 0;
  else if (dy < 0.0)
  {
    dy = -dy;
    if (dy > dist)
      return 0;
  }

  double dz = i_coordinate[2] - j_coordinate[2];
  if (dz > dist)
    return 0;
  else if (dz < 0.0)
  {
    dz = -dz;
    if (dz > dist)
      return 0;
  }

  return ((dx * dx + dy * dy + dz * dz) < (dist * dist));   // avoid sqrt computation
}

double vector_length (const double a [3])
{
  double length_square =0.0;
  for (int i =0; i< 3; i++)
    {
      length_square += (a[i] * a[i]);
    }
  return sqrt (length_square);
}

double dot (const double a [3], const double b [3])
{
  double dot_value = 0.0;
  for (int i=0; i<3; i++)
    {
      dot_value += a[i]*b[i];
    }
  return dot_value;
}

void vector_cross (const double a [3], const double b [3], double c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

void multiply_vector (const double v [3], const double a[3][3], double w[3])
{
  for (int i=0; i<3; i++)
    {
      w[i]=0.0;
      for (int j=0; j<3; j++)
	w[i] += a[i][j] * v[j]; 
    }
}

void vector_normalize ( const double a [3], double b[3])
{
  double a_length = vector_length (a);
  if (a_length !=0.0)
    for (int i=0; i<3; i++)
      {
	b[i] = a[i] / a_length;
      }
  else
    for (int i=0; i<3; i++)
      {
	b[i] = a[i];
      }
}


double dabs (double x)
{
  if (x>=0.0) return x;
  return -x;
} 

void vector_perp(const double a [3], double b [3])
{
  double p[3];
  double small = dabs(a[0]);
  int m = 0;
  for (int i=1; i<3; i++)
    {
      if (dabs(a[i]) >= small)
        continue;
      small = dabs(a[i]);
      m = i;
    }
  for (int i=0; i<3; i++)
    {
      p[i] = 0.0;
      if (i== m) 
        p[i] = 1.0;
    }
  double dt = a[m] / ((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]));

  for (int i=0; i<3; i++)
  {
    p[i] = p[i] - dt * a[i];
  }

  vector_normalize (p, b);
}

void concatenate_matrix (double a [3][3], const double b[3][3])
{
  double temp [3][3];
  for (int i =0; i<3; i++)
    for (int j=0; j<3; j++)
      {
	temp[i][j] = 0.0;
	for (int k=0; k<3; k++)
	  temp[i][j] += a[i][k] * b[k][j];
      }
  for (int i =0; i<3; i++)
    for (int j=0; j<3; j++)
      a[i][j] = temp [i][j];
}

void identity_matrix (double h[3][3])
{
  for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
      {
	h[i][j] =0.0;
      }
      h[i][i] = 1.0;
    }
}

void conjugate_matrix (double h[3][3], double g[3][3], double ghgt[3][3])
{
  double gt [3][3];
  identity_matrix (ghgt);
  concatenate_matrix (ghgt, g);
  concatenate_matrix (ghgt, h);

  for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	gt [i][j] = g [j][i];

  concatenate_matrix (ghgt, gt);
}

double det (const double a [3], const double b [3], const double c [3])
{
  double temp [3];
  vector_cross (a, b, temp);
  return dot (temp, c);
}

// assign coordinate to array 

void
assign_coordinate_to_array (const Molecule &m, double coordinate_array [], int ith_atom)
{
  m.atomi(ith_atom)->getxyz(coordinate_array);
  return;
}


void
rotate_around_axis(const Coordinates & axis, Coordinates & point, angle_t theta)
{
  if (0.0 == theta)
    return;

// The direction cosines for the vector
 
  const coord_t dc1 = axis.x ();
  const coord_t dc2 = axis.y ();
  const coord_t dc3 = axis.z ();
//cerr << "Dc's are " << dc1 << "," << dc2 << "," << dc3 << " sum = " <<
//     dc1 * dc1 + dc2 * dc2 + dc3 * dc3 << ", angle " << theta << endl;
  
// Rather than deal properly with a matrix, just use individual variables

  coord_t rotmat11 = cos(theta) + dc1 * dc1 * (1.0 - cos(theta));
  coord_t rotmat12 = dc1 * dc2 * (1.0 - cos(theta)) - dc3 * sin(theta);
  coord_t rotmat13 = dc1 * dc3 * (1.0 - cos(theta)) + dc2 * sin(theta);
  coord_t rotmat21 = dc1 * dc2 * (1.0 - cos(theta)) + dc3 * sin(theta);
  coord_t rotmat22 = cos(theta) + dc2 * dc2 * (1.0 - cos(theta));
  coord_t rotmat23 = dc2 * dc3 * (1.0 - cos(theta)) - dc1 * sin(theta);
  coord_t rotmat31 = dc3 * dc1 * (1.0 - cos(theta)) - dc2 * sin(theta);
  coord_t rotmat32 = dc3 * dc2 * (1.0 - cos(theta)) + dc1 * sin(theta);
  coord_t rotmat33 = cos(theta) + dc3 * dc3 * (1.0 - cos(theta));

  coord_t x0 = point.x();
  coord_t y0 = point.y();
  coord_t z0 = point.z();
  
  coord_t xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
  coord_t yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
  coord_t zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

  point.setxyz (xx, yy, zz);
  return;
}

void
rotate_coordinates_around_axis (const Coordinates & p1, const Coordinates & p2, 
				Coordinates & h_atom, angle_t theta)
{
  if (0.0==theta) return;

  double distance = p2.distance(p1);

  Coordinates axis = p2-p1;
  axis /= distance;

  Coordinates center_of_rotation = p1 + axis * ((h_atom-p1).dot_product(p2-p1)/distance);

  h_atom -= center_of_rotation;
  rotate_around_axis (axis, h_atom, theta);
  h_atom += center_of_rotation;
} 

void
rotate_coordinates_around_axis(const Coordinates & p1, const Coordinates & p2, 
			       resizable_array_p<Coordinates> &hydrogens_array, angle_t theta)
{
  if (0.0==theta) return;

  int n_of_coord = hydrogens_array.number_elements();
  for (int i=0; i<n_of_coord; i++)
    rotate_coordinates_around_axis(p1, p2, * (hydrogens_array.item(i)), theta);

  return;
} 

