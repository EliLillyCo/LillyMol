/* 
 * Funcitons and subroutines from fortran code (connolly.f) 
 * All of the arrays are 3-dimension
 */
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/coordinates.h"
#include "Molecule_Lib/iwmtypes.h"

double distance_square (const double [3], const double [3]);

double jw_distance (const double [3], const double [3]);

/* 
 * this function is a computational inexpensive version of distance check
 * to see if two points are within distance R   true if <=R
 */

int two_points_within_distance (const double [3], const double [3], double);

// true if distance <R 
int two_points_less_than_distance (const double [3], const double [3], double dist);

double vector_length (const double [3]);

double dot (const double [3], const double [3]);

void vector_cross (const double [3], const double [3], double [3]);

void multiply_vector (const double [3], const double [3][3], double [3]);

void vector_normalize ( const double [3], double [3]);

double dabs (double);

void vector_perp(const double [3], double [3]);

void concatenate_matrix (double [3][3], const double [3][3]);

void identity_matrix (double [3][3]);

void conjugate_matrix (double [3][3], double [3][3], double [3][3]);

double det (const double [3], const double [3], const double [3]);

// assign coordinate to array
void
assign_coordinate_to_array (const Molecule &, double [], int);

void
rotate_around_axis(const Coordinates &, Coordinates &, angle_t);

void
rotate_coordinates_around_axis(const Coordinates &, const Coordinates&, resizable_array_p<Coordinates> &, angle_t); 

void
rotate_coordinates_around_axis(const Coordinates &, const Coordinates&, Coordinates &, angle_t); 
