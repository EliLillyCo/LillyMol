#define SET_OF_COMMA_DESCRIPTOR 4
void comma2_related_descriptor_computation_procedure (Molecule &, int, double [], double [], double **, double [], double **);
void compute_moment_of_inertia_using_center_and_axis (Molecule &, int, const double ** , double [], double **, double [], double []);
void compute_moment_of_inertia (Molecule &, int, const double * const *, const double [], const double [], double **, double []);
void compute_field_center (int, const double * const *, double [], const double []);
void output_ms_comma_result_header_to_buffer (IWString &);
