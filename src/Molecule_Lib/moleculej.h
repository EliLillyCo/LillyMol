// the following class are used to hold all the neighbors of each atom of one molecule
class arrays_holder 
{
public:
  resizable_array<int> one_bond_neighbors;
  resizable_array<int> two_bond_neighbors;
  resizable_array<int> two_bond_neighbors_predecessor;  // keep track where two_bond_neighbors come from
  resizable_array<int> three_bond_neighbors;
  resizable_array<int> three_bond_neighbors_predecessor;  // keep track where three_bond_neighbors come from
};

int _calculate_Huckel_partial_charges (double [], atom_type_t []);
int _calculate_Huckel_partial_charges (double [], int [], int [], atom_type_t [], resizable_array<int> *, double **, double []);
int _compute_charge_for_fragment (double [], int [], int [], double **, double [], int);
int _compute_charge_for_fragment (double [], int [], double **, double [], resizable_array<int> *);
int _compute_charge_for_fragment (double [], int [], double **, double [], resizable_array<int> *, double **, double **, double **, double [], double []);
int _adjust_with_formal_charge (double [], int [], int);

void _setup_conversion_array_between_huckel_matrix_and_atomi (resizable_array<int> *, int []);
int _setup_molecule_huckel_matrix (resizable_array<int> *, int [], atom_type_t [], double **, double []);
int _setup_huckel_element_per_bond (const Bond *, resizable_array<int> *, int [], atom_type_t [], double **);
void _setup_huckel_diagonal_element (int, resizable_array<int> *, atom_type_t [], double **, double []);
int _householder_tri_diagonalization (int, double **, double **);
int _householder_tri_diagonalization (int, double **, double **, double [], double [], double []);
void _compute_eigenvector (int, double **, double **, double [], double **);
void _compute_eigen_value_without_eigen_vector (int, double **, double [], double **);
int _number_of_sign_agreement (int, double, double [], double []);
int _number_of_sign_agreement (int, double, double [], double [], double []);
void _compute_sub_eigen_value(int, int, int, int, double[], double [], double, double, double []);
void _compute_sub_eigen_value_with_low_accuracy(int, int, int, int, double[], double [], double, double, double []);
double _infinite_norm (int, double **);
void _sort_eigenvalue_eigenvector (int, double [], double **);
void _sort_eigenvalue_eigenvector (int, double [], double **, double []);
void _sort_eigenvalue_without_eigenvector (int, double []);

void _triangulation (int, double **);
void _compute_sub_eigenvector(int, double [], double [], double [], double **);
void _compute_sub_eigenvector(int, double [], double [], double [], double **, double **);
void _compute_eigenvector_procedure (int, double **, double **, double [], int, int, double [], double [], double [], double [], double **);
void _compute_eigen_value_without_eigen_vector_procedure (int, double **, double [], int, int, double [], double[], double[]);
void _transform_eigen_vector (int, double **, double **, double **);
void _inverse_iteration (int, double **, double []);
void _inverse_iteration (int, double **, double [], double []);
void _normalize_vector (int, double []);
int _find_conjugated_fragment (int []);
int _find_conjugated_fragment_for_atom(int, int[], int &, int);
void _specific_modification_of_conjugated_fragment (int [], atom_type_t []);

int _assign_huckel_charge_to_atom (resizable_array<int> *, double [], double **, double []);
int _assign_huckel_charge_to_atom (resizable_array<int> *, double [], double **, double [], double []);

double __a_const (atomic_number_t, atomic_number_t);
double _a_const (atom_number_t, atom_number_t);
int _row_number (atom_number_t);

int _find_neighbors_for_molecule (arrays_holder []);
int _calculate_one_bond_partial_charge (arrays_holder [], double[]);
int _calculate_two_and_three_bond_partial_charge(arrays_holder [], double[]);
int _calculate_two_and_three_bond_partial_charge(arrays_holder [], double[], double [], double []);

int _calculate_abraham_partial_charge (double[]);
int _calculate_abraham_partial_charge (double[], arrays_holder []);
int _gasteiger_row_number (atom_number_t);
//int _gasteiger_partial_charge_procedure (double[]);
int _gasteiger_partial_charge_procedure (double[], atom_type_t []);
int _gasteiger_partial_charge_procedure (double [], double [], Set_of_Atoms [], /*resizable_array<int> [], */ double [], double [],double [], double [], double [], atom_type_t []);
int _pi_gasteiger_partial_charge_computation_procedure (double [], atom_type_t []);
void _assign_initial_charge (double [], atom_type_t []);
void _find_number_of_unsaturated_carbon_of_nitrogen_neighbor (atom_number_t, int &, int &, atom_number_t &);
public: 
int find_simplified_sybyl_atom_type (atom_type_t []);
int find_simplified_sybyl_atom_type_sybyl_style (atom_type_t []);
int find_simplified_sybyl_atom_type_sybyl_style (int []);
int find_mcs_atom_type_similar_to_sybyl (int []);
int find_eigen_value_for_matrix (int, double **, double [], double **); 
int find_eigen_value_for_matrix_without_eigen_vector (int, double **, double []);
int compute_Gasteiger_partial_charges(atom_type_t [], int);
int compute_Gasteiger_partial_charges_with_pi_charges ();
int compute_Gasteiger_Huckel_partial_charges (atom_type_t []);
int compute_Huckel_partial_charges(atom_type_t []);
private:
