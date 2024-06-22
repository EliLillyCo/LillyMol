/*
 * This program provide the program for the computation of molconn descriptors
 */

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "e_state_computation_procedure.h"
#include "jw_path_with_values.h"

using std::cerr;

// DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

#define MAX_PATH_NUMBER 100000
// for each path in the computation, these five values are computed at each step
#define NUMBER_OF_VALUES_FOR_PATH 5

static int compute_additional_connectivity_related_descriptors = 0;
static int max_allowed_heavy_atom = 100;

static double simple_linear_path_chi[11];
static double valence_linear_path_chi[11];
static double modified_valence_linear_path_chi[11];

//  When writing to a database we need chemical standardisation

static Chemical_Standardisation chemical_standardisation;

int verbose = 0;

static int molecules_read = 0;

// keep track of number of warning and number of error
static int number_of_error = 0;

#ifdef __GNUG__
template class resizable_array_p<Path_with_values>;
template class resizable_array_base<Path_with_values*>;
#endif

static int output_precision = 3;

static double truncate_to_zero = 1.0e-05;

// Normally we read a smiles file and produce a descriptor file.
// Optionally we can read a descriptor file with a smiles as the first
// column, and optionally, write the same thing.
static int read_descriptor_file_pipeline = 0;
static int write_descriptor_file_pipeline = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -l             compute more connectivity related descriptors besides those in MOLCONNZ program\n";
  cerr << "  -m <n>         maximum allowed heavy atom number for input molecule (default 100)\n";
  cerr << "  -i <type>      input type\n";
  cerr << "  -Y ...         other options, enter '-Y help' for details\n";
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
preprocess_molecule(Molecule& m) {
  //  only deal with the graph of no hydrogen atoms

  m.remove_all(1);
  m.reduce_to_largest_fragment();  // always reduce to largest fragment

  // make sure that there are at least one bond
  if (m.natoms() < 2) {
    return 0;
  }

  // do not process big molecules
  if (m.natoms() > max_allowed_heavy_atom) {
    return 0;
  }

  return m.organic_only();
}

static int
WriteHeader(IWString_and_File_Descriptor& output) {
  // First descriptor does not have leading space.
  output << "jwmc_ncirc";

  for (int i = 0; i < 11; i++) {
    output << " jwmc_xp" << i;
  }

  for (int i = 0; i < 11; i++) {
    output << " jwmc_xpv" << i;
  }

  output << " jwmc_xc3 jwmc_xc4 jwmc_xpc4";

  for (int i = 3; i < 11; i++) {
    output << " jwmc_xch" << i;
  }

  output << " jwmc_xcv3 jwmc_xcv4 jwmc_xpcv4";

  for (int i = 3; i < 11; i++) {
    output << " jwmc_xchv" << i;
  }

  output << " jwmc_dx0";

  for (int i = 1; i < 11; i++) {
    output << " jwmc_dxp" << i;
  }

  output << " jwmc_dxv0";

  for (int i = 1; i < 11; i++) {
    output << " jwmc_dxpv" << i;
  }

  output << " jwmc_ka0";
  for (int i = 1; i < 7; i++) {
    output << " jwmc_kappa" << i;
  }

  output << " jwmc_si jwmc_totops jwmc_sumi jwmc_sumdeli jwmc_tets2 jwmc_phia "
            "jwmc_logidw jwmc_idwbar jwmc_logidc jwmc_idcbar jwmc_w jwmc_logwp jwmc_pf "
            "jwmc_wt jwmc_knotp jwmc_knotvp jwmc_nclass jwmc_numhbd jwmc_numhba";

  for (int i = 1; i < 11; i++) {
    output << " jwmc_nxxp" << i;
  }

  output << " jwmc_nxc3 jwmc_nxc4 jwmc_nxxpc4 jwmc_logntpath jwmc_tm3";

  if (compute_additional_connectivity_related_descriptors) {
    output << " jwmc_jwred jwmc_jwrsymm jwmc_jwlsymm";

    // modified valence path chi
    for (int i = 0; i < 15; i++) {
      output << " jwmc_jwcmvp" << i;
    }

    // modified valence chain chi
    for (int i = 0; i < 8; i++) {
      output << " jwmc_jwcmvc" << i;
    }

    // modified valence difference chi
    for (int i = 0; i < 11; i++) {
      output << " jwmc_jwcmvd" << i;
    }

    // 1st order Zagreb index
    output << " jwmc_jwzas1 jwmc_jwzav1 jwmc_jwzamv1";

    // second order Zagreb index
    output << " jwmc_jwzas2 jwmc_jwzav2 jwmc_jwzamv2";
  }

  output << '\n';

  return 1;
}

static int
atom_valence(atomic_number_t atomic_number_i) {
  int valence = 1;

  // C
  if (6 == atomic_number_i) {
    valence = 4;
  }
  // F Cl Br I
  else if ((9 == atomic_number_i) || (17 == atomic_number_i) || (35 == atomic_number_i) ||
           (53 == atomic_number_i)) {
    valence = 7;
  }
  // O S
  else if ((8 == atomic_number_i) || (16 == atomic_number_i)) {
    valence = 6;
  }
  // N P
  else if ((7 == atomic_number_i) || (15 == atomic_number_i)) {
    valence = 5;
  }
  // H
  else if (1 == atomic_number_i) {
    valence = 1;
  }

  return valence;
}

static int
prime_quantum_number(atomic_number_t atomic_number_i) {
  int pq_number = 1;

  // C,N,O,F
  if ((atomic_number_i > 5) && (atomic_number_i < 10)) {
    pq_number = 2;
  }
  // P,S,Cl
  else if ((atomic_number_i > 14) && (atomic_number_i < 18)) {
    pq_number = 3;
  }
  // Br
  else if (35 == atomic_number_i) {
    pq_number = 4;
  }
  // I
  else if (53 == atomic_number_i) {
    pq_number = 5;
  }
  // H
  // H sould not even be the molecule
  else if (1 == atomic_number_i) {
    pq_number = 2;
  }

  return pq_number;
}

static double
find_alpha_value_for_molecule(Molecule& m, const Atom* const* atoms,
                              const atomic_number_t* z, const int* hcount) {
  double alpha = 0.0;

  int n_atoms = m.natoms();

  for (int i = 0; i < n_atoms; i++) {
    const Atom* atomi = atoms[i];
    int nconi = atomi->ncon();
    int nbondsi = atomi->nbonds();

    int non_sigma_electrons = atom_valence(z[i]) - nconi - hcount[i];

    switch (z[i]) {
      case 6:
        if (nconi == nbondsi - 1) {
          alpha += -0.13;
        } else if (nconi == nbondsi - 2) {
          alpha += -0.22;
        }
        break;
      case 7:
        if ((1 == non_sigma_electrons) || (2 == non_sigma_electrons)) {
          alpha += -0.04;
        } else if (3 == non_sigma_electrons) {
          alpha += -0.20;
        } else if (4 == non_sigma_electrons) {
          alpha += -0.29;
        } else {
          alpha += -0.04;
        }
        break;
      case 8:
        if (nconi == nbondsi) {
          alpha += -0.04;
        } else if (nconi == nbondsi - 1) {
          alpha += -0.20;
        } else {
          alpha += -0.04;
        }
        break;
      case 9:
        alpha += -0.07;
        break;
      case 15:
        if ((3 == nconi) || (4 == nconi)) {
          alpha += 0.43;
        } else if (2 == nconi) {
          alpha += 0.30;
        } else {
          alpha += 0.20;
        }
        break;
      case 16:
        if ((2 == nconi) || (3 == nconi) || (4 == nconi)) {
          alpha += 0.35;
        } else if (1 == nconi) {
          alpha += 0.22;
        } else {
          alpha += 0.15;
        }
        break;
      case 17:
        alpha += 0.29;
        break;
      case 35:
        alpha += 0.48;
        break;
      case 53:
        alpha += 0.73;
        break;
      default:
        break;
    }
  }
  return alpha;
}

static int
number_of_two_bond_count(Molecule& m, const Atom* const* atoms) {
  int n_atoms = m.natoms();

  int number_of_two_bond_fragment = 0;
  for (int i = 0; i < n_atoms; i++) {
    const Atom* atomi = atoms[i];

    int non_hydrogen_atom_connected = atomi->ncon();  // - m.explicit_hydrogens(i);
    number_of_two_bond_fragment +=
        non_hydrogen_atom_connected * (non_hydrogen_atom_connected - 1) / 2;
  }

  return number_of_two_bond_fragment;
}

static int
number_of_three_bond_count(Molecule& m, const Atom* const* atoms) {
  int nb = m.nedges();

  int number_of_three_bond_fragment = 0;

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    atom_number_t atom1 = b->a1();
    atom_number_t atom2 = b->a2();

    const Atom* atom_1 = atoms[atom1];
    const Atom* atom_2 = atoms[atom2];

    int non_hydrogen_connection_1 = atom_1->ncon();  // - m.explicit_hydrogens(atom1);
    int non_hydrogen_connection_2 = atom_2->ncon();  // - m.explicit_hydrogens(atom2);
    number_of_three_bond_fragment +=
        (non_hydrogen_connection_1 - 1) * (non_hydrogen_connection_2 - 1);
  }

  return number_of_three_bond_fragment;
}

static double
shannon_information_content(Molecule& m, int n_classes, int& n_symmetry_class,
                            int& n_lone_symmetry_class, int symmetry_class_count[]) {
  int n_atoms = m.natoms();

  set_vector(symmetry_class_count, n_classes, 0);

  const int* symmetry_classes_array = m.symmetry_classes();

  for (int i = 0; i < n_atoms; i++) {
    symmetry_class_count[symmetry_classes_array[i] - 1]++;
  }

  int non_hydrogen_atom = m.natoms();  // - m.natoms(1);

  double shannon_content = 0.0;

  n_symmetry_class = 0;
  n_lone_symmetry_class = 0;

  for (int i = 0; i < n_classes; i++) {
    if (symmetry_class_count[i] != 0) {
      if (1 == symmetry_class_count[i]) {
        n_lone_symmetry_class++;
      }

      n_symmetry_class++;
      double probability = (double)symmetry_class_count[i] / (double)non_hydrogen_atom;
      shannon_content -= probability * log10(probability);
    }
  }

  return (shannon_content);
}

static double
shannon_information_content(Molecule& m, int& n_symmetry_class,
                            int& n_lone_symmetry_class) {
  int n_classes = m.number_symmetry_classes();

  int* symmetry_class_count = new int[n_classes];
  std::unique_ptr<int[]> free_symmetry_class_count(symmetry_class_count);

  return shannon_information_content(m, n_classes, n_symmetry_class,
                                     n_lone_symmetry_class, symmetry_class_count);
}

static void
compute_kappa_descriptors(Molecule& m, double kappa_simple_descriptors[],
                          double kappa_alpha_descriptors[], const Atom* const* atoms,
                          const atomic_number_t* z, const int* hcount) {
  double alpha = find_alpha_value_for_molecule(m, atoms, z, hcount);

  int one_bond_count = m.nedges();
  int two_bond_count = number_of_two_bond_count(m, atoms);
  int three_bond_count = number_of_three_bond_count(m, atoms);

  int non_hydrogen_atom = m.natoms() - m.natoms(1);

  double one_bcount_adj = one_bond_count + alpha;
  double two_bcount_adj = two_bond_count + alpha;
  double three_bcount_adj = three_bond_count + alpha;

  double non_h_atom_adj = non_hydrogen_atom + alpha;

  if (one_bond_count != 0) {
    kappa_simple_descriptors[0] = non_hydrogen_atom * (non_hydrogen_atom - 1) *
                                  (non_hydrogen_atom - 1) /
                                  ((double)one_bond_count * one_bond_count);
    kappa_alpha_descriptors[0] = (non_h_atom_adj) * (non_h_atom_adj - 1) *
                                 (non_h_atom_adj - 1) / (one_bcount_adj * one_bcount_adj);
  } else {
    kappa_simple_descriptors[0] = 1;
    kappa_alpha_descriptors[0] = 1;
  }
  if (two_bond_count > 1) {
    kappa_simple_descriptors[1] = (non_hydrogen_atom - 1) * (non_hydrogen_atom - 2) *
                                  (non_hydrogen_atom - 2) /
                                  ((double)two_bond_count * two_bond_count);
    kappa_alpha_descriptors[1] = (non_h_atom_adj - 1) * (non_h_atom_adj - 2) *
                                 (non_h_atom_adj - 2) / (two_bcount_adj * two_bcount_adj);
  } else {
    kappa_simple_descriptors[1] = non_hydrogen_atom - 1;
    kappa_alpha_descriptors[1] = non_hydrogen_atom - 1;
  }
  if (three_bond_count != 0) {
    if (non_hydrogen_atom % 2) {
      kappa_simple_descriptors[2] = (non_hydrogen_atom - 1) * (non_hydrogen_atom - 3) *
                                    (non_hydrogen_atom - 3) /
                                    ((double)three_bond_count * three_bond_count);
      kappa_alpha_descriptors[2] = (non_h_atom_adj - 1) * (non_h_atom_adj - 3) *
                                   (non_h_atom_adj - 3) /
                                   (three_bcount_adj * three_bcount_adj);
    } else {
      kappa_simple_descriptors[2] = (non_hydrogen_atom - 3) * (non_hydrogen_atom - 2) *
                                    (non_hydrogen_atom - 2) /
                                    ((double)three_bond_count * three_bond_count);
      kappa_alpha_descriptors[2] = (non_h_atom_adj - 3) * (non_h_atom_adj - 2) *
                                   (non_h_atom_adj - 2) /
                                   (three_bcount_adj * three_bcount_adj);
    }
  } else {
    switch (non_hydrogen_atom) {
      case 2:
        kappa_simple_descriptors[2] = 1.450;
        kappa_alpha_descriptors[2] = 1.450;
        break;
      case 3:
        kappa_simple_descriptors[2] = 2.000;
        kappa_alpha_descriptors[2] = 2.000;
        break;
      case 4:
        kappa_simple_descriptors[2] = 3.378;
        kappa_alpha_descriptors[2] = 3.378;
        break;
      default:
        kappa_simple_descriptors[2] = 0;
        kappa_alpha_descriptors[2] = 0;
        break;
    }
  }

  //  cerr<<"kappa 0="<<kappa_simple_descriptors[0]<<"\tkappa
  //  1="<<kappa_simple_descriptors[1]<<"\tkappa 2="<<kappa_simple_descriptors[2]<<'\n';
  //  cerr<<"kappa 0="<<kappa_alpha_descriptors[0]<<"\tkappa
  //  1="<<kappa_alpha_descriptors[1]<<"\tkappa 3="<<kappa_alpha_descriptors[3]<<'\n';
}

static void
compute_simple_delta_value(const Molecule& m, double simple_delta[],
                           const Atom* const* atoms) {
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    simple_delta[i] = atoms[i]->ncon();  // - m.explicit_hydrogens(i);
  }
}

static void
compute_valence_delta_value(const Molecule& m, double valence_delta[],
                            const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    int atomic_number_i = z[i];

    int valence = atom_valence(z[i]);
    valence_delta[i] = ((double)(valence - hcount[i])) / (atomic_number_i - valence - 1);
  }
}

static void
compute_modified_valence_delta_value(const Molecule& m, double modified_valence_delta[],
                                     const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    int valence = atom_valence(z[i]);
    int pq_number = prime_quantum_number(z[i]) - 1;
    modified_valence_delta[i] =
        ((double)(valence - hcount[i])) / (double)pq_number / (double)pq_number;
  }
}

static double
value_of_zero_order_chi(const Molecule& m, const double delta[]) {
  int n_atoms = m.natoms();
  double value = 0.0;

  for (int i = 0; i < n_atoms; i++) {
    value += 1.0 / sqrt(delta[i]);
  }

  return value;
}


struct ChiDescriptors {
  public:
    double chi_descriptors[15];
    double valence_chi_descriptors[15];
    double difference_modified_v_chi_descriptors[15];

    double chain_chi_descriptors[8];
    double chain_valence_chi_descriptors[8];
    double chain_modified_v_chi_descriptors[8];

    double bt_indexes[5];  // bonchev_trinajstic_information_indexes

    double difference_chi_descriptors[11];
    double difference_valence_chi_descriptors[11];
    double modified_v_chi_descriptors[15];
};

struct DeltaDescriptors {
  double* simple_delta;
  double* valence_delta;
  double* modified_valence_delta;
  double* estate_index;
  double* tp_index;
  double* etp_index;
  double* values;

  DeltaDescriptors(int natoms);
  ~DeltaDescriptors();
};

DeltaDescriptors::DeltaDescriptors(int natoms) {
  simple_delta = new double[natoms];
  valence_delta = new double[natoms];
  modified_valence_delta = new double[natoms];
  estate_index = new double[natoms];
  tp_index = new double[natoms];
  etp_index = new double[natoms];
  values = new double[natoms];
}

struct MoleculeData {
  atomic_number_t* atomic_number;
  const Atom** atoms;
  int* hcount;

  MoleculeData(Molecule& m);
  ~MoleculeData();
};

MoleculeData::MoleculeData(Molecule& m) {
  const int matoms = m.natoms();
  atomic_number = new atomic_number_t[matoms];
  m.atomic_numbers(atomic_number);

  atoms = new const Atom*[matoms];
  m.atoms(atoms);

  hcount = new int[matoms];
  for (int i = 0; i < matoms; ++i) {
    hcount[i] = m.hcount(i);
  }
}

MoleculeData::~MoleculeData() {
  delete [] atomic_number;
  delete [] atoms;
  delete [] hcount;
}

DeltaDescriptors::~DeltaDescriptors() {
  delete [] simple_delta;
  delete [] valence_delta;
  delete [] modified_valence_delta;
  delete [] estate_index;
  delete [] tp_index;
  delete [] etp_index;
  delete [] values;
}

static void
higher_order_chi_descriptor_computation_procedure(
    Molecule& m, double* simple_delta, double* valence_delta,
    double* modified_valence_delta, const double* estate_index,
    resizable_array_p<Path_with_values>& paths, double& chi_descriptors,
    double& valence_chi_descriptors, double& modified_v_chi_descriptors,
    double* topological_index, double* etopological_index, double& chain_chi,
    double& chain_v_chi, double& chain_m_v_chi, int& n_rings, const Atom* const* atoms,
    const atomic_number_t* z, const int* hcount) {
  double values[NUMBER_OF_VALUES_FOR_PATH];

  n_rings = 0;
  chain_chi = 0.0;
  chain_v_chi = 0.0;
  chain_m_v_chi = 0.0;

  chi_descriptors = 0.0;
  valence_chi_descriptors = 0.0;
  modified_v_chi_descriptors = 0.0;

  // Scope these here for efficiency

  double temp_simple_chi_value;
  double temp_valence_chi_value;
  double temp_tp_value;
  double temp_etp_value;
  double temp_modified_v_chi_value;
  double v1, v2, v3;
  double tp_value, etp_value;

  int n_path = paths.number_elements();

  resizable_array_p<Path_with_values> temp_path;
  temp_path.reserve(3 * n_path);

  double length = 1.0;
  for (int j = 0; j < n_path; j++) {
    const Path_with_values* path = paths[j];

    atom_number_t a1 = path->a1();

    const Atom* atoma1 = atoms[a1];
    int ncona1 = atoma1->ncon();
    if (ncona1 < 2) {
      continue;
    }

    atom_number_t a2 = path->a2();

    const Atom* atoma2 = atoms[a2];
    int ncona2 = atoma2->ncon();
    if (ncona2 < 2) {
      continue;
    }

    length = path->length();

    const double* v = path->raw_values();

    for (int k = 0; k < ncona1; k++) {
      atom_number_t a1_next = atoma1->other(a1, k);

      if (path->contains(a1_next)) {
        continue;
      }

      for (int l = 0; l < ncona2; l++) {
        atom_number_t a2_next = atoma2->other(a2, l);

        if (path->contains(a2_next)) {
          continue;
        }

        temp_simple_chi_value = v[0];
        temp_valence_chi_value = v[1];
        temp_tp_value = v[2];
        temp_etp_value = v[3];
        temp_modified_v_chi_value = v[4];

        if (a2_next == a1_next) {
          n_rings++;
          chain_chi += temp_simple_chi_value / sqrt(simple_delta[a1_next]);
          chain_v_chi += temp_valence_chi_value / sqrt(valence_delta[a1_next]);
          chain_m_v_chi +=
              temp_modified_v_chi_value / sqrt(modified_valence_delta[a1_next]);
          continue;
        }

        v1 =
            temp_simple_chi_value / (sqrt(simple_delta[a1_next] * simple_delta[a2_next]));
        v2 = temp_valence_chi_value /
             (sqrt(valence_delta[a1_next] * valence_delta[a2_next]));
        v3 = temp_modified_v_chi_value /
             (sqrt(modified_valence_delta[a1_next] * modified_valence_delta[a2_next]));

        tp_value =
            exp(length / (length + 2.0) * log(temp_tp_value) +
                log(valence_delta[a1_next] * valence_delta[a2_next]) / (length + 2.0));
        if (temp_etp_value > 0.00 && estate_index[a1_next] > 0.0 &&
            estate_index[a2_next] > 0.0) {
          etp_value =
              exp(length / (length + 2.0) * log(temp_etp_value) +
                  log(estate_index[a1_next] * estate_index[a2_next]) / (length + 2.0));
        } else {
          etp_value = 0.0;
        }

        values[0] = v1;
        values[1] = v2;
        values[2] = tp_value;
        values[3] = etp_value;
        values[4] = v3;

        Path_with_values* new_path =
            new Path_with_values(*path, a1_next, a2_next, values);

        temp_path.add(new_path);

        topological_index[a1_next] += tp_value / (length + 2.0) / (length + 2.0);
        topological_index[a2_next] += tp_value / (length + 2.0) / (length + 2.0);

        etopological_index[a1_next] += etp_value / (length + 2.0) / (length + 2.0);
        etopological_index[a2_next] += etp_value / (length + 2.0) / (length + 2.0);

        chi_descriptors += v1;
        valence_chi_descriptors += v2;
        modified_v_chi_descriptors += v3;
      }
    }
  }

  double temp_length = length + 1.0;

  n_rings /= static_cast<int>(temp_length);

  chain_chi /= (temp_length);
  chain_v_chi /= (temp_length);
  chain_m_v_chi /= (temp_length);

  paths.resize_keep_storage(0);

  paths.transfer_in(temp_path);

  return;
}

static void
compute_path_cluster_descriptors(Molecule m, double& c3, double& c4, double& pc4,
                                 double& c3_pc4, double* delta, double& c3v, double& c4v,
                                 double& pc4v, double& c3v_pc4v, double* valence_delta,
                                 double& mvc3, double& mvc4, double& mvpc4,
                                 double& mvc3_mvpc4, double* modified_valence_delta,
                                 int& c3_count, int& c4_count, int& pc4_count,
                                 const Atom* const* atoms, const atomic_number_t* z,
                                 const int* hcount) {
  c3 = 0.0;
  c4 = 0.0;
  pc4 = 0.0;
  c3_pc4 = 0.0;

  c3v = 0.0;
  c4v = 0.0;
  pc4v = 0.0;
  c3v_pc4v = 0.0;

  mvc3 = 0.0;
  mvc4 = 0.0;
  mvpc4 = 0.0;
  mvc3_mvpc4 = 0.0;

  c3_count = 0;
  c4_count = 0;
  pc4_count = 0;

  int n_atoms = m.natoms();
  atom_number_t cluster_end_atom[3];

  for (int i = 0; i < n_atoms; i++) {
    const Atom* atomi = atoms[i];

    int nconi = atomi->ncon();

    // if this atom only connects with hydrogen atoms, continue;
    if (nconi < 3) {
      continue;
    }

    for (int j = 0; j < nconi; j++) {
      cluster_end_atom[0] = atomi->other(i, j);

      for (int k = j + 1; k < nconi; k++) {
        cluster_end_atom[1] = atomi->other(i, k);

        for (int l = k + 1; l < nconi; l++) {
          cluster_end_atom[2] = atomi->other(i, l);

          double component =
              1.0 / sqrt(delta[i] * delta[cluster_end_atom[0]] *
                         delta[cluster_end_atom[1]] * delta[cluster_end_atom[2]]);
          double component_v =
              1.0 / sqrt(valence_delta[i] * valence_delta[cluster_end_atom[0]] *
                         valence_delta[cluster_end_atom[1]] *
                         valence_delta[cluster_end_atom[2]]);
          double component_mv = 1.0 / sqrt(modified_valence_delta[i] *
                                           modified_valence_delta[cluster_end_atom[0]] *
                                           modified_valence_delta[cluster_end_atom[1]] *
                                           modified_valence_delta[cluster_end_atom[2]]);

          c3 += component;
          c3v += component_v;
          mvc3 += component_mv;

          c3_count++;

          for (int ll = l + 1; ll < nconi; ll++) {
            atom_number_t atom_ll = atomi->other(i, ll);

            c4 += component * 1. / sqrt(delta[atom_ll]);
            c4v += component_v * 1. / sqrt(valence_delta[atom_ll]);
            mvc4 += component_mv * 1. / sqrt(modified_valence_delta[atom_ll]);

            c4_count++;
          }

          for (int ll = 0; ll < 3; ll++) {
            const Atom* cluster_end_atom_ll = atoms[cluster_end_atom[ll]];

            int nconll = cluster_end_atom_ll->ncon();

            //  if (nconll-m.explicit_hydrogens(cluster_end_atom[ll]) < 2) continue;

            if (nconll < 2) {
              continue;
            }

            for (int p = 0; p < nconll; p++) {
              atom_number_t atom_p = cluster_end_atom_ll->other(cluster_end_atom[ll], p);

              if (i == atom_p) {
                continue;
              }

              int not_bound_to_common_atom = 1;

              for (int o = 0; o < 3; o++) {
                if (atom_p == cluster_end_atom[o]) {
                  not_bound_to_common_atom = 0;
                  break;
                }
              }

              if (not_bound_to_common_atom) {
                pc4 += component * 1 / sqrt(delta[atom_p]);
                c3_pc4 += component * 1 / sqrt(delta[atom_p]) - component;

                pc4v += component_v * 1. / sqrt(valence_delta[atom_p]);
                c3v_pc4v += component_v * 1. / sqrt(valence_delta[atom_p]) - component_v;

                mvpc4 += component_mv * 1 / sqrt(modified_valence_delta[atom_p]);
                mvc3_mvpc4 += component_mv * 1 / sqrt(modified_valence_delta[atom_p]) -
                              component_mv;

                pc4_count++;
              }
            }
          }
        }
      }
    }
  }
}

static void
setup_terminal_delta(const Molecule& m, double* terminal_delta, const Atom* const* atoms,
                     const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    int valence = 0;
    int terminal_hydrogen = 1;

    int atomic_number_i = z[i];
    // C
    if (6 == atomic_number_i) {
      valence = 4;
      terminal_hydrogen = 3;
    }
    // F Cl Br I
    else if ((9 == atomic_number_i) || (17 == atomic_number_i) ||
             (35 == atomic_number_i) || (53 == atomic_number_i)) {
      valence = 7;
      terminal_hydrogen = 0;
    }
    // O S
    else if ((8 == atomic_number_i) || (16 == atomic_number_i)) {
      valence = 6;
      terminal_hydrogen = 1;
    }
    // N P
    else if ((7 == atomic_number_i) || (15 == atomic_number_i)) {
      valence = 5;
      terminal_hydrogen = 2;
    }
    // H
    if (1 == atomic_number_i) {
      terminal_delta[i] = 0;
      continue;
    }

    terminal_delta[i] =
        ((double)valence - terminal_hydrogen) / (atomic_number_i - valence - 1);
    //  cerr<<"atom i="<<i<<"\t valance_delta ="<<(valance - hcount[i])<<'\n';
    //      return valance - m.hcount(i);
  }
}

static void
setup_terminal_modified_delta(Molecule& m, double* terminal_modified_delta,
                              const atomic_number_t* z) {
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    int valence = 0;
    int terminal_hydrogen = 1;
    int pq_number = prime_quantum_number(z[i]) - 1;

    int atomic_number_i = z[i];
    // C
    if (6 == atomic_number_i) {
      valence = 4;
      terminal_hydrogen = 3;
    }
    // F Cl Br I
    else if ((9 == atomic_number_i) || (17 == atomic_number_i) ||
             (35 == atomic_number_i) || (53 == atomic_number_i)) {
      valence = 7;
      terminal_hydrogen = 0;
    }
    // O S
    else if ((8 == atomic_number_i) || (16 == atomic_number_i)) {
      valence = 6;
      terminal_hydrogen = 1;
    }
    // N P
    else if ((7 == atomic_number_i) || (15 == atomic_number_i)) {
      valence = 5;
      terminal_hydrogen = 2;
    }
    // H
    if (1 == atomic_number_i) {
      terminal_modified_delta[i] = 0;
      continue;
    }

    terminal_modified_delta[i] =
        ((double)valence - terminal_hydrogen) / (double)pq_number / (double)pq_number;
    //  cerr<<"atom i="<<i<<"\t valance_delta ="<<(valance - m.hcount(i))<<'\n';
    //      return valance - m.hcount(i);
  }
}

static void
setup_linear_path_chi(Molecule& m, double* valence_delta, double* simple_linear_path_chi,
                      double* valence_linear_path_chi,
                      double* modified_valence_linear_path_chi, double* terminal_delta,
                      double* terminal_modified_delta, const Atom* const* atoms,
                      const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();
  int non_hydrogen_atoms = n_atoms;  // - m.natoms(1);

  double sqrt2 = sqrt(2.0);

  set_vector(valence_linear_path_chi, 11, 0.0);
  set_vector(modified_valence_linear_path_chi, 11, 0.0);
  set_vector(simple_linear_path_chi, 11, 0.0);
  if (0 == non_hydrogen_atoms) {
    return;
  }

  // setup simple linear path chi
  simple_linear_path_chi[0] = 2.0 + (non_hydrogen_atoms - 2.0) / sqrt2;
  simple_linear_path_chi[1] = sqrt2 + (non_hydrogen_atoms - 3.) / 2.0;
  simple_linear_path_chi[2] = 1.0 + 0.5 * (non_hydrogen_atoms - 4.0) / sqrt2;
  simple_linear_path_chi[3] = 1.0 / sqrt2 + 0.25 * (non_hydrogen_atoms - 5.0);
  simple_linear_path_chi[4] = 0.5 + (non_hydrogen_atoms - 6.) * 0.25 / sqrt2;
  simple_linear_path_chi[5] = 0.5 / sqrt2 + 0.125 * (non_hydrogen_atoms - 7.0);
  simple_linear_path_chi[6] = 0.25 + 0.125 * (non_hydrogen_atoms - 8.) / sqrt2;
  simple_linear_path_chi[7] = 0.25 / sqrt2 + 0.0625 * (non_hydrogen_atoms - 9.);
  simple_linear_path_chi[8] = 0.125 + 0.0625 * (non_hydrogen_atoms - 10.) / sqrt2;
  simple_linear_path_chi[9] = 0.0625 / sqrt2 + 0.03125 * (non_hydrogen_atoms - 11.);
  simple_linear_path_chi[10] = 0.03125 + 0.03125 * (non_hydrogen_atoms - 12.) / sqrt2;

  // setup valence linear path chi

  double temp_chi[11];

  temp_chi[0] = 1.0;

  for (int i = 1; i < 11; i++) {
    temp_chi[i] = temp_chi[i - 1] / sqrt(2.0);
  }

  setup_terminal_delta(m, terminal_delta, atoms, z, hcount);
  setup_terminal_modified_delta(m, terminal_modified_delta, z);

  double valence_delta_sum = 0;
  double modified_valence_delta_sum = 0;

  for (int i = 0; i < n_atoms; i++) {
    valence_delta_sum += terminal_delta[i];
    modified_valence_delta_sum += terminal_modified_delta[i];

    for (int j = 0; j < 11; j++) {
      valence_linear_path_chi[j] += 1.0 / sqrt(terminal_delta[i]) * temp_chi[j];
      modified_valence_linear_path_chi[j] +=
          1.0 / sqrt(terminal_modified_delta[i]) * temp_chi[j];
    }
  }

  for (int i = 0; i < 11; i++) {
    double temp = (sqrt(valence_delta_sum / (double)non_hydrogen_atoms) *
                       (non_hydrogen_atoms - i - 2) / sqrt2 +
                   2) /
                  (double)non_hydrogen_atoms;
    double mv_temp = (sqrt(modified_valence_delta_sum / (double)non_hydrogen_atoms) *
                          (non_hydrogen_atoms - i - 2) / sqrt2 +
                      2) /
                     (double)non_hydrogen_atoms;

    valence_linear_path_chi[i] = valence_linear_path_chi[i] * temp;
    modified_valence_linear_path_chi[i] = modified_valence_linear_path_chi[i] * mv_temp;
  }

  for (int i = non_hydrogen_atoms; i < 11; i++) {
    valence_linear_path_chi[i] = 0;
    modified_valence_linear_path_chi[i] = 0;
  }

  if (non_hydrogen_atoms < 12) {
    for (int i = non_hydrogen_atoms; i < 11; i++) {
      simple_linear_path_chi[i] = 0;
      valence_linear_path_chi[i] = 0;
      modified_valence_linear_path_chi[i] = 0;
    }

    if (non_hydrogen_atoms > 1) {
      switch (non_hydrogen_atoms - 2) {
        case 0:
          simple_linear_path_chi[1] = 1.0;
          break;
        case 1:
          simple_linear_path_chi[2] = 1. / sqrt2;
          break;
        case 2:
          simple_linear_path_chi[3] = 0.5;
          break;
        case 3:
          simple_linear_path_chi[4] = 0.5 / sqrt2;
          break;
        case 4:
          simple_linear_path_chi[5] = 0.25;
          break;
        case 5:
          simple_linear_path_chi[6] = 0.25 / sqrt2;
        case 6:
          simple_linear_path_chi[7] = 0.125;
          break;
        case 7:
          simple_linear_path_chi[8] = 0.125 / sqrt2;
          break;
        case 8:
          simple_linear_path_chi[9] = 0.0625;
          break;
        case 9:
          simple_linear_path_chi[10] = 0.0625 / sqrt2;
          break;
      }
    }
  }
}

static void
setup_linear_path_chi(Molecule& m, double* valence_delta, double* simple_linear_path_chi,
                      double* valence_linear_path_chi,
                      double* modified_valence_linear_path_chi, const Atom* const* atoms,
                      const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();

  double* terminal_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_terminal_delta(terminal_delta);
  double* terminal_modified_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_terminal_modified_delta(terminal_modified_delta);

  setup_linear_path_chi(m, valence_delta, simple_linear_path_chi, valence_linear_path_chi,
                        modified_valence_linear_path_chi, terminal_delta,
                        terminal_modified_delta, atoms, z, hcount);
}

static void
compute_difference_chi_descriptors(double* difference_chi_descriptors,
                                   double* chi_descriptors, double* linear_path_chi) {
  for (int i = 0; i < 11; i++) {
    difference_chi_descriptors[i] = chi_descriptors[i] - linear_path_chi[i];
  }
}

static void
compute_terminal_methyl_group(Molecule& m, resizable_array_p<Path_with_values>& paths,
                              int& terminal_methyl_count, const Atom* const* atoms,
                              const atomic_number_t* z, const int* hcount) {
  int n_path = paths.number_elements();
  for (int j = 0; j < n_path; j++) {
    Path_with_values* path = paths.item(j);
    atom_number_t a1 = path->a1();
    atom_number_t a2 = path->a2();

    const Atom* atoma1 = atoms[a1];
    const Atom* atoma2 = atoms[a2];

    int ncona1 = atoma1->ncon();
    int ncona2 = atoma2->ncon();

    //      if ((ncona1 - m.explicit_hydrogens(a1) < 2) || (ncona2 -
    //      m.explicit_hydrogens(a2) < 2)) continue;
    if ((ncona1 < 2) || (ncona2 < 2)) {
      continue;
    }

    for (int k = 0; k < ncona1; k++) {
      atom_number_t a1_next = atoma1->other(a1, k);

      if ((1 == z[a1_next]) || (path->contains(a1_next))) {
        continue;
      }

      const Atom* atoma1_next = atoms[a1_next];

      for (int l = 0; l < ncona2; l++) {
        atom_number_t a2_next = atoma2->other(a2, l);
        if (a2_next == a1_next || path->contains(a2_next)) {
          continue;
        }

        const Atom* atoma2_next = atoms[a2_next];

        if ((6 == atoma1_next->atomic_number()) && (3 == hcount[a1_next]) &&
            (atoma1_next->ncon() == atoma1_next->nbonds())) {
          if ((6 == atoma2_next->atomic_number()) && (3 == hcount[a2_next]) &&
              (atoma2_next->ncon() == atoma2_next->nbonds())) {
            terminal_methyl_count++;
          }
        }
      }
    }
  }
}

static void
compute_chi_descriptors(
    Molecule& m, double simple_delta[], double valence_delta[],
    double modified_valence_delta[], double estate_index[], double chi_descriptors[],
    double valence_chi_descriptors[], double modified_v_chi_descriptors[],
    int& total_wiener_number, double chain_chi_descriptors[],
    double chain_valence_chi_descriptors[], double chain_modified_v_chi_descriptors[],
    double difference_chi_descriptors[], double difference_valence_chi_descriptors[],
    double difference_modified_v_chi_descriptors[], int bond_count[],
    double& total_topological_index, double& total_e_topological_index,
    double topological_index[], double e_topological_index[], int& n_rings,
    double* values, const Atom* const* atoms, const atomic_number_t* z,
    const int* hcount) {
  int n_atoms = m.natoms();

  set_vector(chi_descriptors, 14, 0.);
  set_vector(valence_chi_descriptors, 14, 0.);
  set_vector(modified_v_chi_descriptors, 15, 0.);

  set_vector(chain_chi_descriptors, 8, 0.);
  set_vector(chain_valence_chi_descriptors, 8, 0.);
  set_vector(chain_modified_v_chi_descriptors, 8, 0.);

  set_vector(bond_count, 15, 0);

  set_vector(topological_index, n_atoms, 0.0);
  set_vector(e_topological_index, n_atoms, 0.0);

  total_wiener_number = 0;

  // compute the delta_value and valence delta value
  compute_simple_delta_value(m, simple_delta, atoms);
  compute_valence_delta_value(m, valence_delta, z, hcount);
  compute_modified_valence_delta_value(m, modified_valence_delta, z, hcount);

  // compute the estate index
  (void)determine_atom_e_state_index(m, estate_index, atoms, z, hcount);

  for (int i = 0; i < n_atoms; i++) {
    estate_index[i] = fabs(estate_index[i]);
    //  cerr<<"atom i="<<i<<"\t"<<m.atomic_symbol(i)<<"t
    //  estate_index["<<i<<"]="<<estate_index[i]<<'\n';
  }

  // take care of the zero order chi values
  chi_descriptors[0] = value_of_zero_order_chi(m, simple_delta);
  valence_chi_descriptors[0] = value_of_zero_order_chi(m, valence_delta);
  modified_v_chi_descriptors[0] = value_of_zero_order_chi(m, modified_valence_delta);

  resizable_array_p<Path_with_values> odd_path;
  resizable_array_p<Path_with_values> even_path;

  odd_path.resize(2000);
  even_path.resize(2000);

  // set up odd number bond path
  int n_bonds = m.nedges();

  chi_descriptors[1] = 0;
  valence_chi_descriptors[1] = 0;
  modified_v_chi_descriptors[1] = 0;

  for (int i = 0; i < n_bonds; i++) {
    const Bond* bondi = m.bondi(i);

    atom_number_t a1 = bondi->a1();
    atom_number_t a2 = bondi->a2();

    double sd1 = simple_delta[a1];
    double sd2 = simple_delta[a2];
    double vd1 = valence_delta[a1];
    double vd2 = valence_delta[a2];
    double mvd1 = modified_valence_delta[a1];
    double mvd2 = modified_valence_delta[a2];

    double temp_chi = 1.0 / sqrt(sd1 * sd2);
    double temp_valence_chi = 1.0 / sqrt(vd1 * vd2);
    double temp_modified_valence_chi = 1.0 / sqrt(mvd1 * mvd2);

    double tp_value = sqrt(vd1 * vd2);
    double etp_value = sqrt(estate_index[a1] * estate_index[a2]);

    values[0] = temp_chi;
    values[1] = temp_valence_chi;
    values[2] = tp_value;
    values[3] = etp_value;
    values[4] = temp_modified_valence_chi;

    odd_path.add(new Path_with_values(a1, a2, values, n_atoms));

    topological_index[a1] += tp_value / 4.0;
    topological_index[a2] += tp_value / 4.0;

    e_topological_index[a1] += etp_value / 4.0;
    e_topological_index[a2] += etp_value / 4.0;

#ifdef DEBUG_HOCD
    cerr << "*e_topological_index updated " << e_topological_index[a1] << " and "
         << e_topological_index[a2] << ", etp_value " << etp_value << '\n';
#endif

    chi_descriptors[1] += temp_chi;
    valence_chi_descriptors[1] += temp_valence_chi;
    modified_v_chi_descriptors[1] += temp_modified_valence_chi;
  }

  total_wiener_number += odd_path.number_elements();
  bond_count[1] = odd_path.number_elements();
  bond_count[0] += bond_count[1];

  compute_terminal_methyl_group(m, odd_path, bond_count[14], atoms, z, hcount);

  // set up even number bond path
  int n_non_hydrogen_atom = n_atoms;  // - m.natoms(1);

  chi_descriptors[2] = 0;
  valence_chi_descriptors[2] = 0;
  modified_v_chi_descriptors[2] = 0;

  for (int i = 0; i < n_atoms; i++) {
    const Atom* ai = atoms[i];

    int nconi = ai->ncon();

    //      if (nconi-m.explicit_hydrogens(i)<2) continue;

    if (nconi < 2) {
      continue;
    }

    for (int j = 0; j < nconi; j++) {
      atom_number_t atomj = ai->other(i, j);

      for (int k = j + 1; k < nconi; k++) {
        atom_number_t atomk = ai->other(i, k);

        double sd1 = simple_delta[i];
        double sd2 = simple_delta[atomj];
        double sd3 = simple_delta[atomk];
        double vd1 = valence_delta[i];
        double vd2 = valence_delta[atomj];
        double vd3 = valence_delta[atomk];
        double mvd1 = modified_valence_delta[i];
        double mvd2 = modified_valence_delta[atomj];
        double mvd3 = modified_valence_delta[atomk];

        double temp_chi = 1.0 / sqrt(sd1 * sd2 * sd3);
        double temp_valence_chi = 1.0 / sqrt(vd1 * vd2 * vd3);
        double temp_modified_valence_chi = 1.0 / sqrt(mvd1 * mvd2 * mvd3);

        double temp = vd1 * vd2 * vd3;
        double etemp = estate_index[i] * estate_index[atomj] * estate_index[atomk];

        double tp_value;
        if (temp > 0.0) {
          tp_value = exp(1.0 / 3.0 * log(temp));
        } else {
          tp_value = 0.0;
        }

        double etp_value;
        if (etemp > 0.0) {
          etp_value = exp(1.0 / 3.0 * log(etemp));
        } else {
          etp_value = 0.0;
        }

        chi_descriptors[2] += temp_chi;
        valence_chi_descriptors[2] += temp_valence_chi;
        modified_v_chi_descriptors[2] += temp_modified_valence_chi;

        values[0] = temp_chi;
        values[1] = temp_valence_chi;
        values[2] = tp_value;
        values[3] = etp_value;
        values[4] = temp_modified_valence_chi;

        even_path.add(new Path_with_values(atomj, atomk, i, values, n_atoms));

        topological_index[atomj] += tp_value / 9.0;
        topological_index[atomk] += tp_value / 9.0;

#ifdef DEBUG_HOCD
        cerr << " updating e_topo " << e_topological_index[atomj] << " and "
             << e_topological_index[atomk] << " etp_value " << etp_value << '\n';
#endif
        e_topological_index[atomj] += etp_value / 9.0;
        e_topological_index[atomk] += etp_value / 9.0;
      }
    }
  }

  total_wiener_number += even_path.number_elements() * 2;
  bond_count[2] = even_path.number_elements();
  bond_count[0] += bond_count[2];

  if (odd_path.number_elements() != 0) {
    for (int i = 0; i < n_non_hydrogen_atom / 2; i++) {
      double temp_chi;
      double temp_v_chi;
      double temp_m_v_chi;

      double temp_chain_chi;
      double temp_chain_v_chi;
      double temp_chain_m_v_chi;

      int temp_n_rings;

      higher_order_chi_descriptor_computation_procedure(
          m, simple_delta, valence_delta, modified_valence_delta, estate_index, odd_path,
          temp_chi, temp_v_chi, temp_m_v_chi, topological_index, e_topological_index,
          temp_chain_chi, temp_chain_v_chi, temp_chain_m_v_chi, temp_n_rings, atoms, z,
          hcount);

      n_rings += temp_n_rings;
      bond_count[0] += odd_path.number_elements();

      if (i < 4) {
        chi_descriptors[i * 2 + 3] = temp_chi;
        valence_chi_descriptors[i * 2 + 3] = temp_v_chi;
        modified_v_chi_descriptors[i * 2 + 3] = temp_m_v_chi;

        chain_chi_descriptors[i * 2] = temp_chain_chi;
        chain_valence_chi_descriptors[i * 2] = temp_chain_v_chi;
        chain_modified_v_chi_descriptors[i * 2] = temp_chain_m_v_chi;
      }

      total_wiener_number += odd_path.number_elements() * (i * 2 + 3);

      if (i * 2 + 3 < 11) {
        bond_count[i * 2 + 3] = odd_path.number_elements();
      }

      if (odd_path.number_elements() == 0) {
        break;
      }
      if (odd_path.number_elements() > MAX_PATH_NUMBER) {
        break;
      }
    }
  }

  odd_path.resize(0);

  if (even_path.number_elements() != 0) {
    for (int i = 0; i < n_non_hydrogen_atom / 2; i++) {
      double temp_chi;
      double temp_v_chi;
      double temp_m_v_chi;

      double temp_chain_chi;
      double temp_chain_v_chi;
      double temp_chain_m_v_chi;

      int temp_n_rings;

      higher_order_chi_descriptor_computation_procedure(
          m, simple_delta, valence_delta, modified_valence_delta, estate_index, even_path,
          temp_chi, temp_v_chi, temp_m_v_chi, topological_index, e_topological_index,
          temp_chain_chi, temp_chain_v_chi, temp_chain_m_v_chi, temp_n_rings, atoms, z,
          hcount);

      n_rings += temp_n_rings;

      bond_count[0] += even_path.number_elements();

      if (i < 4) {
        chi_descriptors[i * 2 + 4] = temp_chi;
        valence_chi_descriptors[i * 2 + 4] = temp_v_chi;
        modified_v_chi_descriptors[i * 2 + 4] = temp_m_v_chi;

        chain_chi_descriptors[i * 2 + 1] = temp_chain_chi;
        chain_valence_chi_descriptors[i * 2 + 1] = temp_chain_v_chi;
        chain_modified_v_chi_descriptors[i * 2 + 1] = temp_chain_m_v_chi;
      }

      total_wiener_number += even_path.number_elements() * (i * 2 + 4);

      if (i * 2 + 4 < 11) {
        bond_count[i * 2 + 4] = even_path.number_elements();
      }

      if (even_path.number_elements() == 0) {
        break;
      }
      if (even_path.number_elements() > MAX_PATH_NUMBER) {
        break;
      }
    }
  }

  // release the memory
  even_path.resize(0);

  total_topological_index = 0.0;
  total_e_topological_index = 0.0;

  for (int i = 0; i < n_atoms; i++) {
    total_topological_index += topological_index[i] / 2 + valence_delta[i];
    //    cerr << " i = " << i << " e_topological_index " << e_topological_index[i] <<
    //    '\n';

    total_e_topological_index += e_topological_index[i] / 2;
  }

  compute_path_cluster_descriptors(
      m, chi_descriptors[11], chi_descriptors[12], chi_descriptors[13],
      chi_descriptors[14], simple_delta, valence_chi_descriptors[11],
      valence_chi_descriptors[12], valence_chi_descriptors[13],
      valence_chi_descriptors[14], valence_delta, modified_v_chi_descriptors[11],
      modified_v_chi_descriptors[12], modified_v_chi_descriptors[13],
      modified_v_chi_descriptors[14], modified_valence_delta, bond_count[11],
      bond_count[12], bond_count[13], atoms, z, hcount);

  //  compute_ring_chi_descriptors(m, chain_chi_descriptors, simple_delta,
  //  chain_valence_chi_descriptors, valence_delta);

  setup_linear_path_chi(m, valence_delta, simple_linear_path_chi, valence_linear_path_chi,
                        modified_valence_linear_path_chi, atoms, z, hcount);

  compute_difference_chi_descriptors(difference_chi_descriptors, chi_descriptors,
                                     simple_linear_path_chi);
  compute_difference_chi_descriptors(difference_valence_chi_descriptors,
                                     valence_chi_descriptors, valence_linear_path_chi);
  compute_difference_chi_descriptors(difference_modified_v_chi_descriptors,
                                     modified_v_chi_descriptors,
                                     modified_valence_linear_path_chi);
}

static void
compute_chi_descriptors(
    Molecule& m, double chi_descriptors[], double valence_chi_descriptors[],
    double modified_v_chi_descriptors[], int& total_wiener_number,
    double* chain_chi_descriptors, double* chain_valence_chi_descriptors,
    double* chain_modified_v_chi_descriptors, double* difference_chi_descriptors,
    double* difference_valence_chi_descriptors,
    double* difference_modified_v_chi_descriptors, int* bond_count,
    double& total_topological_index, double& total_e_topological_index, int& n_rings,
    const Atom* const* atoms, const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();

  double* simple_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_simple_delta(simple_delta);
  double* valence_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_valence_delta(valence_delta);
  double* modified_valence_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_modified_valence_delta(modified_valence_delta);
  double* estate_index = new double[n_atoms];
  std::unique_ptr<double[]> free_estate_index(estate_index);
  double* tp_index = new double[n_atoms];
  std::unique_ptr<double[]> free_tp_index(tp_index);
  double* etp_index = new double[n_atoms];
  std::unique_ptr<double[]> free_etp_index(etp_index);
  double* values = new double[NUMBER_OF_VALUES_FOR_PATH];
  std::unique_ptr<double[]> free_values(values);

  compute_chi_descriptors(
      m, simple_delta, valence_delta, modified_valence_delta, estate_index,
      chi_descriptors, valence_chi_descriptors, modified_v_chi_descriptors,
      total_wiener_number, chain_chi_descriptors, chain_valence_chi_descriptors,
      chain_modified_v_chi_descriptors, difference_chi_descriptors,
      difference_valence_chi_descriptors, difference_modified_v_chi_descriptors,
      bond_count, total_topological_index, total_e_topological_index, tp_index, etp_index,
      n_rings, values, atoms, z, hcount);

  return;
}

static void
compute_Bonchev_Trinajstic_informaiton_indexes_and_wiener_numbers(
    Molecule& m, double* bt_indexes, int* wiener_number, int* bond_distance_count,
    const Atom* const* atoms, const atomic_number_t* z) {
  int n_atoms = m.natoms();

  // initialize the distance count array and the default bt_indexes
  set_vector(wiener_number, 4, 0);
  set_vector(bond_distance_count, n_atoms, 0);
  set_vector(bt_indexes, 5, 0.0);

  // set up the distance count array
  for (int i = 0; i < n_atoms; i++) {
    for (int j = i + 1; j < n_atoms; j++) {
      int dij = m.bonds_between(i, j);

      wiener_number[0] += dij;
      bond_distance_count[dij]++;
    }
  }

  double sumw = 0.0;
  double sum = 0.0;

  for (int i = 1; i < n_atoms; i++) {
    if (0 == bond_distance_count[i]) {
      break;
    }

    double count = bond_distance_count[i];
    double distance = i;

    bt_indexes[4]++;
    sum += count * log10(count) / log10(2.0);
    sumw += count * distance * log10(distance) / log10(2.0);
  }

  // cerr << "Computing 0, wiener_number " << wiener_number[0] << ", sumw " << sumw <<
  // '\n';

  bt_indexes[0] =
      (double)wiener_number[0] * log10((double)wiener_number[0]) / log10(2.0) - sumw;
  if (0 != wiener_number[0]) {
    bt_indexes[1] = bt_indexes[0] / static_cast<double>(wiener_number[0]);
  }

  if (bt_indexes[0] > 0.0) {
    bt_indexes[0] =
        log10(bt_indexes[0]);  // iaw modification to keep from being such large numbers
  }

  //  double temp = (n_atoms - m.natoms(1)) * (n_atoms - m.natoms(1) -1) / 2.0;
  double temp = static_cast<double>((n_atoms) * (n_atoms - 1) / 2);

  bt_indexes[2] = temp * log10(temp) / log10(2.0) - sum;
  if (0.0 != temp) {
    bt_indexes[3] = bt_indexes[2] / temp;
  }

  if (0.0 != bt_indexes[2]) {
    bt_indexes[2] = log10(bt_indexes[2]);  // iaw modification
  }

  wiener_number[1] = number_of_three_bond_count(m, atoms);

  wiener_number[2] = number_of_two_bond_count(m, atoms) * 2;

#ifdef ECHOWBT
  for (int i = 0; i < 4; i++) {
    cerr << "Wiener Number[" << i << "]=" << wiener_number[i] << '\n';
  }
  for (int i = 0; i < 5; i++) {
    cerr << "BT_index[" << i << "]=" << bt_indexes[i] << '\n';
  }
#endif
}

static void
compute_Bonchev_Trinajstic_informaiton_indexes_and_wiener_numbers(
    Molecule& m, double* bt_indexes, int* wiener_number, const Atom* const* atoms,
    const atomic_number_t* z) {
  int n_atoms = m.natoms();
  int* bond_distance_count = new int[n_atoms];
  std::unique_ptr<int[]> free_bond_distance_count(bond_distance_count);

  compute_Bonchev_Trinajstic_informaiton_indexes_and_wiener_numbers(
      m, bt_indexes, wiener_number, bond_distance_count, atoms, z);

  return;
}

static void
determine_number_of_hydrogen_bond_donor_acceptor(Molecule& m, int& h_bond_donor,
                                                 int& h_bond_acceptor,
                                                 const Atom* const* atoms,
                                                 const atomic_number_t* z,
                                                 const int* hcount) {
  int n_atoms = m.natoms();

  for (int i = 0; i < n_atoms; i++) {
    const Atom* atomi = atoms[i];

    int atomic_number_i = z[i];

    int nconi = atomi->ncon();

    int nbondsi = atomi->nbonds();

    int hcounti = hcount[i];

    //      int e_h_counti = m.explicit_hydrogens(i);

    switch (atomic_number_i) {
      case 7:
        if (2 == hcounti) {
          h_bond_donor += 2;
          h_bond_acceptor++;
        } else if ((1 == hcounti) && (m.is_aromatic(i)) && (nconi == nbondsi)) {
          h_bond_donor++;
        } else if ((1 == hcounti) && (1 == nconi)) {
          h_bond_donor++;
          h_bond_acceptor++;
        } else if ((1 == hcounti) && (2 == nconi)) {
          h_bond_donor++;
          h_bond_acceptor++;
        } else if (4 == nconi && 1 == atomi->formal_charge()) {  // quat not an acceptor
          ;
        } else if ((0 == hcounti) && (nbondsi != 5)) {
          h_bond_acceptor++;
        }
        break;
      case 8:
        if (1 == hcounti)  // && (1 == nconi - e_h_counti))
        {
          h_bond_donor++;
          h_bond_acceptor++;
        } else if (0 == hcounti) {  //&& (2 == nbondsi))
          h_bond_acceptor++;
        }
        break;
      case 16:
        if (1 == hcounti) {  //&& (1 == nconi - e_h_counti))
          h_bond_donor++;
        } else if (6 - hcounti + nconi < 9) {  // if ((0 == hcounti)) && (2 == nbondsi))
          h_bond_acceptor++;
        }
        break;
      case 9:
        if (1 == nbondsi) {
          h_bond_acceptor++;
        }
        break;
      case 17:
        if (1 == nbondsi) {
          h_bond_acceptor++;
        }
        break;
      case 6: {
        int f_cl_count = 0;
        int f_cl_br_count = 0;
        for (int j = 0; j < nconi; j++) {
          const Atom* atomj = atoms[atomi->other(i, j)];

          int atomic_number_j = atomj->atomic_number();
          if ((9 == atomic_number_j) || (17 == atomic_number_j)) {
            f_cl_count++;
            f_cl_br_count++;
          } else if (35 == atomic_number_j) {
            f_cl_br_count++;
          }
        }

        if ((f_cl_br_count > 0) && (1 == hcounti)) {
          h_bond_donor++;
        } else if ((f_cl_br_count > 0) && (1 == hcounti)) {
          h_bond_donor++;
          if (2 == f_cl_count) {
            h_bond_donor++;
          }
        } else if ((1 == f_cl_br_count) && (1 == f_cl_count) && (3 == hcounti)) {
          h_bond_donor++;
        }
      } break;
      default:
        break;
    }
  }
}

static double
value_of_kappa_alpha_0(const Molecule& m, double shannon_ic) {
  return ((double)m.natoms() * shannon_ic);
}

static double
value_of_redundancy(const Molecule& m, double shannon_ic) {
  double redundancy = 0.0;
  int n_hydrogen_atom = m.natoms();  // - m.natoms(1);
  if (n_hydrogen_atom > 1) {
    redundancy = 1. - (shannon_ic / log10((double)n_hydrogen_atom));
  }
  return redundancy;
}

static double
value_of_sum_of_intrinsic_state(Molecule& m, const Atom* const* atoms,
                                const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();

  double sumI = 0;
  for (int i = 0; i < n_atoms; i++) {
    sumI += value_of_Kier_and_Hall_atom_intrinsic_state(m, i, atoms, z, hcount);
  }

  return sumI;
}

static double
value_of_sum_of_delta_intrinsic_state(Molecule& m, const Atom* const* atoms,
                                      const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();

  double sumdelI = 0;
  for (int i = 0; i < n_atoms; i++) {
    for (int j = i + 1; j < n_atoms; j++) {
      double i_value_i =
          value_of_Kier_and_Hall_atom_intrinsic_state(m, i, atoms, z, hcount);
      double i_value_j =
          value_of_Kier_and_Hall_atom_intrinsic_state(m, j, atoms, z, hcount);

      double distance = m.bonds_between(i, j) + 1;
      sumdelI += fabs((i_value_i - i_value_j) / distance / distance);
    }
  }

  return sumdelI;
}

static double
value_of_flexible_index_phia(Molecule& m, double* kappa_alpha_descriptors) {
  int n_hydrogen_atom = m.natoms();  // - m.natoms(1);
  if (n_hydrogen_atom > 0) {
    return (kappa_alpha_descriptors[0] * kappa_alpha_descriptors[1] /
            (double)n_hydrogen_atom);
  } else {
    return 0;
  }
}

static double
value_of_first_order_Zagreb_index(Molecule& m, double* delta) {
  double value = 0;
  int n_atoms = m.natoms();
  for (int i = 0; i < n_atoms; i++) {
    value += delta[i] * delta[i];
  }
  return value;
}

static void
compute_values_of_first_order_Zagreb_index(Molecule& m, double* simple_delta,
                                           double* v_delta, double* mv_delta, double& z1,
                                           double& z1v, double& z1mv,
                                           const Atom* const* atoms,
                                           const atomic_number_t* z, const int* hcount) {
  compute_simple_delta_value(m, simple_delta, atoms);
  compute_valence_delta_value(m, v_delta, z, hcount);
  compute_modified_valence_delta_value(m, mv_delta, z, hcount);

  z1 = value_of_first_order_Zagreb_index(m, simple_delta);
  z1v = value_of_first_order_Zagreb_index(m, v_delta);
  z1mv = value_of_first_order_Zagreb_index(m, mv_delta);
}

static void
compute_values_of_first_order_Zagreb_index(Molecule& m, double& z1, double& z1v,
                                           double& z1mv, const Atom* const* atoms,
                                           const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();
  double* simple_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_simple_delta(simple_delta);
  double* valence_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_valence_delta(valence_delta);
  double* modified_valence_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_modified_valence_delta(modified_valence_delta);

  compute_values_of_first_order_Zagreb_index(m, simple_delta, valence_delta,
                                             modified_valence_delta, z1, z1v, z1mv, atoms,
                                             z, hcount);

  return;
}

static double
value_of_second_order_Zagreb_index(Molecule& m, double* delta) {
  double value = 0;
  int n_bonds = m.nedges();
  for (int i = 0; i < n_bonds; i++) {
    const Bond* bondi = m.bondi(i);
    atom_number_t a1 = bondi->a1();
    atom_number_t a2 = bondi->a2();

    value += delta[a1] * delta[a2];
  }
  return value;
}

static void
compute_values_of_second_order_Zagreb_index(Molecule& m, double* simple_delta,
                                            double* v_delta, double* mv_delta, double& z2,
                                            double& z2v, double& z2mv,
                                            const Atom* const* atoms,
                                            const atomic_number_t* z, const int* hcount) {
  compute_simple_delta_value(m, simple_delta, atoms);
  compute_valence_delta_value(m, v_delta, z, hcount);
  compute_modified_valence_delta_value(m, mv_delta, z, hcount);

  z2 = value_of_second_order_Zagreb_index(m, simple_delta);
  z2v = value_of_second_order_Zagreb_index(m, v_delta);
  z2mv = value_of_second_order_Zagreb_index(m, mv_delta);
}

static void
compute_values_of_second_order_Zagreb_index(Molecule& m, double& z2, double& z2v,
                                            double& z2mv, const Atom* const* atoms,
                                            const atomic_number_t* z, const int* hcount) {
  int n_atoms = m.natoms();
  double* simple_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_simple_delta(simple_delta);
  double* valence_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_valence_delta(valence_delta);
  double* modified_valence_delta = new double[n_atoms];
  std::unique_ptr<double[]> free_modified_valence_delta(modified_valence_delta);

  compute_values_of_second_order_Zagreb_index(m, simple_delta, valence_delta,
                                              modified_valence_delta, z2, z2v, z2mv,
                                              atoms, z, hcount);

  return;
}

/*
 * This function is to calculate the molconn descriptors for molecule
 */

static int
jw_molconn(Molecule& m, IWString_and_File_Descriptor& output, const Atom* const* atoms,
           const atomic_number_t* z, const int* hcount) {
  int n_rings = 0;

  double kappa_simple_descriptors[3];
  double kappa_alpha_descriptors[3];

  (void)compute_kappa_descriptors(m, kappa_simple_descriptors, kappa_alpha_descriptors,
                                  atoms, z, hcount);

  int n_symmetry_class = 0;
  int n_lone_symmetry_class = 0;
  double shannon_ic =
      shannon_information_content(m, n_symmetry_class, n_lone_symmetry_class);

  double ratio_n_atoms_vs_symmetry_class = 1;
  double percent_atoms_lone_symmetry_class = 1;
  if (n_symmetry_class > 0) {
    int non_hydrogen_atoms = m.natoms();  // - m.natoms(1);
    ratio_n_atoms_vs_symmetry_class =
        (double)(non_hydrogen_atoms) / (double)n_symmetry_class;
    percent_atoms_lone_symmetry_class =
        (double)n_lone_symmetry_class / (double)(non_hydrogen_atoms);
  }

  // kappa zero is the shannon information content * number of atoms
  double kappa_alpha_0 = value_of_kappa_alpha_0(m, shannon_ic);
  double redundancy = value_of_redundancy(m, shannon_ic);

  double chi_descriptors[15];
  double valence_chi_descriptors[15];
  double difference_modified_v_chi_descriptors[15];

  double chain_chi_descriptors[8];
  double chain_valence_chi_descriptors[8];
  double chain_modified_v_chi_descriptors[8];

  double bt_indexes[5];  // bonchev_trinajstic_information_indexes

  double difference_chi_descriptors[11];
  double difference_valence_chi_descriptors[11];
  double modified_v_chi_descriptors[15];

  int wiener_number[4];
  (void)compute_Bonchev_Trinajstic_informaiton_indexes_and_wiener_numbers(
      m, bt_indexes, wiener_number, atoms, z);

  double total_topological_index = 0.0;
  double total_e_topological_index = 0.0;
  int bond_distance_count[15];

  (void)compute_chi_descriptors(
      m, chi_descriptors, valence_chi_descriptors, modified_v_chi_descriptors,
      wiener_number[3], chain_chi_descriptors, chain_valence_chi_descriptors,
      chain_modified_v_chi_descriptors, difference_chi_descriptors,
      difference_valence_chi_descriptors, difference_modified_v_chi_descriptors,
      bond_distance_count, total_topological_index, total_e_topological_index, n_rings,
      atoms, z, hcount);

  int h_bond_donor = 0;
  int h_bond_acceptor = 0;

  determine_number_of_hydrogen_bond_donor_acceptor(m, h_bond_donor, h_bond_acceptor,
                                                   atoms, z, hcount);

  double sumI = value_of_sum_of_intrinsic_state(m, atoms, z, hcount);
  double sumdelI = value_of_sum_of_delta_intrinsic_state(m, atoms, z, hcount);

  double flexible_index_phia = value_of_flexible_index_phia(m, kappa_alpha_descriptors);

  double zagreb_1 = 0.0;
  double zagreb_1_v = 0.0;
  double zagreb_1_mv = 0.0;

  double zagreb_2 = 0.0;
  double zagreb_2_v = 0.0;
  double zagreb_2_mv = 0.0;

  compute_values_of_first_order_Zagreb_index(m, zagreb_1, zagreb_1_v, zagreb_1_mv, atoms,
                                             z, hcount);
  compute_values_of_second_order_Zagreb_index(m, zagreb_2, zagreb_2_v, zagreb_2_mv, atoms,
                                              z, hcount);

  if (read_descriptor_file_pipeline && write_descriptor_file_pipeline) {
    output << m.smiles() << ' ';
    output << m.name();  // includes all previously calculated descriptors.
  } else if (read_descriptor_file_pipeline) {
    output << m.name();  // includes all previously calculated descriptors.
  } else if (write_descriptor_file_pipeline) {
    output << m.smiles() << ' ';
    append_first_token_of_name(m.name(), output);
  } else {
    append_first_token_of_name(m.name(), output);
  }

  output << ' ' << n_rings;  // col 2

  for (int i = 0; i < 11; i++) {  // col 3-13
    output << ' ' << static_cast<float>(chi_descriptors[i]);
  }

  for (int i = 0; i < 11; i++) {  // col 14-24
    output << ' ' << static_cast<float>(valence_chi_descriptors[i]);
  }

  for (int i = 11; i < 14; i++) {  // col 25-27
    output << ' ' << static_cast<float>(chi_descriptors[i]);
  }

  for (int i = 0; i < 8; i++)  // col 28-35
  {
    if (fabs(chain_chi_descriptors[i]) < truncate_to_zero) {
      output << " 0";
    } else {
      output << ' ' << static_cast<float>(chain_chi_descriptors[i]);
    }
  }

  for (int i = 11; i < 14; i++)  // col 36-38
  {
    if (fabs(valence_chi_descriptors[i]) < truncate_to_zero) {
      output << " 0";
    } else {
      output << ' ' << static_cast<float>(valence_chi_descriptors[i]);
    }
  }

  for (int i = 0; i < 8; i++)  // col 39-46
  {
    if (fabs(chain_valence_chi_descriptors[i]) < truncate_to_zero) {
      output << " 0";
    } else {
      output << ' ' << static_cast<float>(chain_valence_chi_descriptors[i]);
    }
  }

  for (int i = 0; i < 11; i++)  // col 47-57
  {
    if (fabs(difference_chi_descriptors[i]) < truncate_to_zero) {
      output << " 0";
    } else {
      output << ' ' << static_cast<float>(difference_chi_descriptors[i]);
    }
  }

  for (int i = 0; i < 11; i++)  // col 58-68
  {
    if (fabs(difference_valence_chi_descriptors[i]) < truncate_to_zero) {
      output << " 0";
    } else {
      output << ' ' << static_cast<float>(difference_valence_chi_descriptors[i]);
    }
  }

  output << ' ' << static_cast<float>(kappa_alpha_0);  // col 69

  for (int i = 0; i < 3; i++) {  // col 70-72
    output << ' ' << static_cast<float>(kappa_simple_descriptors[i]);
  }

  for (int i = 0; i < 3; i++) {  // coll 73-75
    output << ' ' << static_cast<float>(kappa_alpha_descriptors[i]);
  }

  output << ' ' << static_cast<float>(shannon_ic) << ' '
         << static_cast<float>(total_topological_index) << ' ' << static_cast<float>(sumI)
         << ' ' << static_cast<float>(sumdelI) << ' '
         << static_cast<float>(total_e_topological_index);  // cols 76-80

  output << ' ' << static_cast<float>(flexible_index_phia);  // col 81

  for (int i = 0; i < 4; i++)  // col 82-85
  {
    output << ' ' << static_cast<float>(bt_indexes[i]);
  }

  output << ' '
         << static_cast<float>(log10(static_cast<double>(
                wiener_number[0])));  // col 86  iaw modification, take log

  for (int i = 1; i < 3; i++)  // col 87-88
  {
    output << ' ' << static_cast<float>(wiener_number[i]);
  }

  output << ' '
         << static_cast<float>(log10(static_cast<double>(
                wiener_number[3])));  // col 89 iawp modification, take log

  output << ' ' << static_cast<float>(chi_descriptors[14]) << ' '
         << static_cast<float>(valence_chi_descriptors[14]);  // cols 90-91

  output << ' ' << static_cast<float>(bt_indexes[4]) << ' ' << h_bond_donor << ' '
         << h_bond_acceptor;  // cols 92-94

  for (int i = 1; i < 14; i++) {  // cols 95-107
    output << ' ' << bond_distance_count[i];
  }

  output << ' '
         << static_cast<float>(log10(static_cast<double>(
                bond_distance_count[0])));   // col 108 iaw modification, take log
  output << ' ' << bond_distance_count[14];  // col 109

  if (compute_additional_connectivity_related_descriptors) {
    if (redundancy > 0.001) {  // col 110
      output << ' ' << static_cast<float>(redundancy);
    } else {
      output << " 0";
    }

    output << ' ' << static_cast<float>(ratio_n_atoms_vs_symmetry_class) << ' '
           << static_cast<float>(percent_atoms_lone_symmetry_class);  // col 111-112

    for (int i = 0; i < 15; i++)  // col 113-127
    {
      if (fabs(modified_v_chi_descriptors[i]) < truncate_to_zero) {
        output << " 0";
      } else {
        output << ' ' << static_cast<float>(modified_v_chi_descriptors[i]);
      }
    }

    for (int i = 0; i < 8; i++)  // col 128-135
    {
      if (fabs(chain_modified_v_chi_descriptors[i]) < truncate_to_zero) {
        output << " 0";
      } else {
        output << ' ' << static_cast<float>(chain_modified_v_chi_descriptors[i]);
      }
    }

    for (int i = 0; i < 11; i++)  // col 136-146
    {
      if (fabs(difference_modified_v_chi_descriptors[i]) < truncate_to_zero) {
        output << " 0";
      } else {
        output << ' ' << static_cast<float>(difference_modified_v_chi_descriptors[i]);
      }
    }

    // first order Zagreb index
    output << ' ' << static_cast<float>(zagreb_1) << ' ' << static_cast<float>(zagreb_1_v)
           << ' ' << static_cast<float>(zagreb_1_mv);  // col 147-149

    // second order Zagreb index
    output << ' ' << static_cast<float>(zagreb_2) << ' ' << static_cast<float>(zagreb_2_v)
           << ' ' << static_cast<float>(zagreb_2_mv);  // col 150-152
  }

  /*
  cerr<<"numhbd ="<<h_bond_donor<<'\n';
  cerr<<"numhba ="<<h_bond_acceptor<<'\n';

  cerr<<"sumI ="<< sumI<<'\n';

  cerr<<"sumdelI ="<<sumdelI<<'\n';

  cerr<<"shannon IC="<<shannon_ic<<'\n';
  cerr<<"Kappa_alpha_0="<<kappa_alpha_0<<'\n';
  cerr<<"Redundancy="<<redundancy<<'\n';

  cerr<<"wiener number ="<<wiener_number[3]<<'\n';
  cerr<<"symmetry class number ="<<n_symmetry_class<<'\n';
  cerr<<"flexible index phia = "<<flexible_index_phia<<'\n';
  */

  output << '\n';

  return 1;
}

static int
jw_molconn(Molecule& m, IWString_and_File_Descriptor& output) {
  m.compute_aromaticity_if_needed();

  int matoms = m.natoms();

  if (0 == matoms) {
    cerr << "Skipping empty structure\n";
    return 1;
  }

  if (matoms < 3) {
    cerr << "Can't do molecule with just " << matoms << " atoms\n";
    return 1;
  }

  const Atom** atoms = new const Atom*[matoms];
  std::unique_ptr<const Atom*[]> free_atoms(atoms);
  m.atoms(atoms);

  int* hcount = new int[matoms];
  std::unique_ptr<int[]> free_hcount(hcount);

  // Force storage of nbonds

  for (int i = 0; i < matoms; i++) {
    Atom* ai = const_cast<Atom*>(atoms[i]);
    ai->nbonds();

    hcount[i] = m.hcount(i);
  }

  atomic_number_t* z = new atomic_number_t[matoms];
  std::unique_ptr<atomic_number_t[]> free_z(z);
  m.atomic_numbers(z);

  int rc = jw_molconn(m, output, atoms, z, hcount);

  if (0 == rc) {
    number_of_error++;
  }

  return rc;
}

static int
jw_molconn(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose) {
      cerr << molecules_read << " processing '" << m->name() << "'\n";
    }

    if (!preprocess_molecule(*m)) {
      if (verbose) {
        cerr << "Cannot process '" << m->name() << "'\n";
      }
      continue;
    }

    m->recompute_distance_matrix();

    if (!jw_molconn(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
jw_molconn(const char* fname, FileType input_type, IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jw_molconn(input, output);
}

static int
jw_molconn_descriptor_pipeline_line(const const_IWSubstring& buffer,
                IWString_and_File_Descriptor& output) {
  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "Bad smiles\n";
    return 0;
  }

  return jw_molconn(m, output);
}

static int
jw_molconn_descriptor_pipeline(iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  static int first_call = 1;
  if (first_call) {
    if (! input.next_record(buffer)) {
      cerr << "jw_molconn_descriptor_pipeline:cannot read header\n";
      return 0;
    }

    if (write_descriptor_file_pipeline) {
      output << buffer << ' ';
      WriteHeader(output);
    } else {
      // Remove the smiles and write.
      buffer.remove_leading_words(1);
      output << buffer << ' ';
      WriteHeader(output);
    }

    first_call = 0;
  }

  while (input.next_record(buffer)) {
    if (! jw_molconn_descriptor_pipeline_line(buffer, output)) {
      cerr << "Error processing\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

static int
jw_molconn_descriptor_pipeline(const char* fname,
                IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "jw_molconn_descriptor_pipeline:cannot open '" << fname << "'\n";
    return 0;
  }

  return jw_molconn_descriptor_pipeline(input, output);
}

static void
DisplayDashYOptions() {
  cerr << " -Y rpipe    input is a descriptor file pipeline  smiles id d1 d2 d3 ...\n";
  cerr << " -Y wpipe    write a descriptor file pipeline - smiles id d1 d2 d3 ...\n";
  ::exit(0);
}

static int
jw_molconn(int argc, char** argv) {
  Command_Line cl(argc, argv, "vli:g:m:p:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  //  if (! process_standard_aromaticity_options(cl, verbose))
  //    {
  //      cerr << "Cannot process aromaticity options (-A)\n";
  //      usage (5);
  //    }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "rpipe") {
        read_descriptor_file_pipeline = 1;
        if (verbose) {
          cerr << "Input assumed to be a descriptor file pipeline\n";
        }
      } else if (y == "wpipe") {
        write_descriptor_file_pipeline = 1;
        if (verbose) {
          cerr << "Will write a descriptor file pipeline\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions();
      } else {
        cerr << "Unrecognised -Y directive '" << y << "'\n";
        DisplayDashYOptions();
      }
    }
  }

  if (cl.option_present('l')) {
    compute_additional_connectivity_related_descriptors = 1;
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', output_precision) || output_precision < 1) {
      cerr << "The output precision must be a whole, +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Results written with " << output_precision << " precision\n";
    }
  }

  set_default_iwstring_float_concatenation_precision(output_precision);

  set_global_aromaticity_type(Daylight);
  set_input_aromatic_structures(1);

  if (cl.option_present('g'))  // only recognised with the -d option
  {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g option)\n";
      usage(8);
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_allowed_heavy_atom)) {
      cerr << "Invalid value for max allowed heavy atom (-m) has to be a positive "
              "integar\n";
      usage(2);
    }

    if (max_allowed_heavy_atom < 100) {
      cerr << "The value for max allowed heavy atoms is not valid, should be at least "
              "100\n";
      cerr << "100 (default) is used instead\n";
      max_allowed_heavy_atom = 100;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (read_descriptor_file_pipeline) {
    // do not need input type
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  // In this case, 0 for simple chi, 1 for valence chi, 2 for tp_value, 3 for etp_value
  Path_with_values::set_number_values(NUMBER_OF_VALUES_FOR_PATH);

  IWString_and_File_Descriptor output(1);
  output.resize(36000);

  if (read_descriptor_file_pipeline) {
    // Header generated elsewhere.
  } else {
    if (write_descriptor_file_pipeline) {
      output << "Smiles" << ' ';
    } 

    output << "Name" << ' ';

    if (!WriteHeader(output)) {
      return 5;
    }
  }

  if (read_descriptor_file_pipeline) {
    for (const char* fname : cl) {
      if (! jw_molconn_descriptor_pipeline(fname, output)) {
        cerr << "Fatal error reading '" << fname << "'\n";
        return 1;
      }
    }
  } else {
    for (const char* fname : cl) {
      if (!jw_molconn(fname, input_type, output)) {
        cerr << "Fatal error reading '" << fname << "'\n";
        return 1;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return 0;
}

int
main(int argc, char** argv) {
  int rc = jw_molconn(argc, argv);
  return rc;
}
