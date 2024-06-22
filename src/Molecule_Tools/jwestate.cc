/*
 * This program compute the ESTATE descriptors
 * Reference:
 * 1. JCICS, 2001, 41, 321-329, Liu et al
 * 2. JCICS, 1999, 39, 951-957, Liu et al
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/e_state_computation_procedure_ori.h"

using std::cerr;
using std::endl;

#define UNDEFINED_E_STATE_ATOM_TYPE -1
#define NUMBER_OF_E_STATE_ATOM_TYPE 47
#define NUMBER_OF_H_E_STATE_ATOM_TYPE 8

#define ELOW -5.0
#define EHIGH 15.0
#define ESIZE 0.5

static Chemical_Standardisation chemical_standardisation;

int verbose = 0;

static int molecules_read = 0;

// keep track of number of warning and number of error
static int number_of_error = 0;

static int output_precision = 3;

static IWDigits iwdigits;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString tag;

static int fingerprint_output_is_atom_estate = 0;
static int fingerprint_output_is_hydrogen_estate = 0;

static int output_vector_size = 20;
static double output_vector_std = 1.0;

static int work_as_tdt_filter = 0;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Compute the JWESTATE descriptors for molecules\n";
  cerr << "  -J <tag>        produce a fingerprint\n";
  cerr << "  -J estate       fingerprint is based on atomic EState values\n";
  cerr << "  -J hestate      fingerprint is based on atomic Hydrogen EState values\n";
  cerr << "  -f              work as a TDT filter\n";
  cerr << "  -i <type>      input type\n";
  (void) display_standard_aromaticity_options (cerr);
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
preprocess_molecule(Molecule& m)
{
  m.reduce_to_largest_fragment();  // always reduce to largest fragment

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return 1;
}

/*
  We observe that generally these values tend to fall in the range of -3 to 18,
  a little different from what Jibo implemented initially
*/

static int
write_fingerprint_precise(Molecule& m, const atomic_number_t* z, const double* xstate,
                          Sparse_Fingerprint_Creator& sfc)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    if (1 == z[i]) {
      continue;
    }

    double x = xstate[i] + 3.0;
    if (x < 0.0) {
      x = 0.0;
    }

    int b = static_cast<int>(x + 1.0);
    sfc.hit_bit(b, 5);
    if (b > 0) {
      sfc.hit_bit(b - 1);
    }
    sfc.hit_bit(b + 1);
  }

  return 1;
}

static double
normal_function(double mean, double std, double x)
{
  return exp(-(x - mean) * (x - mean) / (std * std + std * std)) / std / sqrt(2.0 * M_PI);
}

static int
write_fingerprint_dist(Molecule& m, const atomic_number_t* z, const double* xstate,
                       Sparse_Fingerprint_Creator& sfc)
{
  double* ovector = new double[output_vector_size];
  std::unique_ptr<double[]> free_ovector(ovector);
  set_vector(ovector, output_vector_size, 0.0);

  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    double x = xstate[i] + 3.0;
    if (x < 0.0) {
      x = 0.0;
    }

    const int xi = static_cast<int>(x + 0.49999);
    int jstart = xi - 3;
    if (jstart < 0) {
      jstart = 0;
    }
    int jstop = xi + 3;
    if (jstop > output_vector_size) {
      jstop = output_vector_size;
    }

    //  cerr << " i = " << i << " x " << x << " jstart " << jstart << " jstop " << jstop
    //  << endl;
    for (auto j = jstart; j < jstop; ++j) {
      double n = normal_function(x, output_vector_std, static_cast<double>(j));
      ovector[j] += n;
    }
  }

  for (auto i = 0; i < output_vector_size; ++i) {
    //  cerr << " i " << i << " ov " << ovector[i] << endl;

    if (ovector[i] < 0.01) {
      continue;
    }

    sfc.hit_bit(i, static_cast<int>(ovector[i] * 10.0 + 0.4999));
  }

  return 1;
}

static int
write_fingerprint(Molecule& m, const atomic_number_t* z, const double* xstate,
                  IWString_and_File_Descriptor& output)
{
  if (!work_as_tdt_filter) {
    m.remove_all(1);
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator sfc;

  int rc;
  if (output_vector_std >= 0.0) {
    rc = write_fingerprint_dist(m, z, xstate, sfc);
  } else {
    rc = write_fingerprint_precise(m, z, xstate, sfc);
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, output);

  output << tmp << '\n';

  if (!work_as_tdt_filter) {
    output << "|\n";
  }

  return rc;
}

static int
output_result_header(IWString_and_File_Descriptor& output)
{
  output << "Name";

  // print out header for e state descriptor
  for (int i = 0; i < NUMBER_OF_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWESSU" << i;
  }
  for (int i = 0; i < NUMBER_OF_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWESMA" << i;
  }
  for (int i = 0; i < NUMBER_OF_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWESMI" << i;
  }
  for (int i = 0; i < NUMBER_OF_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWESAV" << i;
  }
  for (int i = 0; i < NUMBER_OF_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWESCT" << i;
  }

  // print out header for h_e state descriptor
  for (int i = 0; i < NUMBER_OF_H_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWHESU" << i;
  }
  for (int i = 0; i < NUMBER_OF_H_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWHEMA" << i;
  }
  for (int i = 0; i < NUMBER_OF_H_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWHEMI" << i;
  }
  for (int i = 0; i < NUMBER_OF_H_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWHEAV" << i;
  }
  for (int i = 0; i < NUMBER_OF_H_E_STATE_ATOM_TYPE; i++) {
    output << " estate_JWHECT" << i;
  }

  int n_bin = int(rint((EHIGH - ELOW) / ESIZE)) + 1;
  for (int i = 0; i < n_bin; i++) {
    output << " estate_JWEBIN" << i;
  }

  output << '\n';

  return 1;
}

static int
count_doubly_bonded_oxygens_attached(const Molecule& m, atom_number_t centre,
                                     atom_number_t& o1, atom_number_t& o2)
{
  int doubly_bonded_oxygens_found = 0;

  const Atom* acs = m.atomi(centre);

  for (int i = 0; i < acs->ncon(); i++) {
    const Bond* b = acs->item(i);
    if (!b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(centre);

    if (8 != m.atomic_number(j)) {
      continue;
    }

    if (0 == doubly_bonded_oxygens_found) {
      o1 = j;
    } else if (1 == doubly_bonded_oxygens_found) {
      o2 = j;
    } else {
      //	  cerr << "is_charged_acid: too many doubly bonded oxygens, ignoring atom
      //" << centre << endl;
      return 1;
    }

    doubly_bonded_oxygens_found++;
  }

  //  cerr << "Atom " << centre << " has " << doubly_bonded_oxygens_found << " doubly
  //  bonded oxygens\n";

  return doubly_bonded_oxygens_found;
}

/*
 *  Atom O1 is a charged oxygen atom. Is it part of a carboxyllic or sulf* acid
 * If a carboxyllic acid, we set O2. If Sulphonic, we also set O3
 */

int
is_charged_acid(const Molecule& m, atom_number_t o1, int atom_type[])
{
  const Atom* a1 = m.atomi(o1);

  assert(8 == a1->atomic_number());
  assert(-1 == a1->formal_charge());

  atom_number_t centre =
      a1->other(o1, 0);  // 1 negatively charged oxygen has only one connection

  const Atom* acs = m.atomi(centre);

  if (acs->ncon() < 3) {
    return 0;
  }

  if (acs->ncon() == acs->nbonds()) {  // central must be unsaturated
    return 0;
  }

  if (6 == acs->atomic_number())  // possible carboxyllic acid
  {
    if (3 != acs->ncon()) {
      return 0;
    }
  } else if (16 == acs->atomic_number())  // possible sulph* acid
  {
    if (4 != acs->ncon()) {
      return 0;
    }
  } else if ((7 == acs->atomic_number()) &&
             (1 == acs->formal_charge()))  // possible N oxide
  {
    if (3 != acs->ncon()) {
      return 0;
    }
  } else {
    return 0;
  }

  atom_number_t o2, o3;
  // These initialisations are not needed, just to keep the compiler quiet.
  o2 = kInvalidAtomNumber;
  o3 = kInvalidAtomNumber;
  int no = count_doubly_bonded_oxygens_attached(m, centre, o2, o3);

  //  cerr << "Atom " << centre << " has " << no << " doubly bonded oxygens\n";
  if (0 == no) {
    return 0;
  }

  atom_type[o1] = 30;

  // Carboxyllic or sulph* type acid?

  if (16 == acs->atomic_number()) {
    atom_type[centre] = 44;
  }
  if (7 == acs->atomic_number()) {
    atom_type[centre] = 25;
  }

  atom_type[o2] = 30;
  if (no > 1) {
    atom_type[o3] = 30;
  }

  return 1;
}

/*
 *This funciton determine the e state atom type
 */

static void
determine_e_state_atom_type(Molecule& m, int e_state_atom_type[],
                            const atomic_number_t* z)
{
  int n_atoms = m.natoms();

  set_vector(e_state_atom_type, n_atoms, UNDEFINED_E_STATE_ATOM_TYPE);

  for (int i = 0; i < n_atoms; i++) {
    if (UNDEFINED_E_STATE_ATOM_TYPE != e_state_atom_type[i]) {
      continue;
    }
    int h_count = m.hcount(i);
    int nconi = m.ncon(i);
    int nbondsi = m.nbonds(i);
    switch (z[i]) {
      case 1:
        break;
      case 6:
        if (nconi == nbondsi) {
          if (3 == h_count) {
            e_state_atom_type[i] = 0;
          } else if (2 == h_count) {
            e_state_atom_type[i] = 1;
          } else if (1 == h_count) {
            e_state_atom_type[i] = 2;
          } else if (0 == h_count) {
            e_state_atom_type[i] = 3;
          }
        } else if (nconi == nbondsi - 1) {
          int bonds_in_conjugated_system = number_of_bond_in_conjugated_system(m, i);
          if (2 == h_count) {
            e_state_atom_type[i] = 4;
          } else if (1 == h_count) {
            if (2 == bonds_in_conjugated_system) {
              e_state_atom_type[i] = 10;
            } else {
              e_state_atom_type[i] = 5;
            }
          } else if (0 == h_count) {
            if (2 == bonds_in_conjugated_system) {
              e_state_atom_type[i] = 11;
            } else if (3 == bonds_in_conjugated_system) {
              e_state_atom_type[i] = 12;
            } else {
              e_state_atom_type[i] = 6;
            }
          }
        } else if (nconi == nbondsi - 2) {
          if (!is_triple_bonded(m, i)) {
            e_state_atom_type[i] = 7;
          } else {
            if (1 == h_count) {
              e_state_atom_type[i] = 8;
            } else if (0 == h_count) {
              e_state_atom_type[i] = 9;
            }
          }
        }
        break;
      case 7:
        if (nconi == nbondsi) {
          const formal_charge_t f_charge_i = m.formal_charge(i);
          if (0 == f_charge_i) {
            if (2 == h_count) {
              e_state_atom_type[i] = 13;
            } else if (1 == h_count) {
              e_state_atom_type[i] = 14;
            } else if (0 == h_count) {
              e_state_atom_type[i] = 20;
            }
          } else if (1 == f_charge_i) {
            if (3 == h_count) {
              e_state_atom_type[i] = 15;
            }
            if (2 == h_count) {
              e_state_atom_type[i] = 16;
            } else if (1 == h_count) {
              e_state_atom_type[i] = 17;
            } else if (0 == h_count) {
              e_state_atom_type[i] = 18;
            }
          }
        } else if (nconi == nbondsi - 1) {
          int bonds_in_conjugated_system = number_of_bond_in_conjugated_system(m, i);

          if (1 == h_count) {
            if (2 == bonds_in_conjugated_system) {
              e_state_atom_type[i] = 23;
            } else {
              e_state_atom_type[i] = 19;
            }
          } else if (0 == h_count) {
            if (2 == bonds_in_conjugated_system) {
              e_state_atom_type[i] = 24;
            } else {
              e_state_atom_type[i] = 21;
            }
          }
        } else if ((nconi == nbondsi - 2) && (1 == nconi)) {
          e_state_atom_type[i] = 22;
        } else if ((nconi == nbondsi - 2) && (3 == nconi)) {
          e_state_atom_type[i] = 26;
        }

        break;
      case 8:
        if (-1 == m.formal_charge(i) && is_charged_acid(m, i, e_state_atom_type)) {
          ;
        } else if (nconi == nbondsi) {
          if (1 == h_count) {
            e_state_atom_type[i] = 27;
          } else if (0 == h_count) {
            e_state_atom_type[i] = 28;
          }
        } else if ((1 == nconi) && (2 == nbondsi)) {
          e_state_atom_type[i] = 29;
        }
        break;
      case 9:
        e_state_atom_type[i] = 32;
        break;
      case 15:
        if ((5 == nconi) && (5 == nbondsi)) {
          e_state_atom_type[i] = 38;
        } else if ((4 == nconi) && (5 == nbondsi)) {
          e_state_atom_type[i] = 39;
        } else if (nconi == nbondsi) {
          if (2 == h_count) {
            e_state_atom_type[i] = 35;
          } else if (1 == h_count) {
            e_state_atom_type[i] = 36;
          } else if (0 == h_count) {
            e_state_atom_type[i] = 37;
          }
        }
        break;
      case 16:
        if ((6 == nconi) && (6 == nbondsi)) {
          e_state_atom_type[i] = 45;
        } else if ((4 == nconi) && (6 == nbondsi)) {
          e_state_atom_type[i] = 46;
        } else if ((3 == nconi) && (4 == nbondsi)) {
          e_state_atom_type[i] = 43;
        } else if (nconi == nbondsi) {
          if (1 == h_count) {
            e_state_atom_type[i] = 40;
          } else if (0 == h_count) {
            e_state_atom_type[i] = 41;
          }
        } else if (nconi == nbondsi - 1) {
          e_state_atom_type[i] = 42;
        }
        break;
      case 17:
        e_state_atom_type[i] = 31;
        break;
      case 35:
        e_state_atom_type[i] = 33;
        break;
      case 53:
        e_state_atom_type[i] = 34;
        break;
    }
  }
}

/*
 * compute over-all e_state descriptors
 */

static int
compute_min_max_sum_average(int n_atoms, int NUMBER_OF_TYPE,
                            IWString_and_File_Descriptor& output,
                            const int* e_state_atom_type, const double* e_state_index,
                            double* sum, double* max, double* min, double* avg,
                            int count[])
{
  // initialization
  set_vector(sum, NUMBER_OF_TYPE, 0.0);
  set_vector(max, NUMBER_OF_TYPE, 0.0);
  set_vector(min, NUMBER_OF_TYPE, std::numeric_limits<double>::max());
  set_vector(avg, NUMBER_OF_TYPE, 0.0);
  set_vector(count, NUMBER_OF_TYPE, 0);

  for (int i = 0; i < n_atoms; i++) {
    int atom_type = e_state_atom_type[i];
    if (UNDEFINED_E_STATE_ATOM_TYPE != atom_type) {
      sum[atom_type] += e_state_index[i];
      count[atom_type]++;
      if (max[atom_type] < e_state_index[i]) {
        max[atom_type] = e_state_index[i];
      }
      if (min[atom_type] > e_state_index[i]) {
        min[atom_type] = e_state_index[i];
      }
    }
  }

  for (int i = 0; i < NUMBER_OF_TYPE; i++) {
    if (count[i] > 0) {
      avg[i] = sum[i] / count[i];
    }
  }

  for (int i = 0; i < NUMBER_OF_TYPE; i++) {
    if (count[i] > 0) {
      output << ' ' << static_cast<float>(sum[i]);
    } else {
      output << " 0";
    }
  }

  for (int i = 0; i < NUMBER_OF_TYPE; i++) {
    if (count[i] > 0) {
      output << ' ' << static_cast<float>(max[i]);
    } else {
      output << " 0";
    }
  }

  for (int i = 0; i < NUMBER_OF_TYPE; i++) {
    if (count[i] > 0) {
      output << ' ' << static_cast<float>(min[i]);
    } else {
      output << " 0";
    }
  }

  for (int i = 0; i < NUMBER_OF_TYPE; i++) {
    if (count[i] > 0) {
      output << ' ' << static_cast<float>(avg[i]);
    } else {
      output << " 0";
    }
  }

  for (int i = 0; i < NUMBER_OF_TYPE; i++) {
    iwdigits.append_number(output, count[i]);
  }

  return 1;
}

/*
 * declare and delete some arrays
 */

static int
compute_min_max_sum_average(int n_atoms, const int types,
                            IWString_and_File_Descriptor& output, int atom_type[],
                            double index[], double* d, int* at_count)
{
  return compute_min_max_sum_average(n_atoms, types, output, atom_type, index, d,
                                     d + types, d + 2 * types, d + 3 * types, at_count);
}

/*
 * bined e_state index
 */
static int
e_state_bin_descriptor_computation_procedure(int n_atoms, double low, double high,
                                             double size, const int atom_type[],
                                             const double index[],
                                             IWString_and_File_Descriptor& output)
{
  IWHistogram histogram;
  histogram.initialise(low, high, size);

  for (int i = 0; i < n_atoms; i++) {
    if (UNDEFINED_E_STATE_ATOM_TYPE != atom_type[i]) {
      histogram.extra(index[i]);
    }
  }

  int n_bucket = histogram.nbuckets();
  const unsigned int* count = histogram.raw_counts();

  for (int i = 0; i < n_bucket; i++) {
    iwdigits.append_number(output, count[i]);
  }

  return 1;
}

/*
 * declare and delete some arrays
 */

static int
generate_output(Molecule& m, IWString_and_File_Descriptor& output, double e_state_index[],
                int e_state_atom_type[], double h_e_state_index[],
                int h_e_state_atom_type[])

{
  append_first_token_of_name(m.name(), output);

  const int n_atoms = m.natoms();

  double* d = new double[4 * NUMBER_OF_E_STATE_ATOM_TYPE];
  std::unique_ptr<double[]> free_d(
      d);  // important that NUMBER_OF_E_STATE_ATOM_TYPE []> NUMBER_OF_H_E_STATE_ATOM_TYPE
  int* at_count = new int[NUMBER_OF_E_STATE_ATOM_TYPE];
  std::unique_ptr<int[]> free_at_count(at_count);

  int rc = compute_min_max_sum_average(n_atoms, NUMBER_OF_E_STATE_ATOM_TYPE, output,
                                       e_state_atom_type, e_state_index, d, at_count);

  rc = compute_min_max_sum_average(n_atoms, NUMBER_OF_H_E_STATE_ATOM_TYPE, output,
                                   h_e_state_atom_type, h_e_state_index, d, at_count);

  rc = e_state_bin_descriptor_computation_procedure(
      n_atoms, ELOW, EHIGH, ESIZE, e_state_atom_type, e_state_index, output);

  output << '\n';

  return rc;
}

static int
bonded_to_halogen(Molecule& m, atom_number_t zatom)
{
  const Atom* a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(zatom, i);

    int an = m.atomic_number(j);
    if (6 == an || 7 == an) {  // common cases
      ;
    } else if ((53 == an) || (35 == an) || (9 == an) || (17 == an)) {
      return 1;
    }
  }

  return 0;
}

/*
 * determine H e state atom type
 */

static void
determine_h_e_state_atom_type(Molecule& m, int h_e_state_atom_type[],
                              const atomic_number_t* z)
{
  int n_atoms = m.natoms();

  set_vector(h_e_state_atom_type, n_atoms, UNDEFINED_E_STATE_ATOM_TYPE);

  for (int i = 0; i < n_atoms; i++) {
    int h_count = m.hcount(i);
    int nconi = m.ncon(i);
    int nbondsi = m.nbonds(i);

    switch (z[i]) {
      case 1:
        break;
      case 6:
        if ((nconi == nbondsi - 2) && (1 == h_count)) {
          h_e_state_atom_type[i] = 5;
        } else if ((nconi == nbondsi) && bonded_to_halogen(m, i) && (h_count > 0)) {
          h_e_state_atom_type[i] = 6;
        } else if (h_count > 0) {
          h_e_state_atom_type[i] = 7;
        }
        break;
      case 7:
        if (nconi == nbondsi) {
          if (2 == h_count) {
            h_e_state_atom_type[i] = 3;
          } else if (1 == h_count) {
            h_e_state_atom_type[i] = 4;
          }
        } else if ((nconi == nbondsi - 1) && (1 == h_count)) {
          h_e_state_atom_type[i] = 1;
        } else if (h_count > 0) {
          h_e_state_atom_type[i] = 7;
        }
        break;
      case 8:
        if ((nconi == nbondsi) && (1 == h_count)) {
          h_e_state_atom_type[i] = 0;
        } else if (h_count > 0) {
          h_e_state_atom_type[i] = 7;
        }
        break;
      case 16:
        if ((nconi == nbondsi) && (1 == h_count)) {
          h_e_state_atom_type[i] = 2;
        } else if (h_count > 0) {
          h_e_state_atom_type[i] = 7;
        }
        break;
      default:
        if (h_count > 0) {
          h_e_state_atom_type[i] = 7;
        }
    }
  }
}

/*
 * This function is to calculate the MEDV descriptors for molecule
 */

static int
jw_e_state(Molecule& m, IWString_and_File_Descriptor& output, double e_state_index[],
           double h_e_state[], int atom_type[], int h_atom_type[])
{
  m.recompute_distance_matrix();

  const int matoms = m.natoms();

  atomic_number_t* z = new atomic_number_t[matoms];
  std::unique_ptr<atomic_number_t[]> free_z(z);

  m.atomic_numbers(z);

  determine_atom_e_state_index(m, e_state_index, z);

  if (fingerprint_output_is_atom_estate) {
    return write_fingerprint(m, z, e_state_index, output);
  }

  determine_hydrogen_e_state_index(m, h_e_state, z);

  if (fingerprint_output_is_hydrogen_estate) {
    return write_fingerprint(m, z, h_e_state, output);
  }

  determine_e_state_atom_type(m, atom_type, z);

  determine_h_e_state_atom_type(m, h_atom_type, z);

  generate_output(m, output, e_state_index, atom_type, h_e_state, h_atom_type);

  /*  int n_atoms = m.natoms();
  for (int i=0; i<n_atoms; i++)
  cerr<<"atom i="<<i<<"\testate_index ="<<e_state_index[i]<<endl;*/

  return 1;
}

/*
 * Introduce and delete arrays to separate them from the main program
 */

static int
jw_e_state(Molecule& m, IWString_and_File_Descriptor& output)
{
  m.make_implicit_hydrogens_explicit();  // all hydrogen is needed for the computation of
                                         // connections and bonds

  const int n_atoms = m.natoms();

  double* tmpd = new double[n_atoms + n_atoms];
  std::unique_ptr<double[]> free_tmpd(tmpd);
  int* tmpi = new int[n_atoms + n_atoms];
  std::unique_ptr<int[]> free_tmpi(tmpi);

  int rc = jw_e_state(m, output, tmpd, tmpd + n_atoms, tmpi, tmpi + n_atoms);

  if (0 == rc) {
    number_of_error++;
    if (verbose) {
      cerr << "ERROR"
           << "\tMolecule " << m.name()
           << "\tUnrecognizable molecule structure, notify Jibo" << endl;
    }
  }

  return rc;
}

/*
 * deal with each molecule
 */
static int
jw_e_state(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output)
{
  set_default_iwstring_float_concatenation_precision(output_precision);

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose) {
      cerr << molecules_read << " processing '" << m->name() << "'\n";
    }

    (void)preprocess_molecule(*m);

    if (!jw_e_state(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

/*
 * deal with each files
 */

static int
jw_e_state(const char* fname, FileType input_type, IWString_and_File_Descriptor& output)
{
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jw_e_state(input, output);
}

static int
jw_e_state_filter(const_IWSubstring buffer,  // local copy
                  IWString_and_File_Descriptor& output)
{
  buffer.remove_leading_chars(smiles_tag.length());
  buffer.chop();

  Molecule m;

  if (!m.build_from_smiles(buffer)) {
    cerr << "Cannot interpret smiles '" << buffer << "'\n";
    return 0;
  }

  return jw_e_state(m, output);
}

static int
jw_e_state_filter(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!jw_e_state_filter(buffer, output)) {
      return 0;
    }
  }

  return 1;
}

/*
 * deal with all the commandline options
 */

static int
jw_e_state(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:i:g:J:fl");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  set_global_aromaticity_type(Pearlman);

  if (cl.option_present('p')) {
    if (!cl.value('p', output_precision) || output_precision < 1) {
      cerr << "The output precision must be a whole, +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Results written with " << output_precision << " precision\n";
    }
  }

  if (cl.option_present('g'))  // only recognised with the -d option
  {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g option)\n";
      usage(8);
    }
  }

  if (cl.option_present('J')) {
    const_IWSubstring j;

    for (auto i = 0; cl.value('J', j, i); ++i) {
      if ("estate" == j) {
        fingerprint_output_is_atom_estate = 1;
      } else if ("hestate" == j) {
        fingerprint_output_is_hydrogen_estate = 1;
      } else if (j.starts_with("std=")) {
        j.remove_leading_chars(4);
        if (!j.numeric_value(output_vector_std) || output_vector_std < 0.0) {
          cerr << "The standard deviation on the output values must be positive '" << j
               << "' invalid\n";
          return 2;
        }

        if (verbose) {
          cerr << "Output vector will be a distribution, with standard deviation "
               << output_vector_std << '\n';
        }
      } else if (tag.empty()) {
        tag = j;
      } else {
        cerr << "Unrecognised -J qualifier '" << j << "'\n";
        usage(1);
      }
    }

    if (0 == fingerprint_output_is_atom_estate &&
        0 == fingerprint_output_is_hydrogen_estate) {
      cerr << "Must specify what index to use for fingerprint formation\n";
      usage(1);
    }

    if (tag.length()) {
      if (!tag.ends_with('<')) {
        tag << '<';
      }
    } else if (fingerprint_output_is_atom_estate) {
      tag = "NCESTA<";
    } else if (fingerprint_output_is_hydrogen_estate) {
      tag = "NCESTH<";
    }

    if (verbose) {
      cerr << "Will write output as fingerprint with tag '" << tag << "'\n";
    }

    if (cl.option_present('f')) {
      work_as_tdt_filter = 1;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (work_as_tdt_filter) {
    ;
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

  iwdigits.set_include_leading_space(1);
  iwdigits.initialise(50);

  IWString_and_File_Descriptor output(1);

  if (tag.length()) {
    ;
  } else if (!output_result_header(output)) {
    return 3;
  }

  int rc = 0;
  if (work_as_tdt_filter) {
    iwstring_data_source input(0);

    if (!jw_e_state_filter(input, output)) {
      rc = 2;
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!jw_e_state(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  // output info
  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = jw_e_state(argc, argv);

  return rc;
}
