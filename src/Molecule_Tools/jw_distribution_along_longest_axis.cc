/*
 * This program designed to see distribution of various properties along the longest axis
 * of the molecules The molecule can be 2D or 3D.
 */

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/target.h"

#include "jw_common_alignment.h"

using std::cerr;

// DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

#define NATMTYPE 10  // type of atom

static int output_precision = 4;

static resizable_array_p<Substructure_Hit_Statistics> queries;
static int number_of_regions = 10;
static int number_of_properties = 0;

static Charge_Assigner charge_assigner;

int verbose = 0;

static int align_by_longest_distance = 0;

static int molecules_read = 0;

static int ntest = 0;
static int molecules_failing_region_assignment_test = 0;

// keep track of number of warning and number of error
static int number_of_error = 0;

static Molecule_Output_Object stream_for_oriented_molecules;

static IWDigits iwdigits;

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
  cerr << "Compute the JWDIST descriptors for molecules\n";

  cerr << "  -r <int>       number of regions to run the distribution on (default 10)\n";
  cerr << "  -d             align by longest distance\n";

  cerr << "  -i <type>      input type\n";

  cerr << "  -q <queries>   queries you want to see the distribution of. -q help for detail"<<'\n';

  (void) display_standard_charge_assigner_options (cerr, 'N');
  (void) display_standard_aromaticity_options (cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static void
display_q_directives(int rc) {
  cerr << "  -q QUERY:qfile    query is in <qfile>\n";
  cerr << "  -q QUERY:F:file       queries are in <file>\n";
  cerr << "  -q SMARTS:smt     query is written in smarts\n";

  exit(rc);
}

static int
preprocess_molecule(Molecule& m) {
  m.remove_all(1);  // hydrogens not necessary

  m.remove_all_chiral_centres();

  m.reduce_to_largest_fragment();  // always reduce to largest fragment

  // assign formal charges here
  if (charge_assigner.active()) {
    (void)charge_assigner.process(m);
  }

  return 1;
}

static int
output_result_header(IWString_and_File_Descriptor& output) {
  output << "Name";

  // print out header for bcut descriptor
  for (int i = 0; i < number_of_properties; i++) {
    output << " jwdist_JWD" << i << "FST";
    output << " jwdist_JWD" << i << "LST";
    output << " jwdist_JWD" << i << "SUM";
    output << " jwdist_JWD" << i << "MIN";
    output << " jwdist_JWD" << i << "MAX";
    output << " jwdist_JWD" << i << "AV1";
    output << " jwdist_JWD" << i << "AV2";
    output << " jwdist_JWD" << i << "DE1";
    output << " jwdist_JWD" << i << "DE2";
    output << " jwdist_JWD" << i << "AVP";
  }

  output << '\n';

  return 1;
}

static int
assign_regions(Molecule& m, coord_t xmin, coord_t xmax, int atom_in_region[]) {
  const auto n_atoms = m.natoms();

  // since molecule is now oriented with the longest axis on the x axis, we can assign
  // atoms into differnet region

  float width_of_region = (xmax - xmin) / (float)number_of_regions;
  // cerr << "Range is " << xmin << " to " << xmax << " width_of_region " <<
  // width_of_region << '\n';

  for (int i = 0; i < n_atoms; i++) {
    int region = (int)((m.x(i) - xmin) / width_of_region);
    //    cerr << "Atom " << i << " x " << m.x(i) << " region " << region << ", min " <<
    //    xmin << " width_of_region " << width_of_region << '\n';
    if (region >= number_of_regions) {
      atom_in_region[i] = number_of_regions - 1;
    } else if (region < 0) {
      atom_in_region[i] = 0;
    } else {
      atom_in_region[i] = region;
    }
  }

#define REPORT_REGIONS
#ifdef REPORT_REGIONS
  if (verbose > 1) {
    for (auto i = 0; i < n_atoms; ++i) {
      m.set_isotope(i, atom_in_region[i] + 1);
    }
    cerr << m.unique_smiles() << ' ' << m.name() << '\n';
  }
#endif

  return 1;
}

static int
setup_atom_in_region(Molecule& m, int atom_in_region[], double** coordinates,
                     double** original_matrix, double temp_array[], double** axis) {
  int n_atoms = m.natoms();
  if (n_atoms < 2) {
    return 1;
  }

  // firt order of business -- take care of the arrangement of atoms in various region

  (void)orient_molecule_to_rotational_axis(m, coordinates, original_matrix, temp_array,
                                           axis);

  coord_t xmin, xmax;
  m.spatial_extremeties_x(xmin, xmax);

  assign_regions(m, xmin, xmax, atom_in_region);

  return 1;
}

static int
do_align_by_longest_distance(Molecule& m, int atom_in_region[]) {
  atom_number_t left, right;

  m.rotate_to_longest_distance_along_x(left, right);

  assign_regions(m, m.atomi(left)->x(), m.atomi(right)->x(), atom_in_region);

  return 1;
}

static int
setup_atom_in_region(Molecule& m, int atom_in_region[]) {
  int n_atoms = m.natoms();

  double** coordinates = new double*[n_atoms];
  for (int i = 0; i < n_atoms; i++) {
    coordinates[i] = new double[3];
  }

  double** original_matrix = new double*[3];
  for (int i = 0; i < 3; i++) {
    original_matrix[i] = new double[3];
  }

  double temp_array[3];

  double** axis = new double*[3];
  for (int i = 0; i < 3; i++) {
    axis[i] = new double[3];
  }

  int rc = setup_atom_in_region(m, atom_in_region, coordinates, original_matrix,
                                temp_array, axis);

  for (int i = 0; i < n_atoms; i++) {
    delete[] coordinates[i];
  }

  delete[] coordinates;

  for (int i = 0; i < 3; i++) {
    delete[] original_matrix[i];
    delete[] axis[i];
  }

  delete[] original_matrix;
  delete[] axis;

  return rc;
}

/*
  We return the first non zero element
*/

static int
enforce_canonical_order(int* temp_array, const int n) {
  int first_non_zero = -1;

  for (auto i = 0; i < n; ++i) {
    if (first_non_zero < 0 && temp_array[i] > 0) {
      first_non_zero = i;
    }

    if (temp_array[i] > temp_array[n - 1 - i]) {
      return first_non_zero;
    }

    if (temp_array[i] == temp_array[n - 1 - i]) {
      continue;
    }

    int tmp;
    for (auto j = 0, k = n - 1; j < n / 2; ++j, --k) {
      tmp = temp_array[j];
      temp_array[j] = temp_array[k];
      temp_array[k] = tmp;
      if (first_non_zero < 0 && temp_array[j] > 0) {
        first_non_zero = j;
      }
    }
    for (auto j = 0; j < n; ++j) {
      if (temp_array[j] > 0) {
        return j;
      }
    }
    return -1;
  }

  return first_non_zero;
}

static int
compute_distribution_of_queried_properties(Molecule& m,
                                           IWString_and_File_Descriptor& output,
                                           int** hit_array, int atom_in_region[],
                                           int** already_hit, int temp_array[]) {
  int n_atoms = m.natoms();

  if (n_atoms < 2) {
    return 1;
  }

  for (int i = 0; i < number_of_properties; i++) {
    set_vector(hit_array[i], number_of_regions, 0);
    set_vector(already_hit[i], n_atoms, 0);
  }

  if (align_by_longest_distance) {
    do_align_by_longest_distance(m, atom_in_region);
  } else {
    (void)setup_atom_in_region(m, atom_in_region);
  }

  if (stream_for_oriented_molecules.active()) {
    stream_for_oriented_molecules.write(m);
  }

// #define DEBUG_COMPUTE_QUERIES_PROPERTIES
#ifdef DEBUG_COMPUTE_QUERIES_PROPERTIES
  for (auto i = 0; i < m.natoms(); ++i) {
    cerr << "Atom " << i << " in region " << atom_in_region[i] << '\n';
  }
#endif

  int n_queries = queries.number_elements();

  Molecule_to_Match target(&m);

  for (int i = 0; i < n_queries; i++) {
    Substructure_Results sresults;
    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    double tmp;
    (void)queries[i]->numeric_value(tmp);

    int n = (int)tmp - 1;
    //    cerr << "N set to " << n << '\n';

    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);  // loss of const OK

      int hit_atom = e->item(0);

      // cerr<<"already_hit[n][hit_atom]="<<already_hit [n][hit_atom]<<'\n';
      if (already_hit[n][hit_atom] == 0) {
        // cerr << "hit atom " << hit_atom << " comes from region " <<
        // atom_in_region[hit_atom] << '\n';
        hit_array[n][atom_in_region[hit_atom]]++;

        already_hit[n][hit_atom] = 1;
      }
    }
  }

  for (int i = 0; i < number_of_properties; i++) {
    for (auto j = 0; j < number_of_regions; ++j) {
      temp_array[j] = hit_array[i][j];
    }

    /*  cerr << i <<  " B4 Regions:";
        for (auto j = 0; j < number_of_regions; ++j) {
          cerr << ' ' << temp_array[j];
        }
        cerr << '\n';*/

    int first_non_zero = enforce_canonical_order(temp_array, number_of_regions);

    if (first_non_zero < 0) {
      first_non_zero = number_of_regions;  // will stop the loop below executing
    }

    /*  cerr << i << " AF Regions:";
        for (auto j = 0; j < number_of_regions; ++j) {
          cerr << ' ' << temp_array[j];
        }
        cerr << '\n';*/

    int sum = 0;
    int no_regions_occupied = 0;
    int last_region_occupied = 0;
    int first_region_occupied = -1;
    int max = 0;

    int min;
    if (first_non_zero < 0) {
      min = 0;
    } else {
      min = m.natoms();
    }

    int sum_position = 0;  // this is to compute the average position of the feature

    for (int j = first_non_zero; j < number_of_regions; j++) {
      if (temp_array[j] > 0) {
        sum += temp_array[j];
        sum_position += temp_array[j] * (j + 1);

        // increment number of occupied regions
        no_regions_occupied++;

        // update last region occupied
        last_region_occupied = j;

        // update first region occupied
        if (-1 == first_region_occupied) {
          first_region_occupied = j;
        }

        // update max
        if (temp_array[j] > max) {
          max = temp_array[j];
        }
      }

      // update min
      if (temp_array[j] < min) {
        min = temp_array[j];
      }
    }

    float avg_1;
    float avg_2;
    if (no_regions_occupied > 0) {
      avg_1 = (float)sum / (float)number_of_regions;
      avg_2 = (float)sum / (float)no_regions_occupied;
    } else {
      avg_1 = 0.0f;
      avg_2 = 0.0f;
    }

    float square_deviation_1 = 0.0;
    float square_deviation_2 = 0.0;

    if (no_regions_occupied > 0) {
      for (int j = 0; j < number_of_regions; j++) {
        float temp = temp_array[j] - avg_1;
        square_deviation_1 += temp * temp;

        if (temp_array[j] > 0) {
          square_deviation_2 += temp * temp;
        }
      }
    }

    float dev_1;
    float dev_2;

    if (no_regions_occupied > 0) {
      dev_1 = sqrt(square_deviation_1 / number_of_regions);
      dev_2 = sqrt(square_deviation_2 / no_regions_occupied);
    } else {
      dev_1 = 0.0f;
      dev_2 = 0.0f;
    }

    float average_position;
    if (sum > 0) {
      average_position = (float)sum_position / (float)sum / (float)number_of_regions;
    } else {
      average_position = 0.0f;
    }

    iwdigits.append_number(output, first_region_occupied);
    iwdigits.append_number(output, last_region_occupied);
    iwdigits.append_number(output, sum);
    iwdigits.append_number(output, min);
    iwdigits.append_number(output, max);

    //    output << ' ' << first_region_occupied<< ' ' << last_region_occupied << ' ' <<
    //    sum; output << ' ' << min << ' ' << max << ' ' << avg_1 << ' ' << avg_2;
    output << ' ' << avg_1 << ' ' << avg_2;
    output << ' ' << dev_1 << ' ' << dev_2 << ' ' << average_position;
  }

  /*
  for (int i=0; i<number_of_properties; i++) {
      for (int j=0; j<number_of_regions; j++) {
          cerr<<hit_array[i][j]<<"    ";
        }
      cerr<<'\n';
    }
  */

  return 1;
}

/*
 * This function is to calculate the dist descriptors for molecule
 */

static int
compute_distribution_of_queried_properties(Molecule& m,
                                           IWString_and_File_Descriptor& output) {
  int n_atoms = m.natoms();

  int* temp_array = new int[number_of_regions];
  std::unique_ptr<int[]> free_tmp_array(temp_array);

  int** hit_array = new int*[number_of_properties];
  for (int i = 0; i < number_of_properties; i++) {
    hit_array[i] = new int[number_of_regions];
  }

  int** already_hit = new int*[number_of_properties];
  for (int i = 0; i < number_of_properties; i++) {
    already_hit[i] = new int[n_atoms];
  }

  int* atom_in_region = new int[n_atoms];
  std::unique_ptr<int[]> free_atom_in_region(atom_in_region);

  int rc = compute_distribution_of_queried_properties(
      m, output, hit_array, atom_in_region, already_hit, temp_array);

  for (int i = 0; i < number_of_properties; i++) {
    delete[] hit_array[i];
    delete[] already_hit[i];
  }

  delete[] hit_array;
  delete[] already_hit;

  return rc;
}

static int
jwdist(Molecule& m, IWString_and_File_Descriptor& output) {
  // if we will need the unique smiles to store in the database, generate it now before
  // the charge assigner and such. Since we don't know what state the input molecule is
  // in, make a copy and apply chemical standardisation

  if (m.natoms() < 2) {
    cerr << "Only one atom for this compounds, computation skipped" << '\n';
    return 1;
  }

  (void)preprocess_molecule(m);

  if (m.name().contains(' ')) {
    const_IWSubstring tmp(m.name());
    tmp.truncate_at_first(' ');
    output << tmp;
  } else {
    output << m.name();
  }

  compute_distribution_of_queried_properties(m, output);

  output << '\n';

  return 1;
}

static int
jwdist(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output) {
  output.resize(70000);

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1) {
      cerr << molecules_read << " processing '" << m->name() << "'\n";
    }

    if (!jwdist(*m, output)) {
      number_of_error++;
      return 0;
    }

    output.write_if_buffer_holds_more_than(65536);
  }

  return 1;
}

static int
jwdist(const char* fname, FileType input_type, IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jwdist(input, output);
}

static int
jwdist(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:q:r:N:S:t:dO:");

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

  // charge assigner option
  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage(33);
    }

    charge_assigner.set_min_distance_between_charges(1);
  }

  if (cl.option_present('q')) {
    const_IWSubstring f;
    int i = 0;
    while (cl.value('q', f, i++)) {
      if (f.starts_with("SMARTS:")) {
        f.remove_leading_chars(7);
        Substructure_Hit_Statistics* tmp = new Substructure_Hit_Statistics(f);
        if (!tmp->create_from_smarts(f)) {
          cerr << "Cannot creat query from smarts  '" << f << "'" << '\n';
          delete tmp;
          return 0;
        }

        // mark the query so that we can marge all of them together later on
        tmp->discard_all_numeric_values();
        tmp->add_numeric_value(i);
        queries.add(tmp);
      } else if (f.starts_with("QUERY:")) {
        f.remove_leading_chars(6);
        if (f.starts_with("F:")) {
          f.remove_leading_chars(2);
          resizable_array_p<Substructure_Hit_Statistics> tmp_queries;
          if (!queries_from_file(f, tmp_queries, 1, verbose)) {
            cerr << "Could not process 'F:" << f << "'\n";
            return 0;
          }
          int tmp_no_queries = tmp_queries.number_elements();
          for (int j = 0; j < tmp_no_queries; j++) {
            Substructure_Hit_Statistics* tmp = tmp_queries[j];
            double tmp_value = 0;
            tmp->numeric_value(tmp_value, 0);
            tmp->discard_all_numeric_values();
            tmp->add_numeric_value(i);
            queries.add(tmp);
          }
          tmp_queries.resize_no_delete(0);
        } else {
          Substructure_Hit_Statistics* tmp = new Substructure_Hit_Statistics(f);
          if (!tmp->read(f)) {
            cerr << "fetch_frag_queries: cannot read query from '" << f << "'\n";
            delete tmp;
            return 0;
          }

          tmp->discard_all_numeric_values();
          tmp->add_numeric_value(i);
          queries.add(tmp);
        }
      } else if (f.starts_with("help")) {
        display_q_directives(0);
      }
    }
  } else {
    cerr << "Nothing to be done.  '-q' option has to be present." << '\n';
    cerr << "The Program has ended" << '\n';
    usage(6);
  }

  int n_queries = queries.number_elements();

  if (n_queries == 0) {
    cerr << "Nothing to be done.  '-q' is empty" << '\n';
    cerr << "The Program has ended" << '\n';
    usage(6);
  }

  double temp = 0.0;
  queries[n_queries - 1]->numeric_value(temp);
  number_of_properties = (int)temp;

  if (cl.option_present('r')) {
    if (!cl.value('r', number_of_regions)) {
      cerr << "Invalid value for number_of_regions (-r).  Has to be a integer between 1 "
              "and 10\n";
      usage(2);
    }

    if ((number_of_regions < 1) || (number_of_regions > 10)) {
      cerr << "The value for number of regions is not valid, should be between 1 and 10 "
              "\n";
      cerr << "10 region is used instead\n";
      number_of_regions = 10;
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', ntest) || ntest < 1) {
      cerr << "The number of tests to perform (-t) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will perform " << ntest << " tests\n";
    }
  }

  if (cl.option_present('d')) {
    align_by_longest_distance = 1;

    if (verbose) {
      cerr << "Will align molecules by longest distance\n";
    }
  }

  if (cl.option_present('O')) {
    const_IWSubstring o = cl.string_value('O');
    stream_for_oriented_molecules.add_output_type(FILE_TYPE_SDF);
    if (!stream_for_oriented_molecules.new_stem(o)) {
      cerr << "Cannot open stream for oriented molecules '" << o << "'\n";
      return 2;
    }
  }

  iwdigits.set_include_leading_space(1);
  iwdigits.initialise(30);

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  set_default_iwstring_float_concatenation_precision(output_precision);

  int rc = 0;
  IWString_and_File_Descriptor output(1);

  output.resize(36000);

  output_result_header(output);

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!jwdist(cl[i], input_type, output)) {
      rc = i + 1;
      break;
    }
  }

  // output info
  if (verbose) {
    if (verbose) {
      cerr << "Read " << molecules_read << " molecules\n";
    }
  }

  if (ntest) {
    cerr << molecules_failing_region_assignment_test << " of " << molecules_read
         << " molecules failed the bin assignment test\n";
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = jwdist(argc, argv);
  return rc;
}
