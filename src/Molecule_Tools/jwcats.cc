/*
 * Compute the CATS descriptors according to the following paper
 * Angew. Chem. Int. Ed. 1999, 38, No. 19 by Gisbert Schneider, Werner Neidhart, Thomas
 * Giller, and Gerard Schmid
 */

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <memory>

using std::cerr;
using std::endl;
using std::ostream;

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/allowed_elements.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

// DEBUG_SWITCH is used to print out debugging msgs
#ifndef DEBUG_SWITCH
#define DEBUG_SWITCH 0
#endif

static Donor_Acceptor_Assigner donor_acceptor_assigner;
static Charge_Assigner charge_assigner;

static Chemical_Standardisation chemical_standardisation;

static Allowed_Elements allowed_elements;

static int verbose = 0;

static int molecules_read = 0;

static int reduce_to_largest_fragment = 0;

static int molecules_containing_non_alowed_elements = 0;

static int include_hydrophobic_pairs = 1;

/*
  Output is complicated by the presence or absence of hydrophobic pairs.
  We can make things a little easier by pre-computing an array of those
  values that will be output.
*/

static int* write_array_value = nullptr;

static int min_bond_separation = 0;

static int max_bond_separation = 10;

static int array_size = max_bond_separation * 15;

static int use_queries_to_determine_hydrophobicity = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static int scaling_type = 1;

static int make_implicit_hydrogens_explicit = 0;

// keep track of number of warning and number of error

static int number_of_error = 0;

static IWString fingerprint_tag;
static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int function_as_gfp_filter = 0;

static Accumulator_Int<int> acceptor_acc, donor_acc, negative_acc, positive_acc,
    hydrophobe_acc, features_acc;

static int flush_after_every_molecule = 0;

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

  cerr << "Calculate CAT Search descriptors for molecules\n";

  cerr << "  -p             do not output values related to hydrophobic-hydrophobic pair\n";
  cerr << "  -m <n>         goes up to length of bonds(default 10)\n";
  cerr << "  -z <n>         minimum bond separation (features closer than <n> bonds are ignored\n";
  cerr << "  -q             use queries to define lipophilic agroups\n";
  cerr << "  -s             scaling type\n";
  cerr << "                 0 => no scaling\n";
  cerr << "                 1 => normalized by the number of heavy atom (default)\n";
  cerr << "                 2 => normalized by the sum of the number of two pharmacophore types\n";
  cerr << "                 3 => normalized by the atoms / (sum of the number of two pharmacophore types)\n";
  cerr << "  -J <tag>       create fingerprints\n";
  cerr << "  -H ...         donor acceptor assignment, enter '-H help' for info\n";
  cerr << "  -N ...         charge assigner, enter '-N help' for info\n";
  cerr << "  -h             make implicit hydrogens explicit\n";
  cerr << "  -f             function as a tdt filter\n";
  cerr << "  -X ...         more options\n";
  cerr << "  -i <type>      input type\n";
  cerr << "  -A ...         aromaticity options\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
preprocess_molecule(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation
          .active()) {  // just assume this is getting rid of explicit hydrogens
    chemical_standardisation.process(m);
  } else {
    m.remove_all(1);  // don't want explicit hydrogens here
  }

  return 1;
}

static int
do_fingerprint_output(Molecule& m, const int array_size, const double* scaled_counts,
                      IWString_and_File_Descriptor& output)
{
  Sparse_Fingerprint_Creator sfc;

  for (int i = 0; i < array_size; i++) {
    if (0.0 == scaled_counts[i]) {
      continue;
    }

    if (!write_array_value[i]) {
      continue;
    }

    int j = static_cast<int>(scaled_counts[i] + 0.01);

    sfc.hit_bit(i, j);
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tmp);

  output << fingerprint_tag << tmp << ">\n";

  if (function_as_gfp_filter) {
    ;
  } else {
    output << "|\n";
  }

  return 1;
}

/*
  We give out meaningful descriptor names. Very important that this array
  be kept in sync with what find_property_pair_number returns
*/

// clang-format off
static const char * pair_name[] = {
  "AA",
  "AD",
  "AP",
  "AN",
  "AH",
  "DD",
  "DP",
  "DN",
  "DH",
  "PP",
  "PN",
  "PH",
  "NN",
  "NH",
  "HH" };

// clang-format on

static int
find_property_pair_number(int npi, int npj)
{
  assert((npi > -1) && (npi < 5) && (npj > -1) && (npj < 5));

  int small_no = npi;
  int large_no = npj;

  if (small_no > large_no) {
    small_no = npj;
    large_no = npi;
  }
  return (small_no * 5 - small_no * (small_no - 1) / 2 + large_no - small_no);
}

static int
output_result_header(IWString_and_File_Descriptor& output)
{
  int size = 15;
  // nomally, you print out the whole thing, but, when no hydrophobic-pair flag is on, do
  // not print out them
  if (!include_hydrophobic_pairs) {
    size = 14;
  }

  output << "Name";

  for (int i = 0; i < max_bond_separation; i++) {
    for (int j = 0; j < size; j++) {
      //    output <<" jwc_B"<< (i + 1) <<"P"<<j;
      output << " jwc_B" << (i + 1) << "P" << pair_name[j];
    }
  }

  output << '\n';

  return 1;
}

static int
handle_missing_charge_data_fingerprint(Molecule& m, int array_size,
                                       IWString_and_File_Descriptor& output)
{
  output << fingerprint_tag << ">\n";

  molecules_containing_non_alowed_elements++;

  return 1;
}

static int
handle_missing_charge_data(Molecule& m, int array_size,
                           IWString_and_File_Descriptor& output)
{
  if (fingerprint_tag.length() > 0) {
    return handle_missing_charge_data_fingerprint(m, array_size, output);
  }

  const_IWSubstring mname = m.name();
  mname.truncate_at_first(' ');

  output << mname;

  for (int i = 0; i < array_size; i++) {
    if (write_array_value[i]) {
      output << " .";
    }
  }

  output << '\n';

  return 1;
}

/*
 * This function is to calculate the CAT Search descriptors for molecules
 */

static int
jw_cat_search(Molecule& m, IWString_and_File_Descriptor& output, int donor_acceptor[],
              int** properties, double scaled_counts[])
{
  int property_count[5];
  set_vector(property_count, 5, 0);

  const int n_atoms = m.natoms();

  // initialize vectors
  for (int i = 0; i < 5; i++) {
    set_vector(properties[i], n_atoms, 0);
  }

  set_vector(scaled_counts, array_size, 0.0);

  if (donor_acceptor_assigner.active()) {
    donor_acceptor_assigner.process(m, donor_acceptor);
  }

#ifdef ECHO_DONOR_ACCEPTOR_RESULTS
  for (int i = 0; i < n_atoms; i++) {
    cerr << m.isotope(i) << endl;
  }
#endif

  int number_heavy_atoms = 0;

  for (int i = 0; i < n_atoms; i++) {
    const Atom* atom_i = m.atomi(i);

    if (1 == atom_i->atomic_number()) {
      continue;
    }

    number_heavy_atoms++;

    switch (donor_acceptor[i]) {
      case 0:
        break;

      case 1:
        properties[0][i] = 1;
        break;

      case 2:
        properties[0][i] = 1;
        properties[1][i] = 1;
        break;

      case 3:
        properties[1][i] = 1;
        break;

      default:
        break;
    }

    // positve charge, negative charge

    formal_charge_t fci = atom_i->formal_charge();

    if (0 == fci) {
      ;
    } else if (1 == fci) {
      properties[2][i] = 1;
    } else if (-1 == fci) {
      properties[3][i] = 1;
    }
  }

  if (use_queries_to_determine_hydrophobicity) {
    int nq = queries.number_elements();
    for (int i = 0; i < nq; i++) {
      Substructure_Results sresults;

      int nhits = queries[i]->substructure_search(m, sresults);
      if (verbose > 1) {
        cerr << ' ' << nhits << " hits to query " << i << endl;
      }

      if (0 == nhits) {
        continue;
      }

      for (int j = 0; j < nhits; j++) {
        const Set_of_Atoms* e = sresults.embedding(j);
        e->set_vector(properties[4], 1);
        /*          int n = e->number_elements();
                    for (int k=0; k<n; k++)
                    {
                      properties[4][e->item(k)] = 1;
                    }*/
      }
    }
    //      cerr<<"Hydrophobic"<<endl;
    //      for (int i=0; i<n_atoms; i++)
    //        if (properties[4][i]) cerr<<"atomi="<<i<<endl;
  } else {
    // compute partial atomic charge -- only needed when looking at lipophilicity
    // descriptor

    if (!m.compute_Gasteiger_partial_charges()) {
      if (!function_as_gfp_filter) {
        cerr << m.name() << " ERROR in calculation of partial charge, '"
             << m.name() << "'\n";
      }
      return handle_missing_charge_data(m, array_size, output);
      ;
    }

    for (int i = 0; i < n_atoms; i++) {
      //    cerr << i << ' ' << m.charge_on_atom(i) << endl;
      if (fabs(m.charge_on_atom(i)) <= 0.2) {
        properties[4][i] = 1;
        if (verbose) {
          hydrophobe_acc.extra(1);
        }
      }
      //    cerr << "Atom " << i << " type " << m.atomic_number(i) << " q " <<
      //    m.charge_on_atom(i) << " p " << properties[4][i] << endl;
    }
  }

  // calculate the property count
  if (scaling_type) {
    for (int i = 0; i < n_atoms; i++) {
      for (int j = 0; j < 5; j++) {
        //      cerr << i << ' ' << j << ' ' << properties[j][i] << endl;
        if (properties[j][i]) {
          property_count[j]++;
        }
      }
    }
  }

  if (verbose) {
    int acceptor_this_molecule = 0;
    int donor_this_molecule = 0;
    int positive_this_molecule = 0;
    int negative_this_molecule = 0;
    int hydrophobe_this_molecule = 0;

    for (int i = 0; i < n_atoms; i++) {
      if (properties[0][i]) {
        acceptor_this_molecule++;
      }
      if (properties[1][i]) {
        donor_this_molecule++;
      }
      if (properties[2][i]) {
        positive_this_molecule++;
      }
      if (properties[3][i]) {
        negative_this_molecule++;
      }
      if (properties[4][i]) {
        hydrophobe_this_molecule++;
      }
    }

    acceptor_acc.extra(acceptor_this_molecule);
    donor_acc.extra(donor_this_molecule);
    positive_acc.extra(positive_this_molecule);
    negative_acc.extra(negative_this_molecule);
    hydrophobe_acc.extra(hydrophobe_this_molecule);

    int features_this_molecule = acceptor_this_molecule + donor_this_molecule +
                                 negative_this_molecule + positive_this_molecule +
                                 hydrophobe_this_molecule;
    features_acc.extra(features_this_molecule);
  }

  // for each atom pair, find out if they belong to a properties pair,
  // if yes, increment the count

  // cerr << "HAC " << n_atoms << endl;

  for (int i = 0; i < n_atoms; i++) {
    if (1 == m.atomic_number(i)) {
      continue;
    }

    for (int j = i + 1; j < n_atoms; j++) {
      if (1 == m.atomic_number(j)) {
        continue;
      }

      int bond_distance = m.bonds_between(i, j);

      if (bond_distance > max_bond_separation) {
        continue;
      }

      if (bond_distance < min_bond_separation) {
        continue;
      }

      for (int npi = 0; npi < 5; npi++) {
        if (0 == properties[npi][i]) {
          continue;
        }

        for (int npj = 0; npj < 5; npj++) {
          if (0 == properties[npj][j]) {
            continue;
          }

          int property_pair_no = find_property_pair_number(npi, npj);

          int array_column = (bond_distance - 1) * 15 + property_pair_no;
          //        cerr << " i " << i << " j " << j << " dist " << bond_distance << " npi
          //        " << npi << " npj " << npj << " pair " << property_pair_no << " col "
          //        << array_column << endl;
          if (0 == scaling_type) {
            scaled_counts[array_column]++;
          } else if (1 == scaling_type) {
            scaled_counts[array_column] += 1.0 / static_cast<double>(number_heavy_atoms);
          } else if (2 == scaling_type) {
            scaled_counts[array_column] +=
                1.0 / (property_count[npi] + property_count[npj]);
          } else if (3 == scaling_type) {
            scaled_counts[array_column] += static_cast<double>(number_heavy_atoms) /
                                           (property_count[npi] + property_count[npj]);
          } else {
            cerr << "Not sure what to do with scaling type " << scaling_type << endl;
            return 0;
          }
          //        cerr << "put " << npi << " and " << npj << " into " <<
          //        property_pair_no << " value " << scaled_counts[array_column] << endl;
        }
      }
    }
  }

  if (fingerprint_tag.length() > 0) {
    do_fingerprint_output(m, array_size, scaled_counts, output);

    if (function_as_gfp_filter) {
      ;
    } else {
      output << "|\n";
    }

    if (flush_after_every_molecule) {
      output.flush();
    }

    return 1;
  }

  // Descriptor file output

  append_first_token_of_name(m.name(), output);

  for (int i = 0; i < array_size; i++) {
    if (!write_array_value[i]) {
      continue;
    }

    if (0.0 == scaled_counts[i]) {
      output << " 0";
    } else {
      output << ' ' << static_cast<float>(scaled_counts[i]);
    }
  }

  output << '\n';

  return 1;
}

/*
  This is a little strange.
  Ideally, I should figure out if we can process the molecule first. If not,
  don't write anything.
  But for multi-fragment processing, I want to write the smiles first
  before I preprocess it. That is why this is done in such a seemingly
  strange way.
  And besides, people should not be putting in molecules with disallowed
  atoms anyway.
*/

static int
jw_cat_search(Molecule& m, IWString_and_File_Descriptor& output)
{
  molecules_read++;

  if (function_as_gfp_filter) {
    ;
  } else if (fingerprint_tag.length())  // write smiles before molecule gets changed
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  preprocess_molecule(m);

  if (allowed_elements.contains_non_allowed_atoms(m)) {
    molecules_containing_non_alowed_elements++;

    if (fingerprint_tag.empty()) {  // descriptors, ignore molecule
      return 1;
    }

    output << fingerprint_tag << ">\n";

    if (!function_as_gfp_filter) {
      output << "|\n";
    }

    return 1;
  }

  m.recompute_distance_matrix();

  //    cerr << "Before charge assigner, Number of Atom"<<n_atoms<<endl;

  if (charge_assigner.active()) {
    charge_assigner.process(m);
  }

  //    cerr << "After charge assigner,  Number of Atom"<<m.natoms()<<endl;

  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  int n_atoms = m.natoms();

  int* donor_acceptor = new int[n_atoms];
  std::unique_ptr<int[]> free_donor_acceptor(donor_acceptor);

  int* properties[5];
  for (int i = 0; i < 5; i++) {
    properties[i] = new int[n_atoms];
  }

  int array_size = 15 * max_bond_separation;

  double* scaled_counts = new double[array_size];
  std::unique_ptr<double[]> free_scaled_count(scaled_counts);

  int rc = jw_cat_search(m, output, donor_acceptor, properties, scaled_counts);

  if (0 == rc) {
    number_of_error++;
  }

  for (int i = 0; i < 5; i++) {
    delete[] properties[i];
  }

  return rc;
}

static int
jw_cat_search(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output)
{
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (verbose > 1) {
      cerr << molecules_read << " processing '" << m->name() << "'\n";
    }

    if (!jw_cat_search(*m, output)) {
      cerr << "ERROR in computing descriptors for " << m->name() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
jw_cat_search_record(const const_IWSubstring& buffer,
                     IWString_and_File_Descriptor& output)
{
  const_IWSubstring mybuffer(buffer);
  mybuffer.remove_leading_chars(smiles_tag.length());
  mybuffer.chop();

  // cerr << "Processing '" << buffer << "'\n";

  Molecule m;

  if (!m.build_from_smiles(mybuffer)) {
    cerr << "Invalid smiles '" << mybuffer << "'\n";
    return 0;
  }

  return jw_cat_search(m, output);
}

static int
jw_cat_search(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(32768);

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!jw_cat_search_record(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
jw_cat_search(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return jw_cat_search(input, output);
}

static void
report_accumulator(const Accumulator_Int<int>& acc, const char* atype, std::ostream& os)
{
  if (0 == acc.n()) {
    os << "No molecules had " << atype << " features\n";
    return;
  }

  os << "molecules had between " << acc.minval() << " and " << acc.maxval() << " "
     << atype << " ave " << static_cast<float>(acc.average_if_available_minval_if_not())
     << '\n';

  return;
}

static int
jw_cat_search(const char* fname, FileType input_type,
              IWString_and_File_Descriptor& output)
{
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 2) {
    input.set_verbose(1);
  }

  return jw_cat_search(input, output);
}

static void
DisplayDashXOptions(std::ostream& output)
{
  output << " -X flush      flush output after each molecule\n";

  ::exit(0);
}

static int
jw_cat_search(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:pi:H:q:N:F:m:s:g:J:K:flhz:X:");

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

  if (cl.option_count('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      usage(4);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will recuce to largest fragment\n";
    }
  }

  if (cl.option_present('p')) {
    include_hydrophobic_pairs = 0;

    if (verbose) {
      cerr << "Will xxx hydrophobe-hydrophobe pairs\n";
    }
  }

  if (cl.option_present('s')) {
    int n = cl.option_count('s');

    if (!cl.value('s', scaling_type, n - 1)) {
      cerr << "Invalid value for scaling type (-s) has to be a number\n";
      usage(2);
    }

    if ((scaling_type < 0) || (scaling_type > 3)) {
      cerr << "The value for scaling type is not valid, should be 0, 1 or 2\n";
      cerr << "1 (normalized by number of heavy atoms) is used instead\n";
      scaling_type = 1;
    }
  }

  if (cl.option_present('h')) {
    make_implicit_hydrogens_explicit = 1;

    if (verbose) {
      cerr << "Will make implicit Hydrogens explicit\n";
    }
  }

  if (cl.option_present('z')) {
    if (!cl.value('z', min_bond_separation) || min_bond_separation < 1) {
      cerr << "Minimum bond separation (-z) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Features closer than " << min_bond_separation << " bonds apart ignored\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_bond_separation) || max_bond_separation < 1) {
      cerr << "Maximum bond separation (-m) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will perceive features up to " << max_bond_separation << " bonds apart\n";
    }
  }

  array_size = 15 * max_bond_separation;

  write_array_value = new_int(array_size, 1);

  if (!include_hydrophobic_pairs) {
    for (int i = 0; i < array_size; i++) {
      if (14 == i % 15) {
        write_array_value[i] = 0;
      }
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', fingerprint_tag);

    if (verbose) {
      cerr << "Will produce fingerprints with tag '" << fingerprint_tag << "'\n";
    }

    if (!fingerprint_tag.ends_with('<')) {
      fingerprint_tag.add('<');
    }

    scaling_type = 0;

    if (cl.option_present('f')) {
      function_as_gfp_filter = 1;

      if (verbose) {
        cerr << "Will function as a filter\n";
      }
    }
  }

  if (cl.option_present('q')) {
    queries.resize(cl.option_count('q') + 100);
    if (!process_queries(cl, queries, verbose)) {
      cerr << "Cannot process queries from -q option(s)\n";
      return 6;
    }
    use_queries_to_determine_hydrophobicity = 1;
  }

  int nq = queries.number_elements();
  for (int i = 0; i < nq; i++) {
    queries[i]->set_find_unique_embeddings_only(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (function_as_gfp_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('H')) {
    if (!donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose)) {
      cerr << "Cannot initialise donor/acceptor assignment object\n";
      usage(4);
    }
  }

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, 0, 'N')) {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage(1);
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_after_every_molecule = 1;
        if (verbose) {
          cerr << "Will flush after every molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  set_default_iwstring_float_concatenation_precision(3);

  std::ofstream logfile;
  if (verbose) {
    logfile.open("jwcatsearch_descriptor.log", std::ios::out);
  }

  if (!logfile) {
    cerr << "jwcatsearch_descriptor.log file cannot be opened\n";
    return 0;
  }

  time_t current_time;

  if (verbose) {
    logfile << "This file collect all the info about error during the calculation\n";
    logfile << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
    time(&current_time);
    logfile << "calculation started at " << ctime(&current_time);
    logfile << "The complete Command was:" << endl;
    for (int kk = 0; kk < argc; kk++) {
      logfile << argv[kk] << " ";
    }
    logfile << endl;
  }

  IWString_and_File_Descriptor output(1);

  if (0 == fingerprint_tag.length()) {
    output_result_header(output);
  }

#ifdef DEBUGGING_THING_TO_DETERMINE_PAIRS
  for (int i = 0; i < 5; i++) {
    for (int j = i; j < 5; j++) {
      int k = find_property_pair_number(i, j);
      cerr << " i = " << i << " j = " << j << " pair " << k << endl;
    }
  }
#endif

  int rc = 0;

  if (function_as_gfp_filter) {
    if (cl.number_elements() > 1) {
      cerr << "Fingerprint filter cannot have multiple arguments\n";
      exit(3);
    }

    rc = jw_cat_search(cl[0], output);
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!jw_cat_search(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  // output info
  if (verbose) {
    time(&current_time);
    logfile << "\n\ncalculation ended at " << ctime(&current_time) << endl;
    logfile << "Total Molecules read: " << molecules_read
            << "\tTotal Error: " << number_of_error << endl;
    logfile << "Error Rate: " << (double)number_of_error / (double)molecules_read << endl;
    cerr << "Read " << molecules_read << " molecules\n";

    cerr << "Processed " << molecules_read << " molecules\n";
    report_accumulator(acceptor_acc, "acceptor", cerr);
    report_accumulator(donor_acc, "donor", cerr);
    report_accumulator(positive_acc, "positive", cerr);
    report_accumulator(negative_acc, "negative", cerr);
    report_accumulator(hydrophobe_acc, "hydrophobe", cerr);
    report_accumulator(features_acc, "features", cerr);

    if (molecules_containing_non_alowed_elements) {
      cerr << molecules_containing_non_alowed_elements
           << " molecules with non-allowed atom types\n";
    }
  }

  if (nullptr != write_array_value) {
    delete[] write_array_value;
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = jw_cat_search(argc, argv);
  return rc;
}
