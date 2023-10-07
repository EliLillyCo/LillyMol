/*
  Tester for standardisation.
  Write molecules and their transformed forms
*/

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

static int verbose = 0;

static int molecules_read = 0;

static int ntest = 10;

static int failures = 0;

static int molecules_with_no_changes = 0;

static int molecules_changed = 0;

static int reduce_to_largest_fragment = 0;

static int transform_to_non_isotopic_form = 0;

static int remove_non_chiral_chiral_centres = 0;

static int remove_all_chiral_centres = 0;

static int remove_all_cis_trans_bonds = 0;

static int remove_all_explicit_hydrogens = 0;

static Chemical_Standardisation chemical_standardisation;

static IWString_and_File_Descriptor stream_for_failures;

static int break_on_error = 0;

static Report_Progress report_progress;

static int abort_on_any_failure = 0;

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
  cerr << "Tester for chemical standardisation\n";
  cerr << "  -n <number>    number of tests per molecule\n";
  cerr << "  -b             stop processing any molecule once it encounters a problem\n";
  cerr << "  -c             remove all chirality\n";
  cerr << "  -x             remove all cis-trans bonds\n";
  cerr << "  -u             remove marked chiral centres that aren't chiral\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -I             transform to non-isotopic form\n";
  cerr << "  -F <fname>     write failing molecules to <fname>\n";
  cerr << "  -r <count>     report progress every <count> molecules\n";
  cerr << "  -f             abort on any failure\n";
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static void
do_remove_non_chiral_chiral_centres(Molecule& m)
{
  int nchiral = m.chiral_centres();
  if (0 == nchiral) {
    return;
  }

  int matoms = m.natoms();

  // cerr << "Has " << nchiral << " chiral centres and " << matoms << " atoms\n";
  int centres_processed = 0;

  int rc = 0;
  for (int i = 0; i < matoms; i++) {
    if (nullptr == m.chiral_centre_at_atom(i)) {
      continue;
    }

    if (!is_actually_chiral(m, i)) {
      m.remove_chiral_centre_at_atom(i);
      rc++;
    }

    centres_processed++;

    if (centres_processed == nchiral) {
      break;
    }
  }

  if (verbose > 1 && rc) {
    cerr << "removed " << rc << " false chiral centres\n";
  }

  return;
}

static void
preprocess(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (transform_to_non_isotopic_form) {
    m.transform_to_non_isotopic_form();
  }

  if (remove_all_explicit_hydrogens) {
    m.remove_all(1);
  }

  if (remove_non_chiral_chiral_centres) {
    do_remove_non_chiral_chiral_centres(m);
  }

  if (remove_all_chiral_centres) {
    m.remove_all_chiral_centres();
  }

  if (remove_all_cis_trans_bonds) {
    m.revert_all_directional_bonds_to_non_directional();
  }

  return;
}

#ifdef NO_LONGER_USED_QWEQJWLEKJH
static void
write_atom_numbers(Molecule& m, std::ostream& output)
{
  const int matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    m.set_isotope(i, i);
  }

  output << m.smiles() << '\n';

  m.transform_to_non_isotopic_form();

  return;
}
#endif

static int
tstandardise(Molecule& m, IWString_and_File_Descriptor& output)
{
  preprocess(m);

  assert(chemical_standardisation.active());

#ifdef ECHO_PI_ELECTRONS
  cerr << "Begin processing '" << m.smiles() << "'\n";
  for (auto i = 0; i < m.natoms(); ++i) {
    int pe;
    (void)m.pi_electrons(i, pe);
    cerr << "Atom " << i << " " << m.smarts_equivalent_for_atom(i) << " has " << pe
         << " pi electrons\n";
  }
#endif

  Molecule mcopy(m);

  if (verbose > 1) {
    cerr << "Reference molecule ";
    write_isotopically_labelled_smiles(mcopy, false, cerr);
    cerr << '\n';
  }

  int starting_structure_changed_by_chemical_standardisation =
      chemical_standardisation.process(m);

  if (verbose > 2) {
    cerr << m.name() << " " << starting_structure_changed_by_chemical_standardisation
         << " changes, standardised to " << m.smiles() << '\n';
  }

  if (0 == starting_structure_changed_by_chemical_standardisation) {
    molecules_with_no_changes++;
  } else {
    molecules_changed++;
  }

  if (!m.valence_ok()) {
    cerr << "Abnormal valence in " << m.smiles() << " " << m.name() << ", skipped\n";
    return 1;
  }

  IWString initial_smiles = m.unique_smiles();

  int failures_this_molecule = 0;

  for (int i = 0; i < ntest; i++) {
    IWString smiles = mcopy.random_smiles();
    Molecule tmp;
    if (!tmp.build_from_smiles(smiles)) {
      cerr << "Yipes, cannot build molecule from random smiles '" << smiles << "'\n";
      continue;
    }

    if (verbose > 1) {
      cerr << "Smiles variant " << i << ' ';
      write_isotopically_labelled_smiles(tmp, false, cerr);
      cerr << '\n';
      //    tmp.debug_print(cerr);
      //    cerr << tmp.smiles() << '\n';
    }

    int changes = chemical_standardisation.process(tmp);

    if (changes != starting_structure_changed_by_chemical_standardisation) {
      cerr << "Potential problems for '" << m.name() << "', parent "
           << starting_structure_changed_by_chemical_standardisation << "\nvariant "
           << changes << ", i = " << i << '\n';
      cerr << "Parent after processing " << initial_smiles << '\n';
      cerr << "Random variant smiles " << smiles << ' ' << tmp.name() << '\n';
    }

    if (tmp.unique_smiles() != initial_smiles) {
      cerr << "Smiles mismatch, '" << m.name() << "'\n";
      cerr << initial_smiles << ' ' << m.name() << " INITIAL\n";
      cerr << smiles << ' ' << m.name() << " RANDOM SMILES\n";
      // Save the unique smiles.
      IWString usmi = tmp.unique_smiles();
      tmp.invalidate_smiles();
      cerr << tmp.smiles() << ' ' << m.name() << " STANDARDISED FROM RANDOM "
           << tmp.name() << '\n';
      cerr << usmi << " unique_smiles\n";
      failures_this_molecule++;

      if (break_on_error) {
        break;
      }
    }

    if (!tmp.valence_ok()) {
      cerr << "Valence error in changed molecule " << m.name() << '\n';
      cerr << initial_smiles << ' ' << m.name() << " INITIAL\n";
      cerr << smiles << ' ' << m.name() << " RANDOM SMILES\n";
      cerr << tmp.unique_smiles() << ' ' << m.name()
           << " STANDARDISED FROM RANDOM - VALENCE ERROR\n";
      failures_this_molecule++;

      if (break_on_error) {
        break;
      }
    }

    //  Now see what happens if we standardise the standardised form

    Molecule tmp2;
    (void)tmp2.build_from_smiles(tmp.smiles());

    if (verbose > 1) {
      cerr << tmp.smiles() << " standardised form\n";
      cerr << "Molecule built from standardised form\n";
      cerr << tmp2.smiles() << '\n';
      write_isotopically_labelled_smiles(tmp2, false, cerr);
      cerr << '\n';
      //    tmp2.debug_print(cerr);
    }

    chemical_standardisation.process(tmp2);

    if (tmp2.unique_smiles() != tmp.unique_smiles()) {
      cerr << "Repeated standardisation failure\n";
      cerr << tmp.unique_smiles() << ' ' << m.name()
           << " initial smiles after standardisation\n";
      cerr << tmp2.unique_smiles() << ' ' << m.name()
           << " built from standardised form and standardised\n";
      tmp.invalidate_smiles();
      tmp2.invalidate_smiles();
      cerr << tmp.smiles() << " and " << tmp2.smiles() << " " << tmp.smiles() << '.'
           << tmp2.smiles() << " build from " << smiles << '\n';
      //    tmp2.debug_print(cerr);
      failures_this_molecule++;
      if (break_on_error) {
        break;
      }
    }
  }

  if (failures_this_molecule) {
    if (stream_for_failures.is_open()) {
      stream_for_failures << mcopy.smiles() << ' ' << m.name() << '\n';
      stream_for_failures.write_if_buffer_holds_more_than(4096);
    }

    failures++;
  }

  return 1;
}

static int
tstandardise(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (!tstandardise(*m, output)) {
      return 0;
    }

    if (report_progress()) {
      cerr << "Processed " << molecules_read << " molecules, " << failures
           << " failures\n";
    }

    if (abort_on_any_failure && failures) {
      return 0;
    }
  }

  return 1;
}

static int
tstandardise(const char* fname, FileType input_type, IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return tstandardise(input, output);
}

static int
tstandardise(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vi:A:E:lIg:n:r:uF:bcxs:f");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!cl.option_present('g')) {
    cerr << "The whole point of this programme is testing standardisation. Use the -g "
            "option\n";
    usage(6);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      usage(6);
    }
  }

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose)) {
      usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('b')) {
    break_on_error = 1;

    if (verbose) {
      cerr << "Will break on encountering any error\n";
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('I')) {
    transform_to_non_isotopic_form = 1;

    if (verbose) {
      cerr << "Will transform to non-isotopic form\n";
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', ntest) || ntest < 1) {
      cerr << "The number of random variants to test (-n) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will test " << ntest << " variants\n";
    }
  }

  if (cl.option_present('s')) {
    random_number_seed_t s;
    if (!cl.value('s', s)) {
      cerr << "Cannot determine random number seed (-s)\n";
      usage(2);
    }
    set_smiles_random_number_seed(s);
  } else {
    set_smiles_random_number_seed_random();
  }

  if (cl.option_present('f')) {
    abort_on_any_failure = 1;

    if (verbose) {
      cerr << "Will terminate processing on any failure\n";
    }
  }

  if (cl.option_present('u')) {
    remove_non_chiral_chiral_centres = 1;

    if (verbose) {
      cerr << "False chiral centres will be removed\n";
    }
  }

  if (cl.option_present('c')) {
    remove_all_chiral_centres = 1;

    if (verbose) {
      cerr << "Will remove all chiral centres\n";
    }
  }

  if (cl.option_present('x')) {
    remove_all_cis_trans_bonds = 1;

    if (verbose) {
      cerr << "Will remove all cis trans bonds\n";
    }
  }

  if (cl.option_present('H')) {
    remove_all_explicit_hydrogens = 1;

    if (verbose) {
      cerr << "Explicit hydrogens removed by default\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise report progress option (-r)\n";
      return 3;
    }
  }

  set_copy_name_in_molecule_copy_constructor(1);

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(3);
  }

  if (cl.option_present('F')) {
    IWString f = cl.option_value('F');

    if (!f.ends_with(".smi")) {
      f << ".smi";
    }

    if (!stream_for_failures.open(f.null_terminated_chars())) {
      cerr << "Cannot open stream for failed molecules '" << f << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Failures writtent to '" << f << "'\n";
    }
  }

  set_display_abnormal_valence_messages(0);

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!tstandardise(cl[i], input_type, output)) {
      cerr << "Error processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  if (verbose || failures > 0) {
    cerr << "Processed " << molecules_read << " molecules\n";
    cerr << molecules_changed << " molecules changed, " << molecules_with_no_changes
         << " unchanged\n";
    cerr << failures << " molecules with failures\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  int rc = tstandardise(argc, argv);

  return rc;
}
