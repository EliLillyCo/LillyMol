/*
  Interface to corina
*/

#include <setjmp.h>

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/string_data_source/iwstring_string_data_source.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "corina_lib.h"

#include <sys/signal.h>

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_converted_first_time = 0;

static int molecules_failing_but_no_chiral_centres = 0;

static int molecules_converted_after_chiral_centre_removal = 0;

static int molecules_converted_after_timeout = 0;

static int molecules_failing_for_timeout = 0;

static int last_molecule_failed_for_timeout = 0;

static int current_molecule_failed_for_timeout = 0;

static int molecules_not_converted = 0;

static int computation_not_attempted = 0;

static int random_smiles_to_try = 0;

static int molecules_converted_after_random_smiles = 0;

static Chemical_Standardisation input_chemical_standardisation;

static Chemical_Standardisation from_corina_chemical_standardisation;

static int remove_hydrogens_from_molecules_from_corina = 0;

static int reduce_to_largest_fragment = 1;

static int remove_chiral_centres = 0;

static IWString dash_d_value;

IWString_and_File_Descriptor stream_for_failing_molecules;

static int pass_unique_smiles_to_corina = 1;
static int pass_non_aromatic_unique_smiles_to_corina = 0;

static int generate_3d_smiles = 0;

static int tdt_output = 0;

/*
  Sometimes we want a quick way of identifying those molecules
  that required chiral centre removal or random smiles
*/

IWString_and_File_Descriptor stream_for_problematic_molecules;

static int fix_seed_before_random_smiles = 0;
static unsigned long random_number_seed = 0;

/*
  We can write the .sdf information to TDT output stream
*/

static int function_as_tdt_filter = 0;

static int pass_smiles_as_string_to_corina = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int timeout = 0;

static int max_rings_to_consider = 20;   // numeric_limits<int>::max();
static int max_atoms_to_consider = 100;  // numeric_limits<int>::max();

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
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Interface to corina\n";
  cerr << "  -d ...        options passed to corina\n";
  cerr << "  -x            if a 3D structure is not generated, remove chiral centres\n";
  cerr << "  -r <number>   number of random smiles to try\n";
  cerr << "  -S <fname>    standardise output from corina and write to <fname> in SDF format\n";
  cerr << "  -h            also remove explicit Hydrogens in generated molecule (-S only)\n";
  cerr << "  -F <fname>    output file name for failed molecules - smiles only\n";
  cerr << "  -P <fname>    output file name for problematic molecules - smiles only\n";
  cerr << "  -e <seed>     specify seed for random smiles\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -q            pass normal smiles (not unique smiles) to corina\n";
  cerr << "  -n            pass non-aromatic unique smiles to corina\n";
  cerr << "  -f            function as a tdt filter\n";
  cerr << "  -t            produce TDT ouput\n";
  cerr << "  -u            produce smiles with coordinates\n";
  cerr << "  -s            input is already smiles, pass smiles strings from file to corina\n";
  cerr << "  -R <nrings>   do not process any molecule having more than <nrings> rings (def " << max_rings_to_consider << ")\n";
  cerr << "  -C <natoms>   do not process any molecule having more than <natoms> atoms (def " << max_atoms_to_consider << ")\n";
  cerr << "  -y <seconds>  abandon any conversion after <seconds> seconds\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -o ...        options for corina's -o option (NOT the usual output types for Lilly tools)\n";
  cerr << "                    for example, to get pdb output '-o t=pdb'\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static jmp_buf env;

static void
signal_handler(int zsig) {
  if (verbose) {
    cerr << "Signal handler caught singal type " << zsig << '\n';
  }

  siglongjmp(env, zsig);

  return;
}

static void
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  // cerr << "Reduced to largest frag, now '" << m.smiles() << "'\n";

  if (input_chemical_standardisation.active()) {
    input_chemical_standardisation.process(m);
  }

  return;
}

static void
free_corina_storage() {
  corina_buffer(NULL, NULL, NULL, NULL, NULL);
}

static int
handle_failed_conversion(Molecule& m) {
  molecules_not_converted++;

  if (current_molecule_failed_for_timeout) {
    molecules_failing_for_timeout++;
  }

  if (verbose > 1) {
    cerr << "Could not convert '" << m.name() << "'\n";
  }

  if (stream_for_failing_molecules.is_open()) {
    m.invalidate_smiles();  // don't want unique smiles

    stream_for_failing_molecules << m.smiles() << ' ' << m.name() << '\n';

    stream_for_failing_molecules.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
write_to_problematic_molecule_stream_if_open(Molecule& m,
                                             const char* additional_info = NULL) {
  if (!stream_for_problematic_molecules.is_open()) {
    return 1;
  }

  stream_for_problematic_molecules << m.smiles() << ' ' << m.name();

  if (NULL != additional_info) {
    stream_for_problematic_molecules << ' ' << additional_info;
  }

  stream_for_problematic_molecules << '\n';

  stream_for_problematic_molecules.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
molecule_from_corina_output(const char* from_corina, Molecule& m) {
  String_Data_Source inp(from_corina);

  if (!m.read_molecule_mdl_ds(inp)) {
    cerr << "Cannot parse output from corina\n";
    cerr << from_corina << '\n';
    abort();
    return 0;
  }

  return 1;
}

static int
restandardise(const char* from_corina,
              Chemical_Standardisation& from_corina_chemical_standardisation,
              Molecule_Output_Object& molecule_output) {
  Molecule m;

  if (!molecule_from_corina_output(from_corina, m)) {
    return 0;
  }

  from_corina_chemical_standardisation.process(m);

  if (remove_hydrogens_from_molecules_from_corina) {
    m.remove_all(1);
  }

  return molecule_output.write(m);
}

static int
write_3d_smiles(Molecule& m, std::ostream& text_output) {
  set_append_coordinates_after_each_atom(1);
  m.invalidate_smiles();
  text_output << m.smiles();
  set_append_coordinates_after_each_atom(0);

  free_corina_storage();

  return 1;
}

static int
write_3d_smiles(const char* from_corina,
                Chemical_Standardisation& from_corina_chemical_standardisation,
                std::ostream& text_output) {
  Molecule m;

  if (!molecule_from_corina_output(from_corina, m)) {
    return 0;
  }

  from_corina_chemical_standardisation.process(m);

  return write_3d_smiles(m, text_output);
}

static int
do_output(const IWString& mname, const char* from_corina,
          Chemical_Standardisation& from_corina_chemical_standardisation,
          Molecule_Output_Object& molecule_output, std::ostream& text_output) {
  if (!from_corina_chemical_standardisation
           .active())  // nothing molecular happening here, just write what corina gave us
  {
    text_output << from_corina;
    free_corina_storage();
    return 1;
  }

  Molecule m;

  if (!molecule_from_corina_output(from_corina, m)) {
    return 0;
  }

  from_corina_chemical_standardisation.process(m);

  if (tdt_output) {
    text_output << smiles_tag;
    write_3d_smiles(m, text_output);
    text_output << ">\n";
    text_output << identifier_tag << mname << ">\n";
    text_output << "|\n";
  } else if (generate_3d_smiles) {
    write_3d_smiles(m, text_output);
    text_output << ' ' << mname << '\n';
  } else {
    restandardise(from_corina, from_corina_chemical_standardisation, molecule_output);
  }

  free_corina_storage();

  return 1;
}

static int
write_3d_smiles_zero_coordinates(Molecule& m, std::ostream& text_output) {
  set_append_coordinates_after_each_atom(1);
  m.invalidate_smiles();
  text_output << smiles_tag << m.smiles() << ">\n";
  set_append_coordinates_after_each_atom(0);

  return 1;
}

class Corina_Buffers  // might be a good idea to implement sometime
{
 private:
  char* _out;
  char* _err;

 public:
};

static char* out = NULL;
static char* err = NULL;

static int
rcorina34(const char* corina_input) {
// #define DEBUG_CORINA_BUFFER
#ifdef DEBUG_CORINA_BUFFER
  cerr << "Sending '" << dash_d_value << "'\n";
  cerr << "Input '" << corina_input << "'\n";
#endif

  out = NULL;
  err = NULL;
  unsigned long stats[4];

  corina_buffer(const_cast<char*>(dash_d_value.null_terminated_chars()),
                const_cast<char*>(corina_input), &out, &err, stats);

  if (NULL == out || '\0' == out[0])  // failure
  {
    free_corina_storage();
    return 0;
  }

  return 1;
}

static int
rcorina34_alarm(const char* corina_input) {
  assert(timeout > 0);

  alarm(timeout);

  int returned_from_longjump = sigsetjmp(env, 1);

#ifdef DEBUG_RCORINA_ALARM
  cerr << returned_from_longjump << " returned_from_longjump\n";
#endif

  if (0 == returned_from_longjump) {  // normal, not returning from signal handler
    return rcorina34(corina_input);
  }

  if (SIGALRM == returned_from_longjump) {
    current_molecule_failed_for_timeout++;
  } else {
    cerr << "Unrecognised rc '" << returned_from_longjump << " from longjump\n";
  }

  return 0;
}

static int
rcorina34(Molecule& m) {
  IWString corina_input;
  corina_input.resize(512);

  if (pass_unique_smiles_to_corina) {
    corina_input << m.unique_smiles();
  } else if (pass_non_aromatic_unique_smiles_to_corina) {
    corina_input << m.non_aromatic_unique_smiles();
  } else {
    corina_input << m.smiles();
  }

  corina_input << ' ' << m.name() << '\n';

  if (timeout > 0) {
    return rcorina34_alarm(corina_input.null_terminated_chars());
  } else {
    return rcorina34(corina_input.null_terminated_chars());
  }
}

static int
try_removing_chiral_centres(Molecule& m, const Set_of_Atoms& s) {
  for (int i = 0; i < s.number_elements(); i++) {
    atom_number_t a = s[i];

    m.remove_chiral_centre_at_atom(a);

    if (rcorina34(m)) {
      if (verbose > 1) {
        cerr << m.name() << " converted after " << (i + 1) << " chiral centres removed\n";
      }
      return write_to_problematic_molecule_stream_if_open(m, "RMCHIRAL");
    }
  }

  return 0;
}

static int
try_random_smiles(Molecule& m) {
  if (fix_seed_before_random_smiles) {
    set_smiles_random_number_seed(random_number_seed);
  }

  for (int i = 0; i < random_smiles_to_try; i++) {
    m.invalidate_smiles();

    IWString corina_input;

    corina_input << m.random_smiles() << ' ' << m.name() << '\n';

    //  cerr << "Trying '" << corina_input << "'\n";

    if (rcorina34(corina_input.null_terminated_chars())) {
      if (verbose > 1) {
        cerr << m.name() << " converted after " << (i + 1) << " random smiles\n";
      }
      return write_to_problematic_molecule_stream_if_open(m, "RANDOM");
    }
  }

  return 0;
}

static int
identify_chiral_centres(Molecule& m, Set_of_Atoms& chiral_centres_in_rings,
                        Set_of_Atoms& chiral_centres_in_chains) {
  int nc = m.chiral_centres();

  for (int i = 0; i < nc; i++) {
    const Chiral_Centre* c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    atom_number_t a = c->a();

    if (m.is_ring_atom(a)) {
      chiral_centres_in_rings.add(a);
    } else {
      chiral_centres_in_chains.add(a);
    }
  }

  return nc;
}

static int
rcorina34(Molecule& m, Molecule_Output_Object& molecule_output,
          std::ostream& text_output) {
  if (m.natoms() > max_atoms_to_consider || m.nrings() > max_rings_to_consider) {
    computation_not_attempted++;
    return 1;
  }

  current_molecule_failed_for_timeout = 0;

  if (rcorina34(m)) {
    molecules_converted_first_time++;

    return do_output(m.name(), out, from_corina_chemical_standardisation, molecule_output,
                     text_output);
  }

#ifdef DEBUG_RCORINA_ALARM
  cerr << current_molecule_failed_for_timeout << " current_molecule_failed_for_timeout\n";
#endif

  if (!remove_chiral_centres) {
    return handle_failed_conversion(m);
  }

  if (current_molecule_failed_for_timeout) {
    if (0 == m.chiral_centres()) {
      return handle_failed_conversion(m);
    }

    m.remove_all_chiral_centres();

    if (rcorina34(m)) {
      molecules_converted_after_timeout++;
      return do_output(m.name(), out, from_corina_chemical_standardisation,
                       molecule_output, text_output);
    }

    return handle_failed_conversion(m);
  }

  Set_of_Atoms chiral_centres_in_rings, chiral_centres_in_chains;

  int nc = identify_chiral_centres(m, chiral_centres_in_rings, chiral_centres_in_chains);

  if (0 == nc && 0 == remove_chiral_centres) {
    molecules_failing_but_no_chiral_centres++;

    if (verbose > 1) {
      cerr << m.name() << " failed, no chiral centres\n";
    }

    return handle_failed_conversion(m);
  }

  // First try removing chiral centres in rings

  if (chiral_centres_in_rings.number_elements() &&
      try_removing_chiral_centres(m, chiral_centres_in_rings)) {
    molecules_converted_after_chiral_centre_removal++;

    return do_output(m.name(), out, from_corina_chemical_standardisation, molecule_output,
                     text_output);
  }

  if (chiral_centres_in_chains.number_elements() &&
      try_removing_chiral_centres(m, chiral_centres_in_chains)) {
    molecules_converted_after_chiral_centre_removal++;

    return do_output(m.name(), out, from_corina_chemical_standardisation, molecule_output,
                     text_output);
  }

  if (random_smiles_to_try && try_random_smiles(m)) {
    molecules_converted_after_random_smiles++;

    return do_output(m.name(), out, from_corina_chemical_standardisation, molecule_output,
                     text_output);
  }

  return handle_failed_conversion(m);
}

static int
rcorina34_filter_process(IWString& smiles, std::ostream& text_output) {
  if (rcorina34(smiles.null_terminated_chars())) {
    molecules_converted_first_time++;

    return write_3d_smiles(out, from_corina_chemical_standardisation, text_output);
  }

  Molecule m;
  if (!m.build_from_smiles(smiles)) {
    cerr << "Yipes, cannot interpret smiles of failed conversion '" << smiles << "'\n";
    return 0;
  }

  Set_of_Atoms chiral_centres_in_rings, chiral_centres_in_chains;

  int nc = identify_chiral_centres(m, chiral_centres_in_rings, chiral_centres_in_chains);

  if (0 == nc && 0 == remove_chiral_centres) {
    molecules_failing_but_no_chiral_centres++;

    if (verbose > 1) {
      cerr << m.name() << " failed, no chiral centres, zero coordinates\n";
    }

    return write_3d_smiles_zero_coordinates(m, text_output);
  }

  if (chiral_centres_in_rings.number_elements() &&
      try_removing_chiral_centres(m, chiral_centres_in_rings)) {
    molecules_converted_after_chiral_centre_removal++;

    return write_3d_smiles(out, from_corina_chemical_standardisation, text_output);
  }

  if (chiral_centres_in_chains.number_elements() &&
      try_removing_chiral_centres(m, chiral_centres_in_chains)) {
    molecules_converted_after_chiral_centre_removal++;

    return write_3d_smiles(out, from_corina_chemical_standardisation, text_output);
  }

  if (random_smiles_to_try && try_random_smiles(m)) {
    molecules_converted_after_random_smiles++;

    return write_3d_smiles(out, from_corina_chemical_standardisation, text_output);
  }

  molecules_not_converted++;

  return write_3d_smiles_zero_coordinates(m, text_output);
}

static int
rcorina34_filter(iwstring_data_source& input, std::ostream& text_output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!buffer.starts_with(smiles_tag)) {
      text_output << buffer << '\n';
      continue;
    }

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop(1);

    IWString tmp(buffer);

    if (!rcorina34_filter_process(tmp, text_output)) {
      cerr << "Cannot generate 3D coordinates, 3D TDT will be broken\n";
      molecules_not_converted++;
      return 0;
    }
  }

  return 1;
}

static int
rcorina34_filter(const char* fname, std::ostream& text_output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open stdin '" << fname << "'\n";
    return 0;
  }

  return rcorina34_filter(input, text_output);
}

static int
rcorina34_string(iwstring_data_source& input, std::ostream& output) {
  IWString buffer;

  while (input.next_record(buffer)) {
    if (!rcorina34(buffer.null_terminated_chars())) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output << out;
  }

  return 1;
}

static int
rcorina34_string(const char* fname, std::ostream& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return rcorina34_string(input, output);
}

static int
rcorina34(data_source_and_type<Molecule>& input, Molecule_Output_Object& molecule_output,
          std::ostream& output) {
  Molecule* m;
  while (NULL != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!rcorina34(*m, molecule_output, output)) {
      return 0;
    }
  }

  return 1;
}

static int
rcorina34(const char* fname, FileType input_type, Molecule_Output_Object& molecule_output,
          std::ostream& text_output) {
  assert(NULL != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return rcorina34(input, molecule_output, text_output);
}

static int
rcorina34(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lxd:F:P:S:r:e:qfsnty:R:C:huo:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  optind = 0;

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!input_chemical_standardisation.construct_from_command_line(cl, verbose > 1,
                                                                    'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('x')) {
    remove_chiral_centres = 1;

    if (verbose) {
      cerr << "Failing molecules will have chiral centres removed\n";
    }
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', max_rings_to_consider) || max_rings_to_consider < 1) {
      cerr << "The largest number of rings to consider option (-R) must be a whole +ve "
              "number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will not attempt conversion of molecules having > "
           << max_rings_to_consider << " rings\n";
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_atoms_to_consider) || max_atoms_to_consider < 1) {
      cerr << "The max atoms to consider option (-C) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will not process any molecule having more than " << max_atoms_to_consider
           << " atoms\n";
    }
  }

  if (cl.option_present('y')) {
    if (!cl.value('y', timeout) || timeout < 1) {
      cerr << "The timeout value (-y) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will timeout any calculation taking more than " << timeout << " seconds\n";
    }
    if (SIG_ERR == signal(SIGALRM, signal_handler)) {
      std::cerr << "Cannot set alarm signal handler\n";
      return 0;
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', random_smiles_to_try) || random_smiles_to_try < 0) {
      cerr << "The random smiles to try (-r) option must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will try " << random_smiles_to_try << " random smiles\n";
    }
  }

  if (cl.option_present('e')) {
    if (0 == random_smiles_to_try) {
      cerr << "Must specify random smiles to try via the -r option in order to specify "
              "the -e option\n";
      usage(3);
    }

    if (!cl.value('e', random_number_seed)) {
      cerr << "Invalid random number seed\n";
      return 2;
    }

    fix_seed_before_random_smiles = 1;

    if (verbose) {
      cerr << "Will use seed " << random_number_seed << " before randomised smiles\n";
    }
  }

  if (cl.option_present('q')) {
    pass_unique_smiles_to_corina = 0;
    pass_non_aromatic_unique_smiles_to_corina = 0;

    if (verbose) {
      cerr << "Will pass non-unique smiles to corina\n";
    }
  } else if (cl.option_present('n')) {
    pass_unique_smiles_to_corina = 0;
    pass_non_aromatic_unique_smiles_to_corina = 1;

    if (verbose) {
      cerr << "Will pass non aromatic unique smiles to corina\n";
    }
  }

  if (cl.option_present('t')) {
    tdt_output = 1;
    from_corina_chemical_standardisation.activate_from_corina_transformations();

    if (verbose) {
      cerr << "Will produce TDT output\n";
    }
  } else if (cl.option_present('u')) {
    generate_3d_smiles = 1;
    from_corina_chemical_standardisation.activate_from_corina_transformations();

    if (verbose) {
      cerr << "Will generate smiles with 3d coordinates\n";
    }
  }

  dash_d_value = "corina -i t=smiles";
  if (cl.option_present('o')) {
    dash_d_value << " -o";

    const_IWSubstring o;
    for (int i = 0; cl.value('o', o, i); ++i) {
      if (i > 0) {
        dash_d_value << ',';
      }
      dash_d_value << o;
    }
  } else {
    dash_d_value << " -o mdlcompact";
  }

  if (cl.option_present('d')) {
    int i = 0;
    const_IWSubstring d;
    while (cl.value('d', d, i++)) {
      dash_d_value << " -d " << d;
    }
  }

  if (verbose) {
    cerr << "Dash d option set to '" << dash_d_value << "'\n";
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && cl.option_present('f')) {
    function_as_tdt_filter = 1;
    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }
    from_corina_chemical_standardisation.activate_from_corina_transformations();
  } else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else if (cl.option_present('s')) {
    pass_smiles_as_string_to_corina = 1;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  Molecule_Output_Object molecule_output;

  if (cl.option_present('S')) {
    from_corina_chemical_standardisation.activate_from_corina_transformations();

    const_IWSubstring s = cl.string_value('S');

    molecule_output.add_output_type(FILE_TYPE_SDF);

    if (molecule_output.would_overwrite_input_files(cl, s)) {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 4;
    }

    if (!molecule_output.new_stem(s)) {
      cerr << "Cannot initialise molecule I/O '" << s << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Molecules written to '" << s << "'\n";
    }

    if (cl.option_present('h')) {
      remove_hydrogens_from_molecules_from_corina = 1;

      if (verbose) {
        cerr << "Will remove any explicit hydrogens added by corina\n";
      }
    }
  }

  if (cl.option_present('F')) {
    IWString f = cl.string_value('F');

    if (!f.ends_with(".smi")) {
      f << ".smi";
    }

    if (!stream_for_failing_molecules.open(f.null_terminated_chars())) {
      cerr << "Cannot open stream for failed molecules '" << f << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Smiles of failing molecules written to '" << f << "'\n";
    }
  }

  if (cl.option_present('P')) {
    IWString p = cl.string_value('P');

    if (!p.ends_with(".smi")) {
      p << ".smi";
    }

    if (!stream_for_problematic_molecules.open(p.null_terminated_chars())) {
      cerr << "Cannot open stream for problematic molecules '" << p << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Smiles of problematic molecules written to '" << p << "'\n";
    }
  }

  int rc = 0;
  if (function_as_tdt_filter) {
    if (!rcorina34_filter(cl[0], std::cout)) {
      cerr << "corina3.4 filter failed\n";
      rc = 1;
    }
  } else if (pass_smiles_as_string_to_corina) {
    for (int i = 0; i < cl.number_elements(); ++i) {
      if (!rcorina34_string(cl[i], std::cout)) {
        rc = i + 1;
        break;
      }
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!rcorina34(cl[i], input_type, molecule_output, std::cout)) {
        rc = i + 1;
        break;
      }
    }
  }

  if (molecule_output.active()) {
    molecule_output.do_close();
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    if (computation_not_attempted > 0) {
      cerr << "Computation not attempted on " << computation_not_attempted
           << " molecules due to atom (" << max_atoms_to_consider << ") and ring ("
           << max_rings_to_consider << ") count constraints\n";
    }
    cerr << molecules_converted_first_time << " molecules converted first attempt\n";
    cerr << molecules_failing_but_no_chiral_centres
         << " molecules failed, but no chiral centre(s)\n";
    if (remove_chiral_centres) {
      cerr << molecules_converted_after_chiral_centre_removal
           << " molecules converted after removing chiral centre(s)\n";
    }
    if (random_smiles_to_try > 0) {
      cerr << molecules_converted_after_random_smiles
           << " molecules converted after random smiles\n";
    }
    if (timeout > 0) {
      cerr << molecules_converted_after_timeout
           << " molecules converted after first failing a timeout\n";
    }
    if (timeout > 0 && molecules_failing_for_timeout > 0) {
      cerr << molecules_failing_for_timeout << " did not complete within timeout "
           << timeout << " seconds\n";
    }
    cerr << molecules_not_converted << " molecules not converted\n";
  } else if (computation_not_attempted) {
    cerr << computation_not_attempted
         << " conversions not attempted because of size restriction(s)\n";
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = rcorina34(argc, argv);

  return rc;
}
