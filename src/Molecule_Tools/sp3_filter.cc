/*
  Filter molecules according to SP3 groups
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_written = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int min_carbon_sp3 = 0;
static int min_non_carbon_sp3 = 0;

static int rejected_for_too_few_carbon_sp3 = 0;
static int rejected_for_too_few_non_carbon_sp3 = 0;

static IWString_and_File_Descriptor stream_for_rejected;

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
  cerr << "  -c <number>   min number of Carbon     sp3 atoms\n";
  cerr << "  -x <number>   min number of non-Carbon sp3 atoms\n";
  cerr << "  -U <fname>    stream for rejected molecules\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess(Molecule& m)
{
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
do_write(Molecule& m, IWString_and_File_Descriptor& output)
{
  output << m.smiles() << ' ' << m.name() << "\n";

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

static int
identify_cf3_tbutyl(Molecule& m, int* sp3)
{
  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; i--) {
    if (0 == sp3[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    if (4 != a->ncon()) {
      continue;
    }

    int singly_connected = 0;

    for (int j = 0; j < 4; j++) {
      atom_number_t k = a->other(i, j);

      if (1 != m.ncon(k)) {
        continue;
      }

      atomic_number_t zk = m.atomic_number(k);

      if (6 == zk || 9 == zk || 17 == zk) {
        singly_connected++;
      }
    }

    if (3 == singly_connected) {
      sp3[i] = 0;
      rc++;
    }
  }

  return rc;
}

/*
 */

static int
sp3_filter(Molecule& m, IWString_and_File_Descriptor& output)
{
  int matoms = m.natoms();

  int* sp3 = new_int(matoms);
  std::unique_ptr<int[]> free_sp3(sp3);

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      sp3[i] = 0;
    } else if (1 == m.ncon(i)) {
      sp3[i] = 0;
    } else if (m.nbonds(i) == m.ncon(i)) {  // fully saturated
      sp3[i] = 1;
    } else {
      sp3[i] = 0;
    }
  }

  identify_cf3_tbutyl(m, sp3);

  int carbon_sp3 = 0;
  int non_carbon_sp3 = 0;

  for (int i = 0; i < matoms; i++) {
    if (0 == sp3[i]) {
      continue;
    }

    if (6 == m.atomic_number(i)) {
      carbon_sp3++;
    } else {
      non_carbon_sp3++;
    }
  }

  if (carbon_sp3 < min_carbon_sp3) {
    rejected_for_too_few_carbon_sp3++;
    if (stream_for_rejected.is_open()) {
      return do_write(m, stream_for_rejected);
    }
    return 1;
  }

  if (non_carbon_sp3 < min_non_carbon_sp3) {
    rejected_for_too_few_non_carbon_sp3++;
    if (stream_for_rejected.is_open()) {
      return do_write(m, stream_for_rejected);
    }
    return 1;
  }

  molecules_written++;

  return do_write(m, output);
}

static int
sp3_filter(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!sp3_filter(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
sp3_filter(const char* fname, FileType input_type, IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return sp3_filter(input, output);
}

static int
sp3_filter(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lc:x:U:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
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

  if (!cl.option_present('c') && !cl.option_present('x')) {
    cerr << "Must specify acceptance criteria via either -c and/or -x options\n";
    usage(2);
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', min_carbon_sp3) || min_carbon_sp3 < 1) {
      cerr << "The minimum number of Carbon sp3 atoms (-c) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Molecules will be passed if they have at least " << min_carbon_sp3
           << " sp3 Carbon atoms\n";
    }
  }

  if (cl.option_present('x')) {
    if (!cl.value('x', min_non_carbon_sp3) || min_non_carbon_sp3 < 1) {
      cerr << "The minimum number of non-Carbon sp3 atoms (-x) must be a whole +ve "
              "number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Molecules will be passed if they have at least " << min_non_carbon_sp3
           << " non Carbon sp3 atoms\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('U')) {
    IWString u = cl.string_value('U');
    if (!u.ends_with(".smi")) {
      u << ".smi";
    }

    if (!stream_for_rejected.open(u.null_terminated_chars())) {
      cerr << "Cannot open stream for rejected molecules '" << u << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Rejected molecules written to '" << u << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!sp3_filter(cl[i], input_type, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules, wrote " << molecules_written
         << endl;
    cerr << rejected_for_too_few_carbon_sp3 << " rejected for too few sp3 Carbon atoms\n";
    cerr << rejected_for_too_few_non_carbon_sp3
         << " rejected for too few non Carbon sp3 atoms\n";
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = sp3_filter(argc, argv);

  return rc;
}
