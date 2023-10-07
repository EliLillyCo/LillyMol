/*
  Computes topological torsions
*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <memory>
#include <sstream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "topotorsion_support.h"
#include "torsion_hash.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;
static int molecules_read = 0;

static int function_as_filter = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static Atom_Typing_Specification atom_typing_specification;

/*
  If we are writing our results as a fingerprint, we need a tag
*/

static IWString tag;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int flush_after_every_molecule = 0;

static int
write_no_bits(IWString_and_File_Descriptor& output)
{
  output << tag << ">\n";

  return output.good();
}

/*
  The maximum invariant is 250, so the invariants are stored as base 250 numbers
*/

static void
topotorsion_fingerprints(Molecule& m, const int* ncon, const int* invariant,
                         atom_number_t a0, atom_number_t a1, atom_number_t a2,
                         atom_number_t a3, Sparse_Fingerprint_Creator& fp)
{
  form_canonical_order(invariant, a0, a1, a2, a3);

  unsigned int b = 15625000 * invariant[a0] + 62500 * invariant[a1] +
                   250 * invariant[a2] + invariant[a3];

  fp.hit_bit(b);

  return;
}

static void
topotorsion_fingerprints(Molecule& m, const int* ncon, const int* invariant,
                         atom_number_t a0, atom_number_t a1, atom_number_t a2,
                         Sparse_Fingerprint_Creator& fp)
{
  int nc2 = ncon[a2];

  const Atom* a = m.atomi(a2);

  assert(nc2 == a->ncon());

  for (int i = 0; i < nc2; i++) {
    atom_number_t a3 = a->other(a2, i);
    if (a3 == a1) {
      continue;
    }

    topotorsion_fingerprints(m, ncon, invariant, a0, a1, a2, a3, fp);
  }

  return;
}

static int
topotorsion_fingerprints(Molecule& m, const int* ncon, const int* invariant,
                         const Bond* b, Sparse_Fingerprint_Creator& fp)
{
  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();
  int nc1 = ncon[a1];

  if (1 == nc1 || 1 == ncon[a2]) {  // A1 and A2 are in the middle of a torsion
    return 0;
  }

  const Atom* a = m.atomi(a1);

  for (int i = 0; i < nc1; i++) {
    atom_number_t o = a->other(a1, i);
    if (o == a2) {
      continue;
    }

    topotorsion_fingerprints(m, ncon, invariant, o, a1, a2, fp);
  }

  return 1;
}

static int
topotorsion_fingerprints(Molecule& m, const int* ncon, const int* invariant,
                         IWString_and_File_Descriptor& output)
{
  Sparse_Fingerprint_Creator fp;

  int nb = m.nedges();

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    topotorsion_fingerprints(m, ncon, invariant, b, fp);
  }

  if (verbose > 1) {
    cerr << "Found " << fp.nbits() << " topological torsions\n";
  }

  if (verbose > 1) {
    fp.debug_print(cerr);
  }

  int rc;
  if (0 == fp.nbits()) {
    rc = write_no_bits(output);
  } else {
    rc = fp.write_fingerprint(tag, output);
  }

  return rc;
}

static void
preprocess(Molecule& m)
{
  m.remove_all(1);  // no hydrogens

  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
topotorsion_fingerprints(Molecule& m, IWString_and_File_Descriptor& output)
{
  preprocess(m);

  int matoms = m.natoms();

  if (matoms < 4) {
    cerr << "Not enough atoms, " << matoms << " in molecule '" << m.name()
         << "'\n";
    write_no_bits(output);
    return output.good();
  }

  int* ncon = new int[matoms];
  std::unique_ptr<int[]> free_ncon(ncon);
  m.ncon(ncon);

  int* invariant = new int[matoms];
  std::unique_ptr<int[]> free_invariant(invariant);

  atom_typing_specification.assign_atom_types(m, invariant, ncon);

  return topotorsion_fingerprints(m, ncon, invariant, output);
}

static int
topotorsion_fingerprints(const const_IWSubstring& buffer,
                         IWString_and_File_Descriptor& output)
{
  assert(buffer.ends_with('>'));
  const_IWSubstring smiles(buffer);

  smiles.remove_up_to_first('<');
  smiles.chop();

  Molecule m;

  if (!m.build_from_smiles(smiles)) {
    cerr << "Cannot parse smiles '" << smiles << "'\n";
    return 0;
  }

  molecules_read++;

  return topotorsion_fingerprints(m, output);
}

static int
topotorsion_fingerprints(iwstring_data_source& input,
                         IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!topotorsion_fingerprints(buffer, output)) {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

static int
topotorsion_fingerprints_filter(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return topotorsion_fingerprints(input, output);
}

static int
topotorsion_fingerprints(data_source_and_type<Molecule>& input,
                         IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    output << smiles_tag << m->smiles() << ">\n";
    topotorsion_fingerprints(*m, output);
    output << identifier_tag << m->name() << ">\n";
    output << "|\n";

    if (flush_after_every_molecule) {
      output.flush();
    } else {
      output.write_if_buffer_holds_more_than(32768);
    }
  }

  return output.good();
}

static int
topotorsion_fingerprints(const char* fname, FileType input_type,
                         IWString_and_File_Descriptor& output)
{
  if (function_as_filter) {
    return topotorsion_fingerprints_filter(fname, output);
  }

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  return topotorsion_fingerprints(input, output);
}

static void
DisplayDashXOptions(std::ostream& output)
{
  output << " -X flush      flush output after each molecule\n";

  ::exit(0);
}

static void
usage(int rc = 0)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Usage : " << prog_name << " options file1 file2 file3 ....\n";
  cerr << "  -J <tag>       tag for fingerprints\n";
  cerr << "  -P <atype>     atom typing to use, enter '-P help' for info\n";
  cerr << "  -E <symbol>    create element with symbol\n";
  cerr << "  -f             work as a filter\n";
  cerr << "  -l             strip to largest fragment\n";
  cerr << "  -i <type>      specify input type\n";
  display_standard_aromaticity_options(cerr);
   display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

int
topotorsion(int argc, char** argv)
{
  Command_Line cl(argc, argv, "A:E:vJ:i:fg:lP:X:");

  if (cl.unrecognised_options_encountered()) {
    usage();
  }

  verbose = cl.option_count('v');

  process_elements(cl);

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(7);
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will strip to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g option)\n";
      usage(5);
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Cannot determine atom typying to use '" << p << "'\n";
      return 3;
    }

    if (0 == tag.length()) {
      tag = "NCTO";
      atom_typing_specification.append_to_tag(tag);
      tag << '<';
    }
  } else {
    atom_typing_specification.set_atom_type(IWATTYPE_TT);
    tag = "NCTO<";
  }

  if (cl.option_present('J')) {
    tag = cl.string_value('J');
    if (verbose) {
      cerr << "Torsions written as non-colliding sparse fingerprints, tag '" << tag
           << "'\n";
    }
  } else if (0 == tag.length()) {
    tag = "NCTO<";
  }

  if (!tag.ends_with('<')) {
    tag << '<';
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

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('f')) {
    function_as_filter = 1;
  } else if (!cl.option_present('i')) {
    if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot determine input type(s)\n";
      return 7;
    }
  } else if (!process_input_type(cl, input_type)) {
    cerr << "Cannot determine input type\n";
    usage(6);
  }

  if (cl.empty()) {
    cerr << prog_name << " insufficient arguments\n";
    usage();
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)  // each argument is a file
  {
    const char* fname = cl[i];

    if (verbose) {
      cerr << prog_name << " processing '" << fname << "'\n";
    }

    if (!topotorsion_fingerprints(fname, input_type, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << molecules_read << " molecules read" << endl;
    cerr << endl;
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = topotorsion(argc, argv);

  return rc;
}
