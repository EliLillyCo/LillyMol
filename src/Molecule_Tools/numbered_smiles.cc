/*
  Frequently I need to know the atom numbers of the atoms in a molecule.
  Produce isotopic smiles
*/

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;

static int verbose = 0;
const char* prog_name = nullptr;

static IWString output_file_name;

static Element_Transformations element_transformations;

/*
  The -n switch governs the numbering assigned.
*/

enum class NumberingType {
  kAtomNumber,
  kSymmetryClass,
  kCanonicalRank,
  kRingSystem,
  kAtomType,
  kAtomicNumber
};

// #define NUMBER_BY_ATOM_NUMBER 1
// #define NUMBER_BY_SYMMETRY_CLASS 2
// #define NUMBER_BY_CANONICAL_RANK 3

static NumberingType numbering_type = NumberingType::kAtomNumber;

static int atom_number_offset = 0;

static int number_step = 1;

static int overwrite_existing_isotopic_information = 0;

static int sort_molecule_by_numbering = 0;

static int change_to_graph_form = 0;

static int make_implicit_hydrogens_explicit = 0;

static int compute_bond_symmetry_small = 0;
static int compute_bond_symmetry_large = 0;

static int ring_systems_span_spiro_fusions = 0;

static resizable_array_p<Substructure_Query> queries;

static int apply_isotopic_labels = 1;

static Atom_Typing_Specification atom_typing;

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
  cerr << "  -n atom        number by atom number (the default)\n";
  cerr << "  -n symm        number by symmetry class\n";
  cerr << "  -n canon       number by canonical order\n";
  cerr << "  -n ring        number by ring system\n";
  cerr << "  -n z           number by atomic number\n";
  cerr << "  -n <number>    number by atom number starting with <number>\n";
  cerr << "  -e <step>      number step (default 1)\n";
  cerr << "  -I             overwrite any existing isotopic information\n";
  cerr << "  -P <atype>     number by atom typing (no -n option needed)\n";
  cerr << "  -r             sort molecule by numbering\n";
  cerr << "  -p             change to graph form (all bonds become single etc...)\n";
  cerr << "  -H             make all implicit Hydrogens explicit\n";
  cerr << "  -b ...         write bond symmetry information to stderr\n";
  cerr << "  -q <query>     only label atoms matched by the query/queries\n";
  cerr << "  -s <smarts>    only label atoms matched by the smarts\n";
  cerr << "  -m             apply atom map numbers instead of isotopic labels\n";
  cerr << "  -Y ...         more options, enter '-Y help' for info\n";
  cerr << "  -S <stem>      specify output file name stem\n";
  cerr << "  -K ...         standard smiles control options\n";
  cerr << "  -t ...         element transformation options\n";
  cerr << "  -i <type>      specify input type\n";
  display_standard_aromaticity_options(cerr);
  ;
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int
leave_isotopes_only_on_matched_atoms(Molecule& m)
{
  int matoms = m.natoms();

  int* keep_isotope = new_int(matoms);
  std::unique_ptr<int[]> free_keep_isotope(keep_isotope);

  Molecule_to_Match target(&m);

  int nq = queries.number_elements();

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);

      e->set_vector(keep_isotope, 1);
    }
  }

  for (int i = 0; i < matoms; i++) {
    if (keep_isotope[i]) {
      continue;
    }

    if (apply_isotopic_labels) {
      m.set_isotope(i, 0);
    } else {
      m.set_atom_map_number(i, 0);
    }
  }

  return 1;
}

static int
write_bond_symmetries(Molecule& m, const int* bs, std::ostream& output)
{
  output << m.name() << '\n';

  int nb = m.nedges();

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    output << b->a1() << " (" << m.symmetry_class(b->a1()) << ") " << b->a2() << " ("
           << m.symmetry_class(b->a2()) << ") " << bs[i] << '\n';
  }

  return output.good();
}

static int
do_compute_bond_symmetry_small(Molecule& m, std::ostream& output)
{
  int* bs = new int[m.nedges()];
  std::unique_ptr<int[]> free_bs(bs);

  m.bond_symmetry_class_small_memory(bs);

  return write_bond_symmetries(m, bs, output);
}

static int
do_compute_bond_symmetry_large(Molecule& m, std::ostream& output)
{
  int* bs = new int[m.nedges()];
  std::unique_ptr<int[]> free_bs(bs);

  m.bond_symmetry_class_large_memory(bs);

  return write_bond_symmetries(m, bs, output);
}

static int
adjust_symmetry_numbers(int* isotope, int n, int* new_number)
{
  int next_number_to_assign = 0;

  for (int i = 0; i < n; i++) {
    assert(isotope[i] > 0 && isotope[i] <= n);

    if (new_number[i] >= 0) {  // already processed
      continue;
    }

    new_number[i] = next_number_to_assign;
    next_number_to_assign++;

    for (int j = i + 1; j < n; j++) {
      if (isotope[i] == isotope[j])  // in same symmetry class
      {
        new_number[j] = next_number_to_assign;
        next_number_to_assign++;
      }
    }
  }

  for (int i = 0; i < n; i++) {
    isotope[i] = new_number[i];
  }

  return 1;
}

/*
  We need to make sure we have a set of unique numbers to give to renumber_atoms
*/

static int
adjust_numbers_and_do_sort(Molecule& m, int* isotope) {

  if (numbering_type == NumberingType::kAtomNumber) {
    return m.renumber_atoms(isotope);
  }

  const int matoms = m.natoms();

  // Canonical ranks need to be decremented by 1

  if (numbering_type == NumberingType::kCanonicalRank) {
    for (int i = 0; i < matoms; i++) {
      isotope[i]--;
      assert(isotope[i] >= 0 && isotope[i] < matoms);
    }

    return m.renumber_atoms(isotope);
  }

  // Symmetry classes are more difficult, because there may be multiple atoms
  // with the same symmetry number, but renumber_atoms needs unique numbers

  int* tmp = new_int(matoms, -1);
  std::unique_ptr<int[]> free_tmp(tmp);

  adjust_symmetry_numbers(isotope, matoms, tmp);

  return m.renumber_atoms(isotope);
}

static int
adjust_numbers_and_do_sort(Molecule& m, isotope_t* isotope)
{
  const int matoms = m.natoms();
  std::unique_ptr<int[]> as_int = std::make_unique<int[]>(matoms);
  std::copy_n(isotope, matoms, as_int.get());

  int rc = adjust_numbers_and_do_sort(m, as_int.get());

  std::copy_n(as_int.get(), matoms, isotope);
  return rc;
}

static int
do_apply_isotopic_labels(Molecule& m, const isotope_t* isotope)
{
  if (overwrite_existing_isotopic_information) {
    return m.set_isotopes(isotope);
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    if (m.isotope(i)) {
      continue;
    }

    if (apply_isotopic_labels) {
      m.set_isotope(i, isotope[i]);
    } else {
      m.set_atom_map_number(i, isotope[i]);
    }
  }

  return 1;
}

static int
fetch_atom_numbers(Molecule& m, isotope_t* isotope)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    isotope[i] = atom_number_offset + i;
  }

  return 1;
}

static int
fetch_symmetry_class(Molecule& m, isotope_t* isotope)
{
  std::copy_n(m.symmetry_classes(), m.natoms(), isotope);

  return 1;
}

// The canonical rank is an integer, so a temporary is needed.
static int
fetch_canonical_rank(Molecule& m, isotope_t* isotope)
{
  const int matoms = m.natoms();
  std::unique_ptr<int[]> as_int = std::make_unique<int[]>(matoms);

  m.canonical_ranks(as_int.get());

  std::copy_n(as_int.get(), matoms, isotope);

  return 1;
}

static int
AddDoublyBondedExocyclic(Molecule& m,
                         isotope_t* isotope) {

  // Add doubly bonded exocyclic bonds
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (isotope[i] == 0) {
      continue;
    }
    const Atom* a = m.atomi(i);
    if (a->ncon() == 2 || a->nbonds() == a->ncon()) {
      continue;
    }
    for (const Bond* b : *a) {
      if (!b->is_double_bond()) {
        continue;
      }
      const atom_number_t j = b->other(i);
      if (isotope[j] == isotope[i] || m.atomic_number(j) == 6) {
        continue;
      }
      // We have an exocyclic double bond.
      // If we want to restrict to just =O or =S, uncomment this.
      // if (m.ncon(j) != 1) {
      //   continue;
      // }
      isotope[j] = isotope[i];
    }
  }

  return 1;
}

static int
fetch_ring_system_number(Molecule& m, isotope_t* isotope)
{
  const int nrings = m.nrings();
  if (nrings == 0) {
    return 1;
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> label(new_int(matoms));

  if (ring_systems_span_spiro_fusions) {
    m.label_atoms_by_ring_system_including_spiro_fused(label.get());
  } else {
    m.label_atoms_by_ring_system(label.get());
  }
  std::copy_n(label.get(), matoms, isotope);

  return AddDoublyBondedExocyclic(m, isotope);
}

int
fetch_atomic_number(const Molecule& m, isotope_t* isotope)
{
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    isotope[i] = m.atomic_number(i);
  }

  return matoms;
}

// Needs a temporary array to hold the right int type for
// atom typing.
int
fetch_atom_types(Molecule& m, isotope_t* isotope) {
  const int matoms = m.natoms();
  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);

  atom_typing.assign_atom_types(m, atype.get());

  std::copy_n(atype.get(), matoms, isotope);

  return 1;
}

static int
numbered_smiles(Molecule& m,
                isotope_t* isotope,
                Molecule_Output_Object& output)
{
  switch (numbering_type) {
    case NumberingType::kAtomNumber:
      fetch_atom_numbers(m, isotope);
      break;
    case NumberingType::kSymmetryClass:
      fetch_symmetry_class(m, isotope);
      break;
    case NumberingType::kCanonicalRank:
      fetch_canonical_rank(m, isotope);
      break;
    case NumberingType::kRingSystem:
      fetch_ring_system_number(m, isotope);
      break;
    case NumberingType::kAtomType:
      fetch_atom_types(m, isotope);
      break;
    case NumberingType::kAtomicNumber:
      fetch_atomic_number(m, isotope);
      break;
  }

  do_apply_isotopic_labels(m, isotope);

  if (sort_molecule_by_numbering) {
    adjust_numbers_and_do_sort(m, isotope);
  }

  if (number_step > 1) {
    const int matoms = m.natoms();
    for (int i = 1; i < matoms; i++) {
      isotope_t iso = m.isotope(i);
      if (iso > 0) {
        if (apply_isotopic_labels) {
          m.set_isotope(i, iso * number_step);
        } else {
          m.set_atom_map_number(i, iso * number_step);
        }
      }
    }
  }

  if (queries.number_elements() > 0) {
    leave_isotopes_only_on_matched_atoms(m);
  }

  if (compute_bond_symmetry_small) {
    do_compute_bond_symmetry_small(m, std::cout);
  }
  if (compute_bond_symmetry_large) {
    do_compute_bond_symmetry_large(m, std::cout);
  }

  return output.write(m);
}

static int
numbered_smiles(Molecule& m, Molecule_Output_Object& output)
{
  assert(m.ok());

  if (change_to_graph_form) {
    m.change_to_graph_form();
  }

  int matoms = m.natoms();

  if (0 == matoms) {
    cerr << "numbered_smiles: nostructure encountered '" << m.name() << "'\n";
    return 1;
  }

  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
    matoms = m.natoms();
  }

  std::unique_ptr<isotope_t[]> isotope = std::make_unique<isotope_t[]>(matoms);
  std::fill_n(isotope.get(), matoms, 0);

  return numbered_smiles(m, isotope.get(), output);
}

static int
numbered_smiles(data_source_and_type<Molecule>& input, Molecule_Output_Object& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (verbose) {
      cerr << "Processing '" << input.molecules_read() << " " << m->name() << "'\n";
    }

    if (!numbered_smiles(*m, output)) {
      return 0;
    }
  }

  if (verbose) {
    cerr << "Processed " << input.molecules_read() << " molecules\n";
  }

  return output.good();
}

static int
numbered_smiles(const char* fname, FileType input_type, Molecule_Output_Object& output)
{
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  static int first_call = 1;

  if (0 == output_file_name.length()) {
    IWString ofname = fname;
    int i = ofname.index('.');
    if (i) {
      ofname.iwtruncate(i);
    }

    ofname += "_num.smi";

    if (!output.new_stem(ofname)) {
      cerr << "Cannot set new stem '" << ofname << "'\n";
      return 99;
    }
  } else if (first_call) {
    if (!output.new_stem(output_file_name)) {
      cerr << "Cannot set stem for '" << output_file_name << "'\n";
      return 13;
    }
    first_call = 0;
  }

  return numbered_smiles(input, output);
}

static void
DisplayDashYOptions(std::ostream& output) {
  output << " -Y spiro          ring systems span spiro fusions\n";

  ::exit(0);
}

static int
numbered_smiles(int argc, char** argv)
{
  Command_Line cl(argc, argv, "Ivq:s:q:S:i:o:n:E:A:rpHb:e:mK:P:Y:t:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(6);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl, verbose)) {
    usage(1);
  }

  if (!process_standard_smiles_options(cl, verbose, 'K')) {
    cerr << "Cannot process smiles options\n";
  }

  if (cl.option_present('t')) {
    if (!element_transformations.construct_from_command_line(cl, verbose)) {
      cerr << "Cannot process element transformations\n";
      usage(2);
    }
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process standard aromaticity options\n";
    usage(3);
  }

  if (cl.option_present('I')) {
    overwrite_existing_isotopic_information = 1;
    if (verbose) {
      cerr << "Existing isotopic information will be overwritten\n";
    }
  }

  if (cl.option_present('P')) {
    const IWString p = cl.string_value('P');
    if (! atom_typing.build(p)) {
      cerr << "Cannot build atom typing '" << p << "'\n";
      return 1;
    }
    numbering_type = NumberingType::kAtomType;
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "spiro") {
        ring_systems_span_spiro_fusions = 1;
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('n')) {
    IWString tmp;
    int i = 0;
    while (cl.value('n', tmp, i++)) {
      if (tmp.starts_with("atom")) {
        numbering_type = NumberingType::kAtomNumber;
        if (verbose) {
          cerr << "Atoms will be numbered according to their atom number\n";
        }
      } else if (tmp.starts_with("symm")) {
        numbering_type = NumberingType::kSymmetryClass;
        if (verbose) {
          cerr << "Atoms will be numbered by symmetry class\n";
        }
      } else if (tmp.starts_with("canon")) {
        numbering_type = NumberingType::kCanonicalRank;
        if (verbose) {
          cerr << "Atoms will be numbered by canonical number\n";
        }
      } else if (tmp.starts_with("ring")) {
        numbering_type = NumberingType::kRingSystem;
        if (verbose) {
          cerr << "Atoms numbered by ring system\n";
        }
      } else if (tmp.starts_with("z")) {
        numbering_type = NumberingType::kAtomicNumber;
        if (verbose) {
          cerr << "Atoms numbered by atomic number\n";
        }
      } else if (tmp.is_int(atom_number_offset)) {
        if (verbose) {
          cerr << "Atoms numbered by atom number, offset = " << atom_number_offset
               << '\n';
        }
      } else {
        cerr << "Unrecognised numbering type '" << tmp << "'\n";
        usage(9);
      }
    }
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', number_step) || number_step < 1) {
      cerr << "The number step value (-e) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Gap between numbers set to " << number_step << '\n';
    }
  }

  if (cl.option_present('b')) {
    const_IWSubstring b = cl.string_value('b');

    if ("small" == b) {
      compute_bond_symmetry_small = 1;
      if (verbose) {
        cerr << "Will compute bond symmetry with the small memory model\n";
      }
    } else if ("large" == b) {
      compute_bond_symmetry_large = 1;
      if (verbose) {
        cerr << "Will compute bond symmetry with the large memory model\n";
      }
    } else {
      cerr << "Unrecognised -b qualifier '" << b << "'\n";
      usage(5);
    }
  }

  if (cl.option_present('m')) {
    apply_isotopic_labels = 0;
    if (verbose) {
      cerr << "Will apply atom map numbers\n";
    }
  }

  // By default we echo any chiral info present in the input

  set_include_chiral_info_in_smiles(1);

  if (cl.option_present('H')) {
    make_implicit_hydrogens_explicit = 1;

    if (verbose) {
      cerr << "Implicit hydrogens will be made explicit\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot determine input types\n";
    return 32;
  } else  // good
  {
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output types (-o)\n";
      return 77;
    }
  } else {
    output.add_output_type(FILE_TYPE_SMI);  // our default
  }

  if (cl.option_present('S')) {
    cl.value('S', output_file_name);
    if (0 == output_file_name.length()) {
      cerr << "The -S option must be followed by a stem\n";
      usage(10);
    }

    for (int i = 0; i < cl.number_elements(); i++) {
      if (output.would_use_name(output_file_name.null_terminated_chars(), cl[i])) {
        cerr << "Input and output must be distinct, '" << cl[i] << "'\n";
        return 88;
      }
    }

    if (verbose) {
      cerr << "Stem for output is '" << output_file_name << "'\n";
    }
  }

  if (cl.option_present('r')) {
    sort_molecule_by_numbering = 1;

    if (verbose) {
      cerr << "Molecule will be sorted by numbering scheme\n";
    }
  }

  if (cl.option_present('p')) {
    change_to_graph_form = 1;

    if (verbose) {
      cerr << "Will change to graph form\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(9);
  }

  for (const char* fname : cl) {
    if (!numbered_smiles(fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = numbered_smiles(argc, argv);

  return rc;
}
