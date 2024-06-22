/*
  Identify the scaffold(s) in a molecule
*/

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int add_isotope_label = 0;

static int remove_all_chiral_centres = 0;

static int remove_all_cis_trans_bonds = 0;

static Chemical_Standardisation chemical_standardisation;

static int output_rows = 0;
static int molecules_read = 0;

static int exclude_hcount_from_smarts = 1;

static IW_STL_Hash_Map_int scaffold_found;

static void
preprocess(Molecule& m) {
  m.reduce_to_largest_fragment();

  if (remove_all_chiral_centres) {
    m.remove_all_chiral_centres();
  }

  if (remove_all_cis_trans_bonds) {
    m.revert_all_directional_bonds_to_non_directional();
  }

  //  if (chemical_standardisation.active ())
  //    chemical_standardisation.process(m);

  return;
}

static void
quick_prop(Molecule& m, int& rotbond, int& hetero) {
  // calculate this after ring_membership() is called !!
  (void)m.ring_membership();
  rotbond = 0;
  hetero = 0;
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);
    const atomic_number_t zi = a->atomic_number();
    if (1 != zi && 6 != zi) {
      hetero++;
    }
    const int ncon = a->ncon();
    for (int j = 0; j < ncon; ++j) {
      const Bond* b = a->item(j);
      const atom_number_t k = b->other(i);
      if (k < i) {
        continue;
      }
      if (b->nrings()) {
        continue;
      }
      if (b->is_single_bond() && ncon > 1 && m.atomi(k)->ncon() > 1 &&
          !triple_bond_at_either_end(m, b) &&
          !part_of_otherwise_non_rotabable_entity(m, i, k)) {
        rotbond++;
      }
    }
  }
}

static int
check_csp3(Molecule& m, atom_number_t j) {
  const Atom* a = m.atomi(j);
  if (6 == a->atomic_number() && a->ncon() == a->nbonds()) {
    return a->ncon();
  }
  return 0;
}

static int
add_isotope_to_neighbors(Molecule& m, atom_number_t zatom) {
  if (!add_isotope_label) {
    return 0;
  }

  for (const Bond* b : m[zatom]) {
    atom_number_t k = b->other(zatom);
    m.set_isotope(k, add_isotope_label);
  }

  return 1;
}

static int
remove_spinach(Molecule& m, int* spinach) {
  const int matoms = m.natoms();
  if (0 == matoms) {
    return 0;
  }

  if (0 == m.nrings()) {
    return 0;
  }

  m.identify_spinach(spinach);

  for (int i = 0; i < matoms; i++) {
    if (!spinach[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);
    if (1 != a->ncon()) {
      add_isotope_to_neighbors(m, i);
      continue;
    }

    const Bond* b = a->item(0);
    if (!spinach[b->other(i)] && b->is_double_bond()) {
      spinach[i] = 0;
    }
    if (spinach[i]) {
      add_isotope_to_neighbors(m, i);
    }
  }

  m.remove_atoms(spinach);

  Set_of_Atoms to_be_removed;
  int nr = m.nrings();
  for (int i = 0; i < nr; i++) {
    const Ring* ring = m.ringi(i);
    const int na = ring->number_elements();
    if (3 == na && 12 == (check_csp3(m, ring->item(0)) * check_csp3(m, ring->item(1)) *
                          check_csp3(m, ring->item(2)))) {
      to_be_removed += (*ring);
    }
    if (4 == na && 24 == (check_csp3(m, ring->item(0)) * check_csp3(m, ring->item(1)) *
                          check_csp3(m, ring->item(2)) * check_csp3(m, ring->item(3)))) {
      to_be_removed += (*ring);
    }
  }

  if (to_be_removed.number_elements() > 0) {
    int n = to_be_removed.number_elements();
    for (int i = 0; i < n; i++) {
      add_isotope_to_neighbors(m, to_be_removed[i]);
    }
    m.remove_atoms(to_be_removed);
    if (!remove_spinach(m, spinach)) {
      return 0;
    }
  }
  return m.natoms();
}

template <typename O>
void
output_scaffold(Molecule& m, IWString& original, int number_ring_systems, O& output) {
  int matoms = m.natoms();
  if (matoms < 5 || matoms > 100) {
    return;
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (exclude_hcount_from_smarts) {
    set_include_hcount_in_smiles(0);
  }
  m.invalidate_smiles();
  const IWString usmi = m.unique_smiles();
  if (exclude_hcount_from_smarts) {
    set_include_hcount_in_smiles(1);
  }

  if (scaffold_found.contains(usmi)) {
    return;
  }

  output_rows++;
  scaffold_found[usmi] = 1;

  int rotbond, hetero;
  quick_prop(m, rotbond, hetero);

  output << usmi << ' ' << output_rows << ' ' << original << ' ' << matoms << ' '
         << rotbond << ' ' << hetero << ' ' << (m.nedges() - matoms + 1) << ' '
         << number_ring_systems;
  if (add_isotope_label) {
    output << ' ' << m.number_isotopic_atoms();
  }
  output << "\n";
}

void
UnsetImplicitHydrogenKnown(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
    m.recompute_implicit_hydrogens(i);
  }
}

template <typename O>
int
molecular_scaffold(Molecule& m, O& output) {
  IWString original;
  original << m.name() << ' ' << m.unique_smiles() << ' ' << m.natoms();
  int* spinach = new_int(m.natoms());
  std::unique_ptr<int[]> free_spinach(spinach);

  if (!remove_spinach(m, spinach)) {
    return 1;
  }

  scaffold_found.clear();

  const int matoms = m.natoms();
  int* non_linker = new int[matoms];
  std::unique_ptr<int[]> free_non_linker(non_linker);

  const int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(non_linker);

  original << ' ' << matoms;
  output_scaffold(m, original, number_ring_systems, output);

  if (1 == number_ring_systems) {
    return 1;
  }

  int* subset = new int[matoms];
  std::unique_ptr<int[]> free_subset(subset);

  // Finad all combination of ring_systems
  for (int N = (1 << number_ring_systems) - 2; N > 0; N--) {
    int sub_ringsys = number_ring_systems;
    set_vector(subset, matoms, 0);
    for (int i = 0; i < number_ring_systems; ++i) {
      if (1 == ((N >> i) & 1)) {  // ith ring_system is included in subset
        continue;
      }
      sub_ringsys--;
      for (int j = 0; j < matoms; j++) {
        if (non_linker[j] == (i + 1)) {
          subset[j] = 1;
        }
      }
    }
    Molecule msub(m);

    for (int i = 0; i < matoms; i++) {
      if (subset[i]) {
        add_isotope_to_neighbors(msub, i);
      }
    }

    msub.remove_atoms(subset);
    const int nf = msub.number_fragments();
    int with_ring = 0;
    int k_good = -1;
    for (int k = 0; k < nf; k++) {
      if (msub.rings_in_fragment(k) > 0) {
        with_ring++;
        k_good = k;
      }
      if (with_ring > 1) {
        break;
      }
    }
    if (with_ring != 1) {
      continue;
    }
    if (nf > 1) {
      msub.delete_all_fragments_except(k_good);
    }
    if (remove_spinach(msub, spinach)) {
      UnsetImplicitHydrogenKnown(msub);
      output_scaffold(msub, original, sub_ringsys, output);
    }
  }

  return 1;
}

template <typename O>
int
molecular_scaffold(data_source_and_type<Molecule>& input, O& output)

{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    if (verbose > 1) {
      cerr << "Processing '" << m->name() << "'\n";
    }

    if (!molecular_scaffold(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

template <typename O>
int
molecular_scaffold(const char* fname, FileType input_type, O& output) {
  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Processing '" << fname << "'\n";
  }

  return molecular_scaffold(input, output);
}

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
  cerr << "  -S <string>    create output files with name <string>\n";
  cerr << "  -I <iso>      add isotope label to connection point\n";
  cerr << "  -i <type>      specify input type\n";
  display_standard_aromaticity_options (cerr);
  display_standard_smiles_options (cerr);
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -c            remove all chirality from input molecules\n";
  cerr << "  -t            remove cis-trans bonds from input\n";
  cerr << "  -x            include implicit Hydrogen information in unique smiles\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
molecular_scaffold(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:I:A:K:g:S:ctx");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_smiles_options(cl, verbose)) {
    usage(5);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(6);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (!cl.option_present('i')) {
    if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot discern all file types, use the -i option\n";
      usage(3);
    }
  } else if (!process_input_type(cl, input_type)) {
    usage(3);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('c')) {
    remove_all_chiral_centres = 1;
    if (verbose) {
      cerr << "Chirality will be removed\n";
    }
  }

  if (cl.option_present('t')) {
    remove_all_cis_trans_bonds = 1;
    if (verbose) {
      cerr << "Cis trans bonds will be removed\n";
    }
  }

  if (cl.option_present('x')) {
    exclude_hcount_from_smarts = 0;

    if (verbose) {
      cerr << "Implicit Hydrogen information included in unique smiles\n";
    }
  }

  if (cl.option_present('I')) {
    if (!cl.value('I', add_isotope_label) || add_isotope_label < 0) {
      cerr << "The isotope for connection points (-I) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will add isotope " << add_isotope_label << " to join points\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(87);
  }

  IWString_and_File_Descriptor* output;

  if (cl.option_present('S')) {
    if (1 != cl.option_count('S')) {
      cerr << "Sorry, only one -S option possible\n";
      usage(3);
    }
    IWString output_file_stem;
    cl.value('S', output_file_stem);

    if ("-" == output_file_stem) {
      output = new IWString_and_File_Descriptor(1);
    } else {
      output = new IWString_and_File_Descriptor();
      if (!output->open(output_file_stem)) {
        cerr << "Cannot open '" << output_file_stem << "'\n";
        return 1;
      }
      if (verbose) {
        cerr << "Output to'" << output_file_stem << "'\n";
      }
    }
  } else {
    output = new IWString_and_File_Descriptor(1);
  }

  *output << "FragmentSMILES ROWID ParentId ParentSMILES ParentNatoms ScaffoldNatoms "
             "FragmentNatoms FragmentRotBond FragmentHeteroAtom FragmentNrings "
             "FragmentNringSystem";
  if (add_isotope_label) {
    *output << " FragmentSubstitutions";
  }
  *output << "\n";

  int rc = 0;
  for (int i = 0; i < cl.number_elements() && 0 == rc; i++) {
    if (!molecular_scaffold(cl[i], input_type, *output)) {
      rc = i + 1;
      break;
    }
  }

  output->close();

  delete output;

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = molecular_scaffold(argc, argv);

  return rc;
}
