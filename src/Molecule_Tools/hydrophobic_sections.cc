/*
  Custom descriptors for Dave Cummins work with Phospholipidosis
*/

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/target.h"

using std::cerr;

static int verbose = 0;

static float float_zero = static_cast<float>(0.0);

static int molecules_read = 0;

static int strip_to_largest_fragment = 1;

static resizable_array_p<Substructure_Hit_Statistics> hydrophobic_queries;
static resizable_array_p<Substructure_Hit_Statistics> hydrophillic_queries;

static Accumulator_Int<int> hydrophobic_atom_count;
static Accumulator_Int<int> hydrophillic_atom_count;

static int two_connected_sulphur_is_greasy = 1;

static int cf3_is_hydrophillic = 0;

static int apply_simple_hydrophobic_rules = 0;
static int apply_simple_hydrophillic_rules = 0;

static int non_hydrophobic_atoms_are_hydrophillic = 0;

static charge_t charge_needed_for_hydrophillic = static_cast<charge_t>(0.0);
static charge_t charge_needed_for_hydrophobic = static_cast<charge_t>(0.0);

static int partial_charge_specification_present = 0;  // if either charge is specified

/*
  We can invert things and get hydrophillic info
*/

static int create_hydrophillic_descriptors = 0;

/*
  It can be handy to write out the molecule with the hydrophobic atoms identified
*/

static int hydrophobic_isotope = 0;
static int hydrophillic_isotope = 0;

static Molecule_Output_Object stream_for_labelled_molecules;

/*
  In some contexts we can decide that hydrophobic sections less than a given size
  aren't hydrophobic at all
*/

static int min_hydrophobic_section_size = 1;
static int min_hydrophillic_section_size = 1;

static int min_min_separation = 2;

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
  cerr << "  -G <smarts>    SMARTS query to designate hydrophobic atoms\n";
  cerr << "  -G q:<file>    hydrophobic query file\n";
  cerr << "  -G Q:<file>    hydrophobic queries in <file>\n";
  cerr << "  -G minsize=nn  minimum size for a hydrophobic section\n";
  cerr << "  -G charge=q    max partial charge for a hydrophobic atom\n";
  cerr << "  -G iso=nn      write hydrophillic atoms with isotope <nn>\n";
  cerr << "  -G simple      use simple builtin rules for hydrophobicity\n";
  cerr << "  -L <smarts>    SMARTS query to designate hydrophillic atoms\n";
  cerr << "  -L q:<file>    hydrophillic query file\n";
  cerr << "  -L Q:<file>    hydrophillic queries in <file>\n";
  cerr << "  -L minsize=nn  minimum size for a hydrophillic section\n";
  cerr << "  -L charge=q    min partial charge for a hydrophillic atom\n";
  cerr << "  -L iso=nn      write hydrophillic atoms with isotope <nn>\n";
  cerr << "  -L simple      use simple builtin rules for hydrophillic regions\n";
  cerr << "  -x             hydrophillic_atoms are all non-hydrophobic atoms\n";
  cerr << "  -e <num>       minimum separation between sections (default " << min_min_separation << ")\n";
  cerr << "  -S <fname>     file name for labelled molecules\n";
  cerr << "  -n             suppress writing descriptors\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
do_write_descriptors(const resizable_array<float>& descriptors, IWString& buffer) {
  int n = descriptors.number_elements();

  for (int i = 0; i < n; i++) {
    buffer += ' ';
    buffer.append_number(descriptors[i], 3);
  }

  return 1;
}

static void
compute_separation(Molecule& m, const Set_of_Atoms& s1, const Set_of_Atoms& s2,
                   int& min_separation, int& max_separation) {
  min_separation = m.natoms();
  max_separation = 0;

  int n1 = s1.number_elements();
  int n2 = s2.number_elements();

  for (int i = 0; i < n1; i++) {
    atom_number_t a1 = s1[i];

    for (int j = 0; j < n2; j++) {
      atom_number_t a2 = s2[j];

      int d = m.bonds_between(a1, a2);

      if (d < min_separation) {
        min_separation = d;
      }

      if (d > max_separation) {
        max_separation = d;
      }
    }
  }

  return;
}

static int
get_subset(int n, const int* zarray, int flag, Set_of_Atoms& s) {
  s.resize_keep_storage(0);

  for (int i = 0; i < n; i++) {
    if (flag == zarray[i]) {
      s.add(i);
    }
  }

  return s.number_elements();
}

static void
no_inter_descriptors(resizable_array<float>& inter_descriptors) {
  for (int i = 0; i < 4; i++) {
    inter_descriptors.add(float_zero);
  }

  return;
}

static void
do_create_inter_descriptors(Molecule& m, const int* hyphob, const int* hyphil,
                            resizable_array<float>& inter_descriptors) {
  int matoms = m.natoms();

  int hydrophobic_sections = iwmax_of_array(hyphob, matoms);
  int hydrophillic_sections = iwmax_of_array(hyphil, matoms);

  if (0 == hydrophobic_sections || 0 == hydrophillic_sections) {
    no_inter_descriptors(inter_descriptors);
    return;
  }

  Accumulator_Int<int> separation;

  Set_of_Atoms si, sj;  // scope here for efficiency

  for (int i = 1; i <= hydrophobic_sections; i++) {
    int atoms_in_si = get_subset(matoms, hyphob, i, si);
    if (atoms_in_si < min_hydrophobic_section_size) {
      continue;
    }

    if (0 == atoms_in_si) {  // huh
      continue;
    }

    for (int j = 1; j <= hydrophillic_sections; j++) {
      int atoms_in_sj = get_subset(matoms, hyphil, j, sj);

      if (atoms_in_sj < min_hydrophillic_section_size) {
        continue;
      }

      int min_separation, max_separation;
      compute_separation(m, si, sj, min_separation, max_separation);

      if (min_separation < min_min_separation) {  // not very interesting it would seem
        continue;
      }

      separation.extra(max_separation);
    }
  }

  if (0 == separation.n()) {
    no_inter_descriptors(inter_descriptors);
    return;
  }

  inter_descriptors.add(static_cast<float>(separation.n()));  // DDnumsep
  inter_descriptors.add(separation.minval());                 // DDminsep
  inter_descriptors.add(
      static_cast<float>(separation.average_if_available_minval_if_not()));  // DDavesep
  inter_descriptors.add(static_cast<float>(separation.maxval()));            // DDmaxsep

  return;
}

static void
do_form_ratio_descriptors(const resizable_array<float>& hydrophobic_descriptors,
                          const resizable_array<float>& hydrophillic_descriptors,
                          resizable_array<float>& ratio_descriptors) {
  int n = hydrophobic_descriptors.number_elements();

  if (n != hydrophillic_descriptors.number_elements()) {
    cerr << "Descriptor count mismatch " << hydrophobic_descriptors.number_elements()
         << " hydrophpbic vs " << hydrophillic_descriptors.number_elements()
         << " hydrophillic\n";
    abort();
  }

  assert(ratio_descriptors.empty());

  ratio_descriptors.resize(n);

  for (int i = 0; i < n; i++) {
    float h1 = hydrophobic_descriptors[i];
    float h2 = hydrophillic_descriptors[i];

    if (float_zero == h2) {
      h1 = float_zero;
    } else {
      h1 = h1 / h2;
    }

    ratio_descriptors.add(h1);
  }

  return;
}

static void
copy_isotopic_info(int matoms, int* isotope, const int* phobphil, int iso) {
  for (int i = 0; i < matoms; i++) {
    if (phobphil[i]) {
      isotope[i] = iso;
    }
  }

  return;
}

static int
arom_count(Molecule& m, int flag, const int* hydrophobic) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (flag != hydrophobic[i]) {
      continue;
    }

    if (!m.is_aromatic(i)) {  // found ! aromatic works better than aromatic()
      rc++;
    }
  }

  return rc;
}

static int
mxdst_count(Molecule& m, int flag, const int* hydrophobic) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (flag != hydrophobic[i]) {
      continue;
    }

    for (int j = 0; j < matoms; j++) {
      if (flag != hydrophobic[j]) {
        continue;
      }

      if (j == i) {
        continue;
      }

      int d = m.bonds_between(i, j);

      if (d > rc) {
        rc = d;
      }
    }
  }

  return rc;
}

static int
pi_count(Molecule& m, int flag, const int* hydrophobic) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (flag != hydrophobic[i]) {
      continue;
    }

    int lp;
    if (m.lone_pair_count(i, lp) && lp > 0) {
      rc++;
    }
  }

  return rc;
}

static int
is_terminal(const Molecule& m, int flag, const int* hydrophobic) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (flag != hydrophobic[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    for (int j = 0; j < a->ncon(); j++) {
      atom_number_t k = a->other(i, j);

      if (hydrophobic[k]) {  // another atom in this section
        continue;
      }

      rc++;  // we have a bonded atom, not part of the hydrobic section

      if (rc > 1) {  // if more than 1 connection, we are not terminal
        return 0;
      }
    }
  }

  return 1;  // whole molecule must be greasy
}

static int
identify_contiguous_section(const Molecule& m, atom_number_t zatom, int flag,
                            int* hydrophobic) {
  hydrophobic[zatom] = flag;

  int rc = 1;

  const Atom* a = m.atomi(zatom);

  int nc = a->ncon();
  for (int i = 0; i < nc; i++) {
    atom_number_t j = a->other(zatom, i);
    if (0 == hydrophobic[j]) {
      continue;
    }

    if (flag == hydrophobic[j]) {  // must be a hydrophobic ring somewhere
      continue;
    }

    rc += identify_contiguous_section(m, j, flag, hydrophobic);
  }

  return rc;
}

static int
identify_contiguous_sections(const Molecule& m, resizable_array<int>& section_size,
                             int* hydrophobic) {
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    if (hydrophobic[i]) {
      hydrophobic[i] = matoms + 1;
    }
  }

  int number_sections = 0;

  for (int i = 0; i < matoms; i++) {
    if (0 == hydrophobic[i]) {
      continue;
    }

    if (matoms + 1 != hydrophobic[i]) {  // already assigned
      continue;
    }

    number_sections++;
    hydrophobic[i] = number_sections;
    int size_of_section = identify_contiguous_section(m, i, number_sections, hydrophobic);

    section_size.add(size_of_section);
  }

  return number_sections;
}

static int
do_write_labelled_smiles(Molecule& m, const int* isotope,
                         Molecule_Output_Object& output) {
  m.set_isotopes(isotope);

  return output.write(m);
}

/*
  checks for atoms that are assigned to both types. As a by-product, we also
  sum the number of atoms that are neither
*/

static int
check_atoms_both_hydrophobic_and_hydrophillic(
    const Molecule& m, const int* hyphil, const int* hyphob,
    int& atoms_neither_hydrophobic_nor_hydrophillic) {
  atoms_neither_hydrophobic_nor_hydrophillic = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (hyphil[i] && hyphob[i]) {
      cerr << "\nFatal error, atom " << i << " " << m.smarts_equivalent_for_atom(i)
           << " both hydrophobic and hydrophillic, molecule '" << m.name() << "'\n";
      return 0;
    } else if (hyphil[i] || hyphob[i]) {
      ;
    } else {
      atoms_neither_hydrophobic_nor_hydrophillic++;
    }
  }

  return 1;
}

static int
atoms_from_queries(Molecule& m, int* hydrophobic,
                   const resizable_array_p<Substructure_Hit_Statistics>& queries) {
  Molecule_to_Match target(&m);

  int queries_matching = 0;

  for (int i = 0; i < queries.number_elements(); i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    queries_matching++;
    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);
      e->set_vector(hydrophobic, 1);
    }
  }

  if (verbose > 2) {
    cerr << queries_matching << " queries matched " << m.name() << "'\n";
  }

  return 1;
}

static int
attached_to_Nitrogen(const Molecule& m, atom_number_t zatom) {
  return 0;  // didn't improve things
  const Atom* a = m.atomi(zatom);

  for (int i = 0; i < a->ncon(); i++) {
    const Bond* b = a->item(i);
    if (!b->is_single_bond()) {
      continue;
    }

    atom_number_t j = b->other(i);
    if (7 == m.atomic_number(j)) {
      return 1;
    }
  }

  return 0;
}

/*
  We have an aromatic carbon with 3 connections. Is it aliphatic?
*/

static int
check_3_connected_aromatic(Molecule& m, atom_number_t c) {
  const Atom* a = m.atomi(c);

  assert(3 == a->ncon());

  for (int i = 0; i < a->ncon(); i++) {
    const Bond* b = a->item(i);
    if (b->is_aromatic()) {
      continue;
    }

    atom_number_t j = b->other(c);

    if (m.is_aromatic(j)) {  // maybe...
      continue;
    }

    const Atom* aj = m.atomi(j);

    if (1 == aj->ncon() && 6 != aj->atomic_number()) {  // phenol, aniline
      return 0;
    }

    int j_attached_heteroatom_count = m.attached_heteroatom_count(j);

    if (6 == aj->atomic_number() && j_attached_heteroatom_count &&
        aj->ncon() < aj->nbonds()) {  // maybe an aldehyde or cyano or amide
      return 0;
    }

    if (7 == aj->atomic_number() && 3 == aj->ncon() && 2 == j_attached_heteroatom_count &&
        5 == aj->nbonds()) {  // probably a nitro, can be hydrophobic
      continue;
    }

    if (6 != aj->atomic_number() && j_attached_heteroatom_count &&
        aj->ncon() < aj->nbonds()) {  // sulphonamide, ...
      return 0;
    }
  }

  return 1;
}

static int
part_of_cf3(Molecule& m, atom_number_t f1, int& aromatic_attachment) {
  const Atom* a1 = m.atomi(f1);

  if (1 != a1->ncon()) {  // seems unlikely
    return 0;
  }

  atom_number_t c = a1->other(f1, 0);

  const Atom* ac = m.atomi(c);

  if (6 != ac->atomic_number()) {
    return 0;
  }

  if (4 != ac->ncon()) {
    return 0;
  }

  if (4 != ac->nbonds()) {
    return 0;
  }

  int fluorine_atoms_found = 0;

  for (int i = 0; i < 4; i++) {
    atom_number_t f = ac->other(c, i);

    if (9 == m.atomic_number(f)) {
      fluorine_atoms_found++;
    } else if (m.is_aromatic(f)) {
      aromatic_attachment = 1;
    } else {
      aromatic_attachment = 0;
    }
  }

  return 3 == fluorine_atoms_found;
}

static int
is_aromatic_nitro(Molecule& m, atom_number_t n, int* tmp, int flag) {
  const Atom* a = m.atomi(n);

  assert(3 == a->ncon());

  atom_number_t aromatic_atom = INVALID_ATOM_NUMBER;
  atom_number_t o1 = INVALID_ATOM_NUMBER;
  atom_number_t o2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < a->ncon(); i++) {
    const Bond* b = a->item(i);

    atom_number_t j = b->other(n);

    if (b->is_double_bond()) {
      if (8 != m.atomic_number(j)) {
        return 0;
      }

      if (INVALID_ATOM_NUMBER == o1) {
        o1 = j;
      } else {
        o2 = j;
      }
    } else if (m.is_aromatic(j) && 6 == m.atomic_number(j))  // n-N(=O)=O is hydrophillic
    {
      if (INVALID_ATOM_NUMBER != aromatic_atom) {
        return 0;
      }

      aromatic_atom = j;
    }
  }

  if (INVALID_ATOM_NUMBER == aromatic_atom) {
    return 0;
  }

  tmp[o1] = flag;
  tmp[o2] = flag;
  tmp[n] = flag;
  tmp[aromatic_atom] = flag;

  return 1;
}

static int
part_of_aromatic_nitro(Molecule& m, atom_number_t zatom, int* tmp, int flag) {
  const Atom* a = m.atomi(zatom);

  atomic_number_t z = a->atomic_number();

  if (8 == z) {
    if (1 != a->ncon()) {
      return 0;
    }

    if (2 != a->nbonds()) {
      return 0;
    }

    atom_number_t n = a->other(zatom, 0);

    if (7 != m.atomic_number(n)) {
      return 0;
    }

    return part_of_aromatic_nitro(
        m, n, tmp, flag);  // inefficient. Fix if this ever becomes a problem
  } else if (7 == z) {
    if (3 != a->ncon()) {
      return 0;
    }

    if (5 != a->nbonds()) {
      return 0;
    }

    return is_aromatic_nitro(m, zatom, tmp, flag);
  } else {
    return 0;
  }

  return 1;
}

static void
do_apply_simple_hydrophillic_rules(Molecule& m, int* hydrophillic) {
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    Atom* a = const_cast<Atom*>(m.atomi(i));

    const atomic_number_t z = a->atomic_number();

    if (6 == z) {  // never hydrophillic
      continue;
    }

    if (17 == z || 35 == z || 53 == z) {
      continue;
    }

    if (9 == z) {
      int aromatomic_attachment;

      if (cf3_is_hydrophillic && part_of_cf3(m, i, aromatomic_attachment) &&
          !aromatomic_attachment) {
        hydrophillic[i] = 1;
      } else {
        continue;
      }
    }

    if (16 == z) {
      if (two_connected_sulphur_is_greasy && 2 == a->ncon() && !m.is_aromatic(i)) {
        continue;
      }
    }

    if (part_of_aromatic_nitro(m, i, hydrophillic, 0)) {
      continue;
    }

    hydrophillic[i] = 1;
  }

  return;
}

static void
do_apply_simple_hydrophobic_rules(Molecule& m, int* hydrophobic) {
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (hydrophobic[i]) {
      continue;
    }

    Atom* a = const_cast<Atom*>(m.atomi(i));

    const atomic_number_t z = a->atomic_number();

    if (6 == z) {
      if (attached_to_Nitrogen(m, i)) {
        ;
      } else if (a->implicit_hydrogens()) {
        hydrophobic[i] = 1;
      } else if (4 == a->ncon()) {
        hydrophobic[i] = 1;
      } else if (m.is_aromatic(i) && 3 == a->ncon() && check_3_connected_aromatic(m, i)) {
        hydrophobic[i] = 1;
      }
    } else if (17 == z || 35 == z || 53 == z) {
      hydrophobic[i] = 1;
    } else if (9 == z) {
      int aromatic_attachment;

      if (cf3_is_hydrophillic && part_of_cf3(m, i, aromatic_attachment) &&
          !aromatic_attachment) {
        continue;
      }

      hydrophobic[i] = 1;
    } else if (16 == z) {
      if (two_connected_sulphur_is_greasy && 2 == a->ncon() && !m.is_aromatic(i)) {
        hydrophobic[i] = 1;
      }
    } else if (part_of_aromatic_nitro(m, i, hydrophobic, 1)) {
      hydrophobic[i] = 1;
    }
  }

  return;
}

static int
assign_based_on_partial_charge(Molecule& m, int* hyphob,
                               charge_t charge_needed_for_hydrophobic, int* hyphil,
                               charge_t charge_needed_for_hydrophillic) {
  m.compute_Abraham_partial_charges();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    charge_t q = fabs(m.charge_on_atom(i));

    if (q < charge_needed_for_hydrophobic) {
      hyphob[i] = 1;
    } else if (q > charge_needed_for_hydrophillic) {
      hyphil[i] = 1;
    }
  }

  return 1;
}

static int
invert_hydrophobic(int matoms, int* hyphil, const int* hyphob) {
  for (int i = 0; i < matoms; i++) {
    if (hyphob[i]) {
      hyphil[i] = 0;
    } else {
      hyphil[i] = 1;
    }
  }

  return 1;
}

static int
identify_hydrophillic(Molecule& m, int* hyphil, const int* hyphob) {
  if (non_hydrophobic_atoms_are_hydrophillic) {
    invert_hydrophobic(m.natoms(), hyphil, hyphob);
  }

  if (apply_simple_hydrophillic_rules) {
    do_apply_simple_hydrophillic_rules(m, hyphil);
  }

  if (hydrophillic_queries.number_elements()) {
    atoms_from_queries(m, hyphil, hydrophillic_queries);
  }

  return 1;
}

static int
identify_hydrophobic_atoms(Molecule& m, int* hydrophobic) {
  if (apply_simple_hydrophobic_rules) {
    do_apply_simple_hydrophobic_rules(m, hydrophobic);
  }

  if (hydrophobic_queries.number_elements()) {
    atoms_from_queries(m, hydrophobic, hydrophobic_queries);
  }

  return count_non_zero_occurrences_in_array(hydrophobic, m.natoms());
}

static int
do_create_descriptors(Molecule& m, int hydrophobic_atoms, int* hydrophobic,
                      resizable_array<float>& descriptors) {
  int matoms = m.natoms();

  float hydrophobic_atom_fraction =
      static_cast<float>(hydrophobic_atoms) / static_cast<float>(matoms);

  descriptors.add(static_cast<float>(hydrophobic_atoms));          // DDnumbr
  descriptors.add(static_cast<float>(hydrophobic_atom_fraction));  // DDratio

  resizable_array<int> section_size;

  int contiguous_sections = identify_contiguous_sections(m, section_size, hydrophobic);

  if (verbose > 2) {
    cerr << "Found " << contiguous_sections << " contiguous sections\n";
  }

  assert(contiguous_sections == section_size.number_elements());

  descriptors.add(static_cast<float>(contiguous_sections));  // DDnsect

  if (0 == contiguous_sections) {
    for (int i = 0; i < 14; i++) {
      descriptors.add(float_zero);
    }

    return 1;
  }

  int terminal_hydrophobic_sections = 0;
  int largest_terminal_hydrophobic_section = 0;

  Accumulator_Int<int> size_acc;
  Accumulator_Int<int> pi_acc;
  Accumulator_Int<int> mxdst_acc;
  Accumulator_Int<int> arom_acc;

  for (int i = 0; i < contiguous_sections; i++) {
    size_acc.extra(section_size[i]);
    if (is_terminal(m, i + 1, hydrophobic)) {
      terminal_hydrophobic_sections++;
      if (section_size[i] > largest_terminal_hydrophobic_section) {
        largest_terminal_hydrophobic_section = section_size[i];
      }
    }

    pi_acc.extra(pi_count(m, i + 1, hydrophobic));
    mxdst_acc.extra(mxdst_count(m, i + 1, hydrophobic));
    arom_acc.extra(arom_count(m, i + 1, hydrophobic));
  }

  descriptors.add(static_cast<float>(size_acc.minval()));  // DDminsz
  descriptors.add(
      static_cast<float>(size_acc.average_if_available_minval_if_not()));  // DDavesz
  descriptors.add(static_cast<float>(size_acc.maxval()));                  // DDmaxsz
  descriptors.add(static_cast<float>(size_acc.maxval()) /
                  static_cast<float>(matoms));  // DDmxntr

  descriptors.add(static_cast<float>(pi_acc.minval()));  // DDminpi
  descriptors.add(
      static_cast<float>(pi_acc.average_if_available_minval_if_not()));  // DDavepi
  descriptors.add(static_cast<float>(pi_acc.maxval()));                  // DDmaxpi
  descriptors.add(static_cast<float>(pi_acc.maxval()) /
                  static_cast<float>(matoms));  // DDmpntr

  descriptors.add(static_cast<float>(mxdst_acc.maxval()));  // DDmxmxd
  descriptors.add(
      static_cast<float>(mxdst_acc.average_if_available_minval_if_not()));  // DDavmxd

  descriptors.add(static_cast<float>(arom_acc.maxval()));  // DDmxarm
  descriptors.add(
      static_cast<float>(arom_acc.average_if_available_minval_if_not()));  // DDavarm

  descriptors.add(static_cast<float>(terminal_hydrophobic_sections));         // DDntrml
  descriptors.add(static_cast<float>(largest_terminal_hydrophobic_section));  // DDlgtrm

  return 1;
}

static int
hydrophobic_sections(Molecule& m, int* hyphil, int* isotope, IWString& output_buffer) {
  int matoms = m.natoms();

  if (verbose > 1) {
    cerr << "Processing '" << m.name() << "' with " << matoms << " atoms\n";
  }

  int* hyphob = new_int(matoms);
  std::unique_ptr<int[]> free_hyphob(hyphob);

  if (partial_charge_specification_present) {
    assign_based_on_partial_charge(m, hyphob, charge_needed_for_hydrophobic, hyphil,
                                   charge_needed_for_hydrophillic);
  }

  identify_hydrophobic_atoms(m, hyphob);

  int hydrophobic_atoms = count_non_zero_occurrences_in_array(hyphob, matoms);
  hydrophobic_atom_count.extra(hydrophobic_atoms);

  if (verbose > 2) {
    cerr << hydrophobic_atoms << " hydrophobic atoms\n";
  }

  if (hydrophobic_isotope) {
    copy_isotopic_info(matoms, isotope, hyphob, hydrophobic_isotope);
  }

  const IWString& mname = m.name();
  if (1 == mname.nwords()) {
    output_buffer << mname;
  } else {
    const_IWSubstring tmp(mname);
    tmp.truncate_at_first(' ');
    output_buffer << tmp;
  }

  resizable_array<float> hydrophobic_descriptors;

  do_create_descriptors(m, hydrophobic_atoms, hyphob, hydrophobic_descriptors);

#ifdef DEBUG_WRITE_DESCRIPTORS
  cerr << "Writing " << hydrophobic_descriptors.number_elements()
       << " hydrophobic descriptors\n";
#endif

  do_write_descriptors(hydrophobic_descriptors, output_buffer);

  if (create_hydrophillic_descriptors) {
    identify_hydrophillic(m, hyphil, hyphob);

    if (hydrophillic_isotope) {
      copy_isotopic_info(matoms, isotope, hyphil, hydrophillic_isotope);
    }

    int atoms_neither_hydrophobic_nor_hydrophillic = 0;
    if (!check_atoms_both_hydrophobic_and_hydrophillic(
            m, hyphil, hyphob, atoms_neither_hydrophobic_nor_hydrophillic)) {
      return 0;
    }

    int hydrophillic_atoms = count_non_zero_occurrences_in_array(hyphil, matoms);
    hydrophillic_atom_count.extra(hydrophillic_atoms);

    if (verbose > 2) {
      cerr << hydrophillic_atoms << " hydrophillic atoms\n";
    }

    if (hydrophillic_atoms) {
      output_buffer << ' '
                    << static_cast<float>(hydrophobic_atoms) /
                           static_cast<float>(hydrophillic_atoms);
    } else {
      output_buffer << " 99";
    }

    resizable_array<float> hydrophillic_descriptors;

    do_create_descriptors(m, hydrophillic_atoms, hyphil, hydrophillic_descriptors);

#ifdef DEBUG_WRITE_DESCRIPTORS
    cerr << "Writing " << hydrophillic_descriptors.number_elements()
         << " hydrophillic descriptors\n";
#endif

    do_write_descriptors(hydrophillic_descriptors, output_buffer);

    resizable_array<float> ratio_descriptors;

    do_form_ratio_descriptors(hydrophobic_descriptors, hydrophillic_descriptors,
                              ratio_descriptors);

#ifdef DEBUG_WRITE_DESCRIPTORS
    cerr << "Writing " << ratio_descriptors.number_elements() << " ratio descriptors\n";
#endif

    do_write_descriptors(ratio_descriptors, output_buffer);

    resizable_array<float> inter_descriptors;

    do_create_inter_descriptors(m, hyphob, hyphil, inter_descriptors);

#ifdef DEBUG_WRITE_DESCRIPTORS
    cerr << "Writing " << inter_descriptors.number_elements() << " inter descriptors\n";
#endif

    do_write_descriptors(inter_descriptors, output_buffer);

    if (atoms_neither_hydrophobic_nor_hydrophillic) {
      float tmp = static_cast<float>(atoms_neither_hydrophobic_nor_hydrophillic) /
                  static_cast<float>(matoms);
      output_buffer << ' ' << atoms_neither_hydrophobic_nor_hydrophillic << ' '
                    << tmp;  // DDNTnhphpi DDNTfnhopi
    } else {
      output_buffer << " 0 0";
    }
  }

  output_buffer << '\n';

  return 1;
}

static void
preprocess(Molecule& m) {
  if (strip_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  return;
}

static int
hydrophobic_sections(data_source_and_type<Molecule>& input,
                     IWString_and_File_Descriptor& output) {
  IWString output_buffer;

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    int matoms = m->natoms();

    if (0 == matoms) {
      cerr << "Skipping empty molecule'" << m->name() << "'\n";
      continue;
    }

    int* hyphil;

    if (create_hydrophillic_descriptors) {
      hyphil = new_int(matoms);
    } else {
      hyphil = nullptr;
    }

    int* isotope;

    if (hydrophobic_isotope || hydrophillic_isotope) {
      isotope = new_int(matoms);
    } else {
      isotope = nullptr;
    }

    int rc = hydrophobic_sections(*m, hyphil, isotope, output);

    if (nullptr != hyphil) {
      delete[] hyphil;
    }

    if (nullptr != isotope) {
      do_write_labelled_smiles(*m, isotope, stream_for_labelled_molecules);
      delete[] isotope;
    }

    output.write_if_buffer_holds_more_than(32768);

    if (0 == rc) {
      return 0;
    }
  }

  return 1;
}

static int
hydrophobic_sections(const char* fname, FileType input_type,
                     IWString_and_File_Descriptor& output) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return hydrophobic_sections(input, output);
}

static void
create_header_for_set(const char* prefix, IWString& buffer) {
  if (buffer.length()) {
    buffer << ' ';
  }

  buffer << prefix << "numbr";
  buffer << ' ' << prefix << "ratio";
  buffer << ' ' << prefix << "nsect";
  buffer << ' ' << prefix << "minsz";
  buffer << ' ' << prefix << "avesz";
  buffer << ' ' << prefix << "maxsz";
  buffer << ' ' << prefix << "mxntr";

  buffer << ' ' << prefix << "minpi";
  buffer << ' ' << prefix << "avepi";
  buffer << ' ' << prefix << "maxpi";
  buffer << ' ' << prefix << "mpntr";

  buffer << ' ' << prefix << "mxmxd";
  buffer << ' ' << prefix << "avmxd";
  buffer << ' ' << prefix << "mxarm";
  buffer << ' ' << prefix << "avarm";

  buffer << ' ' << prefix << "ntrml";
  buffer << ' ' << prefix << "lgtrm";

  return;
}

static int
do_create_header(IWString& buffer) {
  create_header_for_set("hpo_HPO", buffer);

  if (!create_hydrophillic_descriptors) {
    return 1;
  }

  buffer << " hpo_OIratio";

  create_header_for_set("hpo_HPI", buffer);

  create_header_for_set("hpo_HOI", buffer);

  buffer << " hpo_OInumsep";
  buffer << " hpo_OIminsep";
  buffer << " hpo_OIavesep";
  buffer << " hpo_OImaxsep";

  buffer << " hpo_NTnhphpi";
  buffer << " hpo_NTfnhopi";

  return 1;
}

static int
do_write_header(IWString_and_File_Descriptor& output) {
  IWString buffer;

  do_create_header(buffer);

  output << "Name " << buffer << '\n';

  return 1;
}

static int
parse_specification(Command_Line& cl, char flag, int& simple, int& minsize, int& isotope,
                    charge_t& q,
                    resizable_array_p<Substructure_Hit_Statistics>& queries) {
  int i = 0;
  const_IWSubstring o;
  while (cl.value(flag, o, i++)) {
    if (o.starts_with("minsize=")) {
      o.remove_leading_chars(8);
      if (!o.numeric_value(minsize) || minsize < 1) {
        cerr << "The minsize directive must be followed by a whole number '" << o
             << "' is invalid\n";
        return 0;
      }

      if (verbose) {
        cerr << "Minsize for -" << flag << " set to " << minsize << '\n';
      }
    } else if (o.starts_with("charge=")) {
      o.remove_leading_chars(7);
      if (!o.numeric_value(q) || q < static_cast<charge_t>(0.0)) {
        cerr << "Partial charge values must be positive real numbers '" << o
             << "' is invalid\n";
        return 0;
      }

      if (verbose) {
        cerr << "Charge for -" << flag << " set to " << q << '\n';
      }
    } else if (o.starts_with("iso=")) {
      o.remove_leading_chars(4);
      if (!o.numeric_value(isotope) || isotope < 1) {
        cerr << "The isotope specification must be a whole positive number\n";
        return 0;
      }

      if (verbose) {
        cerr << "For -" << flag << " matching atoms written with isotope " << isotope
             << '\n';
      }
    } else if ("simple" == o) {
      simple = 1;
    } else if ("def" == o) {
      simple = 1;
    } else if (o.starts_with("q:")) {
      o.remove_leading_chars(2);
      Substructure_Hit_Statistics* tmp = new Substructure_Hit_Statistics;
      if (!tmp->read(o)) {
        cerr << "Cannot initialise query from '" << o << "'\n";
        delete tmp;
        return 0;
      }

      queries.add(tmp);

      if (verbose) {
        cerr << "Read query from '" << o << "'\n";
      }
    } else if (o.starts_with("Q:")) {
      o.remove_leading_chars(2);
      if (!queries_from_file(o, queries, 1, verbose)) {
        cerr << "Error reading queries from file '" << o << "'\n";
        return 0;
      }
    } else  // assume it is a smarts
    {
      Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;
      if (!q->create_from_smarts(o)) {
        delete q;
        return 0;
      }

      if (verbose) {
        cerr << "Created query from smarts '" << o << "'\n";
      }

      queries.add(q);
    }
  }

  if (verbose) {
    cerr << "Created " << queries.number_elements() << " queries\n";
  }

  return 1;
}

static int
hydrophobic_sections(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:E:nG:L:e:xo:S:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  strip_to_largest_fragment = 1;

  if (!process_standard_aromaticity_options(cl, verbose > 1)) {
    cerr << "Cannot process -A option\n";
    usage(11);
  }

  if (!process_elements(cl, verbose > 1, 'E')) {
    cerr << "Cannot initialise elements\n";
    usage(8);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('G')) {
    if (!parse_specification(cl, 'G', apply_simple_hydrophobic_rules,
                             min_hydrophobic_section_size, hydrophobic_isotope,
                             charge_needed_for_hydrophobic, hydrophobic_queries)) {
      cerr << "Cannot parse hydrophobic specification (-G option)\n";
      return 5;
    }

    if (static_cast<charge_t>(0.0) != charge_needed_for_hydrophobic) {
      partial_charge_specification_present = 1;
    }
  } else {
    apply_simple_hydrophobic_rules = 1;
    if (verbose) {
      cerr << "Simple hydrophobic rules applied\n";
    }
  }

  if (cl.option_present('L')) {
    if (!parse_specification(cl, 'L', apply_simple_hydrophillic_rules,
                             min_hydrophillic_section_size, hydrophillic_isotope,
                             charge_needed_for_hydrophillic, hydrophillic_queries)) {
      cerr << "Cannot parse hydrophillic specification (-L option)\n";
      return 4;
    }

    create_hydrophillic_descriptors = 1;

    if (static_cast<charge_t>(0.0) != charge_needed_for_hydrophillic) {
      partial_charge_specification_present = 1;
    }
  }

  if (cl.option_present('x')) {
    non_hydrophobic_atoms_are_hydrophillic = 1;

    if (verbose) {
      cerr << "All non-hydrophobic atoms are hydrophillic\n";
    }

    create_hydrophillic_descriptors = 1;
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', min_min_separation) || min_min_separation < 1) {
      cerr << "The min separation between sections (-e option) must be a whole positive "
              "number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Hydrophobic sections must be separated by at least " << min_min_separation
           << " bonds\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (hydrophillic_isotope || hydrophobic_isotope) {
    if (hydrophillic_isotope == hydrophobic_isotope) {
      cerr << "Warning, hydrophobic and hydrophillic isotope the same\n";
    }

    if (!cl.option_present('S')) {
      cerr << "Must specify file name stem for file of labelled molecules\n";
      usage(3);
    }

    IWString s = cl.string_value('S');

    if (cl.option_present('o')) {
      if (!stream_for_labelled_molecules.determine_output_types(cl)) {
        cerr << "Cannot parse -o option(s)\n";
        usage(8);
      }
    } else {
      stream_for_labelled_molecules.add_output_type(FILE_TYPE_SMI);
    }

    if (stream_for_labelled_molecules.would_overwrite_input_files(cl, s)) {
      cerr << "Sorry, cannot overwrite input file(s) '" << s << "'\n";
      return 4;
    }

    if (!stream_for_labelled_molecules.new_stem(s)) {
      cerr << "Cannot initialise output stream '" << s << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Labelled molecules written to '" << s << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  if (!do_write_header(output)) {
    cerr << "Cannot write header\n";
    return 6;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!hydrophobic_sections(cl[i], input_type, output)) {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      rc = i + 1;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << "Molecules has between " << hydrophobic_atom_count.minval() << " and "
         << hydrophobic_atom_count.maxval() << " hydrophobic atoms\n";
    if (hydrophobic_atom_count.n()) {
      cerr << "ave " << hydrophobic_atom_count.average_if_available_minval_if_not()
           << '\n';
    }
    cerr << "Molecules has between " << hydrophillic_atom_count.minval() << " and "
         << hydrophillic_atom_count.maxval() << " hydrophillic atoms\n";
    if (hydrophillic_atom_count.n()) {
      cerr << "ave " << hydrophillic_atom_count.average_if_available_minval_if_not()
           << '\n';
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = hydrophobic_sections(argc, argv);

  return rc;
}
