/*
  Sort molecules based on a given criterion
*/

#include <iostream>
#include <memory>

#define IWQSORT_FO_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static int nproperties = 0;

static int *comparison_criterion = nullptr;
static int *direction = nullptr;

// Should be an enum.
#define CMP_NATOMS 0
#define CMP_NRINGS 1
#define CMP_AMW 2
#define CMP_NFRAG 3
#define CMP_NCHIRAL 4
#define CMP_NHETERO 5
#define CMP_AROMATIC_ATOMS 6
#define CMP_AROMATIC_RINGS 7
#define CMP_LARGEST_RING_SIZE 8
#define CMP_LARGEST_RING_SYSTEM_SIZE 9
#define CMP_COLUMN_OF_NAME 10
#define CMP_ROTATBLE_BONDS 11
#define CMP_ATOMS_IN_LARGEST_RING_SYSTEM 12
#define CMP_ATOMS_IN_COUNTERION 13
#define CMP_ATOMS_IN_COUNTERION_MP 14
#define CMP_AMW_IN_COUNTERION 15
#define CMP_AMW_IN_COUNTERION_MP 16
#define CMP_ATOMS_IN_LARGEST_FRAGMENT 17
#define CMP_AMW_IN_LARGEST_FRAGMENT 17
#define CMP_AMW_NOH 18
#define CMP_SDF_TAG 19
#define CMP_NUMBER_ISOTOPIC_ATOMS 20
#define CMP_SUBSTRUCTURE_SEARCH 21
#define CMP_UNSATURATION 22
#define CMP_ATOMS_IN_FIRST_FRAGMENT 23
#define CMP_RING_ATOMS 24
#define CMP_ANY_CHARGE 25
#define CMP_SMALLEST_RING_SIZE 26
#define CMP_ATOMIC_NUMBER_TOTAL 27
#define CMP_SP3 27
#define CMP_NBONDS 28

static extending_resizable_array<int> comparison_column;

static atomic_mass_t amw_hydrogen = static_cast<atomic_mass_t>(1.00794);

static IWString *comparison_tag = nullptr;

static int presumed_atoms_in_counterion_if_no_counterion = 0;

static molecular_weight_t presumed_amw_in_counterion_if_no_counterion = 0.0;

static int reduce_to_largest_fragment = 0;

static Chemical_Standardisation chemical_standardisation;

/*
  By default, we allow just one substructure query. The value is the
  total number of matches across all queries.
  We can also run in a mode where there are any number of substructure
  queries. In that case, we need a cross reference between property
  number and query number
*/

static resizable_array_p<Substructure_Query> queries;

static int treat_queries_as_a_single_group = 0;

static int *property_number_to_query_number = nullptr;

/*
  When writing molecules in chunks, we can set a minimum molecules
  per file we'd like to see.
*/

static int min_size_hint = 10000;

static int splits_to_contain_equal_total_atoms = 0;

static unsigned int total_atom_count = 0;

static Report_Progress report_progress;

static int convert_non_numeric_column_contents = 0;

// By default we create a sequence of file names starting at 1.
// Alternatively we can use the first computed property as the name
// index. Most commonly this will be the heavy atom count.
static int file_name_index_is_first_property = 0;

static void
usage(int rc) {
  cerr << "Sorts a structure file by various criteria\n";
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -k <key(s)>    sort specification(s)\n";
  cerr << "                 a preceeding negative sign indicates reverse order for that property\n";
  cerr << "                 'a' 'natoms' number of atoms\n";
  cerr << "                 'r' 'nring'  number of rings\n";
  cerr << "                 'w' 'amw'    molecular weight\n";
  cerr << "                     'amwnH'  molecular weight, excluding Hydrogens - groups tautomers\n";
  cerr << "                 'f' 'nfrag'  number of fragments\n";
  cerr << "                 'c' 'nchiral' number of explicit chiral centres\n";
  cerr << "                 'h' 'hetero' number of heteroatoms\n";
  cerr << "                 'j' 'aroma'  number of aromatic atoms\n";
  cerr << "                 'k' 'aromr'  number of aromatic rings\n";
  cerr << "                 'R' 'lgrsz'  atoms in largest ring\n";
  cerr << "                 'S' 'lgrss'  rings in largest ring system\n";
  cerr << "                 'b' 'rotbond' number rotatable bonds\n";
  cerr << "                     'asr'    atoms in smallest ring\n";
  cerr << "                     'alrss'  atoms in largest ring system\n";
  cerr << "                     'unsat' number of (non aromatic) unsaturated bonds\n";
  cerr << "                     'ailf'   atoms in largest fragment\n";
  cerr << "                     'aiff'   atoms in first fragment\n";
  cerr << "                     'amwlf'  amw in largest fragment\n";
  cerr << "                     'aicm=<n>' atoms in counterions, missing assigned <n>\n";
  cerr << "                     'amwcm=<x>' amw in counterions, missing assigned <x>\n";
  cerr << "                     'col=nn' numeric contents of column <nn> of the name\n";
  cerr << "                     'niso' number of isotopic atoms\n";
  cerr << "                     'rngat' number of ring atoms\n";
  cerr << "                     'Z' sum of atomic numbers in the molecule\n";
  cerr << "                     'sp3' number of sp3 atoms\n";
  cerr << "                     'nbonds' number of bonds (single=1, double=2, triple=3)\n";
  cerr << "                     'charge' number atoms with formal charges\n";
  cerr << "                     'qry=...' hits to substructure query\n";
  cerr << "                     'smt=...' hits to substructure query\n";
  cerr << "                     'sdf=...' values in SDF tag\n";
  cerr << "  -d             descending order\n";
  cerr << "  -y             if multiple queries present, treat as a group\n";
  cerr << "  -D <stem>      write different groups to output files starting with <stem>\n";
  cerr << "  -e <number>    hint for minimum number of molecules per output file\n";
  cerr << "  -M ...         miscellaneous other options, enter '-M help' for info\n";
  cerr << "  -q             quick exit - avoids overhead for deallocation\n";
  display_standard_aromaticity_options (cerr);
  cerr << "  -E <symbol>    create element with symbol\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

class File_Record {
 private:
  off_t _offset;

  int _nlines;

  float *_property;

 public:
  File_Record();
  ~File_Record();

  int echo(iwstring_data_source &, IWString &) const;

  void set_nlines(int n) {
    _nlines = n;
  }

  int nlines() const {
    return _nlines;
  }

  int initialise(Molecule &, off_t, IW_STL_Hash_Map_int *);

  // The compare function includes the position of the record in the file
  // in order for there to be a stable sort.
  int compare(const File_Record &) const;

  // When writing chunked outputs we only care about computed values for
  // assessing same or different.
  int Same(const File_Record& rhs) const;

  // just assumed to be first property
  int natoms() const {
    return static_cast<int>(_property[0] + 0.01F);
  }; 

  float FirstProperty() const  {
    return _property[0];
  }
};

File_Record::File_Record() {
  if (nproperties) {
    _property = new float[nproperties];
  } else {
    _property = nullptr;
  }

  _nlines = -1;

  return;
}

File_Record::~File_Record() {
  if (nullptr != _property) {
    delete[] _property;
  }

  return;
}

static int
identify_largest_fragment(Molecule &m) {
  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  int largest_fragment = 0;

  int nf = m.number_fragments();

  for (int i = 1; i < nf; i++) {
    if (m.atoms_in_fragment(i) > atoms_in_largest_fragment) {
      atoms_in_largest_fragment = m.atoms_in_fragment(i);
      largest_fragment = i;
    }
  }

  return largest_fragment;
}

static int
common_compute_atoms_in_counterions(Molecule &m, int rc_if_just_one_fragment) {
  if (m.number_fragments() < 2) {
    return rc_if_just_one_fragment;
  }

  int largest_fragment = identify_largest_fragment(m);

  int matoms = m.natoms();

  int *fragment_membership = new int[matoms];
  std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  m.fragment_membership(fragment_membership);

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (fragment_membership[i] != largest_fragment) {
      rc++;
    }
  }

  return rc;
}

static atomic_mass_t
common_compute_amw_in_counterions(Molecule &m, atomic_mass_t rc_if_just_one_fragment) {
  if (m.number_fragments() < 2) {
    return rc_if_just_one_fragment;
  }

  int largest_fragment = identify_largest_fragment(m);

  int matoms = m.natoms();

  int *fragment_membership = new int[matoms];
  std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  m.fragment_membership(fragment_membership);

  atomic_mass_t rc = static_cast<molecular_weight_t>(0.0);

  for (int i = 0; i < matoms; i++) {
    if (fragment_membership[i] == largest_fragment) {
      continue;
    }

    rc += m.atomi(i)->atomic_weight();
    rc += m.implicit_hydrogens(i) * amw_hydrogen;
  }

  return rc;
}

static int
compute_ring_atoms(Molecule &m) {
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_ring_atom(i)) {
      rc++;
    }
  }

  return rc;
}

static int
compute_total_atomic_number(Molecule &m) {
  int rc = m.implicit_hydrogens();

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    rc += m.atomic_number(i);
  }

  return rc;
}

static int
compute_sp3(const Molecule &m) {
  int rc = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (m.ncon(i) == m.nbonds(i)) {
      rc++;
    }
  }

  return rc;
}

// Deliberate decision to NOT specifically handle aromatic bonds.
static int
compute_nbonds(const Molecule &m) {
  int rc = 0;
  const int nedges = m.nedges();

  for (int i = 0; i < nedges; ++i) {
    const Bond *b = m.bondi(i);
    if (b->is_single_bond()) {
      rc += 1;
    } else if (b->is_double_bond()) {
      rc += 2;
    } else if (b->is_triple_bond()) {
      rc += 3;
    }
  }

  return rc;
}

static atomic_mass_t
compute_amw_no_h(Molecule &m) {
  double rc = 0.0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Element *e = m.elementi(i);

    rc += e->atomic_mass();
  }

  return rc;
}

/*
  Because molecular weights could be slightly different depending on
  the order of calculation, we truncate all molecular weights to 3
  significant digits
*/

static float
convert_molecular_weight_to_rounded_form(float amw) {
  int tmp = static_cast<int>(amw * 1000.0 + 0.4999);

  return static_cast<float>(tmp) / 1000.0;
}

static float
compute_molecular_weight(const Molecule &m) {
  return convert_molecular_weight_to_rounded_form(m.molecular_weight_ignore_isotopes());
}

static atomic_mass_t
compute_amw_in_largest_fragment(Molecule &m) {
  int matoms = m.natoms();

  if (m.number_fragments() < 2) {
    return compute_molecular_weight(m);
  }

  int largest_fragment = identify_largest_fragment(m);

  int *fragment_membership = new int[matoms];
  std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  m.fragment_membership(fragment_membership);

  atomic_mass_t rc = static_cast<atomic_mass_t>(0.0);

  for (int i = 0; i < matoms; i++) {
    if (largest_fragment != fragment_membership[i]) {
      continue;
    }

    rc += m.atomi(i)->atomic_weight();
    rc += m.implicit_hydrogens(i) * amw_hydrogen;
  }

  return convert_molecular_weight_to_rounded_form(rc);
  ;
}

static int
compute_atoms_in_first_fragment(Molecule &m) {
  int matoms = m.natoms();

  if (m.number_fragments() < 2) {
    return matoms;
  }

  return m.atoms_in_fragment(0);
}

static int
compute_atoms_in_largest_fragment(Molecule &m) {
  int matoms = m.natoms();

  if (m.number_fragments() < 2) {
    return matoms;
  }

  int largest_fragment = identify_largest_fragment(m);

  int *fragment_membership = new int[matoms];
  std::unique_ptr<int[]> free_fragment_membership(fragment_membership);

  m.fragment_membership(fragment_membership);

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (largest_fragment == fragment_membership[i]) {
      rc++;
    }
  }

  return rc;
}

static int
compute_rotatable_bonds(Molecule &m) {
  int rc = 0;

  m.ring_membership();

  int ne = m.nedges();

  for (int i = 0; i < ne; i++) {
    const Bond *b = m.bondi(i);

    if (b->nrings()) {
      continue;
    }

    if (!b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (m.ncon(a1) > 1 && m.ncon(a2) > 1 && !triple_bond_at_either_end(m, b) &&
        !part_of_otherwise_non_rotabable_entity(m, a1, a2)) {
      rc++;
    }
  }

  return rc;
}

static int
compute_aromatic_atoms(Molecule &m) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      rc++;
    }
  }

  return rc;
}

static float
do_convert_non_numeric_column_contents(const IWString &s,
                                       IW_STL_Hash_Map_int &string_to_ndx) {
  const auto f = string_to_ndx.find(s);
  if (f != string_to_ndx.end()) {
    return static_cast<float>(f->second);
  }

  const auto x = string_to_ndx.size();
  string_to_ndx[s] = x;

  return static_cast<float>(x);
}

static int
fetch_sdf_tag_data(const Molecule &m, const IWString &tag, float &v,
                   IW_STL_Hash_Map_int &string_to_ndx) {
  int n = m.number_records_text_info();

  if (0 == n) {
    cerr << "fetch_sdf_tag_data:no text info in '" << m.name() << "'\n";
    return 0;
  }

  for (int i = 0; i < n; i++) {
    const IWString &info = m.text_info(i);

    if (!info.starts_with("> ")) {
      continue;
    }

    if (!info.contains(tag)) {
      continue;
    }

    assert(i != n - 1);

    const IWString &next_record = m.text_info(i + 1);

    if (next_record.numeric_value(v)) {
      return 1;
    }

    if (convert_non_numeric_column_contents) {
      v = do_convert_non_numeric_column_contents(next_record, string_to_ndx);
      return 1;
    }

    cerr << "fetch_sdf_tag_data:invalid value for '" << tag << "' in '" << m.name()
         << "', value '" << next_record << "'\n";
    return 0;
  }

  cerr << "Could not find tag '" << tag << "' in '" << m.name() << "'\n";
  v = static_cast<float>(0.0);
  return 0;
}

static float
fetch_contents_of_column(const Molecule &m, const int col,
                         IW_STL_Hash_Map_int &token_to_ndx) {
  const IWString &mname = m.name();

  const_IWSubstring token;

  if (!mname.word(col, token)) {
    cerr << "Cannot extract column " << (col + 1) << " from '" << mname << "'\n";
    return -99.0;
  }

  float rc;
  if (token.numeric_value(rc)) {
    return rc;
  }

  if (convert_non_numeric_column_contents) {
    return do_convert_non_numeric_column_contents(token, token_to_ndx);
  }

  cerr << "Invalid numeric in column " << (col + 1) << " of '" << mname << "'\n";
  return -99.9;
}

static int
compute_atoms_in_largest_ring_system(Molecule &m) {
  int matoms = m.natoms();

  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  if (1 == nr) {
    return m.ringi(0)->number_elements();
  }

  int *ring_sys = new int[matoms];
  std::unique_ptr<int[]> free_ring_sys(ring_sys);

  m.label_atoms_by_ring_system(ring_sys);

  int *atoms_in_ring_system = new_int(nr + 1);
  std::unique_ptr<int[]> free_atoms_in_ring_system(atoms_in_ring_system);

  for (int i = 0; i < matoms; i++) {
    int r = ring_sys[i];

    if (0 == r) {  // not a ring atom
      continue;
    }

    if (ring_sys[i] > 0) {
      atoms_in_ring_system[r]++;
    }
  }

  int atoms_in_largest_ring_system = 0;

  for (int i = 1; i <= nr; i++) {
    if (atoms_in_ring_system[i] > atoms_in_largest_ring_system) {
      atoms_in_largest_ring_system = atoms_in_ring_system[i];
    }
  }

  return atoms_in_largest_ring_system;
}

static int
compute_largest_ring_system_size(Molecule &m) {
  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  int rc = 1;

  int *ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring *ri = m.ringi(i);

    if (!ri->is_fused()) {
      continue;
    }

    int system_size = 1;

    for (int j = i + 1; j < nr; j++) {
      if (ring_already_done[j]) {
        continue;
      }

      const Ring *rj = m.ringi(j);

      if (rj->fused_system_identifier() == ri->fused_system_identifier()) {
        ring_already_done[j] = 1;
        system_size++;
      }
    }

    if (system_size > rc) {
      rc = system_size;
    }
  }

  return rc;
}

static int
compute_largest_ring_size(Molecule &m) {
  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  int rc = m.ringi(0)->number_elements();

  for (int i = 1; i < nr; i++) {
    const Ring *ri = m.ringi(i);

    int ring_size = ri->number_elements();

    if (ring_size > rc) {
      rc = ring_size;
    }
  }

  return rc;
}

static int
smallest_ring_size(Molecule &m) {
  if (0 == m.nrings()) {
    return 0;
  }

  return m.ringi(0)->number_elements();
}

static int
compute_aromatic_rings(Molecule &m) {
  int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (int i = 0; i < nr; i++) {
    const Ring *ri = m.ringi(i);

    if (ri->is_aromatic()) {
      rc++;
    }
  }

  return rc;
}

static int
perform_substructure_search(Molecule &m, int property_number) {
  Molecule_to_Match target(&m);

  int rc = 0;

  if (treat_queries_as_a_single_group) {
    for (int i = 0; i < queries.number_elements(); i++) {
      Substructure_Results sresults;
      int nhits = queries[i]->substructure_search(target);
      rc += nhits;
    }
  } else {
    int q = property_number_to_query_number[property_number];
    rc = queries[q]->substructure_search(target);
  }

  return rc;
}

static int
count_usaturation(Molecule &m) {
  m.compute_aromaticity_if_needed();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      continue;
    }

    const Atom *a = m.atomi(i);

    if (a->ncon() < a->nbonds()) {
      rc++;
    }
  }

  return rc;
}

int
File_Record::initialise(Molecule &m, const off_t o, IW_STL_Hash_Map_int *string_to_ndx) {
  _offset = o;

  assert(nproperties > 0);
  if (nullptr == _property) {
    _property = new float[nproperties];
  }

  for (int i = 0; i < nproperties; i++) {
    if (CMP_NATOMS == comparison_criterion[i]) {
      _property[i] = static_cast<float>(m.natoms());
    } else if (CMP_NRINGS == comparison_criterion[i]) {
      _property[i] = static_cast<float>(m.nrings());
    } else if (CMP_AMW == comparison_criterion[i]) {
      _property[i] = static_cast<float>(compute_molecular_weight(m));
    } else if (CMP_AMW_NOH == comparison_criterion[i]) {
      _property[i] = compute_amw_no_h(m);
    } else if (CMP_NFRAG == comparison_criterion[i]) {
      _property[i] = static_cast<float>(m.number_fragments());
    } else if (CMP_NCHIRAL == comparison_criterion[i]) {
      _property[i] = static_cast<float>(m.chiral_centres());
    } else if (CMP_NHETERO == comparison_criterion[i]) {
      _property[i] = static_cast<float>(m.natoms() - m.natoms(6));
    } else if (CMP_AROMATIC_ATOMS == comparison_criterion[i]) {
      _property[i] = static_cast<float>(compute_aromatic_atoms(m));
    } else if (CMP_AROMATIC_RINGS == comparison_criterion[i]) {
      _property[i] = static_cast<float>(compute_aromatic_rings(m));
    } else if (CMP_LARGEST_RING_SIZE == comparison_criterion[i]) {
      _property[i] = static_cast<float>(compute_largest_ring_size(m));
    } else if (CMP_LARGEST_RING_SYSTEM_SIZE == comparison_criterion[i]) {
      _property[i] = static_cast<float>(compute_largest_ring_system_size(m));
    } else if (CMP_COLUMN_OF_NAME == comparison_criterion[i]) {
      _property[i] = static_cast<float>(
          fetch_contents_of_column(m, comparison_column[i], string_to_ndx[i]));
    } else if (CMP_ROTATBLE_BONDS == comparison_criterion[i]) {
      _property[i] = compute_rotatable_bonds(m);
    } else if (CMP_ATOMS_IN_LARGEST_RING_SYSTEM == comparison_criterion[i]) {
      _property[i] = compute_atoms_in_largest_ring_system(m);
    } else if (CMP_ATOMS_IN_COUNTERION_MP == comparison_criterion[i]) {
      _property[i] = common_compute_atoms_in_counterions(
          m, presumed_atoms_in_counterion_if_no_counterion);
    } else if (CMP_AMW_IN_COUNTERION_MP == comparison_criterion[i]) {
      _property[i] = common_compute_amw_in_counterions(
          m, presumed_amw_in_counterion_if_no_counterion);
    } else if (CMP_ATOMS_IN_LARGEST_FRAGMENT == comparison_criterion[i]) {
      _property[i] = compute_atoms_in_largest_fragment(m);
    } else if (CMP_ATOMS_IN_FIRST_FRAGMENT == comparison_criterion[i]) {
      _property[i] = compute_atoms_in_first_fragment(m);
    } else if (CMP_AMW_IN_LARGEST_FRAGMENT == comparison_criterion[i]) {
      _property[i] = compute_amw_in_largest_fragment(m);
    } else if (CMP_NUMBER_ISOTOPIC_ATOMS == comparison_criterion[i]) {
      _property[i] = m.number_isotopic_atoms();
    } else if (CMP_SUBSTRUCTURE_SEARCH == comparison_criterion[i]) {
      _property[i] = perform_substructure_search(m, i);
    } else if (CMP_UNSATURATION == comparison_criterion[i]) {
      _property[i] = count_usaturation(m);
    } else if (CMP_SDF_TAG == comparison_criterion[i]) {
      if (!fetch_sdf_tag_data(m, comparison_tag[i], _property[i], string_to_ndx[i])) {
        return 0;
      }
    } else if (CMP_RING_ATOMS == comparison_criterion[i]) {
      _property[i] = compute_ring_atoms(m);
    } else if (CMP_ATOMIC_NUMBER_TOTAL == comparison_criterion[i]) {
      _property[i] = compute_total_atomic_number(m);
    } else if (CMP_SP3 == comparison_criterion[i]) {
      _property[i] = compute_sp3(m);
    } else if (CMP_NBONDS == comparison_criterion[i]) {
      _property[i] = compute_nbonds(m);
    } else if (CMP_ANY_CHARGE == comparison_criterion[i]) {
      _property[i] = m.number_formally_charged_atoms();
    } else if (CMP_SMALLEST_RING_SIZE == comparison_criterion[i]) {
      _property[i] = smallest_ring_size(m);
    } else if (CMP_SMALLEST_RING_SIZE == comparison_criterion[i]) {
      _property[i] = smallest_ring_size(m);
    } else {
      cerr << "Unrecognised property '" << comparison_criterion[i] << endl;
      return 0;
    }

    if (-1 == direction[i]) {
      _property[i] = static_cast<float>(0.0) - _property[i];
    }
  }

// #define ECHO_COMPUTED_PROPERTIES
#ifdef ECHO_COMPUTED_PROPERTIES
  for (int i = 0; i < nproperties; i++) {
    cerr << " i = " << i << " property " << _property[i] << endl;
  }
#endif

  return 1;
}

int
File_Record::compare(const File_Record &rhs) const {
  for (int i = 0; i < nproperties; i++) {
    //cerr << " i = " << i << " com " << _property[i] << " with " << rhs._property[i] << '\n';
    if (_property[i] < rhs._property[i]) {
      return 1;
    } else if (_property[i] > rhs._property[i]) {
      return -1;
    }
  }

  // Make the sort stable.
  if (_offset < rhs._offset) {
    return -1;
  } else {
    return 1;
  }

  return 0;
}

int
File_Record::Same(const File_Record& rhs) const {
  for (int i = 0; i < nproperties; i++) {
    //cerr << " i = " << i << " com " << _property[i] << " with " << rhs._property[i] << '\n';
    if (_property[i] < rhs._property[i]) {
      return 0;
    } else if (_property[i] > rhs._property[i]) {
      return 0;
    }
  }

  return 1;
}

int
File_Record::echo(iwstring_data_source &input, IWString &output) const {
  assert(_nlines > 0);

  // cerr << "Seeking to " << _offset << ", will write " << _nlines << " lines\n";
  if (!input.seekg(_offset)) {
    cerr << "File_Record::echo: cannot seek to " << _offset << endl;
    return 0;
  }

  const_IWSubstring buffer;

  for (int i = 0; i < _nlines; i++) {
    input.next_record(buffer);
    output << buffer << '\n';
  }

  return 1;
}

// #ifdef __GNUG__
template class resizable_array_p<File_Record>;
template class resizable_array_base<File_Record *>;

// #endif

int
file_record_comparitor_ascending(File_Record *const *cv1, File_Record *const *cv2) {
  const File_Record *pfr1 = *cv1;
  const File_Record *pfr2 = *cv2;

  return pfr2->compare(*pfr1);
}

int
file_record_comparitor_descending(File_Record *const *cv1, File_Record *const *cv2) {
  const File_Record *pfr1 = *cv1;
  const File_Record *pfr2 = *cv2;

  return pfr1->compare(*pfr2);
}

// Group all values that are the same into the same file.
static int
identify_next_group_same_value_same_file(const resizable_array<File_Record *> &records,
                        int &istart, int &istop) {
  const auto number_molecules = records.number_elements();

  if (istop >= number_molecules) {
    return 0;
  }

  istart = istop;

  const File_Record *reff = records[istart];

  istop = istart + 1;

  while (istop < number_molecules) {
    // If the same, continue accumulating.
    if (reff->Same(*(records[istop]))) {
      istop++;
      continue;
    }

    return 1;
  }

  return 1;
}

/*
  This is complicated by our desire to not create too many files, so we
  need to group some things together that really should be separate
*/

static int
identify_next_group(const resizable_array<File_Record *> &records, int &istart,
                    int &istop) {
  // cerr << "identify_next_group:on entry istart " << istart << " istop " << istop << ",
  // min_size_hint " << min_size_hint << endl;
  const auto number_molecules = records.number_elements();

  if (istop >= number_molecules) {
    return 0;
  }

  istart = istop;

  const File_Record *reff = records[istart];

  int last_break = -1;

  istop = istart + 1;

  while (istop < number_molecules) {
    // same, just continue accumulating
    if (reff->Same(*(records[istop]))) {
      istop++;
      continue;
    }

    //  We are at a discontinuity

    //  cerr << "At sorting discontinuity istart " << istart << " istop " << istop <<
    //  endl;

    if (istop - istart < min_size_hint) {  // not enough to write
      ;
    // in the first chunk of molecules, we'll check the next group too
    } else if (last_break < 0) {
      ;
    // too much, write previous chunk
    } else if (istop - istart > (min_size_hint + min_size_hint)) {
      istop = last_break;
      return 1;
    // otherwise we get problems when just two chunks migh be written
    } else if (0 == istart) {
      istop++;
      return 1;
    }

    reff = records[istop];
    last_break = istop;
    istop++;
  }

  return 1;
}

static int
echo_set_of_molecules(const resizable_array<File_Record *> &records, int istart,
                      int istop, data_source_and_type<Molecule> &input,
                      IWString_and_File_Descriptor &output) {
  assert(istart <= istop);

  for (int i = istart; i < istop; i++) {
    if (!records[i]->echo(input, output)) {
      cerr << "Cannot write item " << i << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

static int
write_group(const resizable_array<File_Record *> &molecules, int istart, int istop,
            data_source_and_type<Molecule> &input, const IWString &stem, int ndx) {
  IWString fname;

  fname << stem << ndx << ".smi";

  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open group file '" << fname << "'\n";
    return 0;
  }

  return echo_set_of_molecules(molecules, istart, istop, input, output);
}

static int
identify_next_equi_atom_group(const resizable_array<File_Record *> &records, int &istart,
                              int &istop, int atoms_per_chunk, unsigned int &atoms) {
  // cerr << "identify_next_group:on entry istart " << istart << " istop " << istop << ",
  // min_size_hint " << min_size_hint << endl;

  const auto number_molecules = records.number_elements();

  if (istop >= number_molecules) {
    return 0;
  }

  istart = istop;

  const File_Record *reff = records[istart];

  istop = istart + 1;

  int atoms_this_chunk = 0;

  while (istop < number_molecules) {
    atoms_this_chunk += records[istop]->natoms();

    // same, just continue accumulating
    if (reff->Same(*(records[istop]))) {
      istop++;
      continue;
    }

    //  We are at a discontinuity

    if (atoms_this_chunk >= atoms_per_chunk) {
      return 1;
    }

    reff = records[istop];
    istop++;
  }

  return 1;
}

static int
do_same_score_chunked_output(const resizable_array<File_Record *> &molecules,
                                   data_source_and_type<Molecule> &input,
                                   const IWString &stem_for_output_files) {
  int ndx = 1;

  int istart = 0;
  int istop = 0;

  while (identify_next_group_same_value_same_file(molecules, istart, istop)) {
    if (verbose) {
      cerr << "Writing from " << istart << " to " << istop << endl;
    }

    if (file_name_index_is_first_property) {
      ndx = static_cast<int>(molecules[istart]->FirstProperty() + 0.001);
    }
    if (!write_group(molecules, istart, istop, input, stem_for_output_files, ndx)) {
      return 0;
    }

    ndx++;
  }

  return 1;
}

static int
do_equal_atom_count_chunked_output(const resizable_array<File_Record *> &molecules,
                                   data_source_and_type<Molecule> &input,
                                   const IWString &stem_for_output_files) {
  int atoms_per_chunk = total_atom_count / splits_to_contain_equal_total_atoms;

  if (verbose) {
    cerr << "Chunks will contain roughly " << atoms_per_chunk << " atoms\n";
  }

  int ndx = 1;

  unsigned int atoms = 0;

  int istart = 0;
  int istop = 0;

  while (identify_next_equi_atom_group(molecules, istart, istop, atoms_per_chunk, atoms)) {
    if (verbose) {
      cerr << "Writing from " << istart << " to " << istop << endl;
    }

    if (!write_group(molecules, istart, istop, input, stem_for_output_files, ndx)) {
      return 0;
    }

    ndx++;
  }

  return 1;
}

static int
do_chunked_output(const resizable_array<File_Record *> &molecules,
                  data_source_and_type<Molecule> &input,
                  const IWString &stem_for_output_files) {
  int ndx = 1;

  int istart = 0;
  int istop = 0;
  while (identify_next_group(molecules, istart, istop)) {
    if (verbose) {
      cerr << "Writing from " << istart << " to " << istop << endl;
    }

    if (!write_group(molecules, istart, istop, input, stem_for_output_files, ndx)) {
      return 0;
    }

    ndx++;
  }

  if (verbose) {
    cerr << "Created " << (ndx - 1) << " chunked output files\n";
  }

  return 1;
}

static void
preprocess(Molecule &m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
fetch_queries(const_IWSubstring &q) {
  if (!process_cmdline_token(' ', q, queries, verbose)) {
    cerr << "Cannot process query specification '" << q << "'\n";
    return 0;
  }

  return queries.number_elements();
}

static int
fetch_smarts(const const_IWSubstring &s) {
  Substructure_Query *q = new Substructure_Query;
  if (!q->create_from_smarts(s)) {
    cerr << "Invalid smarts '" << s << "'\n";
    return 0;
  }

  queries.add(q);

  return queries.number_elements();
}

static int
display_dash_m_options(std::ostream &os) {
  // clang-format off
  os << " -M ea=<n>         write <n> chunks of output each containing roughly equal total atoms\n";
  os << " -M s2f            convert non numeric tags and column contents to arbitrary\n";
  os << "                   numeric form - grouped rather than sorted\n";
  os << " -M svsf           same values to in the same file (use with the -D option)\n";
  os << " -M firstprop      the file name index of each -D file is the first property\n";
  os << " -M rpt=<n>        report progress every <n> molecules processed\n";
  // clang-format on

  exit(1);
}

static int
msort(int argc, char **argv) {
  Command_Line cl(argc, argv, "vA:E:i:mfcarhonjRSbzZk:dylg:D:e:qM:s:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose)) {
      cerr << "Cannot process standard aromaticity options\n";
      usage(4);
    }
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl)) {
      cerr << "Cannot initialise elements (-E option)\n";
      usage(3);
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
      cerr << "Will reduce molecules to their largest fragment\n";
    }
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

  int splits_contain_equal_values = 0;

  if (cl.option_present('M')) {
    const_IWSubstring m;

    for (int i = 0; cl.value('M', m, i); i++) {
      if (m.starts_with("ea=")) {
        if (!cl.option_present('D')) {
          cerr << "The equal atom count chunk directive (-M ea) only makes sense with "
                  "the -D option\n";
          usage(2);
        }

        m.remove_leading_chars(3);

        if (!m.numeric_value(splits_to_contain_equal_total_atoms) ||
            splits_to_contain_equal_total_atoms < 2) {
          cerr << "The number of equal atom count splits must be a whole +ve number\n";
          exit(3);
        }

        if (verbose) {
          cerr << "Will write " << splits_to_contain_equal_total_atoms
               << " splits containing roughly equal total atom counts\n";
        }
      } else if (m.starts_with("rpt=")) {
        m.remove_leading_chars(4);
        int r;
        if (!m.numeric_value(r) || r < 1) {
          cerr << "Invalid report interval '" << m << "'\n";
          exit(2);
        }

        report_progress.set_report_every(r);

        if (verbose) {
          cerr << "Will report progress every " << r << " molecules processed\n";
        }
      } else if ("s2f" == m) {
        convert_non_numeric_column_contents = 1;
        if (verbose) {
          cerr << "Will convert non numeric tags and column contents to arbitrary "
                  "numbers\n";
        }
      } else if (m == "svsf") {
        splits_contain_equal_values = 1;
        if (verbose) {
          cerr << "The molecules in each file will have the same computed values\n";
        }
      } else if (m == "firstprop") {
        file_name_index_is_first_property = 1;
        if (verbose) {
          cerr << "The -D files will be created with name containing the first property\n";
        }
      } else if ("help" == m) {
        display_dash_m_options(cerr);
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_m_options(cerr);
      }
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, don't know how to do multiple files\n";
    return 8;
  }

  comparison_criterion = new int[100];  // definite overestimate
  std::unique_ptr<int[]> free_comparison_criterion(comparison_criterion);
  direction = new int[100];

  int number_query_specifications = 0;

  if (cl.option_present('k')) {
    if (nproperties) {
      cerr << "Sorry, the -k option cannot be combined with other options\n";
      return 3;
    }

    int nk = cl.option_count('k');

    comparison_tag = new IWString[nk];
    property_number_to_query_number = new_int(nk, -1);

    const_IWSubstring k;
    int i = 0;
    while (cl.value('k', k, i++)) {
      int j = 0;
      const_IWSubstring token;

      while (k.nextword(token, j, ',')) {
        if (token.starts_with('-')) {
          direction[nproperties] = -1;
          token.remove_leading_chars(1);
        } else if (cl.option_present('r')) {
          direction[nproperties] = -1;
        } else {
          direction[nproperties] = 1;

          if (token.starts_with('+')) {
            token.remove_leading_chars(1);
          }
        }

        if ('a' == token || token.starts_with("atom") || "natoms" == token) {
          comparison_criterion[nproperties] = CMP_NATOMS;
        } else if ('r' == token || token.starts_with("nring")) {
          comparison_criterion[nproperties] = CMP_NRINGS;
        } else if ('w' == token || "amw" == token) {
          comparison_criterion[nproperties] = CMP_AMW;
        } else if ("amwnH" == token) {
          comparison_criterion[nproperties] = CMP_AMW_NOH;
        } else if ('f' == token || "nfrag" == token) {
          comparison_criterion[nproperties] = CMP_NFRAG;
        } else if ('c' == token || "nchiral" == token) {
          comparison_criterion[nproperties] = CMP_NCHIRAL;
        } else if ('h' == token || "hetero" == token) {
          comparison_criterion[nproperties] = CMP_NHETERO;
        } else if ('j' == token || "aroma" == token) {
          comparison_criterion[nproperties] = CMP_AROMATIC_ATOMS;
        } else if ('k' == token || "aromr" == token) {
          comparison_criterion[nproperties] = CMP_AROMATIC_RINGS;
        } else if ('R' == token || "lgrsz" == token) {
          comparison_criterion[nproperties] = CMP_LARGEST_RING_SIZE;
        } else if ('S' == token || "lgrss" == token) {
          comparison_criterion[nproperties] = CMP_LARGEST_RING_SYSTEM_SIZE;
        } else if ("alrss" == token) {
          comparison_criterion[nproperties] = CMP_ATOMS_IN_LARGEST_RING_SYSTEM;
        } else if ("asr" == token) {
          comparison_criterion[nproperties] = CMP_SMALLEST_RING_SIZE;
        } else if ('b' == token || "rotbond" == token) {
          comparison_criterion[nproperties] = CMP_ROTATBLE_BONDS;
        } else if ("ailf" == token) {
          comparison_criterion[nproperties] = CMP_ATOMS_IN_LARGEST_FRAGMENT;
        } else if ("aiff" == token) {
          comparison_criterion[nproperties] = CMP_ATOMS_IN_FIRST_FRAGMENT;
        } else if ("unsat" == token) {
          comparison_criterion[nproperties] = CMP_UNSATURATION;
        } else if ("amwlf" == token) {
          comparison_criterion[nproperties] = CMP_AMW_IN_LARGEST_FRAGMENT;
        } else if ("niso" == token) {
          comparison_criterion[nproperties] = CMP_NUMBER_ISOTOPIC_ATOMS;
        } else if ("rngat" == token) {
          comparison_criterion[nproperties] = CMP_RING_ATOMS;
        } else if ("Z" == token) {
          comparison_criterion[nproperties] = CMP_ATOMIC_NUMBER_TOTAL;
        } else if ("sp3" == token) {
          comparison_criterion[nproperties] = CMP_SP3;
        } else if ("nbonds" == token) {
          comparison_criterion[nproperties] = CMP_NBONDS;
        } else if ("charge" == token) {
          comparison_criterion[nproperties] = CMP_ANY_CHARGE;
        } else if (token.starts_with("qry=")) {
          if (queries.number_elements()) {
            cerr << "Sorry, can have only one substructure search criterion\n";
            return 4;
          }
          token.remove_leading_chars(4);
          if (!fetch_queries(token)) {
            cerr << "Invalid query specification(s) '" << token << "'\n";
            return 4;
          }
          comparison_criterion[nproperties] = CMP_SUBSTRUCTURE_SEARCH;
          property_number_to_query_number[nproperties] = number_query_specifications;
          number_query_specifications++;
        } else if (token.starts_with("smt=")) {
          token.remove_leading_chars(4);
          if (!fetch_smarts(token)) {
            cerr << "Invalid smarts '" << token << "'\n";
            return 4;
          }
          comparison_criterion[nproperties] = CMP_SUBSTRUCTURE_SEARCH;
          property_number_to_query_number[nproperties] = number_query_specifications;
          number_query_specifications++;
        } else if (token.starts_with("aicm=")) {
          token.remove_leading_chars(5);
          if (!token.numeric_value(presumed_atoms_in_counterion_if_no_counterion) ||
              presumed_atoms_in_counterion_if_no_counterion < 0) {
            cerr << "INvalid presumed atom count in missing counterion '" << token
                 << "'\n";
            usage(4);
          }
          comparison_criterion[nproperties] = CMP_ATOMS_IN_COUNTERION_MP;
        } else if (token.starts_with("amwcm=")) {
          token.remove_leading_chars(6);
          if (!token.numeric_value(presumed_amw_in_counterion_if_no_counterion) ||
              presumed_amw_in_counterion_if_no_counterion < 0.0) {
            cerr << "INvalid presumed AMW in missing counterion '" << token << "'\n";
            usage(5);
          }
          comparison_criterion[nproperties] = CMP_AMW_IN_COUNTERION_MP;
        } else if (token.starts_with("col=")) {
          token.remove_leading_chars(4);
          int col;
          if (!token.numeric_value(col) || col < 1) {
            cerr << "Invalid column specification '" << token << "'\n";
            return 4;
          }

          comparison_criterion[nproperties] = CMP_COLUMN_OF_NAME;
          comparison_column[nproperties] = (col - 1);
        } else if (token.starts_with("sdf=")) {
          token.remove_leading_chars(4);
          comparison_criterion[nproperties] = CMP_SDF_TAG;
          comparison_tag[nproperties] << "<" << token << '>';
          moleculeio::set_read_extra_text_info(1);
        }
        //      else if (token.starts_with("rx="))
        //      {
        //        token.remove_leading_chars(3);
        //        comparison_tag[nproperties] = CMP_RX;
        //      }
        else {
          cerr << "Unrecognised sort directive '" << token << "'\n";
          usage(3);
        }

        nproperties++;
      }
    }
  } else {
    for (int i = 0; i < argc; i++) {
      const_IWSubstring c = argv[i];
      if ("-a" == c) {
        comparison_criterion[nproperties] = CMP_NATOMS;
        nproperties++;
      } else if ("-m" == c) {
        comparison_criterion[nproperties] = CMP_AMW;
        nproperties++;
      } else if ("-f" == c) {
        comparison_criterion[nproperties] = CMP_NFRAG;
        nproperties++;
      } else if ("-c" == c) {
        comparison_criterion[nproperties] = CMP_NCHIRAL;
        nproperties++;
      } else if ("-r" == c) {
        comparison_criterion[nproperties] = CMP_NRINGS;
        nproperties++;
      } else if ("-h" == c) {
        comparison_criterion[nproperties] = CMP_NHETERO;
        nproperties++;
      } else if ("-j" == c) {
        comparison_criterion[nproperties] = CMP_AROMATIC_ATOMS;
        nproperties++;
      } else if ("-n" == c) {
        comparison_criterion[nproperties] = CMP_AROMATIC_RINGS;
        nproperties++;
      } else if ("-R" == c) {
        comparison_criterion[nproperties] = CMP_LARGEST_RING_SIZE;
        nproperties++;
      } else if ("-S" == c) {
        comparison_criterion[nproperties] = CMP_LARGEST_RING_SYSTEM_SIZE;
        nproperties++;
      } else if ("-b" == c) {
        comparison_criterion[nproperties] = CMP_ROTATBLE_BONDS;
        nproperties++;
      } else if (c.starts_with("col=")) {
        c.remove_leading_chars(4);
        int col;
        if (!c.numeric_value(col) || col < 1) {
          cerr << "Invalid column specification '" << c << "'\n";
          return 4;
        }

        comparison_criterion[nproperties] = CMP_COLUMN_OF_NAME;
        comparison_column[nproperties] = (col - 1);
        nproperties++;
      }
    }
  }

  if (0 == nproperties) {
    cerr << "Must specify one or more sort criteria\n";
    usage(3);
  }

  if (splits_to_contain_equal_total_atoms > 0 && CMP_NATOMS != comparison_criterion[0]) {
    cerr << "When doing equal atom splits, the first comparison property must be number "
            "of atoms\n";
    return 3;
  }

  if (cl.option_present('e') && !cl.option_present('D')) {
    cerr << "The -e option only makes sense with the -D option\n";
    usage(3);
  }

  if (number_query_specifications >
      queries.number_elements())  // don't think this could even happen
  {
    cerr << "Specified " << number_query_specifications << " queries, but only "
         << queries.number_elements() << " queries, impossible\n";
    return 7;
  }

  if (cl.option_present('y')) {
    if (0 == number_query_specifications) {
      cerr << "No query specifications entered, -y doesn't make sense\n";
      usage(5);
    }

    if (number_query_specifications > 1)  // cannot do this
    {
      cerr << "Sorry, don't know how to do -y with more than 1 query specification\n";
      usage(4);
    }

    treat_queries_as_a_single_group = 1;

    if (verbose) {
      cerr << "Will treat " << queries.number_elements() << " queries as a group\n";
    }
  }

  if (queries.number_elements() == number_query_specifications) {  // good
    ;
  } else if (treat_queries_as_a_single_group) {  // > 1 queries but just one
                                                 // query_specification
    ;
  } else {
    cerr << "Entered " << number_query_specifications << " query specifications, and "
         << queries.number_elements() << " queries, must treat as a group\n";
    cerr << "Did you mean to use the -y option also?\n";
    usage(4);
  }

  if (verbose) {
    cerr << "Defined " << nproperties << " sort criteria\n";
  }

  for (int i = 0; i < queries.number_elements(); i++) {
    queries[i]->set_find_unique_embeddings_only(1);
  }

  if (0 == input_type) {
    input_type = discern_file_type_from_name(cl[0]);
  }

  data_source_and_type<Molecule> input(input_type, cl[0]);

  if (!input.good()) {
    cerr << "Cannot open input file '" << cl[0] << "'\n";
    return 4;
  }

  File_Record *records = nullptr;
  int items_in_file = 0;

  if (cl.option_present('s')) {
    if (!cl.value('s', items_in_file) || items_in_file < 2) {
      cerr << "The -s option (records in input) must be a whole +ve number\n";
      usage(2);
    }

  } else {
    items_in_file = input.records_remaining();

    if (0 == items_in_file) {
      cerr << "No structures in input\n";
      return 2;
    }
  }

  IW_STL_Hash_Map_int *string_to_ndx = new IW_STL_Hash_Map_int[nproperties];
  std::unique_ptr<IW_STL_Hash_Map_int[]> free_string_to_ndx(string_to_ndx);

  records = new File_Record[items_in_file];
  std::unique_ptr<File_Record[]> free_records(records);

  off_t offset = static_cast<off_t>(0);
  int lines_read = 0;
  int molecules = 0;

  Molecule *m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (splits_to_contain_equal_total_atoms) {
      total_atom_count += m->natoms();
    }

    if (!records[molecules].initialise(*m, offset, string_to_ndx)) {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 4;
    }

    records[molecules].set_nlines(input.lines_read() - lines_read);

    molecules++;

    if (molecules >= items_in_file) {  // avoid catastrophe
      break;
    }

    lines_read = input.lines_read();

    offset = input.tellg();

    if (report_progress()) {
      cerr << "Read " << molecules << " molecules\n";
    }
  }

  if (verbose) {
    cerr << "Read " << molecules << " molecules\n";
  }

#ifdef MSORT_USE_ARRAY
  File_Record **frp = new File_Record *[molecules];
  std::unique_ptr<File_Record *[]> free_frp(frp);
#else
  resizable_array<File_Record *> frp(molecules);
#endif

  for (auto i = 0; i < molecules; ++i) {
    frp.add(records + i);
#ifdef MSORT_USE_ARRAY
//  frp[i] = records + i;
#endif
  }

  if (cl.option_present('d')) {
    frp.sort(file_record_comparitor_descending);
  } else {
    frp.sort(file_record_comparitor_ascending);
  }

#ifdef MSORT_USE_ARRAY
  if (cl.option_present('d')) {
    iwqsort(frp, molecules, file_record_comparitor_descending);
  } else {
    iwqsort(frp, molecules, file_record_comparitor_ascending);
  }
#endif

  if (verbose) {
    cerr << "Sort complete\n";
  }

  if (cl.option_present('D')) {
    IWString stem_for_output_files = cl.option_value('D');

    if (cl.option_present('e')) {
      if (!cl.value('e', min_size_hint) || min_size_hint < 1) {
        cerr << "The minimum chunk size hint (-e) must be a value +ve integer\n";
        usage(3);
      }

      if (verbose) {
        cerr << "Will try to create chunked output file(s) with at least "
             << min_size_hint << " molecules in each\n";
      }
    }

    if (splits_to_contain_equal_total_atoms) {
      do_equal_atom_count_chunked_output(frp, input, stem_for_output_files);
    } else if (splits_contain_equal_values) {
      do_same_score_chunked_output(frp, input, stem_for_output_files);
    } else {
      do_chunked_output(frp, input, stem_for_output_files);
    }
  } else {
    IWString_and_File_Descriptor output(1);

    if (!echo_set_of_molecules(frp, 0, molecules, input, output)) {
      cerr << "Could not write\n";
      return 3;
    }

    output.flush();
  }

  if (cl.option_present('q')) {
    _exit(0);
  }

  return 0;
}

int
main(int argc, char **argv) {
  int rc = msort(argc, argv);

  return rc;
}
