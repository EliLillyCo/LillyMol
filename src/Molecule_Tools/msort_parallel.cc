/*
  Sort molecules based on a given criterion
  TBB version to initialise the molecular properties
  in parallel
  TODO:ianwatson. Rationalise this and the serial version to share code.
*/

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdlib.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>

//#include "tbb/atomic.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/scalable_allocator.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwqsort/iwqsort_tbb.h"

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

static int keep_temporary_files = 0;

static int nproperties = 0;

static int quick_set_of_properties;

static int* comparison_criterion = nullptr;
static int* direction = nullptr;

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

static extending_resizable_array<int> comparison_column;

static atomic_mass_t amw_hydrogen = static_cast<atomic_mass_t>(1.00794);

static int ascending_order = 1;

static int initial_size_estimate = 1000000;

static IWString* comparison_tag = nullptr;

static int presumed_atoms_in_counterion_if_no_counterion = 0;

static molecular_weight_t presumed_amw_in_counterion_if_no_counterion = 0.0;

static int reduce_to_largest_fragment = 0;

static int scale_first_chunk = 100;

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

static int* property_number_to_query_number = nullptr;

static int do_parallel_sort = 0;

static IWString forked_process_stem("MSORTFRK");

static IWString stem_for_chunked_output;

static Report_Progress file_scope_report_progress;

/*
  When writing molecules in chunks, we can set a minimum molecules
  per file we'd like to see.
*/

static int min_size_hint = 1000;

static int splits_to_contain_equal_total_atoms = 0;

// By default we create a sequence of file names starting at 1.
// Alternatively we can use the first computed property as the name
// index. Most commonly this will be the heavy atom count.
static int file_name_index_is_first_property = 0;

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
  cerr << "Sorts a structure file by various criteria\n";
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "  -a             sort by number of atoms\n";
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
  cerr << "                 'e' 'alrss'  atoms in largest ring system\n";
  cerr << "                 'b' 'rotbond' number rotatable bonds\n";
  cerr << "                     'unsat' number of (non aromatic) unsaturated bonds\n";
  cerr << "                     'ailf'   atoms in largest fragment\n";
  cerr << "                     'aiff'   atoms in first fragment\n";
  cerr << "                     'amwlf'  amw in largest fragment\n";
  cerr << "                     'aicm=<n>' atoms in counterions, missing assigned <n>\n";
  cerr << "                     'amwcm=<x>' amw in counterions, missing assigned <x>\n";
  cerr << "                     'col=nn' numeric contents of column <nn> of the name\n";
  cerr << "                     'niso' number of isotopic atoms\n";
  cerr << "                     'rngat' number of ring atoms\n";
  cerr << "                     'charge' number atoms with formal charges\n";
  cerr << "                     'Z' sum of atomic numbers\n";
//cerr << "                     'qry=...' hits to substructure query\n";

//cerr << "                     'smt=...' hits to substructure query\n";
  cerr << "  -d             descending order\n";
  cerr << "  -y             if multiple queries present, treat as a group\n";
  cerr << "  -D <stem>      write different groups to output files starting with <stem>\n";
  cerr << "  -e <number>    hint for minimum number of molecules per output file\n";
  cerr << "  -q             quick exit - avoids overhead for deallocation\n";
  cerr << "  -h <nthreads>  use Posix threads with <nthreads> threads\n";
  cerr << "  -f             forked process for each input file (recommended)\n";
  cerr << "  -p             use parallel sort for final sort (recommended)\n";
  cerr << "  -Y <x>         decrease size of chunk with highest atom count to <x> pct of rest\n";
  cerr << "  -M ...         other options, enter '-M help' for info\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -E <symbol>    create element with symbol\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

class File_Record
{
 private:
  IWString _smiles;
  IWString _id;

  float* _property;

 public:
  File_Record();
  ~File_Record();

  int echo(IWString&) const;

  const IWString& id() const
  {
    return _id;
  }

  int set_molecule(Molecule*);
  int set_molecule(const const_IWSubstring&);

  int compute_properties();

  // The serial version of this tool considers the offset. This version
  // does not. But for better code similarity we add a Same method.
  int compare(const File_Record&) const;

  // When writing chunked outputs we only care about computed values for
  // assessing same or different.
  int Same(const File_Record& rhs) const;

  const float* zdata() const {
    return _property;
  }

  // just assumed to be the first property.
  int natoms() const {
    return static_cast<int>(_property[0] + 0.01F);
  }

  float FirstProperty() const  {
    return _property[0];
  }
};

File_Record::File_Record()
{
  if (nproperties) {
    _property = new float[nproperties];
  } else {
    _property = nullptr;
  }

  return;
}

File_Record::~File_Record()
{
  if (nullptr != _property) {
    delete[] _property;
  }

  return;
}

static int
identify_largest_fragment(Molecule& m)
{
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
common_compute_atoms_in_counterions(Molecule& m, int rc_if_just_one_fragment)
{
  if (m.number_fragments() < 2) {
    return rc_if_just_one_fragment;
  }

  int largest_fragment = identify_largest_fragment(m);

  int matoms = m.natoms();

  int* fragment_membership = new int[matoms];
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
common_compute_amw_in_counterions(Molecule& m, atomic_mass_t rc_if_just_one_fragment)
{
  if (m.number_fragments() < 2) {
    return rc_if_just_one_fragment;
  }

  int largest_fragment = identify_largest_fragment(m);

  int matoms = m.natoms();

  int* fragment_membership = new int[matoms];
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
compute_ring_atoms(Molecule& m)
{
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
compute_total_atomic_number(Molecule& m)
{
  int rc = m.implicit_hydrogens();

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    rc += m.atomic_number(i);
  }

  return rc;
}

/*static int
compute_absolute_formal_charge (const Molecule & m)
{
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    formal_charge_t fc = m.formal_charge(i);

    if (0 == fc)
      ;
    else if (fc < 0)
      rc -= fc;
    else
      rc += fc;
  }

//cerr << "Rc for " << m.name() << " is " << rc << endl;

  return rc;
}*/

static atomic_mass_t
compute_amw_no_h(Molecule& m)
{
  Molecular_Weight_Control mwc;
  Molecular_Weight_Calculation_Result mwcr;

  mwc._ignore_hydrogens = true;

  m.molecular_weight(mwc, mwcr);

  return mwcr._amw;
}

/*
  Because molecular weights could be slightly different depending on
  the order of calculation, we truncate all molecular weights to 3
  significant digits
*/

static float
convert_molecular_weight_to_rounded_form(double amw)
{
  int tmp = static_cast<int>(amw * 1000.0 + 0.4999);

  return static_cast<float>(tmp) / 1000.0;
}

static float
compute_molecular_weight(const Molecule& m)
{
  Molecular_Weight_Control mwc;
  mwc._ignore_isotopes = true;

  Molecular_Weight_Calculation_Result mwcr;

  m.molecular_weight(mwc, mwcr);

  return convert_molecular_weight_to_rounded_form(mwcr._amw);
}

static atomic_mass_t
compute_amw_in_largest_fragment(Molecule& m)
{
  int matoms = m.natoms();

  if (m.number_fragments() < 2) {
    return compute_molecular_weight(m);
  }

  int largest_fragment = identify_largest_fragment(m);

  int* fragment_membership = new int[matoms];
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
compute_atoms_in_first_fragment(Molecule& m)
{
  int matoms = m.natoms();

  if (m.number_fragments() < 2) {
    return matoms;
  }

  return m.atoms_in_fragment(0);
}

static int
compute_atoms_in_largest_fragment(Molecule& m)
{
  int matoms = m.natoms();

  if (m.number_fragments() < 2) {
    return matoms;
  }

  int largest_fragment = identify_largest_fragment(m);

  int* fragment_membership = new int[matoms];
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
compute_rotatable_bonds(Molecule& m)
{
  int rc = 0;

  m.ring_membership();

  int ne = m.nedges();

  for (int i = 0; i < ne; i++) {
    const Bond* b = m.bondi(i);

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
compute_aromatic_atoms(Molecule& m)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      rc++;
    }
  }

  return rc;
}

static int
fetch_sdf_tag_data(const Molecule& m, const IWString& tag, float& v)
{
  int n = m.number_records_text_info();

  if (0 == n) {
    cerr << "fetch_sdf_tag_data:no text info in '" << m.name() << "'\n";
    return 0;
  }

  for (int i = 0; i < n; i++) {
    const IWString& info = m.text_info(i);

    if (!info.starts_with(tag)) {
      continue;
    }

    assert(i != n - 1);

    const IWString& next_record = m.text_info(i + 1);

    if (!next_record.numeric_value(v)) {
      cerr << "fetch_sdf_tag_data:invalid value for '" << tag << "' in '" << m.name()
           << "', value '" << next_record << "'\n";
      return 0;
    }

    return 1;
  }

  cerr << "Could not find tag '" << tag << "' in '" << m.name() << "'\n";
  v = static_cast<float>(0.0);
  return 0;
}

static float
fetch_contents_of_column(const IWString& mname, int col)
{
  const_IWSubstring token;

  if (!mname.word(col, token)) {
    cerr << "Cannot extract column " << (col + 1) << " from '" << mname << "'\n";
    return -99.0;
  }

  float rc;
  if (token.numeric_value(rc)) {
    return rc;
  }

  cerr << "Invalid numeric in column " << (col + 1) << " of '" << mname << "'\n";
  return -99.9;
}

static int
compute_atoms_in_largest_ring_system(Molecule& m)
{
  int matoms = m.natoms();

  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  if (1 == nr) {
    return m.ringi(0)->number_elements();
  }

  int* ring_sys = new int[matoms];
  std::unique_ptr<int[]> free_ring_sys(ring_sys);

  m.label_atoms_by_ring_system(ring_sys);

  int* atoms_in_ring_system = new_int(nr + 1);
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
compute_largest_ring_system_size(Molecule& m)
{
  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  int rc = 1;

  int* ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    if (!ri->is_fused()) {
      continue;
    }

    int system_size = 1;

    for (int j = i + 1; j < nr; j++) {
      if (ring_already_done[j]) {
        continue;
      }

      const Ring* rj = m.ringi(j);

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
compute_largest_ring_size(Molecule& m)
{
  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  int rc = m.ringi(0)->number_elements();

  for (int i = 1; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    int ring_size = ri->number_elements();

    if (ring_size > rc) {
      rc = ring_size;
    }
  }

  return rc;
}

static int
smallest_ring_size(Molecule& m)
{
  if (0 == m.nrings()) {
    return 0;
  }

  return m.ringi(0)->number_elements();
}

static int
compute_aromatic_rings(Molecule& m)
{
  int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int rc = 0;

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    if (ri->is_aromatic()) {
      rc++;
    }
  }

  return rc;
}

static int
perform_substructure_search(Molecule& m, int property_number)
{
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
count_usaturation(Molecule& m)
{
  m.compute_aromaticity_if_needed();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (m.is_aromatic(i)) {
      continue;
    }

    const Atom* a = m.atomi(i);

    if (a->ncon() < a->nbonds()) {
      rc++;
    }
  }

  return rc;
}

int
File_Record::set_molecule(Molecule* m)
{
  _smiles = m->smiles();
  _id = m->name();

  return 1;
}

int
File_Record::set_molecule(const const_IWSubstring& buffer)
{
  int i = 0;

  buffer.nextword(_smiles, i);
  if (i == buffer.length()) {
    cerr << "File_Record::set_molecule:no id '" << buffer << "'\n";
    return 0;
  }

  buffer.from_to(i + 1, buffer.length() - 1, _id);

  _id.remove_leading_chars(' ');

  // cerr << "smiles '" << _smiles << "' id '" << _id << "'\n";
  if (_smiles.empty() || _id.empty()) {
    cerr << "File_Record:set_molecule:smiles input must have two tokens\n";
    return 0;
  }

  return 1;
}

static int
do_quick_set_of_properties(Molecule& m, float* property)
{
  int matoms = m.natoms();

  property[0] = static_cast<float>(matoms);
  property[1] = static_cast<float>(m.nrings());

  int carbon = 0;
  int nitrogen = 0;
  int oxygen = 0;
  int fluorine = 0;
  int other_element = 0;
  int hcount = 0;
  int terminal = 0;
  int three_connected = 0;

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    atomic_number_t z = a->element()->atomic_number();

    if (6 == z) {
      carbon++;
    } else if (7 == z) {
      nitrogen++;
    } else if (8 == z) {
      oxygen++;
    } else if (9 == z) {
      fluorine++;
    } else {
      other_element++;
    }

    hcount += m.hcount(i);

    if (3 == a->ncon()) {
      three_connected++;
    } else if (1 == a->ncon()) {
      terminal++;
    }
  }

  property[2] = static_cast<float>(carbon);
  property[3] = static_cast<float>(nitrogen);
  property[4] = static_cast<float>(oxygen);
  property[5] = static_cast<float>(fluorine);
  property[6] = static_cast<float>(other_element);
  property[7] = static_cast<float>(hcount);
  property[8] = static_cast<float>(three_connected);
  property[9] = static_cast<float>(terminal);

  return 1;
}

static int
file_scope_compute_properties(Molecule& m, const IWString& id, float* property)
{
  if (quick_set_of_properties) {
    return do_quick_set_of_properties(m, property);
  }

  for (int i = 0; i < nproperties; i++) {
    if (CMP_NATOMS == comparison_criterion[i]) {
      property[i] = static_cast<float>(m.natoms());
    } else if (CMP_NRINGS == comparison_criterion[i]) {
      property[i] = static_cast<float>(m.nrings());
    } else if (CMP_AMW == comparison_criterion[i]) {
      property[i] = static_cast<float>(compute_molecular_weight(m));
    } else if (CMP_AMW_NOH == comparison_criterion[i]) {
      property[i] = compute_amw_no_h(m);
    } else if (CMP_NFRAG == comparison_criterion[i]) {
      property[i] = static_cast<float>(m.number_fragments());
    } else if (CMP_NCHIRAL == comparison_criterion[i]) {
      property[i] = static_cast<float>(m.chiral_centres());
    } else if (CMP_NHETERO == comparison_criterion[i]) {
      property[i] = static_cast<float>(m.natoms() - m.natoms(6));
    } else if (CMP_AROMATIC_ATOMS == comparison_criterion[i]) {
      property[i] = static_cast<float>(compute_aromatic_atoms(m));
    } else if (CMP_AROMATIC_RINGS == comparison_criterion[i]) {
      property[i] = static_cast<float>(compute_aromatic_rings(m));
    } else if (CMP_LARGEST_RING_SIZE == comparison_criterion[i]) {
      property[i] = static_cast<float>(compute_largest_ring_size(m));
    } else if (CMP_LARGEST_RING_SYSTEM_SIZE == comparison_criterion[i]) {
      property[i] = static_cast<float>(compute_largest_ring_system_size(m));
    } else if (CMP_COLUMN_OF_NAME == comparison_criterion[i]) {
      property[i] = static_cast<float>(fetch_contents_of_column(id, comparison_column[i]));
    } else if (CMP_ROTATBLE_BONDS == comparison_criterion[i]) {
      property[i] = compute_rotatable_bonds(m);
    } else if (CMP_ATOMS_IN_LARGEST_RING_SYSTEM == comparison_criterion[i]) {
      property[i] = compute_atoms_in_largest_ring_system(m);
    } else if (CMP_ATOMS_IN_COUNTERION_MP == comparison_criterion[i]) {
      property[i] = common_compute_atoms_in_counterions(
          m, presumed_atoms_in_counterion_if_no_counterion);
    } else if (CMP_AMW_IN_COUNTERION_MP == comparison_criterion[i]) {
      property[i] = common_compute_amw_in_counterions(
          m, presumed_amw_in_counterion_if_no_counterion);
    } else if (CMP_ATOMS_IN_LARGEST_FRAGMENT == comparison_criterion[i]) {
      property[i] = compute_atoms_in_largest_fragment(m);
    } else if (CMP_ATOMS_IN_FIRST_FRAGMENT == comparison_criterion[i]) {
      property[i] = compute_atoms_in_first_fragment(m);
    } else if (CMP_AMW_IN_LARGEST_FRAGMENT == comparison_criterion[i]) {
      property[i] = compute_amw_in_largest_fragment(m);
    } else if (CMP_NUMBER_ISOTOPIC_ATOMS == comparison_criterion[i]) {
      property[i] = m.number_isotopic_atoms();
    } else if (CMP_SUBSTRUCTURE_SEARCH == comparison_criterion[i]) {
      property[i] = perform_substructure_search(m, i);
    } else if (CMP_UNSATURATION == comparison_criterion[i]) {
      property[i] = count_usaturation(m);
    } else if (CMP_SDF_TAG == comparison_criterion[i]) {
      if (!fetch_sdf_tag_data(m, comparison_tag[i], property[i])) {
        return 0;
      }
    } else if (CMP_RING_ATOMS == comparison_criterion[i]) {
      property[i] = compute_ring_atoms(m);
    } else if (CMP_ATOMIC_NUMBER_TOTAL == comparison_criterion[i]) {
      property[i] = compute_total_atomic_number(m);
    } else if (CMP_ANY_CHARGE == comparison_criterion[i]) {
      property[i] = m.number_formally_charged_atoms();
    } else if (CMP_SMALLEST_RING_SIZE == comparison_criterion[i]) {
      property[i] = smallest_ring_size(m);
    } else {
      cerr << "Unrecognised property " << i << " " << comparison_criterion[i] << "\n";
      return 0;
    }

    if (-1 == direction[i]) {
      property[i] = static_cast<float>(0.0) - property[i];
    }
  }

// #define ECHO_COMPUTED_PROPERTIES
#ifdef ECHO_COMPUTED_PROPERTIES
  for (int i = 0; i < nproperties; i++) {
    cerr << " i = " << i << " property " << property[i] << endl;
  }
#endif

  return 1;
}

int
File_Record::compute_properties()
{
  assert(nproperties > 0);
  if (nullptr == _property) {
    _property = new float[nproperties];
  }

  Molecule m;

  m.build_from_smiles(_smiles);

  return file_scope_compute_properties(m, _id, _property);
}

int
File_Record::compare(const File_Record& rhs) const
{
  for (int i = 0; i < nproperties; i++) {
    //  cerr << " i = " << i << " comparing " << _property[i] << " with " <<
    //  rhs._property[i] << endl;
    if (_property[i] < rhs._property[i]) {
      return 1;
    } else if (_property[i] > rhs._property[i]) {
      return -1;
    }
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
File_Record::echo(IWString& output) const
{
  output << _smiles << ' ' << _id << '\n';

  return 1;
}

// #ifdef __GNUG__
template class resizable_array_p<File_Record>;
template class resizable_array_base<File_Record*>;

// #endif

int
file_record_comparitor_ascending(File_Record* const* cv1, File_Record* const* cv2)
{
  const File_Record* pfr1 = *cv1;
  const File_Record* pfr2 = *cv2;

  return pfr2->compare(*pfr1);
}

int
file_record_comparitor_descending(File_Record* const* cv1, File_Record* const* cv2)
{
  const File_Record* pfr1 = *cv1;
  const File_Record* pfr2 = *cv2;

  return pfr1->compare(*pfr2);
}

class File_Record_Comparitor
{
 private:
  const int _ascending;

 public:
  File_Record_Comparitor(int a) : _ascending(a)
  {
  }

  bool
  operator()(const File_Record*, const File_Record*) const;
};

bool
File_Record_Comparitor::operator()(const File_Record* fr1, const File_Record* fr2) const
{
  const float* f1 = fr1->zdata();
  const float* f2 = fr2->zdata();

  if (_ascending) {
    for (int i = 0; i < nproperties; i++) {
      if (f1[i] > f2[i]) {
        return true;
      }
    }
  } else {
    for (int i = 0; i < nproperties; i++) {
      if (f1[i] < f2[i]) {
        return true;
      }
    }
  }

  return false;
}

class Child_Data
{
 private:
  IWString _smiles_id;
  const float* _zdata;

 public:
  void
  set_smiles_id(const const_IWSubstring& s)
  {
    _smiles_id = s;
  }

  void
  set_data(const float* s)
  {
    _zdata = s;
  }

  const float*
  zdata() const
  {
    return _zdata;
  }

  const IWString&
  smiles_id() const
  {
    return _smiles_id;
  }
};

class Child_Data_Descending
{
 private:
 public:
  int
  operator()(const Child_Data*, const Child_Data*) const;
};

int
Child_Data_Descending::operator()(const Child_Data* cd1, const Child_Data* cd2) const
{
  const float* f1 = cd1->zdata();
  const float* f2 = cd2->zdata();

  for (int i = 0; i < nproperties; i++) {
    if (f1[i] < f2[i]) {
      return 1;
    }

    if (f1[i] > f2[i]) {
      return -1;
    }
  }

  return 0;
}

class Child_Data_Ascending
{
 private:
 public:
  int
  operator()(const Child_Data*, const Child_Data*) const;
};

int
Child_Data_Ascending::operator()(const Child_Data* cd1, const Child_Data* cd2) const
{
  const float* f1 = cd1->zdata();
  const float* f2 = cd2->zdata();

  for (int i = 0; i < nproperties; i++) {
    if (f1[i] < f2[i]) {
      return -1;
    }

    if (f1[i] > f2[i]) {
      return 1;
    }
  }

  return 0;
}

class File_Record_Comparator
{
 private:
  const int _ascending_order;

 public:
  File_Record_Comparator(int a) : _ascending_order(a){};

  bool
  operator()(const File_Record*, const File_Record*) const;
};

bool
File_Record_Comparator::operator()(const File_Record* cd1, const File_Record* cd2) const
{
  const float* f1 = cd1->zdata();
  const float* f2 = cd2->zdata();

  if (_ascending_order) {
    for (int i = 0; i < nproperties; i++) {
      if (f1[i] < f2[i]) {
        return true;
      } else if (f1[i] > f2[i]) {
        return false;
      }
    }
  } else {
    for (int i = 0; i < nproperties; i++) {
      if (f1[i] > f2[i]) {
        return true;
      }
    }
  }

  return false;
}

class Child_Data_Comparator
{
 private:
  const int _ascending_order;

 public:
  Child_Data_Comparator(int a) : _ascending_order(a){};

  bool
  operator()(const Child_Data*, const Child_Data*) const;
};

bool
Child_Data_Comparator::operator()(const Child_Data* cd1, const Child_Data* cd2) const
{
  const float* f1 = cd1->zdata();
  const float* f2 = cd2->zdata();

  if (_ascending_order) {
    for (int i = 0; i < nproperties; i++) {
      if (f1[i] < f2[i]) {
        return true;
      } else if (f1[i] > f2[i]) {
        return false;
      }
    }
  } else {
    for (int i = 0; i < nproperties; i++) {
      if (f1[i] > f2[i]) {
        return true;
      }
    }
  }

  return false;
}

// Group all values that are the same into the same file.
int
identify_next_group_same_value_same_file(const resizable_array_p<File_Record> &records,
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
identify_next_group(const resizable_array_p<File_Record>& records, int& istart,
                    int& istop)
{
  int n = records.number_elements();

  // cerr << "identify_next_group:on entry istart " << istart << " istop " << istop << ",
  // min_size_hint " << min_size_hint << endl;

  if (istop >= n) {
    return 0;
  }

  istart = istop;

  const File_Record* reff = records[istart];

  int last_break = -1;

  istop = istart + 1;

  while (istop < n) {
    // same, just continue accumulating
    if (reff->Same(*(records[istop]))) {
      istop++;
      continue;
    }

    //  We are at a discontinuity

    if (istop - istart < min_size_hint) {  // not enough to write
      ;
    } else if (last_break <
               0) {  // in the first chunk of molecules, we'll check the next group too
      ;
    } else if (istop - istart >
               (min_size_hint + min_size_hint))  // too much, write previous chunk
    {
      istop = last_break;
      return 1;
    } else if (0 ==
               istart)  // otherwise we get problems when just two chunks migh be written
    {
      istop++;
      return 1;
    }

    reff = records[istop];
    last_break = istop;
    istop++;
  }

  return 1;
}

template <typename T>
static int
echo_set_of_molecules(const resizable_array_p<T>& records, int istart, int istop,
                      IWString_and_File_Descriptor& output)
{
  assert(istart <= istop);

  for (int i = istart; i < istop; i++) {
    if (!records[i]->echo(output)) {
      cerr << "Cannot write item " << i << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

template <typename T>
int
write_group(const resizable_array_p<T>& molecules, int istart, int istop,
            const IWString& stem, int ndx)
{
  IWString fname;

  fname << stem << ndx << ".smi";

  IWString_and_File_Descriptor output;
  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open group file '" << fname << "'\n";
    return 0;
  }

  return echo_set_of_molecules(molecules, istart, istop, output);
}

/*
  Generally, we are sorting by atom count, and we want to make sure the
  chunk with the largest molecules is not the largest chunk, so identify
  the chunks in descending order
*/

template <typename T>
int
identify_next_equi_atom_group_descending(const resizable_array_p<T>& records, int& istart,
                                         int& istop, const int atoms_per_chunk,
                                         unsigned int& atoms)
{
  // cerr << "identify_next_group:on entry istart " << istart << " istop " << istop << ",
  // min_size_hint " << min_size_hint << endl;

  if (istop <= 0) {
    return 0;
  }

  istart = istop;

  const T* reff = records[istart];

  istop = istart - 1;

  int atoms_this_chunk = 0;

  while (istop >= 0) {
    atoms_this_chunk += records[istop]->natoms();  // + 20;

    // same, just continue accumulating
    if (reff->Same(*(records[istop]))) {
      istop--;
      continue;
    }

    //  We are at a discontinuity

    if (atoms_this_chunk >= atoms_per_chunk) {
      //    cerr << "atoms_this_chunk " << atoms_this_chunk << " per chunk " <<
      //    atoms_per_chunk << endl;

      atoms -= atoms_this_chunk;
      if (atoms < (static_cast<unsigned int>(atoms_per_chunk) * 90 / 100)) {
        istop = 0;
      } else {
        istop++;
      }
      return 1;
    }

    reff = records[istop];
    istop--;
  }

  istop = 0;

  return 1;
}

template <typename T>
int
identify_next_equi_atom_group(const resizable_array_p<T>& records, int& istart,
                              int& istop, int atoms_per_chunk, unsigned int& atoms)
{
  int n = records.number_elements();

  // cerr << "identify_next_group:on entry istart " << istart << " istop " << istop << ",
  // min_size_hint " << min_size_hint << endl;

  if (istop >= n) {
    return 0;
  }

  istart = istop;

  const T* reff = records[istart];

  int last_break = -1;

  istop = istart + 1;

  int atoms_this_chunk = 0;

  while (istop < n) {
    atoms_this_chunk += records[istop]->natoms();

    if (0 == reff->compare(*(records[istop])))  // same, just continue accumulating
    {
      istop++;
      continue;
    }

    //  We are at a discontinuity

    if (atoms_this_chunk >= atoms_per_chunk) {
      atoms -= atoms_this_chunk;
      if (atoms < atoms_per_chunk) {
        istop = n - 1;
      }
      return 1;
    }

    reff = records[istop];
    last_break = istop;
    istop++;
  }

  return 1;
}
template <typename T>
int
do_same_score_chunked_output(const resizable_array_p<T>& molecules,
                  const IWString& stem_for_output_files)
{
  int ndx = 1;

  int istart = 0;
  int istop = 0;
  while (identify_next_group_same_value_same_file(molecules, istart, istop)) {
    //  cerr << "Writing from " << istart << " to " << istop << endl;

    if (file_name_index_is_first_property) {
      ndx = static_cast<int>(molecules[istart]->FirstProperty() + 0.001);
    }
    if (!write_group(molecules, istart, istop, stem_for_output_files, ndx)) {
      return 0;
    }

    ndx++;
  }

  if (verbose) {
    cerr << "Created " << ndx << " chunked output files\n";
  }

  return 1;
}


template <typename T>
int
do_equal_atom_count_chunked_output(const resizable_array_p<T>& molecules,
                                   const unsigned int total_atom_count,
                                   const IWString& stem_for_chunked_output)
{
  int atoms_per_chunk = total_atom_count / splits_to_contain_equal_total_atoms;

  // cerr << "total_atom_count " << total_atom_count << endl;

  if (verbose) {
    cerr << "Chunks will contain roughly " << atoms_per_chunk << " atoms\n";
  }

  int ndx = 1;

  unsigned int atoms = total_atom_count;

  int istart = 0;
  int istop = 0;

  istart = molecules.number_elements() - 1;
  istop = istart;

  int apc = atoms_per_chunk * scale_first_chunk /
            100;  // heuristic to put fewer atoms in first chunk

  while (identify_next_equi_atom_group_descending(molecules, istart, istop, apc, atoms)) {
    if (verbose) {
      cerr << "Writing from " << istart << " to " << istop << "\n";
    }

    if (!write_group(molecules, istop, istart, stem_for_chunked_output, ndx)) {
      return 0;
    }

    ndx++;

    apc = atoms_per_chunk;
  }

  if (verbose) {
    cerr << "Wrote " << (ndx - 1) << " eqi atom count groups, stem '"
         << stem_for_chunked_output << "'\n";
  }

  return 1;
}

template <typename T>
int
do_chunked_output(const resizable_array_p<T>& molecules,
                  const IWString& stem_for_output_files)
{
  int ndx = 1;

  int istart = 0;
  int istop = 0;
  while (identify_next_group(molecules, istart, istop)) {
    //  cerr << "Writing from " << istart << " to " << istop << endl;

    if (!write_group(molecules, istart, istop, stem_for_output_files, ndx)) {
      return 0;
    }

    ndx++;
  }

  if (verbose) {
    cerr << "Created " << ndx << " chunked output files\n";
  }

  return 1;
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

#ifdef NOT_USED_SUBSTRUCTURE_SEARCH_NOT_THREAD_SAFE
static int
fetch_queries(const_IWSubstring& q)
{
  if (!process_cmdline_token(' ', q, queries, verbose)) {
    cerr << "Cannot process query specification '" << q << "'\n";
    return 0;
  }

  return queries.number_elements();
}

static int
fetch_smarts(const const_IWSubstring& s)
{
  Substructure_Query* q = new Substructure_Query;
  if (!q->create_from_smarts(s)) {
    cerr << "Invalid smarts '" << s << "'\n";
    return 0;
  }

  queries.add(q);

  return queries.number_elements();
}
#endif  // NOT_USED_SUBSTRUCTURE_SEARCH_NOT_THREAD_SAFE

static int
read_molecules_text(iwstring_data_source& input, resizable_array_p<File_Record>& records)
{
  Report_Progress report_progress;

  if (file_scope_report_progress.active()) {
    report_progress.initialise(file_scope_report_progress);
  }
  input.set_translate_tabs(1);

  int rc = 0;

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    File_Record* f = new File_Record;

    if (! f->set_molecule(buffer)) {
      cerr << "fatal error reading '" << buffer << "', line " << input.lines_read() << '\n';
      return 0;
    }

    records.add(f);

    if (report_progress()) {
      cerr << "Processed " << records.number_elements() << " molecules\n";
    }

    rc++;
  }

  return rc;
}

static int
read_child_data(int fd, float* zdata, off_t bytes_to_read)
{
  ssize_t rc = IW_FD_READ(fd, zdata, bytes_to_read);

  if (rc == bytes_to_read) {
    return 1;
  }

  cerr << "Gack, cannot read " << bytes_to_read << " bytes of child data rc = " << rc
       << " errbo " << errno << endl;
  perror("read");
  return 0;
}

static int
read_child_data(const char* fname, float* zdata, int& ndx)
{
  int fd = IW_FD_OPEN(fname, O_RDONLY);

  if (fd < 0) {
    cerr << "Gack, cannot open '" << fname << "'\n";
    return 0;
  }

  off_t file_size = dash_s(fname);

  if (verbose > 1) {
    cerr << "Reading file of size " << file_size << " bytes\n";
  }

  if (!read_child_data(fd, zdata + ndx, file_size)) {
    return 0;
  }

  // Increment ndx to account for what we just read

  ndx += file_size / sizeof(float);  // number data points we just read

  // cerr << "ndx incremented to " << ndx << endl;

  return 1;
}

#define WRITE_BLOCK_CHUNK_SIZE 8192

static int
write_some_data(const float* zdata, int nmolecules, int fd)
{
  int bytes_to_write = nmolecules * nproperties * sizeof(float);

  int rc = IW_FD_WRITE(fd, zdata, bytes_to_write);

  if (rc == bytes_to_write) {
    return 1;
  }

  if (rc < 0) {
    cerr << "Cannot write " << bytes_to_write << " bytes\n";
    perror("write");
    cerr << endl;
  } else {
    cerr << "Partial write, wrote " << rc << " of " << bytes_to_write << " bytes\n";
  }

  return 0;
}

static int
do_forked_child_process(data_source_and_type<Molecule>& input, float* zdata, int fd)
{
  int molecules_read = 0;
  int ndx = 0;

  int write_every = WRITE_BLOCK_CHUNK_SIZE / nproperties;

  if (verbose > 1) {
    cerr << "Will write every " << write_every << " molecules\n";
    cerr << "Input contains " << input.molecules_remaining() << " molecules\n";
  }

  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    file_scope_compute_properties(*m, m->name(), zdata + (ndx * nproperties));

    ndx++;

    if (ndx != write_every) {  // no writing this time
      ;
    } else if (write_some_data(zdata, ndx, fd)) {
      ndx = 0;
    } else {
      cerr << "Returning for write failure\n";
      return 0;
    }
  }

  if (0 == ndx) {
    ;
  } else if (write_some_data(zdata, ndx, fd)) {
    ;
  } else {
    return 0;
  }

  IW_FD_CLOSE(fd);

  return 1;
}

static int
do_forked_child_process(const char* fname, FileType input_type, float* zdata, int fd)
{
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  return do_forked_child_process(input, zdata, fd);
}

static void
do_forked_child_process(const char* input_fname, FileType input_type, int ndx,
                        const IWString& forked_process_stem)
{
#ifdef THREADS_OPEN_LOGFILES
  IWString logfilename;
  logfilename << "MSORT" << ndx << ".log";

  IWString_and_File_Descriptor logfile;

  if (!logfile.open(logfilename.null_terminated_chars())) {
    cerr << "Cannot open per thread logfile '" << logfilename << "'\n";
    exit(3);
  }

  logfile << "Child process doing " << nproperties << " properties\n";
  for (int i = 0; i < nproperties; i++) {
    logfile << " i = " << i << " property " << comparison_criterion[i] << '\n';
  }

  logfile.flush();
#endif

  IWString output_fname;
  output_fname << forked_process_stem << ndx << ".dat";

  int mode = O_WRONLY | O_TRUNC | O_CREAT;

#ifdef _WIN32
  int flags = 0;
#else
  int flags = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
#endif

  int fd = IW_FD_OPEN(output_fname, mode, flags);

  if (fd < 0) {
    cerr << "do_forked_child_process:open '" << output_fname << "'\n";
    exit(1);
  }

  int molecules_per_chunk = WRITE_BLOCK_CHUNK_SIZE / nproperties + 1;

  float* zdata = new float[(molecules_per_chunk + 1) * nproperties];
  std::unique_ptr<float[]> free_zdata(zdata);

  if (!do_forked_child_process(input_fname, input_type, zdata, fd)) {
    exit(2);
  }

  // cerr << "Child processing complete, sleeping\n";
  // sleep(ndx + 1);

  exit(0);
}

static int
read_smiles_id(iwstring_data_source& input, resizable_array_p<Child_Data>& cd)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    Child_Data* tmp = new Child_Data;
    tmp->set_smiles_id(buffer);
    cd.add(tmp);
  }

  return 1;
}

static int
read_smiles_id(const char* fname, resizable_array_p<Child_Data>& cd)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles_id(input, cd);
}

static int
do_via_fork(Command_Line& cl, FileType input_type)
{
  int nfiles = cl.number_elements();

  for (int i = 0; i < nfiles; i++) {
    const_IWSubstring fname(cl[i]);

    if (!fname.ends_with(".smi")) {
      cerr << "Sorry, forked processing only works with smiles files '" << fname << "'\n";
      return 0;
    }

    if (!dash_s(cl[i])) {
      cerr << "Missing or empty input file '" << fname << "'\n";
      return 0;
    }
  }

  resizable_array<pid_t> children;
  children.resize(nfiles);

  if (verbose) {
    cerr << "Launching " << nfiles << " children\n";
  }

  for (int i = 0; i < nfiles; i++) {
    pid_t p = fork();
    if (0 == p) {
      do_forked_child_process(cl[i], input_type, i, forked_process_stem);
    } else {
      children.add(p);
      cerr << "kill " << p << endl;
    }
  }

  if (children.number_elements() != nfiles) {
    cerr << "Yipes, only launched " << children.number_elements() << " child processes\n";
    return 0;
  }

  // While the children are off doing work, we assemble the smiles and ID of
  // the participants

  resizable_array_p<Child_Data> cd;
  cd.resize(initial_size_estimate);

  for (int i = 0; i < nfiles; i++) {
    if (!read_smiles_id(cl[i], cd)) {
      cerr << "Cannot read smiles and id data from '" << cl[i] << "'\n";
      return 0;
    }
  }

  if (verbose) {
    cerr << "Parent read data on " << cd.number_elements() << " items from " << nfiles
         << " files\n";
  }

  int fatal_error = 0;
  int molecules_computed = 0;

  for (int i = 0; i < nfiles; i++) {
    IWString fname;
    fname << forked_process_stem << i << ".dat";

    waitpid(children[i], NULL, 0);

    int rc = kill(children[i], 0);
    if (0 == rc) {
      cerr << "Bad news, child " << i << " might still be alive!!\n";
    }

    off_t file_size = dash_s(fname.null_terminated_chars());

    if (0 == file_size) {
      cerr << "Child created file '" << fname << "' missing or empty\n";
      fatal_error = 1;
    }

    if (0 != (file_size % (nproperties * sizeof(float)))) {
      cerr << "Child file '" << fname << "' invalid size " << file_size << endl;
      fatal_error = 1;
    }

    molecules_computed += file_size / (nproperties * sizeof(float));
  }

  if (fatal_error) {
    return 0;
  }

  if (verbose) {
    cerr << "Child processes wrote data for " << molecules_computed << " molecules\n";
  }

  if (0 == molecules_computed) {
    cerr << "No data produced by child processes\n";
    return 0;
  }

  if (molecules_computed != cd.number_elements()) {
    cerr << "Yipes, size mismatch between parent and children\n";
    cerr << "Parent read " << cd.number_elements() << " smiles, children computed "
         << molecules_computed << endl;
    return 0;
  }

  int floats_needed = molecules_computed * (nproperties);

  if (verbose > 1) {
    cerr << "Final sort requires " << floats_needed << " floating point numbers\n";
  }

  float* zdata = new float[floats_needed];
  std::unique_ptr<float[]> free_zdata(zdata);

  if (nullptr == zdata) {
    cerr << "Yipes, cannot allocate " << floats_needed << " floating point numbers\n";
    return 0;
  }

  int ndx = 0;
  for (int i = 0; i < nfiles; i++) {
    IWString fname;
    fname << forked_process_stem << i << ".dat";

    if (!read_child_data(fname.null_terminated_chars(), zdata, ndx)) {
      cerr << "Cannot read child data from '" << fname << "'\n";
      return 0;
    }

    if (!keep_temporary_files) {
      unlink(fname.null_terminated_chars());
    }
  }

  // cerr << "After reading child data ndx " << ndx << endl;

  for (int i = 0; i < molecules_computed; i++) {
    cd[i]->set_data(zdata + i * nproperties);
  }

  time_t tstart = time(NULL);
  if (verbose) {
    cerr << "Begin final sort of " << cd.number_elements() << " molecules\n";
  }

  if (do_parallel_sort) {

    Child_Data** zdata = cd.rawdata();

    Child_Data_Comparator cdc(ascending_order);
    parallel_sort(zdata, zdata + cd.number_elements(), cdc);
  } else if (ascending_order) {
    Child_Data_Ascending cda;
    cd.iwqsort(cda);
  } else {
    Child_Data_Descending cdd;
    cd.iwqsort(cdd);
  }

  if (verbose) {
    time_t tend = time(NULL);

    cerr << "Sorting took " << (tend - tstart) << " seconds\n";
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < molecules_computed; i++) {
    output << cd[i]->smiles_id() << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

struct Thread_Data {
  int _thread_number;
  File_Record** fstart;
  File_Record** fstop;
  int computations_done;
  time_t elapsed_time;
};

static void*
thread_task(void* v)
{
  Thread_Data* td = reinterpret_cast<Thread_Data*>(v);

  time_t tzero = time(NULL);

  int rc = 0;

  if (verbose) {
    cerr << "Thread computing " << (td->fstop - td->fstart + 1) << " items\n";
  }

  for (File_Record** i = td->fstart; i <= td->fstop; i++) {
    File_Record* f = (*i);

    f->compute_properties();

    rc++;
  }

  td->computations_done = rc;

  time_t tend = time(NULL);

  td->elapsed_time = (tend - tzero);

  pthread_exit(NULL);

  return NULL;
}

static int
process_via_pthreads(File_Record** pool, int nmolecules, int nthreads)
{
  Thread_Data* td = new Thread_Data[nthreads];
  std::unique_ptr<Thread_Data[]> free_td(td);

  pthread_t* child_thread = new pthread_t[nthreads];
  std::unique_ptr<pthread_t[]> free_child_thread(child_thread);

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  int items_per_thread = nmolecules / nthreads;

  assert(items_per_thread > 0);

  if (verbose) {
    cerr << nmolecules << " molecules across " << nthreads << " means "
         << items_per_thread << " items per thread\n";
  }

  for (int i = 0; i < nthreads; i++) {
    struct Thread_Data& tdi = td[i];
    tdi._thread_number = i;
    tdi.computations_done = 0;

    tdi.fstart = pool + i * items_per_thread;

    if (i == nthreads - 1) {
      tdi.fstop = pool + nmolecules - 1;
    } else {
      tdi.fstop = tdi.fstart + items_per_thread - 1;
    }

    if (0 != pthread_create(&(child_thread[i]), &attr, thread_task, &tdi)) {
      perror("cannot create child process");
      return 0;
    }
  }

  Accumulator_Int<int> computation;

  int computations_done = 0;

  for (int i = 0; i < nthreads; i++) {
    pthread_join(child_thread[i], NULL);
    computations_done += td[i].computations_done;

    //  cerr << "Joined thread " << i << " computations_done " << td[i].computations_done
    //  <<endl; cerr << "Time " << td[i].elapsed_time << endl;

    if (verbose) {
      cerr << "Thread " << i << " performed " << td[i].computations_done << "\n";
    }

    if (td[i].elapsed_time > 0) {
      int computations_per_second = td[i].computations_done / td[i].elapsed_time;

      computation.extra(computations_per_second);
    }
  }

  if (verbose) {
    cerr << "Performed " << computations_done << " computations\n";
    if (computation.n() > 1) {
      cerr << "Rates between " << computation.minval() << " and " << computation.maxval()
           << " ave " << computation.average() << endl;
    }
  }

  return 1;
}

static int
read_molecules_text(const char* fname, resizable_array_p<File_Record>& records)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_molecules_text(input, records);
}

static int
read_molecules(data_source_and_type<Molecule>& input,
               resizable_array_p<File_Record>& records)
{
  int rc = 0;

  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    File_Record* f = new File_Record;

    if (!f->set_molecule(m)) {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 4;
    }

    records.add(f);

    rc++;

    if (0 == verbose) {
      ;
    } else if (0 == records.number_elements() % 10000) {
      cerr << "Read " << records.number_elements() << " molecules\n";
    }
  }

  return rc;
}

static int
read_molecules(const char* fname, FileType input_type,
               resizable_array_p<File_Record>& records)
{
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return read_molecules(input, records);
}

class TBB_Compute_Properties
{
 private:
 public:
  void
  operator()(const tbb::blocked_range<File_Record**>& r) const;
};

void
TBB_Compute_Properties::operator()(const tbb::blocked_range<File_Record**>& r) const
{
  for (File_Record** i = r.begin(); i != r.end(); i++) {
    File_Record* f = *i;

    f->compute_properties();
  }

  return;
}

static int
display_dash_m_options(std::ostream& os)
{
  os << " -M ea=<n>         write <n> chunks of output each containing roughly equal "
        "total atoms\n";

  exit(1);
}

static unsigned int
count_atoms(const resizable_array_p<File_Record>& records, const int molecules_read,
            const int nproperties)
{
  const File_Record* const* fr = records.rawdata();

  unsigned int rc = 0;

  for (int i = 0; i < molecules_read; i++) {
    float f = fr[i]->natoms();

    rc += static_cast<unsigned int>(f);
  }

  return rc;
}

static int
msort(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:mfcarh:onjRSbzZk:dylg:D:e:qs:pM:x:Y:");

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

  if (cl.option_present('x')) {
    if (!file_scope_report_progress.initialise(cl, 'x', verbose)) {
      cerr << "Cannot initialise report progress (-x) option\n";
      usage(1);
    }
  }

  if (cl.option_present('Y')) {
    if (!cl.value('Y', scale_first_chunk) || scale_first_chunk < 1 ||
        scale_first_chunk > 99) {
      cerr << "The percentage scale largest chunk option (-Y) must be a valid "
              "percentage\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Chunk with largest molecules scaled to " << scale_first_chunk
           << " of regular size\n";
    }
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

  if (cl.option_present('d')) {
    ascending_order = 1;

    if (verbose) {
      cerr << "Descending order sort\n";
    }
  }

  int nthreads = 0;

  if (cl.option_present('h')) {
    if (!cl.value('h', nthreads) || nthreads < 1) {
      cerr << "The number of threads must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will use " << nthreads << " pthreads for processing\n";
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', initial_size_estimate) || initial_size_estimate < 2) {
      cerr << "The initial size estimate (-s) option must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Initial size " << initial_size_estimate << endl;
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  comparison_criterion = new int[1000];  // definite overestimate
  direction = new int[1000];

  int number_query_specifications = 0;

  if (cl.option_present('k')) {
    if (nproperties) {
      cerr << "Sorry, the -k option cannot be combined with other options\n";
      return 3;
    }

    int nk = cl.option_count('k');
    //  cerr << "nk " << nk << endl;

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

        //      if (token.starts_with("amw"))
        //      {
        //        cerr << "Sorry, sorting by any variant of amw not possible in parallel
        //        version, use the serial version\n"; return 1;
        //      }

        if ('a' == token || token.starts_with("atom") || "natoms" == token) {
          comparison_criterion[nproperties] = CMP_NATOMS;
        } else if ('r' == token || "rings" == token || token.starts_with("nring")) {
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
        }
#ifdef QUERIES_ARE_NOT_THREAD_SAFE
        else if (token.starts_with("qry=")) {
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
        }
#endif
        else if (token.starts_with("aicm=")) {
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
        } else if (token.starts_with("tag=")) {
          token.remove_leading_chars(4);
          comparison_criterion[nproperties] = CMP_SDF_TAG;
          comparison_tag[nproperties] << ">  <" << token << '>';
        } else if (token == "charge") {
          comparison_criterion[nproperties] = CMP_ANY_CHARGE;
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

        //      cerr << "Just set comparison_criterion " << nproperties << " to " <<
        //      comparison_criterion[nproperties] << endl;
        nproperties++;
      }
    }
  } else {
    quick_set_of_properties = 1;
    nproperties = 10;
  }

  if (0 == nproperties) {
    cerr << "Must specify one or more sort criteria\n";
    usage(3);
  }

#ifdef ECHO_CRITERIA
  for (int i = 0; i < nproperties; i++) {
    cerr << "Sort criterion " << i << " is " << comparison_criterion[i] << endl;
  }
#endif

  if (cl.option_present('e') && !cl.option_present('D')) {
    cerr << "The -e option only makes sense with the -D option\n";
    usage(3);
  }

  if (cl.option_present('D')) {
    cl.value('D', stem_for_chunked_output);

    if (verbose) {
      cerr << "Chunked output written with stem '" << stem_for_chunked_output << "'\n";
    }
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

  if (cl.option_present('p')) {
    do_parallel_sort = 1;

    if (verbose) {
      cerr << "Will do final sort in parallel\n";
    }
  }

  if (cl.option_present('f')) {
    if (1 == cl.number_elements()) {
      cerr << "Forked processing works by allocating a process to each input file\n";
      cerr << "Use iwsplit to split the file\n";
      return 3;
    }

    if (!do_via_fork(cl, input_type)) {
      cerr << "Forked processing failed\n";
      return 3;
    }

    return 0;
  }

  resizable_array_p<File_Record> records;
  records.resize(1000000);

  int molecules_read = 0;

  for (int i = 0; i < cl.number_elements(); i++) {
    const_IWSubstring fname(cl[i]);

    int tmp;
    if (fname.ends_with(".smi")) {
      tmp = read_molecules_text(cl[i], records);
    } else {
      tmp = read_molecules(cl[i], input_type, records);
    }

    if (0 == tmp) {
      cerr << "Cannot read molecules from '" << cl[i] << "'\n";
      return 3;
    }

    molecules_read += tmp;
  }

  if (0 == molecules_read) {
    cerr << "No molecules\n";
    return 3;
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  // Use parallelism to compute the properties

  File_Record** zdata = (File_Record**)records.rawdata();

  if (nthreads) {
    process_via_pthreads(zdata, molecules_read, nthreads);
  } else {
    TBB_Compute_Properties tbcp;

    tbb::parallel_for(tbb::blocked_range<File_Record**>(zdata, zdata + molecules_read),
                      tbcp);
  }

  time_t tzero = static_cast<time_t>(0);

  if (verbose) {
    cerr << "Data loaded, begin sort\n";
    tzero = time(NULL);
  }

// #define DEBUG_SORTING
#ifdef DEBUG_SORTING
  cerr << "Will compare " << nproperties << " properties\n";
  for (auto i = 0; i < records.number_elements(); ++i) {
    cerr << records[i]->id() << ' ' << records[i]->zdata()[0] << endl;
  }
#endif

  if (do_parallel_sort) {

    File_Record** cd = records.rawdata();

    File_Record_Comparator frc(ascending_order);
    parallel_sort(cd, cd + records.number_elements(), frc);

#ifdef CHECK_RESULTS
    cerr << "Sorting complete, checking...\n";
    for (int i = 1; i < records.number_elements(); i++) {
      cerr << " i = " << i << " compare:\n";
      bool result = frc(cd[i - 1], cd[i]);
      cerr << "Result " << result << endl;
    }
#endif
  } else if (ascending_order) {
    records.sort(file_record_comparitor_ascending);
  } else {
    records.sort(file_record_comparitor_descending);
  }

  if (verbose) {
    time_t tend = time(NULL);
    cerr << "Sort complete, sorting took " << (tend - tzero) << " seconds\n";
  }

#ifdef DEBUG_SORTING
  cerr << "After sort\n";
  for (auto i = 0; i < records.number_elements(); ++i) {
    cerr << records[i]->id() << ' ' << zdata[i]->zdata()[0] << endl;
  }
#endif

  if (cl.option_present('D')) {
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
      unsigned int total_atom_count = count_atoms(records, molecules_read, nproperties);
      do_equal_atom_count_chunked_output(records, total_atom_count, stem_for_chunked_output);
    } else if (splits_contain_equal_values) {
      do_same_score_chunked_output(records, stem_for_chunked_output);
    } else {
      do_chunked_output(records, stem_for_chunked_output);
    }
  } else {
    IWString_and_File_Descriptor output(1);

    if (!echo_set_of_molecules(records, 0, molecules_read, output)) {
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
main(int argc, char** argv)
{
  int rc = msort(argc, argv);

  return rc;
}
