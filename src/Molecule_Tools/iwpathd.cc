// This is an interesting tool. Run-times are not good, but these fingerprints
// do sometimes show up as being strong contributors to a model. This needs some
// attention.

#include <fcntl.h>
#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>

#include <sys/stat.h>
#include <sys/types.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
// #define SUBSTRUCTURE_SEARCH_TSHAPED_PATHS
#ifdef SUBSTRUCTURE_SEARCH_TSHAPED_PATHS
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"
#endif

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int min_path_length = 0;

static int max_path_length = 7;

static int build_disubstituted_ring_paths = 1;

static int build_tshape_paths = 0;

static int discard_tshaped_paths_that_form_rings = 1;

static Atom_Typing_Specification atom_typing_specification;

/*
  We need a unique identifier for each path
*/

static IW_STL_Hash_Map_uint path_to_unique_id;

static resizable_array_p<IWString> id_to_path;

/*
  Cumulative statistics on the number of paths of different types
*/

static int path_type_no_ring = 0;
static int path_type_rings_both_ends = 0;
static int path_type_chain_to_rings = 0;
static int path_type_too_complicated = 0;

/*
  Support for the -Y option
*/

static IWString_and_File_Descriptor stream_for_paths;

static int ntest = 0;

static int skipped_testing_non_sssr_ring_molecules = 0;

static int keep_going_after_test_failure = 0;

static int molecules_failing_tests = 0;

static IWString die_when_formed;

static Element_Transformations element_transformations;

/*
  When running tests, it is seldom useful to keep track of the global counts
*/

static int free_global_structures_after_each_molecule = 0;

/*
  For each path, we keep track of the number of occurrences of each path
*/

static extending_resizable_array<int> number_occurrences;

static extending_resizable_array<int> paths_per_molecule;

static IWDigits iwdigits;

static int min_heteroatoms_needed_in_path = 0;

static int include_hydrogen_info_in_invariant = 0;

static int only_include_hydrogen_info_on_heteroatoms = 0;

static int function_as_tdt_filter = 0;

static IWString smiles_tag("$SMI<");

static IWString fingerprint_tag;

static IWString identifier_tag("PCN<");

static IWString empty_constant_width_fingerprint;
static IWString empty_sparse_fingerprint;

static int create_sparse_fingerprint = 0;

static int create_constant_width_fingerprint = 0;

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
  cerr << "Looks for linear paths in molecules\n";
  cerr << "  -m len        only enumerate paths >= len in size\n";
  cerr << "  -M len        only enumerate paths <= len in size\n";
  cerr << "  -c <n>        suppress paths with <n>  or fewer occurrences\n";
  cerr << "  -c <n%>       suppress paths with <n>% or fewer occurrences\n";
  cerr << "  -C <n>        suppress paths with <n>  or more occurrences\n";
  cerr << "  -C <n%>       suppress paths with <n>% or more occurrences\n";
  cerr << "  -p <string>   prefix for all descriptor names\n";
  cerr << "  -Y fname      output atom paths to filename\n";
  cerr << "  -h <num>      only consider paths containing <num> or more heteroatom\n";
  cerr << "  -y 1          indicate presence or absence of hydrogen on each atom\n";
  cerr << "  -y c          include full hydrogen count\n";
  cerr << "  -y h          only include implicit hydrogen info on heteroatoms\n";
  cerr << "  -t help       testing options, enter '-t help' for details\n";
  cerr << "  -s <fname>    character translation table - translate chars in paths\n";
  cerr << "  -X ...        miscellaneous options, enter '-X help' for details\n";
  cerr << "  -J <tag>      generate fingerprint with tag <tag>\n";
  cerr << "  -f            function as a tdt filter\n";
  cerr << "  -T ele1=ele2  standard element transformation options\n";
  cerr << "  -P <type>     specify atom typing\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static unsigned int
assign_unique_identifier(IW_STL_Hash_Map_uint& path_to_unique_id, const IWString& p)
{
  IW_STL_Hash_Map_uint::const_iterator f = path_to_unique_id.find(p);

  if (f != path_to_unique_id.end()) {
    unsigned int uid = (*f).second;
    number_occurrences[uid]++;

    return uid;
  }

  unsigned int s = path_to_unique_id.size();

  path_to_unique_id[p] = s;

  IWString* tmp = new IWString(p);

  id_to_path.add(tmp);

  number_occurrences[s] = 1;

  return s;
}

static int
write_path_cross_reference(const IW_STL_Hash_Map_uint& path_to_unique_id,
                           IWString_and_File_Descriptor& output)
{
  for (IW_STL_Hash_Map_uint::const_iterator i = path_to_unique_id.begin();
       i != path_to_unique_id.end(); ++i) {
    output << (*i).first << ' ' << (*i).second << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
write_path_cross_reference(const IW_STL_Hash_Map_uint& path_to_unique_id,
                           const char* fname)
{
  IWString_and_File_Descriptor output;
  if (!output.open(fname)) {
    cerr << "write_path_cross_reference:cannot open '" << fname << "'\n";
    return 0;
  }

  return write_path_cross_reference(path_to_unique_id, output);
}

class Molecule_and_Paths : public IWString
{
 private:
  resizable_array<int> _paths_this_molecule;

 public:
  Molecule_and_Paths(const IWString& s) : IWString(s)
  {
  }

  resizable_array<int>&
  paths()
  {
    return _paths_this_molecule;
  }

  int
  get_path_counts(int*) const;
};

static resizable_array_p<Molecule_and_Paths> molecule_and_paths;

int
Molecule_and_Paths::get_path_counts(int* tmp) const
{
  int n = _paths_this_molecule.number_elements();

  for (int i = 0; i < n; i++) {
    int p = _paths_this_molecule[i];

    if (0 == number_occurrences[p]) {
      continue;
    }

    tmp[p]++;
  }

  return 1;
}

template class resizable_array_p<Molecule_and_Paths>;
template class resizable_array_base<Molecule_and_Paths*>;

class PD_Molecule_Data
{
 private:
  Molecule* _m;

  int* _valid_start_atom;

  int _matoms;

  IWString* _atomic_symbol;

  char* _bond_symbol;

  int* _in_same_ring;

  atomic_number_t* _atomic_number;

  int* _aromatic;

  int* _symmetry_class;

  //  We can designate atoms as having or not having an attached hydrogen

  int* _hydrogen_attached;

  //  The computed invariant is a function of the _atom_invariant and adjustments
  //  due to the bonding. We save time by fully computing the atom invariant

  int* _atom_invariant;

  const Atom** _atom;

  //  When adding =O groups, it is handy to have that precomputed

  Set_of_Atoms* _doubly_bonded_singly_connected_neighbour;

  //  And we need to know which atoms are automatically added

  int* _automatically_added;

  int* _are_bonded;

  //  private functions

  int
  _fill_bond_symbols(Molecule& m);
  int
  _identify_automatically_added(const Molecule& m);
  int
  _looks_like_part_of_cf3(const Molecule&, atom_number_t) const;
  int
  _identify_bonded_atoms(const Molecule& m);
  int
  _identify_doubly_bonded_singly_connected_neighbours(Molecule& m);
  int
  _identify_valid_start_atoms(Molecule& m);
  int
  _fill_atom_symbol_array(Molecule& m);
  int
  _fill_atom_symbol_array(const int*);
  int
  _assign_atomic_invariants(Molecule& m);

 public:
  PD_Molecule_Data(Molecule&);
  ~PD_Molecule_Data();

  const Molecule*
  molecule() const
  {
    return _m;
  }

  const Atom*
  atomi(atom_number_t i) const
  {
    return _atom[i];
  }

  const atomic_number_t*
  atomic_number() const
  {
    return _atomic_number;
  }

  int
  natoms() const
  {
    return _matoms;
  }

  const int*
  aromatic() const
  {
    return _aromatic;
  }

  const int*
  symmetry_class() const
  {
    return _symmetry_class;
  }

  const int*
  atom_invariant() const
  {
    return _atom_invariant;
  }

  int
  atom_invariant(atom_number_t a) const
  {
    return _atom_invariant[a];
  }

  int
  aromatic(atom_number_t a) const
  {
    return _aromatic[a];
  }

  int
  nrings()
  {
    return _m->nrings();
  }

  const Ring*
  ringi(int r)
  {
    return _m->ringi(r);
  }

  const int*
  hydrogen_attached() const
  {
    return _hydrogen_attached;
  }

  const int*
  valid_start_atom() const
  {
    return _valid_start_atom;
  }

  int
  valid_start_atom(atom_number_t a) const
  {
    return _valid_start_atom[a];
  }

  int
  is_automatically_added(atom_number_t i) const
  {
    return _automatically_added[i];
  }

  int in_same_ring(atom_number_t, atom_number_t) const;

  int are_bonded(atom_number_t, atom_number_t) const;

  int canonical_bond_number(atom_number_t, atom_number_t) const;

  const Set_of_Atoms&
  doubly_bonded_singly_connected_neighbours(atom_number_t a) const
  {
    return _doubly_bonded_singly_connected_neighbour[a];
  }

  int
  append_bond_between_atoms(atom_number_t, atom_number_t, IWString&) const;

  char bond_symbol_between_atoms(atom_number_t, atom_number_t) const;

  const IWString&
  symbol_for_atom(atom_number_t a) const
  {
    return _atomic_symbol[a];
  }

  int
  ncon(atom_number_t i) const
  {
    return _atom[i]->ncon();
  }
};

int
PD_Molecule_Data::_fill_atom_symbol_array(const int* inv)
{
  for (int i = 0; i < _matoms; i++) {
    _atomic_symbol[i] << '[' << _atom_invariant[i] << ']';
  }

  return 1;
}

int
PD_Molecule_Data::_fill_atom_symbol_array(Molecule& m)
{
  int matoms = m.natoms();

  if (include_hydrogen_info_in_invariant) {
    for (int i = 0; i < matoms; i++) {
      const Element* e = m.elementi(i);

      if (only_include_hydrogen_info_on_heteroatoms && 6 == e->atomic_number()) {
        if (m.is_aromatic(i)) {
          _atomic_symbol[i] = e->aromatic_symbol();
        } else {
          _atomic_symbol[i] = e->symbol();
        }

        continue;
      }

      const_IWSubstring s;

      if (m.is_aromatic(i)) {
        s = e->aromatic_symbol();
      } else {
        s = e->symbol();
      }

      int h = _hydrogen_attached[i];

      if (0 == h) {
        _atomic_symbol[i] = s;
        continue;
      }

      if (1 == include_hydrogen_info_in_invariant) {
        _atomic_symbol[i] << '[' << s << "H]";
      } else {
        _atomic_symbol[i] << '[' << s << 'H' << h << ']';
      }
    }
  } else {
    for (int i = 0; i < matoms; i++) {
      const Element* e = m.elementi(i);

      if (m.is_aromatic(i)) {
        _atomic_symbol[i] = e->aromatic_symbol();
      } else {
        _atomic_symbol[i] = e->symbol();
      }
    }
  }

  return 1;
}

static int
identify_atoms_with_hydrogens(Molecule& m, int* h)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int hi = m.hcount(i);

    if (0 == hi) {
      h[i] = 0;
    } else if (1 == include_hydrogen_info_in_invariant) {
      h[i] = 1;
    } else {
      h[i] = hi;
    }
  }

  return 1;
}

int
PD_Molecule_Data::_looks_like_part_of_cf3(const Molecule& m, atom_number_t f) const
{
  atom_number_t c = _atom[f]->other(f, 0);

  const Atom* ac = _atom[c];

  if (4 != ac->ncon()) {
    return 0;
  }

  int halogen = 1;

  for (int i = 0; i < 4; i++) {
    atom_number_t j = ac->other(c, i);

    if (j == f) {
      continue;
    }

    if (9 == _atomic_number[j]) {
      halogen++;
    } else if (17 == _atomic_number[j]) {
      halogen++;
    }
  }

  return (3 == halogen);
}

int
PD_Molecule_Data::_identify_automatically_added(const Molecule& m)
{
  int rc = 0;

  for (int i = 0; i < _matoms; i++) {
    const Atom* ai = _atom[i];

    if (1 != ai->ncon()) {
      continue;
    }

    if (ai->nbonds() > ai->ncon()) {
      _automatically_added[i] = 1;
    } else if (_looks_like_part_of_cf3(m, i)) {
      _automatically_added[i] = 1;
    }
  }

  return rc;
}

static int
compute_canonical_bond_number(const Bond* b)
{
  if (b->is_aromatic()) {
    return 4;
  }
  if (b->is_single_bond()) {
    return 1;
  }
  if (b->is_double_bond()) {
    return 2;
  }
  if (b->is_triple_bond()) {
    return 3;
  }

  return 1;
}

/*
Slow version using atoms rather than bonds
int
PD_Molecule_Data::_identify_bonded_atoms (const Molecule & m)
{
  for (int i = 0; i < _matoms; i++)
  {
    const Atom * ai = _atom[i];

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = ai->item (j);

      atom_number_t k = b->other(i);

      if (k < i)
        continue;

      int c = compute_canonical_bond_number(b);

      _are_bonded[i * _matoms + k] = c;
      _are_bonded[k * _matoms + i] = c;
    }
  }

#ifdef ECHO_BONDING_INFO
  for (int i = 0; i < _matoms; i++)
  {
    for (int j = i + 1; j < _matoms; j++)
    {
      cerr << "Atom " << i << " and " << j << " bonded? " << _are_bonded[i * _matoms + j]
<< endl;
    }
  }
#endif

  return 1;
}*/

int
PD_Molecule_Data::_identify_bonded_atoms(const Molecule& m)
{
  int ne = m.nedges();

  for (int i = 0; i < ne; i++) {
    const Bond* b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    int c = compute_canonical_bond_number(b);

    _are_bonded[a1 * _matoms + a2] = c;
    _are_bonded[a2 * _matoms + a1] = c;
  }

#ifdef ECHO_BONDING_INFO
  for (int i = 0; i < _matoms; i++) {
    for (int j = i + 1; j < _matoms; j++) {
      cerr << "Atom " << i << " and " << j << " bonded? " << _are_bonded[i * _matoms + j]
           << endl;
    }
  }
#endif

  return 1;
}

static int
identify_in_same_ring(const Ring& r, int matoms, int* in_same_ring)
{
  int n = r.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = r[i];

    for (int k = i + 1; k < n; k++) {
      atom_number_t l = r[k];

      in_same_ring[j * matoms + l] = 1;
      in_same_ring[l * matoms + j] = 1;
    }
  }

  return 1;
}

static int
identify_in_same_ring(Molecule& m, int* in_same_ring)
{
  int matoms = m.natoms();

  int nr = m.nrings();

  for (int i = 0; i < nr; i++) {
    const Ring* ri = m.ringi(i);

    identify_in_same_ring(*ri, matoms, in_same_ring);
  }

  return 1;
}

/*
Slow version using atoms connections
int
PD_Molecule_Data::_fill_bond_symbols (Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi (i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item (j);

      if (b->is_single_bond())
        continue;

      if (b->is_aromatic())
        continue;

      atom_number_t k = b->other (i);

      char s;
      if (b->is_double_bond())
        s = '=';
      else
        s = '#';

      _bond_symbol[i * matoms + k] = s;
      _bond_symbol[k * matoms + i] = s;
    }
  }

  return 1;
}*/

int
PD_Molecule_Data::_fill_bond_symbols(Molecule& m)
{
  int ne = m.nedges();

  for (int i = 0; i < ne; i++) {
    const Bond* b = m.bondi(i);

    if (b->is_single_bond()) {
      continue;
    }

    if (b->is_aromatic()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    char s;
    if (b->is_double_bond()) {
      s = '=';
    } else {
      s = '#';
    }

    _bond_symbol[a1 * _matoms + a2] = s;
    _bond_symbol[a2 * _matoms + a1] = s;
  }

  return 1;
}

int
PD_Molecule_Data::_identify_doubly_bonded_singly_connected_neighbours(Molecule& m)
{
  assert(nullptr != _atom_invariant);

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    if (6 == a->atomic_number()) {
      continue;
    }

    if (1 != a->ncon()) {
      continue;
    }

    if (a->ncon() == a->nbonds()) {
      continue;
    }

    atom_number_t j = a->other(i, 0);

    _doubly_bonded_singly_connected_neighbour[j].add(i);

    _valid_start_atom[i] = 0;

    rc++;
  }

  // When there are multiple attachments, make sure they are in atomic number order

  for (int i = 0; i < matoms; i++) {
    Set_of_Atoms& si = _doubly_bonded_singly_connected_neighbour[i];

    if (si.number_elements() < 2) {
      continue;
    }

    int i0 = si[0];
    int i1 = si[1];

    atomic_number_t z0 = _atomic_number[i0];
    atomic_number_t z1 = _atomic_number[i1];

    if (z0 > z1) {
      si.swap_elements(0, 1);
    }
  }

  // Now adjust the invariant if there are attachments

  for (int i = 0; i < matoms; i++) {
    const Set_of_Atoms& si = _doubly_bonded_singly_connected_neighbour[i];

    if (0 == si.number_elements()) {
      continue;
    }

    atom_number_t si0 = si[0];

    atomic_number_t zsi0 = _atomic_number[si0];

    if (1 == si.number_elements()) {
      if (8 == zsi0) {
        _atom_invariant[i] *= 2;
      } else if (16 == zsi0) {
        _atom_invariant[i] *= 3;
      } else if (7 == zsi0) {
        _atom_invariant[i] *= 4;
      } else if (6 == zsi0) {
        _atom_invariant[i] *= 5;
      } else {
        _atom_invariant[i] *= 6;
      }
    } else if (2 == si.number_elements()) {
      atomic_number_t zsi1 = _atomic_number[si[1]];

      if (8 == zsi0 && 8 == zsi1) {  // O=*=O
        _atom_invariant[i] *= 7;
      } else if (16 == zsi0 && 16 == zsi1) {  // S=*=S
        _atom_invariant[i] *= 8;
      } else if ((8 == zsi0 && 7 == zsi1) || (7 == zsi0 && 8 == zsi1)) {  // O=*=N
        _atom_invariant[i] *= 9;
      } else {
        _atom_invariant[i] *= 10;
      }
    }
  }

  return rc;
}

int
PD_Molecule_Data::_identify_valid_start_atoms(Molecule& m)
{
  int matoms = m.natoms();

  const int* symmetry_class = m.symmetry_classes();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (0 == _valid_start_atom[i]) {
      continue;
    }

    if (min_heteroatoms_needed_in_path > 0 && 6 == m.atomic_number(i)) {
      _valid_start_atom[i] = 0;
      continue;
    }

    if (m.nrings(i) > 1) {
      _valid_start_atom[i] = 0;
      continue;
    }

    for (int j = i + 1; j < matoms; j++) {
      //    cerr << "Symmetryj " << symmetry_class[j] << " vs " << symmetry_class[i] <<
      //    endl;
      if (symmetry_class[j] == symmetry_class[i]) {
        _valid_start_atom[j] = 0;
      }
    }

    rc++;
  }

  // iw_write_array (_valid_start_atom, matoms, "_valid_start_atom", cerr);

  return rc;
}

static int
compute_ring_bonds(const Atom* a)
{
  int acon = a->ncon();

  if (2 == acon) {  // atom A is in a ring
    return 2;
  }

  int rc = 0;

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    if (b->nrings()) {
      rc++;
    }
  }

  return rc;
}

int
PD_Molecule_Data::_assign_atomic_invariants(Molecule& m)
{
  if (atom_typing_specification.active()) {
    return atom_typing_specification.assign_atom_types(m, _atom_invariant);
  }

  for (int i = 0; i < _matoms; i++) {
    if (m.is_aromatic(i)) {
      _aromatic[i] = 1;
    } else {
      _aromatic[i] = 0;
    }

    _atom_invariant[i] = 1000000 * _atomic_number[i] + 900000 * _aromatic[i];

    if (include_hydrogen_info_in_invariant) {
      _atom_invariant[i] += 290000 * _hydrogen_attached[i];
    }

    if (m.is_ring_atom(i)) {
      int r = compute_ring_bonds(_atom[i]);
      if (r > 2) {
        _atom_invariant[i] += 150000;
      }
    }
  }

  return 1;
}

PD_Molecule_Data::PD_Molecule_Data(Molecule& m)
{
  _m = &m;

  m.compute_aromaticity_if_needed();

  _matoms = m.natoms();

  // Must do hydrogen array before atomic symbols

  if (include_hydrogen_info_in_invariant) {
    _hydrogen_attached = new int[_matoms];

    identify_atoms_with_hydrogens(m, _hydrogen_attached);
  } else {
    _hydrogen_attached = nullptr;
  }

  _valid_start_atom = new_int(_matoms, 1);

  _in_same_ring = new_int(_matoms * _matoms);

  assert(nullptr != _in_same_ring);

  identify_in_same_ring(m, _in_same_ring);

  _doubly_bonded_singly_connected_neighbour = new Set_of_Atoms[_matoms];

  _identify_valid_start_atoms(m);

  _bond_symbol = new char[_matoms * _matoms];
  set_vector(_bond_symbol, _matoms * _matoms, ' ');

  assert(nullptr != _bond_symbol);

  _fill_bond_symbols(m);

  _atomic_number = new atom_number_t[_matoms];

  m.atomic_numbers(_atomic_number);

  _aromatic = new int[_matoms];

  _symmetry_class = new int[_matoms];

  copy_vector(_symmetry_class, m.symmetry_classes(), _matoms);

  _atom = new const Atom*[_matoms];

  m.atoms(_atom);

  _automatically_added = new_int(_matoms);

  _identify_automatically_added(m);

  _are_bonded = new_int(_matoms * _matoms);

  _identify_bonded_atoms(m);

  _atom_invariant = new int[_matoms];

  _assign_atomic_invariants(m);

  // Must be done after the atom invariant is computed because it alters
  // the atom invariant

  _identify_doubly_bonded_singly_connected_neighbours(m);

  _atomic_symbol = new IWString[_matoms];

  if (!atom_typing_specification.active()) {
    _fill_atom_symbol_array(m);
  } else {
    _fill_atom_symbol_array(_atom_invariant);
  }

  return;
}

PD_Molecule_Data::~PD_Molecule_Data()
{
  if (nullptr != _atomic_symbol) {
    delete[] _atomic_symbol;
  }

  if (nullptr != _valid_start_atom) {
    delete[] _valid_start_atom;
  }

  if (nullptr != _in_same_ring) {
    delete[] _in_same_ring;
  }

  if (nullptr != _atomic_number) {
    delete[] _atomic_number;
  }

  if (nullptr != _doubly_bonded_singly_connected_neighbour) {
    delete[] _doubly_bonded_singly_connected_neighbour;
  }

  if (nullptr != _bond_symbol) {
    delete[] _bond_symbol;
  }

  if (nullptr != _aromatic) {
    delete[] _aromatic;
  }

  if (nullptr != _atom_invariant) {
    delete[] _atom_invariant;
  }

  if (nullptr != _atom) {
    delete[] _atom;
  }

  if (nullptr != _automatically_added) {
    delete[] _automatically_added;
  }

  if (nullptr != _symmetry_class) {
    delete[] _symmetry_class;
  }

  if (nullptr != _are_bonded) {
    delete[] _are_bonded;
  }

  if (nullptr != _hydrogen_attached) {
    delete[] _hydrogen_attached;
  }

  return;
}

int
PD_Molecule_Data::in_same_ring(atom_number_t a1, atom_number_t a2) const
{
  return _in_same_ring[a1 * _matoms + a2];
}

int
PD_Molecule_Data::are_bonded(atom_number_t a1, atom_number_t a2) const
{
  return _are_bonded[a1 * _matoms + a2];
}

int
PD_Molecule_Data::canonical_bond_number(atom_number_t a1, atom_number_t a2) const
{
  assert(a1 >= 0 && a1 < _matoms);
  assert(a2 >= 0 && a2 < _matoms);

  return _are_bonded[a1 * _matoms + a2];
}

int
PD_Molecule_Data::append_bond_between_atoms(atom_number_t a1, atom_number_t a2,
                                            IWString& s) const
{
  if (' ' != _bond_symbol[a1 * _matoms + a2]) {
    s << _bond_symbol[a1 * _matoms + a2];
  }

  return 1;
}

char
PD_Molecule_Data::bond_symbol_between_atoms(atom_number_t a1, atom_number_t a2) const
{
  return _bond_symbol[a1 * _matoms + a2];
}

class Path_so_Far
{
 private:
  int* _path;

  int* _bond;

  int _atoms_in_path;

  //  Quick indicator of whether or not a given atom has been included in the path before

  int* _in_path;

  //  We can have multiple rings in our path, so item in the path needs to know the
  //  identity of a possible ring closure. Note that because we are tracing linear
  //  paths, there will be a max of 2 ring closures at any atom. Would be 1 without spiro
  //  fusions..

  int* _index_of_ring_closure1;
  int* _index_of_ring_closure2;

  // In order to get unique smiles, we first load into a square array
  // the pairs of atoms where we need to put ring numbers in the smiles

  int* _ring_closure_needed;

  //  At each step along the path, we keep track of any =O type groups that must be added

  Set_of_Atoms* _singly_connected_atoms_to_add;

  //  When we are deciding on canonical orders, we need a set of atomic invariants

  int* _invariant;

  //  When doing canonicaliation of ring systems, we need to know which atoms
  //  are in rings within the path

  int* _in_ring_in_path;

  //  When producing a smiles, we need to know the ring number associated with
  //  each index

  int* _ring_number;

  //  It can be helpful to know the order of the atoms in any smiles produced

  resizable_array<int> _indices_in_smiles;

  // when forming a smiles, ring openings and closings are complicated.
  // Therefore each caller of _produce_smiles is responsible for setting
  // up the ring opening and closing string for each atom

  IWString* _ring_opening_and_closing;

  // When doing ring numbering, we need to know where a given path
  // member is in the _indices_in_smiles array

  int* _position_in_smiles;

  //  private functions

  int
  _add(const PD_Molecule_Data&, const Bond*, atom_number_t);

  int
  _assign_identifier_and_add_to_paths(const IWString& smiles,
                                      resizable_array<int>& paths_this_molecule) const;
  int
  _fill_invariant_array_no_ring_consideration(const PD_Molecule_Data& pdmd);
  int
  _fill_invariant_array_with_ring_consideration(const PD_Molecule_Data& pdmd);

  void
  _append_atom_and_its_connections(const PD_Molecule_Data& pdmd, int ndx,
                                   IWString& smiles) const;

  void
  _tell_atom_ring_no_longer_exists(int ndx, int ndx_being_removed);

  int
  _add_ring_opening_and_closing_characters(const PD_Molecule_Data& pdmd);

  void
  _place_ring_opening_symbols_between_indices(const PD_Molecule_Data& pdmd, int ndx1,
                                              int ndx2, int& ring_number) const;

  int
  _choose_direction_around_ring(int ndx, int start_index_of_ring,
                                int end_index_of_ring) const;
  int
  _produce_smiles(const PD_Molecule_Data& pdmd, IWString& smiles);
  int
  _produce_smiles_with_rings(const PD_Molecule_Data& pdmd,
                             resizable_array<int>& paths_this_molecule);
  int
  _produce_smiles_with_rings(const PD_Molecule_Data& pdmd, int direction,
                             resizable_array<int>& paths_this_molecule);
  int
  _produce_smiles_no_ring(const PD_Molecule_Data& pdmd,
                          resizable_array<int>& paths_this_molecule);

  int
  _produce_smiles_chain_at_both_ends(const PD_Molecule_Data& pdmd,
                                     resizable_array<int>& paths_this_molecule);
  int
  _produce_smiles_chain_at_both_ends(const PD_Molecule_Data& pdmd, int direction,
                                     resizable_array<int>& paths_this_molecule);
  int
  _produce_smiles_with_rings_left_to_right(const PD_Molecule_Data& pdmd,
                                           int index_of_end_of_first_ring,
                                           resizable_array<int>& paths_this_molecule);
  int
  _produce_smiles_with_rings_right_to_left(const PD_Molecule_Data& pdmd,
                                           int index_of_end_of_first_ring,
                                           resizable_array<int>& paths_this_molecule);
  int
  _produce_smiles_with_one_ring(const PD_Molecule_Data& pdmd,
                                resizable_array<int>& paths_this_molecule);
  /*  int _produce_smiles_pure_ring (const PD_Molecule_Data & pdmd,
                                     resizable_array<int> & paths_this_molecule);*/
  /*  int _produce_smiles_chain_at_end (const PD_Molecule_Data & pdmd,
                                        resizable_array<int> & paths_this_molecule,
                                        int start_index_of_ring,
                                        int end_index_of_ring);
      int _produce_smiles_chain_at_beginning (const PD_Molecule_Data & pdmd,
                                        resizable_array<int> & paths_this_molecule,
                                        int start_index_of_ring,
                                        int end_index_of_ring);*/
  int
  _produce_smiles_ring_and_chain(const PD_Molecule_Data& pdmd,
                                 resizable_array<int>& paths_this_molecule, int istart,
                                 int istep);
  int
  _produce_smiles_with_multiple_rings(const PD_Molecule_Data& pdmd, int nr,
                                      resizable_array<int>& paths_this_molecule);

  int
  _determine_direction_chain(const PD_Molecule_Data& pdmd) const;
  int
  _determine_chain_length(int istart, int istep) const;
  int
  _determine_direction_rings_at_both_ends(const PD_Molecule_Data& mdmd,
                                          int& too_complicated) const;
  int
  _determine_direction_around_ring(const PD_Molecule_Data& pdmd, int ndx) const;

  int
  _identify_ring_atoms_in_path();
  void
  _identify_ring_start_and_end(int& start_index_of_ring, int& end_index_of_ring) const;

  int
  _find_corresponding_ring_join(int ndx, const int* a, int avoid) const;

  int
  _process_isolated_rings(const PD_Molecule_Data& pdmd,
                          resizable_array<int>& paths_this_molecule);
  int
  _process_fused_rings(const PD_Molecule_Data& pdmd, int direction,
                       resizable_array<int>& paths_this_molecule);

  int
  _produce_smiles_for_this_part_of_path(const PD_Molecule_Data& pdmd, int istart,
                                        int istop, int istep, IWString& smiles);
  int
  _produce_smiles_starting_at(const PD_Molecule_Data& pdmd, int istart, int istep,
                              IWString& smiles);
  void
  _deal_with_ring_forming_connection(int ifrom, int ito, int& next_ring_number_to_assign,
                                     IWString& smiles);
  int
  _produce_smiles_from_chain_to_rings(const PD_Molecule_Data& pdmd, int direction,
                                      resizable_array<int>& paths_this_molecule);
  int
  _count_ring_openings(int) const;
  int
  _determine_ring_end(int ndx) const;
  int
  _path_is_a_whole_ring_or_ring_system(const PD_Molecule_Data& pdmd) const;

 public:
  Path_so_Far(int matoms);
  ~Path_so_Far();

  int
  debug_print(std::ostream&) const;

  int
  atoms_in_path() const
  {
    return _atoms_in_path;
  }

  int
  rings_present() const;

  int
  last_atom_in_path() const
  {
    return _path[_atoms_in_path - 1];
  }

  int
  add(const PD_Molecule_Data&, const Bond*, atom_number_t);

  int
  is_in_path(atom_number_t s) const
  {
    return _in_path[s];
  }

  int
  pop();

  int
  heteroatoms_in_path(const PD_Molecule_Data&) const;

  int
  produce_smiles(const PD_Molecule_Data& pdmd, resizable_array<int>& paths_this_molecule);
  int
  produce_smiles(const PD_Molecule_Data& pdmd, int direction, IWString& smiles,
                 int istart, int istop) const;

  int
  produce_smiles(const PD_Molecule_Data& pdmd, int direction, IWString& smiles) const;
  int
  produce_smiles(const PD_Molecule_Data& pdmd, int direction, IWString& smiles,
                 int natoms) const;

  int
  write_matched_atoms(int direction, IWString& buffer) const;

  const resizable_array<int>&
  indices_in_smiles() const
  {
    return _indices_in_smiles;
  }

  atom_number_t
  path_member(int i) const
  {
    return _path[i];
  }

  int
  bondi(int i) const
  {
    return _bond[i];
  }
};

Path_so_Far::Path_so_Far(int matoms)
{
  assert(max_path_length > 0);
  assert(matoms > 0);

  _path = new int[max_path_length];

  _bond = new int[max_path_length];

  _atoms_in_path = 0;

  _index_of_ring_closure1 = new int[max_path_length];
  _index_of_ring_closure2 = new int[max_path_length];

  _ring_closure_needed = new int[max_path_length * max_path_length];

  _in_path = new_int(matoms);

  _position_in_smiles = new int[matoms];

  _singly_connected_atoms_to_add = new Set_of_Atoms[max_path_length];

  _invariant = new int[max_path_length];

  _ring_number = new int[max_path_length];

  _ring_opening_and_closing = new IWString[max_path_length];

  _in_ring_in_path = new int[max_path_length];

  return;
}

Path_so_Far::~Path_so_Far()
{
  if (nullptr != _path) {
    delete[] _path;
  }

  if (nullptr != _bond) {
    delete[] _bond;
  }

  if (nullptr != _in_path) {
    delete[] _in_path;
  }

  if (nullptr != _singly_connected_atoms_to_add) {
    delete[] _singly_connected_atoms_to_add;
  }

  if (nullptr != _invariant) {
    delete[] _invariant;
  }

  if (nullptr != _ring_number) {
    delete[] _ring_number;
  }

  if (nullptr != _ring_opening_and_closing) {
    delete[] _ring_opening_and_closing;
  }

  if (nullptr != _index_of_ring_closure1) {
    delete[] _index_of_ring_closure1;
  }

  if (nullptr != _index_of_ring_closure2) {
    delete[] _index_of_ring_closure2;
  }

  delete[] _ring_closure_needed;

  if (nullptr != _in_ring_in_path) {
    delete[] _in_ring_in_path;
  }

  if (nullptr != _position_in_smiles) {
    delete[] _position_in_smiles;
  }

  _atoms_in_path = -128;

  return;
}

int
Path_so_Far::debug_print(std::ostream& os) const
{
  os << _atoms_in_path << " atoms in path\n";

  for (int i = 0; i < _atoms_in_path; i++) {
    os << " i = " << i << " atom " << _path[i];
    if (i != _atoms_in_path - 1) {
      os << " bond " << _bond[i];
    }

    if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[i]) {
      os << " R1 " << _path[_index_of_ring_closure1[i]];
    }
    if (INVALID_ATOM_NUMBER != _index_of_ring_closure2[i]) {
      os << " R2 " << _path[_index_of_ring_closure2[i]];
    }

    os << endl;
  }

  return os.good();
}

int
Path_so_Far::heteroatoms_in_path(const PD_Molecule_Data& pdmd) const
{
  const atomic_number_t* z = pdmd.atomic_number();

  int rc = 0;

  for (int i = 0; i < _atoms_in_path; i++) {
    atom_number_t a = _path[i];

    if (6 != z[a]) {
      rc++;
    }
  }

  return rc;
}

int
Path_so_Far::_add(const PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom)
{
// #define DEBUG_ADD_ATOM_TO_PATH
#ifdef DEBUG_ADD_ATOM_TO_PATH
  cerr << "Adding atom " << zatom << " to path of size " << _atoms_in_path << " atoms\n";
  debug_print(cerr);
#endif

  _index_of_ring_closure1[_atoms_in_path] = INVALID_ATOM_NUMBER;
  _index_of_ring_closure2[_atoms_in_path] = INVALID_ATOM_NUMBER;

  if (nullptr == b) {
    ;
  } else if (_atoms_in_path < 2) {  // smallest ring is 3 membered ring
    ;
  } else if (b->nrings())  // look for any rings formed with previous members of the path
  {
    resizable_array<int> ring_closure_bonds;

    for (int i = 0; i < _atoms_in_path - 1; i++)  // no 2 membered rings!
    {
      atom_number_t j = _path[i];

      if (!pdmd.are_bonded(zatom, j)) {
        continue;
      }

      if (pdmd.in_same_ring(zatom, j)) {
        ring_closure_bonds.add(i);
      }
    }

    int n = ring_closure_bonds.number_elements();

    if (n > 2) {  // ring system is too complex, more that 2 ring closures at ZATOM
      return 0;
    }

    for (int i = 0; i < n; i++) {
      int j = ring_closure_bonds[i];

      // Would form something with too many ring closures at J
      if (_index_of_ring_closure2[j] >= 0) {
        return 0;
      }
    }

    //  It is now safe to add all the ring closure bonds

    for (int i = 0; i < n; i++) {
      int j = ring_closure_bonds[i];

#ifdef DEBUG_ADD_ATOM_TO_PATH
      cerr << "Adding ring bond between atom " << _path[j] << " and " << zatom << endl;
#endif
      if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[j]) {
        _index_of_ring_closure1[j] = _atoms_in_path;
      } else if (INVALID_ATOM_NUMBER == _index_of_ring_closure2[j]) {
        _index_of_ring_closure2[j] = _atoms_in_path;
      } else {  // really should not happen
        return 0;
      }

      if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[_atoms_in_path]) {
        _index_of_ring_closure1[_atoms_in_path] = j;
      } else if (INVALID_ATOM_NUMBER == _index_of_ring_closure2[_atoms_in_path]) {
        _index_of_ring_closure2[_atoms_in_path] = j;
      } else {
        return 0;
      }
    }
  }

  assert(_atoms_in_path < max_path_length);

  _path[_atoms_in_path] = zatom;

  _singly_connected_atoms_to_add[_atoms_in_path] =
      pdmd.doubly_bonded_singly_connected_neighbours(zatom);

  _in_path[zatom] = 1;

  if (0 == _atoms_in_path) {
    ;
  } else if (nullptr != b) {
    _bond[_atoms_in_path - 1] = compute_canonical_bond_number(b);
  }

  _atoms_in_path++;

#ifdef DEBUG_ADD_ATOM_TO_PATH
  cerr << "After adding atom " << zatom << " path is now\n";
  debug_print(cerr);
#endif

  return 1;
}

int
Path_so_Far::add(const PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom)
{
  if (_in_path[zatom]) {
    return 0;
  }

  if (_atoms_in_path >= max_path_length) {
    return 0;
  }

  if (!_add(pdmd, b, zatom)) {
    return 0;
  }

  return _atoms_in_path;
}

/*
  Each path point can have up to two ring closures. If the items
  being removed is our 2nd one, just invalidate it.
  If it is our first one, then copy our 2nd ring closure to our
  first and invalidate number 2
*/

void
Path_so_Far::_tell_atom_ring_no_longer_exists(int ndx, int ndx_being_removed)
{
  if (ndx_being_removed == _index_of_ring_closure2[ndx]) {
    _index_of_ring_closure2[ndx] = INVALID_ATOM_NUMBER;
  } else  // is first ring closure, copy target's 2 to 1
  {
    if (_index_of_ring_closure1[ndx] != ndx_being_removed) {
      cerr << "Yipes, removing item " << ndx_being_removed << ", atom "
           << _path[ndx_being_removed] << " bond to " << ndx << ", atom " << _path[ndx]
           << " not bonded!\n";
      debug_print(cerr);
    }
    assert(_index_of_ring_closure1[ndx] == ndx_being_removed);

    _index_of_ring_closure1[ndx] = _index_of_ring_closure2[ndx];
    _index_of_ring_closure2[ndx] = INVALID_ATOM_NUMBER;
  }

  return;
}

int
Path_so_Far::pop()
{
  assert(_atoms_in_path > 0);

#ifdef DEBUG_POP
  cerr << "Pop operation\n";
  debug_print(cerr);
#endif

  // If this last atom formed any rings, we must let those other atoms know that
  // this ring is now invalid

  const int index_last_atom = _atoms_in_path - 1;

  if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[index_last_atom]) {
    int i = _index_of_ring_closure1[index_last_atom];

    _tell_atom_ring_no_longer_exists(i, index_last_atom);

    _index_of_ring_closure1[index_last_atom] = INVALID_ATOM_NUMBER;
  }

  if (INVALID_ATOM_NUMBER != _index_of_ring_closure2[index_last_atom]) {
    int i = _index_of_ring_closure2[index_last_atom];

    _tell_atom_ring_no_longer_exists(i, index_last_atom);

    _index_of_ring_closure2[index_last_atom] = INVALID_ATOM_NUMBER;
  }

  atom_number_t rc = _path[index_last_atom];

  _atoms_in_path--;

  _in_path[rc] = 0;

  return rc;
}

int
Path_so_Far::_assign_identifier_and_add_to_paths(
    const IWString& smiles, resizable_array<int>& paths_this_molecule) const
{
// #define DEBUG_ASSIGN_IDENTIFIER_AND_ADD_TO_PATHS
#ifdef DEBUG_ASSIGN_IDENTIFIER_AND_ADD_TO_PATHS
  cerr << "_assign_identifier_and_add_to_paths: smiles '" << smiles << "' from\n";
  debug_print(cerr);
  if (die_when_formed.length() && die_when_formed == smiles) {
    cerr << "Formed '" << die_when_formed << "'\n";
    abort();
  }
#endif

  unsigned int s = assign_unique_identifier(path_to_unique_id, smiles);

  paths_this_molecule.add(s);

  return 1;
}

int
Path_so_Far::_fill_invariant_array_no_ring_consideration(const PD_Molecule_Data& pdmd)
{
  const int* atom_invariant = pdmd.atom_invariant();

  for (int i = 0; i < _atoms_in_path; i++) {
    atom_number_t a = _path[i];

    int zextra = _singly_connected_atoms_to_add[i].number_elements();

    _invariant[i] = atom_invariant[a] + 200 * zextra;

    if (0 == i) {
      _invariant[i] += _bond[i];
    } else if (_atoms_in_path - 1 == i) {
      _invariant[i] += _bond[_atoms_in_path - 2];
    } else {
      _invariant[i] += _bond[i - 1] + _bond[i];
    }
  }

  return 1;
}

int
Path_so_Far::_fill_invariant_array_with_ring_consideration(const PD_Molecule_Data& pdmd)
{
  const int* atom_invariant = pdmd.atom_invariant();

  for (int i = 0; i < _atoms_in_path; i++) {
    atom_number_t j = _path[i];

    int zextra = _singly_connected_atoms_to_add[i].number_elements();

    _invariant[i] = atom_invariant[j] + 30000 * zextra;

    int b = 0;
    if (0 == i) {
      b = _bond[0];
      // is first bonded to last?
      int c = pdmd.canonical_bond_number(j, _path[_atoms_in_path - 1]);
      if (c) {
        b += c;
      }
    } else if (_atoms_in_path - 1 == i) {
      b = _bond[_atoms_in_path - 2];
      int c = pdmd.canonical_bond_number(j, _path[0]);  // is last bonded to first?
      if (c) {
        b += c;
      }
    } else {
      b = _bond[i - 1] + _bond[i];
    }

    if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[i]) {
      _invariant[i] += 14000;
      int c = pdmd.canonical_bond_number(j, _path[_index_of_ring_closure1[i]]);
      b += c;
      if (INVALID_ATOM_NUMBER != _index_of_ring_closure2[i]) {
        c = pdmd.canonical_bond_number(j, _path[_index_of_ring_closure2[i]]);
        b += c;
        _invariant[i] += 14000;
      }
    }

    _invariant[i] += b;
  }

  return 1;
}

int
Path_so_Far::rings_present() const
{
  if (_atoms_in_path < 3) {
    return 0;
  }

  for (int i = 0; i < _atoms_in_path; i++) {
    if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[i]) {
      return 1;
    }
  }

  return 0;
}

/*
  Chains are easy to canonicalise

  Rings can be difficult to canonicalise.

  First check, if the two ends are both not involved with any ring stuff,
  then it must be an ortho substituted form of ring or ring system.

  ___  ___    _____       _______
     \/   |__|     |     |
                   |     |
                   | ... |
                   |     |
                   |_____|

  These are fairly easy, because the canonicalisation will be one direction or
  the other.

  If just one end is not involved in a ring, then the canonicalisation
  will start at that end.  The final item in the path will be a ring,
  and we'll need to decide which way around to traverse that ring

                     __
  ___   ____________/  \
     \_/            \__/

  That leaves the case of a ring at either end

  C1CCCCC1CC1CC1CC1CCCCC1

  We first examine the centre chain section.  Can that be
  canonicalised.  If so, that will set the direction of traversal of
  the path.

  We decide on a direction around each ring separately.  As we do
  that, we generate a hopefully unique number for each ring.  If the
  direction isn't set by the central path, we use these numbers to set
  a direction.
*/

int
Path_so_Far::produce_smiles(const PD_Molecule_Data& pdmd,
                            resizable_array<int>& paths_this_molecule)
{
  assert(_atoms_in_path > 0);

  _indices_in_smiles.resize_keep_storage(0);

  for (int i = 0; i < _atoms_in_path; i++) {
    _ring_opening_and_closing[i].resize_keep_storage(0);
  }

  int rc;

  if (!rings_present()) {
    _fill_invariant_array_no_ring_consideration(pdmd);
    rc = _produce_smiles_no_ring(pdmd, paths_this_molecule);
    assert(_indices_in_smiles.number_elements() == _atoms_in_path);
  } else if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[0] &&
             INVALID_ATOM_NUMBER == _index_of_ring_closure1[_atoms_in_path - 1]) {
    _fill_invariant_array_no_ring_consideration(pdmd);
    rc = _produce_smiles_chain_at_both_ends(pdmd, paths_this_molecule);
  } else {
    rc = _produce_smiles_with_rings(pdmd, paths_this_molecule);
  }

  if (0 == rc) {
    return 0;
  }

  if (stream_for_paths.active() && _indices_in_smiles.number_elements() > 0) {
    int p = paths_this_molecule.last_item();

    const IWString* s = id_to_path[p];

    stream_for_paths << (*s);

    for (int i = 0; i < _atoms_in_path; i++) {
      int j = _indices_in_smiles[i];

      atom_number_t k = _path[j];

      iwdigits.append_number(stream_for_paths, k);
    }

    stream_for_paths << '\n';
  }

  return rc;
}

/*
  We are dealing with a path with multiple rings. Figure out which atoms
  are in rings formed by the path
*/

/*int
Path_so_Far::_identify_ring_atoms_in_path()
{
  atom_number_t r1 = INVALID_ATOM_NUMBER;
  atom_number_t r2 = INVALID_ATOM_NUMBER;

  int ring_membership = 0;

  for (int i = 0; i < _atoms_in_path; i++)
  {
    _in_ring_in_path[i] = ring_membership;

    if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[i])
      continue;

    if (_path[i] == r1)    // we have closed r1
    {
      r1 = INVALID_ATOM_NUMBER;     // reset for next use
      ring_membership--;
    }
    else if (_path[i] == r2)    // we have closed r2
    {
      r2 = INVALID_ATOM_NUMBER;
      ring_membership--;
    }
    else if (INVALID_ATOM_NUMBER == r1)    // new ring, can we use r1
    {
      r1 = _path[_index_of_ring_closure1[i]];
      ring_membership++;
//    cerr << "From atom " <<_path[i] << " set ring 1 " << r2 << endl;
    }
    else if (INVALID_ATOM_NUMBER == r2)    // new ring, can we use r2
    {
      r2 = _path[_index_of_ring_closure1[i]];
      ring_membership++;
//    cerr << "From atom " <<_path[i] << " set ring 2 " << r2 << endl;
    }
    else
    {
      cerr << "Path_so_Far::_identify_ring_atoms_in_path:bad news, atom " << _path[i] <<
", r1 " << r1 << " r2 " << r2 << endl; debug_print (cerr); abort();
    }

    if (INVALID_ATOM_NUMBER == _index_of_ring_closure2[i])
      continue;

    if (_path[i] == r1)
    {
      r1 = INVALID_ATOM_NUMBER;
      ring_membership--;
    }
    else if (_path[i] == r2)    // we have closed r2
    {
      r2 = INVALID_ATOM_NUMBER;
      ring_membership--;
    }
    else if (INVALID_ATOM_NUMBER == r1)    // new ring, can we use r1
    {
      r1 = _path[_index_of_ring_closure2[i]];
      ring_membership++;
    }
    else if (INVALID_ATOM_NUMBER == r2)    // new ring, can we use r2
    {
      r2 = _path[_index_of_ring_closure2[i]];
      ring_membership++;
    }
    else
    {
      cerr << "Path_so_Far::_identify_ring_atoms_in_path:bad news2, atom " << _path[i] <<
endl; debug_print (cerr); abort();
    }
  }

  assert (0 == ring_membership);

  return 1;
}*/

/*
  Things get very complicated with 2 rings.

  There are two possibilities. We may have a fused system,
  or something like a biphenyl. From either ring, we may have
  a chain branching out

  We simplify the fused ring case by only processing them
  if there are chain atoms at one or both ends.
*/

/*int
Path_so_Far::_produce_smiles_with_multiple_rings (const PD_Molecule_Data & pdmd,
                                                  int nr,
                                                  resizable_array<int> &
paths_this_molecule)
{
  if (nr > 2)
  {
    cerr << "Sorry, cannot do " << nr << " rings '" << pdmd.molecule()->name() << "'\n";
    return 0;
  }

  _identify_ring_atoms_in_path();

  int index_of_first_ring_atom = -1;
  int index_of_last_ring_atom = -1;
  int fused_system = 0;

  for (int i = 0; i < _atoms_in_path; i++)
  {
    if (_in_ring_in_path[i] > 1)    // must be fused rings present
      fused_system = 1;

    if (0 == _in_ring_in_path[i])
      continue;

    index_of_last_ring_atom = i;

    if (index_of_first_ring_atom < 0)
      index_of_first_ring_atom = i;
  }

  if (! fused_system)
    return _process_isolated_rings (pdmd, paths_this_molecule);

// Only process fused rings if chain at one or both ends

  int lhs_chain_length = index_of_first_ring_atom;
  int rhs_chain_length = _atoms_in_path - index_of_last_ring_atom - 1;

  if (0 == lhs_chain_length && 0 == rhs_chain_length)
    return 0;

  if (lhs_chain_length > rhs_chain_length)
    return _process_fused_rings (pdmd, 1, paths_this_molecule);

  if (lhs_chain_length < rhs_chain_length)
    return _process_fused_rings (pdmd, -1, paths_this_molecule);

  int d = _determine_direction_chain (pdmd);

  return _process_fused_rings (pdmd, d, paths_this_molecule);
}*/

/*
  Fused ring. We know the direction. When we come to the first ring atom,
  we need to figure out the direction around the ring system
*/

/*int
Path_so_Far::_process_fused_rings (const PD_Molecule_Data & pdmd,
                                   int direction,
                                   resizable_array<int> & paths_this_molecule)
{
  int istart, istop, istep;
  if (direction > 0)
  {
    istart = 0;
    istop = _atoms_in_path;
    istep = 1;
  }
  else
  {
    istart = _atoms_in_path - 1;
    istop = -1;
    istep = -1;
  }

// The ring numbers will be assigned in this order

  char * ring_number_to_assign = "1221";

  int rnbndx = 0;   // index into the ring_number_to_assign string

  for (int i = istart; i != istop; i += istep)
  {
    _indices_in_smiles.add (i);

    if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[i])
      continue;

    _ring_opening_and_closing[i] << ring_number_to_assign[rnbndx];
    rnbndx++;

    if (INVALID_ATOM_NUMBER == _index_of_ring_closure2[i])
      continue;

    _ring_opening_and_closing[i] << ring_number_to_assign[rnbndx];
    rnbndx++;
  }

  IWString smiles;

  _produce_smiles (pdmd, smiles);

  return _assign_identifier_and_add_to_paths (smiles, paths_this_molecule);
}*/

/*
 */

/*int
Path_so_Far::_process_isolated_rings (const PD_Molecule_Data & pdmd,
                                      resizable_array<int> & paths_this_molecule)
{
  int index_of_first_ring_atom = -1;
  int index_of_last_ring_atom = -1;

  for (int i = 0; i < _atoms_in_path; i++)
  {
    if (_in_ring_in_path[i] > 1)    // must be fused rings present
      return 0;

    if (0 == _in_ring_in_path[i])
      continue;

    index_of_last_ring_atom = i;

    if (index_of_first_ring_atom < 0)
      index_of_first_ring_atom = i;
  }

// Great, we have two, isolated rings. Must be of the form
//  <optional chain> <ring> <chain of at least one atom> <ring> <optional chain>

  int initial_chain_length = index_of_first_ring_atom;
  int end_chain_length = _atoms_in_path - index_of_last_ring_atom - 1;

  if (initial_chain_length > end_chain_length)
    return _process_isolated_rings (pdmd, 1, paths_this_molecule);

  if (initial_chain_length < end_chain_length)
    return _process_isolated_rings (pdmd, -1, paths_this_molecule);

// We'll never find pure biphenyl, but why would we want to.

  if (0 == initial_chain_length)
    return 0;

  if (_processed_by_differentiating_end_chains (pdmd, paths_this_molecule))
    return 1;

// Cannot discern a direction from the chains.


  for (int i = 0; i < _atoms_in_path; i++)
  {
    _indices_in_smiles.add (i);
  }

  IWString smiles;

  _produce_smiles (pdmd, smiles);

  return _assign_identifier_and_add_to_paths (smiles, paths_this_molecule);
}*/

/*
  We know the direction in which we'll process this path, but we need
  to figure out the canonical order around the two rings
*/

/*int
Path_so_Far::_process_isolated_rings (const PD_Molecule_Data & pdmd,
                                      int direction,
                                      resizable_array<int> & paths_this_molecule)
{
  int istart, istop, istep;
  if (direction > 0)
  {
    istart = 0;
    istop = _atoms_in_path;
    istep = 1;
  }
  else
  {
    istart = _atoms_in_path - 1;
    istop = -1;
    istep = -1;
  }

  for (int i = istart; i != istop; i += istep)
  {
    _indices_in_smiles.add (i);

    if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[i])
      continue;

    _ring_opening_and_closing[i] = '1';    // only one ring number needed
  }

  IWString smiles;

  _produce_smiles (pdmd, smiles)

  return _assign_identifier_and_add_to_paths (smiles, paths_this_molecule);
}*/

void
Path_so_Far::_append_atom_and_its_connections(const PD_Molecule_Data& pdmd, int ndx,
                                              IWString& smiles) const
{
  atom_number_t a = _path[ndx];

  smiles << pdmd.symbol_for_atom(a);

  int n = _singly_connected_atoms_to_add[ndx].number_elements();

  for (int i = 0; i < n; i++) {
    smiles << '(';
    atom_number_t j = _singly_connected_atoms_to_add[ndx].item(i);
    pdmd.append_bond_between_atoms(a, j, smiles);
    smiles << pdmd.symbol_for_atom(j);
    smiles << ')';
  }

  return;
}

/*
  In creating a smiles, we have encountered a ring forming bond between items ifrom and
  ito. If ito already has a ring number assigned, then we use that number. Otherwise
  we assign a new ring number to ourselves - and then when the other end is traversed,
  our ring number will be used
*/

/*void
Path_so_Far::_deal_with_ring_forming_connection (int ifrom,
                                                 int ito,
                                                 int & next_ring_number_to_assign,
                                                 IWString & smiles)
{
  if (_ring_number[ito])   // ring number already assigned
    smiles << _ring_number[ito];
  else
  {
    _ring_number[ifrom] = next_ring_number_to_assign;
    smiles << next_ring_number_to_assign;
    next_ring_number_to_assign++;
  }

  return;
}*/

int
Path_so_Far::_produce_smiles(const PD_Molecule_Data& pdmd, IWString& smiles)
{
  assert(_indices_in_smiles.number_elements() == _atoms_in_path);

  smiles.resize(64);

  // set_vector (_ring_number, _atoms_in_path, 0);

  atom_number_t previous_atom = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _atoms_in_path; i++) {
    int j = _indices_in_smiles[i];

    atom_number_t a = _path[j];

    if (INVALID_ATOM_NUMBER != previous_atom) {
      pdmd.append_bond_between_atoms(previous_atom, a, smiles);
    }

    _append_atom_and_its_connections(pdmd, j, smiles);

    smiles << _ring_opening_and_closing[j];

    previous_atom = a;
  }

// #define DEBUG_PRODUCE_SMILES
#ifdef DEBUG_PRODUCE_SMILES
  cerr << "Formed '" << smiles << "' from";
  debug_print(cerr);

  Molecule tmp;
  if (!tmp.build_from_smiles(smiles)) {
    cerr << "Invalid smiles '" << smiles << endl;
    abort();
  }
#endif

  return 1;
}

// #define DEBUG_ADD_RING_OPENING_AND_CLOSING_CHARACTERS

void
Path_so_Far::_place_ring_opening_symbols_between_indices(const PD_Molecule_Data& pdmd,
                                                         int ndx1, int ndx2,
                                                         int& ring_number) const
{
#ifdef DEBUG_ADD_RING_OPENING_AND_CLOSING_CHARACTERS
  cerr << "Ring " << ring_number << " between atoms " << _path[ndx1] << " and "
       << _path[ndx2] << endl;
#endif

  pdmd.append_bond_between_atoms(_path[ndx1], _path[ndx2],
                                 _ring_opening_and_closing[ndx1]);

  _ring_opening_and_closing[ndx1] << ring_number;
  _ring_opening_and_closing[ndx2] << ring_number;

  ring_number++;

  return;
}

/*
  Doing the ring openings and closings is complex because at each ring,
  we don't know whether we continued along the path, or took the branch

         7_____6_______5
                       |
                       |
  0---1                |
       \               |
        \              |
         \2____3_______4


  In this case, if we take the path 01765432, we need to put a ring
  opening and closing between items 1 and 2, as well as between 3 and 6
*/

int
Path_so_Far::_add_ring_opening_and_closing_characters(const PD_Molecule_Data& pdmd)
{
  assert(_indices_in_smiles.number_elements() == _atoms_in_path);

  set_vector(_ring_closure_needed, _atoms_in_path * _atoms_in_path, 0);

  for (int i = 0; i < _atoms_in_path; i++) {
    int j = _indices_in_smiles[i];

    assert(j >= 0 && j < _atoms_in_path);

#ifdef DEBUG_ADD_RING_OPENING_AND_CLOSING_CHARACTERS
    cerr << "Item " << j << " atom " << _path[j] << " at position " << i
         << " in the smiles\n";
#endif

    _position_in_smiles[j] = i;
  }

  for (int i = 0; i < _atoms_in_path; i++) {
    const int j =
        _indices_in_smiles[i];  // path member J is the I'th item along the canonical path

    int jprev;
    if (i > 0) {
      jprev = _indices_in_smiles[i - 1];
    } else {
      jprev = -1;
    }

    int jnext;
    if (_atoms_in_path - 1 != i) {
      jnext = _indices_in_smiles[i + 1];
    } else {
      jnext = -1;
    }

#ifdef DEBUG_ADD_RING_OPENING_AND_CLOSING_CHARACTERS
    cerr << "Processing j = " << j << " (atom " << _path[j] << ") ";
    if (jprev >= 0) {
      cerr << "Prev " << jprev << " (atom " << _path[jprev] << ") ";
    }
    if (jnext >= 0) {
      cerr << "Next " << jnext << " (atom " << _path[jnext] << ") ";
    }
    cerr << endl;
#endif

    //  See if the previous and the next item in the PATH are adjacent in the smiles

    if (j > 0) {
      int previous_index_in_path = j - 1;
      if (previous_index_in_path == jprev) {
        ;
      } else if (previous_index_in_path == jnext) {
        ;
      } else {
        _ring_closure_needed[j * _atoms_in_path + previous_index_in_path] = 1;
        _ring_closure_needed[previous_index_in_path * _atoms_in_path + j] = 1;
      }
    }

    if (_atoms_in_path - 1 != j) {
      int next_index_in_path = j + 1;
      if (next_index_in_path == jnext) {
        ;
      } else if (next_index_in_path == jprev) {
        ;
      } else {
        _ring_closure_needed[j * _atoms_in_path + next_index_in_path] = 1;
        _ring_closure_needed[next_index_in_path * _atoms_in_path + j] = 1;
      }
    }

    //  Now look for ring closures involving this atom

    int r1 = _index_of_ring_closure1[j];

#ifdef DEBUG_ADD_RING_OPENING_AND_CLOSING_CHARACTERS
    cerr << "Continuing along path from " << j << " (atom " << _path[j] << ") to "
         << jnext << " (atom " << _path[jnext] << "), r1 = " << r1 << "\n";
#endif

    if (r1 < 0) {  // no rings here
      continue;
    }

    if (_position_in_smiles[r1] > (i + 1)) {
      _ring_closure_needed[j * _atoms_in_path + r1] = 1;
      _ring_closure_needed[r1 * _atoms_in_path + j] = 1;
    }

    int r2 = _index_of_ring_closure2[j];

    if (r2 >= 0 && _position_in_smiles[r2] > (i + 1)) {
      _ring_closure_needed[j * _atoms_in_path + r2] = 1;
      _ring_closure_needed[r2 * _atoms_in_path + j] = 1;
    }
  }

  // Now place the ring numbers

  int next_ring_to_open = 1;

  resizable_array<int> ring_closures;

  for (int i = 0; i < _atoms_in_path; i++) {
    int j = _indices_in_smiles[i];

    for (int k = i + 1; k < _atoms_in_path; k++) {
      int l = _indices_in_smiles[k];

      if (_ring_closure_needed[j * _atoms_in_path + l]) {
        ring_closures.add(l);
      }
    }

    int n = ring_closures.number_elements();

    if (0 == n) {
      continue;
    }

    if (1 == n) {
      int r = ring_closures[0];
      _place_ring_opening_symbols_between_indices(pdmd, j, r, next_ring_to_open);
    } else {
      int r1 = ring_closures[0];
      int r2 = ring_closures[1];
      if (_position_in_smiles[r1] > _position_in_smiles[r2]) {
        _place_ring_opening_symbols_between_indices(pdmd, j, r1, next_ring_to_open);
        _place_ring_opening_symbols_between_indices(pdmd, j, r2, next_ring_to_open);
      } else {
        _place_ring_opening_symbols_between_indices(pdmd, j, r2, next_ring_to_open);
        _place_ring_opening_symbols_between_indices(pdmd, j, r1, next_ring_to_open);
      }
    }

    ring_closures.resize_keep_storage(0);
  }

  return 1;
}

/*
  No ring is easy because we will traverse the path either forwards or
  backwards
*/

int
Path_so_Far::_produce_smiles_no_ring(const PD_Molecule_Data& pdmd,
                                     resizable_array<int>& paths_this_molecule)
{
  path_type_no_ring++;

  int direction = _determine_direction_chain(pdmd);

  if (direction > 0) {
    for (int i = 0; i < _atoms_in_path; i++) {
      _indices_in_smiles.add(i);
    }
  } else {
    for (int i = _atoms_in_path - 1; i >= 0; i--) {
      _indices_in_smiles.add(i);
    }
  }

  IWString smiles;

  _produce_smiles(pdmd, smiles);

  return _assign_identifier_and_add_to_paths(smiles, paths_this_molecule);
}

int
Path_so_Far::_determine_direction_chain(const PD_Molecule_Data& pdmd) const
{
  const atomic_number_t* z = pdmd.atomic_number();
  const int* arom = pdmd.aromatic();

  int lndx = 0;
  int rndx = _atoms_in_path - 1;

  while (lndx < rndx) {
    atom_number_t alhs = _path[lndx];
    atom_number_t arhs = _path[rndx];

    //  cerr << "Comparing " << z[alhs] << " and " << z[arhs] << endl;

    if (z[alhs] < z[arhs]) {
      return -1;
    }

    if (z[alhs] > z[arhs]) {
      return 1;
    }

    if (arom[alhs] == arom[arhs]) {
      ;
    } else if (arom[alhs]) {
      return 1;
    } else {
      return -1;
    }

    if (_singly_connected_atoms_to_add[lndx].number_elements() ==
        _singly_connected_atoms_to_add[rndx].number_elements()) {
      ;
    } else if (_singly_connected_atoms_to_add[lndx].number_elements() >
               _singly_connected_atoms_to_add[rndx].number_elements()) {
      return 1;
    } else {
      return -1;
    }

    lndx++;
    rndx--;
  }

  return 1;  // same either direction, just choose one
}

/*
  What is the atom across the ring closure(s) from NDX
*/

int
Path_so_Far::_determine_ring_end(int ndx) const
{
  int r1 = _index_of_ring_closure1[ndx];
  int r2 = _index_of_ring_closure2[ndx];

  if (INVALID_ATOM_NUMBER == r2) {  // just one ring
    return r1;
  }

  // Multiple ring openings at NDX. Determine the one that is furthest away

  int d1;
  if (r1 > ndx) {
    d1 = r1 - ndx;
  } else {
    d1 = ndx - r1;
  }

  int d2;
  if (r2 > ndx) {
    d2 = r2 - ndx;
  } else {
    d2 = ndx - r2;
  }

  if (d1 > d2) {
    return r1;
  }

  return r2;
}

/*
  The most common case is when the first atom closes to the last atom.
  We need special care to deal with things like

  .    ______
  |    |    |
  |____|    |

  (incidentally, things like the above probably won't occur because we
  don't take paths involving bonds in more than 1 ring)

  and spiro fusions like

    /|
  |/

  which are also considered to be completely ring systems.  We are
  trying to differentiate these from things like

  .    _____________________    .
  |    |                   |    |
  |____|                   |____|
*/

int
Path_so_Far::_path_is_a_whole_ring_or_ring_system(const PD_Molecule_Data& pdmd) const
{
  assert(INVALID_ATOM_NUMBER != _index_of_ring_closure1[0]);

  // the most common case of a whole ring system
  if (_atoms_in_path - 1 == _index_of_ring_closure1[0]) {
    return 1;
  }

  int end_of_furthest_ring = _index_of_ring_closure1[0];
  if (_index_of_ring_closure2[0] > end_of_furthest_ring) {
    end_of_furthest_ring = _index_of_ring_closure2[0];
  }

  for (int i = 1; i < _atoms_in_path; i++) {
    //  cerr << "AT i = " << i << " end_of_furthest_ring " << end_of_furthest_ring <<
    //  endl;
    if (_index_of_ring_closure1[i] < 0) {
      if (i > end_of_furthest_ring) {  // have encountered a chain section
        return 0;
      }
    } else {
      if (_index_of_ring_closure1[i] > end_of_furthest_ring) {
        end_of_furthest_ring = _index_of_ring_closure1[i];
      }
      if (_index_of_ring_closure2[i] > end_of_furthest_ring) {
        end_of_furthest_ring = _index_of_ring_closure2[i];
      }

      if (_atoms_in_path - 1 == end_of_furthest_ring) {
        return 1;
      }
    }
  }

  cerr << "Warning, why is the code going to line " << __LINE__ << endl;
  debug_print(cerr);

  return 0;
}

/*
  If either end of the path is non-ring, we will start at that end.

  Then we deal with the case of a ring at both ends

  C1CCCCC1CCC1CCCCC1
*/

int
Path_so_Far::_produce_smiles_with_rings(const PD_Molecule_Data& pdmd,
                                        resizable_array<int>& paths_this_molecule)
{
  if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[0]) {
    return _produce_smiles_from_chain_to_rings(pdmd, 1, paths_this_molecule);
  } else if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[_atoms_in_path - 1]) {
    return _produce_smiles_from_chain_to_rings(pdmd, -1, paths_this_molecule);
  }

  // If every atom in the path is part of a ring, don't process - too expensive to
  // canonicalise

#ifdef ECHO_RING_INFO
  if (_path_is_a_whole_ring_or_ring_system(pdmd)) {
    cerr << "Path is whole ring\n";
  } else {
    cerr << "Not whole ring, continuing\n";
  }
#endif

  if (_path_is_a_whole_ring_or_ring_system(pdmd)) {
    return 0;
  }

  if (_atoms_in_path - 1 == _index_of_ring_closure1[0]) {
    return 0;
  }

  if (_atoms_in_path - 1 == _index_of_ring_closure2[0]) {
    return 0;
  }

  // Now this gets gruesome, we have rings at both ends. First determine the general
  // direction of traversal.

  int index_of_end_of_first_ring = _determine_ring_end(0);
  int index_of_start_of_last_ring = _determine_ring_end(_atoms_in_path - 1);

  int too_complicated = 0;

  int d = _determine_direction_rings_at_both_ends(pdmd, too_complicated);

  if (too_complicated) {
    return 0;
  }

  // cerr << "_produce_smiles_with_rings:direction " << d << endl;

  path_type_rings_both_ends++;

  if (d > 0) {
    return _produce_smiles_with_rings_left_to_right(pdmd, index_of_end_of_first_ring,
                                                    paths_this_molecule);
  } else {
    return _produce_smiles_with_rings_right_to_left(pdmd, index_of_start_of_last_ring,
                                                    paths_this_molecule);
  }
}

/*
  Our path starts with a chain atom. The end of the chain will be a closed ring
*/

int
Path_so_Far::_produce_smiles_from_chain_to_rings(
    const PD_Molecule_Data& pdmd, int direction,
    resizable_array<int>& paths_this_molecule)
{
  path_type_chain_to_rings++;

  // cerr << "Chain to rings, direction " << direction << endl;
  // debug_print(cerr);

  int ring_start = -1;

  if (direction > 0) {
    for (int i = 0; i < _atoms_in_path; i++) {
      _indices_in_smiles.add(i);
      if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[i]) {
        ring_start = i;
        break;
      }
    }
  } else {
    for (int i = _atoms_in_path - 1; i >= 0; i--) {
      _indices_in_smiles.add(i);
      if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[i]) {
        ring_start = i;
        break;
      }
    }
  }

  assert(ring_start >= 0);

  int d = _determine_direction_around_ring(pdmd, ring_start);

  // cerr << "At atom " << ring_start << " ring direction " << d << endl;

  if (d > 0)  // continue along the path
  {
    if (direction > 0) {
      for (int i = ring_start + 1; i < _atoms_in_path; i++) {
        _indices_in_smiles.add(i);
      }
    } else {
      for (int i = ring_start - 1; i >= 0; i--) {
        _indices_in_smiles.add(i);
      }
    }
  } else  // take the branch
  {
    if (direction > 0) {
      for (int i = _atoms_in_path - 1; i > ring_start; i--) {
        _indices_in_smiles.add(i);
      }
    } else {
      for (int i = 0; i < ring_start; i++) {
        _indices_in_smiles.add(i);
      }
    }
  }

  _add_ring_opening_and_closing_characters(pdmd);

  IWString smiles;

  _produce_smiles(pdmd, smiles);

  return _assign_identifier_and_add_to_paths(smiles, paths_this_molecule);
}

int
Path_so_Far::_produce_smiles_with_rings_left_to_right(
    const PD_Molecule_Data& pdmd, int index_of_end_of_first_ring,
    resizable_array<int>& paths_this_molecule)
{
  assert(index_of_end_of_first_ring > 0);

  int d = _determine_direction_around_ring(pdmd, index_of_end_of_first_ring);

  // cerr << "left_to_right:index_of_end_of_first_ring " << index_of_end_of_first_ring <<
  // ", atom " << _path[index_of_end_of_first_ring] << ", direction " << d << endl;

  if (d > 0) {
    for (int i = 0; i < index_of_end_of_first_ring; i++) {
      _indices_in_smiles.add(i);
    }
  } else {
    for (int i = index_of_end_of_first_ring - 1; i >= 0; i--) {
      _indices_in_smiles.add(i);
    }
  }

  _indices_in_smiles.add(index_of_end_of_first_ring);

  int ndx = index_of_end_of_first_ring + 1;

  while (INVALID_ATOM_NUMBER == _index_of_ring_closure1[ndx]) {
    _indices_in_smiles.add(ndx);
    ndx++;
  }

  d = _determine_direction_around_ring(pdmd, ndx);

  // cerr << "2nd ring starts at " << ndx << " atom " << _path[ndx] << " direction " << d
  // << endl;

  if (d > 0) {
    for (int i = ndx; i < _atoms_in_path; i++) {
      _indices_in_smiles.add(i);
    }
  } else {
    _indices_in_smiles.add(ndx);
    for (int i = _atoms_in_path - 1; i > ndx; i--) {
      _indices_in_smiles.add(i);
    }
  }

#ifdef ECHO_SMILES_POSITION
  for (int i = 0; i < _atoms_in_path; i++) {
    int j = _indices_in_smiles[i];

    cerr << "Smiles position " << i << " index " << j << " atom " << _path[j] << endl;
  }
#endif

  _add_ring_opening_and_closing_characters(pdmd);

  IWString smiles;
  _produce_smiles(pdmd, smiles);

  return _assign_identifier_and_add_to_paths(smiles, paths_this_molecule);
}

int
Path_so_Far::_produce_smiles_with_rings_right_to_left(
    const PD_Molecule_Data& pdmd, int index_of_start_of_last_ring,
    resizable_array<int>& paths_this_molecule)
{
  assert(index_of_start_of_last_ring > 0);

  int d = _determine_direction_around_ring(pdmd, index_of_start_of_last_ring);

  // cerr << "right to left, at " << index_of_start_of_last_ring << " atom " <<
  // _path[index_of_start_of_last_ring] << " direction " << d << endl;

  if (d > 0) {
    for (int i = index_of_start_of_last_ring + 1; i < _atoms_in_path; i++) {
      _indices_in_smiles.add(i);
    }
  } else {
    for (int i = _atoms_in_path - 1; i > index_of_start_of_last_ring; i--) {
      _indices_in_smiles.add(i);
    }
  }

  _indices_in_smiles.add(index_of_start_of_last_ring);

  int ndx = index_of_start_of_last_ring - 1;

  while (INVALID_ATOM_NUMBER == _index_of_ring_closure1[ndx]) {
    _indices_in_smiles.add(ndx);
    ndx--;
  }

  d = _determine_direction_around_ring(pdmd, ndx);

  // cerr << "Second ring " << ndx << " atom " << _path[ndx] << " direction " << d <<
  // endl;

  if (d > 0) {
    for (int i = ndx; i >= 0; i--) {
      _indices_in_smiles.add(i);
    }
  } else {
    _indices_in_smiles.add(ndx);

    for (int i = 0; i < ndx; i++) {
      _indices_in_smiles.add(i);
    }
  }

  _add_ring_opening_and_closing_characters(pdmd);

  IWString smiles;

  _produce_smiles(pdmd, smiles);

  return _assign_identifier_and_add_to_paths(smiles, paths_this_molecule);
}

/*
  Someone needs to build the smiles for a ring (or ring system)
  starting with atom NDX.  Which way should they go?  Returning 1
  means continue contiguously along the path, -1 means along the ring
  closure bond.

  We set up paths along the two sides of the ring and look for
  differences.
*/

int
Path_so_Far::_determine_direction_around_ring(const PD_Molecule_Data& pdmd, int ndx) const
{
// #define DEBUG_DETERMINE_DIRECTION_AROUND_RING
#ifdef DEBUG_DETERMINE_DIRECTION_AROUND_RING
  cerr << "NDX " << ndx << endl;
  debug_print(cerr);
#endif

  // The assertion is not necessary, C1[C@@]2([C@H]1[C@@H](C)C(=O)C2)C(C)C MFCD00001313
  // assert (INVALID_ATOM_NUMBER == _index_of_ring_closure2[ndx]);

  int ring_closure_bond_invariant;
  int lhs, rhs;
  int lhs_wins, rhs_wins;
  int easiest;

  if (_index_of_ring_closure1[ndx] < ndx)  // ring closure further back in path
  {
    lhs = _index_of_ring_closure1[ndx];
    rhs = ndx - 1;
    ring_closure_bond_invariant = pdmd.canonical_bond_number(_path[ndx], _path[lhs]);
    lhs_wins = -1;  // take branch
    rhs_wins = 1;   // continue along path
    easiest = rhs_wins;
  } else {
    lhs = ndx + 1;
    rhs = _index_of_ring_closure1[ndx];
    //  cerr << "Bond between " << _path[ndx] << " and " << _path[rhs] << endl;
    ring_closure_bond_invariant = pdmd.canonical_bond_number(_path[ndx], _path[rhs]);
    lhs_wins = 1;
    rhs_wins = -1;
    easiest = lhs_wins;
  }

#ifdef DEBUG_DETERMINE_DIRECTION_AROUND_RING
  cerr << "LHS " << lhs << endl;
  cerr << "Withing ring compare " << _bond[lhs] << " with " << ring_closure_bond_invariant
       << endl;
#endif

  if (_bond[lhs] == ring_closure_bond_invariant) {
    ;
  } else if (_bond[lhs] < ring_closure_bond_invariant) {
    return lhs_wins;
  } else {
    return rhs_wins;
  }

  const int* atom_invariant = pdmd.atom_invariant();

  assert(lhs < rhs);

  int steps_taken = 0;

  while (1) {
#ifdef DEBUG_DETERMINE_DIRECTION_AROUND_RING
    cerr << "Advance consideration to " << lhs << " atom " << _path[lhs] << " and " << rhs
         << " atom " << _path[rhs] << endl;
    cerr << "Compare atom invariants " << atom_invariant[_path[lhs]] << " and "
         << atom_invariant[_path[rhs]] << endl;
#endif

    if (atom_invariant[_path[lhs]] == atom_invariant[_path[rhs]]) {
      ;
    } else if (atom_invariant[_path[lhs]] < atom_invariant[_path[rhs]]) {
      return rhs_wins;
    } else {
      return lhs_wins;
    }

    //   If there is a ring closure on just one atom, then it wins. Make sure it isn't a
    //   ring closure back to NDX.

    if (_index_of_ring_closure1[lhs] ==
        _index_of_ring_closure1[rhs]) {  // most likely neither set
      ;
    } else  // may have as many as two ring closures on each atom
    {
      int lhs1 = _index_of_ring_closure1[lhs];
      int lhs2 = _index_of_ring_closure2[lhs];
      int rhs1 = _index_of_ring_closure1[rhs];
      int rhs2 = _index_of_ring_closure2[rhs];

      if (rhs2 == ndx) {
        rhs2 = -1;
      } else if (rhs1 == ndx) {
        rhs1 = rhs2;
        rhs2 = -1;
      }
      if (lhs2 == ndx) {
        lhs2 = -1;
      } else if (lhs1 == ndx) {
        lhs1 = lhs2;
        lhs2 = -1;
      }

      if (lhs1 < 0 && rhs1 < 0) {  // no extra ring closures here
        ;
      } else if (lhs2 >= 0 && rhs2 < 0) {
        return lhs_wins;
      } else if (lhs2 < 0 && rhs2 >= 0) {
        return rhs_wins;
      } else if (lhs1 >= 0) {
        return lhs_wins;
      } else if (rhs1 >= 0) {
        return rhs_wins;
      }
    }

#ifdef DEBUG_DETERMINE_DIRECTION_AROUND_RING
    cerr << "Compare bonds " << _bond[lhs] << " and " << _bond[rhs - 1] << endl;
#endif

    if (_bond[lhs] == _bond[rhs - 1]) {
      ;
    } else if (_bond[lhs] < _bond[rhs - 1]) {
      return rhs_wins;
    } else {
      return lhs_wins;
    }

    lhs++;
    rhs--;

    if (lhs >= rhs) {
      return easiest;
    }

    steps_taken++;
  }
}

/*
  The array A will be one of _index_of_ring_closure1 or _index_of_ring_closure2
  AVOID will be one of the ring openings for the terminal rings
*/

int
Path_so_Far::_find_corresponding_ring_join(int ndx, const int* a, int avoid) const
{
  if (a[ndx] < 0) {
    return -1;
  }

  if (a[ndx] == avoid) {
    return -1;
  }

  return a[ndx];
}

/*
  During canonicalisation of dumb-bell shaped systems, we need to
  simultaneously compare two items from each side = total of 4 items
*/

#define FOUR_WAY_COMPARISON(Mlhs1, Mlhs2, Mrhs1, Mrhs2) \
  if (Mlhs1 == Mrhs1 && Mlhs2 == Mrhs2)                 \
    ;                                                   \
  else if (Mlhs1 == Mrhs2 && Mlhs2 == Mrhs1)            \
    ;                                                   \
  else if ((Mlhs1 + Mlhs2) < (Mrhs1 + Mrhs2))           \
    return 1;                                           \
  else if ((Mlhs1 + Mlhs2) > (Mrhs1 + Mrhs2))           \
    return -1;

// #define DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS

/*
  We have a path with rings at either end. Some cases will be too complicated
  to do inexpensively.
*/

int
Path_so_Far::_determine_direction_rings_at_both_ends(const PD_Molecule_Data& pdmd,
                                                     int& too_complicated) const
{
  too_complicated = 0;

  int aip1 = _atoms_in_path - 1;

  int lhs_ring_end = _index_of_ring_closure1[0];
  if (_index_of_ring_closure2[0] > lhs_ring_end) {
    lhs_ring_end = _index_of_ring_closure2[0];
  }

  int rhs_ring_start = _index_of_ring_closure1[aip1];
  if (_index_of_ring_closure2[aip1] > 0 &&
      _index_of_ring_closure2[aip1] > rhs_ring_start) {
    rhs_ring_start = _index_of_ring_closure2[aip1];
  }

  // First compare ring sizes

  int lhs_ring_size = lhs_ring_end + 1;
  int rhs_ring_size = _atoms_in_path - rhs_ring_start;

#ifdef DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS
  cerr << "Compare ring sizes " << lhs_ring_size << " and " << rhs_ring_size << endl;
#endif

  if (lhs_ring_size == rhs_ring_size) {
    ;
  } else if (lhs_ring_size < rhs_ring_size) {
    return 1;
  } else {
    return -1;
  }

  // The middle of the collection is determined by whether we have an odd or even number
  // of atoms

  int lhs, rhs;

  if (0 == _atoms_in_path % 2) {
    rhs = _atoms_in_path / 2;
    lhs = rhs - 1;
  } else  // when we have a centre atom, check the bonds on either side
  {
    rhs = _atoms_in_path / 2 + 1;
    lhs = rhs - 2;

    if (_bond[lhs] < _bond[rhs - 1]) {
      return 1;
    }
    if (_bond[lhs] > _bond[rhs - 1]) {
      return 1;
    }
  }

  const int* atom_invariant = pdmd.atom_invariant();

  // The cheapest test is to just sum the atomic invariants on either side

  int sumlhs = 0;
  for (int i = 0; i <= lhs; i++) {
    sumlhs += atom_invariant[_path[i]];
  }

  int sumrhs = 0;
  for (int i = rhs; i < _atoms_in_path; i++) {
    sumrhs += atom_invariant[_path[i]];
  }

#ifdef DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS
  cerr << "lhs " << lhs << " sumlhs " << sumlhs << " rhs " << rhs << " sumrhs " << sumrhs
       << endl;
#endif

  if (sumlhs == sumrhs) {
    ;
  } else if (sumlhs < sumrhs) {
    return -1;
  } else if (sumlhs > sumrhs) {
    return 1;
  }

  // Now we need to scan from the central atom out until we get to the
  // left and right rings.  The ring opening checks are so we can deal
  // with things like C1CCC1C1CC1CNC1CCC1C1CCC1 and
  // C1CCC1CC1CC1NC1CCC1C1CCC1 which might otherwise appear identical.
  // How about C1CCC1C1C2CC12NC1CCC1C1CCC1

#ifdef DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS
  cerr << "Start loop, lhs " << lhs << " lhs_ring_end " << lhs_ring_end
       << " rhs = " << rhs << endl;
#endif

  for (; lhs > lhs_ring_end; lhs--, rhs++) {
#ifdef DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS
    cerr << "lhs " << lhs << " (atom " << _path[lhs] << ") and rhs " << rhs << " (atom "
         << _path[rhs] << ") atom invariants " << atom_invariant[_path[lhs]] << " and "
         << atom_invariant[_path[rhs]] << " bond " << _bond[lhs - 1] << " and "
         << _bond[rhs] << endl;
#endif

    if (atom_invariant[_path[lhs]] == atom_invariant[_path[rhs]]) {
      ;
    } else if (atom_invariant[_path[lhs]] > atom_invariant[_path[rhs]]) {
      return 1;
    } else {
      return -1;
    }

    if (_bond[lhs - 1] == _bond[rhs]) {
      ;
    } else if (_bond[lhs - 1] < _bond[rhs]) {
      return 1;
    } else if (_bond[lhs - 1] > _bond[rhs]) {
      return -1;
    }

#ifdef DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS
    cerr << "lhs " << lhs << " (atom " << _path[lhs] << ") _index_of_ring_closure1[lhs] "
         << _index_of_ring_closure1[lhs] << "(" << _path[_index_of_ring_closure1[lhs]]
         << ") and " << rhs << "(atom " << _path[rhs] << ") _index_of_ring_closure1[rhs] "
         << _index_of_ring_closure1[rhs] << "(" << _path[_index_of_ring_closure1[rhs]]
         << ")\n";
#endif

    if (_index_of_ring_closure1[lhs] < 0 &&
        _index_of_ring_closure1[rhs] < 0) {  // no rings here
      continue;
    }

    if (_index_of_ring_closure1[lhs] > 0 && _index_of_ring_closure1[rhs] < 0) {
      return 1;
    }

    if (_index_of_ring_closure1[lhs] < 0 && _index_of_ring_closure1[rhs] > 0) {
      return -1;
    }

    //  Ring opening/closings on both atoms. Are the ring sizes differentiated

    int lhs_ring_size1 = lhs - _index_of_ring_closure1[lhs] + 1;
    int rhs_ring_size1 = _index_of_ring_closure1[rhs] - rhs + 1;

    if (lhs_ring_size1 < 0 && rhs_ring_size1 < 0) {  // long path has completed a ring
      return 0;
    }

#ifdef DEBUG_DETERMINE_DIRECTION_RINGS_AT_BOTH_ENDS
    cerr << lhs_ring_size1 << " lhs_ring_size1 " << rhs_ring_size1 << " rhs_ring_size1\n";
#endif
    assert(lhs_ring_size1 > 2);
    assert(rhs_ring_size1 > 2);

    if (_index_of_ring_closure2[lhs] < 0 &&
        _index_of_ring_closure2[rhs] < 0)  // just one ring on each size
    {
      if (lhs_ring_size1 == rhs_ring_size1) {
        ;
      } else if (lhs_ring_size1 < rhs_ring_size1) {
        return 1;
      } else {
        return -1;
      }

      continue;
    }

    //  Now this gets ugly. There are multiple ring openings and closings.

    int lhs_ring_size2 = lhs - _index_of_ring_closure2[lhs];
    int rhs_ring_size2 = _index_of_ring_closure2[rhs] - rhs;

    FOUR_WAY_COMPARISON(lhs_ring_size1, lhs_ring_size2, rhs_ring_size1, rhs_ring_size2)
  }

  if (_index_of_ring_closure2[lhs] < 0 && _index_of_ring_closure2[rhs] < 0) {
    ;
  } else if (_index_of_ring_closure2[lhs] > 0 && _index_of_ring_closure2[rhs] < 0) {
    return 1;
  } else {
    return -1;
  }

  // We have scanned all the atoms in between the two end rings.

  lhs = lhs_ring_end;
  rhs = rhs_ring_start;

  // cerr << "lhs " << lhs << " (atom " << _path[lhs] << ") and rhs " << rhs << " (atom "
  // << _path[rhs] << ")\n";

  // This gets very complex here. We need to keep track of 4 atoms, 2 on each side,
  // simultaneously

  int lhsb1 = pdmd.canonical_bond_number(_path[0], _path[lhs]);
  int lhsb2 = _bond[lhs - 1];
  int rhsb1 = pdmd.canonical_bond_number(_path[rhs], _path[_atoms_in_path - 1]);
  int rhsb2 = _bond[rhs];

  // cerr << "lhsb1 " << lhsb1 << " lhsb2 " << lhsb2 << " rhsb1 " << rhsb1 << " rhsb2 " <<
  // rhsb2 << endl;

  FOUR_WAY_COMPARISON(lhsb1, lhsb2, rhsb1, rhsb2);

  int lhs1 = 0;
  int lhs2 = lhs_ring_end - 1;

  int rhs1 = rhs_ring_start + 1;
  int rhs2 = _atoms_in_path - 1;

  for (; lhs1 < lhs2; lhs1++, lhs2--, rhs1++, rhs2--) {
    //  cerr << "lhs1 " << lhs1 << " (atom " << _path[lhs1] << ") lhs2 " << lhs2 << "
    //  (atom " << _path[lhs2] << ") rhs1 " << rhs1 << " (atom " << _path[rhs1] << ") rhs2
    //  " << rhs2 << " (atom " << _path[rhs2] << ")\n"; cerr << "inv: " <<
    //  atom_invariant[_path[lhs1]] << " and " << atom_invariant[_path[lhs2]] << " and "
    //  << atom_invariant[_path[rhs1]] << " and " << atom_invariant[_path[rhs2]] << endl;
    FOUR_WAY_COMPARISON(atom_invariant[_path[lhs1]], atom_invariant[_path[lhs2]],
                        atom_invariant[_path[rhs1]], atom_invariant[_path[rhs2]]);
    //  cerr << "Not resolved by invariants, bonds " << _bond[lhs1] << ", " << _bond[lhs2]
    //  << ", " << _bond[rhs1] << ", " << _bond[rhs2] << endl;
    if (lhs1 > 0) {
      FOUR_WAY_COMPARISON(_bond[lhs1], _bond[lhs2], _bond[rhs1], _bond[rhs2]);
    }
    //  cerr << "Not resolved by bonds\n";

    //  We may have one or more ring openings here

    int lhsr1 =
        _find_corresponding_ring_join(lhs1, _index_of_ring_closure1, lhs_ring_end);
    int lhsr2 =
        _find_corresponding_ring_join(lhs2, _index_of_ring_closure2, lhs_ring_end);
    int rhsr1 =
        _find_corresponding_ring_join(rhs1, _index_of_ring_closure1, rhs_ring_start);
    int rhsr2 =
        _find_corresponding_ring_join(rhs2, _index_of_ring_closure2, rhs_ring_start);

    //  cerr << "lhsr1 " << lhsr1 << " lhsr2 " << lhsr2 << endl;

    if (lhsr1 < 0 && lhsr2 < 0 && rhsr1 < 0 &&
        rhsr2 < 0) {  // hopefully the most common case
      continue;
    }

    if (lhsr1 < 0 && rhsr1 > 0) {  // rings on just rhs
      return 1;
    }

    if (lhsr1 >= 0 && rhsr1 < 0) {  // rings on just lhs
      return -1;
    }

    if (lhsr2 >= 0 && rhsr2 < 0) {  // more rings on lhs
      return 1;
    }

    if (lhsr2 < 0 && rhsr2 > 0) {  // more rings on rhs
      return -1;
    }

    //  Multiple ring openings. Give up

    path_type_too_complicated++;

    too_complicated = 1;

    return 0;
  }

  // cerr << "Not resolved\n";

  return 0;  // not resolved
}

int
Path_so_Far::_count_ring_openings(int ndx) const
{
  if (INVALID_ATOM_NUMBER == _index_of_ring_closure1[ndx]) {
    return 0;
  }

  if (INVALID_ATOM_NUMBER == _index_of_ring_closure2[ndx]) {
    return 1;
  }

  return 2;
}

/*
  Someone has decided that we will start our ring at position ndx. Which direction
  do we take?. Note that
*/

int
Path_so_Far::_choose_direction_around_ring(int ndx, int start_index_of_ring,
                                           int end_index_of_ring) const
{
  int ring_size = end_index_of_ring - start_index_of_ring;

  assert(ring_size > 1);

  int lhs = ndx;
  int rhs = ndx;

  int nsteps = (ring_size - 1) / 2;

  for (int i = 0; i < nsteps; i++) {
    lhs--;
    if (lhs < start_index_of_ring) {
      lhs = end_index_of_ring;
    }

    rhs++;
    if (rhs > end_index_of_ring) {
      rhs = start_index_of_ring;
    }

    if (_invariant[lhs] < _invariant[rhs]) {
      return -1;
    } else if (_invariant[lhs] > _invariant[rhs]) {
      return 1;
    }
  }

  return 1;  // could not differentiate, return anything
}

int
Path_so_Far::_determine_chain_length(int istart, int istep) const
{
  for (int i = 0; i < _atoms_in_path; i++) {
    if (INVALID_ATOM_NUMBER != _index_of_ring_closure1[istart]) {
      return i;
    }

    istart += istep;
  }

  cerr << "Path_so_Far::_determine_chain_length:did not find ring!\n";

  return 0;
}

/*
  The start and end atoms in the path are not involved with rings.
  If one of the chains is longer, choose it. If the chains are the
  same length, then treat the whole path as non-ring
*/

int
Path_so_Far::_produce_smiles_chain_at_both_ends(const PD_Molecule_Data& pdmd,
                                                resizable_array<int>& paths_this_molecule)
{
  int lhs_chain_length = _determine_chain_length(0, 1);
  int rhs_chain_length = _determine_chain_length(_atoms_in_path - 1, -1);

  if (lhs_chain_length < rhs_chain_length) {
    return _produce_smiles_chain_at_both_ends(pdmd, -1, paths_this_molecule);
  } else if (lhs_chain_length > rhs_chain_length) {
    return _produce_smiles_chain_at_both_ends(pdmd, 1, paths_this_molecule);
  }

  // chains same length. Treat this as a chain.

  int direction = _determine_direction_chain(pdmd);

  return _produce_smiles_chain_at_both_ends(pdmd, direction, paths_this_molecule);
}

/*
  DIRECTION signifies whether we are starting at the end or the beginning
*/

int
Path_so_Far::_produce_smiles_chain_at_both_ends(const PD_Molecule_Data& pdmd,
                                                int direction,
                                                resizable_array<int>& paths_this_molecule)
{
  if (direction > 0) {
    for (int i = 0; i < _atoms_in_path; i++) {
      _indices_in_smiles.add(i);
    }
  } else {
    for (int i = _atoms_in_path - 1; i >= 0; i--) {
      _indices_in_smiles.add(i);
    }
  }

  _add_ring_opening_and_closing_characters(pdmd);

  IWString smiles;

  _produce_smiles(pdmd, smiles);

  return _assign_identifier_and_add_to_paths(smiles, paths_this_molecule);
}

int
Path_so_Far::produce_smiles(const PD_Molecule_Data& pdmd, const int istep,
                            IWString& smiles, int istart, const int istop) const
{
  atom_number_t previous_atom = INVALID_ATOM_NUMBER;

  for (;; istart += istep) {
    atom_number_t a = _path[istart];

    if (INVALID_ATOM_NUMBER != previous_atom) {
      pdmd.append_bond_between_atoms(previous_atom, a, smiles);
    }

    _append_atom_and_its_connections(pdmd, istart, smiles);

    previous_atom = a;

    if (istart == istop) {
      return 1;
    }
  }

  return 1;
}

int
Path_so_Far::produce_smiles(const PD_Molecule_Data& pdmd, int direction,
                            IWString& smiles) const
{
  if (direction > 0) {
    return produce_smiles(pdmd, 1, smiles, 0, _atoms_in_path - 1);
  } else {
    return produce_smiles(pdmd, -1, smiles, _atoms_in_path - 1, 0);
  }
}

int
Path_so_Far::write_matched_atoms(int direction, IWString& buffer) const
{
  int istart, istop, istep;
  if (direction > 0) {
    istart = 0;
    istop = _atoms_in_path;
    istep = 1;
  } else {
    istart = _atoms_in_path - 1;
    istop = -1;
    istep = -1;
  }

  for (int i = istart; i != istop; i += istep) {
    iwdigits.append_number(buffer, _path[i]);
  }

  return 1;
}

class TShaped_Path_Generator
{
 private:
  PD_Molecule_Data& _pdmd;
  Path_so_Far& _top_of_t;
  const int _ndx;

  int _canonical_direction_at_top_of_t;

  //  We don't want to form rings, so we keep a list of all the
  //  atoms both in the top of the T, as well as all those bonded

  int* _to_avoid;

  Path_so_Far _vertical;

  //  private functions

  int
  _compute_top_of_t_canonical_direction();

  int
  _compare_divergent_paths(int direction) const;

  int
  _resolve_equal_length_arms(int& top_left, int& top_right, int& vertical) const;
  int
  _form_tshaped_smiles(resizable_array<int>& paths_this_molecule);

  int
  _fill_avoid_array();

  int
  _forms_ring() const;

  int
  _do_build_tshape_paths(const Bond* b, atom_number_t zatom,
                         resizable_array<int>& paths_this_molecule);
  int
  __do_build_tshape_paths(atom_number_t zatom, resizable_array<int>& paths_this_molecule);

 public:
  TShaped_Path_Generator(PD_Molecule_Data&, Path_so_Far&, int);
  ~TShaped_Path_Generator();

  int
  do_build_tshape_paths(resizable_array<int>& paths_this_molecule);
};

TShaped_Path_Generator::TShaped_Path_Generator(PD_Molecule_Data& pdmd, Path_so_Far& t,
                                               int ndx)
    : _pdmd(pdmd), _top_of_t(t), _ndx(ndx), _vertical(pdmd.natoms())
{
  _canonical_direction_at_top_of_t = _compute_top_of_t_canonical_direction();

  if (0 == _canonical_direction_at_top_of_t) {
    _canonical_direction_at_top_of_t = 1;  // just choose something at random
  }

  _to_avoid = new_int(pdmd.natoms());

  _fill_avoid_array();

  return;
}

TShaped_Path_Generator::~TShaped_Path_Generator()
{
  if (nullptr != _to_avoid) {
    delete[] _to_avoid;
  }

  return;
}

int
TShaped_Path_Generator::_fill_avoid_array()
{
  int n = _top_of_t.atoms_in_path();

  for (int i = 0; i < n; i++) {
    atom_number_t j = _top_of_t.path_member(i);

    _to_avoid[j] = 1;

    const Atom* aj = _pdmd.atomi(j);

    int jcon = aj->ncon();

    for (int k = 0; k < jcon; k++) {
      atom_number_t l = aj->other(j, k);

      _to_avoid[l] = 1;
    }
  }

  return 1;
}

/*
  We can save time by one time computing the preferred direction
  across the top of the T shape.
*/

int
TShaped_Path_Generator::_compute_top_of_t_canonical_direction()
{
  int n = _top_of_t.atoms_in_path();

  int l1 = _ndx;
  int l2 = n - _ndx - 1;

  if (l1 < l2) {
    return -1;
  }

  if (l1 > l2) {
    return 1;
  }

  // Need to start scanning outwards

  int lhs_prev = _top_of_t.path_member(_ndx);
  int rhs_prev = _top_of_t.path_member(_ndx);
  int lhs = _ndx - 1;
  int rhs = _ndx + 1;

  while (lhs >= 0) {
    atom_number_t alhs = _top_of_t.path_member(lhs);
    atom_number_t arhs = _top_of_t.path_member(rhs);

    int blhs = _pdmd.canonical_bond_number(lhs, lhs_prev);
    int brhs = _pdmd.canonical_bond_number(rhs, rhs_prev);

    if (blhs < brhs) {
      return 1;
    }
    if (blhs > brhs) {
      return -1;
    }

    int ailhs = _pdmd.atom_invariant(alhs);
    int airhs = _pdmd.atom_invariant(arhs);

    if (ailhs < airhs) {
      return 1;
    }
    if (ailhs > airhs) {
      return -1;
    }

    lhs_prev = alhs;
    rhs_prev = arhs;
    lhs--;
    rhs++;
  }

  return 0;  // unresolved
}

/*
  We are given 3 numbers, V, LHS and RHS. If they are resolved,
  assign the corresponding parameters to a canonical order
*/

static int
is_resolved(int v, int lhs, int rhs, int& vertical, int& top_left, int& top_right)
{
  if (v == lhs && lhs == rhs) {
    return 0;
  }

  if (v == lhs || v == rhs || lhs == rhs) {
    return 0;
  }

  if (v < lhs && lhs < rhs) {
    vertical = -1;
    top_left = 0;
    top_right = 1;
    return 1;
  }

  if (v < rhs && rhs < lhs) {
    vertical = -1;
    top_right = 0;
    top_left = 1;
    return 1;
  }

  if (rhs < v && v < lhs) {
    top_right = -1;
    vertical = 0;
    top_left = 1;
    return 1;
  }

  if (rhs < lhs && lhs < v) {
    top_right = -1;
    top_left = 0;
    vertical = 1;
    return 1;
  }

  if (lhs < rhs && rhs < v) {
    top_left = -1;
    top_right = 0;
    vertical = 1;
    return 1;
  }

  if (lhs < v && v < rhs) {
    top_left = -1;
    vertical = 0;
    top_right = 1;
    return 1;
  }

  cerr << "Strange resolution problem " << lhs << " " << rhs << " " << vertical << endl;
  return 0;
}

int
TShaped_Path_Generator::_resolve_equal_length_arms(int& top_left, int& top_right,
                                                   int& vertical) const
{
  assert(top_left == top_right);
  assert(vertical == top_right);

  int n = _vertical.atoms_in_path();

  atom_number_t vprev = _vertical.path_member(0);

  assert(vprev == _top_of_t.path_member(_ndx));

  atom_number_t lhs_prev = vprev;
  atom_number_t rhs_prev = vprev;

  for (int i = 1; i < n; i++) {
    atom_number_t v = _vertical.path_member(i);
    int bv = _pdmd.canonical_bond_number(vprev, v);

    atom_number_t lhs = _top_of_t.path_member(_ndx - i);
    atom_number_t rhs = _top_of_t.path_member(_ndx + i);

    int vlhs = _pdmd.canonical_bond_number(lhs_prev, lhs);
    int vrhs = _pdmd.canonical_bond_number(rhs_prev, rhs);

    if (is_resolved(bv, vlhs, vrhs, vertical, top_left, top_right)) {
      return 1;
    }

    bv = _pdmd.atom_invariant(v);
    vlhs = _pdmd.atom_invariant(lhs);
    vrhs = _pdmd.atom_invariant(rhs);

    if (is_resolved(bv, vlhs, vrhs, vertical, top_left, top_right)) {
      return 1;
    }

    vprev = v;
    lhs_prev = lhs;
    rhs_prev = rhs;
  }

  return 0;
}

/*
  We are trying to canonicalise a T shaped entity.  The vertical limb
  is the same length as one of the branches within the top of the T.
  Compare them
*/

int
TShaped_Path_Generator::_compare_divergent_paths(int direction) const
{
  int n = _vertical.atoms_in_path();

  atom_number_t tprev;
  atom_number_t vprev;

  int tndx = _ndx;

  for (int vndx = 0; vndx < n; vndx++) {
    atom_number_t t = _top_of_t.path_member(tndx);
    atom_number_t v = _vertical.path_member(vndx);

    if (vndx > 0) {
      int bt = _pdmd.canonical_bond_number(tprev, t);
      int bv = _pdmd.canonical_bond_number(vprev, v);
      if (bt < bv) {
        return -1;
      }
      if (bt > bv) {
        return 1;
      }
    }

    int ti = _pdmd.atom_invariant(t);
    int vi = _pdmd.atom_invariant(v);

    if (ti > vi) {
      return 1;
    }
    if (ti < vi) {
      return -1;
    }

    tndx += direction;

    tprev = t;
    vprev = v;
  }

  return 0;
}

/*
  What if an atom on the vertical path is also in the top of the path
*/

int
TShaped_Path_Generator::_forms_ring() const
{
  int v = _vertical.atoms_in_path();

  for (int i = 1; i < v; i++)  // first atom in vertical is also in top
  {
    atom_number_t v = _vertical.path_member(i);

    if (_top_of_t.is_in_path(v)) {
      return 1;
    }
  }

  return 0;
}

/*
  We have the top of the T, the position from where the vertical branches off.
  We need to decide on a canonical order for the resulting smiles
*/

int
TShaped_Path_Generator::_form_tshaped_smiles(resizable_array<int>& paths_this_molecule)
{
  int top_left = _ndx + 1;
  int top_right = _top_of_t.atoms_in_path() - _ndx;
  int height = _vertical.atoms_in_path();

// #define DEBUG_FORM_TSHAPED_SMILES
#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "_form_tshaped_smiles: top_left " << top_left << " top_right " << top_right
       << " height " << height << endl;
  for (int i = 0; i < _top_of_t.atoms_in_path(); i++) {
    atom_number_t a = _top_of_t.path_member(i);
    cerr << " T" << i << " atom " << a << endl;
  }
  for (int i = 0; i < _vertical.atoms_in_path(); i++) {
    atom_number_t a = _vertical.path_member(i);
    cerr << " V" << i << " atom " << a << endl;
  }
  cerr << "_forms_ring? " << _forms_ring() << endl;
#endif

  if (discard_tshaped_paths_that_form_rings && _forms_ring()) {
    return 1;
  }

  if (top_left != top_right && top_right != height &&
      top_left != height) {  // all different lengths
    ;
  } else if (top_left == top_right && height > top_left) {
    height = -1;
    if (_canonical_direction_at_top_of_t < 0) {
      top_left--;
    } else {
      top_right++;
    }
  } else if (top_left == top_right && height < top_left) {
    height = -1;
    if (_canonical_direction_at_top_of_t < 0) {
      top_left--;
    } else {
      top_right++;
    }
  } else if (top_left == top_right && height == top_right)  // all equal lengths
  {
    _resolve_equal_length_arms(top_left, top_right,
                               height);  // may not actually resolve them
  } else if (top_left == height) {
    top_right = -1;
    if (_compare_divergent_paths(-1) > 0) {
      top_left--;
    } else {
      height++;
    }
  } else if (top_right == height) {
    top_left = -1;
    if (_compare_divergent_paths(1) > 0) {
      top_right--;
    } else {
      height++;
    }
  } else {
    cerr << "Unforseen circumstance: top_left " << top_left << " top_right " << top_right
         << " height " << height << endl;
    abort();
  }

  // All three measures must be distinct. No, the T may be symmetric

  // assert (top_left != top_right);
  // assert (height != top_right);
  // assert (height != top_left);

  // Now that we have a canonical order, we can produce the smiles

  IWString smiles1, smiles2, smiles3;

  // Because of double bonds, we need to know the atoms around the intersection

  atom_number_t tla = _top_of_t.path_member(_ndx - 1);
  atom_number_t tra = _top_of_t.path_member(_ndx + 1);
  atom_number_t centre = _top_of_t.path_member(_ndx);
  atom_number_t v1 = _vertical.path_member(1);

#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "tla " << tla << " centre " << centre << " tra " << tra << " v1 " << v1 << endl;
#endif

#ifdef ECHO_JOINING_SYMBOLS
  IWString junk;
  _pdmd.append_bond_between_atoms(centre, tla, junk);
  _pdmd.append_bond_between_atoms(centre, tra, junk);
  _pdmd.append_bond_between_atoms(centre, v1, junk);

  cerr << "Atoms " << tla << " centre " << centre << " tra " << tra << " v1 " << v1
       << " joining symbols '" << junk << "'\n";
#endif

  assert(centre == _vertical.path_member(0));

  // cerr << "top_left " << top_left << " top_right " << top_right << " height " << height
  // << endl;

  if (top_left < top_right && top_right < height) {
    _top_of_t.produce_smiles(_pdmd, 1, smiles1, 0, _ndx);
    _pdmd.append_bond_between_atoms(centre, tra, smiles2);
    _top_of_t.produce_smiles(_pdmd, 1, smiles2, _ndx + 1, _top_of_t.atoms_in_path() - 1);
    _pdmd.append_bond_between_atoms(centre, v1, smiles3);
    _vertical.produce_smiles(_pdmd, 1, smiles3, 1, _vertical.atoms_in_path() - 1);
  } else if (top_left < height && height < top_right) {
    _top_of_t.produce_smiles(_pdmd, 1, smiles1, 0, _ndx);
    _pdmd.append_bond_between_atoms(centre, v1, smiles2);
    _vertical.produce_smiles(_pdmd, 1, smiles2, 1, _vertical.atoms_in_path() - 1);
    _pdmd.append_bond_between_atoms(centre, tra, smiles3);
    _top_of_t.produce_smiles(_pdmd, 1, smiles3, _ndx + 1, _top_of_t.atoms_in_path() - 1);
  } else if (top_right < top_left && top_left < height) {
    _top_of_t.produce_smiles(_pdmd, -1, smiles1, _top_of_t.atoms_in_path() - 1, _ndx);
    _pdmd.append_bond_between_atoms(centre, tla, smiles2);
    _top_of_t.produce_smiles(_pdmd, -1, smiles2, _ndx - 1, 0);
    _pdmd.append_bond_between_atoms(centre, v1, smiles3);
    _vertical.produce_smiles(_pdmd, 1, smiles3, 1, _vertical.atoms_in_path() - 1);
  } else if (top_right < height && height < top_left) {
    _top_of_t.produce_smiles(_pdmd, -1, smiles1, _top_of_t.atoms_in_path() - 1, _ndx);
    _pdmd.append_bond_between_atoms(centre, v1, smiles2);
    _vertical.produce_smiles(_pdmd, 1, smiles2, 1, _vertical.atoms_in_path() - 1);
    _pdmd.append_bond_between_atoms(centre, tla, smiles3);
    _top_of_t.produce_smiles(_pdmd, -1, smiles3, _ndx - 1, 0);
  } else if (height < top_left && top_left < top_right) {
    _vertical.produce_smiles(_pdmd, -1, smiles1, _vertical.atoms_in_path() - 1, 0);
    _pdmd.append_bond_between_atoms(centre, tla, smiles2);
    _top_of_t.produce_smiles(_pdmd, -1, smiles2, _ndx - 1, 0);
    _pdmd.append_bond_between_atoms(centre, tra, smiles3);
    _top_of_t.produce_smiles(_pdmd, 1, smiles3, _ndx + 1, _top_of_t.atoms_in_path() - 1);
  } else if (height < top_right && top_right < top_left) {
    _vertical.produce_smiles(_pdmd, -1, smiles1, _vertical.atoms_in_path() - 1, 0);
    _pdmd.append_bond_between_atoms(centre, tra, smiles2);
    _top_of_t.produce_smiles(_pdmd, 1, smiles2, _ndx + 1, _top_of_t.atoms_in_path() - 1);
    _pdmd.append_bond_between_atoms(centre, tla, smiles3);
    _top_of_t.produce_smiles(_pdmd, -1, smiles3, _ndx - 1, 0);
  } else  // not resolved, just do something random
  {
    _vertical.produce_smiles(_pdmd, -1, smiles1, _vertical.atoms_in_path() - 1, 0);
    _pdmd.append_bond_between_atoms(centre, tra, smiles2);
    _top_of_t.produce_smiles(_pdmd, 1, smiles2, _ndx + 1, _top_of_t.atoms_in_path() - 1);
    _pdmd.append_bond_between_atoms(centre, tla, smiles3);
    _top_of_t.produce_smiles(_pdmd, -1, smiles3, _ndx - 1, 0);
    //  cerr << "Not resolved: top_left " << top_left << " top_right " << top_right << "
    //  height " << height << endl; abort();
  }

  IWString smiles;

  smiles << smiles1 << '(' << smiles2 << ')' << smiles3;

#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "T-Shaped smiles is '" << smiles << "'\n";
#endif

#ifdef SUBSTRUCTURE_SEARCH_TSHAPED_PATHS
  Substructure_Query q;
  if (!q.create_from_smarts(smiles)) {
    cerr << "Huh, invalid smarts '" << smiles << "'\n";
    return 0;
  }

  int nhits = q.substructure_search(const_cast<Molecule*>(_pdmd.molecule()));
  if (0 == nhits) {
    cerr << "Bad news, path is not a substructure, '" << smiles << "'\n";
    cerr << "smiles1 " << smiles1 << " smiles2 " << smiles2 << " smiles3 " << smiles3
         << endl;
    cerr << "top_left " << top_left << " top_right " << top_right << " height " << height
         << endl;

    cerr << "NDX " << _ndx << endl;
    _top_of_t.debug_print(cerr);
    _vertical.debug_print(cerr);

    return 0;
  }
#endif

  unsigned int s = assign_unique_identifier(path_to_unique_id, smiles);

  paths_this_molecule.add(s);

  if (stream_for_paths.is_open()) {
    stream_for_paths << smiles << "\n";
  }
  // if (file_descriptor_for_paths > 0)
  //   _write_matched_atoms (smiles, stream_for_paths);

  return 1;
}

int
TShaped_Path_Generator::__do_build_tshape_paths(atom_number_t zatom,
                                                resizable_array<int>& paths_this_molecule)
{
  _form_tshaped_smiles(paths_this_molecule);

  const Atom* a = _pdmd.atomi(zatom);

  int acon = a->ncon();

#ifdef DEDBUG_DO_BUILD_TSHAPE_PATHS
  cerr << "Tshape path continues with atom " << zatom << " with " << acon
       << " connections\n";
#endif

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (_to_avoid[j]) {
      continue;
    }

    if (_vertical.is_in_path(j)) {
      continue;
    }

    _do_build_tshape_paths(b, j, paths_this_molecule);
  }

  return 1;
}

int
TShaped_Path_Generator::_do_build_tshape_paths(const Bond* b, atom_number_t zatom,
                                               resizable_array<int>& paths_this_molecule)
{
#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "In _do_build_tshape_paths, _vertical contains " << _vertical.atoms_in_path()
       << endl;
#endif

  if (_vertical.atoms_in_path() >= max_path_length) {
    return 0;
  }

  if (_pdmd.is_automatically_added(zatom)) {
    return 0;
  }

  if (!_vertical.add(_pdmd, b, zatom)) {
    return 0;
  }

  if (_vertical.rings_present()) {
    _vertical.pop();
    return 1;
  }

#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "Vertical part of tshaped path has added atom " << zatom << endl;
#endif

  int rc = __do_build_tshape_paths(zatom, paths_this_molecule);

  _vertical.pop();

  return rc;
}

int
TShaped_Path_Generator::do_build_tshape_paths(resizable_array<int>& paths_this_molecule)
{
  atom_number_t zatom = _top_of_t.path_member(_ndx);

  const Atom* a = _pdmd.atomi(zatom);

  int acon = a->ncon();

  assert(0 == _vertical.atoms_in_path());

#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "Starting T shaped path with atom " << zatom << ", _top_of_t contains "
       << _top_of_t.atoms_in_path() << " atoms\n";
#endif

  _vertical.add(_pdmd, nullptr, zatom);

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    const atom_number_t j = b->other(zatom);

    if (_top_of_t.is_in_path(j)) {
      continue;
    }

    _do_build_tshape_paths(b, j, paths_this_molecule);
  }

  return 1;
}

static int
do_build_tshape_paths(PD_Molecule_Data& pdmd, Path_so_Far& p,
                      resizable_array<int>& paths_this_molecule)
{
  int n = p.atoms_in_path();

#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "Building t shaped path, n = " << n << " start with "
       << paths_this_molecule.number_elements() << endl;
#endif

  for (int i = n - 2; i >= 1; i--)  // avoid end atoms in order to form T's
  {
    atom_number_t j = p.path_member(i);

    if (2 == pdmd.ncon(j)) {
      continue;
    }

    if (pdmd.doubly_bonded_singly_connected_neighbours(j).number_elements() > 0) {
      continue;
    }

    TShaped_Path_Generator tspg(pdmd, p, i);

    tspg.do_build_tshape_paths(paths_this_molecule);
  }

#ifdef DEBUG_FORM_TSHAPED_SMILES
  cerr << "After generating T shaped paths, have "
       << paths_this_molecule.number_elements() << " paths\n";
#endif

  return 1;
}

class Disubstituted_Ring_Builder
{
 private:
  const Set_of_Atoms& _r;

  int _atoms_in_ring;

  Path_so_Far _lhs_path;

  //  Since we don't know which direction the smiles will be formed, we need to form
  //  the lhs smiles both ways, starting at the end and going towards the ring, and
  //  starting at the ring and going out

  IWString _lhs_smiles_starting_at_ring_going_toward_end;
  IWString _lhs_smiles_starting_at_end_going_toward_ring;

  //  We don't know which direction will be canonicalised, so we need two smiles for the
  //  middle

  IWString _smiles_middle_l2r;
  IWString _smiles_middle_r2l;

  Path_so_Far _rhs_path;

  //  private functions

  int
  _write_matched_atoms(const IWString& smiles, IWString& stream_for_paths) const;
  atom_number_t
  _identify_atom_outside_ring(PD_Molecule_Data& pdmd, atom_number_t zatom) const;
  int
  _build_paths(PD_Molecule_Data& pdmd, int i0, int i1,
               resizable_array<int>& paths_this_molecule);
  int
  _build_from_ring_paths(PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
                         atom_number_t other_side_of_ring,
                         resizable_array<int>& paths_this_molecule);
  int
  __build_from_ring_paths(PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
                          atom_number_t other_side_of_ring,
                          resizable_array<int>& paths_this_molecule);
  int
  _build_rhs_from_ring_paths(PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
                             resizable_array<int>& paths_this_molecule);
  int
  __build_rhs_from_ring_paths(PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
                              resizable_array<int>& paths_this_molecule);
  void
  _place_any_ring_closures(PD_Molecule_Data& pdmd, int ndx, int& next_ring_to_assign,
                           int* ring_closures, char* ring_closure_btype, IWString& s);
  char
  _ring_closure_btype_for_start_to_end_of_ring_join(PD_Molecule_Data& pdmd, int i0,
                                                    int direction) const;
  int
  _form_centre_ring_smiles(PD_Molecule_Data& pdmd, const int i0, const int i1,
                           int direction, IWString& s);
  /*  int _form_ring_smiles_longest_to_left (PD_Molecule_Data & pdmd,
                                    const int i0,
                                    const int i1,
                                    int * ring_closures);
      int _form_ring_smiles_longest_to_right (PD_Molecule_Data & pdmd,
                                    const int i0,
                                    const int i1,
                                    int * ring_closures);*/

 public:
  Disubstituted_Ring_Builder(const Set_of_Atoms&, int matoms);

  int
  build_paths(PD_Molecule_Data&, const resizable_array<int>& indices_of_3_connected_atoms,
              resizable_array<int>& paths_this_molecule);
};

Disubstituted_Ring_Builder::Disubstituted_Ring_Builder(const Set_of_Atoms& r, int matoms)
    : _r(r), _lhs_path(matoms), _rhs_path(matoms)
{
  _atoms_in_ring = _r.number_elements();

  return;
}

int
Disubstituted_Ring_Builder::_write_matched_atoms(const IWString& smiles,
                                                 IWString& stream_for_paths) const
{
  stream_for_paths << smiles;

  _lhs_path.write_matched_atoms(-1, stream_for_paths);

  int ring_size = _r.number_elements();

  for (int i = 0; i < ring_size; i++) {
    iwdigits.append_number(stream_for_paths, _r[i]);
  }

  _rhs_path.write_matched_atoms(1, stream_for_paths);

  stream_for_paths << '\n';

  return 1;
}

static int
compare_paths(PD_Molecule_Data& pdmd, const Path_so_Far& p1, const Path_so_Far& p2)
{
  int n1 = p1.atoms_in_path();
  int n2 = p2.atoms_in_path();

  if (n1 == n2) {
    ;
  } else if (n1 < n2) {
    return -1;
  } else {
    return 1;
  }

  for (int i = 0; i < n1; i++) {
    int i1 = pdmd.atom_invariant(p1.path_member(i));
    int i2 = pdmd.atom_invariant(p2.path_member(i));

    if (i1 == i2) {
      ;
    } else if (i1 < i2) {
      return -1;
    } else {
      return 1;
    }
  }

  n1--;
  for (int i = 0; i < n1; i++) {
    int b1 = p1.bondi(i);
    int b2 = p2.bondi(i);

    if (b1 == b2) {
      ;
    } else if (b1 < b2) {
      return -1;
    } else {
      return 1;
    }
  }

  return 0;  // not resolved
}

/*
  Someone has build the lhs of the ring substituted path, we work on the rhs.
  We only produce a smiles if the LHS is longer than the rhs
*/

int
Disubstituted_Ring_Builder::__build_rhs_from_ring_paths(
    PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
    resizable_array<int>& paths_this_molecule)
{
  IWString smiles_here;

  int d = compare_paths(pdmd, _lhs_path, _rhs_path);

  // cerr << " d = " << d << " smiles possibilities " << _smiles_middle_l2r << " and " <<
  // _smiles_middle_r2l << endl;

  if (d > 0) {
    IWString rhs_smiles_starting_at_ring_going_toward_end;
    _rhs_path.produce_smiles(pdmd, 1, rhs_smiles_starting_at_ring_going_toward_end);

    smiles_here << _lhs_smiles_starting_at_end_going_toward_ring;
    smiles_here << _smiles_middle_l2r;
    smiles_here << rhs_smiles_starting_at_ring_going_toward_end;
    //  cerr << _lhs_smiles_starting_at_end_going_toward_ring << '+' << _smiles_middle_l2r
    //  << '+' << rhs_smiles_starting_at_ring_going_toward_end << " -> " << smiles_here <<
    //  '\n';
  } else if (d < 0) {
    IWString rhs_smiles_starting_at_end_going_toward_ring;
    _rhs_path.produce_smiles(pdmd, -1, rhs_smiles_starting_at_end_going_toward_ring);

    smiles_here << rhs_smiles_starting_at_end_going_toward_ring;
    smiles_here << _smiles_middle_r2l;
    smiles_here << _lhs_smiles_starting_at_ring_going_toward_end;
    //  cerr << _lhs_smiles_starting_at_ring_going_toward_end << '+' << _smiles_middle_r2l
    //  << '+' << rhs_smiles_starting_at_end_going_toward_ring << " -> " << smiles_here <<
    //  '\n';
  }

  if (d)  // only if resolved
  {
    unsigned int s = assign_unique_identifier(path_to_unique_id, smiles_here);

    paths_this_molecule.add(s);

    if (stream_for_paths.active()) {
      _write_matched_atoms(smiles_here, stream_for_paths);
    }
  }

  if (_lhs_path.atoms_in_path() + _rhs_path.atoms_in_path() >= max_path_length) {
    return 1;
  }

  const Atom* a = pdmd.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    if (b->nrings() > 1) {
      continue;
    }

    atom_number_t j = b->other(zatom);

    if (_r.contains(j)) {
      continue;
    }

    _build_rhs_from_ring_paths(pdmd, b, j, paths_this_molecule);
  }

  return 1;
}

int
Disubstituted_Ring_Builder::_build_rhs_from_ring_paths(
    PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
    resizable_array<int>& paths_this_molecule)
{
  if (pdmd.is_automatically_added(zatom)) {
    return 0;
  }

  if (!_rhs_path.add(pdmd, b, zatom)) {
    return 0;
  }

  if (_rhs_path.rings_present()) {
    _rhs_path.pop();
    return 1;
  }

  int rc = __build_rhs_from_ring_paths(pdmd, b, zatom, paths_this_molecule);

  _rhs_path.pop();

  return rc;
}

int
Disubstituted_Ring_Builder::__build_from_ring_paths(
    PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
    atom_number_t other_side_of_ring, resizable_array<int>& paths_this_molecule)
{
  _lhs_smiles_starting_at_ring_going_toward_end.resize_keep_storage(0);
  _lhs_smiles_starting_at_end_going_toward_ring.resize_keep_storage(0);

  _lhs_path.produce_smiles(pdmd, 1, _lhs_smiles_starting_at_ring_going_toward_end);
  _lhs_path.produce_smiles(pdmd, -1, _lhs_smiles_starting_at_end_going_toward_ring);

// #define DEBUG_RING_CENTRE_SMILES
#ifdef DEBUG_RING_CENTRE_SMILES
  cerr << "From ring to end '" << _lhs_smiles_starting_at_ring_going_toward_end << "'\n";
  cerr << "From end to ring '" << _lhs_smiles_starting_at_end_going_toward_ring << "'\n";
#endif

  if (_lhs_path.atoms_in_path() + _rhs_path.atoms_in_path() >= max_path_length) {
    return 1;
  }

  atom_number_t start_atom = _identify_atom_outside_ring(pdmd, other_side_of_ring);

  _build_rhs_from_ring_paths(pdmd, nullptr, start_atom, paths_this_molecule);

  const Atom* a = pdmd.atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    if (b->nrings() > 1) {
      continue;
    }

    const atom_number_t j = a->other(zatom, i);

    if (_r.contains(j)) {
      continue;
    }

    if (_lhs_path.is_in_path(j)) {
      continue;
    }

    _build_from_ring_paths(pdmd, b, j, other_side_of_ring, paths_this_molecule);
  }

  return 1;
}

int
Disubstituted_Ring_Builder::_build_from_ring_paths(
    PD_Molecule_Data& pdmd, const Bond* b, atom_number_t zatom,
    atom_number_t other_side_of_ring, resizable_array<int>& paths_this_molecule)
{
  if (pdmd.is_automatically_added(zatom)) {
    return 0;
  }

  if (_r.contains(zatom)) {
    return 0;
  }

  if (!_lhs_path.add(pdmd, b, zatom)) {
    return 0;
  }

  if (_lhs_path.rings_present()) {
    _lhs_path.pop();
    return 1;
  }

  int rc =
      __build_from_ring_paths(pdmd, b, zatom, other_side_of_ring, paths_this_molecule);

  _lhs_path.pop();

  return rc;
}

/*
  We have formed the smiles for the central part of the paths. Now, build up paths
  on either side.

  Remove the check for valid start atoms. I found that it caused problems in cases
  involving consecutive five membered rings. r5-r5-r5-r5, where it may generate 
  different numbers of paths for different random smiles
*/

int
Disubstituted_Ring_Builder::_build_paths(PD_Molecule_Data& pdmd, int i0, int i1,
                                         resizable_array<int>& paths_this_molecule)
{
  atom_number_t a0 = _r[i0];
  atom_number_t a1 = _r[i1];

  // if (! pdmd.valid_start_atom(a0))
  //   return 0;

  // if (! pdmd.valid_start_atom(a1))
  //   return 0;

  // cerr << "Both atom " << a0 << " and " << a1 << " are valid start atoms\n";

  // Find the first atom outside the ring from a0

  atom_number_t start_atom = _identify_atom_outside_ring(pdmd, a0);

  return _build_from_ring_paths(pdmd, nullptr, start_atom, a1, paths_this_molecule);
}

atom_number_t
Disubstituted_Ring_Builder::_identify_atom_outside_ring(PD_Molecule_Data& pdmd,
                                                        atom_number_t zatom) const
{
  const Atom* a = pdmd.atomi(zatom);

  assert(3 == a->ncon());

  for (int i = 0; i < 3; i++) {
    const Bond* b = a->item(i);

    if (b->nrings()) {
      continue;
    }

    return b->other(zatom);
  }

  cerr << "Disubstituted_Ring_Builder:_identify_atom_outside_ring:no outside ring atom "
          "found\n";
  return INVALID_ATOM_NUMBER;
}

// #define DEBUG_DISCERN_CANONICAL_ORDER

/*
  Our ring attachment points are equidistance from each other regardless
  of which direction around the ring we go. Can we discern a canonical order
*/

static int
discern_canonical_order(PD_Molecule_Data& pdmd, const Set_of_Atoms& zring, int i0, int i1)
{
  int ring_size = zring.number_elements();

  assert(2 * (ring_size / 2) == ring_size);

#ifdef DEBUG_DISCERN_CANONICAL_ORDER
  const Molecule& m = *(pdmd.molecule());
  atom_number_t a0 = zring[i0];
  cerr << "discern_canonical_order:ring has " << ring_size << " atoms, initial "
       << m.smarts_equivalent_for_atom(a0) << endl;
#endif

  int lhs1 = i0;
  int lhs2 = i0;

  int nsteps = ring_size / 2 - 1;

  int prev1 = i0;
  int prev2 = i0;

  for (int i = 0; i < nsteps; i++) {
    lhs1--;
    lhs2++;

    if (lhs1 < 0) {
      lhs1 = ring_size - 1;
    }

    if (lhs2 == ring_size) {
      lhs2 = 0;
    }

    int inv1 = pdmd.atom_invariant(zring[lhs1]);
    int inv2 = pdmd.atom_invariant(zring[lhs2]);

#ifdef DEBUG_DISCERN_CANONICAL_ORDER
    cerr << "lhs1 " << zring[lhs1] << " atom invariant " << inv1 << " lhs2 "
         << zring[lhs2] << " atom invariant " << inv2 << endl;
#endif

    if (inv1 == inv2) {
      ;
    } else if (inv1 < inv2) {
      return -1;
    } else {
      return 1;
    }

    int b1 = pdmd.canonical_bond_number(zring[prev1], zring[lhs1]);
    int b2 = pdmd.canonical_bond_number(zring[prev2], zring[lhs2]);

#ifdef DEBUG_DISCERN_CANONICAL_ORDER
    cerr << "Bond invariants " << b1 << " and " << b2 << endl;
#endif

    if (b1 == b2) {
      ;
    } else if (b1 < b2) {
      return -1;
    } else {
      return 1;
    }

    prev1 = lhs1;
    prev2 = lhs2;
  }

  // We haven't yet checked the last bonds

  lhs1--;
  lhs2++;

  if (lhs1 < 0) {
    lhs1 = ring_size - 1;
  }

  if (lhs2 == ring_size) {
    lhs2 = 0;
  }

  int b1 = pdmd.canonical_bond_number(zring[prev1], zring[lhs1]);
  int b2 = pdmd.canonical_bond_number(zring[prev2], zring[lhs2]);

#ifdef DEBUG_DISCERN_CANONICAL_ORDER
  cerr << "Final " << zring[lhs1] << " and " << zring[lhs2] << " bond invariants " << b1
       << " and " << b2 << endl;
#endif

  if (b1 == b2) {
    return 0;
  } else if (b1 < b2) {
    return -1;
  } else {
    return 1;
  }

  return 0;  // ring not resolved
}

/*static int
all_connections_within_ring_system (PD_Molecule_Data & pdmd,
                                   atom_number_t zatom,
                                   const Set_of_Atoms & zring)
{
  const Atom * a = pdmd.atomi (zatom);

  assert (3 == a->ncon());

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = a->item (i);

    if (0 == b->nrings())    // definitely outside the ring
      return 0;

    atom_number_t j = a->other (zatom, i);

    if (! zring.contains (j))
      return 0;
  }

  return 1;   // all connections were in the ring
}*/

/*
  If our atom is bonded to something not adjacent in ZRING, we
  need to put in a ring opening and closing
*/

void
Disubstituted_Ring_Builder::_place_any_ring_closures(PD_Molecule_Data& pdmd, int ndx,
                                                     int& next_ring_to_assign,
                                                     int* ring_closures,
                                                     char* ring_closure_btype,
                                                     IWString& s)
{
  atom_number_t zatom = _r[ndx];

  const Atom* a = pdmd.atomi(zatom);

  if (3 != a->ncon()) {
    return;
  }

  assert(0 == ring_closures[zatom]);

  const int last_index = _r.number_elements() - 1;

  const atom_number_t last_atom_in_ring = _r.last_item();
  const atom_number_t first_atom_in_ring = _r[0];

  for (int i = 0; i < 3; i++) {
    const Bond* b = a->item(i);

    const atom_number_t j = b->other(zatom);

    if (1 == b->nrings()) {  // passing along the outside of the ring system
      continue;
    }

    if (0 == b->nrings()) {  // outside the ring system
      continue;
    }

    assert(2 == b->nrings());

    if (0 == ndx) {
      if (_r[1] == j || j == last_atom_in_ring) {
        continue;
      }
    } else if (ndx == last_index) {
      if (first_atom_in_ring == j || _r[ndx - 1] == j) {
        continue;
      }
    } else if (ndx > 0 && _r[ndx - 1] == j) {
      continue;
    } else if (ndx < last_index && _r[ndx + 1] == j) {
      continue;
    }

    //  If we encounter ring closure 1, ignore that.

#ifdef DEBUG_QQ
    cerr << "Cross ring jump, " << _r << " from atom " << zatom << " to atom " << j
         << endl;
    cerr << ring_closures[zatom] << " and " << ring_closures[j] << endl;
#endif

    assert(0 == ring_closures[j]);

    ring_closures[zatom] = next_ring_to_assign;
    ring_closures[j] = next_ring_to_assign;

    ring_closure_btype[j] = pdmd.bond_symbol_between_atoms(zatom, j);

    if (' ' != ring_closure_btype[j]) {
      s << ring_closure_btype[j];
    }

    s << next_ring_to_assign;
    next_ring_to_assign++;
  }

  return;
}

/*
  We need to set up the ring closure information for where the last part of the ring
  smiles joins up with the start atom
*/

char
Disubstituted_Ring_Builder::_ring_closure_btype_for_start_to_end_of_ring_join(
    PD_Molecule_Data& pdmd, int i0, int direction) const
{
  atom_number_t zatom = _r[i0];

  const Atom* a = pdmd.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    if (b->is_single_bond() || b->is_aromatic()) {
      continue;
    }

    atom_number_t j = b->other(zatom);

    int ndx = _r.index(j);

    if (ndx < 0) {
      continue;
    }

    //  Don't want to process the atom that is second in the traversal

    if (direction > 0) {
      if (i0 + 1 == ndx) {
        continue;
      }
      if (_r.number_elements() - 1 == i0 && 0 == ndx) {
        continue;
      }
    } else {
      if (i0 - 1 == ndx) {
        continue;
      }
      if (0 == i0 && _r.number_elements() - 1 == ndx) {
        continue;
      }
    }

    return pdmd.bond_symbol_between_atoms(zatom, j);
  }

  return ' ';
}

int
Disubstituted_Ring_Builder::_form_centre_ring_smiles(PD_Molecule_Data& pdmd, const int i0,
                                                     const int i1, int direction,
                                                     IWString& s)
{
  int matoms = pdmd.natoms();

  int* ring_closures = new_int(matoms);
  std::unique_ptr<int[]> free_ring_closures(ring_closures);
  char* ring_closure_btype = new char[matoms];
  std::unique_ptr<char[]> free_ring_closure_btype(ring_closure_btype);
  memset(ring_closure_btype, ' ', matoms);

  s.resize_keep_storage(0);

  set_vector(ring_closures, pdmd.natoms(), 0);

  int next_ring_to_assign = 2;

  atom_number_t previous_atom = INVALID_ATOM_NUMBER;

  // The preferred smiles may take a path other than a double bond at the start of its
  // traversal

  char start_end_join =
      _ring_closure_btype_for_start_to_end_of_ring_join(pdmd, i0, direction);

  // cerr << "start_end_join '" << start_end_join << "'\n";

  int ndx = i0;

  while (1) {
    atom_number_t j = _r[ndx];

    if (INVALID_ATOM_NUMBER != previous_atom) {
      pdmd.append_bond_between_atoms(previous_atom, j, s);
    }
    s << pdmd.symbol_for_atom(j);
    if (i0 == ndx) {
      s << '1';
    }

    if (ring_closures[j]) {
      if (' ' != ring_closure_btype[j]) {
        s << ring_closure_btype[j];
      }
      s << ring_closures[j];
    } else {
      _place_any_ring_closures(pdmd, ndx, next_ring_to_assign, ring_closures,
                               ring_closure_btype, s);
    }

    const Set_of_Atoms& dbscn = pdmd.doubly_bonded_singly_connected_neighbours(j);
    for (int k = 0; k < dbscn.number_elements(); k++) {
      s << '(';
      atom_number_t l = dbscn[k];
      pdmd.append_bond_between_atoms(j, l, s);
      s << pdmd.symbol_for_atom(l);
      s << ')';
    }

    if (i1 == ndx) {
      s << '(';
    }

    previous_atom = j;

    if (direction > 0) {
      ndx++;
      if (ndx == _r.number_elements()) {
        ndx = 0;
      }
    } else {
      ndx--;
      if (ndx < 0) {
        ndx = _r.number_elements() - 1;
      }
    }

    //  cerr << "ndx changed to " << ndx << " atom " << _r[ndx] << endl;

    if (i0 == ndx) {  // back to where we started
      break;
    }
  }

  if (' ' != start_end_join) {
    s << start_end_join;
  }

  s << "1)";

  // cerr << "Ring smiles '" << s << "'\n";

  return 1;
}

/*int
Disubstituted_Ring_Builder::_form_ring_smiles_longest_to_right (PD_Molecule_Data & pdmd,
                                  const int i0,
                                  const int i1,
                                  IWString & s,
                                  int * ring_closures)
{
  int ring_size = _r.number_elements();

  int next_ring_to_assign = 2;

  for (int i = i0; i <= i1; i++)
  {
    atom_number_t j = _r[i];

    s << pdmd.symbol_for_atom (j);
    if (i0 == i)
      s << '1';
    if (ring_closures[j])
      s << ring_closures[j];
    else
      _place_any_ring_closures (pdmd, i, next_ring_to_assign, ring_closures);
  }

  s << '(';
  for (int i = i1 + 1; i < ring_size; i++)
  {
    atom_number_t j = _r[i];
    s << pdmd.symbol_for_atom (j);
    if (ring_closures[j])
      s << ring_closures[j];
    else
      _place_any_ring_closures (pdmd, i, next_ring_to_assign, ring_closures);
  }

  for (int i = 0; i < i0; i++)
  {
    atom_number_t j = _r[i];
    s << pdmd.symbol_for_atom (j);
    if (ring_closures[j])
      s << ring_closures[j];
    else
      _place_any_ring_closures (pdmd, i, next_ring_to_assign, ring_closures);
  }

  s << "1)";

//cerr << "Ring smiles right '" << s << "'\n";

  return 1;
}

int
Disubstituted_Ring_Builder::_form_ring_smiles_longest_to_left (PD_Molecule_Data & pdmd,
                                  const int i0,
                                  const int i1,
                                  IWString & s,
                                  int * ring_closures)
{
  int next_ring_to_assign = 2;

  for (int i = i0; i >= 0; i--)
  {
    atom_number_t j = _r[i];
    s << pdmd.symbol_for_atom (j);
    if (i == i0)
      s << '1';
    if (ring_closures[j])
      s << ring_closures[j];
    else
      _place_any_ring_closures (pdmd, i, next_ring_to_assign, ring_closures);
  }

  for (int i = _r.number_elements() - 1; i >= i1; i--)
  {
    atom_number_t j = _r[i];
    s << pdmd.symbol_for_atom (j);
    if (ring_closures[j])
      s << ring_closures[j];
    else
      _place_any_ring_closures (pdmd, i, next_ring_to_assign, ring_closures);
  }

  s << '(';
  for (int i = i1 - 1; i > i0; i--)
  {
    atom_number_t j = _r[i];
    s << pdmd.symbol_for_atom (j);
    if (ring_closures[j])
      s << ring_closures[j];
    else
      _place_any_ring_closures (pdmd, i, next_ring_to_assign, ring_closures);
  }

  s << "1)";

//cerr << "Ring smiles left '" << s << "'\n";

  return 1;
}*/

/*static void
mark_rings_already_done (PD_Molecule_Data & pdmd,
                         const Set_of_Atoms & whole_system,
                         int * ring_already_done)
{
  Molecule & m = *(const_cast<Molecule *> (pdmd.molecule()));

  int nr = m.nrings();

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    if (ri->members_in_common(whole_system) > 1)
      ring_already_done[i] = 1;
  }

  return;
}*/

// #define DEBUG_ASSEMBLE_WHOLE_SYSTEM

static int
assemble_whole_system(PD_Molecule_Data& pdmd, int whichring, Set_of_Atoms& whole_system,
                      resizable_array<int>& indices_with_connections_outside_ring)
{
  Molecule& m = *(const_cast<Molecule*>(pdmd.molecule()));

  const Ring* r0 = m.ringi(whichring);

#ifdef DEBUG_ASSEMBLE_WHOLE_SYSTEM
  cerr << "Building ring system from ring " << whichring << " " << (*r0) << endl;
#endif

  if (r0->strongly_fused_ring_neighbours()) {
    return 0;
  }

  if (m.rings_with_fused_system_identifier(r0->fused_system_identifier()) > 2) {
    return 0;
  }

  const atom_number_t a0 = r0->item(0);
  const atom_number_t a1 = r0->item(1);

  whole_system.add(a1);

  atom_number_t aprev = a0;
  atom_number_t zatom = a1;

  while (1) {
    const Atom* a = pdmd.atomi(zatom);

    int acon = a->ncon();

#ifdef DEBUG_ASSEMBLE_WHOLE_SYSTEM
    cerr << "aprev " << aprev << " zatom " << zatom << ", ncon " << acon << endl;
#endif

    if (acon > 3) {  // cannot process these
      return 0;
    }

    atom_number_t nextatom = INVALID_ATOM_NUMBER;

    for (int i = 0; i < acon; i++) {
      const Bond* b = a->item(i);

      atom_number_t j = b->other(zatom);

      if (j == aprev) {
        continue;
      }

      if (b->nrings() > 1) {  // don't take paths across the ring
        continue;
      }

      if (j == a1) {  // completed the loop
        continue;
      } else if (1 == b->nrings()) {
        nextatom = j;
      } else if (0 == b->nrings() &&
                 pdmd.doubly_bonded_singly_connected_neighbours(zatom).contains(j)) {
        ;
      } else if (0 == b->nrings()) {
        indices_with_connections_outside_ring.add(whole_system.number_elements() - 1);
      }
    }

    if (INVALID_ATOM_NUMBER == nextatom) {
      break;
    }

    if (a1 == nextatom) {
      break;
    }

    whole_system.add(nextatom);

    aprev = zatom;
    zatom = nextatom;
  }

#ifdef DEBUG_ASSEMBLE_WHOLE_SYSTEM
  cerr << "Found " << indices_with_connections_outside_ring.number_elements()
       << " indices with connections outside the ring\n";
#endif

  if (2 != indices_with_connections_outside_ring.number_elements()) {
    return 0;
  }

#ifdef DEBUG_ASSEMBLE_WHOLE_SYSTEM
  int x0 = indices_with_connections_outside_ring[0];
  int x1 = indices_with_connections_outside_ring[1];

  atom_number_t ax0 = whole_system[x0];
  atom_number_t ax1 = whole_system[x1];

  cerr << "Index " << x0 << " atom " << ax0 << " ncon " << pdmd.atomi(ax0)->ncon() << ", "
       << m.smarts_equivalent_for_atom(ax0) << endl;
  cerr << "Index " << x1 << " atom " << ax1 << " ncon " << pdmd.atomi(ax1)->ncon() << ", "
       << m.smarts_equivalent_for_atom(ax1) << endl;
#endif

  return whole_system.number_elements();
}

static int
identify_indices_of_3_connected_atoms(PD_Molecule_Data& pdmd, const Ring& r,
                                      resizable_array<int>& indices_of_3_connected_atoms)
{
  int n = r.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = r[i];

    const Atom* aj = pdmd.atomi(j);

    int jcon = aj->ncon();

    if (jcon > 3) {  // cannot process these
      return 0;
    }

    if (2 == jcon) {
      continue;
    }

    if (pdmd.doubly_bonded_singly_connected_neighbours(j).number_elements()) {
      continue;
    }

    indices_of_3_connected_atoms.add(i);
  }

  return (2 == indices_of_3_connected_atoms.number_elements());
}

/*
  We can do rings and ring systems under some circumstances.
*/

int
Disubstituted_Ring_Builder::build_paths(
    PD_Molecule_Data& pdmd, const resizable_array<int>& indices_of_3_connected_atoms,
    resizable_array<int>& paths_this_molecule)
{
  assert(2 == indices_of_3_connected_atoms.number_elements());

  int n = _r.number_elements();

  int i0 = indices_of_3_connected_atoms[0];
  int i1 = indices_of_3_connected_atoms[1];

  assert(3 == pdmd.atomi(_r[i0])->ncon());
  assert(3 == pdmd.atomi(_r[i1])->ncon());

  //  We don't bother processing ortho substituted things

  if (i0 + 1 == i1) {  // ortho substitution pattern. Processed elsewhere
    return 0;
  }

  if (0 == i0 && n - 1 == i1) {  // ortho substitution
    return 0;
  }

  // We now have a ring or ring system that has 2, well separated, substitution points

  // cerr << "i0 = " << i0 << " atom " << _r[i0] << " and i1 = " << i1 << " atom " <<
  // _r[i1] << endl;

  // the smiles of the ring will look like 'C1CCC(N1)' which allows substituents
  // to be placed at the beginning and the end

  // We need to ensure that we have the canonical path around the ring, when we enter
  // from i0, and when we enter from i1
  // We can measure distance around the ring two ways. Directly, or off the end

  int d1 = i1 - i0;
  int d2 = n + i0 - i1;

#ifdef DEBUG_RING_CENTRE_SMILES
  cerr << "Forming centre ring smiles, i0 " << i0 << " d1 " << d1 << " and i1 " << i1
       << " d2 " << d2 << endl;
#endif

  if (d1 < d2) {
    _form_centre_ring_smiles(pdmd, i1, i0, -1, _smiles_middle_r2l);
    _form_centre_ring_smiles(pdmd, i0, i1, 1, _smiles_middle_l2r);
  } else if (d1 > d2) {
    _form_centre_ring_smiles(pdmd, i1, i0, 1, _smiles_middle_r2l);
    _form_centre_ring_smiles(pdmd, i0, i1, -1, _smiles_middle_l2r);
  } else if (d1 ==
             d2)  // same distance either way, need canonical direction both directions
  {
    int d = discern_canonical_order(pdmd, _r, i0, i1);
    //  cerr << "d1 = " << d2 << " and d2 = " << d2 << " discern direction at " << d <<
    //  endl;
    if (d > 0) {
      _form_centre_ring_smiles(pdmd, i1, i0, -1, _smiles_middle_r2l);
    } else {
      _form_centre_ring_smiles(pdmd, i1, i0, 1, _smiles_middle_r2l);
    }

    d = discern_canonical_order(pdmd, _r, i1, i0);
    if (d < 0) {
      _form_centre_ring_smiles(pdmd, i0, i1, 1, _smiles_middle_l2r);
    } else {
      _form_centre_ring_smiles(pdmd, i0, i1, -1, _smiles_middle_l2r);
    }
  }

#ifdef DEBUG_RING_CENTRE_SMILES
  cerr << "Ring centre smiles '" << _smiles_middle_l2r << "' and '" << _smiles_middle_r2l
       << "'\n";
#endif

  return _build_paths(pdmd, i0, i1, paths_this_molecule);
}

static int
do_build_disubstituted_ring_paths(PD_Molecule_Data& pdmd,
                                  resizable_array<int>& paths_this_molecule)
{
  int nr = pdmd.nrings();

  if (0 == nr) {
    return 0;
  }

  int matoms = pdmd.natoms();

  int* atom_already_done = new_int(matoms);
  std::unique_ptr<int[]> free_atom_already_done(atom_already_done);

  for (int i = 0; i < nr; i++) {
    const Ring* ri = pdmd.ringi(i);

    //  cerr << "What about ring " << i << " count " <<
    //  ri->count_members_set_in_array(atom_already_done, 1) << endl;

    if (ri->count_members_set_in_array(atom_already_done,
                                       1)) {  // for now, we can't handle spiro fusions
      continue;
    }

    resizable_array<int> indices_of_3_connected_atoms;

    //  cerr << "fsid " << ri->fused_system_identifier() << endl;

    if (!ri->is_fused()) {
      if (!identify_indices_of_3_connected_atoms(pdmd, *ri,
                                                 indices_of_3_connected_atoms)) {
        continue;
      }

      assert(2 == indices_of_3_connected_atoms.number_elements());

      assert(3 == pdmd.atomi((*ri)[indices_of_3_connected_atoms[0]])->ncon());
      assert(3 == pdmd.atomi((*ri)[indices_of_3_connected_atoms[1]])->ncon());

      Disubstituted_Ring_Builder drb(*ri, matoms);
      //    cerr << "Disubstituted ring building building from ring " << i << endl;
      drb.build_paths(pdmd, indices_of_3_connected_atoms, paths_this_molecule);
      ri->set_vector(atom_already_done, 1);
    } else {
      Set_of_Atoms whole_system;

      if (!assemble_whole_system(pdmd, i, whole_system, indices_of_3_connected_atoms)) {
        continue;
      }

      assert(3 == pdmd.atomi(whole_system[indices_of_3_connected_atoms[0]])->ncon());
      assert(3 == pdmd.atomi(whole_system[indices_of_3_connected_atoms[1]])->ncon());

      //    cerr << "Disubstituted ring building,. i = " << i << ", system has " <<
      //    whole_system << " atoms\n";
      Disubstituted_Ring_Builder drb(whole_system, matoms);
      drb.build_paths(pdmd, indices_of_3_connected_atoms, paths_this_molecule);
      whole_system.set_vector(atom_already_done, 1);
    }
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

  if (element_transformations.active()) {
    element_transformations.process(m);
  }

  // If writing descriptors, strip to just first token of name

  if (create_sparse_fingerprint || create_constant_width_fingerprint) {
    ;
  } else {
    const IWString& mname = m.name();
    int i = mname.index(' ');
    if (i > 0) {
      const_IWSubstring tmp(mname);
      tmp.iwtruncate(i);
      m.set_name(tmp);
    }
  }

  return;
}

/*
 */

static int
build_paths(PD_Molecule_Data& pdmd, atom_number_t zatom, const Bond* b, Path_so_Far& p,
            resizable_array<int>& paths_this_molecule)

{
// #define DEBUG_BUILD_PATHS
#ifdef DEBUG_BUILD_PATHS
  cerr << "On entry to build_paths, atom " << zatom << ", path has " << p.atoms_in_path()
       << " items, last atom " << p.last_atom_in_path() << "\n";
  p.debug_print(cerr);
#endif

  if (!p.add(pdmd, b, zatom)) {
    return 1;
  }

#ifdef DEBUG_BUILD_PATHS
  cerr << "Atom " << zatom << " added to path, length " << p.atoms_in_path() << endl;
#endif

  assert(zatom == p.last_atom_in_path());

  if (p.atoms_in_path() < min_path_length) {
    ;
  } else if (p.heteroatoms_in_path(pdmd) < min_heteroatoms_needed_in_path) {
    ;
  } else {
#ifdef DEBUG_BUILD_PATHS
    cerr << "Producing smiles from\n";
    p.debug_print(cerr);
#endif

    p.produce_smiles(pdmd, paths_this_molecule);
  }

#ifdef DEBUG_BUILD_PATHS
  cerr << "Path produced, tshape? " << build_tshape_paths << ", " << p.atoms_in_path()
       << " atoms in path\n";
#endif
  if (!build_tshape_paths) {
    ;
  } else if (p.rings_present()) {
    ;
  } else if (p.atoms_in_path() >= build_tshape_paths) {
    do_build_tshape_paths(pdmd, p, paths_this_molecule);
  }

  if (p.atoms_in_path() == max_path_length) {
#ifdef debug_print
    cerr << "Max path length " << max_path_length << ", popping " << p.last_atom_in_path()
         << endl;
#endif

    p.pop();
    return 1;
  }

#ifdef DEBUG_BUILD_PATHS
  cerr << "building on path of length " << p.atoms_in_path() << ", last atom "
       << p.last_atom_in_path() << endl;
#endif

  const Atom* a = pdmd.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);

    if (b->nrings() >
        1) {  // too hard to canonicalise - only go around the outside of naphthalene
      continue;
    }

    atom_number_t j = b->other(zatom);

#ifdef DEBUG_BUILD_PATHS
    cerr << "From atom " << zatom << " to atom " << j << " in path? " << p.is_in_path(j)
         << endl;
#endif

    if (p.is_in_path(j)) {
      continue;
    }

    //  if (j == zatom)
    //    continue;

    if (pdmd.is_automatically_added(j)) {
      continue;
    }

#ifdef DEBUG_BUILD_PATHS
    cerr << " from atom " << zatom << " to atom " << j << " items in path "
         << p.atoms_in_path() << endl;
#endif

    build_paths(pdmd, j, b, p, paths_this_molecule);

    if (zatom != p.last_atom_in_path()) {
      cerr << "NOT last atom " << zatom << endl;
      p.debug_print(cerr);
    }
    assert(zatom == p.last_atom_in_path());
  }

#ifdef DEBUG_BUILD_PATHS
  cerr << "Done building from atom " << zatom << ", popping " << p.last_atom_in_path()
       << endl;
#endif

  p.pop();

  return 1;
}

static int
form_paths(PD_Molecule_Data& pdmd, resizable_array<int>& paths_this_molecule)
{
  int matoms = pdmd.natoms();

  const int* valid_start_atom = pdmd.valid_start_atom();

  for (int i = 0; i < matoms; i++) {
    if (!valid_start_atom[i]) {
      continue;
    }

    //  cerr << "Atom " << i << " is a valid start point\n";

    Path_so_Far p(matoms);

    build_paths(pdmd, i, nullptr, p, paths_this_molecule);
  }

  if (build_disubstituted_ring_paths) {
    do_build_disubstituted_ring_paths(pdmd, paths_this_molecule);
  }

  return paths_this_molecule.number_elements();
}

static int
write_paths(const resizable_array<int>& p1, const resizable_array<int>& p2,
            std::ostream& output)
{
  resizable_array<int> just_p1, just_p2, common;

  assert(static_cast<unsigned int>(id_to_path.number_elements()) ==
         path_to_unique_id.size());

  IW_STL_Hash_Map_uint s1;

  int n = p1.number_elements();

  for (int i = 0; i < n; i++) {
    int pi = p1[i];

    const IWString& s = *(id_to_path[pi]);

    s1[s]++;

    //  cerr << pi << " Count 1 for '" << s << "' now " << s1[s] << endl;

    if (p2.contains(pi)) {
      common.add(pi);
    } else {
      just_p1.add(pi);
    }
  }

  IW_STL_Hash_Map_uint s2;

  n = p2.number_elements();

  for (int i = 0; i < n; i++) {
    int pi = p2[i];

    const IWString& s = *(id_to_path[pi]);

    s2[s]++;

    //  cerr << pi << " count 2 for '" << s << "' now " << s2[s] << endl;

    if (!p1.contains(pi)) {
      just_p2.add(pi);
    }
  }

  cerr << just_p1.number_elements() << " items just 1, " << common.number_elements()
       << " common " << just_p2.number_elements() << " just p2\n";

  for (IW_STL_Hash_Map_uint::const_iterator i = s1.begin(); i != s1.end(); ++i) {
    const IWString& s = (*i).first;

    unsigned int c1 = (*i).second;

    IW_STL_Hash_Map_uint::const_iterator f = s2.find(s);

    if (f == s2.end()) {  // not found in 2nd variant
      output << s << ' ' << c1 << " unmatched occurrences\n";
    } else if (c1 != (*f).second) {
      output << s << ' ' << c1 << " occurrences, vs " << (*f).second << '\n';
    }
  }

  for (IW_STL_Hash_Map_uint::const_iterator i = s2.begin(); i != s2.end(); ++i) {
    const IWString& s = (*i).first;

    unsigned int c2 = (*i).second;

    IW_STL_Hash_Map_uint::const_iterator f = s1.find(s);

    if (f == s1.end())  // not found in initial smiles
    {
      output << s << ' ' << c2 << " extra occurrences\n";
      continue;
    }
  }
  // write_paths ("common\n", common, output);
  // write_paths ("just P1\n", just_p1, output);
  // write_paths ("just P2\n", just_p2, output);

  return 1;
}

/*class Int_Comparator_Larger
{
  private:
  public:
    int operator() (int i1, int i2) const;
};

int
Int_Comparator_Larger::operator() (int i1, int i2) const
{
  if (i1 < i2)
    return -1;

  if (i1 > i2)
    return 1;

  return 0;
}*/

static void
do_free_global_structures_after_each_molecule()
{
  id_to_path.resize_keep_storage(0);

  path_to_unique_id.clear();

  number_occurrences.resize_keep_storage(0);

  paths_per_molecule.resize_keep_storage(0);

  return;
}

static int
perform_tests(Molecule& m)
{
  if (free_global_structures_after_each_molecule) {
    do_free_global_structures_after_each_molecule();
  }

  if (m.non_sssr_rings()) {
    skipped_testing_non_sssr_ring_molecules++;
    cerr << "Skipping tests on molecule with non-sssr rings " << molecules_read << " '"
         << m.name() << "'\n";
    return 1;
  }

  PD_Molecule_Data pdmd(m);

  resizable_array<int> paths_this_molecule;

  form_paths(pdmd, paths_this_molecule);

  Int_Comparator_Larger int_comparator_larger;

  paths_this_molecule.iwqsort(int_comparator_larger);

  IWString initial_smiles = m.smiles();

  if (verbose > 1) {
    cerr << m.name() << " contains " << paths_this_molecule.number_elements()
         << " paths\n";
  }

  IW_STL_Hash_Set smiles_already_tried;

  for (int i = 0; i < ntest; i++) {
    IWString s = m.random_smiles();

    if (smiles_already_tried.contains(s)) {
      continue;
    }

    smiles_already_tried.insert(s);

    Molecule mi;
    if (!mi.build_from_smiles(s)) {
      cerr << "Yipes, cannot build from smiles '" << s << "'\n";
      return keep_going_after_test_failure;
    }

    if (verbose > 2) {
      cerr << " random_smiles '" << s << "'\n";
    }

    PD_Molecule_Data pdmdi(mi);
    resizable_array<int> paths_this_variant;

    form_paths(pdmdi, paths_this_variant);

    if (verbose > 2) {
      cerr << "Random variant " << i << " produced "
           << paths_this_variant.number_elements() << " paths\n";
    }

    paths_this_variant.iwqsort(int_comparator_larger);

    if (paths_this_molecule == paths_this_variant) {
      continue;
    }

    cerr << "Path construction mismatch, '" << m.name() << "', variant " << (i + 1)
         << '\n';
    cerr << "Initial smiles '" << initial_smiles << "'\n";
    cerr << "Variant smiles '" << s << "'\n";
    write_paths(paths_this_molecule, paths_this_variant, cerr);

    molecules_failing_tests++;

    if (keep_going_after_test_failure) {
      return 1;
    }

    return 0;
  }

  if (verbose > 1) {
    cerr << "Generated " << smiles_already_tried.size() << " random smiles variants for '"
         << m.name() << "'\n";
  }

  return ntest;
}

/*
  The hash function from gdbm
*/

// #define CHECK_HASH_COLLISIONS
#ifdef CHECK_HASH_COLLISIONS
static IW_Hash_Map<unsigned int, IWString> hashed_values;
static int hash_collisions_encountered = 0;
static IW_STL_Hash_Set collision_encountered;
#endif

#define PRD
#ifndef PRD
extern unsigned int
MurmurHash2(const void* key, int len, unsigned int seed);

static unsigned murmurhashseed = 1093234123;
#endif

static unsigned int
compute_hash_function(const IWString& p)
{
  int n = p.length();

  const char* dptr = p.rawchars();

  if (1 == n) {
    return dptr[0];
  }

#ifdef PRD
  unsigned int value = (dptr[0] << 3) * ((dptr[n - 1] << 4) + 2) * (dptr[n / 2] << 2) +
                       n + (dptr[0] - dptr[n - 1]);

  for (int i = 0; i < n; i++) {
    value = dptr[i] + (value << 6) + (value << 16) - value;
  }
#else
  unsigned int value = MurmurHash2(dptr, n, murmurhashseed);
#endif

#ifdef CHECK_HASH_COLLISIONS

  IW_Hash_Map<unsigned int, IWString>::const_iterator f = hashed_values.find(value);
  if (f == hashed_values.end()) {
    hashed_values[value] = p;
  } else if (p == (*f).second) {
    ;
  } else if (collision_encountered.contains(p) &&
             collision_encountered.contains((*f).second)) {
    ;
  } else {
    cerr << "Hash collision, bit " << value << " values '" << (*f).second << "' and '"
         << p << "'\n";
    collision_encountered.insert(p);
    collision_encountered.insert((*f).second);
    hash_collisions_encountered++;
  }
#endif

  /* Return the value. */
  return value;
}

static int
do_create_constant_width_fingerprint(Molecule& m, const resizable_array<int>& p,
                                     IWString& output_buffer)
{
  IW_Bits_Base fp(create_constant_width_fingerprint);

  int n = p.number_elements();

  for (int i = 0; i < n; i++) {
    int j = p[i];

    const IWString& pi = *(id_to_path[j]);

    unsigned int b = compute_hash_function(pi);

    fp.set(b % create_constant_width_fingerprint);
  }

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);

  output_buffer << fingerprint_tag << tmp << ">\n";

  return 1;
}

static int
do_create_sparse_fingerprint(Molecule& m, const resizable_array<int>& p,
                             IWString& output_buffer)
{
  Sparse_Fingerprint_Creator sfpc;

  int n = p.number_elements();

  for (int i = 0; i < n; i++) {
    int j = p[i];

    const IWString& pi = *(id_to_path[j]);

    unsigned int b = compute_hash_function(pi);

    sfpc.hit_bit(b);
  }

  IWString dyascii;

  sfpc.daylight_ascii_form_with_counts_encoded(dyascii);

  output_buffer << fingerprint_tag << dyascii << ">\n";

  return 1;
}

/*
  Main stepping off point for computation. Just identify
  the starting atoms and call other functions
*/

static int
path_descriptors(Molecule& m, Molecule_and_Paths& map)
{
  int matoms = m.natoms();

  PD_Molecule_Data pdmd(m);

  resizable_array<int>& paths_this_molecule = map.paths();

  const int* valid_start_atom = pdmd.valid_start_atom();

  for (int i = 0; i < matoms; i++) {
    if (!valid_start_atom[i]) {
      continue;
    }

    Path_so_Far p(matoms);

    build_paths(pdmd, i, nullptr, p, paths_this_molecule);
  }

  if (build_disubstituted_ring_paths) {
    do_build_disubstituted_ring_paths(pdmd, paths_this_molecule);
  }

  paths_per_molecule[paths_this_molecule.number_elements()]++;

  return 1;
}

static int
path_descriptors(Molecule& m, IWString& output_buffer)
{
  if (create_constant_width_fingerprint || create_sparse_fingerprint) {
    Molecule_and_Paths map(m.name());

    if (!path_descriptors(m, map)) {
      return 0;
    }

    if (!function_as_tdt_filter) {
      output_buffer << smiles_tag << m.smiles() << ">\n";
      output_buffer << identifier_tag << m.name() << ">\n";
    }

    if (create_constant_width_fingerprint) {
      do_create_constant_width_fingerprint(m, map.paths(), output_buffer);
    } else if (create_sparse_fingerprint) {
      do_create_sparse_fingerprint(m, map.paths(), output_buffer);
    }

    if (!function_as_tdt_filter) {
      output_buffer << "|\n";
    }
  } else {
    Molecule_and_Paths* map = new Molecule_and_Paths(m.name());

    if (!path_descriptors(m, *map)) {
      return 0;
    }

    molecule_and_paths.add(map);
  }

  return 1;
}

static unsigned int precomputed_hash_values[64];

static void
initialise_precomputed_hash_values()
{
  IWString tmp;
  tmp.resize(4);

  for (int i = 0; i < 64; ++i) {
    tmp.resize_keep_storage(0);

    tmp << i;

    precomputed_hash_values[i] = compute_hash_function(tmp);
  }

  return;
}

static unsigned int
fetch_precomputed_hash_values(int s)
{
  if (s < 64) {
    return precomputed_hash_values[s];
  }

  IWString tmp;
  tmp << s;

  return compute_hash_function(tmp);
}

static int
write_empty_fingerprint(Molecule& m, IWString_and_File_Descriptor& output)
{
  if (create_constant_width_fingerprint) {
    output << empty_constant_width_fingerprint;
  } else if (create_sparse_fingerprint)  // somewhat arbitrary choices here
  {
    Sparse_Fingerprint_Creator sfpc;

    sfpc.hit_bit(fetch_precomputed_hash_values(m.natoms()));
    sfpc.hit_bit(fetch_precomputed_hash_values(m.nrings()));
    sfpc.hit_bit(fetch_precomputed_hash_values(m.natoms(6)));

    IWString dyascii(32);

    sfpc.daylight_ascii_form_with_counts_encoded(dyascii);

    output << fingerprint_tag << dyascii << ">\n";
  }

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
every_bond_in_a_ring(Molecule& m)
{
  (void)m.ring_membership();

  int ne = m.nedges();

  for (int i = 0; i < ne; i++) {
    const Bond* b = m.bondi(i);

    if (0 == b->nrings()) {
      return 0;
    }
  }

  return 1;
}

static int
path_descriptors_filter(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();

    Molecule m;
    if (!m.build_from_smiles(buffer)) {
      cerr << "Invalid smiles '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }

    preprocess(m);

    if (every_bond_in_a_ring(m)) {
      write_empty_fingerprint(m, output);
      continue;
    }

    path_descriptors(m, output);

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
path_descriptors_filter(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return path_descriptors_filter(input, output);
}

static int
handle_molecule_with_every_bond_in_a_ring(Molecule& m,
                                          IWString_and_File_Descriptor& output)
{
  if (verbose) {
    cerr << "Skipping all-ring molecule " << molecules_read << " '" << m.name() << "'\n";
  }

  if (create_constant_width_fingerprint || create_sparse_fingerprint) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
    write_empty_fingerprint(m, output);
    output << "|\n";
  } else {
    Molecule_and_Paths* map = new Molecule_and_Paths(m.name());

    molecule_and_paths.add(map);
  }

  return 1;
}

static int
_path_descriptors(Molecule& m, IWString_and_File_Descriptor& output)
{
  preprocess(m);

  if (0 == m.natoms()) {
    cerr << "Skipping empty molecule '" << m.name() << "'\n";
    return 0;
  }

  if (every_bond_in_a_ring(m)) {
    return handle_molecule_with_every_bond_in_a_ring(m, output);
  }

  if (ntest) {
    if (perform_tests(m)) {
      return 1;
    } else {
      cerr << "Tests failed on '" << m.name() << "'\n";
      return 0;
    }
  }

  if (stream_for_paths.active()) {
    stream_for_paths << m.smiles() << ' ' << m.name() << '\n';
  }

  if (!path_descriptors(m, output)) {
    return 0;
  }

  if (stream_for_paths.active()) {
    stream_for_paths << "|\n";

    stream_for_paths.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
path_descriptors(data_source_and_type<Molecule>& input,
                 IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (!_path_descriptors(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

static int
path_descriptors(const char* fname, FileType input_type,
                 IWString_and_File_Descriptor& output)
{
  assert(nullptr != fname);

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

  // Want to be stdin capable, so remove this
  // molecule_and_paths.make_room_for_extra_items(input.molecules_remaining());

  molecule_and_paths.resize(5000);

  return path_descriptors(input, output);
}

static int
remove_paths_failing_support_criteria(int min_support_level, int max_support_level)
{
  int n = number_occurrences.number_elements();

  int items_remaining = 0;

  cerr << "Support level " << min_support_level << " to " << max_support_level << endl;
  cerr << "Checking " << n << " paths\n";

  int items_below_support_level = 0;
  int items_above_support_level = 0;

  if (0 == max_support_level) {
    max_support_level = n + 1;
  }

  for (int i = 0; i < n; i++) {
    int c = number_occurrences[i];

    if (c < min_support_level) {
      items_below_support_level++;
      number_occurrences[i] = 0;
    } else if (c > max_support_level) {
      items_above_support_level++;
      number_occurrences[i] = 0;
    } else {
      items_remaining++;
    }
  }

  if (verbose) {
    cerr << items_below_support_level << " items below support level, "
         << items_above_support_level << " above\n";
  }

  return items_remaining;
}

static int
remove_paths_failing_support_criteria(int min_support_level, double min_support_fraction,
                                      int max_support_level, double max_support_fraction)
{
  if (0 == min_support_level &&
      min_support_fraction > 0.0)  // convert everything to counts
  {
    min_support_level =
        static_cast<int>(static_cast<double>(molecules_read) * min_support_fraction);
    if (0 == min_support_level) {
      min_support_level = 1;
    }

    if (verbose) {
      cerr << "minimum support level " << min_support_level << endl;
    }
  }

  if (0 == max_support_level && max_support_fraction > 0.0) {
    max_support_level =
        static_cast<int>(static_cast<double>(molecules_read) * max_support_fraction);
    if (0 == max_support_level) {
      max_support_level = 1;
    }

    if (verbose) {
      cerr << "maximum support level " << max_support_level << endl;
    }
  }

  return remove_paths_failing_support_criteria(min_support_level, max_support_level);
}

static int
get_number_or_percent(Command_Line& cl, char flag, int& as_number, double& as_fraction)
{
  const_IWSubstring x = cl.string_value(flag);

  if (0 == x.length()) {
    cerr << "get_number_or_percent: no '" << flag << "' option present\n";
    return 0;
  }

  if (x.ends_with('%')) {
    x.chop(1);
    if (!x.numeric_value(as_fraction)) {
      return 0;
    }

    if (as_fraction <= 0.0 || as_fraction > 100.0) {
      return 0;
    }

    as_fraction = as_fraction / 100.0;

    return 1;
  }

  // Maybe it is an integer or a fraction

  if (x.numeric_value(as_fraction) && as_fraction < 1.0 && as_fraction > 0.0) {
    return 1;
  }

  return x.numeric_value(as_number);
}

static void
append_translated_name(const IWString& s, const IWString* xref, IWString& output_buffer)
{
  int n = s.length();
  for (int i = 0; i < n; i++) {
    int c = static_cast<int>(s[i]);

    //  cerr << "Character '" << s[i] << "' is index " << static_cast<int>(c) << " which
    //  is '" << xref[c] << "'\n";

    if (xref[c].length() > 0) {
      output_buffer << xref[c];
    } else {
      output_buffer << static_cast<char>(c);
    }
  }

  return;
}

static int
read_character_translation_table_record(const const_IWSubstring& buffer, IWString* xref)
{
  if (2 != buffer.nwords()) {
    cerr << "The translation table must have exactly two tokens\n";
    return 0;
  }

  const_IWSubstring cfrom, cto;

  int i = 0;
  buffer.nextword(cfrom, i);
  buffer.nextword(cto, i);

  if (1 != cfrom.length()) {
    cerr << "The character translation table must be for single characters\n";
    return 0;
  }

  int ndx = static_cast<int>(cfrom[0]);

  xref[ndx] = cto;

  if (verbose > 1) {
    cerr << "Will translate '" << cfrom << "' to '" << xref[ndx] << "'\n";
  }

  return 1;
}

static int
read_character_translation_table(iwstring_data_source& input, IWString* xref)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    buffer.strip_trailing_blanks();

    if (buffer.starts_with("##") || 0 == buffer.length()) {
      continue;
    }

    if (!read_character_translation_table_record(buffer, xref)) {
      cerr << "Invalid character translation record '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_character_translation_table(const char* fname, IWString* xref)

{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open character translation file '" << fname << "'\n";
    return 0;
  }

  return read_character_translation_table(input, xref);
}

static int
setup_empty_sparse_fingerprint()
{
  empty_sparse_fingerprint = fingerprint_tag;
  empty_sparse_fingerprint << ">\n";

  return 1;
}

static int
setup_empty_constant_width_fingerprint(int nb)
{
  IW_Bits_Base fp(nb);

  IWString tmp;

  fp.daylight_ascii_representation_including_nset_info(tmp);

  empty_sparse_fingerprint << fingerprint_tag << tmp << ">\n";

  return 1;
}

static int
path_descriptors(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lm:M:p:P:h:c:C:y:Y:t:s:X:J:fT:D:q");

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
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  set_include_chiral_info_in_smiles(0);
  set_input_aromatic_structures(1);

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

  if (cl.option_present('h')) {
    if (!cl.value('h', min_heteroatoms_needed_in_path) ||
        min_heteroatoms_needed_in_path < 0) {
      cerr << "The minimum number of heteroatoms in a path, -h, must be a whole +ve "
              "number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will only produces paths contining at least "
           << min_heteroatoms_needed_in_path << " heteroatoms\n";
    }
  }

  if (cl.option_present('T')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T')) {
      cerr << "Invalid element transformation specifications (-T)\n";
      usage(4);
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', min_path_length) || min_path_length < 1) {
      cerr << "Invalid minimum path length (-m)\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Minimum path lenght " << min_path_length << endl;
    }
  }

  if (cl.option_present('M')) {
    if (!cl.value('M', max_path_length) || max_path_length < 1) {
      cerr << "Invalid maximum path length (-M)\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Maximum path lenght " << max_path_length << endl;
    }
  }

  if (min_path_length > max_path_length) {
    cerr << "Min/max path length inconsistency, " << min_path_length << '/'
         << max_path_length << endl;
    usage(4);
  }

  if (cl.option_present('y')) {
    int i = 0;
    const_IWSubstring y;
    while (cl.value('y', y, i++)) {
      if ('1' == y) {
        include_hydrogen_info_in_invariant = 1;
      } else if ('c' == y) {
        include_hydrogen_info_in_invariant = 2;
      } else if ('h' == y) {
        only_include_hydrogen_info_on_heteroatoms = 1;
      } else {
        cerr << "Unrecognised -y qualifier '" << y << "'\n";
        usage(4);
      }
    }

    if (only_include_hydrogen_info_on_heteroatoms &&
        0 == include_hydrogen_info_in_invariant) {
      include_hydrogen_info_in_invariant = 1;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('f')) {
    function_as_tdt_filter = 1;
    if (verbose) {
      cerr << "Will work as a TDT filter\n";
    }
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  int min_support_level = 0;
  double min_support_fraction = 0.0;

  if (cl.option_present('c')) {
    if (!get_number_or_percent(cl, 'c', min_support_level, min_support_fraction)) {
      cerr << "Cannot parse lower path count option (-c)\n";
      usage(17);
    }

    if (verbose) {
      cerr << "Lower path count support " << min_support_level << ", and "
           << min_support_fraction << " fraction\n";
    }
  }

  int max_support_level = 0;
  double max_support_fraction = 0.0;

  if (cl.option_present('C')) {
    if (!get_number_or_percent(cl, 'C', max_support_level, max_support_fraction)) {
      cerr << "Cannot parse upper path count option (-C)\n";
      usage(17);
    }

    if (verbose) {
      cerr << "Upper path count support " << max_support_level << ", and "
           << max_support_fraction << " fraction\n";
    }
  }

  if (max_support_level > 0 && min_support_level > max_support_level) {
    cerr << "Min support level " << min_support_level << " inconsistent with max "
         << max_support_level << endl;
    return 5;
  }

  if (max_support_fraction > 0.0 && min_support_fraction > max_support_fraction) {
    cerr << "Min support fraction " << min_support_fraction << " inconsistent with max "
         << max_support_fraction << endl;
    return 5;
  }

  number_occurrences.resize(10000);

  if (cl.option_present('t')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++)) {
      if (t.starts_with("seed=")) {
        t.remove_leading_chars(5);

        int s;
        if (!t.numeric_value(s) || s < 0) {
          cerr << "Invalid random number seed '" << t << "'\n";
          return 4;
        }

        set_smiles_random_number_seed(s);
        if (verbose) {
          cerr << "Random number seed " << s << endl;
        }
      } else if ("free" == t) {
        free_global_structures_after_each_molecule = 1;

        if (verbose) {
          cerr << "Will re-initialise global counters after each molecule\n";
        }
      } else if ("kg" == t) {
        keep_going_after_test_failure = 1;

        if (verbose) {
          cerr << "Will keep going after a test failure\n";
        }
      } else if (t.starts_with("dwf=")) {
        t.remove_leading_chars(4);
        die_when_formed = t;
        if (verbose) {
          cerr << "Will die when smiles '" << die_when_formed << "' is formed\n";
        }
      } else if ("help" == t) {
        cerr << " -t seed=<n>        use <n> as random number seed\n";
        cerr << " -t free            free global arrays after each molecule\n";
        cerr << " -t kg              keep going after a test failure\n";
        cerr << " -t <n>             perform <n> random smiles tests on each molecule\n";
        cerr << " -t dwf=smiles      die if smiles <smiles> is formed\n";

        exit(0);
      } else if (t.numeric_value(ntest) && ntest >= 0) {
        if (verbose) {
          cerr << "Will perform " << ntest << " tests on each molecule\n";
        }
      } else {
      }
    }

    if (0 == ntest) {
      cerr << "Test mode, but no tests to be done\n";
      usage(8);
    }
  }

  IWString fingerprint_cross_reference_file_name;

  if (cl.option_present('X')) {
    int i = 0;
    const_IWSubstring x;
    while (cl.value('X', x, i++)) {
      if ("nodis" == x) {
        build_disubstituted_ring_paths = 0;
        if (verbose) {
          cerr << "Will NOT build disubstituted ring paths\n";
        }
      } else if (x.starts_with("tshape=")) {
        x.remove_leading_chars(7);

        if (!x.numeric_value(build_tshape_paths) || build_tshape_paths < 3) {
          cerr << "The tshape= directive must be followed by a whole number >= 3\n";
          return 4;
        }
        if (verbose) {
          cerr << "Will build T shaped paths from linear paths of length "
               << build_tshape_paths << endl;
        }
      } else if ("oktsr" == x) {
        discard_tshaped_paths_that_form_rings = 0;

        if (verbose) {
          cerr << "Tshaped paths forming rings will be allowed\n";
        }
      } else if (x.starts_with("xref=")) {
        x.remove_leading_chars(5);

        fingerprint_cross_reference_file_name = x;

        if (verbose) {
          cerr << "Fingerprint cross reference file '"
               << fingerprint_cross_reference_file_name << "' will be created\n";
        }
      } else if ("help" == x) {
        cerr << "  -X nodis       do NOT build disubstituted ring paths\n";
        cerr << "  -X tshape=nn   produce tshaped subsets for paths at least NN\n";
        cerr << "  -X oktsr       tshaped paths that form rings will be allowed\n";
        cerr << "                 but there will be no ring info in the path!\n";
        cerr << "  -X help        this message\n";
        return 0;
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        return 3;
      }
    }
  }

  int user_specified_tag = 0;

  if (cl.option_present('J')) {
    cl.value('J', fingerprint_tag);

    if (!fingerprint_tag.ends_with('<')) {
      fingerprint_tag << '<';
    }

    if (fingerprint_tag.starts_with("NC")) {
      create_sparse_fingerprint = 1;
      if (verbose) {
        cerr << "Wil create sparse fingerprint with tag '" << fingerprint_tag << "'\n";
      }

      setup_empty_sparse_fingerprint();
    } else if (fingerprint_tag.starts_with("FP")) {
      create_constant_width_fingerprint = 2048;
      if (verbose) {
        cerr << "Wil create constant width fingerprint with tag '" << fingerprint_tag
             << "'\n";
      }

      setup_empty_constant_width_fingerprint(create_constant_width_fingerprint);
    } else {
      cerr << "Fingerprint tags must start with either NC or FP, '" << fingerprint_tag
           << "' is invalid\n";
      return 4;
    }

    user_specified_tag = 1;
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Cannot set atom typing to apply '" << p << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Atom typing initialised '" << p << "'\n";
    }

    if (!user_specified_tag) {
      fingerprint_tag = "NC";
      atom_typing_specification.append_to_tag(fingerprint_tag);
      fingerprint_tag << '<';
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString descriptor_name_prefix;

  if (cl.option_present('p')) {
    cl.value('p', descriptor_name_prefix);

    if (verbose) {
      cerr << "Descriptors will be prefixed with '" << descriptor_name_prefix << "'\n";
    }

    if (!descriptor_name_prefix.ends_with('_')) {
      descriptor_name_prefix << '_';
    }
  }

  if (cl.option_present('Y')) {
    const char* y = cl.option_value('Y');

    if (!stream_for_paths.open(y)) {
      cerr << "Yipes, cannot open -Y file '" << y << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Paths written to '" << y << "'\n";
    }

    stream_for_paths.resize(9000);
  }

  iwdigits.set_include_leading_space(1);
  iwdigits.initialise(512);  // an unlikely number of atoms

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (function_as_tdt_filter) {
      rc = path_descriptors_filter(cl[i], output);
    } else {
      rc = path_descriptors(cl[i], input_type, output);
    }

    if (0 == rc) {
      rc = i + 1;
      break;
    }
  }

  if (output.size()) {
    output.flush();
  }

  if (fingerprint_cross_reference_file_name.length()) {
    write_path_cross_reference(
        path_to_unique_id, fingerprint_cross_reference_file_name.null_terminated_chars());
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules";
    if (!free_global_structures_after_each_molecule) {
      cerr << ", found " << path_to_unique_id.size() << " paths";
    }
    cerr << "\n";
    cerr << path_type_no_ring << " paths with no rings\n";
    cerr << path_type_chain_to_rings << " paths from a chain to a ring\n";
    cerr << path_type_rings_both_ends << " paths with rings at both ends\n";
    cerr << path_type_too_complicated << " paths that were too complex to process\n";

    for (int i = 0; i < paths_per_molecule.number_elements(); i++) {
      if (paths_per_molecule[i]) {
        cerr << paths_per_molecule[i] << " molecules had " << i << " paths\n";
      }
    }
  }

  // 20171003 IL - reversed logic, rc is 0 only when cl.number_elements() is 0,
  // which is probably when we want to bail out
  if (0 == rc) {
    return rc;
  }

  if (ntest > 0) {
    if (skipped_testing_non_sssr_ring_molecules) {
      cerr << "Skipped " << skipped_testing_non_sssr_ring_molecules
           << " molecules containing non sssr rings\n";
    }

    cerr << molecules_failing_tests << " molecules failed testing\n";

    return rc;
  }

#ifdef CHECK_HASH_COLLISIONS
  cerr << hash_collisions_encountered << " hash collisions encountered, "
       << path_to_unique_id.size() << " unique paths\n";
#endif

  if (create_constant_width_fingerprint || create_sparse_fingerprint) {
    if (cl.option_present('q')) {
      _exit(0);
    }
    return 0;
  }

  IWString* xref = nullptr;

  if (cl.option_present('s')) {
    xref = new IWString[std::numeric_limits<unsigned char>::max()];

    const char* s = cl.option_value('s');

    if (!read_character_translation_table(s, xref)) {
      cerr << "Cannot read character translation table from '" << s << "'\n";
      return 4;
    }
  }

  int paths_being_output = path_to_unique_id.size();

  if (0 != min_support_level || 0 != max_support_level || 0.0 != min_support_fraction ||
      0.0 != max_support_fraction) {
    paths_being_output = remove_paths_failing_support_criteria(
        min_support_level, min_support_fraction, max_support_level, max_support_fraction);
    if (0 == paths_being_output) {
      cerr << "Bad news, no paths met support criteria\n";
      return 34;
    }

    if (verbose) {
      cerr << paths_being_output << " of " << number_occurrences.number_elements()
           << " paths met support criteria\n";
    }
  }

  initialise_precomputed_hash_values();

  unsigned int n = path_to_unique_id.size();

  output.resize(paths_being_output * 30);

  output << "Name";

  for (int i = 0; i < id_to_path.number_elements(); i++) {
    if (0 == number_occurrences[i]) {  // suppressed by support considerations
      continue;
    }

    const IWString& s = *(id_to_path[i]);

    //  cerr << s << " is " << i << endl;

    output << ' ' << descriptor_name_prefix;

    if (fingerprint_cross_reference_file_name.length()) {
      output << i;
    } else if (nullptr == xref) {
      output << s;
    } else {
      append_translated_name(s, xref, output);
    }
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  int* tmp = new int[n];
  std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < molecule_and_paths.number_elements(); i++) {
    const Molecule_and_Paths* mpi = molecule_and_paths[i];

    set_vector(tmp, n, 0);

    mpi->get_path_counts(tmp);

    output << (*mpi);
    for (unsigned int j = 0; j < n; j++) {
      if (0 == number_occurrences[j]) {
        ;
      } else if (tmp[j]) {
        iwdigits.append_number(output, tmp[j]);
      } else {
        output << " 0";
      }
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  if (output.length()) {
    output.flush();
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  if (cl.option_present('q')) {
    _exit(0);
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = path_descriptors(argc, argv);

  return rc;
}
