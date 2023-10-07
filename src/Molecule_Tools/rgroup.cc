/*
  Identifies R groups in a molecule based on a substructure query.
*/

#include <ctype.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;
static int molecules_read = 0;
static int molecules_written = 0;
static int molecules_created = 0;
static int molecules_with_no_hits = 0;

static const Element* dummy_atom_element = nullptr;

static Chemical_Standardisation chemical_standardisation;

static isotope_t isotope_for_join_atom = 0;

static int include_attachment_point_in_all_substituents = 0;

static int min_atoms_in_rgroup = 0;

static int max_atoms_in_rgroup = std::numeric_limits<int>::max();

static int make_implicit_hydrogens_explicit = 0;

/*
  When H atoms are made explicit, we probably want to remove
  most of them before output.
*/

static int remove_explicit_hydrogens_not_attached_to_isotope = 0;

static int do_all_embeddings = 0;

static int combine_substituents_on_same_atom = 1;

static int process_ring_substituents = 0;

//  What do we do when a molecule does not match the query

static int ignore_molecules_which_dont_match_query = 0;

// What do we do when a molecule has more than one match to the query

static int ignore_molecules_with_multiple_hits = 0;

static extending_resizable_array<int>* query_stats = nullptr;
static IW_STL_Hash_Map_int** rgroups_found = nullptr;

static void
usage(int rc = 0) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Usage: " << prog_name << " -q/-s query <options> <input file>\n";
  cerr << "Identifies substituents based on what is attached to query atoms\n";
  cerr << "  -q <file>      specify query(s)\n";
  cerr << "  -s <smarts>    queries as smarts\n";
  cerr << "  -r q.n         find R group from atom N in query Q\n";
  cerr << "  -z i           ignore molecules which don't match a query\n";
  cerr << "  -u             unique embeddings only\n";
  cerr << "  -d <ele>       add an embedding dummy atom, type <ele>, to all groups\n";
  cerr << "  -a             process all embeddings\n";
  cerr << "  -b             process symmetry equivalent embeddings\n";
  cerr << "  -e             do multiple substituents on same atom separately\n";
  cerr << "  -k             include attachment atom in R groups\n";
  cerr << "  -I <iso>       label join points with isotope <iso>\n";
  cerr << "  -X <symbol>    after building, remove all elements of type <symbol>\n";
  cerr << "  -c <natoms>    min atoms in an R group\n";
  cerr << "  -C <natoms>    max atoms in an R group\n";
  cerr << "  -H             make implicit hydrogens explicit\n";
  cerr << "  -h             after processing, remove explicit Hydrogens not attached to isotope (-I)\n";
  cerr << "  -y             process substituents involved in ring bonds\n";
  cerr << "  -i <type>      input type\n";
  display_standard_aromaticity_options (cerr);
  display_standard_smiles_options (cerr);
  cerr << "  -E <symbol>    create element with symbol <symbol>\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  Needed a class for carrying around a bunch of arguments
*/

class RGroup_Construction_Current_State {
 private:
  int _query_number;
  int _embedding_number;
  int _query_atom_number;
  int _parent_molecule_written;

  int* _atom_already_done;

 public:
  RGroup_Construction_Current_State(int);
  ~RGroup_Construction_Current_State();

  int
  query_number() const {
    return _query_number;
  }

  int
  embedding_number() const {
    return _embedding_number;
  }

  int
  query_atom_number() const {
    return _query_atom_number;
  }

  void
  set_query_number(int s) {
    _query_number = s;
  }

  void
  set_embedding_number(int s) {
    _embedding_number = s;
  }

  void
  set_query_atom_number(int s) {
    _query_atom_number = s;
  }

  void
  set_parent_molecule_written(int s) {
    _parent_molecule_written = s;
  }

  int write_parent_molecule_if_needed(Molecule&, IWString_and_File_Descriptor&);

  int do_output(Molecule& r, int ndx, IWString_and_File_Descriptor& output) const;

  int add_dummy_element_and_do_output(Molecule& r, int ndx,
                                      IWString_and_File_Descriptor& output) const;

  int* atom_already_done() {
    return _atom_already_done;
  }
};

RGroup_Construction_Current_State::RGroup_Construction_Current_State(int matoms) {
  _atom_already_done = new_int(matoms);

  _parent_molecule_written = 0;

  return;
}

RGroup_Construction_Current_State::~RGroup_Construction_Current_State() {
  if (nullptr != _atom_already_done) {
    delete[] _atom_already_done;
  }

  return;
}

int
RGroup_Construction_Current_State::write_parent_molecule_if_needed(
    Molecule& m, IWString_and_File_Descriptor& output) {
  if (_parent_molecule_written) {
    return 0;
  }

  output << m.smiles() << ' ' << m.name() << '\n';

  output.write_if_buffer_holds_more_than(32768);

  _parent_molecule_written = 1;

  return 1;
}

/*
  Remove any explicit Hydrogen NOT bonded to an atom with
  isotope ISO
*/

static int
do_remove_explicit_hydrogens_not_attached_to_isotope(Molecule& m, isotope_t iso) {
  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; i--) {
    const Atom* a = m.atomi(i);

    if (1 != a->atomic_number()) {
      continue;
    }

    if (1 != a->ncon()) {
      continue;
    }

    atom_number_t j = a->other(i, 0);

    if (iso == m.isotope(j)) {
      continue;
    }

    m.remove_atom(i);
    rc++;
  }

  return rc;
}

int
RGroup_Construction_Current_State::add_dummy_element_and_do_output(
    Molecule& r, int ndx, IWString_and_File_Descriptor& output) const {
  if (nullptr != dummy_atom_element) {
    Atom* aa = new Atom(dummy_atom_element);
    r.add(aa);
    r.add_bond(0, r.natoms() - 1, SINGLE_BOND);
  }

  if (isotope_for_join_atom) {
    r.set_isotope(0, isotope_for_join_atom);

    if (remove_explicit_hydrogens_not_attached_to_isotope) {
      do_remove_explicit_hydrogens_not_attached_to_isotope(r, isotope_for_join_atom);
    }
  }

  return do_output(r, ndx, output);
}

int
RGroup_Construction_Current_State::do_output(Molecule& r, int ndx,
                                             IWString_and_File_Descriptor& output) const {
  output << r.unique_smiles() << " RG " << _query_number << '.' << _embedding_number
         << '.' << _query_atom_number << '.' << ndx << '\n';

  molecules_written++;

  output.write_if_buffer_holds_more_than(32768);

  if (verbose) {
    query_stats[_query_number][_query_atom_number]++;
    rgroups_found[_query_number][_query_atom_number][r.unique_smiles()]++;
    output.flush();
  }

  return 1;
}

// #define DEBUG_BUILD_MOLECULE

/*
  We need to keep track of atom number cross references between molecule
  M and molecule R.
  To do that, we use the array atom_already_done. For those atoms in R,
  the value of ATOM_ALREADY_DONE[i] will be their atom number (PLUS ONE)
  in the molecule R. Yes, somewhat of a kludge, but it works...
*/

static int
build_molecule(const Molecule& m, atom_number_t parent_atom,
               bond_type_t bt, atom_number_t anchor,
               const Set_of_Atoms& embedding, int* atom_already_done, Molecule& r) {
  assert(0 == atom_already_done[parent_atom]);
  assert(parent_atom != anchor);

  int initial_natoms = r.natoms();

  atom_already_done[parent_atom] = initial_natoms + 1;

  Atom* aa = new Atom(m.atomic_number(parent_atom));
  r.add(aa);

  if (INVALID_ATOM_NUMBER != anchor) {
    r.add_bond(initial_natoms, atom_already_done[anchor] - 1, bt);
  }

  int rc = 1;
  for (const Bond* b : m[parent_atom]) {
    atom_number_t j = b->other(parent_atom);

    if (j == anchor) {
      continue;
    }

#ifdef DEBUG_BUILD_MOLECULE
    cerr << "From atom " << parent_atom << " examine " << j
         << " already_done = " << atom_already_done[j]
         << " contains = " << embedding.contains(j) << '\n';
#endif

    //  Be careful placing bonds, otherwise problems with cyclic structures

    if (atom_already_done[j]) {
      if (atom_already_done[j] > initial_natoms) {
        r.add_bond(initial_natoms, atom_already_done[j] - 1, b->btype());
      }
      continue;
    }

    if (embedding.contains(j)) {
      continue;
    }

    rc += build_molecule(m, j, b->btype(), parent_atom, embedding, atom_already_done, r);
  }

  return rc;
}

/*
  For each query, we need a list of the atoms to be processed.
*/

class Set_of_Query_Atoms : public resizable_array<int> {
 private:
 public:
};

static resizable_array_p<Set_of_Query_Atoms> query_atoms;

template class resizable_array_p<Set_of_Query_Atoms>;
template class resizable_array_base<Set_of_Query_Atoms*>;

static int
size_constraints_met(Molecule& m) {
  int matoms = m.natoms();

  if (matoms < min_atoms_in_rgroup) {
    return 0;
  }

  if (matoms > max_atoms_in_rgroup) {
    return 0;
  }

  return 1;
}

/*
  Build a molecule that includes the attachment point and all
  attachments
*/

static int
rgroup_including_substituents(Molecule& m, RGroup_Construction_Current_State& rgrcs,
                              atom_number_t a, const Set_of_Atoms& atoms_in_substituent,
                              const Set_of_Atoms& embedding,
                              IWString_and_File_Descriptor& output) {
  // cerr << "Processing embedding " << embedding << " from atom " << a << '\n';
  Molecule r;

  if (nullptr != dummy_atom_element) {
    Atom* aa = new Atom(dummy_atom_element);
    r.add(aa);
  } else {
    Atom* aa = new Atom(m.atomic_number(a));
    r.add(aa);
  }

  if (isotope_for_join_atom) {
    r.set_isotope(0, isotope_for_join_atom);
  }

  int* atom_already_done = rgrcs.atom_already_done();

  for (int i = 0; i < atoms_in_substituent.number_elements(); i++) {
    atom_number_t j = atoms_in_substituent[i];

    const Bond* b = m.bond_between_atoms(a, j);

    if (atom_already_done[j])  // ring
    {
      r.add_bond(0, atom_already_done[j] - 1, b->btype());
      continue;
    }

    int natoms_before_addition = r.natoms();

    build_molecule(m, j, NOT_A_BOND, INVALID_ATOM_NUMBER, embedding, atom_already_done, r);

    r.add_bond(0, natoms_before_addition, b->btype());
  }

  if (remove_explicit_hydrogens_not_attached_to_isotope) {
    do_remove_explicit_hydrogens_not_attached_to_isotope(r, isotope_for_join_atom);
  }

  rgrcs.write_parent_molecule_if_needed(m, output);

  return rgrcs.do_output(r, 0, output);
}

static int
identify_atoms_in_substituent(const Molecule& m, atom_number_t ok_to_hit,
                              atom_number_t zatom, int* atom_already_done) {
#ifdef DEBUG_IDENTIFY_ATOMS_IN_SUBSTITUENT
  cerr << "identify_atoms_in_substituent:on to atom " << zatom << '\n';
#endif

  atom_already_done[zatom] = 1;

  const Atom* a = m.atomi(zatom);

  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);

    if (1 == atom_already_done[j]) {  // been this way before
      continue;
    }

    if (0 == atom_already_done[j]) {  // great, keep going
      ;
    } else if (ok_to_hit == j) {
      continue;
    } else if (2 == atom_already_done[j])  // bad, looped back to initial embedding
    {
      if (zatom == ok_to_hit) {  // coming off the attachment point, of course we hit
                                 // other matched atoms
        continue;
      }

#ifdef DEBUG_IDENTIFY_ATOMS_IN_SUBSTITUENT
      cerr << "From atom " << zatom << " get to " << j << " is loopback\n";
#endif
      return 0;
    }

    if (!identify_atoms_in_substituent(m, ok_to_hit, j, atom_already_done)) {
      return 0;
    }
  }

  return 1;
}

static int
add_attachment_point(Molecule& m,  // not const because the is_aromatic is non const
                     atom_number_t zatom, Molecule& r) {
  Atom* aa;
  if (!m.is_aromatic(zatom)) {
    aa = new Atom(m.atomic_number(zatom));
  } else if (6 == m.atomic_number(zatom)) {
    aa = new Atom(20);
  } else if (7 == m.atomic_number(zatom)) {
    aa = new Atom(11);
  } else {
    aa = new Atom(20);  // huh?
  }

  r.add(aa);

  return 1;
}

/*
  Atom zatom is in a ring. We must first ensure that there is no
  path back to other matched atoms
*/

static int
rgroup_do_ring(Molecule& m, RGroup_Construction_Current_State& rgrcs,
               const Set_of_Atoms& embedding, atom_number_t zatom,
               IWString_and_File_Descriptor& output) {
  if (verbose > 1) {
    cerr << "Trying to find ring R group from atom " << zatom
         << " ncon = " << m.ncon(zatom) << '\n';
  }

  int* atom_already_done = rgrcs.atom_already_done();

  int matoms = m.natoms();

  set_vector(atom_already_done, m.natoms(), 0);

  embedding.set_vector(atom_already_done, 2);  // 2 means in original embedding

  // If this loops back, we cannot process it.
  if (0 == identify_atoms_in_substituent(m, zatom, zatom, atom_already_done)) {
    return 1;
  }

  if (verbose > 1) {
    cerr << "Found OK substituent from " << zatom << '\n';
  }

  Molecule r;

  // When doing ring substituents, the atachment point is always included

  add_attachment_point(m, zatom, r);

  atom_already_done[zatom] = 0;  // so it does not get processed below

  int* xref = new_int(matoms);
  std::unique_ptr<int[]> free_xref(xref);
  int ndx = 1;

  for (int i = 0; i < matoms; i++) {
    if (1 != atom_already_done[i]) {
      continue;
    }

    r.add(new Atom(m.atomi(i)));

    xref[i] = ndx;
    ndx++;
  }

  if (isotope_for_join_atom) {
    r.set_isotope(xref[zatom], isotope_for_join_atom);
  }

  atom_already_done[zatom] = 1;  // now we want it processed

  for (const Bond* b : m.bond_list()) {
    if (1 != atom_already_done[b->a1()]) {
      continue;
    }

    if (1 != atom_already_done[b->a2()]) {
      continue;
    }

    atom_number_t r1 = xref[b->a1()];
    atom_number_t r2 = xref[b->a2()];

    r.add_bond(r1, r2, b->btype());
  }

  molecules_created++;

  if (!size_constraints_met(r)) {
    return 1;
  }

  rgrcs.write_parent_molecule_if_needed(m, output);

  return rgrcs.do_output(r, 0, output);
}

// Return the set of atoms attached to `zatom` which are not
// in `embedding` and which are not associated with a ring.
// Matched atoms are marked in `atom_already_done`.
Set_of_Atoms
UnmatchedConnections(Molecule& m,
                     const Set_of_Atoms& embedding,
                     const int* atom_already_done,
                     atom_number_t zatom) {
  Set_of_Atoms result;

  for (const Bond* b : m[zatom]) {
    const atom_number_t j = b->other(zatom);

    if (verbose > 1) {
      cerr << "From atom " << zatom << " atom " << j << " type "
           << m.smarts_equivalent_for_atom(j)
           << " already_done = " << atom_already_done[j]
           << " in embedding = " << embedding.contains(j) << '\n';
    }

    if (atom_already_done[j]) {
      cerr << "From atom " << zatom << " we find atom " << j
           << " already done. Possible cyclic\n";
      continue;
    }

    if (embedding.contains(j)) {
      continue;
    }

    if (m.in_same_ring(zatom, j)) {
      continue;
    }

    // Do not break biphenyl type linkages.
    if (m.nrings(zatom) && m.nrings(j) && b->nrings()) {
      continue;
    }

    result << j;
  }

  return result;
}

/*
  Slightly smarter about how to define the R group.

  This is designed to deal with the case where there may be
  multiple attachments to the atom. We don't want to produce
  unconnected fragments, but want to produce a fragment
  with the central atom included
*/

static int
rgroup_combine_same_atom_substituents(Molecule& m,
                                      RGroup_Construction_Current_State& rgrcs,
                                      const Set_of_Atoms& embedding, atom_number_t zatom,
                                      IWString_and_File_Descriptor& output) {
  if (verbose > 1) {
    cerr << "Trying to find R group from atom " << zatom << " ncon = " << m.ncon(zatom) << '\n';
  }

  int* atom_already_done = rgrcs.atom_already_done();
  std::fill_n(atom_already_done, m.natoms(), 0);

  const Set_of_Atoms atoms_in_substituent = UnmatchedConnections(m, embedding,
                atom_already_done, zatom);

  // nothing sprouting off here this time
  if (atoms_in_substituent.empty()) {
    return 0;
  }

  int n = atoms_in_substituent.number_elements();
  // cerr << "At atom " << zatom << " unmatched atoms are " << atoms_in_substituent << '\n';

  // Jun 2023. I am not sure why hcount should make a difference. If Hydrogens are
  // explicit, then they will be captured in `atoms_in_substituent`. If they are implicit
  // then they do not matter.
  // I have changed this, but am a little nervous about the effects.
  // if (n > 1 || m.hcount(zatom)) {
  if (n > 1 || (make_implicit_hydrogens_explicit && m.hcount(zatom))) {
    return rgroup_including_substituents(m, rgrcs, zatom, atoms_in_substituent, embedding,
                                         output);
  }

  // cerr << "Normal processing atom atom " << zatom << '\n';

  Molecule r;

  if (include_attachment_point_in_all_substituents) {
    if (!m.is_aromatic(zatom)) {
      m.add(new Atom(m.atomic_number(zatom)));
    } else if (6 == m.atomic_number(zatom)) {
      m.add(new Atom(20));
    } else if (7 == m.atomic_number(zatom)) {
      m.add(new Atom(11));
    } else {
      m.add(new Atom(20));  // huh?
    }

    atom_already_done[zatom] = 1;

    atom_number_t j = atoms_in_substituent[0];

    // const Bond * b = m.atomi(zatom)->bond_to_atom(j);
    const Bond* b = m.atomi(zatom)->bond_to_atom(zatom, j);

    build_molecule(m, j, b->btype(), zatom, embedding, atom_already_done, r);
  } else {
    build_molecule(m, atoms_in_substituent[0], NOT_A_BOND, INVALID_ATOM_NUMBER, embedding,
                   atom_already_done, r);
  }

  molecules_created++;

  if (!size_constraints_met(r)) {
    return 1;
  }

  rgrcs.write_parent_molecule_if_needed(m, output);

  return rgrcs.add_dummy_element_and_do_output(r, 0, output);
}

static int
rgroup(Molecule& m, RGroup_Construction_Current_State& rgrcs,
       const Set_of_Atoms& embedding, atom_number_t parent_atom,
       IWString_and_File_Descriptor& output) {
  int rc = 0;

  if (verbose > 1) {
    cerr << "Trying to find R group from atom " << parent_atom << " ncon = " << m.ncon(parent_atom) << '\n';
  }

  int* atom_already_done = rgrcs.atom_already_done();

  int ndx = 0;

  for (const Bond* b : m[parent_atom]) {
    atom_number_t j = b->other(parent_atom);

    if (verbose > 1) {
      cerr << "From atom " << parent_atom << " atom " << j
           << " already_done = " << atom_already_done[j]
           << " in embedding = " << embedding.contains(j) << '\n';
    }

    if (atom_already_done[j]) {
      cerr << "From atom " << parent_atom << " we find atom " << j
           << " already done. Possible cyclic\n";
      continue;
    }

    if (embedding.contains(j)) {
      continue;
    }

    if (m.in_same_ring(parent_atom, j)) {  // cannot process these
      continue;
    }

    Molecule r;

    rc += build_molecule(m, j, NOT_A_BOND, INVALID_ATOM_NUMBER, embedding,
                         atom_already_done, r);

    molecules_created++;

    if (!size_constraints_met(r)) {
      continue;
    }

    rgrcs.write_parent_molecule_if_needed(m, output);

    rgrcs.add_dummy_element_and_do_output(r, ndx, output);
    ++ndx;
  }

  return rc;
}

// Return true if `zatom` is attached via a ring bond to something
// outside `s`.
static int
is_within_ring(Molecule& m, const Set_of_Atoms& s, atom_number_t zatom) {
  (void)m.ring_membership();

  if (0 == m.nrings(zatom)) {
    return 0;
  }

  const Atom* a = m.atomi(zatom);

#ifdef DEBUG_IS_WITHIN_RING
  cerr << "is_within_ring:examining atom " << zatom << " acon " << acon << '\n';
#endif

  for (const Bond* b : *a) {
    if (0 == b->nrings()) {
      continue;
    }

    atom_number_t j = b->other(zatom);

    // ring bond to something outside the embedding, bad
    if (!s.contains(j)) {
      return 1;
    }
  }

#ifdef DEBUG_IS_WITHIN_RING
  cerr << "is_within_ring:atom " << zatom << " not part of ring\n";
#endif
  return 0;
}

// REturn true if any atom in `s` is bonded, via a ring bond, to
// an atom outside `s`.
static int
is_within_ring(Molecule& m, const Set_of_Atoms& s) {
  (void)m.ring_membership();

  for (const atom_number_t j : s) {
    if (0 == m.nrings(j)) {
      continue;
    }

    if (is_within_ring(m, s, j)) {
      return 1;
    }
  }

  return 0;
}

// #define DEBUG_START_RGROUP

/*
  We have an embedding in the molecule.
  For each of the query atoms in ATOMS, we need to find the R group
*/

static int
rgroups(Molecule& m, RGroup_Construction_Current_State& rgrcs,
        const Set_of_Atoms& embedding, const Set_of_Query_Atoms* atoms,
        IWString_and_File_Descriptor& output) {
  int rc = 0;

  int na = atoms->number_elements();
  for (int i = 0; i < na; i++) {
    int j = atoms->item(i);

#ifdef DEBUG_START_RGROUP
    cerr << "In query " << query_number << ", embedding " << embedding_number << " item "
         << i << " query atom " << j << " matches atom " << embedding->item(j) << '\n';
    cerr << " embedding " << *atoms << " atom " << j << 
#endif

    rgrcs.set_query_atom_number(i);

    if (is_within_ring(m, embedding, embedding[j])) {
      if (process_ring_substituents) {
        rc += rgroup_do_ring(m, rgrcs, embedding, embedding[j], output);
      }
    } else if (combine_substituents_on_same_atom) {
      rc += rgroup_combine_same_atom_substituents(m, rgrcs, embedding, embedding[j],
                                                  output);
    } else {
      rc += rgroup(m, rgrcs, embedding, embedding[j], output);
    }
  }

  return rc;
}

/*
  Even if just one atom is in a ring, we will discard the embedding
*/

static int
remove_embeddings_within_rings(Molecule& m, Substructure_Results& sresults) {
  for (int i = sresults.number_embeddings() - 1; i >= 0; i--) {
    const Set_of_Atoms* e = sresults.embedding(i);

    if (is_within_ring(m, *e)) {
      sresults.remove_embedding(i);
    }
  }

  return sresults.number_embeddings();
}

static int
rgroups(Molecule& m, Substructure_Hit_Statistics* q,
        RGroup_Construction_Current_State& rgrcs, const Set_of_Query_Atoms* atoms,
        IWString_and_File_Descriptor& output) {
  Substructure_Results sresults;

  int nhits = q->substructure_search(m, sresults);

  if (verbose > 1) {
    cerr << nhits << " hits to query '" << q->comment() << "'\n";
  }

  if (0 == nhits) {
    if (verbose) {
      cerr << "Only matched " << sresults.max_query_atoms_matched_in_search()
           << " query atoms when searching '" << m.name() << "'\n";
    }
    return 0;
  }

  if (!process_ring_substituents) {
    nhits = remove_embeddings_within_rings(m, sresults);
    if (nhits == 0) {
      if (verbose) {
        cerr << m.name() << "after removing matches in rings, no hits\n";
      }
      return 0;
    }
  }

  if (1 == nhits) {
    ;
  } else if (ignore_molecules_with_multiple_hits) {
    return 0;
  } else if (verbose > 1) {
    cerr << m.name() << ", " << nhits << " hits to query '" << q->comment() << "'\n";
  }

  int istop = 1;
  if (do_all_embeddings) {
    istop = nhits;
  }

  for (int i = 0; i < istop; i++) {
    const Set_of_Atoms* p = sresults.embedding(i);

    rgrcs.set_embedding_number(i);

    rgroups(m, rgrcs, *p, atoms, output);
  }

  return 1;
}

static resizable_array_p<Substructure_Hit_Statistics> queries;

static int
rgroups(Molecule& m, IWString_and_File_Descriptor& output) {
  int nq = queries.number_elements();

  int rc = 0;

  RGroup_Construction_Current_State rgrcs(m.natoms());

  for (int i = 0; i < nq; i++) {
    Substructure_Hit_Statistics* q = queries[i];

    const Set_of_Query_Atoms* sqa = query_atoms[i];

    rgrcs.set_query_number(i);

    if (0 == rgroups(m, q, rgrcs, sqa, output)) {
      continue;
    }

    rc++;
  }

  if (rc) {
    return rc;
  }

  if (verbose > 1) {
    cerr << m.name() << ", no hits\n";
  }

  molecules_with_no_hits++;

  return ignore_molecules_which_dont_match_query;
}

static void
preprocess(Molecule& m) {
  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  return;
}

static int
rgroups(data_source_and_type<Molecule>& input, IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (verbose > 1) {
      cerr << "Processing " << input.molecules_read() << " '" << m->name() << "'\n";
    }

    preprocess(*m);

    (void)rgroups(*m, output);

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

/*
  Looks like

  <q>.<m>
  or just
  <q>
*/

static int
parse_query_and_matched_atom_specification(const const_IWSubstring& r) {
  int nq = queries.number_elements();

  int query_number;

  if (r.numeric_value(query_number)) {
    if (query_number < 0 || query_number >= nq) {
      cerr << "Invalid query number '" << r << "' must be between 0 and " << (nq - 1)
           << "\n";
      return 0;
    }

    Substructure_Hit_Statistics* q = queries[query_number];
    int n = q->max_atoms_in_query();

    Set_of_Query_Atoms* s = query_atoms[query_number];
    for (int i = 0; i < n; i++) {
      s->add(i);
    }

    return 1;
  }

  // Must be of the form Q.N

  const_IWSubstring q, n;

  if (!r.split(q, '.', n)) {
    cerr << "Query atom specifications must be of the form 'Q.A' where Q is the query "
            "number and A is the atom\n";
    return 0;
  }

  if (!q.numeric_value(query_number) || query_number < 0 || query_number >= nq) {
    cerr << "Invalid query number '" << r << "' must be between 0 and " << (nq - 1)
         << "\n";
    return 0;
  }

  Substructure_Hit_Statistics* qry = queries[query_number];
  int aiq = qry->max_atoms_in_query();

  int a;
  if (!n.numeric_value(a) || a < 0 || a >= aiq) {
    cerr << "Query " << query_number << " has " << aiq << " atoms only, so '" << n
         << "' is an invalid matched atom specifier\n";
    return 0;
  }

  Set_of_Query_Atoms* s = query_atoms[query_number];
  s->add(a);

  if (verbose) {
    cerr << "Will report R Groups from atom " << a << " in query " << query_number
         << '\n';
  }

  return 1;
}

static int
rgroups(const char* fname, FileType input_type, IWString_and_File_Descriptor& output) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 2;
  }

  return rgroups(input, output);
}

static int
rgroups(int argc, char** argv) {
  Command_Line cl(argc, argv, "A:E:vi:d:r:z:q:s:c:C:HhaeI:ukyb");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(4);
  }

  if (!process_elements(cl)) {
    cerr << "Cannot process -E specifiers\n";
    usage(5);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (!cl.option_present('i')) {
    if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot determine input type(s)\n";
      return 6;
    }
  } else if (!process_input_type(cl, input_type)) {
    usage(2);
  }

  if (cl.option_present('z')) {
    int i = 0;
    const_IWSubstring z;

    while (cl.value('z', z, i++)) {
      if ('i' == z) {
        ignore_molecules_which_dont_match_query = 1;
        if (verbose) {
          cerr << "Molecules which don't match substructure queries will be ignored\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        usage(3);
      }
    }
  }

  if (cl.option_present('d')) {
    const_IWSubstring d = cl.string_value('d');

    int isave = auto_create_new_elements();
    set_auto_create_new_elements(1);

    dummy_atom_element = get_element_from_symbol_no_case_conversion(d);
    if (nullptr == dummy_atom_element) {
      dummy_atom_element = create_element_with_symbol(d);
    }

    set_auto_create_new_elements(isave);

    if (verbose) {
      cerr << "An embedding dummy atom of type '" << dummy_atom_element->symbol()
           << "' will be added to all groups\n";
    }
  }

  if (cl.option_present('I')) {
    if (!cl.value('I', isotope_for_join_atom) || isotope_for_join_atom <= 0) {
      cerr << "The isotope to apply option (-I) must be a valid isotope number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will label the join point with isotope " << isotope_for_join_atom << '\n';
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', min_atoms_in_rgroup) || min_atoms_in_rgroup < 1) {
      cerr << "The minimum number of atoms in an rgroup (-c) option must be a whole +ve "
              "number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will discard rgroups with fewer than " << min_atoms_in_rgroup
           << " atoms\n";
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_atoms_in_rgroup) || max_atoms_in_rgroup < 1 ||
        max_atoms_in_rgroup < min_atoms_in_rgroup) {
      cerr << "The maximum number of atoms in an rgroup (-C) option must be a whole +ve "
              "number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will discard rgroups with more than " << max_atoms_in_rgroup << " atoms\n";
    }
  }

  if (cl.option_present('H')) {
    make_implicit_hydrogens_explicit = cl.option_count('H');

    if (verbose) {
      cerr << "Explicit hydrogens added\n";
    }
  }

  if (cl.option_present('h')) {
    if (0 == isotope_for_join_atom) {
      cerr << "The remove explicit hydrogens option (-h) option must be combined with "
              "the -I option\n";
      usage(3);
    }

    remove_explicit_hydrogens_not_attached_to_isotope = isotope_for_join_atom;

    if (verbose) {
      cerr << "After processing, will remove explicit H not attached to join points\n";
    }
  }

  if (cl.option_present('a')) {
    do_all_embeddings = 1;

    if (verbose) {
      cerr << "Will process all embeddings\n";
    }
  }

  if (cl.option_present('e')) {
    combine_substituents_on_same_atom = 0;

    if (verbose) {
      cerr << "Multiple substituents on one atom done separately\n";
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, queries, verbose, 'q')) {
      cerr << "Cannot discern queries from -q switch\n";
      usage(21);
    }
  }

  if (cl.option_present('s')) {
    int i = 0;
    const_IWSubstring s;
    while (cl.value('s', s, i++)) {
      Substructure_Hit_Statistics* q = new Substructure_Hit_Statistics;
      if (!q->create_from_smarts(s)) {
        cerr << "Cannot parse smarts '" << s << "'\n";
        return 5;
      }

      queries.add(q);
    }
  }

  int nq = queries.number_elements();
  if (0 == nq) {
    cerr << "Must specify one or more queries via -q or -s options\n";
    usage(24);
  }

  if (verbose) {
    cerr << "Processing " << nq << " queries\n";
  }

  int process_symmetry_equivalent_matches = 0;

  if (cl.option_present('b')) {
    process_symmetry_equivalent_matches = 1;

    if (verbose) {
      cerr << "Will process symmetry equivalent matches\n";
    }
  }

  query_stats = new extending_resizable_array<int>[nq];
  rgroups_found = new IW_STL_Hash_Map_int*[nq];
  for (int i = 0; i < nq; i++) {
    Substructure_Hit_Statistics* q = queries[i];
    if (!process_symmetry_equivalent_matches) {
      q->set_do_not_perceive_symmetry_equivalent_matches(1);
    }
    rgroups_found[i] = new IW_STL_Hash_Map_int[q->max_atoms_in_query()];
  }

  if (cl.option_present('u')) {
    for (int i = 0; i < nq; i++) {
      queries[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('k')) {
    include_attachment_point_in_all_substituents = 1;

    if (verbose) {
      cerr << "Will include attachment point in all R groups\n";
    }
  }

  if (cl.option_present('y')) {
    process_ring_substituents = 1;

    if (verbose) {
      cerr << "Will process substituents involved in rings\n";
    }
  }

  // For each query, we must assemble the associated query atom numbers

  if (!cl.option_present('r')) {  // will get all sites on all query atoms
    ;
  } else if (cl.option_count('r') < cl.option_count('q')) {
    cerr << "Must specify which query atoms via -r option(s)\n";
    usage(25);
  }

  query_atoms.resize(nq);
  for (int i = 0; i < nq; i++) {
    Set_of_Query_Atoms* tmp = new Set_of_Query_Atoms;
    query_atoms.add(tmp);
  }

  if (cl.option_present('r')) {
    int i = 0;
    const_IWSubstring r;
    while (cl.value('r', r, i++)) {
      if (!parse_query_and_matched_atom_specification(r)) {
        return 3;
      }
    }
  } else {
    int natoms =
        queries[0]->highest_initial_atom_number();  // or should we use max_atoms_in_query

    for (int i = 1; i < nq; i++) {
      if (queries[i]->highest_initial_atom_number() != natoms) {
        cerr << "Without the -r option, all queries must have same number of atoms "
             << queries[i]->highest_initial_atom_number() << " vs " << natoms << '\n';
        return i;
      }
    }

    for (int i = 0; i < nq; i++) {
      IWString r;
      r << i;
      parse_query_and_matched_atom_specification(r);
    }
  }
  // cerr << "Got to rerqwer\n";

  if (0 == cl.number_elements()) {
    cerr << "No input files specified\n";
    return 2;
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!rgroups(cl[i], input_type, output)) {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << molecules_read << " molecules read\n";
    cerr << molecules_created << " molecules created\n";
    cerr << molecules_written << " molecules written\n";
    if (molecules_with_no_hits) {
      cerr << molecules_with_no_hits << " molecules matched no queries\n";
    }

    for (int i = 0; i < nq; i++) {
      const extending_resizable_array<int>& qi = query_stats[i];

      Substructure_Hit_Statistics* qry = queries[i];

      int jstop = qi.number_elements();
      if (jstop > qry->max_atoms_in_query()) {
        jstop = qry->max_atoms_in_query();
      }

      const IW_STL_Hash_Map_int* h = rgroups_found[i];

      for (int j = 0; j < jstop; j++) {
        cerr << "Query " << i << " matched atom " << j << " hit " << qi[j] << " times\n";

        const IW_STL_Hash_Map_int& hj = h[j];
        for (IW_STL_Hash_Map_int::const_iterator l = hj.begin(); l != hj.end(); ++l) {
          cerr << i << '.' << j << ' ' << (*l).second << " instances of '" << (*l).first
               << "'\n";
        }
      }
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = rgroups(argc, argv);

  return rc;
}
