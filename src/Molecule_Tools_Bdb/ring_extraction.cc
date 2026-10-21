#include "ring_extraction.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/primes.h"

#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/target.h"

namespace smi2rings_bdb {
using std::cerr;

static bool
AllValuesAre(const int* v, const int n, int value) {
  for (int i = 0; i < n; ++i) {
    if (v[i] != value) {
      return false;
    }
  }

  return true;
}

/*
  Atom OUTSIDE_RING is an unsaturated atom outside the ring (attached to
  atom IN_RING). Add its neighbours
*/

static int
do_continue_chain_neighbours_through_multiple_bonds(Molecule& m, atom_number_t in_ring,
                                                    atom_number_t outside_ring,
                                                    int* in_set, int uid) {
  in_set[outside_ring] = uid;

  const Atom* a = m.atomi(outside_ring);

  for (const Bond* b : *a) {
    atom_number_t j = b->other(outside_ring);

    if (j == in_ring) {
      continue;
    }

    in_set[j] = uid;
  }

  return 1;
}

RingExtractionData::RingExtractionData(Molecule& mol,
                                       const RingExtractionOptions& options)
    : _m(mol), _options(options) {
  const int matoms = _m.natoms();
  _in_system = new_int(matoms, -1);
  _atype = nullptr;
  _ring_system_size = _m.nrings() ? new int[_m.nrings()] : nullptr;

  _xref.reset(new int[matoms]);

  // Original version had isotopic atoms saved here, but that breaks the
  // case of the atom typing including the isotopic label.
  // So isotopes in the incoming molecule are saved and removed in
  // AssignAtomTypes.
}

RingExtractionData::~RingExtractionData() {
  if (_in_system != nullptr) {
    delete[] _in_system;
  }
  if (_atype != nullptr) {
    delete[] _atype;
  }
  if (_ring_system_size != nullptr) {
    delete[] _ring_system_size;
  }
}

int
RingExtractionData::AssignAtomTypes(Atom_Typing_Specification& atom_typing) {
  _atype = new uint32_t[_m.natoms()];
  const int rc = atom_typing.assign_atom_types(_m, _atype);

  if (_m.number_isotopic_atoms() > 0) {
    _starting_isotopes = _m.GetIsotopes();
    _m.transform_to_non_isotopic_form();
  }

  return rc;
}

int
RingExtractionData::MaybeRevertInitialIsotopes(Molecule& m) {
  if (!_starting_isotopes) {
    return 0;
  }

  return m.set_isotopes(_starting_isotopes.get());
}

// Modified to also detect S+ and O+
// C12=C(ON=C1OC)CC[S+](C)C2 CHEMBL278618
int
RingExtractionData::IdentifyThreeConnectedAromaticSulphur(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    const Atom& a = m[i];

    if (a.formal_charge() == 1 && (a.atomic_number() == 8 || a.atomic_number() == 16)) {
    } else if (a.atomic_number() == 16 && a.ncon() == 3 && m.is_aromatic(i)) {
    } else {
      continue;
    }

    _aromatic_sulphur << i;
  }

  return _aromatic_sulphur.size();
}

int
RingExtractionData::AddConnectedAtomToPositiveAromaticNitrogen(Molecule& m,
                                                               int uid) const {
  int rc = 0;

  for (const atom_number_t n : _positive_aromatic_nitrogen) {
    //  cerr << "Atom " << j << " +ve arom N, flag " << uid << " val " << in_set[j] <<
    //  '\n';
    if (uid != _in_system[n]) {
      continue;
    }

    const Atom* a = m.atomi(n);
    if (a->ncon() == 2) {
      continue;
    }

    for (const Bond* b : *a) {
      atom_number_t l = b->other(n);

      if (uid == _in_system[l]) {
        continue;
      }

      if (m.is_ring_atom(l)) {  // hmmm, is this correct?
        continue;
      }

      _in_system[l] = uid;

      rc++;
    }
  }

  return rc;
}

int
RingExtractionData::ConvertPositiveAroamticNitrogenToNeutral(Molecule& m, int uid) const {
  int rc = 0;

  for (atom_number_t j : _positive_aromatic_nitrogen) {
    if (uid != _in_system[j]) {
      continue;
    }

    m.set_formal_charge(j, 0);

    rc++;
  }

  return rc;
}

static uint32_t
compute_environment_isotope(Molecule& m, atom_number_t zatom, const Bond* b) {
  atom_number_t j = b->other(zatom);

  int jcon = m.ncon(j);

  atomic_number_t jz = m.atomic_number(j);

  if (1 == jcon && b->is_single_bond()) {
    if (6 == jz) {
      return 1;
    }
    if (7 == jz) {
      return 2;
    }
    if (8 == jz) {
      return 3;
    }
    if (9 == jz) {
      return 4;
    }
    if (16 == jz) {
      return 5;
    }
    if (17 == jz) {
      return 6;
    }
    if (35 == jz) {
      return 7;
    }
    if (53 == jz) {
      return 8;
    }

    return 9;
  }

  if (1 == jcon)  // double bond
  {
    if (6 == jz) {
      return 10;
    }
    if (7 == jz) {
      return 11;
    }
    if (8 == jz) {
      return 12;
    }
    if (16 == jz) {
      return 13;
    }

    return 14;
  }

  if (!b->is_single_bond())  // some kind of multiply connected double bond
  {
    if (6 == jz) {
      return 15;
    }
    if (7 == jz) {
      return 16;
    }

    return 17;
  }

  // All the single bond cases

  if (8 == jz) {
    return 20;
  }

  if (16 == jz) {
    return 21;
  }

  if (6 == jz) {
    if (4 == jcon) {
      return 22;
    }

    if (2 == jcon) {
      return 23;
    }

    if (m.nbonds(j) == jcon) {
      return 24;
    }

    return 25;
  }

  if (7 == jz) {
    if (2 == jcon) {
      return 30;
    }

    if (m.nbonds(j) == jcon) {
      return 31;
    }

    return 32;
  }

  return 33;  // huh, what is this?
}

static int
do_label_attachment_points_with_existing_isotope(Molecule& m, const int* in_system,
                                                 int uid) {
  m.transform_to_non_isotopic_form();

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (uid != in_system[i]) {
      continue;
    }

    const Atom& a = m[i];

    if (2 == a.ncon()) {
      continue;
    }

    isotope_t to_assign = 0;
    for (const Bond* b : a) {
      atom_number_t k = b->other(i);

      if (uid == in_system[k]) {
        continue;
      }

      if (isotope_t iso = m.isotope(k); iso > to_assign) {
        to_assign = iso;
      }
    }

    m.set_isotope(i, to_assign);
  }

  return 1;
}

static int
do_label_attachment_points_with_environment(Molecule& m, const int* in_system, int uid,
                                            const uint32_t* atype) {
  m.transform_to_non_isotopic_form();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    if (uid != in_system[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    int acon = a->ncon();

    if (2 == acon) {
      continue;
    }

    for (int j = 0; j < acon; j++) {
      const Bond* b = a->item(j);

      atom_number_t k = b->other(i);

      if (uid == in_system[k]) {
        continue;
      }

      isotope_t iso;
      if (nullptr != atype) {
        iso = atype[k];
      } else {
        iso = compute_environment_isotope(m, i, b);
      }

      m.set_isotope(i, iso);
    }
  }

  return 1;
}

/*
  No differentiation as to the number of connections outside the system
*/

static int
do_label_attachment_points_1(Molecule& m, const int* in_system, int uid,
                             const RingExtractionOptions& _options) {
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (uid != in_system[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    for (const Bond* b : *a) {
      atom_number_t k = b->other(i);

      if (uid == in_system[k]) {
        continue;
      }

      if (b->is_single_bond()) {
        m.set_isotope(i, 1);
      } else if (_options._add_double_bonded_neighbors_to_rings) {
      } else {
        m.set_isotope(i, 2);
      }

      rc++;
    }
  }

  return rc;
}

static int
do_label_attachment_points_2(Molecule& m, const int* in_system, int uid) {
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++) {
    if (uid != in_system[i]) {
      continue;
    }

    const Atom* a = m.atomi(i);

    int acon = a->ncon();

    if (2 == acon && m.is_ring_atom(i)) {
      continue;
    }

    int connections_outside_ring = 0;

    for (int j = 0; j < acon; j++) {
      atom_number_t k = a->other(i, j);

      if (uid != in_system[k]) {
        connections_outside_ring++;
      }
    }

    if (connections_outside_ring) {
      m.set_isotope(i, connections_outside_ring);
      rc++;
    }
  }

  return rc;
}

static int
do_label_attachment_points(Molecule& m, const int* in_system, int uid,
                           const RingExtractionOptions& _options) {
  m.transform_to_non_isotopic_form();

  if (1 == _options._label_ring_atoms) {
    return do_label_attachment_points_1(m, in_system, uid, _options);
  } else {
    return do_label_attachment_points_2(m, in_system, uid);
  }
}

/*static int
is_spiro_fused (Molecule & m,
                const Ring & r,
                const int * in_system)
{
  int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
  }
}*/

static int
embedding_contains_ring_atom(Molecule& m, const Set_of_Atoms& e) {
  int n = e.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = e[i];

    if (m.is_ring_atom(j)) {
      return 1;
    }
  }

  return 0;
}

//  Identifying functional groups is hard, because an atom can be part of
//  several functional groups. We use prime numbers to allow this.

int
RingExtractionData::IdentifyFunctionalGroups(
    Molecule& m, resizable_array_p<Substructure_Query>& queries) {
  _functional_group = std::make_unique<uint32_t[]>(m.natoms());

  Molecule_to_Match target(&m);

  int rc = 0;

  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      if (!embedding_contains_ring_atom(m, *embedding)) {
        continue;
      }

      int p = primes[rc];
      rc++;

      for (atom_number_t l : *embedding) {
        if (0 == _functional_group[l]) {
          _functional_group[l] = p;
        } else {
          _functional_group[l] *= p;
        }
      }
    }
  }

  return rc;
}

//  In order to avoid passing around lots of arguments, create a class

class Ring_System_Info {
 private:
  const RingExtractionOptions& _options;
  int _rings_in_system;
  int _keep_this_ring_system;
  int _uid_to_assign;
  int* _ring_already_done;

 public:
  Ring_System_Info(int, const RingExtractionOptions& options);
  ~Ring_System_Info();

  int
  ring_already_done(int r) const {
    return _ring_already_done[r];
  }

  int*
  ring_already_done() {
    return _ring_already_done;
  }

  int
  uid_to_assign() const {
    return _uid_to_assign;
  }

  void initialise(int u);

  int add_ring(const Ring* r);
  int already_processed(const Ring* r);

  int
  keep_this_ring_system() const {
    return _keep_this_ring_system;
  }

  void
  set_keep_this_ring_system(int s) {
    _keep_this_ring_system = s;
  }

  int
  rings_in_system() const {
    return _rings_in_system;
  }
};

Ring_System_Info::Ring_System_Info(int nr, const RingExtractionOptions& options)
    : _options(options) {
  _ring_already_done = new_int(nr);

  return;
}

Ring_System_Info::~Ring_System_Info() {
  delete[] _ring_already_done;

  return;
}

void
Ring_System_Info::initialise(int u) {
  _rings_in_system = 0;
  _keep_this_ring_system = 1;
  _uid_to_assign = u;

  return;
}

int
Ring_System_Info::add_ring(const Ring* r) {
  _rings_in_system++;

  if (r->number_elements() > _options._max_ring_size_to_consider) {
    _keep_this_ring_system = 0;
  }

  assert(0 == _ring_already_done[r->ring_number()]);

  _ring_already_done[r->ring_number()] = 1;

  return _keep_this_ring_system;
}

/*
  Are any of the atoms in R already set in the IN_SYSTEM array?
*/

static int
any_members_non_negative_in_array(const int* in_system, const Ring& ring) {
  for (atom_number_t r : ring) {
    if (in_system[r] >= 0) {
      return 1;
    }
  }

  return 0;
}

/*int
Ring_System_Info::already_processed(const Ring * r)
{
  return _ring_already_done[r->ring_number()];
}*/

static atom_number_t
identify_all_double_bonds_outside_ring(Molecule& m, const Ring* r, atom_number_t zatom,
                                       const Ring_System_Info& rsi, int* in_system) {
  const Atom* a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++) {
    const Bond* b = a->item(i);
    if (!b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(zatom);

    if (r->contains(j)) {  // not interested within ring
      continue;
    }

    in_system[j] = rsi.uid_to_assign();
    rc++;
  }

  return rc;
}

/*
  We are bridging across spiro rings, or adding all double bonds to rings.
  Either of these actions may bring other rings into a group. Since we
  don't know how large this may make a ring system, we need an infinite loop
*/

static int
assign_atoms_to_ring_system(Molecule& m, Ring_System_Info& rsi, int zring, int* in_system,
                            const RingExtractionOptions& _options) {
  const Ring* r = m.ringi(zring);

  r->set_vector(in_system, rsi.uid_to_assign());

  rsi.add_ring(r);

  if (_options._add_double_bonded_neighbors_to_rings) {
    int ring_size = r->number_elements();

    for (int i = 0; i < ring_size; i++) {
      atom_number_t j = r->item(i);

      const Atom* aj = m.atomi(j);

      if (2 == aj->ncon()) {
        continue;
      }

      if (aj->nbonds() == aj->ncon()) {  // fully saturated
        continue;
      }

      identify_all_double_bonds_outside_ring(m, r, j, rsi, in_system);
    }
  }

  int nr = m.nrings();

  // Now that we have added the doubly bonded atoms, any ring
  // that touches any set atom, gets added

  for (int i = 0; i < nr; i++) {
    if (i == zring) {
      continue;
    }

    if (rsi.ring_already_done(i)) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    if (any_members_non_negative_in_array(in_system, *ri)) {
      assign_atoms_to_ring_system(m, rsi, i, in_system, _options);
    }
  }

  return rsi.keep_this_ring_system();
}

// For atoms that are part of the scaffold, set their `in_system` value to 0,
// other atoms get -1 - they probably have this value already due to defaults in
// PerMoleculeData
static int
assign_atoms_to_scaffold(Molecule& m, int* in_system, int* ring_system_size) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> spinach = std::make_unique<int[]>(matoms);
  m.identify_spinach(spinach.get());

  for (int i = 0; i < matoms; ++i) {
    if (spinach[i] == 0) {
      in_system[i] = 0;
    } else if (in_system[i] >= 0) {
    } else {
      in_system[i] = -1;
    }
  }

  ring_system_size[0] = m.nrings();

  return 1;  // One "ring system" found
}

static int
assign_atoms_to_ring_systems(Molecule& m, int* in_system, int* ring_system_size,
                             const RingExtractionOptions& _options) {
  int matoms = m.natoms();

  int nr = m.nrings();

  Ring_System_Info rsi(nr, _options);

  int uid_to_assign = 0;

  for (int i = 0; i < nr; i++) {
    if (rsi.ring_already_done(i)) {
      continue;
    }

    rsi.initialise(uid_to_assign);

    if (assign_atoms_to_ring_system(m, rsi, i, in_system, _options)) {
      ring_system_size[uid_to_assign] = rsi.rings_in_system();
      uid_to_assign++;
    } else  // some reason to not keep that ring system
    {
      for (int j = 0; j < matoms; j++) {
        if (uid_to_assign == in_system[j]) {
          in_system[j] = -1;
        }
      }
    }
  }

  return uid_to_assign;  // the number of ring systems we have identified
}

int
RingExtractionData::AddImplicitHydrogenToAromaticSulphur(Molecule& m,
                                                         int ring_system) const {
  if (_aromatic_sulphur.empty()) {
    return 0;
  }

  int rc = 0;
  for (atom_number_t s : _aromatic_sulphur) {
    if (_in_system[s] != ring_system) {
      continue;
    }
    m.set_implicit_hydrogens(s, 1, 1);
    ++rc;
  }

  return rc;
}

// Do we need to modify any charged aromatic nitrogen atoms.
// `positive_aromatic_nitrogen` contains a list of '[n+]' atoms.
// `ring_system` is the ring system we are processing.
int
RingExtractionData::MaybeConvertPositiveAromaticNitrogen(Molecule& m,
                                                         int ring_system) const {
  if (_positive_aromatic_nitrogen.empty()) {
    return 0;
  }

  if (_options._add_connected_atom_to_positive_aromatic_nitrogen) {
    AddConnectedAtomToPositiveAromaticNitrogen(m, ring_system);
  } else if (_options._convert_positive_aromatic_nitrogen_to_neutral) {
    ConvertPositiveAroamticNitrogenToNeutral(m, ring_system);
  }

  return 1;
}

int
RingExtractionData::MaybeAddOutsideRingAtoms(Molecule& m, int ring_system_number) {
  if (!_functional_group) {
    return 0;
  }

  return AddExtraRingSingleBondedNeighgours(m, ring_system_number);
}

int
RingExtractionData::IdentifyPositiveAromaticNitrogen(Molecule& m) {
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    if (1 != a->formal_charge()) {
      continue;
    }

    if (3 != a->ncon()) {
      continue;
    }

    if (7 != a->atomic_number()) {
      continue;
    }

    if (!m.is_aromatic(i)) {
      continue;
    }

    if (1 == m.nrings(i)) {
      _positive_aromatic_nitrogen.add(i);
    }
  }

  return _positive_aromatic_nitrogen.number_elements();
}

// #define DEBUG_SUBSET_CREATION

//  Create the unique smiles for a given subset `uid` of the molecule.
int
RingExtractionData::CreateSubset(Molecule& m, int uid, IWString& unique_smiles) const {
  Molecule tmp;

  m.create_subset(tmp, _in_system, uid, _xref.get());

  // We could be more precise about this - change just the atoms
  // that have changed bonding.
  // O[C@@H]1C[C@H]2[N+](C)(C)[C@@H](C1)[C@H]1[C@@H]2O1 CHEMBL3706465
  const int matoms = tmp.natoms();
  for (int i = 0; i < matoms; ++i) {
    tmp.set_implicit_hydrogens_known(i, 0);
  }

  // cerr << tmp.unique_smiles() << " after formation\n";
  // tmp.invalidate_smiles();
  tmp.recompute_implicit_hydrogens();
  // cerr << tmp.unique_smiles() << " after implicit H " << m.valence_ok() << '\n';

  for (atom_number_t s : _aromatic_sulphur) {
    if (_xref[s] >= 0) {
      tmp.set_implicit_hydrogens(_xref[s], 1, 1);
    }
  }

#ifdef DEBUG_SUBSET_CREATION
  tmp.debug_print(cerr);
  cerr << "Smiles for subset '" << tmp.smiles() << "' unique '";
  cerr << tmp.unique_smiles() << "'\n";
#endif

  if (_options._remove_invalid_chiral_centres) {
    lillymol::do_remove_invalid_chiral_centres(tmp);
  }

  // Restore the caller's setting even if SMILES generation throws.
  class CoordinateScope {
   private:
    int _saved;
    bool _changed;

   public:
    explicit CoordinateScope(bool include)
        : _saved(lillymol::include_coordinates_with_smiles()),
          _changed(_saved != include) {
      if (_changed) {
        lillymol::set_include_coordinates_with_smiles(include);
      }
    }

    ~CoordinateScope() {
      if (_changed) {
        lillymol::set_include_coordinates_with_smiles(_saved);
      }
    }
  } scope(_options._write_3d_smiles_if_needed &&
          tmp.highest_coordinate_dimensionality() == 3);

  unique_smiles = tmp.unique_smiles();

  return 1;
}

static int
get_next_prime_factor(int znumber, int& pindex, int istop) {
  while (1) {
    int p = primes[pindex++];

    if (0 == znumber % p) {
      return p;
    }

    if (p > istop) {
      return 0;
    }
  }

  return 0;  // should not come to here
}

/*
  Add all the atoms which were hit by a functional group query to the set.
  Atom K is part of the set. Any atom whose functional_group value is the
  same as K's will be added
*/

static void
add_functional_group(int matoms, int k, const uint32_t* functional_group, int* in_set,
                     int uid) {
  int istop = static_cast<int>(sqrt(static_cast<double>(functional_group[k])));

  int i = 0;
  int p;
  while (0 != (p = get_next_prime_factor(functional_group[k], i, istop))) {
    for (int i = 0; i < matoms; i++) {
      if (0 == functional_group[i]) {
        continue;
      }

      if (0 == (functional_group[i] % p)) {
        in_set[i] = uid;
      }
    }
  }

  return;
}

/*
  Must make sure we don't accidentially add too many things
  We process the atoms for which UID == IN_SYSTEM[i]
*/

static constexpr int kTmpid = -55512;

int
RingExtractionData::AddExtraRingSingleBondedNeighgours(Molecule& m, int uid) {
  int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++) {
    if (uid != _in_system[i]) {  // just examine our subset
      continue;
    }

    Atom* a = const_cast<Atom*>(m.atomi(i));  // loss of const OK

    if (a->ncon() == 2) {  // cannot have branches outside the system
      continue;
    }

    for (const Bond* b : *a) {
      if (!b->is_single_bond()) {
        continue;
      }

      atom_number_t k = b->other(i);
      if (uid == _in_system[k]) {  // we are looking for atoms outside the ring system
        continue;
      }

      // adding all singly bonded attachments
      if (2 == _options._add_single_bonded_neighbours_to_rings) {
        _in_system[k] = kTmpid;
        if (_functional_group[i] && _functional_group[k]) {
          add_functional_group(matoms, i, _functional_group.get(), _in_system, kTmpid);
        }

        rc++;
      } else if (_functional_group[i] &&
                 _functional_group[k])  // just adding extra's which are part of a
                                        // functional group
      {
        _in_system[k] = kTmpid;
        add_functional_group(matoms, i, _functional_group.get(), _in_system, kTmpid);
        rc++;
      } else if (_options._add_all_single_atom_ring_substituents && 1 == m.ncon(k)) {
        _in_system[k] = kTmpid;
        rc++;
      }

      if (_options._continue_chain_neighbours_through_multiple_bonds &&
          0 == m.nrings(k)) {
        Atom* ak = const_cast<Atom*>(m.atomi(k));

        if (ak->nbonds() > ak->ncon()) {
          do_continue_chain_neighbours_through_multiple_bonds(m, i, k, _in_system,
                                                              kTmpid);
        }
        rc++;
      }
    }
  }

  if (rc) {
    for (int i = 0; i < matoms; i++) {
      //    cerr << "Atom " << i << " set " << in_system[i] << '\n';
      if (kTmpid == _in_system[i]) {
        _in_system[i] = uid;
      }
    }
  }

  // a zero return code is OK - just means no extra-ring double bonds found
  return rc;
}

bool
RingExtractionData::AlreadySeen(const IWString& unique_smiles) {
  if (_found_this_molecule.contains(unique_smiles)) {
    return true;
  }

  _found_this_molecule.insert(unique_smiles);

  return false;
}

int
RingExtractor::Prepare(Molecule& m, RingExtractionData& data) {
  if (_atom_typing.active() && !data.AssignAtomTypes(_atom_typing)) {
    return 0;
  }
  if (!_queries.empty()) {
    data.IdentifyFunctionalGroups(m, _queries);
  }
  if (_options._add_connected_atom_to_positive_aromatic_nitrogen ||
      _options._convert_positive_aromatic_nitrogen_to_neutral) {
    data.IdentifyPositiveAromaticNitrogen(m);
  }
  data.IdentifyThreeConnectedAromaticSulphur(m);
  return 1;
}

int
RingExtractor::IdentifySystems(Molecule& m, RingExtractionData& data) const {
  if (_options._use_scaffolds) {
    return assign_atoms_to_scaffold(m, data._in_system, data._ring_system_size);
  }
  return assign_atoms_to_ring_systems(m, data._in_system, data._ring_system_size,
                                      _options);
}

bool
RingExtractor::PrepareSystem(Molecule& m, RingExtractionData& data, int uid) const {
  if (data._ring_system_size[uid] > _options._max_ring_system_size_to_consider) {
    return false;
  }
  data.MaybeConvertPositiveAromaticNitrogen(m, uid);
  data.AddImplicitHydrogenToAromaticSulphur(m, uid);
  data.MaybeAddOutsideRingAtoms(m, uid);
  return !_options._discard_molecules_where_ring_is_whole_molecule ||
         !AllValuesAre(data._in_system, m.natoms(), uid);
}

void
RingExtractor::LabelEnvironment(Molecule& m, RingExtractionData& data, int uid) const {
  do_label_attachment_points_with_environment(m, data._in_system, uid, data._atype);
}

void
RingExtractor::LabelExistingIsotopes(Molecule& m, RingExtractionData& data,
                                     int uid) const {
  do_label_attachment_points_with_existing_isotope(m, data._in_system, uid);
}

void
RingExtractor::LabelAttachments(Molecule& m, RingExtractionData& data, int uid) const {
  do_label_attachment_points(m, data._in_system, uid, _options);
}

RingKeyTypes
RingExtractor::KeyTypes(bool building) const {
  RingKeyTypes result;
  if (_options._label_attachment_points_with_environment) {
    result._types[result._size++] = RingKeyType::kContextual;
  }
  // The legacy builder does not store the isoatt lookup representation.
  if (!building && _options._label_attachment_points_with_existing_isotope) {
    result._types[result._size++] = RingKeyType::kExistingIsotopes;
  }
  if (_options._label_attachment_points) {
    result._types[result._size++] = RingKeyType::kAttachment;
  }
  if (_options._always_process_unsubstituted_ring) {
    result._types[result._size++] = RingKeyType::kBare;
  }
  return result;
}

int
RingExtractor::CreateKey(Molecule& m, RingExtractionData& data, int uid, RingKeyType type,
                         IWString& key, bool building) const {
  switch (type) {
    case RingKeyType::kContextual:
      LabelEnvironment(m, data, uid);
      break;
    case RingKeyType::kExistingIsotopes:
      LabelExistingIsotopes(m, data, uid);
      break;
    case RingKeyType::kAttachment:
      LabelAttachments(m, data, uid);
      break;
    case RingKeyType::kBare:
      // Preserve the legacy distinction: bare store always removes isotopes;
      // bare-only lookup can intentionally use input isotopes in its key.
      if (building) {
        m.transform_to_non_isotopic_form();
      }
      break;
  }
  return data.CreateSubset(m, uid, key);
}

}  // namespace smi2rings_bdb
