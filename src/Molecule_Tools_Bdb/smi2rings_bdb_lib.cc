#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "Molecule_Lib/molecule.h"

#include "smi2rings_bdb_lib.h"

namespace mol2rings {

using std::cerr;

// Return true if all members of the array `v` are `value`.
// There is probably a library function to do this...
static bool
AllValuesAre(const int* v,
             const int n,
             int value) {
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

  for (const Bond * b : *a) {
    atom_number_t j = b->other(outside_ring);

    if (j == in_ring) {
      continue;
    }

    in_set[j] = uid;
  }

  return 1;
}
PerMoleculeData::PerMoleculeData(Molecule& mol) : m(mol) {
  const int matoms = m.natoms();
  in_system = new_int(matoms, -1);
  atype = nullptr;
  if (m.nrings() == 0) {
    rsu = nullptr;
    ring_system_size = 0;
  } else {
    rsu = new Ring_System_Uniqueness[m.nrings()];
    ring_system_size = new int[m.nrings()];
  }

  if (m.number_isotopic_atoms() > 0) {
    starting_isotopes = m.GetIsotopes();
    m.transform_to_non_isotopic_form();
  }

  xref.reset(new int[matoms]);
}

PerMoleculeData::~PerMoleculeData() {
  if (in_system != nullptr) {
    delete[] in_system;
  }
  if (atype != nullptr) {
    delete[] atype;
  }
  if (rsu != nullptr) {
    delete[] rsu;
  }
  if (ring_system_size != nullptr) {
    delete[] ring_system_size;
  }
}

int
PerMoleculeData::AssignAtomTypes(Atom_Typing_Specification& atom_typing) {
  atype = new int[m.natoms()];
  return atom_typing.assign_atom_types(m, atype);
}

int
PerMoleculeData::MaybeRevertInitialIsotopes(Molecule& m) {
  if (! starting_isotopes) {
    return 0;
  }

  return m.set_isotopes(starting_isotopes.get());
}

int
PerMoleculeData::IdentifyThreeConnectedAromaticSulphur(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() != 16) {
      continue;
    }
    if (a.ncon() != 3) {
      continue;
    }
    if (! m.is_aromatic(i)) {
      continue;
    }

    aromatic_sulphur << i;
  }

  return aromatic_sulphur.size();
}


int
PerMoleculeData::AddConnectedAtomToPositiveAromaticNitrogen(Molecule& m, int uid) const { 
  int rc = 0;

  for (const atom_number_t n : positive_aromatic_nitrogen) {
    //  cerr << "Atom " << j << " +ve arom N, flag " << uid << " val " << in_set[j] << '\n';
    if (uid != in_system[n]) {
      continue;
    }

    const Atom* a = m.atomi(n);
    if (a->ncon() == 2) {
      continue;
    }

    for (const Bond* b : *a) {
      atom_number_t l = b->other(n);

      if (uid == in_system[l]) {
        continue;
      }

      if (m.is_ring_atom(l)) {  // hmmm, is this correct?
        continue;
      }

      in_system[l] = uid;

      rc++;
    }
  }

  return rc;
}

int
PerMoleculeData::ConvertPositiveAroamticNitrogenToNeutral(Molecule& m, int uid) const {
  int rc = 0;

  for (atom_number_t j : positive_aromatic_nitrogen) {
    if (uid != in_system[j]) {
      continue;
    }

    m.set_formal_charge(j, 0);

    rc++;
  }

  return rc;
}

int
all_fragments_too_small(Molecule& m, int min_atoms_per_fragment) {
  int matoms = m.natoms();

  int* fm = new_int(matoms);
  std::unique_ptr<int[]> free_fm(fm);

  m.fragment_membership(fm);

  Set_of_Atoms to_remove;
  to_remove.resize(matoms);

  int nf = m.number_fragments();

  for (int i = 0; i < nf; i++) {
    if (m.atoms_in_fragment(i) < min_atoms_per_fragment) {
      for (int j = 0; j < matoms; j++) {
        if (i == fm[j]) {
          to_remove.add(j);
        }
      }
    }
  }

  if (to_remove.empty()) {
    return 0;
  }

  if (matoms == to_remove.number_elements()) {
    return 1;
  }

  m.remove_atoms(to_remove);

  return 0;
}

static void
WriteSmiles(Molecule& m, IWString_and_File_Descriptor& output) {
  if (write_3d_smiles_if_needed && m.highest_coordinate_dimensionality() == 3) {
    set_append_coordinates_after_each_atom(1);
  }

  output << m.smiles();

  set_append_coordinates_after_each_atom(0);
}

static int
write_ring(Molecule& m, const PerMoleculeData& mdata, int uid,
           IWString_and_File_Descriptor& output) {
  IWString unique_smiles;

#ifdef DEBUG_WRITE_RING
  cerr << "Creating subset id " << uid << " :";
  for (int i = 0; i < m.natoms(); i++) {
    if (uid == in_set[i]) {
      cerr << ' ' << i;
    }
  }
  cerr << '\n';
#endif

  mdata.CreateSubset(m, uid, unique_smiles);

  output << unique_smiles;

  if (use_parent_name_for_rings) {
    output << ' ' << m.name() << 'R';
  }

  output << '\n';

  return output.good();
}

static int
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
do_label_attachment_points_with_environment(Molecule& m, const int* in_system, int uid,
                                            const int* atype) {
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

      int iso;
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
do_label_attachment_points_1(Molecule& m, const int* in_system, int uid) {
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
      } else if (add_double_bonded_neighbors_to_rings) {
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
do_label_attachment_points(Molecule& m, const int* in_system, int uid) {
  m.transform_to_non_isotopic_form();

  if (1 == label_ring_atoms) {
    return do_label_attachment_points_1(m, in_system, uid);
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
PerMoleculeData::IdentifyFunctionalGroups(Molecule& m) {
  functional_group = std::make_unique<uint32_t[]>(m.natoms());
  std::fill_n(functional_group.get(), m.natoms(), 0);

  Molecule_to_Match target(&m);

  int rc = 0;

  for (Substructure_Query* q : singly_bonded_queries) {
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
        if (0 == functional_group[l]) {
          functional_group[l] = p;
        } else {
          functional_group[l] *= p;
        }
      }
    }
  }

  return rc;
}

//  In order to avoid passing around lots of arguments, create a class

class Ring_System_Info {
 private:
  int _rings_in_system;
  int _keep_this_ring_system;
  int _uid_to_assign;
  int* _ring_already_done;

 public:
  Ring_System_Info(int);
  ~Ring_System_Info();

  int ring_already_done(int r) const {
    return _ring_already_done[r];
  }

  int* ring_already_done() {
    return _ring_already_done;
  }

  int uid_to_assign() const {
    return _uid_to_assign;
  }

  void initialise(int u);

  int add_ring(const Ring* r);
  int already_processed(const Ring* r);

  int keep_this_ring_system() const {
    return _keep_this_ring_system;
  }

  void set_keep_this_ring_system(int s) {
    _keep_this_ring_system = s;
  }

  int rings_in_system() const {
    return _rings_in_system;
  }
};

Ring_System_Info::Ring_System_Info(int nr) {
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

  if (r->number_elements() > max_ring_size_to_consider) {
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
assign_atoms_to_ring_system(Molecule& m, Ring_System_Info& rsi, int zring,
                            int* in_system) {
  const Ring* r = m.ringi(zring);

  r->set_vector(in_system, rsi.uid_to_assign());

  rsi.add_ring(r);

  if (add_double_bonded_neighbors_to_rings) {
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
      assign_atoms_to_ring_system(m, rsi, i, in_system);
    }
  }

  return rsi.keep_this_ring_system();
}

static int
MaybeWriteNotFoundRings(Molecule& m,
                        const IWString& unique_smiles) {
  if (! stream_for_missing_ring_systems.is_open()) {
    return 0;
  }

  stream_for_missing_ring_systems << m.smiles() << ' ' << m.name() << '\n';
  stream_for_missing_ring_systems << unique_smiles << '\n';
  stream_for_missing_ring_systems.write_if_buffer_holds_more_than(32768);

  return 1;
}

//#define DEBUG_DO_DATABASE_LOOKUP

static int
do_database_lookup_multiple_dbs(Molecule& m,
                                resizable_array_p<Db>& databases,
                                const IWString& unique_smiles,
                                Ring_System_Uniqueness& rsu) {
  Dbt dkey((void*)unique_smiles.rawchars(), unique_smiles.length());

  // The number of databases in which this if found.
  int ndb_found = 0;
  uint32_t highest_count = 0;
  IWString exemplar_with_highest_count;

  for (auto* database : databases) {
    Dbt fromdb;

    if (0 != database->get(NULL, &dkey, &fromdb, 0)) {  // not in this db.
      continue;
    }

    ++ndb_found;

    google::protobuf::io::ArrayInputStream input(fromdb.get_data(), fromdb.get_size());
    Smi2Rings::Ring ring;
    if (!google::protobuf::TextFormat::Parse(&input, &ring)) {
      cerr << "do_database_lookup_multiple_dbs:invalid db contents '";
      cerr.write(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
      cerr << "'\n";
      return 0;
    }

    if (ring.n() < highest_count) {
      continue;
    }
    highest_count = ring.n();
    exemplar_with_highest_count = ring.ex();
  }

#ifdef DEBUG_DO_DATABASE_LOOKUP
  cerr << "Found in database, dcount " << ndb_found << '\n';
#endif

  rsu.set_uniqueness_measure(highest_count);
  rsu.set_name_from_database(exemplar_with_highest_count);
  rsu.set_ndb_found(ndb_found);

  if (ndb_found == 0) {
    MaybeWriteNotFoundRings(m, unique_smiles);
  }

  return 0;
}

static int
do_database_lookup_single_db(Molecule& m,
                             Db& database, const IWString& unique_smiles,
                             Ring_System_Uniqueness& rsu) {
#ifdef DEBUG_DO_DATABASE_LOOKUP
  cerr << "do_database_lookup_single_db key " << unique_smiles << '\n';
#endif

  Dbt dkey((void*)unique_smiles.rawchars(), unique_smiles.length());

  Dbt fromdb;

  if (0 != database.get(NULL, &dkey, &fromdb, 0)) {  // found new ring
#ifdef DEBUG_DO_DATABASE_LOOKUP
    cerr << "Not found in database '" << unique_smiles << "'\n";
#endif
    rsu.set_uniqueness_measure(0);
    ring_systems_not_found_in_database++;
    MaybeWriteNotFoundRings(m, unique_smiles);
    return 0;
  }

#ifdef DEBUG_DO_DATABASE_LOOKUP
  cerr << "FOUND\n";
#endif

  google::protobuf::io::ArrayInputStream input(fromdb.get_data(), fromdb.get_size());
  Smi2Rings::Ring ring;
  if (!google::protobuf::TextFormat::Parse(&input, &ring)) {
    cerr << "do_database_lookup_single_db:invalid db contents '";
    cerr.write(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
    cerr << "'\n";
    return 0;
  }

#ifdef DEBUG_DO_DATABASE_LOOKUP
  cerr << "Found in database, " << ring.ShortDebugString() << '\n';
#endif

  rsu.set_uniqueness_measure(ring.n());
  rsu.set_name_from_database(ring.ex());

  return 1;
}

static int
do_database_lookup(Molecule& m,
                   PerMoleculeData& mdata,
                   int uid,
                   resizable_array_p<Db>& databases) {
  IWString unique_smiles;
  mdata.CreateSubset(m, uid, unique_smiles);

#ifdef DEBUG_DO_DATABASE_LOOKUP
  cerr << "Unique smiles of subset '" << unique_smiles << "' from " << m.name() << '\n';
#endif

  if (databases.size() == 1) {
    return do_database_lookup_single_db(m, *databases[0], unique_smiles, mdata.rsu[uid]);
  } else {
    return do_database_lookup_multiple_dbs(m, databases, unique_smiles, mdata.rsu[uid]);
  }
}

// These operations are ordered in a specific way, only for doing lookups.
// The most specific lookup is done first.

static int
DoLookups(Molecule& m, int ring_number, PerMoleculeData& mdata,
          resizable_array_p<Db>& databases, IWString_and_File_Descriptor& output) {

#ifdef DEBUG_DO_DATABASE_LOOKUP
  cerr << "DoLookups " << m.smiles() << '\n';
#endif

  if (label_attachment_points_with_environment) {
    do_label_attachment_points_with_environment(m, mdata.in_system, ring_number,
                                                mdata.atype);
    if (write_rings) {
      write_ring(m, mdata, ring_number, output);
    }

    if (do_database_lookup(m, mdata, ring_number, databases)) {
      return 1;
    }

    m.transform_to_non_isotopic_form();
  }

  if (label_attachment_points) {
    do_label_attachment_points(m, mdata.in_system, ring_number);
    if (write_rings) {
      write_ring(m, mdata, ring_number, output);
    }

    if (do_database_lookup(m, mdata, ring_number, databases)) {
      return 1;
    }

    m.transform_to_non_isotopic_form();
  }

  if (always_process_unsubstituted_ring) {
    if (do_database_lookup(m, mdata, ring_number, databases)) {
      return 1;
    }
  }

  return 1;
}

// Print `proto` and place in `destination`.
// Variable `in_scope` is used to hold the newly
// allocated data, and must remain in scope during
// the lifetime of `destination`.
template <typename P>
void
ProtoToDbt(const P& proto, Dbt& destination, std::string& in_scope) {
  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(1);
  printer.PrintToString(proto, &in_scope);

  destination.set_data(in_scope.data());
  destination.set_size(in_scope.size());
}

static int
do_store(Db& database, Dbt& dkey, Dbt& zdata) {
  int rc = database.put(NULL, &dkey, &zdata, 0);
  if (0 == rc) {
    return 1;
  }

  database.err(rc, "Cannot store");
  return 0;
}

static int
store_new_ring(Db& database, Dbt& dkey, const IWString& mname) {
  new_ring_systems_stored++;

  Smi2Rings::Ring ring;
  ring.set_n(1);
  ring.set_ex(mname.AsString());

  Dbt to_store;
  std::string in_scope;
  ProtoToDbt(ring, to_store, in_scope);
  return do_store(database, dkey, to_store);
}

// A ring, `dkey` has been retrieved from the database, with the string
// representation of the proto in `fromdb`.
// Convert `fromdb` to proto form, increment count and write back to `database`.
static int
increment_count(Db& database, Dbt& dkey, Dbt& fromdb) {
  google::protobuf::io::ArrayInputStream input(fromdb.get_data(), fromdb.get_size());
  Smi2Rings::Ring ring;
  if (!google::protobuf::TextFormat::Parse(&input, &ring)) {
    cerr << "increment_count:invalid db contents '";
    cerr.write(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
    cerr << "'\n";
    return 0;
  }

  ring.set_n(ring.n() + 1);

  std::string in_scope;
  Dbt value;
  ProtoToDbt(ring, value, in_scope);
  return do_store(database, dkey, value);
}

static int
DoDatabaseStore(Molecule& m, PerMoleculeData& mdata, int uid,
                Db& database) {
  IWString unique_smiles;
  mdata.CreateSubset(m, uid, unique_smiles);

#ifdef DEBUG_STORE
  cerr << "do_database_store " << unique_smiles << '\n';
#endif
  if (mdata.AlreadySeen(unique_smiles)) {
    return 1;
  }

  Dbt dkey((void*)unique_smiles.rawchars(), unique_smiles.length());
  // dkey.set_data((void *)unique_smiles.rawchars());   // loss of const OK
  // dkey.set_size(unique_smiles.nchars());

  Dbt fromdb;
  fromdb.set_data(NULL);
  fromdb.set_size(0);

  if (0 != database.get(NULL, &dkey, &fromdb, 0)) {
    return store_new_ring(database, dkey, m.name());
  }

  int rc = increment_count(database, dkey, fromdb);

  // delete reinterpret_cast<char *>(fromdb.get_data());

  return rc;
}


static int
DoStore(Molecule& m, int ring_number, PerMoleculeData& mdata,
        Db& database, IWString_and_File_Descriptor& output) {
  if (label_attachment_points_with_environment) {
    do_label_attachment_points_with_environment(m, mdata.in_system, ring_number,
                                                mdata.atype);
    if (write_rings) {
      write_ring(m, mdata, ring_number, output);
    }

    DoDatabaseStore(m, mdata, ring_number, database);

    m.transform_to_non_isotopic_form();
  }

  if (label_attachment_points) {
    do_label_attachment_points(m, mdata.in_system, ring_number);
    if (write_rings) {
      write_ring(m, mdata, ring_number, output);
    }

    DoDatabaseStore(m, mdata, ring_number, database);

    m.transform_to_non_isotopic_form();
  }

  if (always_process_unsubstituted_ring) {
    m.transform_to_non_isotopic_form();

    DoDatabaseStore(m, mdata, ring_number, database);
  }

  return 1;
}

static int
DoLookupsOrStore(Molecule& m, int ring_number, PerMoleculeData& mdata,
                 resizable_array_p<Db>& databases, IWString_and_File_Descriptor& output) {

  if (lookup_in_database) {
    return DoLookups(m, ring_number, mdata, databases, output);
  } else {
    return DoStore(m, ring_number, mdata, *databases[0], output);
  }
}

int
PerMoleculeData::AddImplicitHydrogenToAromaticSulphur(Molecule& m,
                int ring_system) const {
  if (aromatic_sulphur.empty()) {
    return 0;
  }

  int rc = 0;
  for (atom_number_t s : aromatic_sulphur) {
    if (in_system[s] != ring_system) {
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
PerMoleculeData::MaybeConvertPositiveAromaticNitrogen(Molecule& m,
                int ring_system) const {
  if (positive_aromatic_nitrogen.empty()) {
    return 0;
  }

  if (add_connected_atom_to_positive_aromatic_nitrogen) {
    AddConnectedAtomToPositiveAromaticNitrogen(m, ring_system);
  } else if (convert_positive_aromatic_nitrogen_to_neutral) {
     ConvertPositiveAroamticNitrogenToNeutral(m, ring_system);
  }

  return 1;
}

// Update various file scope statistics

int
PerMoleculeData::MaybeAddOutsideRingAtoms(Molecule& m,
                         int ring_system_number) {
  if (! functional_group) {
    return 0;
  }

  return AddExtraRingSingleBondedNeighgours(m, ring_system_number);
}


int
PerMoleculeData::IdentifyPositiveAromaticNitrogen(Molecule& m) {
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
      positive_aromatic_nitrogen.add(i);
    }
  }

  return positive_aromatic_nitrogen.number_elements();
}


//  Create the unique smiles for a given subset of the molecule

int
PerMoleculeData::CreateSubset(Molecule& m, int uid,
                IWString& unique_smiles) const {
  Molecule tmp;

  m.create_subset(tmp, in_system, uid, xref.get());

  int matoms = tmp.natoms();
  for (int i = 0; i < matoms; i++) {
    tmp.set_implicit_hydrogens_known(i, 0);

    if (0 == tmp.isotope(i)) {
      continue;
    }

    tmp.recompute_implicit_hydrogens(i);
  }

  for (atom_number_t s : aromatic_sulphur) {
    if (xref[s] >= 0) {
      tmp.set_implicit_hydrogens(xref[s], 1, 1);
    }
  }

// #define DEBUG_SUBSET_CREATION
#ifdef DEBUG_SUBSET_CREATION
  tmp.debug_print(cerr);
  cerr << "Smiles for subset '" << tmp.smiles() << "' unique '";
  cerr << tmp.unique_smiles() << "'\n";
#endif

  if (remove_invalid_chiral_centres) {
    lillymol::do_remove_invalid_chiral_centres(subset);
  }

  if (write_3d_smiles_if_needed && tmp.highest_coordinate_dimensionality() == 3) {
    set_append_coordinates_after_each_atom(1);
  }

  unique_smiles = tmp.unique_smiles();
  set_append_coordinates_after_each_atom(0);

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
PerMoleculeData::AddExtraRingSingleBondedNeighgours(Molecule& m, int uid) {
  int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++) {
    if (uid != in_system[i]) {  // just examine our subset
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
      if (uid == in_system[k]) {  // we are looking for atoms outside the ring system
        continue;
      }

      // adding all singly bonded attachments
      if (2 == add_single_bonded_neighbours_to_rings) {
        in_system[k] = kTmpid;
        if (functional_group[i] && functional_group[k]) {
          add_functional_group(matoms, i, functional_group.get(), in_system, kTmpid);
        }

        rc++;
      } else if (functional_group[i] &&
                 functional_group[k])  // just adding extra's which are part of a
                                       // functional group
      {
        in_system[k] = kTmpid;
        add_functional_group(matoms, i, functional_group.get(), in_system, kTmpid);
        rc++;
      } else if (add_all_single_atom_ring_substituents && 1 == m.ncon(k)) {
        in_system[k] = kTmpid;
        rc++;
      }

      if (continue_chain_neighbours_through_multiple_bonds && 0 == m.nrings(k)) {
        Atom* ak = const_cast<Atom*>(m.atomi(k));

        if (ak->nbonds() > ak->ncon()) {
          do_continue_chain_neighbours_through_multiple_bonds(m, i, k, in_system, kTmpid);
        }
        rc++;
      }
    }
  }

  if (rc) {
    for (int i = 0; i < matoms; i++) {
      //    cerr << "Atom " << i << " set " << in_system[i] << '\n';
      if (kTmpid == in_system[i]) {
        in_system[i] = uid;
      }
    }
  }

  // a zero return code is OK - just means no extra-ring double bonds found
  return rc;
}

bool
PerMoleculeData::AlreadySeen(const IWString& unique_smiles) {
  if (found_this_molecule.contains(unique_smiles)) {
    return true;
  }

  found_this_molecule.insert(unique_smiles);

  return false;
}

Smi2Rings::Smi2Rings() {
}

int
Smi2Rings::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  return 1;
}

}  // namespace mol2rings
