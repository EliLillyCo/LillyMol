/*
  Produces unique smiles for ring systems within a molecule
  Can store these ring systems in a Berkeley database, or can lookup
  ring systems in a Berkeley database
*/

#include <math.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/primes.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "db_cxx.h"

#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools_Bdb/smi2rings.pb.h"

using std::cerr;

static Chemical_Standardisation chemical_standardisation;

static Atom_Typing_Specification atom_typing;

const char* prog_name = nullptr;

static int verbose = 0;

// If working as a filter, the minimum number of exemplars required.
static int min_urs_needed = 0;

static int discarded_for_urs_too_low = 0;

// We can optionally add ALL exocyclic double bonds.
static int add_double_bonded_neighbors_to_rings = 0;

static int add_single_bonded_neighbours_to_rings = 0;

// This will add CH3, OH, NH2, etc groups

static int add_all_single_atom_ring_substituents = 0;

static int continue_chain_neighbours_through_multiple_bonds = 0;

static int label_attachment_points = 0;

// If we are processing isotopically labelled attachment points, we can also
// process the unsubstuted ring.
static int always_process_unsubstituted_ring = 0;

static int label_ring_atoms = 0;

static int label_attachment_points_with_environment = 0;

static int report_all_ring_uniqueness_measurements_looked_up = 1;

/*
  In talking with Michal Vieth, he wants to capture functional groups
*/

static resizable_array_p<Substructure_Query> singly_bonded_queries;

static int max_ring_size_to_consider = std::numeric_limits<int>::max();

static int discard_molecules_where_ring_is_whole_molecule = 0;

static int molecules_read = 0;

static resizable_array_p<Db> databases;

static int new_ring_systems_stored = 0;

static int lookup_in_database = 0;
static int store_in_database = 0;

static int molecules_with_unique_ring_systems_found = 0;

static int ring_systems_not_found_in_database = 0;

static int write_rings = 1;

// Rather than storing rings in the database as they are found, we can
// cache the results and then copy that hash to the database when done.
// This is very desirable for performance.
static int hash_store = 0;

static IWString smiles_tag("$SMI<");

static int write_parent_structure = 0;

static int use_parent_name_for_rings = 0;

static int write_non_ring_atoms = 0;

static Element_Transformations element_transformations;

static int max_ring_system_size_to_consider = std::numeric_limits<int>::max();

//  Various things we can do if we have positively charged aromatic Nitrogens

static int convert_positive_aromatic_nitrogen_to_neutral = 0;
static int add_connected_atom_to_positive_aromatic_nitrogen = 0;

//  Various statistics on the ring systems found

static Report_Progress report_progress;

// The text appended to each molecule to give the count of the lowest count ring system

static IWString urs_string(" URS:");

//  A count of the number of times a ring system with a given number of
//  component rings is found

static extending_resizable_array<int> ring_system_size_array;

//  An accumulator of the number of rings in the systems

static Accumulator_Int<uint64_t> ring_system_ring_stats;

//  An accumulator of the number of atoms in the systems

static Accumulator_Int<uint64_t> ring_system_atom_stats;

// Accumulator on the number of exemplars.
static Accumulator_Int<uint64_t> acc_exemplar;
// And for low numbers, we keep track of the number of instances of
// each number of exemplars.
static int max_exemplar_count = 20;

static extending_resizable_array<uint64_t> exemplar_count;

static int bridge_across_spiro_groups = 0;

static int remove_all_chiral_centres = 0;

static int remove_invalid_chiral_centres = 0;

static IWString_and_File_Descriptor stream_for_rejected_molecules;

static int write_3d_smiles_if_needed = 0;

// We can do two kinds of proto output.
// Note that currently we are not filling in the ring_type attribute.
// Not sure it is necessary.
// smiles id textproto
static int smiles_and_textproto = 0;
// smiles and id are part of the proto.
static int full_proto_output = 0;

static IWString_and_File_Descriptor stream_for_missing_ring_systems;

// The programme attempts to re-use code for both the store
// and lookup tasks as common as possible.
// WHile it would be better to refactor this into a library,
// we can make things a little smoother with a class that
// holds all the output information.
class Output {
  private:
    // If we are writing all rings to an output stream.
    IWString_and_File_Descriptor* _text_output = nullptr;

    // If doing lookups there will be one or more databases being
    // queried.
    // If we are doing a store, there will be a single database that
    // will be updated.
    resizable_array_p<Db> _database;

    // As an optimisation during database loads, we can cache all the rings
    // in memory and write them when done.
    // Output will be to the first (and only) database in `_database`.
    int _accumulate_rings = 0;
    absl::flat_hash_map<std::string, Smi2Rings::Ring> _hash;

  public:
    Output();
    ~Output();

    void set_accumulate_rings(int s) {
      _accumulate_rings = s;
    }
};

/*
  When we are applying isotopic labels to the rings, we need to
  keep track of the uniqueness of each ring system
*/

class Ring_System_Uniqueness {
 private:
  int _uid;
  int _uniqueness_measure;

  //  If we are looking up in a database with names, we can return
  //  the name of the molecule that first saw that ring in the database

  IWString _name_from_database;

  // If we are doing lookups across multiple databases, we keep track
  // of the number of databases in which a ring is found.
  int _ndb_found;

 public:
  Ring_System_Uniqueness();

  void set_uid(int s) {
    _uid = s;
  }

  int uid() const {
    return _uid;
  }

  void set_uniqueness_measure(int s) {
    _uniqueness_measure = s;
  }

  int uniqueness_measure() const {
    return _uniqueness_measure;
  }

  void set_name_from_database(const const_IWSubstring& s) {
    _name_from_database = s;
  }

  const IWString& name_from_database() const {
    return _name_from_database;
  }

  void set_ndb_found(int s) {
    _ndb_found = s;
  }

  int ndb_found() const {
    return _ndb_found;
  }
};

Ring_System_Uniqueness::Ring_System_Uniqueness() {
  _uid = -1;
  _uniqueness_measure = std::numeric_limits<int>::max();
  _ndb_found = 0;
}

class Ring_System_Uniqueness_Comparitor {
 private:
 public:
  int operator()(const Ring_System_Uniqueness&, const Ring_System_Uniqueness&) const;
};

int
Ring_System_Uniqueness_Comparitor::operator()(const Ring_System_Uniqueness& rs1,
                                              const Ring_System_Uniqueness& rs2) const {
  int n1 = rs1.uniqueness_measure();
  int n2 = rs2.uniqueness_measure();

  if (n1 < n2) {
    return -1;
  }

  if (n1 > n2) {
    return 1;
  }

  return 0;
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
  cerr << R"(Identifies rings and ring systems. Stores or retrieves by unique smiles of the ring.
Data is stored in a BerkeleyDb database with key unique smiles and value a textproto
Smi2Rings::Ring proto.
smi2rings_bdb -d dbname -d STORE|LOOKUP <options> file.smi
  -j ...         options for adding atoms to rings, enter '-j help' for info
  -N add         add non-ring connection to positive aromatic nitrogen
  -N neut        convert positive aromatic nitrogen to neutral form
  -y             put isotopes on ring atoms to indicate ring system membership
  -z <size>      discard rings or ring systems with a ring of size > <size>
  -R <max>       max ring system size (number rings) to consider
  -d <dbname>    specify database for lookup or storage
  -d STORE       store ring systems found in the database
  -n             suppress output if writing to a database
  -d LOOKUP      lookup each ring system in the database
  -U <string>    string to append with the lowest ring count (default 'URS')
  -p w           write the parent structure - all atoms
  -p n           use the name of the parent for the name of the rings
  -s <min>       write non-ring fragments with at least <min> atoms
  -f <cutoff>    function as filter (lookup), discard molecules below <cutoff> URS
                 -f file=<fname> to write rejected molecules to <fname>
  -S <query>     one or more queries for functional groups. Typically single atom ring substituents
  -r <number>    report progress every <number> molecules processed
  -a             only report rarest in lookups
  -e             discard rings where the ring is the entire molecule
  -c             remove all chiral centres
  -Y ...         obscure options, enter '-Y help' for info
  -E ...         standard element options, enter '-E help' for info
  -A ...         standard aromaticity options
  -g ...         chemical standardisation options
  -v             verbose output
Typical usage
  smi2rings_bdb -d rings.bdb -d STORE -j ring -j iso -j env -j double -N add -n -c collection.smi
  smi2rings_bdb -d rings.bdb -d LOOKUP -n -f 10 -v -c -j double -N add -j iso query.smi
)";
// clang-format on

  exit(rc);
}

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

struct PerMoleculeData {
  Molecule& m;
  int* in_system;
  uint32_t* atype;

  // For each ring system discovered, we need to know the number of rings in
  // it. This array will only be partially filled if there are fused rings.
  int* ring_system_size;

  // For each ring system found, a measure of uniqueness.
  Ring_System_Uniqueness* rsu;

  // the unique smiles of the rings identified in this molecule.
  IW_STL_Hash_Set found_this_molecule;

  // If the molecule has isotopes, we capture them here.
  std::unique_ptr<isotope_t[]> starting_isotopes;

  // If we are adding functional groups.
  std::unique_ptr<isotope_t[]> functional_group;

  // We may need to adjust the smiles written for certain atom types.
  Set_of_Atoms positive_aromatic_nitrogen;
  // Three connected only;
  Set_of_Atoms aromatic_sulphur;

  // create_subset needs a cross reference array.
  std::unique_ptr<int[]> xref;

 public:
  PerMoleculeData(Molecule& m);
  ~PerMoleculeData();

  int natoms() const {
    return m.natoms();
  }

  int AssignAtomTypes(Atom_Typing_Specification& atom_typing);

  int IdentifyThreeConnectedAromaticSulphur(Molecule& m);
  int IdentifyPositiveAromaticNitrogen(Molecule& m);
  int IdentifyFunctionalGroups(Molecule& m);

  int MaybeAddOutsideRingAtoms(Molecule& m,
                         int ring_system_number);
  int MaybeConvertPositiveAromaticNitrogen(Molecule& m,
                int ring_system) const;

  int AddExtraRingSingleBondedNeighgours(Molecule& m, int uid);
  int AddConnectedAtomToPositiveAromaticNitrogen(Molecule& m, int uid) const;
  int ConvertPositiveAroamticNitrogenToNeutral(Molecule& m, int uid) const;
  int AddImplicitHydrogenToAromaticSulphur(Molecule& m, int ring_system) const;

  int CreateSubset(Molecule& m, int uid,
                IWString& unique_smiles) const;

  // Return true if found_this_molecule contains m.unique_smiles().
  bool AlreadySeen(const IWString& unique_smiles);

  // If starting_isotopes is set, set those isotopes in `m`.
  int MaybeRevertInitialIsotopes(Molecule& m);
};

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
  atype = new uint32_t[m.natoms()];
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

static int
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

static int
do_write_non_ring_atoms(Molecule& m, int* in_system,
                        IWString_and_File_Descriptor& output) {
  int matoms = m.natoms();
  int atoms_in_subset = 0;
  for (int i = 0; i < matoms; i++) {
    if (in_system[i])  // was a ring atom, don't want those
    {
      in_system[i] = 0;
    } else {
      in_system[i] = 1;
      atoms_in_subset++;
    }
  }

  if (atoms_in_subset < write_non_ring_atoms) {
    return 1;
  }

  Molecule subset;
  m.create_subset(subset, in_system, 1);

  if (subset.number_fragments() > 1) {
    if (all_fragments_too_small(subset, write_non_ring_atoms)) {
      return 1;
    }
  }

  if (remove_invalid_chiral_centres) {
    lillymol::do_remove_invalid_chiral_centres(subset);
  }

  output << subset.smiles() << ' ' << m.name() << " non ring\n";

  return output.good();
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
assign_atoms_to_ring_systems(Molecule& m, int* in_system, int* ring_system_size) {
  int matoms = m.natoms();

  int nr = m.nrings();

  Ring_System_Info rsi(nr);

  int uid_to_assign = 0;

  for (int i = 0; i < nr; i++) {
    if (rsi.ring_already_done(i)) {
      continue;
    }

    rsi.initialise(uid_to_assign);

    if (assign_atoms_to_ring_system(m, rsi, i, in_system)) {
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

static Ring_System_Uniqueness_Comparitor ring_system_uniqueness_comparitor;

// Molecule `m` has been disarded for having a lowest ring count value
// too low. Increment global counter and maybe write.
static int
HandleDiscardedForMinUrs(Molecule& m, int lowest_ring_count) {
  discarded_for_urs_too_low++;
  if (stream_for_rejected_molecules.active()) {
    stream_for_rejected_molecules << m.smiles() << ' ' << m.name() << urs_string
                                  << lowest_ring_count << '\n';
    stream_for_rejected_molecules.write_if_buffer_holds_more_than(4096);
  }

  return 1;
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
StoreInHash(Molecule& m,
            const IWString& unique_smiles,
            std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash) {
  auto iter = dbhash->find(unique_smiles);
  if (iter != dbhash->end()) {
    auto n = iter->second.n();
    iter->second.set_n(n + 1);
    return 1;
  }

  ++new_ring_systems_stored;

  Smi2Rings::Ring proto;
  proto.set_n(1);
  proto.set_ex(m.name().data(), m.name().length());

  dbhash->emplace(std::make_pair(unique_smiles, std::move(proto)));

  return 1;
}
static int
DoDatabaseStore(Molecule& m, PerMoleculeData& mdata, int uid,
                std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
                Db& database) {
  IWString unique_smiles;
  mdata.CreateSubset(m, uid, unique_smiles);

#ifdef DEBUG_STORE
  cerr << "do_database_store " << unique_smiles << '\n';
#endif
  if (mdata.AlreadySeen(unique_smiles)) {
    return 1;
  }

  if (dbhash) {
    return StoreInHash(m, unique_smiles, dbhash);
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
        Db& database,
        std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
        IWString_and_File_Descriptor& output) {
  if (label_attachment_points_with_environment) {
    do_label_attachment_points_with_environment(m, mdata.in_system, ring_number,
                                                mdata.atype);
    if (write_rings) {
      write_ring(m, mdata, ring_number, output);
    }

    DoDatabaseStore(m, mdata, ring_number, dbhash, database);

    m.transform_to_non_isotopic_form();
  }

  if (label_attachment_points) {
    do_label_attachment_points(m, mdata.in_system, ring_number);
    if (write_rings) {
      write_ring(m, mdata, ring_number, output);
    }

    DoDatabaseStore(m, mdata, ring_number, dbhash, database);

    m.transform_to_non_isotopic_form();
  }

  if (always_process_unsubstituted_ring) {
    m.transform_to_non_isotopic_form();

    DoDatabaseStore(m, mdata, ring_number, dbhash, database);
  }

  return 1;
}

static int
DoLookupsOrStore(Molecule& m, int ring_number, PerMoleculeData& mdata,
                 resizable_array_p<Db>& databases,
                 std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
                 IWString_and_File_Descriptor& output) {

  if (lookup_in_database) {
    return DoLookups(m, ring_number, mdata, databases, output);
  } else {
    return DoStore(m, ring_number, mdata, *databases[0], dbhash, output);
  }
}

static int
WriteProto(Molecule& m, const PerMoleculeData& mdata,
           int number_ring_systems,
           IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  Smi2Rings::Results results;
  int nprocess = report_all_ring_uniqueness_measurements_looked_up ?
        number_ring_systems : 1;
  for (int i = 0; i < nprocess; ++i) {
    Smi2Rings::Result* r = results.add_result();
    const auto& u = mdata.rsu[i];
    r->set_urs(u.uniqueness_measure());
    if (u.uniqueness_measure() > 0) {
      const IWString& s = u.name_from_database();
      r->set_ex(s.data(), s.length());
    }
  }

  if (full_proto_output) {
    const IWString& s = m.smiles();
    results.set_smiles(s.data(), s.length());
    const IWString& n = m.name();
    results.set_id(n.data(), n.length());
  } else {
    output << m.smiles() << ' ' << m.name() << ' ';
  }

  std::string buffer;
  printer.PrintToString(results, &buffer);

  output << buffer << '\n';

  return 1;
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

static void
UpdateGlobalStatistics(const PerMoleculeData& mdata,
                int ring_system_number) {
  if (verbose == 0) {
    return;
  }

  int rings_in_system = mdata.ring_system_size[ring_system_number];

  ring_system_size_array[rings_in_system]++;
  ring_system_ring_stats.extra(rings_in_system);

  int atoms_in_system = std::count(mdata.in_system, mdata.in_system + mdata.natoms(),
                ring_system_number);

  ring_system_atom_stats.extra(atoms_in_system);
}

int
PerMoleculeData::MaybeAddOutsideRingAtoms(Molecule& m,
                         int ring_system_number) {
  if (! functional_group) {
    return 0;
  }

  return AddExtraRingSingleBondedNeighgours(m, ring_system_number);
}


static int
smi2rings(Molecule& m, PerMoleculeData& mdata,
          std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
          IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  // The number of ring systems.
  const int n = assign_atoms_to_ring_systems(m, mdata.in_system, mdata.ring_system_size);

  for (int i = 0; i < n; i++) {
    mdata.rsu[i].set_uid(i);

    UpdateGlobalStatistics(mdata, i);

    if (mdata.ring_system_size[i] > max_ring_system_size_to_consider) {
      continue;
    }

    mdata.MaybeConvertPositiveAromaticNitrogen(m, i);
    mdata.AddImplicitHydrogenToAromaticSulphur(m, i);

    mdata.MaybeAddOutsideRingAtoms(m, i);

    if (discard_molecules_where_ring_is_whole_molecule &&
        AllValuesAre(mdata.in_system, matoms, i)) {
      continue;
    }

    DoLookupsOrStore(m, i, mdata, databases, dbhash, output);
  }

  if (write_non_ring_atoms) {
    do_write_non_ring_atoms(m, mdata.in_system, output);
  }

  if (!lookup_in_database) {
    return output.good();
  }

  // If we were doing a lookup, check the results

  if (label_attachment_points || label_attachment_points_with_environment ||
      always_process_unsubstituted_ring) {
    m.transform_to_non_isotopic_form();
  }

  if (n > 1) {
    iwqsort(mdata.rsu, n, ring_system_uniqueness_comparitor);
  }

#ifdef ECHO_SORTED_LIST
  if (i > 1) {
    for (int i = 0; i < n; i++) {
      cerr << i << ' ' << rsu[i].uniqueness_measure() << '\n';
    }
  }
#endif

  int lowest_ring_count = mdata.rsu[0].uniqueness_measure();

  // none of the ring systems were looked up
  if (std::numeric_limits<int>::max() == lowest_ring_count) { 
    return output.good();
  }

  if (verbose > 1) {
    cerr << m.name() << " of " << n << " ring systems, least common ring system count "
         << lowest_ring_count << '\n';
  }

  acc_exemplar.extra(lowest_ring_count);
  if (lowest_ring_count <= max_exemplar_count) {
    ++exemplar_count[lowest_ring_count];
  }

  if (0 == lowest_ring_count) {
    molecules_with_unique_ring_systems_found++;
  }

  if (lowest_ring_count < min_urs_needed) {
    return HandleDiscardedForMinUrs(m, lowest_ring_count);
  }

  if (smiles_and_textproto || full_proto_output) {
    return WriteProto(m, mdata, n, output);
  }

  IWString tmp(m.name());

  tmp << urs_string << mdata.rsu[0].uniqueness_measure();
  if (0 == mdata.rsu[0].uniqueness_measure()) {
    tmp << " .";
  } else {
    tmp << ' ' << mdata.rsu[0].name_from_database();
  }

  if (report_all_ring_uniqueness_measurements_looked_up) {
    for (int i = 1; i < n; i++) {
      tmp << ' ' << mdata.rsu[i].uniqueness_measure();

      if (i != (n - 1)) {
        tmp << ',';
      }
    }
  }

#ifdef DEBUG_WRITE_RING
  cerr << label_ring_atoms << " label_ring_atoms, n = " << n << '\n';
  for (int i = 0; i < matoms; ++i) {
    cerr << i << " in_system " << mdata.in_system[i] << ' '
         << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

  mdata.MaybeRevertInitialIsotopes(m);

  if (label_ring_atoms) {
    int uid = mdata.rsu[0].uid();
    for (int i = 0; i < matoms; ++i) {
      if (mdata.in_system[i] == uid) {
        m.set_isotope(i, 1);
      }
    }
    IWString name(m.name());
    name << urs_string << mdata.rsu[0].uniqueness_measure();
    m.set_name(name);
#ifdef OLD_VERSION_THAT_IS_NOT_VERY_USEFUL
    for (int i = 0; i < n; i++) {
      const int u = mdata.rsu[i].uid();

      for (int j = 0; j < matoms; j++) {
        if (u == mdata.in_system[j]) {
          m.set_isotope(j, i + 1);
        }
      }
    }
#endif
  }

  WriteSmiles(m, output);
  output << ' ' << tmp << '\n';

  return output.good();
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
    lillymol::do_remove_invalid_chiral_centres(tmp);
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

static int
smi2rings(Molecule& m,
          std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
          IWString_and_File_Descriptor& output) {
  if (report_progress()) {
    cerr << "Read " << molecules_read << " molecules, found " << new_ring_systems_stored
         << " new ring systems\n";
  }

  if (write_parent_structure || write_non_ring_atoms) {
    WriteSmiles(m, output);
    output << ' ' << m.name() << '\n';
  }

  int nr = m.nrings();

  if (0 == nr) {  // that's OK, we just have nothing to do
    return 1;
  }

  PerMoleculeData mdata(m);
  if (atom_typing.active()) {
    mdata.AssignAtomTypes(atom_typing);
  }

  if (singly_bonded_queries.number_elements()) {
    mdata.IdentifyFunctionalGroups(m);
  }

  if (add_connected_atom_to_positive_aromatic_nitrogen ||
      convert_positive_aromatic_nitrogen_to_neutral) {
    mdata.IdentifyPositiveAromaticNitrogen(m);
  }

  mdata.IdentifyThreeConnectedAromaticSulphur(m);

  return smi2rings(m, mdata, dbhash, output);
}

/*
  Until we get directional bonds in create_subset working, just strip
  directional bonds. Besides there aren't any in rings
*/

static void
Preprocess(Molecule& m) {
  if (element_transformations.active()) {
    (void)element_transformations.process(&m);
  }

  if (chemical_standardisation.active()) {
    (void)chemical_standardisation.process(m);
  }

  if (remove_all_chiral_centres) {
    m.remove_all_chiral_centres();
  }

  m.revert_all_directional_bonds_to_non_directional();

  return;
}

static int
smi2rings(data_source_and_type<Molecule>& input,
          std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
          IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule()) && output.good()) {
    std::unique_ptr<Molecule> free_m(m);

    Preprocess(*m);

    molecules_read++;

    if (!smi2rings(*m, dbhash, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(16384);
  }

  return output.good();
}

static int
smi2rings(const char* fname, FileType input_type,
          std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>>& dbhash,
          IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return smi2rings(input, dbhash, output);
}

static int
DoStore(const IWString& usmi,
        const Smi2Rings::Ring& proto,
        Db& database) {
  Dbt dbkey((void*)usmi.data(), usmi.length());

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  printer.PrintToString(proto, &buffer);

  Dbt dbdata((void*)buffer.data(), buffer.length());

  return do_store(database, dbkey, dbdata);
}

static int
WriteHash(absl::flat_hash_map<IWString, Smi2Rings::Ring>& dbhash,
          Db* database) {
  for (const auto& [usmi, proto] : dbhash) {
    DoStore(usmi, proto, *database);
  }

  return 1;
}

static void
display_dash_j_options(std::ostream& os) {
  // clang-format off
  os << R"(
 -j double      add doubly bonded neighbours to rings - always do this.
 -j all         add all singly bonded neighbours to rings
 -j fg          add only singly bonded functional groups to rings
 -j unsat       add all atoms bonded to an unsaturated neighbour
 -j spiro       include spiro fusions in ring systems
 -j safg        include all single atom functional groups, CH3, NH2, OH, F etc
 -j iso         isotopically label the attachment points
 -j isocon      isotope number is number of outside-system attachments
 -j env         isotope number reflects environment of attachment
                use '-j env=UST:achry' to specify an atom type for the env
 -j ring        always process the bare ring - in addition to isotopic forms
When building a database (-d STORE) use all of '-j ring', '-j iso' and '-j env'.
When doing lookups use one or more of them.
)";
  // clang-format on

  exit(0);
}

static void
DisplayDashYOptions(std::ostream& output) {
  // clang-format off
  cerr << R"( -Y smi3d        write 3D smiles when 3D information present
 -Y smiproto     write smiles, id followed by textproto output
 -Y proto        write as Smi2Rings::Results textproto
 -Y writemissing=<fname> during lookup write any missing rings to <fname>.smi
 -Y hash         when storing, hash all results and write database at end (recommended)
)";
  // clang-format on

  usage(1);
}

static int
smi2rings(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:T:i:d:j:D:g:nN:s:R:f:z:S:r:p:U:uycCaeY:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "The report progress option (-r) must be whole positive number\n";
      usage(5);
    }
  }

  if (cl.option_present('U')) {
    cl.value('U', urs_string);

    if (verbose) {
      cerr << "URS string set to '" << urs_string << "'\n";
    }

    if (urs_string.empty()) {
      cerr << "Warning, zero length Unique Ring System string\n";
    }
  }

  if (cl.option_present('T')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T')) {
      cerr << "Cannot parse -T option\n";
      usage(11);
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      usage(6);
    }
  }

  if (cl.option_present('j')) {
    int i = 0;
    const_IWSubstring j;
    while (cl.value('j', j, i++)) {
      if ("all" == j) {
        add_single_bonded_neighbours_to_rings = 2;
        if (verbose) {
          cerr << "Will include singly bonded atoms outside the ring system\n";
        }
      } else if ("iso" == j) {
        label_attachment_points = 1;
        if (verbose) {
          cerr << "Isotopic labels applied at attachment points\n";
        }
      } else if ("isocon" == j) {
        label_attachment_points = 2;
        if (verbose) {
          cerr << "The isotopic labels applied will be the number of connections\n";
        }
      } else if ("spiro" == j) {
        bridge_across_spiro_groups = 1;
        if (verbose) {
          cerr << "Will bridge across spiro fusions for ring systems\n";
        }
      } else if ("fg" == j) {
        add_single_bonded_neighbours_to_rings = 1;
        if (verbose) {
          cerr << "Will add singly bonded, identified functional groups to the ring "
                  "system\n";
        }
      } else if ("double" == j) {
        add_double_bonded_neighbors_to_rings = 1;
        if (verbose) {
          cerr << "Multiple bonds outside the ring system will be included\n";
        }
      } else if ("unsat" == j) {
        continue_chain_neighbours_through_multiple_bonds = 1;
        if (verbose) {
          cerr << "Will continue chain neighbours through multiple bonds\n";
        }
      } else if ("glaxo" == j) {
        label_attachment_points = 1;
        add_double_bonded_neighbors_to_rings = 1;
      } else if ("safg" == j) {
        add_all_single_atom_ring_substituents = 1;
        if (verbose) {
          cerr << "Will add all single atom functional groups\n";
        }
      } else if ("env" == j) {
        label_attachment_points_with_environment = 1;
        if (verbose) {
          cerr << "isotopic labels will reflect environment\n";
        }
      } else if (j.starts_with("env=")) {
        j.remove_leading_chars(4);

        if (!atom_typing.build(j)) {
          cerr << "Cannot determing atom typing for environment '" << j << "'\n";
          return 3;
        }
        label_attachment_points_with_environment = 1;
        if (verbose) {
          cerr << "Isotopic labels will reflect environment with type '" << j << "'\n";
        }
      } else if (j == "ring") {
        always_process_unsubstituted_ring = 1;
        if (verbose) {
          cerr << "Will always examine the non-isotopically labelled ring\n";
        }
      } else if ("help" == j) {
        display_dash_j_options(cerr);
      } else {
        cerr << "Unrecognised -j qualifier '" << j << "'\n";
        display_dash_j_options(cerr);
      }
    }
  }

  if (label_attachment_points_with_environment) {
  } else if (label_attachment_points) {
  } else if (always_process_unsubstituted_ring) {
  } else {
    cerr << "Default action, store labelled attachment points\n";
    label_attachment_points = 1;
  }

  if (cl.option_present('S')) {
    if (!process_queries(cl, singly_bonded_queries, verbose, 'S')) {
      cerr << "Cannot parse functional group queries\n";
      usage(4);
    }

    add_single_bonded_neighbours_to_rings = 1;

    if (verbose) {
      cerr << "Defined " << singly_bonded_queries.number_elements()
           << " queries for singly bonded functional groups\n";
    }

    for (int i = 0; i < singly_bonded_queries.number_elements(); i++) {
      singly_bonded_queries[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('a')) {
    report_all_ring_uniqueness_measurements_looked_up = 0;
    if (verbose) {
      cerr << "Will only report the rarest ring system during lookups\n";
    }
  }

  if (cl.option_present('z')) {
    if (!cl.value('z', max_ring_size_to_consider) || max_ring_size_to_consider < 3) {
      cerr << "The max ring size to consider option (-z) must be followed by a whole "
              "number > 2\n";
      usage(41);
    }

    if (verbose) {
      cerr << "Will not consider rings larger than " << max_ring_size_to_consider << '\n';
    }
  }

  if (cl.option_present('e')) {
    discard_molecules_where_ring_is_whole_molecule = 1;

    if (verbose) {
      cerr << "Discard molecules where the ring is the whole molecule\n";
    }
  }

  if (! cl.option_present('d')) {
    cerr << "Must specify database(s) to process via the -d option\n";
    usage(1);
  }

  if (cl.option_present('d')) {
    int oflags = 0;
    DBTYPE dbtype = DB_UNKNOWN;

    resizable_array<IWString> dbnames;

    IWString d;
    for (int i = 0; cl.value('d', d, i); ++i) {
      if ("STORE" == d) {
        oflags = DB_CREATE;
        dbtype = DB_BTREE;
        store_in_database = 1;
      } else if (d.starts_with("LOOK")) {
        oflags = DB_RDONLY;
        dbtype = DB_UNKNOWN;
        lookup_in_database = 1;
      } else if (d == "NAME") {
        // obsolete option now ignored.
      } else {
        dbnames << d;
      }
    }

    if (lookup_in_database && store_in_database) {
      cerr << "Cannot specify both store and lookup operations\n";
      usage(4);
    }

    if (lookup_in_database && !cl.option_present('n')) {
      cerr << "When doing database lookups, must add the -n option\n";
      usage(17);
    }

    if (0 == store_in_database && 0 == lookup_in_database) {
      cerr << "Must specify either '-d STORE' or '-d LOOKUP'\n";
      usage(17);
    }

    databases.resize(dbnames.number_elements());
    for (IWString& dbname : dbnames) {
      std::unique_ptr<Db> db(new Db(NULL, DB_CXX_NO_EXCEPTIONS));
      if (int rc =
              db->open(NULL, dbname.null_terminated_chars(), 0, dbtype, oflags, 0644);
          rc != 0) {
        cerr << "Cannot open database '" << dbname << "'\n";
        db->err(rc, "");
        return 1;
      }
      if (verbose) {
        cerr << "Opened database '" << dbname << "'\n";
      }
      databases << db.release();
    }
  }

  if (cl.option_present('y')) {
    if (!lookup_in_database) {
      cerr << "The label atoms by ring system option (-y) only makes sense when doing "
              "lookups\n";
      usage(5);
    }

    label_ring_atoms = 1;

    if (verbose) {
      cerr << "Atoms labelled by ring system membership\n";
    }
  }

  if (cl.option_present('u')) {
    cerr << "The -u option is now ignored\n";
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', max_ring_system_size_to_consider) ||
        max_ring_system_size_to_consider < 1) {
      cerr << "The -R option must be followed by a whole positive number\n";
      usage(39);
    }

    if (verbose) {
      cerr << "Ring systems with more than " << max_ring_system_size_to_consider
           << " rings will be ignored\n";
    }
  }

  if (cl.option_present('n')) {
    write_rings = 0;
    if (verbose) {
      cerr << "Default output to stdout suppressed\n";
    }
  }

  if (cl.option_present('N')) {
    int i = 0;
    const_IWSubstring n;
    while (cl.value('N', n, i++)) {
      if ("add" == n) {
        add_connected_atom_to_positive_aromatic_nitrogen = 1;
        if (verbose) {
          cerr << "Will add connected atoms to positive aromatic nitrogens\n";
        }
      } else if ("neut" == n) {
        convert_positive_aromatic_nitrogen_to_neutral = 1;
        if (verbose) {
          cerr << "Will convert positive aromatic nitrogen to neutral form\n";
        }
      } else {
        cerr << "Unrecognised -N qualifier '" << n << "'\n";
        usage(5);
      }
    }
  }

  if (cl.option_present('p')) {
    int i = 0;
    const_IWSubstring p;

    while (cl.value('p', p, i++)) {
      if ('w' == p) {
        write_parent_structure = 1;

        if (verbose) {
          cerr << "Will write the parent structure\n";
        }
      } else if ('n' == p) {
        use_parent_name_for_rings = 1;

        if (verbose) {
          cerr << "Each ring inherits the name of its parent structure\n";
        }
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(5);
      }
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', write_non_ring_atoms) || write_non_ring_atoms < 1) {
      cerr << "The minimum number of atoms in non-ring fragments to write (-s) must be a "
              "whole positive number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will write non ring fragments with " << write_non_ring_atoms
           << " or more atoms\n";
    }
  }

  if (cl.option_present('f')) {
    lookup_in_database = 1;

    IWString fname;
    const_IWSubstring f;
    for (int i = 0; cl.value('f', f, i); ++i) {
      if (f.starts_with("file=") || f.starts_with("FILE=")) {
        f.remove_leading_chars(5);
        fname = f;
      } else if (!cl.value('f', min_urs_needed) || min_urs_needed < 0) {
        cerr << "The filter option requires a whole non-negative cutoff (-f)\n";
        usage(2);
      }
    }

    if (verbose) {
      cerr << "Will discard molecules with URS values below " << min_urs_needed << "\n";
    }

    if (!fname.empty()) {
      if (!fname.ends_with(".smi")) {
        fname << ".smi";
      }
      if (!stream_for_rejected_molecules.open(fname)) {
        cerr << "Cannot open stream for discarded molecules '" << fname << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Discarded molecules written to '" << fname << "'\n";
      }
    }
  }

  if (cl.option_present('c')) {
    remove_all_chiral_centres = 1;
    if (verbose) {
      cerr << "Chirality will be discarded\n";
    }
  } else if (cl.option_present('C')) {
    remove_invalid_chiral_centres = 1;
    if (verbose) {
      cerr << "Will remove any invalid chiral centres in rings produced\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "smi3d") {
        write_3d_smiles_if_needed = 1;
        if (verbose) {
          cerr << "Will write 3D smiles if coordinates present\n";
        }
      } else if (y == "proto") {
        full_proto_output = 1;
        if (verbose) {
          cerr << "Will write Smi2Rings::Results textproto output\n";
        }
      } else if (y == "smiproto") {
        smiles_and_textproto = 1;
        if (verbose) {
          cerr << "Will write smiles + Smi2Rings::Results textproto\n";
        }
      } else if (y.starts_with("writemissing=")) {
        y.remove_leading_chars(13);
        IWString fname(y);
        if (! fname.ends_with(".smi")) {
          fname << ".smi";
        }
        if (! stream_for_missing_ring_systems.open(fname)) {
          cerr << "Cannot open stream for missing rings '" << fname << "'\n";
          return 1;
        }
        if (verbose) {
          cerr << "Missing rings written to '" << fname << "'\n";
        }
      } else if (y == "hash") {
        hash_store = 1;
        if (verbose) {
          cerr << "Will cache rings and flush at the end\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  std::unique_ptr<absl::flat_hash_map<IWString, Smi2Rings::Ring>> dbhash;
  if (hash_store) {
    dbhash = std::make_unique<absl::flat_hash_map<IWString, Smi2Rings::Ring>>();
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!smi2rings(cl[i], input_type, dbhash, output)) {
      rc = i + 1;
      break;
    }
  }

  if (hash_store) {
    if (verbose) {
      cerr << "Writing " << dbhash->size() << " ring systems to the database\n";
    }
    WriteHash(*dbhash, databases[0]);
  }

  output.flush();

  for (auto* database : databases) {
    database->close(0);
  }

  if (verbose) {
    cerr << "Produced rings for " << molecules_read << " molecules\n";

    for (int i = 0; i < ring_system_size_array.number_elements(); i++) {
      if (ring_system_size_array[i]) {
        cerr << ring_system_size_array[i] << " ring systems had " << i << " rings\n";
      }
    }

    cerr << "Ring systems had between " << ring_system_ring_stats.minval() << " and "
         << ring_system_ring_stats.maxval() << " rings";
    if (ring_system_ring_stats.n() > 1) {
      cerr << " average " << ring_system_ring_stats.average();
    }
    cerr << '\n';

    cerr << "Ring systems had between " << ring_system_atom_stats.minval() << " and "
         << ring_system_atom_stats.maxval() << " atoms";
    if (ring_system_atom_stats.n() > 1) {
      cerr << " average " << ring_system_atom_stats.average();
    }
    cerr << '\n';

    if (store_in_database) {
      cerr << "Stored " << new_ring_systems_stored << " new ring systems\n";
    } else {
      cerr << "Found " << molecules_with_unique_ring_systems_found
           << " molecules with unique ring systems, "
           << ring_systems_not_found_in_database << " ring systems not in database\n";
      cerr << "Exemplar count values btw " << acc_exemplar.minval() << " and "
           << acc_exemplar.maxval() << " average "
           << static_cast<float>(acc_exemplar.average()) << "\n";
      for (int i = 0; i < exemplar_count.number_elements(); ++i) {
        if (exemplar_count[i]) {
          cerr << exemplar_count[i] << " rings had exemplar count " << i << '\n';
        }
      }
    }

    if (min_urs_needed > 0) {
      cerr << discarded_for_urs_too_low << " molecules discarded for URS below "
           << min_urs_needed << '\n';
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = smi2rings(argc, argv);

  return rc;
}
