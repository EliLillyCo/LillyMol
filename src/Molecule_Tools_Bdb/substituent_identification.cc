/*
  Scans a set of molecules and identifies all the substituents that can appear
  off a given atom type. Those are stored in a Berkeley DB database

  Fundamental flaw with including unsaturation in the atom typing.

  When we build the database, the atoms at the end of a double bond that is to
  be broken will be typed as unsaturated. But, if we want to have new molecules
  be able to accept a double bond, they will not match because their anchor
  point is NOT unsaturated - we would like to make it unsaturated.
  For now, this is just a known limitation
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "db_cxx.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char *prog_name = nullptr;

#define PROCESSING_FINISHED 1
#define READY_TO_PROCESS 2
#define NEXT_TIME 3

class Molecule_Specific_Temporary_Arrays {
public:
  const int _matoms;
  int *_atype;
  int *_aromatic_bond;
  int *_processing_status;
  int *_to_remove;
  IWString _molecule_name;

public:
  Molecule_Specific_Temporary_Arrays(int matoms);
  ~Molecule_Specific_Temporary_Arrays();

  void initialise_aromatic_bond_array(Molecule &m);

  void set_processing_status(int ndx, int s) { _processing_status[ndx] = s; }
  void set_processing_status(int s) {
    set_vector(_processing_status, _matoms, s);
  }

  int *to_remove();

  int processing_status(int s) const { return _processing_status[s]; }

  int atom_type(int a) const { return _atype[a]; }

  void temporarily_saturate(atom_number_t a1, atom_number_t a2);
  void undo_temporary_saturation(atom_number_t a1, atom_number_t a2);

  int bond_constant(const Bond &b, const atom_number_t,
                    const atom_number_t) const;

  const IWString &molecule_name() const { return _molecule_name; }
};

Molecule_Specific_Temporary_Arrays::Molecule_Specific_Temporary_Arrays(
    int matoms)
    : _matoms(matoms) {
  _atype = new int[matoms];
  _aromatic_bond = new_int(matoms * matoms);
  _processing_status = new int[matoms];
  _to_remove = new int[matoms];

  return;
}

Molecule_Specific_Temporary_Arrays::~Molecule_Specific_Temporary_Arrays() {
  delete[] _atype;
  delete[] _aromatic_bond;
  delete[] _processing_status;
  delete[] _to_remove;

  return;
}

int Molecule_Specific_Temporary_Arrays::bond_constant(
    const Bond &b, const atom_number_t a1, const atom_number_t a2) const {
  if (_aromatic_bond[a1 * _matoms + a2])
    return 11;

  if (b.is_single_bond())
    return 3;

  if (b.is_double_bond())
    return 5;

  return 7;
}

void Molecule_Specific_Temporary_Arrays::initialise_aromatic_bond_array(
    Molecule &m) {
  m.compute_aromaticity_if_needed();

  const auto ne = m.nedges();

  for (auto i = 0; i < ne; ++i) {
    const Bond *b = m.bondi(i);

    if (b->is_aromatic()) {
      const auto a1 = b->a1();
      const auto a2 = b->a2();
      _aromatic_bond[a1 * _matoms + a2] = 1;
      _aromatic_bond[a2 * _matoms + a1] = 1;
    }
  }

  _molecule_name = m.name();

  return;
}

int *Molecule_Specific_Temporary_Arrays::to_remove() {
  std::fill_n(_to_remove, _matoms, 0);

  return _to_remove;
}

/*
  This is a big hack and quite unstable. By looking into the atom typing code we
  observe that the numeric difference between a saturated atom and an
  unsaturated atom is 1. Very bad for stability, but I cannot think of a better
  way of doing this. In case of problems with multiple bonds, look here first
*/

void Molecule_Specific_Temporary_Arrays::temporarily_saturate(
    atom_number_t a1, atom_number_t a2) {
  _atype[a1]++;
  _atype[a2]++;

  return;
}

void Molecule_Specific_Temporary_Arrays::undo_temporary_saturation(
    atom_number_t a1, atom_number_t a2) {
  _atype[a1]--;
  _atype[a2]--;
}
/*
  We generate a bunch of different circular fingerprint values in the anchor
  molecule. For each value, that is, each anchor environment, we need a hash map
  of the substituents that are associated with that attachment environment

  Jul 2022. Needed to change the name of the class to avoid colliding with the
  one in substructure. For the ASubstituent class, the unique_smiles will be the
  key in the unordered_map
*/

static char _inter_molecule_name_separator = '.';

class ASubstituent {
private:
  IWString _first_molecule;
  int _number_instances;

public:
  ASubstituent(const IWString &s) : _first_molecule(s) {
    _number_instances = 1;
  }

  void extra() { _number_instances++; }
  void extra(const IWString &);

  const IWString &first_molecule() const { return _first_molecule; }
  void set_first_molecule(const IWString &s) { _first_molecule = s; }
  int number_instances() const { return _number_instances; }

  void write_first_and_instances(IWString &) const;
};

/*
this causes all kinds of problems generating ambiguous name re4solutions
template <typename T>
T &
operator << (T & os, const ASubstituent & s)
{
  os << s.first_molecule() << ':' << s.number_instances();

  return os;
}*/

void ASubstituent::write_first_and_instances(IWString &s) const {
  s << _first_molecule << ':' << _number_instances;

  // cerr << "First '" << _first_molecule << "' " << _number_instances << endl;

  return;
}

void ASubstituent::extra(const IWString &mname) {
  _first_molecule << _inter_molecule_name_separator << mname;
  _number_instances++;

  return;
}

class Substituent_Hash {
private:
public:
#if defined(IW_INTEL_COMPILER)

  static const size_t bucket_size = 4;
  static const size_t min_buckets = 8;
  bool operator()(const ASubstituent &, const IWString &) const;

#endif

  size_t operator()(const ASubstituent &) const;
};

size_t Substituent_Hash::operator()(const ASubstituent &s) const { return 1; }

typedef std::unordered_map<IWString, ASubstituent, IWStringHash>
    Usmi2Substituent;

/*
  When doing MMP studies, there are a number of things we need to keep track of
*/

class MMP_Related {
private:
public:
  MMP_Related();
};

MMP_Related::MMP_Related() { return; }

class SubstituentIdentification {
private:
  int _verbose;
  int _molecules_read;
  Chemical_Standardisation _chemical_standardisation;
  int _reduce_to_largest_fragment;
  int _remove_chirality;

  //  When we chop off a substituent, what is the smallest fragment size
  //  that can be left. For example, we would not want to consider benzene as
  //  a substituent of oxygen in phenol

  int _min_residual_size;

  int _max_substituent_size;

  int _process_hydrogen_substituents;

  //  When we are building a molecule, we can impose a limit on the number of
  //  atoms that are stripped off. Otherwise we get non useful results

  int _max_atoms_lost_during_build;

  Db _db;

  int _shell_radius;

  //  Do next shell atoms keep track of where they are attached

  int _precise_fingerprints;

  int _make_implicit_hydrogens_explicit;

  Atom_Typing_Specification _atom_typing_specification;

  std::unordered_map<unsigned int, Usmi2Substituent> *_bit;

  Report_Progress _report_progress;

  //  Some characteristics of the database being created

  unsigned int _pairs_stored;
  unsigned int _keys_stored;

  //  If we are now doing replacements, we need a means of identifying the
  //  anchor points

  resizable_array_p<Substructure_Query> _anchor_query;

  //  Or we can just take a default, which is all atoms with Hydrogens, and all
  //  non-ring single bonds

  int _default_new_molecule_starting_points;

  //  We can interpret the anchor query as defining a bond to be broken.

  int _break_molecule_at_first_two_matched_atoms;

  //  And during lookups, we may want to impose a minimum shell radius

  int _min_shell_radius;

  int _only_produce_molecules_at_biggest_radius;

  int _min_substituent_size;

  unsigned int _only_add_bond;

  resizable_array_p<Substructure_Query> _substituents_must_contain;

  int _invalid_valences_ignored;

  int _max_atoms_in_product;

  int _molecules_hitting_queries;
  int _molecules_written;

  int _write_parent_molecule;

  int _remove_isotopes_from_product;

  int _apply_atom_map_labels;

  int _max_molecules_per_input_molecule;
  int _molecules_producing_too_many_new_molecules;

  Accumulator_Int<int> _acc_new_molecules_produced;

  int _min_examples_needed_for_addition;

  int _write_fragments_added;

  int _concatenate_all_examples;

  int _strip_leading_zeros;

  int _matched_atoms_to_process;

  Db **_dbs;
  int _ndb;

  Numeric_Data_From_File<float> _expt;

  // private functions

  void _default_values();
  void _usage(int rc);
  void _preprocess(Molecule &m);

  int _ok_atom_count(const int matoms) const;
  int _matches_one_of_substituents_must_contain(Molecule &m);

  int _store_radius();
  int _get_radius();

  int _get_matched_atoms_to_process(const Set_of_Atoms *);

  int _read_smiles_already_found(const char *fname,
                                 IW_STL_Hash_Set &smiles_already_found);
  int _read_smiles_already_found(data_source_and_type<Molecule> &input,
                                 IW_STL_Hash_Set &smiles_already_found);

  //  Functions used during replacement

  int _check_already_made(const Molecule &m,
                          IW_STL_Hash_Set &alread_processed) const;

  int _enough_examples(const const_IWSubstring &fromdb) const;

  int _look_for_new_substituents(Molecule &m,
                                 IW_STL_Hash_Set &smiles_already_found,
                                 IWString_and_File_Descriptor &output);
  int _look_for_new_substituents(Molecule &m, atom_number_t zatom,
                                 const IWString &smiles_of_fragment_lost,
                                 IW_STL_Hash_Set &smiles_already_found,
                                 Molecule_Specific_Temporary_Arrays &msta,
                                 IWString_and_File_Descriptor &output);
  int _look_for_new_substituents(Molecule &m, atom_number_t zatom,
                                 const IWString &smiles_of_fragment_lost,
                                 IW_STL_Hash_Set &smiles_already_found,
                                 struct DBKey &dbkey,
                                 int &new_molecules_produced,
                                 IWString_and_File_Descriptor &output);
  int _look_for_new_substituents_db(Molecule &m, atom_number_t zatom,
                                    const IWString &smiles_of_fragment_lost,
                                    IW_STL_Hash_Set &smiles_already_found,
                                    const struct DBKey &dbkey, Db &db,
                                    Dbt &zkey,
                                    IW_STL_Hash_Set &alread_processed,
                                    IWString_and_File_Descriptor &output);
  int _form_new_molecules(Molecule &m, atom_number_t zatom,
                          const IWString &smiles_of_fragment_lost,
                          IW_STL_Hash_Set &smiles_already_found,
                          const struct DBKey &rad_and_bit,
                          const_IWSubstring &fromdb,
                          IW_STL_Hash_Set &alread_processed,
                          IWString_and_File_Descriptor &output);
  int _form_new_molecule(Molecule &m, atom_number_t zatom,
                         const IWString &smiles_of_fragment_lost,
                         IW_STL_Hash_Set &smiles_already_found,
                         const struct DBKey &rad_and_bit,
                         const const_IWSubstring &buffer,
                         IW_STL_Hash_Set &alread_processed,
                         IWString_and_File_Descriptor &output);

  int _matched_pairs_qsar(Molecule &m, const atom_number_t zatom,
                          const IWString &smiles_of_fragment_lost,
                          IW_STL_Hash_Set &smiles_already_found,
                          Molecule_Specific_Temporary_Arrays &msta,
                          MMP_Related &mmp);
  int _matched_pairs_qsar_across_bond(Molecule &m, const atom_number_t a1,
                                      const atom_number_t a2,
                                      IW_STL_Hash_Set &smiles_already_found,
                                      Molecule_Specific_Temporary_Arrays &msta,
                                      MMP_Related &mmp);

  int _matched_pairs_qsar(Molecule &m, IWString_and_File_Descriptor &output);
  int _matched_pairs_qsar(Molecule &m, atom_number_t zatom,
                          const IWString &smiles_of_fragment_lost,
                          IW_STL_Hash_Set &smiles_already_found,
                          const struct DBKey &rad_and_bit,
                          const_IWSubstring &fromdb,
                          IW_STL_Hash_Set &alread_processed,
                          IWString_and_File_Descriptor &output);
  int _matched_pairs_qsar_db(Molecule &m, const atom_number_t zatom,
                             const IWString &smiles_of_fragment_lost,
                             IW_STL_Hash_Set &smiles_already_found, Dbt &dbkey,
                             MMP_Related &mmp);

  atom_number_t _do_break_across_first_two_matched_atoms(
      Molecule &m, const atom_number_t a1, const atom_number_t a2,
      Molecule_Specific_Temporary_Arrays &msta,
      IWString &smiles_of_fragment_lost);
  int _identify_fragment_being_removed(Molecule &m, const atom_number_t a0,
                                       const int *to_remove,
                                       IWString &smiles_of_fragment_lost) const;
  int _identify_non_ring_single_bonds(Molecule &m,
                                      int *process_these_pairs) const;
  int _try_bond_breaking_across(const Molecule &m, const atom_number_t a0,
                                const atom_number_t a1, const int *to_remove,
                                IW_STL_Hash_Set &smiles_already_found,
                                Molecule_Specific_Temporary_Arrays &msta,
                                IWString_and_File_Descriptor &output);
  int _try_bond_breaking_across(const Molecule &m, const atom_number_t a0,
                                const atom_number_t a1,
                                IW_STL_Hash_Set &smiles_already_found,
                                Molecule_Specific_Temporary_Arrays &msta,
                                IWString_and_File_Descriptor &output);

  //  Functions used during database builds

  void _do_build_database_report(std::ostream &os) const;
  void _do_create_molecules_report(std::ostream &os) const;

  void _associate_substituent_with_bit(
      int r, unsigned int b, bond_type_t bt, Molecule &substituent,
      const Molecule_Specific_Temporary_Arrays &msta);
  int _compute_environment2x(Molecule &m, int radius, int rmax,
                             unsigned int *rc,
                             Molecule_Specific_Temporary_Arrays &msta) const;
  int _compute_environment01(Molecule &m, const atom_number_t zatom, int rmax,
                             unsigned int *rc,
                             Molecule_Specific_Temporary_Arrays &msta) const;

  void _initialise_msta(Molecule &m, int *number,
                        Molecule_Specific_Temporary_Arrays &msta);

  int _id_attch_pt_and_make_substituent_associations(
      Molecule &m1, bond_type_t, Molecule &m2,
      Molecule_Specific_Temporary_Arrays &msta);
  int _divide_molecule(Molecule &m, const atom_number_t a1,
                       const atom_number_t a2,
                       Molecule_Specific_Temporary_Arrays &msta);
  int _build_database(Molecule &m, IWString_and_File_Descriptor &output);
  int _process_molecules(data_source_and_type<Molecule> &input,
                         IW_STL_Hash_Set &smiles_already_found,
                         IWString_and_File_Descriptor &output);
  int _process_molecules(const char *fname, FileType input_type,
                         IW_STL_Hash_Set &smiles_already_found,
                         IWString_and_File_Descriptor &output);

  int _write_in_memory_hashes_to_database();
  int _write_in_memory_hash_to_database(
      int radius, const std::unordered_map<unsigned int, Usmi2Substituent> &b);
  int _write_substituents_associated_with_bit(const struct DBKey &dbkey,
                                              const Usmi2Substituent &u);

public:
  SubstituentIdentification();
  ~SubstituentIdentification();

  int operator()(int argc, char **argv);
};

SubstituentIdentification::SubstituentIdentification()
    : _db(NULL, DB_CXX_NO_EXCEPTIONS) {
  _default_values();

  return;
}

void SubstituentIdentification::_default_values() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 1; // always
  _remove_chirality = 0;

  _min_residual_size = 3;
  _max_substituent_size = 1;

  _process_hydrogen_substituents = 0;

  _max_atoms_lost_during_build = std::numeric_limits<int>::max();

  _precise_fingerprints = 1;

  _shell_radius = 0;

  _make_implicit_hydrogens_explicit = 0;

  _bit = nullptr;

  _min_shell_radius = 0;

  _only_produce_molecules_at_biggest_radius = 0;

  _min_substituent_size = 0;

  _only_add_bond = SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND;

  _invalid_valences_ignored = 0;

  _max_atoms_in_product = std::numeric_limits<int>::max();

  _remove_isotopes_from_product = 0;

  _apply_atom_map_labels = 0;

  _pairs_stored = 0;
  _keys_stored = 0;

  _default_new_molecule_starting_points = 0;

  _molecules_hitting_queries = 0;
  _molecules_written = 0;

  _write_parent_molecule = 0;

  _max_molecules_per_input_molecule = 10000;

  _molecules_producing_too_many_new_molecules = 0;

  _min_examples_needed_for_addition = 1;

  _write_fragments_added = 0;

  _break_molecule_at_first_two_matched_atoms = 0;

  _concatenate_all_examples = 0;

  _strip_leading_zeros = 0;

  _matched_atoms_to_process = 1;

  _dbs = nullptr;
  _ndb = 0;

  return;
}

SubstituentIdentification::~SubstituentIdentification() {
  _db.close(0);

  if (nullptr != _bit)
    delete[] _bit;

  if (nullptr != _dbs) {
    for (auto i = 0; i < _ndb; ++i) {
      _dbs[i]->close(0);
    }
    delete[] _dbs;
  }

  return;
}

void
SubstituentIdentification::_usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Identifies substituents in a set of molecules. Then looks for replacements\n";
  cerr << "  -d <dbname>   Berkeley database for substituent data\n";
  cerr << "  -R <bonds>    radius for circular fingerprint from attachment point\n";
  cerr << "  -m <natoms>   min number of atoms in a substituent\n";
  cerr << "  -M <natoms>   max number of atoms in a substituent\n";
  cerr << "  -B            flag to enable building\n";
  cerr << "  -w <natoms>   during building, min number of atoms left after removing substituent\n";
  cerr << "  -e            during building, concatenate names of all example molecules\n";
  cerr << "  -h            during building, discern Hydrogen as a substituent\n";
  cerr << "  -q <query>    during addition, query  describing anchor atoms needing new substituents\n";
  cerr << "  -s <smarts>   during addition, smarts describing anchor atoms needing new substituents\n";
  cerr << "  -y            during addition, process all atoms with H's and all non ring single bonds (no -q or -s)\n";
  cerr << "  -k            during addition, break bond btw first two matched atoms, discard first frag, and join\n";
  cerr << "                replacement requires -k option e.g. to replace an aniline, try -s '-k -s c-[NH2]'\n";
  cerr << "  -L <natoms>   during addition, max atoms lost from parent with -k option\n";
  cerr << "  -r <bonds>    during addition, min radius to be considered\n";
  cerr << "  -a            during addition, only produce molecules at the largest radius found\n";
  cerr << "  -b S,D,T      during addition, only add via the kinds of bond(s) specified\n";
  cerr << "  -C <natoms>   during addition, discard any products with more than <natoms> atoms\n";
  cerr << "  -H ...        during addition, queries that substituents must contain\n";
  cerr << "  -p            during addition, write starting molecule\n";
  cerr << "  -I            during addition, remove isotopes before writing\n";
  cerr << "  -x <number>   during addition, max number of molecules produced per input molecule (def 10000)\n";
  cerr << "  -u <number>   during addition, min number of example molecules needed for a fragment to be added\n";
  cerr << "  -j            during addition, write the fragment added\n";
  cerr << "  -n <number>   during addition, process <number> matched atoms in each query (default 1)\n";
  cerr << "  -V <fname>    during addition, avoid forming any of the molecules in <fname>\n";
  cerr << "  -t <number>   report progress every <number> molecules processed\n";
  cerr << "  -z            strip leading zero's from identifiers before storing\n";
  cerr << "  -l            reduce to largest fragment (automatically turned on)\n";
  cerr << "  -c            remove all chirality (automatically turned on during building)\n";
  cerr << "  -U            label via atom map numbers rather than isotopes\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

void SubstituentIdentification::_preprocess(Molecule &m) {
  if (_reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  m.revert_all_directional_bonds_to_non_directional();

  if (_remove_chirality)
    m.remove_all_chiral_centres();

  if (_chemical_standardisation.active())
    _chemical_standardisation.process(m);

  if (_make_implicit_hydrogens_explicit)
    m.make_implicit_hydrogens_explicit();

  if (_strip_leading_zeros) {
    const_IWSubstring mname = m.name();
    if (mname.starts_with('0')) {
      mname.remove_leading_chars('0');
      m.set_name(mname);
    }
  }

  return;
}

int SubstituentIdentification::_ok_atom_count(const int matoms) const {
  // cerr << "_ok_atom_count checking " << matoms << " vs " <<
  // _min_substituent_size << " to " << _max_substituent_size << endl;

  if (matoms > _max_substituent_size)
    return 0;

  if (matoms < _min_substituent_size)
    return 0;

  return 1;
}

struct DBKey {
  int radius;
  unsigned int b;
};

int SubstituentIdentification::_write_in_memory_hashes_to_database() {
  for (auto i = 0; i <= _shell_radius; ++i) {
    if (!_write_in_memory_hash_to_database(i, _bit[i]))
      return 0;
  }

  return 1;
}

int SubstituentIdentification::_write_in_memory_hash_to_database(
    int radius, const std::unordered_map<unsigned int, Usmi2Substituent> &b) {
  struct DBKey dbkey {
    radius, 0
  };

  for (auto i = b.cbegin(); i != b.cend(); ++i) {
    dbkey.b = (*i).first;

    if (!_write_substituents_associated_with_bit(dbkey, (*i).second))
      return 0;
  }

  return 1;
}

int SubstituentIdentification::_write_substituents_associated_with_bit(
    const struct DBKey &dbkey, const Usmi2Substituent &u) {
  IWString to_store(32768);

  for (auto i = u.cbegin(); i != u.cend(); ++i) {
    const IWString &usmi = (*i).first;
    const ASubstituent &s = (*i).second;

    _pairs_stored += s.number_instances();

    if (i != u.cbegin())
      to_store << '|';

    to_store << usmi << ':';
    s.write_first_and_instances(to_store);
  }

  Dbt zkey{(void *)(&dbkey), sizeof(dbkey)};
  Dbt zdata{to_store.rawdata(), static_cast<u_int32_t>(to_store.length())};

  // cerr << "Writing radius " << dbkey.radius << " " << dbkey.b << endl;

  int s = _db.put(NULL, &zkey, &zdata, 0);

  if (0 != s) {
    cerr << "Berkeley database Put operation failed\n";
    _db.err(s, "");
    return 0;
  }

  _keys_stored++;

  return 1;
}

void SubstituentIdentification::_initialise_msta(
    Molecule &m, int *number, Molecule_Specific_Temporary_Arrays &msta) {
  _atom_typing_specification.assign_atom_types(m, msta._atype);

  msta.initialise_aromatic_bond_array(m);

  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    number[i] = i;
    m.set_user_specified_atom_void_ptr(i, number + i);
  }

  return;
}

/*
  We need to remove the fragment containing a0, but we need to capture
  the unique smiles of the fragment being removed
*/

int SubstituentIdentification::_identify_fragment_being_removed(
    Molecule &m, const atom_number_t attachment_point, const int *to_remove,
    IWString &smiles_of_fragment_lost) const {
  if (_apply_atom_map_labels)
    m.set_atom_map_number(attachment_point, 1);
  else
    m.set_isotope(attachment_point, 1);

  Molecule being_discarded;
  m.create_subset(being_discarded, to_remove, 1);

  smiles_of_fragment_lost = being_discarded.unique_smiles();

  return 1;
}

/*
  This can be made a lot more efficient
*/

atom_number_t
SubstituentIdentification::_do_break_across_first_two_matched_atoms(
    Molecule &m, const atom_number_t a0, const atom_number_t a1,
    Molecule_Specific_Temporary_Arrays &msta,
    IWString &smiles_of_fragment_lost) {
  if (!m.are_bonded(a0, a1)) {
    cerr << "SubstituentIdentification::_do_break_across_first_two_matched_"
            "atoms:atoms not bonded\n";
    return INVALID_ATOM_NUMBER;
  }

#ifdef DEBUG_DO_BREAK_ACROSS_FIRST_TWO_MATCHED_ATOMS
  Molecule mcopy(m);
  for (int i = 0; i < m.natoms(); ++i) {
    mcopy.set_atom_map_number(i, i);
  }
  cerr << mcopy.smiles() << " atoms " << a0 << " and " << a1 << endl;
#endif

  int *to_remove = msta.to_remove();

  const auto atoms_being_removed =
      m.identify_side_of_bond(to_remove, a1, 1, a0);

  if (1 == m.ncon(a1)) //  a common case
  {
    const int delta = a1 < a0 ? 1 : 0;
    m.remove_atom(a1);
    return a0 - delta;
  }

  // cerr << "Found " << atoms_being_removed << " atoms to be removed\n";

  if (0 == atoms_being_removed) // should not happen
    return INVALID_ATOM_NUMBER;

  if (atoms_being_removed > _max_atoms_lost_during_build)
    return INVALID_ATOM_NUMBER;

  // We are going to remove all the atoms connected to a1, but we then need to
  // work out the new atom number of atom a0 which will be retained in the
  // resulting molecule.

  const int delta = std::count_if(to_remove, to_remove + a0,
                                  [](const int x) { return 1 == x; });

  if (_write_fragments_added) // we need to capture the unique smiles of what is
                              // discarded
    _identify_fragment_being_removed(m, a1, to_remove, smiles_of_fragment_lost);

  m.remove_atoms(to_remove);

  return a0 - delta;
}

static int symmetrically_equivalent_to_something_already_found(
    const Bond &b, const int *process_these_pairs, const int ndx,
    const int *symmetry_classes) {
  const auto sa1 = symmetry_classes[b.a1()];
  const auto sa2 = symmetry_classes[b.a2()];

  for (int i = 0; i < ndx; i += 2) {
    const auto si1 = symmetry_classes[process_these_pairs[i]];
    const auto si2 = symmetry_classes[process_these_pairs[i + 1]];

    if (si1 == sa1 && si2 == sa2)
      return 1;
    if (si1 == sa2 && si2 == sa1)
      return 1;
  }

  return 0;
}

int SubstituentIdentification::_identify_non_ring_single_bonds(
    Molecule &m, int *process_these_pairs) const {
  (void)m.ring_membership();

  const int *symmetry_classes = m.symmetry_classes();

  int ndx = 0;

  const int nb = m.nedges();

  for (int i = 0; i < nb; ++i) {
    const Bond *b = m.bondi(i);

    if (!b->is_single_bond())
      continue;

    if (b->nrings())
      continue;

    if (symmetrically_equivalent_to_something_already_found(
            *b, process_these_pairs, ndx, symmetry_classes))
      continue;

    process_these_pairs[ndx] = b->a1();
    ndx++;
    process_these_pairs[ndx] = b->a2();
    ndx++;
  }

  return ndx / 2;
}

// We are going to remove all the atoms connected to a0, so we need to work out
// the new atom number of atom a1 which will be retained in the resulting
// molecule.

int SubstituentIdentification::_try_bond_breaking_across(
    const Molecule &m, const atom_number_t a0, const atom_number_t a1,
    const int *to_remove, IW_STL_Hash_Set &smiles_already_found,
    Molecule_Specific_Temporary_Arrays &msta,
    IWString_and_File_Descriptor &output) {
  IWString smiles_of_fragment_lost;

  Molecule mcopy(m);

  const int delta = std::count_if(to_remove, to_remove + a1,
                                  [](const int x) { return 1 == x; });

  if (_write_fragments_added) // we need to capture the unique smiles of what is
                              // discarded
    _identify_fragment_being_removed(mcopy, a0, to_remove,
                                     smiles_of_fragment_lost);

  mcopy.remove_atoms(to_remove);

  return _look_for_new_substituents(mcopy, a1 - delta, smiles_of_fragment_lost,
                                    smiles_already_found, msta, output);
}

/*
  We try both ways across the bond
*/

int SubstituentIdentification::_try_bond_breaking_across(
    const Molecule &m, const atom_number_t a0, const atom_number_t a1,
    IW_STL_Hash_Set &smiles_already_found,
    Molecule_Specific_Temporary_Arrays &msta,
    IWString_and_File_Descriptor &output) {
  if (!m.are_bonded(a0, a1)) {
    cerr << "SubstituentIdentification::_try_bond_breaking_across:atoms not "
            "bonded '"
         << m.name() << "' atoms " << a0 << " " << a1 << endl;
    return 0;
  }

  const int matoms = m.natoms();

  int *to_remove = msta.to_remove(); // comes out zero filled

  const auto atoms_being_removed =
      m.identify_side_of_bond(to_remove, a0, 1, a1);

  if (0 == atoms_being_removed) // should not happen
    return 0;

  if (atoms_being_removed <= _max_atoms_lost_during_build)
    _try_bond_breaking_across(
        m, a0, a1, to_remove, smiles_already_found, msta,
        output); // kind of dangerous to pass both to_remove and msta...

  // Now try things the other way around

  if (matoms - atoms_being_removed <= _max_atoms_lost_during_build) {
    std::transform(to_remove, to_remove + matoms, to_remove,
                   [](const int x) { return !x; });

    _try_bond_breaking_across(m, a1, a0, to_remove, smiles_already_found, msta,
                              output);
  }

  return 1;
}

/*
  When doing matched pairs, we have a lot of situations.
  There may be multiple sites of matched pairs.
  There may be different numbers of examples at each of those sites.
  The shifts will be different within a site
*/

int SubstituentIdentification::_matched_pairs_qsar(
    Molecule &m, IWString_and_File_Descriptor &output) {
  const auto matoms = m.natoms();

  // cerr << "Looking for new substituents to " << m.smiles() << endl;

  int *numbers = new int[matoms];
  std::unique_ptr<int[]> free_numbers(numbers);

  Molecule_Specific_Temporary_Arrays msta(matoms);

  _initialise_msta(m, numbers, msta);

  int queries_hitting_this_molecule = 0;

  IW_STL_Hash_Set smiles_already_found;

  smiles_already_found.insert(
      m.unique_smiles()); // do not want to recreate the parent

  MMP_Related mmp;

  if (_default_new_molecule_starting_points) // non ring single bonds and atoms
                                             // with Hydrogen
  {
    IWString smiles_of_fragment_lost("H");

    for (int i = 0; i < matoms; ++i) {
      if (0 == m.hcount(i))
        continue;

      _matched_pairs_qsar(m, i, smiles_of_fragment_lost, smiles_already_found,
                          msta, mmp);
      queries_hitting_this_molecule++;
    }

    int *process_these_pairs = new int[matoms * matoms];
    std::unique_ptr<int[]> free_process_these_pairs(process_these_pairs);
    const auto nb = _identify_non_ring_single_bonds(m, process_these_pairs);
    if (nb > 0) {
      queries_hitting_this_molecule++;

      for (int i = 0; i < nb; ++i) {
        atom_number_t a1 = process_these_pairs[i + i];
        atom_number_t a2 = process_these_pairs[i + i + 1];

        _matched_pairs_qsar_across_bond(m, a1, a2, smiles_already_found, msta,
                                        mmp);
      }
    }
  } else // based on queries
  {
    Molecule_to_Match target(&m);

    const int nq = _anchor_query.number_elements();

    for (auto i = 0; i < nq; ++i) {
      Substructure_Results sresults;

      int nhits = _anchor_query[i]->substructure_search(target, sresults);

      if (0 == nhits)
        continue;

      queries_hitting_this_molecule++;

      if (_break_molecule_at_first_two_matched_atoms) {
        for (auto j = 0; j < nhits; ++j) {
          Molecule mcopy(m);
          IWString smiles_of_fragment_lost;
          const Set_of_Atoms *e = sresults.embedding(j);

          if (e->number_elements() < 2) // should never happen
            continue;

          _matched_pairs_qsar_across_bond(mcopy, e->item(0), e->item(1),
                                          smiles_already_found, msta, mmp);
        }
      } else {
        IWString smiles_of_fragment_lost("H");
        for (auto j = 0; j < nhits; ++j) {
          const atom_number_t k = sresults.embedding(j)->item(0);

          if (0 == m.hcount(k)) // will not happen if the person constructed the
                                // query properly
            continue;

          if (!_matched_pairs_qsar(m, k, smiles_of_fragment_lost,
                                   smiles_already_found, msta, mmp))
            continue;
        }
      }
    }
  }

  if (queries_hitting_this_molecule)
    _molecules_hitting_queries++;

  return 1;
}

int SubstituentIdentification::_get_matched_atoms_to_process(
    const Set_of_Atoms *e) {
  const int n = e->number_elements();

  if (n > _matched_atoms_to_process)
    return _matched_atoms_to_process;

  return n;
}

int SubstituentIdentification::_look_for_new_substituents(
    Molecule &m, IW_STL_Hash_Set &smiles_already_found,
    IWString_and_File_Descriptor &output) {
  const auto matoms = m.natoms();

  if (matoms >= _max_atoms_in_product) {
    cerr << m.name() << " already has " << matoms
         << " atoms, but max in products is " << _max_atoms_in_product
         << " skipped\n";
    return 1;
  }

  // cerr << "Looking for new substituents to " << m.smiles() << endl;

  int *numbers = new int[matoms];
  std::unique_ptr<int[]> free_numbers(numbers);

  Molecule_Specific_Temporary_Arrays msta(matoms);

  _initialise_msta(m, numbers, msta);

  int queries_hitting_this_molecule = 0;

  int parent_written = 0;

  smiles_already_found.insert(
      m.unique_smiles()); // do not want to recreate the parent

  int molecules_created = 0;

  if (_default_new_molecule_starting_points) // non ring single bonds and atoms
                                             // with Hydrogen atoms
  {
    IWString smiles_of_fragment_lost("H");

    for (int i = 0; i < matoms; ++i) {
      if (0 == m.hcount(i))
        continue;

      molecules_created += _look_for_new_substituents(
          m, i, smiles_of_fragment_lost, smiles_already_found, msta, output);
      queries_hitting_this_molecule++;
      if (molecules_created > _max_molecules_per_input_molecule)
        break;
    }

    int *process_these_pairs = new int[matoms * matoms];
    std::unique_ptr<int[]> free_process_these_pairs(process_these_pairs);
    const auto nb = _identify_non_ring_single_bonds(m, process_these_pairs);
    if (nb > 0) {
      queries_hitting_this_molecule++;

      for (int i = 0; i < nb; ++i) {
        atom_number_t a1 = process_these_pairs[i + i];
        atom_number_t a2 = process_these_pairs[i + i + 1];

        if (_write_parent_molecule && 0 == i)
          output << m.smiles() << ' ' << m.name() << '\n';

        molecules_created += _try_bond_breaking_across(
            m, a1, a2, smiles_already_found, msta, output);
        if (molecules_created > _max_molecules_per_input_molecule)
          break;
      }
    }
  } else {
    Molecule_to_Match target(&m);

    const int nq = _anchor_query.number_elements();

    for (auto i = 0; i < nq; ++i) {
      Substructure_Results sresults;

      int nhits = _anchor_query[i]->substructure_search(target, sresults);

      if (0 == nhits)
        continue;

      //    cerr << nhits << " hits to " << m.name() << endl;

      queries_hitting_this_molecule++;

      if (_write_parent_molecule && !parent_written) {
        output << m.smiles() << ' ' << m.name() << '\n';
        parent_written = 1;
      }

      if (_break_molecule_at_first_two_matched_atoms) {
        for (auto j = 0; j < nhits; ++j) {
          Molecule mcopy(m);
          IWString smiles_of_fragment_lost;
          const Set_of_Atoms *e = sresults.embedding(j);

          if (e->number_elements() < 2) // should never happen
            continue;

          const auto attachment_point =
              _do_break_across_first_two_matched_atoms(
                  mcopy, e->item(0), e->item(1), msta, smiles_of_fragment_lost);
          //        cerr << m.name() << ' ' << j << " atoms " << e->item(0) << "
          //        and " << e->item(1) << " attachment_point " <<
          //        attachment_point << endl;
          if (INVALID_ATOM_NUMBER == attachment_point)
            continue;

          //        cerr << "From " << m.smiles() << "\nget  " << mcopy.smiles()
          //        << " atom " << attachment_point << " type " <<
          //        mcopy.smarts_equivalent_for_atom(attachment_point) << endl;

          molecules_created += _look_for_new_substituents(
              mcopy, attachment_point, smiles_of_fragment_lost,
              smiles_already_found, msta, output);
          if (molecules_created > _max_molecules_per_input_molecule)
            break;
        }
      } else {
        IWString smiles_of_fragment_lost("H");
        for (auto j = 0; j < nhits; ++j) {
          const int nprocess =
              _get_matched_atoms_to_process(sresults.embedding(j));

          for (int x = 0; x < nprocess; ++x) {
            const atom_number_t k = sresults.embedding(j)->item(x);

            if (INVALID_ATOM_NUMBER == k)
              continue;

            if (0 == m.hcount(k)) // will not happen if the person constructed
                                  // the query properly
              continue;

            molecules_created +=
                _look_for_new_substituents(m, k, smiles_of_fragment_lost,
                                           smiles_already_found, msta, output);

            if (molecules_created > _max_molecules_per_input_molecule)
              break;
          }
        }
      }
    }
  }

  if (queries_hitting_this_molecule)
    _molecules_hitting_queries++;

  if (molecules_created > _max_molecules_per_input_molecule) {
    _molecules_producing_too_many_new_molecules++;
    if (_verbose)
      cerr << m.smiles() << ' ' << m.name() << " would have produced more than "
           << _max_molecules_per_input_molecule << " new molecules\n";
  }

  return 1;
}

int SubstituentIdentification::_look_for_new_substituents(
    Molecule &m, atom_number_t zatom, const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found,
    Molecule_Specific_Temporary_Arrays &msta,
    IWString_and_File_Descriptor &output) {
  unsigned int *b = new unsigned int[_shell_radius + 1];
  std::unique_ptr<unsigned int[]> free_b(b);

  auto max_shell_radius_formed =
      _compute_environment01(m, zatom, _shell_radius, b, msta);

// #define ECHO_BITS_FORMED
#ifdef ECHO_BITS_FORMED
  cerr << "Looking for additions to " << m.smiles() << " atom " << zatom << ' '
       << m.smarts_equivalent_for_atom(zatom) << " bits";
  for (auto i = 0; i <= max_shell_radius_formed; ++i) {
    cerr << ' ' << b[i];
  }
  cerr << endl;
#endif

  int new_molecules_produced = 0;

  for (auto r = max_shell_radius_formed; r >= _min_shell_radius; --r) {
    DBKey dbkey{r, b[r]};

    const auto tmp = _look_for_new_substituents(
        m, zatom, smiles_of_fragment_lost, smiles_already_found, dbkey,
        new_molecules_produced, output);

    //  cerr << tmp << " molecules generated at radois " << r << endl;

    if (0 == tmp)
      continue;

    new_molecules_produced += tmp;

    if (new_molecules_produced > _max_molecules_per_input_molecule)
      return new_molecules_produced;

    if (_only_produce_molecules_at_biggest_radius)
      return new_molecules_produced;
  }

  // cerr << m.name() << " generated " << new_molecules_produced << endl;

  return new_molecules_produced;
}

int SubstituentIdentification::_look_for_new_substituents(
    Molecule &m, atom_number_t zatom, const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found, struct DBKey &dbkey,
    int &new_molecules_produced, IWString_and_File_Descriptor &output) {
  Dbt zkey{(void *)(&dbkey), sizeof(dbkey)};

  // cerr << "Looking for bit " << dbkey.b << " radius " << dbkey.radius <<
  // endl;

  IW_STL_Hash_Set
      alread_processed; // names of reagents when dealing with multiple
                        // databases, don't re-create molecules

  if (0 == _ndb)
    return _look_for_new_substituents_db(m, zatom, smiles_of_fragment_lost,
                                         smiles_already_found, dbkey, _db, zkey,
                                         alread_processed, output);

  int rc = 0;
  for (auto i = 0; i < _ndb; ++i) {
    int tmp = _look_for_new_substituents_db(
        m, zatom, smiles_of_fragment_lost, smiles_already_found, dbkey,
        *(_dbs[i]), zkey, alread_processed, output);

    new_molecules_produced += tmp;
    rc += tmp;

    if (new_molecules_produced > _max_molecules_per_input_molecule)
      return rc;
  }

  return rc;
}

int SubstituentIdentification::_look_for_new_substituents_db(
    Molecule &m, atom_number_t zatom, const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found, const struct DBKey &rad_and_bit,
    Db &db, Dbt &dbkey, IW_STL_Hash_Set &alread_processed,
    IWString_and_File_Descriptor &output) {
  Dbt zdata;

  int s = db.get(NULL, &dbkey, &zdata, 0);
  if (0 != s)
    return 0;

  // cerr << "Found bit in database\n";

  const_IWSubstring ss(reinterpret_cast<const char *>(zdata.get_data()),
                       zdata.get_size());

  return _form_new_molecules(m, zatom, smiles_of_fragment_lost,
                             smiles_already_found, rad_and_bit, ss,
                             alread_processed, output);
}

/*
  A typical database record will look like
  [1OH2]:36473:129|[1OH2]=:36473:1
*/

int SubstituentIdentification::_form_new_molecules(
    Molecule &m, atom_number_t zatom, const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found, const struct DBKey &rad_and_bit,
    const_IWSubstring &fromdb, IW_STL_Hash_Set &alread_processed,
    IWString_and_File_Descriptor &output) {
  const_IWSubstring token;
  int i = 0;

  int rc = 0;

  while (fromdb.nextword(token, i, '|')) {
    _form_new_molecule(m, zatom, smiles_of_fragment_lost, smiles_already_found,
                       rad_and_bit, token, alread_processed, output);
    rc++;
  }

  return rc;
}

/*
  If all example molecules are stored in the database, a typica entry will look
  like

  |O=C(Nc1s[1cH]c([n]1)C)C:1158204.1164662:2|
*/

/*
  We are branching from a single atom.
*/

int SubstituentIdentification::_matched_pairs_qsar(
    Molecule &m, const atom_number_t zatom,
    const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found,
    Molecule_Specific_Temporary_Arrays &msta, MMP_Related &mmp) {
  unsigned int *b = new unsigned int[_shell_radius + 1];
  std::unique_ptr<unsigned int[]> free_b(b);

  auto max_shell_radius_formed =
      _compute_environment01(m, zatom, _shell_radius, b, msta);

// #define ECHO_BITS_FORMED
#ifdef ECHO_BITS_FORMED
  cerr << "Looking for additions to " << m.smiles() << " atom " << zatom << ' '
       << m.smarts_equivalent_for_atom(zatom) << " bits";
  for (auto i = 0; i <= max_shell_radius_formed; ++i) {
    cerr << ' ' << b[i];
  }
  cerr << endl;
#endif

  int new_molecules_produced = 0;

  for (auto r = max_shell_radius_formed; r >= _min_shell_radius; --r) {
    DBKey dbkey{r, b[r]};
    Dbt zkey{(void *)(&dbkey), sizeof(dbkey)};

    const auto tmp = _matched_pairs_qsar_db(m, zatom, smiles_of_fragment_lost,
                                            smiles_already_found, zkey, mmp);

    //  cerr << tmp << " molecules generated at radois " << r << endl;

    if (0 == tmp)
      continue;

    new_molecules_produced += tmp;

    if (new_molecules_produced > _max_molecules_per_input_molecule)
      return new_molecules_produced;

    if (_only_produce_molecules_at_biggest_radius)
      return new_molecules_produced;
  }

  // cerr << m.name() << " generated " << new_molecules_produced << endl;

  return new_molecules_produced;
}

int SubstituentIdentification::_matched_pairs_qsar_db(
    Molecule &m, const atom_number_t zatom,
    const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found, struct Dbt &dbkey,
    MMP_Related &mmp) {
  Dbt fromdb;

  for (int i = 0; i < _ndb; ++i) {
    int rc = _dbs[i]->get(NULL, &dbkey, &fromdb, 0);
    if (0 != rc)
      continue;

    const_IWSubstring s(reinterpret_cast<const char *>(fromdb.get_data()),
                        fromdb.get_size());
  }

  return 1;
}

int SubstituentIdentification::_matched_pairs_qsar_across_bond(
    Molecule &m, const atom_number_t a1, const atom_number_t a2,
    IW_STL_Hash_Set &smiles_already_found,
    Molecule_Specific_Temporary_Arrays &msta, MMP_Related &mmp) {
  m.remove_bond_between_atoms(a1, a2);
  resizable_array_p<Molecule> c;
  m.create_components(c);

  int adj = 0;

  const int fma1 = m.fragment_membership(a1);

  for (int i = 0; i < a1; ++i) {
    if (m.fragment_membership(i) != fma1)
      adj++;
  }

  if (m.fragment_membership(0) == m.fragment_membership(a1))
    return _matched_pairs_qsar(*c[0], a1 - adj, c[1]->unique_smiles(),
                               smiles_already_found, msta, mmp);
  else
    return _matched_pairs_qsar(*c[1], a1 - adj, c[0]->unique_smiles(),
                               smiles_already_found, msta, mmp);

  // here

  return 1;
}

int SubstituentIdentification::_matched_pairs_qsar(
    Molecule &m, atom_number_t zatom, const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found, const struct DBKey &rad_and_bit,
    const_IWSubstring &fromdb, IW_STL_Hash_Set &alread_processed,
    IWString_and_File_Descriptor &output) {
  const_IWSubstring token;
  int i = 0;

  int rc = 0;

  while (fromdb.nextword(token, i, '|')) {
    _form_new_molecule(m, zatom, smiles_of_fragment_lost, smiles_already_found,
                       rad_and_bit, token, alread_processed, output);
    rc++;
  }

  return rc;
}

int SubstituentIdentification::_matches_one_of_substituents_must_contain(
    Molecule &m) {
  Molecule_to_Match target(&m);

  const auto nq = _substituents_must_contain.number_elements();

  for (auto i = 0; i < nq; ++i) {
    if (_substituents_must_contain[i]->substructure_search(target))
      return 1;
  }

  return 0;
}

int SubstituentIdentification::_check_already_made(
    const Molecule &m, IW_STL_Hash_Set &alread_processed) const {
  Molecule mcopy(m);
  mcopy.transform_to_non_isotopic_form();

  const IWString &usmi = mcopy.unique_smiles();

  // cerr << "Checking " << usmi << endl;

  const auto f = alread_processed.find(usmi);

  if (f != alread_processed.end())
    return 1;

  alread_processed.insert(usmi);

  return 0; // never seen this before
}

int SubstituentIdentification::_form_new_molecule(
    Molecule &m, atom_number_t zatom, const IWString &smiles_of_fragment_lost,
    IW_STL_Hash_Set &smiles_already_found, const struct DBKey &rad_and_bit,
    const const_IWSubstring &buffer, IW_STL_Hash_Set &alread_processed,
    IWString_and_File_Descriptor &output) {
// #define DEBUG_FORM_NEW_MOLECULE
#ifdef DEBUG_FORM_NEW_MOLECULE
  cerr << "Building from database contents '" << buffer << "'\n";
#endif

  const_IWSubstring smiles;
  int i = 0;
  if (!buffer.nextword(smiles, i, ':')) {
    cerr << "SubstituentIdentification::_form_new_molecule:invalid database "
            "contents '"
         << buffer << "'\n";
    return 0;
  }

  // cerr << "Looking for bonding information in '" << token << "'\n";

  bond_type_t bt = SINGLE_BOND;
  if (smiles.ends_with('=')) {
    if (m.hcount(zatom) < 2)
      return 0;

    bt = DOUBLE_BOND;
    smiles.chop();
  } else if (smiles.ends_with('#')) {
    if (m.hcount(zatom) < 3)
      return 0;

    bt = TRIPLE_BOND;
    smiles.chop();
  } else if (0 == m.hcount(zatom))
    return 0;

  if (0 == (_only_add_bond & bt))
    return 0;

  if (_ndb > 0) {
    if (alread_processed.contains(smiles))
      return 0;

    alread_processed.insert(smiles);
  }

  if (smiles ==
      smiles_of_fragment_lost) // pulling off something and putting it back on
    return 0;

  Molecule f;
  if (!f.build_from_smiles(smiles)) {
    cerr << "SubstituentIdentification::_form_new_molecule:invalid smiles '"
         << buffer << "'\n";
    return 0;
  }

  if (!_ok_atom_count(f.natoms()))
    return 0;

  const auto initial_matoms = m.natoms();

  // cerr << "Fragment contains " << f.natoms() << " atoms\n";
  if (initial_matoms + f.natoms() > _max_atoms_in_product)
    return 0;

  const auto fragment_attachment_point = f.atom_with_isotope(1);
  if (INVALID_ATOM_NUMBER == fragment_attachment_point) {
    cerr << "SubstituentIdentification::_form_new_molecule:no isotopic atom in "
            "stored fragment '"
         << buffer << "'\n";
    return 0;
  }

  if (_min_examples_needed_for_addition > 1 && !_enough_examples(buffer))
    return 0;

  if (_substituents_must_contain.number_elements() &&
      !_matches_one_of_substituents_must_contain(f))
    return 0;

  m.add_molecule(&f);

  m.add_bond(zatom, initial_matoms + fragment_attachment_point, bt);

  m.unset_all_implicit_hydrogen_information(zatom);
  m.set_implicit_hydrogens_known(zatom, 0);
  m.unset_all_implicit_hydrogen_information(initial_matoms +
                                            fragment_attachment_point);
  m.set_implicit_hydrogens_known(initial_matoms + fragment_attachment_point, 0);

  if (_remove_isotopes_from_product)
    ;
  else if (_apply_atom_map_labels)
    m.set_atom_map_number(zatom, 1);
  else
    m.set_isotope(zatom, 1);

  m.recompute_implicit_hydrogens(zatom);
  m.recompute_implicit_hydrogens(initial_matoms + fragment_attachment_point);

  const_IWSubstring token;

  buffer.nextword(token, i); // will get everything
  token++;                   // get rid of leading ':'

  if (token == m.name()) // molecule being merged with itself
    return 0;

  if (!m.valence_ok()) {
    if (_verbose > 1) {
      cerr << "Warning, invalid valence '" << m.name() << " with '" << token
           << "' " << m.smiles() << "\n";
      for (auto i = 0; i < m.natoms(); ++i) {
        if (m.valence_ok(i))
          continue;

        cerr << " problem with atom " << i << ' '
             << m.smarts_equivalent_for_atom(i) << "\n";
      }
    }
    _invalid_valences_ignored++;
  } else if (_check_already_made(m, smiles_already_found))
    ;
  else {
    if (_remove_isotopes_from_product)
      m.set_isotope(initial_matoms + fragment_attachment_point, 0);

    output << m.smiles() << ' ' << m.name() << " %% ";

    output << token << ":R" << rad_and_bit.radius << " B" << rad_and_bit.b;

    if (_write_fragments_added)
      //    output << ' ' << smiles_of_fragment_lost << '%' <<
      //    rad_and_bit.radius << '%' << smiles;
      output << ' ' << smiles_of_fragment_lost << '%' << smiles;

    output << '\n';

    _molecules_written++;

    output.write_if_buffer_holds_more_than(32768);
  }

  m.resize(initial_matoms);
  m.set_isotope(zatom, 0);

  return 1;
}

/*
  db contents looks like smiles:id:examples

  We need to see if 'examples' is >= _min_examples_needed_for_addition
*/

int SubstituentIdentification::_enough_examples(
    const const_IWSubstring &fromdb) const {
  const auto i = fromdb.rindex(':');

  if (i + 2 == fromdb.length()) // most common case, single digit
    return (fromdb.last_item() - '0') >= _min_examples_needed_for_addition;

  if (i + 3 == fromdb.length()) // two digit number of examples
    return (10 * (fromdb[i + 1] - '0') + (fromdb[i + 2] - '0')) >=
           _min_examples_needed_for_addition;

  const_IWSubstring s;
  s = fromdb.from_to(i + 1, fromdb.length() - 1);

  int u;

  if (!s.numeric_value(u) || u < 1) {
    cerr << "SubstituentIdentification::_enough_examples:invalid count stored '"
         << fromdb << "'\n";
    return 0;
  }

  return u >= _min_examples_needed_for_addition;
}

void SubstituentIdentification::_associate_substituent_with_bit(
    int r, unsigned int b, bond_type_t bt, Molecule &substituent,
    const Molecule_Specific_Temporary_Arrays &msta) {
// #define DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
  cerr << "At radius " << r << ", b = " << b << " processing "
       << substituent.unique_smiles() << endl;
#endif

  IWString usmi(
      substituent
          .unique_smiles()); // note that we add the attachment type to the end

  if (SINGLE_BOND == bt)
    ;
  else if (DOUBLE_BOND == bt)
    usmi << '=';
  else
    usmi << '#';

  auto &y = _bit[r]; // appropriate radius hash

  auto f = y.find(b);

  if (f == y.end()) // never seen this bit before
  {
    y[b] = std::unordered_map<IWString, ASubstituent, IWStringHash>();

    y[b].emplace(usmi, msta.molecule_name());
    auto qq = y[b].find(usmi);
    qq->second.set_first_molecule(msta.molecule_name());
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
    cerr << "After starting new bit " << y[b].size() << " items stored\n";
#endif
  } else // bit has been seen before
  {
    auto &usmi2subs = (*f).second;

    //  unordered_map<IWString, ASubstituent, IWStringHash>::iterator f2 =
    //  usmi2subs.find(substituent.unique_smiles());

    auto f2 = usmi2subs.find(usmi);

    if (f2 == usmi2subs.end()) {
      usmi2subs.emplace(usmi, msta.molecule_name());
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
      cerr << "For bit " << b << " new item " << substituent.unique_smiles()
           << ", name " << msta.molecule_name() << endl;
      auto qq = usmi2subs.find(usmi);
      qq->second.set_first_molecule(msta.molecule_name());
      cerr << "new key '" << qq->first << "' new value "
           << qq->second.first_molecule() << endl;
#endif
    } else {
      if (_concatenate_all_examples)
        f2->second.extra(msta.molecule_name());
      else
        (*f2).second.extra();
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
      cerr << "Extra instance of " << substituent.unique_smiles() << " now has "
           << (*f2).second.number_instances() << endl;
#endif
    }
  }

  return;
}

int SubstituentIdentification::_id_attch_pt_and_make_substituent_associations(
    Molecule &anchor, bond_type_t bt, Molecule &substituent,
    Molecule_Specific_Temporary_Arrays &msta) {
  unsigned int *b = new unsigned int[_shell_radius + 1];
  std::unique_ptr<unsigned int[]> free_b(b);

  const auto a1 = anchor.atom_with_isotope(1);

  if (INVALID_ATOM_NUMBER == a1) {
    cerr << "SubstituentIdentification::no atom with isotope 1 in fragment '"
         << anchor.smiles() << "'\n";
    return 0;
  }

  auto max_radius_formed =
      _compute_environment01(anchor, a1, _shell_radius, b, msta);

  for (auto i = 0; i <= max_radius_formed; ++i) {
    //  cerr << "Radius " << i << " bit " << b[i] << endl;
    _associate_substituent_with_bit(i, b[i], bt, substituent, msta);
  }

  return 1;
}

int convert_to_atom_number_in_parent(Molecule &m, const atom_number_t zatom) {
  // cerr << "Atom " << zatom << " was atom " << *(reinterpret_cast<const
  // atom_number_t *>(m.user_specified_atom_void_ptr(zatom))) << " in parent\n";

  return *(reinterpret_cast<const atom_number_t *>(
      m.user_specified_atom_void_ptr(zatom)));
}

int SubstituentIdentification::_compute_environment01(
    Molecule &m, const atom_number_t zatom, int rmax, unsigned int *rc,
    Molecule_Specific_Temporary_Arrays &msta) const {
  const Atom *a = m.atomi(zatom);

  const atom_number_t zip = convert_to_atom_number_in_parent(m, zatom);

  rc[0] = msta.atom_type(zip);

  if (0 == rmax)
    return 0;

  unsigned int s = rc[0] + 127516;

  msta.set_processing_status(0);
  msta.set_processing_status(zatom, PROCESSING_FINISHED);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i) {
    const Bond *b = a->item(i);

    const auto a1 = b->other(zatom);

    const auto a1p = convert_to_atom_number_in_parent(m, a1);

    int bc = msta.bond_constant(*b, zip, a1p);

    s += bc * msta.atom_type(a1p);

    msta.set_processing_status(a1, NEXT_TIME); // may not be needed
  }

  rc[1] = s;

  if (1 == rmax) // we are done
    return 1;

  return _compute_environment2x(m, 2, rmax, rc, msta);
}

int SubstituentIdentification::_compute_environment2x(
    Molecule &m, int radius, int rmax, unsigned int *rc,
    Molecule_Specific_Temporary_Arrays &msta) const {
  const auto matoms = m.natoms();

  unsigned int s = 0;

  for (auto i = 0; i < matoms; ++i) {
    if (NEXT_TIME != msta.processing_status(i))
      continue;

    msta.set_processing_status(i, PROCESSING_FINISHED);

    const auto iip = convert_to_atom_number_in_parent(m, i);

    const Atom *a = m.atomi(i);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j) {
      const Bond *b = a->item(j);

      const auto k = b->other(i);

      if (0 != msta.processing_status(k))
        continue;

      msta.set_processing_status(k, READY_TO_PROCESS);

      const auto kip = convert_to_atom_number_in_parent(m, k);

      if (_precise_fingerprints)
        s += msta.atom_type(iip) *
             (msta.bond_constant(*b, iip, kip) + 17 * msta.atom_type(kip));
      else
        s += msta.bond_constant(*b, iip, kip) * msta.atom_type(kip);
    }
  }

  if (0 == s)
    return radius - 1;

  rc[radius] = rc[radius - 1] * radius + s;

  if (radius == rmax)
    return radius;

  for (auto i = 0; i < matoms; ++i) {
    if (READY_TO_PROCESS == msta.processing_status(i))
      msta.set_processing_status(i, NEXT_TIME);
  }

  return _compute_environment2x(m, radius + 1, rmax, rc, msta);
}

int SubstituentIdentification::_divide_molecule(
    Molecule &m, const atom_number_t a1, const atom_number_t a2,
    Molecule_Specific_Temporary_Arrays &msta) {
  const bond_type_t put_back = m.btype_between_atoms(a1, a2);

  if (SINGLE_BOND != put_back)
    msta.temporarily_saturate(a1, a2);

  m.remove_bond_between_atoms(a1, a2);
  m.recompute_implicit_hydrogens(a1);
  m.recompute_implicit_hydrogens(a2);
  m.set_isotope(a1, 1);
  m.set_isotope(a2, 1);
  resizable_array_p<Molecule> components;
  m.create_components(components);
  // cerr << "From " << x << " created " << components[0]->smiles() << " " <<
  // components[1]->smiles() << endl;

  const auto atom_count_0 = components[0]->natoms();
  const auto atom_count_1 = components[1]->natoms();

  if (_ok_atom_count(atom_count_0) && atom_count_1 >= _min_residual_size)
    _id_attch_pt_and_make_substituent_associations(*(components[1]), put_back,
                                                   *(components[0]), msta);

  if (_ok_atom_count(atom_count_1) && atom_count_0 >= _min_residual_size)
    _id_attch_pt_and_make_substituent_associations(*(components[0]), put_back,
                                                   *(components[1]), msta);

  m.set_isotope(a1, 0);
  m.set_isotope(a2, 0);
  m.add_bond(a1, a2, put_back);

  if (SINGLE_BOND != put_back)
    msta.undo_temporary_saturation(a1, a2);

  return 1;
}

int SubstituentIdentification::_build_database(
    Molecule &m, IWString_and_File_Descriptor &output) {
  const auto matoms = m.natoms();

  Molecule_Specific_Temporary_Arrays msta(matoms);

  int *numbers = new int[matoms];
  std::unique_ptr<int[]> free_numbers(numbers);

  _initialise_msta(m, numbers, msta);

  Molecule mcopy(m);

  for (int i = 0; i < matoms; i++) {
    if (0 == msta.atom_type(i))
      continue;

    const Atom *a = m.atomi(i);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j) {
      const Bond *b = a->item(j);

      if (b->nrings())
        continue;

      const auto k = b->other(i);

      if (k < i)
        continue;

      _divide_molecule(mcopy, i, k, msta);
    }
  }

  static const Element *hydrogen = get_element_from_atomic_number(1);

  if (_process_hydrogen_substituents) {
    for (int i = 0; i < matoms; ++i) {
      if (0 == msta.atom_type(i))
        continue;

      const Atom *a = m.atomi(i);

      const auto acon = a->ncon();

      if (1 == acon)
        continue;

      if (0 == m.hcount(i))
        continue;

      mcopy.add(hydrogen);
      mcopy.add_bond(i, matoms, SINGLE_BOND);
      _divide_molecule(mcopy, i, matoms, msta);
      mcopy.remove_atom(matoms);
    }
  }

  return 1;
}

int SubstituentIdentification::_process_molecules(
    data_source_and_type<Molecule> &input,
    IW_STL_Hash_Set &smiles_already_found,
    IWString_and_File_Descriptor &output) {
  Molecule *m;
  while (nullptr != (m = input.next_molecule())) {
    _molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    int rc;
    if (_anchor_query.number_elements() ||
        _default_new_molecule_starting_points)
      rc = _look_for_new_substituents(*m, smiles_already_found, output);
    else
      rc = _build_database(*m, output);

    if (0 == rc)
      return 0;

    output.write_if_buffer_holds_more_than(32768);

    if (_report_progress()) {
      if (_anchor_query.number_elements() ||
          _default_new_molecule_starting_points)
        _do_create_molecules_report(cerr);
      else
        _do_build_database_report(cerr);
    }
  }

  return 1;
}

void SubstituentIdentification::_do_build_database_report(
    std::ostream &os) const {
  os << "Processed " << _molecules_read << " molecules,";
  for (auto i = 0; i <= _shell_radius; ++i) {
    const auto &bi = _bit[i];

    int s = 0;
    for (auto j = bi.cbegin(); j != bi.cend(); ++j) {
      s += (*j).second.size();
    }
    os << " rad " << i << ' ' << s;
  }
  os << "\n";
}

void SubstituentIdentification::_do_create_molecules_report(
    std::ostream &os) const {
  cerr << "Read " << _molecules_read << " molecules, "
       << _molecules_hitting_queries << " hit queries. Produced "
       << _molecules_written << ", " << _invalid_valences_ignored
       << " invalid valences ignored\n";
  if (_molecules_producing_too_many_new_molecules)
    cerr << _molecules_producing_too_many_new_molecules
         << " would have produced more than "
         << _max_molecules_per_input_molecule << " new molecules\n";

  return;
}

int SubstituentIdentification::_process_molecules(
    const char *fname, FileType input_type,
    IW_STL_Hash_Set &smiles_already_found,
    IWString_and_File_Descriptor &output) {
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

  if (_verbose > 1)
    input.set_verbose(1);

  return _process_molecules(input, smiles_already_found, output);
}

static int opendb(Db &db, const char *dbname, DBTYPE dbtype, int flags,
                  int mode) {
  int rc = db.open(NULL, dbname, NULL, dbtype, flags, mode);

  if (0 != rc) {
    cerr << "Cannot open database '" << dbname << "'\n";
    db.err(rc, "");
    return 0;
  }

  return 1;
}

static int opendb_read(Db &db, const char *dbname) {
  return opendb(db, dbname, DB_UNKNOWN, DB_RDONLY, S_IREAD);
}

static int opendb_write(Db &db, const char *dbname) {
  return opendb(db, dbname, DB_BTREE, DB_CREATE,
                S_IREAD | S_IWRITE | S_IRGRP | S_IROTH);
}

#define RADIUS_KEY "_RADIUS"

int SubstituentIdentification::_store_radius() {
  Dbt dbkey((void *)RADIUS_KEY, strlen(RADIUS_KEY));
  IWString tmp;
  tmp << _shell_radius;

  Dbt zdata(const_cast<char *>(tmp.rawdata()), tmp.length());

  const auto s = _db.put(NULL, &dbkey, &zdata, 0);
  if (0 == s)
    return 1;

  cerr << "SubstituentIdentification::_store_radius:Berkeley database Put "
          "operation failed\n";
  _db.err(s, "");
  return 0;
}

int SubstituentIdentification::_get_radius() {
  Dbt dbkey((void *)RADIUS_KEY, strlen(RADIUS_KEY));
  Dbt zdata;

  int s;
  if (_ndb > 0)
    s = _dbs[0]->get(NULL, &dbkey, &zdata, 0);
  else
    s = _db.get(NULL, &dbkey, &zdata, 0);

  if (0 != s) {
    cerr << "SubstituentIdentification::_get_radius:Berkeley database get "
            "operation failed\n";
    _db.err(s, "");
    return 0;
  }

  const_IWSubstring tmp((const char *)zdata.get_data(), zdata.get_size());

  tmp.truncate_at_first(':'); // in case databases have been concatenated

  int rc;
  if (!tmp.numeric_value(rc) || rc <= 0) {
    cerr << "SubstituentIdentification::_get_radius:invalid radius '" << tmp
         << "'\n";
    return 0;
  }

  return rc;
}

int SubstituentIdentification::_read_smiles_already_found(
    data_source_and_type<Molecule> &input,
    IW_STL_Hash_Set &smiles_already_found) {
  Molecule *m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    smiles_already_found.insert(m->unique_smiles());
  }

  return smiles_already_found.size();
}

int SubstituentIdentification::_read_smiles_already_found(
    const char *fname, IW_STL_Hash_Set &smiles_already_found) {
  const FileType input_type = discern_file_type_from_name(fname);
  if (input_type == FILE_TYPE_INVALID) {
    cerr << "read_smiles_already_found:do not know how to process '" << fname
         << "'\n";
    return 0;
  }

  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open already formed molecules file '" << fname << "'\n";
    return 0;
  }

  return _read_smiles_already_found(input, smiles_already_found);
}

int SubstituentIdentification::operator()(int argc, char **argv) {
  Command_Line cl(
      argc, argv,
      "vA:E:i:g:ld:R:P:M:m:ct:s:q:ab:r:C:H:fw:pIBx:u:kYL:jyezX:hUn:V:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    _usage(1);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, _verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      _usage(5);
    }
  } else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E')) {
    if (!process_elements(cl, _verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1,
                                                               'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      _usage(32);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;

    if (_verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;

    if (_verbose)
      cerr << "Will remove all chirality\n";
  }

  if (cl.option_present('a')) {
    _only_produce_molecules_at_biggest_radius = 1;

    if (_verbose)
      cerr << "Will only produce molecules at the largest radius for which a "
              "match is found\n";
  }

  if (cl.option_present('b')) {
    _only_add_bond = 0;

    IWString b;
    for (int i = 0; cl.value('b', b, i); ++i) {
      b.to_lowercase();

      if ('s' == b[0])
        _only_add_bond |= SINGLE_BOND;
      else if ('d' == b[0])
        _only_add_bond |= DOUBLE_BOND;
      else if ('t' == b[0])
        _only_add_bond |= TRIPLE_BOND;
      else {
        cerr << "Unrecognised bond type specification '" << b << "'\n";
        _usage(1);
      }
    }

    if (_verbose)
      cerr << "Bond mask " << _only_add_bond << endl;
  }

  if (cl.option_present('t')) {
    if (!_report_progress.initialise(cl, 't', _verbose))
      return 1;
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', _shell_radius) || _shell_radius < 0) {
      cerr << "The shell radius option (-R) must be a valid radius\n";
      _usage(2);
    }

    if (_verbose)
      cerr << "Shell radius set to " << _shell_radius << endl;
  } else
    _shell_radius = 1;

  if (cl.option_present('r')) {
    if (!cl.value('r', _min_shell_radius) || _min_shell_radius < 0) {
      cerr << "The minimum shell radius (-r) must be a non negative whole "
              "number\n";
      _usage(1);
    }

    //  no further checking of the shell radius until we can read what is in the
    //  database
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', _max_atoms_in_product) || _max_atoms_in_product < 2) {
      cerr << "The max atoms in product molecule option (-C) must be a whole "
              "+ve number\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "Will discard product molecules with more than "
           << _max_atoms_in_product << " atoms\n";
  }

  if (cl.option_present('p')) {
    _write_parent_molecule = 1;

    if (_verbose)
      cerr << "Will write the parent molecule\n";
  }

  if (cl.option_present('I')) {
    _remove_isotopes_from_product = 1;

    if (_verbose)
      cerr << "Will remove isotopic labels from product molecules\n";
  }

  _bit =
      new std::unordered_map<unsigned int, Usmi2Substituent>[_shell_radius + 1];

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!_atom_typing_specification.build(p)) {
      cerr << "Cannot initialise attachment atom type specification '" << p
           << "'\n";
      return 1;
    }
  } else {
    _atom_typing_specification.build("UST:AFUY");
  }

  if (cl.option_present('M')) {
    if (!cl.value('M', _max_substituent_size) || _max_substituent_size < 1) {
      cerr << "The max substituent size option (-M) must be a whole +ve "
              "number\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "Will only examine substituents with " << _max_substituent_size
           << " or fewer atoms\n";
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', _min_substituent_size) || _min_substituent_size < 1) {
      cerr << "The minimum substituent size (-m) must be a whole +ve number\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "Will ignore substituents with fewer than "
           << _min_substituent_size << " atoms\n";
  }

  if (cl.option_present('L')) {
    if (!cl.value('L', _max_atoms_lost_during_build) ||
        _max_atoms_lost_during_build < 0) {
      cerr << "The max atoms lost during build option (-L) must be a whole +ve "
              "number\n";
      _usage(1);
    }
    if (_verbose)
      cerr << "Parent molecules can lose a max of "
           << _max_atoms_lost_during_build << " atoms during building\n";
  }

  if (cl.option_present('j')) {
    _write_fragments_added = 1;
    if (_verbose)
      cerr << "Will write the fragment added during molecule creation\n";
  }

  if (cl.option_present('y')) {
    _default_new_molecule_starting_points = 1;
    if (_verbose)
      cerr << "Will use default atom specifications for growing molecules\n";
  }

  if (cl.option_present('X')) {
    const_IWSubstring e = cl.string_value('X');

    if (!_expt.read_data(e)) {
      cerr << "Cannot read experimental data data from '" << e << "'\n";
      return 1;
    }

    if (_verbose)
      cerr << _expt.size() << " experimental values read from '" << e << "'\n";
  }

  if (cl.option_present('e')) {
    _concatenate_all_examples = 1;

    if (_verbose)
      cerr << "Will store all examples of a substituent - not just first\n";
  }

  if (cl.option_present('z')) {
    _strip_leading_zeros = 1;

    if (_verbose)
      cerr << "Will strip leading zeros from identifiers\n";
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, _anchor_query, _verbose, 'q')) {
      cerr << "Cannot process anchor point queries (-q)\n";
      return 2;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      Substructure_Query *t = new Substructure_Query;
      if (!t->create_from_smarts(s)) {
        cerr << "Cannot parse smarts '" << s << "'\n";
        delete t;
        return 2;
      }
      _anchor_query.add(t);
    }
  }

  for (int i = 0; i < _anchor_query.number_elements(); ++i) {
    _anchor_query[i]->set_find_one_embedding_per_atom(1);
  }

  if (cl.option_present('k')) {
    _break_molecule_at_first_two_matched_atoms = 1;

    if (_verbose)
      cerr << "Will break the input molecule at the first two matched atoms\n";
  }

  if (cl.option_present('x')) {
    if (!cl.value('x', _max_molecules_per_input_molecule) ||
        _max_molecules_per_input_molecule < 1) {
      cerr << "The maximum number of molecules produced per input molecule "
              "must be a valid +ve number\n";
      return 2;
    }

    if (_verbose)
      cerr << "Will produce a maximum of " << _max_molecules_per_input_molecule
           << " molecules per input molecule\n";
  }

  if (cl.option_present('u')) {
    if (!cl.value('u', _min_examples_needed_for_addition) ||
        _min_examples_needed_for_addition < 1) {
      cerr << "The minimum number of example structures needed for addition "
              "(-u) must be a whole +ve number\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "Will only add a substituent if there are "
           << _min_examples_needed_for_addition
           << " or more examples in the db\n";
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', _matched_atoms_to_process) ||
        _matched_atoms_to_process < 1) {
      cerr << "The number of matched atoms to process (-n) must be a whole +ve "
              "number\n";
      _usage(1);
    }
    if (_verbose)
      cerr << "Will process the first " << _matched_atoms_to_process
           << " matched atoms in each query\n";
  }

  if (cl.option_present('U')) {
    _apply_atom_map_labels = 1;
    if (_verbose)
      cerr << "Will apply atom map labels\n";
  }

  if (cl.option_present('B'))
    _remove_chirality = 1;

  if (0 == _anchor_query.number_elements() && !cl.option_present('B') &&
      0 == _default_new_molecule_starting_points) {
    cerr << "No queries specified (make new molecules), but the -B option for "
            "building not specified, cannot continue\n";
    _usage(2);
  }

  // Not building a database, but did not specify a max substituent size, just
  // give them whatever is in the database

  if (!cl.option_present('B') && !cl.option_present('M'))
    _max_substituent_size = std::numeric_limits<int>::max();

  if (cl.option_present('h')) {
    _process_hydrogen_substituents = 1;
    if (_verbose)
      cerr << "Will also generate Hydrogen substituents\n";
  }

  if (cl.option_present('H')) {
    if (!process_queries(cl, _substituents_must_contain, _verbose, 'H')) {
      cerr << "Cannot process queries substituents must contain (-H)\n";
      return 2;
    }
  }

  if (_verbose && _anchor_query.number_elements())
    cerr << "Defined " << _anchor_query.number_elements()
         << " queries to identify attachment points\n";

  if (cl.option_present('f')) {
    _precise_fingerprints = 0;

    if (_verbose)
      cerr << "Will generate fingerprints that allow imprecise matching\n";
  }

  if (cl.option_present('Y')) {
    _make_implicit_hydrogens_explicit = 1;

    if (_verbose)
      cerr << "Will make implicit Hydrogens explicit\n";
  }

  if (cl.option_present('w')) {
    if (!cl.value('w', _min_residual_size) || _min_residual_size < 1) {
      cerr << "The minimum residual atom count option (-w) must be a whole +ve "
              "number\n";
      _usage(2);
    }

    if (_verbose)
      cerr << "Will not produce a fragment if the residual contains fewer than "
           << _min_residual_size << " atoms\n";
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify database to build via the -d option\n";
    _usage(1);
  }

  auto dcount = cl.option_count('d');

  if (dcount > 1 && 0 == _anchor_query.number_elements()) {
    cerr << "Can only build one database at a time\n";
    return 1;
  }

  if (1 == dcount && 0 == _anchor_query.number_elements() &&
      0 == _default_new_molecule_starting_points) // great, building a database
                                                  // and no queries
  {
    const char *dbname = cl.option_value('d');
    if (!opendb_write(_db, dbname))
      return 1;

    if (!_store_radius())
      return 1;
  } else if (1 == dcount) // lookup, but just one database
  {
    const char *dbname = cl.option_value('d');

    if (!opendb_read(_db, dbname))
      return 1;
  } else // must be a lookup from multiple databases
  {
    _ndb = dcount;
    _dbs = new Db *[_ndb];

    for (auto i = 0; i < _ndb; ++i) {
      const char *dbname = cl.option_value('d', i);

      _dbs[i] = new Db(0, DB_CXX_NO_EXCEPTIONS);

      if (!opendb_read(*_dbs[i], dbname))
        return 1;
    }

    if (_verbose)
      cerr << "Opened " << _ndb << " databases for lookup\n";
  }

  if (cl.option_present('r')) {
    auto t = _get_radius();
    if (t < 1) {
      cerr << "SubstituentIdentification:no cannot retrieve radius from DB, "
              "beware\n";
      t = _min_shell_radius; // just use what the user entered
    }

    if (_min_shell_radius > t) {
      cerr << "Min shell radius " << _min_shell_radius
           << " too large, database built with max radius " << t << endl;
      return 1;
    }

    _shell_radius = t;

    if (_verbose)
      cerr << "Will only look at shells of radius " << _min_shell_radius
           << " or more\n";

    if (_shell_radius < _min_shell_radius)
      _shell_radius = _min_shell_radius;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      _usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (!all_files_recognised_by_suffix(cl))
    return 4;

  set_copy_name_in_molecule_copy_constructor(1);
  set_copy_user_specified_atom_void_ptrs_during_create_subset(1);

  IW_STL_Hash_Set smiles_already_found;

  if (cl.option_present('V')) {
    const char *v = cl.option_value('V');

    if (!_read_smiles_already_found(v, smiles_already_found)) {
      cerr << "Cannot read already formed smiles '" << v << "'\n";
      return 0;
    }

    if (_verbose)
      cerr << "Read " << smiles_already_found.size()
           << " already formed molecules to avoid from '" << v << "'\n";
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    _usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!_process_molecules(cl[i], input_type, smiles_already_found, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (_verbose) {
    cerr << "Read " << _molecules_read << " molecules\n";
    if (_anchor_query.number_elements() ||
        _default_new_molecule_starting_points)
      _do_create_molecules_report(cerr);
    else {
      for (auto i = 0; i <= _shell_radius; ++i) {
        const auto &b = _bit[i];

        cerr << "At radius " << i << " found " << b.size() << " items\n";
      }
      _do_build_database_report(cerr);
      cerr << "Stored " << _keys_stored << " bits, " << _pairs_stored
           << " pairs\n";
    }
  }

  if (0 == _anchor_query.number_elements() &&
      0 == _default_new_molecule_starting_points)
    _write_in_memory_hashes_to_database();

  return rc;
}

int main(int argc, char **argv) {
  prog_name = argv[0];

  SubstituentIdentification SubstituentIdentification;

  return SubstituentIdentification(argc, argv);
}
