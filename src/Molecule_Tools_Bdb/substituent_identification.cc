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
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream_impl_lite.h"

#define IWQSORT_IMPLEMENTATION 1
#define IWQSORT_FO_IMPLEMENTATION 1

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/mformula.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools_Bdb/substituent_identification.pb.h"
#else
#include "substituent_identification.pb.h"
#endif

#include "db_cxx.h"

using std::cerr;

const char* prog_name = nullptr;

#define PROCESSING_FINISHED 1
#define READY_TO_PROCESS 2
#define NEXT_TIME 3

class Molecule_Specific_Temporary_Arrays {
 public:
  const int _matoms;
  int* _atype;
  int* _aromatic_bond;
  int* _processing_status;
  int* _to_remove;
  IWString _molecule_name;

 public:
  Molecule_Specific_Temporary_Arrays(int matoms);
  ~Molecule_Specific_Temporary_Arrays();

  void initialise_aromatic_bond_array(Molecule& m);

  void
  set_processing_status(int ndx, int s) {
    _processing_status[ndx] = s;
  }

  void
  set_processing_status(int s) {
    set_vector(_processing_status, _matoms, s);
  }

  int* to_remove();

  int
  processing_status(int s) const {
    return _processing_status[s];
  }

  int
  atom_type(int a) const {
    return _atype[a];
  }

  void temporarily_saturate(atom_number_t a1, atom_number_t a2);
  void undo_temporary_saturation(atom_number_t a1, atom_number_t a2);

  int bond_constant(const Bond& b, const atom_number_t, const atom_number_t) const;

  const IWString&
  molecule_name() const {
    return _molecule_name;
  }
};

Molecule_Specific_Temporary_Arrays::Molecule_Specific_Temporary_Arrays(int matoms)
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

int
Molecule_Specific_Temporary_Arrays::bond_constant(const Bond& b, const atom_number_t a1,
                                                  const atom_number_t a2) const {
  if (_aromatic_bond[a1 * _matoms + a2]) {
    return 11;
  }

  if (b.is_single_bond()) {
    return 3;
  }

  if (b.is_double_bond()) {
    return 5;
  }

  return 7;
}

void
Molecule_Specific_Temporary_Arrays::initialise_aromatic_bond_array(Molecule& m) {
  m.compute_aromaticity_if_needed();

  const auto ne = m.nedges();

  for (auto i = 0; i < ne; ++i) {
    const Bond* b = m.bondi(i);

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

int*
Molecule_Specific_Temporary_Arrays::to_remove() {
  std::fill_n(_to_remove, _matoms, 0);

  return _to_remove;
}

/*
  This is a big hack and quite unstable. By looking into the atom typing code we
  observe that the numeric difference between a saturated atom and an
  unsaturated atom is 1. Very bad for stability, but I cannot think of a better
  way of doing this. In case of problems with multiple bonds, look here first
*/

void
Molecule_Specific_Temporary_Arrays::temporarily_saturate(atom_number_t a1,
                                                         atom_number_t a2) {
  _atype[a1]++;
  _atype[a2]++;

  return;
}

void
Molecule_Specific_Temporary_Arrays::undo_temporary_saturation(atom_number_t a1,
                                                              atom_number_t a2) {
  _atype[a1]--;
  _atype[a2]--;
}

// There are certain common things used throughout lookups.
class UsedDuringLookups {
  private:
    // The unique smiles of product molecules already formed.
    IW_STL_Hash_Set _product_already_formed;

    // The unique smiles of fragments that have already been added to the
    // current molecule.
    // This is for when dealing with multiple databases do avoid forming
    // duplicates;
    absl::flat_hash_set<std::string> _reagent_already_used_this_molecule;

    // The names of substituents that have been added to the current
    // molecule. Avoids forming duplicates in the case of concatenated
    // databases.
    IW_STL_Hash_Set _name_already_found;

    uint32_t _molecules_formed_from_current_molecule;

    IWString _smiles_of_fragment_lost;

  public:
    UsedDuringLookups();

    // Note that we implement the idea of products to NOT form by reading
    // these molecules into _product_already_formed.
    int ReadFragmentsToIgnore(const char* fname);
    int ReadFragmentsToIgnore(data_source_and_type<Molecule>& input);

    void PrepareForNewMolecule();

    void ClearReagentNames() {
      _reagent_already_used_this_molecule.clear();
    }

    IWString& smiles_of_fragment_lost() {
      return _smiles_of_fragment_lost;
    }
    const IWString& smiles_of_fragment_lost() const {
      return _smiles_of_fragment_lost;
    }

    uint32_t molecules_formed_from_current_molecule() const {
      return _molecules_formed_from_current_molecule;
    }
    uint32_t& molecules_formed_from_current_molecule() {
      return _molecules_formed_from_current_molecule;
    }
    void AnotherMoleculeFormed() {
      ++_molecules_formed_from_current_molecule;
    }

    IW_STL_Hash_Set& product_already_formed() {
      return _product_already_formed;
    }

    // Return true if the non isotopic unique smiles is in
    // _product_already_formed
    bool ProductMoleculeSeenBefore(Molecule& m);

    bool ReagentAlreadyUsedThisMolecule(const std::string& smiles);
};

UsedDuringLookups::UsedDuringLookups() {
  _molecules_formed_from_current_molecule = 0;
}

void
UsedDuringLookups::PrepareForNewMolecule() {
  _molecules_formed_from_current_molecule = 0;
  _name_already_found.clear();
}

int
UsedDuringLookups::ReadFragmentsToIgnore(const char* fname) {
  const FileType input_type = discern_file_type_from_name(fname);
  if (input_type == FILE_TYPE_INVALID) {
    cerr << "UsedDuringLookups::ReadFragmentsToIgnore:do not know how to process '" << fname << "'\n";
    return 0;
  }

  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open already formed molecules file '" << fname << "'\n";
    return 0;
  }

  return ReadFragmentsToIgnore(input);
}

int
UsedDuringLookups::ReadFragmentsToIgnore(data_source_and_type<Molecule>& input) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    // Cannot call preprocess here.

    _product_already_formed.insert(m->unique_smiles());
  }

  return _product_already_formed.size();
}

bool
UsedDuringLookups::ProductMoleculeSeenBefore(Molecule& m) {
  Molecule mcopy(m);
  mcopy.transform_to_non_isotopic_form();

  const IWString& usmi = mcopy.unique_smiles();

  // cerr << "Checking " << usmi << '\n';

  if (_product_already_formed.contains(usmi)) {
    return true;
  }

  _product_already_formed.insert(usmi);

  return false;  // never seen this before
}

bool
UsedDuringLookups::ReagentAlreadyUsedThisMolecule(const std::string& smiles) {
  auto iter = _reagent_already_used_this_molecule.find(smiles);
  if (iter != _reagent_already_used_this_molecule.end()) {
    return true;
  }

  _reagent_already_used_this_molecule.insert(smiles);

  return false;
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
  uint32_t _number_instances;

  // If we are concatenating the names of the examplars, we can limit the number
  // of names concatenated. Otherwise the database can get too large.
  static uint32_t _max_instances;

 public:
  ASubstituent(const IWString& s) : _first_molecule(s) {
    _number_instances = 1;
  }

  void SetMaxInstancesToConcatenate(uint32_t s) {
    _max_instances = s;
  }

  void
  extra() {
    _number_instances++;
  }

  void extra(const IWString&);

  const IWString&
  first_molecule() const {
    return _first_molecule;
  }

  void
  set_first_molecule(const IWString& s) {
    _first_molecule = s;
  }

  uint32_t
  number_instances() const {
    return _number_instances;
  }

  void write_first_and_instances(IWString&) const;
};

uint32_t ASubstituent::_max_instances = 0;

int
SetMaxInstancesToConcatenate(const const_IWSubstring& buffer, int verbose) {
  static ASubstituent asubstituent("NOTUSED");

  uint32_t s;
  if (! buffer.numeric_value(s)) {
    cerr << "SetMaxInstancesToConcatenate:the max number of instances to concatenate must be a whole +ve number\n";
    return 0;
  }

  asubstituent.SetMaxInstancesToConcatenate(s);
  if (verbose) {
    cerr << "Will store a max of " << s << " examplars in concatenated names\n";
  }

  return 1;
}

/*
this causes all kinds of problems generating ambiguous name re4solutions
template <typename T>
T &
operator << (T & os, const ASubstituent & s)
{
  os << s.first_molecule() << ':' << s.number_instances();

  return os;
}*/

void
ASubstituent::write_first_and_instances(IWString& s) const {
  s << _first_molecule << ':' << _number_instances;

  // cerr << "First '" << _first_molecule << "' " << _number_instances << '\n';

  return;
}

void
ASubstituent::extra(const IWString& mname) {
  _number_instances++;
  if (_max_instances > 0 && _number_instances > _max_instances) {
    return;
  }

  _first_molecule << _inter_molecule_name_separator << mname;

  return;
}

class Substituent_Hash {
 private:
 public:
#if defined(IW_INTEL_COMPILER)

  static const size_t bucket_size = 4;
  static const size_t min_buckets = 8;
  bool operator()(const ASubstituent&, const IWString&) const;

#endif

  size_t operator()(const ASubstituent&) const;
};

size_t
Substituent_Hash::operator()(const ASubstituent& s) const {
  return 1;
}

typedef std::unordered_map<IWString, ASubstituent, IWStringHash> Usmi2Substituent;

/*
  When doing MMP studies, there are a number of things we need to keep track of
*/

class MMP_Related {
 private:
 public:
  MMP_Related();
};

MMP_Related::MMP_Related() {
  return;
}

// When keeping track of a fragment that has been lost we may need
// to keep track of various attributes of it.
struct FragmentLost {
  IWString smiles;
  int natoms;
  // If this is needed, it will be allocated and computed.
  std::unique_ptr<mformula::MFormula> formula;

  public:
    FragmentLost() {
      natoms = -1;
    }

  void EnsureNatomsComputed() {
    if (natoms <= 0) {
      natoms = lillymol::count_atoms_in_smiles(smiles);
    }
  }

  void SetHydrogen() {
    smiles = "H";
    natoms = 1;
  }

  void Set(Molecule& m) {
    smiles = m.unique_smiles();
    natoms = m.natoms();
  }
};

class SubstituentIdentification {
 private:
  int _verbose;

  uint64_t _molecules_read;

  Chemical_Standardisation _chemical_standardisation;

  int _reduce_to_largest_fragment;

  int _remove_chirality;

  // During building.
  //  When we chop off a substituent, what is the smallest fragment size
  //  that can be left. For example, we would not want to consider benzene as
  //  a substituent of oxygen in phenol

  int _min_residual_size;

  // During database building, a limit on what gets stored.
  // During replacement, the maximum size fragment size that will be added.
  int _max_substituent_size;

  // During replacemet, the maximum number of extra atoms allowed in going from
  // an existing substituent to the replacement.
  // And a limit on the number of atoms by which a replacement can shrink from the
  // starting number of atoms.
  uint32_t _max_extra_atoms_added = 0;
  uint32_t _max_atoms_lost_during_addition = 0;

  // The maximum number of bonds from the attachment point to the end of the molecule.
  uint32_t _max_length_substituent = 0;

  int _process_hydrogen_substituents;

  //  When we are building a molecule, we can impose a limit on the number of
  //  atoms that are stripped off. Otherwise we get non useful results
  int _max_atoms_lost_during_build;

  // We can optionally store serialized protos.
  int _store_protos;

  // We can optionally write as textproto.
  IWString_and_File_Descriptor _stream_for_textproto;

  int _shell_radius;

  //  Do next shell atoms keep track of where they are attached

  int _precise_fingerprints;

  int _make_implicit_hydrogens_explicit;

  Atom_Typing_Specification _atom_typing_specification;

  std::unordered_map<uint32_t, Usmi2Substituent>* _bit;

  Report_Progress _report_progress;

  //  Some characteristics of the database being created

  uint32_t _pairs_stored;
  uint32_t _keys_stored;

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

  uint32_t _only_add_bond;

  // If any of these queries match, the product is OK.
  // but see all_must_have_queries_must_match below.
  resizable_array_p<Substructure_Query> _replacements_must_contain;
  resizable_array_p<Substructure_Query> _replacements_must_not_contain;


  // If set then all the queries in _replacements_must_contain array must
  // match - rather than any of them matching which is the default.
  int _all_must_have_queries_must_match;

  int _invalid_valences_ignored;

  int _max_atoms_in_product;

  int _molecules_hitting_queries;
  int _molecules_written;

  int _write_parent_molecule;

  int _remove_isotopes_from_product;

  int _apply_atom_map_labels;

  uint32_t _max_molecules_per_input_molecule;
  uint32_t _molecules_producing_too_many_new_molecules;

  Accumulator_Int<int> _acc_new_molecules_produced;

  uint32_t _min_examples_needed_for_addition;

  int _write_fragments_added;

  int _concatenate_all_examples;

  int _strip_leading_zeros;

  int _matched_atoms_to_process;

  // We can open any number of databases.
  Db** _dbs;
  int _ndb;

  Numeric_Data_From_File<float> _expt;

  int _output_is_textproto;
  int _prepend_smiles_to_textproto;

  // set if we are processing matched pairs.
  int _do_matched_pairs;

  // private functions

  void _default_values();
  void _usage(int rc);
  void _preprocess(Molecule& m);

  int SetupTextProtoStream(const const_IWSubstring fname);

  int _ok_atom_count(const int matoms) const;
  int OkSubstituentSize(const std::string& smiles, int& natoms, int atoms_in_fragment_lost) const;
  int _matches_one_of_substituents_must_contain(Molecule& m);
  int _matches_replacement_must_not_have_query(Molecule& m);
  int OkSubstructureQueries(Molecule& m);

  int _store_radius();
  int _get_radius();

  int _get_matched_atoms_to_process(const Set_of_Atoms*);

  int _read_smiles_already_found(const char* fname,
                                 IW_STL_Hash_Set& smiles_already_found);
  int _read_smiles_already_found(data_source_and_type<Molecule>& input,
                                 IW_STL_Hash_Set& smiles_already_found);

  //  Functions used during replacement

  int _check_already_made(const Molecule& m, IW_STL_Hash_Set& already_processed) const;

  uint32_t _enough_examples(const const_IWSubstring& fromdb) const;

  int MaybeRemoveTrailingBondOkH(Molecule& m,
                           atom_number_t zatom,
                           const_IWSubstring& smiles,
                           bond_type_t& bt) const;

  int _look_for_new_substituents(Molecule& m, IW_STL_Hash_Set& smiles_already_found,
                                 IWString_and_File_Descriptor& output);
  int _look_for_new_substituents(Molecule& m, 
                    UsedDuringLookups& lookup_data,
                    IWString_and_File_Descriptor& output);
  int _look_for_new_substituents(Molecule& m, atom_number_t zatom,
                                 FragmentLost& fragment_lost,
                                 IW_STL_Hash_Set& smiles_already_found,
                                 Molecule_Specific_Temporary_Arrays& msta,
                                 IWString_and_File_Descriptor& output);
  int _look_for_new_substituents(Molecule& m, atom_number_t zatom,
                                 FragmentLost& fragment_lost,
                                 IW_STL_Hash_Set& smiles_already_found,
                                 struct DBKey& dbkey, uint32_t& new_molecules_produced,
                                 IWString_and_File_Descriptor& output);
  int _look_for_new_substituents(
    Molecule& m, atom_number_t zatom, 
    UsedDuringLookups& lookup_data,
    struct DBKey& dbkey,
    IWString_and_File_Descriptor& output);
  int _look_for_new_substituents(
    Molecule& m, atom_number_t zatom, 
    UsedDuringLookups& lookup_data,
    Molecule_Specific_Temporary_Arrays& msta,
    IWString_and_File_Descriptor& output);
  int _look_for_new_substituents_db(Molecule& m, atom_number_t zatom,
                                    FragmentLost& fragment_lost,
                                    IW_STL_Hash_Set& smiles_already_found,
                                    const struct DBKey& dbkey, Db& db, Dbt& zkey,
                                    IW_STL_Hash_Set& already_processed,
                                    IWString_and_File_Descriptor& output);
  int _look_for_new_substituents_db(
    Molecule& m, atom_number_t zatom,
    UsedDuringLookups& lookup_data, 
    const struct DBKey& rad_and_bit, Db& db,
    Dbt& dbkey, IWString_and_File_Descriptor& output);

  int OkFragmentSize(const std::string& replacement_smiles, FragmentLost& fragment_lost);

  int _form_new_molecule(Molecule& m, atom_number_t zatom,
                         const FragmentLost& fragment_lost,
                         IW_STL_Hash_Set& smiles_already_found,
                         const struct DBKey& rad_and_bit, const const_IWSubstring& buffer,
                         IW_STL_Hash_Set& already_processed,
                         IWString_and_File_Descriptor& output);
  int AddFragment(Molecule& m, atom_number_t zatom,
                  Molecule& frag, atom_number_t fragment_attachment_point,
                  bond_type_t bt) const;
  int FormNewMolecules(Molecule& m, atom_number_t zatom,
                                               FragmentLost& fragment_lost,
                                               IW_STL_Hash_Set& smiles_already_found,
                                               const struct DBKey& rad_and_bit,
                                               const substituent_identification::Replacements& proto,
                                               IW_STL_Hash_Set& already_processed,
                                               IWString_and_File_Descriptor& output);
  int FormNewMolecules(Molecule& m, atom_number_t zatom,
                       const struct DBKey& rad_and_bit,
                       UsedDuringLookups& lookup_data,
                       const substituent_identification::Replacements& proto,
                       IWString_and_File_Descriptor& output);
  int FormNewMoleculeAtAtomWithReplacement(Molecule& m, atom_number_t zatom,
                                               FragmentLost& fragment_lost,
                                               IW_STL_Hash_Set& smiles_already_found,
                                               const struct DBKey& rad_and_bit,
                                               const substituent_identification::Replacement& proto,
                                               IW_STL_Hash_Set& already_processed,
                                               IWString_and_File_Descriptor& output);
  int FormNewMoleculeAtAtom(Molecule& m, atom_number_t zatom,
                      const struct DBKey& rad_and_bit,
                      UsedDuringLookups& lookup_data,
                      const substituent_identification::Replacement& proto,
                      IWString_and_File_Descriptor& output);

  int _matched_pairs_qsar_by_radius(Molecule& m, const atom_number_t zatom,
                          FragmentLost& fragment_lost,
                          IW_STL_Hash_Set& smiles_already_found,
                          Molecule_Specific_Temporary_Arrays& msta, MMP_Related& mmp,
                          IWString_and_File_Descriptor& output);
  int _matched_pairs_qsar_across_bond(Molecule& m, const atom_number_t a1,
                                      const atom_number_t a2,
                                      IW_STL_Hash_Set& smiles_already_found,
                                      Molecule_Specific_Temporary_Arrays& msta,
                                      MMP_Related& mmp,
                                      IWString_and_File_Descriptor& ouput);

  int _matched_pairs_qsar(Molecule& m, IWString_and_File_Descriptor& output);
  int _matched_pairs_qsar(Molecule& m, atom_number_t zatom,
                          FragmentLost& fragment_lost,
                          IW_STL_Hash_Set& smiles_already_found,
                          const struct DBKey& rad_and_bit, const_IWSubstring& fromdb,
                          IW_STL_Hash_Set& already_processed,
                          IWString_and_File_Descriptor& output);
  int _matched_pairs_qsar_db(Molecule& m, const atom_number_t zatom,
                             FragmentLost& fragment_lost,
                             IW_STL_Hash_Set& smiles_already_found, Dbt& dbkey,
                             int radius, MMP_Related& mmp,
                             IWString_and_File_Descriptor& output);
  int ListMatchedPairs(Molecule& m,
                atom_number_t zatom,
                FragmentLost& fragment_lost,
                const substituent_identification::Replacement& replacement,
                int radius,
                MMP_Related& mmp,
                IWString_and_File_Descriptor& output);

  atom_number_t _do_break_across_first_two_matched_atoms(
      Molecule& m, atom_number_t a1, atom_number_t a2,
      Molecule_Specific_Temporary_Arrays& msta, FragmentLost& fragment_lost);
  int _identify_fragment_being_removed(Molecule& m, const atom_number_t a0,
                                       const int* to_remove,
                                       FragmentLost& fragment_lost) const;
  int _identify_non_ring_single_bonds(Molecule& m, int* process_these_pairs) const;
  int _try_bond_breaking_across(const Molecule& m, const atom_number_t a0,
                                const atom_number_t a1, const int* to_remove,
                                IW_STL_Hash_Set& smiles_already_found,
                                Molecule_Specific_Temporary_Arrays& msta,
                                IWString_and_File_Descriptor& output);
  int _try_bond_breaking_across(const Molecule& m, const atom_number_t a0,
                                const atom_number_t a1,
                                IW_STL_Hash_Set& smiles_already_found,
                                Molecule_Specific_Temporary_Arrays& msta,
                                IWString_and_File_Descriptor& output);
  int MoleculeIsRejected(Molecule& m, Molecule& frag,
                UsedDuringLookups& lookup_data);

  //  Functions used during database builds

  void _do_build_database_report(std::ostream& os) const;
  void _do_create_molecules_report(std::ostream& os) const;

  void _associate_substituent_with_bit(int r, uint32_t b, bond_type_t bt,
                                       Molecule& substituent,
                                       const Molecule_Specific_Temporary_Arrays& msta);
  int _compute_environment2x(Molecule& m, int radius, int rmax, uint32_t* rc,
                             Molecule_Specific_Temporary_Arrays& msta) const;
  int _compute_environment01(Molecule& m, const atom_number_t zatom, int rmax,
                             uint32_t* rc,
                             Molecule_Specific_Temporary_Arrays& msta) const;

  void _initialise_msta(Molecule& m, int* number,
                        Molecule_Specific_Temporary_Arrays& msta);

  int _id_attch_pt_and_make_substituent_associations(
      Molecule& m1, bond_type_t, Molecule& m2, Molecule_Specific_Temporary_Arrays& msta);
  int _divide_molecule(Molecule& m, const atom_number_t a1, const atom_number_t a2,
                       Molecule_Specific_Temporary_Arrays& msta);
  int _build_database(Molecule& m, IWString_and_File_Descriptor& output);
  int _process_molecules(data_source_and_type<Molecule>& input,
                         IW_STL_Hash_Set& smiles_already_found,
                         IWString_and_File_Descriptor& output);
  int _process_molecules(data_source_and_type<Molecule>& input,
                                              UsedDuringLookups& lookup_data,
                                              IWString_and_File_Descriptor& output);
  int _process_molecules(const char* fname, FileType input_type,
                                              UsedDuringLookups& lookup_data,
                                              IWString_and_File_Descriptor& output);
  int _process_molecules(const char* fname, FileType input_type,
                         IW_STL_Hash_Set& smiles_already_found,
                         IWString_and_File_Descriptor& output);

  int _write_in_memory_hashes_to_database();
  int _write_in_memory_hash_to_database(
      int radius, const std::unordered_map<uint32_t, Usmi2Substituent>& b);
  int _write_substituents_associated_with_bit(const struct DBKey& dbkey,
                                              const Usmi2Substituent& u);
  int StoreProto(const struct DBKey& dbkey, const Usmi2Substituent& u);

  int WriteTextProto(Molecule& m,
                const substituent_identification::Replacement& proto,
                Molecule& frag,
                const struct DBKey& rad_and_bit,
                FragmentLost& fragment_list,
                IWString_and_File_Descriptor& output) const;
  int WriteTextProto(const struct DBKey& dbkey,
               const substituent_identification::Replacements& proto,
               IWString_and_File_Descriptor & output);

  int AlreadyMadeEnoughMolecules(const UsedDuringLookups& lookup_data) const;

 public:
  SubstituentIdentification();
  ~SubstituentIdentification();

  int operator()(int argc, char** argv);
};

SubstituentIdentification::SubstituentIdentification() {
  _default_values();
}

void
SubstituentIdentification::_default_values() {
  _verbose = 0;
  
  _molecules_read = 0;

  _reduce_to_largest_fragment = 1;  // always

  _remove_chirality = 0;

  _min_residual_size = 3;

  _max_substituent_size = 1;

  _max_extra_atoms_added = 0;
  _max_atoms_lost_during_addition = 0;

  _process_hydrogen_substituents = 0;

  _max_atoms_lost_during_build = std::numeric_limits<int>::max();

  _precise_fingerprints = 1;

  _bit = nullptr;

  _shell_radius = 0;

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

  _make_implicit_hydrogens_explicit = 0;

  _all_must_have_queries_must_match = 0;

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

  _store_protos = 1;

  _output_is_textproto = 0;
  _prepend_smiles_to_textproto = 0;

  _do_matched_pairs = 0;

  return;
}

SubstituentIdentification::~SubstituentIdentification() {
  if (nullptr != _bit) {
    delete[] _bit;
  }

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
  cerr << R"(Identifies substituents and their atomic contexts as EC type fingerprints.
First step is to scan existing collections and build a database of substituents and contexts.

Once that database is built, new molecules can be generated by looking up atomic contexts and replacing or substituting
substituents that have been found in the same context.
 -d <dbname>   Berkeley database for substituent data. Just one during builds, any number during replacement.
 -B            flag to enable building
 -R <bonds>    radius for circular fingerprint from attachment point
 -w <natoms>   during building, min number of atoms left after removing substituent
 -e <max>      during building, concatenate as many as <max> names of exemplar molecules (generate large db)
 -h            during building, discern Hydrogen as a substituent

Options used during addition/replacement
 -m <natoms>   min number of atoms in a substituent. Use -n for as many as 'n' fewer atoms than an existing sidechain.
 -M <natoms>   max number of atoms in a substituent. Use +n for as many as 'n' more atoms than an existing sidechain.
 -q <query>    during addition, query  describing anchor atoms needing new substituents
 -s <smarts>   during addition, smarts describing anchor atoms needing new substituents
 -k            during addition, break bond btw first two matched atoms, discard the fragment
               containing the second atom, and then look for replacements.
               For example to replace an aniline try    -k -s 'c-[NH2]'
 -y            during addition, process all atoms with H's and all non ring single bonds (no -q or -s)
 -L <natoms>   during addition, max atoms lost from parent with -k option
 -r <bonds>    during addition, min radius to be considered
 -a            during addition, only produce molecules at the largest radius found
 -b S,D,T      during addition, only add via the kinds of bond(s) specified
 -C <natoms>   during addition, discard any products with more than <natoms> atoms
 -H ...        during addition, queries that substituents must     contain, SMARTS:ccc for smarts. See -Y allHmatch.
 -K ...        during addition, queries that substituents must NOT contain, SMARTS:ccc for smarts.
 -p            during addition, write starting molecule
 -u <number>   during addition, min number of example molecules needed for a fragment to be added
 -V <fname>    during addition, avoid forming any of the molecules in <fname>
 -I            during addition, remove isotopes before writing
 -Y ...        more options, enter '-Y help' for info
 -l            reduce to largest fragment (automatically turned on)
 -c            remove all chirality (automatically turned on during building)
 -i <type>     input specification
 -g ...        chemical standardisation options
 -E ...        standard element specifications
 -A ...        standard aromaticity specifications
 -v            verbose output

Suggested usage for Building a database:
identify_substituents -d /path/to/db.bdb -B -R 5 -w 10 -M 12 -v -Y dbproto collection.smi
And then for lookups:
identify_substituents -d /path/to/db.bdb  -M 4 -s '[cD3x2]!@*' -k -r 2  file.smi > new.smi
which looks up max 4 atoms substituents, that are identical to radius 2 contexts.
)";
  // clang-format on

  exit(rc);
}

void
SubstituentIdentification::_preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  m.revert_all_directional_bonds_to_non_directional();

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  if (_strip_leading_zeros) {
    const_IWSubstring mname = m.name();
    if (mname.starts_with('0')) {
      mname.remove_leading_chars('0');
      m.set_name(mname);
    }
  }

  return;
}

int
SubstituentIdentification::_ok_atom_count(const int matoms) const {
  // cerr << "_ok_atom_count checking " << matoms << " vs " <<
  // _min_substituent_size << " to " << _max_substituent_size << '\n';

  if (matoms > _max_substituent_size) {
    return 0;
  }

  if (matoms < _min_substituent_size) {
    return 0;
  }

  return 1;
}

// Return true of `smiles` satisfies both _min_substituent_size and
// _max_substituent_size.
// `natoms` is determined if it is zero or less.
int
SubstituentIdentification::OkSubstituentSize(const std::string& replacement_smiles,
                        int& atoms_in_replacement,
                        int atoms_in_fragment_lost) const {

  if (_max_substituent_size > 0) {
    if (atoms_in_replacement <= 0) {
      const_IWSubstring s(replacement_smiles);
      atoms_in_replacement = lillymol::count_atoms_in_smiles(s);
    }
    if (atoms_in_replacement > _max_substituent_size) {
      return 0;
    }
  }

  if (_min_substituent_size > 0) {
    if (atoms_in_replacement <= 0) {
      const_IWSubstring s(replacement_smiles);
      atoms_in_replacement = lillymol::count_atoms_in_smiles(s);
    }
    if (atoms_in_replacement < _min_substituent_size) {
      return 0;
    }
  }

  if (_max_extra_atoms_added > 0 && atoms_in_fragment_lost > 0) {
    if (atoms_in_replacement <= 0) {
      atoms_in_replacement = lillymol::count_atoms_in_smiles(replacement_smiles);
    }

    if (atoms_in_replacement <= atoms_in_fragment_lost) {
        // replacement is smaller, no need to check.
    } else if ((atoms_in_replacement - atoms_in_fragment_lost) > static_cast<int>(_max_extra_atoms_added)) {
      return 0;
    }
  }

  if (_max_atoms_lost_during_addition > 0 && atoms_in_fragment_lost > 0) {
    if (atoms_in_replacement <= 0) {
      atoms_in_replacement = lillymol::count_atoms_in_smiles(replacement_smiles);
    }

    if (atoms_in_replacement >= atoms_in_fragment_lost) {
      // no need to check
    } else if ((atoms_in_fragment_lost - atoms_in_replacement) > static_cast<int>(_max_atoms_lost_during_addition)) {
      return 0;
    }
  }

  return 1;
}

struct DBKey {
  int radius;
  uint32_t b;
};

int
SubstituentIdentification::_write_in_memory_hashes_to_database() {
  for (auto i = 0; i <= _shell_radius; ++i) {
    if (!_write_in_memory_hash_to_database(i, _bit[i])) {
      return 0;
    }
  }

  return 1;
}

int
SubstituentIdentification::_write_in_memory_hash_to_database(
    int radius, const std::unordered_map<uint32_t, Usmi2Substituent>& b) {
  struct DBKey dbkey {
    radius, 0
  };

  for (auto i = b.cbegin(); i != b.cend(); ++i) {
    dbkey.b = (*i).first;

    if (!_write_substituents_associated_with_bit(dbkey, (*i).second)) {
      return 0;
    }
  }

  return 1;
}

int
SubstituentIdentification::_write_substituents_associated_with_bit(
    const struct DBKey& dbkey, const Usmi2Substituent& u) {
  if (_store_protos) {
    return StoreProto(dbkey, u);
  }

  return 1;
}

// Pointers to keys and values for sorting.
struct SmiSubstituent{
  const IWString* smi;
  const ASubstituent* substituent;
};

struct SmiSubstituentComparitor {
  int operator() (const SmiSubstituent& sc1, const SmiSubstituent& sc2) {
    uint32_t n1 = sc1.substituent->number_instances();
    uint32_t n2 = sc2.substituent->number_instances();

    if (n1 > n2) {
      return -1;
    }
    if (n1 == n2) {
      return 0;
    }

    return 1;
  }
};

int
SubstituentIdentification::WriteTextProto(const struct DBKey& dbkey,
               const substituent_identification::Replacements& proto,
               IWString_and_File_Descriptor & output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  static constexpr char kSep = ' ';
  output << dbkey.radius << kSep << dbkey.b << kSep;

  std::string buffer;
  if (! printer.PrintToString(proto, &buffer)) {
    cerr << "WriteTextProto:cannot write '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  output << buffer;
  output << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
SubstituentIdentification::StoreProto(
    const struct DBKey& dbkey, const Usmi2Substituent& u) {
  substituent_identification::Replacements proto;

  const uint32_t n = u.size();
  std::unique_ptr<SmiSubstituent[]> tmp = std::make_unique<SmiSubstituent[]>(n);
  uint32_t ndx = 0;
  for (const auto& [usmi, substituent] : u) {
    tmp[ndx].smi = &usmi;
    tmp[ndx].substituent = &substituent;
    ++ndx;
  }
  // Decision about how to sort. Should it be by number of atoms or prevalence.
  // Use prevalence since it is already computed.
  SmiSubstituentComparitor cmp;

  iwqsort(tmp.get(), n, cmp);

  for (uint32_t i = 0; i < n; ++i) {
    substituent_identification::Replacement* r = proto.mutable_replacement()->Add();
    const SmiSubstituent& ss = tmp[i];
    r->set_smiles(ss.smi->AsString());
    r->add_id(ss.substituent->first_molecule().AsString());
    r->set_n(ss.substituent->number_instances());

    _pairs_stored += ss.substituent->number_instances();
  }

  std::string serialized;
  proto.SerializeToString(&serialized);

  Dbt zkey{(void*)(&dbkey), sizeof(dbkey)};
  Dbt zdata(serialized.data(), serialized.size());

  if (const int s = _dbs[0]->put(NULL, &zkey, &zdata, 0); s != 0) {
    cerr << "Berkeley database Put operation failed\n";
    _dbs[0]->err(s, "");
    return 0;
  }

  ++_keys_stored;

  if (_stream_for_textproto.active()) {
    WriteTextProto(dbkey, proto, _stream_for_textproto);
  }

  return 1;
}

void
SubstituentIdentification::_initialise_msta(Molecule& m, int* number,
                                            Molecule_Specific_Temporary_Arrays& msta) {
  _atom_typing_specification.assign_atom_types(m, msta._atype);

  msta.initialise_aromatic_bond_array(m);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    number[i] = i;
    m.set_user_specified_atom_void_ptr(i, number + i);
  }

  return;
}

/*
  We need to remove the fragment containing a0, but we need to capture
  the unique smiles of the fragment being removed
*/

int
SubstituentIdentification::_identify_fragment_being_removed(
    Molecule& m, const atom_number_t attachment_point, const int* to_remove,
    FragmentLost& fragment_lost) const {
  if (_apply_atom_map_labels) {
    m.set_atom_map_number(attachment_point, 1);
  } else {
    m.set_isotope(attachment_point, 1);
  }

  Molecule being_discarded;
  m.create_subset(being_discarded, to_remove, 1);

  fragment_lost.smiles = being_discarded.unique_smiles();
  fragment_lost.natoms = being_discarded.natoms();

  return 1;
}

/*
  This can be made a lot more efficient
  Break the bond between `a0` and `a1`.

  Returns the atom number of the `a0` - which may change because of
  the removal of the fragment containing `a1`.
*/

atom_number_t
SubstituentIdentification::_do_break_across_first_two_matched_atoms(
    Molecule& m, atom_number_t a0, atom_number_t a1,
    Molecule_Specific_Temporary_Arrays& msta, FragmentLost& fragment_lost) {
  if (!m.are_bonded(a0, a1)) {
    cerr << "SubstituentIdentification::_do_break_across_first_two_matched_"
            "atoms:atoms not bonded\n";
    return kInvalidAtomNumber;
  }

#ifdef DEBUG_DO_BREAK_ACROSS_FIRST_TWO_MATCHED_ATOMS
  Molecule mcopy(m);
  for (int i = 0; i < m.natoms(); ++i) {
    mcopy.set_atom_map_number(i, i);
  }
  cerr << mcopy.smiles() << " atoms " << a0 << " and " << a1 << '\n';
#endif

  int* to_remove = msta.to_remove();

  const auto atoms_being_removed = m.identify_side_of_bond(to_remove, a1, 1, a0);

  if (1 == m.ncon(a1)) {  //  a common case
    if (a1 < a0) {  // Adjust for loss of a1
      --a0;
    }
    if (_write_fragments_added) {  
      fragment_lost.smiles = m.atomic_symbol(a1);
      fragment_lost.natoms = 1;
    }
    m.remove_atom(a1);
    return a0;
  }

  // cerr << "Found " << atoms_being_removed << " atoms to be removed\n";

  if (0 == atoms_being_removed) [[unlikely]] {  // should not happen
    return kInvalidAtomNumber;
  }

  if (atoms_being_removed > _max_atoms_lost_during_build) {
    return kInvalidAtomNumber;
  }

  // We are going to remove all the atoms connected to a1, but we then need to
  // work out the new atom number of atom a0 which will be retained in the
  // resulting molecule.

  const int delta =
      std::count_if(to_remove, to_remove + a0, [](const int x) { return 1 == x; });

  // Do we need to capture the unique smiles of what is discarded?
  if (_write_fragments_added) {  
    _identify_fragment_being_removed(m, a1, to_remove, fragment_lost);
  }

  m.remove_atoms(to_remove);

  return a0 - delta;
}

static int
symmetrically_equivalent_to_something_already_found(const Bond& b,
                                                    const int* process_these_pairs,
                                                    const int ndx,
                                                    const int* symmetry_classes) {
  const int sa1 = symmetry_classes[b.a1()];
  const int sa2 = symmetry_classes[b.a2()];

  for (int i = 0; i < ndx; i += 2) {
    const int si1 = symmetry_classes[process_these_pairs[i]];
    const int si2 = symmetry_classes[process_these_pairs[i + 1]];

    if (si1 == sa1 && si2 == sa2) {
      return 1;
    }
    if (si1 == sa2 && si2 == sa1) {
      return 1;
    }
  }

  return 0;
}

int
SubstituentIdentification::_identify_non_ring_single_bonds(Molecule& m,
                int* process_these_pairs) const {
  (void)m.ring_membership();

  const int* symmetry_classes = m.symmetry_classes();

  int ndx = 0;

  const int nb = m.nedges();

  for (int i = 0; i < nb; ++i) {
    const Bond* b = m.bondi(i);

    if (!b->is_single_bond()) {
      continue;
    }

    if (b->nrings()) {
      continue;
    }

    if (symmetrically_equivalent_to_something_already_found(*b, process_these_pairs, ndx,
                                                            symmetry_classes)) {
      continue;
    }

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

int
SubstituentIdentification::_try_bond_breaking_across(
    const Molecule& m, const atom_number_t a0, const atom_number_t a1,
    const int* to_remove, IW_STL_Hash_Set& smiles_already_found,
    Molecule_Specific_Temporary_Arrays& msta, IWString_and_File_Descriptor& output) {
  FragmentLost fragment_lost;

  Molecule mcopy(m);

  const int delta =
      std::count_if(to_remove, to_remove + a1, [](const int x) { return 1 == x; });

  // If writing fragments added, capture the unique smiles of what is discarded.
  if (_write_fragments_added) {
    _identify_fragment_being_removed(mcopy, a0, to_remove, fragment_lost);
  }

  mcopy.remove_atoms(to_remove);

  return _look_for_new_substituents(mcopy, a1 - delta, fragment_lost,
                                    smiles_already_found, msta, output);
}

/*
  We try both ways across the bond
*/

int
SubstituentIdentification::_try_bond_breaking_across(
    const Molecule& m, const atom_number_t a0, const atom_number_t a1,
    IW_STL_Hash_Set& smiles_already_found, Molecule_Specific_Temporary_Arrays& msta,
    IWString_and_File_Descriptor& output) {
  if (!m.are_bonded(a0, a1)) {
    cerr << "SubstituentIdentification::_try_bond_breaking_across:atoms not "
            "bonded '"
         << m.name() << "' atoms " << a0 << " " << a1 << '\n';
    return 0;
  }

  const int matoms = m.natoms();

  int* to_remove = msta.to_remove();  // comes out zero filled

  const auto atoms_being_removed = m.identify_side_of_bond(to_remove, a0, 1, a1);

  if (0 == atoms_being_removed) [[unlikely]] {  // should not happen
    return 0;
  }

  if (atoms_being_removed <= _max_atoms_lost_during_build) {
    _try_bond_breaking_across(
        m, a0, a1, to_remove, smiles_already_found, msta,
        output);  // kind of dangerous to pass both to_remove and msta...
  }

  // Now try things the other way around

  if (matoms - atoms_being_removed <= _max_atoms_lost_during_build) {
    std::transform(to_remove, to_remove + matoms, to_remove,
                   [](const int x) { return !x; });

    _try_bond_breaking_across(m, a1, a0, to_remove, smiles_already_found, msta, output);
  }

  return 1;
}

/*
  When doing matched pairs, we have a lot of situations.
  There may be multiple sites of matched pairs.
  There may be different numbers of examples at each of those sites.
  The shifts will be different within a site
*/

int
SubstituentIdentification::_matched_pairs_qsar(Molecule& m,
                                               IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  // cerr << m.smiles() << " _matched_pairs_qsar \n";

  // cerr << "Looking for new substituents to " << m.smiles() << '\n';

  std::unique_ptr<int[]> numbers = std::make_unique<int[]>(matoms);

  Molecule_Specific_Temporary_Arrays msta(matoms);

  _initialise_msta(m, numbers.get(), msta);

  int queries_hitting_this_molecule = 0;

  IW_STL_Hash_Set smiles_already_found;

  smiles_already_found.insert(m.unique_smiles());  // do not want to recreate the parent

  MMP_Related mmp;

  // If we are losing an H atom, keep one copy of the smiles
  FragmentLost hydrogen;
  hydrogen.SetHydrogen();

  // non ring single bonds and atoms with Hydrogen
  if (_default_new_molecule_starting_points) {

    for (int i = 0; i < matoms; ++i) {
      if (0 == m.hcount(i)) {
        continue;
      }

      _matched_pairs_qsar_by_radius(m, i, hydrogen, smiles_already_found, msta, mmp, output);
      queries_hitting_this_molecule++;
    }

    std::unique_ptr<int[]> process_these_pairs = std::make_unique<int[]>(matoms + matoms);

    const auto nb = _identify_non_ring_single_bonds(m, process_these_pairs.get());
    if (nb > 0) {
      queries_hitting_this_molecule++;

      for (int i = 0; i < nb; ++i) {
        atom_number_t a1 = process_these_pairs[i + i];
        atom_number_t a2 = process_these_pairs[i + i + 1];

        _matched_pairs_qsar_across_bond(m, a1, a2, smiles_already_found, msta, mmp, output);
      }
    }
  } else {  // based on queries
    Molecule_to_Match target(&m);

    for (Substructure_Query* q : _anchor_query) {
      Substructure_Results sresults;

      if (! q->substructure_search(target, sresults)) {
        continue;
      }

      queries_hitting_this_molecule++;

      if (_break_molecule_at_first_two_matched_atoms) {
        for (const Set_of_Atoms* e : sresults.embeddings()) {
          Molecule mcopy(m);

          if (e->size() < 2) [[unlikely]] {  // should never happen
            continue;
          }

          _matched_pairs_qsar_across_bond(mcopy, e->item(0), e->item(1),
                                          smiles_already_found, msta, mmp, output);
        }
      } else {
        for (const Set_of_Atoms* e : sresults.embeddings()) {
          const atom_number_t k = e->front();

          // will not happen if the person constructed the query properly
          if (0 == m.hcount(k)) [[unlikely]] {
            continue;
          }

          if (!_matched_pairs_qsar_by_radius(m, k, hydrogen, smiles_already_found,
                                   msta, mmp, output)) {
            continue;
          }
        }
      }
    }
  }

  // cerr << "_matched_pairs_qsar complete\n";
  if (queries_hitting_this_molecule) {
    _molecules_hitting_queries++;
  }

  return 1;
}

int
SubstituentIdentification::_get_matched_atoms_to_process(const Set_of_Atoms* e) {
  const int n = e->number_elements();

  if (n > _matched_atoms_to_process) {
    return _matched_atoms_to_process;
  }

  return n;
}

// Return true if we have already made more then _max_molecules_per_input_molecule molecules
// from the current parent.
int
SubstituentIdentification::AlreadyMadeEnoughMolecules(const UsedDuringLookups& lookup_data) const {
  if (lookup_data.molecules_formed_from_current_molecule() > _max_molecules_per_input_molecule) {
    return 1;
  }

  return 0;
}

int
SubstituentIdentification::_look_for_new_substituents(Molecule& m, 
                    UsedDuringLookups& lookup_data,
                    IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  if (matoms >= _max_atoms_in_product) {
    cerr << m.name() << " already has " << matoms << " atoms, but max in products is "
         << _max_atoms_in_product << " skipped\n";
    return 1;
  }

  // cerr << "Looking for new substituents to " << m.smiles() << '\n';

  std::unique_ptr<int[]> numbers = std::make_unique<int[]>(matoms);

  Molecule_Specific_Temporary_Arrays msta(matoms);

  _initialise_msta(m, numbers.get(), msta);

  int queries_hitting_this_molecule = 0;

  int parent_written = 0;

  // Do not want to recreate the parent.
  lookup_data.ProductMoleculeSeenBefore(m);

  FragmentLost fragment_lost;

  // Note that we do not check _max_molecules_per_input_molecule as often
  // as we could, for efficiency and avoiding clutter in the code.

  // non ring single bonds and atoms with Hydrogen atoms
  if (_default_new_molecule_starting_points) {
    fragment_lost.smiles = "H";
    fragment_lost.natoms = 1;

    for (int i = 0; i < matoms; ++i) {
      if (0 == m.hcount(i)) {
        continue;
      }

      _look_for_new_substituents(m, i, lookup_data, msta, output);
      queries_hitting_this_molecule++;
    }

    if (AlreadyMadeEnoughMolecules(lookup_data)) {
      return 1;
    }

    std::unique_ptr<int[]> process_these_pairs = std::make_unique<int[]>(matoms + matoms);

    const auto nb = _identify_non_ring_single_bonds(m, process_these_pairs.get());
    if (nb > 0) {
      queries_hitting_this_molecule++;

      for (int i = 0; i < nb; ++i) {
        atom_number_t a1 = process_these_pairs[i + i];
        atom_number_t a2 = process_these_pairs[i + i + 1];

        if (_write_parent_molecule && 0 == i) {
          output << m.smiles() << ' ' << m.name() << '\n';
        }

        _try_bond_breaking_across(m, a1, a2, lookup_data.product_already_formed(), msta, output);
        if (AlreadyMadeEnoughMolecules(lookup_data)) {
          break;
        }
      }
    }
  } else {
    Molecule_to_Match target(&m);

    for (Substructure_Query* q : _anchor_query) {
      Substructure_Results sresults;

      if (int nhits = q->substructure_search(target, sresults); nhits == 0) {
        continue;
      }

      queries_hitting_this_molecule++;

      if (_write_parent_molecule && !parent_written) {
        output << m.smiles() << ' ' << m.name() << '\n';
        parent_written = 1;
      }

      if (_break_molecule_at_first_two_matched_atoms) {
        for (const Set_of_Atoms* e : sresults.embeddings()) {
          // cerr << "Matched atoms " << *e << '\n';

          if (e->number_elements() < 2) [[ unlikely ]] {  // should never happen
            continue;
          }

          Molecule mcopy(m);
          const atom_number_t attachment_point = _do_break_across_first_two_matched_atoms(
              mcopy, e->item(0), e->item(1), msta, fragment_lost);
          // cerr << m.name() << ' ' << " atoms " << e->item(0) << " and " << e->item(1) << " attachment_point " << attachment_point << '\n';
          if (kInvalidAtomNumber == attachment_point) {
            continue;
          }

          //        cerr << "From " << m.smiles() << "\nget  " << mcopy.smiles()
          //        << " atom " << attachment_point << " type " <<
          //        mcopy.smarts_equivalent_for_atom(attachment_point) << '\n';

          _look_for_new_substituents(mcopy, attachment_point, lookup_data, msta, output);
          if (AlreadyMadeEnoughMolecules(lookup_data)) {
            break;
          }
        }
      } else {
        lookup_data.smiles_of_fragment_lost() = "H";
        for (const Set_of_Atoms* e : sresults.embeddings()) {
          const int nprocess = _get_matched_atoms_to_process(e);

          for (int x = 0; x < nprocess; ++x) {
            const atom_number_t k = e->item(x);

            if (INVALID_ATOM_NUMBER == k) {
              continue;
            }

            // Will not happen if the query was properly constructed.
            if (0 == m.hcount(k)) {
              continue;
            }

            _look_for_new_substituents(m, k, lookup_data, msta, output);
          }
        }

        if (AlreadyMadeEnoughMolecules(lookup_data)) {
          break;
        }
      }
    }
  }

  if (queries_hitting_this_molecule) {
    _molecules_hitting_queries++;
  }

  if (lookup_data.molecules_formed_from_current_molecule() >
        _max_molecules_per_input_molecule) {
    _molecules_producing_too_many_new_molecules++;
    if (_verbose) {
      cerr << m.smiles() << ' ' << m.name() << " would have produced more than "
           << _max_molecules_per_input_molecule << " new molecules\n";
      cerr << "Increase limit with '-Y maxgen=nnnn' if desired, or use a larger radius (-r) or increase support requirement (-u)\n";
    }
  }

  return 1;
}


int
SubstituentIdentification::_look_for_new_substituents(
    Molecule& m, atom_number_t zatom, FragmentLost& fragment_lost,
    IW_STL_Hash_Set& smiles_already_found, Molecule_Specific_Temporary_Arrays& msta,
    IWString_and_File_Descriptor& output) {
  std::unique_ptr<uint32_t[]> b = std::make_unique<uint32_t[]>(_shell_radius + 1);

  auto max_shell_radius_formed = _compute_environment01(m, zatom, _shell_radius, b.get(), msta);

// #define ECHO_BITS_FORMED
#ifdef ECHO_BITS_FORMED
  cerr << "Looking for additions to " << m.smiles() << " atom " << zatom << ' '
       << m.smarts_equivalent_for_atom(zatom) << " bits";
  for (auto i = 0; i <= max_shell_radius_formed; ++i) {
    cerr << ' ' << b[i];
  }
  cerr << '\n';
#endif

  uint32_t new_molecules_produced = 0;

  cerr << "Radii " << max_shell_radius_formed << " and _min_shell_radius " << _min_shell_radius << '\n';
  for (int r = max_shell_radius_formed; r >= _min_shell_radius; --r) {
    cerr << "Looking for " << b[r] << " at radius " << r << '\n';
    DBKey dbkey{r, b[r]};

    const int tmp = _look_for_new_substituents(m, zatom, fragment_lost,
                                               smiles_already_found, dbkey,
                                               new_molecules_produced, output);

    //  cerr << tmp << " molecules generated at radius " << r << '\n';

    if (0 == tmp) {
      continue;
    }

    new_molecules_produced += tmp;

    if (new_molecules_produced > _max_molecules_per_input_molecule) {
      return new_molecules_produced;
    }

    if (_only_produce_molecules_at_biggest_radius) {
      return new_molecules_produced;
    }
  }

  // cerr << m.name() << " generated " << new_molecules_produced << '\n';

  return new_molecules_produced;
}

int
SubstituentIdentification::_look_for_new_substituents(
    Molecule& m, atom_number_t zatom, 
    UsedDuringLookups& lookup_data,
    Molecule_Specific_Temporary_Arrays& msta,
    IWString_and_File_Descriptor& output) {
  std::unique_ptr<uint32_t[]> b = std::make_unique<uint32_t[]>(_shell_radius + 1);

  int max_shell_radius_formed = _compute_environment01(m, zatom, _shell_radius, b.get(), msta);

// #define ECHO_BITS_FORMED
#ifdef ECHO_BITS_FORMED
  cerr << "Looking for additions to " << m.smiles() << " atom " << zatom << ' '
       << m.smarts_equivalent_for_atom(zatom) << " bits";
  for (auto i = 0; i <= max_shell_radius_formed; ++i) {
    cerr << ' ' << b[i];
  }
  cerr << '\n';
#endif

  // cerr << "Radii " << max_shell_radius_formed << " and _min_shell_radius " << _min_shell_radius << '\n';
  for (int r = max_shell_radius_formed; r >= _min_shell_radius; --r) {
    //cerr << "Looking at radius " << r << '\n';
    DBKey dbkey{r, b[r]};

    if (! _look_for_new_substituents(m, zatom, lookup_data, dbkey, output)) {
      continue;
    }

    if (AlreadyMadeEnoughMolecules(lookup_data)) {
      return 1;
    }

    if (_only_produce_molecules_at_biggest_radius) {
      return 1;
    }
  }

  // cerr << m.name() << " generated " << new_molecules_produced << '\n';

  return 1;
}


int
SubstituentIdentification::_look_for_new_substituents(
    Molecule& m, atom_number_t zatom, FragmentLost& fragment_lost,
    IW_STL_Hash_Set& smiles_already_found,
    struct DBKey& dbkey,
    uint32_t& new_molecules_produced,
    IWString_and_File_Descriptor& output) {
  Dbt zkey{(void*)(&dbkey), sizeof(dbkey)};

  // cerr << "Looking for bit " << dbkey.b << " rad " << dbkey.radius << " lost " << fragment_lost.smiles << '\n';

  IW_STL_Hash_Set already_processed;  // names of reagents when dealing with multiple
                                     // databases, don't re-create molecules

  int rc = 0;
  for (int i = 0; i < _ndb; ++i) {
    int tmp = _look_for_new_substituents_db(m, zatom, fragment_lost,
                                            smiles_already_found, dbkey, *(_dbs[i]), zkey,
                                            already_processed, output);

    if (tmp == 0) {
      continue;
    }

    new_molecules_produced += tmp;
    rc += tmp;

    if (new_molecules_produced > _max_molecules_per_input_molecule) {
      return rc;
    }
  }

  return rc;
}

int
SubstituentIdentification::_look_for_new_substituents(
    Molecule& m, atom_number_t zatom, 
    UsedDuringLookups& lookup_data,
    struct DBKey& dbkey,
    IWString_and_File_Descriptor& output) {
  Dbt zkey{(void*)(&dbkey), sizeof(dbkey)};

  // cerr << "Looking for bit " << dbkey.b << " rad " << dbkey.radius << " ndb " << _ndb << '\n';

  // If there are multiple databases, we don't want to add the same reagent
  // multiple times.
  lookup_data.ClearReagentNames();

  // cerr << "LINE " << __LINE__ << '\n';

  int rc = 0;
  for (int i = 0; i < _ndb; ++i) {
    int tmp = _look_for_new_substituents_db(m, zatom, lookup_data, dbkey,
                                            *(_dbs[i]), zkey, output);
    if (tmp == 0) {
      continue;
    }

    rc += tmp;

    if (lookup_data.molecules_formed_from_current_molecule() > _max_molecules_per_input_molecule) {
      return rc;
    }
  }

  return rc;
}

int
SubstituentIdentification::_look_for_new_substituents_db(
    Molecule& m, atom_number_t zatom, FragmentLost& fragment_lost,
    IW_STL_Hash_Set& smiles_already_found, const struct DBKey& rad_and_bit, Db& db,
    Dbt& dbkey, IW_STL_Hash_Set& already_processed, IWString_and_File_Descriptor& output) {
  Dbt zdata;

  cerr << "Looking up bit q " << rad_and_bit.b << '\n';
  if (int s = db.get(NULL, &dbkey, &zdata, 0); s != 0) {
    return 0;
  }

  // cerr << "Found bit in database\n";

  const_IWSubstring fromdb(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());

  if (_store_protos) {
    const absl::string_view sv(fromdb.data(), fromdb.length());
    substituent_identification::Replacements proto;
    proto.ParseFromString(sv);

    return FormNewMolecules(m, zatom, fragment_lost, smiles_already_found,
                        rad_and_bit, proto, already_processed, output);
  }

  return 0;
}

// Parse the contents of `fromdb` as a proto.
static int
DbToProto(const Dbt& fromdb, substituent_identification::Replacements& proto) {
  const absl::string_view sv(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());

  if (! proto.ParseFromString(sv)) {
    cerr << "DbToProto::_matched_pairs_qsar_db:invalid proto data\n";
    return 0;
  }

  return 1;
}

// Proto only implementation.
int
SubstituentIdentification::_look_for_new_substituents_db(
    Molecule& m, atom_number_t zatom,
    UsedDuringLookups& lookup_data, 
    const struct DBKey& rad_and_bit, Db& db,
    Dbt& dbkey, IWString_and_File_Descriptor& output) {
  assert(_store_protos);

  Dbt fromdb;
  if (db.get(NULL, &dbkey, &fromdb, 0) != 0) {
    return 0;
  }

  substituent_identification::Replacements proto;
  if (! DbToProto(fromdb, proto)) {
    return 0;
  }

  return FormNewMolecules(m, zatom, rad_and_bit, lookup_data, proto, output);
}

#ifdef OLD_VERSION_WITHOUT_PROTOS
/*
  A typical database record will look like
  [1OH2]:36473:129|[1OH2]=:36473:1
*/

int
SubstituentIdentification::_form_new_molecules(Molecule& m, atom_number_t zatom,
                                               FragmentLost& fragment_lost,
                                               IW_STL_Hash_Set& smiles_already_found,
                                               const struct DBKey& rad_and_bit,
                                               const_IWSubstring& fromdb,
                                               IW_STL_Hash_Set& already_processed,
                                               IWString_and_File_Descriptor& output) {
  const_IWSubstring token;
  int i = 0;

  uint32_t rc = 0;

  while (fromdb.nextword(token, i, '|')) {
    rc += _form_new_molecule(m, zatom, fragment_lost, smiles_already_found,
                       rad_and_bit, token, already_processed, output);
    if (rc > _max_molecules_per_input_molecule) {
      return rc;
    }
  }

  return rc;
}
#endif

int
SubstituentIdentification::FormNewMolecules(Molecule& m, atom_number_t zatom,
                                               FragmentLost& fragment_lost,
                                               IW_STL_Hash_Set& smiles_already_found,
                                               const struct DBKey& rad_and_bit,
                                               const substituent_identification::Replacements& proto,
                                               IW_STL_Hash_Set& already_processed,
                                               IWString_and_File_Descriptor& output) {
  fragment_lost.EnsureNatomsComputed();

  uint32_t rc = 0;
  for (const substituent_identification::Replacement& replacement: proto.replacement()) {
    if (replacement.n() < _min_examples_needed_for_addition) {
      continue;
    }
    int natoms = 0;
    if (! OkSubstituentSize(replacement.smiles(), natoms, fragment_lost.natoms)) {
      continue;
    }
    // Check for size of product here? Complicated by possible loss of atoms...
    // TODO:ianwatson figure this out
    rc += FormNewMoleculeAtAtomWithReplacement(m, zatom, fragment_lost, smiles_already_found,
                rad_and_bit, replacement, already_processed, output);
    if (rc > _max_molecules_per_input_molecule) {
      return rc;
    }
  }

  return rc;
}

int
SubstituentIdentification::FormNewMolecules(Molecule& m, atom_number_t zatom,
                                            const struct DBKey& rad_and_bit,
                                            UsedDuringLookups& lookup_data,
                                            const substituent_identification::Replacements& proto,
                                            IWString_and_File_Descriptor& output) {
  static constexpr int kZero = 0;

  int rc = 0;
  for (const substituent_identification::Replacement& replacement: proto.replacement()) {
    if (replacement.n() < _min_examples_needed_for_addition) {
      continue;
    }
    int natoms = 0;
    if (! OkSubstituentSize(replacement.smiles(), natoms, kZero)) {
      continue;
    }
    // Check for size of product here? Complicated by possible loss of atoms...
    // TODO:ianwatson figure this out, I think it is OK to check here...
    rc += FormNewMoleculeAtAtom(m, zatom, rad_and_bit, lookup_data, replacement, output);
    if (lookup_data.molecules_formed_from_current_molecule() > _max_molecules_per_input_molecule) {
      return rc;
    }
  }

  return rc;
}


int
SubstituentIdentification::AddFragment(Molecule& m, atom_number_t zatom,
                                       Molecule& frag, atom_number_t fragment_attachment_point,
                                       bond_type_t bt) const {
  const int initial_matoms = m.natoms();

  m.add_molecule(&frag);

  m.add_bond(zatom, initial_matoms + fragment_attachment_point, bt);

  m.unset_all_implicit_hydrogen_information(zatom);
  m.set_implicit_hydrogens_known(zatom, 0);
  m.unset_all_implicit_hydrogen_information(initial_matoms + fragment_attachment_point);
  m.set_implicit_hydrogens_known(initial_matoms + fragment_attachment_point, 0);

  if (_remove_isotopes_from_product) {
    ;
  } else if (_apply_atom_map_labels) {
    m.set_atom_map_number(zatom, 1);
  } else {
    m.set_isotope(zatom, 1);
  }

  m.recompute_implicit_hydrogens(zatom);
  m.recompute_implicit_hydrogens(initial_matoms + fragment_attachment_point);

  return 1;
}

static void
IssueInvalidValenceWarning(Molecule& m, Molecule& frag, std::ostream& output) {
  output << "Warning, invalid valence '" << m.name() << " with " << frag.smiles() << ' '
           << m.smiles() << "\n";
  for (int i = 0; i < m.natoms(); ++i) {
    if (m.valence_ok(i)) {
      continue;
    }

    output << " problem with atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << "\n";
  }
}

int
SubstituentIdentification::FormNewMoleculeAtAtomWithReplacement(Molecule& m, atom_number_t zatom,
                                               FragmentLost& fragment_lost,
                                               IW_STL_Hash_Set& smiles_already_found,
                                               const struct DBKey& rad_and_bit,
                                               const substituent_identification::Replacement& proto,
                                               IW_STL_Hash_Set& already_processed,
                                               IWString_and_File_Descriptor& output) {
  const_IWSubstring smiles(proto.smiles());
  bond_type_t bt;
  if (! MaybeRemoveTrailingBondOkH(m, zatom, smiles, bt)) {
    return 0;
  }

  if (_ndb > 0) {
    if (already_processed.contains(smiles)) {
      return 0;
    }

    already_processed.insert(smiles);
  }

  // pulling off something and putting it back on.
  // TODO:ianwatson does this actually work?
  if (smiles == fragment_lost.smiles) {
    return 0;
  }

  // string comparisons not fully implemented.
  if (iwstring::Equals(m.name(), proto.id(0))) {  // maybe molecule being merged with itself
    return 0;
  }

  Molecule frag;
  if (!frag.build_from_smiles(smiles)) {
    cerr << "SubstituentIdentification::FormNewMolecule:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  // Atom count in frag does not need to be checked here.

  const int initial_matoms = m.natoms();

  // cerr << "Fragment contains " << frag.natoms() << " atoms\n";
  if (initial_matoms + frag.natoms() > _max_atoms_in_product) {
    return 0;
  }

  const atom_number_t fragment_attachment_point = frag.atom_with_isotope(1);
  if (INVALID_ATOM_NUMBER == fragment_attachment_point) {
    cerr << "SubstituentIdentification::FormNewMolecule:no isotopic atom in "
            "stored fragment '" << smiles << '\n';
    return 0;
  }

  if (! OkSubstructureQueries(frag)) {
    return 0;
  }

  AddFragment(m, zatom, frag, fragment_attachment_point, bt);

  if (!m.valence_ok()) {
    if (_verbose > 1) {
      IssueInvalidValenceWarning(m, frag, cerr);
    }
    _invalid_valences_ignored++;
  } else if (_check_already_made(m, smiles_already_found)) {
    ;
  } else {
    if (_remove_isotopes_from_product) {
      m.set_isotope(initial_matoms + fragment_attachment_point, 0);
    }

    if (_output_is_textproto) {
      WriteTextProto(m, proto, frag, rad_and_bit, fragment_lost, output);
    } else {
      output << m.smiles() << ' ' << m.name() << " %% ";

      output << proto.id(0) << ":R" << rad_and_bit.radius << " N:" << proto.n();

      if (_write_fragments_added) {
        output << ' ' << fragment_lost.smiles << '%' << smiles;
      }

      output << '\n';
    }

    _molecules_written++;

    output.write_if_buffer_holds_more_than(32768);
  }

  m.resize(initial_matoms);
  m.set_isotope(zatom, 0);

  return 1;
}

int
SubstituentIdentification::WriteTextProto(Molecule& m,
                const substituent_identification::Replacement& proto,
                Molecule& frag,
                const struct DBKey& rad_and_bit,
                FragmentLost& fragment_lost,
                IWString_and_File_Descriptor& output) const {

  substituent_identification::Result result;
  const IWString& s = m.smiles();
  result.set_smi(s.data(), s.length());
  const IWString& id = m.name();
  result.set_id(id.data(), id.length());
  result.set_donor(proto.id(0));
  result.set_rad(rad_and_bit.radius);
  result.set_n(proto.n());
  if (_write_fragments_added) {
    result.set_frag(fragment_lost.smiles.data(), fragment_lost.smiles.length());
  }

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (! printer.PrintToString(result, &buffer)) {
    cerr << "SubstituentIdentification::WriteTextProto:cannot write '"
         << result.ShortDebugString() << "'\n";
    return 0;
  }

  if (_prepend_smiles_to_textproto) {
    output << m.smiles() << ' ';
  }

  output << buffer << '\n';

  return 1;
}

int
SubstituentIdentification::FormNewMoleculeAtAtom(Molecule& m, atom_number_t zatom,
                                               const struct DBKey& rad_and_bit,
                                               UsedDuringLookups& lookup_data,
                                               const substituent_identification::Replacement& proto,
                                               IWString_and_File_Descriptor& output) {
  const_IWSubstring smiles(proto.smiles());
  bond_type_t bt;
  if (! MaybeRemoveTrailingBondOkH(m, zatom, smiles, bt)) {
    return 0;
  }

  if (_ndb > 0) {
    if (lookup_data.ReagentAlreadyUsedThisMolecule(proto.smiles())) {
      return 0;
    }
  }

  // pulling off something and putting it back on.
  // TODO:ianwatson does this actually work?
  if (smiles == lookup_data.smiles_of_fragment_lost()) {
    return 0;
  }

  // string comparisons not fully implemented.
  if (iwstring::Equals(m.name(), proto.id(0))) {  // maybe molecule being merged with itself
    return 0;
  }

  Molecule frag;
  if (!frag.build_from_smiles(smiles)) {
    cerr << "SubstituentIdentification::FormNewMolecule:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  // Atom count in frag does not need to be checked here.

  const int initial_matoms = m.natoms();

  // cerr << "Fragment contains " << frag.natoms() << " atoms\n";
  if (initial_matoms + frag.natoms() > _max_atoms_in_product) {
    return 0;
  }

  const atom_number_t fragment_attachment_point = frag.atom_with_isotope(1);
  if (INVALID_ATOM_NUMBER == fragment_attachment_point) {
    cerr << "SubstituentIdentification::FormNewMolecule:no isotopic atom in "
            "stored fragment '" << smiles << '\n';
    return 0;
  }

  // Support level for the fragment does not need to be checked.

  if (! OkSubstructureQueries(frag)) {
    return 0;
  }

  AddFragment(m, zatom, frag, fragment_attachment_point, bt);

  if (MoleculeIsRejected(m, frag, lookup_data)) {
    m.resize(initial_matoms);
    m.set_isotope(zatom, 0);
    return 0;
  }

  if (_remove_isotopes_from_product) {
    m.set_isotope(initial_matoms + fragment_attachment_point, 0);
  }

  if (_output_is_textproto) {
    FragmentLost fragment_lost;
    fragment_lost.smiles = lookup_data.smiles_of_fragment_lost();
    WriteTextProto(m, proto, frag, rad_and_bit, fragment_lost, output);
  } else {
    output << m.smiles() << ' ' << m.name() << " %% ";

    output << proto.id(0) << ":R" << rad_and_bit.radius << " N:" << proto.n();

    if (_write_fragments_added) {
      //    output << ' ' << smiles_of_fragment_lost << '%' <<
      //    rad_and_bit.radius << '%' << smiles;
      output << ' ' << lookup_data.smiles_of_fragment_lost() << '%' << smiles;
    }

    output << '\n';
  }

  _molecules_written++;

  lookup_data.AnotherMoleculeFormed();

  output.write_if_buffer_holds_more_than(32768);

  m.resize(initial_matoms);
  m.set_isotope(zatom, 0);

  return 1;
}

int
SubstituentIdentification::MoleculeIsRejected(Molecule& m,
                        Molecule& frag,
                        UsedDuringLookups& lookup_data) {
  if (!m.valence_ok()) {
    if (_verbose > 1) {
      IssueInvalidValenceWarning(m, frag, cerr);
    }
    _invalid_valences_ignored++;
    return 1;
  } 

  if (lookup_data.ProductMoleculeSeenBefore(m)) {
    return 1;
  }

  return 0;  // not rejected
}

/*
  If all example molecules are stored in the database, a typical entry will look
  like

  |O=C(Nc1s[1cH]c([n]1)C)C:1158204.1164662:2|
*/

/*
  We are branching from a single atom.
*/

int
SubstituentIdentification::_matched_pairs_qsar_by_radius(Molecule& m, const atom_number_t zatom,
                                               FragmentLost& fragment_lost,
                                               IW_STL_Hash_Set& smiles_already_found,
                                               Molecule_Specific_Temporary_Arrays& msta,
                                               MMP_Related& mmp,
                                               IWString_and_File_Descriptor& output) {
  std::unique_ptr<uint32_t[]> b = std::make_unique<uint32_t[]>(_shell_radius + 1);

  auto max_shell_radius_formed = _compute_environment01(m, zatom, _shell_radius, b.get(), msta);

// #define ECHO_BITS_FORMED
#ifdef ECHO_BITS_FORMED
  cerr << "Looking for additions to " << m.smiles() << " atom " << zatom << ' '
       << m.smarts_equivalent_for_atom(zatom) << " bits";
  for (auto i = 0; i <= max_shell_radius_formed; ++i) {
    cerr << ' ' << b[i];
  }
  cerr << '\n';
#endif

  uint32_t new_molecules_produced = 0;

  for (auto r = max_shell_radius_formed; r >= _min_shell_radius; --r) {
    DBKey dbkey{r, b[r]};
    Dbt zkey{(void*)(&dbkey), sizeof(dbkey)};

    const auto tmp = _matched_pairs_qsar_db(m, zatom, fragment_lost,
                                            smiles_already_found, zkey, r, mmp, output);

    //  cerr << tmp << " molecules generated at radius " << r << '\n';

    if (0 == tmp) {
      continue;
    }

    new_molecules_produced += tmp;

    if (new_molecules_produced > _max_molecules_per_input_molecule) {
      return new_molecules_produced;
    }

    if (_only_produce_molecules_at_biggest_radius) {
      return new_molecules_produced;
    }
  }

  // cerr << m.name() << " generated " << new_molecules_produced << '\n';

  return new_molecules_produced;
}

int
SubstituentIdentification::_matched_pairs_qsar_db(Molecule& m, const atom_number_t zatom,
                                                  FragmentLost& fragment_lost,
                                                  IW_STL_Hash_Set& smiles_already_found,
                                                  struct Dbt& dbkey, int radius, MMP_Related& mmp,
                                                  IWString_and_File_Descriptor& output) {
  Dbt fromdb;

  fragment_lost.EnsureNatomsComputed();

  // cerr << "_matched_pairs_qsar_db checking " << _ndb << " databases\n";
  for (int i = 0; i < _ndb; ++i) {
    int rc = _dbs[i]->get(NULL, &dbkey, &fromdb, 0);
    if (0 != rc) {
      cerr << "No match in db " << i << '\n';
      continue;
    }
    cerr << "Found bit in database\n";

    substituent_identification::Replacements proto;
    if (! DbToProto(fromdb, proto)) {
      return 0;
    }

    for (const substituent_identification::Replacement& replacement : proto.replacement()) {
      if (replacement.n() < _min_examples_needed_for_addition) {
        continue;
      }

      int natoms = 0;
      if (! OkSubstituentSize(replacement.smiles(), natoms, fragment_lost.natoms)) {
        continue;
      }

      ListMatchedPairs(m, zatom, fragment_lost, replacement,
                      radius, mmp, output);
    }
  }

  return 1;
}

int
SubstituentIdentification::ListMatchedPairs(Molecule& m,
                atom_number_t zatom,
                FragmentLost& fragment_lost,
                const substituent_identification::Replacement& replacement,
                int radius,
                MMP_Related& mmp,
                IWString_and_File_Descriptor& output) {
  if (! OkFragmentSize(replacement.smiles(), fragment_lost)) {
    return 1;
  }

  m.set_isotope(zatom, 1);
  output << m.smiles() << ' ' << m.name() << " radius: " << radius;
  output << " lost: " << fragment_lost.smiles << '\n';

  // The databases I am dealing with only have 1 exemplar stored....
  output << replacement.smiles() << " new sidechain n: " << replacement.n() << " exemplified by";
  for (const std::string& id : replacement.id()) {
    output << ' ' << id << '\n';
  }

  m.set_isotope(zatom, 0);

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
SubstituentIdentification::OkFragmentSize(const std::string& replacement_smiles,
                FragmentLost& fragment_lost) {
  const const_IWSubstring tmp(replacement_smiles);
  int atoms_in_replacement = -1;
  cerr << "OkFragmentSize:checking losing " << fragment_lost.smiles << " replace with " << replacement_smiles << '\n';

  if (_max_substituent_size == 0) {
    // Nothing to check.
  } else  {
    atoms_in_replacement = lillymol::count_atoms_in_smiles(tmp);
    if (atoms_in_replacement > _max_substituent_size) {
      return 0;
    }
  }

  if (_max_extra_atoms_added > 0) {
    if (atoms_in_replacement < 0) {
      atoms_in_replacement = lillymol::count_atoms_in_smiles(replacement_smiles);
    }

    if (atoms_in_replacement <= fragment_lost.natoms) {
      // new fragment is smaller, no need to check extra atoms
    } else if ((atoms_in_replacement - fragment_lost.natoms) > static_cast<int>(_max_extra_atoms_added)) {
      return 0;
    }
  }

  if (_max_atoms_lost_during_addition > 0) {
    if (atoms_in_replacement < 0) {
      atoms_in_replacement = lillymol::count_atoms_in_smiles(replacement_smiles);
    }

    if (atoms_in_replacement >= fragment_lost.natoms) {
      // new fragment is larger, no need to check for atoms lost.
    } else if ((fragment_lost.natoms - atoms_in_replacement) > static_cast<int>(_max_atoms_lost_during_addition)) {
      return 0;
    }
  }

  cerr << "OK with " << fragment_lost.natoms << " lost atoms, gaining " << atoms_in_replacement << " as replacement\n";
  return 1;
}

int
SubstituentIdentification::_matched_pairs_qsar_across_bond(
    Molecule& m, const atom_number_t a1, const atom_number_t a2,
    IW_STL_Hash_Set& smiles_already_found, Molecule_Specific_Temporary_Arrays& msta,
    MMP_Related& mmp, IWString_and_File_Descriptor& output) {
  m.remove_bond_between_atoms(a1, a2);
  resizable_array_p<Molecule> c;
  m.create_components(c);

  const int fma1 = m.fragment_membership(a1);

  int adj = 0;

  for (int i = 0; i < a1; ++i) {
    if (m.fragment_membership(i) != fma1) {
      adj++;
    }
  }

  c[0]->set_name(m.name());
  c[1]->set_name(m.name());

  FragmentLost fragment_lost;
  if (m.fragment_membership(0) == m.fragment_membership(a1)) {
    fragment_lost.Set(*c[1]);
    return _matched_pairs_qsar_by_radius(*c[0], a1 - adj, fragment_lost,
                               smiles_already_found, msta, mmp, output);
  } else {
    fragment_lost.Set(*c[0]);
    return _matched_pairs_qsar_by_radius(*c[1], a1 - adj, fragment_lost,
                               smiles_already_found, msta, mmp, output);
  }

  return 1;
}

int
SubstituentIdentification::OkSubstructureQueries(Molecule& m) {
  if (_replacements_must_contain.empty()) {
    // Do nothing
  } else if (_all_must_have_queries_must_match) {
    if (! lillymol::AllQueriesMatch(m, _replacements_must_contain)) {
      return 0;
    }
  } else if (!lillymol::AnyQueryMatches(m, _replacements_must_contain)) {
    return 0;
  }

  if (_replacements_must_not_contain.size() > 0) {
    if (lillymol::AnyQueryMatches(m, _replacements_must_not_contain)) {
      return 0;
    }
  }

  return 1;
}

int
SubstituentIdentification::_check_already_made(const Molecule& m,
                                               IW_STL_Hash_Set& already_processed) const {
  Molecule mcopy(m);
  mcopy.transform_to_non_isotopic_form();

  const IWString& usmi = mcopy.unique_smiles();

  // cerr << "Checking " << usmi << '\n';

  const auto f = already_processed.find(usmi);

  if (f != already_processed.end()) {
    return 1;
  }

  already_processed.insert(usmi);

  return 0;  // never seen this before
}

int
SubstituentIdentification::MaybeRemoveTrailingBondOkH(Molecule& m,
                           atom_number_t zatom,
                           const_IWSubstring& smiles,
                           bond_type_t& bt) const {
  bt = SINGLE_BOND;

  if (smiles.ends_with('=')) {
    if (m.hcount(zatom) < 2) {
      return 0;
    }

    bt = DOUBLE_BOND;
    smiles.chop();
  } else if (smiles.ends_with('#')) {
    if (m.hcount(zatom) < 3) {
      return 0;
    }

    bt = TRIPLE_BOND;
    smiles.chop();
  } else if (0 == m.hcount(zatom)) {
    return 0;
  }

  if (0 == (_only_add_bond & bt)) {
    return 0;
  }

  return 1;
}

int
SubstituentIdentification::_form_new_molecule(Molecule& m, atom_number_t zatom,
                                              const FragmentLost& fragment_lost,
                                              IW_STL_Hash_Set& smiles_already_found,
                                              const struct DBKey& rad_and_bit,
                                              const const_IWSubstring& buffer,
                                              IW_STL_Hash_Set& already_processed,
                                              IWString_and_File_Descriptor& output) {
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

  bond_type_t bt;
  if (! MaybeRemoveTrailingBondOkH(m, zatom, smiles, bt)) {
    return 0;
  }

  if (_ndb > 0) {
    if (already_processed.contains(smiles)) {
      return 0;
    }

    already_processed.insert(smiles);
  }

  // pulling off something and putting it back on
  if (smiles == fragment_lost.smiles) {
    return 0;
  }

  Molecule f;
  if (!f.build_from_smiles(smiles)) {
    cerr << "SubstituentIdentification::_form_new_molecule:invalid smiles '" << buffer << "'\n";
    return 0;
  }

  if (!_ok_atom_count(f.natoms())) {
    return 0;
  }

  const auto initial_matoms = m.natoms();

  // cerr << "Fragment contains " << f.natoms() << " atoms\n";
  if (initial_matoms + f.natoms() > _max_atoms_in_product) {
    return 0;
  }

  const auto fragment_attachment_point = f.atom_with_isotope(1);
  if (INVALID_ATOM_NUMBER == fragment_attachment_point) {
    cerr << "SubstituentIdentification::_form_new_molecule:no isotopic atom in "
            "stored fragment '"
         << buffer << "'\n";
    return 0;
  }

  if (_min_examples_needed_for_addition > 1 && !_enough_examples(buffer)) {
    return 0;
  }

  if (! OkSubstructureQueries(f)) {
    return 0;
  }

  m.add_molecule(&f);

  m.add_bond(zatom, initial_matoms + fragment_attachment_point, bt);

  m.unset_all_implicit_hydrogen_information(zatom);
  m.set_implicit_hydrogens_known(zatom, 0);
  m.unset_all_implicit_hydrogen_information(initial_matoms + fragment_attachment_point);
  m.set_implicit_hydrogens_known(initial_matoms + fragment_attachment_point, 0);

  if (_remove_isotopes_from_product) {
    ;
  } else if (_apply_atom_map_labels) {
    m.set_atom_map_number(zatom, 1);
  } else {
    m.set_isotope(zatom, 1);
  }

  m.recompute_implicit_hydrogens(zatom);
  m.recompute_implicit_hydrogens(initial_matoms + fragment_attachment_point);

  const_IWSubstring token;

  buffer.nextword(token, i);  // will get everything
  token++;                    // get rid of leading ':'

  if (token == m.name()) {  // molecule being merged with itself
    return 0;
  }

  if (!m.valence_ok()) {
    if (_verbose > 1) {
      cerr << "Warning, invalid valence '" << m.name() << " with '" << token << "' "
           << m.smiles() << "\n";
      for (auto i = 0; i < m.natoms(); ++i) {
        if (m.valence_ok(i)) {
          continue;
        }

        cerr << " problem with atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << "\n";
      }
    }
    _invalid_valences_ignored++;
  } else if (_check_already_made(m, smiles_already_found)) {
    ;
  } else {
    if (_remove_isotopes_from_product) {
      m.set_isotope(initial_matoms + fragment_attachment_point, 0);
    }

    output << m.smiles() << ' ' << m.name() << " %% ";

    output << token << ":R" << rad_and_bit.radius << " B" << rad_and_bit.b;

    if (_write_fragments_added) {
      //    output << ' ' << smiles_of_fragment_lost << '%' <<
      //    rad_and_bit.radius << '%' << smiles;
      output << ' ' << fragment_lost.smiles << '%' << smiles;
    }

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

uint32_t
SubstituentIdentification::_enough_examples(const const_IWSubstring& fromdb) const {
  const auto i = fromdb.rindex(':');

  if (i + 2 == fromdb.length()) {  // most common case, single digit
    return (uint32_t)(fromdb.last_item() - '0') >= _min_examples_needed_for_addition;
  }

  if (i + 3 == fromdb.length()) {  // two digit number of examples
    return (uint32_t)(10 * (fromdb[i + 1] - '0') + (uint32_t)(fromdb[i + 2] - '0')) >=
           _min_examples_needed_for_addition;
  }

  const_IWSubstring s;
  s = fromdb.from_to(i + 1, fromdb.length() - 1);

  uint32_t u;

  if (!s.numeric_value(u) || u < 1) {
    cerr << "SubstituentIdentification::_enough_examples:invalid count stored '" << fromdb
         << "'\n";
    return 0;
  }

  return u >= _min_examples_needed_for_addition;
}

// We have formed a new substituent, and the join point is
// associated with bit `b` and radius `r`.
// Add this information to the hashed data for radius `r`.
void
SubstituentIdentification::_associate_substituent_with_bit(
    int r, uint32_t b, bond_type_t bt, Molecule& substituent,
    const Molecule_Specific_Temporary_Arrays& msta) {
// #define DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
  cerr << "At radius " << r << ", b = " << b << " processing "
       << substituent.unique_smiles() << '\n';
#endif

  // note that we add the attachment type to the end
  IWString usmi(substituent.unique_smiles());

  if (SINGLE_BOND == bt) {
    ;
  } else if (DOUBLE_BOND == bt) {
    usmi << '=';
  } else {
    usmi << '#';
  }

  auto& y = _bit[r];  // appropriate radius hash

  auto f = y.find(b);

  if (f == y.end())  // never seen this bit before
  {
    y[b] = std::unordered_map<IWString, ASubstituent, IWStringHash>();

    y[b].emplace(usmi, msta.molecule_name());
    auto qq = y[b].find(usmi);
    qq->second.set_first_molecule(msta.molecule_name());
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
    cerr << "After starting new bit " << y[b].size() << " items stored\n";
#endif
  } else  // bit has been seen before
  {
    auto& usmi2subs = (*f).second;

    //  unordered_map<IWString, ASubstituent, IWStringHash>::iterator f2 =
    //  usmi2subs.find(substituent.unique_smiles());

    auto f2 = usmi2subs.find(usmi);

    if (f2 == usmi2subs.end()) {
      usmi2subs.emplace(usmi, msta.molecule_name());
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
      cerr << "For bit " << b << " new item " << substituent.unique_smiles() << ", name "
           << msta.molecule_name() << '\n';
      auto qq = usmi2subs.find(usmi);
      qq->second.set_first_molecule(msta.molecule_name());
      cerr << "new key '" << qq->first << "' new value " << qq->second.first_molecule()
           << '\n';
#endif
    } else {
      if (_concatenate_all_examples) {
        f2->second.extra(msta.molecule_name());
      } else {
        (*f2).second.extra();
      }
#ifdef DEBUG_ASSOCIATE_SUBSTITUENT_WITH_BIT
      cerr << "Extra instance of " << substituent.unique_smiles() << " now has "
           << (*f2).second.number_instances() << '\n';
#endif
    }
  }

  return;
}

int
SubstituentIdentification::_id_attch_pt_and_make_substituent_associations(
    Molecule& anchor, bond_type_t bt, Molecule& substituent,
    Molecule_Specific_Temporary_Arrays& msta) {
  std::unique_ptr<uint32_t[]> b = std::make_unique<uint32_t[]>(_shell_radius + 1);

  const atom_number_t a1 = anchor.atom_with_isotope(1);

  if (INVALID_ATOM_NUMBER == a1) {
    cerr << "SubstituentIdentification::no atom with isotope 1 in fragment '"
         << anchor.smiles() << "'\n";
    return 0;
  }

  int max_radius_formed = _compute_environment01(anchor, a1, _shell_radius, b.get(), msta);

  for (int i = 0; i <= max_radius_formed; ++i) {
    //  cerr << "Radius " << i << " bit " << b[i] << '\n';
    _associate_substituent_with_bit(i, b[i], bt, substituent, msta);
  }

  return 1;
}

int
convert_to_atom_number_in_parent(Molecule& m, const atom_number_t zatom) {
  // cerr << "Atom " << zatom << " was atom " << *(reinterpret_cast<const
  // atom_number_t *>(m.user_specified_atom_void_ptr(zatom))) << " in parent\n";

  return *(reinterpret_cast<const atom_number_t*>(m.user_specified_atom_void_ptr(zatom)));
}

int
SubstituentIdentification::_compute_environment01(
    Molecule& m, const atom_number_t zatom, int rmax, uint32_t* rc,
    Molecule_Specific_Temporary_Arrays& msta) const {
  const Atom* a = m.atomi(zatom);

  const atom_number_t zip = convert_to_atom_number_in_parent(m, zatom);

  rc[0] = msta.atom_type(zip);

  if (0 == rmax) {
    return 0;
  }

  uint32_t s = rc[0] + 127516;

  msta.set_processing_status(0);
  msta.set_processing_status(zatom, PROCESSING_FINISHED);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i) {
    const Bond* b = a->item(i);

    const auto a1 = b->other(zatom);

    const auto a1p = convert_to_atom_number_in_parent(m, a1);

    int bc = msta.bond_constant(*b, zip, a1p);

    s += bc * msta.atom_type(a1p);

    msta.set_processing_status(a1, NEXT_TIME);  // may not be needed
  }

  rc[1] = s;

  if (1 == rmax) {  // we are done
    return 1;
  }

  return _compute_environment2x(m, 2, rmax, rc, msta);
}

int
SubstituentIdentification::_compute_environment2x(
    Molecule& m, int radius, int rmax, uint32_t* rc,
    Molecule_Specific_Temporary_Arrays& msta) const {
  const auto matoms = m.natoms();

  uint32_t s = 0;

  for (auto i = 0; i < matoms; ++i) {
    if (NEXT_TIME != msta.processing_status(i)) {
      continue;
    }

    msta.set_processing_status(i, PROCESSING_FINISHED);

    const auto iip = convert_to_atom_number_in_parent(m, i);

    const Atom* a = m.atomi(i);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      const auto k = b->other(i);

      if (0 != msta.processing_status(k)) {
        continue;
      }

      msta.set_processing_status(k, READY_TO_PROCESS);

      const auto kip = convert_to_atom_number_in_parent(m, k);

      if (_precise_fingerprints) {
        s += msta.atom_type(iip) *
             (msta.bond_constant(*b, iip, kip) + 17 * msta.atom_type(kip));
      } else {
        s += msta.bond_constant(*b, iip, kip) * msta.atom_type(kip);
      }
    }
  }

  if (0 == s) {
    return radius - 1;
  }

  rc[radius] = rc[radius - 1] * radius + s;

  if (radius == rmax) {
    return radius;
  }

  for (auto i = 0; i < matoms; ++i) {
    if (READY_TO_PROCESS == msta.processing_status(i)) {
      msta.set_processing_status(i, NEXT_TIME);
    }
  }

  return _compute_environment2x(m, radius + 1, rmax, rc, msta);
}

int
SubstituentIdentification::_divide_molecule(Molecule& m, const atom_number_t a1,
                                            const atom_number_t a2,
                                            Molecule_Specific_Temporary_Arrays& msta) {
  const bond_type_t put_back = m.btype_between_atoms(a1, a2);

  if (SINGLE_BOND != put_back) {
    msta.temporarily_saturate(a1, a2);
  }

  m.remove_bond_between_atoms(a1, a2);
  m.recompute_implicit_hydrogens(a1);
  m.recompute_implicit_hydrogens(a2);
  m.set_isotope(a1, 1);
  m.set_isotope(a2, 1);
  resizable_array_p<Molecule> components;
  m.create_components(components);
  // cerr << "From " << x << " created " << components[0]->smiles() << " " <<
  // components[1]->smiles() << '\n';

  const auto atom_count_0 = components[0]->natoms();
  const auto atom_count_1 = components[1]->natoms();

  if (_ok_atom_count(atom_count_0) && atom_count_1 >= _min_residual_size) {
    _id_attch_pt_and_make_substituent_associations(*(components[1]), put_back,
                                                   *(components[0]), msta);
  }

  if (_ok_atom_count(atom_count_1) && atom_count_0 >= _min_residual_size) {
    _id_attch_pt_and_make_substituent_associations(*(components[0]), put_back,
                                                   *(components[1]), msta);
  }

  m.set_isotope(a1, 0);
  m.set_isotope(a2, 0);
  m.add_bond(a1, a2, put_back);

  if (SINGLE_BOND != put_back) {
    msta.undo_temporary_saturation(a1, a2);
  }

  return 1;
}

int
SubstituentIdentification::_build_database(Molecule& m,
                                           IWString_and_File_Descriptor& output) {
  const auto matoms = m.natoms();

  Molecule_Specific_Temporary_Arrays msta(matoms);

  std::unique_ptr<int[]> numbers = std::make_unique<int[]>(matoms);

  _initialise_msta(m, numbers.get(), msta);

  Molecule mcopy(m);

  for (int i = 0; i < matoms; i++) {
    if (0 == msta.atom_type(i)) {
      continue;
    }

    const Atom* a = m.atomi(i);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      if (b->nrings()) {
        continue;
      }

      const auto k = b->other(i);

      if (k < i) {
        continue;
      }

      _divide_molecule(mcopy, i, k, msta);
    }
  }

  static const Element* hydrogen = get_element_from_atomic_number(1);

  if (_process_hydrogen_substituents) {
    for (int i = 0; i < matoms; ++i) {
      if (0 == msta.atom_type(i)) {
        continue;
      }

      const Atom* a = m.atomi(i);

      const auto acon = a->ncon();

      if (1 == acon) {
        continue;
      }

      if (0 == m.hcount(i)) {
        continue;
      }

      mcopy.add(hydrogen);
      mcopy.add_bond(i, matoms, SINGLE_BOND);
      _divide_molecule(mcopy, i, matoms, msta);
      mcopy.remove_atom(matoms);
    }
  }

  return 1;
}


int
SubstituentIdentification::_process_molecules(data_source_and_type<Molecule>& input,
                                              UsedDuringLookups& lookup_data,
                                              IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    _molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    lookup_data.PrepareForNewMolecule();

    // cerr << "_do_matched_pairs " << _do_matched_pairs << '\n';
    int rc;
    if (_do_matched_pairs) {
      rc = _matched_pairs_qsar(*m, output);
    } else if (_anchor_query.size() || _default_new_molecule_starting_points) {
      rc = _look_for_new_substituents(*m, lookup_data, output);
    } else {
      rc = _build_database(*m, output);
    }

    if (0 == rc) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);

    if (_report_progress()) {
      if (_anchor_query.number_elements() || _default_new_molecule_starting_points) {
        _do_create_molecules_report(cerr);
      } else {
        _do_build_database_report(cerr);
      }
    }
  }

  return 1;
}

void
SubstituentIdentification::_do_build_database_report(std::ostream& os) const {
  os << "Processed " << _molecules_read << " molecules,";
  for (auto i = 0; i <= _shell_radius; ++i) {
    const auto& bi = _bit[i];

    int s = 0;
    for (auto j = bi.cbegin(); j != bi.cend(); ++j) {
      s += (*j).second.size();
    }
    os << " rad " << i << ' ' << s;
  }
  os << "\n";
}

void
SubstituentIdentification::_do_create_molecules_report(std::ostream& os) const {
  cerr << "Read " << _molecules_read << " molecules, " << _molecules_hitting_queries
       << " hit queries. Produced " << _molecules_written << ", "
       << _invalid_valences_ignored << " invalid valences ignored\n";
  if (_molecules_producing_too_many_new_molecules) {
    cerr << _molecules_producing_too_many_new_molecules
         << " would have produced more than " << _max_molecules_per_input_molecule
         << " new molecules\n";
  }

  return;
}


int
SubstituentIdentification::_process_molecules(const char* fname, FileType input_type,
                                              UsedDuringLookups& lookup_data,
                                              IWString_and_File_Descriptor& output) {
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

  if (_verbose > 1) {
    input.set_verbose(1);
  }

  return _process_molecules(input, lookup_data, output);
}

static int
opendb(Db& db, const char* dbname, DBTYPE dbtype, int flags, int mode) {
  int rc = db.open(NULL, dbname, NULL, dbtype, flags, mode);

  if (0 != rc) {
    cerr << "Cannot open database '" << dbname << "'\n";
    db.err(rc, "");
    return 0;
  }

  return 1;
}

static int
opendb_read(Db& db, const char* dbname) {
  return opendb(db, dbname, DB_UNKNOWN, DB_RDONLY, S_IREAD);
}

static int
opendb_write(Db& db, const char* dbname) {
  return opendb(db, dbname, DB_BTREE, DB_CREATE, S_IREAD | S_IWRITE | S_IRGRP | S_IROTH);
}

#define RADIUS_KEY "_RADIUS"

int
SubstituentIdentification::_store_radius() {
  Dbt dbkey((void*)RADIUS_KEY, strlen(RADIUS_KEY));
  IWString tmp;
  tmp << _shell_radius;

  Dbt zdata(const_cast<char*>(tmp.rawdata()), tmp.length());

  const auto s = _dbs[0]->put(NULL, &dbkey, &zdata, 0);
  if (0 == s) {
    return 1;
  }

  cerr << "SubstituentIdentification::_store_radius:Berkeley database Put "
          "operation failed\n";
  _dbs[0]->err(s, "");

  return 0;
}

int
SubstituentIdentification::_get_radius() {
  Dbt dbkey((void*)RADIUS_KEY, strlen(RADIUS_KEY));
  Dbt zdata;

  int s = _dbs[0]->get(NULL, &dbkey, &zdata, 0);

  if (0 != s) {
    cerr << "SubstituentIdentification::_get_radius:Berkeley database get "
            "operation failed\n";
    _dbs[0]->err(s, "");
    return 0;
  }

  const_IWSubstring tmp((const char*)zdata.get_data(), zdata.get_size());

  tmp.truncate_at_first(':');  // in case databases have been concatenated

  int rc;
  if (!tmp.numeric_value(rc) || rc <= 0) {
    cerr << "SubstituentIdentification::_get_radius:invalid radius '" << tmp << "'\n";
    return 0;
  }

  return rc;
}

int
SubstituentIdentification::_read_smiles_already_found(
    data_source_and_type<Molecule>& input, IW_STL_Hash_Set& smiles_already_found) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    smiles_already_found.insert(m->unique_smiles());
  }

  return smiles_already_found.size();
}

int
SubstituentIdentification::_read_smiles_already_found(
    const char* fname, IW_STL_Hash_Set& smiles_already_found) {
  const FileType input_type = discern_file_type_from_name(fname);
  if (input_type == FILE_TYPE_INVALID) {
    cerr << "read_smiles_already_found:do not know how to process '" << fname << "'\n";
    return 0;
  }

  data_source_and_type<Molecule> input(input_type, fname);

  if (!input.good()) {
    cerr << "Cannot open already formed molecules file '" << fname << "'\n";
    return 0;
  }

  return _read_smiles_already_found(input, smiles_already_found);
}

static void
DisplayDashYOptions(std::ostream& output) {
  output << R"(
 -Y dbproto        during building store serialized protos, then needed during read.
 -Y textproto=<fname> during building also write database contents as textproto form to <fname>.
 -Y exph           make implicit Hydrogens explicit.
 -Y x0             during building strip leading 0's from identifiers - may help with db size.
 -Y amap           during addition, label connections points with atom map numbers - rather than isotopes.
 -Y write_frag     during addition, write the fragment added.
 -Y maxgen=<n>     during addition, max number of molecules produced per input molecule (def 10000).
 -Y qat=<n>        during addition, number of query atoms to process - default the first query atom match.
 -Y rpt=<n>        report progress every <n> molecules processed.
 -Y textproto      write substituent_identification::Result textproto as output.
 -Y smitextproto   write smiles + substituent_identification::Result textproto as output.
 -Y matchedpairs   during lookup, identify local matched pairs.
 -Y allHmatch      all -H queries must match. By default, a product is OK if any of the -H queries match.
)";

  ::exit(0);
}

int
SubstituentIdentification::SetupTextProtoStream(const const_IWSubstring fname) {
  IWString tmp(fname);
  tmp.EnsureEndsWith(".textproto");
  if (! _stream_for_textproto.open(tmp)) {
    cerr << "Cannot open '" << tmp << "'\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Database contents written as textproto to '" << tmp << "'\n";
  }

  return 1;
}

int
SubstituentIdentification::operator()(int argc, char** argv) {
  Command_Line cl(argc, argv,
                  "vA:E:i:g:ld:R:P:M:m:cs:q:ab:r:C:H:K:fw:pIBu:kY:L:ye:X:hV:");

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
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, _verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      _usage(32);
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;

    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;

    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('a')) {
    _only_produce_molecules_at_biggest_radius = 1;

    if (_verbose) {
      cerr << "Will only produce molecules at the largest radius for which a "
              "match is found\n";
    }
  }

  if (cl.option_present('b')) {
    _only_add_bond = 0;

    IWString b;
    for (int i = 0; cl.value('b', b, i); ++i) {
      b.to_lowercase();

      if ('s' == b[0]) {
        _only_add_bond |= SINGLE_BOND;
      } else if ('d' == b[0]) {
        _only_add_bond |= DOUBLE_BOND;
      } else if ('t' == b[0]) {
        _only_add_bond |= TRIPLE_BOND;
      } else {
        cerr << "Unrecognised bond type specification '" << b << "'\n";
        _usage(1);
      }
    }

    if (_verbose) {
      cerr << "Bond mask " << _only_add_bond << '\n';
    }
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', _shell_radius) || _shell_radius < 0) {
      cerr << "The shell radius option (-R) must be a valid radius\n";
      _usage(2);
    }

    if (_verbose) {
      cerr << "Shell radius set to " << _shell_radius << '\n';
    }
  } else {
    _shell_radius = 1;
  }

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

    if (_verbose) {
      cerr << "Will discard product molecules with more than " << _max_atoms_in_product
           << " atoms\n";
    }
  }

  if (cl.option_present('p')) {
    _write_parent_molecule = 1;

    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('I')) {
    _remove_isotopes_from_product = 1;

    if (_verbose) {
      cerr << "Will remove isotopic labels from product molecules\n";
    }
  }

  _bit = new std::unordered_map<uint32_t, Usmi2Substituent>[_shell_radius + 1];

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!_atom_typing_specification.build(p)) {
      cerr << "Cannot initialise attachment atom type specification '" << p << "'\n";
      return 1;
    }
  } else {
    _atom_typing_specification.build("UST:AFUY");
  }

  if (cl.option_present('M')) {
    IWString m = cl.string_value('M');
    int value_set = 0;
    if (m.starts_with('+')) {
      m.remove_leading_chars(1);
      if (! m.numeric_value(_max_extra_atoms_added)) {
        cerr << "The maximum number of extra atoms added must be a whole +ve number\n";
        return 1;
      }
      if (_verbose) {
        cerr << "During addition replacement fragments can increas size by a max of " << 
                        _max_extra_atoms_added << " atoms\n";
      }
    } else {
      if (! m.numeric_value(_max_substituent_size) || _max_substituent_size < 1) {
        cerr << "The max substituent size option (-M) must be a whole +ve number\n";
        _usage(1);
      }
      if (_verbose) {
        cerr << "Will only examine substituents with " << _max_substituent_size
             << " or fewer atoms\n";
      }
      value_set = 1;
    }

    // If just +3 is specified, we need to turn off the checking for number of atoms.
    if (! value_set) {
      _max_substituent_size = 0;
    }
  }

  if (cl.option_present('m')) {
    IWString m = cl.string_value('m');
    int value_set = 0;
    if (m.starts_with('-')) {
      m.remove_leading_chars(1);
      if (! m.numeric_value(_max_atoms_lost_during_addition)) {
        cerr << "The max atoms lost during addition must be a whole +ve number\n";
        _usage(1);
      }
      if (_verbose) {
        cerr << "During addition, can only lose " << _max_atoms_lost_during_addition << 
                " atoms\n";
      }
    } else {
      if (! m.numeric_value(_min_substituent_size) || _min_substituent_size < 1) {
        cerr << "The minimum substituent size (-m) must be a whole +ve number\n";
        _usage(1);
      }

      if (_verbose) {
        cerr << "Will ignore substituents with fewer than " << _min_substituent_size
             << " atoms\n";
      }
      value_set = 1;
    }
    if (! value_set) {
      _min_substituent_size = 0;
    }
  }

  if (cl.option_present('L')) {
    if (!cl.value('L', _max_atoms_lost_during_build) ||
        _max_atoms_lost_during_build < 0) {
      cerr << "The max atoms lost during build option (-L) must be a whole +ve "
              "number\n";
      _usage(1);
    }
    if (_verbose) {
      cerr << "Parent molecules can lose a max of " << _max_atoms_lost_during_build
           << " atoms during building\n";
    }
  }

  if (cl.option_present('y')) {
    _default_new_molecule_starting_points = 1;
    if (_verbose) {
      cerr << "Will use default atom specifications for growing molecules\n";
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring e = cl.string_value('X');

    if (!_expt.read_data(e)) {
      cerr << "Cannot read experimental data data from '" << e << "'\n";
      return 1;
    }

    if (_verbose) {
      cerr << _expt.size() << " experimental values read from '" << e << "'\n";
    }
  }

  if (cl.option_present('e')) {
    const_IWSubstring e = cl.string_value('e');
    if (e == '.') {
      _concatenate_all_examples = 1;

      if (_verbose) {
        cerr << "Will store all examples of a substituent - not just first\n";
      }
    } else if (! SetMaxInstancesToConcatenate(e, _verbose)) {
      cerr << "The maximum number of instances to store (-e) must be a whole +ve number\n";
      return 1;
    }
  }

  if (cl.option_present('z')) {
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
      std::unique_ptr<Substructure_Query> t = std::make_unique<Substructure_Query>();
      if (!t->create_from_smarts(s)) {
        cerr << "Cannot parse smarts '" << s << "'\n";
        return 2;
      }
      _anchor_query << t.release();
    }
  }

  for (int i = 0; i < _anchor_query.number_elements(); ++i) {
    _anchor_query[i]->set_find_one_embedding_per_atom(1);
  }

  if (cl.option_present('k')) {
    _break_molecule_at_first_two_matched_atoms = 1;

    if (_verbose) {
      cerr << "Will break the input molecule at the first two matched atoms\n";
    }
  }

  if (cl.option_present('u')) {
    if (!cl.value('u', _min_examples_needed_for_addition) ||
        _min_examples_needed_for_addition < 1) {
      cerr << "The minimum number of example structures needed for addition "
              "(-u) must be a whole +ve number\n";
      _usage(1);
    }

    if (_verbose) {
      cerr << "Will only add a substituent if there are "
           << _min_examples_needed_for_addition << " or more examples in the db\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "dbproto") {
        _store_protos = 1;
        if (_verbose) {
          cerr << "Database will contain serialized protos\n";
        }
      } else if (y.starts_with("textproto=")) {
        y.remove_leading_chars(10);
        if (! SetupTextProtoStream(y)) {
          cerr << "Cannot initialise textproto file '" << y << "'\n";
          return 1;
        }
      } else if (y == "amap") {
        _apply_atom_map_labels = 1;
        if (_verbose) {
          cerr << "Will apply atom map labels\n";
        }
      } else if (y == "write_frag") {
        _write_fragments_added = 1;
        if (_verbose) {
          cerr << "Will write the fragment added during molecule creation\n";
        }
      } else if (y == "x0") {
        _strip_leading_zeros = 1;

        if (_verbose) {
          cerr << "Will strip leading zeros from identifiers\n";
        }
      } else if (y == "exph") {
        _make_implicit_hydrogens_explicit = 1;

        if (_verbose) {
          cerr << "Will make implicit Hydrogens explicit\n";
        }
      } else if (y.starts_with("maxgen=")) {
        y.remove_leading_chars(7);
        if (! y.numeric_value(_max_molecules_per_input_molecule) ||
            _max_molecules_per_input_molecule < 1) {
          cerr << "The maximum number of molecules produced per input molecule "
                  "must be a valid +ve number\n";
          return 2;
        }

        if (_verbose) {
          cerr << "Will produce a maximum of " << _max_molecules_per_input_molecule
               << " molecules per input molecule\n";
        }
      } else if (y.starts_with("qat=")) {
        y.remove_leading_chars(4);
        if (! y.numeric_value(_matched_atoms_to_process) || _matched_atoms_to_process < 1) {
          cerr << "The number of matched atoms to process qat= must be a whole +ve "
                  "number\n";
          return 1;
        }
        if (_verbose) {
          cerr << "Will process the first " << _matched_atoms_to_process
               << " matched atoms in each query\n";
        }
      } else if (y.starts_with("rpt=")) {
        y.remove_leading_chars(4);
        uint32_t rpt;
        if (! y.numeric_value(rpt)) {
          cerr << "Invalid report progress directive '" << y << "'\n";
          return 1;
        }
        _report_progress.set_report_every(rpt);
        if (_verbose) {
          cerr << "Will eport progress every " << rpt << " molecules\n";
        }
      } else if (y == "textproto") {
        _output_is_textproto = 1;
        if (_verbose) {
          cerr << "Will write textproto output\n";
        }
      } else if (y == "smitextproto") {
        _output_is_textproto = 1;
        _prepend_smiles_to_textproto = 1;
        if (_verbose) {
          cerr << "Will write product smiles and textproto\n";
        }
      } else if (y == "matchedpairs") {
        _do_matched_pairs = 1;
        if (_verbose) {
          cerr << "Will identify local mached pairs\n";
        }
      } else if (y.starts_with("maxlen=")) {
        y.remove_leading_chars(7);
        if (! y.numeric_value(_max_length_substituent) || _max_length_substituent < 1) {
          cerr << "The laxlen directive must be a whole +ve number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will only use fragments as long as " << _max_length_substituent << 
                  " bonds in length\n";
        }
      } else if (y == "allHmatch") {
        _all_must_have_queries_must_match = 1;
        if (_verbose) {
          cerr << "All must have queries (-H) must match, rather than any of them\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  if (cl.option_present('B')) {
    _remove_chirality = 1;
  }

  if (_anchor_query.empty() && !cl.option_present('B') &&
      0 == _default_new_molecule_starting_points) {
    cerr << "No queries specified (make new molecules), but the -B option for "
            "building not specified, cannot continue\n";
    _usage(2);
  }

  // Not building a database, but did not specify a max substituent size, just
  // give them whatever is in the database

  if (!cl.option_present('B') && !cl.option_present('M')) {
    _max_substituent_size = std::numeric_limits<int>::max();
  }

  if (cl.option_present('h')) {
    _process_hydrogen_substituents = 1;
    if (_verbose) {
      cerr << "Will also generate Hydrogen substituents\n";
    }
  }

  if (cl.option_present('H')) {
    if (!process_queries(cl, _replacements_must_contain, _verbose, 'H')) {
      cerr << "Cannot process queries substituents must contain (-H)\n";
      return 2;
    }
  }
  if (cl.option_present('K')) {
    if (!process_queries(cl, _replacements_must_not_contain, _verbose, 'K')) {
      cerr << "Cannot process queries substituents must NOT contain (-K)\n";
      return 2;
    }
  }

  if (_verbose && _anchor_query.size()) {
    cerr << "Defined " << _anchor_query.size()
         << " queries to identify attachment points\n";
  }

  if (cl.option_present('f')) {
    _precise_fingerprints = 0;

    if (_verbose) {
      cerr << "Will generate fingerprints that allow imprecise matching\n";
    }
  }

  if (cl.option_present('w')) {
    if (!cl.value('w', _min_residual_size) || _min_residual_size < 1) {
      cerr << "The minimum residual atom count option (-w) must be a whole +ve "
              "number\n";
      _usage(2);
    }

    if (_verbose) {
      cerr << "Will not produce a fragment if the residual contains fewer than "
           << _min_residual_size << " atoms\n";
    }
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify database to build via the -d option\n";
    _usage(1);
  }

  auto dcount = cl.option_count('d');

  if (dcount > 1 && _anchor_query.empty()) {
    cerr << "Can only build one database at a time\n";
    return 1;
  }

  if (1 == dcount && _anchor_query.empty() &&
      0 == _default_new_molecule_starting_points)  // great, building a database
                                                   // and no queries
  {
    _dbs = new Db*[1];
    _ndb = 1;
    _dbs[0] = new Db(0, DB_CXX_NO_EXCEPTIONS);

    const char* dbname = cl.option_value('d');
    if (!opendb_write(*_dbs[0], dbname)) {
      return 1;
    }

    if (!_store_radius()) {
      return 1;
    }
  } else { // lookup
    _ndb = dcount;
    _dbs = new Db*[_ndb];

    for (auto i = 0; i < _ndb; ++i) {
      const char* dbname = cl.option_value('d', i);

      _dbs[i] = new Db(0, DB_CXX_NO_EXCEPTIONS);

      if (!opendb_read(*_dbs[i], dbname)) {
        return 1;
      }
    }

    if (_verbose) {
      cerr << "Opened " << _ndb << " databases for lookup\n";
    }
  }

  if (cl.option_present('r')) {
    auto stored = _get_radius();
    if (stored < 1) {
      cerr << "SubstituentIdentification:cannot retrieve radius from DB, beware\n";
      stored = _min_shell_radius;  // just use what the user entered
    }

    if (_min_shell_radius > stored) {
      cerr << "Min shell radius " << _min_shell_radius
           << " too large, database built with max radius " << stored << '\n';
      return 1;
    }

    _shell_radius = stored;

    if (_verbose) {
      cerr << "Will only look at shells of radius " << _min_shell_radius << " or more\n";
    }

    if (_shell_radius < _min_shell_radius) {
      _shell_radius = _min_shell_radius;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      _usage(6);
    }
  } else if (1 == cl.size() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  set_copy_name_in_molecule_copy_constructor(1);
  set_copy_user_specified_atom_void_ptrs_during_create_subset(1);

  IW_STL_Hash_Set smiles_already_found;

  if (cl.option_present('V')) {
    const char* v = cl.option_value('V');

    if (!_read_smiles_already_found(v, smiles_already_found)) {
      cerr << "Cannot read already formed smiles '" << v << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << smiles_already_found.size()
           << " already formed molecules to avoid from '" << v << "'\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    _usage(2);
  }

  UsedDuringLookups lookup_data;

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!_process_molecules(cl[i], input_type, lookup_data, output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (_anchor_query.empty() &&
      0 == _default_new_molecule_starting_points) {
    _write_in_memory_hashes_to_database();
  }

  if (_verbose) {
    cerr << "Read " << _molecules_read << " molecules\n";
    if (_anchor_query.size() || _default_new_molecule_starting_points) {
      _do_create_molecules_report(cerr);
    } else {
      for (auto i = 0; i <= _shell_radius; ++i) {
        const auto& b = _bit[i];

        cerr << "At radius " << i << " found " << b.size() << " items\n";
      }
      _do_build_database_report(cerr);
      cerr << "Stored " << _keys_stored << " bits, " << _pairs_stored << " pairs\n";
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  SubstituentIdentification SubstituentIdentification;

  return SubstituentIdentification(argc, argv);
}
