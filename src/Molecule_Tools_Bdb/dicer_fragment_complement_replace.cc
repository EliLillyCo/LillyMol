// Read texproto output from dicer, including the complementary fragments.
// For those fragments that meet the specified constraints, lookup
// other complementary fragments in a database and genrate new molecules.

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "db_cxx.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Tools/mformula.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/mformula.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools_Bdb/dicer_fragment_complement.pb.h"
#else
#include "Molecule_Tools/dicer_fragments.pb.h"
#include "dicer_fragment_complement.pb.h"
#endif

namespace dicer_fragment_replace {

using std::cerr;

struct Attachment {
  atom_number_t atom;
  isotope_t isotope;
};

using AttachmentList = std::vector<Attachment>;

struct AttachmentPair {
  atom_number_t frag_atom;
  atom_number_t comp_atom;
};

using AttachmentPairing = std::vector<AttachmentPair>;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(
 -d <dbname>    one or more BerkeleyDB databases containing mappings from unique smiles to complementary fragments
 -v             verbose output
)";
  // clang-format on

  ::exit(rc);
}

class Replacement {
  private:
    // As read during statup.
    std::unique_ptr<dicer_data::DicerFragment> _proto;

    Molecule _m;

    mformula::MFormula _mformula;

    int _aromatic_atom_count;
    int _nrings;

  public:
    Replacement();

    int Build(dicer_data::DicerFragment* proto);
};

Replacement::Replacement() {
  _aromatic_atom_count = 0;
  _nrings = 0;
}

int
Replacement::Build(dicer_data::DicerFragment* proto) {
  _proto.reset(proto);

  if (! _m.build_from_smiles(_proto->smi())) {
    cerr << "Replacement::Build:invalid smiles " << _proto->ShortDebugString() << '\n';
    return 0;
  }

  _nrings = _m.nrings();
  _aromatic_atom_count = _m.aromatic_atom_count();

  return 1;
}

class SetOfReplacements {
  private:
    resizable_array_p<Replacement> _replacements;

  public:
    int Build(IWString& fname);
    int Build(iw_tf_data_record::TFDataReader& reader);
};

int
SetOfReplacements::Build(IWString& fname) {
  iw_tf_data_record::TFDataReader reader(fname);

  if (! reader.good()) {
    cerr << "SetOfReplacements::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(reader);
}

int
SetOfReplacements::Build(iw_tf_data_record::TFDataReader& reader) {

  while (true) {
    std::unique_ptr<dicer_data::DicerFragment> proto =
                reader.ReadProtoPtr<dicer_data::DicerFragment>();
    if (! proto) {
      return _replacements.size();
    }
  }
}

#ifdef PMD_NOT_USED_HERE
struct PerAtomData {
  Molecule_to_Match target;
};

class PerMoleculeData {
  private:
    Molecule& _m;
    int _natoms;

    // A cross reference between atom numbers in a fragment and atom numbers in the parent.
    int * _xref_f2p;
    // A cross reference between atom numbers in the parent and atom numbers in a fragment.
    int * _xref_p2f;

    uint32_t* _atype;

  public:
    PerMoleculeData(Molecule& m);
    ~PerMoleculeData();

    int AssignAtomTypes(Atom_Typing_Specification& ats);

    const uint32_t* atype() const {
      return _atype;
    }

    void EstablishXref(const Set_of_Atoms& embedding);

    int* xref_f2p() const {
      return _xref_f2p;
    }
    int* xref_p2f() const {
      return _xref_p2f;
    }
};

PerMoleculeData::PerMoleculeData(Molecule& m) : _m(m) {
  _natoms = _m.natoms();

  _atype = nullptr;

  _xref_f2p = new int[_natoms];
  _xref_p2f = new int[_natoms];
}

PerMoleculeData::~PerMoleculeData() {
  delete [] _xref_f2p;
  delete [] _xref_p2f;

  if (_atype != nullptr) {
    delete [] _atype;
  }
}

int
PerMoleculeData::AssignAtomTypes(Atom_Typing_Specification& ats) {
  if (_atype == nullptr) {
    _atype = new uint32_t[_natoms];
  }

  return ats.assign_atom_types(_m, _atype);
}

void
PerMoleculeData::EstablishXref(const Set_of_Atoms& embedding) {
  std::fill_n(_xref_f2p, _natoms, 0);
  std::fill_n(_xref_p2f, _natoms, 0);

  const int n = embedding.number_elements();
  for (int i = 0; i < n; ++i) {
    atom_number_t j = embedding[i];
    _xref_f2p[i] = j;
    _xref_p2f[j] = i;
  }
}
#endif // PMD_NOT_USED_HERE

struct PerMoleculeData {
  // If we are writing parent molecules, high up in the call tree
  // this will be initialised with what should be written as the
  // parent. When MaybeWriteProduct writes a product, it will check to see
  // if this is non-empty. If non-empty, it will write it and resize it
  // to empty.
  // That way the parent molecule information gets written once - and only
  // if products are generated.
  IWString write_parent_molecule;

  uint32_t products_generated = 0;

  std::unique_ptr<uint32_t[]> atype;
};

class SetOfDatabases {
  private:
    int _verbose;

    dicer_replace_complement::Options _options;

    // As dicer fragments are read from the input, we only process those fragments
    // with at least this many atoms.
    int _min_starting_fragment_natoms = 1;
    uint32_t _fragments_too_small = 0;
    uint32_t _duplicates_generated = 0;

    resizable_array_p<Substructure_Query> _fragment_must_contain;
    resizable_array_p<Substructure_Query> _fragment_must_not_contain;
    resizable_array_p<Substructure_Query> _complement_must_contain;
    resizable_array_p<Substructure_Query> _complement_must_not_contain;

    // If the input was generated with '-B nousmi' then we can match atoms without
    // needing to do any substructure searching.
    int _input_has_atom_maps;

    resizable_array_p<Db> _database;

    uint64_t _fragments_found_in_db;
    uint64_t _fragments_not_found_in_db;;
    uint64_t _ok_fragment_failed = 0;
    uint64_t _ok_new_complement_failed = 0;
    uint64_t _ok_complementary_failed = 0;

    uint64_t _molecules_processed;

    Atom_Typing_Specification _atom_typing;

    uint32_t _max_products_per_parent;

    int _remove_isotopes;

    // THe -p <string> sets this value.
    IWString _parent_molecule_string;

    char _write_initial_complementary_fragment;

    extending_resizable_array<uint32_t> _complements_retrieved;

    absl::flat_hash_set<IWString> _seen;

    // private functions
    int ProcessFragment(Molecule& parent,
                        PerMoleculeData& pmd,
                        const dicer_data::DicerFragment& fragment,
                        IWString_and_File_Descriptor& output);
    int AddFragments(PerMoleculeData& pmd,
                Molecule& frag,
                Molecule& initial_comp,
                const resizable_array_p<Molecule>& complements,
                const std::vector<int>& counts,
                IWString_and_File_Descriptor& output);
    int MakeMolecules(PerMoleculeData& pmd,
                const Molecule& frag,
                Molecule& initial_comp,
                const Molecule& comp,
                int count,
                const AttachmentList& frag_attachments,
                const AttachmentList& comp_attachments,
                IWString_and_File_Descriptor& output);
    int MaybeWriteProduct(PerMoleculeData& pmd, Molecule& m,
                        int count,
                        Molecule& initial_comp,
                        IWString_and_File_Descriptor& output);

    bool Seen(Molecule& m);

    int OkComplementary(const Molecule& parent,
                        const Molecule& frag,
                        Molecule& comp);
    int OkFragment(const Molecule& parent, Molecule& fragment);
    int GetComplements(const IWString& usmi,
                       resizable_array_p<Molecule>& complements,
                       std::vector<int>& counts);
    int GetComplementsInner(const const_IWSubstring& fromdb,
               resizable_array_p<Molecule>& complements,
               std::vector<int>& counts);
    int OkNewComplement(const Molecule& initial_comp,
                std::unique_ptr<mformula::MFormula>& initial_formula,
                Molecule& new_comp) const;
    int  OkProduct(Molecule& m);

    int InitialiseQueries();

  public:
    SetOfDatabases();
    ~SetOfDatabases();

    int Initialise(Command_Line& cl);

    int Process(const dicer_data::DicedMolecule& proto, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

SetOfDatabases::SetOfDatabases() {
  _verbose = 0;
  _molecules_processed = 0;
  _input_has_atom_maps = 0;
  _min_starting_fragment_natoms = 0;
  _fragments_too_small = 0;
  _remove_isotopes = 0;
  _write_initial_complementary_fragment = '\0';
  _max_products_per_parent = std::numeric_limits<uint32_t>::max();
  _fragments_found_in_db = 0;
  _fragments_not_found_in_db = 0;
}

static void
DisplayDashBOptions() {
  cerr << R"( -B nousmi         input is from dicer with the -B nousmi option
)";

  ::exit(0);
}

SetOfDatabases::~SetOfDatabases() {
  for (Db* db : _database) {
    db->close(0);
  }
}

int
SetOfDatabases::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');
  if (! cl.option_present('d')) {
    cerr << "Must specify one or more databases via the -d option\n";
    return 0;
  }

  set_copy_name_in_molecule_copy_constructor(1);

  int flags = DB_RDONLY;
  DBTYPE dbtype = DB_UNKNOWN;
  int mode = 0;
  DbEnv* env = NULL;

  IWString dbname;
  for (int i = 0; cl.value('d', dbname, i); ++i) {
    std::unique_ptr<Db> db = std::make_unique<Db>(env, DB_CXX_NO_EXCEPTIONS);
    if (int rc = db->open(NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);
        rc != 0) {
      cerr << "Cannot open '" << dbname << "' ";
      db->err(rc, "");
    }
    _database << db.release();
  }

  if (_verbose) {
    cerr << "Opened " << _database.size() << " dicer precedent databases\n";
  }

  if (cl.option_present('B')) {
    const_IWSubstring b;
    for (int i = 0; cl.value('B', b, i); ++i) {
      if (b == "nousmi") {
        _input_has_atom_maps = 1;
        if (_verbose) {
          cerr << "Assume that input has atom map information\n";
        }
      } else if (b == "help") {
        DisplayDashBOptions();
      } else {
        cerr << "Unrecognised -B qualifier '" << b << "'\n";
        DisplayDashBOptions();
      }
    }
  }

  if (cl.option_present('P')) {
    const const_IWSubstring p = cl.string_value('P');
    if (! _atom_typing.build(p)) {
      cerr << "Invalid atom typing '" << p << "'\n";
      return 0;
    }
  }

  if (cl.option_present('I')) {
    _remove_isotopes = 1;
    if (_verbose) {
      cerr << "Will remove isotopes from products\n";
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    std::optional<dicer_replace_complement::Options> opts = 
         iwmisc::ReadTextProto<dicer_replace_complement::Options>(fname);
    if (! opts) {
      cerr << "SetOfDatabases::Initialise:cannot read configuration file '" << fname << "'\n";
      return 0;
    }

    _options = *opts;  // maybe move would also be ok, but will not matter.

    if (! InitialiseQueries()) {
      cerr << "SetOfDatabases::Initialise:cannot initialise queries\n";
      return 0;
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_starting_fragment_natoms) || _min_starting_fragment_natoms < 1) {
      cerr << "The minimum starting fragment size (-m) must be a whole +ve number\n";
      return 1;
    }
    if (_verbose) {
      cerr << "Will skip dicer fragments with fewer than " << _min_starting_fragment_natoms << " atoms\n";
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_products_per_parent) || _max_products_per_parent < 1) {
      cerr << "SetOfDatabases::Initialise:the max products per parent (-x) option must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will generate a max of " << _max_products_per_parent <<
              " products per starting molecule\n";
    }
  }

  if (cl.option_present('c')) {
    IWString c = cl.string_value('c');
    if (! char_name_to_char(c)) {
      cerr << "Invalid initial fragment separator specification '" << c << "'\n";
      return 0;
    }

    _write_initial_complementary_fragment = c[0];
    if (_verbose) {
      cerr << "Will include the initial complementary fragment in the output\n";
    }
  }

  if (cl.option_present('p')) {
    _parent_molecule_string = cl.string_value('p');
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  return 1;
}

int
InitialiseQueries(const google::protobuf::RepeatedPtrField<std::string>& smarts, resizable_array_p<Substructure_Query>& queries) {
  for (const std::string& s : smarts) {
    std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();

    if (! q->create_from_smarts(s)) {
      cerr << "InitialiseQueries:invalid smarts '" << s << "'\n";
      return 0;
    }

    queries << q.release();
  }

  return queries.number_elements();
}

// We have read a configuration proto and need to convert anyof the smarts in
// there to queries.
int
SetOfDatabases::InitialiseQueries() {
  if (_options.fragment_conditions().smarts_must_have_size() > 0) {
    if (! dicer_fragment_replace::InitialiseQueries(_options.fragment_conditions().smarts_must_have(),
                _fragment_must_contain)) {
      return 0;
    }
  }

  if (_options.fragment_conditions().smarts_must_not_have_size() > 0) {
    if (! dicer_fragment_replace::InitialiseQueries(_options.fragment_conditions().smarts_must_not_have(),
                _fragment_must_not_contain)) {
      return 0;
    }
  }

  if (_options.complement_conditions().smarts_must_have_size() > 0) {
    if (! dicer_fragment_replace::InitialiseQueries(_options.complement_conditions().smarts_must_have(),
                _complement_must_contain)) {
      return 0;
    }
  }

  if (_options.complement_conditions().smarts_must_not_have_size() > 0) {
    if (! dicer_fragment_replace::InitialiseQueries(_options.complement_conditions().smarts_must_not_have(),
                _complement_must_not_contain)) {
      return 0;
    }
  }

  return 1;
}

int
SetOfDatabases::Report(std::ostream& output) const {
  output << "SetOfDatabases processed " << _molecules_processed << '\n';
  output << "Found " << _fragments_found_in_db << " fragments, did not find " <<
            _fragments_not_found_in_db << '\n';
  output << _ok_complementary_failed << " complements failed prerequisites\n";
  output << _ok_new_complement_failed << " ok new complement failed\n";
  output << _ok_fragment_failed << " ok fragment failed\n";

  if (_min_starting_fragment_natoms > 0) {
    output << _fragments_too_small << " fragments skipped for fewer than " <<
              _min_starting_fragment_natoms << " atoms\n";
  }

  for (int i = 0; i < _complements_retrieved.number_elements(); ++i) {
    if (_complements_retrieved[i]) {
      cerr << _complements_retrieved[i] << " generated " << i << " variants\n";
    }
  }
  return output.good();
}


int
SetOfDatabases::Process(const dicer_data::DicedMolecule& proto,
                        IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  Molecule m;
  if (! m.build_from_smiles(proto.smiles())) {
    cerr << "SetOfDatabases::Process:invalid smiles '" << proto.smiles() << "'\n";
    return 0;
  }

  m.set_name(proto.name());
  cerr << " read " << proto.name() << " with " << proto.fragment_size() << " fragmetns\n";

//if (!_atom_typing.active()) {
//  pmd.AssignAtomTypes(_atom_typing);
//}

  PerMoleculeData pmd;

  if (!_parent_molecule_string.empty()) {
    pmd.write_parent_molecule << m.smiles() << ' ' << proto.name()
                              << ' ' << _parent_molecule_string;
  }

  if (_atom_typing.active()) {
    pmd.atype = std::make_unique<uint32_t[]>(m.natoms());
    _atom_typing.assign_atom_types(m, pmd.atype.get());
  }

  int rc = 0;
  for (const dicer_data::DicerFragment& fragment : proto.fragment()) {
    rc += ProcessFragment(m, pmd, fragment, output);
  }

  ++_complements_retrieved[rc];

  return 1;
}


// `fragment` is a fragment in `parent`. `fragment` and the
// complement in `fragment` comprise `parent`.
// Look for complements in the database.
int
SetOfDatabases::ProcessFragment(Molecule& parent,
                        PerMoleculeData& pmd,
                        const dicer_data::DicerFragment& fragment,
                        IWString_and_File_Descriptor& output) {
  if (! fragment.has_comp()) {
    cerr << "No complementary fragment\n";
    return 0;
  }

  Molecule frag;
  if (! frag.build_from_smiles(fragment.smi())){
    cerr << "SetOfDatabases::Process:invalid smiles " << fragment.ShortDebugString() << '\n';
    return 0;
  }

  if (frag.natoms() < _min_starting_fragment_natoms) {
    cerr << "Too few atoms " << frag.natoms() << '\n';
    ++_fragments_too_small;
    return 1;
  }

  cerr << parent.smiles() << " parent smiles\n";
  cerr << frag.smiles() << " fragment smiles\n";
  if (! OkFragment(parent, frag)) {
    cerr << "OkFragment failed\n";
    ++_ok_fragment_failed;
    return 0;
  }

  Molecule comp;
  if (! comp.build_from_smiles(fragment.comp())) {
    cerr << "SetOfDatabases::Process:invalid complementary smiles\n";
    cerr << fragment.ShortDebugString() << '\n';
    return 0;
  }

  // Sanity check.
  if ((frag.natoms() + comp.natoms() != parent.natoms())) [[unlikely]] {
    cerr << "Atom count mismatch " << frag.natoms() << '+' << comp.natoms() << " != "
         << parent.natoms() << '\n';
    return 0;
  }

  cerr << comp.smiles() << " complementary smiles\n";
  if (! OkComplementary(parent, frag, comp)) {
    cerr << "OkComplementary failed\n";
    ++_ok_complementary_failed;
    return 0;
  }
  cerr << "OkComplementary worked\n";

  // Look up complements of this fragments in the database.
  resizable_array_p<Molecule> complements;
  std::vector<int> counts;
  if (! GetComplements(frag.unique_smiles(), complements, counts)) {
    cerr << "GetComplements did not return complements\n";
    return 0;
  }
  cerr << "GetComplements worked, have " << complements.size() << '\n';

  frag.set_name(parent.name());

  return AddFragments(pmd, frag, comp, complements, counts, output);
}

int
GetComplement(const const_IWSubstring& fromdb,
              resizable_array_p<Molecule>& complements) {
  const_IWSubstring id, smi;
  if (! fromdb.split(id, ':', smi) || id.empty() || smi.empty()) {
    return 0;
  }

  std::unique_ptr<Molecule> m = std::make_unique<Molecule>();
  if (! m->build_from_smiles(smi)) {
    cerr << "GetComplement:invalid smiles '" << smi << "'\n";
    return 0;
  }

  m->set_name(id);

  complements << m.release();

  return 1;
}

int
SetOfDatabases::GetComplementsInner(const const_IWSubstring& fromdb,
               resizable_array_p<Molecule>& complements,
               std::vector<int>& counts) {

  dicer_data::ComplementaryFragments proto;
  if (! proto.ParseFromArray(fromdb.data(), fromdb.nchars())) {
    cerr << "SetOfDatabases::GetComplementsInner:cannot decode proto\n";
    return 0;
  }

  for (const auto& complement : proto.frag()) {
    std::unique_ptr<Molecule> m = std::make_unique<Molecule>();
    if (! m->build_from_smiles(complement.smi())) {
      cerr << "SetOfDatabases::GetComplementsInner:invalid smiles '" << complement.smi() << "'\n";
      return 0;
    }
    m->set_name(complement.id());
    if (complement.n() == 0) {
      cerr << "Read complement with zero count " << m->aromatic_smiles() << '\n';
    }
    complements << m.release();
    counts.push_back(complement.n());
    if (complements.size() > _max_products_per_parent) {
      break;
    }
  }

  return complements.size();
}

int
SetOfDatabases::GetComplements(const IWString& usmi,
                               resizable_array_p<Molecule>& complements,
                               std::vector<int>& counts) {
  cerr << "Looking up " << usmi << "'\n";
  Dbt dbkey((void* ) (usmi.data()), usmi.length());

  for (Db* db : _database) {
    Dbt fromdb;
    if (db->get(NULL, &dbkey, &fromdb, 0) != 0) {
      continue;
    }

    cerr << usmi << " found in db\n";
    ++_fragments_found_in_db;
    const_IWSubstring tmp(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
    return GetComplementsInner(tmp, complements, counts);
  }

  ++_fragments_not_found_in_db;
  cerr << usmi << " not found\n";

  return 0;
}

// `atom_conditions` contains limits on the differences between the atoms
// in the parent and in a frgment - which might be either side of a split.
// Return 1 if the atom counts are OK wrt the constraints.
int
OkAtomDifferences(const dicer_replace_complement::Conditions& atom_conditions,
        uint32_t atoms_in_parent,
        uint32_t atoms_in_fragment) {
  cerr << atoms_in_parent << " atoms_in_parent " << atoms_in_fragment << " atoms_in_fragment\n";
  if (! atom_conditions.has_min_atoms_in_fragment()) {
  } else if (atoms_in_fragment < atom_conditions.min_atoms_in_fragment()) {
    return 0;
  }

  if (! atom_conditions.has_max_atoms_in_fragment()) {
  } else if (atoms_in_fragment > atom_conditions.max_atoms_in_fragment()) {
    return 0;
  }

  if (! atom_conditions.has_min_atoms_lost_parent_to_fragment()) {
  } else if (atoms_in_parent - atoms_in_fragment < atom_conditions.min_atoms_lost_parent_to_fragment()) {
    return 0;
  }

  if (! atom_conditions.has_max_atoms_lost_parent_to_fragment()) {
  } else if (atoms_in_parent - atoms_in_fragment > atom_conditions.max_atoms_lost_parent_to_fragment()) {
    return 0;
  }

  return 1;
}

// We have just read the diced form of `parent` and `fragment` is a fragment.
// Is it OK to start looking for complementary fragments, given that `fragment`
// will be the conserved part of `parent`.
int
SetOfDatabases::OkFragment(const Molecule& parent,
                           Molecule& fragment) {
  return 1;



  if (! _options.has_fragment_conditions()) {
    return 1;
  }

  const dicer_replace_complement::Conditions& atom_conditions = _options.fragment_conditions();

  const uint32_t atoms_in_parent = parent.natoms();
  const uint32_t atoms_in_fragment = fragment.natoms();

  if (! OkAtomDifferences(atom_conditions, atoms_in_parent, atoms_in_fragment)) {
    cerr << "OkAtomDifferences failed\n";
    return 0;
  }

  if (_fragment_must_contain.empty()) {
  } else if (! lillymol::AnyQueryMatches(fragment, _fragment_must_contain)) {
    return 0;
  }

  if (_fragment_must_not_contain.empty()) {
  } else if (lillymol::AnyQueryMatches(fragment, _fragment_must_not_contain)) {
    return 0;
  }

  return 1;
}

// Given that `frag` is a fragment in `parent` is it OK to
// begin searching for different complementary fragments.
// This is the first test.
int
SetOfDatabases::OkComplementary(const Molecule& parent,
                        const Molecule& frag,
                        Molecule& comp) {
  return 1;



  if (! _options.has_complement_conditions()) {
    return 1;
  }

  const dicer_replace_complement::Conditions& atom_conditions = _options.complement_conditions();

  const uint32_t atoms_in_parent = parent.natoms();
  const uint32_t atoms_in_complement = comp.natoms();

  if (! OkAtomDifferences(atom_conditions, atoms_in_parent, atoms_in_complement)) {
    return 0;
  }

  if (_complement_must_contain.empty()) {
  } else if (! lillymol::AnyQueryMatches(comp, _complement_must_contain)) {
    return 0;
  }

  if (_complement_must_not_contain.empty()) {
  } else if (lillymol::AnyQueryMatches(comp, _complement_must_not_contain)) {
    return 0;
  }

  return 1;
}

// `initial_comp` is the complementary fragment that was found in the
// starting molecule and we are proposing that it be switched with `new_comp`.
// Check to see if that is an OK change.
int
SetOfDatabases::OkNewComplement(const Molecule& initial_comp,
                std::unique_ptr<mformula::MFormula>& initial_formula,
                Molecule& new_comp) const {
  return 1;

  if (! _options.has_min_comp_formula_difference() &&
      ! _options.has_max_comp_formula_difference()) {
    return 1;
  }

  mformula::MFormula comp_mf;
  comp_mf.Build(new_comp);

  uint32_t diff = comp_mf.Diff(*initial_formula);
  if (! _options.has_min_comp_formula_difference()) {
  } else if (diff < _options.min_comp_formula_difference()) {
    return 0;
  }

  if (! _options.has_max_comp_formula_difference()) {
  } else if (diff > _options.max_comp_formula_difference()) {
    return 0;
  }

  return 1;
}

int
SetOfDatabases::OkProduct(Molecule& m) {
  return 1;
}


int
CompareAttachment(const Attachment& a1, const Attachment& a2) {
  if (a1.isotope < a2.isotope) {
    return -1;
  }
  if (a1.isotope > a2.isotope) {
    return 1;
  }
  if (a1.atom < a2.atom) {
    return -1;
  }
  if (a1.atom > a2.atom) {
    return 1;
  }
  return 0;
}

bool
AttachmentLess(const Attachment& a1, const Attachment& a2) {
  return CompareAttachment(a1, a2) < 0;
}

AttachmentList
GetAttachments(const Molecule& m) {
  AttachmentList result;

  const int matoms = m.natoms();
  result.reserve(matoms);

  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso = m.isotope(i);
    if (iso) {
      result.push_back(Attachment{static_cast<atom_number_t>(i), iso});
    }
  }

  std::sort(result.begin(), result.end(), AttachmentLess);

  return result;
}

IWString
AttachmentsToString(const AttachmentList& attachments) {
  IWString result;

  for (const Attachment& attachment : attachments) {
    if (! result.empty()) {
      result << ' ';
    }
    result << attachment.atom << ':' << attachment.isotope;
  }

  return result;
}

bool
SameIsotopeMultiset(const AttachmentList& frag_attachments,
                    const AttachmentList& comp_attachments) {
  if (frag_attachments.size() != comp_attachments.size()) {
    return false;
  }

  for (uint32_t i = 0; i < frag_attachments.size(); ++i) {
    if (frag_attachments[i].isotope != comp_attachments[i].isotope) {
      return false;
    }
  }

  return true;
}

void
GenerateAttachmentPairings(const AttachmentList& frag_attachments,
                           const AttachmentList& comp_attachments,
                           std::vector<AttachmentPairing>& pairings) {
  pairings.clear();

  AttachmentList comp_permutation = comp_attachments;
  do {
    AttachmentPairing pairing;
    pairing.reserve(frag_attachments.size());

    bool matched = true;
    for (uint32_t i = 0; i < frag_attachments.size(); ++i) {
      if (frag_attachments[i].isotope != comp_permutation[i].isotope) {
        matched = false;
        break;
      }

      pairing.push_back(AttachmentPair{frag_attachments[i].atom, comp_permutation[i].atom});
    }

    if (matched) {
      pairings.push_back(std::move(pairing));
    }
  } while (std::next_permutation(comp_permutation.begin(), comp_permutation.end(), AttachmentLess));
}

// We have looked up complementary fragments from the database. These are in
// `complements`. The parent molecule consists of two fragments, `frag` and
// `initial_comp`.
int
SetOfDatabases::AddFragments(
                PerMoleculeData& pmd,
                Molecule& frag,
                Molecule& initial_comp,
                const resizable_array_p<Molecule>& complements,
                const std::vector<int>& counts,
                IWString_and_File_Descriptor& output) {
  std::unique_ptr<mformula::MFormula> initial_formula;
  if (_options.has_min_comp_formula_difference() ||
      _options.has_max_comp_formula_difference()) {
    initial_formula = std::make_unique<mformula::MFormula>();
    initial_formula->Build(initial_comp);
  }

  AttachmentList frag_attachments = GetAttachments(frag);
  cerr << "AddFragments frag " << frag.smiles() << " attachments " << AttachmentsToString(frag_attachments) << '\n';

  const int initial_comp_number_fragments = initial_comp.number_fragments();

  int rc = 0;
  for (int i = 0; i < complements.number_elements(); ++i) {
    Molecule* new_comp = complements[i];
    cerr << "new_comp " << new_comp->smiles() << " initial_comp_number_fragments " << initial_comp_number_fragments << '\n';
    if (new_comp->number_fragments() != initial_comp_number_fragments) {
      cerr << "Fragment count mismatch\n";
      continue;
    }

    if (! OkNewComplement(initial_comp, initial_formula, *new_comp)) {
      cerr << "OkNewComplement failed\n";
      ++_ok_new_complement_failed;
      continue;
    }

    AttachmentList comp_attachments = GetAttachments(*new_comp);
    cerr << new_comp->smiles() << " attachments " << AttachmentsToString(comp_attachments) << '\n';

    if (frag_attachments.size() != comp_attachments.size()) {
      cerr << "Attachment size mismatch\n";
      continue;
    }
    if (! SameIsotopeMultiset(frag_attachments, comp_attachments)) {
      cerr << "Attachment isotope mismatch\n";
      continue;
    }

    rc += MakeMolecules(pmd, frag, initial_comp, *new_comp, counts[i], frag_attachments, comp_attachments, output);
  }

  return rc;
}

void
UnfixImplicitHydrogens(Molecule& m,
                       atom_number_t a1,
                       atom_number_t a2) {
  m.unset_all_implicit_hydrogen_information(a1);
  m.unset_all_implicit_hydrogen_information(a2);
}

// Fragments `frag` and `initial_comp` are the starting molecule.
// We are going to add `comp` to `frag` to make new molecules.
int
SetOfDatabases::MakeMolecules(PerMoleculeData& pmd,
                const Molecule& frag,
                Molecule& initial_comp,
                const Molecule& comp,
                int count,
                const AttachmentList& frag_attachments,
                const AttachmentList& comp_attachments,
                IWString_and_File_Descriptor& output) {
  if (frag_attachments.empty()) {
    return 0;
  }

  if (frag_attachments.size() > 3) {
    cerr << "SetOfDatabases::MakeMolecules: " << frag_attachments.size() <<
            " attachment points not implemented\n";
    return 0;
  }

// #define DEBUG_MAKE_MOLECULES
#ifdef DEBUG_MAKE_MOLECULES
  Molecule t1(frag);
  cerr << t1.smiles() << frag.name() << " FRAGMENT " << frag.natoms() << " starting atoms\n";
  cerr << initial_comp.smiles() << " STARTING COMPLEMENTARY\n";
  Molecule t2(comp);
  cerr << t2.aromatic_smiles() << " proposed replacement, count " << count << '\n';
  cerr << "frag attachments " << AttachmentsToString(frag_attachments) << '\n';
  cerr << "comp attachments " << AttachmentsToString(comp_attachments) << '\n';
#endif

  std::vector<AttachmentPairing> pairings;
  GenerateAttachmentPairings(frag_attachments, comp_attachments, pairings);
  if (pairings.empty()) {
    cerr << "SetOfDatabases::MakeMolecules:no isotope-compatible pairings\n";
    return 0;
  }

  int rc = 0;
  const int frag_natoms = frag.natoms();

  for (const AttachmentPairing& pairing : pairings) {
    Molecule product(frag);
    product.add_molecule(&comp);
    product << " %% " << comp.name() << ' ' << count;

#ifdef DEBUG_MAKE_MOLECULES
    cerr << product.aromatic_smiles() << " after adding\n";
    cerr << "have " << pairing.size() << " attachments\n";
#endif

    for (const AttachmentPair& pair : pairing) {
      atom_number_t comp_atom = frag_natoms + pair.comp_atom;
      product.add_bond(pair.frag_atom, comp_atom, SINGLE_BOND);
      UnfixImplicitHydrogens(product, pair.frag_atom, comp_atom);
    }

    if (product.number_fragments() > 1) {
      cerr << product.smiles() << ' ' << product.name() <<
              " multiple attachment multiple fragments\n";
    }

    rc += MaybeWriteProduct(pmd, product, count, initial_comp, output);
  }

  return rc;
}

// Return true if m.unique_smiles() is part of the _seen hash.
// Update the hash of not.
bool
SetOfDatabases::Seen(Molecule& m) {
  if (_seen.contains(m.unique_smiles())) {
    ++_duplicates_generated;
    return 1;
  }

  _seen.insert(m.unique_smiles());

  return 0;
}

int
SetOfDatabases::MaybeWriteProduct(PerMoleculeData& pmd,
                        Molecule& m,
                        int count,
                        Molecule& initial_comp,
                        IWString_and_File_Descriptor& output) {
  if (m.number_fragments() > 1) {
    cerr << m.smiles() << ' ' << m.name() << " multi fragment product\n";
  }
  // cerr << m.aromatic_smiles() << " MaybeWriteProduct\n";

  if (! OkProduct(m)) {
    return 0;
  }

  if (_remove_isotopes) {
    m.unset_isotopes();
    if (Seen(m)) {
      return 0;
    }
  } else {
    Molecule mcopy(m);
    mcopy.unset_isotopes();
    if (Seen(mcopy)) {
      return 0;
    }
  }

  if (! pmd.write_parent_molecule.empty()) {
    output << pmd.write_parent_molecule << '\n';
    pmd.write_parent_molecule.resize(0);
  }

  // Avoid writing unique smiles.
  m.invalidate_smiles();

  output << m.smiles() << ' ' << m.name();

  if (m.number_fragments() > 1) {
    cerr << m.smiles() << ' ' << m.name() << " multi fragment product\n";
  }

  if (_write_initial_complementary_fragment != '\0') {
    output << _write_initial_complementary_fragment << initial_comp.smiles();
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
DicerFragmentReplace(const dicer_data::DicedMolecule& proto,
                     SetOfDatabases& options,
                     IWString_and_File_Descriptor& output) {
  return options.Process(proto, output);
}

int
DicerFragmentReplaceLine(const const_IWSubstring& buffer,
                         SetOfDatabases& options,
                         IWString_and_File_Descriptor& output) {
  google::protobuf::io::ArrayInputStream zero_copy_array(buffer.data(), buffer.nchars());
  dicer_data::DicedMolecule proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "DicerFragmentReplaceLine:cannot parse proto " << buffer << '\n';
    return 0;
  }

  return DicerFragmentReplace(proto, options, output);
}

int
DicerFragmentReplace(iwstring_data_source& input,
                     SetOfDatabases& options,
                     IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    cerr << "processing " << buffer << '\n';
    if (! DicerFragmentReplaceLine(buffer, options, output)) {
      cerr << "DicerFragmentReplace:fatal error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
DicerFragmentReplace(const char* fname, 
                     SetOfDatabases& options,
                     IWString_and_File_Descriptor& output) {

  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DicerFragmentReplace:cannot open '" << fname << "'\n";
    return 0;
  }

  return DicerFragmentReplace(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:d:B:IC:x:c:m:p:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (! cl.option_present('d')) {
    cerr << "Must specify one or more databases via the -d option\n";
    Usage(1);
  }

  SetOfDatabases options;
  if (! options.Initialise(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
//  Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    cerr << "Processing '" << fname << "'\n";
    if (! DicerFragmentReplace(fname, options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_fragment_replace

int
main(int argc, char ** argv) {

  int rc = dicer_fragment_replace::Main(argc, argv);

  return rc;
}
