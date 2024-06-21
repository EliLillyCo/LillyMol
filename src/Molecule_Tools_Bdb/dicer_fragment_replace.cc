// Read texproto output from dicer.
// Use the output from dicer_to_topological_types as replacement fragments.
// For each fragment found in the new molecule, switch in plausible replacement
// fragments from previous diced collections to create new molecules.
// Note this needs to be combined with the tool of the same name in Molecule_Tools
// TODO: ianwatson implement ... 

#include <iostream>
#include <memory>
#include <optional>

#include "google/protobuf/text_format.h"

#include "db_cxx.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecular_formula.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_fragment_replace {

using std::cerr;

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

    molecular_formula::MolecularFormula<uint32_t> _mformula;

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

class SetOfDatabases {
  private:
    int _verbose;

    // If the input was generated with '-B nousmi' then we can match atoms without
    // needing to do any substructure searching.
    int _input_has_atom_maps;

    resizable_array_p<Db> _database;

    uint64_t _molecules_processed;

    Atom_Typing_Specification _atom_typing;

    extending_resizable_array<uint32_t> _variants_generated;

    // private functions
    int Process(Molecule& m,
                        const dicer_data::DicerFragment& fragment,
                        const Set_of_Atoms& embedding,
                        IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                        Molecule_to_Match& target,
                        const dicer_data::DicerFragment& fragment,
                        IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                        Molecule_to_Match& target,
                        const dicer_data::DicerFragment& fragment,
                        Substructure_Query& query,
                        resizable_array_p<dicer_data::FragmentAndComplements>& fromdb,
                        IWString_and_File_Descriptor& output);

  public:
    SetOfDatabases();

    int Initialise(Command_Line& cl);

    std::optional<int> Process(dicer_data::DicedMolecule& proto, IWString_and_File_Descriptor& output);
};

SetOfDatabases::SetOfDatabases() {
  _verbose = 0;
  _input_has_atom_maps = 0;
}

static void
DisplayDashBOptions() {
  cerr << R"( -B nousmi         input is from dicer with the -B nousmi option
)";

  ::exit(0);
}

int
SetOfDatabases::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');
  if (! cl.option_present('d')) {
    cerr << "Must specify one or more databases via the -d option\n";
    return 0;
  }

  int flags = 0;
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

  if (_verbose) {
    cerr << "Opened " << _database.size() << " dicer precedent databases\n";
  }

  return 1;
}

std::optional<int>
SetOfDatabases::Process(dicer_data::DicedMolecule& proto,
                        IWString_and_File_Descriptor& output) {
  Molecule m;
  if (! m.build_from_smiles(proto.smiles())) {
    cerr << "SetOfDatabases::Process:invalid smiles '" << proto.smiles() << "'\n";
    return std::nullopt;
  }

  m.set_name(proto.name());

  PerMoleculeData pmd(m);

  if (!_atom_typing.active()) {
    pmd.AssignAtomTypes(_atom_typing);
  }

  Molecule_to_Match target(&m);

  int rc = 0;
  for (const dicer_data::DicerFragment& fragment : proto.fragment()) {
    rc += Process(m, target, fragment, output);
  }

  return rc;
}

#ifdef NO_LONGER_NEEDED_QQ
int
FragmentFoundInDatabase(const IWString& smi, 
                        resizable_array_p<dicer_data::FragmentAndComplements>& retrieved) {
  Dbt dbkey, fromdb;
  dbkey.set_data((void*) smi.data());
  dbkey.set_size(smi.length());

  for (const Db* db : _databases) {
    if (0 != db->get(NULL, &dbkey, &fromdb, 0)) {
      continue;
    }
  }

  return 1;
}
#endif

int
SetOfDatabases::Process(Molecule& m,
                        Molecule_to_Match& target,
                        const dicer_data::DicerFragment& fragment,
                        IWString_and_File_Descriptor& output) {
  resizable_array_p<dicer_data::FragmentAndComplements> fromdb;
  Molecule frag;
  if (! frag.build_from_smiles(fragment.smi())){
    cerr << "SetOfDatabases::Process:invalid smiles " << fragment.ShortDebugString() << '\n';
    return 0;
  }

  Molecule f;
  if (! f.build_from_smiles(fragment.smi())) {
    cerr << "SetOfDatabases::Process:invalid fragment smiles " <<
            fragment.ShortDebugString() << '\n';
    return 0;
  }

  static Molecule_to_Query_Specifications mqs;
  mqs.set_substituents_only_at_isotopic_atoms(1);

  Substructure_Query qry;
  if (! qry.create_from_molecule(f, mqs)) {
    cerr << "SetOfDatabases::Process:cannot convert molecule to query " <<
            fragment.ShortDebugString() << '\n';
    return 0;
  }

  return Process(m, target, fragment, qry, fromdb, output);
}

int
SetOfDatabases::Process(Molecule& m,
                        Molecule_to_Match& target,
                        const dicer_data::DicerFragment& fragment,
                        Substructure_Query& query,
                        resizable_array_p<dicer_data::FragmentAndComplements>& fromdb,
                        IWString_and_File_Descriptor& output) {
  Substructure_Results sresults;
  if (query.substructure_search(target, sresults) == 0) {
    return 0;
  }

  int rc = 0;
  for (const Set_of_Atoms* e : sresults.embeddings()) {
    rc += Process(m, fragment, *e, output);
  }

  return rc;
}

int
SetOfDatabases::Process(Molecule& m,
                        const dicer_data::DicerFragment& fragment,
                        const Set_of_Atoms& embedding,
                        IWString_and_File_Descriptor& output) {
  return 1;
}

int
DicerFragmentReplace(iwstring_data_source& input,
                     SetOfDatabases& options,
                     IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    cerr << "finish implementation\n";
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
  Command_Line cl(argc, argv, "-vE:A:d:B:P:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

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
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    if (! DicerFragmentReplace(fname, options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  return 0;
}

}  // namespace dicer_fragment_replace

int
main(int argc, char ** argv) {

  int rc = dicer_fragment_replace::Main(argc, argv);

  return rc;
}
