// Use existing QupKake results to offer predictions for unknown molecules.
// First phase is to ingest QupKake data and build a database.
// When I originally formulated this I imagined doing an EC shell type
// implementation, but as I worked through that I decided that better
// would be to focus more on precise identification of the fragments
// of interest. Code to support a future EC fingerprint type implementation
// remains, but has not been tested at all and needs more work.

#include <algorithm>
#include <iostream>
#include <limits>

#include "absl/container/flat_hash_map.h"
#include "absl/hash/hash.h"
#include "absl/strings/string_view.h"

#include "db_cxx.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

// We use the same database Key as the iwecfp databases.
#include "iwecfp_database.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools_Bdb/qupkake_transfer.pb.h"
#else
#include "qupkake_transfer.pb.h"
#endif


namespace iwecfp_database {
// Make sure a DBKey struct can be hashed.
template <typename H>
H
AbslHashValue(H state, const iwecfp_database::DBKey& s) {
  return H::combine(std::move(state), (s._bit + s._acca + 127 * s._radius));
}

}  // namespace iwecfp_database

namespace qupkake_transfer {

using std::cerr;
using iwecfp_database::DBKey;

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
  cerr << R"(Transference of precomputed QupKake micro pKa results to new molecules.
Qupkake https://pubs.acs.org/doi/10.1021/acs.jctc.4c00328 
Could be used for any atomic property where the atomic property can be encoded as an isotopic value.
Currently being used for both Qupkake and EpiK pKa computations.

This tool performs two tasks:
  1. Read a series of Qupkake computations and store them in a database.
  2. Read an unknown molecule and lookup sites in a database from step 1.

 -d <dbname>    the name of the database.
 -d STORE       Build the database - the database file must not exist.
 -d LOOKUP      Lookup in the database - the database file must have been created by this tool.

The following options apply to building.
 -S <fname>     Generate a summary of all fragments stored in the database.

The following options apply to lookups.
 -s <smarts>    Only look for atoms matching <smarts>.
 -z i           Ignore molecules not matching any of the query atoms.
 -G <dbname>    The id's of computed fragments are in the Db. If you want to see structures, add
                one of more selimsteg databases providing id->smiles relationships.
 -G <dbname>,ZPAD=<n>    This selimsteg database is zero padded to <n> characters.
 -n <nsample>   Only write <nsample> values from the database - helps when many values are stored.
 -U <fname>     Write fragments not found to <fname> - useful for identifying needed computations.

The following options apply to both building and lookups.
 -m <natoms>    the minimum fragment size. If a fragment is too small, not enough context will be known.
 -c             remove chirality (recommended).
 -l             reduce to largest fragment (recommended).
 -g ...         chemical standardisation.
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

// A simplistic atom typing used for labeling attachment points.
enum CustomAtomType {
  kCarbonAromatic = 10,
  kCarbonSaturatedRing = 20,
  kCarbonSaturatedChain = 40,
  kCarbonUnsaturated = 60,

  kNitrogenAromatic = 80,
  kNitrogenSaturated = 90,
  kNitrogenUnSaturated = 110,

  kOxygenAromatic = 130,
  kOxygenSaturated = 140,
  kOxygenUnsaturated = 150,

  kFluorine = 160,
  kPhosphorus = 170,

  kSulphurAromatic = 180,
  kSulphurSaturated = 190,
  kSulpuhrUnsaturated = 210,

  KChlorine = 220,
  kBromine = 230,
  kIodine = 240,

  kOther = 250
};

// Number of Aryl connections to `zatom`
int
ArylCount(Molecule& m, atom_number_t zatom) {
  int rc = 0;
  for (const Bond* b : m[zatom]) {
    if (b->is_aromatic()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (m.is_aromatic(o)) {
      ++rc;
    }
  }

  return rc;
}

isotope_t
AtomType(Molecule& m, atom_number_t zatom) {
  const Atom& a = m[zatom];
  const atomic_number_t z = a.atomic_number();

  const int ahc = m.attached_heteroatom_count(zatom);
  const int aryl = ArylCount(m, zatom);
  const int acon = a.ncon() - 1;

  const int delta = 4 * acon + 2 * ahc * aryl;

  if (z == 6) {
    if (m.is_aromatic(zatom)) {
      return kCarbonAromatic + delta;
    }
    if (m.saturated(zatom)) {
      if (m.ring_bond_count(zatom)) {
        return kCarbonSaturatedRing + delta;
      } else {
        return kCarbonSaturatedChain + delta;
      }
    }
    return kCarbonUnsaturated + delta;
  }

  if (z == 7) {
    if (m.is_aromatic(zatom)) {
      return kNitrogenAromatic + delta;
    }
    if (m.saturated(zatom)) {
      return kNitrogenSaturated + delta;
    }
    return kNitrogenUnSaturated + delta;
  }

  if (z == 8) {
    if (m.is_aromatic(zatom)) {
      return kOxygenAromatic + delta;
    }
    if (m.saturated(zatom)) {
      return kOxygenSaturated + delta;
    }
    return kOxygenUnsaturated + delta;
  }

  if (z == 9) {
    return kFluorine + delta;
  }

  if (15 == z) {
    return kPhosphorus + delta;
  }

  if (16 == z) {
    if (m.is_aromatic(zatom)) {
      return kSulphurAromatic + delta;
    }
    if (m.saturated(zatom)) {
      return kSulphurSaturated + delta;
    }
    return kSulpuhrUnsaturated + delta;
  }

  if (z == 17) {
    return KChlorine + delta;
  }

  if (z == 35) {
    return kBromine + delta;
  }

  if (z == 53) {
    return kIodine + delta;
  }

  return kOther + delta;
}



struct PerMoleculeData {
  int* processing_status;
  atomic_number_t* z;
  int* aromatic;
  int* attached_heteroatom_count;
  int* saturated;
  // Used by whatever needs it.
  int* tmp;

  uint32_t* atype;

  IWString parent_smiles;

  PerMoleculeData();
  ~PerMoleculeData();

  int Build(Molecule& m);

  bool IsSaturatedCarbom(atom_number_t zatom) const;

  // If we have already allocated `atype`, just return it.
  // If not, allocate it and return the newly allocated array.
  // NOT initialised.
  uint32_t* get_atype(int n);
};

class MoleculeData;

// In order to avoid passing around a lot of arguments...
struct ShellData {
  atom_number_t centre_atom;
  const uint32_t* atype;
  int radius;
  int max_radius;
  uint32_t bit;

  ShellData();
};

enum class AcidBase {
  kAcidic,
  kBasic
};

// THis is constructed from a smiles that looks like
//   smiles id acid/basic
// And somewhere in the molecule is an isotope that is 1000 times the pka value.
class MoleculeData {
  private:
    Molecule _m;
    // The atom of interest will be the first atom with an isotope.
    atom_number_t _zatom;
    // these will be extracted from the name.
    float _pka;
    AcidBase _acid_base;

  // Private Functions

  public:
    MoleculeData();

    // during database builds the molecule name contains the experimental values.
    int BuildWithKnownValue(const Molecule& m);

    // during lookups, we don't have that extra information.
    int BuildNoKnownValue(const Molecule& m);

    Molecule& molecule() {
      return _m;
    }
    const Molecule& molecule() const {
      return _m;
    }
    atom_number_t centre_atom() const {
      return _zatom;
    }
    void set_centre_atom(atom_number_t s) {
      _zatom = s;
    }

    float value() const {
      return _pka;
    }

    bool IsAcid() const {
      return _acid_base == qupkake_transfer::AcidBase::kAcidic;
    }
};

class Storage {
  private:
    // Accumulate the data in memory before committing it at the end.
    absl::flat_hash_map<DBKey, qupkake_transfer::QupKakeData, absl::Hash<DBKey>> _data;
    absl::flat_hash_map<IWString, qupkake_transfer::QupKakeData> _usmi_to_data;

    std::unique_ptr<Db> _database;

  public:
    Storage();
    ~Storage();

    int OpenForWriting(IWString& dbname);
    int OpenForReading(IWString& dbname);

    // Write the contents of _data to _database.
    int WriteToDatabase();

    // Write the hash contents to a file.
    int WriteHash(IWString_and_File_Descriptor& output);

    int StoreUsmi(Molecule& m, const MoleculeData& moldata);

    std::optional<qupkake_transfer::QupKakeData> Lookup(const IWString& usmi);

    int StoreShell(Molecule& m, const MoleculeData& moldata, PerMoleculeData& pmd, int radius, uint32_t sum_so_far);
};

Storage::Storage() {
  _database.reset(new Db(NULL, DB_CXX_NO_EXCEPTIONS));
}

Storage::~Storage() {
  _database->close(0);
}

int
OpenDbForReading(IWString& dbname, Db* db) {
  if (! dash_s(dbname)) {
    cerr << "OpenForReading:db not found '" << dbname << "'\n";
    return 0;
  }

  int flags = DB_RDONLY;
  DBTYPE dbtype = DB_UNKNOWN;
  int mode = 0;

  int rc = db->open(NULL, dbname, NULL, dbtype, flags, mode);

  if (rc == 0) {
    return 1;
  }

  cerr << "OpenForReading:cannot open database '" << dbname << "'\n";
  db->err(rc, "");
  return 0;
}

int
Storage::OpenForReading(IWString& dbname) {
  return OpenDbForReading(dbname, _database.get());

  if (! dash_s(dbname)) {
    cerr << "Storage::OpenForReading:db not found '" << dbname << "'\n";
    return 0;
  }

  int flags = 0;
  DBTYPE dbtype = DB_UNKNOWN;
  int mode = 0;

  int rc = _database->open(NULL, dbname, NULL, dbtype, flags, mode);

  if (rc == 0) {
    return 1;
  }

  cerr << "Storage::OpenForReading:cannot open database '" << dbname << "'\n";
  _database->err(rc, "");
  return 0;
}

int
Storage::OpenForWriting(IWString& dbname) {
  if (dash_s(dbname)) {
    cerr << "Storage::OpenForWriting:db already exists '" << dbname << "'\n";
    return 0;
  }

  int flags = DB_CREATE;
  DBTYPE dbtype = DB_BTREE;
  int mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

  int rc = _database->open(NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);
  if (rc == 0) {
    return 1;
  }
  cerr << "Storage::OpenForWriting:cannot open database '" << dbname << "'\n";
  _database->err(rc, "");
  return 0;
}

// Add another ComputedValue to `proto`.
int
AddExtraValue(const MoleculeData& moldata, qupkake_transfer::QupKakeData& proto) {
  IWString name = moldata.molecule().name();
  name.truncate_at_first(' ');

  if (moldata.IsAcid()) {
    qupkake_transfer::ParentValue* pv = proto.mutable_acid()->mutable_values()->Add();
    pv->set_par(name.data(), name.length());
    pv->set_value(moldata.value());
  } else {
    qupkake_transfer::ParentValue* pv = proto.mutable_base()->mutable_values()->Add();
    pv->set_par(name.data(), name.length());
    pv->set_value(moldata.value());
  }

  return 1;
}

int
Storage::StoreShell(Molecule& m, const MoleculeData& moldata, PerMoleculeData& pmd, int radius, uint32_t sum_so_far) {

  cerr << "StoreShell: radius " << radius << " sum_so_far " << sum_so_far << '\n';
  iwecfp_database::DBKey dbkey;
  iwecfp_database::form_key(sum_so_far, radius, pmd.atype[moldata.centre_atom()], dbkey);

  auto iter = _data.find(dbkey);
  if (iter != _data.end()) {
    return AddExtraValue(moldata, iter->second);
  }

  // Build new proto and store.

  const int matoms = m.natoms();

  const int* processing_status = pmd.processing_status;
  int* tmp = pmd.tmp;
  for (int i = 0; i < matoms; ++i) {
    if (processing_status[i] == 0) {
      tmp[i] = 0;
    } else {
      tmp[i] = 1;
    }
  }

  Molecule subset;

  m.create_subset(subset, pmd.tmp, 1);

  const IWString& usmi = subset.unique_smiles();

  qupkake_transfer::QupKakeData proto;
  proto.set_usmi(usmi.data(), usmi.size());

  AddExtraValue(moldata, proto);

  _data.emplace(dbkey, std::move(proto));

  return 1;
}

int
Storage::StoreUsmi(Molecule& m, const MoleculeData& moldata) {
  const IWString& usmi = m.unique_smiles();

  auto iter = _usmi_to_data.find(usmi);
  if (iter != _usmi_to_data.end()) {
    AddExtraValue(moldata, iter->second);
    return 1;
  }

  qupkake_transfer::QupKakeData proto;
  proto.set_usmi(usmi.data(), usmi.length());
  AddExtraValue(moldata, proto);

  _usmi_to_data[usmi] = std::move(proto);

  // cerr << usmi << " added now have " << _usmi_to_data.size() << " relationships\n";

  return 1;
}

int
Store(const DBKey& dbkey, const qupkake_transfer::QupKakeData& proto, Db& db) {
  std::string serialized;
  proto.SerializeToString(&serialized);

  Dbt key;
  key.set_data((void*) &dbkey);
  key.set_size(sizeof(dbkey));

  Dbt to_store;

  to_store.set_data(serialized.data());
  to_store.set_size(serialized.size());

  if (db.put(NULL, &key, &to_store, 0) == 0) {
    return 1;
  }

  cerr << "Cannot store " << proto.ShortDebugString() << '\n';
  return 0;
}

int
Store(const IWString& usmi, const qupkake_transfer::QupKakeData& proto, Db& db) {
  std::string serialized;
  proto.SerializeToString(&serialized);

  Dbt key;
  key.set_data((void*)usmi.data());
  key.set_size(usmi.length());

  Dbt to_store;

  to_store.set_data(serialized.data());
  to_store.set_size(serialized.size());

  if (db.put(NULL, &key, &to_store, 0) == 0) {
    return 1;
  }

  cerr << "Cannot store " << proto.ShortDebugString() << '\n';
  return 0;
}

void
CalculateSummaryStats(qupkake_transfer::ComputedValues& proto) {
  Accumulator<double> acc;

  for (const auto& value : proto.values()) {
    acc.extra(value.value());
  }
  proto.set_minval(acc.minval());
  proto.set_maxval(acc.maxval());
  proto.set_mean(acc.average());
}

int
Storage::WriteToDatabase() {
  cerr << "Storage::WriteToDatabase:writing " << _usmi_to_data.size() << " items to db\n";
  for (auto& [k, v] : _usmi_to_data) {
    CalculateSummaryStats(*v.mutable_acid());
    CalculateSummaryStats(*v.mutable_base());
    Store(k, v, *_database);
  }

  return 1;
}

struct ValueAndId {
  float value;
  IWString id;
};

// Write a sorted list of all the ComputedValue items in `proto`.
// Only write computations of type `target`.
// `name` will likely be one of 'acid' or 'base'
// `value_id` is a temporary array we use.
int
WriteAcidBase(const qupkake_transfer::ComputedValues& proto,
              const char* name,
              ValueAndId* value_id,
              IWString_and_File_Descriptor& output) {
  Accumulator<double> acc;
  resizable_array<float> values;
  values.reserve(proto.values_size());

  int ndx = 0;
  for (const auto& value : proto.values()) {
    acc.extra(value.value());
    value_id[ndx].value = value.value();
    value_id[ndx].id = value.par();
    values << value.value();
    ++ndx;
  }

  if (ndx == 0) {
    return 0;
  }

  if (ndx > 1) {
    std::sort(value_id, value_id + ndx, [](const ValueAndId& vid1, const ValueAndId& vid2) {
      return vid1.value < vid2.value;
    });
    values.iwqsort_lambda([](float v1, float v2) { return v1 < v2;});
  }

  static constexpr char kSep = ' ';

  for (int i = 0; i < ndx; ++i) {
    output << value_id[i].id << kSep << name << kSep << value_id[i].value << '\n';
  }

  output << name << kSep << acc.n() << " values btw " << static_cast<float>(acc.minval()) << " ave " << 
            static_cast<float>(acc.average()) << " " << static_cast<float>(acc.maxval()) << 
            " range " << static_cast<float>(acc.range()) << " median ";
  if (ndx == 1) {
    output << values[0];
  } else if (ndx == 2) {
    output << static_cast<float>(acc.average());
  } else if (ndx % 2 == 0) {
    int i1 = ndx / 2 - 1;
    int i2 = ndx / 2;
    output << ((values[i1] + values[i2]) * 0.5f);
  } else {
    output << values[ndx / 2];
  }
  output << '\n';

  return 1;
}

int
WriteHashMember(const IWString& usmi, const qupkake_transfer::QupKakeData& proto,
          IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';  // not fully implemented...
  output << usmi << kSep << "acid" << kSep << proto.acid().values_size() << kSep << 
            "base" << kSep << proto.base().values_size() << '\n';

  int n = std::max(proto.acid().values_size(), proto.base().values_size());

  std::unique_ptr<ValueAndId[]> value_id = std::make_unique<ValueAndId[]>(n);
  WriteAcidBase(proto.acid(), "acid", value_id.get(), output);

  output.write_if_buffer_holds_more_than(4096);

  WriteAcidBase(proto.base(), "base", value_id.get(), output);

  output.write_if_buffer_holds_more_than(4096);
  return 1;
}

int
WriteHashMember(const DBKey& key, const qupkake_transfer::QupKakeData& proto,
          IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';  // not fully implemented...

  output << key._radius << kSep << key._acca << kSep << key._bit << '\n';
  output << proto.usmi() << '\n';

  int n = std::max(proto.acid().values_size(), proto.base().values_size());

  std::unique_ptr<ValueAndId[]> value_id = std::make_unique<ValueAndId[]>(n);

  WriteAcidBase(proto.acid(), "acid", value_id.get(), output);

  output.write_if_buffer_holds_more_than(4096);

  WriteAcidBase(proto.base(), "base", value_id.get(), output);

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Storage::WriteHash(IWString_and_File_Descriptor& output) {
  for (const auto& [k, v] : _usmi_to_data) {
    WriteHashMember(k, v, output);
  }

//for (const auto& [k, v] : _data) {
//  WriteHashMember(k, v, output);
//}

  return 1;
}

std::optional<qupkake_transfer::QupKakeData>
Storage::Lookup(const IWString& usmi) {
  Dbt key;
  key.set_data((void*) usmi.data());
  key.set_size(usmi.length());

  Dbt fromdb;

  if (_database->get(NULL, &key, &fromdb, 0) != 0) {
    return std::nullopt;
  }

  const absl::string_view sv((const char*) fromdb.get_data(), fromdb.get_size());
  qupkake_transfer::QupKakeData proto;
  proto.ParseFromString(sv);

  return std::move(proto);
}

PerMoleculeData::PerMoleculeData() {
  processing_status = nullptr;
  z = nullptr;
  aromatic = nullptr;
  attached_heteroatom_count = nullptr;
  saturated = nullptr;
  atype = nullptr;

  tmp = nullptr;
}

PerMoleculeData::~PerMoleculeData() {
  delete [] processing_status;
  delete [] z;
  delete [] aromatic;
  delete [] attached_heteroatom_count;
  delete [] saturated;
  delete [] atype;

  delete [] tmp;
}

int
PerMoleculeData::Build(Molecule& m) {
  const int matoms = m.natoms();
  processing_status = new int[matoms];
  z = new atomic_number_t[matoms];
  aromatic = new int[matoms];
  attached_heteroatom_count = new int[matoms];
  saturated = new int[matoms];
  tmp = new int[matoms];

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < matoms; ++i) {
    z[i] = m.atomic_number(i);
    aromatic[i] = m.is_aromatic(i);
    attached_heteroatom_count[i] = m.attached_heteroatom_count(i);
    saturated[i] = m.saturated(i);
  }

  parent_smiles = m.smiles();

  return 1;
}

// Return true if `zatom` is a fully saturated carbon atom.
bool
PerMoleculeData::IsSaturatedCarbom(atom_number_t zatom) const {
  if (z[zatom] != 6) {
    return 0;
  }

  return saturated[zatom];
}

uint32_t* 
PerMoleculeData::get_atype(int n) {
  if (atype == nullptr) {
    atype = new uint32_t[n];
  }

  return atype;
}

MoleculeData::MoleculeData() {
  _zatom = kInvalidAtomNumber;
  _pka = 0.0;
}

// Any value less than 1000 (1.0) will be ignored. I think that is OK.
constexpr float kDivide = 1000.0f;
constexpr isotope_t kDefaultIsotope = 9;


int
MoleculeData::BuildWithKnownValue(const Molecule& m) {
  // In order to make unique smiles match, need to impose a single
  // isotope at each centre atom - once we have processed the incoming value.
  _m = m;  // make out own copy.
  const int matoms = _m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso = _m.isotope(i);
    if (iso == 0) {
      continue;
    }

    if (iso <= kDivide) {  // previous test is redundant, kept for clarity.
      continue;
    }

    _pka = static_cast<float>(iso) / kDivide;
    _zatom = i;
    break;
  }

  if (_zatom == kInvalidAtomNumber) {
    cerr << "Data::BuildWithKnownValue:no isotope\n";
    cerr << _m.smiles() << ' ' << _m.name() << '\n';
    return 0;
  }

  // the name contains the assignment type.
  // 'id basic' or 'id acidic'

  const IWString& name = m.name();
  if (name.empty()) {
    cerr << "MoleculeData::BuildWithKnownValue:name is empty\n";
    return 0;
  }
  int i = 0;
  const_IWSubstring token;
  name.nextword(token, i);
  if (! name.nextword(token, i)) {
    cerr << "MoleculeData::BuildWithKnownValue:invalid name '" << name << "'\n";
    return 0;
  }

  if (token == "acidic") {
    _acid_base = AcidBase::kAcidic;
  } else if (token == "basic") {
    _acid_base = AcidBase::kBasic;
  } else {
    cerr << "MoleculeData::BuildWithKnownValue:unknown acidic/basic '" << name << "'\n";
    return 0;
  }

  return 1;
}

int
MoleculeData::BuildNoKnownValue(const Molecule& m) {
  _m = m;

  return 1;
}

uint32_t
BondConstant(const Bond* bondi) {
//if (all_bonds_same_type) {
//  return 1;
//}

  if (bondi->is_aromatic()) {
    return 11;
  }
  if (bondi->is_triple_bond()) {
    return 7;
  }
  if (bondi->is_double_bond()) {
    return 5;
  }

  return 3;
}

int additive = 1;

void
Increment(const int bc, const uint32_t atom_constant, uint32_t& sum_so_far) {
  // cerr << "Before increment " << sum_so_far << " bc " << bc << " atom_constant " <<
  // atom_constant << '\n';

  if (additive) {
    sum_so_far += bc * atom_constant;
  } else {
    sum_so_far *= bc * atom_constant;
  }

  return;
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    // We can either build a database or do lookups.
    int _build_database;

    int _reduce_to_largest_fragment = 0;

    extending_resizable_array<uint32_t> _atoms_in_subset;
    extending_resizable_array<uint32_t> _atoms_lost;

    Atom_Typing_Specification _atom_typing_specification;

    Storage _storage;

    // The minimum radius is used when forming fragments.
    // In chains, if an atom is within _min_radius of the atom of
    // interest, it will be included in the subset.
    int _min_radius;

    int _max_radius;

    // the -d option.
    IWString _dbname;

    // Find that radius 1 shells are not very specific and have lots
    // of collisions. Optional as to whether they are stored or not.
    int _store_radius_1_shells;

    // Once we identify the fragment if interest, we can force it to be
    // larger if we want.
    int _min_fragment_size;

    // If a large number of fragments are in the database, we can overwhelm the user
    // with data. Impose a limit on the number of calculated values displayed.
    int _nwrite;

    // For each molecule, we keep track of the number times we find a match
    // to a fragment.
    extending_resizable_array<int> _fragments_found;

    int _remove_chirality = 0;

    // during lookups, the default is to lookup all N and O atoms.
    // That search can be limited via substructure queries.
    resizable_array_p<Substructure_Query> _queries;
    int _ignore_molecules_not_matching_queries;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    uint64_t _molecules_read = 0;

    // Only identifiers are stored in the database. It can be convenient
    // to be able to look up smiles when presenting results.
    resizable_array_p<Db> _selimsteg;
    // Some databases might be stored with leading zero's.
    std::vector<int> _zero_pad_identifiers;

    // The -U option.
    IWString_and_File_Descriptor _stream_for_fragments_not_found;

    // Private functions

    int BuildDatabase(Molecule& m,
                 MoleculeData& mdata,
                 PerMoleculeData& pmd);
    int DoLookups(Molecule& m,
                 MoleculeData& mdata,
                 PerMoleculeData& pmd,
                 IWString_and_File_Descriptor& output);
    int DoLookups(Molecule& m,
                   atom_number_t zatom,
                   IWString_and_File_Descriptor& output);
    int IdentifyAtomsToProcess(Molecule& m, int* process_atom);
    int LookupViaUniqueSmiles(Molecule& parent, Molecule& m, atom_number_t zatom,
                IWString_and_File_Descriptor& output);
    int ReduceToFragmentOfInterest(Molecule& m, MoleculeData& mdata, atom_number_t zatom);
    int PropagateThroughChains(Molecule& m, atom_number_t anchor,
                atom_number_t previous_atom, atom_number_t zatom, 
                isotope_t* iso, int* tmp);
    int MaybeInitialiseAtomTypes(Molecule& m, PerMoleculeData& pmd);
    int WriteComputedValue(const qupkake_transfer::ParentValue& c,
                            const char* acid_base,
                            IWString_and_File_Descriptor& output);
    int IdentifyShells(Molecule& m,
                        MoleculeData& mdata,
                        PerMoleculeData& pmd,
                        int radius,
                        uint32_t& sum_so_far);
    int MaybeWriteUnknownFragment(const Molecule & parent, Molecule& m);
    int ComputeIstep(int nvalues) const;

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    int build_database() const {
      return _build_database;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int BuildDatabase(Molecule& m);

    // Write the contents of `_data` to `_dbname`.
    int WriteToDatabase();
    // When doing database builds, write the hash to a stream.
    int WriteHash(IWString& fname);

    int DoLookups(Molecule& m, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _build_database = 0;
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _min_radius = 3;
  _max_radius = 3;
  _store_radius_1_shells = 1;
  _ignore_molecules_not_matching_queries = 0;
  _min_fragment_size = 0;
  _nwrite = std::numeric_limits<int32_t>::max();

  // We do a bunch of copies.
  set_copy_name_in_molecule_copy_constructor(1);
}

void
DisplayDashMOptions(std::ostream& output) {
  output << " -M no1            do NOT store radius 1 bits\n";
  exit(0);
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
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

  if (cl.option_present('P')) {
    const IWString p = cl.string_value('P');
    if (! _atom_typing_specification.build(p)) {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 0;
    }
  } else {
    _atom_typing_specification.build("UST:AY");
  }

  if (! cl.option_present('d')) {
    cerr << "Must specify name of database to read/write via the -d option\n";
    Usage(1);
  }

  IWString d;
  for(int i = 0; cl.value('d', d, i); ++i) {
    if (d == "STORE") {
      _build_database = 1;
    } else if (d == "LOOKUP") {
      _build_database = 0;
    } else {
      _dbname = d;
    }
  }

  if (_dbname.empty()) {
    cerr << "Must specify name of database to read/write via the -d option\n";
    Usage(1);
  }

  if (_build_database) {
    _dbname.EnsureEndsWith(".bdb");
  }

  if (_build_database && dash_s(_dbname.null_terminated_chars())) {
    cerr << "Options::Initialise:building a database, cannot write to existing database '" << _dbname << "'\n";
    return 0;
  }

  if (_build_database) {
    if (! _storage.OpenForWriting(_dbname)) {
      cerr << "Options::Initialise:cannot open '" << _dbname << "' for writing\n";
      return 0;
    }
  } else {
    if (! _storage.OpenForReading(_dbname)) {
      cerr << "Options::Initialise:cannot open '" << _dbname << "' for reading\n";
      return 0;
    }
  }

  if (cl.option_present('M')) {
    IWString m;
    for (int i = 0; cl.value('M', m, i); ++i) {
      if (m == "no1") {
        _store_radius_1_shells = 0;
        if (_verbose) {
          cerr << "Radius 1 shells will NOT be stored\n";
        }
      } else if (m == "help") {
        DisplayDashMOptions(cerr);
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        DisplayDashMOptions(cerr);
      }
    }
  }

  if (cl.option_present('s')) {
    IWString smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Options::Initialise:invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _queries << qry.release();
    }

    if (_verbose) {
      cerr << "Defined " << _queries.size() << " queries for matched atoms\n";
    }

    if (cl.option_present('z')) {
      _ignore_molecules_not_matching_queries = 1;
      if (_verbose) {
        cerr << "Will ignore molecules not matching any queries\n";
      }
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_fragment_size) || _min_fragment_size < 1) {
      cerr << "The min fragment size (-m) option must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Minimum fragment size " << _min_fragment_size << '\n';
    }
  }

  if (cl.option_present('G')) {
    IWString goption;
    for (int i = 0; cl.value('G', goption, i); ++i) {
      IWString dbname, zpad_directive;
      goption.split(dbname, ',', zpad_directive);  // zpad might be empty.

      int zpad = 0;
      if (zpad_directive.size() > 0) {
        zpad_directive.remove_leading_chars(5);
        if (! zpad_directive.numeric_value(zpad) || zpad <  6) {
          cerr << "The zero pad identifiers option must be a whole +ve number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Database '" << dbname << " zero padded to " << zpad << " width before selimisteg lookup\n";
        }
      }

      std::unique_ptr<Db> db(new Db(NULL, DB_CXX_NO_EXCEPTIONS));
      if (! OpenDbForReading(dbname, db.get())) {
        cerr << "Options::Initialise:cannot open selimsteg db '" << dbname << "'\n";
        return 0;
      }
      _selimsteg << db.release();
      _zero_pad_identifiers.push_back(zpad);
    }

    if (_verbose) {
      cerr << "Defined " << _selimsteg.size() << " selimsteg databases\n";
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _nwrite) || _nwrite < 1) {
      cerr << "The number of database entries to write (-n) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only write " << _nwrite << " entries from the database\n";
    }
  }

  if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    fname.EnsureEndsWith(".smi");
    if (! _stream_for_fragments_not_found.open(fname)) {
      cerr << "Options::Initialise:cannot open -U file '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Fragments not found in db written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";

  for (int i = 0; i < _atoms_lost.number_elements(); ++i) {
    if (_atoms_lost[i]) {
      output << _atoms_lost[i] << " molecule lost " << i << " atoms in making fragment of interest\n";
    }
  }

  for (int i = 0; i < _fragments_found.number_elements(); ++i) {
    if (_fragments_found[i]) {
      output << _fragments_found[i] << " molecules found " << i << " fragments\n";
    }
  }

  for (int i = 0; i < _atoms_in_subset.number_elements(); ++i) {
    if (_atoms_in_subset[i]) {
      output << _atoms_in_subset[i] << " fragments had " << i << " atoms\n";
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

// Exocyclic double bonds not handled.
void
PropagateThroughRingSystem(Molecule& m, atom_number_t zatom, int* tmp) {
  // cerr << "PropagateThroughRingSystem atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << '\n';
  tmp[zatom] = 1;
  for (const Bond* b : m[zatom]) {
    if (b->nrings() == 0) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (tmp[o]) {
      continue;
    }
    if (m.ring_bond_count(o)) {
      PropagateThroughRingSystem(m, o, tmp);
    }
  }
}

void
PropagateThroughChains(Molecule& m, atom_number_t zatom, int* tmp) {
  tmp[zatom] = 1;
  for (const Bond* b : m[zatom]) {
    if (b->nrings() > 0) [[unlikely]] {  // not sure it can happen at all.
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (tmp[o]) {
      continue;
    }
    if (m.ring_bond_count(o) == 0) {
      PropagateThroughChains(m, o, tmp);
    }
  }
}

int
Options::PropagateThroughChains(Molecule& m, atom_number_t anchor,
                atom_number_t previous_atom, atom_number_t zatom, 
                isotope_t* iso, int* tmp) {
  tmp[zatom] = 1;
  for (const Bond* b : m[zatom]) {
    if (b->nrings() > 0) {  // can happen if we are expanding from a substituted ring atom.
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (tmp[o]) {
      continue;
    }

    if (! b->is_single_bond()) {
      PropagateThroughChains(m, anchor, zatom, o, iso, tmp);
      continue;
    }

    if (m.ring_bond_count(o) > 0) {
      if (m.bonds_between(anchor, o) <= _min_radius) {
        PropagateThroughRingSystem(m, o, tmp);
      } else {
        tmp[o] = 1;
        iso[zatom] = AtomType(m, o);
      }
      continue;
    }

    atomic_number_t z = m.atomic_number(zatom);
    if (z != 6) {
      PropagateThroughChains(m, anchor, zatom, o, iso, tmp);
      continue;
    }

    if (z != m.atomic_number(previous_atom)) {
      PropagateThroughChains(m, anchor, zatom, o, iso, tmp);
      continue;
    }

    // We are left with a C-C single bond.
    if (m.bonds_between(anchor, o) <= _min_radius) {
      PropagateThroughChains(m, anchor, zatom, o, iso, tmp);
    } else {
      tmp[o] = 1;
      iso[o] = AtomType(m, o);
    }
  }

  return 1;
}

// If we can find a singly bonded neighbour of `zatom` that is set in `tmp`, return it.
std::optional<atom_number_t>
BondedToTmp(const Molecule& m, atom_number_t zatom, const int* tmp) {
  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (tmp[o]) {
      return o;
    }
  }

  return std::nullopt;
}

// `tmp` holds the ids of the chain atoms which include `zatom`. 
// zatom is non ring.
// Add rings that connect to the marked atoms, and remove everything else.
void
RemoveNonAdjacent(Molecule& m, atom_number_t zatom, int* tmp) {
  if (m.nrings() == 0) {
    return;
  }

  const int matoms = m.natoms();
  // int nonzero = std::count_if(tmp, tmp + matoms, [](int x) { return x == 1;});
  // cerr << "Start with " << nonzero << " nonzero atoms\n";

  // Find ring atoms adjacent to an atom in `tmp` and propagate through that ring system.
  for (int i = 0; i < matoms; ++i) {
    if (tmp[i] == 0) {
      continue;
    }
    for (const Bond* b : m[i]) {
      atom_number_t o = b->other(i);
      if (tmp[o]) {
        continue;
      }
      if (m.ring_bond_count(o)) {
        PropagateThroughRingSystem(m, o, tmp);
      }
    }
  }

  // nonzero = std::count_if(tmp, tmp + matoms, [](int x) { return x == 1;});
  //cerr << "before adding doubly bonded " << nonzero << " nonzero atoms\n";
  //cerr << "Attaching doubly bonded atoms\n";
  // Attach any atoms doubly bonded to an included atom.
  for (int i = 0; i < matoms; ++i) {
    if (tmp[i]) {
      continue;
    }
    const Atom& a = m[i];
    if (a.ncon() != 1) {
      continue;
    }
    const Bond* b = a[0];
    if (! b->is_double_bond()) {
      continue;
    }
    if (a.atomic_number() == 7 || a.atomic_number() == 8) {
    } else {
      continue;
    }
    atom_number_t o = b->other(i);
    if (tmp[o]) {  // If `o` is in the system add `i`
      tmp[i] = 1;
    }
  }

  for (int i = matoms - 1; i >= 0; --i) {
    if (tmp[i] == 0) {
      m.remove_atom(i);
    }
  }
}
// `tmp` holds the identify of the atoms in the ring system holding `zatom`.
// Remove all atoms not bonded to an atom in that system
void
RemoveBeyondRingOfInterest(Molecule& m, atom_number_t zatom, const int* tmp) {
  Set_of_Atoms to_remove;
  // Make sure atoms are added largest atom number first.
  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (tmp[i]) {
      continue;
    }
    std::optional<atom_number_t> bonded_to_tmp = BondedToTmp(m, i, tmp);
    if (bonded_to_tmp) {
      m.set_isotope(*bonded_to_tmp, AtomType(m, i));
    } else {
      to_remove << i;
    }
  }

  for (atom_number_t a : to_remove) {
    m.remove_atom(a);
  }
}

// Identify any unmatched heteratoms atoms in `tmp` that are doubly bonded to
// an atom in `tmp` and include those atoms.
int
IncludeAdjacentDoublyBonded(Molecule& m, int* tmp) {
  int rc = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (tmp[i]) {
      continue;
    }
    if (m.atomic_number(i) == 6) {
      continue;
    }

    // No allowance for =S groups - who cares...
    if (m.atomic_number(i) > 8) {
      continue;
    }

    if (m.ring_bond_count(i)) {
      continue;
    }

    bool bonded_to_tmp = false;
    for (const Bond* b : m[i]) {
//    Turns out that for many ring systems, external heteratoms are important.
//    if (b->is_single_bond()) {
//      continue;
//    }
      if (b->is_aromatic()) {
        continue;
      }
      atom_number_t o = b->other(i);
      if (tmp[o]) {
        bonded_to_tmp = true;
        break;
      }
    }

    if (! bonded_to_tmp) {
      continue;
    }

    tmp[i] = 1;
    ++rc;
  }

  return rc;
}

// Examine bonds where one side is in the system and the other is not.
int
IncludeInfoAboutJoined(Molecule& m,
                       isotope_t* iso,
                       int* tmp) {

  Set_of_Atoms added_here;

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (tmp[a1] && tmp[a2]) {
      continue;
    }
    if (tmp[a1] == 0 && tmp[a2] == 0) {
      continue;
    }

    // Make sure that `a1` is in the system.
    if (tmp[a1] == 0) {
      std::swap(a1, a2);
    }

    // A biphenyl like connection.
    if (m.is_aromatic(a1) && m.is_aromatic(a2)) {
      iso[a1] = AtomType(m, a2);
      continue;
    }

    // Add all doubly bonded neighbours.
    if (b->is_double_bond() && b->nrings() == 0) {
      added_here << a2;
      iso[a2] = AtomType(m, a2);
      continue;
    }

    if (m.is_aromatic(a1)) {
      iso[a1] = AtomType(m, a2);
      continue;
    } else if (m.is_aromatic(a2)) {
      iso[a1] = AtomType(m, a2);
      continue;
    }

    if (b->is_double_bond()) {
      added_here << a2;
      iso[a2] = AtomType(m, a2);
      continue;
    }
    // It is possible this could include a single atom of a ring. OK.
    if (m.atomic_number(a2) != 6) {
      added_here << a2;
      iso[a2] = AtomType(m, a2);
      continue;
    }
    if (m.isotope(a1) == 0) {
      iso[a1] = AtomType(m, a2);
    }
  }

  added_here.set_vector(tmp, 1);

  return 1;
}
int
AddHeteratomsToRingAtoms(Molecule& m, atom_number_t ignore, isotope_t* iso, int* tmp) {
  Set_of_Atoms added_here;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (tmp[i] == 0) {
      continue;
    }
    if (! m.is_aromatic(i)) {
      continue;
    }
    const Atom& a = m[i];
    if (a.ncon() == 2) {
      continue;
    }

    bool isanchor = a.isotope() == kDefaultIsotope;

    for (const Bond* b : a) {
      if (b->is_aromatic()) {
        continue;
      }
      atom_number_t o = b->other(i);
      if (tmp[o] || o == ignore) {
        continue;
      }
      if (isanchor) {
        added_here << o;
      } else if (m.atomic_number(o) == 6) {
        iso[i] = AtomType(m, o);
      } else {
        added_here << o;
      }
    }
  }

  added_here.set_vector(tmp, 1);

  return 1;
}

// We have formed the fragment, but it does not have enough atoms.
// Add atoms that are ajdacent to the matched atoms until we have `extra_atoms_needed`
// atoms added.
// Note that this may be impossible, since we only look at the first layer of
// attached atoms.
void
ExpandBeyondCurrent(Molecule& m, int* in_fragment,
                    int extra_atoms_needed) {
  const int matoms = m.natoms();

  // Gather the types of atoms that are adjacent to already included atoms.
  Set_of_Atoms heteroatoms, aromatic, unsaturated, other;

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (in_fragment[i] == 0) {
      continue;
    }

    for (const Bond* b : m[i]) {
      atom_number_t o = b->other(i);
      if (in_fragment[o]) {
        continue;
      }

      if (m.is_aromatic(o)) {
        aromatic << o;
      } else if (! m.saturated(o)) {
        unsaturated << o;
      } else if (m.atomic_number(o) != 6) {
        heteroatoms << o;
      } else {
        other << o;
      }
    }
  }

  for (atom_number_t a : aromatic) {
    in_fragment[a] = 1;
  }

  rc += aromatic.number_elements();

  if (rc >= extra_atoms_needed) {
    return;
  }

  for (atom_number_t a : unsaturated) {
    in_fragment[a] = 1;
  }

  rc += unsaturated.number_elements();

  if (rc >= extra_atoms_needed) {
    return;
  }

  for (atom_number_t a : heteroatoms) {
    in_fragment[a] = 1;
  }

  rc += heteroatoms.number_elements();

  if (rc >= extra_atoms_needed) {
    return;
  }

  for (atom_number_t a : other) {
    in_fragment[a] = 1;
  }

  rc += other.number_elements();

  if (rc >= extra_atoms_needed) {
    return;
  }

  return;
}


// `zatom` is an atom of interest. Remove atoms from `m` so we are left with just
// relevant atoms joined to `zatom`.
int
Options::ReduceToFragmentOfInterest(Molecule& m, MoleculeData& mdata, atom_number_t zatom) {
  const int initial_matoms = m.natoms();

  std::unique_ptr<int[]> tmp = std::make_unique<int[]>(initial_matoms);
  std::fill_n(tmp.get(), initial_matoms, 0);
  std::unique_ptr<isotope_t[]> iso = std::make_unique<isotope_t[]>(initial_matoms);
  std::fill_n(iso.get(), initial_matoms, 0);

  m.compute_aromaticity_if_needed();

  // cerr << "Starting molecule\n";
  // cerr << m.smiles() << ' ' << m.name() << '\n';
  if (m.ring_bond_count(zatom)) {
    PropagateThroughRingSystem(m, zatom, tmp.get());
    AddHeteratomsToRingAtoms(m, zatom, iso.get(), tmp.get());
    // If substitued, expand through the chain.
    if (m.ncon(zatom) > 2 && m.ring_bond_count(zatom) == 2) {
      PropagateThroughChains(m, zatom, kInvalidAtomNumber, zatom, iso.get(), tmp.get());
    }
  } else {
    PropagateThroughChains(m, zatom, kInvalidAtomNumber, zatom, iso.get(), tmp.get());
  }

  // Attach any adjacent =O and =N
  IncludeInfoAboutJoined(m, iso.get(), tmp.get());

  if (_min_fragment_size > 0) {
    int atoms_in_frag = std::count_if(tmp.get(), tmp.get() + m.natoms(), [](int t) { return t > 0;});
    if (atoms_in_frag < _min_fragment_size) {
      ExpandBeyondCurrent(m, tmp.get(), _min_fragment_size - atoms_in_frag);
    }
  }

  // This could be folded into the next loop, but this is clearer.  
  // Cannot use m.set_isotopes
  for (int i = 0; i < initial_matoms; ++i) {
    if (iso[i]) {
      m.set_isotope(i, iso[i]);
    }
  }

  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (tmp[i] == 0) {
      m.remove_atom(i);
    }
  }

  int final_matoms = m.natoms();

  ++_atoms_lost[initial_matoms - final_matoms];
  ++_atoms_in_subset[final_matoms];

  // Re-identify the atom of interest.
  for (int i = 0; i < final_matoms; ++i) {
    if (m.isotope(i) == kDefaultIsotope) {
      mdata.set_centre_atom(i);
      return 1;
    }
  }

  cerr << "Options::ReduceToFragmentOfInterest:Hmmm, did not find atom with isotope " <<  kDefaultIsotope << '\n';
  mdata.set_centre_atom(kInvalidAtomNumber);
  cerr << m.smiles() << ' ' << m.name() << '\n';
  return 0;
}

int
Options::BuildDatabase(Molecule& m) {
  ++_molecules_read;

  MoleculeData mdata;
  if (! mdata.BuildWithKnownValue(m)) {
    cerr << "Cannot initialise molecule data\n";
    return 1;  // ignore errors
  }

  m.set_isotope(mdata.centre_atom(), kDefaultIsotope);

  if ( !ReduceToFragmentOfInterest(m, mdata, mdata.centre_atom())) {
    return 1;
  }

  _storage.StoreUsmi(m, mdata);
  return 1;

  PerMoleculeData pmd;
  pmd.Build(m);

  MaybeInitialiseAtomTypes(m, pmd);

  return BuildDatabase(m, mdata, pmd);
}

int
Options::MaybeInitialiseAtomTypes(Molecule& m, PerMoleculeData& pmd) {
  if (! _atom_typing_specification.active()) {
    return 0;
  }

  uint32_t* atype = pmd.get_atype(m.natoms());
  _atom_typing_specification.assign_atom_types(m, atype);

  return 1;
}

constexpr int kComplete = 2;
constexpr int kNextShell = 3;
constexpr int kProcessedHere = 4;

int
ExpandRadiusZero(Molecule& m,
                 atom_number_t zatom,
                 PerMoleculeData& pmd,
                 uint32_t& sum_so_far) {
//cerr << m.smiles() << " ExpandRadiusZero\n";
  sum_so_far = 0;
  for (const Bond* b : m[zatom]) {
    uint32_t bc = BondConstant(b);
    const atom_number_t o = b->other(zatom);

    Increment(bc, pmd.atype[o], sum_so_far);
    pmd.processing_status[o] = kNextShell;
  }

  return 1;
}

int
Options::BuildDatabase(Molecule& m,
                 MoleculeData& mdata,
                 PerMoleculeData& pmd) {
  const int matoms = m.natoms();

  std::fill_n(pmd.processing_status, matoms, 0);

  const atom_number_t zatom = mdata.centre_atom();

  pmd.processing_status[zatom] = kComplete;

  // We always include the first shell;
  uint32_t sum_so_far = 0;
  ExpandRadiusZero(m, zatom, pmd, sum_so_far);

  if (_store_radius_1_shells) {
    _storage.StoreShell(m, mdata, pmd, 1, sum_so_far);
  }

  return IdentifyShells(m, mdata, pmd, 2, sum_so_far);
}

#ifdef QQWE

  // The IWString_and_File_Descriptor object needs to be flushed.
  output.write_if_buffer_holds_more_than(4092);

  return 1;
}
#endif

int
Options::IdentifyShells(Molecule& m,
                        MoleculeData& mdata,
                        PerMoleculeData& pmd,
                        int radius,
                        uint32_t& sum_so_far) {
  int* processing_status = pmd.processing_status;
  const int matoms = m.natoms();

  int h = std::count_if(processing_status, processing_status + matoms, [](int p) {return p == kNextShell;});
  cerr << "Begin radius " << radius << " processing " << h << " atoms\n";

  // Keep track of what is added during this shell.
  Set_of_Atoms connections;
  connections.reserve(10);

  for (int i = 0; i < matoms; ++i) {
    if (processing_status[i] != kNextShell) {
      continue;
    }

    uint32_t added_this_atom = 0;
    for (const Bond* b : m[i]) {
      atom_number_t o = b->other(i);
      if (processing_status[o] == kComplete) {
        continue;
      }
      uint32_t bc = BondConstant(b);

      Increment(bc, pmd.atype[o], added_this_atom);
      connections << o;

      if (b->is_single_bond() && pmd.IsSaturatedCarbom(o)) {
        // AddTail();
      }
    }
    sum_so_far += added_this_atom;
    for (atom_number_t a : connections) {
      processing_status[a] = kProcessedHere;
    }
    connections.resize_keep_storage(0);
  }

  _storage.StoreShell(m, mdata, pmd, radius, sum_so_far);

  if (radius == _max_radius) {
    return 1;
  }

  for (int i = 0; i < matoms; ++i) {
    if (processing_status[i] == kNextShell) {
      processing_status[i] = kComplete;
    } else if (processing_status[i] == kProcessedHere) {
      processing_status[i] = kNextShell;
    }
  }

  return IdentifyShells(m, mdata, pmd, radius + 1, sum_so_far);
}

int
Options::WriteToDatabase() {
  if (! _build_database) {
    return 0;
  }

  return _storage.WriteToDatabase();
}

int
Options::WriteHash(IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteHash:cannot open '" << fname << "'\n";
    return 0;
  }

  return _storage.WriteHash(output);
}

// Use _queries to set the matched atoms in `process_atom`.
int
Options::IdentifyAtomsToProcess(Molecule& m, int* process_atom) {
  int rc = 0;

  Molecule_to_Match target(&m);

  for (Substructure_Query* q : _queries) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    sresults.each_embedding_set_vector(process_atom, 1);
    ++rc;
  }

  return rc;
}

int
Options::MaybeWriteUnknownFragment(const Molecule & parent, Molecule& m) {
  if (! _stream_for_fragments_not_found.is_open()) {
    return 0;
  }

  _stream_for_fragments_not_found << m.smiles() << ' ' << parent.name() << '\n';
  _stream_for_fragments_not_found.write_if_buffer_holds_more_than(4096);

  return 1;
}

// Return true if `zatom` is the Nitrogen atom of a Nitro group.
// Also works as the Sulphur in a O=S=O grouping.
int
IsNNitro(Molecule& m, atom_number_t zatom) {
  const Atom& a = m[zatom];
  if (a.atomic_number() == 7) {
    if (a.ncon() != 3) {
      return 0;
    }
  }

  int doubly_bond_oxygen_count = 0;
  for (const Bond* b : a) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (m.atomic_number(o) != 8) {
      continue;
    }
    ++doubly_bond_oxygen_count;
    if (doubly_bond_oxygen_count == 2) {
      return 1;
    }
  }

  return 0;
}

// By default we only consider N and O and S atoms.
int
DefaultAtomsToProces(Molecule& m, int* process_atom) {
  int rc = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {

    atomic_number_t z = m.atomic_number(i);
    if (z == 6) {
      continue;
    }
    if (z == 7) {
      if (IsNNitro(m, i)) {
        continue;
      }
    } else if (z == 8) {
    } else if (z == 16) {
      if (IsNNitro(m, i)) {   // sulfonamide
        continue;
      }
    } else {
      continue;
    }

    process_atom[i] = 1;
    ++rc;
  }

  return rc;
}

int
Options::DoLookups(Molecule& m,
                   IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  Preprocess(m);

  const int matoms = m.natoms();

  std::unique_ptr<int[]> process_atom = std::make_unique<int[]>(matoms);
  std::fill_n(process_atom.get(), matoms, 0);

  if (_queries.empty()) {
    DefaultAtomsToProces(m, process_atom.get());
  } else {
    if (! IdentifyAtomsToProcess(m, process_atom.get())) {
      return _ignore_molecules_not_matching_queries;
    }
  }

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (! process_atom[i]) {
      continue;
    }

    Molecule mcopy(m);
    mcopy.set_isotope(i, kDefaultIsotope);
    rc += DoLookups(mcopy, i, output);
  }

  ++_fragments_found[rc];

  return 1;
}

void
DoZeroPad(int width, const std::string& s, IWString& result) {
  int extra_chars = width - s.length();
  if (extra_chars <= 0) {
    return;
  }

  result.reserve(width);
  for (int i = 0; i < extra_chars; ++i) {
    result << '0';
  }

  result << s;
}

// `id` is the name of a molecule.
// If we are NOT zero padding identifiers, set `result` to point to the data in `id`.
// If we are zero padding, use `tmp` to hold permanent storage of the zero padded
// identifier, and place the contents of `tmp` into `result`.
int
FormKey(const std::string& id, IWString& tmp,
        int zero_pad_identifiers,
        Dbt& result) {
  if (zero_pad_identifiers == 0) {
    result.set_data((void*)id.data());
    result.set_size(id.size());
    return 1;
  }

  DoZeroPad(zero_pad_identifiers, id, tmp);
  result.set_data((void*)tmp.data());
  result.set_size(tmp.length());

  return 1;
}

int
Options::WriteComputedValue(const qupkake_transfer::ParentValue& c,
                            const char* acid_base,
                            IWString_and_File_Descriptor& output) {
  if (_selimsteg.empty()) {
    return 1;
  }

  const std::string& id = c.par();

  // If we zero pad the identifier, we need something with scope here.
  IWString tmp;

  static constexpr char kSep = ' ';

  if (_selimsteg.size() != _zero_pad_identifiers.size()) {
    cerr << "Options::WriteComputedValue:selimsteg/zero pad size mismatch\n";
    return 0;
  }

  Dbt fromdb;
  for (uint32_t i = 0; i < _selimsteg.size(); ++i) {
    Dbt dbkey;
    FormKey(id, tmp, _zero_pad_identifiers[i], dbkey);

    if (int rc = _selimsteg[i]->get(NULL, &dbkey, &fromdb, 0); rc != 0) {
      // cerr << "Not found ";
      // db->err(rc, "");
      continue;
    }
    const_IWSubstring s((const char*) fromdb.get_data(), fromdb.get_size());
    output << s << kSep << id << kSep << c.value() << kSep << acid_base << '\n';
    return 1;
  }

  output << c.par() << kSep << c.value() << '\n';

  return 1;
}

// Return istep if `nvalues` have been found in the database.
int
Options::ComputeIstep(int nvalues) const {
  // all values can be shown.
  if (nvalues <= _nwrite) {
    return 1;
  }

  return nvalues / _nwrite + 1;
}

int
Options::LookupViaUniqueSmiles(Molecule& parent, Molecule& m, atom_number_t zatom,
                IWString_and_File_Descriptor& output) {
  const IWString& usmi = m.unique_smiles();

  Dbt key;
  key.set_data((void*) usmi.data());
  key.set_size(usmi.length());

  std::optional<qupkake_transfer::QupKakeData> maybe_proto = _storage.Lookup(usmi);
  if (! maybe_proto) {
    MaybeWriteUnknownFragment(parent, m);
    return 0;
  }

  static constexpr char kSep = ' ';  // Not fully implemented

  output << parent.smiles() << kSep << parent.name() << kSep << parent.smarts_equivalent_for_atom(zatom) << '\n';
  if (maybe_proto->acid().values_size() > 0) {
    const ComputedValues& acid = maybe_proto->acid();
    output << m.smiles() << " acid " << acid.values_size() << kSep << acid.minval() << kSep << acid.mean() << kSep << acid.maxval() << '\n';

    int istep = ComputeIstep(acid.values_size());

    for (int i = 0; i < acid.values_size(); i += istep) {
      WriteComputedValue(acid.values()[i], "acid", output);
    }
  }
  if (maybe_proto->base().values_size() > 0) {
    const ComputedValues& base = maybe_proto->base();
    output << m.smiles() << " base " << base.values_size() << kSep << base.minval() << kSep << base.mean() << kSep << base.maxval() << '\n';

    int istep = ComputeIstep(base.values_size());
    for (int i = 0; i < base.values_size(); i += istep) {
      WriteComputedValue(base.values()[i], "base", output);
    }
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::DoLookups(Molecule& m,
                   atom_number_t zatom,
                   IWString_and_File_Descriptor& output) {
  Molecule parent(m);

  MoleculeData mdata;
  mdata.BuildNoKnownValue(m);

  if (!ReduceToFragmentOfInterest(m, mdata, zatom)) {
    return 0;
  }

  return LookupViaUniqueSmiles(parent, m, zatom, output);

  PerMoleculeData pmd;
  pmd.Build(m);

  MaybeInitialiseAtomTypes(m, pmd);

  return DoLookups(m, mdata, pmd, output);
}

int
Options::DoLookups(Molecule& m, MoleculeData& mdata, PerMoleculeData& pmd,
                   IWString_and_File_Descriptor& output) {
  uint32_t sum_so_far = 0;
  ExpandRadiusZero(m, mdata.centre_atom(), pmd, sum_so_far);

  return 0;
}

int
ApplicationName(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  if (options.build_database()) {
    return options.BuildDatabase(m);
  }

  return options.DoLookups(m, output);
}

int
ApplicationName(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! ApplicationName(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ApplicationName(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:R:P:d:M:S:s:z:G:U:n:m:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! ApplicationName(options, fname, input_type, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (options.build_database()) {
    options.WriteToDatabase();
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    options.WriteHash(fname);
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace qupkake_transfer

int
main(int argc, char ** argv) {

  int rc = qupkake_transfer::Main(argc, argv);

  return rc;
}
