#include <iostream>

#include "Foundational/iwmisc/misc.h"

#include "structure_database.h"

namespace structure_database {

using std::cerr;

StructureDatabase::StructureDatabase() {
  _database = nullptr;
  _graph_database = nullptr;

  _mol2graph.TurnOnMostUsefulOptions();

  _chemical_standardisation.activate_all();

  _echo_lookup_key = 0;
}

StructureDatabase::~StructureDatabase() {
  static constexpr uint32_t kFlags = 0;
  if (_database != nullptr) {
    _database->close(kFlags);  // probably not necessary
    delete _database;
  }
  if (_graph_database != nullptr) {
    _graph_database->close(kFlags);  // probably not necessary
    delete _graph_database;
  }
}

int
StructureDatabase::TurnOffChemicalStandardisation() {
  _chemical_standardisation.deactivate();
  return 1;
}

int
StructureDatabase::OpenStructureDatabaseForReading(IWString& dbname) {
  return CommonOpenDatabaseForReading(dbname, _database);
}

int
StructureDatabase::OpenGraphDatabaseForReading(IWString& dbname) {
  return CommonOpenDatabaseForReading(dbname, _graph_database);
}

bool
StructureDatabase::CommonOpenDatabaseForReading(IWString& dbname,
                                Db*& database) {
  if (database) {
    cerr << "StructureDatabase::CommonOpenDatabaseForReading:already open "
         << dbname << '\n';
    return false;
  }
  
  database = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

  int rc = database->open(NULL, dbname.null_terminated_chars(), NULL, DB_UNKNOWN,
                          DB_RDONLY, 0);

  if (0 != rc) {
    cerr << "StructureDatabase::CommonOpenDatabaseForReading::cannot open '" << dbname << "'\n";
    database->err(rc, "");
    delete database;
    database = nullptr;
    return false;
  }

  return true;
}

int
StructureDatabase::OpenStructureDatabaseForWriting(IWString& dbname) {
  return CommonOpenDatabaseForWriting(dbname, _database);
}

int
StructureDatabase::OpenGraphDatabaseForWriting(IWString& dbname) {
  return CommonOpenDatabaseForWriting(dbname, _graph_database);
}

int
StructureDatabase::OpenForReading(const IWString& stem) {
  IWString dbname;
  dbname << stem << ".smi.bdb";
  if (! OpenStructureDatabaseForReading(dbname)) {
    cerr << "StructureDatabase::OpenForReading:cannot open '" << dbname << "'\n";
    return 0;
  }

  dbname.resize_keep_storage(0);
  dbname << stem << ".graph.bdb";
  if (! OpenGraphDatabaseForReading(dbname)) {
    cerr << "StructureDatabase::OpenForReading:cannot open '" << dbname << "'\n";
    return 0;
  }

#ifdef USEFUL_FOR_DEBUGGING
  Dbc* cursor = NULL;

  int rc = _graph_database->cursor(NULL, &cursor, 0);
  if (0 != rc) {
    _graph_database->err(rc, "cannot acquire cursor");
    return 0;
  }

  Dbt zkey, zdata;

  int records_in_database = 0;
  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    records_in_database++;
    if (records_in_database > 10) {
      break;
    }
    IWString key((const char*)zkey.get_data(), zkey.get_size());
    IWString data((const char*)zdata.get_data(), zdata.get_size());
    // cerr << "key '" << key << "' data '" << data << "'\n";
  }
#endif


  return 1;
}

int
StructureDatabase::OpenForWriting(const IWString& stem) {
  IWString dbname;
  dbname << stem << ".smi.bdb";
  if (! OpenStructureDatabaseForWriting(dbname)) {
    cerr << "StructureDatabase::OpenForWriting:cannot open '" << dbname << "'\n";
    return 0;
  }

  dbname.resize_keep_storage(0);
  dbname << stem << ".graph.bdb";
  if (! OpenGraphDatabaseForWriting(dbname)) {
    cerr << "StructureDatabase::OpenForWriting:cannot open '" << dbname << "'\n";
    return 0;
  }
  return 1;
}

bool
StructureDatabase::CommonOpenDatabaseForWriting(IWString& dbname,
                                Db*& database) {
  if (database) {
    cerr << "StructureDatabase::CommonOpenDatabaseForReading:already open "
         << dbname << '\n';
    return false;
  }

  int flags;
  DBTYPE dbtype;
  int mode;

  if (dash_s(dbname)) {
    dbtype = DB_UNKNOWN;
    flags = 0;
    mode = 0;
  } else {
    dbtype = DB_BTREE;
    flags = DB_CREATE;
    mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;
  }

  database = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

  int rc = database->open(NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);
  if (rc != 0) {
    cerr << "StructureDatabase::CommonOpenDatabaseForWriting::cannot open '" << dbname << "'\n";
    database->err(rc, "");
    delete database;
    database = nullptr;
    return false;
  }

  return true;
}

void
Prepend(const char* prefix, IWString& destination) {
  IWString tmp;
  tmp.reserve(strlen(prefix) + destination.length());
  tmp << prefix << destination;
  destination = tmp;
}

// Place in `key` the unique smiles + formula key used in graph databases.
void
StructureDatabase::FormGraphDatabaseKey(Molecule& m, IWString& key) {
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();

  m.compute_aromaticity_if_needed();
  IWString formula;
  m.formula_distinguishing_aromatic(formula);

  m.change_to_graph_form(_mol2graph);
  const IWString& usmi = m.unique_smiles();

  key.resize_keep_storage(0);

  key.reserve(usmi.size() + 1 + formula.size());
  key << m.unique_smiles() << ':' << formula;
}

// things that must be done both on storing and lookup.
void
StructureDatabase::Preprocess(Molecule& m) {
  // LillyMol does not properly handle cis trans bonding.
  m.revert_all_directional_bonds_to_non_directional();

  // We could allow isotopes if there was a need.
  m.unset_isotopes();
}

int
StructureDatabase::Store(Molecule& m) {
  Preprocess(m);

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  IWString mname = m.name();
  mname.truncate_at_first(' ');

  if (! DoStore(m.unique_smiles(), mname, _database)) {
    cerr << "StructureDatabase::Store:cannot store '" << m.name() << "'\n";
    return 0;
  }

  if (m.number_fragments() > 1) {
    m.reduce_to_largest_fragment_carefully();
    Prepend("F%", mname);
    DoStore(m.unique_smiles(), mname, _database);
  }

  if (m.chiral_centres()) {
    m.remove_all_chiral_centres();
    Prepend("@%", mname);
    DoStore(m.unique_smiles(), mname, _database);
  }

  if (_graph_database) {
    IWString key;
    FormGraphDatabaseKey(m, key);
    DoStore(key, m.name(), _graph_database);
  }

  return 1;
}

int
StructureDatabase::DoStore(const IWString& key, const IWString& value, Db* database) {
  Dbt dbkey;
  dbkey.set_data(key.rawdata());
  dbkey.set_size(key.length());

  Dbt fromdb;
  if (database->get(NULL, &dbkey, &fromdb, 0) == 0) {
    return AppendToExisting(dbkey, fromdb, value, database);
  }

  Dbt dbvalue;
  dbvalue.set_data(value.rawdata());
  dbvalue.set_size(value.length());

  const int rc = database->put(NULL, &dbkey, &dbvalue, 0);
  if (rc == 0) {
    return 1;
  }

  cerr << "StructureDatabase::DoStore:cannot store '" << key << "' value '" << value << "'\n";
  database->err(rc, "");
  return 0;
}
    
int
StructureDatabase::AppendToExisting(Dbt& dbkey, const Dbt& fromdb, 
                const IWString& value, Db* database) {
  IWString to_store(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
  to_store.append_with_spacer(value, ':');

  Dbt dbvalue;
  dbvalue.set_data((void*)to_store.data());
  dbvalue.set_size(to_store.size());

  const int rc = database->put(NULL, &dbkey, &dbvalue, 0);
  if (rc == 0) {
    return 1;
  }

  IWString tmp(reinterpret_cast<const char*>(dbkey.get_data()), dbkey.get_size());
  cerr << "StructureDatabase::append_with_spacer:cannot store '" << tmp << "' value '" << value << "'\n";
  database->err(rc, "");
  return 0;
}

int
StructureDatabase::Lookup(Molecule& m, uint32_t params, IWString& ids_matched) {
  // cerr << "params " <<params << '\n';
  ids_matched.resize_keep_storage(0);

  Molecule mcopy(m);

  Preprocess(mcopy);

  if (! _chemical_standardisation.active()) {
  } else if (params & Lookup::kNoStandardise) {
  } else {
    _chemical_standardisation.process(mcopy);
  }

  int rc = 0;
  if (params & Lookup::kExact) {
    rc += DoLookup(_database, mcopy.unique_smiles(), ids_matched);
  }

  if (params & Lookup::kStrip && mcopy.number_fragments()) {
    mcopy.reduce_to_largest_fragment_carefully();
    rc += DoLookup(_database, mcopy.unique_smiles(), ids_matched);
  }

  if (params & Lookup::kNoChiral && mcopy.chiral_centres()) {
    mcopy.remove_all_chiral_centres();
    rc += DoLookup(_database, mcopy.unique_smiles(), ids_matched);
  }

  if ((params & Lookup::kGraph) == 0) {
    return rc;
  }

  if (! _graph_database) {
    cerr << "StructureDatabase::Lookup:graph database not open\n";
    return rc;
  }

  IWString key;
  FormGraphDatabaseKey(mcopy, key);

  return rc + DoLookup(_graph_database, key, ids_matched);
}
  
int
StructureDatabase::DoLookup(Db* database, const IWString& key, IWString& ids_matched) {
  Dbt dbkey;
  dbkey.set_data(key.rawdata());
  dbkey.set_size(key.length());

  if (_echo_lookup_key) {
    cerr << "StructureDatabase::DoLookup:Lookup '" << key << "' " << database << '\n';
  }

  Dbt fromdb;
  if (0 != database->get(NULL, &dbkey, &fromdb, 0)) {
    return 0;
  }

  const_IWSubstring tmp((const char *)(fromdb.get_data()), fromdb.get_size());
  ids_matched.append_with_spacer(tmp, ':');

  return 1;
}

}  // namespace structure_database
