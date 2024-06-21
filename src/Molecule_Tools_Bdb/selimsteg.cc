#include <iostream>

#include "selimsteg.h"

namespace selimsteg {

using std::cerr;

Selimsteg::Selimsteg() {
  _database = nullptr;
}

Selimsteg::~Selimsteg() {
  if (_database) {
    delete _database;
  }
}

bool
Selimsteg::OpenDatabase(const std::string& dbname) {
  if (_database) {
    cerr << "Selimsteg::OpenDatabase:already open\n";
    return false;
  }

  _database = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

  int rc = _database->open(NULL, dbname.c_str(), NULL, DB_UNKNOWN,
                                 DB_RDONLY, 0);

  if (0 != rc) {
    cerr << "Selimsteg::OpenDatabase::cannot open '" << dbname << "'\n";
    _database->err(rc, "");
    delete _database;
    _database = nullptr;
    return false;
  }

  return true;
}

std::optional<std::string>
Selimsteg::Lookup(const std::string& key) {
  Dbt dbkey;

  dbkey.set_data((void*)(key.data()));  // loss of const OK
  dbkey.set_size(key.size());

  Dbt fromdb;
  if (0 != _database->get(NULL, &dbkey, &fromdb, 0)) {
    return std::nullopt;
  }

  return std::string(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
}

std::optional<Molecule>
Selimsteg::GetMolecule(const std::string& key) {
  std::optional<std::string> maybe_smiles = Lookup(key);
  if (! maybe_smiles) {
    return std::nullopt;
  }

  Molecule m;
  if (! m.build_from_smiles(*maybe_smiles)) {
    cerr << "Selimsteg::GetSmiles:invalid smiles '" << *maybe_smiles
         << " key " << key << '\n';
    return std::nullopt;
  }

  m.set_name(key);

  return m;
}

std::vector<Molecule>
Selimsteg::GetMolecules(const std::vector<std::string>& keys) {
  const uint32_t nmols = keys.size();
  std::vector<Molecule> result(nmols);
  for (uint32_t i = 0; i < nmols; ++i) {
    std::optional<std::string> maybe_smiles = Lookup(keys[i]);
    if (! maybe_smiles) {
      continue;
    }
    if (result[i].build_from_smiles(*maybe_smiles)) {
      result[i].set_name(keys[i]);
    } else {
      cerr << "Selimsteg::GetMolecules:invalid smiles " << *maybe_smiles 
           << " key " << keys[i] << '\n';
    }
  }

  return result;
}

}  // namespace selimsteg
