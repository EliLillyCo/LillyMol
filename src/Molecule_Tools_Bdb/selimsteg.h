#ifndef MOLECULE_TOOLS_BDB_H
#define MOLECULE_TOOLS_BDB_H

#include <memory>
#include <optional>
#include <string>

#include "db_cxx.h"

#include "Molecule_Lib/molecule.h"

namespace selimsteg {

class Selimsteg {
  private:
    //This just never worked, not sure why.
    //std::unique_ptr<Db> _database;

    Db* _database;

  public:
    Selimsteg();
    ~Selimsteg();

    bool OpenDatabase(const std::string& dbname);

    std::optional<std::string> Lookup(const std::string& key);
    std::optional<Molecule> GetMolecule(const std::string& key);
    std::vector<Molecule> GetMolecules(const std::vector<std::string>& keys);
};

}  // selimsteg
#endif // MOLECULE_TOOLS_BDB_H
