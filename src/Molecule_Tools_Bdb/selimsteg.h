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

    // This owns _database. Copying it would give two objects the same Db* and
    // both would close and delete it. The implicitly generated copy constructor
    // did exactly that, and it was reachable - a python binding returning
    // Selimsteg& with the default return value policy copies, so a with block
    // ended up closing a different object than the one that had opened the
    // database, leaving the file open.
    Selimsteg(const Selimsteg&) = delete;
    Selimsteg& operator=(const Selimsteg&) = delete;

    bool OpenDatabase(const std::string& dbname);

    // Close the database, if one is open. Idempotent, and called by the
    // destructor, so there is one release path.
    //
    // Worth having explicitly rather than leaving it to the destructor. The open
    // database holds a file handle, and from python the only way to release it
    // was to drop the last reference to the object. That is invisible when it
    // goes wrong: unlinking the database file while it is still open succeeds on
    // a local filesystem but on NFS becomes a silly rename to .nfsXXXX, so a
    // directory holding it cannot be removed. Also lets a caller open a
    // different database on the same object, which OpenDatabase otherwise
    // refuses.
    void Close();

    // Returns nullopt when no database is open, rather than dereferencing a null
    // pointer as this used to.
    std::optional<std::string> Lookup(const std::string& key);
    std::optional<Molecule> GetMolecule(const std::string& key);
    std::vector<Molecule> GetMolecules(const std::vector<std::string>& keys);
};

}  // selimsteg
#endif // MOLECULE_TOOLS_BDB_H
