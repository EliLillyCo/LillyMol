#ifndef MOLECULE_TOOLS_BDB_STRUCTURE_DATABASE_H_
#define MOLECULE_TOOLS_BDB_STRUCTURE_DATABASE_H_

#include "db_cxx.h"

#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

namespace structure_database {

// During lookups we can control what gets returned via the bitwise OR of these
enum Lookup {
  kNothing = 0,
  // No fragment stripping, chirality preserved
  kExact = 1,
  // Strip to largest fragment
  kStrip = 2,
  // Remove chirality
  kNoChiral = 4,
  // Graph
  kGraph = 8,
  // Whether or not to apply standardisation. If not set, molecules are
  // standardised, but if set that is skipped.
  kNoStandardise = 16
};

class StructureDatabase {
  private:
    // We have two databases, one for smiles,
    //   all fragments
    //   largest fragment
    //   largest fragment, no chirality
    // and
    //   a graph database.
    Db* _database;
    Db* _graph_database;

    Mol2Graph _mol2graph;

    Chemical_Standardisation _chemical_standardisation;

    int _remove_invalid_chirality;

    // Useful during debugging.
    int _echo_lookup_key = 0;

  // private functions.
    bool CommonOpenDatabaseForReading(IWString& dbname, Db*& database);
    bool CommonOpenDatabaseForWriting(IWString& dbname, Db*& database);
    void FormGraphDatabaseKey(Molecule& m, IWString& key);
    void Preprocess(Molecule& m);
    int DoStore(const IWString& key, const IWString& value, Db* database);
    int DoLookup(Db* database, const IWString& key, IWString& ids_matched);
    int AppendToExisting(Dbt& dbkey, const Dbt& fromdb, 
                const IWString& value, Db* database);

  public:
    StructureDatabase();
    ~StructureDatabase();

    // Structure databases come in pairs, one db for smiles and one for
    // graph representations. These functions will open both databases.
    int OpenForReading(const IWString& stem);
    int OpenForWriting(const IWString& stem);

    int OpenStructureDatabaseForReading(IWString& dbname);
    int OpenStructureDatabaseForWriting(IWString& dbname);
    int OpenGraphDatabaseForReading(IWString& dbname);
    int OpenGraphDatabaseForWriting(IWString& dbname);

    int TurnOffChemicalStandardisation();

    int Store(Molecule& m);

    int Lookup(Molecule& m, enum Lookup params, IWString& ids_matched);
    int Lookup(Molecule& m, uint32_t mask, IWString& ids_matched);
};

}  // namespace structure_database

#endif // MOLECULE_TOOLS_BDB_STRUCTURE_DATABASE_H_
