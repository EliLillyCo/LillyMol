#ifndef MOLECULE_TOOLS_MOLECULE_DATABASE_OPTIONS_H_
#define MOLECULE_TOOLS_MOLECULE_DATABASE_OPTIONS_H_

#include <optional>

#include "leveldb/db.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/molecule_database_options.pb.h"

namespace molecule_database {

inline const char* kStructureTypeKey = "_StructureTtypeKey";

void Preprocess(Molecule& m, const LLYMol::MoleculeDatabaseOptions& options,
          Chemical_Standardisation& chemical_standardisation);

// Store the structure specifications in kStructureTypeKey.
int StoreStructureTypeInformation(const LLYMol::MoleculeDatabaseOptions& options,
        const leveldb::WriteOptions& write_options,
        leveldb::DB* database);

// Retrieve the structure specifications stored in kStructureTypeKey.
std::optional<LLYMol::MoleculeDatabaseOptions> GetStructureTypeInformation(const leveldb::ReadOptions& read_options,
        leveldb::DB* database);

std::optional<LLYMol::MoleculeDatabaseOptions> OptionsFromCmdline(Command_Line& cl, int verbose);
}  // namespace molecule_database  
#endif  // MOLECULE_TOOLS_MOLECULE_DATABASE_OPTIONS_H_
