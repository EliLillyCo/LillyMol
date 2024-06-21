#include <iostream>

#include <google/protobuf/util/message_differencer.h>

#include "Molecule_Lib/mol2graph_proto.h"

#include "molecule_database_options.h"

namespace  molecule_database {

using std::cerr;

void Preprocess(Molecule& m, const LLYMol::MoleculeDatabaseOptions& options,
                Chemical_Standardisation& chemical_standardisation) {
  if (options.reduce_to_largest_fragment()) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (!options.include_chirality()) {
    m.remove_all_chiral_centres();
  }

  if (!options.include_directional_bonds()) {
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }
}

int StoreStructureTypeInformation(const LLYMol::MoleculeDatabaseOptions& options,
        const leveldb::WriteOptions& write_options,
        leveldb::DB* database) {
  leveldb::Slice key(kStructureTypeKey);
  std::string existing_value;
  const leveldb::ReadOptions read_options;
  leveldb::Status status = database->Get(read_options, key, &existing_value);
  if (status.ok()) {  // Existing value present.
    LLYMol::MoleculeDatabaseOptions current;
    if (! current.ParseFromString(existing_value)) {
      cerr << "StoreStructureTypeInformation:could not parse existing " << kStructureTypeKey << "\n";
      return 0;
    }
    if (! google::protobuf::util::MessageDifferencer::Equivalent(options, current)) {
      cerr << "StoreStructureTypeInformation:Incompatible storage conditions\n";
      cerr << "proposed\n" << options.DebugString() << "\n";
      cerr << "existing\n" << current.DebugString() <<  "\n";
      return 0;
    }

    return 1;  // New is same as current, no need to store.
  }


  std::string serialized;
  options.SerializeToString(&serialized);
  leveldb::Slice value(serialized);
  status = database->Put(write_options, key, value);
  if (status.ok()) {
    return 1;
  }

  cerr << "Did not store structure information key\n";
  return 0;
}

std::optional<LLYMol::MoleculeDatabaseOptions> OptionsFromCmdline(Command_Line& cl, int verbose) {
  LLYMol::MoleculeDatabaseOptions to_be_returned;

  if (cl.option_present('l')) {
    to_be_returned.set_reduce_to_largest_fragment(true);
    if (verbose)
      cerr << "Will strip to largest fragment\n";
  }

  if (cl.option_present('c')) {
    to_be_returned.set_include_chirality(false);
    if (verbose)
      cerr << "WIll exclude chirality\n";
  }

  if (!cl.option_present('G')) {
    return to_be_returned;
  }

  Mol2Graph mol2graph;

  if (! mol2graph.construct(cl, 'G', verbose)) {
    cerr << "Cannot initialize mol2graph (-G)\n";
    return std::nullopt;
  }

  *to_be_returned.mutable_mol2graph() = Mol2GraphToProto(mol2graph);

  return to_be_returned;
}

std::optional<LLYMol::MoleculeDatabaseOptions> GetStructureTypeInformation(const leveldb::ReadOptions& read_options,
        leveldb::DB* database) {
  leveldb::Slice key(kStructureTypeKey);
  std::string value;
  const leveldb::Status status = database->Get(read_options, key, &value);
  if (! status.ok()) {
    cerr << "GetStructureTypeInformation:cannot retrieve " << kStructureTypeKey << "\n";
    return std::nullopt;
  }

  LLYMol::MoleculeDatabaseOptions to_be_returned;
  if (! to_be_returned.ParseFromString(value)) {
    cerr << "GetStructureTypeInformation:cannot parse proto " << kStructureTypeKey << "\n";
    return std::nullopt;
  }

  return to_be_returned;
}

}  // namespace molecule_database
