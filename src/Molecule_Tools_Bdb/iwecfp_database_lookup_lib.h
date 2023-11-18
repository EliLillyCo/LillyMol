#ifndef MOLECULE_TOOLS_BDB_IWECF_DATABASE_LOOKUP_LIB_
#define MOLECULE_TOOLS_BDB_IWECF_DATABASE_LOOKUP_LIB_

#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>

#include "db_cxx.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/istream_and_type.h"

#include "iwecfp_database.h"

namespace iwecfp_database_lookup {

using iwecfp_database::Count_Radius;
using iwecfp_database::DBKey;
using iwecfp_database::IWdbkeyHash;

// A class describing a BerkeleyDB database built by iwecfp_database_load.
class SP_Database
{
 private:
  Db* _db;
  DbEnv* _env;

  int _lookups_done;
  int _bits_found;

  // We need a mutex to protect access to the _hash related things
  // Initial implementation (cilk) did not help, so currently not active.

  std::mutex _mutex;

  // We hash lookups for speed.
  std::unordered_map<DBKey, Count_Radius, IWdbkeyHash> _hash;
  int _hash_size;
  int _cache_hit;

  // The database may, or may not contain example structures.
  // Usually not, storing examples results in databases that are
  // too large.
  int _contains_examples;

  IWString _dbname;

  int _logfile_open;

  std::ofstream _logfile;

  int _need_to_swap_bytes;  // not used

  //  private functions

  int _do_open(const char*, int multi_threaded);
  int _parse_database_record(const const_IWSubstring& fromdb, const DBKey& dbkey, int& radius,
                         int& count, IWString& example);
  int _in_cache(Dbt& zkey, const DBKey& dbkey, int& radius, int& count);
  int _in_cache_mutex_protected(Dbt& zkey, const DBKey& dbkey, int& radius, int& count);

 public:
  SP_Database();
  ~SP_Database();

  int report(std::ostream& output) const;
  int debug_print(std::ostream& output) const;

  int open_logfile(const char* fname);

  int determine_if_examples_present();

  int contains_example_structures() const {
    return _contains_examples;
  }

  void set_env(DbEnv* e) {
    _env = e;
  }

  int open(const char*, int multi_threaded);
  //  int open(const const_IWSubstring & s, int ndx);

  unsigned int slurp_to_cache(int);

  int check_exists(Dbt* zkey);

  int do_lookup(Dbt& zkey, const DBKey& dbkey, int& radius, int& count, IWString& example);
  int get(DbTxn*, Dbt&, Dbt&, u_int32_t);
};

class Count_Example {
 private:
  const int _number_examples;

  const int _radius;

  IWString _first_example_smiles;
  IWString _first_example_name;

 public:
  //  Count_Example(Molecule &, atom_number_t centre_atom, int r);
  Count_Example(const IWString&, int r, int c);

  int
  radius() const {
    return _radius;
  }

  void
  set_first_example_smiles(const const_IWSubstring& s) {
    _first_example_smiles = s;
  }

  int
  number_examples() const {
    return _number_examples;
  }

  //  const IWString & first_example() const {return _first_example;}

  int do_output(IWString_and_File_Descriptor&) const;
};

/*
  When doing retrospective studies, it can be interesting to look
  at the database without certain molecules
*/

typedef std::unordered_map<DBKey, int, IWdbkeyHash> Subtract;

class Subtraction_Set {
 private:
  Subtract _subtract;
  int _molecules_in_subtraction_set;

  std::atomic<int> _subtraction_bits_hit;

  //  private functions

  int
  _transfer_fingerprints_to_subtract(const iwecfp_database::Set_of_Bits<iwecfp_database::Bit_Produced>& sob);
  int
  _adjust_count_for_subtraction(const DBKey& dbkey, int& count) const;

 public:
  Subtraction_Set();

  size_t
  nbits() const
  {
    return _subtract.size();
  }

  int
  report(std::ostream&) const;

  typedef int (*Preprocessor)(Molecule& m);

  int
  build(const Command_Line& cl, const char flag, iwecfp_database::Fingerprint_Characteristics& rpc,
        Preprocessor preprocess_molecule);
  int
  build(IWString& fname, iwecfp_database::Fingerprint_Characteristics& fpc,
        Preprocessor preprocess_molecule);
  int
  build(data_source_and_type<Molecule>& input, iwecfp_database::Fingerprint_Characteristics& fpc,
        Preprocessor preprocess_molecule);

  int
  nmolecules() const
  {
    return _molecules_in_subtraction_set;
  }

  int
  active() const
  {
    return _molecules_in_subtraction_set;
  }

  int
  adjust_count_for_subtraction(const DBKey& dbkey, int& count);

  int
  adjust_count_for_subtraction_threaded(
      const DBKey& dbkey, int& count);  // does not adjust _subtraction_bits_hit
};

class Set_of_Databases {
 private:
  resizable_array_p<SP_Database> _db;
  DbEnv _dbenv;

  int _user_wants_example_structures;

  unsigned char _userdata[512];  // user memory for retrievals

  Subtraction_Set _subtraction;

  //  private functions

  int
  _determine_if_examples_present();

 public:
  Set_of_Databases();
  ~Set_of_Databases();

  int number_databases() const {
    return _db.number_elements();
  }

  int debug_print(std::ostream& output) const;

  int open_logfile(const char* stem);

  int build(Command_Line& cl, char flag_env, char flag_db, const int verbose);
  int build(const Command_Line& cl, const char flag_db, const int verbose, DbEnv* envptr);

  // Used by the python api.
  int AddDatabase(const std::string& fname);

  int BuildSubtractionDatabase(IWString& fname, iwecfp_database::Fingerprint_Characteristics& fpc,
                int(*preprocess_molecule)(Molecule& m));
  int ReportSubtractDatabaseStatus(std::ostream& output) const;
  int ReportSubtractionDatabaseImpact(std::ostream& output) const;

  int report(std::ostream&) const;

  void set_user_wants_example_structures(int s) {
    _user_wants_example_structures = s;
  }

  int user_wants_example_structures() const {
    return _user_wants_example_structures;
  }

  int determine_atom_typing_in_use(IWString& string_atype, int& iwecfp_atom_type);
  int determine_max_search_radius(int& min_radius);

  int slurp_to_cache(const int s);

  int lookup_bit(const DBKey& dbkey, int& count, IWString& example_structure);
  int lookup_bit(const DBKey& dbkey, resizable_array_p<Count_Example>& ce);
  int lookup_bit_threaded(const DBKey& dbkey);
};

// Needed for the python interface, and a more generally useful class.
class SP_Set_of_Databases {
  private:
    Set_of_Databases _dbs;
    // As we open databases, we determine the atom typing. That
    // must be consistent across databases.
    const IWString _ats;

    // And using _ats, instantiate an actual typing.
    Atom_Typing_Specification _atom_typing_specification;

    iwecfp_database::Fingerprint_Characteristics _fingerprint_characteristics;

  public:
    int AddDatabase(const std::string& dbname);

    // Note that currently there is no checking whether max_radius might be
    // larger than what is in any of the databases.
    int set_max_radius(int max_radius) {
      _fingerprint_characteristics.set_max_shell_radius(max_radius);
      return 1;
    }

    int PerShellData(Molecule& m, std::vector<int>& result);
};


void set_max_hash_size(int s);
void set_verbose(int s);
void set_slurp_examples(int s);

} // namespace iwecfp_database_lookup

#endif // MOLECULE_TOOLS_BDB_IWECF_DATABASE_LOOKUP_LIB_
