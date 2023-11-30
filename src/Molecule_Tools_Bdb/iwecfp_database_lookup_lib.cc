// Functions supporting lookups in synthetic precedent databases

#include <limits>
#include <memory>
#include <string>

#include "iwecfp_database.h"
#include "iwecfp_database_lookup_lib.h"

namespace iwecfp_database_lookup {

using std::cerr;
using iwecfp_database::Fingerprint_Characteristics;
using iwecfp_database::Set_of_Bits;
using iwecfp_database::Bit_Produced;

int max_hash_size = 0;

void
set_max_hash_size(int s) {
  max_hash_size = s;
}

int verbose = 0;
void
set_verbose(int s) {
  verbose = s;
}

int slurp_examples = 0;
void
set_slurp_examples(int s) {
  slurp_examples = s;
}

SP_Database::SP_Database()
{
  _db = nullptr;
  _env = nullptr;

  _lookups_done = 0;
  _bits_found = 0;
  _hash_size = 0;
  _cache_hit = 0;

  _contains_examples = 1;  // until proven otherwise

  _logfile_open = 0;

  _need_to_swap_bytes = 0;

  return;
}

SP_Database::~SP_Database()
{
  if (nullptr != _db) {
    delete _db;
  }

  return;
}

int
SP_Database::debug_print(std::ostream& output) const
{
  output << "SP_Database::debug_print:database '" << _dbname << "' _contains_examples "
         << _contains_examples << " _hash_size " << _hash_size << '\n';

  return 1;
}

int
SP_Database::report(std::ostream& os) const
{
  os << "Database '" << _dbname << "' looked up " << _lookups_done << " bits, found "
     << _bits_found << " bits\n";
  if (_hash.size() > 0) {
    os << "Bit hash contains " << _hash.size() << " bits, " << _cache_hit
       << " cache hits\n";
  }

  return 1;
}

int
SP_Database::open(const char* dbname, int multi_threaded)
{
  return _do_open(dbname, multi_threaded);
}

int
SP_Database::open_logfile(const char* fname)
{
  _logfile.open(fname, std::ios::out);

  if (!_logfile.good()) {
    cerr << "SP_Database::open_logfile:cannot open '" << fname << "'\n";
    return 0;
  }

  _logfile_open = 1;

  return 1;
}

int
SP_Database::_do_open(const char* dbname, int multi_threaded)
{
  if (nullptr != _db) {
#ifdef DEBUG_SEQUENCE_DATABASE_STUFF
    cerr << "Deleting old database\n";
#endif

    delete _db;
  }

  _db = new Db(_env, DB_CXX_NO_EXCEPTIONS);

#ifdef DEBUG_SEQUENCE_DATABASE_STUFF
  cerr << "Opening '" << dbname << "'\n";
#endif

  int oflags = DB_RDONLY;
  if (multi_threaded) {
    oflags |= DB_THREAD;
  }

  int rc = _db->open(NULL, dbname, NULL, DB_UNKNOWN, oflags, 0);

  if (0 != rc) {
    cerr << "Cannot open database '" << dbname << "' ";
    _db->err(rc, "");

    return 0;
  }

  _dbname = dbname;

  return 1;
}

unsigned int
SP_Database::slurp_to_cache(int min_examples)
{
  Dbc* cursor = nullptr;

  int rc = _db->cursor(NULL, &cursor, 0);
  if (0 != rc) {
    _db->err(rc, "cannot acquire cursor");
    return 0;
  }

  int records_in_database = 0;

  Dbt zkey, zdata;

  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    records_in_database++;

    if (sizeof(Count_Radius) != zdata.get_size()) {
      continue;
    }

    const Count_Radius* cr = reinterpret_cast<const Count_Radius*>(zdata.get_data());

    if (cr->_count < min_examples) {
      continue;
    }

    DBKey* dbkey = reinterpret_cast<DBKey*>(zkey.get_data());

    Count_Radius cr2;
    cr2._count = cr->_count;
    cr2._radius = cr->_radius;
    _hash[*dbkey] = cr2;
  }

  if (DB_NOTFOUND != rc) {
    _db->err(rc, "Strange error at end of cursor");
  }

  if (verbose) {
    cerr << "Slurped " << _hash.size() << " of " << records_in_database
         << " database records into hash "
         << static_cast<float>(_hash.size()) / static_cast<float>(records_in_database)
         << '\n';
  }

  return _hash.size();
}

int
SP_Database::_parse_database_record(
    const const_IWSubstring& fromdb,
    const DBKey& dbkey,  // just for informational messages
    int& radius, int& count, IWString& example)
{
  if (fromdb.nwords() < 3) {
    cerr << "SP_Database::_parse_database_record:corrupted database contents for bit "
         << dbkey._bit << " found '" << fromdb << "'\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  fromdb.nextword(token, i);

  if (!token.numeric_value(count) || count < 1) {
    cerr << "Invalid count in database for bit " << dbkey._bit << " found '" << fromdb
         << "'\n";
    return 0;
  }

  fromdb.nextword(token, i);

  if ('4' == token) {
    radius = 4;
  } else if ('3' == token) {
    radius = 3;
  } else if ('2' == token) {
    radius = 2;
  } else if ('1' == token) {
    radius = 1;
  } else if ('0' == token) {
    radius = 0;
  } else if (!token.numeric_value(radius) || radius < 0)  // Strange
  {
    cerr << "Invalid radius in database for bit " << dbkey._bit << " found '" << fromdb
         << "'\n";
    return 0;
  }

  // What follows is the smiles and name of the example

  example = fromdb.substr(i + 1);

  // cerr << "Example set to '" << example << "'\n";

  return 1;
}

int
SP_Database::_in_cache_mutex_protected(Dbt& zkey, const DBKey& dbkey, int& radius,
                                       int& count)
{
  // const DBKey * dbkey = reinterpret_cast<const DBKey *>(zkey.get_data());

  const auto f = _hash.find(dbkey);

  if (f == _hash.end()) {
    return 0;
  }

  _cache_hit++;

  count = (*f).second._count;
  radius = (*f).second._radius;
  _bits_found++;
  return 1;
}

int
SP_Database::_in_cache(Dbt& zkey, const DBKey& dbkey, int& radius, int& count)
{
  _mutex.lock();
  int rc = _in_cache_mutex_protected(zkey, dbkey, radius, count);
  _mutex.unlock();
  return rc;
}

/*
  The DBKey object is not really used, just for informational messages
*/

int
SP_Database::do_lookup(Dbt& zkey, const DBKey& dbkey, int& radius, int& count,
                       IWString& example)
{
// #define DEBUG_SP_DATABASE_DO_LOOKUP
#ifdef DEBUG_SP_DATABASE_DO_LOOKUP
  cerr << "SP_Database::do_lookup, bit " << dbkey._bit << " examples? "
       << _contains_examples << '\n';
#endif

  example.resize_keep_storage(0);

  _lookups_done++;

  if (max_hash_size > 0 && _in_cache(zkey, dbkey, radius, count)) {
    return 1;
  }

  if (slurp_examples > 0) {  // no longer doing lookups, bit not present
    return 0;
  }

  Dbt zdata;
  char* user_memory;
  if (!_contains_examples) {
    user_memory = new char[8];
    zdata.set_ulen(8);
  } else {
    user_memory = new char[1024];
    zdata.set_ulen(1024);
  }

  std::unique_ptr<char[]> free_user_memory(user_memory);

  zdata.set_data(user_memory);
  zdata.set_flags(DB_DBT_USERMEM);

  //_logfile << "key " << zkey.get_size() << " bytes, flags " << zkey.get_flags() << '\n';

  int dbrc = _db->get(NULL, &zkey, &zdata, 0);

  if (0 == dbrc) {  // great, found
    ;
  } else if (DB_NOTFOUND == dbrc) {  // not here
    return 0;
  } else {
    _db->err(dbrc, "Unspecified database error, see Ian");
    return 0;
  }

  _bits_found++;

#ifdef DEBUG_SP_DATABASE_DO_LOOKUP
  cerr << "Db contents size " << zdata.get_size() << '\n';
#endif

  if (8 == zdata.get_size())  // no examples here. Should be a #defined symbol
  {
    const int* iptr = reinterpret_cast<const int*>(zdata.get_data());

    count = iptr[0];
    radius = iptr[1];
    if (_logfile_open) {
      _mutex.lock();
      _logfile << "bit " << dbkey._bit << " centre " << dbkey._acca << " found " << count
               << " count, radius " << radius << '\n';
      _mutex.unlock();
    }
  } else {
    const_IWSubstring fromdb(reinterpret_cast<const char*>(zdata.get_data()),
                             zdata.get_size());

    if (!_parse_database_record(fromdb, dbkey, radius, count, example)) {
      cerr << "Invalid database contents, ignored\n";
      return 0;
    }
    //  cerr << "Example is '" << example << "'\n";
  }

  // slightly dangerous access to unguarded variable _hash_size
  if (_hash_size < static_cast<int>(max_hash_size)) {
    DBKey* dbkey = reinterpret_cast<DBKey*>(zkey.get_data());
    Count_Radius cr;
    cr._count = count;
    cr._radius = radius;

    _mutex.lock();
    _hash[*dbkey] = cr;
    _hash_size++;
    _mutex.unlock();
  }

  return 1;
}

int
SP_Database::check_exists(Dbt* zkey)
{
  return 0 == _db->exists(NULL, zkey, 0);
}

int
SP_Database::get(DbTxn* env, Dbt& dkey, Dbt& zdata, u_int32_t flags)
{
  int rc = _db->get(env, &dkey, &zdata, flags);

  if (0 == rc) {
    return 1;
  }

  if (DB_NOTFOUND) {
    return 0;
  }

  _db->err(rc, "Trying to retrieve ");
  cerr.write(reinterpret_cast<const char*>(dkey.get_data()), dkey.get_size());
  cerr << '\n';

  return 0;
}

int
SP_Database::determine_if_examples_present()
{
  Dbc* cursor = nullptr;

  int rc = _db->cursor(NULL, &cursor, 0);
  if (0 != rc) {
    _db->err(rc, "SP_Database::determine_if_examples_present:cannot acquire cursor");
    return 0;
  }

  Dbt zkey, zdata;

  _contains_examples = 0;

  for (int records_checked = 0; 0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT));
       records_checked++) {
    if (records_checked > 100) {
      break;
    }

    //  cerr << records_checked << " got " << zdata.get_size() << " bytes\n";

    if (8 == zdata.get_size()) {  // most common case
      continue;
    }

    if (zdata.get_size() > 8) {
      _contains_examples = 1;
      break;
    }
  }

  cursor->close();

  return _contains_examples;
}

Set_of_Databases::Set_of_Databases() : _dbenv(0u)
{
  _user_wants_example_structures = 1;

  return;
}

Set_of_Databases::~Set_of_Databases()
{
  return;
}

int
Set_of_Databases::debug_print(std::ostream& output) const
{
  output << "Set_of_Databases::debug_print:set contains " << _db.size() << " databases\n";
  cerr << " _user_wants_example_structures " << _user_wants_example_structures << '\n';
  for (const SP_Database* db : _db) {
    db->debug_print(output);
  }

  return 1;
}

int
Set_of_Databases::report(std::ostream& os) const
{
  os << "Report on " << _db.size() << " databases\n";

  for (const SP_Database* db : _db) {
    db->report(os);
  }

  return 1;
}

int
Set_of_Databases::open_logfile(const char* stem)
{
  IWString s(stem);

  for (int i = 0; i < _db.number_elements(); ++i) {
    IWString fname(stem);
    fname << '.' << i << ".log";
    if (! _db[i]->open_logfile(fname.null_terminated_chars())) {
      cerr << "Set_of_Databases::open_logfile:cannot set logfile '" << fname << "'\n";
      return 0;
    }
  }

  return _db.number_elements();
}

int
Set_of_Databases::build(Command_Line& cl, char flag_env, char flag_db, const int verbose)
{
  if (!cl.option_present(flag_env)) {
    return build(cl, flag_db, verbose, nullptr);
  }

  const char* e = cl.option_value(flag_env);

  u_int32_t env_flags = DB_CREATE | DB_INIT_MPOOL;

  int rc = _dbenv.open(e, env_flags, 0);

  if (0 != rc) {
    cerr << "Set_of_Databases::build:cannot initialise DB environment '" << e << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Using database environment '" << e << "'\n";
  }

  return build(cl, flag_db, verbose, &_dbenv);
}

int
Set_of_Databases::build(const Command_Line& cl, const char flag_db, const int verbose,
                        DbEnv* envptr) {
  if (! cl.option_present(flag_db)) {
    return 0;
  }

  IWString fname;
  for (int i = 0; cl.value(flag_db, fname, i); ++i) {
    if (verbose) {
      cerr << "Opening database '" << fname << "'" << '\n';
    }
    
    std::unique_ptr<SP_Database> db = std::make_unique<SP_Database>();

    if (nullptr != envptr) {
      db->set_env(envptr);
    }

    if (!db->open(fname.null_terminated_chars(), 0 /* no threading */)) {
      cerr << "Set_of_Databases::build:cannot open '" << fname << "'\n";
      return 0;
    }

    _db << db.release();
  }

  if (verbose) {
    cerr << "Set_of_Databases::build:opened " << _db.size() << " databases\n";
  }

  if (_user_wants_example_structures) {
    _determine_if_examples_present();
  }

  return _db.number_elements();
}

int
Set_of_Databases::AddDatabase(const std::string& dbname) {
  std::unique_ptr<SP_Database> db = std::make_unique<SP_Database>();

  // No allowances for database environments here.

  static constexpr int kNoThreading = 0;
  if (!db->open(dbname.c_str(), kNoThreading)) {
    cerr << "Set_of_Databases::AddDatabase:cannot open '" << dbname << "'\n";
    return 0;
  }

  _db << db.release();

  if (_user_wants_example_structures) {
    _determine_if_examples_present();
  }

  return _db.number_elements();
}

int
Set_of_Databases::slurp_to_cache(int min_examples) {
  for (int i = 0; i < _db.number_elements(); ++i) {
    if (!_db[i]->slurp_to_cache(min_examples)) {
      cerr << "Set_of_Databases::slurp_to_cache:cannot slurp database " << i << " for "
           << min_examples << " examples\n";
    }
  }

  return _db.number_elements();
}

int
Set_of_Databases::_determine_if_examples_present()
{
  for (SP_Database* db : _db) {
    db->determine_if_examples_present();

    if (!db->contains_example_structures()) {
      _user_wants_example_structures = 0;
    }
  }

  // cerr <<
  // "Set_of_Databases::_determine_if_examples_present:_user_wants_example_structures " <<
  // _user_wants_example_structures << '\n';

  return 1;
}

// #define DEBUG_BIT_LOOKUP

int
Set_of_Databases::lookup_bit(const DBKey& dbkey, int& count, IWString& example_structure)
{
#ifdef DEBUG_BIT_LOOKUP
  cerr << "Set_of_Databases::lookup_bit:looking for bit " << dbkey << " across " << _db.size()
       << " databases\n";
#endif

  int nexamples = 0;
  example_structure.resize_keep_storage(0);

  Dbt zkey((void*)(&dbkey), sizeof(dbkey));  // loss of const OK

  IWString s;  // scope here for efficiency
  for (SP_Database* db : _db) {
    int radius;

    if (!db->do_lookup(zkey, dbkey, radius, count, s)) {
      //    cerr << "Did not find bit " << dbkey._bit << " (radius " <<
      //    static_cast<int>(dbkey._radius) << ") in database " << i << '\n';
      continue;
    }

    if (radius != dbkey._radius) {  // collision
      continue;
    }

#ifdef DEBUG_BIT_LOOKUP
    cerr << "Found bit " << dbkey._bit << " (radius " << radius << ") in database " << i
         << ", count " << count << " example '" << s << "'\n";
#endif

    if (_subtraction.active()) {
      _subtraction.adjust_count_for_subtraction(dbkey, count);
      if (0 == count) {
        continue;
      }
    }

    nexamples += count;

    if (_user_wants_example_structures && example_structure.empty()) {
      example_structure = s;
    }
  }

  return nexamples;
}

int
Set_of_Databases::lookup_bit_threaded(const DBKey& dbkey)
{
  if (1 == _db.size()) {  // no need for multi threaded execution
    int notused1;
    IWString notused2;
    return lookup_bit(dbkey, notused1, notused2);
  }

  int nexamples = 0;

  for (SP_Database* db : _db) {
    int radius, count;
    IWString s;

    DBKey tmp(dbkey);
    // loss of const OK. Each thread needs its own key
    Dbt zkey((void*)(&tmp), sizeof(tmp));
    assert(tmp._bit == dbkey._bit);

    if (!db->do_lookup(zkey, dbkey, radius, count, s)) {
      continue;
    }

    if (radius != dbkey._radius) {  // collision
      continue;
    }

#ifdef DEBUG_BIT_LOOKUP
    cerr << "Found bit " << dbkey._bit << " (radius " << radius << ") in database " << i
         << ", count " << count << '\n';
#endif

    if (_subtraction.active()) {
      _subtraction.adjust_count_for_subtraction_threaded(dbkey, count);  // count may be zero
    }

    nexamples += count;
  }

  return nexamples;
}

int
Set_of_Databases::lookup_bit(const DBKey& dbkey, resizable_array_p<Count_Example>& ce)
{
  int count;
  IWString zexample;
  const auto rc = lookup_bit(dbkey, count, zexample);

  if (0 == rc) {
    return 0;
  }

  ce.add(new Count_Example(zexample, dbkey._radius, count));

  return rc;
}

// Determine the atom typing in use across the databases and
// set both arguments.
int
Set_of_Databases::determine_atom_typing_in_use(IWString& string_atype, int& iwecfp_atom_type)
{
  Dbt dkey((void*)(ATYPE_KEY), ::strlen(ATYPE_KEY));
  Dbt d;

  char buffer[32];
  d.set_data(buffer);
  d.set_ulen(sizeof(buffer));
  d.set_flags(DB_DBT_USERMEM);

  for (int i = 0; i < _db.number_elements(); ++i) {
    SP_Database* db = _db[i];

    if (!db->get(NULL, dkey, d, 0)) {
      cerr << "Database " << i << " lacks atom type key '" << ATYPE_KEY << "'\n";
      return 0;
    }

    const_IWSubstring tmp(reinterpret_cast<const char*>(d.get_data()), d.get_size());

    //  cerr << "type '" << tmp << "' fetched\n";

    if (0 == i) {
      string_atype = tmp;
    } else if (string_atype == tmp) {  // great, consistent
      ;
    } else {
      cerr << "INconsistent atom types in use, stored " << tmp << "', using '"
           << string_atype << "', cannot continue\n";
      return 0;
    }
  }

  iwecfp_atom_type = determine_atom_type(string_atype);

  if (0 == iwecfp_atom_type) {
    cerr << "INvalid atom type stored in database\n";
    return 0;
  }

  if (verbose) {
    cerr << "Atom type '" << string_atype << "', numeric " << iwecfp_atom_type << '\n';
  }

  return 1;
}

/*
  The max search radius will be the shortest radius stored in any of our databases
*/

int
Set_of_Databases::determine_max_search_radius(int& min_radius)
{
  min_radius = 99;

  Dbt dkey((void*)(RADIUS_KEY), ::strlen(RADIUS_KEY));
  Dbt d;

  char buffer[32];
  d.set_data(buffer);
  d.set_ulen(sizeof(buffer));
  d.set_flags(DB_DBT_USERMEM);

  for (int i = 0; i < _db.number_elements(); ++i) {
    if (!_db[i]->get(NULL, dkey, d, 0)) {
      cerr << "Database " << i << " lacks atom type key '" << RADIUS_KEY << "'\n";
      continue;
    }

    const_IWSubstring tmp(reinterpret_cast<const char*>(d.get_data()), d.get_size());

    //  cerr << "type '" << tmp << "' fetched\n";

    int r;
    if (!tmp.numeric_value(r) || r < 1) {
      cerr << "Set_of_Databases::determine_max_radius:invalid '" << RADIUS_KEY
           << "' value '" << tmp << "', ignored\n";
      continue;
    }

    if (r < min_radius) {
      min_radius = r;
    }
  }

  if (verbose) {
    cerr << "Set_of_Databases::determine_max_radius:across " << _db.size()
         << " databases, min radius " << min_radius << '\n';
  }

  return 1;
}

int
Set_of_Databases::ReportSubtractDatabaseStatus(std::ostream& output) const {
  output << "Read " << _subtraction.nbits() << " bits produced by " <<
            _subtraction.nmolecules() << " molecules for subtraction set\n";
  return 1;
}

int
Set_of_Databases::ReportSubtractionDatabaseImpact(std::ostream& output) const {
  return _subtraction.report(output);
}

Subtraction_Set::Subtraction_Set()
{
  _molecules_in_subtraction_set = 0;

  _subtraction_bits_hit = 0;

  return;
}

int
Subtraction_Set::report(std::ostream& output) const
{
  if (_subtract.empty()) {
    output << "Subtraction_Set not active\n";
    return 1;
  }

  output << "Subtraction_Set::report:contains " << _subtract.size() << " bits, "
     << _subtraction_bits_hit << " adjustments made\n";

  return 1;
}

int
Subtraction_Set::build(const Command_Line& cl, const char flag,
                       Fingerprint_Characteristics& fpc,
                       Preprocessor preprocess_molecule) {
  IWString fname;
  for (int i = 0; cl.value(flag, fname, i); ++i) {
    if (!build(fname, fpc, preprocess_molecule)) {
      cerr << "Subtraction_Set::build:cannot build from '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Subtraction_Set::build(data_source_and_type<Molecule>& input,
                       Fingerprint_Characteristics& fpc,
                       Preprocessor preprocess_molecule) {
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    _molecules_in_subtraction_set++;

    preprocess_molecule(*m);

    std::unique_ptr<Molecule> free_m(m);

    Set_of_Bits<iwecfp_database::Bit_Produced> sob;

    compute_fingerprints(*m, fpc, sob);

    _transfer_fingerprints_to_subtract(sob);
  }

  return _subtract.size();
}

int
Subtraction_Set::build(IWString& fname, Fingerprint_Characteristics& fpc,
                Preprocessor preprocess_molecule)
{
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);

  if (!input.good()) {
    cerr << "build_subtraction_set:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input, fpc, preprocess_molecule);
}

int
Subtraction_Set::_transfer_fingerprints_to_subtract(const Set_of_Bits<Bit_Produced>& sob)
{
  for (auto i : sob) {
    const auto& k = i.first;

    //  unordered_map<DBKey, int, IWdbkeyHash>::const_iterator f = sub.find(k);
    const auto f = _subtract.find(k);

    //  cerr << "Loaded subtraction for " << k._bit << '\n';

    if (f == _subtract.end()) {
      _subtract[k] = 1;
    } else {
      _subtract[k]++;
    }
  }

  return 1;
}

int
Set_of_Databases::BuildSubtractionDatabase(IWString& fname,
                iwecfp_database::Fingerprint_Characteristics& fpc,
                int(*preprocess_molecule)(Molecule& m)) {
  return _subtraction.build(fname, fpc, preprocess_molecule);
}

int
Subtraction_Set::_adjust_count_for_subtraction(const DBKey& dbkey, int& count) const
{
  Subtract::const_iterator f = _subtract.find(dbkey);

  if (f == _subtract.end()) {
    return 0;
  }

  int s = (*f).second;

  if (count >= s) {
    count = count - s;
  } else {
    count = 0;
  }

  return 1;
}

int
Subtraction_Set::adjust_count_for_subtraction_threaded(const DBKey& dbkey, int& count)
{
  return _adjust_count_for_subtraction(dbkey, count);
}

int
Subtraction_Set::adjust_count_for_subtraction(const DBKey& dbkey, int& count)
{
  const auto rc = _adjust_count_for_subtraction(dbkey, count);
  if (rc) {
    _subtraction_bits_hit++;
  }

  return rc;
}

Count_Example::Count_Example(const IWString& s, int r, int c)
    : _number_examples(c), _radius(r)
{
  if (s.nwords() < 2) {
    cerr << "Count_Example::Count_Example:invalid db contents '" << s << "'\n";
    return;
  }

  int i = 0;

  s.nextword(_first_example_smiles, i);
  s.nextword(_first_example_name, i);

  return;
}

int
Count_Example::do_output(IWString_and_File_Descriptor& output) const
{
  if (0 == _first_example_smiles.length()) {
    return 1;
  }

  output << _first_example_smiles << ' ';

  output << _first_example_name << ' ';

  if (_radius >= 0) {
    output << _radius << ' ';
  }

  output << _number_examples << " EX\n";

  return 1;
}

int
SP_Set_of_Databases::AddDatabase(const std::string& dbname) {
  if (! _dbs.AddDatabase(dbname)) {
    return 0;
  }

  int iwecfp_atom_type;
  IWString string_atype;
  if (! _dbs.determine_atom_typing_in_use(string_atype, iwecfp_atom_type)) {
    cerr << "SP_Set_of_Databases::AddDatabase:cannot determine atom typing in use\n";
    return 0;
  }

  _fingerprint_characteristics.set_atype(iwecfp_atom_type);

  if (_atom_typing_specification.active()) {
    return 1;
  }

  if (! _atom_typing_specification.build(string_atype)) {
    cerr << "SP_Set_of_Databases:AddDatabase:cannot initialise atom typing " <<
            string_atype << '\n';
    return 0;
  }

  int mr = 0;
  _dbs.determine_max_search_radius(mr);

  _fingerprint_characteristics.set_max_shell_radius(mr);

  return 1;
}

int
SP_Set_of_Databases::PerShellData(Molecule& m, std::vector<int>& result) {
  result.resize(_fingerprint_characteristics.max_shell_radius() + 1);
  std::fill(result.begin(), result.end(), std::numeric_limits<int>::max());
  cerr <<"Examining " << m.aromatic_smiles() << '\n';

  Set_of_Bits<iwecfp_database::Bit_Produced> sob;
  compute_fingerprints(m, _fingerprint_characteristics, sob);

  for (auto i : sob) {
    const auto& k = i.first;

    const int r = k._radius;

    const int nexamples = _dbs.lookup_bit_threaded(k);
    cerr << "Bit at radius " << r << " with " << nexamples << " examples\n";

    if (nexamples < result[r]) {
      result[r] = nexamples;
    }
  }

  for (int& c: result) {
    cerr << " value " << c << '\n';
    if (c == std::numeric_limits<int>::max()) {
      c = 0;
    }
  }

  return 1;
}

}  // namespace iwecfp_database_lookup
