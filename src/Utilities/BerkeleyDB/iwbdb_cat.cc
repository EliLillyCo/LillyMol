/*
  Concatenate multiple Berkeley databases into one
*/

#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "db_cxx.h"

#include "re2/re2.h"
#include <sys/stat.h>

using std::cerr;

static int verbose = 0;

static int store_flag = DB_NOOVERWRITE;

static DBTYPE dbtype = DB_BTREE;

static int duplicate_entries = 0;

static uint64_t items_read = 0;

static uint64_t items_stored = 0;

static uint64_t append_data_to_stored_data = 0;

static Report_Progress report_progress;

/*
  When loading descriptor databases, we can impose conditions on the number
  of tokens per record
*/

static int count_records_to_determine_tokens_per_record = 0;
static int count_records_to_determine_bytes_per_record = 0;

static int tokens_per_record = 0;

static unsigned int bytes_per_record = 0;

static int records_suppressed_by_size_or_token_mismatch = 0;

/*
  We have several kinds of databases where there is a count followed
  by some number of identifiers. When these are concatenated, we must
  sum the first tokens.
  If this value is 1, we just keep the first identifier encountered,
  If this value is 2, we concatenate all the identifiers
*/

static int special_processing_for_dicer_type_databases = 0;

static int special_processing_for_synthetic_precedent_databases = 0;

// From the -N option
static int numeric_value_column = -1;

static int store_items_not_matching_pattern = 0;

static int dicer_items_summed = 0;

// Interesting to keep statistics on the numbers encountered in the count field

static Accumulator_Int<int> count_stats;

/*
  When concatenating new data onto the end of old data, we may have
  a separator
*/

static IWString separator(':');

static int items_appended_to_existing_data = 0;

/*
  Sometimes we are trying to update databases in use. We can optionally
  wait until we get a write lock
*/

static int sleep_before_next_attempt = 0;

static int max_open_attempts = 0;

static std::unique_ptr<re2::RE2> data_rx_must_match, data_rx_non_match;

static int skipped_for_matching_data_rx = 0;
static int skipped_for_not_matching_data_rx = 0;

/*
  We can speed things up a lot by doing a background wc
  as we open each database
*/

static int perform_wc = 0;

/*
  Sometimes when combining databases we insist that a certain key be
  common to the 2 databases.
  I included this for dealing with the case of descriptor databases,
  and ensuring that the _HEADER records matched
*/

class Key_to_Match {
 private:
  IWString _dkey;

  Dbt _key_datum;

  int _value_known;

  IWString _data;

 public:
  Key_to_Match();

  const IWString&
  dkey() const {
    return _dkey;
  }

  const Dbt&
  key_datum() const {
    return _key_datum;
  }

  void
  set_key(const const_IWSubstring& k);

  int
  matches_key_in(Db& db);
};

Key_to_Match::Key_to_Match() {
  _value_known = 0;

  return;
}

void
Key_to_Match::set_key(const const_IWSubstring& k) {
  _dkey = k;

  _key_datum.set_data(const_cast<char*>(k.rawchars()));
  _key_datum.set_size(k.length());

  return;
}

int
Key_to_Match::matches_key_in(Db& db) {
  Dbt fromdb;

  int rc = db.get(NULL, &_key_datum, &fromdb, 0);

  if (DB_NOTFOUND == rc) {
    return 1;
  }

  if (0 == rc) {
    ;
  } else if (DB_NOTFOUND == rc) {
    ;
  } else {
    cerr << "Key_to_Match::matches_key_in: error fetching '" << _dkey << "' ";
    db.err(rc, "");
    return 0;
  }

  // If nothing is known and there is nothing in the database, we don't know anything
  // Note that this leaves open the bug possibility of joining 3 databases, the
  // first of which doesn't have the key and the two others have different values...

  if (DB_NOTFOUND == rc && !_value_known) {
    return 1;
  }

  if (!_value_known) {
    _data.strncpy(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());

    _value_known = 1;

    return 1;
  }

  // The value is known

  if (DB_NOTFOUND == rc) {
    cerr << "Key '" << _dkey << "' mismatch, stored '" << _data
         << "' but missing from new database\n";
    return 0;
  }

  IWString tmp;
  tmp.strncpy(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());

  if (tmp == _data) {
    if (verbose > 1) {
      cerr << " key '" << _dkey << "' matches\n";
    }

    return 1;
  }

  cerr << "Key '" << _dkey << "' mismatch: stored '" << _data << "' got '" << tmp
       << "'\n";
  return 0;
}

static Key_to_Match* keys_to_match = NULL;
static int number_keys_to_match = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << DB_VERSION_STRING << '\n';

  // clang-format off
  cerr << "Combines multiple Berkeley DB databases into one\n";
  cerr << '\n';

  cerr << "usage: iwbdb_cat target database1 database2 database3 ... \n";
  cerr << "        the contents of database1, database2, ... are put into TARGET\n";
  cerr << " -d <target>  specify target database via different syntax\n";
  cerr << " -o           overwrite existing entries in TARGET\n";
//cerr << " -c <size>    value for GDBM_CACHESIZE (sometimes smaller is better here)\n";
  cerr << " -M <key>     key(s) whose value must be consistent in all databases\n";
  cerr << "              note this is done sequentially\n";
  cerr << " -a           append new data to existing data\n";
  cerr << " -S <string>  append and use <string> as the separator\n";
  cerr << " -r <number>  report progress every <number> items\n";
  cerr << " -L <seconds> if the database is write-locked, wait <seconds> seconds\n";
  cerr << "              and try again\n";
  cerr << " -L maxtries=nn  the maximum number of times to try to get a write-lock\n";
  cerr << " -C <number>  only load records that have <number> tokens\n";
  cerr << " -C count     read records from the existing database to determine how\n";
  cerr << "              many tokens per record\n";
  cerr << " -p           special processing for databases with 'count ...', increments count\n";
  cerr << " -N <col>     values have a numeric value in column <col>, increment\n";
  cerr << " -q           special processing for synthetic precedent databases (two int's)\n";
  cerr << " -b           store items that do not match the dicer pattern\n";
  cerr << " -y <type>    access type (btree, hash, recno, queue)\n";
  cerr << " -x <rx>      only process data values that        match <rx> (use FILE:<rx> if needed)\n";
  cerr << " -X <rx>      only process data values that do not match <rx>\n";
  cerr << " -v           verbose output\n";
  // clang-format on

  exit(rc);
}

static int
check_key_matches(Db& db) {
  for (int i = 0; i < number_keys_to_match; i++) {
    Key_to_Match& k = keys_to_match[i];

    if (!k.matches_key_in(db)) {
      cerr << "Key mismatch\n";
      return 0;
    }
  }

  return 1;
}

static int
_matches_tokens_or_bytes(const Dbt& zdata) {
  if (bytes_per_record) {
    return zdata.get_size() == bytes_per_record;
  }

  if (tokens_per_record) {
    const_IWSubstring d(reinterpret_cast<const char*>(zdata.get_data()),
                        zdata.get_size());

    return tokens_per_record == d.nwords();
  }

  abort();  // should never come to here

  return 0;
}

static int
matches_tokens_or_bytes(const Dbt& zdata) {
  if (_matches_tokens_or_bytes(zdata)) {  // great
    return 1;
  }

  records_suppressed_by_size_or_token_mismatch++;

  if (0 == verbose) {  // don't say anything
    ;
  } else if (bytes_per_record) {
    cerr << "Record size mismatch " << zdata.get_size() << " vs " << bytes_per_record
         << '\n';
  } else if (tokens_per_record) {
    const_IWSubstring d(reinterpret_cast<const char*>(zdata.get_data()),
                        zdata.get_size());

    cerr << "Token count mismatch " << d.nwords() << " vs " << tokens_per_record << '\n';
    cerr << "Data '" << d << "'\n";
  }

  return 0;
}

static int
find_max(const extending_resizable_array<int>& era) {
  int zmax = 0;   // the actual maximum value
  int imax = -1;  // the array index that has the max value - we return it
  int nmax = 0;   // the number of instances ofthe maximum value

  for (int i = 0; i < era.number_elements(); i++) {
    if (era[i] == zmax) {
      nmax++;
    } else if (era[i] > zmax) {
      zmax = era[i];
      imax = i;
    }
  }

  if (1 == nmax) {
    return imax;
  }

  cerr << "Warning, " << nmax << " records each had " << zmax
       << " items, possibly ambiguous sizing\n";

  return imax;
}

/*
  We will only be loading records that have the same number of tokens,
  or are of the same size as records currently in the database.
  Fetch some records from the database to determine what these
  values should be
*/

static int
determine_tokens_per_record_in_target_database(Db& target) {
  Dbc* cursor;

  int rc = target.cursor(NULL, &cursor, 0);
  if (0 != rc) {
    target.err(rc, "cannot acquire cursor");
    return 0;
  }

  Dbt zdata, zkey;

  int nfetch = 0;

  extending_resizable_array<int> tokens(1000);
  extending_resizable_array<int> record_length(1000);

  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    record_length[zdata.get_size()]++;

    const_IWSubstring d(reinterpret_cast<const char*>(zdata.get_data()),
                        zdata.get_size());

    int t = d.nwords();

    tokens[t]++;

    nfetch++;

    if (nfetch > 100) {
      break;
    }
  }

  cursor->close();

  if (count_records_to_determine_bytes_per_record) {
    bytes_per_record = find_max(record_length);
  } else if (count_records_to_determine_tokens_per_record) {
    tokens_per_record = find_max(tokens);
  } else {
    abort();  // should not have been called
  }

  return 1;
}

static int
do_the_store(Db& target, Dbt& dkey, Dbt& zdata) {
  int rc = target.put(NULL, &dkey, &zdata, store_flag);

  if (0 == rc) {
    items_stored++;
    return 1;
  }

  if (DB_KEYEXIST == rc) {
    if (verbose > 2) {
      cerr.write((const char*)dkey.get_data(), dkey.get_size());
      cerr << " already in target database\n";
    }

    duplicate_entries++;

    return 1;
  }

  target.err(rc, "Fatal error on store");
  return 0;
}

static int
determine_existing_count(IWString& zdata, int& count, int do_shift) {
  count = 0;

  int n = zdata.length();

  const char* c = zdata.rawchars();

  for (int i = 0; i < n; i++) {
    if (isdigit(c[i])) {
      count = 10 * count + c[i] - '0';
    } else if (' ' == c[i]) {
      if (do_shift) {
        zdata.remove_leading_chars(i + 1);
      }
      return 1;
    }
  }

  cerr << "Cannot extract count and identifier '" << zdata << "'\n";
  return 0;
}

// These special purpose processings should be updated to use protos.

struct Count_Radius {
  int _count;
  int _radius;
};

typedef struct Count_Radius Count_Radius;

static int
do_special_processing_for_synthetic_precedent_databases(Db& target, Dbt& dkey,
                                                        Dbt& zdata,
                                                        Dbt& fromdb) {
  if (sizeof(Count_Radius) != fromdb.get_size()) {
    cerr << "Incorrect size for stored data " << fromdb.get_size() << " bytes, ignored\n";
    return do_the_store(target, dkey, zdata);
  }

  items_appended_to_existing_data++;

  const Count_Radius* cr1 = reinterpret_cast<const Count_Radius*>(fromdb.get_data());
  const Count_Radius* cr2 = reinterpret_cast<const Count_Radius*>(zdata.get_data());

  if (cr1->_radius != cr2->_radius) {
    cerr << "do_special_processing_for_synthetic_precedent_databases:radius mismatch "
         << cr1->_radius << " vs " << cr2->_radius << " ignored!\n";
    return do_the_store(target, dkey, zdata);
  }

  Count_Radius crs;
  crs._count = cr1->_count + cr2->_count;
  crs._radius = cr1->_radius;

  Dbt to_store(&crs, sizeof(crs));

  return do_the_store(target, dkey, to_store);
}

// We know that `dbkey` is not in the database.
// `fromdb` is not const because we re-use it to store.
static int
do_special_processing_for_dicer_type_databases(Db& target, Dbt& dkey, Dbt& zdata,
                Dbt& fromdb) {
  items_appended_to_existing_data++;

  IWString already_stored(reinterpret_cast<const char*>(fromdb.get_data()),
                          fromdb.get_size());

  int existing_count;
  if (determine_existing_count(already_stored, existing_count, 1)) {
    ;
  } else if (store_items_not_matching_pattern) {
    return do_the_store(target, dkey, fromdb);
  } else {
    cerr << "Cannot determine count '" << already_stored << "', ignored\n";
    return 1;
  }

  // cerr << "STORED: from '" << already_stored << "' count is " << existing_count <<
  // '\n';

  IWString newdata(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());

  int strip_flag = 0;
  if (special_processing_for_dicer_type_databases > 1) {
    strip_flag = 1;
  }

  int new_count;
  if (!determine_existing_count(newdata, new_count, strip_flag)) {
    cerr << "Cannot determine count '" << already_stored << "', ignored\n";
    return 1;
  }

  // cerr << "NEWDATA: from '" << newdata << "' count is " << new_count << '\n';

  IWString to_store;

  if (verbose) {
    count_stats.extra(existing_count + new_count);
  }

  to_store << (existing_count + new_count) << ' ' << already_stored;
  if (special_processing_for_dicer_type_databases > 1) {
    to_store << separator << newdata;
  }

  dicer_items_summed++;

  // re-use fromdb

  fromdb.set_data((void*)(to_store.rawchars()));  // horrible old style cast
  fromdb.set_size(to_store.length());

  return do_the_store(target, dkey, fromdb);
}

std::optional<int>
GetCount(const const_IWSubstring& buffer, int column) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i); ++col) {
    if (col != column) {
      continue;
    }

    int v;
    if (! token.numeric_value(v)) {
      return std::nullopt;
    }
    return v;
  }

  return std::nullopt;
}

static int
do_special_processing_for_numeric_column(Db& target,
                Dbt& dkey,
                Dbt& zdata,
                int numeric_value_column,
                Dbt& fromdb) {
  const_IWSubstring buffer((char*) zdata.get_data(), zdata.get_size());
  std::optional<int> extra = GetCount(buffer, numeric_value_column);
  if (! extra) {
    cerr << "do_special_processing_for_numeric_column:invalid numeric '" << buffer << "'\n";
    return do_the_store(target, dkey, zdata);
  }

  buffer.set((char*)fromdb.get_data(), static_cast<size_t>(fromdb.get_size()));
  std::optional<int> initial_count = GetCount(buffer, numeric_value_column);
  if (! extra) {
    cerr << "do_special_processing_for_numeric_column:invalid DB count '" << buffer << "'\n";
    return do_the_store(target, dkey, zdata);
  }

  int new_value = *initial_count + *extra;
  IWString new_data;
  new_data.resize(zdata.get_size() + 3);  // 3 is arbitrary guess.

  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i); ++col) {
    if (col > 0) {
      new_data << ' ';
    }
    if (col == numeric_value_column) {
      new_data << new_value;
    } else {
      new_data << token;
    }
  }

  items_appended_to_existing_data++;

  fromdb.set_data((void*) new_data.data());
  fromdb.set_size(new_data.length());
  return do_the_store(target, dkey, fromdb);
}

static int
do_append_data_to_stored_data(Db& target, Dbt& dkey, Dbt& zdata) {

  Dbt already_stored;

  int rc = target.get(NULL, &dkey, &already_stored, 0);

  if (DB_NOTFOUND == rc) {
    return do_the_store(target, dkey, zdata);
  }

  if (special_processing_for_dicer_type_databases) {
    return do_special_processing_for_dicer_type_databases(target, dkey, zdata, already_stored);
  }

  if (special_processing_for_synthetic_precedent_databases) {
    return do_special_processing_for_synthetic_precedent_databases(target, dkey, zdata, already_stored);
  }

  if (numeric_value_column >= 0) {
    return do_special_processing_for_numeric_column(target, dkey, zdata, numeric_value_column, already_stored);
  }

  items_appended_to_existing_data++;

  IWString newdata;
  newdata.resize(already_stored.get_size() + separator.length() + zdata.get_size());

  newdata.strncat(reinterpret_cast<const char*>(already_stored.get_data()),
                  already_stored.get_size());
  newdata += separator;
  newdata.strncat(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());

  // re-use already_stored

  already_stored.set_data(const_cast<char*>(newdata.rawchars()));
  already_stored.set_size(newdata.length());

  return do_the_store(target, dkey, already_stored);
}

static int
store_or_append(Db& target, Dbt& dkey, Dbt& zdata) {
  if (bytes_per_record || tokens_per_record) {
    if (!matches_tokens_or_bytes(zdata)) {
      return 1;
    }
  }

  if (append_data_to_stored_data) {
    return do_append_data_to_stored_data(target, dkey, zdata);
  }

  return do_the_store(target, dkey, zdata);
}

static int
matches(re2::RE2& rx, const Dbt& d) {
  re2::StringPiece tmp(static_cast<const char*>(d.get_data()), d.get_size());

  return RE2::PartialMatch(tmp, rx);
}

/*
  We can often speed things up by scanning the target database for all
  the items in the second database. This forces things into RAM, so
  we hold the database write lock for a shorter period of time

  Implement this sometime
*/

/*static int
prefetch (GDBM_FILE & target, GDBM_FILE & source)
{
}*/

static int
iwbdb_cat(Db& target, Db& source) {
  Dbc* cursor;

  int rc = source.cursor(NULL, &cursor, 0);
  if (0 != rc) {
    source.err(rc, "cannot acquire cursor");
    return 0;
  }

  Dbt zdata, zkey;

  int items_fetched_this_database = 0;

  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    items_fetched_this_database++;
    items_read++;

    if (report_progress()) {
      cerr << "Read " << items_read << " items, stored " << items_stored << '\n';
    }

    if (data_rx_must_match && !matches(*data_rx_must_match, zdata)) {
      skipped_for_not_matching_data_rx++;
      continue;
    }

    if (data_rx_non_match && matches(*data_rx_non_match, zdata)) {
      skipped_for_matching_data_rx++;
      continue;
    }

    if (!store_or_append(target, zkey, zdata)) {
      return 0;
    }
  }

  if (verbose) {
    cerr << items_fetched_this_database << " items fetched, " << items_stored
         << " items stored\n";
  }

  cursor->close();

  return 1;
}

static int
iwbdb_cat(Db& target, const char* dbname);  // forward declaration

static int
iwbdb_cat_list_of_db(Db& target, const char* fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "iwbdb_cat_list_of_db:cannot open list of db's '" << fname << "'\n";
    return 0;
  }

  IWString buffer;

  while (input.next_record(buffer)) {
    if (!iwbdb_cat(target, buffer.null_terminated_chars())) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
iwbdb_cat(Db& target, const char* dbname) {
  if (0 == strncmp(dbname, "F:", 2)) {
    return iwbdb_cat_list_of_db(target, dbname);
  }

  Db source(0, DB_CXX_NO_EXCEPTIONS);

  int rc = source.open(NULL, dbname, NULL, DB_UNKNOWN, DB_RDONLY, 0);

  if (0 != rc) {
    cerr << "Cannot open source database '" << dbname << "' ";
    source.err(rc, "");
    return 0;
  }

  if (perform_wc) {
    if (0 == fork()) {
      IWString cmd;
      execl("/usr/bin/wc", dbname, NULL);
      cerr << "Yipesl execl for wc failed\n";
      _exit(2);
    }
  }

  if (verbose) {
    cerr << "Processing '" << dbname << "'\n";
  }

  if (number_keys_to_match) {
    if (!check_key_matches(source)) {
      cerr << "Required key mismatch in database '" << dbname << "'\n";
      return 0;
    }
  }

  rc = iwbdb_cat(target, source);

  source.close(0);

  return rc;
}

static int
setup_regular_expression_from_file(std::unique_ptr<re2::RE2>& rx,
                                   iwstring_data_source& input) {
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Empty regular expression file, cannot continue\n";
    return 0;
  }

  if (iwre2::RE2Reset(rx, buffer)) {
    return 1;
  }

  cerr << "Cannot initialise regular expression '" << buffer << "'\n";
  return 0;
}

static int
setup_regular_expression_from_file(std::unique_ptr<re2::RE2>& rx, const char* fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open regular expression file '" << fname << "'\n";
    return 0;
  }

  return setup_regular_expression_from_file(rx, input);
}

static int
setup_regular_expression(std::unique_ptr<re2::RE2>& rx, const char* s) {
  IWString tmp(s);

  if (tmp.starts_with("FILE:")) {
    tmp.remove_leading_chars(5);
    return setup_regular_expression_from_file(rx, tmp.null_terminated_chars());
  }

  if (iwre2::RE2Reset(rx, tmp)) {
    return 1;
  }

  cerr << "Cannot initialise regular expression '" << s << "'\n";
  return 0;
}

static int
open_target_database(const char* target_fname, Db& target) {
  int attempts_made = 0;

  int flags = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

  while (1) {
    int rc;

    if (dash_s(target_fname)) {
      dbtype = DB_UNKNOWN;
      rc = target.open(NULL, target_fname, NULL, dbtype, 0, flags);
    } else {
      rc = target.open(NULL, target_fname, NULL, dbtype, DB_CREATE, flags);
    }

    if (0 == rc) {
      if (verbose && attempts_made > 0) {
        cerr << "Obtained write lock\n";
      }

      return 1;
    }

    attempts_made++;

    if (max_open_attempts > 0 && attempts_made >= max_open_attempts) {
      return 0;
    }

    if (0 == sleep_before_next_attempt) {
      return 0;
    }

    if (verbose) {
      cerr << "Could not obtain write lock, sleeping for " << sleep_before_next_attempt
           << " seconds\n";
    }

    sleep(sleep_before_next_attempt);
  }

  assert(NULL == "Should not come to here");

  return 0;
}

static int
iwbdb_cat(int argc, char** argv) {
  Command_Line cl(argc, argv, "voc:M:L:aS:r:C:pqd:by:x:X:wN:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('o') && cl.option_present('a')) {
    cerr << "The overwrite (-o) and append (-a) options are mutually incompatible\n";
    usage(3);
  }

  if (cl.option_present('o')) {
    store_flag = 0;

    if (verbose) {
      cerr << "All entries in TARGET will be overwritten\n";
    }
  }

  if (cl.option_present('a')) {
    append_data_to_stored_data = 1;

    if (verbose) {
      cerr << "Will concatenate new data onto the end of existing data\n";
    }
  }

  if (cl.option_present('S')) {
    separator = cl.string_value('S');

    if ("vbar" == separator) {
      separator = "|";
    } else if ("space" == separator) {
      separator = " ";
    } else if ("tab" == separator) {
      separator = "\t";
    } else if ("semic" == separator) {
      separator = ";";
    } else if ("comma" == separator) {
      separator = ",";
    } else if ("colon" == separator) {
      separator = ":";
    } else if ("dollar" == separator) {
      separator = "$";
    } else if ("squote" == separator) {
      separator = "'";
    } else if ("dquote" == separator) {
      separator = '"';
    }

    append_data_to_stored_data = 1;

    if (verbose) {
      cerr << "Will concatenate new data onto end of existing data, separator is '"
           << separator << "'\n";
    }
  }

  if (append_data_to_stored_data) {
    store_flag = 0;
  }

  if (cl.option_present('x')) {
    const char* x = cl.option_value('x');

    if (!setup_regular_expression(data_rx_must_match, x)) {
      cerr << "Invalid data match regular expression '" << x << "'\n";
      exit(3);
    }

    if (verbose) {
      cerr << "Will only concatenate cases where data matches '"
           << data_rx_must_match->pattern() << "'\n";
    }
  }

  if (cl.option_present('X')) {
    const char* x = cl.option_value('X');

    if (!setup_regular_expression(data_rx_non_match, x)) {
      cerr << "Cannot initialise data non match regular expression '" << x << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will only concatenate cases where data does not match '"
           << data_rx_non_match->pattern() << "'\n";
    }
  }

  if (cl.option_present('w')) {
    perform_wc = 1;

    if (verbose) {
      cerr << "Will perform background wc before loading each database\n";
    }
  }

  if (cl.option_present('y')) {
    const_IWSubstring y = cl.string_value('y');

    if ("btree" == y) {
      ;
    } else if ("hash" == y) {
      dbtype = DB_HASH;
    } else if ("recno" == y) {
      dbtype = DB_RECNO;
    } else if ("queue" == y) {
      dbtype = DB_QUEUE;
    } else {
      cerr << "Unrecognised access method (-y) '" << y << "'\n";
    }
  }

  // the cache size idea never seemed to make much difference. Deprecated.
  int cache_size = 0;

  if (cl.option_present('c')) {
    if (!cl.value('c', cache_size) || cache_size < 1) {
      cerr << "The cache size option (-c) must be a whole positive number\n";
      usage(8);
    }

    if (verbose) {
      cerr << "cache size set to " << cache_size << '\n';
    }
  }

  if (cl.option_present('C')) {
    const_IWSubstring c = cl.string_value('C');

    if ("tokens" == c) {
      count_records_to_determine_tokens_per_record = 1;
    }
    if ("bytes" == c) {
      count_records_to_determine_bytes_per_record = 1;
    } else if (!c.numeric_value(tokens_per_record) || tokens_per_record < 1) {
      cerr << "The only load tokens per record (-C option) must be a whole positive "
              "number\n";
      usage(11);
    }
  }

  number_keys_to_match = cl.option_count('M');

  if (number_keys_to_match) {
    keys_to_match = new Key_to_Match[number_keys_to_match];

    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i)) {
      keys_to_match[i].set_key(m);

      i++;
    }

    if (verbose) {
      cerr << "For databases to be merged, the following keys must match\n";
      for (int i = 0; i < number_keys_to_match; i++) {
        cerr << " '" << keys_to_match[i].dkey() << "'\n";
      }
    }
  }

  if (cl.option_present('L')) {
    int i = 0;
    const_IWSubstring l;
    while (cl.value('L', l, i++)) {
      if (l.starts_with("maxtries=")) {
        l.remove_leading_chars(9);
        if (!l.numeric_value(max_open_attempts) || max_open_attempts < 1) {
          cerr << "The max open attempts value must be a whole positive number\n";
          usage(6);
        }

        if (verbose) {
          cerr << "Will try " << max_open_attempts
               << " times to open the database for writing\n";
        }

        continue;
      }

      //    Anything else must be the sleep time

      if (!l.numeric_value(sleep_before_next_attempt) || sleep_before_next_attempt < 1) {
        cerr << "The sleep between open attempts option (-L) must be a whole positive "
                "number\n";
        usage(5);
      }

      if (verbose) {
        cerr << "Will sleep " << sleep_before_next_attempt
             << " seconds between database open attempts\n";
      }
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting mechanism (-r)\n";
      usage(3);
    }
  }

  if (cl.option_present('p')) {
    special_processing_for_dicer_type_databases = cl.option_count('p');

    if (verbose) {
      cerr << "Special processing for databases with 'count ...' type data\n";
    }

    store_flag = 0;
    append_data_to_stored_data = 1;

    if (cl.option_present('b')) {
      store_items_not_matching_pattern = 1;

      if (verbose) {
        cerr << "Will store items even if they don't match the 'count ...' pattern\n";
      }
    }
  } else if (cl.option_present('q')) {
    special_processing_for_synthetic_precedent_databases = cl.option_count('q');

    if (verbose) {
      cerr << "Special processing for synthetic precedent databases\n";
    }

    store_flag = 0;
    append_data_to_stored_data = 1;
  } else if (cl.option_present('N')) {
    if (! cl.value('N', numeric_value_column) || numeric_value_column < 1) {
      cerr << "The numeric value column (-N) option must be a valid column number\n";
      return 1;
    }
    if (verbose) {
      cerr << "Numeric value expected in column " << numeric_value_column << '\n';
    }
    --numeric_value_column;
    store_flag = 0;
    append_data_to_stored_data = 1;
  }

  // We can either support iwbdb_cat <target> <db1> <db2> etc
  // or iwbdb_cat -d <target> <db1> <db2> etc..

  if (cl.option_present('d')) {
    ;
  } else if (cl.number_elements() < 2) {
    cerr << "Requires at least two databases\n";
    usage(1);
  }

  const char* target_fname;
  int istart;
  if (cl.option_present('d')) {
    target_fname = cl.option_value('d');
    istart = 0;
  } else {
    target_fname = cl[0];
    istart = 1;
  }

  Db target(NULL, DB_CXX_NO_EXCEPTIONS);

  if (cache_size > 0) {
    int rc = target.set_cachesize(0, cache_size, 1);
    if (0 != rc) {
      target.err(rc, "Cannot set cache size");
      exit(3);
    } else {
      cerr << "Database cache set to " << cache_size << " bytes\n";
    }
  }

  if (!open_target_database(target_fname, target)) {
    cerr << "Cannot open database '" << target_fname << "'\n";
    return 6;
  }

  if (number_keys_to_match) {
    (void)check_key_matches(target);
  }

  if (count_records_to_determine_tokens_per_record ||
      count_records_to_determine_bytes_per_record) {
    determine_tokens_per_record_in_target_database(target);
  }

  int rc = 0;
  for (int i = istart; i < cl.number_elements(); i++) {
    if (!iwbdb_cat(target, cl[i])) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << items_read << " from " << (cl.number_elements() - 1)
         << " databases\n";
    cerr << "Stored " << items_stored << " items\n";
    if (append_data_to_stored_data) {
      cerr << "Appended " << items_appended_to_existing_data
           << " items to existing stored data\n";
    }

    if (count_records_to_determine_tokens_per_record ||
        count_records_to_determine_bytes_per_record) {
      cerr << "Suppressed " << records_suppressed_by_size_or_token_mismatch
           << " records due to token count or size mismatch\n";
    }

    if (special_processing_for_dicer_type_databases && count_stats.n() > 1) {
      cerr << dicer_items_summed << " items incremented to existing data\n";
      cerr << "Final counts between " << count_stats.minval() << " and "
           << count_stats.maxval() << " ave " << static_cast<float>(count_stats.average())
           << '\n';
    }

    if (data_rx_must_match) {
      cerr << skipped_for_not_matching_data_rx
           << " records skipped for not matched data rx '"
           << data_rx_must_match->pattern() << "'\n";
    }

    if (data_rx_non_match) {
      cerr << skipped_for_matching_data_rx
           << " records skipped for matching non match rx '"
           << data_rx_non_match->pattern() << "'\n";
    }
  }

  target.close(0);

  return rc;
}

int
main(int argc, char** argv) {
  int rc = iwbdb_cat(argc, argv);

  return rc;
}
