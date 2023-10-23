/*
  Load a BerkeleyDB database
*/

#include <ctype.h>
#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"

#include "db_cxx.h"
#include "zlib.h"

#include "re2/re2.h"
#include <sys/stat.h>
#include <sys/time.h>

using std::cerr;

static int verbose = 0;

static const char* prog_name = NULL;

static std::unique_ptr<re2::RE2> key_rx;

static int column_for_identifier = 0;

static int overwrite_existing_records = 0;

static int append_to_stored_data = 0;

static IWString append_to_stored_data_separator = ':';

static int append_only_if_different_data = 0;

static resizable_array<int> columns_for_data;

static int records_read = 0;

static int records_stored = 0;

static int records_overwritten = 0;

static int records_appended = 0;

static int ignore_bad_input = 0;

static int invalid_records_ignored = 0;

static int report_progress = 0;

static IWString_and_File_Descriptor stream_for_items_loaded;

static struct timeval tprev;
static struct timeval tzero;

static int compress_with_zlib = 0;

static int lowercase_key = 0;
static int lowercase_data = 0;
static int uppercase_key = 0;
static int uppercase_data = 0;

static int need_to_change_case = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << DB_VERSION_STRING << '\n';

  cerr << "Loads records from a iwbdb database\n";
  cerr << prog_name << " -d dbname identifier_file\n";
  cerr << " -d <db>         database\n";
  cerr << " -a <char>       append new data to existing data, separated by <char>\n";
  cerr
      << " -A              only append new data if it is different from already stored\n";
  cerr << " -o              overwrite existing database entries\n";
  cerr << " -B <file>       file for identifiers already in the database\n";
  cerr << " -c <column>     column for identifiers (default 1)\n";
  cerr << " -C <column>     column(s) for the data to be stored\n";
  cerr << " -p <pattern>    pattern for identifiers\n";
  cerr << " -K <key>        key  value for command line items\n";
  cerr << " -V <key>        data value for command line items\n";
  cerr << " -r <number>     report progress every <number> stores\n";
  cerr << " -b              fault tolerant mode, just ignore seemingly invalid input\n";
  cerr << " -Z              compress data with zlib\n";
  cerr << " -y <type>       access type (btree, hash, recno, queue)\n";
  cerr << " -T <t>          UK UD lK lD Uppercase or Lowercase the Key and/or Data\n";
  cerr << " -v              verbose output\n";
  // clang-format on

  exit(rc);
}

// clang-format off
/*
  Make sure we handle the case of trying to store 'foo' when what is stored
  contains
  ...:foo:...
  foo:...
  ...:foo
*/
// clang-format on

static int
data_is_the_same(const Dbt& fromdb, const IWString& s) {
  if (0 ==
      s.strncmp(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size())) {
    return 1;
  }

  const_IWSubstring f(reinterpret_cast<const char*>(fromdb.get_data()),
                      fromdb.get_size());

  // There must be at least foo: in the found string

  if (f.length() + 1 < s.length()) {
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  while (f.nextword(token, i, append_to_stored_data_separator[0])) {
    if (s == token) {
      return 1;
    }
  }

  return 0;
}

/*
  Compute milliseconds since previous time structure. Put current
  values into prev
*/

static long
milliseconds_since(struct timeval& prev) {
  struct timeval tv;

  (void)gettimeofday(&tv, static_cast<struct timezone*>(NULL));

  long rc = tv.tv_sec - prev.tv_sec;
  assert(rc >= 0);

  rc = rc * 1000000;

  rc += tv.tv_usec - prev.tv_usec;

  prev.tv_sec = tv.tv_sec;
  prev.tv_usec = tv.tv_usec;

  return rc / 1000;
}

static int
do_zlib_compression(Dbt& zdata) {
  uLongf destlen = zdata.get_size() * 100 + 12;
  Bytef* dest = new Bytef[destlen];

  int rc = compress(dest, &destlen, reinterpret_cast<const Bytef*>(zdata.get_data()),
                    static_cast<uLongf>(zdata.get_size()));

  if (Z_OK == rc) {
    zdata.set_data(reinterpret_cast<char*>(dest));
    zdata.set_size(static_cast<int>(destlen));

    return 1;
  }

  delete[] dest;

  if (Z_MEM_ERROR == rc) {
    cerr << "Not enough memory for compress " << destlen << '\n';
  } else if (Z_BUF_ERROR == rc) {
    cerr << "Buffer too small for compression " << destlen << '\n';
  } else {
    cerr << "Unrecognised error from compress " << rc << '\n';
  }

  return 0;
}

static int
iwbdb_load(const const_IWSubstring& identifier, Db& database, Dbt& dbkey, Dbt& zdata) {
  if (compress_with_zlib) {
    do_zlib_compression(zdata);
  }

  int rc = database.put(NULL, &dbkey, &zdata, 0);

  if (0 == rc) {
    records_stored++;

    if (report_progress > 0 && 0 == records_stored % report_progress) {
      cerr << "Read " << records_read << " records, stored " << records_stored;
      long ms = milliseconds_since(tprev);
      if (ms > 0) {
        cerr << " rate " << (static_cast<long>(report_progress) * 1000 / ms)
             << " per sec";
      }
      cerr << '\n';
    }

    if (stream_for_items_loaded.is_open()) {
      stream_for_items_loaded.strncat((const char*)dbkey.get_data(), dbkey.get_size());
      stream_for_items_loaded << ' ';
      stream_for_items_loaded.strncat((const char*)zdata.get_data(), zdata.get_size());
      stream_for_items_loaded << '\n';
      stream_for_items_loaded.write_if_buffer_holds_more_than(32768);
    }

    if (compress_with_zlib) {
      delete[] static_cast<char*>(zdata.get_data());  // Cast to keep compiler quiet.
    }

    return 1;
  }

  cerr << "Cannot store: key '" << identifier << "' data " << zdata.get_size()
       << " bytes\n";
  database.err(rc, "");

  return 0;  // no fault tolerance for this, this is really bad
}

static int
iwbdb_load(const const_IWSubstring& identifier, const const_IWSubstring& to_store,
           Db& database, IWString_and_File_Descriptor& stream_for_already_stored) {
  Dbt dbkey((void*)identifier.rawchars(), identifier.length());

  // We only need to check if this key is present if we are writing the
  // already stored records, or cannot overwrite

  if (stream_for_already_stored.is_open() || (!overwrite_existing_records) ||
      append_to_stored_data || (overwrite_existing_records && verbose)) {
    Dbt fromdb;
    int rc = database.get(NULL, &dbkey, &fromdb, 0);

    if (0 == rc) {
      if (stream_for_already_stored.is_open()) {
        stream_for_already_stored << identifier << ' ' << to_store << '\n';
      }

      if (!append_to_stored_data) {
        ;
      } else if (append_only_if_different_data && data_is_the_same(fromdb, to_store)) {
        ;
      } else {
        IWString tmp;
        tmp.resize(fromdb.get_size() + append_to_stored_data_separator.length() +
                   to_store.length());

        tmp.strncat(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());
        tmp << append_to_stored_data_separator << to_store;

        records_appended++;

        Dbt zdata((void*)tmp.rawchars(), tmp.length());

        return iwbdb_load(identifier, database, dbkey, zdata);
      }

      if (!overwrite_existing_records) {
        return 1;
      }
    }

    records_overwritten++;
  }

  Dbt zdata((void*)to_store.rawchars(), to_store.length());

  return iwbdb_load(identifier, database, dbkey, zdata);
}

static int
determine_identifier(const const_IWSubstring& buffer, const_IWSubstring& identifier,
                     IWString& to_write) {
  if (!buffer.word(column_for_identifier, identifier)) {
    return 0;
  }

  to_write.resize(buffer.length());

  if (columns_for_data.empty()) {
    to_write = buffer;
    to_write.remove_word(column_for_identifier);

    return 1;
  }

  for (int i = 0; i < columns_for_data.number_elements(); i++) {
    if (i > 0) {
      to_write += ' ';
    }

    const_IWSubstring c;
    if (!buffer.word(columns_for_data[i], c)) {
      cerr << "Yipes, cannot extract column " << (columns_for_data[i] + 1) << " from '"
           << buffer << "'\n";
      return 0;
    }

    to_write += c;
  }

  return 1;
}

static int
load_after_changing_case(const const_IWSubstring& identifier,
                         const const_IWSubstring& to_store, Db& database,
                         IWString_and_File_Descriptor& stream_for_already_stored) {
  IWString newid(identifier);
  if (lowercase_key) {
    newid.to_lowercase();
  } else if (uppercase_key) {
    newid.to_uppercase();
  }

  IWString newdata(to_store);

  if (lowercase_data) {
    newdata.to_lowercase();
  } else if (uppercase_data) {
    newdata.to_uppercase();
  }

  return iwbdb_load(newid, newdata, database, stream_for_already_stored);
}

static int
iwbdb_load(iwstring_data_source& input, Db& database,
           IWString_and_File_Descriptor& stream_for_already_stored) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    records_read++;

    const_IWSubstring identifier;
    IWString to_store;

    if (!determine_identifier(buffer, identifier, to_store)) {
      cerr << "Cannot determine identifier from '" << buffer << "'\n";
      if (ignore_bad_input) {
        invalid_records_ignored++;
        continue;
      }

      return 0;
    }

    if (key_rx && !iwre2::RE2PartialMatch(identifier, *key_rx)) {
      cerr << "Invalid key '" << identifier << "'\n";

      if (ignore_bad_input) {
        continue;
      }

      return 0;
    }

    if (verbose > 2) {
      cerr << "Loading '" << identifier << "' value '" << to_store << "'\n";
    }

    if (need_to_change_case) {
      if (!load_after_changing_case(identifier, to_store, database,
                                    stream_for_already_stored)) {
        return 0;
      }
    } else if (!iwbdb_load(identifier, to_store, database, stream_for_already_stored)) {
      return 0;
    }
  }

  return 1;
}

// forward declaration

static int
iwbdb_load(const char* fname, Db& database,
           IWString_and_File_Descriptor& stream_for_already_stored);

static int
iwbdb_list_of_files(iwstring_data_source& input, Db& database,
                    IWString_and_File_Descriptor& stream_for_already_stored) {
  IWString buffer;

  while (input.next_record(buffer)) {
    if (!iwbdb_load(buffer.null_terminated_chars(), database,
                    stream_for_already_stored)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
iwbdb_load(const char* fname, Db& database,
           IWString_and_File_Descriptor& stream_for_already_stored) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  if (0 == strncmp(fname, "F:", 2)) {
    return iwbdb_list_of_files(input, database, stream_for_already_stored);
  }

  return iwbdb_load(input, database, stream_for_already_stored);
}

int
iwbdb_load(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:p:c:C:oa:AK:V:br:B:h:Zy:L:T:N:O:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('r')) {
    if (!cl.value('r', report_progress) || report_progress < 1) {
      cerr << "The report progress option (-r) must be a whole positive number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will report progress every " << report_progress << " items stored\n";
    }

    milliseconds_since(tzero);
    milliseconds_since(tprev);
  }

  if (cl.option_present('Z')) {
    compress_with_zlib = 1;

    if (verbose) {
      cerr << "Will compress data with zlib\n";
    }
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify the database to use via the -d option\n";
    usage(31);
  }

  DbEnv dbenv(DB_CXX_NO_EXCEPTIONS);

  DbEnv* dbenv_ptr;

  if (cl.option_present('h')) {
    IWString zenv = cl.string_value('h');

    int rc = dbenv.open(zenv.null_terminated_chars(),
                        DB_CREATE | DB_INIT_MPOOL | DB_INIT_CDB | DB_USE_ENVIRON, 0644);
    if (0 != rc) {
      cerr << "Yipes, cannot open the environment, '" << zenv << "', rc = " << rc << " "
           << dbenv.strerror(rc) << '\n';
      dbenv.close(0);
      return 4;
    }

    if (verbose) {
      cerr << "Opened environment '" << zenv << "'\n";
    }

    dbenv_ptr = &dbenv;
  } else {
    dbenv_ptr = NULL;
  }

  DBTYPE dbtype = DB_BTREE;

  if (cl.option_present('y')) {
    const_IWSubstring t = cl.string_value('y');

    if ("btree" == t) {
      ;
    } else if ("hash" == t) {
      dbtype = DB_HASH;
    } else if ("recno" == t) {
      dbtype = DB_RECNO;
    } else if ("queue" == t) {
      dbtype = DB_QUEUE;
    } else {
      cerr << "Unrecognised access method (-y) '" << t << "'\n";
    }
  }

  Db database(dbenv_ptr, DB_CXX_NO_EXCEPTIONS);

  if (cl.option_present('O')) {
    const_IWSubstring o;
    for (int i = 0; cl.value('O', o, i); i++) {
      if (o.starts_with("h_ffactor=")) {
        o.remove_leading_chars(10);
        int h;
        if (!o.numeric_value(h) || h < 1) {
          cerr << "Invalid hash fill factor '" << o << "'\n";
          return 2;
        }

        int rc = database.set_h_ffactor(h);
        if (0 == rc) {
          cerr << "Hash fill factor set to " << h << '\n';
        } else {
          cerr << "Unable to set hash fill factor to '" << h << "'\n";
          database.err(rc, "");
          return 2;
        }
      } else if (o.starts_with("h_nelem=")) {
        o.remove_leading_chars(8);
        int h;
        if (!o.numeric_value(h) || h < 1) {
          cerr << "Invalid hash nelem directive '" << h << '\n';
          return 2;
        }

        int rc = database.set_h_nelem(h);
        if (0 == rc) {
          cerr << "Set hash nelem " << h << '\n';
        } else {
          cerr << "Cannot set h_nelem " << h << ", ";
          database.err(rc, "");
          return 2;
        }
      } else if (o.starts_with("cachesize=")) {
        o.remove_leading_chars(10);
        int h;
        if (!o.numeric_value(h) || h < 1) {
          cerr << "Invalid cachesize directive '" << h << '\n';
          return 2;
        }

        int rc = database.set_cachesize(0, h, 1);
        if (0 == rc) {
          cerr << "Set cache size " << h << '\n';
        } else {
          cerr << "Cannot set cachesize " << h << ", ";
          database.err(rc, "");
          return 2;
        }
      } else {
        cerr << "Unrecognised db option specification '" << o << "'\n";
        return 3;
      }
    }
  }

  if (cl.option_present('d')) {
    IWString dbname = cl.string_value('d');

    u_int32_t flags;

    if (dash_s(dbname))  // if the file already exists, don't specify a type
    {
      dbtype = DB_UNKNOWN;
      flags = 0;
    } else {
      flags = DB_CREATE;
    }

    // do any cache size setting here

    int mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

    int rc =
        database.open(NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);

    if (0 != rc) {
      cerr << "Cannot open database '" << dbname << "' '";
      database.err(rc, "");
      return 0;
    }

    if (verbose) {
      cerr << "Using database '" << dbname << "'\n";
    }
  }

  if (cl.option_present('b')) {
    ignore_bad_input = 1;

    if (verbose) {
      cerr << "Will ignore records for which no identifier can be determined\n";
    }
  }

  if (cl.option_present('p')) {
    const_IWSubstring p;
    cl.value('p', p);

    if (!iwre2::RE2Reset(key_rx, p)) {
      cerr << "Could not compile key regular expression '" << p << "'\n";
      return 13;
    }

    if (verbose) {
      cerr << "Will only display records for which key matches '" << key_rx->pattern()
           << "'\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', column_for_identifier) || column_for_identifier < 1) {
      cerr << "INvalid identifier column (-c option)\n";
      usage(41);
    }

    if (verbose) {
      cerr << "Identifiers found in column " << column_for_identifier << '\n';
    }

    column_for_identifier--;  // we start with column 0
  }

  if (cl.option_present('C')) {
    int i = 0;
    const_IWSubstring C;
    while (cl.value('C', C, i++)) {
      int c;
      if (!C.numeric_value(c) || c < 1) {
        cerr << "INvalid column '" << C << "'\n";
        return 3;
      }

      columns_for_data.add(c - 1);
    }
  }

  if (cl.option_present('o') && (cl.option_present('a') || cl.option_present('A'))) {
    cerr << "The overwrite (-o) and append (-a) options are mutually incompatible\n";
    usage(3);
  }

  if (cl.option_present('o')) {
    overwrite_existing_records = 1;
    if (verbose) {
      cerr << "Existing database contents can be overwritten\n";
    }
  }

  if (cl.option_present('a')) {
    append_to_stored_data = 1;

    const_IWSubstring a = cl.string_value('a');

    if (a.length() > 1) {
      cerr << "Sorry, the inter-item database separator is restricted to one character\n";
      return 6;
    }

    append_to_stored_data_separator = a[0];

    if (verbose) {
      cerr << "Will append new data to existing stored data, separator '"
           << append_to_stored_data_separator << "'\n";
    }
  }

  if (cl.option_present('A')) {
    append_to_stored_data = 1;
    append_only_if_different_data = 1;

    if (verbose) {
      cerr << "Will only append data to stored items if different\n";
    }
  }

  if (cl.option_present('T')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('T', t, i++)) {
      if ("UK" == t) {
        uppercase_key = 1;
      } else if ("lK" == t) {
        lowercase_key = 1;
      } else if ("UD" == t) {
        uppercase_data = 1;
      } else if ("lD" == t) {
        lowercase_data = 1;
      } else {
        cerr << "Unrecognised -T qualifier '" << t << "'\n";
        return 3;
      }
    }

    need_to_change_case = 1;
  }

  // Items to load must come from either a file on the command line,
  // or pairs in -K -V pairs

  int nk = cl.option_count('K');

  if (nk != cl.option_count('V')) {
    cerr << "There must be the same number of key (-K) and value (-V) options\n";
    usage(7);
  }

  if (nk) {  // good, items on the command line
    ;
  } else if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  IWString_and_File_Descriptor stream_for_already_stored;

  if (cl.option_present('B')) {
    IWString b = cl.string_value('B');
    if (!stream_for_already_stored.open(b.null_terminated_chars())) {
      cerr << "Cannot open stream for already stored identifiers '" << b << "'\n";
      return 15;
    }

    if (verbose) {
      cerr << "Keys already stored written to '" << b << "'\n";
    }
  }

  if (cl.option_present('L')) {
    IWString l = cl.string_value('L');
    if (!stream_for_items_loaded.open(l.null_terminated_chars())) {
      cerr << "Cannot open stream for loaded data '" << l << "'\n";
      return 15;
    }

    if (verbose) {
      cerr << "Data stored into the database written to '" << l << "'\n";
    }
  }

  int rc = 0;

  // first do the items on the command line

  for (int i = 0; i < nk; i++) {
    const_IWSubstring k = cl.string_value('K', i);
    const_IWSubstring v = cl.string_value('V', i);

    if (!iwbdb_load(k, v, database, stream_for_already_stored)) {
      rc = i + 1;
      break;
    }
  }

  // If one of the command line items failed, don't do the file

  for (int i = 0; i < cl.number_elements() && 0 == rc; i++) {
    if (!iwbdb_load(cl[i], database, stream_for_already_stored)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << records_read << " records, stored " << records_stored << '\n';
    if (ignore_bad_input && invalid_records_ignored) {
      cerr << "Ignored " << invalid_records_ignored << " invalid input records\n";
    }
    long ms = milliseconds_since(tzero);

    if (ms > 0 && records_stored > 100) {
      cerr << "Rate "
           << (static_cast<double>(records_stored) * 1000.0 / static_cast<double>(ms))
           << " records per second\n";
    }

    if (overwrite_existing_records) {
      cerr << records_overwritten << " existing records overwritten\n";
    }
    if (append_to_stored_data) {
      cerr << records_appended << " items appended to existing stored data\n";
    }
  }

  database.close(0);

  if (NULL != dbenv_ptr) {
    dbenv.close(0);
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = iwbdb_load(argc, argv);

  return rc;
}
