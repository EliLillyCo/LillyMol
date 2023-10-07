// Delete specific keys in a BerkeleyDB database

#include <iostream>
#include <memory>

#if (__GNUC__ >= 3)
#define HAVE_CXX_STDHEADERS 1
#endif

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"

#include "db_cxx.h"

#include "re2/re2.h"

using std::cerr;

static int verbose = 0;

static int items_deleted = 0;

static int ignore_keys_not_in_database = 0;

static int report_keys_not_in_database = 1;

static Db database(NULL, DB_CXX_NO_EXCEPTIONS);

static IWString_and_File_Descriptor stream_for_deleted_records;

static int identifier_column = -1;

static int remove_leading_zeros = 0;

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

  cerr << "Deletes entries from a Berkeley DB database\n";
  cerr << "usage: <options> file1 file2 ...\n";
  cerr << " -d <dbname>     specify database from which to remove records\n";
  cerr << " -I <identifier> specify individual identifier(s) to remove\n";
  cerr << " -p <pattern>    delete all records with KEY  which match <pattern>\n";
  cerr << " -P <pattern>    delete all records with DATA which matches <pattern>\n";
  cerr << " -g              ignore keys not in the database\n";
  cerr << " -q              don't write messages about keys not in database\n";
  cerr << " -X <fname>      write deleted records to <fname>\n";
  cerr << " -c <col>        identifiers to delete are in column <col>\n";
  cerr << " -z              remove leading zero's from identifiers\n";
  cerr << " -v              verbose output\n";
// clang-format on

  exit(rc);
}

static int
write_existing_record(const Dbt& zkey, const Dbt& zdata) {
  stream_for_deleted_records.write(reinterpret_cast<const char*>(zkey.get_data()),
                                   zkey.get_size());
  stream_for_deleted_records << ' ';
  stream_for_deleted_records.write(reinterpret_cast<const char*>(zdata.get_data()),
                                   zdata.get_size());
  stream_for_deleted_records << '\n';

  stream_for_deleted_records.write_if_buffer_holds_more_than(32768);

  return stream_for_deleted_records.good();
}

static int
write_existing_record(Dbt& zkey) {
  Dbt zdata;

  int rc = database.get(NULL, &zkey, &zdata, 0);

  if (0 != rc)  // huh, couldn't fetch the data
  {
    if (ignore_keys_not_in_database) {
      return 1;
    }

    cerr << "Cannot fetch existing data for '";
    cerr.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size()) << "'\n";
    database.err(rc, "");
    return 0;
  }

  return write_existing_record(zkey, zdata);
}

static int
iwbdb_delete_pattern(std::unique_ptr<re2::RE2>& key_rx,
                     std::unique_ptr<re2::RE2>& data_rx) {
  Dbc* cursor = NULL;
  if (int rc = database.cursor(NULL, &cursor, 0); rc != 0) {
    database.err(rc, "cannot acquire cursor");
    return 0;
  }

  int records_in_database = 0;
  int items_deleted_this_database = 0;

  Dbt zdata, zkey;

  int rc;
  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    records_in_database++;

    int matched = 0;
    int written = 0;  // to the deleted records stream

    if (key_rx) {
      re2::StringPiece tmp(reinterpret_cast<const char*>(zkey.get_data()),
                           zkey.get_size());
      matched = RE2::PartialMatch(tmp, *key_rx);
      if (matched && verbose > 1) {
        cerr << "Key '";
        cerr.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size())
            << "' to be deleted\n";
      }
    }

    if (0 == matched && data_rx) {
      re2::StringPiece tmp(reinterpret_cast<const char*>(zdata.get_data()),
                           zdata.get_size());
      matched = RE2::PartialMatch(tmp, *data_rx);

      if (matched && stream_for_deleted_records.is_open()) {
        (void)write_existing_record(zkey, zdata);
        written = 1;
      }

      if (matched && verbose) {
        cerr << "Key '";
        cerr.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size())
            << "' data '";
        cerr.write(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size())
            << "' to be deleted\n";
      }
    }

    if (matched) {
      int rc = cursor->del(0);
      if (0 != rc) {
        database.err(rc, "Cannot delete ");
        cerr.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size());
        cerr << '\n';
        return 0;
      }

      items_deleted_this_database++;
      if (!written && stream_for_deleted_records.is_open()) {
        if (!write_existing_record(zkey)) {
          return 0;
        }
      }
    }
  }

  cursor->close();

  if (verbose) {
    cerr << "Scanned " << records_in_database << " records, "
         << items_deleted_this_database << " match ";
    if (key_rx) {
      cerr << " key pattern '" << key_rx->pattern() << "' ";
    }
    if (data_rx) {
      cerr << " data pattern '" << data_rx->pattern() << "'";
    }
    cerr << '\n';
  }

  items_deleted += items_deleted_this_database;

  return 1;
}

static int
iwbdb_delete_pattern(const const_IWSubstring& key_pattern,
                     const const_IWSubstring& data_pattern) {
  std::unique_ptr<re2::RE2> key_rx;

  if (key_pattern.length() && !iwre2::RE2Reset(key_rx, key_pattern)) {
    cerr << "Cannot parse key regular expression pattern '" << key_pattern << "'\n";
    return 0;
  }

  std::unique_ptr<re2::RE2> data_rx;

  if (data_pattern.length() && !iwre2::RE2Reset(data_rx, data_pattern)) {
    cerr << "Cannot parse data regular expression pattern '" << data_pattern << "'\n";
    return 0;
  }

  return iwbdb_delete_pattern(key_rx, data_rx);
}

static int
iwbdb_delete(const const_IWSubstring& zkey, Dbt& dkey) {
  // If we are supposed to be writing the deleted records, we need to fetch that now

  if (stream_for_deleted_records.is_open()) {
    if (!write_existing_record(dkey)) {
      return 0;
    }
  }

  int rc = database.del(NULL, &dkey, 0);

  if (0 == rc) {
    if (verbose > 1) {
      cerr << "Key '" << zkey << "' deleted\n";
    }

    items_deleted++;

    return 1;
  }

  // Why couldn't we delete this record? Was it because it isn't in the database?

  Dbt fromdb;

  rc = database.get(NULL, &dkey, &fromdb, 0);

  if (0 == rc) {  // great, found it in the database
    ;
  } else if (!report_keys_not_in_database &&
             ignore_keys_not_in_database) {  // not there and we don't say anything and
                                             // ignore it
    return 1;
  } else {
    cerr << "Key '" << zkey << " not in database";
    if (ignore_keys_not_in_database) {
      cerr << ", ignored\n";
      return 1;
    }
    cerr << '\n';
    return 0;
  }

  database.err(rc, "Delete failed");

  return 0;
}

static int
iwbdb_delete(const const_IWSubstring& zkey) {
  Dbt dkey;

  if (identifier_column < 0) {
    const_IWSubstring mykey(zkey);
    if (remove_leading_zeros) {
      mykey.remove_leading_chars('0');
    }

    dkey.set_data(const_cast<char*>(mykey.rawchars()));  // loss of const OK
    dkey.set_size(mykey.length());

    return iwbdb_delete(mykey, dkey);
  }

  const_IWSubstring id;
  if (!zkey.word(identifier_column, id)) {
    cerr << "Cannot extract column " << identifier_column << " from '" << zkey << "'\n";
    return 0;
  }

  if (remove_leading_zeros) {
    id.remove_leading_chars('0');
  }

  dkey.set_data(const_cast<char*>(id.rawchars()));
  dkey.set_size(id.length());

  return iwbdb_delete(zkey, dkey);
}

static int
iwbdb_delete(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (iwbdb_delete(buffer)) {
      ;
    } else {
      cerr << "Yipes, cannot delete key on line '" << input.lines_read() << '\n';
      return 0;
    }
  }

  return 1;
}

static int
iwbdb_delete(const char* fname) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open file of keys to delete '" << fname << "'\n";
    return 0;
  }

  return iwbdb_delete(input);
}

static int
iwbdb_delete(int argc, char** argv) {
  Command_Line cl(argc, argv, "d:vI:gqp:P:X:c:z");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('d')) {
    cerr << "Must specify database via the -d option\n";
    usage(3);
  }

  IWString dbname;
  cl.value('d', dbname);

  int rc = database.open(NULL, dbname.null_terminated_chars(), NULL, DB_UNKNOWN, 0, 0);
  if (0 != rc) {
    cerr << "Sorry, cannot open database '" << dbname << "'\n";
    database.err(rc, "");
    cerr << '\n';
    return 5;
  }

  if (cl.option_present('g')) {
    ignore_keys_not_in_database = 1;
    if (verbose) {
      cerr << "Will ignorekeys not in the database\n";
    }
  }

  if (cl.option_present('q')) {
    report_keys_not_in_database = 0;

    if (verbose) {
      cerr << "No messages about keys not in database\n";
    }
  }

  if (cl.option_present('X')) {
    IWString x;
    cl.value('X', x);

    if (!stream_for_deleted_records.open(x.null_terminated_chars())) {
      cerr << "Yipes, cannot open stream for deleted records '" << x << "'\n";
      return 28;
    }
  }

  if (cl.option_present('I')) {
    int i = 0;
    const_IWSubstring zkey;
    while (cl.value('I', zkey, i++)) {
      if (!iwbdb_delete(zkey)) {
        return i;
      }
    }

    if (cl.empty()) {
      return 0;
    }
  }

  if (cl.option_present('p') || cl.option_present('P')) {
    const_IWSubstring key_pattern, data_pattern;

    cl.value('p', key_pattern);
    cl.value('P', data_pattern);

    return !iwbdb_delete_pattern(key_pattern, data_pattern);
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', identifier_column) || identifier_column < 1) {
      cerr << "Invalid column for identifiers (-c option)\n";
      usage(7);
    }

    if (verbose) {
      cerr << "Identifiers assumed to be in column " << identifier_column << '\n';
    }

    identifier_column--;
  }

  if (cl.option_present('z')) {
    remove_leading_zeros = 1;

    if (verbose) {
      cerr << "Will remove leading zero's from identifiers\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!iwbdb_delete(cl[i])) {
      return i + 1;
    }
  }

  if (verbose) {
    cerr << "Deleted " << items_deleted << " items\n";
  }

  database.close(0);

  return 0;
}

int
main(int argc, char** argv) {
  int rc = iwbdb_delete(argc, argv);

  return rc;
}
