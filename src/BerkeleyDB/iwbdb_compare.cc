/*
  Compare the contents of Berkeley databases.
  Note that we only check keys in the first database
*/

#include <ctype.h>

#include <iostream>
#include <memory>
#include <numeric>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "db_cxx.h"
#include "unpack_data.h"

using std::cerr;

static int verbose = 0;

static int strip_blanks = 0;

static int items_matching_only_with_blank_removal = 0;

static int keys_with_no_data_in_other_database = 0;

static IWString_and_File_Descriptor stream_for_keys_only_in_reference_database;

static IWString_and_File_Descriptor stream_for_new_and_changed_keys;

static int write_key_of_updates_first = 1;

static int display_no_match_message = 0;

static int perhaps_data_appended_in_different_order = 0;

static char concatenated_data_separator;

static int same_data_different_order = 0;

static Unpack_Binary_Data key_unpack_format, data_unpack_format;

static Report_Progress report_progress;

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
  cerr << "Compares the contents of Berkeley databases. Keys from the first database are\n";
  cerr << "   compared with the other databases\n";
  cerr << DB_VERSION_STRING << '\n';
  cerr << " -b              break on first difference\n";
  cerr << " -P <rx>         only process keys that match <rx>\n";
  cerr << " -w              if values don't match, remove leading and trailing blanks\n";
  cerr << "                 and try the comparison again\n";
  cerr << " -R <fname>      write keys only found in reference database to <fname>\n";
  cerr << " -D <fname>      write new and changed keys to <fname>\n";
  cerr << " -u <pattern>    unpack pattern for the data\n";
  cerr << " -U <pattern>    unpack pattern for the key\n";
  cerr << " -r <number>     report progess every <number> comparisons done\n";
  cerr << " -o <char>       maybe data is in different order, separated by <char> (normally :)\n";
  cerr << " -s              delta file is smiles, write data then key\n";
  cerr << " -q              suppress no match messages\n";
  cerr << " -v              verbose output\n";
// clang-format on

  exit(rc);
}

static int records_in_database = 0;

static int records_matching = 0;

static int break_on_first_difference = 0;

/*
  the number if different records in each database
*/

static int* different_records = NULL;

static int
contains_binary_characters(const Dbt& zdata) {
  const char* c = reinterpret_cast<const char*>(zdata.get_data());
  int n = zdata.get_size();

  for (int i = 0; i < n; i++) {
    if (isprint(c[i])) {
      continue;
    }

    if (isspace(c[i])) {  // tab, space or newline
      continue;
    }

    return 1;
  }

  return 0;
}

static void
write_datum(const char* header, const Dbt& d, const Unpack_Binary_Data& unpack_format,
            IWString_and_File_Descriptor& os) {
  os << header << ' ';

  if (NULL == d.get_data()) {
    os << "NULL\n";
    return;
  }

  os << "'";

  if (!contains_binary_characters(d)) {
    os.write(reinterpret_cast<const char*>(d.get_data()), d.get_size());
  } else if (unpack_format.active()) {
    unpack_format.write_unpacked_data(reinterpret_cast<const char*>(d.get_data()),
                                      d.get_size(), os);
  } else {
    os << d.get_size() << " binary";
  }

  os << "'\n";

  return;
}

static int
report_difference(const Dbt& dkey, const Dbt& reference_value, const Dbt& fromdb,
                  IWString_and_File_Descriptor& output) {
  write_datum("key", dkey, key_unpack_format, output);

  write_datum("reference", reference_value, data_unpack_format, output);

  write_datum("found    ", fromdb, data_unpack_format, output);

  output.write_if_buffer_holds_more_than(32768);

  return 0;
}

static int
just_different_order(const Dbt& reference_value, const Dbt& fromdb) {
  if (reference_value.get_size() != fromdb.get_size()) {
    return 0;
  }

  const_IWSubstring r(reinterpret_cast<const char*>(reference_value.get_data()),
                      reference_value.get_size());
  const_IWSubstring f(reinterpret_cast<const char*>(fromdb.get_data()),
                      fromdb.get_size());

  const int n = r.ccount(concatenated_data_separator);

  if (0 == n) {
    return 0;
  }

  if (n != f.ccount(concatenated_data_separator)) {
    return 0;
  }

  if (1 == n) {
    const_IWSubstring r1, r2, f1, f2;
    r.split(r1, concatenated_data_separator, r2);
    f.split(f1, concatenated_data_separator, f2);

    return r1 == f2 && r2 == f1;
  }

  IW_STL_Hash_Set rs, fs;

  int i = 0;
  const_IWSubstring s;
  while (r.nextword(s, i, concatenated_data_separator)) {
    rs.insert(s);
  }

  i = 0;
  while (f.nextword(s, i, concatenated_data_separator)) {
    fs.insert(s);
  }

  // for (auto i : rs)
  for (auto i = rs.begin(); i != rs.end(); ++i) {
    if (!rs.contains(*i)) {
      return 0;
    }
  }

  return 1;
}

static int
iwbdb_compare(const Dbt& dkey, const Dbt& reference_value, const Dbt& fromdb,
              IWString_and_File_Descriptor& output) {
  if (NULL == fromdb.get_data()) {  // reference_value is NOT null
    return report_difference(dkey, reference_value, fromdb, output);
  }

  if (0 == memcmp(reinterpret_cast<const char*>(reference_value.get_data()),
                  reinterpret_cast<const char*>(fromdb.get_data()),
                  fromdb.get_size())) {  // they are the same
    return 1;
  }

  if (strip_blanks) {
    const_IWSubstring r(reinterpret_cast<const char*>(reference_value.get_data()),
                        reference_value.get_size());
    const_IWSubstring f(reinterpret_cast<const char*>(fromdb.get_data()),
                        fromdb.get_size());

    r.strip_leading_blanks();
    r.strip_trailing_blanks();

    f.strip_leading_blanks();
    f.strip_trailing_blanks();

    if (r == f) {
      items_matching_only_with_blank_removal++;
      if (verbose > 1) {
        const_IWSubstring k(reinterpret_cast<const char*>(dkey.get_data()),
                            dkey.get_size());
        cerr << "Key " << k << " only matches with blanks removed\n";
      }
      return 1;
    }
  }

  if (perhaps_data_appended_in_different_order &&
      just_different_order(reference_value, fromdb)) {
    same_data_different_order++;
    return 1;
  }

  return report_difference(dkey, reference_value, fromdb, output);
}

static int
write_new_or_changed_value(const Dbt& zkey, const Dbt& zdata,
                           IWString_and_File_Descriptor& os) {
  if (write_key_of_updates_first) {
    os.strncat(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size());
    os << ' ';
    os.strncat(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());
  } else {
    os.strncat(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());
    os << ' ';
    os.strncat(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size());
  }

  os << '\n';
  os.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
process_no_match_for(const Dbt& zkey, const Dbt& zdata, std::ostream& os) {
  keys_with_no_data_in_other_database++;

  if (display_no_match_message) {
    os << "No match for key '";

    os.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size());

    os << "' in comparator databases\n";
  }

  if (stream_for_keys_only_in_reference_database.is_open()) {
    stream_for_keys_only_in_reference_database.strncat(
        reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size());
    stream_for_keys_only_in_reference_database << '\n';
    stream_for_keys_only_in_reference_database.write_if_buffer_holds_more_than(32768);
  }

  if (!stream_for_new_and_changed_keys.is_open()) {
    return 1;
  }

  return write_new_or_changed_value(zkey, zdata, stream_for_new_and_changed_keys);
}

static int
iwbdb_compare(int ndb, Db** db, IWString_and_File_Descriptor& output) {
  Db* first_database = db[0];

  Dbc* cursor = NULL;

  int rc = first_database->cursor(NULL, &cursor, 0);
  if (0 != rc) {
    first_database->err(rc, "cannot acquire cursor");
    return 0;
  }

  Dbt zkey, reference_value;

  while (0 == (rc = cursor->get(&zkey, &reference_value, DB_NEXT))) {
    records_in_database++;

    int found_key_in_another_database = 0;

    for (int i = 1; i < ndb; i++) {
      Dbt fromdb;

      if (0 != db[i]->get(NULL, &zkey, &fromdb, 0)) {
        continue;
      }

      found_key_in_another_database++;

      int c = iwbdb_compare(zkey, reference_value, fromdb, output);

      if (0 == c)  // data values are different
      {
        different_records[i]++;
        if (stream_for_new_and_changed_keys.is_open()) {
          write_new_or_changed_value(zkey, reference_value,
                                     stream_for_new_and_changed_keys);
        }
        if (break_on_first_difference) {
          return 1;
        }
      } else {
        ++records_matching;
      }
    }

    if (0 == found_key_in_another_database) {
      process_no_match_for(zkey, reference_value, cerr);
    }

    if (report_progress()) {
      cerr << "Processed " << records_in_database << " records";
      for (int i = 0; i < ndb; ++i) {
        cerr << ' ' << i << " diff " << different_records[i];
      }
      cerr << " no value in other db " << keys_with_no_data_in_other_database << '\n';
    }
  }

  return 1;
}

static int
iwbdb_compare(int argc, char** argv) {
  Command_Line cl(argc, argv, "vbwR:D:squ:U:r:o:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('b')) {
    break_on_first_difference = 1;

    if (verbose) {
      cerr << "Will stop processing once the first difference is found\n";
    }
  }

  if (cl.option_present('w')) {
    strip_blanks = 1;

    if (verbose) {
      cerr << "Will strip leading and trailing blanks when doing comparisons\n";
    }
  }

  int ndb = cl.number_elements();

  if (ndb < 2) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('u')) {
    const char* u = cl.option_value('u');
    if (!data_unpack_format.initialise(u)) {
      cerr << "Invalid data unpack format '" << u << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Data unpack format set to '" << u << "'\n";
    }
  }

  if (cl.option_present('U')) {
    const char* u = cl.option_value('U');

    if (!key_unpack_format.initialise(u)) {
      cerr << "Invalid key unpack format '" << u << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Key unpack format set to '" << u << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 2;
    }
  }

  if (cl.option_present('o')) {
    const_IWSubstring tmp = cl.string_value('o');
    if (1 != tmp.length()) {
      cerr << "The concatenated data separator option (-o) must be a single character, '"
           << tmp << "' invalid\n";
      usage(1);
    }

    concatenated_data_separator = tmp[0];

    if (verbose) {
      cerr << "Will check to see if data differs just by concatenation order\n";
    }

    perhaps_data_appended_in_different_order = 1;
  }

  Db** db = new Db*[ndb];

  for (int i = 0; i < ndb; i++) {
    db[i] = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

    int rc = db[i]->open(NULL, cl[i], NULL, DB_UNKNOWN, DB_RDONLY, 0);

    if (0 != rc) {
      cerr << "Cannot open database :";
      db[i]->err(rc, cl[i]);
      return i + 1;
    }
  }

  different_records = new_int(ndb);
  std::unique_ptr<int[]> free_different_records(different_records);

  if (cl.option_present('R') && cl.option_present('D')) {
    cerr << "Sorry, can't use the -R and -D options\n";
    usage(3);
  }

  if (cl.option_present('R')) {
    const char* r = cl.option_value('R');

    if (!stream_for_keys_only_in_reference_database.open(r)) {
      cerr << "Cannot open stream for keys only in reference db '" << r << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Keys only in reference database written to '" << r << "'\n";
    }
  } else if (cl.option_present('D')) {
    const char* d = cl.option_value('D');

    if (!stream_for_new_and_changed_keys.open(d)) {
      cerr << "Cannot open stream for new and updated keys '" << d << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "New keys and keys with different data written to '" << d << "'\n";
    }

    if (cl.option_present('s')) {
      write_key_of_updates_first = 0;
      if (verbose) {
        cerr << "In delta file will write data then key\n";
      }
    }
  }

  if (cl.option_present('q')) {
    display_no_match_message = 0;

    if (verbose) {
      cerr << "Will suppress display of no match messages\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = iwbdb_compare(ndb, db, output);

  output.flush();

  cerr << "Read " << records_in_database << " database records, " << 
          records_matching << " records matched across databases\n";
  for (int i = 1; i < ndb; i++) {
    cerr << "Database '" << cl[i] << "' " << different_records[i]
         << " records differing\n";

    db[i]->close(0);
    delete db[i];
  }

  delete[] db;

  if (verbose && items_matching_only_with_blank_removal) {
    cerr << items_matching_only_with_blank_removal
         << " items only matched with blank removal\n";
  }

  if (verbose && same_data_different_order) {
    cerr << same_data_different_order
         << " data values were the same, but in different concatenation order\n";
  }

  if (keys_with_no_data_in_other_database) {
    cerr << keys_with_no_data_in_other_database
         << " in reference database were unmatched in any comparator db's\n";
  }

  if (stream_for_keys_only_in_reference_database.is_open()) {
    stream_for_keys_only_in_reference_database.close();
  }

  if (stream_for_new_and_changed_keys.is_open()) {
    stream_for_new_and_changed_keys.close();
  }

  if (0 == rc) {  // some kind of failure
    return rc;
  }

  const auto d = std::accumulate(different_records, different_records + ndb, 0);
  if (0 == d && 0 == keys_with_no_data_in_other_database) {
    return 0;
  }

  return 1;
}

int
main(int argc, char** argv) {
  int rc = iwbdb_compare(argc, argv);

  return rc;
}
