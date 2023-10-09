/*
  Looks for items in a Berkeley DB database
*/

#include <stdlib.h>

#include <fstream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/report_progress.h"

#include "db_cxx.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static int remove_leading_zeros = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  cerr << DB_VERSION_STRING << '\n';
  // clang-format on

  // clang-format off
  cerr << "Checks for existence of keys in a Berkeley database\n";
  cerr << " -d <dbname>     database to scan\n";
  cerr << " -r <number>     report progress\n";
  cerr << " -c <col>        identifiers in column <col>\n";
  cerr << " -B <file>       file for identifiers not in database\n";
  cerr << " -z              remove leading 0's from identifiers before lookup\n";
  cerr << " -v              verbose output\n";
  // clang-format on

  exit(rc);
}

static int keys_looked_up = 0;

static int keys_found = 0;

static Report_Progress report_progress;

static int identifier_column = -1;

static IWString_and_File_Descriptor stream_for_not_found;

static int
iwbdb_exists3(const const_IWSubstring& id, const const_IWSubstring& buffer,
              Db& database) {
  Dbt dkey;

  dkey.set_data(const_cast<char*>(id.rawchars()));
  dkey.set_size(id.length());

  keys_looked_up++;

  int rc = database.exists(NULL, &dkey, 0);

  // cerr << "Looking up '" << id << "' from '" << buffer << "', rc " << rc << '\n';

  if (0 == rc) {
    keys_found++;
  } else if (DB_NOTFOUND == rc) {
    if (stream_for_not_found.is_open()) {
      stream_for_not_found << buffer << '\n';
      stream_for_not_found.write_if_buffer_holds_more_than(32768);
    }
  } else {
    cerr << "Unusual error from database.exists '";
    database.err(rc, "");
    cerr << '\n';
    return 0;
  }

  if (report_progress()) {
    cerr << "Looked up " << keys_looked_up << " found " << keys_found << " items\n";
  }

  return 1;
}

static int
iwbdb_exists2(const const_IWSubstring& id, const const_IWSubstring& buffer,
              Db& database) {
  // cerr << "remove_leading_zeros? " << remove_leading_zeros << " id '" << id << "'\n";
  if (remove_leading_zeros && id.starts_with('0')) {
    const_IWSubstring myid(id);
    myid.remove_leading_chars('0');
    return iwbdb_exists3(myid, buffer, database);
  } else {
    return iwbdb_exists3(id, buffer, database);
  }
}

static int
iwbdb_exists(const const_IWSubstring& buffer, Db& database) {
  if (identifier_column < 0) {
    return iwbdb_exists2(buffer, buffer, database);
  }

  const_IWSubstring token;

  if (!buffer.word(identifier_column, token)) {
    cerr << "Cannot extract identifier from '" << buffer << "'\n";
    return 0;
  }

  // cerr << "Got '" << token << "' from '" << buffer << "'\n";

  return iwbdb_exists2(token, buffer, database);
}

static int
iwbdb_exists(iwstring_data_source& input, Db& database) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!iwbdb_exists(buffer, database)) {
      return 0;
    }
  }

  return 1;
}

static int
iwbdb_exists(const char* fname, Db& database) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return iwbdb_exists(input, database);
}

static int
iwbdb_exists(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:r:c:B:z");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (!cl.option_present('d')) {
    cerr << "Must specify the database to use via the -d option\n";
    usage(31);
  }

  Db database(NULL, DB_CXX_NO_EXCEPTIONS);

  if (cl.option_present('d')) {
    const char* dbname = cl.option_value('d');

    int rc = database.open(NULL, dbname, NULL, DB_UNKNOWN, DB_RDONLY, 0);

    if (0 != rc) {
      cerr << "Cannot open database '" << dbname << "'\n";
      database.err(rc, "Didn't open");
      return 0;
    }

    if (verbose) {
      cerr << "Using database '" << dbname << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "The report progress option (-r) must be a whole positive number\n";
      usage(3);
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', identifier_column) || identifier_column < 1) {
      cerr << "INvalid identifier column (-c option)\n";
      usage(41);
    }

    if (verbose) {
      cerr << "Identifiers found in column " << identifier_column << '\n';
    }

    identifier_column--;  // we start with column 0
  }

  if (cl.option_present('z')) {
    remove_leading_zeros = 1;

    if (verbose) {
      cerr << "Will remove leading 0's from identifiers\n";
    }
  }

  if (cl.option_present('B')) {
    const char* b = cl.option_value('B');
    if (!stream_for_not_found.open(b)) {
      cerr << "Cannot open stream for not found identifiers '" << b << "'\n";
      return 15;
    }

    if (verbose) {
      cerr << "Identifiers not found written to '" << b << "'\n";
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!iwbdb_exists(cl[i], database)) {
      rc = i + 1;
      break;
    }
  }

  std::cout << "Looked up " << keys_looked_up << " found " << keys_found << " keys\n";

  database.close(0);

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = iwbdb_exists(argc, argv);

  return rc;
}
