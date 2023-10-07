// Load a BerkeleyDB database from the contents of a TDT file.

#include <stdlib.h>
#include <zlib.h>

#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "db_cxx.h"

#include <sys/stat.h>

using std::cerr;
using std::endl;

static int verbose = 0;

static Report_Progress report_progress;

static int tdts_read = 0;

static int items_stored = 0;

static IWString identifier_tag("PCN<");

static std::unique_ptr<re2::RE2> rx;

static int store_flag = DB_NOOVERWRITE;

static int ignore_tdts_with_no_identifier = 0;

static int gzip_data = 0;

static int fault_tolerant = 0;

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
  cerr << DB_VERSION_STRING << endl;

  cerr << "Builds a Berkeley database from a tdt\n";
  cerr << " -d <dbname>      name of database to load\n";
  cerr << " -o               overwrite existing database entries\n";
  cerr << " -I <tag>         identifier tag - becomes the key\n";
  cerr << " -X <rx>          regular expression for identifiers\n";
  cerr << " -r <number>      report progress every <number> tdts processed\n";
  cerr << " -k               ignore TDT's with no identifier - by default we fail\n";
  cerr << " -f               fault tolerant mode - ignore most errors\n";
  cerr << " -z               compress data with zlib\n";
  cerr << " -y <type>        access type (btree, hash, recno, queue)\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
convert_to_gzipd_form(const const_IWSubstring& descriptors, Dbt& to_store) {
  uLongf destlen = descriptors.length() * 100 + 12;
  Bytef* dest = new Bytef[destlen];

  int rc =
      compress(dest, &destlen, reinterpret_cast<const Bytef*>(descriptors.rawchars()),
               static_cast<uLongf>(descriptors.length()));

  if (Z_OK == rc) {
    to_store.set_data(reinterpret_cast<char*>(dest));
    to_store.set_size(static_cast<int>(destlen));

    return 1;
  }

  delete[] dest;

  if (Z_MEM_ERROR == rc) {
    cerr << "Not enough memory for compress " << destlen << endl;
  } else if (Z_BUF_ERROR == rc) {
    cerr << "Buffer too small for compression " << destlen << endl;
  } else {
    cerr << "Unrecognised error from compress " << rc << endl;
  }

  return 0;
}

static int
iwbdb_from_tdt(const IWString& id, const const_IWSubstring& tostore, Db& database) {
  Dbt dkey(const_cast<char*>(id.rawchars()), id.length());

  Dbt zdata;

  if (gzip_data) {
    convert_to_gzipd_form(tostore, zdata);
  } else {
    zdata.set_data(const_cast<char*>(tostore.rawchars()));
    zdata.set_size(tostore.length());
  }

  int rc = database.put(NULL, &dkey, &zdata, store_flag);

  if (gzip_data) {
    delete reinterpret_cast<char*>(zdata.get_data());
  }

  if (0 == rc) {
    items_stored++;
    return 1;
  }

  if (DB_KEYEXIST == rc && DB_NOOVERWRITE == store_flag) {
    cerr << "Cannot overwrite existing data for '" << id << "'\n";
    return 1;
  }

  cerr << "Fatal error storing '" << id << "' ";
  database.err(rc, "");
  return 0;
}

static int
iwbdb_from_tdt(const IW_TDT& tdt, const IWString& id, Db& database) {
  const_IWSubstring zdata = tdt.rawdata();

  if (zdata.ends_with('\n')) {
    zdata.chop();
  }

  return iwbdb_from_tdt(id, zdata, database);
}

static int
iwbdb_from_tdt(const IW_TDT& tdt, Db& database) {
  IWString id;
  if (!tdt.dataitem_value(identifier_tag, id)) {
    cerr << "Cannot extract '" << identifier_tag << "' from tdt\n";
    return ignore_tdts_with_no_identifier;
  }

  if (rx && !iwre2::RE2FullMatch(id, *rx)) {
    cerr << "Invalid identifier '" << id << "' does not match '" << rx->pattern()
         << "'\n";
    return fault_tolerant;
  }

  return iwbdb_from_tdt(tdt, id, database);
}

static int
iwbdb_from_tdt(iwstring_data_source& input, Db& database) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    tdts_read++;

    if (report_progress()) {
      cerr << "Processing " << tdts_read << " tdt's\n";
    }

    if (!iwbdb_from_tdt(tdt, database)) {
      cerr << "Fatal error processing tdt " << tdts_read << endl;
      cerr << tdt;

      return 0;
    }
  }

  return 1;
}

static int
iwbdb_from_tdt(const char* fname, Db& database) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return iwbdb_from_tdt(input, database);
}

static int
iwbdb_from_tdt(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:oI:r:X:kzfy:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('d')) {
    cerr << "Must specify database via the -d option\n";
    usage(3);
  }

  if (cl.option_present('I')) {
    identifier_tag = cl.string_value('I');

    if (verbose) {
      cerr << "Identifier tag '" << identifier_tag << "'\n";
    }

    if (!identifier_tag.ends_with('<')) {
      identifier_tag += '<';
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "The report value (-r option) must be a whole number > 0\n";
      usage(5);
    }
  }

  if (cl.option_present('z')) {
    gzip_data = 1;

    if (verbose) {
      cerr << "Will gzip stored data\n";
    }
  }

  if (cl.option_present('f')) {
    fault_tolerant = 1;

    if (verbose) {
      cerr << "Will operate in a fault tolerant mode\n";
    }
  }

  if (cl.option_present('X')) {
    IWString x = cl.string_value('X');

    if (!iwre2::RE2Reset(rx, x)) {
      cerr << "Invalid regular expression pattern '" << x << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Regular expression for identifiers '" << rx->pattern() << "'\n";
    }
  }

  if (cl.option_present('o')) {
    store_flag = 0;

    if (verbose) {
      cerr << "Will overwrite existing entries\n";
    }
  }

  if (cl.option_present('k')) {
    ignore_tdts_with_no_identifier = 1;

    if (verbose) {
      cerr << "Will ignore TDT's from which no identifier can be extracted\n";
    }
  }

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    usage(4);
  }

  Db database(NULL, DB_CXX_NO_EXCEPTIONS);

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

  if (cl.option_present('d')) {
    IWString dbname = cl.string_value('d');

    int mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

    u_int32_t flags;

    if (dash_s(dbname.null_terminated_chars())) {
      dbtype = DB_UNKNOWN;
      flags = 0;
    } else {
      flags = DB_CREATE;
    }

    int rc =
        database.open(NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);

    if (0 != rc) {
      cerr << "Cannot open database '" << dbname << "' '";
      database.err(rc, "");
      return 0;
    }

    if (verbose) {
      cerr << "opened database '" << dbname << "' for writing\n";
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!iwbdb_from_tdt(cl[i], database)) {
      rc = i + 1;
      break;
    }
  }

  database.close(0);

  if (verbose) {
    cerr << "Processed " << tdts_read << " TDT's, " << items_stored << " items stored\n";
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = iwbdb_from_tdt(argc, argv);

  return rc;
}
