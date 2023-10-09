/*
  Lists a Berkeley database
*/

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/report_progress.h"

#include "db_cxx.h"
#include "unpack_data.h"
#include "zlib.h"

#include "re2/re2.h"

using std::cerr;

static int verbose = 0;

static int output_needed = 1;

static int print_binary_data = 0;

static IWString print_binary_data_as_daylight_fingerprint_tag;

static int byte_swap_binary_data = 0;

static int write_database_contents_first = 0;

static int write_database_key = 1;

static int write_the_data = 1;

static int records_to_skip = 0;

static int records_to_process = 0;

static std::unique_ptr<re2::RE2> key_rx, data_rx;

static Unpack_Binary_Data key_unpack_format, data_unpack_format;

static int need_to_uncompress_with_zlib = 0;

static int flush_when = 32768;

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
  cerr << DB_VERSION_STRING << '\n';

  cerr << "Lists the contents of a Berkeley DB database\n";
  cerr << " -b              print binary data\n";
  cerr << " -p <pattern>    only print data for keys which match <pattern>\n";
  cerr << " -P <pattern>    only print data for keys with data which matches <pattern>\n";
  cerr << " -u <pattern>    unpack pattern for the data\n";
  cerr << " -U <pattern>    unpack pattern for the key\n";
  cerr << " -w              write database contents before key values (smiles for example)\n";
  cerr << " -n              no output - useful for counting items in a database\n";
  cerr << " -f <number>     skip the first <number> records\n";
  cerr << " -o <number>     only write <number> records\n";
  cerr << " -k              suppress writing the key\n";
  cerr << " -D <tag>        write data in Daylight ASCII bit represenation form\n";
  cerr << " -a              byte swap before forming Daylight ASCII representation\n";
  cerr << " -Z              stored data is compressed by gzip, uncompress\n";
  cerr << " -h <size>       flush output buffers when <size> bytes accumulated\n";
  cerr << " -r <records>    report progress every <records> records retrieved\n";
  cerr << " -v              verbose output\n";
  // clang-format on

  exit(rc);
}

static int
uncompress_and_write(const Dbt& fromdb, IWString_and_File_Descriptor& output) {
  int buffer_size = fromdb.get_size() * 120;

  for (int i = 0; i < 10; i++)  // try 10 times
  {
    char* expanded = new char[buffer_size];
    std::unique_ptr<char[]> free_expanded(expanded);
    uLongf expanded_size = buffer_size;

    int rc = uncompress(reinterpret_cast<Bytef*>(expanded), &expanded_size,
                        reinterpret_cast<const Bytef*>(fromdb.get_data()),
                        static_cast<uLong>(fromdb.get_size()));
    if (Z_OK == rc) {
      output.write(expanded, expanded_size);
      return 1;
    }

    if (Z_MEM_ERROR == rc) {
      cerr << "Not enough memory for expansion, tried " << i << " times\n";
      return 0;
    } else if (Z_BUF_ERROR == rc) {
      cerr << "Buffer not large enough, input " << fromdb.get_size() << " tried "
           << buffer_size << '\n';
    } else if (Z_DATA_ERROR == rc) {
      cerr << "Data corrupted after " << i << " times\n";
      return 0;
    } else {
      cerr << "Unrecognised return code from uncompress " << rc << '\n';
    }

    buffer_size = buffer_size + buffer_size;
  }

  cerr << "decode_gzipd:cannot expand " << fromdb.get_size() << " bytes\n";
  return 0;
}

static int
contains_binary_characters(const Dbt& zdata) {
  const char* c = reinterpret_cast<const char*>(zdata.get_data());
  for (unsigned int i = 0; i < zdata.get_size(); i++) {
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



/*
  This function is made complex by the possibility that p may not be
  aligned properly
*/

template <typename T>
int
write_as(const char* p, IWString_and_File_Descriptor& output) {
  T tmp;

  char* ptmp = reinterpret_cast<char*>(&tmp);

  memcpy(ptmp, p, sizeof(T));

  output << tmp;

  return output.good();
}

template int
write_as<unsigned int>(const char* p, IWString_and_File_Descriptor& output);
template int
write_as<int>(const char* p, IWString_and_File_Descriptor& output);
template int
write_as<unsigned char>(const char* p, IWString_and_File_Descriptor&);

template <typename T>
int
write_as(const char*& p, int count, IWString_and_File_Descriptor& output) {
  for (int i = 0; i < count; i++) {
    write_as<T>(p, output);
    p += sizeof(T);
  }

  return 1;
}

template int
write_as<unsigned int>(const char*&, int, IWString_and_File_Descriptor&);
template int
write_as<int>(const char*&, int, IWString_and_File_Descriptor&);
template int
write_as<unsigned char>(const char*&, int, IWString_and_File_Descriptor&);

static int
write_data(const Dbt& zdata, const Unpack_Binary_Data& unpack_format, int writing_key,
           IWString_and_File_Descriptor& output) {
  if (writing_key) {
    if (!contains_binary_characters(zdata)) {
      output.write(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());
    } else if (unpack_format.active()) {
      unpack_format.write_unpacked_data(reinterpret_cast<const char*>(zdata.get_data()),
                                        zdata.get_size(), output);
    } else {
      output << zdata.get_size() << " binary";
    }
    return 1;
  }

  if (unpack_format.active()) {
    return unpack_format.write_unpacked_data(
        reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size(), output);
  }

  if (need_to_uncompress_with_zlib) {
    return uncompress_and_write(zdata, output);
  }

  if (print_binary_data || !contains_binary_characters(zdata)) {
    output.write(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());
  } else {
    output << zdata.get_size() << " bytes of binary data";
  }

  return output.good();
}

static int
list_berkeley_database(const Dbt& zkey,
                       Dbt& zdata,  // not const because it may get byte swapped
                       IWString_and_File_Descriptor& output) {
  if (verbose > 2) {
    cerr << "Processing '";
    cerr.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size());
    cerr << "'\n";
  }

  if (write_database_contents_first) {
    write_data(zdata, data_unpack_format, 0, output);
  } else if (write_database_key) {
    write_data(zkey, key_unpack_format, 1, output);
  }

  if (write_database_key && write_the_data) {
    output << ' ';
  }

  if (write_database_contents_first) {
    write_data(zkey, key_unpack_format, 1, output);
  } else if (write_the_data) {
    write_data(zdata, data_unpack_format, 0, output);
  }

  output << '\n';

  return output.good();
}

/*
  Does the associated data match the pattern for the data_rx - why is this a separate
  function, check it when it is looked up!!
*/

static int
data_rx_matches(Db& database, re2::RE2& data_rx, Dbt& zkey) {
  Dbt fromdb;
  fromdb.set_flags(DB_DBT_MALLOC);

  int rc = database.get(NULL, &zkey, &fromdb, 0);
  if (0 != rc) {
    cerr << "Yipes, could not fetch data for key '";
    cerr.write(reinterpret_cast<const char*>(zkey.get_data()), zkey.get_size()) << "'\n";
    database.err(rc, "");
    return 0;
  }

  re2::StringPiece tmp(reinterpret_cast<const char*>(fromdb.get_data()),
                       fromdb.get_size());
  rc = RE2::PartialMatch(tmp, data_rx);

  delete static_cast<char*>(fromdb.get_data());

  return rc;
}

static int
list_berkeley_database(Db& database, IWString_and_File_Descriptor& output) {
  Dbc* cursor = NULL;

  int rc = database.cursor(NULL, &cursor, 0);
  if (0 != rc) {
    database.err(rc, "cannot acquire cursor");
    return 0;
  }

  int records_in_database = 0;
  int records_matching = 0;

  Accumulator_Int<int64_t> bytes;
  extending_resizable_array<int> bytes_per_record;
  bytes_per_record.resize(200000);

  Dbt zkey, zdata;

  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    records_in_database++;

    if (report_progress()) {
      cerr << "Scanned " << records_in_database << " records\n";
    }

    if (records_in_database <= records_to_skip) {
      continue;
    }
    if (records_to_process > 0 && records_matching >= records_to_process) {
      break;
    }
    if (key_rx) {
      re2::StringPiece tmp(reinterpret_cast<const char*>(zkey.get_data()),
                           zkey.get_size());
      if (!RE2::PartialMatch(tmp, *key_rx)) {
        continue;
      }
    }

    if (data_rx && !data_rx_matches(database, *data_rx, zkey)) {
      continue;
    }

    bytes.extra(zdata.get_size());
    bytes_per_record[zdata.get_size()]++;

    records_matching++;

    if (output_needed) {
      (void)list_berkeley_database(zkey, zdata, output);
      output.write_if_buffer_holds_more_than(flush_when);
    }
  }

  if (records_to_process > 0) {
    ;
  } else if (DB_NOTFOUND != rc) {
    database.err(rc, "Strange error at end of cursor");
  }

  if (verbose) {
    cerr << records_in_database << " records in database, " << records_matching
         << " records matched\n";
    if (bytes.n() > 1) {
      cerr << "Records between " << bytes.minval() << " and " << bytes.maxval()
           << " bytes, ave " << static_cast<float>(bytes.average()) << "\n";
    }
  }

  cursor->close();

  return output.good();
}

static int
list_berkeley_database(const char* dbname, IWString_and_File_Descriptor& output) {
  Db database(NULL, DB_CXX_NO_EXCEPTIONS);

  // Set cache size here if needed

  int rc = database.open(NULL, dbname, NULL, DB_UNKNOWN, DB_RDONLY, 0);

  if (0 != rc) {
    cerr << "Cannot open database '" << dbname << "' '";
    database.err(rc, "");
    return 0;
  }

  if (verbose) {
    cerr << "Database '" << dbname << "'\n";
  }

  rc = list_berkeley_database(database, output);

  database.close(0);

  return rc;
}

int
list_berkeley_database(int argc, char** argv) {
  Command_Line cl(argc, argv, "vbp:P:wu:U:nf:o:kxD:aZh:r:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('b')) {
    print_binary_data = 1;
    if (verbose) {
      cerr << "Will print binary data\n";
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

  if (cl.option_present('P')) {
    const_IWSubstring p;
    cl.value('P', p);

    if (!iwre2::RE2Reset(data_rx, p)) {
      cerr << "Could not compile data regular expression '" << p << "'\n";
      return 13;
    }

    if (verbose) {
      cerr << "Will only display records for which data matches '" << data_rx->pattern()
           << "'\n";
    }
  }

  if (cl.option_present('w')) {
    write_database_contents_first = 1;
    if (verbose) {
      cerr << "Database contents written before key values\n";
    }
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

  if (cl.option_present('n')) {
    output_needed = 0;

    if (verbose) {
      cerr << "Output suppressed\n";
    }

    if (0 == verbose) {
      verbose = 1;
    }
  }

  if (cl.option_present('Z')) {
    need_to_uncompress_with_zlib = 1;

    if (verbose) {
      cerr << "Data will be uncompressed with zlib\n";
    }
  }

  if (cl.option_present('f')) {
    if (!cl.value('f', records_to_skip) || records_to_skip < 0) {
      cerr << "The records to skip option (-f) must be a whole positive number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will skip the first " << records_to_skip << " records\n";
    }
  }

  if (cl.option_present('o')) {
    if (!cl.value('o', records_to_process) || records_to_process < 0) {
      cerr << "The records to process option (-o) must be a whole positive number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will only print the first " << records_to_process << " records\n";
    }
  }

  if (cl.option_present('k')) {
    write_database_key = 0;

    if (verbose) {
      cerr << "Will suppress writing the database key\n";
    }
  }

  if (cl.option_present('x')) {
    write_the_data = 0;

    if (verbose) {
      cerr << "WIll suppress writing the data\n";
    }
  }

  if (!write_database_key && !write_the_data) {
    cerr << "No output\n";
    usage(8);
  }

  if (cl.option_present('D')) {
    print_binary_data_as_daylight_fingerprint_tag = cl.string_value('D');

    if (verbose) {
      cerr << "Data written as Daylight fingerprint form\n";
    }

    if (cl.option_present('a')) {
      byte_swap_binary_data = 1;

      if (verbose) {
        cerr << "Data will be byte swapped before being written in Daylight form\n";
      }
    }
  }

  if (cl.option_present('h')) {
    if (!cl.value('h', flush_when) || flush_when < 1) {
      cerr << "The max output buffer size option (-h) must be a whole ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will flush output buffers when they contain " << flush_when << " bytes\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "The report progress option (-r) must be a whole +ve number\n";
      usage(3);
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!list_berkeley_database(cl[i], output)) {
      rc = i + 1;
      break;
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = list_berkeley_database(argc, argv);

  return rc;
}
