/*
  Fetch records from a BerkeleyDb database
*/

#include <ctype.h>
#include <stdlib.h>

#include <algorithm>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"

#include "db_cxx.h"
#include "zlib.h"

#include "re2/re2.h"

using std::cerr;

static int verbose = 0;

static const char* prog_name = NULL;

static int print_binary_data = 0;

static std::unique_ptr<re2::RE2> key_rx;

static int column_for_identifier = 0;

static int write_other_data = 0;

static int write_identifiers = 1;

static int single_file_for_each_item_retrieved = 0;

static IWString single_file_stem;

static int need_to_uncompress_with_zlib = 0;

static int translate_tabs = 0;

//  Sometimes we want everything fetched, but with a special string for
//  those items that are absent.

static IWString write_items_not_found;

static int write_database_contents_first = 0;

static int records_tried = 0;

static int records_retrieved = 0;

//  Special processing for some selimsteg applications.
//  Key is identifier, data is smiles + other name tokens

static int smiles_id_with_multiple_tokens = 0;

// One or more databases can be opened.
static Db** database = NULL;
static int number_databases = 0;

static int* found_in_database = NULL;

//  When dealing with 6 and 12 digit Lilly numbers, it may be convenient to
//  left pad identifiers read in

static int left_pad_length = 0;

static char left_pad_character = '0';

static int strip_leading_zeros = 0;

//  Horrible hack for ASL problem

static int strip_zeros_after_prefix_and_before_other_numbers = 0;

//  For descriptor databases, we may have records that are either bits or bytes

static int data_is_bytes = 0;

static int data_is_bits = 0;

static int write_newline_character = 1;

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

  cerr << "fetches records from a Berkeley DB database\n";
  cerr << prog_name << " -d dbname identifier_file\n";
  cerr << " -d <db>         database\n";
  cerr << " -B <file>       file for identifiers not in database\n";
  cerr << " -M <string>     write missing identifiers with normal output, but\n";
  cerr << "                 append <string>\n";
  cerr << " -c <column>     column for identifiers (default 1)\n";
  cerr << " -p <pattern>    pattern for identifiers\n";
  cerr << " -o              write complete records from identifier file (s)\n";
  cerr << " -K <key>        fetch individual key(s)\n";
  cerr << " -k              suppress writing the key\n";
  cerr << " -b              print binary data\n";
  cerr << " -w              write database contents before key values (smiles for example)\n";
  cerr << " -l <num>        left pad identifiers to <num> chars with '0' chars\n";
  cerr << " -z              strip leading 0's from identifiers on input\n";
  cerr << " -D              strip 0's from identifiers like 'ASL0000nnnn'\n";
  cerr << " -Z              stored data is compressed by gzip, uncompress\n";
  cerr << " -f <size>       user specified buffer for fetches\n";
  cerr << " -h <dir>        database environment home directory\n";
  cerr << " -S <stem>       individual file for each item retrieved\n";
  cerr << " -G sel          selimsteg type application, data contains extra name tokens\n";
  cerr << " -T              translate tabs to spaces on input\n";
  cerr << " -u              suppress addition of newline after data\n";
  cerr << " -v              verbose output\n";
  // clang-format on

  exit(rc);
}

static int
uncompress_and_write(IWString_and_File_Descriptor& output, const Dbt& fromdb) {
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
      cerr << "Not enough memory for expansion\n";
      return 0;
    } else if (Z_BUF_ERROR == rc) {
      cerr << "Buffer not large enough, input " << fromdb.get_size() << " tried "
           << buffer_size << '\n';
    } else if (Z_DATA_ERROR == rc) {
      cerr << "Data corrupted, tried " << (i + 1) << " times\n";
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
do_write_smiles_and_identifiers(const Dbt& fromdb, const const_IWSubstring& to_write,
                                IWString_and_File_Descriptor& output) {
  const_IWSubstring myfromdb(reinterpret_cast<const char*>(fromdb.get_data()),
                             fromdb.get_size());

  myfromdb.remove_leading_chars(' ');

  int i = myfromdb.index(' ');

  if (i < 0)  // data field just contains smiles
  {
    output << myfromdb << ' ' << to_write << '\n';
    return 1;
  }

  // The data field contains 'smiles token token token...'

  output.strncat(myfromdb, i);

  output << ' ' << to_write;

  myfromdb += i;

  output << myfromdb << '\n';

  return 1;
}

static int
do_single_file_for_each_item_retrieved(const Dbt& fromdb,
                                       const const_IWSubstring& to_write) {
  IWString fname;

  if (single_file_stem.length()) {
    fname << single_file_stem << records_retrieved;
  } else {
    fname = to_write;
    fname.truncate_at_first(' ');
  }

  IWString_and_File_Descriptor output;

  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  output.write(reinterpret_cast<const char*>(fromdb.get_data()), fromdb.get_size());

  return output.good();
}

static int
contains_binary_characters(const Dbt& zdata) {
  const char* c = reinterpret_cast<const char*>(zdata.get_data());
  for (unsigned int i = 0; i < zdata.get_size(); i++) {
    if (isprint(c[i])) {
      ;
    } else if ('\n' == c[i]) {
      ;
    } else if ('\t' == c[i]) {
      ;
    } else {
      return 1;
    }
  }

  return 0;
}

static int
write_database_contents_as_bytes(IWString_and_File_Descriptor& output,
                                 const Dbt& from_database) {
  const unsigned char* b =
      reinterpret_cast<const unsigned char*>(from_database.get_data());

  for (unsigned int i = 0; i < from_database.get_size(); i++) {
    output << ' ' << static_cast<int>(b[i]);
  }

  return output.good();
}

static int
write_database_contents_as_bits(IWString_and_File_Descriptor& output,
                                const Dbt& from_database) {
  IW_Bits_Base b;

  if (!b.construct_from_array_of_bits((const unsigned char*)from_database.get_data(),
                                      from_database.get_size() * IW_BITS_PER_BYTE)) {
    cerr << "Huh, could not construct bits object from array!!!\n";
    return 0;
  }

  IWString buffer;

  if (!b.append_string_form(buffer, '1', '0',
                            1))  // last 1 means include spaces between bits
  {
    cerr << "Yipes, b.append_string_form failed\n";
    b.debug_print(cerr);
    return 0;
  }

  output << buffer;

  return output.good();
}

static int
write_database_contents(const Dbt& from_database,
                        const const_IWSubstring& id,  // only needed for uncompress errors
                        IWString_and_File_Descriptor& output) {
  if (data_is_bits) {
    return write_database_contents_as_bits(output, from_database);
  }

  if (data_is_bytes) {
    return write_database_contents_as_bytes(output, from_database);
  }

  if (need_to_uncompress_with_zlib) {
    if (!uncompress_and_write(output, from_database)) {
      cerr << "Failed to uncompress data for '" << id << "'\n";
      return 0;
    }
    return 1;
  }

  // Otherwise just write it as text

  // cerr << "Writing " << from_database.get_size() << " bytes\n";

  output.write(reinterpret_cast<const char*>(from_database.get_data()),
               from_database.get_size());

  return output.good();
}

static int
write_data(const Dbt& d, const const_IWSubstring& to_write,
           IWString_and_File_Descriptor& output) {
  if (single_file_for_each_item_retrieved) {
    return do_single_file_for_each_item_retrieved(d, to_write);
  }

  if (smiles_id_with_multiple_tokens) {
    return do_write_smiles_and_identifiers(d, to_write, output);
  }

  if (data_is_bits || data_is_bytes || need_to_uncompress_with_zlib) {
    ;
  } else if (!print_binary_data && contains_binary_characters(d)) {
    if (verbose) {
      cerr << "Data for '" << to_write << "' " << d.get_size()
           << " bytes containing binary data\n";
    }
    return 1;
  }

  if (write_database_contents_first) {
    write_database_contents(d, to_write, output);
    if (to_write.length()) {
      output << ' ' << to_write;
    }
  } else {
    if (to_write.length()) {
      output << to_write << ' ';
    }
    write_database_contents(d, to_write, output);
  }

  // Be careful of newlines - the IW_DY_Fingerprint method adds a newline

  if (write_newline_character) {
    output << '\n';
  }

  return output.good();
}

static int
write_missing_identifier_info(const const_IWSubstring& identifier,
                              const const_IWSubstring& to_write,
                              IWString_and_File_Descriptor& stream_for_not_found) {
  if (to_write.length()) {
    stream_for_not_found << to_write << '\n';
  } else {
    stream_for_not_found << identifier << '\n';
  }

  stream_for_not_found.write_if_buffer_holds_more_than(32768);

  return stream_for_not_found.good();
}

/*
  Can we speed things up with a user specified buffer
*/

static unsigned char* mybuffer = NULL;
static int lenbuf = 0;

static int
iwdb_fetch_3(const const_IWSubstring& identifier, const const_IWSubstring& to_write,
             IWString_and_File_Descriptor& output,
             IWString_and_File_Descriptor& stream_for_not_found) {
  // cerr << "Looking up '" << identifier << "'\n";

  Dbt dbkey((void*)identifier.rawchars(), identifier.length());

  Dbt fromdb;
  if (NULL != mybuffer) {
    assert(lenbuf > 0);
    fromdb.set_data(mybuffer);
    fromdb.set_ulen(lenbuf);
  }

  int flags = 0;

  int found = 0;
  for (int i = 0; i < number_databases; i++) {
    int tmp = database[i]->get(NULL, &dbkey, &fromdb, flags);
    if (0 == tmp) {
      found_in_database[i]++;
      found = 1;
      break;
    }
  }

  if (!found) {
    if (verbose > 1) {
      cerr << "Did not find '" << identifier << "' in database\n";
    }

    if (write_items_not_found.length()) {
      if (write_database_contents_first) {
        output << write_items_not_found << ' ' << to_write << '\n';
      } else {
        output << to_write << ' ' << write_items_not_found << '\n';
      }

      output.write_if_buffer_holds_more_than(32768);

      return output.good();
    }

    if (stream_for_not_found.is_open()) {
      return write_missing_identifier_info(identifier, to_write, stream_for_not_found);
    }

    return 1;
  }

  records_retrieved++;

  return write_data(fromdb, to_write, output);
}

static int
iwdb_fetch_2(const const_IWSubstring& identifier, const const_IWSubstring& to_write,
             IWString_and_File_Descriptor& output,
             IWString_and_File_Descriptor& stream_for_not_found) {
  if (key_rx && !iwre2::RE2FullMatch(identifier, *key_rx)) {
    cerr << "Invalid key '" << identifier << "'\n";
    return 0;
  }

  return iwdb_fetch_3(identifier, to_write, output, stream_for_not_found);
}

/*
  These classes are all about efficiency, and these days, it may not even matter.
  Generally, the identifier will flow through as a const_IWSubstring, but if we want
  to add leading zero's, or remove zero's from the middle, we need to change the
  string, so we cannot use a const_IWSubstring.
  So, each of these objects retains a buffer that the const_IWSubstring can point to.
  Definitely not thread safe.
*/

class Zero_Padder {
 private:
  int _width;
  char* _zeros;
  char* _mybuffer;

 public:
  Zero_Padder();
  ~Zero_Padder();

  int
  initialise(int n, char c);

  int
  process(const_IWSubstring& s);
};

Zero_Padder::Zero_Padder() {
  _width = 0;
  _zeros = nullptr;
  _mybuffer = nullptr;

  return;
}

Zero_Padder::~Zero_Padder() {
  if (nullptr != _zeros) {
    delete[] _zeros;
  }
  if (nullptr != _mybuffer) {
    delete[] _mybuffer;
  }

  return;
}

int
Zero_Padder::initialise(int n, char c) {
  _width = n;

  if (nullptr != _zeros) {  // should not happen
    delete[] _zeros;
  }

  _zeros = new char[_width];
  std::fill_n(_zeros, _width, c);

  std::fill_n(_zeros, _width, '0');

  _mybuffer = new char[_width];

  return 1;
}

int
Zero_Padder::process(const_IWSubstring& s) {
  const int nzero = _width - s.length();

  strncpy(_mybuffer, _zeros, nzero);
  strncpy(_mybuffer + nzero, s.rawchars(), s.length());

  s.set(_mybuffer, _width);

  return nzero;
}

static Zero_Padder zero_padder;

class In_The_Middle_Zero_Remover {
 private:
  char* _buffer;

 public:
  In_The_Middle_Zero_Remover();
  ~In_The_Middle_Zero_Remover();

  int
  process(const_IWSubstring& s);
};

In_The_Middle_Zero_Remover::In_The_Middle_Zero_Remover() {
  _buffer = nullptr;

  return;
}

In_The_Middle_Zero_Remover::~In_The_Middle_Zero_Remover() {
  if (nullptr != _buffer) {
    delete[] _buffer;
  }

  return;
}

int
In_The_Middle_Zero_Remover::process(const_IWSubstring& id) {
  const auto n = id.length();

  for (int i = 0; i < n; ++i) {
    if (!isdigit(id[i])) {
      continue;
    }

    //  got our first digit

    if (0 == i) {  // string starts with 0, that's not what we process
      return 0;
    }

    if ('0' != id[i]) {  // not a zero, cannot process
      return 0;
    }

    //  Need to make sure there is a non-zero number following

    const auto f = std::find_if(id.data() + i + 1, id.cend(),
                                [](char c) { return c > '1' && c <= '9'; });
    cerr << "end? " << (f == id.cend()) << '\n';
    if (f == id.cend()) {  // no non-zero digit following
      return 0;
    }

    //  Now remove all these zero's

    if (nullptr != _buffer) {
      delete[] _buffer;
    }

    const int newlength = i + (id.cend() - f);

    _buffer = new char[newlength];

    std::copy_n(id.cbegin(), i, _buffer);
    std::copy(f, id.cend(), _buffer + i);

    id.set(_buffer, newlength);

    return 1;
  }

  return 0;
}

static In_The_Middle_Zero_Remover in_the_middle_zero_remover;

static int
iwdb_fetch_1(const const_IWSubstring& identifier, const const_IWSubstring& to_write,
             IWString_and_File_Descriptor& output,
             IWString_and_File_Descriptor& stream_for_not_found) {
  if (0 == left_pad_length) {
    return iwdb_fetch_2(identifier, to_write, output, stream_for_not_found);
  }

  if (identifier.length() >= left_pad_length) {
    return iwdb_fetch_2(identifier, to_write, output, stream_for_not_found);
  }

  const_IWSubstring tmp(identifier);

  zero_padder.process(tmp);

  return iwdb_fetch_2(tmp, to_write, output, stream_for_not_found);
}

static int
determine_identifier(const const_IWSubstring& buffer, const_IWSubstring& identifier,
                     const_IWSubstring& to_write) {
  if (!buffer.word(column_for_identifier, identifier)) {
    return 0;
  }

  // cerr << "Determining the identifier, result '" << identifier << "'\n";

  if (write_other_data) {
    to_write = buffer;  // write everything
  } else if (write_identifiers) {
    to_write = identifier;
  }

  if (strip_leading_zeros) {
    identifier.remove_leading_chars('0');
    if (0 == identifier.length()) {
      cerr << "Yipes, identifier was all leading 0's\n";
      return 0;
    }
  } else if (strip_zeros_after_prefix_and_before_other_numbers) {
    in_the_middle_zero_remover.process(identifier);
  }

  return 1;
}

static int
iwdb_fetch(iwstring_data_source& input, IWString_and_File_Descriptor& output,
           IWString_and_File_Descriptor& stream_for_not_found) {
  if (translate_tabs) {
    input.set_translate_tabs(1);
  }

  input.set_dos(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    records_tried++;

    const_IWSubstring identifier, to_write;

    if (!determine_identifier(buffer, identifier, to_write)) {
      cerr << "Cannot determine identifier from '" << buffer << "'\n";
      return 0;
    }

    if (!iwdb_fetch_1(identifier, to_write, output, stream_for_not_found)) {
      return 0;
    }
  }

  return 1;
}

static int
iwdb_fetch(const char* fname, IWString_and_File_Descriptor& output,
           IWString_and_File_Descriptor& stream_for_not_found) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return iwdb_fetch(input, output, stream_for_not_found);
}

int
iwdb_fetch(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:bwp:okK:c:aB:t:l:OM:S:h:C:zZf:G:TuD");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  number_databases = cl.option_count('d');

  if (0 == number_databases) {
    cerr << "Must specify the database to use via the -d option\n";
    usage(31);
  }

  int cache_size = 0;
  if (cl.option_present('C')) {
    if (!cl.value('C', cache_size) || cache_size < 100) {
      cerr << "Must specify a good cache size (-C)\n";
      usage(3);
    }
  }

  int rc;

  DbEnv dbenv(DB_CXX_NO_EXCEPTIONS);

  DbEnv* dbenv_ptr;

  if (cl.option_present('h')) {
    IWString zenv = cl.string_value('h');

    rc = dbenv.open(zenv.null_terminated_chars(),
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

  database = new Db*[number_databases];

  try {
    int env_flags = DB_CXX_NO_EXCEPTIONS;

    for (int i = 0; i < number_databases; i++) {
      database[i] = new Db(dbenv_ptr, env_flags);

      if (cache_size > 0) {
        int rc = database[i]->set_cachesize(0, cache_size, 1);
        if (0 != rc) {
          cerr << "Cannot set cache size to " << cache_size << " " << dbenv.strerror(rc)
               << "'\n";
          return 3;
        }

        cerr << "Cache set to " << cache_size << '\n';
      }

      IWString dbname = cl.string_value('d', i);

      int rc = database[i]->open(NULL, dbname.null_terminated_chars(), NULL, DB_UNKNOWN,
                                 DB_RDONLY, 0);

      //    cerr << "Database '" << dbname << "' opened, rc = " << rc << '\n';

      if (0 != rc) {
        cerr << "Cannot open database '" << dbname << "' '" << dbenv.strerror(rc)
             << "'\n";
        database[i]->err(rc, "Didn't open");
        dbenv.close(0);
        return 0;
      }

      if (verbose) {
        cerr << "database " << i << " is '" << dbname << "'\n";
      }
    }
  } catch (const char* s) {
    cerr << "Caught exception '" << s << "'\n";
    return 1;
  } catch (void* s) {
    cerr << "Caught default exception\n";
    return 3;
  }

  found_in_database = new_int(number_databases);
  std::unique_ptr<int[]> free_found_in_database(found_in_database);

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

  if (cl.option_present('l') && (cl.option_present('z') || cl.option_present('D'))) {
    cerr << "The -l and -z options are mutually exclusive\n";
    usage(4);
  }

  if (cl.option_present('l')) {
    if (!cl.value('l', left_pad_length) || left_pad_length < 1) {
      cerr << "The left pad length option (-l) must be a whole positive number\n";
      usage(7);
    }

    if (verbose) {
      cerr << "Will left pad identifiers input to " << left_pad_length << " characters\n";
    }

    zero_padder.initialise(left_pad_length, left_pad_character);
  }

  if (cl.option_present('z')) {
    strip_leading_zeros = 1;

    if (verbose) {
      cerr << "Leading 0's will be stripped from identifiers\n";
    }
  }

  if (cl.option_present('D')) {
    strip_zeros_after_prefix_and_before_other_numbers = 1;

    if (verbose) {
      cerr << "Embedded 0's after prefix and before numbers removed (ASL0000nnn)\n";
    }
  }

  if (cl.option_present('o')) {
    write_other_data = 1;
    if (verbose) {
      cerr << "Will print all data in the identifier file\n";
    }
  }

  if (cl.option_present('k')) {
    write_identifiers = 0;

    if (verbose) {
      cerr << "Will not write key values\n";
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

  if (cl.option_present('w')) {
    write_database_contents_first = 1;
    if (verbose) {
      cerr << "Database contents written before key values\n";
    }
  }

  if (cl.option_present('K')) {
    ;
  } else if (cl.number_elements()) {
    ;
  } else {
    cerr << "Specify either a file of keys or individual keys via the -K option\n";
    usage(2);
  }

  if (cl.option_present('f')) {
    if (!cl.value('f', lenbuf) || lenbuf < 2) {
      cerr << "The length of the user supplied buffer must be > 1 bytes\n";
      usage(5);
    }

    mybuffer = new unsigned char[lenbuf];
    if (NULL == mybuffer) {
      cerr << "Yipes, could not allocate " << lenbuf << " bytes for user buffer\n";
      return 4;
    }

    if (verbose) {
      cerr << "Using user supplied buffer " << lenbuf << " butes\n";
    }
  }

  if (cl.option_present('G')) {
    int i = 0;
    const_IWSubstring g;
    while (cl.value('G', g, i++)) {
      if ("sel" == g) {
        smiles_id_with_multiple_tokens = 1;
        if (verbose) {
          cerr << "Special processing for selimsteg type databases\n";
        }
      } else {
        cerr << "Unrecognised -G qualifier '" << g << "'\n";
        usage(4);
      }
    }
  }

  if (cl.option_present('Z')) {
    need_to_uncompress_with_zlib = 1;

    if (verbose) {
      cerr << "Data will be uncompressed with zlib\n";
    }
  }

  if (cl.option_present('M') && cl.option_present('B')) {
    cerr << "The write missing values (-M) and missing values to file (-B) options\n";
    cerr << "are mutually inconsistent\n";
    usage(5);
  }

  IWString_and_File_Descriptor stream_for_not_found;

  if (cl.option_present('B')) {
    IWString b = cl.string_value('B');
    if (!stream_for_not_found.open(b.null_terminated_chars())) {
      cerr << "Cannot open stream for not found identifiers '" << b << "'\n";
      return 15;
    }

    if (verbose) {
      cerr << "Identifiers not found written to '" << b << "'\n";
    }
  }

  if (cl.option_present('M')) {
    write_items_not_found = cl.string_value('M');

    if (verbose) {
      cerr << "Identifiers not in database written to output - append '"
           << write_items_not_found << "'\n";
    }
  }

  if (cl.option_present('t')) {
    const_IWSubstring t = cl.string_value('t');

    if ("bits" == t) {
      data_is_bits = 1;
      if (verbose) {
        cerr << "Will treat the data as bits\n";
      }
    } else if ("bytes" == t) {
      data_is_bytes = 1;
      if (verbose) {
        cerr << "Will treat the data as bytes\n";
      }
    } else {
      cerr << "Unrecognised type (-t) qualifier '" << t << "'\n";
      usage(7);
    }
  }

  if (cl.option_present('T')) {
    translate_tabs = 1;

    if (verbose) {
      cerr << "Will translate tabs in the input file\n";
    }
  }

  if (cl.option_present('u')) {
    write_newline_character = 0;

    if (verbose) {
      cerr << "Will suppress addition of newline characters on output\n";
    }
  }

  if (cl.option_present('S')) {
    single_file_for_each_item_retrieved = 1;

    single_file_stem = cl.string_value('S');

    if (verbose) {
      cerr << "File for each item retrieved\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  rc = 0;

  if (cl.option_present('K')) {
    int i = 0;
    const_IWSubstring k;
    while (cl.value('K', k, i++)) {
      const_IWSubstring to_write;
      if (write_identifiers) {
        to_write = k;
      }

      records_tried++;

      if (!iwdb_fetch_1(k, to_write, output, stream_for_not_found)) {
        rc = i;
        break;
      }
    }
  }

  if (cl.number_elements()) {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!iwdb_fetch(cl[i], output, stream_for_not_found)) {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose) {
    cerr << "Retrieved " << records_retrieved << " of " << records_tried << " records\n";

    if (number_databases > 1) {
      for (int i = 0; i < number_databases; i++) {
        cerr << found_in_database[i] << " items found in database " << i << '\n';
      }
    }
  }

  for (int i = 0; i < number_databases; i++) {
    database[i]->close(0);
    if (verbose) {
      //    database[i]->stat(stats,     // implement sometime
    }

    delete database[i];
  }

  delete database;

  dbenv.close(0);

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = iwdb_fetch(argc, argv);

  return rc;
}
