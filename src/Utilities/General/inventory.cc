/*
  Scans one or more inventory files for matches
*/

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

// static int negative_amount_skipped = 0;

static int lines_read = 0;

static int lines_written = 0;

static int include_amount = 0;

static int duplicates_suppressed = 0;

static int check_hash_set = 0;

static IW_STL_Hash_Set restricted;

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
  cerr << "Scans one or more inventory files for items meeting min\n";
  cerr << " -I <fname>     name of inventory file\n";
  cerr << " -c <amt>       minimum amount for corresponding -I file\n";
  cerr << " -u <units>     units for corresponding -I file\n";
  cerr << " -C <col>       column with amount\n";
  cerr << " -M <text>      must have text on each line\n";
  cerr << " -a             append the amount present to the output\n";
  cerr << " -h             check for duplicates even if only one file\n";
  cerr << " -g             input is fingerprint file\n";
  cerr << " -R <fname>     restricted identifier file\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
read_restricted_items(iwstring_data_source& input, IW_STL_Hash_Set& restricted) {
  IWString buffer;

  while (input.next_record(buffer)) {
    restricted.insert(buffer);
  }

  return restricted.size();
}

static int
read_restricted_items(const char* fname, IW_STL_Hash_Set& restricted) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open restricted file '" << fname << "'\n";
    return 1;  // we can operate without it
  }

  return read_restricted_items(input, restricted);
}

class Inventory_File {
 private:
  const char* _fname;
  iwstring_data_source _input;
  const float _min_amount;
  const int _col;
  const const_IWSubstring _units;
  const IWString _must_have;

  int _lines_read;
  int _lines_written;
  int _duplicates_suppressed;
  int _restricted_items_skipped;

  // private functions

  int _filter_record(const const_IWSubstring& buffer, IW_STL_Hash_Set& written,
                     IWString_and_File_Descriptor& output);

 public:
  Inventory_File(const char* fname, float min_amount, int col,
                 const const_IWSubstring& units, const IWString& must_have);

  int
  lines_read() const {
    return _lines_read;
  }

  int
  lines_written() const {
    return _lines_written;
  }

  int
  duplicates_suppressed() const {
    return _duplicates_suppressed;
  }

  int report(std::ostream& os) const;

  int
  good() const {
    return _input.is_open();
  };

  int filter(IW_STL_Hash_Set& written, IWString_and_File_Descriptor&);
};

Inventory_File::Inventory_File(const char* fname, float min_amount, int col,
                               const const_IWSubstring& units, const IWString& must_have)
    : _fname(fname),
      _min_amount(min_amount),
      _col(col),
      _units(units),
      _must_have(must_have) {
  _lines_read = 0;
  _lines_written = 0;
  _duplicates_suppressed = 0;
  _restricted_items_skipped = 0;

  if (!_input.open(fname)) {
    cerr << "Inventory_File::Inventory_File:cannot open '" << fname << "'\n";
  }

  return;
}

int
Inventory_File::filter(IW_STL_Hash_Set& written, IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (_input.next_record(buffer)) {
    if (!_filter_record(buffer, written, output)) {
      cerr << "Fatal error processing '" << buffer << "', line " << _input.lines_read()
           << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
Inventory_File::_filter_record(const const_IWSubstring& buffer, IW_STL_Hash_Set& written,
                               IWString_and_File_Descriptor& output) {
  _lines_read++;

  if (0 == _must_have.length()) {
    ;
  } else if (!buffer.contains(_must_have)) {
    return 1;
  }

  IWString id;
  int i = 0;

  if (!buffer.nextword(id, i)) {
    cerr << "Cannot extract identifier from '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring token;
  if (1 == _col) {
    if (!buffer.nextword(token, i)) {
      cerr << "Cannot extract amount from '" << buffer << "'\n";
      return 0;
    }
  } else if (!buffer.word(_col, token)) {
    cerr << "Cannot extract column " << (_col + 1) << " from '" << buffer << "'\n";
    return 0;
  }

  float amt;

  if (!token.numeric_value(amt) || amt < 0.0) {  // just ignore these
    return 1;
  }

  if (amt < _min_amount) {
    return 1;
  }

  if (!check_hash_set) {
    ;
  } else if (written.contains(id)) {
    _duplicates_suppressed++;
    return 1;
  } else if (restricted.contains(id)) {
    _restricted_items_skipped++;
    return 1;
  } else {
    written.insert(id);
  }

  output << id;

  if (include_amount) {
    output << ' ' << token << ' ' << _units;
  }

  output << '\n';

  _lines_written++;

  return 1;
}

int
Inventory_File::report(std::ostream& os) const {
  os << "Inventory_File::report:file '" << _fname << "' read " << _lines_read
     << " records, wrote " << _lines_written << ", duplicates_suppressed "
     << _duplicates_suppressed << '\n';

  if (verbose > 1) {
    cerr << " Column " << (_col + 1) << " must have '" << _must_have << "'\n";
  }

  return os.good();
}

static int
inventory(int argc, char** argv) {
  Command_Line cl(argc, argv, "vI:c:C:u:M:ahR:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('a')) {
    include_amount = 1;

    if (verbose) {
      cerr << "Will include the amount with any output\n";
    }
  }

  int nfiles = cl.option_count('I');

  if (0 == nfiles) {
    cerr << "Must specify one or more inventory files via the -I option\n";
    usage(3);
  }

  if (nfiles != cl.option_count('c') || nfiles != cl.option_count('u') ||
      nfiles != cl.option_count('C') || nfiles != cl.option_count('M')) {
    cerr << "Must specify same number of -u and -c and -C options as -I options "
         << nfiles << '\n';
    usage(5);
  }

  if (1 == nfiles) {
    check_hash_set = 0;
  } else {
    check_hash_set = 1;
  }

  if (cl.option_present('h')) {
    check_hash_set = 1;
    if (verbose) {
      cerr << "Duplicate check performed\n";
    }
  }

  if (cl.number_elements()) {
    cerr << "command line arguments ignored\n";
    usage(2);
  }

  if (cl.option_present('R')) {
    const char* r = cl.option_value('R');

    if (!read_restricted_items(r, restricted)) {
      cerr << "Cannot read restricted items from '" << r << "'\n";
    }

    if (verbose) {
      cerr << "Read " << restricted.size() << " restricted items from '" << r << "'\n";
    }
  }

  IW_STL_Hash_Set written;

  IWString_and_File_Descriptor output(1);

  int i = 0;
  IWString fname;
  while (cl.value('I', fname, i)) {
    if (!dash_s(fname.null_terminated_chars())) {
      cerr << "Missing or empty file '" << fname << "', cannot continue\n";
      return i + 1;
    }

    const_IWSubstring u = cl.string_value('u', i);

    float c;
    if (!cl.value('c', c, i) || c < static_cast<float>(0.0)) {
      cerr << "Invalid amount\n";
      usage(4);
    }

    int col;
    if (!cl.value('C', col, i) || col < 2) {
      cerr << "The column option (-C) must be a valid column number\n";
      usage(3);
    }

    col--;

    IWString must_have;
    cl.value('M', must_have, i);
    if ('.' == must_have) {
      must_have.resize(0);  // these aren't really regular expressions
    } else {
      must_have.gsub('.', ' ');
    }

    Inventory_File ifile(fname, c, col, u, must_have);
    if (!ifile.good()) {
      cerr << "Cannot initialise inventory file '" << fname << "'\n";
      return 3;
    }

    if (!ifile.filter(written, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return i + 1;
    }

    if (verbose) {
      ifile.report(cerr);
      lines_read += ifile.lines_read();
      lines_written += ifile.lines_written();
      duplicates_suppressed += ifile.duplicates_suppressed();
    }

    i++;
  }

  if (verbose && nfiles > 1) {
    cerr << "Read " << lines_read << " lines from " << nfiles << " files, wrote "
         << lines_written << '\n';
    if (check_hash_set) {
      cerr << "Suppressed " << duplicates_suppressed << " duplicates\n";
    }
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = inventory(argc, argv);

  return rc;
}
