/*
  Concatentate several files
  Sept 99, add option for dealing with TDT files
*/

#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>

using std::cerr;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "iwtokeniser.h"

static IWString missing_value = '.';

static int* missing_dataitem = nullptr;

/*
  It can be convenient to give a wildcard and have the programme know
  of various files to NOT include
*/

static std::unique_ptr<RE2> suffix_exclusion_list;

static char input_separator = ' ';
static char output_separator = ' ';

// We can impose a requirement that the output be tabular.
// Note that this is not robust in the case of files with differening
// input separators. The test is an nwords check on the output.
static uint32_t output_must_be_tabular = 0;

/*
  since the records may be ordered differently in each file, we need to determine
  up front, for each file, where is each identifier.
*/

class AFile : public iwstring_data_source, public IW_STL_Hash_Map<IWString, off_t> {
 private:
  IWString _fname;

  int _columns_in_file;

  int _records_in_file;

  IWString _header;

  int _identifier_column;

  //  What we write when we don't have data on an identifier

  IWString _missing_values;

  //  We can keep track of the number of times each file doesn't have a record

  int _missing_records;

  //  We also track the number of duplicate identifiers

  int _duplicate_identifiers;

  int _quoted_fields;

  char _input_separator;

  // If we are suppressing duplicate column names, we may only
  // be writing a subset of columns int he input.
  int* _write_column;

  //  private functions

  int _append_everything_but_identifier_column(const_IWSubstring& buffer,
                                               IWString& output_buffer) const;
  int _identify_duplicate_descriptor_names(
      IWString_STL_Hash_Set& descriptor_names_encountered);
  template <typename T>
  IWString _column_subset(const T& s, char output_separator) const;

 public:
  AFile();
  ~AFile();

  int active() const {
    return _fname.length();
  }

  const IWString& fname() const {
    return _fname;
  }

  void set_file_name(const IWString& f) {
    _fname = f;
  }

  void set_identifier_column(int c) {
    _identifier_column = c;
  }

  void set_input_separator(char c) {
    _input_separator = c;
  }

  int initialise(const char*);

  int columns_in_file() const {
    return _columns_in_file;
  }

  int records_in_file() const {
    return _records_in_file;
  }

  int duplicate_identifiers() const {
    return _duplicate_identifiers;
  }

  const IWString& header() const {
    return _header;
  }

  const IWString& filename() const {
    return _fname;
  }

  const IWString& missing_values() const {
    return _missing_values;
  }

  int missing_records() const {
    return _missing_records;
  }

  int& missing_records() {
    return _missing_records;
  }

  int append_header(IWString& output, const char output_separator) const;

  int echo(const IWString&, IWString&);

  int write_unprocessed_identifiers(std::ostream& os) const;

  int descriptor_names_are_unique(IWString_STL_Hash_Set& descriptor_names_encountered);
  int disambiguate_duplicate_descriptors(
      IWString_STL_Hash_Set& descriptor_names_encountered, char output_separator);
};

static int verbose = 0;

/*
  Sometimes identifiers will be of the form
  R(12345)_xxxxx
  and the really important part is the R(12345)
*/

static int chop_ = 0;

static int skip_blank_lines = 0;

/*
  If we keep track of the identifiers encountered in the first file, we
  can skip duplicates in the first file
*/

static int skip_duplicate_ids_in_first_file = 0;

/*
  We can ignore duplicate identifiers in subsequent files too
*/

static int ignore_duplicate_identifiers = 0;

/*
  By default, we only write records for those identifiers which are in the first
  file.
*/

static int all_identifiers = 0;

/*
  So this can be used with files that aren't descriptor files
  By default, header records are present and we process descriptor files
*/

static int header_records_present = 1;

/*
  Sometimes we want to only write records for which every file has data
*/

static int only_write_records_when_all_files_have_data = 0;

static std::ofstream stream_for_identifiers_not_written;

static int records_discarded = 0;

static int trim_leading_zeros_from_identifiers = 0;

static int translate_identifiers_to_lowercase = 0;

static int die_if_duplicate_descriptor_names = 0;

static int first_file_may_contain_duplicate_ids = 0;

static int suppress_duplicate_descriptors = 0;

static int disambiguate_duplicate_descriptors = 0;

void
do_trim_leading_zeros_from_identifiers(const_IWSubstring& buffer) {
  int nzero = 0;
  for (int i = 0; i < buffer.length(); i++) {
    if ('0' != buffer[i]) {
      break;
    }

    nzero++;
  }

  if (0 == nzero) {
    return;
  }

  if (buffer.length() == nzero) {
    cerr << "Cannot strip zeros from '" << buffer << "'\n";
    return;
  }

  if (nzero) {
    buffer.remove_leading_chars(nzero);
  }

  return;
}

/*
  There are two places in the code where we need to get an identifier
  from a buffer, and then maybe massaging the identifier
*/

static int
preprocess_to_identifier(const_IWSubstring& buffer, const int identifier_column,
                         const int quot, const char sep) {
  // cerr << "Preprocessing, identifier_column " << identifier_column << '\n';
  // cerr << "Initially '" << buffer << "'\n";

  if (identifier_column > 0) {
    IWTokeniser iwt(buffer);
    iwt.set_sep(sep);
    if (' ' != sep) {
      iwt.set_empty_fields_valid(1);
    }

    const_IWSubstring token;
    for (int c = 0; c <= identifier_column; ++c) {
      if (!iwt.next_token(token)) {
        return 0;
      }
    }
    //  cerr << "From '" << buffer << "' col " << identifier_column << " get '" << token
    //  << "'\n";
    buffer = token;
  } else {
    buffer.truncate_at_first(sep);  // just the identifier
  }

  if (chop_ && buffer.contains('_')) {
    buffer.truncate_at_first('_');
  }

  if (trim_leading_zeros_from_identifiers) {
    do_trim_leading_zeros_from_identifiers(buffer);
  }

  if (0 == buffer.length()) {
    cerr << "preprocess_to_identifier:zero length identifier, column "
         << identifier_column << " sep '" << sep << "'\n";
    return 0;
  }

  assert(0 != buffer.length());

  return 1;
}

AFile::AFile() {
  _missing_records = 0;
  _records_in_file = 0;
  _duplicate_identifiers = 0;

  _identifier_column = 0;

  _quoted_fields = 0;

  _input_separator = input_separator;

  _write_column = nullptr;

  return;
}

AFile::~AFile() {
  if (_write_column != nullptr) {
    delete[] _write_column;
  }
}

/*
  Recognise things like

  fname,sep=,col=
*/

static int
recognise_file_qualifiers(const char* s, char& sep, int& col, int& quot,
                          IWString& fname) {
  const_IWSubstring buffer(s);

  int i = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i, ',')) {
    //  cerr << "recognise_file_qualifiers:examining '" << token << "'\n";

    if (0 == fname.length()) {
      fname = token;
    } else if (token.starts_with("sep=")) {
      token.remove_leading_chars(4);
      IWString tmp = token;
      if (!char_name_to_char(tmp)) {
        cerr << "recognise_file_qualifiers:unrecognised input separator '" << token
             << "'\n";
        return 0;
      }

      sep = tmp[0];
    } else if (token.starts_with("col=")) {
      token.remove_leading_chars(4);
      if (!token.numeric_value(col) || col <= 0) {
        cerr << "recognise_file_qualifiers:invalid identifier column '" << token << "'\n";
        return 0;
      }
      col--;
    } else if (token.starts_with("quot")) {
      quot = 1;
    } else {
      cerr << "recognise_file_qualifiers:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  if (0 == fname.length())  // should never happen
  {
    cerr << "recognise_file_qualifiers:file name not specified\n";
    return 0;
  }

  return 1;
}

int
AFile::initialise(const char* fname) {
  IWString myfname;

  if (!recognise_file_qualifiers(fname, _input_separator, _identifier_column,
                                 _quoted_fields, myfname)) {
    return 0;
  }

  if (!iwstring_data_source::open(myfname.null_terminated_chars())) {
    cerr << "AFile::initialise: cannot open '" << myfname << "'\n";
    return 0;
  }

  // cerr << "AFile::initialise:from '" << fname << "' column " << _identifier_column << "
  // sep '" << _input_separator << "'\n";

  iwstring_data_source::set_dos(1);

  off_t offset = iwstring_data_source::tellg();  // probably 0

  if (!iwstring_data_source::next_record(_header)) {
    cerr << "AFile::initialise: cannot read header record\n";
    return 0;
  }

  IWTokeniser iwt(_header);
  iwt.set_sep(_input_separator);
  if (' ' != _input_separator) {
    iwt.set_empty_fields_valid(1);
  }
  if (_quoted_fields) {
    iwt.set_quoted_tokens(1);
  }

  const_IWSubstring token;

  _columns_in_file = 0;
  IWString new_header;

  for (; iwt.next_token(token); ++_columns_in_file) {
    if (_columns_in_file == _identifier_column) {
      continue;
    }

    if (new_header.length() > 0) {
      new_header += output_separator;
    }

    new_header << token;
  }

  _header = new_header;

  _columns_in_file--;

  if (header_records_present) {
    offset = iwstring_data_source::tellg();
  } else {
    push_record();
  }

  const_IWSubstring buffer;
  while (iwstring_data_source::next_record(buffer)) {
    buffer.strip_leading_blanks();
    buffer.strip_trailing_blanks();

    if (0 == buffer.length()) {
      cerr << "Blank line at line " << iwstring_data_source::lines_read() << " in '"
           << fname << "'\n";

      if (skip_blank_lines) {
        offset = iwstring_data_source::tellg();
        continue;
      }

      return 0;
    }

    if (!preprocess_to_identifier(buffer, _identifier_column, _quoted_fields,
                                  _input_separator)) {
      cerr << "Cannot extract identifier from '" << buffer << "' col "
           << (_identifier_column + 1) << '\n';
      return 0;
    }

    if (contains(buffer)) {
      if (0 == _duplicate_identifiers) {
        cerr << "File '" << fname << "'\n";
      }

      _duplicate_identifiers++;

      if (ignore_duplicate_identifiers > 1) {  // multiple -g options, don't even warn
        ;
      } else if (1 == ignore_duplicate_identifiers) {  // just warn
        cerr << "AFile::initialise: duplicate identifier '" << buffer << "'\n";
      } else  // fatal
      {
        cerr << "AFile::initialise: duplicate identifier '" << buffer << "'\n";
        return 0;
      }
    } else if (translate_identifiers_to_lowercase) {
      IWString tmp(buffer);
      tmp.to_lowercase();

      (*this)[tmp] = offset;
    } else {
      (*this)[buffer] = offset;

      if (verbose < 2) {
        ;
      } else if (2 == verbose && 0 == size() % 1000) {
        cerr << "Identifier '" << buffer << "' found at " << offset << ", size " << size()
             << '\n';
      } else {
        cerr << "Identifier '" << buffer << "' found at " << offset << '\n';
      }
    }

    offset = iwstring_data_source::tellg();
  }

  _fname = fname;
  _records_in_file = iwstring_data_source::lines_read();

  _missing_values.resize(2 * _columns_in_file);
  for (int i = 0; i < _columns_in_file; i++) {
    _missing_values += output_separator;
    _missing_values += missing_value;
  }

  return 1;
}

/*
  Append the columns in `buffer` to `output`.
  Skip over _identifier_column and any columns that are unset in _write_column.
*/

int
AFile::_append_everything_but_identifier_column(const_IWSubstring& buffer,
                                                IWString& output_buffer) const {
  IWTokeniser iwt(buffer);
  iwt.set_sep(_input_separator);
  if (_quoted_fields) {
    iwt.set_quoted_tokens(1);
  }
  if (' ' != _input_separator) {
    iwt.set_empty_fields_valid(1);
  }

  const_IWSubstring token;

  for (int col = 0; iwt.next_token(token); ++col) {
    if (col == _identifier_column) {
      continue;
    }

    if (_write_column != nullptr && _write_column[col] == 0) {
      continue;
    }

    output_buffer += output_separator;

    output_buffer << token;
  }

  return 1;
}

template <typename T>
int
append_possibly_translated(const T& s, const char inp_sep, const char out_sep,
                           const int quot, IWString& output) {
  if (inp_sep == out_sep) {
    output << s;
    return 1;
  }

  IWString tmp(s);

  if (!quot) {
    tmp.gsub(inp_sep, out_sep);

    output << tmp;

    return 1;
  }

  static const char dquote = '"';

  int in_quote = s[0] == dquote;

  const int n = s.length();

  for (int i = 0; i < n; ++i) {
    const char c = s[i];

    if (dquote == c) {
      in_quote = !in_quote;
    } else if (in_quote) {
      ;
    } else if (inp_sep == c) {
      tmp[i] = out_sep;
    }
  }

  output << tmp;

  return 1;
}

// Make a subset of the columns as governed by _write_column
template <typename T>
IWString
AFile::_column_subset(const T& buffer, char output_separator) const {
  IWString result;
  result.reserve(buffer.nwords_single_delimiter(_input_separator));

  int i = 0;
  IWString token;
  for (int col = 0; buffer.nextword_single_delimiter(token, i, _input_separator); ++col) {
    if (!_write_column[col]) {
      continue;
    }
    if (result.length() > 0) {
      result << output_separator;
    }
    result << token;
  }

  return result;
}

int
AFile::append_header(IWString& output, const char output_separator) const {
  // If no column subsetting, write everything.
  if (_write_column == nullptr) {
    return append_possibly_translated(_header, _input_separator, output_separator,
                                      _quoted_fields, output);
  }

  // Write a subset of columns.
  // Not dealing with quoted fields here
  output << _column_subset(_header, output_separator);
  return 1;
}

int
AFile::echo(const IWString& id, IWString& output_buffer) {
  // cerr << "AFile::echo '" << id << "'\n";

  IW_STL_Hash_Map<IWString, off_t>::const_iterator f = find(id);

  if (f != end()) {  // great, found it
    ;
  } else if (translate_identifiers_to_lowercase) {
    IWString tmp(id);
    tmp.to_lowercase();
    f = find(tmp);
  }

  if (f == end()) {
    _missing_records++;
    return 0;
  }

  off_t offset = f->second;

  // If we will later come back and do rows not yet processed, remove references to ID
  // Get rid of it now, because if processing fails now, it will also fail later

  if (first_file_may_contain_duplicate_ids) {
    ;
  } else if (all_identifiers || stream_for_identifiers_not_written.rdbuf()->is_open()) {
    erase(id);
  }

  if (!iwstring_data_source::seekg(offset)) {
    cerr << "AFile::echo: very bad news, cannot seek to " << offset << " for id '" << id
         << "'\n";
    return 0;
  }

  const_IWSubstring buffer;
  if (!iwstring_data_source::next_record(buffer)) {
    cerr << "AFile::echo: very bad news, cannot read from offset " << offset << '\n';
    return 0;
  }

  buffer.strip_trailing_blanks();

  // If _write_column is active, output must be done by looking at columns.

  if (_write_column != nullptr) {
    _append_everything_but_identifier_column(buffer, output_buffer);
  } else if (0 == _identifier_column) {                // the most common case.
    buffer.remove_leading_words(1, _input_separator);  // don't write the identifier again
    output_buffer << output_separator;
    append_possibly_translated(buffer, _input_separator, output_separator, _quoted_fields,
                               output_buffer);
  } else {
    _append_everything_but_identifier_column(buffer, output_buffer);
  }

  // cerr << "After appending data '" << buffer << "'\n";
  // cerr << "output_buffer is '" << output_buffer << "'\n";

  return 1;
}

static int
write_unwritten_identifier(const const_IWSubstring& id, const IWString& fname,
                           std::ostream& os) {
  os << id << " from " << fname << " not written\n";

  return os.good();
}

int
AFile::write_unprocessed_identifiers(std::ostream& os) const {
  IW_STL_Hash_Map<IWString, off_t>::const_iterator i;
  for (i = begin(); i != end(); i++) {
    const IWString& id = (*i).first;

    write_unwritten_identifier(id, _fname, os);
  }

  return os.good();
}

static int
all_descriptors_unique(const const_IWSubstring& buffer,
                       IWString_STL_Hash_Set& descriptor_names_encountered,
                       const char sep) {
  IWTokeniser iwt(buffer);
  iwt.set_sep(sep);
  if (' ' != sep) {
    iwt.set_empty_fields_valid(1);
  }

  IWString dname;

  while (iwt.next_token(dname)) {
    if (descriptor_names_encountered.contains(dname)) {
      cerr << "Duplicate descriptor '" << dname << "'\n";
      return 0;
    }

    descriptor_names_encountered.insert(dname);
  }

  return 1;
}

// Examine the contents of `_header` and determine if there are duplicate
// descriptor names, as governed by `descriptor_names_encountered` which is
// updated. Note that a duplicate within this file is not distinguished from
// a duplicate feature name that came from a previous file.
int
AFile::descriptor_names_are_unique(IWString_STL_Hash_Set& descriptor_names_encountered) {
  if (suppress_duplicate_descriptors) {
    return _identify_duplicate_descriptor_names(descriptor_names_encountered);
  }

  if (all_descriptors_unique(_header, descriptor_names_encountered, _input_separator)) {
    return 1;
  }

  cerr << "AFile::descriptor_names_are_unique:duplicate descriptor found in '" << _fname
       << "'\n";
  return 0;
}

int
AFile::_identify_duplicate_descriptor_names(
    IWString_STL_Hash_Set& descriptor_names_encountered) {
  _write_column = new_int(_header.nwords_single_delimiter(_input_separator));
  int duplicates_encountered = 0;
  int i = 0;
  IWString token;
  for (int col = 0; _header.nextword_single_delimiter(token, i, _input_separator);
       ++col) {
    if (descriptor_names_encountered.contains(token)) {
      _write_column[col] = 0;
      ++duplicates_encountered;
      continue;
    }
    descriptor_names_encountered.insert(token);
    _write_column[col] = 1;
  }

  // if no duplicates found, we wil be writing whole records.
  if (duplicates_encountered == 0) {
    delete[] _write_column;
    _write_column = nullptr;
  }

  if (verbose) {
    cerr << "Identified " << duplicates_encountered << " duplicate column names\n";
  }

  return 1;
}

// Return a variant on `name` that is not in `descriptor_names_encountered`.
// `descriptor_names_encountered` is updated once we generate an unused name.
IWString
NewName(const IWString& name, IWString_STL_Hash_Set& descriptor_names_encountered) {
  IWString new_name;
  for (int j = 2;; ++j) {
    new_name = name;
    new_name << '_' << j;
    if (!descriptor_names_encountered.contains(new_name)) {
      descriptor_names_encountered.insert(new_name);
      return new_name;
    }
  }
}

// Examine the contents of `_header` and if there are duplicate feature names in there
// reconstruct `_header` with updated names.
// returns the number of name changes made.
int
AFile::disambiguate_duplicate_descriptors(
    IWString_STL_Hash_Set& descriptor_names_encountered, char output_separator) {
  int rc = 0;

  int i = 0;
  IWString token;
  IWString new_header;

  new_header.reserve(_header.length() + 20);  // just a guess

  for (int col = 0; _header.nextword_single_delimiter(token, i, _input_separator);
       ++col) {
    if (!descriptor_names_encountered.contains(token)) {
      descriptor_names_encountered.insert(token);
      new_header.append_with_spacer(token, output_separator);
      continue;
    }

    new_header.append_with_spacer(NewName(token, descriptor_names_encountered),
                                  output_separator);
    ++rc;
  }

  _header = new_header;

  return 1;
}

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
  cerr << "Concatenates descriptor files by joining on identifiers\n";
  cerr << "           set per file settings with 'fname,sep=comma,col=4'\n";
  cerr << " -u               truncate identifiers at first '_' char\n";
  cerr << " -a               write all identifiers (includes those not in first file)\n";
  cerr << " -M <missing>     missing value string (default " << missing_value << ")\n";
  cerr << " -d               skip duplicate identifiers in the first file\n";
  cerr << " -f               first file may contain duplicate ID's\n";
  cerr << " -g               ignore duplicate identifiers in files\n";
  cerr << " -c <column>      identifier column(s) (default 1)\n";
  cerr << " -z               trim leading zero's from identifiers\n";
  cerr << " -I               only write records for which identifier is present in every file\n";
  cerr << " -K <fname>       write identifiers discarded by -I option to <fname>\n";
  cerr << " -n               input files are NOT descriptor files - header records not special\n";
  cerr << " -k               skip blank lines in all files\n";
  cerr << " -s               ignore case when comparing identifiers\n";
  cerr << " -D die           stop processing if duplicate descriptor names are encountered\n";
  cerr << " -D rm            remove duplicate descriptors\n";
  cerr << " -D disambiguate  assign new unique names to duplicate descriptors\n";
  cerr << " -i <sep>         input  file separator (default space)\n";
  cerr << " -o <sep>         output file separator (default space)\n";
  cerr << " -q               input consists of quoted fields\n";
  cerr << " -Y ...           other options, enter '-Y help' for info\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
write_unprocessed_identifiers_to_K_file_if_requested(
    AFile* files, const int nfiles, std::ofstream& stream_for_identifiers_not_written) {
  if (!stream_for_identifiers_not_written.rdbuf()->is_open()) {
    return 1;
  }

  for (int i = 1; i < nfiles; i++) {
    if (!files[i].active()) {
      continue;
    }

    files[i].write_unprocessed_identifiers(stream_for_identifiers_not_written);
  }

  return stream_for_identifiers_not_written.good();
}

/*
  In processing the first file, for some reason we didn't write a record.
  Fill the -K file with the information
*/

static int
write_identifiers_not_written(AFile* files, const int nfiles, const const_IWSubstring& id,
                              std::ostream& os) {
  write_unwritten_identifier(id, files[0].fname(), os);

  for (int i = 1; i < nfiles; i++) {
    if (!files[i].active()) {
      continue;
    }

    if (files[i].contains(id)) {
      write_unwritten_identifier(id, files[i].fname(), os);
      files[i].erase(id);  // don't want that file to try writing it later
    }
  }

  return os.good();
}

// There is a requirement that the output be tabular.
// `output` holds the last record written, starting at `start`.
// If `output_must_be_tabular` is uint32_t::max() then this must be the first
// call and we set it. Otherwise we make sure it is the same as what was first
// set.

static int
IsTabular(const IWString& output, uint32_t start,
          uint32_t& output_must_be_tabular) {
  const_IWSubstring s(output.data() + start, output.size() - start);
  const uint32_t nw = s.nwords_single_delimiter(output_separator);

  if (output_must_be_tabular == std::numeric_limits<uint32_t>::max()) {
    output_must_be_tabular = nw;
    return 1;
  } else if (output_must_be_tabular == nw) {
    return 1;
  }

  cerr << "Non tabular output found. Got " << nw << " words, expected " <<
          output_must_be_tabular << '\n';
  cerr << s;

  return 0;
}

/*
  We are processing a remaining file, MISSING will be a bunch of missing
  values for all the items from earlier files not present
*/

static int
do_all_identifiers(const IWString& missing, int& records_written, AFile* files,
                   const int nfiles, const int z, IWString_and_File_Descriptor& output) {
  AFile& f = files[z];

  // If we are checking for tabular output.
  uint32_t initial_size = output.size();

  IW_STL_Hash_Map<IWString, off_t>::const_iterator i;
  for (i = f.begin(); i != f.end(); i++) {
    const IWString& id = (*i).first;

    output << id << missing;

    //  We actively mess with all_identifiers to ensure that we don't delete entries
    //  from file Z. That would mess up our iterator. Yes, this is ugly.

    all_identifiers = 0;
    for (int j = z; j < nfiles; j++) {
      AFile& jfile = files[j];

      if (!jfile.active()) {
        continue;
      }

      if (!jfile.echo(id, output)) {
        if (verbose > 1) {
          cerr << "No offset for '" << id << "' in '" << jfile.filename() << "'\n";
        }
        output += jfile.missing_values();
        missing_dataitem[j]++;
      }

      all_identifiers = 1;  // set now for all files (z+1)....
    }

    output << '\n';
    records_written++;

    output.write_if_buffer_holds_more_than(4096);

    if (output_must_be_tabular && ! IsTabular(output, initial_size, output_must_be_tabular)) {
      return 0;
    }
  }

  return 1;
}

static int
do_all_identifiers(int columns_so_far, int& records_written, AFile* files,
                   const int nfiles, int z, IWString_and_File_Descriptor& output) {
  IWString tmp;
  tmp += output_separator;
  tmp += missing_value;

  IWString missing;
  missing.append(columns_so_far, tmp);

  // cerr << "Missing for file " << z << " is '" << missing << "'\n";

  return do_all_identifiers(missing, records_written, files, nfiles, z, output);
}

/*
  We have read the entire first file and matched up everything from the
  other files. Now echo any identifiers which were in other files
*/

static int
do_all_identifiers(int columns_so_far, int records_written, AFile* files,
                   const int nfiles, IWString_and_File_Descriptor& output) {
  for (int i = 1; i < nfiles; i++) {
    const AFile& f = files[i];

    if (!f.active()) {
      continue;
    }

    if (verbose) {
      cerr << "File '" << f.fname() << "' has " << f.size() << " unwritten ecords\n";
    }

    if (0 == f.size()) {  // all identifiers in this file written
      ;
    } else if (!do_all_identifiers(columns_so_far, records_written, files, nfiles, i,
                                   output)) {
      return 0;
    }

    columns_so_far += f.columns_in_file();
  }

  return 1;
}

static int
concat_files(iwstring_data_source& input, const char input_separator,
             const int identifier_column, const int quot, AFile* files, const int nfiles,
             IWString_and_File_Descriptor& output) {
  output.resize(30000);

  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Very bad news, cannot read first header\n";
    return 0;
  }

  int columns;

  if (all_identifiers) {
    columns = buffer.nwords(input_separator) - 1;
  } else {
    columns = 0;  // keep the compiler quiet about possibly uninitialised variables
  }

  // Do we need to examine the contents of the header records?
  if (die_if_duplicate_descriptor_names || suppress_duplicate_descriptors ||
      disambiguate_duplicate_descriptors) {
    IWString_STL_Hash_Set descriptor_names_encountered;

    if (all_descriptors_unique(buffer, descriptor_names_encountered, input_separator)) {
    } else if (die_if_duplicate_descriptor_names) {
      cerr << "Duplicate descriptor names present in first file header\n";
      return 0;
    }

    for (int i = 1; i < nfiles; ++i) {
      if (disambiguate_duplicate_descriptors) {
        files[i].disambiguate_duplicate_descriptors(descriptor_names_encountered,
                                                    output_separator);
        continue;
      }

      if (files[i].descriptor_names_are_unique(descriptor_names_encountered)) {
        continue;
      }

      if (die_if_duplicate_descriptor_names) {
        return 0;
      }
    }
  }

  // If there are no header records present, then push this record back.

  if (!header_records_present) {
    input.push_record();
  } else {
    append_possibly_translated(buffer, input_separator, output_separator, quot, output);

    for (int i = 1; i < nfiles; i++) {
      if (!files[i].active()) {
        continue;
      }

      output << output_separator;
      files[i].append_header(output, output_separator);
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  // now process the bulk of the file

  IWString_STL_Hash_Set ids_encountered;

  while (input.next_record(buffer)) {
    //  buffer.strip_trailing_blanks();     no, breaks with tab separated input
    //  buffer.strip_leading_blanks();

    if (0 == buffer.length()) {
      if (skip_blank_lines) {
        continue;
      }

      cerr << "Blank line in first file, cannot continue\n";
      return 0;
    }

    const_IWSubstring id = buffer;

    preprocess_to_identifier(id, identifier_column, quot, input_separator);

    if (skip_duplicate_ids_in_first_file) {
      if (ids_encountered.contains(id)) {
        cerr << "Skipping duplicate identifier '" << id << "' in first file\n";
        continue;
      }

      ids_encountered.insert(id);
    }

    uint32_t size_of_buffer_before_processing_id = output.size();

    append_possibly_translated(buffer, input_separator, output_separator, quot, output);

    int nmissing = 0;
    for (int i = 1; i < nfiles; i++) {
      if (!files[i].active()) {
        continue;
      }

      if (!files[i].echo(id, output)) {
        if (verbose > 1) {
          cerr << "No offset for '" << id << "' in '" << files[i].filename() << "'\n";
        }
        output << files[i].missing_values();
        nmissing++;
        missing_dataitem[i]++;
      }
    }

    if (nmissing) {
      if (stream_for_identifiers_not_written.rdbuf()->is_open()) {
        write_identifiers_not_written(files, nfiles, id,
                                      stream_for_identifiers_not_written);
      }

      if (only_write_records_when_all_files_have_data) {
        records_discarded++;
        output.resize_keep_storage(size_of_buffer_before_processing_id);
        continue;
      }
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(8192);

    if (output_must_be_tabular && ! IsTabular(output, size_of_buffer_before_processing_id,
                output_must_be_tabular)) {
      return 0;
    }
  }

  if (only_write_records_when_all_files_have_data) {
    return write_unprocessed_identifiers_to_K_file_if_requested(
        files, nfiles, stream_for_identifiers_not_written);
  }

  if (all_identifiers) {
    return do_all_identifiers(columns, input.lines_read() - 1, files, nfiles, output);
  }

  return 1;
}

static int
concat_files(const char* fname, int identifier_column, int quot, AFile* files,
             const int nfiles, IWString_and_File_Descriptor output) {
  IWString myfname;

  if (!recognise_file_qualifiers(fname, input_separator, identifier_column, quot,
                                 myfname)) {
    return 0;
  }

  iwstring_data_source input(myfname.null_terminated_chars());

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  files[0].set_file_name(fname);  // files[0] never calls its initialise() method, so we
                                  // need to tell it its file name

  return concat_files(input, input_separator, identifier_column, quot, files, nfiles,
                      output);
}

static void
DisplayDashYOptions() {
  cerr << " -Y tabular          fail if non-tabular output is generated\n";
  cerr << "                     Note that this uses the value of the -o option to count tokens\n";

  ::exit(0);
}

static int
concat_files(int argc, char** argv) {
  Command_Line cl(argc, argv, "vuM:adgzc:IK:nksX:qQbi:o:fD:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options present\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('z')) {
    trim_leading_zeros_from_identifiers = 1;

    if (verbose) {
      cerr << "Leading zero's will be trimmed from identifiers\n";
    }
  }

  if (cl.option_present('s')) {
    translate_identifiers_to_lowercase = 1;

    if (verbose) {
      cerr << "Will ignore case when comparing identifiers\n";
    }
  }

  resizable_array<int> identifier_columns;

  if (cl.option_present('c')) {
    const_IWSubstring c;
    int i = 0;
    while (cl.value('c', c, i++)) {
      int col;
      if (!c.numeric_value(col) || col < 1) {
        cerr << "Invalid identifier column specification '" << c << "'\n";
        usage(7);
      }

      if (verbose) {
        cerr << "Identifiers in column " << col << '\n';
      }

      col--;

      identifier_columns.add(col);
    }
  } else {
    identifier_columns.add(0);
  }

  if (cl.option_present('u')) {
    chop_ = 1;
    if (verbose) {
      cerr << "Will truncate identifiers at first '_' char\n";
    }
  }

  if (cl.option_present('k')) {
    skip_blank_lines = 1;

    if (verbose) {
      cerr << "Will skip blank lines\n";
    }
  }

  if (cl.option_present('a')) {
    all_identifiers = 1;
    if (verbose) {
      cerr << "Will write all identifiers\n";
    }
  }

  if (cl.option_present('b')) {
    die_if_duplicate_descriptor_names = 1;

    if (verbose) {
      cerr << "Will check for duplicate descriptor names\n";
    }
  }

  if (cl.option_present('D')) {
    const_IWSubstring d = cl.string_value('D');
    if (d == "die") {
      die_if_duplicate_descriptor_names = 1;
      if (verbose) {
        cerr << "Will terminate processing on encountering duplicate descriptor names\n";
      }
    } else if (d == "rm") {
      suppress_duplicate_descriptors = 1;
      if (verbose) {
        cerr << "Duplicate descriptors suppressed\n";
      }
    } else if (d == "disambiguate") {
      disambiguate_duplicate_descriptors = 1;
      if (verbose) {
        cerr << "Duplicate descriptors will be disambiguated\n";
      }
    } else {
      cerr << "Unrecognised -D qualifier '" << d << "'\n";
      usage(1);
    }
  }

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');
    if (!char_name_to_char(i)) {
      cerr << "Invalid input separtor '" << i << "'\n";
      return 1;
    }

    input_separator = i[0];
  }

  if (cl.option_present('o')) {
    IWString i = cl.string_value('o');
    if (!char_name_to_char(i)) {
      cerr << "Invalid output separtor '" << i << "'\n";
      return 1;
    }

    output_separator = i[0];
  }

  if (cl.option_present('f')) {
    if (all_identifiers) {
      cerr << "Sorry, cannot use the -a and -f options together\n";
      return 1;
    }

    first_file_may_contain_duplicate_ids = 1;

    if (verbose) {
      cerr << "First file may contain duplicate ID's\n";
    }
  }

  if (cl.option_present('M')) {
    missing_value = cl.string_value('M');

    if (verbose) {
      cerr << "Missing value string '" << missing_value << "'\n";
    }
  }

  if (cl.option_present('d')) {
    skip_duplicate_ids_in_first_file = 1;
    if (verbose) {
      cerr << "Will skip duplicate identifiers in the first file\n";
    }
  }

  if (cl.option_present('g')) {
    ignore_duplicate_identifiers = cl.option_count('g');

    if (verbose) {
      cerr << "Will ignore duplicate identifiers\n";
    }
  }

  if (cl.option_present('I')) {
    only_write_records_when_all_files_have_data = 1;

    if (verbose) {
      cerr
          << "Will only write records when each input file has data for the identifier\n";
    }
  }

  if (all_identifiers && only_write_records_when_all_files_have_data) {
    cerr << "The do all identifiers (-a) and only do identifiers with complete data (-I) "
            "options are mutually exclusive\n";
    usage(6);
  }

  if (cl.option_present('n')) {
    header_records_present = 0;

    if (verbose) {
      cerr << "Header records not present\n";
    }
  }

  if (cl.option_present('K')) {
    IWString fname = cl.string_value('K');

    stream_for_identifiers_not_written.open(fname.null_terminated_chars(), std::ios::out);

    if (!stream_for_identifiers_not_written.good()) {
      cerr << "Cannot open stream for discarded identifiers '" << fname << "'\n";
      return 3;
    }

    only_write_records_when_all_files_have_data = 1;

    if (verbose) {
      cerr << "Will only write complete records, discarded identifiers written to '"
           << fname << "'\n";
    }
  }

  int nfiles = cl.number_elements();

  if (0 == nfiles) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (nfiles < 2) {
    cerr << "Just one file specified, will function like cat\n";
  }

  if (cl.option_present('X')) {
    const_IWSubstring x = cl.string_value('X');

    if ("." == x) {
      suffix_exclusion_list.reset(new RE2("\\.(smi|log)$"));
    } else {
      if (!iwre2::RE2Reset(suffix_exclusion_list, x)) {
        cerr << "Invalid suffix exclusion regular expression '" << x << "'\n";
        return 4;
      }
    }

    if (verbose) {
      cerr << "Files matching '" << suffix_exclusion_list->pattern()
           << "' will be skipped\n";
    }

    int valid_files_on_command_line = 0;
    int files_excluded = 0;
    for (int i = 0; i < cl.number_elements(); i++) {
      if (RE2::PartialMatch(cl[i], *suffix_exclusion_list)) {
        files_excluded++;
      } else {
        valid_files_on_command_line++;
      }
    }

    if (valid_files_on_command_line < 2) {
      cerr << "Sorry, only " << valid_files_on_command_line
           << " files can be processed\n";
      return 4;
    }

    if (verbose) {
      cerr << files_excluded << " files excluded, " << valid_files_on_command_line
           << " files remain\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "tabular") {
        output_must_be_tabular = std::numeric_limits<uint32_t>::max();
        if (verbose) {
          cerr << "Will fail if non tabular output is encountered\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions();
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions();
      }
    }
  }

  missing_dataitem = new_int(nfiles);
  std::unique_ptr<int> free_missing_dataitem(missing_dataitem);

  // Make sure we have an identifier column for every input file

  if (identifier_columns.number_elements() < nfiles) {
    identifier_columns.extend(nfiles, identifier_columns.last_item());
  }

  // We will echo the first file and append the others. Establish offsets in the others

  AFile* files = new AFile[nfiles];  // number 0 is not used

  for (int i = 1; i < nfiles; i++) {
    if (!suffix_exclusion_list) {
      ;
    } else if (RE2::PartialMatch(cl[i], *suffix_exclusion_list)) {
      continue;
    }

    files[i].set_identifier_column(identifier_columns[i]);
    files[i].set_input_separator(input_separator);

    if (!files[i].initialise(cl[i])) {
      cerr << "Cannot initialise file '" << cl[i] << "'\n";
      return i;
    }

    if (verbose) {
      cerr << "File '" << cl[i] << "' contains " << files[i].records_in_file()
           << " records of " << files[i].columns_in_file() << " columns each";
      if (files[i].duplicate_identifiers()) {
        cerr << ", " << files[i].duplicate_identifiers() << " duplicates";
      }
      cerr << '\n';
    }
  }

  int quot = 0;

  if (cl.option_present('q')) {
    quot = 1;

    if (verbose) {
      cerr << "Input may contain quoted fields\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  if (!concat_files(cl[0], identifier_columns[0], quot, files, nfiles, output)) {
    return 5;
  }

  output.flush();

  if (verbose) {
    int items_missing = 0;

    for (int i = 1; i < nfiles; i++) {
      if (!files[i].active()) {
        continue;
      }

      items_missing += files[i].missing_records();
    }

    if (items_missing) {
      for (int i = 1; i < nfiles; i++) {
        if (!files[i].active()) {
          continue;
        }

        cerr << "File " << i << " '" << files[i].filename() << "' had "
             << files[i].missing_records() << " missing records\n";
      }
    }

    if (records_discarded) {
      cerr << "Discarded " << records_discarded << " records with incomplete data\n";
    }

    items_missing = 0;

    for (int i = 0; i < nfiles; i++) {
      if (!files[i].active()) {
        continue;
      }

      items_missing += missing_dataitem[i];
    }

    if (items_missing) {
      for (int i = 0; i < nfiles; i++) {
        if (!files[i].active()) {
          continue;
        }

        cerr << missing_dataitem[i] << " dataitems missing from file " << i << ", '"
             << files[i].fname() << "'\n";
      }
    } else {
      cerr << "All items found in all files\n";
    }
  }

  if (cl.option_present('Q')) {
    return 0;
  }

  delete[] files;

  return 0;
}

int
main(int argc, char** argv) {
  int rc = concat_files(argc, argv);

  return rc;
}
