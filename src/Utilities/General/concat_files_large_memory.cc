/*
  Concatentate several files
  Sept 99, add option for dealing with TDT files

  This version reads the content of all files into memory - fast, but
  memory intensive, but it can process any number of files at once.
*/

#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

static IWString missing_value = '.';

static int* missing_dataitem = NULL;

/*
  since the records may be ordered differently in each file, we need to determine
  up front, for each file, where is each identifier.
*/

class AFile : public IW_STL_Hash_Map<IWString, IWString> {
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

  int _active;

 public:
  AFile();

  const IWString&
  fname() const {
    return _fname;
  }

  void
  set_file_name(const IWString& f) {
    _fname = f;
  }

  void
  set_identifier_column(int c) {
    _identifier_column = c;
  }

  int initialise(const char*);

  int
  active() const {
    return _active;
  }

  void
  deactivate() {
    _active = 0;
  }

  int
  columns_in_file() const {
    return _columns_in_file;
  }

  int
  records_in_file() const {
    return _records_in_file;
  }

  int
  duplicate_identifiers() const {
    return _duplicate_identifiers;
  }

  const IWString&
  header() const {
    return _header;
  }

  const IWString&
  filename() const {
    return _fname;
  }

  const IWString&
  missing_values() const {
    return _missing_values;
  }

  int
  missing_records() const {
    return _missing_records;
  }

  int&
  missing_records() {
    return _missing_records;
  }

  int echo(const IWString&, IWString&);

  int write_unprocessed_identifiers(IWString_and_File_Descriptor& os) const;
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

static IWString_and_File_Descriptor stream_for_identifiers_not_written;

static int records_discarded = 0;

static int trim_leading_zeros_from_identifiers = 0;

static int translate_identifiers_to_lowercase = 0;

static int suppress_warning_messages = 0;

static int ignore_zero_length_files = 0;

static char input_separator = ' ';
static char output_separator = ' ';

/*
  We can make things easy for ourselves by skipping certain
  files - like log files and smiles files
*/

static std::unique_ptr<re2::RE2> files_to_skip;

void
do_trim_leading_zeros_from_identifiers(IWString& buffer) {
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

template <typename T>
int
get_next_token(const const_IWSubstring& buffer, T& token, int& i, const char sep) {
  if (' ' == sep) {
    return buffer.nextword(token, i);
  }

  return buffer.nextword_single_delimiter(token, i, sep);
}

/*
  There are two places in the code where we need to get an identifier
  from a buffer, and then maybe massaging the identifier
*/

static int
preprocess_to_identifier_and_data(const_IWSubstring& buffer, int identifier_column,
                                  IWString& zid, IWString& zdata) {
  if (identifier_column > 0) {
    int i = 0;
    const_IWSubstring token;
    for (int c = 0; c <= identifier_column; ++c) {
      if (!get_next_token(buffer, token, i, input_separator)) {
        return 0;
      }

      if (c > 0) {
        zdata << output_separator;
      }

      zdata << token;
    }
    //  cerr << "From '" << buffer << "' col " << identifier_column << " get '" << token
    //  << "'\n";
    zid = token;

    while (get_next_token(buffer, token, i, input_separator)) {
      zdata << output_separator;
      zdata << token;
    }
  } else {
    const int ndx = buffer.index(input_separator);

    if (ndx < 0) {
      return 0;
    }

    buffer.from_to(0, ndx - 1, zid);
    buffer.from_to(ndx + 1, buffer.length() - 1, zdata);

    //  cerr << "ID '" << zid << "' data '" << zdata << "'\n";
  }

  if (chop_ && zid.contains('_')) {
    zid.truncate_at_first('_');
  }

  if (trim_leading_zeros_from_identifiers) {
    do_trim_leading_zeros_from_identifiers(zid);
  }

  if (translate_identifiers_to_lowercase) {
    zid.to_lowercase();
  }

  assert(0 != zid.length());

  return 1;
}

AFile::AFile() {
  _missing_records = 0;
  _records_in_file = 0;
  _duplicate_identifiers = 0;

  _identifier_column = 0;

  _active = 1;

  return;
}

int
AFile::initialise(const char* fname) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    _active = 0;
    cerr << "AFile::initialise: cannot open '" << fname << "'\n";
    return 0;
  }

  if (!input.next_record(_header)) {
    _active = 0;
    if (ignore_zero_length_files > 1) {
      return 1;
    }

    cerr << "AFile::initialise: cannot read header record\n";
    return ignore_zero_length_files;
  }

  const_IWSubstring tmp(_header);
  IWString zid, zdata;
  preprocess_to_identifier_and_data(tmp, _identifier_column, zid, zdata);
  // cerr << "Split into '" << zid << "' and '" << zdata << "'\n";

  _header = zdata;

  if (' ' == input_separator) {
    _columns_in_file = _header.nwords();
  } else {
    _columns_in_file = _header.nwords_single_delimiter(input_separator);
  }

  if (!header_records_present) {
    input.push_record();
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    buffer.strip_leading_blanks();
    buffer.strip_trailing_blanks();

    if (0 == buffer.length()) {
      if (skip_blank_lines) {
        continue;
      }

      cerr << "Blank line at line " << input.lines_read() << " in '" << _fname << "'\n";
      return 0;
    }

    IWString zid, zdata;
    if (!preprocess_to_identifier_and_data(buffer, _identifier_column, zid, zdata)) {
      cerr << "AFile::initialise: error on line " << input.lines_read() << endl;
      return 0;
    }

    if (contains(zid)) {
      if (ignore_duplicate_identifiers > 1) {
        ;
      } else if (0 == _duplicate_identifiers) {
        cerr << "File '" << fname << "'\n";
      }

      _duplicate_identifiers++;

      if (ignore_duplicate_identifiers > 1) {  // multiple -g options, don't even warn
        ;
      } else if (1 == ignore_duplicate_identifiers) {  // just warn
        cerr << "AFile::initialise: duplicate identifier '" << zid << "'\n";
      } else  // fatal
      {
        cerr << "AFile::initialise:duplicate identifier '" << zid << "'\n";
        return 0;
      }
    } else {
      zdata.strip_trailing_blanks();
      (*this)[zid] = zdata;

      //    cerr << "Identifier '" << buffer << "' contains " << zdata <<endl;
    }
  }

  _fname = fname;
  _records_in_file = input.lines_read();

  _missing_values.resize(2 * _columns_in_file);
  for (int i = 0; i < _columns_in_file; i++) {
    _missing_values += ' ';
    _missing_values += missing_value;
  }

  // input.close();

  return 1;
}

int
AFile::echo(const IWString& id, IWString& output_buffer) {
  const IW_STL_Hash_Map<IWString, IWString>::const_iterator f = find(id);

  if (f == end()) {
    _missing_records++;
    return 0;
  }

  if (!output_buffer.ends_with(output_separator)) {
    output_buffer += output_separator;
  }

  const_IWSubstring zdata =
      (*f).second;  // not const, and not a reference - we may change it below

  output_buffer += zdata;

  // If we will later come back and do rows not yet processed, remove references to ID
  // Get rid of it now, because if processing fails now, it will also fail later

  if (all_identifiers || stream_for_identifiers_not_written.is_open()) {
    erase(id);
  }

  return 1;
}

static int
write_unwritten_identifier(IWString_and_File_Descriptor& os, const const_IWSubstring& id,
                           const IWString& fname) {
  os << id << " from " << fname << " not written\n";

  return os.good();
}

int
AFile::write_unprocessed_identifiers(IWString_and_File_Descriptor& os) const {
  IW_STL_Hash_Map<IWString, IWString>::const_iterator i;
  for (i = begin(); i != end(); i++) {
    const IWString& id = (*i).first;

    write_unwritten_identifier(os, id, _fname);

    os.write_if_buffer_holds_more_than(32768);
  }

  return os.good();
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
  cerr << " -u               truncate identifiers at first '_' char\n";
  cerr << " -a               write all identifiers (includes those not in first file)\n";
  cerr << " -M <missing>     missing value string (default " << missing_value << ")\n";
  cerr << " -d               skip duplicate identifiers in the first file (repeat to suppress warning)\n";
  cerr << " -g               ignore duplicate identifiers in files (repeat to suppress warning)\n";
  cerr << " -c <column>      identifier column (default 1)\n";
  cerr << " -z               trim leading zero's from identifiers\n";
  cerr << " -I               only write records for which identifier is present in every file\n";
  cerr << " -K <fname>       write identifiers discarded by -I option to <fname>\n";
  cerr << " -n               input files are NOT descriptor files - header records not special\n";
  cerr << " -k               skip blank lines in all files\n";
  cerr << " -s               ignore case when comparing identifiers\n";
  cerr << " -H <rx>          skip files whose names match regular expression <rx>. Use 'def' for default\n";
  cerr << " -q               quiet operation, suppress most warning messages\n";
  cerr << " -e               skip empty files (repeat to silently skip)\n";
  cerr << " -i <sep>         input file separator\n";
  cerr << " -L <fname>       process a list of files read from <fname>\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int nfiles = 0;

static int
write_unprocessed_identifiers_to_K_file_if_requested(AFile* files) {
  if (!stream_for_identifiers_not_written.is_open()) {
    return 1;
  }

  for (int i = 1; i < nfiles; i++) {
    if (files[i].active()) {
      files[i].write_unprocessed_identifiers(stream_for_identifiers_not_written);
    }
  }

  return stream_for_identifiers_not_written.good();
}

/*
  In processing the first file, for some reason we didn't write a record.
  Fill the -K file with the information
*/

static int
write_identifiers_not_written(IWString_and_File_Descriptor& os, AFile* files,
                              const IWString& id) {
  write_unwritten_identifier(os, id, files[0].fname());

  for (int i = 1; i < nfiles; i++) {
    if (!files[i].active()) {
      continue;
    }

    if (files[i].contains(id)) {
      write_unwritten_identifier(os, id, files[i].fname());
      files[i].erase(id);  // don't want that file to try writing it later
    }
  }

  return os.good();
}

/*
  We are processing a remaining file, MISSING will be a bunch of missing
  values for all the items from earlier files not present
*/

static int
do_all_identifiers(const IWString& missing, int& records_written, AFile* files,
                   const int z, IWString_and_File_Descriptor& output) {
  AFile& f = files[z];

  // IW_STL_Hash_Map<IWString, IWString>::const_iterator i;
  for (auto i = f.begin(); i != f.end(); i++) {
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
        if (verbose && !suppress_warning_messages) {
          cerr << "No data for '" << id << "' in '" << jfile.filename() << "'\n";
        }
        output += jfile.missing_values();
        missing_dataitem[j]++;
      }

      all_identifiers = 1;  // set now for all files (z+1)....
    }

    output << '\n';
    records_written++;

    output.write_if_buffer_holds_more_than(32768);
  }

  if (output.length()) {
    output.flush();
  }

  return 1;
}

static int
do_all_identifiers(int columns_so_far, int& records_written, AFile* files, int z,
                   IWString_and_File_Descriptor& output) {
  IWString tmp;
  tmp += output_separator;
  tmp += missing_value;

  IWString missing;
  missing.append(columns_so_far, tmp);

  // cerr << "Missing for file " << z << " is '" << missing << "'\n";

  return do_all_identifiers(missing, records_written, files, z, output);
}

/*
  We have read the entire first file and matched up everything from the
  other files. Now echo any identifiers which were in other files
*/

static int
do_all_identifiers(int columns_so_far, int records_written, AFile* files,
                   IWString_and_File_Descriptor& output) {
  for (int i = 1; i < nfiles; i++) {
    const AFile& f = files[i];

    if (verbose) {
      cerr << "File '" << f.fname() << "' has " << f.size() << " unwritten ecords\n";
    }

    if (0 == f.size()) {  // all identifiers in this file written
      ;
    } else if (!do_all_identifiers(columns_so_far, records_written, files, i, output)) {
      return 0;
    }

    columns_so_far += f.columns_in_file();
  }

  return 1;
}

static int
concat_files(iwstring_data_source& input, int identifier_column, AFile* files,
             IWString_and_File_Descriptor& output) {
  output.resize(36000);

  const_IWSubstring input_buffer;
  if (!input.next_record(input_buffer)) {
    cerr << "Very bad news, cannot read first header\n";
    return 0;
  }

  int columns = 0;

  if (all_identifiers) {
    if (' ' == input_separator) {
      columns = input_buffer.nwords() - 1;
    } else {
      columns = input_buffer.nwords_single_delimiter(input_separator) - 1;
    }
  } else {
    columns = 0;  // keep the compiler quiet about possibly uninitialised variables
  }

  // IWString output_buffer;
  // output_buffer.resize(30000);

  // If there are no header records present, then push this record back.

  if (!header_records_present) {
    input.push_record();
  } else {
    output << input_buffer;

    for (int i = 1; i < nfiles; i++) {
      //    cerr << "File " << i << " active " << files[i].active() << endl;
      //    cerr << files[i].header() << endl;

      if (files[i].active()) {
        output << output_separator << files[i].header();
      }
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(8196);
  }

  IWString_STL_Hash_Set ids_encountered;

  while (input.next_record(input_buffer)) {
    input_buffer.strip_trailing_blanks();
    input_buffer.strip_leading_blanks();

    if (0 == input_buffer.length()) {
      if (skip_blank_lines) {
        continue;
      }

      cerr << "Blank line in first file, cannot continue\n";
      return 0;
    }

    IWString zid, zdata;
    preprocess_to_identifier_and_data(input_buffer, identifier_column, zid, zdata);

    if (skip_duplicate_ids_in_first_file) {
      if (ids_encountered.contains(zid)) {
        if (skip_duplicate_ids_in_first_file < 2) {
          cerr << "Skipping duplicate identifier '" << zid << "' in first file\n";
        }
        continue;
      }

      ids_encountered.insert(zid);
    }

    int size_of_buffer_before_starting_current_record = output.length();

    output << zid << output_separator << zdata;

    int nmissing = 0;
    for (int i = 1; i < nfiles; i++) {
      if (!files[i].active()) {
        continue;
      }

      if (!files[i].echo(zid, output)) {
        if (verbose && !suppress_warning_messages) {
          cerr << "No data for '" << zid << "' in '" << files[i].filename() << "'\n";
        }
        output << files[i].missing_values();
        nmissing++;
      }
    }

    if (nmissing && only_write_records_when_all_files_have_data) {
      if (stream_for_identifiers_not_written.is_open()) {
        write_identifiers_not_written(stream_for_identifiers_not_written, files, zid);
      }

      output.resize_keep_storage(size_of_buffer_before_starting_current_record);

      records_discarded++;
      continue;
    }

    output << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  if (output.length()) {
    output.flush();
  }

  if (only_write_records_when_all_files_have_data) {
    return write_unprocessed_identifiers_to_K_file_if_requested(files);
  }

  if (all_identifiers) {
    return do_all_identifiers(columns, input.lines_read() - 1, files, output);
  }

  return 1;
}

static int
concat_files(const char* fname, int identifier_column, AFile* files,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  files[0].set_file_name(fname);  // files[0] never calls its initialise() method, so we
                                  // need to tell it its file name

  return concat_files(input, identifier_column, files, output);
}

static int
read_file_names(iwstring_data_source& input, resizable_array_p<IWString>& file_names) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (0 == buffer.length()) {
      continue;
    }

    file_names.add(new IWString(buffer));
  }

  return file_names.number_elements();
}

static int
read_file_names(const char* fname, resizable_array_p<IWString>& file_names) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot read list of files from '" << fname << "'\n";
    return 0;
  }

  return read_file_names(input, file_names);
}

static int
concat_files(int argc, char** argv) {
  Command_Line cl(argc, argv, "vuM:adgzc:IK:nkH:sqeL:i:o:");

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
        cerr << "Identifiers in column " << col << endl;
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

  if (cl.option_present('M')) {
    missing_value = cl.string_value('M');

    if (verbose) {
      cerr << "Missing value string '" << missing_value << "'\n";
    }
  }

  if (cl.option_present('d')) {
    skip_duplicate_ids_in_first_file = cl.option_count('d');
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

    if (!stream_for_identifiers_not_written.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for discarded identifiers '" << fname << "'\n";
      return 3;
    }

    only_write_records_when_all_files_have_data = 1;

    if (verbose) {
      cerr << "Will only write complete records, discarded identifiers written to '"
           << fname << "'\n";
    }
  }

  if (cl.option_present('H')) {
    const_IWSubstring h = cl.string_value('H');

    if ("def" == h) {
      const_IWSubstring pattern("\\.{smi,log,sdf}$");
      iwre2::RE2Reset(files_to_skip, pattern);

      if (verbose) {
        cerr << "Will skip files matching '" << files_to_skip->pattern() << "'\n";
      }
    } else if (!iwre2::RE2Reset(files_to_skip, h)) {
      cerr << "Invalid regular expression '" << h << "'\n";
      return 6;
    } else {
      if (verbose) {
        cerr << "Will skip files matching '" << h << "'\n";
      }
    }
  }

  if (cl.option_present('q')) {
    suppress_warning_messages = 1;

    if (verbose) {
      cerr << "Will suppress most warning messages\n";
    }
  }

  if (cl.option_present('e')) {
    ignore_zero_length_files = cl.option_count('e');

    if (verbose) {
      cerr << "Will ignore zero length files\n";
    }
  }

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');
    if (!char_name_to_char(i)) {
      cerr << "Invalid input separtor '" << i << "'\n";
      return 1;
    }

    input_separator = i[0];
    output_separator = input_separator;
  }

  if (cl.option_present('o')) {
    IWString i = cl.string_value('o');
    if (!char_name_to_char(i)) {
      cerr << "Invalid output separtor '" << i << "'\n";
      return 1;
    }

    output_separator = i[0];
  }

  nfiles = cl.number_elements();

  if (nfiles > 0 && cl.option_present('L')) {
    cerr << "Cannot specify both files on the command line and files via the -L option\n";
    usage(2);
  }

  resizable_array_p<IWString> file_names;

  if (cl.option_present('L')) {
    const char* l = cl.option_value('L');

    if (!read_file_names(l, file_names) || file_names.number_elements() < 2) {
      cerr << "Cannot determine list of files to read from '" << l << "'\n";
      return 3;
    }

    nfiles = file_names.number_elements();
  } else if (cl.number_elements() < 2) {
    cerr << "Must specify at least two files on the command line\n";
    usage(2);
  } else {
    file_names.resize(nfiles);
    for (int i = 0; i < nfiles; i++) {
      file_names.add(new IWString(cl[i]));
    }
  }

  missing_dataitem = new_int(nfiles);
  std::unique_ptr<int> free_missing_dataitem(missing_dataitem);

  // Make sure we have an identifier column for every input file

  if (identifier_columns.number_elements() < nfiles) {
    identifier_columns.extend(nfiles, identifier_columns.last_item());
  }

  // We will echo the first file and append the others. Establish data in the others

  AFile* files = new AFile[nfiles];  // number 0 is not used

  int active_files = 0;

  for (int i = 1; i < nfiles; i++) {
    if (files_to_skip && iwre2::RE2PartialMatch(*(file_names[i]), *files_to_skip)) {
      files[i].deactivate();
      if (verbose) {
        cerr << "Skipping file '" << *(file_names[i]) << "' due to -H match\n";
      }
      continue;
    }

    files[i].set_identifier_column(identifier_columns[i]);

    if (!files[i].initialise(file_names[i]->null_terminated_chars())) {
      cerr << "Cannot initialise file '" << file_names[i] << "'\n";
      return i;
    }

    if (verbose) {
      cerr << "File '" << *(file_names[i]) << "' contains " << files[i].records_in_file()
           << " records of " << files[i].columns_in_file() << " columns each";
      if (files[i].duplicate_identifiers()) {
        cerr << ", " << files[i].duplicate_identifiers() << " duplicates";
      }
      cerr << endl;
    }

    active_files++;
  }

  if (0 == active_files) {
    cerr << "No files to process!\n";
    return 4;
  }

  int rc = 0;

  IWString_and_File_Descriptor output(1);

  if (!concat_files(file_names[0]->null_terminated_chars(), identifier_columns[0], files,
                    output)) {
    rc = 1;
  }

  if (verbose) {
    for (int i = 1; i < nfiles; i++) {
      if (files[i].active()) {
        cerr << "File " << i << " '" << files[i].filename() << "' had "
             << files[i].missing_records() << " missing records\n";
      }
    }

    if (records_discarded) {
      cerr << "Discarded " << records_discarded << " records with incomplete data\n";
    }

    for (int i = 0; i < nfiles; i++) {
      cerr << missing_dataitem[i] << " dataitems missing from file " << i << '\n';
    }
  }

  output.flush();

  if (stream_for_identifiers_not_written.is_open()) {
    stream_for_identifiers_not_written.flush();
  }

  delete[] files;

  return rc;
}

int
main(int argc, char** argv) {
  int rc = concat_files(argc, argv);

  return rc;
}
