/*
  Does a line-by-line diff on a two or more files
  Designed for large files where diff cannot handle them.
  Designed only for the case where a small number of
  lines are missing from each file.
  We maintain a read-ahead buffer, with a hash to
  quickly tell which identifiers are in the cache.
*/

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static int break_on_first_diff = 0;

static Report_Progress report_progress;

static int identifier_column = 0;

static int include_lines_numbers_in_output = 0;

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
  cerr << "Does a line by line diff across two or more files\n";
  cerr << "Useful for comparing very large files, when the main diffs are lines missing\n";
  cerr << " -s <n>         number of records to buffer - look-ahead buffer\n";
  cerr << " -c <col>       identifier in column <col>\n";
  cerr << " -b             break upon finding first difference\n";
  cerr << " -r <n>         report progress every <n> records\n";
  cerr << " -l             include line numbers in output\n";
  cerr << " -B <stem>      file name stem for records missing from each file\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  Basically a file with a read-ahead cache
*/

class File_Being_Compared : public iwstring_data_source {
 private:
  int _nbuffers;
  IWString* _buffer;
  IWString* _id;
  int _buf_ptr;
  IW_STL_Hash_Map_int _id_to_line;
  IWString _fname;
  int _read_ahead_records;
  bool _eof;

  IWString_and_File_Descriptor _stream_for_skipped;
  int _items_written_to_skipped_stream;

  //  private functions

  int _write_current_to_stream_for_skipped();
  void _increment_buf_ptr();

 public:
  File_Being_Compared();
  ~File_Being_Compared();

  int debug_print(std::ostream&) const;

  int allocate_buffer(int);

  const IWString&
  current_id() const {
    return _id[_buf_ptr];
  }

  const IWString&
  current_buffer() const {
    return _buffer[_buf_ptr];
  }

  const IWString&
  fname() const {
    return _fname;
  }

  int open_stream_for_skipped(const char*);

  bool in_buffer(const IWString& s) const;

  int fill_buffer();

  int mark_current_id_skipped();
  int all_remaining_ids_skipped();

  int advance();
  int advance_to_identifier(const IWString& s);

  void
  set_fname(const char* s) {
    _fname = s;
  }

  int report(std::ostream&) const;
};

File_Being_Compared::File_Being_Compared() {
  _nbuffers = 0;
  _buffer = NULL;
  _id = NULL;
  _buf_ptr = 0;
  _read_ahead_records = 0;
  _eof = false;

  _items_written_to_skipped_stream = 0;

  return;
}

File_Being_Compared::~File_Being_Compared() {
  if (NULL != _buffer) {
    delete[] _buffer;
  }

  if (NULL != _id) {
    delete[] _id;
  }

  return;
}

int
File_Being_Compared::allocate_buffer(int n) {
  if (NULL != _buffer) {
    delete[] _buffer;
  }

  if (NULL != _id) {
    delete[] _id;
  }

  _buffer = new IWString[n];
  _id = new IWString[n];

  _nbuffers = n;

  _buf_ptr = 0;

  return 1;
}

int
File_Being_Compared::report(std::ostream& os) const {
  os << "Report on '" << _fname << "', read " << lines_read() << " records.";
  if (_items_written_to_skipped_stream) {
    os << ' ' << _items_written_to_skipped_stream << " records not in other files.";
  }
  os << "\n";

  return 1;
}

void
File_Being_Compared::_increment_buf_ptr() {
  _buf_ptr++;
  if (_buf_ptr >= _nbuffers) {
    _buf_ptr = 0;
  }

  return;
}

int
File_Being_Compared::debug_print(std::ostream& output) const {
  output << "File_Being_Compared::debug_print:contains " << _nbuffers
         << " read-ahead buffers, id to line " << _id_to_line.size() << " items\n";

  int b = _buf_ptr;
  for (int i = 0; i < _nbuffers; ++i) {
    output << _id[b] << " : " << _buffer[b] << "\n";

    b++;
    if (b >= _nbuffers) {
      b = 0;
    }
  }

  for (auto i = _id_to_line.begin(); i != _id_to_line.end(); ++i) {
    output << (*i).first << " is item " << (*i).second << endl;
  }

  return 1;
}

int
File_Being_Compared::open_stream_for_skipped(const char* f) {
  if (!_stream_for_skipped.open(f)) {
    cerr << "File_Being_Compared::open_stream_for_skipped:cannot open '" << f << "'\n";
    return 0;
  }

  return 1;
}

bool
File_Being_Compared::in_buffer(const IWString& s) const {
  // auto f = _id_to_line.find(s);
  // cerr << "Looking for '" << s << "' in buffers, " << (f != _id_to_line.end()) << endl;
  return _id_to_line.find(s) != _id_to_line.end();
}

/*
  _buf_ptr points to the current record. The read-ahead buffer is always full, so
  when we fetch a new record, it must replace _buf_ptr
*/

int
File_Being_Compared::advance() {
  if (_eof)  // already encountered eof, do we still have buffered records to offer
  {
    if (_read_ahead_records <= 0) {
      return 0;
    }

    _increment_buf_ptr();
    _read_ahead_records--;
    return 1;
  }

  if (_id[_buf_ptr].length() > 0) {
    _id_to_line.erase(_id[_buf_ptr]);
  }

  if (!next_record(_buffer[_buf_ptr])) {
    _eof = true;
    return 0;
  }

  IWString id;
  if (!_buffer[_buf_ptr].word(identifier_column, id)) {
    cerr << "File_Being_Compared::advance:cannot extract identifier from column "
         << identifier_column << " in '" << _buffer[_buf_ptr] << "'\n";
    return 0;
  }

  _id[_buf_ptr] = id;

  _id_to_line[id] = _buf_ptr;

  // cerr << "Just set identifier '" << id << "' to be in slot " << _buf_ptr << " hash
  // contains " << _id_to_line.size() << " items\n";

  _increment_buf_ptr();

  return 1;
}

int
File_Being_Compared::fill_buffer() {
  _buf_ptr = 0;

  _read_ahead_records = _nbuffers;  // unless we get an eof below

  for (int i = 0; i < _nbuffers; i++) {
    if (advance()) {  // great, looking good..
      ;
    } else if (_eof)  // fewer records in file than read-ahead buffer
    {
      if (verbose > 1) {
        cerr << "File_Being_Compared::fill_buffer:only read " << i << " of " << _nbuffers
             << " read-ahead records\n";
      }

      _nbuffers = i;
      _read_ahead_records = i;
    } else  // how can this happen?
    {
      cerr << "File_Being_Compared::fill_buffer:only read " << i << " records from '"
           << _fname << "'\n";
      return 0;
    }
  }

  _buf_ptr = 1;  // we were populating the read-ahead buffer

  return 1;
}

int
File_Being_Compared::_write_current_to_stream_for_skipped() {
  _stream_for_skipped << _buffer[_buf_ptr] << "\n";
  _stream_for_skipped.write_if_buffer_holds_more_than(4096);

  _items_written_to_skipped_stream++;

  return 1;
}

int
File_Being_Compared::advance_to_identifier(const IWString& s) {
  while (s != _id[_buf_ptr]) {
    if (_stream_for_skipped.is_open()) {
      _write_current_to_stream_for_skipped();
    }

    if (!advance()) {
      return 0;
    }
  }

  return 1;
}

int
File_Being_Compared::mark_current_id_skipped() {
  if (!_stream_for_skipped.is_open()) {
    return 0;
  }

  return _write_current_to_stream_for_skipped();
}

int
File_Being_Compared::all_remaining_ids_skipped() {
  if (!_stream_for_skipped.is_open()) {
    return 0;
  }

  int rc = 0;
  while (_read_ahead_records > 0) {
    mark_current_id_skipped();
    _increment_buf_ptr();
    _read_ahead_records--;
    rc++;
  }

  return rc;
}

static bool
all_buffers_the_same(const File_Being_Compared* file, int nfiles) {
  for (int i = 1; i < nfiles; ++i) {
    if (file[0].current_buffer() != file[i].current_buffer()) {
      return false;
    }
  }

  return true;
}

static bool
all_files_pointing_to_same_identifier(File_Being_Compared* file, int nfiles, bool& eof) {
  eof = false;

  // Try what is hopefully the most common case

  IWString id0 = file[0].current_id();  // not a reference since we may change it below

  bool same_ids = 1;
  for (int i = 0; i < nfiles; ++i) {
    if (file[i].current_id() != id0) {
      same_ids = false;
      break;
    }
  }

  if (same_ids) {
    return true;
  }

  // Now things get difficult, the files are pointing at different identifiers

#ifdef DEBUG_ALL_FILES_PTNG_TO_SAME_IDENTIFIER
  cerr << "Files not aligned at record level\n";

  for (int i = 0; i < nfiles; i++) {
    file[i].debug_print(cerr);
  }
#endif

  bool all_have_id;  // scope here for efficiency

  while (1) {
    all_have_id = true;
    for (int i = 1; i < nfiles; ++i) {
      if (!file[i].in_buffer(id0)) {
        all_have_id = false;
        break;
      }
    }

#ifdef DEBUG_ALL_FILES_PTNG_TO_SAME_IDENTIFIER
    cerr << "Do all files have '" << id0 << "'? " << all_have_id << endl;
#endif

    if (all_have_id) {
#ifdef DEBUG_ALL_FILES_PTNG_TO_SAME_IDENTIFIER
      cerr << "All files had id '" << id0 << "'\n";
#endif

      for (int i = 1; i < nfiles; ++i) {
        file[i].advance_to_identifier(id0);
      }

      return 1;
    }

    //  Not all the other files have id0, advance file 0 and try again

    file[0].mark_current_id_skipped();

    if (!file[0].advance()) {
      eof = true;
      return 0;
    }

    id0 = file[0].current_id();
#ifdef DEBUG_ALL_FILES_PTNG_TO_SAME_IDENTIFIER
    cerr << "id0 advanced to '" << id0 << "'\n";
#endif
  }

  return 1;
}

static int
diff_line_by_line(File_Being_Compared* file, int nfiles, std::ostream& output) {
  for (int i = 0; i < nfiles; ++i) {
    if (!file[i].fill_buffer()) {
      cerr << "File '" << file[i].fname() << "' cannot fill buffer\n";
      return 0;
    }

    if (verbose > 2) {
      file[i].debug_print(cerr);
    }
  }

  int rc = 0;  // means no differences

  while (1) {
    bool eof;
    if (all_files_pointing_to_same_identifier(file, nfiles, eof)) {
      ;
    } else if (eof) {
      return 1;
    } else {
      return 0;
    }

    if (verbose > 2) {
      cerr << "Found common id '" << file[0].current_id() << "'\n";
    }

    if (report_progress()) {
      cerr << "Processed " << file[0].lines_read() << " records, found " << rc
           << " different records\n";
    }

    if (!all_buffers_the_same(file, nfiles)) {
      for (int i = 0; i < nfiles; i++) {
        if (i > 0) {
          output << '|';
        }
        if (include_lines_numbers_in_output) {
          output << file[i].lines_read() << '|';
        }
        output << file[i].current_buffer();
        ;
      }
      output << "\n";

      rc++;

      if (break_on_first_diff) {
        return rc;
      }
    }

    for (int i = 0; i < nfiles; ++i) {
      if (!file[i].advance()) {
        return rc;
      }
    }
  }
}

static int
diff_line_by_line(int argc, char** argv) {
  Command_Line cl(argc, argv, "vbr:s:c:lB:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('b')) {
    break_on_first_diff = 1;

    if (verbose) {
      cerr << "Will break on first different record encountered\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      return 2;
    }
  }

  int buffer_size = 1000;

  if (cl.option_present('s')) {
    if (!cl.value('s', buffer_size) || buffer_size < 1) {
      cerr << "The number of records to cache (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Each file will buffer '" << buffer_size << " records\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', identifier_column) || identifier_column < 1) {
      cerr << "The identifier column (-c) option must be a valid column number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Identifiers in column " << identifier_column << endl;
    }

    identifier_column--;
  }

  if (cl.option_present('l')) {
    include_lines_numbers_in_output = 1;

    if (verbose) {
      cerr << "Will include line numbers in output\n";
    }
  }

  int nfiles = cl.number_elements();

  if (0 == nfiles) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (nfiles < 2) {
    cerr << "Requires more than one file\n";
    usage(2);
  }

  IWString bstem;
  if (cl.option_present('B')) {
    cl.value('B', bstem);

    if (verbose) {
      cerr << "Missing records written to files with stem '" << bstem << "'\n";
    }
  }

  File_Being_Compared* file = new File_Being_Compared[nfiles];
  std::unique_ptr<File_Being_Compared[]> free_file(file);

  for (int i = 0; i < nfiles; i++) {
    if (!file[i].open(cl[i])) {
      cerr << "Cannot open '" << cl[i] << "'\n";
      return i + 1;
    }

    file[i].set_fname(cl[i]);
    file[i].allocate_buffer(buffer_size);

    if (bstem.length()) {
      IWString tmp(bstem);
      tmp << i;

      if (!file[i].open_stream_for_skipped(tmp.null_terminated_chars())) {
        cerr << "Cannot open stream for missed items '" << tmp << "'\n";
        return i + 1;
      }
    }
  }

  int differences_found = diff_line_by_line(file, nfiles, std::cout);

  std::cout.flush();

  if (verbose) {
    cerr << "Across " << nfiles << " files found " << differences_found
         << " differing records\n";
    for (int i = 0; i < nfiles; i++) {
      file[i].report(cerr);
    }
  }

  if (differences_found > 0) {
    return 1;
  } else {
    return 0;
  }
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = diff_line_by_line(argc, argv);

  return rc;
}
