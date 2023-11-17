/*
  Ensures that the rows in two descriptor files match up
  Most common use will be to produce input files for various learners,
  that need the data from a descriptor file, but not the identifiers in
  column 1
*/

#include <stdlib.h>

#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static char input_separator = ' ';

static int remove_identifier = 0;

static int remove_header_records = 0;

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
  cerr << "Ensures the rows in two descriptor files match up\n";
  cerr << " -L <fname>     file with the desired row ordering\n";
  cerr << " -X <fname>     write the data (identifier column removed) from the -L file\n";
  cerr << " -r             remove the identifier from the output file\n";
  cerr << " -f             remove the header records\n";
  cerr << " -i <sep>       input separator (space, tab...)\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
nextword(const const_IWSubstring& buffer, int& i, IWString& s) {
  if (' ' == input_separator) {
    return buffer.nextword(s, i);
  } else {
    return buffer.nextword_single_delimiter(s, i, input_separator);
  }
}

static int
extract_id(const const_IWSubstring& buffer, int col, IWString& id) {
  int i = 0;
  if (!nextword(buffer, i, id)) {
    return 0;
  }

  if (0 == col) {
    return id.length();
  }

  while (col > 0) {
    if (!nextword(buffer, i, id)) {
      return 0;
    }
  }

  return id.length();
}

static int
store_data(IW_STL_Hash_Map_String& id_to_data,
           const_IWSubstring buffer,  // note local copy
           const int col) {
  IWString id;
  int i = 0;
  buffer.nextword(id, i, input_separator);

  buffer.remove_leading_chars(id.length() + 1);

  id_to_data.emplace(id, buffer);

  return 1;
}

#ifdef IMPLEMENT_SOMETIME
static int
do_write(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output) {
  return 1;
}
#endif

/*
  this is all kind of broken. If COL is anything but the first column, it will break
  For now, that is the only use, so we can ignore for now. Fix if it ever becomes
  problematic
*/

static int
descriptor_file_same_row_order(iwstring_data_source& input, const int col,
                               IW_STL_Hash_Map_int& id_to_ndx, IWString* ndx_to_id,
                               IWString_and_File_Descriptor& output) {
  IWString id;  // scope here for efficiency

  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Empty file\n";
    return 0;
  }

  if (!remove_header_records) {
    if (remove_identifier) {
      buffer.remove_leading_words(1, input_separator);
    }

    output << buffer << '\n';
  }

  const int n = id_to_ndx.size();

  int same_order = 1;

  while (input.next_record(buffer)) {
    if (!extract_id(buffer, col, id)) {
      return 0;
    }

    if (ndx_to_id[input.lines_read()] != id) {
      same_order = 0;
      break;
    }

    if (remove_identifier) {
      buffer.remove_leading_words(1, input_separator);
      output << buffer << '\n';
    } else {
      output << buffer << '\n';
    }
  }

  if (same_order) {
    if (verbose) {
      cerr << "Files in same order\n";
    }
    return 1;
  }

  if (verbose) {
    cerr << "Out of order detected at " << input.lines_read() << " records read\n";
  }

  int istart = input.lines_read();

  // Now things get hard. We need to read the rest of the file

  IW_STL_Hash_Map_String id_to_data;

  if (!store_data(id_to_data, buffer, col)) {
    return 0;
  }

  while (input.next_record(buffer)) {
    if (!store_data(id_to_data, buffer, col)) {
      return 0;
    }
  }

  for (int i = istart; i <= n; ++i) {
    const IWString& id = ndx_to_id[i];

    const auto f = id_to_data.find(id);

    if (f == id_to_data.end()) {
      cerr << "Unmatched identifier in input file '" << id << "'\n";
      return 0;
    }

    output << f->second << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  return output.good();
}

static int
descriptor_file_same_row_order(const char* fname, const int col,
                               IW_STL_Hash_Map_int& id_to_ndx, IWString* ndx_to_id,
                               IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptor_file_same_row_order(input, col, id_to_ndx, ndx_to_id, output);
}

static int
read_desired_row_ordering(iwstring_data_source& input, IW_STL_Hash_Map_int& id_to_ndx,
                          IWString* ndx_to_id, const int col,
                          IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Cannot read header\n";
    return 0;
  }

  if (remove_header_records) {
    ;
  } else if (output.is_open()) {
    buffer.remove_leading_words(1, input_separator);
    output << buffer << '\n';
  }

  ndx_to_id[0] = "ID";
  id_to_ndx["ID"] = 0;

  while (input.next_record(buffer)) {
    IWString id;
    if (!extract_id(buffer, col, id)) {
      cerr << "Cannot extract id '" << buffer << "' col " << (col + 1) << '\n';
      return 0;
    }

    id_to_ndx[id] = input.lines_read();
    ndx_to_id[input.lines_read()] = id;

    if (output.is_open()) {
      buffer.remove_leading_words(1, input_separator);
      output << buffer << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  }

  if (output.is_open()) {
    output.flush();
  }

  return id_to_ndx.size();
}

static int
read_desired_row_ordering(const char* fname, IW_STL_Hash_Map_int& id_to_ndx,
                          IWString*& ndx_to_id, const int col,
                          IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  const int n = input.records_remaining();  // don't worry about header

  if (0 == n) {
    cerr << "Empty file\n";
    return 0;
  }

  ndx_to_id = new IWString[n + 1];  // don't forget the header

  return read_desired_row_ordering(input, id_to_ndx, ndx_to_id, col, output);
}

static int
descriptor_file_same_row_order(int argc, char** argv) {
  Command_Line cl(argc, argv, "vL:c:C:i:X:rf");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('f')) {
    remove_header_records = 1;

    if (verbose) {
      cerr << "Will remove header records from output file(s)\n";
    }
  }

  if (!cl.option_present('L')) {
    cerr << "Must specify file with desired line ordering via the -L option\n";
    usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i)) {
      cerr << "Invalid input separator (-i) '" << i << "'\n";
      return 1;
    }
    input_separator = i[0];
  }

  if (cl.size() > 1) {
    cerr << "Do not know how to handle multiple input files\n";
    return 2;
  }

  IW_STL_Hash_Map_int id_to_ndx;
  IWString* ndx_to_id = nullptr;

  if (cl.option_present('L')) {
    int lcol = 0;

    if (cl.option_present('C')) {
      if (!cl.value('C', lcol) || lcol < 1) {
        cerr << "Invalid like file column number (-C)\n";
        usage(1);
      }

      if (verbose) {
        cerr << "Identifiers in like file in column " << lcol << '\n';
      }

      lcol--;
    }

    IWString_and_File_Descriptor stream_for_predictors;

    if (cl.option_present('X')) {
      if (0 != lcol) {
        cerr << "Sorry, the -X option only works when the identifier is in column 1 - "
                "see Ian\n";
        return 1;
      }

      const char* x = cl.option_value('X');

      if (!stream_for_predictors.open(x)) {
        cerr << "Cannot open stream for predictor values '" << x << "'\n";
      }

      if (verbose) {
        cerr << "Predictor values written to '" << x << "'\n";
      }
    }

    const char* l = cl.option_value('L');

    if (!read_desired_row_ordering(l, id_to_ndx, ndx_to_id, lcol,
                                   stream_for_predictors)) {
      cerr << "Cannot read desired row ordering from '" << l << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Read " << id_to_ndx.size() << " row orderings from '" << l << "'\n";
    }
  }

  int icol = 0;

  if (cl.option_present('c')) {
    if (!cl.value('c', icol) || icol < 1) {
      cerr << "Invalid input file column number (-c)\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Identifiers in input file in column " << icol << '\n';
    }
  }

  if (cl.option_present('r')) {
    remove_identifier = 1;

    if (verbose) {
      cerr << "Identifier removed from output file\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!descriptor_file_same_row_order(cl[i], icol, id_to_ndx, ndx_to_id, output)) {
      rc = i + 1;
      break;
    }
  }

  if (nullptr != ndx_to_id) {
    delete[] ndx_to_id;
  }

  output.flush();

  if (verbose) {
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = descriptor_file_same_row_order(argc, argv);

  return rc;
}
