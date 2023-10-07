/*
  Selects rows from a descriptor file
*/

#include <stdlib.h>

#include <fstream>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static int identifier_column = 0;

static int strip_leading_zeros_from_identifiers = 0;

static int records_read = 0;

static int records_written = 0;

static int invert_selection = 0;

static IWString_and_File_Descriptor stream_for_records_not_written;
static IWString_and_File_Descriptor stream_for_identifiers_not_written;

static int allow_multiple_record_dataitems = 0;

static char token_separator = ' ';

static int identifier_column_in_descriptor_file = 0;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Selects rows from a descriptor file: <idfile> <descriptor file>\n";
  cerr << " -c <col>       identifiers are in column <col> of identifier file\n";
  cerr << " -C <col>       identifiers are in column <col> of descriptor file (default 1)\n";
  cerr << " -K <id>        fetch individual ID's\n";
  cerr << " -z             strip leading zero's from identifiers\n";
  cerr << " -Y <fname>     file for descriptor file records not written\n";
  cerr << " -X <fname>     file for identifiers not written\n";
  cerr << " -m             write all instances of identifiers, not just first\n";
  cerr << " -x             invert selection, write records NOT in identifier file\n";
  cerr << " -i <char>      input token separator (space by default)\n";
  cerr << " -t             input is tab separated\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
process_items_in_hash(const IWString_STL_Hash_Set& ids_to_select,
                      IWString_and_File_Descriptor& output)
{
  if (allow_multiple_record_dataitems) {  // ids_to_select will be unchanged
    return 1;
  }

  if (0 == ids_to_select.size()) {
    return 1;
  }

  if (verbose) {
    cerr << ids_to_select.size() << " identifiers not found in descriptor file\n";
  }
  if (stream_for_identifiers_not_written.is_open()) {
    for (IWString_STL_Hash_Set::const_iterator f = ids_to_select.begin();
         f != ids_to_select.end(); ++f) {
      stream_for_identifiers_not_written << (*f) << '\n';
      stream_for_identifiers_not_written.write_if_buffer_holds_more_than(32768);
    }
  } else if (verbose) {
    for (IWString_STL_Hash_Set::const_iterator f = ids_to_select.begin();
         f != ids_to_select.end(); ++f) {
      cerr << (*f) << endl;
    }
  }

  return 1;
}

static int
write_header_if_needed_and_record(const IWString& header, const const_IWSubstring& buffer,
                                  int& header_written,
                                  IWString_and_File_Descriptor& output)
{
  if (!header_written) {
    output << header << '\n';
    header_written = 1;
  }

  output << buffer << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
descriptor_file_select_rows(iwstring_data_source& input,
                            IWString_STL_Hash_Set& ids_to_select,
                            IWString_and_File_Descriptor& output)
{
  IWString header;

  if (!input.next_record(header) || 0 == header.length()) {
    cerr << "Cannot read first record of descriptor file\n";
    return 0;
  }

  int header_written_to_output = 0;
  int header_written_to_dash_B_file = 0;

  records_read++;

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    records_read++;

    const_IWSubstring id;
    if (0 == identifier_column_in_descriptor_file) {
      id = buffer;
      id.truncate_at_first(token_separator);
    } else if (!buffer.word(identifier_column_in_descriptor_file, id, token_separator)) {
      cerr << "Cannot extract identifier from column "
           << (identifier_column_in_descriptor_file + 1) << " in '" << buffer << "'\n";
      return 0;
    }

    if (strip_leading_zeros_from_identifiers) {
      id.remove_leading_chars('0');
    }

    int found_in_hash = ids_to_select.contains(id);
    if (invert_selection) {
      found_in_hash = !found_in_hash;
    }

    if (!found_in_hash) {
      if (stream_for_records_not_written.active()) {
        write_header_if_needed_and_record(header, buffer, header_written_to_dash_B_file,
                                          stream_for_records_not_written);
      }
      continue;
    }

    write_header_if_needed_and_record(header, buffer, header_written_to_output, output);

    records_written++;

    if (!allow_multiple_record_dataitems) {
      ids_to_select.erase(id);
    }
  }

  if (records_written) {
    records_written++;  // we also wrote the header
  }

  return process_items_in_hash(ids_to_select, output);

  if (allow_multiple_record_dataitems) {  // ids_to_select will be unchanged
    ;
  } else if (ids_to_select.size()) {
    cerr << ids_to_select.size() << " identifiers not found in descriptor file\n";
    if (stream_for_identifiers_not_written.is_open()) {
      for (IWString_STL_Hash_Set::const_iterator f = ids_to_select.begin();
           f != ids_to_select.end(); ++f) {
        stream_for_identifiers_not_written << (*f) << '\n';
        stream_for_identifiers_not_written.write_if_buffer_holds_more_than(32768);
      }
    } else if (verbose) {
      for (IWString_STL_Hash_Set::const_iterator f = ids_to_select.begin();
           f != ids_to_select.end(); ++f) {
        cerr << (*f) << endl;
      }
    }
  }

  return output.good();
}

static int
descriptor_file_select_rows(const char* fname, IWString_STL_Hash_Set& ids_to_fetch,
                            IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptor_file_select_rows(input, ids_to_fetch, output);
}

static int
read_ids(const const_IWSubstring& buffer, IWString_STL_Hash_Set& ids_to_fetch)
{
  IWString id;

  if (!buffer.word(identifier_column, id, token_separator)) {
    cerr << "Cannot extract column " << (identifier_column + 1) << " from input\n";
    return 0;
  }

  if (strip_leading_zeros_from_identifiers) {
    id.remove_leading_chars('0');
  }

  ids_to_fetch.insert(id);

  return 1;
}

static int
read_ids(iwstring_data_source& input, IWString_STL_Hash_Set& ids_to_fetch)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!read_ids(buffer, ids_to_fetch)) {
      cerr << "Fatal error on line " << input.lines_read() << " '" << buffer << "'\n";
      return 0;
    }
  }

  if (verbose) {
    cerr << "Read " << ids_to_fetch.size() << " identifiers to select\n";
  }

  return ids_to_fetch.size();
}

static int
read_ids(const char* fname, IWString_STL_Hash_Set& ids_to_fetch)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_ids(input, ids_to_fetch);
}

static int
read_identifiers(const const_IWSubstring& k, IWString_STL_Hash_Set& ids_to_fetch)
{
  int i = 0;
  IWString token;
  while (k.nextword(token, i, ',')) {
    //  cerr << "Fetching '" << token << "'\n";
    ids_to_fetch.insert(token);
  }

  return ids_to_fetch.size();
}

static int
descriptor_file_select_rows(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vzi:tc:X:Y:mxK:C:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('z')) {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose) {
      cerr << "Will strip leading 0's from identifiers\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', identifier_column) || identifier_column < 1) {
      cerr << "The identifier column must be a positive whole number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Identifiers in column " << identifier_column << endl;
    }

    identifier_column--;
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', identifier_column_in_descriptor_file) ||
        identifier_column_in_descriptor_file < 1) {
      cerr << "The identifier column in descriptor file option (-C) must be a positive "
              "whole number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Identifiers in column " << identifier_column
           << " of the descriptor file\n";
    }

    identifier_column_in_descriptor_file--;
  }

  if (cl.option_present('i') && cl.option_present('t')) {
    cerr << "The -i and -t options cannot both be used\n";
    usage(1);
  }

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');

    if (!char_name_to_char(i)) {
      cerr << "Unrecognised token separator '" << i << "'\n";
      return 1;
    }

    token_separator = i[0];

    if (verbose) {
      cerr << "Token separator set to '" << token_separator << "'\n";
    }
  } else if (cl.option_present('t')) {
    token_separator = '\t';
  }

  IWString_STL_Hash_Set ids_to_fetch;

  int file_to_process = 1;

  if (cl.option_present('K')) {
    if (1 != cl.number_elements()) {
      cerr
          << "When specifying items on command line (-K), there must be one input file\n";
      usage(2);
    }

    int i = 0;
    const_IWSubstring k;
    while (cl.value('K', k, i++)) {
      read_identifiers(k, ids_to_fetch);
    }

    file_to_process = 0;
  } else if (2 != cl.number_elements()) {
    cerr << "Insufficient arguments, must have identifier file and descriptor file\n";
    usage(2);
  }

  if (ids_to_fetch.size() > 0) {  // read from the -K option
    ;
  } else if (!read_ids(cl[0], ids_to_fetch) || 0 == ids_to_fetch.size()) {
    cerr << "Cannot read identifiers to be fetched '" << cl[0] << "'\n";
    return 5;
  }

  if ((cl.option_present('Y') || cl.option_present('X')) && cl.option_present('m')) {
    cerr << "Cannot use the -m option with either the -X or -Y options\n";
    return 6;
  }

  if (cl.option_present('Y')) {
    const char* b = cl.option_value('Y');

    if (!stream_for_records_not_written.open(b)) {
      cerr << "Cannot open rejected record file '" << b << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Discarded records written to '" << b << "'\n";
    }
  }

  if (cl.option_present('X')) {
    const char* b = cl.option_value('X');

    if (!stream_for_identifiers_not_written.open(b)) {
      cerr << "Cannot open unwritten identifier file '" << b << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Identifiers not written will be written to '" << b << "'\n";
    }
  }

  if (cl.option_present('m')) {
    allow_multiple_record_dataitems = 1;

    if (verbose) {
      cerr << "Multi record dataitems permitted\n";
    }
  }

  if (cl.option_present('x')) {
    invert_selection = 1;

    if (verbose) {
      cerr << "Will invert selection action\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = descriptor_file_select_rows(cl[file_to_process], ids_to_fetch, output);

  output.flush();

  if (verbose) {
    cerr << "Wrote " << records_written << " of " << records_read << " records read\n";
  }

  if (rc) {
    return 0;
  } else {
    return 1;
  }
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_select_rows(argc, argv);

  return rc;
}
