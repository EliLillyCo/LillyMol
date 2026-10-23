// Fetches identifiers from another file in the order specified
// by the list of identifiers.
// Different from fetch_smiles_quick, that retrieves things based on
// the order in the data file.

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

namespace fetch_smiles {

using std::cerr;

// what to do if duplicate identifiers are encountered.

enum DuplicateIdentifer {
  kFail,
  kTakeFirst,
  kTakeLast
};

struct JobParams {
  // Column where the identifier is found.
  int query_file_column = 0;

  // Default set up for smiles.
  int data_file_column = 1;

  // What to do if  an identifier is not found.
  int ignore_missing_identifiers = 0;

  bool identifier_file_is_descriptor_file = false;

  // What to do if we encounter a duplicate identifier.
  DuplicateIdentifer duplicate_identifier;

  char data_file_token_separator = ' ';
  char query_file_token_separator = ' ';
};

// The number of identifiers not found.
int missing_identifier = 0;

void
Usage(int rc) {
  cerr << R"(Fetches records from one file based on the identifiers in another.
The key to this tool is that the records fetched from the donor file will be in the same order that
they appear in the identifier file. It works by reading the whole file, determining where each identifier
is, and then reading the file of identifiers to be written one at a time. It then writes the requested
identifier line.
This is fundamentally different from fetch_smiles_quick which ingests the identifiers, and then as
records are read, a record is written if it has been requested. That is fundamentally more efficient
and should be used unless there is a specific reason to preserve the ordering.
 -c <col>       identifier column in identifier file.
 -C <col>       identifier column in smiles file.
 -j             skip the identifier header record - useful for descriptor files.
 -x             ignore missing identifiers.
 -v             verbose output.
)";

  exit(rc);
}

int
FetchSmilesRecord(const const_IWSubstring& line,
            const JobParams& params,
            IW_STL_Hash_Map_String& data,
            IWString_and_File_Descriptor& output) { 
  const_IWSubstring token;
  int i = 0;
  for (int col = 0;
       line.nextword_single_delimiter(token, i, params.data_file_token_separator);
       ++col) {
    if (col != params.query_file_column) {
      continue;
    }
    const auto iter = data.find(token);
    if (iter == data.end()) {
      if (params.ignore_missing_identifiers) {
        missing_identifier++;
        continue;
      }
      cerr << "No data for " << token << '\n';
      return 0;
    }
    output << iter->second << '\n';
    output.write_if_buffer_holds_more_than(4096);
    return 1;
  }

  cerr << "Did not find column " << params.data_file_token_separator << " in " << line << '\n';
  return 0;
}

int
FetchSmiles(iwstring_data_source& input,
            const JobParams& params,
            IW_STL_Hash_Map_String& data,
            IWString_and_File_Descriptor& output) {
  const_IWSubstring line;
  if (params.identifier_file_is_descriptor_file) {
    input.next_record(line);
  }

  while (input.next_record(line)) {
    if (! FetchSmilesRecord(line, params, data, output)) {
      cerr << "Error processing '" << line << "'\n";
      return 0;
    }
  }

  return 1;
}

int
FetchSmiles(const char * fname,
            const JobParams& params,
            IW_STL_Hash_Map_String& data,
            IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.ok()) {
    cerr << "Cannot open identifier file '" << fname << "'\n";
    return 0;
  }

  return FetchSmiles(input, params, data, output);
}

int
ReadRecord(const const_IWSubstring& line,
           const JobParams& params,
           IW_STL_Hash_Map_String& data) {
  const_IWSubstring token;
  int i = 0;
  for (int col = 0; line.nextword_single_delimiter(token, i, params.query_file_token_separator); ++col) {
    if (col == params.data_file_column) {
      const auto iter = data.find(token);
      if (iter == data.end()) {
        data.emplace(IWString(token), IWString(line));
        return 1;
      }
      switch (params.duplicate_identifier) {
        case kFail:
          cerr << "Duplicate identifier '" << token << "'\n";
          return 0;
        case kTakeFirst:
          return 1;
        case kTakeLast:
          data.emplace(IWString(token), IWString(line));
          return 1;
      }
    }
  }

  cerr << "Did not find column " << params.data_file_column << " in '" << line << "'\n";
  return 0;
}

int
ReadData(iwstring_data_source& input,
         const JobParams& params,
         IW_STL_Hash_Map_String& data) {
  const_IWSubstring line;

  while (input.next_record(line)) {
    if (! ReadRecord(line, params, data)) {
      cerr << "Unable to process " << line << '\n';
      return 0;
    }
  }

  return data.size();
}

int
ReadData(const char * fname,
         const JobParams& params,
         IW_STL_Hash_Map_String& data) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadData(input, params, data);
}

int
FetchSmiles(int argc, char** argv) {
  Command_Line cl(argc, argv, "vc:C:xj");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  JobParams params;

  params.data_file_column = 1;
  params.query_file_column = 0;

  if (cl.option_present('c')) {
    if (! cl.value('c', params.query_file_column) || params.query_file_column < 1) {
      cerr << "Invalid query file column\n";
      Usage(1);
    }
    if (verbose)
      cerr << "Query file identifiers in column " << params.query_file_column << '\n';
    params.query_file_column--;
  }

  if (cl.option_present('C')) {
    if (! cl.value('c', params.data_file_column) || params.data_file_column < 1) {
      cerr << "Invalid data file column\n";
      Usage(1);
    }
    if (verbose)
      cerr << "data file identifiers in column " << params.data_file_column << '\n';
    params.data_file_column--;
  }

  if (cl.option_present('x')) {
    params.ignore_missing_identifiers = 1;
    if (verbose) {
      cerr << "Missing identifiers will be ignored\n";
    }
  }

  if (cl.option_present('j')) {
    params.identifier_file_is_descriptor_file = true;
    if (verbose) {
      cerr << "Identifier file is a descriptor file - first record skipped\n";
    }
  }

  if (cl.size() < 2) {
    cerr << "Must specify at least 2 files\n";
    Usage(1);
  }

  IW_STL_Hash_Map_String data;
  for (int i = 1; i < cl.number_elements(); ++i) {
    if (! ReadData(cl[i], params, data)) {
      cerr << "Cannot read data from '" << cl[i] << "'\n";
      return 1;
    }
  }

  if (verbose)
    cerr << "Read " << data.size() << " items\n";

  IWString_and_File_Descriptor output(1);
  if (! FetchSmiles(cl[0], params, data, output)) {
    return 1;
  }

  return 0;
}

}  // namespace fetch_smiles


int
main (int argc, char ** argv)
{
  int rc = fetch_smiles::FetchSmiles(argc, argv);

  return rc;
}
