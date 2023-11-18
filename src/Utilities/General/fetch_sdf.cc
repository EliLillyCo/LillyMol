/*
  Fetches sd file records based on identifiers in another file
*/

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static int identifier_column_in_identifier_file = 0;

static int identifier_column_in_sdf_name = -1;

static IWString identifier_tag;  // the angle brackets only

static IWString identifier_tag_in_sd_file;  // with the > ..*< stuff

static int ignore_duplicate_identifiers_in_identifier_file = 0;

static int ignore_duplicate_identifiers_in_sdfile = 0;

static int duplicate_identifiers_in_identifier_file = 0;

static IWString_and_File_Descriptor stream_for_not_in_identifier_file;
static IWString_and_File_Descriptor stream_for_not_in_sd_file;

/*
  RDF files are all over the place, so we ask the user to specify
  a regexp for the start of each entry
*/
static std::unique_ptr<re2::RE2> rdf_start_of_entry_rx;

static int swap_identifier_and_data = 0;

static int sdfile_records_read = 0;
static int items_written = 0;
static int not_in_identifier_file = 0;

static int erase_identifiers_written = 1;

static int include_identifier_file_info = 1;

static IWString string_before_identifier_file_info(' ');

static int strip_leading_zeros_from_identifiers = 0;

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
  cerr << "Fetches records from one file based on identifiers in another file\n";
  cerr << prog_name << " identifier_file sd_file > newfile\n";
  cerr << " -c <col>       identifier column in identifier file\n";
  cerr << " -C <tag>       identifier tag in sd file - default is 1st line\n";
  cerr << " -f <col>       identifier column in name field of .sdf file (default whole record)\n";
  cerr << " -d             ignore duplicate identifiers in identifier file\n";
  cerr << " -a             write all instances of identifiers in sd file\n";
  cerr << "                by default, only the first is written\n";
  cerr << " -X <fname>     write sd records not in <identifier_file> to <fname>\n";
  cerr << " -Y <fname>     write identifiers not in <sd_file> to <fname>\n";
  cerr << " -k             suppress addition of info from identifier file\n";
  cerr << " -n <string>    string to insert between record and info from identifier file\n";
  cerr << " -w             swap id and info in the -Y file (useful for smiles)\n";
  cerr << " -R <regexp>    input file is an RDF file rather than SDF\n";
  cerr << "                <regexp> is a regular expression to define start of record\n";
  cerr << " -z             strip leading zero's from identifiers\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
separate_into_id_and_data_2(const const_IWSubstring& buffer, int col, IWString& id,
                            IWString& zdata)
{
  if (0 == col)  // the most common case
  {
    //  cerr << "Extracting column " << col << " from '" << buffer << "'\n";
    int ispace = buffer.index(' ');

    if (ispace < 0)  // just one token on the line
    {
      id = buffer;
      return 1;
    }

    const char* s = buffer.rawchars();
    id.set(s, ispace);

    zdata.set(s + ispace + 1, buffer.length() - ispace - 1);

    //  cerr << "Split into '" << id << "' and '" << zdata << "'\n";

    return 1;
  }

  if (!buffer.word(col, id)) {
    cerr << "Cannot extract column " << col << " from '" << buffer << "'\n";
    return 0;
  }

  zdata = buffer;
  zdata.remove_word(col);

  return 1;
}

static int
separate_into_id_and_data(const const_IWSubstring& buffer, int col, IWString& id,
                          IWString& zdata)
{
  if (!separate_into_id_and_data_2(buffer, col, id, zdata)) {
    return 0;
  }

  if (strip_leading_zeros_from_identifiers) {
    id.remove_leading_chars('0');
  }

  return 1;
}

#ifdef NOT_BEING_USED
static int
insert_identifier_file_info(IWString& sd, const IWString& zinfo)
{
  if (0 == zinfo.length()) {
    return 1;
  }

  cerr << "Insert identifier info feature not implemented, see Ian\n";
  return 0;
}
#endif

/*static int
read_next_rdf (iwstring_data_source & input,
              IWString & id)
{
  assert (rdf_start_of_entry_rx);
  id.resize_keep_storage(0);

  const_IWSubstring buffer;

  if (! input.next_record (buffer))   // normal eof
    return 0;

  if (! iwre2::RE2FullMatch(buffer, *rdf_start_of_entry_rx))
  {
    cerr << "Bad news, first record in RDF does not match pattern '" <<
rdf_start_of_entry_rx->pattern() << "'\n"; return 0;
  }

  int next_record_is_id = 0;
  int lines_read = 1;

  while (input.next_record (buffer))
  {
    lines_read++;

    if (iwre2::RE2FullMatch(buffer, *rdf_start_of_entry_rx))
    {
      input.push_record ();
      return 1;
    }

    if (buffer == identifier_tag_in_sd_file)
    {
      if (id.length ())
        cerr << "yipes, duplicate " << identifier_tag_in_sd_file << "' in file\n";
      next_record_is_id = 1;
    }
    else if (next_record_is_id)
    {
      id = buffer;
      if (! buffer.starts_with ("$DATUM"))
      {
        cerr << "Should be '$DATUM', but got '" << buffer << "'\n";
        return 0;
      }
      id.remove_leading_words (1);

      next_record_is_id = 0;
    }
  }

  return 1;
}*/

static int
looks_like_sdf_id(const const_IWSubstring& buffer, const IWString& tag)
{
  if (!buffer.starts_with("> ")) {  // must be at least that much present
    return 0;
  }

  if (buffer.ends_with(tag)) {  // good enough
    return 1;
  }

  if (!buffer.contains(tag)) {
    return 0;
  }

  // Maybe try some more things some time later

  return 0;
}

static int
read_next_sd(iwstring_data_source& input, IWString& id)
{
  id.resize_keep_storage(0);

  const_IWSubstring buffer;

  int got_dollars = 0;
  int next_record_is_id = 0;
  int lines_read = 0;

  while (input.next_record(buffer)) {
    lines_read++;

    if (buffer.starts_with("$$$$")) {
      got_dollars = 1;
      break;
    }

    if (0 == identifier_tag_in_sd_file.length()) {
      if (1 == lines_read) {
        id = buffer;
      }
    } else if (buffer.starts_with(identifier_tag_in_sd_file) ||
               looks_like_sdf_id(buffer, identifier_tag)) {
      if (id.length() > 0) {
        cerr << "yipes, duplicate " << identifier_tag_in_sd_file << "' in file, id '"
             << id << "', line " << input.lines_read() << endl;
      }
      next_record_is_id = 1;
    } else if (next_record_is_id) {
      id = buffer;
      next_record_is_id = 0;
    }
  }

  if (0 == lines_read) {  // normal eof
    return 0;
  }

  if (!got_dollars) {
    cerr << "Premature EOF\n";
    return 0;
  }

  if (identifier_column_in_sdf_name < 0) {
    return 1;
  }

  if (0 == identifier_column_in_sdf_name) {
    id.truncate_at_first(' ');
    return 1;
  }

  if (id.nwords() <= identifier_column_in_sdf_name) {
    cerr << "Identifier column " << (identifier_column_in_sdf_name + 1) << " but id '"
         << id << "' does not have enough tokens\n";
    return 0;
  }

  id.remove_leading_words(identifier_column_in_sdf_name);

  return 1;
}

/*static int
read_next_entity (iwstring_data_source & input,
                  IWString & id)
{
  id.resize_keep_storage (0);

  int rc;
  if (rdf_start_of_entry_rx)
    rc = read_next_rdf(input, id);
  else
    rc = read_next_sd(input, id);

  if (0 == rc)
    return 0;

  if (0 == id.length ())
  {
    cerr << "No identifier '" << identifier_tag_in_sd_file << "'\n";
    return 0;
  }

  return rc;
}*/

static int
determine_identifier_offsets(iwstring_data_source& input,
                             IW_STL_Hash_Map_long& id_to_offset)
{
  off_t offset = static_cast<off_t>(0);

  if (rdf_start_of_entry_rx) {
    if (!input.skip_records(*rdf_start_of_entry_rx, 1)) {
      cerr << "Cannot find first rdf entry\n";
      return 0;
    }

    input.push_record();

    offset = input.tellg();
  }

  // cerr << "File starts at " << offset << endl;

  IWString id;
  while (read_next_sd(input, id)) {
    if (id_to_offset.contains(id)) {
      cerr << "Duplicate identifier in SD file '" << id << "'\n";
      if (!ignore_duplicate_identifiers_in_sdfile) {
        return 0;
      }
    }
    id_to_offset[id] = offset;
    offset = input.tellg();
  }

  return id_to_offset.size();
}

static int
echo_sdfile(iwstring_data_source& sdfile, IWString& output_buffer)
{
  const_IWSubstring buffer;
  while (sdfile.next_record(buffer)) {
    output_buffer << buffer << '\n';
    if ("$$$$" == buffer) {
      return 1;
    }
  }

  cerr << "EOF without encountering $$$$\n";
  return 0;
}

static int
fetch_sdf(const const_IWSubstring& buffer, IW_STL_Hash_Set& identifier_encountered,
          iwstring_data_source& sdfile, IW_STL_Hash_Map_long& id_to_offset,
          IWString& output_buffer)
{
  IWString id, zdata;

  if (!separate_into_id_and_data(buffer, identifier_column_in_identifier_file, id,
                                 zdata)) {
    return 0;
  }

  if (!identifier_encountered.contains(id)) {
    identifier_encountered.insert(id);
  } else if (ignore_duplicate_identifiers_in_identifier_file) {
    if (verbose) {
      cerr << "Ignoring duplicate identifer '" << id << "'\n";
    }
    duplicate_identifiers_in_identifier_file++;
  } else {
    cerr << "Duplicate identifier in identifier file '" << id << "'\n";
    return 0;
  }

  IW_STL_Hash_Map_long::const_iterator f = id_to_offset.find(id);

  if (f == id_to_offset.end())  // Entry not in sd file
  {
    if (verbose > 1) {
      cerr << "Identifier '" << id << "' not found in sd file\n";
    }

    if (stream_for_not_in_sd_file.is_open()) {
      if (swap_identifier_and_data && zdata.length()) {
        stream_for_not_in_sd_file << zdata << ' ' << id << '\n';
      } else {
        stream_for_not_in_sd_file << id;
        if (zdata.length()) {
          stream_for_not_in_sd_file << ' ' << zdata;
        }
        stream_for_not_in_sd_file << '\n';
      }

      stream_for_not_in_sd_file.write_if_buffer_holds_more_than(32768);
    }

    return 1;
  }

  if (!sdfile.seekg((*f).second)) {
    cerr << "Cannot seek to " << (*f).second << " for '" << id << "'\n";
    return 0;
  }

  if (erase_identifiers_written) {
    id_to_offset.erase(id);
  }

  items_written++;

  return echo_sdfile(sdfile, output_buffer);
}

static int
fetch_sdf(iwstring_data_source& identifier_file, iwstring_data_source& sdfile,
          IW_STL_Hash_Map_long& id_to_offset, IWString_and_File_Descriptor& output)
{
  IW_STL_Hash_Set identifier_encountered;

  identifier_file.set_dos(1);

  const_IWSubstring buffer;
  while (identifier_file.next_record(buffer)) {
    sdfile_records_read++;

    if (!fetch_sdf(buffer, identifier_encountered, sdfile, id_to_offset, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
fetch_sdf(const char* identifier_fname, iwstring_data_source& sdfile,
          IW_STL_Hash_Map_long& id_to_offset, IWString_and_File_Descriptor& output)
{
  iwstring_data_source identifier_stream(identifier_fname);
  if (!identifier_stream.good()) {
    cerr << "Cannot open identifier file '" << identifier_fname << "'\n";
    return 0;
  }

  identifier_stream.set_strip_trailing_blanks(1);

  return fetch_sdf(identifier_stream, sdfile, id_to_offset, output);
}

static int
fetch_sdf(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vc:C:X:Y:kdzn:wzaR:f:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('n')) {
    if (cl.option_present('k')) {
      cerr << "The -n and -k options are mutually exclusive\n";
      usage(5);
    }

    cl.value('n', string_before_identifier_file_info);

    if (verbose) {
      cerr << "Will put '" << string_before_identifier_file_info
           << " before identifier file info\n";
    }
  }

  if (cl.option_present('a')) {
    if (cl.option_present('Y')) {
      cerr << "Sorry, the -a and -Y options do not work together\n";
      return 6;
    }

    erase_identifiers_written = 0;

    if (verbose) {
      cerr << "Will write all instances of identifiers in the smiles file\n";
    }
  }

  if (cl.option_present('d')) {
    ignore_duplicate_identifiers_in_identifier_file = 1;

    if (verbose) {
      cerr << "Will ignore duplicates in identifier file\n";
    }
  }

  if (cl.option_present('z')) {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose) {
      cerr << "Will strip leading zero's from identifiers\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', identifier_column_in_identifier_file) ||
        identifier_column_in_identifier_file < 1) {
      cerr << "The column in identifier file option (-c) must be a whole +ve number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Identifiers in identifier file in column "
           << identifier_column_in_identifier_file << endl;
    }

    identifier_column_in_identifier_file--;
  }

  if (cl.option_present('R')) {
    const_IWSubstring r = cl.string_value('R');

    if (!iwre2::RE2Reset(rdf_start_of_entry_rx, r)) {
      cerr << "Invalid RDF regular expression '" << r << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Will treat as an RDF file\n";
    }

    if (!cl.option_present('C')) {
      cerr << "When processing an rdf file, the -C option is mandatory\n";
      usage(4);
    }
  }

  if (cl.option_present('C')) {
    const_IWSubstring c;

    cl.value('C', c);

    identifier_tag << '<' << c << '>';

    if (rdf_start_of_entry_rx) {
      identifier_tag_in_sd_file << "$DTYPE " << c;
    } else {
      identifier_tag_in_sd_file << ">  <" << c << ">";
    }

    if (verbose) {
      cerr << "Identifiers in sd file in '" << identifier_tag_in_sd_file << "' tag\n";
    }
  }

  if (cl.option_present('f')) {
    if (!cl.value('f', identifier_column_in_sdf_name) ||
        identifier_column_in_sdf_name < 1) {
      cerr << "The identifier column in .sdf file name (-f) must be a valid column "
              "number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Identifier column in .sdf file in column " << identifier_column_in_sdf_name
           << endl;
    }

    identifier_column_in_sdf_name--;
  }

  if (cl.number_elements() < 2) {
    cerr << "Must specify at least two files\n";
    usage(2);
  }

  if (cl.option_present('k')) {
    include_identifier_file_info = 0;

    if (verbose) {
      cerr << "Suppress addition of identifier file info\n";
    }
  }

  if (cl.option_present('X')) {
    const char* fname = cl.option_value('X');

    if (!stream_for_not_in_identifier_file.open(fname)) {
      cerr << "Cannot open stream for not in identifier file (-X option), '" << fname
           << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Smiles not in identifier file written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('Y')) {
    const char* fname = cl.option_value('Y');

    if (!stream_for_not_in_sd_file.open(fname)) {
      cerr << "Cannot open stream for not in smiles file (-Y option), '" << fname
           << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "Identifiers not in smiles file written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('w')) {
    if (!cl.option_present('Y')) {
      cerr << "The -w option only makes sense with the -Y option\n";
      usage(3);
    }

    swap_identifier_and_data = 1;

    if (verbose) {
      cerr << "Will swap positions of identifier and data in -Y file\n";
    }
  }

  if (cl.number_elements() < 2) {
    cerr << "Must specify both identifier file and sd file\n";
    usage(1);
  } else if (cl.number_elements() > 2) {
    cerr << "Sorry, can only do one sd file at a time\n";
    usage(2);
  }

  iwstring_data_source sdfile(cl[1]);
  if (!sdfile.good()) {
    cerr << "Cannot open input sd file '" << cl[1] << "'\n";
    return 5;
  }

  sdfile.set_dos(1);
  sdfile.set_strip_trailing_blanks(1);

  IW_STL_Hash_Map_long id_to_offset;

  if (!determine_identifier_offsets(sdfile, id_to_offset)) {
    cerr << "Cannot determine offsets\n";
    return 3;
  }

  if (verbose) {
    cerr << "SDfile contains " << id_to_offset.size() << " items\n";
  }

  if (0 == id_to_offset.size()) {
    cerr << "No structures in sdfile '" << cl[1] << "'\n";
    return 5;
  }

  IWString_and_File_Descriptor output(1);

  if (!fetch_sdf(cl[0], sdfile, id_to_offset, output)) {
    return 3;
  }

  if (0 == items_written) {
    cerr << "Warning, nothing written!\n";
  }

  if (verbose) {
    cerr << sdfile_records_read << " sd file records read, " << items_written
         << " items written\n";
    if (not_in_identifier_file) {
      cerr << not_in_identifier_file << " items not in identifier file\n";
    }
  }

  output.close();

  if (stream_for_not_in_identifier_file.is_open()) {
    for (IW_STL_Hash_Map_long::const_iterator i = id_to_offset.begin();
         i != id_to_offset.end(); i++) {
      off_t o = (*i).second;

      if (!sdfile.seekg(o)) {
        cerr << "Cannot seek to " << o << " for '" << (*i).first << "'\n";
        return 4;
      }

      echo_sdfile(sdfile, stream_for_not_in_identifier_file);

      stream_for_not_in_identifier_file.write_if_buffer_holds_more_than(32768);
    }

    stream_for_not_in_identifier_file.close();
  }

  if (stream_for_not_in_sd_file.is_open()) {
    stream_for_not_in_sd_file.close();
  }

  if (cl.option_present('q')) {
    _exit(0);
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = fetch_sdf(argc, argv);

  return rc;
}
