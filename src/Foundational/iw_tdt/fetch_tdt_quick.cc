/*
  Fetches specific TDT's from a file of TDT's
*/

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "cmdline.h"
#include "iw_stl_hash_map.h"
#include "misc.h"

#include "iw_tdt.h"
#include "fetch_via_regexp.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString identifier_tag ("PCN<");

static int identifier_column_in_identifier_file = -1;

static int identifier_column_in_TDT = -1;

static int ignore_duplictes_in_identifier_file = 0;

static int ignore_tdts_missing_the_tag = 0;

static int tdts_read = 0;

static int tdts_written = 0;

static int identifiers_written = 0;

static int write_all_instances_of_duplicates = 0;

static int which_one_to_retrieve = 0;

static IWString_and_File_Descriptor stream_for_tdts_not_in_identifier_file;

static int write_buffer_size = 32768;

static int invert_fetching_operation = 0;

/*
  Feb 2004.
  Allow replacing the data in the TDT with what is in the identifier file
*/

static int replace_tdt_data_with_data_from_identifier_file = 0;

static IW_STL_Hash_Map<IWString, IWString> replacement_info;

static int strip_leading_zeros = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Fetches a subset of TDT's from a file of TDT's\n";
  cerr << prog_name << " options <identifier_file> <tdt_file>\n";
  cerr << " -T <tag>       tag for identifiers\n";
  cerr << " -w <number>    which instance of <tag> to process (default 1)\n";
  cerr << " -p <rx>        identifiers must match regular expression <rx>\n";
  cerr << " -i             ignore identifiers not matching <rx>\n";
  cerr << " -d ignore      ignore duplicates in <identifier_file>\n";
  cerr << " -D all         fetch all instances of an identifier in <tdt_file>\n";
  cerr << " -c <col>       identifier in column <col> in the identifier file\n";
  cerr << " -C <col>       identifier in column <col> in the tdt file\n";
  cerr << " -X <fname>     write tdts not in <identifier file> to <fname>\n";
  cerr << " -Y <fname>     write identifiers not in <tdt_file> to <fname>\n";
  cerr << " -I <id>        identifiers to fetch\n";
  cerr << " -R <rx>        identifiers to fetch - specified by regular expressions\n";
  cerr << " -r             replace identifier in TDT with data from identifier file\n";
//cerr << " -b <bytes>     write buffer size\n";
  cerr << " -x             invert behaviour, selection becomes deselection\n";
  cerr << " -z             strip leading 0's from identifiers\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
extract_identifier (const IWString & buffer,
                    int col,
                    IWString & id,
                    int & fatal)
{
  if (col >= 0)
  {
    if (buffer.nwords () < col + 1)
    {
      cerr << "Sorry, cannot extract column " << (col + 1) << " from identifier record\n";
      fatal = 1;
      return 0;
    }

    (void) buffer.word (col, id);
  }
  else
    id = buffer;

  if (strip_leading_zeros)
    id.remove_leading_chars('0');

  return 1;
}

static int
handle_tdts_not_in_identifier_file (const IW_TDT & tdt)
{
  if (! stream_for_tdts_not_in_identifier_file.is_open())
    return 1;

  stream_for_tdts_not_in_identifier_file << tdt;

  stream_for_tdts_not_in_identifier_file.write_if_buffer_holds_more_than(write_buffer_size);

  return stream_for_tdts_not_in_identifier_file.good ();
}

static int
write_while_replacing_the_identifier (const IW_TDT & tdt,
                                      const const_IWSubstring & id,
                                      IWString & output)
{
  int i = 0;
  const_IWSubstring buffer;

  int replacement_info_written = 0;

  while (tdt.next_dataitem (buffer, i))
  {
    if (replacement_info_written || ! buffer.starts_with (identifier_tag))
    {
      output << buffer;
      continue;
    }

    IW_STL_Hash_Map<IWString, IWString>::const_iterator f = replacement_info.find (id);

    assert (f != replacement_info.end ());

    output << identifier_tag << (*f).second << ">\n";

    replacement_info_written = 1;
  }

  output << "|\n";

  return 1;
}

#ifdef FETCH_VIA_IDENTIFIERS

Maybe this is a good idea, maybe not...

class Fetch_via_Identifiers
{
  private:
    IW_STL_Hash_Map_int _identifiers_to_fetch;
    IW_STL_Hash_Map_String _replacement_info;

  public:
    int build (const char * fname, const int col, const int ignore_duplictes_in_identifier_file);
    int build (iwstring_data_source &, const int col, const int ignore_duplictes_in_identifier_file);

    int matches (const IWString & s);
};

int
Fetch_via_Identifiers::build (const char * fname,
                              const int col,
                              const int ignore_duplictes_in_identifier_file)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Fetch_via_Identifiers::cannot open '" << fname << "'\n";
    return 0;
  }

  return build (input, col, ignore_duplictes_in_identifier_file);
}

int
Fetch_via_Identifiers::build (iwstring_data_source & input,
                              const int col,
                              const int ignore_duplictes_in_identifier_file) 
{
  IWString buffer;

  int fatal;
  while (input.next_record(buffer))
  {
    IWString id;
    if (! extract_identifier(buffer, col, id, fatal))
    {
      cerr << "Fetch_via_Identifiers::build:cannot extract identifier from '" << buffer << "'\n";
      return 0;
    }

    if (_identifiers_to_fetch.contains(id))
    {
      if (ignore_duplictes_in_identifier_file)
        continue;

      cerr << "Duplicate identifier in identifier file '" << id << "'\n";
      return 0;
    }

    _identifiers_to_fetch[id] = 0;           // the number of times fetched
    if (replace_tdt_data_with_data_from_identifier_file)
      _replacement_info[id] = buffer;
  }

  return _identifiers_to_fetch.size();
}

int
Fetch_via_Identifiers::matches (const IWString & s)
{
  return _identifiers_to_fetch.contains(s);
}
#endif

static Fetch_via_Regexp fetch_via_regexp;

static int
fetch_tdt_grep (const IW_TDT & tdt,
                const IWString & id,
                IWString_and_File_Descriptor & output)
{
  bool is_match = fetch_via_regexp.matches(id);

  if (invert_fetching_operation)
    is_match = ! is_match;

  if (! is_match)
    return handle_tdts_not_in_identifier_file(tdt);

  output << tdt;

  identifiers_written++;
  tdts_written++;

  return 1;
}

int
fetch_tdt (const IW_TDT & tdt,
           IW_STL_Hash_Map_int & identifiers_to_fetch,
           IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (tdt.dataitem_value (identifier_tag, buffer, which_one_to_retrieve))
    ;
  else if (ignore_tdts_missing_the_tag)
  {
    cerr << "Ignoring TDT without value for '" << identifier_tag << "'\n";
    return 1;
  }
  else
  {
    cerr << "Cannot retrieve '" << identifier_tag << "' from TDT\n";
    return 0;
  }

  IWString id;
  int fatal;

  if (! extract_identifier (buffer, identifier_column_in_TDT, id, fatal))
  {
    cerr << "No identifier\n";
    if (! fatal)
      return 1;

    cerr << "Missing or invalid identifier in '" << buffer << "'\n";
    return 0;
  }

  if (identifiers_to_fetch.empty())
    return fetch_tdt_grep (tdt, id, output);

//cerr << "Checking id '" << id << "'\n";

  const auto f = identifiers_to_fetch.find(id);

  bool is_match;

  if (f == identifiers_to_fetch.end())
    is_match = false;
  else
    is_match = true;

  if (invert_fetching_operation)
    is_match = ! is_match;

//cerr << "Do we want to retrieve '" << id << "'\n";

  if (! is_match)   // not interested in this one
    return handle_tdts_not_in_identifier_file (tdt);

  if (f == identifiers_to_fetch.end())    // we are writing a non-match
    ;
  else if (0 == (*f).second)    // first time we've encountered this ID
    identifiers_written++;
  else if (write_all_instances_of_duplicates)
    ;
  else
  {
    cerr << "Ignoring duplicate '" << id << "' in TDT file\n";
    return 1;
  }

  if (! invert_fetching_operation)
    (*f).second = (*f).second + 1;

  if (replace_tdt_data_with_data_from_identifier_file)
    write_while_replacing_the_identifier (tdt, id, output);
  else
    output << tdt;

  tdts_written++;

  return 1;
}

static int
fetch_tdt (iwstring_data_source & input,
           IW_STL_Hash_Map_int & identifiers_to_fetch,
           IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;
  while (tdt.next (input))
  {
    tdts_read++;

    if (! fetch_tdt (tdt, identifiers_to_fetch, output))
    {
      cerr << "Fatal error\n";
      cerr << tdt;
      return 0;
    }

    output.write_if_buffer_holds_more_than(write_buffer_size);
  }

  return 1;
}

static int
fetch_tdt (const char * fname,
           IW_STL_Hash_Map_int & identifiers_to_fetch,
           IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return fetch_tdt (input, identifiers_to_fetch, output);
}

static int
read_identifiers_to_fetch (iwstring_data_source & input,
                           IW_STL_Hash_Map_int & identifiers_to_fetch)
{
  IWString buffer;

  while (input.next_record (buffer))
  {
    buffer.strip_trailing_blanks ();

    if (0 == buffer.length ())
      continue;

    IWString id;
    int fatal;

    if (! extract_identifier (buffer, identifier_column_in_identifier_file, id, fatal))
    {
      if (! fatal)
        continue;

      cerr << "Missing or invalid identifier in '" << buffer << "'\n";
      return 0;
    }

    if (! identifiers_to_fetch.contains (id))
    {
      identifiers_to_fetch[id] = 0;    // set to the number fetched
      if (replace_tdt_data_with_data_from_identifier_file)
        replacement_info[id] = buffer;
    }
    else if (ignore_duplictes_in_identifier_file)
      cerr << "Ignoring duplicate identifier '" << id << "'\n";
    else
    {
      cerr << "Duplicate identifier in ID file '" << id << "'\n";
      return 0;
    }
  }

  return identifiers_to_fetch.size ();
}

static int
read_identifiers_to_fetch (const char * fname,
                           IW_STL_Hash_Map_int & identifiers_to_fetch)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_identifiers_to_fetch (input, identifiers_to_fetch);
}

static int
report_occurrences_of_identifiers (const extending_resizable_array<int> & occurrences_of_identifier,
                                   ostream & output)
{
  for (int i = 0; i < occurrences_of_identifier.number_elements (); i++)
  {
    if (0 == i || occurrences_of_identifier[i])
      output << occurrences_of_identifier[i] << " identifiers found " << i << " times in TDT\n";
  }

  return output.good ();
}

static int
fetch_tdt (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vc:C:X:Y:p:I:R:id:D:T:vrw:qb:zx");    // -p and -i options no longer do anything

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('c'))
  {
    if (cl.option_present('I') && cl.number_elements() < 2)
    {
      cerr << "The -c option (identifier column for id file) does not make sense without an identifier file\n";
      return 1;
    }

    if (! cl.value ('c', identifier_column_in_identifier_file) || identifier_column_in_identifier_file < 1)
    {
      cerr << "The identifier column in identifier file (-c) must be a valid column number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Identifiers found in column " << identifier_column_in_identifier_file << " of identifier file\n";

    identifier_column_in_identifier_file--;
  }

  if (cl.option_present ('C'))
  {
    if (! cl.value ('C', identifier_column_in_TDT) || identifier_column_in_TDT < 1)
    {
      cerr << "The identifier column in TDT file (-C) must be a valid column number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Identifiers found in column " << identifier_column_in_TDT << " of TDT file\n";

    identifier_column_in_TDT--;
  }

  if (cl.option_present ('T'))
  {
    identifier_tag = cl.string_value ('T');

    if (verbose)
      cerr << "Will look for identifiers in the '" << identifier_tag << "' tag\n";

    if (! identifier_tag.ends_with ('<'))
      identifier_tag += '<';
  }

  if (cl.option_present ('w'))
  {
    if (! cl.value ('w', which_one_to_retrieve) || which_one_to_retrieve < 1)
    {
      cerr << "The which tag item to sort (-w option) must be a whole +ve number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will sort on the " << which_one_to_retrieve << " instance of '" << identifier_tag << "'\n";

    which_one_to_retrieve--;
  }

  if (cl.option_present ('d'))
  {
    int i = 0;
    const_IWSubstring d;
    while (cl.value ('d', d, i++))
    {
      if (d.starts_with ("ign"))
      {
        ignore_duplictes_in_identifier_file = 1;

        if (verbose)
          cerr << "Will ignore duplicates in the identifier file\n";
      }
      else
      {
        cerr << "Unrecognsed -d qualifier '" << d << "'\n";
        usage (5);
      }
    }
  }

  if (cl.option_present ('D'))
  {
    int i = 0;
    const_IWSubstring d;
    while (cl.value ('D', d, i++))
    {
      if ("all" == d)
      {
        write_all_instances_of_duplicates = 1;

        if (verbose)
          cerr << "Will write all instances of duplicate identifiers in TDT file\n";
      }
      else
      {
        cerr << "Unrecognsed -d qualifier '" << d << "'\n";
        usage (5);
      }
    }
  }

  if (cl.option_present ('r'))
  {
    replace_tdt_data_with_data_from_identifier_file = 1;

    if (verbose)
      cerr << "Will replace identifier in TDT with data from identifier file\n";
  }

  if (cl.option_present('b'))
  {
    if (! cl.value('b', write_buffer_size) || write_buffer_size < 2)
    {
      cerr << "The write buffer size must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will flush write buffers when the have " << write_buffer_size << " bytes or more\n";
  }

  if (cl.option_present('z'))
  {
    strip_leading_zeros = 1;

    if (verbose)
      cerr << "Will strip leading 0's from identifiers\n";
  }

  if (cl.option_present('x'))
  {
    invert_fetching_operation = 1;

    if (verbose)
      cerr << "Will discard rather than fetch identifiers\n";
  }

  if (cl.option_present('p') || cl.option_present('i'))
    cerr << "Options -p and -i are deprecated, remove from invocation. Ignored for now\n";

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IW_STL_Hash_Map_int identifiers_to_fetch;

// We need to know where in the command line to start

  int clstart = -1;

  if (cl.option_present ('I'))
  {
    int i = 0;
    const_IWSubstring ii;
    while (cl.value ('I', ii, i++))
    {
      identifiers_to_fetch[ii] = 0;         // number times fetched
    }

    clstart = 0;
  }
  else if (cl.option_present('R'))
  {
    if (! fetch_via_regexp.build(cl, 'R'))
    {
      cerr << "Cannot initialise regular expressions to fetch\n";
      return 1;
    }
    clstart = 0;
  }
  else if (2 == cl.number_elements ())
  {
    if (! read_identifiers_to_fetch (cl[0], identifiers_to_fetch))
    {
      cerr << "Cannot determine identifiers to fetch from '" << cl[0] << "'\n";
      return 4;
    }
    clstart = 1;

    if (! replace_tdt_data_with_data_from_identifier_file)
      ;
    else if (identifiers_to_fetch.size () != replacement_info.size ())
    {
      cerr << "Bad news, read " << identifiers_to_fetch.size () << " identifiers, but " << replacement_info.size () << " replacement info\n";
      return 4;
    }
  }
  else
  {
    cerr << "Must specify file of identifiers and TDT file, or list of identifiers (-I)\n";
    usage (5);
  }

  if (fetch_via_regexp.active())
    ;
  else if (0 == identifiers_to_fetch.size ())
  {
    cerr << "No identifiers to fetch\n";
    return 5;
  }

  if (verbose && ! identifiers_to_fetch.empty())
    cerr << "Looking for " << identifiers_to_fetch.size () << " identifiers\n";

  if (cl.option_present ('X'))
  {
    IWString x = cl.string_value ('X');

    stream_for_tdts_not_in_identifier_file.open (x.null_terminated_chars ());

    if (! stream_for_tdts_not_in_identifier_file.good ())
    {
      cerr << "Cannot open stream for TDT's not in identifier file '" << x << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "TDT's not in identifier file written to '" << x << "'\n";
  }

  IWString_and_File_Descriptor stream_for_identifiers_not_in_tdt_file;

  if (cl.option_present ('Y'))
  {
    const char * y = cl.option_value ('Y');

    if (! stream_for_identifiers_not_in_tdt_file.open (y))
    {
      cerr << "Cannot open stream for identifiers not in TDT file '" << y << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Identifiers not in TDT file written to '" << y << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  assert (clstart >= 0);

  int rc = 0;
  for (int i = clstart; i < cl.number_elements (); i++)
  {
    if (! fetch_tdt (cl[i], identifiers_to_fetch, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (stream_for_identifiers_not_in_tdt_file.is_open () || verbose)
  {
    extending_resizable_array<int> occurrences_of_identifier;

    for (IW_STL_Hash_Map_int::const_iterator f = identifiers_to_fetch.begin (); f != identifiers_to_fetch.end (); f++)
    {
      int times_written = (*f).second;

      occurrences_of_identifier[times_written]++;

      if (times_written)
        continue;

      if (stream_for_identifiers_not_in_tdt_file.is_open ())
      {
        stream_for_identifiers_not_in_tdt_file << (*f).first << '\n';
        stream_for_identifiers_not_in_tdt_file.write_if_buffer_holds_more_than(write_buffer_size);
      }
    }

    if (stream_for_identifiers_not_in_tdt_file.is_open())
      stream_for_identifiers_not_in_tdt_file.close();

    if (verbose)
      report_occurrences_of_identifiers (occurrences_of_identifier, cerr);
  }

  if (verbose)
  {
    cerr << "read " << tdts_read << " TDT's, wrote " << identifiers_written << " identifiers and " << tdts_written << " tdts\n";
    if (fetch_via_regexp.active())
      fetch_via_regexp.report(cerr);
  }

  if (cl.option_present('q'))
  {
    if (verbose)
      cerr << "Quick exit via -q option\n";
    _exit(0);
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fetch_tdt (argc, argv);

  return rc;
}
