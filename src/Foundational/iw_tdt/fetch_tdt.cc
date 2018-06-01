/*
  Fetches specific TDT's from a file of TDT's
  Output in the same order that they are in the identifier file
*/

#include <stdlib.h>
#include <fstream>
using namespace std;

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_stl_hash_map.h"
#include "iwcrex.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString identifier_tag ("PCN<");

static IW_Regular_Expression identifier_rx;

static int ignore_identifiers_not_matching_rx = 0;

static int identifier_column_in_identifier_file = -1;

static int gsub_underscores_in_identifier_file = 0;

static int identifier_column_in_TDT = -1;

static int ignore_duplictes_in_identifier_file = 0;

static int ignore_tdts_missing_the_tag = 0;

static int identifiers_read = 0;

static int tdts_written = 0;

static int identifiers_written = 0;

static int write_all_instances_of_duplicates = 0;

static IWString_and_File_Descriptor stream_for_identifiers_not_in_tdt;

static int replace_tdt_data_with_data_from_identifier_file = 0;

static int remove_leading_zeros = 0;

static int invert_fetching_operation = 0;

/*
  We need to signal that an identifier has been written. Byte offset
  1 is an impossible offset in a TDT file
*/

#define IMPOSSIBLE_OFFSET static_cast<off_t> (1)

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Fetches a subset of TDT's from a file of TDT's\n";
  cerr << prog_name << " options <identifier_file> <tdt_file>\n";
  cerr << " -T <tag>       tag for identifiers\n";
  cerr << " -p <rx>        identifiers must match regular expression <rx>\n";
  cerr << " -i             ignore identifiers not matching <rx>\n";
  cerr << " -d ignore      ignore duplicates in <identifier_file>\n";
  cerr << " -D all         fetch all instances of an identifier in <tdt_file>\n";
  cerr << " -c <col>       identifier in column <col> in the identifier file\n";
  cerr << " -C <col>       identifier in column <col> in the tdt file\n";
  cerr << " -u             gsub _ to space in identifier file\n";
  cerr << " -X <fname>     write tdts not in <identifier file> to <fname>\n";
  cerr << " -Y <fname>     write identifiers not in <tdt_file> to <fname>\n";
  cerr << " -r             replace identifier in TDT with data from identifier file\n";
  cerr << " -x             invert behaviour, selection becomes deselection\n";
  cerr << " -z             remove leading 0's from identifiers\n";
  cerr << " -q             quick exit, don't wait for hash deallocation\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
handle_identifiers_not_in_tdt_file (const IWString & id)
{
  if (stream_for_identifiers_not_in_tdt.is_open ())
  {
    stream_for_identifiers_not_in_tdt << id << '\n';
    stream_for_identifiers_not_in_tdt.write_if_buffer_holds_more_than(16000);
    return stream_for_identifiers_not_in_tdt.good ();
  }

  return 1;
}

static int
extract_identifier (const IWString & buffer,
                    IWString & id,
                    int & fatal)
{
  if (identifier_column_in_identifier_file >= 0)
  {
    if (buffer.nwords () < identifier_column_in_identifier_file + 1)
    {
      cerr << "Sorry, cannot extract column " << (identifier_column_in_identifier_file + 1) << " from identifier record\n";
      fatal = 1;
      return 0;
    }

    (void) buffer.word (identifier_column_in_identifier_file, id);
  }
  else
    id = buffer;

  if (remove_leading_zeros)
    id.remove_leading_chars('0');

  if (! identifier_rx.active ())
    ;
  else if (identifier_rx.matches (id))
    ;
  else if (ignore_identifiers_not_matching_rx)
  {
    cerr << "Ignoring identifier '" << id << "' not matching '" << identifier_rx.source () << "'\n";
    fatal = 0;
    return 0;
  }
  else
  {
    cerr << "Invalid identifier '" << id << "' does not match '" << identifier_rx.source () << "'\n";
    fatal = 1;
    return 0;
  }

  return 1;
}

template <typename O> int
echo_tdt (iwstring_data_source & tdt_file,
          off_t o,
          const IWString & id,
          const const_IWSubstring & from_identifier_file,
          O & output_buffer)
{
  if (! tdt_file.seekg (o))
  {
    cerr << "Yipes, cannot seek to " << o << " for '" <<id << "'\n";
    return 0;
  }

  const_IWSubstring buffer;

  int identifier_replaced = 0;

  while (tdt_file.next_record (buffer))
  {
    if (replace_tdt_data_with_data_from_identifier_file && ! identifier_replaced && buffer.starts_with (identifier_tag))
    {
      output_buffer << identifier_tag << from_identifier_file << ">\n";
      identifier_replaced = 1;
    }
    else
      output_buffer << buffer << '\n';

    if ('|' == buffer)
      break;
  }

  return 1;
}

static int
fetch_tdt (const const_IWSubstring & buffer,
           iwstring_data_source & tdt_file,
           IW_STL_Hash_Map<IWString, off_t> & offset,
           IWString_and_File_Descriptor & output)
{
  IWString id;
  int fatal;

  if (! extract_identifier (buffer, id, fatal))
  {
    if (! fatal)
      return 1;

    cerr << "Missing or invalid identifier in '" << buffer << "'\n";
    return 0;
  }

  IW_STL_Hash_Map<IWString, off_t>::iterator f = offset.find (id);

  bool is_match;

  if (f == offset.end())
    is_match = false;
  else
    is_match = true;

  if (invert_fetching_operation)
  {
    if (is_match)
      is_match = false;
    else
      is_match = true;
  }

//cerr << "Do we want to retrieve '" << id << "'\n";

  if (! is_match)   // not interested in this one
    return handle_identifiers_not_in_tdt_file (id);

  off_t o = (*f).second;

  if (IMPOSSIBLE_OFFSET != o)   // never seen this one before
    identifiers_written++;
  else if (write_all_instances_of_duplicates)
    ;
  else
  {
    if (verbose)
      cerr << "Ignoring duplicate '" << id << "' in identifier file\n";
    return 1;
  }

  if (! echo_tdt (tdt_file, o, id, buffer, output))
    return 0;

  tdts_written++;

  if (write_all_instances_of_duplicates)     // keep the offset
    ;
  else
    offset[id] = IMPOSSIBLE_OFFSET;    // make sure it won't be written again

  return 1;
}

static int
fetch_tdt (iwstring_data_source & identifier_file,
           iwstring_data_source & tdt_file,
           IW_STL_Hash_Map<IWString, off_t> & offset,
           IWString_and_File_Descriptor & output)
{
  IWString buffer;

  identifier_file.set_dos (1);

  while (identifier_file.next_record (buffer))
  {
    if (gsub_underscores_in_identifier_file)
      buffer.gsub ('_', ' ');

    if (! fetch_tdt (buffer, tdt_file, offset, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(16000);

    identifiers_read++;
  }

  return 1;
}

static int
extract_identifier_from_tdt (const_IWSubstring buffer,    // not passed by reference, we have a copy
                             IWString & id,
                             int & fatal)
{
  buffer.remove_leading_chars (identifier_tag.length ());
  buffer.chop ();

  if (identifier_column_in_TDT >= 0)
  {
    if (! buffer.word (identifier_column_in_TDT, id))
    {
      cerr << "Cannot extract word " << identifier_column_in_TDT << " from '" << buffer << "'\n";
      fatal = 1;
      return 0;
    }

    if (remove_leading_zeros)
      id.remove_leading_chars('0');

    return 1;
  }

  id = buffer;

  if (remove_leading_zeros)
    id.remove_leading_chars('0');

  if (! identifier_rx.active ())
    ;
  else if (identifier_rx.matches (id))
    ;
  else if (ignore_identifiers_not_matching_rx)
  {
    cerr << "Ignoring identifier '" << id << "' not matching '" << identifier_rx.source () << "'\n";
    fatal = 0;
    return 0;
  }
  else
  {
    cerr << "Invalid identifier '" << id << "' does not match '" << identifier_rx.source () << "'\n";
    fatal = 1;
    return 0;
  }

  return 1;
}

static int
read_offsets (iwstring_data_source & tdt_file,
              IW_STL_Hash_Map<IWString, off_t> & offset)
{
  IWString buffer;

  off_t o = static_cast<off_t> (0);

  int got_id_this_tdt = 0;

  while (tdt_file.next_record (buffer))
  {
    if ('|' == buffer)
    {
      if (got_id_this_tdt)
        ;
      else if (ignore_tdts_missing_the_tag)
        ;
      else
      {
        cerr << "No '" << identifier_tag << "' tag in TDT, line " << tdt_file.lines_read () << endl;
        return 0;
      }

      o = tdt_file.tellg ();
      got_id_this_tdt = 0;

      continue;
    }

    if (got_id_this_tdt)
      continue;

//  if (0 == buffer.length())   // PD2 problem
//    continue;

    if (! buffer.starts_with (identifier_tag))
      continue;

    IWString id;
    int fatal;

    if (! extract_identifier_from_tdt (buffer, id, fatal))
    {
      if (! fatal)
        continue;

      cerr << "Missing or invalid identifier in '" << buffer << "'\n";
      return 0;
    }

    if (! offset.contains (id))
      offset[id] = o;    // set to the number fetched
    else if (ignore_duplictes_in_identifier_file)
      cerr << "Ignoring duplicate identifier '" << id << "'\n";
    else
    {
      cerr << "Duplicate identifier in ID file '" << id << "'\n";
      return 0;
    }

    got_id_this_tdt = 1;
  }

  return offset.size ();
}

static int
fetch_tdt (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vc:C:X:Y:p:id:D:T:vurzqx");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('c'))
  {
    if (! cl.value ('c', identifier_column_in_identifier_file) || identifier_column_in_identifier_file < 1)
    {
      cerr << "The identifier column in identifier file (-c) must be a valid column number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Identifiers found in column " << identifier_column_in_identifier_file << " of identifier file\n";

    identifier_column_in_identifier_file--;

    if (cl.option_present ('u'))
    {
      gsub_underscores_in_identifier_file = 1;

      if (verbose)
        cerr << "Will gsub underscores in identifier file\n";
    }
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

  if (cl.option_present ('p'))
  {
    const_IWSubstring p = cl.string_value ('p');

    if (! identifier_rx.set_pattern (p))
    {
      cerr << "Invalid identifier regexp '" << p << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Identifiers must match '" << p << "'\n";

    if (cl.option_present ('i'))
    {
      ignore_identifiers_not_matching_rx = 1;

      if (verbose)
        cerr << "Will ignore identifiers not matching RX\n";
    }
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

  if (cl.option_present('z'))
  {
    remove_leading_zeros = 1;

    if (verbose)
      cerr << "Will remove leading 0's from identifers\n";
  }

  if (cl.option_present('x'))
  {
    invert_fetching_operation = 1;

    if (verbose)
      cerr << "Will discard rather than fetch identifiers\n";
  }

  if (cl.option_present ('r'))
  {
    replace_tdt_data_with_data_from_identifier_file = 1;

    if (verbose)
      cerr << "Will replace identifier in TDT with data from identifier file\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (2 != cl.number_elements ())
  {
    cerr << "Must specify both an identifier file and a TDT file\n";
    usage (4);
  }

  iwstring_data_source identifier_file (cl[0]);

  if (! identifier_file.good ())
  {
    cerr << "Cannot open identifier file '" << cl[0] << "'\n";
    return 5;
  }

  iwstring_data_source tdt_file (cl[1]);

  if (! tdt_file.good ())
  {
    cerr << "Cannot open input file '" << cl[1] << "'\n";
    return 4;
  }

  IW_STL_Hash_Map<IWString, off_t> offset;

  if (! read_offsets (tdt_file, offset))
  {
    cerr << "Cannot determine TDT offsets in '" << cl[1] << "'\n";
    return 5;
  }

  if (verbose)
    cerr << "Input contains " << offset.size () << " tdts\n";

#ifdef DEBUG_FETCH_TDT
  for (IW_STL_Hash_Map<IWString, off_t>::const_iterator i = offset.begin(); i != offset.end(); ++i)
  {
    cerr << "Item '" << (*i).first << "' at " << (*i).second << endl;
  }
#endif

  if (cl.option_present ('Y'))
  {
    IWString y = cl.string_value ('Y');

    if (! stream_for_identifiers_not_in_tdt.open (y.null_terminated_chars ()))
    {
      cerr << "Cannot open stream for identifiers not in TDT file '" << y << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Identifiers's not in TDT file written to '" << y << "'\n";
  }

  IWString_and_File_Descriptor stream_for_tdts_not_in_identifier_file;

  if (cl.option_present ('X'))
  {
    IWString x = cl.string_value ('X');

    if (! stream_for_tdts_not_in_identifier_file.open (x.null_terminated_chars()))
    {
      cerr << "Cannot open stream for TDTs not in identifier file '" << x << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "TDT's not in identifier file written to '" << x << "'\n";
  }

  IWString_and_File_Descriptor output(1);
  if (! fetch_tdt (identifier_file, tdt_file, offset, output))
  {
    return 5;
  }

  output.flush();

  replace_tdt_data_with_data_from_identifier_file = 0;   // doesn't make sense now

  if (stream_for_tdts_not_in_identifier_file.is_open() || verbose)
  {
    int tdts_not_written = 0;

    for (IW_STL_Hash_Map<IWString, off_t>::const_iterator f = offset.begin (); f != offset.end (); f++)
    {
      off_t o = (*f).second;

      if (IMPOSSIBLE_OFFSET == o)    // means it was written
        continue;

      tdts_not_written++;

      if (! stream_for_tdts_not_in_identifier_file.is_open ())
        ;
      else if (echo_tdt (tdt_file, o, (*f).first, "", stream_for_tdts_not_in_identifier_file))
        ;
      else
      {
        cerr << "Cannot echo TDT for '" << (*f).first << "' at offset " << o << endl;
        return 6;
      }
    }

    if (stream_for_tdts_not_in_identifier_file.is_open())
      stream_for_tdts_not_in_identifier_file.close();

    if (verbose)
      cerr << tdts_not_written << " tdts were not written\n";
  }

  if (verbose)
    cerr << "read " << identifiers_read << " identifiers, wrote " << identifiers_written << " identifiers and " << tdts_written << " tdts\n";

  if (stream_for_identifiers_not_in_tdt.is_open())
    stream_for_identifiers_not_in_tdt.close();

  if (cl.option_present('q'))
    _exit(0);

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fetch_tdt (argc, argv);

  return rc;
}
