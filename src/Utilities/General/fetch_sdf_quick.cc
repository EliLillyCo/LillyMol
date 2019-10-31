/*
  Fetches sd file records based on identifiers in another file
*/

#include <stdlib.h>
#include <fstream>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_stl_hash_map.h"
#include "fetch_via_regexp.h"

const char * prog_name = NULL;

static int verbose = 0;

static int identifier_column_in_identifier_file = 0;

static int identifier_column_in_sdf_name = -1;

static IWString identifier_tag;    //the angle brackets only

static IWString identifier_tag_in_sd_file;   // with the > ..*< stuff

static int ignore_duplicate_identifiers_in_identifier_file = 0;

static int duplicate_identifiers_in_identifier_file = 0;

static IWString_and_File_Descriptor stream_for_not_in_identifier_file;
static IWString_and_File_Descriptor stream_for_not_in_sd_file;

static int sdfile_records_read = 0;
static int items_written = 0;
static int not_in_identifier_file = 0;

static int include_identifier_file_info = 0;

static IWString tag_for_identifier_file_info(' ');

static int strip_leading_zeros_from_identifiers = 0;

static int erase_identifiers_written = 1;

static int swap_identifier_and_data = 0;

static int block_size = 32768;

static IWString record_separator ("$$$$");

static int input_is_sdf = 1;

static int suppress_normal_output = 0;

static int write_nostructs = 1;

static int nostructs_encountered = 0;

static int ignore_sdfs_with_no_identifier = 0;

/*
  Conformations might come in as ID_n
*/

static int input_contains_conformations = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Fetches records from one file based on identifiers in another file\n";
  cerr << prog_name << " identifier_file sd_file > newfile\n";
  cerr << " -c <col>       identifier column in identifier file\n";
  cerr << " -C <tag>       identifier tag in sd file - default is 1st line\n";
  cerr << " -f <col>       identifier column in name field of .sdf file (default whole record)\n";
  cerr << " -d             ignore duplicate identifiers in identifier file\n";
  cerr << " -g             ignore possibly invalid sd records with no identifier\n";
  cerr << " -a             write all instances of identifiers in sd file\n";
  cerr << "                by default, only the first is written\n";
  cerr << " -X <fname>     write sd records not in <identifier_file> to <fname>\n";
  cerr << " -Y <fname>     write identifiers not in <sd_file> to <fname>\n";
  cerr << " -o             do NOT write nostructs to the -X file\n";
  cerr << " -K <id>        single identifiers on command line\n";
  cerr << " -R <rx>        match via regular expressions\n";
  cerr << " -T <string>    tag to insert between record and info from identifier file\n";
  cerr << " -w             swap id and info in the -Y file (useful for smiles)\n";
  cerr << " -z             strip leading zero's from identifiers\n";
  cerr << " -p <string>    record separator (default '" << record_separator << "')\n";
  cerr << " -e             input contains conformers with names ID_n. Remove _n part to match names\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
separate_into_id_and_data_2 (const const_IWSubstring & buffer,
                             int col,
                             IWString & id,
                             IWString & zdata)
{
  if (0 == col)   // the most common case
  {
//  cerr << "Extracting column " << col << " from '" << buffer << "'\n";
    int ispace = buffer.index (' ');

    if (ispace < 0)   // just one token on the line
    {
      id = buffer;
      return 1;
    }

    const char * s = buffer.rawchars ();
    id.set (s, ispace);

    zdata.set (s + ispace +1, buffer.length () - ispace - 1);

//  cerr << "Split into '" << id << "' and '" << zdata << "'\n";

    return 1;
  }

  if (! buffer.word (col, id))
  {
    cerr << "Cannot extract column " << col << " from '" << buffer << "'\n";
    return 0;
  }

  zdata = buffer;
  zdata.remove_word (col);

  return 1;
}

static int
separate_into_id_and_data (const const_IWSubstring & buffer,
                           int col,
                           IWString & id,
                           IWString & zdata)
{
  if (! separate_into_id_and_data_2 (buffer, col, id, zdata))
    return 0;

  if (strip_leading_zeros_from_identifiers)
    id.remove_leading_chars ('0');

  return 1;
}

/*
  There are two very distinct functionings at work here.

  In the case where we are adding any extra info from the
  identifier file to the sdfiles written, the hash holds
  that extra information. In the case where we are not
  processing that information, we keep track of how
  many times that ID has been written
*/

class Identifier_Info_String
{
  private:
    IW_STL_Hash_Map_String _identifier_to_fetch;

  public:
    int initialise (const IWString &, const IWString &);

    int contains (const IWString & s) const { return _identifier_to_fetch.contains (s);}

    int size () const { return _identifier_to_fetch.size ();}

    void set_has_been_written (const IWString &);

    int write_unwritten_items (IWString_and_File_Descriptor & output) const;

    const IWString & data_for_id (const IWString & s) const;
};

int
Identifier_Info_String::initialise (const IWString & s1,
                                    const IWString & s2)
{
  _identifier_to_fetch[s1] = s2;

  return 1;
}

void
Identifier_Info_String::set_has_been_written (const IWString & s)
{
  if (erase_identifiers_written)
    _identifier_to_fetch.erase (s); 

  return;
}

int
Identifier_Info_String::write_unwritten_items (IWString_and_File_Descriptor & output) const
{
  if (verbose)
    cerr << _identifier_to_fetch.size () << " identifiers not in sd file\n";

  for (IW_STL_Hash_Map_String::const_iterator i = _identifier_to_fetch.begin (); i != _identifier_to_fetch.end (); i++)
  {
    const IWString & zdata = (*i).second;

    if (swap_identifier_and_data && zdata.length ())
      output << zdata << ' ' << (*i).first;
    else
    {
      output << (*i).first;
      if (zdata.length ())
        output << ' ' << zdata;
    }
    output << '\n';

    output.write_if_buffer_holds_more_than(block_size);
  }

  return output.good ();
}

const IWString &
Identifier_Info_String::data_for_id (const IWString & s) const
{
  IW_STL_Hash_Map_String::const_iterator f = _identifier_to_fetch.find (s);

  if (f == _identifier_to_fetch.end ())
  {
    cerr << "Identifier_Info_String::data_for_id:no data for '" << s << "'\n";
    abort ();
  }

  return (*f).second;
}

static IWString unused_srring;

class Identifier_Info_Count
{
  private:
    IW_STL_Hash_Map_int _identifier_to_fetch;

//  private functions

    void _set_has_been_written (const IWString & id);

  public:
    int initialise (const IWString &, const IWString &);

    int contains (const IWString & s) const { return _identifier_to_fetch.contains (s);}

    int size () const { return _identifier_to_fetch.size ();}

    void set_has_been_written (const IWString &);

    int write_unwritten_items (IWString_and_File_Descriptor & output) const;

    const IWString & data_for_id (const IWString & id) const { return unused_srring;}   // is never called

    int add(const IWString & s) { _identifier_to_fetch[s] = 1; return _identifier_to_fetch.size();}
};

int
Identifier_Info_Count::initialise (const IWString & s1, const IWString & s2)
{
  _identifier_to_fetch[s1] = 0;

  return 1;
}

void
Identifier_Info_Count::set_has_been_written (const IWString & id)
{
  if (input_contains_conformations)
  {
    IWString tmp(id);
    tmp.truncate_at_last('_');
    _set_has_been_written(tmp);
  }
  else
    _set_has_been_written(id);
}

void
Identifier_Info_Count::_set_has_been_written (const IWString & id)
{
  IW_STL_Hash_Map_int::iterator f = _identifier_to_fetch.find (id);

  if (f == _identifier_to_fetch.end ())
  {
    cerr << "Identifier_Info_Count::_set_has_been_written:no info on '" << id << "'\n";
    return;
  }

  if (erase_identifiers_written)
    _identifier_to_fetch.erase (id);
  else
    (*f).second++;

  return;
}

int
Identifier_Info_Count::write_unwritten_items (IWString_and_File_Descriptor & output) const
{
  int rc = 0;

  for (IW_STL_Hash_Map_int::const_iterator i = _identifier_to_fetch.begin (); i != _identifier_to_fetch.end (); i++)
  {
    if ((*i).second > 0)   // has been written
      continue;

    output << (*i).first << '\n';

    output.write_if_buffer_holds_more_than(block_size);

    rc++;
  }

  if (verbose)
    cerr << rc << " items in identifier file not in sd file\n";

  return output.good ();
}

template <typename T>
int
read_identifiers_to_fetch (const const_IWSubstring & buffer,
                           int record_number,
                           T & identifiers_to_fetch)
{
  IWString id, zdata;
  if (! separate_into_id_and_data (buffer, identifier_column_in_identifier_file, id, zdata))
    return 0;

  if (! identifiers_to_fetch.contains (id))
    return identifiers_to_fetch.initialise (id, zdata);

  if (ignore_duplicate_identifiers_in_identifier_file)
  {
    cerr << "Ignoring duplicate identifer '" << id << "', line " << record_number << "\n";
    duplicate_identifiers_in_identifier_file++;
    return 1;
  }

  cerr << "Duplicate identifier '" << id << "'\n";
  return 0;
}

template <typename T>
int
read_identifiers_to_fetch (iwstring_data_source & input,
                           T & identifiers_to_fetch)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! read_identifiers_to_fetch (buffer, input.lines_read (), identifiers_to_fetch))
    {
      cerr << "Fatal error processing identifier file '" << buffer << "', line " << input.lines_read () << endl;
      return 0;
    }
  }

  return identifiers_to_fetch.size ();
}

template <typename T>
int
read_identifiers_to_fetch (const char * fname,
                           T & identifiers_to_fetch)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_strip_trailing_blanks(1);

  return read_identifiers_to_fetch (input, identifiers_to_fetch);
}

static int
insert_identifier_file_info(IWString & sd,
                            const IWString & tag,
                            const IWString & zinfo)
{
  if (0 == zinfo.length ())
    return 1;

  if (! sd.ends_with("$$$$\n"))
  {
    cerr << "insert_identifier_file_info:sdf data does not end in $$$$\n";
    return 1;
  }

  IWString tmp;
  tmp << ">  <" << tag << ">\n";
  tmp << zinfo << '\n';
  tmp << '\n';

  sd.insert(tmp, sd.length() - 5);

  return 0;
}

static int
looks_like_sdf_id (const const_IWSubstring & buffer,
                   const IWString & tag)
{
  if (! buffer.starts_with("> "))   // must be at least that much present
    return 0;

  if (buffer.ends_with(tag))    // good enough
    return 1;

  if (! buffer.contains(tag))
    return 0;

// Maybe try some more things some time later

  return 0;
}

static int
read_next_sd(iwstring_data_source & input,
             IWString & sd,
             IWString & id,
             int & nostruct)
{
  sd.resize_keep_storage(0);
  id.resize_keep_storage(0);

  nostruct = 0;

  const_IWSubstring buffer;

  int got_record_separator = 0;
  int next_record_is_id = 0;
  int lines_read = 0;

  while (input.next_record(buffer))
  {
    lines_read++;

    sd << buffer << '\n';

    if (record_separator == buffer)
    {
      got_record_separator = 1;
      break;
    }

    if (! input_is_sdf)
      ;
    else if (4 == lines_read && buffer.starts_with("  0  0"))
    {
      nostruct = 1;
      nostructs_encountered++;
      continue;
    }

    if (0 == identifier_tag_in_sd_file.length())
    {
      if (1 == lines_read)
        id = buffer;
    }
    else if (buffer.starts_with(identifier_tag_in_sd_file) || looks_like_sdf_id(buffer, identifier_tag))
    {
      if (id.length ())
        cerr << "yipes, duplicate " << identifier_tag_in_sd_file << "', id '" << id << "'\n";
      next_record_is_id = 1;
    }
    else if (next_record_is_id)
    {
      id = buffer;
      next_record_is_id = 0;
    }
  }

  if (0 == lines_read)   // normal eof
    return 0;

  if (! got_record_separator)
  {
    cerr << "Premature EOF\n";
    return 0;
  }

  if (0 == id.length())
  {
    cerr << "No identifier '" << identifier_tag_in_sd_file << "'\n";
    cerr << sd;
    return ignore_sdfs_with_no_identifier;
  }

  if (identifier_column_in_sdf_name < 0)
    return 1;

  if (0 == identifier_column_in_sdf_name)
  {
    id.truncate_at_first(' ');
    return 1;
  }

  if (id.nwords() <= identifier_column_in_sdf_name)
  {
    cerr << "Identifier column " << (identifier_column_in_sdf_name+1) << " but id '" << id << "' does not have enough tokens\n";
    return 0;
  }

  id.remove_leading_words(identifier_column_in_sdf_name);

  return 1;
}

int
fetch_sdf_quick_rx (IWString & sd,
                    const IWString & id,
                    int nostruct,
                    Fetch_via_Regexp & fetch_via_regexp,
                    IWString_and_File_Descriptor & output)
{
  int to_be_output = fetch_via_regexp.matches(id);

  if (! to_be_output)   // identifier not requested
  {
    if (verbose > 2)
      cerr << "SD file identifier '" << id << "' not requested\n";

    if (nostruct && ! write_nostructs)
      ;
    else
    {
      not_in_identifier_file++;
      if (stream_for_not_in_identifier_file.is_open())
      {
        stream_for_not_in_identifier_file << sd;
        stream_for_not_in_identifier_file.write_if_buffer_holds_more_than(block_size);
      }
    }

    return 1;
  }

  if (suppress_normal_output)
    ;
  else if (nostruct && ! write_nostructs)
    ;
  else
    output << sd;

  items_written++;

  return 1;
}

template <typename T>
int
fetch_sdf_quick(IWString & sd,
                const IWString & id,
                int nostruct,
                T & identifiers_to_fetch,
                IWString_and_File_Descriptor & output)
{
  int to_be_output = 0;

  if (input_contains_conformations)
  {
    IWString tmp(id);
    tmp.truncate_at_last('_');
    to_be_output = identifiers_to_fetch.contains(tmp);
  }
  else
    to_be_output = identifiers_to_fetch.contains(id);

  if (! to_be_output)   // identifier not requested
  {
    if (verbose > 2)
      cerr << "SD file identifier '" << id << "' not requested\n";

    if (nostruct && ! write_nostructs)
      ;
    else
    {
      not_in_identifier_file++;
      if (stream_for_not_in_identifier_file.is_open())
      {
        stream_for_not_in_identifier_file << sd;
        stream_for_not_in_identifier_file.write_if_buffer_holds_more_than(block_size);
      }
    }

    return 1;
  }

  if (suppress_normal_output)
    ;
  else if (nostruct && ! write_nostructs)
    ;
  else
  {
    if (tag_for_identifier_file_info.length() > 0)
      insert_identifier_file_info(sd, tag_for_identifier_file_info, identifiers_to_fetch.data_for_id(id));

    output << sd;
  }

  items_written++;

  identifiers_to_fetch.set_has_been_written(id);

  return 1;
}

int
fetch_sdf_quick_rx (iwstring_data_source & input,
                    Fetch_via_Regexp & fetch_via_regexp,
                    IWString_and_File_Descriptor & output)
{
  input.set_dos(1);
  input.set_strip_trailing_blanks(1);

  IWString sd;
  IWString id;
  int nostruct;

  while (read_next_sd (input, sd, id, nostruct))
  {
    sdfile_records_read++;

    if (0 == id.length ())
    {
      cerr << "Possibly erroneous structure, ending line " << input.lines_read () << endl;
      continue;
    }

    if (! fetch_sdf_quick_rx (sd, id, nostruct, fetch_via_regexp, output))
      return 0;

    output.write_if_buffer_holds_more_than(block_size);
  }

  return 1;
}

template <typename T>
int
fetch_sdf_quick(iwstring_data_source & input,
                T & identifiers_to_fetch,
                IWString_and_File_Descriptor & output)
{
  input.set_dos(1);
  input.set_strip_trailing_blanks(1);

  IWString sd;
  IWString id;
  int nostruct;

  while (read_next_sd(input, sd, id, nostruct))
  {
    sdfile_records_read++;

    if (0 == id.length ())
    {
      cerr << "Possibly erroneous structure, ending line " << input.lines_read () << endl;
      continue;
    }

    if (! fetch_sdf_quick(sd, id, nostruct, identifiers_to_fetch, output))
      return 0;

    output.write_if_buffer_holds_more_than(block_size);

    if (0 == identifiers_to_fetch.size())
      return 1;
  }

  return 1;
}

template <typename T>
int
fetch_sdf_quick (const char * fname,
                 T & identifiers_to_fetch,
                 IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return fetch_sdf_quick (input, identifiers_to_fetch, output);
}

static int
fetch_sdf_quick_rx (const char * fname,
                    Fetch_via_Regexp & fetch_via_regexp,
                    IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return fetch_sdf_quick_rx(input, fetch_via_regexp, output);
}

static void
summary_after_fetching (std::ostream & output)
{
  if (verbose)
  {
    output << sdfile_records_read << " sd file records read, " << items_written << " items written\n";
    if (not_in_identifier_file)
      output << not_in_identifier_file << " items not in identifier file\n";
    if (nostructs_encountered)
      output << "Encountered " << nostructs_encountered << " nostructs\n";
  }

  return;
}

static int
fetch_sdf_quick_rx (const Command_Line & cl,
                    Fetch_via_Regexp & fetch_via_regexp,
                    IWString_and_File_Descriptor & output)
{
  if (cl.option_present ('X'))
  {
    const char * fname = cl.option_value ('X');

    if (! stream_for_not_in_identifier_file.open (fname))
    {
      cerr << "Cannot open stream for not in identifier file (-X option), '" << fname << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Smiles not in identifier file written to '" << fname << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); ++i)
  {
    if (! fetch_sdf_quick_rx(cl[i], fetch_via_regexp, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (0 == items_written)
    cerr << "Warning, nothing written!\n";

  if (verbose)
    summary_after_fetching(std::cerr);

  return rc;
}


template <typename T>
int
fetch_sdf_quick (const Command_Line & cl,
                 T & identifiers_to_fetch,
                 IWString_and_File_Descriptor & output)
{
  if (cl.option_present ('X'))
  {
    const char * fname = cl.option_value ('X');

    if (! stream_for_not_in_identifier_file.open (fname))
    {
      cerr << "Cannot open stream for not in identifier file (-X option), '" << fname << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Smiles not in identifier file written to '" << fname << "'\n";
  }

  if (cl.option_present ('Y'))
  {
    const char * fname = cl.option_value ('Y');

    if (! stream_for_not_in_sd_file.open (fname))
    {
      cerr << "Cannot open stream for not in smiles file (-Y option), '" << fname << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Identifiers not in smiles file written to '" << fname << "'\n";
  }

  int rc = 0;
  for (int i = 1; i < cl.number_elements (); i++)
  {
    if (! fetch_sdf_quick (cl[i], identifiers_to_fetch, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (0 == items_written)
    cerr << "Warning, nothing written!\n";

  if (verbose)
    summary_after_fetching(std::cerr);

  if (stream_for_not_in_sd_file.is_open () && identifiers_to_fetch.size ())
    identifiers_to_fetch.write_unwritten_items (stream_for_not_in_sd_file);

  return rc;
}

static int
fetch_sdf_quick (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vc:C:X:Y:dzB:aT:wzogK:p:ef:R:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present ('a'))
  {
    erase_identifiers_written = 0;

    if (cl.option_present ('Y'))
    {
      cerr << "Sorry, the -a and -Y options do not work together\n";
      return 6;
    }

    if (verbose)
      cerr << "Will write all instances of identifiers in the smiles file\n";
  }

  if (cl.option_present ('T'))
  {
    cl.value('T', tag_for_identifier_file_info);

    if (verbose)
      cerr << "Will put '" << tag_for_identifier_file_info << " before identifier file info\n";

    include_identifier_file_info = 1;
  }

  if (cl.option_present ('B'))
  {
    if (! cl.value ('B', block_size) || block_size < 1)
    {
      cerr << "INvalid block size (-B option)\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will use " << block_size << " as the output block size\n";
  }

  if (cl.option_present('p'))
  {
    cl.value('p', record_separator);

    input_is_sdf = 0;

    if (verbose)
      cerr << "Record separator set to '" << record_separator << "'\n";
  }

  if (cl.option_present ('d'))
  {
    ignore_duplicate_identifiers_in_identifier_file = 1;

    if (verbose)
      cerr << "Will ignore duplicates in identifier file\n";
  }

  if (cl.option_present ('g'))
  {
    ignore_sdfs_with_no_identifier = 1;

    if (verbose)
      cerr << "Will ignore possibly invalid sdf records with missing or empty identifier\n";
  }

  if (cl.option_present ('z'))
  {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Will strip leading zero's from identifiers\n";
  }

  if (cl.option_present ('o'))
  {
    write_nostructs = 0;
    if (verbose)
      cerr << "Will not write nostructs\n";
  }

  if (cl.option_present ('c'))
  {
    if (! cl.value ('c', identifier_column_in_identifier_file) || identifier_column_in_identifier_file < 1)
    {
      cerr << "The column in identifier file option (-c) must be a whole +ve number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Identifiers in identifier file in column " << identifier_column_in_identifier_file << endl;

    identifier_column_in_identifier_file--;
  }

  if (cl.option_present ('C'))
  {
    const_IWSubstring c;

    cl.value ('C', c);

    identifier_tag << '<' << c << '>';

    identifier_tag_in_sd_file << ">  <" << c << ">";

    if (verbose)
      cerr << "Identifiers in sd file in '" << identifier_tag_in_sd_file << "' tag\n";
  }

  if (cl.option_present('f'))
  {
    if (! cl.value('f', identifier_column_in_sdf_name) || identifier_column_in_sdf_name < 1)
    {
      cerr << "The identifier column in .sdf file name (-f) must be a valid column number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Identifier column in .sdf file in column " << identifier_column_in_sdf_name << endl;

    identifier_column_in_sdf_name--;
  }

  if (cl.option_present ('w'))
  {
    if (! cl.option_present ('Y'))
    {
      cerr << "The -w option only makes sense with the -Y option\n";
      usage (3);
    }

    swap_identifier_and_data = 1;

    if (verbose)
      cerr << "Will swap positions of identifier and data in -Y file\n";
  }

  if (cl.option_present('e'))
  {
    input_contains_conformations = 1;

    if (verbose)
      cerr << "Input assumed to be conformers\n";

    if (! cl.option_present('a'))
      cerr << "The -a option not specified, only first conformation of each molecule written\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('K'))
  {
    Identifier_Info_Count identifiers_to_fetch;
    int i = 0;
    IWString k;
    while (cl.value('K', k, i++))
    {
      identifiers_to_fetch.add(k);
    }

    if (verbose)
      cerr << "Will fetch " << identifiers_to_fetch.size () << " identifiers\n";

//  the underlying method starts processing at the 2nd input file. Insert
//  a dummy thing in there

    const char * notused = "notused";
    cl.insert_at_beginning(notused);

    rc = fetch_sdf_quick (cl, identifiers_to_fetch, output);
  }
  else if (cl.option_present('R'))
  {
    Fetch_via_Regexp fvr;

    if (! fvr.build(cl, 'R'))
    {
      cerr << "Cannot build regular expressions to fetch (-R)\n";
      return 1;
    }

    rc = fetch_sdf_quick_rx(cl, fvr, output);
  }
  else if (cl.number_elements() < 2)
  {
    cerr << "Must specify at least two files\n";
    usage (2);
  }
  else if (include_identifier_file_info)
  {
    Identifier_Info_String identifiers_to_fetch;

    if (! read_identifiers_to_fetch(cl[0], identifiers_to_fetch))
    {
      cerr << "Cannot read identifiers to fetch from '" << cl[0] << "'\n";
      return 6;
    }

    if (verbose)
      cerr << "Will fetch " << identifiers_to_fetch.size() << " identifiers\n";

    rc = fetch_sdf_quick(cl, identifiers_to_fetch, output);
  }
  else
  {
    Identifier_Info_Count identifiers_to_fetch;

    if (! read_identifiers_to_fetch (cl[0], identifiers_to_fetch))
    {
      cerr << "Cannot read identifiers to fetch from '" << cl[0] << "'\n";
      return 6;
    }

    if (verbose)
      cerr << "Will fetch " << identifiers_to_fetch.size () << " identifiers\n";

    rc = fetch_sdf_quick (cl, identifiers_to_fetch, output);
  }

  output.flush();

  if (stream_for_not_in_sd_file.is_open())
    stream_for_not_in_sd_file.close();

  if (stream_for_not_in_identifier_file.is_open())
    stream_for_not_in_identifier_file.close();

  return ! rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fetch_sdf_quick (argc, argv);

  return rc;
}
