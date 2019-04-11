/*
  
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iw_stl_hash_map.h"
#include "iwstring_data_source.h"

using std::cerr;
using std::cout;
using std::ostream;

const char * prog_name = NULL;

static int verbose = 0;

static IWString descriptor_tag ("DSC<");

static int header_records_to_skip = 1;

static int strip_leading_zeros = 0;

static IWString identifier_tag ("PCN<");

static int truncate_descriptor_file_identifiers_at_first_underscore = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "gfp_add_descriptors.cc,v 1.1 2004/03/22 13:05:43 \n";
  cerr << "Adds a descriptor file to a gfp file\n";
  cerr << prog_name << " gfp_file descriptor_file > newfile\n";
  cerr << " -D <tag>       tag to use for descriptors ( default '" << descriptor_tag << "')\n";
  cerr << " -z             strip leading '0' characters when matching identifiers\n";
  cerr << " -I <tag>       identifier tag (default '" << identifier_tag << "'\n";
  cerr << " -u             truncate descriptor file identifiers at first underscore\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
fetch_descriptors (const const_IWSubstring & id, 
                   IW_STL_Hash_Map_String & descriptors,
                   ostream & output)
{
  IWString tmp (id);

  if (strip_leading_zeros)
    tmp.remove_leading_chars ('0');

  tmp.truncate_at_first (' ');

  IW_STL_Hash_Map_String::const_iterator f = descriptors.find (tmp);

  if (f == descriptors.end ())
  {
    cerr << "No descriptors for '" << tmp << "' from '" << id << "'\n";
    return 0;
  }

  output << descriptor_tag << (*f).second << ">\n";

  return output.good ();
}

static int
gfp_add_descriptors (iwstring_data_source & input,
                     IW_STL_Hash_Map_String & descriptors,
                     ostream & output)
{
  const_IWSubstring buffer;

  int got_id = 0;

  int tdts_read = 0;

  while (input.next_record (buffer))
  {
    output << buffer << endl;

    if ('|' == buffer)
    {
      if (! got_id)
      {
        cerr << "TDT with no identifier, impossible\n";
        return 0;
      }

      got_id = 0;

      tdts_read++;
    }
    else if (buffer.starts_with (identifier_tag))
    {
      if (got_id)
        cerr << "HUH, duplicate identifiers in TDT!!\n";

      got_id = 1;

      buffer.remove_leading_chars (identifier_tag.length ());
      buffer.chop ();    // remove >

      if (! fetch_descriptors (buffer, descriptors, output))
      {
        cerr << "Fatal error fetching descriptors for '" << buffer << "'\n";
        return 0;
      }
    }
  }

  return output.good ();
}

static int
gfp_add_descriptors (const char * fname,
                     IW_STL_Hash_Map_String & descriptors,
                     ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_add_descriptors (input, descriptors, output);
}

static int
read_descriptors (const const_IWSubstring & buffer,
                  IW_STL_Hash_Map_String & descriptors)
{
  const_IWSubstring id, zdata;

  if (! buffer.split (id, ' ', zdata))
  {
    cerr << "Invalid descriptor record\n";
    return 0;
  }

  if (strip_leading_zeros)
    id.remove_leading_chars ('0');

  if (truncate_descriptor_file_identifiers_at_first_underscore)
    id.truncate_at_first ('_');

  if (descriptors.contains (id))
    cerr << "Ignoring duplicate id in descriptor file '" << id << "'\n";
  else
    descriptors[id] = zdata;

  return 1;
}

static int
read_descriptors (iwstring_data_source & input,
                  IW_STL_Hash_Map_String & descriptors)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input.next_record (buffer))
    {
      cerr << "EOF reading descriptors\n";
      return 0;
    }
  }

  while (input.next_record (buffer))
  {
    if (! read_descriptors (buffer, descriptors))
    {
      cerr << "Fatal error reading line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return descriptors.size ();
}

static int
read_descriptors (const char * fname,
                  IW_STL_Hash_Map_String & descriptors)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_descriptors (input, descriptors);
}


static int
gfp_add_descriptors (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vD:zI:u");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('D'))
  {
    cl.value ('D', descriptor_tag);

    if (verbose)
      cerr << "Descriptors written with tag '" << descriptor_tag << "'\n";

    if (! descriptor_tag.ends_with ('<'))
      descriptor_tag.add ('<');
  }

  if (cl.option_present ('I'))
  {
    cl.value ('I', identifier_tag);

    if (verbose)
      cerr << "Identifiers in gfp file '" << identifier_tag << "'\n";

    if (! identifier_tag.ends_with ('<'))
      identifier_tag.add ('<');
  }

  if (cl.option_present ('z'))
  {
    strip_leading_zeros = 1;

    if (verbose)
      cerr << "Will remove leading '0' characters when matching id's\n";
  }

  if (cl.option_present ('u'))
  {
    truncate_descriptor_file_identifiers_at_first_underscore = 1;

    if (verbose)
      cerr << "Will truncate descriptor file identifiers at first underscore\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (2 != cl.number_elements ())
  {
    cerr << "Must specify exactly two files, gfp file and descriptor file\n";
    usage (4);
  }

  IW_STL_Hash_Map_String descriptors;

  if (! read_descriptors (cl[1], descriptors))
  {
    cerr << "Cannot read descriptors from '" << cl[0] << "'\n";
    return 6;
  }

  if (verbose)
    cerr << "Read " << descriptors.size () << " descriptor values\n";

  if (! gfp_add_descriptors (cl[0], descriptors, cout))
  {
    cerr << "Fatal error merging fingerprints and descriptors\n";
    return 4;
  }

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_add_descriptors (argc, argv);

  return rc;
}
