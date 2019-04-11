/*
  Converts a TDT file to smiles
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_tdt.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = NULL;

static IWString smiles_tag ("$SMI");
static IWString identifier_tag ("PCN");

static IW_Regular_Expression also_match;

static int number_smiles_per_tdt_to_write = 0;

static int verbose = 0;

static int tdts_read = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "tdt2smi.cc,v 1.1.1.1 2003/09/25 12:35:55\n";
  cerr << "Converts TDT file to smiles\n";
  cerr << " -S <tag>       smiles tag (default '" << smiles_tag << "')\n";
  cerr << " -I <tag>       identifier tag (default '" << identifier_tag << "')\n";
  cerr << " -R <rx>        also write tags that match <rx>\n";
  cerr << " -n             each TDT may have multiple smiles in it\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static IWString output_buffer;

/*static int
append_data (IWString & buffer,
             int chars_to_strip,
             const_IWSubstring & zdata)
{
  if (buffer.length () > 0)
    buffer.add (' ');

  zdata.remove_leading_chars (chars_to_strip);
  zdata.chop ();

  buffer << zdata;

  return 1;
}

static int
append_data (IWString & buffer,
             const_IWSubstring & zdata)
{
  int chars_to_strip = zdata.index ('<') + 1;

  return append_data (buffer, chars_to_strip, zdata);
}*/


static int
tdt2smi (IW_TDT & tdt,
         ostream & output)
{
  output_buffer.resize_keep_storage (0);

  const_IWSubstring ztag, zdata;

  int smiles_found = 0;

  int i = 0;

  while (tdt.next_dataitem_value (ztag, zdata, i))
  {
//  cerr << "Examining '" << ztag << "' value '" << zdata << "'\n";
    if (smiles_tag == ztag)
    {
      if (output_buffer.length () > 0)
      {
        output_buffer.add ('\n');
        output << output_buffer;

        smiles_found++;
        if (number_smiles_per_tdt_to_write > 0 && smiles_found >= number_smiles_per_tdt_to_write)
          return output.good ();

        output_buffer.resize_keep_storage (0);
      }
      output_buffer.append_with_spacer (zdata);
    }
    else if (identifier_tag == ztag)
      output_buffer.append_with_spacer (zdata);
    else if (also_match.active () && also_match.matches (ztag))
      output_buffer.append_with_spacer (zdata);
  }

  if (output_buffer.length () > 0)
  {
    output_buffer.add ('\n');

    output << output_buffer;
  }

  return output.good ();
}


static int 
tdt2smi (iwstring_data_source & input,
         ostream & output)
{
  output_buffer.resize_keep_storage (0);

  IW_TDT tdt;

  while (tdt.next (input))
  {
    tdts_read++;

    if (! tdt2smi (tdt, output))
      return 0;
  }

  return output.good ();
}

static int 
tdt2smi (const char * fname,
         ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tdt2smi (input, output);
}

static int
tdt2smi (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vS:I:R:n:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('S'))
  {
    smiles_tag = cl.string_value ('S');

    if (verbose)
      cerr << "Smiles tag '" << smiles_tag << "'\n";
  }

  if (cl.option_present ('I'))
  {
    identifier_tag = cl.string_value ('I');

    if (verbose)
      cerr << "Identifier tag '" << identifier_tag << "'\n";
  }

  if (cl.option_present ('R'))
  {
    IWString r = cl.string_value ('R');

    if (! also_match.set_pattern (r))
    {
      cerr << "Cannot parse also match regular expression '" << r << "'\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will also match records that match '" << also_match.source () << "'\n";
  }

  if (cl.option_present ('n'))
  {
    const_IWSubstring n = cl.string_value ('n');

    if ("all" == n)
    {
      number_smiles_per_tdt_to_write = -4;
      if (verbose)
        cerr << "All smiles in a TDT will be written\n";
    }
    else if (! n.numeric_value (number_smiles_per_tdt_to_write) || number_smiles_per_tdt_to_write < 1)
    {
      cerr << "Invalid number of smiles per tdt to write '" << n << "'\n";
    }
    else if (verbose)
    {
      cerr << "Will write a max of " << number_smiles_per_tdt_to_write << " smiles per input TDT\n";
    }
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! tdt2smi (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << tdts_read << " tdt's\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tdt2smi (argc, argv);

  return rc;
}
