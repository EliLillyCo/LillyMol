/*
  Converts a .rpt file to a .nn file
*/

#include <stdlib.h>
#include <iostream>
#include <limits>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Utilities/GFP_Tools/smiles_id_dist.h"

using std::cerr;
using std::endl;
using std::numeric_limits;

const char * prog_name = nullptr;

static int verbose = 0;

static int neighbours_to_write = numeric_limits<int>::max();

static int column_to_process = -1;

static IWString tag;

static int comboscore = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Converts a ROCS .rpt file to a .nn file\n";
  cerr << " -S <fname>     file with smiles\n";
  cerr << " -d <tag>       header tag to process\n";
  cerr << " -n <nbrs>      number neighbours to write\n";
  cerr << " -h             discard self neighbours\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static IW_STL_Hash_Map_String smiles;

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  IWString s, id;

  while (input.next_record(buffer))
  {
    if (buffer.nwords() < 2)
    {
      cerr << "Smiles file must have at least two tokens\n";
      return 0;
    }

    int i = 0;
    buffer.nextword(s, i);
    buffer.nextword(id, i);

    smiles[id] = s;
  }

  return smiles.size();
}

static int
read_smiles (const char * fname,
             IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open smiles file '" << fname << "'\n";
    return 0;
  }

  return read_smiles (input, smiles);
}

static int
rocs_rpt2nn_record (const const_IWSubstring & buffer,
                    resizable_array_p<Smiles_ID_Dist> & sid)
{
  if (buffer.nwords() < 12)
  {
    cerr << "ROCS .rpt files have 12 tokens\n";
    return 0;
  }

  IWString haystack_id, needle_id;

  int i = 0;

  buffer.nextword(haystack_id, i);
  buffer.nextword(needle_id, i);

  IWString haystack_smiles;

  if (! smiles.contains(haystack_id))
  {
    cerr << "No smiles for '" << haystack_id << "'\n";
    return 0;
  }

  haystack_smiles = smiles[haystack_id];

  const_IWSubstring token;

  for (int col = 2 ; buffer.nextword(token, i); col++)
  {
    if (col != column_to_process)
      continue;

    float d;
    if (! token.numeric_value(d))
    {
      cerr << "Invalid numeric '" << token << "'\n";
      return 0;
    }

    if (comboscore)
      d = (2.0 - d) * 0.5;
    else
      d = 1.0 - d;

    Smiles_ID_Dist * s = new Smiles_ID_Dist(haystack_smiles, haystack_id, d);
    sid.add(s);
  }

  return sid.size();
}

class Smiles_ID_Distance_Comparator
{
  private:
  public:
    int operator() (const Smiles_ID_Dist *, const Smiles_ID_Dist *) const;
};

int
Smiles_ID_Distance_Comparator::operator() (const Smiles_ID_Dist * s1, const Smiles_ID_Dist * s2) const
{
  float d1 = s1->distance();
  float d2 = s2->distance();

  if (d1 < d2)
    return -1;

  if (d1 > d2)
    return 1;

  return 0;
}

static Smiles_ID_Distance_Comparator sid_comparator;

static int
write_needle_smiles_and_id (const const_IWSubstring & buffer,
                            const IW_STL_Hash_Map_String & smiles,
                            IWString_and_File_Descriptor & output)
{
  int i = 0;
  const_IWSubstring token;

  if (! buffer.nextword(token, i))   // haystack molecule
    return 0;

  IWString needle_id;
  if (! buffer.nextword(needle_id, i))
    return 0;

  IW_STL_Hash_Map_String::const_iterator f = smiles.find(needle_id);

  if (f == smiles.end())
  {
    cerr << "No smiles for needle '" << needle_id << "'\n";
    return 0;
  }

  output << smiles_tag << (*f).second << ">\n";
  output << identifier_tag << needle_id << ">\n";

  return 1;
}

static int
determine_column_to_process (const const_IWSubstring & buffer,
                             int & column_to_process,
                             const IWString & tag)
{
  const_IWSubstring token;

  int col = 0;
  for (int i = 0; buffer.nextword(token, i); col++)
  {
    if (token == tag)
    {
      column_to_process = col;
      return 1;
    }
  }

  return 0;
}


static int
rocs_rpt2nn (iwstring_data_source & input,
             IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  if (column_to_process > 0)
    ;
  else if (! determine_column_to_process(buffer, column_to_process, tag))
  {
    cerr << "Cannot find '" << tag << "' in header record\n";
    return 0;
  }

  resizable_array_p<Smiles_ID_Dist> sid;

  sid.resize(20000);

  while (input.next_record(buffer))
  {
    if (sid.number_elements())   // look for needle on first item read
      ;
    else if (! write_needle_smiles_and_id(buffer, smiles, output))
    {
      cerr << "Cannot fetch data about needle\n";
      return 0;
    }

    if (! rocs_rpt2nn_record (buffer, sid))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  int n = sid.number_elements();

  if (0 == n)
    return 0;

  if (verbose)
    cerr << "Read " << n << " records, sorting...\n";

  if (n > 1)
    sid.iwqsort (sid_comparator);

  if (verbose > 1)
    cerr << "Sort complete\n";

  if (n > neighbours_to_write)
    n = neighbours_to_write;

  for (int i = 0; i < n; i++)
  {
    const Smiles_ID_Dist * s = sid[i];

    output << smiles_tag << s->smiles() << ">\n";
    output << identifier_tag << s->id() << ">\n";
    output << distance_tag << s->distance() << ">\n";

    output.write_if_buffer_holds_more_than(32768);
  }

  output << "|\n";

  return 1;
}

static int
rocs_rpt2nn (const char * fname,
             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return rocs_rpt2nn(input, output);
}

static int
rocs_rpt2nn (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vS:d:n:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('d'))
  {
    cerr << "Must specify which tag to process via the -d option\n";
    usage(3);
  }

  if (cl.option_present('d'))
  {
    cl.value('d', tag);

    if ("ComboScore" == tag)
      comboscore = 1;

    if (verbose)
      cerr << "Examining tag '" << tag << "\n";
  }

  if (! cl.option_present('S'))
  {
    cerr << "Must specify smiles file via the -S option\n";
    usage(3);
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', neighbours_to_write) || neighbours_to_write < 1)
    {
      cerr << "The number of neighbours_to_write must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << " Will write " << neighbours_to_write << " neighbours\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('S'))
  {
    const char * s = cl.option_value('S');
    if (! read_smiles(s, smiles))
    {
      cerr << "Cannot read smiles from '" << s << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Read " << smiles.size() << " smiles from '" << s << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! rocs_rpt2nn(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rocs_rpt2nn(argc, argv);

  return rc;
}
