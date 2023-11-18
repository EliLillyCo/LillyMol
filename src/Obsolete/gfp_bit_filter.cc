/*
  Finds those fingerprints that match a given set
*/

#include <stdlib.h>

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString_and_File_Descriptor stream_for_not_selected;

static int write_tdt = 0;
static int write_smiles = 0;

static int tdts_read = 0;

static int tdts_matching = 0;

static int all_queries_must_match = 0;

static IWString smiles_tag("$SMI<");

class Match_Criterion
{
  private:
    IWString _tag;
    int _fixed;          // or sparse
    unsigned int _fingerprint;
    unsigned int _bit;

    int _comparisons_done;
    int _successful_matches;

//  private functions

    int _determine_fingerprint_numbers (const IW_General_Fingerprint &);

  public:
    Match_Criterion();

    int initialise (const const_IWSubstring &);

    int matches (const IW_General_Fingerprint &);
};

Match_Criterion::Match_Criterion()
{
  _fixed = -1;    // flag to examine this on first call

  _comparisons_done = 0;
  _successful_matches = 0;

  return;
}

int
Match_Criterion::initialise (const const_IWSubstring & s)
{
  const_IWSubstring tag, b;

  if (! s.split(tag, ',', b))
  {
    cerr << "Match_Criterion::initialise:invalid tag,bit specification '" << s << "'\n";
    return 0;
  }

  _tag = tag;
  if (_tag.ends_with('<'))
    _tag.chop();

  if (! b.numeric_value(_bit))
  {
    cerr << "Match_Criterion::initialise:invalid bit number '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
Match_Criterion::_determine_fingerprint_numbers (const IW_General_Fingerprint & fp)
{
  int n = number_fingerprints();

  for (int i = 0; i < n; i++)
  {
    const IWString & t = fixed_fingerprint_tag(i);

    if (_tag == t)
    {
      _fixed = 1;
      _fingerprint = i;
      return 1;
    }
  }

  n = number_sparse_fingerprints();

  for (int i = 0; i < n; i++)
  {
    const IWString & t = sparse_fingerprint_tag(i);
    cerr << "Checking tag '" << t << "'\n";

    if (_tag == t)
    {
      _fixed = 0;
      _fingerprint = i;
      return 1;
    }
  }

  cerr << "Match_Criterion::initialise:did not match '" << _tag << "' in input fingerprints\n";
  return 0;
}

int
Match_Criterion::matches (const IW_General_Fingerprint & fp)
{
  if (_fixed < 0)   // need initialising
  {
    if (! _determine_fingerprint_numbers(fp))
    {
      cerr << "Match_Criterion::matches:cannot initialise fingerprints\n";
      abort();
    }
  }

  _comparisons_done++;

  if (_fixed)
  {
    const IWDYFP & f = fp[_fingerprint];
    if (f.is_set(_bit))
    {
      _successful_matches++;
      return 1;
    }
  }
  else
  {
    const Sparse_Fingerprint & f = fp.sparse_fingerprint(_fingerprint);

    if (f.is_set(_bit))
    {
      _successful_matches++;
      return 1;
    }
  }

  return 0;
}

static int nq = 0;
static Match_Criterion * query = nullptr;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Filters a .gfp file by bits matched\n";
  cerr << " -s <tag,number>  selection criteria\n";
  cerr << " -a               all queries must match\n";
  cerr << " -W tdt           write tdt rather than smiles output\n";
  cerr << " -W smi           write tdt rather than smiles output\n";
  cerr << " -B <fname>       write items not selected to <fname>\n";
  cerr << " -F ...           standard gfp options\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
handle_non_match (const IW_TDT & tdt,
                  const IW_General_Fingerprint & fp)
{
  if (! stream_for_not_selected.is_open())
    return 1;

  if (write_tdt)
    stream_for_not_selected << tdt.rawdata();
  else
    stream_for_not_selected << fp.id() << '\n';

  return 1;
}

static int
handle_successful_match (const IW_TDT & tdt,
                         const IW_General_Fingerprint & fp,
                         IWString_and_File_Descriptor & output)
{
  tdts_matching++;

  if (write_tdt)
    output << tdt.rawdata();
  else if (write_smiles)
  {
    const_IWSubstring tmp;
    tdt.dataitem_value(smiles_tag, 0, tmp);
    if (tmp.starts_with(smiles_tag))
      tmp.remove_leading_chars(smiles_tag.length());
    output << tmp;
    output << ' ' << fp.id() << '\n';
  }
  else
    output << fp.id() << '\n';

  return 1;
}

static int
gfp_filter_by_bits (const IW_TDT & tdt,
                    IW_General_Fingerprint & fp,
                    IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < nq; i++)
  {
    if (! query[i].matches(fp))
    {
      if (all_queries_must_match)
        return handle_non_match(tdt, fp);
    }
    else if (! all_queries_must_match)    // first match is all we need
      return handle_successful_match(tdt, fp, output);
    else if (i == nq - 1)    // just matched last query
      return handle_successful_match(tdt, fp, output);
  }

  return handle_non_match(tdt, fp);
}

static int
gfp_filter_by_bits (IW_TDT & tdt,
                    IWString_and_File_Descriptor & output)
{
  int fatal;
  IW_General_Fingerprint fp;

  if (fp.construct_from_tdt(tdt, fatal))
    return gfp_filter_by_bits(tdt, fp, output);

  if (fatal)
  {
    cerr << "Fatal error converting from tdt to fingerprint\n";
    return 0;
  }
  else
  {
    cerr << "Non fatal error converting to fingerprint, ignored\n";
    return 1;
  }
}

static int
gfp_filter_by_bits (iwstring_data_source & input,
                    IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    tdts_read++;

    if (! gfp_filter_by_bits (tdt, output))
    {
      cerr << "Fatal error processing TDT, line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
    stream_for_not_selected.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
gfp_filter_by_bits (const char * fname,
                    IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_filter_by_bits(input, output);
}


static int
gfp_filter_by_bits (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:W:F:P:a");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('W'))
  {
    const_IWSubstring w = cl.string_value('W');

    if ("tdt" == w)
    {
      write_tdt = 1;

      if (verbose)
        cerr << "Will write tdt output\n";
    }
    else if ("smi" == w)
    {
      write_smiles = 1;

      if (verbose)
        cerr << "Will write smiles\n";
    }
    else
    {
      cerr << "Unrecognised -W qualifier '" << w << "'\n";
      return 2;
    }

  }

  if (cl.option_present('a'))
  {
    all_queries_must_match = 1;

    if (verbose)
      cerr << "For a match all query specifications must match\n";
  }

  if (! cl.option_present('s'))
  {
    cerr << "Must specify one or more selection criteria via the -s option\n";
    usage(2);
  }

  nq = cl.option_count('s');

  query = new Match_Criterion[nq];

  if (cl.option_present('s'))
  {
    const_IWSubstring s;
    for (int i = 0; i < nq; i++)
    {
      cl.value('s', s, i);

      if (! query[i].initialise(s))
      {
        cerr << "Cannot initialise bit filter query '" << s << "'\n";
        return 3;
      }
    }

    if (verbose)
      cerr << "Initialised " << nq << " bit filter queries\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! gfp_filter_by_bits(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << tdts_read << " tdts, " << tdts_matching << " matched\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_filter_by_bits(argc, argv);

  return rc;
}
