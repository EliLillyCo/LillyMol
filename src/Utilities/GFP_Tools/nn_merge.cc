/*
  We have a series of .nn files and need to merge them
*/

#include <stdlib.h>
#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/smiles_id_dist.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString distance_tag("DIST<");
static IWString identifier_tag("PCN<");
static IWString smiles_tag("$SMI<");

static IWString distance_tag_no_open_bracket("DIST");
static IWString identifier_tag_no_open_bracket("PCN");
static IWString smiles_tag_no_open_bracket("$SMI");

static int neighbours_to_keep = std::numeric_limits<int>::max();

/*
  The default output is to include the neighbour number with
  each neighbour
*/

static int write_number_neighbours_with_parent = 1;

/*
  If we are running in a mode where we have one input file run against
  a number of different chunks of the haystack, then it will be an
  error if one file doesn't contain info about each identifier
*/

static int display_missing_data_message = 1;

static float upper_distance_threshold = std::numeric_limits<float>::max();

static int brief_output_tabular_form = 0;

static Accumulator<float> shortest_distance;

static Fraction_as_String fraction_as_string;

static int needles_processed = 0;

static int filter_to_shortest_distance = 0;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Merges a series of .nn files\n";
  cerr << " -n <number>    max number of neighbours to keep\n";
  cerr << " -T <dist>      upper distance threshold\n";
  cerr << " -b             brief tabular output form of just NN distance\n";
  cerr << " -w             suppress number neighbours in output\n";
  cerr << " -q             do NOT warn about missing identifiers in files\n";
  cerr << " -u             unique neighbours only, keep shortest NN distance\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

template class resizable_array_p<Smiles_ID_Dist>;
template class resizable_array_base<Smiles_ID_Dist *>;

typedef IW_STL_Hash_Map<IWString, off_t> IW_STL_Hash_Map_off_t;

class NNFile
{
  private:
    iwstring_data_source _input;

    IW_STL_Hash_Map_off_t _offset;

  public:
    int initialise (const char *);

    unsigned int number_identifiers_stored() const { return _offset.size();}

    int tdt_for_item (const IWString &, IW_TDT & tdt);

    int next_unprocessed_identifier (IWString &);
};

int
NNFile::initialise (const char * fname)
{
  if (! _input.open(fname))
  {
    cerr << "NNFile::initialise:cannot open '" << fname << "'\n";
    return 0;
  }

  IW_TDT tdt;

  IWString pcn;   // scope here for efficiency

  off_t o = static_cast<off_t>(0);

  while (tdt.next(_input))
  {
    if (! tdt.dataitem_value(identifier_tag, pcn))
    {
      cerr << "NNFile::initialise:no tag in tdt, '" << fname << "', line " << _input.lines_read() << "\n";
      return 0;
    }

    pcn.truncate_at_last(' ');

    _offset[pcn] = o;
    o = _input.tellg();
  }

  return static_cast<int>(_offset.size());
}

int
NNFile::tdt_for_item (const IWString & id, 
                      IW_TDT & tdt)
{
  IW_STL_Hash_Map_off_t::iterator f = _offset.find(id);

  if (f == _offset.end())
    return 0;

  if (! _input.seekg ((*f).second))
  {
    cerr << "NNFile::tdt_for_item:cannot seek to " << (*f).second << " for '" << id << "'\n";
    return 0;
  }

  _offset.erase(f);

  return tdt.next(_input);
}

int
NNFile::next_unprocessed_identifier (IWString & id)
{
  if (0 == _offset.size())
    return 0;

  IW_STL_Hash_Map_off_t::iterator f = _offset.begin();

  id = (*f).first;
//id = (*(_offset.begin())).first;

  _offset.erase(f);

  return 1;
}

/*
  Get all the neighbours into a resizable array. Beware of indices because
  the target will have a smiles and id, but no distance
*/

static int
tdt_to_sid (const IW_TDT & tdt,
            resizable_array_p<Smiles_ID_Dist> & sid)
{
  IWString smiles;
  IWString id;

  const_IWSubstring ztag, zdata;

  int i = 0;
  while (tdt.next_dataitem_value (ztag, zdata, i))
  {
//  cerr << "Examining tag '" << ztag << "' data '" << zdata << " i = " << i << endl;
    if (identifier_tag_no_open_bracket == ztag)
    {
      id = zdata;
      continue;
    }

    if (smiles_tag_no_open_bracket == ztag)
    {
      smiles = zdata;
      continue;
    }

    if (distance_tag_no_open_bracket != ztag)
      continue;

    float d;

    if (! zdata.numeric_value(d) || d < 0.0f || d > 1.0f)
    {
      cerr << "Invalid distance '" << zdata << "'\n";
      return 0;
    }

    if (d > upper_distance_threshold)
      continue;

    Smiles_ID_Dist * s = new Smiles_ID_Dist(smiles, id, d);
    sid.add(s);

    smiles.resize_keep_storage(0);
    id.resize_keep_storage(0);
  }

  return 1;
}

class SID_Comparator
{
  private:
  public:
    int operator() (const Smiles_ID_Dist *, const Smiles_ID_Dist *) const;
};

int
SID_Comparator::operator() (const Smiles_ID_Dist * sid1, const Smiles_ID_Dist * sid2) const
{
  float d1 = sid1->distance();
  float d2 = sid2->distance();

  if (d1 < d2)
    return -1;

  if (d1 > d2)
    return 1;

  return 0;
}

/*
  If we are doing brief tabular output, then we just need to find
  the shortest distance
*/

static int
do_brief_output_tabular_form (const IWString & id,    // ID of needle
                              resizable_array_p<Smiles_ID_Dist> & sid,
                              IWString_and_File_Descriptor & output)
{
  int n = sid.number_elements();

  output << id;
  if (0 == n)
  {
    output << " .\n";
    return 1;
  }

// Need to find the shortest distance

  int id_with_shortest_distance = 0;
  similarity_type_t mindist = sid[0]->distance();

  for (int i = 1; i < n; i++)
  {
    const Smiles_ID_Dist * sidi = sid[i];

    if (sidi->distance() >= mindist)
      continue;

    mindist = sidi->distance();
    id_with_shortest_distance = i;
  }

  fraction_as_string.append_number(output, mindist);

  if (verbose)
    shortest_distance.extra(mindist);

  return 1;
}

static SID_Comparator sidc;

static NNFile * nfile = nullptr;
static int number_files = 0;

static int
do_filter_to_shortest_distance(resizable_array_p<Smiles_ID_Dist> & sid)
{
  IW_STL_Hash_Set seen_before;

  return sid.remove_items_fn([& seen_before] (const Smiles_ID_Dist * s) 
                                    {
                                      if (seen_before.contains(s->id()))
                                        return 1;
                                      seen_before.emplace(s->id());
                                      return 0;
                                    });
}

/*static int
do_filter_to_shortest_distance(resizable_array_p<Smiles_ID_Dist> & sid)
{
  int n = sid.number_elements();

  int rc = 0;

  IW_STL_Hash_Set seen_before;

  for (int i = 0; i < n; ++i)
  {
    const IWString & id = sid[i]->id();

    if (seen_before.contains(id))
    {
      sid.remove_item(i);
      i--;
      n--;
      rc++;
    }
    else
      seen_before.insert(id);
  }

  return rc;
}*/

static int
gather_neighbours_sort_and_write (const IW_TDT & tdt,
                                  int fstart,
                                  IWString_and_File_Descriptor & output)
{
  IWString id;
  if (! tdt.dataitem_value(identifier_tag, id))
  {
    cerr << "gather_neighbours_sort_and_write:tdt without identifier\n";
    cerr << "TDT has " << tdt.number_elements() << " items\n";
    cerr << tdt << '\n';
    return 0;
  }

  id.truncate_at_last(' ');

  IWString smiles;

  (void) tdt.dataitem_value(smiles_tag, smiles);

  resizable_array_p<Smiles_ID_Dist> sid;

  if (! tdt_to_sid (tdt, sid))
  {
    cerr << "Cannot parse tdt\n";
    return 0;
  }

  for (int i = fstart; i < number_files; i++)
  {
    IW_TDT tdt;
    if (! nfile[i].tdt_for_item(id, tdt))
    {
      if (display_missing_data_message)
        cerr << "No data for '" << id << "' in file " << i << endl;
      continue;
    }

    if (! tdt_to_sid(tdt, sid))
    {
      cerr << "Cannot parse TDT from file " << i << " for '" << id << "'\n";
      return 0;
    }
  }

  if (verbose > 2)
    cerr << "Got " << sid.number_elements() << " nbrs for '" << id << "'\n";

  if (brief_output_tabular_form)   // no need to sort
    return do_brief_output_tabular_form(id, sid, output);

  sid.iwqsort(sidc);

  if (filter_to_shortest_distance)
    do_filter_to_shortest_distance(sid);

  int nwrite = sid.number_elements();
  if (nwrite > neighbours_to_keep)
    nwrite = neighbours_to_keep;

  if (verbose && nwrite > 0)
    shortest_distance.extra(sid[0]->distance());

  output << smiles_tag << smiles << ">\n";
  output << identifier_tag << id;
  if (write_number_neighbours_with_parent)
    output << ' ' << nwrite;
  output << ">\n";

  for (int i = 0; i < nwrite; i++)
  {
    const Smiles_ID_Dist * sidi = sid[i];

    output << smiles_tag << sidi->smiles() << ">\n";
    output << identifier_tag << sidi->id() << ">\n";
    fraction_as_string.append_number(output, sidi->distance());

    output.write_if_buffer_holds_more_than(8192);
  }

  output << "|\n";

  return 1;
}

static int
nn_merge (iwstring_data_source & input,
          NNFile * nfile,
          const int number_files,
          IWString_and_File_Descriptor & output)
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    if (0 == tdt.number_elements())
    {
      cerr << "Skipping empty tdt\n";
      continue;
    }

    needles_processed++;

    if (! gather_neighbours_sort_and_write(tdt, 1, output))   // 1 means start with file 1
      return 0;
  }

// At this stage, we have written all nn lists that were in the first file.

  IWString id;  // scope here for efficiency

  for (int i = 1; i < number_files; i++)
  {
    while (nfile[i].next_unprocessed_identifier(id))
    {
      if (! nfile[i].tdt_for_item(id, tdt))
      {
        cerr << "Huh, file " << i << " cannot return inprocessed identifier '" << id << "'\n";
        return 0;
      }

      needles_processed++;

      if (! gather_neighbours_sort_and_write(tdt, i + 1, output))
        return 0;

      output.write_if_buffer_holds_more_than(32768);
    }
  }

  return 1;
}

static int
nn_merge (const char * fname,
          NNFile * nfile,
          const int number_files,
          IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return nn_merge (input, nfile, number_files, output);
}

static int
nn_merge (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:qT:wbu");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    if (! cl.value('n', neighbours_to_keep) || neighbours_to_keep < 1)
    {
      cerr << "The neighbours to keep option (-n) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will keep a maximum of " << neighbours_to_keep << " neighbours\n";
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', upper_distance_threshold) || upper_distance_threshold < 0.0)
    {
      cerr << "The upper distance threshold (-T) option must be a valid distance\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard neighbours with distances beyond " << upper_distance_threshold << endl;
  }

  if (cl.option_present('q'))
  {
    display_missing_data_message = 1;

    if (verbose)
      cerr << "No messages about missing identifiers\n";
  }

  if (cl.option_present('b'))
  {
    brief_output_tabular_form = 1;

    if (verbose)
      cerr << "Output is just ID and NNDist\n";
  }

  if (cl.option_present('u'))
  {
    filter_to_shortest_distance = 1;

    if (verbose)
      cerr << "Will discard neighbours other than the shortest\n";
  }

  if (cl.option_present('w'))
  {
    write_number_neighbours_with_parent = 0;

    if (verbose)
      cerr << "Neighbour numbers not written\n";
  }

  number_files = cl.number_elements();

  if (0 == number_files)
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (1 == number_files)
  {
    cerr << "programme merges multiple .nn files, really should specify more than 1 input\n";
  }

  nfile = new NNFile[number_files]; std::unique_ptr<NNFile[]> free_nfile(nfile);

  for (int i = 1; i < number_files; i++)
  {
    if (! nfile[i].initialise(cl[i]))
    {
      cerr << "Cannot read neighbour info from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  unsigned int items_per_file = nfile[1].number_identifiers_stored();

  if (verbose)
    cerr << "Second file contains " << items_per_file << " items\n";

  for (int i = 2; i < number_files; i++)
  {
    unsigned int j = nfile[i].number_identifiers_stored();

    if (j == items_per_file)
      continue;

    cerr << "Possible size mismatch, first file contains " << items_per_file << " items, file " << i << " contains " << j << endl;
  }

  if (brief_output_tabular_form)
    fraction_as_string.set_leading_string(" ");
  else
    fraction_as_string.set_leading_string(distance_tag);

  fraction_as_string.initialise(0.0, 1.0, 3);

  if (brief_output_tabular_form)
    fraction_as_string.append_to_each_stored_string("\n");
  else
    fraction_as_string.append_to_each_stored_string(">\n");
  
  IWString_and_File_Descriptor output(1);

  if (brief_output_tabular_form)
    output << "ID NNDist\n";

  if (! nn_merge (cl[0], nfile, number_files, output))
  {
    cerr << "Merging failed\n";
    return 3;
  }

  output.flush();

  if (verbose)
  {
    cerr << "Processed " << cl.number_elements() << " files\n";
    if (shortest_distance.n() > 1)
      cerr << "Shortest distances between " << shortest_distance.minval() << " and " << shortest_distance.maxval() << " ave " << shortest_distance.average_if_available_minval_if_not() << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = nn_merge(argc, argv);

  return rc;
}
