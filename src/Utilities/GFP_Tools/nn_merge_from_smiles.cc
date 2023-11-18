/*
  Merges a series of near neighbour files. These must have been processed by nplotnn
*/

#include <stdlib.h>
#include <iostream>
#include <limits>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Utilities/GFP_Tools/smiles_id_dist.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static float upper_distance_threshold = std::numeric_limits<float>::max();

static int nbrs_to_write = std::numeric_limits<int>::max();

static int needles_processed = 0;

static int ignore_self_neighbours = 0;

static extending_resizable_array<int> nbr_count;

static Fraction_as_String fraction_as_string;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Merges near neighbour smiles files that have been processed by nplotnn\n";
  cerr << " -n <number>    max number of neighbours to keep\n";
  cerr << " -T <dist>      upper distance threshold\n";
  cerr << " -h             ignore self neighbours\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
write_nbr_list(const IWString & smiles0,
               const IWString & id0,
               const resizable_array_p<Smiles_ID_Dist> & nbrs,
               IWString_and_File_Descriptor & output)
{
  int n = nbrs.number_elements();

  if (n > nbrs_to_write)
    n = nbrs_to_write;

  output << smiles0 << ' ' << id0 << ' ' << n << "\n";

  for (int i = 0; i < n; i++)
  {
    const Smiles_ID_Dist * sidi = nbrs[i];

    output << sidi->smiles() << ' ' << sidi->id();
    fraction_as_string.append_number(output, sidi->distance());

    output.write_if_buffer_holds_more_than(8192);
  }

  output.flush();

  nbr_count[n]++;

  return 1;
}

class Distance_Comparator
{
  private:
  public:
    int operator () (const Smiles_ID_Dist *, const Smiles_ID_Dist *) const;
};

int
Distance_Comparator::operator() (const Smiles_ID_Dist * sid1, const Smiles_ID_Dist * sid2) const
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
  smiles id <skip> distanc3
*/

template <typename T>
int
parse_smiles_id_dist (const const_IWSubstring & buffer,
                      IWString & smiles,
                      IWString & id,
                      T & x)
{
  int i = 0;

  buffer.nextword(smiles, i);
  buffer.nextword(id, i);

  const_IWSubstring s;

  buffer.nextword(s, i);
  buffer.nextword(s, i);

  if (! s.numeric_value(x))
  {
    cerr << "Invalid distance '" << s << "'\n";
    return 0;
  }

  return 1;
}

/*
  smiles id and number nbrs
*/

template <typename T>
int
parse_smiles_id_numeric (const const_IWSubstring & buffer,
                         IWString & smiles,
                         IWString & id,
                         T & x)
{
  int i = 0;

  buffer.nextword(smiles, i);
  buffer.nextword(id, i);

  const_IWSubstring s;

  buffer.nextword(s, i);

  if (! s.numeric_value(x))
  {
    cerr << "Invalid numeric '" << s << "'\n";
    return 0;
  }

  return 1;
}

static int
read_nbrs (iwstring_data_source & input,
           const IWString & id0,
           int number_nbrs,
           resizable_array_p<Smiles_ID_Dist> & nbrs,
           int & fatal)
{
  nbrs.make_room_for_extra_items(number_nbrs);

//cerr << "Requiest to read " << number_nbrs << " nbrs\n";

  const_IWSubstring buffer;

  for (int i = 0; i < number_nbrs; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Premature eof, expected " << number_nbrs << " nbrs, only read " << i << endl;
      return 0;
    }

    IWString smiles, id;     // yes these shadow a function parameter
    float d;

    if (! parse_smiles_id_dist(buffer, smiles, id, d) || d < 0.0f)
    {
      cerr << "Invalid nbr specification '" << buffer << "'\n";
      fatal = 1;
      return 0;
    }

    if (0.0f == d && ignore_self_neighbours && id == id0)
      continue;

    if (d <= upper_distance_threshold)
      nbrs.add(new Smiles_ID_Dist(smiles, id, d));
  }

  return 1;
}

static int
get_next_needle (iwstring_data_source & input,
                 const IWString & id,
                 int & fatal,
                 resizable_array_p<Smiles_ID_Dist> & nbrs)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    fatal = 0;
    return 0;
  }

  int number_nbrs;

  IWString mysmiles, my_id;

  if (! parse_smiles_id_numeric (buffer, mysmiles, my_id, number_nbrs) || number_nbrs < 0)
  {
    cerr << "get_next_needle:invalid input '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  if (my_id != id)
  {
    cerr << "Identifier mismatch, looking for nbrs of '" << id << "' got '" << my_id << "'\n";
    fatal = 1;
    return 0;
  }

  return read_nbrs (input, my_id, number_nbrs, nbrs, fatal);
}

static int
get_next_needle0 (iwstring_data_source & input,
                  IWString & smiles,
                  IWString & id,
                  int & fatal,
                  resizable_array_p<Smiles_ID_Dist> & nbrs)
{
  nbrs.resize_keep_storage(0);

  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    fatal = 0;
    return 0;
  }

//cerr << "First record read '" << buffer << "'\n";

  int number_nbrs;

  if (! parse_smiles_id_numeric (buffer, smiles, id, number_nbrs) || number_nbrs < 0)
  {
    cerr << "get_next_needle:invalid input '" << buffer << "'\n";
    return 0;
  }

  return read_nbrs(input, id, number_nbrs, nbrs, fatal);
}

static int
nn_merge_from_smiles (iwstring_data_source * input,
                      int nfiles,
                      IWString_and_File_Descriptor & output)
{
  Distance_Comparator dc;

  IWString smiles0, id0;

  resizable_array_p<Smiles_ID_Dist> nbrs;

  int fatal;

  while (get_next_needle0 (input[0], smiles0, id0, fatal, nbrs))
  {
//  cerr << "Processing '" << id0 << "'\n";
    needles_processed++;

    for (int i = 1; i < nfiles; i++)
    {
      if (! get_next_needle(input[i], id0, fatal, nbrs))
      {
        cerr << "Cannot fetch nbrs for '" << id0 << "' from file " << i << endl;
        return 0;
      }
    }

    if (nbrs.number_elements() > 1)
      nbrs.iwqsort(dc);

    if (! write_nbr_list(smiles0, id0, nbrs, output))
      return 0;
  }

  if (fatal)
  {
    cerr << "Error exit\n";
    return 0;
  }

  return 1;
}

static int
nn_merge_from_smiles (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vT:n:h");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('T'))
  {
    if (! cl.value('T', upper_distance_threshold) || upper_distance_threshold < 0.0f || upper_distance_threshold > 1.0f)
    {
      cerr << "The upper distance threshold must be a valid distance\n";
      usage(2);
    }

    if (verbose)
      cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', nbrs_to_write) || nbrs_to_write < 0)
    {
      cerr << "The number of nbrs to write (-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will write a max of " << nbrs_to_write << " nbrs for each molecule\n";
  }

  if (cl.option_present('h'))
  {
    ignore_self_neighbours = 1;

    if (verbose)
      cerr << "Will ignore self neighbours\n";
  }

  int nfiles = cl.number_elements();

  if (0 == nfiles)
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (1 == nfiles)
  {
    cerr << "Only works with multiple files\n";
    usage(2);
  }

  iwstring_data_source * input = new iwstring_data_source[nfiles];

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! input[i].open (cl[i]))
    {
      cerr << "Cannot open '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  fraction_as_string.set_leading_string(" ");

  fraction_as_string.initialise(0.0, 1.0, 3);

  fraction_as_string.append_to_each_stored_string("\n");

  IWString_and_File_Descriptor output(1);

  if (! nn_merge_from_smiles (input, nfiles, output))
  {
    cerr << "Cannot merge nn lists\n";
    return 2;
  }

  if (verbose)
  {
    cerr << "Processed nbrs for " << needles_processed << " needles\n";
    for (int i = 0; i < nbr_count.number_elements(); i++)
    {
      if (nbr_count[i])
        cerr << nbr_count[i] << " molecules had " << i << " neighbours\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = nn_merge_from_smiles(argc, argv);

  return rc;
}
