/*
  Produces a .nn file from a distance matrix
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/iwdigits.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "IWDistanceMatrixBase.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

int neighbours_to_find = 0;

static int molecules_processed = 0;

static float lower_distance_threshold = static_cast<float> (-1.0);
static float upper_distance_threshold = static_cast<float> (0.0);

static IW_STL_Hash_Map_String smiles;

static IWString default_smiles;

static int write_molecules_with_no_neighbours = 0;

static int molecules_with_no_neighbours = 0;

static Fraction_as_String fraction_as_string;

static int fast_io = 0;

static int write_neighbours_as_indices = 0;

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");
static IWString nbr_tag("NBR<");

/*
  Various summary statistics
*/

static Accumulator<float> closest_neighbour;
static Accumulator<float> distance_stats;
static Accumulator<float> number_neighbours;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id: distance_matrix_nn.cc,v 1.1 2003/05/28 17:05:57 ian Exp $\n";
  cerr << "Scans a distance matrix file and produces output for nplotnn\n";
  cerr << " -n <number>    number of neighbours for each molecule\n";
  cerr << " -t <dist>      lower distance threshold\n";
  cerr << " -T <dist>      upper distance threshold\n";
  cerr << " -D <fname>     the distance matrix file\n";
  cerr << " -S <fname>     smiles file for molecules in distance matrix\n";
  cerr << " -S default=XX  no smiles data available, use XX as smiles for all id's\n";
  cerr << " -o             write neighbours as indices\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
read_smiles_2 (const const_IWSubstring & buffer,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring smi, id;

  if (! buffer.split (smi, ' ', id))
  {
    cerr << "Cannot extract smiles and ID\n";
    return 0;
  }

  id.truncate_at_first (' ');

  smiles[id] = smi;

  return 1;
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! read_smiles_2 (buffer, smiles))
      return 0;
  }

  return smiles.size ();
}

static int
read_smiles (const const_IWSubstring & fname,
             IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles (input, smiles);
}

static int
write_smiles_and_id (const IWString & id,
                     int ndx,
                     const IW_STL_Hash_Map_String & smiles,
                     IWString & output)
{
  if (smiles.size ())
  {
    output << smiles_tag;

    IW_STL_Hash_Map_String::const_iterator f = smiles.find (id);
    if (f == smiles.end ())
    {
      cerr << "No smiles data for '" << id << "'\n";
      output << "C";
    }
    else
      output << (*f).second;
    output << ">\n";
  }
  else if (default_smiles.length () > 0)
    output << default_smiles;

  if (write_neighbours_as_indices && ndx >= 0)
    output << nbr_tag << ndx << ">\n";
  else
    output << identifier_tag << id << ">\n";

  return 1;
}

template <typename T, typename C>
int
distance_matrix_nn (IWDistanceMatrixBase<T> & dm,
                    ID_Distance<T> * idd,
                    C & comparitor,
                    int t,             // the target molecule
                    IWString & output_buffer)
{
  molecules_processed++;

  int neighbours_written = 0;

  int n = dm.number_molecules ();

  for (int i = 0; i < n; i++)
  {
    if (i == t)     // no self neighbours
      continue;

    T d = dm.zvalue (t, i);

    if (lower_distance_threshold >= static_cast<float> (0.0) && d <= static_cast<T> (lower_distance_threshold))
      continue;

    if (upper_distance_threshold > static_cast<float> (0.0) && d > static_cast<T> (upper_distance_threshold))
      continue;

    if (verbose)
      distance_stats.extra (static_cast<float> (d));

    idd[neighbours_written].set_id (i);
    idd[neighbours_written].set_distance (d);
    neighbours_written++;
  }

  if (verbose)
    number_neighbours.extra (neighbours_written);

  if (0 == neighbours_written && ! write_molecules_with_no_neighbours)
  {
    molecules_with_no_neighbours++;
    return 1;
  }

  if (neighbours_written > 1)
    iwqsort (idd, neighbours_written, comparitor);

  if (verbose)
    closest_neighbour.extra (static_cast<float> (idd[0].distance ()));

  const IWString & id = dm.id (t);

  write_smiles_and_id (id, -1, smiles, output_buffer);

  if (neighbours_to_find > 0 && neighbours_written > neighbours_to_find)
    neighbours_written = neighbours_to_find;

  for (int i = 0; i < neighbours_written; i++)
  {
    const ID_Distance<T> & iddi = idd[i];

    int j = iddi.id ();

    write_smiles_and_id (dm.id (j), j, smiles, output_buffer);
  
    output_buffer << distance_tag;

    if (fast_io)
      fraction_as_string.append_number (output_buffer, iddi.distance ());
    else
      output_buffer << iddi.distance ();

    output_buffer << ">\n";
  }

  output_buffer << "|\n";

  return 1;
}

template <typename T, typename C>
int
distance_matrix_nn (IWDistanceMatrixBase<T> & dm,
                    ID_Distance<T> * idd,
                    C & comparitor,
                    int output_fd)
{
  int n = dm.number_molecules ();

  if (verbose)
    cerr << "Processing " << n << " molecules\n";

  IWString output_buffer;

  for (int i = 0; i < n; i++)
  {
    distance_matrix_nn (dm, idd, comparitor, i, output_buffer);

    if (output_buffer.length() > 32768)
      output_buffer.write_whole_blocks_shift_unwritten(output_fd);
  }

  if (output_buffer.length())
    output_buffer.write(output_fd);

  return 1;
}

template <typename T, typename C>
int
distance_matrix_nn_2 (const const_IWSubstring & buffer,
                      IWDistanceMatrixBase<T> & dm,
                      ID_Distance<T> * idd,
                      C & comparitor,
                      IWString & output_buffer)
{
  int i = dm.which_item_has_id (buffer);

  if (i < 0)
  {
    cerr << "Distance matrix has no data on '" << buffer << "'\n";
    return 0;
  }

  return distance_matrix_nn (dm, idd, comparitor, i, output_buffer);
}

template <typename T, typename C>
int
distance_matrix_nn (iwstring_data_source & input,
                    IWDistanceMatrixBase<T> & dm,
                    ID_Distance<T> * idd,
                    C & comparitor,
                    int output_fd)
{
  const_IWSubstring buffer;

  IWString output_buffer;

  while (input.next_record (buffer))
  {
    buffer.truncate_at_first (' ');

    if (! distance_matrix_nn_2 (buffer, dm, idd, comparitor, output_buffer))
    {
      cerr << "Failure processing '" << buffer << "'\n";
      return 0;
    }

    if (output_buffer.length() > 32768)
      output_buffer.write_whole_blocks_shift_unwritten(output_fd);
  }

  if (output_buffer.length())
    output_buffer.write(output_fd);

  return 1;
}

template <typename T, typename C>
int
distance_matrix_nn (const const_IWSubstring & fname,
                    IWDistanceMatrixBase<T> & dm,
                    ID_Distance<T> * idd,
                    C & comparitor,
                    int output_fd)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return distance_matrix_nn (input, dm, idd, comparitor, output_fd);
}

template <typename T, typename C>
int
distance_matrix_nn (const Command_Line & cl,
                    IWDistanceMatrixBase<T> & dm,
                    ID_Distance<T> * idd,
                    C & comparitor,
                    int output_fd)
{
  if (! cl.option_present ('D'))
    return distance_matrix_nn (dm, idd, comparitor, output_fd);

  int i = 0;
  const_IWSubstring d;
  while (cl.value ('D', d, i++))
  {
    if (! distance_matrix_nn (d, dm, idd, comparitor, output_fd))
      return 0;
  }

  return 1;
}

template <typename T, typename C>
int
distance_matrix_nn (const Command_Line & cl,
                    IWDistanceMatrixBase<T> & dm,
                    const const_IWSubstring & dm_name,
                    C & comparitor,
                    int output_fd)
{
  cerr << "Opening distance matrix '" << dm_name << "'\n";
  if (! dm.do_read (dm_name))
  {
    cerr << "Cannot initialise distance matrix '" << dm_name << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Distance matrix '" << dm_name << "' contains data on " << dm.number_molecules () << " molecules\n";

  if (fast_io > 0)
  {
    T m = dm.maxval ();

    if (! fraction_as_string.initialise (static_cast<float> (0.0), static_cast<float> (m), fast_io))
      return 0;
  }

  ID_Distance<T> * idd = new ID_Distance<T>[dm.number_molecules ()];

  if (NULL == idd)
  {
    cerr << "Cannot allocate " << dm.number_molecules () << " comparison objects\n";
    return 0;
  }

  if (neighbours_to_find > dm.number_molecules () - 1)
  {
    cerr << "Requested " << neighbours_to_find << " but only " << dm.number_molecules () << " molecules in matrix\n";
    neighbours_to_find = dm.number_molecules () - 1;
  }

  int rc = distance_matrix_nn (cl, dm, idd, comparitor, output_fd);

  delete [] idd;

  return rc;
}

template <typename T>
class ID_Distance_Comparitor
{
  private:
  public:
    int operator () (const ID_Distance<T> &, const ID_Distance<T> &);
};

template <typename T>
int
ID_Distance_Comparitor<T>::operator () (const ID_Distance<T> & idd1,
                                        const ID_Distance<T> & idd2)
{
  T d1 = idd1.distance ();
  T d2 = idd2.distance ();

  if (d1 < d2)
    return -1;

  if (d1 > d2)
    return 1;

  return 0;
}

static int
distance_matrix_nn (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:t:T:D:k:S:f:o");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present ('n') && ! cl.option_present ('t') && ! cl.option_present ('T'))
  {
    neighbours_to_find = 1;
  }

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', neighbours_to_find) || neighbours_to_find < 1)
    {
      cerr << "The number of neighbours to find (-n) must be a whole positive number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will find " << neighbours_to_find << " neighbours for each item\n";
  }

  if (cl.option_present('o'))
  {
    write_neighbours_as_indices = 1;
    if (verbose)
      cerr << "Neighbours written as indices\n";
  }

  if (cl.option_present ('t'))
  {
    if (! cl.value ('t', lower_distance_threshold) || lower_distance_threshold < static_cast<float> (0.0))
    {
      cerr << "The lower distance threshold must be non-negative\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will discard distances of " << lower_distance_threshold << " or shorter\n";
  }

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', upper_distance_threshold) || upper_distance_threshold < lower_distance_threshold)
    {
      cerr << "The upper distance threshold must be greater than " << lower_distance_threshold << endl;
      usage (5);
    }

    if (verbose)
      cerr << "Will discard distances larger than " << upper_distance_threshold << endl;
  }

  if (cl.option_present ('S'))
  {
    const_IWSubstring s = cl.string_value ('S');

    if (s.starts_with ("default="))
    {
      s.remove_leading_chars (8);
      default_smiles << smiles_tag << s << ">\n";
    }
    else if (! read_smiles (s, smiles))
    {
      cerr << "Cannot read smiles data '" << s << "'\n";
      return 3;
    }
    else if (verbose)
      cerr << "Read " << smiles.size () << " id->smiles relationships from '" << s << "'\n";
  }

  if (cl.option_present ('f'))
  {
    if (! cl.value ('f', fast_io) || fast_io < 2)
    {
      cerr << "The fast I/O option needs a precision of > 1 digit\n";
      usage (6);
    }

    if (verbose)
      cerr << "Will use fast I/O, " << fast_io << " decimal digits\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Must specify the distance matrix on the command line\n";
    usage (5);
  }

  const_IWSubstring d = cl[0];

  const_IWSubstring k;
  if (cl.option_present ('k'))
    k = cl.string_value ('k');
  else
    k = "float";

  int rc;

  if ("float" == k)
  {
    ID_Distance_Comparitor<float> comparitor;
    IWDistanceMatrixBase<float> dm;
    rc = distance_matrix_nn (cl, dm, d, comparitor, 1);
  }
  else if ("double" == k)
  {
    ID_Distance_Comparitor<double> comparitor;
    IWDistanceMatrixBase<double> dm;
    rc = distance_matrix_nn (cl, dm, d, comparitor, 1);
  }
  else if ("int" == k)
  {
    IWDistanceMatrixBase<int> dm;
    ID_Distance_Comparitor<int> comparitor;
    rc = distance_matrix_nn (cl, dm, d, comparitor, 1);
  }
  else if ("uchar" == k)
  {
    IWDistanceMatrixBase<unsigned char> dm;
    ID_Distance_Comparitor<unsigned char> comparitor;
    rc = distance_matrix_nn (cl, dm, d, comparitor, 1);
  }
  else
  {
    cerr << "Not sure what to do with distance matrix type '" << k << "'\n";
    return 8;
  }

  if (0 == rc)
    return 5;

  if (verbose)
  {
    cerr << "Processed " << molecules_processed << " molecules\n";
    cerr << "Had between " << number_neighbours.minval () << " and " << number_neighbours.maxval () << " neighbours";
    if (number_neighbours.n () > 1 && number_neighbours.minval () != number_neighbours.maxval ())
      cerr << " ave " << number_neighbours.average ();
    cerr << endl;

    cerr << "Shortest distances between " << closest_neighbour.minval () << " and " << closest_neighbour.maxval ();
    if (closest_neighbour.n () > 1)
      cerr << " ave " << closest_neighbour.average ();
    cerr << endl;
     
    cerr << "All distances between " << distance_stats.minval () << " and " << distance_stats.maxval ();
    if (distance_stats.n () > 1)
      cerr << " ave " << distance_stats.average ();
    cerr << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_nn (argc, argv);

  return rc;
}
