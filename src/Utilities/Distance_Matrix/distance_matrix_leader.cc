/*
  Produces a .ldr file from a distance matrix
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

using std::cerr;
using std::endl;
using std::cout;

const char * prog_name = NULL;

static int verbose = 0;

static int number_molecules = 0;

static IW_STL_Hash_Map_String smiles;

static Fraction_as_String fraction_as_string;

static int fast_io = 0;

/*
  Enable a per-molecule max cluster size and threshold
*/

static int * max_cluster_size = NULL;

static int global_max_cluster_size = 0;

static float global_threshold;

static float * threshold = NULL;

static int threshold_column = -1;

static IWString threshold_file_name;

static int * already_selected = NULL;

static extending_resizable_array<int> cluster_size;

/*
  We may have preferences for what gets selected first
*/

static float * score = NULL;

static Accumulator<float> distance_stats;

static int max_clusters_to_find = 0;

static int clusters_found = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static IWString default_smiles;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id: distance_matrix_leader.cc,v 1.1 2003/05/28 17:06:05 ian Exp $\n";
  cerr << "Does a leader computation on a distance matrix\n";
  cerr << " -t <dist>      distance threshold\n";
  cerr << " -m <number>    maximum cluster size\n";
  cerr << " -S <fname>     smiles file for molecules in distance matrix\n";
  cerr << " -S default=XX  no smiles data available, use XX as smiles for all id's\n";
  cerr << " -H col=<col>   threshold for each item in column <col> of name\n";
  cerr << " -H <fname>     thresholds in file <fname>\n";
//cerr << " -k <type>      type of distance matrix\n";
  cerr << " -f <digits>    fast I/O, use <digits> precision\n";
  cerr << " -s <fname>     score or selection order for each molecule in <fname>\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
read_smiles_2(const const_IWSubstring & buffer,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring smi, id;

  if (! buffer.split(smi, ' ', id))
  {
    cerr << "Cannot extract smiles and ID\n";
    return 0;
  }

  id.truncate_at_first(' ');

  smiles[id] = smi;

  return 1;
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! read_smiles_2(buffer, smiles))
      return 0;
  }

  return smiles.size();
}

static int
read_smiles (const const_IWSubstring & fname,
             IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles(input, smiles);
}

static int
write_smiles_and_id (const IWString & id,
                     const IW_STL_Hash_Map_String & smiles,
                     std::ostream & output)
{
  if (smiles.size())
  {
    output << smiles_tag;

    IW_STL_Hash_Map_String::const_iterator f = smiles.find(id);
    if (f == smiles.end())
    {
      cerr << "No smiles data for '" << id << "'\n";
      output << "C";
    }
    else
      output << (*f).second;
    output << ">\n";
  }
  else if (default_smiles.length() > 0)
    output << default_smiles;

  output << identifier_tag << id << ">\n";

  return output.good();
}

template <typename T, typename C>
int
distance_matrix_leader (const IWDistanceMatrixBase<T> & dm,
                        ID_Distance<T> * idd,
                        C & comparitor,
                        int leader,
                        T threshold,
                        int max_cluster_size,
                        std::ostream & output)
{
  int items_this_cluster = 0;

  int n = dm.number_molecules();

  for (int i = 0; i < n; i++)
  {
    if (already_selected[i])
      continue;

    T d = dm.zvalue(leader, i);

    if (d > threshold)
      continue;

    if (verbose)
      distance_stats.extra(static_cast<float>(d));

    idd[items_this_cluster].set_id(i);
    idd[items_this_cluster].set_distance(d);
    items_this_cluster++;
  }

  if (items_this_cluster > 1) 
    iwqsort(idd, items_this_cluster, comparitor);

  if (max_cluster_size > 0 && (items_this_cluster + 1) > max_cluster_size)
    items_this_cluster = max_cluster_size - 1;

  if (verbose)
    cluster_size[items_this_cluster + 1]++;

  const IWString & id = dm.id(leader);
  write_smiles_and_id(id, smiles, output);

  output << "CLUSTER<" << clusters_found << ">\n";
  output << "CSIZE<" << items_this_cluster << ">\n";

  IWString output_buffer(distance_tag);

  for (int i = 0; i < items_this_cluster; i++)
  {
    const ID_Distance<T> & iddi = idd[i];

    int j = iddi.id();

    already_selected[j] = 1;

    write_smiles_and_id(dm.id(j), smiles, output);
  
    output_buffer.resize_keep_storage(distance_tag.length());
    if (fast_io)
      fraction_as_string.append_number(output_buffer, iddi.distance());
    else
      output_buffer << iddi.distance ();

    output_buffer << ">\n";

    output << output_buffer;
  }

  output << "|\n";

  return output.good();
}

static int
find_next_leader (const int * already_selected,
                  int n)
{
  if (NULL == score)
    return locate_item_in_array(0, n, already_selected);

  int rc = -1;
  float highest_score = score[0];

  for (int i = 0; i < n; i++)
  {
    if (already_selected[i])
      continue;

    if (rc < 0 || score[i] > highest_score)
    {
      highest_score = score[i];
      rc = i;
    }
  }

  return rc;
}

template <typename T, typename C>
int
distance_matrix_leader (const IWDistanceMatrixBase<T> & dm,
                    ID_Distance<T> * idd,
                    C & comparitor,
                    std::ostream & output)
{
  int n = dm.number_molecules();

  if (max_clusters_to_find <= 0)
    max_clusters_to_find = n;

  while (clusters_found < max_clusters_to_find)
  {
    int leader = find_next_leader(already_selected, n);

    if (leader < 0)
      break;

    already_selected[leader] = 1;

    if (! distance_matrix_leader(dm, idd, comparitor, leader, static_cast<T>(threshold[leader]), max_cluster_size[leader], output))
      break;

    clusters_found++;
  }

  return output.good();
}

/*
  If there is no scoring data in the file, we want to infer highest score at
  the top of the file. So, we have a global counter of next score to be allocated
*/

static int next_default_score_to_be_allocated = 100000;   // arbitrary large number

template <typename T, typename D>
int
read_data_2 (const const_IWSubstring & buffer,
             const IWDistanceMatrixBase<T> & dm,
             D * v)
{
  const_IWSubstring id, string_value;

  D numeric_value;

  if (buffer.split(id, ' ', string_value))
  {
    if (! string_value.numeric_value(numeric_value))
    {
      cerr << "Invalid score '" << buffer << "'\n";
      return 0;
    }
  }
  else
  {
    numeric_value = static_cast<D>(next_default_score_to_be_allocated);
    next_default_score_to_be_allocated--;
    id = buffer;
  }

  int i = dm.which_item_has_id(id);

  if (i < 0)
  {
    cerr << "Distance matrix has no data on '" << id << "'\n";
    return 0;
  }

  v[i] = numeric_value;

  return 1;
}

template <typename T, typename D>
int
read_data (iwstring_data_source & input,
           const IWDistanceMatrixBase<T> & dm,
           D * v)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! read_data_2(buffer, dm, v))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

template <typename T, typename D>
int
read_data (const const_IWSubstring & fname,
           const IWDistanceMatrixBase<T> & dm,
           D * v)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return read_data(input, dm, v);
}

template <typename T>
int
read_scores (const const_IWSubstring & fname,
             const IWDistanceMatrixBase<T> & dm)
{
  score = new_float(dm.number_molecules(), static_cast<float>(0.0));

  if (NULL == score)
  {
    cerr << "Cannot allocate " << dm.number_molecules() << " float items\n";
    return 0;
  }

  return read_data(fname, dm, score);
}

template <typename T>
int
determine_thresholds (IWDistanceMatrixBase<T> & dm,
                      iwstring_data_source & input,
                      float * threshold)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    const_IWSubstring id, d;

    if (! buffer.split(id, ' ', d))
    {
      cerr << "The threshold file must contain two tokens, id and distance\n";
      return 0;
    }

    int j = dm.which_item_has_id(id);

    if (j < 0)
    {
      cerr << "Yipes, no information on '" << id << "' in distance matrix\n";
      return 0;
    }

    if (! d.numeric_value(threshold[j]) || threshold[j] < static_cast<float>(0.0))
    {
      cerr << "Invalid threshold '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

template <typename T>
int
determine_thresholds (IWDistanceMatrixBase<T> & dm,
                      const IWString & fname,
                      float * threshold)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open threshold file '" << fname << "'\n";
    return 0;
  }

  return determine_thresholds(dm, input, threshold);
}

template <typename T>
int
determine_thresholds (IWDistanceMatrixBase<T> & dm,
                      int threshold_column,
                      float * threshold)
{
  int n = dm.number_molecules();

  for (int i = 0; i < n; i++)
  {
    const IWString & id = dm.id(i);

    const_IWSubstring token;

    if (! id.word(threshold_column, token))
    {
      cerr << "Cannot extract threshold column from '" << id << "'\n";
      return 0;
    }

    if (! token.numeric_value(threshold[i]) || threshold[i] <= static_cast<float>(0.0))
    {
      cerr << "Invalid threshold '" << token << "' in '" << id << "'\n";
      return 0;
    }
  }

  return 1;
}

template <typename T, typename C>
int
distance_matrix_leader (Command_Line & cl,
                        IWDistanceMatrixBase<T> & dm,
                        const const_IWSubstring & dm_name,
                        C & comparitor,
                        std::ostream & output)
{
  if (! dm.do_read(dm_name))
  {
    cerr << "Cannot initialise distance matrix '" << dm_name << "'\n";
    return 0;
  }

  int n = dm.number_molecules();

  number_molecules = n;

  if (verbose)
    cerr << "Distance matrix '" << dm_name << "' contains data on " << number_molecules << " molecules\n";

  if (fast_io > 0)
  {
    T m = dm.maxval();

    if (! fraction_as_string.initialise(static_cast<float>(0.0), static_cast<float>(m), fast_io))
      return 0;
  }

  if (threshold_column >= 0)
  {
    if (! determine_thresholds(dm, threshold_column, threshold))
      return 0;
  }
  else if (threshold_file_name.length() > 0)
  {
    if (! determine_thresholds(dm, threshold_file_name, threshold))
      return 0;
  }

  ID_Distance<T> * idd = new ID_Distance<T>[number_molecules];

  already_selected = new_int(number_molecules); std::unique_ptr<int[]> free_already_selected(already_selected);

  threshold = new_float(number_molecules, global_threshold); std::unique_ptr<float[]> free_threshold(threshold);

  max_cluster_size = new_int(number_molecules, global_max_cluster_size); std::unique_ptr<int[]> free_max_cluster_size(max_cluster_size);

  if (NULL == idd || NULL == already_selected || NULL == threshold || NULL == max_cluster_size)
  {
    cerr << "Cannot allocate " << dm.number_molecules() << " objects\n";
    return 0;
  }

  if (cl.option_present('s'))
  {
    const_IWSubstring s = cl.string_value('s');

    if (! read_scores(s, dm))
    {
      cerr << "Cannot read scoring data from '" << s << "'\n";
      return 0;
    }
  }

  int rc = distance_matrix_leader(dm, idd, comparitor, output);

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
  T d1 = idd1.distance();
  T d2 = idd2.distance();

  if (d1 < d2)
    return -1;

  if (d1 > d2)
    return 1;

  return 0;
}

static int
distance_matrix_leader (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vt:k:S:f:m:s:H:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('t'))
  {
    cerr << "Must specify a threshold via the -t option\n";
    usage(4);
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', global_threshold) || global_threshold <= static_cast<float>(0.0))
    {
      cerr << "The global distance threshold must be non-negative\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will cluster at " << global_threshold << endl;
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (s.starts_with("default="))
    {
      s.remove_leading_chars(8);
      default_smiles << smiles_tag << s << ">\n";
    }
    else if (! read_smiles(s, smiles))
    {
      cerr << "Cannot read smiles data '" << s << "'\n";
      return 3;
    }
    else if (verbose)
      cerr << "Read " << smiles.size() << " id->smiles relationships from '" << s << "'\n";
  }

  if (cl.option_present('f'))
  {
    if (! cl.value('f', fast_io) || fast_io < 2)
    {
      cerr << "The fast I/O option needs a precision of > 1 digit\n";
      usage(6);
    }

    if (verbose)
      cerr << "Will use fast I/O, " << fast_io << " decimal digits\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', global_max_cluster_size) || global_max_cluster_size < 1)
    {
      cerr << "The max cluster size (-m) option must be a whole positive number\n";
      usage(6);
    }

    if (verbose)
      cerr << "Max cluster size " << global_max_cluster_size << endl;
  }

  if (cl.option_present('H'))
  {
    const_IWSubstring h = cl.string_value('H');

    if (h.starts_with("col="))
    {
      h.remove_leading_chars(4);

      if (! h.numeric_value(threshold_column) || threshold_column < 1)
      {
        cerr << "Invalid threshold column '" << h << "'\n";
        usage(4);
      }

      if (verbose)
        cerr << "Thresholds in column " << threshold_column << endl;

      threshold_column--;
    }
    else
    {
      threshold_file_name = h;

      if (verbose)
        cerr << "Thresholds in '" << threshold_file_name << endl;
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(5);
  }

  const_IWSubstring d = cl[0];    // the name of the distance matrix

  const_IWSubstring k;
  if (cl.option_present('k'))
    k = cl.string_value('k');
  else
    k = "float";

  int rc;

  if ("float" == k)
  {
    ID_Distance_Comparitor<float> comparitor;
    IWDistanceMatrixBase<float> dm;
    rc = distance_matrix_leader(cl, dm, d, comparitor, cout);
  }
  else if ("double" == k)
  {
    ID_Distance_Comparitor<double> comparitor;
    IWDistanceMatrixBase<double> dm;
    rc = distance_matrix_leader(cl, dm, d, comparitor, cout);
  }
  else if ("int" == k)
  {
    IWDistanceMatrixBase<int> dm;
    ID_Distance_Comparitor<int> comparitor;
    rc = distance_matrix_leader(cl, dm, d, comparitor, cout);
  }
  else if ("uchar" == k)
  {
    IWDistanceMatrixBase<unsigned char> dm;
    ID_Distance_Comparitor<unsigned char> comparitor;
    rc = distance_matrix_leader(cl, dm, d, comparitor, cout);
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
    cerr << "Clustered " << number_molecules << " molecules into " << clusters_found << " clusters\n";
    int isum = 0;
    for (int i = 0; i < cluster_size.number_elements(); i++)
    {
      int j = cluster_size[i];
      if (0 == j)
        continue;

      cerr << j << " clusters were of size " << i << " members\n";

      isum += j * i;
    }

    cerr << "In clusters " << isum << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_leader(argc, argv);

  return rc;
}
