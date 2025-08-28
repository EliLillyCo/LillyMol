/*
  Spread implementation with distances from a Distance Matrix file
*/

#include <stdlib.h>
#include <memory>
#include <random>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "IWDistanceMatrixBase.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int items_to_select = 0;

static int items_selected = 0;

static int behave_stochastically = 0;

static float distance_fuzz = static_cast<float>(0.0);

static IW_STL_Hash_Map_String id_to_smiles;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static IWString default_smiles;

static int distance_scale_specified = 0;

static float * distance_scale = NULL;

static int every_object_must_have_a_scale_factor = 1;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Spread implementation, distances from a distance matrix file\n";
  cerr << " -n <number>    number of items to select\n";
  cerr << " -d             choose between identical distances randomly\n";
  cerr << " -s <id>        first select <id>\n";
  cerr << " -r             select first item randomly\n";
  cerr << " -S <fname>     file of smiles\n";
  cerr << " -S default=XX  no smiles data available, use XX as smiles for all id's\n";
  cerr << " -A <fname>     file name with ID's already selected\n";
  cerr << " -p FILE=<fname> file of distance scaling factors\n";
  cerr << " -P COL=<col>    scale factors are in column <col>\n";
  cerr << " -k <type>      type of database 'byte', 'double', 'int', 'float', 'fmqb' (float masquerading byte\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

template <typename T>
int
process_distance_scaling_record (const IWDistanceMatrixBase<T> & dm,
                                 const const_IWSubstring & buffer,
                                 float * distance_scale)
{
  const_IWSubstring id, s;

  if (! buffer.split(id, ' ', s))
    return 0;

  int i = dm.which_item_has_id(id);

  if (i < 0)
  {
    cerr << "No distance matrix data for '" << id << "'\n";
    return 0;
  }

  float v;

  if (! s.numeric_value(v) || v <= 0.0F)
  {
    cerr << "Invalid scaling value '" << s << "'\n";
    return 0;
  }

  distance_scale[i] = v;

//cerr << "Scaling factor for '" << id << "' (" << i << ") set to " << v << endl;

  return 1;
}

template <typename D, typename T>
int
assign_scaling_from_name_token (const D & dm,
                                const T * distance_to_previously_selected,
                                int col)
{
  int n = dm.number_molecules();

  int items_with_no_scale_factor = 0;

  for (int i = 0; i < n; i++)
  {
    if (static_cast<T>(0) == distance_to_previously_selected[i])    // already selected, no need for a scaling factor
      continue;

    const IWString & id = dm.id(i);
    const_IWSubstring token;
    if (! id.word(col, token))
    {
      if ( ! every_object_must_have_a_scale_factor)
      {
        items_with_no_scale_factor++;
        continue;
      }

      if (! token.numeric_value(distance_scale[i]) || distance_scale[i] < 0.0f)
      {
        cerr << "Invalid column " << (col+1) << " for distance scale in '" << id << "'\n";
        return 0;
      }
    }
  }

  if (items_with_no_scale_factor)
    cerr << items_with_no_scale_factor << " items had no scale factor\n";

  return 1;
}

template <typename D, typename T>
int
read_distance_scaling (const D & dm,
                       const T * distance_to_previously_selected,
                       iwstring_data_source & input)
{
  const_IWSubstring buffer;

  int n = dm.number_molecules();

  input.next_record(buffer);    // discard header

  while (input.next_record(buffer))
  {
    if (! process_distance_scaling_record(dm, buffer, distance_scale))
    {
      cerr << "Fatal error processing distance scaling data '" << buffer << "'\n";
      return 0;
    }
  }

  int items_with_no_scaling_value = 0;

  for (int i = 0; i < n; i++)
  {
    if (distance_scale[i] > 0.0f)
      continue;

    if (static_cast<T>(0) == distance_to_previously_selected[i])    // already selected, no need for a scaling factor
      continue;

    if (every_object_must_have_a_scale_factor)
    {
      items_with_no_scaling_value++;
      cerr << "No scaling data for '" << dm.id(i) << "'\n";
    }
    else
      distance_scale[i] = 1.0f;
  }

  if (items_with_no_scaling_value && every_object_must_have_a_scale_factor)
    return 0;

  return 1;
}

template <typename D, typename T>
int
read_distance_scaling (const D & dm,
                       const T * distance_to_previously_selected,
                       const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open distance scaling file '" << fname << "'\n";
    return 0;
  }

  return read_distance_scaling(dm, distance_to_previously_selected, input);
}

static int
process_smiles (const const_IWSubstring & smiles,
                const_IWSubstring & id,
                IW_STL_Hash_Map_String & id_to_smiles)
{
  id.truncate_at_first(' ');

  if (id_to_smiles.contains(id))
  {
    cerr << "Duplicate identifier '" << id << "'\n";
    return 0;
  }

  id_to_smiles[id] = smiles;

  return 1;
}

static int
process_smiles (const const_IWSubstring & buffer,
                IW_STL_Hash_Map_String & id_to_smiles)
{
  const_IWSubstring smiles, id;

  if (! buffer.split(smiles, ' ', id))
  {
    cerr << "Cannot tokenise buffer into smiles and id\n";
    return 0;
  }

  return process_smiles(smiles, id, id_to_smiles);
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & id_to_smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! process_smiles(buffer, id_to_smiles))
    {
      cerr << "Fatal error processing smiles info '" << buffer << "'\n";
      return 0;
    }
  }

  return id_to_smiles.size();
}

static int
read_smiles (const const_IWSubstring & fname,
             IW_STL_Hash_Map_String & id_to_smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open smiles file '" << fname << "'\n";
    return 0;
  }

  return read_smiles(input, id_to_smiles);
}

template <typename T>
int
special_processing_for_first_item_selected (const IWDistanceMatrixBase<T> & dm,
                                            int zitem,
                                            T * distance_to_previously_selected,
                                            int * nearest_previously_selected_neighbour)
{
  int n = dm.number_molecules();

  for (int i = 0; i < n; i++)
  {
    if (i == zitem)
    {
      distance_to_previously_selected[i] = static_cast<T>(0);
      continue;
    }

    if (static_cast<T>(0.0) == distance_to_previously_selected[i])
      continue;

    T d = dm.zvalue(zitem, i);

    distance_to_previously_selected[i] = d * distance_scale[i];

    nearest_previously_selected_neighbour[i] = zitem;
  }

  return 1;
}

static int
write_smiles_and_id (const IWString & id,
                     IWString_and_File_Descriptor & output)
                     
{
  IW_STL_Hash_Map_String::const_iterator f = id_to_smiles.find(id);

  if (f != id_to_smiles.end())
    output << smiles_tag << (*f).second << ">\n";
  else if (default_smiles.length() > 0)
    output << default_smiles;

  output << identifier_tag << id << ">\n";

  return output.good();
}

template <typename D, typename T>
int
item_is_selected (const D & dm,
                  int zitem,
                  T * distance_to_previously_selected,
                  int * nearest_previously_selected_neighbour,
                  IWString_and_File_Descriptor & output,
                  bool do_output = true)
{
  if (verbose > 1)
  {
    cerr << "Selected item " << zitem << " '" << dm.id(zitem) << "'";
    if (items_selected)
      cerr << " dist " << distance_to_previously_selected[zitem];
    cerr << endl;
  }

  if (do_output)
  {
    write_smiles_and_id(dm.id(zitem), output);

    int nbr = nearest_previously_selected_neighbour[zitem];
    if (nbr >= 0)
    {
      write_smiles_and_id(dm.id(nbr), output);
      if (distance_scale_specified)
        output << "SCALE<" << distance_scale[zitem] << ">\n";
      output << distance_tag << distance_to_previously_selected[zitem] << ">\n";
    }
    output << "|\n";

    output.write_if_buffer_holds_more_than(4096);
  }

//distance_to_previously_selected[zitem] = static_cast<T>(-1.0f);

  items_selected++;

  nearest_previously_selected_neighbour[zitem] = -1;   // very important, mark as selected

//if (1 == items_selected)
//  return special_processing_for_first_item_selected(dm, zitem, distance_to_previously_selected, nearest_previously_selected_neighbour);

  int n = dm.number_molecules();

  static int first_call = 1;

  if (first_call)
  {
    for (int i = 0; i < n; i++)
    {
      if (i == zitem)
        continue;

      T d = dm.zvalue(zitem, i) * distance_scale[i];

      distance_to_previously_selected[i] = d;
      nearest_previously_selected_neighbour[i] = zitem;
//    cerr << "From '" << dm.id(zitem) << "' update " << dm.id(i) << " to " << d << endl;
    }

    first_call = 0;
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      if (nearest_previously_selected_neighbour[i] < 0)   // already selected
        continue;

      T d = dm.zvalue(zitem, i) * distance_scale[i];

      if (d > distance_to_previously_selected[i])
        continue;

      distance_to_previously_selected[i] = d;
      nearest_previously_selected_neighbour[i] = zitem;
//    cerr << "From '" << dm.id(zitem) << "' update " << dm.id(i) << " to " << d << endl;
    }
  }

  return output.good();
}

template <typename T>
int
choose_next_item (const int * nearest_previously_selected_neighbour,
                  const T * distance_to_previously_selected,
                  int n)
{
  int rc = -1;

  T largest_distance = static_cast<T>(-1.0);

  std::random_device rd;
  std::uniform_int_distribution<int> u01(0, 1);

  for (int i = 0; i < n; i++)
  {
    if (nearest_previously_selected_neighbour[i] < 0)   // already selected
      continue;

    if (rc < 0 || distance_to_previously_selected[i] > largest_distance)
    {
      largest_distance = distance_to_previously_selected[i];
      rc = i;
      continue;
    }

    if (! behave_stochastically)
      continue;

    if (largest_distance == distance_to_previously_selected[i])
    {
      if (1 == u01(rd))
        continue;

      largest_distance = distance_to_previously_selected[i];
      rc = i;
      continue;
    }

    if (largest_distance - distance_to_previously_selected[i] < distance_fuzz)
    {
      if (1 == u01(rd))
        continue;

      largest_distance = distance_to_previously_selected[i];
      rc = i;
      continue;
    }
  }

  assert (rc >= 0);

  return rc;
}

template <typename D, typename T>
int
process_previously_selected (D & dm,
                             int * nearest_previously_selected_neighbour,
                             T * distance_to_previously_selected,
                             iwstring_data_source & input,
                             int column,
                             IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (column > 0)
      buffer.remove_leading_words(1);

    buffer.truncate_at_first(' ');

    int s = dm.which_item_has_id(buffer);

    if (s < 0)
    {
      cerr << "process_previously_selected:unrecognised identifier '" << buffer << "'\n";
      return 0;
    }

    item_is_selected(dm, s, distance_to_previously_selected, nearest_previously_selected_neighbour, output, false);   // no output for previously selected
  }

  return 1;
}

template <typename D, typename T>
int
process_previously_selected (D & dm,
                             int * nearest_previously_selected_neighbour,
                             T * distance_to_previously_selected,
                             const char * fname,
                             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring tmp(fname);

  int column = 0;
  if (tmp.ends_with(".smi"))
    column = 1;

  return process_previously_selected(dm, nearest_previously_selected_neighbour, distance_to_previously_selected, input, column, output);
}

template <typename D, typename T>
int
distance_matrix_spread (D & dm,
                        int * nearest_previously_selected_neighbour,
                        T * distance_to_previously_selected,
                        IWString_and_File_Descriptor & output)
{
  int n = dm.number_molecules();

  while (items_selected < items_to_select)
  {
//  cerr << items_selected << " items selected\n";
    int i = choose_next_item(nearest_previously_selected_neighbour, distance_to_previously_selected, n);

#ifdef DEBUG_DM_SPREAD
    for (int i = 0; i < dm.number_molecules(); i++)
    {
      if (distance_to_previously_selected[i] > 0.0f)
        cerr << items_selected << " distance " << i << " '" << dm.id(i) << "' is " << distance_to_previously_selected[i] << endl;
    }
    cerr << "next " << i << " '" << dm.id(i) << "'\n";
#endif

    item_is_selected(dm, i, distance_to_previously_selected, nearest_previously_selected_neighbour, output);
  }

  return output.good();
}

template <typename T>
void
initialise_arrays (T * distance_to_previously_selected,
                   int * nearest_previously_selected_neighbour, 
                   int n)
{
  for (int i = 0; i < n; i++)
  {
    nearest_previously_selected_neighbour[i] = -1;
    distance_to_previously_selected[i] = static_cast<T>(-1.0);
  }

  return;
}

template <typename D, typename T>
int
distance_matrix_spread (Command_Line & cl,
                        D & dm,
                        T notused,
                        IWString_and_File_Descriptor & output)
{
  if (! dm.do_read(cl[0]))
  {
    cerr << "Cannot initialise distance matrix from '" << cl[0] << "'\n";
    return 6;
  }

  int n = dm.number_molecules();

  if (verbose)
    cerr << "Distance matrix built on " << n << " molecules\n";

  if (cl.option_present('S'))
  {
    const_IWSubstring s;
    for (int i = 0; cl.value('S', s, i); i++)
    {
      if (s.starts_with("default="))
      {
        s.remove_leading_chars(8);
        default_smiles << smiles_tag << s << ">\n";
      }
      else if (! read_smiles(s, id_to_smiles))
      {
        cerr << "Cannot read smiles data '" << s << "'\n";
        return 3;
      }
      else if (verbose)
        cerr << "Read " << id_to_smiles.size() << " id->smiles relationships from '" << s << "'\n";
    }
  }

  if (items_to_select > n)
  {
    cerr << "Requested " << items_to_select << " items, but distance matrix contains only " << n << endl;
    items_to_select = n;
  }

  if (0 == items_to_select)
    items_to_select = n;

  T * tmp = new T[n]; std::unique_ptr<T[]> free_tmp(tmp);
  int * itmp = new int[n]; std::unique_ptr<int[]> free_itmp(itmp);

  set_vector(tmp, n, static_cast<T>(1.0));

  distance_scale = new float[n];

  set_vector(distance_scale, n, -1.0f);

  if (cl.option_present('p'))
  {
    const_IWSubstring p;
    for (int i = 0; cl.value('p', p, i); i++)
    {
      if ("oknoscale" == p)
      {
        every_object_must_have_a_scale_factor = 0;

        if (verbose)
          cerr << "Not all objects need to have a scale factor\n";
      }
      else if (p.starts_with("FILE="))
      {
        p.remove_leading_chars(5);

        IWString fname(p);

        if (! read_distance_scaling(dm, tmp, fname.null_terminated_chars()))
        {
          cerr << "Cannot read distance scaling from '" << p << "'\n";
          return 0;
        }
      }
      else if (p.starts_with("COL="))
      {
        p.remove_leading_chars(4);
        int c;
        if (! p.numeric_value(c) || c < 2)
        {
          cerr << "The scaling column must be a +ve whole number > 1, 'COL=" << p << "' invalid\n";
          return 2;
        }

        if (verbose)
          cerr << "Distance scaling value in column " << c << endl;

        if (! assign_scaling_from_name_token(dm, tmp, c--))
        {
          cerr << "Could not assign scaling factors from name token\n";
          return 2;
        }
      }
      else
      {
        cerr << "Unrecognised scale specification '" << p << "'\n";
        usage(23);
      }
    }

    int n = dm.number_molecules();

    for (int i = 0; i < n; i++)
    {
      if (tmp[i] > static_cast<T>(0.0f))
        tmp[i] *= distance_scale[i];
    }

    distance_scale_specified = 1;
  }
  else
    set_vector(distance_scale, dm.number_molecules(), 1.0f);

  initialise_arrays(tmp, itmp, n);

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring s;

    while (cl.value('s', s, i++))
    {
      int j = dm.which_item_has_id(s);
      if (j < 0)
      {
        cerr << "Sorry, no item '" << s << "' in the set\n";
        return 5;
      }

      item_is_selected(dm, j, tmp, itmp, output);
    }
  }
  else if (cl.option_present('r'))
  {
    std::random_device rd;
    std::uniform_int_distribution<int> u(0, n - 1);
    int i = u(rd);
    item_is_selected(dm, i, tmp, itmp, output);
  }
  else if (cl.option_present('A'))
  {
    const char * a = cl.option_value('A');

    if (! process_previously_selected(dm, itmp, tmp, a, output))
    {
      cerr << "Cannot process previously selected file '" << a << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Processed " << items_selected << " previously selected items from '" << a << "'\n";
  }
  else
    item_is_selected(dm, 0, tmp, itmp, output);

  return distance_matrix_spread(dm, itmp, tmp, output);
}

static int
distance_matrix_spread(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:n:rS:f:k:dA:p:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (1 != cl.number_elements())
  {
    cerr << "Can process one, and only one, distance matrix file\n";
    usage(4);
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', items_to_select) || items_to_select < 1)
    {
      cerr << "The number of items to select (-n) option must be followed by a whole positive number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will select " << items_to_select << " items\n";
  }

  if (cl.option_present('d'))
  {
    behave_stochastically = 1;

    if (verbose)
      cerr << "Will randomly choose among equal distances\n";
  }

  if (cl.option_present('z'))
  {
    if (! cl.value('z', distance_fuzz) || distance_fuzz < static_cast<float>(0.0))
    {
      cerr << "The distance 'fuzz' value must be a valid non-negative floating point number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Distances within " << distance_fuzz << " will be decided randomly\n";
  }

  const_IWSubstring k;

  if (cl.option_present('k'))
    k = cl.string_value('k');
  else
    k = "float";

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if ("float" == k)
  {
    IWDistanceMatrixBase<float> dm;
    rc = distance_matrix_spread(cl, dm, 1.0f, output);
  }
  else if ("double" == k)
  {
    IWDistanceMatrixBase<double> dm;
    rc = distance_matrix_spread(cl, dm, 1.0, output);
  }
  else if ("int" == k)
  {
    IWDistanceMatrixBase<int> dm;
    rc = distance_matrix_spread(cl, dm, 1, output);
  }
  else if ("byte" == k)
  {
    IWDistanceMatrixBase<unsigned char> dm;
    rc = distance_matrix_spread(cl, dm, static_cast<unsigned char>(0), output);
  }
  else if ("fmqb" == k)
  {
    IWDistanceMatrixMasquerading_as_Byte<float> dm;
    rc = distance_matrix_spread(cl, dm, 1.0f, output);
  }
  else
  {
    cerr << "What am I supposed to do with type '" << k << "'\n";
    return 5;
  }

  output.flush();

  if (0 == rc)
    return 4;

  if (verbose)
  {
    cerr << "Selected " << items_selected << " molecules\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_spread(argc, argv);

  return rc;
}
