/*
  Takes distances from an ASCII file and produces a descriptor matrix
*/

#include <stdlib.h>
#include <math.h>

#include "cmdline.h"
#include "misc.h"
#include "iw_stl_hash_map.h"
#include "iw_stl_hash_set.h"
#include "IWDistanceMatrixBase.h"
#include "iwstring_data_source.h"
#include "iwmmap.h"
#include "report_progress.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString missing_value;

static int missing_values_encountered = 0;

static int square_matrix_contains_header = 0;

static int input_is_similarity = 0;

static int distance_matrix_size = 0;

static int updating_existing_distance_matrix = 0;

static int values_stored = 0;

static Report_Progress report_progress;

static int rescale_distances = 0;

static int three_column_data_has_header = 0;

static char input_separator = ' ';

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Takes a file of ASCII distances and produces a distance matrix file\n";
  cerr << " -S <fname>     distance matrix file to produce\n";
  cerr << " -t <form>      what kind of input\n";
  cerr << " -t 3           3 column format 'id1 id2 dist'\n";
  cerr << " -t 3h          3 column format including a header record (discarded)\n";
  cerr << " -t sq          square data - also works for lower triangular\n";
  cerr << " -t sqh         square data with a header record\n";
  cerr << " -t rh          square proximity data, identifiers in -I file\n";
  cerr << " -I <file>      file of identifiers, used with '-t rh'\n";
  cerr << " -u <float>     default value for distances - unset value\n";
  cerr << " -M <string>    missing value\n";
  cerr << " -m             input is a similarity measure, change to a distance before storing\n";
  cerr << " -s <size>      size of new distance matrix\n";
  cerr << " -a             automatically determine size of distance matrix (must read input file twice)\n";
  cerr << " -r <number>    report progress every <number> items processed\n";
  cerr << " -i             incremental update of existing database\n";
  cerr << " -y             rescale input distances to the range 0-1\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

template <typename T>
int
resize_distance_matrix (IWDistanceMatrixBase<T> & dm,
                        int n)
{
  if (n < 2)
  {
    cerr << "Not enough data in file, N = " << n << endl;
    return 0;
  }

  if (! dm.resize(n))
  {
    cerr << "Memory failure, N = " << n << endl;
    return 0;
  }

  if (verbose)
    cerr << "Problem sized for " << dm.number_molecules() << " molecules\n";

  distance_matrix_size = dm.number_molecules();

  return 1;
}

template <typename T>
int
read_header_record (iwstring_data_source & input,
                    IWDistanceMatrixBase<T> & dm)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  int columns_in_input = buffer.nwords() - 1;

  if (! resize_distance_matrix(dm, columns_in_input))
    return 0;

  return 1;
}

/*static int
read_row_of_three_column_data (const const_IWSubstring & buffer,
                               IW_STL_Hash_Map_int & id_to_index,
                               int row,
                               IWDistanceMatrixBase<float> & dm)
{

  int i = 0;
  IWString token;
  int col = 0;

  buffer.nextword(token, i);     // discard first token

  while (buffer.nextword(token, i))
  {
    if (id_to_index.contains(token))
    {
      cerr << "Duplicate identifier '" << token << "'\n";
      return 0;
    }

    int i = id_to_index.size();
    id_to_index[token] = i;
  }

  if (id_to_index.size() != n)
  {
    cerr << "Internal error, expected " << n << " identifiers, got " << id_to_index.size() << endl;
    return 0;
  }

  return 1;
}*/

template <typename T, typename S>
int
store_value (IWDistanceMatrixBase<T> & dm,
             int row,
             int col,
             S & token)
{
  if (token == missing_value)
  {
    missing_values_encountered++;
    return 1;
  }

  T d;
  if (! token.numeric_value(d))
  {
    cerr << "Invalid numeric value '" << token << "'\n";
    return 0;
  }

// These next lines are to avoid compiler warnings. If we do 'if (d < 0)' we get warnings when using unsigned variables

  if (d > static_cast<T>(0))
    ;
  else if (static_cast<T>(0) == d)
    ;
  else
  {
    cerr << "Negative values not allowed '" << token << "'\n";
    return 0;
  }

  if (input_is_similarity)
  {
    assert (d >= static_cast<float>(0.0) && d <= static_cast<float>(1.0));
    d = static_cast<T>(static_cast<float>(1.0) - static_cast<float>(d));
  }

  if (row >= dm.number_molecules() || col >= dm.number_molecules() || row == col)
  {
    cerr << "Invalid row/col combination " << row << '/' << col << " number molecules " << dm.number_molecules() << endl;
    return 0;
  }

//cerr << "Distance between " << row << " and " << col << " is '" << d << endl;

  values_stored++;

  return dm.set(row, col, d);
}

template <typename T>
int
read_square_data_record (const const_IWSubstring & buffer,
                         int row,
                         IWDistanceMatrixBase<T> & dm)
{
  if (0 == dm.number_molecules())
  {
    if (! resize_distance_matrix(dm, buffer.nwords() - 1))
      return 0;
  }

  assert (row < dm.number_molecules());

  int i = 0;
  const_IWSubstring token;

  (void) buffer.nextword(token, i);
  dm.set_id(row, token);

  int col = 0;

  while (buffer.nextword(token, i))
  {
    if (col < row)
    {
      if (! store_value(dm, row, col, token))
        return 0;
    }

    col++;

    if (col >= row)
      break;
  }

  if (col < row)
  {
    cerr << "Not enough columns, got " << col << " expected " << row << endl;
    return 0;
  }

  return 1;
}

template <typename T>
int
read_square_data (iwstring_data_source & input,
                  IWDistanceMatrixBase<T> & dm)
{
  if (square_matrix_contains_header)
  {
    if (! read_header_record(input, dm))
    {
      cerr << "Cannot read header record\n";
      return 0;
    }
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! read_square_data_record(buffer, input.lines_read() - 1 - square_matrix_contains_header, dm))
    {
      cerr << "Fatal error reading line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    if (report_progress())
      cerr << "Processed " << input.lines_read() << " lines, stored " << values_stored << " values\n";
  }

  if (input.lines_read() - square_matrix_contains_header != dm.number_molecules())
  {
    cerr << "Invalid data. Only read " << input.lines_read() << " records for " << dm.number_molecules() << " molecules\n"; 
    return 0;
  }

  return 1;
}

template <typename T>
int
read_square_data (const char * fname,
                  IWDistanceMatrixBase<T> & dm)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_square_data(input, dm);
}

template <typename T>
int
assign_index (IW_STL_Hash_Map_int & id_to_index,
              const IWString & id,
              IWDistanceMatrixBase<T> & dm)
{
  IW_STL_Hash_Map_int::const_iterator f = id_to_index.find(id);

  if (f != id_to_index.end())
    return (*f).second;

  int rc = id_to_index.size();

  if (rc >= distance_matrix_size)
  {
    cerr << "Gack, distance matrix sized to " << distance_matrix_size << ", too many unique identifiers\n";
    return -1;
  }

  id_to_index[id] = rc;

  dm.set_id(rc, id);

  return rc;
}

template <typename T>
int
read_proximity (const const_IWSubstring & buffer,
                int zrow,
                IWDistanceMatrixBase<T> & dm)
{
  const_IWSubstring token;

  int col = 0;
  for (int i = 0; buffer.nextword(token, i); col++)
  {
    T v;

    if (token.numeric_value(v))
      ;
    else if (missing_value == token)
    {
      missing_values_encountered++;
      continue;
    }

    if (col == zrow)
      continue;

    if (input_is_similarity)
      v = static_cast<T>(1.0) - v;

    dm.set(zrow, col, v);

    values_stored++;
  }

  if (col != dm.number_molecules())
  {
    cerr << "Incomplete data for row " << zrow << " expected " << dm.number_molecules() << " items\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
read_proximity (iwstring_data_source & input,
                IWDistanceMatrixBase<T> & dm)
{
  const_IWSubstring buffer;

  int zrow;
  for (zrow = 0; input.next_record(buffer); zrow++)
  {
    if (zrow >= dm.number_molecules())
    {
      cerr << "Too much data in input file, matrix sized for " << dm.number_molecules() << " molecules\n";
      return 0;
    }

    if (! read_proximity(buffer, zrow, dm))
    {
      cerr << "Fatal error processing line " << (zrow + 1) << endl;
      return 0;
    }

    if (report_progress())
      cerr << "Read " << input.lines_read() << " lines, stored " << values_stored << " values\n";
  }

  if (zrow != dm.number_molecules())
  {
    cerr << "Incomplete data, only read " << zrow << " rows of data, expected " << dm.number_molecules() << endl;
    return 0;
  }

  return 1;
}

template <typename T>
int
read_proximity (const char * fname,
                IWDistanceMatrixBase<T> & dm)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open proximity data file '" << fname << "'\n";
    return 0;
  }

  return read_proximity(input, dm);
}

template <typename T>
int
read_identifier_and_proximity (iwstring_data_source & identifier_input,
                               const char * fname,
                               IWDistanceMatrixBase<T> & dm)
{
  int n = identifier_input.records_remaining();

  if (n < 2)
  {
    cerr << "Identifier file contains " << n << " records, impossible\n";
    return 0;
  }

  if (! dm.resize(n))
  {
    cerr << "Cannot allocate " << n << " records in distance matrix\n";
    return 0;
  }

  if (verbose)
    cerr << "Problem sized for " << n << " records\n";

  const_IWSubstring buffer;

  for (int i = 0; identifier_input.next_record(buffer); i++)
  {
    dm.set_id(i, buffer);
  }

  return read_proximity(fname, dm);
}

template <typename T>
int
read_identifier_and_proximity (const IWString & identifier_fname,
                               const char * fname,
                               IWDistanceMatrixBase<T> & dm)
{
  iwstring_data_source input(identifier_fname);

  if (! input.good())
  {
    cerr << "Cannot open identifier file '" << identifier_fname << "'\n";
    return 0;
  }

  return read_identifier_and_proximity(input, fname, dm);
}

template <typename T>
int
read_three_column_data (const const_IWSubstring & buffer,
                        IW_STL_Hash_Map_int & id_to_index,
                        IWDistanceMatrixBase<T> & dm)
{
  int i = 0;
  IWString token;

  if (! buffer.nextword(token, i, input_separator))
  {
    cerr << "empty record\n";
    return 0;
  }

  int id1 = assign_index(id_to_index, token, dm);

  if (id1 < 0)
    return 0;

//cerr << "Name '" << token << "' assigned " << id1 << endl;

  if (! buffer.nextword(token, i, input_separator))
  {
    cerr << "Cannot extract ID2\n";
    return 0;
  }

  int id2 = assign_index(id_to_index, token, dm);

  if (id2 < 0)
    return 0;

//cerr << "Name '" << token << "' assigned " << id2 << endl;

  if (id1 == id2)   // don't care about diagonal items
    return 1;

  if (! buffer.nextword(token, i, input_separator))
  {
    cerr << "Cannot extract distance\n";
    return 0;
  }

  return store_value(dm, id1, id2, token);
}

/*
  When reading three column input, we have various possibilities for the number of
  records in the input
*/

static int
try_various_things (int nr)
{
  int n = static_cast<int>(1.0 + sqrt(static_cast<double>(1 + 8 * nr)) / 2.0 + 0.0001);  // NR = n * (n - 1) / 2

  if (n * (n - 1) / 2 == nr)
    return n;

  n = static_cast<int>(1.0 + sqrt(static_cast<double>(1 + 4 * nr)) / 2.0 + 0.0001);      // NR = n * (n - 1)

  if (n * (n - 1) == nr)
    return n;

  n = static_cast<int>(sqrt(static_cast<double>(nr)) + 0.001);                           // NR = n * n

  if (n * n == nr)
    return n;

  return 0;
}

/*
  We are reading three column data, how big is the corresponding distance matrix?
*/

template <typename T>
int
resize_based_on_records_in_input (IWDistanceMatrixBase<T> & dm,
                                  int nr)
{
  if (nr < 2)
  {
    cerr << "Three column input must contain at least two records\n";
    return 0;
  }

  if (verbose)
    cerr << "Input contains " << nr << " records\n";

  if (three_column_data_has_header)
    nr--;

/* nr is actually n* (n - 1) and so must be even

  if (nr != (nr / 2) * 2)
  {
    cerr << "Number of records in file must be an even number, " << nr << " is invalid\n";
    return 0;
  }*/

  int n = try_various_things(nr);

  if (0 == n)
  {
    cerr << "Cannot be " << nr << " records in a distance matrix file, n = " << n << "\n";
    return 0;
  }

  if (verbose)
    cerr << nr << " records in three column distance file corresponds to " << n << " molecules\n";

  return resize_distance_matrix(dm, n);
}

template <typename T>
int
read_identifiers_from_database (const IWDistanceMatrixBase<T> & dm,
                                IW_STL_Hash_Map_int & id_to_index)
{
  int n = dm.number_molecules();

  for (int i = 0; i < n; i++)
  {
    const IWString & id = dm.id(i);

    if (id.length() > 0)
      id_to_index[id] = i;
  }

  return 1;    // even if no valid identifiers found
}

template <typename I, typename T>
int
read_three_column_data (I & input,
                        IWDistanceMatrixBase<T> & dm)
{
//input.set_translate_tabs(1);

  if (0 == distance_matrix_size)
  {
    int nr = input.records_remaining();

    if (! resize_based_on_records_in_input(dm, nr))
      return 0;
  }

  IW_STL_Hash_Map_int id_to_index;

  if (updating_existing_distance_matrix)
  {
    if (! read_identifiers_from_database(dm, id_to_index))
      return 0;
  }

  const_IWSubstring buffer;

  if (three_column_data_has_header)
    input.next_record(buffer);

  while (input.next_record(buffer))
  {
    if (! read_three_column_data(buffer, id_to_index, dm))
    {
      cerr << "Invalid 3 column record, line " << input.lines_read() << " '" << buffer << "'\n";
      return 0;
    }

    if (report_progress())
      cerr << "Stored " << values_stored << " values\n";
  }

  if (static_cast<int>(id_to_index.size()) != dm.number_molecules())
    cerr << "Warning, distance matrix sized for " << dm.number_molecules() << " but only " << id_to_index.size() << " unique identifiers found\n";

  return 1;
}

template <typename T>
int
read_three_column_data (const char * fname,
                        IWDistanceMatrixBase<T> & dm)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_three_column_data(input, dm);
}

static int
gather_ids (IW_STL_Hash_Set & ids,
            const const_IWSubstring & buffer)
{
  int i = 0;
  IWString token;
  if (! buffer.nextword(token, i, input_separator))
    return 0;

  ids.insert(token);

  if (! buffer.nextword(token, i, input_separator))
    return 0;

  ids.insert(token);
  
  return 1;
}

template <typename T>
int
read_three_column_data_auto_size (IWString_Data_Source_MMAP & input,
                                  IWDistanceMatrixBase<T> & dm)
{
  IW_STL_Hash_Set ids;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! gather_ids(ids, buffer))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  if (verbose)
    cerr << "Found " << ids.size() << " identifiers in three column input\n";

  if (! resize_distance_matrix(dm, ids.size()))
  {
    cerr << "Cannot resize distance matrix " << ids.size() << endl;
    return 0;
  }

  if (! input.seekg(0))
  {
    cerr << "read_three_column_data_auto_size:cannot seek back to start of file\n";
    return 0;
  }

  return read_three_column_data(input, dm);
}

template <typename T>
int
read_three_column_data_auto_size (const char * fname,
                                  IWDistanceMatrixBase<T> & dm)
{
  IWString_Data_Source_MMAP input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_three_column_data_auto_size(input, dm);
}

template <typename T>
int
read_three_column_data (Command_Line & cl,
                        IWDistanceMatrixBase<T> & dm)
{
  if (cl.option_present('s'))
  {
    if (! cl.value('s', distance_matrix_size) || distance_matrix_size < 2)
    {
      cerr << "A distance matrix must have at least two items in it (-s option)\n";
      return 0;
    }

    if (! resize_distance_matrix(dm, distance_matrix_size)) return 0;
    return read_three_column_data(cl[0], dm);
  }
  else if (cl.option_present('a'))
  {
    return read_three_column_data_auto_size(cl[0], dm);
  }
  else
  {
	  //GH Added to fix the bug
    return read_three_column_data(cl[0], dm);
  }

  return 1;    // keep the compiler quiet
}

template <typename T>
class Scale_Distances
{
  private:
    const T _maxval;
  public:
    Scale_Distances (T m) : _maxval(m) {}

    void operator () (T &) const;
};

template <typename T>
void
Scale_Distances<T>::operator () (T & v) const
{
  if (static_cast<T>(0.0) == v)
    return;

  v = v / _maxval;

  return;
}

template <typename T>
int
do_rescale_distances (IWDistanceMatrixBase<T> & dm)
{
  T m = dm.maxval();

  Scale_Distances<T> sd(m);

  dm.each(sd);

  return 1;

  for (typename IWDistanceMatrixBase<T>::iterator i = dm.begin(); i != dm.end(); i++)
  {
    sd(*i);
  }

  return 1;
}

template <typename T>
int
distance_matrix_from_distances (Command_Line & cl,
                                IWDistanceMatrixBase<T> & dm)
{
  IWString s = cl.string_value('S');

  if (! updating_existing_distance_matrix)    // no need to check on existing file
    ;
  else if (! dash_s(s.null_terminated_chars()))
  {
    cerr << "Existing distance matrix '" << s << "' not found\n";
    return 0;
  }
  else
  {
    if (! dm.do_read(s))
    {
      cerr << "Cannot read existing distance matrix '" << s << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Read existing database with " << dm.number_molecules() << " molecules\n";

    distance_matrix_size = dm.number_molecules();
  }

  if (! updating_existing_distance_matrix)
  {
    T unset_value = static_cast<T>(1);

    if (cl.option_present('u'))
    {
      IWString tmp = cl.string_value('u');
      if (! tmp.numeric_value(unset_value))
      {
        cerr << "The unset value must be a valid number\n";
        usage(5);
      }
  
      if (verbose)
        cerr << "Matrix initialised with " << unset_value << endl;
    }

    dm.set_initialiser(unset_value);
  }

  const_IWSubstring t = cl.string_value('t');
  
  int rc;

  if ("3" == t)
    rc = read_three_column_data(cl, dm);
  else if ("3h" == t)
  {
    three_column_data_has_header = 1;
    rc = read_three_column_data(cl, dm);
  }
  else if ("sq" == t)
    rc = read_square_data(cl[0], dm);
  else if ("sqh" == t)
  {
    square_matrix_contains_header = 1;

    rc = read_square_data(cl[0], dm);
  }
  else if ("rh" == t)
  {
    if (! cl.option_present('I'))
    {
      cerr << "Must also specify the file of identifiers via the -I option\n";
      usage(5);
    }

    IWString identifier_fname = cl.string_value('I');

    rc = read_identifier_and_proximity(identifier_fname, cl[0], dm);
  }
  else
  {
    cerr << "Unrecognised -t qualifier '" << t << "'\n";
    return 0;
  }

  if (0 == rc)
    return 0;

  if (verbose)
    cerr << "Stored " << values_stored << " items, distance matrix values between " << dm.minval() << " and " << dm.maxval() << endl;

  if (rescale_distances)
    do_rescale_distances(dm);

  return dm.do_write(s);
}

static int
distance_matrix_from_distances (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vS:t:u:M:ms:iI:r:yad:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('S'))
  {
    cerr << "Must specify name of distance matrix file via the -S option\n";
    usage(1);
  }

  if (cl.option_present('s') && cl.option_present('i'))
  {
    cerr << "The -i and -s options are mutually incompatible\n";
    usage(5);
  }

  if (cl.option_present('d'))
  {
    const_IWSubstring d = cl.string_value('d');
    if ("tab" == d)
      input_separator = '\t';
    else if ("comma" == d)
      input_separator = ',';
    else if ("space" == d)
      input_separator = ' ';
    else if (1 == d.length())
      input_separator = d[0];
    else
    {
      cerr << "Unrecognised input separator (-d) '" << d << "'\n";
      return 0;
    }
  }

  if (cl.option_present('i'))
  {
    updating_existing_distance_matrix = 1;

    if (verbose)
      cerr << "Update of existing distance matrix\n";
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 1;
    }
  }

  if (cl.option_present('m'))
  {
    input_is_similarity = 1;

    if (verbose)
      cerr << "Input is similarity, will convert to distances\n";
  }

  if (cl.option_present('y'))
  {
    rescale_distances = 1;

    if (verbose)
      cerr << "Will rescale input distances\n";
  }

  if (cl.option_present('M'))
  {
    missing_value = cl.string_value('M');

    if (verbose)
      cerr << "Missing value string '" << missing_value << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (! cl.option_present('t'))
  {
    cerr << "Must specify type input via the -t option\n";
    usage(1);
  }

  const_IWSubstring k;

  if (cl.option_present('k'))
    k = cl.string_value('k');
  else
    k = "float";

  int rc = 0;

  if ("float" == k)
  {
    IWDistanceMatrixBase<float> dm;
    rc = distance_matrix_from_distances(cl, dm);
  }
  else if ("double" == k)
  {
    IWDistanceMatrixBase<double> dm;
    rc = distance_matrix_from_distances(cl, dm);
  }
  else if ("int" == k)
  {
    IWDistanceMatrixBase<int> dm;
    rc = distance_matrix_from_distances(cl, dm);
  }
  else if ("byte" == k)
  {
    IWDistanceMatrixBase<unsigned char> dm;
    rc = distance_matrix_from_distances(cl, dm);
  }
  else
  {
    cerr << "Not sure what to do with distance matrix type '" << k << "'\n";
    return 6;
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

  int rc = distance_matrix_from_distances(argc, argv);

  return rc;
}
