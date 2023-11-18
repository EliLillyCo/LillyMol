/*
  Takes a descriptor file and produces a gfp file with molecular
  properties This variant computes all pair-wise distances and scales
  things so the longest distance will be 1.0
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;


static int verbose = 0;

static IWString tag ("DSC<");
static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static IWString dummy_fingerprint ("zU..1;8;7;8;7;1");
static IWString dummy_fingerprint_tag ("FPJUNK<");

static IWString missing_value ('.');

static int produce_descriptor_file = 0;

static IW_STL_Hash_Map_String id_to_smiles;

/*
  We want to scan a record and look for instances of missing values.
*/

static IWString start_of_record_missing_value;
static IWString space_delimited_missing_value;
static IWString end_of_record_missing_value;

/*
  Sometimes we have numbers in the form 123E+4
*/

static int big_e = 0;

static int records_containing_missing_values = 0;

static float mahalanobis_scaling_factor = 0.0;

static IWHistogram histogram;

static Fraction_as_String fraction_as_string;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;

  cerr << "Takes a descriptor file and produces a gfp file with molecular properties\n";
  cerr << " -T <tag>    tag for molecular properties (default '" << tag << "')\n";
  cerr << " -M <string> ignore records with missing values - missing values are represented by <string>\n";
  cerr << " -E          input data contains numeric forms that include E\n";
  cerr << " -H <maxd>   initialise a histogram to report distances within <maxd>\n";
  cerr << " -d          produce a descriptor file\n";
  cerr << " -S <fname>  corresponding smiles file\n";
  cerr << " -v          verbose output\n";

  exit (rc);
}

class Descriptors
{
  private:
    IWString _id;

    float * _descriptor;

//  private functions

    int _build_from_descriptor_record (const const_IWSubstring &, int, int, int &);

  public:
    Descriptors ();
    ~Descriptors ();

    int build_from_descriptor_record (const const_IWSubstring &, int &);

    void scan  (Accumulator<float> * acc) const;
    void shift (float);

    float cartesian_distance (const Descriptors &) const;

    int write_gfp_record (int, float, IWString_and_File_Descriptor &) const;
};

Descriptors::Descriptors ()
{
  _descriptor = nullptr;

  return;
}

Descriptors::~Descriptors ()
{
  if (NULL != _descriptor)
    delete _descriptor;

  return;
}

static int
convert_from_text (const const_IWSubstring & token,
                   float & d,
                   int & fatal)
{
  if (token.numeric_value (d))
    return 1;

  if (missing_value == token)
  {
    fatal = 1;
    cerr << "Sorry, don't know how to handle missing values\n";
    return 0;
  }

  IWString copy_of_token (token);
  copy_of_token.gsub ('E', 'e');

  double tmp;

  if (! token.numeric_value (tmp))
  {
    cerr << "Invalid numeric '" << token << "'\n";
    fatal = 1;
    return 0;
  }

  d = static_cast<float> (tmp);

  return 1;
}

static int
contains_missing_values (const const_IWSubstring & buffer)
{
  if (0 == space_delimited_missing_value.length ())    // not checking
    return 0;

  if (buffer.contains (space_delimited_missing_value))
  {
    records_containing_missing_values++;
    return 1;
  }

  if (buffer.starts_with (start_of_record_missing_value) || buffer.ends_with (end_of_record_missing_value))
  {
    records_containing_missing_values++;
    return 1;
  }

  return 0;
}

static int number_descriptors = 0;

static int
contains_E_in_numeric_fields (const const_IWSubstring & buffer,
                              int istart)
{
  const char * s = buffer.data ();

  for (int i = istart; i < buffer.length (); i++)
  {
    if ('E' == s[i])
      return 1;
  }

  return 0;
}

int
Descriptors::build_from_descriptor_record (const const_IWSubstring & buffer,
                                           int & fatal)
{
  if (contains_missing_values (buffer))
  {
    fatal = 0;
    return 0;
  }

  if (0 == number_descriptors)
  {
    cerr << "Descriptors::build_from_descriptor_record: unknown number of descriptors\n";
    return 0;
  }

  assert (NULL == _descriptor);

  _descriptor = new float[number_descriptors];

  int i = 0;
  const_IWSubstring token;

  (void) buffer.nextword (_id, i);

  if (big_e && contains_E_in_numeric_fields (buffer, i))
  {
    IWString tmp (buffer);
    tmp.gsub ('E', 'e');
//  cerr << "Transformed to '" << tmp << "'\n";
    return _build_from_descriptor_record (tmp, i, number_descriptors, fatal);
  }

  return _build_from_descriptor_record (buffer, i, number_descriptors, fatal);
}

int
Descriptors::_build_from_descriptor_record (const const_IWSubstring & buffer,
                                           int i,
                                           int number_descriptors,
                                           int & fatal)
{
  int col = 0;
  const_IWSubstring token;

  while (buffer.nextword (token, i))
  {
    if (col >= number_descriptors)
    {
      cerr << "Yipes, too many columns " << col << endl;
      return 0;
    }

    int fatal;
    if (! convert_from_text (token, _descriptor[col], fatal))
    {
      cerr << "Descriptors::build_from_descriptor_record: cannot build\n";
      cerr << buffer << endl;
      return 0;
    }

    col++;
  }

  if (col != number_descriptors)
  {
    cerr << "Number descriptors mismatch " << col << " vs " << number_descriptors << endl;
    return 0;
  }

  return 1;
}

void
Descriptors::scan (Accumulator<float> * acc) const
{
  for (int i = 0; i < number_descriptors; i++)
  {
    acc[i].extra (_descriptor[i]);
  }

  return;
}

void
Descriptors::shift (float offset)
{
  for (int i = 0; i < number_descriptors; i++)
  {
    _descriptor[i] = _descriptor[i] - offset;
  }

  return;
}

float 
Descriptors::cartesian_distance (const Descriptors & rhs) const
{
  assert (number_descriptors > 0);

  float rc = static_cast<float> (0.0);

  for (int i = 0; i < number_descriptors; i++)
  {
    float diff = _descriptor[i] - rhs._descriptor[i];

    rc += (diff * diff);
  }

  return sqrt (rc);
}

/*
  We can save time with I/O if we pre-compute the string representations of the values found in the descriptors
*/

static int use_precomputed_descriptor_values_if_possible = 0;

static IWString * dist = nullptr;

static int
initialise_distances (float maxrange)
{
  int n = static_cast<int> (maxrange / 0.01) + 1;

  if (verbose)
    cerr << "Fast I/O requires " << n << " precomputed numbers\n";

  if (n > 5000)
    return 0;

  dist = new IWString[n + 1];

  dist[0] = '0';
  for (int i = 1; i <= n; i++)
  {
    float f = static_cast<float> (i) * static_cast<float> (0.01);

    if (f < 1.0)
      dist[i].append_number (f, 2);
    else if (f < 10.0)
      dist[i].append_number (f, 3);
    else
      dist[i].append_number (f, 4);

//  cerr << "i = " << i << " distance '" << dist[i] << "'\n";
  }

  return 1;
}

int
Descriptors::write_gfp_record (int seq,
                               float scale_factor,
                               IWString_and_File_Descriptor & output) const
{
  IWString molecular_properties;
  molecular_properties.resize (number_descriptors * 8);

  if (NULL != dist)
  {
    for (int i = 0; i < number_descriptors; i++)
    {
      int j = static_cast<int> (rint (_descriptor[i] * scale_factor * 100.0));

      molecular_properties.append_with_spacer (dist[j]);
    }
  }
  else
  {
    for (int i = 0; i < number_descriptors; i++)
    {
      if (i > 0)
        molecular_properties += ' ';
  
      float f = _descriptor[i] * scale_factor;
      molecular_properties << f;
    }
  }

  if (produce_descriptor_file)
  {
    output << _id << ' ' << molecular_properties << '\n';
    return output.good();
  }

  if (0 == id_to_smiles.size())
  {
    output << smiles_tag << "[" << seq << "C]>\n";
    output << identifier_tag << _id << ">\n";
    output << dummy_fingerprint_tag << dummy_fingerprint << ">\n";
    output << tag << molecular_properties << ">\n";
    output << "|\n";

    return 1;
  }

  if (! id_to_smiles.contains(_id))
  {
    cerr << "GACK, no smiles for '" << _id << "', cannot continue\n";
    return 0;
  }

  output << smiles_tag << id_to_smiles[_id] << ">\n";
  output << identifier_tag << _id << ">\n";
  output << tag << molecular_properties << ">\n";
  output << "|\n";

  return output.good ();
}

static Descriptors * descriptors = nullptr;
static int number_records = 0;

static int
write_gfp_file (IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < number_records; i++)
  {
    if (! descriptors[i].write_gfp_record (i, mahalanobis_scaling_factor, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good ();
}

/*
  within gfp_nearneighbours, we want all distances to be within the range [0,1]. Find the longest possible
  distance in the dataset.
*/

static int
determine_longest_distance()
{
  Accumulator<float> * acc = new Accumulator<float>[number_descriptors];
  std::unique_ptr<Accumulator<float>[]> free_acc(acc);

  for (int i = 0; i < number_records; i++)
  {
    descriptors[i].scan (acc);
  }

  float minmin = acc[0].minval ();
  float maxrange = acc[0].maxval () - acc[0].minval ();
  for (int i = 1; i < number_descriptors; i++)
  {
    float d = acc[i].minval ();
    if (d < minmin)
      minmin = d;

    d = acc[i].maxval () - d;
    if (d > maxrange)
      maxrange = d;
  }

  if (verbose)
    cerr << "Minimum coordinate " << minmin << endl;

  for (int i = 0; i < number_records; i++)
  {
    descriptors[i].shift (minmin);
  }

  float maxd = 0.0;

  for (int i = 0; i < number_records; i++)
  {
    const Descriptors & di = descriptors[i];

    if (verbose && i > 0 && 0 == i % 1000)
      cerr << "Computing closest distances " << i << endl;

    for (int j = i + 1; j < number_records; j++)
    {
      float d = di.cartesian_distance (descriptors[j]);

      if (d > maxd)
        maxd = d;

      if (histogram.active () && d <= histogram.maxval ())
        histogram.extra (d);
    }
  }

  if (verbose)
  {
    float max_range = acc[0].maxval () - acc[0].minval ();
    for (int i = 1; i < number_descriptors; i++)
    {
      if (0 == acc[i].n ())
        continue;

      float ri = acc[i].maxval () - acc[i].minval ();
      if (ri > max_range)
        max_range = ri;

      if (verbose > 1)
        cerr << "Descriptor " << i << " " << acc[i].n () << " values between " << acc[i].minval () << " and " << acc[i].maxval () << endl;
    }
    cerr << "Longest distance " << maxd << ", max range " << max_range << endl;

    histogram.debug_print (cerr);
  }

// We want to be able to multiply the factor rather than dividing it

  mahalanobis_scaling_factor = 1.0 / maxd;

  if (verbose > 1)
    cerr << "Multiplicative scaling factor " << mahalanobis_scaling_factor << endl;

  if (use_precomputed_descriptor_values_if_possible)
    initialise_distances (maxrange * mahalanobis_scaling_factor);

  return 1;
}

static int
build_pool (iwstring_data_source & input)
{
  int items_in_pool = 0;
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    int fatal;
    if (! descriptors[items_in_pool].build_from_descriptor_record (buffer, fatal))
    {
      if (! fatal)
        continue;

      cerr << "Cannot parse descriptor file record at line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }

    items_in_pool++;
  }

  number_records = items_in_pool;

  if (verbose)
    cerr << "Pool contains " << number_records << " items\n";

  return number_records;
}

static int
write_header (IWString & header,
              IWString_and_File_Descriptor & output)
{
  output << header << '\n';

  return output.good ();
}

static int
descriptor_file_to_molecular_properties (const Command_Line & cl,
                                         iwstring_data_source & input,
                                         IWString_and_File_Descriptor & output)
{
  IWString header;

  if (! input.next_record (header))
  {
    cerr << "Cannot fetch header record\n";
    return 0;
  }

  int columns_in_input = header.nwords ();

  if (verbose)
    cerr << "Input contains " << columns_in_input << " columns\n";

  number_descriptors = columns_in_input - 1;

  number_records = input.records_remaining ();
  if (0 == number_records)
  {
    cerr << "Yipes, only the header record in the file\n";
    return 0;
  }

  descriptors = new Descriptors[number_records];
  if (NULL == descriptors)
  {
    cerr << "Yipes, cannot allocate " << number_records << " objects\n";
    return 0;
  }

  std::unique_ptr<Descriptors[]> free_descriptors(descriptors);

  if (verbose)
    cerr << "There are " <<number_records << " records in the file\n";

  if (! build_pool (input))
  {
    cerr << "Could not read all the data\n";
    return 0;
  }

  (void) determine_longest_distance ();

  if (produce_descriptor_file)
    write_header (header, output);

  return write_gfp_file (output);
}

static int
descriptor_file_to_molecular_properties (const Command_Line & cl,
                                         IWString_and_File_Descriptor & output)
{
  const char * fname = cl[0];

  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptor_file_to_molecular_properties (cl, input, output);
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & id_to_smiles)
{
  const_IWSubstring buffer;

  IWString smiles, id;
  int i;

  while (input.next_record(buffer))
  {
    i = 0;
    buffer.nextword(smiles, i);
    buffer.nextword(id, i);

    if (id_to_smiles.contains(id))
      cerr << "Ignoring duplicate smiles identifier '" << id << "'\n";
    else
      id_to_smiles[id] = smiles;
  }

  return id_to_smiles.size();
}

static int
read_smiles (const char * fname,
             IW_STL_Hash_Map_String & id_to_smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open smiles file '" << fname << "'\n";
    return 0;
  }

  return read_smiles (input, id_to_smiles);
}

static int
descriptor_file_to_molecular_properties (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vT:M:EfH:dS:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('H'))
  {
    double h;
    if (! cl.value ('H', h) || h <= 0.0)
    {
      cerr << "INvalid histogram maximum distance (-H option)\n";
      usage (5);
    }

    histogram.initialise (0.0, h, h / 100.0);

    if (verbose)
      histogram.debug_print (cerr);
  }

  if (cl.option_present ('T'))
  {
    tag = cl.string_value ('T');

    if (verbose)
      cerr << "Molecular properties written as '" << tag << "' dataitem\n";

    if (! tag.ends_with ('<'))
      tag += '<';
  }

  if (cl.option_present ('E'))
  {
    big_e = 1;

    if (verbose)
      cerr << "Will translate E to e for numeric input\n";
  }

  if (cl.option_present ('f'))
  {
    use_precomputed_descriptor_values_if_possible = 1;

    if (verbose)
      cerr << "If possible, we will use precomputed fast I/O\n";
  }

  if (cl.option_present ('d'))
  {
    produce_descriptor_file = 1;

    if (verbose)
      cerr << "A descriptor file will be produced\n";
  }

  if (cl.option_present ('M'))
  {
    missing_value = cl.string_value ('M');

    if (verbose)
      cerr << "Missing value string '" << missing_value << "'\n";

    space_delimited_missing_value << ' ' << missing_value << ' ';
    start_of_record_missing_value << missing_value << ' ';
    end_of_record_missing_value << ' ' << missing_value;
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.number_elements () > 1)
  {
    cerr << "Sorry, only does one file at a time\n";
    return 7;
  }

  if (cl.option_present('S'))
  {
    const char * s = cl.option_value('S');

    if (! read_smiles (s, id_to_smiles))
    {
      cerr << "Cannot read smiles from '" << s << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Read " << id_to_smiles.size() << " smiles from '" << s << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! descriptor_file_to_molecular_properties (cl, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    if (records_containing_missing_values)
      cerr << records_containing_missing_values << " records contained missing values - discared\n";
  }

  return rc;
}


int
main (int argc, char ** argv)
{
  int rc = descriptor_file_to_molecular_properties (argc, argv);

  return rc;
}
