/*
  Takes a descriptor file and produces a gfp file with molecular properties
*/

#include <stdlib.h>
#include <math.h>

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

static int verbose = 0;

static IWString tag ("DSC<");

static IWString dummy_fingerprint ("zU..1;8;7;8;7;1");
static IWString dummy_fingerprint_tag ("FPJUNK<");

static IWString missing_value ('.');

/*
  We want to scan a record and look for instances of missing values.
*/

static IWString space_delimited_missing_value;
static IWString start_of_record_missing_value;
static IWString end_of_record_missing_value;

static int records_containing_missing_values = 0;

static double mahalanobis_scaling_factor = 0.0;

static int records_read = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;

  cerr << "Takes a descriptor file and produces a gfp file with molecular properties\n";
  cerr << " -T <tag>    tag for molecular properties (default '" << tag << "')\n";
  cerr << " -M <string> ignore records with missing values - missing values are represented by <string>\n";
  cerr << " -m          Mahalanobis data - normalise the results\n";
  cerr << " -v          verbose output\n";

  exit (rc);
}

static int
report_extremeties (int number_descriptors,
                    const Accumulator<float> * acc,
                    ostream & output)
{
  int columns_needing_shifting = 0;

  for (int i = 0; i < number_descriptors; i++)
  {
    output << " descriptor " << i << " between " << acc[i].minval () << " and " << acc[i].maxval () << '\n';
    if (acc[i].minval () < 0.0)
      columns_needing_shifting++;
  }

  output << columns_needing_shifting << " of " << number_descriptors << " descriptors need shifting\n";

  return output.good ();
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

static int
descriptor_file_to_molecular_properties (const const_IWSubstring & buffer,
                                         int number_descriptors,
                                         const Accumulator<float> * acc,
                                         IWString_and_File_Descriptor & output)
{
  if (contains_missing_values (buffer))
    return 1;

  int i = 0;
  const_IWSubstring id;
  buffer.nextword (id, i);

  IWString molecular_properties;
  molecular_properties.resize (number_descriptors * 7);

  int col = 0;
  const_IWSubstring token;

  while (buffer.nextword (token, i))
  {
    float d;
    
    token.numeric_value (d);

    if (0.0 != mahalanobis_scaling_factor)
      d = (d - acc[col].minval ()) * mahalanobis_scaling_factor;
    else if (acc[col].minval () < 0.0)
      d = d - acc[col].minval () + static_cast<float> (1.0);

    if (d < static_cast<float> (0.0))
      cerr << "Yipes, negative descriptor value " << d << " from '" << token << "'\n";

    if (col > 0)
      molecular_properties.add (' ');

    molecular_properties.append_number (d, 3);
    col++;
  }

  output << "$SMI<[" << records_read << "C]>\n";
  output << "PCN<" << id << ">\n";
  output << dummy_fingerprint_tag << dummy_fingerprint << ">\n";
  output << tag << molecular_properties << ">\n";
  output << "|\n";

  return output.good ();
}

static int
descriptor_file_to_molecular_properties (iwstring_data_source & input,
                                         int number_descriptors,
                                         const Accumulator<float> * acc,
                                         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  input.next_record (buffer);
   
  records_read = 0;

  while (input.next_record (buffer))
  {
    records_read++;

    if (! descriptor_file_to_molecular_properties (buffer, number_descriptors, acc, output))
    {
      cerr << "Fatal error writing " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return output.good ();
}

/*
  within gfp_nearneighbours, we want all distances to be within the range [0,1]. Find the longest possible
  distance in the dataset.
*/

static int
determine_mahalanobis_scaling_factor (const Accumulator<float> * acc,
                                      int number_descriptors)
{
  double d = 0.0;

  for (int i = 0; i < number_descriptors; i++)
  {
    float diff = acc[i].maxval () - acc[i].minval ();

    d += (diff * diff);
  }

  mahalanobis_scaling_factor = sqrt (d);

  if (verbose)
    cerr << "Mahalanobis scaling factor " << mahalanobis_scaling_factor << endl;

// We want to be able to multiply the factor rather than dividing it

  mahalanobis_scaling_factor = 1.0 / mahalanobis_scaling_factor;

  return 1;
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
scan_for_extremeties (const const_IWSubstring & buffer,
                      int number_descriptors,
                      Accumulator<float> * acc)
{
  if (contains_missing_values (buffer))
    return 1;

  int i = 0;
  const_IWSubstring token;

  (void) buffer.nextword (token, i);    // identifier in column 1

  int col = 0;
  while (buffer.nextword (token, i))
  {
    if (col >= number_descriptors)
    {
      cerr << "Yipes, too many columns " << col << endl;
      return 0;
    }

    float d;
    int fatal;
    if (! convert_from_text (token, d, fatal))
    {
      if (fatal)
        return 0;
    }
    else
      acc[col].extra (d);

    col++;
  }

  if (col != number_descriptors)
  {
    cerr << "Number descriptors mismatch " << col << " vs " << number_descriptors << endl;
    return 0;
  }

  return 1;
}

static int
scan_for_extremeties (iwstring_data_source & input,
                      int number_descriptors,
                      Accumulator<float> * acc)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! scan_for_extremeties (buffer, number_descriptors, acc))
    {
      cerr << "Fatal error on line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
descriptor_file_to_molecular_properties (const Command_Line & cl,
                                         iwstring_data_source & input,
                                         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (! input.next_record (buffer))
  {
    cerr << "Cannot fetch header record\n";
    return 0;
  }

  int columns_in_input = buffer.nwords ();

  if (verbose)
    cerr << "Input contains " << columns_in_input << " columns\n";

  int number_descriptors = columns_in_input - 1;

  Accumulator<float> * acc = new Accumulator<float>[number_descriptors];

  if (! scan_for_extremeties (input, number_descriptors, acc))
  {
    cerr << "Cannot determine extremeties in descriptors\n";
    return 0;
  }

  if (verbose > 1)
    report_extremeties (number_descriptors, acc, cerr);

  if (! input.seekg (0))
  {
    cerr << "Cannot seek back to beginning of file\n";
    return 0;
  }

  if (cl.option_present ('m'))
  {
    determine_mahalanobis_scaling_factor (acc, number_descriptors);
  }

  for (int i = 0; i < number_descriptors; i++)
  {
    if (0 == acc[i].n ())
    {
      cerr << "Sorry no values for column " << i << endl;
      return 3;
    }
  }

  int rc = descriptor_file_to_molecular_properties (input, number_descriptors, acc, output);

  delete [] acc;

  return rc;
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
descriptor_file_to_molecular_properties (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vT:M:m");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('T'))
  {
    tag = cl.string_value ('T');

    if (verbose)
      cerr << "Molecular properties written as '" << tag << "' dataitem\n";

    if (! tag.ends_with ('<'))
      tag += '<';
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
    cerr << "Read " << records_read << " records\n";
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
