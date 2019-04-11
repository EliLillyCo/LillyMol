/*
  Scans a descriptor file for similarity to a given vector
*/

#include <stdlib.h>
#include <memory>
#include <limits>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_stl_hash_map.h"
#include "iwbits.h"
#include "misc.h"
#include "sparse_fp_creator.h"
using std::numeric_limits;

const char * prog_name = NULL;

static int verbose = 0;

static int records_read = 0;

static IWString tag("FPDSC<");

static IW_STL_Hash_Map_String id_to_smiles;

typedef IW_STL_Hash_Map<IWString, int *> ID_and_Vector;

static int columns_in_input = 0;

static int nbits = 0;

static int * numeric_value = NULL;

static IWString missing_value('.');

static int truncate_large_values = numeric_limits<int>::max();

static int strip_leading_zeros_from_identifiers = 0;

static int input_file_has_a_header = 1;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static IWString insert_dummy_smiles;

static int work_as_tdt_filter = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Converts a integer descriptor file to fingerprints, either fixed 0/1 or non colliding counted\n";
  cerr << " -F <tag>       tag to create (default '" << tag << "', use NC... to get non colliding, counted fingerprints\n";
  cerr << " -S <fname>     fetch smiles from <fname>\n";
  cerr << " -s             insert a dummy smiles with each fingerprint\n";
  cerr << " -f             work as a TDT filter\n";
  cerr << " -P <fname>     when working as TDT filter, fetch descriptor values from <fname>\n";
  cerr << " -z             strip leading zero's from identifiers\n";
//cerr << " -t <int>       numeric values above <int> will be set to 1\n";   // not sure what this was designed to do
  cerr << " -k             input file does NOT have a header record\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}


static int
fill_array_with_bit_values(const const_IWSubstring & buffer,
                           int & i,
                           int * numeric_value)
{
  set_vector(numeric_value, nbits, 0);

  const_IWSubstring token;

  for (int col = 0; buffer.nextword(token, i); col++)
  {
    if (1 == token.length())
    {
      char c = token[0];
      if (c >= '0' && c <= '9')
      {
        numeric_value[col] = c - '0';
        continue;
      }
    }

    if (missing_value == token)
      numeric_value[col] = 0;      // not sure this makes sense...
    else if (! token.numeric_value(numeric_value[col]))
    {
      cerr << "Invalid numeric '" << token << "' in '" << buffer << "'\n";
      return 0;
    }
    else if (numeric_value[col] < 0)
    {
      cerr << "Negative values not allowed\n";
      return 0;
    }
  }

  return 1;
}

static int
get_smiles(const IWString & id,
           IWString & smiles)
{
  IW_STL_Hash_Map_String::const_iterator f = id_to_smiles.find(id);

  if (f != id_to_smiles.end())
  {
    smiles = (*f).second;
    return 1;
  }

  if (! strip_leading_zeros_from_identifiers)
    return 0;

  if (! id.starts_with('0'))
    return 0;

  IWString tmp(id);
  tmp.remove_leading_chars('0');
  f = id_to_smiles.find(tmp);

  if (f == id_to_smiles.end())
    return 0;

  smiles = (*f).second;

  return 1;
}

static int
write_smiles (const IWString & id,
              IWString_and_File_Descriptor & output)
{
  IWString smiles;

  if (! get_smiles(id, smiles))
    return 0;

  output << smiles_tag << smiles << ">\n";
  return 1;
}

static int
descriptor_file_to_01_fingerprints_record (const const_IWSubstring & buffer,
                                           IWString_and_File_Descriptor & output)
{
  if (columns_in_input != buffer.nwords())
  {
    cerr << "Column count mismatch, got" << buffer.nwords() << " expected " << columns_in_input << endl;
    return 0;
  }

  IWString id;

  int i = 0;

  if (! buffer.nextword(id, i))
  {
    cerr << "Cannot extract identfier\n";
    return 0;
  }

  if (! fill_array_with_bit_values(buffer, i, numeric_value))
    return 0;

  if (insert_dummy_smiles.length())
    output << smiles_tag << insert_dummy_smiles << ">\n";
  else if (0 == id_to_smiles.size())
    ;
  else if (! write_smiles(id, output))
  {
    cerr << "No smiles for '" << id << "'\n";
    return 0;
  }

  output << identifier_tag << id << ">\n";

  IW_Bits_Base fp(nbits);

  fp.construct_from_array_of_ints(numeric_value, nbits);

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);

  output << tag << tmp << ">\n";

  output << "|\n";

  return 1;
}

static int
descriptor_file_to_01_fingerprints_main (iwstring_data_source & input,
                                         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    records_read++;

    if (! descriptor_file_to_01_fingerprints_record (buffer, output))
    {
      cerr << "Fatal error processing line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }

    if (output.length() > 32768)
      output.write_whole_blocks_shift_unwritten();
  }

  return output.good();
}

static int
descriptor_file_to_01_fingerprints (iwstring_data_source & input,
                                    IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  int c = buffer.nwords();

  if (c <= 0)
    return 0;

  if (0 == columns_in_input)   // first file processed, set nbits
  {
    columns_in_input = c;

    c--;
    if (0 == c % 8)
      nbits = c;
    else
      nbits = ( (c / 8) * 8 + 8);

    if (verbose)
      cerr << "Found " << c << " columns in input, nbits = " << nbits << endl;

    numeric_value = new int[nbits];
  }
  else if (c == columns_in_input)
    ;
  else
  {
    cerr << "Column count mismatch, got " << c << " expected " << columns_in_input << endl;
    return 0;
  }

  if (! input_file_has_a_header)
    input.push_record();

  return descriptor_file_to_01_fingerprints_main(input, output);
}

static int
descriptor_file_to_01_fingerprints (const char * fname,
     IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptor_file_to_01_fingerprints(input, output);
}

static int
descriptor_file_to_01_fingerprints_filter_record_sparse (const IWString & tag, 
                                        const int * c,
                                        IWString_and_File_Descriptor & output)
{

  Sparse_Fingerprint_Creator sfp;

  sfp.create_from_array_of_ints(c, nbits);

  IWString tmp;

  sfp.daylight_ascii_form_with_counts_encoded (tag, tmp);

  output << tmp << "\n";

  return 1;
}

static int
descriptor_file_to_01_fingerprints_filter_record_dense (const IWString & tag, 
                                        const int * c,
                                        IWString_and_File_Descriptor & output)
{
  IW_Bits_Base fp(nbits);

  fp.construct_from_array_of_ints(c, nbits);

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);

  output << tag << tmp << ">\n";

  return 1;
}

static int
descriptor_file_to_01_fingerprints_filter_record (const_IWSubstring & buffer,
                                           const ID_and_Vector & id_and_vector,
                                           IWString_and_File_Descriptor & output)
{
  buffer.remove_leading_chars(identifier_tag.length());
  buffer.chop();

  ID_and_Vector::const_iterator f = id_and_vector.find(buffer);

  if (f == id_and_vector.end())
  {
    cerr << "Canot find data for '" << buffer << "'\n";
    return 0;
  }

  if (tag.starts_with("FP"))
    return descriptor_file_to_01_fingerprints_filter_record_dense (tag, (*f).second, output);
  else
    return descriptor_file_to_01_fingerprints_filter_record_sparse (tag, (*f).second, output);
}

static int
descriptor_file_to_01_fingerprints_filter (iwstring_data_source & input,
                                           const ID_and_Vector & id_and_vector,
                                           IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    output << buffer << "\n";

    output.write_if_buffer_holds_more_than(32768);

    if (! buffer.starts_with(identifier_tag))
      continue;

    if (! descriptor_file_to_01_fingerprints_filter_record (buffer, id_and_vector, output))
    {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
descriptor_file_to_01_fingerprints_filter (const char * fname,
                                           const ID_and_Vector & id_and_vector,
                                           IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);   // should be stdin

  if (! input.good())
  {
    cerr << "Cannot open '" <<fname << "'\n";
    return 0;
  }

  return descriptor_file_to_01_fingerprints_filter (input, id_and_vector, output);
}

static int
fetch_id_and_smiles(const const_IWSubstring & buffer,
                    IW_STL_Hash_Map_String & id_to_smiles)
{
  if (buffer.nwords() < 2)
  {
    cerr << "No smiles info\n";
    return 0;
  }

  IWString smiles, id;
  int i = 0;

  buffer.nextword(smiles, i);
  buffer.nextword(id, i);

  if (id_to_smiles.contains(id))
  {
    cerr << "Warning, duplicate id '" << id << "', ignored\n";
    return 1;
  }

  id_to_smiles[id] = smiles;

  return id_to_smiles.size();
}

static int
read_smiles(iwstring_data_source & input,
            IW_STL_Hash_Map_String & id_to_smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! fetch_id_and_smiles(buffer, id_to_smiles))
    {
      cerr << "Cannot extract id and smiles from '" << buffer << "'\n";
      return 0;
    }
  }

  return id_to_smiles.size();
}

static int
read_descriptor_file_record (const const_IWSubstring & buffer,
                             ID_and_Vector & id_and_vector)
{
  if (columns_in_input != buffer.nwords())
  {
    cerr << "Column count mismatch, got" << buffer.nwords() << " expected " << columns_in_input << endl;
    return 0;
  }

  IWString id;

  int i = 0;

  if (! buffer.nextword(id, i))
  {
    cerr << "Cannot extract identfier\n";
    return 0;
  }

  if (strip_leading_zeros_from_identifiers)
    id.remove_leading_chars('0');

  if (! fill_array_with_bit_values(buffer, i, numeric_value))
    return 0;

  int * tmp = new int[nbits];

  for (int i = 0; i < nbits; i++)
  {
    if (numeric_value[i] > truncate_large_values)
      tmp[i] = truncate_large_values;
    else
      tmp[i] = numeric_value[i];
  }

  id_and_vector[id] = tmp;

  return 1;
}

static int
read_descriptors(iwstring_data_source & input,
                 ID_and_Vector & id_and_vector)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record from descriptor file\n";
    return 0;
  }

  columns_in_input = buffer.nwords();

  nbits = columns_in_input - 1;

  if (nbits < 1)
  {
    cerr << "No descriptors in header record '" << buffer << "'\n";
    return 0;
  }

  if (0 != nbits % 8)
    nbits = (nbits / 8 + 1) * 8;

  numeric_value = new int[nbits];

  while (input.next_record (buffer))
  {
    if (! read_descriptor_file_record (buffer, id_and_vector))
    {
      cerr << "Fatal error processing descriptor file record '" << buffer << "'\n";
      return 0;
    }
  }

  return id_and_vector.size();
}

static int
read_smiles(const const_IWSubstring & fname,
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

static int
read_descriptors (const const_IWSubstring & fname,
                  ID_and_Vector & id_and_vector)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open descriptor file '" << fname << "'\n";
    return 0;
  }

  return read_descriptors (input, id_and_vector);
}

static int
descriptor_file_to_01_fingerprints (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:S:st:zkfP:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('F'))
  {
    cl.value('F', tag);

    if (verbose)
      cerr << "Fingerprints produced with tag '" << tag << "'\n";

    if (! tag.ends_with('<'))
      tag << '<';

    if (tag.starts_with("FP"))
      truncate_large_values = 1;
  }

  if (cl.option_present('s') && cl.option_present('S'))
  {
    cerr << "Sorry, the -s (dummy smiles) and -S (file of smiles) options cannot both be used\n";
    usage(4);
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (! read_smiles (s, id_to_smiles))
    {
      cerr << "Cannot read identifier->smiles relationship from '" << s << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Read " << id_to_smiles.size() << " id->smiles relationships from '" << s << "'\n";
  }

  if (cl.option_present('s'))
  {
    insert_dummy_smiles = "C";
    if (verbose)
      cerr << "Will insert a dummy smiles\n";
  }

  if (cl.option_present('z'))
  {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Will strip leading zero's from identifiers\n";
  }

  if (cl.option_present('k'))
  {
    input_file_has_a_header = 0;

    if (verbose)
      cerr << "Input file does not have a header record\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', truncate_large_values) || truncate_large_values < 0)
    {
      cerr << "Invalid threshold (-t)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Values above " << truncate_large_values << " will be set to 1\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('f'))
  {
    if (cl.number_elements() > 1)
    {
      cerr << "When working as a tdt filter, can be only one argument\n";
      return 2;
    }

    work_as_tdt_filter = 1;

    if (verbose)
      cerr << "Will work as a TDT filter\n";

    if (! cl.option_present('P'))
    {
      cerr << "When working as a TDT filter must specify descriptor file via the -P option\n";
      usage(2);
    }
    if (cl.option_present('s') || cl.option_present('S'))
    {
      cerr << "The smiles options (-s or -S) options and -f (work as tdt filter) do not make sense together\n";
      return 2;
    }

    const_IWSubstring p = cl.string_value('P');

    ID_and_Vector id_and_vector;

    if (! read_descriptors (p, id_and_vector))
    {
      cerr << "Cannot read id's and numeric values from '" << p << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Read " << id_and_vector.size() << " id and numeric values from '" << p << "'\n";

    if (! descriptor_file_to_01_fingerprints_filter (cl[0], id_and_vector, output))
    {
      cerr << "Fatal error processing filter\n";
      rc = 1;
    }
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! descriptor_file_to_01_fingerprints(cl[i], output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose)
  {
    cerr << "Wrote " << records_read << " records\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_to_01_fingerprints(argc, argv);

  return rc;
}
