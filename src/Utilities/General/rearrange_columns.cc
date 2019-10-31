#include <stdlib.h>

using namespace std;

/*
  Rearrange the columns in one file to be the same as in another
*/

#include "iw_stl_hash_map.h"

#include "cmdline.h"
#include "iwstring_data_source.h"

#include "valid_descriptor_name.h"

static int verbose = 0;

static int ignore_case_when_comparing_descriptor_names = 0;

static int is_descriptor_file = 1;

static int check_for_valid_descriptors = 0;   // silly feature no longer needed

/*
  What do we do when descriptors are not in the reference file, but
  do appear in the input file(s)
*/

static int delete_columns_not_in_reference_file = 0;

static int fail_unless_all_descriptors_present = 0;

/*
  What do we do when descriptors are in the reference file, but not
  in any of the input files
*/

static IWString echo_columns_in_likefile_but_not_in_fromfile;

static int only_write_records_when_all_files_have_data = 0;

static int only_write_descriptors_present_in_all_files = 0;

static IWString_and_File_Descriptor  tream_for_discarded_identifiers;

static int identifiers_not_in_every_file = 0;

static int records_read = 0;
static int records_written = 0;

static IWString_and_File_Descriptor stream_for_discarded_identifiers;

/*
  sometimes people want to merge
  000010 with 10
*/

static int trim_leading_zeros_from_identifiers = 0;

static Valid_Descriptor_Name valid_descriptor_name;

static int function_as_cat_if_just_one_file_present = 0;

static int dash_p_option_present = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;

  cerr << "Rearranges the columns in one descriptor file to look like another descriptor file\n";
  cerr << " rearrange_columns -l reference_file <other options> <file1> <file2> ... > rearranged_file\n";
  cerr << " -l <file>      the reference descriptor file\n";
  cerr << " -o             ignore case when comparing descriptor names\n";
  cerr << " -c <column>    identifier column (default 1)\n";
  cerr << " -t             delete columns not in the reference file\n";
  cerr << " -f <string>    include columns in reference file that aren't in <file1>...\n";
  cerr << " -I row         only write records when every input file contains data\n";
  cerr << " -I col         only write descriptors where every reference file contains the descriptor\n";
  cerr << " -K <fname>     write identifiers discarded by -I option to <fname>\n";
  cerr << " -z             trim leading zero's from identifiers\n";
  cerr << " -J             input file(s) are NOT descriptor files\n";
  cerr << " -y             do the same as 'cat' if just one file on the command line\n";
  cerr << " -p             do NOT automatically apply a default -f value when multiple files present\n";
  cerr << " -b             fail unless all descriptors in the likefile are present\n";
  display_standard_valid_descriptor_name_options (cerr, 'm', 'V');
  cerr << " -C             allow any form of descriptor name\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
do_cat (iwstring_data_source & input,
        int output)
{
  IWString output_buffer;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    output_buffer << buffer << '\n';

    if (output_buffer.length () > 32768)
    {
      if (output_buffer.length () != output_buffer.write_whole_blocks_shift_unwritten (output))
        return 0;
    }
  }

  if (output_buffer.length ())
    output_buffer.write (output);

  return 1;
}

static int
do_cat (const char * fname,
        int output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return do_cat (input, output);
}

template <typename T>
void
do_trim_leading_zeros_from_identifiers (T & buffer)
{
  int nzero = 0;
  for (int i = 0; i < buffer.length (); i++)
  {
    if ('0' != buffer[i])
      break;

    nzero++;
  }

  if (0 == nzero)
    return;

  if (buffer.length () == nzero)
  {
    cerr << "Cannot strip zeros from '" << buffer << "'\n";
    return;
  }

  if (nzero)
    buffer.remove_leading_chars (nzero);

  return;
}

template void do_trim_leading_zeros_from_identifiers (IWString &);
template void do_trim_leading_zeros_from_identifiers (const_IWSubstring &);

class Input_File
{
  private:
    IWString _fname;

    IWString _header;

    int * _column_cross_reference;

    int _columns_in_file;

    int _identifier_column;

    iwstring_data_source _input;

    IW_STL_Hash_Map_off_t _offset;

//  private functions

    void _establish_column_cross_reference (const IWString & descriptor,
                                            int col,
                                            IW_STL_Hash_Map_int & global_cross_reference);
    int _echo (const const_IWSubstring & buffer,
                   const_IWSubstring * output_tokens) const;

  public:
    Input_File ();
    ~Input_File ();

    const IWString & fname () const { return _fname;}

    void set_identifier_column (int c) { _identifier_column = c;}

    int open_and_read_header (const char * n);

    int position_just_after_header_record ();

    int gather_descriptor_names (resizable_array_p<IWString> &, IW_STL_Hash_Map_int &) const;

    int print_cross_reference (ostream &) const;

    int establish_offsets ();

//  The first file will be asked to cycle through its input

    int next_record (const_IWSubstring &);

    int echo (const IWString & id, const_IWSubstring * output_buffer, int & missing);
    int echo_next_record_in_file (IWString &, const_IWSubstring * output_buffer, int &);

    int establish_column_cross_reference (IW_STL_Hash_Map_int &);
    int adjust_cross_reference_for_removed_descriptors (IW_STL_Hash_Map_int & global_cross_reference);
};

Input_File::Input_File ()
{
  _column_cross_reference = NULL;

  _columns_in_file = 0;

  if (is_descriptor_file)
    _identifier_column = 0;
  else
    _identifier_column = -1;

  return;
}

Input_File::~Input_File ()
{
  if (NULL != _column_cross_reference)
    delete _column_cross_reference;

  return;
}

int
Input_File::open_and_read_header (const char * n)
{
  assert (0 == _fname.length ());

  if (! _input.open (n))
  {
    cerr << "Input_File::open: cannot open '" << n << "'\n";
    return 0;
  }

  if (! _input.next_record (_header))
  {
    cerr << "Input_File::open_and_read_header: cannot read header\n";
    return 0;
  }

  _columns_in_file = _header.nwords ();

  _fname = n;

  return 1;
}

int
Input_File::position_just_after_header_record ()
{
  if (! _input.seekg (0))
  {
    cerr << "Input_File::position_just_after_header_record: cannot seek back to beginning of file\n";
    return 0;
  }

  const_IWSubstring buffer;

  return _input.next_record (buffer);
}

int
Input_File::gather_descriptor_names (resizable_array_p<IWString> & descriptors,
                                     IW_STL_Hash_Map_int & descriptor_names) const
{
  assert (_header.length () > 0);

  IW_STL_Hash_Map_int descriptors_found_this_file;

  int i = 0;
  IWString descriptor;
  int col = 0;

  while (_header.nextword (descriptor, i))
  {
    if (! check_for_valid_descriptors)
      ;
    else if (! valid_descriptor_name.valid (descriptor))
    {
      cerr << "Input_File::gather_descriptor_names: invalid descriptor name '" << descriptor << "'\n";
      return 0;
    }

//  cerr << "Examinining descriptor '" << descriptor << "' in column " << col << " id in " << _identifier_column << endl;

    if (col == _identifier_column)
    {
      col++;
      continue;
    }

    col++;

    if (ignore_case_when_comparing_descriptor_names)
      descriptor.to_lowercase ();

    if (descriptors_found_this_file.contains (descriptor))
    {
      cerr << "Input_File::establish_column_cross_reference: duplicate descriptor '" << descriptor << "'\n";
//    if (check_for_valid_descriptors)
//      return 0;
      continue;
//    return 0;
    }

    descriptors_found_this_file[descriptor] = 1;     // should use a hash set since we don't use the value '1' - just don't want to have to compile two headers

    if (descriptor_names.contains (descriptor))
    {
      descriptor_names[descriptor]++;
      continue;
    }

    IWString * tmp = new IWString (descriptor);

    descriptors.add (tmp);

    descriptor_names[descriptor] = 1;
  }

  return 1;
}

int
Input_File::establish_offsets ()
{
  assert (_fname.length () > 0);

  const_IWSubstring buffer;

  off_t offset = _input.tellg ();

  while (_input.next_record (buffer))
  {
    if (_identifier_column > 0)
      buffer.remove_leading_words (_identifier_column);

    buffer.truncate_at_first (' ');

    if (trim_leading_zeros_from_identifiers)
      do_trim_leading_zeros_from_identifiers (buffer);

    _offset[buffer] = offset;

    offset = _input.tellg ();
  }

  if (verbose)
    cerr << "Input_File::establish_offsets: file '" << _fname << "' contains " << _offset.size () << " offsets\n";

  return 1;
}


/*
  A common operation is to assign a unique number to a descriptor name
  That unique number is the current size of the hash
*/

static void
assign_next_number (IW_STL_Hash_Map_int & hash,
                    const IWString & zindex)
{
  int s = hash.size ();

  hash[zindex] = s;

  return;
}

void
Input_File::_establish_column_cross_reference (const IWString & descriptor,
                                               int col,
                                               IW_STL_Hash_Map_int & global_cross_reference)
{
  if (col == _identifier_column)
  {
    _column_cross_reference[col] = -1;
    return;
  }

  if (global_cross_reference.contains (descriptor))
  {
#ifdef DEBUG_ESTABLISH_COLUMN_CROSS_REFERENCE
    cerr << "Global cross reference puts '" << descriptor << "' our column " << col << " in column " << global_cross_reference[descriptor] << endl;
#endif

    _column_cross_reference[col] = global_cross_reference[descriptor];
    return;
  }

  if (delete_columns_not_in_reference_file)
  {
    _column_cross_reference[col] = -1;
    return;
  }

#ifdef DEBUG_ESTABLISH_COLUMN_CROSS_REFERENCE
  cerr << "New descriptor '" << descriptor << "'\n";
#endif

//  need to add an extra column

  assign_next_number (global_cross_reference, descriptor);

  return;
}

int
Input_File::establish_column_cross_reference (IW_STL_Hash_Map_int & global_cross_reference)
{
  assert (_header.length () > 0);
  assert (_columns_in_file > 0);
  assert (NULL == _column_cross_reference);

  _column_cross_reference = new int[_columns_in_file];

  int i = 0;
  IWString descriptor;
  int col = 0;

  while (_header.nextword (descriptor, i))
  {
    if (ignore_case_when_comparing_descriptor_names)
      descriptor.to_lowercase ();

    _establish_column_cross_reference (descriptor, col, global_cross_reference);

    col++;
  }

  return 1;
}


int
Input_File::print_cross_reference (ostream & os) const
{
  os << "Cross references for '" << _fname << "', " << _columns_in_file << " columns\n";

  for (int i = 0; i < _columns_in_file; i++)
  {
    os << "  column " << i;
    const_IWSubstring dname;
    _header.word (i, dname);
    cerr << ' ' << dname;
    cerr << " goes to column " << _column_cross_reference[i] << endl;
  }

  return os.good ();
}

/*
  Echo the data for a given identifier
*/

int
Input_File::echo (const IWString & id,
                  const_IWSubstring * output_buffer,
                  int & missing)
{
  IW_STL_Hash_Map_off_t::const_iterator f = _offset.find (id);

  if (f == _offset.end ())
  {
    if (verbose > 1)
      cerr << "Input_File::echo: no offset for '" << id << "' in '" << _fname << "'\n";

    missing = 1;

    return 1;
  }

  off_t o = (*f).second;

  if (! _input.seekg (o))
  {
    cerr << "Input_File::echo: cannot seek to " << o << " for '" << id << "' in '" << _fname << "'\n";
    return 0;
  }

//cerr << "Seeked to " << o << " for '" << id << "'\n";

  const_IWSubstring buffer;

  if (! _input.next_record (buffer))
  {
    cerr << "Input_File::echo: cannot read data at " << o << " for '" << id << "' in '" << _fname << "'\n";
    return 0;
  }

//cerr << "Fetched raw data '" << buffer << "'\n";

  return _echo (buffer, output_buffer);
}

/*
  
*/

int
Input_File::echo_next_record_in_file (IWString & id,
                                      const_IWSubstring * output_tokens,
                                      int & fatal)
{
  const_IWSubstring buffer;

  if (! _input.next_record (buffer))     // normal eof
  {
    fatal = 0;
    return 0;
  }

  if (_identifier_column >= 0)
    buffer.word (_identifier_column, id);

  if (!_echo (buffer, output_tokens))
  {
    fatal = 1;
    return 0;
  }

  return 1;
}

int
Input_File::_echo (const const_IWSubstring & buffer,
                   const_IWSubstring * output_tokens) const
{
  int i = 0;
  const_IWSubstring token;
  int col = 0;
  while (buffer.nextword (token, i))
  {
    if (col == _identifier_column)
    {
      col++;
      continue;
    }

//  cerr << "In column " << col << " got '" << token << "'\n";
    if (col >= _columns_in_file)
    {
      cerr << "Input_File::echo:invalid column count, only " << _columns_in_file << " columns\n";
      return 0;
    }

    int j = _column_cross_reference[col];

//  cerr << "Column " << col << " goes to column " << j << endl;

    if (j >= 0)
      output_tokens[j] = token;

    col++;
  }

  return 1;
}

static Input_File * input_file = NULL;
static int number_input_files = 0;

/*
  We have a request to only write descriptors present in all likefiles
*/

static void
discard_descriptors_not_in_all_files (resizable_array_p<IWString> & likefile_descriptors,
                       IW_STL_Hash_Map_int & likefile_descriptors_column_numbers,
                       IW_STL_Hash_Map_int & likefile_descriptor_count,
                       int number_likefiles)
{
  int j = 0;     // new column number to be assigned
  for (int i = 0; i < likefile_descriptors.number_elements (); i++)
  {
    const IWString & d = *(likefile_descriptors[i]);

    if (number_likefiles == likefile_descriptor_count[d])
    {
      likefile_descriptors_column_numbers[d] = j;
      j++;
    }
    else
    {
      if (verbose)
        cerr << "Descriptor '" << d << "' in only " << likefile_descriptor_count[d] << " reference files, discarded\n";

      likefile_descriptors_column_numbers.erase (d);
      likefile_descriptors.remove_item (i);
      i--;
    }
  }

  return;
}

static int
establish_column_cross_reference (IW_STL_Hash_Map_int & global_cross_reference)
{
  for (int i = 0; i < number_input_files; i++)
  {
    if (! input_file[i].establish_column_cross_reference (global_cross_reference))
    {
      cerr << "Cannot establish column cross reference for '" << input_file[i].fname () << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  Just read all the unique descriptors in the header
*/

static int
collect_likefile_descriptor_names (const const_IWSubstring & header,
                      resizable_array_p<IWString> & likefile_descriptors,
                      IW_STL_Hash_Map_int & likefile_descriptors_column_numbers,
                      IW_STL_Hash_Map_int & likefile_descriptor_count)
{
  IWString token;
  int i = 0;

  IW_STL_Hash_Map_int in_this_file;    // should be a hash_set

  while (header.nextword (token, i))
  {
    if (! check_for_valid_descriptors)
      ;
    else if (! valid_descriptor_name.valid (token))
      cerr << "Invalid descriptor name in reference file '" << token << "', ignored!\n";

    if (ignore_case_when_comparing_descriptor_names)
      token.to_lowercase ();

    if (in_this_file.contains (token))
    {
      cerr << "Warning, duplicate descriptor in reference file '" << token << "', ignored\n";
      continue;
    }

    in_this_file[token] = 1;

    if (likefile_descriptors_column_numbers.contains (token))   // great, exact match
    {
      likefile_descriptor_count[token]++;
    }
    else
    {
      IWString * tmp = new IWString (token);
      likefile_descriptors.add (tmp);

      assign_next_number (likefile_descriptors_column_numbers, token);

      likefile_descriptor_count[token] = 1;
    }
  }

  return 1;
}

static int
read_likefile_header (IWString & fname,
                      resizable_array_p<IWString> & likefile_descriptors,
                      IW_STL_Hash_Map_int & likefile_descriptors_column_numbers,
                      IW_STL_Hash_Map_int & likefile_descriptor_count)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open reference file '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring header;

  if (! input.next_record (header))
  {
    cerr << "read_likefile_header: cannot read header\n";
    return 0;
  }

  if (is_descriptor_file)
    header.remove_leading_words (1);      // get rid of the NAME or ID field

  return collect_likefile_descriptor_names (header, likefile_descriptors, likefile_descriptors_column_numbers, likefile_descriptor_count);
}

/*
  If there is no likefile, we assume that the input files themselves constitute
  the likefiles
*/

static int
read_likefile_headers (Command_Line & cl,
                       resizable_array_p<IWString> & likefile_descriptors,
                       IW_STL_Hash_Map_int & likefile_descriptors_column_numbers,
                       IW_STL_Hash_Map_int & likefile_descriptor_count)
{
  assert (number_input_files > 1);
  assert (cl.number_elements () == number_input_files);

  for (int i = 0; i < number_input_files; i++)
  {
    IWString f = cl[i];

    if (! read_likefile_header (f, likefile_descriptors, likefile_descriptors_column_numbers, likefile_descriptor_count))
    {
      cerr << "Fatal error reading reference file '" << f << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_likefile_headers (Command_Line & cl,
                       char flag,
                       resizable_array_p<IWString> & likefile_descriptors,
                       IW_STL_Hash_Map_int & likefile_descriptors_column_numbers,
                       IW_STL_Hash_Map_int & likefile_descriptor_count)
{

  int i = 0;
  IWString f;

  while (cl.value (flag, f, i++))
  {
    if (! read_likefile_header (f, likefile_descriptors, likefile_descriptors_column_numbers, likefile_descriptor_count))
    {
      cerr << "Fatal error reading reference file '" << f << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  Keep at file scope for efficiency
*/

static IWString output_buffer;

static int
write_output (const IWString & id,
              const const_IWSubstring * token,
              int columns_in_output,
              IWString_and_File_Descriptor & output)
{
  output << id;

  for (int i = 0; i < columns_in_output; i++)
  {
    const const_IWSubstring & ti = token[i];

    if (ti.length())
      output.append_with_spacer(ti);
    else if (! dash_p_option_present)
      cerr << "Huh, column " << i << " is blank\n";
    else if (echo_columns_in_likefile_but_not_in_fromfile.length ())
      output.append_with_spacer(echo_columns_in_likefile_but_not_in_fromfile);
    else
      output << " .";
  }

//cerr << "Output buffer contains " << output_buffer.nwords () << " words, columns_in_output = " << columns_in_output << endl;
//cerr << output_buffer << endl;

  output += '\n';

  records_written++;

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}


static void
fill_output_with_missing_values (const_IWSubstring * output_buffer,
                                 int ncol,
                                 const IWString & fill_value)
{
  for (int i = 0; i < ncol; i++)
  {
    output_buffer[i] = fill_value;
  }

  return;
}

static void
handle_missing_identifiers (const IWString & id,
                            const resizable_array<int> & missing_from_file)
{
  if (stream_for_discarded_identifiers.is_open ())
  {
    stream_for_discarded_identifiers << id << " missing from";
    for (int i = 0; i < missing_from_file.number_elements (); i++)
    {
      const Input_File & f = input_file[missing_from_file[i]];

      stream_for_discarded_identifiers << ' ' << f.fname ();
    }

    stream_for_discarded_identifiers << '\n';

    stream_for_discarded_identifiers.write_if_buffer_holds_more_than(32768);
  }

  return;
}

static int
rearrange_columns (IW_STL_Hash_Map_int & likefile_column_hash,
                   const_IWSubstring * output_tokens,
                   int columns_in_output,
                   int zfile,
                   IW_STL_Hash_Map_int & identifiers_processed,
                   IWString_and_File_Descriptor & output)
{
  Input_File & fi = input_file[zfile];

  if (! fi.position_just_after_header_record ())
    return 0;

// we want to be able to report the file names that are missing certain identifiers. Keep the
// list at this scope for efficiency

  resizable_array<int> missing_from_file;

  int fatal;    // scope here for efficiency

  while (1)
  {
    if (echo_columns_in_likefile_but_not_in_fromfile.length ())
      fill_output_with_missing_values (output_tokens, columns_in_output, echo_columns_in_likefile_but_not_in_fromfile);

    IWString id;

    if (! fi.echo_next_record_in_file (id, output_tokens, fatal))
    {
      if (fatal)
      {
        cerr << "Failure echo next record for '" << id << "'\n";
        return 0;
      }

      return 1;
    }

    if (trim_leading_zeros_from_identifiers)
      do_trim_leading_zeros_from_identifiers (id);

    if (identifiers_processed.contains (id))
      continue;

// If we are processing anything but the last of multiple files, we
// need to keep track of which identifiers we have done.

    if (zfile < number_input_files - 1)
      identifiers_processed[id] = 1;

    records_read++;

    for (int i = zfile + 1; i < number_input_files; i++)
    {
      int missing = 0;
      if (! input_file[i].echo (id, output_tokens, missing))
      {
        cerr << "Fatal error processing '" << id << "'\n";
        return 0;
      }

      if (missing)
        missing_from_file.add (i);
    }

    if (missing_from_file.number_elements ())
    {
      handle_missing_identifiers (id, missing_from_file);

      missing_from_file.resize_keep_storage (0);
      if (only_write_records_when_all_files_have_data)
        continue;
    }

    if (! write_output (id, output_tokens, columns_in_output, output))
    {
      cerr << "write_output failed\n";
      return 0;
    }
  }

  return 1;
}

int
rearrange_columns (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vl:tc:of:I:K:zm:V:JCypb");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised options encountered\n";
    usage (5);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('z'))
  {
    trim_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Leading zero's will be trimmed from identifiers\n";
  }

  if (cl.option_present ('o'))     // process this before the -l option
  {
    ignore_case_when_comparing_descriptor_names = 1;

    if (verbose)
      cerr << "Will ignore case when comparing descriptor names\n";
  }

  if (cl.option_present ('C'))
  {
    if (cl.option_present ('m') || cl.option_present ('V'))
    {
      cerr << "The -m and -V options are mutually inconsistent with the -C option\n";
      usage (4);
    }

    check_for_valid_descriptors = 0;

    if (verbose)
      cerr << "All checking for valid descriptor names suppressed\n";
  }

  if (! valid_descriptor_name.construct_from_command_line (cl, 'm', 'V', verbose))
  {
    cerr << "Cannot initialise valid descriptor name conditions (-m and -V options)\n";
    usage (7);
  }

  if (cl.option_present ('J'))
  {
    is_descriptor_file = 0;

    if (verbose)
      cerr << "Not descriptor file(s)\n";
  }

  resizable_array<int> identifier_columns;

  if (cl.option_present ('c'))
  {
    int i = 0;
    const_IWSubstring c;
    while (cl.value ('c', c, i++))
    {
      int col;
      if (! c.numeric_value (col) || col < 1)
      {
        cerr << "Invalid identifier column specification '" << c << "'\n";
        usage (7);
      }

      if (verbose)
        cerr << "Identifiers in column " << col << endl;

      col--;

      identifier_columns.add (col);
    }
  }
  else if (! is_descriptor_file)
    ;
  else
    identifier_columns.add (0);

  if (cl.option_present ('t'))
  {
    delete_columns_not_in_reference_file = 1;

    if (verbose)
      cerr << "Columns not in the reference file will be deleted\n";
  }

  if (cl.option_present('b'))
  {
    fail_unless_all_descriptors_present = 1;

    if (verbose)
      cerr << "Will fail unless all descriptors present\n";
  }

  if (cl.option_present ('f'))
  {
    echo_columns_in_likefile_but_not_in_fromfile = cl.string_value ('f');

    if (verbose)
      cerr << "Columns in reference file but not in fromfile echo'd as '" << echo_columns_in_likefile_but_not_in_fromfile << endl;
  }

  if (cl.option_present ('I'))
  {
    const_IWSubstring ival;
    int i = 0;
    while (cl.value ('I', ival, i++))
    {
      if (ival.starts_with ("row"))
      {
        only_write_records_when_all_files_have_data = 1;

        if (verbose)
          cerr << "Will only write records when each input file has data for the identifier\n";
      }
      else if (ival.starts_with ("col"))
      {
        only_write_descriptors_present_in_all_files = 1;
        delete_columns_not_in_reference_file = 1;

        if (verbose)
          cerr << "Will only write descriptors present in all files\n";
      }
      else
      {
        cerr << "Unrecognised -I qualifier '" << ival << "'\n";
        usage (11);
      }
    }
  }

  number_input_files = cl.number_elements ();

  if (0 == number_input_files)
  {
    cerr << "Insufficient arguments\n";
    usage (5);
  }

  if (cl.option_present ('K'))
  {
    const char * fname = cl.option_value ('K');

    if (! stream_for_discarded_identifiers.open (fname))
    {
      cerr << "Cannot open stream for discarded identifiers '" << fname << "'\n";
      return 3;
    }

    only_write_records_when_all_files_have_data = 1;

    if (verbose)
      cerr << "Will only write complete records, discarded identifiers written to '" << fname << "'\n";
  }

  if (cl.option_present ('y'))
  {
    function_as_cat_if_just_one_file_present = 1;

    if (verbose)
      cerr << "Will do the same thing as 'cat' if just one file present\n";
  }

// If they want data from all input files, they need a missing value fill string

  if (cl.option_present('p'))
    dash_p_option_present = 1;
  else if (number_input_files > 1 && 
      ! only_write_records_when_all_files_have_data &&
      0 == echo_columns_in_likefile_but_not_in_fromfile.length ())
  {
    echo_columns_in_likefile_but_not_in_fromfile = ".";
  }

  resizable_array_p<IWString> likefile_descriptors;
  IW_STL_Hash_Map_int likefile_descriptors_column_numbers;
  IW_STL_Hash_Map_int likefile_descriptor_count;

  int number_likefiles = cl.option_count ('l');

  if (number_likefiles)
  {
    if (! read_likefile_headers (cl, 'l', likefile_descriptors, likefile_descriptors_column_numbers, likefile_descriptor_count))
    {
      cerr << "Cannot read reference file column header(s)\n";
      return 12;
    }
  }
  else
  {
    if (number_input_files > 1)
      ;
    else if (function_as_cat_if_just_one_file_present)
    {
      if (verbose)
        cerr << "Just one input file present, same as 'cat'\n";

      do_cat (cl[0], 1);

      return 0;
    }
    else
    {
      cerr << "Without a reference file (-l option), you must have > 1 input file\n";
      usage (7);
    }

    number_likefiles = number_input_files;

    if (! read_likefile_headers (cl, likefile_descriptors, likefile_descriptors_column_numbers, likefile_descriptor_count))
    {
      cerr << "Cannot read reference file column header(s)\n";
      return 13;
    }
  }

  if (verbose)
    cerr << number_likefiles << " like files contain " << likefile_descriptors.number_elements () << " descriptors, " << likefile_descriptors_column_numbers.size () << " unique\n";

  if (1 == number_input_files && only_write_records_when_all_files_have_data)
    cerr << "The -a option doesn't do anything when there is only one input file\n";

  if (only_write_descriptors_present_in_all_files && 1 == number_likefiles)
    cerr << "Only writing descriptors present in all files, but only one reference file\n";

  if (only_write_descriptors_present_in_all_files && number_likefiles > 1)
    discard_descriptors_not_in_all_files (likefile_descriptors, likefile_descriptors_column_numbers, likefile_descriptor_count, number_likefiles);

//#define ECHO_DESCRIPTORS
#ifdef ECHO_DESCRIPTORS
  cerr << "After trimming descriptors\n";
  for (IW_STL_Hash_Map_int::const_iterator i = likefile_descriptors_column_numbers.begin (); i != likefile_descriptors_column_numbers.end (); i++)
  {
    cerr << "Reference file descriptor '" << (*i).first << "' in column " << (*i).second << endl;
  }
#endif

  if (is_descriptor_file && identifier_columns.number_elements () < number_input_files)
    identifier_columns.extend (number_input_files, identifier_columns.last_item ());

  input_file = new Input_File[number_input_files];

  for (int i = 0; i < number_input_files; i++)
  {
    if (is_descriptor_file)
      input_file[i].set_identifier_column (identifier_columns[i]);

    if (! input_file[i].open_and_read_header (cl[i]))
    {
      cerr << "Cannot initialise input file '" << cl[i] << "'\n";
      return i + 1;
    }
  }

// We need to get the unique descriptor names out, in order

  resizable_array_p<IWString> input_file_descriptors;
  IW_STL_Hash_Map_int input_file_descriptor_hash;

  for (int i = 0; i < number_input_files; i++)
  {
    if (! input_file[i].gather_descriptor_names (input_file_descriptors, input_file_descriptor_hash))
    {
      cerr << "Input file " << i << " contains an invalid header\n";
      return i + 1;
    }
  }

#ifdef ECHO_DESCRIPTORS
  for (IW_STL_Hash_Map_int::const_iterator i = input_file_descriptor_hash.begin (); i != input_file_descriptor_hash.end (); i++)
  {
    cerr << "Input file descriptor '" << (*i).first << "' found " << (*i).second << " times\n";
  }
#endif

// At this stage, we know the descriptors in the likefiles,
// and the descriptors in the input files
// Figure out all the cross referencing...

  IW_STL_Hash_Map_int column_cross_reference;

  IWString header;
  header.resize ((likefile_descriptors.number_elements () + 1) * 9);
  if (is_descriptor_file)
    header = "Name";

  for (int i = 0; i < likefile_descriptors.number_elements (); i++)
  {
    const IWString & d = *(likefile_descriptors[i]);
//  cerr << "Examining descriptor '" << d << "' seen in inputs? " << input_file_descriptor_hash.contains(d) << " echo '" << echo_columns_in_likefile_but_not_in_fromfile << "'\n";

    if (input_file_descriptor_hash.contains (d) ||        // the easy case, in both likefile and input files
        echo_columns_in_likefile_but_not_in_fromfile.length ())
    {
      header.append_with_spacer (d);

      assign_next_number (column_cross_reference, d);
    }
  }

// Now deal with descriptors that are in the input files, but not in the likefile(s)

  if (! delete_columns_not_in_reference_file)
  {
    for (int i = 0; i < input_file_descriptors.number_elements (); i++)
    {
      const IWString & d = *(input_file_descriptors[i]);

      if (likefile_descriptors_column_numbers.contains (d))
        continue;

      header.append_with_spacer (d);

      assign_next_number (column_cross_reference, d);
    }
  }

  if (verbose)
  {
    IW_STL_Hash_Map_int::const_iterator i;
    for (i = column_cross_reference.begin (); i != column_cross_reference.end (); ++i)
    {
      cerr << "Descriptor '" << (*i).first << "' goes in column " << (*i).second << endl;
    }
  }

  if (! establish_column_cross_reference (column_cross_reference))
  {
    cerr << "Cannot build column cross reference\n";
    return 4;
  }

  if (verbose)
  {
    for (int i = 0; i < number_input_files; i++)
    {
      input_file[i].print_cross_reference (cerr);
    }
  }

  if (fail_unless_all_descriptors_present)
  {
    if (column_cross_reference.size() != likefile_descriptors_column_numbers.size())
    {
      cerr << "Likefile contains " << likefile_descriptors_column_numbers.size() << " descriptors, input file(s) contain " << column_cross_reference.size() << ", cannot continue(-b)\n";
      return 5;

    }
  }

  int columns_in_output = 0;

  for (IW_STL_Hash_Map_int::const_iterator i = column_cross_reference.begin (); i != column_cross_reference.end (); i++)
  {
    int col = (*i).second;

    if (col >= 0)
      columns_in_output++;
  }

  if (0 == verbose)
    ;
  else if (is_descriptor_file)
    cerr << "Output will contain " << columns_in_output << " descriptors\n";
  else
    cerr << "Output will contain " << columns_in_output << " columns\n";

  const_IWSubstring * output_tokens = new const_IWSubstring[columns_in_output];
  assert (NULL != output_tokens);

  IWString_and_File_Descriptor output(1);

  output << header << '\n';
  output.flush();

  if (number_input_files > 1)
  {
    for (int i = 0; i < number_input_files; i++)
    {
      if (! input_file[i].establish_offsets ())
      {
        cerr << "Yipes, file " << i << " cannot establish identifier offsets\n";
        return i + 1;
      }
    }
  }

  int rc = 0;

  IW_STL_Hash_Map_int identifiers_processed;

  for (int i = 0; i < number_input_files; i++)
  {
    if (! rearrange_columns (column_cross_reference, output_tokens, columns_in_output,
                             i, identifiers_processed, output))
    {
      cerr << "Fatal error processing file " << i << endl;
      rc = i + 1;
      break;
    }

    if (only_write_records_when_all_files_have_data)
      break;
  }

  output.close();

  delete [] output_tokens;

  if (verbose)
  {
    cerr << "Read " << records_read << " records";
    if (records_written != records_read)
      cerr << ", wrote " << records_written << endl;
    else
      cerr << endl;

    if (identifiers_not_in_every_file)
      cerr << identifiers_not_in_every_file << " records not written because of missing data\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = rearrange_columns (argc, argv);

  return rc;
}
