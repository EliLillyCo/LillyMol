#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <memory>
#include <fstream>

/*
  Compare 2 descriptor files
*/

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::cout;
using std::endl;

static int verbose = 0;

static int ignore_case_when_comparing_descriptor_names = 0;

/*
  Sept 03. We will be inserting prefixes at the start of
  descriptor names. Allow the trimming of a string from
  the front of descriptor names
*/

static IWString ignore_prefix;

static double tolerance = 0.0;

static int interpret_tolerance_as_absolute_difference = 0;

static double noise = 0.0;

static int differences_within_noise = 0;

static IWString missing_value('.');

static int records_read = 0;

static int records_differing = 0;

static int columns_in_file1 = 0;

static int same_row_order = 1;

static int assume_columns_in_correct_order = 0;

static IW_STL_Hash_Map_off_t identifier_offset_cross_reference;

static IWString * column_title = NULL;

static int * differences_in_column = NULL;

static int differences_less_than_tolerance = 0;

static extending_resizable_array<int> differences_per_record;

static int ignore_missing_value_mismatches = 0;

static int missing_value_mismatches = 0;

static int records_with_missing_value_mismatches = 0;

static int ignore_leading_zeros = 0;

static int remove_chars_after_underscore = 0;

// the number of columns int he first, but not in subsequent inputs.
static int missing_columns = 0;

//static ofstream stream_for_identifiers_with_different_data;
static IWString_and_File_Descriptor stream_for_identifiers_with_different_data;

/*
  Normally when an identifier is missing, we just abort
*/

static int ignore_missing_identifiers = 0;

static int report_missing_identifiers = 1;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Compare 2 descriptor files, doing comparisons by numeric differences\n";
  cerr << " -c             ignore case when comparing descriptor names\n";
  cerr << " -t <tol>       tolerance for floating point comparisons (ratio by default)\n";
  cerr << " -b             interpret the tolerance as an absolute difference\n";
  cerr << " -n <noise>     ignore differences if both values less than <noise>\n";
  cerr << " -M <char>      missing value string (default '" << missing_value << "')\n";
  cerr << " -m             ignore missing value mis-matches\n";
  cerr << " -x             files are not in the same row order\n";
  cerr << " -z             ignore leading zero's when comparing identifiers\n";
  cerr << " -u             strip identifiers at first _ character\n";
  cerr << " -G             ignore missing identifiers\n";
  cerr << " -D <file>      write identifiers with differences to <file>\n";
  cerr << " -q mid         no messages about missing identifiers\n";
  cerr << " -P <prefix>    strip <prefix> from all descriptor names\n";
  cerr << " -o             assume columns in same order - ignore descriptor names\n";
  cerr << " -e             stderr messages written to stdout - convenience\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static void
do_remove_chars_after_underscore(const_IWSubstring & s)
{
  int u = s.index('_');

//cerr << "INdex of underscore in '" << s << "' is " << u << '\n';
  if (u < 1)    // cannot do anything if first char is underscore
    return;

  s.iwtruncate(u);

  return;
}

static int
tokens_compare_after_exclusions(const const_IWSubstring & t1,
                                        const const_IWSubstring & t2)
{
  const_IWSubstring cpt1(t1);
  const_IWSubstring cpt2(t2);

  if (ignore_leading_zeros)
  {
    cpt1.remove_leading_chars('0');
    cpt2.remove_leading_chars('0');
  }

  if (remove_chars_after_underscore)
  {
    do_remove_chars_after_underscore(cpt1);
    do_remove_chars_after_underscore(cpt2);
  }

  return cpt1 == cpt2;
}

static int
setup_default_cross_reference(const IWString & header1,
                               const IWString & header2,
                               int * column_xref)
{
  int n = header2.nwords();

  for (int i = 0; i < n; i++)
  {
    column_xref[i] = i;
  }

  return 1;
}

/*
  
*/

static int
establish_column_xref(const IWString & file1_header,
                       const IWString & file2_header,
                       int * column_xref)
{
  IW_STL_Hash_Map_int file1_column_hash;
  IW_STL_Hash_Map_int file1_column_hash_lowercase;

  int i = 0;
  IWString token;
  int col = 0;

  int prefix_stripped = 0;

  while (file1_header.nextword(token, i))
  {
    if (ignore_prefix.length() > 0 && token.starts_with(ignore_prefix) && token.length() > ignore_prefix.length())
    {
      token.remove_leading_chars(ignore_prefix.length());
      prefix_stripped++;
    }

    file1_column_hash[token] = col;

    column_title[col] = token;

    if (ignore_case_when_comparing_descriptor_names)
    {
      token.to_lowercase();
      file1_column_hash_lowercase[token] = col;
    }

    col++;
  }

  if (verbose && ignore_prefix.length())
    cerr << prefix_stripped << " of " << (col - 1) << " columns in file 1 stripped prefix '" << ignore_prefix << "'\n";

  prefix_stripped = 0;

  int columns_matching = 0;
  // The number of colums that are in the same position in both files.
  int same_column = 0;

  i = 0;
  col = 0;
  for (;file2_header.nextword(token, i); ++col)
  {
    if (ignore_prefix.length() > 0 && token.starts_with(ignore_prefix) && token.length() > ignore_prefix.length())
    {
      token.remove_leading_chars(ignore_prefix.length());
      prefix_stripped++;
    }

    IW_STL_Hash_Map_int::const_iterator f = file1_column_hash.find(token);

    int to_col = -1;
    
    if (f != file1_column_hash.end()) {   // great, found it
      to_col = (*f).second;
    } else if (ignore_case_when_comparing_descriptor_names)
    {
      token.to_lowercase();
      f = file1_column_hash_lowercase.find(token);

      if (f != file1_column_hash_lowercase.end()) {
        to_col = (*f).second;
      }
    }

    if (to_col >= 0)
    {
      column_xref[col] = to_col;

      if (verbose > 1) {
        cerr << "Column " << (col + 1) << " '" << token << "' found in column " << (to_col + 1);
        if (col != to_col) {
          cerr << " *";
        }
        cerr << '\n';
      }

      columns_matching++;
      if (to_col == col) {
        ++same_column;
      }
    }
    else
    {
      column_xref[col] = -1;

      if (verbose)
        cerr << "Column " << (col + 1) << " '" << token << "' not in rhs. Ignored\n";
      ++missing_columns;
    }
  }

  if (verbose) {
    cerr << columns_matching << " columns match\n";
    cerr << same_column << " columns in the same position\n";
    if (ignore_prefix.length())
      cerr << prefix_stripped << " of " << (col - 1) << " columns in file 2 stripped prefix '" << ignore_prefix << "'\n";
  }

  if (columns_matching < 2)   // name will always match
  {
    cerr << "Yipes, no columns match\n";
    return 0;
  }

  return 1;
}

static int
tokenise(const const_IWSubstring & buffer,
          const_IWSubstring * token,
          const int * column_xref,
          int columns_in_input)
{
  int i = 0;
  const_IWSubstring v;
  int col = 0;

  while (buffer.nextword(v, i))
  {
    if (col >= columns_in_input)
    {
      cerr << "Too many columns in input, max is " << columns_in_input << '\n';
      return 0;
    }

    int to_col = column_xref[col];

    if (to_col >= 0)
    {
      token[to_col] = v;
    }

    col++;
  }

  return 1;
}

static void
report_mismatch(const const_IWSubstring & token1,
                 int col,
                 const const_IWSubstring & token2,
                 int & differences_this_record,
                 const const_IWSubstring & id)
{
  if (0 == differences_this_record)
  {
    cout << "ID: " << id << '\n';
    if (stream_for_identifiers_with_different_data.is_open())
      stream_for_identifiers_with_different_data << id << '\n';
  }

  differences_this_record++;

  cout << " mismatch in column " << (col + 1) << " '" << column_title[col] << "' '" << token1 << "' vs '" << token2 << "'\n";

  differences_in_column[col]++;
  
  return;
}

static int
jfilecompare(const const_IWSubstring & token1,
              int col,
              const const_IWSubstring & token2,
              int & differences_this_record,
              int & missing_value_mismatches_this_record,
              const const_IWSubstring & id)      // the id for the row
{
#ifdef DEBUG_JFILECOMPARE
  cerr << "Token1 '" << token1 << "'\n";
  cerr << "Token2 '" << token2 << "'\n";
#endif

  if (0 == token2.length())     // that column not present in RHS
    return 1;

  if (token1 == token2)
    return 1;

  if (missing_value == token1 || missing_value == token2)
  {
    if (ignore_missing_value_mismatches)
    {
      missing_value_mismatches_this_record++;
      return 1;
    }

    report_mismatch(token1, col, token2, differences_this_record, id);

    return 1;
  }

  if (0 == col)
  {
    if ((ignore_leading_zeros || remove_chars_after_underscore) && tokens_compare_after_exclusions(token1, token2))
      return 1;

    cout << "Identifier mismatch\n";
    report_mismatch(token1, col, token2, differences_this_record, id);
    cout << "Fatal\n";

    return 0;
  }

  if (0.0 == tolerance && 0.0 == noise)
  {
    cerr.flush();
//  cerr << "No tolerance or noise '" << token2 << "', id '" << id << "'\n";
    report_mismatch(token1, col, token2, differences_this_record, id);

    return 1;
  }

  double v1, v2;

  if (! token1.numeric_value(v1))
  {
    cerr << "INvalid numeric value in lhs : '" << token1 << "'\n";
    return 1;
  }

  if (! token2.numeric_value(v2))
  {
    cerr << "INvalid numeric value in rhs : '" << token2 << "'\n";
    return 1;
  }

  if (v1 == v2)
    return 1;

  if (interpret_tolerance_as_absolute_difference)
  {
    if (fabs(v1 - v2) <= tolerance)
    {
      differences_less_than_tolerance++;
      return 1;
    }
  }
  else
  {
    if (fabs( (v1 - v2) / (v1 + v2) * 0.5) <= tolerance)
    {
      differences_less_than_tolerance++;
      return 1;
    }
  }

  if (fabs(v1) <= noise && fabs(v2) <= noise)
  {
    differences_within_noise++;
    return 1;
  }

  report_mismatch(token1, col, token2, differences_this_record, id);

  return 1;
}

static int
jfilecompare(const const_IWSubstring & buffer1,
              int columns_in_buffer1,
              const const_IWSubstring & buffer2,
              int columns_in_buffer2,
              const int * column_xref,
              const_IWSubstring * token,
              int & differences_this_record)
{
  for (int i = 0; i < columns_in_buffer1; i++)
  {
    token[i].make_empty();
  }

  if (! tokenise(buffer2, token, column_xref, columns_in_buffer2))
  {
    cerr << "Yipes, could not tokenise RHS\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token1;
  int col = 0;

  differences_this_record = 0;

  int missing_value_mismatches_this_record = 0;

  const_IWSubstring id;     // we take this from buffer1

  while (buffer1.nextword(token1, i))
  {
    if (0 == col)
      id = token1;
    else if (col >= columns_in_buffer1)
    {
      cerr << "Too many columns in file1\n";
      return 0;
    }

//  cerr << "Fetching column " << col << '\n';
    if (! jfilecompare(token1, col, token[col], differences_this_record, missing_value_mismatches_this_record, id))
    {
      cerr << "Fatal error in comparison\n";
      return 0;
    }

    col++;
  }

  if (missing_value_mismatches_this_record)
  {
    missing_value_mismatches += missing_value_mismatches_this_record;
    records_with_missing_value_mismatches++;
  }

  return 1;
}

static int
fetch_second_record(const const_IWSubstring & buffer1,
                     iwstring_data_source & input2,
                     const_IWSubstring & buffer2,
                     int & eof)
{
  eof = 0;

  if (same_row_order)
  {
    if (! input2.next_record(buffer2))
    {
      cerr << "Premature EOF reading second file\n";
      eof = 1;
      return 0;
    }

    return 1;
  }

  const_IWSubstring id1 = buffer1;
  id1.truncate_at_first(' ');     // get just the first word

  if (ignore_leading_zeros)
    id1.remove_leading_chars('0');

  if (remove_chars_after_underscore)
      do_remove_chars_after_underscore(id1);

  if (identifier_offset_cross_reference.contains(id1))
  {
    off_t o = identifier_offset_cross_reference[id1];

    if (! input2.seekg(o))
    {
      cerr << "Yipes, cannot seek to " << o << " for '" << id1 << "'\n";
      return 0;
    }

    return input2.next_record(buffer2);
  }

  if (report_missing_identifiers)
    cerr << "Identifier '" << id1 << "' not found in second file\n";

  return ignore_missing_identifiers;
}

static int
jfilecompare(iwstring_data_source & input1,
              int columns_in_file1,
              iwstring_data_source & input2,
              int columns_in_file2,
              const int * column_xref,
              const_IWSubstring * token)
{
  const_IWSubstring buffer1;

  while (input1.next_record(buffer1))
  {
    const_IWSubstring buffer2;

    int eof;
    if (! fetch_second_record(buffer1, input2, buffer2, eof))
    {
      if (eof)
      {
        cerr << "Premature EOF reading rhs\n";
        return 0;
      }

      if (ignore_missing_identifiers)
        continue;

      return 0;
    }

    records_read++;

    int differences_this_record = 0;

    if (! jfilecompare(buffer1, columns_in_file1, buffer2, columns_in_file2, column_xref, token, differences_this_record))
    {
      cerr << "Error on line " << input1.lines_read() << '\n';
      return 0;
    }

    differences_per_record[differences_this_record]++;

    if (verbose > 1)
      cerr << differences_this_record << " differences\n";

    if (differences_this_record)
    {
      records_differing++;
    }
  }

  return 1;
}

static int
jfilecompare(iwstring_data_source & input1,
              const IWString & header1,
              int columns_in_file1,
              iwstring_data_source & input2,
              const IWString & header2,
              int columns_in_file2,
              int * column_xref,
              const_IWSubstring * token)
{
  if (assume_columns_in_correct_order)
    setup_default_cross_reference(header1, header2, column_xref);
  else if (! establish_column_xref(header1, header2, column_xref))
  {
    cerr << "Cannot establish column cross reference\n";
    cerr << header1 << '\n';
    cerr << header2 << '\n';

    return 0;
  }

  return jfilecompare(input1, columns_in_file1, input2, columns_in_file2, column_xref, token);
}

static int
jfilecompare(iwstring_data_source & input1,
              iwstring_data_source & input2)
{
  IWString header1;
  if (! input1.next_record(header1))
  {
    cerr << "Cannot fetch header record from lhs\n";
    return 0;
  }

  IWString header2;
  if (! input2.next_record(header2))
  {
    cerr << "Cannot fetch header record from rhs\n";
    return 0;
  }

  columns_in_file1 = header1.nwords();

  differences_per_record.resize(columns_in_file1);

  if (verbose)
    cerr << "Lhs contains " << columns_in_file1 << " columns\n";

  if (0 == columns_in_file1)
  {
    cerr << "Yipes, no collumns in lhs\n";
    return 0;
  }

  column_title = new IWString[columns_in_file1];
  differences_in_column = new int[columns_in_file1];

  for (int i = 0; i < columns_in_file1; i++)
  {
    differences_in_column[i] = 0;
  }

  int columns_in_file2 = header2.nwords();

  if (0 == columns_in_file2)
  {
    cerr << "Yipes, no columns in rhs header\n";
    return 0;
  }

  int * column_xref = new int[columns_in_file2]; std::unique_ptr<int[]> free_column_xref(column_xref);

  int nalloc = columns_in_file1;
  if (columns_in_file2 > columns_in_file1)
    nalloc = columns_in_file2;

  const_IWSubstring * token = new const_IWSubstring[nalloc]; std::unique_ptr<const_IWSubstring[]> free_token(token);

  return jfilecompare(input1, header1, columns_in_file1, input2, header2, columns_in_file2, column_xref, token);
}

static int
establish_identifier_cross_reference(iwstring_data_source & input2)
{
  const_IWSubstring buffer;

  if (! input2.next_record(buffer))
  {
    cerr << "Cannot read header record of second file\n";
    return 0;
  }

  off_t o = input2.tellg();

  while (input2.next_record(buffer))
  {
    buffer.truncate_at_first(' ');

    if (ignore_leading_zeros)
      buffer.remove_leading_chars('0');

    if (remove_chars_after_underscore)
      do_remove_chars_after_underscore(buffer);

    if (identifier_offset_cross_reference.end() != identifier_offset_cross_reference.find(buffer))
      cerr << "Warning, duplicate record for '" << buffer << "', ignored\n";
    else
      identifier_offset_cross_reference[buffer] = o;

    o = input2.tellg();
  }

  if (verbose)
    cerr << "Identifier/offset cross reference contains " << identifier_offset_cross_reference.size() << " items\n";

  return input2.seekg(0);
}

static int
jfilecompare(const char * fname1, const char * fname2)
{
  iwstring_data_source input1(fname1);

  if (! input1.ok())
  {
    cerr << "Cannot open input file '" << fname1 << "'\n";
    return 0;
  }

  iwstring_data_source input2(fname2);

  if (! input2.ok())
  {
    cerr << "Cannot open input file '" << fname2 << "'\n";
    return 0;
  }

  if (! same_row_order)
  {
    if (! establish_identifier_cross_reference(input2))
    {
      cerr << "Cannot establish identifier cross reference in file2 '" << fname2 << "'\n";
      return 0;
    }
  }

  return jfilecompare(input1, input2);
}

static int
final_verbose_report(std::ostream & os)
{
  os << "Read " << records_read << " records. " << records_differing << " records contained differences\n";
  if (differences_less_than_tolerance)
    os << differences_less_than_tolerance << " differences less than " << tolerance << '\n';
  if (differences_within_noise)
    os << differences_within_noise << " differences in values less than noise value " << noise << '\n';

  if (missing_value_mismatches)
    os << missing_value_mismatches << " missing value mismatches on " << records_with_missing_value_mismatches << " records ignored\n";

  for (int i = 0; i < columns_in_file1; i++)
  {
    if (differences_in_column[i])
    {
      os << differences_in_column[i] << " differences in column " << (i + 1) << " '" << column_title[i] << "'\n";
    }
  }

  for (int i = 0; i < differences_per_record.number_elements(); i++)
  {
    if (differences_per_record[i])
      os << differences_per_record[i] << " records had " << i << " different values\n";
  }

  return os.good();
}

int
jfilecompare(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vt:cM:mxzGul:V:D:n:q:P:eob");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised options encountered\n";
    usage(5);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('t'))
  {
    if (! cl.value('t', tolerance) || tolerance < 0.0)
    {
      cerr << "The comparison tolerance must be non-negative\n";
      usage(5);
    }

    if (verbose)
      cerr << "Differences less than " << tolerance << " will be ignored\n";
  }

  if (cl.option_present('b'))
  {
    interpret_tolerance_as_absolute_difference = 1;

    if (verbose)
      cerr << "Tolerance interpreted as an absolute difference\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', noise) || noise < 0.0)
    {
      cerr << "The noise level must be non-negative\n";
      usage(5);
    }

    if (verbose)
      cerr << "Values less than " << noise << " will be considered noise\n";
  }

  if (cl.option_present('q'))
  {
    int i = 0;
    const_IWSubstring q;
    while (cl.value('q', q, i++))
    {
      if ("mid" == q)
      {
        report_missing_identifiers = 0;

        if (verbose)
          cerr << "Will not report missing identifiers\n";
      }
      else if ("dup" == q)
      {
//      report_duplicate_identifiers = 0;    implement this sometime

        if (verbose)
          cerr << "Will not report duplicate identifiers\n";
      }
      else
      {
        cerr << "Unrecognised -q qualifier '" << q << "'\n";
        usage(5);
      }
    }
  }

  if (cl.option_present('c'))
  {
    ignore_case_when_comparing_descriptor_names = 1;

    if (verbose)
      cerr << "Will ignore case when comparing descriptor names\n";
  }

  if (cl.option_present('P'))
  {
    ignore_prefix = cl.string_value('P');

    if (verbose)
      cerr << "Will strip leading '" << ignore_prefix << "' from descriptor names\n";
  }

  if (cl.option_present('M'))
  {
    missing_value = cl.string_value('M');

    if (verbose)
      cerr << "Missing value '" << missing_value << "'\n";
  }

  if (cl.option_present('m'))
  {
    ignore_missing_value_mismatches = 1;

    if (verbose)
      cerr << "Will not report values that differ because of missing values\n";
  }

  if (cl.option_present('z'))
  {
    ignore_leading_zeros = 1;

    if (verbose)
      cerr << "Identifiers compared without regard to leading zero's\n";
  }

  if (cl.option_present('u'))
  {
    remove_chars_after_underscore = 1;

    if (verbose)
      cerr << "Will discard characters after first underscore in identifiers\n";
  }

  if (cl.option_present('x'))
  {
    same_row_order = 0;

    if (verbose)
      cerr << "Files not necessarily in same row order\n";
  }

  if (cl.option_present('G'))
  {
    ignore_missing_identifiers = 1;

    if (verbose)
      cerr << "Will ignore missing identifiers in the 2nd file\n";
  }

  if (cl.option_present('o'))
  {
    assume_columns_in_correct_order = 1;

    if (verbose)
      cerr << "Will assume columns in same order - ignore header\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(5);
  }

  if (2 != cl.number_elements())
  {
    cerr << "Must specify two files to compare\n";
    usage(3);
  }

  if (cl.option_present('D'))
  {
    IWString d = cl.string_value('D');

    stream_for_identifiers_with_different_data.open(d.null_terminated_chars());

    if (! stream_for_identifiers_with_different_data.good())
    {
      cerr << "Cannot open stream for identifiers with different data '" << d << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Will write identifiers with different values to '" << d << "'\n";
  }

  jfilecompare(cl[0], cl[1]);

  if (verbose)
  {
    cout.flush();
    if (cl.option_present('e'))
      final_verbose_report(cout);
    else
      final_verbose_report(cerr);
  } else if (missing_columns) {
    cerr << "Warning, " << missing_columns << " columns in first file, not in subsequent\n";
  }

  if (NULL != column_title)
    delete [] column_title;

  if (NULL != differences_in_column)
    delete [] differences_in_column;

  return records_differing;
}

int
main(int argc, char ** argv)
{
  int rc = jfilecompare(argc, argv);

  return rc;
}
