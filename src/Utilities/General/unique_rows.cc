/*
  identify the unique rows in a file
*/

#include <random>
#include <fstream>
#include <algorithm>
#include <memory>

#include "re2/re2.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "iwtokeniser.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static std::unique_ptr<RE2> column_rx;

static int strip_leading_zeros = 0;

static int ignore_case = 0;

static int show_counts = 0;

static char truncate_at = '\0';

static int do_write = 1;

static int number_instances_to_allow = 1;

static char token_separator = ' ';

static int process_all_columns = 0;

static int sort_dash_o_output = 0;

static double probability_reject_duplicate = 1.0;

static int is_descriptor_file = 0;

std::default_random_engine generator;

std::uniform_real_distribution<double> distribution(0.0,1.0);

std::mt19937_64 rng;

static int input_is_quoted_tokens = 0;

static resizable_array_p<IWString> descriptors_to_process;

static int header_records_to_skip = 0;

// Want the ability to only examine part of a token. If this is set,
// then only the characters that match the regex are added to the hash
// of already seen things.
// If this does not match, it is a fatal error.
static std::unique_ptr<RE2> text_subset;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Identifies the unique rows in a file\n";
  cerr << " -c <col>       only consider column(s) <col>\n";
  cerr << " -x <col>       consider all column(s) except <col>\n";
  cerr << " -D <file>      file for duplicate records\n";
  cerr << " -o             display counts of how many times each identifier encountered\n";
  cerr << " -O <fname>     write counts and times each identifier found to <fname>\n";
  cerr << " -s             sort output from both -o and -O. Simlified output too, just token and count\n";
//cerr << " -t <tol>       for columns that can be interpreted as numeric, allow a tolerance\n";
  cerr << " -z             remove leading 0's from comparisons\n";
  cerr << " -n             no output, just collect statistics\n";
  cerr << " -nc            ignore case during comparisons\n";
  cerr << " -r <number>    number of instances of each unique item to write (default 1)\n";
  cerr << " -p <prob>      discard duplicate records with probability <prob>\n";
  cerr << " -trunc <char>  truncate fields at first <char> - useful for truncating decimals\n";
  cerr << " -whash         wait for hash de-allocation - very slow, don't use\n";
  cerr << " -tab,-csv      deal with differently formatted files\n";
  cerr << " -i <...>       input column separator\n";
  cerr << " -nwz           No Warnings about Zero length text comparisons\n";
  cerr << " -subset        regexp specifying which part of the column to consider '(..*)...' would be all except last 3 chars\n";
  cerr << " -h <n>         first <n> records in file\n";
  cerr << " -j             treat as descriptor file\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

/*
  We function in one of two ways. If all column specifications are positive,
  we fill the columns_to_process array.
  But if there are negative column numbers present, then we have a list
  of columns, some of which will be positive, some will be negative

*/

static extending_resizable_array<int> columns_to_process;

static resizable_array<int> list_of_column_numbers;

static int negative_column_numbers_present = 0;

static extending_resizable_array<int> columns_to_ignore;

/*
  If both arrays only deal with column 1, then no reason to check column 2
*/

static int highest_column_number_to_check = 0;

static IWString_STL_Hash_Set previously_seen;

static IW_STL_Hash_Map_int times_seen;

static int lines_read = 0;

static int duplicate_lines_found = 0;

static int lines_written = 0;

static IWString_and_File_Descriptor stream_for_duplicate_lines;

static int warn_zero_length_comparisons = 1;

static int
read_previously_found_identifiers(iwstring_data_source & input)
{
  IWString buffer;

  while (input.next_record(buffer))
  {
    buffer.truncate_at_first(' ');  // hard coded to be first column

    previously_seen.insert(buffer);
  }

  return 1;
}

static int
read_previously_found_identifiers(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open previously seen identifier file '" << fname << "'\n";
    return 0;
  }

  return read_previously_found_identifiers(input);
}

// Alter `comparison_text` to be only the parts that are matched by `text_subset`.
int
TruncateToSubset(RE2* text_subset, IWString& text) {
  const re2::StringPiece tmp(text.data(), text.length());
  std::string submatch;
  if (! RE2::FullMatch(tmp, *text_subset, &submatch)) {
    cerr << "TruncateToSubset:no match to " << text_subset->pattern() << " in '" << text << "'\n";
    return 0;
  }

  text = submatch;

  return 1;
}

/*
  We have negative column numbers, fill the per column array
*/

static int
fill_my_columns_to_check(const const_IWSubstring & buffer,
                         const resizable_array<int> & columns_to_process,
                         extending_resizable_array<int> & my_columns_to_process)
{
  int ncol = buffer.nwords();

  int n = columns_to_process.number_elements();

  int rc = -1;    // highest column to check

  for (int i = 0; i < n; i++)
  {
    int c = columns_to_process[i];

    if (c >= 0)
    {
      my_columns_to_process[c] = 1;
      if (c > rc)
        rc = c;
      continue;
    }

    if (c + ncol < 0)
    {
      cerr << "Record has only " << ncol << " columns, cannot fetch column " << c << endl;
      return -1;
    }

    my_columns_to_process[c + ncol] = 1;

    if (c + ncol > rc)
      rc = c + ncol;
  }

  return rc;
}

static int
nextword_quoted_tokens(const const_IWSubstring & buffer,
                       const_IWSubstring & token,
                       int & i,
                       const char token_separator)
{
#ifdef DEBUG_NEXTWORD_QUOTED_TOKENS
  cerr << "Ftch from '" << buffer << "', i = " << i << endl;
#endif

  const char dquote = '"';

  const int n = buffer.length();

  if (i >= n)
    return 0;

#ifdef DEBUG_NEXTWORD_QUOTED_TOKENS
  cerr << "nextword_quoted_tokens begin at " << i << " char '" << buffer[i] << "'\n";
#endif

  int in_quote;
  if (dquote == buffer[i])
  {
    in_quote = 1;
    i++;
  }
  else
    in_quote = 0;

  const int istart = i;

  for ( ; i < n; ++i)
  {
//  cerr << "Checking '" << buffer[i] << "' quote " << in_quote << endl;
    if (dquote == buffer[i])
      in_quote = ! in_quote;
    else if (in_quote)
      ;
    else if (token_separator == buffer[i])
    {
      token.set(buffer.rawchars() + istart, i - istart);
      if (token.ends_with(dquote))
        token.chop();

#ifdef DEBUG_NEXTWORD_QUOTED_TOKENS
      cerr << "Found '" << token << "', i = " << i << ", just examined '" << buffer[i] << "'\n";
#endif
      i++;
      return 1;
    }
  }

  token.set(buffer.rawchars() + istart, n - istart);
#ifdef DEBUG_NEXTWORD_QUOTED_TOKENS
  cerr << "Final token '" << token << "'\n";
#endif
  if (token.ends_with(dquote))
    token.chop();

  return 1;
}

static int
get_next_token(const const_IWSubstring & buffer,
               const_IWSubstring & token,
               int & i,
               char token_separator)
{
  if (' ' == token_separator)
    return buffer.nextword(token, i);

  if (input_is_quoted_tokens)
    return nextword_quoted_tokens(buffer, token, i, token_separator);

  return buffer.nextword_single_delimiter(token, i, token_separator);
}

static int
unique_rows(const const_IWSubstring & buffer,
            const int record_number,
            IWString & output_buffer)
{
  IWString comparison_text;

  extending_resizable_array<int> my_columns_to_process;

  if (negative_column_numbers_present)
  {
    highest_column_number_to_check = fill_my_columns_to_check(buffer, list_of_column_numbers, my_columns_to_process);

    if (highest_column_number_to_check < 0)
      return 0;
  }

  comparison_text.resize(buffer.length());

//#define DEBUG_UNIQUE_ROWS
#ifdef DEBUG_UNIQUE_ROWS
  cerr << "Looking at " << columns_to_process.number_elements() << " columns to process\n";
  for (int i = 0; i < columns_to_process.number_elements(); i++)
  {
    cerr << " process column " << i << " " << columns_to_process[i] << endl;
  }
  for (int i = 0; i < columns_to_ignore.number_elements(); i++)
  {
    cerr << " ignore column " << i << " " << columns_to_ignore[i] << endl;
  }
#endif

  IWTokeniser iwt(buffer);
  iwt.set_sep(token_separator);
  if (' ' != token_separator)
    iwt.set_empty_fields_valid(1);
  if (input_is_quoted_tokens)
    iwt.set_quoted_tokens(1);

  const_IWSubstring token;

//cerr << "Line " << __LINE__ << endl;

  if (process_all_columns)
    comparison_text = buffer;
  else
  {
    for (int col = 0; iwt.next_token(token); col++)
    {
//    cerr << " col " << col << " token '" << token << "'\n";
      int append_token = 0;    // if this column will be part of the matching
  
      if (column_rx)
      {
        re2::StringPiece tmp(token.data(), token.length());
        if (RE2::PartialMatch(tmp, *column_rx))
          append_token = 1;
      }
      else if (negative_column_numbers_present)
      {
        if (my_columns_to_process[col])
          append_token = 1;
      }
      else if (columns_to_process.number_elements() > 0)
      {
        if (columns_to_process[col])
          append_token = 1;
      }
      else if (0 == columns_to_ignore[col])
        append_token = 1;
  
      if (! append_token)
        continue;
  
      if (strip_leading_zeros)
        token.remove_leading_chars('0');
  
      comparison_text.append_with_spacer(token);

      if (highest_column_number_to_check > 0 && col >= highest_column_number_to_check)
        break;
    }
  }

  if (0 == comparison_text.length())
  {
    if (warn_zero_length_comparisons)
      cerr << "Warning, zero length comparison text\n";
    return 1;
  }

  if (ignore_case)
    comparison_text.to_lowercase();

  if ('\0' == truncate_at)
    ;
  else if (comparison_text.starts_with(truncate_at)) 
    cerr << "Cannot truncate '" << comparison_text << "' at first '" << truncate_at << "'\n";
  else
    comparison_text.truncate_at_first(truncate_at);

  if (text_subset) {
    if (! TruncateToSubset(text_subset.get(), comparison_text)) {
      return 0;
    }
  }

  if (! show_counts)
    ;
  else if (1 == record_number && is_descriptor_file)
    ;
  else
    times_seen[comparison_text] += 1;

#ifdef DEBUG_UNIQUE_ROWS
  cerr << "From '" << buffer << "' comparison_text '" << comparison_text << "'\n";
  if (previously_seen.contains(comparison_text))
    cerr << "Seen " << times_seen[comparison_text] << endl;
#endif

  if (! previously_seen.contains(comparison_text))   // definitely not a dup
  {
    if (1 == record_number && is_descriptor_file)
      ;
    else
    {
      previously_seen.insert(comparison_text);
      if (number_instances_to_allow > 1)    // we need to count
        times_seen[comparison_text] = 1;
    }
  }
  else if (number_instances_to_allow > 1 && (times_seen[comparison_text]+show_counts) < number_instances_to_allow)  // not considered a duplicate, strangeness with show_counts because incremented below
  {
    if (! show_counts)     // was incremented above if show_counts
      times_seen[comparison_text]++;
  }
  else if (probability_reject_duplicate < 1.0 && distribution(generator) < probability_reject_duplicate)
  {
    if (! show_counts)
      times_seen[comparison_text]++;
  }
  else       // treated as a duplicate
  {
    duplicate_lines_found++;

    if (stream_for_duplicate_lines.is_open())
    {
      stream_for_duplicate_lines << buffer << '\n';
      stream_for_duplicate_lines.write_if_buffer_holds_more_than(32768);

      return stream_for_duplicate_lines.good();
    }

    if (verbose > 2)
      cerr << "Duplicate: '" << buffer << "'\n";

    return 1;
  }

  if (do_write)
  {
    output_buffer << buffer << '\n';

    lines_written++;
  }

//cerr << "After insertion hash contains " << previously_seen.size() << endl;

return 1;
  if (show_counts)
    times_seen[comparison_text] = 1;
  else if (number_instances_to_allow > 1)
  {
    if (! times_seen.contains(comparison_text))
      times_seen[comparison_text] = 1;
    else
      times_seen[comparison_text] += 1;
  }
  
  return 1;
}

static int
determine_columns_to_process (const const_IWSubstring & buffer,
                              const resizable_array_p<IWString> & descriptors_to_process,
                              extending_resizable_array<int> & columns_to_process)
{
  const int nd = descriptors_to_process.number_elements();

  int * dfound = new_int(nd); std::unique_ptr<int[]> free_dfound(dfound);

  const_IWSubstring token;
  int i = 0;
  for (int col = 0; get_next_token(buffer, token, i, token_separator); ++col)
  {
    for (int j = 0; j < nd; ++j)
    {
      if (token == *descriptors_to_process[j])
      {
        columns_to_process[col] = 1;
        dfound[j] = 1;

        if (col > highest_column_number_to_check)
          highest_column_number_to_check = col;

        if (verbose)
          cerr << "Descriptor '" << *descriptors_to_process[j] << "' found in column " << (col+1) << endl;
        break;
      }
    }
  }

  int rc = 1;
  for (int i = 0; i < nd; ++i)
  {
    if (! dfound[i])
    {
      cerr << "No match for descriptor '" << *descriptors_to_process[i] << "'\n";
      rc = 0;
    }
  }

  return rc;
}

static int
unique_rows (iwstring_data_source & input,
             IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (header_records_to_skip > 0)
  {
    for (int i = 0; i < header_records_to_skip; ++i)
    {
      if (! input.next_record(buffer))
      {
        cerr << "Cannot read header record " << i << endl;
        return 0;
      }
    }
  }

  while (input.next_record(buffer))
  {
    if (0 == lines_read && descriptors_to_process.number_elements() > 0)
    {
      if (! determine_columns_to_process(buffer, descriptors_to_process, columns_to_process))
      {
        cerr << "Cannot discern columns to process '" << buffer << "'\n";
        return 0;
      }
    }

    lines_read++;

    if (! unique_rows(buffer, lines_read, output))
      return 0;

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int
unique_rows (const char * fname,
             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return unique_rows(input, output);
}

static int
extract_column_range (const const_IWSubstring & cs,
                      resizable_array<int> & col)
{
  const_IWSubstring s1, s2;
  if (! cs.split(s1, '-', s2) || 0 == s1.length() || 0 == s2.length())
    return 0;

  int c1;
  if (! s1.numeric_value(c1) || c1 < 1)
  {
    cerr << "Invalid column range start '" << s1 << "'\n";
    return 0;
  }

  int c2;
  if (! s2.numeric_value(c2) || c2 < c1)
  {
    cerr << "Invalid column range end '" << s1 << "'\n";
    return 0;
  }

  for (int i = c1; i <= c2; i++)
  {
    col.add(i-1);
  }

  return 1;
}

static int
extract_comma_separated_columns (const const_IWSubstring & cs,
                                 resizable_array<int> & col)
{
  int i = 0;
  const_IWSubstring token;
  while (cs.nextword(token, i, ','))
  {
    int c;

    if (! token.numeric_value(c) || c < 1)
    {
      cerr << "Invalid column number '" << token << "'\n";
      return 0;
    }

    col.add(c-1);
  }

  return 1;
}

static int
fetch_columns (Command_Line_v2 & cl,
               char flag,
               resizable_array<int> & col)
{
  int i = 0;
  const_IWSubstring cs;
  while (cl.value(flag, cs, i++))
  {
    if ("all" == cs)
    {
      process_all_columns = 1;

      return 1;
    }

//  cerr << "Examining '" << cs << "'\n";

    if (cs.index('-') > 0)    // don't grab negative column numbers here
    {
      if (! extract_column_range(cs, col))
      {
        cerr << "Invalid range of columns '" << cs << "'\n";
        return 0;
      }
      continue;
    }

    if (cs.contains(','))
    {
      if (! extract_comma_separated_columns(cs, col))
      {
        cerr << "Invalid set of columns '" << cs << "'\n";
        return 0;
      }

      continue;
    }

    int c;
    if (! cs.numeric_value(c))
    {
       cerr << "Invalid column specification'" << cs << "'\n";
       return 0;
    }

    if (c > 0)
    {
      c--;

      col.add(c);
    }
    else
    {
      col.add(c);
      negative_column_numbers_present = 1;
    }
  }

  return col.number_elements();
}

static int
write_count_data_sorted (const IW_STL_Hash_Map_int & times_seen,
                         std::ostream & os)
{
  const auto n = times_seen.size();

  std::pair<IWString, int> * string_and_count = new std::pair<IWString, int>[n]; std::unique_ptr<std::pair<IWString, int>[]> free_string_and_count(string_and_count);

  Accumulator_Int<int> acc;

  int ndx = 0;
  for (IW_STL_Hash_Map_int::const_iterator i = times_seen.begin(); i != times_seen.end(); ++i)
  {
    string_and_count[ndx].first = (*i).first;
    const auto c = (*i).second;
    string_and_count[ndx].second = (*i).second;

    ndx++;

    if (verbose)
      acc.extra(c);
  }

  std::sort(string_and_count, string_and_count + n, [](const std::pair<IWString, int> & sc1, const std::pair<IWString, int> & sc2) {if (sc1.second > sc2.second) return 1;
                                                                                                                                      return 0;});

  for (unsigned int i = 0; i < n; ++i)
  {
    os << string_and_count[i].first << ' ' << string_and_count[i].second << '\n';
  }

  if (verbose)
    cerr << "Items seen between " << acc.minval() << " and " << acc.maxval() << " times, ave " << static_cast<float>(acc.average_if_available_minval_if_not()) << '\n';

  return 1;
}

static int
write_count_data(const IW_STL_Hash_Map_int & times_seen,
                 std::ostream & os)
{
  if (sort_dash_o_output)
    return write_count_data_sorted(times_seen, os);

  Accumulator_Int<int> acc;

  for (IW_STL_Hash_Map_int::const_iterator i = times_seen.begin(); i != times_seen.end(); ++i)
  {
      const int c = (*i).second;

      if (c >= show_counts)
        os << (*i).first << " seen " << c << " times\n";

      if (verbose)
        acc.extra(c);
    }

    if (verbose)
      cerr << "Items seen between " << acc.minval() << " and " << acc.maxval() << " times, ave " << static_cast<float>(acc.average_if_available_minval_if_not()) << '\n';

  return 1;
}

static int
unique_rows (int argc, char ** argv)
{
//Command_Line cl (argc, argv, "vc:x:D:R:z");
  Command_Line_v2 cl(argc, argv, "-v-c=s-x=ipos-D=s-R=s-z-nc-o-trunc=s-n-r=ipos-whash-P=sfile-O=s-tab-csv-vbar-comma-s-all-p=f-q-nwz-d=s-h=ipos-j=ipos-i=s-subset=s");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('z'))
  {
    strip_leading_zeros = 1;

    if (verbose)
      cerr << "Leading zero's removed during comparisons\n";
  }

  if (cl.option_present('n'))
  {
    do_write = 0;

    if (verbose)
      cerr << "Normal output suppressed\n";
  }

  int specifications = (0 != cl.option_present('c')) + (0 != cl.option_present('x')) + (0 != cl.option_present('R')) + (cl.option_count('d') > 0);

  if (0 == specifications)
    columns_to_process[0] = 1;
  else if (1 != specifications)
  {
    cerr << "Must have just one of -c, -d, -x or -R options\n";
    usage(5);
  }

  if (cl.option_present("trunc"))
  {
    const const_IWSubstring t = cl.string_value("trunc");

    if (t.length() > 1)
    {
      cerr << "Sorry, the truncate value can be just a single character '" << t << "' not valid\n";
      usage(4);
    }

    if (' ' == t[0])
    {
      cerr << "Sorry, cannot truncate at space\n";
      return 4;
    }

    truncate_at = t[0];

    if (verbose)
      cerr << "Will truncate at first '" << truncate_at << "' in each identifier\n";
  }

  if (cl.option_present("nc"))
  {
    ignore_case = 1;

    if (verbose)
      cerr << "Will ignore case when comparing records\n";
  }

  if (cl.option_present('j'))
  {
    is_descriptor_file = 1;
    if (verbose)
      cerr << "Will treat as a descriptor file\n";
  }

  if (cl.option_present("h"))
  {
    if (! cl.value("h", header_records_to_skip) || header_records_to_skip < 0)
    {
      cerr << "The number of header records to skip (-h) option must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will skip the first " << header_records_to_skip << " header records\n";
  }

  if (cl.option_present('c'))
  {
    if (! fetch_columns(cl, 'c', list_of_column_numbers))
      return 4;

    if (negative_column_numbers_present)
    {
    }
    else if (process_all_columns)
      ;
    else
    {
      for (int i = 0; i < list_of_column_numbers.number_elements(); i++)
      {
        int c = list_of_column_numbers[i];

        columns_to_process[c] = 1;
        if (c > highest_column_number_to_check)
          highest_column_number_to_check = c;
      }
    }

    if (verbose > 2)
    {
      for (int i = 0; i < columns_to_process.number_elements(); i++)
      {
        if (columns_to_process[i])
          cerr << "Will examine column " << i << endl;
      }
    }
  }
  else if (cl.option_present('x'))
  {
    resizable_array<int> tmp;
    if (! fetch_columns(cl, 'x', tmp))
      return 4;

    for (int i = 0; i < tmp.number_elements(); i++)
    {
      int c = tmp[i];
      columns_to_ignore[c] = 1;
    }
  }
  else if (cl.option_present('R'))
  {
    const_IWSubstring r = cl.string_value('R');

    re2::StringPiece tmp(r.data(), r.length());
    column_rx.reset(new RE2(tmp));
    if (! column_rx->ok())
    {
      cerr << "Invalid column regexp '" << r << "'\n";
      return 6;
    }

    if (verbose)
      cerr << "Will check columns that match '" << column_rx->pattern() << "'\n";
  }
  else if (cl.option_present('d'))
  {
    IWString s;
    for (int i = 0; cl.value('d', s, i); ++i)
    {
      s.unhtml();
      descriptors_to_process.add(new IWString(s));
    }
    is_descriptor_file = 1;
  }
  else if (0 == specifications)
    ;
  else
  {
    cerr << "Must specify either the -c or -x options\n";
    usage(6);
  }

  if (cl.option_present("csv") || cl.option_present("comma"))
  {
    token_separator = ',';
    if (verbose)
      cerr << "Will treat as comma separated\n";
  }
  else if (cl.option_present("tab"))
  {
    token_separator = '\t';
    if (verbose)
      cerr << "Will treat as tab separated\n";
  }
  else if (cl.option_present("vbar"))
  {
    token_separator = '|';
    if (verbose)
      cerr << "Input fields separated by | characters\n";
  }
  else if (cl.option_present("i"))
  {
    IWString t = cl.string_value("i");
    if (! char_name_to_char(t))
    {
      cerr << "Unrecognised input token separator '" << t << "'\n";
      return 1;
    }
    token_separator = t[0];
  }

  if (cl.option_present('o'))
  {
    show_counts = cl.option_count('o');
    if (verbose)
      cerr << "Will display counts of how many times each identifier seen\n";
  }

  if (cl.option_present('s'))
  {
    sort_dash_o_output = 1;

    if (verbose)
      cerr << "Summary statistics of frequncies will be sorted\n";

    if (0 == show_counts)
      show_counts = 1;
  }

  IWString count_fname;

  if (cl.option_present('O'))
  {
    cl.value('O', count_fname);

    if (verbose)
      cerr << "Count data written to '" << count_fname << "'\n";

    show_counts = 1;
  }

  if (! do_write && ! show_counts && 0 == verbose && ! cl.option_present('D'))
    cerr << "Will there be any output???\n";

  if (cl.option_present('r'))
  {
    cl.value('r', number_instances_to_allow);

    if (verbose)
      cerr << "Will allow at most " << number_instances_to_allow << " instances of each unique item\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', probability_reject_duplicate) || probability_reject_duplicate <= 0.0 || probability_reject_duplicate > 1.0)
    {
      cerr << "The probability of rejecting a duplicate record (-p) must be a valid fraction\n";
      usage(1);
    }

    if (verbose)
      cerr << "Duplicate records rejected with probability " << probability_reject_duplicate << endl;
  }

  if (cl.option_present('P'))
  {
    const char * p = cl.option_value('P');

    if (! read_previously_found_identifiers(p))
    {
      cerr << "Cannot read previously encountered identifiers file '" << p << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Read " << previously_seen.size() << " previously seen identifiers from '" << p << "'\n";
  }

  if (cl.option_present('q'))
  {
    input_is_quoted_tokens = 1;
    if (verbose)
      cerr << "INput is quoted tokens\n";
  }

  if (cl.option_present("nwz"))
  {
    warn_zero_length_comparisons = 0;

    if (verbose)
      cerr << "Will not report occurrences of zero length comparisons\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('D'))
  {
    IWString d (cl.string_value('D'));

    if (! stream_for_duplicate_lines.open(d.null_terminated_chars()))
    {
      cerr << "Cannot open stream for duplicates '" << d << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Duplicates written to '" << d << "'\n";
  }

  if (cl.option_present("subset")) {
    const_IWSubstring r = cl.string_value("subset");

    re2::StringPiece tmp(r.data(), r.length());
    text_subset.reset(new RE2(tmp));
    if (! text_subset->ok())
    {
      cerr << "Invalid text subset regexp '" << r << "'\n";
      return 6;
    }

    if (verbose)
      cerr << "Only the parts of a column that match " << text_subset->pattern() << " will be considered for uniqueness\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (verbose)
    {
      cerr << "Begin processing '" << cl[i] << "'";
      if (previously_seen.size() > 0)
        cerr << ", found " << previously_seen.size() << " unique items";
      cerr << endl;
    }

    if (! unique_rows(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
     cerr << "Read " << lines_read << ", " << duplicate_lines_found << " duplicates, wrote " << lines_written << endl;
     cerr << previously_seen.size() << " unique items\n";
  }

  if (! show_counts)
    ;
  else if (0 == times_seen.size())
    cerr << "No count data\n";
  else if (count_fname.length())
  {
    std::ofstream ofile(count_fname.null_terminated_chars(), std::ios::out);

    if (! ofile.good())
    {
      cerr << "Cannot open count data file '" << count_fname << "'\n";
      return 2;
    }

    write_count_data(times_seen, ofile);
  }
  else
    write_count_data(times_seen, cerr);

#ifdef SHOW_SET
  for (IWString_STL_Hash_Set::const_iterator f = previously_seen.begin(); f != previously_seen.end(); ++f)
  {
    cerr << "Saw '" << (*f) << "'\n";
  }
#endif

  if (stream_for_duplicate_lines.is_open())
    stream_for_duplicate_lines.close();

  if (! cl.option_present("whash"))
    _exit(0);

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = unique_rows(argc, argv);

  return rc;
}
