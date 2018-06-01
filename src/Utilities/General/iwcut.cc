/*
  My own cut utility
*/

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                                          + __GNUC_PATCHLEVEL__)

#include <iostream>
#include <memory>
#include <algorithm>
#include <random>
#include <cctype>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "misc.h"
#include "iwstring_data_source.h"

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Extracts columns - like cut\n";
  cerr << " -f <col>         extract single column\n";
  cerr << " -f <col1,col2>   extract multiple columns\n";
  cerr << " -f <col1-col2>   extract range of columns\n";
  cerr << " -F <file>        column numbers are in <file>\n";
  cerr << " -m               descriptors in -F file are on multiple lines *\n";
  cerr << " -d <descriptor>  extract single descriptor\n";
  cerr << " -d <d1,d2,d3>    extract multiple descriptors\n";
  cerr << " -D <file>        descriptors are in <file>\n";
  cerr << " -b               do NOT always add the first column when specifying column headings\n";
  cerr << " -m               descriptors in -D file are on multiple lines *\n";
  cerr << " -K <char>        KLUDGE, remove descriptor name prefixes up to <char>\n";
  cerr << " -R <regexp>      print column(s) matching <regexp> - do match on each line\n";
  cerr << " -x               prints the columns that do not match\n";
  cerr << " -r               match descriptor names as regular expressions\n";
  cerr << " -c               ignore case when comparing descriptor names\n";
  cerr << " -i <char>        input  token separator (default ' ')\n";
  cerr << " -t               input is tab separated\n";
  cerr << " -q <k|rm>        input may contain quoted tokens, keep or remove the quotes\n";
  cerr << " -o <char>        output token separator (default ' ')\n";
  cerr << " -z <string>      pad missing descriptors/columns with <string>\n";
  cerr << " -u               consecutive delimiters mean empty words in between\n";
  cerr << " -E <col>         columns that must NOT be empty - record not written if empty. Use . for any\n";
#if GCC_VERSION > 50000
  cerr << " -a <n>           select <n> columns at random\n";
#endif
  cerr << " -g <char>        in tokens with whitespace, replace with <char>\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int verbose = 0;

static int descriptor_names_on_one_line = 1;

static int suppress_selected_columns = 0;

static int ignore_case = 0;

static int match_descriptor_names_as_regular_expressions = 0;

static IWString missing_column_fill;

static char input_token_separator = ' ';

static char output_token_separator = ' ';

static char descriptor_name_remove_to = ' ';

static int random_columns = 0;

/*
  If we are filling missing columns, we may have a constant thing to add to every output record

  Note that what gets appended to the header is different from what gets appended to each record
*/

static IWString extra_stuff_header_record;
static IWString extra_stuff_each_record;

/*
  The regular expression associated with the -R option
*/

static IW_Regular_Expression token_regexp;

static int consecutive_delimiters_mean_empty_words_between = 0;

static int input_is_quoted_tokens = 0;

static extending_resizable_array<int> columns_that_must_not_be_empty;
static int all_fields_must_have_content = 0;

static int add_first_token_when_processing_descriptors = 1;

static int records_discarded_for_zero_field = 0;

static char gsub_spaces_in_tokens = ' ';

static int
invert_selections(resizable_array<int> & columns_requested, int ncol,
                  int is_descriptor_file, int * tmp)
{
  int nr = columns_requested.number_elements();

  for (int i = 0; i < nr; i++)
  {
    int c = columns_requested[i];

    tmp[c] = 0;
  }

  if (is_descriptor_file)
    tmp[0] = 1;

  columns_requested.resize(0);
  columns_requested.resize(ncol);

  for (int i = 0; i < ncol; i++)
  {
    if (tmp[i])
      columns_requested.add(i);
  }

  return 1;
}

//#define DEBUG_IWCUT

static int
invert_selections(resizable_array<int> & columns_requested,
                  int ncol,
                  int is_descriptor_file)
{
#ifdef DEBUG_IWCUT
  cerr << "Inverting selections, input contains " << ncol << " columns\n";
#endif

  int * tmp = new_int(ncol, 1); std::unique_ptr<int> free_tmp(tmp);

  return invert_selections(columns_requested, ncol, is_descriptor_file, tmp);
}

static int
parse_dash_D_option_multiple_lines(iwstring_data_source & input,
                     resizable_array_p<IWString> & descriptors_to_get)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.contains(input_token_separator))
      buffer.truncate_at_first(input_token_separator);

    if (0 == buffer.length())
      continue;

    IWString * tmp = new IWString(buffer);

    descriptors_to_get.add(tmp);
  }

  return 1;
}

static int
parse_dash_D_option_one_line(const const_IWSubstring & buffer,
                     resizable_array_p<IWString> & descriptors_to_get)
{
  int i = 0;
  const_IWSubstring dname;

  while (buffer.nextword(dname, i))
  {
    IWString * tmp = new IWString(dname);

    descriptors_to_get.add(tmp);
  }

  return 1;
}

static int
parse_dash_D_option_one_line(iwstring_data_source & input,
                     resizable_array_p<IWString> & descriptors_to_get)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read descriptor names from file\n";
    return 0;
  }

  return parse_dash_D_option_one_line(buffer, descriptors_to_get);
}

static int
parse_dash_D_option(iwstring_data_source & input,
                    resizable_array_p<IWString> & descriptors_to_get)
{
  if (descriptor_names_on_one_line)
    return parse_dash_D_option_one_line(input, descriptors_to_get);
  else
    return parse_dash_D_option_multiple_lines(input, descriptors_to_get);
}

static int
parse_dash_D_option(const_IWSubstring buffer,
                    resizable_array_p<IWString> & descriptors_to_get)
{
  iwstring_data_source input(buffer);
  if (! input.ok())
  {
    cerr << "Cannot open '" << buffer << "'\n";
    return 0;
  }

  return parse_dash_D_option(input, descriptors_to_get);
}

static int
parse_dash_D_option(Command_Line & cl,
                    char flag,
                    resizable_array_p<IWString> & descriptors_to_get)
{
  int i = 0;
  const_IWSubstring d;
  while (cl.value(flag, d, i++))
  {
    if (! parse_dash_D_option(d, descriptors_to_get))
    {
      cerr << "INvalid -" << flag << " combination '" << d << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
parse_dash_d_option(const_IWSubstring buffer,      // our own copy
                     resizable_array_p<IWString> & descriptors_to_get)
{

  int i = 0;
  const_IWSubstring token;

  if (match_descriptor_names_as_regular_expressions)    // do NOT tokenise on comma, that can be part of the RX
  {
    if (buffer.contains(','))
      cerr << "Assuming regular expression containing comma\n";

    IWString * tmp = new IWString(buffer);
    descriptors_to_get.add(tmp);

    return 1;
  }

  while (buffer.nextword(token, i, ','))
  {
    IWString * tmp = new IWString(token);

    descriptors_to_get.add(tmp);
  }

  return 1;
}

static int
parse_dash_d_option(Command_Line & cl,
                     char flag,
                     resizable_array_p<IWString> & descriptors_to_get)
{
  int i = 0;
  const_IWSubstring d;
  while (cl.value(flag, d, i++))
  {
    if (! parse_dash_d_option(d, descriptors_to_get))
    {
      cerr << "Invalid -" << flag << " option '" << d << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  Consume digits from the front of BUFFER
*/

static int
get_number(const_IWSubstring & buffer,
            int & zresult)
{
  zresult = 0;

  int chars_consumed = 0;

  while (buffer.length())
  {
    int j = buffer[0] - '0';

    if (j > 9 || j < 0)
      return chars_consumed;

    zresult = zresult * 10 + j;
    chars_consumed++;
    buffer.remove_leading_chars(1);
  }

  return chars_consumed;
}

static int
get_range(const_IWSubstring & buffer,
           int rstart,
           resizable_array<int> & columns_requested)
{
  int rend;

  if (! get_number(buffer, rend))
  {
    cerr << "Invalid range specifier '" << buffer << "'\n";
    return 0;
  }

  if (rend < rstart)
  {
    cerr << "Invalid range specifier - must increase\n";
    return 0;
  }

  for (int i = rstart + 1; i <= rend; i++)
  {
    columns_requested.add_if_not_already_present(i - 1);
  }

  return 1;
}

static int
add_column_from_text_representation(const_IWSubstring buffer,
                     resizable_array<int> & columns_requested)
{
  int col;
  if (! buffer.numeric_value(col) || col < 1)
  {
    cerr << "INvalid column specifier '" << buffer << "'\n";
    return 0;
  }

  col--;

  columns_requested.add(col);

  return 1;
}

static int
parse_dash_F_option_one_line(const_IWSubstring & buffer,
                     resizable_array<int> & columns_requested)
{
  int i = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i))
  {
    if (! add_column_from_text_representation(token, columns_requested))
      return 0;
  }

  return 1;
}

static int
parse_dash_F_option_one_line(iwstring_data_source & input,
                     resizable_array<int> & columns_requested)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read first record of -F file\n";
    return 0;
  }

  return parse_dash_F_option_one_line(buffer, columns_requested);
}

static int
parse_dash_F_option_multiple_lines(iwstring_data_source & input,
                     resizable_array<int> & columns_requested)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! add_column_from_text_representation(buffer, columns_requested))
      return 0;
  }

  return 1;
}

static int
parse_dash_F_option(const const_IWSubstring & fname,
                     resizable_array<int> & columns_requested)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open -F file '" << fname << "'\n";
    return 0;
  }

  if (descriptor_names_on_one_line)
    return parse_dash_F_option_one_line(input, columns_requested);
  else
    return parse_dash_F_option_multiple_lines(input, columns_requested);
}

static int
parse_dash_F_option(Command_Line & cl,
                     char flag,
                     resizable_array<int> & columns_requested)
{
  int i = 0;
  const_IWSubstring f;

  while (cl.value(flag, f, i++))
  {
    if (! parse_dash_F_option(f, columns_requested))
    {
      cerr << "Invalid -" << flag << " qualifier '" << f << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  -f qualifiers look like

   -f 3
   -f 1,2,3
   -f 4-5
*/

static int
parse_dash_f_option(const_IWSubstring buffer,     // our own copy, not passed by reference
                     resizable_array<int> & columns_requested)
{
  const_IWSubstring f;
  int i = 0;

  while (buffer.nextword(f, i, ','))
  {
    int col;

    if (f.starts_with('-'))
    {
      if (! f.numeric_value(col) || 0 == col)
      {
        cerr << "Invalid negative token offset '" << f << "'\n";
        return 0;
      }
      columns_requested.add_if_not_already_present(col);
      continue;
    }

    if (! get_number(f, col) || 0 == col)
    {
      cerr << "Invalid column number '" << f << "'\n";
      return 0;
    }

    if (col > 0)
      columns_requested.add_if_not_already_present(col - 1);
    else
      columns_requested.add_if_not_already_present(col);

    if (0 == f.length())     // just a single column number
      continue;

    if (f.starts_with('-'))
    {
      f.remove_leading_chars(1);
      if (! get_range(f, col, columns_requested))
      {
        cerr << "Invalid range specifier -'" << f << "'\n";
        return 0;
      }
    }
    else
    {
      cerr << "Invalid -f qualifier '" << f << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
parse_dash_f_option (Command_Line & cl,
                     char flag,
                     resizable_array<int> & columns_requested)
{
  int i = 0;
  const_IWSubstring f;

  while (cl.value(flag, f, i++))
  {
    if (! parse_dash_f_option(f, columns_requested))
    {
      cerr << "Invalid -" << flag << " qualifier '" << f << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
append_possible_translated(const const_IWSubstring & token,
                           const char rep,
                           IWString_and_File_Descriptor & output)
{
  const int n = token.length();

  output.make_room_for_extra_items(n);

  int rc = 0;

  for (int i = 0; i < n; ++i)
  {
    const char c = token[i];
    if (isspace(c))
      output << rep;
    else
    {
      output << c;
      rc++;
    }
  }
  return rc;
}

static int
iwcut_token_regexp(const const_IWSubstring & buffer,
                   const resizable_array<int> & columns_requested,
                   IWString_and_File_Descriptor & output)
{
  int i = 0;
  const_IWSubstring token;

  int col = 0;

  int need_separator = false;

  while (buffer.nextword(token, i, input_token_separator))
  {
    int m;

    if (columns_requested.contains(col))
      m = 1;
    else if (token_regexp.matches(token))
      m = 1;
    else
      m = 0;

    col++;

    if (suppress_selected_columns)
      m = !m;

    if (! m)
      continue;

    if (need_separator)
      output << output_token_separator;

    need_separator = true;
    append_possible_translated(token, gsub_spaces_in_tokens, output);
  }

  output += '\n';

  return 1;
}

static int
next_quoted_token_inner(const const_IWSubstring & buffer,
                  const int s,
                  const char input_token_separator,
                  const_IWSubstring & token)
{
  if (input_token_separator == buffer[s])   // empty field
  {
    token.make_empty();
    return 1;
  }

  const int n = buffer.length();

  const char dquote = '"';

  int in_quote = buffer[s] == dquote;

  for (int i = s+1; i < n; ++i)
  {
    if (dquote == buffer[i])
      in_quote = ! in_quote;
    else if (in_quote)
      ;
    else if (input_token_separator == buffer[i])
    {
      token.set(buffer.rawchars() + s, i - s);
      return 1;
    }
  }

  token.set(buffer.rawchars() + s, n - s);
  return 1;
}

static int
next_quoted_token(const const_IWSubstring & buffer,
                  const int s,
                  const char input_token_separator,
                  const_IWSubstring & token)
{
  const int n = buffer.length();

  if (s == n)    // last token is empty
  {
    token.make_empty();
    return 1;
  }

  assert(s < n);

  next_quoted_token_inner(buffer, s, input_token_separator, token);

  if (0 == token.length())
    return 1;

  if (2 != input_is_quoted_tokens)
    return 1;

  const char dquote('"');

  if (token.starts_with(dquote) && token.ends_with(dquote))
  {
    token.remove_leading_chars(1);
    token.chop();
  }

  return 1;
}

static int
iwcut_token_regexp(iwstring_data_source & input,
                    const resizable_array<int> & columns_requested,
                    IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! iwcut_token_regexp(buffer, columns_requested, output))
      return 0;

    if (output.size() > 32768)
      output.write_whole_blocks_shift_unwritten();
  }

  return 1;
}

static int
iwcut(const const_IWSubstring & buffer,
      const resizable_array<int> & word_beginnings,
      const resizable_array<int> & columns_requested,
      IWString_and_File_Descriptor & output)
{
  const int initial_size_output = output.number_elements();

  int nr = columns_requested.number_elements();

#ifdef DEBUG_IWCUT
  for (int i = 0; i < word_beginnings.number_elements(); i++)
  {
    cerr << ' ' << word_beginnings[i];
  }
  cerr << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    int c = columns_requested[i];

    if (c < 0)
      c = word_beginnings.number_elements() + c;

    int j;

    if (word_beginnings.ok_index(c))
      j = word_beginnings[c];
    else if (missing_column_fill.length() > 0)
    {
      if (i > 0)
        output += output_token_separator;
      output += missing_column_fill;
      continue;
    }
    else
    {
      cerr << "Column " << c << " out of range, only " << word_beginnings.number_elements() << " columns in input (use the -z option to fix)\n";
      return 0;
    }

    const_IWSubstring token;

    if (input_is_quoted_tokens)
      next_quoted_token(buffer, j, input_token_separator, token);
    else if (consecutive_delimiters_mean_empty_words_between)
      (void) buffer.nextword_single_delimiter(token, j, input_token_separator);
    else
      (void) buffer.nextword(token, j, input_token_separator);

    if (0 == token.length() && (all_fields_must_have_content || columns_that_must_not_be_empty[c]))
    {
      output.resize_keep_storage(initial_size_output);
      records_discarded_for_zero_field++;
      return 1;
    }

    if (i > 0)
      output += output_token_separator;

    if (' ' != gsub_spaces_in_tokens)
      append_possible_translated(token, gsub_spaces_in_tokens, output);
    else
      output += token;
  }

  if (extra_stuff_each_record.length())
    output << extra_stuff_each_record;

  output += '\n';

  return 1;
}

//#define DEBUG_QUOTED_WB

static int
locate_quoted_tokens_word_beginnings(const const_IWSubstring & buffer,
                                     resizable_array<int> & word_beginnings,
                                     const char input_token_separator)
{
  const int n = buffer.length();

  if (0 == n)
    return 0;

  const char dquote = '"';

  word_beginnings.add(0);

  int in_quote;
  if (input_token_separator == buffer[0])
  {
    word_beginnings.add(1);
    in_quote = false;
  }
  else
    in_quote = buffer[0] == dquote;

  for (int i = 1; i < n; ++i)
  {
//  cerr << " char " << i << " '" << buffer[i] << "' quote " << in_quote << endl;
    if (dquote == buffer[i])
      in_quote = ! in_quote;
    else if (in_quote)
      ;
    else if (input_token_separator == buffer[i])
      word_beginnings.add(i+1);
  }

#ifdef DEBUG_QUOTED_WB
  for (int i = 0; i < word_beginnings.size(); ++i)
  {
    cerr << " wb " << i << ' ' << word_beginnings[i] << endl;
  }
#endif

  return word_beginnings.number_elements();
}

static int
iwcut(const const_IWSubstring & buffer,
      const resizable_array<int> & columns_requested,
      IWString_and_File_Descriptor & output)
{
  static int columns_in_input = 0;

  resizable_array<int> word_beginnings;

  if (columns_in_input > 0)
    word_beginnings.resize(columns_in_input);

//cerr << "Line " << __LINE__ << " iqt " << input_is_quoted_tokens << endl;

  int ncol;
  if (input_is_quoted_tokens)
    ncol = locate_quoted_tokens_word_beginnings(buffer, word_beginnings, input_token_separator);
  else if (consecutive_delimiters_mean_empty_words_between)
    ncol = buffer.locate_word_beginnings_single_delimiter(word_beginnings, input_token_separator);
  else
    ncol = buffer.locate_word_beginnings(word_beginnings, input_token_separator);

#ifdef DEBUG_IWCUT
  cerr << "Processing '" << buffer << "'\n";
  cerr << "ncol " << ncol << " count " << buffer.ccount(input_token_separator) << endl;
#endif

  if (ncol > columns_in_input)
    columns_in_input = ncol;

  return iwcut(buffer, word_beginnings, columns_requested, output);
}

#if GCC_VERSION > 50000

static int
identify_random_columns(const const_IWSubstring & buffer,
                        resizable_array<int> & columns_requested)
{
  int nw;
  if (' ' == input_token_separator)
    nw = buffer.nwords();
  else
    nw = buffer.nwords_single_delimiter(input_token_separator);

  if (nw <= random_columns)
  {
    cerr << "identify_random_columns::request " << random_columns << " random columns, but input only contains " << nw << " impossible\n";
    return 0;
  }

  int * tmp = new int[nw]; std::unique_ptr<int[]> free_tmp(tmp);

  std::fill_n(tmp, nw, 0);

  for (int i = 0; i < random_columns; ++i)
  {
    tmp[i] = 1;
  }

  std::random_device rd;
  std::mt19937 rng(rd());

  std::shuffle(tmp + 1, tmp + nw - 1, rng);   // leave the first one set to 1 - presumably most common case

  columns_requested.resize_keep_storage(0);

  for (int i = 0; i < nw; ++i)
  {
    if (tmp[i])
      columns_requested.add(i);
  }

  return 1;
}
#endif

static int
matches_except_for_quotes(const IWString & descriptor,
                          const const_IWSubstring & d)
{
  if (! d.starts_with('"'))
    return 0;

  if (! d.ends_with('"'))
    return 0;

  const int istop = d.length() - 1;
  for (int i = 1; i < istop; ++i)
  {
    if (descriptor[i-1] != d[i])
      return 0;
  }

  return 1;
}

static int
find_column_number(const IWString & descriptor,
                   const resizable_array_p<const_IWSubstring> & header,
                   resizable_array<int> & columns_requested)
{
  int n = header.number_elements();

  const_IWSubstring noquotes;
  if (',' == input_token_separator && noquotes.starts_with('"') && descriptor.ends_with('"'))
  {
    descriptor.from_to(1, descriptor.length() - 2, noquotes);
  }

  for (int i = 0; i < n; i++)
  {
    const const_IWSubstring & d = *(header[i]);

//  cerr << "Comparing '" << descriptor << "' with '" << d << "'\n";

    if (descriptor == d)
      ;
    else if (ignore_case && descriptor.equals_ignore_case(d))
      ;
    else if (',' == input_token_separator && matches_except_for_quotes(descriptor, d))
      ;
    else
      continue;

    columns_requested.add_if_not_already_present(i);

    if (verbose > 1)
      cerr << "Descriptor '" << d << " in column " << (i + 1) << endl;

    return 1;
  }

  cerr << "Strange, no column matches '" << descriptor << "'\n";
  if (verbose)
  {
    cerr << "Found these descriptors\n";
    for (int i = 0; i < n; ++i)
    {
      cerr << " col " << i << " desc '" << *(header[i]) << "'\n";
    }
  }

  return 0;
}

static int
identify_column (const IWString & descriptor,
                 const resizable_array_p<const_IWSubstring> & header,
                 resizable_array<int> & columns_requested)
        
{
  if (! match_descriptor_names_as_regular_expressions)
    return find_column_number(descriptor, header, columns_requested);

  IW_Regular_Expression rx;

  if (! rx.set_pattern(descriptor))
  {
    cerr << "INvalid regular expression '" << descriptor << "'\n";
    return 0;
  }

  int n = header.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    const const_IWSubstring & d = *(header[i]);

//  cerr << "Comparing '" << descriptor << "' with '" << d << "'\n";

    if (rx.matches(d))
    {
      columns_requested.add_if_not_already_present(i);
      rc++;

      if (verbose > 1)
        cerr << "Descriptor '" << d << " in column " << (i + 1) << endl;
    }
  }

  return rc;
}

static int
do_split(const const_IWSubstring & buffer,
         const char input_token_separator,
         resizable_array_p<const_IWSubstring> & word_beginnings)
{
  const int n = buffer.length();

  if (0 == n)
    return 0;

  const char dquote = '"';

//cerr << "Looking for wb in '" << buffer << endl;

  int previous_delimiter = -1;    // 

  int in_quote;
  if (input_token_separator == buffer[0])
  {
    word_beginnings.add(new const_IWSubstring());
    in_quote = 0;
    previous_delimiter = 0;
  }
  else
    in_quote = dquote == buffer[0];

  for (int i = 1; i < n; ++i)
  {
    if (dquote == buffer[i])
      in_quote = ! in_quote;
    else if (in_quote)
      ;
    else if (input_token_separator == buffer[i])
    {
      if (in_quote)
        ;
      else 
      {
        word_beginnings.add(new const_IWSubstring(buffer.rawchars() + previous_delimiter + 1, i - previous_delimiter - 1));
        previous_delimiter = i;
      }
    }
  }

  word_beginnings.add(new const_IWSubstring(buffer.rawchars() + previous_delimiter + 1, n - previous_delimiter - 1));

  if (2 == input_is_quoted_tokens)
  {
    const char dquote = '"';

    for (int i = 0; i < word_beginnings.number_elements(); ++i)
    {
      auto * x = word_beginnings[i];

      if (x->starts_with(dquote) && x->ends_with(dquote))
      {
        x->remove_leading_chars(1);
        x->chop();
      }
    }
  }

  return word_beginnings.number_elements();
}

/*
  Convert descriptor names to columns
*/

static int
determine_descriptors_to_be_output(const const_IWSubstring & buffer,
                                   resizable_array_p<IWString> & descriptors_requested,
                                   resizable_array<int> & columns_requested)
{
  resizable_array_p<const_IWSubstring> header;

  if (input_is_quoted_tokens)
    do_split(buffer, input_token_separator, header);
  else
    buffer.split(header, input_token_separator);

#ifdef DEBUG_IWCUT
  cerr << "header split into " << header.size() << " items\n";
  for (int i = 0; i < header.number_elements(); ++i)
  {
    cerr << " col " << i << " dname " << *header[i] << endl;
  }
#endif

  if (' ' != descriptor_name_remove_to)
  {
    for (int i = 0; i < header.number_elements(); i++)
    {
      const_IWSubstring & d = *(header[i]);

      if (d.contains(descriptor_name_remove_to))
        d.remove_up_to_first(descriptor_name_remove_to);
    }
  }

  if (add_first_token_when_processing_descriptors)
    columns_requested.add_if_not_already_present(0);

  int nr = descriptors_requested.number_elements();

  int rc = nr;

  for (int i = 0; i < nr; i++)
  {
    const IWString & d = *(descriptors_requested[i]);

    if (identify_column(d, header, columns_requested))
      continue;

    cerr << "Cannot find descriptor matching '" << d << "' in header\n";

// Note that we temporarily reverse the meanings of each and header. They will be swapped after the first record in the main loop

    if (missing_column_fill.length())
    {
      extra_stuff_each_record << output_token_separator << d;
      extra_stuff_header_record << missing_column_fill;
      continue;
    }

    cerr << buffer << endl;
    rc = 0;
  }

  if (verbose)
    cerr << "Selected " << columns_requested.number_elements() << " descriptors for output\n";

  return rc;
}

static int
iwcut(iwstring_data_source & input,
      resizable_array<int> & columns_requested,
      resizable_array_p<IWString> & descriptors_requested,
      IWString_and_File_Descriptor & output)
{
  static int first_record = 1;

  const_IWSubstring buffer;

  IWString tmp;     // needs to be kept in scope. Could go in the whole loop, but very inefficient to instantiate it each iteration but only use it on the first record

  while (input.next_record(buffer))
  {
    int first_record_this_file = (1 == input.lines_read());

    if (first_record_this_file && descriptors_requested.number_elements())
    {
      if (1 == ignore_case)
      {
        tmp = buffer;
        tmp.to_lowercase();
        buffer = tmp;
      }

      if (! determine_descriptors_to_be_output(buffer, descriptors_requested, columns_requested))
        return 0;

      if (first_record)   // first file being processed
        first_record = 0;
      else                // don't write out the header again
        continue;
    }
#if GCC_VERSION > 50000
    else if (random_columns)
    {
      if (! identify_random_columns(buffer, columns_requested))
        return 0;
    }
#endif

    if (suppress_selected_columns && first_record_this_file && columns_requested.number_elements())
    {
      if (' ' == input_token_separator)
        invert_selections(columns_requested, buffer.nwords(), descriptors_requested.number_elements());
      else
        invert_selections(columns_requested, buffer.nwords_single_delimiter(input_token_separator), descriptors_requested.number_elements());
    }

    if (! iwcut(buffer, columns_requested, output))
    {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      return 0;
    }

    if (first_record_this_file && extra_stuff_header_record.length())
    {
      extra_stuff_each_record = extra_stuff_header_record;
      extra_stuff_header_record.resize(0);
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();

  return 1;
}

static int
iwcut(const char * fname,
      resizable_array<int> & columns_requested,
      resizable_array_p<IWString> & descriptors_requested,
      IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  if (token_regexp.active())
    return iwcut_token_regexp(input, columns_requested, output);

  return iwcut(input, columns_requested, descriptors_requested, output);
}

static int
iwcut (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vf:F:d:D:mxcrz:R:i:o:K:utq:a:E:bg:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('a'))
  {
    if (! cl.value('a', random_columns) || random_columns < 1)
    {
      cerr << "The number of random columns (-a) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will select " << random_columns << " random columns\n";
  }
  else if (! cl.option_present('f') && ! cl.option_present('d') && ! cl.option_present('D') && ! cl.option_present('F') && ! cl.option_present('R'))
  {
    cerr << "Must specify which columns to omit via the -f, -d or -D options\n";
    usage(4);
  }

  if (cl.option_present('i') && cl.option_present('t'))
  {
    cerr << "Sorry, cannot use both -i and -t options\n";
    usage(3);
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');

    if (! char_name_to_char(i))
    {
      cerr << "Unrecognised input token specifier '" << i << "'\n";
      return 1;
    }

    input_token_separator = i[0];

    if (verbose)
      cerr << "Input token separator '" << input_token_separator << "'\n";
  }
  else if (cl.option_present('t'))
  {
    input_token_separator = '\t';

    if (verbose)
      cerr << "Input token separator set to tab\n";
  }

  if (cl.option_present('q'))
  {
    const_IWSubstring q = cl.string_value('q');

    if ("rm" == q)
    {
      input_is_quoted_tokens = 2;
    }
    else if ("k" == q || "keep" == q)
    {
      input_is_quoted_tokens = 1;
      if (verbose)
        cerr << "Input is quoted tokens\n";
    }
    else
    {
      cerr << "Unrecognised -q qualifier '" << q << "'\n";
      usage(1);
    }
  }

  if (cl.option_present('o'))
  {
    IWString o = cl.string_value('o');
    if (! char_name_to_char(o))
    {
      cerr << "Sorry, the output token separator can only be a single character\n";
      usage(5);
    }

    output_token_separator = o[0];

    if (verbose)
      cerr << "Output token separator '" << output_token_separator << "'\n";
  }

  if (cl.option_present('u'))
  {
    consecutive_delimiters_mean_empty_words_between = 1;
    if (verbose)
      cerr << "Consecutive delimiters mean empty words in between adjacent delimiters\n";
  }

  if (' ' != input_token_separator)
    consecutive_delimiters_mean_empty_words_between = 1;

  if (cl.option_present('g'))
  {
    const_IWSubstring g = cl.string_value('g');

    if (1 != g.length())
    {
      cerr << "The space replace option (-g) must specify a single character, '" << g << "' invalid\n";
      return 1;
    }

    gsub_spaces_in_tokens = g[0];

    if (verbose)
      cerr << "Whilespace inside tokens translated to '" << gsub_spaces_in_tokens << "'\n";
  }

  if (cl.option_present('c'))
  {
    ignore_case = cl.option_count('c');
    if (verbose)
      cerr << "Will ignore case when comparing descriptor names\n";
  }

  if (cl.option_present('r'))
  {
    match_descriptor_names_as_regular_expressions = 1;
    if (verbose)
      cerr << "Will match descriptors as regular expressions\n";
  }

  resizable_array<int> columns_requested;

  if (cl.option_present('f'))
  {
    if (! parse_dash_f_option(cl, 'f', columns_requested))
    {
      cerr << "Invalid -f specifiers\n";
      usage(3);
    }

    if (verbose)
    {
      cerr << "Will extract these columns\n";
      for (int i = 0; i < columns_requested.number_elements(); i++)
      {
        cerr << ' ' << (columns_requested[i] + 1);
      }
      cerr << endl;
    }
  }

  if (cl.option_present('E'))
  {
    if (1 == cl.option_count('E'))
    {
      const_IWSubstring e = cl.string_value('E');

      if ('.' == e)
        all_fields_must_have_content = 1;

      if (verbose)
        cerr << "Will discard any record with an empty selected field\n";
    }

    if (! all_fields_must_have_content)
    {
      resizable_array<int> tmp;
      if (! parse_dash_f_option(cl, 'E', tmp))
      {
        cerr << "Invalid -E specifiers\n";
        usage(1);
      }

      for (int i = 0; i < tmp.number_elements(); ++i)
      {
        columns_that_must_not_be_empty[tmp[i]] = 1;
      }

      if (verbose)
        cerr << "Defined " << tmp.number_elements() << " columns that must not be empty\n";
    }
  }

// Make sure the -m option is examined before the -F and -D options

  if (cl.option_present('m'))
  {
    descriptor_names_on_one_line = 0;

    if (verbose)
      cerr << "Descriptor names one per line in -D file(s), columns one per line in -F file(s)\n";
  }

  if (cl.option_present('F'))
  {
    if (! parse_dash_F_option(cl, 'F', columns_requested))
    {
      cerr << "Invalid -F combinations\n";
      usage(6);
    }
  }

  resizable_array_p<IWString> descriptors_requested;
  if (cl.option_present('d'))
  {
    if (! parse_dash_d_option(cl, 'd', descriptors_requested))
    {
      cerr << "INvalid -d combinations\n";
      usage(5);
    }
  }

  if (cl.option_present('D'))
  {
    if (! parse_dash_D_option(cl, 'D', descriptors_requested))
    {
      cerr << "Invalid -D combinations\n";
      usage(5);
    }
  }

  if (cl.option_present('b'))
  {
    add_first_token_when_processing_descriptors = 0;

    if (verbose)
      cerr << "Will not automatically add first column\n";
  }

  if (cl.option_present('K'))
  {
    IWString tmp;

    cl.value('K', tmp);

    assert (tmp.length() > 0);

    descriptor_name_remove_to = tmp[0];

    if (verbose)
      cerr << "Will remove characters up to " << descriptor_name_remove_to << " from descriptor names\n";

    for (int i = 0; i < descriptors_requested.number_elements(); i++)
    {
      IWString & d = *(descriptors_requested[i]);

      if (d.contains(descriptor_name_remove_to))
        d.remove_up_to_first(descriptor_name_remove_to);
    }
  }

  if (cl.option_present('R'))
  {
    if (cl.option_count('R') > 1)
    {
      cerr << "Only one -R option is supported. Create a component regexp '(XXX|YYY|ZZZ)'\n";
      return 4;
    }

    const_IWSubstring r = cl.string_value('R');

    if (! token_regexp.set_pattern(r))
    {
      cerr << "Invalid token regular expression '" << r << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Columns matching '" << token_regexp.source() << "' will be output\n";
  }

  if (! ignore_case)
    ;
  else if (0 == descriptors_requested.number_elements())
  {
    cerr << "Very strange, ignoring case, but no descriptors requested for selection\n";
  }
  else    // convert the descriptor names to lowercase
  {
//  for (int i = 0; i < descriptors_requested.number_elements (); i++)
//  {
//    descriptors_requested[i]->to_lowercase();
//  }
  }

  if (verbose && descriptors_requested.number_elements())
  {
    cerr << "Will extract these descriptors\n";
    for (int i = 0; i < descriptors_requested.number_elements(); i++)
    {
      cerr << ' ' << *(descriptors_requested[i]) << endl;
    }
  }

  if (verbose && columns_requested.number_elements())
  {
    cerr << "Will extract these columns\n";
    for (int i = 0; i < columns_requested.number_elements(); i++)
    {
      cerr << ' ' << (columns_requested[i] + 1) << endl;
    }
  }

  if (cl.option_present('x'))
  {
    suppress_selected_columns = 1;

    if (verbose)
      cerr << "Will suppress selected columns\n";
  }

  if (cl.option_present('z'))
  {
    missing_column_fill = cl.string_value('z');

    if (verbose)
      cerr << "Missing column(s) filled with '" << missing_column_fill << "'\n";

    missing_column_fill.insert_before(0, output_token_separator);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! iwcut(cl[i], columns_requested, descriptors_requested, output))
    {
      cerr << "iwcut: error processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    if (records_discarded_for_zero_field)
      cerr << records_discarded_for_zero_field << " records discarded for empty field\n";
  }

  return rc;
}
int
main (int argc, char ** argv)
{
  int rc = iwcut(argc, argv);

  return rc;
}
