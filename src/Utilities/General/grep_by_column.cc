/*
  Scans a descriptor file for similarity to a given vector
*/

#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static char input_separator = ' ';
static char output_separator = ' ';

static int just_matched_columns = 0;

static int lines_read = 0;
static int lines_written = 0;

static int first_column_contains_identifiers = 0;

static int min_matches_needed = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "grep like utility but with range of search restricted to specific columns\n";
  cerr << " -e <col>:string   search for <string> in <col>.\n";
  cerr << "                     <col> can be column number of name. Separate list with commas\n";
  cerr << "                     -e 3,foo,bar:hello    search for 'hello' in columns 3 or 'foo' or 'bar'\n";
  cerr << " -E <col>:<rx>     search for regexp <rx> in col <col>\n";
  cerr << " -W <col>          extra columns to write without checking\n";
  cerr << " -m <matches>      will write record if <matches> of the queries match\n";
  cerr << " -I                first column contains identifiers\n";
  cerr << " -s                only write the selected columns\n";
  cerr << " -i <sep>          input  column separator\n";
  cerr << " -o <sep>          output column separator\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}


#define SPECIAL_VALUE_ALWAYS_WRITE -2

class Column_Specification
{
  protected:
    resizable_array<int> _column_numbers;
    resizable_array_p<IWString> _column_names;
  public:
    Column_Specification() {}

    int build(const const_IWSubstring &);

    int  column_described_by_name() const { return _column_names.number_elements();}

    int identify_column_numbers(int * process_column, const int ncols, const int ndx) const;
    int identify_column_names(const resizable_array_p<IWString> & header, int * process_column, const int ncols, const int ndx) const;


    int set_process_columns(int * t, const int flag) const;
};

int
Column_Specification::build(const const_IWSubstring & buffer)
{
  int i = 0;
  const_IWSubstring token;

  while(buffer.nextword(token, i, ','))
  {
    int c;
    if (! token.numeric_value(c) || c < 1)
      _column_names.add(new IWString(token));
    else
      _column_numbers.add(c - 1);
  }
  
  return 1;
}

int
Column_Specification::set_process_columns(int * t, const int flag) const
{
  return 1;
}

/*
  For each -e and -E, we build a column specification object. they can then examine a header
  and determine the columns
*/

class Column_Match_Specification : public Column_Specification
{
  private:

    IWString _column_contents;
    int _is_regexp;
    std::unique_ptr<re2::RE2> _rx;

  public:
    Column_Match_Specification();

    void debug_print(std::ostream &) const;

    int  build(const const_IWSubstring & a, const const_IWSubstring & b);

    int set_is_regexp(const int s);

    int matches(const const_IWSubstring & token);
};

Column_Match_Specification::Column_Match_Specification()
{
  _is_regexp = 0;

  return;
}

int
Column_Match_Specification::build(const const_IWSubstring & a,
                                  const const_IWSubstring & b)
{
  if (! this->Column_Specification::build(a))
  {
    cerr << "Column_Match_Specification::build:invalid column specification '" << a << "'\n";
    return 0;
  }

  _column_contents = b;

  return 1;
}

class Set_of_Column_Match_Specification
{
  private:
    resizable_array_p<Column_Match_Specification> _column_specifications;

    int _ncols;    // will be set to the smallest number of columns we need to check

    int * _process_column;

    resizable_array_p<Column_Specification> _extra_columns_to_write;

//  private functions

    int _identify_columns_to_process(const IWString & first_record);
    int _identify_columns_to_write(const IWString & first_record);
    int _write_selected_columns_header(const IWString & headr, IWString_and_File_Descriptor & output) const;
    int _any_column_specifications_are_names() const;
    int _process(iwstring_data_source &, IWString_and_File_Descriptor & output);
    int _process_record(const const_IWSubstring & line, IWString_and_File_Descriptor & output);

  public:
    Set_of_Column_Match_Specification();
    ~Set_of_Column_Match_Specification();

    void debug_print(std::ostream & output) const;

    int build(const Command_Line & cl, const char flag);
    int build_extra_columns(const Command_Line & cl, const char flag);

    int process(iwstring_data_source & input, IWString_and_File_Descriptor & output);
};

Set_of_Column_Match_Specification::Set_of_Column_Match_Specification()
{
  _ncols = 0;

  _process_column = nullptr;

  return;
}

Set_of_Column_Match_Specification::~Set_of_Column_Match_Specification()
{
  if (nullptr != _process_column)
    delete [] _process_column;
}


void
Set_of_Column_Match_Specification::debug_print(std::ostream & output) const
{
  for (int i = 0; i < _ncols; ++i)
  {
    output << " _process_column[" << i << "] " << _process_column[i] << endl;
  }

  return;
}

static int
build_column_specifications(const Command_Line & cl,
                            const char flag,
                            resizable_array_p<Column_Match_Specification> & column_specifications)
{
  const_IWSubstring e;
  for (int i = 0; cl.value(flag, e, i); ++i)
  {
    const_IWSubstring a, b;

    if (! e.split(a, ':', b) || 0 == a.length())
    {
      cerr << "Set_of_Column_Match_Specification::build:invalid column/contents specification '" << e << "' must be of the form 'col:contents'\n";
      return 0;
    }

    Column_Match_Specification * t = new Column_Match_Specification();

    if (! t->build(a, b))
    {
      cerr << "Set_of_Column_Match_Specification::build:invalid column specification '" << e << "'\n";
      delete t;
      return 0;
    }

    if ('E' == flag)
      t->set_is_regexp(1);

    column_specifications.add(t);
   }

   return 1;
}

int
Set_of_Column_Match_Specification::build(const Command_Line & cl,
                                    const char flag)
{
  return build_column_specifications(cl, flag, _column_specifications);
}

int
Set_of_Column_Match_Specification::build_extra_columns(const Command_Line & cl,
                                    const char flag)
{
  const_IWSubstring s;
  for(int i = 0; cl.value(flag, s, i); ++i)
  {
    Column_Specification * t = new Column_Specification();
    if (! t->build(s))
    {
      delete t;
      cerr << "Set_of_Column_Match_Specification::build_extra_columns:invalid column specifier '" << s << "'\n";
      return 0;
    }

    _extra_columns_to_write.add(t);
  }

  return 1;
}

int
Set_of_Column_Match_Specification::process(iwstring_data_source & input,
                                      IWString_and_File_Descriptor & output)
{
  IWString first_record;

  if (! input.next_record(first_record))
  {
    cerr << "Set_of_Column_Match_Specification::build:empty file\n";
    return 0;
  }

  if (' ' == input_separator)
    _ncols = first_record.nwords();
  else
    _ncols = first_record.nwords_single_delimiter(input_separator);

  _process_column = new int[_ncols];
  std::fill_n(_process_column, _ncols, -1);

  if (! _identify_columns_to_process(first_record))
  {
    cerr << "Set_of_Column_Match_Specification::build:cannot identify columns to process\n";
    cerr << first_record << "\n";
    return 0;
  }

  if (_extra_columns_to_write.number_elements())
  {
    if (! _identify_columns_to_write(first_record))
    {
      cerr << "Set_of_Column_Match_Specification::process:cannot discern extra columns to write\n";
      return 0;
    }
  }

  int last_column = 0;
  for (int i = 0; i < _ncols; ++i)
  {
    if (_process_column[i] >= 0)
      last_column = i;
    else if (SPECIAL_VALUE_ALWAYS_WRITE == _process_column[i])
      last_column = i;
  }

  _ncols = last_column + 1;

  if (verbose > 1)
    debug_print(cerr);

  if (_any_column_specifications_are_names())
  {
    if (just_matched_columns)
      _write_selected_columns_header(first_record, output);
    else
      output << first_record << '\n';
  }
  else
    input.push_record();

  return _process(input, output);
}

static int 
get_next_token(const const_IWSubstring & buffer,
               const char input_separator,
               const_IWSubstring & token,
               int & i)
{
  if (' ' == input_separator)
    return buffer.nextword(token, i);
  else
    return buffer.nextword_single_delimiter(token, i, input_separator);
}

int
Set_of_Column_Match_Specification::_write_selected_columns_header(const IWString & header, 
                                            IWString_and_File_Descriptor & output) const
{
  int i = 0;
  const_IWSubstring token;

  int need_separator = 0;

  for (int c = 0; get_next_token(header, input_separator, token, i) && c < _ncols; ++c)
  {
    if (0 == c && first_column_contains_identifiers && -1 == _process_column[0])
    {
      output << token;
      need_separator = 1;
      continue;
    }

    if (_process_column[c] >= 0)
      ;
    else if (SPECIAL_VALUE_ALWAYS_WRITE == _process_column[c])
      ;
    else
      continue;

    if (need_separator)
      output << output_separator;
     else
       need_separator = 1;

    output << token;
  }

  output << '\n';

  return 1;
}

void
Column_Match_Specification::debug_print(std::ostream & output) const
{
  output << "Column_Match_Specification::debug_print:";
  if (_column_numbers.number_elements())
  {
    output << "column numbers";
    for (int i = 0; i < _column_numbers.number_elements(); ++i)
    {
      output << ' ' << _column_numbers[i];
    }
  }
  else
  {
    output << "columns";
    for (int i = 0; i < _column_names.number_elements(); ++i)
    {
      output << ' ' << *_column_names[i];
    }
  }
  output << "\n";

  output << "Contents ";
  if (_is_regexp)
    output << "regexp '" << _rx->pattern();
  else
    output << "string '" << _column_contents;
    
  output << "'\n";

  return;
}

int
Column_Match_Specification::matches(const const_IWSubstring & s)
{
//cerr << "Column_Match_Specification::matches what about " << s << " contents " << _column_contents << "'\n";
  if (_is_regexp)
    return iwre2::RE2PartialMatch(s, *_rx);

  return s == _column_contents;
}

int
Column_Match_Specification::set_is_regexp(const int s)
{
  _is_regexp = s;

  if (! _is_regexp)     // shuold check if _rx is initialised. Will not happen here...
    return 1;

  if (! iwre2::RE2Reset(_rx, _column_contents))
  {
    cerr << "Column_Match_Specification::set_is_regexp:invalid regexp " << _column_contents << "'\n";
    return 0;
  }

  return 1;
}

int
Column_Specification::identify_column_numbers(int * process_column, const int ncols,
                                              const int ndx) const
{
  for (int i = 0; i < _column_numbers.number_elements(); ++i)
  {
    const auto c = _column_numbers[i];
//  cerr << "What about column " << c << endl;

    if (c >= ncols)
    {
      cerr << "Column_Match_Specification::identify_numeric_column_numbers:out of range " << c << " max is " << ncols << endl;
      return 0;
    }

    process_column[c] = ndx;
  }

  return 1;
}

int
Column_Specification::identify_column_names(const resizable_array_p<IWString> & header,
                                            int * process_column,
                                            const int ncols,
                                            const int ndx) const
{
  for (int c = 0; c < _column_names.number_elements(); ++c)
  {
    int found_match = 0;
    const IWString & si = *_column_names[c];
//  cerr << "looking for match to '" << si << "'\n";

    for (int j = 0; j < header.number_elements(); ++j)
    {
//    cerr << "   does it match " << *header[j] << endl;
      if (si != *header[j])
        continue;

      process_column[j] = ndx;
      found_match = 1;
      break;
    }

    if (! found_match)
    {
      cerr << "Column_Match_Specification::identify_column_names:no match for " << si << endl;
      if (verbose)
      {
        for (int i = 0; i < header.number_elements(); ++i)
        {
          cerr << " col " << (i+1) << ' ' << *header[i] << endl;
        }
      }
      return 0;
    }
  }

  return 1;
}

static void
tokenise(const IWString & s,
         const char input_separator,
         resizable_array_p<IWString> & f)
{
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; get_next_token(s, input_separator, token, i); ++col)
  {
    f.add(new IWString(token));
  }
}

int
Set_of_Column_Match_Specification::_identify_columns_to_process(const IWString & first_record)
{
  int need_to_check_names = 0;

//cerr << "Checking " << _column_specifications.number_elements() << " column specifications\n";

  for (int i = 0; i < _column_specifications.number_elements(); ++i)
  {
    if (! _column_specifications[i]->identify_column_numbers(_process_column, _ncols, i))
    {
      cerr << "Set_of_Column_Match_Specification::_identify_columns_to_process:specifier " << i << " failed numeric identification, ncols " << _ncols << endl;
      return 0;
    }

    if (_column_specifications[i]->column_described_by_name())
      need_to_check_names = 1;
  }

//cerr << need_to_check_names << " need_to_check_names\n";

  if (! need_to_check_names)
    return 1;

  resizable_array_p<IWString> header;
  tokenise(first_record, input_separator, header);

  for (int i = 0; i < _column_specifications.number_elements(); ++i)
  {
    if (! _column_specifications[i]->identify_column_names(header, _process_column, _ncols, i))
    {
      return 0;
    }
  }

  return 1;
}

int
Set_of_Column_Match_Specification::_identify_columns_to_write(const IWString & first_record)
{
  int need_to_check_names = 0;

//cerr << "Checking " << _column_specifications.number_elements() << " column specifications\n";

  int * tmp = new int[_ncols]; std::unique_ptr<int[]> free_tmp(tmp);
  std::fill_n(tmp, _ncols, -1);

  for (int i = 0; i < _extra_columns_to_write.number_elements(); ++i)
  {
    if (! _extra_columns_to_write[i]->identify_column_numbers(tmp, _ncols, i))
    {
      cerr << "Set_of_Column_Match_Specification::_identify_columns_to_process:specifier " << i << " failed numeric identification, ncols " << _ncols << endl;
      return 0;
    }

    if (_extra_columns_to_write[i]->column_described_by_name())
      need_to_check_names = 1;
  }

//cerr << need_to_check_names << " need_to_check_names\n";

  if (need_to_check_names)
  {
    resizable_array_p<IWString> header;
    tokenise(first_record, input_separator, header);

    for (int i = 0; i < _extra_columns_to_write.number_elements(); ++i)
    {
      if (! _extra_columns_to_write[i]->identify_column_names(header, tmp, _ncols, i))
      {
        return 0;
      }
    }
  }

  for (int i = 0; i < _ncols; ++i)
  {
    if (tmp[i] < 0)
      continue;

    if (_process_column[i] >= 0)
    {
      cerr << "Set_of_Column_Match_Specification::_identify_columns_to_write:column " << i << " set for searching\n";
      continue;
    }

    _process_column[i] = SPECIAL_VALUE_ALWAYS_WRITE;
  }

  return 1;
}

int
Set_of_Column_Match_Specification::_any_column_specifications_are_names() const
{
  for (int i = 0; i < _column_specifications.number_elements(); ++i)
  {
    if (_column_specifications[i]->column_described_by_name())
      return 1;
  }

  return 0;
}

static void
append_to_buffer(const const_IWSubstring & token,
                 bool & need_separator,
                 IWString & buffer)
{
  if (! need_separator)
    need_separator = true;
  else
    buffer << output_separator;

  buffer << token;

  return;
}

int
Set_of_Column_Match_Specification::_process(iwstring_data_source & input,
                                       IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    lines_read++;

    if (! _process_record(buffer, output))
    {
      cerr << "Set_of_Column_Match_Specification::_process_record:fatal error processing\n";
      cerr << buffer << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Set_of_Column_Match_Specification::_process_record(const const_IWSubstring & line,
                                              IWString_and_File_Descriptor & output)
{
  bool need_separator = false;

  IWString buffer;
  buffer.resize(line.length());

  const_IWSubstring id_col_0;

  int matches_found = 0;

  int i = 0;
  const_IWSubstring token;
  for (int c = 0; get_next_token(line, input_separator, token, i) && c < _ncols; ++c)
  {
    if ((0 == c && first_column_contains_identifiers && -1 == _process_column[0]) || (SPECIAL_VALUE_ALWAYS_WRITE == _process_column[c]))
    {
      if (just_matched_columns)
        append_to_buffer(token, need_separator, buffer);
      continue;
    }

    if (_process_column[c] < 0)
      continue;

    if (! _column_specifications[_process_column[c]]->matches(token))
    {
      if (min_matches_needed > 0)
      {
        if (just_matched_columns)
          append_to_buffer(token, need_separator, buffer);
        continue;
      }
      return 1;
    }

    matches_found++;
    if (just_matched_columns)
      append_to_buffer(token, need_separator, buffer);
  }

  if (min_matches_needed > 0 && matches_found < min_matches_needed)
    return 1;

  if (just_matched_columns)
    output << buffer;
  else
    output << line;

  output << '\n';

  lines_written++;

  return 1;
}

static int
grep_by_column(const char * fname,
               Set_of_Column_Match_Specification & socs,
               IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return socs.process(input, output);
}

static int
grep_by_column (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "ve:E:Isi:o:W:m:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  Set_of_Column_Match_Specification socs;

  if (! cl.option_present('e') && ! cl.option_present('E'))
  {
    cerr << "Must specify column(s) to fetch via the -e/-E option(s)\n";
    usage(1);
  }

  if (! socs.build(cl, 'e') || ! socs.build(cl, 'E'))
  {
    cerr << "Invalid column specifications\n";
    return 1;
  }

  if (cl.option_present('W'))
  {
    if (! socs.build_extra_columns(cl, 'W'))
    {
      cerr << "Cannot parse extra columns to write (-W) specifications\n";
      return 1;
    }

    just_matched_columns = 1;
  }

  if (cl.option_present('I'))
  {
    first_column_contains_identifiers = 1;

    if (verbose)
      cerr << "First column contains identifiers\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', min_matches_needed) || min_matches_needed < 1)
    {
      cerr << "The minimum number of matches needed (-m) must be a +ve whole number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will write a record if at least " << min_matches_needed << " conditions matched\n";
  }

  if (cl.option_present('s'))
  {
    just_matched_columns = 1;

    if (verbose)
      cerr << "Will just write selected columns\n";
  }

  if (cl.option_present('i'))
  {
    IWString tmp;
    cl.value('i', tmp);

    if (! char_name_to_char(tmp))
    {
      cerr << "The word delimiter must be a single character, '" << tmp << "' is invalid\n";
      usage(4);
    }

    input_separator = tmp[0];

    if (verbose)
      cerr << "input delimiter '" << input_separator << "'\n";
  }

  if (cl.option_present('o'))
  {
    IWString tmp;
    cl.value('o', tmp);

    if (! char_name_to_char(tmp))
    {
      cerr << "The output delimiter must be a single character, '" << tmp << "' is invalid\n";
      usage(4);
    }

    output_separator = tmp[0];

    if (verbose)
      cerr << "Output delimiter '" << output_separator << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! grep_by_column(cl[i], socs, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "wrote " << lines_written << " of " << lines_read << " records\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = grep_by_column(argc, argv);

  return rc;
}
