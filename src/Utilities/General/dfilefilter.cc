#include <stdlib.h>
#include <ctype.h>

#include "cmdline.h"
#include "set_or_unset.h"
#include "misc.h"
#include "iwstring_data_source.h"
#include "logical_expression.h"

static int verbose = 0;

static int records_read = 0;

static int records_written = 0;

static IWString missing_value('.');

static int missing_value_passes = 0;

static int missing_column_ok = 0;

static IW_Logical_Expression logical_expression;

static int just_write_identifiers = 0;

static IWString_and_File_Descriptor stream_for_rejected_records;

/*
  We can save time by only looking at columns that will be needed
*/

static int highest_column_needed = 0;

/*
  When processing multiple files, we only want to write the header
  record one time
*/

static int header_record_written = 0;

/*
  May 2007. Chemexplorer tie-in. The first N columns are reagent types,
  but the last M are molecular properties. We want to be able to filter
  the file to all rows that have just one reagent present - that is, zero
  in all the other reagent rows
*/

static IWString first_non_reagent_column_header;
static int number_columns_reagent_types = 0;

static int all_other_columns_zero = 0;

static char input_separator = ' ';

static int write_passing_records = 1;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Filters a descriptor file based on expressions\n";
  cerr << " -e <expression>  the expression for filtering\n";
  cerr << " -m <text>        missing value string\n";
  cerr << " -p               missing values pass all conditions (by default they fail all)\n";
  cerr << " -c               missing column fails condition, but is non-fatal\n";
  cerr << " -i               just write the identifiers, rather than the whole descriptor file\n";
  cerr << " -B <fname>       write rejected records to <fname>\n";
  cerr << " -K <dname>       Kludge for trxn - first non-reagent type column\n";
  cerr << " -z               all non-selected columns must be zero\n";
  cerr << " -n               no output, just an informational scan\n";
  cerr << " -s <delim>       input delimiter (default ' ')\n";
  cerr << " -v               verbose output\n";

  cerr << endl;
  cerr << "Expressions can be simple: 'clogp<=5.3' or 'nrings !=3' or 'amw>250', or complex\n";
  cerr << "'clogp<=5.3 && amw>200 && numhba<=5 && numhbd<=10' or\n";
  cerr << "'clogp>2.0&&clogp<=5.0||amw<600'\n";
  cerr << endl;
  cerr << "The relational operators '=', '!=', '>', '>=', '<', '<=' are supported\n";
  cerr << "The logical operators '&&' (AND), '||' (OR), '^^' (XOR) and ';;' (low priority AND) are supported\n";
  cerr << endl;
  cerr << "Regular expressions are supported with the =~ operator\n";
  cerr << " 'amw>150 && natoms=~5' which means amw > 150 or natoms containing a '5'\n";

  exit(rc);
}

/*
  For each column, we keep track of the text representation as well as the corresponding numeric
  value
*/

class Column_Data : public const_IWSubstring, public Set_or_Unset<double>
{
  private:
  public:
    void reset () { Set_or_Unset<double>::unset (); const_IWSubstring::make_empty ();}

    int numeric_value (double &);

    void new_value (const const_IWSubstring & rhs);

    int is_zero ();
};

int
Column_Data::numeric_value (double & v)
{
  assert (const_IWSubstring::length() > 0);

  if (Set_or_Unset<double>::is_set())
    return Set_or_Unset<double>::value(v);

  if (! const_IWSubstring::numeric_value(v))
  {
    const_IWSubstring & string_ref = *this;
    cerr << "Column_Data::numeric_value: invalid numeric '" << string_ref << "'\n";
    return 0;
  }

  Set_or_Unset<double>::set(v);

  return 1;
}

int 
Column_Data::is_zero ()
{
  double v;
  if (! Set_or_Unset<double>::is_set())
  {
    if (! const_IWSubstring::numeric_value(v))
    {
      const_IWSubstring & string_ref = *this;
      cerr << "Column_Data::numeric_value: invalid numeric '" << string_ref << "'\n";
      return 0;
    }

    Set_or_Unset<double>::set(v);
  }
  else
    Set_or_Unset<double>::value(v);

  return 0.0 == v;
}

/*
  We are processing a new record and each column gets new data
*/

void
Column_Data::new_value (const const_IWSubstring & rhs)
{
  const_IWSubstring::operator=(rhs);

  Set_or_Unset<double>::unset();

  return;
}

static int
get_next_token (const const_IWSubstring & buffer,
                const_IWSubstring & token,
                int & i)
{
  if (' ' == input_separator)
    return buffer.nextword(token, i);
  else
    return buffer.nextword_single_delimiter(token, i, input_separator);
}

static int
tokenise (const const_IWSubstring & buffer,
          int ncol,
          Column_Data * cd)
{
  const_IWSubstring token;
  int i = 0;
  int col = 0;

  while (get_next_token(buffer, token, i))
  {
    cd[col].new_value(token);

    col++;

    if (col > highest_column_needed)
      return 1;
  }

  if (col != ncol)
  {
    if (!missing_column_ok) 
      cerr << "Column count mismatch, found " << col << " expected " << ncol << endl;
    return 0;
  }

  return 1;
}

/*
  The expression being parsed looks like
  'clogp <= 5&& amw < 500&&xvc>=4'

  We need a class to hold the individual tokens
*/

#define CC_OPERATOR_EQUAL 1
#define CC_OPERATOR_NE 2
#define CC_OPERATOR_LE 3
#define CC_OPERATOR_LT 4
#define CC_OPERATOR_GE 5
#define CC_OPERATOR_GT 6
#define CC_OPERATOR_RX 7

/*
  Unary conditions on columns
*/

class Column_Condition
{
  private:
    IWString _column_name;

//  Once we are passed a header, we can convert column name to column number

    int      _column;

//  The operator

    int _operator;

//  The numeric rhs

    double _value;

//  implement regular expressions sometime

    IW_Regular_Expression _rx;

  public:
    Column_Condition ();

    int debug_print (std::ostream &) const;

    void set_column_name (const_IWSubstring & s) { _column_name = s;}
    const IWString & column_name () const { return _column_name;}

    int column () const { return _column; }

    int build (const const_IWSubstring &);

    int determine_column_number (const const_IWSubstring &);

    int passes (Column_Data * cd);
};

Column_Condition::Column_Condition ()
{
  _column = -1;

  _operator = -1;

  return;
}

/*
  Parse an expression of the kind 'foo == 3'

  Actually this is kind of ugly. There are many different kinds of valid descriptor file
  headers. For example

  Name SVMFP:cRPH_LC50 SVMFP:cRPH_LC50.nonlethal SVMFP:cRPH_LC50.lethal
  Name SVMFP:y SVMFP:y.>80uM SVMFP:y.<80uM

  This makes parsing a column selection very difficult.
*/

int
Column_Condition::build (const const_IWSubstring & e)
{
  IWString expr(e);    // make a copy we can change

  expr.strip_leading_blanks();

  if (0 == expr.length())
  {
    cerr << "Column_Condition::build:empty token, cannot build\n";
    return 0;
  }

  int i = 0;    // our index into expr

  _column_name.resize(expr.length());

  while (isalnum(expr[i]) || '_' == expr[i] || '-' == expr[i] || ':' == expr[i] || '.' == expr[i])
  {
    _column_name += expr[i];
    i++;

    if (i == expr.length())
    {
      cerr << "Column_Condition::build: no operator found in '" << e << "'\n";
      return 0;
    }
  }

  if (expr.length() < 2)    // must be operator and numeric left
  {
    cerr << "Column_Condition::build: invalid expression (rhs too short) '" << e << "'\n";
    return 0;
  }

// what about the FORTRAN type operator specifications

  const int four_character_operator_possible = (expr.length() > 4);

// Now determine the operator

  int nchars = 1;    // the number of chars in the operaor

  if ('=' == expr[i])
  {
    _operator = CC_OPERATOR_EQUAL;

    if ('~' == expr[i + 1])
    {
      _operator = CC_OPERATOR_RX;
      nchars = 2;
    }
    else if ('=' == expr[i + 1])
      nchars = 2;
  }
  else if ('<' == expr[i])
  {
    if ('=' == expr[i + 1])
    {
      nchars = 2;
      _operator = CC_OPERATOR_LE;
    }
    else
      _operator = CC_OPERATOR_LT;
  }
  else if ('>' == expr[i])
  {
    if ('=' == expr[i + 1])
    {
      nchars = 2;
      _operator = CC_OPERATOR_GE;
    }
    else
      _operator = CC_OPERATOR_GT;
  }
  else if ('!' == expr[i] && '=' == expr[i + 1])
  {
    nchars = 2;
    _operator = CC_OPERATOR_NE;
  }
  else if (four_character_operator_possible && expr.matches_ignore_case(".lt.", 4, i))
  {
    nchars = 4;
    _operator = CC_OPERATOR_LT;
  }
  else if (four_character_operator_possible && expr.matches_ignore_case(".le.", 4, i))
  {
    nchars = 4;
    _operator = CC_OPERATOR_LE;
  }
  else if (four_character_operator_possible && expr.matches_ignore_case(".gt.", 4, i))
  {
    nchars = 4;
    _operator = CC_OPERATOR_GT;
  }
  else if (four_character_operator_possible && expr.matches_ignore_case(".ge.", 4, i))
  {
    nchars = 4;
    _operator = CC_OPERATOR_GE;
  }
  else if (four_character_operator_possible && expr.matches_ignore_case(".eq.", 4, i))
  {
    nchars = 4;
    _operator = CC_OPERATOR_EQUAL;
  }
  else if (four_character_operator_possible && expr.matches_ignore_case(".ne.", 4, i))
  {
    nchars = 4;
    _operator = CC_OPERATOR_NE;
  }
  else
  {
    cerr << "Column_Condition::build: unrecognised operator '" << e << "'\n";
    return 0;
  }

  i += nchars;

  if (i >= expr.length())
  {
    cerr << "Column_Condition::build: missing rhs '" << e << "'\n";
    return 0;
  }

  const_IWSubstring v = expr.substr(i);

  if (CC_OPERATOR_RX == _operator)
  {
    if (! _rx.set_pattern(v))
    {
      cerr << "Invalid regular expression '" << v << "'\n";
      return 0;
    }
  }
  else if (! v.numeric_value(_value))
  {
    cerr << "Column_Condition::build: invalid numeric '" << e << "'\n";
    return 0;
  }

  return 1;
}

int
Column_Condition::determine_column_number (const const_IWSubstring & header)
{
  int i = 0;
  int c = 0;
  const_IWSubstring token;
  while (get_next_token(header, token, i))
  {
    if (token == _column_name)
    {
      _column = c;
      return 1;
    }

    c++;
  }

  cerr << "Column_Condition::determine_column_number: bad news, no '" << _column_name << "' in header\n";
  cerr << header << endl;
  return 0;
}

/*
  The only reason this method isn't const is because regular expression
  matching is non-const
*/

int
Column_Condition::passes (Column_Data * cd)
{
  Column_Data & mycd = cd[_column];

  if (missing_value == mycd)
    return missing_value_passes;

  if (CC_OPERATOR_RX == _operator)
  {
    const_IWSubstring p = mycd;

#ifdef DEBUG_RX_MATCHES
    cerr << "Does '" << p << "' match '" << _rx.source() << "'\n";
#endif

    return _rx.matches(p);
  }

  if (0 == mycd.length())   // should we warn about this...
  {
    return missing_value_passes;
  }

  double in_file;
  if (! mycd.numeric_value(in_file))
  {
    cerr << "Column_Condition::passes: cannot extract numeric value\n";
    abort();
  }

//cerr << "Checking " << in_file << " vs " << _value << endl;

  switch (_operator)
  {
    case CC_OPERATOR_EQUAL:
      return in_file == _value;

    case CC_OPERATOR_NE:
      return in_file != _value;

    case CC_OPERATOR_LT:
      return in_file < _value;

    case CC_OPERATOR_LE:
      return in_file <= _value;

    case CC_OPERATOR_GT:
      return in_file > _value;

    case CC_OPERATOR_GE:
      return in_file >= _value;

    default:
      abort();
  }

  return 1;
}

static Column_Condition * column_condition;
static int nc = 0;

static int
_passes_all_filters (Column_Data * cd, int ncol)
{
  logical_expression.reset();

  for (int i = 0; i < nc; i++)
  {
    if (! logical_expression.result_needed(i))    // coupled via an OR with something previously
      continue;

    Column_Condition & c = column_condition[i];

    int r = c.passes(cd);

    logical_expression.set_result(i, r);

    int rc;
    if (logical_expression.evaluate(rc))
      return rc;
  }

  return 1;    // must be OK
}

/*
  If we want the operation "all other columns zero", then we need
  to know which columns to check. These will be the ones that are
  NOT part of a column condition
*/

int * check_column_for_zero = NULL;

static int
initialise_check_column_for_zero (int ncol)
{
  assert (NULL == check_column_for_zero);

//cerr << "check_column_for_zero contains " << ncol << " items\n";

  check_column_for_zero = new_int(ncol, 1);

  check_column_for_zero[0] = 0;     // identifier column

  for (int i = 0; i < nc; i++)
  {
    const Column_Condition & ci = column_condition[i];

    int c = ci.column();

    check_column_for_zero[c] = 0;
  }

  return 1;
}

static int
passes_all_filters (Column_Data * cd, int ncol)
{
  if (! _passes_all_filters(cd, ncol))
    return 0;

  if (NULL == check_column_for_zero)
    return 1;

  int istop = ncol;
  if (number_columns_reagent_types > 0)
    istop = number_columns_reagent_types;

//cerr << "istop set to " << istop << endl;

  for (int i = 1; i < istop; i++)
  {
//  cerr << "Check column " << i << "? " << check_column_for_zero[i] << " zero? " << cd[i].is_zero() << endl;

    if (! check_column_for_zero[i])
      continue;

    if (! cd[i].is_zero())
      return 0;
  }

//cerr << "Passes filter\n";
  return 1;
}

static int
identify_first_non_reagent_column (const const_IWSubstring & header)
{
  int i = 0;
  const_IWSubstring token;

  get_next_token(header, token, i);      // skip over identifier column

  for (int col = 1; get_next_token(header, token, i); col++)
  {
    if (token == first_non_reagent_column_header)
    {
      number_columns_reagent_types = col;

      if (1 == col)
        cerr << "Warning, no non-reagent type columns\n";

      if (verbose)
        cerr << "Column " << col << " contains descriptor '" << first_non_reagent_column_header << "'\n";
      return 1;
    }
  }

  cerr << "Did not find '" << first_non_reagent_column_header << "' descriptor in header\n";
  return 0;
}

static int
write_record (const_IWSubstring & buffer,
              IWString_and_File_Descriptor & output)
{
  if (! write_passing_records)
    return 1;

  if (just_write_identifiers)
    buffer.truncate_at_first(' ');

  output << buffer << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return output.good();
}

static int
dfilefilter (iwstring_data_source & input,
             int ncol,
             Column_Data * cd,
             IWString_and_File_Descriptor & output)
{
  input.set_dos(1);

  const_IWSubstring buffer;
  while (input.next_record(buffer) && output.good())
  {
    records_read++;

//  cerr << "Just read '" << buffer << "'\n";

    if (! tokenise(buffer, ncol, cd))
    {
      if (missing_column_ok) {
        if (stream_for_rejected_records.is_open())
          write_record(buffer, stream_for_rejected_records);
        continue;
      }
      else {
        cerr << "Cannot tokenise line " << input.lines_read() << endl;
        cerr << buffer << endl;
        return 0;
      }
    }

    if (! passes_all_filters(cd, ncol))
    {
      if (stream_for_rejected_records.is_open())
        write_record(buffer, stream_for_rejected_records);
    }
    else
    {
      write_record(buffer, output);

      records_written++;
    }
  }

  return output.good();
}

static int 
dfilefilter (iwstring_data_source & input, 
             IWString_and_File_Descriptor & output)
{
  const_IWSubstring header;

  if (! input.next_record(header))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

// identify those columns for which numeric values are needed

  int ncol;
  if (' ' == input_separator)
    ncol = header.nwords();
  else
    ncol = header.nwords_single_delimiter(input_separator);

  if (verbose)
    cerr << "Input file contains " << ncol << " columns\n";

  highest_column_needed = 0;

  for (int i = 0; i < nc; i++)
  {
    Column_Condition & c = column_condition[i];

    if (! c.determine_column_number(header))
    {
      return 0;
    }

    int col = c.column();

    if (verbose)
      cerr << "Column " << (col + 1) << " is descriptor '" << c.column_name() << "'\n";

    assert (col >= 0 && col < ncol);

    if (col > highest_column_needed)
      highest_column_needed = col;
  }

  if (0 == first_non_reagent_column_header.length())
    ;
  else if (! identify_first_non_reagent_column(header))
  {
    cerr << "Cannot identify descriptor for trxn processing\n";
    return 0;
  }
  else if (number_columns_reagent_types > highest_column_needed)
    highest_column_needed = number_columns_reagent_types;

  if (all_other_columns_zero)
    initialise_check_column_for_zero(ncol);

  if (just_write_identifiers)     // ignore header
    ;
  else if (header_record_written)
    ;
  else if (! write_passing_records)
    ;
  else
  {
    output << header << '\n';
    header_record_written = 1;
  }

  Column_Data * cd = new Column_Data[ncol];

  int rc = dfilefilter(input, ncol, cd, output);

  delete [] cd;

  return rc;
}

/*
  When parsing the expression, we need these classes
*/

class Operator;

class Condition
{
  private:
    IWString _text;

    Operator * _op;

  public:
    Condition ();
    ~Condition ();

    int debug_print (std::ostream &) const;

    int build (const const_IWSubstring &);
    int build (const IWString &, int);

    const IWString & text () const { return _text;}

    int operator_count () const;

    Condition * next_condition () const;

    Operator * op () const { return _op;}
};

#define OPERATOR_TYPE_AND 1
#define OPERATOR_TYPE_OR 2
#define OPERATOR_TYPE_XOR 3
#define OPERATOR_TYPE_LOW_PRIORITY_AND 4

class Operator
{
  private:
    int _type;

    Condition _next;

  public:
    Operator ();
    ~Operator ();

    int debug_print (std::ostream &) const;

    int build (const IWString &, int);

    int optype () const { return _type;}

    int operator_count () const;

    Condition * next_condition () { return & _next;}
};

Condition::Condition ()
{
  _op = NULL;

  return;
}

Condition::~Condition ()
{
  if (NULL != _op)
    delete _op;

  return;
}

int
Condition::debug_print (std::ostream & os) const
{
  os << "Condition '" << _text << "'\n";

  if (NULL == _op)
    return os.good();

  return _op->debug_print(os);
}

int
Condition::build (const const_IWSubstring & e)
{
  IWString expr(e);

  expr.remove_all(' ');    // get rid of blanks

  if (0 == expr.length())
  {
    cerr << "Condition::build: empty expression '" << e << "'\n";
    return 0;
  }

  return build(expr, 0);
}

/*
  Our part of the expression will be from where we start until we encounter
  a double operator '&&', etc..
*/

int
Condition::build (const IWString & expr, int istart)
{
  assert (NULL == _op);

  int iptr = istart;

  while (iptr < expr.length())
  {
    char c = expr[iptr];

    if ('|' == c || '&' == c || '^' == c || ';' == c)
    {
      if (iptr == expr.length() - 1)
      {
        cerr << "Condition::build: expressions cannot end in operators\n";    // well, they probably could...
        return 0;
      }

      if (c == expr[iptr + 1])
        break;
    }

    iptr++;
  }

  iptr--;
  if (iptr == istart)
  {
    cerr << "Condition::build: empty expression. Impossible\n";
    return 0;
  }

  expr.from_to(istart, iptr, _text);

  if (iptr >= expr.length() - 1)    // we are done
    return 1;

  _op = new Operator;

  return _op->build(expr, iptr + 1);
}

int
Condition::operator_count () const
{
  if (NULL == _op)
    return 0;

  return _op->operator_count();
}

Condition *
Condition::next_condition () const
{
  if (NULL == _op)
    return NULL;

  return _op->next_condition();
}

Operator::Operator ()
{
  _type = -1;

  return;
}

Operator::~Operator ()
{
  return;
}

int
Operator::debug_print (std::ostream & os) const
{
  if (_type < 0)
  {
    os << "Operator::debug_print: operator not set\n";
  }
  else if (OPERATOR_TYPE_AND == _type)
    os << "Operator AND";
  else if (OPERATOR_TYPE_OR == _type)
    os << "Operator OR";
  else if (OPERATOR_TYPE_XOR == _type)
    os << "Operator XOR";
  else if (OPERATOR_TYPE_LOW_PRIORITY_AND == _type)
    os << "Operator LOW_PRIORITY_AND";

  os << ' ';

  return _next.debug_print(os);
}

int
Operator::build (const IWString & expr, int iptr)
{
  if (iptr + 2 >= expr.length())
  {
    cerr << "Operator::build: expression too short '" << expr << "' iptr = " << iptr << endl;
    return 0;
  }

  char firstchar = expr[iptr];

  if (firstchar != expr[iptr + 1])
  {
    cerr << "Operator::build: operators must be repeated '" << expr << "' iptr = " << iptr << endl;
    return 0;
  }

  if ('&' == firstchar)
    _type = OPERATOR_TYPE_AND;
  else if ('|' == firstchar)
    _type = OPERATOR_TYPE_OR;
  else if ('^' == firstchar)
    _type = OPERATOR_TYPE_XOR;
  else if (';' == firstchar)
    _type = OPERATOR_TYPE_LOW_PRIORITY_AND;
  else
  {
    cerr << "Operator::build: unrecognised operator '" << expr << "'\n";
    return 0;
  }

  return _next.build(expr, iptr + 2);
}

int
Operator::operator_count () const
{
  return 1 + _next.operator_count();
}

static int
build_logical_expression (Condition * c)
{
  if (1 == nc)    // just one expression (foo > 3), no logical expression
    return 1;

  for (int i = 0; i < nc - 1; i++)
  {
    const Operator * op = c->op();

    int t = op->optype();

    if (OPERATOR_TYPE_AND == t)
      logical_expression.add_operator(IW_LOGEXP_AND);
    else if (OPERATOR_TYPE_OR == t)
      logical_expression.add_operator(IW_LOGEXP_OR);
    else if (OPERATOR_TYPE_XOR == t)
      logical_expression.add_operator(IW_LOGEXP_XOR);
    else if (OPERATOR_TYPE_LOW_PRIORITY_AND == t)
      logical_expression.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
    else
    {
      cerr << "Huh, what kind of operator is this " << t << endl;
      return 0;
    }

    c = c->next_condition();
  }

  if (verbose)
  logical_expression.debug_print(cerr);

  return 1;
}

static int
build_column_conditions (Condition * c)
{
  for (int i = 0; i < nc; i++)
  {
    Column_Condition & cc = column_condition[i];

    const IWString & s = c->text();

    if (verbose)
      cerr << "Building column condition " << i << " from '" << s << "'\n";

    if (! cc.build(s))
    {
      cerr << "Cannot parse expression '" << c->text() << "'\n";
      return 0;
    }

    c = c->next_condition();
    if (NULL == c)
      return 1;
  }

  return 1;
}

static int
dfilefilter (const char * fname, 
             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return dfilefilter(input, output);
}

static int
parse_expression (const const_IWSubstring & e,
                  Condition & c)
{
  if (! c.build(e))
  {
    cerr << "Cannot parse expression '" << e << "'\n";
    return 0;
  }

  return 1;
}

static int
dfilefilter (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "ve:m:pciB:K:zs:n");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('m'))
  {
    cl.value('m', missing_value);
    if (verbose)
      cerr << "The missing value string is '" << missing_value << "'\n";
  }

  if (cl.option_present('p'))
  {
    missing_value_passes = 1;
    if (verbose)
      cerr << "Missing values pass all filters\n";
  }

  if (cl.option_present('c'))
  {
    missing_column_ok = 1;
    if (verbose)
      cerr << "Missing column is non-terminal failure\n";
  }

  if (cl.option_present('i'))
  {
    just_write_identifiers = 1;

    if (verbose)
      cerr << "Will just write identifiers (column 1)\n";
  }

  if (cl.option_present('s'))
  {
    IWString s = cl.string_value('s');
    if (! char_name_to_char(s))
    {
      cerr << "Unrecognised input separator (-s ) '" << s << "'\n";
      usage(1);
    }

    input_separator = s[0];
  }

  if (cl.option_present('n'))
  {
    write_passing_records = 0;

    if (verbose)
      cerr << "Normal output suppressed\n";
  }

  if (1 != cl.option_count('e'))
  {
    cerr << "Must specify one (and only one) -e option with the selection expression\n";
    usage(2);
  }

  if (cl.option_present('e'))
  {
    const_IWSubstring e;
    cl.value('e', e);

    Condition c;
    if (! parse_expression(e, c))
    {
      cerr << "Cannot process -e option\n";
      usage(15);
    }

    nc = c.operator_count() + 1;

    if (verbose)
      c.debug_print(cerr);

    column_condition = new Column_Condition[nc];

    if (verbose)
      cerr << "Expression contains " << nc << " conditions\n";

    if (! build_column_conditions(&c))
    {
      cerr << "Bad -e specifier '" << e << "'\n";
      usage(18);
    }

    if (! build_logical_expression(&c))
    {
      cerr << "Cannot build logical expression\n";
      usage(7);
    }
  }

  if (cl.option_present('K'))
  {
    first_non_reagent_column_header = cl.string_value('K');

    if (verbose)
      cerr << "Special processing for trxn, first non reagent column '" << first_non_reagent_column_header << "'\n";

  }

  if (cl.option_present('z'))
  {
    all_other_columns_zero = 1;

    if (verbose)
      cerr << "Selected records written only if all other columns zero\n";
  }

  if (cl.option_present('B'))
  {
    IWString b = cl.string_value('B');

    if (! stream_for_rejected_records.open(b.null_terminated_chars()))
    {
      cerr << "Cannot open -B file '" << b << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Rejected records written to '" << b << "'\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(3);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! dfilefilter(cl[i], output))
    {
      cerr << "Error processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (stream_for_rejected_records.is_open())
    stream_for_rejected_records.flush();

  if (verbose || ! write_passing_records)
  {
    cerr << "Read " << records_read << " records, wrote " << records_written;
    if (records_read > 0)
      cerr << ", fraction " << static_cast<float>(records_written) / static_cast<float>(records_read) << "\n";
    else
      cerr << endl;
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = dfilefilter(argc, argv);

  return rc;
}
