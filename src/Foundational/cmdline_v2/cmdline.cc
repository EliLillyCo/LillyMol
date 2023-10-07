#include <stdlib.h>


#include "cmdline_v2.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;

/*
  All the different qualifiers that can follow an option

  We use a special flag to indicate that the option is
  interpretable as a numeric
*/

#define IWCMDLINE_NUMERIC 4096

#define IWCMDLINE_NONE 0
#define IWCMDLINE_STRING 1
#define IWCMDLINE_INTEGER (2 | IWCMDLINE_NUMERIC)
#define IWCMDLINE_POSITIVE_INTEGER (3 | IWCMDLINE_NUMERIC)
#define IWCMDLINE_UNSIGNED_INTEGER (4 | IWCMDLINE_NUMERIC)
#define IWCMDLINE_FLOAT (5 | IWCMDLINE_NUMERIC)
#define IWCMDLINE_FRACTION (6 | IWCMDLINE_NUMERIC)

#define IWCMDLINE_XFILE 7
#define IWCMDLINE_SFILE 8
#define IWCMDLINE_DIRECTORY 9

#define IWCMDLINE_CLOSE 10

CmdLineDefault::CmdLineDefault ()
{
  return;
}

CmdLineDefault cmdlinedefault;

CmdLine_Option_and_Value::CmdLine_Option_and_Value ()
{
  _opt = nullptr;
  _valid_as = IWCMDLINE_NONE;

  return;
}

int
CmdLine_Option_and_Value::is_option (const char * o,
                                     int leno) const
{
  if (nullptr == _opt)   // happens with erroneous initialisations
    return 0;

  if (static_cast<size_t>(leno) != strlen(_opt))
    return 0;

  return 0 == strncmp(_opt, o, leno);
}

int
CmdLine_Option_and_Value::initialise (const char * opt)
{
  _opt = opt;

  return 1;
}

int
CmdLine_Option_and_Value::initialise (const char * opt,
                                      int valid_as,
                                      const IWString & v)
{
  _opt = opt;
  _value = v;
  _value.null_terminate();

  if (IWCMDLINE_STRING == valid_as)
    ;
  else if (IWCMDLINE_INTEGER == valid_as)
  {
    long int tmp;
    if (! _value.numeric_value(tmp))
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' invalid integer value '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_POSITIVE_INTEGER == valid_as)
  {
    long int tmp;
    if (! _value.numeric_value(tmp) || tmp < 0)    // should this be <= ?
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' invalid positive integer value '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_UNSIGNED_INTEGER == valid_as)
  {
    unsigned long tmp;
    if (! _value.numeric_value(tmp) || tmp < 0)
    {
      cerr << "CmdLine_Option_and_Value:set_value: '-" << _opt << "' invalid unsigned integer value '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_FLOAT == valid_as)
  {
    double tmp;
    if (! _value.numeric_value(tmp))
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' invalid float value '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_FRACTION == valid_as)
  {
    double tmp;
    if (! _value.numeric_value(tmp) || tmp < 0.0 || tmp > 1.0)
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' invalid fraction value '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_SFILE == valid_as)
  {
    if (! dash_s(_value.null_terminated_chars()))
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' missing or empty file '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_XFILE == valid_as)
  {
    if (! dash_x(_value.null_terminated_chars()))
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' missing inaccessible executable  '" << _value << "'\n";
      return 0;
    }
  }
  else if (IWCMDLINE_DIRECTORY == valid_as)
  {
    if (! dash_d(_value.null_terminated_chars()))
    {
      cerr << "CmdLine_Option_and_Value:set_value: '" << _opt << "' not a valid directory  '" << _value << "'\n";
      return 0;
    }
  }
  else
  {
    cerr << "CmdLine_Option_and_Value::set_value: unrecognised type " << valid_as << endl;
    return 0;
  }

  return 1;
}

template <typename T>
int
CmdLine_Option_and_Value::value(T & v) const
{
//cerr << "CmdLine_Option_and_Value:converting '" << _value << "' to type " << sizeof(T) << endl;
  return _value.numeric_value(v);
}

template int CmdLine_Option_and_Value::value(int &) const;
template int CmdLine_Option_and_Value::value(unsigned int &) const;
template int CmdLine_Option_and_Value::value(long &) const;
template int CmdLine_Option_and_Value::value(unsigned long &) const;
template int CmdLine_Option_and_Value::value(float &) const;
template int CmdLine_Option_and_Value::value(double &) const;
template int CmdLine_Option_and_Value::value(long long &) const;

/*
  For each of the options, we need a structure to hold the
  string for the option, and its type.
  We could do this with virtual functions, but I've had
  problems with the SGI compiler doing resizable_array's
  of objects with virtual functions - sigh...
*/

class Option_Type
{
  private:
    IWString _s;   // needs permanent storage
    int _type;

  public:
    Option_Type ();

    const IWString & flag () const { return _s;}

    int construct (const const_IWSubstring &);

    int takes_qualifier () const { return (IWCMDLINE_NONE != _type);}

    int valid_as () const { return _type;}
};

Option_Type::Option_Type ()
{
  _type = IWCMDLINE_NONE;

  return;
};

int
Option_Type::construct (const const_IWSubstring & s)
{
  if (! s.contains('='))
  {
    _s = s;
    return 1;
  }

  const_IWSubstring opt, d;

  if (! s.split(opt, '=', d))
  {
    cerr << "Option_Type::construct:invalid specification '" << s << "'\n";
    return 0;
  }

  _s = opt;

  if ("s" == d)
    _type = IWCMDLINE_STRING;
  else if ("i" == d || "int" == d)
    _type = IWCMDLINE_INTEGER;
  else if ("ipos" == d)
    _type = IWCMDLINE_POSITIVE_INTEGER;
  else if ("u" == d || "uint" == d)
    _type = IWCMDLINE_UNSIGNED_INTEGER;
  else if ("f" == d || "float" == d)
    _type = IWCMDLINE_FLOAT;
  else if ("fraction" == d)
    _type = IWCMDLINE_FRACTION;
  else if ("xfile" == d)
    _type = IWCMDLINE_XFILE;
  else if ("sfile" == d)
    _type = IWCMDLINE_SFILE;
  else if ("dir" == d)
    _type = IWCMDLINE_DIRECTORY;
  else if ("close" == d)
    _type = IWCMDLINE_CLOSE;
  else
  {
    cerr << "Option_Type::construct:unrecognised directive '" << s << "'\n";

    return 0;
  }

  return 1;
}

Command_Line_v2::Command_Line_v2 (int argc,
                  char * const * argv,
                  const char * schar,
                  const CmdLineDefault & c)
{
  _number_option_value_pairs = 0;
  _clov = nullptr;

  _unrecognised_options_encountered = 0;

#ifdef DEBUG_COMMAND_LINE
  cerr << "Command_Line::on entry\n";
  for (int i = 0; i < argc; i++)
  {
    cerr << "'" << argv[i] << "'\n";
  }
#endif

  resize_keep_storage(0);

  const const_IWSubstring s(schar);

  if (0 == s.length())
  {
    _load_rest_into_array(argc, argv, 0);
    return;
  }

  int number_options = s.ccount('-');

// Cannot have something like "v-foo", all options must be the same form

  if (s.starts_with('-'))
    ;
  else if (number_options)
  {
    cerr << "Command_Line:cannot mix old and new forms '" << s << "'\n";
    _good = 0;
    return;
  }
  else
  {
    _initialise_old_form(argc, argv, s);
    return;
  }

  if (0 == number_options)
  {
#ifdef DEBUG_COMMAND_LINE
    cerr << "Command_Line::no options specified\n";
#endif

    _load_rest_into_array(argc, argv, 0);

    return;
  }

  Option_Type * ot = new Option_Type[number_options];

  _constructor(argc, argv, s, ot, number_options, c);

  delete [] ot;

  return;
}

Command_Line_v2::~Command_Line_v2 ()
{
  if (nullptr != _clov)
    delete [] _clov;

  if (-5 == _number_option_value_pairs)
    cerr << "Command_Line:deleting already deleted Command_Line\n";

  _number_option_value_pairs = -5;

  return;
}

static const Option_Type *
find_corresponding_option_and_type (const Option_Type * ot,
                                    int n,
                                    const char * opt)
{
  if ('-' == *opt)
    opt++;

//cerr << "Testing " << n << " options for match to '" << opt << "'\n";

  for (int i = 0; i < n; i++)
  {
//  cerr << "Compare '" << ot[i].flag() << "' and '" << opt << "'\n";

    if (ot[i].flag() == opt)
      return &(ot[i]);
  }

  return nullptr;
}

int
Command_Line_v2::_load_rest_into_array (int argc,
                                char * const * argv,
                                int istart)
{
//cerr << "Loading from " << istart << " to " << argc << endl;
  for (int i = istart; i < argc; i++)
  {
    resizable_array<const char *>::add(argv[i]);
  }

  return _number_elements;
}

int
Command_Line_v2::_constructor (const int argc,
                       char * const * argv,
                       const const_IWSubstring & s,
                       Option_Type * ot,
                       const int number_options,
                       const CmdLineDefault & c)
{
  _good = 1;

  int i = 0;
  const_IWSubstring option_and_directive;
  int nopt = 0;

  while (s.nextword(option_and_directive, i, '-'))
  {
    if (! ot[nopt].construct(option_and_directive))
    {
      cerr << "Command_Line::_constructor:invalid specification '" << option_and_directive << "'\n";
      _good = 0;
      return 0;
    }

    nopt++;
  }

// Determine an upper bound on the number of - options 

  _number_option_value_pairs = 0;
  for (int i = 0; i < argc; i++)
  {
    const char * s = argv[i];

    if (0 == strlen(s))   // huh?
      ;
    else if ('-' == *s)
      _number_option_value_pairs++;
  }

  if (_number_option_value_pairs)
    _clov = new CmdLine_Option_and_Value[_number_option_value_pairs];

  _number_option_value_pairs = 0;

  assert (number_options == nopt);

  int ndx;
  for (ndx = 1; ndx < argc; ndx++)
  {
    const char * opt = argv[ndx];

#ifdef DEBUG_CMDLINE
    cerr << "ndx " << ndx << " considering '" << opt << "'\n";
#endif

    int len_opt = static_cast<int>(strlen(opt));

    if (0 == len_opt)   // not sure what this means
    {
      cerr << "CmdLine:warning, zero length item\n";
      break;
    }

    if ('-' != *opt)
      break;

    if (1 == len_opt)
      break;

    const Option_Type * o = find_corresponding_option_and_type(ot, number_options, opt);

    if (nullptr == o)
    {
      _good = 0;
      _unrecognised_options_encountered++;
      continue;
    }

    CmdLine_Option_and_Value & clov = _clov[_number_option_value_pairs];
    _number_option_value_pairs++;

#ifdef DEBUG_CMDLINE
    if (o->takes_qualifier())
      cerr << "Takes qualifier\n";
    else
      cerr << "no qualifier\n";
#endif

    if (! o->takes_qualifier())
      clov.initialise(opt + 1);
    else if (IWCMDLINE_CLOSE == o->valid_as())
    {
//    cerr << "Option '" << opt << "' takes closing pair\n";
      int got_closing_option = 0;

      IWString v;
      while (ndx < argc)
      {
        ndx++;

        const char * tmp = argv[ndx];
        if (0 == strcmp(opt, tmp))
        {
          got_closing_option = 1;
          break;
        }

        v.append_with_spacer(tmp);
      }

      if (! got_closing_option)
      {
        cerr << "Command_Line::unclosed -" << opt << " option\n";
        _good = 0;
        return 0;
      }

//    cerr << "Qualifier for '" << opt << "' is '" << v << "'\n";

      v.null_terminate();

      clov.initialise(opt, IWCMDLINE_STRING, v);
    }
    else
    {
      ndx++;
      if (argc == ndx)
      {
        cerr << "Command_Line::_constructor:the '" << o->flag() << "' option requires a qualifier\n";
        _good = 0;
        return 0;
      }

      const char * qualifier = argv[ndx];

      if (! clov.initialise(opt + 1, o->valid_as(), qualifier))   // opt+1 to skip over the leading -
      {
        cerr << "Command_Line:invalid option value\n";
        _good = 0;
      }
    }
  }

  if (ndx < argc)
      return _load_rest_into_array(argc, argv, ndx);

  return 1;
}

int
Command_Line_v2::option_present(char c) const
{
  return option_present(&c, 1);
}

int
Command_Line_v2::option_present(const char * s, int lens) const
{
#ifdef DEBUG_OPTION_PRESENT
  cerr << "Search " << _number_option_value_pairs << " option value pairs for '";
  cerr.write(s, lens);
  cerr << "'\n";
#endif

  for (int i = 0; i < _number_option_value_pairs; i++)
  {
    const CmdLine_Option_and_Value & c = _clov[i];

    if (c.is_option(s, lens))
      return 1;
  }

#ifdef DEBUG_OPTION_PRESENT
  cerr << "Option not present\n";
#endif

  return 0;
}

int
Command_Line_v2::option_count (const char * opt, int lens) const
{
  int rc = 0;

  for (int i = 0; i < _number_option_value_pairs; i++)
  {
    const CmdLine_Option_and_Value & c = _clov[i];

    if (c.is_option(opt, lens))
      rc++;
  }

  return rc;
}

template <typename T>
int
Command_Line_v2::value (int len_opt,
                     const char * opt,
                     T & v,
                     int which_one) const
{
  int number_found = 0;

  for (int i = 0; i < _number_option_value_pairs; i++)
  {
    const CmdLine_Option_and_Value & c = _clov[i];

    if (! c.is_option(opt, len_opt))
      continue;

    if (which_one == number_found)
      return c.value(v);

    number_found++;
  }

  return 0;
}

const char *
Command_Line_v2::option_value (int len_opt,
                            const char * opt,
                            int which_one) const
{
  int number_found = 0;

  for (int i = 0; i < _number_option_value_pairs; i++)
  {
    const CmdLine_Option_and_Value & c = _clov[i];

    if (! c.is_option(opt, len_opt))
      continue;

    if (which_one == number_found)
      return c.value().rawchars();    // may be NULL if never set. Hmmm, should we distinguish NULL from the zero length string. Hmmm, should we distinguish NULL from the zero length string? Actually, it would have been NULL terminated, so maybe oK

    number_found++;
  }

  return nullptr;
}

const const_IWSubstring 
Command_Line_v2::string_value (int len_opt,
                            const char * opt, 
                            int which_one) const
{
  int number_found = 0;

  for (int i = 0; i < _number_option_value_pairs; i++)
  {
    const CmdLine_Option_and_Value & c = _clov[i];

    if (! c.is_option(opt, len_opt))
      continue;

    if (which_one == number_found)
    {
      const IWString & v = c.value();
      if (0 == v.length())
        return "";

      return v;
    }

    number_found++;
  }

  return "";

}
template int Command_Line_v2::value(int, const char *, IWString &, int) const;
template int Command_Line_v2::value(int, const char *, const_IWSubstring &, int) const;
template int Command_Line_v2::value(int, const char *, int &, int) const;
template int Command_Line_v2::value(int, const char *, unsigned int &, int) const;
template int Command_Line_v2::value(int, const char *, long &, int) const;
template int Command_Line_v2::value(int, const char *, float &, int) const;
template int Command_Line_v2::value(int, const char *, double &, int) const;
template int Command_Line_v2::value(int, const char *, long long &, int) const;
template int Command_Line_v2::value(int, const char *, unsigned long &, int) const;

#ifdef IW_STD_STRING_DEFINED
template int Command_Line_v2::value(int, const char *, std::string &, int) const;
#endif

/*template <typename T>
int
Command_Line_v2::value (char o, T & s, int which_one) const
{
  return value (1, &o, s, which_one);
}*/

template int Command_Line_v2::value(char, IWString &, int) const;
template int Command_Line_v2::value(char, const_IWSubstring &, int) const;
template int Command_Line_v2::value(char, int &, int) const;
template int Command_Line_v2::value(char, unsigned int &, int) const;
template int Command_Line_v2::value(char, float &, int) const;
template int Command_Line_v2::value(char, double &, int) const;

int
Command_Line_v2::_initialise_old_form(int argc,
                                    char * const * argv,
                                    const const_IWSubstring & s)
{
// Work out how many options there will be

  int colons = s.ccount(':');

  int n = s.length() - colons;

  Option_Type * ot = new Option_Type[n];

  int rc = _initialise_old_form(argc, argv, s, ot);

  delete [] ot;

  return rc;
}

int
Command_Line_v2::_initialise_old_form (int argc,
                                    char * const * argv,
                                    const const_IWSubstring & s,
                                    Option_Type * ot)
{
  _good = 1;

  IWString c;

  int nopt = 0;

  for (int i = 0; i < s.length(); i++)
  {
    c = s[i];

    if (i < s.length() - 1 && ':' == s[i + 1])
    {
      i++;
      c << "=s";
    }

//  cerr << "Building option from '" << c << "'\n";
    (void) ot[nopt].construct(c);
//  cerr << "Built flag '" << ot[nopt].flag() << "'\n";

    nopt++;
  }

// Determine an upper bound on the number of - options 

  _number_option_value_pairs = 0;
  for (int i = 0; i < argc; i++)
  {
    const char * s = argv[i];

    if (0 == strlen(s))   // huh?
      ;
    else if ('-' == *s)
      _number_option_value_pairs++;
  }

  if (_number_option_value_pairs)
    _clov = new CmdLine_Option_and_Value[_number_option_value_pairs];

  _number_option_value_pairs = 0;

  int ndx;
  for (ndx = 1; ndx < argc; ndx++)
  {
    const char * opt = argv[ndx];

//#define DEBUG_CMDLINE
#ifdef DEBUG_CMDLINE
    cerr << "ndx " << ndx << " considering '" << opt << "'\n";
#endif

    int len_opt = static_cast<int>(strlen(opt));

    if (0 == len_opt)   // not sure what this means
    {
      cerr << "CmdLine:warning, zero length item\n";
      break;
    }

    if ('-' != *opt)
      break;

    if (1 == len_opt)
      break;

    const Option_Type * o = find_corresponding_option_and_type(ot, nopt, opt);

    if (nullptr == o)
    {
    cerr << "No option for '" << opt << "'\n";
      _good = 0;
      _unrecognised_options_encountered++;
      continue;
    }

    CmdLine_Option_and_Value & clov = _clov[_number_option_value_pairs];
    _number_option_value_pairs++;

#ifdef DEBUG_CMDLINE
    if (o->takes_qualifier())
      cerr << "Takes qualifier\n";
    else
      cerr << "no qualifier\n";
#endif

    if (! o->takes_qualifier())
      clov.initialise(opt + 1);
    else
    {
      ndx++;
      if (argc == ndx)
      {
        cerr << "Command_Line::_constructor:the '" << o->flag() << "' option requires a qualifier\n";
        _good = 0;
        return 0;
      }

      const char * qualifier = argv[ndx];

      if (! clov.initialise(opt + 1, o->valid_as(), qualifier))   // opt+1 to skip over the leading -
      {
        cerr << "Command_Line:invalid option value\n";
        _good = 0;
      }
    }
  }

  if (ndx < argc)
      return _load_rest_into_array(argc, argv, ndx);

  return 1;
}
