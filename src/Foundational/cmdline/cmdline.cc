//#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <limits>
#include "iwconfig.h"

#undef _GETOPT_H

// Unclear when we need to switch between C and C++ bindings for optxxx
#if defined(_WIN32) || defined(NEED_EXTERN_OPT)
/* Global Exportable */
 extern "C" int optind;
 extern "C" char *optarg;
 extern "C" int opterr;

#else
#include <unistd.h>
#endif

using std::cerr;
using std::endl;

#include "cmdline.h"
#include "iwstring.h"

#define CL_MAGIC 97531
#define OV_MAGIC 97532

#define UNKNOWN_VALIDITY -1

#define VALID_AS_INT 1
#define VALID_AS_UNSIGNED_INT 2

Option_and_Value::Option_and_Value (int o, const char * val) : _o(o), _value(val)
{
  _magic = OV_MAGIC;

  return;
}

int
Option_and_Value::ok () const
{
  if (OV_MAGIC != _magic)
    return 0;

  return 1;
}

Option_and_Value::~Option_and_Value ()
{
  assert (ok());

  _magic = -5;
}

int
Option_and_Value::value (IWString & result)
{
  if (NULL == _value)
    return 0;

  result = _value;

  return 1;
}

int
Option_and_Value::value (const_IWSubstring & result)
{
  if (NULL == _value)
    return 0;

  result = _value;

  return 1;
}

int
Option_and_Value::value (char * buffer)
{
  if (NULL == _value)
    return 0;

  IW_STRCPY(buffer, _value);

  return 1;
}

std::ostream &
operator << (std::ostream & os, const Option_and_Value & ov)
{
  return os << "Option '" << ov.option() << "', value '" << ov.value() << "'";
}

Command_Line::Command_Line (int argc, char ** argv, const char * options)
{
  _magic = CL_MAGIC;

#if defined(_WIN32)
// unsure about this, definitely needed in some circumstances, but behaviour
// keeps changing with compiler version and platform
  argptr = NULL;

#endif

  optarg = NULL;

#if defined(_WIN32) || defined(NEED_EXTERN_OPT)
  optind = 0;
#else
  optind = 1;   // reinitialise in case of multiple invocations
  opterr = 0;     // suppress error messages
#endif

  _options.resize(argc);     // yes, this is too many

  _some_options_start_with_dash = 0;
  _unrecognised_options_encountered = 0;

  int o;


  while ((o = getopt(argc, argv, options)) != EOF)
  {
    if ('?' == o)
    {
      cerr << "Command_Line: unrecognised option '" << argv[optind - 1] << "'\n";
      _unrecognised_options_encountered++;
    }
    else
    {
      Option_and_Value * tmp = new Option_and_Value(o, optarg);
      _options.add(tmp);
    }
  }

  resize (argc - optind);

  for (int i = optind; i < argc; i++)
  {
    if ('-' == *(argv[i]))
      _some_options_start_with_dash++;
    add(argv[i]);
  }

  return;
}

Command_Line::~Command_Line ()
{
  assert (ok());

  _magic = -9;
}

int
Command_Line::ok () const
{
  if (_magic != CL_MAGIC)
    return 0;

  for (int i = 0; i < _options.number_elements(); i++)
    if (! _options[i]->ok())
      return 0;

  return resizable_array<const char *>::ok();
}

int
Command_Line::debug_print (std::ostream & os) const
{
  assert (os.good());

  os << "Command line object contains " << _options.number_elements() << " options and " <<
        _number_elements << " values\n";

  for (int i = 0; i < _options.number_elements(); i++)
    os << *(_options[i]) << endl;

  for (int i = 0; i < _number_elements; i++)
    os << "Value '" << _things[i] << "'\n";

  return 1;
}

int
Command_Line::option_present (const char c) const
{
  for (int i = 0; i < _options.number_elements(); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option())
      return i + 1;
  }

  return 0;
}

const char *
Command_Line::option_value (const char c, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option() && occurrence == nfound++)
      return oo->value();
  }

  return 0;
}

int
Command_Line::all_values (const char c, resizable_array<const char *> & values) const
{
  int rc = 0;

  for (int i = 0; i < _options.number_elements(); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      values.add (oo->value());
      rc++;
    }
  }
  return rc;
}

/*
  When returning all instances of an option, we can optionally
  further tokenise these instances. 

  prog -r a -r 'c d'

  would return 3 elements (if split_tokens is specified)
*/

int
Command_Line::all_values (const char c, resizable_array_p<IWString> & values,
                          int split_tokens) const
{
  int rc = 0;

  for (int i = 0; i < _options.number_elements(); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c != oo->option())
      continue;


    if (0 == split_tokens || 0 == ::ccount(oo->value(), ' '))    // no need to split, either just one token, or not requested
    {
      rc++;
      IWString * tmp = new IWString(oo->value());
      values.add(tmp);
      continue;
    }

    const_IWSubstring v(oo->value());
    int j = 0;
    const_IWSubstring token;
    while (v.nextword(token, j))
    {
      IWString * tmp = new IWString(token);
      values.add(tmp);
      rc++;
    }
  }

  return rc;
}

int
Command_Line::option_count (const char c) const
{
  int rc = 0;

  for (int i = 0; i < _options.number_elements(); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option())
      rc++;
  }

  return rc;
}

template <typename T>
int
Command_Line::value (const char c, T & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->value(result);

      nfound++;
    }
  }

  return 0;
}

template int Command_Line::value(const char, int &, int) const;
template int Command_Line::value(const char, unsigned int &, int) const;
template int Command_Line::value(const char, long int &, int) const;
template int Command_Line::value(const char, unsigned long int &, int) const;
template int Command_Line::value(const char, long long int &, int) const;
template int Command_Line::value(const char, unsigned long long int &, int) const;
template int Command_Line::value(const char, double &, int) const;
template int Command_Line::value(const char, float &, int) const;


/*int
Command_Line::double_val (const char c, double & value, int occurrence)
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->double_val(value);
      nfound++;
    }
  }

  return 0;
}*/

#ifdef NOT_USED_QWREKJQHWELK
int
Command_Line::value (const char c, double & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->value(result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, float & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->value(result);
      nfound++;
    }
  }

  return 0;
}
#endif

int
Command_Line::value (const char c, IWString & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->value(result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, const_IWSubstring & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->value(result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, char * result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
        return oo->value(result);
      nfound++;
    }
  }

  return 0;
}

const const_IWSubstring
Command_Line::string_value (const char c, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements(); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option())
    {
      if (nfound == occurrence)
      {
        const char * v = oo->value();
        if (NULL == v)
          return "";
        else
          return oo->value();
      }

      nfound++;
    }
  }

  return "";
}

std::ostream &
operator << (const Command_Line & cl, std::ostream & os)
{
  return os;
}
