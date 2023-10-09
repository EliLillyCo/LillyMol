#include <stdlib.h>
#include <assert.h>
//#include <ctype.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "molecule.h"
#include "numass.h"

using std::cerr;
using std::endl;

Number_Assigner::Number_Assigner()
{
  _prefix_string = 'R';
  _next_number_to_assign = -1;
  _replace_name = 0;

  _include_parentheses = 1;

  _number_digits = 0;

  _only_apply_if_existing_name_is_empty = 0;

  _separator_to_existing = ' ';

  return;
}

int
Number_Assigner::initialise(Command_Line & cl,
                            char nflag,
                            int verbose)
{
  if (! cl.option_present(nflag))
    return 1;

  int i = 0;
  const_IWSubstring n;

  int rc = 0;

  while (cl.value(nflag, n, i++))
  {
    int tmp;
    if (n.numeric_value(tmp))
    {
      if (tmp < 0)
      {
        cerr << "Number_Assigner::initialise: the '" << nflag << "' option must be followed by a positive whole number\n";
        return 0;
      }

      _next_number_to_assign = tmp;
      rc++;
      continue;
    }

    if ("replace" == n)
    {
      _replace_name = 1;
      if (verbose)
        cerr << "Will replace name with sequence number\n";

      continue;
    }

    if ("rpfirst" == n)
    {
      _replace_first_token = 1;
      if (verbose)
        cerr << "Will replace only the first token of an existing name\n";
      continue;
    }

    if ("noparen" == n)
    {
      _include_parentheses = 0;

      if (verbose)
        cerr << "No parentheses used\n";

      continue;
    }

    if (n.starts_with("ndigits="))
    {
      n.remove_leading_chars(8);

      if (! n.numeric_value(_number_digits) || _number_digits < 1)
      {
         cerr << "Invalid number of digits specification '" << n << "'\n";
         return 0;
      }

      if (verbose)
        cerr << "Numeric component padded to " << _number_digits << " digits with leading 0's\n";

      continue;
    }

    if ("nochange" == n)
    {
      _only_apply_if_existing_name_is_empty = 1;

      if (verbose)
        cerr << "Will not change any existing name\n";

      continue;
    }

    if (n.starts_with("nsep="))
    {
      n.remove_leading_chars(5);
      _separator_to_existing = n;

      if (verbose)
        cerr << "Sepaarator to existing name '" << _separator_to_existing << "'\n";

      continue;
    }

    if (n.starts_with("xrf="))
    {
      IWString fname = n;
      fname.remove_leading_chars(4);


      _cross_reference_file.open(fname.null_terminated_chars(), std::ios::out);

      if (! _cross_reference_file.good())
      {
        cerr << "Number_Assigner::initialise: cannot open cross reference file '" << fname << "'\n";
        return 0;
      }

      if (verbose)
        cerr << "Cross reference file '" << fname << "' opened\n";

      continue;
    }

    if ("help" == n)
    {
      display_standard_number_assigner_options(cerr, nflag);
      exit(1);
    }

    _prefix_string = n;
    rc++;
  }

// If they just entered a letter, start them with 1

  if (rc && -1 == _next_number_to_assign)
    _next_number_to_assign = 1;

  if (verbose)
    cerr << "Molecules will be assigned numbers starting with '" << _prefix_string << '(' << _next_number_to_assign << ")\n";

  return 1;
}

int
Number_Assigner::initialise(int i)
{
  if (i < 0)
  {
    cerr << "Number_Assigner::initialise: the initial value must be >= 0, " << i << " not allowed\n";
    return 0;
  }

  _next_number_to_assign = i;

  return i;
}

int
Number_Assigner::process(Molecule & m)
{
  IWString tmp(m.name());

  process(tmp);

  m.set_name(tmp);

  return 1;
}

int
Number_Assigner::process(IWString & s)
{
  if (_next_number_to_assign < 0)
    return 1;

  if (_only_apply_if_existing_name_is_empty && s.length() > 0)
    return 1;

  IWString buffer;
  buffer.resize(256);    // just a guess
  buffer = _prefix_string;
  if (_include_parentheses)
    buffer += '(';

  if (_number_digits > 0)
  {
    IWString tmp;
    tmp.append_number(_next_number_to_assign++);

    if (tmp.length() < _number_digits)
      tmp.shift(_number_digits - tmp.length(), '0');

    buffer << tmp;
  }
  else
    buffer.append_number(_next_number_to_assign++);

  if (_include_parentheses)
    buffer += ')';

  if (_cross_reference_file.rdbuf()->is_open())
  {
    _cross_reference_file << buffer << ' ' << s << '\n';
    if (! _cross_reference_file.good())
      cerr << "Number_Assigner::process: cross reference file corrupted\n";
  }

//cerr << "LIne " <<__LINE__ << " buffer '" <<  buffer << "' s '" << s << "'\n";
  if (_replace_name)
    ;
  else if (s.length())
    buffer << _separator_to_existing << s;

  s = buffer;

  return 1;
}

int
display_standard_number_assigner_options(std::ostream & os, char flag)
{
  os << "  -" << flag << " <number>    assign sequential numbers R(%d) starting with <number>\n";
  os << "  -" << flag << " <string>    use <string> rather than 'R'\n";
  os << "  -" << flag << " 'replace'   replace the entire name with R(nnn)\n";
  os << "  -" << flag << " xrf=fname   create a cross reference file\n";
  os << "  -" << flag << " noparen     omit parentheses when forming the names\n";
  os << "  -" << flag << " digits=<n>  left pad numbers with leading 0's to yield <n> digits\n";
  os << "  -" << flag << " nochange    only apply if the existing name is empty\n";
  os << "  -" << flag << " nsep='c'    separator between number and existing name (default ' ')\n";

  return os.good();
}
