#include <stdlib.h>
#include <string>

#include "cmdline.h"

template <>
int
Command_Line::value<std::string>(const char c, std::string& zvalue, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
      {
        if (nullptr == oo->value())
          zvalue = "";
        else
          zvalue = oo->value();
        return 1;
      }

      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value_as_std_string (const char c, std::string & zvalue, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
      {
        if (nullptr == oo->value())
          zvalue = "";
        else
          zvalue = oo->value();
        return 1;
      }

      nfound++;
    }
  }

  return 0;
}

std::string
Command_Line::std_string_value(const char c, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
      {
        const char * v = oo->value ();
        if (nullptr == v)
          return "";
        else
          return oo->value ();
      }

      nfound++;
    }
  }

  return "";
}
