#include <stdlib.h>

#include "valid_descriptor_name.h"

static int default_max_name_length = 32;

int
display_standard_valid_descriptor_name_options (std::ostream & os,
                        char mflag, 
                        char vflag)
{
  os << " -" << mflag << " <number>    maximum length of the name of a descriptor (default " << default_max_name_length << ")\n";
  os << " -" << vflag << " <char>      <char> is valid in a descriptor name\n";

  return os.good ();
}

Valid_Descriptor_Name::Valid_Descriptor_Name ()
{
  _max_name_length = default_max_name_length;

  for (int i = 0; i < 256; i++)
  {
    _valid_character[i] = 0;
  }

  for (int i = 0; i <= 9; i++)
  {
    int j = static_cast<int> ('0' + i);
    _valid_character[j] = 1;
  }

  for (int i = 0; i < 26; i++)
  {
    int j = static_cast<int> ('a' + i);
    _valid_character[j] = 1;
    j = static_cast<int> ('A' + i);
    _valid_character[j] = 1;
  }

  int i = static_cast<int> ('_');
  _valid_character[i] = 1;

  i = static_cast<int> (':');    // Jan 2008
  _valid_character[i] = 1;

  return;
}

int
Valid_Descriptor_Name::construct_from_command_line (Command_Line & cl,
                                                    char mflag,
                                                    char vflag,
                                                    int verbose)
{
  if (cl.option_present (mflag))
  {
    if (! cl.value (mflag, _max_name_length) || _max_name_length < 1)
    {
      cerr << "The maximum length of a descriptor name must be a whole positive number\n";
      return 0;
    }

    if (verbose)
      cerr << "Descriptor names must be no longer than " << _max_name_length << " characters\n";
  }

  if (cl.option_present (vflag))
  {
    int i = 0;
    const_IWSubstring v;
    while (cl.value (vflag, v, i++))
    {
      int j = 0;
      const_IWSubstring token;
      while (v.nextword (token, j, ','))
      {
        if (1 != token.length ())
        {
          cerr << "Valid descriptor name characters must be of length 1, '" << token << "' not valid\n";
          return 0;
        }

        int k = static_cast<int> (token[0]);
        _valid_character[k] = 1;

        if (verbose)
          cerr << "'" << token << "' valid in descriptor names\n";
      }
    }
  }

  return 1;
}

int
Valid_Descriptor_Name::valid (const IWString & d) const
{
  IWString not_used;

  return valid (d, not_used);
}

int
Valid_Descriptor_Name::valid (const IWString & d,
                              IWString & reason) const
{
  int rc = 1;

  if (d.length () > _max_name_length)
  {
    reason << "name too long";
    rc = 0;
  }

  for (int j = 0; j < d.length (); j++)
  {
    char c = d[j];
    if (0 == _valid_character[static_cast<unsigned int> (c)])
    {
      reason << " invalid character '" << c << "'";
      rc = 0;
    }
  }

  return rc;
}
