#include <stdlib.h>

/*
  Tversky specifications
*/

#include "cmdline.h"

#include "tversky.h"

static const char * cvsid = "tversky.cc,v 1.5 2002/06/21 19:05:32";

Tversky::Tversky ()
{
  (void) cvsid;

  _a = -99.9;
  _b = -99.9;

  _active = 0;

  _optimistic_mode = 0;

  _treat_non_colliding_as_01 = 0;

  _nset_sensitive_zero_bit_similarity = 0;

  return;
}

int
display_standard_tversky_options (std::ostream & os)
{
  os << " -a <tva>         Tversky A coefficient\n";
  os << " -b <tva>         Tversky B coefficient\n";

  return os.good ();
}

int
display_standard_tversky_options (std::ostream & os, char flag)
{
  os << " -" << flag << " a=nn          Tversky A coefficient\n";
  os << " -" << flag << " b=nn          Tversky B coefficient\n";
  os << " -" << flag << " optimistic    Optimistic hybrid Tversky/Tanimoto distances\n";
  os << " -" << flag << " noscale       normally parameters are scaled to a+b=2\n";
  os << " -" << flag << " nc01          treat non-colliding counted fingerprints as 0/1 fingerprints\n";
  os << " -" << flag << " nset0s        distance when 0 bits set depends on nset in other fingerprint\n";

  return os.good ();
}

int
Tversky::parse_command_line (Command_Line & cl,
                             char aflag, char bflag, 
                             int verbose)
{
  if (! cl.option_present (aflag) && ! cl.option_present (bflag))    // nothing specified
  {
    _active = 0;
    return 1;
  }

  if (! cl.value (aflag, _a) || _a < 0.0)
  {
    cerr << "Missing or invalid value of Tversky A parameter\n";
    return 0;
  }

  if (! cl.value (bflag, _b) || _b < 0.0)
  {
    cerr << "Missing or invalid value of Tversky B parameter\n";
    return 0;
  }

  if (verbose)
    cerr << "Tversky parameters set to " << _a << " and " << _b << endl;

  _active = 1;

  return 1;
}

int
Tversky::parse_command_line (Command_Line & cl,
                             char flag,
                             int verbose)
{
  if (! cl.option_present (flag))
  {
    _active = 0;
    return 1;
  }

  int scale_to_2 = 1;

  int i = 0;
  const_IWSubstring f;
  while (cl.value (flag, f, i++))
  {
    if (f.starts_with ("opt"))
    {
      _optimistic_mode = 1;

      if (verbose)
        cerr << "Optimistic Tversky mode\n";
    }
    else if (f.starts_with ("a="))
    {
      f.remove_leading_chars (2);

      if (! f.numeric_value (_a) || _a < 0.0)
      {
        cerr << "INvalid Tversky A parameter '" << f << "'\n";
        return 0;
      }
    }
    else if (f.starts_with ("b="))
    {
      f.remove_leading_chars (2);

      if (! f.numeric_value (_b) || _b < 0.0)
      {
        cerr << "INvalid Tversky B parameter '" << f << "'\n";
        return 0;
      }
    }
    else if ("noscale" == f)
      scale_to_2 = 0;
    else if ("nc01" == f)
    {
      _treat_non_colliding_as_01 = 1;
      if (verbose)
        cerr << "Non colliding sparse fingerprints treated as 0/1 Tversky\n";
    }
    else if ("nset0s" == f)
    {
      _nset_sensitive_zero_bit_similarity = 1;
      if (verbose)
        cerr << "Zero bit fingerprints treated specially for Tversky\n";
    }
    else if ("help" == f)
    {
      display_standard_tversky_options (cerr, flag);
      exit (1);
    }
    else
    {
      cerr << "Unrecognised tversky qualifier '" << f << "'\n";
      return 0;
    }
  }

  if (_a >= 0.0 && _b >= 0.0)     // great, both specified
  {
    if (scale_to_2)
    {
      double sum = _a + _b;
      _a = _a / sum * 2.0;
      _b = _b / sum * 2.0;
    }
  }
  else if (_a < 0.0 && _b >= 0.0 && _b <= 2.0)
  {
    _a = 2.0 - _b;
  }
  else if (_a >= 0.0 && _a <= 2.0 && _b < 0.0)
  {
    _b = 2.0 - _a;
  }
  else
  {
    cerr << "tversky::parse_command_line: incomplete specification a = " << _a << " b = " << _b << "\n";
    return 0;
  }

  if (verbose)
  {
    cerr << "Tversky parameters set to " << _a << " and " << _b;
    if (_optimistic_mode)
      cerr << ", optimistic distances.\n";

    cerr << endl;
  }

  _active = 1;

  return 1;
}

similarity_type_t
Tversky::tanimoto (int nset1, int nset2,
                   int bic) const
{
  if (0 == nset1 || 0 == nset2)
  {
    if (0 == nset1 && 0 == nset2)
      return static_cast<similarity_type_t> (1.0);

    if (! _nset_sensitive_zero_bit_similarity)
      return static_cast<similarity_type_t> (0.0);

    if (nset1)
      return static_cast<similarity_type_t> (1.0) / static_cast<similarity_type_t> (nset1 + 1);
    else
      return static_cast<similarity_type_t> (1.0) / static_cast<similarity_type_t> (nset2 + 1);
  }

// Use Bradshaw's notation

  int a = nset1 - bic;
  int b = nset2 - bic;

//cerr << "a = " << a << " b = " << b << " bic = " << bic << endl;
  return (float (bic) / (_a * float (a) + _b * float (b) + float (bic)));
}
