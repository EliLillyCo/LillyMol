#include <stdlib.h>

#include "Foundational/mtrand/iwrandom.h"

#include "iw_tdt_filter.h"
#include "iw_tdt.h"

#if defined(SKIP_TDT_FILTER)
int
display_tdt_filter_syntax (std::ostream & os)
{
  os << "TDT Filter functionality not implemented this platform\n";

  return os.good ();
}
#else

void
IW_TDT_Filter::_default_values ()
{
  _tdts_examined = 0;

  _tdts_passing = 0;

  _tdts_with_tag_present = 0;

  _tag_absent_means_match = 0;

  return;
}

IW_TDT_Filter::IW_TDT_Filter ()
{
  _default_values ();

  return;
}

IW_TDT_Filter::~IW_TDT_Filter ()
{
  assert (ok ());

  return;
}

int
IW_TDT_Filter::ok () const
{
  if (0 == _f.number_elements () && 0 == _logexp.number_operators ())    // empty
    return 1;

  if (_f.number_elements () != _logexp.number_operators () + 1)
    return 0;

  return 1;
}

int
IW_TDT_Filter::debug_print (std::ostream & os) const
{
  os << "IW_TDT_Filter with " << _f.number_elements () << " components\n";

  for (int i = 0; i < _f.number_elements (); i++)
  {
    const IW_TDT_Filter_Base * f = _f[i];
    f->debug_print (os);
  }

  return os.good ();
}

int
IW_TDT_Filter::report (std::ostream & os) const
{
  os << "Filter with " << _f.number_elements () << " components examined " << _tdts_examined;
  if (_tdts_examined)
    os << ", " << _tdts_passing << " passed";
  os << endl;

  if (_tdts_with_tag_present < _tdts_examined)
    os << _tdts_with_tag_present << " TDT's had tag(s) present\n";

  return os.good ();
}

/*
  The specifier for a tdt filter will be a logical expression. Tokens
  look like
    FOO     dataitem FOO must be present
    !FOO    dataitem FOO must not be present
    3:FOO   the 3rd instance of FOO must be present
    BAR=3   the 1st BAR must be equal 3 (numeric test)
    BAR>4   the 1st BAR must be greater than 4 (numeric test)
    2:FOO<5 the 2nd FOO must be less than 5
    2:QQ=~string  the 2nd QQ must match 'string' as a regular expression

  These are combined in logical expressions as

    FOO&&2:BAR<4,,3:FOO=~^F

  Note that all operators must be repeated.
*/


/*
  Function to extract a grouping from a composite expression

  ZOP will be set to any operator found following the token
*/

static int
extract_next_component (const const_IWSubstring & s,
                        int & istart,
                        const_IWSubstring & zresult,
                        int & zop)
{
//cerr << "Looking for next component in '" << s << "' istart = " << istart << endl;

  if (istart >= s.length ())
    return 0;

  zop = IW_LOGEXP_UNDEFINED;

  int nchars = 0;

  int istop = istart;
  while (istop < s.length ())
  {
    if (s.matches_at_position (istop, "&&", 2))
    {
      zop = IW_LOGEXP_AND;
      nchars = 2;
      break;
    }
    else if (s.matches_at_position (istop, ",,", 2))
    {
      zop = IW_LOGEXP_OR;
      nchars = 2;
      break;
    }
    else if (s.matches_at_position (istop, "||", 2))
    {
      zop = IW_LOGEXP_OR;
      nchars = 2;
      break;
    }
    else if (s.matches_at_position (istop, "^^", 2))
    {
      zop = IW_LOGEXP_XOR;
      nchars = 2;
      break;
    }
    else if (s.matches_at_position (istop, ";;", 2))
    {
      zop = IW_LOGEXP_LOW_PRIORITY_AND;
      nchars = 2;
      break;
    }
    else if (s.matches_at_position (istop, ".AND.", 5))
    {
      zop = IW_LOGEXP_AND;
      nchars = 5;
      break;
    }
    else if (s.matches_at_position (istop, ".or.", 4) || s.matches_at_position (istop, ".OR.", 4))
    {
      zop = IW_LOGEXP_OR;
      nchars = 4;
      break;
    }
    else if (s.matches_at_position (istop, ".and.", 5))
    {
      zop = IW_LOGEXP_LOW_PRIORITY_AND;
      nchars = 5;
      break;
    }

    istop++;
  }

  s.from_to (istart, istop - 1, zresult);
  istart = istop + nchars;
  return 1;
}

/*
  Look for something like 'nn:XXX'
  and parse the nn into ITEM_NUMBER
*/

static int
extract_item_number (const_IWSubstring & s,
                     int & item_number)
{
  item_number = 0;     // the default value

  int istop = s.index (':');

  if (istop < 0)
    return 1;

  if (0 == istop)
  {
    cerr << "Token '" << s << "' starts with ':', this is surprising\n";
    return 1;
  }

  if (s.length () - 1 == istop)
  {
    cerr << "Token '" << s << "' ends with ':', this is surprising\n";
    return 1;
  }

  if (1 == istop && '*' == s[0])
  {
    item_number = -1;
    s.remove_leading_chars(2);
    return 1;
  }

  for (int i = 0; i < istop; i++)
  {
    if (s[0] < '0' || s[0] > '9')
    {
      cerr << "Non numeric character in number '" << s << "'\n";
      return 0;
    }

    item_number = item_number * 10 + s[0] - '0';
    s++;
  }

  if (item_number == 0)
  {
    cerr << "Item numbers start at 1, not 0\n";
    return 0;
  }

  s++;     // skip over the ':' character

  item_number--;     // our numbers start at 0

  return 1;
}

static int
extract_column_number (const_IWSubstring & s,
                       int & column)
{
  if (! s.starts_with (":COLUMN"))
    return 1;

  s.remove_leading_chars (7);

  if (0 == s.length ())
  {
    cerr << "Column specification cannot be end\n";
    return 0;
  }

  int chars_consumed = 0;

  while (s.length ())
  {
    char c = s[0];
    if (c >= '0' && c <= '9')
    {
      if (0 == chars_consumed)
        column = c - '0';
      else
        column = 10 * column + c - '0';

      chars_consumed++;
      s++;
    }
  }

  if (0 == chars_consumed)
  {
    cerr << "Invalid column specification\n";
    return 0;
  }

  return 1;
}

static int 
looks_like_regexp (const const_IWSubstring & s)
{
  return s.contains ("=~");
}

static int
looks_like_ne (const const_IWSubstring & s)
{
  return s.contains ("!=");
}

static int
looks_like_eq (const const_IWSubstring & s)
{
  return s.contains ('=');
}

static int
looks_like_lt (const const_IWSubstring & s)
{
  if (s.contains ('<'))
    return 1;

  if (s.contains (".lt."))
    return 1;

  if (s.contains (".LT."))
    return 1;

  return 0;
}

static int
looks_like_gt (const const_IWSubstring & s)
{
  if (s.contains ('>'))
    return 1;

  if (s.contains (".gt."))
    return 1;

  if (s.contains (".GT."))
    return 1;

  return 0;
}

static int
looks_like_random (const const_IWSubstring & s)
{
  return s.starts_with ("%RAND");
}

//#define DEBUG_BUILD_FROM_STRING

int
IW_TDT_Filter::build_from_string (const const_IWSubstring & s)
{
  assert (0 == _f.number_elements ());    // cannot already have been built - maybe we should just delete it...

  _f.resize (10);

  int istart = 0;

  if (s.starts_with ("TMMM"))
  {
    _tag_absent_means_match = 1;
    istart += 4;
  }

  int zop;
  const_IWSubstring token;

  int i = 0;
  while (extract_next_component (s, istart, token, zop))
  {
#ifdef DEBUG_BUILD_FROM_STRING
    cerr << "Examining token " << i << " '" << token << "'\n";
#endif

    if (token.starts_with ('!'))
    {
      _logexp.set_unary_operator (i, 0);
      token++;
    }

    if (IW_LOGEXP_UNDEFINED != zop)
    {
      if (! _logexp.add_operator (zop))
      {
        cerr << "IW_TDT_Filter::build_from_string: cannot add operator\n";
        return 0;
      }
    }

//  Does it look like '2:XXX...'

    int item_number = 0;
    if (! extract_item_number (token, item_number))
    {
      cerr << "IW_TDT_Filter::build_from_string: invalid numeric '" << s << "'\n";
      return 0;
    }

    int column = 1;
    if (! extract_column_number (token, column) || column < 1)
    {
      cerr << "IW_TDT_Filter::build_from_string: invalid column specification '" << s << "'\n";
      return 0;
    }

    column--;

    IW_TDT_Filter_Base * f = NULL;

    if (looks_like_regexp (token))
    {
      f = new IW_TDT_Filter_regexp ();
    }
    else if (looks_like_random (token))
    {
      f = new IW_TDT_Filter_random ();
    }
    else if (looks_like_eq (token))
    {
      f = new IW_TDT_Filter_eq ();
    }
    else if (looks_like_lt (token))
    {
      f = new IW_TDT_Filter_lt ();
    }
    else if (looks_like_gt (token))
    {
      f = new IW_TDT_Filter_gt ();
    }
    else
    {
      f = new IW_TDT_Filter_present ();
    }

    if (! f->construct (token))
    {
      cerr << "IW_TDT_Filter::build_from_string: cannot parse '" << token << "'\n";
      return 0;
    }

    f->set_item_to_process (item_number);

    _f.add (f);
    i++;
  }


  return 1;
}

int
IW_TDT_Filter::matches (const IW_TDT & tdt)
{
  _tdts_examined++;

  _logexp.reset ();

  for (int i = 0; i < _f.number_elements (); i++)
  {
    if (! _logexp.result_needed (i))
      continue;

    IW_TDT_Filter_Base * f = _f[i];

    int zresult = f->matches (tdt, _tag_absent_means_match);

    _logexp.set_result (i, zresult);

    int rc;
    if (! _logexp.evaluate (rc))    // need more results
      continue;

//  We got a result

    if (rc)
      _tdts_passing++;

    return rc;
  }

  return 0;
}

IW_TDT_Filter_Base::IW_TDT_Filter_Base ()
{
  _item_to_process = 0;

  _column_to_process = -1;

  _truncate_before_first = '\0';

  return;
}

IW_TDT_Filter_Base::~IW_TDT_Filter_Base ()
{
  return;
}

int
IW_TDT_Filter_Base::ok () const
{
  return IW_Regular_Expression::ok ();
}

int
IW_TDT_Filter_Base::_extract_data (const IW_TDT & tdt,
                                   const_IWSubstring & tdtvalue,
                                   int which_item_to_retrieve,
                                   int print_error_message)
{
  const_IWSubstring tag;

  if (! tdt.dataitem_value (*this, tag, tdtvalue, which_item_to_retrieve))
  {
    if (print_error_message)
      cerr << "IW_TDT_Filter_Base::matches: no " << _item_to_process << " '" << source () << "' item in tdt\n";
    return 0;
  }

  if (_column_to_process < 0)
    return 1;

  if (0 == _column_to_process)
  {
    tdtvalue.truncate_at_first (' ');
    return 1;
  }

  const_IWSubstring tmp;

  if (! tdtvalue.word (_column_to_process, tmp))
  {
    cerr << "IW_TDT_Filter_Base::_extract_data:cannot extract column " << (_column_to_process + 1) << "\n";
    return 0;
  }

  tdtvalue = tmp;

  return 1;
}

int
IW_TDT_Filter_Base::_extract_data (const IW_TDT & tdt,
                                   double & tdtvalue,
                                   int & tag_missing_from_tdt,
                                   int which_item_to_retrieve)
{
  const_IWSubstring stringval;

  if (! _extract_data (tdt, stringval, which_item_to_retrieve, 0))
  {
    tag_missing_from_tdt = 1;
    return 0;
  }

  if (! stringval.numeric_value (tdtvalue))
  {
    cerr << "IW_TDT_Filter_Base::_extract_data: non numeric value '" << stringval << "'\n";
    tag_missing_from_tdt = 0;
    return 0;
  }

  return 1;
}

int
IW_TDT_Filter_regexp::matches (const IW_TDT & tdt,
                      int return_code_if_absent)
{
  if (_item_to_process >= 0)
  {
    const_IWSubstring tdtvalue;
    if (! _extract_data (tdt, tdtvalue, _item_to_process))
    return return_code_if_absent;

    return _regexp.matches (tdtvalue);
  }

  if (-1 == _item_to_process)
  {
    int n = tdt.count_dataitems(*this);

    if (0 == n)
      return return_code_if_absent;

    for (int i = 0; i < n; i++)
    {
      const_IWSubstring tdtvalue;
      _extract_data(tdt, tdtvalue, i);
      if (_regexp.matches(tdtvalue))
        return 1;
    }
  }

  return 0;
}

int
IW_TDT_Filter_eq::matches (const IW_TDT & tdt,
                      int return_code_if_absent)
{
  if (_item_to_process >= 0)
  {
    double tdtvalue;
    int tag_missing_from_tdt;
    if (! _extract_data (tdt, tdtvalue, tag_missing_from_tdt, _item_to_process))
    {
      if (tag_missing_from_tdt)
        return return_code_if_absent;
    
      return 0;
    }

    return tdtvalue == _value;
  }

  if (-1 == _item_to_process)
  {
    int n = tdt.count_dataitems(*this);

    if (0 == n)
      return return_code_if_absent;

    for (int i = 0; i < n; i++)
    {
      double tdtvalue;
      int tag_missing_from_tdt;
      if (! _extract_data(tdt, tdtvalue, tag_missing_from_tdt, i))
        return 0;

      if (tdtvalue == _value)
        return 1;
    }
  }

  return 0;
}

int
IW_TDT_Filter_lt::matches (const IW_TDT & tdt,
                      int return_code_if_absent)
{
  if (_item_to_process >= 0)
  {
    double tdtvalue;
    int tag_missing_from_tdt;
    if (! _extract_data (tdt, tdtvalue, tag_missing_from_tdt, _item_to_process))
    {
      if (tag_missing_from_tdt)
        return return_code_if_absent;
    
      return 0;
    }

    return tdtvalue < _value;
  }

  if (-1 == _item_to_process)
  {
    int n = tdt.count_dataitems(*this);

    if (0 == n)
      return return_code_if_absent;

    for (int i = 0; i < n; i++)
    {
      double tdtvalue;
      int tag_missing_from_tdt;
      if (! _extract_data(tdt, tdtvalue, tag_missing_from_tdt, i))
        return 0;

      if (tdtvalue < _value)
        return 1;
    }
  }

  return 0;
}

int
IW_TDT_Filter_gt::matches (const IW_TDT & tdt,
                      int return_code_if_absent)
{
  if (_item_to_process >= 0)
  {
    double tdtvalue;
    int tag_missing_from_tdt;
    if (! _extract_data (tdt, tdtvalue, tag_missing_from_tdt, _item_to_process))
    {
      if (tag_missing_from_tdt)
        return return_code_if_absent;
    
      return 0;
    }

    return tdtvalue > _value;
  }

  if (-1 == _item_to_process)
  {
    int n = tdt.count_dataitems(*this);

    if (0 == n)
      return return_code_if_absent;

    for (int i = 0; i < n; i++)
    {
      double tdtvalue;
      int tag_missing_from_tdt;
      if (! _extract_data(tdt, tdtvalue, tag_missing_from_tdt, i))
        return 0;

      if (tdtvalue > _value)
        return 1;
    }
  }

  return 0;
}

int
IW_TDT_Filter_present::matches (const IW_TDT & tdt,
                      int return_code_if_absent)    // doesn't make sense here
{
  const_IWSubstring tdtvalue;
  return _extract_data (tdt, tdtvalue, 0);
}

int
IW_TDT_Filter_random::matches (const IW_TDT & tdt,
                      int return_code_if_absent)    // not needed here
{
  if (iwrandom () < _ratio)
    return 1;

  return 0;
}

int
IW_TDT_Filter_regexp::construct (const const_IWSubstring & s)
{
  int iop = s.index ("=~");

  if (iop <= 0 || iop == s.length () - 2)
  {
    cerr << "IW_TDT_Filter_regexp::construct: operator =~ badly positioned '" << s << "'\n";
    return 0;
  }

  const_IWSubstring tag;
  s.from_to (0, iop - 1, tag);
  tag.strip_trailing_blanks ();

  if (! IW_Regular_Expression::set_pattern (tag))
  {
    cerr << "IW_TDT_Filter_regexp::construct: cannot parse tag regexp '" << tag << "'\n";
    return 0;
  }

  const_IWSubstring m;
  s.from_to (iop + 2, s.length () - 1, m);

  cerr << "IW_TDT_Filter_regexp::construct: tag is '" << tag << "' matcher is '" << m << "'\n";

  if (! _regexp.set_pattern (m))
  {
    cerr << "IW_TDT_Filter_regexp::construct: cannot parse regexp '" << m << "'\n";
    return 0;
  }

  return 1;
}

int
IW_TDT_Filter_present::construct (const const_IWSubstring & s)
{
  if (! IW_Regular_Expression::set_pattern (s))
  {
    cerr << "IW_Regular_Expression::construct: cannot parse tag rx '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
IW_TDT_Filter_Base::_extract_numeric (const const_IWSubstring & s,
                                      const char * op,
                                      double & zresult)
{
  unsigned int iop = s.find (op);

  if (iop <= 0 || s.length () - ::strlen (op) == iop)
  {
    cerr << "IW_TDT_Filter_Base::_extract_numeric: operator '" << op << "' missing or incorrectly placed in '" << s << "'\n";
    return 0;
  }

  const_IWSubstring tag;
  s.from_to (0, iop - 1, tag);
  tag.strip_trailing_blanks ();

  if (! IW_Regular_Expression::set_pattern (tag))
  {
    cerr << "IW_TDT_Filter_Base::_extract_numeric: cannot parse tag rx '" << tag << "'\n";
    return 0;
  }

  const_IWSubstring v;
  s.from_to (iop + ::strlen (op), s.length () - 1, v);

  if (! v.numeric_value (zresult))
  {
    cerr << "IW_TDT_Filter_Base::_extract_numeric: invalid numeric '" << v << "' in '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
IW_TDT_Filter_Base::_extract_numeric (const const_IWSubstring & s,
                                      char op,
                                      double & zresult)
{
  int iop = s.find (op);

  if (iop <= 0 || s.length () - 1 == iop)
  {
    cerr << "IW_TDT_Filter_Base::_extract_numeric: operator '" << op << "' missing or incorrectly placed in '" << s << "'\n";
    return 0;
  }

  const_IWSubstring tag;
  s.from_to (0, iop - 1, tag);
  tag.strip_trailing_blanks ();

  if (! IW_Regular_Expression::set_pattern (tag))
  {
    cerr << "IW_TDT_Filter_Base::_extract_numeric: cannot parse tag rx '" << tag << "'\n";
    return 0;
  }

  const_IWSubstring v;
  s.from_to (iop + 1, s.length () - 1, v);
  v.strip_leading_blanks ();

  if (! v.numeric_value (zresult))
  {
    cerr << "IW_TDT_Filter_Base::_extract_numeric: invalid numeric '" << v << "' in '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
IW_TDT_Filter_eq::construct (const const_IWSubstring & s)
{
  return _extract_numeric (s, '=', _value);
}

int
IW_TDT_Filter_gt::construct (const const_IWSubstring & s)
{
  if (s.contains ('>'))
    return _extract_numeric (s, '>', _value);

  if (s.contains (".gt."))
    return _extract_numeric (s, ".gt.", _value);

  if (s.contains (".GT."))
    return _extract_numeric (s, ".GT.", _value);

  abort ();

  return 0;
}

int
IW_TDT_Filter_lt::construct (const const_IWSubstring & s)
{
  if (s.contains ('<'))
    return _extract_numeric (s, '<', _value);

  if (s.contains (".lt."))
    return _extract_numeric (s, ".lt.", _value);

  if (s.contains (".LT."))
    return _extract_numeric (s, ".LT.", _value);

  abort ();

  return 0;
}

/*
  The random specifier looks like
    %RAND=0.5
*/

int
IW_TDT_Filter_random::construct (const const_IWSubstring & s)
{
  const_IWSubstring mys = s;
  mys.remove_leading_chars (6);    // get rid of '%RAND='

  if (! mys.numeric_value (_ratio) || _ratio <= 0.0 || _ratio >= 1.0)
  {
    cerr << "IW_TDT_Filter_random::construct: invalid ratio '" << s << "'\n";
    return 0;
  }

  iw_random_seed ();

  return 1;
}

int
IW_TDT_Filter_regexp::ok () const
{
  if (!  _regexp.ok ())
    return ok ();

  return IW_TDT_Filter_Base::ok ();
}

int
IW_TDT_Filter_present::ok () const
{
  return IW_TDT_Filter_Base::ok ();
}

int
IW_TDT_Filter_eq::ok () const
{
  return IW_TDT_Filter_Base::ok ();
}

int
IW_TDT_Filter_lt::ok () const
{
  return IW_TDT_Filter_Base::ok ();
}

int
IW_TDT_Filter_gt::ok () const
{
  return IW_TDT_Filter_Base::ok ();
}

int
IW_TDT_Filter_random::ok () const
{
  return (_ratio > 0.0 && _ratio < 1.0);
}

int
IW_TDT_Filter_present::debug_print (std::ostream & os) const
{
  cerr << "filter for presence of " << _item_to_process << " instance of '" << IW_Regular_Expression::source () << "'\n";
  return os.good ();
}
int
IW_TDT_Filter_regexp::debug_print (std::ostream & os) const
{
  cerr << "filter " << _item_to_process << " instance of '" << IW_Regular_Expression::source () << "' matches '" << _regexp.source () << "'\n";

  return os.good ();
}
int
IW_TDT_Filter_eq::debug_print (std::ostream & os) const
{
  os << "filter " << _item_to_process << " instance of '" << IW_Regular_Expression::source () << "' = " << _value << endl;

  return os.good ();
}
int
IW_TDT_Filter_lt::debug_print (std::ostream & os) const
{
  os << "filter " << _item_to_process << " instance of '" << IW_Regular_Expression::source () << "' < " << _value << endl;

  return os.good ();
}
int
IW_TDT_Filter_gt::debug_print (std::ostream & os) const
{
  os << "filter " << _item_to_process << " instance of '" << IW_Regular_Expression::source () << "' > " << _value << endl;

  return os.good ();
}

int
IW_TDT_Filter_random::debug_print (std::ostream & os) const
{
  os << "random filter, ratio = " << _ratio << endl;

  return os.good ();
}

int
display_tdt_filter_syntax (std::ostream & os)
{
  os << "A TDT filter is used to select some TDT's from a set of TDT's\n";
  os << "Each entry in a TDT (Thor DataTree) file looks like\n";
  os << "PCN<FOO>\n";
  os << "In this case the TAG is 'PCN' and the value is 'FOO'\n";

  os << "You must provide the name of the TAG, and optionally\n";
  os << "some condition for the value\n";

  os << "Entering 'FOO' is a simple match, that is, a dataitem of\n";
  os << "type 'FOO' must be present\n";
  os << "Since this is a regular expression match, FOO matches 'XXFOOYYY'\n";
  os << "Use '^FOO$' to match FOO and just FOO\n";
  os << "Something like '2:FOO' means that a second occurrence of\n";
  os << "something matching FOO must be present\n";
  os << "FOO>3.1 matches only TDT's for which the VALUE of dataitem\n";
  os << "matching FOO is greater than 3.1\n";
  os << "while '3:FOO<5' means that the 3rd dataitem to match FOO must be less than 5\n";
  os << "Regular expressions are also available\n";
  os << "'2:FOO=~^BAR' will match those TDT's whose TAG's match 'FOO' and\n";
  os << "whose VALUE starts with 'BAR'\n";
  os << endl;

  os << "A specific column can be designated by following the tag with :COLUMN\n";

  os << "These can be combined into regular expressions\n";
  os << "For example 'FOO&&BAR' means that a TAG which matches FOO and a\n";
  os << "tag which matches BAR must both be present\n";
  os << "While '2:FOO<3,,BAR=~^XX' means the 2nd occurrence of something\n";
  os << "whose tag matches FOO must be less than 3, OR (,,) something whose\n";
  os << "TAG matches BAR must start with 'XX'\n";


  os << endl;

  return os.good ();
}

#endif
