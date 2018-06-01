#include <stdlib.h>
#include <ctype.h>
#include <iostream>

#include "ematch.h"
#include "iwstring.h"

void
Element_Matcher::_default_values ()
{
  _e = NULL;
  _isotope = -1;

  _match_organic_only = 0;
  _match_non_organic_only = 0;
  _match_non_periodic_only = 0;

   return;
}

Element_Matcher::Element_Matcher ()
{
  _default_values();

  return;
}

Element_Matcher::Element_Matcher (const Element * e)
{
  _default_values();

  set_element(e);

  return;
}

/*
  If the leading character is *, then we do an atomic number
  match. So, '*Ra' matches any isotopic variant of Radium
*/

Element_Matcher::Element_Matcher (const char * s)
{
  _default_values();

  (void) construct_from_string(s, static_cast<int>(strlen(s)));

  return;
}

Element_Matcher::Element_Matcher (const IWString & s)
{
  _default_values();

  (void) construct_from_string(s.rawchars(), s.nchars());

  return;
}

std::ostream &
operator<< (std::ostream & os, const Element_Matcher & em)
{
  em.operator_less_less(os);

  return os;
}

int
Element_Matcher::operator_less_less (std::ostream & os) const
{
  os << "Element Matcher:match";

  if (_isotope >= 0)
    os << " isotope " << _isotope;
     
  if (_e)
  {
    os << " element '" << _e->symbol() << "'";
    return os.good();
  }

  if (_match_organic_only)
  {
    os << " organic";
    return os.good();
  }
  
  if (_match_non_organic_only)
  {
    os << " nonorganic";
    return os.good();
  }

  if (_match_non_periodic_only)
  {
    os << " nonperiodic";
    return os.good();
  }

  os << " not sure what I'm supposed to match";

  return os.good();
}

void
display_element_matcher_syntax (std::ostream & os)
{
  os << "Element_Matcher recognised syntax options\n";
  os << " RX=<regexp>      elements whose symbols match <regexp>, e.g. 'RX=^[C,O]$'\n";
  os << " norganic         organic elements\n";
  os << " nonorganic       non organic elements\n";
  os << " nonperiodic      elements not in the periodic table\n";
  os << " <sym>            match the element with symbol <sym>\n";
  os << " <n><sym>         match isotope <n> of <sym>\n";
//os << " *<sym>           match all isotopic variants of <sym>\n";   not sure this works properly, it just converts to atomic number match, which isn't what should happen

  return;
}

//#define DEBUG_EMATCH_CONSTRUCT_FROM_STRING

/*
  Need to be careful with the meaning of *
  The problem is that it is not only the wildcard for all isotopic
  variants, but also an element by itself
*/

int
Element_Matcher::construct_from_string (const const_IWSubstring & directive)
{
  const_IWSubstring s(directive);

  if (0 == s.length())
    return 1;

  _isotope = -1;

  while (s.length() > 0 && isdigit(s[0]))
  {
    if (_isotope < 0)
      _isotope = s[0] - '0';
    else
      _isotope = 10 * _isotope + s[0] - '0';

    s++;
  }

  if (0 == s.length())   // just an isotope specification
    return 1;

  if ("organic" == s)
  {
    _match_organic_only = 1;
    return 1;
  }

  if ("nonorganic" == s)
  {
    _match_non_organic_only = 1;
    return 1;
  }

  if ("nonperiodic" == s)
  {
    _match_non_periodic_only = 1;
    return 1;
  }

// We need a way of specifying the * element itself

  if ('*' == directive)
  {
    const Element * e = get_element_from_symbol_no_case_conversion("*");
    set_element(e);

    return 1;
  }

  if (s.starts_with("RX="))
  {
    s += 3;
    if (! _symbol_rx.set_pattern(s))
    {
      cerr << "Element_Matcher::construct_from_string:invalid symbol regular expression '" << s << "'\n";
      return 0;
    }

    return 1;
  }

#ifdef DEBUG_EMATCH_CONSTRUCT_FROM_STRING
  cerr << "Element matcher making from '" << s << "'\n";
#endif

  if (! isalnum(s[0]))
  {
    cerr << "Element_Matcher::construct_from_string: elements must start with a letter or number '" << s << "'\n";
    return 0;
  }

  const Element * e = get_element_from_symbol(s, _isotope);

  if (NULL == e)
    e = create_element_with_symbol(s);

  set_element(e);

  return 1;
}

int
Element_Matcher::construct_from_string (const char * s, int lens)
{
  const const_IWSubstring tmp(s, lens);

  return construct_from_string(tmp);
}

int
Element_Matcher::construct_from_string (const char * s)
{
  const_IWSubstring tmp(s);

  return construct_from_string(tmp);
}

int
Element_Matcher::construct_from_string (const IWString & s)
{
  const const_IWSubstring tmp(s);

  return construct_from_string(tmp);
}

Element_Matcher::Element_Matcher (atomic_number_t z)
{
  _default_values();

  const Element * e = get_element_from_atomic_number(z);
  assert(e);

  set_element(e);

  return;
}

int
Element_Matcher::ok() const
{
  if (NULL == _e)
    ;
  else if (! _e->ok())
    return 0;

  return 1;
}

int
Element_Matcher::debug_print (std::ostream & os) const
{
  assert (os.good());

  os << "Element matcher ";
  if (_isotope >= 0)
    os << "isotope " << _isotope << ' ';
  if (NULL != _e)
    os << "element '" << _e->symbol() << "'";

  os << endl;

  if (_match_non_periodic_only)
    os << "Non periodic table elements\n";
  if (_match_organic_only)
    os << "Organic only\n";
  if (_match_non_organic_only)
    os << "Non organic only\n";

  if (_symbol_rx.active())
    os << " matched symbol rx '" << _symbol_rx.source() << "'\n";

  return os.good();
}

void
Element_Matcher::set_element (const Element * e)
{
  assert (NULL != e);

  assert (e->ok());

  _e = e;

  return;
}

//#define DEBUG_ELEMENT_MATCHER_MATCHES

int
Element_Matcher::matches (const Element * e, int iso)
{
  assert (e->ok());

#ifdef DEBUG_ELEMENT_MATCHER_MATCHES
  cerr << "Trying to match '" << e->symbol() << "' iso " << iso << endl;
  debug_print(cerr);
#endif

  int isotope_matched = 0;   // we need to record the fact that something matched

  if (_isotope < 0)    // not active, do not check
    ;
  else if (iso != _isotope)
    return 0;
  else if (NULL == _e)    // need to check other attributes
    return isotope_matched = 1;

#ifdef DEBUG_ELEMENT_MATCHER_MATCHES
  if (_e == e)
    cerr << "Element matcher returning 1 on exact match\n";
#endif

  if (NULL == _e)    // do not check
    ;
  else if (_e == e)
    return 1;

  if (_match_organic_only)
    return e->organic();

  if (_match_non_organic_only)
    return ! e->organic();

  if (_match_non_periodic_only)
    return ! e->is_in_periodic_table();

  if (_symbol_rx.active())
    return _symbol_rx.matches(e->symbol());

  return isotope_matched;
}

#include "cmdline.h"

/*
  Note special treatment for the zero length string.
  This for fileconv's -O switch. -O must take an argument, but if
  the regular "organic" set is to be used, one invokes it as -O ''
  Mar 2003. Because of problems with quotes and such, recognise the
  special case of 'none'
*/

int
Set_of_Element_Matches::construct_from_command_line (Command_Line & cl,
                               int verbose,
                               char mflag)
{
  IWString ele;
  int i = 0;
  while (cl.value(mflag, ele, i++))
  {
    if (0 == ele.length() || "none" == ele)    // special case(s) for fileconv
      continue;

    if (1 == ele.nwords())
    {
      Element_Matcher * em = new Element_Matcher(ele);

      if (verbose)
        cerr << "Element '" << ele << "' included in element matches\n";

      add(em);
    }
    else
    {
      int j = 0;
      const_IWSubstring token;
      while (ele.nextword(token, j))
      {
        Element_Matcher * em = new Element_Matcher(token);
        if (verbose)
          cerr << "Element '" << (*em) << "' included in element matches\n";

        add(em);
      }
    }
  }

  return _number_elements + 1;
}

int
Set_of_Element_Matches::matches (const Element * e, int iso)
{
  for (int i = 0; i < _number_elements; i++)
  {
    Element_Matcher * em = _things[i];

    if (em->matches(e, iso))
      return i + 1;
  }

  return 0;
}
