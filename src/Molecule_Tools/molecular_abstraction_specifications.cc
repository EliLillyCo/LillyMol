#include <stdlib.h>

#include "molecular_abstraction_specifications.h"

#define MA_OPEN_PAREN '('
#define MA_CLOSE_PAREN ')'

Molecular_Abstraction_Directives_Node::Molecular_Abstraction_Directives_Node()
{
  _type = 0;

  _next = NULL;

  return;
}

Molecular_Abstraction_Directives_Node::~Molecular_Abstraction_Directives_Node()
{
  if (NULL != _next)
    delete _next;

  return;
}

/*
  Input looks like

  xxxx.yyyyy
  xxxx().yyyy()
  xxxx()yyyy.zzz
*/

static int
next_directive(const const_IWSubstring & s,
               IWString & zdirective,
               IWString & args)
{
//cerr << "Extracting next directive from '" << s << "'\n";

  int nchars = s.length();

  if (0 == nchars)
    return 0;

  int rc = 0;    // number characters consumed

  int ndx = 0;

  while (ndx < nchars)
  {
    char c = s[ndx];

    ndx++;
    rc++;

    if (MA_CLOSE_PAREN == c)
    {
      cerr << "Cannot have closing paren in directive\n";
      return 0;
    }

    if (MA_OPEN_PAREN == c)
      break;
    else if ('.' == c)
      return rc;

    if (isalnum(c))
      ;
    else if ('_' == c)
      ;
    else
    {
      cerr << "Directives can be alphanumeric only '" << c << "' not allowed\n";
      return 0;
    }

    zdirective.add(c);
  }

  if (0 == zdirective.length())
  {
    cerr << "No directive\n";
    return 0;
  }

  if (ndx == nchars)
    return rc;

  int paren_level = 1;

  while (ndx < nchars)
  {
    char c = s[ndx];

    ndx++;
    rc++;

    if (MA_OPEN_PAREN == c)
      paren_level++;
    else if (MA_CLOSE_PAREN == c)
    {
      paren_level--;
      if (0 == paren_level)
        return rc;
    }

    args.add(c);
  }

  cerr << "No closing paren\n";
  return 0;
}

int
Molecular_Abstraction_Directives_Node::build(const const_IWSubstring & s)
{
  if (0 == s.length())
  {
    cerr << "Molecular_Abstraction_Directives_Node::build:empty specification, cannot process\n";
    return 0;
  }

  const_IWSubstring tmp(s);
  tmp.strip_trailing_blanks();

  int rc =  _build(tmp);

  if (! rc)
  {
    cerr << "Molecular_Abstraction_Directives_Node::build:cannot parse '" << tmp << "'\n";
    return 0;
  }

  return rc;
}

int
Molecular_Abstraction_Directives_Node::_build(const const_IWSubstring & s)
{
  assert (s.length() > 0);

#ifdef DEBUG_MADN_BUILD
  cerr << "_build from '" << s << "'\n";
#endif

  int chars_consumed = next_directive(s, _directive, _args);

  if (0 == chars_consumed)
    return 0;

  if (chars_consumed == s.length())
    return 1;

  const_IWSubstring tmp(s);

  tmp.remove_leading_chars(chars_consumed);

  if (tmp.starts_with('.'))
  {
    if (1 == tmp.length())
      return 1;

    tmp.remove_leading_chars(1);
  }

  _next = new Molecular_Abstraction_Directives_Node;

//cerr << "tmp is '" << tmp << "'\n";
  return _next->_build(tmp);
}

int
Molecular_Abstraction_Directives_Node::_extract_directive_and_args(const const_IWSubstring & s)
{
  int open_paren = -1;

  for (int i = 0; i < s.length(); i++)
  {
    if (MA_CLOSE_PAREN == s[i])
    {
      cerr << "Molecular_Abstraction_Directives_Node::_extract_directive_and_args:close paren found '" << s << "'\n";
      return 0;
    }

    if (MA_OPEN_PAREN == s[i])
    {
      open_paren = i;
      break;
    }
  }

  if (open_paren < 0)
  {
    cerr << "Molecular_Abstraction_Directives_Node::_extract_directive_and_args:no opening paren '" << s << "'\n";
    return 0;
  }

  s.from_to(0, open_paren - 1, _directive);

  int close_paren = -1;

  for (int i = open_paren + 1; i < s.length(); i++)
  {
    if (MA_OPEN_PAREN == s[i])
    {
      cerr << "Molecular_Abstraction_Directives_Node::_extract_directive_and_args:multiple open parens '" << s << "'\n";
      return 0;
    }

    if (MA_CLOSE_PAREN == s[i])
    {
      close_paren = i;
      break;
    }
  }

  if (close_paren < 0)
  {
    cerr << "Molecular_Abstraction_Directives_Node::_extract_directive_and_args:no closing paren '" << s << "'\n";
    return 0;
  }

  if (open_paren + 1 < close_paren)
    s.from_to(open_paren + 1, close_paren - 1, _args);

  if (close_paren == s.length() - 1)
    return s.length();

  const_IWSubstring tmp(s);

  tmp.remove_leading_chars(close_paren);

  return close_paren + 1;
}

int
Molecular_Abstraction_Directives_Node::debug_print(ostream & os) const
{
  os << "Directive '" << _directive << "', args '" << _args << "'\n";
  if (NULL != _next)
  {
    os << " next\n";
    return _next->debug_print(os);
  }

  return 1;
}

int
Molecular_Abstraction_Directives_Node::number_abstractions() const
{
  if (NULL == _next)
    return 1;

  return 1 + _next->number_abstractions();
}

/*
  We have read the abstraction specification, now we need to know if
  it contains recognised terms
*/

int
Molecular_Abstraction_Directives_Node::directive_recognised()
{
  if (MAD_RMSPINACH == _directive)
    _type = MAD_TYPE_RMSPINACH;
  else if (MAD_RMSCAFFOLD == _directive)
    _type = MAD_TYPE_RMSCAFFOLD;
  else if (MAD_SCAFFOLD == _directive)
    _type = MAD_TYPE_SCAFFOLD;
  else if (MAD_RINGS == _directive)
    _type = MAD_TYPE_RINGS;
  else if (MAD_BIGRING == _directive)
    _type = MAD_TYPE_BIGRING;
  else if (MAD_TRANSLATE == _directive)
    _type = MAD_TYPE_TRANSLATE;
  else if (MAD_REMOVE_ATOM == _directive)
    _type = MAD_TYPE_REMOVE_ATOM;
  else if (MAD_CHANGE_BOND_TYPE == _directive)
    _type = MAD_TYPE_CHANGE_BOND_TYPE;
  else if (MAD_REPLACE_LINKER == _directive)
    _type = MAD_TYPE_REPLACE_LINKER;
  else if (MAD_ABSTRACT_RING_FORM == _directive)
    _type = MAD_TYPE_ABSTRACT_RING_FORM;
  else if (MAD_FRAGMENT == _directive)
    _type = MAD_TYPE_FRAGMENT;
  else if (MAD_REMOVE_ATOMS == _directive)
    _type = MAD_TYPE_REMOVE_ATOMS;
  else if (MAD_PLACE_ISOTOPE == _directive)
    _type = MAD_TYPE_PLACE_ISOTOPE;
  else if (MAD_PLACE_CHARGE == _directive)
    _type = MAD_TYPE_PLACE_CHARGE;
  else if (MAD_ALL_ATOMS_TRANSFORM == _directive)
    _type = MAD_TYPE_ALL_ATOMS_TRANSFORM;
  else if (MAD_ALL_BONDS_TRANSFORM == _directive)
    _type = MAD_TYPE_ALL_BONDS_TRANSFORM;
  else if (MAD_COMPRESS_CONSECUTIVE == _directive)
    _type = MAD_TYPE_COMPRESS_CONSECUTIVE;
  else if (MAD_RINGSYS == _directive)
    _type = MAD_TYPE_RINGSYS;
  else if (MAD_RMBOND == _directive)
    _type = MAD_TYPE_RMBOND;
  else if (MAD_RMRD2 == _directive)
    _type = MAD_TYPE_RMRD2;
  else if (MAD_INVSCAF == _directive)
    _type = MAD_TYPE_INVSCAF;
  else if (MAD_SSS == _directive)
    _type = MAD_TYPE_SSS;
  else if (MAD_SPINACH == _directive)
    _type = MAD_TYPE_SPINACH;
  else
  {
    cerr << "Molecular_Abstraction_Directives_Node::directive_recognised:unrecognised directive '" << _directive << "'\n";
    return 0;
  }

  if (NULL != _next)
    return _next->directive_recognised();

  return 1;
}

//#define TEST_MOLECULAR_ASBTRACTION_SPECIFICATION
#ifdef TEST_MOLECULAR_ASBTRACTION_SPECIFICATION

int
main (int argc, char ** argv)
{
  if (1 == argc)
  {
    cerr << "Must specify argument(s)\n";
    return 4;
  }

  IWString s(argv[1]);

  cerr << "Building from '" << s << "'\n";

  Molecular_Abstraction_Directives_Node m;

  if (! m.build(s))
  {
    cerr << "Cannot parse '" << s << "'\n";
    return 23;
  }

  m.debug_print(cout);

  return 0;
}

#endif
