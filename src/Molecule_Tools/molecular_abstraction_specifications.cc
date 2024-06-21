#include <stdlib.h>
#include <iostream>

#include "molecular_abstraction_specifications.h"

using std::cerr;

#define MA_OPEN_PAREN '('
#define MA_CLOSE_PAREN ')'

void
DisplayAbstractionNames(std::ostream& output) {
  output << "Recognised abstraction functions\n";
  output << " allatoms    transformations applied to all atoms\n";
  output << " allbonds    transformations applied to all bonds\n";
  output << " bigring     remove all atoms except those in the largest ring\n";
  output << " arf         create an abstract ring form\n";
  output << " cbt         change bond types\n";
  output << " charge      placement, removal and changing formal charges\n";
  output << " comprconsec compress consecutive groups\n";
  output << " frag        fragment removal\n";
  output << " invscaf     invert scaffold\n";
  output << " isotope     placement, removal and changing isotopic labels\n";
  output << " rings       remove all non ring atoms, keep just the rings\n";
  output << " ringsys     generate ring systems\n";
  output << " rmat        remove specified atoms\n";
  output << " rmatoms     remove specified atoms\n";
  output << " rmbond      remove specified bonds\n";
  output << " rmrd2       remove two connected ring atoms - make rings smaller\n";
  output << " rmscaffold  remove the scaffold atoms\n";
  output << " rmspinach   same as scaffold\n";
  output << " rplink      replace a linker\n";
  output << " scaffold    reduce the molecule to scaffold atoms\n";
  output << " sss         perform substructure search\n";
  output << " translate   translate all atoms of a type to another element\n";
  return;
}

void
DisplayUsageExamples(std::ostream& output) {
  output << "Specify a set of operations that are performed in a pipeline\n";
  output << "allbonds(<1,2,3>)    change all bonds to type <1,2,3>\n";
  output << "arf(<lrs,lhc,ELE=,AROM=,ALIPH=>)\n";
  output << "bigring(spiro)\n";
  output << "cbt(<smarts>=<1,2,3>)  change the bond between first two matched atoms to <1,2,3>\n";
  output << "charge(<smarts, CHARGE=, SMARTS=, n=q, *>)\n";
  output << "compress(<>)\n";
  output << "fragment(KEEP=<smt> REMOVE=<smt>) fragment filtering\n";
  output << "isotope(<ISO=<i> SMARTS=<i=n> <i>)\n";
  output << "invscaf(SCH=<ele>)  keep the atoms that are NOT the scaffold\n";
  output << "rings()\n";
  output << "ringsys(<>)\n";
  output << "rmrd2(<e>)  remove [eRD2] atoms in ring\n";
  output << "rplink(ELE=<e>)\n";
  output << "scafold(keepfirst)  keepfirst means keep first ring attachment\n";
  output << "spinach(RMDBSC, AROM=<e>, ALIPH=<e>, CHAIN=<e>\n";
  output << "sss(SMARTS=<smt>, NONM, ...)   substructure search, ... is a query file, NONM pass non matches\n";
  output << "translate()\n";
  output << "\n";
  output << "Recognised by all operators...\n";
  output << " WRITE     write the current molecule\n";
  output << " FP        generate fingerprint for the current molecule\n";
  output << " WRITEC    append number of changes to output\n";
  output << " WRITEIF   write if the molecule is changed\n";
  output << " WRITE_MIN_ATOMS=<nn>   write if the molecule contains at least nn atoms\n";
  output << " WRITE_MAX_ATOMS=<nn>   write if the molecule contains at most  nn atoms\n";
  output << " WRITE_MIN_PARENT_ATOM_RATIO=<frac>   write if ratio of atoms now to parent > frac\n";
  output << " WRITE_MAX_PARENT_ATOM_RATIO=<frac>   write if ratio of atoms now to parent < frac\n";
  output << " ISO=<n>       apply atom typing\n";
  output << " AT=<atype>       apply atom typing\n";
  return;
}

Molecular_Abstraction_Directives_Node::Molecular_Abstraction_Directives_Node()
{
  _type = 0;

  _next = nullptr;

  return;
}

Molecular_Abstraction_Directives_Node::~Molecular_Abstraction_Directives_Node()
{
  if (nullptr != _next)
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
Molecular_Abstraction_Directives_Node::debug_print(std::ostream & os) const
{
  os << "Directive '" << _directive << "', args '" << _args << "'\n";
  if (nullptr != _next)
  {
    os << " next\n";
    return _next->debug_print(os);
  }

  return 1;
}

int
Molecular_Abstraction_Directives_Node::number_abstractions() const
{
  if (nullptr == _next)
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
  else if (MAD_GRAPH == _directive)
    _type = MAD_TYPE_GRAPH;
  else
  {
    cerr << "Molecular_Abstraction_Directives_Node::directive_recognised:unrecognised directive '" << _directive << "'\n";
    return 0;
  }

  if (nullptr != _next)
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
