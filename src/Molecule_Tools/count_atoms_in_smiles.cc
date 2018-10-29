#include <cctype>

#include "iwstring.h"

#include "smiles.h"

int
count_atoms_in_smiles(const const_IWSubstring & smiles)
{
  int rc = 0;

  int n = smiles.length();

//cerr << "Counting atoms in '" << smiles << "'\n";

  const char * s = smiles.rawchars();

  for (int i = 0; i < n; i++)
  {
    char c = s[i];

    if (isupper(c))
      rc++;
    else if ('l' == c)    // end of a Cl
      ;
    else if ('r' == c)    // end of a Br
      ;
    else if (islower(c))
      rc++;
    else if ('[' == c)
    {
      rc++;
      i++;
      while (']' != s[i])   // assume smiles is OK
      {
        i++;
      }
    }
    else if (' ' == c)
      break;
  }

//cerr << "Count " << rc << " atoms in '" << smiles << "'\n";

  return rc;
}


int
count_atoms_in_smiles(const const_IWSubstring & smiles,
                      int & nrings)
{
  int rc = 0;

  nrings = 0;

  int n = smiles.length();

  for (int i = 0; i < n; i++)
  {
    char c = smiles[i];

    if (isupper(c))
      rc++;
    else if ('l' == c)    // end of a Cl
      ;
    else if ('r' == c)    // end of a Br
      ;
    else if (islower(c))
      rc++;
    else if ('[' == c)
    {
      rc++;
      i++;
      while (']' != smiles[i])   // assume smiles is OK
      {
        i++;
      }
    }
    else if (isdigit(c))
      nrings++;
    else if (' ' == c)
      break;
  }

//cerr << "Count " << rc << " atoms in '" << smiles << "'\n";

  if (0 != nrings % 2)
  {
    cerr << "Ring estimation error '" << smiles << "' got " << nrings << '\n';
    nrings = 0;
  }
  else
    nrings = nrings / 2;

  return rc;
}

Smiles_Text::Smiles_Text()
{
  _ntokens = 0;
  _token = nullptr;
}

void
Smiles_Text::_free_array()
{
  if (nullptr != _token)
  {
    delete [] _token;
    _token = nullptr;
  }

  _ntokens = 0;

  return;
}

Smiles_Text::~Smiles_Text()
{
  _free_array();
}

int
Smiles_Text::_error(const int pos, const char * msg)
{
  cerr << msg << endl;
  for (int i = 0; i <= pos; ++i)
  {
    cerr << _smiles[i];
  }
  cerr << '\n';

  for (int i = 0; i < pos; ++i)
  {
    cerr << ' ';
  }

  cerr << "^\n";

  return 0;
}

template <typename T>
int
Smiles_Text::debug_print(T & output) const
{
  output << "Smiles_Text:debug_print:_ntokens " << _ntokens << '\n';

  for (int i = 0; i < _ntokens; ++i)
  {
    output << i << ' ';
    _token[i].debug_print(output);
  }

  return 1;
}

SmilesSmarts_Component::SmilesSmarts_Component()
{
  _isa = 0;
  _specific = 0;
  _atom_number = -1;
}

template <typename T>
int
SmilesSmarts_Component::debug_print(T & output) const
{
  output << "SmilesSmarts_Component::ISA " << _isa << ' ';
  if (SSC_ATOM == _isa)
    output << "atom " << _specific;
  else if (SSC_BOND == _isa)
  {
    if (SSC_SINGLE_BOND == _specific)
      output << "single bond";
    else if (SSC_DOUBLE_BOND == _specific)
      output << "double bond";
    else if (SSC_TRIPLE_BOND == _specific)
      output << "triple bond";
    else if (SSC_AROMATIC_BOND == _specific)
      output << "aromatic bond";
    else if (SSC_DIRECTIONAL_UP == _specific)
      output << "directional up bond";
    else if (SSC_DIRECTIONAL_DOWN == _specific)
      output << "directional down bond";
    else if (this->length() > 1)
      output << "smarts bond";
  }
  else if (SSC_OPAREN == _isa)
    output << "open paren";
  else if (SSC_CPAREN == _isa)
    output << "close paren";
  else if (SSC_DOT == _isa)
  {
    if (1 == this->length())
      output << "dot";
    else
      output << "dots";
  }
  else if (SSC_LEADING_NUMERIC == _isa)
    output << "leading numeric";
  else if (SSC_RING == _isa)
    output << "ring " << _specific;
  else if (0 == _isa)
    output << "NOT SET";
  else
  {
    cerr << "SmilesSmarts_Component::debug_print:unrecognised form " << _isa << endl;
    return 0;
  }

  const const_IWSubstring & s = (*this);
  output << ' ' << s << '\n';

  return 1;
}

int
SmilesSmarts_Component::set_ring(const char * s, const int nchars, const int r)
{
  const_IWSubstring & me = *this;
  me.set(s, nchars);

  _isa = SSC_RING;

  _specific = r;

  return 1;
}

int
SmilesSmarts_Component::set_element(const char * s, const int nchars, const int e)
{
  assert(e > 0 && e < 1000);

  _isa = SSC_ATOM;
  _specific = e;

  const_IWSubstring & me = *this;
  me.set(s, nchars);

  return 1;
}

int
SmilesSmarts_Component::set(const char * s, const int nchars, const int x)
{
  const_IWSubstring & me = *this;
  me.set(s, nchars);

  _isa = x;

  return 1;
}

int
SmilesSmarts_Component::set_bond(const char * s, const int nchars, const int x)
{
  _isa = SSC_BOND;

  _specific = x;

  const_IWSubstring & me = *this;
  me.set(s, nchars);

  return 1;
}

int
SmilesSmarts_Component::set_leading_numeric(const char * s, const int nchars)
{
  _isa = SSC_LEADING_NUMERIC;

  const_IWSubstring & me = *this;
  me.set(s, nchars);

  return 1;
}

int
Smiles_Text::_index_of_closing(const int x, const char copen, const char cclose) const
{
  const int nchars = _smiles.length();

  int level = 1;

  for (int i = (x+1); i < nchars; ++i)
  {
    if (_smiles[i] == copen)
      level++;
    else if (_smiles[i] == cclose)
    {
      level--;
      if (0 == level)
        return i;
    }
  }

  return -1;
}

int
Smiles_Text::build(const char * s,
                   const int nchars,
                   const int is_smiles)
{
  _smiles.set(s, nchars);

  _free_array();

  _token = new SmilesSmarts_Component[nchars];     // significant over-estimate

  _ntokens = 0;

  const char c0 = _smiles[0];

  int i = 0;

  if ('<' == c0 || '>' == c0 || isdigit(c0))
  {
    if (is_smiles)
      return _error(0, "not valid with smiles");
    if (! _consume_leading_numeric(i))
      return 0;
  }

  int paren_level = 0;

  while (i < nchars)
  {
    const char c = s[i];
//  cerr << " i = " << i << " '" << c << "'\n";

    if (']' == c)
      return _error(i, "mismatched close square bracket");

    if ('[' == c)
    {
      const int ccb = _index_of_closing(i, '[', ']');
      if (ccb < 0)
        return _error(i, "no closing square bracket");
      _token[_ntokens++].set(s + i, ccb - i, SSC_ATOM);
      i = ccb+1;
    }
    else if ('(' == c)
    {
      paren_level++;
      _token[_ntokens++].set(s + i++, 1, SSC_OPAREN);
    }
    else if (')' == c)
    {
      paren_level--;
      if (paren_level < 0)
        return _error(i, "mismatched parentheses");
      _token[_ntokens++].set(s + i++, 1, SSC_CPAREN);
    }
    else if ('%' == c)
    {
      if (! _consume_digits(i))
        return 0;
    }
    else if (isdigit(c))
      _consume_digit(i);
    else if ('.' == c)
    {
      if (! _handle_dot(i))
        return 0;
    }
    else if (_consume_element(i))
      ;
    else if (_consume_bond(i, is_smiles))
      ;
    else
      return _error(i, "unrecognised smiles character");
  }

  if (0 != paren_level)
    _error(i-1, "mismatched parentheses");

  return 1;
}

int
Smiles_Text::_consume_element(int & i)
{
  const char c = _smiles[i];
  char nextchar;
  if (i == _smiles.length() - 1)
    nextchar = ' ';
  else
    nextchar = _smiles[i+1];

  if ('C' == c && 'l' == nextchar)
  {
    _token[_ntokens++].set_element(_smiles.rawdata() + i, 2, 17);
    i += 2;
  }
  else if ('B' == c && 'r' == nextchar)
  {
    _token[_ntokens++].set_element(_smiles.rawdata() + i, 2, 35);
    i += 2;
  }
  else if ('C' == c || 'c' == c) 
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 6);
  else if ('O' == c || 'o' == c) 
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 8);
  else if ('N' == c || 'n' == c) 
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 7);
  else if ('P' == c || 'p' == c) 
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 15);
  else if ('S' == c || 's' == c) 
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 16);
  else if ('H' == c)
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 1);
  else if ('F' == c)
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 9);
  else if ('I' == c)
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 53);
  else if ('B' == c || 'b' == c) 
    _token[_ntokens++].set_element(_smiles.rawdata() + i++, 1, 5);
  else
    return 0;

  return 1;
}

/*
  Get the digits following a %
*/

int
Smiles_Text::_consume_digits(int & i)
{
  assert ('%' == _smiles[i]);

  const int chars_remaining = _smiles.length() - i;
  if (chars_remaining < 3)     //   %nn
    return _error(i, "Invalid ring specification");

  int digits = 0;
  int r = 0;

  for (int j = i + 1; j < _smiles.length(); ++j)
  {
    if (isdigit(_smiles[j]))
    {
      r = 10 * r + (_smiles[j] - '0');
      digits++;
    }
    else
      break;
  }

  if (0 == digits)
    return _error(i, "Invalid % ring specification");

  _token[_ntokens++].set_ring(_smiles.rawdata() + i, 1 + digits, r);

  i += 1 + digits;

  return 1;
}
int
Smiles_Text::_consume_digit(int & i)
{
  const int r = _smiles[i] - '0';

  assert (r >= 0 && r <= 9);

  _token[_ntokens++].set_ring(_smiles.rawdata() + i, 1, r);

  i++;

  return 1;
}

int
Smiles_Text::_handle_dot(int & cstart)
{
  assert ('.' == _smiles[cstart]);

  if (cstart == _smiles.length() - 1)
    return _error(cstart, "Smiles cannot end in .");

  if ('.' != _smiles[cstart+1])    // single dot
  {
    _token[_ntokens++].set(_smiles.rawdata() + cstart, 1, SSC_DOT);
    cstart++;
    return 1;
  }

// must be ... or ...{..}

  const int chars_remaining = _smiles.length() - cstart;
  if (chars_remaining < 4)    // cannot end in ...
    return _error(cstart, "Invalid . construct");

  if ('.' != _smiles[cstart+2])      // we checked +1 above
    return _error(cstart, "Invalid . construct");

  if ('{' != _smiles[cstart + 3])     // just three dots
  {
    _token[_ntokens++].set(_smiles.rawdata() + cstart, 3, SSC_DOT);
    cstart += 3;
    return 1;
  }

  const int close_brace = _index_of_closing(cstart + 3, '{', '}');
  if (close_brace < 0)
    return _error(cstart, "Mismatched {}\n");

  _token[_ntokens++].set(_smiles.rawdata() + cstart, close_brace - cstart, SSC_DOT);
  cstart = close_brace + 1;
  
  return 1;
}

int
Smiles_Text::_consume_bond(int & i, const int is_smiles)
{
  if (! is_smiles)
    return _consume_smarts_bond(i);

  const char c = _smiles[i];

  if ('-' == c)
   _token[_ntokens++].set_bond(_smiles.rawdata() + i, 1, SSC_SINGLE_BOND);
  else if ('=' == c)
   _token[_ntokens++].set_bond(_smiles.rawdata() + i, 1, SSC_DOUBLE_BOND);
  else if ('#' == c)
   _token[_ntokens++].set_bond(_smiles.rawdata() + i, 1, SSC_TRIPLE_BOND);
  else if ('/' == c)
   _token[_ntokens++].set_bond(_smiles.rawdata() + i, 1, SSC_DIRECTIONAL_UP);
  else if ('\\' == c)
   _token[_ntokens++].set_bond(_smiles.rawdata() + i, 1, SSC_DIRECTIONAL_DOWN);
  else if (':' == c)
   _token[_ntokens++].set_bond(_smiles.rawdata() + i, 1, SSC_AROMATIC_BOND);
  else    // not a bond
    return 0;

  i++;
  return 1;
}

int
Smiles_Text::_consume_smarts_bond(int & i)
{
  int consumed = 0;

  for ( ;i < _smiles.length(); ++i)
  {
    const char c = _smiles[i];

    if ('-' == c)
      consumed++;
    else if ('=' == c)
      consumed++;
    else if ('#' == c)
      consumed++;
    else if ('/' == c)
      consumed++;
    else if ('\\' == c)
      consumed++;
    else if (':' == c)
      consumed++;
    else if ('!' == c)
      consumed++;
    else if ('@' == c)
      consumed++;
    else if (',' == c)
      consumed++;
    else if ('&' == c)
      consumed++;
    else
      break;
  }

  if (0 == consumed)
    return 0;

  _token[_ntokens++].set_bond(_smiles.rawdata() + i - consumed, consumed, SSC_SMARTS_BOND);

  i += consumed;

  return 1;
}

int
Smiles_Text::_consume_leading_numeric(int & i)
{
  int consumed = 0;

  if ('<' == _smiles[i] || '>' == _smiles[i])
  {
    consumed++;
    i++;
  }

  if (i == _smiles.length())
    return _error(i, "incomplete");

  if (! isdigit(_smiles[i]))    // does not start with a numeric
  {
    if (consumed)
      return _error(i, "No numeric qualifier");

    return 1;
  }
  
  for ( ; i < _smiles.length(); ++i)
  {
    if (isdigit(_smiles[i]))
    {
      consumed++;
      i++;
    }
  }

  _token[_ntokens++].set_leading_numeric(_smiles.rawdata() + i - consumed, consumed);

  return 1;
}

int
Smiles_Text::establish_atom_numbers()
{
  int ndx = 0;

  for (int i = 0; i < _ntokens; ++i)
  {
    if (_token[i].is_atom())
    {
      _token[i].set_atom_number(ndx);
      ndx++;
    }
  }

  return ndx;
}

/*
  THE ORDER MUST BE THE SAME AS IN THE ENUM

  SSC_ATOM,
  SSC_BOND,
  SSC_RING,
  SSC_OPAREN,
  SSC_CPAREN,
  SSC_DOT,
*/

static int OK_ADJACENCY[6*6] = 
{
  1, 1, 1, 1, 1, 1,         // anything after an atom is OK
  1, 0, 1, 0, 0, 0,         // only atom or ring after bond
  1, 1, 1, 1, 1, 1,         // anything OK after ring
  1, 1, 1, 0, 0, 0,         // atom bond or ring after open paren
  1, 1, 1, 1, 0, 1,         // atom, bond ring, oparen, dot after close paran
  1, 0, 0, 1, 0, 0
};

int
Smiles_Text::ok_adjacency(const int issue_warning) const
{
  int istart = 1;

  if (SSC_LEADING_NUMERIC == _token[0].isa())
    istart = 2;

  int prev = _token[istart - 1].isa();

  for (int i = istart; i < _ntokens; ++i)
  {
    const int curr = _token[i].isa();
    if (! OK_ADJACENCY[6 * (prev-1) + (curr-1)])
    {
      if (issue_warning)
      {
        cerr << "Smiles_Text::ok_adjacency:invalid adjacency.\n";
        _token[i-1].debug_print(cerr);
        _token[i].debug_print(cerr);
      }

      return 0;
    }
  }

  return 1;
}

//#define TEST_SMILES_TEXT
#ifdef TEST_SMILES_TEXT

#include "cmdline.h"
#include "iwstring_data_source.h"

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "What does this programme do?\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
test_smiles_tokeniser_buffer(const const_IWSubstring & buffer,
                             const int is_smiles,
                             IWString_and_File_Descriptor & output)
{
  Smiles_Text stext;

  const_IWSubstring smiles, id;
  if (! buffer.split(smiles, ' ', id))
  {
    cerr << "Invalid form, not smiles and id '" << buffer << "'\n";
    return 0;
  }

  if (! stext.build(smiles.rawchars(), smiles.length(), is_smiles))
  {
    cerr << "Invalid input " << buffer << endl;
    return 0;
  }

  stext.establish_atom_numbers();

  if (! stext.ok_adjacency(1))
  {
    cerr << "Invalid adjacency\n";
    return 0;
  }

  output << buffer << '\n';
  stext.debug_print(output);

  return 1;
}

static int
test_smiles_tokeniser(iwstring_data_source & input,
                      const int is_smiles,
                      IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    buffer.strip_trailing_blanks();

    if (! test_smiles_tokeniser_buffer(buffer, is_smiles, output))
    {
      cerr << "Test failed\n";
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
test_smiles_tokeniser(const char * fname,
                      const int is_smiles,
                      IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_smiles_tokeniser(input, is_smiles, output);
}

static int
test_smiles_tokeniser(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  const int is_smiles = cl.option_present('s');

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! test_smiles_tokeniser(cl[i], is_smiles, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_smiles_tokeniser(argc, argv);

  return rc;
}
#endif
