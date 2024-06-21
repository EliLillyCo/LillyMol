#include <stdlib.h>
#include <ctype.h>

#include "Foundational/iwmisc/logical_expression.h"
#include "tokenise_atomic_smarts.h"

using std::cerr;

static constexpr const char* kDollarOpenParen = "$(";

Atomic_Smarts_Component::Atomic_Smarts_Component ()
{
  _unary_operator = 1;         // the default
  _op = IW_LOGEXP_UNDEFINED;   

  _next = nullptr;
      
  _from_primitive = 0;
}

Atomic_Smarts_Component::~Atomic_Smarts_Component()
{
  if (nullptr != _next)
    delete _next;

  return;
}

int
Atomic_Smarts_Component::ok() const
{
  if (_next != nullptr && _next->empty()) {
    return 0;
  }

  return 1;
}

static int
write_operator(std::ostream & os, int op)
{
  os << ' ';

  if (IW_LOGEXP_AND == op)
    os << '&';
  else if (IW_LOGEXP_OR == op)
    os << ',';
  else if (IW_LOGEXP_XOR == op)
    os << '^';
  else if (IW_LOGEXP_LOW_PRIORITY_AND == op)
    os << ';';
  else if (IW_LOGEXP_UNDEFINED == op)
    os << "undefined operator";
  else
    os << "What operator is this " << op;

  os << ' ';

  return os.good();
}

int
Atomic_Smarts_Component::debug_print(std::ostream & os) const
{
  if (! ok()) {
    cerr << "NOT OK\n";
  }

  write_operator(os, _op);

  if (0 == _unary_operator) {
    os << '!';
  }

  os << '\'';

  os.write(rawchars(), nchars());

  os << '\'';

  if (nullptr == _next) {
    os << '\n';
    return 1;
  }

  return _next->debug_print(os);
}

std::ostream &
operator << (std::ostream & os, const Atomic_Smarts_Component & rhs)
{
  write_operator(os, rhs.op());

  if (0 == rhs.unary_operator()) {
    os << '!';
  }

  os.write (rhs.rawchars(), rhs.nchars());

  if (nullptr == rhs.next()) {
    os << " _";
    return os;
  }

  return os << " (next)";
}

int
Atomic_Smarts_Component::ConvertToEnvironment() {
  if (IWString::starts_with(kDollarOpenParen)) {
    return 0;
  }

  IWString& me = *this;

  IWString new_value;
  new_value << kDollarOpenParen << me << ')';
  me = new_value;

  _from_primitive = 1;

  return 1;
}

static int
characters_in_environment(const_IWSubstring & smarts)
{
  if (smarts.nchars() <= 3) {
    cerr << "Atomic environment too short, must be at least '$(*)'\n";
    return 0;
  }

  assert (smarts.starts_with("$("));

  int paren_level = 1;
  for (int i = 2; i < smarts.nchars(); i++) {
/// cerr << "in env '" << smarts[i] << "'\n";
    if ('(' == smarts[i]) {
      paren_level++;
    } else if (')' == smarts[i]) {
      paren_level--;
      if (0 == paren_level){
        return i + 1;
      }
    }
  }

  cerr << "Mismatched parentheses in environment\n";
  return 0;
}

static int
number_repeated_characters(const const_IWSubstring & smarts,
                           char c)
{
  assert (c == smarts[0]);

  int rc = 1;

  for (int i = 1; i < smarts.length(); i++)
  {
    if (c != smarts[i])
      return rc;

    rc++;
  }

  return rc;
}

/*
  We are looking at something like 'R2' or '#35' and need to know how many
  numeric characters there are after 
*/

static int
number_numeric_characters (const const_IWSubstring & smarts,
                           int istart)
{
//cerr << "Counting numeric characters '" << smarts << "' starting at " << istart << '\n';

  int rc = 0;

  char c = smarts[istart];
  if ('<' == c)
  {
    rc++;
    istart++;
  }
  else if ('>' == c)
  {
    rc++;
    istart++;
  }

  for (int i = istart; i < smarts.length(); i++)
  {
    char c = smarts[i];
    if (c >= '0' && c <= '9')
      rc++;
    else
      break;
  }

//cerr << "Found " << rc << " numeric characters\n";
  return rc;
}

/*
  May 2002. The ! operator binds tightly. Identify the number of characters in the
  specification following the ! character
*/

static int
characters_in_next_primitive(const const_IWSubstring & smarts,
                             const int istart)
{
  // Just one character.
  if (istart == smarts.length() - 1) {
    return 1;
  }

  // "Cl" for example.
  if (isupper (smarts[istart]) && islower(smarts[istart + 1])) {
    return 2;
  }

// All the single character elements

  char c = smarts[istart];

  if ('a' == c || 'A' == c || 'B' == c || 'C' == c || 'F' == c || 'I' == c ||
      'N' == c || 'O' == c || 'P' == c || 'S' == c || 'U' == c || 'K' == c ||
      'V' == c || 'W' == c || 'c' == c || 'o' == c || 'n' == c || 's' == c || 'p' == c)
    return 1;

  int numeric_characters = number_numeric_characters(smarts, istart + 1);

//cerr << " contains " << numeric_characters << " numeric characters\n";

  if ('#' == c) {
    return 1 + numeric_characters;
  }

  if ('D' == c) {
    return 1 + numeric_characters;
  }

  if ('G' == c) {
    return 1 + numeric_characters;
  }

  if ('H' == c) {
    return 1 + numeric_characters;
  }

  if ('R' == c)
    return 1 + numeric_characters;

  if ('r' == c)
    return 1 + numeric_characters;

  if ('T' == c)
    return 1 + numeric_characters;

  if ('X' == c)
    return 1 + numeric_characters;

  if ('v' == c)
    return 1 + numeric_characters;

  if ('*' == c)     // what are they thinking    !*
    return 1;

  if (c >= '0' && c <= '9')
    return number_numeric_characters(smarts, 0);

  if ('+' == c)
  {
    if (numeric_characters)
      return 1 + numeric_characters;

    return number_repeated_characters(smarts, '+');
  }

  if ('-' == c)
  {
    if (numeric_characters)
      return 1 + numeric_characters;

    return number_repeated_characters(smarts, '-');
  }

  return 0;
}

static int
characters_in_next_token(const_IWSubstring & smarts)
{
  if ('$' == smarts[0]) {
    return characters_in_environment(smarts);
  }

  int square_bracket_level = 0;
  int curly_brace_level = 0;
  for (int i = 0; i < smarts.nchars(); i++)
  {
    char c = smarts[i];

//  cerr << "Examining '" << c << "' sqbrklvl = " << square_bracket_level << '\n';

    if ('[' == c)
    {
      square_bracket_level++;
      continue;
    }

    if (']' == c)
    {
      square_bracket_level--;
      if (0 == square_bracket_level)
        return i + 1;

      continue;
    }

    if ('{' == c)
    {
      curly_brace_level++;
      continue;
    }

    if ('}' == c)
    {
      curly_brace_level--;
      continue;
    }

    if (square_bracket_level) {
      continue;
    }

    if (',' == c && curly_brace_level) {
      continue;
    }

    if (',' == c || ';' == c || '!' == c || '^' == c || '$' == c || ':' == c) {
      return i;
    }
  }

  if (square_bracket_level) {
    cerr << "Mismatched square brackets\n";
    abort();
  }

  return smarts.nchars();
}

//#define DEBUG_PARSE

int
Atomic_Smarts_Component::parse(const_IWSubstring smarts)      // our own copy
{
//cerr << "Last char to process is '" << smarts[characters_to_process - 1] << "'\n";

#ifdef DEBUG_PARSE
  cerr << "Begin atomis smarts component parse '" << smarts << "'\n";
#endif
  
  int rc = _parse(smarts);

  if (0 == rc) {
    return 0;
  }

// If there is a mixture of primitive and $() tokens, conver them to the same form.
  int dollar_count = 0;
  int primitive_count = 0;
  for (const Atomic_Smarts_Component* asc = this; asc != nullptr; asc = asc->next()) {
    if (asc->starts_with(kDollarOpenParen)) {
      ++dollar_count;
    } else {
      ++primitive_count;
    }
  }

  // If either one is zero, this is not mixed.
  if (dollar_count == 0 || primitive_count == 0) {
    return rc;
  }

  // Mixed, convert all to environments.
  for (Atomic_Smarts_Component* asc = this; asc != nullptr; asc = asc->next()) {
    asc->ConvertToEnvironment();
  }

  // THis is no longer needed.

  return rc;
// Because of limitations of the implentation, we need to "fix" some special
// cases. Mainly we need to remove any operator which follows an atom and
// which preceeds an environment. For example, '[N;$(C-O)]'

  Atomic_Smarts_Component * a = this;
  while (a->_next)
  {
    if (a->starts_with("$("))    // we assume that all environments follow all the atoms specifications :-(
      return rc;

    if (a->_next->starts_with("$(")) {
      a->_op = IW_LOGEXP_UNDEFINED;
      return rc;
    }

    a = a->_next;
  }

  return rc;
}

int
Atomic_Smarts_Component::_parse(const_IWSubstring & smarts)
{
  assert (smarts.length());

#ifdef DEBUG_PARSE
  cerr << "Atomic smarts component parsing '" << smarts << "'\n";
#endif

  int characters_processed = 0;

  int length_of_our_token;

  // First extract any operator specification.

  if ('&' == smarts[0])
  {
    _op = IW_LOGEXP_AND;
    smarts++;
    characters_processed++;
  }
  else if (',' == smarts[0])
  {
    _op = IW_LOGEXP_OR;
    smarts++;
    characters_processed++;
  }
  else if ('^' == smarts[0])
  {
    _op = IW_LOGEXP_XOR;
    smarts++;
    characters_processed++;
  }
  else if (';' == smarts[0])
  {
    _op = IW_LOGEXP_LOW_PRIORITY_AND;
    smarts++;
    characters_processed++;
  }
  else           // no operator present
  {
    _op = IW_LOGEXP_UNDEFINED;
  }

#ifdef DEBUG_PARSE
  cerr << " op = " << _op << '\n';

  cerr << "Smarts now '" << smarts << "'\n";
#endif

  if (smarts.empty()) {
    cerr << "Atomic_Smarts_Component::_parse:just operator with no smarts\n";
    return 0;
  }

  if ('!' == smarts[0]) {
    _unary_operator = 0;

#ifdef DEBUG_PARSE
    cerr << "unary operator found in front of '" << smarts << "'\n";
#endif

    ++characters_processed;

    smarts++;
    if (smarts.empty()) {
      cerr << "Atomic_Smarts_Component::_parse: must have something to negate\n";
      return 0;
    }

// The ! operator binds tightly to a primitive    !#6H   gets treated differently from !$(CC)

    if (smarts.starts_with("$(")) {
      length_of_our_token = characters_in_next_token(smarts);
    } else {
      length_of_our_token = characters_in_next_primitive(smarts, 0);
      if (0 == length_of_our_token) {
        cerr << "Atomic_Smarts_Component::_parse: no recognised primitive following negation - contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
        return 0;
      }
    }
  } else {
    length_of_our_token = characters_in_next_token(smarts);
  }

#ifdef NOLONGERNEEDEEE
  cerr << "length_of_our_token " << length_of_our_token << " from '" << smarts << "'\n";
  if (smarts.starts_with("$(")) {
    length_of_our_token = characters_in_next_token(smarts);
  } else {
    length_of_our_token = characters_in_next_primitive(smarts, 0);
    if (0 == length_of_our_token) {
      cerr << "Atomic_Smarts_Component::_parse: no recognised primitive following negation - contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
      return 0;
    }
    cerr << "NEw determination of length_of_our_token " << length_of_our_token << '\n';
  }
#endif

  characters_processed = characters_processed + length_of_our_token;

  IWString::operator=(smarts);    // copy the smarts

  if (0 == length_of_our_token) {
    return 0;
  }

  iwtruncate(length_of_our_token);     // keep the number we used

#ifdef DEBUG_PARSE
  cerr << "Consumed " << length_of_our_token << " characters '";
  cerr.write(rawchars(), nchars());
  cerr << '\'';
  if (0 == _unary_operator)
    cerr << " unary op = " << _unary_operator;
  cerr << '\n';
#endif

  smarts += (length_of_our_token);      // get rid of the characters we consume

  characters_processed += length_of_our_token;

  // If empty, we are done.
  if (smarts.empty()) {
    return 1;
  }

  assert (nullptr == _next);

  _next = new Atomic_Smarts_Component;

  return _next->_parse(smarts);
}

#ifdef TEST_TKATSMARTS

static int verbose = 0;

static int parse_as_composite_query = 0;

static int parse_as_atomic_smarts = 0;

#include "cmdline.h"
#include "iwstring_data_source.h"

static int
test_composite_query (const_IWSubstring & buffer, std::ostream & output)
{
  assert (buffer.starts_with('[') && buffer.ends_with(']'));

  return 1;
}
  
static int
test_atomic_smarts (const_IWSubstring & buffer, std::ostream & output)
{
//assert (buffer.starts_with('[') && buffer.ends_with(']'));

  Atomic_Smarts_Component asc;

  if (! asc.parse(buffer))
  {
    cerr << "Cannot parse smarts\n";
    return 0;
  }

  asc.debug_print(output);

  return 1;
}

static int
test_tokenise_atomic_smarts (const_IWSubstring & buffer, std::ostream & output)
{
  if (parse_as_composite_query)
    return test_composite_query(buffer, output);
  else if (parse_as_atomic_smarts)
    return test_atomic_smarts(buffer, output);
}

static int
test_tokenise_atomic_smarts (iwstring_data_source & input, std::ostream & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! test_tokenise_atomic_smarts(buffer, output))
    {
      cerr << "Cannot parse '" << buffer << "', line " << input.lines_read() << '\n';
      return 0;
    }
  }

  return output.good();
}

static int
test_tokenise_atomic_smarts (const char * fname, std::ostream & output)
{
  if (::strlen(fname) > 2 && 'F' == fname[0] && ':' == fname[1])
  {
    fname += 2;
    iwstring_data_source input(fname);
    if (! input.ok())
    {
      cerr << "Cannot open '" << fname << "'\n";
      return 0;
    }

    return test_tokenise_atomic_smarts(input, output);
  }

  const_IWSubstring buffer = fname;

  return test_tokenise_atomic_smarts(buffer, output);
}

static void
usage (int rc)
{
  cerr << "Tester for smarts tokeniser\n";

  cerr << "  -m                  tokenise as composite query 'C#N&&CO'\n";
  cerr << "  -a                  tokenise as atomic smarts   '[CD2,R2;X2]'\n";

  cerr << "  -v                  verbose output\n";

  exit(rc);
}

static int
test_tokenise_atomic_smarts (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "mav");

  verbose = cl.option_count('v');

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  if (cl.number_elements().empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  parse_as_atomic_smarts = 1;

  if (cl.option_present('m'))
  {
    parse_as_composite_query = 1;
    if (verbose)
      cerr << "Will parse as a composite query\n";
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! test_tokenise_atomic_smarts(cl[i], cout))
      return i + 1;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_tokenise_atomic_smarts(argc, argv);

  return rc;
}
#endif

