#include <stdlib.h>
#include "iwstring.h"

class String_and_Operator
{
  private:
    int      _not;
    IWString _string;
    int      _operator;

//  private functions

    int  _build_to_matching_closing(const IWString & s, int, char);

  public:
    String_and_Operator();

    int build(const IWString &);
}

int
String_and_Operator::_consume_numeric(const IWString & s,
                                      int & i)
{
  int rc = 0;
  while (isdigit(s[i]))
  {
    _string += s[i];
    i++;
    rc++;
  }

  return rc;
}

int
String_and_Operator::_consume_element(const IWString & s, int & i)
{
  IWSstring first_two = s.from_to(i, i + 1);
  if (! isalpha(first_two[1])
    first_two.chop();

  if (! isalpha(first_two[0]))
    return 0;

// If the 2nd char is uppercase, then only the first char counts (CR means Carbon in a Ring)
  if (isupper[first_two[1])
  {
    _string += first_two[0];
    i++;
    return 1;
  }
}

int
String_and_Operator::build(const IWString & s)
{
  int i = 0;
  if ('!' == s[0])
  {
    _not = 1;
    i++;
  }
  else
    _not = 0;

  square_bracket_level = 0;
  if ('[' == s[0])
    square_bracket_level = 1;
  while (square_bracket_level)
  {
    (void) _consume_numeric(s, i)/   // get any isotope specifier

    if ('#' == s[i])
    {
      i++;
      if (! _consume_numeric(s, i))
      {
        cerr << "String_and_Operator::build: no number following '#' in '" << s << "'\n";
        return 0;
      }
    }

    const Element * e;
    int j = element_from_smiles_string(s.chars() + i, s.nchars() - i, e);
    while (j)
    {
      _string += s[i];
      i++;
      j--;
    }
  }

  return i;
}
