#include "Foundational/iwstring/iwstring.h"

#include "tokens_in_quoted_string.h"

/*
  NOTE THIS IS NOT ROBUST.
  while it works for something like

  ||"<img src=""/img/s.gif"" alt=""Special Protocol"" height=""20"" width=""20"" border=""0""/>"||||

  it does not really check that properly...

  Make sure it works for something like


1,2,3,4,5,6,7,8,9
1,"2,x",3,4,5,6,7,8,9

  where some of the tokens are quoted
*/

int
tokens_in_quoted_string(const const_IWSubstring & buffer,
                        const char token_delimiter,
                        const int line_no)
{
  int rc = 0;

  const int n = buffer.length();

  if (0 == n)
    return 0;

  const char dquote = '"';

  int in_quote = buffer[0] == dquote;

  rc = 1;

  for (int i = 1; i < n; ++i)
  {
    if (dquote == buffer[i])
      in_quote = ! in_quote;
    else if (in_quote)
      ;
    else if (token_delimiter == buffer[i])
    {
//    if (dquote != buffer[i-1])
//      cerr << "Warning '" << buffer << "' column " << i << endl;
      rc++;
    }

//  cerr << " i = " << i << " rc " << rc << " char '" << buffer[i] << "'\n";
  }

  if (in_quote)
    std::cerr << "tokens_in_quoted_string:warning still quoted '" << buffer << "' line " << line_no << '\n';

  return rc;
}
