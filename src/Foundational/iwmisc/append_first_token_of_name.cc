#include <stdlib.h>

#include "Foundational/iwstring/iwstring.h"
#include "misc.h"

int
append_first_token_of_name(const IWString & mname,
                           IWString & output)
{
  int space_position = mname.index(' ');

  if (space_position < 0)
  {
    output << mname;
    return output.size();
  }

  if (space_position > 0)
    return output.strncat(mname, space_position);

// The unusual case of leading spaces in the name

  const_IWSubstring tmp(mname);

  tmp.remove_leading_chars(' ');

  if (0 == tmp.length())
  {
    std::cerr << "append_first_token_of_name:blank name!\n";
    output << "NONAME";
    return 0;
  }

  tmp.truncate_at_first(' ');

  output << tmp;

  return output.size();
}
