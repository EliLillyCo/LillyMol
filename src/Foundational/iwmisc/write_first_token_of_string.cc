
#include "Foundational/iwstring/iwstring.h"

#include "misc.h"

std::ostream &
write_first_token_of_string(const IWString & zstring,
                            std::ostream & os)
{
  if (zstring.index(' ') < 0)
  {
    os << zstring;
    return os;
  }

  const_IWSubstring tmp(zstring);

  tmp.strip_leading_blanks();
  tmp.truncate_at_first(' ');

  os << tmp;

  return os;
}
