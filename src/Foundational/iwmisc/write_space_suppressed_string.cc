#include "misc.h"
#include "iwstring.h"

std::ostream &
write_space_suppressed_string (const IWString & zstring,
                               std::ostream & os,
                               char fill_char)
{
  if (zstring.contains (' '))
  {
    IWString tmp (zstring);
    tmp.gsub (' ', fill_char);

    os << tmp;
  }
  else
    os << zstring;

  return os;
}

