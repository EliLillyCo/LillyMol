#include <stdlib.h>

#include "iwstring_data_source.h"
//#include "iwcrex.h"

#include "gfp.h"

int
count_tdts_in_file (iwstring_data_source & input,
                    const IWString & identifier_tag)
{
  IWString tmp;
  tmp << '^' << identifier_tag;

  IW_Regular_Expression pcn (tmp);
  int rc = input.grep (pcn);

  if (0 == rc)
  {
    cerr << "No occurrences of " << pcn.source () << "' in input\n";
    return 0;
  }

  return rc;
}
