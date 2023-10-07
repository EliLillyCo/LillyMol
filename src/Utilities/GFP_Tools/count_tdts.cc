
#include "re2/re2.h"

#include "Foundational/data_source/iwstring_data_source.h"

#include "gfp.h"

using std::cerr;
using std::endl;

int
count_tdts_in_file(iwstring_data_source & input,
                   const IWString & identifier_tag)
{
  IWString tmp;
  tmp << '^' << identifier_tag;

  re2::StringPiece string_piece(tmp.data(), tmp.length());

  re2::RE2 pcn(string_piece);
  int rc = input.grep(pcn);

  if (0 == rc)
  {
    cerr << "No occurrences of " << pcn.pattern () << "' in input\n";
    return 0;
  }

  return rc;
}
