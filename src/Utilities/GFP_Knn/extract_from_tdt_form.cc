#include <stdlib.h>
#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#include "extract_from_tdt_form.h"

using std::cerr;
using std::endl;

int
extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       IWString & zresult)
{
  assert (tag.ends_with ('<'));

  if (! buffer.starts_with (tag))
  {
    cerr << "Tag mismatch. Expecting '" << tag << "' got '" << buffer << "'\n";
    return 0;
  }

  if (! buffer.ends_with ('>'))
  {
    cerr << "Malformed TDT item '" << buffer << "'\n";
    return 0;
  }

  zresult = buffer;
  zresult.remove_leading_chars (tag.length ());
  zresult.chop ();

  return 1;
}

template <typename T>
int
extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       T & zresult)
{
  IWString tmp;

  if (! extract_from_tdt_form (buffer, tag, tmp))
    return 0;

  if (! tmp.numeric_value (zresult))
  {
    cerr << "extract_from_tdt_form: invalid numeric '" << tmp << "'\n";
    return 0;
  }

  return 1;
}

#ifdef __GNUG__

template int extract_from_tdt_form (const_IWSubstring const &, IWString const &, int &);
template int extract_from_tdt_form (const_IWSubstring const &, IWString const &, float &);

#endif

// For some reason, the templates aren't being handled properly

int
extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       int & zresult)
{
  IWString tmp;

  if (! extract_from_tdt_form (buffer, tag, tmp))
    return 0;

  if (! tmp.numeric_value (zresult))
  {
    cerr << "extract_from_tdt_form: invalid numeric '" << tmp << "'\n";
    return 0;
  }

  return 1;
}

int
extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       float & zresult)
{
  IWString tmp;

  if (! extract_from_tdt_form (buffer, tag, tmp))
    return 0;

  if (! tmp.numeric_value (zresult))
  {
    cerr << "extract_from_tdt_form: invalid numeric '" << tmp << "'\n";
    return 0;
  }

  return 1;
}
