#include <stdlib.h>

#include "iwstring.h"
#include "iwbits.h"
#include "iwhistogram.h"

int
IWHistogram::write_as_daylight_fingerprint (const IWString & tag,
                                            std::ostream & os) const
{
  const int * c = reinterpret_cast<const int *> (_count);

  IW_Bits_Base fp;

  fp.construct_from_array_of_ints (c, _nbuckets);

  return fp.write_daylight_ascii_representation (os, tag);
}

