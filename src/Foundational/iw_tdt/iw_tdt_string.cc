#include <stdlib.h>

#include "iw_tdt.h"

int
IW_TDT::build (const char * s, int lens)
{
  return build (const_IWSubstring(s, lens));
}
int
IW_TDT::build (const const_IWSubstring & s)
{
  _zdata = s;

  int lens = s.length();

  _end.resize_keep_storage(0);

  lens--;   // we don't want the last newline character

  for (int i = 0; i < lens; i++)
  {
    if ('\n' == s[i])
      _end.add(i+1);
  }

  return _end.number_elements();
}
