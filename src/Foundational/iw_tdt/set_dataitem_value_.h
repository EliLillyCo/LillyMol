#include <stdlib.h>

#include "iw_tdt.h"

template <typename S>
int
IW_TDT::set_dataitem_value (const char * tag,
                            int len_tag,
                            const S & new_data,
                            int which_one)
{
  int i = _find_index_in_end_array (tag, len_tag, which_one);
  if (i < 0)
  {
    cerr << "IW_TDT::set_dataitem_value:no " << which_one << " instance of " << cerr.write (tag, len_tag) << "' found\n";
    return 0;
  }

  int e = _end[i];   // end of the current data

  int s;             // start of current data
  if (0 == i)
    s = 0;
  else
    s = _end[i - 1];

  return _zdata.change (s, e, new_data);
}
