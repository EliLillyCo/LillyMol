#ifndef ADD_DATAITEM_H
#define ADD_DATAITEM_H

#include "iw_tdt.h"

template <typename S>
int
IW_TDT::add_dataitem (const char * tag,
                      int len_tag,
                      const S & s,
                      int where_to_put)
{
  IWString to_insert;
  to_insert.strncat (tag, len_tag);
  to_insert << s;
  to_insert << '>';
  if (include_newlines_in_tdt ())
    to_insert.add ('\n');

  if (where_to_put < 0)
    where_to_put = _end.number_elements () - where_to_put;
   
  assert (where_to_put >= 0);

  if (where_to_put > 0)
    where_to_put = _end[where_to_put];

  return _zdata.insert (to_insert, where_to_put);
}

#endif
