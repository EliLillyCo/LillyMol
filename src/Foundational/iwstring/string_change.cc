#include <stdlib.h>

#include "iwstring.h"

int
IWString::change (int istart,
                  int istop,
                  const char * new_data,
                  int len_new_data)
{
  assert (istart <= istop);
  assert (istart >= 0);
  assert (istop < _number_elements);

  if (0 == len_new_data)
    return remove_from_to (istart, istop);

  int size_difference = len_new_data - (istop - istart + 1);

#ifdef DEBUG_CHANGE
  cerr << (*this) << " change from " << istart << " to " << istop << " diff " << size_difference << endl;
#endif

  if (size_difference < 0)    // shrinking
  {
    ::memmove (_things + istart + len_new_data, _things + istop + 1, _number_elements - istop);
  }
  else if (size_difference > 0)    // growing
  {
#ifdef DEBUG_CHANGE
    cerr << "Increasing " << _number_elements << " diff " << size_difference << " alloc " << _elements_allocated << endl;
#endif
    if (_number_elements + size_difference > _elements_allocated)
      make_room_for_extra_items (size_difference);


    ::memmove (_things + istart + len_new_data, _things + istop + 1, _number_elements - istop - 1);
  }

  ::memcpy (_things + istart, new_data, len_new_data);

  _number_elements += size_difference;

#ifdef DEBUG_CHANGE
  cerr << "N set to " << _number_elements << endl;
#endif

  assert (ok ());

  return 1;
}
