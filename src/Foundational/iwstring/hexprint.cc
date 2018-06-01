#include "iwstring.h"

const static char * h = "0123456789abcdef";

int
// IWString::append_hex(const unsigned char * v, const int n)
IWString::append_hex(unsigned char const * v, const int n)
{
  make_room_for_extra_items(n+n);

  char * dest = _things + _number_elements;

#ifdef DEBUG_APPEND_HX
  for (int i = 0; i < n; ++i)
  {
    cerr << "append_hex writing i = " << i << ' ' << static_cast<int>(v[i]) << endl;
  }
#endif

  for (int i = 0; i < n; ++i)
  {
    int j = v[i] >> 4;
    *dest = h[j];
    dest++;

    j = v[i] & 0x0f;
    *dest = h[j];
    dest++;
  }

  _number_elements = dest - _things;

  return _number_elements;
}
