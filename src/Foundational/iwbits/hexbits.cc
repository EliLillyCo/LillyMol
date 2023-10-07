#include <stdlib.h>

#include "iwbits.h"
#include "Foundational/iwstring/iwstring.h"

using std::cerr;
using std::endl;

/*
  Array is laid out so we can use hex characters as array indices
*/

static unsigned char hex2bit [] = {
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        0,    /* 0 */
        1,    /* 1 */
        2,    /* 2 */
        3,    /* 3 */
        4,    /* 4 */
        5,    /* 5 */
        6,    /* 6 */
        7,    /* 7 */
        8,    /* 8 */
        9,    /* 9 */
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        10,    /* a */
        11,    /* b */
        12,    /* c */
        13,    /* d */
        14,    /* e */
        15,    /* f */
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255,
        255};

int
IW_Bits_Base::construct_from_hex(const const_IWSubstring & zhex)
{
#ifdef ECHO_ARRAY
  for (int i = 0; i < 256; i++)
  {
    if (255 == hex2bit[i])
      continue;

    cerr << " i = " << i << ' ' << static_cast<char>(i) << " value " << static_cast<int>(hex2bit[i]) << endl;
  }
#endif

  int nchars = zhex.length();

  int nb = nchars * 4;     // each hex char is 4 bits

  if (0 == nb)     // maybe should do something else
    return 0;

  allocate_space_for_bits(nb);

  return HexToBinary(zhex, _bits);

#ifdef NOWITSOWNFUNCTION
  unsigned char * s = reinterpret_cast<unsigned char *>(_bits);

  for (int i = 0; i < nchars; i += 2)
  {
    unsigned char c = static_cast<unsigned char>(zhex[i]);

    unsigned char d1 = hex2bit[c];

    if (255 == d1)
    {
      cerr << "Invalid hex character '" << c << "'\n";
      return 0;
    }

//  cerr << "hex for " << c << " is " << static_cast<int>(d1) << endl;

    unsigned char d2;

    if (i + 1 < nchars)
    {
      c = zhex[i + 1];

      d2 = hex2bit[c];

      if (255 == d2)
      {
        cerr << "Invalid hex character '" << c << "'\n";
        return 0;
      }

//    cerr << "hex for " << c << " is " << static_cast<int>(d1) << endl;
    }
    else
      d2 = 0;

    d1 = (d1 << 4) | d2;

    *s = d1;
    s++;
  }

  return 1;
#endif
}

int
HexToBinary(const const_IWSubstring& zhex,
            unsigned char * destination) {
  const int nchars = zhex.length();
  for (int i = 0; i < nchars; i += 2) {
    unsigned char c = static_cast<unsigned char>(zhex[i]);

    unsigned char d1 = hex2bit[c];

    if (255 == d1)
    {
      cerr << "Invalid hex character '" << c << "'\n";
      return 0;
    }

//  cerr << "hex for " << c << " is " << static_cast<int>(d1) << endl;

    unsigned char d2;

    if (i + 1 < nchars)
    {
      c = zhex[i + 1];

      d2 = hex2bit[c];

      if (255 == d2)
      {
        cerr << "Invalid hex character '" << c << "'\n";
        return 0;
      }

//    cerr << "hex for " << c << " is " << static_cast<int>(d1) << endl;
    }
    else
      d2 = 0;

    d1 = (d1 << 4) | d2;

    *destination = d1;
    destination++;
  }

  return 1;
}

static const char bit2ascii [] = {
            '0',
            '1',
            '2',
            '3',
            '4',
            '5',
            '6',
            '7',
            '8',
            '9',
            'a',
            'b',
            'c',
            'd',
            'e',
            'f'};

namespace iwbits {
void
InternalHexForm(const unsigned char * s,
                int nbytes,
                IWString & destination)
{
  for (int i = 0; i < nbytes; i++)
  {
    const unsigned char ds = s[i];

//#define DEBUG_HEX_FORM
#ifdef DEBUG_HEX_FORM
    cerr << hex << static_cast<int>(ds) << " becomes " << hex << static_cast<int>((ds & 0xf0) >> 4) << " and " << hex << static_cast<int>(ds & 0x0f) << dec << endl;
#endif

    destination += bit2ascii[(ds & 0xf0) >> 4];
    destination += bit2ascii[ds & 0x0f];
  }

  return;
}

}  // namespace iwbits

int
IW_Bits_Base::hex_form(IWString & destination) const
{
  destination.resize_keep_storage(nbits() / 4);    // each hex character encodes 4 bits

  iwbits::InternalHexForm(_bits, _whole_bytes, destination);

  if (_extra_bits)
  {
    int extra_bytes = _extra_bits / IW_BITS_PER_BYTE;

    iwbits::InternalHexForm(_bits + _whole_bytes, extra_bytes, destination);
  }

  return ok();
}
