#include <stdlib.h>

#include "iwstring.h"

#include "iwuuencode.h"

/*
  I want something to create an ASCII representation of binary data.

  Take code from uuencode/decode

  We always use the base64 encoding
*/

const char uu_base64[64] =
{
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
  'w', 'x', 'y', 'z', '0', '1', '2', '3',
  '4', '5', '6', '7', '8', '9', '+', '/'
};

/* Pointer to the translation table we currently use.  */
const char *trans_ptr = uu_base64;

/* ENC is the basic 1 character encoding function to make a char printing.  */
#define ENC(Char) (trans_ptr[(Char) & 077])


int
IWuuencode_append (const void * v,
                   int n,
                   IWString & destination)
{
  const char * p = reinterpret_cast<const char *> (v);

  char ch;

  for (; n > 2; n -= 3, p += 3)
  {
    ch = *p >> 2;
    ch = ENC (ch);
    destination.add (ch);
    ch = ((*p << 4) & 060) | ((p[1] >> 4) & 017);
    ch = ENC (ch);
    destination.add (ch);
    ch = ((p[1] << 2) & 074) | ((p[2] >> 6) & 03);
    ch = ENC (ch);
    destination.add (ch);
    ch = p[2] & 077;
    ch = ENC (ch);
    destination.add (ch);
  }

  if (0 == n)
    return 1;

  char c1 = *p;
  char c2 = n == 1 ? 0 : p[1];

  ch = c1 >> 2;
  ch = ENC (ch);
  destination.add (ch);

  ch = ((c1 << 4) & 060) | ((c2 >> 4) & 017);
  ch = ENC (ch);
  destination.add (ch);

  if (2 == n)
  {
    ch = (c2 << 2) & 074;
    ch = ENC (ch);
    destination.add (ch);
  }

  return 1;
}
