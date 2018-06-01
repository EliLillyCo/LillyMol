#include <stdlib.h>
#ifdef _WIN32
  #include <winsock2.h>
#else
  #include <netinet/in.h>
#endif

#include "misc.h"

void
htonl_unsigned_long (void * z, unsigned int nw)
{
  uint32_t * l = reinterpret_cast<uint32_t *>(z);

  for (unsigned int i = 0; i < nw; i++)
  {
    l[i] = htonl(l[i]);
  }

  return;
}

void
htons_unsigned_short (void * z, unsigned int nw)
{
  uint16_t * s = reinterpret_cast<uint16_t *>(z);

  for (unsigned int i = 0; i < nw; i++)
  {
    s[i] = htons(s[i]);
  }

  return;
}

void
ntohl_unsigned_long (void * z, unsigned int nw)
{
  uint32_t * l = reinterpret_cast<uint32_t *>(z);

//cerr << "ntohl_unsigned_long " << nw << " words\n";
  for (unsigned int i = 0; i < nw; i++)
  {
//  cerr << "Convert " << l[i];
    l[i] = ntohl(l[i]);
//  cerr << " to " << l[i] << endl;
  }

  return;
}

void
ntohs_unsigned_short (void * z, unsigned int nw)
{
  uint16_t * s = reinterpret_cast<uint16_t *>(z);

  for (unsigned int i = 0; i < nw; i++)
  {
    s[i] = ntohs(s[i]);
  }

  return;
}
