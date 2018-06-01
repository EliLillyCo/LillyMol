#include <stdlib.h>

#include "iwbits.h"

#include "precompbit8.h"

int
bits_in_common (const unsigned char * b1, const unsigned char * b2, int nb)
{
  int whole_words = nb / IW_BYTES_PER_WORD;

  int rc;

  if (whole_words)
  {
    const unsigned int * ib1 = (const unsigned int *) b1;
    const unsigned int * ib2 = (const unsigned int *) b2;

    rc = bits_in_common_word ((const unsigned int *) b1, (const unsigned int *) b2, whole_words);
  }
  else
    rc = 0;

  int extra_bytes = nb % IW_BITS_PER_WORD;

  if (extra_bytes)
  {
  }

  return rc;
} 

