#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "precompbit.h"

/*#define CHECK_RESULTS*/
#ifdef CHECK_RESULTS
static const unsigned char one_bit_8[8] = {128, 64, 32, 16, 8, 4, 2, 1};

static const unsigned int one_bit_32[32] = {
  0x80000000,
  0x40000000,
  0x20000000,
  0x10000000,
  0x08000000,
  0x04000000,
  0x02000000,
  0x01000000,
  0x00800000,
  0x00400000,
  0x00200000,
  0x00100000,
  0x00080000,
  0x00040000,
  0x00020000,
  0x00010000,
  0x00008000,
  0x00004000,
  0x00002000,
  0x00001000,
  0x00000800,
  0x00000400,
  0x00000200,
  0x00000100,
  0x00000080,
  0x00000040,
  0x00000020,
  0x00000010,
  0x00000008,
  0x00000004,
  0x00000002,
  0x00000001 };


static int
check_bits_in_common (unsigned int b1, unsigned int b2,
                      int other_result)
{
  int bic = 0;
  unsigned int band = b1 & b2;
  int i;

  for (i = 0; i < 32; i++)
  {
    if (one_bit_32[i] & band)
      bic++;
  }

  if (bic == other_result)
    return 1;

  fprintf (stderr, "Warning, mismatch %u (%x) %u (%x) and %x, bic = %d, other = %d\n",
           b1, b1, b2, b2, band, bic, other_result);

  return 0;
}

#endif

#ifdef __sun

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  int ns = nwords * 2;     /* 2 shorts per word */
  const unsigned short * sb1 = (const unsigned short *) b1;
  const unsigned short * sb2 = (const unsigned short *) b2;

  for (i = 0; i < ns; i++)
  {
    unsigned short c = *sb1 & *sb2;
    sb1++;
    sb2++;
    rc += sixteen_bit_count[c];
  }

  return rc;
}

#else

struct Uc2
{
  unsigned short s0;
  unsigned short s1;
};

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int j;
  int rc = 0;

  union And_Target {
    unsigned int uint;
    struct Uc2 uc2;
  } and_target;

  for (j = 0; j < nwords; j++)
  {
#ifdef CHECK_RESULTS
    int tmp;
#endif

    and_target.uint = (*b1 & *b2);
    b1++;
    b2++;

#ifndef __i386__
    /* found this slowed things down on Intel */
    if (0 == and_target.uint)
      continue;
#endif

#ifdef CHECK_RESULTS
    tmp = sixteen_bit_count[and_target.uc2.s0] + 
          sixteen_bit_count[and_target.uc2.s1];
#elifdef __i386__
    ;
#else
    rc += sixteen_bit_count[and_target.uc2.s0] + 
          sixteen_bit_count[and_target.uc2.s1];
#endif

#ifdef CHECK_RESULTS
    check_bits_in_common (b1[j], b2[j], tmp);
    rc += tmp;
#endif
  }

  return rc;
}
#endif
