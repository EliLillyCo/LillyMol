#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*#define CHECK_RESULTS*/

#if defined (CHECK_RESULTS)

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

/*
   Sun   METHOD1    sun compiler
   Linux METHOD4    Intel compiler

   Jun 2015. Some experiments using the Intel compiler with profile guided feedback
   time ./tbic  -n 10000 -b 2048 -d 0.2

   1 1.396
   2 1.462
   3 1.508
   4 1.471
   5 1.35
   6 1.637
   7 1.308
   1L 1.288

*/


// if we have sse3 use this

#define BIC_METHOD1L

//#define BIC_METHOD7
//#define BIC_METHOD5


#if defined (BIC_METHOD1)

#include "precompbit8.h"

#define IW_LAST_BYTE 0x000000ff

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  unsigned int c;
  const unsigned char * cp = (const unsigned char *)(&c);

  for (i = 0; i < nwords; i++)
  {
#if defined(CHECK_RESULTS)
    int tmp;
#endif

    c = *b1 & *b2;

#if defined(CHECK_RESULTS)
    tmp = eight_bit_count[c & IW_LAST_BYTE] + eight_bit_count[(c >> 16) & IW_LAST_BYTE] + eight_bit_count[(c >> 8) & IW_LAST_BYTE] + eight_bit_count[c >> 24];
    check_bits_in_common (*b1, *b2, tmp);
    rc += tmp;
#else
//  rc += eight_bit_count[c & IW_LAST_BYTE] + eight_bit_count[(c >> 16) & IW_LAST_BYTE] + eight_bit_count[(c >> 8) & IW_LAST_BYTE] + eight_bit_count[c >> 24];
//  rc += eight_bit_count[cp[0]] + eight_bit_count[cp[1]] + eight_bit_count[cp[2]] + eight_bit_count[cp[3]];
    rc += __builtin_popcount (c);
#endif

    b1++;
    b2++;
  }

  return rc;
}

#elif defined (BIC_METHOD1L)

static int
bits_in_common_32bit (const unsigned int * b1, const unsigned int * b2,
                      const int nwords)
{
  int i;
  int rc = 0;
  unsigned int c;

  for (i = 0; i < nwords; i++)
  {
    c = *b1 & *b2;

#if defined(CHECK_RESULTS)
    tmp = eight_bit_count[c & IW_LAST_BYTE] + eight_bit_count[(c >> 16) & IW_LAST_BYTE] + eight_bit_count[(c >> 8) & IW_LAST_BYTE] + eight_bit_count[c >> 24];
    check_bits_in_common (*b1, *b2, tmp);
    rc += tmp;
#else
    rc += __builtin_popcount (c);
#endif

    b1++;
    b2++;
  }

  return rc;
}

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;

#ifdef CHECK_RESULTS
  int j;
  int tmp;
#endif

  int nwords2 = nwords / 2;    // should check that it is even....

//if (2 * nwords2 != nwords)
//  return bits_in_common_32bit(b1, b2, nwords);

  const unsigned long * lb1 = (const unsigned long *) b1;
  const unsigned long * lb2 = (const unsigned long *) b2;

  for (i = 0; i < nwords2; i++)
  {
#if defined(CHECK_RESULTS)
    tmp = 0;
    for (j = 0; j < 64; ++j)
    {
      if (c & 0x1)
        tmp++;
      c >>= 1;
    }
    if (tmp != __builtin_popcountl(*lb1&*lb2))
    {
      fprintf(stderr, "mismatch on bic, compute %d, builtin %d\n", tmp, __builtin_popcountl(*lb1&*lb2));
      exit(1);
    }
    rc += tmp;
#else
//  rc += eight_bit_count[c & IW_LAST_BYTE] + eight_bit_count[(c >> 16) & IW_LAST_BYTE] + eight_bit_count[(c >> 8) & IW_LAST_BYTE] + eight_bit_count[c >> 24];
//  rc += eight_bit_count[cp[0]] + eight_bit_count[cp[1]] + eight_bit_count[cp[2]] + eight_bit_count[cp[3]];

    rc += __builtin_popcountl (*lb1&*lb2);  // gcc
//  rc += _mm_popcnt_u64(*lb1&*lb2);    // intel compiler
#endif

    lb1++;
    lb2++;
  }

  if (2 * nwords2 != nwords)
    rc += __builtin_popcount(*(b1 + nwords - 1) & *(b2 + nwords - 1));

  return rc;
}

#elif defined (BIC_METHOD2)

#include "precompbit.h"

static unsigned int andtmp[512];

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int nsw = nwords * 2;
  int rc = 0;
  const unsigned short * s = (const unsigned short *) &andtmp;

  for (i = 0; i < nwords; i++)
  {
    andtmp[i] = (b1[i]) & (b2[i]);
  }

#ifdef CHECK_AND
  for (i = 0; i < nwords; i++)
  {
    unsigned int j = b1[i] & b2[i];
    if (j != andtmp[i])
    {
      fprintf (stderr, "Vectorisation failure, i = %d, vector %u actually %u\n", i, andtmp[i], j);
    }
  }
#endif

  for (i = 0; i < nsw; i++)
  {
    rc += sixteen_bit_count[s[i]];
  }

  return rc;
}

#elif defined (BIC_METHOD3)

#include "precompbit.h"

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  int ns = nwords * 2;     /* 2 shorts per word */
  const unsigned short * sb1 = (const unsigned short *) b1;
  const unsigned short * sb2 = (const unsigned short *) b2;
  unsigned short c;

  for (i = 0; i < ns; i++)
  {
    c = *sb1 & *sb2;
    sb1++;
    sb2++;
    rc += sixteen_bit_count[c];
  }

  return rc;
}

#elif defined (BIC_METHOD4)

#include "precompbit.h"

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
#if defined (CHECK_RESULTS)
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

#if defined(CHECK_RESULTS)
    tmp = sixteen_bit_count[and_target.uc2.s0] + 
          sixteen_bit_count[and_target.uc2.s1];
    check_bits_in_common (*(b1 - 1), *(b2 - 1), tmp);
    rc += tmp;
#else
    rc += sixteen_bit_count[and_target.uc2.s0] + 
          sixteen_bit_count[and_target.uc2.s1];
#endif
  }

  return rc;
}

#elif defined (BIC_METHOD5)

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  unsigned int x;

  for (i = 0; i < nwords; i++)
  {
    x = b1[i] & b2[i];

    x -= (x >>1) & 0x55555555;
    x  = ((x >> 2) & 0x33333333) + (x & 0x33333333);
    x  = ((x >> 4) + x) & 0x0f0f0f0f0f;
    x *= 0x01010101;

    rc += x >> 24;
  }

  return rc;
}


#elif defined (BIC_METHOD6)

#include "precompbit8.h"

struct Uc4
{
  unsigned char c0;
  unsigned char c1;
  unsigned char c2;
  unsigned char c3;
};

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int j;
  int rc = 0;

  union And_Target {
    unsigned int uint;
    struct Uc4 uc4;
  } and_target;

  for (j = 0; j < nwords; j++)
  {
#if defined (CHECK_RESULTS)
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

#if defined(CHECK_RESULTS)
    tmp = eight_bit_count[and_target.uc4.c0] + 
          eight_bit_count[and_target.uc4.c1] +
          eight_bit_count[and_target.uc4.c2] +
          eight_bit_count[and_target.uc4.c3];
    check_bits_in_common (*(b1 - 1), *(b2 - 1), tmp);
    rc += tmp;
#else
    rc += eight_bit_count[and_target.uc4.c0] + 
          eight_bit_count[and_target.uc4.c1] +
          eight_bit_count[and_target.uc4.c2] +
          eight_bit_count[and_target.uc4.c3];
#endif
  }

  return rc;
}

#elif defined(BIC_METHOD7)

/*
  Taken from http://wm.ite.pl/articles/sse-popcount.html, author Wojciech Muła
*/

#include<stdint.h>

#ifdef ALIGN_DATA
#	define __aligned__ __attribute__((aligned(16)))
#else
#	define __aligned__
#endif

// lookup for SSE
static uint8_t POPCOUNT_4bit[16] __aligned__ = {
	/* 0 */ 0,
	/* 1 */ 1,
	/* 2 */ 1,
	/* 3 */ 2,
	/* 4 */ 1,
	/* 5 */ 2,
	/* 6 */ 2,
	/* 7 */ 3,
	/* 8 */ 1,
	/* 9 */ 2,
	/* a */ 2,
	/* b */ 3,
	/* c */ 2,
	/* d */ 3,
	/* e */ 3,
	/* f */ 4
};

// ---- SSSE3 - better alorithm, inner loop unrolled ----------------------
static int ssse3_popcount3(uint8_t* buffer, int chunks16) {
	static char MASK_4bit[16] = {0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf};

	uint32_t result;

	assert(chunks16 % 4 == 0);

	__asm__ volatile ("movdqu (%%eax), %%xmm7" : : "a" (POPCOUNT_4bit));
	__asm__ volatile ("movdqu (%%eax), %%xmm6" : : "a" (MASK_4bit));
	__asm__ volatile ("pxor    %%xmm5, %%xmm5" : : ); // xmm5 -- global accumulator

	result = 0;

	int k, n, i;

	i = 0;
	while (chunks16 > 0) {
		// max(POPCOUNT_8bit) = 8, thus byte-wise addition could be done
		// for floor(255/8) = 31 iterations
#define MAX (7*4)
		if (chunks16 > MAX) {
			k = MAX;
			chunks16 -= MAX;
		}
		else {
			k = chunks16;
			chunks16 = 0;
		}
#undef MAX
		__asm__ volatile ("pxor %xmm4, %xmm4"); // xmm4 -- local accumulator
		for (n=0; n < k; n+=4) {
#define body(index) \
			__asm__ volatile( \
				"movdqa	  (%%eax), %%xmm0	\n" \
				"movdqa    %%xmm0, %%xmm1	\n" \
				"psrlw         $4, %%xmm1	\n" \
				"pand      %%xmm6, %%xmm0	\n" \
				"pand      %%xmm6, %%xmm1	\n" \
				"movdqa    %%xmm7, %%xmm2	\n" \
				"movdqa    %%xmm7, %%xmm3	\n" \
				"pshufb    %%xmm0, %%xmm2	\n" \
				"pshufb    %%xmm1, %%xmm3	\n" \
				"paddb     %%xmm2, %%xmm4	\n" \
				"paddb     %%xmm3, %%xmm4	\n" \
				: : "a" (&buffer[index]));

			body(i);
			body(i + 1*16);
			body(i + 2*16);
			body(i + 3*16);
#undef body
			i += 4*16;
		}

		// update global accumulator (two 32-bits counters)
		__asm__ volatile (
			"pxor	%xmm0, %xmm0		\n"
			"psadbw	%xmm0, %xmm4		\n"
			"paddd	%xmm4, %xmm5		\n"
		);
	}

	// finally add together 32-bits counters stored in global accumulator
	__asm__ volatile (
		"movhlps   %%xmm5, %%xmm0	\n"
		"paddd     %%xmm5, %%xmm0	\n"
		"movd      %%xmm0, %%eax	\n"
		: "=a" (result)
	);

	return result;
}

int
bits_in_commonq (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  static unsigned int * tmp = NULL;
  static int ntmp = 0;

  int i;

  if (ntmp >= nwords)
    ;
  else if (NULL == tmp)
  {
    tmp = (unsigned int *) malloc(nwords * 4);
    ntmp = nwords;
  }
  else
  {
    free(tmp);
    tmp = (unsigned int *) malloc(nwords * 4);
    ntmp = nwords;
  }

  for (i = 0; i < nwords; ++i)
  {
    tmp[i] = b1[i] & b2[i];
  }

  return ssse3_popcount3((uint8_t *) tmp, nwords / 4);
}

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i, j;

  if (0 != nwords % 32)
  {
    fprintf(stderr, "Bad nwords %d\n", nwords);
    return -1;
  }

  static unsigned int * tmp = NULL;

  if (NULL == tmp)
    tmp = (unsigned int *) malloc(64);

#ifdef DEBUG_BIC7
  for (i = 0; i < 64; ++i)
  {
    tmp[i] = -1;
//  fprintf(stderr, "%u\n", tmp[i]);
  }
  fprintf(stderr, "%u sizeof(uint8_t) %lu\n", tmp[0], sizeof(uint8_t));
  for (i = 4; i <= 16; i += 4)
  {
    fprintf(stderr, "i = %d bic %d\n", i, ssse3_popcount3((uint8_t *) tmp, i));
  }
#endif

  int rc = 0;

//fprintf(stderr, "nwords %d\n", nwords);
  for (i = 0; i < nwords; i += 16)
  {
    for (j = 0; j < 16; ++j)
    {
      tmp[j] = b1[i+j] & b2[i+j];
    }
    rc += ssse3_popcount3((uint8_t *) tmp, 4);
//  fprintf(stderr, " i = %d, rc %d\n", i, rc);
  }

  return rc;
}

#else

#error "Must define a BIC_METHOD symbol"

#endif
