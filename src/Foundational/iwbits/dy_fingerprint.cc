#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>


    /**********************************************************************
    * This source code is public domain, and may be freely distributed,	  *
    * modified, and used for any purpose.  IT COMES WITH ABSOLUTELY NO    *
    * WARRANTY OF ANY KIND.						  *
    **********************************************************************/

/*===========================================================================
* CONVERTING DAYLIGHT'S ASCII FORMAT TO BINARY AND BACK
*
* Daylight uses a simple 6-bit-to-8-bit encoding to convert binary data
* into a printable ASCII form.  This for is used, for example, in Thor
* Datatrees (TDTs) to represent fingerprints.
*
* The following macros and array provide fast lookup to translate binary
* data into printable ascii, and vica versa.  Each 6 bits of binary (range
* 0-63) is converted to one of 64 characters in the set [,.0-9A-Za-z]; each
* 3-byte triplet thus converts to a 4-byte ASCII string.
*
* The table is a simple lookup table, but serves both as a "to" and "from"
* table.  The left column is index-to-ascii, the right column is
* ascii-to-binary.
*
* Every binary array is padded to a multiple of 3 bytes for the conversion;
* once the conversion is done you can't tell whether the last two bytes are
* pad bytes or real bytes containing zero.  To remedy this, an extra
* character is tacked on the ASCII representation; it will always be one of
* the characters '3', '2', or '1', indicating how many of the bytes in the
* last triplet are genuine.  That is, an ASCII-to-binary conversion will
* always produce an array whose length is a multiple of 3, but the last one
* or two bytes might just be pad bytes; the last ascii character indicates
* this.
===========================================================================*/


//#define ERROR(f,e) {fprintf(stderr, "%s: %e\n"); return nullptr;}

// Change definition, IAW

#define ERROR(f,e) {fprintf(stderr, "%s: %s\n", "DEATH", (e)); return nullptr;}


#define BIN2ASCII(s,c,d,e) \
{s[0] = lookup[ (c & 0xFC) >> 2].ascii;				\
 s[1] = lookup[((c & 0x03) << 4) | ((d & 0xF0) >> 4)].ascii;	\
 s[2] = lookup[((d & 0x0F) << 2) | ((e & 0xC0) >> 6)].ascii;	\
 s[3] = lookup[ (e & 0x3F)].ascii;				\
}

#define ASCII2BIN(s,c,d,e) \
{c = (lookup[s[0]].bin << 2) | ((lookup[s[1]].bin & 0x30) >> 4);	 \
 d = ((lookup[s[1]].bin & 0x0F) << 4) | ((lookup[s[2]].bin & 0x3C) >> 2);\
 e = ((lookup[s[2]].bin & 0x03) << 6) | lookup[s[3]].bin;		 \
}
  
static struct {
  char ascii, bin;
} lookup[] = {
/*  0   */ {'.', 0},
/*  1   */ {',', 0},
/*  2   */ {'0', 0},
/*  3   */ {'1', 0},
/*  4   */ {'2', 0},
/*  5   */ {'3', 0},
/*  6   */ {'4', 0},
/*  7   */ {'5', 0},
/*  8   */ {'6', 0},
/*  9   */ {'7', 0},
/* 10   */ {'8', 0},
/* 11   */ {'9', 0},
/* 12   */ {'A', 0},
/* 13   */ {'B', 0},
/* 14   */ {'C', 0},
/* 15   */ {'D', 0},
/* 16   */ {'E', 0},
/* 17   */ {'F', 0},
/* 18   */ {'G', 0},
/* 19   */ {'H', 0},
/* 20   */ {'I', 0},
/* 21   */ {'J', 0},
/* 22   */ {'K', 0},
/* 23   */ {'L', 0},
/* 24   */ {'M', 0},
/* 25   */ {'N', 0},
/* 26   */ {'O', 0},
/* 27   */ {'P', 0},
/* 28   */ {'Q', 0},
/* 29   */ {'R', 0},
/* 30   */ {'S', 0},
/* 31   */ {'T', 0},
/* 32   */ {'U', 0},
/* 33 ! */ {'V', 0},
/* 34 " */ {'W', 0},
/* 35 # */ {'X', 0},
/* 36 $ */ {'Y', 0},
/* 37 % */ {'Z', 0},
/* 38 & */ {'a', 0},
/* 39 ' */ {'b', 0},
/* 40 ( */ {'c', 0},
/* 41 ) */ {'d', 0},
/* 42 * */ {'e', 0},
/* 43 + */ {'f', 0},
/* 44 , */ {'g',  1},
/* 45 - */ {'h', 0},
/* 46 . */ {'i',  0},
/* 47 / */ {'j', 0},
/* 48 0 */ {'k',  2},
/* 49 1 */ {'l',  3},
/* 50 2 */ {'m',  4},
/* 51 3 */ {'n',  5},
/* 52 4 */ {'o',  6},
/* 53 5 */ {'p',  7},
/* 54 6 */ {'q',  8},
/* 55 7 */ {'r',  9},
/* 56 8 */ {'s', 10},
/* 57 9 */ {'t', 11},
/* 58 : */ {'u', 0},
/* 59 ; */ {'v', 0},
/* 60 < */ {'w', 0},
/* 61 = */ {'x', 0},
/* 62 > */ {'y', 0},
/* 63 ? */ {'z', 0},
/* 64 @ */ {'-', 0},
/* 65 A */ {'-', 12},
/* 66 B */ {'-', 13},
/* 67 C */ {'-', 14},
/* 68 D */ {'-', 15},
/* 69 E */ {'-', 16},
/* 70 F */ {'-', 17},
/* 71 G */ {'-', 18},
/* 72 H */ {'-', 19},
/* 73 I */ {'-', 20},
/* 74 J */ {'-', 21},
/* 75 K */ {'-', 22}, 
/* 76 L */ {'-', 23}, 
/* 77 M */ {'-', 24}, 
/* 78 N */ {'-', 25}, 
/* 79 O */ {'-', 26}, 
/* 80 P */ {'-', 27}, 
/* 81 Q */ {'-', 28}, 
/* 82 R */ {'-', 29}, 
/* 83 S */ {'-', 30},
/* 84 T */ {'-', 31},
/* 85 U */ {'-', 32},
/* 86 V */ {'-', 33},
/* 87 W */ {'-', 34},
/* 88 X */ {'-', 35},
/* 89 Y */ {'-', 36},
/* 90 Z */ {'-', 37},
/* 91 [ */ {'-', 0},
/* 92 \ */ {'-', 0},
/* 93 ] */ {'-', 0},
/* 94 ^ */ {'-', 0},
/* 95 _ */ {'-', 0},
/* 96 ` */ {'-', 0},
/* 97 a */ {'-', 38},
/* 98 b */ {'-', 39},
/* 99 c */ {'-', 40},
/*100 d */ {'-', 41},
/*101 e */ {'-', 42},
/*102 f */ {'-', 43},
/*103 g */ {'-', 44},
/*104 h */ {'-', 45},
/*105 i */ {'-', 46},
/*106 j */ {'-', 47},
/*107 k */ {'-', 48},
/*108 l */ {'-', 49},
/*109 m */ {'-', 50},
/*110 n */ {'-', 51},
/*111 o */ {'-', 52},
/*112 p */ {'-', 53},
/*113 q */ {'-', 54},
/*114 r */ {'-', 55},
/*115 s */ {'-', 56},
/*116 t */ {'-', 57},
/*117 u */ {'-', 58},
/*118 v */ {'-', 59},
/*119 w */ {'-', 60},
/*120 x */ {'-', 61},
/*121 y */ {'-', 62},
/*122 z */ {'-', 63},
/*123 { */ {'-', 0},
/*124 | */ {'-', 0},
/*125 } */ {'-', 0},
/*126 ~ */ {'-', 0},
/*127 */ {'-', 0}
};

/***************************************************************************
* FUNCTION: du_bin2ascii
*
* DESCRIPTION:
*	Converts binary to ASCII using the 6-bits-to-8 expansion.
*	Returns a newly-malloc'ed array containing the ASCII; the calling
*	program is responsible for deallocating the array.
*
*	On entry, the parameter blen contains the binary string's length;
*
***************************************************************************/

char *
du_bin2ascii(int *palen, int blen, const char *b)
{
  int ntriples, nleftover, i, j;
  char *ascii, *p, leftover[3];

  ntriples = blen / 3;
  nleftover = blen - (ntriples * 3);

  /**** Allocate the return array; use multiples of 4 bytes plus
        1 more for "pad" indication */
  *palen = ntriples * 4 + 1;
  if (nleftover != 0)
    *palen += 4;
  if (nullptr == (ascii = (char *) malloc(*palen + 1)))		/* 1 more for '\0' */
    ERROR("du_bin2ascii", "malloc failed\n");
  p = ascii;
  
  /**** convert every 3 bytes to 4 ascii chars ****/
  for (i = 0; i < ntriples * 3; i += 3) {
    BIN2ASCII(p,b[i],b[i+1],b[i+2]);
    p += 4;
  }

  /**** convert the trailing 1 or 2 bytes ****/
  if (nleftover > 0) {
    leftover[0] = leftover[1] = leftover[2] = 0;
    j = 0;
    while(i < blen)
      leftover[j++] = b[i++];
    BIN2ASCII(p,leftover[0],leftover[1],leftover[2]);
    p += 4;
  }

  /**** append with 3, 2, or 1 indicate # valid bytes in last triplet ****/
  if (nleftover == 0) nleftover = 3;
  *(p++) = '0' + nleftover;
  *p = '\0';

  return ascii;
}

/***************************************************************************
* FUNCTION: du_bin2ascii_len
*
* DESCRIPTION:
*	Given a binary array of data, returns the length the ascii string
*	would be, but without actually creating the ascii string.
***************************************************************************/

/*int du_bin2ascii_len(int blen, char *bin)
{
  int ntriples, nleftover, alen;

  ntriples = blen / 3;
  nleftover = blen - (ntriples * 3);

  alen = ntriples * 4 + 1;
  if (nleftover != 0)
    alen += 4;

  return alen;
}*/

/***************************************************************************
* FUNCTION: du_ascii2bin
*
* DESCRIPTION:
*	The converse of the above function: Converts each 4 bytes of
*	ASCII to 3 bytes of binary.  The last byte of ASCII will be either
*	0, 1, or 2, indicating how many of the last 3 bytes of binary
*	to chop off to get the correct length.
*
*	Returns a pointer to a newly-malloc'ed array, and *pblen contains
*	the number of bytes in that array.  The calling program is 
*	responsible for deallocating the returned array.
***************************************************************************/

/*char *
du_ascii2bin(int *pblen, int alen, char *ascii)
{
  char *binary, *p;
  int i;
  
  if ((alen % 4) != 1)
    ERROR("du_ascii2bin", "Invalid ASCII string (length wrong)");

  *pblen = (alen - 1) / 4;
  *pblen *= 3;
  if (nullptr == (binary = malloc(*pblen)))
    ERROR("du_ascii2bin", "malloc failed");

  p = ascii;
  for (i = 0; i < *pblen; i += 3) {
    ASCII2BIN(ascii, binary[i],  binary[i+1], binary[i+2]);
    ascii += 4;
  }

  i = *ascii - '0';
  *pblen -= 3 - i;

  return binary;
}*/

#include "assert.h"

int
du_ascii2bin(const char *ascii, const int nchars,
             unsigned char * binary, unsigned int & nbytes)
{
  assert (nullptr != binary);

  if ((nchars % 4) != 1)
    ERROR("du_ascii2bin", "Invalid ASCII string (length wrong)");

  /**** Compute the binary array's length ****/
  nbytes = (nchars - 1) / 4;
  nbytes *= 3;

  /**** Convert the ascii to binary ****/
//const char * p = ascii;
  for (unsigned int i = 0; i < nbytes; i += 3) {
    ASCII2BIN(ascii, binary[i],  binary[i+1], binary[i+2]);
    ascii += 4;
  }

  /**** Trim the pad zeros off the binary string ****/
  int i = *ascii - '0';
  nbytes -= 3 - i;

  return 1;
}
