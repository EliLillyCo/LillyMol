#ifndef IWUUENCODE_H
#define IWUUENCODE_H

class IWString;
class const_IWSubstring;

extern int IWuuencode_append (const void * p,
                   int nchars,
                   IWString & destination);

/*
  If you have an encoded form, you need to know how many bytes are needed for decoding
*/

extern int IWuudecode_bytes_needed (int);  // length of string holding encoded form

extern int IWuudecode (const unsigned char * encoded, int nchars, unsigned char * destination);
extern int IWuudecode (const IWString & encoded, unsigned char * destination);
extern int IWuudecode (const const_IWSubstring & encoded, unsigned char * destination);

#endif
