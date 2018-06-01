#ifndef IWBITS_SUPPORT_H
#define IWBITS_SUPPORT_H
extern void
word_string_form (IWString & buffer, const unsigned int zword,
                  int bits_to_print,
                  const char t, const char f,
                  int include_space);
extern void
byte_string_form (IWString & buffer, const unsigned char zword,
                  int bits_to_print,
                  const char t, const char f,
                  int include_space);
#endif

