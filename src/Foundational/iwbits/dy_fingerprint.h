#ifndef FOUNDATIONAL_IWBITS_DY_FINGERPRINT_H_
#define FOUNDATIONAL_IWBITS_DY_FINGERPRINT_H_

extern int du_ascii2bin(const char *ascii, int nchars,
             unsigned char * binary, unsigned int & nbytes);

extern char * du_bin2ascii (int *palen, int blen, char *b);

#endif  // FOUNDATIONAL_IWBITS_DY_FINGERPRINT_H_
