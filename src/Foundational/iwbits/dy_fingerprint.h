#ifndef IW_DAYLIGHT_FINGERPRINT_H
#define IW_DAYLIGHT_FINGERPRINT_H

extern int du_ascii2bin(const char *ascii, int nchars,
             unsigned char * binary, unsigned int & nbytes);

extern char * du_bin2ascii (int *palen, int blen, char *b);

#endif       
