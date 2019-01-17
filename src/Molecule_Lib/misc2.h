#ifndef IW_MISC2_H
#define IW_MISC2_H

#include "iwmtypes.h"

extern void iwabort ();

extern int iw_rename (const char *, const char *);

extern int iw_getpid ();

extern int int_comparitor_larger (const int *, const int *);

extern int uint64_comparitor_smaller (const iw_uint64_t *, const iw_uint64_t *);

extern void iwxor (const int *, int *, int);

class const_IWSubstring;

extern int fetch_numeric (const const_IWSubstring & string, int & value, int max_chars = 0);
extern int fetch_numeric_char (const char * string, int & value, int max_chars);

/*
  Sometimes we need to compute combinatorial permutations and we
  may be dealing with numbers larger than can be held in an int
*/

extern iw_uint64_t iw_combinatorial_combinations (int n, int k);

template <typename T> int skip_to_string (T & input, const char * target, int report_discard);

// identify the + characters in a reaction smiles. Complicated by the presence of + signs inside square brackets

extern int identify_plus_positions (const const_IWSubstring & buffer, resizable_array<int> & pos);
extern int splitOnPlusses(const const_IWSubstring & buffer,
                        resizable_array_p<const_IWSubstring> & parts);

#endif
