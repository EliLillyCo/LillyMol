/*
  Various supporting functions for sparse fingerprints.

  The reason for splitting these out into separate files is to avoid having
  all the hash_map stuff that is part of the Sparse_Fingerprint_Creator object
  dragged in if it is not needed
*/

#include <stdlib.h>
#ifdef _WIN32
  #include <winsock2.h>
  #pragma comment(lib, "Ws2_32.lib")
#else
  #include <netinet/in.h>
#endif

#include <iostream>
#include <memory>

#include "Foundational/iwbits/iwbits.h"

#include "misc.h"
#include "sparse_fp_creator.h"

using std::cerr;
using std::endl;

/*
  When writing out fingerprints that include counts, we write out groups of
  4 fingerprints as whole words, followed by a word containing the counts of
  the four previous fingerprints

  aaaabbbbccccdddd1234

  Note that when we get to the last word, we leave random values in the last bytes
  of the count word.
*/

int
words_needed_for_counted_form(int ndx)
{
  int words_needed = ndx + ndx / 4;
  if (0 != ndx % 4)
    words_needed++;

  return words_needed;
}

/*
  The final stage of forming ascii representations
*/

int
form_sparse_fingerprint(int ndx,
                         unsigned int * tmp,
                         IWString & dyascii)
{
  IW_Bits_Base fp;
  fp.construct_from_array_of_bits(reinterpret_cast<const unsigned char *>(tmp), ndx * IW_BITS_PER_WORD);

  return fp.daylight_ascii_representation(dyascii);
}

/*
  Sometimes we will have a fixed width fingerprint that just needs to be written
  We pack the array TMP
*/

int
pack_into_sparse_counted_form_v1(int nb,
                                  const int * k,
                                  unsigned int * tmp)
{
  union 
  {
    unsigned int counts;
    unsigned char c[IW_BYTES_PER_WORD];
  } counts;

  counts.counts = 0;

  int ndx = 0;
  int j = 0;          // index into the counts.c array

  for (int i = 0; i < nb; i++)
  {
    if (0 == k[i])
      continue;

    tmp[ndx] = htonl(i);
    ndx++;
    counts.c[j] = static_cast<unsigned char>(k[i]);
    j++;

    if (IW_BYTES_PER_WORD == j)
    {
      tmp[ndx] = counts.counts;
      ndx++;
      j = 0;
      counts.counts = 0;
    }
  }

// We didn't have 4 bits at the end

  if (j > 0)
  {
    tmp[ndx] = counts.counts;
    ndx++;
  }

#ifdef DEBUG_ENCODE_WITH_COUNTS
  cerr << "After building array, ndx = " << ndx << endl;
  for (int i = 0; i < ndx; i++)
  {
    cerr << tmp[i] << endl;
  }
#endif

  return 1;
}

int
non_colliding_counted_fingerprint_daylight_representation(int nb, const int * k, IWString & dyascii)
{
  int non_zero_bits = count_non_zero_occurrences_in_array(k, nb);

  if (0 == non_zero_bits)
  {
    cerr << "non_colliding_counted_fingerprint_daylight_representation: no bits set\n";
    dyascii.resize(0);

    return 1;
  }

  int words_needed = words_needed_for_counted_form(non_zero_bits);

//cerr << nb << " bits, but only " << non_zero_bits << " non zero ones, need " << words_needed << " words\n";
  
  unsigned int * tmp = new unsigned int[words_needed]; std::unique_ptr<unsigned int[]> free_tmp(tmp);

  pack_into_sparse_counted_form_v1(nb, k, tmp);

  IW_Bits_Base fp;

  fp.construct_from_array_of_bits(reinterpret_cast<const unsigned char *>(tmp), words_needed * IW_BITS_PER_WORD);

  return fp.daylight_ascii_representation(dyascii);
}
