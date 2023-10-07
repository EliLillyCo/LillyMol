#include <stdlib.h>

#include "Foundational/iwmisc/misc.h"

#include "iwbits.h"

// I didn't finish this because I could not figure out what to do if
// there are extra bytes. Our bit vector may have 3 bytes allocated,
// but swapping those would reuquire 4 bytes. Revisit this if it ever
// matters

#ifdef NOT_IMPLEMENTED
void
IW_Bits_Base::swap_byte_order()
{
  if (0 == _nbits)
    return;

  int whole_words = _whole_bytes / IW_BYTES_PER_WORD;

  if (whole_words)
    rick_higgs_byte_swap(whole_words, _bit);

  int extra_bytes = _whole_bytes % IW_BYTES_PER_WORD;
  if (_extra_bits)
    extra_bytes++;

  if (extra_bytes)

  if (0 == _extra_bits)
  {

  return;
}
#endif
