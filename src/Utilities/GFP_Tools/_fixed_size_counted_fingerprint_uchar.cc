#include <stdlib.h>

/*
  We only need to include iwbits so we get the definition for similarity_type_t
*/

#include "Foundational/iwbits/iwbits.h"

#define FIXED_SIZE_COUNTED_FINGERPRINT_IMPLEMENTATION 
#include "fixed_size_counted_fingerprint.h"

template class Fixed_Size_Counted_Fingerprint_Base<unsigned char>;
