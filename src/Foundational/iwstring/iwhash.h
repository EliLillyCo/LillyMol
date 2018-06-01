#ifndef IW_HASH_H
#define IW_HASH_H

#include "iwconfig.h"

class IWStringHash
{
  private:
  public:

#if defined (IW_INTEL_COMPILER)

    static const size_t bucket_size = 4;
    static const size_t min_buckets = 8;
    bool  operator () (const IWString &, const IWString &) const;

#endif

    size_t operator () (const IWString &) const;
};

#endif
