#ifndef IW_STL_HASH_MAP_H
#define  IW_STL_HASH_MAP_H

#include <iostream>
#include <stdlib.h>

#if (__GNUC__ >= 3)
#include <ext/hash_set>
using namespace __gnu_cxx;
#else
#include <hash_set>
using namespace stdext;
#endif

#include "iwstring.h"
#include "iwhash.h"

class IW_STL_Hash_Multiset_IWString : public hash_multiset<IWString, IWStringHash>
{
  private:
  public:
};

#endif
