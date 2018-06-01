#ifndef IW_STL_HASH_SET_H
#define  IW_STL_HASH_SET_H

#include <stdlib.h>
#include <iostream>

#include "iwconfig.h"

#if (GCC_VERSION >= 40405)
#include <unordered_set>
#define IW_Hash_Set std::unordered_set
#elif defined(__clang__)
#include <unordered_set>
#define IW_Hash_Set std::unordered_set
#elif (__GNUC__ >= 3)
#include <ext/hash_set>
using namespace __gnu_cxx;
#define IW_Hash_Set std::hash_set
#else
#include <hash_set>
using namespace stdext;
#define IW_Hash_Set hash_set
#endif

using std::cerr;
using std::endl;

#include "iwstring.h"

#include "iwhash.h"

class IWString_STL_Hash_Set : public IW_Hash_Set<IWString, IWStringHash>
{
  private:
  public:
    int contains (const IWString & t) const 
    {
      IW_Hash_Set<IWString, IWStringHash >::const_iterator f = find (t);

      return f != IW_Hash_Set<IWString, IWStringHash>::end ();
    }
};

typedef IWString_STL_Hash_Set IW_STL_Hash_Set;

#endif
