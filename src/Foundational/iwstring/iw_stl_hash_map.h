#ifndef IW_STL_HASH_MAP_H
#define  IW_STL_HASH_MAP_H

#include <stdlib.h>
#include <iostream>

#include "Foundational/iwmisc/iwconfig.h"

#if (GCC_VERSION >= 40405)
#include <unordered_map>
#define IW_Hash_Map std::unordered_map
#elif defined(__clang__)
#include <unordered_map>
#define IW_Hash_Map std::unordered_map
#elif (__GNUC__ >= 3)
#include <ext/hash_map>
using namespace __gnu_cxx;
#define IW_Hash_Map hash_map
#else
#include <hash_map>
#define IW_Hash_Map hash_map
#endif

#include "iwstring.h"
#include "iwhash.h"

template <typename T, typename V>
class IW_STL_Hash_Map : public IW_Hash_Map<T, V, IWStringHash>
{
  private:
  public:
    int contains (const T & t) const 
    {
      typename IW_Hash_Map<T, V, IWStringHash>::const_iterator f = IW_Hash_Map<T, V, IWStringHash>::find (t);

      return f != IW_Hash_Map<T, V, IWStringHash>::end ();
    }

    typedef typename IW_STL_Hash_Map<T, V>::const_iterator const_iterator;
};

typedef IW_STL_Hash_Map<IWString, int> IW_STL_Hash_Map_int;
typedef IW_STL_Hash_Map<IWString, unsigned int> IW_STL_Hash_Map_uint;
typedef IW_STL_Hash_Map<IWString, float> IW_STL_Hash_Map_float;
typedef IW_STL_Hash_Map<IWString, double> IW_STL_Hash_Map_double;
typedef IW_STL_Hash_Map<IWString, long> IW_STL_Hash_Map_long;
typedef IW_STL_Hash_Map<IWString, off_t> IW_STL_Hash_Map_off_t;
typedef IW_STL_Hash_Map<IWString, IWString> IW_STL_Hash_Map_String;

#endif
