#ifndef IW_STL_HASH_MAP_H
#define  IW_STL_HASH_MAP_H

#include <stdlib.h>

#if (GCC_VERSION >= 40405)
#include <unordered_map>
#define IW_Hash_Multimap std::unordered_multimap
#elif (__GNUC__ >= 3)
#include <ext/hash_map>
using namespace __gnu_cxx;
#define IW_Hash_Multimap hash_map
#else
#include <hash_map>
#define IW_Hash_Multimap hash_map
#endif

#include "iwstring.h"
#include "iwhash.h"

template <typename K, typename V>
class IW_STL_Hash_Multimap : public IW_Hash_Multimap<K, V, IWStringHash>
{
  private:
  public:
};

typedef IW_STL_Hash_Multimap<IWString, int> IW_STL_Hash_Multimap_int;
typedef IW_STL_Hash_Multimap<IWString, float> IW_STL_Hash_Multimap_float;
typedef IW_STL_Hash_Multimap<IWString, double> IW_STL_Hash_Multimap_double;
typedef IW_STL_Hash_Multimap<IWString, IWString> IW_STL_Hash_Multimap_IWString;

#endif
