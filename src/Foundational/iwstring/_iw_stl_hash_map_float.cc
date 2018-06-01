/*
  Template instantiation stuff for hash maps
*/

#include <stdlib.h>

#include "iw_stl_hash_map.h"

#ifdef __GNUG__

//template class IW_STL_Hash_Map<IWString, float>;
//template class hash_map<IWString, float, IWStringHash >;
//template void hashtable<pair<IWString const, float>, IWString, IWStringHash, _Select1st<pair<IWString const, float> >, equal_to<IWString>, allocator<float> >::clear(void);
//template void hashtable<pair<IWString const, float>, IWString, IWStringHash, _Select1st<pair<IWString const, float> >, equal_to<IWString>, allocator<float> >::find_or_insert(pair<IWString const, float> const &);

#else   /* not __GNUG__ */

#endif  /* __GNUG__ */
