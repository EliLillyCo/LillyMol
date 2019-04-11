#include <stdlib.h>

#if (__GNUC__ == 3)

#include <ext/hash_set>
using namespace __gnu_cxx;

#else

#include <hash_set>
using namespace std;

#endif

#include "iwstring.h"
#include "iw_stl_hash_set.h"

#ifdef __GNUG__

/*   gcc-2.95
template class hashtable<pair<IWString const, long>, IWString, IWStringHash, select1st<pair<IWString const, long> >, equal_to<IWString>, __default_alloc_template<false, 0> >;

template class hashtable<IWString, IWString, IWStringHash, _Identity<IWString>, equal_to<IWString>, allocator<IWString> >;
template class vector<_Hashtable_node<IWString> *, allocator<IWString> >;
template void fill<_Hashtable_node<IWString> **, _Hashtable_node<IWString> *>(_Hashtable_node<IWString> **, _Hashtable_node<IWString> **, _Hashtable_node<IWString> * const &);
template _Hashtable_node<IWString> ** fill_n<_Hashtable_node<IWString> **, unsigned int, _Hashtable_node<IWString> *>(_Hashtable_node<IWString> **, unsigned int, _Hashtable_node<IWString> * const &);
template class _Hashtable_const_iterator<IWString, IWString, IWStringHash, _Identity<IWString>, equal_to<IWString>, allocator<IWString> >;
*/

//template class std::hashtable<IWString, IWString, IWStringHash, std::_Identity<IWString>, std::equal_to<IWString>, std::allocator<IWString> >;
//template class std::vector<std::_Hashtable_node<IWString>*, std::allocator<IWString> >;

#else

#endif
