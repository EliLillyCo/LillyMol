#ifndef FOUNDATIONAL_IWSTRING_IW_STL_HASH_SET_H_
#define FOUNDATIONAL_IWSTRING_IW_STL_HASH_SET_H_

#include <initializer_list>

#include "Foundational/iwmisc/iwconfig.h"

#include <unordered_set>
#define IW_Hash_Set std::unordered_set

#include "iwstring.h"

#include "iwhash.h"

class IWString_STL_Hash_Set : public IW_Hash_Set<IWString, IWStringHash>
{
  private:
  public:
    IWString_STL_Hash_Set() {
    };
    IWString_STL_Hash_Set(std::initializer_list<IWString> ilist) {
      for (const auto& item : ilist) {
        insert(item);
      }
    }

    int contains (const IWString & t) const 
    {
      IW_Hash_Set<IWString, IWStringHash >::const_iterator f = find (t);

      return f != IW_Hash_Set<IWString, IWStringHash>::end ();
    }
};

typedef IWString_STL_Hash_Set IW_STL_Hash_Set;

#endif  // FOUNDATIONAL_IWSTRING_IW_STL_HASH_SET_H_
