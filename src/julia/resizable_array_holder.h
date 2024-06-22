#ifndef JULIA_RESIZABLE_ARRAY_HOLDER_H_
#define JULIA_RESIZABLE_ARRAY_HOLDER_H_

#include "Foundational/iwaray/iwaray.h"

namespace lillymol_julia {
// Several times we need to use a resizable_array_p of something

template <typename T>
class ResizableArrayHolder {
  protected:
    const resizable_array_p<T>& _ref;

  public:
    ResizableArrayHolder(const resizable_array_p<T>& rhs) : _ref(rhs) {
    }

    uint32_t length() const {
      return _ref.size();
    }
    const T& operator[](int ndx) const {
      return *_ref[ndx];
    }
    const T& item(int ndx) const {
      return *_ref[ndx];
    }

    uint32_t size() const {
      return _ref.size();
    }
};

}  // lillymol_julia
#endif // JULIA_RESIZABLE_ARRAY_HOLDER_H_
