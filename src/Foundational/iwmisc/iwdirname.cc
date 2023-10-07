#include "Foundational/iwmisc/misc.h"

namespace iwmisc {

IWString
IWDirname(const IWString& fname) {
  IWString result(fname);
  int ndx = fname.rindex('/');
  if (ndx < 0) {
    return ".";
  }
  result.iwtruncate(ndx);
  return result;
}

}  // namespace iwmisc
