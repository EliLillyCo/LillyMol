// AFile implementation functions.

#include "proto_support.h"

namespace iwmisc {

AFile::AFile(IWString& fname, int mode) {
  int flags = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  _fd = IW_FD_OPEN(fname.null_terminated_chars(), mode, flags);
}

AFile::~AFile() {
  IW_FD_CLOSE(_fd);
}

}  // namespace iwmisc
