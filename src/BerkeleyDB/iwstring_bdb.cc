// Support functions for IWString/BerkeleyDB
// Note that the operator itself is declared in iwstring.h, we
// don't bother with a separate header file.

#include "db_cxx.h"

#include "Foundational/iwstring/iwstring.h"

IWString&
operator<< (IWString& output, const Dbt& dbt) {
  output.strncat(reinterpret_cast<const char*>(dbt.get_data()),
                 dbt.get_size());
  return output;
}
