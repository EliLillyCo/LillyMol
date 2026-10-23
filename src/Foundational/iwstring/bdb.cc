#include "db_cxx.h"

#include "iwstring.h"

IWString&
operator<< (IWString& output, const Dbt& dbt) {
  output.strncat(reinterpret_cast<const char*>(dbt.get_data()),
                 dbt.get_size());
  return output;
}
