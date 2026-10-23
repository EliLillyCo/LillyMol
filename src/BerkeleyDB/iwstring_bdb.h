#ifndef BERKELEYDB_IWSTRING_BDB_H_
#define BERKELEYDB_IWSTRING_BDB_H_

#include "db_cxx.h"

#include "Foundational/iwstring/iwstring.h"

extern IWString & operator << (IWString & output, const Dbt& dbt);

#endif // BERKELEYDB_IWSTRING_BDB_H_

