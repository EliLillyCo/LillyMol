#ifndef UNPACK_DATA_H
#define UNPACK_DATA_H

#include "Foundational/iwstring/iwstring.h"

class Unpack_Binary_Data
{
  private:
    IWString _unpack_format;

//  private functions

  public:
    int initialise (const char * s) { _unpack_format = s; return 1;}   // sometime make this do some checking

    int active () const { return _unpack_format.length();}

    template <typename T> int write_unpacked_data (const char * s, const int len, T & output) const;
};

#endif
