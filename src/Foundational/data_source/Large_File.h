#ifndef IW_LARGE_FILE_H
#define IW_LARGE_FILE_H

#include "iwstring.h"

class IW_Large_File
{
  private:
    int _fd;

//  We have an array of characters we have read from the file.


    char _pending[4096];

//  We also need an index into the _pending array for where to start

    int  _pstart;

//  as well as the location of the last character in _pending

    int _pend;

    int _records_read;

//  We need to keep track of our state if the constructor with a file name is used

    int _ok;

//  private functions

    int _index_of_newline (char newline) const;

    void _default_values ();

  public:
    IW_Large_File ();
    IW_Large_File (const char * fname);
    ~IW_Large_File ();

    int open (const char * fname);
    int is_open () const { return _fd >= 0;}

    int ok   () const { return _ok;}
    int good () const { return _ok;}

    int next_record (IWString &, char = '\n');

    int records_read () const { return _records_read; }
};

#endif
