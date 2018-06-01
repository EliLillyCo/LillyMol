#ifndef IWZLIB_H
#define IWZLIB_H

/*
  I need a wrapper for zlib so we can do record oriented I/O
*/

class IWString;

// zlib not needed on MS Windows - John, what's the proper preprocessor symbol?

#ifdef NO_ZLIB

typedef off_t z_off_t;

class IW_ZLib_Wrapper
{
  private:
    void * _gzfile;

    char * _buffer;

    int _start_of_next_record;
    int _chars_in_buffer;

    char _record_delimiter;

  public:
    IW_ZLib_Wrapper () { abort();}
    ~IW_ZLib_Wrapper () {};

    int open_file (const char *) { return 0;}
    int close_file () { return 0;}

    int active () const { return 0;}

    int eof () const { return 0;}

    int read_bytes (void * destination, int len_destination) { return 0;};

    void set_record_delimiter (char s) { _record_delimiter = s;}

    size_t next_record (IWString &) { return 0;};

    int seekg (z_off_t) { return 0;}
    z_off_t tellg () const { return 0;}
};
#else
#include "zlib.h"

class IW_ZLib_Wrapper
{
  private:
    gzFile _gzfile;

    char * _buffer;

    int _start_of_next_record;
    int _chars_in_buffer;

    char _record_delimiter;

  public:
    IW_ZLib_Wrapper ();
    ~IW_ZLib_Wrapper ();

    int open_file (const char *);
    int close_file ();

    int active () const { return NULL != _gzfile;}

    int eof () const;

    int read_bytes (void * destination, int len_destination);

    void set_record_delimiter (char s) { _record_delimiter = s;}

    size_t next_record (IWString &);

    int seekg (z_off_t);
    z_off_t tellg () const;
};

#endif

#endif
