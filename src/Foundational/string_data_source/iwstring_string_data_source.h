#ifndef STRING_DATA_SOURCE_H
#define STRING_DATA_SOURCE_H

#include "iwstring.h"
#include "iwcrex.h"

/*
  We want to read data from a string rather than a file descriptor
*/

class String_Data_Source
{
  private:
    const char * _src;

    int _iptr;

    int _lines_read;

    int _strip_trailing_blanks;

    int _dos;

//  private functions

    void _default_values();

  public:
    String_Data_Source(const char *);
    String_Data_Source(const char *, int);

    template <typename T> int next_record(T &);

    int lines_read() const { return _lines_read;}

    off_t tellg() const { return _iptr;}
    int seekg (off_t s);
    int  records_remaining (int = 0);
    off_t  file_size () const { return strlen(_src);}

    int push_record();
    int record_buffered() const;

    int skip_records();
    int skip_records (IW_Regular_Expression & rx, int nskip);

    int set_strip_trailing_blanks (int s);
    void set_dos (int s);

    int good () const;
    int ok () const;
    int eof () const { return '\0' == _src[_iptr];}
    int at_eof () const { return '\0' == _src[_iptr];}   // backwards compatability
};

#endif
