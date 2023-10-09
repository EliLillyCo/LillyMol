#ifndef STRING_DATA_SOURCE_H
#define STRING_DATA_SOURCE_H

#include "re2/re2.h"

#include "Foundational/iwstring/iwstring.h"

/*
  Read data from a string rather than a file descriptor.
*/

class String_Data_Source
{
  private:
    const char * _src;

    int _iptr;

    int _lines_read;

    int _strip_trailing_blanks;

    int _dos;

    // What was returned most recently by next_record().
    const_IWSubstring _most_recent_record;

//  private functions

    void _default_values();

  public:
    String_Data_Source(const char *);
    String_Data_Source(const char *, int);

    template <typename T> int next_record(T &);

    int lines_read() const { return _lines_read;}

    off_t tellg() const { return _iptr;}
    int seekg(off_t s);
    int records_remaining() const;
    off_t file_size() const { return strlen(_src);}

    int push_record();
    int record_buffered() const;

    // Skip over `nskip` records. Note that _lines_read is incremented.
    int skip_records(int nskip);
    int skip_records(RE2 & rx, int nskip);

    void set_strip_trailing_blanks(int s) {
      _strip_trailing_blanks = s;
    }
    // Not implemented.
    void set_dos(int s);

    int good() const;
    int ok() const { return 1;}
    int eof() const { return '\0' == _src[_iptr];}
    int at_eof() const { return '\0' == _src[_iptr];}   // backwards compatability

    // For some reason, this would not link with cmake (bazel OK), so create
    // two separate methods.
    // template <typename T> int most_recent_record(T& buffer) const;
    // Return the previously returned record
    int most_recent_record(IWString& buffer) const;
    int most_recent_record(const_IWSubstring& buffer) const;
};

#endif
