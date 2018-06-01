#ifndef String_Data_Source_H
#define String_Data_Source_H 1

#include <std/string.h>
#include <fstream.h>

class String_Data_Source : public istream
{
  private:
    string  _buffer;

    int     _but_size;
    char *  _input_buffer;

    int     _buf_ptr;
    int     _record_buffered;
    int     _record_delimiter;
    int     _longest_record;
    int     _lines_read;
    int     _lines_which_are_returned;
    int     _strip_leading_blanks;
    int     _strip_trailing_blanks;
    int     _compress_spaces;
    int     _skip_blank_lines;
    string  _ignore_pattern_string;
    char * * _ignore_pattern_regexp;
    char * * _filter_pattern_regexp;
    string  _filter_pattern_string;
    int     _filter_match_position;
    int     _convert_to_lowercase;
    int     _convert_to_uppercase;

// private functions

    void  _default_values ();
    int   _apply_all_filters ();
    int   _fetch_complete_record (int);
    int   _fetch_record ();

  public:
    String_Data_Source ();
    ~String_Data_Source ();

    int ok () const;

    int debug_print (ostream &) const;

    int next_record (string &);
    int most_recent_record (string &);

    int lines_read () const { return _lines_read; }
    int lines_returned () const { return _lines_returned;}

    int next_record_matches (const char *);
    int push_record ();
    int skip_to   (const char *);
    int skip_past (const char *);

    void set_strip_leading_blanks () { _strip_leading_blanks = 1; }
    void set_strip_trailing_blanks () { _strip_trailing_blanks = 1; }
    void set_skip_blank_lines () { _skip_blank_lines = 1; }
    void set_compress_spaces  () { _compress_spaces = 1;}
    void set_ignore_pattern (const char *);
    void set_filter_pattern (const char *);
    void set_convert_to_lowercase () { _convert_to_lowercase = 1; _convert_to_uppercase = 0;}
    void set_convert_to_uppercase () { _convert_to_uppercase = 1; _convert_to_lowercase = 0;}
    void set_no_case_conversion   () { _convert_to_uppercase = 0; _convert_to_lowercase = 0;}

    int filter_match_position () const { return _filter_match_position;}

    int  longest_record () const { return _longest_record; }
};

#endif
