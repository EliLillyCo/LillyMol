#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include <stdio.h>

const int DEFAULT_BUF_SIZE = 512;

class data_source 
{
  friend
    data_source * new_data_source (const char *, int);
  friend
    data_source * new_data_source (FILE *, int);

  private:
    FILE *_fp;
    char *_buffer;
    int   _buf_size;
    int   _record_buffered;
    int   _longest_record;
    int   _lines_read;
    int   _at_eof;
    int   _strip_newlines;
    int   _strip_leading_blanks;
    int   _strip_trailing_blanks;
    int   _skip_blank_lines;
    char  *_ignore_pattern;
    int   _strlen_ignore_pattern;
    char  *_filter_pattern;
    int   _strlen_filter_pattern;
    int   _match_filter_at_start;
    char  *_filter_match_position;
    int   _created_from_fp;
    int   _convert_to_lowercase;
    int   _convert_to_uppercase;

    void  _default_values ();
    int   _line_is_blank () const;
    int   _apply_all_filters ();
    int   _got_complete_record (char *);
    int   _fetch_record ();
    int   _handle_eof_processing ();

  public:
    data_source (int = DEFAULT_BUF_SIZE);
    ~data_source ();

    const char * next_record ();
    const char * most_recent_record ();

    int lines_read () const { return _lines_read; }
    int at_eof () const     { return _at_eof; }
    int not_at_eof () const { return ! _at_eof; }
    int next_record_matches (const char *);
    int push_record ();
    int skip_to   (const char *);
    int skip_past (const char *);

    void set_strip_newlines () { _strip_newlines = 1; }
    void set_strip_leading_blanks ()  { _strip_leading_blanks = 1; }
    void set_strip_trailing_blanks () { _strip_trailing_blanks = 1; }
    void set_skip_blank_lines () { _skip_blank_lines = 1; }
    void set_ignore_pattern (const char *);
    void set_filter_pattern (const char *);
    void set_convert_to_lowercase () { _convert_to_lowercase = 1; _convert_to_uppercase = 0;}
    void set_convert_to_uppercase () { _convert_to_uppercase = 1; _convert_to_lowercase = 0;}
    const char * filter_match_position () const { return _filter_match_position;}

    int  longest_record () const { return _longest_record; }
};

extern data_source * new_data_source (const char *, int = DEFAULT_BUF_SIZE);
extern data_source * new_data_source (FILE *, int = DEFAULT_BUF_SIZE);

#endif
