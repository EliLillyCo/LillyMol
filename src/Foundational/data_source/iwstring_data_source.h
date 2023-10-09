#ifndef IWSTRING_DATA_SOURCE_H
#define IWSTRING_DATA_SOURCE_H 1

#include <sys/types.h>
#include <iostream>
#include <memory>

#include "re2/re2.h"

#include "Foundational/iwstring/iwzlib.h"
// IL #include "lzma.h"
#include "Foundational/iwstring/iwstring.h"

//#define STRING_DEFAULT_BUF_SIZE 4096
#define STRING_DEFAULT_BUF_SIZE 8192

class IWSDS_State;       // for saving and restoring state info

/*
  Jun 2004. Found that seek operations on gnu iostream objects have
  become horribly, horribly slow. Therefore re-implement this to
  exclude c++ i/o
*/

class iwstring_data_source
{
  private:
    int      _fd;
    int      _open;
    int      _good;
    int      _eof;

    IWString _fname;

//  _buffer is what is often used to communicate with arguments.

    IWString _buffer;

//  _read_buffer is used for reading bytes from our file descriptor.
//  What this means is that if you attempt to use the file descriptor
//  while it is being used by an iwstring_data_source object, that will
//  produce unpredictable results

    char *   _read_buffer;

    int      _lrecl;

//  As we read data into _read_buffer, we copy it record at a time into _buffer.
//  So we need to know two positions in _read_buffer

    int      _next_char_in_read_buffer_to_transfer_to_buffer;
    int      _chars_in_read_buffer;

    int      _record_buffered;
    int      _record_delimiter;
    int      _longest_record;
    int      _lines_read;
    int      _lines_which_are_returned;
    int      _strip_leading_blanks;
    int      _strip_trailing_blanks;
    int      _compress_spaces;
    int      _skip_blank_lines;

    std::unique_ptr<RE2> _ignore_pattern;
    std::unique_ptr<RE2> _filter_pattern;

    int      _convert_to_lowercase;
    int      _convert_to_uppercase;
    int      _dos;
    int      _translate_tabs;

//  Feb 2000. I want programmes to be able to work in a pipeline,
//  so let's have the ability to echo every record returned by next_record ()
//  to a stream

    std::ostream * _echo_returned_records_stream;

//  Aug 2003. I want to be able to read gzip'd files

    IW_ZLib_Wrapper _gzfile;

// private functions

    void  _default_values (int);
    void  _setup_stream (const char *);
    int   _apply_all_filters ();
    int   _fetch_record ();

    int   _copy_next_record_from_read_buffer_to_buffer ();
    int   _read_more_data_into_read_buffer ();
    int   _fetch_record_into_buffer ();
    int   _write_read_buffer (std::ostream & output, size_t & nbytes);
    int   _write_read_buffer (IWString_and_File_Descriptor & output, size_t & nbytes);
    int   _copy_read_buffer_to_destination (void * destination, int & nbytes);
    int   _matches_ignore_pattern(const const_IWSubstring& buffer) const;
    int   _matches_filter_pattern(const const_IWSubstring& buffer) const;

    int   _save_state (IWSDS_State &);
    int   _restore_state (IWSDS_State &);

    size_t _copy_raw_bytes (void *, const size_t); 

protected:
	
	// for stringbuffer.
	bool _isstringbuffer;
	const char *_stringbuffer;
	int _stringbuffer_size;
	bool _isstringbuffer_loaded;


  public:
	  bool get_isstringbuffer()
	  {
		  return _isstringbuffer;
	  }

    iwstring_data_source ();
    iwstring_data_source (int);    // an already opened file descriptor
    iwstring_data_source (const char *,   int = STRING_DEFAULT_BUF_SIZE);
    iwstring_data_source (const IWString &, int = STRING_DEFAULT_BUF_SIZE);
		
	iwstring_data_source (bool isstringbuffer, const char *stringbuffer, int stringbuffer_size);

    ~iwstring_data_source ();

//  If created without a file name, use open () to open

    int open (const char *);
    int open (IWString & f) { return open (f.null_terminated_chars ());}

    int do_close ();

    int fd() const { return _fd;}

    int is_open () const { return _open;}

    int is_pipe () const;

    const IWString & fname () const { return _fname;}

    void set_dos (int d) { _dos = d;}
    void set_record_delimiter (char d);

    int record_buffered() const { return _record_buffered;}

    void set_echo_returned_records_stream (std::ostream & os) { _echo_returned_records_stream = &os;}
    void turn_off_echoing () { _echo_returned_records_stream = NULL;}
    std::ostream * echo_stream () const { return _echo_returned_records_stream;}

    int good () const;
    int ok () const;
    int eof () const { return _eof;}
    int at_eof () const { return _eof;}   // backwards compatability

    int debug_print (std::ostream &) const;

    off_t    tellg () const;
    int      seekg (off_t, int = 0);
    off_t    file_size ();

    void reset_record_counter () { _lines_read = 0;}

    template <typename T> int next_record (T &);

    int most_recent_record (IWString &);

    int lines_read () const { return _lines_read; }
    int next_record_matches (const char *);
    int push_record ();
    int skip_to   (const char *);
    int skip_past (const char *);

    void set_strip_leading_blanks (int s=1) { _strip_leading_blanks = s; }
    void set_strip_trailing_blanks (int s=1) { _strip_trailing_blanks = s; }
    void set_skip_blank_lines (int s=1) { _skip_blank_lines = s; }
    void set_compress_spaces  (int s=1) { _compress_spaces = s;}
    int  set_filter_pattern (const const_IWSubstring &);
    int  set_ignore_pattern (const const_IWSubstring &);
    void set_convert_to_lowercase () { _convert_to_lowercase = 1; _convert_to_uppercase = 0;}
    void set_convert_to_uppercase () { _convert_to_uppercase = 1; _convert_to_lowercase = 0;}
    void set_no_case_conversion   () { _convert_to_uppercase = 0; _convert_to_lowercase = 0;}
    void set_translate_tabs (int s) { _translate_tabs = s;}

    int  longest_record () const { return _longest_record; }
    int  records_remaining (int = 0);

    int  at_least_X_records_remaining (int x);

    int  grep (char);
    int  grep (const const_IWSubstring &);
    int  grep (RE2 &);
    int  grep (int n, RE2 *, int *);   // look for N regular expressions at once

    int  count_records_starting_with(const const_IWSubstring &);   // the most common thing we do with regexps

    int  echo (std::ostream &, size_t);    // echo's bytes
    int  echo (IWString_and_File_Descriptor &, size_t);    // echo's bytes

    int echo_records (std::ostream & os, int necho);    // echo's records
    int echo_records (IWString_and_File_Descriptor & os, int necho);    // echo's records

    int skip_records (int nskip);

    int skip_records (RE2& rx, int nskip);

    int read_bytes (void *, size_t);

    size_t copy_raw_bytes (void *, const size_t);    // does a save and restore of the state, so it will not advance the file pointer - just got too complicated...
};

#endif
