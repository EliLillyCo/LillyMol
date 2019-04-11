#ifndef IWSTRING_H
#define IWSTRING_H 1

#if (__GNUC__ == 3)
#if defined(linux)
#include <sys/cdefs.h>
#endif
#endif

#include <cstring>
#include <iostream>

#include <assert.h>

#include "zlib.h"

#if defined(__GNUC__)
#include <sstream>
#else
#include <sstream>
#endif
#include <vector>
#include <string>

using std::cerr;
using std::endl;

#if defined (__STRING__) || defined (_CPP_STRING) || defined (_IOSTREAM_) || defined (__STD_IOSTREAM__) || defined (_STLP_IOSTREAM) || defined (_GLIBCXX_STRING)
#define IW_STD_STRING_DEFINED 1
#endif

#include "iwaray.h"

extern void strip_leading_blanks (char *);
extern void strip_trailing_blanks (char *);
extern void to_lowercase (char *, int = 0);
//extern void to_lowercase (resizable_array_p<char> *);
extern void to_uppercase (char *, int = 0);
extern void to_uppercase (resizable_array_p<char> *);
extern void no_newline            (char *);
extern int  is_double (const char *, double *);
extern int  is_int (const char *, int *);
extern int  ccount (const char *, char);
extern int  remove_blanks (char *);
extern int  compress_blanks (char *);
extern int  compress_blanks (char *, int);
//extern int find_string (const resizable_array_p<char> *, const char *);
extern const char * iwbasename (const char *);

extern char * make_copy (const char *);

extern int  words (const char *, char = ' ');

extern int iwstrncasecmp (const char *, const char *, int);

/*
  An IWString object is very much like a standard C string.
  In general, the string representation is NOT null terminated.
  If you wish to access a NULL terminated const version, use chars ().
  chars () will apply the null terminator and pass you the string.
  Note that you are getting the actual data, not a copy, so it
  will change as the underlying object changes.

  May 96: New calls, null_terminated_chars (), with the intention of
  having chars () NOT do null termination
  For now, rawchars () returns the chars.
*/

class IWString;

class const_IWSubstring
{
  private:
    const char * _data;
    int          _nchars;

  friend
    int operator == (const char * lhs, const const_IWSubstring & rhs);
  friend
    int operator != (const char * lhs, const const_IWSubstring & rhs);
  friend
    int operator == (char, const const_IWSubstring &);
  friend
    const_IWSubstring substr (const const_IWSubstring &, int, int);
  friend
    class IWString;

  public:
    const_IWSubstring ();
    const_IWSubstring (const char *);
    const_IWSubstring (const char *, int);
    const_IWSubstring (const IWString &);

#if defined (IW_STD_STRING_DEFINED)
    const_IWSubstring (const std::string &);
#endif

    ~const_IWSubstring ();

    const_IWSubstring & operator = (const char &);
    const_IWSubstring & operator = (const char *);
    const_IWSubstring & operator = (const const_IWSubstring &);
    const_IWSubstring & operator = (const IWString &);

#if defined (IW_STD_STRING_DEFINED)
    const_IWSubstring & operator = (const std::string &);
#endif

    int debug_print (std::ostream &) const;

    void make_empty () { _data = NULL; _nchars = 0;};
    void set (const char * s, int l) { _data = s; _nchars = l;}
    void set (const char * s, size_t l) { _data = s; _nchars = static_cast<int>(l);}

    int length () const { return _nchars;}
    int nchars () const { return _nchars;}

    const char * cend () const { return _data + _nchars;}
    const char * cbegin () const { return _data;}

//  Note that rawchars () does NOT null terminate the string. Caveat user!

    const char * rawchars () const { return _data;}
    const char * data () const { return _data;}    // for compatability with std::string

    const_IWSubstring substr (int, int = -1) const;

    const_IWSubstring from_to (int, int) const;
    void from_to (int, int, char *) const;
    void from_to (int, IWString &) const;
    const_IWSubstring from_to (int, const char *) const;

    const_IWSubstring before (char) const;
    const_IWSubstring after  (char) const;

    void from_to (int, int, const_IWSubstring &) const;
    void from_to (int, int, IWString &) const;

    char operator [] (int i) const { return _data[i];}

    int is_int (int &) const;

    int numeric_value (unsigned char &) const;
    int numeric_value (unsigned short &) const;
    int numeric_value (int &) const;
    int numeric_value (unsigned int &) const;
    int numeric_value (long &) const;
    int numeric_value (float &) const;
    int numeric_value (double &) const;
    int numeric_value (long long &) const;
    int numeric_value (unsigned long &) const;
    int numeric_value (unsigned long long &) const;

    int starts_with (char) const;
    int starts_with (const char *) const;
    int starts_with (const const_IWSubstring &) const;
    int starts_with (const IWString &) const;

    int ends_with (char) const;
    int ends_with (const char *) const;
    int ends_with (const const_IWSubstring &) const;
    int ends_with (const IWString &) const;

    int index  (const char) const;
    int index  (const char * s) const { return find (s);}
    int rindex (const char) const;

    int find  (const char c) const { return index (c);}
    int find  (const char *) const;
    int find  (const IWString &) const;
    int find  (const const_IWSubstring &) const;

    int ccount (char) const;

    int contains (char c) const { return index (c) + 1;}
    int contains (const char * s) const { return find (s) >= 0;}
    int contains (const IWString & s) const { return find (s) >= 0;}
    int contains (const const_IWSubstring & s) const { return find (s) >= 0;}

    int equals_ignore_case (const char *, size_t) const;
    int equals_ignore_case (const char * s) const { return equals_ignore_case (s, static_cast<int>(std::strlen (s)));}
    int equals_ignore_case (const const_IWSubstring & s) const { return equals_ignore_case (s.rawchars (), s.length ()); }
    int equals_ignore_case (const IWString & s) const;

    int strncmp (const char *, int) const;

    int strcmp (const const_IWSubstring &) const;
    int strcasecmp (const const_IWSubstring &) const;

    int matches_ignore_case (const char *, size_t, int offset = 0) const;
    int matches_ignore_case (const char * s, int offset = 0) const { return matches_ignore_case (s, static_cast<int>(std::strlen (s)), offset);}
    int matches_ignore_case (const const_IWSubstring & s, int offset = 0) const { return matches_ignore_case (s.rawchars (), s.length (), offset);}
    int matches_ignore_case (const IWString & s, int offset = 0) const;

//  In the matches_at_position methods, the first argument is the offset

    int matches_at_position (int o, const char * s, int lens) const;
    int matches_at_position (int o, const char * s) const { return matches_at_position (o, s, static_cast<int>(std::strlen (s)));}
    int matches_at_position (int o, const IWString & s) const;
    int matches_at_position (int o, const const_IWSubstring & s) const { return matches_at_position (o, s.rawchars (), s.length ());}
     
#ifdef INST_RXPOSIXH
    int matches (regex_t &) const;
#endif

    int  next (const char, int &) const;

    int  remove_leading_chars (int);
    int  remove_leading_chars (char);
    int  remove_leading_words (int, char = ' ');
    int  remove_up_to_first (char);
    void  remove_line_terminators();
    
    int chop (int = 1);

    void strip_leading_blanks ();
    void strip_trailing_blanks ();

    char operator++ ();         // prefix version
    char operator++ (int);      // postfix version

    void lengthen (int i) { _nchars += i;}

    void operator += (int);

    int iwtruncate (int);
    int shorten (int i) { return iwtruncate (i);}
    int truncate_at_first (char);     // the matching character will be lost
    int truncate_at_last  (char);     // the matching character will be lost

    void iwbasename (const_IWSubstring &, char = '/') const;
    void iwbasename (IWString &, char = '/') const;

//  Common operation to remove a suffix. We can't just chop at the last period because what
//  happens if we have something like '/path.to/file.smi'

    int remove_suffix (char = '.', char = '/');

    int nwords () const;
    int nwords (char) const;

    int nwords_single_delimiter(char) const;   // consecutive delimiters means empty word in between

    int locate_word_beginnings (resizable_array<int> &, char = ' ') const;
    int locate_word_beginnings_single_delimiter(resizable_array<int> &, char = '	') const;

    const_IWSubstring word (int, char = ' ') const;
    int word (int, IWString &, char = ' ') const;
    int word (int, const_IWSubstring &, char = ' ') const;
    int whitespace_delimited_word (int which_word, IWString & result) const;
    int whitespace_delimited_word (int which_word, const_IWSubstring & result) const;
    int wordindex (int, char = ' ') const;
    int nextword (IWString &, int &, char = ' ') const;
    int nextword (const_IWSubstring &, int &, char = ' ') const;
    int prevword (IWString &, int &, char = ' ') const;
    int prevword (const_IWSubstring &, int &, char = ' ') const;

    template <typename T> int nextword_single_delimiter(T &, int &, char = '	') const;

    int expand_environment_variables(IWString & destination) const;

    int copy_to_char_array (char *) const;

    int split (resizable_array_p<const_IWSubstring> &, char = ' ') const;
    int split (iwaray<const_IWSubstring> &, char = ' ') const;
    int split (const_IWSubstring &, char, const_IWSubstring &) const;
    int split (IWString &, char, IWString &) const;
    int split (IWString *, char = ' ') const;
    int split (std::vector<std::string> &, std::string) const;

    char last_item () const { return *(_data + _nchars - 1);}

    int balance (char, char) const;   // reports the numeric imbalance - returns 0 if balanced

    int write (int) const;     // write to an open file descriptor

    int operator == (const char * rhs) const;
    int operator != (const char * rhs) const;
    int operator == (const const_IWSubstring &) const;
    int operator != (const const_IWSubstring &) const;
    int operator == (const IWString &) const;
    int operator != (const IWString & rhs) const;
    int operator == (char) const;

//  The relational operators are implemented using strncmp

    int operator < (const const_IWSubstring &) const;
    int operator < (const IWString &) const;
    int operator > (const const_IWSubstring &) const;
    int operator > (const IWString &) const;

    int operator < (int) const;
    int operator <= (int) const;
    int operator > (int) const;
    int operator >= (int) const;
    int operator == (int) const;
    int operator != (int) const;

    int operator < (float) const;
    int operator <= (float) const;
    int operator > (float) const;
    int operator >= (float) const;
    int operator == (float) const;
    int operator != (float) const;

    int operator < (double) const;
    int operator <= (double) const;
    int operator > (double) const;
    int operator >= (double) const;
    int operator == (double) const;
    int operator != (double) const;

//  A common operation is to parse a string into 'directive=value'

    template <typename T> int split_into_directive_and_value (const_IWSubstring & d, char separator, T & dvalue) const;
};

class IWString : public resizable_array<char>
{
  private:
//  using resizable_array_base<char>::_number_elements;
//  using resizable_array_base<char>::_things;

// private functions

    void _default_values ();
    void _const_appearing_null_terminate () const;    // used for ostream <<
    void _recompute_length ();
    int  _locate_word_boundaries (int, char, int &, int &) const;
    int  _insert (int, const char *, int);
    int  _overwrite (const char *, int, int);
    void _append (int, const char *, int);

    int _common_gsub_getting_smaller (const char * zfrom, int len_from, const char * zto, int len_to);
    int _common_gsub_getting_larger (const char * zfrom, int len_from, const char * zto, int len_to);
    int _common_gsub_same_size (const char * zfrom, int len_from, const char * zto);
    int _common_gsub (const char * zfrom, int len_from, const char * zto, int len_to);

    template <typename I> void _append_int_form(I);

  friend
    int operator == (char, const IWString &);
  friend
    int operator != (char, const IWString &);
  friend
    int operator == (const char * lhs, const IWString & rhs);
  friend
    int operator != (const char * lhs, const IWString & rhs);
  friend
    IWString operator + (const char *, const IWString &);
  friend
    IWString operator + (const IWString &, const char *);
  friend
    IWString operator + (const IWString &, const IWString &);
  friend
    class const_IWSubstring;

  public:
    IWString ();
    IWString (char);
    IWString (const char *);
    IWString (const char *, int);
    IWString (const IWString &);
    IWString (const const_IWSubstring &);
    IWString (int);                // resize operation
    IWString (unsigned int);                // resize operation

#if defined (IW_STD_STRING_DEFINED)
    IWString (const std::string &);
#endif

    int ok () const;
    int good () const { return ok ();}

    int null_terminate ();

    int set (const char *, int);
    int set (const char * s, size_t n);
    int set (const char * s, u_int32_t n);

    int set_and_assume_ownership (char *, int);

    void make_empty() { if (_number_elements > 0) resize(0);}
     
    IWString & operator = (char);
    IWString & operator = (const char *);
    IWString & operator = (const const_IWSubstring &);
    IWString & operator = (const IWString &);
#if defined(__GNUC__)
    IWString & operator = (std::ostringstream &);
#endif

#if defined (IW_STD_STRING_DEFINED)
    IWString & operator = (const std::string &);
#endif

    IWString & operator = (IWString &&);

    int operator == (char) const;
    int operator != (char) const;
    int operator == (const const_IWSubstring & rhs) const;
    int operator != (const const_IWSubstring & rhs) const;
    int operator == (const char *) const;
    int operator != (const char *) const;
    int operator == (const IWString &) const;
    int operator != (const IWString &) const;

//  The relational operators are implemented using strncmp

    int operator < (const const_IWSubstring &) const;
    int operator < (const IWString &) const;
    int operator > (const const_IWSubstring &) const;
    int operator > (const IWString &) const;

    int strncpy (const char *, int);
    int strncpy (const IWString &, int);
    int strncpy (const const_IWSubstring &, int);

    int strncat (const char *, int);
    int strncat (const IWString &, int);
    int strncat (const const_IWSubstring &, int);

    int append_with_spacer (const const_IWSubstring, char = ' ');
    int append_with_spacer (const const_IWSubstring, const IWString &);

    int copy_to_char_array (char *) const;

    int strspn (const char *);

    int length () const { return _number_elements;}
    int nchars () const { return _number_elements;}

    int  next (const char, int &) const;

    void strip_leading_blanks ();
    void strip_trailing_blanks ();

    void chop (int = 1);
    void remove_leading_chars (int);
    int  remove_leading_chars (char);
    void remove_leading_chars (int, char);
    int  remove_leading_words (int, char = ' ');

    int  remove_word (int, char = ' ');
    void iwtruncate (int);
    int  truncate_at_first (char);     // the matching character will be lost
    int  truncate_at_last  (char);     // the matching character will be lost

    int remove_chars (int zfrom, int nchars);
    int remove_from_to (int zfrom, int zto);

    void iwbasename (const_IWSubstring &, char = '/') const;
    void iwbasename (IWString &, char = '/') const;

    int remove_all_these (const char *);
    int remove_up_to_first (char);

//  Common operation to remove a suffix. We can't just chop at the last period because what
//  happens if we have something like '/path.to/file.smi'

    int remove_suffix (char = '.', char = '/');

    int shift (int, char);

    int nwords () const;
    int nwords (char) const;

    int nwords_single_delimiter(char) const;   // consecutive delimiters means empty word in between

    int locate_word_beginnings (resizable_array<int> &, char = ' ') const;
    int locate_word_beginnings_single_delimiter(resizable_array<int> &, char = '	') const;

    const_IWSubstring word (int, char = ' ') const;
    int word (int, IWString &, char = ' ') const;
    int word (int, const_IWSubstring &, char = ' ') const;
    int wordindex (int, char = ' ') const;
    int nextword (IWString &, int &, char = ' ') const;
    int nextword (const_IWSubstring &, int &, char = ' ') const;
    int prevword (IWString &, int &, char = ' ') const;

    template <typename T> int nextword_single_delimiter(T &, int &, char = '	') const;

    const_IWSubstring substr (int, int = -1) const;

    const_IWSubstring from_to (int, int) const;
    void from_to (int, int, char *) const;
    const_IWSubstring from_to (int, const char *) const;

    void from_to (int, int, const_IWSubstring &) const;
    void from_to (int, int, IWString &) const;

    const_IWSubstring before (char) const;
    const_IWSubstring after  (char) const;

    int shorten (int);

    int is_int (int &) const;
    int is_double (double &) const;
    int is_hex (unsigned int &) const;

    int numeric_value (unsigned char &) const;
    int numeric_value (unsigned short &) const;
    int numeric_value (int &) const;
    int numeric_value (unsigned int &) const;
    int numeric_value (long &) const;
    int numeric_value (float &) const;
    int numeric_value (double &) const;
    int numeric_value (long long &) const;
    int numeric_value (unsigned long &) const;

    int numeric_value_fast(int &) const;

    const char * chars ();
    const char * null_terminated_chars ();
    const char * c_str() { return null_terminated_chars();}
    operator const char * () { if (0 == _number_elements) return NULL; else return null_terminated_chars ();}
    const char * rawchars () const { return _things;}

    const char * data () const { return _things;}    // for compatability with std::string

    int equals_ignore_case (const char *, size_t) const;
    int equals_ignore_case (const char * s) const { return equals_ignore_case (s, static_cast<int>(std::strlen (s)));}
    int equals_ignore_case (const const_IWSubstring & s) const { return equals_ignore_case (s.rawchars (), s.length ()); }
    int equals_ignore_case (const IWString & s) const { return equals_ignore_case (s.rawchars (), s.length ()); }

    int strcmp (const IWString &) const;
    int strcmp (const char *) const;
    int strcasecmp (const IWString &) const;
    int strcasecmp (const char *) const;

    int strncmp (const IWString &, int) const;
    int strncmp (const char *, int) const;

    int matches_ignore_case (const char *, size_t, int offset = 0) const;
    int matches_ignore_case (const char * s, int offset = 0) const { return matches_ignore_case (s, static_cast<int>(std::strlen (s)), offset);}
    int matches_ignore_case (const const_IWSubstring & s, int offset = 0) const { return matches_ignore_case (s.rawchars (), s.length (), offset);}
    int matches_ignore_case (const IWString & s, int offset = 0) const;

    int starts_with (char) const;
    int starts_with (const char *, int = -1) const;
    int starts_with (const IWString &) const;
    int starts_with (const const_IWSubstring &) const;

    int ends_with (char) const;
    int ends_with (const char *, int = -1) const;
    int ends_with (const IWString &) const;
    int ends_with (const const_IWSubstring &) const;

    int looks_like (const char *, int) const;
    int looks_like (const IWString &, int) const;
    int looks_like (const const_IWSubstring &, int) const;

//  In the matches_at_position methods, the first argument is the offset

    int matches_at_position (int o, const char * s, int lens) const;
    int matches_at_position (int o, const char * s) const { return matches_at_position (o, s, static_cast<int>(std::strlen (s)));}
    int matches_at_position (int o, const IWString & s) const { return matches_at_position (o, s.rawchars (), s.length ());}
    int matches_at_position (int o, const const_IWSubstring & s) const { return matches_at_position (o, s.rawchars (), s.length ());}

#ifdef INST_RXPOSIXH
    int matches (regex_t &) const;
#endif

    void to_lowercase ();
    void to_uppercase ();

    int compress_blanks ();

    int gsub (char, char, int = -1);
    int gsub (const char *, const char *);
    int gsub (const const_IWSubstring &, const const_IWSubstring &);
    int gsub (const IWString &, const IWString &);
    int gsub (char, const char *);
    int gsub (char, const const_IWSubstring &);
    int gsub (char, const IWString &);

    int unhtml();

    void operator += (const char * rhs);
    void operator += (char cc) {(void) add (cc);}
    void operator += (const const_IWSubstring &);
    void operator += (const IWString &);

#if defined (IW_STD_STRING_DEFINED)
    void operator += (const std::string &);
#endif

    void append_number (int);
    void append_number (unsigned int);
    void append_number (long);
    void append_number (long long);
    void append_number (unsigned long long);
    void append_number (unsigned long);
    void append_number (float);
    void append_number (double);
    void append_number (float, int);       // 2nd arg is the precision
    void append_number (double, int);      // 2nd arg is the precision
    void append_number (float, const char *);

    // Moved const around to match calls in tsclass and descriptor_file_to_01_fingerprints
    // int append_hex(const unsigned char * v, const int n);
    int append_hex(unsigned char const * v, const int n);

    void append (int, const char *);     // appends N copies of S
    void append (int, const IWString &);     // appends N copies of S

//  Often we want to build up a space separated list

    IWString & append_to_space_separated_list (const IWString & rhs);

    IWString & operator += (int f) { append_number (f); return *this;}
    IWString & operator += (float f) { append_number (f); return *this;}
    IWString & operator += (double f) { append_number (f); return *this;}

    int getline (std::istream &, char = '\n');
    int getline (int fd, char = '\n');

#ifdef ZLIB_H
    int getline (gzFile *, char = '\n');
#endif

    int write (int) const;     // write to an open file descriptor

//  Write blocks to an open file descriptor. Only write whole blocks
//  and move any unwritten data to the beginning

    int write_whole_blocks_shift_unwritten (int);

    int find (const char *) const;
    int find (const IWString &) const;
    int find (const const_IWSubstring &) const;

    int index (const char * s) const { return find (s);}
    int index (char c) const { return resizable_array<char>::index (c);}

    int ccount (char) const;

    int contains (char c) const { return resizable_array<char>::contains (c);}
    int contains (const IWString & c) const { return find (c) >= 0;}
    int contains (const char * c)     const { return find (c) >= 0;}

    int insert   (char, int);
    int insert   (const char *, int);
    int insert   (const IWString &, int);
    int insert   (const const_IWSubstring &, int);

    int overwrite   (char, int);
    int overwrite   (const char *, int);
    int overwrite   (const IWString &, int);
    int overwrite   (const const_IWSubstring &, int);

    typedef const_IWSubstring IWSubstring;

    typedef int word_iterator;

    int translate (char, char);

    int split (const_IWSubstring &, char, const_IWSubstring &) const;
    int split (IWString &, char, IWString &) const;
    int split (iwaray<const_IWSubstring> &, char = ' ') const;
    int split (iwaray<IWString> &, char = ' ') const;
    int split (resizable_array_p<const_IWSubstring> &, char = ' ') const;
    int split (resizable_array_p<IWString> &, char = ' ') const;
    int split (IWString *, char = ' ') const;

    int balance (char, char) const;    // like open and close parens, or [], or <>

    IWString & operator << (const IWString &);
    IWString & operator << (const const_IWSubstring &);
    IWString & operator << (const char *);
    IWString & operator << (char);
    IWString & operator << (int);
    IWString & operator << (long);
    IWString & operator << (long long);
    IWString & operator << (unsigned long long);
    IWString & operator << (unsigned int);
    IWString & operator << (long unsigned int);
    IWString & operator << (float);
    IWString & operator << (double);
#ifdef IW_STD_STRING_DEFINED
    IWString & operator << (const std::string & s) { this->operator+=(s); return *this;}
#endif

    int operator < (int) const;
    int operator <= (int) const;
    int operator > (int) const;
    int operator >= (int) const;
    int operator == (int) const;
    int operator != (int) const;

    int operator < (float) const;
    int operator <= (float) const;
    int operator > (float) const;
    int operator >= (float) const;
    int operator == (float) const;
    int operator != (float) const;

    int operator < (double) const;
    int operator <= (double) const;
    int operator > (double) const;
    int operator >= (double) const;
    int operator == (double) const;
    int operator != (double) const;

//  Change characters between (and including) istart to istop to a new string. Grow or shorten as needed

    int change (int istart, int istop, const char * s, int len_s);
    int change (int istart, int istop, const char * s) { return change (istart, istop, s, static_cast<int>(std::strlen (s)));}
    int change (int istart, int istop, const IWString & s) { return change (istart, istop, s.rawchars (), s.length ());}
    int change (int istart, int istop, const const_IWSubstring & s) { return change (istart, istop, s.rawchars (), s.length ());}

    template <typename T> int split_into_directive_and_value (const_IWSubstring & d, char separator, T & dvalue) const;
};

inline int
IWString::operator == (const char * rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  if (_number_elements != (int) std::strlen (rhs))    // cannot be the same if different lengths
    return 0;

  return 0 == ::strncmp (_things, rhs, _number_elements);    // loss of const OK for strcmp
}

inline int
operator != (char lhs, const const_IWSubstring & rhs)
{
  if (1 != rhs.length ())
    return 1;

  return lhs != rhs[0];
}

inline int
const_IWSubstring::operator != (const const_IWSubstring & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;
  int l2 = rhs._nchars;
  if (_nchars != l2)
    return 1;

  if (0 == ::strncmp (_data, rhs._data, _nchars))
    return 0;
  else
    return 1;
}

inline int
const_IWSubstring::operator != (const IWString & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  if (_nchars != rhs.nchars ())
    return 1;

  if (0 == ::strncmp (_data, rhs.rawchars (), _nchars))
    return 0;
  else
    return 1;
}

inline int
IWString::operator == (const IWString & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (_number_elements != rhs._number_elements)
    return 0;

//return (0 == ::strncmp (_things, rhs._things, _number_elements));
  return (0 == ::memcmp (_things, rhs._things, _number_elements));
}

inline int
IWString::operator != (const IWString & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  if (_number_elements != rhs._number_elements)
    return 1;

  return (0 != ::strncmp (_things, rhs._things, _number_elements));
}


inline int
const_IWSubstring::operator == (const const_IWSubstring & rhs) const
{
  if (_nchars != rhs._nchars)
    return 0;

  return (0 == ::strncmp (_data, rhs._data, _nchars));
}

inline int
const_IWSubstring::operator == (const IWString & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (_nchars != rhs._number_elements)
    return 0;

  return (0 == ::strncmp (_data, rhs._things, _nchars));
}

inline int
IWString::operator == (char rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  if (1 != _number_elements)
    return 0;

  return rhs == _things[0];
}

inline int
operator == (const char * lhs, const IWString & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  int lchars = static_cast<int>(std::strlen (lhs));
  if (lchars != rhs._number_elements)
    return 0;

  return 0 == ::strncmp (lhs, rhs._things, lchars);
}

inline int
operator != (const char * lhs, const IWString & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  int lchars = static_cast<int>(std::strlen (lhs));
  if (lchars != rhs._number_elements)
    return 1;

  return 0 != ::strncmp (lhs, rhs._things, lchars);
}

inline int
operator == (const char * lhs, const const_IWSubstring & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  int lchars = static_cast<int>(std::strlen (lhs));
  if (lchars != rhs._nchars)
    return 0;

  return (0 == ::strncmp (lhs, rhs._data, lchars));
}

inline int
operator != (const char * lhs, const const_IWSubstring & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  int lchars = static_cast<int>(std::strlen (lhs));
  if (lchars != rhs._nchars)
    return 1;

  return (0 != ::strncmp (lhs, rhs._data, lchars));
}

inline int 
operator != (char lhs, const IWString & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (1 != rhs._number_elements)
    return 1;

  return lhs != rhs._things[0];
}

inline int
const_IWSubstring::operator == (const char * rhs) const
{
  int l2 = static_cast<int>(std::strlen (rhs));

  if (_nchars != l2)
    return 0;

  return 0 == ::strncmp (_data, rhs, _nchars);
}

inline int 
IWString::operator == (const const_IWSubstring & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  int l2 = rhs.nchars ();
  if (_number_elements != l2)
    return 0;

  return 0 == ::strncmp (_things, rhs.rawchars (), _number_elements);
}

inline int 
IWString::operator != (const const_IWSubstring & rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  if (_number_elements != rhs.nchars ())
    return 1;

  return 0 != ::strncmp (_things, rhs.rawchars (), _number_elements);
}

inline int
operator == (char lhs, const IWString & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (1 != rhs._number_elements)
    return 0;

  return lhs == rhs._things[0];
}

inline int
operator == (char lhs, const const_IWSubstring & rhs)
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (1 != rhs._nchars)
    return 0;

  return lhs == rhs._data[0];
}

inline int
const_IWSubstring::operator == (char rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (1 != _nchars)
    return 0;

  return rhs == _data[0];
}

inline int
IWString::operator != (char rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;
  if (1 != _number_elements)
    return 1;

  return rhs != _things[0];
}

inline int
IWString::operator != (const char * rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  if (_number_elements != (int) std::strlen (rhs))    // indeed they are different if different lengths
    return 1;

  return 0 != ::strncmp (_things, rhs, _number_elements);
}

inline int
const_IWSubstring::operator != (const char * rhs) const
{
//cerr << "invoked at line " << __LINE__ << endl;

  int l2 = static_cast<int>(std::strlen (rhs));
  if (l2 != _nchars)    // the two strings are definitely not equal
    return 1;

  return 0 != ::strncmp (_data, rhs, l2);
}

inline std::ostream &
operator << (std::ostream & os, const IWString & s)
{
  assert (s.ok ());

  if (0 == s.length ())
    return os;

  return os.write (s.rawchars (), s.length ());
}

inline std::ostream &
operator << (std::ostream & os, const IWString * s)
{
  assert (s->ok ());

  if (0 == s->nchars ())
    return os;

  return os.write (s->rawchars (), s->length ());
}

inline std::ostream &
operator << (std::ostream & os, const const_IWSubstring & s)
{
  if (0 == s.nchars ())
    return os;

  return os.write (s.rawchars (), s.length ());
}

/*
  By default, the flush operator resizes the IWString object down to 0.
  If we want it to keep its allocated storage, pass flush() this arg
*/

#define IWSFD_KEEP_STORAGE_ON_FLUSH -7232

class IWString_and_File_Descriptor : public IWString
{
  private:
    int _fd;

//  If the file name ends in .gz, we gzip our output

    gzFile _gzfile;

//  private functions

    int _compress_and_write();
    int _do_resize(int);
    int _open_gzipd_stream(const char * fname);

  public:
    IWString_and_File_Descriptor();
    IWString_and_File_Descriptor(int);
    ~IWString_and_File_Descriptor();

    int write_whole_blocks_shift_unwritten();
    int write_if_buffer_holds_more_than (int s) { if (_number_elements > s) write_whole_blocks_shift_unwritten(); return 1;}
    int flush(int keep_storage = 0);   // by default, we size our buffer back to 0

    int good () const { return (-1 != _fd && IWString::ok());}  // failed open sets _fd to -1

    ssize_t write(const char * s, size_t nchars);

    int fd() const { return _fd;}

    int active () const { return _fd > 0 || NULL != _gzfile;}
    int is_open() const { return _fd > 0 || NULL != _gzfile;}

    int open(const char *);
    int close();
};

//extern IWString_and_File_Descriptor & operator<<(int i, IWString_and_File_Descriptor & os) { return os.append_int(i);}
//extern IWString_and_File_Descriptor & operator<<(unsigned int i, IWString_and_File_Descriptor & os) { return os.append_unsigned_int(i);}

#if defined (IW_STD_STRING_DEFINED)

/*inline std::string &
operator = (std::string & lhs, const const_IWSubstring & rhs)
{
  lhs.assign (rhs.rawchars (), rhs.length ();

  return lhs;
}*/

inline void
operator += (std::string & lhs, const const_IWSubstring & rhs)
{
  lhs.append (rhs.rawchars (), rhs.length ());

  return;
}

#endif

#ifdef _DB_CXX_H_
extern IWString & operator << (const Dbt &, IWString_and_File_Descriptor &);
#endif

//extern IWString & operator << (IWString &, const IWString &);
//extern IWString & operator << (IWString &, const char *);
//extern IWString & operator << (IWString &, const const_IWSubstring &);
//extern IWString & operator << (IWString &, char);
//extern IWString & operator << (IWString &, int);
//extern IWString & operator << (IWString &, float);

extern IWString operator + (const char *, const IWString &);
extern IWString operator + (const IWString &, const char *);
extern IWString operator + (const IWString &, const IWString &);

#if defined(__GNUC__)
extern IWString & operator << (IWString &, std::ostringstream &);
#endif

extern int operator >> (std::istream &, IWString &);

extern const_IWSubstring substr (const IWString &, int, int = -1);
extern const_IWSubstring substr (const const_IWSubstring &, int, int = -1);

extern const_IWSubstring from_to (const IWString &, int, int);

extern const_IWSubstring from_to (const const_IWSubstring &, int, int);

extern void append_digit (IWString &, int);
extern int change_suffix (IWString &, const IWString &);

extern int is_int (const IWString &, int &);
extern int is_int (const const_IWSubstring &, int &);

extern int expand_environment_variables (const char *, IWString &);

extern int compress_descriptor_file_record (const const_IWSubstring & buffer,
                                 IWString & destination);

extern void set_default_iwstring_float_concatenation_precision (int s);
extern void set_default_iwstring_double_concatenation_precision (int s);

extern int char_name_to_char(IWString & s);
extern void char_name_to_char_usage(const IWString & s);

/*
  Append a number to an ostream or IWString with a given width
  Initial implementation just for ints
*/

template <typename O, typename T> int append_number(O &, T, int);
template <> int append_number(std::ostream &, int, int);
template <> int append_number(IWString &, int, int);

#if defined(IW_IMPLEMENTATIONS_EXPOSED) || defined(APPEND_INT_FORM_IMPLEMENTATION)

#ifndef APPEND_INT_FORM_DONE
#define APPEND_INT_FORM_DONE

template <typename I>
void
IWString::_append_int_form (I znumber)
{
  if (0 == znumber)
  {
    resizable_array<char>::add('0');
    return;
  }

  char buffer[24];
  int ndigits = 0;
  while (znumber)
  {
    I j = znumber % 10;
    buffer[ndigits++] = '0' + j;  // digits[j];
    znumber = znumber / 10;
  }

  int chars_needed = ndigits;

  if (_number_elements + chars_needed > _elements_allocated)
    resize (_number_elements + chars_needed);

  for (int i = 0; i < ndigits; i++)
  {
    _things[_number_elements++] = buffer[ndigits - i - 1];
  }

  return;
}
#endif
#endif

#if defined(SPLIT_DV_IMPLEMENTATION) || defined(IW_IMPLEMENTATIONS_EXPOSED)

template <typename T>
int
IWString::split_into_directive_and_value (const_IWSubstring & directive,
                                          char separator,
                                          T & v) const
{
  const_IWSubstring tmp (*this);

  return tmp.split_into_directive_and_value (directive, separator, v);
}

template <typename T>
int
const_IWSubstring::split_into_directive_and_value (const_IWSubstring & directive,
                                        char separator,
                                        T & v) const
{
  int i = 0;

  if (! nextword (directive, i, separator))
  {
    cerr << "IWString::split_into_directive_and_value:cannot split directive on '" << separator << "'\n";
    return 0;
  }

  const_IWSubstring tmp;

  if (! nextword (tmp, i, separator))
  {
    cerr << "const_IWSubstring::split_into_directive_and_value:cannot split value on '" << separator << "'\n";
    return 0;
  }

  if (tmp.numeric_value (v))
    return 1;

  cerr << "const_IWSubstring::split_into_directive_and_value:invalid numeric '" << tmp << "'\n";
  return 0;
}

#endif
#endif

/* arch-tag: a3c38e2b-3093-495e-a598-8ad0d05cc122 */
