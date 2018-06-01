#ifndef IW_CMDLINE_H
#define IW_CMDLINE_H

#include <iostream>

#include "iwstring.h"

typedef int clov_magic_number_t;

class Option_and_Value
{
  friend
    std::ostream &
    operator << (std::ostream &, const Option_and_Value &);

  private:
    int _o;
    const char * _value;
    int _valid_as;

    long int _long_int_val;
    unsigned long int _unsigned_long_int_val;

    double _double_val;
    int _valid_as_double;
    clov_magic_number_t _magic;

//  private functions

    int _discern_numeric_values ();

  public:
    Option_and_Value (int, const char * = NULL);
    ~Option_and_Value ();

    int ok () const;

    char   option () const { return _o;}
    const char * value  () const { return _value;}

    template <typename T> int value (T &) const;

    int value (IWString &);
    int value (const_IWSubstring &);
    int value (char *);
};


class Command_Line : public resizable_array<const char *>
{
  friend
    std::ostream &
      operator << (const Command_Line &, std::ostream &);

  private:
    resizable_array_p<Option_and_Value> _options;
    int _some_options_start_with_dash;    // perhaps indicative of an error
    int _unrecognised_options_encountered;
    clov_magic_number_t _magic;

  public:
    Command_Line (int, char **, const char *);
    ~Command_Line ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int some_options_start_with_dash () const { return _some_options_start_with_dash;}
    int unrecognised_options_encountered () const { return _unrecognised_options_encountered;}

    Option_and_Value * ov (const char, int = 0);

    int option_present (const char) const;
    int option_count (const char) const;

    const char * option_value (const char, int = 0) const;

    template <typename T> int value (const char, T &, int = 0) const;

//  int value (const char, double &, int = 0) const;
//  int value (const char, float &, int = 0) const;
    int value (const char, IWString &, int = 0) const;
    int value (const char, const_IWSubstring &, int = 0) const;
    int value (const char, char *, int = 0) const;

//#ifdef IW_STD_STRING_DEFINED
    int value (const char, std::string &, int = 0) const;
    std::string std_string_value(const char, int = 0) const;
//#endif

    const const_IWSubstring string_value (const char, int = 0) const;

    int all_values (const char, resizable_array<const char *> &) const;
    int all_values (const char, resizable_array_p<IWString> &, int = 0) const;

//  int position (char c, int f = 0) const { return option_present (c, f) - 1;}
};

template <typename T>
int
Option_and_Value::value (T & rc) const
{
  if (NULL == _value)
    return 0;

  const_IWSubstring tmp(_value);

  return tmp.numeric_value(rc);
}

#endif
