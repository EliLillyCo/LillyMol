#ifndef IWDIGITS_H
#define IWDIGITS_H

#include "Foundational/iwstring/iwstring.h"

class IWDigits : public iwaray<IWString>
{
  private:

    IWString _leading_space;    // need not be a space
    IWString _appended_string;

//  private functions

    void _default_values ();
    void _fill_in_the_digits ();
    template <typename T, typename I> int _append_number (T & buffer, I zdigit) const;
    template <typename T, typename I> int _append_negative_number (T & buffer, I zdigit) const;

  public:
    IWDigits ();
    IWDigits (int);

    void debug_print (std::ostream &) const;

    int set_include_leading_space (int);

    int initialise (int);

//  template <typename T> int append_digit  (T & s, int v) const { return append_number (s, v);}
//  template <typename T> int append_number (T &, int) const;

    template <typename T, typename I> int append_number (T & buffer, I zdigit) const;

    const IWString & string_for_digit (int) const;

    int set_leading_string (const const_IWSubstring &);
    int set_leading_string (char s);

    int append_to_each_stored_string (const const_IWSubstring &);
    int append_to_each_stored_string (char s);
};

/*
  We can also speed things up when dealing with fractions by pre-computing the
  string representation of all fractions
*/

class Fraction_as_String
{
  private:
    int _digits;
    double _minval;
    double _maxval;
    double _dx;

    int _nbuckets;

    IWString * _fraction;

    IWString _leading_space;    // need not be a space

//  potential problems if non-zero numbers get written as 0

    int _small_non_zero_numbers_written_as_zero;

//  private functions

   void _append_number_no_string_rep(IWString & s, float f) const;
   int  _fill_string_data();

  public:
    Fraction_as_String();
    ~Fraction_as_String();

    int initialise(float minval, float maxval, int digits);

    int nbuckets() const { return _nbuckets;}

    int active() const { return nullptr != _fraction;}

    int set_include_leading_space(int);

    int append_to_each_stored_string(const const_IWSubstring &);

    int set_leading_string(const const_IWSubstring &);
    int set_leading_string(char c);

    const IWString & string_for_fraction(float f) const;

    void append_number(IWString &, float) const;

    void set_small_non_zero_numbers_written_as_zero(int s) { _small_non_zero_numbers_written_as_zero = s;}
};

#endif
