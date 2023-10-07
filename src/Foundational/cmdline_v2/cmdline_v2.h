#ifndef IWNCMDLINE_H
#define IWNCMDLINE_H

#include "Foundational/iwstring/iwstring.h"

class CmdLineDefault
{
  private:
  public:
    CmdLineDefault ();
};

extern CmdLineDefault cmdlinedefault;

class CmdLine_Option_and_Value
{
  private:
    const char * _opt;
    int _valid_as;
    IWString _value;

  public:
    CmdLine_Option_and_Value ();

    int initialise (const char * opt, int vas, const IWString & v);
    int initialise (const char * opt);

    int is_option (const char * s, int lens) const;

    template <typename T> int value (T &) const;
    int value (IWString & v) const { v = _value; return 1;}
    int value (const_IWSubstring & v) const { v = _value; return 1;}

#ifdef IW_STD_STRING_DEFINED
    int value (std::string & v) const { v = _value.rawchars (); return 1;}
#endif

    const IWString & value () const { return _value;}
};

class Option_Type;

class Command_Line_v2 : public resizable_array<const char *>
{
  private:
    int _unrecognised_options_encountered;

    int _good;

    int _number_option_value_pairs;

    CmdLine_Option_and_Value * _clov;

// private functions

  int _constructor (const int argc,
                       char * const * argv,
                       const const_IWSubstring & s,
                       Option_Type * ot,
                       const int number_options,
                       const CmdLineDefault & c);
    int _load_rest_into_array (int argc,
                                char * const * argv,
                                int istart);
    int _initialise_old_form (int argc,
                                    char * const * argv,
                                    const const_IWSubstring & s);
    int _initialise_old_form (int argc,
                                    char * const * argv,
                                    const const_IWSubstring & s,
                                    Option_Type * ot);

  public:
    Command_Line_v2 (int, char * const *, const char *, const CmdLineDefault & = cmdlinedefault);
    ~Command_Line_v2 ();

    int good () const { return _good;}

    int unrecognised_options_encountered () const { return _unrecognised_options_encountered;}

    int option_present (const char * s, int lens) const;
    int option_present (char c) const;
    int option_present (const char * s) const { return option_present (s, static_cast<int> (strlen (s)) );}
    int option_present (const IWString & s) const {return option_present (s.rawchars (), s.length ());}
    int option_present (const const_IWSubstring & s) const {return option_present (s.rawchars (), s.length ());}

#ifdef IW_STD_STRING_DEFINED
    int option_present (const std::string & s) const {return option_present (s.data (), static_cast<int> (s.length ()) );}
#endif

    int option_count (const char * s, int lens) const;
    int option_count (char c) const { return option_count (&c, 1);}
    int option_count (const char * s) const { return option_count (s, static_cast<int> (strlen (s)) );}
    int option_count (const IWString & s) const { return option_count (s.rawchars (), s.length ());}
    int option_count (const const_IWSubstring & s) const { return option_count (s.rawchars (), s.length ());}

#ifdef IW_STD_STRING_DEFINED
    int option_count (const std::string & s) const { return option_count (s.data (), static_cast<int> (s.length()) );}
#endif

    const char * option_value (int lens, const char * s, int which_one) const;
    const char * option_value (char c, int which_one = 0) const { return option_value (1, &c, which_one);}
    const char * option_value (const char * s, int which_one = 0) const { return option_value ( static_cast<int> (strlen (s)), s, which_one);}
    const char * option_value (const IWString & s, int which_one = 0) const { return option_value (s.length (), s.rawchars (), which_one);}
    const char * option_value (const const_IWSubstring & s, int which_one = 0) const { return option_value (s.length (), s.rawchars (), which_one);}

#ifdef IW_STD_STRING_DEFINED
    const char * option_value (const std::string & s, int which_one = 0) const { return option_value ( static_cast<int> (s.length()), s.data (), which_one);}
#endif

    template <typename T> int value (int lens, const char * s, T &, int which_one = 0) const;
    template <typename T> int value (const char * s, T & v, int which_one = 0) const { return value ( static_cast<int> (strlen (s)), s, v, which_one);}
    template <typename T> int value (char c, T & v, int which_one = 0) const { return value (1, &c, v, which_one);}

    const const_IWSubstring string_value (int lens, const char * s, int which_one = 0) const;
    const const_IWSubstring string_value (const char c, int which_one = 0) const { return string_value (1, &c, which_one);}
    const const_IWSubstring string_value (const char * s, int which_one = 0) const { return string_value ( static_cast<int> (strlen (s)), s, which_one);}
    const const_IWSubstring string_value (const IWString & s, int which_one = 0) const { return string_value (s.length (), s.rawchars (), which_one);}
    const const_IWSubstring string_value (const const_IWSubstring & s, int which_one = 0) const { return string_value (s.length (), s.rawchars (), which_one);}

#ifdef IW_STD_STRING_DEFINED
    const const_IWSubstring string_value (const std::string & s, int which_one = 0) const { return string_value ( static_cast<int> (s.length ()), s.data (), which_one);}
#endif
};

#endif
