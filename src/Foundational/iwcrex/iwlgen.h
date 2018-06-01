#ifndef IW_RX_LGEN_H
#define IW_RX_LGEN_H

#include "iwstring.h"

/*
  This class is used with IW_Regular_Expression_Template
*/

class IW_lgen
{
  private:
    char * _compiled_regular_expression;

  public:
    IW_lgen ();
    ~IW_lgen ();

    int ok () const;

    int set_pattern (const char *, int);
    int set_pattern (const char *);
    int set_pattern (const IWString &);
    int set_pattern (const const_IWSubstring &);

    int matches (const char *, int);
    int matches (const char *);
    int matches (const IWString &);
    int matches (const const_IWSubstring &);

//  The following are undefined for our lgen interface

    int matches_save_subexpressions (const char *, int);
    int matches_save_subexpressions (const char *);
    int matches_save_subexpressions (const IWString &);
    int matches_save_subexpressions (const const_IWSubstring &);

    int number_subexpressions () const;

    int dollar (int whichdollar, const_IWSubstring & zresult) const;

    const const_IWSubstring dollar (int d) const;

    int dollar (int whichdollar, int &);
    int dollar (int whichdollar, float &);
    int dollar (int whichdollar, double &);

    int s (IWString &, const const_IWSubstring &);
    int s (const IWString &, const const_IWSubstring &, IWString &);

    void set_use_extended (int);
};

#endif
