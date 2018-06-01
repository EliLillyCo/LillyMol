#ifndef IW_REGEXP_H
#define IW_REGEXP_H

#include "RegExp.h"

class IW_RegExp
{
  private:
    CRegExp _cregexp;

  public:
    IW_RegExp ();
    IW_RegExp (const char *);
    IW_RegExp (const char *, int);
    IW_RegExp (const const_IWSubstring &);
    IW_RegExp (const IWString &);
    ~IW_RegExp ();

    int ok () const;

    int set_pattern (const char *, int);
    int set_pattern (const char *);
    int set_pattern (const IWString &);
    int set_pattern (const const_IWSubstring &);

    int matches (const char *);
    int matches (const char *, int);
    int matches (const IWString &);
    int matches (const const_IWSubstring &);

    int matches_save_subexpressions (const char *);
    int matches_save_subexpressions (const char *, int);
    int matches_save_subexpressions (const IWString &);
    int matches_save_subexpressions (const const_IWSubstring &);

    const regmatch_t * regmatch () const;
    int number_subexpressions () const;

    int dollar (int whichdollar, const_IWSubstring & zresult) const;
    int dollar (int whichdollar, IWString & zresult) const;

    const const_IWSubstring dollar (int) const;

    int dollar (int whichdollar, int &);
    int dollar (int whichdollar, float &);
    int dollar (int whichdollar, double &);

    int s (IWString &, const const_IWSubstring &);
    int s (const IWString &, const const_IWSubstring &, IWString &);

    void set_use_extended (int);
    void set_ignore_case (int s) { _regcomp_flags |= REG_ICASE;}
};

// In place conversion

extern int iwcrex_s (IWString &, IW_RegExp &, const const_IWSubstring &);

// Result put in a new variable

extern int iwcrex_s (const IWString &, IW_RegExp &, const const_IWSubstring &, IWString &);


#endif
