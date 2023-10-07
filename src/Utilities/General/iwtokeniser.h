#ifndef IWTOKENISER_H
#define IWTOKENISER_H

#include "Foundational/iwstring/iwstring.h"

/*
  Want a common interface for fetching tokens from a line.
  We need to do different things depending on whether there are
  quoted tokens, or how consecutive delimiters are treated.
  Possibly this should be a base class and derived classes, but too lazy for that.
*/

class IWTokeniser
{
  private:
    const const_IWSubstring _buffer;
    int _i;    //used by nextword

    int _empty_fields_valid;

    int _quoted_tokens;
    int _in_quote;    // as we scan the buffer

    char _sep;

//  private functions

    template <typename S> int _next_token_with_quotes(S & token);

  public:
    IWTokeniser(const const_IWSubstring & s);

    void set_sep(char c) { _sep = c;}
    void set_quoted_tokens(const int s) { _quoted_tokens = s;}
    void set_empty_fields_valid(const int s) { _empty_fields_valid = s;}

    template <typename S> int next_token(S &);
};

#endif
