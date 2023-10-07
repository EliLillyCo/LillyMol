// Tester for the string class

#include <stdlib.h>

#include <fstream>

#include "iwstring.h"

#include "should_match.h"

using std::cerr;
using std::endl;

extern int verbose;

template <typename Char_Type>
int
_should_match(const Char_Type & result,
               int nchars,
               const char * expected,
               const char * caller)
{
  if (verbose > 1)
    cerr << "Testing, result = '" << result << "' (" << nchars << "), expected = '" << expected << 
            "', origin " << caller << endl;

  if (result == expected)
    return 1;

  cerr << "Failure, result = '" << result << "', expected = '" << expected << 
          "', origin " << caller << endl;

  abort();
  return 0;
}

template int _should_match(const IWString&,int,const char*,const char*);
template int _should_match(const const_IWSubstring&,int,const char*,const char*);
template int _should_match(const char *const&,int,const char*,const char*);

int
should_match(const IWString & result, const char * expected,
              const char * caller)
{
  assert (result.ok ());

  return _should_match(result, result.nchars (), expected, caller);
}

int
should_match(const const_IWSubstring & result, const char * expected,
              const char * caller)
{
  return _should_match(result, result.nchars (), expected, caller);
}

int
should_match(const char * result, const char * expected,
              const char * caller)
{
  return _should_match(result, strlen (result), expected, caller);
}
