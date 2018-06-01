#include "iwstring.h"

static void 
foo (const IWString & bar)
{
  cerr << "Arg is '" << bar << "'\n";

  return;
}

int
main ()
{
   foo ("hello world");

   return 0;
}
