#include <stdlib.h>

#include "iwstring.h"

static void
foo (const char * s)
{
  cout << "Passed '" << s << "'\n";

  return;
}

int
main ()
{
  IWString s ("hello world");

  foo (s);

  return 0;
}
