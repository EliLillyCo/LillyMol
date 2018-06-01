#include "iwstring.h"

static int
foo (const char * s)
{
  cerr << "Called char * version\n";

  return 1;
}

static int
foo (const IWString & s)
{
  cerr << "Called iwstring version\n";

  return 1;
}

int
main ()
{
  char buffer[9];

  (void) foo (buffer);

  return 0;
}
