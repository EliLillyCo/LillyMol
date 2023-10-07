#include <stdlib.h>
#include <iostream>

#include "iwaray.h"
#include "iwrandom.h"

static int
test_extending (int argc, char ** argv)
{
  extending_resizable_array<int> foo;

  for (int i = 0; i < 10; i++)
  {
    int j = intbtwij (0, 100);
    foo[j]++;
  }

  int nf = foo.number_elements ();
  for (int i = 0; i < nf; i++)
    if (foo[i])
      cerr << "foo[" << i << "] = " << foo[i] << endl;
  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_extending (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
