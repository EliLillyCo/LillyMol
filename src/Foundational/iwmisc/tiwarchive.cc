#include <stdlib.h>
#include <iostream>
using namespace std;

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwarchive.h"

static int
tiwarchive (int argc, char ** argv)
{
  iwarchive<int> foo;
  
  cerr << "After construction, foo is " << foo << endl;

  foo.add (1);
  foo.add (2);

  cerr << "After adding, foo is " << foo << endl;

  iwarchive<int> bar = foo;

  foo.add (3);
  bar.add (4);

  cerr << "After construction, bar is " << bar << endl;

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = tiwarchive (argc, argv);

#ifdef USE_IWMALLOC
  iwmalloc_terse_malloc_status (stderr);
#endif

  return rc;
}
