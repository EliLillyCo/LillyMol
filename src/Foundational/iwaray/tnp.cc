/*
  Tester for the iterator
*/

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

int
main (int argc, char ** argv)
{
  const int fsize = 10;
  resizable_array_np<int> foo (fsize);

  for (int i = 0; i < fsize; i++)
    foo.add (i);

  int j;
  while (foo.next (j))
  {
    cerr << "Got value " << j << " from foo\n";
  }

  return 0;
}
