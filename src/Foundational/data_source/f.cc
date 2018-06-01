#include <fstream.h>
#include <std/string.h>
#include <assert.h>
#include "string_data_source.h"

int 
bar (istream & input)
{
  assert (input.good ());

  const int maxchar = 10;
  char buffer[maxchar];

  int i = 0;
  while (input.good ())
  {
    input.get (buffer, maxchar);
    cerr << "Getline returned " << input.gcount () << " characters\n";

    cerr << "Read record '" << buffer << "'\n";
    if ('\n' == input.peek ())
      input.seekg (input.tellg () + 1);
  }

  return 1;
}

int
bb (string_data_source & input)
{
  string buffer;

  while (input.next_record (buffer))
  {
    cerr << "Read record '" << buffer << "'\n";
  }

  return 1;
}

int
main (int argc, char **argv)
{
//fstream foo (argv[1], ios::in);
  string_data_source foo (argv[1]);

  if (! foo.good ())
  {
    cerr << "Failed to create string data source from '" << argv[1] << "'\n";
    return 1;
  }

  bb (foo);

  return 0;
}
