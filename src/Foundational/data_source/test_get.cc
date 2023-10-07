#include <stdlib.h>
#include <iostream>
#include <fstream>

int
main ()
{
  ifstream input ("foo.bar", ios::in);

  if (! input.good ())
  {
    cerr << "Cannot open input file\n";
    return 1;
  }

  char buffer[256];

  char newline = '\n';

  while (1)
  {
    input.get (buffer, sizeof (buffer), newline);
    int nchars = input.gcount ();

    if (nchars > 0)
    {
      cerr << "Read " << nchars << " good is " << input.good () << endl;
      cerr.write (buffer, nchars);
      cerr << endl;

      input.get ();

      continue;
    }

    cerr << "Read 0 characters, eof " << input.eof () << " good " << input.good () << endl;
    if (input.eof ())
    {
      cerr << "Endfile, done\n";
      return 0;
    }

    if (! input.good ())
    {
      cerr << "Stream is not good\n";
      input.clear ();
      cerr << "After clearing eof " << input.eof () << endl;
    }

    char nextchar = input.get ();
    if (newline == nextchar)
    {
      cerr << "Looks like a blank line came in\n";
      continue;
    }
    else
    {
      cerr << "Next char is '" << nextchar << "'\n";
    }

    if (! input.good ())
    {
      cerr << "Input stream not good, breaking\n";
      input.clear ();
    }

    if (input.eof ())
    {
      cerr << "Encountered EOF\n";
      return 0;
    }
  }

  return 0;
}
