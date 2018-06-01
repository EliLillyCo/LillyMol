#include <String.h>

int
main ()
{
  String foo;

  readline (cin, foo);
  foo.gsub ("^  *", "");

  cerr << "Foo is '" << foo << "'\n";
  cerr << "Length of foo is " << foo.length () << endl;

  return 0;
}
