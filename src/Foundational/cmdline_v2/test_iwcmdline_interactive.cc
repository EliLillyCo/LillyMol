/*
  Tester for command line
*/

#include <stdlib.h>

#include "cmdline_v2.h"
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Interactive tester for command line options\n";
  cerr << " -int <int>     integer option\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

template <typename T>
int
echo_option_values (const Command_Line_v2 & cl,
                    const char * opt,
                    T & v)
{
  cerr << cl.option_count (opt) << " instances of '" << opt << "'\n";
  int rc = 0;

  int i = 0;
  while (cl.value (opt, v, i++))
  {
    cerr << " " << (i - 1) << " instance of '" << opt << "', value '" << v << "'\n";
    rc++;
  }

  return rc;
}

template int echo_option_values (const Command_Line_v2 &, const char *, const_IWSubstring &);
template int echo_option_values (const Command_Line_v2 &, const char *, int &);
template int echo_option_values (const Command_Line_v2 &, const char *, IWString &);
template int echo_option_values (const Command_Line_v2 &, const char *, float &);
template int echo_option_values (const Command_Line_v2 &, const char *, unsigned int &);
template int echo_option_values (const Command_Line_v2 &, const char *, double &);

static int
test_iwcmdline_interactive (int argc, char ** argv)
{
  IWString initialiser;

  const char * s = "-v";
  cerr << "Single letter option '" << s << "'\n";
  initialiser << s;

  s = "-foo";
  cerr << "Multi character option '" << s << "'\n";
  initialiser << s;

  s = "-string=s";
  cerr << "string qualified option '" << s << "'\n";
  initialiser << s;

  s = "-dir=dir";
  cerr << "directory option '" << s << "'\n";
  initialiser << s;

  s = "-sfile=sfile";
  cerr << "file option '" << s << "'\n";
  initialiser << s;

  s = "-xfile=xfile";
  cerr << "executable file option '" << s << "'\n";
  initialiser << s;

  s = "-int=int";
  cerr << "integer option '" << s << "'\n";
  initialiser << s;

  s = "-ipos=ipos";
  cerr << "positive integer option '" << s << "'\n";
  initialiser << s;

  s = "-ipos=ipos";
  cerr << "positive integer option '" << s << "'\n";
  initialiser << s;

  s = "-uint=uint";
  cerr << "unsigned integer option '" << s << "'\n";
  initialiser << s;

  s = "-float=float";
  cerr << "float option '" << s << "'\n";
  initialiser << s;

  s = "-fraction=fraction";
  cerr << "double option '" << s << "'\n";
  initialiser << s;

  s = "-close=close";
  cerr << "double option '" << s << "'\n";
  initialiser << s;

  Command_Line_v2 cl (argc, argv, initialiser.null_terminated_chars ());

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    return 3;
  }

  if (! cl.good ())
  {
    cerr << "Error, one or more options invalid\n";
    return 4;
  }

  cerr << cl.option_count ('v') << " instances of -v option\n";

  if (cl.option_count ('v') != cl.option_count ("v"))
  {
    cerr << "Huh, mismatch between quoted forms!, " << cl.option_count ('v') << " vs " << cl.option_count ("v") << endl;
    return (3);
  }

  cerr << cl.option_count ("foo") << " instances of -foo option\n";

  cerr << cl.option_count ("dir") << " instances of -dir option\n";

  const_IWSubstring tmp;
  echo_option_values (cl, "dir", tmp);
  echo_option_values (cl, "sfile", tmp);
  echo_option_values (cl, "xfile", tmp);

  float f;
  echo_option_values (cl, "float", f);
  echo_option_values (cl, "fraction", f);

  double d;
  echo_option_values (cl, "double", d);

  int i;
  echo_option_values (cl, "int", i);
  echo_option_values (cl, "ipos", i);

  unsigned int u;
  echo_option_values (cl, "uint", u);

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_iwcmdline_interactive (argc, argv);

  return rc;
}
